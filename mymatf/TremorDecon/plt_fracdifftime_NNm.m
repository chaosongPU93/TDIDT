function [f,Nn,fracm,ampbincntm,Nnn,fracnm,ampbinncntm]=...
  plt_fracdifftime_NNm(f,ampplt,dtarvlplt,amppltn,dtarvlpltn,sps,typepltnoi,m,nbin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,Nn,fracm,ampbincntm,Nnn,fracnm,ampbinncntm]=...
%   plt_fracdifftime_NNm(f,ampplt,dtarvlplt,amppltn,dtarvlpltn,dtcut,sps,typepltnoi,m)
%
% Similar purpose as 'plt_fracdifftime_NNm_mmax', This function is for the certain m,
% to plot the diff time distribution of srcs, binned by amp
% first, for each m smaller than mmax, of N & N-m source pairs. 
% Amp here is the median amp of the 
% cluster defined by consecutive events within N & N-m source pairs.
% Deal with data and synthetic noise only, NOT for synthetics.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/11/04
% Last modified date:   2023/11/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('m',1);
defval('nbin',5);

dtcut = 0.25*m+0.125;
if m < 3
  xran = [0 2];
else
  xran = [0 2*ceil(0.25*m+0.125)];
end

ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% %%%%%%%%%% if bin by amp with a equal width
% nbin = h.NumBins;
% color = jet(nbin-1);
% binwdt = 0.1;
% iplt = 0;
% for i = 1: nbin-1 
%   ind = find(log10(ampplt)>=h.BinEdges(i) & log10(ampplt)<h.BinEdges(i+1));
%   dtarvli = dtarvlplt(ind);
%   if length(dtarvli)>=minnum
%     iplt = iplt+1;
%     dtarvlimed(iplt) = median(dtarvli)/sps;
%     ampbincnt(iplt) = (h.BinEdges(i)+h.BinEdges(i+1))/2;
%     [N,edges]=histcounts(dtarvli/sps,'binwidth',binwdt,'normalization','count');
%     edges = edges-binwdt/2;
%   %   edges = -binwdt/2: binwdt: 2+binwdt/2;
%     [N,edges]=histcounts(dtarvli/sps,edges,'normalization','count');
%     Nn = N./length(dtarvli);
% %     p(iplt)=stairs(ax,edges,[Nn Nn(end)],'color',color(iplt,:),'LineWidth',1);
%     cnt = (edges(1:end-1)+edges(2:end))/2;
%     p(iplt)=plot(ax,cnt,Nn,'color',color(iplt,:),'LineWidth',1);
%     label{iplt} = sprintf('%.2f',ampbincnt(iplt));
%   end  
% end
% %%%%%%%%%% if bin by amp with a equal width

%%%%%%%%%% if bin by amp with a equal number
% nbin = 5;
[ampbin,indbin,n] = binxeqnum(ampplt,nbin);
color = gradientblue(nbin);
% color = plasma(nbin);
% color = gray(nbin+1);
% color = flipud(color(1:end-1,:));
% color = flipud(kelicmap(nbin));
binwdt = 0.05;
% dtcut1 = 0.5*nsep*sps;
% dtcut2 = dtcut1+0.125*sps;
% dtcut = 0.25*nsep+0.125;
nx = round(xran(2)/binwdt)+1;
Nn = zeros(nx, nbin);
edges = xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2;
cnt = xran(1): binwdt: xran(2);
for i = 1: nbin
% for i = 1: 1
  ind = indbin{i};
  dtarvli = dtarvlplt(ind);
  % dtarvlimed1(i) = median(dtarvli(dtarvli<=dtcut1))/sps;
  % dtarvlimed2(i) = median(dtarvli(dtarvli<=dtcut2))/sps;
  dtarvlimed(i) = median(dtarvli)/sps;
  dtarvlimode(i) = mode(dtarvli)/sps;
  ampbincntm(i) = median(ampbin{i});
  % [N,edges]=histcounts(dtarvli/sps,'binwidth',binwdt,'normalization','count');
  % edges = edges-binwdt/2;
  N=histcounts(dtarvli/sps,edges,'normalization','count');
  Nn(:,i) = N./mode(n);
  fracm(i) = sum(dtarvli/sps<=dtcut)/mode(n);
  % cnt = (edges(1:end-1)+edges(2:end))/2;
  p(i)=plot(ax,cnt,Nn(:,i),'color',color(i,:),'LineWidth',1);
  label{i} = sprintf('amp of %.1f',ampbincntm(i));
  % label{i} = sprintf('%d/%dth amp',i,nbin);
  % keyboard
end
%%%%%%%%%% if bin by amp with a equal number
Nnm = mean(Nn,2);
% Nnm = median(Nn,2);
% for i = 1: nbin
%   p(i)=plot(ax,cnt,Nn(:,i)-Nnm,'color',color(i,:),'LineWidth',1);
%   label{i} = sprintf('amp of %.1f - mean',ampbincnt(i));
% end
p(nbin+1)=plot(ax,cnt,Nnm,'k-','LineWidth',1.5);
label{nbin+1} = sprintf('mean');  %median
legend(ax,p,label);
ylabel(ax,'Normalized count');
xlabel(ax,sprintf('Arrival time difference N to N-%d (s)',m));
xlim(ax,xran);
if m == 1
  yran = [0 0.2];
elseif m == 2
  yran = [0 0.06];
else
  yran = ax.YLim;
end
ylim(ax,yran);
longticks(ax,2);
hold(ax,'off');
title(ax,'Data');
% keyboard

ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% %%%%%%%%%% if bin by amp with a equal width
% iplt = 0;
% for i = 1: h.NumBins-1 
%   ind = find(log10(amppltn)>=h.BinEdges(i) & log10(amppltn)<h.BinEdges(i+1));
%   dtarvlin = dtarvlpltn(ind);  
%   if length(dtarvlin)>=minnumn
%     iplt = iplt+1;
%     dtarvlinmed(iplt) = median(dtarvlin)/sps;
%     ampbinncnt(iplt) = (h.BinEdges(i)+h.BinEdges(i+1))/2;
%     [N,edges]=histcounts(dtarvlin/sps,'binwidth',binwdt,'normalization','count');
%     edges = edges-binwdt/2;
%   %   edges = -binwdt/2: binwdt: 2+binwdt/2;
%     [N,edges]=histcounts(dtarvlin/sps,edges,'normalization','count');
%     Nn = N./length(dtarvlin);
% %     p(iplt)=stairs(ax,edges,[Nn Nn(end)],'color',color(iplt,:),'LineWidth',1);
%     cnt = (edges(1:end-1)+edges(2:end))/2;
%     p(iplt)=plot(ax,cnt,Nn,'color',color(iplt,:),'LineWidth',1);
%     label{iplt} = sprintf('%.2f',ampbinncnt(iplt));
%   end  
% end
% %%%%%%%%%% if bin by amp with a equal width

%%%%%%%%%% if bin by amp with a equal number
[ampbinn,indbinn,nn] = binxeqnum(amppltn,nbin);
Nnn = zeros(nx, nbin);
for i = 1: nbin
% for i = 1: 1
  ind = indbinn{i};
  dtarvlin = dtarvlpltn(ind);
  % dtarvlinmed1(i) = median(dtarvlin(dtarvlin<=dtcut1))/sps;
  % dtarvlinmed2(i) = median(dtarvlin(dtarvlin<=dtcut2))/sps;
  dtarvlinmed(i) = median(dtarvlin)/sps;
  dtarvlinmode(i) = mode(dtarvlin)/sps;
  ampbinncntm(i) = median(ampbinn{i});
  % [N,edges]=histcounts(dtarvlin/sps,'binwidth',binwdt,'normalization','count');
  % edges = edges-binwdt/2;
  [N,edges]=histcounts(dtarvlin/sps,edges,'normalization','count');
  Nnn(:,i) = N./mode(n);
  % cnt = (edges(1:end-1)+edges(2:end))/2;
  if typepltnoi == 1 
    p(i)=plot(ax,cnt,Nnn(:,i),'color',color(i,:),'LineWidth',1);
    fracnm(i) = sum(dtarvlin/sps<=dtcut)/mode(nn);
  elseif typepltnoi == 2 
    p(i)=plot(ax,cnt,Nn(:,i)-Nnn(:,i),'color',color(i,:),'LineWidth',1);
    fracnm(i) = (fracm(i)*mode(n)-sum(dtarvlin/sps<=dtcut)) / (mode(n)-mode(nn));  
    % fracn(i) = (frac(i)*mode(n)-sum(dtarvlin/sps<=dtcut)) / mode(n);  
  end
  label{i} = sprintf('amp of %.1f',ampbinncntm(i));
end
%%%%%%%%%% if bin by amp with a equal number
Nnnm = mean(Nnn,2);
% Nnm = median(Nn,2);
% for i = 1: nbin
%   p(i)=plot(ax,cnt,Nnn(:,i)-Nnnm,'color',color(i,:),'LineWidth',1);
%   label{i} = sprintf('amp of %.1f - mean',ampbinncnt(i));
% end
label{nbin+1} = sprintf('mean');  %median
if typepltnoi == 1 
  title(ax,'Synthetic noise');
  p(nbin+1)=plot(ax,cnt,Nnnm,'k-','LineWidth',1.5);
  legend(ax,p,label);
elseif typepltnoi == 2 
  title(ax,'Data - Synthetic noise');
  p(nbin+1)=plot(ax,cnt,Nnm-Nnnm,'k-','LineWidth',1.5);
end
ylabel(ax,'Normalized count');
xlabel(ax,sprintf('Arrival time difference N to N-%d (s)',m));
xlim(ax,xran);
ylim(ax,yran);
longticks(ax,2);
hold(ax,'off');
% keyboard

ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% plot(ax,ampbincnt,dtarvlimed,'k-'); %if bin by amp with a equal width
% p1=plot(ax,log10(ampbincnt),dtarvlimed1,'k-');  %if bin by amp with a equal number
% p2=plot(ax,log10(ampbincnt),dtarvlimed2,'k--'); 
% p3=plot(ax,log10(ampbincnt),dtarvlimed,'k-.'); 
% if nsep==1
%   yran=[0.25 0.75];
%   % yran=[0.1 0.5];
% %   legend(ax,[p1 p2 p3],'w/i 0.5 s','w/i 0.625 s','all');
% elseif nsep==2
%   yran=[0.1 0.6]; 
%   % yran=[0 0.4];
% %   legend(ax,[p1 p2 p3],'w/i 1 s','w/i 1.125 s','all');
% elseif nsep==3
%   yran=[0.5 4.5];
% %   legend(ax,[p1 p2 p3],'w/i 1.5 s','w/i 1.625 s','all');
% end
yran=[0 1];
plot(ax,log10(ampbincntm),fracm,'k-','linew',1,'marker','o','markersize',4,'markerfacec','k');
xlabel(ax,'Median log_{10}{amp}');
% ylabel(ax,'Median diff. arrival (s)');
ylabel(ax,sprintf('Fraction w/i %.3f s',dtcut));
ylim(ax,yran);
title(ax,'Data');

ax=f.ax(4); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% plot(ax,ampbinncnt,dtarvlinmed,'k-');
% plot(ax,log10(ampbinncnt),dtarvlinmed1,'k-');
% plot(ax,log10(ampbinncnt),dtarvlinmed2,'k--'); 
% plot(ax,log10(ampbinncnt),dtarvlinmed,'k-.');
if typepltnoi == 1 
  plot(ax,log10(ampbinncntm),fracnm,'k-','linew',1,'marker','o','markersize',4,'markerfacec','k');
  title(ax,'Synthetic noise');
elseif typepltnoi == 2 
  plot(ax,log10(ampbincntm),fracnm,'k-','linew',1,'marker','o','markersize',4,'markerfacec','k');
  title(ax,'Data - Synthetic noise');
end
xlabel(ax,'Median log_{10}{amp}');
% ylabel(ax,'Median diff. arrival (s)');
ylabel(ax,sprintf('Fraction w/i %.3f s',dtcut));
ylim(ax,yran);

% keyboard
