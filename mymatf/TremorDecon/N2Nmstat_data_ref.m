% N2Nmstat_data_ref.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the script to carry out a bunch of N to N-m statisical analysis
% to deconvoluted catalog from the real data. The codes are cut from 
% 'deconv_ref_4s_exp_4thsta.m' because they are getting more and more 
% complicated and different from the original scope of 
% 'deconv_ref_4s_exp_4thsta.m'
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/10/23
% Last modified date:   2023/10/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
% format short e   % Set the format to 5-digit floating point
clear
clc
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, resol] = pixelperinch(1);

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');

cutout = 'ellipse';
ttol = 35;
ntol = 3;
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.pgc002.',cutout(1:4)));
% trange = trange(1:end-1,:);
tlen = trange(:,3)-trange(:,2);
tlensum = sum(tlen);

modname = 'timeoff_plfit_4thsta_160sps.mat';
planefit = load(strcat(rstpath, '/MAPS/',modname));
rmse = planefit.gof{4}.rmse;
offmax = round(2.0*rmse);

sps = 160;

%%%load data
savefile = 'deconv_stats4th_allbstsig.mat';
load(strcat(rstpath, '/MAPS/',savefile));
savefile = 'deconv_stats4th_allbstnoi.mat';
load(strcat(rstpath, '/MAPS/',savefile));
 
keyboard

%%
%%%param for secondary sources removed
locxyprojall = allbstsig.locxyprojall;
tarvlsplstall = allbstsig.impindepall(:,1);
nsrc = allbstsig.nsrcraw;
imp = allbstsig.impindepall;
locxyprojalln = allbstnoi.locxyprojall;
tarvlsplstalln = allbstnoi.impindepall(:,1);
nsrcn = allbstnoi.nsrcraw;
impn = allbstnoi.impindepall;
supertstr = 'Secondary sources removed';
fnsuffix = [];
% dtarvlnn1 = allbstsig.dtarvlnn1all;
% dtarvlnn2 = allbstsig.dtarvlnn2all;
% dtarvlnn3 = allbstsig.dtarvlnn3all;
% dtarvlnn4 = allbstsig.dtarvlnn4all;
% dtarvlnn5 = allbstsig.dtarvlnn5all;
% dtarvlnn6 = allbstsig.dtarvlnn6all;
% dtarvlnn7 = allbstsig.dtarvlnn7all;
% dtarvlnn8 = allbstsig.dtarvlnn8all;
% dtarvlnn1n = allbstnoi.dtarvlnn1all;
% dtarvlnn2n = allbstnoi.dtarvlnn2all;
% dtarvlnn3n = allbstnoi.dtarvlnn3all;
% dtarvlnn4n = allbstnoi.dtarvlnn4all;
% dtarvlnn5n = allbstnoi.dtarvlnn5all;
% dtarvlnn6n = allbstnoi.dtarvlnn6all;
% dtarvlnn7n = allbstnoi.dtarvlnn7all;
% dtarvlnn8n = allbstnoi.dtarvlnn8all;

% %%%param for further checked at KLNB
% locxyprojall = allbstsig.locxyproj4thall;
% tarvlsplstall = allbstsig.impindep4thall(:,1);
% nsrc = allbstsig.nsrc;
% imp = allbstsig.impindep4thall;
% locxyprojalln = allbstnoi.locxyproj4thall;
% tarvlsplstalln = allbstnoi.impindep4thall(:,1);
% nsrcn = allbstnoi.nsrc;
% impn = allbstnoi.impindep4thall;
% supertstr = 'Further checked at KLNB';
% fnsuffix = '4th';

typepltnoi = 1; %plot noise
% typepltnoi = 2; %plot data -noise

% dtarvlplt = dtarvlnn1;
% dtarvlpltn = dtarvlnn1n;
% m = 5;

%% for any m, fraction of unique events inside clusters w/i a diff time cut, wrt. all catalog for data
fracsrc2all = [];
mmax = 15;
% nsep = 14;
for m = 1:15
  
  dcut = 0.25*m+0.125;
  nsrcsep = nsrc-m;
  nsrcsep(nsrcsep<0) = 0;
  ndcutpair = zeros(size(trange,1), 1);
  ndcutsrc = zeros(size(trange,1), 1);
  ndcutsrc2 = zeros(size(trange,1), 1);
  imppair = cell(size(trange,1), 1);
  imppairuni = cell(size(trange,1), 1);
  impcont = cell(size(trange,1), 1);
  impcontuni = cell(size(trange,1), 1);
  indcont = cell(size(trange,1), 1);
  indcontuni = cell(size(trange,1), 1);
  
  for i = 1: size(trange,1)
    if nsrc(i) == 0
      continue
    end
    ist = sum(nsrc(1:i-1))+1;
    ied = ist+nsrc(i)-1;
    impi = imp(ist:ied,:);
    %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
    dtarvl = srcdistNtoNm(impi(:,1),impi(:,7:8),mmax);
    %   dtarvlnnsepall = [dtarvlnnsepall; dtarvl{nsep}];
    if isempty(dtarvl{m})
      continue
    end
    impbf = impi(1:end-m,:);
    impaf = impi(1+m:end,:);
    if ~isequal(size(impbf,1),length(dtarvl{m}))
      disp('Check');
    end
    ind = find(dtarvl{m}/sps <= dcut);
    ndcutpair(i,1) = length(ind);
    tmp = [];
    for j = 1: length(ind)
      tmp = [tmp; impbf(ind(j),:); impaf(ind(j),:)];
    end
    imppair{i,1} = tmp;
    % impsort = sortrows(imppair,1);
    imppairuni{i,1} = unique(tmp,'rows','stable');
    % impall = [impall; impuni];
    ndcutsrc(i,1) = size(imppairuni{i,1},1);
    
    tmp = [];
    for j = 1: length(ind)
      tmp = [tmp; impi(ind(j):ind(j)+m,:)];
    end
    impcont{i,1} = tmp;
    impcontuni{i,1} = unique(tmp,'rows','stable');
    % impall = [impall; impuni];
    ndcutsrc2(i,1) = size(impcontuni{i,1},1);
    
    tmp = [];
    for j = 1: length(ind)
      tmp = [tmp; reshape(ind(j):ind(j)+m, [], 1)];
    end
    indcont{i,1} = tmp;
    indcontuni{i,1} = unique(tmp,'rows','stable');
  end
  
  imppairunia = cat(1,imppairuni{:});
  ndcutsrcall = sum(ndcutsrc);
  fracsrcall = ndcutsrcall/sum(nsrc)*100;
  
  impcontunia = cat(1,impcontuni{:});
  ndcutsrc2all(m,1) = sum(ndcutsrc2);
  fracsrc2all(m,1) = ndcutsrc2all(m,1)/sum(nsrc)*100;
  
  ndcutpairall(m,1) = sum(ndcutpair);
  fracpairall = ndcutpairall(m,1)/sum(nsrcsep)*100;
  
end

for m = 1:15
  dcut = 0.25*m+0.125;
  fprintf('%d clusters of %d consecutive events (%d/%d unique) w/i %.3f s \n',...
    ndcutpairall(m,1), m+1, ndcutsrc2all(m,1), sum(nsrc), dcut);
end

fracsrc2all = [100; fracsrc2all];
fprintf('%.3f \n',fracsrc2all);
dfracsrc2all = fracsrc2all(1:end-1) - fracsrc2all(2:end);
fprintf('%.3f \n',dfracsrc2all);

f=initfig;
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
plot(ax,0:1:15,fracsrc2all,'ko-','linew',1,'markersize',4);
plot(ax,1:15,dfracsrc2all,'ro-','linew',1,'markersize',4);
xlabel(ax,'m');
ylabel(ax,'Frac of srcs');
title(ax,'Frac of srcs in cluster of N & N-m whose diff time w/i 0.25*m+0.125 s');
legend(ax,'inclusive','exclusive');

keyboard

%% for any m, bin by median amp for all events within N&N-m, fraction of diff time measurements w/i dtcut  
f1 = initfig(10,4.5,1,2); %initialize fig
tit=supertit(f1.ax,strcat({'Diff. arrival between N & N-m, '},supertstr));
movev(tit,0.3);

f2 = initfig(10,9,2,2); %initialize fig
tit=supertit(f2.ax,strcat({'Diff. arrival between N & N-m, binned by amp, '},supertstr));
movev(tit,0.2);

mmax = 15;
color = jet(mmax);
for m = 1:15
  dtcut = 0.25*m+0.125;
  if m < 3
    xran = [0 2];
  else
    xran = [0 2*ceil(dtcut)];
  end
  
  nsrcsep = nsrc-m;
  nsrcsep(nsrcsep<0) = 0;
  ampdt = []; % the median amp for all sources in between N & N-m pair
  dtarvlplt = [];
  for i = 1: size(trange,1)
    if nsrc(i) == 0
      continue
    end
    ist = sum(nsrc(1:i-1))+1;
    ied = ist+nsrc(i)-1;
    impi = imp(ist:ied,:);
    %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
    dtarvl = srcdistNtoNm(impi(:,1),impi(:,7:8),mmax);
    %   dtarvlnnsepall = [dtarvlnnsepall; dtarvl{nsep}];
    if isempty(dtarvl{m})
      continue
    end
    dtarvlplt = [dtarvlplt; dtarvl{m}];
    ampcont = [];
    for j = 1: nsrc(i)-m
      impcont = impi(j: j+m,:);
      ampcont(j,1) = median(mean(impcont(:,[2 4 6]),2));
    end
    ampdt = [ampdt; ampcont];
  end
  
  ampdtn = [];
  dtarvlpltn = [];
  for i = 1: size(trange,1)
    if nsrcn(i) == 0
      continue
    end
    ist = sum(nsrcn(1:i-1))+1;
    ied = ist+nsrcn(i)-1;
    impi = impn(ist:ied,:);
    %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
    dtarvln = srcdistNtoNm(impi(:,1),impi(:,7:8),mmax);
    %   dtarvlnnsepall = [dtarvlnnsepall; dtarvl{nsep}];
    if isempty(dtarvln{m})
      continue
    end
    dtarvlpltn = [dtarvlpltn; dtarvln{m}];
    ampcont = [];
    for j = 1: nsrcn(i)-m
      impcont = impi(j: j+m,:);
      ampcont(j,1) = median(mean(impcont(:,[2 4 6]),2));
    end
    ampdtn = [ampdtn; ampcont];
  end
  ampplt = ampdt;
  amppltn = ampdtn;

  %%%distribution of diff time
  binwdt = 0.05;
  % dtcut = 0.25*nsep+0.125;
  if m < 3
    xran = [0 2];
  else 
    xran = [0 2*ceil(dtcut)];
  end
  nx = round(xran(2)/binwdt)+1;
  edges = xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2;
  cnt = xran(1): binwdt: xran(2);
  N=histcounts(dtarvlplt/sps,edges,'normalization','count');
  Nn = N./length(dtarvlplt);
  frac(m) = sum(dtarvlplt/sps<=dtcut)/length(dtarvlplt);
  N=histcounts(dtarvlpltn/sps,edges,'normalization','count');
  Nnn = N./length(dtarvlplt);
  ax=f1.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  plot(ax,cnt,Nn,'-','Color',color(m,:),'linew',1);
  text(ax,0.95,0.7-0.04*m,sprintf('Frac w/i %.3f s: %.2f \n',dtcut,frac(m)),'HorizontalAlignment','right',...
    'Units','normalized','FontSize',8);
  ylabel(ax,'Normalized count');
  xlabel(ax,sprintf('Diff. arrival between sources N and N-m (s)'));
  xlim(ax,xran);
  % ylim(ax,yran);
  yran = ax.YLim;
  longticks(ax,2);
  title(ax,'Data');
  hold(ax,'off');

  ax=f1.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  if typepltnoi == 1
    fracn(m) = sum(dtarvlpltn/sps<=dtcut)/length(dtarvlpltn);
    plot(ax,cnt,Nnn,'-','Color',color(m,:),'linew',1);
    title(ax,'Synthetic noise');
  elseif typepltnoi == 2
    fracn(m) = (sum(dtarvlplt/sps<=dtcut)-sum(dtarvlpltn/sps<=dtcut)) / ...
      (length(dtarvlplt)-length(dtarvlpltn));
    plot(ax,cnt,Nn-Nnn,'-','Color',color(m,:),'linew',1);%
    title(ax,'Data - Synthetic noise');
  end
  text(ax,0.95,0.7-0.04*m,sprintf('Frac w/i %.3f s: %.2f \n',dtcut,fracn(m)),'HorizontalAlignment','right',...
    'Units','normalized','FontSize',8);
  ylabel(ax,'Normalized count');
  xlabel(ax,sprintf('Diff. arrival between sources N and N-m (s)'));
  xlim(ax,xran);
  ylim(ax,yran);
  longticks(ax,2);
  hold(ax,'off');
  

  %%%fraction of diff time w/i 'dtcut' time
  ax=f2.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  %%%%%%%%%% if bin by amp with a equal number
  nbin = 5;
  [ampbin,indbin,n] = binxeqnum(ampplt,nbin);
  binwdt = 0.05;
  nx = round(xran(2)/binwdt)+1;
  Nbn = zeros(nx, nbin);
  edges = xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2;
  cnt = xran(1): binwdt: xran(2);
  for i = 1: nbin
    ind = indbin{i};
    dtarvli = dtarvlplt(ind);
    ampbincnt(i,m) = median(ampbin{i});
    N=histcounts(dtarvli/sps,edges,'normalization','count');
    Nbn(:,i) = N./mode(n);
    fracb(i,m) = sum(dtarvli/sps<=dtcut)/mode(n);
  end
  %%%%%%%%%% if bin by amp with a equal number
  Nbnm = mean(Nbn,2);
  p(m) = plot(ax,log(ampbincnt(:,m)),fracb(:,m),'-','Color',color(m,:),'linew',1,...
    'marker','o','markersize',4,'markerfacec',color(m,:));
  label{m} = sprintf('m=%d',m);
  xlabel(ax,'Median log_{10}{amp}');
  % ylabel(ax,'Median diff. arrival (s)');
  ylabel(ax,'Frac. of diff. arrival w/i 0.25*m+0.125 s');
  yran=[0 1];
  ylim(ax,yran);
  title(ax,'Data');
  
  ax=f2.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  plot(ax,log(ampbincnt(:,m)),fracb(:,m),'-','Color',color(m,:),'linew',1,...
    'marker','o','markersize',4,'markerfacec',color(m,:));
  yran=[0 0.1];
  ylim(ax,yran);

  ax=f2.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  %%%%%%%%%% if bin by amp with a equal number
  [ampbinn,indbinn,nn] = binxeqnum(amppltn,nbin);
  Nbnn = zeros(nx, nbin);
  for i = 1: nbin
    ind = indbinn{i};
    dtarvlin = dtarvlpltn(ind);
    ampbinncnt(i,m) = median(ampbinn{i});
    N=histcounts(dtarvlin/sps,edges,'normalization','count');
    Nbnn(:,i) = N./mode(n);
    if typepltnoi == 1
      fracbn(i,m) = sum(dtarvlin/sps<=dtcut)/mode(nn);
    elseif typepltnoi == 2
      fracbn(i,m) = (fracb(i,m)*mode(n)-sum(dtarvlin/sps<=dtcut)) / (mode(n)-mode(nn));
      % fracn(i) = (frac(i)*mode(n)-sum(dtarvlin/sps<=dtcut)) / mode(n);
    end
  end
  %%%%%%%%%% if bin by amp with a equal number
  Nbnnm = mean(Nbnn,2);
  if typepltnoi == 1
    plot(ax,log(ampbinncnt(:,m)),fracbn(:,m),'-','Color',color(m,:),'linew',1,...
    'marker','o','markersize',4,'markerfacec',color(m,:));
    title(ax,'Synthetic noise');
  elseif typepltnoi == 2
    plot(ax,log(ampbincnt(:,m)),fracbn(:,m),'-','Color',color(m,:),'linew',1,...
    'marker','o','markersize',4,'markerfacec',color(m,:));%
    title(ax,'Data - Synthetic noise');
  end
  xlabel(ax,'Median log_{10}{amp}');
  ylabel(ax,'Frac. of diff. arrival w/i 0.25*m+0.125 s');
  yran=[0 1];
  ylim(ax,yran);
  
  ax=f2.ax(4); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  if typepltnoi == 1
    plot(ax,log(ampbinncnt(:,m)),fracbn(:,m),'-','Color',color(m,:),'linew',1,...
    'marker','o','markersize',4,'markerfacec',color(m,:));
  elseif typepltnoi == 2
    plot(ax,log(ampbincnt(:,m)),fracbn(:,m),'-','Color',color(m,:),'linew',1,...
    'marker','o','markersize',4,'markerfacec',color(m,:));
  end
  yran=[0 0.1];
  ylim(ax,yran);
  
end
legend(f1.ax(1),p,label,'NumColumns',3);
legend(f2.ax(1),p,label,'NumColumns',3);

keyboard

%% for a certain m, similar plots for N & N-m pairs 
%%%Median amplitude for srcs w/i the cluster
m = 1;
dtcut = 0.25*m+0.125;
nsrcsep = nsrc-m;
nsrcsep(nsrcsep<0) = 0; 
mmax = 15;
ampdt = []; % the median amp for all sources in between N & N-m pair
dtarvlplt = [];
for i = 1: size(trange,1)
  if nsrc(i) == 0
    continue
  end
  ist = sum(nsrc(1:i-1))+1;
  ied = ist+nsrc(i)-1;
  impi = imp(ist:ied,:);
  %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
  dtarvl = srcdistNtoNm(impi(:,1),impi(:,7:8),mmax);
%   dtarvlnnsepall = [dtarvlnnsepall; dtarvl{nsep}];
  if isempty(dtarvl{m})
    continue
  end
  dtarvlplt = [dtarvlplt; dtarvl{m}];
  ampcont = [];
  for j = 1: nsrc(i)-m
    impcont = impi(j: j+m,:);
    ampcont(j,1) = median(mean(impcont(:,[2 4 6]),2));
  end
  ampdt = [ampdt; ampcont];
end

ampdtn = [];
dtarvlpltn = [];
for i = 1: size(trange,1)
  if nsrcn(i) == 0
    continue
  end
  ist = sum(nsrcn(1:i-1))+1;
  ied = ist+nsrcn(i)-1;
  impi = impn(ist:ied,:);
  %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
  dtarvln = srcdistNtoNm(impi(:,1),impi(:,7:8),mmax);
%   dtarvlnnsepall = [dtarvlnnsepall; dtarvl{nsep}];
  if isempty(dtarvln{m})
    continue
  end
  dtarvlpltn = [dtarvlpltn; dtarvln{m}]; 
  ampcont = [];
  for j = 1: nsrcn(i)-m
    impcont = impi(j: j+m,:);
    ampcont(j,1) = median(mean(impcont(:,[2 4 6]),2));
  end
  ampdtn = [ampdtn; ampcont];
end
ampplt = ampdt;
amppltn = ampdtn;

%%%summarize the whole catalog, diff arrival time and fractions
p = [];
label = [];
binwdt = 0.05;
% dtcut = 0.25*nsep+0.125;
if m < 3
  xran = [0 2];
else
  xran = [0 2*ceil(dtcut)];
end
nx = round(xran(2)/binwdt)+1;
edges = xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2;
cnt = xran(1): binwdt: xran(2);
N=histcounts(dtarvlplt/sps,edges,'normalization','count');
Nn = N./length(dtarvlplt);
frac = sum(dtarvlplt/sps<=dtcut)/length(dtarvlplt);
N=histcounts(dtarvlpltn/sps,edges,'normalization','count');
Nnn = N./length(dtarvlplt);
fracn = sum(dtarvlpltn/sps<=dtcut)/length(dtarvlpltn);
f = initfig(5.5,4,1,1); %initialize fig
tit=supertit(f.ax,supertstr);
movev(tit,0.2);
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
p(1)=plot(ax,cnt,Nn,'color','b','LineWidth',1);
label{1}='Data';
p(2)=plot(ax,cnt,Nnn,'color','r','LineWidth',1);
label{2}='Synthetic noise';
fracdif = (sum(dtarvlplt/sps<=dtcut)-sum(dtarvlpltn/sps<=dtcut)) / ...
  (length(dtarvlplt)-length(dtarvlpltn));
p(3)=plot(ax,cnt,Nn-Nnn,'color','k','LineWidth',1.5);
label{3}='Data - Synthetic noise';
text(ax,0.7,0.7,sprintf('%.2f',frac),'Color','b','HorizontalAlignment','left','Units','normalized');
text(ax,0.7,0.63,sprintf('%.2f',fracn),'Color','r','HorizontalAlignment','left','Units','normalized');
text(ax,0.7,0.56,sprintf('%.2f',fracdif),'Color','k','HorizontalAlignment','left','Units','normalized');
legend(ax,p,label);
ylabel(ax,'Normalized count');
xlabel(ax,sprintf('Diff. arrival between sources N and N-%d (s)',m));
xlim(ax,xran);
% if nsep == 1
%   yran = [0 0.2];
% elseif nsep == 2
%   yran = [0 0.06];
% else
%   yran = [0 0.04];
% end
% ylim(ax,yran);
longticks(ax,2);
hold(ax,'off');
% keyboard

%%
%%%Bin by amp, then plot diff time distribution, and frac w/i some 'dtcut'
f = initfig(12,8,2,2); %initialize fig
tit=supertit(f.ax,supertstr);
movev(tit,0.2);
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
nbin = 5;
[ampbin,indbin,n] = binxeqnum(ampplt,nbin);
color = jet(nbin);
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
xlabel(ax,sprintf('Diff. arrival between sources N and N-%d (s)',m));
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
xlabel(ax,sprintf('Diff. arrival between sources N and N-%d (s)',m));
xlim(ax,xran);
ylim(ax,yran);
longticks(ax,2);
hold(ax,'off');
% keyboard

ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% plot(ax,ampbincnt,dtarvlimed,'k-'); %if bin by amp with a equal width
% p1=plot(ax,log(ampbincnt),dtarvlimed1,'k-');  %if bin by amp with a equal number
% p2=plot(ax,log(ampbincnt),dtarvlimed2,'k--'); 
% p3=plot(ax,log(ampbincnt),dtarvlimed,'k-.'); 
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
plot(ax,log(ampbincntm),fracm,'k-','linew',1,'marker','o','markersize',4,'markerfacec','k');
xlabel(ax,'Median log_{10}{amp}');
% ylabel(ax,'Median diff. arrival (s)');
ylabel(ax,sprintf('Frac. of diff. arrival w/i %.3f s',dtcut));
ylim(ax,yran);
title(ax,'Data');

ax=f.ax(4); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% plot(ax,ampbinncnt,dtarvlinmed,'k-');
% plot(ax,log(ampbinncnt),dtarvlinmed1,'k-');
% plot(ax,log(ampbinncnt),dtarvlinmed2,'k--'); 
% plot(ax,log(ampbinncnt),dtarvlinmed,'k-.');
if typepltnoi == 1 
  plot(ax,log(ampbinncntm),fracnm,'k-','linew',1,'marker','o','markersize',4,'markerfacec','k');
  title(ax,'Synthetic noise');
elseif typepltnoi == 2 
  plot(ax,log(ampbincntm),fracnm,'k-','linew',1,'marker','o','markersize',4,'markerfacec','k');
  title(ax,'Data - Synthetic noise');
end
xlabel(ax,'Median log_{10}{amp}');
% ylabel(ax,'Median diff. arrival (s)');
ylabel(ax,sprintf('Frac. of diff. arrival w/i %.3f s',dtcut));
ylim(ax,yran);

keyboard

orient(f.fig,'landscape');
fname = strcat(sprintf('nn%d%dbinsdifftime',m,nbin),fnsuffix,'.pdf');
print(f.fig,'-dpdf','-fillpage',strcat('/home/chaosong/Pictures/',fname));

%%
% %%%abs distance along min-rmse direction between each source and all others whose arrival 
% %%%separation is <=2 s
% widin = 11;  % maximum width allowed is 8.5 inches
% htin = 5;   % maximum height allowed is 11 inches
% nrow = 1;
% ncol = 2;
% f = initfig(widin,htin,nrow,ncol); %initialize fig
% 
% ax=f.ax(1);
% hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% binwdist = 0.1;
% binedge = (0: binwdist: 50*binwdist)';
% %if looking at distance in map view, use plain count or normalize by area
% %%%projection along the min-rmse direction
% %%%for data
% nsrcprop = nsrc;
% nsrcprop(nsrcprop<=2)=0;
% dlocproj2allbst = [];
% dlocprojnn1bst = [];
% dlocprojnn2bst = [];
% for i = 1: size(trange,1)
%   ist = sum(nsrc(1:i-1))+1;
%   ied = ist+nsrc(i)-1;
%   tarvlsplst = tarvlsplstall(ist:ied);
%   ist = sum(nsrcprop(1:i-1))+1;
%   ied = ist+nsrcprop(i)-1;
%   locxyproj = locxyprojall(ist:ied, :);
%   if ~isempty(tarvlsplst) && ~isempty(locxyproj) 
%     [~,dlocproj2all] = srcdistall(tarvlsplst,locxyproj,[0 2*sps]);
%     dlocproj2allbst = [dlocproj2allbst; dlocproj2all];
%     [~,dlocproj] = srcdistNtoNm(tarvlsplst,locxyproj,2);
%     dlocprojnn1bst = [dlocprojnn1bst; dlocproj{1}];
%     dlocprojnn2bst = [dlocprojnn2bst; dlocproj{2}];
%   end
% end
% dist2all = abs(dlocproj2allbst(:,1));
% [bincnt,binhgt,count,normalizer] = histbinbyarea(dist2all,binedge,'countdensity');
% % binhgt = binhgt./length(dist2all);
% binhgt = count./length(dist2all);
% p1=bar(ax,bincnt,binhgt,1,'stacked','b','facea',0.6);
% % stairs(ax,binedge,[binhgt; binhgt(end)],'k','LineWidth',1);
% plot(ax,[median(dist2all) median(dist2all)],ax.YLim,'b--','linew',1.5);
% text(ax,median(dist2all)+0.1,ax.YLim(2)-0.1*range(ax.YLim),...
%   sprintf('%.2f',median(dist2all)),'HorizontalAlignment','left');
% %%%for noise
% nsrcpropn = nsrcn;
% nsrcpropn(nsrcpropn<=2)=0;
% dlocproj2allbstn = [];
% for i = 1: size(trange,1)
% %   i
%   ist = sum(nsrcn(1:i-1))+1;
%   ied = ist+nsrcn(i)-1;
%   tarvlsplst = tarvlsplstalln(ist:ied);
%   ist = sum(nsrcpropn(1:i-1))+1;
%   ied = ist+nsrcpropn(i)-1;
%   locxyproj = locxyprojalln(ist:ied, :);
%   if ~isempty(tarvlsplst) && ~isempty(locxyproj) 
%     [~,dlocproj2all] = srcdistall(tarvlsplst,locxyproj,[0 2*sps]);
%     dlocproj2allbstn = [dlocproj2allbstn; dlocproj2all];
%   end
% end
% dist2alln = abs(dlocproj2allbstn(:,1));
% [bincnt,binhgt,count,normalizer] = histbinbyarea(dist2alln,binedge,'countdensity');
% % binhgt = binhgt./length(dist2all);
% binhgt = count./length(dist2all);
% p2=bar(ax,bincnt,binhgt,1,'stacked','r','facea',0.6);
% % stairs(ax,binedge,[binhgt; binhgt(end)],'k','LineWidth',1);
% plot(ax,[median(dist2alln) median(dist2alln)],ax.YLim,'r--','linew',1.5);
% text(ax,median(dist2alln)+0.1,ax.YLim(2)-0.2*range(ax.YLim),...
%   sprintf('%.2f',median(dist2alln)),'HorizontalAlignment','left');
% xlabel(ax,'Dist. (km) along min-scatter direc. between each source and all others with diff. arrival time \leq 2 s');
% ylabel(ax,'Normalized count');
% legend(ax,[p1,p2],'Data','Synthetic noise','location','east');
% xlim(ax,[0 4]);
% longticks(ax,2);
% hold(ax,'off');
% 
% ax=f.ax(2);
% hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% binwdt = 0.05;
% [N,edges]=histcounts(dtarvlplt/sps,'binwidth',binwdt,'normalization','count');
% N = N./length(dtarvlplt);
% bincnt = (edges(1:end-1)+edges(2:end))/2;
% p1=bar(ax,bincnt,N,1,'stacked','b','facea',0.6,'edgecolor','k');
% % p1=histogram(ax,dtarvlnn1/sps,'binwidth',binwdt,'normalization','count',...
% %   'facecolor','b','EdgeColor','k','facea',0.6);
% % p1.Values = p1.Values ./ length(dtarvlnn1);
% plot(ax,[median(dtarvlplt/sps) median(dtarvlplt/sps)],ax.YLim, '--', ...
%   'Color', 'b', 'linew', 1.5);
% text(ax,median(dtarvlplt/sps)+0.02,ax.YLim(2)-0.1*range(ax.YLim),...
%   sprintf('%.2f',median(dtarvlplt/sps)),'HorizontalAlignment','left');
% 
% [N,edges]=histcounts(dtarvlpltn/sps,'binwidth',binwdt,'normalization','count');
% N = N./length(dtarvlplt);
% bincnt = (edges(1:end-1)+edges(2:end))/2;
% p2=bar(ax,bincnt,N,1,'stacked','r','facea',0.6,'edgecolor','k');
% plot(ax,[median(dtarvlpltn/sps) median(dtarvlpltn/sps)],ax.YLim, '--', ...
%   'Color', 'r', 'linew', 1.5);
% text(ax,median(dtarvlpltn/sps)+0.02,ax.YLim(2)-0.2*range(ax.YLim),...
%   sprintf('%.2f',median(dtarvlpltn/sps)),'HorizontalAlignment','left');
% ylabel(ax,'Normalized count');
% xlabel(ax,sprintf('Diff. arrival time between sources N and N-%d (s)',nsep));
% legend(ax,[p1,p2],'Data','Synthetic noise','location','east');
% xlim(ax,[0 1.5]);
% % ylim(ax,[0 3.5]);
% longticks(ax,2);
% hold(ax,'off');
% tit=supertit(f.ax,supertstr);
% movev(tit,0.2);
% 
% orient(f.fig,'landscape');
% fname = strcat(sprintf('nn%ddiff',nsep),fnsuffix,'.pdf');
% print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));

% keyboard

%% amp, but only for N & N-m pair, ignoring sources in between
% ampafdt = [];
% ampbfdt = [];
% for i = 1: size(trange,1)
%   ist = sum(nsrc(1:i-1))+1;
%   ied = ist+nsrc(i)-1;
%   impi = imp(ist:ied,:);
%   impaf = impi(1+nsep:end,:);
%   ampaf = mean(impaf(:,[2 4 6]),2);
%   ampafdt = [ampafdt; ampaf];
%   impbf = impi(1:end-nsep,:);
%   ampbf = mean(impbf(:,[2 4 6]),2);
%   ampbfdt = [ampbfdt; ampbf];
% end
% ampdt = mean([ampbfdt ampafdt],2);
% minnum = 500;
% 
% ampafdtn = [];
% ampbfdtn = [];
% for i = 1: size(trange,1)
%   ist = sum(nsrcn(1:i-1))+1;
%   ied = ist+nsrcn(i)-1;
%   impi = impn(ist:ied,:);
%   impaf = impi(1+nsep:end,:);
%   ampaf = mean(impaf(:,[2 4 6]),2);
%   ampafdtn = [ampafdtn; ampaf];
%   impbf = impi(1:end-nsep,:);
%   ampbf = mean(impbf(:,[2 4 6]),2);
%   ampbfdtn = [ampbfdtn; ampbf];
% end
% ampdtn = mean([ampbfdtn ampafdtn],2);
% minnumn = 250;
% 
% ampplt = ampdt;
% amppltn = ampdtn;

% figure
% subplot(311)
% h=histogram(log10(ampbfdt),'BinWidth',0.25); hold on; %,'NumBins',5
% ax=gca;
% plot(ax,ax.XLim,[200 200],'k--');
% xlabel('log_{10}{Amp}');
% ylabel('Count');
% xlim([-1.5 1]);
% text(0.95,0.9,'Earlier one of the pair','Units','normalized','HorizontalAlignment','right');
% title(sprintf('Between sources N and N-%d (s)',nsep));
% subplot(312)
% h=histogram(log10(ampafdt),'BinWidth',0.25); hold on;
% ax=gca;
% plot(ax,ax.XLim,[200 200],'k--');
% xlim([-1.5 1]);
% text(0.95,0.9,'Later one of the pair','Units','normalized','HorizontalAlignment','right');
% subplot(313)
% h=histogram(log10(ampdt),'BinWidth',0.25); hold on;
% ax=gca;
% plot(ax,ax.XLim,[200 200],'k--');
% xlim([-1.5 1]);
% text(0.95,0.9,'Mean of the pair','Units','normalized','HorizontalAlignment','right');
% 
% %amp distribution of source pair 
% figure
% subplot(121)
% h=histogram(log10(ampplt),'BinWidth',0.25); hold on;
% ax=gca;
% xlim([-1.5 1]);
% title('Mean of the pair for data');
% ylabel(ax,'Count');
% xlabel(ax,sprintf('log_{10}{amp} between sources N and N-%d (s)',nsep));
% 
% subplot(122)
% h=histogram(log10(amppltn),'BinWidth',0.25); hold on;
% ax=gca;
% xlim([-1.5 1]);
% title('Mean of the pair for syn noise');
% ylabel(ax,'Count');
% xlabel(ax,sprintf('log_{10}{amp} between sources N and N-%d (s)',nsep));


%%
% figure
% subplot(121)
% eucdistnn1 = allbstsig.distarvlnn1all(:,1);
% eucdistnn2 = allbstsig.distarvlnn2all(:,1);
% 
% eucdistplt = eucdistnn1;
% dlocprojplt = dlocprojnn1bst;
% 
% ax = gca;
% hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% tmp = [dtarvlplt/sps eucdistplt log10(ampplt)];
% tmp = sortrows(tmp,3,'descend');
% scatter(ax,tmp(:,1),tmp(:,2),15,tmp(:,3),'o','filled');
% oldc = colormap(ax,'kelicol');
% newc = flipud(oldc);
% colormap(ax,newc);
% c=colorbar(ax,'SouthOutside');
% c.Label.String = strcat({'log_{10}(amp)'});
% ylabel(ax,'Distance (km)');
% xlabel(ax,sprintf('Diff. arrival time between sources N and N-%d (s)',nsep));
% xlim(ax,[0 2]);
% ylim(ax,[0 8]);
% hold(ax,'off');
% 
% subplot(122)
% ax = gca;
% hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% tmp = [dtarvlplt/sps abs(dlocprojplt(:,1)) log10(ampplt)];
% tmp = sortrows(tmp,3,'descend');
% scatter(ax,tmp(:,1),tmp(:,2),15,tmp(:,3),'o','filled');
% oldc = colormap(ax,'kelicol');
% newc = flipud(oldc);
% colormap(ax,newc);
% c=colorbar(ax,'SouthOutside');
% c.Label.String = strcat({'log_{10}(amp)'});
% ylabel(ax,'Distance along min-rmse (km)');
% xlabel(ax,sprintf('Diff. arrival time between sources N and N-%d (s)',nsep));
% xlim(ax,[0 2]);
% ylim(ax,[0 4]);
% hold(ax,'off');

