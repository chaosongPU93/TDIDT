function [f,ampbincnt,Nbnall,fracdtb,Nbnnall,fracdtbn]=...
  plt_fracdt_NNm_mmax(f,nbst,imp,nsrc,impn,nsrcn,mmax,sps,typepltnoi,scale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,ampbincnt,Nbnall,fracb,Nbnnall,fracbn]=...
%   plt_fracdt_NNm_mmax(f,nbst,imp,nsrc,impn,nsrcn,mmax,sps,typepltnoi)
%
% This function is to plot the Fraction of diff time distribution of srcs, 
% binned by amp first, for each m smaller than mmax, of N & N-m source pairs. 
% Amp here is the median amp of the cluster defined by consecutive events
% within N & N-m source pairs.
% Since binning by amp first, there is no way to show the diff time distri
% itself for each m in one figure, so here only the fraction is summarized.
% Deal with data and synthetic noise only, NOT for synthetics.
% Deal with a series of m < mmax
%
% --for noise comparison, there is no need to bin by amplitude. For comparison
% to the data, just divide the total number by the number of bins used in the 
% corresponding data m. You might want the number of bins to decrease as m
% increases, since your sample is getting smaller. It's possible that at some
% value of m you might want to lump all values of m, m+1, m+2, ... into a 
% single pot to get decent statistics.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/11/04
% Last modified date:   2023/11/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('scale','linear');

if mmax < 3
  xran = [0 2];
else
  xran = [0 2*ceil(0.25*mmax+0.125)];
end

%%%for every m, bin diff time by amp, then compute the norm. count in each
%%%time bin, then save all
Nbnall = cell(mmax,1);  
Nbnnall = cell(mmax,1);
ampbincnt = cell(mmax,1);
fracdtb = cell(mmax,1); %frac of diff time measurements w/i dtcut in each amp bin
fracdtbn = cell(mmax,1); %same as above, but for noise, or, data-noise
ndtb = cell(mmax,1); %num of unique event wrt the total # evts w/i dtcut in each amp bin

for m=1:mmax
  dtcut = 0.25*m+0.125;

  %for larger m, use fewer bins for useful statistics
  if m<=2
    nbin = 10;
  elseif m>2 && m<=6
    nbin = 5;
  elseif m>6 && m<=10
    nbin = 3;
  elseif m>10 && m<=12  
    nbin = 2;
  elseif m>12
%     nbin = 1;
    break   %skip the rest, as we want to see the variation with amp, 1 bin is not enough
  end

  %%%DATA
  %%%%%%%%%% if bin by amp with a equal number
  [ampplt,dtarvlplt] = med_amp_incluster(nbst,imp,nsrc,m);  
  [ampbin,indbin,n] = binxeqnum(ampplt(:,1),nbin);  %bin by amp with same number
  binwdt = 0.05;  %width of time bin
  nx = round(xran(2)/binwdt)+1; %num of total x coords
  Nbn = zeros(nx, nbin);  %normalized counts in each time bin and amp bin
  edges = xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2; %edge of time bin
%   cnt = xran(1): binwdt: xran(2);
  medampbin = zeros(nbin,1);  %median amp of each amp bin
  fracdtbin = zeros(nbin,1); %frac of diff time measurements w/i dtcut
  ndtbin = zeros(nbin,1); %num of unique events in each amp bin
  for i = 1: nbin
    indi = indbin{i};
    dtarvli = dtarvlplt(indi,1);
    medampbin(i) = median(ampbin{i});
    N=histcounts(dtarvli/sps,edges,'normalization','count');
    Nbn(:,i) = reshape(N,[],1)./mode(n);
    fracdtbin(i) = sum(dtarvli/sps<=dtcut)/mode(n);
    ndtbin(i) = mode(n);
  end
  ampbincnt{m} = medampbin;
  fracdtb{m} = fracdtbin;  %store for each m
  ndtb{m} = ndtbin;
  Nbnall{m} = Nbn;
  %%%%%%%%%% if bin by amp with a equal number
  Nbnm = mean(Nbn,2);
  

  %%%NOISE, NO need to bin by amp
  [~,dtarvlpltn]=med_amp_incluster(nbst,impn,nsrcn,m);
  %since no amp binning, for each 'bin', simply divide by nbin  
  Nbnn = zeros(nx, nbin);
  Nn=histcounts(dtarvlpltn(:,1)/sps,edges,'normalization','count');
  fracdtbinn = zeros(nbin,1);
  for i = 1: nbin
    if typepltnoi == 1
      Nbnn(:,i) = reshape(Nn,[],1)./nbin./mode(n);
      fracdtbinn(i) = sum(dtarvlpltn(:,1)/sps<=dtcut)/length(dtarvlpltn);
    elseif typepltnoi == 2
      Nbnn(:,i) = (Nbn(:,i)*mode(n)-reshape(Nn,[],1)./nbin)./(mode(n)-length(dtarvlpltn)/nbin);
      fracdtbinn(i) = (fracdtbin(i)*mode(n)-sum(dtarvlpltn(:,1)/sps<=dtcut)/nbin) / ...
        (mode(n)-length(dtarvlpltn)/nbin);
    end
  end
  fracdtbn{m} = fracdtbinn;
  Nbnnall{m} = Nbnn;
  Nbnnm = mean(Nbnn,2);

  if ~(length(fracdtb{m})==1 &&  sum(fracdtbn{m})==0)
    mmaxnonzero = m;
  end

end

color = jet(mmaxnonzero);
for m = 1:mmaxnonzero
  %plot data
  ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  tmp1 = ampbincnt{m};
  tmp2 = fracdtb{m};
  tmp3 = fracdtbn{m};
  tmp4 = ndtb{m};
  tmp1 = tmp1(tmp2>0);
  tmp2 = tmp2(tmp2>0);
  tmp4 = tmp4(tmp2>0);
  if strcmp(scale,'log')
    tmp2 = log10(tmp2);
%     ylblstr = 'log_{10}{Frac. of diff. arrival w/i 0.25*m+0.125 s}';
    ylblstr = 'log_{10}{Fraction}';
    tmp3 = log10(tmp3);
    ylim(ax,[-2.5 0]);
%     ylim(ax,[-3 0]);
  else
%     ylblstr = 'Frac. of diff. arrival w/i 0.25*m+0.125 s';
    ylblstr = 'Fraction';
    ylim(ax,[0 1]);
  end
%   p(m) = plot(ax,log10(tmp1),tmp2,'-','Color',color(m,:),'linew',1,...
%     'marker','o','markersize',4,'markerfacec',color(m,:));
  p(m) = plot(ax,log10(tmp1),tmp2,'-','Color',color(m,:),'linew',1);
  scatter(ax,log10(tmp1),tmp2,log10(tmp4)*10,color(m,:),'filled');
  label{m} = sprintf('m=%d',m);
  xlabel(ax,'log_{10}{Median amp.}');
  ylabel(ax,ylblstr);
  xlim(ax,[-0.9 0.1]);
%   title(ax,'Data');
  longticks(ax,2);
%   if typepltnoi == 1
%     scatter(ax,abs_loc_on_axis(ax.XLim,0.98),tmp3(1),20,color(m,:),'^','LineWidth',1);
%   end
%   text(ax,abs_loc_on_axis(ax.XLim,0.02),tmp2(1),label{m});

%   if strcmp(scale,'linear')
%     %plot data zoom-in
%     ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%     p(m) = plot(ax,log10(tmp1),tmp2,'-','Color',color(m,:),'linew',1,...
%       'marker','o','markersize',4,'markerfacec',color(m,:));
%     ylim(ax,[0 0.15]);
%     xlim(ax,[-0.9 0.1]);
%     longticks(ax,2);
% %     if typepltnoi == 1
% %       scatter(ax,abs_loc_on_axis(ax.XLim,0.98),tmp3(1),20,color(m,:),'^','LineWidth',1);
% %     end
% %     text(ax,abs_loc_on_axis(ax.XLim,0.02),tmp2(1),label{m},'fontsize',8);
%     
%     % convert data units to global figure units
%     [xs,ys] = ds2nfu(f.ax(2),f.ax(2).XLim(1),f.ax(2).YLim(2));
%     [xe,ye] = ds2nfu(f.ax(1),f.ax(1).XLim(2),f.ax(2).YLim(2));
%     % plot two dashed lines denoting the zoom-in effect
%     annotation('line',[xs,xe],[ys,ye],'color',[.5 .5 .5],'linestyle','--',...
%       'LineWidth',.5);
%     % convert data units to global figure units
%     [xs,ys] = ds2nfu(f.ax(2),f.ax(2).XLim(1),f.ax(2).YLim(1));
%     [xe,ye] = ds2nfu(f.ax(1),f.ax(1).XLim(2),f.ax(2).YLim(1));
%     % plot two dashed lines denoting the zoom-in effect
%     annotation('line',[xs,xe],[ys,ye],'color',[.5 .5 .5],'linestyle','--',...
%       'LineWidth',.5);
% 
%   end

% %   if typepltnoi == 2
%     ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%     if typepltnoi == 1
%       plot(ax,log10(ampbincnt(:,m)),log10(fracbn(:,m)),'-','Color',color(m,:),'linew',1,...
%       'marker','o','markersize',4,'markerfacec',color(m,:));
%       title(ax,'Synthetic noise');
%     elseif typepltnoi == 2
%       plot(ax,log10(ampbincnt(:,m)),log10(fracbn(:,m)),'-','Color',color(m,:),'linew',1,...
%       'marker','o','markersize',4,'markerfacec',color(m,:));%
%       title(ax,'Data - Synthetic noise');
%     end
%     xlabel(ax,'Median log_{10}{amp}');
%     ylabel(ax,'Frac. of diff. arrival w/i 0.25*m+0.125 s');
% %     ylim(ax,[0 1]);
%   xlim(ax,[-1 0.2]);
%     
%     ax=f.ax(4); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%     if typepltnoi == 1
%       plot(ax,log10(ampbincnt(:,m)),log10(fracbn(:,m)),'-','Color',color(m,:),'linew',1,...
%       'marker','o','markersize',4,'markerfacec',color(m,:));
%     elseif typepltnoi == 2
%       plot(ax,log10(ampbincnt(:,m)),log10(fracbn(:,m)),'-','Color',color(m,:),'linew',1,...
%       'marker','o','markersize',4,'markerfacec',color(m,:));
%     end
% %     ylim(ax,[0 0.15]);
%   xlim(ax,[-1 0.2]);
% %   end
end

if strcmp(scale,'log')
  legend(f.ax(1),p,label,'NumColumns',4,'Location','south');
else
  legend(f.ax(1),p,label,'NumColumns',4,'Location','north');
end
  
%  keyboard 

