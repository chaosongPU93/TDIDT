function [f,ampbincnt,Nbn,fracb,ampbinncnt,Nbnn,fracbn]=...
  plt_fracdifftime_NNm_mmax(f,nbst,imp,nsrc,impn,nsrcn,mmax,sps,typepltnoi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,cnt,Nn,Nnn,frac,fracn,ampdt,dtarvlplt,ampdtn,dtarvlpltn]=...
%    plt_srcdifftime_NNm_mmax(f,nbst,imp,nsrc,impn,nsrcn,mmax,sps,typepltnoi)
%
% This function is to plot the diff time distribution of srcs, binned by amp
% first, for each m smaller than mmax, of N & N-m source pairs. 
% Amp here is the median amp of the 
% cluster composed by consecutive events within N & N-m source pairs.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/11/04
% Last modified date:   2023/11/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

color = jet(mmax);
if mmax < 3
  xran = [0 2];
else
  xran = [0 2*ceil(0.25*mmax+0.125)];
end
for m = 1:mmax
  dtcut = 0.25*m+0.125;
  [ampplt,dtarvlplt]=med_amp_incluster(nbst,imp,nsrc,mmax,m);
  [amppltn,dtarvlpltn]=med_amp_incluster(nbst,impn,nsrcn,mmax,m);

  %%%fraction of diff time w/i 'dtcut' time
  ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
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
  
  ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  plot(ax,log(ampbincnt(:,m)),fracb(:,m),'-','Color',color(m,:),'linew',1,...
    'marker','o','markersize',4,'markerfacec',color(m,:));
  yran=[0 0.1];
  ylim(ax,yran);

  ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
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
  
  ax=f.ax(4); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
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
legend(f.ax(1),p,label,'NumColumns',3);