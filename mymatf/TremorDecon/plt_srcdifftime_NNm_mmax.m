function [f,cnt,Nn,Nnn,frac,fracn,mmaxnonzero,mmaxnonzeron]=...
  plt_srcdifftime_NNm_mmax(f,nbst,imp,nsrc,impn,nsrcn,mmax,sps,typepltnoi,pltflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,cnt,Nn,Nnn,frac,fracn,ampplt,dtarvlplt,amppltn,dtarvlpltn]=...
%   plt_srcdifftime_NNm_mmax(f,nbst,imp,nsrc,impn,nsrcn,mmax,sps,typepltnoi)
%
% This function is to plot the diff time distribution of srcs, binned by amp
% first, for each m smaller than mmax, of N & N-m source pairs. 
% Amp here is the median amp of the cluster defined by consecutive events
% within N & N-m source pairs. 
% Deal with data and synthetic noise only, NOT for synthetics.
% Deal with a series of m < mmax
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/11/04
% Last modified date:   2023/11/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('pltflag',1);

if mmax < 3
  xran = [0 2];
else
  xran = [0 2*ceil(0.25*mmax+0.125)];
end
binwdt = 0.05;
nx = round(xran(2)/binwdt)+1;
edges = xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2;
cnt = xran(1): binwdt: xran(2);

for m = 1:mmax
  dtcut(m) = 0.25*m+0.125;

  [~,dtarvlplt]=med_amp_incluster(nbst,imp,nsrc,m);
  N=histcounts(dtarvlplt(:,1)/sps,edges,'normalization','count');
  Nn(:,m) = N./length(dtarvlplt(:,1));
  [~,dtarvlpltn]=med_amp_incluster(nbst,impn,nsrcn,m);
  N=histcounts(dtarvlpltn(:,1)/sps,edges,'normalization','count');
  Nnn(:,m) = N./length(dtarvlplt(:,1));
  
  frac(m) = sum(dtarvlplt(:,1)/sps<=dtcut(m))/length(dtarvlplt(:,1));
  if typepltnoi == 1
    fracn(m) = sum(dtarvlpltn(:,1)/sps<=dtcut(m))/length(dtarvlpltn(:,1));
  elseif typepltnoi == 2
    fracn(m) = (sum(dtarvlplt(:,1)/sps<=dtcut(m))-sum(dtarvlpltn(:,1)/sps<=dtcut(m))) / ...
      (length(dtarvlplt(:,1))-length(dtarvlpltn(:,1)));
  end
end
if isempty(find(frac==0,1))  %if the trial mmax is not big enough to exhaust the list
  mmaxnonzero = mmax;
  mmaxnonzeron = mmax;
  warning('mmax is not big enough, better increase it');
else
  mmaxnonzero = find(frac==0,1)-1;
  mmaxnonzeron = find(fracn==0,1)-1;
end

% frac = frac(1: mmaxnonzero);
% fracn = fracn(1: mmaxnonzeron);

if pltflag
  color = jet(mmaxnonzero);

  %for data
  for m = 1:mmaxnonzero
    %%%distribution of diff time
    ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    p(m) = plot(ax,cnt,Nn(:,m),'-','Color',color(m,:),'linew',1);
    label{m} = sprintf('m=%d',m);
    text(ax,0.95,0.7-0.04*m,sprintf('Frac w/i %.3f s: %.3e \n',...
      dtcut(m),frac(m)),'HorizontalAlignment','right','Units','normalized',...
      'FontSize',8);
    ylabel(ax,'Normalized count');
    xlabel(ax,sprintf('Arrival time difference N to N-m (s)'));
    xlim(ax,xran);
    % ylim(ax,yran);
    yran = ax.YLim;
    longticks(ax,2);
    title(ax,'Data');
    hold(ax,'off');
    legend(ax,p,label,'NumColumns',ceil(mmaxnonzero/5));
  end

  %for noise, or data-noise
  p = []; label = [];
  for m = 1:mmaxnonzeron
    ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    if typepltnoi == 1
      p(m) = plot(ax,cnt,Nnn(:,m),'-','Color',color(m,:),'linew',1);
      title(ax,'Synthetic noise');
    elseif typepltnoi == 2
      p(m) = plot(ax,cnt,Nn(:,m)-Nnn(:,m),'-','Color',color(m,:),'linew',1);%
      title(ax,'Data - Synthetic noise');
    end
    label{m} = sprintf('m=%d',m);
    text(ax,0.95,0.7-0.04*m,sprintf('Frac w/i %.3f s: %.3f \n',...
      dtcut(m),fracn(m)),'HorizontalAlignment','right','Units','normalized',...
      'FontSize',8);
    ylabel(ax,'Normalized count');
    xlabel(ax,sprintf('Arrival time difference N to N-m (s)'));
    xlim(ax,xran);
    ylim(ax,yran);
    longticks(ax,2);
    hold(ax,'off');
    legend(ax,p,label,'NumColumns',ceil(mmaxnonzeron/5));
  end
end

