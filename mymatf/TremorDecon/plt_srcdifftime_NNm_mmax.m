function [f,cnt,Nn,Nnn,frac,fracn,ampplt,dtarvlplt,amppltn,dtarvlpltn]=...
  plt_srcdifftime_NNm_mmax(f,nbst,imp,nsrc,impn,nsrcn,mmax,sps,typepltnoi)
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

color = jet(mmax);
if mmax < 3
  xran = [0 2];
else
  xran = [0 2*ceil(0.25*mmax+0.125)];
end
for m = 1:mmax
  dtcut = 0.25*m+0.125;
  [ampplt,dtarvlplt]=med_amp_incluster(nbst,imp,nsrc,m);
  [amppltn,dtarvlpltn]=med_amp_incluster(nbst,impn,nsrcn,m);

  %%%distribution of diff time
  binwdt = 0.05;
  nx = round(xran(2)/binwdt)+1;
  edges = xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2;
  cnt = xran(1): binwdt: xran(2);
  N=histcounts(dtarvlplt/sps,edges,'normalization','count');
  Nn = N./length(dtarvlplt);
  frac(m) = sum(dtarvlplt/sps<=dtcut)/length(dtarvlplt);
  N=histcounts(dtarvlpltn/sps,edges,'normalization','count');
  Nnn = N./length(dtarvlplt);
  ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  p(m) = plot(ax,cnt,Nn,'-','Color',color(m,:),'linew',1);
  label{m} = sprintf('m=%d',m);
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

  ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
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

end
legend(f.ax(1),p,label,'NumColumns',3);
