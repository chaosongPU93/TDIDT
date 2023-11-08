function [f1,f2,Nn,frac,Nnn,fracn,fracdif]=...
  plt_difftime_NNm(f1,f2,dtarvlplt,dtarvlpltn,dtcut,sps,typepltnoi,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,Nn,frac,Nnn,fracn,fracdif]=...
%   plt_difftime_NNm(f,dtarvlplt,dtarvlpltn,dtcut,sps,m)
%
% This function is to plot the diff time distribution of srcs, for a certain
% m, for data and synthetics noise, plotting them together in one panel or
% separately.
% Text shown the fraction of diff time measurements w/i a dtcut.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/11/04
% Last modified date:   2023/11/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('m',1);

%%%summarize the whole catalog, diff arrival time and fractions
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

ax=f1.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
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


ax=f2.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
plot(ax,cnt,Nn,'color','k','LineWidth',1);
text(ax,0.7,0.7,sprintf('%.2f',frac),'Color','b','HorizontalAlignment','left','Units','normalized');
title(ax,'Data');
yran = ax.YLim;
ylabel(ax,'Normalized count');
xlabel(ax,sprintf('Diff. arrival between sources N and N-%d (s)',m));
xlim(ax,xran);
longticks(ax,2);
hold(ax,'off');

ax=f2.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
if typepltnoi == 1
  plot(ax,cnt,Nnn,'color','k','LineWidth',1);
  text(ax,0.7,0.7,sprintf('%.2f',fracn),'Color','k','HorizontalAlignment','left','Units','normalized');
  title(ax,'Synthetic noise');
elseif typepltnoi == 2
  plot(ax,cnt,Nn-Nnn,'color','k','LineWidth',1);
  text(ax,0.7,0.7,sprintf('%.2f',fracdif),'Color','k','HorizontalAlignment','left','Units','normalized');
  title(ax,'Data - Synthetic noise');
end
xlim(ax,xran);
ylim(ax,yran);
longticks(ax,2);
hold(ax,'off');
