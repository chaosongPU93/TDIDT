function [f] = plt_rcccatvsrcc(irccran,ircccat,rcccat,nwin,sps,ircc,rcc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_rcccatvsrcc(irccran,ircccat,rcccat,nwin,sps,ircc,rcc)
%
% Plot the comparison between concatenated rcc based on alignments 
% of each short windows and original rcc based on a single 
% alignment of the whole window. 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/03/11
% Last modified date:   2024/03/11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

color=jet(nwin);

widin = 8.3;  % maximum width allowed is 8.5 inches
htin = 3;   % maximum height allowed is 11 inches
nrow = 2;
ncol = 1;
f = initfig(widin,htin,nrow,ncol); %initialize fig
xran = [0.08 0.98]; yran = [0.15 0.95];
xsep = 0.05; ysep = 0.04;
optaxpos(f,nrow,ncol,xran,yran,xsep,ysep);

ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
xlim(ax,[0 ircccat(end)]/sps);
ylim(ax,[-1 1]);
for iwin = 1: nwin
  ist = irccran(iwin,1);
  ied = irccran(iwin,2);
  plot(ax,ircccat(ist:ied)/sps,rcccat(ist:ied),'Color',color(iwin,:),'linew',1);
  plot(ax,[ied ied]/sps,ax.YLim,'--','Color',[.3 .3 .3]);
end
longticks(ax,4);
nolabels(ax,1);

ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
p1=plot(ax,ircc/sps,rcc,'k-','linew',1);
p2=plot(ax,ircccat/sps,rcccat,'r-','linew',1);
legend(ax,[p1 p2],'Whole-window','Concatenated','Orientation','horizontal',...
  'Location','south');
xlim(ax,[0 ircccat(end)]/sps);
ylim(ax,[-1 1]);
xlabel(ax,'Time (s)');
ylabel(ax,'Running CC');
longticks(ax,4);

fname = 'rcccatvsrcc.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

