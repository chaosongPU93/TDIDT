function f=plt_002lfe(green,greenf,sps,saveflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the function to plot the broadband 'green' and bandpassed 'greenf'
% templates, similar to figure 1 of Song&Rubin2024, but only templates w/o
% the seismograms created by 'plt_002lfe_signalzoom.m'.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/08/07
% Last modified date:   2024/08/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('saveflag',1);  %default is short-win detections

nsta = size(greenf,2);
% lwlet = size(greenf,1);
lwlet = 12*sps;

%%%compose a summary figure for each data win and save it, so that we will have a feeling for
%%%all migrations
widin = 8.4;  % maximum width allowed is 8.5 inches
htin = 2;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f=initfig(widin,htin,nrow,ncol);

figxran = [0.08 0.94]; figyran = [0.25 0.95];
figxsep = 0.03; figysep = 0.04;
axpos=optaxpos(f,nrow,ncol,figxran,figyran,figxsep,figysep);

color = ['r';'b';'k';'c'];

%%%LFE templates, broadband
ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on'); 
for i = 1: nsta
  p(i)=plot(ax,(1:lwlet)/sps, green(1:lwlet,i), '-','Color',color(i,:),'linew',1); 
end
% legend(ax,p,{'PGC','SSIB','SILB'},'Location','southeast','NumColumns',3,...
%   'fontsize',8);
text(ax,0.02,0.1,'Broadband','Units','normalized','HorizontalAlignment','left',...
  'FontSize',10);
text(ax,0.98,0.9,'LFE templates at family 002','Units','normalized',...
  'HorizontalAlignment','right','FontSize',10);
text(ax,0.02,0.90,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1);
xlim(ax,[0 lwlet/sps]);  
% ylim(ax,yran/2);
ylim(ax,[-0.5 0.5]);
% shrink(ax,range(xzoom)*sps/lwlet,1);
% ax.Position(1)=axpos(1,1);
xticks(ax,0: 2: lwlet/sps);
longticks(ax,2); 
% nolabels(ax,1);
ylabel(ax,'Amplitude','FontSize',10);
xlabel(ax,'Time (s)','FontSize',10);

%%%LFE templates, bandpassed
ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on'); 
for i = 1: nsta
  p(i)=plot(ax,(1:lwlet)/sps, greenf(1:lwlet,i), '-','Color',color(i,:),'linew',1); 
end
xlim(ax,[0 lwlet/sps]);  
% ylim(ax,yran/2);
ylim(ax,[-0.5 0.5]);
legend(ax,p,{'PGC','SSIB','SILB','KLNB'},'Location','southeast','NumColumns',3,...
  'fontsize',8);
text(ax,0.02,0.1,'1.8-18 Hz','Units','normalized','HorizontalAlignment','left',...
  'FontSize',10);
% text(ax,0.98,0.9,'LFE templates at family 002','Units','normalized',...
%   'HorizontalAlignment','right','FontSize',10);
text(ax,0.02,0.90,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1);
% shrink(ax,range(xzoom)*sps/lwlet,1);
% ax.Position(1)=axpos(1,1);
xticks(ax,0: 2: lwlet/sps);
longticks(ax,2); 
nolabels(ax,2);
xlabel(ax,'Time (s)','FontSize',10);

if saveflag
  fname = 'lfetemp002.pdf';
  % print(f.fig,'-dpdf',...
  %   strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
  print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));
end