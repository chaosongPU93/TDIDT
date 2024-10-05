function [f] = plt_intertimeandbursts(hfinter03,hfinter04,hfinter05,newt03,newt04,newt05,ymax,ttol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_intertimeandbursts(hfinter03,hfinter04,hfinter05,newt03,newt04,newt05,ymax,ttol)
% This function is to plot the separation in time between itself and its preceding detection,
% i.e., inter-event time, for all 3 ETS episodes, on top of the automatically-grouped tremor
% burst windows in gray.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/02/13
% Last modified date:   2024/02/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

defval('ymax',2e-3);
defval('ttol',[]);

nrow = 3;
ncol = 1;
widin = 10;  % maximum width allowed is 8.5 inches
htin = 5.5;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.08 0.95]; pltyran = [0.1 0.95];
pltxsep = 0.08; pltysep = 0.08; 
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

msizehf = 2;

ax = f.ax(1); hold(ax,'on'); ax.Box = 'on'; 
for j = 1: size(newt03,1)
  patarea = [newt03(j,1) -3e+4;
             newt03(j,2) -3e+4;
             newt03(j,2) 3e+4;
             newt03(j,1) 3e+4;
             newt03(j,1) -3e+4];
  patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.8,'edgecolor','none');
end
scatter(ax, hfinter03(:,1), hfinter03(:,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w')
plot(ax, ax.XLim, [ttol(1) ttol(1)], 'b--');
text(ax,ax.XLim(1)+0.99*(range(ax.XLim)),1.2*ttol(1),sprintf('%d s',...
  round(ttol(1)*86400)),'FontSize',11,'horizontalalignment','right');
text(ax,0.99,0.85,'2003','FontSize',11,'unit','normalized','horizontalalignment',...
  'right','EdgeColor','k','Margin',2);
% ylabel(ax, 'Time (day) to the preceding detection');
xlabel(ax, 'Time (day) since 01 Mar. 2003');
ylim(ax,[0, ymax]);
longticks(ax,4);
hold(ax,'off');

ax = f.ax(2); hold(ax,'on'); ax.Box = 'on'; 
for j = 1: size(newt04,1)
  patarea = [newt04(j,1) -3e+4;
             newt04(j,2) -3e+4;
             newt04(j,2) 3e+4;
             newt04(j,1) 3e+4;
             newt04(j,1) -3e+4];
  patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.8,'edgecolor','none');
end
scatter(ax, hfinter04(:,1), hfinter04(:,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w')
plot(ax, ax.XLim, [ttol(2) ttol(2)], 'b--');
text(ax,ax.XLim(1)+0.99*(range(ax.XLim)),1.2*ttol(2),sprintf('%d s',...
  round(ttol(2)*86400)),'FontSize',11,'horizontalalignment','right');
text(ax,0.99,0.85,'2004','FontSize',11,'unit','normalized','horizontalalignment',...
  'right','EdgeColor','k','Margin',2);
ylabel(ax, 'Time (day) to the preceding detection');
xlabel(ax, 'Time (day) since 12 Jul. 2004');
ylim(ax,[0, ymax]);
longticks(ax,4);
hold(ax,'off');

ax = f.ax(3); hold(ax,'on'); ax.Box = 'on'; 
for j = 1: size(newt05,1)
  patarea = [newt05(j,1) -3e+4;
             newt05(j,2) -3e+4;
             newt05(j,2) 3e+4;
             newt05(j,1) 3e+4;
             newt05(j,1) -3e+4];
  patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.8,'edgecolor','none');
end
scatter(ax, hfinter05(:,1), hfinter05(:,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w')
plot(ax,ax.XLim, [ttol(3) ttol(3)], 'b--');
text(ax,ax.XLim(1)+0.99*(range(ax.XLim)),1.2*ttol(3),sprintf('%d s',...
  round(ttol(3)*86400)),'FontSize',11,'horizontalalignment','right');
text(ax,0.99,0.85,'2005','FontSize',11,'unit','normalized','horizontalalignment',...
  'right','EdgeColor','k','Margin',2);
% ylabel(ax, 'Time (day) to the preceding detection');
xlabel(ax, 'Time (day) since 11 Sep. 2005');
ylim(ax,[0, ymax]);
longticks(ax,4);
hold(ax,'off');

orient(f.fig,'landscape');
fname = '4stremor_intertandbursts.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));



