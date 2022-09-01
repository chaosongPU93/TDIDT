function [f,barsw1,barne1,barsw2,barne2,lgd1,lgd2] = ...
    plt_hist_combined_RTMs_v2(indne,indsw,ranvechf98,medallhf,angbest,xran,yran,binw,vdistsw1,wtsw1,...
                           vdistne1,wtne1,vdistsw2,wtsw2,vdistne2,wtne2,conf)
% this function is to simplify the plotting of RTM coverage, LF-HF offset distribution
% of combined RTMs and its corresponding confidence interval and its zoom-in, for two
% opposite propagation groups at one region,
%
% different from 'plt_hist_combined_RTMs', this v2 deal with the specific case that we
% use a line to divide the study region into 2 spatial portions, so that each portion
% would contain data points (LF-HF) from the opposite propagating RTMs
%
% 
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/11/22
% Last modified date:   2020/11/22 
[scrsz, res] = pixelperinch(1);

f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 3;
for isub = 1:nrow*ncol-1
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
xwid = 0.27;
set(f.ax(1), 'position', [ 0.05, 0.08, xwid, 0.85]);

% subplot 1
hold(f.ax(1),'on');
ax = f.ax(1);
ax.FontSize = 8;
% xran = [-15 5];
% yran = [-15 15];
% indplt = union(indsw,indne);    % index to plot
% indplt = indne;    % index to plot
[ax,c] = plt_mig_prop_onmap(ax, indne, indsw, ranvechf98, medallhf, angbest, xran,yran);
% text(f.ax(1),0.5,0.06,'All detections','FontSize',12,'unit','normalized',...
%      'horizontalalignment','center');
% text(f.ax(1),0.5,0.12,'Main LZB region,','FontSize',12,'unit','normalized',...
%      'horizontalalignment','center'); 
text(ax,0.04,0.95,'a','FontSize',10,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.95,0.95,'LZB','FontSize',10,'unit','normalized','EdgeColor','k','Margin',2,...
    'horizontalalignment','right');
hold(ax,'off');
f.ax(1).Position(2) = f.ax(1).Position(2)-0.05;
% c.Position(2) = c.Position(2)-0.1;

% subplot 2 
% binw = 1;

ywid = f.ax(1).Position(4)+f.ax(1).Position(2)-(c.Position(2)) - 0.05;
% ywid = f.ax(1).Position(4) + 0.092 - 0.05;
% ywid = f.ax(1).Position(4)-f.ax(1).Position(2)+(c.Position(2)) - 0.06;
set(f.ax(2), 'position', [ 0.40, c.Position(2), xwid, ywid*2/3 ]);

meansw = wt_mean(vdistsw1,wtsw1);
meanne = wt_mean(vdistne1,wtne1);
sigmasw = sqrt(wt_var(vdistsw1,wtsw1,2));
neff = sum(wtsw1);
% conf = 99;
CIsw = confidence_interval(meansw,sigmasw,neff,conf);
sigmane = sqrt(wt_var(vdistne1,wtne1,2));
neff = sum(wtne1);
% conf = 99;
CIne = confidence_interval(meanne,sigmane,neff,conf);

hold(f.ax(2),'on');
ax = f.ax(2);
ax.FontSize = 8;
[ax, barsw1, ~, ~, ~, ~] = plt_weighted_dist(ax, vdistsw1, wtsw1, binw,'dec');
barsw1(1).FaceAlpha = 1;
[ax, barne1, ~, ~, ~, ~] = plt_weighted_dist(ax, vdistne1, wtne1, binw,'dec');
barne1(1).FaceColor = [0 1 1];
plot(ax,[meansw meansw],[-100 100],'r--','linew',1.5);    % median of wt dist
plot(ax,[meanne meanne],[-100 100],'b--','linew',1.5);    % median of wt dist
errorbar(ax,meansw,0.26,CIsw(1)-meansw,CIsw(2)-meansw,'horizontal','o','markersize',3,'color',...
         'r','linew',1,'MarkerEdgeColor','r','MarkerFaceColor','r');
errorbar(ax,meanne,0.26,CIne(1)-meanne,CIne(2)-meanne,'horizontal','o','markersize',3,'color',...
         'b','linew',1,'MarkerEdgeColor','b','MarkerFaceColor','b');
ax.XLim = [-12.5 12.5];
xx = ax.XLim;
xvect = [-7 xx(1)+0.5];
yvect = [0.015 0.015];
drawArrow(ax,xvect,yvect,ax.XLim,ax.YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
xvect = [7 xx(2)-0.5];
yvect = [0.015 0.015];
drawArrow(ax,xvect,yvect,ax.XLim,ax.YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
text(ax,0.2,0.6,num2str(length(vdistne1)),'fontsize',10,'unit','normalized',...
     'horizontalalignment','center','color','b');
text(ax,0.2,0.54,strcat({'detections'}),'fontsize',10,'unit','normalized',...
     'horizontalalignment','center','color','k');
text(ax,0.85,0.6,num2str(length(vdistsw1)),'fontsize',10,'unit','normalized',...
     'horizontalalignment','center','color','r');
text(ax,0.85,0.54,strcat({'detections'}),'fontsize',10,'unit','normalized',...
     'horizontalalignment','center','color','k');
% text(f.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% lgd1 = legend(ax,[barne1,barsw1],{'ENE propagation','WSW propagation'},'fontsize',7,...
%        'numcolumns',2,...
%        'Position',[0.55  c.Position(2)+ywid*2/3*0.93  xwid*0.9  ywid*2/3*0.05]);
lgd1 = [];
hold(ax,'off');

% subplot 3
set(f.ax(3), 'position', [ 0.40, c.Position(2)+ywid*2/3+0.05, xwid, ywid*1/3 ]);
ax = f.ax(3);
ax.FontSize = 8;
hold(ax,'on');
plot(ax,[0 0],[-100 100],'--','color',[0.6 0.6 0.6],'linew',1);
plot(ax,[meansw meansw],[-100 100],'r--','linew',2);    % median of wt dist
plot(ax,[meanne meanne],[-100 100],'b--','linew',2);    % median of wt dist
errorbar(ax,meansw,0.26,CIsw(1)-meansw,CIsw(2)-meansw,'horizontal','o','markersize',5,'color',...
         'r','linewidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','CapSize',10);
errorbar(ax,meanne,0.26,CIne(1)-meanne,CIne(2)-meanne,'horizontal','o','markersize',5,'color',...
         'b','linewidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','CapSize',10);
% text(ax,0.95,0.86,strcat({'10X'}),'fontsize',10,'unit','normalized','EdgeColor','k',...
%      'Margin',2,'horizontalalignment','right');
% text(ax,0.05,0.15,strcat({'99% CI'}),'fontsize',12,'unit','normalized','horizontalalignment','left'); 
text(ax,0.04,0.86,'b','FontSize',10,'unit','normalized','EdgeColor','k','Margin',2);
ylim(ax,[0.2 0.3]);
xticks(ax,-3: 0.5: 3);
yticks(ax,[0.2 0.25 0.3]);
ylabel(ax,'PDF estimate','fontsize',11);
ax.Box = 'on';
grid(ax, 'on');
ax.GridLineStyle = '--';
hold(ax,'off');

xlimit = f.ax(2).XLim/5;
xlim(f.ax(3),xlimit); 
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(3),xlimit(1),0.2);
[xe,ye] = ds2nfu(f.ax(2),xlimit(1),0.3);
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(3),xlimit(2),0.2);
[xe,ye] = ds2nfu(f.ax(2),xlimit(2),0.3);
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');

% subplot 4
set(f.ax(4), 'position', [ 0.7, c.Position(2), xwid, ywid*2/3 ]);

meansw = wt_mean(vdistsw2,wtsw2);
sigmasw = sqrt(wt_var(vdistsw2,wtsw2,2));
neff = sum(wtsw2);
% conf = 99;
CIsw = confidence_interval(meansw,sigmasw,neff,conf);
meanne = wt_mean(vdistne2,wtne2);
sigmane = sqrt(wt_var(vdistne2,wtne2,2));
neff = sum(wtne2);
% conf = 99;
CIne = confidence_interval(meanne,sigmane,neff,conf);
hold(f.ax(4),'on');
ax = f.ax(4);
ax.FontSize = 8;
[ax, barsw2, ~, ~, ~, ~] = plt_weighted_dist(ax, vdistsw2, wtsw2, binw,'dec');
barsw2(1).FaceAlpha = 1;
[ax, barne2, ~, ~, ~, ~] = plt_weighted_dist(ax, vdistne2, wtne2, binw,'dec');
barne2(1).FaceColor = [0 1 1];
plot(ax,[meansw meansw],[-100 100],'r--','linew',1.5);    % median of wt dist
plot(ax,[meanne meanne],[-100 100],'b--','linew',1.5);    % median of wt dist
errorbar(ax,meansw,0.26,CIsw(1)-meansw,CIsw(2)-meansw,'horizontal','o','markersize',3,'color',...
         'r','linew',1,'MarkerEdgeColor','r','MarkerFaceColor','r');
errorbar(ax,meanne,0.26,CIne(1)-meanne,CIne(2)-meanne,'horizontal','o','markersize',3,'color',...
         'b','linew',1,'MarkerEdgeColor','b','MarkerFaceColor','b');
ax.XLim = [-12.5 12.5];
xx = ax.XLim;
xvect = [-7 xx(1)+0.5];
yvect = [0.015 0.015];
drawArrow(ax,xvect,yvect,ax.XLim,ax.YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
xvect = [7 xx(2)-0.5];
yvect = [0.015 0.015];
drawArrow(ax,xvect,yvect,ax.XLim,ax.YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
text(ax,0.2,0.6,num2str(length(vdistne2)),'fontsize',10,'unit','normalized',...
     'horizontalalignment','center','color','b');
text(ax,0.2,0.54,strcat({'detections'}),'fontsize',10,'unit','normalized',...
     'horizontalalignment','center','color','k');
text(ax,0.85,0.6,num2str(length(vdistsw2)),'fontsize',10,'unit','normalized',...
     'horizontalalignment','center','color','r');
text(ax,0.85,0.54,strcat({'detections'}),'fontsize',10,'unit','normalized',...
     'horizontalalignment','center','color','k');
% text(f.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
lgd2 = legend(ax,[barne2,barsw2],{'ENE propagation','WSW propagation'},'fontsize',7,...
       'numcolumns',2,...
       'Position',[0.570  c.Position(2)+ywid*2/3*0.915  xwid*0.9  ywid*2/3*0.05]);
% lgd2=[];
ax.YLabel.String = [];
ax.YAxisLocation = 'right';
hold(ax,'off');

set(f.ax(5), 'position', [ 0.7, c.Position(2)+ywid*2/3+0.05, xwid, ywid*1/3 ]);
ax = f.ax(5);
ax.FontSize = 8;
hold(ax,'on');
plot(ax,[0 0],[-100 100],'--','color',[0.6 0.6 0.6],'linew',1);
plot(ax,[meansw meansw],[-100 100],'r--','linew',2);    % median of wt dist
plot(ax,[meanne meanne],[-100 100],'b--','linew',2);    % median of wt dist
errorbar(ax,meansw,0.26,CIsw(1)-meansw,CIsw(2)-meansw,'horizontal','o','markersize',5,'color',...
         'r','linewidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','CapSize',10);
errorbar(ax,meanne,0.26,CIne(1)-meanne,CIne(2)-meanne,'horizontal','o','markersize',5,'color',...
         'b','linewidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','CapSize',10);
% text(ax,0.95,0.86,strcat({'10X'}),'fontsize',10,'unit','normalized','EdgeColor','k',...
%      'Margin',2,'horizontalalignment','right');
% text(ax,0.05,0.15,strcat({'99% CI'}),'fontsize',12,'unit','normalized','horizontalalignment','left'); 
text(ax,0.04,0.86,'c','FontSize',10,'unit','normalized','EdgeColor','k','Margin',2);
ylim(ax,[0.2 0.3]);
xticks(ax,-3: 0.5: 3);
yticks(ax,[0.2 0.25 0.3]);
% ylabel(ax,'PDF estimate','fontsize',11);
ax.Box = 'on';
grid(ax, 'on');
ax.GridLineStyle = '--';
ax.YAxisLocation = 'right';
hold(ax,'off');

xlimit = f.ax(4).XLim/5;
xlim(f.ax(5),xlimit); 
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(5),xlimit(1),0.2);
[xe,ye] = ds2nfu(f.ax(4),xlimit(1),0.3);
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(5),xlimit(2),0.2);
[xe,ye] = ds2nfu(f.ax(4),xlimit(2),0.3);
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');

% keyboard