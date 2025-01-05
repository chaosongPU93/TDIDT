function f=plt_shortwin_summary(sigsta,imp,timevec,locxy,locxyproj,stats,...
  imp4th,timevec4th,locxy4th,locxyproj4th,stats4th,tstbuf,dy,mo,yr,sps,ttype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_shortwin_summary(timevec,locxy,locxyproj,lsig,stats,sps,...
%   ttype,wt,timevec4th,locxy4th,locxyproj4th,stats4th,wt4th)
%
% This function is to plot a summary figure for short-win detections of DATA.
% It contains the signal seismogram, map view of sources after 2ndary removed,
% projection of sources after 2ndary removed, map view of sources checked at
% KLNB,  projection of sources checked at KLNB. Sources that fail the 4th-sta
% check is shown in different color/symbol in the projection plot of successful
% ones. 
% --For DATA detection, this function is expected to be called for each burst
% to create a summary plot for each burst, then all summary figures will be
% merged into a single pdf file in the driving script, 
% 'deconv_ref_4s_exp_4thsta_fn.m'. 
% --Parts of the code are similar and consistent with
% 'plt_wholewin_mapandproj4th', and 'plt_wholewin_mapandproj' for the whole-win
% detection.
% --Added the option of 'purpose' so that is called as 'summary', you will plot
% the max. slope direction too
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/03/18
% Last modified date:   2024/03/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('sps',160);
defval('ttype','tori');  % default is sort the sources by origin time

%the cut-out boundary of 4-s detections
cutout = 'ellipse';
x0 = 0.2; 
y0 = 0.2;
semia = 1.75;
semib = 1.25;
angrot = 45;
[xcut, ycut] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);

color = ['r';'b';'k'];

widin = 8.3;  % maximum width allowed is 8.5 inches
htin = 8.5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 5;
f = initfig(widin,htin,nrow,ncol); %initialize fig

[~,indremove] = setdiff(locxy,locxy4th,'rows','stable');

%%% seismograms
set(f.ax(1), 'position', [0.08 0.88 0.85 0.1]);
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% yyaxis(ax,'right');
% % plot(ax,ircccat/sps,rcccat,'o','color',[.7 .7 .7],'markersize',1);
% plot(ax,ircccat/sps,rcccat,'color',[.7 .7 .7]);
% ylim(ax,[-1 1]);
% yyaxis(ax,'left');
ym = max(abs(sigsta(:)));
yran=1.5*[-ym ym];
% %patch a gray shaded part for zoomed-in in other plots
% xzoom = [5 30];
% patarea = [xzoom(1) yran(2);
%            xzoom(2) yran(2);
%            xzoom(2) yran(1);
%            xzoom(1) yran(1);
%            xzoom(1) yran(2)];
% patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.15,'edgecolor','none');
lsig = size(sigsta,1); 
nsta = size(sigsta,2); 
for i = 1: nsta
  p(i)=plot(ax,(1:lsig)/sps, sigsta(:,i), '-','Color',color(i,:));
end
xlim(ax,[0,lsig/sps]);  ylim(ax,yran); 
%plot the deconvolved sources
yloc = (yran(1)+range(yran)*0.05);
if ~isempty(imp)
  scatter(ax,imp(:,1)/sps,yloc*ones(size(imp(:,1))),4,[.6 .6 .6],'filled');
end
% text(ax,0.5,0.9,num2str(size(impindepst,1)),'fontsize',10,...
%   'Units','normalized');
text(ax,0.99,0.9,'Signal','Units','normalized','HorizontalAlignment','right',...
  'FontSize',10);
text(ax,0.01,0.8,'a','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1);
legend(ax,p,'PGC','SSIB','SILB','Location','north','Orientation','horizontal',...
  'fontsize',8);
xlabel(ax,sprintf('Time (s) since %.1f s on %s %s %s',tstbuf,dy,mo,yr),...
  'FontSize',10);
ylabel(ax,'Amplitude','FontSize',10);
longticks(ax,5); 
hold(ax,'off');

%%% 2ndary removed, MAP VIEW
%mannually set the locations for each axis
set(f.ax(2), 'position', [0.08 0.48 0.42 0.35]);
ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
plot(ax,xcut,ycut,'k-','linew',2);
cran = [0 lsig/sps];
xran = [-4 4];
yran = [-4 4];
msize = 40;
if ~isempty(imp)
  wt = mean(imp(:,[2 4 6]),2);
  wtmax = prctile(wt,95); %use percentile in case
  refscl = wt./wtmax;
  refscl(refscl>=1) = 1;  %force the larger amp to be plotted as the same size in case of saturation
  scatter(ax,locxy(:,1),locxy(:,2),msize*refscl,timevec/sps,'filled','o',...
    'MarkerEdgeColor',[.5 .5 .5]);
end
text(ax,0.99,0.05,sprintf('%d events',size(locxy,1)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',8);
% colormap(ax,flipud(colormap(ax,'kelicol')));
colormap(ax,'viridis');
c1=colorbar(ax);
caxis(ax,cran);
if isequal(ttype,'tori')
  c1.Label.String = sprintf('Relative origin time (s)');
elseif isequal(ttype,'tarvl')
  c1.Label.String = sprintf('Arrival time (s)');
end
c1.Label.FontSize = 9;
scatter(ax,xran(1)+0.1*range(xran),yran(2)-0.05*range(yran),msize,'w','filled',...
  'MarkerEdgeColor',[.5 .5 .5],'linew',1);
text(ax,0.02,0.9,strcat({'Amplitude '},'$\geq$',{' 95th prctile'}) ,'Units',...
  'normalized','HorizontalAlignment','left','FontSize',8,'interpreter','latex');
axis(ax,'equal');
xticks(ax,xran(1):1:xran(2));
yticks(ax,yran(1):1:yran(2));
axis(ax,[xran yran]);
if ~isempty(stats)
  angrmse = stats.angrmse;
  [rotx, roty] = complex_rot(0,1,-angrmse);
  xarrow = [0.5-rotx 0.5+rotx];
  yarrow = [-2.5-roty -2.5+roty];
  p=annotation('arrow','color','k','linestyle','-','linewidth',1.5);
  p.Parent = ax;
  p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
  % text(ax,0.6,0.2,strcat(num2str(angrmse),'$^{\,\circ}$'),'FontSize',11,...
  %   'unit','normalized','interpreter','latex','HorizontalAlignment','left');
  text(ax,0.6,0.2,sprintf('%d%c',angrmse,char(176)),...
    'Units','normalized','HorizontalAlignment','left','FontSize',10);
%   angslope = stats.angslope;
%   [rotx, roty] = complex_rot(0,1,-angslope);
%   xarrow = [0.5-rotx 0.5+rotx];
%   yarrow = [-2.5-roty -2.5+roty];
%   p1=annotation('arrow','color','k','linestyle','--','linewidth',1.5,'color',[.6 .6 .6]);
%   p1.Parent = ax;
%   p1.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
%   % text(ax,0.6,0.25,strcat(num2str(angslope),'$^{\,\circ}$'),'FontSize',11,...
%   %   'unit','normalized','interpreter','latex','HorizontalAlignment','left',...
%   %   'color',[.6 .6 .6]);
%   text(ax,0.6,0.25,sprintf('%d%c',angslope,char(176)),'color',[.6 .6 .6],...
%     'Units','normalized','HorizontalAlignment','left','FontSize',10);    
end
% text(ax,0.99,0.95,'Secondary removed','HorizontalAlignment','right','Units',...
%   'normalized','FontSize',8);  
text(ax,0.99,0.95,'3-station','HorizontalAlignment','right','Units',...
  'normalized','FontSize',10);  
text(ax,0.02,0.05,'b','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
longticks(ax,2);
hold(ax,'off');

%%% 2ndary removed, PROJECTION
%mannually set the locations for each axis
c1pos=c1.Position;
ax2pos=f.ax(2).Position;
set(f.ax(3), 'position', [c1pos(1)+c1pos(3)+0.12 c1pos(2) ax2pos(3) c1pos(4)]);
ax=f.ax(3);
hold(ax,'on');
ax.Box='on'; grid(ax,'on');
if ~isempty(locxyproj) && ~isempty(stats)
  %   scatter(ax,timevec/sps,locxyproj(:,2),15,[.6 1 1],'filled','o','MarkerEdgeColor','k');
  % scatter(ax,timevec/sps,locxyproj(:,1),30*refscl,[.5 .5 .5],'filled','o',...
  %   'MarkerEdgeColor','k');
  scatter(ax,timevec/sps,locxyproj(:,1),msize*refscl,[.4 .4 .4],'filled','o',...
    'MarkerEdgeColor','w');
  % linear robust least square
  fttpfree = fittype( @(a,b,x) a*x+b);
  fitobjproj = fit(timevec/sps, locxyproj(:,1),fttpfree,'Robust','Bisquare',...
    'StartPoint',[1 1]);
  % output fit parameters
  fitproj = feval(fitobjproj,timevec/sps);
  plot(ax,timevec/sps,fitproj,'-','linewidth',2,'color','k');
  fitproj = reshape(fitproj,[],1);
  % text(ax,0.99,0.10,sprintf('Speed: %.1f m/s',range(fitproj)*1e3/range(timevec/sps)),...
  %   'HorizontalAlignment','right','Units','normalized','FontSize',8);
  text(ax,0.99,0.10,sprintf('Speed: %.1f m/s',stats.slope*1e3),...
    'HorizontalAlignment','right','Units','normalized','FontSize',8);
  text(ax,0.99,0.05,sprintf('Pearson: %.2f',stats.pearwt),...
    'HorizontalAlignment','right','Units','normalized','FontSize',8);
else
  text(ax,0.99,0.10,'No projection can be made',...
    'HorizontalAlignment','right','Units','normalized','FontSize',8);
  text(ax,0.99,0.05,'due to 2 or fewer events',...
    'HorizontalAlignment','right','Units','normalized','FontSize',8);  
end
text(ax,0.99,0.95,'3-station','HorizontalAlignment','right','Units',...
  'normalized','FontSize',10);  
text(ax,0.02,0.05,'c','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
% legend(ax,[p1,p2],'Along prop. direction', 'Orthogonal to prop. dir.','FontSize',8);
if isequal(ttype,'tori')
  xlabel(ax,'Relative origin time (s)');
elseif isequal(ttype,'tarvl')
  xlabel(ax,'Arrival time (s)');
end
ylabel(ax,'Projected location (km)');
% ym = ceil(max(abs(locxyproj(:,1))));
ym = 3;  % a value that is suitable for all bursts
ylim(ax,[-ym ym]);
yticks(ax,-ym: 1: ym);
xlim(ax,cran);
longticks(ax,2);
hold(ax,'off');
% keyboard

%%% 4th-sta checked, MAP VIEW
%mannually set the locations for each axis
set(f.ax(4), 'position', [0.08 0.08 0.42 0.35]);
ax=f.ax(4); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
plot(ax,xcut,ycut,'k-','linew',2);
if ~isempty(imp4th) && ~isempty(imp)
  wt4th = mean(imp4th(:,[2 4 6]),2);
  % wt4thmax = prctile(wt4th,95); %use percentile in case
  refscl4th = wt4th./wtmax;
  refscl4th(refscl4th>=1) = 1;  %force the larger amp to be plotted as the same size in case of saturation
  scatter(ax,locxy4th(:,1),locxy4th(:,2),msize*refscl4th,timevec4th/sps,'filled',...
    'o','MarkerEdgeColor',[.5 .5 .5]);
end
text(ax,0.99,0.05,sprintf('%d events',size(locxy4th,1)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',8);
% colormap(ax,flipud(colormap(ax,'kelicol')));
colormap(ax,'viridis');
c2=colorbar(ax);
caxis(ax,cran);
if isequal(ttype,'tori')
  c2.Label.String = sprintf('Relative origin time (s)');
elseif isequal(ttype,'tarvl')
  c2.Label.String = sprintf('Arrival time (s)');
end
c2.Label.FontSize = 9;
% scatter(ax,xran(1)+0.1*range(xran),yran(2)-0.05*range(yran),msize,'w','filled',...
%   'MarkerEdgeColor',[.5 .5 .5]);
% text(ax,0.02,0.9,'Amplitude\geq 95th prctile','Units','normalized',...
%   'HorizontalAlignment','left','FontSize',8);
axis(ax,'equal');
xticks(ax,xran(1):1:xran(2));
yticks(ax,yran(1):1:yran(2));
axis(ax,[xran yran]);
if ~isempty(stats4th)
  angrmse4th = stats4th.angrmse;
  [rotx, roty] = complex_rot(0,1,-angrmse4th);
  xarrow = [0.5-rotx 0.5+rotx];
  yarrow = [-2.5-roty -2.5+roty];
  % [xarrown,yarrown] = ds2nfu(ax,xarrow,yarrow);
  % p=annotation('arrow',xarrown,yarrown,'color','k','linestyle','-','linewidth',1.5);
  p=annotation('arrow','color','k','linestyle','-','linewidth',1.5);
  p.Parent = ax;
  p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
%   text(ax,0.6,0.2,strcat(num2str(angrmse4th),{'$\,^{\circ}$'}),'FontSize',11,...
%     'unit','normalized','interpreter','latex','HorizontalAlignment','left');
  text(ax,0.6,0.2,sprintf('%d%c',angrmse4th,char(176)),...
    'Units','normalized','HorizontalAlignment','left','FontSize',10);
%   angslope4th = stats4th.angslope;
%   [rotx, roty] = complex_rot(0,1,-angslope4th);
%   xarrow = [0.5-rotx 0.5+rotx];
%   yarrow = [-2.5-roty -2.5+roty];
%   p1=annotation('arrow','color','k','linestyle','--','linewidth',1.5,'color',[.6 .6 .6]);
%   p1.Parent = ax;
%   p1.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
%   % text(ax,0.6,0.25,strcat(num2str(angslope4th),'$^{\,\circ}$'),'FontSize',11,...
%   %   'unit','normalized','interpreter','latex','HorizontalAlignment','left',...
%   %   'color',[.6 .6 .6]);
%   text(ax,0.6,0.25,sprintf('%d%c',angslope4th,char(176)),'color',[.6 .6 .6],...
%     'Units','normalized','HorizontalAlignment','left','FontSize',10);    
end
% text(ax,0.99,0.95,'Checked at KLNB','HorizontalAlignment','right','Units',...
%   'normalized','FontSize',8);  
text(ax,0.99,0.95,'4-station','HorizontalAlignment','right','Units',...
  'normalized','FontSize',10);  
text(ax,0.02,0.05,'d','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
longticks(ax,2);
% %%%add the indication of the reasonable physical source size
% x0 = -3;
% y0 = 3.5;
% Vs = 3;   % S-wave speed
% %       Vprop = 0.8*Vs;   % reasonable rupture propagation speed
% Vprop = Vs;   % max rupture propagation speed
% td = 40/sps;  % estimated from the bin where the most tsep falls in, 32-48-sample bin
% %         td = median(tsep(:));
% srcsz = Vprop*td;
% radi = srcsz/2;
% scatter(ax,x0,y0,2,'ko','filled');
% [x, y] = circle_chao(x0,y0,radi,0.1);
% plot(ax,x,y,'-','Color','k','linew',1);
% text(ax,0.02,0.86,'Max. reasonable','HorizontalAlignment','left',...
%   'Units','normalized','FontSize',8,'interpreter','latex');
% text(ax,0.02,0.81,'physical source size','HorizontalAlignment','left',...
%   'Units','normalized','FontSize',8,'interpreter','latex');
% % text(ax,0.02,0.76,strcat('$V_{s} \cdot t_{d} = 3 \cdot 0.25$',{' km'}),'HorizontalAlignment','left',...
% %   'Units','normalized','FontSize',8,'interpreter','latex');
% text(ax,0.02,0.76,strcat('$V_{s} \cdot t_{d} \approx 750$',{' m'}),...
%   'HorizontalAlignment','left','Units','normalized','FontSize',8,...
%   'interpreter','latex');
% %%%add the indication of the reasonable physical source size
hold(ax,'off');

%%% 4th-sta checked, PROJECTION
%mannually set the locations for each axis
c2pos=c2.Position;
ax4pos=f.ax(4).Position;
set(f.ax(5), 'position', [c2pos(1)+c2pos(3)+0.12 c2pos(2) ax4pos(3) c2pos(4)]);
ax=f.ax(5);
hold(ax,'on');
ax.Box='on'; grid(ax,'on');
if ~isempty(locxyproj4th) && ~isempty(stats4th) && ~isempty(locxy) && ~isempty(locxyproj)
  %   scatter(ax,timevec/sps,locxyproj(:,2),15,[.6 1 1],'filled','o','MarkerEdgeColor','k');
  % scatter(ax,timevec/sps,locxyproj(:,1),30*refscl,[.5 .5 .5],'filled','o',...
  %   'MarkerEdgeColor','k');
  p1=scatter(ax,timevec4th/sps,locxyproj4th(:,1),msize*refscl4th,[.4 .4 .4],...
    'filled','o','MarkerEdgeColor','w');
  %%%project the removed sources along the same direction as the preserved ones
  [~,~,locxyprojrem] = customprojection(locxy(indremove,1:2),angrmse4th);
  p2=scatter(ax,timevec(indremove)/sps,locxyprojrem(:,1),...
    msize*refscl(indremove),'r','linewidth',0.75);
  % linear robust least square
  fttpfree = fittype( @(a,b,x) a*x+b);
  fitobjproj = fit(timevec4th/sps, locxyproj4th(:,1),fttpfree,'Robust','Bisquare',...
    'StartPoint',[1 1]);
  % output fit parameters
  fitproj = feval(fitobjproj,timevec4th/sps);
  plot(ax,timevec4th/sps,fitproj,'-','linewidth',2,'color','k');
  fitproj = reshape(fitproj,[],1);
  % text(ax,0.99,0.10,sprintf('Speed: %.1f m/s',range(fitproj)*1e3/range(timevec4th/sps)),...
  %   'HorizontalAlignment','right','Units','normalized','FontSize',8);
  text(ax,0.99,0.10,sprintf('Speed: %.1f m/s',stats4th.slope*1e3),...
    'HorizontalAlignment','right','Units','normalized','FontSize',8);
  text(ax,0.99,0.05,sprintf('Pearson: %.2f',stats4th.pearwt),...
    'HorizontalAlignment','right','Units','normalized','FontSize',8);
  if length(indremove)>=2
    %%%linear fit to removed srcs along the same prop direc as 4th srcs
    [fitobj2,gof2,out2] = fit(timevec(indremove)/sps,locxyprojrem(:,1),fttpfree,...
      'Robust','Bisquare','StartPoint',[1 1]);
    %%%Some statistics
    statsrem = statsofrobustlnfit(fitobj2,gof2,out2,timevec(indremove)/sps,...
      locxyprojrem(:,1));
    % output fit parameters
    fitproj2 = feval(fitobj2,timevec(indremove)/sps);
    plot(ax,timevec(indremove)/sps,fitproj2,'--','linewidth',1.5,'color','r');
    text(ax,0.99,0.20,sprintf('Speed: %.1f m/s',statsrem.slope*1e3),...
      'HorizontalAlignment','right','Units','normalized','FontSize',8,'color','r');
    text(ax,0.99,0.15,sprintf('Pearson: %.2f',statsrem.pearwt),...
      'HorizontalAlignment','right','Units','normalized','FontSize',8,'color','r');
  end
  lgd=legend(ax,[p1 p2],'Preserved','Discarded','Location','northwest'); %,'Orientation','horizontal'
  set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));
else
  text(ax,0.99,0.10,'No projection can be made',...
    'HorizontalAlignment','right','Units','normalized','FontSize',8);
  text(ax,0.99,0.05,'due to 2 or fewer events',...
    'HorizontalAlignment','right','Units','normalized','FontSize',8);  
end
text(ax,0.99,0.95,'4-station','HorizontalAlignment','right','Units',...
  'normalized','FontSize',10);  
text(ax,0.02,0.05,'e','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
if isequal(ttype,'tori')
  xlabel(ax,'Relative origin time (s)');
elseif isequal(ttype,'tarvl')
  xlabel(ax,'Arrival time (s)');
end
ylabel(ax,'Projected location (km)');
ylim(ax,[-ym ym]);
yticks(ax,-ym: 1: ym);
xlim(ax,cran);
longticks(ax,2);
hold(ax,'off');

% keyboard


