function f=plt_wholewin_mapandproj(timevec,locxy,locxyproj,stats,wt,lsig,sps,...
  ttype,str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_wholewin_mapandproj(timevec,locxy,locxyproj,stats,wt,lsig,sps,ttype,str)
%
% For whole window detections, we don't care much details, at least not in
% terms of various plots. A simple plot here is the map-view of deconvolved
% sources scaled by the amplitude and colorcoded by origin time, annotated
% with the propagation direction. Besides that, the projection of sources 
% along the propagation direction with some fitting statistics is also shown.
% Could be viewed as an complement to 'plt_deconlfit.m'. 
% --This script can plot detections EITHER after 2ndary removed OR after 
% 4th-sta checked.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/03/11
% Last modified date:   2024/03/11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('ftrans','interpchao');  % default is sort the sources by origin time
% defval('str','Secondary removed');  % default is sort the sources by origin time
defval('str','3-station');  % default is sort the sources by origin time

%the cut-out boundary of 4-s detections
cutout = 'ellipse';
x0 = 0.2; 
y0 = 0.2;
semia = 1.75;
semib = 1.25;
angrot = 45;
[xcut, ycut] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);

widin = 8.3;  % maximum width allowed is 8.5 inches
htin = 4;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f = initfig(widin,htin,nrow,ncol); %initialize fig

%%% 2ndary removed, map view
%mannually set the locations for each axis
set(f.ax(1), 'position', [0.08 0.1 0.42 0.8]);
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
plot(ax,xcut,ycut,'k-','linew',2);
cran = [0 lsig/sps];
xran = [-4 4];
yran = [-4 4];
msize = 40;
if ~isempty(locxy)
  wtmax = prctile(wt,95); %use percentile in case
  refscl = wt./wtmax;
  refscl(refscl>=1) = 1;  %force the larger amp to be plotted as the same size in case of saturation
  scatter(ax,locxy(:,1),locxy(:,2),msize*refscl,timevec/sps,'filled','o',...
    'MarkerEdgeColor',[.5 .5 .5]);
end
text(ax,0.99,0.05,sprintf('%d events',size(locxy,1)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
% colormap(ax,flipud(colormap(ax,'kelicol')));
colormap(ax,'viridis');
% colormap(ax,'plasma');
c=colorbar(ax);
caxis(ax,cran);
if isequal(ttype,'tori')
  c.Label.String = sprintf('Relative origin time (s)');
elseif isequal(ttype,'tarvl')
  c.Label.String = sprintf('Arrival time (s)');
end
c.Label.FontSize = 9;
scatter(ax,xran(1)+0.1*range(xran),yran(2)-0.05*range(yran),msize,'w','filled',...
  'MarkerEdgeColor',[.5 .5 .5]);
text(ax,0.02,0.9,'Amplitude\geq 95th prctile','Units','normalized',...
  'HorizontalAlignment','left','FontSize',8);
axis(ax,'equal');
xticks(ax,xran(1):1:xran(2));
yticks(ax,yran(1):1:yran(2));
axis(ax,[xran yran]);
if ~isempty(locxy)
  angrmse = stats.angrmse;
  [rotx, roty] = complex_rot(0,1,-angrmse);
  xarrow = [0.5-rotx 0.5+rotx];
  yarrow = [-2.5-roty -2.5+roty];
  % [xarrown,yarrown] = ds2nfu(ax,xarrow,yarrow);
  % p=annotation('arrow',xarrown,yarrown,'color','k','linestyle','-','linewidth',1.5);
  p=annotation('arrow','color','k','linestyle','-','linewidth',1.5);
  p.Parent = ax;
  p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
  % text(ax,0.6,0.2,strcat(num2str(angrmse),'$^{\circ}$'),'FontSize',11,...
  %     'unit','normalized','interpreter','latex','HorizontalAlignment','left');
  text(ax,0.6,0.2,sprintf('%d%c',angrmse,char(176)),...
    'Units','normalized','HorizontalAlignment','left','FontSize',10);
end
text(ax,0.99,0.95,str,'HorizontalAlignment','right','Units',...
  'normalized','FontSize',9);  
text(ax,0.02,0.05,'a','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
longticks(ax,2);

%%% 2ndary removed, projection
%mannually set the locations for each axis
cpos=c.Position;
ax1pos=f.ax(1).Position;
set(f.ax(2), 'position', [cpos(1)+cpos(3)+0.12 cpos(2) ax1pos(3) cpos(4)]);
ax=f.ax(2);
hold(ax,'on');
ax.Box='on'; grid(ax,'on');
if ~isempty(locxy)
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
  %   'HorizontalAlignment','right','Units','normalized','FontSize',9);
  text(ax,0.99,0.10,sprintf('Speed: %.1f m/s',stats.slope*1e3),...
    'HorizontalAlignment','right','Units','normalized','FontSize',9);
  text(ax,0.99,0.05,sprintf('Pearson: %.2f',stats.pearwt),...
    'HorizontalAlignment','right','Units','normalized','FontSize',9);
end
text(ax,0.99,0.95,str,'HorizontalAlignment','right','Units',...
  'normalized','FontSize',9);  
text(ax,0.02,0.05,'b','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
% legend(ax,[p1,p2],'Along prop. direction', 'Orthogonal to prop. dir.','FontSize',9);
if isequal(ttype,'tori')
  xlabel(ax,'Relative origin time (s)');
elseif isequal(ttype,'tarvl')
  xlabel(ax,'Arrival time (s)');
end
ylabel(ax,'Projected location (km)');
ym = ceil(max(abs(locxyproj(:,1))));
ylim(ax,[-ym ym]);
xlim(ax,cran);
longticks(ax,2);
% keyboard


