function f=plt_shortwin_summary_noise(sigsta,imp,timevec,locxy,imp4th,timevec4th,...
  locxy4th,tstbuf,dy,mo,yr,sps,ttype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_shortwin_summary_noise(sigsta,imp,timevec,locxy,imp4th,timevec4th,...
%   locxy4th,tstbuf,dy,mo,yr,sps,ttype)
%
% This function is to plot a summary figure for short-win detections of SYNTHETIC
% NOISE.vIt contains the noise seismogram, map view of sources after 2ndary 
% removed, and map view of sources checked at KLNB, projection of sources checked
% at KLNB. 
% --For NOISE detection, this function is expected to be called for each burst
% to create a summary plot for each burst, then all summary figures will be
% merged into a single pdf file in the driving script, 
% 'deconv_ref_4s_exp_4thsta_fn.m'. 
% --Note that this function do not plot the projection along the propagation
% direction as there should be no expected migration from noise. Something might
% be there because you force individual alignment for each short window.
% --Most of the code are similar and consistent with
% 'plt_shortwin_summary' for short-win detections of DATA.
% --Unlike 'plt_shortwin_summary' that has the option 'purpose', there is NO
% need to plot the projection, so NO 'purpose' here.
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
htin = 4.5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 3;
f = initfig(widin,htin,nrow,ncol); %initialize fig

%%% seismograms
set(f.ax(1), 'position', [0.08 0.79 0.85 0.18]);
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% yyaxis(ax,'right');
% % plot(ax,ircccat/sps,rcccat,'o','color',[.7 .7 .7],'markersize',1);
% plot(ax,ircccat/sps,rcccat,'color',[.7 .7 .7]);
% ylim(ax,[-1 1]);
% yyaxis(ax,'left');
ym = max(abs(sigsta(:)));
yran=1.5*[-ym ym];
% yran=[-1.325 1.325];
lsig = size(sigsta,1); 
nsta = size(sigsta,2); 
for i = 1: nsta
  p(i)=plot(ax,(1:lsig)/sps, sigsta(:,i), '-','Color',color(i,:));
end
xlim(ax,[0,lsig/sps]);  ylim(ax,yran); 
%plot the deconvolved sources
yloc = (yran(1)+range(yran)*0.05);
if ~isempty(imp)
  scatter(ax,imp(:,1)/sps,yloc*ones(size(imp(:,1))),4,[.7 .7 .7],'filled');
end
% text(ax,0.5,0.9,num2str(size(impindepst,1)),'fontsize',10,...
%   'Units','normalized');
text(ax,0.99,0.9,'Synthetic noise','Units','normalized','HorizontalAlignment',...
  'right','FontSize',10);
text(ax,0.01,0.85,'a','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1);
legend(ax,p,'PGC','SSIB','SILB','Location','north','Orientation','horizontal',...
  'fontsize',8);
xlabel(ax,sprintf('Time (s) since %.1f s on %s %s %s',tstbuf,dy,mo,yr),...
  'FontSize',10);
ylabel(ax,'Amplitude','FontSize',10);
longticks(ax,5); 
hold(ax,'off');

%%% 2ndary removed, map view
%mannually set the locations for each axis
set(f.ax(2), 'position', [0.08 0.1 0.38 0.6]);
ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
plot(ax,xcut,ycut,'k-','linew',2);
cran = [0 lsig/sps];
xran = [-4 4];
yran = [-4 4];
msize = 35;
if ~isempty(imp)
  wt = median(imp(:,[2 4 6]),2);
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
% text(ax,0.99,0.95,'Secondary removed','HorizontalAlignment','right','Units',...
%   'normalized','FontSize',8);  
text(ax,0.99,0.95,'3-station','HorizontalAlignment','right','Units',...
  'normalized','FontSize',10);  
text(ax,0.02,0.05,'b','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
longticks(ax,2);

%%% 4th-sta checked, map view
%mannually set the locations for each axis
% c1pos=c1.Position;
% ax2pos=f.ax(2).Position;
% set(f.ax(3), 'position', [c1pos(1)+c1pos(3)+0.12 0.1 ax2pos(3) ax2pos(4)]);
set(f.ax(3), 'position', [0.54 0.1 0.38 0.6]);
ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
plot(ax,xcut,ycut,'k-','linew',2);
if ~isempty(imp4th) && ~isempty(imp)
  wt4th = median(imp4th(:,[2 4 6]),2);
%   wt4thmax = prctile(wt4th,95); %use percentile in case
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
% text(ax,0.99,0.95,'Checked at KLNB','HorizontalAlignment','right','Units',...
%   'normalized','FontSize',8);  
text(ax,0.99,0.95,'4-station','HorizontalAlignment','right','Units',...
  'normalized','FontSize',10);  
text(ax,0.02,0.05,'c','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
longticks(ax,2);

% keyboard

