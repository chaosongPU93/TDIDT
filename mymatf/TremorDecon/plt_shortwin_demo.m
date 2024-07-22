function f=plt_shortwin_demo(imp,timevec,locxy,locxyproj,stats,...
  imp4th,timevec4th,locxy4th,locxyproj4th,stats4th,lsig,sps,ttype,saveflag,hfall)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_shortwin_demo(imp,timevec,locxy,locxyproj,stats,...
%   imp4th,timevec4th,locxy4th,locxyproj4th,stats4th,lsig,sps,ttype)
%
% Similar to 'plt_wholewin_demo.m', this script would plot the detections
% BOTH after 2ndary removed AND after 4th-sta checked, for DATA detection, but
% without the seismogram, since in the paper, that would be duplicate
% --This is to plot a summary figure for one burst window for the purpose of
% demonstration, should be called in 'deconv_ref_4s_exp_4thsta_181demo.m'. 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/03/11
% Last modified date:   2024/04/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('sps',160);
defval('ttype','tori');  % default is sort the sources by origin time
defval('saveflag',1);  %default is short-win detections
defval('hfall',[]);  %4-s tremor catalog

%the cut-out boundary of 4-s detections
cutout = 'ellipse';
x0 = 0.2; 
y0 = 0.2;
semia = 1.75;
semib = 1.25;
angrot = 45;
[xcut, ycut] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);

widin = 8.3;  % maximum width allowed is 8.5 inches
htin = 7;   % maximum height allowed is 11 inches
nrow = 2;
ncol = 2;
f = initfig(widin,htin,nrow,ncol); %initialize fig

[~,indremove] = setdiff(locxy,locxy4th,'rows','stable');

%%% 2ndary removed, map view
%mannually set the locations for each axis
set(f.ax(1), 'position', [0.08 0.55 0.42 0.4]);
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
plot(ax,xcut,ycut,'k-','linew',2);
cran = [0 lsig/sps];
xran = [-4 4];
yran = [-4 4];
msize = 40;
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
text(ax,0.02,0.9,strcat({'Amplitude '},'$\geq$',{' 95th prctile'}) ,'Units','normalized',...
  'HorizontalAlignment','left','FontSize',8,'interpreter','latex');
axis(ax,'equal');
xticks(ax,xran(1):1:xran(2));
yticks(ax,yran(1):1:yran(2));
axis(ax,[xran yran]);
if ~isempty(stats)
  angrmse = stats.angrmse;
  [rotx, roty] = complex_rot(0,1,-angrmse);
  xarrow = [0.5-rotx 0.5+rotx];
  yarrow = [-2.5-roty -2.5+roty];
  % [xarrown,yarrown] = ds2nfu(ax,xarrow,yarrow);
  % p=annotation('arrow',xarrown,yarrown,'color','k','linestyle','-','linewidth',1.5);
  p=annotation('arrow','color','k','linestyle','-','linewidth',1.5);
  p.Parent = ax;
  p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
  % text(ax,0.6,0.2,strcat(num2str(angrmse),'$^{\,\circ}$'),'FontSize',11,...
  %     'unit','normalized','interpreter','latex','HorizontalAlignment','left');
  text(ax,0.6,0.2,sprintf('%d%c',angrmse,char(176)),...
    'Units','normalized','HorizontalAlignment','left','FontSize',10);
end
% text(ax,0.99,0.95,'Secondary removed','HorizontalAlignment','right','Units',...
%   'normalized','FontSize',9);  
text(ax,0.99,0.95,'3-station','HorizontalAlignment','right','Units',...
  'normalized','FontSize',10);  
text(ax,0.02,0.05,'a','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
longticks(ax,2);

%%% 2ndary removed, projection
%mannually set the locations for each axis
c1pos=c1.Position;
ax1pos=f.ax(1).Position;
set(f.ax(2), 'position', [c1pos(1)+c1pos(3)+0.12 c1pos(2) ax1pos(3) c1pos(4)]);
ax=f.ax(2);
hold(ax,'on');
ax.Box='on'; grid(ax,'on');
if ~isempty(locxyproj)&& ~isempty(stats)
  %   scatter(ax,timevec/sps,locxyproj(:,2),15,[.6 1 1],'filled','o','MarkerEdgeColor','k');
  % scatter(ax,timevec/sps,locxyproj(:,1),30*refscl,[.5 .5 .5],'filled','o',...
  %   'MarkerEdgeColor','k');
  p1=scatter(ax,timevec/sps,locxyproj(:,1),msize*refscl,[.4 .4 .4],'filled','o',...
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
text(ax,0.02,0.05,'b','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);

%%%%%%%%%%%%%%% projections of loc of tremors within the same win
ttol = 35;
tranbstbuf = load(strcat(rstpath, '/MAPS/tdec.bstranbuf',...
  num2str(ttol),'s.pgc002.',cutout(1:4)));
ttol = 1e-3*86400;
tranmig = load(strcat(rstpath, '/MAPS/migran',num2str(round(ttol)),'s.pgc002'),'w+');
indmig = 109;
ind = find(hfall(:,daycol)==tranmig(indmig,1) & ...
    hfall(:,seccol)>=tranmig(indmig,2) & ...
    hfall(:,seccol)<=tranmig(indmig,3));
mig = hfall(ind,:);
mig(:,seccol) = mig(:,seccol)-tranmig(indmig,2);
indtmr = find(mig(:,seccol)+tranmig(indmig,2)>=tranbstbuf(181,2) & ...
  mig(:,seccol)+tranmig(indmig,2)<=tranbstbuf(181,3));
migtmr = sortrows(mig(indtmr,:),seccol);
timevecbst = migtmr(:,seccol)+tranmig(indmig,2)-tranbstbuf(181,2);
migtmrdum = migtmr;
for j = 1: size(migtmr,1)
  x0 = migtmrdum(j,1);
  y0 = migtmrdum(j,2);
  [newx,newy] = coordinate_rot(x0,y0,-(stats.angrmse-90),[0 0]);
  migtmrdum(j,1) = newx;
  migtmrdum(j,2) = newy;
end
p2=scatter(ax,timevecbst,migtmrdum(:,1),msize*3/4,'rs','linewidth',0.75);
% create fit object  
[fitobjtmr,goftmr,outtmr] = fit(timevecbst,migtmrdum(:,1),fttpfree,'Robust',...
  'Bisquare','StartPoint',[1 1]);
%%%Some statistics
statstmr = statsofrobustlnfit(fitobjtmr,goftmr,outtmr,timevecbst,migtmrdum(:,1));
fittmr = feval(fitobjtmr,timevecbst);
plot(ax,timevecbst,fittmr,'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
text(ax,0.99,0.22,sprintf('Speed: %.1f m/s',statstmr.slope*1e3),...
  'HorizontalAlignment','right','Units','normalized','FontSize',8,'color','r');
text(ax,0.99,0.17,sprintf('Pearson: %.2f',statstmr.pearwt),...
  'HorizontalAlignment','right','Units','normalized','FontSize',8,'color','r');
%%%%%%%%%%%%%%% projections of loc of tremors within the same win

lgd=legend(ax,[p1,p2],'LFE', '4-s tremor','Location','northwest');
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));
if isequal(ttype,'tori')
  xlabel(ax,'Relative origin time (s)');
elseif isequal(ttype,'tarvl')
  xlabel(ax,'Arrival time (s)');
end
ylabel(ax,'Projected location (km)');
% ym = ceil(max(abs(locxyproj(:,1))));
% ym = 3;  % a value that is suitable for all bursts
ylim(ax,[-5 2]);
yticks(ax,-5: 1: 2);
xlim(ax,cran);
longticks(ax,2);
% keyboard

%%% 4th-sta checked, map view
%mannually set the locations for each axis
set(f.ax(3), 'position', [0.08 0.08 0.42 0.4]);
ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
plot(ax,xcut,ycut,'k-','linew',2);
if ~isempty(imp4th) && ~isempty(imp)
  wt4th = median(imp4th(:,[2 4 6]),2);
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
  % text(ax,0.6,0.2,strcat(num2str(angrmse4th),'$^{\,\circ}$'),'FontSize',11,...
  %     'unit','normalized','interpreter','latex','HorizontalAlignment','left');
  text(ax,0.6,0.2,sprintf('%d%c',angrmse4th,char(176)),...
    'Units','normalized','HorizontalAlignment','left','FontSize',10);
end
% text(ax,0.99,0.95,'Checked at KLNB','HorizontalAlignment','right','Units',...
%   'normalized','FontSize',8);  
text(ax,0.99,0.95,'4-station','HorizontalAlignment','right','Units',...
  'normalized','FontSize',10);  
text(ax,0.02,0.05,'c','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
longticks(ax,2);

%%% 4th-sta checked, PROJECTION
%mannually set the locations for each axis
c2pos=c2.Position;
ax3pos=f.ax(3).Position;
set(f.ax(4), 'position', [c2pos(1)+c2pos(3)+0.12 c2pos(2) ax3pos(3) c2pos(4)]);
ax=f.ax(4);
hold(ax,'on');
ax.Box='on'; grid(ax,'on');
if ~isempty(locxyproj4th) && ~isempty(stats4th) && ~isempty(locxy) && ~isempty(locxyproj)
  %   scatter(ax,timevec/sps,locxyproj(:,2),15,[.6 1 1],'filled','o','MarkerEdgeColor','k');
  % scatter(ax,timevec/sps,locxyproj(:,1),30*refscl,[.5 .5 .5],'filled','o',...
  %   'MarkerEdgeColor','k');
  p1=scatter(ax,timevec4th/sps,locxyproj4th(:,1),msize*refscl4th,[.4 .4 .4],'filled','o',...
    'MarkerEdgeColor','w');
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
  if length(indremove) >=2
    %%%linear fit to removed srcs along the same prop direc as 4th srcs
    [fitobj2,gof2,out2] = fit(timevec(indremove)/sps,locxyprojrem(:,1),fttpfree,...
      'Robust','Bisquare','StartPoint',[1 1]);
    %%%Some statistics
    statsrem = statsofrobustlnfit(fitobj2,gof2,out2,timevec(indremove)/sps,...
      locxyprojrem(:,1));
    % output fit parameters
    fitproj2 = feval(fitobj2,timevec(indremove)/sps);
    plot(ax,timevec(indremove)/sps,fitproj2,'--','linewidth',1.5,'color','r');
    text(ax,0.99,0.22,sprintf('Speed: %.1f m/s',statsrem.slope*1e3),...
      'HorizontalAlignment','right','Units','normalized','FontSize',8,'color','r');
    text(ax,0.99,0.17,sprintf('Pearson: %.2f',statsrem.pearwt),...
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
text(ax,0.02,0.05,'d','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
if isequal(ttype,'tori')
  xlabel(ax,'Relative origin time (s)');
elseif isequal(ttype,'tarvl')
  xlabel(ax,'Arrival time (s)');
end
ylabel(ax,'Projected location (km)');
ylim(ax,[-5 2]);
yticks(ax,-5: 1: 2);
xlim(ax,cran);
longticks(ax,2);
hold(ax,'off');

if saveflag
  fname = strcat('shortwinsum_demo.pdf');
  print(f.fig,'-dpdf',...
    strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
end

% keyboard
