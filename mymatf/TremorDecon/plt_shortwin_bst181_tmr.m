function f=plt_shortwin_bst181_tmr(imp,timevec,locxy,locxyproj,stats,...
  imp4th,timevec4th,locxy4th,locxyproj4th,stats4th,lsig,sps,ttype,saveflag,hfall)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_shortwin_bst181_tmr(imp,timevec,locxy,locxyproj,stats,...
%   imp4th,timevec4th,locxy4th,locxyproj4th,stats4th,lsig,sps,ttype)
%
% This function is a combination of 'plt_shortwin_demo.m' and
% 'plt_rtmlfit002_demo.m', in order plot the short-win detections of
% LFE burst win #181, and the 4-s tremor detections inside the same
% window, to make an easier comparison between them.
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


widin = 7.2;  % maximum width allowed is 8.5 inches
htin = 9.5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 7;
f = initfig(widin,htin,nrow,ncol); %initialize fig

[~,indremove] = setdiff(locxy,locxy4th,'rows','stable');

%%% 2ndary removed, map view
%mannually set the locations for each axis
set(f.ax(1), 'position', [0.08 0.70 0.42 0.3]);
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
% c1pos(3)=0.15;
set(f.ax(2), 'position', [c1pos(1)+c1pos(3)+0.13 c1pos(2) ax1pos(3) c1pos(4)]);
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
  [fitobjproj,gof,output] = fit(timevec/sps, locxyproj(:,1),fttpfree,'Robust','Bisquare',...
    'StartPoint',[1 1]);
  reslfe=stats.output.residuals;
%   stats = statsofrobustlnfit(fitobjproj,gof,output,timevec/sps,locxyproj(:,1));
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

%%%%%%%%%%%%%%% below, projections of loc of tremors within the same win
%%%%%%% copy from 'plt_rtmlfit002_demo.m'
workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/PGCtrio');

%original burst win ranges
ttol = 35;
ntol = 3;
tranbst = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.',...
  num2str(ntol),'.pgc002.',cutout(1:4)));

%actually-used in decon, ranges with a little buffer
tranbstbuf = load(strcat(rstpath, '/MAPS/tdec.bstranbuf',...
  num2str(ttol),'s.pgc002.',cutout(1:4)));
tlen = tranbstbuf(:,3)-tranbstbuf(:,2);

ttol = 1e-3*86400;
tranmig = load(strcat(rstpath, '/MAPS/migran',num2str(round(ttol)),'s.pgc002'),'w+');
% tranmig(:,2:3) = tranmig(:,2:3)/3600;

%%%corresponding
indplt = [3; 12; 15; 23; 24; 27; 37; 57; 58; 69; 70; 71; 72; 73; 74;
  92; 93; 94; 95; 98; 100; 103; 109];

%%%4-s migration
indmig = 109;
daycol = 14;
seccol = 16;  % 2 choices here, one is the center of detecting window, other is the start of strongest .5 s
ind = find(hfall(:,daycol)==tranmig(indmig,1) & ...
  hfall(:,seccol)>=tranmig(indmig,2) & ...
  hfall(:,seccol)<=tranmig(indmig,3));
mig = hfall(ind,:);
mig(:,seccol) = mig(:,seccol)-tranmig(indmig,2);
indtmr = find(mig(:,seccol)+tranmig(indmig,2)>=tranbstbuf(181,2) & ...
  mig(:,seccol)+tranmig(indmig,2)<=tranbstbuf(181,3));
tlenmig = tranmig(indmig,3)-tranmig(indmig,2);

%%%4-s detections inside the burst win 181
migtmr = sortrows(mig(indtmr,:),seccol);
timevecbst = migtmr(:,seccol)+tranmig(indmig,2)-tranbstbuf(181,2);
%%%force the projection onto prop. direction of LFEs in burst 181
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
restmr=statstmr.output.residuals;
fittmr = feval(fitobjtmr,timevecbst);
plot(ax,timevecbst,fittmr,'--','linewidth',1.5,'color','r');
text(ax,0.99,0.22,sprintf('Speed: %.1f m/s',statstmr.slope*1e3),...
  'HorizontalAlignment','right','Units','normalized','FontSize',8,'color','r');
text(ax,0.99,0.17,sprintf('Pearson: %.2f',statstmr.pearwt),...
  'HorizontalAlignment','right','Units','normalized','FontSize',8,'color','r');
%%%%%%%%%%%%%%% above, projections of loc of tremors within the same win

figure; hold on; 
binw=0.1;
histogram(reslfe,'binw',binw,'normalization','probability');
histogram(restmr,'binw',binw,'normalization','probability');
[mu1,sigma1]=normfit(reslfe);
x=-2.5:binw:2.5; 
yfit1=normpdf(x,mu1,sigma1)*binw;
plot(x,yfit1,'b-','LineWidth',2);
std1=std(reslfe);
text(0.5,0.5,sprintf('%.1f',std1),'Units','normalized','Color','b');

[mu2,sigma2]=normfit(restmr(abs(restmr)<3));
yfit2=normpdf(x,mu2,sigma2)*binw;
plot(x,yfit2,'r-','LineWidth',2);
std2=std(restmr(abs(restmr)<3));
text(0.5,0.4,sprintf('%.1f',std2),'Units','normalized','Color','r');


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
set(f.ax(3), 'position', [0.08 0.395 0.42 0.3]);
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
%   'normalized','FontSize',9);
text(ax,0.99,0.95,'4-station','HorizontalAlignment','right','Units',...
  'normalized','FontSize',10);
text(ax,0.02,0.05,'c','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
% xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
longticks(ax,2);

%%% 4th-sta checked, PROJECTION
%mannually set the locations for each axis
c2pos=c2.Position;
ax3pos=f.ax(3).Position;
set(f.ax(4), 'position', [c2pos(1)+c2pos(3)+0.13 c2pos(2) ax3pos(3) c2pos(4)]);
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
  %   'HorizontalAlignment','right','Units','normalized','FontSize',9);
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
% if isequal(ttype,'tori')
%   xlabel(ax,'Relative origin time (s)');
% elseif isequal(ttype,'tarvl')
%   xlabel(ax,'Arrival time (s)');
% end
ylabel(ax,'Projected location (km)');
ylim(ax,[-5 2]);
yticks(ax,-5: 1: 2);
xlim(ax,cran);
longticks(ax,2);
hold(ax,'off');


%%%%%% Comparison LFE locations with contemponeous tremors
axpos = [0.08 0.055 0.415 0.3;
  c2pos(1)+c2pos(3)+0.13 0.078 ax3pos(3) 0.256;
  %  0.58 0.04 0.36 0.28;
  0.62 0.195 0.13 0.12];
for isub = 5:7
  set(f.ax(isub), 'position', axpos(isub-4,:));
end

%range as 4-s tremor density
xran = [-6 4];
yran = [-4 4];

%%%find the propagation direction for 4-s RTM, and project onto its own direction
[locxyprojmig,~,statsmig] = srcprojdistNtoNm(mig(:,seccol),mig(:,1:2),1,1);


%%%map locations of tremor
ax=f.ax(5); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
plot(ax,xcut,ycut,'k-','linew',2);
migdum = sortrows(mig,seccol);
scatter(ax,migdum(:,1),migdum(:,2), 20, migdum(:,seccol), 'filled','o',...
  'MarkerEdgeColor',[.5 .5 .5]);
%   colormap(ax,'jet');
colormap(ax,flipud(colormap(ax,'kelicol')));
%   colormap(ax,'viridis');
c1=colorbar(ax,'SouthOutside');
c1.Position = [axpos(1,1) axpos(1,2)+0.005 axpos(1,3) 0.01];
%     c.TickLabels=[];
juldate = num2str(tranmig(indmig,1));
yr = str2double(juldate(1:4));
date = str2double(juldate(5:end));
a = jul2dat(yr,date);
if a(1) == 9
  mo = 'Sep.';
elseif a(1) == 7
  mo = 'Jul.';
else
  mo = 'Mar.';
end
day = num2str(a(2));
yr = num2str(a(3));
c1.Label.String = sprintf('Time (s) since %.1f s on %s %s %s', ...
  tranmig(indmig,2),day, mo, yr);
c1.Label.FontSize = 9;
%     c.Label.String = strcat(num2str(trange(i,1)),' of HF',' (hr)');
caxis(ax,[0 tlenmig]);
% text(f.ax(1),0.85,0.1,'HF','FontSize',12,'unit','normalized');
%   text(ax,0.04,0.1,'PGC','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
axis(ax, 'equal');
axis(ax,[xran yran]);
xticks(ax,xran(1):1:xran(2));
yticks(ax,yran(1):1:yran(2));
% ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
%   medxhf = median(migdum(:,1));
%   medyhf = median(migdum(:,2));
medxhf = -2;
medyhf = -1;
[rotx, roty] = complex_rot(0,1,-statsmig.angrmse);
xarrow = [medxhf-rotx medxhf+rotx];
yarrow = [medyhf-roty medyhf+roty];
%   drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5);
p=annotation('arrow','color',[0.5 0.5 0.5],'linestyle','-','linewidth',1.5);
p.Parent = ax;
p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
% text(f.ax(1),0.5,0.1,num2str(i),'FontSize',12,'unit','normalized');
text(ax,0.98,0.15,strcat(num2str(size(migdum,1)),{' events'}),'FontSize',8,...
  'unit','normalized','horizontalalignment','right');
text(ax,0.98,0.10,strcat({'in '},num2str(tlenmig),{' s'}),'FontSize',8,...
  'unit','normalized','horizontalalignment','right');
rate = sprintf('%.3f',size(migdum,1)/tlenmig);
text(ax,0.98,0.05,strcat({'rate: '},rate),'FontSize',8,...
  'unit','normalized','horizontalalignment','right');
% text(ax,0.98,0.05,strcat(num2str(statsmig.angrmse),'$^{\,\circ}$'),'interpreter',...
%   'latex','Units','normalized','FontSize',12,'HorizontalAlignment','right');
text(ax,0.30,0.34,sprintf('%d%c',statsmig.angrmse,char(176)),...
  'Units','normalized','HorizontalAlignment','left','FontSize',10);
text(ax,0.98,0.95,'4-s tremor','FontSize',10,'unit','normalized','horizontalalignment','right',...
  'fontweight','bold');
text(ax,0.02,0.94,'e','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
hold(ax,'off');

%%%projections of tremor and those in particular within burst 181
ax=f.ax(6); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%%%plot other tremor detections along the general prop. direction
indelse=setdiff((1:length(mig))', indtmr);
scatter(ax,mig(indelse,seccol),locxyprojmig(indelse,1),...
  20,[.6 .6 .6],'filled','o','MarkerEdgeColor','w');
fitprop = feval(statsmig.fitobj,mig(:,seccol));
plot(ax,mig(:,seccol),fitprop,'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
%%%plot tremor inside 181 along its own prop. direction
timevecbst = migtmr(:,seccol);
[locxyprojtmr,~,statstmr] = srcprojdistNtoNm(timevecbst,migtmr(:,1:2),1,1);
scatter(ax,migtmr(:,seccol),locxyprojtmr(:,1),20,timevecbst,'filled','o',...
  'MarkerEdgeColor',[.5 .5 .5]);
fitprop1 = feval(statstmr.fitobj,migtmr(:,seccol));
plot(ax,migtmr(:,seccol),fitprop1,'-','linewidth',2,'color','k');
colormap(ax,'viridis');
c2=colorbar(ax,'SouthOutside');
pos = ax.Position;
% c2.Position = [pos(1), 0.11, pos(3), 0.02];
c2st = pos(1)+(tranbstbuf(181,2)-tranmig(indmig,2))/tlenmig*pos(3);
c2len = tlen(181)/tlenmig*pos(3);
c2.Position = [c2st, c1.Position(2), c2len, c1.Position(4)];
% caxis(ax,[0 tlen(181)]);
% c2.Label.String = sprintf('Time (s) since %.1f s on %s %s %s', ...
%   tranbstbuf(181,2),day, mo, yr);
caxis(ax,minmax(timevecbst'));
c2.Ticks = [1050 1250];
c2.TickLabels = ['1050'; '1250'];
c2.TickLength = 0.05;
c2.Label.FontSize = 9;
% c2.Label.String = sprintf('Time (s) since %.1f s on %s %s %s', ...
%   tranmig(indmig,2),day, mo, yr);
xran1 = [0 tlenmig];
yran1 = [-5 2];
xlim(ax,xran1);
ylim(ax,yran1);
xlabel(ax,sprintf('Time (s) since %.1f s on %s %s %s',tranmig(indmig,2),day, mo, yr));
ylabel(ax,'Projected location (km)');
text(ax,0.99,0.95,sprintf('Speed: %.1f m/s',statstmr.slope*1e3),'FontSize',8,...
  'unit','normalized','horizontalalignment','right');
%   text(ax,0.98,0.1,sprintf('SE: %.2f',stats.slopese),'FontSize',8,...
%     'unit','normalized','horizontalalignment','right');
text(ax,0.99,0.9,sprintf('Pearson: %.2f',statstmr.pearwt),'FontSize',8,...
  'unit','normalized','horizontalalignment','right');
text(ax,0.99,0.1,sprintf('Speed: %.1f m/s',statsmig.slope*1e3),'FontSize',8,...
  'unit','normalized','horizontalalignment','right');
%   text(ax,0.98,0.1,sprintf('SE: %.2f',stats.slopese),'FontSize',8,...
%     'unit','normalized','horizontalalignment','right');
text(ax,0.99,0.05,sprintf('Pearson: %.2f',statsmig.pearwt),'FontSize',8,...
  'unit','normalized','horizontalalignment','right');
text(ax,0.02,0.94,'f','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
% ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
hold(ax,'off');

%an inset on top, map view of tremors that correspond to burst win #181
timevecbst = migtmr(:,seccol)+tranmig(indmig,2)-tranbstbuf(181,2);
ax=f.ax(7); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
plot(ax,xcut,ycut,'k-','linew',1.5);
scatter(ax,migtmr(:,1),migtmr(:,2), 12, timevecbst, 'filled','o',...
  'MarkerEdgeColor',[.5 .5 .5]);
colormap(ax,'viridis');
caxis(ax,[0 tlen(181)]);
[rotx, roty] = complex_rot(0,1,-statstmr.angrmse);
xarrow = [1-rotx 1+rotx];
yarrow = [-1-roty -1+roty];
p=annotation('arrow','color','k','linestyle','-','linewidth',1.5);
p.Parent = ax;
p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
% text(ax,0.98,0.1,strcat(num2str(statstmr.angrmse),'$^{\,\circ}$'),'interpreter',...
%   'latex','Units','normalized','HorizontalAlignment','right','FontSize',11);
text(ax,0.98,0.1,sprintf('%d%c',statstmr.angrmse,char(176)),...
  'Units','normalized','HorizontalAlignment','right','FontSize',10);
text(ax,0.05,0.13,'g','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
axis(ax, 'equal');
axis(ax,[-2.5 2.5 -2.5 2.5]);
xticks(ax,xran(1):1:xran(2));
yticks(ax,yran(1):1:yran(2));
% ax.GridLineStyle = '--';
nolabels(ax,3);
longticks(ax,1.5);
hold(ax,'off');


if saveflag
  fname = strcat('shortwinsum_demo_4stmr.pdf');
  print(f.fig,'-dpdf',...
    strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
end

