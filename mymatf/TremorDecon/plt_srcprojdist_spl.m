function [f] = plt_srcprojdist_spl(timevec,locxy,dtime,dist,locxyproj,...
  dlocxyproj,stats,sps,ttype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_srcprojdist_spl(timevec,locxy,dtime,dist,locxyproj,...
%   dlocxyproj,stats,sps,ttype)
%
% Plot the projected distance between the source N and source N-nsep, in
% the sequential order of origin time or arrival time (specified by 'ttype'),
% along the proposed propagation direction (direction that has the min RMSE
% of the robust linear regression), in time offset space.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/06/22
% Last modified date:   2023/03/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('ttype','tori');  % default is sort the sources by origin time

%%%number of sources has to be larger than 2, otherwise there can be only 1 possible line to fit 
if size(locxy,1) <= 2
  f=[];
  return

else
    
  angrmse = stats.angrmse;
  angslope = stats.angslope;

  widin = 10;  % maximum width allowed is 8.5 inches
  htin = 6.5;   % maximum height allowed is 11 inches
  nrow = 2;
  ncol = 3;
  f = initfig(widin,htin,nrow,ncol); %initialize fig
  xran = [0.06 0.96]; yran = [0.06 0.98];
  xsep = 0.05; ysep = 0.07;
  optaxpos(f,nrow,ncol,xran,yran,xsep,ysep);

  ax=f.ax(1);
  hold(ax,'on');
  ax.Box='on'; grid(ax,'on');
  % scatter(implocst(1+nsep:end,1),implocst(1+nsep:end,2),25,dist,'filled',...
  %   'MarkerEdgeColor',[.5 .5 .5]);
  % text(0.9,0.9,sprintf('%d',size(implocst,1)-nsep),'Units','normalized',...
  %   'HorizontalAlignment','right');
  scatter(ax,locxy(:,1),locxy(:,2),15,timevec/sps,'filled',...
    'MarkerEdgeColor',[.5 .5 .5]);
  text(ax,0.9,0.95,sprintf('%d',size(locxy,1)),'Units','normalized',...
    'HorizontalAlignment','right','FontSize',9);
  oldc = colormap(ax,'kelicol');
  newc = flipud(oldc);
  colormap(ax,newc);
  c=colorbar(ax);
  %       c.Label.String = sprintf('Distance between consecutive sources (km)');
  %       caxis([0 3.5]);
  if isequal(ttype,'tori')
    c.Label.String = sprintf('Relative origin time (s)');
  elseif isequal(ttype,'tarvl')
    c.Label.String = sprintf('Arrival time (s)');
  end
  axis(ax,'equal');
  xran = [-max(abs(locxy(:)))-1 max(abs(locxy(:)))+1];
  yran = [-max(abs(locxy(:)))-1 max(abs(locxy(:)))+1];
  [rotx, roty] = complex_rot(0,6,-angrmse);
  xvect = [4-rotx 4+rotx];
  yvect = [-18-roty -18+roty];
  drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5);
  [rotx, roty] = complex_rot(0,6,-angslope);
  xvect = [4-rotx 4+rotx];
  yvect = [-18-roty -18+roty];
  drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[.5 .5 .5]);
  text(ax,0.7,0.2,strcat(num2str(angrmse),{'{\circ}'}),'FontSize',9,...
    'unit','normalized','horizontalalignment','left');
  text(ax,0.7,0.1,strcat(num2str(angslope),{'{\circ}'}),'FontSize',9,...
    'unit','normalized','horizontalalignment','left','color',[.5 .5 .5]);
  xlabel(ax,sprintf('Offset 12 (samples at %d Hz)',sps));
  ylabel(ax,sprintf('Offset 13 (samples at %d Hz)',sps));
  xticks(ax,xran(1):8:xran(2));
  yticks(ax,yran(1):8:yran(2));
  axis(ax,[xran yran]);
  plot(ax,[xran(1) xran(2)],[yran(1) yran(2)],'k--','linew',2);
  longticks(ax,2);

  ax=f.ax(2);
  hold(ax,'on');
  ax.Box='on'; grid(ax,'on');
  % dt = diffcustom(timevec,nsep,'forward');
  % dist = sqrt(diffcustom(locxy(:,1),nsep,'forward').^2 + ...
  %   diffcustom(locxy(:,2),nsep,'forward').^2 );
  scatter(ax,dtime/sps,dist,15,'k');
  text(ax,0.98,0.95,strcat({'med. dist of dt\leq1.0s: '},sprintf('%.1f spls',median(dist(dtime/sps<=1)))),...
    'Units','normalized','HorizontalAlignment','right','FontSize',9);
  text(ax,0.98,0.9,strcat({'med. dist of dt\leq0.5s: '},sprintf('%.1f spls',median(dist(dtime/sps<=0.5)))),...
    'Units','normalized','HorizontalAlignment','right','FontSize',9);
  ylim(ax,[0 45]);
  if isequal(ttype,'tori')
    xlabel(ax,'Diff. relative origin time (s)');
  elseif isequal(ttype,'tarvl')
    xlabel(ax,'Diff. arrival time (s)');
  end
  ylabel(ax,'Dist. between consecutive sources (spls)');
  ax.YAxisLocation = 'right';
  xlim(ax,[0 2]);
  longticks(ax,2);

  ax=f.ax(3);
  hold(ax,'on');
  ax.Box='on'; grid(ax,'on');
  histogram(ax,dist,'BinWidth',1,'FaceColor','w','EdgeColor','k','Orientation','horizontal');
  text(ax,0.98,0.95,sprintf('med. dist: %.1f spls',median(dist)),'Units','normalized',...
    'HorizontalAlignment','right','FontSize',9);
  ylim(ax,[0 45]);
  % xlim(ax,[0 50]);
  % ylabel(ax,'Dist. between consecutive sources (km)');
  xlabel(ax,'Counts');
  nolabels(ax,2);
  longticks(ax,2);

  ax=f.ax(4);
  hold(ax,'on');
  ax.Box='on'; grid(ax,'on');
  p2=scatter(ax,timevec/sps,locxyproj(:,2),15,[.6 1 1],'filled','o','MarkerEdgeColor','k');
  p1=scatter(ax,timevec/sps,locxyproj(:,1),15,[.5 .5 .5],'filled','o','MarkerEdgeColor','k');
  % linear robust least square
  fttpfree = fittype( @(a,b,x) a*x+b);
  fitobjproj = fit(timevec/sps, locxyproj(:,1),fttpfree,'Robust','Bisquare',...
    'StartPoint',[1 1]);
  % output fit parameters
  fitproj = feval(fitobjproj,timevec/sps);
  plot(ax,timevec/sps,fitproj,'-','linewidth',2,'color','k');
  text(ax,0.98,0.95,sprintf('speed: %.1f spls / %.1f s',range(fitproj),range(timevec/sps)),...
    'HorizontalAlignment','right','Units','normalized','FontSize',9);
  text(ax,0.98,0.88,sprintf('Pearson: %.2f',stats.pearwt),...
    'HorizontalAlignment','right','Units','normalized','FontSize',9);
  % legend(ax,[p1,p2],'Along prop. direction', 'Orthogonal to prop. dir.','FontSize',9);
  if isequal(ttype,'tori')
    xlabel(ax,'Relative origin time (s)');
  elseif isequal(ttype,'tarvl')
    xlabel(ax,'Arrival time (s)');
  end
  ylabel(ax,'Projected location (spls)');
  ylim(ax,yran);
  longticks(ax,2);

  ax=f.ax(5);
  hold(ax,'on');
  ax.Box='on'; grid(ax,'on');
  % dlocxyproj = abs(dlocxyproj);
  % distort = abs(diffcustom(impindepdum(:,2),nsep,'forward'));
  scatter(ax,dtime/sps,dlocxyproj(:,2),15,[.6 1 1],'filled','o','MarkerEdgeColor','k');
  % distproj = abs(diffcustom(impindepdum(:,1),nsep,'forward'));
  scatter(ax,dtime/sps,dlocxyproj(:,1),15,[.5 .5 .5],'filled','o','MarkerEdgeColor','k');
  text(ax,0.98,0.95,strcat({'med. |dist| of dt\leq1.0s: '},sprintf('%.1f spls',...
    median(abs(dlocxyproj(dtime/sps<=1,1))))),'Units','normalized','HorizontalAlignment','right',...
    'FontSize',9);
  text(ax,0.98,0.9,strcat({'med. |dist| of dt\leq0.5s: '},sprintf('%.1f spls',...
    median(abs(dlocxyproj(dtime/sps<=0.5,1))))),'Units','normalized','HorizontalAlignment','right',...
    'FontSize',9);
  % ylim(ax,[0 45]);
  if isequal(ttype,'tori')
    xlabel(ax,'Diff. relative origin time (s)');
  elseif isequal(ttype,'tarvl')
    xlabel(ax,'Diff. arrival time (s)');
  end
  ylabel(ax,'Diff. proj. location (spls)');
  xlim(ax,[0 2]);
  ax.YAxisLocation = 'right';
  longticks(ax,2);
  % nolabels(ax,1);

  ax=f.ax(6);
  hold(ax,'on');
  ax.Box='on'; grid(ax,'on');
  p2=histogram(ax,dlocxyproj(:,2),'BinWidth',1,'FaceColor',[.6 1 1],'Orientation','horizontal');
  p1=histogram(ax,dlocxyproj(:,1),'BinWidth',1,'FaceColor','k','Orientation','horizontal');
  text(ax,0.98,0.95,sprintf('med. |dist|: %.1f spls',median(abs(dlocxyproj(:,1)))),...
    'Units','normalized','HorizontalAlignment','right','FontSize',9);
  ylim(ax,f.ax(5).YLim);
  % xlim(ax,[0 50]);
  % ylabel(ax,'Proj. dist. between consecutive sources (km)');
  xlabel(ax,'Counts');
  legend(ax,[p1,p2],'Along proj. direction', 'Orthogonal to proj. dir.','FontSize',9,...
    'location','southeast');
  longticks(ax,2);
  nolabels(ax,2);


end

% keyboard



