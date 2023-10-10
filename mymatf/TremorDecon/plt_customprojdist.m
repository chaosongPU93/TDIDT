function [f] = plt_customprojdist(timevec,locxy,dtime,dist,locxyproj,...
  dlocxyproj,projang,sps,ttype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_customprojdist(timevec,locxy,dtime,dist,locxyproj,...
%   dlocxyproj,projang,sps,ttype)
%
% Plot the projected distance between the source N and source N-nsep, in
% the sequential order of origin time or arrival time (specified by 'ttype'),
% along the CUSTOM direction (projang), in map view.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/04/27
% Last modified date:   2023/04/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('ttype','tori');  % default is sort the sources by origin time

%%%number of sources has to be larger than 2, otherwise there can be only 1 possible line to fit 
if size(locxy,1) <= 2
  f=[];
  return

else

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
  xran = [-4 4];
  yran = [-4 4];
  [rotx, roty] = complex_rot(0,1,-projang);
  xvect = [0.5-rotx 0.5+rotx];
  yvect = [-2.5-roty -2.5+roty];
  drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5);
  text(ax,0.7,0.2,strcat(num2str(projang),{'{\circ}'}),'FontSize',9,...
    'unit','normalized','horizontalalignment','left');
  xlabel(ax,'E (km)');
  ylabel(ax,'N (km)');
  xticks(ax,xran(1):1:xran(2));
  yticks(ax,yran(1):1:yran(2));
  axis(ax,[xran yran]);
  plot(ax,[-4 4],[-4 4],'k--','linew',2);
  longticks(ax,2);

  ax=f.ax(2);
  hold(ax,'on');
  ax.Box='on'; grid(ax,'on');
  % dt = diffcustom(tsplst,nsep,'forward');
  % dist = sqrt(diffcustom(implocst(:,1),nsep,'forward').^2 + ...
  %   diffcustom(implocst(:,2),nsep,'forward').^2 );
  scatter(ax,dtime/sps,dist,15,'k');
  text(ax,0.98,0.95,strcat({'med. dist of dt\leq1.0s: '},sprintf('%.2f km',median(dist(dtime/sps<=1)))),...
    'Units','normalized','HorizontalAlignment','right','FontSize',9);
  text(ax,0.98,0.9,strcat({'med. dist of dt\leq0.5s: '},sprintf('%.2f km',median(dist(dtime/sps<=0.5)))),...
    'Units','normalized','HorizontalAlignment','right','FontSize',9);
  ylim(ax,[0 6]);
  if isequal(ttype,'tori')
    xlabel(ax,'Diff. relative origin time (s)');
  elseif isequal(ttype,'tarvl')
    xlabel(ax,'Diff. arrival time (s)');
  end
  ylabel(ax,'Dist. between consecutive sources (km)');
  ax.YAxisLocation = 'right';
  xlim(ax,[0 2]);
  longticks(ax,2);

  ax=f.ax(3);
  hold(ax,'on');
  ax.Box='on'; grid(ax,'on');
  histogram(ax,dist,'BinWidth',0.1,'FaceColor','w','EdgeColor','k','Orientation','horizontal');
  text(ax,0.98,0.95,sprintf('med. dist: %.2f km',median(dist)),'Units','normalized',...
    'HorizontalAlignment','right','FontSize',9);
  ylim(ax,[0 6]);
  % xlim(ax,[0 50]);
  % ylabel(ax,'Dist. between consecutive sources (km)');
  xlabel(ax,'Counts');
  nolabels(ax,2);
  longticks(ax,2);

  ax=f.ax(4);
  hold(ax,'on');
  ax.Box='on'; grid(ax,'on');
  % p2=scatter(ax,timevec/sps,locxyproj(:,2),15,[.6 1 1],'filled','o','MarkerEdgeColor','k');
  p1=scatter(ax,timevec/sps,locxyproj(:,1),15,[.5 .5 .5],'filled','o','MarkerEdgeColor','k');
  % linear robust least square
  fttpfree = fittype( @(a,b,x) a*x+b);
  fitobjproj = fit(timevec/sps, locxyproj(:,1),fttpfree,'Robust','Bisquare',...
    'StartPoint',[1 1]);
  % output fit parameters
  fitproj = feval(fitobjproj,timevec/sps);
  plot(ax,timevec/sps,fitproj,'-','linewidth',2,'color','k');
  text(ax,0.98,0.95,sprintf('speed: %.1f km / %.1f s',range(fitproj),range(timevec/sps)),...
    'HorizontalAlignment','right','Units','normalized','FontSize',9);
  if isequal(ttype,'tori')
    xlabel(ax,'Relative origin time (s)');
  elseif isequal(ttype,'tarvl')
    xlabel(ax,'Arrival time (s)');
  end
  ylabel(ax,'Projected location (km)');
  ylim(ax,yran);
  longticks(ax,2);

  ax=f.ax(5);
  hold(ax,'on');
  ax.Box='on'; grid(ax,'on');
  % scatter(ax,dtime/sps,dlocxyproj(:,2),15,[.6 1 1],'filled','o','MarkerEdgeColor','k');
  scatter(ax,dtime/sps,dlocxyproj(:,1),15,[.5 .5 .5],'filled','o','MarkerEdgeColor','k');
  text(ax,0.98,0.95,strcat({'med. |dloc| of dt\leq1.0s: '},sprintf('%.2f km',...
    median(abs(dlocxyproj(dtime/sps<=1,1))))),'Units','normalized','HorizontalAlignment','right',...
    'FontSize',9);
  text(ax,0.98,0.9,strcat({'med. |dloc| of dt\leq0.5s: '},sprintf('%.2f km',...
    median(abs(dlocxyproj(dtime/sps<=0.5,1))))),'Units','normalized','HorizontalAlignment','right',...
    'FontSize',9);
  % ylim(ax,[0 6]);
  if isequal(ttype,'tori')
    xlabel(ax,'Diff. relative origin time (s)');
  elseif isequal(ttype,'tarvl')
    xlabel(ax,'Diff. arrival time (s)');
  end
  ylabel(ax,'Diff. proj. location (km)');
  xlim(ax,[0 2]);
  ax.YAxisLocation = 'right';
  longticks(ax,2);
  % nolabels(ax,1);

  ax=f.ax(6);
  hold(ax,'on');
  ax.Box='on'; grid(ax,'on');
  % p2=histogram(ax,dlocxyproj(:,2),'BinWidth',0.1,'FaceColor',[.6 1 1],'Orientation','horizontal');
  p1=histogram(ax,dlocxyproj(:,1),'BinWidth',0.1,'FaceColor','k','Orientation','horizontal');
  text(ax,0.98,0.95,sprintf('med. |dloc|: %.2f km',median(abs(dlocxyproj(:,1)))),...
    'Units','normalized','HorizontalAlignment','right','FontSize',9);
  ylim(ax,f.ax(5).YLim);
  % xlim(ax,[0 50]);
  % ylabel(ax,'Proj. dist. between consecutive sources (km)');
  xlabel(ax,'Counts');
  % legend(ax,[p1,p2],'Along proj. direction', 'Orthogonal to proj. dir.','FontSize',9,...
  %   'location','southeast');
  legend(ax,p1,'Along proj. direction','FontSize',9,'location','southeast');    
  longticks(ax,2);
  nolabels(ax,2);

end
