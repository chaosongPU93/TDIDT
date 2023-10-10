function [f] = plt_srcprojdist(timevec,locxy,dtime,dist,locxyproj,...
  dlocxyproj,stats,sps,ttype,wt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_srcprojdist(timevec,locxy,dtime,dist,locxyproj,...
%   dlocxyproj,stats,sps,ttype)
%
% Plot the projected distance between the source N and source N-nsep, in
% the sequential order of origin time or arrival time (specified by 'ttype'),
% along the proposed propagation direction (direction that has the min RMSE
% of the robust linear regression), in map view.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/06/22
% Last modified date:   2023/03/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('ttype','tori');  % default is sort the sources by origin time
defval('wt',ones(size(timevec))); 

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
  xran = [-4 4];
  yran = [-4 4];
  [rotx, roty] = complex_rot(0,1,-angrmse);
  xvect = [0.5-rotx 0.5+rotx];
  yvect = [-2.5-roty -2.5+roty];
  drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5);
  [rotx, roty] = complex_rot(0,1,-angslope);
  xvect = [0.5-rotx 0.5+rotx];
  yvect = [-2.5-roty -2.5+roty];
  drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[.5 .5 .5]);
  text(ax,0.7,0.2,strcat(num2str(angrmse),{'{\circ}'}),'FontSize',9,...
    'unit','normalized','horizontalalignment','left');
  text(ax,0.7,0.1,strcat(num2str(angslope),{'{\circ}'}),'FontSize',9,...
    'unit','normalized','horizontalalignment','left','color',[.5 .5 .5]);
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
  wtmax = prctile(wt,95); %use percentile in case
  refscl = wt./wtmax;
  refscl(refscl>=1) = 1;  %force the larger amp to be plotted as the same size in case of saturation
%   scatter(ax,timevec/sps,locxyproj(:,2),15,[.6 1 1],'filled','o','MarkerEdgeColor','k');
  scatter(ax,timevec/sps,locxyproj(:,1),30*refscl,[.5 .5 .5],'filled','o','MarkerEdgeColor','k');  
  % linear robust least square
  fttpfree = fittype( @(a,b,x) a*x+b);
  fitobjproj = fit(timevec/sps, locxyproj(:,1),fttpfree,'Robust','Bisquare',...
    'StartPoint',[1 1]);
  % output fit parameters
  fitproj = feval(fitobjproj,timevec/sps);
  plot(ax,timevec/sps,fitproj,'-','linewidth',2,'color','k');
  fitproj = reshape(fitproj,[],1);
%   binw = stats.gof.rmse;
%   npropbin = 4;
%   propbineg = -npropbin/2: 1: npropbin/2;
%   fitprojbin = repmat(fitproj,1,npropbin+1)+propbineg*binw.*ones(size(fitproj));
%   for i = 1: npropbin+1
%     plot(ax,timevec/sps,fitprojbin(:,i),'r--');
%   end
%   for i = 1: npropbin
%     cnt(i) = median((fitprojbin(:,i)+fitprojbin(:,i+1))/2);
%     ind = find(locxyproj(:,1)>=fitprojbin(:,i) & locxyproj(:,1)<fitprojbin(:,i+1));
%     amp(i) = median(wt(ind));
%     sprintf('%.2f \n',sort(wt(ind)))
%   end
%   ntbin = 5;
%   [~,indbin] = binxeqnum(timevec,ntbin);
%   for i = 1: ntbin
%     cnt(i) = median(timevec(indbin{i}));
%     amp(i) = median(wt(indbin{i}));
%     ampp = median(wt((locxyproj(indbin{i},1)-fitproj(indbin{i}))>=0));
%     ampn = median(wt((locxyproj(indbin{i},1)-fitproj(indbin{i}))<0));
%     ampdiff(i) = ampp-ampn;
%   end
%   text(ax,0.98,0.85,sprintf('ampsum diff: %.1f',ampp-ampn),...
%     'HorizontalAlignment','right','Units','normalized','FontSize',9);
  text(ax,0.98,0.95,sprintf('speed: %.1f km / %.1f s',range(fitproj),range(timevec/sps)),...
    'HorizontalAlignment','right','Units','normalized','FontSize',9);
  text(ax,0.98,0.90,sprintf('Pearson: %.2f',stats.pearwt),...
    'HorizontalAlignment','right','Units','normalized','FontSize',9);
% legend(ax,[p1,p2],'Along prop. direction', 'Orthogonal to prop. dir.','FontSize',9);
  if isequal(ttype,'tori')
    xlabel(ax,'Relative origin time (s)');
  elseif isequal(ttype,'tarvl')
    xlabel(ax,'Arrival time (s)');
  end
  ylabel(ax,'Projected location (km)');
  ylim(ax,yran);
  longticks(ax,2);
  
%   f1 = initfig(5,5,1,1); %initialize fig
%   ax=f1.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%   plot(ax,cnt/sps,amp,'ko-','MarkerSize',4);
%   plot(ax,cnt/sps,ampdiff,'ro-','MarkerSize',4);
%   if isequal(ttype,'tori')
%     xlabel(ax,'Relative origin time (s)');
%   elseif isequal(ttype,'tarvl')
%     xlabel(ax,'Arrival time (s)');
%   end
%   ylabel(ax,'Amp');
%   legend(ax,'median amp','median amp diff wrt prop front');
% 
%   f1 = initfig(5,5,1,1); %initialize fig
%   ax=f1.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%   plot(ax,cnt,amp,'ko-','MarkerSize',4);
%   xlabel(ax,'Projected location (km)');
%   ylabel(ax,'Med amp');
   
  ax=f.ax(5);
  hold(ax,'on');
  ax.Box='on'; grid(ax,'on');
%   scatter(ax,dtime/sps,dlocxyproj(:,2),15,[.6 1 1],'filled','o','MarkerEdgeColor','k');
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
%   p2=histogram(ax,dlocxyproj(:,2),'BinWidth',0.1,'FaceColor',[.6 1 1],'Orientation','horizontal');
  p1=histogram(ax,dlocxyproj(:,1),'BinWidth',0.1,'FaceColor','k','Orientation','horizontal');
  text(ax,0.98,0.95,sprintf('med. |dloc|: %.2f km',median(abs(dlocxyproj(:,1)))),...
    'Units','normalized','HorizontalAlignment','right','FontSize',9);
  ylim(ax,f.ax(5).YLim);
  % xlim(ax,[0 50]);
  % ylabel(ax,'Proj. dist. between consecutive sources (km)');
  xlabel(ax,'Counts');
  legend(ax,[p1],'Along proj. direction','FontSize',9,...
    'location','southeast'); %, 'Orthogonal to proj. dir.'
  longticks(ax,2);
  nolabels(ax,2);


% subplot(3,3,7); hold on
% %%% Actually used is the pre-determined best prop direc to do the fitting
% implocdum = implocst;
% for jj = 1: size(implocst,1)
%   x0 = implocdum(jj,1);
%   y0 = implocdum(jj,2);
%   [newx,newy] = coordinate_rot(x0,y0,-(angslope-90),[0 0]);
%   implocdum(jj,1) = newx;
%   implocdum(jj,2) = newy;
% end
% p2=scatter(tsplst/sps,implocdum(:,2),15,[.6 1 1],'filled','o','MarkerEdgeColor','k');
% p1=scatter(tsplst/sps,implocdum(:,1),15,[.5 .5 .5],'filled','o','MarkerEdgeColor','k');
% text(0.9,0.9,sprintf('%d',size(implocst,1)-nsep),'Units','normalized',...
%   'HorizontalAlignment','right');
% % linear robust least square
% [fitobjprop,gof,~] = fit(tsplst/sps, implocdum(:,1),fttpfree,'Robust','Bisquare',...
%   'StartPoint',[1 1]);
% % output fit parameters
% coefprop = coeffvalues(fitobjprop);
% slopeprop = coefprop(1);
% intcptprop = coefprop(2);
% fitprop = feval(fitobjprop,tsplst/sps);
% plot(tsplst/sps,fitprop,'-','linewidth',2,'color','k');
% legend([p1,p2],'Along prop. (max speed)', 'Along ort. to prop.');
% if isequal(ttype,'tori')
%   xlabel('Relative origin time (s)');
% elseif isequal(ttype,'tarvl')
%   xlabel('Arrival time (s)');
% end
% ylabel('Projected location (km)');
% grid on; box on;
% 
% subplot(3,3,8); hold on
% distort = abs(diffcustom(implocdum(:,2),nsep,'forward'));
% p2=scatter(dt/sps,distort,15,[.6 1 1],'filled','o','MarkerEdgeColor','k');
% distprop = abs(diffcustom(implocdum(:,1),nsep,'forward'));
% p1=scatter(dt/sps,distprop,15,[.5 .5 .5],'filled','o','MarkerEdgeColor','k');
% text(0.9,0.8,strcat({'med of dt\leq1.0s: '},sprintf('%.2f km',median(distprop(dt/sps<=1)))),...
%   'Units','normalized','HorizontalAlignment','right');
% text(0.9,0.7,strcat({'med of dt\leq0.5s: '},sprintf('%.2f km',median(distprop(dt/sps<=0.5)))),...
%   'Units','normalized','HorizontalAlignment','right');
% % ylim([0 3.5]);
% if isequal(ttype,'tori')
%   xlabel('Diff. relative origin time (s)');
% elseif isequal(ttype,'tarvl')
%   xlabel('Diff. arrival time (s)');
% end
% ylabel('Projected dist. between consecutive sources (km)');
% xlim([0 2]);
% % legend([p1,p2],'Along prop. (max speed)', 'Along ort. to prop.');
% grid on; box on;
% 
% subplot(3,3,9); hold on
% p2=histogram(distort,'BinWidth',0.1,'FaceColor',[.6 1 1]);
% p1=histogram(distprop,'BinWidth',0.1,'FaceColor','k');
% text(0.9,0.7,sprintf('med: %.2f km',median(distprop)),'Units','normalized',...
%   'HorizontalAlignment','right');
% % xlim([0 3.5]);
% xlabel('Projected dist. between consecutive sources (km)');
% ylabel('Counts');
% % legend([p1,p2],'Along prop. (max speed)', 'Along ort. to prop.');
% grid on; box on;

end
