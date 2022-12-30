function f = plt_absdistvsdt(f,dtdist,nsep,sps,ttype,color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f = plt_absdistvsdt(f,dtdist,sps,ttype,color)
%
% Plot the absolute distance between the source N and source N-nsep, in
% the sequential order of origin time or arrival time (specified by 'ttype'),
% usually called to plot for many bursts combined together, and compare
% data with noise case.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/12/08
% Last modified date:   2022/12/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('ttype','tori');  % default is sort the sources by origin time

nrow = round(size(dtdist,2)/2);
for i = 1: nrow
  
  dt = dtdist(:,(i-1)*2+1);
  dist = dtdist(:,i*2);
  
%   ax=f.ax((i-1)*3+1);
%   hold(ax,'on');
%   ax.Box='on'; grid(ax,'on');
%   scatter(ax,dt/sps,dist,15,color,'filled','o','MarkerEdgeColor','k');
%   ylim(ax,[0 8]);
%   xlim(ax,[0 1.5]);
% %   text(ax,0.98,0.95,strcat({'med. dist of dt\leq1.0s: '},sprintf('%.2f km',median(dist(dt/sps<=1)))),...
% %     'Units','normalized','HorizontalAlignment','right','FontSize',9);
% %   text(ax,0.98,0.9,strcat({'med. dist of dt\leq0.5s: '},sprintf('%.2f km',median(dist(dt/sps<=0.5)))),...
% %     'Units','normalized','HorizontalAlignment','right','FontSize',9);
%   if isequal(ttype,'tori')
%     xlabel(ax,sprintf('Diff. relative origin time between sources N and N-%d (s)',nsep));
%   elseif isequal(ttype,'tarvl')
%     xlabel(ax,sprintf('Diff. arrival time between sources N and N-%d (s)',nsep));
%   end
%   ylabel(ax,sprintf('Dist. between sources N and N-%d (km)',nsep));
%   ax.YAxisLocation = 'right';
%   longticks(ax,2);
  
%   ax=f.ax((i-1)*3+2);
%   hold(ax,'on');
%   ax.Box='on'; grid(ax,'on');
%   scatter(ax,dt/sps,dist,15,color,'filled','o','MarkerEdgeColor','k');
%   text(ax,0.98,0.95,strcat({'med. dist of dt\leq1.0s: '},sprintf('%.2f km',median(dist(dt/sps<=1)))),...
%     'Units','normalized','HorizontalAlignment','right','FontSize',9);
%   text(ax,0.98,0.9,strcat({'med. dist of dt\leq0.5s: '},sprintf('%.2f km',median(dist(dt/sps<=0.5)))),...
%     'Units','normalized','HorizontalAlignment','right','FontSize',9);
%   ylim(ax,[0 10]);
%   if isequal(ttype,'tori')
%     xlabel(ax,'Diff. relative origin time (s)');
%   elseif isequal(ttype,'tarvl')
%     xlabel(ax,'Diff. arrival time (s)');
%   end
%   ylabel(ax,sprintf('Dist. between sources N and N-%d (km)',nsep));
%   ax.YAxisLocation = 'right';
%   xlim(ax,[0 1.5]);
%   longticks(ax,2);

  ax=f.ax((i-1)*2+1);
  hold(ax,'on');
  ax.Box='on'; grid(ax,'on');
  histogram(ax,dist,'BinWidth',0.1,'FaceColor',color,'EdgeColor','k',...
    'Normalization','pdf'); %,'Orientation','horizontal'
  xlim(ax,[0 8]);
  ylim(ax,[0 0.7]);
%   text(ax,0.98,0.95,sprintf('med. dist: %.2f km',median(dist)),'Units','normalized',...
%     'HorizontalAlignment','right','FontSize',9);
  plot(ax,[median(dist) median(dist)],ax.YLim, '--', 'Color', color, 'linew', 1.5);
  xlabel(ax,sprintf('Dist. between sources N and N-%d (km)',nsep));
  ylabel(ax,'Probability density function');
%   ax.YAxisLocation = 'right';
%   nolabels(ax,2);
  longticks(ax,2);
  
  ax=f.ax(i*2);
  hold(ax,'on');
  ax.Box='on'; grid(ax,'on');
  histogram(ax,dt/sps,'BinWidth',4/sps,'FaceColor',color,'EdgeColor','k','Normalization','pdf');
  xlim(ax,[0 2]);
  ylim(ax,[0 3.5]);
%   text(ax,0.98,0.95,sprintf('med. dt: %.2f s',median(dt/sps)),'Units','normalized',...
%     'HorizontalAlignment','right','FontSize',9);
  plot(ax,[median(dt/sps) median(dt/sps)],ax.YLim, '--', 'Color', color, 'linew', 1.5);
  % ylabel(ax,'Dist. between consecutive sources (km)');
  ylabel(ax,'Probability density function');
  if isequal(ttype,'tori')
    xlabel(ax,sprintf('Diff. relative origin time between sources N and N-%d (s)',nsep));
  elseif isequal(ttype,'tarvl')
    xlabel(ax,sprintf('Diff. arrival time between sources N and N-%d (s)',nsep));
  end
%   nolabels(ax,2);
%   ax.YAxisLocation = 'right';
  longticks(ax,2);

end