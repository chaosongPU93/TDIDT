function f = plt_propdistvsdt(f,dtdist,nsep,sps,ttype,color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f = plt_propdistvsdt(f,dtdist,sps,color)
%
% Plot the along-propagation distance between the source N and source N-nsep, in
% the sequential order of origin time or arrival time (specified by 'ttype'),
% usually called to plot for many bursts combined together, and compare
% data with noise case.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/12/08
% Last modified date:   2022/12/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrow = round(size(dtdist,2)/2);
for i = 1: nrow
  
  dt = dtdist(:,(i-1)*2+1);
  distprop = dtdist(:,i*2);
  
  ax=f.ax((i-1)*2+1);
  hold(ax,'on');
  ax.Box='on'; grid(ax,'on');
  scatter(ax,dt/sps,distprop,15,color,'filled','o','MarkerEdgeColor','k');
  xlim(ax,[0 2]);
  ylim(ax,[0 4]);
%   text(ax,0.98,0.95,strcat({'med. dist of dt\leq1.0s: '},sprintf('%.2f km',median(distprop(dt/sps<=1)))),...
%     'Units','normalized','HorizontalAlignment','right','FontSize',9);
%   text(ax,0.98,0.9,strcat({'med. dist of dt\leq0.5s: '},sprintf('%.2f km',median(distprop(dt/sps<=0.5)))),...
%     'Units','normalized','HorizontalAlignment','right','FontSize',9);
  if isequal(ttype,'tori')
    xlabel(ax,sprintf('Diff. relative origin time between sources N and N-%d (s)',nsep));
  elseif isequal(ttype,'tarvl')
    xlabel(ax,sprintf('Diff. arrival time between sources N and N-%d (s)',nsep));
  end
  ylabel(ax,'Proj. dist. along prop. dir (km)');
  ax.YAxisLocation = 'right';
  longticks(ax,2);
%   nolabels(ax,1);
  
  ax=f.ax(i*2);
  hold(ax,'on');
  ax.Box='on'; grid(ax,'on');
  histogram(ax,distprop,'BinWidth',0.1,'FaceColor',color,'EdgeColor','k','Orientation','horizontal',...
    'Normalization','pdf');
  ylim(ax,[0 4]);
%   xlim(ax,[0 50]);
%   text(ax,0.98,0.95,sprintf('med. dist: %.2f km',median(distprop)),'Units','normalized',...
%     'HorizontalAlignment','right','FontSize',9);
  plot(ax,ax.XLim,[median(distprop) median(distprop)], '--', 'Color', color, 'linew', 1.5);
  % ylabel(ax,'Proj. dist. between consecutive sources (km)');
  xlabel(ax,'Probability density function');
  nolabels(ax,2);

  longticks(ax,2);
  
end