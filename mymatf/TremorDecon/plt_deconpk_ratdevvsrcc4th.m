function f=plt_deconpk_ratdevvsrcc4th(f,lgdevsrcamprall,rccpairsrcall,rcccatsrcall,color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_deconpk_ratdevvsrcc4th(f,lgdevsrcamprall,rccpairsrcall,rcccatsrcall,color)
%
% This function is to plot the direct source peak ratios between all 
% station pairs (12,13, and 23) including 14. VS. the RCC value at 
% the arrival time. Number of 14 depends 
% on number of 4th stations that are actually useful. 
% Considering every source as a data point, this function combines all 
% burst windows together and show them as scatter dots.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/11/17
% Last modified date:   2022/11/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1: size(lgdevsrcamprall,2)
  ax=f.ax(i);
  hold(ax,'on');  
  grid(ax,'on');
  if i<4
    ydata = rccpairsrcall(:,i);
  else
    ydata = rcccatsrcall(:,2);
  end
  scatter(ax,lgdevsrcamprall(:,i),ydata,40,'MarkerFaceColor',color,'MarkerEdgeColor',...
    'none','MarkerFaceAlpha',.2);
%   text(ax,0.95,0.1,sprintf('med=%.2f; %.2f',median(srcamprall(:,i)), median(rccpairsrcall(:,i))),...
%     'Units','normalized','HorizontalAlignment','right');
  axis(ax,[-1 1 -1 1]);
  plot(ax,ax.XLim,[median(ydata) median(ydata)],'--',...
    'color',color,'linew',1);
  if i == 1
    ylabel(ax,'RCC between the same pair');
  end
  longticks(ax,2);
  
  ax=f.ax(4+i);
  hold(ax,'on');  
  grid(ax,'on');
  if i<4
    ydata = rcccatsrcall(:,1);
  else
    ydata = rcccatsrcall(:,2);
  end
  scatter(ax,lgdevsrcamprall(:,i),ydata,40,'MarkerFaceColor',color,'MarkerEdgeColor','none',...
    'MarkerFaceAlpha',.2);
%   text(ax,0.95,0.1,sprintf('med=%.2f; %.2f',median(srcamprall(:,i)), median(rcccatsrcall)),...
%     'Units','normalized','HorizontalAlignment','right');
  axis(ax,[-1 1 -1 1]);
  plot(ax,ax.XLim,[median(ydata) median(ydata)],'--',...
    'color',color,'linew',1);
  if i == 1
    xlabel(ax,'deviation from median amp ratio 1/2 (log)');
  elseif i == 2
    xlabel(ax,'deviation from median amp ratio 1/3 (log)');
  elseif i == 3
    xlabel(ax,'deviation from median amp ratio 2/3 (log)');
  else
    xlabel(ax,'deviation from median amp ratio 1/4 (log)');
  end
  if i == 1
    ylabel(ax,'Mean RCC');
  end
  longticks(ax,2);
end
