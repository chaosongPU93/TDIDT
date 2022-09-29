function plt_deconpk_ratdevvsrcc(f,lgdevsrcamprall,rccpairsrcall,rcccatsrcall,color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plt_deconpk_rat_comb(srcamprall,psrcamprsall,nsrcamprsall)
%
% This function is to plot the direct and scaled deconvolved positive and 
% negative source peak ratios between all station pairs (12,13, and 23). 
% Considering every source as a data point, this function combines all 
% burst windows together and show them as histograms.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/09/28
% Last modified date:   2022/09/28 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1: 3
  ax=f.ax(i);
  hold(ax,'on');  
  grid(ax,'on');
  scatter(ax,lgdevsrcamprall(:,i),rccpairsrcall(:,i),40,'MarkerFaceColor',color,'MarkerEdgeColor',...
    'none','MarkerFaceAlpha',.2);
%   text(ax,0.95,0.1,sprintf('med=%.2f; %.2f',median(srcamprall(:,i)), median(rccpairsrcall(:,i))),...
%     'Units','normalized','HorizontalAlignment','right');
  axis(ax,[-1 1 -1 1]);
  plot(ax,ax.XLim,[median(rccpairsrcall(:,i)) median(rccpairsrcall(:,i))],'--',...
    'color',color,'linew',1);
  if i == 1
    ylabel(ax,'RCC between the same pair');
  end
  longticks(ax,2);
  
  ax=f.ax(3+i);
  hold(ax,'on');  
  grid(ax,'on');
  scatter(ax,lgdevsrcamprall(:,i),rcccatsrcall,40,'MarkerFaceColor',color,'MarkerEdgeColor','none',...
    'MarkerFaceAlpha',.2);
%   text(ax,0.95,0.1,sprintf('med=%.2f; %.2f',median(srcamprall(:,i)), median(rcccatsrcall)),...
%     'Units','normalized','HorizontalAlignment','right');
  axis(ax,[-1 1 -1 1]);
  plot(ax,ax.XLim,[median(rcccatsrcall) median(rcccatsrcall)],'--',...
    'color',color,'linew',1);
  if i == 1
    xlabel(ax,'deviation from median amp ratio 1/2 (log)');
  elseif i == 2
    xlabel(ax,'deviation from median amp ratio 1/3 (log)');
  else
    xlabel(ax,'deviation from median amp ratio 2/3 (log)');
  end
  if i == 1
    ylabel(ax,'Mean RCC of the 2 best pairs');
  end
  longticks(ax,2);
end


