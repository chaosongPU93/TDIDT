function f=plt_deconpk_sclrat_comb(f,psrcamprsall,nsrcamprsall,color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_deconpk_sclrat_comb(f,psrcamprsall,nsrcamprsall,color)
%
% This function is to plot the scaled deconvolved positive and 
% negative source peak ratios between all station pairs (12,13, and 23). 
% Considering every source as a data point, this function combines all 
% burst windows together and show them as histograms.
%
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/11/20
% Last modified date:   2022/11/20 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1: 3
  %hist of scaled positive amp ratios
  ax=f.ax(i);
  hold(ax,'on');  
  grid(ax,'on');
  histogram(ax,log10(psrcamprsall(:,i)),'FaceColor',color);
  plot(ax,[median(log10(psrcamprsall(:,i))) median(log10(psrcamprsall(:,i)))],ax.YLim,'--',...
    'color',color,'linew',1);
  errorbar(ax,median(log10(psrcamprsall(:,i))), 0.9*ax.YLim(2), mad(log10(psrcamprsall(:,i)),1),...
    mad(log10(psrcamprsall(:,i)),1),'horizontal','o',...
    'color',color,'linewidth',0.8,'CapSize',5,'MarkerSize',0.5);
%   text(ax,0.95,0.9,sprintf('med=%.2f; MAD=%.2f',median(psrcamprsall(:,i)),mad(psrcamprsall(:,i), 1)),...
%     'Units','normalized','HorizontalAlignment','right');  %note here is linear scale
  if i ==1
    xlabel(ax,'log_{10}{src scaled pos amp ratio 1/2}');
  elseif i ==2
    xlabel(ax,'log_{10}{src scaled pos amp ratio 1/3}');
  else
    xlabel(ax,'log_{10}{src scaled pos amp ratio 2/3}');
  end
  xlim(ax,[-1 1]);

  %hist of scaled negative amp ratios
  ax=f.ax(3+i);
  hold(ax,'on');  
  grid(ax,'on');
  histogram(ax,log10(nsrcamprsall(:,i)),'FaceColor',color);
  plot(ax,[median(log10(nsrcamprsall(:,i))) median(log10(nsrcamprsall(:,i)))],ax.YLim,'--',...
    'color',color,'linew',1);
  errorbar(ax,median(log10(nsrcamprsall(:,i))), 0.9*ax.YLim(2), mad(log10(nsrcamprsall(:,i)),1),...
    mad(log10(nsrcamprsall(:,i)),1),'horizontal','o',...
    'color',color,'linewidth',0.8,'CapSize',5,'MarkerSize',0.5);
%   text(ax,0.95,0.9,sprintf('med=%.2f; MAD=%.2f',median(nsrcamprsall(:,i)),mad(nsrcamprsall(:,i), 1)),...
%     'Units','normalized','HorizontalAlignment','right');  %note here is linear scale
  if i ==1
    xlabel(ax,'log_{10}{src scaled neg amp ratio 1/2}');
  elseif i ==2
    xlabel(ax,'log_{10}{src scaled neg amp ratio 1/3}');
  else
    xlabel(ax,'log_{10}{src scaled neg amp ratio 2/3}');
  end
  xlim(ax,[-1 1]);
  
end





