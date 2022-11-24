function f=plt_deconpk_rat_comb4th(f,srcamprall,impindepstall,color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_deconpk_rat_comb4th(f,srcamprall,color)
%
% This function is to plot the direct source peak ratios between all 
% station pairs (12,13, and 23)and include 14. Number of 14 depends 
% on number of 4th stations that are actually useful. 
% Considering every source as a data point, this function combines all 
% burst windows together and show them as histograms.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/11/17
% Last modified date:   2022/11/17 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1: round(length(f.ax)/2)
  %hist of direct amp ratios
  ax=f.ax(i);
  hold(ax,'on');  
  grid(ax,'on');
  histogram(ax,log10(srcamprall(:,i)),'FaceColor',color,'Normalization','pdf','BinWidth',0.05);
  plot(ax,[median(log10(srcamprall(:,i))) median(log10(srcamprall(:,i)))],ax.YLim,'--',...
    'color',color,'linew',1);
  errorbar(ax,median(log10(srcamprall(:,i))), 0.9*ax.YLim(2), mad(log10(srcamprall(:,i)),1),...
    mad(log10(srcamprall(:,i)),1),'horizontal','o',...
    'color',color,'linewidth',0.8,'CapSize',5,'MarkerSize',0.5);
%   text(ax,0.95,0.9,sprintf('med=%.2f; MAD=%.2f',median(srcamprall(:,i)),mad(srcamprall(:,i), 1)),...
%     'Units','normalized','HorizontalAlignment','right');  %note here is linear scale
  if i==1
%     ylabel(ax,'# of source');
    ylabel(ax,'PDF');
    xlabel(ax,'log_{10}{src amp ratio 1/2}');
  elseif i==2
    xlabel(ax,'log_{10}{src amp ratio 1/3}');
  elseif i==3
    xlabel(ax,'log_{10}{src amp ratio 2/3}');
  elseif i==4
    xlabel(ax,'log_{10}{src amp ratio 1/4}');
  end
  xlim(ax,[-1 1]);
  
  %scatter of direct amp VS. amp ratios
  ax = f.ax(4+i);
  hold(ax,'on'); grid(ax, 'on');
  if i ==1
    scatter(ax,log10(impindepstall(:,2)),log10(srcamprall(:,i)),20,color);
    xlabel(ax,'log_{10}{src amp at 1}');
    ylabel(ax,'log_{10}{src amp ratio 1/2}');
  elseif i ==2
    scatter(ax,log10(impindepstall(:,2)),log10(srcamprall(:,i)),20,color);
    xlabel(ax,'log_{10}{src amp at 1}');
    ylabel(ax,'log_{10}{src amp ratio 1/3}');
  elseif i==3
    scatter(ax,log10(impindepstall(:,4)),log10(srcamprall(:,i)),20,color);
    xlabel(ax,'log_{10}{src amp at 2}');
    ylabel(ax,'log_{10}{src amp ratio 2/3}');
  elseif i==4
    scatter(ax,log10(impindepstall(:,2)),log10(srcamprall(:,i)),20,color);
    xlabel(ax,'log_{10}{src amp at 1}');
    ylabel(ax,'log_{10}{src amp ratio 1/4}');
  end
  axis(ax,'equal');

end