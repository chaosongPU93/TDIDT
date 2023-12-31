function axall=plt_deconpk_rat_comb(axall,srcamprall,impindepstall,color,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axall=plt_deconpk_rat_comb(axall,srcamprall,impindepstall,color,flag)
%
% This function is to plot the direct and scaled deconvolved positive and 
% negative source peak ratios between all station pairs (12,13, and 23). 
% Considering every source as a data point, this function combines all 
% burst windows together and show them as histograms.
%
% --now separate the direct and scaled parts into 2 codes. This code now 
%   plots the direct ratio only, and 'plt_deconpk_sclrat' will plots the 
%   scaled ratio, positive and negative
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/09/28
% Last modified date:   2022/09/28 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('flag','both');

for i = 1: size(srcamprall,2)
  %hist of direct amp ratios
  % ax = f.ax(i);
  ax = axall(i);
  hold(ax,'on');  
  grid(ax,'on');
  histogram(ax,log10(srcamprall(:,i)),'FaceColor',color,'Normalization','count','BinWidth',0.05);
  plot(ax,[median(log10(srcamprall(:,i))) median(log10(srcamprall(:,i)))],ax.YLim,'--',...
    'color',color,'linew',1);
  errorbar(ax,median(log10(srcamprall(:,i))), 0.9*ax.YLim(2), mad(log10(srcamprall(:,i)),1),...
    mad(log10(srcamprall(:,i)),1),'horizontal','o',...
    'color',color,'linewidth',0.8,'CapSize',5,'MarkerSize',0.5);
%   text(ax,0.95,0.9,sprintf('med=%.2f; MAD=%.2f',median(srcamprall(:,i)),mad(srcamprall(:,i), 1)),...
%     'Units','normalized','HorizontalAlignment','right');  %note here is linear scale
  if i==1
    ylabel(ax,'Count');
    % ylabel(ax,'PDF');
    xlabel(ax,'log_{10}{src amp ratio PGC/SSIB}');
  elseif i==2
    xlabel(ax,'log_{10}{src amp ratio PGC/SILB}');
  elseif i==3
    xlabel(ax,'log_{10}{src amp ratio SSIB/SILB}');
  end
  xlim(ax,[-1 1]);
  
  if strcmp(flag,'both')
    %scatter of direct amp VS. amp ratios
    ax = f.ax(size(srcamprall,2)+i);
    hold(ax,'on'); grid(ax, 'on');
    if i ==1
      scatter(ax,log10(impindepstall(:,2)),log10(srcamprall(:,i)),20,color);
      xlabel(ax,'log_{10}{src amp at PGC}');
      ylabel(ax,'log_{10}{src amp at SSIB}');
    elseif i ==2
      scatter(ax,log10(impindepstall(:,2)),log10(srcamprall(:,i)),20,color);
      xlabel(ax,'log_{10}{src amp at PGC}');
      ylabel(ax,'log_{10}{src amp at SILB}');
    elseif i==3
      scatter(ax,log10(impindepstall(:,4)),log10(srcamprall(:,i)),20,color);
      xlabel(ax,'log_{10}{src amp at SSIB}');
      ylabel(ax,'log_{10}{src amp at SILB}');
    end
    axis(ax,'equal');
  end
  
end
  
  
