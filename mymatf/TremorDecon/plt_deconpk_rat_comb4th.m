function [axall,mampr,madampr]=plt_deconpk_rat_comb4th(axall,srcamprall,impindepstall,color,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axall=plt_deconpk_rat_comb4th(axall,srcamprall,impindepstall,color,flag)
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

defval('flag','both');

for i = 1: size(srcamprall,2)

  %hist of direct amp ratios
  % ax = f.ax(i);
  ax = axall(i); hold(ax,'on'); grid(ax,'on'); ax.Box = 'on';
  mampr(i,1) = median(log10(srcamprall(:,i)));
  madampr(i,1) =  mad(log10(srcamprall(:,i)),1);
  histogram(ax,log10(srcamprall(:,i)),'FaceColor',color,'Normalization','count','BinWidth',0.05);
  yran = ax.YLim;
  plot(ax,[median(log10(srcamprall(:,i))) median(log10(srcamprall(:,i)))],yran,'--',...
    'color',color,'linew',1.5);
  errorbar(ax,median(log10(srcamprall(:,i))), abs_loc_on_axis(ax.YLim,0.95), mad(log10(srcamprall(:,i)),1),...
    mad(log10(srcamprall(:,i)),1),'horizontal','o',...
    'color',color,'linewidth',1.5,'CapSize',5,'MarkerSize',0.5);
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
  elseif i==4
    xlabel(ax,'log_{10}{src amp ratio PGC/KLNB}');
  end
  xlim(ax,[-1 1]);
  ylim(ax,yran);
  ax.YAxis.Exponent = 2;
  longticks(ax,2);
  hold(ax,'off');

%   keyboard
  if strcmp(flag,'both')
    %scatter of direct amp VS. amp ratios
    % ax=f.ax(size(srcamprall,2)+i);
    ax=axall(size(srcamprall,2)+i);
    hold(ax,'on'); grid(ax, 'on');
    if i ==1
      scatter(ax,log10(impindepstall(:,2)),log10(srcamprall(:,i)),20,color); %log10(srcamprall(:,i))
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
    elseif i==4
      scatter(ax,log10(impindepstall(:,2)),log10(srcamprall(:,i)),20,color);
      xlabel(ax,'log_{10}{src amp at PGC}');
      ylabel(ax,'log_{10}{src amp at KLNB}');
    end
    axis(ax,'equal');
    longticks(ax,2);
    hold(ax,'off');

  end

end