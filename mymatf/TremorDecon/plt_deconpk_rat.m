function f=plt_deconpk_rat(f,msrcampr,madsrcampr,nsrc,color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_deconpk_rat(f,msrcampr,madsrcampr,nsrc,color)
%
% This function is to plot the direct deconvolved source peak ratios 
% between all station pairs (12,13, and 23). 
% Each burst win is represented by one dot, whose size is proportional to
% the number of sources, along with an errorbar.
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

for i = 1: length(f.ax)
  ax=f.ax(i);
  hold(ax,'on');  
  grid(ax,'on');
  for j = 1:size(msrcampr,1)
    if nsrc(j)>0
      errorbar(ax,msrcampr(j,i), j, madsrcampr(j,i),madsrcampr(j,i),'horizontal','o',...
        'color',[.5 .5 .5],'linewidth',0.8,'CapSize',5,'MarkerSize',0.5);
      scatter(ax,msrcampr(j,i),j,nsrc(j)/max(nsrc)*40,color,'o','filled');
    end
  end
  ylim(ax,[0 size(msrcampr,1)+1]);
  if size(mnsrcamprs,1) >1
    plot(ax,[wt_median(msrcampr(:,i),nsrc) wt_median(msrcampr(:,i),nsrc)],ax.YLim,'--','color',color,'linew',1);
%     text(ax,0.95,0.1,sprintf('wtmed=%.2f',wt_median(msrcampr(:,i),nsrc)),'Units','normalized',...
%       'HorizontalAlignment','right');
  end
  if i ==1
    ylabel(ax,'Burst win index');
    xlabel(ax,'med src amp ratio 1/2');
  elseif i ==2
    xlabel(ax,'med src amp ratio 1/3');
  elseif i==3
    xlabel(ax,'med src amp ratio 2/3');
  end
  

end

