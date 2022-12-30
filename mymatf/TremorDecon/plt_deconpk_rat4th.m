function f=plt_deconpk_rat4th(f,msrcampr,madsrcampr,nsrc,color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_deconpk_rat4th(f,msrcampr,madsrcampr,nsrc,color)
%
% Similar to 'plt_deconpk_rat', this function is to plot the direct source
% peak ratios between all station pairs (12,13, and 23) and includes 14. 
% Number of 14 depends on number of 4th stations that are actually useful.
% Each burst win is represented by one dot, whose size is proportional to
% the number of sources, along with an errorbar. NOTE that as of 2022/12/22,
% the median amp ratio and its std is in log scale!
% 
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/11/17
% Last modified date:   2022/11/17
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
  if size(msrcampr,1) >1
    plot(ax,[wt_median(msrcampr(:,i),nsrc) wt_median(msrcampr(:,i),nsrc)],ax.YLim,'--','color',color,'linew',1);
  end
  if i==1
    ylabel(ax,'Burst win index');
    xlabel(ax,'log_{10}{med src amp ratio 1/2}');
  elseif i==2
    xlabel(ax,'log_{10}{med src amp ratio 1/3}');
  elseif i==3
    xlabel(ax,'log_{10}{med src amp ratio 2/3}');
  elseif i==4
    xlabel(ax,'log_{10}{med src amp ratio 1/4}');
  end
  
end


