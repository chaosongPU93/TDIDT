function f=plt_deconpk_sclrat(f,mpsrcamprs,madpsrcamprs,mnsrcamprs,madnsrcamprs,nsrc,color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_deconpk_rat(msrcampr,madsrcampr,nsrc,symbol,mpsrcamprs,madpsrcamprs,...
%   mnsrcamprs,madnsrcamprs)
%
% Similar to 'plt_deconpk_rat', this function is to plot the SCALED 
% deconvolved positive and negative source peak ratios between all
% station pairs (12,13, and 23). 
% Each burst win is represented by one dot, whose size is proportional to
% the number of sources, along with an errorbar.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/11/20
% Last modified date:   2022/11/20 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1: 3
  ax=f.ax(i);
  hold(ax,'on');  
  grid(ax,'on');
  for j = 1:size(mpsrcamprs,1)
    if nsrc(j)>0
      errorbar(ax,mpsrcamprs(j,i), j, madpsrcamprs(j,i),madpsrcamprs(j,i),'horizontal','o',...
        'color',[.5 .5 .5],'linewidth',0.8,'CapSize',5,'MarkerSize',0.5);
      scatter(ax,mpsrcamprs(j,i),j,nsrc(j)/max(nsrc)*40,color,'o','filled');
    end
  end
  ylim(ax,[0 size(msrcampr,1)+1]);
  if size(mnsrcamprs,1) >1
    plot(ax,[wt_median(mpsrcamprs(:,i),nsrc) wt_median(mpsrcamprs(:,i),nsrc)],ax.YLim,'--','color',color,'linew',1);
%     text(ax,0.95,0.1,sprintf('wtmed=%.2f',wt_median(mpsrcamprs(:,i),nsrc)),'Units','normalized',...
%       'HorizontalAlignment','right');
  end
  if i ==1
    xlabel(ax,'med src scaled pos amp ratio 1/2');
  elseif i ==2
    xlabel(ax,'med src scaled pos amp ratio 1/3');
  else
    xlabel(ax,'med src scaled pos amp ratio 2/3');
  end

  ax=f.ax(3+i);
  hold(ax,'on');  
  grid(ax,'on');
  for j = 1:size(mnsrcamprs,1)
    if nsrc(j)>0
      errorbar(ax,mnsrcamprs(j,i), j, madnsrcamprs(j,i),madnsrcamprs(j,i),'horizontal','o',...
        'color',[.5 .5 .5],'linewidth',0.8,'CapSize',5,'MarkerSize',0.5);
      scatter(ax,mnsrcamprs(j,i),j,nsrc(j)/max(nsrc)*40,color,'o','filled');
    end
  end
  ylim(ax,[0 size(msrcampr,1)+1]);
  if size(mnsrcamprs,1) >1
    plot(ax,[wt_median(mnsrcamprs(:,i),nsrc) wt_median(mnsrcamprs(:,i),nsrc)],ax.YLim,'--','color',color,'linew',1);
%     text(ax,0.95,0.1,sprintf('wtmed=%.2f',wt_median(mnsrcamprs(:,i),nsrc)),'Units','normalized',...
%       'HorizontalAlignment','right');
  end
  if i ==1
    xlabel(ax,'med src scaled neg amp ratio 1/2');
  elseif i ==2
    xlabel(ax,'med src scaled neg amp ratio 1/3');
  else
    xlabel(ax,'med src scaled neg amp ratio 2/3');
  end

end



