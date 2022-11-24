function f=plt_errorof4thtarvlpred(f,pred14offtrall,offmax,color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_errorof4thtarvlpred(f,pred14offtrall,offmax,color)
%
% This function is to plot the difference in time (samples) between predicted
% arrival and finally selected peak at a 4th sta for each given source 
% location. This plot aims to lump all found sources from all burst windows 
% into the same histogram.
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
  yyaxis(ax,'left');
  histogram(ax,abs(pred14offtrall),'Normalization','probability','BinWidth',1,'FaceColor',color);
  p1=plot(ax,[offmax offmax],ax.YLim,'k--');
  %       [muHat,sigmaHat] = normfit(pred14offtr);
  %       pdfhat = normpdf(-offmax:offmax,muHat,sigmaHat);
  %       plot(ax,-offmax:offmax,pdfhat,'r','linew',2);
  %       plot(ax,[muHat muHat],ax.YLim,'r--','linew',1);
  xlabel(ax,'Abs diff in samples between predicted arrival and selected peak');
  ylabel(ax,'Probability');
  xlim(ax,[0 offmax]);
  
  yyaxis(ax,'right');
  [cdfval,x] = ecdf(abs(pred14offtrall)); %between Nth and (N-1)th source
  plot(ax,x,cdfval,'-','linew',1,'color',color);
  ylabel(ax,'Empirical CDF');
  legend(ax,[p1],'max allowed diff','Location','north');
  
end