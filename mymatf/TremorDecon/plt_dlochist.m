function [ax,MUHAT1,SIGMAHAT1,MUHAT2,SIGMAHAT2]=plt_dlochist(ax,dloc,xran,binw,legendstr,normopt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ax,MUHAT1,SIGMAHAT1,MUHAT2,SIGMAHAT2]=plt_dlochist(ax,dloc,xran,binw,legendstr,normopt)
%
% This function is to plot the location difference between events (sign 
% included) in terms of histogram, and Gaussian fits, in two orthogonal 
% directions that are stored in different columns in data 'dloc'. The
% function reads the binwidth 'binw' and normalization option 'normopt'.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/01/23
% Last modified date:   2024/01/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('normopt','countdensity');

ndim = size(dloc,2);

%Gaussian fit
X=xran(1):binw:xran(2);
[MUHAT1,SIGMAHAT1] = normfit(dloc(:,1));
if ndim==2
  [MUHAT2,SIGMAHAT2] = normfit(dloc(:,2));
else
  MUHAT2=[]; SIGMAHAT2=[];
end

if ~isempty(ax)
  %plot histogram
  hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  p(1)=histogram(ax,dloc(:,1),'binw',binw,'normalization',normopt,'Facec','b');
  if ndim==2
    p(2)=histogram(ax,dloc(:,2),'binw',binw,'normalization',normopt,'Facec','r');
  end

  if strcmp(normopt,'countdensity')
    Y1=normpdf(X,MUHAT1,SIGMAHAT1)*size(dloc,1);
    if ndim==2
      Y2=normpdf(X,MUHAT2,SIGMAHAT2)*size(dloc,1);
    end
    ylabel(ax,'Count per km');
  elseif strcmp(normopt,'count')
    Y1=normpdf(X,MUHAT1,SIGMAHAT1)*size(dloc,1)*binw;
    if ndim==2
      Y2=normpdf(X,MUHAT2,SIGMAHAT2)*size(dloc,1)*binw;
    end
    ylabel(ax,'Count');
  elseif strcmp(normopt,'pdf')
    Y1=normpdf(X,MUHAT1,SIGMAHAT1);
    if ndim==2
      Y2=normpdf(X,MUHAT2,SIGMAHAT2);
    end
    ylabel(ax,'PDF');
  elseif strcmp(normopt,'probability')
    Y1=normpdf(X,MUHAT1,SIGMAHAT1)*binw;
    if ndim==2
      Y2=normpdf(X,MUHAT2,SIGMAHAT2)*binw;
    end
    ylabel(ax,'Probability');
  end
  plot(ax,X,Y1,'Color','b','linew',2);
  text(ax,0.02,0.95,sprintf('\\mu=%.2f, \\sigma=%.2f, med(|x|)=%.2f',...
    MUHAT1,SIGMAHAT1,median(abs(dloc(:,1)))),...
    'Units','normalized','HorizontalAlignment','left','color','b');
  if ndim==2
    plot(ax,X,Y2,'Color','r','linew',2);
    text(ax,0.02,0.9,sprintf('\\mu=%.2f, \\sigma=%.2f, med(|x|)=%.2f',...
    MUHAT2,SIGMAHAT2,median(abs(dloc(:,2)))),...
      'Units','normalized','HorizontalAlignment','left','color','r');
  end
  xlim(ax,xran);
  xlabel(ax,'Location difference (km)');
  legend(ax,p,legendstr,'Location','east');
  hold(ax,'off');
end