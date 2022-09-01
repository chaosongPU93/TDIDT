function [f,fitobj12,fitobj13,dyfit12,dyfit13] = plt_deconlfit(impindepst,torispl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,fitobj12,fitobj13,dyfit12,dyfit13] = plt_deconlfit(impindepst,torispl)
%
% Plot the offset 23 from different perspectives. For example, 
% 1. map view in the space set by off12 and off13, or relative location N 
%   and E. 
% 2. Histogram.
% 3. offset vs. origin time (you can also use the arrival time)
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/07/01
% Last modified date:   2022/07/01 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

fttpfree = fittype( @(a,b,x) a*x+b);
[fitobj12,~,~] = fit(torispl,impindepst(:,7),fttpfree,'Robust','Bisquare',...
  'StartPoint',[1 1]);
coef12 = coeffvalues(fitobj12);
coefci12 = confint(fitobj12);
xplt = min(torispl): 1: max(torispl);
yfit12 = coef12(1)*xplt + coef12(2);
dyfit12 = yfit12(end)-yfit12(1);
ylow12 = coefci12(1,1)*xplt + coefci12(1,2);
yupp12 = coefci12(2,1)*xplt + coefci12(2,2);

[fitobj13,~,~] = fit(torispl,impindepst(:,8),fttpfree,'Robust','Bisquare',...
  'StartPoint',[1 1]);
coef13 = coeffvalues(fitobj13);
coefci13 = confint(fitobj13);
yfit13 = coef13(1)*xplt + coef13(2);
dyfit13 = yfit13(end)-yfit13(1);
ylow13 = coefci13(1,1)*xplt + coefci13(1,2);
yupp13 = coefci13(2,1)*xplt + coefci13(2,2);

f.fig=figure;
f.fig.Renderer = 'painters';
hold on
p1=scatter(torispl,impindepst(:,7),10,'r','filled');  %off12
plot(xplt,yfit12,'r-','linew',2);
patch([xplt fliplr(xplt)],[ylow12 fliplr(yupp12)],'k','Facealpha',0.2,'edgecolor','none');
p2=scatter(torispl,impindepst(:,8),10,'b','filled');  %off13
plot(xplt,yfit13,'b-','linew',2);
patch([xplt fliplr(xplt)],[ylow13 fliplr(yupp13)],'k','Facealpha',0.2,'edgecolor','none');
text(0.98,0.15,sprintf('d(off12): %.1f',dyfit12),'Units','normalized','HorizontalAlignment',...
  'right');
text(0.98,0.05,sprintf('d(off13): %.1f',dyfit13),'Units','normalized','HorizontalAlignment',...
  'right');
legend([p1,p2],'offset 12','offset 13','Location','south');
box on; grid on
xlabel('Relative origin time (samples)');
ylabel('Offset (samples)');


