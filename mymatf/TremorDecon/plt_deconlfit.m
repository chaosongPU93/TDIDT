function [f,fitobj12,fitobj13,dyfit12,dyfit13] = plt_deconlfit(impindepst,imploc,torispl,sps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,fitobj12,fitobj13,dyfit12,dyfit13] = plt_deconlfit(impindepst,imploc,torispl,sps)
%
% Plot the source movement (change over time) from different perspectives. 
% 1. source off12 and off13 vs. origin time (you can also use the arrival time)
% 2. source relative locations N and E in map view
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/03/11
% Last modified date:   2024/03/11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
defval('sps',160);

%sort by 'tori'
[torispl, indsort] = sortrows(torispl,1);
imploc = imploc(indsort, :);
impindepst = impindepst(indsort, :);

nrow = 1;
ncol = 2;
widin = 7;
htin = 4;
pltxran = [0.1 0.95]; pltyran = [0.12 0.95]; % optimal axis location
pltxsep = 0.07; pltysep = 0.02;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

%%% in sample space
fttpfree = fittype( @(a,b,x) a*x+b);
[fitobj12,~,~] = fit(torispl/sps,impindepst(:,7)/sps,fttpfree,'Robust','Bisquare',...
  'StartPoint',[1 1]);
coef12 = coeffvalues(fitobj12);
coefci12 = confint(fitobj12);
xplt = (min(torispl): 1: max(torispl))/sps;
yfit12 = coef12(1)*xplt + coef12(2);
dyfit12 = yfit12(end)-yfit12(1);
ylow12 = coefci12(1,1)*xplt + coefci12(1,2);
yupp12 = coefci12(2,1)*xplt + coefci12(2,2);

[fitobj13,~,~] = fit(torispl/sps,impindepst(:,8)/sps,fttpfree,'Robust','Bisquare',...
  'StartPoint',[1 1]);
coef13 = coeffvalues(fitobj13);
coefci13 = confint(fitobj13);
yfit13 = coef13(1)*xplt + coef13(2);
dyfit13 = yfit13(end)-yfit13(1);
ylow13 = coefci13(1,1)*xplt + coefci13(1,2);
yupp13 = coefci13(2,1)*xplt + coefci13(2,2);

ax = f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
p1=scatter(ax,torispl/sps,impindepst(:,7)/sps,8,'r','filled','markeredgec',[.5 .5 .5]);  %off12
plot(ax,xplt,yfit12,'r-','linew',2);
patch(ax,[xplt fliplr(xplt)],[ylow12 fliplr(yupp12)],'k','Facealpha',0.2,'edgecolor','none');
p2=scatter(ax,torispl/sps,impindepst(:,8)/sps,8,'b','filled','markeredgec',[.5 .5 .5]);  %off13
plot(ax,xplt,yfit13,'b-','linew',2);
patch(ax,[xplt fliplr(xplt)],[ylow13 fliplr(yupp13)],'k','Facealpha',0.2,'edgecolor','none');
text(ax,0.98,0.15,sprintf('d(\\Delta{t}_{12}): %.1f',dyfit12),'Units','normalized',...
  'HorizontalAlignment','right');
text(ax,0.98,0.05,sprintf('d(\\Delta{t}_{13}): %.1f',dyfit13),'Units','normalized',...
  'HorizontalAlignment','right');
legend(ax,[p1,p2],sprintf('\\Delta{t}_{12}'),sprintf('\\Delta{t}_{13}'),'Location','southwest');
xlabel(ax,'Relative origin time (s)');
ylabel(ax,sprintf('\\Delta{t} (s)'));

%%% in map view
fttpfree = fittype( @(a,b,x) a*x+b);
[fitobj12,~,~] = fit(torispl/sps,imploc(:,1),fttpfree,'Robust','Bisquare',...
  'StartPoint',[1 1]);
coef12 = coeffvalues(fitobj12);
coefci12 = confint(fitobj12);
xplt = (min(torispl): 1: max(torispl))/sps;
yfit12 = coef12(1)*xplt + coef12(2);
dyfit12 = yfit12(end)-yfit12(1);
ylow12 = coefci12(1,1)*xplt + coefci12(1,2);
yupp12 = coefci12(2,1)*xplt + coefci12(2,2);

[fitobj13,~,~] = fit(torispl/sps,imploc(:,2),fttpfree,'Robust','Bisquare',...
  'StartPoint',[1 1]);
coef13 = coeffvalues(fitobj13);
coefci13 = confint(fitobj13);
yfit13 = coef13(1)*xplt + coef13(2);
dyfit13 = yfit13(end)-yfit13(1);
ylow13 = coefci13(1,1)*xplt + coefci13(1,2);
yupp13 = coefci13(2,1)*xplt + coefci13(2,2);

ax = f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
p1=scatter(ax,torispl/sps,imploc(:,1),8,'r','filled','markeredgec',[.5 .5 .5]);  %E loc
plot(ax,xplt,yfit12,'r-','linew',2);
patch(ax,[xplt fliplr(xplt)],[ylow12 fliplr(yupp12)],'k','Facealpha',0.2,'edgecolor','none');
p2=scatter(ax,torispl/sps,imploc(:,2),8,'b','filled','markeredgec',[.5 .5 .5]);  %N loc
plot(ax,xplt,yfit13,'b-','linew',2);
patch(ax,[xplt fliplr(xplt)],[ylow13 fliplr(yupp13)],'k','Facealpha',0.2,'edgecolor','none');
text(ax,0.98,0.15,sprintf('d(E): %.1f',dyfit12),'Units','normalized',...
  'HorizontalAlignment','right');
text(ax,0.98,0.05,sprintf('d(N): %.1f',dyfit13),'Units','normalized',...
  'HorizontalAlignment','right');
legend(ax,[p1,p2],'E','N','Location','southwest');
xlabel(ax,'Relative origin time (s)');
ylabel(ax,sprintf('Location (km)'));

% orient(f.fig,'landscape');
fname = 'lfemigproof.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

