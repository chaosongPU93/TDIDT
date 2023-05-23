function [f] = plt_wholewinvssubwin(f,impindep1,impindep2,sps,angproj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_wholewinvssubwin(f,impindep1,impindep2,sps,angproj)
%
% This is to plot the comparison between results from 25-s subwins and 
% whole-win. The most direct way is the map view of two colorcoded by time, 
% but time can be blurry when the time extent of the burst is long. So this
% script simply plot the comparison of time-distance along a custom direction 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/01/26
% Last modified date:   2023/01/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ista = 1;
%convert time offset to relative loc
imploc1 = off2space002(impindep1(:,7:8),sps,'interpchao',0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[impindepst1, indsort] = sortrows(impindep1, (ista-1)*2+1);
implocst1 = imploc1(indsort, :);
tsplst1 = impindepst1(:,(ista-1)*2+1);

imploc2 = off2space002(impindep2(:,7:8),sps,'interpchao',0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[impindepst2, indsort] = sortrows(impindep2, (ista-1)*2+1);
implocst2 = imploc2(indsort, :);
tsplst2 = impindepst2(:,(ista-1)*2+1);

%project loc along the SAME custom direction
[propx1,orty1,nlocxy1] = customprojection(implocst1(:,1:2),angproj);
[propx2,orty2,nlocxy2] = customprojection(implocst2(:,1:2),angproj);

fttpfree = fittype( @(a,b,x) a*x+b);

%scatter the time-distance for whole-win detection results VS subwin results
% f = initfig(6,6,1,1);

ax=f.ax(1);
hold(ax,'on');
ax.Box='on'; grid(ax,'on');
p1=scatter(ax,tsplst1/sps,propx1,15,[.5 .5 .5],'filled','o','MarkerEdgeColor','k');
p2=scatter(ax,tsplst2/sps,propx2,15,[.6 1 1],'filled','o','MarkerEdgeColor','k');
% linear robust least square
[fitobj1,gof1,output1] = fit(tsplst1/sps,propx1,fttpfree,'Robust','Bisquare','StartPoint',[1 1]);
[fitobj2,gof2,output2] = fit(tsplst2/sps,propx2,fttpfree,'Robust','Bisquare','StartPoint',[1 1]);
%fit parameters and statistics
stats1 = statsofrobustlnfit(fitobj1,gof1,output1,tsplst1/sps,propx1);
stats2 = statsofrobustlnfit(fitobj2,gof2,output2,tsplst2/sps,propx2);
%best-fit y
propx1fit = feval(fitobj1,tsplst1/sps);
propx2fit = feval(fitobj2,tsplst2/sps);

plot(ax,tsplst1/sps,propx1fit,'-','linewidth',2,'color','k');
plot(ax,tsplst2/sps,propx2fit,'-','linewidth',2,'color','c');

text(ax,0.98,0.95,sprintf('speed: %.1f km/%.1f s; Pear: %.2f',range(propx1fit),range(tsplst1/sps),...
  stats1.pearwt),'HorizontalAlignment','right','Units','normalized','FontSize',9,'color','k');
text(ax,0.98,0.88,sprintf('speed: %.1f km/%.1f s; Pear: %.2f',range(propx2fit),range(tsplst2/sps),...
  stats2.pearwt),'HorizontalAlignment','right','Units','normalized','FontSize',9,'color','c');
legend(ax,[p1,p2],'Sub-window', 'Whole-window','FontSize',9);
% xlabel(ax,'Relative origin time (s)');
xlabel(ax,'Arrival time (s)');
ylabel(ax,sprintf('Proj. loc. along %d (km)',angproj));
% ylim(ax,[-4 4]);
longticks(ax,2);












