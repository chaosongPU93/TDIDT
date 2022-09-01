function [f] = plt_srcprojdist(implocst,nsep,sps,dist,dt,tsplst,ttype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_srcprojdist(implocst,impindepstst,nsep,sps,dist,dt,torisplst,ttype,mapview)
%
% Plot the projected distance between the source N and source N-nsep, in
% the sequential order of origin time or arrival time (specified by 'ttype'),
% along the proposed propagation direction (direction that has the min RMSE
% of the robust linear regression). Plot the map view of sources in
% 'relaloc' space (N and E in km).
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/06/22
% Last modified date:   2022/06/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('ttype','tori');  % default is sort the sources by origin time

angle = 0:5:355;
slope = zeros(length(angle),1);
rmse = zeros(length(angle),1);
fttpfree = fittype( @(a,b,x) a*x+b);

%%% find best angle, now is mainly to get the variation of se and slope with the trial angle
for iang = 1: length(angle)
  %%% propagation trial
  implocdum = implocst;
  for jj = 1: size(implocst,1)
    x0 = implocdum(jj,1);
    y0 = implocdum(jj,2);
    [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),[0 0]);
    implocdum(jj,1) = newx;
    implocdum(jj,2) = newy;
  end
  % linear robust least square
  [fitobj,gof,~] = fit(tsplst/sps, implocdum(:,1),fttpfree,'Robust','Bisquare',...
    'StartPoint',[1 1]);
  % output fit parameters
  coef = coeffvalues(fitobj);
  slope(iang) = coef(1);
  rmse(iang) = gof.rmse;
end

%%% best angle estimate from hf
ind = find(slope>0);
ind3 = find(rmse(ind)==min(rmse(ind)));     % rmse, Root Mean Squared Error, the smaller, the better
if length(ind3) > 1
  disp(strcat({'multiple angles have same rmse: '},num2str(angle(ind3))));
end
angrmse = angle(ind(ind3(1)));
ind6 = find(slope==max(slope)); % one with the largest slope, i.e., migrating speed
if length(ind6) > 1
  disp(strcat({'multiple angles have same slope:'},num2str(angle(ind6))));
end
angslope = angle(ind6(1));

angrmse = 140;

f.fig = figure;
f.fig.Renderer = 'painters';

subplot(3,3,1); hold on
% scatter(implocst(1+nsep:end,1),implocst(1+nsep:end,2),25,dist,'filled',...
%   'MarkerEdgeColor',[.5 .5 .5]);
% text(0.9,0.9,sprintf('%d',size(implocst,1)-nsep),'Units','normalized',...
%   'HorizontalAlignment','right');
scatter(implocst(:,1),implocst(:,2),15,tsplst/sps,'filled',...
  'MarkerEdgeColor',[.5 .5 .5]);
text(0.9,0.9,sprintf('%d',size(implocst,1)-nsep),'Units','normalized',...
  'HorizontalAlignment','right');
colormap jet
c=colorbar;
%       c.Label.String = sprintf('Distance between consecutive sources (km)');
%       caxis([0 3.5]);
if isequal(ttype,'tori')
  c.Label.String = sprintf('Relative origin time (s)');
elseif isequal(ttype,'tarvl')
  c.Label.String = sprintf('Arrival time (s)');
end
axis equal
xran = [-4 4];
yran = [-4 4];
[rotx, roty] = complex_rot(0,1,-angrmse);
xvect = [0.5-rotx 0.5+rotx];
yvect = [-2.5-roty -2.5+roty];
ax=gca;
drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5);
[rotx, roty] = complex_rot(0,1,-angslope);
xvect = [0.5-rotx 0.5+rotx];
yvect = [-2.5-roty -2.5+roty];
drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[.5 .5 .5]);
text(0.84,0.2,strcat(num2str(angrmse),{'{\circ}'}),'FontSize',10,...
  'unit','normalized','horizontalalignment','center');
text(0.84,0.1,strcat(num2str(angslope),{'{\circ}'}),'FontSize',10,...
  'unit','normalized','horizontalalignment','center','color',[.5 .5 .5]);
xlabel('E (km)');
ylabel('N (km)');
xticks(xran(1):1:xran(2));
yticks(yran(1):1:yran(2));
axis([xran yran]);
plot([-4 4],[-4 4],'k--','linew',2);
grid on; box on;

subplot(3,3,2); hold on
scatter(dt/sps,dist,15);
text(0.9,0.9,strcat({'med of dt\leq1.0s: '},sprintf('%.2f km',median(dist(dt/sps<=1)))),...
  'Units','normalized','HorizontalAlignment','right');
text(0.9,0.75,strcat({'med of dt\leq0.5s: '},sprintf('%.2f km',median(dist(dt/sps<=0.5)))),...
  'Units','normalized','HorizontalAlignment','right');
% ylim([0 3.5]);
if isequal(ttype,'tori')
  xlabel('Diff. relative origin time (s)');
elseif isequal(ttype,'tarvl')
  xlabel('Diff. arrival time (s)');
end
ylabel('Dist. between consecutive sources (km)');
xlim([0 2]);
grid on; box on;

subplot(3,3,3); hold on
histogram(dist,'BinWidth',0.1);
text(0.9,0.9,sprintf('med: %.2f km',median(dist)),'Units','normalized',...
  'HorizontalAlignment','right');
% xlim([0 3.5]);
xlabel('Dist. between consecutive sources (km)');
ylabel('Counts');
grid on; box on;

subplot(3,3,4); hold on
%%% Actually used is the pre-determined best prop direc to do the fitting
implocdum = implocst;
for jj = 1: size(implocst,1)
  x0 = implocdum(jj,1);
  y0 = implocdum(jj,2);
  [newx,newy] = coordinate_rot(x0,y0,-(angrmse-90),[0 0]);
  implocdum(jj,1) = newx;
  implocdum(jj,2) = newy;
end
p2=scatter(tsplst/sps,implocdum(:,2),15,'k');
p1=scatter(tsplst/sps,implocdum(:,1),15,[.5 .5 .5],'filled','o',...
  'MarkerEdgeColor','k');
text(0.9,0.9,sprintf('%d',size(implocst,1)-nsep),'Units','normalized',...
  'HorizontalAlignment','right');
% linear robust least square
[fitobjprop,gof,~] = fit(tsplst/sps, implocdum(:,1),fttpfree,'Robust','Bisquare',...
  'StartPoint',[1 1]);
% output fit parameters
coefprop = coeffvalues(fitobjprop);
slopeprop = coefprop(1);
intcptprop = coefprop(2);
fitprop = feval(fitobjprop,tsplst/sps);
plot(tsplst/sps,fitprop,'-','linewidth',2,'color','k');
legend([p1,p2],'Along prop. (min rmse)', 'Along ort. to prop.');
if isequal(ttype,'tori')
  xlabel('Relative origin time (s)');
elseif isequal(ttype,'tarvl')
  xlabel('Arrival time (s)');
end
ylabel('Projected location (km)');
grid on; box on;

subplot(3,3,5); hold on
distort = abs(diffcustom(implocdum(:,2),nsep,'forward'));
p2=scatter(dt/sps,distort,15,'k');
distprop = abs(diffcustom(implocdum(:,1),nsep,'forward'));
p1=scatter(dt/sps,distprop,15,[.5 .5 .5],'filled','o','MarkerEdgeColor','k');
text(0.9,0.8,strcat({'med of dt\leq1.0s: '},sprintf('%.2f km',median(distprop(dt/sps<=1)))),...
  'Units','normalized','HorizontalAlignment','right');
text(0.9,0.7,strcat({'med of dt\leq0.5s: '},sprintf('%.2f km',median(distprop(dt/sps<=0.5)))),...
  'Units','normalized','HorizontalAlignment','right');
% ylim([0 3.5]);
if isequal(ttype,'tori')
  xlabel('Diff. relative origin time (s)');
elseif isequal(ttype,'tarvl')
  xlabel('Diff. arrival time (s)');
end
ylabel('Projected dist. between consecutive sources (km)');
xlim([0 2]);
% legend([p1,p2],'Along prop. (min rmse)', 'Along ort. to prop.');
grid on; box on;

subplot(3,3,6); hold on
p2=histogram(distort,'BinWidth',0.1,'FaceColor','w');
p1=histogram(distprop,'BinWidth',0.1,'FaceColor','k');
text(0.9,0.7,sprintf('med: %.2f km',median(distprop)),'Units','normalized',...
  'HorizontalAlignment','right');
% xlim([0 3.5]);
xlabel('Projected dist. between consecutive sources (km)');
ylabel('Counts');
% legend([p1,p2],'Along prop. (min rmse)', 'Along ort. to prop.');
grid on; box on;

subplot(3,3,7); hold on
%%% Actually used is the pre-determined best prop direc to do the fitting
implocdum = implocst;
for jj = 1: size(implocst,1)
  x0 = implocdum(jj,1);
  y0 = implocdum(jj,2);
  [newx,newy] = coordinate_rot(x0,y0,-(angslope-90),[0 0]);
  implocdum(jj,1) = newx;
  implocdum(jj,2) = newy;
end
p2=scatter(tsplst/sps,implocdum(:,2),15,'k');
p1=scatter(tsplst/sps,implocdum(:,1),15,[.5 .5 .5],'filled','o',...
  'MarkerEdgeColor','k');
text(0.9,0.9,sprintf('%d',size(implocst,1)-nsep),'Units','normalized',...
  'HorizontalAlignment','right');
% linear robust least square
[fitobjprop,gof,~] = fit(tsplst/sps, implocdum(:,1),fttpfree,'Robust','Bisquare',...
  'StartPoint',[1 1]);
% output fit parameters
coefprop = coeffvalues(fitobjprop);
slopeprop = coefprop(1);
intcptprop = coefprop(2);
fitprop = feval(fitobjprop,tsplst/sps);
plot(tsplst/sps,fitprop,'-','linewidth',2,'color','k');
legend([p1,p2],'Along prop. (max speed)', 'Along ort. to prop.');
if isequal(ttype,'tori')
  xlabel('Relative origin time (s)');
elseif isequal(ttype,'tarvl')
  xlabel('Arrival time (s)');
end
ylabel('Projected location (km)');
grid on; box on;

subplot(3,3,8); hold on
distort = abs(diffcustom(implocdum(:,2),nsep,'forward'));
p2=scatter(dt/sps,distort,15,'k');
distprop = abs(diffcustom(implocdum(:,1),nsep,'forward'));
p1=scatter(dt/sps,distprop,15,[.5 .5 .5],'filled','o','MarkerEdgeColor','k');
text(0.9,0.8,strcat({'med of dt\leq1.0s: '},sprintf('%.2f km',median(distprop(dt/sps<=1)))),...
  'Units','normalized','HorizontalAlignment','right');
text(0.9,0.7,strcat({'med of dt\leq0.5s: '},sprintf('%.2f km',median(distprop(dt/sps<=0.5)))),...
  'Units','normalized','HorizontalAlignment','right');
% ylim([0 3.5]);
if isequal(ttype,'tori')
  xlabel('Diff. relative origin time (s)');
elseif isequal(ttype,'tarvl')
  xlabel('Diff. arrival time (s)');
end
ylabel('Projected dist. between consecutive sources (km)');
xlim([0 2]);
% legend([p1,p2],'Along prop. (max speed)', 'Along ort. to prop.');
grid on; box on;

subplot(3,3,9); hold on
p2=histogram(distort,'BinWidth',0.1,'FaceColor','w');
p1=histogram(distprop,'BinWidth',0.1,'FaceColor','k');
text(0.9,0.7,sprintf('med: %.2f km',median(distprop)),'Units','normalized',...
  'HorizontalAlignment','right');
% xlim([0 3.5]);
xlabel('Projected dist. between consecutive sources (km)');
ylabel('Counts');
% legend([p1,p2],'Along prop. (max speed)', 'Along ort. to prop.');
grid on; box on;


