function f=plt_srcsharpness(sharp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,distprop,distort,angrmse,fitobjprop,gof,output] = plt_srcprojdist(implocst,nsep,sps,dist,dt,tsplst,ttype)
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
widin = 10;
htin = 5;
nrow = 1;
ncol = 2;
f = initfig(widin,htin,nrow,ncol);

xran = [0.1 0.9]; yran = [0.1 0.9];
xsep = 0.08; ysep = 0.08;
optaxpos(f,nrow,ncol,xran,yran,xsep,ysep);

nsrc = size(sharp,1);
ax = f.ax(1); ax.Box='on'; grid(ax,'on');
hold(ax,'on');
scatter(ax,1:nsrc,sharp(:,1),20,'r','filled');
scatter(ax,1:nsrc,sharp(:,2),20,'b','filled');
scatter(ax,1:nsrc,sharp(:,3),20,'k','filled');
legend(ax,'PGC','SSIB','SILB','Location','southeast');
xlabel(ax,'Source # (as order in decon)');
ylabel(ax,'Focal length from a parabola fitting');
hold(ax,'off');

ax = f.ax(2); ax.Box='on'; grid(ax,'on');
hold(ax,'on');
histogram(ax,sharp(:,1),'FaceColor','r','Orientation','horizontal','binw',0.0005);
histogram(ax,sharp(:,2),'FaceColor','b','Orientation','horizontal','binw',0.0005);
histogram(ax,sharp(:,3),'FaceColor','k','Orientation','horizontal','binw',0.0005);
plot(ax,ax.XLim,[median(sharp(:,1)) median(sharp(:,1))],'r--','linew',2);
plot(ax,ax.XLim,[median(sharp(:,2)) median(sharp(:,2))],'b--','linew',2);
plot(ax,ax.XLim,[median(sharp(:,3)) median(sharp(:,3))],'k--','linew',2);
xlabel(ax,'Count');
% ylabel(ax,'Focal length from a parabola fitting');
ylim(ax,f.ax(1).YLim);
hold(ax,'off');

hdl = supertit(f.ax,'sharpness of grouped peaks of res-wlet CC');
movev(hdl,0.2);

% keyboard




