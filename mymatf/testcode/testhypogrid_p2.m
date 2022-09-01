% testhypogrid_p2
%
% this script will test a few weird behaviors of Hypoinverse, as a
% continuity of 'testhypogrid.m'
% 
% 1. How is it different by carrying out 40 iterations directly from
%    20 iterations then a refinement of another 20 iterations whose
%    first iteration is based on the inverted depth of the former 20
%    iterations? Ideally, if Hypoinverse is stable and reliable 
%    during the inversion, the 2 ways should produce the same result.
% --Found the problem: the endline of each event in 'chat00p' contains
%   a starting guess of the event origin time, lat, lon and depth. In
%   our case, we fix the depth. But this line is read and updated for
%   each round of the iteration. When we do the refinement, we only 
%   set the depth to be the same as the inversion from the 20 iterations,
%   but the rest of this line remains the same as the 1st iteration,
%   makes a difference in 'chat00p' for the refinement and 'chat20p' if
%   extend to 40 iterations. To clarify, start inverting lon and lat from
%   'chat00p', and use the output 'chat00o' (simplify it to get chat00s)
%   to update the phase file 'chat01p', which is different from 'chat00p'
%   ONLY at the end line. And the information of that line is the inversion
%   result same as 'chat00o' or 'chat00s', except the depth which comes
%   from interpolating the inverted lat and lon onto the slab interface.
% --So if I did the refinement by using the same information from the last
%   inversion, it is no different from just doing 40 iterations.
% --But I thought the end line is not important, but looks like it needs a 
%   starting model to get the travel time estimate, derivatives and so on,
%   so that line is actually read and used in the 1st iteration, and that's
%   why it is also updated.
% --But looking like the two inverted locations (lat, lon, dep) for the same
%   source from 40 iterations or refinement after 20 iterations, the lat, lon
%   and origin time changes significantly. However, both solution meets 
%   the criteria: 2 observed time offsets and source on the slab (can be seen
%   by the small diff between slab interface and inverted depth), so they
%   are all 'correct'. 
% --This indicates the fact that, the solution is NOT unique. There is a
%   tradeoff between location on the slab and origin time. Image if sources
%   are all on the same plane as the stations, the 2 constant offsets would
%   outline a hypobola on the plane. Similarly, when the plane is tilted, and
%   the velocity is layered, there can be a more complex structure on which
%   all sources can meet your restrictions roughly within your resolution
%   (error) limit. Therefore, you may end up with non-unique solutions for
%   the same detection represented by 2 time offsets. In other words, the 
%   transformation between time offset and locations could be NON-unique
%   sometimes.
% --The reason why refinement sometimes did a better job is because it happens
%   to force the starting model loc and origin time to be more 'consistent'
%   with the rest.
% --Final lesson, the strategy I adopted, ie., explore the max range of a grid
%   with no 'jumping' to another solution is still sensable, as a result,
%   the transformed locaton grid looks regular in space. +/-13 samples at 20 sps 
%
%
% 2. Plot several cross-section profiles of the slab, and check the 
%    gradient of the depth and estimate the dip angle of the slab.
%    With that information in mind, is plotting the map view as we 
%    always do is acceptable instead of a fault view? If not, how much
%    is the distortion?
% --Looks like the distortion is not big, either using the 'true' depth from 
%   interpolation at the slab interface, or using a plane model with a median
%   dipping angle.
% --The byproduct of this experiment is the cross-section profile, and the spatial
%   derivative of the depth of slab interface. And slab1.0 is much smoother than 
%   the old slab grid used by John Armb.
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/03/16
% Last modified date:   2022/03/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
% close all;

%% STEP 1
%fix depth to 37 km, 20 iterations, directly read the 'chat??s', more detailed result
n=20;
dep20 = zeros(625,n);
for i = 0:1:n-1
  if i<10
    suf = strcat('0',num2str(i),'s');
  else
    suf = strcat(num2str(i),'s');
  end
  iter = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/20iter/chat',suf));
  dep20(:,i+1) = iter(:,8);
end

%fix depth to 37 km, 40 iterations, directly read the 'chat??s', more detailed result
n=40;
dep40 = zeros(625,n);
for i = 0:1:n-1
  if i<10
    suf = strcat('0',num2str(i),'s');
  else
    suf = strcat(num2str(i),'s');
  end
  iter = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/40iter/chat',suf));
  dep40(:,i+1) = iter(:,8);
end

%fix depth to 37 km, 20 iterations, then refine by 20 more iterations directly read the
%'chat??s', more detailed result
n=20;
dep20ref = zeros(625,n);
for i = 0:1:n-1
  if i<10
    suf = strcat('0',num2str(i),'s');
  else
    suf = strcat(num2str(i),'s');
  end
  iter = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/ref20iter/chat',suf));
  dep20ref(:,i+1) = iter(:,8);
end

%randomly choose several sources to plot the depth variation with iteration
iplt = randi(625,10,1);
figure
subplot(1,2,1)
for i = 1: 10
  plot(dep20(iplt(i), :),'linew',2); hold on
  plot(dep40(iplt(i), :),'linew',1);
end

subplot(1,2,2)
for i = 1: 10
  plot(1:20,dep20(iplt(i), :),'linew',2); hold on
  plot(21:40,dep20ref(iplt(i), :),'linew',2);
  plot(dep40(iplt(i), :),'linew',1);
end


%% STEP 2
%%%read in 2 slab grids that are actively using
%Old one, used by John Armb
oslab = load('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/xyz.grid');
oloc0 = off2space002([0 0],20,'interpArmb',0); % loc has 8 cols, format: dx,dy,lon,lat,dep,tori,off12,off13
[dx, dy] = absloc2relaloc(-oslab(:,1),oslab(:,2),oloc0(:,3),oloc0(:,4));
oslab = [dx dy oslab];

%create a rectangle grid for interpolation, necessary for the later gradient computation
xgrdran = [-40 40];
ygrdran = [-40 40];
dxgrd=1;
dygrd=1;
x = xgrdran(1):dxgrd:xgrdran(2);
y = ygrdran(2):-dygrd:ygrdran(1);
nx = length(x);
ny = length(y);
[xgrd, ygrd] = meshgrid(x,y);
xmesh = reshape(xgrd,[],1);
ymesh = reshape(ygrd,[],1);
% depnew = interp2(grd(:,1),grd(:,2),grd(:,3),lonmesh,latmesh,'spline');
%interpolate for the depth within the rect grid
Fo = scatteredInterpolant(oslab(:,1),oslab(:,2),oslab(:,5),'linear','linear');
zomesh = Fo(xmesh,ymesh);
zogrd = reshape(zomesh,nx,ny);

%New one, used by Chao
nslab = load('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/slab1.0.grid');
nloc0 = off2space002([0 0],20,'interpchao',0);
[dx, dy] = absloc2relaloc(-nslab(:,1),nslab(:,2),nloc0(:,3),nloc0(:,4));
nslab = [dx dy nslab];
%interpolate for the depth within the rect grid
Fn = scatteredInterpolant(nslab(:,1),nslab(:,2),nslab(:,5),'linear','linear');
znmesh = Fn(xmesh,ymesh);
zngrd = reshape(znmesh,nx,ny);

%All detections around 002, just to have a sense on the spatial extent that detections cover
hftime= load(fullfile('/home/data2/chaosong/matlab/allan/data-no-resp/PGCtrio/MAPS',...
  'tdectimeori_002.up.lo2.1.cc0.44.22.ms26_1.25-6.5_4s40sps4add'));
sps=40; iup=4;
offint = round(hftime(:, 1:2)*iup);
nloc = off2space002(offint,sps*iup,'interpchao',0); % loc has 8 cols, format: dx,dy,lon,lat,dep,tori,off12,off13
[k1n,~] = boundary(nloc(:,3),nloc(:,4),0.5);

ind = find(offint(:,1)>=-23*iup & offint(:,1)<=25*iup & ...
           offint(:,2)>=-24*iup & offint(:,2)<=24*iup);
oloc = off2space002(offint(ind,:),sps*iup,'interpArmb',0);
[k1o,~] = boundary(oloc(:,3),oloc(:,4),0.5);

%highest density region outlined by an ellipse
x0 = 0.2; 
y0 = 0.2;
rotang = 45;
semia = 1.75;
semib = 1.0;
[xell, yell] = ellipse_chao(x0,y0,semia,semib,0.01,rotang,[x0,y0]);

%highest density region outlined by a rectangle
ns = [-3 3];
ew = [-6 2];
x0 = mean(ew);
y0 = mean(ns);
wid = 8;
hgt = 6;
[xrec, yrec] = rectangle_chao(x0,y0,wid,hgt,0.01);

%a cross-section profile with angle of 45 deg
xprof=-20:0.1:15;
rotang=45;
a=tand(rotang);
yprof=linefcn(xprof,a,0);


%%
%what is the slab model like, in terms of lon, lat, dep
figure
ax=subplot(1,2,1);
scatter(ax,-oslab(:,3),oslab(:,4),20,oslab(:,5),'s','filled'); hold(ax,'on');
scatter(ax,oloc0(:,3),oloc0(:,4),30,'k','filled');
plot(ax,oloc(k1o,3),oloc(k1o,4),'k');
colormap('jet');
c=colorbar;
c.Label.String = strcat('Depth (km) to interface');
% cran = c.Limits;
cran = [30 50];
caxis(ax,cran);
yran = [48.2 48.8];
xran = [-124 -123];
xlim(ax,xran);
ylim(ax,yran);
title(ax,'Slab model xyz.grid in /LOCpgsssitest');
xlabel(ax,'Lon');
ylabel(ax,'Lat');
hold(ax,'off');

ax=subplot(1,2,2);
scatter(ax,-nslab(:,3),nslab(:,4),20,nslab(:,5),'s','filled'); hold(ax,'on');
scatter(ax,nloc0(:,3),nloc0(:,4),30,'k','filled');
plot(ax,nloc(k1n,3),nloc(k1n,4),'k');
colormap('jet');
c=colorbar;
c.Label.String = strcat('Depth (km) to interface');
caxis(ax,cran);
xlim(ax,xran);
ylim(ax,yran);
title(ax,'Slab model slab1.0.grid in /LOCpgsssitest');
xlabel(ax,'Lon');
ylabel(ax,'Lat');
hold(ax,'off');


%what is the slab model like, in terms of relative distance N, E, dep
figure
ax=subplot(1,2,1);
scatter(ax,oslab(:,1),oslab(:,2),20,oslab(:,5),'s','filled'); hold(ax,'on');
scatter(ax,oloc0(:,1),oloc0(:,2),30,'k','filled');
plot(ax,oloc(k1o,1),oloc(k1o,2),'k');
colormap('jet');
c=colorbar;
c.Label.String = strcat('Depth (km) to interface');
% cran = c.Limits;
cran = [30 50];
caxis(ax,cran);
axis equal
xlim(ax,xgrdran);
ylim(ax,ygrdran);
title(ax,'Slab model xyz.grid in /LOCpgsssitest');
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
hold(ax,'off');

ax=subplot(1,2,2);
scatter(ax,nslab(:,1),nslab(:,2),20,nslab(:,5),'s','filled'); hold(ax,'on');
scatter(ax,nloc0(:,1),nloc0(:,2),30,'k','filled');
plot(ax,nloc(k1n,1),nloc(k1n,2),'k');
colormap('jet');
c=colorbar;
c.Label.String = strcat('Depth (km) to interface');
caxis(ax,cran);
axis equal
xlim(ax,xgrdran);
ylim(ax,ygrdran);
title(ax,'Slab model slab1.0.grid in /LOCpgsssitest');
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
hold(ax,'off');


%how does my interpolation work, I need this rect grid to compute gradient at any direction.
%what is the slab model like, in terms of relative distance N, E, dep
figure
ax=subplot(1,2,1);
scatter(ax,xmesh,ymesh,20,zomesh,'s','filled'); hold(ax,'on');
% imagesc(ax,[-40 40], [40 -40],zogrd);
% ax.YDir = 'reverse';
scatter(ax,oloc0(:,1),oloc0(:,2),30,'k','filled');
plot(ax,oloc(k1o,1),oloc(k1o,2),'k');
colormap('jet');
c=colorbar;
c.Label.String = strcat('Depth (km) to interface');
% cran = c.Limits;
cran = [30 50];
caxis(ax,cran);
axis equal
xlim(ax,xgrdran);
ylim(ax,ygrdran);
title(ax,'Slab model xyz.grid in /LOCpgsssitest');
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
hold(ax,'off');

ax=subplot(1,2,2);
scatter(ax,xmesh,ymesh,20,znmesh,'s','filled'); hold(ax,'on');
scatter(ax,nloc0(:,1),nloc0(:,2),30,'k','filled');
plot(ax,nloc(k1n,1),nloc(k1n,2),'k');
colormap('jet');
c=colorbar;
c.Label.String = strcat('Depth (km) to interface');
caxis(ax,cran);
axis equal
xlim(ax,xgrdran);
ylim(ax,ygrdran);
title(ax,'Slab model slab1.0.grid in /LOCpgsssitest');
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
hold(ax,'off');

%% spatial gradient of depth 
%%%gradient works for both a grid matrix or a single vector
[dzdxo, dzdyo] = gradient(zogrd,dxgrd,-dygrd);
[dzdxn, dzdyn] = gradient(zngrd,dxgrd,-dygrd);

%median of the gradient vector, which is the ortogonal of the strike, clockwise from N
strikeort = 90-median(median(atand(dzdyn./dzdxn)));  
strike = strikeort+270;

%gradient vector, note that gradient points to the acsending direction
figure
ax=subplot(1,2,1);
contour(ax,xgrd,ygrd,zogrd,10); hold(ax,'on');
colormap('jet');
quiver(ax,xgrd,ygrd,dzdxo,dzdyo,2,'color',[.6 .6 .6]);
scatter(ax,oloc0(:,1),oloc0(:,2),30,'k','filled');
plot(ax,oloc(k1o,1),oloc(k1o,2),'k');
plot(ax,xell,yell,'k-','linew',1);
plot(ax,xrec,yrec,'k--','linew',1);
axis equal
xlim(ax,xgrdran);
ylim(ax,ygrdran);
title(ax,'Gradient vector, xyz.grid');
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
hold(ax,'off');

ax=subplot(1,2,2);
contour(ax,xgrd,ygrd,zngrd,10); hold(ax,'on');
colormap('jet');
quiver(ax,xgrd,ygrd,dzdxn,dzdyn,2,'color',[.6 .6 .6]);
scatter(ax,nloc0(:,1),nloc0(:,2),30,'k','filled');
plot(ax,nloc(k1n,1),nloc(k1n,2),'k');
plot(ax,xell,yell,'k-','linew',1);
plot(ax,xrec,yrec,'k--','linew',1);
axis equal
xlim(ax,xgrdran);
ylim(ax,ygrdran);
title(ax,'Gradient vector, slab1.0.grid');
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
hold(ax,'off');

moveh(ax,-0.05)

%%
%Absoulte magnitude of the gradient 
magogrd = sqrt(dzdxo.^2+dzdyo.^2);
magomesh = reshape(magogrd,[],1);

magngrd = sqrt(dzdxn.^2+dzdyn.^2);
magnmesh = reshape(magngrd,[],1);


figure
ax=subplot(1,2,1);
scatter(ax,xmesh,ymesh,40,magomesh,'s','filled'); hold(ax,'on');
scatter(ax,oloc0(:,1),oloc0(:,2),30,'k','filled');
plot(ax,oloc(k1o,1),oloc(k1o,2),'k');
plot(ax,xell,yell,'k-','linew',1);
plot(ax,xrec,yrec,'k--','linew',1);
plot(ax,xprof,yprof,'r-','linew',2);
colormap('jet');
c=colorbar;
c.Label.String = strcat('Magnitude of spatial gradient');
axis equal
xlim(ax,xgrdran);
ylim(ax,ygrdran);
title(ax,'Gradient magnitude, xyz.grid');
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
hold(ax,'off');

ax=subplot(1,2,2);
scatter(ax,xmesh,ymesh,40,magnmesh,'s','filled'); hold(ax,'on');
scatter(ax,nloc0(:,1),nloc0(:,2),30,'k','filled');
plot(ax,nloc(k1n,1),nloc(k1n,2),'k');
plot(ax,xell,yell,'k-','linew',1);
plot(ax,xrec,yrec,'k--','linew',1);
plot(ax,xprof,yprof,'r-','linew',2);
colormap('jet');
c=colorbar;
c.Label.String = strcat('Magnitude of spatial gradient');
axis equal
xlim(ax,xgrdran);
ylim(ax,ygrdran);
title(ax,'Gradient magnitude, slab1.0.grid');
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
hold(ax,'off');


%% cross-section profile
%interpolate the depth on the cross-section profile
zoprof = (Fo(xprof,yprof))';
znprof = (Fn(xprof,yprof))';

%the cross-section passes 0,0, so you have multiple ways to get a relative distance measure along
%the the profile. For instance, 'cart2pol', sqrt
[xprofrot,yprofrot] = coordinate_rot(xprof,yprof,rotang,[0,0]);

%get the 1-D derivative of depth relative to the distance along profile
dzoprof = dfxdx(zoprof, xprofrot);
dznprof = dfxdx(znprof, xprofrot);

%dip angle at each point
angdipo = atand(dzoprof);
angdipn = atand(dznprof);

%fit the profile with a single line, and get the dip angle from the slope
fttpfree = fittype( @(a,b,x) a*x+b);
[fitobj,~,~] = fit(xprofrot, zoprof,fttpfree,'Robust','Bisquare','StartPoint',[1 1]);
coef = coeffvalues(fitobj);
angdipoave = atand(coef(1));
[fitobj,~,~] = fit(xprofrot, znprof,fttpfree,'Robust','Bisquare','StartPoint',[1 1]);
coef = coeffvalues(fitobj);
angdipnave = atand(coef(1));

%get a rough sense on the range of the ellipse and rectangle along the profile 
[xellrot,yellrot] = coordinate_rot(xell,yell,rotang,[0,0]);
[xrecrot,yrecrot] = coordinate_rot(xrec,yrec,rotang,[0,0]);

%%%
figure
ax=subplot(2,1,1);
yyaxis(ax,'left'); ax.YAxis(2).Color
plot(ax,xprofrot,zoprof,'r','linew',1); hold(ax,'on');
% plot(ax,xprofrot(1:end-1),dzoprof+mean(cran),'b');
scatter(ax,oloc0(:,1),oloc0(:,5),30,'k','filled');
plot(ax,minmax(xellrot'),[cran(1)+4 cran(1)+4],'k-','linew',2);
plot(ax,minmax(xrecrot'),[cran(1)+2 cran(1)+2],'k--','linew',2);
text(ax,0.1,0.9,sprintf('Dip angle from linear fitting: %.1f deg',angdipoave),'Units','normalized');
axis equal
xlim(ax,[floor(min(xprofrot)) ceil(max(xprofrot))]);
ylim(ax,cran);
xlabel(ax,'Distance along profile SW-NE (km)');
ylabel(ax,'Depth (km)');
ax.YDir = 'reverse';
yyaxis(ax,'right');
plot(ax,xprofrot(1:end-1),angdipo,'k');
ylabel(ax,'Dip angle (deg)');
ylim(ax,[5 25]);
ax.YAxis(1).Color='r';
ax.YAxis(2).Color='k';
hold(ax,'off');

ax=subplot(2,1,2);
yyaxis(ax,'left');
plot(ax,xprofrot,znprof,'r','linew',1); hold(ax,'on');
% plot(ax,xprofrot(1:end-1),dznprof+mean(cran),'b');
scatter(ax,nloc0(:,1),nloc0(:,5),30,'k','filled');
plot(ax,minmax(xellrot'),[cran(1)+4 cran(1)+4],'k-','linew',2);
plot(ax,minmax(xrecrot'),[cran(1)+2 cran(1)+2],'k--','linew',2);
text(ax,0.1,0.9,sprintf('Dip angle from linear fitting: %.1f',angdipnave),'Units','normalized');
axis equal
xlim(ax,[floor(min(xprofrot)) ceil(max(xprofrot))]);
ylim(ax,cran);
xlabel(ax,'Distance along profile SW-NE (km)');
ylabel(ax,'Depth (km)');
ax.YDir = 'reverse';
yyaxis(ax,'right');
plot(ax,xprofrot(1:end-1),angdipn,'k');
ylabel(ax,'Dip angle (deg)');
ylim(ax,[5 25]);
ax.YAxis(1).Color='r';
ax.YAxis(2).Color='k';
hold(ax,'off');


%% Now let's do the coordinate rotation
%%%After rotation, what is a horizontal rectangle in X-Y-Z system like in X''-Y''-Z''
%1, Rotate E(X) and N(Y), U(Z) to the strike orthogonal and strike direction, done via rotating
%counter-clockwise by an angle of 90-strikeort by Z axis with the bottom left corner as the rotation
%center
rotcnt = [-6 -3 0];
[xr,yr,zr] = coordinate_rot3d(xrec,yrec,zeros(length(xrec),1),90-strikeort,3,rotcnt);
[x0r,y0r,z0r] = coordinate_rot3d(0,0,0,90-strikeort,3,rotcnt);

%2, Project the object on a plane that is parallel to the slab interface, which has a angle of dip 
%relative to the horizontal plane. Obtain the location of projection in X'-Y'-Z' system. You don't 
%have to project to the real slab, as the difference is a constant, the depth of bottom left corner. 
xrp = xr;
yrp = yr;
zrppln = -xrp.*tand(angdipnave);
zrpint = -(Fn(xrec,yrec)-Fn(rotcnt(1),rotcnt(2)));
zrp = zrpint;

%3. Rotate the object on the plane counter-clockwise by the dip angle by Y' axis. Obtain the
%location in X'-Y'-Z' system.
[xrpr,zrpr] = complex_rot(xrp,zrp,angdipnave);
yrpr = yrp;

%4. Rotate the X'-Y'-Z' coordinate system back to X-Y-Z system with rotation center being the 
%current loc in X'-Y'-Z' of (0,0,0) in X-Y-Z.
[xrprr,yrprr,zrprr] = coordinate_rot3d(xrpr,yrpr,zrpr,-(90-strikeort),3,[x0r,y0r,z0r]);

figure
%In original X-Y-Z coordinate
ax=subplot(2,2,1);
plot(ax,xrec,yrec,'k-','linew',1); hold(ax,'on');
scatter(ax,0,0,15,'ko','filled'); % origin in X-Y-Z
scatter(ax,rotcnt(1),rotcnt(2),15,[0.4 0.4 0.4],'filled'); % rotation center in X-Y-Z
xran = [-10 10];
yran = [-10 10];
axis equal
[rotx, roty] = complex_rot(-2,-3,90-strikeort,rotcnt(1:2));
xvect = [rotcnt(1) rotx];
yvect = [rotcnt(2) roty];
drawArrow(ax,xvect,yvect,xran,yran,'linewidth',0.5,'color',[0.4 0.4 0.4]);  %down-dip direction
[rotx, roty] = complex_rot(-6,1,90-strikeort,rotcnt(1:2));
xvect = [rotcnt(1) rotx];
yvect = [rotcnt(2) roty];
drawArrow(ax,xvect,yvect,xran,yran,'linewidth',0.5,'color',[0.4 0.4 0.4]);  %strike direction
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
title(ax,'Original object');

%Intermediate rotation in X'-Y'-Z' system.
ax=subplot(2,2,2);
plot(ax,xr,yr,'k-','linew',1); hold(ax,'on');
scatter(ax,x0r,y0r,15,'ko','filled'); % origin in X'-Y'-Z'
scatter(ax,0,0,15,[0.4 0.4 0.4],'filled'); % rotation center in X'-Y'-Z'
plot(ax,xrpr,yrpr,'r-','linew',1);  % projected object in X'-Y'-Z'
axis equal
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Down-dip (km)');
ylabel(ax,'Strike (km)');
title(ax,'Intermediate rotation');

%Comparison between projection after rotation and original object in X-Y-Z system
ax=subplot(2,2,3);
plot(ax,xrec,yrec,'k-','linew',1); hold(ax,'on');
scatter(ax,0,0,15,'ko','filled'); % origin in X-Y-Z
scatter(ax,rotcnt(1),rotcnt(2),15,[0.4 0.4 0.4],'filled'); % rotation center in X-Y-Z
plot(ax,xrprr,yrprr,'r-','linew',1);  % projected object in X-Y-Z
axis equal
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
title(ax,'Final comparison');
hold(ax,'off');

ax=subplot(2,2,4);
plot(ax,zrppln,'k-','linew',1); hold(ax,'on');
plot(ax,zrpint,'k--','linew',1);
ylabel(ax,'Relative depth (km)');
title(ax,'Depth comparison between plane approx. and truth');
hold(ax,'off');





