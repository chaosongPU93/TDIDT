% testhypogrid
%
% this script systematically tests the behavior of hypoinverse in different ways.
% 1. 
% --Plot the used slab model by John's package, compared to the inversion of a grid
%   directly from John. Get a difference in depth between what is interpolated based
%   on John's lon and lat using the slab model, and the inverted one in John's grid.
% --The lesson is, John's is also not perfect, but very regular 
% 2.
% --Invert a grid of offset using John's package, with a fixed depth of 35 for every 
%   source, for 10 iteration; Then try fix depth to 37 km; Furthermore, try 20 
%   iterations.
% --The lesson is, result is slightly getting more regularly slightly, indicating that
%   the whole package is doing fine, but the starting depth is not good for all 
%   sources, so 10 iterations is not enough to calibrate locations for all. 
% 3.
% --Decrease the sampling rate to 20 Hz, try 2 difference sizes of grid, weird things
%   still occur, indicating that the intial depth is not proper. Then we tried to 
%   invert the same grid as John's, of course it is worse. And it can not be fixed by
%   simply add more iterations. But we can take the depth of each source from the 
%   raw result as the starting depth for each source and do a refinement. As expected,
%   the result is improved significantly. Then it is obvious to repeat this refinement
%   a few times, and then compare John's grid and my refined grid, and also look at
%   the depth difference as done in step 1.
% --The lesson is, this method works. An approprite regime is, for the 1st run of 
%   Hypoinv, fix depth to 37 km, run 20 iterations, and then run 20 rounds of 
%   refinement, each of which has 20 iterations.
%   More iterations that is rather meaningless.
% 4. 
% --Based on the findings through benchmarking in 1-3, we can now invert our own grids.
%   The lesson is, you are demanding too much given the current resolution of the grid,
%   so that there is no obvious improvement in the locations of sources. The OLD grid
%   model has a comparable spatial resolution to 2-sample spacing at 40 hz, meaning that
%   temporal resolution higher than 0.05 s (1 smaple at 20 Hz) might give you random 
%   error in location because of interpolation.
% --It also tells you that, the larger the grid, the larger inaccuracy in locations
%   using a initial fixed depth to all sources, so that you need more rounds of refinement. 
% --All results in steps 1-4 are in /LOCpgsssitest.
% 5.
% --Now let's try our new slab grid, and switch folder to /forsummary. Its turns out
%   that the slab1.0.grid in /forsummary has the same spatial resolution as the old one,
%   it just that may be the depth is more accurate. But the max spatial resolution 
%   determines that the inversion result won't be very different from the old slab,
%   especially in lat and lon, but the depth changes due to the slab depth variation.
% --It also tells you that for this smaller grid, a fixed depth to all sources are not
%   too bad, so you don't need so many rounds of refinements
% 6. 
% --Create a denser grid of slab1.0 by interpolation > 0.01 deg in lat or lon, see
%   'edit_slabgridv2.m' for details. For example, we create a 0.002*0.002 deg of slab
%   model from spline interpolation, and create 'iterchaov2.c' to read slab model and
%   find the depth of source from the lat and lon from hypoinverse. 
% --It turns out that there is almost no improvement. The reason behind is, the resolution
%   of depth is only accurate to 0.01 km. A higher resolution in lat and lon cannot cause
%   difference in depth larger than 0.01 km, so the locations are not becoming more regular.
% --Naturally, you can try to increase the resolution of depth to 0.001 km in 'chat00p'
%   and output in 'iterchaov2.c', but hypoinverse only takes in the first 2 digits of 
%   decimals in depth for inversion, so trying that does NOT help. There maybe somewhere
%   to change the parameter setting of hypo, but hard to know without the manual.
% 7.
% --Now the only hope is, can we obtain a smoother location grid using the better slab 
%   model, either the 0.01*0.01 deg slab1.0 (sam resolution as the old one), or the 
%   0.002*0.002 deg slab1.0, from the same 2-sample spacing at 40 hz (0.05 s, 1 smaple
%   at 20 Hz) as John was using?
% --Hard to say if it is better, but definitely the new slab model is smoother in depth,
%   but depth is something we care the least, nevertheless, I think it is better to use
%   the new slab model. By the way, 0.002*0.002 deg slab does not help make the locations
%   more regular. 
% --At this point, I think the best strategy now is that, 
%   either use John's grid 'p085p020' to interpolate to any higher resolution than 2 
%   samples at 40 Hz, 
%   or, use new slab model 'slab1.0.grid' to do an inversion to 2-sample spacing grids
%   and add 5 rounds of extra refinement, take which as the starting point of 
%   further interpolation. A direct inversion of locations using the time offsets in 
%   samples finer than 20 sps is NOT recommended.
% 8.
% --Through trial and error, it looks like +-13 samples at 20 sps is the maximum range that
%   won't cause any wrapping in locations (feel like some kind of cycle-skipping), which 
%   is clear from the correspondance of corners of the grid.
% --So now, when the transformation between the time offset to spatial location is needed,
%   We can either use John's grid 'p085p020' to interpolate to any higher resolution than 
%   2 samples at 40 Hz,
%   Or use the inverted grid of +-13 samples at 20 sps using the 0.01*0.01 deg slab1.0 
%   model for further interpolation. Detailed inversion is, fix the depth of each source
%   to 33 km, run 20 iterations. Then we fix the depth of each source to the depth 
%   computed from the raw inversion, then run 5 or 10 rounds of refinements. 
%   Note that a finer slab model like 0.002*0.002 deg is a overkill, so are more iterations
%   and more refinements.
%   
% 
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/02/27
% Last modified date:   2022/03/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clc;
clear;
close all;

%% STEP 1
%old slab grid file from John Armbruster
olds = load('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssi/xyz.grid');

%Inverted locations from John Armbruster from a grid of offsets
A=load(strcat(getenv('ALLAN'), '/matfils/','p085p020'));
mm=25; nn=25;
dlat=reshape(A(:,3),mm,nn); %degrees
mlat=reshape(A(:,4),mm,nn); %minutes
lato=dlat+mlat/60.;
dlon=reshape(A(:,5),mm,nn); %degrees
mlon=reshape(A(:,6),mm,nn); %minutes
lono=-(dlon+mlon/60.);
depo=reshape(A(:,7),mm,nn);
off12o=reshape(A(:,1),mm,nn); %PGSS is uniform in columns; changes across rows.
off13o=reshape(A(:,2),mm,nn); %PGSI is uniform in rows; changes down columns.
% x2 = floor(off12o(1,1)): -1: ceil(off12o(1,end));   % avoid extrapolation
% y2 = floor(off13o(1,1)): -1: ceil(off13o(end,1));
% [off12n, off13n] = meshgrid(x2,y2); % new denser grid
lon = reshape(lono,[],1);
lat = reshape(lato,[],1);
dep = reshape(depo,[],1);

%this is the time offset 0,0 location relative to fam 002 from Allan
lat0=48.0+26.32/60.;    % lat of the event inverted from (0,0), off12,off13
lon0=-(123.0+35.07/60.);   % lon of the event inverted from (0,0), off12,off13

%use the old grid file to interpolate the location of the origin
F = scatteredInterpolant(-olds(:,1),olds(:,2),olds(:,3),'linear','linear');
dep0si = F(lon0,lat0);

%use the old grid file to interpolate the depth for John's locations
deposi = F(lono,lato);
depsi = reshape(deposi,[],1);

%find the location of the orgin by interpolating the existed grid
lat0i = interp2(off12o,off13o,lato,86,20,'linear');   % do each separately
lon0i = interp2(off12o,off13o,lono,86,20,'linear');
dep0i = interp2(off12o,off13o,depo,86,20,'linear');


figure
%what is the slab model like?
ax=subplot(2,2,1);
scatter(ax,-olds(:,1),olds(:,2),20,olds(:,3),'filled'); hold(ax,'on');
c=colorbar;
c.Label.String = strcat('Depth to interface');
% cran = c.Limits;
cran = [30 40];
caxis(ax,cran);

[k1,v1] = boundary(lon,lat,0.5);
% scatter(ax,lon(k1),lat(k1),20,'ko');
plot(ax,lon(k1),lat(k1),'k');
% scatter(ax,lon0,lat0,20,'ro','filled');
yran = [48.2 48.7];
xran = [-124 -123];
xlim(ax,xran);
ylim(ax,yran);
title(ax,'Slab model xyz.grid in /LOCpgsssi');
xlabel(ax,'Lon');
ylabel(ax,'Lat');
hold(ax,'off');

%what is the John's grid of inversions like?
ax=subplot(2,2,2);
scatter(ax,lon,lat,20,dep,'filled'); hold(ax,'on');
plot(ax,lon(k1),lat(k1),'k');
% scatter(ax,lon0,lat0,20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(ax,cran);
xlim(ax,xran);
ylim(ax,yran);
title(ax,'Grid file p085p020');
xlabel(ax,'Lon');
ylabel(ax,'Lat');
hold(ax,'off');

%if John's inversion of lon and lat are correct, is depth correct? ie, are sources indeed on the
%slab interface? Obtain the depth difference
ax=subplot(2,2,3);
scatter(ax,lon,lat,20,dep-depsi,'filled'); hold(ax,'on');
% scatter(ax,lon0,lat0,20,'ro','filled');
colormap(ax,'redblue');
% oldc = colormap(ax,'kelicol');
% newc = flipud(oldc);
% colormap(ax,newc);
c=colorbar;
c.Label.String = strcat('Depth difference');
caxis(ax,[-max(abs(dep-depsi)) max(abs(dep-depsi))]);
xlim(ax,xran);
ylim(ax,yran);
title(ax,"Depth difference between true slab and p085p020");
xlabel(ax,'Lon');
ylabel(ax,'Lat');
hold(ax,'off');

depdifstd = std(dep-depsi);
disp(depdifstd);
depdifmed = median(dep-depsi);
disp(depdifmed);

%%%By checking, John's slab file is a rectangle, while the lon decreases first while the lat stays 
%%%the same, then the lat decreases by one grid point. Basically the order is tracking each row of
%%%the 'grid', which includes 401 longitudes and 301 latitudes
nlon = length(unique(olds(:,1)));
nlat = length(unique(olds(:,2)));
if nlon*nlat ~= size(olds,1)
  disp('The slab is not a rectangle grid');
end

keyboard

%% STEP 2
%inverion results using LOCpgsssi by fixing the initial trial depth to 35 km
loc = load('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/evtloc.offset_002_rectgrid_6_40.Armb');

%similar to loc, but fix the initial trial depth to 37 km
loc1 = load('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/evtloc.offset_002_rectgrid_6_40.Armb1');

%similar to loc1, but undergoes 20 iterations instead of 10
loc2 = load('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/evtloc.offset_002_rectgrid_6_40.Armb2');

figure
subplot(2,2,1)
scatter(loc(:,1),loc(:,2),20,loc(:,3),'filled'); hold on
% scatter(loc(85,1),loc(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis([36.3 38.3]);
xlabel('Lon');
ylabel('Lat');
title('fixed depth to 35 km, 10 iterations');

subplot(2,2,2)
scatter(loc1(:,1),loc1(:,2),20,loc1(:,3),'filled'); hold on
% scatter(loc1(85,1),loc1(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis([36.3 38.3]);
xlabel('Lon');
ylabel('Lat');
title('fixed depth to 37 km, 10 iterations');

subplot(2,2,3)
scatter(loc2(:,1),loc2(:,2),20,loc2(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis([36.3 38.3]);
xlabel('Lon');
ylabel('Lat');
title('fixed depth to 37 km, 20 iterations');

keyboard

%% STEP 3
%old slab grid file from John Armbruster
olds = load('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssi/xyz.grid');
%use the old grid file to interpolate the location of the origin
F = scatteredInterpolant(-olds(:,1),olds(:,2),olds(:,3),'linear','linear');

%fix the initial trial depth to 37 km, 20 iterations, +-3 samples at 20 sps
loc3 = load('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/evtloc.offset_002_rectgrid_3_20.Armb2');

%fix the initial trial depth to 37 km, 20 iterations, +-20 samples at 20 sps
loc20 = load('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/evtloc.offset_002_rectgrid_20_20.Armb2');

%John's grid, -23:2:25 in off12, 24:2:24 in off13 at 40 sps
locarmb = load(strcat(getenv('ALLAN'), '/matfils/','p085p020'));
dlat=locarmb(:,3); %degrees
mlat=locarmb(:,4); %minutes
lato=dlat+mlat/60.;
dlon=locarmb(:,5); %degrees
mlon=locarmb(:,6); %minutes
lono=-(dlon+mlon/60.);
depo=locarmb(:,7);

%fix the initial trial depth to 37 km, 20 iterations, -23:2:25 in off12, 24:2:24 in off13 at 40 sps
loc12 = load('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/evtloc.offset_002_rectgrid_armb_40.Armb2');

%fix the initial trial depth to 37 km, 40 iterations, -23:2:25 in off12, 24:2:24 in off13 at 40 sps
loc123 = load('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/evtloc.offset_002_rectgrid_armb_40.Armb3');

%fixed depth to 1st-round inversion, 20 iterations, -23:2:25 in off12, 24:2:24 in off13 at 40 sps
loc12ref1 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/',...
  'evtloc.offset_002_rectgrid_armb_40.Armb2.ref1'));

%fixed depth to 1st-round refinement, 20 iterations, -23:2:25 in off12, 24:2:24 in off13 at 40 sps
loc12ref2 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/',...
  'evtloc.offset_002_rectgrid_armb_40.Armb2.ref2'));

%fixed depth to 2nd-round refinement, 20 iterations, -23:2:25 in off12, 24:2:24 in off13 at 40 sps
loc12ref3 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/',...
  'evtloc.offset_002_rectgrid_armb_40.Armb2.ref3'));

%fixed depth to 3rd-round refinement, 20 iterations, -23:2:25 in off12, 24:2:24 in off13 at 40 sps
loc12ref4 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/',...
  'evtloc.offset_002_rectgrid_armb_40.Armb2.ref4'));

%fixed depth to 4th-round refinement, 20 iterations, -23:2:25 in off12, 24:2:24 in off13 at 40 sps
loc12ref5 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/',...
  'evtloc.offset_002_rectgrid_armb_40.Armb2.ref5'));

%fixed depth to 19th-round refinement, 20 iterations, -23:2:25 in off12, 24:2:24 in off13 at 40 sps
loc12ref20 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/',...
  'evtloc.offset_002_rectgrid_armb_40.Armb2.ref20'));

figure
%what is the slab model like?
subplot(3,4,1);
scatter(-olds(:,1),olds(:,2),20,olds(:,3),'filled');  hold on
c=colorbar;
c.Label.String = strcat('Depth to interface');
% cran = c.Limits;
caxis([36.3 38.3]);
xran=[-124 -123.4];
yran=[48.22 48.56];
xlim(xran);
ylim(yran);
title('Slab model xyz.grid in /LOCpgsssi');
xlabel('Lon');
ylabel('Lat');

subplot(3,4,2)
scatter(loc20(:,1),loc20(:,2),20,loc20(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis([36.3 38.3]);
xlabel('Lon');
ylabel('Lat');
title('+/-20 samples at 20 Hz, fixed depth to 37 km, 20 iterations');

subplot(3,4,3);
scatter(lono,lato,20,depo,'filled'); hold on
% scatter(ax,lon0,lat0,20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis([36.3 38.3]);
xran=[-124 -123.4];
yran=[48.22 48.56];
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('Grid file p085p020');

subplot(3,4,4)
scatter(loc12(:,1),loc12(:,2),20,loc12(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis([36.3 38.3]);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('Same grid as p085p020, fixed depth to 37 km, 20 iterations');

subplot(3,4,5)
scatter(loc123(:,1),loc123(:,2),20,loc123(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis([36.3 38.3]);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('fixed depth to 37 km, 40 iterations');

subplot(3,4,6)
scatter(loc12ref1(:,1),loc12ref1(:,2),20,loc12ref1(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis([36.3 38.3]);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('fixed depth to 1st-round inversion, 20 iterations');

subplot(3,4,7)
scatter(loc12ref2(:,1),loc12ref2(:,2),20,loc12ref2(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis([36.3 38.3]);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('fixed depth to 1st-round refinement');

subplot(3,4,8)
scatter(loc12ref3(:,1),loc12ref3(:,2),20,loc12ref3(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis([36.3 38.3]);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('fixed depth to 2nd-round refinement');

subplot(3,4,9)
scatter(loc12ref4(:,1),loc12ref4(:,2),20,loc12ref4(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis([36.3 38.3]);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('fixed depth to 3rd-round refinement');

subplot(3,4,10)
scatter(loc12ref5(:,1),loc12ref5(:,2),20,loc12ref5(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis([36.3 38.3]);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('fixed depth to 4th-round refinement');

subplot(3,4,11)
scatter(loc12ref20(:,1),loc12ref20(:,2),20,loc12ref20(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis([36.3 38.3]);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('fixed depth to 19th-round refinement');

ax=subplot(3,4,12);
deprefsi = F(loc12ref20(:,1),loc12ref20(:,2));
scatter(ax,loc12ref20(:,1),loc12ref20(:,2),20,loc12ref20(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
caxis(ax,[-max(abs(loc12ref20(:,3)-deprefsi)) max(abs(loc12ref20(:,3)-deprefsi))]);
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference between true slab and newest inversion");

depdifstd = std(loc12ref20(:,3)-deprefsi);
disp(depdifstd);
depdifmed = median(loc12ref20(:,3)-deprefsi);
disp(depdifmed);

%%%
figure; histogram(dep-depsi,'BinWidth',0.002); hold on;
histogram(loc12ref20(:,3)-deprefsi,'BinWidth',0.002);
xlim([-0.2 0.2]);

keyboard

%% STEP 4
%fix the initial trial depth to 37 km, 20 iterations, +-15 samples in off12 & off13 at 100 sps
locraw = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/',...
  'evtloc.offset_002_rectgrid_15_100.Armb2'));

%5th round of refinement based on the raw locations
locref5 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/',...
  'evtloc.offset_002_rectgrid_15_100.Armb2.ref5'));

%10th round of refinement based on the raw locations
locref10 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/',...
  'evtloc.offset_002_rectgrid_15_100.Armb2.ref10'));

%15th round of refinement based on the raw locations
locref15 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/',...
  'evtloc.offset_002_rectgrid_15_100.Armb2.ref15'));

%20th round of refinement based on the raw locations
locref20 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/',...
  'evtloc.offset_002_rectgrid_15_100.Armb2.ref20'));

figure
%what is the slab model like?
subplot(3,4,1);
scatter(-olds(:,1),olds(:,2),20,olds(:,3),'filled'); hold on
c=colorbar;
c.Label.String = strcat('Depth to interface');
% cran = c.Limits;
cran = [36 38];
caxis(cran);
yran = [48.38 48.48];
xran = [-123.64 -123.54];
xlim(xran);
ylim(yran);
title('Slab model xyz.grid in /LOCpgsssi');
xlabel('Lon');
ylabel('Lat');

subplot(3,4,3)
scatter(locraw(:,1),locraw(:,2),15,locraw(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('Raw inversion, fixed depth to 37 km, 20 iterations');

ax=subplot(3,4,4);
deprefsi = F(locraw(:,1),locraw(:,2));
scatter(ax,locraw(:,1),locraw(:,2),15,locraw(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
caxis(ax,[-0.2 0.2]);
depdifmed = median(locraw(:,3)-deprefsi);
depdifstd = std(locraw(:,3)-deprefsi);
text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference between true slab and Raw inversion");


subplot(3,4,5)
scatter(locref5(:,1),locref5(:,2),15,locref5(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('5th refinement');

ax=subplot(3,4,6);
deprefsi = F(locref5(:,1),locref5(:,2));
scatter(ax,locref5(:,1),locref5(:,2),15,locref5(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
% caxis(ax,[-max(abs(locref5(:,3)-deprefsi)) max(abs(locref5(:,3)-deprefsi))]);
caxis(ax,[-0.2 0.2]);
depdifmed = median(locref5(:,3)-deprefsi);
depdifstd = std(locref5(:,3)-deprefsi);
text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference");


subplot(3,4,7)
scatter(locref10(:,1),locref10(:,2),15,locref10(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('10th refinement');

ax=subplot(3,4,8);
deprefsi = F(locref10(:,1),locref10(:,2));
scatter(ax,locref10(:,1),locref10(:,2),15,locref10(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
% caxis(ax,[-max(abs(locref5(:,3)-deprefsi)) max(abs(locref5(:,3)-deprefsi))]);
caxis(ax,[-0.2 0.2]);
depdifmed = median(locref10(:,3)-deprefsi);
depdifstd = std(locref10(:,3)-deprefsi);
text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference");

subplot(3,4,9)
scatter(locref15(:,1),locref15(:,2),15,locref15(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('15th refinement');

ax=subplot(3,4,10);
deprefsi = F(locref15(:,1),locref15(:,2));
scatter(ax,locref15(:,1),locref15(:,2),15,locref15(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
% caxis(ax,[-max(abs(locref5(:,3)-deprefsi)) max(abs(locref5(:,3)-deprefsi))]);
caxis(ax,[-0.2 0.2]);
depdifmed = median(locref15(:,3)-deprefsi);
depdifstd = std(locref15(:,3)-deprefsi);
text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference");

subplot(3,4,11)
scatter(locref20(:,1),locref20(:,2),15,locref20(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('20th refinement');

ax=subplot(3,4,12);
deprefsi = F(locref20(:,1),locref20(:,2));
scatter(ax,locref20(:,1),locref20(:,2),15,locref20(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
% caxis(ax,[-max(abs(locref5(:,3)-deprefsi)) max(abs(locref5(:,3)-deprefsi))]);
caxis(ax,[-0.2 0.2]);
depdifmed = median(locref20(:,3)-deprefsi);
depdifstd = std(locref20(:,3)-deprefsi);
text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference");

keyboard

%% STEP 5
%new slab grid file slab1.0 with a resolution of 0.01*0.01 deg
news = load('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/slab1.0.grid');
Fn = scatteredInterpolant(-news(:,1),news(:,2),news(:,3),'linear','linear');

%fix the initial trial depth to 37 km, 20 iterations, +-15 samples in off12 & off13 at 100 sps
locraw = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/',...
  'evtloc.offset_002_rectgrid_15_100.mit'));

%5th round of refinement based on the raw locations
locref5 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/',...
  'evtloc.offset_002_rectgrid_15_100.mit.ref5'));

%10th round of refinement based on the raw locations
locref10 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/',...
  'evtloc.offset_002_rectgrid_15_100.mit.ref10'));

%15th round of refinement based on the raw locations
locref15 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/',...
  'evtloc.offset_002_rectgrid_15_100.mit.ref15'));

%20th round of refinement based on the raw locations
locref20 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/',...
  'evtloc.offset_002_rectgrid_15_100.mit.ref20'));

figure
%what is the slab model like?
subplot(3,4,1);
scatter(-news(:,1),news(:,2),20,news(:,3),'filled'); hold on
c=colorbar;
c.Label.String = strcat('Depth to interface');
% cran = c.Limits;
cran = [36 38];
caxis(cran);
yran = [48.38 48.48];
xran = [-123.64 -123.54];
xlim(xran);
ylim(yran);
title('New slab model slab1.0.grid in /forsummary');
xlabel('Lon');
ylabel('Lat');

subplot(3,4,3)
scatter(locraw(:,1),locraw(:,2),15,locraw(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('Raw inversion, fixed depth to 37 km, 20 iterations');

ax=subplot(3,4,4);
deprefsi = Fn(locraw(:,1),locraw(:,2));
scatter(ax,locraw(:,1),locraw(:,2),15,locraw(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
caxis(ax,[-0.2 0.2]);
depdifmed = median(locraw(:,3)-deprefsi);
depdifstd = std(locraw(:,3)-deprefsi);
text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference between true slab and Raw inversion");


subplot(3,4,5)
scatter(locref5(:,1),locref5(:,2),15,locref5(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('5th refinement');

ax=subplot(3,4,6);
deprefsi = Fn(locref5(:,1),locref5(:,2));
scatter(ax,locref5(:,1),locref5(:,2),15,locref5(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
% caxis(ax,[-max(abs(locref5(:,3)-deprefsi)) max(abs(locref5(:,3)-deprefsi))]);
caxis(ax,[-0.2 0.2]);
depdifmed = median(locref5(:,3)-deprefsi);
depdifstd = std(locref5(:,3)-deprefsi);
text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference");


subplot(3,4,7)
scatter(locref10(:,1),locref10(:,2),15,locref10(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('10th refinement, fixed depth to last refinement, 20 iterations');

ax=subplot(3,4,8);
deprefsi = Fn(locref10(:,1),locref10(:,2));
scatter(ax,locref10(:,1),locref10(:,2),15,locref10(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
% caxis(ax,[-max(abs(locref5(:,3)-deprefsi)) max(abs(locref5(:,3)-deprefsi))]);
caxis(ax,[-0.2 0.2]);
depdifmed = median(locref10(:,3)-deprefsi);
depdifstd = std(locref10(:,3)-deprefsi);
text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference");

subplot(3,4,9)
scatter(locref15(:,1),locref15(:,2),15,locref15(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('15th refinement');

ax=subplot(3,4,10);
deprefsi = Fn(locref15(:,1),locref15(:,2));
scatter(ax,locref15(:,1),locref15(:,2),15,locref15(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
% caxis(ax,[-max(abs(locref5(:,3)-deprefsi)) max(abs(locref5(:,3)-deprefsi))]);
caxis(ax,[-0.2 0.2]);
depdifmed = median(locref15(:,3)-deprefsi);
depdifstd = std(locref15(:,3)-deprefsi);
text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference");

subplot(3,4,11)
scatter(locref20(:,1),locref20(:,2),15,locref20(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('20th refinement');

ax=subplot(3,4,12);
deprefsi = Fn(locref20(:,1),locref20(:,2));
scatter(ax,locref20(:,1),locref20(:,2),15,locref20(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
% caxis(ax,[-max(abs(locref5(:,3)-deprefsi)) max(abs(locref5(:,3)-deprefsi))]);
caxis(ax,[-0.2 0.2]);
depdifmed = median(locref20(:,3)-deprefsi);
depdifstd = std(locref20(:,3)-deprefsi);
text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference");

keyboard

%% STEP 6
%new slab grid file slab1.0 with a resolution of 0.002*0.002 deg
news = load('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/slab1.0_2e-3.grid');
Fn2 = scatteredInterpolant(-news(:,1),news(:,2),news(:,3),'linear','linear');

%fix the initial trial depth to 37 km, 20 iterations, +-15 samples in off12 & off13 at 100 sps
locraw = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/',...
  'evtloc.offset_002_rectgrid_15_100.mit'));

%5th round of refinement based on the raw locations
locref5 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/',...
  'evtloc.offset_002_rectgrid_15_100.mit.ref5'));

% %10th round of refinement based on the raw locations
% locref10 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/',...
%   'evtloc.offset_002_rectgrid_15_100.mit.ref10'));
% 
% %15th round of refinement based on the raw locations
% locref15 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/',...
%   'evtloc.offset_002_rectgrid_15_100.mit.ref15'));
% 
%20th round of refinement based on the raw locations
locref20 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/',...
  'evtloc.offset_002_rectgrid_15_100.mit.ref20'));

figure
%what is the slab model like?
subplot(3,4,1);
scatter(-news(:,1),news(:,2),20,news(:,3),'filled'); hold on
c=colorbar;
c.Label.String = strcat('Depth to interface');
% cran = c.Limits;
cran = [36 38];
caxis(cran);
yran = [48.38 48.48];
xran = [-123.64 -123.54];
xlim(xran);
ylim(yran);
title('New slab model slab1.0_2e-3.grid in /forsummary');
xlabel('Lon');
ylabel('Lat');

subplot(3,4,3)
scatter(locraw(:,1),locraw(:,2),15,locraw(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('Raw inversion, fixed depth to 37 km, 20 iterations');

clear new

ax=subplot(3,4,4);
deprefsi = Fn2(locraw(:,1),locraw(:,2));
scatter(ax,locraw(:,1),locraw(:,2),15,locraw(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
caxis(ax,[-0.2 0.2]);
depdifmed = median(locraw(:,3)-deprefsi);
depdifstd = std(locraw(:,3)-deprefsi);
text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference between true slab and Raw inversion");


subplot(3,4,5)
scatter(locref5(:,1),locref5(:,2),15,locref5(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('5th refinement');

ax=subplot(3,4,6);
deprefsi = Fn2(locref5(:,1),locref5(:,2));
scatter(ax,locref5(:,1),locref5(:,2),15,locref5(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
% caxis(ax,[-max(abs(locref5(:,3)-deprefsi)) max(abs(locref5(:,3)-deprefsi))]);
caxis(ax,[-0.2 0.2]);
depdifmed = median(locref5(:,3)-deprefsi);
depdifstd = std(locref5(:,3)-deprefsi);
text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference");

subplot(3,4,7)
scatter(locref20(:,1),locref20(:,2),15,locref20(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('20th refinement');

ax=subplot(3,4,8);
deprefsi = Fn2(locref20(:,1),locref20(:,2));
scatter(ax,locref20(:,1),locref20(:,2),15,locref20(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
% caxis(ax,[-max(abs(locref5(:,3)-deprefsi)) max(abs(locref5(:,3)-deprefsi))]);
caxis(ax,[-0.2 0.2]);
depdifmed = median(locref20(:,3)-deprefsi);
depdifstd = std(locref20(:,3)-deprefsi);
text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference");

keyboard

%% STEP 7
%new slab grid file slab1.0 with a resolution of 0.01*0.01 deg
news = load('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/slab1.0.grid');
% %new slab grid file slab1.0 with a resolution of 0.002*0.002 deg
% news = load('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/slab1.0_2e-3.grid');
Fn3 = scatteredInterpolant(-news(:,1),news(:,2),news(:,3),'linear','linear');

%fix the initial trial depth to 37 km, 20 iterations, same grid as Armb's grid 'p085p020'
locraw = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/',...
  'evtloc.offset_002_rectgrid_armb_40.mit'));

%5th round of refinement based on the raw locations
locref5 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/',...
  'evtloc.offset_002_rectgrid_armb_40.mit.ref5'));

%10th round of refinement based on the raw locations
locref10 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/',...
  'evtloc.offset_002_rectgrid_armb_40.mit.ref10'));

%15th round of refinement based on the raw locations
locref15 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/',...
  'evtloc.offset_002_rectgrid_armb_40.mit.ref15'));

%20th round of refinement based on the raw locations
locref20 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/',...
  'evtloc.offset_002_rectgrid_armb_40.mit.ref20'));

figure
%what is the slab model like?
subplot(3,4,1);
scatter(-news(:,1),news(:,2),15,news(:,3),'filled'); hold on
c=colorbar;
c.Label.String = strcat('Depth to interface');
% cran = c.Limits;
cran = [36 38];
caxis(cran);
xran=[-124 -123.4];
yran=[48.22 48.56];
xlim(xran);
ylim(yran);
title('New slab model slab1.0.grid in /forsummary');
xlabel('Lon');
ylabel('Lat');

subplot(3,4,3)
scatter(locraw(:,1),locraw(:,2),15,locraw(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('Raw inversion, fixed depth to 37 km, 20 iterations');

clear new

ax=subplot(3,4,4);
deprefsi = Fn3(locraw(:,1),locraw(:,2));
scatter(ax,locraw(:,1),locraw(:,2),15,locraw(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
caxis(ax,[-0.2 0.2]);
depdifmed = median(locraw(:,3)-deprefsi);
depdifstd = std(locraw(:,3)-deprefsi);
text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference between true slab and Raw inversion");

subplot(3,4,5)
scatter(locref5(:,1),locref5(:,2),15,locref5(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('5th refinement');

ax=subplot(3,4,6);
deprefsi = Fn3(locref5(:,1),locref5(:,2));
scatter(ax,locref5(:,1),locref5(:,2),15,locref5(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
% caxis(ax,[-max(abs(locref5(:,3)-deprefsi)) max(abs(locref5(:,3)-deprefsi))]);
caxis(ax,[-0.2 0.2]);
depdifmed = median(locref5(:,3)-deprefsi);
depdifstd = std(locref5(:,3)-deprefsi);
text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference");


subplot(3,4,7)
scatter(locref10(:,1),locref10(:,2),15,locref10(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('10th refinement');

ax=subplot(3,4,8);
deprefsi = Fn3(locref10(:,1),locref10(:,2));
scatter(ax,locref10(:,1),locref10(:,2),15,locref10(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
% caxis(ax,[-max(abs(locref5(:,3)-deprefsi)) max(abs(locref5(:,3)-deprefsi))]);
caxis(ax,[-0.2 0.2]);
depdifmed = median(locref10(:,3)-deprefsi);
depdifstd = std(locref10(:,3)-deprefsi);
text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference");

subplot(3,4,9)
scatter(locref15(:,1),locref15(:,2),15,locref15(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('15th refinement');

ax=subplot(3,4,10);
deprefsi = Fn3(locref15(:,1),locref15(:,2));
scatter(ax,locref15(:,1),locref15(:,2),15,locref15(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
% caxis(ax,[-max(abs(locref5(:,3)-deprefsi)) max(abs(locref5(:,3)-deprefsi))]);
caxis(ax,[-0.2 0.2]);
depdifmed = median(locref15(:,3)-deprefsi);
depdifstd = std(locref15(:,3)-deprefsi);
text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference");

subplot(3,4,11)
scatter(locref20(:,1),locref20(:,2),15,locref20(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('20th refinement');

ax=subplot(3,4,12);
deprefsi = Fn3(locref20(:,1),locref20(:,2));
scatter(ax,locref20(:,1),locref20(:,2),15,locref20(:,3)-deprefsi,'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
% caxis(ax,[-max(abs(locref5(:,3)-deprefsi)) max(abs(locref5(:,3)-deprefsi))]);
caxis(ax,[-0.2 0.2]);
depdifmed = median(locref20(:,3)-deprefsi);
depdifstd = std(locref20(:,3)-deprefsi);
text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference");

keyboard

%% STEP 8
%new slab grid file slab1.0 with a resolution of 0.01*0.01 deg
news = load('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/slab1.0.grid');
% %new slab grid file slab1.0 with a resolution of 0.002*0.002 deg
% news = load('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/slab1.0_2e-3.grid');
Fn3 = scatteredInterpolant(-news(:,1),news(:,2),news(:,3),'linear','linear');

%fix the initial trial depth to 37 km, 20 iterations, same grid as Armb's grid 'p085p020'
locraw = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/',...
  'evtloc.offset_002_rectgrid_13_20.mit'));

%5th round of refinement based on the raw locations
locref5 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/',...
  'evtloc.offset_002_rectgrid_13_20.mit.ref5'));

%10th round of refinement based on the raw locations
locref10 = load(strcat('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/',...
  'evtloc.offset_002_rectgrid_13_20.mit.ref10'));

figure
%what is the slab model like?
subplot(3,4,1);
scatter(-news(:,1),news(:,2),8,news(:,3),'filled'); hold on
c=colorbar;
c.Label.String = strcat('Depth to interface');
% cran = c.Limits;
cran = [31 41];
caxis(cran);
xran=[-124.1 -123.3];
yran=[47.9 48.7];
xlim(xran);
ylim(yran);
title('New slab model slab1.0.grid in /forsummary');
xlabel('Lon');
ylabel('Lat');

subplot(3,4,3)
scatter(locraw(:,1),locraw(:,2),8,locraw(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(cran);
xlim(xran);
ylim(yran);
xlabel('Lon');
ylabel('Lat');
title('Raw inversion, fixed depth to 33 km, 20 iterations');

clear new

ax=subplot(3,4,4);
deprefsi = Fn3(locraw(:,1),locraw(:,2));
scatter(ax,locraw(:,1),locraw(:,2),8,locraw(:,3)-deprefsi,'filled'); hold on
off12c = 7;
off13c = 8;
ran12 = [min(locraw(:,off12c)) max(locraw(:,off12c))];
ran13 = [min(locraw(:,off13c)) max(locraw(:,off13c))];
indtr = find(locraw(:,off12c)==ran12(2) & locraw(:,off13c)==ran13(2));  % top right
scatter(ax,locraw(indtr,1),locraw(indtr,2),20,'r^','filled');
indbr = find(locraw(:,off12c)==ran12(2) & locraw(:,off13c)==ran13(1)); % bottom right
scatter(ax,locraw(indbr,1),locraw(indbr,2),20,'b^','filled');
indbl = find(locraw(:,off12c)==ran12(1) & locraw(:,off13c)==ran13(1));  % bottom left
scatter(ax,locraw(indbl,1),locraw(indbl,2),20,'k^','filled');
indtl = find(locraw(:,off12c)==ran12(1) & locraw(:,off13c)==ran13(2)); % top left
scatter(ax,locraw(indtl,1),locraw(indtl,2),20,'g^','filled');
indc = find(locraw(:,off12c)==0 & locraw(:,off13c)==0);  % center
scatter(ax,locraw(indc,1),locraw(indc,2),20,'kd','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
caxis(ax,[-0.1 0.1]);
depdifmed = median(locraw(:,3)-deprefsi);
depdifstd = std(locraw(:,3)-deprefsi);
% text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
% text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference between true slab and Raw inversion");

ax=subplot(3,4,5);
scatter(ax,locref5(:,1),locref5(:,2),8,locref5(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(ax,cran);
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,'5th refinement');

ax=subplot(3,4,6);
deprefsi = Fn3(locref5(:,1),locref5(:,2));
scatter(ax,locref5(:,1),locref5(:,2),8,locref5(:,3)-deprefsi,'filled'); hold on
scatter(ax,locref5(indtr,1),locref5(indtr,2),20,'r^','filled');
scatter(ax,locref5(indbr,1),locref5(indbr,2),20,'b^','filled');
scatter(ax,locref5(indbl,1),locref5(indbl,2),20,'k^','filled');
scatter(ax,locref5(indtl,1),locref5(indtl,2),20,'g^','filled');
scatter(ax,locref5(indc,1),locref5(indc,2),20,'kd','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
% caxis(ax,[-max(abs(locref5(:,3)-deprefsi)) max(abs(locref5(:,3)-deprefsi))]);
caxis(ax,[-0.1 0.1]);
depdifmed = median(locref5(:,3)-deprefsi);
depdifstd = std(locref5(:,3)-deprefsi);
% text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
% text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference");

ax=subplot(3,4,7);
scatter(ax,locref10(:,1),locref10(:,2),8,locref10(:,3),'filled'); hold on
% scatter(loc2(85,1),loc2(85,2),20,'ro','filled');
c=colorbar;
c.Label.String = strcat('Depth to interface');
caxis(ax,cran);
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,'10th refinement');

ax=subplot(3,4,8);
deprefsi = Fn3(locref10(:,1),locref10(:,2));
scatter(ax,locref10(:,1),locref10(:,2),8,locref10(:,3)-deprefsi,'filled'); hold on
scatter(ax,locref10(indtr,1),locref10(indtr,2),20,'r^','filled');
scatter(ax,locref10(indbr,1),locref10(indbr,2),20,'b^','filled');
scatter(ax,locref10(indbl,1),locref10(indbl,2),20,'k^','filled');
scatter(ax,locref10(indtl,1),locref10(indtl,2),20,'g^','filled');
scatter(ax,locref10(indc,1),locref10(indc,2),20,'kd','filled');
colormap(ax,'redblue');
c=colorbar;
c.Label.String = strcat('Depth difference');
% caxis(ax,[-max(abs(locref10(:,3)-deprefsi)) max(abs(locref10(:,3)-deprefsi))]);
caxis(ax,[-0.1 0.1]);
depdifmed = median(locref10(:,3)-deprefsi);
depdifstd = std(locref10(:,3)-deprefsi);
% text(0.8,0.2,sprintf('Med: %.4f',depdifmed),'Units','normalized','HorizontalAlignment','right');
% text(0.8,0.1,sprintf('Std: %.4f',depdifstd),'Units','normalized','HorizontalAlignment','right');
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Lon');
ylabel(ax,'Lat');
title(ax,"Depth difference");

