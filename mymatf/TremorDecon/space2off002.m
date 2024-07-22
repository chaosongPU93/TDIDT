function offset = space2off002(loc,sps,ftrans,fplt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% offset = space2off002(loc,sps,ftrans,fplt)
%
% This function aims to quickly transform the locations on the slab
% interface loc=[dx, dy], into travel time offset between  station 1 and 2,
% off12, and the offset between station 1 and 3, off13, in samples
% (sps information needed) to, using the flag of the ways to transform.
% Essentially, it is a reverse transform of function 'off2space002' which
% tranforms offset to spatial locations. The key algorithm is interpolation 
% using 'griddata'.
%
% --'loc' has 2 colomns for dx and dy relative to eg. 002 family location
% --You can also use the 'interpArmb' to choose a sparse grid of locations
%   inverted by John Armbruster (2-sample spacing at 40 sps) to interpolate
%   at the desired locations to obtain the offset in samples at the assigned
%   sampling rate. You can open the plotting flag 'fplt' to validate the
%   result. Note a bunch of bad locations at the SW corner, could be
%   problematic if your desired locations are close to them.
% --The adopted way in 'grid_cat.m' is limited to interpolation to 2^n
%   times of the original spacing. Vq = interp2(V,K) returns the interpolated
%   values on a refined grid formed by repeatedly halving the intervals K
%   times in each dimension.This results in 2^K-1 interpolated points
%   between sample values.
% --Via the comprehensive tests in 'testhypogrid.m', we learned that the
%   max. spatial resolution reported by Hypoinverse (determined by 2-digit
%   decimal of depth) leads to a max temporal resolution about 1 sample at 20
%   sps, no matter how fine your slab model is (old slab with 0.01^2 deg^2;
%   slab1.0 with 0.01^2 deg^2; slab1.0 with 0.002^2 deg^2).
%   Therefore, now we offer a slightly better option: use the inverted grid
%   of +-13 samples at 20 sps using the 0.01*0.01 deg slab1.0 model for
%   further interpolation, called 'interpchao'.
% --This new option also recognize the fact that, because the source location
%   varies, the origin time for each source can NOT be viewed as the arrival
%   time, because the travel time changes. Therefore, the new location grid
%   file contains the orgin time estimate, which can be used to get a
%   calibration for a relative origin time to the center [0,0], or an absolute
%   origin time.
% --2022/04/14, Added a new option 'interpArmbreloc' which take John's existing 
%   solution's depth as the fixed depth in hypo, and redo the inversion for 
%   a few rounds, ideally the result should not change too much relative to 
%   the input.  This is done by '/hyposhell/pgc/relocArmbgrid.sh'. The result
%   is stored in file '/LOCpgsssitest/evtloc.p085p020'. In testing script
%   'testrelocArmbgrid.m', we can see that the agreement between the original 
%   grid 'p085p020' and the relocated grid is not too bad. The benefit of this
%   option is, the new result has the travel time information that can be 
%   used to interpolate at the current grid.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/04/12
% Last modified date:   2022/04/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% default value for easy debugging or no input
defval('loc',[]);
defval('sps',100);
defval('ftrans','interpArmb');
defval('fplt',1);

%%%Interpolation to any resolution requested from John Armb's 2-sample spacing, 40-sps grid
%%%'p085p020'.
if isequal(ftrans, 'interpArmb')
  %set path
  gridpath = strcat(getenv('ALLAN'), '/matfils');
  gridfile = 'p085p020';  % refer to 'gridcat.m' from Allan
  offgrido = load(fullfile(gridpath, gridfile)); % 2-sample spacing at 40 Hz
  dlat=offgrido(:,3); %degrees
  mlat=offgrido(:,4); %minutes
  lat=dlat+mlat/60.;
  dlon=offgrido(:,5); %degrees
  mlon=offgrido(:,6); %minutes
  lon=-(dlon+mlon/60.);
  dep=offgrido(:,7);
  centPGSS=86;  % at 40 Hz
  centPGSI=20;
  spsscale = sps/40;  % scale to the requested sps
  off12o=spsscale*(centPGSS-offgrido(:,1)); %PGSS is uniform in columns; changes across rows.
  off13o=spsscale*(centPGSI-offgrido(:,2)); %PGSI is uniform in rows; changes down columns.
  %center, location of 002
  lon0 = interp2(reshape(off12o,25,25), reshape(off13o,25,25), reshape(lon,25,25),0,0,'spline');  % interpolate as well
  lat0 = interp2(reshape(off12o,25,25), reshape(off13o,25,25), reshape(lat,25,25),0,0,'spline');  
  %convert to rela loc around 002
  [dx,dy] = absloc2relaloc(lon,lat,lon0,lat0);
  
%%%Interpolation to any resolution requested from the relocation of John Armb's 2-sample spacing,
%%%40-sps grid 'p085p020'.
elseif isequal(ftrans, 'interpArmbreloc')
  %set path
  gridpath = strcat(getenv('MHOME'), '/Seisbasics/hypoinverse/LOCpgsssitest');
  gridfile = 'evtloc.p085p020';
  offgrido = load(fullfile(gridpath, gridfile)); % 2-sample spacing at 40 Hz
  mm=25; nn=25;
  lon=reshape(offgrido(:,1),mm,nn); %degrees
  lat=reshape(offgrido(:,2),mm,nn); %degrees
  dep=reshape(offgrido(:,3),mm,nn); %km
  ttrvl=reshape(30-offgrido(:,4),mm,nn); %travel time, + origin time = 30 s, due to setting in Hypoinverse
  centPGSS=86;  % at 40 Hz
  centPGSI=20;
  spsscale = sps/40;  % scale to the requested sps
  off12o=spsscale*reshape(centPGSS-offgrido(:,5),mm,nn); %PGSS is uniform in columns; changes across rows.
  off13o=spsscale*reshape(centPGSI-offgrido(:,6),mm,nn); %PGSI is uniform in rows; changes down columns.
  %center, location of 002
  lon0 = interp2(reshape(off12o,25,25), reshape(off13o,25,25), reshape(lon,25,25),0,0,'spline');  % interpolate as well
  lat0 = interp2(reshape(off12o,25,25), reshape(off13o,25,25), reshape(lat,25,25),0,0,'spline');  
  %convert to rela loc around 002
  [dx,dy] = absloc2relaloc(lon,lat,lon0,lat0);
  
%%%Interpolation to any resolution requested from Chao's inverted +-13 samples, 20 sps grid with
%%%Hypoinverse, using the 0.01*0.01-deg slab1.0
elseif isequal(ftrans, 'interpchao')
  %set path
  gridpath = strcat(getenv('MHOME'), '/Seisbasics/hypoinverse/forsummary');
  gridfile = 'evtloc.offset_002_rectgrid_13_20.mit.ref5';
  offgrido = load(fullfile(gridpath, gridfile)); % 2-sample spacing at 40 Hz
  lon=offgrido(:,1); %degrees
  lat=offgrido(:,2); %degrees
  dep=offgrido(:,3); %km
  ttrvl=30-offgrido(:,4); %travel time, + origin time = 30 s, due to setting in Hypoinverse
  spsscale = sps/20;  % scale to the requested sps
  off12o=spsscale*offgrido(:,7); %PGSS is uniform in columns; changes across rows.
  off13o=spsscale*offgrido(:,8); %PGSI is uniform in rows; changes down columns.
  %center, location of 002
  lon0 = interp2(reshape(off12o,27,27), reshape(off13o,27,27), reshape(lon,27,27),0,0,'spline');  % interpolate as well
  lat0 = interp2(reshape(off12o,27,27), reshape(off13o,27,27), reshape(lat,27,27),0,0,'spline');
  %convert to rela loc around 002
  [dx,dy] = absloc2relaloc(lon,lat,lon0,lat0);
  
end

%use grid data to interpolate at requested locations
off12n = griddata(dx,dy,off12o,loc(:,1),loc(:,2),'cubic'); % do each separately
off13n = griddata(dx,dy,off13o,loc(:,1),loc(:,2),'cubic');
lonn = griddata(dx,dy,lon,loc(:,1),loc(:,2),'cubic');
latn = griddata(dx,dy,lat,loc(:,1),loc(:,2),'cubic');
depn = griddata(dx,dy,dep,loc(:,1),loc(:,2),'cubic');
if isequal(ftrans, 'interpArmb')  % John's original grid does not have origin time 
  ttrvl = zeros(size(lat));
  ttrvln = zeros(size(latn));
else
  ttrvln = griddata(dx,dy,ttrvl,loc(:,1),loc(:,2),'cubic');
end
  
% %8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
% locgrid = [loc(:,1) loc(:,2) lonn latn depn ttrvln off12n off13n];
%2 cols, format: off12,off13
offset = [off12n off13n]; 

if fplt
  color1 = ttrvl;
  color2 = ttrvln;
%   color1 = dep;
%   color2 = depn;
  figure
  subplot(321)
  scatter(dx,dy,10,color1,'filled');
  xlabel('E (km)'); ylabel('N (km)'); axis equal tight; box on
  subplot(322)
  scatter(off12o,off13o,10,color1,'filled');
  xlabel('off12'); ylabel('off13'); axis equal tight; box on
  subplot(323)
  scatter(locgrid(:,3),locgrid(:,4),10,color2,'filled'); hold on
  scatter(lon0,lat0,20,'r','filled');
  xlabel('lon'); ylabel('lat');
  subplot(324)
  scatter(locgrid(:,7),locgrid(:,8),10,color2,'filled'); hold on
  scatter(0,0,20,'r','filled');
  xlabel('off12'); ylabel('off13'); axis equal tight; box on
  subplot(325)
  scatter(locgrid(:,1),locgrid(:,2),10,color2,'filled'); hold on
  scatter(0,0,20,'r','filled');
  xlabel('E (km)'); ylabel('N (km)'); axis equal tight; box on
end

















