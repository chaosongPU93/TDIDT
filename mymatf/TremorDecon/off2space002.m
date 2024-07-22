function [locuse, indinput] = off2space002(offset,sps,ftrans,fplt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [locuse, indinput] = off2space002(offset,sps,ftrans,fplt)
%
% This function aims to quickly transform the travel time offset between 
% station 1 and 2, off12, and the offset between station 1 and 3, off13, 
% in samples (sps information needed) to locations on the slab interface,
% using the flag of the ways to transform. 
%
% --'offset' has 2 colomns for off12 and off13
% --You can either use the 'directhypo' to choose the inverted locations
%   using hypoinverse from a grid a offsets with several options of sps,
%   like 20, 40, 80, 100, then find a match. NOT recommended any more,
%   especially for the higher sampling rate than 20 sps.
% --You can also use the 'interpArmb' to choose a sparse grid of locations
%   inverted by John Armbruster (2-sample spacing at 40 sps) to interpolate
%   to the desired sampling rate, then find a match. The interpolation has
%   been benchmarked to be correct with that in 'grid_cat.m'. Note a bunch
%   of bad locations at the SW corner, could be problematic.
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
% First created date:   2022/02/16
% Last modified date:   2022/04/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% default value for easy debugging or no input
defval('offset',[]);
defval('sps',100);
defval('ftrans','interpchao');
defval('fplt',0);

%%%Find the exact match directly from the inversion using a grid of time offsets with the same
%%%sampling rate, NOT recommended any more!
if isequal(ftrans, 'directhypo')
  %set path
  gridpath = strcat(getenv('MHOME'), '/Seisbasics/hypoinverse/forsummary');
  if sps==20
    gridfile = 'evtloc.offset_002_rectgrid_3_20';
  elseif sps==40
    gridfile = 'evtloc.offset_002_rectgrid_6_40';
  elseif sps==80
    gridfile = 'evtloc.offset_002_rectgrid_12_80';
  elseif sps==100
    gridfile = 'evtloc.offset_002_rectgrid_15_100';
  else
    warning("The requested sps doesn't have a precomputed grid");
  end
  offgrid = load(fullfile(gridpath, gridfile));
  if isempty(offset)
    offset = offgrid(:,7:8);
  end
  [~, induniq] = unique(offset(:,1:2),'rows','stable');
  [~,indmatch,~] = intersect(offset(:,1:2),offgrid(:,7:8),'rows','stable');
  if length(indmatch) < length(induniq)
    warning("Some inputs did not find a match, used grid is too small");
    return
%   elseif length(indmatch) < length(induniq)  
%     warning("Some inputs find more than one match, needs check");
%     return
  end 
  for i =1:size(offset,1)
    tmp = find(offset(i,1)==offgrid(:,7) & offset(i,2)==offgrid(:,8));
    if length(tmp)>1
      disp(i);
      warning("this input offset has more than one match in the loc grid");
    else
      indgrid(i,1) = tmp;
      indinput(i,1) = i;
    end
  end
  tmp = offgrid(indgrid, :); % only use the matched points
  lon0 = -123.585000;
  lat0 = 48.436667;
  [dx,dy] = absloc2relaloc(tmp(:,1),tmp(:,2),lon0,lat0); % convert to rela loc around 002
  tmp = [dx dy tmp];
  locuse = [tmp(:, 1:6) tmp(:, 9:10)];  % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

  if fplt
    figure
    subplot(221)
    scatter(offgrid(:,1),offgrid(:,2),10,offgrid(:,1)+offgrid(:,2)); hold on
    scatter(locuse(:,3),locuse(:,4),10,'k');
    xlabel('lon'); ylabel('lat');
    subplot(222)
    scatter(offgrid(:,7),offgrid(:,8),10,offgrid(:,7)+offgrid(:,8)); hold on
    scatter(locuse(:,7),locuse(:,8),10,'k');
    xlabel('off12'); ylabel('off13');
    subplot(223)
    [dx,dy] = absloc2relaloc(offgrid(:,1),offgrid(:,2),lon0,lat0);
    scatter(dx,dy,10,dx+dy); hold on
    scatter(locuse(:,1),locuse(:,2),10,'k');
    xlabel('E (km)'); ylabel('N (km)');
  end

else
  %%%Interpolation to any resolution requested from John Armb's 2-sample spacing, 40-sps grid
  %%%'p085p020'.
  if isequal(ftrans, 'interpArmb')
    %set path
    gridpath = strcat(getenv('ALLAN'), '/matfils');
    gridfile = 'p085p020';  % refer to 'gridcat.m' from Allan
    offgrido = load(fullfile(gridpath, gridfile)); % 2-sample spacing at 40 Hz
    mm=25; nn=25;
    dlat=reshape(offgrido(:,3),mm,nn); %degrees
    mlat=reshape(offgrido(:,4),mm,nn); %minutes
    lat=dlat+mlat/60.;
    dlon=reshape(offgrido(:,5),mm,nn); %degrees
    mlon=reshape(offgrido(:,6),mm,nn); %minutes
    lon=-(dlon+mlon/60.);
    dep=reshape(offgrido(:,7),mm,nn);
    centPGSS=86;  % at 40 Hz
    centPGSI=20;
    spsscale = sps/40;  % scale to the requested sps
    off12o=spsscale*reshape(centPGSS-offgrido(:,1),mm,nn); %PGSS is uniform in columns; changes across rows.
    off13o=spsscale*reshape(centPGSI-offgrido(:,2),mm,nn); %PGSI is uniform in rows; changes down columns.
    x2 = floor(off12o(1,1)): -1: ceil(off12o(1,end));   % avoid extrapolation
    y2 = floor(off13o(1,1)): -1: ceil(off13o(end,1)); 
    [off12n, off13n] = meshgrid(x2,y2); % create a new denser grid

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
    x2 = floor(off12o(1,1)): -1: ceil(off12o(1,end));   % avoid extrapolation
    y2 = floor(off13o(1,1)): -1: ceil(off13o(end,1)); 
    [off12n, off13n] = meshgrid(x2,y2); % create a new denser grid  %%%Interpolation to any resolution requested from Chao's inverted +-13 samples, 20 sps grid with
  
  %%%Interpolation to any resolution requested from Chao's inverted +-13 samples, 20 sps grid with
  %%%Hypoinverse, using the 0.01*0.01-deg slab1.0
  elseif isequal(ftrans, 'interpchao')
    %set path
    gridpath = strcat(getenv('MHOME'), '/Seisbasics/hypoinverse/forsummary');
    gridfile = 'evtloc.offset_002_rectgrid_13_20.mit.ref5';
    offgrido = load(fullfile(gridpath, gridfile)); % 2-sample spacing at 40 Hz
    mm=27; nn=27;
    lon=reshape(offgrido(:,1),mm,nn); %degrees
    lat=reshape(offgrido(:,2),mm,nn); %degrees
    dep=reshape(offgrido(:,3),mm,nn); %km
    ttrvl=reshape(30-offgrido(:,4),mm,nn); %travel time, + origin time = 30 s, due to setting in Hypoinverse
    spsscale = sps/20;  % scale to the requested sps
    off12o=spsscale*reshape(offgrido(:,7),mm,nn); %PGSS is uniform in columns; changes across rows.
    off13o=spsscale*reshape(offgrido(:,8),mm,nn); %PGSI is uniform in rows; changes down columns.
    x2 = ceil(off12o(1,1)): 1: floor(off12o(1,end));   % avoid extrapolation
    y2 = ceil(off13o(1,1)): 1: floor(off13o(end,1)); 
    [off12n, off13n] = meshgrid(x2,y2); % create a new denser grid

  end
  
  %%%For interpolation, you can either do lat, lon or dep independently, or construct a complex
  %%%using the lon and lat as the real and imag parts
  %%%
  latn = interp2(off12o,off13o,lat,off12n,off13n,'spline');   % do each separately
  lonn = interp2(off12o,off13o,lon,off12n,off13n,'spline');
  depn = interp2(off12o,off13o,dep,off12n,off13n,'spline');
  if isequal(ftrans, 'interpArmb')  % John's original grid does not have origin time 
    ttrvl = zeros(size(lat));
    ttrvln = zeros(size(latn));
  else
    ttrvln = interp2(off12o,off13o,ttrvl,off12n,off13n,'spline');
  end
%   loc = complex(lon,lat); % construct a complex grid    % construct a complex of lon and lat
%   loc2 = interp2(off12o,off13o,loc,off12n,off13n,'spline');
%   lon2 = real(loc2);
%   lat2 = imag(loc2);

  %center, location of 002 
  lon0 = interp2(off12o,off13o,lon,0,0,'spline');  % interpolate as well
  lat0 = interp2(off12o,off13o,lat,0,0,'spline');
  
  %reshape to column vector
  lonn = reshape(lonn,[],1);
  latn = reshape(latn,[],1);
  depn = reshape(depn,[],1);
  ttrvln = reshape(ttrvln,[],1);
  off12n = reshape(off12n,[],1);
  off13n = reshape(off13n,[],1);
  offgrid = [lonn latn depn ttrvln off12n off13n];
  if isempty(offset)
    offset = offgrid(:,5:6);
  end
  [~, induniq] = unique(offset(:,1:2),'rows','stable');
  [~,indmatch,~] = intersect(offset(:,1:2),offgrid(:,5:6),'rows','stable');
  if length(indmatch) < length(induniq)
    warning("Some inputs did not find a match, used grid is too small");
    return
  end 
  for i =1:size(offset,1)
    tmp = find(offset(i,1)==offgrid(:,5) & offset(i,2)==offgrid(:,6));
    if length(tmp)>1
      disp(i);
      warning("this input offset has more than one match in the loc grid");
    else
      indgrid(i,1) = tmp;
      indinput(i,1) = i;
    end
  end
  tmp = offgrid(indgrid, :); % only use the matched points
  [dx,dy] = absloc2relaloc(tmp(:,1),tmp(:,2),lon0,lat0); % convert to rela loc around 002
  tmp = [dx dy tmp];
  locuse = tmp;   % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%   figure
%   subplot(411)
%   scatter(-lon,lat);
%   subplot(412)
%   scatter(-lon2,lat2);  
%   subplot(413)
%   scatter(-lon21,lat21);
%   subplot(414)
%   scatter(-lon22,lat22);

  if fplt
    figure
    subplot(221)
    scatter(reshape(lon,[],1),reshape(lat,[],1),10,reshape(ttrvl,[],1),'filled');
    xlabel('lon'); ylabel('lat'); colormap jet; colorbar;
    subplot(222)
    scatter(reshape(off12o,[],1),reshape(off13o,[],1),10,reshape(ttrvl,[],1),'filled');
    xlabel('off12'); ylabel('off13');
    subplot(223)
    scatter(offgrid(:,1),offgrid(:,2),10,offgrid(:,4),'filled'); hold on
    scatter(lon0,lat0,10,'r','filled');
    scatter(locuse(:,3),locuse(:,4),10,'k','filled');  
    xlabel('lon'); ylabel('lat'); colormap jet; colorbar;
    subplot(224)
    scatter(offgrid(:,5),offgrid(:,6),10,offgrid(:,4),'filled'); hold on
    scatter(0,0,10,'r','filled');
    scatter(locuse(:,7),locuse(:,8),10,'k','filled');
    xlabel('off12'); ylabel('off13');
%     subplot(325)
%     [dx,dy] = absloc2relaloc(offgrid(:,1),offgrid(:,2),lon0,lat0);
%     scatter(dx,dy,10,offgrid(:,4),'filled'); hold on
%     scatter(0,0,10,'r','filled');
%     scatter(locuse(:,1),locuse(:,2),10,'k','filled');
%     xlabel('E (km)'); ylabel('N (km)'); colormap jet; colorbar;
  end
  
end

% % Compare my hypoinverse with Armbruster's version using diff slab models
% ind = find(offgrid(:,4)>=-6 & offgrid(:,4)<=6 & ...
%   offgrid(:,5)>=-6 & offgrid(:,5)<=6);
% [dx,dy] = absloc2relaloc(offgrid(ind,1),offgrid(ind,2),-(123.0+35.07/60),48.0+26.32/60);
% offgridcmp1 = [dx dy offgrid(ind,:)];	% format, dx,dy,lon,lat,dep,off12,off13
% 
% gridpath = strcat(getenv('MHOME'), '/Seisbasics/hypoinverse/LOCpgsssitest');
% gridfile = 'evtloc.offset_002_rectgrid_6_40.Armb';
% offgrid = load(fullfile(gridpath, gridfile));
% [dx,dy] = absloc2relaloc(offgrid(:,1),offgrid(:,2),-123.588833,48.433833);
% offgrid = [dx dy offgrid];
% offgridcmp2 = [offgrid(:, 1:5) offgrid(:, 8:9)];  % format, dx,dy,lon,lat,dep,off12,off13
% offgridcmp2 = sortrows(offgridcmp2, [6 7],'descend');
% figure
% subplot(131)
% scatter(offgridcmp1(:,1),offgridcmp1(:,2));
% title('Direct interpolation from p085p020, accurate to 1 sample at 40 sps');
% axis equal
% axis([-6 4 -6 4]);
% subplot(132)
% scatter(offgridcmp2(:,1),offgridcmp2(:,2));
% title("Inversion with Armb's version of hypo");
% axis equal
% axis([-6 4 -6 4]);
% 
% gridpath = strcat(getenv('MHOME'), '/Seisbasics/hypoinverse/forsummary');
% gridfile = 'evtloc.offset_002_rectgrid_6_40';
% offgrid = load(fullfile(gridpath, gridfile));
% [dx,dy] = absloc2relaloc(offgrid(:,1),offgrid(:,2),-123.585000,48.436667);
% offgrid = [dx dy offgrid];
% offgridcmp3 = [offgrid(:, 1:5) offgrid(:, 8:9)];  % format, dx,dy,lon,lat,dep,off12,off13
% offgridcmp3 = sortrows(offgridcmp3, [6 7],'descend');
% subplot(133)
% scatter(offgridcmp3(:,1),offgridcmp3(:,2));
% title('Inversion with my version (diff. slab)');
% axis equal
% axis([-6 4 -6 4]);


% keyboard
