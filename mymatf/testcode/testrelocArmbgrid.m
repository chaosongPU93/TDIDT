% testrelocArmbgrid.m
%
% --John Armbruster's original catalog 'p085p020' does not give travel
% time or origin time, but hypoinverse does report the origin time
% given the arrival time at each station. 
% 
% --One way to do it is, take John's existing solution's depth as
% the fixed depth in hypo, and redo the inversion for a few rounds, ideally
% the result should not change too much relative to the input. But the new 
% result has the travel time information that can be used to interpolate at
% the current grid. This is done by '/hyposhell/pgc/relocArmbgrid.sh'
%
% --This test is trying to see how different the two grids are, for the same
% source, which is stored in the result file '/LOCpgsssitest/evtloc.p085p020'
%
% --looks like the agreement is not too bad.
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/04/14
% Last modified date:   2022/04/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%load the relocation result file
gridpath = strcat(getenv('MHOME'), '/Seisbasics/hypoinverse/LOCpgsssitest');
gridfile = 'evtloc.p085p020';
offgrido = load(fullfile(gridpath, gridfile)); % 2-sample spacing at 40 Hz

%relocation
lonn = offgrido(:,1);
latn = offgrido(:,2);
depn = offgrido(:,3);
ttrvln = 30-offgrido(:,4);

%original location
dlat = offgrido(:,7); %degrees
mlat = offgrido(:,8); %minutes
lato = dlat+mlat/60.;
dlon = offgrido(:,9); %degrees
mlon = offgrido(:,10); %minutes
lono = -(dlon+mlon/60.);
depo = offgrido(:,11);

lon0 = interp2(reshape(86-offgrido(:,5),25,25), reshape(20-offgrido(:,6),25,25), ...
  reshape(lono,25,25),0,0,'spline');  % interpolate as well
lat0 = interp2(reshape(86-offgrido(:,5),25,25), reshape(20-offgrido(:,6),25,25), ...
  reshape(lato,25,25),0,0,'spline');

%looks like the agreement is not too bad
figure
scatter(lon0,lat0,40,'b','filled'); hold on
scatter(lono,lato,20,'k','filled');
scatter(lonn,latn,6,'r','filled');

figure
scatter(lonn,latn,30,ttrvln,'filled');
colormap jet;
colorbar;

off2space002([86-offgrido(:,5) 20-offgrido(:,6)],40,'interpchao',1);

keyboard



















