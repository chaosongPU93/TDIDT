function [lon,lat] = relaloc2absloc(dx,dy,lon0,lat0) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [lon,lat] = relaloc2absloc(dx,dy,lon0,lat0) 
% 
% This function is used to convert the relative location in km (dx, dy) wrt
% a reference location in longitude and latitude (lon0,lat0), to absolute
% location in longitude and latitude
% 
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2024/03/13
% Last modified date:   2024/03/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rads=pi/180.;
erad=6372.028;
srad=erad*cos(lat0*rads);
% N--S, only depends on lat
lat = dy/rads/erad+lat0;
% W--E, depends on both lon and lat; high lat leads to smaller dx
lon = dx/rads/srad+lon0;

% keyboard