function [projx,projy,projxy] = customprojection(locxy,projang)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [projx,projy,projxy] = CUSTOMPROJECTION(locxy,projang)
%
% Sometimes you might have the need to project the location coordinate 'locxy'
% of sources to along a custom direction 'projang'. Basically it is a coordinate 
% frame transformation, ie., get the location in a custom frame. Call
% the 'coordinate_rot' for actual transform.
% --Here 'projang' follows the geographic convention, ie, measure from N 
%   clockwise, range is [0 360]. Be sure to convert from four-quadrant 
%   convention which measures from E cclockwise to 180, cclockwise to -180.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/01/26
% Last modified date:   2023/01/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projxy = locxy;
for i = 1: size(locxy,1)
  [projxy(i,1), projxy(i,2)] = coordinate_rot(locxy(i,1),locxy(i,2),-(projang-90),[0 0]);
end
projx = projxy(:,1);
projy = projxy(:,2);