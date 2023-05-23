function [projx,orty,nlocxy] = customprojection(olocxy,projang)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [propx,orty,nlocxy] = customprojection(olocxy,propang)
%
% Sometimes you might have the need to project the location coordinate
% of sources to along a custom direction. Basically it is a coordinate 
% frame transformation, ie., get the location in a custom frame. Call
% the 'coordinate_rot' for actual transform.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/01/26
% Last modified date:   2023/01/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nlocxy = olocxy;
for i = 1: size(olocxy,1)
  [nlocxy(i,1), nlocxy(i,2)] = coordinate_rot(olocxy(i,1),olocxy(i,2),-(projang-90),[0 0]);
end
projx = nlocxy(:,1);
orty = nlocxy(:,2);