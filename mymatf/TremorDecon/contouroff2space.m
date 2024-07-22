function conmatxy = contouroff2space(conmat,sps,ftrans)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conmatxy = contouroff2space(conmat,sps,ftrans)
% 
% This function converts the contour lines in offsets in time (s) to map view 
% (km), by mapping all data points on the contour lines expressed in arrival 
% time offset in seconds to map locations in km using the same mapping 
% function 'off2space002'. Linear interpolation is adopted.
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/03/23
% Last modified date:   2024/03/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('sps',160);
defval('ftrans','interpchao');

locall = off2space002([],sps,ftrans); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
Fx = scatteredInterpolant(locall(:,7),locall(:,8),locall(:,1),'linear','none');
Fy = scatteredInterpolant(locall(:,7),locall(:,8),locall(:,2),'linear','none');

grp = unique(conmat(:,2));
ngrp = max(grp);
conmatxy = conmat;
for igrp = 1:ngrp
  ind = find(conmat(:,2)==igrp);  
  conmatxy(ind,3) = Fx(conmat(ind,3)*sps,conmat(ind,4)*sps);
  conmatxy(ind,4) = Fy(conmat(ind,3)*sps,conmat(ind,4)*sps);
end
