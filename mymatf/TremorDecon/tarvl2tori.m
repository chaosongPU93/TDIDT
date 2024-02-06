function [torispl,impsort,indsort]=tarvl2tori(imp,sps,ftrans,sortflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [torispl,impsort,indsort]=tarvl2tori(imp,sps,ftrans,sortflag)
%
% This function is to compute the relative time origin time of each event  
% based on the arrival time, and the location in terms of time offset that
% is store in different columns in 'imp'. If the 'sortflag' is true or 1,
% then the sorted events or the order based on the origin time are also
% returned.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/11/04
% Last modified date:   2024/01/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('sps',160);
defval('ftrans','interpchao');
defval('sortflag',0);

[imploc0, ~] = off2space002([0 0],sps,ftrans,0);  % a ref source at 0,0
[imploc, ~] = off2space002(imp(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
tcor = round((imploc(:,6)-imploc0(6))*sps);   % travel time difference between each source and a ref source at 0,0
torispl = imp(:,1)-tcor;

if sortflag
  [torispl, indsort] = sortrows(torispl,1);
  impsort = imp(indsort, :);
else
  impsort = [];
  indsort = [];
end
