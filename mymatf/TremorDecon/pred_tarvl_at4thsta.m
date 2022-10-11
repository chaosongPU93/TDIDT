function tarvl4=pred_tarvl_at4thsta(stanm,off12,off13,tarvl1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tries to ease the predication of arrival times at the 4th sta
% of the same sources based on the source locations and arrival time at 
% the 1st station (or the aligned 2nd or 3rd stas). The idea is, we get the
% empirical time off14 for source location (off12,off13), then 'tarvl4' is
% simply the sum of 'tarvl1' and this 'off14'.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/10/04
% Last modified date:   2022/10/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the plane fitting model 
off14mod = load(strcat(rstpath, '/MAPS/timeoff_planefitparam_4thsta'),'w+');

%candidate new 4th stations
stasnew=['LZB  '
         'TWKB '
         'MGCB '
         'KLNB '];  % twkb lzb mgcb
nstasnew = size(stasnew,1); 

%find which 4th sta is being requested
[~,idx]=ismember(stanm,stasnew,'rows');

%the simple plane fitting model is:  z = a.*x + b.*y + c 
off14 = off14mod(idx,1).*off12 + off14mod(idx,2).*off13 + off14mod(idx,3);

%the arrival time at 4th sta is just arrival time at 1st sta plus the off14
tarvl4 = tarvl1 + off14;



