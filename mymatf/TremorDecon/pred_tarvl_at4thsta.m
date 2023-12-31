function [tarvl4,off14]=pred_tarvl_at4thsta(stanm,off12,off13,tarvl1,off1i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tries to ease the predication of arrival times at the 4th sta
% of the same sources based on the source locations and arrival time at 
% the 1st station (or the aligned 2nd or 3rd stas). The idea is, we get the
% empirical time off14 caused by the deviation of source location (off12,
% off13) from (0,0), then 'tarvl4' is simply 'tarvl1' subtracted by
% 'off14'. 'off1i' is the prealignment between 4th sta and 1st sta to the 
% entire window that needs to corrected as well.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/10/04
% Last modified date:   2022/10/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('tarvl1',[]);
defval('off1i',0);

%load the plane fitting model 
workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/PGCtrio');
% modname = 'timeoff_planefitparam_4thsta_160sps';
% off14mod = load(strcat(rstpath, '/MAPS/',modname));
modname = 'timeoff_plfit_4thsta_160sps.mat';
planefit = load(strcat(rstpath, '/MAPS/',modname));
off14mod = planefit.modparam;
% rmse = modelfit.gof.rmse;

%candidate new 4th stations
stasnew=['LZB  '
         'TWKB '
         'MGCB '
         'KLNB '];  % twkb lzb mgcb

%find which 4th sta is being requested
[~,idx]=ismember(stanm,stasnew,'rows');

% keyboard
%the simple plane fitting model is:  z = a.*x + b.*y + c 
off14 = round(off14mod(idx,1).*off12 + off14mod(idx,2).*off13 + off14mod(idx,3));
% off14 = round(off14mod(idx,:)*[off12; off13; 1]);

%calibrate for the waveform prealignment, if any
off14 = off14 - off1i*ones(size(off14)); 

if ~isempty(tarvl1)
  %the arrival time at 4th sta is just arrival time at 1st sta minus the off14
  %because, off14 = tarvl1-tarvl4, so if off14>0, sta4 needs to be shifted to the right to align with
  %sta1
  %This also implies that, if 'off1i' is not zero, say >0, ie., the 4th sta trace has been prealigned
  %with 1st sta, 'tarvl4' would be larger, ie., tarvl4 = tarvl4 + off1i
  tarvl4 = tarvl1 - off14;
else
  tarvl4 = [];
end

% keyboard

