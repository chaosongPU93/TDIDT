function [tarvl4,off14]=pred_tarvl_at4thsta(stanm,off12,off13,tarvl1,off1i,sps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tries to ease the predication of arrival times at the 4th sta
% of the same sources based on the source locations and arrival time at 
% the 1st station (or the aligned 2nd or 3rd stas). 
% Suppose that ALL seismograms have be aligned to image a source located at 
% the specific LFE family (eg, 002), the plane fit model give the empirical 
% off14 for each pair of off12 and off13 (the total offset). 'off1i' is the 
% EXTRA offset if for some window, station 4 is additionally best-aligned 
% relative to station 1, then the resulting off14 needs to corrected, so is
% the arrival time at sta 4. 
% The idea is, we get the
% empirical time off14 caused by the deviation of source location (off12,
% off13) from (0,0), then 'tarvl4' is simply 'tarvl1' subtracted by
% 'off14'. 'off1i' is the prealignment between 4th sta and 1st sta to the 
% entire window that needs to corrected as well.
%
% INPUT:  
%     off12: total offset PGC-SSIB, after alignment to 0,0, ie, loc of lfe fam
%     off13: total offset PGC-SILB, after alignment to 0,0, ie, loc of lfe fam
%     tarvl1: arrival time at PGC
%     off1i: ADDITIONAL offset to alignment to 0,0, ie, loc of lfe fam
%     sps: sampling rate of off12, off13, tarvl1, off1i
%
% OUTPUT:
%     off14: APParent offset PGC-4th sta after alignment to 0,0, ie, loc of lfe fam
%     tarvl4: APParent arrival at 4th sta after alignment to 0,0, ie, loc of lfe fam
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/10/04
% Last modified date:   2024/07/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('off12',0);
defval('off13',0);
defval('tarvl1',[]);
defval('off1i',0);
defval('sps',160);

%load the plane fitting model that is based on 160 hz data
workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/PGCtrio');
% modname = 'timeoff_planefitparam_4thsta_160sps';
% off14mod = load(strcat(rstpath, '/MAPS/',modname));
modname = 'timeoff_plfit_4thsta_160sps.mat';  %based on 160 hz data.
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

off14 = round(off14mod(idx,1).*(off12*160/sps) + ...
              off14mod(idx,2).*(off13*160/sps) + ...
              off14mod(idx,3));   %in 160 hz
off14 = off14/(160/sps);  %convert back to desired sps            
% off14 = round(off14mod(idx,:)*[off12; off13; 1]);

%%%2024/05/01
%plane fit is based on 40 sps data and has error, techinically, a source at 0,0
%should have been aligned at all stations, such that off14 for at 0,0 is zero
off14 = off14-round(off14mod(idx,3)/(160/sps));

%calibrate for the waveform prealignment, if any, double-checked on 2024/07/19
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

