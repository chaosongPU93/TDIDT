function [STAEd,STANd,timepola] = readpols_norotsv2(datafnm,POLSTA,idx,sps,fact,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Read the one day of Polaris station data, rotate according to split
% correction angle, polarization angle, to the optimal & orthogonal
% direction. version 2 for data removed from station response
% 
% USE data removed from station response 
%
% Modified from readpols.m by Allan Rubin.
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/06/29
% Last modified date:   2023/07/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% decide which component to read based on the number of inputs, if ==14, it means that 'COMPNAME' 
%%% is not defined, so do the default, ie., read the horizontal components
if ~isempty(varargin)
  COMPNAME = (varargin{1});
else
  COMPNAME = [];
end

if isempty(COMPNAME)
  STAEdat = [datafnm,'.',POLSTA(idx,1:4),'..HHE.D.SAC']; %HHE for polstas.
  STANdat = [datafnm,'.',POLSTA(idx,1:4),'..HHN.D.SAC'];
  [STAE, ~, ~, ~, timepola] = readsac(STAEdat, 0, 'l');   %%% Function 'readsac'
  [STAN, ~, ~, ~, ~] = readsac(STANdat,0, 'l');    %%% 'l' means linux

  %%% SCALE
  STAE = fact* detrend(STAE(:)); 
  STAN = fact* detrend(STAN(:));

  % resample if needed, by Chao
  [num, denom] = rat(sps/100);
  STAEd = resample(STAE,num,denom);   % times = num/denom
  STANd = resample(STAN,num,denom);
  timepola = linspace(round(timepola(1)), round(timepola(end)),...
          round(length(timepola)*sps/40));

%%% this means the 'COMPNAME' is defined as 'Z'
elseif strcmp(COMPNAME, 'Z')
  STAZdat = [datafnm,'.',POLSTA(idx,1:4),'..HHZ.D.SAC']; %HHE for polstas.
  [STAZ, ~, ~, ~, timepola] = readsac(STAZdat, 0, 'l');   %%% Function 'readsac'

  %%% SCALE
  STAZ = fact* detrend(STAZ(:)); 
  
  % resample if needed, by Chao
  [num, denom] = rat(sps/100);
  STAZd = resample(STAZ,num,denom);   % times = num/denom
  timepola = linspace(round(timepola(1)), round(timepola(end)),...
          round(length(timepola)*sps/40));
  STAEd=STAZd;
  STANd=0*STAEd;
end

