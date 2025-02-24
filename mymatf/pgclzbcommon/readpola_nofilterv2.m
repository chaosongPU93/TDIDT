function [opt, ort, timepola] = readpola_nofilterv2(datafnm, POLSTA, ...
    POLROTS, idx, sps, fact, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Read the one day of Polaris station data, rotate according to split
% correction angle, polarization angle, to the optimal & orthogonal
% direction. version 2 for data removed from station response
% 
% Modified from readpols.m by Allan Rubin.
%
%
% By Chao Song, chaosong@princeton.edu
% Last modified date:   2019/06/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% decide which component to read. If not defined, it means that 'COMPNAME' 
%%% is empty, so do the default, ie., read the horizontal components
% if(POLSTA(idx,5:5)==' ')   %%% The name has 4 charcters
if isempty(varargin)
  STAEdat = [datafnm,'.',POLSTA(idx,1:4),'..HHE.D.SAC']; %HHE for polstas.
  STANdat = [datafnm,'.',POLSTA(idx,1:4),'..HHN.D.SAC'];
  [STAE, ~, ~, ~, timepola] = readsac(STAEdat, 0, 'l');   %%% Function 'readsac'
  [STAN, ~, ~, ~, ~] = readsac(STANdat,0, 'l');    %%% 'l' means linux

  %%% SCALE
  STAE = fact* detrend(STAE(:)); 
  STAN = fact*STAN(:);

  % resample if needed, by Chao
  [num, denom] = rat(sps/100);
  STAEd = resample(STAE,num,denom);   % times = num/denom
  STANd = resample(STAN,num,denom);
  timepola = linspace(round(timepola(1)), round(timepola(end)),...
          round(length(timepola)*sps/100));

  tracelen = length(STAEd); %after decimation/upsampling

  % convert the offset in 40 sps unit to the sps specified.
  %%% 1st column of PERMROTS == offset in samples between fast and slow components
  fastslowoff = round(POLROTS(idx,1)*sps/40);  % this offset is based on 40sps data, even for polaris stations
  %%% 4th column of PERMROTS == offset in samples between the arrival times at different stations
  STAsoff = round(POLROTS(idx,4)*sps/40);  %input offsets given at 40 sps.

  %Rotate & split-correct
  STA = STAEd+ 1i* STANd;
  STAfastslow = STA* exp(-1i* POLROTS(idx,2));   %%% 2th column of PERMROTS, == rotation angle to get fast/slow direction
  STAslow = real(STAfastslow);
  STAfast = imag(STAfastslow);
  len = length(STA);
  off = round(10*sps/40);
  if off < fastslowoff 
      off = round(15*sps/40);
  end
  %%% 1st column of PERMROTS == offset in samples between fast/slow direction
  STAslow(off: len-off) = STAslow(off+ fastslowoff: len- off+ fastslowoff); 
  STAsplitcorrected = (STAslow+ 1i* STAfast)* exp(1i* POLROTS(idx,2));
  STAscrot = STAsplitcorrected* exp(-1i* POLROTS(idx,3));

  %Timeshift
  if STAsoff > -1
      STAscrot(1: tracelen- STAsoff) = STAscrot(STAsoff+ 1: tracelen);
      STAscrot(tracelen- STAsoff+ 1: tracelen) = 0;
  else
      STAscrot(-STAsoff+ 1: tracelen) = STAscrot(1: tracelen+ STAsoff);
      STAscrot(1: -STAsoff) = 0;
  end

  opt = real(STAscrot);
  ort = imag(STAscrot);
  
% else    %%% The name has 5 characters??? ONLY read Z component
%%% this means the 'COMPNAME' is defined as 'Z'
else
  COMPNAME = (varargin{1});
  if strcmp(COMPNAME, 'Z')
    
    STAZdat=[datafnm,'.',POLSTA(idx,1:4),'..HHZ.D.SAC']; %BHE for permstas.
    [STAZ,~,~,~,timepola]=readsac(STAZdat,0,'l');
    
    STAZ = fact* detrend(STAZ(:));
    
    % resample if needed, by Chao
    [num, denom] = rat(sps/100);
    STAZd=resample(STAZ,num,denom);   % times = num/denom
    tracelen=length(STAZd); %after decimation
    timepola = linspace(round(timepola(1)), round(timepola(end)),...
      round(length(timepola)*sps/100));
    
    %Timeshift
    STAsoff=round(POLROTS(idx,4)*sps/40);  %input offsets given at 40 sps
    if STAsoff > -1
      STAZd(1:tracelen-STAsoff)=STAZd(STAsoff+1:tracelen);
      STAZd(tracelen-STAsoff+1:tracelen)=0;
    else
      STAZd(-STAsoff+1:tracelen)=STAZd(1:tracelen+STAsoff);
      STAZd(1:-STAsoff)=0;
    end
    opt=STAZd;
    %   ort = [];
    ort=0*opt;
  end

end




