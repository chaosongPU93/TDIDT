%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Read the one day of Polaris station data, rotate according to split
% correction angle, polarization angle, to the optimal & orthogonal
% direction.
% 
% Modified from readpols.m by Allan Rubin.
%
%
% By Chao Song, chaosong@princeton.edu
% Last modified date:   2019/03/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [opt, ort, timepola] = readpola_1daynofilter(datafnm, POLSTA, ...
    POLROTS, idx, sps, fact)

STAEdat = [datafnm,'.',POLSTA(idx,1:4),'..HHE.D.SAC']; %HHE for polstas.
STANdat = [datafnm,'.',POLSTA(idx,1:4),'..HHN.D.SAC'];
[STAE, ~, ~, ~, timepola] = readsac(STAEdat, 0, 'l');   %%% Function 'readsac'
[STAN, ~, ~, ~, ~] = readsac(STANdat,0, 'l');    %%% 'l' means linux
tracelen=length(STAE);   % trace length of N/E comp

%%% TAPER, cosine taper over 1 second (at 40 sps) before filtering:
x = (0: pi/200: pi/2- pi/200)';
STAE(1: 100) = sin(x).* STAE(1: 100); %Only at start of day on E component!
STAN(1: 100) = sin(x).* STAN(1: 100);
x = flipud(x);    % OR, could use cos(x) instead of sin
STAE(tracelen- 99: tracelen) = sin(x).* STAE(tracelen- 99: tracelen);
STAN(tracelen- 99: tracelen) = sin(x).* STAN(tracelen- 99: tracelen);

%%% SCALE and detrend
STAE = fact* detrend(STAE(:));  
STAN = fact* detrend(STAN(:));  

% resample if needed, by Chao
[num, denom] = rat(sps/100);
STAEd = resample(STAE,num,denom);   % times = num/denom
STANd = resample(STAN,num,denom);
timepola = linspace(round(timepola(1)), round(timepola(end)),...
        round(length(timepola)*sps/40));

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
