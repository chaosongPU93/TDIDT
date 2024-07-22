% synthshift_onespot_multibsts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the script to call 'synthshift_onespot_fn.m' to generate synthetics 
% from a single spot, different noise levels and saturation levels. This script
% sets up some parameters for the synthetics.
% --If you want to generate a single length of synthetics of the same kind,
% run 'synthshift_onespot.m'.
% --See also 'synthshift_chao_multibsts.m' for synthetics from diff region
% sizes and saturation levels with no noise.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/05/05
% Last modified date:   2024/05/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
format short e   % Set the format to 5-digit floating point
clear
clc
close all

% set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');

%whether to normalize templates
normflag = 0;

%which templates to use
tempflag = 'chao';

%specify the amplitude-frequency (counts) distribution as uniform
distr='UN';  % uniform distribution

sps = 160;      
greenlen = pow2(9)*sps/40;

bufsec = 1;
msftaddm = bufsec*sps;  %buffer range for later CC alignment, +1 for safety
rccmwsec = 0.5;
rccmwlen = rccmwsec*sps;  %window length for computing RCC
overshoot = rccmwlen/2; %RCC points to the center of computing window, so extra room is needed

%length of each simulation
% Twin=0.5*3600+(greenlen+msftaddm*2+overshoot*2-2)/sps; %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
% nrun = 6;

Twin=3*3600+(greenlen+msftaddm*2+overshoot*2-2)/sps; %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
nrun = 1;

%%%for the duration of templates, there are several options
%%%1. (ppeak-npeak)*2 of the BB template: 44;54;38;38 0.275s, 0.3375s, 0.2375s, 0.2375s
%%%2. direct eyeballing for between zerocrossings: ~65, 0.4s
%%%3. binned peak-to-peak separation for decently saturated unfiltered synthetics: ~37, 0.25s
%%%4. similar to 3, but synthetics are filtered first: ~37, 0.25s
% tdura = 0.4;  % duration from Chao's broadband template, width is about 795-730=65 spls at 160 hz
tdura = 0.25; % start to use on 2023/12/07

%%%specify which time is uniform in time
% timetype = 'tarvl';
timetype = 'tori';

%%%specify regime for transformation from time offset to map location 
% ftrans = 'interpArmb';
% ftrans = 'interpArmbreloc';
ftrans = 'interpchao';

%%%whether to plot to check the synthetics
pltsynflag = 0;
% pltsynflag = 1;

%%%whether to plot to check the new synthetics with added noise
pltnewsynflag = 0;
% % pltnewsynflag = 1;

%%%flag for validating if ground truth of sources can recover the record
testsrcflag = 0;
% testsrcflag = 1;

%different percent of noise
perctrial = 0.1*(0:2:16)';
ntrial = length(perctrial);

%variation of saturation level
nsat=[0.1 0.4 1 2 4 10 20 40 100];  % times of saturation, # of arrivals per duration 

parfor irun = 1: nrun
%   for iperc = 1: ntrial
%     perc = perctrial(iperc);  %ratio of noise to add into synthetics
%     
    %seed for random number generator
    % seed=3; 
    seed=irun*ntrial; 
%     seed=(irun-1)*ntrial+iperc; 

    synthshift_onespot_fn(perctrial,nsat,distr,Twin,tdura,seed,timetype,...
      normflag,tempflag,ftrans,pltsynflag,pltnewsynflag,testsrcflag,irun);
%     keyboard
%   end
end








