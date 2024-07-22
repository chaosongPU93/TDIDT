% synthshift_regandnoi_multibsts.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the script to call 'synthshift_regandnoi_fn.m' to generate synthetics
% with diff. saturation levels made by sources from the combination of a small 
% set of region sizes and with a small set of noise levels. 
% For generating synthetics EITHER from various region sizes with no added 
% noise, refer to 'synthshift_chao_fn' and 'synthshift_chao_multibsts',
% OR from a single spot with various noise levels, refer to
% 'synthshift_onespot_fn' and 'synthshift_onespot_multibsts'.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/07/20
% Last modified date:   2024/07/20
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

%%%specify distribution for source location  
% distrloc = 'custompdf'; %using a custom PDF function
distrloc = 'uniform'; %uniformly random in a specified region,

%%%specify if considering the physical size of each source
% physicalsize = 1;
physicalsize = 0;

%%%specify shape of the source region
srcregion='ellipse';
% srcregion='rectangle';
% srcregion='circle';

%%%specify if forcing a min speration of arrival time for 2 events from the same spot 
forcesep = 0;
% forcesep = 1;

%whether to show a plot of raw synthetics
pltflag = 0;
% pltflag = 1;

%whether to test if synthetics can be reproduced by convolution
testsrcflag = 0;
% testsrcflag = 1;

%variation of source region size
if strcmp(srcregion,'ellipse')
  semia = 1.75*(0.8:0.2:1.2);
  semib = 1.25*(0.8:0.2:1.2);
  nreg = length(semia);
%   semia = 3;
%   semib = 1.5;
%   nreg = length(semia);
end

%variation of saturation level
% nsat=[0.1 0.4 1 2 4 10 20 40 100];  % times of saturation, # of arrivals per duration 
% nsat=[0.4 1 2 4 10 20 40 100];  % times of saturation, # of arrivals per duration 
nsat=[0.4 1 2 4 10 20];  % times of saturation, # of arrivals per duration 

%%%whether to plot to check the synthetics
pltsynflag = 0;
% pltsynflag = 1;

%%%whether to plot to check the new synthetics with added noise
pltnewsynflag = 0;
% % pltnewsynflag = 1;

%different percent of noise
perctrial = 0.1*(4:2:8)';
ntrial = length(perctrial);

for irun = 1: nrun
%   parfor ireg = 1: nreg
  for ireg = 1: nreg
%     ireg = 3;
    xaxis = semia(ireg); %axis length of the same ellipse of my 4-s catalog
    yaxis = semib(ireg);

    %seed for random number generator
    % seed=3; 
%     seed=ireg; 
    seed=(irun-1)*nreg+ireg; 

    synthshift_regandnoi_fn(srcregion,xaxis,yaxis,perctrial,nsat,distr,Twin,tdura,...
      seed,distrloc,timetype,physicalsize,forcesep,testsrcflag,normflag,...
      tempflag,ftrans,pltsynflag,pltnewsynflag,pltflag,irun);

  end
end






