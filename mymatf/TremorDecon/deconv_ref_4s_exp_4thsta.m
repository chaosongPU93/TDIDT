% deconv_ref_4s_exp_4thsta.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the 'driver' script to call 'deconv_ref_4s_exp_4thsta_fn' to do the
% analysis based on the desired 'indofburst'. 
% so that outputs from either from quieter burst windows (presumbly night
% times) or from noisier windows presumbly day times), could be directly
% compared in one single plot or more. Here i am more interested into how
% the signal waveform coherency between 4th stations and trio stations changes
% wrt to the background noise level (night or day). Maybe the prediction at
% a 4th station would be better if the burst was during night times given a 
% set of deconvolved sources using trio stations.
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/10/13
% Last modified date:   2022/10/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
format short e   % Set the format to 5-digit floating point
clear
clc
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, resol] = pixelperinch(1);

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');


%% some options for burst windows depending on the assumed noise level 
%empirically determined indices of bursts fall into local day times (noisier)
%or night times (quieter)
inbst = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,...
  34,35,56,57,58,59,60,61,62,70,71,72,73,74,75,76,77,78,79,80,81,82,110,111,112,113,114,115,116,117,...
  118,119,120,142,143,144,145,146,147,148,149,150,151,152,153,154,172,173,174,175];
idbst = [36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,63,64,65,66,67,68,69,83,84,85,...
  86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,121,122,123,124,...
  125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,155,156,157,158,159,160,161,...
  162,163,164,165,166,167,168,169,170,171,176,177,178,179,180,181,182,183,184,185,186,187,188,189,...
  190,191,192,193,194,195];

%indices of bursts whose wlet-sig cc between the opt and ort components are very similar, above 75
%percentile
ihighoo = [20,23,59,60,80,113,116,120,134,189,194];

%indices of bursts whose sig cc between sta 1/2/3 are very high, above 90 percentile
ihicc123 = [1,3,6,7,8,24,56,71,75,77,81,83,93,102,114,116,132,145,149,185];
ihicc123n = intersect(inbst,ihicc123);

indtest = [18,21,22,23];

%% call function 'deconv_ref_4s_exp_rand_fn'
normflag = 0; %do not normalize the templates
pltflag = 0;  %do not create summary plots for each choice of inputs

%%% quiet burst wins using real data
nbst = deconv_ref_4s_exp_4thsta_fn(inbst,normflag,pltflag);

dbst = deconv_ref_4s_exp_4thsta_fn(idbst,normflag,pltflag);








