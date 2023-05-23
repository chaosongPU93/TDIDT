% lfecaterror.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We already know that through many different ways to change the order of 
% deconvolution, the resulting catalogs can be a little different, in terms
% of number of sources, source location, etc. We might not care that much 
% about the exact number, since any catalog might overestimate it. But we 
% want to know the location difference, which may enlight us on the error
% of the catalog, in particular the arrival time error. The way to do it
% is to evaluate the location difference in samples for the sources that 
% are associated with the same waveform peak.
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/02/10
% Last modified date:   2023/02/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
format short e   % Set the format to 5-digit floating point
clear
clc
close all

% set(0,'DefaultFigureVisible','on');
set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, resol] = pixelperinch(1);

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');

cutout = 'ellipse';
ttol = 35;
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.pgc002.',cutout(1:4)));

modname = 'timeoff_plfit_4thsta_160sps.mat';
planefit = load(strcat(rstpath, '/MAPS/',modname));
rmse = planefit.gof{4}.rmse;
offmax = round(2.0*rmse);

%% call function to carry out deconvolution in different order
normflag = 0; %do not normalize the templates
pltflag = 0;  %do not create summary plots for each choice of inputs
% rccmwsec = 0.25; %use 0.5s or 0.25s
rccmwsec = 0.5; %use 0.5s or 0.25s
noiseflag = 0;

%%%use the current scheme, start from peak with the largest weighted height by RCC
allbstsig = deconv_ref_4s_exp_4thsta_fn(181,normflag,noiseflag,pltflag,rccmwsec); %

%%%use the alternative scheme, start from chronologically the earliest peak above the threshold
allbstnoi = deconv_ref_4s_exp_4thsta_fn(181,normflag,noiseflag,pltflag,rccmwsec);  %

%% associate sources with wavefrom peaks, identify those linked with same peaks

%% what is location difference of them in samples as a vector?

%% move to a common origin, what's the shape like formed by all vectors?






















