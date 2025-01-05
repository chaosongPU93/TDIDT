% deconv_4s_exp_4thsta_suppledemo.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the script to call the function "deconv_4s_exp_4thsta_fn.m"
% to run the whole-win deconvolution to a few particular burst windows 
% that are shown in the supplement of the paper.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/11/23
% Last modified date:   2024/11/23
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
ntol = 3;
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.',...
  num2str(ntol),'.pgc002.',cutout(1:4)));
tlen = trange(:,3)-trange(:,2);

modname = 'timeoff_plfit_4thsta_160sps.mat';
planefit = load(strcat(rstpath, '/MAPS/',modname));
rmse = planefit.gof{4}.rmse;
offmax = round(2.0*rmse);

sps = 160;

indhi = [99; 100; 115; 164; 170; 174; 181]; 

%% call function 'deconv_ref_4s_exp_rand_fn', all bursts, DATA VS NOISE
%%%Flag to indicate if it is necessary to recalculate everything
% flagrecalc = 0;
flagrecalc = 1;

if flagrecalc
  normflag = 0; %do not normalize the templates
  pltflag = 0;  %do not create summary plots for each choice of inputs
  % rccmwsec = 0.25; %use 0.5s or 0.25s 
  rccmwsec = 0.5; %use 0.5s or 0.25s
  
  %%%all bursts using real data
  noiseflag = 0;
  allbstsig = deconv_4s_exp_4thsta_fn(indhi,normflag,noiseflag,pltflag,rccmwsec); %

  savefile = 'deconv1win_stats4th_no23_allbstsig_indhi.mat';
  save(strcat(rstpath, '/MAPS/',savefile), 'allbstsig');
  
end