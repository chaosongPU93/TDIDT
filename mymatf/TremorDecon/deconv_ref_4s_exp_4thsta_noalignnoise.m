% deconv_ref_4s_exp_4thsta_noalignnoise.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the 'driver' script to call 'deconv_ref_4s_exp_4thsta_fn' to 
% re-run the deconvolution to original synthetic noise but DO NOT align the
% both the whole-win and 25-s wins
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/11/05
% Last modified date:   2024/11/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
% format short e   % Set the format to 5-digit floating point
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
% trange = trange(1:end-1,:);
tlen = trange(:,3)-trange(:,2);
tlensum = sum(tlen);

modname = 'timeoff_plfit_4thsta_160sps.mat';
planefit = load(strcat(rstpath, '/MAPS/',modname));
rmse = planefit.gof{4}.rmse;
offmax = round(2.0*rmse);

sps = 160;

%% call function 'deconv_ref_4s_exp_4thsta_fn_diffseed', all bursts NOISE
%%%Flag to indicate if it is necessary to recalculate everything
% flagrecalc = 0;
flagrecalc = 1;

if flagrecalc
  normflag = 0; %do not normalize the templates
  pltflag = 0;  %do not create summary plots for each burst
%   rccmwsec = 0.25; %use 0.5s or 0.25s 
  rccmwsec = 0.5; %use 0.5s or 0.25s
  
  %%%all bursts using synthetic noise
  noiseflag = 1; %
  alignflag = 0; %do not align the trio stations fror 25-s windows
  tmp = setdiff(1:size(trange,1),45);
  allbstnoi = deconv_ref_4s_exp_4thsta_fn(tmp,normflag,noiseflag,pltflag,rccmwsec,alignflag);  %
  savefile = 'deconv_stats4th_no23_allbstnoi_noalign.mat';
%   savefile = 'deconv_stats4th_no23_allbstnoi0.25s_noalign.mat';
  save(strcat(rstpath, '/MAPS/',savefile), 'allbstnoi');
  
else
  savefile = 'deconv_stats4th_no23_allbstnoi_noalign.mat';
  load(strcat(rstpath, '/MAPS/',savefile));
end


% keyboard