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

cutout = 'ellipse';
ttol = 35;
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.pgc002.',cutout(1:4)));

modname = 'timeoff_plfit_4thsta_160sps.mat';
planefit = load(strcat(rstpath, '/MAPS/',modname));
rmse = planefit.gof{4}.rmse;
offmax = round(2.0*rmse);

%% some options for burst windows depending on the assumed noise level 
%empirically determined indices of bursts fall into local day times (noisier)
%or night times (quieter)
inbst = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,...
  34,35,56,57,58,59,60,61,62,70,71,72,73,74,75,76,77,78,79,80,81,82,110,111,112,113,114,115,116,117,...
  118,119,120,142,143,144,145,146,147,148,149,150,151,152,153,154,172,173,174,175]';
idbst = [36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,63,64,65,66,67,68,69,83,84,85,...
  86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,121,122,123,124,...
  125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,155,156,157,158,159,160,161,...
  162,163,164,165,166,167,168,169,170,171,176,177,178,179,180,181,182,183,184,185,186,187,188,189,...
  190,191,192,193,194,195]';

%indices of bursts whose CC of sig-wlet cc between sta 1 and 4 are high, above 80th percentile
ind41 = [1,4,6,10,12,15,21,22,23,24,29,46,52,56,74,77,78,82,83,102,111,113,114,121,124,125,...
  129,145,151,152,155,166,167,169,174,175,180,183,187]';

%indices of bursts whose CC of sig-wlet cc between sta 14, 24, 34 are all high, above 80th percentile
ind4123 = [1,12,15,46,56,77,102,114,125,129,155,180]';

%indices of bursts whose CC of envelope between sta 1,2,3 are high, above 75th percentile
indenv123 = [1,3,8,10,31,67,78,81,82,91,102,108,114,116,121,129,146,153,167,175];

%indices of bursts whose CC of sig-wlet cc between sta 1,2,3 are high, above 75th percentile
indswcc123 = [1,2,3,6,8,18,21,46,56,77,78,82,83,84,102,107,114,125,127,145,169];

%% call function 'deconv_ref_4s_exp_rand_fn', high-correlation VS low-correlation bursts, DATA VS NOISE
flagrecalc = 0;
% flagrecalc = 1;

if flagrecalc
  normflag = 0; %do not normalize the templates
  pltflag = 0;  %do not create summary plots for each choice of inputs
  % rccmwsec = 0.25; %use 0.5s or 0.25s 
  rccmwsec = 0.5; %use 0.5s or 0.25s
  
  %%%high-correlation bursts using real data
  noiseflag = 0;
  hibstsig = deconv_ref_4s_exp_4thsta_fn(indenv123,normflag,noiseflag,pltflag,rccmwsec);
  savefile = 'deconv_stats4th_hienvccbstsig.mat';
  save(strcat(rstpath, '/MAPS/',savefile), 'hibstsig');
  
  %%%high-correlation bursts using synthetic noise
  noiseflag = 1;
  hibstnoi = deconv_ref_4s_exp_4thsta_fn(indenv123,normflag,noiseflag,pltflag,rccmwsec);
  savefile = 'deconv_stats4th_hienvccbstnoi.mat';
  save(strcat(rstpath, '/MAPS/',savefile), 'hibstnoi');

else
  savefile = 'deconv_stats4th_hienvccbstsig.mat';
  load(strcat(rstpath, '/MAPS/',savefile));
  savefile = 'deconv_stats4th_hienvccbstnoi.mat';
  load(strcat(rstpath, '/MAPS/',savefile));
end

bstsig = hibstsig;
bstnoi = hibstnoi;

% keyboard

%% call function 'deconv_ref_4s_exp_rand_fn', all bursts, DATA VS NOISE
%%%Flag to indicate if it is necessary to recalculate everything
flagrecalc = 0;
% flagrecalc = 1;

if flagrecalc
  normflag = 0; %do not normalize the templates
  pltflag = 0;  %do not create summary plots for each choice of inputs
  % rccmwsec = 0.25; %use 0.5s or 0.25s 
  rccmwsec = 0.5; %use 0.5s or 0.25s
  
  %%%all bursts using real data
  noiseflag = 0;
  allbstsig = deconv_ref_4s_exp_4thsta_fn(1:size(trange,1),normflag,noiseflag,pltflag,rccmwsec); %
  savefile = 'deconv_stats4th_allbstsig.mat';
%   save(strcat(rstpath, '/MAPS/',savefile), 'allbstsig');
  
  %%%all bursts using synthetic noise
  noiseflag = 1;
  allbstnoi = deconv_ref_4s_exp_4thsta_fn(1:size(trange,1),normflag,noiseflag,pltflag,rccmwsec);  %
  savefile = 'deconv_stats4th_allbstnoi.mat';
%   save(strcat(rstpath, '/MAPS/',savefile), 'allbstnoi');
else
  savefile = 'deconv_stats4th_allbstsig.mat';
  load(strcat(rstpath, '/MAPS/',savefile));
  savefile = 'deconv_stats4th_allbstnoi.mat';
  load(strcat(rstpath, '/MAPS/',savefile));
end

bstsig = allbstsig;
bstnoi = allbstnoi;
% keyboard

%% comparison plots, syn noise vs data
%%%the direct/scaled deconvolved pos/neg source peak ratio between all station pairs, for each
%%%burst win separately
nsrc = bstsig.nsrc;
msrcampr = bstsig.msrcampr;
madsrcampr = bstsig.madsrcampr;
f1 = initfig(16,5,1,size(msrcampr,2));
f1 = plt_deconpk_rat4th(f1,msrcampr,madsrcampr,nsrc,'b');

nsrc = bstnoi.nsrc;
msrcampr = bstnoi.msrcampr;
madsrcampr = bstnoi.madsrcampr;
f1 = plt_deconpk_rat4th(f1,msrcampr,madsrcampr,nsrc,'r');

%%
%%%combine the direct/scaled deconvolved pos/neg source peak ratio between all station pairs of
%%%all burst wins, and summarize into one histogram
srcamprall = bstsig.srcamprall;
impindepstall = bstsig.impindepstall;
f2 = initfig(16,9,2,size(srcamprall,2));
f2 = plt_deconpk_rat_comb4th(f2,srcamprall,impindepstall,'b');

srcamprall = bstnoi.srcamprall;
impindepstall = bstnoi.impindepstall;
f2 = plt_deconpk_rat_comb4th(f2,srcamprall,impindepstall,'r');

%%
%%%deviation of source amp ratio from some median vs. RCC
lgdevsrcamprall = bstsig.lgdevsrcamprall;
rccpairsrcall = bstsig.rccpairsrcall;
rcccatsrcall = bstsig.rcccatsrcall;
f3 = initfig(15,7,2,size(lgdevsrcamprall,2)); %initialize fig
f3 = plt_deconpk_ratdevvsrcc4th(f3,lgdevsrcamprall,rccpairsrcall,rcccatsrcall,'b');

lgdevsrcamprall = bstnoi.lgdevsrcamprall;
rccpairsrcall = bstnoi.rccpairsrcall;
rcccatsrcall = bstnoi.rcccatsrcall;
f3 = plt_deconpk_ratdevvsrcc4th(f3,lgdevsrcamprall,rccpairsrcall,rcccatsrcall,'r');

%%
%%%diff between predicted arrival and selected peak at 4th sta
pred4offtrall = bstsig.pred4offtrall;
f4 = initfig(5,5,1,size(pred4offtrall,2)); %initialize fig
f4 = plt_errorof4thtarvlpred(f4,pred4offtrall,offmax,'b');

pred4offtrall = bstnoi.pred4offtrall;
f4 = plt_errorof4thtarvlpred(f4,pred4offtrall,offmax,'r');

%%
%%%preserved sources' amp ratio between 4th and 1st stas
srcamprall = bstsig.srcamprall;
impindepstall = bstsig.impindepstall;
f5 = initfig(12,5,1,3); %initialize fig
f5 = plt_deconpk_rat14(f5,impindepstall,srcamprall,'b');

srcamprall = bstnoi.srcamprall;
impindepstall = bstnoi.impindepstall;
f5 = plt_deconpk_rat14(f5,impindepstall,srcamprall,'r');

%%
%%%cloest waveform peak ratio between sta pairs VS. decon src amp ratio between same pairs
clppkhtwfall = bstsig.clppkhtwfall;
psrcampsall = bstsig.psrcampsall;
clnpkhtwfall = bstsig.clnpkhtwfall;
nsrcampsall = bstsig.nsrcampsall;
f6 = initfig(15,7,2,size(clppkhtwfall,2)); %initialize fig
f6 = plt_deconpkratvswfpkrat4th(f6,clppkhtwfall,psrcampsall,clnpkhtwfall,nsrcampsall);
% f6 = initfig(15,7,2,4); %initialize fig
% clppkhtwfall = [clppkhtwfall(:,1:3) clppkhtwfall(:,end)];  %I only want KLNB
% clnpkhtwfall = [clnpkhtwfall(:,1:3) clnpkhtwfall(:,end)];  %I only want KLNB
% f6 = plt_deconpkratvswfpkrat4th(f6,clppkhtwfall,psrcampsall,clnpkhtwfall,nsrcampsall);

clppkhtwfall = bstnoi.clppkhtwfall;
psrcampsall = bstnoi.psrcampsall;
clnpkhtwfall = bstnoi.clnpkhtwfall;
nsrcampsall = bstnoi.nsrcampsall;
f7 = initfig(15,7,2,size(clppkhtwfall,2)); %initialize fig
f7 = plt_deconpkratvswfpkrat4th(f7,clppkhtwfall,psrcampsall,clnpkhtwfall,nsrcampsall);

%%
%%%cloest waveform peak ratio between sta pairs VS. decon src amp ratio between same pairs
clppkhtwfall = bstsig.clppkhtwfall;
clnpkhtwfall = bstsig.clnpkhtwfall;
clpkspanwfall = clppkhtwfall-clnpkhtwfall;
psrcampsall = bstsig.psrcampsall;
nsrcampsall = bstsig.nsrcampsall;
srcampsspanall = psrcampsall-nsrcampsall;
f7 = initfig(15,4,1,size(clppkhtwfall,2)); %initialize fig
f7 = plt_deconpkratvswfpkrat4th(f7,clpkspanwfall,srcampsspanall);


%% how many sources from noise compared to data, after checked at 4th stations
nsrcs = bstsig.nsrc;  %number of sources from signal
nsrcn = bstnoi.nsrc;  %number of sources from synthetic noise

frac = nsrcn./nsrcs; %fraction during day (noisy) times

f8 = initfig(5,4,1,1); %initialize fig
ax = f8.ax(1);
hold(ax,'on');
grid(ax,'on');
histogram(ax,frac,'BinWidth',0.05);
plot(ax,[median(frac) median(frac)],ax.YLim,'--','color','r','linew',1);
ylabel(ax,'Count of burst wins');
xlabel(ax,'Fraction');
% title(ax,'Burst wins with high-correlation between 14');
title(ax,'All burst wins');

keyboard


