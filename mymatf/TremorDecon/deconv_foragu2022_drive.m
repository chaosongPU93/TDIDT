% deconv_foragu2022_drive.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the 'driver' script to call 'deconv_foragu2022' to do the
% analysis based on the desired 'indofburst'. 'deconv_foragu2022' is a
% special version of 'deconv_ref_4s_exp_4thsta_fn' for AGU 2022, contains lots
% of minor modifications of plotting. We need this driver for more statistics
% especially the distance vs. time for the sources from both data and synthetic
% noise.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/12/08
% Last modified date:   2022/12/08
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

%% call function 'deconv_ref_4s_exp_rand_fn', high-correlation bursts, DATA VS NOISE
% flagrecalc = 0;
% % flagrecalc = 1;
% 
% if flagrecalc
%   normflag = 0; %do not normalize the templates
%   pltflag = 0;  %do not create summary plots for each choice of inputs
%  % rccmwlen = sps/4; %use 0.5s or 0.25s 
%  % rccmwlen = sps/2; %use 0.5s or 0.25s
%   
%   %%%high-correlation bursts using real data
%   noiseflag = 0;
%   hibstsig = deconv_foragu2022(indenv123,normflag,noiseflag,pltflag,rccmwlen);
%   savefile = 'deconv_agu2022stats4th_hienvccbstsig.mat';
%   save(strcat(rstpath, '/MAPS/',savefile), 'hibstsig');
%   
%   %%%high-correlation bursts using synthetic noise
%   noiseflag = 1;
%   hibstnoi = deconv_foragu2022(indenv123,normflag,noiseflag,pltflag,rccmwlen);
%   savefile = 'deconv_agu2022stats4th_hienvccbstnoi.mat';
%   save(strcat(rstpath, '/MAPS/',savefile), 'hibstnoi');
% 
% else
%   savefile = 'deconv_agu2022stats4th_hienvccbstsig.mat';
%   load(strcat(rstpath, '/MAPS/',savefile));
%   savefile = 'deconv_agu2022stats4th_hienvccbstnoi.mat';
%   load(strcat(rstpath, '/MAPS/',savefile));
% end
% 
% bstsig = hibstsig;
% bstnoi = hibstnoi;
% sps = 160;
% 
% % keyboard

%% call function 'deconv_ref_4s_exp_rand_fn', all bursts, DATA VS NOISE
%%%Flag to indicate if it is necessary to recalculate everything
flagrecalc = 0;
% flagrecalc = 1;

if flagrecalc
  normflag = 0; %do not normalize the templates
  pltflag = 0;  %do not create summary plots for each choice of inputs
  % rccmwlen = sps/4; %use 0.5s or 0.25s 
  rccmwlen = sps/2; %use 0.5s or 0.25s

  %%%all bursts using real data
  noiseflag = 0;
  allbstsig = deconv_foragu2022(1:size(trange,1),normflag,noiseflag,pltflag,rccmwlen); %
  savefile = 'deconv_agu2022stats4th_allbstsig.mat';
%   save(strcat(rstpath, '/MAPS/',savefile), 'allbstsig');
  
  %%%all bursts using synthetic noise
  noiseflag = 1;
  allbstnoi = deconv_foragu2022(1:size(trange,1),normflag,noiseflag,pltflag,rccmwlen);  %
  savefile = 'deconv_agu2022stats4th_allbstnoi.mat';
%   save(strcat(rstpath, '/MAPS/',savefile), 'allbstnoi');
else
  savefile = 'deconv_agu2022stats4th_allbstsig.mat';
  load(strcat(rstpath, '/MAPS/',savefile));
  savefile = 'deconv_agu2022stats4th_allbstnoi.mat';
  load(strcat(rstpath, '/MAPS/',savefile));
end

bstsig = allbstsig;
bstnoi = allbstnoi;
sps = 160;
% keyboard

%% comparison plots, syn noise vs data
%%%the separation in arrival time vs. distance between sources N--N-1, N--N-2, and N--N-3,
%%%for data and noise
widin = 7;  % maximum width allowed is 8.5 inches
htin = 3.5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f1 = initfig(widin,htin,nrow,ncol); %initialize fig
xran = [0.1 0.94]; yran = [0.12 0.96];
xsep = 0.08; ysep = 0.05;
optaxpos(f1,nrow,ncol,xran,yran,xsep,ysep);

dtarvlnn1all = bstsig.dtarvlnn1all;
distarvlnn1all = bstsig.distarvlnn1all;
dtarvlnn2all = bstsig.dtarvlnn2all;
distarvlnn2all = bstsig.distarvlnn2all;
dtarvlnn3all = bstsig.dtarvlnn3all;
distarvlnn3all = bstsig.distarvlnn3all;
f1 = plt_absdistvsdt(f1,[dtarvlnn1all,distarvlnn1all],1,sps,'tarvl','b');
% median(distarvlnn1all(dtarvlnn1all/sps<=1))
% median(distarvlnn1all(dtarvlnn1all/sps<=0.5))
% median(distarvlnn1all)
% median(dtarvlnn1all/sps)
% legend(f1.ax(3),'Data','location','east');

dtarvlnn1all = bstnoi.dtarvlnn1all;
distarvlnn1all = bstnoi.distarvlnn1all;
dtarvlnn2all = bstnoi.dtarvlnn2all;
distarvlnn2all = bstnoi.distarvlnn2all;
dtarvlnn3all = bstnoi.dtarvlnn3all;
distarvlnn3all = bstnoi.distarvlnn3all;
f1 = plt_absdistvsdt(f1,[dtarvlnn1all,distarvlnn1all],1,sps,'tarvl','r');
ylim(f1.ax(1),[0 0.6]);
ylim(f1.ax(2),[0 3]);

% orient(f1.fig,'landscape');
% print(f1.fig,'-dpdf','/home/chaosong/Pictures/allbsttarvlnn1dist.pdf');
% print(f1.fig,'-dpdf','/home/chaosong/Pictures/hienvccbsttarvlnn1dist.pdf');

%%
%%%similar as above, but only for sources whose time separation from neighbors is within 1s
widin = 7;  % maximum width allowed is 8.5 inches
htin = 3.5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f1 = initfig(widin,htin,nrow,ncol); %initialize fig
xran = [0.1 0.94]; yran = [0.12 0.96];
xsep = 0.08; ysep = 0.05;
optaxpos(f1,nrow,ncol,xran,yran,xsep,ysep);

dtarvlnn1all = bstsig.dtarvlnn1all;
distarvlnn1all = bstsig.distarvlnn1all;
ind=find(dtarvlnn1all<=1.0*sps);
f1 = plt_absdistvsdt(f1,[dtarvlnn1all(ind),...
  distarvlnn1all(ind)],1,sps,'tarvl','b');

dtarvlnn1all = bstnoi.dtarvlnn1all;
distarvlnn1all = bstnoi.distarvlnn1all;
ind=find(dtarvlnn1all<=1.0*sps);
f1 = plt_absdistvsdt(f1,[dtarvlnn1all(ind),...
  distarvlnn1all(ind)],1,sps,'tarvl','r');
ylim(f1.ax(1),[0 0.6]);
ylim(f1.ax(2),[0 4.1]);


%%
%%%Same as above, but how about for sources after the 4th sta check
widin = 7;  % maximum width allowed is 8.5 inches
htin = 3.5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f2 = initfig(widin,htin,nrow,ncol); %initialize fig
xran = [0.1 0.94]; yran = [0.12 0.96];
xsep = 0.08; ysep = 0.05;
optaxpos(f2,nrow,ncol,xran,yran,xsep,ysep);

dtarvlnn14thall = bstsig.dtarvlnn14thall;
distarvlnn14thall = bstsig.distarvlnn14thall;
dtarvlnn24thall = bstsig.dtarvlnn24thall;
distarvlnn24thall = bstsig.distarvlnn24thall;
dtarvlnn34thall = bstsig.dtarvlnn34thall;
distarvlnn34thall = bstsig.distarvlnn34thall;
f2 = plt_absdistvsdt(f2,[dtarvlnn14thall,distarvlnn14thall],1,sps,'tarvl','b');

dtarvlnn14thall = bstnoi.dtarvlnn14thall;
distarvlnn14thall = bstnoi.distarvlnn14thall;
dtarvlnn24thall = bstnoi.dtarvlnn24thall;
distarvlnn24thall = bstnoi.distarvlnn24thall;
dtarvlnn34thall = bstnoi.dtarvlnn34thall;
distarvlnn34thall = bstnoi.distarvlnn34thall;
f2 = plt_absdistvsdt(f2,[dtarvlnn14thall,distarvlnn14thall],1,sps,'tarvl','r');
ylim(f2.ax(1),[0 0.6]);
ylim(f2.ax(2),[0 3]);
xlim(f2.ax(2),[0 5]);
orient(f2.fig,'landscape');
print(f2.fig,'-dpdf','/home/chaosong/Pictures/allbsttarvlnn1dist4th.pdf');
% print(f2.fig,'-dpdf','/home/chaosong/Pictures/hienvccbsttarvlnn1dist4th.pdf');

%%
%%%the separation in arrival time vs. distance along 'propagation' between sources N--N-1, 
%%%for data and noise. The caveat is, not all bursts of data is migrating, let alone
%%%noise. So this plot is NOT very reliable
widin = 6.5;  % maximum width allowed is 8.5 inches
htin = 3.5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f3 = initfig(widin,htin,nrow,ncol); %initialize fig
xran = [0.05 0.96]; yran = [0.12 0.96];
xsep = 0.08; ysep = 0.05;
optaxpos(f3,nrow,ncol,xran,yran,xsep,ysep);

dtarvlnn1all = bstsig.dtarvlnn1all;
distarvlpropall = bstsig.distarvlpropall;
f3 = plt_propdistvsdt(f3,[dtarvlnn1all,distarvlpropall],1,sps,'tarvl','b');

dtarvlpropall = bstnoi.dtarvlpropall;
distarvlpropall = bstnoi.distarvlpropall;
f3 = plt_propdistvsdt(f3,[dtarvlpropall,distarvlpropall],1,sps,'tarvl','r');
orient(f3.fig,'landscape');
print(f3.fig,'-dpdf','/home/chaosong/Pictures/allbsttarvlnn1propdist.pdf');

%%
%%%Same as above, but how about for sources after the 4th sta check
widin = 6.5;  % maximum width allowed is 8.5 inches
htin = 3.5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f4 = initfig(widin,htin,nrow,ncol); %initialize fig
xran = [0.05 0.96]; yran = [0.12 0.96];
xsep = 0.08; ysep = 0.05;
optaxpos(f4,nrow,ncol,xran,yran,xsep,ysep);

dtarvlnn14thall = bstsig.dtarvlnn14thall;
distarvlprop4thall = bstsig.distarvlprop4thall;
f4 = plt_propdistvsdt(f4,[dtarvlnn14thall,distarvlprop4thall],1,sps,'tarvl','b');

dtarvlprop4thall = bstnoi.dtarvlprop4thall;
distarvlprop4thall = bstnoi.distarvlprop4thall;
f4 = plt_propdistvsdt(f4,[dtarvlprop4thall,distarvlprop4thall],1,sps,'tarvl','r');
orient(f4.fig,'landscape');
print(f4.fig,'-dpdf','/home/chaosong/Pictures/allbsttarvlnn1propdist4th.pdf');

%%
%%%how does the distribution of off14 (LFE caalog) compare to that from the tremor catalog?
%%%NOTE as of 2022/12/22, due to the prealignment between traces 'off1i', it is not simply the direct
%%%subtraction between the arrival time difference at station pairs. A recalculation may be required
%%%depending on if 'off1i' has been accounted. 
impindepstall = bstsig.impindepstall;
nsrc = bstsig.nsrc;
% off1ic = zeros(length(nsrc),2);
% for i = 1: length(nsrc)
%   off12a = [impindepstall(sum(nsrc(1:i)),1)-impindepstall(sum(nsrc(1:i)),3) ...
%     impindepstall(sum(nsrc(1:i)),1)-impindepstall(sum(nsrc(1:i)),5)];
%   off12 = [impindepstall(sum(nsrc(1:i)),7) impindepstall(sum(nsrc(1:i)),8)];
%   off1ic(i,1:2) = off12-off12a;
% end
f1.fig = figure;
f1.fig.Renderer = 'painters';
ax=gca; hold(ax,'on'); grid(ax,'on'); ax.Box = 'on';
scatter3(ax,impindepstall(:,7),impindepstall(:,8),impindepstall(:,end-1),15,[.5 .5 .5],'filled',...
  'MarkerEdgeColor','k');
xlabel(ax,sprintf('off12 at %d sps',sps), 'Interpreter', 'none' );
ylabel(ax,sprintf('off13 at %d sps',sps), 'Interpreter', 'none' );
zlabel(ax,sprintf('off14 at %d sps',sps), 'Interpreter', 'none' );
view(-30, 5);


