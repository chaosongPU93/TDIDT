% plt_cc_synth.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --This script aims to plot the variation of CC of the synthetic sesimograms,
% combining the noise-free case as one panel and the single-spot case as
% the other panel. We emphasize only the zero-lag CC, so that 2 types can
% be combined as one plot.
% --This is different from 'cc_synth' or 'cc_syn_onespot' which do more
% complicated tasks and create a plot of RCC and CC together.
%
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/08/07
% Last modified date:   2024/08/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e   % Set the format to 5-digit floating point
clear
clc
close all

%% for easy testing
defval('normflag',0); %whether to normalize templates
defval('rccmwsec',0.5); %moving win len in sec for computing RCC

tdura = 0.25;  %must be consistent with what synthetics actually used
% tdura = 0.4;

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
% if pltflag
set(0,'DefaultFigureVisible','on');
% else
%     set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
% end

[scrsz, resol] = pixelperinch(1);

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');

freqflag='hf';  % flag to indicate whether to do hf or lf;

FLAG = 'PGC'; % detector

fam = '002';   % family number

[timoffrot,~] = GetDays4Stack(fam);
nday = size(timoffrot, 1);

%Use new LFE catalog
CATA = 'new';

% get permanent and polaris station rotation parameters, based on 40-sps data
sft2=0;     % centroid shift of station 2
sft3=0;     % centroid shift of station 3
[PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,sft2,sft3);
if ~isequal(CATA, 'fixed')
  reftime = PERMROTS(1,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
  PERMROTS(:,4) = PERMROTS(:,4)-reftime;    % to make sure that 1st station is 0
  POLROTS(:,4) = POLROTS(:,4)-reftime;
end
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;     % convert 2-3 columns to rad
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
  'LZB'];
POLSTA=['SSIB '           % polaris station names
  'SILB '
  'KLNB '
  'MGCB '
  'TWKB '];

stas=['PGC  '
  'SSIB '
  'SILB '
  %   'LZB  '
  %   'TWKB '
  %   'MGCB '
  'KLNB '
  ]; % determine the trio and order, here the 1st sta is PGC
nsta=size(stas,1);         %  number of stations

sps = 40;

iup = 4;  % upsample 4 times

cutout = 'ellipse';

%load detections
if isequal(fam,'002')
  cyclskip = 0;
  mshift=26+cyclskip; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
  loopoffmax=2.1; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
  xcmaxAVEnmin=0.44; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?
end
hi=6.5;    % frequency band
lo=1.25;
npo=2;     % poles, passes of filters
npa=2;
winlensec=4;
winoffsec=1;        % window offset in sec, which is the step of a moving window
winlen=winlensec*sps;      % length in smaples

%load detections inside the cutout boundary of interest, an output of 'locinterp002_4s.m'
PREFIX = strcat(fam,'.up.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),'.',...
  int2str(npo),int2str(npa),'.ms', int2str(mshift));
hfbnd = load(strcat(rstpath, '/MAPS/eloc.pgc.tdec_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',...
  num2str(winlen/sps),'s',num2str(sps),'sps4add_',cutout(1:4)));
daycol = 14;
seccol = 16;  % 2 choices here, one is the center of detecting window, other is the start of strongest .5 s
hfbnd = sortrows(hfbnd, [daycol, seccol]);
off12ran = minmax(hfbnd(:,9)');
off13ran = minmax(hfbnd(:,10)');

%load all detections, an output of 'locinterp002_4s.m'
hfall = load(strcat(rstpath, '/MAPS/eloc.pgc.tdec_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',...
  num2str(winlen/sps),'s',num2str(sps),'sps4add_'));
hfall = sortrows(hfall, [daycol, seccol]);
%get ones outside your boundary
hfout = setdiff(hfall,hfbnd,'rows');
hfout = sortrows(hfout, [daycol, seccol]);
%format as follows:
%%% 8+34+4*nstanew cols, if 4 new stas, then will be 58 cols
%%% UPDATED at 2021/06/23
%%%   dx,dy,lon,lat,dep,tori,off12,off13 (integer samples at upsampled sps)
%%%   1:off12(n) 2:off13(n) 3:off12sec 4:off13sec 5:fam 6:date 7:timswin(n)
%%%   8:timswin(n)-winlensec/2+idiff/sps 9:xcmaxAVEnbang(nin) 10:loopoff(n)
%%%   11:cumsumtrdiff  12:cumsumtrdiff/cumsumtr(winlen)
%%%   for the rest, n=n+4;
%%%   9:Ampsq(nin) 10:Prior(nin) 11:Post(nin) 12:Prior4(nin) 13:Post4(nin) 14:Prior8(nin)
%%%   15:Post8(nin) 16:Prior16(nin) 17:Post16(nin) 18:cc(nin) 19:ccprior(nin)
%%%   20:ccpost(nin) 21:ccprior4(nin) 22:ccpost4(nin) 23:STAamp(nin,2) 24:STAamp(nin,3)
%%%   25:sumsSTA12n(n,iSTA12bang) 26:sumsSTA13n(n,iSTA13bang)
%%%   27:sumsSTA32n(n,iSTA32bang) 28:sigma(nin) 29:sigma12(nin) 30:sigma13(nin)
%%%   in(1:nstanew) loff(1:nstanew) ioff(1:nstanew) ccmaxave(1:nstanew)

ttol = 35;
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.pgc002.',cutout(1:4)));
tlen = trange(:,3)-trange(:,2);
nbst = size(trange,1);

%%%load the empirical model param of time offset 14 for each addtional stations
off14mod = load(strcat(rstpath, '/MAPS/timeoff_planefitparam_4thsta_160sps'),'w+');

%% load synthetic seismograms
%%%specify distribution for source location
distr='UN';  % uniform distribution

%The ratio of elsewhere events to local events.  Zero for Discrete Ide.
fracelsew=0; %0.25 %0.5; %0.6;

%%%specify regime for transformation from time offset to map location
% ftrans = 'interpArmb';
% ftrans = 'interpArmbreloc';
ftrans = 'interpchao';

%times of saturation
% nsat=[0.1 0.4 1 2 4 10 20 40 100];
nsat=[0.4 1 2 4 10 20 40 100];
% nsat=[0.4 2 10];
nnsat = length(nsat);

%length of each simulation
sps = 160;
greenlen = pow2(9)*sps/40;
bufsec = 1;
msftaddm = bufsec*sps;  %buffer range for later CC alignment, +1 for safety
rccmwsec = 0.5;
rccmwlen = rccmwsec*sps;  %window length for computing RCC
overshoot = rccmwlen/2; %RCC points to the center of computing window, so extra room is needed
% Twin=0.5*3600+(greenlen+msftaddm*2+overshoot*2-2)/sps; %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
% nrun = 6;

Twin=3*3600+(greenlen+msftaddm*2+overshoot*2-2)/sps; %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
nrun = 1;

%% synthetics from different region sizes and saturation levels
%%%specify if considering the physical size of each source
% physicalsize = 1;
physicalsize = 0;

%%%diameter of physical size
if physicalsize
  diam=0.15;% 0.5; %0.6; %
else
  diam=0;
end

%%%specify shape of the source region
srcregion='ellipse';
% srcregion='rectangle';
% srcregion='circle';

%variation of source region size
if strcmp(srcregion,'ellipse')
  semia = 1.75*(0.6:0.2:2.0);
  semib = 1.25*(0.6:0.2:2.0);
  nreg = length(semia);
end
nround = nreg;

% ttstr1 = {'Noise-free syn, '};
% fnsuffix1 = '';

for iround = 1: nround
  for insat = 1: nnsat
    sat = nsat(insat);
    
    xaxis = semia(iround);
    yaxis = semib(iround);
    savefile = strcat('rst_decon_synth_reg',num2str(xaxis),...
      '-',num2str(yaxis),'_nsat',num2str(sat),'_td',num2str(tdura),'.mat'); %,srcregion(1:3),'_'
    load(strcat(workpath,'/synthetics/',savefile));
    
    rcccat{insat,iround} = allsyn.rcccat;
    ccwpair{insat,iround} = allsyn.ccwpair;
    
  end
end

%% synthetics from different noise and saturation levels, sources at a single spot
%different percent of noise
perctrial = 0.1*(0:2:16)';
ntrial = length(perctrial);
nround = ntrial;

% ttstr1 = {'Single-spot syn, '};
% fnsuffix1 = '_onespot';

for iround = 1: nround
  for insat = 1: nnsat
    sat = nsat(insat);
    
    perc = perctrial(iround);
    savefile = strcat('rst_decon_synth_onespot','_noi',num2str(perc),'_nsat',...
      num2str(sat),'_td',num2str(tdura),'.mat');
    load(strcat(workpath,'/synthetics/',savefile));
    
    rcccat1s{insat,iround} = allsyn.rcccat;
    ccwpair1s{insat,iround} = allsyn.ccwpair;
  end
end

%% RCC vs. saturation rate & region size
%%%load the LFE catalog, 25-s-win
% savefile = 'deconv_stats4th_allbstsig.mat';
savefile = 'deconv_stats4th_no23_allbstsig.mat';
allsig = load(strcat(rstpath, '/MAPS/',savefile));
rccbst = allsig.allbstsig.rccbst; %cat rcc; whole-win rcc; same for noise
ccwpairk = allsig.allbstsig.ccwpairk; %zero-lag CC of 25-s-window
mrccwpairk = allsig.allbstsig.mrccwpairk; %median RCC of 25-s-windows

%%%for RCC, noise catalog contains data and noise, 25-s-win, noise
% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
allnoi = load(strcat(rstpath, '/MAPS/',savefile));
rccbstn = allnoi.allbstnoi.rccbst; %cat rcc; whole-win rcc; same for noise
ccwpairkn = allnoi.allbstnoi.ccwpairk; %zero-lag CC of 25-s-window
mrccwpairkn = allnoi.allbstnoi.mrccwpairk; %median RCC of 25-s-windows

%%%load the LFE catalog, whole-win
% savefile = 'deconv1win_stats4th_allbstsig.mat';
savefile = 'deconv1win_stats4th_no23_allbstsig.mat';
allsig1win = load(strcat(rstpath, '/MAPS/',savefile));
rccbst1win = allsig1win.allbstsig.rccbst; %whole-win rcc; same for noise
ccpair1win = allsig1win.allbstsig.ccpair; %whole-win 0-lag cc paris
mcc1win = allsig1win.allbstsig.mcc; %mean whole-win 0-lag cc

%%%load the RCC and CC value from real data, all bursts
% load(strcat('rcc',num2str(rccmwsec),'.mat'));
templump = [];  % lumped RCC from all bursts, from 25-s windows
templump1 = [];  % lumped RCC from all bursts, from whole-windows
templump2 = []; % lumped zero-lag CC of whole-window from all bursts
templump3 = []; % lumped zero-lag CC of 25-s-window from all bursts
templump4 = [];  % lumped RCC from all bursts, from 25-s windows, NOISE
templump5 = [];  % lumped zero-lag CC from all bursts, from 25-s windows, NOISE
for i = 1: length(trange)
  temp = rccbst{i};
  %%%cat rcc from 25-s-wins
  temp1 = sort(temp(:,1),'ascend');
  temp1(:,2) = (1:size(temp1,1))/size(temp1,1);
  templump = [templump; temp1];  
  
  %%%0-lag cc from 25-s-wins
  temp3 = ccwpairk{i};  %all 25-s subwins, pairs 12, 13, 23, as col
  temp3 = (temp3(:,1)+temp3(:,2))./2; %retain the mean of 12 and 13
  templump3 = [templump3; temp3]; %lump them
  
  %%%rcc from the whole-win
  temp1 = sort(temp(:,2),'ascend');
  temp1(:,2) = (1:size(temp1,1))/size(temp1,1);
  templump1 = [templump1; temp1];
  
  %%%0-lag cc from the whole-win
  temp2 = mcc1win(i); % the mean of 12 and 13 of the whole win,
  templump2 = [templump2; temp2];

  %%%%%%%% for noise
  %%%cat rcc from 25-s-wins
  temp = rccbstn{i};
  temp4 = sort(temp(:,3),'ascend');
  temp4(:,2) = (1:size(temp4,1))/size(temp4,1);
  templump4 = [templump4; temp4];  
  
  %%%0-lag cc from 25-s-wins
  temp5 = ccwpairkn{i};  %all 25-s subwins, pairs 12, 13, 23, as col
  temp5 = (temp5(:,1)+temp5(:,2))./2; %retain the mean of 12 and 13
  templump5 = [templump5; temp5]; %lump them
  %%%%%%%% for noise
end
%%%bin based on the lumped RCC for some 'median' value, 25-s-win version
rccplt = templump(:,1:2); %data
% %%%1. bin x by equal width, then get the median of y
% [xcnt,ycnt] = ranybinx(rccplt(:,1),rccplt(:,2),'median',[],[],-1:0.01:1);
% rccreal = [xcnt,ycnt];
%%%2. bin x by equal number
nbin = 200;
[xbin,indbin,n] = binxeqnum(rccplt(:,1),nbin);
for i = 1: nbin
  ind = indbin{i};
  rccreal(i,1) = median(rccplt(ind,1));
  rccreal(i,2) = median(rccplt(ind,2));
end

%%%same as above, but for noise
rccplt = templump4(:,1:2); %data
%%%2. bin x by equal number
nbin = 200;
[xbin,indbin,n] = binxeqnum(rccplt(:,1),nbin);
for i = 1: nbin
  ind = indbin{i};
  rccnoi(i,1) = median(rccplt(ind,1));
  rccnoi(i,2) = median(rccplt(ind,2));
end

%median lumped rcc from all data bursts, 25-s-win version
mrccreal = median(templump(:,1));
%median lumped rcc from all data bursts, whole-win version
mrccreal1win = median(templump1(:,1));
%median lumped 0-lag CC from all data bursts, 25-s-win version
mccreal = median(templump3);
%median lumped 0-lag CC from all data bursts, whole-win version
mccreal1win = median(templump2);
%median lumped rcc from all data bursts, 25-s-win version, NOISE
mrccnoi = median(templump4(:,1));
%median lumped 0-lag CC from all data bursts, 25-s-win version
mccnoi = median(templump5);


%% plot of median RCC/CC value wrt saturation & region size
% rcccat = allsyn.rcccatk;
% ccwpair = allsyn.ccwpairk;
nrow = 1; ncol = 2;
widin = 6;  % maximum width allowed is 8.5 inches
htin = 3.5;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);
pltxran = [0.08 0.98]; pltyran = [0.12 0.98]; % optimal axis location
pltxsep = 0.06; pltysep = 0.05;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log';
p=[]; label=[];
nround = nreg;
color = gradientblue(nround);
for iround = 1: nround
  for insat = 1: nnsat
    aa = ccwpair{insat,iround}; %0-lag CC between 3 pairs from all 25-s wins
    bb = mean(aa(:,[1 2]), 2); %mean pof pairs 12 and 13
    mcc(insat,iround) = median(bb);  %use median
  end
  % p(iround) = plot(ax,nsat,mcc(:,iround),'-','Color',color(iround,:),'linew',1);
  % scatter(ax,nsat,mcc(:,iround),20,color(iround,:),'filled','MarkerEdgeColor','k');
  p(iround) = plot(ax,nsat,mcc(:,iround),'-',...
    'linew',1,'color',color(iround,:),'marker','o','markersize',4,...
    'markerfacec',color(iround,:));
  label{iround} = sprintf('%.1fx%.1f',2*semia(iround),2*semib(iround));
end
p(nround+1) = plot(ax,ax.XLim,[mccreal mccreal],'k--','linew',1.5);
label{nround+1} = 'Data';
lgd=legend(ax,p,label,'NumColumns',2,'Location','best','fontsize',6);
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
lgdtit = 'Region size (km)';
title(lgd,lgdtit,'fontsize',7);
xlabel(ax,'Saturation');
% ylabel(ax,'Zero-lag CC');
ylabel(ax,'Median of the maximum CC');
longticks(ax,2);
ylim(ax,[0 1]);
xticks(ax,nsat);
text(ax,0.02,0.95,'a','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');


ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log';
p=[]; label=[];
nround = ntrial;
color = gradientblue(nround);
for iround = 1: nround
  for insat = 1: nnsat
    aa = ccwpair1s{insat,iround}; %0-lag CC between 3 pairs from all 25-s wins
    bb = mean(aa(:,[1 2]), 2); %mean pof pairs 12 and 13
    mcc1s(insat,iround) = median(bb);  %use median
  end
  p(iround) = plot(ax,nsat,mcc1s(:,iround),'-',...
    'linew',1,'color',color(iround,:),'marker','o','markersize',4,...
    'markerfacec',color(iround,:));
  label{iround} = sprintf('%.1f',perctrial(iround));
end
p(nround+1) = plot(ax,ax.XLim,[mccreal mccreal],'k--','linew',1.5);
label{nround+1} = 'Data';
p(nround+2) = plot(ax,ax.XLim,[mccnoi mccnoi],'r--','linew',1.5);
label{nround+2} = 'Noise';
lgd=legend(ax,p,label,'NumColumns',3,'Location','best','fontsize',6);
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
lgdtit = 'Noise level';
title(lgd,lgdtit,'fontsize',7);
xlabel(ax,'Saturation');
% ylabel(ax,'Median of the maximum CC');
longticks(ax,2);
ylim(ax,[0 1]);
xticks(ax,nsat);
text(ax,0.02,0.95,'b','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
% nolabels(ax,2);
hold(ax,'off');

fname = strcat('cc_syn_comb.pdf');
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));
