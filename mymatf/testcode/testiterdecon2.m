% testiterdecon2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to test the iterative deconvolution methods written by me, Chao Song,
% using some real data and real template, segments similar from 
% 'plt_tremor_spectra_ETS'.
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/11/19
% Last modified date:   2021/11/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e   % Set the format to 5-digit floating point
clc
clear
close all


%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, resol] = pixelperinch(1);

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
%%% choose the data path here
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');
figpath = strcat(workpath, '/project2021/PGCtrio');
hypopath = strcat(getenv('MHOME'), '/Seisbasics/hypoinverse/forproj21');

fam = '002';

FLAG = 'PGC';

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
  'KLNB '
  ];
nsta=size(stas,1);         %  number of stations

%load the template
sps = 40;
wlen = 60*sps;
for ista = 1:3
  fname = strcat(temppath, '/',fam,'_',stas(ista,:),'_',num2str(sps),'sps','_',num2str(wlen/sps),...
    's_BBCCS_opt_Nof_Non_Chao');
  tracetemp(:,ista) = load(fname);
end

%% use a real window of seismogram
TYPE='HF';
bndflag='bnd';
  
nfampool = [
  '002';
  '243';
  '240';  % the most recent catalog of PGC is from fam 002, 243 and 240
  ];
nfam = size(nfampool,1);

% read the rotation parameters accroding to the fam
% Note: here I guess we have to choose a controlling fam and use its rotation para to read the data,
% etc.
% get permanent and polaris station rotation parameters
if isequal(fam,'002')
  %     [timoffrot,stas,bostname,tempoffs] = GetDays(fam,freqflag);
  %     nsta=size(stas,1);         %  number of stations
  % get permanent and polaris station rotation parameters
  sft2=0;     % centroid shift of station 2
  sft3=0;     % centroid shift of station 3
  CATA = 'fixed';
  [PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,sft2,sft3);
else
  sft2=0;     % centroid shift of station 2
  sft3=0;     % centroid shift of station 3
  CATA = 'new';
  [PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,sft2,sft3);
  reftime = PERMROTS(1,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
  PERMROTS(:,4) = PERMROTS(:,4)-reftime;    % to make sure that 1st station is 0
  POLROTS(:,4) = POLROTS(:,4)-reftime;
end
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;     % convert 2-3 columns to rad
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;

%%% load the time range of the tremor bursts
%%% option 1, windows from automatic clustering algorithm
fname = strcat(hypopath, '/tremor_burst_ranges_',FLAG,TYPE,bndflag);
trange = load(fname);

%%% load the merged catalog of hf detections, within dist range, no double counting, of tremor
%%% bursts
winlensechf = 4;
winlensec = 16;
loopoffmax = 4;
xcmaxAVEnmin = 0.45;
SUFFIXhf = strcat('up.hf.time.',num2str(winlensechf),'_',num2str(winlensec),'.',...
  num2str(loopoffmax),'.',num2str(xcmaxAVEnmin));
distmaxhf = 10;
distmaxlf = 10;
hffname = strcat(hypopath, '/eloc.pgc.pj21.',num2str(distmaxhf),'dcutnodb.bst',...
  TYPE,bndflag,'.',SUFFIXhf);
hftime = load(hffname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: if using the burst catalog, there are 59 cols where the last col is the occurrence time in
% days w.r.t. the starting day of each ETS, 2003060, 2004194, 2005254
%
% 58 cols, format is:
%   E(002) N(002) E(own) N(own) dist(own) lon lat dep  (8 + previous 50 cols)
%%% UPDATED at 2021/06/23
%%%   n=n+8;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
daycol = 14;
seccol = 16;
hftime = sortrows(hftime, [daycol, seccol]);

% for fairly broadband filtering, will filter again before deconvolution
hi=15;
lo=0.1;
npo=2;
npa=2;

% dates and ets
%%% NOTE: 'dates' is the dates that the tremor is active at the region of interest, so that the
%%% waveform at the stations on these dates you are looking at mainly comes from the sources at the
%%% region of interest. We want to see if there is a noticable change in spectra during the burst
%%% windows on these dates
dates = unique(trange(:,1));
years = unique(floor(dates/1000));
nets = length(years);

% keyboard
iets = 3;
% dates in each ets
year = years(iets);
datesets = dates(floor(dates/1000)==year);

if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
  stas(4,:)='KELB ';
else
  stas(4,:)='KLNB ';  % remember to change it back
end

k = 0;  % for counting
optcatdat = [];   % concataneted array for all windows in the same ETS, optimal component
ortcatdat = [];   % orthogonal component
optsnrcatdat = [];  % optimal component within a freq range of a high SNR
i = 2;
%%% read daily tremor data
date = datesets(i);
jday = floor(date-year*1000);

%read horizontal optimal and orthogonal components
[STAopt,STAort,~,fileflag] = rd_daily_bpdata(year,jday,datapath,stas,PERMSTA,POLSTA,...
  PERMROTS,POLROTS,sps,lo,hi,npo,npa,[],[],[],[]);

%% load the LFE catalog of Michael Bostock on the same date
%format: [fam yyyy mm dd sec dx dy lon lat dep magnitude number-of-stations], 12 cols
catalog = ReformBostock;
temp = jul2dat(year,jday);
mm = temp(1);
dd = temp(2);
bostcat = catalog(catalog(:,2)==year & catalog(:,3)==mm & catalog(:,4)==dd,:);
bostcat = sortrows(bostcat, 5);

%% load the tremor catalog of John Armbruster on the same date
%format: [yyyy mm dd sec dx dy lon lat dep], 9 cols;
catalog = ReformArmbruster;
armcat = catalog(catalog(:,1)==year & catalog(:,2)==mm & catalog(:,3)==dd,:);
armcat = sortrows(armcat, 4);
clear catalog

%%
%%% We do 2 things here:
%%% 1. simply get the spectra of all windows and store them
%%% 2. concatenate all windows from the same ETS as if they are continuous, then get the
%%% spectrogram of it
% get the windows (time ranges) of the same day, for efficiency
rangetemp = trange(trange(:,1)==datesets(i), :);
%     keyboard
j = 4;  % this range comes from Rubin&Armbruster2013 Fig. 6 continued lower left, 6g 
% j = 11;  % this range comes from Rubin&Armbruster2013 Fig. 6 continued lower right, 6h
rangeplt = rangetemp(j,:);

%%%%%%%%%%%%%%%%%%%%% for testing purpose %%%%%%%%%%%%%%%%%%%%%5
% tst = 16.37*3600; % this range comes from Rubin&Armbruster2013 Fig. 6 continued lower right 
% ted = 16.47*3600;
% tst = 32600; % this range does not contain any burst from merged catalog of 002, 243, and 240
% ted = 33000;
% tst = 25600; % this range does not contain any burst from merged catalog of 002, 243, and 240
% ted = 26000;
% rangeplt = [2005255 tst ted];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%%% plot the bursts in map       	
f.fig=figure;
widin = 5;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*resol htin*resol]);
nrow = 1;
ncol = 1;
for isub = 1:nrow*ncol
  f.ax(isub) = subplot(nrow,ncol,isub);
end
colnum = [daycol seccol 1 2];
[f.ax(1),hf] = plt_detections_inmap(f.ax(1),hftime,rangeplt,rangeplt(2:3),colnum);

% choose the start and end time of data
% tst = rangetemp(j,2);
% ted = rangetemp(j,3);
timetap = 10; % buffer time intended for tapering
tst = 9.9*3600-timetap; % this range focus on 002 surroundings
ted = 10.0*3600+timetap;
ist = floor(tst*sps);
ied = floor(ted*sps-1);
cla(f.ax(1));
[f.ax(1),~] = plt_detections_inmap(f.ax(1),hftime,[rangeplt(1) tst ted],rangeplt(2:3),colnum);

%time offsets between stations
% off12 = 0;
% off13 = 0;
hfobj = hf(hf(:,seccol)>=tst & hf(:,seccol)<=ted, :); % object hf tremor detections

%get a rough sense on 'noise' from other LFE/tremor from other regions to the target area 
bnduse = [min(hfobj(:,1)) min(hfobj(:,2));
          max(hfobj(:,1)) min(hfobj(:,2));
          max(hfobj(:,1)) max(hfobj(:,2));
          min(hfobj(:,1)) max(hfobj(:,2));
          min(hfobj(:,1)) min(hfobj(:,2));
          ];
hold(f.ax(1),'on');
plot(f.ax(1),bnduse(:,1),bnduse(:,2),'k-','linew',2);
hold(f.ax(1),'off');

[iin,ion] = inpolygon(armcat(:,5),armcat(:,6),bnduse(:,1),bnduse(:,2));
isinbnd = iin | ion;
% armobj = armcat(armcat(:,4)>=tst & armcat(:,4)<=ted, :);
armobj = armcat(armcat(:,4)>=tst & armcat(:,4)<=ted & isinbnd~=1, :);  % object armbruster's tremor detections

[iin,ion] = inpolygon(bostcat(:,6),bostcat(:,7),bnduse(:,1),bnduse(:,2));
isinbnd = iin | ion;
% bostobj = bostcat(bostcat(:,5)>=tst & bostcat(:,5)<=ted, :);
bostobj = bostcat(bostcat(:,5)>=tst & bostcat(:,5)<=ted & isinbnd~=1, :);  % object bostock's lfe detections

%should account for time offsets between stations here
%flag for how to pre-align your records before computing running CC, deconvolution, etc
PREALIGN = 'median';  % use the median offset from detections within the time range
if strcmp(PREALIGN, 'median')
  off12med = median(hfobj(:, 11))*sps; % median offset from detections
  off13med = median(hfobj(:, 12))*sps;
  off12 = round(off12med);
  off13 = round(off13med);
  
  optdat = [];
  ortdat = [];
  optdat(:, 1) = STAopt(max(ist,1): min(ied,86400*sps), 2);
  ortdat(:, 1) = STAort(max(ist,1): min(ied,86400*sps), 2);
  optdat(:, 2) = STAopt(max(ist-off12,1): min(ied-off12,86400*sps), 3);
  ortdat(:, 2) = STAort(max(ist-off12,1): min(ied-off12,86400*sps), 3);
  optdat(:, 3) = STAopt(max(ist-off13,1): min(ied-off13,86400*sps), 4);
  ortdat(:, 3) = STAort(max(ist-off13,1): min(ied-off13,86400*sps), 4);

elseif strcmp(PREALIGN, 'constrained')
  %maybe a median offset for a long migration is too rough, better to separate into shorter segments
  %and plot, align them separately, and compute the running CC
  %also, the median offset may not be enclosed
  %%%%%%%%%%% independent enclosed alignment %%%%%%%%%%%%%%%%%%
  %first align the broadband trace
  mshiftadd = sps/8;    % maximum allowed shift between 2 traces
  buffer = 2*mshiftadd; % need some buffer for shifting
  opttmp = STAopt(max(ist-buffer,1): min(ied+buffer,86400*sps), 2:end); 
  orttmp = STAort(max(ist-buffer,1): min(ied+buffer,86400*sps), 2:end); 
  mid = ceil(size(opttmp,1)/2);
  fixlen = ied-ist+1;
  loffmax = 5;
  ccmin = 0.01;  % depending on the length of trace, cc could be very low
  iup = 1;    % times of upsampling
  [off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(opttmp',mid,...
    fixlen,mshiftadd,loffmax,ccmin,iup);
  off12 = round(off12con);
  off13 = round(off13con);
  
  % align the records
  istart = buffer+1;
  iend = istart+ied-ist;
  optdat = [];
  ortdat = [];
  optdat(:, 1) = opttmp(istart: iend, 1);
  ortdat(:, 1) = orttmp(istart: iend, 1);
  optdat(:, 2) = opttmp(istart-round(off12): iend-round(off12), 2);
  ortdat(:, 2) = orttmp(istart-round(off12): iend-round(off12), 2);
  optdat(:, 3) = opttmp(istart-round(off13): iend-round(off13), 3);
  ortdat(:, 3) = orttmp(istart-round(off13): iend-round(off13), 3);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% log-scale envelope of the high-freq seismogram, like fig.2 in Tape_etal2018NatGeo
hisig=8;
losig=2;
temp = Bandpass(detrend(optdat(:,2)), sps, losig, hisig, npo, npa, 'butter');
%detrend again for caution
temp=detrend(temp);

[envup, envlo] = envelope(temp);
envup = moving(envup,40);
timemax = (hfobj(:, seccol) -tst)*sps +1; % time of max power rate of half sec
figure
subplot(2,1,1);
plot(temp); hold on;
scatter(timemax, 0.18*ones(size(timemax)),20,'y','filled','MarkerEdgeColor','k');
xlim([6000 16000]);%10000
subplot(2,1,2);
plot(log10(envup)); hold on
% plot(envup); hold on
scatter(timemax, 0.18*ones(size(timemax)),20,'y','filled','MarkerEdgeColor','k');
xlim([6000 16000]);
  
%%
stasig = [];
stawlet = [];
for ista = 1: 3
  
  sigbb = optdat(:,ista);
  % sigbb = sigbb(1:50*sps);
  lsig = length(sigbb);
  
%   figure
%   plot(tracetemp(:, ista));
  
  lwlet = pow2(9);
  zerocs = 1198;
  %according analysis to the template, the zero-crossing to the end of coda is safe to set as 8 s
  wletbb = tracetemp(zerocs+8*sps-lwlet+1: zerocs+8*sps, ista);
  
  %% remove mean and linear trend and taper
  %of signal
  sigbbd = detrend(sigbb);
  %and taper with tukeywin, which is actually a tapered cosine window
  %tapered length is adaptative to frequency, maybe at least longer than one full period length of
  %the lowest frequency
  ptstap = timetap*sps; % points tapered at the start or end 
  fractap = 2*ptstap/size(sigbbd,1); % if fractap is >=1, n-point von Hann window is returned
  w = tukeywin(size(sigbbd,1),fractap);
  sigbbdt = w.* sigbbd;
  %detrend again for caution
  sigbbdt=detrend(sigbbdt);

  %romve mean, linear trend of template
  wletbbd = detrend(wletbb);
  %and taper with tukeywin, which is actually a tapered cosine window
  w = tukeywin(size(wletbbd,1),fractap);
  wletbbdt = w.* wletbbd;
  %detrend again for caution
  wletbbdt=detrend(wletbbdt);
  
  %% filter the signal and template to reach a similar spectra shape
  %%% 1.8-4.5 for signal and 1.8-6.0 for template seem like a good pair
  %%% OR, 1.8-5.4 for signal and 1.8-18 for template to keep the most of the template at the high end
  %%% OR, use 1.8-4.5 for both, but use more poles (eg., 3) for signal to achieve a faster decay
  %filter the signal
  hisig=5.4;
  losig=1.8;
  sig = Bandpass(sigbbdt, sps, losig, hisig, npo, npa, 'butter');
  %detrend again for caution
  sig=detrend(sig);
  
  %filter the template
  hiwlet=18;
  lowlet=1.8;
  wlet = Bandpass(wletbbdt, sps, lowlet, hiwlet, npo, npa, 'butter');
  %detrend again for caution
  wlet=detrend(wlet);

  %store them for each station
  stasig(:,ista) = sig;
  stawlet(:,ista) = wlet;
  
  %%% plot the original and bandpassed waveforms and spectra
  f.fig = figure;
  widin = 5;  % maximum width allowed is 8.5 inches
  htin = 8.5;   % maximum height allowed is 11 inches
  set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);
  
  subplot(611)
  hold on
  box on
  grid on
  plot(1:lwlet,wletbb,'r','linew',1);
  xlim([0 lwlet]);
  legend('Broadband template','location','southeast');
  text(0.9,0.9,stas(ista,:),'unit','normalized');
  xlabel('Samples');
  ylabel('Amplitude');
  
  subplot(612)
  hold on
  box on
  grid on
  plot(1:lsig,sigbb,'k','linew',1);
  xlim([0 lsig]);
  legend(sprintf('Broader-band signal, %.1f-%.1f Hz',lo,hi),'location','southeast');
  text(0.9,0.9,stas(ista,:),'unit','normalized');
  xlabel('Samples');
  ylabel('Amplitude');
  
  subplot(613)
  hold on
  box on
  grid on
  plot(1:lwlet,wlet,'r','linew',1);
  xlim([0 lwlet]);
  legend(sprintf('Bandpassed template, %.1f-%.1f Hz',lowlet,hiwlet),'location','southeast');
  text(0.9,0.9,stas(ista,:),'unit','normalized');
  xlabel('Samples');
  ylabel('Amplitude');
  
  subplot(614)
  hold on
  box on
  grid on
  plot(1:lsig,sig,'k','linew',1);
  xlim([0 lsig]);
  legend(sprintf('Bandpassed signal, %.1f-%.1f Hz',losig,hisig),'location','southeast');
  text(0.9,0.9,stas(ista,:),'unit','normalized');
  xlabel('Samples');
  ylabel('Amplitude');
  
  subplot(6,1,[5 6]);
  %%% plot the spectrum, using 'pchave', which is capable of time-averaging a longer segment
  wlen = lwlet;
  polap = 50;
  nfft = wlen;
  Fs = sps;
  %note that i request for only 1 window, so there is no robust estimate for the template, it is
  %pretty much the same as 'periodogram', CHECK IT
  [~,~,~,~,~,~,pcwletbb]=pchave(wletbb,wlen,polap,nfft,Fs,'MAD','dpss');
  [pcsigbb,pcft]=pchave(sigbb,wlen,polap,nfft,Fs,'MAD','dpss');
  
  ax = gca;
  loglog(ax,pcft,sqrt(pcwletbb),'r--','linew',1);
  hold(ax,'on');
  % grid(ax,'on');
  loglog(ax,pcft,sqrt(pcsigbb),'k--','linew',1);
  
  [~,~,~,~,~,~,pcwlet]=pchave(wlet,wlen,polap,nfft,Fs,'MAD','dpss');
  [pcsig,pcft]=pchave(sig,wlen,polap,nfft,Fs,'MAD','dpss');
  
  loglog(ax,pcft,sqrt(pcwlet),'color','r','linew',1.5);
  hold(ax,'on');
  % grid(ax,'on');
  loglog(ax,pcft,sqrt(pcsig),'color','k','linew',1.5);
  yran = [1e-5 1e1];
  plot(ax,[1 1],yran,':','color',[0.7 0.7 0.7]);
  plot(ax,[2 2],yran,':','color',[0.7 0.7 0.7]);
  plot(ax,[4 4],yran,':','color',[0.7 0.7 0.7]);
  plot(ax,[8 8],yran,':','color',[0.7 0.7 0.7]);
  plot(ax,[losig losig],yran,'-.','color',[0.5 0.5 0.5]);
  plot(ax,[hisig hisig],yran,'-.','color',[0.5 0.5 0.5]);
  plot(ax,[lowlet lowlet],yran,'-.','color',[1 0.65 0]);
  plot(ax,[hiwlet hiwlet],yran,'-.','color',[1 0.65 0]);
  xticks(ax,[0.1 1 10]);
  xlim(ax,[0.1 10]);
  ylim(ax,yran);
  %     ax.GridLineStyle = '-';
  lgd = legend(ax,'Broadband template',sprintf('Broader-band signal, %.1f-%.1f Hz',lo,hi),...
    sprintf('Bandpassed template, %.1f-%.1f Hz',lowlet,hiwlet),...
    sprintf('Bandpassed signal, %.1f-%.1f Hz',losig,hisig),'location','northeast','fontsize',8);
  % make background transparent
  set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));
  xlabel(ax,'Frequency');
  ylabel(ax,'Amplitude');
  hold(ax,'off');

end


%% get the running average CC between 2 stations
cclen=sps*4;
[ircc,rcc12] = RunningCC(stasig(:,1), stasig(:,2), cclen);
[~,rcc13] = RunningCC(stasig(:,1), stasig(:,3), cclen);
[~,rcc23] = RunningCC(stasig(:,2), stasig(:,3), cclen);
rcc = (rcc12+rcc13+rcc23)/3;

%% convolve the original signal with template of other station, and get the running average of it
dt = 1/sps;  % sampling interval
twlet = [190/sps; 191/sps; 191/sps];  % Time shift of main arrival of wavelet
s1w2 = conv(stasig(:,1),stawlet(:,2),'full'); % the order doesn't matter
s1w2 = s1w2(1+round(twlet(2)/dt):lsig+round(twlet(2)/dt));  % cut accordingly
s1w3 = conv(stasig(:,1),stawlet(:,3),'full');
s1w3 = s1w3(1+round(twlet(3)/dt):lsig+round(twlet(3)/dt));
s2w1 = conv(stasig(:,2),stawlet(:,1),'full');
s2w1 = s2w1(1+round(twlet(1)/dt):lsig+round(twlet(1)/dt));  % cut accordingly
s2w3 = conv(stasig(:,2),stawlet(:,3),'full');
s2w3 = s2w3(1+round(twlet(3)/dt):lsig+round(twlet(3)/dt));
s3w1 = conv(stasig(:,3),stawlet(:,1),'full');
s3w1 = s3w1(1+round(twlet(1)/dt):lsig+round(twlet(1)/dt));  % cut accordingly
s3w2 = conv(stasig(:,3),stawlet(:,2),'full');
s3w2 = s3w2(1+round(twlet(2)/dt):lsig+round(twlet(2)/dt));  % cut accordingly
cclen=sps*4;
[ircccon,rcc12c12] = RunningCC(s1w2, s2w1, cclen);  % running cc of 12 between conv of 1 and 2
[~,rcc13c13] = RunningCC(s1w3, s3w1, cclen);
[~,rcc23c23] = RunningCC(s2w3, s3w2, cclen);
rcccon = (rcc12c12+rcc13c13+rcc23c23)/3;

%%%If you want to compare between 3 stations, you have to convolve the record at one station with 
%%%both templates at the other 2 stations
s1w23 = conv(s1w2,stawlet(:,3),'full'); 
s1w23 = s1w23(1+round(twlet(3)/dt):lsig+round(twlet(3)/dt));  % cut accordingly
s2w13 = conv(s2w1,stawlet(:,3),'full');
s2w13 = s2w13(1+round(twlet(3)/dt):lsig+round(twlet(3)/dt));  % cut accordingly
s3w12 = conv(s3w1,stawlet(:,2),'full');
s3w12 = s3w12(1+round(twlet(2)/dt):lsig+round(twlet(2)/dt));  % cut accordingly
cclen=sps*4;
[irccc123,rcc12c123] = RunningCC(s1w23, s2w13, cclen); % running cc of 12 between conv of 1, 2 and 3
[~,rcc13c123] = RunningCC(s1w23, s3w12, cclen);
[~,rcc23c123] = RunningCC(s2w13, s3w12, cclen);
rccc123 = (rcc12c123+rcc13c123+rcc23c123)/3;


%% plot the processed signal in time domain in detail
f.fig = figure;
f.fig.Renderer = 'painters';
widin = 18;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);

ncol = 2;
sublen = 30;  % window length for each subplot 
nrow = ceil((ted-tst)/ncol/sublen);
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
%     grid(f.ax(isub),'on');
end

pltxran = [0.03 0.98]; pltyran = [0.05 0.98];
pltxsep = 0.02; pltysep = 0.04; 
axpos = optaxposition(nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
for isub = 1:nrow*ncol
  set(f.ax(isub), 'position', axpos(isub,:));
end

timemax = (hfobj(:, seccol) -tst)*sps +1; % time of max power rate of half sec
timecnt = (hfobj(:, seccol-1) -tst)*sps +1; % time of center of the 4-s detecting window
off12all = round(hfobj(:, 11)*sps); % offsets of 4-s tremor detections
off13all = round(hfobj(:, 12)*sps);

timearm = round((armobj(:,4) -tst)*sps +1); % time of armbruster's tremor detections
timebost = round((bostobj(:,5) -tst)*sps +1); % time of armbruster's tremor detections

for isub = 1:nrow*ncol
  ax = f.ax(isub);
  hold(ax,'on');
  plot(ax,stasig(:,1),'r-');
  plot(ax,stasig(:,2),'b-');
  plot(ax,stasig(:,3),'k-'); %,'linew',0.5
  xran = [(isub-1)*sublen isub*sublen]*sps;
  yran(2) = max(max(stasig));
  yran(1) = min(min(stasig));
  yloc = yran(1)+(yran(2)-yran(1))*14/15;
  ind = find(timemax>=xran(1) & timemax<=xran(2));
  for i = 1: length(ind)
    barst = timecnt(ind(i))-sps*2+1;
    bared = timecnt(ind(i))+sps*2;
    plot(ax,[barst bared],[yloc yloc],'c-','linew',3);
  end
  scatter(ax,timemax(ind), yloc*ones(size(timemax(ind))),20,'y','filled','MarkerEdgeColor','k');
  yloc12 = yran(1)+(yran(2)-yran(1))*1/6;
  yloc13 = yran(1)+(yran(2)-yran(1))*1/12;
  text(ax,timemax(ind),yloc12*ones(size(timemax(ind))),num2str(off12all(ind)),...
    'HorizontalAlignment','right','FontSize',7,'color','b');
  text(ax,timemax(ind),yloc13*ones(size(timemax(ind))),num2str(off13all(ind)),...
    'HorizontalAlignment','right','FontSize',7,'color','k');
  ind = find(timearm>=xran(1) & timearm<=xran(2));
  if ~isempty(ind)
    yloc = yran(1)+(yran(2)-yran(1))*13/15;
    scatter(ax,timearm(ind), yloc*ones(size(timearm(ind))),15,'g','filled','MarkerEdgeColor','k');
  end
  ind = find(timebost>=xran(1) & timebost<=xran(2));
  if ~isempty(ind)
    yloc = yran(1)+(yran(2)-yran(1))*12/15;
    scatter(ax,timebost(ind), yloc*ones(size(timebost(ind))),15,'m','filled','MarkerEdgeColor','k');
  end

  xlim(ax,xran);
  ylim(ax,yran);
end
ax = f.ax((nrow-1)*ncol+1);
hold(ax,'on');
xlabel(ax,'Samples','FontSize',10);
ylabel(ax,'Amplitude','FontSize',10);
% print(f.fig,'-depsc2',strcat('/home/chaosong/Pictures/untitled.eps'));


%% take the portion after the end of coda as the 'noise', 
noisebb = tracetemp(zerocs+8*sps+1: end);
%romve mean, linear trend of template
noisebb = detrend(noisebb);
%and taper with tukeywin, which is actually a tapered cosine window
w = tukeywin(size(noisebb,1),fractap);
noisebb = w.* noisebb;
%filter the template
noise = Bandpass(noisebb, sps, losig, hisig, npo, npa, 'butter');
[coef, lag] = xcorr(noise, wlet, size(noise, 1), 'none'); % master raw CC
% figure;
% plot(coef);
%find the max. of the master raw CC for its location and amplitude
[amp, idx] = max(coef);   % max of raw cc
lagsamp = lag(idx);   % offset in samples
amp = amp/(sum(wlet.^2));   % convert max raw CC to amplitude

%% iterative deconvolution
ista = 3; % which station are we focusing on?
sig = stasig(:,ista);
wlet = stawlet(:,ista);

dt = 1/sps;  % sampling interval
if ista == 1
  twlet = 190/sps;  % Time shift of main arrival of wavelet
else
  twlet = 191/sps;  % Time shift of main arrival of wavelet
end
width = 2.5;  % width for Gaussian filter
dres_min = 0.5;  % tolerance, percentage change in residual per iteration
mfit_min = 5e-1;  % tolerance, norm of the residual/misfit
% nit_max = 2000;  % max numer of iterations 
tdura = 0.4;
nit_max = round((ted-tst)/tdura);  
fpltit = 1;  % plot flag for each iteration
fpltend = 1;  % plot flag for the final iteration
[sigdecon(:,ista),pred,res,dresit,mfitit,ampit,nit,fighdl] = ...
  iterdecon(sig,wlet,noise,dt,twlet,width,dres_min,mfit_min,nit_max,fpltit,fpltend);

ax = fighdl{2}.ax(1);
hold(ax,'on');
text(ax,0.8,0.85,stas(ista,:),'unit','normalized');
hold(ax,'off');

ax = fighdl{2}.ax(3);
hold(ax,'on');
yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*1/4;
scatter(ax,timemax, yloc*ones(size(timemax)),10,'y','filled','MarkerEdgeColor','k');
yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*1/6;
scatter(ax,timearm, yloc*ones(size(timearm)),10,'g','filled','MarkerEdgeColor','k');
yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*1/8;
scatter(ax,timebost, yloc*ones(size(timebost)),10,'m','filled','MarkerEdgeColor','k');
legend(ax,'Signal','Prediction','4-s Tremor','Tremor(Arm)','LFE(Bost)');
hold(ax,'off');

%%
f.fig = figure;
f.fig.Renderer = 'painters';
widin = 12;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);

ncol = 1;
nrow = 3;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
%     grid(f.ax(isub),'on');
end

pltxran = [0.06 0.94]; pltyran = [0.06 0.98];
pltxsep = 0.02; pltysep = 0.04; 
axpos = optaxposition(nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
for isub = 1:nrow*ncol
  set(f.ax(isub), 'position', axpos(isub,:));
end

for isub = 1:nrow*ncol
ax = f.ax(isub);
hold(ax,'on');
ax.Box = 'on';
yyaxis(ax,'left');
ind = find(sigdecon(:,isub)>0);
stem(ax, ind, sigdecon(ind,isub), 'k','MarkerSize',4); xlim(ax, [0,length(sigdecon(:,isub))]);
yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*14/15;
scatter(ax,timemax, yloc*ones(size(timemax)),20,'y','filled','MarkerEdgeColor','k');
yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*13/15;
scatter(ax,timearm, yloc*ones(size(timearm)),20,'g','filled','MarkerEdgeColor','k');
yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*12.5/15;
scatter(ax,timebost, yloc*ones(size(timebost)),20,'m','filled','MarkerEdgeColor','k');
text(ax,0.3,0.93,strtrim(stas(isub,:)),'FontSize',10,'unit','normalized','EdgeColor','k','Margin',2);
if isub == nrow
  xlabel(ax,'Samples');
  ylabel(ax,'Amplitude of deconvolved non-zero impulses');
end
yyaxis(ax,'right');
cctemp = rcc;
cctemp(cctemp<0)=-10;
% plot(ax,irunccn,cctemp,'co','markersize',2);
plot(ax,ircc,rcc,'co','markersize',2);
plot(ax,ax.XLim,[0 0],'--','color',[0.7 0.7 0.7],'linew',1.5);
if isub == nrow
  ylabel(ax,'Average short-win CC between 3 stations');
end
ylim(ax,[-1 1]);
hold(ax,'off');
end


% %% Is there a good way to smooth it? Tried envelope, gaussian filter, not ideal
% %%% maybe check out the paper on sliding average
% gf = gauss_time(20,5);
% sigdecongf = conv(sigdecon,gf,'same');
% figure
% stem(sigdecon); hold on
% plot(sigdecongf);
% [sigdeconenv,~] = envelope(sigdecon,20,'peak');
% plot(sigdeconenv);

% for i = 1: length(sigdecon)
%   temp = zeros(length(sigdecon),1);
%   temp(i) = sigdecon(i);
%   sigdecongf() = sigdecon*
% end
%%
itmat = 1: nit;
figure;
yyaxis left;
plot(itmat,ampit(:,2));
ylabel('Amplitude of deconvolved impulses in each iteration');
yyaxis right;
plot(itmat,rcc(ampit(:,1)-cclen/2),'co','markersize',3);
ylabel('Average short-win CC at the index updated in each iteration');
xlabel('Iteration number');

%%
figure
box on
scatter(ampit(:,2),rcc(ampit(:,1)-cclen/2),20,'k','filled');
ylabel('Average short-win CC at the index updated in each iteration');
xlabel('Amplitude of deconvolved impulses in each iteration');

%%
ampitsort = sort(ampit);
for i = 1: length(ampitsort)
  nle(i,1) = length(ampitsort)-i+1;
end
lognle = log10(nle);
logampit = log10(ampitsort);
f.fig=figure;
nrow = 3;
ncol = 1;
for isub = 1:nrow*ncol
  f.ax(isub) = subplot(nrow,ncol,isub);
end
ax=f.ax(1);
scatter(ax,ampitsort,nle,10,'k','filled');
% xlabel('Amplitude of deconvolved impulses in each iteration');
% ylabel('Number of impulses whose amplitude is no smaller than it');
title(ax,'Linear vs. linear');
ax=f.ax(2);
scatter(ax,ampitsort,lognle,10,'k','filled'); hold(ax,'on');
[fitobj,gof,outp,lognlefit] = linear_bisquare_fit_free(ampitsort,lognle);
plot(ax,ampitsort,lognlefit,'-','linewidth',1.5,'color','b');
% xlabel('Amplitude of deconvolved impulses in each iteration');
% ylabel('Number of impulses whose amplitude is no smaller than it');
title(ax,'logarithmic vs. Linear, indicates a exponential law if a linear trend'); hold(ax,'off');
ax=f.ax(3);
scatter(ax,logampit,lognle,10,'k','filled'); hold(ax,'on');
[fitobj,gof,outp,lognlefit] = linear_bisquare_fit_free(logampit,lognle);
plot(ax,logampit,lognlefit,'-','linewidth',1.5,'color','b'); 
xlabel(ax,'Sorted amplitude of deconvolved impulses');
ylabel(ax,'Number of impulses whose amplitude is no smaller than it');
title(ax,'logarithmic vs. logarithmic, indicates a power law if a linear trend'); hold(ax,'off');















