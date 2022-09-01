% testiterdecon_groupm20220224
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to test the iterative deconvolution methods written by me, Chao Song,
% using some real data and real template, segments similar from 
% 'plt_tremor_spectra_ETS'.
%
% --Different from 'testiterdecon3' which looks at a segment of a RTM,
%   This version we are focusing on a even shorter record, found by Allan, in 
%   which the start and end of the signal have a low amplitude, indicating a
%   low noise leve (possibly no interference of tremor from elsewhere) and decent
%   SNR. Check the record from Allan's sending. 
% --What we want to know is, what is the minimum number of LFEs (templates) 
%   could possibly be in this segment of record.
% --Use the start and end as the proxy of noise; Use various kinds of 
%   deconvolution techniques.
% --If do not want a taper for a short segment of record, given the potential
%   artifact of filtering due to not tapering, one strategy is to filter the daily
%   data directly to desired passband after tapering, but this would lose the
%   the flexibility of adjusting the corners to ensure a similar spectral shape
%   between the record and template. The other one is to taper and filter a slightly
%   longer segment containing the target record, and cut the record and noise
% --Now is compatible with 40, 80 and 100 sps. The main difference related to 
%   the sampling rate is the zero-crossing time of templates computed. Might be 
%   unified later.
% --80 Hz is still higher than the temporary station's original sampling rate 
%   40 Hz, assuming the corresponding templates are computed.
% --The broader-band signal has a plenty of long-period energy, so the 
%   prealignment based on the enclosed CC between the chosen record segment
%   should always be done AFTER the filtering.            
% --The prealignment between all stations has a rule for the sign. 
%   if off12 >0, sta 1 arrives later than sta 2, so you need to move 2 to
%   the right --> to align with 1. 'ampit' contains the indices of the 
%   impulses at each station and the offsets between each triplets. The
%   sign of the offset is set so if off12 >0, it means sta 1 arrives later
%   than sta 2, so you need to move 2 to the right --> to align with 1.
%   This means, it has the SAME sign rule as the prealignment between the
%   signals at 3 stations. Therefore, the net offset relative to the
%   origin (0,0) is the sum of the 2 parts.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/02/24
% Last modified date:   2022/02/24
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

FLAG = 'PGC'; % use pgc trio

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
%   'KLNB '
  ];
nsta=size(stas,1);         %  number of stations

%load the template
sps = 100;

templensec = 60;
if sps==40 || sps==80
  CATA = 'fixed';
elseif sps == 100
  CATA = 'new';
end

if isequal(fam,'002') && isequal(CATA,'old')
  suffix = '_catold';
elseif isequal(fam,'002') && isequal(CATA,'new')
  suffix = '_catnew';
else
  suffix = [];
end
for ista = 1:nsta
  fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
    num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao',suffix);
  tracetemp(:,ista) = load(fname);
end
figure
for ista = 1:nsta
  plot(tracetemp(:,ista)+0.3*ista); hold on
  text(size(tracetemp,1)*0.45, 0.3*ista, stas(ista,:));
end

%% use a real window of seismogram
TYPE='HF';
bndflag='bnd';
  
% read the rotation parameters accroding to the fam
% Note: here I guess we have to choose a controlling fam and use its rotation para to read the data,
% etc.
% get permanent and polaris station rotation parameters
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
seccol = 16;  % 2 choices here, one is the center of detecting window, other is the start of strongest .5 s 
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
iets = 1;
% dates in each ets
year = years(iets);
datesets = dates(floor(dates/1000)==year);

EXAMP = 2;  % which example to look at

if EXAMP == 1
  i = 2;
elseif EXAMP == 2
  i = 3;
end
%%% read daily tremor data
date = datesets(i);
jday = floor(date-year*1000);

%read horizontal optimal and orthogonal components
JDAY = num2zeropadstr(jday,3);
MO=day2month(jday,year);     % EXTERNAL function, day2month, get the month of one particular date
yr=num2str(year); 
direc=[datapath, '/arch', yr,'/',MO,'/'];     % directory name
prename=[direc,yr,'.',JDAY,'.00.00.00.0000.CN'];    %  path plus prefix of data file,
disp(prename);
[STAopt,STAort,~,fileflag] = rd_daily_bpdata(year,jday,prename,stas,PERMSTA,POLSTA,...
  PERMROTS,POLROTS,sps,lo,hi,npo,npa,[],[],[],[]);

%read vertical components too, as we want to have a sense of the SNR at the orthogonal and vertical
%components
[STAvert,~,~,fileflag] = rd_daily_bpdata(year,jday,prename,stas,PERMSTA,POLSTA,...
  PERMROTS,POLROTS,sps,lo,hi,npo,npa,[],[],[],[],'Z');

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

%% See if the target time windows are inside any tremor burst
%%% Allan's examples:
%%% EX1: 2003062, 85058-85082 s, loud and short-duration signal with a quiet start and end,
%%% according to the amplitude of the waveform, and the running CC variation, it seems reasonable to
%%% deconvolve the 85066-85078 s, and take 85062-85066 s as the noise
%%% EX2: 2003063, 65634-65654 s, medium amp and long-duration signal with a quiet start and end
rangetemp = trange(trange(:,1)==datesets(i), :);
%     keyboard
if EXAMP == 1 
  j = 21;  % EX1 seems inside the 21st burst of day 2003062
elseif EXAMP == 2
  j = 11;  % EX2 seems inside the 11st burst of day 2003063
end
rangeplt = rangetemp(j,:);

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
if EXAMP == 1
  timetap = 5; % buffer time intended for tapering
  tlst = 85058+timetap-timetap; % time range for a longer record, see why in the start
  tled = 85082-timetap+timetap;
  tnoi = [85062 85066];   % time range for the noise
  tsig = [85066 85078];   % time range for the target signal
elseif EXAMP == 2
  timetap = 5; % buffer time intended for tapering
  tlst = 65630+timetap-timetap; % time range for a longer record, see why in the start
  tled = 65658-timetap+timetap;
  tnoi = [65636 65639];   % time range for the noise
  tsig = [65639 65648];   % time range for the target signal

end
ilst = floor(tlst*sps+1);
iled = floor(tled*sps);
isst = floor(tsig(1)*sps+1);
ised = floor(tsig(2)*sps);
inst = floor(tnoi(1)*sps+1);
ined = floor(tnoi(2)*sps);

cla(f.ax(1));
[f.ax(1),~] = plt_detections_inmap(f.ax(1),hftime,[rangeplt(1) tlst tled],rangeplt(2:3),colnum);

hfobj = hf(hf(:,seccol)>=tlst & hf(:,seccol)<=tled, :); % object hf tremor detections

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
armobj = armcat(armcat(:,4)>=tlst & armcat(:,4)<=tled & isinbnd~=1, :);  % object armbruster's tremor detections

[iin,ion] = inpolygon(bostcat(:,6),bostcat(:,7),bnduse(:,1),bnduse(:,2));
isinbnd = iin | ion;
% bostobj = bostcat(bostcat(:,5)>=tst & bostcat(:,5)<=ted, :);
bostobj = bostcat(bostcat(:,5)>=tlst & bostcat(:,5)<=tled & isinbnd~=1, :);  % object bostock's lfe detections


%% simple processing to the different components of the signal, then plot them
opttmp = [];
orttmp = [];
verttmp = [];
for ista = 1: nsta
%   hisig=6.5; % this is the previously adopted HF detecting band
%   losig=1.25;
  hisig=6.3; % this will give a similar spectral shape between template and signal 
  losig=1.8;
  
  mshiftadd = round(sps/8);    % maximum allowed shift between 2 traces
  buffer = 2*mshiftadd;
  
  %optimal component
%   tmp = optdat(:,ista);
  tmp = STAopt(max(ilst-buffer,1): min(iled+buffer,86400*sps), ista+1);
  tmp = detrend(tmp); % remove mean and linear trend
  fractap = sps/size(tmp,1); % if fractap is >=1, n-point von Hann window is returned
  tmp = tukeywin(size(tmp,1),fractap).* tmp; % taper
  tmp = Bandpass(tmp, sps, losig, hisig, npo, npa, 'butter');
  tmp = detrend(tmp);   %detrend again for caution
  opttmp(:,ista) = tmp;
  
  %orthogonal component
%   tmp = ortdat(:,ista);
  tmp = STAort(max(ilst-buffer,1): min(iled+buffer,86400*sps), ista+1);
  tmp = detrend(tmp); % remove mean
  tmp = tukeywin(size(tmp,1),fractap).* tmp; % taper
  tmp = Bandpass(tmp, sps, losig, hisig, npo, npa, 'butter');
  tmp = detrend(tmp);   %detrend again for caution
  orttmp(:,ista) = tmp;
  
  %vertical component
%   tmp = vertdat(:,ista);
  tmp = STAvert(max(ilst-buffer,1): min(iled+buffer,86400*sps), ista+1);
  tmp = detrend(tmp); % remove mean
  tmp = tukeywin(size(tmp,1),fractap).* tmp; % taper
  tmp = Bandpass(tmp, sps, losig, hisig, npo, npa, 'butter');
  tmp = detrend(tmp);   %detrend again for caution
  verttmp(:,ista) = tmp;

end

%% Pre-alignment of all stations based on filtered data
%%%As a reminder, the sign of off12 means, if off12 >0, sta 1 arrives later than sta 2, so you need
%%%to move 2 to the right --> to align with 1
off1i = zeros(nsta,1);
PREALIGN = 'constrained';
if strcmp(PREALIGN, 'median')
  off12med = median(hfobj(:, 11))*sps; % median offset from detections
  off13med = median(hfobj(:, 12))*sps;
  off1i(2) = round(off12med);
  off1i(3) = round(off13med);
  if EXAMP == 1
    off1i(2) = -1*sps/40;
    off1i(3) = -1*sps/40;
  elseif EXAMP == 2
    off1i(2) = 1*sps/40;
    off1i(3) = 1*sps/40; %
  end
    
elseif strcmp(PREALIGN, 'constrained')
  %maybe a median offset for a long migration is too rough, better to separate into shorter segments
  %and plot, align them separately, and compute the running CC
  %also, the median offset may not be enclosed
  %%%%%%%%%%% independent enclosed alignment %%%%%%%%%%%%%%%%%%
  optcc = opttmp(buffer+1+isst-ilst: buffer+1+isst-ilst+ ised-isst, 1:nsta);
%   ortcc = orttmp(buffer+1+isst-ilst: end-buffer-(ised-iled), 1:nsta);
%   vertcc = verttmp(buffer+1+isst-ilst: end-buffer-(ised-iled), 1:nsta);
  ccmid = ceil(size(optcc,1)/2);
  ccwlen = round(size(optcc,1)*0.75);
  loffmax = 5*sps/40;
  ccmin = 0.01;  % depending on the length of trace, cc could be very low
  iup = 1;    % times of upsampling
  [off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(optcc',ccmid,...
    ccwlen,mshiftadd,loffmax,ccmin,iup);
  off1i(2) = round(off12con);
  off1i(3) = round(off13con);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%Align records
optdat = [];
ortdat = [];
vertdat = [];
optdat(:, 1) = opttmp(buffer+1: end-buffer, 1); % sta 1
ortdat(:, 1) = orttmp(buffer+1: end-buffer, 1);
vertdat(:, 1) = verttmp(buffer+1: end-buffer, 1);
optdat(:, 2) = opttmp(buffer+1-off1i(2): end-buffer-off1i(2), 2); % sta 2
ortdat(:, 2) = orttmp(buffer+1-off1i(2): end-buffer-off1i(2), 2);
vertdat(:, 2) = verttmp(buffer+1-off1i(2): end-buffer-off1i(2), 2);
optdat(:, 3) = opttmp(buffer+1-off1i(3): end-buffer-off1i(3), 3); % sta 3
ortdat(:, 3) = orttmp(buffer+1-off1i(3): end-buffer-off1i(3), 3);
vertdat(:, 3) = verttmp(buffer+1-off1i(3): end-buffer-off1i(3), 3);

%%% get the running average CC between 3 stations
cclen=sps/2;
[ircc,rcc12] = RunningCC(optdat(:,1), optdat(:,2), cclen);  % optimal comp.
[~,rcc13] = RunningCC(optdat(:,1), optdat(:,3), cclen);
[~,rcc23] = RunningCC(optdat(:,2), optdat(:,3), cclen);
rccopt = (rcc12+rcc13+rcc23)/3;

[ircc,rcc12] = RunningCC(ortdat(:,1), ortdat(:,2), cclen);  % orthogonal comp.
[~,rcc13] = RunningCC(ortdat(:,1), ortdat(:,3), cclen);
[~,rcc23] = RunningCC(ortdat(:,2), ortdat(:,3), cclen);
rccort = (rcc12+rcc13+rcc23)/3;

[ircc,rcc12] = RunningCC(vertdat(:,1), vertdat(:,2), cclen);  % vertical comp.
[~,rcc13] = RunningCC(vertdat(:,1), vertdat(:,3), cclen);
[~,rcc23] = RunningCC(vertdat(:,2), vertdat(:,3), cclen);
rccvert = (rcc12+rcc13+rcc23)/3;

%% plot the different component
f.fig = figure;
f.fig.Renderer = 'painters';
widin = 16;  % maximum width allowed is 8.5 inches
htin = 2.5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);

ncol = 1;
nrow = 1;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
%     grid(f.ax(isub),'on');
end

pltxran = [0.05 0.95]; pltyran = [0.2 0.9];
pltxsep = 0.02; pltysep = 0.04; 
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
for isub = 1:nrow*ncol
  set(f.ax(isub), 'position', axpos(isub,:));
end

timemax = (hfobj(:, seccol) -tlst)*sps +1; % starting time of max power rate of half sec
timecnt = (hfobj(:, seccol-1) -tlst)*sps +1; % central time of the 4-s detecting window
off12all = round(hfobj(:, 11)*sps); % offsets of 4-s tremor detections
off13all = round(hfobj(:, 12)*sps);

timearm = round((armobj(:,4) -tlst)*sps +1); % time of armbruster's tremor detections
timebost = round((bostobj(:,5) -tlst)*sps +1); % time of armbruster's tremor detections

ax = f.ax(1);
hold(ax,'on');
% yyaxis(ax,'left');
plot(ax,optdat(:,1),'r-');
plot(ax,optdat(:,2),'b-');
plot(ax,optdat(:,3),'k-'); %,'linew',0.5
text(ax,0.05,0.1,'Optimal','unit','normalized',...
  'HorizontalAlignment','left','FontSize',12);
xran = [0 tled-tlst]*sps;
yran(2) = 1.5*max(max(abs(optdat)));
yran(1) = -yran(2);
yloc = yran(1)+(yran(2)-yran(1))*14/15;
ind = find(timemax>=xran(1) & timemax<=xran(2));
for i = 1: length(ind)
  barst = timecnt(ind(i))-sps*2;
  bared = timecnt(ind(i))+sps*2-1;
  plot(ax,[barst bared],[yloc yloc],'c-','linew',3);
end
for i = 1: length(ind)
  barst = timemax(ind(i));
  bared = timemax(ind(i))+sps/2-1;
  plot(ax,[barst bared],[yloc yloc],'b-','linew',4);
end
% scatter(ax,timemax(ind), yloc*ones(size(timemax(ind))),20,'y','filled','MarkerEdgeColor','k');
% yloc12 = yran(1)+(yran(2)-yran(1))*1/6;
% yloc13 = yran(1)+(yran(2)-yran(1))*1/12;
% text(ax,timemax(ind),yloc12*ones(size(timemax(ind))),num2str(off12all(ind)),...
%   'HorizontalAlignment','right','FontSize',7,'color','b');
%   text(ax,timemax(ind),yloc13*ones(size(timemax(ind))),num2str(off13all(ind)),...
%     'HorizontalAlignment','right','FontSize',7,'color','k');
text(ax,0.45,0.2,num2str(off1i(2)),'Units','normalized',...
  'HorizontalAlignment','left','FontSize',10,'color','b');
text(ax,0.45,0.1,num2str(off1i(3)),'Units','normalized',...
  'HorizontalAlignment','left','FontSize',10,'color','k');

% ind = find(timearm>=xran(1) & timearm<=xran(2));
% if ~isempty(ind)
%   yloc = yran(1)+(yran(2)-yran(1))*13/15;
%   scatter(ax,timearm(ind), yloc*ones(size(timearm(ind))),15,'g','filled','MarkerEdgeColor','k');
% end
% ind = find(timebost>=xran(1) & timebost<=xran(2));
% if ~isempty(ind)
%   yloc = yran(1)+(yran(2)-yran(1))*12/15;
%   scatter(ax,timebost(ind), yloc*ones(size(timebost(ind))),15,'m','filled','MarkerEdgeColor','k');
% end
noiran = (tnoi-tlst)*sps;
plot(ax,[noiran(1) noiran(1)],yran,'k--');
plot(ax,[noiran(2) noiran(2)],yran,'k--');
sigran = (tsig-tlst)*sps;
plot(ax,[sigran(1) sigran(1)],yran,'c--');
plot(ax,[sigran(2) sigran(2)],yran,'c--');
text(ax,0.95,0.1,sprintf('Bandpassed by %.2f-%.2f Hz',losig,hisig),'unit','normalized',...
  'HorizontalAlignment','right','FontSize',12);
xlim(ax,xran);
xticks(ax,xran(1): 4*sps: xran(2))
xticklabels(ax,(xran(1)/sps+tlst): 4: (xran(2)/sps+tlst));
ylim(ax,yran);
a = jul2dat(year,jday);
mo = a(1);
if mo == 9
  mo = {' Sep. '};
elseif mo == 7
  mo = {' Jul. '};
else
  mo = {' Mar. '};
end
day = num2str(a(2));
yr = num2str(a(3));
xlabel(ax,strcat({'Time (sec) on '}, day, mo, yr),'FontSize',12);
ylabel(ax,'Amplitude','FontSize',12);
% yyaxis(ax,'right');
% plot(ax,ircc,rccopt,'o','color',[.5 .5 .5],'markersize',2);
% ylabel(ax,'Running CC','FontSize',12);
% ylim(ax,[-1 1]);


%% pre-process the signal and template for a similar spectral shape essential for deconvolution
stasig = [];
% stawlet = [];
stanoi = [];
inoi = (tnoi-tlst)*sps+1;
isig = (tsig-tlst)*sps+1;
  
%%%broader-band signal, optimally shifted, no rtr, no taper 
sigbb = [];
for ista = 1: nsta
  sigbb(:,ista) = STAopt(max(isst-off1i(ista),1): min(ised-off1i(ista),86400*sps), ista+1);
end
lsig = size(sigbb,1);

%%% choose a portion as the proxy of noise, like signal, cut directly
noibb = [];
for ista = 1: nsta
  noibb(:,ista) = STAopt(max(inst-off1i(ista),1): min(ined-off1i(ista),86400*sps), ista+1);
end

%%% use a short segment centered at the main dipole to align the filtered template by constrained CC
tmpwlet = tracetemp;
for ista = 1: nsta
  %romve mean, linear trend of template
  tmpwlet(:,ista) = detrend(tmpwlet(:,ista));
  %and taper with tukeywin, which is actually a tapered cosine window
  w = tukeywin(size(tmpwlet(:,ista),1),fractap);
  tmpwlet(:,ista) = w.* tmpwlet(:,ista);
  %detrend again for caution
  tmpwlet(:,ista)=detrend(tmpwlet(:,ista));
  %filter the template
  hiwlet=18;
  lowlet=1.8;
  tmpwlet(:,ista) = Bandpass(tmpwlet(:,ista), sps, lowlet, hiwlet, npo, npa, 'butter');
  %detrend again for caution
  tmpwlet(:,ista)=detrend(tmpwlet(:,ista));
end
ccmid = round(templensec*sps/2);
ccwlen = 10*sps;
loffmax = 5*sps/40;
ccmin = 0.01;  % depending on the length of trace, cc could be very low
iup = 1;    % times of upsampling
[off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(tmpwlet',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
offwlet1i(1) = 0;
offwlet1i(2) = round(off12con);
offwlet1i(3) = round(off13con);
% figure
% for ista = 1:nsta
%   plot(tmpwlet(:,ista)+0.3*ista); hold on
%   text(size(tmpwlet,1)*0.45, 0.3*ista, stas(ista,:));
% end


%automatically find the rough zero-crossing time, whose abs. value is closest to 0, whether + or -
[~,imin] = min(tmpwlet(:,1));
[~,imax] = max(tmpwlet(:,1));
[~,zcsta1] = min(abs(tmpwlet(imin:imax,1)));
zcsta1 = zcsta1+imin-1;
stawlet = [];
lwlet = pow2(9)*sps/40;
for ista = 1: nsta
  stawlet(:,ista) = tmpwlet(zcsta1+8*sps-lwlet+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  %detrend again for caution
  stawlet(:,ista)=detrend(stawlet(:,ista));

  [~,imin] = min(stawlet(:,ista));
  [~,imax] = max(stawlet(:,ista));
  [~,zcsta(ista)] = min(abs(stawlet(imin:imax,ista)));
  zcsta(ista) = zcsta(ista)+imin-1;
end
%the following is just a check, because now the templates must be best aligned 
ccmid = round(size(stawlet,1)/2);
ccwlen = 4*sps;
loffmax = 5*sps/40;
ccmin = 0.01;  % depending on the length of trace, cc could be very low
iup = 1;    % times of upsampling
[off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(stawlet',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
if ~(off12con==0 && off13con==0)
  disp('Filtered templates are NOT best aligned');
end

figure
for ista = 1:nsta
  plot(stawlet(:,ista)+0.3*ista); hold on
  scatter(zcsta(ista),stawlet(zcsta(ista),ista)+0.3*ista, 20,'k','filled');
  text(size(stawlet,1)*0.1, 0.3*ista, stas(ista,:));
end
xlabel(sprintf('Samples at %d Hz',sps));
ylabel('Amplitude');

figure
plot(stawlet(:,1),'r-'); hold on
plot(stawlet(:,2),'b-');
plot(stawlet(:,3),'k-');
legend(stas);
ylim([-0.3 0.3]);
xlabel(sprintf('Samples at %d Hz',sps),'FontSize',12);
ylabel('Amplitude','FontSize',12);

for ista = 1: nsta
  
  %according analysis to the template, the zero-crossing to the end of coda is safe to set as 8 s
  wletbb(:,ista) = tracetemp(zcsta1+8*sps-lwlet+1-offwlet1i(ista): ...
                             zcsta1+8*sps-offwlet1i(ista), ista);
%   wletbb(:,ista) = tracetemp(zerocs(ista)+8*sps-lwlet+1: zerocs(ista)+8*sps, ista);
  
  %%%instead of detrend, taper and bandpass the desired segment, chop directly from the longer 
  %%%segment that is already bandpassed to the desired freq range
  sig = detrend(optdat(isig(1): isig(2)-1, ista));
  
%   %filter the template
%   hiwlet=18;
%   lowlet=1.8;
%   wlet = Bandpass(wletbbdt, sps, lowlet, hiwlet, npo, npa, 'butter');
%   %detrend again for caution
%   wlet=detrend(wlet);
  wlet = stawlet(:,ista);

  %%%bandpassed noise
  noi = detrend(optdat(inoi(1): inoi(2)-1, ista));

  [f] = plt_wletsig_timefreq(wletbb(:,ista),sigbb(:,ista),wlet,sig,noibb(:,ista),noi,sps,...
                             [lowlet hiwlet],[losig hisig]);
%   [f] = plt_wletsig_timefreq(wletbb,sigbb,wlet,sig,[],[],sps,[lowlet hiwlet],[losig hisig]);
  text(f.ax(1),0.9,0.9,stas(ista,:),'unit','normalized');

  %store them for each station
  stasig(:,ista) = sig;
%   stawlet(:,ista) = wlet;
  stanoi(:,ista) = noi;
  
end

dt = 1/sps;  % sampling interval
statwlet = zcsta/sps; % take the zero-crossing as the time shift of main arrival of wavelet

%%
figure
plot(stawlet(:,2),'r-'); hold on
text(size(stawlet,1)*0.02, 0.15, stas(2,:));
ylim([-0.2 0.2]); xlim([0 1400]);
xlabel(sprintf('Samples at %d Hz',sps),'FontSize',10);
ylabel('Amplitude','FontSize',10);

figure
plot(stasig(:,2),'k-'); hold on
text(size(stawlet,1)*0.02, 0.2, stas(2,:));
ylim([-0.3 0.3]); xlim([0 1400]);
xlabel(sprintf('Samples at %d Hz',sps),'FontSize',10);
ylabel('Amplitude','FontSize',10);


%% running CC
%%%Get the running average CC between 3 stations
cclen=sps/2;
[ircc,rcc12] = RunningCC(stasig(:,1), stasig(:,2), cclen);
[~,rcc13] = RunningCC(stasig(:,1), stasig(:,3), cclen);
[~,rcc23] = RunningCC(stasig(:,2), stasig(:,3), cclen);
rcc = (rcc12+rcc13+rcc23)/3;

%%%Convolve the original signal with template of other station, and get the running average of it
s1w2 = conv(stasig(:,1),stawlet(:,2),'full'); % the order doesn't matter
s1w2 = s1w2(1+round(statwlet(2)/dt):lsig+round(statwlet(2)/dt));  % cut accordingly
s1w3 = conv(stasig(:,1),stawlet(:,3),'full');
s1w3 = s1w3(1+round(statwlet(3)/dt):lsig+round(statwlet(3)/dt));
s2w1 = conv(stasig(:,2),stawlet(:,1),'full');
s2w1 = s2w1(1+round(statwlet(1)/dt):lsig+round(statwlet(1)/dt));  % cut accordingly
s2w3 = conv(stasig(:,2),stawlet(:,3),'full');
s2w3 = s2w3(1+round(statwlet(3)/dt):lsig+round(statwlet(3)/dt));
s3w1 = conv(stasig(:,3),stawlet(:,1),'full');
s3w1 = s3w1(1+round(statwlet(1)/dt):lsig+round(statwlet(1)/dt));  % cut accordingly
s3w2 = conv(stasig(:,3),stawlet(:,2),'full');
s3w2 = s3w2(1+round(statwlet(2)/dt):lsig+round(statwlet(2)/dt));  % cut accordingly
cclen=sps/2;
[ircccon,rcc12c12] = RunningCC(s1w2, s2w1, cclen);  % running cc of 12 between conv of 1 and 2
[~,rcc13c13] = RunningCC(s1w3, s3w1, cclen);
[~,rcc23c23] = RunningCC(s2w3, s3w2, cclen);
rcccon = (rcc12c12+rcc13c13+rcc23c23)/3;

%%%If you want to compare between 3 stations, you have to convolve the record at one station with 
%%%both templates at the other 2 stations
s1w23 = conv(s1w2,stawlet(:,3),'full'); 
s1w23 = s1w23(1+round(statwlet(3)/dt):lsig+round(statwlet(3)/dt));  % cut accordingly
s2w13 = conv(s2w1,stawlet(:,3),'full');
s2w13 = s2w13(1+round(statwlet(3)/dt):lsig+round(statwlet(3)/dt));  % cut accordingly
s3w12 = conv(s3w1,stawlet(:,2),'full');
s3w12 = s3w12(1+round(statwlet(2)/dt):lsig+round(statwlet(2)/dt));  % cut accordingly
cclen=sps/2;
[irccc123,rcc12c123] = RunningCC(s1w23, s2w13, cclen);
[~,rcc13c123] = RunningCC(s1w23, s3w12, cclen);
[~,rcc23c123] = RunningCC(s2w13, s3w12, cclen);
rccc123 = (rcc12c123+rcc13c123+rcc23c123)/3;


%% plot the pre-processed signal in more detail
f.fig = figure;
f.fig.Renderer = 'painters';
% widin = 12;  % maximum width allowed is 8.5 inches
% htin = 2.5;   % maximum height allowed is 11 inches
widin = 8.5;  % maximum width allowed is 8.5 inches
htin = 4;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);

ncol = 1;
nrow = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
%     grid(f.ax(isub),'on');
end

pltxran = [0.1 0.9]; pltyran = [0.15 0.95];
pltxsep = 0.02; pltysep = 0.06; 
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

timemax = (hfobj(:, seccol) -tsig(1))*sps +1; % starting time of max power rate of half sec
timecnt = (hfobj(:, seccol-1) -tsig(1))*sps +1; % central time of the 4-s detecting window
off12all = round(hfobj(:, 11)*sps); % offsets of 4-s tremor detections
off13all = round(hfobj(:, 12)*sps);

timearm = round((armobj(:,4) -tsig(1))*sps +1); % time of armbruster's tremor detections
timebost = round((bostobj(:,5) -tsig(1))*sps +1); % time of armbruster's tremor detections

ax = f.ax(1);
hold(ax,'on');
plot(ax,stawlet(:,1),'r-');
plot(ax,stawlet(:,2),'b-');
plot(ax,stawlet(:,3),'k-');
legend(ax,stas);
text(ax,0.98,0.1,sprintf('Bandpassed %.2f-%.2f Hz',lowlet,hiwlet),'unit','normalized',...
  'HorizontalAlignment','right','FontSize',10);
text(ax,0.02,0.1,'Templates','unit','normalized',...
  'HorizontalAlignment','left','FontSize',10);
ylim(ax,[-0.3 0.3]);
xlim(ax,[0 lwlet]);
% xlabel(sprintf('Samples at %d Hz',sps),'FontSize',11);
ylabel(ax,'Amplitude','FontSize',11);
hold(ax,'off');

%   isub = (ista-1)*2+1;
ax = f.ax(2);
hold(ax,'on');
yyaxis(ax,'left');
plot(ax,stasig(:,1),'r-');
plot(ax,stasig(:,2),'b-');
plot(ax,stasig(:,3),'k-'); %,'linew',0.5
text(ax,0.02,0.1,'Optimal data','unit','normalized',...
  'HorizontalAlignment','left','FontSize',10);
xran = [0 tsig(2)-tsig(1)]*sps;
% yran(2) = 1.6*max(max(abs(stasig)));
yran(2) = 0.6;
yran(1) = -yran(2);
ind = find(timemax>=xran(1) & timemax<=xran(2));
for i = 1: length(ind)
  barst = timecnt(ind(i))-sps*2;
  bared = timecnt(ind(i))+sps*2-1;
  yloc = yran(1)+(yran(2)-yran(1))*(20-i)/20;

  plot(ax,[barst bared],[yloc yloc],'k-','linew',2.5);
end
for i = 1: length(ind)
  barst = timemax(ind(i));
  bared = timemax(ind(i))+sps/2-1;
  yloc = yran(1)+(yran(2)-yran(1))*(20-i)/20;
  plot(ax,[barst bared],[yloc yloc],'b-','linew',4);
end
% scatter(ax,timemax(ind), yloc*ones(size(timemax(ind))),20,'y','filled','MarkerEdgeColor','k');
% yloc12 = yran(1)+(yran(2)-yran(1))*1/6;
% yloc13 = yran(1)+(yran(2)-yran(1))*1/12;
% text(ax,timemax(ind),yloc12*ones(size(timemax(ind))),num2str(off12all(ind)),...
%   'HorizontalAlignment','right','FontSize',7,'color','b');
%   text(ax,timemax(ind),yloc13*ones(size(timemax(ind))),num2str(off13all(ind)),...
%     'HorizontalAlignment','right','FontSize',7,'color','k');
text(ax,0.45,0.2,num2str(off1i(2)),'Units','normalized',...
  'HorizontalAlignment','left','FontSize',10,'color','b');
text(ax,0.45,0.1,num2str(off1i(3)),'Units','normalized',...
  'HorizontalAlignment','left','FontSize',10,'color','k');

% ind = find(timearm>=xran(1) & timearm<=xran(2));
% if ~isempty(ind)
%   yloc = yran(1)+(yran(2)-yran(1))*13/15;
%   scatter(ax,timearm(ind), yloc*ones(size(timearm(ind))),15,'g','filled','MarkerEdgeColor','k');
% end
% ind = find(timebost>=xran(1) & timebost<=xran(2));
% if ~isempty(ind)
%   yloc = yran(1)+(yran(2)-yran(1))*12/15;
%   scatter(ax,timebost(ind), yloc*ones(size(timebost(ind))),15,'m','filled','MarkerEdgeColor','k');
% end
text(ax,0.98,0.1,sprintf('Bandpassed %.2f-%.2f Hz',losig,hisig),'unit','normalized',...
  'HorizontalAlignment','right','FontSize',10);
xlim(ax,xran);
ylim(ax,yran);
xlabel(sprintf('Samples at %d Hz',sps),'FontSize',11);
ylabel(ax,'Amplitude','FontSize',11);

yyaxis(ax,'right');
plot(ax,ircc,rcc,'o','color',[.5 .5 .5],'markersize',2);
ylabel(ax,'Running CC','FontSize',11);
ylim(ax,[-1 1]);
hold(ax,'off');

if lwlet > lsig
  shrink(f.ax(2),lwlet/lsig,1);
  f.ax(2).Position(1)=axpos(1,1);
end

print(f.fig,'-dpdf','/home/chaosong/Downloads/1.pdf');

%% independent iterative deconvolution
for ista = 1:nsta
  % ista = 3 % which station are we focusing on?
  sig = stasig(:,ista);
  wlet = stawlet(:,ista);
  noi = stanoi(:, ista);
  twlet = statwlet(ista);
  
  dt = 1/sps;  % sampling interval
  % if ista == 1
  %   twlet = 190/sps;  % Time shift of main arrival of wavelet
  % else
  %   twlet = 191/sps;  % Time shift of main arrival of wavelet
  % end
  width = 2.5;  % width for Gaussian filter
  dres_min = 0.05;  % tolerance, percentage change in residual per iteration
  mfit_min = 5e-1;  % tolerance, norm of the residual/misfit
  nit_max = 200;  % max numer of iterations
  % nit_max = 10000;  % max numer of iterations
  tdura = 0.4;
  npul_max = round((tsig(2)-tsig(1))/tdura);  % assume the max distinguishable number allows half overlapping
  nit_max;
  fpltit = 0;  % plot flag for each iteration
  fpltend = 1;  % plot flag for the final iteration
  fcheck = 0;  % plot flag for intermediate computations
  [sigdecon(:,ista),pred(:,ista),res(:,ista),dresit,mfitit,ampit,nit,fighdl] = ...
    iterdecon(sig,wlet,rcc,noi,dt,twlet,width,dres_min,mfit_min,nit_max,npul_max,fpltit,fpltend,fcheck);
%   [sigdecon(:,ista),pred(:,ista),res(:,ista),dresit,mfitit,ampit,nit,fighdl] = ...
%     iterdecon(sig,wlet,ones(size(sig)),noi,dt,twlet,width,dres_min,mfit_min,nit_max,npul_max,fpltit,fpltend,fcheck);
  ax = fighdl{2}.ax(1);
  hold(ax,'on');
  text(ax,0.05,0.85,stas(ista,:),'unit','normalized','FontSize',12)
  hold(ax,'off');
  
  ampitsta{ista} = ampit;
  nit
  
  % %check the spectral shape between wavelet, signal and prediction
  % [f] = plt_wletsigpred_freq(wlet,sig,pred(:,ista),sps);
  % ax = f.ax(1);
  % hold(ax,'on');
  % text(ax,0.8,0.85,stas(ista,:),'unit','normalized');
  % hold(ax,'off');
  
end

% %%
% figure
% [irccr,rccr12] = RunningCC(res(:,1), res(:,2), cclen);
% [~,rccr13] = RunningCC(res(:,1), res(:,3), cclen);
% [~,rccr23] = RunningCC(res(:,2), res(:,3), cclen);
% rccr = (rccr12+rccr13+rccr23)/3;
% p1=plot(ircc,rcc,'-','color',[.5 .5 .5],'LineWidth',1); hold on; grid on
% plot([0 350],[median(rcc) median(rcc)],'r--');
% p2=plot(irccr,rccr,'k-','LineWidth',1);
% legend([p1,p2],'Signal rcc','Residual rcc','Location','south');

%% direct CC between the deconvolved impulses
% figure;
% ax = gca;
% yyaxis(ax,'left');
% hold(ax,'on');
% ax.Box = 'on';
% imp1 = [];
% imp1(:,1) = find(sigdecon(:,1)>0);
% imp1(:,2) = sigdecon(sigdecon(:,1)>0, 1);
% stem(ax,imp1(:,1),imp1(:,2), 'r-','MarkerSize',5);
% imp2 = [];
% imp2(:,1) = find(sigdecon(:,2)>0);
% imp2(:,2) = sigdecon(sigdecon(:,2)>0, 2);
% stem(ax,imp2(:,1),imp2(:,2), 'b-','MarkerSize',5);
% imp3 = [];
% imp3(:,1) = find(sigdecon(:,3)>0);
% imp3(:,2) = sigdecon(sigdecon(:,3)>0, 3);
% stem(ax,imp3(:,1),imp3(:,2), 'k-','MarkerSize',5);
% 
% [coef12,lag12] = xcorr(sigdecon(:,1), sigdecon(:,2), sps/2, 'coeff');
% [mcoef12, idx] = max(coef12);   % max of master raw cc
% lagsamp12 = lag12(idx);   % offset in samples
% [coef13,lag13] = xcorr(sigdecon(:,1), sigdecon(:,3), sps/2, 'coeff');
% [mcoef13, idx] = max(coef13);   % max of master raw cc
% lagsamp13 = lag13(idx);   % offset in samples
% [coef23,lag23] = xcorr(sigdecon(:,2), sigdecon(:,3), sps/2, 'coeff');
% [mcoef23, idx] = max(coef23);   % max of master raw cc
% lagsamp23 = lag23(idx);   % offset in samples
% text(ax,0.1,0.9,sprintf('lag12: %d, max coef12: %.2f', lagsamp12,mcoef12),'Units','normalized',...
%   'HorizontalAlignment','left');
% text(ax,0.1,0.8,sprintf('lag13: %d, max coef13: %.2f', lagsamp13,mcoef13),'Units','normalized',...
%   'HorizontalAlignment','left');
% text(ax,0.1,0.7,sprintf('lag23: %d, max coef23: %.2f', lagsamp23,mcoef23),'Units','normalized',...
%   'HorizontalAlignment','left');
% 
% yyaxis(ax,'right');
% cctemp = rcc;
% cctemp(cctemp<0)=0;
% % plot(ax,ircc,cctemp,'co','markersize',4);
% plot(ax,ircc,rcc,'co','markersize',4);
% ylim(ax,[-1 1]);
% % xlim(ax,[50 300]);
% hold(ax,'off');

%% Group nearest impulses from different stations into pairs
%different stations has different number of non-zero impulses
for i = 1: nsta
  nimp(i) = sum(sigdecon(:,i)>0);
end
[npair, ista] = max(nimp);
imp = [];
imp(:,1) = find(sigdecon(:,ista)>0);
imp(:,2) = sigdecon(sigdecon(:,ista)>0, ista);
imp(:,3) = rcc(imp(:,1)-cclen/2);
imp(:,4) = imp(:,3).*imp(:,2);
%sort them based on the product of amplitude and rcc?
imp = sortrows(imp,4,'descend');

imppair = zeros(npair, 4*3);
for ip = 1: npair
  imppair(ip,(ista-1)*4+1:ista*4) = imp(ip,1:4);
end

spsscale = sps/40;
loff_max = 4*spsscale;
for i = 1:nsta
  if i == ista
    continue
  end
  imp1 = [];
  imp1(:,1) = find(sigdecon(:,i)>0);
  imp1(:,2) = sigdecon(sigdecon(:,i)>0, i);
  imp1(:,3) = rcc(imp1(:,1)-cclen/2);
  imp1(:,4) = imp1(:,3).*imp1(:,2);
  for ip = 1: npair
    [loff,itemp] = min(abs(imp1(:,1)-imp(ip,1))); % find the closest impulse
    if loff>loff_max
      continue
    end
    imppair(ip,(i-1)*4+1:i*4) = imp1(itemp,1:4); % pair them
    imp1(itemp,:) = [];   % mute the used ones
  end
end
indpair = find(sum(imppair==0,2)==0); %pairs that have corresponding peaks at all stations
imppairf = imppair(indpair,:);
impindep = imppairf(:, [1 2 5 6 9 10]);
impindep(:,7:9) = [impindep(:,1)-impindep(:,3) impindep(:,1)-impindep(:,5) ...
                   impindep(:,3)-impindep(:,5)];

figure;
colorn = {'r','b','k'};
for i = 1: nsta
  subplot(4,1,i)
  ax = gca;
  indtemp = find(sigdecon(:,i)>0);  % independent impulses from other stations
  stem(ax,indtemp,sigdecon(indtemp,i),'color',[.6 .6 .6],'MarkerSize',4); hold(ax,'on');
  text(ax,0.95,0.9,num2str(length(indtemp)),'unit','normalized','HorizontalAlignment','right',...
    'FontSize',10,'color',[.6 .6 .6]);
  indtemp = find(impindep(:,(i-1)*2+1)>0);  % that can be individually paired with 3-sta among independent
  stem(ax,impindep(indtemp,(i-1)*2+1),impindep(indtemp,(i-1)*2+2),'color',colorn{i},...
    'MarkerSize',4);
  text(ax,0.95,0.75,num2str(length(indtemp)),'unit','normalized','HorizontalAlignment','right',...
    'FontSize',10,'color',colorn{i});
  text(ax,0.05,0.9,strcat(stas(i,:)),'unit','normalized','HorizontalAlignment','left',...
    'FontSize',12);
  xlim(ax,[0 lsig]); ylim(ax,[0 2.5]); hold(ax,'off');
  ylabel(ax,'Amplitude','FontSize',12);
end

subplot(4,1,4)
ax = gca;
yyaxis(ax,'left');
hold(ax,'on');
ax.Box = 'on';
stem(ax,impindep(:,1),impindep(:,2), 'r-','MarkerSize',5);
stem(ax,impindep(:,3),impindep(:,4), 'b-','MarkerSize',5);
stem(ax,impindep(:,5),impindep(:,6), 'k-','MarkerSize',5);
text(ax,0.95,0.9,sprintf('Number of grouped triplets: %d', length(indpair)),'fontsize',10,...
  'Units','normalized','HorizontalAlignment','right');
xlabel(ax,sprintf('Samples at %d Hz',round(1/dt)),'FontSize',12);
ylabel(ax,'Amplitude','FontSize',12);
ylim(ax,[0 2.5]);
yyaxis(ax,'right');
plot(ircc,rcc,'o','color',[.5 .5 .5],'markersize',2);
ylim(ax,[-1 1]);
xlim(ax,[0 lsig]);
hold(ax,'off');

%%
%%%scatter of offsets in samples, and account for prealignment offset 
xran = [-loff_max loff_max]; 
yran = [-loff_max loff_max];
cran = [0 lsig];
impindept = impindep;
impindept(:,7:8) = impindept(:,7:8)+repmat([off1i(2) off1i(3)],size(impindept,1),1);
[f] = plt_decon_imp_scatter(impindept,xran,yran,cran,sps,1);
scatter(gca,off1i(2),off1i(3),30,'ks','filled','MarkerEdgeColor','k');
title(gca,sprintf('Independent, grouped, using data of %d Hz',sps));

%%% scatter of rela locations, and account for prealignment offset 
xran = [-2 2]; 
yran = [-2 2];
cran = [0 lsig];
impindept = impindep;
impindept(:,7:8) = impindept(:,7:8)+repmat([off1i(2) off1i(3)],size(impindept,1),1);
[f] = plt_decon_imp_scatter_space(impindept,xran,yran,cran,sps,'directhypo',0);
title(gca,sprintf('Independent, grouped, using data of %d Hz',sps));


%% Joint deconvoluion between 3 stations
spsscale = sps/40;

loff_max = 4*spsscale;
nit_max = 200;  % max numer of iterations 
fpltit = 0;  % plot flag for each iteration
fpltend = 1;  % plot flag for the final iteration
fpltchk = 0;  % plot flag for intermediate computations
[sigdecon,pred,res,dresit,mfitit,ampit,nit,fighdl] = ...
  iterdecon_joint(stasig,stawlet,stanoi,dt,statwlet,loff_max,dres_min,mfit_min,nit_max,npul_max,fpltit,fpltend,fpltchk);
for i = 1:nsta
  ax = fighdl{2}.ax(i+3);
  hold(ax,'on');
  text(ax,0.98,0.85,strtrim(stas(i,:)),'unit','normalized','HorizontalAlignment','right','FontSize',10)
  hold(ax,'off');
end

%%%'ampit' contains the indices of the impulses at each station and the offsets between each
%%%triplets. The sign of the offset is set so if off12 >0, it means sta 1 arrives later than sta 2,
%%%so you need to move 2 to the right --> to align with 1. This means, it has the SAME sign rule as
%%%the prealignment between the signals at 3 stations. Therefore, the net offset relative to the
%%%origin (0,0) is the sum of the 2 parts.

%%
%%%scatter of offsets in samples, and account for prealignment offset 
xran = [-loff_max loff_max]; 
yran = [-loff_max loff_max];
cran = [0 lsig];
ampitt = ampit;
ampitt(:,7:8) = ampitt(:,7:8)+repmat([off1i(2) off1i(3)],size(ampitt,1),1);
[f] = plt_decon_imp_scatter(ampitt,xran,yran,cran,sps,1);
scatter(gca,off1i(2),off1i(3),30,'ks','filled','MarkerEdgeColor','k');
title(gca,sprintf('Joint, using data of %d Hz',sps));

%%%scatter of rela locations, and account for prealignment offset 
xran = [-2 2]; 
yran = [-2 2];
cran = [0 lsig];
ampitt = ampit;
ampitt(:,7:8) = ampitt(:,7:8)+repmat([off1i(2) off1i(3)],size(ampitt,1),1);
[f] = plt_decon_imp_scatter_space(ampitt,xran,yran,cran,sps,'directhypo',0);
title(gca,sprintf('Joint, using data of %d Hz',sps));

%% Imagine you used some stopping criteria, and wonder if you have missed some important peaks
% %%%because you might have stopped too early
% %Let's check if the early iterations are EXACTLY the same if you carry out a bit more iterations
% %This needs at least two rounds and save the result of the 1st round
% % sz = min(size(ampit,1),size(ampit1,1));
% % isequaln(ampit(1:sz,:), ampit1(1:sz,:))   % this should be true
% %Let's check the trend of amplitude when you do a bit more iterations
% for i = 1: nsta
%   [coefnoi, lagnoi] = xcorr(stanoi(:,i), stawlet(:,i), size(stanoi, 1), 'none'); % unnormalized cc
%   coefnoieff = coefnoi(lagnoi+statwlet(i)/dt >= 1 & lagnoi+statwlet(i)/dt <= size(stanoi, 1));
%   [pkhgt, pkind] = findpeaks(coefnoieff);
%   ampnoi(i) = median(pkhgt);
% end
% ampnoin = ampnoi./sum(stawlet.^2, 1);   % normalized amp of CC between noise and wavelet 
% 
% figure
% subplot(311)
% plot(ampit(1:27,2),'linew',2); hold on
% plot(ampit(:,2),'linew',1); plot([0 nit],[ampnoin(1) ampnoin(1)],'--');
% plot([49 49],[0 2.5],'k--'); text(50, 2, 'min. AIC'); ylim([0 2.5]);
% 
% subplot(312)
% plot(ampit(1:27,4),'linew',2); hold on
% plot(ampit(:,4),'linew',1); plot([0 nit],[ampnoin(2) ampnoin(2)],'--');
% plot([49 49],[0 2.5],'k--'); ylim([0 2.5]);
% 
% subplot(313)
% plot(ampit(1:27,6),'linew',2); hold on
% plot(ampit(:,6),'linew',1); plot([0 nit],[ampnoin(3) ampnoin(3)],'--');
% plot([49 49],[0 2.5],'k--'); ylim([0 2.5]);
% xlabel('Iteration number');
% ylabel('Amplitude');
% 

%% compare independent results with joint results
figure;
subplot(211)
yyaxis left;
stem(impindep(:,1),impindep(:,2), 'r-','MarkerSize',5); hold on; box on
stem(impindep(:,3),impindep(:,4), 'b-','MarkerSize',5);
stem(impindep(:,5),impindep(:,6), 'k-','MarkerSize',5);
text(0.95,0.9,sprintf('Number of triplets: %d', size(impindep,1)),'fontsize',10,...
  'Units','normalized','HorizontalAlignment','right');
xlim([0 lsig]);
ylim([0 2.5]);
title(sprintf('Independent, grouped, using data of %d Hz',sps));
ylabel('Amplitude','FontSize',12);
yyaxis right;
plot(ircc,rcc,'o','color',[.5 .5 .5],'markersize',2);
ylim([-1 1]);

subplot(212)
stem(ampit(:,1),ampit(:,2), 'r-','MarkerSize',5); hold on; box on
stem(ampit(:,3),ampit(:,4), 'b-','MarkerSize',5);
stem(ampit(:,5),ampit(:,6), 'k-','MarkerSize',5);
text(0.95,0.9,sprintf('Number of triplets: %d', size(ampit,1)),'fontsize',10,...
  'Units','normalized','HorizontalAlignment','right');
title(sprintf('Joint, using data of %d Hz',sps));
xlim([0 lsig]); ylim([0 2.5]);
xlabel(sprintf('Samples at %d Hz',round(1/dt)),'FontSize',12);
ylabel('Amplitude','FontSize',12);

%%%scatter of offsets in samples, and account for prealignment offset 
xran = [-loff_max loff_max]; 
yran = [-loff_max loff_max];
cran = [0 lsig];
impindept = impindep;
impindept(:,7:8) = impindept(:,7:8)+repmat([off1i(2) off1i(3)],size(impindept,1),1);
[f] = plt_decon_imp_scatter(impindept,xran,yran,cran,sps,1);
title(gca,sprintf('Independent, grouped, using data of %d Hz',sps));

ampitt = ampit;
ampitt(:,7:8) = ampitt(:,7:8)+repmat([off1i(2) off1i(3)],size(ampitt,1),1);
[f] = plt_decon_imp_scatter(ampitt,xran,yran,cran,sps,1);
title(gca,sprintf('Joint, using data of %d Hz',sps));


%% Refinement, redo deconvlution starting from the first resolved impulse
%%%Looks like the result from refinement does not make it obviously better 
fpltit = 0;  % plot flag for each iteration
fpltend = 1;  % plot flag for the final iteration
fpltchk = 0;  % plot flag for intermediate computations
% impraw = median(ampit(:,[1 3 5]),2);
impraw = ampit;
[sigdecon,pred,res,dresit,mfitit,ampitrf,nit,fighdl] = ...
  iterdecon_joint_refine(stasig,stawlet,stanoi,impraw,dt,statwlet,loff_max,dres_min,mfit_min,nit_max,...
  npul_max,fpltit,fpltend,fpltchk);

ampitrft = ampitrf; ampitrft(:,[1 3 5 7 8 9]) = ampitrft(:,[1 3 5 7 8 9])/spsscale;
ampitrft = sortrows(ampitrft,1);
xran = [-max(abs(ampitrft(:,7))) max(abs(ampitrft(:,7)))];
yran = [-max(abs(ampitrft(:,8))) max(abs(ampitrft(:,8)))];
cran = [0 lsig/spsscale];
[f] = plt_decon_imp_scatter(ampitrft,xran,yran,cran,sps,spsscale);
reg = [-loff_max/spsscale -loff_max/spsscale;
       -loff_max/spsscale loff_max/spsscale;
       loff_max/spsscale loff_max/spsscale;
       loff_max/spsscale -loff_max/spsscale;
       -loff_max/spsscale -loff_max/spsscale;]; 
plot(reg(:,1),reg(:,2),'-','Color',[.5 .5 .5],'linew',1.5);
title(gca,sprintf('Joint, from 1st impulse in time, using data of %d Hz, normalized to 40 Hz',sps));


%% spectral division, Chao's codes
cctemp = [zeros(cclen/2,1); rcc; zeros(cclen/2,1)];
stdec1(:,ista) = specdiv_damp(sig,wlet,noi,dt,twlet);
figure
subplot(311)
stem(imp3(:,1),imp3(:,2), 'b-','MarkerSize',5); hold on
plot(10*stdec1(:,ista),'r');
plot([0 500],[10*rms(noi) 10*rms(noi)],'--','color',[.5 .5 .5]);
legend('Iterdecon','Specdiv\_damp');

stdec2(:,ista) = specdiv_water(sig,wlet,noi,dt,twlet);
subplot(312)
stem(imp3(:,1),imp3(:,2), 'b-','MarkerSize',5); hold on
plot(10*stdec2(:,ista),'r');
plot([0 500],[10*rms(noi) 10*rms(noi)],'--','color',[.5 .5 .5]);
legend('Iterdecon','Specdiv\_water');

stdec3(:,ista) = specdiv_arraycond(sig,wlet,noi,dt,twlet);
subplot(313)
stem(imp3(:,1),imp3(:,2), 'b-','MarkerSize',5); hold on
plot(10*stdec3(:,ista),'r');
plot([0 500],[10*rms(noi) 10*rms(noi)],'--','color',[.5 .5 .5]);
legend('Iterdecon','Specdiv\_arraycond');




