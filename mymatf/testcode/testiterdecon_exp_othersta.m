% testiterdecon_exp_othersta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to test the iterative deconvolution methods written by me, Chao Song,
% using some real data and real template, segments similar from 
% 'plt_tremor_spectra_ETS'.
%
% --The main difference between this version and 'testiterdecon_exp' is that
%   we want to look at the deconvolution result at the stations which have
%   similar envelope of seismogram to that of the trio stations in focus.
%   If the envelope is close, it might be useful to have more constraints
%   and/or validation to the 3-sta joint deconvolution and the paired
%   impulses from independently determined results.
% --In this script, we foucs on 100 Hz, as we already see the 
%   slight improvement from 80 Hz in 'testiterdecon_exp'.
% --May have to obtain new templates at current interested stations, if 
%   using 100 Hz, and at new stations. The rots parameters have already
%   been computed for all stations for many fams including 002, so no worry
%   about that.
% --After checking, the existed 100-sps templates are based on 'fixed' 
%   rots params from Yajun/Allan; Now, I added templates that are based
%   on rots params from 'old' and 'new' days and lfe catalog;
%   Existed 40-sps-60-s, and 80-sps-60-s templates have all 3 types, with 
%   suffixes indicating if using the 'old' or 'new' lfe catalog, but the
%   newly made 'fixed' ones have a subtle difference from the existed ones,
%   slightly more than 1-sample shift, and i can't tell why. At present, 
%   i am keeping the existed ones.
% --Among the 3 types of templates, ones based on 'new' lfe catalog is the
%   best, as it brings down the orthogonal components the most, this is most
%   obvious at station KLNB; ones based on 'old' is close to 'new'; 'fixed'
%   is the worst, as the orthogonal component at KLNB has comparable amplitude
%   to optimal components. SO, choose 'new' templates from now on!!!
% --At this point, we assume the time shift that are needed to prealign
%   the trace at an additional station is the same as the current 
%   interested stations. So you just align first based on the station 
%   location difference, then assume the source is at the prealignment of
%   3 stations and use their running CC.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/02/08
% Last modified date:   2022/02/08
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
FLAG = 'PGC';   % use pgc trio

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
      'LZB  '
      'TWKB '
      'MGCB '
      'KLNB '
%       'VGZ  ' % starting from VGZ, these stations are what we didn't pay much attention to
%       'GOWB '
%       'SNB  '
%       'TWBB '
      ];
nsta=size(stas,1);         %  number of stations

%% make template or read template
remake = 0;  % re-make the template or not, 1/0
dstack = [];
ccstack = [];
sps = 100;
templensec = 60;
% templensec = 32*4;
% CATA = 'fixed';  % use fixed rots params from Yajun/Allan   
CATA = 'new';  % use rots computed from new lfe catalog
% CATA = 'old';  % use rots computed from old lfe catalog
ccbp = [2 8];   % bandpass in CC the raw templates

if isequal(fam,'002') && isequal(CATA,'old')
  suffix = '_catold';
elseif isequal(fam,'002') && isequal(CATA,'new')
  suffix = '_catnew';
else
  suffix = [];
end

if remake   % if requested templates do not exist, recompute them
  ccmethod = 2; % only using lfes passed the cc thresholds (ccmethod=1); using all lfes (ccmethod=2)
  plflag = 0;
  if strcmp(FLAG, 'TWKB')
%     ccbp = [2 8];
    [dstack, ccstack] = mk_bbtemp_LZB(fam,sps,templensec,ccmethod,ccbp,plflag);
    ind = [5 4 6];
    stack = ccstack(ind, :);
    %write into files, NOTE the remade stacks contain all 7 stations
    allstas=['PGC  '
      'SSIB '
      'SILB '
      'LZB  '
      'TWKB '
      'MGCB '
      'KLNB '];
    allnsta = size(allstas,1);
    for ista = 1: allnsta
      fid = fopen(strcat(temppath, fam, '_', strtrim(allstas(ista, :)), '_', ...
        num2str(sps), 'sps_', num2str(templensec), 's_', ...
        'BBDS_', 'opt_Nof_Non_Chao'), 'w+');  % bandpassed direct stack, no filter, no norm
      
      fprintf(fid, '%f \n', dstack(ista, :)');
      fclose(fid);
      
      fid = fopen(strcat(temppath, fam, '_', strtrim(allstas(ista, :)), '_', ...
        num2str(sps), 'sps_', num2str(templensec), 's_', ...
        'BBCCS_', 'opt_Nof_Non_Chao'), 'w+');  % bandpassed cc stack, no filter, no norm
      fprintf(fid, '%f \n', ccstack(ista, :)');
      fclose(fid);
    end
    
  elseif strcmp(FLAG, 'PGC')
    if isequal(fam,'002')
%       CATA = 'fixed';
%       ccbp = [2 8];
      [dstack, ccstack] = mk_bbtemp_PGC(fam,CATA,sps,templensec,ccmethod,ccbp,plflag);
    else
%       ccbp = [2 8];
      [dstack, ccstack] = mk_bbtemp_PGC(fam,CATA,sps,templensec,ccmethod,ccbp,plflag);
    end
    ind = [1 2 3];
    stack = ccstack(ind, :);
    
    %write into files, NOTE the remade stacks contain all 7 stations
    allstas=['PGC  '
      'SSIB '
      'SILB '
      'LZB  '
      'TWKB '
      'MGCB '
      'KLNB '];
    nsta = size(allstas,1);
    
    %write into files
    for ista = 1: nsta
      fidds = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
        num2str(sps), 'sps_', num2str(templensec), 's_', ...
        'BBDS_', 'opt_Nof_Non_Chao',suffix), 'w+');  % bandpassed direct stack, no filter, no norm
      fidccs = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
        num2str(sps), 'sps_', num2str(templensec), 's_', ...
        'BBCCS_', 'opt_Nof_Non_Chao',suffix), 'w+');  % bandpassed cc stack, no filter, no norm
      fprintf(fidds, '%f \n', dstack(ista, :)');
      fclose(fidds);
      fprintf(fidccs, '%f \n', ccstack(ista, :)');
      fclose(fidccs);
    end
  end
  
else    % if exist, load them directly
  for ista = 1: nsta
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
      num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao',suffix);
    ccstack(ista,:) = load(fname);
  end
  tracetemp = ccstack';
end

figure
for ista = 1: nsta
  plot(tracetemp(:,ista)+0.3*ista); hold on
  text(size(tracetemp,1)*0.45, 0.3*ista, stas(ista,:));
  xlim([size(tracetemp,1)*0.4 size(tracetemp,1)*0.6] );
end

%% Double check the if existed and newly-made templates are the same
%load the template
% wlen = 60*sps;
% for ista = 1:3
%   fname = strcat(temppath, '/',fam,'_',stas(ista,:),'_',num2str(sps),'sps','_',num2str(wlen/sps),...
%     's_BBCCS_opt_Nof_Non_Chao_catold');
%   tracetemp(:,ista) = load(fname);
% end
% 
%  
% figure
% subplot(311)
% plot(tracetemp(:,1)); hold on;
% plot(stack(1,:));
% % xlim([2500 3500]);
% ylim([-0.3 0.3]);
% 
% subplot(312)
% plot(tracetemp(:,2)); hold on
% plot(stack(2,:));
% % xlim([2500 3500]);
% ylim([-0.3 0.3]);
% 
% subplot(313)
% plot(tracetemp(:,3)); hold on
% plot(stack(3,:));
% % xlim([2500 3500]);
% ylim([-0.3 0.3]);

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

if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
  stas(7,:)='KELB ';
else
  stas(7,:)='KLNB ';  % remember to change it back
end

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
[STAopt,STAort,~,fileflag] = rd_daily_bpdata(year,jday,datapath,stas,PERMSTA,POLSTA,...
  PERMROTS,POLROTS,sps,lo,hi,npo,npa,[],[],[],[]);

%read vertical components too, as we want to have a sense of the SNR at the orthogonal and vertical
%components
[STAvert,~,~,fileflag] = rd_daily_bpdata(year,jday,datapath,stas,PERMSTA,POLSTA,...
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


hfobj = hf(hf(:,seccol)>=tlst & hf(:,seccol)<=tled, :); % object hf tremor detections

%get a rough sense on 'noise' from other LFE/tremor from other regions to the target area 
bnduse = [min(hfobj(:,1)) min(hfobj(:,2));
          max(hfobj(:,1)) min(hfobj(:,2));
          max(hfobj(:,1)) max(hfobj(:,2));
          min(hfobj(:,1)) max(hfobj(:,2));
          min(hfobj(:,1)) min(hfobj(:,2));
          ];

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
  optcc = opttmp(buffer+1+isst-ilst: buffer+1+isst-ilst+ ised-isst, 1:3);
%   ortcc = orttmp(buffer+1+isst-ilst: end-buffer-(ised-iled), 1:3);
%   vertcc = verttmp(buffer+1+isst-ilst: end-buffer-(ised-iled), 1:3);
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
% 
for ista = 4: 7
  optref = opttmp(buffer+1+isst-ilst: buffer+1+isst-ilst+ ised-isst, 1);
  optcc = opttmp(buffer+1+isst-ilst: buffer+1+isst-ilst+ ised-isst, ista);
  [coef,lag] = xcorr(optref, optcc, mshiftadd, 'coeff');
  [mcoef, idx] = max(coef);   % max of master raw cc
  off1i(ista) = lag(idx);   % offset in samples
  off1i(ista) = 0;
  optdat(:, ista) = opttmp(buffer+1-off1i(ista): end-buffer-off1i(ista), ista);
  ortdat(:, ista) = orttmp(buffer+1-off1i(ista): end-buffer-off1i(ista), ista);
  vertdat(:, ista) = verttmp(buffer+1-off1i(ista): end-buffer-off1i(ista), ista);
end


%%% get the running average CC between 3 stations
cclen=sps/2;
[ircc,rcc12] = RunningCC(optdat(:,1), optdat(:,2), cclen);  % optimal comp.
[~,rcc13] = RunningCC(optdat(:,1), optdat(:,3), cclen);
[~,rcc23] = RunningCC(optdat(:,2), optdat(:,3), cclen);
rccopt = (rcc12+rcc13+rcc23)/3;
% 
% [ircc,rcc12] = RunningCC(ortdat(:,1), ortdat(:,2), cclen);  % orthogonal comp.
% [~,rcc13] = RunningCC(ortdat(:,1), ortdat(:,3), cclen);
% [~,rcc23] = RunningCC(ortdat(:,2), ortdat(:,3), cclen);
% rccort = (rcc12+rcc13+rcc23)/3;
% 
% [ircc,rcc12] = RunningCC(vertdat(:,1), vertdat(:,2), cclen);  % vertical comp.
% [~,rcc13] = RunningCC(vertdat(:,1), vertdat(:,3), cclen);
% [~,rcc23] = RunningCC(vertdat(:,2), vertdat(:,3), cclen);
% rccvert = (rcc12+rcc13+rcc23)/3;

%%% Optimal for all stations
f.fig = figure;
f.fig.Renderer = 'painters';
widin = 18;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);

ncol = 1;
nrow = 5;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
%     grid(f.ax(isub),'on');
end

pltxran = [0.03 0.95]; pltyran = [0.05 0.95];
pltxsep = 0.02; pltysep = 0.04; 
axpos = optaxposition(nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
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
yyaxis(ax,'left');
plot(ax,optdat(:,1),'r-');
plot(ax,optdat(:,2),'b-');
plot(ax,optdat(:,3),'k-'); %,'linew',0.5
text(ax,0.05,0.95,'Optimal','unit','normalized',...
  'HorizontalAlignment','left','FontSize',12);
xran = [0 tled-tlst]*sps;
yran(2) = 1.5*max(max(abs(optdat(:,1:3))));
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
noiran = (tnoi-tlst)*sps;
plot(ax,[noiran(1) noiran(1)],yran,'k--');
plot(ax,[noiran(2) noiran(2)],yran,'k--');
sigran = (tsig-tlst)*sps;
plot(ax,[sigran(1) sigran(1)],yran,'c--');
plot(ax,[sigran(2) sigran(2)],yran,'c--');
text(ax,0.95,0.95,sprintf('Bandpassed %.2f-%.2f Hz',losig,hisig),'unit','normalized',...
  'HorizontalAlignment','right','FontSize',12);
[envup, envlo] = envelope(optdat(:, 1:3));
medenvup = median(envup,2);
medenvlo = median(envlo,2);
plot(ax,medenvup,'-','color',[0.5 0.5 0.5],'linew',1.5);
plot(ax,medenvlo,'-','color',[0.5 0.5 0.5],'linew',1.5);
xlim(ax,xran);
xticks(ax,xran(1): 4*sps: xran(2))
xticklabels(ax,(xran(1)/sps+tlst): 4: (xran(2)/sps+tlst));
ylim(ax,yran);
xlabel(ax,'Time (s)','FontSize',10);
ylabel(ax,'Amplitude','FontSize',10);
yyaxis(ax,'right');
plot(ax,ircc,rccopt,'o','color',[.7 .7 .7],'markersize',2);
ylabel(ax,'Running CC','FontSize',10);
ylim(ax,[-1 1]);
hold(ax,'off');

for i = 2:5
  ax = f.ax(i);
  hold(ax,'on');
  plot(ax,optdat(:,i+2),'k-');
  text(ax,0.05,0.95,strcat({'Optimal '},stas(i+2,:)),'unit','normalized',...
    'HorizontalAlignment','left','FontSize',12);
  [envup, envlo] = envelope(optdat(:, i+2));
  plot(ax,envup,'-','color',[0.5 0.5 0.5],'linew',1.5);
  plot(ax,envlo,'-','color',[0.5 0.5 0.5],'linew',1.5);
  xlim(ax,xran);
  xticks(ax,xran(1): 4*sps: xran(2))
  xticklabels(ax,(xran(1)/sps+tlst): 4: (xran(2)/sps+tlst));
  ylim(ax,yran);
  hold(ax,'off');
end
% 
% %%% Orthogogal for all stations
% f.fig = figure;
% f.fig.Renderer = 'painters';
% widin = 18;  % maximum width allowed is 8.5 inches
% htin = 8;   % maximum height allowed is 11 inches
% set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);
% 
% ncol = 1;
% nrow = 5;
% for isub = 1:nrow*ncol
%     f.ax(isub) = subplot(nrow,ncol,isub);
%     f.ax(isub).Box = 'on';
% %     grid(f.ax(isub),'on');
% end
% 
% pltxran = [0.03 0.95]; pltyran = [0.05 0.95];
% pltxsep = 0.02; pltysep = 0.04; 
% axpos = optaxposition(nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
% for isub = 1:nrow*ncol
%   set(f.ax(isub), 'position', axpos(isub,:));
% end
% ax = f.ax(1);
% hold(ax,'on');
% yyaxis(ax,'left');
% plot(ax,ortdat(:,1),'r-');
% plot(ax,ortdat(:,2),'b-');
% plot(ax,ortdat(:,3),'k-'); %,'linew',0.5
% text(ax,0.05,0.95,'Orthogonal','unit','normalized',...
%   'HorizontalAlignment','left','FontSize',12);
% xlim(ax,xran);
% xticks(ax,xran(1): 4*sps: xran(2))
% xticklabels(ax,(xran(1)/sps+tlst): 4: (xran(2)/sps+tlst));
% ylim(ax,yran);
% 
% for i = 2:5
%   ax = f.ax(i);
%   hold(ax,'on');
%   plot(ax,ortdat(:,i+2),'k-');
%   text(ax,0.05,0.95,strcat({'Orthogonal '},stas(i+2,:)),'unit','normalized',...
%     'HorizontalAlignment','left','FontSize',12);
%   xlim(ax,xran);
%   xticks(ax,xran(1): 4*sps: xran(2))
%   xticklabels(ax,(xran(1)/sps+tlst): 4: (xran(2)/sps+tlst));
%   ylim(ax,yran);
%   hold(ax,'off');
% end
%
% 
% %%% Vertical for all stations
% f.fig = figure;
% f.fig.Renderer = 'painters';
% widin = 18;  % maximum width allowed is 8.5 inches
% htin = 8;   % maximum height allowed is 11 inches
% set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);
% 
% ncol = 1;
% nrow = 5;
% for isub = 1:nrow*ncol
%     f.ax(isub) = subplot(nrow,ncol,isub);
%     f.ax(isub).Box = 'on';
% %     grid(f.ax(isub),'on');
% end
% 
% pltxran = [0.03 0.95]; pltyran = [0.05 0.95];
% pltxsep = 0.02; pltysep = 0.04; 
% axpos = optaxposition(nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
% for isub = 1:nrow*ncol
%   set(f.ax(isub), 'position', axpos(isub,:));
% end
% ax = f.ax(1);
% hold(ax,'on');
% yyaxis(ax,'left');
% plot(ax,vertdat(:,1),'r-');
% plot(ax,vertdat(:,2),'b-');
% plot(ax,vertdat(:,3),'k-'); %,'linew',0.5
% text(ax,0.05,0.95,'Vertical','unit','normalized',...
%   'HorizontalAlignment','left','FontSize',12);
% xlim(ax,xran);
% xticks(ax,xran(1): 4*sps: xran(2))
% xticklabels(ax,(xran(1)/sps+tlst): 4: (xran(2)/sps+tlst));
% ylim(ax,yran);
% 
% for i = 2:5
%   ax = f.ax(i);
%   hold(ax,'on');
%   plot(ax,vertdat(:,i+2),'k-');
%   text(ax,0.05,0.95,strcat({'Vertical '},stas(i+2,:)),'unit','normalized',...
%     'HorizontalAlignment','left','FontSize',12);
%   xlim(ax,xran);
%   xticks(ax,xran(1): 4*sps: xran(2))
%   xticklabels(ax,(xran(1)/sps+tlst): 4: (xran(2)/sps+tlst));
%   ylim(ax,yran);
%   hold(ax,'off');
% end

%% pre-process the signal and template for a similar spectral shape essential for deconvolution
stasig = [];
stawlet = [];
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
[off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(tmpwlet(:,1:3)',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
offwlet1i = zeros(nsta,1);
offwlet1i(2) = round(off12con);
offwlet1i(3) = round(off13con);

for ista = 4: nsta
  [coef,lag] = xcorr(tmpwlet(:,1), tmpwlet(:,ista), mshiftadd, 'coeff');
  [mcoef, idx] = max(coef);   % max of master raw cc
  offwlet1i(ista) = lag(idx);   % offset in samples  
end

figure
for ista = 1:nsta
  plot(tmpwlet(:,ista)+0.3*ista); hold on
  text(size(tmpwlet,1)*0.1, 0.3*ista, stas(ista,:));
end

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

%Plot spectra to see if they are similar in shape
for ista = 1: nsta
  
  %according analysis to the template, the zero-crossing to the end of coda is safe to set as 8 s
  wletbb(:,ista) = tracetemp(zcsta1+8*sps-lwlet+1-offwlet1i(ista): ...
                             zcsta1+8*sps-offwlet1i(ista), ista);
%   wletbb(:,ista) = tracetemp(zerocs(ista)+8*sps-lwlet+1: zerocs(ista)+8*sps, ista);
  
%   %%% remove mean and linear trend and taper
%   %of signal
%   sigbbd = detrend(sigbb);
%   %and taper with tukeywin, which is actually a tapered cosine window
%   %tapered length is adaptative to frequency, maybe at least longer than one full period length of
%   %the lowest frequency
%   fractap = sps*4/size(sigbbd,1); % if fractap is >=1, n-point von Hann window is returned
%   ptstap = fractap/2*size(sigbbd,1); % points tapered at the start/end 
%   w = tukeywin(size(sigbbd,1),fractap);
% %   sigbbdt = w.* sigbbd;
%   sigbbdt = ones(size(sigbbd)).* sigbbd;  % don't taper
%   %detrend again for caution
%   sigbbdt=detrend(sigbbdt);

%   %romve mean, linear trend of template
%   wletbbd = detrend(wletbb(:,ista));
%   %and taper with tukeywin, which is actually a tapered cosine window
%   w = tukeywin(size(wletbbd,1),fractap);
% %   wletbbdt = w.* wletbbd;
%   wletbbdt = ones(size(wletbbd)).* wletbbd;  % don't taper
%   %detrend again for caution
%   wletbbdt=detrend(wletbbdt);
  
  %%% filter the signal and template to reach a similar spectra shape
  %%% 1.8-4.5 for signal and 1.8-6.0 for template seem like a good pair
  %%% OR, 1.8-5.4 for signal and 1.8-18 for template to keep the most of the template at the high end
  %%% OR, use 1.8-4.5 for both, but use more poles (eg., 3) for signal to achieve a faster decay
%   %filter the signal
%   hisig=6.3;    % seems like signal has higher amp, so 6.3 leads to a closer decay rate
%   losig=1.8;
%   sig = Bandpass(sigbbdt, sps, losig, hisig, npo, npa, 'butter');
%   %detrend again for caution
%   sig=detrend(sig);

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

% if sps==40
%   statwlet = [191 191 191]/sps;  % Time shift of main arrival of wavelet
% elseif sps==80
%   statwlet = [384 384 384]/sps;  % Time shift of main arrival of wavelet
% elseif sps==100
%   statwlet = [494 494 494]/sps;  % Time shift of main arrival of wavelet
% end


%% Before doing all the rest, we need to ensure the seismogram has a 'similar' shape
%%%1. Similarity in shape could be reflected by the running CC between station i and the reference
cclen=sps/2;
rcc1i(:,1) = ones(lsig-cclen,1);
[ircc,rcc1i(:,2)] = RunningCC(stasig(:,1), stasig(:,2), cclen);  % optimal comp.
[~,rcc1i(:,3)] = RunningCC(stasig(:,1), stasig(:,3), cclen);
[~,rcc23] = RunningCC(stasig(:,2), stasig(:,3), cclen);
rcc = (rcc1i(:,2)+rcc1i(:,3)+rcc23)/3;
for ista = 4:nsta
  [~,rcc1i(:,ista)] = RunningCC(stasig(:,1), stasig(:,ista), cclen);
end

%%%2. Similarity in shape could be also reflected by the envelope CC between station i and the 
%%%reference
[envup,~] = envelope(stasig);
medenvup = median(envup(:,1:3),2);

ftobj = fittype('gauss1');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Robust = 'Bisquare';
opts.StartPoint = [0.87498036315248 -16 186.736840939868];

figure
mu = [];
for ista = 4: nsta
  [coef,lag] = xcorr(medenvup, envup(:,ista), length(medenvup), 'coeff');
  [mcoef, idx] = max(coef);   % max of master raw cc
  off1i(ista) = lag(idx);   % offset in samples
  subplot(nsta-3,1,ista-3);
  plot(lag,coef); hold on
  scatter(lag(idx),mcoef,20,'k','filled');
  % Fit model to data.
  [xData, yData] = prepareCurveData(lag',coef);
  [fitrst, gof,~] = fit(xData,yData,ftobj,opts);
  param = coeffvalues(fitrst);
  mu(ista) = param(2);
  sigma = param(3);
  text(-800,0.8,num2str(lag(idx)),'Color','k');
  text(-800,0.6,num2str(mu(ista)),'Color','b');
end

%% plot the pre-processed signal in more detail
f.fig = figure;
f.fig.Renderer = 'painters';
widin = 12;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);

ncol = 1;
nrow = 5;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
%     grid(f.ax(isub),'on');
end

pltxran = [0.05 0.95]; pltyran = [0.15 0.95];
pltxsep = 0.02; pltysep = 0.04; 
axpos = optaxposition(nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
for isub = 1:nrow*ncol
  set(f.ax(isub), 'position', axpos(isub,:));
end

timemax = (hfobj(:, seccol) -tsig(1))*sps +1; % starting time of max power rate of half sec
timecnt = (hfobj(:, seccol-1) -tsig(1))*sps +1; % central time of the 4-s detecting window
off12all = round(hfobj(:, 11)*sps); % offsets of 4-s tremor detections
off13all = round(hfobj(:, 12)*sps);

timearm = round((armobj(:,4) -tsig(1))*sps +1); % time of armbruster's tremor detections
timebost = round((bostobj(:,5) -tsig(1))*sps +1); % time of armbruster's tremor detections

%   isub = (ista-1)*2+1;
ax = f.ax(1);
hold(ax,'on');
yyaxis(ax,'left');
plot(ax,stasig(:,1),'r-');
plot(ax,stasig(:,2),'b-');
plot(ax,stasig(:,3),'k-'); %,'linew',0.5
text(ax,0.05,0.95,'Optimal','unit','normalized',...
  'HorizontalAlignment','left','FontSize',12);
xran = [0 tsig(2)-tsig(1)]*sps;
yran(2) = 1.5*max(max(abs(stasig(:,1:3))));
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
text(ax,0.95,0.95,sprintf('Bandpassed %.2f-%.2f Hz',losig,hisig),'unit','normalized',...
  'HorizontalAlignment','right','FontSize',12);
% medenvup = median(envup(:,1:3),2);
% medenvlo = median(envlo,2);
plot(ax,medenvup,'-','color',[0.5 0.5 0.5],'linew',1.5);
% plot(ax,medenvlo,'-','color',[0.5 0.5 0.5],'linew',1.5);
xlim(ax,xran);
ylim(ax,yran);
xlabel(ax,'Samples','FontSize',10);
ylabel(ax,'Amplitude','FontSize',10);

yyaxis(ax,'right');
plot(ax,ircc,rcc,'o','color',[.7 .7 .7],'markersize',2);
plot(ax,[0 lsig],[0 0],'--','color',[0.8 0.8 0.8]);
ylabel(ax,'Running CC','FontSize',10);
ylim(ax,[-1 1]);
hold(ax,'off');

for i = 2: 5
  ax = f.ax(i);
  hold(ax,'on');
  plot(ax,stasig(:,i+2),'k-');
  text(ax,0.05,0.95,strcat({'Optimal '},stas(i+2,:)),'unit','normalized',...
    'HorizontalAlignment','left','FontSize',12);
  text(ax,0.25,0.95,num2str(off1i(i+2)),'unit','normalized',...
    'HorizontalAlignment','left','FontSize',10);
  plot(ax,envup(:,i+2),'-','color',[0.5 0.5 0.5],'linew',1.5);
%   plot(ax,envlo,'-','color',[0.5 0.5 0.5],'linew',1.5);
  xlim(ax,xran);
  % ylim(ax,yran);
  
  yyaxis(ax,'right');
  plot(ax,ircc,rcc1i(:,i+2),'o','color',[.7 .7 .7],'markersize',2);
  plot(ax,[0 lsig],[0 0],'--','color',[0.8 0.8 0.8]);
  ylabel(ax,'Running CC','FontSize',10);
  ylim(ax,[-1 1]);
  hold(ax,'off');
end


%% independent iterative deconvolution
for ista = 1: nsta
%   ista = 4 % which station are we focusing on?
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
  nit_max = 50;  % max numer of iterations 
  % nit_max = 10000;  % max numer of iterations 
  tdura = 0.4;
  npul_max = round((tsig(2)-tsig(1))/tdura);  % assume the max distinguishable number allows half overlapping
  nit_max
  fpltit = 0;  % plot flag for each iteration
  fpltend = 1;  % plot flag for the final iteration
  fcheck = 0;  % plot flag for intermediate computations
  [sigdecon(:,ista),pred(:,ista),res(:,ista),dresit,mfitit,ampit,nit,fighdl] = ...
    iterdecon(sig,wlet,rcc,noi,dt,twlet,width,dres_min,mfit_min,nit_max,npul_max,fpltit,fpltend,fcheck);
  % [sigdecon(:,ista),pred,res,dresit,mfitit,ampit,nit,fighdl] = ...
  %   iterdecon(sig,wlet,ones(size(sig)),noi,dt,twlet,width,dres_min,mfit_min,nit_max,npul_max,fpltit,fpltend);
  ax = fighdl{2}.ax(1);
  hold(ax,'on');
  text(ax,0.8,0.85,stas(ista,:),'unit','normalized');
  hold(ax,'off');

  ampitsta{ista} = ampit;

end

%% Group nearest impulses from the 3 controling stations into pairs
%different stations has different number of non-zero impulses
for i = 1: 3
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
for i = 1:3
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
  stem(ax,indtemp,sigdecon(indtemp,i),'color',[.8 .8 .8],'MarkerSize',4); hold(ax,'on');
  text(ax,0.9,0.9,num2str(length(indtemp)),'unit','normalized','HorizontalAlignment','right',...
    'FontSize',10,'color',[.8 .8 .8]);
  indtemp = find(impindep(:,(i-1)*2+1)>0);  % that can be individually paired with 3-sta among independent
  stem(ax,impindep(indtemp,(i-1)*2+1),impindep(indtemp,(i-1)*2+2),'color',colorn{i},...
    'MarkerSize',4);
  text(ax,0.9,0.75,num2str(length(indtemp)),'unit','normalized','HorizontalAlignment','right',...
    'FontSize',10,'color',colorn{i});
  text(ax,0.05,0.9,strcat(stas(i,:)),'unit','normalized','HorizontalAlignment','left',...
    'FontSize',12);
  xlim(ax,[0 lsig]); hold(ax,'off');
end

subplot(4,1,4)
ax = gca;
yyaxis(ax,'left');
hold(ax,'on');
ax.Box = 'on';
stem(ax,impindep(:,1),impindep(:,2), 'r-','MarkerSize',5);
stem(ax,impindep(:,3),impindep(:,4), 'b-','MarkerSize',5);
stem(ax,impindep(:,5),impindep(:,6), 'k-','MarkerSize',5);
text(ax,0.1,0.9,sprintf('Number of triplets: %d', length(indpair)),'fontsize',10,...
  'Units','normalized');
yyaxis(ax,'right');
plot(ircc,rcc,'o','color',[.5 .5 .5],'markersize',2);
ylim(ax,[-1 1]);
% xlim(ax,[50 300]);
hold(ax,'off');

impindept = impindep; 
impindept(:,7:8) = impindept(:,7:8)+repmat([off1i(2) off1i(3)],size(impindept,1),1);
impindept(:,[1 3 5 7 8 9]) = impindept(:,[1 3 5 7 8 9])/spsscale;
xran = [-loff_max/spsscale loff_max/spsscale]; 
yran = [-loff_max/spsscale loff_max/spsscale];
cran = [0 lsig/spsscale];
[f] = plt_decon_imp_scatter(impindept,xran,yran,cran,sps,spsscale);
title(gca,sprintf('Independent, grouped, using data of %d Hz, normalized to 40 Hz',sps));


%% identify the impulses that are close enough to the paired triplets 
imppairf(:,13) = median(imppairf(:,[4 8 12]), 2);  % median of the weighted coefs
imppairf = sortrows(imppairf,13,'descend');  % use this median to rank the triplets to decide which to pair first
imppairelse = zeros(size(imppairf,1), 4*4);
for i = 4:7
  imp1 = [];
  imp1(:,1) = find(sigdecon(:,i)>0);
  imp1(:,2) = sigdecon(sigdecon(:,i)>0, i);
  imp1(:,3) = rcc(imp1(:,1)-cclen/2);
  imp1(:,4) = imp1(:,3).*imp1(:,2);
  for ip = 1: size(imppairf,1)
    [loff,itemp] = min(abs(imp1(:,1)-median(imppairf(ip,[1 5 9])))); % find the closest impulse
    if loff>loff_max
      continue
    end
    imppairelse(ip,(i-4)*4+1:(i-3)*4) = imp1(itemp,1:4); % pair them
    imp1(itemp,:) = [];   % mute the used ones
  end

end
% indpairelse = find(sum(imppairelse==0,2)==0); %pairs that have corresponding peaks at all stations
% imppairelsef = imppairelse(indpairelse,:);
% impindep = imppairelsef(:, [1 2 5 6 9 10]);
% impindep(:,7:9) = [impindep(:,1)-impindep(:,3) impindep(:,1)-impindep(:,5) ...
%                    impindep(:,3)-impindep(:,5)];

%% make a summary plot of other stations
figure;
subplot(5,1,1)
ax = gca;
yyaxis(ax,'left');
hold(ax,'on');
ax.Box = 'on';
stem(ax,imppairf(:,1),imppairf(:,2), 'r-','MarkerSize',5);  % paired triplets from 3 stations
stem(ax,imppairf(:,5),imppairf(:,6), 'b-','MarkerSize',5);
stem(ax,imppairf(:,9),imppairf(:,10), 'k-','MarkerSize',5);
text(ax,0.1,0.9,sprintf('Number of triplets: %d', length(indpair)),'fontsize',10,...
  'Units','normalized');
yyaxis(ax,'right');
plot(ircc,rcc,'o','color',[.5 .5 .5],'markersize',2);
plot([0 lsig],[0 0],'--','color',[.3 .3 .3]);
ylim(ax,[-1 1]);
xlim(ax,[0 lsig]);
hold(ax,'off');

for i = 2:5
  subplot(5,1,i)
  ax = gca;
  indtemp = find(sigdecon(:,i+2)>0);  % independent impulses from other stations
  stem(ax,indtemp,sigdecon(indtemp,i+2),'color',[.8 .8 .8],'MarkerSize',4); hold on
  text(ax,0.9,0.9,num2str(length(indtemp)),'unit','normalized','HorizontalAlignment','right',...
    'FontSize',10,'color',[.8 .8 .8]);
  indtemp = find(imppairelse(:,(i-2)*4+1)>0);  % that can be individually paired with 3-sta among independent
  stem(ax,imppairelse(indtemp,(i-2)*4+1),imppairelse(indtemp,(i-2)*4+2),'color',[.2 .2 .2],...
    'MarkerSize',4);
  text(ax,0.9,0.75,num2str(length(indtemp)),'unit','normalized','HorizontalAlignment','right',...
    'FontSize',10,'color',[.2 .2 .2]);
  indtemp = find(sum(imppairelse==0,2)==0); % that can be  paired with 3-sta among ALL other stations
  stem(ax,imppairelse(indtemp,(i-2)*4+1),imppairelse(indtemp,(i-2)*4+2),'color','r','MarkerSize',4);
  text(ax,0.9,0.6,num2str(length(indtemp)),'unit','normalized','HorizontalAlignment','right',...
    'FontSize',10,'color','r');
  text(ax,0.05,0.9,strcat(stas(i+2,:)),'unit','normalized','HorizontalAlignment','left',...
    'FontSize',12);
  xlim(ax,[0 lsig]);
end


%% Joint deconvoluion between 3 stations
spsscale = sps/40;

loff_max = 4*spsscale;
nit_max = 200;  % max numer of iterations 
fpltit = 0;  % plot flag for each iteration
fpltend = 1;  % plot flag for the final iteration
fpltchk = 0;  % plot flag for intermediate computations
[sigdecon,pred,res,dresit,mfitit,ampit,nit,fighdl] = ...
  iterdecon_joint(stasig(:,1:3),stawlet(:,1:3),stanoi(:,1:3),dt,statwlet(1:3),loff_max,...
                  dres_min,mfit_min,nit_max,npul_max,fpltit,fpltend,fpltchk);

%%%'ampit' contains the indices of the impulses at each station and the offsets between each
%%%triplets. The sign of the offset is set so if off12 >0, it means sta 1 arrives later than sta 2,
%%%so you need to move 2 to the right --> to align with 1. This means, it has the SAME sign rule as
%%%the prealignment between the signals at 3 stations. Therefore, the net offset relative to the
%%%origin (0,0) is the sum of the 2 parts.
ampitt = ampit; 
ampitt(:,7:8) = ampitt(:,7:8)+repmat([off1i(2) off1i(3)],size(ampitt,1),1);
ampitt(:,[1 3 5 7 8 9]) = ampitt(:,[1 3 5 7 8 9])/spsscale;
ampitt = sortrows(ampitt,1);
xran = [-loff_max/spsscale loff_max/spsscale]; 
yran = [-loff_max/spsscale loff_max/spsscale];
cran = [0 lsig/spsscale];
[f] = plt_decon_imp_scatter(ampitt,xran,yran,cran,sps,spsscale);
title(gca,sprintf('Using data of %d Hz, normalized to 40 Hz',sps));

xran = [-3 3]; 
yran = [-3 3];
cran = [0 lsig];
[f] = plt_decon_imp_scatter_space(ampit,xran,yran,cran,sps,'interpArmb',0);
title(gca,sprintf('Joint, using data of %d Hz',sps));

ampitt = ampit; 
ampitt(:,7:8) = ampitt(:,7:8)+repmat([off1i(2) off1i(3)],size(ampitt,1),1);
[f] = plt_decon_imp_scatter_space(ampitt,xran,yran,cran,sps,'interpArmb',0);
title(gca,sprintf('Joint, using data of %d Hz',sps));

%% compare independent results with joint results
figure;
subplot(211)
yyaxis left;
stem(impindep(:,1),impindep(:,2), 'r-','MarkerSize',5); hold on; box on
stem(impindep(:,3),impindep(:,4), 'b-','MarkerSize',5);
stem(impindep(:,5),impindep(:,6), 'k-','MarkerSize',5);
text(0.1,0.9,sprintf('Number of triplets: %d', size(impindep,1)),'fontsize',10,...
  'Units','normalized');
text(0.9,0.9,sprintf('Independent, grouped'),'fontsize',10,...
  'Units','normalized','HorizontalAlignment','right');
xlim([0 lsig]);
ylim([0 2.5]);
yyaxis right;
plot(ircc,rcc,'o','color',[.5 .5 .5],'markersize',2);
ylim([-1 1]);

subplot(212)
stem(ampit(:,1),ampit(:,2), 'r-','MarkerSize',5); hold on; box on
stem(ampit(:,3),ampit(:,4), 'b-','MarkerSize',5);
stem(ampit(:,5),ampit(:,6), 'k-','MarkerSize',5);
text(0.1,0.9,sprintf('Number of triplets: %d', size(ampit,1)),'fontsize',10,...
  'Units','normalized');
text(0.9,0.9,sprintf('Joint'),'fontsize',10,...
  'Units','normalized','HorizontalAlignment','right');
xlim([0 lsig]); ylim([0 2.5]);


xran = [-loff_max loff_max]; 
yran = [-loff_max loff_max];
cran = [0 lsig];
[f] = plt_decon_imp_scatter(impindep,xran,yran,cran,sps,1);
title(gca,sprintf('Independent, grouped, using data of %d Hz',sps));

[f] = plt_decon_imp_scatter(ampit,xran,yran,cran,sps,1);
title(gca,sprintf('Joint, using data of %d Hz',sps));














