% cc_synth_onespot.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --This script aims to do some analysis to the synthetic seismograms other
% than deconvolution that is done particularly in 'decon_synth.m', to be
% compared to real data.
%
% Analysis would include:
% 1. tracking RCC variation with different saturation level and region size;
% 2. tracking RCC vs. envelope (?)
% --This script itself should be similar to 'synthrccvssatur.m' and
% 'synthrccvssrcregsize.m', which were designed to analyze synthetics
% generated by older versions of codes.
% --some part of comparing the RCC w/ the real data is similar to that in
% 'rccexperiment.m'
%
% NOTE: use this code rather than 'rccexperiment', 'synthrccvssrcregsize',
% 'synthrccvssatur', IF the source distri is uniform in space. If it is custom
% PDF like 2D Gaussian or from 4-s tremors, then go to these individual codes!
%
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/05/26
% Last modified date:   2023/09/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e   % Set the format to 5-digit floating point
clear
clc
close all

%% for easy testing
defval('normflag',0); %whether to normalize templates
defval('rccmwsec',0.5); %moving win len in sec for computing RCC

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

%% prepare templates (Green's functions), from 'lfetemp002_160sps.m' or Allan's templates
% tempflag = 'allan';
tempflag = 'chao';

adatapath = '/home/data2/chaosong/matlab/allan/matfils/';  %path for Allan's data

stas=['PGC  '
  'SSIB '
  'SILB '
  %   'LZB  '
  %   'TWKB '
  %   'MGCB '
  'KLNB '
  ]; % determine the trio and order, here the 1st sta is PGC
nsta=size(stas,1);         %  number of stations

if strcmp(tempflag,'chao')
  fam = '002';   % family number
  sps = 160;
  templensec = 60;
  ccstack = [];
  for ista = 1: nsta
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
      num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao_catnew');
    ccstack(:,ista) = load(fname);
  end
  STA = detrend(ccstack);
elseif strcmp(tempflag,'allan')
  sps = 100;
  templensec = 120;
  fname = ['PGCopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh_0.5-15Hz_CC1.25-6.5Hz_le20shift';
    'SSIopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh_0.5-15Hz_CC1.25-6.5Hz_le20shift';
    'SILopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh_0.5-15Hz_CC1.25-6.5Hz_le20shift'];
  %   templensec = 60;
  %   fname = ['PGCopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh0.5-15Hz_le14sh';
  %            'SSIopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh0.5-15Hz_le14sh';
  %            'SILopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh0.5-15Hz_le14sh'];
  for ista=1:nsta
    temp=load(strcat(adatapath,fname(ista,:)));
    STA(:,ista)=detrend(temp(:,1))/350;
  end
end

% %plot the raw templates, not filtered, not best aligned
% figure
% ltemp = size(STA,1);
% subplot(111)
% hold on
% for ista = 1: nsta
%   plot((1:ltemp)/sps,STA(:,ista));
% end

%%%The below aligns the templates by x-correlation
if strcmp(tempflag,'chao')
  ist = templensec*sps*4/10;  %not using whole window in case any station has very long-period energy
  ied = templensec*sps*6/10;
elseif strcmp(tempflag,'allan')
  ist = templensec*sps*6/10;
  ied = templensec*sps*8/10;
end
[maxses,imaxses]=max(STA(ist:ied,:),[],1);
[minses,iminses]=min(STA(ist:ied,:),[],1);
spread=maxses-minses;
imaxses = imaxses+ist-1;  %convert to global indices
iminses = iminses+ist-1;
% zcrosses=round(0.5*(imaxses+iminses));  % rough, assuming symmetry, Chao 2021/07/16
%automatically find the zero-crossings
zcrosses = zeros(nsta,1);
for ista = 1:nsta
  seg = detrend(STA(iminses(ista): imaxses(ista),ista));  % for zero-crossing timing, only use the main station
  [~,zcrosses(ista)] = min(abs(seg));
  zcrosses(ista) = zcrosses(ista)-1+iminses(ista);  % convert to global index
end
%now we want to cut a segment around the zero-crossing at each station
sampbef=6*sps;
sampaft=10*sps;
is=zcrosses-sampbef;
ie=zcrosses+sampaft;
for ista=1:nsta
  STAtmp(:,ista)=detrend(STA(is(ista):ie(ista),ista));  % this means templates are 'aligned' at zero-crossings
end
%x-correlation independently between each station pair
mshiftadd=10*sps/40;
tempxc(:,1)=xcorr(STAtmp(:,2),STAtmp(:,3),mshiftadd,'coeff');
tempxc(:,2)=xcorr(STAtmp(:,1),STAtmp(:,2),mshiftadd,'coeff'); %shift STAtmp(3,:) to right for positive values
tempxc(:,3)=xcorr(STAtmp(:,1),STAtmp(:,3),mshiftadd,'coeff'); %shift STAtmp(2,:) to right for positive values
for ista = 4: nsta
  tempxc(:,ista)=xcorr(STAtmp(:,1),STAtmp(:,ista),mshiftadd,'coeff'); %shift STAtmp(2,:) to right for positive values
end
[~,imax]=max(tempxc,[],1);
imax=imax-(mshiftadd+1); %This would produce a slightly different shift, if filtered seisms were used.
imax(2)-imax(3)+imax(1);   %enclosed if it equals to 0
for ista=2:nsta
  STAtmp(mshiftadd+1:end-(mshiftadd+1),ista)=detrend(STAtmp(mshiftadd+1-imax(ista):end-(mshiftadd+1)-imax(ista),ista));
end
%normalization
if normflag
  for ista=1:nsta
    STAtmp(:,ista)=STAtmp(:,ista)/spread(ista); % now templates are 'aligned' indeoendently by x-corr wrt. sta 1
  end
end
% figure
% ltemp = size(STAtmp,1);
% subplot(111)
% hold on
% for ista = 1: nsta
%   plot((1:ltemp)/sps,STAtmp(:,ista));
% end
%%%The above aligns the templates by x-correlation

%%%detrend, taper and bandpass templates
tmpwlet = STAtmp; % no bandpass
tmpwletf = STAtmp;  % bandpassed version
fractap = sps/size(tmpwlet,1);
for ista = 1: nsta
  %romve mean, linear trend of template
  tmpwlet(:,ista) = detrend(tmpwlet(:,ista));
  %and taper with tukeywin, which is actually a tapered cosine window
  w = tukeywin(size(tmpwlet(:,ista),1),fractap);
  tmpwlet(:,ista) = w.* tmpwlet(:,ista);
  %detrend again for caution
  tmpwlet(:,ista) = detrend(tmpwlet(:,ista));
  %filter the template
  hiwlet=18;
  lowlet=1.8;
  tmpwletf(:,ista) = Bandpass(tmpwlet(:,ista), sps, lowlet, hiwlet, 2, 2, 'butter');
  %detrend again for caution
  tmpwletf(:,ista) = detrend(tmpwletf(:,ista));
end

%%%constrained CC, so that only 2 offsets are independent
ccmid = round(size(tmpwletf,1)/2);
ccwlen = 10*sps;
loffmax = 5*sps/40;
ccmin = 0.01;  % depending on the length of trace, cc could be very low
iup = 1;    % times of upsampling
[off12conf,off13conf,ccf] = constrained_cc_interp(tmpwletf(:,1:3)',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
offwlet1if(1) = 0;
offwlet1if(2) = round(off12conf);
offwlet1if(3) = round(off13conf);
if nsta>3
  for ista = 4: nsta
    [mcoef,offwlet1if(ista)] = xcorrmax(tmpwletf(:,1), tmpwletf(:,ista), mshiftadd, 'coeff');
  end
end
%%%for the broadband templates as well
[off12con,off13con,cc] = constrained_cc_interp(tmpwlet(:,1:3)',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
offwlet1i(1) = 0;
offwlet1i(2) = round(off12con);
offwlet1i(3) = round(off13con);
if nsta>3
  for ista = 4: nsta
    [mcoef,offwlet1i(ista)] = xcorrmax(tmpwlet(:,1), tmpwlet(:,ista), mshiftadd, 'coeff');
  end
end

%%%automatically find the rough zero-crossing time, whose abs. value is closest to 0, whether + or -
[~,imin] = min(tmpwletf(:,1));
[~,imax] = max(tmpwletf(:,1));
[~,zcsta1] = min(abs(tmpwletf(imin:imax,1)));
zcsta1 = zcsta1+imin-1;
[~,imin] = min(tmpwlet(:,1));
[~,imax] = max(tmpwlet(:,1));
[~,zcsta2] = min(abs(tmpwlet(imin:imax,1)));
zcsta2 = zcsta2+imin-1;
greenlen = pow2(9)*sps/40;
green = zeros(greenlen,nsta); % no bandpass
greenf = zeros(greenlen,nsta);  % bandpassed version
ppeaks = zeros(nsta,1); % positive peaks
npeaks = zeros(nsta,1); % negative peaks
zcrosses = zeros(nsta,1);
ppeaksf = zeros(nsta,1); % positive peaks
npeaksf = zeros(nsta,1); % negative peaks
zcrossesf = zeros(nsta,1);
for ista = 1: nsta
  %cut according to the zero-crossing and the time shift from the constrained CC
  green(:,ista) = tmpwlet(zcsta2+8*sps-greenlen+1-offwlet1i(ista): zcsta2+8*sps-offwlet1i(ista), ista);
  greenf(:,ista) = tmpwletf(zcsta1+8*sps-greenlen+1-offwlet1if(ista): zcsta1+8*sps-offwlet1if(ista), ista);
  %detrend again for caution
  green(:,ista) = detrend(green(:,ista));
  greenf(:,ista) = detrend(greenf(:,ista));
  
  %re-find the zero-crossing as the template length has changed
  [~,imin] = min(green(:,ista));
  [~,imax] = max(green(:,ista));
  [~,zcrosses(ista)] = min(abs(green(imin:imax,ista)));
  zcrosses(ista) = zcrosses(ista)+imin-1;
  ppeaks(ista) = imax;
  npeaks(ista) = imin;
  
  %for bandpassed templates
  [~,imin] = min(greenf(:,ista));
  [~,imax] = max(greenf(:,ista));
  [~,zcrossesf(ista)] = min(abs(greenf(imin:imax,ista)));
  zcrossesf(ista) = zcrossesf(ista)+imin-1;
  ppeaksf(ista) = imax;
  npeaksf(ista) = imin;
  
end
%the following is just a check, because now the templates must be best aligned
ccmid = round(size(greenf,1)/2);
ccwlen = 4*sps;
loffmax = 5*sps/40;
ccmin = 0.01;  % depending on the length of trace, cc could be very low
iup = 1;    % times of upsampling
[off12con,off13con,cc] = constrained_cc_interp(greenf(:,1:3)',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
if ~(off12con==0 && off13con==0)
  disp('Filtered templates are NOT best aligned \n');
end
if nsta>3
  for ista = 4: nsta
    [mcoef,mlag] = xcorrmax(greenf(:,1), greenf(:,ista), mshiftadd, 'coeff');
    if mlag~=0   % offset in samples
      fprintf('Filtered templates are NOT best aligned at %s \n',stas(ista,:));
    end
  end
end
[off12con,off13con,cc] = constrained_cc_interp(green(:,1:3)',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
if ~(off12con==0 && off13con==0)
  disp('Broadband templates are NOT best aligned \n');
end
if nsta>3
  for ista = 4: nsta
    [mcoef,mlag] = xcorrmax(green(:,1), green(:,ista), mshiftadd, 'coeff');
    if mlag~=0   % offset in samples
      fprintf('Broadband templates are NOT best aligned at %s \n',stas(ista,:));
    end
  end
end

amprat(1,:) = minmax(greenf(:,1)')./minmax(greenf(:,2)');	% amp ratio between max at sta 3 and 2 or min
amprat(2,:) = minmax(greenf(:,1)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min
amprat(3,:) = minmax(greenf(:,2)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min
spread = range(greenf);   % range of the amp of template

%%%plot the unfiltered and filtered templates
% plt_templates(green,greenf,stas,[],[],lowlet,hiwlet,sps);

%just the filtered templates
% plt_templates_bp(greenf,stas,lowlet,hiwlet,sps);

% plt_templates(green(:,1:3),greenf(:,1:3),stas(1:3,:),[],[],lowlet,hiwlet,sps);
% zcrosses
% 2*(ppeaks-npeaks)/sps
%
% zcrossesf
% 2*(ppeaksf-npeaksf)/sps
% keyboard

%% load synthetic seismograms
%%%Specify the amplitude-frequency (counts) distribution
distr='UN'  % uniform distribution

fracelsew=0; %0.25 %0.5; %0.6; %The ratio of elsewhere events to local events.  Zero for Discrete Ide.

%%%specify regime for transformation from time offset to map location
% ftrans = 'interpArmb';
% ftrans = 'interpArmbreloc';
ftrans = 'interpchao';

%different percent of noise
perctrial = 0.1*(0:2:16)';
ntrial = length(perctrial);

%times of saturation
% nsat=[0.1 0.4 1 2 4 10 20 40 100];
nsat=[0.4 1 2 4 10 20 40 100];
nnsat = length(nsat);

nrun = 6;

%length of each simulation
sps = 160;
greenlen = pow2(9)*sps/40;
bufsec = 1;
msftaddm = bufsec*sps;  %buffer range for later CC alignment, +1 for safety
rccmwsec = 0.5;
rccmwlen = rccmwsec*sps;  %window length for computing RCC
overshoot = rccmwlen/2; %RCC points to the center of computing window, so extra room is needed
%   Twin=3*3600+(greenlen+msftaddm*2+overshoot*2-2)/sps; %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
Twin=0.5*3600+(greenlen+msftaddm*2+overshoot*2-2)/sps; %Early and late portions will be deleted. Twin includes only the early portion. In seconds.

tdura = 0.25;  %must be consistent with what synthetics actually used
% tdura = 0.4;

if tdura == 0.25
  %     savefile = strcat('rst_decon_synth_onespot','_td',num2str(tdura),'.mat');
  %     savefile = strcat('rst_decon_synth_onespotnitm8500','_td',num2str(tdura),'.mat'));
  savefile = strcat('rst_decon_synth_onespotmedwtcoef','_td',num2str(tdura),'.mat');
elseif tdura == 0.4
  %     savefile = ''rst_synth_onespot.mat';
  %     savefile = ''rst_synth_onespotnitm8500.mat';
  savefile = 'rst_decon_synth_onespotmedwtcoef.mat';
end
load(strcat(workpath,'/synthetics/',savefile));

% flagrecalc = 0;
flagrecalc = 1;
if flagrecalc

%%%loop for noise level
for iperc = 1: ntrial
  perc = perctrial(iperc);
  disp(perc);
  
  fname = ['/synthetics/STAS.',distr,'.',int2str(sps),'sps.',...
      'onespot','.noi',num2str(perc),'.else',num2str(fracelsew,2),...
      'nsat'];
  
  %%%loop for saturation level
  for insat = 1: nnsat
    %   insat = 1;
    disp(nsat(insat));

    sigstacat = [];   %cat preprocessed seismograms from all runs
    
    %loop for each run (each segment of synthetics to be combined later)
    for irun = 1: nrun
      irunstr = num2str(irun);
      
    if tdura == 0.25
      %%%load synthetics of certain saturation level, (time, station, sat)
      STAopt = load(strcat(workpath,fname,num2str(nsat(insat)),'tdura',...
        num2str(tdura),'T',num2str(round(Twin)),'p',irunstr));
      %%%load sources
      synsrc = load(strcat(workpath,fname,num2str(nsat(insat)),'tdura',...
        num2str(tdura),'T',num2str(round(Twin)),'p',irunstr,'_sources'));
      %%%load starting indices of added sources at sta 1
      synsrcstind = load(strcat(workpath,fname,num2str(nsat(insat)),'tdura',...
        num2str(tdura),'T',num2str(round(Twin)),'p',irunstr,'_stind'));
    elseif tdura == 0.4
      STAopt = load(strcat(workpath,fname,num2str(nsat(insat))));
      synsrc = load(strcat(workpath,fname,num2str(nsat(insat)),'_sources'));
      synsrcstind = load(strcat(workpath,fname,num2str(nsat(insat)),'_stind'));
    end 
    
    %% filter data (and templates)
    %%%filter data
    hisig=6.3; % this will give a similar spectral shape between template and signal
    losig=1.8;
    optseg = [];
    for ista = 1:nsta
      optseg(:,ista) = Bandpass(STAopt(:,ista), sps, losig, hisig, 2, 2, 'butter');
    end
    
    
    %% break into short windows
    %some params
    bufsec = 1;
    msftaddm = bufsec*sps;  %buffer range for later CC alignment, +1 for safety
    rccmwsec = 0.5;
    rccmwlen = rccmwsec*sps;  %window length for computing RCC
    overshoot = rccmwlen/2; %RCC points to the center of computing window, so extra room is needed
    
    indst = 1+msftaddm;
    inded = size(optseg,1)-msftaddm;
    subwsec = 25;
    subwlen = subwsec*sps;
    %since the rcc would lose rccmwlen/2 at both ends, this results in overlapping of 'rccmwlen' in rcc
    %across consecutive windows; if use 'ovlplen' of rccmwlen, then rcc has no overlapping at all
    %       ovlplen = rccmwlen*2;
    ovlplen = rccmwlen;
    windows = movingwins(indst,inded,subwlen,ovlplen,0);
    
    nwin =  size(windows,1);
    
    off1iw = zeros(nwin,nsta);  % the best alignment between sta2, sta3 wrt sta1 for each subwin
    ccaliw = zeros(nwin,1+nsta-3);  % CC value using the best alignment, including 4th stas
    
    ircccat = [];   % concatenated indices of RCC
    irccran = zeros(nwin,2);  % start and end indices (range) of RCC of all subwins
    rcccat = [];  % average concatenated RCC
    rcc1icat = [];  % concatenated RCC between sta 1 and 4
    rccpaircat = [];  % concatenated RCC between each station pair, order is 12, 13, 23
    mrccwpair = []; % median of concatenated RCC between each station pair
    ccwpair = []; % 0-lag overall cc of each subwin, between each station pair, order is 12, 13, 23
    ccw1i = []; % same as above, but between sta 1 and 4
    
    for iwin = 1: nwin
      isubwst = windows(iwin,1);
      isubwed = windows(iwin,2);
      
      %align records
      optcc = detrend(optseg(isubwst: isubwed,:));
      msftadd = (round(max(abs([off12ran off13ran])))+1)*sps/40;  %+1 for safety
      loffmax = 4*sps/40;
      ccmid = ceil(size(optcc,1)/2);
      ccwlen = round(size(optcc,1)-2*(msftadd+1));  % minus ensures successful shifting of records
      ccmin = 0.01;  % depending on the length of trace, cc could be very low
      iup = 1;    % times of upsampling
      [off12con,off13con,ccaliw(iwin,1)] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
        ccwlen,msftadd,loffmax,ccmin,iup);
      % if a better alignment cannot be achieved, use 0,0
      if off12con == msftadd+1 && off13con == msftadd+1
        fprintf('Short window %d cannot be properly aligned, double-check needed \n',iwin);
        off12con = 0;
        off13con = 0;
        cc12 = xcorr(optcc(:,1), optcc(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
        cc13 = xcorr(optcc(:,1), optcc(:,3),0,'normalized');
        cc23 = xcorr(optcc(:,2), optcc(:,3),0,'normalized');
        ccaliw(iwin,1) = (cc12+cc13+cc23)/3;
      end
      off1iw(iwin,1) = 0;
      off1iw(iwin,2) = round(off12con);
      off1iw(iwin,3) = round(off13con);
      
      for ista = 4: nsta
        [ccaliw(iwin,ista-2),off1iw(iwin,ista)] = xcorrmax(optcc(:,1),optcc(:,ista), 1.5*msftadd, 'coeff');
      end
      
      %Align records
      optdat = [];  % win +/-3 s, segment of interest,first 1s will be tapered
      optdat(:, 1) = optseg(isubwst: isubwed, 1); % time column
      for ista = 1: nsta
        optdat(:, ista) = optseg(isubwst-off1iw(iwin,ista): isubwed-off1iw(iwin,ista), ista);
      end
      
      subw = zeros(size(optdat,1), nsta);
      for ista = 1: nsta
        tmp = optdat(:,ista); %best aligned, filtered
        %detrend and taper only the data, NOT the noise
        tmp = detrend(tmp);
        subw(:,ista) = tmp;
      end
      %compute running CC between 3 stations
      [irccw,rccw12,rccw13,rccw23] = RunningCC3sta(subw,rccmwlen);
      rccw = (rccw12+rccw13+rccw23)/3;
      irccw = irccw + windows(iwin,1) - windows(1,1);   %convert to global index
      %         if iwin == 1
      irccw = irccw - overshoot;
      %         end
      %         plot(irccw, rccw);
      
      ircccat = [ircccat; irccw];
      irccran(iwin,:) = [irccw(1) irccw(end)];
      rcccat = [rcccat; rccw];
      rccpaircat = [rccpaircat; rccw12 rccw13 rccw23];
      mrccwpair = [mrccwpair; median(rccw12) median(rccw13) median(rccw23)];
      
      rccw1i = [];  % between 4th stas and 1st sta
      for ista = 4:nsta
        [~,rccw1i(:,ista-3)] = RunningCC(subw(:,1), subw(:,ista), rccmwlen);
      end
      rcc1icat = [rcc1icat; rccw1i];
      
      ccw12 = xcorr(subw(:,1), subw(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
      ccw13 = xcorr(subw(:,1), subw(:,3),0,'normalized');
      ccw23 = xcorr(subw(:,2), subw(:,3),0,'normalized');
      ccwpair = [ccwpair; ccw12 ccw13 ccw23];
      tmp = zeros(1,nsta-3);
      for ista = 4:nsta
        tmp(1,ista-3) = xcorr(subw(:,1), subw(:,ista),0,'normalized');
      end
      ccw1i = [ccw1i; tmp];
    end
    
    %if only use the mean RCC from pair 12 and 13
    rcccat = mean(rccpaircat(:,[1 2]), 2);
    
    %% Best alignment for the whole window
    %%%obtain a single best alignment based on the entire win
    optcc = detrend(optseg(1+msftaddm: end-msftaddm, :));
    %       msftadd = 10*sps/40;
    ccmid = ceil(size(optcc,1)/2);
    ccwlen = round(size(optcc,1)-2*(msftadd+1));  % minus ensures successful shifting of records
    ccmin = 0.01;  % depending on the length of trace, cc could be very low
    iup = 1;    % times of upsampling
    [off12con,off13con,ccali,iloopoff,loopoff] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
      ccwlen,msftadd,loffmax,ccmin,iup);
    % if a better alignment cannot be achieved, use 0,0
    if off12con == msftadd+1 && off13con == msftadd+1
      fprintf('Whole window cannot be properly aligned, double-check needed \n');
      off12con = 0;
      off13con = 0;
      cc12 = xcorr(optcc(:,1), optcc(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
      cc13 = xcorr(optcc(:,1), optcc(:,3),0,'normalized');
      cc23 = xcorr(optcc(:,2), optcc(:,3),0,'normalized');
      ccali = (cc12+cc13+cc23)/3;
    end
    off1i = zeros(nsta,1);
    off1i(2) = round(off12con);
    off1i(3) = round(off13con);
    
    %%%obtain a single 'observational' alignment between 4th and 1st station, note that this might
    %%%be very different from the empirical prediction from 'empioffset4thsta002'
    %%%--there are also different options here, using sig, or sig-wlet CC as in decon routine
    mcoef = zeros(nsta-3, 1);
    mlag = zeros(nsta-3, 1);
    for ista = 4: nsta
      [mcoef(ista-3),off1i(ista)] = xcorrmax(optcc(:,1), optcc(:,ista), msftadd, 'coeff');
      mlag(ista-3) = off1i(ista);
    end
    
    %       %%%%%%%%%%%%%%%%%%%%%%%%%%
    %       %if you want to avoid the case when the alignment is way off the centroid
    %       %by chance while the saturation level is low, you can force it to be the
    %       %an average location, this reference value is from the abs location of the
    %       %the centroid (0.2,0.2), and prediction of off14 from plane fit model
    %       off1i(2:3) = [2 2];
    %       [~,off1i(4)] = pred_tarvl_at4thsta(stas(4,:),off1i(2),off1i(3));
    %       %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%Align and compute the RCC based on the entire win, and take that as the input signal!
    optdat = [];  % win segment of interest
    for ista = 1: nsta
      optdat(:, ista) = optseg(1+msftaddm-off1i(ista): end-msftaddm-off1i(ista), ista);
    end
    
    %location of the whole-win best alignment
    loc0 = off2space002([off1i(2) off1i(3)],sps,ftrans,0);
    % loc has 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
    
    %%%2022/06/06, do NOT taper whatsoever!!
    sigsta = zeros(size(optdat,1), nsta);
    for ista = 1:nsta
      tmp = optdat(:,ista); %best aligned, filtered
      %detrend and taper only the data, NOT the noise
      tmp = detrend(tmp);
      sigsta(:,ista) = tmp;
    end
    %compute running CC between 3 stations
    [ircc,rcc12,rcc13,rcc23] = RunningCC3sta(sigsta,rccmwlen);
    ircc = ircc-overshoot;
    rccpair = [rcc12 rcc13 rcc23];
    %if only use the mean RCC from pair 12 and 13
    rcc = mean(rccpair(:,[1 2]), 2);
    mrcc = median(rcc);
    
    rcc1i = zeros(length(rcc),nsta-3);
    for ista = 4:nsta
      [~,rcc1i(:,ista-3)] = RunningCC(sigsta(:,1), sigsta(:,ista), rccmwlen);
    end
    
    sigsta = detrend(sigsta(overshoot+1:end-overshoot, :));  %excluding the overshoot
    
    %zero-lag max CC
    cc12 = xcorr(sigsta(:,1), sigsta(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
    cc13 = xcorr(sigsta(:,1), sigsta(:,3),0,'normalized');
    cc23 = xcorr(sigsta(:,2), sigsta(:,3),0,'normalized');
    cc1i = zeros(nsta-3,1);
    for ista = 4:nsta
      cc1i(ista-3) = xcorr(sigsta(:,1), sigsta(:,ista),0,'normalized');
    end
    ccpair = [cc12 cc13 cc23];
    mcc = (cc12+cc13)/2;  %if only use the pair 12 and 13
        
    sigstacat = [sigstacat; sigsta];
    
    end
    
    %%%envelope of the trace
    for i = 1:nsta
      [envup,~] = envelope(detrend(sigstacat(:,i)));
      prc = 95;
      envhi(i,insat,ireg) = prctile(envup,prc);
      envlo(i,insat,ireg) = prctile(envup,100-prc);
      envmed(i,insat,ireg) = prctile(envup,50);
    end
    
    % %a qucik look of the data
    % figure
    % lsig = size(sigsta,1);
    % plot((1:lsig)/sps,sigsta(:,1),'r'); hold on
    % plot((1:lsig)/sps,sigsta(:,2),'b');
    % plot((1:lsig)/sps,sigsta(:,3),'k');
    % ax = gca;
    % axsym(ax);
    % plot((1:lsig)/sps,rcc*ax.YLim(2),'o','color',[.6 .6 .6],'markersize',2);
    % text(0.95,0.1,sprintf('Saturation= %.1f; Semia= %.2f, Semib= %.2f',...
    %   nsat(insat),xaxis,yaxis),'Units','normalized','HorizontalAlignment','right');
    % xlim([0 20]);
    
  end
  % toc
end

end %if need to recalculate


%% RCC vs. saturation rate & region size
%%%load the LFE catalog, 25-s-win
% savefile = 'deconv_stats4th_allbstsig.mat';
savefile = 'deconv_stats4th_no23_allbstsig.mat';
allsig = load(strcat(rstpath, '/MAPS/',savefile));
rccbst = allsig.allbstsig.rccbst; %cat rcc; whole-win rcc; same for noise
ccwpairk = allsig.allbstsig.ccwpairk; %zero-lag CC of 25-s-window
mrccwpairk = allsig.allbstsig.mrccwpairk; %median RCC of 25-s-windows

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
for i = 1: length(trange)
  temp = rccbst{i};
  %%%cat rcc from 25-s-wins
  temp1 = sort(temp(:,1),'ascend');
  temp1(:,2) = (1:size(temp1,1))/size(temp1,1);
  templump = [templump; temp1];
  
  %%%rcc from the whole-win
  temp1 = sort(temp(:,2),'ascend');
  temp1(:,2) = (1:size(temp1,1))/size(temp1,1);
  templump1 = [templump1; temp1];
  
  %%%0-lag cc from the whole-win
  temp2 = mcc1win(i); % the mean of 12 and 13 of the whole win,
  templump2 = [templump2; temp2];
  
  %%%0-lag cc from 25-s-wins
  temp3 = ccwpairk{i};  %all 25-s subwins, pairs 12, 13, 23, as col
  temp4 = (temp3(:,1)+temp3(:,2))./2; %retain the mean of 12 and 13
  templump3 = [templump3; temp4]; %lump them
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

%%%same as above, whole-win version
rccplt = templump1(:,1:2); %data
% %%%1. bin x by equal width, then get the median of y
% [xcnt,ycnt] = ranybinx(rccplt(:,1),rccplt(:,2),'median',[],[],-1:0.01:1);
% rccreal1win = [xcnt,ycnt];
%%%2. bin x by equal number
nbin = 200;
[xbin,indbin,n] = binxeqnum(rccplt(:,1),nbin);
for i = 1: nbin
  ind = indbin{i};
  rccreal1win(i,1) = median(rccplt(ind,1));
  rccreal1win(i,2) = median(rccplt(ind,2));
end

%median lumped rcc from all data bursts, 25-s-win version
mrccreal = median(templump(:,1));
%median lumped rcc from all data bursts, whole-win version
mrccreal1win = median(templump1(:,1));
%median lumped 0-lag CC from all data bursts, 25-s-win version
mccreal = median(templump3);
%median lumped 0-lag CC from all data bursts, whole-win version
mccreal1win = median(templump2);

%% plot of median RCC/CC value wrt saturation & region size
rcccat = allsyn.rcccatk;
ccwpair = allsyn.ccwpairk;

f = initfig(8,4.5,1,2); %initialize fig
color = gradientblue(ntrial);
ax=f.ax(1);
hold(ax,'on'); ax.Box='on'; grid(ax,'on');
for iperc = 1: ntrial
  for insat = 1: nnsat
    mrcc(insat,iperc) = median(rcccat{insat,iperc});  %median of rcc
  end%   plot(ax,log10(nsat),mrcc,'-o','markersize',4,'color',color(iperc,:));
  plot(ax,log10(nsat),mrcc(:,iperc),'-','Color',color(iperc,:),'linew',1);
  scatter(ax,log10(nsat),mrcc(:,iperc),20,color(iperc,:),'filled','MarkerEdgeColor','k');
end
plot(ax,ax.XLim,[mrccreal mrccreal],'k--','linew',1);
xlabel(ax,'log_{10}(Saturation)');
ylabel(ax,'Median RCC');
longticks(ax,2);
ylim(ax,[0 1]);

ax=f.ax(2);
hold(ax,'on'); ax.Box='on'; grid(ax,'on');
p=[]; label=[];
for iperc = 1: ntrial
  for insat = 1: nnsat
    aa = ccwpair{insat,iperc}; %0-lag CC between 3 pairs from all 25-s wins
    bb = mean(aa(:,[1 2]), 2); %mean pof pairs 12 and 13
    mcc(insat,iperc) = median(bb);  %use median
  end
  %   plot(ax,log10(nsat),mcc(:,iperc),'-o','markersize',4,'color',color(iperc,:));
  p(iperc) = plot(ax,log10(nsat),mcc(:,iperc),'-','Color',color(iperc,:),'linew',1);
  scatter(ax,log10(nsat),mcc(:,iperc),20,color(iperc,:),'filled','MarkerEdgeColor','k');
  label{iperc} = sprintf('Noise=%.1f',perctrial(iperc));
end
p(ntrial+1) = plot(ax,ax.XLim,[mccreal mccreal],'k--','linew',1);
% plot(ax,ax.XLim,[mccreal1win mccreal1win],'b--');
label{ntrial+1} = 'Data';
% legend(ax,p,label,'NumColumns',2,'Location','north');
xlabel(ax,'log_{10}(Saturation)');
ylabel(ax,'Zero-lag CC');
longticks(ax,2);
ylim(ax,[0 1]);

keyboard

%% plot of the cdf of RCC wrt saturation & region size
% f = initfig(12,8,3,3); %initialize fig
% color = jet(nnsat);
% k=0;
% for iperc =  1: ntrial
%   k=k+1;
%   ax=f.ax(k);
%   hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%   p(1)=plot(ax,rccreal(:,1),rccreal(:,2),'k-','LineWidth',1);
%   label2{1}='Med of data';
%   for i = 1: nnsat
%     temp = sort(rccsat(:,i,iperc),'ascend');
%     temp(:,2) = (1:size(temp,1))/size(temp,1);
%     p(i+1) = plot(ax,temp(:,1),temp(:,2),'-','color',color(i,:));
%     label2{i+1} = sprintf('Satur=%.1f',nsat(i));
%   end
%   title(ax,sprintf('noise=%.1f',perctrial(iperc)));
%   if k==1
%     legend(ax,p,label2);
%   end
%   xlabel(ax,'RCC');
%   ylabel(ax,'CDF');
% end



