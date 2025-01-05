%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is highly identical to "deconv_ref_4s_exp_4thsta_181demo.m"
% which is designed specifically to burst 181. This script is to deal with
% ANY other single burst to reveal the detailed intermediate result and
% plots.
% Especially, I like to see the zoom-in seismograms in addition to the LFE
% clusters defined by the deconvolved 3- or 4-sta LFE catalog. One benefit
% is to see the correspondence. Another is to found the segments containing
% larger event clusters, which is needed for a plot showing more panels like
% the highlighted 25-s win of burst 181.
%
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/11/03
% Last modified date:   2024/11/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e   % Set the format to 5-digit floating point
clear
clc
close all

%% for easy testing
%%%Candidate bursts to run: indhi = [99; 100; 115; 164; 170; 174; 181];
%%%these are chosen according to amp binning from small to large
%bin 1: 079, 094, 011, 013
%bin 2: 085, 174,
%bin 3: 049, 174, 057
%bin 4: 164, 067, 158
%bin 5: 060, 181, 193
%%%finally chosen for paper: 079, 011, 085, 174, 060

defval('idxbst',79); %global indices of bursts to run
defval('normflag',0); %whether to normalize templates
defval('noiseflag',0);  %whether to use synthetic noises
defval('pltflag',0);  %whether to plot figs for each burst
defval('rccmwsec',0.5); %moving win len in sec for computing RCC

rccflag = 1; %1 means RCC weighting is used
whichrcc = 0; %if rcc weighting, which pair is used, 0 is best 2 pairs; 1 is 12; 2 is 13; 3 is 23

saveflag = 1; %0 means don't save figures; 1 means yes

%Choice to make upon the actual-used alignment at 4th stations
if noiseflag
  align14flag = 0;  %do NOT align sta 4 wrt. 1 if using noise
else
  align14flag = 1;
end

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

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
  'LZB  '
  'TWKB '
  'MGCB '
  'KLNB ']; % determine the trio and order, here the 1st sta is PGC
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

% dates and ets
%%% NOTE: 'dates' is the dates that the tremor is active at the region of interest, so that the
%%% waveform at the stations on these dates you are looking at mainly comes from the sources at the
%%% region of interest. We want to see if there is a noticable change in spectra during the burst
%%% windows on these dates
dates = unique(trange(:,1));
years = unique(floor(dates/1000));
nets = length(years);

%%%load the empirical model param of time offset 14 for each addtional stations
off14mod = load(strcat(rstpath, '/MAPS/timeoff_planefitparam_4thsta_160sps'),'w+');

%% load other catalogs
%%% load the LFE catalog of Michael Bostock, inside and outside the rectangle in 'locinterp002_4s.m'
%obtain the location of fam 002, lon0 and lat0
ftrans = 'interpchao';
% loc0 = off2space002([0 0],sps*iup,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
%format: [fam yyyy mm dd sec dx dy lon lat dep magnitude number-of-stations], 12 cols
%time should point to the peak, zero-crossing?
% bostcat = ReformBostock(loc0(3),loc0(4),0);

bosdir = ('/home/data2/chaosong/matlab/allan/BOSTOCK');
lfefnm = ('newlfeloc');
lfeloc = load(fullfile(bosdir, lfefnm));
loc0 = lfeloc(lfeloc(:,1)==2,:);
% bostcat = ReformBostock(loc0(3),loc0(2),1);
%%%if use the lumped catalog that combine unique events from 002 and 246
bostcat = load(fullfile(bosdir, '002-246_lumped.2003-2005_cull_NEW_chao'));

%all MB's LFEs of the SAME dates that I have been using
bostdayi = cell(length(dates),1);
for i = 1: length(dates)
  date = dates(i);
  year = floor(date/1000);
  jday = floor(date-year*1000);
  a = jul2dat(year,jday);
  
  temp = bostcat(bostcat(:,2)==year & bostcat(:,3)==a(1) & bostcat(:,4)==a(2),:);
  bostdayi{i} = sortrows(temp, 5);
  nbostdayi(i) = size(temp,1);
  bomsumday(i) = sum(temp(:,end));
  
end
bostdayia = cat(1,bostdayi{:});
bostcat = bostdayia;

%the cut-out boundary of 4-s detections
cutout = 'ellipse';
x0 = 0.2;
y0 = 0.2;
semia = 1.75;
semib = 1.25;
angrot = 45;
[xcut, ycut] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);

%same rectangle in 'locinterp002_4s.m'
EW = [-7 3];
NS = [-3 4];
wid = range(EW);
hgt = range(NS);
x0 = mean(EW);
y0 = mean(NS);
[x, y] = rectangle_chao(x0,y0,wid,hgt,0.01);

%%%2022/06/29, use the same ellipse to exclude fam 047
% bnd = [x y];
bnd = [xcut ycut];
[iin,ion] = inpolygon(bostcat(:,6),bostcat(:,7),bnd(:,1),bnd(:,2));
isinbnd = iin | ion;
bostcati = bostcat(isinbnd == 1, :);
bostcato = bostcat(isinbnd ~= 1, :);
clear bostcat

% %if you choose only fam 002 regardless
% bostcati = bostcati(bostcati(:,1)==2,:);
%convert moment mag to moment, Mw = (2/3)*log_10(M0)-10.7, where M0 has the unit of dyne.cm
%(10^-7 N.m), so Mw = (2/3)*log_10(M0)-6 if M0 has the unit of N.m
bostcati(:,13) = 10.^(1.5*(6+bostcati(:,11)));

% bostcato = bostcato(bostcato(:,1)~=2 & bostcato(:,1)~=47 & bostcato(:,1)~=246,:);
bostcato(:,13) = 10.^(1.5*(6+bostcato(:,11)));


%%%load the tremor catalog of John Armbruster,
%%%2022/06/29, not really very useful as the detecting window could be 128-s long (not sure)
%format: [yyyy mm dd sec dx dy lon lat dep], 9 cols;
armcat = ReformArmbrusterv2(loc0(3),loc0(4),0);
[~,ind] = min(abs(armcat(:,5))+abs(armcat(:,6)));
% armcat1 = ReformArmbruster(loc0(3),loc0(4),1);
% [~,ind1] = min(abs(armcat1(:,5))+abs(armcat1(:,6)));
[iin,ion] = inpolygon(armcat(:,5),armcat(:,6),bnd(:,1),bnd(:,2));
isinbnd = iin | ion;
armcati = armcat(isinbnd == 1, :);
armcato = armcat(isinbnd ~= 1, :);
clear armcat

%% prepare templates (Green's functions), from 'lfetemp002_160sps.m'
sps = 160;
templensec = 60;

ccstack = [];
for ista = 1: nsta
  fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
    num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao_catnew');
  ccstack(:,ista) = load(fname);
end
STA = detrend(ccstack);

ccstackort = [];
for ista = 1: nsta
  fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
    num2str(templensec), 's_', 'BBCCS_', 'ort_Nof_Non_Chao_catnew');
  ccstackort(:,ista) = load(fname);
end
STAort = detrend(ccstackort);

ccstackvert = [];
for ista = 1: nsta
  fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
    num2str(templensec), 's_', 'BBCCS_', 'vert_Nof_Non_Chao_catnew');
  ccstackvert(:,ista) = load(fname);
end
STAvert = detrend(ccstackvert);

%flag of normalization
% normflg = 0;

% %plot the raw templates, not filtered, not best aligned
% figure
% subplot(211)
% hold on
% for ista = 1: nsta
%   plot(STA(:,ista));
% end
% subplot(212)
% hold on
% for ista = 1: nsta
%   plot(STAort(:,ista));
% end

%%%The below aligns the templates by x-correlation
ist = templensec*sps*4/10;  %not using whole window in case any station has very long-period energy
ied = templensec*sps*6/10;
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
  STAtmport(:,ista)=detrend(STAort(is(ista):ie(ista),ista));
  STAtmpvert(:,ista)=detrend(STAvert(is(ista):ie(ista),ista));
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
  STAtmport(mshiftadd+1:end-(mshiftadd+1),ista)=detrend(STAtmport(mshiftadd+1-imax(ista):end-(mshiftadd+1)-imax(ista),ista));
  STAtmpvert(mshiftadd+1:end-(mshiftadd+1),ista)=detrend(STAtmpvert(mshiftadd+1-imax(ista):end-(mshiftadd+1)-imax(ista),ista));
end
%normalization
if normflag
  for ista=1:nsta
    STAtmp(:,ista)=STAtmp(:,ista)/spread(ista); % now templates are 'aligned' indeoendently by x-corr wrt. sta 1
    STAtmport(:,ista)=STAtmport(:,ista)/spread(ista);
    STAtmpvert(:,ista)=STAtmpvert(:,ista)/spread(ista);
  end
end
% figure
% subplot(211)
% hold on
% for ista = 1: nsta
%   plot(STAtmp(:,ista));
% end
% subplot(212)
% hold on
% for ista = 1: nsta
%   plot(STAtmport(:,ista));
% end
%%%The above aligns the templates by x-correlation

%%%detrend, taper and bandpass templates
tmpwlet = STAtmp; % no bandpass
tmpwletf = STAtmp;  % bandpassed version
fractap = sps/size(tmpwlet,1);
tmpwletort = STAtmport; % no bandpass
tmpwletfort = STAtmport;  % bandpassed version
tmpwletvert = STAtmpvert; % no bandpass
tmpwletfvert = STAtmpvert;  % bandpassed version
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
  
  %same process for orthogonal
  tmpwletort(:,ista) = detrend(tmpwletort(:,ista));
  tmpwletort(:,ista) = w.* tmpwletort(:,ista);
  tmpwletort(:,ista) = detrend(tmpwletort(:,ista));
  tmpwletfort(:,ista) = Bandpass(tmpwletort(:,ista), sps, lowlet, hiwlet, 2, 2, 'butter');
  tmpwletfort(:,ista) = detrend(tmpwletfort(:,ista));
  
  %same process for vertical
  tmpwletvert(:,ista) = detrend(tmpwletvert(:,ista));
  tmpwletvert(:,ista) = w.* tmpwletvert(:,ista);
  tmpwletvert(:,ista) = detrend(tmpwletvert(:,ista));
  tmpwletfvert(:,ista) = Bandpass(tmpwletvert(:,ista), sps, lowlet, hiwlet, 2, 2, 'butter');
  tmpwletfvert(:,ista) = detrend(tmpwletfvert(:,ista));
  
end

%%%constrained CC, so that only 2 offsets are independent
ccmid = round(size(tmpwletf,1)/2);
ccwlen = 10*sps;
loffmax = 5*sps/40;
ccmin = 0.01;  % depending on the length of trace, cc could be very low
iup = 1;    % times of upsampling
[off12con,off13con,cc] = constrained_cc_interp(tmpwletf(:,1:3)',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
offwlet1i(1) = 0;
offwlet1i(2) = round(off12con);
offwlet1i(3) = round(off13con);

for ista = 4: nsta
  [mcoef,offwlet1i(ista)] = xcorrmax(tmpwletf(:,1), tmpwletf(:,ista), mshiftadd, 'coeff');
end

%%%automatically find the rough zero-crossing time, whose abs. value is closest to 0, whether + or -
[~,imin] = min(tmpwletf(:,1));
[~,imax] = max(tmpwletf(:,1));
[~,zcsta1] = min(abs(tmpwletf(imin:imax,1)));
zcsta1 = zcsta1+imin-1;
greenlen = pow2(9)*sps/40;
green = zeros(greenlen,nsta); % no bandpass
greenf = zeros(greenlen,nsta);  % bandpassed version
ppeaks = zeros(nsta,1); % positive peaks
npeaks = zeros(nsta,1); % negative peaks
greenort = zeros(greenlen,nsta); % no bandpass
greenfort = zeros(greenlen,nsta);  % bandpassed version
greenvert = zeros(greenlen,nsta); % no bandpass
greenfvert = zeros(greenlen,nsta);  % bandpassed version
for ista = 1: nsta
  %cut according to the zero-crossing and the time shift from the constrained CC
  green(:,ista) = tmpwlet(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenf(:,ista) = tmpwletf(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  %detrend again for caution
  green(:,ista) = detrend(green(:,ista));
  greenf(:,ista) = detrend(greenf(:,ista));
  if normflag
    %normalize by max amp
    green(:,ista) = green(:,ista)/max(abs(green(:,ista)));    % normalize
    greenf(:,ista) = greenf(:,ista)/max(abs(green(:,ista)));    % normalize
  end
  
  %same process for orthogonal
  greenort(:,ista) = tmpwletort(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenfort(:,ista) = tmpwletfort(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenort(:,ista) = detrend(greenort(:,ista));
  greenfort(:,ista) = detrend(greenfort(:,ista));
  if normflag
    greenort(:,ista) = greenort(:,ista)/max(abs(green(:,ista)));    % normalize
    greenfort(:,ista) = greenfort(:,ista)/max(abs(green(:,ista)));    % normalize
  end
  
  %same process for vertical
  greenvert(:,ista) = tmpwletvert(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenfvert(:,ista) = tmpwletfvert(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenvert(:,ista) = detrend(greenvert(:,ista));
  greenfvert(:,ista) = detrend(greenfvert(:,ista));
  if normflag
    greenvert(:,ista) = greenvert(:,ista)/max(abs(green(:,ista)));    % normalize
    greenfvert(:,ista) = greenfvert(:,ista)/max(abs(green(:,ista)));    % normalize
  end
  
  %re-find the zero-crossing as the template length has changed
  [~,imin] = min(greenf(:,ista));
  [~,imax] = max(greenf(:,ista));
  [~,zcrosses(ista)] = min(abs(greenf(imin:imax,ista)));
  zcrosses(ista) = zcrosses(ista)+imin-1;
  ppeaks(ista) = imax;
  npeaks(ista) = imin;
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
for ista = 4: nsta
  [mcoef,mlag] = xcorrmax(greenf(:,1), greenf(:,ista), mshiftadd, 'coeff');
  if mlag~=0   % offset in samples
    fprintf('Filtered templates are NOT best aligned at %s \n',stas(ista,:));
  end
end

amprat(1,:) = minmax(greenf(:,1)')./minmax(greenf(:,2)');	% amp ratio between max at sta 3 and 2 or min
amprat(2,:) = minmax(greenf(:,1)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min
amprat(3,:) = minmax(greenf(:,2)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min
spread = range(greenf);   % range of the amp of template

%%%plot the unfiltered templates at diff stas in one panel, and filtered in another
% plt_templates(green,greenf,stas,greenort,greenfort,lowlet,hiwlet,sps);
% plt_templates(green,greenf,stas,[],[],lowlet,hiwlet,sps);

%just the filtered templates
% plt_templates_bp(greenf,stas,lowlet,hiwlet,sps);

%plot the filtered templates at all components
staplt = 1:7;
winnoi = [1 3]; %empirical window for pre-arrival noise
winsig = [4 7]; %empirical window for main arrival
% [f,snrf,snrfort,snrfvert]=plt_templates_bp_bysta(greenf,greenfort,greenfvert,stas,staplt,lowlet,hiwlet,sps,...
%   [winnoi;winsig]);


%% prepare the signal and noise windows
%filtering passband for reading data, confirmed by 'spectrabursts002_4s.m'
hisig=6.3; % this will give a similar spectral shape between template and signal
losig=1.8;

%%%moving window length in samples for running CC, envelope, etc.
%standard window length is about 0.5s, this is about the visual duration of the filtered and unfiltered
%template, although in fact to include as least one cycle of the main dipole of template
rccmwlen=rccmwsec*sps;
% rccmwlen=sps/2;
% rccmwlen=sps;

windowsk = cell(size(trange,1),1); %window division
irccrank = cell(size(trange,1),1); %indices of rcc concatenation
off1iwk = cell(size(trange,1),1);  % the best alignment between sta 2, 3 wrt 1 for all subwins and all burst wins
off1ic = zeros(size(trange,1),nsta);  % single best alignment 'computed' between ALL stas wrt 1 for entire win
off1i = zeros(size(trange,1),nsta);  % single best alignment 'actually used' between ALL stas wrt 1 for entire win
off14pred = zeros(size(trange,1),nsta-3); %empirical pred of off14 from plane fit given single best alignment
ccali = zeros(size(trange,1),1);  % CC value using the best alignment
ccaliwk = cell(size(trange,1),1); % CC value using the best alignment for all subwins and all burst wins
ccwpairk = cell(size(trange,1),1); % 0-lag CC value for all subwins and all burst wins
mrccwpairk = cell(size(trange,1),1); % median RCC value for all subwins and all burst wins
subwsec = zeros(size(trange,1),1);  % subwin length in sec used in practice
subwseclfit = zeros(size(trange,1),1);  % subwin length from linear fitting, ie, time for 1-sample offset change
ndetwink = cell(size(trange,1),1); % num of detections in each sub win after 2ndary removed
ndetwin4thk = cell(size(trange,1),1); % num of detections in each sub win after 4th sta check

for iii = 1: length(idxbst)
  
  [iets,i,j] = indofburst(trange,idxbst(iii));
  
  % for iets = 3: nets
  % dates in each ets
  year = years(iets);
  datesets = dates(floor(dates/1000)==year);
  
  %   for i = 2: length(datesets)
  
  date = datesets(i);
  jday = floor(date-year*1000);
  a = jul2dat(year,jday);
  if a(1) == 9
    mo = 'Sep.';
  elseif a(1) == 7
    mo = 'Jul.';
  else
    mo = 'Mar.';
  end
  dy = num2str(a(2));
  yr = num2str(a(3));
  
  %Bostock's LFE catalog on the same date
  bostdayi = bostcati(bostcati(:,2)==year & bostcati(:,3)==a(1) & bostcati(:,4)==a(2),:);
  bostdayi = sortrows(bostdayi, 5);
  bostdayo = bostcato(bostcato(:,2)==year & bostcato(:,3)==a(1) & bostcato(:,4)==a(2),:);
  bostdayo = sortrows(bostdayo, 5);
  
  %Armbruster's tremor catalog on the same date
  armdayi = armcati(armcati(:,1)==year & armcati(:,2)==a(1) & armcati(:,3)==a(2),:);
  armdayi = sortrows(armdayi, 4);
  armdayo = armcato(armcato(:,1)==year & armcato(:,2)==a(1) & armcato(:,3)==a(2),:);
  armdayo = sortrows(armdayo, 4);
  
  %bursts and 4-s detections of the same day
  rangetemp = trange(trange(:,1)==datesets(i), :);
  hfdayi = hfbnd(hfbnd(:,daycol)==datesets(i), :);  % inside bound of the day
  hfdayo = hfout(hfout(:,daycol)==datesets(i), :);  % outside bound of the day
  
  %read horizontal optimal and orthogonal components
  JDAY = num2zeropadstr(jday,3);
  MO=day2month(jday,year);     % EXTERNAL function, day2month, get the month of one particular date
  direc=[datapath, '/arch', yr,'/',MO,'/'];     % directory name
  prename=[direc,yr,'.',JDAY,'.00.00.00.0000.CN'];    %  path plus prefix of data file,
  %     disp(prename);
  [STAopt,STAort,~,fileflag] = rd_daily_bpdata(year,jday,prename,stas,PERMSTA,POLSTA,...
    PERMROTS,POLROTS,sps,losig,hisig,npo,npa,[],[],[],[]);
  
  if fileflag == 0    % means there are missing files
    fprintf('Day %s / %s will be omitted because of missing files. \n', yr, JDAY);
    continue    % continue to the next day
  end
  
  %read vertical components too, as we want to have a sense of the SNR at the orthogonal and vertical
  %components
  [STAvert,~,~,fileflag] = rd_daily_bpdata(year,jday,prename,stas,PERMSTA,POLSTA,...
    PERMROTS,POLROTS,sps,losig,hisig,npo,npa,[],[],[],[],'Z');
  
  %     keyboard
  %%
  %     for j = 14: size(rangetemp,1)
  %       close all
  %       k = k+1;
  k = idxbst(iii);
  disp(k);
  
  tmaxi = hfdayi(:, seccol); % starting time of max power rate of half sec inside the ellipse
  tmaxo = hfdayo(:, seccol); % starting time of max power rate of half sec outside the ellipse
  tcnti = hfdayi(:, 15);  % the center of detecting win is the 15th col
  tcnto = hfdayo(:, 15);  % the center of detecting win is the 15th col
  tbosti = bostdayi(:,5); % (arrival) time of Bostock's LFE catalog inside the rectangle
  tbosto = bostdayo(:,5); % (arrival) time of Bostock's LFE catalog outside the rectangle
  tarmi = armdayi(:,4); % (arrival) time of Armbruster's tremor catalog inside the rectangle
  tarmo = armdayo(:,4); % (arrival) time of Armbruster's tremor catalog outside the rectangle
  
  tst = rangetemp(j,2); % start and end time of bursts
  ted = rangetemp(j,3);
  
  %how many 4-s detections fall into the burst range
  indtmaxi = find(tmaxi>=tst-0.1 & tmaxi<=ted+0.1);
  ninbst(k,1) = length(indtmaxi); %+-0.1 s in case of resolution while saving to file
  
  %       %%%%Use a fixed range of time before and after the 0.5-s strongest arrival
  %       tbuffer = 3;   % buffer time to include some coherent precursor or coda, first 1s will be tapered
  %       tstbuf = tst-tbuffer; % start and end time of bursts, buffer added
  %       tedbuf = ted+tbuffer;
  %%%%Use the start and end of the 4-s detecting window
  tstbuf = min(tcnti(indtmaxi)-2);
  tedbuf = max(tcnti(indtmaxi)+2);
  tlenbuf = tedbuf-tstbuf;
  trangenew(iii,2)=tstbuf;
  trangenew(iii,3)=tedbuf;
  
  %plot the bostock's LFE catalog inside ellipse
  indbi = find(tbosti>=tstbuf & tbosti<=tedbuf);
  
  %max allowable shift in best alignment
  %       msftaddm = (round(max(abs([off12ran off13ran])))+1)*sps/40;  %+1 for safety
  %       msftaddm = sps+1;  %+1 for safety
  msftaddm = 1.5*sps+1;  %+1 for safety
  %       msftaddm = round(sps/8);    % maximum allowed shift between 2 traces
  
  %have some overshoot, so that the resulted rcc would have the same length as the signal
  overshoot = rccmwlen/2;
  %       overshoot = 0;
  
  %%%%2022/09/26, obtain all information for data before decon, so that you know the threshold
  %%%%being used, although this was used mainly for noise experiment, we still use it here
  %chop a record segment
  optseg = STAopt(max(floor(tstbuf*sps+1-overshoot-msftaddm),1): ...
    min(floor(tedbuf*sps+overshoot+msftaddm),86400*sps), :); % sta 1
  ortseg = STAort(max(floor(tstbuf*sps+1-overshoot-msftaddm),1): ...
    min(floor(tedbuf*sps+overshoot+msftaddm),86400*sps), :);
  vertseg = STAvert(max(floor(tstbuf*sps+1-overshoot-msftaddm),1): ...
    min(floor(tedbuf*sps+overshoot+msftaddm),86400*sps), :);
  
  %%%according the linear fitting result of off12 and off13 VS. origin time, we sort of know how
  %%%much it needs to change the overall offset to change by 1 sample
  %generate overlapping windows of the same length
  indst = floor(tstbuf*sps+1);
  inded = floor(tedbuf*sps);
  if subwseclfit(k)< 25 && subwseclfit(k)>0
    subwsec(k) = subwseclfit(k);
  elseif tlenbuf<25
    subwsec(k) = tlenbuf;
  else
    subwsec(k) = 25;   % determined from fitting
  end
  subwlen = subwsec(k)*sps;
  %since the rcc would lose rccmwlen/2 at both ends, this results in overlapping of 'rccmwlen' in rcc
  %across consecutive windows; if use 'ovlplen' of rccmwlen, then rcc has no overlapping at all
  %       ovlplen = rccmwlen*2;
  ovlplen = rccmwlen;
  windows = movingwins(indst,inded,subwlen,ovlplen,0);
  %       %if some portions are not included
  %       if windows(end,2) < inded
  %         %make it an separate window if its length is at least half of that of the others
  %         if inded-(windows(end,1)+(subwlen-ovlplen))+1 >= subwlen/2
  %           windows = [windows; [windows(end,1)+(subwlen-ovlplen) inded]];
  %         %otherwise just combine the portion to the last window
  %         else
  %           windows(end,2) = inded;
  %         end
  %       end
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
  
  rccpaircatx = [];
  
  for iwin = 1: nwin
    %         isubwst = windows(iwin,1);
    %         isubwed = windows(iwin,2);
    isubwst = windows(iwin,1)-max(floor(tstbuf*sps+1-overshoot-msftaddm),1);
    isubwed = windows(iwin,2)-max(floor(tstbuf*sps+1-overshoot-msftaddm),1);
    if iwin == 1
      isubwst = isubwst-overshoot;
    end
    if iwin == nwin
      isubwed = isubwed+overshoot;
    end
    
    %align records
    %         optcc = [];
    %         optcc(:,1) = STAopt(max(isubwst,1): min(isubwed,86400*sps), 2);
    %         optcc(:,2) = STAopt(max(isubwst*sps,1): min(isubwed*sps,86400*sps), 3);
    %         optcc(:,3) = STAopt(max(isubwst*sps,1): min(isubwed*sps,86400*sps), 4);
    optcc = detrend(optseg(isubwst: isubwed, 2:end));
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
      off12con = 0;
      off13con = 0;
      cc12 = xcorr(optcc(:,1), optcc(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
      cc13 = xcorr(optcc(:,1), optcc(:,3),0,'normalized');
      cc23 = xcorr(optcc(:,2), optcc(:,3),0,'normalized');
      ccaliw(iwin,1) = (cc12+cc13+cc23)/3;
      fprintf('Tremor burst %d cannot be properly aligned, double-check needed \n',k);
    end
    off1iw(iwin,1) = 0;
    off1iw(iwin,2) = round(off12con);
    off1iw(iwin,3) = round(off13con);
    
    for ista = 4: nsta
      [ccaliw(iwin,ista-2),off1iw(iwin,ista)] = xcorrmax(optcc(:,1),optcc(:,ista), 1.5*msftadd, 'coeff');
    end
    
    %Align records
    optdat = [];  % win +/-3 s, segment of interest,first 1s will be tapered
    ortdat = [];
    optdat(:, 1) = optseg(isubwst: isubwed, 1); % time column
    ortdat(:, 1) = ortseg(isubwst: isubwed, 1);
    for ista = 1: nsta
      optdat(:, ista+1) = optseg(isubwst-off1iw(iwin,ista): isubwed-off1iw(iwin,ista), ista+1);
      ortdat(:, ista+1) = ortseg(isubwst-off1iw(iwin,ista): isubwed-off1iw(iwin,ista), ista+1);
    end
    
    subw = zeros(size(optdat,1), nsta);
    for ista = 1: nsta
      tmp = optdat(:,ista+1); %best aligned, filtered
      %detrend and taper only the data, NOT the noise
      tmp = detrend(tmp);
      %%%2022/06/06, do NOT taper whatsoever!!
      %           %taper the very start or end subwin and obtain the new rcc between tapered signals
      %           ltmp = length(tmp);
      %           fractap = sps/ltmp; % if fractap is >=1, n-point von Hann window is returned
      %           ptstap = fractap*ltmp; % if fractap is >=1, n-point von Hann window is returned
      %           if iwin==1
      %             w = coswinlh(ltmp,fractap);
      %             tmp = w.* tmp;
      %             tmp = detrend(tmp); %detrend again for caution
      %           elseif iwin==nwin
      %             w = coswinrh(ltmp,fractap);
      %             tmp = w.* tmp;
      %             tmp = detrend(tmp); %detrend again for caution
      %           end
      subw(:,ista) = tmp;
    end
    %compute running CC between 3 stations
    [irccw,rccw12,rccw13,rccw23] = RunningCC3sta(subw,rccmwlen);
    rccw = (rccw12+rccw13+rccw23)/3;
    irccw = irccw + windows(iwin,1) - windows(1,1);   %convert to global index
    if iwin == 1
      irccw = irccw - overshoot;
    end
    %         plot(irccw, rccw);
    
    ircccat = [ircccat; irccw];
    irccran(iwin,:) = [irccw(1) irccw(end)];
    rcccat = [rcccat; rccw];
    rccpaircat = [rccpaircat; rccw12 rccw13 rccw23];
    mrccwpair = [mrccwpair; median(rccw12) median(rccw13) median(rccw23)];
    
    %         %%%%%%%%%%%%%%%%%% if detrend is needed for each rcc win
    %         tic
    %         [irccwx,rccw12x,rccw13x,rccw23x] = RunningCC3sta(subw,rccmwlen,'xcorr');
    %         toc
    %         irccwx = irccwx + windows(iwin,1) - windows(1,1);   %convert to global index
    %         rccpaircatx = [rccpaircatx; rccw12x rccw13x rccw23x];
    %         %%%%%%%%%%%%%%%%%% if detrend is needed for each rcc win
    
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
    
    
    %%%
    subwtmp = detrend(subw(overshoot+1:end-overshoot, :));
    subwsave{iwin} = subwtmp;
  end
  irccrank{k} = irccran;
  windowsk{k} = windows-windows(1,1)+1;
  off1iwk{k} = off1iw;
  ccaliwk{k} = ccaliw;
  ccwpairk{k} = ccwpair;
  mrccwpairk{k} = mrccwpair;
  % keyboard
  %       %if only use the mean RCC from the 2 pairs that have the highest overall CC
  %       [~,ind] = min(sum(ccwpair,1));
  %       rcccat = mean(rccpaircat(:,setdiff(1:3,ind)), 2);
  
  %if only use the mean RCC from pair 12 and 13
  rcccat = mean(rccpaircat(:,[1 2]), 2);
  
  %       %%
  %       rcccatx = mean(rccpaircatx(:,[1 2]), 2);
  %       figure
  %       for ii=1:5
  %         subplot(5,1,ii); hold on
  %         plot(ircccat/sps,rcccat,'k-');
  %         plot(ircccat/sps,rcccatx,'r-');
  %         xlim(50*[(ii-1) ii]);
  %       end
  
  %%
  %if choose not to use RCC weighting; for easier comparison
  if ~rccflag
    rcccat = ones(size(rcccat));
    rcc1icat = ones(size(rcc1icat));
  else
    if whichrcc == 1
      rcccat = rccpaircat(:,1);
    elseif whichrcc == 2
      rcccat = rccpaircat(:,2);
    elseif whichrcc == 3
      rcccat = rccpaircat(:,3);
    end
  end
  
  rcccomb = [];
  rcccomb(:,1) = rcccat;
  
  %       figure
  %       subplot(231); hold on; ax=gca;
  %       histogram(rccpaircat(:,1)); title('rcc12'); [MUHAT,SIGMAHAT] = normfit(rccpaircat(:,1));
  %       plot([MUHAT MUHAT],ax.YLim,'r--');
  %       subplot(232); hold on; ax=gca;
  %       histogram(rccpaircat(:,2)); title('rcc13'); [MUHAT,SIGMAHAT] = normfit(rccpaircat(:,2));
  %       plot([MUHAT MUHAT],ax.YLim,'r--');
  %       subplot(233); hold on; ax=gca;
  %       histogram(rccpaircat(:,3)); title('rcc23'); [MUHAT,SIGMAHAT] = normfit(rccpaircat(:,3));
  %       plot([MUHAT MUHAT],ax.YLim,'r--');
  %       subplot(234); hold on; ax=gca;
  %       histogram(sum(rccpaircat(:,[1 2]),2)); title('rcc12+rcc13'); [MUHAT,SIGMAHAT] = normfit(sum(rccpaircat(:,[1 2]),2));
  %       plot([MUHAT MUHAT],ax.YLim,'r--');
  %       subplot(235); hold on; ax=gca;
  %       histogram(sum(rccpaircat(:,[2 3]),2)); title('rcc13+rcc23'); [MUHAT,SIGMAHAT] = normfit(sum(rccpaircat(:,[2 3]),2));
  %       plot([MUHAT MUHAT],ax.YLim,'r--');
  %       subplot(236); hold on; ax=gca;
  %       histogram(sum(rccpaircat(:,[1 3]),2)); title('rcc12+rcc23'); [MUHAT,SIGMAHAT] = normfit(sum(rccpaircat(:,[1 3]),2));
  %       plot([MUHAT MUHAT],ax.YLim,'r--');
  
  %%%Is it true that the coherence between 2-3 is the highest among 3 pairs?
  %%%---Yes, for the concatenated rcc
  %       [f] = plt_rcccat(rccpaircat,sps);
  
  %%%obtain a single best alignment based on the entire win
  %       optcc = optseg(:, 2:end);
  optcc = detrend(optseg(1+msftaddm: end-msftaddm, 2:end));
  ccmid = ceil(size(optcc,1)/2);
  ccwlen = round(size(optcc,1)-2*(msftadd+1));
  ccmin = 0.01;  % depending on the length of trace, cc could be very low
  iup = 1;    % times of upsampling
  [off12con,off13con,ccali(k)] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
    ccwlen,msftadd,loffmax,ccmin,iup);
  % if a better alignment cannot be achieved, use 0,0
  if off12con == msftadd+1 && off13con == msftadd+1
    off12con = 0;
    off13con = 0;
    cc12 = xcorr(optcc(:,1), optcc(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
    cc13 = xcorr(optcc(:,1), optcc(:,3),0,'normalized');
    cc23 = xcorr(optcc(:,2), optcc(:,3),0,'normalized');
    ccali(k) = (cc12+cc13+cc23)/3;
    fprintf('Tremor burst %d cannot be properly aligned, double-check needed \n',k);
  end
  off1ic(k,1) = 0;
  off1ic(k,2) = round(off12con);
  off1ic(k,3) = round(off13con);
  
  %%%obtain a single 'observational' alignment between 4th and 1st station, note that this might
  %%%be very different from the empirical prediction from 'empioffset4thsta002'
  %%%--there are also different options here, using sig, or sig-wlet CC as in decon routine
  mcoef = zeros(nsta-3, 3);
  mlag = zeros(nsta-3, 3);
  envrat = zeros(nsta-3, 3);
  for ista = 4: nsta
    [mcoef(ista-3, 1),off1ic(k,ista)] = xcorrmax(optcc(:,1), optcc(:,ista), msftadd, 'coeff');
    mlag(ista-3, 1) = off1ic(k,ista);
    envrat(ista-3, 1) = median(envelope(optcc(:,ista)))./median(envelope(optcc(:,1)));
    %do an overall CC between 4th and 2nd/3rd stas, to see which one they are most coehrent with
    for jjj = 2:3
      [mcoef(ista-3, jjj),mlag(ista-3, jjj)] = xcorrmax(optcc(:,1), optcc(:,ista), msftadd, 'coeff');
      envrat(ista-3, jjj) = median(envelope(optcc(:,ista)))./median(envelope(optcc(:,jjj)));
    end
    %empirical prediction from plane fitting in 'empioffset4thsta002'
    off14pred(k,ista-3) = round(off14mod(ista-3,1).*off1ic(k,2) + off14mod(ista-3,2).*off1ic(k,3) + ...
      off14mod(ista-3,3));
  end
  
  %for real data, USE the best whole-win alignment before decon
  off1i(k,1:3) = off1ic(k,1:3);
  %       off1i(k,1:3) = [0 0 0];
  
  %Choice to make upon the actual-used alignment at 4th stations
  %for data case, DO align!
  %       align14flag = 1;
  if align14flag
    off1i(k,4:end) = off1ic(k,4:end); %if you trust overall alignment at 4th stations
  else
    off1i(k,4:end) = zeros(1,nsta-3); %if you don't
  end
  
  %       off1i = zeros(1,nsta);  %is that necessary to align the whole trace?
  
  %%%Align and compute the RCC based on the entire win, and take that as the input signal!
  optdat = [];  % win segment of interest
  ortdat = [];
  vertdat = [];
  optdat(:, 1) = optseg(1+msftaddm: end-msftaddm, 1); % time column
  ortdat(:, 1) = ortseg(1+msftaddm: end-msftaddm, 1);
  vertdat(:, 1) = vertseg(1+msftaddm: end-msftaddm, 1);
  for ista = 1: nsta
    optdat(:, ista+1) = optseg(1+msftaddm-off1i(k,ista): end-msftaddm-off1i(k,ista), ista+1);
    ortdat(:, ista+1) = ortseg(1+msftaddm-off1i(k,ista): end-msftaddm-off1i(k,ista), ista+1);
    vertdat(:, ista+1) = vertseg(1+msftaddm-off1i(k,ista): end-msftaddm-off1i(k,ista), ista+1);
  end
  
  %Align the noise using the same offset
  noidat = [];  % 4-s prior to signal win
  noidat(:, 1:2) = STAopt(max(floor((tstbuf-4)*sps+1),1): min(floor((tstbuf-0)*sps),86400*sps), 1:2); % sta 1
  noidat(:, 3) = STAopt(max(floor((tstbuf-4)*sps+1)-off1i(k,2),1): ...
    min(floor((tstbuf-0)*sps)-off1i(k,2),86400*sps), 3); % sta 2
  noidat(:, 4) = STAopt(max(floor((tstbuf-4)*sps+1)-off1i(k,3),1): ...
    min(floor((tstbuf-0)*sps)-off1i(k,3),86400*sps), 4); % sta 3
  
  %%%taper the signal and obtain the new rcc between tapered signals
  %%%2022/06/06, do NOT taper whatsoever!!
  sigsta = zeros(size(optdat,1), nsta);
  for ista = 1:nsta
    tmp = optdat(:,ista+1); %best aligned, filtered
    %detrend and taper only the data, NOT the noise
    tmp = detrend(tmp);
    %         ltmp = length(tmp);
    %         fractap = 2*sps/ltmp; % if fractap is >=1, n-point von Hann window is returned
    %         ptstap = fractap/2*size(tmp,1); % if fractap is >=1, n-point von Hann window is returned
    %         w = tukeywin(size(tmp,1),fractap);
    %         tmp = w.* tmp;
    %         tmp = detrend(tmp); %detrend again for caution
    sigsta(:,ista) = tmp;
  end
  %compute running CC between 3 stations
  [ircc,rcc12,rcc13,rcc23] = RunningCC3sta(sigsta,rccmwlen);
  ircc = ircc-overshoot;
  rcc = (rcc12+rcc13+rcc23)/3;
  rcc1i = zeros(length(rcc),nsta-3);
  for ista = 4:nsta
    [~,rcc1i(:,ista-3)] = RunningCC(sigsta(:,1), sigsta(:,ista), rccmwlen);
  end
  sigsta = detrend(sigsta(overshoot+1:end-overshoot, :));  %excluding the overshoot
  
  %if only use pairs 12 and 13
  rcc = (rcc12+rcc13)/2;
  rcccomb(:,2) = rcc;
  
  %%%for ort. comp
  sigstaort = zeros(size(ortdat,1), nsta);
  for ista = 1:nsta
    tmp = ortdat(:,ista+1); %best aligned, filtered
    tmp = detrend(tmp);
    sigstaort(:,ista) = tmp;
  end
  [irccort,rcc12,rcc13,rcc23] = RunningCC3sta(sigstaort,rccmwlen);
  irccort = irccort-overshoot;
  rccort = (rcc12+rcc13)/3;
  sigstaort = detrend(sigstaort(overshoot+1:end-overshoot, :));  %excluding the overshoot
  
  %%%for vertical comp
  sigstavert = zeros(size(vertdat,1), nsta);
  for ista = 1:nsta
    tmp = vertdat(:,ista+1); %best aligned, filtered
    tmp = detrend(tmp);
    sigstavert(:,ista) = tmp;
  end
  sigstavert = detrend(sigstavert(overshoot+1:end-overshoot, :));  %excluding the overshoot
  
  %plot the comparison between rcccat and rcc
  %       if ~noiseflag
  %         [f] = plt_rcccatvsrcc(irccran,ircccat,rcccat,nwin,sps,ircc,rcc);
  %       end
  %       median(rcccat)
  %       median(rcc)
  
  %       f1 = plt_traceindetail(sigsta,stas,25*sps);
  %       title(f1.ax(1),'Signal');
  
  %       keyboard
  
  %%%Doing this here is mainly to note down the thresholds used in the real data
  fixthresh = zeros(nsta,1);  %fixed thresholds used in real data to stop iteration
  mswccpksep = zeros(nsta,1); %Med sep of peaks in sig-wlet CC
  for ista = 1:nsta
    wlet = greenf(:,ista);  %template here is best aligned, tapered, linear trend removed, filtered
    lwlet = length(wlet);
    sig = sigsta(:,ista); %best aligned, filtered, tapered
    lsig = length(sig);
    dt = 1/sps;  % sampling interval
    twlet = zcrosses(ista)*dt;
    %get the master CC between sig and wlet, to know what is range of weighted master CC
    nfft = lsig;
    [coef, lag] = xcorr(sig, wlet, nfft, 'none'); % unnormalized master raw CC
    lrcc = length(rcccat); % length of running CC
    ldiff = nfft-lrcc;  % difference in length
    %effective raw CC that corresponds to the index of the overlapping portion between signal (or
    %sigdecon) and running CC, same length as rcc
    itwlet = round(twlet/dt); % time shift of main arrival of wavelet in samples
    coefeff = coef(lag+itwlet >= 1+round(ldiff/2) & ...
      lag+itwlet <= nfft-round(ldiff/2));
    %find all peaks in the effective master raw CC
    [pkhgt, pkind] = findpeaks(coefeff);
    mswccpksep(ista) = round(median(diff(pkind)));
    %rcc serves the weight as the peak height, aka the master raw cc value at the peak
    if ista <=3
      wtcoef = rcccat(pkind).* pkhgt;
    else
      wtcoef = rcc1icat(pkind, ista-3).* pkhgt;
    end
    fixthresh(ista) = median(wtcoef);  % median of the weighted master CC, could be percentile?
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %FLAG to simulate the behavior of noise
  %       noiseflag = 1;
  
  %       seedmat = randi(1000,200,1);
  %       for k = 1: length(seedmat)
  %       seed = seedmat(k);
  
  if noiseflag
    %obtain the amp and phase spectra of records via fft
    nfft = size(optseg,1); % number of points in fft
    [xf,ft,amp,pha] = fftspectrum(optseg(:,2:end), nfft, sps,'twosided');
    
    %uniform, random phase with the same span [-pi,pi];
    mpharan = minmax(pha');
    seed = idxbst(iii);
    %         seed = k;
    %         seed = 3;
    rng(seed);
    pharand = (rand(nfft,nsta)-0.5)*2*pi;  %make the phases span from -pi to pi
    
    %         figure
    %         for ista = 1: nsta
    %           subplot(1,nsta,ista)
    %           histogram(pharand(:,ista));
    %         end
    %         keyboard
    %
    %         xf2 = amp*nfft.*exp(1i*pha);
    %         opt2 = optseg;
    %         opt2(:,2:end) = real(ifft(xf2,nfft));
    %         figure
    %         plot(optseg(:,end),'k','linew',1); hold on
    %         plot(opt2(:,end),'r','linew',0.5);
    %         xlabel('Samples at 160 sps');
    %         ylabel('Amplitude');
    %         xlim([0 lsig]);
    %
    %         figure
    %         plot((1:nfft)/sps, optseg(:,end),'k','linew',1);
    %         xlabel('Time (s)');
    %         ylabel('Amplitude');
    %         xlim([0 lsig/sps]);
    %
    %         figure
    %         subplot(3,1,1);
    %         plot(ft,amp(:,end));
    %         xlabel('Frequency (Hz)');
    %         ylabel('Amplitude');
    %         subplot(3,1,2)
    %         plot(ft,pha(:,end));
    %         xlabel('Frequency (Hz)');
    %         ylabel('Phase (rad)');
    %         ylim([-4 4]);
    %         subplot(3,1,3)
    %         plot(ft,pharand(:,end));
    %         xlabel('Frequency (Hz)');
    %         ylabel('Randomized phase (rad)');
    %         ylim([-4 4]);
    
    %construct record with the same amplitude but random phase
    xfrand = amp.*nfft.*exp(1i.*pharand);
    optseg(:,2:end) = real(ifft(xfrand,nfft));
    
    %         figure
    %         plot((1:nfft)/sps, optseg(:,end),'k','linew',1);
    %         xlabel('Time (s)');
    %         ylabel('Amplitude');
    %         xlim([0 lsig/sps]);
    %         ylim([-0.3 0.3]);
    
    %         figure
    %         subplot(131)
    %         histogram(optseg(:,2)); [MUHAT,SIGMAHAT] = normfit(optseg(:,2))
    %         subplot(132)
    %         histogram(optseg(:,3)); [MUHAT,SIGMAHAT] = normfit(optseg(:,3))
    %         subplot(133)
    %         histogram(optseg(:,4)); [MUHAT,SIGMAHAT] = normfit(optseg(:,4))
    
    %%%for orthogonal components, with the same phase??
    [xf,ft,amp,pha] = fftspectrum(ortseg(:,2:end), nfft, sps,'twosided');
    %         pharand = (rand(nfft,3)-0.5)*2*pi;
    xfrand = amp.*nfft.*exp(1i.*pharand);
    ortseg(:,2:end) = real(ifft(xfrand,nfft));
    
    %%%for vertical components, with the SAME phase!!
    [xf,ft,amp,pha] = fftspectrum(vertseg(:,2:end), nfft, sps,'twosided');
    %         pharand = (rand(nfft,3)-0.5)*2*pi;
    xfrand = amp.*nfft.*exp(1i.*pharand);
    vertseg(:,2:end) = real(ifft(xfrand,nfft));
    
    off1iw = zeros(nwin,nsta);  % the best alignment between sta2, sta3 wrt sta1 for each subwin
    ccaliw = zeros(nwin,1+nsta-3);  % CC value using the best alignment, including 4th stas
    ircccat = [];   % concatenated indices of RCC
    irccran = zeros(nwin,2);  % start and end indices (range) of RCC of all subwins
    rcccat = [];  % average concatenated RCC
    rcc1icat = [];  % concatenated RCC between sta 1 and 4
    rccpaircat = [];  % concatenated RCC between each station pair, order is 12, 13, 23
    ccwpair = []; % 0-lag overall cc of each subwin, between each station pair, order is 12, 13, 23
    ccw1i = []; % same as above, but between sta 1 and 4
    
    for iwin = 1: nwin
      isubwst = windows(iwin,1)-max(floor(tstbuf*sps+1-overshoot-msftaddm),1);
      isubwed = windows(iwin,2)-max(floor(tstbuf*sps+1-overshoot-msftaddm),1);
      if iwin == 1
        isubwst = isubwst-overshoot;
      end
      if iwin == nwin
        isubwed = isubwed+overshoot;
      end
      
      %align records
      %         optcc = [];
      optcc = detrend(optseg(isubwst: isubwed, 2:end));
      %%%for each short win, allow the same shift for noise as for data
      %         if noiseflag
      %           msftadd = 50;
      %           loffmax = 20*sps/40;
      %         else
      msftadd = (round(max(abs([off12ran off13ran])))+1)*sps/40;  %+1 for safety
      loffmax = 4*sps/40;
      %         end
      ccmid = ceil(size(optcc,1)/2);
      ccwlen = round(size(optcc,1)-2*(msftadd+1));  % minus ensures successful shifting of records
      ccmin = 0.01;  % depending on the length of trace, cc could be very low
      iup = 1;    % times of upsampling
      
      %%%For short win, do we align noise *exactly* the same way as data?
      %ask how well can you align the noise, does not matter if the location is
      %the most reliable, but the resulting CC is the highest possible
      [off12con,off13con,ccaliw(iwin,1)] = constrained_cc_loose(optcc(:,1:3)',ccmid,...
        ccwlen,msftadd,ccmin,iup);
      % if a better alignment cannot be achieved, use 0,0
      if off12con == msftadd+1 && off13con == msftadd+1
        off12con = 0;
        off13con = 0;
        fprintf('Tremor burst %d cannot be properly aligned, double-check needed \n',k);
      end
      off1iw(iwin,1) = 0;
      off1iw(iwin,2) = round(off12con);
      off1iw(iwin,3) = round(off13con);
      for ista = 4: nsta
        [ccaliw(iwin,ista-2),off1iw(iwin,ista)] = xcorrmax(optcc(:,1), optcc(:,ista), 1.5*msftadd, 'coeff');
      end
      %Align records
      optdat = [];  % win +/-3 s, segment of interest,first 1s will be tapered
      ortdat = [];
      optdat(:, 1) = optseg(isubwst: isubwed, 1); % time column
      ortdat(:, 1) = ortseg(isubwst: isubwed, 1);
      for ista = 1: nsta
        optdat(:, ista+1) = optseg(isubwst-off1iw(iwin,ista): isubwed-off1iw(iwin,ista), ista+1);
        ortdat(:, ista+1) = ortseg(isubwst-off1iw(iwin,ista): isubwed-off1iw(iwin,ista), ista+1);
      end
      
      subw = zeros(size(optdat,1), nsta);
      for ista = 1:nsta
        tmp = optdat(:,ista+1); %best aligned, filtered
        %detrend and taper only the data, NOT the noise
        tmp = detrend(tmp);
        %%%2022/06/06, do NOT taper whatsoever!!
        %           %taper the very start or end subwin and obtain the new rcc between tapered signals
        %           ltmp = length(tmp);
        %           fractap = sps/ltmp; % if fractap is >=1, n-point von Hann window is returned
        %           ptstap = fractap*ltmp; % if fractap is >=1, n-point von Hann window is returned
        %           if iwin==1
        %             w = coswinlh(ltmp,fractap);
        %             tmp = w.* tmp;
        %             tmp = detrend(tmp); %detrend again for caution
        %           elseif iwin==nwin
        %             w = coswinrh(ltmp,fractap);
        %             tmp = w.* tmp;
        %             tmp = detrend(tmp); %detrend again for caution
        %           end
        subw(:,ista) = tmp;
      end
      %compute running CC between 3 stations
      [irccw,rccw12,rccw13,rccw23] = RunningCC3sta(subw,rccmwlen);
      rccw = (rccw12+rccw13+rccw23)/3;
      irccw = irccw + windows(iwin,1) - windows(1,1);   %convert to global index
      if iwin == 1
        irccw = irccw - overshoot;
      end
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
    irccrank{k} = irccran;
    windowsk{k} = windows-windows(1,1)+1;
    off1iwk{k} = off1iw;
    ccaliwk{k} = ccaliw;
    ccwpairk{k} = ccwpair;
    mrccwpairk{k} = mrccwpair;
    % keyboard
    %         %if only use the mean RCC from the 2 pairs that have the highest overall CC
    %         [~,ind] = min(sum(ccwpair,1));
    %         rcccat = mean(rccpaircat(:,setdiff(1:3,ind)), 2);
    
    %if only use the mean RCC from pair 12 and 13
    rcccat = mean(rccpaircat(:,[1 2]), 2);
    
    %if choose NOT to use RCC weighting; for easier comparison
    if ~rccflag
      rcccat = ones(size(rcccat));
      rcc1icat = ones(size(rcc1icat));
    else
      if whichrcc == 1
        rcccat = rccpaircat(:,1);
      elseif whichrcc == 2
        rcccat = rccpaircat(:,2);
      elseif whichrcc == 3
        rcccat = rccpaircat(:,3);
      end
    end
    
    rcccomb(:,3) = rcccat;
    
    %         figure
    %         subplot(231); hold on; ax=gca;
    %         histogram(rccpaircat(:,1)); title('rcc12'); [MUHAT,SIGMAHAT] = normfit(rccpaircat(:,1));
    %         plot([MUHAT MUHAT],ax.YLim,'r--');
    %         subplot(232); hold on; ax=gca;
    %         histogram(rccpaircat(:,2)); title('rcc13'); [MUHAT,SIGMAHAT] = normfit(rccpaircat(:,2));
    %         plot([MUHAT MUHAT],ax.YLim,'r--');
    %         subplot(233); hold on; ax=gca;
    %         histogram(rccpaircat(:,3)); title('rcc23'); [MUHAT,SIGMAHAT] = normfit(rccpaircat(:,3));
    %         plot([MUHAT MUHAT],ax.YLim,'r--');
    %         subplot(234); hold on; ax=gca;
    %         histogram(mean(rccpaircat(:,[1 2]),2)); title('rcc12+rcc13'); [MUHAT,SIGMAHAT] = normfit(sum(rccpaircat(:,[1 2]),2));
    %         plot([MUHAT MUHAT],ax.YLim,'r--');
    %         subplot(235); hold on; ax=gca;
    %         histogram(mean(rccpaircat(:,[2 3]),2)); title('rcc13+rcc23'); [MUHAT,SIGMAHAT] = normfit(sum(rccpaircat(:,[2 3]),2));
    %         plot([MUHAT MUHAT],ax.YLim,'r--');
    %         subplot(236); hold on; ax=gca;
    %         histogram(mean(rccpaircat(:,[1 3]),2)); title('rcc12+rcc23'); [MUHAT,SIGMAHAT] = normfit(sum(rccpaircat(:,[1 3]),2));
    %         plot([MUHAT MUHAT],ax.YLim,'r--');
    
    %%%Is it true that the coherence between 2-3 is the highest among 3 pairs?
    %%%---Yes, for the concatenated rcc
    %       [f] = plt_rcccat(rccpaircat,sps);
    
    %%%obtain a single best alignment based on the entire win
    %       optcc = optseg(:, 2:end);
    optcc = detrend(optseg(1+msftaddm: end-msftaddm, 2:end));
    %         msftadd = 0.5*sps;
    ccmid = ceil(size(optcc,1)/2);
    ccwlen = round(size(optcc,1)-2*(msftadd+1));
    ccmin = 0.01;  % depending on the length of trace, cc could be very low
    iup = 1;    % times of upsampling
    %for the whole win, ask how well can you align the noise, does not matter if the location is
    %the most reliable, but the resulting CC is the highest possible
    [off12con,off13con,ccali(k)] = constrained_cc_loose(optcc(:,1:3)',ccmid,...
      ccwlen,msftadd,ccmin,iup);
    % if a better alignment cannot be achieved, use 0,0
    if off12con == msftadd+1 && off13con == msftadd+1
      off12con = 0;
      off13con = 0;
      fprintf('Tremor burst %d cannot be properly aligned, double-check needed \n',k);
    end
    off1ic(k,1) = 0;
    off1ic(k,2) = round(off12con);
    off1ic(k,3) = round(off13con);
    
    %%%obtain a single 'observational' alignment between 4th and 1st station, note that this might
    %%%be very different from the empirical prediction from 'empioffset4thsta002'
    %%%--there are also different options here, using sig, or sig-wlet CC as in decon routine
    mcoef = zeros(nsta-3, 3);
    mlag = zeros(nsta-3, 3);
    envrat = zeros(nsta-3, 3);
    for ista = 4: nsta
      [mcoef(ista-3, 1),off1ic(k,ista)] = xcorrmax(optcc(:,1), optcc(:,ista), msftadd, 'coeff');
      mlag(ista-3, 1) = off1ic(k,ista);
      envrat(ista-3, 1) = median(envelope(optcc(:,ista)))./median(envelope(optcc(:,1)));
      %do an overall CC between 4th and 2nd/3rd stas, to see which one they are most coehrent with
      for jjj = 2:3
        [mcoef(ista-3, jjj),mlag(ista-3, jjj)] = xcorrmax(optcc(:,1), optcc(:,ista), msftadd, 'coeff');
        envrat(ista-3, jjj) = median(envelope(optcc(:,ista)))./median(envelope(optcc(:,jjj)));
      end
      %empirical prediction from plane fitting in 'empioffset4thsta002'
      off14pred(k,ista-3) = round(off14mod(ista-3,1).*off1ic(k,2) + off14mod(ista-3,2).*off1ic(k,3) + ...
        off14mod(ista-3,3));
    end
    
    %for synthetic noise case, do NOT align traces
    off1i(k,1:3) = [0 0 0];  % technically unnecessary given it has been predefined
    
    %for synthetic noise case, do NOT align traces
    %Choice to make upon the actual-used alignment at 4th stations
    %         align14flag = 0;
    if align14flag
      off1i(k,4:end) = off1ic(k,4:end); %if you trust overall alignment at 4th stations
    else
      off1i(k,4:end) = zeros(1,nsta-3); %if you don't
    end
    
    %%%Align and compute the RCC based on the entire win, and take that as the input signal!
    optdat = [];  % win segment of interest
    ortdat = [];
    vertdat = [];
    optdat(:, 1) = optseg(1+msftaddm: end-msftaddm, 1); % time column
    ortdat(:, 1) = ortseg(1+msftaddm: end-msftaddm, 1);
    vertdat(:, 1) = vertseg(1+msftaddm: end-msftaddm, 1);
    for ista = 1: nsta
      optdat(:, ista+1) = optseg(1+msftaddm-off1i(k,ista): end-msftaddm-off1i(k,ista), ista+1);
      ortdat(:, ista+1) = ortseg(1+msftaddm-off1i(k,ista): end-msftaddm-off1i(k,ista), ista+1);
      vertdat(:, ista+1) = vertseg(1+msftaddm-off1i(k,ista): end-msftaddm-off1i(k,ista), ista+1);
    end
    
    %%%taper the signal and obtain the new rcc between tapered signals
    %%%2022/06/06, do NOT taper whatsoever!!
    sigsta = zeros(size(optdat,1), nsta);
    for ista = 1:nsta
      tmp = optdat(:,ista+1); %best aligned, filtered
      %detrend and taper only the data, NOT the noise
      tmp = detrend(tmp);
      sigsta(:,ista) = tmp;
    end
    %compute running CC between 3 stations
    [ircc,rcc12,rcc13,rcc23] = RunningCC3sta(sigsta,rccmwlen);
    ircc = ircc-overshoot;
    rcc = (rcc12+rcc13+rcc23)/3;
    rcc1i = zeros(length(rcc),nsta-3);
    for ista = 4:nsta
      [~,rcc1i(:,ista-3)] = RunningCC(sigsta(:,1), sigsta(:,ista), rccmwlen);
    end
    sigsta = detrend(sigsta(overshoot+1:end-overshoot, :));  %excluding the overshoot
    
    rcc = (rcc12+rcc13)/2;
    rcccomb(:,4) = rcc;
    
    %%%for ort. comp
    sigstaort = zeros(size(ortdat,1), nsta);
    for ista = 1:nsta
      tmp = ortdat(:,ista+1); %best aligned, filtered
      tmp = detrend(tmp);
      sigstaort(:,ista) = tmp;
    end
    [irccort,rcc12,rcc13,rcc23] = RunningCC3sta(sigstaort,rccmwlen);
    irccort = irccort-overshoot;
    rccort = (rcc12+rcc13)/3;
    sigstaort = detrend(sigstaort(overshoot+1:end-overshoot, :));  %excluding the overshoot
    
    %%%for vertical comp
    sigstavert = zeros(size(vertdat,1), nsta);
    for ista = 1:nsta
      tmp = vertdat(:,ista+1); %best aligned, filtered
      tmp = detrend(tmp);
      sigstavert(:,ista) = tmp;
    end
    sigstavert = detrend(sigstavert(overshoot+1:end-overshoot, :));  %excluding the overshoot
    
    %       figure
    %       hold on
    %       box on
    %       plot(ircc,rcc,'k-')
    %       plot(ircccat,rcccat,'r-');
    %       ylim([-1 1]);
    
  end
  
  %% load LFE catalogs
  % savefile = 'deconv_stats4th_allbstsig.mat';
  savefile = 'deconv_stats4th_no23_allbstsig.mat';
  load(strcat(rstpath, '/MAPS/',savefile));
  % savefile = 'deconv_stats4th_allbstnoi.mat';
  savefile = 'deconv_stats4th_no23_allbstnoi.mat';
  load(strcat(rstpath, '/MAPS/',savefile));
  
  %%
  %%%param for secondary sources removed
  locxyprojall = allbstsig.locxyprojall;
  tarvlsplstall = allbstsig.impindepall(:,1);
  nsrc = allbstsig.nsrc;
  imp = allbstsig.impindepall;
  locxyprojalln = allbstnoi.locxyprojall;
  tarvlsplstalln = allbstnoi.impindepall(:,1);
  nsrcn = allbstnoi.nsrc;
  impn = allbstnoi.impindepall;
  supertstr = 'Secondary sources removed';
  fnsuffix = [];
  
  impuse = imp;
  nsrcuse = nsrc;
  fnsuffix2 = [];
  
  %% signal + zoom-in + map locations + some reference symbols
  if ~noiseflag
    % %%%plot all short windows, and LFE detections
    % htin = 0.6*(nwin+1);   % maximum height allowed is 11 inches
    % widin = 8.4;  % maximum width allowed is 8.5 inches
    % nrow = 1;
    % ncol = 1;
    % f = initfig(widin,htin,nrow,ncol); %initialize fig
    % figxran = [0.08 0.98]; figyran = [0.08 0.98];
    % figxsep = 0.05; figysep = 0.02;
    % optaxpos(f,nrow,ncol,figxran,figyran,figxsep,figysep);
    % ax=f.ax(1); hold(ax,'on'); ax.Box = 'off'; grid(ax,'off');
    % %       set(ax, 'color', 'none');
    % %   xlim(ax,[]);
    % %   nolabels(ax,2);
    % longticks(ax,4);
    % axis(ax,'tight');
    % ylabel(f.ax(end),sprintf('Amplitude'),'FontSize',10);
    % xlabel(f.ax(end),sprintf('Time (s) on %s %s %s',dy,mo,yr),'FontSize',10);
    % vsep = 1;
    % ylim(ax,[-vsep*(nwin+1) -0.5*vsep]);
    
    % [clusibst,nclusibst,tclusibst]=clusters_in_burst(idxbst,impuse,nsrcuse);
    
    % cumutime=0;
    % for iwin = 1:nwin
    %   subwtmp = subwsave{iwin};
    %   plot(ax,(1:length(subwtmp))/sps,subwtmp(:,1)-vsep*iwin,'r','linew',1);
    %   plot(ax,(1:length(subwtmp))/sps,subwtmp(:,2)-vsep*iwin,'b','linew',1);
    %   plot(ax,(1:length(subwtmp))/sps,subwtmp(:,3)-vsep*iwin,'k','linew',1);
    %   %%%plot the deconvolved sources
    %   yloc = -vsep*(iwin+0.45);
    %   tst = cumutime;
    %   ted = cumutime+length(subwtmp)/sps;
    %   %           ind = find(impindepst(:,1)/sps >= tst & impindepst(:,1)/sps < ted);
    %   %           scatter(ax,impindepst(ind,1)/sps-tst,yloc*ones(size(impindepst(ind,1))),...
    %   %             6,[.6 .6 .6],'filled');
    %   ind = find(tclusibst(:,1)>=tst & tclusibst(:,2)<ted);
    %   clusibstxzoom = clusibst(ind);
    %   tclusibstxzoom = tclusibst(ind);
    %   nclus = length(ind);
    %   for iclus=1:2:nclus
    %     impclus=clusibstxzoom{iclus};
    %     scatter(ax,impclus(:,1)/sps-tst,yloc*ones(size(impclus(:,1))),10,'^','k');
    %   end
    %   for iclus=2:2:nclus
    %     impclus=clusibstxzoom{iclus};
    %     scatter(ax,impclus(:,1)/sps-tst,yloc*ones(size(impclus(:,1))),10,'v','k');
    %   end
    %   impuniclus = unique(cat(1,clusibstxzoom{:}),'rows');
    %   ist = sum(nsrcuse(1:idxbst-1))+1;
    %   ied = ist+nsrcuse(idxbst)-1;
    %   impindepst = impuse(ist:ied,:);
    %   ind = find(impindepst(:,1)/sps >= tst & impindepst(:,1)/sps <ted);
    %   impiso = setdiff(impindepst(ind,:),impuniclus,'rows');
    %   scatter(ax,impiso(:,1)/sps-tst,yloc*ones(size(impiso(:,1))),10,'k');
      
    %   cumutime = ted;
    % end
    % keyboard
    % fname = sprintf('seisbst%saligned.pdf',num2zeropadstr(idxbst, 3));
    % print(f.fig,'-dpdf',...
    %   strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/supfigures/',fname));
    %         keyboard

    %%
%     %%%these are outstanding migrations
%     if idxbst == 99
%       xzoom = [0 25]+52;
%     elseif idxbst == 100
%       xzoom = [0 25]+210;
%     elseif idxbst == 115
%       xzoom = [0 25]+129;
%     elseif idxbst == 164
%       xzoom = [0 25]+146;   % or 146, 49
%     elseif idxbst == 170
%       xzoom = [0 25]+27;
%     elseif idxbst == 174
%       xzoom = [0 25]+123;  % or 123, 75
%     elseif idxbst == 181
%       xzoom = [0 25]+5;
%     end
    
    %%%these are chosen according to amp binning from small to large
    %bin 1: 079, 094, 011, 013
    %bin 2: 085, 174,
    %bin 3: 049, 174, 057
    %bin 4: 164, 067, 158
    %bin 5: 060, 181, 193
    if idxbst == 79
      xzoom = [0 17];
      panelnm = 'b';
      yran2l = [-0.1 0.1];
    % elseif idxbst == 94
    %   xzoom = [0 25]+130;
    %   panelnm = 'b';
    %   yran2l = [-0.4 0.4];
    elseif idxbst == 11
      xzoom = [0 25]+66;
      panelnm = 'c';
      yran2l = [-0.2 0.2];
    % elseif idxbst == 13
    %   xzoom = [0 25]+20;
    %   panelnm = 'b';
    %   yran2l = [-0.2 0.2];
    elseif idxbst == 85
      xzoom = [0 25]+25;
      panelnm = 'd';
      yran2l = [-0.3 0.3];
    elseif idxbst == 174
      xzoom = [0 25]+75;   % or 123, 75
      panelnm = 'e';
      yran2l = [-0.5 0.5];
    % elseif idxbst == 49
    %   xzoom = [0 25]+119;
    %   panelnm = 'd';
    %   yran2l = [-0.5 0.5];
    % elseif idxbst == 57
    %   xzoom = [0 25]+225;
    %   panelnm = 'd';
    %   yran2l = [-0.4 0.4];
    % elseif idxbst == 164
    %   xzoom = [0 25]+49;   % or 146, 49
    %   panelnm = 'd';
    %   yran2l = [-0.8 0.8];
    elseif idxbst == 60
      xzoom = [0 25]+25;
      panelnm = 'f';
      yran2l = [-1.2 1.2];
    end
    %%
    tbostplt=tbosti(indbi); %bostock's LFE detections inside the burst, pointing to positive peaks
    tbostzcplt=tbostplt-(ppeaks(1)-zcrosses(1))/sps; %shift it to zero-crossings
    winnoi = [1 3]; %empirical window for pre-arrival noise
    winsig = [4 7]; %empirical window for main arrival
    patchwins = [winnoi;winsig];
    ist = sum(nsrcuse(1:idxbst-1))+1;
    ied = ist+nsrcuse(idxbst)-1;
    impindepst = impuse(ist:ied,:);
    if idxbst == 79
      f=plt_002lfe_signalzoomv2bst79(greenf,sigsta,[],impindepst,sps,xzoom,...
        tstbuf,dy,mo,yr,tbostzcplt,idxbst,yran2l,saveflag);
      text(f.ax(2),0.01,0.85,panelnm,'FontSize',11,'unit','normalized','EdgeColor','k','Margin',1);
    else
      f=plt_002lfe_signalzoomv2(greenf,sigsta,[],impindepst,sps,xzoom,...
        tstbuf,dy,mo,yr,tbostzcplt,idxbst,yran2l,saveflag);
      text(f.ax(2),0.01,0.85,panelnm,'FontSize',11,'unit','normalized','EdgeColor','k','Margin',1);
      delete(f.ax(1));
    end
%     f.ax(2).YAxis(1).Limits = yran;
    keyboard
%     if saveflag
      % orient(f.fig,'landscape');
      fname = sprintf('lfetempvstremorv2_%s_%d-%ds.pdf',num2zeropadstr(idxbst, 3),...
        xzoom(1),xzoom(2));
      print(f.fig,'-dpdf',...
        strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
%     end
    
    % %         [f] = plt_agu2022abstractv3(greenf(:,1:3),sigsta(:,1:3),impindepst,sps,xzoom,off1iw,loff_max,...
    % %           tstbuf,dy,mo,yr,ftrans,'spl');
  end
  
end

