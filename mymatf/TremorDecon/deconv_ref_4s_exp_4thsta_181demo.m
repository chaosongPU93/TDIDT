% deconv_ref_4s_exp_4thsta_181demo.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on 'deconvbursts002_ref_4s_exp' and 'deconv_ref_4s_exp_rand_fn', but
% this version tries to involve the 4th station to check the deconvolution
% result from the trio stations
%
% Directly from 'deconv_ref_4s_exp_4thsta_fn.m'. Serve as the script to run 
% all demonstrative figures needed in the paper using burst #181 as an 
% example. This makes it easy to run 'deconv_ref_4s_exp_4thsta_fn.m' solely
% as a function to compute the statistics, or summary figures that are needed
% for every burst.
%
% --The key of making this work is to obtain the observational/empirical/
%   theoretical time offset between the 4th sta and 1st sta (or equivalently
%   the aligned 2nd or 3rd stas), which has already been estimated from the 
%   script 'empioffset4thsta002', for each deconvolved source location.
% --As long as you know the time offset off14 for each source, you would know
%   the arrival time at sta 4 based on the arrival time at sta 1 for all
%   sources.
% --Convolving the arrival times of sources at 4th sta and its own templates
%   at the optimal component results in the prediction at 4th sta. The residual
%   of subtracting pred from the opt data, if this is smaller than the data,
%   it means the deconvolved sources are explaining the opt data at 4th sta
%   as well.
% --2022/11/21, add the 'noise' option similar to that in 'deconv_ref_4s_exp_rand_fn'
% --2022/11/30, to resolve the amp issue at KLNB wrt. PGC, I modified 'fact' in 
%   'rd_daily_bpdata', I had to recompute templates in 'lfetemp002_160sps.m',
%   and then all results were recalculated!
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/03/06
% Last modified date:   2024/03/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e   % Set the format to 5-digit floating point
clear
clc
close all

%% for easy testing
defval('idxbst',181); %global indices of bursts to run 
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
[f,snrf,snrfort,snrfvert]=plt_templates_bp_bysta(greenf,greenfort,greenfvert,stas,staplt,lowlet,hiwlet,sps,...
  [winnoi;winsig]);


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
      
      %%  
      %%%finalize the signal, noise, and template (Green's function)
      sigdecon = [];
      pred = [];
      ampit = [];
      for ista = 1:3
        wlet = greenf(:,ista);  %template here is best aligned, tapered, linear trend removed, filtered
        lwlet = length(wlet);
        sig = sigsta(:,ista); %best aligned, filtered, tapered
        lsig = length(sig);
        noi = noidat(:,ista+1);
        noi = detrend(noi); %best aligned, filtered, linear trend removed
        
        dt = 1/sps;  % sampling interval
        twlet = zcrosses(ista)*dt;
        width = 2.5;  % width for Gaussian filter
        dres_min = 0.5;  % tolerance, percentage change in residual per iteration, in terms of variance reduction
        mfit_min = 5e-1;  % tolerance, norm of the residual/misfit
        tdura = 0.5;  % estimate from the broadband template from fam 002
        nit_max = round(1.5*1/tdura*(tlenbuf));  % max numer of iterations
        nimp_max = round(1/tdura*(tlenbuf));%a single peak is ~20 samples wide; maybe a little less (at 100 sps). ~0.4s, 1/0.4=2.5
        fpltit = 0;  % plot flag for each iteration
        fpltend = 0;  % plot flag for the final iteration
        fpltchk = 0; % plot flag for intermediate computations
        
        if noiseflag
          %%%As of 2022/09/27, 'threshold' would be auto-computed based on data then feed to the
          %%%decon if the 'noiseflag' is on
          [sigdecon(:,ista),pred(:,ista),res,dresit,mfitit{ista},ccchgit,ampit{ista},nit,fighdl] = ...
            iterdecon(sig,wlet,rcccat,noi,fixthresh(ista),dt,twlet,width,dres_min,...
            mfit_min,nit_max,nimp_max,fpltit,fpltend,fpltchk);
        else
          [sigdecon(:,ista),pred(:,ista),res,dresit,mfitit{ista},ccchgit,ampit{ista},nit,fighdl] = ...
            iterdecon(sig,wlet,rcccat,noi,[],dt,twlet,width,dres_min,mfit_min,nit_max,nimp_max,...
            fpltit,fpltend,fpltchk);
%           [sigdecon(:,ista),pred(:,ista),res,dresit,mfitit,ampit{ista},nit,fighdl] = ...
%             iterdecon_chron(sig,wlet,rcccat,noi,[],dt,twlet,width,dres_min,mfit_min,nit_max,nimp_max,...
%             fpltit,fpltend,fpltchk);
%         [sigdecon(:,ista),pred(:,ista),res,dresit,mfitit,ampit,nit,fighdl] = ...
%           iterdecon_rcc(sig,wlet,rcc,medrcc,dt,twlet,width,dres_min,mfit_min,nit_max,nimp_max,fpltit,fpltend,fcheck);

          %%% do you want to enforce the chronological order?
          chronoflag = 0;
          if chronoflag
            x = ampit{ista};
%             [~,ikeep,idiscard] = ascendvector(x(:,1));
            [~,ikeep,idiscard] = ascendvector(x(:,1),x(:,2),prctile(x(:,2),50));
            
%             % %%
%             figure;
%             subplot(221);
%             scatter(1:size(x,1),x(:,1),15,log10(x(:,2)),'filled');
%             c=colorbar;
%             colormap('jet');
%             c.Label.String = 'log_{10}{Amp.}';
%             xlabel('Iteration #');
%             ylabel('Decon index');
%             title('Original');
%             
%             subplot(222);
%             scatter(1:size(x(ikeep,:),1),x(ikeep,1),15,log10(x(ikeep,2)),'filled');
%             c=colorbar;
%             colormap('jet');
%             c.Label.String = 'log_{10}{Amp.}';
%             xlabel('Iteration #');
%             ylabel('Decon index');
%             title('Force to maintain chrono. order');
%             
%             subplot(223);
%             histogram(log10(x(:,2)),'FaceColor',[.5 .5 .5],'BinWidth',0.05); hold on
%             histogram(log10(x(ikeep,2)),'FaceColor','r','BinWidth',0.05);
%             histogram(log10(x(idiscard,2)),'FaceColor','b','BinWidth',0.05);
%             ax=gca;
%             plot([log10(prctile(x(:,2),50)) log10(prctile(x(:,2),50))],ax.YLim,'k--','linew',1);
%             xlabel('log_{10}{Amp.}');
%             ylabel('Count');
%           
            sigdecon(:,ista) = zeros(lsig,1);
            for jj = 1:length(ikeep)
              sigdecon(x(ikeep(jj),1),ista) = sigdecon(x(ikeep(jj),1),ista)+x(ikeep(jj),2);
            end
            ampit{ista} = x(ikeep,:);
            dresit = dresit(ikeep);
            mfitit = mfitit(ikeep);
          end
        end
       
        if fpltend
          ax = fighdl{2}.ax(1);
          hold(ax,'on');
          text(ax,0.05,0.85,stas(ista,:),'unit','normalized');
          hold(ax,'off');
        end
        
        nit
        
      end

      detecttype = [];
      xzoom = [5 30];

      if ~noiseflag
        ista = 1;
        %%%orthogonal comp.
        mfitort = msfitvsiter(ampit,sigstaort,greenfort,zcrosses,ista,sps);
        %%%vertical comp.
        mfitvert = msfitvsiter(ampit,sigstavert,greenfvert,zcrosses,ista,sps);
        %%%optimal comp.
        varsig = var(sigsta(:,ista));
        mfit = mfitit{ista};
        varsigort = var(sigstaort(:,ista));
        varsigvert = var(sigstavert(:,ista));
        [f1,f2] = plt_msfitvsiter_allcomp(mfit,varsig,mfitort,varsigort,...
          mfitvert,varsigvert,stas,ista,detecttype,0);
      end

%       keyboard
%       %%%For paper demo, plot the impulses, signal, prediction and res right after deconvolution
%       if ~noiseflag
%         ista = 1;
        % saveflag = 0; %don't save figures
        % xzoom = [5 30]; 
%         f = plt_impsigpredres_1sta(sigdecon,sigsta,pred,stas,ista,sps,xzoom,detecttype,saveflag);
%       end
      
      %% Group nearest impulses from different stations into pairs, using moving searching range
      spsscale = sps/40;
      loff_max = 6*spsscale;  % maximum allowed shift is 4 samples at 40 Hz, it is also true from the synthetic distribution
      %note the output 'impindep' gives the arrival index of impulse at each station, after
      %alignment based upon the entire window 'off1i', and the last three cols are the arrival time
      %difference, NOT the true location yet! 
      refsta = 1;
%       [impindep,imppairf,indpair,sharp] = groupimptripdecon_ref(sigdecon,ampit,irccran,rcccat,...
%         off1i(k,1:3),off1iw(:,1:3),loff_max,refsta);
      [impindep,imppairf,indpair,sharp] = groupimptripdecon_ref_nooff23(sigdecon,ampit,irccran,rcccat,...
        off1i(k,1:3),off1iw(:,1:3),loff_max,refsta);
      
%       %plot the sharpness of grouped peaks in res-wlet CC 
%       f=plt_srcsharpness(sharp);
%       
%       nsrc(k) = size(impindep,1);
%       end
%     
%     figure
%     histogram(nsrc);
%     text(0.05,0.9,sprintf('med=%d',median(nsrc)),'Units','normalized');
%     xlabel('number of sources');
%     ylabel('counts');
%     
      
%       if ~noiseflag
%         %%%plot the individually deconvolved impulses and the grouping result
%         [f] = plt_groupimptriplets(sigdecon,impindep,stas,ircccat,rcccat,sps,xzoom,detecttyp,saveflag);
%       end
      
      %note here 'impindepst' inherits the first 6 cols from 'impindep', but the last three cols 
      %are adjusted from arrival time difference to the true location offset accounting for the best
      %alignment upon each subwin that is also used in grouping!
      impindep(:,7:8) = impindep(:,7:8)+repmat([off1i(k,2) off1i(k,3)],size(impindep,1),1); %account for prealignment
      impindepst = sortrows(impindep,1);
      
      if ~noiseflag
        %%%prediction via convolution between grouped impulses and template at each station
        [~,predgrp,resgrp,predgrpl,resgrpl,l2normred1,varred]= ...
          predsig_conv_imptemp(sigsta(:,1:3),optdat(:,1:4),impindepst,...
          greenf,zcrosses,overshoot,stas,0,sps);
        %%%simply plot the trio stations signal, prediction and residual
        pltsta = [1 2 3];
%         f1 = plt_selsigpredres(sigsta,predgrp,resgrp,varred1,stas,pltsta,sps,...
%           xzoom,detecttype,saveflag);
        f1 = plt_selsigpred(sigsta,predgrp,varred,stas,pltsta,sps,...
          xzoom,detecttype,saveflag);
        %%
        %%%plot the predictions at ort and vert. comp. together in same plot,
        %%%as well as the templates
        [~,predgrport,resgrport,~,~,l2normred1ort,varredort]= ...
          predsig_conv_imptemp(sigstaort(:,1:3),ortdat(:,1:4),impindepst,...
          greenfort,zcrosses,overshoot,stas,0,sps);
        %%%prediction at vert. comp.
        [~,predgrpvert,resgrpvert,~,~,l2normred1vert,varredvert]= ...
          predsig_conv_imptemp(sigstavert(:,1:3),vertdat(:,1:4),impindepst,...
          greenfvert,zcrosses,overshoot,stas,0,sps);
        winnoi = [1 3]; %empirical window for pre-arrival noise
        winsig = [4 7]; %empirical window for main arrival    
        patchwins = [winnoi;winsig];
        %%
        normflag=1;
        f = plt_selsigpred_allcomp(sigsta,predgrp,varred,sigstaort,predgrport,...
          varredort,sigstavert,predgrpvert,varredvert,greenfort,greenfvert,...
          stas,pltsta,sps,xzoom,patchwins,normflag,detecttype,saveflag);

        f = plt_selsigpred_ortvert(sigstaort,predgrport,varredort,...
          sigstavert,predgrpvert,varredvert,greenfort,greenfvert,...
          stas,pltsta,sps,xzoom,patchwins,detecttype,saveflag);

      end
%   keyboard

      %%
%       %%%plot the scatter of offsets, accounting for prealignment offset, == true offset
%       span = max(range(off1iw(:,2))+2*loff_max, range(off1iw(:,3))+2*loff_max);
%       xran = [round(mean(minmax(off1iw(:,2)'))-span/2)-1, round(mean(minmax(off1iw(:,2)'))+span/2)+1];
%       yran = [round(mean(minmax(off1iw(:,3)'))-span/2)-1, round(mean(minmax(off1iw(:,3)'))+span/2)+1];
%       cran = [0 lsig];
%       f1.fig = figure;
%       f1.fig.Renderer = 'painters';
%       ax1=gca;
%       [ax1,torispl,mamp] = plt_decon_imp_scatter_ref(ax1,impindepst,xran,yran,cran,off1iw,loff_max,...
%         sps,50,'mean','tori','comb');
%       scatter(ax1,off1i(k,2),off1i(k,3),20,'ks','filled','MarkerEdgeColor','k');
%       title(ax1,'Independent, grouped');
%       keyboard
      
%       %%%plot the scatter of offsets, shifted to the same origin, the best alignment
%       xran = [-loff_max+off1i(k,2)-1 loff_max+off1i(k,2)+1];
%       yran = [-loff_max+off1i(k,3)-1 loff_max+off1i(k,3)+1];
%       f1.fig = figure;
%       f1.fig.Renderer = 'painters';
%       ax1=gca;
%       ax1 = plt_decon_imp_scatter_ref_sft(ax1,impindepst,xran,yran,cran,off1i(k,:),off1iw,...
%         loff_max,irccran,sps,50,'mean','tori');
%       scatter(ax1,off1i(k,2),off1i(k,3),20,'ks','filled','MarkerEdgeColor','k');
%       title(ax1,'Independent, grouped, shifted');
      

      %% plot amp ratio 12 and 13, and 23, right after grouping
%       %convert time offset to relative loc
%       imploc = off2space002(impindepst(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
%       mapview = 'offset';
%       [f,srcamprat,psrcampratscl,nsrcampratscl] = plt_srcampratio(imploc,impindepst,greenf,sps,mapview);
%       figure; 
%       subplot(211); hold on; ax=gca;
%       histogram(log10(srcamprat(:,3))); 
%       msrcampr(k) = median(log10(srcamprat(:,3)));
%       plot([msrcampr(k) msrcampr(k)],ax.YLim,'r--');      
%       xlabel('log_{10}{Amp ratio 2/3}');
%       ylabel('Counts');
%       subplot(212); hold on; ax=gca;
%       histogram(log10(psrcampratscl(:,3))); 
%       mpsrcamprs(k) = median(log10(psrcampratscl(:,3)));
%       plot([mpsrcamprs(k) mpsrcamprs(k)],ax.YLim,'r--');      
%       xlabel('log_{10}{Positive scaled amp ratio 2/3}');
%       ylabel('Counts');

      %% with 2ndary srcs, prediction of impulse tarvl at 4th sta given sources and empirical off14-src relation
%       %%%predict the tarvl and amp at 4th sta 
%       for ista = 4:nsta
%         impindepst(:,9+(ista-4)*2+1) = pred_tarvl_at4thsta(stas(ista,:),impindepst(:,7),...
%           impindepst(:,8),impindepst(:,1),off1i(k,ista));
%         %%%there is decision making here about what impulse amp should be given to 4th sta
%         %%%1st way is to use the amp SAME as of the sta that is most similar in waveform to 4th
%         [~,ind] = max(mcoef(ista-3,:)); % which of the trio stas is most similar to 4th sta?
%         impindepst(:,9+(ista-3)*2) = impindepst(:,ind*2); % assign to 4th sta the same amp from that sta
% %         %%%2nd way is to use the mean amp from the trio station
% %         impindepst(:,9+(ista-3)*2) = mean(impindepst(:,[2 4 6]),2);
%       end
%       
%       %%%given zero-crossing indices, obtain corresponding positive and negative peak indices
%       ppkindep4th = impindepst(:,10:end);
%       npkindep4th = impindepst(:,10:end);
%       for ista = 4:nsta
%         ppkindep4th(:,(ista-4)*2+1) = ppkindep4th(:,(ista-4)*2+1)+ppeaks(ista)-zcrosses(ista);
%         npkindep4th(:,(ista-4)*2+1) = npkindep4th(:,(ista-4)*2+1)+npeaks(ista)-zcrosses(ista);
%       end
%          
%       %%%plot the predicted tarvl and amp of decon impulses at 4th sta vs. waveform
%       [f1] = plt_deconpk_sigpk_comp_4thsta(sigsta(:,4:nsta),stas(4:nsta,:),...
%         impindepst(:,10:end),ppkindep4th,npkindep4th,greenf(:,4:nsta)); 


      %% check the difference by grouping using different stations as the reference station
%       spsscale = sps/40;
%       loff_max = 6*spsscale;  % maximum allowed shift is 4 samples at 40 Hz, it is also true from the synthetic distribution
%       %note the output 'impindep' gives the arrival index of impulse at each station, after
%       %alignment based upon the entire window 'off1i', and the last three cols are the arrival time
%       %difference, NOT the true location yet! 
% %       aa = groupimptripdecon_refOLD(sigdecon,ampit,irccran,rcccat,off1i(k,:),off1iw,loff_max);
%       aa = groupimptriplets_ref(sigdecon,irccran,rcccat,off1i(k,:),off1iw,loff_max,'wtamp',1);
% 
%       for refsta = 1:3
% %         bb = groupimptripdecon_ref(sigdecon,ampit,irccran,rcccat,off1i(k,:),off1iw,loff_max,refsta);
%         bb = groupimptriplets_ref(sigdecon,irccran,rcccat,off1i(k,:),off1iw,loff_max,'wtamp',refsta);
%         size(bb,1)
%         indsame = find(ismember(bb,aa,'rows')); %which sources are the same?
%         inddiff = find(~ismember(bb,aa,'rows'));
%         perca = length(indsame)/size(aa,1)*100
%         percb = length(indsame)/size(bb,1)*100
%       end      


      %% Remove the small-amplitude, secondary triplets from the grouped result
      %convert the sources in terms of the arrival time of zero-crossing to positve peaks' indices
      ppkindep = impindep;  %positive peaks
      for ista = 1: 3
        ppkindep(:,(ista-1)*2+1) = ppkindep(:,(ista-1)*2+1)+ppeaks(ista)-zcrosses(ista);
      end
      npkindep = impindep;  %negative peaks 
      for ista = 1: 3
        npkindep(:,(ista-1)*2+1) = npkindep(:,(ista-1)*2+1)+npeaks(ista)-zcrosses(ista);
      end
      
%       %are there any grouped sources are reversed in arrival order at different sta? 
%       pkist = [];
%       indpk = [];
%       tsep = [];
%       for ista = 1:nsta
%         pki = ppkindep(:,(ista-1)*2+1);
%         [pkist(:,ista),indpk(:,ista)] = sort(pki,'ascend');
%         tsep(:,ista) = diff(pkist(:,ista));
%       end      
%       ind = find(indpk(:,1)~=indpk(:,2) | indpk(:,1)~=indpk(:,3) | indpk(:,2)~=indpk(:,3));
%       disp(ind)
%       length(ind)
      
      %use function 'removesecondarysrc' to remove the smaller amplitude, secondary triplets that 
      %are too close in time to large triplets
      minsrcsep = 20;   % min allowed separation between deconvolved peaks
      [ppkindepsave,indremove] = removesecondarysrc(ppkindep,sigsta(:,1:3));
%       [pkindepsave,indremove] = removesecondarysrc(pkindep,sigsta(:,1:3),minsrcsep);
      
%       if ~noiseflag
%         %%%plot 2ndary sources
%         normflag = 1; %whether to normalize seismograms for plotting
%         f = plt_secondarysrc(ppkindep,sigsta(:,1:3),ppkindepsave,indremove,stas,...
%           xzoom,sps,normflag,detecttype,saveflag);
%       end

      %REMOVE the secondary sources from the grouped result
      impindep(indremove, :) = [];
      ppkindep(indremove, :) = [];
      npkindep(indremove, :) = [];
      sharp(indremove, :) = [];
      impindepst = sortrows(impindep,1);
      nsrc(k,1) = size(impindepst,1);  % number of sources AFTER removing 2ndary

      if ~isempty(impindepst)

%       %plot the sharpness of grouped peaks in res-wlet CC
%       f=plt_srcsharpness(sharp);

      %% plot the scatter of sources in terms of offsets, accounting for prealignment offset
%       span = max(range(off1iw(:,2))+2*loff_max, range(off1iw(:,3))+2*loff_max);
%       xran = [round(mean(minmax(off1iw(:,2)'))-span/2)-1, round(mean(minmax(off1iw(:,2)'))+span/2)+1];
%       yran = [round(mean(minmax(off1iw(:,3)'))-span/2)-1, round(mean(minmax(off1iw(:,3)'))+span/2)+1];
%       cran = [0 lsig];
%       %%%plot the scatter of offsets, accounting for prealignment offset, == true offset
%       f1.fig = figure;
%       f1.fig.Renderer = 'painters';
%       ax1=gca;
%       [ax1,torispl,mamp,xbndcvhl,ybndcvhl] = plt_decon_imp_scatter_ref(ax1,impindepst,xran,yran,...
%         cran,off1iw,loff_max,sps,50,'mean','tori','comb');
%       scatter(ax1,off1i(k,2),off1i(k,3),20,'ks','filled','MarkerEdgeColor','k');
%       title(ax1,'Secondary sources removed');
%       keyboard

%       %%%plot the scatter of offsets, shifted to the same origin, the best alignment
%       xran = [-loff_max+off1i(k,2)-1 loff_max+off1i(k,2)+1];
%       yran = [-loff_max+off1i(k,3)-1 loff_max+off1i(k,3)+1];
%       cran = [0 lsig];
%       f1.fig = figure;
%       f1.fig.Renderer = 'painters';
%       ax1=gca;
%       ax1 = plt_decon_imp_scatter_ref_sft(ax1,impindepst,xran,yran,cran,off1i(k,:),off1iw,...
%         loff_max,irccran,sps,50,'mean','tori');
%       scatter(ax1,off1i(k,2),off1i(k,3),20,'ks','filled','MarkerEdgeColor','k');
%       title(ax1,'Independent, grouped, no secondary sources, shifted');
% keyboard
      %% plot the scatter of sources in terms of rela locations
%       xran = [-5 5];
%       yran = [-5 5];
%       cran = [0 lsig/sps];
%       f2.fig = figure;
%       f2.fig.Renderer = 'painters';
%       ax2=gca;
%       [ax2] = plt_decon_imp_scatter_space_ref(ax2,impindepst,xran,yran,cran,off1iw,loff_max,...
%         sps,50,ftrans,'mean','tori','comb');
%       plot(ax2,xcut,ycut,'k-','linew',2);
%       title(ax2,'Independent, grouped, no secondary sources');
      
      %% num of detections checked in each sub win
      %find which subwin this lfe's time belongs to
      % aa = windows-windows(1,1)+1; %max(floor(tstbuf*sps+1-overshoot-msftaddm),1)
      aa = irccran;
      iwin = findwhichrange(impindepst(:,1),aa);
      ndetwin = zeros(nwin,1);
      for j = 1:nwin
        ndetwin(j) = sum(iwin == j);
      end
      ndetwink{k} = ndetwin;

      %% separation in arrival time between deconvolved positive peaks
      %%Plot the separation in time between these preserved positive peaks after removing the
      %%secondary ones, to see if they can be too close to each other
      if size(ppkindepsave,1) > 1
        [f,tsep] = plt_tsep_deconpk(ppkindepsave,sps);
        close(f.fig);
      end
%       median(tsep)
% keyboard

      %% in sample space, distance for consecutive sourcces, in terms of arrival time?
      ista=1;
      impindepstst = sortrows(impindepst, (ista-1)*2+1);
      tarvlsplst = impindepstst(:,(ista-1)*2+1);
      
      %For each LFE source, get its distance to all other LFEs, maybe we don't care that long separation in time 
      [dt2all,dloc2all,dist2all] = srcdistall(tarvlsplst,impindepstst(:,7:8),[0 2*sps]);
      
%       %plot euclidean distance between each LFE source to all others
%       tlensec = 14981;
%       f = plt_srcdistall(dt2all,dist2all,sps,40/sps,tlensec,1,'spl');
%       %plot the loc diff between each LFE source to all others
%       f = plt_srcdlocall(dloc2all,1,'spl');  
            
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%if instead using other 2 stations as reference, would diff arrival time significantly
      %%%different?
      ista=2;
      aaa = sortrows(impindepst, (ista-1)*2+1);
      bbb = aaa(:,(ista-1)*2+1);
      [dt2all,~,~] = srcdistall(bbb,aaa(:,7:8),[0 2*sps]);
      ista=3;
      aaa = sortrows(impindepst, (ista-1)*2+1);
      bbb = aaa(:,(ista-1)*2+1);
      [dt2all,~,~] = srcdistall(bbb,aaa(:,7:8),[0 2*sps]);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
      m = 5;
      [dtarvl,doffset,eucdist] = srcdistNtoNm(tarvlsplst,impindepstst(:,7:8),m);
      
      % %plot the loc diff between above source pairs
      % f = plt_srcdlocNtoNm(doffset,1,'spl');
      % %plot the diff time and distance between above source pairs
      % f = plt_srcdistNtoNm(dtarvl,eucdist,sps,40/sps,1,'spl'); 
   
      %%%Projected distance along specific directions, eg., propagation and its orthogonal, in terms of arrival time
      [locxyproj,dlocxyproj,stats] = srcprojdistNtoNm(tarvlsplst,impindepstst(:,7:8),m,sps);
      if ~isempty(dlocxyproj)
        nsep = 1;
        ttype = 'tarvl';        
        % [f] = plt_srcprojdist_spl(tarvlsplst,impindepstst(:,7:8),dtarvl{nsep},eucdist{nsep},...
        %   locxyproj,dlocxyproj{nsep},stats,sps,ttype);
%         close(f.fig);
      end            
% keyboard

      %% in map view, the distance for consecutive sourcces, in terms of arrival time?
      %note the 'tsep' obtained from the deconvolved positive peaks should be identical to that if
      %obtained from the deconvolved impulses themselves, which represent the arrival indices of the
      %zero-crossing
      ista=1;
      %convert time offset to relative loc
      [imploc, indinput] = off2space002(impindepst(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
      [impindepstst, indsort] = sortrows(impindepst, (ista-1)*2+1);
      implocst = imploc(indsort, :);
      tarvlsplst = impindepstst(:,(ista-1)*2+1);
      
      %For each LFE source, get its distance to all other LFEs, maybe we don't care that long separation in time 
      [dt2all,dloc2all,dist2all] = srcdistall(tarvlsplst,implocst,[0 2*sps]);

      %in terms of origin time?
      [imploc0, ~] = off2space002([0 0],sps,ftrans,0);  % a ref source at 0,0
      tcor = round((imploc(:,6)-imploc0(6))*sps);   % travel time difference between each source and a ref source at 0,0
      torispl = impindepst(:,1)-tcor;    
      [torisplst, indsort] = sortrows(torispl,1);
      implocorist = imploc(indsort, :);
      impindeporist = impindepst(indsort, :);
      [dto2all,dloco2all,disto2all] = srcdistall(torisplst,implocorist,[0 2*sps]);
      
      % %plot euclidean distance between each LFE source to all others
      % f = plt_srcdistall(dt2all,dist2all,sps,40/sps,0.1,'km');  
      % %plot the loc diff between each LFE source to all others
      % f = plt_srcdlocall(dloc2all,0.1,'km');

      %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
      m = 8;
      [dtarvl,dneloc,eucdist] = srcdistNtoNm(tarvlsplst, implocst, m);

      % %plot the loc diff between above source pairs
      % f = plt_srcdlocNtoNm(dneloc,0.1,'km');
      % %plot the diff time and distance between above source pairs
      % f = plt_srcdistNtoNm(dtarvl,eucdist,sps,40/sps,0.1,'km'); 
        
      %%%Projected distance along specific directions, eg., propagation and its orthogonal, in terms of arrival time
      % [locxyproj,dlocxyproj,stats] = srcprojdistNtoNm(tarvlsplst,implocst,m,sps);
      % if ~isempty(dlocxyproj)
      %   nsep = 1;
      %   ttype = 'tarvl';
      %   wt = median(impindepst(:,[2 4 6]),2);
      %   [f] = plt_srcprojdist(tarvlsplst,implocst,dtarvl{nsep},eucdist{nsep},...
      %     locxyproj,dlocxyproj{nsep},stats,sps,ttype,wt);
%         close(f.fig);
      % end  

      if isempty(impindepst)
        continue
      else
      srcampr = [impindepst(:,2)./impindepst(:,4) impindepst(:,2)./impindepst(:,6) ...
                 impindepst(:,4)./impindepst(:,6)];
               
      psrcamprs = [impindepst(:,2)*max(greenf(:,1))./(impindepst(:,4)*max(greenf(:,2))) ...
                   impindepst(:,2)*max(greenf(:,1))./(impindepst(:,6)*max(greenf(:,3))) ...
                   impindepst(:,4)*max(greenf(:,2))./(impindepst(:,6)*max(greenf(:,3)))];
      psrcamps = [impindepst(:,2)*max(greenf(:,1)) impindepst(:,4)*max(greenf(:,2)) ...
                  impindepst(:,6)*max(greenf(:,3))];
      
      nsrcamprs = [impindepst(:,2)*min(greenf(:,1))./(impindepst(:,4)*min(greenf(:,2))) ...
                   impindepst(:,2)*min(greenf(:,1))./(impindepst(:,6)*min(greenf(:,3))) ...
                   impindepst(:,4)*min(greenf(:,2))./(impindepst(:,6)*min(greenf(:,3)))];
      nsrcamps = [impindepst(:,2)*min(greenf(:,1)) impindepst(:,4)*min(greenf(:,2))...
                  impindepst(:,6)*min(greenf(:,3))];
      
      msrcampr(k,:) = median(log10(srcampr), 1);
      madsrcampr(k,:) = mad(log10(srcampr), 1, 1);
      mpsrcamprs(k,:) = median(log10(psrcamprs), 1);
      madpsrcamprs(k,:) = mad(log10(psrcamprs), 1, 1);
      mnsrcamprs(k,:) = median(log10(nsrcamprs), 1);
      madnsrcamprs(k,:) = mad(log10(nsrcamprs), 1, 1);

      %what is the deviation of amp ratio from the median for each source?
      lndevsrcampr = srcampr-median(srcampr, 1); % in linear scale
      lgdevsrcampr = log10(srcampr)-median(log10(srcampr), 1); % in log scale, note that log2+log5=log10, so this means a ratio

      %%%what are the corresponding RCC at each source
      rccpairsrc = [];
      rccpairsrc(:,1) = rccpaircat(round(mean(impindepst(:,[1 3]),2)),1);
      rccpairsrc(:,2) = rccpaircat(round(mean(impindepst(:,[1 5]),2)),2);
      rccpairsrc(:,3) = rccpaircat(round(mean(impindepst(:,[3 5]),2)),3);
      
      %use the concatenated rcc at the average arrival time of each source
      rcccatsrc = [];
      rcccatsrc(:,1) = rcccat(round(mean(impindepst(:,[1 3 5]),2)));

      %% signal + zoom-in + map locations + some reference symbols 
      if ~noiseflag
        
        %%%plot all short windows, and LFE detections
        htin = 6.3;   % maximum height allowed is 11 inches
        widin = 8.4;  % maximum width allowed is 8.5 inches
        nrow = 1;
        ncol = 1;
        f = initfig(widin,htin,nrow,ncol); %initialize fig
        figxran = [0.08 0.98]; figyran = [0.08 0.98];
        figxsep = 0.05; figysep = 0.02;
        optaxpos(f,nrow,ncol,figxran,figyran,figxsep,figysep);
        ax=f.ax(1); hold(ax,'on'); ax.Box = 'off'; grid(ax,'off');
        %       set(ax, 'color', 'none');
        %   xlim(ax,[]);
        %   nolabels(ax,2);
        longticks(ax,4);
        axis(ax,'tight');
        ylabel(f.ax(end),sprintf('Amplitude'),'FontSize',10);
        xlabel(f.ax(end),sprintf('Time (s) on %s %s %s',dy,mo,yr),'FontSize',10);
        vsep = 1;
%         ylim(ax,[-vsep*(nwin+1) 0]);
        
        [clusibst,nclusibst,tclusibst]=clusters_in_burst(181);
        
        cumutime=0;
        for iwin = 1:nwin
          subwtmp = subwsave{iwin};
          plot(ax,(1:length(subwtmp))/sps,subwtmp(:,1)-vsep*iwin,'r','linew',1);
          plot(ax,(1:length(subwtmp))/sps,subwtmp(:,2)-vsep*iwin,'b','linew',1);
          plot(ax,(1:length(subwtmp))/sps,subwtmp(:,3)-vsep*iwin,'k','linew',1);
%           %%%plot the deconvolved sources
%           yloc = -vsep*(iwin+0.45);
%           tst = cumutime;
%           ted = cumutime+length(subwtmp)/sps;
% %           ind = find(impindepst(:,1)/sps >= tst & impindepst(:,1)/sps < ted);
% %           scatter(ax,impindepst(ind,1)/sps-tst,yloc*ones(size(impindepst(ind,1))),...
% %             6,[.6 .6 .6],'filled');
%           ind = find(tclusibst(:,1)>=tst & tclusibst(:,2)<ted);
%           clusibstxzoom = clusibst(ind);
%           tclusibstxzoom = tclusibst(ind);
%           nclus = length(ind);
%           for iclus=1:2:nclus
%             impclus=clusibstxzoom{iclus};
%             scatter(ax,impclus(:,1)/sps-tst,yloc*ones(size(impclus(:,1))),10,'^','k');
%           end
%           for iclus=2:2:nclus
%             impclus=clusibstxzoom{iclus};
%             scatter(ax,impclus(:,1)/sps-tst,yloc*ones(size(impclus(:,1))),10,'v','k');
%           end
%           impuniclus = unique(cat(1,clusibstxzoom{:}),'rows');
%           ind = find(impindepst(:,1)/sps >= tst & impindepst(:,1)/sps <ted);
%           impiso = setdiff(impindepst(ind,:),impuniclus,'rows');
%           scatter(ax,impiso(:,1)/sps-tst,yloc*ones(size(impiso(:,1))),10,'k');
          
          cumutime = ted;
        end
        
%         xzoom = [0 25];
%         [f] = plt_agu2022abstractv4(greenf(:,1:3),optdat(:,2:4),impindepst,sps,xzoom,off1iw,loff_max,...
%           rcccat,overshoot,tstbuf,dy,mo,yr,ftrans,'spl');
%         text(f.ax(3),0.98,0.88,'Using synthetic noise','HorizontalAlignment','right',...
%           'Units','normalized','FontSize',9,'FontWeight','bold');
%       else        
        xzoom = [5 30];
        tbostplt=tbosti(indbi); %bostock's LFE detections inside the burst, pointing to positive peaks
        tbostzcplt=tbostplt-(ppeaks(1)-zcrosses(1))/sps; %shift it to zero-crossings
        winnoi = [1 3]; %empirical window for pre-arrival noise
        winsig = [4 7]; %empirical window for main arrival    
        patchwins = [winnoi;winsig];
        f=plt_002lfe_signalzoom(green,greenf,sigsta,[],impindepst,sps,xzoom,...
          tstbuf,dy,mo,yr,tbostzcplt,patchwins,saveflag);
% %         [f] = plt_agu2022abstractv3(greenf(:,1:3),sigsta(:,1:3),impindepst,sps,xzoom,off1iw,loff_max,...
% %           tstbuf,dy,mo,yr,ftrans,'spl');
      end
      % keyboard
%       text(f.ax(3),0.98,0.95,'Secondary removed','HorizontalAlignment','right',...
%         'Units','normalized','FontSize',10);
% %       text(f.ax(3),0.98,0.95,'Further checked at KLNB','HorizontalAlignment','right',...
% %         'Units','normalized','FontSize',10);
%       text(f.ax(4),0.98,0.95,'Zoom-in','HorizontalAlignment','right',...
%         'Units','normalized','FontSize',10);
% 
%       if ~noiseflag
%         ax=f.ax(3);
%         hold(ax,'on');
%         xran1 = [-4 4];
%         yran1 = [-4 4];
%         wtmax = prctile(mean(impindepst(:,[2 4 6]),2),95); %use percentile in case
%         scatter(ax,xran1(1)+0.1*range(xran1),yran1(2)-0.05*range(yran1),35,'w','filled',...
%           'MarkerEdgeColor',[.5 .5 .5]);
%         text(ax,0.02,0.9,'amplitude ','Units','normalized',...
%           'HorizontalAlignment','left','FontSize',8);
%         text(ax,0.02,0.85,'\geq 95th prctile','Units','normalized',...
%           'HorizontalAlignment','left','FontSize',8);
% %         errloc = errloc(:,1:2)+[-0.5 -0.5];
% %         plot(ax,errloc(:,1),errloc(:,2),'-','Color',[.4 .4 .4],'linew',1);
%         angrmse = stats.angrmse;
%         [rotx, roty] = complex_rot(0,1,-angrmse);
%         xvect = [0.5-rotx 0.5+rotx];
%         yvect = [-2.5-roty -2.5+roty];
%         drawArrow(ax,xvect,yvect,xran1,yran1,'linewidth',1);
%         text(ax,0.62,0.18,strcat(num2str(angrmse),'$^{\circ}$'),'FontSize',10,...
%           'unit','normalized','interpreter','latex');
%         hold(ax,'off');
%       
%         %%%add the indication of the reasonable physical source size
%         ax=f.ax(4);
%         hold(ax,'on');
%         x0 = 2;
%         y0 = -1;
%         Vs = 3;   % S-wave speed
% %       Vprop = 0.8*Vs;   % reasonable rupture propagation speed
%         Vprop = Vs;   % max rupture propagation speed
%         td = 40/sps;  % estimated from the bin where the most tsep falls in, 32-48-sample bin
% %         td = median(tsep(:));
%         srcsz = Vprop*td;
%         radi = srcsz/2;
%         scatter(ax,x0,y0,2,'ko','filled');
%         [x, y] = circle_chao(x0,y0,radi,0.1);
%         plot(ax,x,y,'-','Color','k','linew',1);
%         text(ax,0.55,0.25,'max. reasonable','HorizontalAlignment','left',...
%           'Units','normalized','FontSize',10,'interpreter','latex');
%         text(ax,0.55,0.2,'source size','HorizontalAlignment','left',...
%           'Units','normalized','FontSize',10,'interpreter','latex');
%         text(ax,0.55,0.15,strcat('$V_{s} \cdot t_{d} = 3 \cdot 0.25$'),'HorizontalAlignment','left',...
%           'Units','normalized','FontSize',10,'interpreter','latex');
% 
%         %%%add the indication of the error ellipse
%         mmmax = 6;
%         nnmax = 6;
%         erroff1 = zeros((mmmax*2+1)^2,2);
%         kk = 0;
%         for mm = -mmmax:1:mmmax
%           for nn = -nnmax:1:nnmax
%             kk = kk+1;
%             erroff1(kk,:) = [mm nn];
%           end
%         end
%         [errloc1, ~] = off2space002(erroff1,sps,ftrans,0);
%         
%         F1 = scatteredInterpolant(erroff1(:,1),erroff1(:,2),errloc1(:,1),'linear');
%         F2 = scatteredInterpolant(erroff1(:,1),erroff1(:,2),errloc1(:,2),'linear');
%         
% %         erroff = 3*[0 2; 1 sqrt(3); sqrt(2) sqrt(2); sqrt(3) 1;
% %           2 0; sqrt(3) -1; sqrt(2) -sqrt(2); 1 -sqrt(3);
% %           0 -2; -1 -sqrt(3); -sqrt(2) -sqrt(2); -sqrt(3) -1;
% %           -2 0; -sqrt(3) 1; -sqrt(2) sqrt(2); -1 sqrt(3);
% %           0 2;];
% %         errloc = [];
% %         errloc(:,1) = F1(erroff(:,1),erroff(:,2));
% %         errloc(:,2) = F2(erroff(:,1),erroff(:,2));
%         
%         [off1iloc, ~] = off2space002(off1i(k,2:3),sps,ftrans,0);
%         errloc = errloc(:,1:2)+off1iloc(1:2);
%         
%         plot(ax,errloc(:,1),errloc(:,2),'-','Color',[.4 .4 .4],'linew',1);
%         plot(f.ax(3),errloc(:,1),errloc(:,2),'-','Color',[.4 .4 .4],'linew',1);
%         text(ax,0.05,0.9,strcat('$\pm6$','-sample'),'HorizontalAlignment','left',...
%           'Units','normalized','FontSize',10,'interpreter','latex');
%         text(ax,0.05,0.85,'error contour','HorizontalAlignment','left',...
%           'Units','normalized','FontSize',10,'interpreter','latex');
%       end
%       keyboard
      %% test a bunch of offset max to see residual reduction VS. # of sources left
%       %%%According to the resulting plot, ~2* RMSE (15 samples at 160 sps) seems proper
%       modname = 'timeoff_plfit_4thsta_160sps.mat';
%       planefit = load(strcat(rstpath, '/MAPS/',modname));
%       ista = 7;
%       [nsrctest(k,:),nsrcnorm(k,:),rr(k,:,:),vr(k,:,:)] = ...
%         testvarredvsoffmax(ista,sigsta,optdat,greenf,sps,zcrosses,planefit,irccran,...
%         rcc1icat,impindep,stas,off1i(k,:),off1iw,overshoot,mswccpksep,1);      
%       keyboard
      %% 2ndary src removed, prediction of impulse tarvl at 4th sta given sources and empirical off14-src relation
%       %%%predict the tarvl and amp at 4th sta 
%       for ista = 4:nsta
%         impindepst(:,9+(ista-4)*2+1) = pred_tarvl_at4thsta(stas(ista,:),impindepst(:,7),...
%           impindepst(:,8),impindepst(:,1),off1i(k,ista));
%         %%%there is decision making here about what impulse amp should be given to 4th sta
%         [~,ind] = max(mcoef(ista-3,:)); % which of the trio stas is most similar to 4th sta?
% 
%         %1st way is to use the amp SAME as of the sta that is most similar in waveform to 4th
% %         impindepst(:,9+(ista-3)*2) = impindepst(:,ind*2); %assign to 4th sta the same amp from that sta
% %         impindepst(:,9+(ista-3)*2) = impindepst(:,9+(ista-3)*2)*envrat(ista-3,ind); %account for envelope ratio
%         
% %         %2nd way is to use the mean amp from the trio station
%         impindepst(:,9+(ista-3)*2) = mean(impindepst(:,[2 4 6]),2);
% %         impindepst(:,9+(ista-3)*2) = impindepst(:,9+(ista-3)*2)*envrat(ista-3,ind); %account for envelope ratio
%       end

      %%%carry out 'deconvolution' at 4th stations as well for the tarvl and amp
      modname = 'timeoff_plfit_4thsta_160sps.mat';
      planefit = load(strcat(rstpath, '/MAPS/',modname));
%       rmse = planefit.gof.rmse;

      pred4diff = [];
      off14pred = [];
      for ista = 4:nsta
        wlet = greenf(:,ista);  %template here is best aligned, tapered, linear trend removed, filtered
        lwlet = length(wlet);
        sig = sigsta(:,ista); %best aligned, filtered, tapered
        lsig = length(sig);
        
        dt = 1/sps;  % sampling interval
        twlet = zcrosses(ista)*dt;
        fpltit = 0;  % plot flag for each iteration
        fpltend = 0;  % plot flag for the final iteration
        fpltchk = 0; % plot flag for intermediate computations
        rmse = planefit.gof{ista-3}.rmse;
        %1.6*rmse seems proper from whole-win RCC, 2 from concatenated RCC
        offmax = round(2.0*rmse);
        
        if noiseflag
          %%%As of 2022/09/27, 'threshold' would be auto-computed based on data then feed to the
          %%%decon if the 'noiseflag' is on
          [sigdecon(:,ista),pred,res,dresit,mfitit{ista},ampit{ista},fighdl] = ...
            iterdecon_4thsta(sig,wlet,irccran,rcc1icat(:,ista-3),fixthresh(ista),...
            dt,twlet,impindep,stas(ista,:),off1i(k,ista),off1iw(:,ista),offmax,...
            fpltit,fpltend,fpltchk);
        else
%         [sigdecon(:,ista),pred,res,dresit,mfitit,ampit{ista},fighdl] = iterdecon_4thsta(sig,wlet,...
%           irccran,rcc1i(:,ista-3),dt,twlet,impindep,stas(ista,:),off1i(k,ista),[],offmax,...
%           fpltit,fpltend,fpltchk);
          [sigdecon(:,ista),pred,res,dresit,mfitit{ista},ampit{ista},fighdl] = ...
            iterdecon_4thsta(sig,wlet,irccran,rcc1icat(:,ista-3),[],...
            dt,twlet,impindep,stas(ista,:),off1i(k,ista),off1iw(:,ista),offmax,...
            fpltit,fpltend,fpltchk);
        end
        
        ampiti = ampit{ista};
        impindep(:,9+(ista-4)*2+1) = ampiti(:,1);
        impindep(:,9+(ista-3)*2) = ampiti(:,2);
        if ista == nsta
          pred4diff(:,ista-3) = ampiti(:,end-1);  %difference between found peak and predicted arrival
          off14pred(:,ista-3) = ampiti(:,end);  %predicted off14
        end
      end
      
      %%%given zero-crossing indices, obtain corresponding positive and negative peak indices
      ppkindep = impindep;
      npkindep = impindep;
      for ista = 1: 3
        ppkindep(:,(ista-1)*2+1) = ppkindep(:,(ista-1)*2+1)+ppeaks(ista)-zcrosses(ista);
        npkindep(:,(ista-1)*2+1) = npkindep(:,(ista-1)*2+1)+npeaks(ista)-zcrosses(ista);
      end
      for ista = 4:nsta
        ppkindep(:,9+(ista-4)*2+1) = ppkindep(:,9+(ista-4)*2+1)+ppeaks(ista)-zcrosses(ista);
        npkindep(:,9+(ista-4)*2+1) = npkindep(:,9+(ista-4)*2+1)+npeaks(ista)-zcrosses(ista);
      end
         
%       %%%plot the predicted tarvl and amp of decon impulses at 4th sta vs. waveform
%       [f1] = plt_deconpk_sigpk_comp_4thsta(sigsta(:,4:nsta),stas(4:nsta,:),...
%         impindep(:,10:end),ppkindep4th,npkindep4th,greenf(:,4:nsta));
      
      %%%further ELIMINATE sources that fail the check at 4th stations
      trust4th = 7; % trust KLNB the most among all 4th stations
      indremove = find(impindep(:,9+(trust4th-4)*2+1)==0 & impindep(:,9+(trust4th-3)*2)==0);
      pred4difftr = pred4diff(setdiff(1:size(pred4diff,1),indremove),trust4th-3);
      off14predtr = off14pred(setdiff(1:size(off14pred,1),indremove),trust4th-3);
      impindep(indremove,:) = [];
      ppkindep(indremove, :) = [];
      npkindep(indremove, :) = [];
      [impindepst,indsort] = sortrows(impindep,1);
      nsrc4th(k,1) = size(impindepst,1);  % number of sources AFTER 4th-sta check
      pred4difftr = pred4difftr(indsort);
      
      if ~isempty(impindepst)
%       %%%plot the scatter of sources in terms of offsets, accounting for prealignment offset
%       span = max(range(off1iw(:,2))+2*loff_max, range(off1iw(:,3))+2*loff_max);
%       xran = [round(mean(minmax(off1iw(:,2)'))-span/2)-1, round(mean(minmax(off1iw(:,2)'))+span/2)+1];
%       yran = [round(mean(minmax(off1iw(:,3)'))-span/2)-1, round(mean(minmax(off1iw(:,3)'))+span/2)+1];
%       cran = [0 lsig];
%       %%%plot the scatter of offsets, accounting for prealignment offset, == true offset
%       f1.fig = figure;
%       f1.fig.Renderer = 'painters';
%       ax1=gca;
%       [ax1,torispl,mamp,xbndcvhl,ybndcvhl] = plt_decon_imp_scatter_ref(ax1,impindepst,xran,yran,cran,off1iw,loff_max,...
%         sps,50,'mean','tori','comb');
%       scatter(ax1,off1i(k,2),off1i(k,3),20,'ks','filled','MarkerEdgeColor','k');
%       title(ax1,'Secondary sources & 4th-sta checked');

      %% num of detections checked in each sub win
      %find which subwin this lfe's time belongs to
      iwin = findwhichrange(impindepst(:,1),aa);
      ndetwin4th = zeros(nwin,1);
      for j = 1:nwin
        ndetwin4th(j) = sum(iwin == j);
      end
      ndetwin4thk{k} = ndetwin4th;

      %% error in off14 between prediction and observation
      %off14 prediction for decon srcs using plane fit model WITH alignment 
%       [~,off14pred] = pred_tarvl_at4thsta(stas(trust4th,:),impindepst(:,7),impindepst(:,8),...
%         impindepst(:,1),0);
      %off14 computed from actually-matched arrivals at stas 1 and 4 from decon
      % +repmat(off1i(trust4th),size(impindepst,1),1)
%       off14 = impindepst(:,1)-impindepst(:,9+(trust4th-4)*2+1); %after prealignment
%       diffoff14tr = off14pred-off14;
      
      off14 = impindep(:,1)-impindep(:,9+(trust4th-4)*2+1); %after prealignment
      diffoff14tr = off14predtr-off14;
      diffoff14tr = diffoff14tr(indsort);

      f = initfig(8,5,1,2); %initialize fig
      ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); 
      histogram(ax,pred4difftr);
      plot(ax,[median(pred4difftr) median(pred4difftr)],ax.YLim,'r--','LineWidth',1);
      plot(ax,[-offmax -offmax],ax.YLim,'k--');
      plot(ax,[offmax offmax],ax.YLim,'k--');
      xlabel(ax,'diff. in 4th arrival between pred and decon');
      ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');  
      histogram(ax,diffoff14tr);
      plot(ax,[median(diffoff14tr) median(diffoff14tr)],ax.YLim,'r--','LineWidth',1);
      plot(ax,[-offmax -offmax],ax.YLim,'k--');
      plot(ax,[offmax offmax],ax.YLim,'k--');
      xlabel(ax,'diff. in off14 between plane-fit and decon');
      ylabel(ax,'count');

      if ~noiseflag
        %%%final prediction via convolution between grouped impulses and template at each station
        [f2,predgrp,resgrp,predgrpl,resgrpl,l2normred1,varred4th]=...
          predsig_conv_imptemp(sigsta,optdat,impindepst,greenf,zcrosses,...
          overshoot,stas,0,sps);
        %%%simply plot the trio stations and KLNB signal, prediction and residual
        pltsta = [1 2 3 7];
%         f1 = plt_selsigpredres(sigsta,predgrp,resgrp,varred1,stas,pltsta,sps,...
%           xzoom,detecttype,saveflag);
        f1 = plt_selsigpred(sigsta,predgrp,varred4th,stas,pltsta,sps,...
          xzoom,detecttype,saveflag);
%%
        %%%final prediction at ort. comp.
        [~,predgrport,resgrport,~,~,l2normred1ort,varredort4th]= ...
          predsig_conv_imptemp(sigstaort,ortdat,impindepst,...
          greenfort,zcrosses,overshoot,stas,0,sps);
        %%%prediction at vert. comp.
        [~,predgrpvert,resgrpvert,~,~,l2normred1vert,varredvert4th]= ...
          predsig_conv_imptemp(sigstavert,vertdat,impindepst,...
          greenfvert,zcrosses,overshoot,stas,0,sps);
        normflag=1;
        f = plt_selsigpred_allcomp(sigsta,predgrp,varred4th,sigstaort,predgrport,...
          varredort4th,sigstavert,predgrpvert,varredvert4th,greenfort,greenfvert,...
          stas,pltsta,sps,xzoom,patchwins,normflag,detecttype,saveflag);

        f = plt_selsigpred_ortvert(sigstaort,predgrport,varredort4th,...
          sigstavert,predgrpvert,varredvert4th,greenfort,greenfvert,...
          stas,pltsta,sps,xzoom,patchwins,detecttype,saveflag);

      end
      
% %%
% impindepst = imp4thrmst;

      %% recompute time separation and distance after 4th sta check  
      %%%Plot the separation in time between these preserved positive peaks after removing the
      %%secondary ones, to see if they can be too close to each other
      %discard the sources that are determined to be too close and secondary compared to a major source
      ppkindepsave = ppkindep;
      if size(ppkindepsave,1) > 1
        [f,tsep] = plt_tsep_deconpk(ppkindepsave,sps);
        close(f.fig);
      end
%       median(tsep)

      %%%in sample space, distance for consecutive sourcces in terms of arrival time
      ista=1;
      impindepstst = sortrows(impindepst, (ista-1)*2+1);
      tarvlsplst = impindepstst(:,(ista-1)*2+1);
      
      %%
%       %%%obtain a single best alignment based on the zoom-in segment
%       xzoom = [5 30];
%       msftadd = 20;
%       optcc = sigsta(xzoom(1)*sps+1: xzoom(2)*sps, :);
%       ccmid = ceil(size(optcc,1)/2);
%       ccwlen = round(size(optcc,1)-2*(msftadd+1));
%       loffmax = 4*sps/40;
%       ccmin = 0.01;  % depending on the length of trace, cc could be very low
%       iup = 1;    % times of upsampling
%       [off12con,off13con] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
%         ccwlen,msftadd,loffmax,ccmin,iup);
%       % if a better alignment cannot be achieved, use 0,0
%       if off12con == msftadd+1 && off13con == msftadd+1
%         off12con = 0;
%         off13con = 0;
%         fprintf('This segment cannot be properly aligned, double-check needed \n');
%       end
%       sigseg(:, 1) = sigsta(xzoom(1)*sps+1: xzoom(2)*sps, 1);
%       sigseg(:, 2) = sigsta(xzoom(1)*sps+1-off12con: xzoom(2)*sps-off12con, 2); % sta 2
%       sigseg(:, 3) = sigsta(xzoom(1)*sps+1-off13con: xzoom(2)*sps-off13con, 3); % sta 3
%       
%       [mcoef,off17] = xcorrmax(optcc(:,1), optcc(:,nsta), msftadd, 'coeff');
%       sigseg(:, 4) = sigsta(xzoom(1)*sps+1-off17: xzoom(2)*sps-off17, nsta); % sta 3
% 
%       [ircc,rcc12] = RunningCC(sigseg(:,1), sigseg(:,2), rccmwlen);
%       [~,rcc13] = RunningCC(sigseg(:,1), sigseg(:,3), rccmwlen);
%       [~,rcc23] = RunningCC(sigseg(:,2), sigseg(:,3), rccmwlen);
%       % ircc = ircc+*sps;
%       rcc = (rcc12+rcc13+rcc23)/3;
% %       rcc = (rcc12+rcc13)/2;
      
%       %plot a zoomed-in view of seismograms
%       timp = impindepall(:,1)+ppeaks(1)-zcrosses(1);
%       timp4th = impindep4thall(:,1)+ppeaks(1)-zcrosses(1);
%       f=plt_seiszoomedin(sigseg,xzoom,sps,tstbuf,tbosto,...
%         tbosti,timp,timp4th,ircc,rcc);        
%       keyboard

      %For each LFE source, get its distance to all other LFEs, maybe we don't care that long separation in time 
      [dt2all,dloc2all,dist2all] = srcdistall(tarvlsplst,impindepstst(:,7:8),[0 2*sps]);
      
      % %plot euclidean distance between each LFE source to all others
      % f = plt_srcdistall(dt2all,dist2all,sps,40/sps,1,'spl');  
      % %plot the loc diff between each LFE source to all others
      % f = plt_srcdlocall(dloc2all,1,'spl');

      %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
      m = 5;
      [dtarvl,doffset,eucdist] = srcdistNtoNm(tarvlsplst, impindepstst(:,7:8), m);

      % %plot the loc diff between above source pairs
      % f = plt_srcdlocNtoNm(doffset,1,'spl');
      % %plot the diff time and distance between above source pairs
      % f = plt_srcdistNtoNm(dtarvl,eucdist,sps,40/sps,1,'spl'); 
    
      %%%Projected distance along specific directions, eg., propagation and its orthogonal, in terms of arrival time
      [locxyproj,dlocxyproj,stats] = srcprojdistNtoNm(tarvlsplst,impindepstst(:,7:8),m,sps);
      if ~isempty(dlocxyproj)
        nsep = 1;
        ttype = 'tarvl';        
        % [f] = plt_srcprojdist_spl(tarvlsplst,impindepstst(:,7:8),dtarvl{nsep},eucdist{nsep},...
        %   locxyproj,dlocxyproj{nsep},stats,sps,ttype);
%         close(f.fig);
      end
% keyboard

      %%%map view, distance in terms of arrival time 
      ista=1;      
      [imploc, indinput] = off2space002(impindepst(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
      [impindepstst, indsort] = sortrows(impindepst, (ista-1)*2+1);
      implocst = imploc(indsort, :);
      tarvlsplst = impindepstst(:,(ista-1)*2+1);
            
      %For each LFE source, get its distance to all other LFEs, maybe we don't care that long separation in time 
      [dt2all,dloc2all,dist2all] = srcdistall(tarvlsplst,implocst,[0 2*sps]);

      %in terms of origin time?
      [imploc0, ~] = off2space002([0 0],sps,ftrans,0);  % a ref source at 0,0
      tcor = round((imploc(:,6)-imploc0(6))*sps);   % travel time difference between each source and a ref source at 0,0
      torispl = impindepst(:,1)-tcor;    
      [torisplst4th, indsort] = sortrows(torispl,1);
      implocorist4th = imploc(indsort, :);
      impindeporist4th = impindepst(indsort, :);
      [dto2all,dloco2all,disto2all] = srcdistall(torisplst4th,implocorist4th,[0 2*sps]);
      
      % %plot euclidean distance between each LFE source to all others
      % f = plt_srcdistall(dt2all,dist2all,sps,40/sps,0.1,'km');  
      % %plot the loc diff between each LFE source to all others
      % f = plt_srcdlocall(dloc2all,0.1,'km');
      
      %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
      m = 8;
      [dtarvl,dneloc,eucdist] = srcdistNtoNm(tarvlsplst, implocst, m);

      % %plot the loc diff between above source pairs
      % f = plt_srcdlocNtoNm(dneloc,0.1,'km');
      % %plot the diff time and distance between above source pairs
      % f = plt_srcdistNtoNm(dtarvl,eucdist,sps,40/sps,0.1,'km'); 
      
      %%%Projected distance along specific directions, eg., propagation and its orthogonal, in terms of arrival time
      % [locxyproj,dlocxyproj,stats] = srcprojdistNtoNm(tarvlsplst,implocst,m,sps);
      % if ~isempty(dlocxyproj)
      %   nsep = 1;
      %   ttype = 'tarvl';
      %   wt = median(impindepst(:,[2 4 6]),2);
      %   [f] = plt_srcprojdist(tarvlsplst,implocst,dtarvl{nsep},eucdist{nsep},...
      %     locxyproj,dlocxyproj{nsep},stats,sps,ttype,wt);
%         close(f.fig);
      % end

      %%      
      %%%Both 3-sta and 4-sta detections, map locations and projections
      ttype = 'tori';
      [locxyproj,~,stats] = srcprojdistNtoNm(torisplst,implocorist,1,sps);      
      [locxyproj4th,~,stats4th] = srcprojdistNtoNm(torisplst4th,implocorist4th,1,sps);
      
      %%%independent search for prop direc, etc. to the removed sources
      [~,indremove] = setdiff(implocorist,implocorist4th,'rows','stable');
      implocoristrem = implocorist(indremove,:);
      torisplstrem = torisplst(indremove,:);
      [locxyprojrem,~,statsrem] = srcprojdistNtoNm(torisplstrem,implocoristrem,1,sps);

      %%
      if noiseflag
        % %if the noise seismogram is NOT needed, use 'plt_shortwin_demo_noise'
        % f=plt_shortwin_demo_noise(impindeporist,torisplst,implocorist,...
        %   impindeporist4th,torisplst4th,implocorist4th,lsig,sps,ttype,saveflag);
        %if the noise seismogram is needed as well, use 'plt_wholewin_summary_noise'
        f=plt_wholewin_summary_noise(sigsta(:,1:3),impindeporist,torisplst,...
          implocorist,impindeporist4th,torisplst4th,implocorist4th,...
          tstbuf,dy,mo,yr,sps,ttype,'demo',0);  
        fname = strcat('shortwinsum_noi_demo.pdf');
        print(f.fig,'-dpdf',...
          strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));          
      else
        %%%summary figure for short-win data of the demo burst, w/o waveforms
        % f=plt_shortwin_demo(impindeporist,torisplst,implocorist,...
        %   locxyproj,stats,impindeporist4th,torisplst4th,implocorist4th,...
        %   locxyproj4th,stats4th,lsig,sps,ttype,saveflag);
        %%%particularly for burst #181, plus tremor locations
        f=plt_shortwin_bst181_tmr(impindeporist,torisplst,implocorist,...
          locxyproj,stats,impindeporist4th,torisplst4th,implocorist4th,...
          locxyproj4th,stats4th,lsig,sps,ttype,saveflag,hfall);
          
      end
  

% keyboard

      end
      end
      %% diff between predicted arrival and selected peak at 4th sta
%       figure; 
%       ax = gca;
%       hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%       yyaxis(ax,'left');
%       histogram(ax,abs(pred14offtr),'Normalization','probability','BinWidth',1*sps/40);
%       p1=plot(ax,[offmax offmax],ax.YLim,'k--');
% %       [muHat,sigmaHat] = normfit(pred14offtr);
% %       pdfhat = normpdf(-offmax:offmax,muHat,sigmaHat);
% %       plot(ax,-offmax:offmax,pdfhat,'r','linew',2);
% %       plot(ax,[muHat muHat],ax.YLim,'r--','linew',1);
%       xlabel(ax,'Abs diff in samples between predicted arrival and selected peak');
%       ylabel(ax,'Probability');
%       xlim(ax,[0 offmax]);
%       
%       yyaxis(ax,'right');
%       [cdfval,x] = ecdf(abs(pred14offtr)); %between Nth and (N-1)th source
%       plot(ax,x,cdfval,'linew',1,'color','r');
%       ylabel(ax,'Empirical CDF');
%       legend(ax,[p1],'max allowed diff','Location','east');
      
      %% preserved sources' amp ratio between 4th and 1st stas
%       figure;
%       set(gcf,'Renderer','painters');
%       ax = subplot(1,3,1);
%       hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%       scatter(ax,impindep(:,2),impindep(:,9+(trust4th-3)*2),20,'k');
%       axis(ax,'equal');
% %       axis(ax,[0 1.2 0 1.2]);
%       plot(ax,[0 max([ax.YLim(2) ax.XLim(2)])],[0 max([ax.YLim(2) ax.XLim(2)])],'r--','linew',1);
%       xlabel(ax,'Amp at sta 1');
%       ylabel(ax,'Amp at sta 4');
% 
%       ax = subplot(1,3,2);
%       hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%       histogram(ax,log10(impindep(:,9+(trust4th-3)*2)./impindep(:,2)),'Normalization','pdf',...
%         'BinWidth',0.1);
%       quant = log10(impindep(:,9+(trust4th-3)*2)./impindep(:,2));
%       [muHat,sigmaHat] = normfit(quant);
%       pdfhat = normpdf(-1:0.02:1,muHat,sigmaHat);
%       plot(ax,-1:0.02:1,pdfhat,'r-','linew',2);
%       plot(ax,[muHat muHat],ax.YLim,'r--','linew',1);      
%       plot(ax,[median(quant) median(quant)],ax.YLim,'b--','linew',1);
%       text(ax,0.05,0.9,sprintf('1*sigma=%.2f',sigmaHat),'Units','normalized');
%       xlabel(ax,'log_{10}{Amp ratio sta 4/1}');
%       ylabel(ax,'PDF');
% 
%       ax = subplot(1,3,3);
%       hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%       scatter(ax,log10(impindep(:,2)),log10(impindep(:,9+(trust4th-3)*2)./impindep(:,2)),20,'k');
%       xlabel(ax,'log_{10}{Amp ratio at 1}');
%       ylabel(ax,'log_{10}{Amp ratio 4/1}');
%       axis(ax,'equal');
      
%       ax = subplot(1,3,3);
%       hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%       scatter(ax,impindep(:,2),impindep(:,9+(trust4th-3)*2)./impindep(:,2),20,'k');
%       xlabel(ax,'Amp ratio at 1');
%       ylabel(ax,'Amp ratio 4/1');
%       axis(ax,'equal');

%       keyboard


      %% plot amp ratio 12 and 13, and 23, after 4th-sta check
%       %%ideally, when templates and data are not normalized, and there is no particular noise or any
%       %%other factors causing the amplitude scaling between temp and data and each station to be
%       %%vastly different, then for each deconvolved source, the direct impulse amp should be
%       %%~identical at all stations, ie., the ratio between station pairs should be ~1
      srcampr4th = [impindepst(:,2)./impindepst(:,4) impindepst(:,2)./impindepst(:,6) ...
                 impindepst(:,4)./impindepst(:,6) impindepst(:,2)./impindepst(:,17)];
               
      psrcamprs4th = [impindepst(:,2)*max(greenf(:,1))./(impindepst(:,4)*max(greenf(:,2))) ...
                   impindepst(:,2)*max(greenf(:,1))./(impindepst(:,6)*max(greenf(:,3))) ...
                   impindepst(:,4)*max(greenf(:,2))./(impindepst(:,6)*max(greenf(:,3))) ...
                   impindepst(:,2)*max(greenf(:,1))./(impindepst(:,17)*max(greenf(:,7)))];
      psrcamps4th = [impindepst(:,2)*max(greenf(:,1)) impindepst(:,4)*max(greenf(:,2)) ...
                  impindepst(:,6)*max(greenf(:,3)) impindepst(:,17)*max(greenf(:,7))];
      
      nsrcamprs4th = [impindepst(:,2)*min(greenf(:,1))./(impindepst(:,4)*min(greenf(:,2))) ...
                   impindepst(:,2)*min(greenf(:,1))./(impindepst(:,6)*min(greenf(:,3))) ...
                   impindepst(:,4)*min(greenf(:,2))./(impindepst(:,6)*min(greenf(:,3))) ...
                   impindepst(:,2)*min(greenf(:,1))./(impindepst(:,17)*min(greenf(:,7)))];
      nsrcamps4th = [impindepst(:,2)*min(greenf(:,1)) impindepst(:,4)*min(greenf(:,2))...
                  impindepst(:,6)*min(greenf(:,3)) impindepst(:,17)*min(greenf(:,7))];
      
      msrcampr4th(k,:) = median(log10(srcampr4th), 1);
      madsrcampr4th(k,:) = mad(log10(srcampr4th), 1, 1);
      mpsrcamprs4th(k,:) = median(log10(psrcamprs4th), 1);
      madpsrcamprs4th(k,:) = mad(log10(psrcamprs4th), 1, 1);
      mnsrcamprs4th(k,:) = median(log10(nsrcamprs4th), 1);
      madnsrcamprs4th(k,:) = mad(log10(nsrcamprs4th), 1, 1);
      nsrc(k,1) = size(srcampr,1);

      %what is the deviation of amp ratio from the median for each source?
      lndevsrcampr4th = srcampr4th-median(srcampr4th, 1); % in linear scale
      lgdevsrcampr4th = log10(srcampr4th)-median(log10(srcampr4th), 1); % in log scale, note that log2+log5=log10, so this means a ratio
      
      %%%what are the corresponding RCC at each source
      rccpairsrc4th = [];
      rccpairsrc4th(:,1) = rccpaircat(round(mean(impindepst(:,[1 3]),2)),1);
      rccpairsrc4th(:,2) = rccpaircat(round(mean(impindepst(:,[1 5]),2)),2);
      rccpairsrc4th(:,3) = rccpaircat(round(mean(impindepst(:,[3 5]),2)),3);
      
      %use the concatenated rcc at the average arrival time of each source
      rcccatsrc4th = [];
      rcccatsrc4th(:,1) = rcccat(round(mean(impindepst(:,[1 3 5]),2)));
      rcccatsrc4th(:,2) = rcc1icat(impindepst(:,9+(trust4th-4)*2+1),trust4th-3);
        
%       keyboard

      %% source amp ratio vs. RCC
%       f.fig=figure;
%       f.fig.Renderer='Painters';
%       for jj = 1:4
%         subplot(2,4,jj); hold on; box on; grid on;
%         if jj<4
%           scatter(log10(srcampr(:,jj)),rccpairsrc(:,jj),50,'MarkerFaceColor','k','MarkerEdgeColor',...
%             'none','MarkerFaceAlpha',.2);
%         else
%           scatter(log10(srcampr(:,jj)),rcccatsrc(:,2),50,'MarkerFaceColor','k','MarkerEdgeColor',...
%             'none','MarkerFaceAlpha',.2);
%         end
%         axis([-1 1 -1 1]);
%         plot([log10(median(srcampr(:,jj))) log10(median(srcampr(:,jj)))], [-1 1], 'r--','linew',1);
%         if jj == 1
%           ylabel('RCC between the same pair');
%         end
%         longticks(gca,2);
%       end
%       for jj = 1:4
%         subplot(2,4,4+jj); hold on; box on; grid on;
%         if jj<4
%           scatter(log10(srcampr(:,jj)),rcccatsrc(:,1),50,'MarkerFaceColor','k','MarkerEdgeColor','none',...
%             'MarkerFaceAlpha',.2);
%         else
%           scatter(log10(srcampr(:,jj)),rcccatsrc(:,2),50,'MarkerFaceColor','k','MarkerEdgeColor',...
%             'none','MarkerFaceAlpha',.2);
%         end
%         axis([-1 1 -1 1]);
%         plot([log10(median(srcampr(:,jj))) log10(median(srcampr(:,jj)))], [-1 1], 'r--','linew',1);
%         if jj == 1
%           xlabel('log_{10}{Amp ratio 1/2}');
%         elseif i==2
%           xlabel('log_{10}{Amp ratio 1/3}');
%         elseif i==3
%           xlabel('log_{10}{Amp ratio 2/3}');
%         else
%           xlabel('log_{10}{Amp ratio 1/4}');
%         end

%         if jj == 1
%           ylabel('Mean RCC');
%         end
%         longticks(gca,2);
%       end      
%       
%       end
%       
%       keyboard
      %% deviation of source amp ratio from some median vs. RCC; shift the plot above to center
%       f.fig=figure;
%       f.fig.Renderer='Painters';
%       for jj = 1:3
%         subplot(2,3,jj); hold on; box on; grid on;
%         scatter(lgdevsrcampr(:,jj),rccpairsrc(:,jj),50,'MarkerFaceColor','k','MarkerEdgeColor',...
%           'none','MarkerFaceAlpha',.2);
%         axis([-1 1 -1 1]);
%         if jj == 1
%           ylabel('RCC between the same pair');
%         end
%         longticks(gca,2);
%       end
%       for jj = 1:3
%         subplot(2,3,3+jj); hold on; box on; grid on;
%         scatter(lgdevsrcampr(:,jj),rcccatsrc,50,'MarkerFaceColor','k','MarkerEdgeColor','none',...
%           'MarkerFaceAlpha',.2);
%         axis([-1 1 -1 1]);
%         if jj == 1
%           xlabel('deviation from median amp ratio 1/2 (log)');
%         elseif jj == 2
%           xlabel('deviation from median amp ratio 1/3 (log)');
%         else
%           xlabel('deviation from median amp ratio 2/3 (log)');
%         end
%         if jj == 1
%           ylabel('Mean RCC of the 2 best pairs');
%         end
%         longticks(gca,2);
%       end      
      
      %% histograms of source amp
%       figure;
%       for ista = 1: nsta
%         subplot(3,3,ista); hold on; box on; grid on; ax=gca;
%         histogram(log10(srcampr(:,ista))); 
%         plot([log10(msrcampr(k,ista)) log10(msrcampr(k,ista))],ax.YLim,'r--');      
% %         errorbar(mclppkhtwfr(i),mpsrcamprs(i),madpsrcamprs(i),madpsrcamprs(i),...
% %           madclppkhtwfr(i),madclppkhtwfr(i),'color',[.5 .5 .5],'linewidth',0.8,'CapSize',5);
%         text(0.05,0.9,sprintf('med=%.2f',msrcampr(k,ista)),'Units','normalized');
%         text(0.9,0.9,'imp','Units','normalized','HorizontalAlignment','right');
%         ylabel('Counts');
%         if ista == 1
%           xlabel('log_{10}{Amp ratio 1/2}');
%         elseif ista == 2
%           xlabel('log_{10}{Amp ratio 1/3}');
%         else
%           xlabel('log_{10}{Amp ratio 2/3}');
%         end
%         xlim([-1 1]);
%         
%         subplot(3,3,3+ista); hold on; box on; grid on; ax=gca;
%         histogram(log10(psrcamprs(:,ista))); 
%         plot([log10(mpsrcamprs(k,ista)) log10(mpsrcamprs(k,ista))],ax.YLim,'r--');      
%         text(0.05,0.9,sprintf('med=%.2f',mpsrcamprs(k,ista)),'Units','normalized');
%         text(0.9,0.9,'imp*temp max','Units','normalized','HorizontalAlignment','right');
%         ylabel('Counts');      
%         if ista == 1
%           xlabel('log_{10}{Scaled positive amp ratio 1/2}');
%         elseif ista == 2
%           xlabel('log_{10}{Scaled positive amp ratio 1/3}');
%         else
%           xlabel('log_{10}{Scaled positive amp ratio 2/3}');
%         end
%         xlim([-1 1]);
%         
%         subplot(3,3,6+ista); hold on; box on; grid on; ax=gca;
%         histogram(log10(nsrcamprats(:,ista))); 
%         plot([log10(mnsrcamprs(k,ista)) log10(mnsrcamprs(k,ista))],ax.YLim,'r--');      
%         text(0.05,0.9,sprintf('med=%.2f',mnsrcamprs(k,ista)),'Units','normalized');
%         text(0.9,0.9,'imp*temp min','Units','normalized','HorizontalAlignment','right');
%         ylabel('Counts');      
%         if ista == 1
%           xlabel('log_{10}{Scaled negative amp ratio 1/2}');
%         elseif ista == 2
%           xlabel('log_{10}{Scaled negative amp ratio 1/3}');
%         else
%           xlabel('log_{10}{Scaled negative amp ratio 2/3}');
%         end
%         xlim([-1 1]);
% 
%       end
%       
%       keyboard
      %% plot deconvolved positive peaks with waveform peaks
%       %%plot the positive peaks indicated by the grouped triplets, see if they indeed match the
%       %%peaks of the data at each station, could be helpful which triplets are minor
%       fpltremv = 1;   % plot or not the secondary sources that are removed as well
%       [f,ppk,ppkhgt] = plt_deconpk_sigpk(sigsta,ppkindep,indremove,fpltremv);
        
      %% Compare impulse amp with waveform peaks 
      %%%compare the deconvolved positive peaks and negative peaks with the closest waveform peaks, 
      %%%ideally they should very similar in amplitude 
      ppkisave = [ppkindep(:,1:6) ppkindep(:,10:17)];  %deconvolved pos peaks of saved sources
      ppkisave = sortrows(ppkisave,1);
      npkisave = [npkindep(:,1:6) npkindep(:,10:17)];  %deconvolved neg peaks of saved sources
      npkisave = sortrows(npkisave,1);
      zcrsisave = [impindep(:,1:6) impindep(:,10:17)]; %deconvolved zero-crossings of saved sources
      zcrsisave = sortrows(zcrsisave,1);
      [f1,f2,clppk,clnpk] = plt_deconpk_sigpk_comp(sigsta,zcrsisave,ppkisave,npkisave,greenf);
      close(f1.fig);
      close(f2.fig);
      clppkhtwf = clppk.clppkhtwf;  %waveform peak height
      clppkhtwf = [clppkhtwf(:,1:3) clppkhtwf(:,end)];  %I only want KLNB
%       mclppkhtwf(k, :) = median(clppkhtwf, 1);
%       madclppkhtwf(k, :) = mad(clppkhtwf,1, 1);
      clppkwf = clppk.clppkwf;  %also store the waveform peak separation
      if size(clppkwf,1)>1
        clppkwfsep = diff([clppkwf(:,1:3) clppkwf(:,end)]);  %I only want KLNB
      else
        clppkwfsep = [];
      end
      
      clnpkhtwf = clnpk.clnpkhtwf;
      clnpkhtwf = [clnpkhtwf(:,1:3) clnpkhtwf(:,end)];  %I only want KLNB
%       mclnpkhtwf(k, :) = median(clnpkhtwf, 1);
%       madclnpkhtwf(k, :) = mad(clnpkhtwf,1, 1);
      clnpkwf = clnpk.clnpkwf;
      if size(clnpkwf,1)>1
        clnpkwfsep = diff([clnpkwf(:,1:3) clnpkwf(:,end)]);  %I only want KLNB
      else
        clnpkwfsep = [];
      end
      
      %also get the waveform peak separation for the full trace, since the number of peaks is not
      %the same for each station, keep 
      for ista = [1 2 3 nsta]
        [~, pk] = findpeaks(sigsta(:,ista));
        ppkwfsepmed4th(k,ista) = median(diff(pk));
        ppkwfsepmod4th(k,ista) = mode(diff(pk));
        [~, pk] = findpeaks(-sigsta(:,ista));
        npkwfsepmed4th(k,ista) = median(diff(pk));
        npkwfsepmod4th(k,ista) = mode(diff(pk));
      end
      
%       keyboard
      end
            
      %% How much area is slipping?
%       %%%Assume a Vs in the LVL, Vs * Tdura == source region size
%       Vs = 3;   % S-wave speed
%       Vprop = 0.8*Vs;   % rupture propagation speed
% %       td = 0.8*tdura;  % a rough median duration of the seismogram, can be 80% wide of the template
%       td = 40/sps;  % estimated from the bin where the most tsep falls in, 32-48-sample bin
%       srcsz = Vprop*td;
%       radi = srcsz/2;
%       xran = [-5 5];
%       yran = [-5 5];
%       cran = [0 lsig/sps];
%       f3.fig = figure;
%       fsty = 'shade';
% %       fsty = 'line';
%       ax3=gca;
%       ax3 = plt_srcsliparea(ax3,impindepst,xran,yran,cran,radi,sps,ftrans,fsty);
%       if isequal(fsty,'shade')
%         print(f3.fig,'-dpdf','-opengl','-r600','/home/chaosong/Pictures/abc.pdf');
%       elseif isequal(fsty,'line')
%         print(f3.fig,'-dpdf','-painters','/home/chaosong/Pictures/abc.pdf');
%       end
       
      %% summary the distance using diff 'nsep'
%       for nsep = 1: 10
%         dt = diffcustom(torisplst,nsep,'forward');
%         dist = sqrt(diffcustom(implocst(:,1),nsep,'forward').^2 + ...
%           diffcustom(implocst(:,2),nsep,'forward').^2 );
%         
%         med05sdist(nsep) = median(dist(dt/sps<=0.5));
%         med10sdist(nsep) = median(dist(dt/sps<=1));
%         med20sdist(nsep) = median(dist(dt/sps<=2));
%       end
%       figure
%       plot(1:10,med05sdist,'o-'); hold on
%       plot(1:10,med10sdist,'^-');
%       plot(1:10,med20sdist,'d-');
%       grid on; box on
%       legend('Source pairs whose origin time diff. is <=0.5 s', ...
%         'Source pairs whose origin time diff. is <=1 s', ...
%         'Source pairs whose origin time diff. is <=2 s');
%       xlabel('N samples of separation');
%       ylabel('Median distance (km)')
      
      %% analyze offset 23
%       mapview = 'offset';
%       [f] = plt_srcoffset23(imploc,impindepst,torispl,indsort,mapview);
            
      %% what's the corresponding rcc at source arrival, associated with the offset
%       [f] = plt_srcrccwithoffset(impindepst,rccpaircat);
     
      %% linear regression of offset change to know how fast the centroid of sources are moving
%       [f,fitobj12,fitobj13,dyfit12,dyfit13] = plt_deconlfit(impindepst,torispl);
      
      %% plot amp ratio 12 and 13, and 23
%       imploc = off2space002(impindepst(:,7:8),sps,ftrans,0);
%       mapview = 'offset';
%       [f,srcamprat,srcampratscl] = plt_srcampratio(imploc,impindepst,greenf,sps,mapview);
%       figure; 
%       subplot(211);
%       histogram(srcamprat(:,3),'binw',0.1); 
%       msrcampr(k) = median(srcamprat(:,3))
%       subplot(212);
%       histogram(srcampratscl(:,3),'binw',0.1); 
%       msrcamprs(k) = median(srcampratscl(:,3))
%       xlabel('Amp ratio 2/3');
%       ylabel('Counts');
      
      %% what is the distribution of 4-s detections inside the same window?
%       hfbsto = hfdayo(tmaxo>=tstbuf & tmaxo<=tedbuf, :);
%       hfbsti = hfdayi(tmaxi>=tstbuf & tmaxi<=tedbuf, :);
%       hfbst = sortrows([hfbsti; hfbsto], [daycol, seccol]);
%       colnum = [daycol seccol 1 2];
%       rangeplt = rangetemp(j,:);
%       figure
%       ax = gca;
%       [ax,hf] = plt_detections_inmap(ax,hfbst,[rangetemp(j,1) tstbuf tedbuf],[tstbuf tedbuf],colnum);

      %% final prediction via convolution between grouped impulses and template at each station 
%       [predgrp,resgrp,predgrpl,resgrpl]=predsig_conv_imptemp(sigsta,optdat,impindepst,...
%         greenf,zcrosses,overshoot,stas,1);
%       
%       %%%obtain an concatenated RCC for the residual as well, you can run this part after removing
%       %%%the secondary sources, but to justify the use of coda in templates if sources are widely-
%       %%%distributed, you can run it right after grouping. You should see the rcc of the residual is
%       %%%pretty low (eg., averaged to 0) compared to signal
%       %%%after grouping
%       irccrcat = [];   % concatenated indices of RCC for residual
%       rccrcat = [];  % concatenated RCC for residual
%       reslwin = windows-windows(1,1)+1+overshoot; % segmentation of windows for residual
%       for iwin = 1: nwin
%         ist = reslwin(iwin,1);
%         ied = reslwin(iwin,2);
%         resdat = [];
%         resdat(:, 1) = resgrpl(ist: ied, 1); % sta 1
%         resdat(:, 2) = resgrpl(ist-off1iw(iwin,2): ied-off1iw(iwin,2), 2); % sta 2
%         resdat(:, 3) = STAopt(ist-off1iw(iwin,3): ied-off1iw(iwin,3), 3); % sta 3
%         %compute running CC between 3 stations
%         [irccr,rccr12] = RunningCC(resdat(:,1), resdat(:,2), rccmwlen);
%         [~,rccr13] = RunningCC(resdat(:,1), resdat(:,3), rccmwlen);
%         [~,rccr23] = RunningCC(resdat(:,2), resdat(:,3), rccmwlen);
%         rccr = (rccr12+rccr13+rccr23)/3;  % use averge
%         irccr = irccr + reslwin(iwin,1) - reslwin(1,1);   %convert to global index
%         irccrcat = [irccrcat; irccr];
%         rccrcat = [rccrcat; rccr];
%       end
      % %without aligning to the best for each subwin, the whole-win rcc of residuals
      % %compute running CC between 3 stations
      % [irccr,rcc12,rcc13,rcc23] = RunningCC3sta(resgrpl,rccmwlen);
      % irccr = irccr-overshoot;
      % rccr = (rcc12+rcc13)/2;     
% 
%       %%% what about the orthogonal component?
%       %CC between sig and wlet at opt
%       nfft = lsig;
%       coef = [];
%       lag = [];
%       for ista = 1:3
%         [coef(:,ista), lag(:,ista)] = xcorr(sigsta(:,ista), greenf(:,ista), nfft, 'none'); % unnormalized master raw CC
%       end
%       %CC between sig and wlet at ort
%       coefort = [];
%       lagort = [];
%       for ista = 1:3
%         [coefort(:,ista), lagort(:,ista)] = xcorr(sigstaort(:,ista), greenfort(:,ista), nfft, 'none');
%       end
%       %%%for stas 1, 2, 3, cc between opt and ort
%       maxlag = 2*sps;
%       for ista = 1:3
%         %CC opt with ort
%         [cccoef, lagcoef] = xcorr(detrend(coef(:,ista)), detrend(coefort(:,ista)), maxlag, 'coeff');
%         [ccboo(k,ista), mind] = max(cccoef);
%         lagboo(k,ista) = lagcoef(mind);
%       end
% 
%       predgrport = zeros(size(sigstaort));  % predefine the prediction array
%       impfull = zeros(size(sigsta));  % a full array of impulses containing zeros
%       for ista = 1: nsta
%         impfull(impindepst(:,(ista-1)*2+1), ista) = impindepst(:,ista*2); % the non-zero index has the amp 
%         predtmp = conv(greenfort(:,ista), impfull(:,ista), 'full');
%         twlet = zcrosses(ista)*dt;
%         itwlet = round(twlet/dt);
%         predgrport(:,ista) = predtmp(itwlet:lsig+itwlet-1);  % cut accordingly
%         
%       end
%       resgrport = sigstaort - predgrport;  % the residual, different between signal and prediction
%       
%       figure
%       subplot(221)
%       hold on
%       box on; grid on
%       plot(sigsta(:,1),'r');
%       plot(sigsta(:,2),'b');
%       plot(sigsta(:,3),'k');
%       ym = max(abs(sigsta(:)));
%       yran=1.2*[-ym ym];
%       ylim(yran);
%       text(0.95,0.9,sprintf('%.2f; %.2f; %.2f',norm(sigsta(:,1)),norm(sigsta(:,2)),...
%         norm(sigsta(:,3))),'Units','normalized','HorizontalAlignment','right');
%       subplot(222)
%       hold on
%       box on; grid on
%       plot(resgrp(:,1),'r');
%       plot(resgrp(:,2),'b');
%       plot(resgrp(:,3),'k');
%       ylim(yran);
%       text(0.95,0.9,sprintf('%.2f; %.2f; %.2f',norm(resgrp(:,1)),norm(resgrp(:,2)),...
%         norm(resgrp(:,3))),'Units','normalized','HorizontalAlignment','right');
%       text(0.95,0.8,sprintf('%.1f%%; %.1f%%; %.1f%%',...
%         (norm(sigsta(:,1))-norm(resgrp(:,1)))/norm(sigsta(:,1))*100,...
%         (norm(sigsta(:,2))-norm(resgrp(:,2)))/norm(sigsta(:,2))*100,...
%         (norm(sigsta(:,3))-norm(resgrp(:,3)))/norm(sigsta(:,3))*100),...
%         'Units','normalized','HorizontalAlignment','right');
%       text(0.95,0.1,sprintf('%.2f; %.2f; %.2f',ccboo(:,1),ccboo(:,2),...
%         ccboo(:,3)),'Units','normalized','HorizontalAlignment','right');
%       subplot(223)
%       hold on
%       box on; grid on
%       plot(sigstaort(:,1),'r');
%       plot(sigstaort(:,2),'b');
%       plot(sigstaort(:,3),'k');
%       ym = max(abs(sigstaort(:)));
%       yran=1.2*[-ym ym];
%       ylim(yran);
%       text(0.95,0.9,sprintf('%.2f; %.2f; %.2f',norm(sigstaort(:,1)),norm(sigstaort(:,2)),...
%         norm(sigstaort(:,3))),'Units','normalized','HorizontalAlignment','right');
%       subplot(224)
%       hold on
%       box on; grid on
%       plot(resgrport(:,1),'r');
%       plot(resgrport(:,2),'b');
%       plot(resgrport(:,3),'k');
%       ylim(yran);
%       text(0.95,0.9,sprintf('%.2f; %.2f; %.2f',norm(resgrport(:,1)),norm(resgrport(:,2)),...
%         norm(resgrport(:,3))),'Units','normalized','HorizontalAlignment','right');
%       text(0.95,0.8,sprintf('%.1f%%; %.1f%%; %.1f%%',...
%         (norm(sigstaort(:,1))-norm(resgrport(:,1)))/norm(sigstaort(:,1))*100,...
%         (norm(sigstaort(:,2))-norm(resgrport(:,2)))/norm(sigstaort(:,2))*100,...
%         (norm(sigstaort(:,3))-norm(resgrport(:,3)))/norm(sigstaort(:,3))*100),...
%         'Units','normalized','HorizontalAlignment','right');



%     end
%     
%   end
%   
% end

end

% keyboard

%% if 'pltflag' is on, then summary plots for each choice of inputs would be made 
if pltflag && ~isempty(impindepst)
%   %%%the direct deconvolved pos/neg source peak ratio between all station pairs, for each
%   %%%burst win separately
%   f1 = initfig(16,5,1,size(msrcampr,2));
%   plt_deconpk_rat4th(f1,msrcampr,madsrcampr,nsrc,'k');
% 
%   %%%combine the direct/scaled deconvolved pos/neg source peak ratio between all station pairs of
%   %%%all burst wins, and summarize into one histogram
%   f2 = initfig(16,9,2,size(srcamprall,2));
%   plt_deconpk_rat_comb4th(f2,srcamprall,impindepstall,'k');
% 
%   %%%deviation of source amp ratio from some median vs. RCC
%   f3 = initfig(15,7,2,size(lgdevsrcamprall,2)); %initialize fig
%   plt_deconpk_ratdevvsrcc4th(f3,lgdevsrcamprall,rccpairsrcall,rcccatsrcall,'k');

  %%%diff between predicted arrival and selected peak at 4th sta
  f4 = initfig(5,5,1,size(pred4difftrall,2)); %initialize fig
  plt_errorof4thtarvlpred(f4,pred4difftrall,offmax,'k');

  %%%preserved sources' amp ratio between 4th and 1st stas
  f5 = initfig(12,5,1,3); %initialize fig
  plt_deconpk_rat14(f5,impindep4thall,srcamprall,'k');

%   %%%histogram of RCC itslef
%   figure
%   for i = 1:3
%     subplot(1,4,i)
%     hold on; box on; grid on; ax=gca;
%     histogram(rccpairsrcall(:,i));
%     plot([median(rccpairsrcall(:,i)) median(rccpairsrcall(:,i))],ax.YLim,'r--','linew',1);  
%     text(0.05,0.9,sprintf('med=%.2f',median(rccpairsrcall(:,i))),...
%       'Units','normalized','HorizontalAlignment','left');
%     if i ==1
%       ylabel('# of source');
%       xlabel('RCC_{12}');
%     elseif i ==2
%       xlabel('RCC_{13}');
%     else
%       xlabel('RCC_{23}');
%     end
%   end
%   subplot(1,4,4)
%   hold on; box on; grid on; ax=gca;
%   histogram(rcccatsrcall);
%   plot([median(rcccatsrcall) median(rcccatsrcall)],ax.YLim,'r--','linew',1);
%   text(0.05,0.9,sprintf('med=%.2f',median(rcccatsrcall)),...
%     'Units','normalized','HorizontalAlignment','left');
%   xlabel('Mean RCC of 2 best pairs');

end


% keyboard






