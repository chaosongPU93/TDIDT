% deconvbursts002_ref_4s_exp_random.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on 'deconvbursts002_ref_4s_exp', but this script tries to simulate what 
% the current deconvolution algorithm would behave if the records are 
% random noise. 
%
% --Preivously in 'deconvbursts002_ref_4s_exp', we try to mimic the
%   noise by shifting records up to a few seconds, but looks like there is
%   still some implicit coherence between stations that is hard to explain.
% --This code would simulate noise by using the same amplitude spectrum from
%   the real records and random phase spectrum with a uniform distribution 
%   that has the same range as the records. This is done by 'fft' and 'ifft',
%   and the signal in freq domain is a complex expressed by 
%   xf = abs(xf)*exp(i*angle(xf)). The amplitude spectrum is abs(xf)/nfft, 
%   and the phase spectrum is angle(xf).
%
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/08/09
% Last modified date:   2022/08/09
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
  'SILB '];     % determine the trio and order, here the 1st sta is PGC
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


%% load other catalogs
%%% load the LFE catalog of Michael Bostock, inside and outside the rectangle in 'locinterp002_4s.m'
%obtain the location of fam 002, lon0 and lat0
ftrans = 'interpchao';
loc0 = off2space002([0 0],sps*iup,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
%format: [fam yyyy mm dd sec dx dy lon lat dep magnitude number-of-stations], 12 cols
%time should point to the peak, zero-crossing?
bostcat = ReformBostock(loc0(3),loc0(4),0);

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
STA = ccstack;

ccstackort = [];
for ista = 1: nsta
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
        num2str(templensec), 's_', 'BBCCS_', 'ort_Nof_Non_Chao_catnew');
    ccstackort(:,ista) = load(fname);
end
STAort = ccstackort;

%flag of normalization
normflg = 0;

% %plot the raw templates, not filtered, not best aligned
% figure
% subplot(211)
% hold on
% plot(STA(:,1),'r')
% plot(STA(:,2),'b')
% plot(STA(:,3),'k')
% subplot(212)
% hold on
% plot(STAort(:,1),'r')
% plot(STAort(:,2),'b')
% plot(STAort(:,3),'k')

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
    seg = STA(iminses(ista): imaxses(ista),ista);  % for zero-crossing timing, only use the main station
    [~,zcrosses(ista)] = min(abs(seg));
    zcrosses(ista) = zcrosses(ista)-1+iminses(ista);  % convert to global index
end
%now we want to cut a segment around the zero-crossing at each station
sampbef=6*sps;
sampaft=10*sps;
is=zcrosses-sampbef;
ie=zcrosses+sampaft;
for ista=1:nsta
    STAtmp(:,ista)=STA(is(ista):ie(ista),ista);  % this means templates are 'aligned' at zero-crossings
    STAtmport(:,ista)=STAort(is(ista):ie(ista),ista); 
end
%x-correlation independently between each station pair 
mshiftadd=10*sps/40;
tempxc(:,1)=xcorr(STAtmp(:,2),STAtmp(:,3),mshiftadd,'coeff');
tempxc(:,2)=xcorr(STAtmp(:,1),STAtmp(:,2),mshiftadd,'coeff'); %shift STAtmp(3,:) to right for positive values
tempxc(:,3)=xcorr(STAtmp(:,1),STAtmp(:,3),mshiftadd,'coeff'); %shift STAtmp(2,:) to right for positive values
[~,imax]=max(tempxc,[],1);
imax=imax-(mshiftadd+1); %This would produce a slightly different shift, if filtered seisms were used.
imax(2)-imax(3)+imax(1);   %enclosed if it equals to 0
for ista=2:nsta
    STAtmp(mshiftadd+1:end-(mshiftadd+1),ista)=STAtmp(mshiftadd+1-imax(ista):end-(mshiftadd+1)-imax(ista),ista);
    STAtmport(mshiftadd+1:end-(mshiftadd+1),ista)=STAtmport(mshiftadd+1-imax(ista):end-(mshiftadd+1)-imax(ista),ista);
end
%normalization
if normflg 
  for ista=1:nsta
      STAtmp(:,ista)=STAtmp(:,ista)/spread(ista); % now templates are 'aligned' indeoendently by x-corr wrt. sta 1
      STAtmport(:,ista)=STAtmport(:,ista)/spread(ista);
  end
end
% figure
% subplot(211)
% hold on
% plot(STAtmp(:,1),'r')
% plot(STAtmp(:,2),'b')
% plot(STAtmp(:,3),'k')
% subplot(212)
% hold on
% plot(STAtmport(:,1),'r')
% plot(STAtmport(:,2),'b')
% plot(STAtmport(:,3),'k')
%%%The above aligns the templates by x-correlation

%%%detrend, taper and bandpass templates
tmpwlet = STAtmp; % no bandpass
tmpwletf = STAtmp;  % bandpassed version
fractap = sps/size(tmpwlet,1);
tmpwletort = STAtmport; % no bandpass
tmpwletfort = STAtmport;  % bandpassed version
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
  tmpwletf(:,ista) = Bandpass(tmpwlet(:,ista), sps, lowlet, hiwlet, 2, 2, 'butter');
  %detrend again for caution
  tmpwletf(:,ista)=detrend(tmpwletf(:,ista));
  
  %same process for orthogonal
  tmpwletort(:,ista) = detrend(tmpwletort(:,ista));
  tmpwletort(:,ista) = w.* tmpwletort(:,ista);
  tmpwletort(:,ista)=detrend(tmpwletort(:,ista));
  tmpwletfort(:,ista) = Bandpass(tmpwletort(:,ista), sps, lowlet, hiwlet, 2, 2, 'butter');
  tmpwletfort(:,ista)=detrend(tmpwletfort(:,ista));
end

%%%constrained CC, so that only 2 offsets are independent
ccmid = round(size(tmpwletf,1)/2);
ccwlen = 10*sps;
loffmax = 5*sps/40;
ccmin = 0.01;  % depending on the length of trace, cc could be very low
iup = 1;    % times of upsampling
[off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(tmpwletf',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
offwlet1i(1) = 0;
offwlet1i(2) = round(off12con);
offwlet1i(3) = round(off13con);

%%%automatically find the rough zero-crossing time, whose abs. value is closest to 0, whether + or -
[~,imin] = min(tmpwletf(:,1));
[~,imax] = max(tmpwletf(:,1));
[~,zcsta1] = min(abs(tmpwletf(imin:imax,1)));
zcsta1 = zcsta1+imin-1;
greenlen = pow2(9)*sps/40;
green = zeros(greenlen,nsta); % no bandpass
greenf = zeros(greenlen,nsta);  % bandpassed version
ppeaks = zeros(nsta,1);
greenort = zeros(greenlen,nsta); % no bandpass
greenfort = zeros(greenlen,nsta);  % bandpassed version
for ista = 1: nsta
  %cut according to the zero-crossing and the time shift from the constrained CC
  green(:,ista) = tmpwlet(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenf(:,ista) = tmpwletf(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  %detrend again for caution
  green(:,ista)=detrend(green(:,ista));
  greenf(:,ista)=detrend(greenf(:,ista));
  if normflg
    %normalize by max amp
    green(:,ista)=green(:,ista)/max(abs(green(:,ista)));    % normalize
    greenf(:,ista)=greenf(:,ista)/max(abs(green(:,ista)));    % normalize
  end
  
  %same process for orthogonal
  greenort(:,ista) = tmpwletort(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenfort(:,ista) = tmpwletfort(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenort(:,ista)=detrend(greenort(:,ista));
  greenfort(:,ista)=detrend(greenfort(:,ista));
  if normflg
    greenort(:,ista)=greenort(:,ista)/max(abs(green(:,ista)));    % normalize
    greenfort(:,ista)=greenfort(:,ista)/max(abs(green(:,ista)));    % normalize
  end
  
  %re-find the zero-crossing as the template length has changed
  [~,imin] = min(greenf(:,ista));
  [~,imax] = max(greenf(:,ista));
  [~,zcrosses(ista)] = min(abs(greenf(imin:imax,ista)));
  zcrosses(ista) = zcrosses(ista)+imin-1;
  ppeaks(ista) = imax;
end
%the following is just a check, because now the templates must be best aligned 
ccmid = round(size(greenf,1)/2);
ccwlen = 4*sps;
loffmax = 5*sps/40;
ccmin = 0.01;  % depending on the length of trace, cc could be very low
iup = 1;    % times of upsampling
[off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(greenf',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
if ~(off12con==0 && off13con==0)
  disp('Filtered templates are NOT best aligned');
end

amprat(1,:) = minmax(greenf(:,1)')./minmax(greenf(:,2)');	% amp ratio between max at sta 3 and 2 or min
amprat(2,:) = minmax(greenf(:,1)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min  
amprat(3,:) = minmax(greenf(:,2)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min  
spread = range(greenf);   % range of the amp of template

%%%plot the unfiltered and filtered templates
% plt_templates(green,greenf,greenort,greenfort,lowlet,hiwlet,sps);


%% prepare the signal and noise windows
% dates and ets
%%% NOTE: 'dates' is the dates that the tremor is active at the region of interest, so that the
%%% waveform at the stations on these dates you are looking at mainly comes from the sources at the
%%% region of interest. We want to see if there is a noticable change in spectra during the burst
%%% windows on these dates
dates = unique(trange(:,1));
years = unique(floor(dates/1000));
nets = length(years);

%filtering passband for reading data, confirmed by 'spectrabursts002_4s.m'
hisig=6.3; % this will give a similar spectral shape between template and signal
losig=1.8;

%moving window length in samples for running CC, envelope, etc.
mwlen=sps/2;
% mwlen=sps;

k = 0;  % counting the burst windows

off1iwk = cell(size(trange,1),1);  % the best alignment between sta 2, 3 wrt 1 for all subwins and all burst wins
off1i = zeros(size(trange,1),3);  % single best alignment between sta 2, 3 wrt 1 for entire win
ccali = zeros(size(trange,1),1);  % CC value using the best alignment
subwsec = zeros(size(trange,1),1);  % subwin length in sec used in practice
subwseclfit = zeros(size(trange,1),1);  % subwin length from linear fitting, ie, time for 1-sample offset change

for iets = 3: nets
  % dates in each ets
  year = years(iets);
  datesets = dates(floor(dates/1000)==year);
    
  for i = 2: length(datesets)
    
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
    disp(prename);
    [STAopt,STAort,~,fileflag] = rd_daily_bpdata(year,jday,prename,stas,PERMSTA,POLSTA,...
      PERMROTS,POLROTS,sps,losig,hisig,npo,npa,[],[],[],[]);

    if fileflag == 0    % means there are missing files
      fprintf('Day %s / %s will be omitted because of missing files. \n', yr, JDAY);
      continue    % continue to the next day
    end
    
%     %read vertical components too, as we want to have a sense of the SNR at the orthogonal and vertical
%     %components
%     [STAvert,~,~,fileflag] = rd_daily_bpdata(year,jday,prename,stas,PERMSTA,POLSTA,...
%       PERMROTS,POLROTS,sps,losig,hisig,npo,npa,[],[],[],[],'Z');

    
%     keyboard
%%
    for j = 14: size(rangetemp,1)  
%       close all
      k = k+1;  
      disp(k);

      tmaxi = hfdayi(:, seccol); % starting time of max power rate of half sec inside the ellipse
      tmaxo = hfdayo(:, seccol); % starting time of max power rate of half sec outside the ellipse  
      tcnti = hfdayi(:, 15);  % the center of detecting win is the 15th col
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
                          
      %max allowable shift in best alignment
%       msftaddm = (round(max(abs([off12ran off13ran])))+1)*sps/40;  %+1 for safety
%       msftaddm = sps+1;  %+1 for safety
      msftaddm = 1.5*sps+1;  %+1 for safety
%       msftaddm = round(sps/8);    % maximum allowed shift between 2 traces

      %have some overshoot, so that the resulted rcc would have the same length as the signal
      overshoot = mwlen/2;
%       overshoot = 0;
      
      %FLAG to simulate the behavior of noise
      noiseflag = 1;

%       seedmat = randi(1000,200,1);
%       for iii = 1: length(seedmat)
%       seed = seedmat(iii);

      %chop a record segment
      optseg = STAopt(max(floor(tstbuf*sps+1-overshoot-msftaddm),1): ...
        min(floor(tedbuf*sps+overshoot+msftaddm),86400*sps), :); % sta 1
      ortseg = STAort(max(floor(tstbuf*sps+1-overshoot-msftaddm),1): ...
        min(floor(tedbuf*sps+overshoot+msftaddm),86400*sps), :);

      if noiseflag
        %obtain the amp and phase spectra of records via fft
        nfft = size(optseg,1); % number of points in fft
        [xf,ft,amp,pha] = fftspectrum(optseg(:,2:4), nfft, sps,'twosided');

        %uniform, random phase with the same span [-pi,pi];
        mpharan = minmax(pha');
        seed = k;
        rng(seed);
        pharand = (rand(nfft,3)-0.5)*2*pi;  %make the phases span from -pi to pi

%         figure
%         subplot(131)
%         histogram(pharand(:,1));
%         subplot(132)
%         histogram(pharand(:,2));
%         subplot(133)
%         histogram(pharand(:,3));

        %construct record with the same amplitude but random phase
        xfrand = amp*nfft.*exp(1i*pharand);
        optseg(:,2:4) = real(ifft(xfrand,nfft));

%         figure
%         subplot(311)
%         plot(optseg(:,2));
%         subplot(312)
%         plot(optseg(:,3));
%         subplot(313)
%         plot(optseg(:,4));
%         
%         figure
%         subplot(131)
%         histogram(optseg(:,2)); [MUHAT,SIGMAHAT] = normfit(optseg(:,2))
%         subplot(132)
%         histogram(optseg(:,3)); [MUHAT,SIGMAHAT] = normfit(optseg(:,3))
%         subplot(133)
%         histogram(optseg(:,4)); [MUHAT,SIGMAHAT] = normfit(optseg(:,4))
        
        %%%for orthogonal components, with the same phase??
        [xf,ft,amp,pha,power,psd] = fftspectrum(ortseg(:,2:4), nfft, sps,'twosided');
%         pharand = (rand(nfft,3)-0.5)*2*pi;
        xfrand = amp*nfft.*exp(1i*pharand);
        ortseg(:,2:4) = real(ifft(xfrand,nfft));
                
      end

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
      %since the rcc would lose mwlen/2 at both ends, this results in overlapping of 'mwlen' in rcc 
      %across consecutive windows; if use 'ovlplen' of mwlen, then rcc has no overlapping at all
%       ovlplen = mwlen*2;
      ovlplen = mwlen;
      [windows] = movingwins(indst,inded,subwlen,ovlplen,0);
      %if some portions are not included
      if windows(end,2) < inded
        %make it an separate window if its length is at least half of that of the others
        if inded-(windows(end,1)+(subwlen-ovlplen))+1 >= subwlen/2
          windows = [windows; [windows(end,1)+(subwlen-ovlplen) inded]];
        %otherwise just combine the portion to the last window
        else
          windows(end,2) = inded;
        end
      end
      nwin =  size(windows,1);
      
      off1iw = zeros(nwin,3);  % the best alignment between sta2, sta3 wrt sta1 for each subwin
      
      ircccat = [];   % concatenated indices of RCC
      irccran = zeros(nwin,2);  % start and end indices (range) of RCC of all subwins
      rcccat = [];  % average concatenated RCC
      rccpaircat = [];  % concatenated RCC between each station pair, order is 12, 13, 23
      ccwpair = []; % 0-lag overall cc of each subwin, between each station pair, order is 12, 13, 23 
      
%       figure
%       hold on
%       box on
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
        optcc = optseg(isubwst: isubwed, 2:end);
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
        if noiseflag
          %ask how well can you align the noise, does not matter if the location is
          %the most reliable, but the resulting CC is the highest possible
          [off12con,off13con,ccali(k),iloopoff] = constrained_cc_loose(optcc',ccmid,...
            ccwlen,msftadd,ccmin,iup);
        else
          [off12con,off13con,ccali(k),iloopoff,loopoff] = constrained_cc_interp(optcc',ccmid,...
            ccwlen,msftadd,loffmax,ccmin,iup);
        end
        % if a better alignment cannot be achieved, use 0,0
        if off12con == msftadd+1 && off13con == msftadd+1
          off12con = 0;
          off13con = 0;
          fprintf('Tremor burst %d cannot be properly aligned, double-check needed \n',k);
        end
        off1iw(iwin,1) = 0;
        off1iw(iwin,2) = round(off12con);
        off1iw(iwin,3) = round(off13con);
        
        %Align records
        optdat = [];  % win +/-3 s, segment of interest,first 1s will be tapered 
        ortdat = [];
        optdat(:, 1:2) = optseg(isubwst: isubwed, 1:2); % sta 1
        ortdat(:, 1:2) = ortseg(isubwst: isubwed, 1:2);      
        optdat(:, 3) = optseg(isubwst-off1iw(iwin,2): isubwed-off1iw(iwin,2), 3); % sta 2
        ortdat(:, 3) = ortseg(isubwst-off1iw(iwin,2): isubwed-off1iw(iwin,2), 3);
        optdat(:, 4) = optseg(isubwst-off1iw(iwin,3): isubwed-off1iw(iwin,3), 4); % sta 3
        ortdat(:, 4) = ortseg(isubwst-off1iw(iwin,3): isubwed-off1iw(iwin,3), 4);
        
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
        [irccw,rccw12] = RunningCC(subw(:,1), subw(:,2), mwlen);
        [~,rccw13] = RunningCC(subw(:,1), subw(:,3), mwlen);
        [~,rccw23] = RunningCC(subw(:,2), subw(:,3), mwlen);
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
        
        ccw12 = xcorr(subw(:,1), subw(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
        ccw13 = xcorr(subw(:,1), subw(:,3),0,'normalized');
        ccw23 = xcorr(subw(:,2), subw(:,3),0,'normalized');
        ccwpair = [ccwpair; ccw12 ccw13 ccw23];
      end
      off1iwk{k} = off1iw;
      
      %if only use the mean RCC from the 2 pairs that have the highest overall CC
      [~,ind] = min(sum(ccwpair,1));
      rcccat = sum(rccpaircat(:,setdiff(1:3,ind)), 2) / 2;
      
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
      optcc = optseg(1+msftaddm: end-msftaddm, 2:end);
      msftadd = 0.5*sps;
      ccmid = ceil(size(optcc,1)/2);
      ccwlen = round(size(optcc,1)-2*(msftadd+1));
      ccmin = 0.01;  % depending on the length of trace, cc could be very low
      iup = 1;    % times of upsampling
      %for the whole win, ask how well can you align the noise, does not matter if the location is
      %the most reliable, but the resulting CC is the highest possible
      [off12con,off13con,ccali(k),iloopoff] = constrained_cc_loose(optcc',ccmid,...
        ccwlen,msftadd,ccmin,iup);
      % if a better alignment cannot be achieved, use 0,0
      if off12con == msftadd+1 && off13con == msftadd+1
        off12con = 0;
        off13con = 0;
        fprintf('Tremor burst %d cannot be properly aligned, double-check needed \n',k);
      end
      off1i(k,1) = 0;
      off1i(k,2) = round(off12con);
      off1i(k,3) = round(off13con);
      
      %%%Align and compute the RCC based on the entire win, and take that as the input signal!      
      optdat = [];  % win segment of interest
      ortdat = [];
      optdat(:, 1:2) = optseg(1+msftaddm: end-msftaddm, 1:2); % sta 1
      ortdat(:, 1:2) = ortseg(1+msftaddm: end-msftaddm, 1:2);
      optdat(:, 3) = optseg(1+msftaddm-off1i(k,2): end-msftaddm-off1i(k,2), 3); % sta 2
      ortdat(:, 3) = ortseg(1+msftaddm-off1i(k,2): end-msftaddm-off1i(k,2), 3);
      optdat(:, 4) = optseg(1+msftaddm-off1i(k,3): end-msftaddm-off1i(k,3), 4); % sta 3
      ortdat(:, 4) = ortseg(1+msftaddm-off1i(k,3): end-msftaddm-off1i(k,3), 4);

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
      [ircc,rcc12] = RunningCC(sigsta(:,1), sigsta(:,2), mwlen);
      [~,rcc13] = RunningCC(sigsta(:,1), sigsta(:,3), mwlen);
      [~,rcc23] = RunningCC(sigsta(:,2), sigsta(:,3), mwlen);
      ircc = ircc-overshoot;
      rcc = (rcc12+rcc13+rcc23)/3;
      sigsta = sigsta(overshoot+1:end-overshoot, :);  %excluding the overshoot
      
      %for ort. comp
      sigstaort = zeros(size(ortdat,1), nsta);
      for ista = 1:nsta
        tmp = ortdat(:,ista+1); %best aligned, filtered
        tmp = detrend(tmp);
        sigstaort(:,ista) = tmp;
      end
      [irccort,rcc12] = RunningCC(sigstaort(:,1), sigstaort(:,2), mwlen);
      [~,rcc13] = RunningCC(sigstaort(:,1), sigstaort(:,3), mwlen);
      [~,rcc23] = RunningCC(sigstaort(:,2), sigstaort(:,3), mwlen);
      irccort = irccort-overshoot;
      rccort = (rcc12+rcc13+rcc23)/3;
      sigstaort = sigstaort(overshoot+1:end-overshoot, :);  %excluding the overshoot
      
%       figure
%       hold on
%       box on
%       plot(ircc,rcc,'k-')
%       plot(ircccat,rcccat,'r-');
%       ylim([-1 1]);
            
      %%%finalize the signal, noise, and template (Green's function)
      sigdecon = [];
      pred = [];
      ampit = [];
      for ista = 1:nsta
        wlet = greenf(:,ista);  %template here is best aligned, tapered, linear trend removed, filtered
        lwlet = length(wlet);
        sig = sigsta(:,ista); %best aligned, filtered, tapered
        lsig = length(sig);
        noi = noidat(:,ista+1);
        noi = detrend(noi); %best aligned, filtered, linear trend removed
        
        dt = 1/sps;  % sampling interval
        twlet = zcrosses(ista)*dt;
        width = 2.5;  % width for Gaussian filter
        dres_min = 0.5;  % tolerance, percentage change in residual per iteration
        mfit_min = 5e-1;  % tolerance, norm of the residual/misfit
        tdura = 0.5;  % estimate from the broadband template from fam 002
        nit_max = round(1.5*1/tdura*(tlenbuf));  % max numer of iterations
        nimp_max = round(1/tdura*(tlenbuf));%a single peak is ~20 samples wide; maybe a little less (at 100 sps). ~0.4s, 1/0.4=2.5
        fpltit = 0;  % plot flag for each iteration
        fpltend = 0;  % plot flag for the final iteration
        fpltchk = 0; % plot flag for intermediate computations
        if noiseflag
          if ista == 1
%             fixthresh = 2.3280e-01; % if use 3-pair rcc and normalize templates, sta 1
%             fixthresh = 2.7255e-01; % if use 2-pair rcc and normalize templates
%             fixthresh = 1.2045e-01; % if use 3-pair rcc and do not normalize templates
            fixthresh = 1.4102e-01; % if use 2-pair rcc and do not normalize templates
          elseif ista == 2
%             fixthresh = 1.4667e-01; % if use 3-pair rcc and normalize templates, sta 2
%             fixthresh = 1.5914e-01; % if use 2-pair rcc and normalize templates
%             fixthresh = 5.8671e-02; % if use 3-pair rcc and do not normalize templates
            fixthresh = 6.3658e-02; % if use 2-pair rcc and do not normalize templates
          elseif ista == 3
%             fixthresh = 1.9222e-01; % if use 3-pair rcc and normalize templates, sta 3
%             fixthresh = 2.2233e-01; % if use 2-pair rcc and normalize templates
%             fixthresh = 7.4332e-02; % if use 3-pair rcc and do not normalize templates
            fixthresh = 8.5975e-02; % if use 2-pair rcc and do not normalize templates
          end
          [sigdecon(:,ista),pred(:,ista),res,dresit,mfitit,ampit{ista},nit,fighdl] = ...
            iterdecon_fixthresh(sig,wlet,rcccat,noi,fixthresh,dt,twlet,width,dres_min,mfit_min,nit_max,nimp_max,fpltit,fpltend,fpltchk);
        else
          [sigdecon(:,ista),pred(:,ista),res,dresit,mfitit,ampit{ista},nit,fighdl] = ...
            iterdecon(sig,wlet,rcccat,noi,dt,twlet,width,dres_min,mfit_min,nit_max,nimp_max,fpltit,fpltend,fpltchk);
%         [sigdecon(:,ista),pred(:,ista),res,dresit,mfitit,ampit,nit,fighdl] = ...
%           iterdecon_rcc(sig,wlet,rcc,medrcc,dt,twlet,width,dres_min,mfit_min,nit_max,nimp_max,fpltit,fpltend,fcheck);
        end
        
        if fpltend
          ax = fighdl{2}.ax(1);
          hold(ax,'on');
          text(ax,0.05,0.85,stas(ista,:),'unit','normalized');
          hold(ax,'off');
        end
        
        nit
        
      end
            
      %% Group nearest impulses from different stations into pairs, using moving searching range
      spsscale = sps/40;
      loff_max = 6*spsscale;  % maximum allowed shift is 4 samples at 40 Hz, it is also true from the synthetic distribution
      %note the output 'impindep' gives the arrival index of impulse at each station, after
      %alignment based upon the entire window 'off1i', and the last three cols are the arrival time
      %difference, NOT the true location yet! 
      refsta = 1;
%       [impindep,imppairf,indpair] = groupimptriplets_ref(sigdecon,irccran,rcccat,...
%         off1i(k,:),off1iw,loff_max,'wtamp',refsta);
      [impindep,imppairf,indpair] = groupimptripdecon_ref(sigdecon,ampit,irccran,rcccat,...
        off1i(k,:),off1iw,loff_max,refsta);
      
%       nsrc(iii) = size(impindep,1);
%       end
%     
%     figure
%     histogram(nsrc);
%     text(0.05,0.9,sprintf('med=%d',median(nsrc)),'Units','normalized');
%     xlabel('number of sources');
%     ylabel('counts');
%     
%     keyboard
      
%       %%%plot the individually deconvolved impulses and the grouping result
%       [f] = plt_groupimptriplets(sigdecon,impindep,stas,ircccat,rcccat);
      
      %note here 'impindepst' inherits the first 6 cols from 'impindep', but the last three cols 
      %are adjusted from arrival time difference to the true location offset accounting for the best
      %alignment upon each subwin that is also used in grouping!
      impindepst = sortrows(impindep,1);
      impindepst(:,7:8) = impindepst(:,7:8)+repmat([off1i(k,2) off1i(k,3)],size(impindepst,1),1); %account for prealignment

      %%%plot the scatter of offsets, accounting for prealignment offset, == true offset
      span = max(range(off1iw(:,2))+2*loff_max, range(off1iw(:,3))+2*loff_max);
      xran = [round(mean(minmax(off1iw(:,2)'))-span/2)-1, round(mean(minmax(off1iw(:,2)'))+span/2)+1];
      yran = [round(mean(minmax(off1iw(:,3)'))-span/2)-1, round(mean(minmax(off1iw(:,3)'))+span/2)+1];
      cran = [0 lsig];
      f1.fig = figure;
      f1.fig.Renderer = 'painters';
      ax1=gca;
      [ax1,torispl,mamp] = plt_decon_imp_scatter_ref(ax1,impindepst,xran,yran,cran,off1iw,loff_max,...
        sps,50,'mean','tori','comb');
      scatter(ax1,off1i(k,2),off1i(k,3),20,'ks','filled','MarkerEdgeColor','k');
      title(ax1,'Independent, grouped');
      
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

      
      %% check the difference by grouping using different stations as the reference station
      spsscale = sps/40;
      loff_max = 6*spsscale;  % maximum allowed shift is 4 samples at 40 Hz, it is also true from the synthetic distribution
      %note the output 'impindep' gives the arrival index of impulse at each station, after
      %alignment based upon the entire window 'off1i', and the last three cols are the arrival time
      %difference, NOT the true location yet! 
%       aa = groupimptripdecon_refOLD(sigdecon,ampit,irccran,rcccat,off1i(k,:),off1iw,loff_max);
      aa = groupimptriplets_ref(sigdecon,irccran,rcccat,off1i(k,:),off1iw,loff_max,'wtamp',1);

      for refsta = 1:3
%         bb = groupimptripdecon_ref(sigdecon,ampit,irccran,rcccat,off1i(k,:),off1iw,loff_max,refsta);
        bb = groupimptriplets_ref(sigdecon,irccran,rcccat,off1i(k,:),off1iw,loff_max,'wtamp',refsta);
        size(bb,1)
        indsame = find(ismember(bb,aa,'rows')); %which sources are the same?
        inddiff = find(~ismember(bb,aa,'rows'));
        perca = length(indsame)/size(aa,1)*100
        percb = length(indsame)/size(bb,1)*100
      end      


      %% Remove the small-amplitude, secondary triplets from the grouped result
      %convert the sources in terms of the arrival time of zero-crossing to positve peaks' indices
      pkindep = impindep;
      for ista = 1: nsta
        pkindep(:,(ista-1)*2+1) = pkindep(:,(ista-1)*2+1)+ppeaks(ista)-zcrosses(ista);
      end
      
      %are there any grouped sources are reversed in arrival order at different sta? 
      pkist = [];
      indpk = [];
      tsep = [];
      for ista = 1:nsta
        pki = pkindep(:,(ista-1)*2+1);
        [pkist(:,ista),indpk(:,ista)] = sort(pki,'ascend');
        tsep(:,ista) = diff(pkist(:,ista));
      end      
      ind = find(indpk(:,1)~=indpk(:,2) | indpk(:,1)~=indpk(:,3) | indpk(:,2)~=indpk(:,3));
%       disp(ind)
%       length(ind)
      
      %use function 'removesecondarysrc' to remove the smaller amplitude, secondary triplets that 
      %are too close in time to large triplets
      minsrcsep = 20;   % min allowed separation between deconvolved peaks
      [pkindepsave,indremove] = removesecondarysrc(pkindep,sigsta);
%       [pkindepsave,indremove] = removesecondarysrc(pkindep,sigsta,minsrcsep);
      
      %remove the secondary sources from the grouped result
      impindep(indremove, :) = [];
      impindepst = sortrows(impindep,1);
      impindepst(:,7:8) = impindepst(:,7:8)+repmat([off1i(k,2) off1i(k,3)],size(impindepst,1),1); %account for prealignment

      %% plot the scatter of sources in terms of offsets, accounting for prealignment offset
      span = max(range(off1iw(:,2))+2*loff_max, range(off1iw(:,3))+2*loff_max);
      xran = [round(mean(minmax(off1iw(:,2)'))-span/2)-1, round(mean(minmax(off1iw(:,2)'))+span/2)+1];
      yran = [round(mean(minmax(off1iw(:,3)'))-span/2)-1, round(mean(minmax(off1iw(:,3)'))+span/2)+1];
      cran = [0 lsig];
      %%%plot the scatter of offsets, accounting for prealignment offset, == true offset
      f1.fig = figure;
      f1.fig.Renderer = 'painters';
      ax1=gca;
      [ax1,torispl,mamp,xbndcvhl,ybndcvhl] = plt_decon_imp_scatter_ref(ax1,impindepst,xran,yran,cran,off1iw,loff_max,...
        sps,50,'mean','tori','comb');
      scatter(ax1,off1i(k,2),off1i(k,3),20,'ks','filled','MarkerEdgeColor','k');
      title(ax1,'Independent, grouped, no secondary sources');
      
%       %%%plot the scatter of offsets, shifted to the same origin, the best alignment
%       xran = [-loff_max+off1i(k,2)-1 loff_max+off1i(k,2)+1];
%       yran = [-loff_max+off1i(k,3)-1 loff_max+off1i(k,3)+1];
%       f1.fig = figure;
%       f1.fig.Renderer = 'painters';
%       ax1=gca;
%       ax1 = plt_decon_imp_scatter_ref_sft(ax1,impindepst,xran,yran,cran,off1i(k,:),off1iw,...
%         loff_max,irccran,sps,50,'mean','tori');
%       scatter(ax1,off1i(k,2),off1i(k,3),20,'ks','filled','MarkerEdgeColor','k');
%       title(ax1,'Independent, grouped, no secondary sources, shifted');

      %% plot the scatter of sources in terms of rela locations
      xran = [-5 5];
      yran = [-5 5];
      cran = [0 lsig/sps];
      f2.fig = figure;
      f2.fig.Renderer = 'painters';
      ax2=gca;
      [ax2] = plt_decon_imp_scatter_space_ref(ax2,impindepst,xran,yran,cran,off1iw,loff_max,...
        sps,50,ftrans,'mean','tori','comb');
      plot(ax2,xcut,ycut,'k-','linew',2);
      title(ax2,'Independent, grouped, no secondary sources');
      
      %% Compare deconvolved positive peaks with waveform peaks 
%       %%plot the positive peaks indicated by the grouped triplets, see if they indeed match the
%       %%peaks of the data at each station, could be helpful which triplets are minor
%       fpltremv = 1;   % plot or not the secondary sources that are removed as well
%       [f] = plt_deconpk_sigpk(sigsta,pkindep,indremove,fpltremv);
        
      %% separation in arrival time between deconvolved positive peaks
      %%%Plot the separation in time between these preserved positive peaks after removing the
      %%%secondary ones, to see if they can be too close to each other
      [f,tsep,pkist,indpk,indreverse] = plt_tsep_deconpk(pkindepsave,sps);
      
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
       
      %% what is the distance for consecutive soures, in terms of origin time?      
      %convert time offset to relative loc
      [imploc, indinput] = off2space002(impindepst(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
      [torisplst, indsort] = sortrows(torispl,1);
      dtorinn1 = diff(torisplst);
      implocst = imploc(indsort, :);
      impindepstst = impindepst(indsort, :);
      
      %%%For each LFE source, get its distance to all other LFEs
      [f] = plt_srcalldistCDF(impindepstst,implocst,torisplst,[],sps);
      [f] = plt_srcalldistCDF(impindepstst,implocst,torisplst,60*sps,sps);
      
      %between Nth and (N-1)th source
      distnn1 = sqrt(diffcustom(implocst(:,1),1,'forward').^2 + ...
        diffcustom(implocst(:,2),1,'forward').^2 );
      
      %between Nth and (N-2)th source
      dtorinn2 = diffcustom(torisplst,2,'forward');
      distnn2 = sqrt(diffcustom(implocst(:,1),2,'forward').^2 + ...
        diffcustom(implocst(:,2),2,'forward').^2 );
      
      %between Nth and (N-3)th source
      dtorinn3 = diffcustom(torisplst,3,'forward');
      distnn3 = sqrt(diffcustom(implocst(:,1),3,'forward').^2 + ...
        diffcustom(implocst(:,2),3,'forward').^2 );
      
      dist = distnn1;
      dift = dtorinn1;
      nsep = 1;
            
      %%%ECDF of sources in terms of offset 12, 13 and 23
      nsepvect = [1 2 3];
      [f] = plt_srcoffCDF(impindepstst,torisplst,nsepvect);
      
      %%%Absolute distance, in terms of origin time
      ttype = 'tori';
      mapview = 'offset';
      [f] = plt_srcdist(implocst,impindepstst,nsep,sps,dist,dift,torisplst,ttype,mapview);
      
      %%%ECDF of sources in terms of absolute distance
      nsepvect = 1:2:15;
      [f] = plt_srcdistCDF(implocst,torisplst,nsepvect);
            
      %%%Projected distance along specific directions, eg., propagation and its orthogonal, in terms of origin time
      [f] = plt_srcprojdist(implocst,nsep,sps,dist,dift,torisplst,ttype);
      
      %%%ECDF of sources in terms of projected distance in min rmse and max speed direction
      [f] = plt_srcprojdistCDF(implocst,nsepvect,sps,torisplst);  

      %% abnormal sources
%       ind = find(dift/sps<=0.1 & dist > median(dist(dift/sps<=1)));
%       nabn = length(ind);
%       abndist = dist(ind);
%       abndt = dift(ind);
%       abnlocl = implocst(ind+nsep, :);   %it should be a pair, earlier and later
%       abnloce = implocst(ind, :);
%       abntoril = torisplst(ind+nsep);
%       abntorie = torisplst(ind);
%       abnl = impindepst(indsort(ind+nsep),:);
%       abne = impindepst(indsort(ind),:);
%       abndtarvl = abnl(:,[1 3 5])-abne(:,[1 3 5]);
%       
%       mapview = 'offset';
%       [f] = plt_abnsrcdist_p1(abnloce,abnlocl,abne,abnl,nsep,sps,abndist,abndt,abntoril,dist,...
%         torisplst,mapview);
      
      %% abnormal sources, continued
%       [f] = plt_abnsrcdist_p2(abne,abnl,stas,abndtarvl,abndist);
      
      %% Are these abnormal sources consecutive peaks in arrival?
%       %%% ie. are there any other peaks in between?
%       for jj = 1:nabn
%         for ista = 1: nsta
%           ind = [];
%           if abnl(jj,(ista-1)*2+1) > abne(jj,(ista-1)*2+1)
%             ind = find(impindepst(:,(ista-1)*2+1)>abne(jj,(ista-1)*2+1) & ...
%               impindepst(:,(ista-1)*2+1)<abnl(jj,(ista-1)*2+1));
%           else
%             ind = find(impindepst(:,(ista-1)*2+1)>abnl(jj,(ista-1)*2+1) & ...
%               impindepst(:,(ista-1)*2+1)<abne(jj,(ista-1)*2+1));
%           end
%           disp(ind)
%         end
%       end

      %% what is the distance for consecutive sourcces, but in terms of arrival time?
      %note the 'tsep' obtained from the deconvolved positive peaks should be identical to that if
      %obtained from the deconvolved impulses themselves, which represent the arrival indices of the
      %zero-crossing
      dtarvlnn1 = [];
      distnn1 = [];
      for ista = 1: nsta
        [impindepstst, indsort] = sortrows(impindepst, (ista-1)*2+1);
        implocst = imploc(indsort, :);
        tarvlsplst = impindepstst(:,(ista-1)*2+1);

        nsep = 1;
        dtarvlnn1(:,ista) = diffcustom(tarvlsplst, nsep,'forward');
        distnn1(:,ista) = sqrt(diffcustom(implocst(:,1),nsep,'forward').^2 + ...
          diffcustom(implocst(:,2),nsep,'forward').^2 );
        dift = dtarvlnn1(:,ista);
        dist = distnn1(:,ista);
        
        %%%Absolute distance, in terms of arrival time 
        ttype = 'tarvl';
        mapview = 'offset';
        [f] = plt_srcdist(implocst,impindepstst,nsep,sps,dist,dift,tarvlsplst,ttype,mapview);
        
        %%%Projected distance along specific directions, eg., propagation and its orthogonal, in terms of arrival time
        [f] = plt_srcprojdist(implocst,nsep,sps,dist,dift,tarvlsplst,ttype);

      end
      
      
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
      mapview = 'offset';
      [f] = plt_srcoffset23(imploc,impindepst,torispl,indsort,mapview);
            
      %% what's the corresponding rcc at source arrival, associated with the offset
      [f] = plt_srcrccwithoffset(impindepst,rccpaircat);
     
      %% linear regression of offset change to know how fast the centroid of sources are moving
      [f,fitobj12,fitobj13,dyfit12,dyfit13] = plt_deconlfit(impindepst,torispl);
      
      %% plot amp ratio 12 and 13, and 23
      mapview = 'offset';
      [f] = plt_srcampratio(imploc,impindepst,greenf,sps,mapview);
      
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
      predgrp = zeros(size(sigsta));  % predefine the prediction array
      impfull = zeros(size(sigsta));  % a full array of impulses containing zeros
      predgrpl = zeros(lsig+overshoot*2, nsta);   % longer including an overshoot on both sides for alignment
      impfulll = zeros(lsig+overshoot*2, nsta); 
      for ista = 1: nsta
        impfull(impindepst(:,(ista-1)*2+1), ista) = impindepst(:,ista*2); % the non-zero index has the amp 
        predtmp = conv(greenf(:,ista), impfull(:,ista), 'full');
        twlet = zcrosses(ista)*dt;
        itwlet = round(twlet/dt);
        predgrp(:,ista) = predtmp(itwlet:lsig+itwlet-1);  % cut accordingly
        
        %obtain a longer predication and residual for alignment in each subwin
        impfulll(impindepst(:,(ista-1)*2+1)+overshoot, ista) = impindepst(:,ista*2);
        predtmpl = conv(greenf(:,ista), impfulll(:,ista), 'full');
        predgrpl(:,ista) = predtmpl(itwlet:lsig+overshoot*2+itwlet-1);  % cut accordingly
      end
      resgrp = sigsta - predgrp;  % the residual, different between signal and prediction
      resgrpl = detrend(optdat(:,2:end)) - predgrpl;  % longer residual

%       figure;
%       subplot(211)
%       plot(predtmp,'b');hold on;
%       plot(predtmpl,'r');
%       subplot(212)
%       plot(predgrp(:,3),'b');hold on;
%       plot(predgrpl(:,3),'r');
      
      %%%obtain an concatenated RCC for the residual as well, you can run this part after removing
      %%%the secondary sources, but to justify the use of coda in templates if sources are widely-
      %%%distributed, you can run it right after grouping. You should see the rcc of the residual is
      %%%pretty low (eg., averaged to 0) compared to signal
      %%%after grouping
      irccrcat = [];   % concatenated indices of RCC for residual
      rccrcat = [];  % concatenated RCC for residual
      reslwin = windows-windows(1,1)+1+overshoot; % segmentation of windows for residual
      for iwin = 1: nwin
        ist = reslwin(iwin,1);
        ied = reslwin(iwin,2);
        resdat = [];
        resdat(:, 1) = resgrpl(ist: ied, 1); % sta 1
        resdat(:, 2) = resgrpl(ist-off1iw(iwin,2): ied-off1iw(iwin,2), 2); % sta 2
        resdat(:, 3) = STAopt(ist-off1iw(iwin,3): ied-off1iw(iwin,3), 3); % sta 3
        %compute running CC between 3 stations
        [irccr,rccr12] = RunningCC(resdat(:,1), resdat(:,2), mwlen);
        [~,rccr13] = RunningCC(resdat(:,1), resdat(:,3), mwlen);
        [~,rccr23] = RunningCC(resdat(:,2), resdat(:,3), mwlen);
        rccr = (rccr12+rccr13+rccr23)/3;  % use averge
        irccr = irccr + reslwin(iwin,1) - reslwin(1,1);   %convert to global index
        irccrcat = [irccrcat; irccr];
        rccrcat = [rccrcat; rccr];
      end

      %%% what about the orthogonal component? 
      nfft = lsig;
      coef = [];
      lag = [];
      for ista = 1:nsta
        [coef(:,ista), lag(:,ista)] = xcorr(sigsta(:,ista), greenf(:,ista), nfft, 'none'); % unnormalized master raw CC
      end
      maxlag = 2*sps;
      %%%for stas 1, 2, 3, between opt and ort
      coefort = [];
      lagort = [];
      for ista = 1:nsta
        [coefort(:,ista), lagort(:,ista)] = xcorr(sigstaort(:,ista), greenfort(:,ista), nfft, 'none');
      end
      for ista = 1:nsta
        %CC opt with ort
        [cccoef, lagcoef] = xcorr(coef(:,ista), coefort(:,ista), maxlag, 'coeff');
        [ccboo(k,ista), mind] = max(cccoef);
        lagboo(k,ista) = lagcoef(mind);
      end

      predgrport = zeros(size(sigstaort));  % predefine the prediction array
      impfull = zeros(size(sigsta));  % a full array of impulses containing zeros
      for ista = 1: nsta
        impfull(impindepst(:,(ista-1)*2+1), ista) = impindepst(:,ista*2); % the non-zero index has the amp 
        predtmp = conv(greenfort(:,ista), impfull(:,ista), 'full');
        twlet = zcrosses(ista)*dt;
        itwlet = round(twlet/dt);
        predgrport(:,ista) = predtmp(itwlet:lsig+itwlet-1);  % cut accordingly
        
      end
      resgrport = sigstaort - predgrport;  % the residual, different between signal and prediction
      
      figure
      subplot(221)
      hold on
      box on; grid on
      plot(sigsta(:,1),'r');
      plot(sigsta(:,2),'b');
      plot(sigsta(:,3),'k');
      ym = max(abs(sigsta(:)));
      yran=1.2*[-ym ym];
      ylim(yran);
      text(0.95,0.9,sprintf('%.2f; %.2f; %.2f',norm(sigsta(:,1)),norm(sigsta(:,2)),...
        norm(sigsta(:,3))),'Units','normalized','HorizontalAlignment','right');
      subplot(222)
      hold on
      box on; grid on
      plot(resgrp(:,1),'r');
      plot(resgrp(:,2),'b');
      plot(resgrp(:,3),'k');
      ylim(yran);
      text(0.95,0.9,sprintf('%.2f; %.2f; %.2f',norm(resgrp(:,1)),norm(resgrp(:,2)),...
        norm(resgrp(:,3))),'Units','normalized','HorizontalAlignment','right');
      text(0.95,0.8,sprintf('%.1f%%; %.1f%%; %.1f%%',...
        (norm(sigsta(:,1))-norm(resgrp(:,1)))/norm(sigsta(:,1))*100,...
        (norm(sigsta(:,2))-norm(resgrp(:,2)))/norm(sigsta(:,2))*100,...
        (norm(sigsta(:,3))-norm(resgrp(:,3)))/norm(sigsta(:,3))*100),...
        'Units','normalized','HorizontalAlignment','right');
      text(0.95,0.1,sprintf('%.2f; %.2f; %.2f',ccboo(:,1),ccboo(:,2),...
        ccboo(:,3)),'Units','normalized','HorizontalAlignment','right');
      subplot(223)
      hold on
      box on; grid on
      plot(sigstaort(:,1),'r');
      plot(sigstaort(:,2),'b');
      plot(sigstaort(:,3),'k');
      ym = max(abs(sigstaort(:)));
      yran=1.2*[-ym ym];
      ylim(yran);
      text(0.95,0.9,sprintf('%.2f; %.2f; %.2f',norm(sigstaort(:,1)),norm(sigstaort(:,2)),...
        norm(sigstaort(:,3))),'Units','normalized','HorizontalAlignment','right');
      subplot(224)
      hold on
      box on; grid on
      plot(resgrport(:,1),'r');
      plot(resgrport(:,2),'b');
      plot(resgrport(:,3),'k');
      ylim(yran);
      text(0.95,0.9,sprintf('%.2f; %.2f; %.2f',norm(resgrport(:,1)),norm(resgrport(:,2)),...
        norm(resgrport(:,3))),'Units','normalized','HorizontalAlignment','right');
      text(0.95,0.8,sprintf('%.1f%%; %.1f%%; %.1f%%',...
        (norm(sigstaort(:,1))-norm(resgrport(:,1)))/norm(sigstaort(:,1))*100,...
        (norm(sigstaort(:,2))-norm(resgrport(:,2)))/norm(sigstaort(:,2))*100,...
        (norm(sigstaort(:,3))-norm(resgrport(:,3)))/norm(sigstaort(:,3))*100),...
        'Units','normalized','HorizontalAlignment','right');

      %%
      %%%plot agu 2022 abstract
%       [f] = plt_agu2022abstract(sigsta,resgrp,impindepst,implocst,torisplst,nsep,sps,tmaxi,...
%         tmaxo,tbosti,tbosto,ircccat,rcccat,irccrcat,rccrcat,off1iw,off1i,loff_max,tstbuf,tedbuf,...
%         dy,mo,yr,k,ftrans,xcut,ycut);
%       print(f.fig,'-dpdf',fullfile(rstpath, '/FIGS','agu2022abstract.pdf'));
      
      xzoom = [225 245];
%       xzoom = [20 35];
      [f] = plt_agu2022abstractv2(greenf,sigsta,impindepst,sps,xzoom,off1iw,loff_max,tstbuf,...
        dy,mo,yr,ftrans);
      print(f.fig,'-dpdf',fullfile(rstpath, '/FIGS','agu2022abstractv2.pdf'));
      
      xzoom = [231 235];
      ind = find(impindepst(:,1)/sps >= xzoom(1) & impindepst(:,1)/sps <= xzoom(2));
      impindepstaa = impindepst(ind,:);
      torisplaa = torispl(ind);
      %convert time offset to relative loc
      [imploc, indinput] = off2space002(impindepstaa(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
      [torisplst, indsort] = sortrows(torisplaa,1);
      dtorinn1 = diff(torisplst);
      implocst = imploc(indsort, :);
      impindepstst = impindepstaa(indsort, :);

      %between Nth and (N-1)th source
      distnn1 = sqrt(diffcustom(implocst(:,1),1,'forward').^2 + ...
        diffcustom(implocst(:,2),1,'forward').^2 );
      
      dist = distnn1;
      dift = dtorinn1;
      nsep = 1;
      ttype = 'tori';
      [f] = plt_srcprojdist(implocst,nsep,sps,dist,dift,torisplst,ttype);


    end
    
  end
  
end








