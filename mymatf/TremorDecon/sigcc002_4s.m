% sigcc002_4s.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To some extent quantify if the 4th or more station could be useful into
% the deconvolution together or just use it as a check, another option, and
% possibly the simplest option is to directly cross-correlate the seismograms
% between the 4th station and trio stations. If, the deconvolution works
% better than counting the peaks directly by using the templates that has 
% coda beside the main dipole, then this CC between sigs should be lower in 
% general than the CC between sig-wlet CCs (done in 'sigwletcc002_4s.m').
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/08/04
% Last modified date:   2022/08/04
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
  'SILB '
  'KLNB '
  'LZB  '
  'TWKB '
  'MGCB '
  ];     % determine the trio and order, here the 1st sta is PGC
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
tranbst = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(round(ttol)),'s.pgc002.',cutout(1:4)));
tlen = tranbst(:,3)-tranbst(:,2);
nbst = size(tranbst,1);

ttol = 1e-3*86400;
tranmig = load(strcat(rstpath, '/MAPS/migran',num2str(round(ttol)),'s.pgc002'),'w+');
tranmig(:,2:3) = tranmig(:,2:3)/3600;
nmig = size(tranmig,1);

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
normflag = 0;

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
for ista = 4: nsta
  tempxc(:,ista)=xcorr(STAtmp(:,1),STAtmp(:,ista),mshiftadd,'coeff'); %shift STAtmp(2,:) to right for positive values
end
[~,imax]=max(tempxc,[],1);
imax=imax-(mshiftadd+1); %This would produce a slightly different shift, if filtered seisms were used.
imax(2)-imax(3)+imax(1);   %enclosed if it equals to 0
for ista=2:nsta
    STAtmp(mshiftadd+1:end-(mshiftadd+1),ista)=STAtmp(mshiftadd+1-imax(ista):end-(mshiftadd+1)-imax(ista),ista);
    STAtmport(mshiftadd+1:end-(mshiftadd+1),ista)=STAtmport(mshiftadd+1-imax(ista):end-(mshiftadd+1)-imax(ista),ista);
end
%normalization
if normflag 
  for ista=1:nsta
      STAtmp(:,ista)=STAtmp(:,ista)/spread(ista); % now templates are 'aligned' indeoendently by x-corr wrt. sta 1
      STAtmport(:,ista)=STAtmport(:,ista)/spread(ista);
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
[off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(tmpwletf(:,1:3)',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
offwlet1i(1) = 0;
offwlet1i(2) = round(off12con);
offwlet1i(3) = round(off13con);

for ista = 4: nsta
  [coef,lag] = xcorr(tmpwletf(:,1), tmpwletf(:,ista), mshiftadd, 'coeff');
  [mcoef, idx] = max(coef);   % max of master raw cc
  offwlet1i(ista) = lag(idx);   % offset in samples  
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
for ista = 1: nsta
  %cut according to the zero-crossing and the time shift from the constrained CC
  green(:,ista) = tmpwlet(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenf(:,ista) = tmpwletf(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  %detrend again for caution
  green(:,ista)=detrend(green(:,ista));
  greenf(:,ista)=detrend(greenf(:,ista));
  if normflag
    %normalize by max amp
    green(:,ista)=green(:,ista)/max(abs(green(:,ista)));    % normalize
    greenf(:,ista)=greenf(:,ista)/max(abs(green(:,ista)));    % normalize
  end
  
  %same process for orthogonal
  greenort(:,ista) = tmpwletort(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenfort(:,ista) = tmpwletfort(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenort(:,ista)=detrend(greenort(:,ista));
  greenfort(:,ista)=detrend(greenfort(:,ista));
  if normflag
    greenort(:,ista)=greenort(:,ista)/max(abs(green(:,ista)));    % normalize
    greenfort(:,ista)=greenfort(:,ista)/max(abs(green(:,ista)));    % normalize
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
[off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(greenf(:,1:3)',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
if ~(off12con==0 && off13con==0)
  disp('Filtered templates are NOT best aligned \n');
end
for ista = 4: nsta
  [coef,lag] = xcorr(greenf(:,1), greenf(:,ista), mshiftadd, 'coeff');
  [mcoef, idx] = max(coef);   % max of master raw cc
  if lag(idx)~=0   % offset in samples
    fprintf('Filtered templates are NOT best aligned at %s \n',stas(ista,:));
  end
end

amprat(1,:) = minmax(greenf(:,1)')./minmax(greenf(:,2)');	% amp ratio between max at sta 3 and 2 or min
amprat(2,:) = minmax(greenf(:,1)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min  
amprat(3,:) = minmax(greenf(:,2)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min  
spread = range(greenf);   % range of the amp of template

%%%plot the unfiltered and filtered templates
% plt_templates(green,greenf,stas,greenort,greenfort,lowlet,hiwlet,sps);

%just the filtered templates
% plt_templates_bp(greenf,stas,lowlet,hiwlet,sps);


%%
% dates and ets
%%% NOTE: 'dates' is the dates that the tremor is active at the region of interest, so that the
%%% waveform at the stations on these dates you are looking at mainly comes from the sources at the
%%% region of interest. We want to see if there is a noticable change in spectra during the burst
%%% windows on these dates
dates = unique(tranbst(:,1));
years = unique(floor(dates/1000));
nets = length(years);

%filtering passband for reading data, confirmed by 'spectrabursts002_4s.m'
hisig=6.3; % this will give a similar spectral shape between template and signal
losig=1.8;

k = 0;  %deconvolution burst win count
n = 0;  %auto-determined migration win count

sps = 160;

%empirically determined indices of bursts fall into local day times (noisier)
%or night times (quieter)
inbst = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,...
  34,35,56,57,58,59,60,61,62,70,71,72,73,74,75,76,77,78,79,80,81,82,110,111,112,113,114,115,116,117,...
  118,119,120,142,143,144,145,146,147,148,149,150,151,152,153,154,172,173,174,175];
idbst = [36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,63,64,65,66,67,68,69,83,84,85,...
  86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,121,122,123,124,...
  125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,155,156,157,158,159,160,161,...
  162,163,164,165,166,167,168,169,170,171,176,177,178,179,180,181,182,183,184,185,186,187,188,189,...
  190,191,192,193,194,195];

%indices of bursts whose wlet-sig cc between the opt and ort components at sta 1/2/3 are very
%similar, above 75 percentile
ihioo123 = [20,23,59,60,80,113,116,120,134,189,194];
ihioo123n = intersect(inbst,ihioo123);

%indices of bursts whose sig cc between sta 1/2/3 are very similar, above 90 percentile
ihicc123 = [1,3,6,7,8,24,56,71,75,77,81,83,93,102,114,116,132,145,149,185];
ihicc123n = intersect(inbst,ihicc123);

idxburst = 1:size(tranbst,1);

nibst = size(idxburst,1);
[~,~,~,idate,ibst] = indofburst(tranbst,idxburst);

off1i = zeros(nibst,nsta);

ccmij = zeros(3, nmig, nsta-3);   %coef of CC at 4th sta with that at sta 1/2/3 for migs
lagmij = zeros(3, nmig, nsta-3);  %lag of CC at 4th sta with that at sta 1/2/3 for migs
ccm123 = zeros(nmig, 3);  %coef of CC within sta 1,2,3 for migs
lagm123 = zeros(nmig, 3); %lag of CC within sta 1,2,3 for migs
ccm45 = zeros(nmig, 6);  %coef of CC within sta 4,5,6,7 for migs
lagm45 = zeros(nmig, 6); %lag of CC within sta 4,5,6,7 for migs

ccbij = zeros(3, nibst, nsta-3); %coef of CC at 4th sta with that at sta 1/2/3 for bursts
lagbij = zeros(3, nibst, nsta-3);  %lag of CC at 4th sta with that at sta 1/2/3 for bursts
ccb123 = zeros(nibst, 3);  %coef of CC within sta 1,2,3 for bursts
lagb123 = zeros(nibst, 3); %lag of CC within sta 1,2,3 for bursts
ccb45 = zeros(nibst, 6);  %coef of CC within sta 4,5,6,7 for bursts
lagb45 = zeros(nibst, 6); %lag of CC within sta 4,5,6,7 for bursts
ccboo = zeros(nibst, 3); %coef of CC opt comp. with that of ort within sta 1,2,3 for bursts
lagboo = zeros(nibst, 3);  %lag of CC opt comp. with that of ort within sta 1,2,3 for bursts

%%%Flag to indicate if it is necessary to recalculate everything
flagrecalc = 0;
% flagrecalc = 1;

if flagrecalc
  for iii = 1: length(idxburst)
    iii
    % for i = 1: length(dates)  % dates in each ets
    i = idate(iii);
    date = dates(i);
    year = floor(date/1000);
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
    
    %bursts and 4-s detections of the same day
    rangetemp = tranbst(tranbst(:,1)==date, :);
    hfdayi = hfbnd(hfbnd(:,daycol)==date, :);  % inside bound of the day
    hfdayo = hfout(hfout(:,daycol)==date, :);  % outside bound of the day
    k = k+size(rangetemp,1);
    
    %migrations of the same day
    mig = tranmig(tranmig(:,1)==date, :);
    n = n+size(mig,1);
    
    tmaxi = hfdayi(:, seccol); % starting time of max power rate of half sec inside the ellipse
    tmaxo = hfdayo(:, seccol); % starting time of max power rate of half sec outside the ellipse
    tcnti = hfdayi(:, 15);  % the center of detecting win is the 15th col
    tcnto = hfdayo(:, 15);  % the center of detecting win is the 15th col
    
    %read horizontal optimal and orthogonal components
    JDAY = num2zeropadstr(jday,3);
    MO=day2month(jday,year);     % EXTERNAL function, day2month, get the month of one particular date
    direc=[datapath, '/arch', yr,'/',MO,'/'];     % directory name
    prename=[direc,yr,'.',JDAY,'.00.00.00.0000.CN'];    %  path plus prefix of data file,
    %     disp(prename);
    [STAopt,STAort,~,fileflag] = rd_daily_bpdata(year,jday,prename,stas,PERMSTA,POLSTA,...
      PERMROTS,POLROTS,sps,losig,hisig,npo,npa,[],[],[],[]);
    
    %   if fileflag == 0    % means there are missing files
    %     fprintf('Day %s / %s will be omitted because of missing files. \n', yr, JDAY);
    %     continue    % continue to the next day
    %   end
    
    %   %read vertical components too, as we want to have a sense of the SNR at the orthogonal and vertical
    %   %components
    %   [STAvert,~,~,fileflag] = rd_daily_bpdata(year,jday,prename,stas,PERMSTA,POLSTA,...
    %     PERMROTS,POLROTS,sps,losig,hisig,npo,npa,[],[],[],[],'Z');
    
    %%%CC between stations for decon windows
    %   for j = 1: size(rangetemp,1)
    j = ibst(iii);
    tst = rangetemp(j,2); % start and end time of bursts
    ted = rangetemp(j,3);
    
    icount = iii;
    %     icount = j+k-size(rangetemp,1);
    
    %how many 4-s detections fall into the burst range
    indtmaxi = find(tmaxi>=tst-0.1 & tmaxi<=ted+0.1);
    
    %%%%Use the start and end of the 4-s detecting window
    %     tstbuf = min(tcnti(indtmaxi)-2)-4;  %with some buffer
    %     tedbuf = max(tcnti(indtmaxi)+2)+4;
    tstbuf = min(tcnti(indtmaxi)-2);  %no buffer
    tedbuf = max(tcnti(indtmaxi)+2);
    
    optcc = STAopt(max(floor((tstbuf+1)*sps+1),1): min(floor((tedbuf-1)*sps),86400*sps), 2:nsta+1);
    msftadd = 10/40*sps;
    ccmid = ceil(size(optcc,1)/2);
    ccwlen = round(size(optcc,1)-2*(msftadd+1));
    loffmax = 4*sps/40;
    ccmin = 0.01;  % depending on the length of trace, cc could be very low
    iup = 1;    % times of upsampling
    [off12con,off13con,ccali(icount),iloopoff,loopoff] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
      ccwlen,msftadd,loffmax,ccmin,iup);
    % if a better alignment cannot be achieved, use 0,0
    if off12con == msftadd+1 && off13con == msftadd+1
      off12con = 0;
      off13con = 0;
      fprintf('Tremor burst %d cannot be properly aligned, double-check needed \n',k);
    end
    off1i(icount,1) = 0;
    off1i(icount,2) = round(off12con);
    off1i(icount,3) = round(off13con);
    
    for ista = 4: nsta
      [coef,lag] = xcorr(optcc(:,1), optcc(:,ista), msftadd, 'coeff');
      [mcoef, idx] = max(coef);   % max of master raw cc
      off1i(icount,ista) = lag(idx);   % offset in samples
    end
    
    %Align records
    optdat = [];  % win segment of interest
    ortdat = [];
    optdat(:, 1:2) = STAopt(max(floor(tstbuf*sps+1),1): min(floor(tedbuf*sps),86400*sps), 1:2); % sta 1
    ortdat(:, 1:2) = STAort(max(floor(tstbuf*sps+1),1): min(floor(tedbuf*sps),86400*sps), 1:2);
    for ista = 2: nsta
      optdat(:, ista+1) = STAopt(max(floor(tstbuf*sps+1)-off1i(icount,ista),1): ...
        min(floor(tedbuf*sps)-off1i(icount,ista),86400*sps), ista+1); % sta 2
      ortdat(:, ista+1) = STAort(max(floor(tstbuf*sps+1)-off1i(icount,ista),1): ...
        min(floor(tedbuf*sps)-off1i(icount,ista),86400*sps), ista+1);
    end
    
    sigsta = zeros(size(optdat,1), nsta);
    for ista = 1:nsta
      tmp = optdat(:,ista+1); %best aligned, filtered
      tmp = detrend(tmp);
      sigsta(:,ista) = tmp;
    end
    sigstaort = zeros(size(ortdat,1), nsta);
    for ista = 1:nsta
      tmp = ortdat(:,ista+1); %best aligned, filtered
      tmp = detrend(tmp);
      sigstaort(:,ista) = tmp;
    end
    
    %     figure
    %     hold on;
    % %     plot(optdat(:,1),sigsta(:,1),'r');
    %     plot(optdat(:,1),sigsta(:,2),'b');
    %     plot(optdat(:,1),sigsta(:,3),'k');
    %     xlabel('Time (s)');
    %     ylabel('Amplitude');
    
    maxlag = 2*sps;
    
    %%%among stas 4, 5, 6, 7, between opt and opt, in order 45 46 47 56 57 67
    tmp = 0;
    for mm = 4: nsta-1
      for nn = mm+1: nsta
        tmp = tmp+1;
        [coef, lag] = xcorr(sigsta(:,mm), sigsta(:,nn), maxlag, 'coeff');
        [ccb45(icount,tmp), mind] = max(coef);
        lagb45(icount,tmp) = lag(mind);
      end
    end
    
    %%%for 4th and more stations, relative to sta 1/2/3, between opt and opt
    for mm = 4: nsta
      for nn = 1: 3
        [coef, lag] = xcorr(sigsta(:,nn), sigsta(:,mm), maxlag, 'coeff');
        [ccbij(nn,icount,mm-3), mind] = max(coef);
        lagbij(nn,icount,mm-3) = lag(mind);
      end
    end
    
    
    %%%among stas 1, 2, 3, between opt and opt, in order 12 13 23
    tmp = 0;
    for mm = 1: 3-1
      for nn = mm+1: 3
        tmp = tmp+1;
        [coef, lag] = xcorr(sigsta(:,mm), sigsta(:,nn), maxlag, 'coeff');
        [ccb123(icount,tmp), mind] = max(coef);
        lagb123(icount,tmp) = lag(mind);
      end
    end    
    
    %   end
    % end
    
  end
  
  %%% save some variables
  savefile = 'rst_sigcc_dtr.mat';
  save(strcat(rstpath, '/MAPS/',savefile), 'off1i','ccbij','lagbij','ccb123','lagb123','ccb45',...
    'lagb45');
  
else
  maxlag = 2*sps;
  savefile = 'rst_sigcc_dtr.mat';
  load(strcat(rstpath, '/MAPS/',savefile));
end


%% target some high-correlation bursts
% ind = find(ccb123(:,1)>=prctile(ccb123(:,1),75) & ccb123(:,2)>=prctile(ccb123(:,2),75) & ...
%   ccb123(:,3)>=prctile(ccb123(:,3),75));
% ind = find(ccb123(:,3)>=prctile(ccb123(:,3),75));
% ind = [20,23,59,60,80,113,116,120,134,189,194];

%% burst windows for stas 4/5/6/7 vs. 1/2/3
%%%scatter of lag and CC 
widin = 12;
htin = 9;
nrow = 3;
ncol = nsta-3;
pltxran = [0.06 0.96]; pltyran = [0.06 0.96]; % optimal axis location
pltxsep = 0.03; pltysep = 0.03;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

for ii = 1:nrow
  for jj = 1:ncol
    isub = (ii-1)*ncol+jj;
    ax = f.ax(isub);
    hold(ax,'on');
    ax.Box = 'on';
    grid(ax,'on');
    scatter(ax,lagbij(ii,:,jj)/sps,ccbij(ii,:,jj),16,...
      'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.15);
    ylim(ax,[0.0 0.7]);
    xlim(ax,[-maxlag,maxlag]/sps);
    text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(jj+3,:)),strtrim(stas(ii,:))),'Units',...
      'normalized');
    text(ax,0.98,0.05,sprintf('%.2f, %.2f',median(lagbij(ii,:,jj)/sps),median(ccbij(ii,:,jj))),...
      'Units','normalized','HorizontalAlignment','right');
    if jj ~= 1
      nolabels(ax,2);
    end
    if ii ~= nrow
      nolabels(ax,1);
    end
    if ii == nrow && jj==1
      xlabel(ax,'Lag (s) of max CC');
      ylabel(ax,'Max CC');
    end
    longticks(ax,2);
  end
end
supertit(f.ax(1:ncol),'cc of sig; bursts; 4th stas vs. trio stas');

% %%%CDF of CC 
% figure
% nrow = 3;
% ncol = nsta-3;
% color=jet(ncol);
% p = [];
% for ii = 1:nrow
%   ax=subplot(1,nrow,ii);
%   hold on
%   for jj = 1:ncol
%     [cdfval,x] = ecdf(ccbij(ii,:,jj)); %between Nth and (N-1)th source
%     if ii == 1
%       p(jj)=plot(ax,x,cdfval,'linew',1,'color', color(jj,:));
%     else
%       plot(ax,x,cdfval,'linew',1,'color', color(jj,:));
%     end
%   end
%   xlim(ax,[0.0 0.7]);
%   grid on; box on;
%   if ii == 1
%     ylabel(ax,'Empirical CDF');
%     xlabel(ax,'Max CC');
%   else
%     nolabels(ax,2);
%   end
%   text(ax,0.98,0.6,sprintf('wrt. %s',strtrim(stas(ii,:))),'Units',...
%       'normalized','HorizontalAlignment','right');
%   if ii == nrow
%     legend(p,stas(ncol:end, :),'Location','southeast');
%   end
% end
% 
% keyboard

%% burst windows for stas 4/5/6/7 
%%%scatter of lag and CC 
widin = 9;
htin = 6;
nrow = 2;
ncol = 3;
pltxran = [0.06 0.96]; pltyran = [0.08 0.96]; % optimal axis location
pltxsep = 0.03; pltysep = 0.03;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

isub = 0;
for ii = 4: nsta-1
  for jj = ii+1: nsta
    isub = isub+1;
    ax = f.ax(isub);
    hold(ax,'on');
    ax.Box = 'on';
    grid(ax,'on');
    scatter(ax,lagb45(:,isub)/sps,ccb45(:,isub),16,...
    'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.15);
    ylim(ax,[0.0 0.7]);
    xlim(ax,[-maxlag,maxlag]/sps);
    text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(ii,:)),strtrim(stas(jj,:))),'Units',...
      'normalized');
    text(ax,0.98,0.05,sprintf('%.2f, %.2f',median(lagb45(:,isub)/sps),median(ccb45(:,isub))),...
      'Units','normalized','HorizontalAlignment','right');
    if rem(isub,ncol) ~= 1
      nolabels(ax,2);
    end
    if isub <= (nrow-1)*ncol
      nolabels(ax,1);
    end
    if isub == (nrow-1)*ncol+1
      xlabel(ax,'Lag (s) of max CC');
      ylabel(ax,'Max CC');
    end
    longticks(ax,2);
  end
end
supertit(f.ax(1:ncol),'cc of sig; bursts; among 4th stas');


%% burst windows for stas 1/2/3
%%%scatter of lag and CC 
widin = 9;
htin = 3.5;
nrow = 1;
ncol = 3;
pltxran = [0.06 0.96]; pltyran = [0.12 0.94]; % optimal axis location
pltxsep = 0.03; pltysep = 0.03;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

for isub = 1:3
  ax = f.ax(isub);
  hold(ax,'on');
  ax.Box = 'on';
  grid(ax,'on');
  scatter(ax,lagb123(:,isub)/sps,ccb123(:,isub),16,...
    'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.15);
  ylim(ax,[0.1 0.8]);
  xlim(ax,[-maxlag,maxlag]/sps);
  if isub ==1
    text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(1,:)),strtrim(stas(2,:))),'Units',...
      'normalized');
  elseif isub ==2
    text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(1,:)),strtrim(stas(3,:))),'Units',...
      'normalized');
  else
    text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(2,:)),strtrim(stas(3,:))),'Units',...
      'normalized');
  end
  text(ax,0.98,0.05,sprintf('%.2f, %.2f',median(lagb123(:,isub)/sps),median(ccb123(:,isub))),...
    'Units','normalized','HorizontalAlignment','right');
  if isub~=1
    nolabels(ax,2);
  end
  if isub==1
    xlabel(ax,'Lag (s) of max CC');
    ylabel(ax,'Max CC');
  end
  longticks(ax,2);
end
supertit(f.ax(1:ncol),'cc of sig; bursts; among trio stas');

% %%%CDF of CC 
% figure
% ncol = 3;
% color=jet(ncol);
% p = [];
% hold on
% ax = gca;
% for jj = 1:ncol
%   [cdfval,x] = ecdf(ccb123(:,jj)); %between Nth and (N-1)th source
%   p(jj)=plot(ax,x,cdfval,'linew',1,'color', color(jj,:));
% end
% xlim(ax,[0.1 0.8]);
% grid on; box on;
% ylabel(ax,'Empirical CDF');
% xlabel(ax,'Max CC');
% 
% legend(ax,p,'PGC-SSIB','PGC-SILB','SSIB-SILB','Location','southeast');

% %% burst windows for stas 4/5/6/7 vs. 1/2/3, minus reference
% widin = 12;
% htin = 9;
% nrow = 3;
% ncol = nsta-3;
% pltxran = [0.06 0.96]; pltyran = [0.06 0.96]; % optimal axis location
% pltxsep = 0.03; pltysep = 0.03;
% f = initfig(widin,htin,nrow,ncol);
% optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
% 
% for ii = 1:nrow
%   for jj = 1:ncol
%     isub = (ii-1)*ncol+jj;
%     ax = f.ax(isub);
%     hold(ax,'on');
%     ax.Box = 'on';
%     grid(ax,'on');
%     if ii == 1
%       ref = mean(ccb123(:,[1 2]),2)';
%       text(ax,0.98,0.95,'ref=mean(12+13)','Units','normalized','HorizontalAlignment','right');
%     elseif ii == 2
%       ref = mean(ccb123(:,[1 3]),2)';
%       text(ax,0.98,0.95,'ref=mean(12+23)','Units','normalized','HorizontalAlignment','right');
%     else
%       ref = mean(ccb123(:,[2 3]),2)';
%       text(ax,0.98,0.95,'ref=mean(13+23)','Units','normalized','HorizontalAlignment','right');
%     end
%     scatter(ax,lagbij(ii,:,jj)/sps,ccbij(ii,:,jj)-ref,8,1:size(ccbij,2),'filled');
%     colormap(ax,'jet');
%     ylim(ax,[-0.5 0.2]);
%     xlim(ax,[-maxlag,maxlag]/sps);
%     text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(jj+3,:)),strtrim(stas(ii,:))),'Units',...
%       'normalized');
%     text(ax,0.98,0.05,sprintf('%.2f, %.2f',median(lagbij(ii,:,jj)/sps),median(ccbij(ii,:,jj)-ref)),...
%       'Units','normalized','HorizontalAlignment','right');
%     if jj ~= 1
%       nolabels(ax,2);
%     end
%     if ii ~= nrow
%       nolabels(ax,1);
%     end
%     if ii == nrow && jj==1
%       xlabel(ax,'Lag (s) of max CC');
%       ylabel(ax,'Max CC - ref CC');
%       colorbar(ax,'Location','west');
%     end
%   end
% end
% supertit(f.ax(1:ncol),'Bursts (-reference)');











