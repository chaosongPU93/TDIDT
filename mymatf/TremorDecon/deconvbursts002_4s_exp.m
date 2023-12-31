% deconvbursts002_4s_exp.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to 'deconvbursts002_4s', but now 'deconvbursts002_4s' is 
% used to carry out all USEFUL operations that are common to all burst
% windows. And this script serves as experiments could be useful to some
% specific examples. 
% 
% Check 'deconvbursts002_4s' for more instructions and tips for carrying out
% deconvolution to all windows in general.
%
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/07/08
% Last modified date:   2022/07/08
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

%%%load the tremor catalog of John Armbruster
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

off1i = zeros(size(trange,1),3);  % the best alignment between sta2, sta3 wrt sta1
ccali = zeros(size(trange,1),1);  % CC value using the best alignment
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
%     disp(prename);
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
    for j = 14: size(rangetemp,1)  
%       close all
      k = k+1;  
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
            
      %have some overshoot, so that the resulted rcc would have the same length as the signal
      overshoot = mwlen/2;
%       overshoot = 0;

      %FLAG to simulate the behavior of noise                    
      noiseflag = 1;
      
      randsftlm = 5*(rand(200,1)-0.5)+7.5; % uniform between 5 and 10
      randsftrm = 5*(rand(200,1)-0.5)+7.5; % uniform between 5 and 10
%       figure; histogram(randsftlm);
%       hold on; histogram(randsftrm);
      
      for iii = 1: length(randsftlm)

      if noiseflag
%         randsftl = 5;   % randomly shift records by a few sec so that CC between them should behave like noise
%         randsftr = 5;
        randsftl = randsftlm(iii);
        randsftr = randsftrm(iii);
      else
        randsftl = 0;
        randsftr = 0;
      end
      %align records
      optcc = [];
      optcc(:,1) = STAopt(max(floor((tstbuf+1)*sps+1),1): min(floor((tedbuf-1)*sps),86400*sps), 2);
      optcc(:,2) = STAopt(max(floor((tstbuf-randsftl+1)*sps+1),1): min(floor((tedbuf-randsftl-1)*sps),86400*sps), 3);
      optcc(:,3) = STAopt(max(floor((tstbuf+randsftr+1)*sps+1),1): min(floor((tedbuf+randsftr-1)*sps),86400*sps), 4);      
      if noiseflag
        msftadd = 50;
        loffmax = 20*sps/40;
      else
        msftadd = (round(max(abs([off12ran off13ran])))+1)*sps/40;  %+1 for safety
        loffmax = 4*sps/40;
      end
      ccmid = ceil(size(optcc,1)/2);
      ccwlen = round(size(optcc,1)-2*(msftadd+1));
      ccmin = 0.01;  % depending on the length of trace, cc could be very low
      iup = 1;    % times of upsampling
      [off12con,off13con,ccali(k),iloopoff,loopoff] = constrained_cc_interp(optcc',ccmid,...
        ccwlen,msftadd,loffmax,ccmin,iup);
      % if a better alignment cannot be achieved, use 0,0
      if off12con == msftadd+1 && off13con == msftadd+1
        off12con = 0;
        off13con = 0;
        fprintf('Tremor burst %d cannot be properly aligned, double-check needed \n',k);
      end
      off1i(k,1) = 0;
      off1i(k,2) = round(off12con);
      off1i(k,3) = round(off13con);
      
      %%%what if the alignment is random or as extreme to the search boundary?
      ALIGN = 'random';
      ALIGN = 'extreme';
      ALIGN = 'best';
      if isequal(ALIGN,'random')
        spsscale = sps/40;
        loff_max = 6*spsscale;
        itry = 0;
        for i12 = -loff_max+off1i(k,2): 1: loff_max+off1i(k,2)
          for i13 = -loff_max+off1i(k,3): 1: loff_max+off1i(k,3)
            if abs(i13-i12-(off1i(k,3)-off1i(k,2))) > loff_max
              continue
            end
            itry = itry+1;
            off1ik(itry,1) = 0;
            off1ik(itry,2) = i12;
            off1ik(itry,3) = i13;
          end
        end
      elseif isequal(ALIGN,'extreme')
        spsscale = sps/40;
        loff_max = 6*spsscale;
        off1i(k,1) = 0;
        off1i(k,2) = -loff_max+off1i(k,2);
        off1i(k,3) = off1i(k,3);
      end
      
%       for itry = 1: size(off1ik,1)
%                 
%         off1i(k,:) =  off1ik(itry,:);

      %Align records
      optdat = [];  % win segment of interest
      ortdat = [];
      optdat(:, 1:2) = STAopt(max(floor(tstbuf*sps+1-overshoot),1): ...
                              min(floor(tedbuf*sps+overshoot),86400*sps), 1:2); % sta 1
      ortdat(:, 1:2) = STAort(max(floor(tstbuf*sps+1-overshoot),1): ...
                              min(floor(tedbuf*sps+overshoot),86400*sps), 1:2);      
      optdat(:, 3) = STAopt(max(floor((tstbuf-randsftl)*sps+1-overshoot)-off1i(k,2),1): ...
                            min(floor((tedbuf-randsftl)*sps+overshoot)-off1i(k,2),86400*sps), 3); % sta 2
      ortdat(:, 3) = STAort(max(floor((tstbuf-randsftl)*sps+1-overshoot)-off1i(k,2),1): ...
                            min(floor((tedbuf-randsftl)*sps+overshoot)-off1i(k,2),86400*sps), 3);
      optdat(:, 4) = STAopt(max(floor((tstbuf+randsftr)*sps+1-overshoot)-off1i(k,3),1): ...
                            min(floor((tedbuf+randsftr)*sps+overshoot)-off1i(k,3),86400*sps), 4); % sta 3
      ortdat(:, 4) = STAort(max(floor((tstbuf+randsftr)*sps+1-overshoot)-off1i(k,3),1): ...
                            min(floor((tedbuf+randsftr)*sps+overshoot)-off1i(k,3),86400*sps), 4);

      %Align the noise using the same offset
      noidat = [];  % 4-s prior to signal win
      noidat(:, 1:2) = STAopt(max(floor((tstbuf-4)*sps+1),1): min(floor((tstbuf-0)*sps),86400*sps), 1:2); % sta 1
      noidat(:, 3) = STAopt(max(floor((tstbuf-4)*sps+1)-off1i(k,2),1): ...
                            min(floor((tstbuf-0)*sps)-off1i(k,2),86400*sps), 3); % sta 2
      noidat(:, 4) = STAopt(max(floor((tstbuf-4)*sps+1)-off1i(k,3),1): ...
                            min(floor((tstbuf-0)*sps)-off1i(k,3),86400*sps), 4); % sta 3                          
                                
%       %%%segment including both the signal and noise together
%       pltdat = [];  % 4-s noise, signal win
%       pltdat(:, 1:2) = STAopt(max(floor((tstbuf-4)*sps+1-overshoot),1): min(floor(tedbuf*sps+overshoot),86400*sps), 1:2); % sta 1
%       pltdat(:, 3) = STAopt(max(floor((tstbuf-4)*sps+1-overshoot)-off1i(k,2),1): ...
%                             min(floor(tedbuf*sps+overshoot)-off1i(k,2),86400*sps), 3); % sta 2
%       pltdat(:, 4) = STAopt(max(floor((tstbuf-4)*sps+1-overshoot)-off1i(k,3),1): ...
%                             min(floor(tedbuf*sps+overshoot)-off1i(k,3),86400*sps), 4); % sta 3                          
%       %compute running CC between 3 stations
%       [ircc,rcc12] = RunningCC(pltdat(:,2), pltdat(:,3), mwlen);
%       [~,rcc13] = RunningCC(pltdat(:,2), pltdat(:,4), mwlen);
%       [~,rcc23] = RunningCC(pltdat(:,3), pltdat(:,4), mwlen);
%       ircc = ircc-overshoot;
%       rcc = (rcc12+rcc13+rcc23)/3;
%       pltdat = pltdat(overshoot+1:end-overshoot, :);  %excluding the overshoot
% 
%       %plot the overview of both the signal and noise together
%       f.fig = figure;
%       f.fig.Renderer = 'painters';
%       widin = 12;  % maximum width allowed is 8.5 inches
%       htin = 3;   % maximum height allowed is 11 inches
%       set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);
%       ncol = 1;
%       nrow = 1;
%       for isub = 1:nrow*ncol
%         f.ax(isub) = subplot(nrow,ncol,isub);
%         f.ax(isub).Box = 'on';
%         %     grid(f.ax(isub),'on');
%       end
%       pltxran = [0.08 0.92]; pltyran = [0.15 0.95];
%       pltxsep = 0.02; pltysep = 0.04;
%       axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
%       for isub = 1:nrow*ncol
%         set(f.ax(isub), 'position', axpos(isub,:));
%       end
%       
%       ax = f.ax(1); hold(ax,'on');
%       yyaxis(ax,'left');
%       plot(ax,pltdat(:,1),pltdat(:,2),'r-','linew',0.5);
%       plot(ax,pltdat(:,1),pltdat(:,3),'b-','linew',0.5);
%       plot(ax,pltdat(:,1),pltdat(:,4),'k-','linew',0.5); %,'linew',0.5
%       xran = [pltdat(1,1) pltdat(end,1)]; % x axis range in sec
%       sclseis = 1.4;
%       yran = sclseis*[-max(max(abs(optdat(:,2:end)))) max(max(abs(optdat(:,2:end))))];  % left y axis range
%       xlim(ax,xran);
%       ylim(ax,yran);
%       longticks(ax,5);
%       %plot the strongest 0.5-s arrival outside ellipse
%       ind = find(tmaxo>=xran(1) & tmaxo<=xran(2));
%       for ii = 1: length(ind)
%         barst = tmaxo(ind(ii)); % arrival of strongest 0.5-s, in sec
%         bared = tmaxo(ind(ii))+0.5; % last for 0.5 s
%         yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.03;
%         plot(ax,[barst bared],[yloc yloc],'-','linew',2.5,'color',[.6 .6 .6]);
%       end
%       %plot the strongest 0.5-s arrival inside ellipse
%       ind = find(tmaxi>=xran(1) & tmaxi<=xran(2));
%       for ii = 1: length(ind)
%         barst = tmaxi(ind(ii)); % arrival of strongest 0.5-s, in sec
%         bared = tmaxi(ind(ii))+0.5; % last for 0.5 s
%         plot(ax,[barst bared],[yloc yloc],'k-','linew',2.5);
%       end
% %       %plot the armbruster's tremor catalog outside rectangle
% %       ind = find(tarmo>=xran(1) & tarmo<=xran(2));
% %       yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.97;
% %       scatter(ax,tarmo(ind), yloc*ones(size(tarmo(ind))),10,'g','linew',1); % armbruster
% %       %plot the armbruster's tremor catalog inside rectangle
% %       ind = find(tarmi>=xran(1) & tarmi<=xran(2));
% %       scatter(ax,tarmi(ind), yloc*ones(size(tarmi(ind))),10,'g','filled');
%       %plot the bostock's LFE catalog outside rectangle
%       ind = find(tbosto>=xran(1) & tbosto<=xran(2));
%       scatter(ax,tbosto(ind), yloc*ones(size(tbosto(ind))),10,'m','linew',1); % bostock
%       %plot the bostock's LFE catalog inside rectangle
%       ind = find(tbosti>=xran(1) & tbosti<=xran(2));
%       scatter(ax,tbosti(ind), yloc*ones(size(tbosti(ind))),10,'m','filled');
%       yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.18;
%       plot(ax,[tstbuf tstbuf],ax.YLim,'--','Color',[.5 .5 .5],'linew',1);
%       text(ax,tstbuf,yloc,'start','fontsize',8);
%       plot(ax,[tedbuf tedbuf],ax.YLim,'--','Color',[.5 .5 .5],'linew',1);
%       text(ax,tedbuf,yloc,'end','fontsize',8,'HorizontalAlignment','right');
%       plot(ax,[tstbuf-4 tstbuf-4],ax.YLim,':','Color',[.5 .5 .5],'linew',1);
%       text(ax,tstbuf-4,yloc,'Automated choice of noise','fontsize',8);
%       text(ax,0.01,0.1,num2str(off1i(k,2)),'unit','normalized',...
%         'HorizontalAlignment','left','fontsize',10,'color','b'); % offset 12
%       text(ax,0.04,0.1,num2str(off1i(k,3)),'unit','normalized',...
%         'HorizontalAlignment','left','fontsize',10,'color','k'); % offset 13
%       text(ax,0.07,0.1,sprintf('%.3f',ccali(k)),'unit','normalized',...
%         'HorizontalAlignment','left','fontsize',10,'color','k'); % offset 13
%       text(ax,0.01,0.95,sprintf('Win%d--N%d',k,ninbst),'unit','normalized',...
%         'HorizontalAlignment','left','fontsize',12); % burst number
%       xlabel(sprintf('Time (s) on %s %s %s',dy,mo,yr),'FontSize',12);
%       ylabel(ax,'Amplitude','FontSize',12);
%       yyaxis(ax,'right');
%       plot(ax,ircc/sps+xran(1),rcc,'co','markersize',0.5);  % scale it to the amplitude range
% %       plot(ax,ir/sps+tstbuf,renv*1.2+0.3*yran(2),'Color',[.6 .6 .6],'linew',0.5);  % scale it to the amplitude range
%       ylim(ax,[-1 1]);
%       ylabel(ax,'Running CC','FontSize',12);
%       hold(ax,'off');
% %       orient(f.fig,'landscape');
% %       print(f.fig,'-dpdf','-bestfit',fullfile('/home/chaosong/Downloads/23.pdf'));


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
      rccpair = [rcc12 rcc13 rcc23];
      sigsta = sigsta(overshoot+1:end-overshoot, :);  %excluding the overshoot
      
%       figure
%       subplot(231)
%       histogram(rcc12); title('rcc12'); [MUHAT,SIGMAHAT] = normfit(rcc12)
%       subplot(232)
%       histogram(rcc13); title('rcc13'); [MUHAT,SIGMAHAT] = normfit(rcc13)
%       subplot(233)
%       histogram(rcc23); title('rcc23'); [MUHAT,SIGMAHAT] = normfit(rcc23)
%       subplot(234)
%       histogram(rcc12+rcc13); title('rcc12+rcc13'); [MUHAT,SIGMAHAT] = normfit(rcc12+rcc13)
%       subplot(235)
%       histogram(rcc13+rcc23); title('rcc13+rcc23'); [MUHAT,SIGMAHAT] = normfit(rcc13+rcc23)
%       subplot(236)
%       histogram(rcc12+rcc23); title('rcc12+rcc23'); [MUHAT,SIGMAHAT] = normfit(rcc12+rcc23)

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

      cc12 = xcorr(sigsta(:,1), sigsta(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
      cc13 = xcorr(sigsta(:,1), sigsta(:,3),0,'normalized');
      cc23 = xcorr(sigsta(:,2), sigsta(:,3),0,'normalized');
      if isequal(ALIGN, 'best')
        mrcc(k,1) = median(rcc);
        mcc(k,1) = (cc12+cc13+cc23)/3;
      else
        mrcc(itry,1) = median(rcc);
        mcc(itry,1) = (cc12+cc13+cc23)/3;
      end

      %if only use the mean RCC from the 2 pairs that have the highest overall CC
      [~,ind] = min([cc12 cc13 cc23]);
      rcc = sum(rccpair(:,setdiff(1:3,ind)), 2) / 2;

%       end
      
      if isequal(ALIGN,'random')
        figure
        subplot(121)
        hold on
        scatter(off1ik(:,2),off1ik(:,3),40,mrcc,'filled');
        scatter(off1ik(mrcc==max(mrcc),2),off1ik(mrcc==max(mrcc),3),30,'k','filled');
        scatter(round(off12con),round(off13con),15,'k');
        plot([off1ik(mrcc==max(mrcc),2)-10 off1ik(mrcc==max(mrcc),2)+10], ...
          [off1ik(mrcc==max(mrcc),3)-10 off1ik(mrcc==max(mrcc),3)+10],'k--');
        colormap jet
        c=colorbar;
        c.Label.String = 'Median of RCC';
        xlabel(sprintf('PGC-SSIB offset (samples at %d Hz)',sps),'FontSize',11);
        ylabel(sprintf('PGC-SILB offset (samples at %d Hz)',sps),'FontSize',11);
        axis equal
        xran = [-loff_max+round(off12con)-1 loff_max+round(off12con)+1];
        yran = [-loff_max+round(off13con)-1 loff_max+round(off13con)+1];
        xlim(xran);
        ylim(yran);
        off12m = round(range(xran)/2);
        off13m = round(range(yran)/2);
        xbndlow = [offxran(1):offxran(2), offxran(2)*ones(1,off13m)];
        xbndupp = [offxran(2)-1:-1:offxran(1), offxran(1)*ones(1,off13m)];
        ybndlow = [offyran(1)*ones(1, off12m), offyran(1):offyran(2)];
        ybndupp = [offyran(2)*ones(1, off12m-1), offyran(2):-1:offyran(1)];
        xbnd = [xbndlow xbndupp]';
        ybnd = [ybndlow ybndupp]';
        plot(ax, xbnd,ybnd,'-','Color',[.4 .4 .4],'linew',2);

        box on

        subplot(122)
        hold on
        scatter(off1ik(:,2),off1ik(:,3),40,mcc,'filled');
        scatter(off1ik(mcc==max(mcc),2),off1ik(mcc==max(mcc),3),30,'k','filled');
        scatter(round(off12con),round(off13con),15,'k');
        plot([off1ik(mcc==max(mcc),2)-10 off1ik(mcc==max(mcc),2)+10], ...
          [off1ik(mcc==max(mcc),3)-10 off1ik(mcc==max(mcc),3)+10],'k--');
        colormap jet
        c=colorbar;
        c.Label.String = 'Overall CC';
        xlabel(sprintf('PGC-SSIB offset (samples at %d Hz)',sps),'FontSize',11);
        ylabel(sprintf('PGC-SILB offset (samples at %d Hz)',sps),'FontSize',11);
        axis equal
        xlim(xran);
        ylim(yran);
        plot(ax, xbnd,ybnd,'-','Color',[.4 .4 .4],'linew',2);
        box on
      end
        
      %%%Is it true that the coherence between 2-3 is the highest among 3 pairs?
      %%%---Not really, or align with 1 offset upon the whole long window might not support it
%       [f] = plt_rcccat(rccpair,sps);

      %%%what about the indices of times when there are 4-s in-bound detections; no detections;
      %%%out-bound detections, etc?
      %in-bound detections
      inddeti = [];
      for jj = 1:ninbst(k,1)
        inddeti = [inddeti round((max(tcnti(indtmaxi(jj))-2, tst-2)-tstbuf)*sps+1 -mwlen/2) : ...
                           round((min(tcnti(indtmaxi(jj))+2, ted+2)-tstbuf)*sps -mwlen/2) ];
      end
      inddeti = unique(inddeti);  % in case of duplicates due to overlapping windows
      %out-bound detections
      indtmaxo = find(tmaxo>=tst-0.1 & tmaxo<=ted+0.1); %+-0.1 s in case of resolution while saving to file
      inddeto = [];
      if ~isempty(indtmaxo)
        tcnto = hfdayo(:, 15);  % the center of detecting win is the 15th col
        for jj = 1: length(indtmaxo)
          inddeto = [inddeto round((max(tcnto(indtmaxo(jj))-2, tst-2)-tstbuf)*sps+1 -mwlen/2) : ...
                             round((min(tcnto(indtmaxo(jj))+2, ted+2)-tstbuf)*sps -mwlen/2) ];
        end
        inddeto = unique(inddeto);  % in case of duplicates due to overlapping windows
      end
      %no detections
      indnodet = setdiff(setdiff(1:length(rcc), inddeti), inddeto);
    
      %the median rcc as the stopping threshold in iterations
%       medrcc = median(rcc([indnodet inddeto]));
      medrcc = median(rcc(indnodet));
      
      %%%what if we don't use rcc as the weighting, ie., set to a constant of 1
%       rcc = ones(size(sigsta,1),1);
      
      %%
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
%             fixthresh = 1.2352e-01; % if use 3-pair rcc and normalize templates, sta 1
%             fixthresh = 1.7766e-01; % if use 2-pair rcc and normalize templates
%             fixthresh = 6.3908e-02; % if use 3-pair rcc and do not normalize templates
            fixthresh = 9.1922e-02; % if use 2-pair rcc and do not normalize templates
          elseif ista == 2
%             fixthresh = 7.3467e-02; % if use 3-pair rcc and normalize templates, sta 2
%             fixthresh = 9.5637e-02; % if use 2-pair rcc and normalize templates
%             fixthresh = 2.9388e-02; % if use 3-pair rcc and do not normalize templates
            fixthresh = 3.8256e-02; % if use 2-pair rcc and do not normalize templates
          elseif ista == 3
%             fixthresh = 1.0873e-01; % if use 3-pair rcc and normalize templates, sta 3
%             fixthresh = 1.3470e-01; % if use 2-pair rcc and normalize templates
%             fixthresh = 4.2047e-02; % if use 3-pair rcc and do not normalize templates
            fixthresh = 5.2089e-02; % if use 2-pair rcc and do not normalize templates
          end
          [sigdecon(:,ista),pred(:,ista),res,dresit,mfitit,ampit{ista},nit,fighdl] = ...
            iterdecon_fixthresh(sig,wlet,rcc,noi,fixthresh,dt,twlet,width,dres_min,mfit_min,nit_max,nimp_max,fpltit,fpltend,fpltchk);
        else
          [sigdecon(:,ista),pred(:,ista),res,dresit,mfitit,ampit{ista},nit,fighdl] = ...
            iterdecon(sig,wlet,rcc,noi,dt,twlet,width,dres_min,mfit_min,nit_max,nimp_max,fpltit,fpltend,fpltchk);
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
      
      %% Group nearest impulses from different stations into pairs
      spsscale = sps/40;
      loff_max = 6*spsscale;  % maximum allowed shift is 4 samples at 40 Hz, it is also true from the synthetic distribution
      %note the output 'impindep' gives the arrival index of impulse at each station, after
      %alignment based upon the entire window 'off1i', and the last three cols are the arrival time
      %difference, NOT the true location yet! 
      refsta = 1;
%       [impindep,imppairf,indpair] = groupimptriplets(sigdecon,rcc,loff_max,'wtamp',refsta);
      [impindep,imppairf,indpair] = groupimptripdecon(sigdecon,ampit,rcc,loff_max,refsta);
      
      nsrc(iii) = size(impindep,1);
      end
    
    figure
    histogram(nsrc);
    text(0.05,0.9,sprintf('med=%d',median(nsrc)),'Units','normalized');
    xlabel('number of sources');
    ylabel('counts');
    
    keyboard

      %%%plot the individually deconvolved impulses and the grouping result                 
%       [f] = plt_groupimptriplets(sigdecon,impindep,stas,ircc,rcc);
     
      %note here 'impindepst' inherits the first 6 cols from 'impindep', but the last three cols 
      %are adjusted from arrival time difference to the true location offset accounting for the best
      %alignment upon the entire window that is also used in grouping!
      impindepst = sortrows(impindep,1);
      impindepst(:,7:8) = impindepst(:,7:8)+repmat([off1i(k,2) off1i(k,3)],size(impindepst,1),1); %account for prealignment
      
      %%%plot the scatter of offsets, accounting for prealignment offset, == true offset
      xran = [-loff_max+off1i(k,2)-1 loff_max+off1i(k,2)+1];
      yran = [-loff_max+off1i(k,3)-1 loff_max+off1i(k,3)+1];
      offxran = [-loff_max+off1i(k,2) loff_max+off1i(k,2)];
      offyran = [-loff_max+off1i(k,3) loff_max+off1i(k,3)];
      cran = [0 lsig];
      f1.fig = figure;
      f1.fig.Renderer = 'painters';
      ax1=gca;
      [ax1,torispl,mamp] = plt_decon_imp_scatter(ax1,impindepst,xran,yran,cran,offxran,offyran,...
        sps,75,'mean','tori');
      scatter(ax1,off1i(k,2),off1i(k,3),20,'ks','filled','MarkerEdgeColor','k');
      title(ax1,'Independent, grouped');
      
keyboard
      %% check the difference by grouping using different stations as the reference station
%       spsscale = sps/40;
%       loff_max = 4*spsscale;  % maximum allowed shift is 4 samples at 40 Hz, it is also true from the synthetic distribution
%       %note the output 'impindep' gives the arrival index of impulse at each station, after
%       %alignment based upon the entire window 'off1i', and the last three cols are the arrival time
%       %difference, NOT the true location yet! 
% %       aa = groupimptripdeconOLD(sigdecon,ampit,rcc,loff_max);
%       aa = groupimptriplets(sigdecon,rcc,loff_max,'wtamp',1);
%       
%       for refsta = 1:3
% %         bb = groupimptripdecon(sigdecon,ampit,rcc,loff_max,refsta);
%         bb = groupimptriplets(sigdecon,rcc,loff_max,'wtamp',refsta);
%         size(bb,1)
%         indsame = find(ismember(bb,aa,'rows'));
%         inddiff = find(~ismember(bb,aa,'rows'));
%         perca = length(indsame)/size(aa,1)*100
%         percb = length(indsame)/size(bb,1)*100
%       end      
            
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
      length(ind)

      %use function 'removesecondarysrc' to remove the smaller amplitude, secondary triplets that 
      %are too close in time to large triplets
      minsrcsep = 20;   % min allowed separation between deconvolved peaks
      [pkindepsave,indremove] = removesecondarysrc(pkindep,sigsta);
%       [pkindepsave,indremove] = removesecondarysrc(pkindep,sigsta,minsrcsep);
      
      %remove the secondary sources from the grouped result
      impindep(indremove, :) = [];

      %% plot the scatter of sources in terms of offsets, accounting for prealignment offset
      impindepst = sortrows(impindep,1);
      impindepst(:,7:8) = impindepst(:,7:8)+repmat([off1i(k,2) off1i(k,3)],size(impindepst,1),1); %account for prealignment
      xran = [-loff_max+off1i(k,2)-1 loff_max+off1i(k,2)+1];
      yran = [-loff_max+off1i(k,3)-1 loff_max+off1i(k,3)+1];
      offxran = [-loff_max+off1i(k,2) loff_max+off1i(k,2)];
      offyran = [-loff_max+off1i(k,3) loff_max+off1i(k,3)];
      cran = [0 lsig];
      %%%plot the scatter of offsets, accounting for prealignment offset, == true offset
      f1.fig = figure;
      f1.fig.Renderer = 'painters';
      ax1=gca;
      [ax1,torispl,mamp] = plt_decon_imp_scatter(ax1,impindepst,xran,yran,cran,offxran,offyran,...
        sps,75,'mean','tori');
      scatter(ax1,off1i(k,2),off1i(k,3),20,'ks','filled','MarkerEdgeColor','k');
      title(ax1,'Independent, grouped, no secondary sources');
      
      %% plot the scatter of sources in terms of rela locations
      xran = [-3 3];
      yran = [-3 3];
      cran = [0 lsig/sps];
      f2.fig = figure;
      f2.fig.Renderer = 'painters';
      ax2=gca;
      [ax2] = plt_decon_imp_scatter_space(ax2,impindepst,xran,yran,cran,offxran,...
        offyran,sps,75,ftrans,'mean','tori');
      plot(ax2,xcut,ycut,'k-','linew',2);
      title(ax2,'Independent, grouped, no secondary sources');
      
      %% Compare deconvolved positive peaks with waveform peaks 
      %%plot the positive peaks indicated by the grouped triplets, see if they indeed match the
      %%peaks of the data at each station, could be helpful which triplets are minor
      fpltremv = 1;   % plot or not the secondary sources that are removed as well
      [f] = plt_deconpk_sigpk(sigsta,pkindep,indremove,fpltremv);
      xlim(f.ax(1),[0 2000]);
      xlim(f.ax(2),[0 2000]);
      xlim(f.ax(3),[0 2000]);
            
      %% separation in arrival time between deconvolved positive peaks
      %%%Plot the separation in time between these preserved positive peaks after removing the
      %%%secondary ones, to see if they can be too close to each other
      [f,tsep,pkist,indpk,indreverse] = plt_tsep_deconpk(pkindepsave,sps);   
      
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
%       [f] = plt_srcalldistCDF(impindepstst,implocst,torisplst,1*sps,sps);

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

      %% linear regression of offset change to know how fast the centroid of sources are moving
      [f,fitobj12,fitobj13,dyfit12,dyfit13] = plt_deconlfit(impindepst,torispl);

      %%%estimate a preliminary subwin length that the offset could change by 1 sample
      madoff = max(abs([dyfit12; dyfit13]));
      subwseclfit(k) = tlenbuf/madoff;

      %% what's the corresponding rcc at source arrival, associated with the offset
      [f] = plt_srcrccwithoffset(impindepst,rccpair);
      
%       %% what is the distribution of 4-s detections inside the same window?
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
      
      %% what about the orthogonal component?
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
      text(0.95,0.75,sprintf('%.1f%%; %.1f%%; %.1f%%',...
        (norm(sigsta(:,1))-norm(resgrp(:,1)))/norm(sigsta(:,1))*100,...
        (norm(sigsta(:,2))-norm(resgrp(:,2)))/norm(sigsta(:,2))*100,...
        (norm(sigsta(:,3))-norm(resgrp(:,3)))/norm(sigsta(:,3))*100),...
        'Units','normalized','HorizontalAlignment','right');
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
      text(0.95,0.75,sprintf('%.1f%%; %.1f%%; %.1f%%',...
        (norm(sigstaort(:,1))-norm(resgrport(:,1)))/norm(sigstaort(:,1))*100,...
        (norm(sigstaort(:,2))-norm(resgrport(:,2)))/norm(sigstaort(:,2))*100,...
        (norm(sigstaort(:,3))-norm(resgrport(:,3)))/norm(sigstaort(:,3))*100),...
        'Units','normalized','HorizontalAlignment','right');


      %% plot the comparision between the real data and prediction using the grouped impulses
      %compute running CC between 3 stations
      [irccp,rccp12] = RunningCC(predgrp(:,1), predgrp(:,2), mwlen);
      [~,rccp13] = RunningCC(predgrp(:,1), predgrp(:,3), mwlen);
      [~,rccp23] = RunningCC(predgrp(:,2), predgrp(:,3), mwlen);
      irccp = irccp-overshoot;
      rccp = (rccp12+rccp13+rccp23)/3;  % use averge
      
      [f] = plt_sigpred(sigsta,predgrp,impindepst,sps,tmaxi,tmaxo,tbosti,tbosto,...
        ircc,rcc,irccp,rccp,mwlen,tstbuf,tedbuf,dy,mo,yr);
      
      %% A summary figure
      %%%compose a summary figure for each data win and save it, so that we will have a feeling for
      %%%all migrations
      [f3] = plt_decon_summary(sigsta,resgrp,pkindepsave,impindepst,sps,tmaxi,tmaxo,tbosti,tbosto,...
        ircc,rcc,mwlen,off1i,loff_max,tstbuf,tedbuf,dy,mo,yr,k,ftrans,xcut,ycut);
      h=supertit(f3.ax(1), sprintf('Window %s, %d, %.1f-%.1f s, Passband: %.1f-%.1f Hz, sps: %d',...
        num2zeropadstr(k, 3),date,tstbuf,tedbuf,losig,hisig,sps),10);
      movev(h,-0.05);

      %save figure
      fignm{k,1} = sprintf('decon002_%.1fsFig%s_%.1f-%.1f.pdf',...
        mwlen/sps,num2zeropadstr(k,3),losig,hisig);
      print(f3.fig,'-dpdf',fullfile(rstpath, '/FIGS',fignm{k,1}));


%       keyboard
    end
    
  end
  
end








