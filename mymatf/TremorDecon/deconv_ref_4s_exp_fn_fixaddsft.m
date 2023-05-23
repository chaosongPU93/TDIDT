% deconv_ref_4s_exp_4thsta_fn.m
% function rststruct = deconv_ref_4s_exp_fn_fixaddsft(idxbst,normflag,rccmwsec,addsft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A modified version from 'deconv_ref_4s_exp_4thsta_fn' to handle the case 
% where a customized additional shift defined in 'addsft' can be applied to 
% data traces after aligning to our best. However, the decon stopping threshold
% remains still the same. Sp basically, you align the data as before, then 
% invoke recompute RCC, etc after applying the additional shifts.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/01/23
% Last modified date:   2023/01/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e   % Set the format to 5-digit floating point
clear
clc
close all

%% for easy testing
defval('idxbst',181); %global indices of bursts to run 
defval('normflag',0); %whether to normalize templates
defval('rccmwsec',0.5); %moving win len in sec for computing RCC
defval('addsft',[-6 0 0]); %moving win len in sec for computing RCC

rccflag = 1; %1 means RCC weighting is used

if addsft(1)~= 0
  addsft = addsft-addsft(1); 
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

%%%load the empirical model param of time offset 14 for each addtional stations
off14mod = load(strcat(rstpath, '/MAPS/timeoff_planefitparam_4thsta_160sps'),'w+');

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
STA = detrend(ccstack);

ccstackort = [];
for ista = 1: nsta
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
        num2str(templensec), 's_', 'BBCCS_', 'ort_Nof_Non_Chao_catnew');
    ccstackort(:,ista) = load(fname);
end
STAort = detrend(ccstackort);

%flag of normalization
% normflag = 0;

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

%%%plot the unfiltered and filtered templates
% plt_templates(green,greenf,stas,greenort,greenfort,lowlet,hiwlet,sps);

%just the filtered templates
% plt_templates_bp(greenf,stas,lowlet,hiwlet,sps);


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

%%%moving window length in samples for running CC, envelope, etc.
%standard window length is about 0.5s, this is about the visual duration of the filtered and unfiltered
%template, although in fact to include as least one cycle of the main dipole of template
rccmwlen=rccmwsec*sps;
% rccmwlen=sps/2;
% rccmwlen=sps;

k = 0;  % counting the burst windows

off1iwk = cell(size(trange,1),1);  % the best alignment between sta 2, 3 wrt 1 for all subwins and all burst wins
off1ic = zeros(size(trange,1),nsta);  % single best alignment 'computed' between ALL stas wrt 1 for entire win
off1i = zeros(size(trange,1),nsta);  % single best alignment 'actually used' between ALL stas wrt 1 for entire win
off14pred = zeros(size(trange,1),nsta-3); %empirical pred of off14 from plane fit given single best alignment
ccali = zeros(size(trange,1),1);  % CC value using the best alignment
ccaliwk = cell(size(trange,1),1); % CC value using the best alignment for all subwins and all burst wins 
subwsec = zeros(size(trange,1),1);  % subwin length in sec used in practice
subwseclfit = zeros(size(trange,1),1);  % subwin length from linear fitting, ie, time for 1-sample offset change

% %empirically determined indices of bursts fall into local day times (noisier)
% %or night times (quieter)
% inbst = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,...
%   34,35,56,57,58,59,60,61,62,70,71,72,73,74,75,76,77,78,79,80,81,82,110,111,112,113,114,115,116,117,...
%   118,119,120,142,143,144,145,146,147,148,149,150,151,152,153,154,172,173,174,175];
% idbst = [36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,63,64,65,66,67,68,69,83,84,85,...
%   86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,121,122,123,124,...
%   125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,155,156,157,158,159,160,161,...
%   162,163,164,165,166,167,168,169,170,171,176,177,178,179,180,181,182,183,184,185,186,187,188,189,...
%   190,191,192,193,194,195];
% 
% %indices of bursts whose wlet-sig cc between the opt and ort components are very similar, above 75
% %percentile
% ihighoo = [20,23,59,60,80,113,116,120,134,189,194];
% 
% %indices of bursts whose sig cc between sta 1/2/3 are very high, above 90 percentile
% ihicc123 = [1,3,6,7,8,24,56,71,75,77,81,83,93,102,114,116,132,145,149,185];
% ihicc123n = intersect(inbst,ihicc123);
% 
% indtest = [18,21,22,23];

%store all target windows
srcamprall = [];  %src amp ratio 
lndevsrcamprall = []; %linear deviation from median src amp ratio
lgdevsrcamprall = []; %log deviation from median src amp ratio
rcccatsrcall = [];  %mean concat RCC among trio and RCC14 at src arrival 
rccpairsrcall = []; %concat RCC for trio sta pairs at src arrival 
psrcampsall = []; %positive scaled src amp  
nsrcampsall = []; %negative scaled src amp  
psrcamprsall = [];  %positive scaled src amp ratio  
nsrcamprsall = [];  %negative scaled src amp ratio  
clppkhtwfall = [];  %closest pos peak height of waveform
clnpkhtwfall = [];  %closest neg peak height of waveform
clppkwfsepall = [];  %closest pos peak separation of waveform
clnpkwfsepall = [];  %closest neg peak separation of waveform
ppkwfsepmed = []; %median of the pos peak separation of waveform
ppkwfsepmod = []; %mode of the pos peak separation of waveform
npkwfsepmed = []; %median of the neg peak separation of waveform
npkwfsepmod = []; %mode of the neg peak separation of waveform
pred4offtrall = [];  %difference in arrival from prediction at 4th sta  
impindepstall = [];  %all sources' arrivals and amplitudes

%secondary arrivals removed, decon impulse tarvl separation, spatial distance, etc.
tsepall = []; 
dtorinn1all = [];
distorinn1all = [];
dtorinn2all = [];
distorinn2all = [];
dtorinn3all = [];
distorinn3all = [];
dtoripropall = [];
distoripropall = [];
distoriortall = [];
dtarvlnn1all = [];
distarvlnn1all = [];
dtarvlnn2all = [];
distarvlnn2all = [];
dtarvlnn3all = [];
distarvlnn3all = [];
dtarvlpropall = [];
distarvlpropall = [];
distarvlortall = [];

%4th station checked, decon impulse tarvl separation, spatial distance, etc.
tsep4thall = []; 
dtorinn14thall = [];
distorinn14thall = [];
dtorinn24thall = [];
distorinn24thall = [];
dtorinn34thall = [];
distorinn34thall = [];
dtoriprop4thall = [];
distoriprop4thall = [];
distoriort4thall = [];
dtarvlnn14thall = [];
distarvlnn14thall = [];
dtarvlnn24thall = [];
distarvlnn24thall = [];
dtarvlnn34thall = [];
distarvlnn34thall = [];
dtarvlprop4thall = [];
distarvlprop4thall = [];
distarvlort4thall = [];

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
    
%     %read vertical components too, as we want to have a sense of the SNR at the orthogonal and vertical
%     %components
%     [STAvert,~,~,fileflag] = rd_daily_bpdata(year,jday,prename,stas,PERMSTA,POLSTA,...
%       PERMROTS,POLROTS,sps,losig,hisig,npo,npa,[],[],[],[],'Z');

    
%     keyboard
%%
%     for j = 14: size(rangetemp,1)  
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
      overshoot = rccmwlen/2;
%       overshoot = 0;
      
      %%%%2022/09/26, obtain all information for data before decon, so that you know the threshold
      %%%%being used, although this was used mainly for noise experiment, we still use it here 
      %chop a record segment
      optseg = STAopt(max(floor(tstbuf*sps+1-overshoot-msftaddm),1): ...
        min(floor(tedbuf*sps+overshoot+msftaddm),86400*sps), :); % sta 1
      ortseg = STAort(max(floor(tstbuf*sps+1-overshoot-msftaddm),1): ...
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
      off1iwk{k} = off1iw;
      ccaliwk{k} = ccaliw;

      %if only use the mean RCC from the 2 pairs that have the highest overall CC
      [~,ind] = min(sum(ccwpair,1));
      rcccat = sum(rccpaircat(:,setdiff(1:3,ind)), 2) / 2;
      
      %if choose not to use RCC weighting; for easier comparison
      if ~rccflag
        rcccat = ones(size(rcccat));
        rcc1icat = ones(size(rcc1icat));
      end


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
%       off1ia(k,4:end) = off14pred(k,:);
      
      %for real data, USE the best whole-win alignment before decon
      off1i(k,1:3) = off1ic(k,1:3); 
      
      %Choice to make upon the actual-used alignment at 4th stations
      %for data case, DO align!
      align14flag = 1;  
      if align14flag
        off1i(k,4:end) = off1ic(k,4:end); %if you trust overall alignment at 4th stations
      else
        off1i(k,4:end) = zeros(1,nsta-3); %if you don't 
      end

      %%%Align and compute the RCC based on the entire win, and take that as the input signal!      
      optdat = [];  % win segment of interest
      ortdat = [];
      optdat(:, 1) = optseg(1+msftaddm: end-msftaddm, 1); % time column
      ortdat(:, 1) = ortseg(1+msftaddm: end-msftaddm, 1);
      for ista = 1: nsta 
        optdat(:, ista+1) = optseg(1+msftaddm-off1i(k,ista): end-msftaddm-off1i(k,ista), ista+1);
        ortdat(:, ista+1) = ortseg(1+msftaddm-off1i(k,ista): end-msftaddm-off1i(k,ista), ista+1);
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
      
      %for ort. comp
      sigstaort = zeros(size(ortdat,1), nsta);
      for ista = 1:nsta
        tmp = ortdat(:,ista+1); %best aligned, filtered
        tmp = detrend(tmp);
        sigstaort(:,ista) = tmp;
      end
      [irccort,rcc12,rcc13,rcc23] = RunningCC3sta(sigstaort,rccmwlen);
      irccort = irccort-overshoot;
      rccort = (rcc12+rcc13+rcc23)/3;
      sigstaort = detrend(sigstaort(overshoot+1:end-overshoot, :));  %excluding the overshoot
      
%       figure
%       hold on
%       box on
%       plot(ircc,rcc,'k-')
%       plot(ircccat,rcccat,'r-');
%       ylim([-1 1]);      
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
%           wtcoef = rcc1i(pkind, ista-3).* pkhgt;
          wtcoef = rcc1icat(pkind, ista-3).* pkhgt;
        end
        fixthresh(ista) = median(wtcoef);  % median of the weighted master CC, could be percentile?
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%Apply additional shift to traces as requested in input
      ircccat = [];   % concatenated indices of RCC
      irccran = zeros(nwin,2);  % start and end indices (range) of RCC of all subwins
      rcccat = [];  % average concatenated RCC
      rcc1icat = [];  % concatenated RCC between sta 1 and 4
      rccpaircat = [];  % concatenated RCC between each station pair, order is 12, 13, 23
      ccwpair = []; % 0-lag overall cc of each subwin, between each station pair, order is 12, 13, 23
      ccw1i = []; % same as above, but between sta 1 and 4
      
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
        
        off1iw(iwin,1) = off1iw(iwin,1)+addsft(1);
        off1iw(iwin,2) = off1iw(iwin,2)+addsft(2);
        off1iw(iwin,3) = off1iw(iwin,3)+addsft(3);
        
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
      off1iwk{k} = off1iw;
      ccaliwk{k} = ccaliw;
      
      %if only use the mean RCC from the 2 pairs that have the highest overall CC
      [~,ind] = min(sum(ccwpair,1));
      rcccat = sum(rccpaircat(:,setdiff(1:3,ind)), 2) / 2;
      
      %if choose NOT to use RCC weighting; for easier comparison
      if ~rccflag
        rcccat = ones(size(rcccat));
        rcc1icat = ones(size(rcc1icat));
      end
      
      %%%Is it true that the coherence between 2-3 is the highest among 3 pairs?
      %%%---Yes, for the concatenated rcc
      %       [f] = plt_rcccat(rccpaircat,sps);
      
      %%%obtain a single best alignment based on the entire win
      off1ic(k,1) = off1ic(k,1)+addsft(1);
      off1ic(k,2) = off1ic(k,2)+addsft(2);
      off1ic(k,3) = off1ic(k,3)+addsft(3);
      
      %for real data, USE the best whole-win alignment before decon
      off1i(k,1:3) = off1ic(k,1:3);
      
      %Choice to make upon the actual-used alignment at 4th stations
      %for data case, DO align!
      %       align14flag = 1;
      if align14flag
        off1i(k,4:end) = off1ic(k,4:end); %if you trust overall alignment at 4th stations
      else
        off1i(k,4:end) = zeros(1,nsta-3); %if you don't
      end
      
      %%%Align and compute the RCC based on the entire win, and take that as the input signal!
      optdat = [];  % win segment of interest
      ortdat = [];
      optdat(:, 1) = optseg(1+msftaddm: end-msftaddm, 1); % time column
      ortdat(:, 1) = ortseg(1+msftaddm: end-msftaddm, 1);
      for ista = 1: nsta
        optdat(:, ista+1) = optseg(1+msftaddm-off1i(k,ista): end-msftaddm-off1i(k,ista), ista+1);
        ortdat(:, ista+1) = ortseg(1+msftaddm-off1i(k,ista): end-msftaddm-off1i(k,ista), ista+1);
      end
      
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
      
      %for ort. comp
      sigstaort = zeros(size(ortdat,1), nsta);
      for ista = 1:nsta
        tmp = ortdat(:,ista+1); %best aligned, filtered
        tmp = detrend(tmp);
        sigstaort(:,ista) = tmp;
      end
      [irccort,rcc12,rcc13,rcc23] = RunningCC3sta(sigstaort,rccmwlen);
      irccort = irccort-overshoot;
      rccort = (rcc12+rcc13+rcc23)/3;
      sigstaort = detrend(sigstaort(overshoot+1:end-overshoot, :));  %excluding the overshoot
      
      %       figure
      %       hold on
      %       box on
      %       plot(ircc,rcc,'k-')
      %       plot(ircccat,rcccat,'r-');
      %       ylim([-1 1]);
      
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
        
        [sigdecon(:,ista),pred(:,ista),res,dresit,mfitit,ampit{ista},nit,fighdl] = ...
          iterdecon(sig,wlet,rcccat,noi,fixthresh(ista),dt,twlet,width,dres_min,...
          mfit_min,nit_max,nimp_max,fpltit,fpltend,fpltchk);
       
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
      [impindep,imppairf,indpair,sharp] = groupimptripdecon_ref(sigdecon,ampit,irccran,rcccat,...
        off1i(k,1:3),off1iw(:,1:3),loff_max,refsta);
      
      %plot the sharpness of grouped peaks in res-wlet CC 
      f=plt_srcsharpness(sharp);
%       
%       nsrc(iii) = size(impindep,1);
%       end
%     
%     figure
%     histogram(nsrc);
%     text(0.05,0.9,sprintf('med=%d',median(nsrc)),'Units','normalized');
%     xlabel('number of sources');
%     ylabel('counts');
%     
      
%       %%%plot the individually deconvolved impulses and the grouping result
%       [f] = plt_groupimptriplets(sigdecon,impindep,stas,ircccat,rcccat);
      
      %note here 'impindepst' inherits the first 6 cols from 'impindep', but the last three cols 
      %are adjusted from arrival time difference to the true location offset accounting for the best
      %alignment upon each subwin that is also used in grouping!
      impindep(:,7:8) = impindep(:,7:8)+repmat([off1i(k,2) off1i(k,3)],size(impindep,1),1); %account for prealignment
      impindepst = sortrows(impindep,1);

      %%
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
%           impindepst(:,8),impindepst(:,1),off1ia(k,ista));
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
%       
%       %%%waveform prediction via conv of impulses and template at each sta, and res reduction
%       [f2,predgrp,resgrp,predgrpl,resgrpl]=predsig_conv_imptemp(sigsta,optdat,impindepst,...
%         greenf,zcrosses,overshoot,stas,1);


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
      
      %REMOVE the secondary sources from the grouped result
      impindep(indremove, :) = [];
      ppkindep(indremove, :) = [];
      npkindep(indremove, :) = [];
      sharp(indremove, :) = [];
      impindepst = sortrows(impindep,1);
      nsrcraw(iii,1) = size(impindepst,1);  % number of sources AFTER removing 2ndary 

      if ~isempty(impindepst)

      %plot the sharpness of grouped peaks in res-wlet CC
      f=plt_srcsharpness(sharp);

      %% plot the scatter of sources in terms of offsets, accounting for prealignment offset
      span = max(range(off1iw(:,2))+2*loff_max, range(off1iw(:,3))+2*loff_max);
      xran = [round(mean(minmax(off1iw(:,2)'))-span/2)-1, round(mean(minmax(off1iw(:,2)'))+span/2)+1];
      yran = [round(mean(minmax(off1iw(:,3)'))-span/2)-1, round(mean(minmax(off1iw(:,3)'))+span/2)+1];
      cran = [0 lsig];
      %%%plot the scatter of offsets, accounting for prealignment offset, == true offset
      f1.fig = figure;
      f1.fig.Renderer = 'painters';
      ax1=gca;
      [ax1,torispl,mamp,xbndcvhl,ybndcvhl] = plt_decon_imp_scatter_ref(ax1,impindepst,xran,yran,...
        cran,off1iw,loff_max,sps,50,'mean','tori','comb');
      scatter(ax1,off1i(k,2),off1i(k,3),20,'ks','filled','MarkerEdgeColor','k');
      title(ax1,'Independent, grouped, no secondary sources');
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
%       title(ax1,'Independent, grouped, no secondary sources, shifted');

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
            
      %% separation in arrival time between deconvolved positive peaks
      %%Plot the separation in time between these preserved positive peaks after removing the
      %%secondary ones, to see if they can be too close to each other
      if size(ppkindepsave,1) > 1
        [f,tsep,pkist,indpk,indreverse] = plt_tsep_deconpk(ppkindepsave,sps);
        tsepall = [tsepall; tsep];
        close(f.fig);
      end
%       median(tsep)
% keyboard

      %% signal + zoom-in + map locations + some reference symbols 
      xzoom = [5 30];
      [f] = plt_agu2022abstractv3(greenf(:,1:3),sigsta(:,1:3),impindepst,sps,xzoom,off1iw,loff_max,...
        tstbuf,dy,mo,yr,ftrans);
      text(f.ax(3),0.98,0.95,'Secondary removed','HorizontalAlignment','right',...
        'Units','normalized','FontSize',10);
      text(f.ax(4),0.98,0.95,'Zoom-in','HorizontalAlignment','right',...
        'Units','normalized','FontSize',10);

      ax=f.ax(3);
      hold(ax,'on');
      xran1 = [-4 4];
      yran1 = [-4 4];
      wtmax = prctile(mean(impindepst(:,[2 4 6]),2),95); %use percentile in case
      scatter(ax,xran1(1)+0.1*range(xran1),yran1(2)-0.05*range(yran1),35,'w','filled',...
        'MarkerEdgeColor',[.5 .5 .5]);
      %       text(ax,0.02,0.9,strcat({'amp.\geq'},sprintf('%.1f',wtmax)),'Units','normalized',...
      %         'HorizontalAlignment','left','FontSize',8);
      text(ax,0.02,0.9,'amplitude ','Units','normalized',...
        'HorizontalAlignment','left','FontSize',8);
      text(ax,0.02,0.85,'\geq 95th prctile','Units','normalized',...
        'HorizontalAlignment','left','FontSize',8);
      angrmse = 140;
      [rotx, roty] = complex_rot(0,1,-angrmse);
      xvect = [0.5-rotx 0.5+rotx];
      yvect = [-2.5-roty -2.5+roty];
      drawArrow(ax,xvect,yvect,xran1,yran1,'linewidth',1);
      text(ax,0.62,0.18,strcat(num2str(angrmse),'$^{\circ}$'),'FontSize',10,...
        'unit','normalized','interpreter','latex');
      hold(ax,'off');
      
      %%%add the indication of the reasonable physical source size
      ax=f.ax(4);
      hold(ax,'on');
      x0 = 2;
      y0 = -1;
      Vs = 3;   % S-wave speed
      %       Vprop = 0.8*Vs;   % reasonable rupture propagation speed
      Vprop = Vs;   % max rupture propagation speed
      td = 40/sps;  % estimated from the bin where the most tsep falls in, 32-48-sample bin
      %         td = median(tsep(:));
      srcsz = Vprop*td;
      radi = srcsz/2;
      scatter(ax,x0,y0,2,'ko','filled');
      [x, y] = circle_chao(x0,y0,radi,0.1);
      plot(ax,x,y,'-','Color','k','linew',1);
      text(ax,0.55,0.25,'max. reasonable','HorizontalAlignment','left',...
        'Units','normalized','FontSize',10,'interpreter','latex');
      text(ax,0.55,0.2,'source size','HorizontalAlignment','left',...
        'Units','normalized','FontSize',10,'interpreter','latex');
      text(ax,0.55,0.15,strcat('$V_{s} \cdot t_{d} = 3 \cdot 0.25$'),'HorizontalAlignment','left',...
        'Units','normalized','FontSize',10,'interpreter','latex');
      
      %%%add the indication of the error ellipse
      mmmax = 6;
      nnmax = 6;
      erroff1 = zeros((mmmax*2+1)^2,2);
      kk = 0;
      for mm = -mmmax:1:mmmax
        for nn = -nnmax:1:nnmax
          kk = kk+1;
          erroff1(kk,:) = [mm nn];
        end
      end
      [errloc1, ~] = off2space002(erroff1,sps,ftrans,0);
      
      F1 = scatteredInterpolant(erroff1(:,1),erroff1(:,2),errloc1(:,1),'linear');
      F2 = scatteredInterpolant(erroff1(:,1),erroff1(:,2),errloc1(:,2),'linear');
      
      erroff = 3*[0 2; 1 sqrt(3); sqrt(2) sqrt(2); sqrt(3) 1;
        2 0; sqrt(3) -1; sqrt(2) -sqrt(2); 1 -sqrt(3);
        0 -2; -1 -sqrt(3); -sqrt(2) -sqrt(2); -sqrt(3) -1;
        -2 0; -sqrt(3) 1; -sqrt(2) sqrt(2); -1 sqrt(3);
        0 2;];
      errloc = [];
      errloc(:,1) = F1(erroff(:,1),erroff(:,2));
      errloc(:,2) = F2(erroff(:,1),erroff(:,2));
      
      errloc = errloc(:,1:2)+[-2.5 1.5];
      plot(ax,errloc(:,1),errloc(:,2),'-','Color',[.4 .4 .4],'linew',1);
      text(ax,0.05,0.9,strcat('$\pm6$','-sample'),'HorizontalAlignment','left',...
        'Units','normalized','FontSize',10,'interpreter','latex');
      text(ax,0.05,0.85,'error contour','HorizontalAlignment','left',...
        'Units','normalized','FontSize',10,'interpreter','latex');
      
%       if noiseflag
% %         print(f.fig,'-dpdf','/home/chaosong/Pictures/sig+map+rm2ndarynoi.pdf');
%         print(f.fig,'-dpdf','/home/chaosong/Pictures/sig+map+rm2ndarynoiv2.pdf');
%       else
% %         print(f.fig,'-dpdf','/home/chaosong/Pictures/sig+map+rm2ndary.pdf');
%         print(f.fig,'-dpdf','/home/chaosong/Pictures/sig+map+rm2ndaryv2.pdf');
%       end

%       keyboard

      %% what is the distance for consecutive sourcces, but in terms of arrival time?
      %note the 'tsep' obtained from the deconvolved positive peaks should be identical to that if
      %obtained from the deconvolved impulses themselves, which represent the arrival indices of the
      %zero-crossing
      ista=1;
      %convert time offset to relative loc
      [imploc, indinput] = off2space002(impindepst(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
      [impindepstst, indsort] = sortrows(impindepst, (ista-1)*2+1);
      implocst = imploc(indsort, :);
      tarvlsplst = impindepstst(:,(ista-1)*2+1);
      
      dtarvlnn1 = diffcustom(tarvlsplst, 1,'forward');
      distnn1 = sqrt(diffcustom(implocst(:,1),1,'forward').^2 + ...
        diffcustom(implocst(:,2),1,'forward').^2 );
      dtarvlnn1all = [dtarvlnn1all; dtarvlnn1];
      distarvlnn1all = [distarvlnn1all; distnn1];

      dtarvlnn2 = diffcustom(tarvlsplst, 2,'forward');
      distnn2 = sqrt(diffcustom(implocst(:,1),2,'forward').^2 + ...
        diffcustom(implocst(:,2),2,'forward').^2 );
      dtarvlnn2all = [dtarvlnn2all; dtarvlnn2];
      distarvlnn2all = [distarvlnn2all; distnn2];

      dtarvlnn3 = diffcustom(tarvlsplst, 3,'forward');
      distnn3 = sqrt(diffcustom(implocst(:,1),3,'forward').^2 + ...
        diffcustom(implocst(:,2),3,'forward').^2 );
      dtarvlnn3all = [dtarvlnn3all; dtarvlnn3];
      distarvlnn3all = [distarvlnn3all; distnn3];
            
      dift = dtarvlnn1;
      dist = distnn1;      
      nsep = 1;
      ttype = 'tarvl';
      %%%Projected distance along specific directions, eg., propagation and its orthogonal, in terms of arrival time
      [f,distprop,distort,stats] = plt_srcprojdist(implocst,nsep,sps,dist,dift,tarvlsplst,ttype);
      if ~isempty(distprop) && ~isempty(distort)
%         close(f.fig);
        distarvlpropall = [distarvlpropall; distprop];
        distarvlortall = [distarvlortall; distort];
        dtarvlpropall = [dtarvlpropall; dift];
        propang(iii) = stats.angrmse;
        proppear(iii) = stats.pearwt;
      end
%       orient(f.fig,'landscape');
%       if noiseflag
%         print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',ttype,'nn1distnoi.pdf'));
%       else
%         print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',ttype,'nn1dist.pdf'));
%       end
% keyboard
      
      end
      
end


