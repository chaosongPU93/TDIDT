% deconv_4s_exp_4thsta_fn.m
% function rststruct = deconv_4s_exp_4thsta_fn(idxbst,normflag,noiseflag,pltflag,rccmwlen)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to 'deconvbursts002_4s_exp', 'deconv_ref_4s_exp_4thsta_fn', this 
% script aims to combine all useful operations based on our experience up
% to now, to revise the deconvolution to each burst window based on the 
% whole window, rather than breaking into 25-s long windows, then concatenate
% them. 
%
% --It is useful because if you see a migration in some bursts using short
% windows, is that due to the artifical movement of alignment based on each
% short window due to the enforcement, or is intrinsic. In addition, the 
% this movement might be different depending on how large the scale in time
% you are looking at. So it is critical to convince others that a migrating
% pattern is even evident using the whole window, which then motivate you 
% to use short windows to image the migration better, as you might lose some
% resolution at the edge if the grouping/detecting window is invariant. 
%
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/12/28
% Last modified date:   2022/12/28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e   % Set the format to 5-digit floating point
clear
clc
% close all

%% for easy testing
defval('idxbst',181);   % the 250-s example 
defval('normflag',0);
defval('noiseflag',0);
defval('pltflag',1);

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
if pltflag
    set(0,'DefaultFigureVisible','on');
else
    set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
end

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
STA = detrend(ccstack);

ccstackort = [];
for ista = 1: nsta
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
        num2str(templensec), 's_', 'BBCCS_', 'ort_Nof_Non_Chao_catnew');
    ccstackort(:,ista) = load(fname);
end
STAort = detrend(ccstackort);

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
% rccmwlen=sps/2;
% rccmwlen=sps;

k = 0;  % counting the burst windows

off1ic = zeros(size(trange,1),nsta);  % single best alignment 'computed' between ALL stas wrt 1 for entire win
off1i = zeros(size(trange,1),nsta);  % single best alignment 'actually used' between ALL stas wrt 1 for entire win
off14pred = zeros(size(trange,1),nsta-3); %empirical pred of off14 from plane fit given single best alignment
ccali = zeros(size(trange,1),1);  % CC value using the best alignment
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
rccsrcall = [];  %mean RCC among trio and RCC14 at src arrival 
rccpairsrcall = []; %RCC for trio sta pairs at src arrival 
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
    % for j = 14: size(rangetemp,1)  
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
        
      %%%obtain a single best alignment based on the entire win 
      %       optcc = optseg(:, 2:end);
      optcc = detrend(optseg(1+msftaddm: end-msftaddm, 2:end));
      msftadd = (round(max(abs([off12ran off13ran])))+1)*sps/40;  %+1 for safety
      ccmid = ceil(size(optcc,1)/2);
      ccwlen = round(size(optcc,1)-2*(msftadd+1));
      ccmin = 0.01;  % depending on the length of trace, cc could be very low
      iup = 1;    % times of upsampling
      [off12con,off13con,ccali(k),iloopoff,loopoff] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
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
        [coef,lag] = xcorr(optcc(:,1), optcc(:,ista), msftadd, 'coeff');
        [mcoef(ista-3, 1), idx] = max(coef);   % max of master raw cc
        off1ic(k,ista) = lag(idx);   % offset in samples
        mlag(ista-3, 1) = lag(idx);
        envrat(ista-3, 1) = median(envelope(optcc(:,ista)))./median(envelope(optcc(:,1)));
        %do an overall CC between 4th and 2nd/3rd stas, to see which one they are most coehrent with
        for jjj = 2:3
          [coef,lag] = xcorr(optcc(:,jjj), optcc(:,ista), msftadd, 'coeff');
          [mcoef(ista-3, jjj), idx] = max(coef);
          mlag(ista-3, jjj) = lag(idx);
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
      flgalign14 = 1;  %do you want to actually align them?
%       flgalign14 = 0;
      if flgalign14
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
                                
%       %%%segment including both the signal and noise together
%       pltdat = [];  % 4-s noise, signal win
%       pltdat(:, 1:2) = STAopt(max(floor((tstbuf-4)*sps+1-overshoot),1): min(floor(tedbuf*sps+overshoot),86400*sps), 1:2); % sta 1
%       pltdat(:, 3) = STAopt(max(floor((tstbuf-4)*sps+1-overshoot)-off1i(k,2),1): ...
%                             min(floor(tedbuf*sps+overshoot)-off1i(k,2),86400*sps), 3); % sta 2
%       pltdat(:, 4) = STAopt(max(floor((tstbuf-4)*sps+1-overshoot)-off1i(k,3),1): ...
%                             min(floor(tedbuf*sps+overshoot)-off1i(k,3),86400*sps), 4); % sta 3                          
%       %compute running CC between 3 stations
%       [ircc,rcc12] = RunningCC(pltdat(:,2), pltdat(:,3), rccmwlen);
%       [~,rcc13] = RunningCC(pltdat(:,2), pltdat(:,4), rccmwlen);
%       [~,rcc23] = RunningCC(pltdat(:,3), pltdat(:,4), rccmwlen);
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
      [ircc,rcc12] = RunningCC(sigsta(:,1), sigsta(:,2), rccmwlen);
      [~,rcc13] = RunningCC(sigsta(:,1), sigsta(:,3), rccmwlen);
      [~,rcc23] = RunningCC(sigsta(:,2), sigsta(:,3), rccmwlen);
      ircc = ircc-overshoot;
      rcc = (rcc12+rcc13+rcc23)/3;
      rccpair = [rcc12 rcc13 rcc23];
      rcc1i = zeros(length(rcc),nsta-3);
      for ista = 4:nsta
        [~,rcc1i(:,ista-3)] = RunningCC(sigsta(:,1), sigsta(:,ista), rccmwlen);
      end
      sigsta = detrend(sigsta(overshoot+1:end-overshoot, :));  %excluding the overshoot
      
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
      [irccort,rcc12] = RunningCC(sigstaort(:,1), sigstaort(:,2), rccmwlen);
      [~,rcc13] = RunningCC(sigstaort(:,1), sigstaort(:,3), rccmwlen);
      [~,rcc23] = RunningCC(sigstaort(:,2), sigstaort(:,3), rccmwlen);
      irccort = irccort-overshoot;
      rccort = (rcc12+rcc13+rcc23)/3;
      sigstaort = detrend(sigstaort(overshoot+1:end-overshoot, :));  %excluding the overshoot

      cc12 = xcorr(sigsta(:,1), sigsta(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
      cc13 = xcorr(sigsta(:,1), sigsta(:,3),0,'normalized');
      cc23 = xcorr(sigsta(:,2), sigsta(:,3),0,'normalized');
      ccpair = [cc12 cc13 cc23];
      mrcc(k,1) = median(rcc);
      mcc(k,1) = (cc12+cc13+cc23)/3;

      %if only use the mean RCC from the 2 pairs that have the highest overall CC
      [~,ind] = min(ccpair);
      rcc = sum(rccpair(:,setdiff(1:3,ind)), 2) / 2;

%       end
      
      %%%plot to show how does the RCC/CC change wrt. the diff alignment between sta pairs 12 and 13
      %%%note that the result is gonna be burst dependent, as data is changing
      %%%in contrast, if using templates for similar analysis, result would be invariant 
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
      
      %%%what if we don't use rcc as the weighting, ie., set to a constant of 1
%       rcc = ones(size(sigsta,1),1);
      
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
        lrcc = length(rcc); % length of running CC
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
          wtcoef = rcc(pkind).* pkhgt;
        else
          wtcoef = rcc1i(pkind, ista-3).* pkhgt;
        end
        fixthresh(ista) = median(wtcoef);  % median of the weighted master CC, could be percentile?
      end
      
      %FLAG to simulate the behavior of noise
%       noiseflag = 1;

%       seedmat = randi(1000,200,1);
%       for iii = 1: length(seedmat)
%       seed = seedmat(iii);

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
%         plot(optseg(:,end),'k','linew',1);
%         xlabel('Samples at 160 sps');
%         ylabel('Amplitude');
%         xlim([0 lsig]);
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
%         subplot(3,1,3)
%         plot(ft,pharand(:,end));
%         xlabel('Frequency (Hz)');
%         ylabel('Randomized phase (rad)');
%                 
        %construct record with the same amplitude but random phase
        xfrand = amp.*nfft.*exp(1i.*pharand);
        optseg(:,2:end) = real(ifft(xfrand,nfft));

%         figure
%         plot(optseg(:,end),'k','linew',1);
%         xlabel('Samples at 160 sps');
%         ylabel('Amplitude');
%         xlim([0 lsig]);
%         ylim([-0.3 0.3]);
%                
        %%%for orthogonal components, with the same phase??
        [xf,ft,amp,pha] = fftspectrum(ortseg(:,2:end), nfft, sps,'twosided');
%         pharand = (rand(nfft,3)-0.5)*2*pi;
        xfrand = amp.*nfft.*exp(1i.*pharand);
        ortseg(:,2:end) = real(ifft(xfrand,nfft));

        ircc = [];   % indices of RCC
        rcc = [];  % average RCC
        rcc1i = [];  % RCC between sta 1 and 4
        rccpair = [];  % RCC between each station pair, order is 12, 13, 23
        ccpair = []; % 0-lag overall cc of each subwin, between each station pair, order is 12, 13, 23 
        cc1i = []; % same as above, but between sta 1 and 4

        %%%obtain a single best alignment based on the entire win 
%       optcc = optseg(:, 2:end);
        optcc = detrend(optseg(1+msftaddm: end-msftaddm, 2:end));
        msftadd = 0.5*sps;
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
          [coef,lag] = xcorr(optcc(:,1), optcc(:,ista), msftadd, 'coeff');
          [mcoef(ista-3, 1), idx] = max(coef);   % max of master raw cc
          off1ic(k,ista) = lag(idx);   % offset in samples
          mlag(ista-3, 1) = lag(idx);
          envrat(ista-3, 1) = median(envelope(optcc(:,ista)))./median(envelope(optcc(:,1)));
          %do an overall CC between 4th and 2nd/3rd stas, to see which one they are most coehrent with
          for jjj = 2:3
            [coef,lag] = xcorr(optcc(:,jjj), optcc(:,ista), msftadd, 'coeff');
            [mcoef(ista-3, jjj), idx] = max(coef);
            mlag(ista-3, jjj) = lag(idx);
            envrat(ista-3, jjj) = median(envelope(optcc(:,ista)))./median(envelope(optcc(:,jjj)));
          end
          %empirical prediction from plane fitting in 'empioffset4thsta002'
          off14pred(k,ista-3) = round(off14mod(ista-3,1).*off1ic(k,2) + off14mod(ista-3,2).*off1ic(k,3) + ...
            off14mod(ista-3,3));
        end
        %       off1ia(k,4:end) = off14pred(k,:);
        
        %for synthetic noise case, do NOT align traces
        off1i(k,1:3) = [0 0 0];  % technically unnecessary given it has been predefined 
        
        %for synthetic noise case, do NOT align traces
        %Choice to make upon the actual-used alignment at 4th stations
%         flgalign14 = 1;  %do you want to actually align them?
        flgalign14 = 0;
        if flgalign14
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
          sigsta(:,ista) = tmp;
        end
        %compute running CC between 3 stations
        [ircc,rcc12] = RunningCC(sigsta(:,1), sigsta(:,2), rccmwlen);
        [~,rcc13] = RunningCC(sigsta(:,1), sigsta(:,3), rccmwlen);
        [~,rcc23] = RunningCC(sigsta(:,2), sigsta(:,3), rccmwlen);
        ircc = ircc-overshoot;
        rcc = (rcc12+rcc13+rcc23)/3;
        rccpair = [rcc12 rcc13 rcc23];
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
        [irccort,rcc12] = RunningCC(sigstaort(:,1), sigstaort(:,2), rccmwlen);
        [~,rcc13] = RunningCC(sigstaort(:,1), sigstaort(:,3), rccmwlen);
        [~,rcc23] = RunningCC(sigstaort(:,2), sigstaort(:,3), rccmwlen);
        irccort = irccort-overshoot;
        rccort = (rcc12+rcc13+rcc23)/3;
        sigstaort = detrend(sigstaort(overshoot+1:end-overshoot, :));  %excluding the overshoot

        cc12 = xcorr(sigsta(:,1), sigsta(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
        cc13 = xcorr(sigsta(:,1), sigsta(:,3),0,'normalized');
        cc23 = xcorr(sigsta(:,2), sigsta(:,3),0,'normalized');
        ccpair = [cc12 cc13 cc23];
        mrcc(k,1) = median(rcc);
        mcc(k,1) = (cc12+cc13+cc23)/3;
  
        %if only use the mean RCC from the 2 pairs that have the highest overall CC
        [~,ind] = min(ccpair);
        rcc = sum(rccpair(:,setdiff(1:3,ind)), 2) / 2;
          
      end  

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
          %%%As of 2022/09/27, 'threshold' would be auto-computed based on data then feed to the
          %%%decon if the 'noiseflag' is on
          [sigdecon(:,ista),pred(:,ista),res,dresit,mfitit,ampit{ista},nit,fighdl] = ...
            iterdecon(sig,wlet,rcc,noi,fixthresh(ista),dt,twlet,width,dres_min,...
            mfit_min,nit_max,nimp_max,fpltit,fpltend,fpltchk);
        else
          [sigdecon(:,ista),pred(:,ista),res,dresit,mfitit,ampit{ista},nit,fighdl] = ...
            iterdecon(sig,wlet,rcc,noi,[],dt,twlet,width,dres_min,mfit_min,nit_max,nimp_max,...
            fpltit,fpltend,fpltchk);
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
      
    %   nsrc(iii) = size(impindep,1);
    %   end
    
    % figure
    % histogram(nsrc);
    % text(0.05,0.9,sprintf('med=%d',median(nsrc)),'Units','normalized');
    % xlabel('number of sources');
    % ylabel('counts');
    
    % keyboard

      %%%plot the individually deconvolved impulses and the grouping result                 
%       [f] = plt_groupimptriplets(sigdecon,impindep,stas,ircc,rcc);
     
      %note here 'impindepst' inherits the first 6 cols from 'impindep', but the last three cols 
      %are adjusted from arrival time difference to the true location offset accounting for the best
      %alignment upon the entire window that is also used in grouping!
      impindepst(:,7:8) = impindepst(:,7:8)+repmat([off1i(k,2) off1i(k,3)],size(impindepst,1),1); %account for prealignment
      impindepst = sortrows(impindep,1);
      
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
        sps,50,'mean','tori');
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
      ppkindep = impindep;
      for ista = 1: nsta
        ppkindep(:,(ista-1)*2+1) = ppkindep(:,(ista-1)*2+1)+ppeaks(ista)-zcrosses(ista);
      end
      npkindep = impindep;  %negative peaks 
      for ista = 1: 3
        npkindep(:,(ista-1)*2+1) = npkindep(:,(ista-1)*2+1)+npeaks(ista)-zcrosses(ista);
      end

      %use function 'removesecondarysrc' to remove the smaller amplitude, secondary triplets that 
      %are too close in time to large triplets
      minsrcsep = 20;   % min allowed separation between deconvolved peaks
      [ppkindepsave,indremove] = removesecondarysrc(ppkindep,sigsta(:,1:3));
%       [pkindepsave,indremove] = removesecondarysrc(pkindep,sigsta,minsrcsep);
      
      %REMOVE the secondary sources from the grouped result
      impindep(indremove, :) = [];
      ppkindep(indremove, :) = [];
      npkindep(indremove, :) = [];
      impindepst = sortrows(impindep,1);
      nsrcraw(iii,1) = size(impindepst,1);  % number of sources AFTER removing 2ndary 

      if ~isempty(impindepst)

      %% plot the scatter of sources in terms of offsets, accounting for prealignment offset
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
        sps,50,'mean','tori');
      scatter(ax1,off1i(k,2),off1i(k,3),20,'ks','filled','MarkerEdgeColor','k');
      title(ax1,'Independent, grouped, no secondary sources');

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

      %% what is the distance for consecutive soures, in terms of origin time?      
      %convert time offset to relative loc
      [imploc, indinput] = off2space002(impindepst(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
      [torisplst, indsort] = sortrows(torispl,1);
      dtorinn1 = diff(torisplst);
      implocst = imploc(indsort, :);
      impindepstst = impindepst(indsort, :);
            
      %between Nth and (N-1)th source
      distnn1 = sqrt(diffcustom(implocst(:,1),1,'forward').^2 + ...
        diffcustom(implocst(:,2),1,'forward').^2 );
      dtorinn1all = [dtorinn1all; dtorinn1];
      distorinn1all = [distorinn1all; distnn1];
      
      %between Nth and (N-2)th source
      dtorinn2 = diffcustom(torisplst,2,'forward');
      distnn2 = sqrt(diffcustom(implocst(:,1),2,'forward').^2 + ...
        diffcustom(implocst(:,2),2,'forward').^2 );
      dtorinn2all = [dtorinn2all; dtorinn2];
      distorinn2all = [distorinn2all; distnn2];
      
      %between Nth and (N-3)th source
      dtorinn3 = diffcustom(torisplst,3,'forward');
      distnn3 = sqrt(diffcustom(implocst(:,1),3,'forward').^2 + ...
        diffcustom(implocst(:,2),3,'forward').^2 );
      dtorinn3all = [dtorinn3all; dtorinn3];
      distorinn3all = [distorinn3all; distnn3];
      
      dist = distnn1;
      dift = dtorinn1;
      nsep = 1;
      %%%Projected distance along specific directions, eg., propagation and its orthogonal, in terms of origin time
      ttype = 'tori';
      [f,distprop,distort] = plt_srcprojdist(implocst,nsep,sps,dist,dift,torisplst,ttype);
      if ~isempty(distprop) && ~isempty(distort)
        close(f.fig);
        distoripropall = [distoripropall; distprop];
        distoriortall = [distoriortall; distort];
        dtoripropall = [dtoripropall; dift];
      end

            %% what is the distance for consecutive sourcces, but in terms of arrival time?
      %note the 'tsep' obtained from the deconvolved positive peaks should be identical to that if
      %obtained from the deconvolved impulses themselves, which represent the arrival indices of the
      %zero-crossing
      ista=1;
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
      [f,distprop,distort] = plt_srcprojdist(implocst,nsep,sps,dist,dift,tarvlsplst,ttype);
      if ~isempty(distprop) && ~isempty(distort)
        close(f.fig);
        distarvlpropall = [distarvlpropall; distprop];
        distarvlortall = [distarvlortall; distort];
        dtarvlpropall = [dtarvlpropall; dift];
      end

    
      %% plot the scatter of sources in terms of rela locations
      xran = [-3 3];
      yran = [-3 3];
      cran = [0 lsig/sps];
      f2.fig = figure;
      f2.fig.Renderer = 'painters';
      ax2=gca;
      [ax2] = plt_decon_imp_scatter_space(ax2,impindepst,xran,yran,cran,offxran,...
        offyran,sps,50,ftrans,'mean','tori');
      plot(ax2,xcut,ycut,'k-','linew',2);
      title(ax2,'Independent, grouped, no secondary sources');
      

      %% test a bunch of offset max to see residual reduction VS. # of sources left
%       %%%According to the resulting plot, ~1.6* RMSE (12 samples at 160 sps) seems proper
%       modname = 'timeoff_plfit_4thsta_160sps.mat';
%       planefit = load(strcat(rstpath, '/MAPS/',modname));
%       ista = 7;
%       wlet = greenf(:,ista);  %template here is best aligned, tapered, linear trend removed, filtered
%       lwlet = length(wlet);
%       sig = sigsta(:,ista); %best aligned, filtered, tapered
%       lsig = length(sig);
%       
%       dt = 1/sps;  % sampling interval
%       twlet = zcrosses(ista)*dt;
%       fpltit = 0;  % plot flag for each iteration
%       fpltend = 0;  % plot flag for the final iteration
%       fpltchk = 0; % plot flag for intermediate computations
%       rmse = planefit.gof{ista-3}.rmse;
%       %+-1-sigma ~0.68, 1.5-sigma ~0.87, 2-sigma ~0.95, 3-sigma ~0.997
%       %         sigma1 = round(rmse);
%       %         sigma2 = round(2*rmse);
%       offmaxtry = 4:1:40;
%       for ioff = 1: length(offmaxtry)
%         offmax = offmaxtry(ioff);
% %         [~,pred,res,dresit,mfitit,ampiti,fighdl] = iterdecon_4thsta(sig,wlet,...
% %           irccran,rcc1i(:,ista-3),dt,twlet,impindep,stas(ista,:),off1i(k,ista),[],offmax,...
% %           fpltit,fpltend,fpltchk);
%         [~,pred,res,dresit,mfitit,ampiti,fighdl] = iterdecon_4thsta(sig,wlet,...
%           irccran,rcc1icat(:,ista-3),dt,twlet,impindep,stas(ista,:),off1i(k,ista),off1iw(:,ista),offmax,...
%           fpltit,fpltend,fpltchk);
%         imptemp = impindep;
%         imptemp(:,9+(ista-4)*2+1) = ampiti(:,1);
%         imptemp(:,9+(ista-3)*2) = ampiti(:,2);
%         pred14off(:,ista-3) = ampiti(:,end);  %difference between found peak and predicted arrival
% 
%               %%% further eliminate sources that fail the check at 4th stations
%         trust4th = 7; % trust KLNB the most among all 4th stations
%         indremove = find(imptemp(:,9+(trust4th-4)*2+1)==0 & imptemp(:,9+(trust4th-3)*2)==0);
%         pred14offtr = pred14off(setdiff(1:size(pred14off,1),indremove),trust4th-3);
%         imptemp(indremove,:) = [];
%         imptempst = sortrows(imptemp,1);
%         nsrc(k,ioff) = size(imptempst,1);
%         nsrcnorm(k,ioff) = size(imptempst,1)/ size(impindep,1)* 100;  %normalized by # after removing 2ndary srcs
% 
%         %%%final prediction via convolution between grouped impulses and template at each station 
%         [~,predgrp,resgrp,predgrpl,resgrpl,resred]=predsig_conv_imptemp(sigsta,optdat,imptempst,...
%           greenf,zcrosses,overshoot,stas,0);
%         rr(k,ioff,1) = resred(trust4th,3);
%         rr(k,ioff,2) = mean(resred(1:3,3));
%       end
%       
%       %%
%       figure
%       ax = subplot(1,2,1);
%       hold(ax,'on');
%       ax.Box='on'; grid(ax,'on');
%       yyaxis(ax,'left');
%       plot(ax,offmaxtry/rmse,nsrc(k,:),'b-','linew',1);
%       p1=plot(ax,[mswccpksep(ista) mswccpksep(ista)]/2/rmse,ax.YLim,'k--');
%       xlabel(ax,sprintf('Allowed diff from tarvl pred (* RMSE %.1f spls)',rmse));
%       ylabel(ax,'# of sources left');
%       yyaxis(ax,'right');
%       plot(ax,offmaxtry/rmse,rr(k,:,1),'r-','linew',1);
%       p2=plot(ax,offmaxtry/rmse,rr(k,:,2),'k:','linew',1);
%       xticks(ax,0:0.5:6);
%       ylabel(ax,'Change in L-2 norm (%)');
%       legend(ax,[p1,p2],'Half of med sep of peaks in sig-wlet CC','Mean MR at trio stas',...
%         'Location','south');
%       ylim(ax,[16 34]);
% 
%       ax = subplot(1,2,2);
%       hold(ax,'on');
%       ax.Box='on'; grid(ax,'on');
%       scatter(ax,nsrcnorm(k,:),rr(k,:,1),25,'k','filled');
%       xlabel('# sources left (%)');
%       ylabel('Change in L-2 norm (%)');
%       ylim(ax,[16 34]);
%       
% 
%       keyboard

      %%%carry out 'deconvolution' at 4th stations as well for the tarvl and amp
      modname = 'timeoff_plfit_4thsta_160sps.mat';
      planefit = load(strcat(rstpath, '/MAPS/',modname));
%       rmse = planefit.gof.rmse;

      pred4off = [];
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
          [sigdecon(:,ista),pred,res,dresit,mfitit,ampit{ista},fighdl] = ...
            iterdecon_4thsta(sig,wlet,irccran,rcc1icat(:,ista-3),fixthresh(ista),...
            dt,twlet,impindep,stas(ista,:),off1i(k,ista),off1iw(:,ista),offmax,...
            fpltit,fpltend,fpltchk);
        else
%         [sigdecon(:,ista),pred,res,dresit,mfitit,ampit{ista},fighdl] = iterdecon_4thsta(sig,wlet,...
%           irccran,rcc1i(:,ista-3),dt,twlet,impindep,stas(ista,:),off1i(k,ista),[],offmax,...
%           fpltit,fpltend,fpltchk);
          [sigdecon(:,ista),pred,res,dresit,mfitit,ampit{ista},fighdl] = ...
            iterdecon_4thsta(sig,wlet,irccran,rcc1icat(:,ista-3),[],...
            dt,twlet,impindep,stas(ista,:),off1i(k,ista),off1iw(:,ista),offmax,...
            fpltit,fpltend,fpltchk);
        end
        
        ampiti = ampit{ista};
        impindep(:,9+(ista-4)*2+1) = ampiti(:,1);
        impindep(:,9+(ista-3)*2) = ampiti(:,2);
        if ista == nsta
          pred4off(:,ista-3) = ampiti(:,end);  %difference between found peak and predicted arrival
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
      
      %%% further ELIMINATE sources that fail the check at 4th stations
      trust4th = 7; % trust KLNB the most among all 4th stations
      indremove = find(impindep(:,9+(trust4th-4)*2+1)==0 & impindep(:,9+(trust4th-3)*2)==0);
      pred4offtr = pred4off(setdiff(1:size(pred4off,1),indremove),trust4th-3);
      impindep(indremove,:) = [];
      ppkindep(indremove, :) = [];
      npkindep(indremove, :) = [];
      impindepst = sortrows(impindep,1);
      impindepstall = [impindepstall; impindepst];
      pred4offtrall = [pred4offtrall; pred4offtr];

      if ~isempty(impindepst)
      %%% plot the scatter of sources in terms of offsets, accounting for prealignment offset
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
      close(f1.fig);
%       title(ax1,'Independent, grouped, no secondary sources, after 4th-sta check');
%       print(f1.fig,'-dpdf','/home/chaosong/Pictures/checkat4th.pdf');
%       print(f1.fig,'-dpdf','/home/chaosong/Pictures/checkat4thnoi.pdf');

      %%%final prediction via convolution between grouped impulses and template at each station 
      [f2,predgrp,resgrp,predgrpl,resgrpl,l2normred(iii,:,:)]=predsig_conv_imptemp(sigsta,optdat,impindepst,...
        greenf,zcrosses,overshoot,stas,1);
      close(f2.fig);  

      %% recompute time separation and distance after 4th sta check  
      %%%Plot the separation in time between these preserved positive peaks after removing the
      %%secondary ones, to see if they can be too close to each other
      %discard the sources that are determined to be too close and secondary compared to a major source
      ppkindepsave = ppkindep;
      if size(ppkindepsave,1) > 1
        [f,tsep,pkist,indpk,indreverse] = plt_tsep_deconpk(ppkindepsave,sps);
        tsep4thall = [tsep4thall; tsep];
        close(f.fig);
      end
%       median(tsep)

      %%%distance in terms of arrival time 
      ista=1;
      [impindepstst, indsort] = sortrows(impindepst, (ista-1)*2+1);
      implocst = imploc(indsort, :);
      tarvlsplst = impindepstst(:,(ista-1)*2+1);
            
      dtarvlnn1 = diffcustom(tarvlsplst, 1,'forward');
      distnn1 = sqrt(diffcustom(implocst(:,1),1,'forward').^2 + ...
        diffcustom(implocst(:,2),1,'forward').^2 );
      dtarvlnn14thall = [dtarvlnn14thall; dtarvlnn1];
      distarvlnn14thall = [distarvlnn14thall; distnn1];

      dtarvlnn2 = diffcustom(tarvlsplst, 2,'forward');
      distnn2 = sqrt(diffcustom(implocst(:,1),2,'forward').^2 + ...
        diffcustom(implocst(:,2),2,'forward').^2 );
      dtarvlnn24thall = [dtarvlnn24thall; dtarvlnn2];
      distarvlnn24thall = [distarvlnn24thall; distnn2];

      dtarvlnn3 = diffcustom(tarvlsplst, 3,'forward');
      distnn3 = sqrt(diffcustom(implocst(:,1),3,'forward').^2 + ...
        diffcustom(implocst(:,2),3,'forward').^2 );
      dtarvlnn34thall = [dtarvlnn34thall; dtarvlnn3];
      distarvlnn34thall = [distarvlnn34thall; distnn3];

      dift = dtarvlnn1;
      dist = distnn1;
      nsep = 1;      
      ttype = 'tarvl';
      %%%Projected distance along specific directions, eg., propagation and its orthogonal, in terms of arrival time
      [f,distprop,distort] = plt_srcprojdist(implocst,nsep,sps,dist,dift,tarvlsplst,ttype);
      if ~isempty(distprop) && ~isempty(distort)
        close(f.fig);
        distarvlprop4thall = [distarvlprop4thall; distprop];
        distarvlort4thall = [distarvlort4thall; distort];
        dtarvlprop4thall = [dtarvlprop4thall; dift];
      end

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
      
      %% plot amp ratio 12 and 13, and 23, without 2ndary sources
%       %%ideally, when templates and data are not normalized, and there is no particular noise or any
%       %%other factors causing the amplitude scaling between temp and data and each station to be
%       %%vastly different, then for each deconvolved source, the direct impulse amp should be
%       %%~identical at all stations, ie., the ratio between station pairs should be ~1
if isempty(impindepst)
  continue
else
srcampr = [impindepst(:,2)./impindepst(:,4) impindepst(:,2)./impindepst(:,6) ...
           impindepst(:,4)./impindepst(:,6) impindepst(:,2)./impindepst(:,17)];
psrcamprs = [impindepst(:,2)*max(greenf(:,1))./(impindepst(:,4)*max(greenf(:,2))) ...
             impindepst(:,2)*max(greenf(:,1))./(impindepst(:,6)*max(greenf(:,3))) ...
             impindepst(:,4)*max(greenf(:,2))./(impindepst(:,6)*max(greenf(:,3))) ...
             impindepst(:,2)*max(greenf(:,1))./(impindepst(:,17)*max(greenf(:,7)))];
psrcamps = [impindepst(:,2)*max(greenf(:,1)) impindepst(:,4)*max(greenf(:,2)) ...
            impindepst(:,6)*max(greenf(:,3)) impindepst(:,17)*max(greenf(:,7))];

nsrcamprs = [impindepst(:,2)*min(greenf(:,1))./(impindepst(:,4)*min(greenf(:,2))) ...
             impindepst(:,2)*min(greenf(:,1))./(impindepst(:,6)*min(greenf(:,3))) ...
             impindepst(:,4)*min(greenf(:,2))./(impindepst(:,6)*min(greenf(:,3))) ...
             impindepst(:,2)*min(greenf(:,1))./(impindepst(:,17)*min(greenf(:,7)))];

nsrcamps = [impindepst(:,2)*min(greenf(:,1)) impindepst(:,4)*min(greenf(:,2))...
            impindepst(:,6)*min(greenf(:,3)) impindepst(:,17)*min(greenf(:,7))];

msrcampr(iii,:) = median(srcampr, 1);
madsrcampr(iii,:) = mad(srcampr, 1, 1);
mpsrcamprs(iii,:) = median(psrcamprs, 1);
madpsrcamprs(iii,:) = mad(psrcamprs, 1, 1);
mnsrcamprs(iii,:) = median(nsrcamprs, 1);
madnsrcamprs(iii,:) = mad(nsrcamprs, 1, 1);
nsrc(iii,1) = size(srcampr,1);
srcamprall = [srcamprall; srcampr];
psrcampsall = [psrcampsall; psrcamps];
nsrcampsall = [nsrcampsall; nsrcamps];
psrcamprsall = [psrcamprsall; psrcamprs];
nsrcamprsall = [nsrcamprsall; nsrcamprs];

%what is the deviation of amp ratio from the median for each source?
lndevsrcampr = srcampr-median(srcampr, 1); % in linear scale
lgdevsrcampr = log10(srcampr)-log10(median(srcampr, 1)); % in log scale, note that log2+log5=log10, so this means a ratio
lndevsrcamprall = [lndevsrcamprall; lndevsrcampr];
lgdevsrcamprall = [lgdevsrcamprall; lgdevsrcampr];

%%%what are the corresponding RCC at each source
rccpairsrc = [];
rccpairsrc(:,1) = rccpair(round(mean(impindepst(:,[1 3]),2)),1);
rccpairsrc(:,2) = rccpair(round(mean(impindepst(:,[1 5]),2)),2);
rccpairsrc(:,3) = rccpair(round(mean(impindepst(:,[3 5]),2)),3);
rccpairsrcall = [rccpairsrcall; rccpairsrc];

%use the concatenated rcc at the average arrival time of each source
rccsrc = [];
rccsrc(:,1) = rcc(round(mean(impindepst(:,[1 3 5]),2)));
rccsrc(:,2) = rcc1i(impindepst(:,9+(trust4th-4)*2+1),trust4th-3);
rccsrcall = [rccsrcall; rccsrc];
  
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
%         plot([log10(msrcampr(iii,ista)) log10(msrcampr(iii,ista))],ax.YLim,'r--');      
% %         errorbar(mclppkhtwfr(i),mpsrcamprs(i),madpsrcamprs(i),madpsrcamprs(i),...
% %           madclppkhtwfr(i),madclppkhtwfr(i),'color',[.5 .5 .5],'linewidth',0.8,'CapSize',5);
%         text(0.05,0.9,sprintf('med=%.2f',msrcampr(iii,ista)),'Units','normalized');
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
%         plot([log10(mpsrcamprs(iii,ista)) log10(mpsrcamprs(iii,ista))],ax.YLim,'r--');      
%         text(0.05,0.9,sprintf('med=%.2f',mpsrcamprs(iii,ista)),'Units','normalized');
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
%         plot([log10(mnsrcamprs(iii,ista)) log10(mnsrcamprs(iii,ista))],ax.YLim,'r--');      
%         text(0.05,0.9,sprintf('med=%.2f',mnsrcamprs(iii,ista)),'Units','normalized');
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
%       mclppkhtwf(iii, :) = median(clppkhtwf, 1);
%       madclppkhtwf(iii, :) = mad(clppkhtwf,1, 1);
clppkhtwfall = [clppkhtwfall; clppkhtwf];
clppkwf = clppk.clppkwf;  %also store the waveform peak separation
clppkwfsep = diff([clppkwf(:,1:3) clppkwf(:,end)]);  %I only want KLNB
clppkwfsepall = [clppkwfsepall; clppkwfsep];

clnpkhtwf = clnpk.clnpkhtwf;
clnpkhtwf = [clnpkhtwf(:,1:3) clnpkhtwf(:,end)];  %I only want KLNB
%       mclnpkhtwf(iii, :) = median(clnpkhtwf, 1);
%       madclnpkhtwf(iii, :) = mad(clnpkhtwf,1, 1);
clnpkhtwfall = [clnpkhtwfall; clnpkhtwf];
clnpkwf = clnpk.clnpkwf;
clnpkwfsep = [clnpkwf(:,1:3) clnpkwf(:,end)];  %I only want KLNB
clnpkwfsepall = [clnpkwfsepall; clnpkwfsep];

%also get the waveform peak separation for the full trace, since the number of peaks is not
%the same for each station, keep 
for ista = [1 2 3 nsta]
  [~, pk] = findpeaks(sigsta(:,ista));
  ppkwfsepmed(iii,ista) = median(diff(pk));
  ppkwfsepmod(iii,ista) = mode(diff(pk));
  [~, pk] = findpeaks(-sigsta(:,ista));
  npkwfsepmed(iii,ista) = median(diff(pk));
  npkwfsepmod(iii,ista) = mode(diff(pk));
end

%       keyboard
end


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
      [predgrp,resgrp,predgrpl,resgrpl]=predsig_conv_imptemp(sigsta,optdat,impindepst,...
        greenf,zcrosses,overshoot,stas,1);

      %% plot the comparision between the real data and prediction using the grouped impulses
      %compute running CC between 3 stations
      [irccp,rccp12] = RunningCC(predgrp(:,1), predgrp(:,2), rccmwlen);
      [~,rccp13] = RunningCC(predgrp(:,1), predgrp(:,3), rccmwlen);
      [~,rccp23] = RunningCC(predgrp(:,2), predgrp(:,3), rccmwlen);
      irccp = irccp-overshoot;
      rccp = (rccp12+rccp13+rccp23)/3;  % use averge
      
      [f] = plt_sigpred(sigsta,predgrp,impindepst,sps,tmaxi,tmaxo,tbosti,tbosto,...
        ircc,rcc,irccp,rccp,rccmwlen,tstbuf,tedbuf,dy,mo,yr);
      

%       keyboard
    % end
    
  % end
  
end

%% Ouput everything in the form of a structure array
rststruct.srcamprall = srcamprall;
rststruct.lndevsrcamprall = lndevsrcamprall;
rststruct.lgdevsrcamprall = lgdevsrcamprall;
rststruct.rcccatsrcall = rccsrcall;
rststruct.rccpairsrcall = rccpairsrcall;
rststruct.psrcampsall = psrcampsall;
rststruct.nsrcampsall = nsrcampsall;
rststruct.psrcamprsall = psrcamprsall;
rststruct.nsrcamprsall = nsrcamprsall;
rststruct.clppkhtwfall = clppkhtwfall;
rststruct.clnpkhtwfall = clnpkhtwfall;
rststruct.clppkhtwfall = clppkwfsepall;
rststruct.clnpkhtwfall = clnpkwfsepall;
rststruct.ppkwfsepmed = ppkwfsepmed;
rststruct.ppkwfsepmod = ppkwfsepmod;
rststruct.npkwfsepmed = npkwfsepmed;
rststruct.npkwfsepmod = npkwfsepmod;


rststruct.pred4offtrall = pred4offtrall;
rststruct.impindepstall = impindepstall;

rststruct.tsepall = tsepall;
rststruct.dtorinn1all = dtorinn1all;
rststruct.distorinn1all = distorinn1all;
rststruct.dtorinn2all = dtorinn2all;
rststruct.distorinn2all = distorinn2all;
rststruct.dtorinn3all = dtorinn3all;
rststruct.distorinn3all = distorinn3all;
rststruct.dtoripropall = dtoripropall;
rststruct.distoripropall = distoripropall;
rststruct.distoriortall = distoriortall;
rststruct.dtarvlnn1all = dtarvlnn1all;
rststruct.distarvlnn1all = distarvlnn1all;
rststruct.dtarvlnn2all = dtarvlnn2all;
rststruct.distarvlnn2all = distarvlnn2all;
rststruct.dtarvlnn3all = dtarvlnn3all;
rststruct.distarvlnn3all = distarvlnn3all;
rststruct.dtarvlpropall = dtarvlpropall;
rststruct.distarvlpropall = distarvlpropall;
rststruct.distarvlortall = distarvlortall;

rststruct.tsep4thall = tsep4thall; 
rststruct.dtorinn14thall = dtorinn14thall;
rststruct.distorinn14thall = distorinn14thall;
rststruct.dtorinn24thall = dtorinn24thall;
rststruct.distorinn24thall = distorinn24thall;
rststruct.dtorinn34thall = dtorinn34thall;
rststruct.distorinn34thall = distorinn34thall;
rststruct.dtoriprop4thall = dtoriprop4thall;
rststruct.distoriprop4thall = distoriprop4thall;
rststruct.distoriort4thall = distoriort4thall;
rststruct.dtarvlnn14thall = dtarvlnn14thall;
rststruct.distarvlnn14thall = distarvlnn14thall;
rststruct.dtarvlnn24thall = dtarvlnn24thall;
rststruct.distarvlnn24thall = distarvlnn24thall;
rststruct.dtarvlnn34thall = dtarvlnn34thall;
rststruct.distarvlnn34thall = distarvlnn34thall;
rststruct.dtarvlprop4thall = dtarvlprop4thall;
rststruct.distarvlprop4thall = distarvlprop4thall;
rststruct.distarvlort4thall = distarvlort4thall;

if ~isempty(impindepst)
  rststruct.nsrcraw = nsrcraw;
  rststruct.nsrc = nsrc;
  rststruct.msrcampr = msrcampr;
  rststruct.madsrcampr = madsrcampr;
  rststruct.mpsrcamprs = mpsrcamprs;
  rststruct.madpsrcamprs = madpsrcamprs;
  rststruct.mnsrcamprs = mnsrcamprs;
  rststruct.madnsrcamprs = madnsrcamprs;
  rststruct.l2normred = l2normred;
end

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
  f4 = initfig(5,5,1,size(pred4offtrall,2)); %initialize fig
  plt_errorof4thtarvlpred(f4,pred4offtrall,offmax,'k');

  %%%preserved sources' amp ratio between 4th and 1st stas
  f5 = initfig(12,5,1,3); %initialize fig
  plt_deconpk_rat14(f5,impindepstall,srcamprall,'k');

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






