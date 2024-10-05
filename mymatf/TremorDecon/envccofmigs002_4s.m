% envccofmigs002_4s.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From the plots of envelopes at different stations from 'envofmigs002_4s',
% we could visionly see how similar the envelopes are. This scripts aims
% to carry out a cross-correlation between the envelopes around the auto-
% detected migrations, and those windows for deconvolution. Comparing the
% CC with that of the PGC-SSIB-SILB, we will have a quantative sense of
% the similarity of seismograms, and possibly know if any station could
% function as the 4th station to check the deconvolution result from the
% trio stations.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/07/26
% Last modified date:   2022/07/26
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

% median(tranmig(:,3)-tranmig(:,2))

%% load HF decon results for the alignment of 25-s windows
% savefile = 'deconv_stats4th_allbstsig.mat';
savefile = 'deconv_stats4th_no23_allbstsig.mat';
load(strcat(rstpath, '/MAPS/',savefile));

% off1i = allbstsig.off1i;
% off1i(1:3,:)


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

idxburst = (1:size(tranbst,1))';

nibst = size(idxburst,1);
[~,~,~,idate,ibst] = indofburst(tranbst,idxburst);

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

% prct = [10 50 90];
prct = [5 50 95];
envprct = zeros(nibst, 3, nsta);  % percentiles of the envelope
envortprct = zeros(nibst, 3, nsta);  % percentiles of the envelope, ort comp.
envvertprct = zeros(nibst, 3, nsta);  % percentiles of the envelope, vert comp.

%%%Flag to indicate if it is necessary to recalculate everything
flagrecalc = 0;
% flagrecalc = 1;

if flagrecalc
  for iii = 1: length(idxburst)
    iii
    %   for i = 1: length(dates)  % dates in each ets
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
    
    if fileflag == 0    % means there are missing files
      fprintf('Day %s / %s will be omitted because of missing files. \n', yr, JDAY);
      continue    % continue to the next day
    end
    
    %read vertical components too, as we want to have a sense of the SNR at the orthogonal and vertical
    %components
    [STAvert,~,~,fileflag] = rd_daily_bpdata(year,jday,prename,stas,PERMSTA,POLSTA,...
      PERMROTS,POLROTS,sps,losig,hisig,npo,npa,[],[],[],[],'Z');
    
    %   %%%CC between stations for migration windows
    %   for j = 1: size(mig,1)
    %     tst = mig(j,2)*3600;  %each migration is long enough, median is ~7 mins, no extra time needed
    %     ted = mig(j,3)*3600;
    %     envseg = envup(max(floor(tst*sps)+1, 1): min(floor(ted*sps),86400*sps), :);
    %
    %     maxlag = 2*sps;
    %
    %     for jj = 4: nsta
    %       [coef, lag] = xcorr(envseg(:,jj), envseg(:,1), maxlag, 'coeff');
    %       [ccmij(1,j+n-size(mig,1),jj-3),mind] = max(coef);
    %       lagmij(1,j+n-size(mig,1),jj-3) = lag(mind);
    %
    %       [coef, lag] = xcorr(envseg(:,jj), envseg(:,2), maxlag, 'coeff');
    %       [ccmij(2,j+n-size(mig,1),jj-3),mind] = max(coef);
    %       lagmij(2,j+n-size(mig,1),jj-3) = lag(mind);
    %
    %       [coef, lag] = xcorr(envseg(:,jj), envseg(:,3), maxlag, 'coeff');
    %       [ccmij(3,j+n-size(mig,1),jj-3),mind] = max(coef);
    %       lagmij(3,j+n-size(mig,1),jj-3) = lag(mind);
    %     end
    %
    %     [coef, lag] = xcorr(envseg(:,1), envseg(:,2), maxlag, 'coeff');
    %     [ccm123(j+n-size(mig,1),1),mind] = max(coef);
    %     lagm123(j+n-size(mig,1),1) = lag(mind);
    %     [coef, lag] = xcorr(envseg(:,1), envseg(:,3), maxlag, 'coeff');
    %     [ccm123(j+n-size(mig,1),2),mind] = max(coef);
    %     lagm123(j+n-size(mig,1),2) = lag(mind);
    %     [coef, lag] = xcorr(envseg(:,2), envseg(:,3), maxlag, 'coeff');
    %     [ccm123(j+n-size(mig,1),3),mind] = max(coef);
    %     lagm123(j+n-size(mig,1),3) = lag(mind);
    %
    %   end
    
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
    tstbuf = min(tcnti(indtmaxi)-2);
    tedbuf = max(tcnti(indtmaxi)+2);
    tlenbuf = tedbuf - tstbuf;
%     %add a buffer of the same length at the start and end
%     tstbuf = tstbuf - tlenbuf;
%     tedbuf = tedbuf + tlenbuf;
    
    %max allowable shift in best alignment
    msftaddm = 1.5*sps+1;  %+1 for safety

    %have some overshoot, so that the resulted rcc would have the same length as the signal
    rccmwsec=0.5;
    rccmwlen=rccmwsec*sps;
    overshoot = rccmwlen/2;

    %%%%2022/09/26, obtain all information for data before decon, so that you know the threshold
    %%%%being used, although this was used mainly for noise experiment, we still use it here
    %chop a record segment
    optseg = STAopt(max(floor(tstbuf*sps+1-overshoot-msftaddm),1): ...
      min(floor(tedbuf*sps+overshoot+msftaddm),86400*sps), :); % sta 1
    ortseg = STAort(max(floor(tstbuf*sps+1-overshoot-msftaddm),1): ...
      min(floor(tedbuf*sps+overshoot+msftaddm),86400*sps), :);
    vertseg = STAvert(max(floor(tstbuf*sps+1-overshoot-msftaddm),1): ...
      min(floor(tedbuf*sps+overshoot+msftaddm),86400*sps), :);

    %save room for memory
    clear STAopt STAort STAvert

    %%%obtain a single best alignment based on the entire win
    %       optcc = optseg(:, 2:end);
    optcc = detrend(optseg(1+msftaddm: end-msftaddm, 2:end));
    msftadd = 10/40*sps;
    loffmax = 4*sps/40;
    ccmid = ceil(size(optcc,1)/2);
    ccwlen = round(size(optcc,1)-2*(msftadd+1));
    ccmin = 0.01;  % depending on the length of trace, cc could be very low
    iup = 1;    % times of upsampling
    [off12con,off13con,ccali(icount),iloopoff,loopoff] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
      ccwlen,msftadd,loffmax,ccmin,iup);
    % if a better alignment cannot be achieved, use 0,0
    if off12con == msftadd+1 && off13con == msftadd+1
      off12con = 0;
      off13con = 0;
      fprintf('Tremor burst %d cannot be properly aligned, double-check needed \n',icount);
    end
    off1i(icount,1) = 0;
    off1i(icount,2) = round(off12con);
    off1i(icount,3) = round(off13con);
    
    for ista = 4: nsta
      [coef,lag] = xcorr(detrend(optcc(:,1)), detrend(optcc(:,ista)), msftadd, 'coeff');
      [mcoef, idx] = max(coef);   % max of master raw cc
      off1i(icount,ista) = lag(idx);   % offset in samples
    end

    %%%Align and compute the RCC based on the entire win, and take that as the input signal!
    optdat = [];  % win segment of interest
    ortdat = [];
    vertdat = [];
    optdat(:, 1) = optseg(1+msftaddm: end-msftaddm, 1); % time column
    ortdat(:, 1) = ortseg(1+msftaddm: end-msftaddm, 1);
    vertdat(:, 1) = vertseg(1+msftaddm: end-msftaddm, 1);
    for ista = 1: nsta
      optdat(:, ista+1) = optseg(1+msftaddm-off1i(icount,ista): ...
        end-msftaddm-off1i(icount,ista), ista+1);
      ortdat(:, ista+1) = ortseg(1+msftaddm-off1i(icount,ista): ...
        end-msftaddm-off1i(icount,ista), ista+1);
      vertdat(:, ista+1) = vertseg(1+msftaddm-off1i(icount,ista): ...
        end-msftaddm-off1i(icount,ista), ista+1);
    end
    
    %%%for optimal comp
    sigsta = zeros(size(optdat,1), nsta);
    for ista = 1:nsta
      tmp = optdat(:,ista+1); %best aligned, filtered
      tmp = detrend(tmp);
      sigsta(:,ista) = tmp;
    end
    sigsta = detrend(sigsta(overshoot+1:end-overshoot, :));  %excluding the overshoot
    
    %%%for orthogonal comp
    sigstaort = zeros(size(ortdat,1), nsta);
    for ista = 1:nsta
      tmp = ortdat(:,ista+1); %best aligned, filtered
      tmp = detrend(tmp);
      sigstaort(:,ista) = tmp;
    end
    sigstaort = detrend(sigstaort(overshoot+1:end-overshoot, :));  %excluding the overshoot
    
    %%%for vertical comp
    sigstavert = zeros(size(vertdat,1), nsta);
    for ista = 1:nsta
      tmp = vertdat(:,ista+1); %best aligned, filtered
      tmp = detrend(tmp);
      sigstavert(:,ista) = tmp;
    end
    sigstavert = detrend(sigstavert(overshoot+1:end-overshoot, :));  %excluding the overshoot

%     optcc = STAopt(max(floor((tstbuf+1)*sps+1),1): min(floor((tedbuf-1)*sps),86400*sps), 2:nsta+1);
%     msftadd = 10/40*sps;
%     ccmid = ceil(size(optcc,1)/2);
%     ccwlen = round(size(optcc,1)-2*(msftadd+1));
%     loffmax = 4*sps/40;
%     ccmin = 0.01;  % depending on the length of trace, cc could be very low
%     iup = 1;    % times of upsampling
%     [off12con,off13con,ccali(icount),iloopoff,loopoff] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
%       ccwlen,msftadd,loffmax,ccmin,iup);
%     % if a better alignment cannot be achieved, use 0,0
%     if off12con == msftadd+1 && off13con == msftadd+1
%       off12con = 0;
%       off13con = 0;
%       fprintf('Tremor burst %d cannot be properly aligned, double-check needed \n',icount);
%     end
%     off1i(icount,1) = 0;
%     off1i(icount,2) = round(off12con);
%     off1i(icount,3) = round(off13con);
%     
%     for ista = 4: nsta
%       [coef,lag] = xcorr(detrend(optcc(:,1)), detrend(optcc(:,ista)), msftadd, 'coeff');
%       [mcoef, idx] = max(coef);   % max of master raw cc
%       off1i(icount,ista) = lag(idx);   % offset in samples
%     end
%     
%     %Align records
%     optdat = [];  % win segment of interest
%     ortdat = [];
%     optdat(:, 1:2) = STAopt(max(floor(tstbuf*sps+1),1): min(floor(tedbuf*sps),86400*sps), 1:2); % sta 1
%     ortdat(:, 1:2) = STAort(max(floor(tstbuf*sps+1),1): min(floor(tedbuf*sps),86400*sps), 1:2);
%     for ista = 2: nsta
%       optdat(:, ista+1) = STAopt(max(floor(tstbuf*sps+1)-off1i(icount,ista),1): ...
%         min(floor(tedbuf*sps)-off1i(icount,ista),86400*sps), ista+1); % sta 2
%       ortdat(:, ista+1) = STAort(max(floor(tstbuf*sps+1)-off1i(icount,ista),1): ...
%         min(floor(tedbuf*sps)-off1i(icount,ista),86400*sps), ista+1);
%     end
%     
%     %%%obtain envelope directly from signal
%     [envseg,~] = envelope(optdat(:,2:end));


    %%%obtain envelope directly from signal
    [envseg,~] = envelope(sigsta);
    [envortseg,~] = envelope(sigstaort);
    [envvertseg,~] = envelope(sigstavert);

    %get some estimate of the amplitude range/variation, using envelope
    for ip = 1: 3
      envprct(icount,ip,:) = prctile(envseg,prct(ip));
      envortprct(icount,ip,:) = prctile(envortseg,prct(ip));
      envvertprct(icount,ip,:) = prctile(envvertseg,prct(ip));
    end
    
    maxlag = 2*sps;
    
    %%%among stas 4, 5, 6, 7, between opt and opt, in order 45 46 47 56 57 67
    tmp = 0;
    for mm = 4: nsta-1
      for nn = mm+1: nsta
        tmp = tmp+1;
        [coef, lag] = xcorr(detrend(envseg(:,mm)), detrend(envseg(:,nn)), maxlag, 'coeff');
        [ccb45(icount,tmp), mind] = max(coef);
        lagb45(icount,tmp) = lag(mind);
      end
    end
    
    %%%for 4th and more stations, relative to sta 1/2/3, between opt and opt
    for mm = 4: nsta
      for nn = 1: 3
        [coef, lag] = xcorr(detrend(envseg(:,nn)), detrend(envseg(:,mm)), maxlag, 'coeff');
        [ccbij(nn,icount,mm-3), mind] = max(coef);
        lagbij(nn,icount,mm-3) = lag(mind);
      end
    end
    
    %%%among stas 1, 2, 3, between opt and opt, in order 12 13 23
    tmp = 0;
    for mm = 1: 3-1
      for nn = mm+1: 3
        tmp = tmp+1;
        [coef, lag] = xcorr(detrend(envseg(:,mm)), detrend(envseg(:,nn)), maxlag, 'coeff');
        [ccb123(icount,tmp), mind] = max(coef);
        lagb123(icount,tmp) = lag(mind);
      end
    end
    
%     f1 = plt_traceindetail(envseg,stas,25*sps);
%     title(f1.ax(1),'Envelope');
%     keyboard

    %   end
    % end
    
  end
  
  %%% save some variables
%   savefile = 'rst_envcc_dtr.mat';
%   savefile = strcat('rst_envcc',num2str(prct(end)),'_dtr.mat');
  savefile = strcat('rst_envcc_3comps',num2str(prct(end)),'_dtr.mat');
  save(strcat(rstpath, '/MAPS/',savefile), 'off1i','ccbij','lagbij','ccb123','lagb123','ccb45',...
    'lagb45','envprct','envortprct','envvertprct');
  
else
  maxlag = 2*sps;
%   savefile = strcat('rst_envcc',num2str(prct(end)),'_dtr.mat');
  savefile = strcat('rst_envcc_3comps',num2str(prct(end)),'_dtr.mat');
  load(strcat(rstpath, '/MAPS/',savefile));
end

keyboard

%% target some high-correlation bursts
% % ind = find(ccboo(:,1)>=prctile(ccboo(:,1),75) & ccboo(:,2)>=prctile(ccboo(:,2),75) & ...
% %   ccboo(:,3)>=prctile(ccboo(:,3),75));
% % ind = [20,23,59,60,80,113,116,120,134,189,194];
% % 
% ind123 = find(ccb123(:,1)>=prctile(ccb123(:,1),75) & ccb123(:,2)>=prctile(ccb123(:,2),75) & ...
%   ccb123(:,3)>=prctile(ccb123(:,3),75));
% % ind123 = [1,3,8,10,31,67,78,81,82,91,102,108,114,116,121,129,146,153,167,175];

%% scatter between envelope (amplitude) range and median
widin = 6.5;
htin = 6.8;
nrow = 2;
ncol = 2;
pltxran = [0.1 0.96]; pltyran = [0.1 0.96]; % optimal axis location
pltxsep = 0.05; pltysep = 0.04;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
color=['r';'b';'k';'c'];
sybl = ['o';'^';'v';'d'];
p=[];
xran=[-2.2 -0.2];
yran=[0 2];
for i = 1: 4
  ax=f.ax(i); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); axis(ax,'equal');
  envmed(:,i) = envprct(:,2,i);
  envran(:,i) = envprct(:,end,i)./envprct(:,1,i);
  axis(ax,[xran yran]);
  plot(ax, ax.XLim, [median(log10(envran(:,i))) median(log10(envran(:,i)))],...
    '--','Color',color(i,:),'Linew',2);
  plot(ax, [median(log10(envmed(:,i))) median(log10(envmed(:,i)))], ax.YLim,...
    '--','Color',color(i,:),'Linew',2);
  p(i)=scatter(ax,log10(envmed(:,i)),log10(envran(:,i)),30,...
      color(i,:),sybl(i,:),'filled','markeredgec','w');
  text(ax,0.98,0.95,stas(i,:),'Units','normalized','HorizontalAlignment',...
    'right','FontSize',12);
end
ylabel(f.ax(3),'log_{10}{Envelope range}');
ylabel(f.ax(1),'log_{10}{Envelope range}');
xlabel(f.ax(3),'log_{10}{Median envelope}');
xlabel(f.ax(4),'log_{10}{Median envelope}');
% legend(f.ax(3),p,stas(1:4,:),'Location','northwest');

fname = strcat('bstenvran',num2str(prct(end)),num2str(prct(1)),'vsmed.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

%% ort comp. scatter between envelope (amplitude) range and median
widin = 6.5;
htin = 6.8;
nrow = 2;
ncol = 2;
pltxran = [0.1 0.96]; pltyran = [0.1 0.96]; % optimal axis location
pltxsep = 0.05; pltysep = 0.04;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
color=['r';'b';'k';'c'];
sybl = ['o';'^';'v';'d'];
p=[];
xran=[-2.2 -0.2];
yran=[0 2];
for i = 1: 4
  ax=f.ax(i); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); axis(ax,'equal');
  envortmed(:,i) = envortprct(:,2,i);
  envortran(:,i) = envortprct(:,end,i)./envortprct(:,1,i);
  axis(ax,[xran yran]);
  plot(ax, ax.XLim, [median(log10(envortran(:,i))) median(log10(envortran(:,i)))],...
    '--','Color',color(i,:),'Linew',2);
  plot(ax, [median(log10(envortmed(:,i))) median(log10(envortmed(:,i)))], ax.YLim,...
    '--','Color',color(i,:),'Linew',2);
  p(i)=scatter(ax,log10(envortmed(:,i)),log10(envortran(:,i)),30,...
      color(i,:),sybl(i,:),'filled','markeredgec','w');
  text(ax,0.98,0.95,stas(i,:),'Units','normalized','HorizontalAlignment',...
    'right','FontSize',12);
  text(ax,0.02,0.95,'Orthogonal','Units','normalized','HorizontalAlignment',...
    'left','FontSize',12);
end
ylabel(f.ax(3),'log_{10}{Envelope range}');
ylabel(f.ax(1),'log_{10}{Envelope range}');
xlabel(f.ax(3),'log_{10}{Median envelope}');
xlabel(f.ax(4),'log_{10}{Median envelope}');
% legend(f.ax(3),p,stas(1:4,:),'Location','northwest');

fname = strcat('bstenvortran',num2str(prct(end)),num2str(prct(1)),'vsmed.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));


%% vert comp. scatter between envelope (amplitude) range and median
widin = 6.5;
htin = 6.8;
nrow = 2;
ncol = 2;
pltxran = [0.1 0.96]; pltyran = [0.1 0.96]; % optimal axis location
pltxsep = 0.05; pltysep = 0.04;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
color=['r';'b';'k';'c'];
sybl = ['o';'^';'v';'d'];
p=[];
xran=[-2.2 -0.2];
yran=[0 2];
for i = 1: 4
  ax=f.ax(i); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); axis(ax,'equal');
  envvertmed(:,i) = envvertprct(:,2,i);
  envvertran(:,i) = envvertprct(:,end,i)./envvertprct(:,1,i);
  axis(ax,[xran yran]);
  plot(ax, ax.XLim, [median(log10(envvertran(:,i))) median(log10(envvertran(:,i)))],...
    '--','Color',color(i,:),'Linew',2);
  plot(ax, [median(log10(envvertmed(:,i))) median(log10(envvertmed(:,i)))], ax.YLim,...
    '--','Color',color(i,:),'Linew',2);
  p(i)=scatter(ax,log10(envvertmed(:,i)),log10(envvertran(:,i)),30,...
      color(i,:),sybl(i,:),'filled','markeredgec','w');
  text(ax,0.98,0.95,stas(i,:),'Units','normalized','HorizontalAlignment',...
    'right','FontSize',12);
  text(ax,0.02,0.95,'Vertical','Units','normalized','HorizontalAlignment',...
    'left','FontSize',12);
end
ylabel(f.ax(3),'log_{10}{Envelope range}');
ylabel(f.ax(1),'log_{10}{Envelope range}');
xlabel(f.ax(3),'log_{10}{Median envelope}');
xlabel(f.ax(4),'log_{10}{Median envelope}');
% legend(f.ax(3),p,stas(1:4,:),'Location','northwest');

fname = strcat('bstenvvertran',num2str(prct(end)),num2str(prct(1)),'vsmed.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

%% histogram of envelope (amplitude) range
% widin = 6.5;
% htin = 7;
% nrow = 2;
% ncol = 2;
% pltxran = [0.1 0.96]; pltyran = [0.1 0.96]; % optimal axis location
% pltxsep = 0.05; pltysep = 0.05;
% f = initfig(widin,htin,nrow,ncol);
% optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
% color=['r';'b';'k';'c'];
% binw=0.05;
% p=[];
% for i = 1: 4
%   ax=f.ax(i); hold(ax,'on'); ax.Box='on'; grid(ax,'on');%  set(ax,'XScale','log');
%   envran(:,i) = envprct(:,end,i)./envprct(:,1,i);
%   histogram(ax,log10(envran(:,i)),'binw',binw,'Facec',color(i,:),'edgec','none');
% %   [N, binedge] = histcounts(log10(envran(:,i)),'binw',binw,'normalization','pdf');
% %   N=[N N(end)];
% %   p(i)=stairs(ax,binedge,N,'color',color(i,:),'LineWidth',1);
%   plot(ax,[median(log10(envran(:,i))) median(log10(envran(:,i)))], ax.YLim, ...
%     '--','Color',color(i,:),'Linew',2);
% %   xlim(ax,[0.5 1]);
%   xlim(ax,[0 2]);
%   text(ax,0.98,0.95,stas(i,:),'Units','normalized','HorizontalAlignment',...
%     'right','FontSize',12);
% end
% xlabel(f.ax(3),'log_{10}{Envelope range}');
% ylabel(f.ax(3),'Count');
% % legend(f.ax(3),p,stas(1:4,:),'Location','northwest');
% 
% fname = strcat('bstenvran',num2str(prct(end)),num2str(prct(1)),'hist.pdf');
% print(f.fig,'-dpdf',...
%   strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));


%% envelope (amplitude) range VS. burst #
% widin = 12;
% htin = 9;
% nrow = 3;
% ncol = 1;
% pltxran = [0.06 0.96]; pltyran = [0.06 0.96]; % optimal axis location
% pltxsep = 0.03; pltysep = 0.03;
% f = initfig(widin,htin,nrow,ncol);
% optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
% itime = [1 69; 70 138; 139 195];
% iyr = ['2003';'2004';'2005'];
% color = jet(nsta);
% sybl = ['o';'x';'s';'^';'d';'v';'p'];
% for isub = 1: nrow*ncol
%   ax = f.ax(isub);
%   hold(ax,'on');
%   ax.Box = 'on';
%   grid(ax,'on');
%   
% %   %bar plot range bar
% %   for ii = 1: nsta
% %     for jj = 1: nibst
% %       if jj == 1
% %         p(ii)=plot(ax,[jj jj],[log10(envprct(jj,1,ii)) log10(envprct(jj,2,ii))],'-','Marker',sybl(ii,:),...
% %           'markersize',4,'color',color(ii,:));
% %       else
% %         plot(ax,[jj jj],[log10(envprct(jj,1,ii)) log10(envprct(jj,2,ii))],'-','Marker',sybl(ii,:),'markersize',4,...
% %           'color',color(ii,:));
% %       end
% %     end
% %   end
%   
%   %line plot
%   for ii = 1: nsta
%     p(ii)=plot(ax,log10(envprct(:,3,ii)),'o-','Color',color(ii,:),'linew',1.5,'markersize',2.5);
%   end
%   for ii = 1: nsta
%     plot(ax,log10(envprct(:,1,ii)),':','Color',color(ii,:),'linew',1);   
%   end
% 
%   xlim(ax,[itime(isub,1) itime(isub,2)]);
% %   ylim(ax,[0 0.8]);
%   ylim(ax,[-2.5 0]);
%   text(ax,0.5,0.95,iyr(isub,:),'Units','normalized','FontSize',12);
%   if isub == 1
%     legend(ax,p,stas,'Location','northwest');
%   end
%   if isub == nrow
%     xlabel(ax,'Burst #');
%     ylabel(ax,'log_{10}{10--90 prctile of env}');
%   end
% end

%% envelope (amplitude) range VS. med env.
% widin = 9;
% htin = 9;
% nrow = 3;
% ncol = 1;
% pltxran = [0.08 0.96]; pltyran = [0.08 0.96]; % optimal axis location
% pltxsep = 0.03; pltysep = 0.04;
% f = initfig(widin,htin,nrow,ncol);
% optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
% color = jet(nsta);
% sybl = ['o';'^';'v';'s';'d';'p';'h'];
% 
% %just scatter plot, env median vs 90-10 percentile
% ax = f.ax(1);
% hold(ax,'on');
% ax.Box = 'on';
% grid(ax,'on');
% for ii = 1: nsta
%   p(ii)=scatter(ax,log10(envprct(:,3,ii)),log10(envprct(:,2,ii))-log10(envprct(:,1,ii)),20,...
%     color(ii,:),sybl(ii,:),'filled');
% end
% %   text(ax,0.5,0.95,iyr(isub,:),'Units','normalized','FontSize',12);
% legend(ax,p,stas,'Location','northeast');
% xlabel(ax,'log_{10}{Median env}');
% ylabel(ax,'log_{10}{90th prctile}-log_{10}{10th prctile}');
% nolabels(ax,1);
% 
% %bin according to median env, median + error bar
% ax = f.ax(2);
% hold(ax,'on');
% ax.Box = 'on';
% grid(ax,'on');
% nbin = 10;
% xcnt = zeros(nbin,nsta);
% ycnt = zeros(nbin,nsta);
% y1sig = zeros(nbin,nsta);
% x = reshape(log10(envprct(:,3,:)),[],1);
% binw = (max(x)-min(x))/nbin;
% binEdges = min(x): binw: min(x)+nbin*binw;
% for ii = 1: nsta
%   [xcnt(:,ii),ycnt(:,ii),y1sig(:,ii)] = ranybinx(log10(envprct(:,3,ii)),...
%     log10(envprct(:,2,ii))-log10(envprct(:,1,ii)),'median',[],[],binEdges);
%   p(ii)=errorbar(ax,xcnt(:,ii),ycnt(:,ii),-y1sig(:,ii),y1sig(:,ii),'vertical',sybl(ii,:),...
%     'markersize',5,'color','k','linewidth',0.8,'MarkerEdgeColor','none',...
%     'MarkerFaceColor',color(ii,:),'CapSize',4);
% end
% xlabel(ax,'Bin center');
% ylabel(ax,'Median of each bin');
% nolabels(ax,1);
% 
% %bin according to median env, median, connect with line
% ax = f.ax(3);
% hold(ax,'on');
% ax.Box = 'on';
% grid(ax,'on');
% for ii = 1: nsta
%   plot(ax,xcnt(:,ii),ycnt(:,ii),'-','Color',color(ii,:),'linew',1.5);
% end
% xlabel(ax,'Bin center');
% ylabel(ax,'Median of each bin');
% 
% for isub = 1:nrow*ncol
%   ax = f.ax(isub);
%   axis(ax,'equal');
%   xlim(ax,[-2.2 -0.2]);
%   ylim(ax,[0.4 1.2]);
%   xticks(ax,-2.2:0.1:-0.2);
%   yticks(ax,0.4:0.1:1.2);
%   longticks(ax,4);
% end
% 
% keyboard


%% burst windows for stas 1/2/3
%load the plane fit model and gradient computation between the trio stations
savefile = 'timeoff_plfitmap_160sps.mat';
load(strcat(rstpath, '/MAPS/',savefile));

nrow = 1;
ncol = 3;
widin = ncol*2.2;
htin = nrow*2.4;
pltxran = [0.10 0.96]; pltyran = [0.15 0.96]; % optimal axis location
pltxsep = 0.025; pltysep = 0.02;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

for isub = 1:3
  ax = f.ax(isub);
  hold(ax,'on');
  ax.Box = 'on';
  grid(ax,'on');
  scatter(ax,lagb123(:,isub)/sps,ccb123(:,isub),10,...
    'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
  scatter(ax,median(lagb123(:,isub)/sps),median(ccb123(:,isub)),40,'^',...
    'MarkerFaceColor','r','MarkerEdgeColor','k');
  ylim(ax,[0.0 0.8]);
  yticks(ax,0:0.2:0.8);
  xlim(ax,[-maxlag,maxlag]/sps);
  xticks(ax,-maxlag/sps: 1: maxlag/sps);
  if isub ==1
    text(ax,0.98,0.95,sprintf('%s-%s',strtrim(stas(1,:)),strtrim(stas(2,:))),'Units',...
      'normalized','HorizontalAlignment','right');
%     text(ax,0.02,0.13,strcat("$\partial\Delta{t}_{12} / \partial x'=$",...
%       sprintf(' %.3f',doffdprojloc(isub,1))),'FontSize',9,...
%       'unit','normalized','interpreter','latex','HorizontalAlignment','left');
%     text(ax,0.02,0.05,strcat("$\partial\Delta{t}_{12} / \partial y'=$",...
%       sprintf(' %.3f',doffdprojloc(isub,2))),'FontSize',9,...
%       'unit','normalized','interpreter','latex','HorizontalAlignment','left');
  elseif isub ==2
    text(ax,0.98,0.95,sprintf('%s-%s',strtrim(stas(1,:)),strtrim(stas(3,:))),'Units',...
      'normalized','HorizontalAlignment','right');
%     text(ax,0.02,0.13,strcat("$\partial\Delta{t}_{13} / \partial x'=$",...
%       sprintf(' %.3f',doffdprojloc(isub,1))),'FontSize',9,...
%       'unit','normalized','interpreter','latex','HorizontalAlignment','left');
%     text(ax,0.02,0.05,strcat("$\partial\Delta{t}_{13} / \partial y'=$",...
%       sprintf(' %.3f',doffdprojloc(isub,2))),'FontSize',9,...
%       'unit','normalized','interpreter','latex','HorizontalAlignment','left');
  else
    text(ax,0.98,0.95,sprintf('%s-%s',strtrim(stas(2,:)),strtrim(stas(3,:))),'Units',...
      'normalized','HorizontalAlignment','right');
%     text(ax,0.02,0.13,strcat("$\partial\Delta{t}_{23} / \partial x'=$",...
%       sprintf(' %.3f',doffdprojloc(isub,1))),'FontSize',9,...
%       'unit','normalized','interpreter','latex','HorizontalAlignment','left');
%     text(ax,0.02,0.05,strcat("$\partial\Delta{t}_{23} / \partial y'=$",...
%       sprintf(' %.3f',doffdprojloc(isub,2))),'FontSize',9,...
%       'unit','normalized','interpreter','latex','HorizontalAlignment','left');
  end
%   text(ax,0.98,0.05,sprintf('%.2f, %.2f',median(lagb123(:,isub)/sps),median(ccb123(:,isub))),...
%     'Units','normalized','HorizontalAlignment','right','FontSize',8);
  if isub~=1
    nolabels(ax,2);
  end
  if isub==1
    ylabel(ax,'Max CC');
  end
  xlabel(ax,'Time lag to reach max CC (s)');
  longticks(ax,1.5);
end
% supertit(f.ax(1:ncol),'cc of sig env; extended bursts; among trio stas');
fname = strcat('sigenvcctrio.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

% keyboard


%% burst windows for stas 4/5/6/7 vs. 1/2/3
%%%if off14 = a*x+b*y+c, then d(off14)/dx = a, d(off14)/dy = b, where x and y
%%%are map locations in N and E in km. If you instead want to know the gradient
%%%along a custom coordinate frame x' and y' where x' is rotated from x by
%%%theta, so [x; y] = [cos(theta) -sin(theta); sin(theta) cos(theta)]* [x'; y'].
%%%Then, the new gradient is d(off14)/dx' = a*cos(theta)+b*sin(theta);
%%%d(off14)/dy' = b*cos(theta)-a*sin(theta);
%load the plane fit model and gradient computation between the trio stations
savefile = 'timeoff_plfitmap_4thsta_160sps.mat';
load(strcat(rstpath, '/MAPS/',savefile));
modparammap = [modparammap(4,:); modparammap(1:3,:)];  %the station order is a bit different 
doff14dloc = [doff14dloc(4,:); doff14dloc(1:3,:)];  %the station order is a bit different  
doff14dprojloc = [doff14dprojloc(4,:); doff14dprojloc(1:3,:)];  %the station order is a bit different  

%%%scatter of lag and CC 
nrow = 3;
ncol = nsta-3;
widin = ncol*2.1;
htin = nrow*2.1;
pltxran = [0.06 0.96]; pltyran = [0.06 0.96]; % optimal axis location
pltxsep = 0.02; pltysep = 0.02;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

for ii = 1:nrow
  for jj = 1:ncol
    isub = (ii-1)*ncol+jj;
    ax = f.ax(isub); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
    scatter(ax,lagbij(ii,:,jj)/sps,ccbij(ii,:,jj),10,...
      'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
    scatter(ax,median(lagbij(ii,:,jj)/sps),median(ccbij(ii,:,jj)),40,'^',...
      'MarkerFaceColor','r','MarkerEdgeColor','k');
    ylim(ax,[0.0 0.7]);
    yticks(ax,0:0.1:0.7);
    xlim(ax,[-maxlag,maxlag]/sps);
    xticks(ax,-maxlag/sps: 1: maxlag/sps);
    text(ax,0.98,0.95,sprintf('%s-%s',strtrim(stas(ii,:)),strtrim(stas(jj+3,:))),'Units',...
      'normalized','HorizontalAlignment','right');
%     text(ax,0.98,0.05,sprintf('%.2f, %.2f',median(lagbij(ii,:,jj)/sps),median(ccbij(ii,:,jj))),...
%       'Units','normalized','HorizontalAlignment','right','FontSize',8);
%     text(ax,0.02,0.2,'$\frac{\partial\Delta{t}_{14}}{\partial x}=0.02$','FontSize',11,...
%       'unit','normalized','interpreter','latex','HorizontalAlignment','left');
%     if ii == 1  
%       text(ax,0.02,0.29,strcat('$\partial\Delta{t}_{14} / \partial x=$',...
%         sprintf(' %.3f',doff14dloc(jj,1))),'FontSize',9,...
%         'unit','normalized','interpreter','latex','HorizontalAlignment','left');
%       text(ax,0.02,0.21,strcat('$\partial\Delta{t}_{14} / \partial y=$',...
%         sprintf(' %.3f',doff14dloc(jj,2))),'FontSize',9,...
%         'unit','normalized','interpreter','latex','HorizontalAlignment','left');
%       text(ax,0.02,0.14,strcat("$\partial\Delta{t}_{14} / \partial x'=$",...
%         sprintf(' %.3f',doff14dprojloc(jj,1))),'FontSize',10,...
%         'unit','normalized','interpreter','latex','HorizontalAlignment','left');
%       text(ax,0.02,0.05,strcat("$\partial\Delta{t}_{14} / \partial y'=$",...
%         sprintf(' %.3f',doff14dprojloc(jj,2))),'FontSize',10,...
%         'unit','normalized','interpreter','latex','HorizontalAlignment','left');
%     elseif ii == 2  
%       text(ax,0.02,0.14,strcat("$\partial\Delta{t}_{24} / \partial x'=$",...
%         sprintf(' %.3f',doff14dprojloc(jj,1)-doffdprojloc(ii-1,1))),'FontSize',10,...
%         'unit','normalized','interpreter','latex','HorizontalAlignment','left');
%       text(ax,0.02,0.05,strcat("$\partial\Delta{t}_{24} / \partial y'=$",...
%         sprintf(' %.3f',doff14dprojloc(jj,2)-doffdprojloc(ii-1,2))),'FontSize',10,...
%         'unit','normalized','interpreter','latex','HorizontalAlignment','left');
%     elseif ii == 3  
%       text(ax,0.02,0.14,strcat("$\partial\Delta{t}_{34} / \partial x'=$",...
%         sprintf(' %.3f',doff14dprojloc(jj,1)-doffdprojloc(ii-1,1))),'FontSize',10,...
%         'unit','normalized','interpreter','latex','HorizontalAlignment','left');
%       text(ax,0.02,0.05,strcat("$\partial\Delta{t}_{34} / \partial y'=$",...
%         sprintf(' %.3f',doff14dprojloc(jj,2)-doffdprojloc(ii-1,2))),'FontSize',10,...
%         'unit','normalized','interpreter','latex','HorizontalAlignment','left');
%     end

    if jj ~= 1
      nolabels(ax,2);
    end
    if jj == 1
      ylabel(ax,'Max CC');
    end
    if ii ~= nrow
      nolabels(ax,1);
    end
    if ii == nrow
      xlabel(ax,'Time lag to reach max CC (s)');
    end
    longticks(ax,1.5);
  end
end
% supertit(f.ax(1:ncol),'cc of sig env; extended bursts; 4th stas vs. trio stas');
fname = strcat('sigenvcc4thvstrio.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

%%
for i = 1: nsta-3
  for j = 1:2
    gradratx(i,j) = doff14dprojloc(i,1)/doffdprojloc(j,1);
    gradraty(i,j) = doff14dprojloc(i,2)/doffdprojloc(j,2);
  end
end


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
supertit(f.ax(1:ncol),'cc of sig env; extended bursts; among 4th stas');




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
%     ylim(ax,[-0.3 0.2]);
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


% %% migration windows for stas 4/5/6/7 vs. 1/2/3
% %%%scatter of lag and CC 
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
%     scatter(ax,lagmij(ii,:,jj)/sps,ccmij(ii,:,jj),8,'k');
%     ylim(ax,[0.2 0.9]);
%     xlim(ax,[-maxlag,maxlag]/sps);
%     text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(jj+3,:)),strtrim(stas(ii,:))),'Units',...
%       'normalized');
%     text(ax,0.98,0.05,sprintf('%.2f, %.2f',median(lagmij(ii,:,jj)/sps),median(ccmij(ii,:,jj))),...
%       'Units','normalized','HorizontalAlignment','right');
%     if jj ~= 1
%       nolabels(ax,2);
%     end
%     if ii ~= nrow
%       nolabels(ax,1);
%     end
%     if ii == nrow && jj==1
%       xlabel(ax,'Lag (s) of max CC');
%       ylabel(ax,'Max CC');
%     end
%   end
% end
% supertit(f.ax(1:ncol),'Migrations');
% 
% %% migration windows for stas 4/5/6/7 vs. 1/2/3, minus reference
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
%       ref = mean(ccm123(:,[1 2]),2)';
%       text(ax,0.98,0.95,'ref=mean(12+13)','Units','normalized','HorizontalAlignment','right');
%     elseif ii == 2
%       ref = mean(ccm123(:,[1 3]),2)';
%       text(ax,0.98,0.95,'ref=mean(12+23)','Units','normalized','HorizontalAlignment','right');
%     else
%       ref = mean(ccm123(:,[2 3]),2)';
%       text(ax,0.98,0.95,'ref=mean(13+23)','Units','normalized','HorizontalAlignment','right');
%     end
%     scatter(ax,lagmij(ii,:,jj)/sps,ccmij(ii,:,jj)-ref,8,1:size(ccmij,2),'filled');
%     colormap(ax,'jet');
%     ylim(ax,[-0.4 0.3]);
%     xlim(ax,[-maxlag,maxlag]/sps);
%     text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(jj+3,:)),strtrim(stas(ii,:))),'Units',...
%       'normalized');
%     text(ax,0.98,0.05,sprintf('%.2f, %.2f',median(lagmij(ii,:,jj)/sps),median(ccmij(ii,:,jj)-ref)),...
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
% supertit(f.ax(1:ncol),'Migrations (-reference)');
% 
% 
% %% migration windows for stas 1/2/3
% widin = 9;
% htin = 3;
% nrow = 1;
% ncol = 3;
% pltxran = [0.06 0.96]; pltyran = [0.12 0.95]; % optimal axis location
% pltxsep = 0.05; pltysep = 0.03;
% f = initfig(widin,htin,nrow,ncol);
% optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
% 
% for isub = 1:3
%   ax = f.ax(isub);
%   hold(ax,'on');
%   ax.Box = 'on';
%   grid(ax,'on');
%   scatter(ax,lagm123(:,isub)/sps,ccm123(:,isub),8,'k');
%   ylim(ax,[0.2 0.9]);
%   xlim(ax,[-maxlag,maxlag]/sps);
%   if isub ==1
%     text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(1,:)),strtrim(stas(2,:))),'Units',...
%       'normalized');
%   elseif isub ==2
%     text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(1,:)),strtrim(stas(3,:))),'Units',...
%       'normalized');
%   else
%     text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(2,:)),strtrim(stas(3,:))),'Units',...
%       'normalized');
%   end
%   text(ax,0.98,0.05,sprintf('%.2f, %.2f',median(lagm123(:,isub)/sps),median(ccm123(:,isub))),...
%     'Units','normalized','HorizontalAlignment','right');
%   if isub~=1
%     nolabels(ax,2);
%   end
%   if isub==1
%     xlabel(ax,'Lag (s) of max CC');
%     ylabel(ax,'Max CC');
%   end
% end
% supertit(f.ax(1:ncol),'Migrations');




  
  
  
  
  
  
  