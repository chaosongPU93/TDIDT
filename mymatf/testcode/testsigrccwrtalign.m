% testsigrccwrtalign.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to 'testwletrccwrtalign.m', this is to test how CC/RCC changes wrt. 
% to different alignments between
% signals at stations. However, since this 
% would be burst dependent, as data is changing, we choose specific bursts 
% for a qualitative understanding.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/01/05
% Last modified date:   2023/01/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format short e   % Set the format to 5-digit floating point
clear
clc
close all

%% for easy testing
defval('idxbst',181); %global indices of bursts to run 
defval('normflag',0); %whether to normalize templates
defval('noiseflag',0);  %whether to use synthetic noises
defval('pltflag',1);  %whether to plot figs for each burst
defval('rccmwsec',0.5); %moving win len in sec for computing RCC

alignflg = 'wholewin';
% alignflg = 'subwin';

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
set(0,'DefaultFigureVisible','on');

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
%   'KLNB '
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

%% prepare the signal and noise windows
% dates and ets
%%% NOTE: 'dates' is the dates that the tremor is active at the region of interest, so that the
%%% waveform at the stations on these dates you are looking at mainly comes from the sources at the
%%% region of interest. We want to see if there is a noticable change in spectra during the burst
%%% windows on these dates
dates = unique(trange(:,1));
years = unique(floor(dates/1000));
nets = length(years);

sps = 160;

%filtering passband for reading data, confirmed by 'spectrabursts002_4s.m'
hisig=6.3; % this will give a similar spectral shape between template and signal
losig=1.8;

%moving window length in samples for running CC, envelope, etc.
%standard window length is about 0.5s, this is about the visual duration of the filtered and unfiltered
%template, although in fact to include as least one cycle of the main dipole of template
rccmwlen=rccmwsec*sps;
% rccmwlen=sps/2;
% rccmwlen=sps;

iii = 1;
[iets,i,j] = indofburst(trange,idxbst(iii));

year = years(iets);
datesets = dates(floor(dates/1000)==year);

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
end

tmaxi = hfdayi(:, seccol); % starting time of max power rate of half sec inside the ellipse
tmaxo = hfdayo(:, seccol); % starting time of max power rate of half sec outside the ellipse
tcnti = hfdayi(:, 15);  % the center of detecting win is the 15th col
tcnto = hfdayo(:, 15);  % the center of detecting win is the 15th col

tst = rangetemp(j,2); % start and end time of bursts
ted = rangetemp(j,3);

%how many 4-s detections fall into the burst range
indtmaxi = find(tmaxi>=tst-0.1 & tmaxi<=ted+0.1);

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

offmax = 28;
off12 = -offmax: 1: offmax;
off13 = -offmax: 1: offmax;

[X1,X2] = meshgrid(off12,off13);
tmp = [X1(:) X2(:)];
npts = length(tmp);
  
if strcmp(alignflg,'wholewin')
  %%%obtain a single best alignment based on the entire win
  %       optcc = optseg(:, 2:end);
  optcc = detrend(optseg(1+msftaddm: end-msftaddm, 2:end));
  msftadd = (round(max(abs([off12ran off13ran])))+1)*sps/40;  %+1 for safety
  ccmid = ceil(size(optcc,1)/2);
  ccwlen = round(size(optcc,1)-2*(msftadd+1));
  loffmax = 5*sps/40;
  ccmin = 0.01;  % depending on the length of trace, cc could be very low
  iup = 1;    % times of upsampling
  [off12con,off13con,ccali] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
    ccwlen,msftadd,loffmax,ccmin,iup);
  % if a better alignment cannot be achieved, use 0,0
  if off12con == msftadd+1 && off13con == msftadd+1
    off12con = 0;
    off13con = 0;
    fprintf('Tremor burst %d cannot be properly aligned, double-check needed \n',k);
  end
  off1ic(1) = 0;
  off1ic(2) = round(off12con);
  off1ic(3) = round(off13con);

  off1i = zeros(npts,3);
  off1i(:,2) = tmp(:,1)+off1ic(2);
  off1i(:,3) = tmp(:,2)+off1ic(3);

  rccmwlen=rccmwsec*sps;

  for ipt = 1: npts
    ipt
    sigsta = [];
    for ista = 1: nsta
      %cut according to the zero-crossing and the time shift from the constrained CC
      sigsta(:, ista) = optseg(1+msftaddm-off1i(ipt,ista): end-msftaddm-off1i(ipt,ista), ista+1);
      %detrend again for caution
      sigsta(:,ista) = detrend(sigsta(:,ista));
    end

    %compute running CC between 3 stations
    [ircc,rcc12,rcc13,rcc23] = RunningCC3sta(sigsta,rccmwlen);
    rccpair = [rcc12 rcc13 rcc23];
    meanmedrcc123(ipt) = mean(median(rccpair,1));
    meanmeanrcc123(ipt) = mean(mean(rccpair,1));

    cc12 = xcorr(sigsta(:,1), sigsta(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
    cc13 = xcorr(sigsta(:,1), sigsta(:,3),0,'normalized');
    cc23 = xcorr(sigsta(:,2), sigsta(:,3),0,'normalized');
    ccpair = [cc12 cc13 cc23];
    meancc123(ipt) = mean(ccpair);

%     [~,ind] = min(ccpair);
%     meancc12(ipt) = mean(ccpair(:,setdiff(1:3,ind)));
%     meanmedmedrcc12(ipt) = mean(median(rccpair(:,setdiff(1:3,ind)),1)); 
%     meanmeanrcc12(ipt) = mean(mean(rccpair(:,setdiff(1:3,ind)),1));

    meancc12(ipt) = mean(ccpair(:,[1 2]));
    meanmedrcc12(ipt) = mean(median(rccpair(:,[1 2]),1)); 
    meanmeanrcc12(ipt) = mean(mean(rccpair(:,[1 2]),1));

    meancc13(ipt) = mean(ccpair(:,[1 3]));
    meanmedrcc13(ipt) = mean(median(rccpair(:,[1 3]),1)); 
    meanmeanrcc13(ipt) = mean(mean(rccpair(:,[1 3]),1));

    meancc23(ipt) = mean(ccpair(:,[2 3]));
    meanmedrcc23(ipt) = mean(median(rccpair(:,[2 3]),1)); 
    meanmeanrcc23(ipt) = mean(mean(rccpair(:,[2 3]),1));

  end
  
elseif strcmp(alignflg,'subwin')
  %%%according the linear fitting result of off12 and off13 VS. origin time, we sort of know how
  %%%much it needs to change the overall offset to change by 1 sample
  %generate overlapping windows of the same length
  indst = floor(tstbuf*sps+1);
  inded = floor(tedbuf*sps);

  subwsec = 25;   % determined from fitting

  subwlen = subwsec*sps;
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
  
  off1iw = zeros(nwin,3);  % the best alignment between sta2, sta3 wrt sta1 for each subwin
  ccaliw = zeros(nwin,1);  % CC value using the best alignment, including 4th stas
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
    
  end
  
  off1i = zeros(npts,3);
  off1i(:,2) = tmp(:,1);
  off1i(:,3) = tmp(:,2);

  for ipt = 1: npts
    ipt
    
    rcccat = [];  % average concatenated RCC
    rccpair = [];  % concatenated RCC between each station pair, order is 12, 13, 23
    ccpair = []; % 0-lag overall cc of each subwin, between each station pair, order is 12, 13, 23
    
    for iwin = 1: nwin
%       iwin
      isubwst = windows(iwin,1)-max(floor(tstbuf*sps+1-overshoot-msftaddm),1);
      isubwed = windows(iwin,2)-max(floor(tstbuf*sps+1-overshoot-msftaddm),1);
      if iwin == 1
        isubwst = isubwst-overshoot;
      end
      if iwin == nwin
        isubwed = isubwed+overshoot;
      end
      
      %Align records
      subw = [];  % win +/-3 s, segment of interest,first 1s will be tapered
      for ista = 1: 3
        subw(:, ista) = optseg(isubwst-(off1i(ipt,ista)+off1iw(iwin,ista)): ...
          isubwed-(off1i(ipt,ista)+off1iw(iwin,ista)), ista+1);
        subw(:, ista) = detrend(subw(:, ista));
      end
      
      %compute running CC between 3 stations
      [irccw,rccw12,rccw13,rccw23] = RunningCC3sta(subw,rccmwlen);
      rccw = (rccw12+rccw13+rccw23)/3;
      rcccat = [rcccat; rccw];
      rccpair = [rccpair; rccw12 rccw13 rccw23];
        
      ccw12 = xcorr(subw(:,1), subw(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
      ccw13 = xcorr(subw(:,1), subw(:,3),0,'normalized');
      ccw23 = xcorr(subw(:,2), subw(:,3),0,'normalized');
      ccpair = [ccpair; ccw12 ccw13 ccw23];
      
    end
    meanmedrcc123(ipt) = mean(median(rccpair,1));
    meanmeanrcc123(ipt) = mean(mean(rccpair,1));
    meancc123(ipt) = mean(mean(ccpair,1));

%     [~,ind] = min(mean(ccpair,1));
%     cc(ipt) = mean(mean(ccpair(:,setdiff(1:3,ind)), 1));
%     medrcc(ipt) = mean(median(rccpair(:,setdiff(1:3,ind)),1)); 
%     meanrcc(ipt) = mean(mean(rccpair(:,setdiff(1:3,ind)),1)); 
    
    meancc12(ipt) = mean(mean(ccpair(:,[1 2]), 1));
    meanmedrcc12(ipt) = mean(median(rccpair(:,[1 2]),1));
    meanmeanrcc12(ipt) = mean(mean(rccpair(:,[1 2]),1));
    
    meancc13(ipt) = mean(mean(ccpair(:,[1 3]), 1));
    meanmedrcc13(ipt) = mean(median(rccpair(:,[1 3]),1)); 
    meanmeanrcc13(ipt) = mean(mean(rccpair(:,[1 3]),1));

    meancc23(ipt) = mean(mean(ccpair(:,[2 3]), 1));
    meanmedrcc23(ipt) = mean(median(rccpair(:,[2 3]),1)); 
    meanmeanrcc23(ipt) = mean(mean(rccpair(:,[2 3]),1));

  end

end

meanmedrcc123 = reshape(meanmedrcc123,length(off13),length(off12));
meanmeanrcc123 = reshape(meanmeanrcc123,length(off13),length(off12));
meanmedrcc12 = reshape(meanmedrcc12,length(off13),length(off12));
meanmeanrcc12 = reshape(meanmeanrcc12,length(off13),length(off12));
meanmedrcc13 = reshape(meanmedrcc13,length(off13),length(off12));
meanmeanrcc13 = reshape(meanmeanrcc13,length(off13),length(off12));
meanmedrcc23 = reshape(meanmedrcc23,length(off13),length(off12));
meanmeanrcc23 = reshape(meanmeanrcc23,length(off13),length(off12));
meancc123 = reshape(meancc123,length(off13),length(off12));
meancc12 = reshape(meancc12,length(off13),length(off12));
meancc13 = reshape(meancc13,length(off13),length(off12));
meancc23 = reshape(meancc23,length(off13),length(off12));

%the cut-out boundary of 4-s detections
ftrans = 'interpchao';
x0 = 0.2; 
y0 = 0.2;
semia = 1.75;
semib = 1.25;
angrot = 45;
[xcut, ycut] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);
%what's the shape of it in sample space?
offcut = space2off002([xcut-x0, ycut-y0],sps,ftrans,0);

%%%range that encloses about 90 percent of LFE detections aftering shifting back
%%%to the origin
load('90thprcrangeoflfes.mat');

%% load the LFE catalog, whole-win
if strcmp(alignflg,'wholewin')
  savefile = 'deconv1win_stats4th_allbstsig.mat';
elseif strcmp(alignflg,'subwin')
  savefile = 'deconv_stats4th_allbstsig.mat';
end
allsig= load(strcat(rstpath, '/MAPS/',savefile));
imp = allsig.allbstsig.impindepall;
%convert time offset to relative loc
ftrans = 'interpchao';
[imploc, ~] = off2space002(imp(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

nsrc = allsig.allbstsig.nsrc;
% propang = allsig.allbstsig.propang;
% proppear = allsig.allbstsig.proppear;
if strcmp(alignflg,'wholewin')
  rccsrc = allsig.allbstsig.rccsrcall;
elseif strcmp(alignflg,'subwin') 
  rccsrc = allsig.allbstsig.rcccatsrcall;
end

%time and loc of LFE of that burst, from 25-win detection
ist = sum(nsrc(1:(idxbst(iii)-1)))+1;
ied = ist+nsrc(idxbst(iii))-1;
lfe = imp(ist:ied, :);
off1ibst = allsig.allbstsig.off1i;
lfe(:,7:8) = lfe(:,7:8)-repmat([off1ibst(idxbst(iii),2) off1ibst(idxbst(iii),3)],size(lfe,1),1);

lfeloc = imploc(ist:ied, :);
rcc = rccsrc(ist:ied, 1); %cat RCC at the average src arrival of 3 stas

%%
%%%plot to show how does the RCC/CC change wrt. the diff alignment between sta pairs 12 and 13
%%%note that the result is gonna be burst dependent, as data is changing
%%%in contrast, if using templates for similar analysis, result would be invariant
widin = 6.5;  % maximum width allowed is 8.5 inches
htin = 6.5;   % maximum height allowed is 11 inches
nrow = 2;
ncol = 2;
f1=initfig(widin,htin,nrow,ncol);

figxran = [0.1 0.98]; figyran = [0.1 0.97];
figxsep = 0.05; figysep = 0.08;
optaxpos(f1,nrow,ncol,figxran,figyran,figxsep,figysep);

matplt = meancc12;
ax=f1.ax(1);
hold(ax,'on'); ax.Box = 'on'; 
imagesc(ax,off12/sps,off13/sps,matplt);
plot(ax,offcut(:,1)/sps,offcut(:,2)/sps,'k-','linew',1.5);
loffm = 6*4;  %at 160 Hz
detoffm = [-loffm loffm; -loffm -loffm; loffm -loffm; loffm loffm; -loffm loffm];
plot(ax,detoffm(:,1)/sps,detoffm(:,2)/sps,'-','linew',1.5,'color',[.5 .5 .5]);
%plot the range of 90 percent of detected LFEs
plot(ax,conmat10th(:,3),conmat10th(:,4),'r-','LineWidth',1.5);
[C,h] = contour(ax,X1/sps,X2/sps,matplt,'k-','ShowText','on','LineWidth',0.5);
clabel(C,h,'FontSize',8);
% scatter(ax,lfe(:,7)/sps,lfe(:,8)/sps,10,[.7 .7 .7],'filled','MarkerEdgeColor','k');
% plot(ax,minmax(off12)/sps,minmax(off13)/sps,'--','Color',[.3 .3 .3],'linew',1);
ind = find(matplt==max(matplt,[],'all'));
[sub13, sub12] = ind2sub(size(matplt),ind);
scatter(ax,off12(sub12)/sps,off13(sub13)/sps,30,'k^','linew',1.5);
% colormap(ax,'jet');
colormap(ax, flipud(colormap(ax,'kelicol')));
% caxis(ax,[-1 1]);
% c=colorbar(ax);
% c.Label.String = 'mean overall CC of pairs 12 and 13';
title(ax,'Average CC of pairs 12 and 13','FontSize',10,'FontWeight','bold');  %,'FontWeight','normal'
text(ax,0.02,0.95,'a','FontSize',11,'unit','normalized','EdgeColor','k',...
  'Margin',1,'color','k','backgroundcolor','w');
axis(ax,'equal','tight');
xlim(ax,minmax(off12)/sps);
ylim(ax,minmax(off13)/sps);
xticks(ax,-0.1: 0.1: 0.1);
yticks(ax,-0.1: 0.1: 0.1);
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
% xlabel(ax,sprintf('\\Delta{t}_{12} (s)'));
ylabel(ax,sprintf('\\Delta{t}_{13} (s)'));
longticks(ax,1);

matplt = meancc13;
ax=f1.ax(2); hold(ax,'on'); ax.Box = 'on'; 
imagesc(ax,off12/sps,off13/sps,matplt);
% plot(ax,offcut(:,1)/sps,offcut(:,2)/sps,'b-','linew',1.5);
[C,h] = contour(ax,X1/sps,X2/sps,matplt,'k-','ShowText','on','LineWidth',0.5);
clabel(C,h,'FontSize',8);
% plot(ax,minmax(off12)/sps,minmax(off13)/sps,'--','Color',[.3 .3 .3],'linew',1);
ind = find(matplt==max(matplt,[],'all'));
[sub13, sub12] = ind2sub(size(matplt),ind);
scatter(ax,off12(sub12)/sps,off13(sub13)/sps,30,'k^','linew',1.5);
% colormap(ax,'jet');
colormap(ax, flipud(colormap(ax,'kelicol')));
% caxis(ax,[-1 1]);
% c=colorbar(ax);
% c.Label.String = 'median of mean RCC of 3 pairs';
title(ax,'Average CC of pairs 12 and 23','FontSize',10,'FontWeight','normal');
text(ax,0.02,0.95,'b','FontSize',11,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w'); %,'backgroundcolor','w'
axis(ax,'equal','tight');
xlim(ax,minmax(off12)/sps);
ylim(ax,minmax(off13)/sps);
xticks(ax,-0.1: 0.1: 0.1);
yticks(ax,-0.1: 0.1: 0.1);
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
% xlabel(ax,sprintf('PGC-SSIB offset (samples at %d Hz)',sps),'FontSize',11);
% ylabel(ax,sprintf('PGC-SILB offset (samples at %d Hz)',sps),'FontSize',11);
% xlabel(ax,sprintf('\\Delta{t}_{12} (s)'));
% ylabel(ax,sprintf('\\Delta{t}_{13} (s)'));
longticks(ax,1);

matplt = meancc23;
ax=f1.ax(3); hold(ax,'on'); ax.Box = 'on'; 
imagesc(ax,off12/sps,off13/sps,matplt);
% plot(ax,offcut(:,1)/sps,offcut(:,2)/sps,'b-','linew',1.5);
[C,h] = contour(ax,X1/sps,X2/sps,matplt,'k-','ShowText','on','LineWidth',0.5);
clabel(C,h,'FontSize',8);
% scatter(ax,lfe(:,7)/sps,lfe(:,8)/sps,10,[.7 .7 .7],'filled','MarkerEdgeColor','k');
% plot(ax,minmax(off12)/sps,minmax(off13)/sps,'--','Color',[.3 .3 .3],'linew',1);
ind = find(matplt==max(matplt,[],'all'));
[sub13, sub12] = ind2sub(size(matplt),ind);
scatter(ax,off12(sub12)/sps,off13(sub13)/sps,30,'k^','linew',1.5);
% colormap(ax,'jet');
colormap(ax, flipud(colormap(ax,'kelicol')));
% caxis(ax,[-1 1]);
% c=colorbar(ax);
% c.Label.String = 'median of mean RCC of pairs 12 and 13';
title(ax,'Average CC of pairs 13 and 23','FontSize',10,'FontWeight','normal');
text(ax,0.02,0.95,'c','FontSize',11,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
axis(ax,'equal','tight');
xlim(ax,minmax(off12)/sps);
ylim(ax,minmax(off13)/sps);
xticks(ax,-0.1: 0.1: 0.1);
yticks(ax,-0.1: 0.1: 0.1);
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
xlabel(ax,sprintf('\\Delta{t}_{12} (s)'));
ylabel(ax,sprintf('\\Delta{t}_{13} (s)'));
longticks(ax,1);

matplt = meancc123;
ax=f1.ax(4);
hold(ax,'on'); ax.Box = 'on'; 
imagesc(ax,off12/sps,off13/sps,matplt);
% plot(ax,offcut(:,1)/sps,offcut(:,2)/sps,'b-','linew',1.5);
[C,h] = contour(ax,X1/sps,X2/sps,matplt,'k-','ShowText','on','LineWidth',0.5);
clabel(C,h,'FontSize',8);
% conmat = contour(ax,X1,X2,matplt,-0.15:0.01:0.3,'k-','ShowText','on'); %
% plot(ax,minmax(off12)/sps,minmax(off13)/sps,'--','Color',[.3 .3 .3],'linew',1);
ind = find(matplt==max(matplt,[],'all'));
[sub13, sub12] = ind2sub(size(matplt),ind);
scatter(ax,off12(sub12)/sps,off13(sub13)/sps,30,'k^','linew',1.5);
% colormap(ax,'jet');
colormap(ax, flipud(colormap(ax,'kelicol')));
% caxis(ax,[-1 1]);
% c=colorbar(ax);
% c.Label.String = 'mean overall CC of 3 pairs';
title(ax,'Average CC of all 3 pairs','FontSize',10,'FontWeight','normal');
text(ax,0.02,0.95,'d','FontSize',11,'unit','normalized','EdgeColor','k',...
  'Margin',1,'color','k','backgroundcolor','w');
axis(ax,'equal','tight');
xlim(ax,minmax(off12)/sps);
ylim(ax,minmax(off13)/sps);
xticks(ax,-0.1: 0.1: 0.1);
yticks(ax,-0.1: 0.1: 0.1);
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
xlabel(ax,sprintf('\\Delta{t}_{12} (s)'));
% ylabel(ax,sprintf('\\Delta{t}_{13} (s)'));
longticks(ax,1);

% orient(f1.fig,'landscape');
fname = 'sigccwrtoffset.pdf';
print(f1.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

keyboard

%% choose a set of offset12, 13 with equal CC/RCC
ccval = 0.1;
[~,indcol] = find(abs(conmat(1,:)-ccval)<=1e-6);
nind = length(indcol);
for i = 1: nind
  erroff = conmat(:,indcol(i)+1: indcol(i)+conmat(2,indcol(i)))';
end
%fit an ellipse to it
ellfit = fit_ellipse(erroff(:,1),erroff(:,2));

%%%add the indication of the error ellipse
mmmax = offmax+1;
nnmax = offmax+1;
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

errloc = [];
errloc(:,1) = F1(erroff(:,1),erroff(:,2));
errloc(:,2) = F2(erroff(:,1),erroff(:,2));


widin = 8; htin = 4;
nrow = 1; ncol = 2;
[f2] = initfig(widin,htin,nrow,ncol);
xran = [0.1 0.95]; yran = [0.1 0.95];
xsep = 0.08; ysep = 0.05;
optaxpos(f2,nrow,ncol,xran,yran,xsep,ysep);

ax=f2.ax(1);
hold(ax,'on'); ax.Box = 'on'; 
plot(ax,[-mmmax mmmax],[-nnmax nnmax],'--','Color',[.3 .3 .3],'linew',1);
plot(ax,erroff(:,1),erroff(:,2),'k','linew',1.5);
x0 = ellfit.X0_in; 
y0 = ellfit.Y0_in;
semia = ellfit.a; 
semib = ellfit.b;
rotang = -rad2deg(ellfit.phi);
rotcnt = [x0, y0];
[xe, ye] = ellipse_chao(x0,y0,semia,semib,0.01,rotang,rotcnt);
plot(ax,xe,ye,'r-');
[xa,ya] = complex_rot(x0,y0+semib,rotang,rotcnt);
plot(ax,[x0 xa],[y0 ya],'b--','linew',1);
[xb,yb] = complex_rot(x0+semia,y0,rotang,rotcnt);
plot(ax,[x0 xb],[y0 yb],'b--','linew',1);
text(ax,0.95,0.1,strcat(sprintf('%.1f; %.1f; %.1f',semia,semib,rotang),' {\circ}'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',10);
xlabel(ax,sprintf('PGC-SSIB offset (samples at %d Hz)',sps),'FontSize',11);
ylabel(ax,sprintf('PGC-SILB offset (samples at %d Hz)',sps),'FontSize',11);
axis(ax,'equal','tight');
axis(ax,[-mmmax mmmax -nnmax nnmax]);
title(ax,sprintf('Contour of CC=%.2f',ccval));

ax=f2.ax(2);
hold(ax,'on'); ax.Box = 'on'; 
plot(ax,errloc(:,1),errloc(:,2),'k','linew',1.5);
xlabel(ax,'E (km)','FontSize',11);
ylabel(ax,'N (km)','FontSize',11);
axis(ax,'equal','tight');
% plot(ax,[-1 1],[-1 1],'--','Color',[.3 .3 .3],'linew',1);
% axis(ax,[-1 1 -1 1]);

keyboard
%% what is the shape like for sqrt(off12^2+off13^2+off23^2)==const
const = 8.1;

widin = 8; htin = 4;
nrow = 1; ncol = 2;
[f3] = initfig(widin,htin,nrow,ncol);
xran = [0.1 0.95]; yran = [0.1 0.95];
xsep = 0.08; ysep = 0.05;
optaxpos(f3,nrow,ncol,xran,yran,xsep,ysep);

ax=f3.ax(1);
hold(ax,'on'); ax.Box = 'on'; 
plot(ax,[-mmmax mmmax],[-nnmax nnmax],'--','Color',[.3 .3 .3],'linew',1);
syms x y
func1 = @(x,y) x.^2+y.^2-const^2;
fimplicit(ax,func1,'color',[.3 .3 .3],'linew',1);
func2 = @(x,y) x.^2+y.^2-x.*y-const^2/2;
fp=fimplicit(ax,func2,'k','linew',1.5);
%fit an ellipse to it
erroff = [reshape(fp.XData,[],1) reshape(fp.YData,[],1)];
ellfit = fit_ellipse(erroff(:,1),erroff(:,2));
x0 = ellfit.X0_in; 
y0 = ellfit.Y0_in;
semia = ellfit.a; 
semib = ellfit.b;
rotang = -rad2deg(ellfit.phi);
rotcnt = [x0, y0];
[xe, ye] = ellipse_chao(x0,y0,semia,semib,0.01,rotang,rotcnt);
plot(ax,xe,ye,'r-');
[xa,ya] = complex_rot(x0,y0+semib,rotang,rotcnt);
plot(ax,[x0 xa],[y0 ya],'b--','linew',1);
[xb,yb] = complex_rot(x0+semia,y0,rotang,rotcnt);
plot(ax,[x0 xb],[y0 yb],'b--','linew',1);
text(ax,0.95,0.1,strcat(sprintf('%.1f; %.1f; %.1f',semia,semib,rotang),' {\circ}'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',10);
xlabel(ax,sprintf('PGC-SSIB offset (samples at %d Hz)',sps),'FontSize',11);
ylabel(ax,sprintf('PGC-SILB offset (samples at %d Hz)',sps),'FontSize',11);
axis(ax,'equal','tight');
axis(ax,[-mmmax mmmax -nnmax nnmax]);
xticks(ax,-mmmax: 2: mmmax);
yticks(ax,-nnmax: 2: nnmax);
title(ax,strcat('$\sqrt(x^2+y^2-(y-x)^2)=$',sprintf(' %.1f',const)),'Interpreter','latex');

ax=f3.ax(2);
hold(ax,'on'); ax.Box = 'on';
errloc = [];
errloc(:,1) = F1(erroff(:,1),erroff(:,2));
errloc(:,2) = F2(erroff(:,1),erroff(:,2));
plot(ax,errloc(:,1),errloc(:,2),'k','linew',1.5);
plot(ax,[-1 1],[-1 1],'--','Color',[.3 .3 .3],'linew',1);
xlabel(ax,'E (km)','FontSize',11);
ylabel(ax,'N (km)','FontSize',11);
axis(ax,'equal','tight');
axis(ax,[-1 1 -1 1]);







