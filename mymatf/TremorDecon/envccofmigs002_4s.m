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

ccmij = zeros(3, nmig, nsta-3);
lagmij = zeros(3, nmig, nsta-3);
ccm123 = zeros(nmig, 3);
lagm123 = zeros(nmig, 3);

ccbij = zeros(3, nbst, nsta-3);
lagbij = zeros(3, nbst, nsta-3);
ccb123 = zeros(nbst, 3);
lagb123 = zeros(nbst, 3);

for i = 1: length(dates)  % dates in each ets
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
 
  %%%obtain the mean envelope
  tmp = detrend(STAopt(:,2:end));
%   fractap = 20*sps/size(STAopt,1); % if fractap is >=1, n-point von Hann window is returned
%   ptstap = fractap/2*size(STAopt,1); % if fractap is >=1, n-point von Hann window is returned
%   w = tukeywin(size(STAopt,1),fractap);
%   tmp = w.* tmp;
  tmp = detrend(tmp); %detrend again for caution
  [envup,~] = envelope(tmp);

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
  for j = 1: size(rangetemp,1)
    tst = rangetemp(j,2); % start and end time of bursts
    ted = rangetemp(j,3);
    
    %how many 4-s detections fall into the burst range
    indtmaxi = find(tmaxi>=tst-0.1 & tmaxi<=ted+0.1);
    
    %%%%Use the start and end of the 4-s detecting window 
%     tstbuf = min(tcnti(indtmaxi)-2)-4;  %with some buffer
%     tedbuf = max(tcnti(indtmaxi)+2)+4;
    tstbuf = min(tcnti(indtmaxi)-2);  %no buffer
    tedbuf = max(tcnti(indtmaxi)+2);
    
    envseg = envup(max(floor(tstbuf*sps)+1, 1): min(floor(tedbuf*sps),86400*sps), :);
    
    maxlag = 2*sps;
    
    for jj = 4: nsta
      [coef, lag] = xcorr(envseg(:,jj), envseg(:,1), maxlag, 'coeff');
      [ccbij(1,j+k-size(rangetemp,1),jj-3),mind] = max(coef);
      lagbij(1,j+k-size(rangetemp,1),jj-3) = lag(mind);
      
      [coef, lag] = xcorr(envseg(:,jj), envseg(:,2), maxlag, 'coeff');
      [ccbij(2,j+k-size(rangetemp,1),jj-3),mind] = max(coef);
      lagbij(2,j+k-size(rangetemp,1),jj-3) = lag(mind);
      
      [coef, lag] = xcorr(envseg(:,jj), envseg(:,3), maxlag, 'coeff');
      [ccbij(3,j+k-size(rangetemp,1),jj-3),mind] = max(coef);
      lagbij(3,j+k-size(rangetemp,1),jj-3) = lag(mind);      
    end
    
    [coef, lag] = xcorr(envseg(:,1), envseg(:,2), maxlag, 'coeff');
    [ccb123(j+k-size(rangetemp,1),1),mind] = max(coef);
    lagb123(j+k-size(rangetemp,1),1) = lag(mind);
    [coef, lag] = xcorr(envseg(:,1), envseg(:,3), maxlag, 'coeff');
    [ccb123(j+k-size(rangetemp,1),2),mind] = max(coef);
    lagb123(j+k-size(rangetemp,1),2) = lag(mind);
    [coef, lag] = xcorr(envseg(:,2), envseg(:,3), maxlag, 'coeff');
    [ccb123(j+k-size(rangetemp,1),3),mind] = max(coef);
    lagb123(j+k-size(rangetemp,1),3) = lag(mind);  
      
  end
      
end

%% migration windows for stas 4/5/6/7 vs. 1/2/3
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
    scatter(ax,lagmij(ii,:,jj)/sps,ccmij(ii,:,jj),8,'k');
    ylim(ax,[0.2 0.9]);
    xlim(ax,[-maxlag,maxlag]/sps);
    text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(jj+3,:)),strtrim(stas(ii,:))),'Units',...
      'normalized');
    text(ax,0.98,0.05,sprintf('%.2f, %.2f',median(lagmij(ii,:,jj)/sps),median(ccmij(ii,:,jj))),...
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
  end
end
supertit(f.ax(1:ncol),'Migrations');

%% migration windows for stas 4/5/6/7 vs. 1/2/3, minus reference
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
    if ii == 1
      ref = mean(ccm123(:,[1 2]),2)';
      text(ax,0.98,0.95,'ref=mean(12+13)','Units','normalized','HorizontalAlignment','right');
    elseif ii == 2
      ref = mean(ccm123(:,[1 3]),2)';
      text(ax,0.98,0.95,'ref=mean(12+23)','Units','normalized','HorizontalAlignment','right');
    else
      ref = mean(ccm123(:,[2 3]),2)';
      text(ax,0.98,0.95,'ref=mean(13+23)','Units','normalized','HorizontalAlignment','right');
    end
    scatter(ax,lagmij(ii,:,jj)/sps,ccmij(ii,:,jj)-ref,8,1:size(ccmij,2),'filled');
    colormap(ax,'jet');
    ylim(ax,[-0.4 0.3]);
    xlim(ax,[-maxlag,maxlag]/sps);
    text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(jj+3,:)),strtrim(stas(ii,:))),'Units',...
      'normalized');
    text(ax,0.98,0.05,sprintf('%.2f, %.2f',median(lagmij(ii,:,jj)/sps),median(ccmij(ii,:,jj)-ref)),...
      'Units','normalized','HorizontalAlignment','right');
    if jj ~= 1
      nolabels(ax,2);
    end
    if ii ~= nrow
      nolabels(ax,1);
    end
    if ii == nrow && jj==1
      xlabel(ax,'Lag (s) of max CC');
      ylabel(ax,'Max CC - ref CC');
      colorbar(ax,'Location','west');
    end
  end
end
supertit(f.ax(1:ncol),'Migrations (-reference)');


%% migration windows for stas 1/2/3
widin = 9;
htin = 3;
nrow = 1;
ncol = 3;
pltxran = [0.06 0.96]; pltyran = [0.12 0.95]; % optimal axis location
pltxsep = 0.05; pltysep = 0.03;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

for isub = 1:3
  ax = f.ax(isub);
  hold(ax,'on');
  ax.Box = 'on';
  grid(ax,'on');
  scatter(ax,lagm123(:,isub)/sps,ccm123(:,isub),8,'k');
  ylim(ax,[0.2 0.9]);
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
  text(ax,0.98,0.05,sprintf('%.2f, %.2f',median(lagm123(:,isub)/sps),median(ccm123(:,isub))),...
    'Units','normalized','HorizontalAlignment','right');
  if isub~=1
    nolabels(ax,2);
  end
  if isub==1
    xlabel(ax,'Lag (s) of max CC');
    ylabel(ax,'Max CC');
  end
end
supertit(f.ax(1:ncol),'Migrations');


%% burst windows for stas 4/5/6/7 vs. 1/2/3
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
    scatter(ax,lagbij(ii,:,jj)/sps,ccbij(ii,:,jj),8,'k');
    ylim(ax,[0.4 0.9]);
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
  end
end
supertit(f.ax(1:ncol),'Bursts');

%% burst windows for stas 4/5/6/7 vs. 1/2/3, minus reference
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
    if ii == 1
      ref = mean(ccb123(:,[1 2]),2)';
      text(ax,0.98,0.95,'ref=mean(12+13)','Units','normalized','HorizontalAlignment','right');
    elseif ii == 2
      ref = mean(ccb123(:,[1 3]),2)';
      text(ax,0.98,0.95,'ref=mean(12+23)','Units','normalized','HorizontalAlignment','right');
    else
      ref = mean(ccb123(:,[2 3]),2)';
      text(ax,0.98,0.95,'ref=mean(13+23)','Units','normalized','HorizontalAlignment','right');
    end
    scatter(ax,lagbij(ii,:,jj)/sps,ccbij(ii,:,jj)-ref,8,1:size(ccbij,2),'filled');
    colormap(ax,'jet');
    ylim(ax,[-0.3 0.2]);
    xlim(ax,[-maxlag,maxlag]/sps);
    text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(jj+3,:)),strtrim(stas(ii,:))),'Units',...
      'normalized');
    text(ax,0.98,0.05,sprintf('%.2f, %.2f',median(lagbij(ii,:,jj)/sps),median(ccbij(ii,:,jj)-ref)),...
      'Units','normalized','HorizontalAlignment','right');
    if jj ~= 1
      nolabels(ax,2);
    end
    if ii ~= nrow
      nolabels(ax,1);
    end
    if ii == nrow && jj==1
      xlabel(ax,'Lag (s) of max CC');
      ylabel(ax,'Max CC - ref CC');
      colorbar(ax,'Location','west');
    end
  end
end
supertit(f.ax(1:ncol),'Bursts (-reference)');


%% burst windows for stas 1/2/3
widin = 9;
htin = 3;
nrow = 1;
ncol = 3;
pltxran = [0.06 0.96]; pltyran = [0.12 0.95]; % optimal axis location
pltxsep = 0.05; pltysep = 0.03;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

for isub = 1:3
  ax = f.ax(isub);
  hold(ax,'on');
  ax.Box = 'on';
  grid(ax,'on');
  scatter(ax,lagb123(:,isub)/sps,ccb123(:,isub),8,'k');
  ylim(ax,[0.4 0.9]);
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
end
supertit(f.ax(1:ncol),'Bursts');




  
  
  
  
  
  
  