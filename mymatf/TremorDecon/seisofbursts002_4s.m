% seisofbursts002_4s.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The next script to run after 'tremorbursts002_4s.m', although technically
% 'tremorbursts002_4s.m' is the end of finding the target windows containing
% the high-density tremor detections around 002 region. However, we need this
% script to plot and see how is the seismogram like for these windows, in 
% order to know, eg., need to break into subwindows for later 
% deconvolution? How to approximate noise for that window? Do you need diff.
% strategy of deconvolution for diff. windows with diff. amplitude variation
% characteristics?
% 
% --it's guaranteed that if the requested time before data win is < t_max=15s 
%   used for groupingburst, there will be no 4-s detections that occurred in
%   the region of interest with that time,
%   otherwise it will be grouped into the same burst.
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/03/21
% Last modified date:   2022/03/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
format short e   % Set the format to 5-digit floating point
clear
clc
close all

% set(0,'DefaultFigureVisible','on');
set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

%%%Flag to indicate if it is necessary to recalculate everything
flagrecalc = 0;
% flagrecalc = 1;

if flagrecalc

[scrsz, res] = pixelperinch(1);

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

% dates and ets
%%% NOTE: 'dates' is the dates that the tremor is active at the region of interest, so that the
%%% waveform at the stations on these dates you are looking at mainly comes from the sources at the
%%% region of interest. We want to see if there is a noticable change in spectra during the burst
%%% windows on these dates
dates = unique(trange(:,1));
years = unique(floor(dates/1000));
nets = length(years);


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


%% read daily data, break into windows of segments, plot
%filtering passband for reading data
hisig=6.3; % this will give a similar spectral shape between template and signal
losig=1.8;

sps = 160;

%moving window length in samples for running CC, envelope, etc.
mwlen=0.5*sps;
% mwlen=sps;

%settings for each seismogram figure f1
nrow1 = 8; % rows and cols of subplots in each figure
ncol1 = 2; 
widin1 = 16; % size of each figure
htin1 = 10;
pltxran1 = [0.06 0.96]; pltyran1 = [0.06 0.96]; % optimal axis location
pltxsep1 = 0.03; pltysep1 = 0.03;
pltwlen1 = 25;   % length of record to plot for each subplot, in sec 
ifig1 = 0; % figure number

%settings for running CC VS running envelope for each burst f2
nrow2 = 3; % rows and cols of subplots in each figure
ncol2 = 5; 
widin2 = 15; % size of each figure
htin2 = 9;
pltxran2 = [0.06 0.96]; pltyran2 = [0.06 0.96]; % optimal axis location
pltxsep2 = 0.03; pltysep2 = 0.03;
ifig2 = 100; % figure number
%initialize the figure of rcc VS renv
ifig2 = ifig2+1;
f2 = initfig(widin2,htin2,nrow2,ncol2,ifig2);
axpos2 = optaxpos(f2,nrow2,ncol2,pltxran2,pltyran2,pltxsep2,pltysep2);
axtit = [f2.ax(1) f2.ax(2) f2.ax(3) f2.ax(4) f2.ax(5)];
supertit(axtit, sprintf('Fig %s, Passband: %.1f-%.1f Hz, cclen: %.1f s, rcc vs. renv',...
  num2zeropadstr(ifig2-100, 2),losig,hisig,mwlen/sps),10);
xlabel(f2.ax((nrow2-1)*ncol2+1),'Running envelope','FontSize',10);
ylabel(f2.ax((nrow2-1)*ncol2+1),'Running CC','FontSize',10);
nsubtot2 = 0;

%settings for the empirical CDF of running envelope for each burst f3, same as f2
ifig3 = 200; % figure number
%initialize the figure of rcc VS renv
ifig3 = ifig3+1;
f3 = initfig(widin2,htin2,nrow2,ncol2,ifig3);
axpos3 = optaxpos(f3,nrow2,ncol2,pltxran2,pltyran2,pltxsep2,pltysep2);
axtit = [f3.ax(1) f3.ax(2) f3.ax(3) f3.ax(4) f3.ax(5)];
supertit(axtit, sprintf('Fig %s, Passband: %.1f-%.1f Hz, cclen: %.1f s, renv CDF',...
  num2zeropadstr(ifig3-200, 2),losig,hisig,mwlen/sps),10);
xlabel(f3.ax((nrow2-1)*ncol2+1),'Running envelope','FontSize',10);
ylabel(f3.ax((nrow2-1)*ncol2+1),'Empirical CDF','FontSize',10);

k = 0;  % counting the burst windows

off1i = zeros(size(trange,1),3);  % the best alignment between sta2, sta3 wrt sta1
ccali = zeros(size(trange,1),1);  % CC using the best alignment

% keyboard
runall = cell(nets,1);  % store the results of running CC, envelope and abs amplitude for all ETS 

for iets = 1: nets
  % dates in each ets
  year = years(iets);
  datesets = dates(floor(dates/1000)==year);
  
  runtmp = [];   % for each ETS, store the results of running CC, envelope and abs amplitude
  
  for i = 1: length(datesets)
    
    date = datesets(i);
    jday = floor(date-year*1000);
    a = jul2dat(year,jday);
    mo = a(1);
    if mo == 9
      mo = 'Sep.';
    elseif mo == 7
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
    
    ifig1 = ifig1+1;
    f1 = initfig(widin1,htin1,nrow1,ncol1,ifig1);
    axpos1 = optaxpos(f1,nrow1,ncol1,pltxran1,pltyran1,pltxsep1,pltysep1);
    axtit = [f1.ax(1) f1.ax(2)];
    supertit(axtit, sprintf('Fig %s, %d, %s %s, %s, Passband: %.1f-%.1f Hz, wlen: %d s',...
      num2zeropadstr(ifig1, 2),date,mo,dy,yr,losig,hisig,pltwlen1),10);
    xlabel(f1.ax((nrow1-1)*ncol1+1),sprintf('Time (sec) on %s %s, %s',mo,dy,yr),'FontSize',10);
    ylabel(f1.ax((nrow1-1)*ncol1+1),'Amplitude','FontSize',10);

    nsubtot1 = 0;  % count of the total subplots on the current figure f1
    
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
    
    %bursts and 4-s detections of the same day
    rangetemp = trange(trange(:,1)==datesets(i), :);
    hfdayi = hfbnd(hfbnd(:,daycol)==datesets(i), :);
    hfdayo = hfout(hfout(:,daycol)==datesets(i), :);
    
%     keyboard
    for j = 1: size(rangetemp,1)  
      k = k+1;  
      disp(k);

      tst = rangetemp(j,2); % start and end time of bursts
      ted = rangetemp(j,3);
      %it's guaranteed that if the requested time before data win is < t_max=15s used for grouping
      %burst, there will be no 4-s detections that occurred in the region of interest with that time,
      %otherwise it will be grouped into the same burst
      tbuffer = 10;   % buffer time to have a sense of noise level
      tstbuf = tst-tbuffer; % start and end time of bursts, buffer added
      tedbuf = ted+tbuffer;
      tlenbuf = tedbuf-tstbuf;
      
      nsub = tlenbuf/pltwlen1; % how mnay subplots needed for a single burst?
      if nsub < ceil(nsub)  % for most cases this is true
        tstbuf = tedbuf-ceil(nsub)*pltwlen1; % give the extra time to the start, pre-event is better for noise proxy 
      end
      nsub = ceil(nsub);  % the actual number of subplots needed for a single burst
      
      %align records
      msftadd = (round(max(abs([off12ran off13ran])))+1)*sps/40;  %+1 for safety
%       msftadd = round(sps/8);    % maximum allowed shift between 2 traces
      optcc = STAopt(max(floor((tst-1)*sps+1),1): min(floor((ted+1)*sps),86400*sps), 2:end);
      ccmid = ceil(size(optcc,1)/2);
      ccwlen = round(size(optcc,1)-2*(msftadd+1));
      loffmax = 4*sps/40;
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
      
      %Align records
      optdat = [];
      ortdat = [];
      optdat(:, 1:2) = STAopt(max(floor(tstbuf*sps+1),1): min(floor(tedbuf*sps),86400*sps), 1:2); % sta 1
      ortdat(:, 1:2) = STAort(max(floor(tstbuf*sps+1),1): min(floor(tedbuf*sps),86400*sps), 1:2);      
      optdat(:, 3) = STAopt(max(floor(tstbuf*sps+1)-off1i(k,2),1): ...
                            min(floor(tedbuf*sps)-off1i(k,2),86400*sps), 3); % sta 2
      ortdat(:, 3) = STAort(max(floor(tstbuf*sps+1)-off1i(k,2),1): ...
                            min(floor(tedbuf*sps)-off1i(k,2),86400*sps), 3);
      optdat(:, 4) = STAopt(max(floor(tstbuf*sps+1)-off1i(k,3),1): ...
                            min(floor(tedbuf*sps)-off1i(k,3),86400*sps), 4); % sta 3
      ortdat(:, 4) = STAort(max(floor(tstbuf*sps+1)-off1i(k,3),1): ...
                            min(floor(tedbuf*sps)-off1i(k,3),86400*sps), 4);
      
      tmaxi = hfdayi(:, seccol); % starting time of max power rate of half sec inside the ellipse
      tmaxo = hfdayo(:, seccol); % starting time of max power rate of half sec outside the ellipse
      tcnti = hfdayi(:, 15);  % the center of detecting win is the 15th col
      tbosti = bostdayi(:,5); % (arrival) time of Bostock's LFE catalog inside the rectangle
      tbosto = bostdayo(:,5); % (arrival) time of Bostock's LFE catalog outside the rectangle
      tarmi = armdayi(:,4); % (arrival) time of Armbruster's tremor catalog inside the rectangle
      tarmo = armdayo(:,4); % (arrival) time of Armbruster's tremor catalog outside the rectangle
      
      %how many 4-s detections fall into the burst range 
      indtmaxi = find(tmaxi>=tst-0.1 & tmaxi<=ted+0.1);
      ninbst(k,1) = length(indtmaxi); %+-0.1 s in case of resolution while saving to file
      
      %%%What is the max separation in time between these detections?
      %%%We need this, as we want to know how CC varies with this separation, rather than the data
      %%%window length, because data win could be long but temporally-clustered arrivals
      msept(k,1) = max(diff(tmaxi(indtmaxi)));
      
      %what is the max separation in time if use the center of the detecting window? this has to
      %larger than the step 1s of detecting win, overlapping unless >=4 s
      septcnt = diff(tcnti(indtmaxi));
      mseptcnt(k,1) = max(septcnt); 
      
      %what is the fraction of time of the data win that contains 4-s in-bound detections?
      tleni = 0;
      for jj = 1:ninbst(k,1)
        tleni = tleni+ range([max(tcnti(indtmaxi(jj))-2, tst-2) min(tcnti(indtmaxi(jj))+2, ted+2)]);
      end
      tlenovlp = 0; % length of time that is double-counted due to overlapping
      for jj = 1: length(septcnt)
        if septcnt(jj) < 4  % detecting win length, meanning there is overlapping
          tlenovlp = tlenovlp + 4-septcnt(jj);
        end
      end
      tleni = tleni - tlenovlp;  % subtract the double-counted overlapping length
      wlensec(k,1) = ted-tst+4;   % total length of data win in sec
      tfraci(k,1) = tleni / wlensec(k,1);
      
      %compute running CC between 3 stations
      [ircc,rcc12] = RunningCC(optdat(:,2), optdat(:,3), mwlen);
      [~,rcc13] = RunningCC(optdat(:,2), optdat(:,4), mwlen);
      [~,rcc23] = RunningCC(optdat(:,3), optdat(:,4), mwlen);
      rcc = (rcc12+rcc13+rcc23)/3;
      
      %compute the median absolute amplitude and envelope of the same moving window
      %for the moving window at the same station, sensable to use median
      [ir,ramp1,renv1] = Runningampenv(optdat(:,2),mwlen,mwlen-1,'median');
      [~,ramp2,renv2] = Runningampenv(optdat(:,3),mwlen,mwlen-1,'median');
      [~,ramp3,renv3] = Runningampenv(optdat(:,4),mwlen,mwlen-1,'median');
      %looks like using the amplitude and envelope are pretty similar
      %maybe better to mean, otherwise stations PGC will be underrepresented due to a larger amp
      %variation
      ramp = mean([ramp1 ramp2 ramp3], 2);  % use mean or median??
      renv = mean([renv1 renv2 renv3], 2);
%       ramp = median([ramp1 ramp2 ramp3], 2);  % use mean or median??
%       renv = median([renv1 renv2 renv3], 2);
            
      %target indices in terms of index of 'rcc', --> need -mwlen/2 
      ind = round((tst-2-tstbuf)*sps+1 -mwlen/2): round((ted+2-tstbuf)*sps -mwlen/2); % win +/-2 s, segment of interest
      
      %overall average max CC based on the current best alignment
      indopt = round((tst-tstbuf-2)*sps+1): round((ted-tstbuf+2)*sps);
      mcc12 = sum(optdat(indopt,2).*optdat(indopt,3))./ ...
        (sqrt(sum(optdat(indopt,2).^2)).*sqrt(sum(optdat(indopt,3).^2)));
      mcc13 = sum(optdat(indopt,2).*optdat(indopt,4))./ ...
        (sqrt(sum(optdat(indopt,2).^2)).*sqrt(sum(optdat(indopt,4).^2)));
      mcc23 = sum(optdat(indopt,3).*optdat(indopt,4))./ ...
        (sqrt(sum(optdat(indopt,3).^2)).*sqrt(sum(optdat(indopt,4).^2)));
      mcc(k,1) = (mcc12+mcc13+mcc23)/3;

      %store the result of the interested portion, normalize to its maximum
      runtmp = [runtmp; renv(ind)/max(renv(ind)) ramp(ind)/max(ramp(ind)) rcc(ind)];
      
      %for all the rest segment excluding win +/-2 s
      dumrcc = rcc; dumrcc(ind)=[];   
      dumrenv = renv; dumrenv(ind)=[];
      dumramp = ramp; dumramp(ind)=[];
      dumir = ir; dumir(ind)=[];
      
      %As explained above, the 4-s win below definitely doesn't have 4-s detections, in terms of index of 'rcc'
      ind1 = round((tst-7-tstbuf)*sps+1 -mwlen/2): round((tst-3-tstbuf)*sps -mwlen/2);  % 4-s prior to win +/-2 s, 1-s buffer
      
      rcomp(k,1:4) = [median(renv(ind)), median(renv(ind1)), median(rcc(ind)), median(rcc(ind1))];
      if median(renv(ind)) >= median(renv(ind1))
        rcomp(k,5) = k;
      else
        rcomp(k,5) = -k;
      end
      if median(rcc(ind)) >= median(rcc(ind1))
        rcomp(k,6) = k;
      else
        rcomp(k,6) = -k;
      end
      
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
      indnodet = setdiff(setdiff(ind, inddeti), inddeto);

      
      nsubtot2 = nsubtot2+1;  % count of the total subplots on the current figure f2
      if nsubtot2 > nrow2*ncol2
          %print the current figure f2 
          orient(f2.fig,'landscape');
          f2name{ifig2-100,1} = sprintf('rccvsrenv002_%.1fsFig%s_%.1f-%.1f.pdf',...
            mwlen/sps,num2zeropadstr(ifig2-100,2),losig,hisig);
          print(f2.fig,'-dpdf','-fillpage',fullfile(rstpath, '/FIGS',f2name{ifig2-100,1}));
          %move on to next figure f2
          ifig2 = ifig2+1;
          f2 = initfig(widin2,htin2,nrow2,ncol2,ifig2);
          axpos2 = optaxpos(f2,nrow2,ncol2,pltxran2,pltyran2,pltxsep2,pltysep2);
          axtit = [f2.ax(1) f2.ax(2) f2.ax(3) f2.ax(4) f2.ax(5)];
          supertit(axtit, sprintf('Fig %s, Passband: %.1f-%.1f Hz, cclen: %.1f s, rcc vs. renv',...
            num2zeropadstr(ifig2-100,2),losig,hisig,mwlen/sps),10);
          xlabel(f2.ax((nrow2-1)*ncol2+1),'Running envelope','FontSize',10);
          ylabel(f2.ax((nrow2-1)*ncol2+1),'Running CC','FontSize',10);
          
          %print the current figure f3, same as to f2 
          orient(f3.fig,'landscape');
          f3name{ifig3-200,1} = sprintf('renvCDF002_%.1fsFig%s_%.1f-%.1f.pdf',...
            mwlen/sps,num2zeropadstr(ifig3-200,2),losig,hisig);
          print(f3.fig,'-dpdf','-fillpage',fullfile(rstpath, '/FIGS',f3name{ifig3-200,1}));
          %move on to next figure f2
          ifig3 = ifig3+1;
          f3 = initfig(widin2,htin2,nrow2,ncol2,ifig3);
          axpos3 = optaxpos(f3,nrow2,ncol2,pltxran2,pltyran2,pltxsep2,pltysep2);
          axtit = [f3.ax(1) f3.ax(2) f3.ax(3) f3.ax(4) f3.ax(5)];
          supertit(axtit, sprintf('Fig %s, Passband: %.1f-%.1f Hz, cclen: %.1f s, renv CDF',...
            num2zeropadstr(ifig3-200,2),losig,hisig,mwlen/sps),10);
          xlabel(f3.ax((nrow2-1)*ncol2+1),'Running envelope','FontSize',10);
          ylabel(f3.ax((nrow2-1)*ncol2+1),'Empirical CDF','FontSize',10);

          nsubtot2 = nsubtot2-nrow2*ncol2;  %refresh
      end

      %%%plot elements for f2, number: 1??
      ax = f2.ax(nsubtot2); hold(ax,'on');      
%       scatter(ax,dumrenv,dumrcc,6,[.6 .6 .6],'filled'); % segment excluding win +/-2 s
      scatter(ax,renv(ind1),rcc(ind1),6,[.6 .6 .6],'filled'); % 4-s prior
      scatter(ax,renv(indnodet),rcc(indnodet),6,'c','filled'); % no detections
      scatter(ax,renv(inddeto),rcc(inddeto),6,'r','filled'); % out-bound 4-s detections
      scatter(ax,renv(inddeti),rcc(inddeti),6,'k','filled'); % in-bound 4-s detections
%       scatter(ax,renv(ind),rcc(ind),10,ir(ind)/sps,'filled'); % win +/-2 s
%       scatter(ax,median(dumrenv),median(dumrcc),8,'ko','linew',1);
      scatter(ax,median(renv(ind1)),median(rcc(ind1)),15,'ko','linew',1,'MarkerFaceColor',...
        [.6 .6 .6]);
      scatter(ax,median(renv(ind)),median(rcc(ind)),30,'ko','linew',1.5);
%       scatter(ax,dumramp,dumrcc/yran(2),6,[.6 .6 .6],'filled'); hold on
%       scatter(ax,ramp(ind),rcc(ind)/yran(2),10,ir(ind),'filled');
%       [xcnt,ycnt,y1sig] = ranybinx(renv(ind),rcc(ind),'median',10);
%       errorbar(ax,xcnt,ycnt,-y1sig,y1sig,'vertical','o','markersize',3,'color',...
%          'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor','k','CapSize',4);
      text(ax,0.98,0.35,sprintf('%d',date),'unit','normalized',...
        'HorizontalAlignment','right','fontsize',8.5); % ETS
      text(ax,0.98,0.3,sprintf('Win %d',k),'unit','normalized',...
        'HorizontalAlignment','right','fontsize',8.5); % burst number
      text(ax,0.98,0.25,sprintf('(-2)%.1f - %.1f(+2) s',tst-tstbuf,ted-tstbuf),'unit','normalized',...
        'HorizontalAlignment','right','fontsize',8.5);  %start and end of window
      [spear(k,1),~] = corr(renv(ind),rcc(ind),'Type','Spearman'); % Spearman coeff
      [ken(k,1),~] = corr(renv(ind),rcc(ind),'Type','Kendall');  % Kendall coeff
      text(ax,0.98,0.55,sprintf('S: %.2f',spear(k,1)),'unit','normalized',...
        'HorizontalAlignment','right','fontsize',8);
      text(ax,0.98,0.50,sprintf('K: %.2f',ken(k,1)),'unit','normalized',...
        'HorizontalAlignment','right','fontsize',8);
      colormap(ax,'jet');
      caxis(ax,[0 size(optdat,1)/sps]);
      c=colorbar(ax,'south');
      c.Position(2) = c.Position(2)-0.005;
%       c.Position(3) = axpos2(1,3)*0.75;
      c.Position(4) = 0.01;
      c.Label.String = strcat({'Time (s) since the start'});
      c.Label.FontSize = 8;
      ylim(ax,[-1 1]);
      longticks(ax,2);
      hold(ax,'off');
      
      %%%plot elements for f3, number: 2??
      ax = f3.ax(nsubtot2); hold(ax,'on');
%       [f,x] = ecdf(dumrenv);  % segment excluding win +/-2 s
      [f,x] = ecdf(renv(ind1)); % 4-s prior to segment of interest 
      plot(ax,x,f,'color',[.6 .6 .6],'linew',1);
      ind2 = round((tst-tstbuf-5)*sps+1): round((ted-tstbuf+5)*sps);  % win +/-5 s
      [f,x] = ecdf(renv(ind2));
      plot(ax,x,f,'color','c','linew',1);      
      [f,x] = ecdf(renv(ind));  % win +/-2 s, segment of interest
      plot(ax,x,f,'color','k','linew',2);  
%       pctl = [0.13 2.28 15.87 50 84.13 97.72 99.87];  % -3\sigma: \sigma: 3\sigma
      pctl = [2.28 50]; % -2\sigma, 0(median)
      pctlv = prctile(dumrenv,pctl);
      text(ax,0.98,0.7,sprintf('%.3f; %.3f',pctlv(1),pctlv(2)),'unit','normalized',...
        'HorizontalAlignment','right','fontsize',8.5,'color',[.6 .6 .6]); 
      pctlv = prctile(renv(ind2),pctl);
      text(ax,0.98,0.65,sprintf('%.3f; %.3f',pctlv(1),pctlv(2)),'unit','normalized',...
        'HorizontalAlignment','right','fontsize',8.5,'color','c'); 
      pctlv = prctile(renv(ind),pctl);
      text(ax,0.98,0.60,sprintf('%.3f; %.3f',pctlv(1),pctlv(2)),'unit','normalized',...
        'HorizontalAlignment','right','fontsize',8.5,'color','k'); 
      plot(ax,ax.XLim, [pctl(1)/100 pctl(1)/100], 'k--');
      plot(ax,ax.XLim, [pctl(2)/100 pctl(2)/100], 'k--');
      text(ax,0.98,0.35,sprintf('%d',date),'unit','normalized',...
        'HorizontalAlignment','right','fontsize',8.5); % ETS
      text(ax,0.98,0.3,sprintf('Win %d',k),'unit','normalized',...
        'HorizontalAlignment','right','fontsize',8.5); % burst number
      text(ax,0.98,0.25,sprintf('(-2)%.1f - %.1f(+2) s',tst-tstbuf,ted-tstbuf),'unit','normalized',...
        'HorizontalAlignment','right','fontsize',8.5);  %start and end of window
      ylim(ax,[0 1]);
      longticks(ax,2);
      hold(ax,'off');

      %scale seismogram for nicer plotting
      sclseis = 1.4;
      yran = sclseis*[-max(max(abs(optdat(:,2:end)))) max(max(abs(optdat(:,2:end))))];  % left y axis range
      
      for jj = 1: nsub
        isub = nsubtot1+jj;
        if isub > nrow1*ncol1
          
          %print the current figure f1
          orient(f1.fig,'landscape');
          f1name{ifig1,1} = sprintf('seisbst002_%dsFig%s_%d_%.1f-%.1f.pdf',...
            pltwlen1,num2zeropadstr(ifig1,2),date,losig,hisig);
          print(f1.fig,'-dpdf','-fillpage',fullfile(rstpath, '/FIGS',f1name{ifig1,1}));
          %move on to next figure f1
          ifig1 = ifig1+1;  %initialize another figure
          f1 = initfig(widin1,htin1,nrow1,ncol1,ifig1);
          axpos1 = optaxpos(f1,nrow1,ncol1,pltxran1,pltyran1,pltxsep1,pltysep1);
          axtit = [f1.ax(1) f1.ax(2)];
          supertit(axtit, sprintf('Fig %s, %d, %s %s, %s, Passband: %.1f-%.1f Hz, wlen: %d s',...
            num2zeropadstr(ifig1,2),date,mo,dy,yr,losig,hisig,pltwlen1),10);
          xlabel(f1.ax(nrow1*ncol1-1),sprintf('Time (s) on %s %s, %s',mo,dy,yr),'FontSize',10);
          ylabel(f1.ax(nrow1*ncol1-1),'Amplitude','FontSize',10);
          
          nsubtot1 = nsubtot1-nrow1*ncol1;  %refresh
          isub = nsubtot1+jj;
        end

        %%%plot elements for f1, number: ??
        ax = f1.ax(isub); hold(ax,'on');
%         yyaxis(ax,'left');
        plot(ax,optdat(:,1),optdat(:,2),'r-','linew',0.5);
        plot(ax,optdat(:,1),optdat(:,3),'b-','linew',0.5);
        plot(ax,optdat(:,1),optdat(:,4),'k-','linew',0.5); %,'linew',0.5
        xran = [tstbuf+(jj-1)*pltwlen1 tstbuf+jj*pltwlen1]; % x axis range in sec
        xlim(ax,xran); 
        ylim(ax,yran);
        longticks(ax,5);
        %plot the strongest 0.5-s arrival outside ellipse
        ind = find(tmaxo>=xran(1) & tmaxo<=xran(2));        
        for ii = 1: length(ind)
          barst = tmaxo(ind(ii)); % arrival of strongest 0.5-s, in sec 
          bared = tmaxo(ind(ii))+0.5; % last for 0.5 s
          yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.97;
          plot(ax,[barst bared],[yloc yloc],'-','linew',2.5,'color',[.6 .6 .6]);
        end
        %plot the strongest 0.5-s arrival inside ellipse
        ind = find(tmaxi>=xran(1) & tmaxi<=xran(2));        
        for ii = 1: length(ind)
          barst = tmaxi(ind(ii)); % arrival of strongest 0.5-s, in sec 
          bared = tmaxi(ind(ii))+0.5; % last for 0.5 s
          yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.97;
          plot(ax,[barst bared],[yloc yloc],'k-','linew',2.5);
        end
        %plot the armbruster's tremor catalog outside rectangle
        ind = find(tarmo>=xran(1) & tarmo<=xran(2));
        yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.97;
        scatter(ax,tarmo(ind), yloc*ones(size(tarmo(ind))),10,'g','linew',1); % armbruster
        %plot the armbruster's tremor catalog inside rectangle
        ind = find(tarmi>=xran(1) & tarmi<=xran(2));
        scatter(ax,tarmi(ind), yloc*ones(size(tarmi(ind))),10,'g','filled');
        %plot the bostock's LFE catalog outside rectangle
        ind = find(tbosto>=xran(1) & tbosto<=xran(2));
        scatter(ax,tbosto(ind), yloc*ones(size(tbosto(ind))),10,'m','linew',1); % bostock
        %plot the bostock's LFE catalog inside rectangle
        ind = find(tbosti>=xran(1) & tbosti<=xran(2));
        scatter(ax,tbosti(ind), yloc*ones(size(tbosti(ind))),10,'m','filled');
        yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.1;
        if tst>=xran(1) && tst<=xran(2)
          plot(ax,[tst tst],ax.YLim,'--','Color',[.5 .5 .5],'linew',1);
          text(ax,tst,yloc,'start','fontsize',8);
        end
        if ted>=xran(1) && ted<=xran(2)
          plot(ax,[ted ted],ax.YLim,'--','Color',[.5 .5 .5],'linew',1);
          text(ax,ted,yloc,'end','fontsize',8);
        end
        text(ax,0.01,0.1,num2str(off1i(k,2)),'unit','normalized',...
          'HorizontalAlignment','left','fontsize',8,'color','b'); % offset 12
        text(ax,0.04,0.1,num2str(off1i(k,3)),'unit','normalized',...
          'HorizontalAlignment','left','fontsize',8,'color','k'); % offset 13
        text(ax,0.07,0.1,sprintf('%.3f',ccali(k)),'unit','normalized',...
          'HorizontalAlignment','left','fontsize',8,'color','k'); % offset 13
        text(ax,0.01,0.85,sprintf('Win%d--N%d--P%d/%d',k,ninbst(k,1),jj,nsub),'unit','normalized',...
          'HorizontalAlignment','left','fontsize',8); % burst number
%         plot(ax,ircc/sps+tstbuf,rcc*yran(2),'co','markersize',0.5);  % scale it to the amplitude range
%         plot(ax,ir/sps+tstbuf,renv*1.2+0.3*yran(2),'Color',[.6 .6 .6],'linew',0.5);  % scale it to the amplitude range
        %%%%%%%%%%% Below is colorcode rcc by renv
        sclrcc = 1.2;
        scatter(ax,ircc/sps+tstbuf,rcc*yran(2)/sclseis*sclrcc,1,renv);  % scale it to the amplitude range
        colormap(ax,'jet');
        caxis(ax,[min(renv) max(renv)]);
        if jj==nsub
          c=colorbar(ax,'south');
          c.Position(3) = axpos1(isub,3)*0.25;
          c.Position(1) = (axpos1(isub,3)+axpos1(isub,1))-c.Position(3)-0.012;
          c.Position(4) = 0.005;
          c.Ticks = [];
          c.Label.FontSize = 6;
          text(ax,0.72,0.2,sprintf('%.3f',min(renv)),...
            'HorizontalAlignment','center','FontSize',8,'unit','normalized');
          text(ax,0.96,0.2,sprintf('%.3f',max(renv)),...
            'HorizontalAlignment','center','FontSize',8,'unit','normalized');
          text(ax,0.84,0.2,'renv',...
            'HorizontalAlignment','center','FontSize',8,'unit','normalized');
        end
        %%%%%%%%%%% Above is colorcode rcc by renv
%         %%%%%%%%%%% Below is colorcode renv by rcc
%         sclrenv = 2;
%         yzero = 0.2*yran(2);
%         plot(ax,[xran(1) xran(1)+0.5],[yzero yzero],'k-','linew',1);
%         scatter(ax,ir/sps+tstbuf,renv*sclrenv+yzero,2,rcc);
%         colormap(ax,'jet');
%         caxis(ax,[min(rcc) max(rcc)]);
%         if jj==nsub
%           c=colorbar(ax,'south');
%           c.Position(3) = axpos1(isub,3)*0.25;
%           c.Position(1) = (axpos1(isub,3)+axpos1(isub,1))-c.Position(3)-0.012;
%           c.Position(4) = 0.005;
%           c.Ticks = [];
%           c.Label.FontSize = 6;
%           text(ax,0.72,0.2,sprintf('%.3f',min(rcc)),...
%             'HorizontalAlignment','center','FontSize',8,'unit','normalized');
%           text(ax,0.96,0.2,sprintf('%.3f',max(rcc)),...
%             'HorizontalAlignment','center','FontSize',8,'unit','normalized');
%           text(ax,0.84,0.2,'rcc',...
%             'HorizontalAlignment','center','FontSize',8,'unit','normalized');
%           yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.3;
%           xloc = ax.XLim(1)+(ax.XLim(2)-ax.XLim(1))*0.72;
%           plot(ax,[xloc xloc+0.5],[yloc yloc],'-','linew',2,'color','c');
%           xloc = ax.XLim(1)+(ax.XLim(2)-ax.XLim(1))*0.77;
%           plot(ax,[xloc xloc+0.5*max(renv)/min(renv)],[yloc yloc],'-','linew',2,'color','c');
%         end
%         %%%%%%%%%%% Above is colorcode renv by rcc
        
%         [irr,rcorr,~] = Runningcorr(renv,rcc,mwlen,mwlen/2,'Spearman');
%         plot(ax,(irr+mwlen/2)/sps+tstbuf,(rcorr*0.5+0.5)*yran(2),'Color',[.6 .6 .6],'linew',0.5);
        hold(ax,'off');
      end
      
      nsubtot1 = nsubtot1 + nsub; % counting number subplots for each figure
      
    end
    
    %if the same date is finished, print it
    %print the current figure f1
    orient(f1.fig,'landscape');
    f1name{ifig1,1} = sprintf('seisbst002_%dsFig%s_%d_%.1f-%.1f.pdf',...
      pltwlen1,num2zeropadstr(ifig1,2),date,losig,hisig);
    print(f1.fig,'-dpdf','-fillpage',fullfile(rstpath, '/FIGS',f1name{ifig1,1}));
  end
  
  %store the result of the interested portion
  runall{iets} = runtmp;
  
  %if this is already the last burst, print it
  if  k==size(trange,1)
    %print the current figure f1
    orient(f1.fig,'landscape');
    f1name{ifig1,1} = sprintf('seisbst002_%dsFig%s_%d_%.1f-%.1f.pdf',...
      pltwlen1,num2zeropadstr(ifig1,2),date,losig,hisig);
    print(f1.fig,'-dpdf','-fillpage',fullfile(rstpath, '/FIGS',f1name{ifig1,1}));
  end

  %if this is already the last burst, print it
  if k==size(trange,1)
    %print the current figure f2
    orient(f2.fig,'landscape');
    f2name{ifig2-100,1} = sprintf('rccvsrenv002_%.1fsFig%s_%.1f-%.1f.pdf',...
      mwlen/sps,num2zeropadstr(ifig2-100,2),losig,hisig);
    print(f2.fig,'-dpdf','-fillpage',fullfile(rstpath, '/FIGS',f2name{ifig2-100,1}));
    %print the current figure f3, same as to f2
    orient(f3.fig,'landscape');
    f3name{ifig3-200,1} = sprintf('renvCDF002_%.1fsFig%s_%.1f-%.1f.pdf',...
      mwlen/sps,num2zeropadstr(ifig3-200,2),losig,hisig);
    print(f3.fig,'-dpdf','-fillpage',fullfile(rstpath, '/FIGS',f3name{ifig3-200,1}));

  end


%   close all
end
%save all workspace variables
save('seisbursts_allvari.mat');

else
  load('seisbursts_allvari.mat');
end
keyboard

%% A few more plots
% %%%Spearman vs. data win length
% figure; 
% scatter(wlensec,spear,20,1:length(wlensec),'filled');
% colorbar; colormap jet; box on; grid on;
% xlabel('Length of selected window (s)');
% ylabel("Spearman's coeff");

%%%Median RCC or overall CC of the data win vs. data win length
figure
ax=subplot(121);
scatter(ax,wlensec,rcomp(:,3),20,1:length(wlensec),'filled');hold(ax,'on');
colorbar; colormap jet; box on; grid on;
[speartmp,~] = corr(wlensec(:),rcomp(:,3),'Type','Spearman'); % Spearman coeff
text(ax,0.95,0.9,sprintf('S: %.2f',speartmp),'unit','normalized',...
  'HorizontalAlignment','right','fontsize',10);
xlabel(ax,'Length of selected window (s)');
ylabel(ax,'Median RCC of selected window');
hold(ax,'off');
ax=subplot(122);
scatter(ax,wlensec,mcc,20,1:length(wlensec),'filled');hold(ax,'on');
[speartmp,~] = corr(wlensec(:),mcc(:),'Type','Spearman'); % Spearman coeff
text(ax,0.95,0.9,sprintf('S: %.2f',speartmp),'unit','normalized',...
  'HorizontalAlignment','right','fontsize',10);
colorbar; colormap jet; box on; grid on;
xlabel(ax,'Length of selected window (s)');
ylabel(ax,'Overall CC selected window');
hold(ax,'off');

%%%Median RCC or overall CC of the data win vs. max. separation between strongest 0.5-s arrivals
figure
ax=subplot(221);
scatter(ax,msept,rcomp(:,3),20,1:length(msept),'filled'); hold(ax,'on');
colorbar; colormap jet; box on; grid on;
[speartmp,~] = corr(msept(:),rcomp(:,3),'Type','Spearman'); % Spearman coeff
text(ax,0.95,0.9,sprintf('S: %.2f',speartmp),'unit','normalized',...
  'HorizontalAlignment','right','fontsize',10);
% create fit object with free constraints
fttpfree = fittype( @(a,b,x) a*x+b);
[fitobj,gof,outp] = fit(msept(:),rcomp(:,3),fttpfree,'Robust','Bisquare',...
  'StartPoint',[1 1]);
coef = coeffvalues(fitobj);
coefci = confint(fitobj);
xplt = min(msept(:)): 0.01: max(msept(:));
yfit = coef(1)*xplt + coef(2);
ylow = coefci(1,1)*xplt + coefci(1,2);
yupp = coefci(2,1)*xplt + coefci(2,2);
plot(ax,xplt,yfit,'k-','linew',2);
patch(ax,[xplt fliplr(xplt)],[ylow fliplr(yupp)],'k','Facealpha',0.2,'edgecolor','none');
xlabel(ax,'Max. separation between strongest 0.5-s arrivals (s)');
ylabel(ax,'Median RCC of selected window');
hold(ax,'off');

ax=subplot(222);
scatter(ax,msept,mcc,20,1:length(msept),'filled'); hold(ax,'on');
colorbar; colormap jet; box on; grid on;
[speartmp,~] = corr(msept(:),mcc(:),'Type','Spearman'); % Spearman coeff
text(ax,0.95,0.9,sprintf('S: %.2f',speartmp),'unit','normalized',...
  'HorizontalAlignment','right','fontsize',10);
[fitobj,~,~] = fit(msept(:),mcc(:),fttpfree,'Robust','Bisquare',...
  'StartPoint',[1 1]);
coef = coeffvalues(fitobj);
coefci = confint(fitobj);
yfit = coef(1)*xplt + coef(2);
ylow = coefci(1,1)*xplt + coefci(1,2);
yupp = coefci(2,1)*xplt + coefci(2,2);
plot(ax,xplt,yfit,'k-','linew',2);
patch(ax,[xplt fliplr(xplt)],[ylow fliplr(yupp)],'k','Facealpha',0.2,'edgecolor','none');
xlabel(ax,'Max. separation between strongest 0.5-s arrivals (s)');
ylabel(ax,'Overall CC selected window');
hold(ax,'off');

ax=subplot(223);
scatter(ax,mseptcnt,rcomp(:,3),20,1:length(mseptcnt),'filled'); hold(ax,'on');
colorbar; colormap jet; box on; grid on;
[speartmp,~] = corr(mseptcnt(:),rcomp(:,3),'Type','Spearman'); % Spearman coeff
text(ax,0.95,0.9,sprintf('S: %.2f',speartmp),'unit','normalized',...
  'HorizontalAlignment','right','fontsize',10);
[fitobj,~,~] = fit(mseptcnt(:),rcomp(:,3),fttpfree,'Robust','Bisquare',...
  'StartPoint',[1 1]);
coef = coeffvalues(fitobj);
coefci = confint(fitobj);
xplt = min(mseptcnt(:)): 0.01: max(mseptcnt(:));
yfit = coef(1)*xplt + coef(2);
ylow = coefci(1,1)*xplt + coefci(1,2);
yupp = coefci(2,1)*xplt + coefci(2,2);
plot(ax,xplt,yfit,'k-','linew',2);
patch(ax,[xplt fliplr(xplt)],[ylow fliplr(yupp)],'k','Facealpha',0.2,'edgecolor','none');
xlabel(ax,'Max. separation between 4-s detecting windows (s)');
ylabel(ax,'Median RCC of selected window');
hold(ax,'off');

ax=subplot(224);
scatter(ax,mseptcnt,mcc,20,1:length(mseptcnt),'filled'); hold(ax,'on');
colorbar; colormap jet; box on; grid on;
[speartmp,~] = corr(mseptcnt(:),mcc(:),'Type','Spearman'); % Spearman coeff
text(ax,0.95,0.9,sprintf('S: %.2f',speartmp),'unit','normalized',...
  'HorizontalAlignment','right','fontsize',10);
[fitobj,~,~] = fit(mseptcnt(:),mcc(:),fttpfree,'Robust','Bisquare',...
  'StartPoint',[1 1]);
coef = coeffvalues(fitobj);
coefci = confint(fitobj);
yfit = coef(1)*xplt + coef(2);
ylow = coefci(1,1)*xplt + coefci(1,2);
yupp = coefci(2,1)*xplt + coefci(2,2);
plot(ax,xplt,yfit,'k-','linew',2);
patch(ax,[xplt fliplr(xplt)],[ylow fliplr(yupp)],'k','Facealpha',0.2,'edgecolor','none');
xlabel(ax,'Max. separation between 4-s detecting windows (s)');
ylabel(ax,'Overall CC selected window');
hold(ax,'off');

%%%Median RCC or overall CC of the data win vs. fraction of time that has 4-s in-bound detections
figure
ax=subplot(121);
scatter(ax,tfraci,rcomp(:,3),20,1:length(tfraci),'filled'); hold(ax,'on');
colorbar; colormap jet; box on; grid on;
[speartmp,~] = corr(tfraci(:),rcomp(:,3),'Type','Spearman'); % Spearman coeff
text(ax,0.95,0.9,sprintf('S: %.2f',speartmp),'unit','normalized',...
  'HorizontalAlignment','right','fontsize',10);
[fitobj,~,~] = fit(tfraci(:),rcomp(:,3),fttpfree,'Robust','Bisquare',...
  'StartPoint',[1 1]);
coef = coeffvalues(fitobj);
coefci = confint(fitobj);
xplt = min(tfraci(:)): 0.01: max(tfraci(:));
yfit = coef(1)*xplt + coef(2);
ylow = coefci(1,1)*xplt + coefci(1,2);
yupp = coefci(2,1)*xplt + coefci(2,2);
plot(ax,xplt,yfit,'k-','linew',2);
patch(ax,[xplt fliplr(xplt)],[ylow fliplr(yupp)],'k','Facealpha',0.2,'edgecolor','none');
xlabel(ax,'Fraction of time with 4-s in-bound detections');
ylabel(ax,'Median RCC of selected window');
hold(ax,'off');

ax=subplot(122);
scatter(ax,tfraci,mcc,20,1:length(tfraci),'filled'); hold(ax,'on');
colorbar; colormap jet; box on; grid on;
[speartmp,~] = corr(tfraci(:),mcc(:),'Type','Spearman'); % Spearman coeff
text(ax,0.95,0.9,sprintf('S: %.2f',speartmp),'unit','normalized',...
  'HorizontalAlignment','right','fontsize',10);
[fitobj,~,~] = fit(tfraci(:),mcc(:),fttpfree,'Robust','Bisquare',...
  'StartPoint',[1 1]);
coef = coeffvalues(fitobj);
coefci = confint(fitobj);
yfit = coef(1)*xplt + coef(2);
ylow = coefci(1,1)*xplt + coefci(1,2);
yupp = coefci(2,1)*xplt + coefci(2,2);
plot(ax,xplt,yfit,'k-','linew',2);
patch(ax,[xplt fliplr(xplt)],[ylow fliplr(yupp)],'k','Facealpha',0.2,'edgecolor','none');
xlabel(ax,'Fraction of time with 4-s in-bound detections');
ylabel(ax,'Overall CC selected window');
hold(ax,'off');


%% summary plot of rc VS renv for all bursts
pctl = [2.28 50];
dnpts=round(mwlen/4);
for iets = 1: nets
  runtmp = runall{iets};
  runplt = runtmp(1:dnpts:end,:);
  nbin=10;
  [xbin,indbin] = binxeqnum(runplt(:,1),nbin);  %bin by env with same number
  for i = 1: nbin
    indi = indbin{i};
    xcnt(i,iets) = median(xbin{i});
    ycnt(i,iets) = median(runplt(indi,3));
    y1sig(i,iets) = std(runplt(indi,3));  
  end
%   [xcnt(:,iets),ycnt(:,iets),y1sig(:,iets)] = ranybinx(runplt(:,1),runplt(:,3),'median',10,[],[]);
  [spearets(iets),~] = corr(runplt(:,1),runplt(:,3),'Type','Spearman');
  [kenets(iets),~] = corr(runplt(:,1),runplt(:,3),'Type','Kendall');
%   pctlv(1:2,iets) = prctile(runplt(:,1),pctl);
end

%%
nrow = 1; % rows and cols of subplots in each figure
ncol = 3; 
widin = 8.3; % size of each figure
htin = 4;
pltxran = [0.07 0.98]; pltyran = [0.18 0.88];
pltxsep = 0.02; pltysep = 0.03;
f = initfig(widin,htin,nrow,ncol);
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
axtit = [f.ax(1) f.ax(2) f.ax(3)];
% supertit(axtit, sprintf('Fig %d, Passband: %.1f-%.1f Hz, cclen: %.1f s, rcc vs. renv',...
%   ifig2,losig,hisig,cclen/sps),10);
xlabel(f.ax(1),'Normalized running envelope','FontSize',10);
ylabel(f.ax(1),'Running CC','FontSize',10);

for iets = 1: nets
  ax = f.ax(iets); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  runtmp = runall{iets};
  runplt = runtmp(1:dnpts:end,:);
  
  scale='log10';
  binmethod='grid';
  xran=[0 1]; yran=[-1 1];
  dx=0.01; dy=0.01;
  den1d = density_matrix(runplt(:,1),runplt(:,3),xran,yran,dx,dy);
  den1d = den1d(den1d(:,3)>0, :);
  den1d = sortrows(den1d,3);
  
  dum = den1d;
  dum(dum(:,3)>1, :) = [];
  if strcmp(scale,'log10')
    dum(:,3) = log10(dum(:,3));
  end
  scatter(ax,dum(:,1),dum(:,2),3,dum(:,3),'o','linew',0.2);  %, 'MarkerEdgeColor', 'w')
  
  dum = sortrows(den1d,3);
  dum(dum(:,3)==1, :) = [];
  if strcmp(scale,'log10')
    dum(:,3) = log10(dum(:,3));
  end  
  scatter(ax,dum(:,1),dum(:,2),3,dum(:,3),'o','filled','MarkerEdgeColor','none');
%   colormap(ax,flipud(colormap(ax,'kelicol')));
  % colormap(ax,jet);
  colormap(ax,'plasma');
  c=colorbar(ax,'SouthOutside');
  pos = ax.Position;
  c.Position = [pos(1), pos(2)-0.05, pos(3), 0.03];
  if strcmp(scale,'log10')
    cstr = strcat({'log_{10}(# points / '},binmethod,')');
  elseif strcmp(scale,'linear')
    cstr = strcat({'# points / '},binmethod);  
  end
  c.Label.String = cstr;
  
%   scatter(ax,runplt(:,1),runplt(:,3),6,'MarkerFaceColor',[.2 .2 .2],'MarkerEdgeColor',...
%     'none','MarkerFaceAlpha',.2); hold on
  scatter(ax,median(runplt(:,1)),median(runplt(:,3)),40,'k^','filled');
%   text(ax,median(runplt(:,1))*1.05,median(runplt(:,3))*1.05,sprintf('(%.3f, %.3f)',...
%     median(runplt(:,1)),median(runplt(:,3))),'HorizontalAlignment',...
%     'left','fontsize',10);
  text(ax,0.98,0.05,sprintf('%d',years(iets)),'unit','normalized',...
    'HorizontalAlignment','right','fontsize',10);
%   text(ax,0.2,0.1,sprintf('2.5 prctile: %.3f',prctile(runtmp(:,1),2.5)),'unit','normalized',...
%     'HorizontalAlignment','left','fontsize',10);
  errorbar(ax,xcnt(:,iets),ycnt(:,iets),-y1sig(:,iets),y1sig(:,iets),'vertical','o',...
    'markersize',4,'color','k','linewidth',0.8,'MarkerEdgeColor','k',...
    'MarkerFaceColor','none','CapSize',4);
  text(ax,0.02,0.1,sprintf('S: %.2f',spearets(iets)),'unit','normalized',...
    'HorizontalAlignment','left','fontsize',8);
  text(ax,0.02,0.05,sprintf('K: %.2f',kenets(iets)),'unit','normalized',...
    'HorizontalAlignment','left','fontsize',8);
%   text(ax,0.02,0.05,sprintf('%.3f; %.3f',pctlv(1,iets),pctlv(2,iets)),'unit','normalized',...
%     'HorizontalAlignment','left','fontsize',8.5,'color','k');
  ax.GridLineStyle = '--';
  ax.XAxisLocation = 'top';
  xlim(ax,xran);
  ylim(ax,yran);
  longticks(ax,2);
  xticks(ax,xran(1):0.2:xran(2));
  yticks(ax,yran(1):0.2:yran(2));
  
  if iets ~= 1
    nolabels(ax,2);
  end
  hold(ax,'off');

end

% orient(f.fig,'landscape');
fname = 'rccvsrenv.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

% clear runall runtmp

keyboard

%% merge all figures into a single pdf file
status = system('rm -f /home/data2/chaosong/matlab/allan/data-no-resp/PGCtrio/FIGS/aaa.pdf');
for i = 1:size(f1name,1)
  fname{i} = fullfile(rstpath, '/FIGS/',f1name{i});
end
append_pdfs(fullfile(rstpath, '/FIGS/aaa.pdf'),fname);


