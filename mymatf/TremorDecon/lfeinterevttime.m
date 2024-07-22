% lfeinterevttime.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the script to obtain the inter-event time of the deconvolved LFE
% catalog, either 3-station or 4-station catalog, and compare with noise.
% A further step is to first bin them by amplitude.
% --Instead of using the direct amp of the source itself, we try to link the
% the detection with seismogram. But it is not the instant seismic amp at
% the same time, but the median amp of the time window which is centered at
% the detection. The length of the time window should be chosen carefully.
% --In addition to the amplitude which can be positive or negative although
% the median may be more stable, the envelope might be better.
% --Also, note that some detections are very close in time, therefore the 
% windows chosen would be highly overlapping, thus the median amplitude can 
% be very simialr.
% --The final amplitude is the mean of the two. 
%
% --The whole analysis is another way to address whether is the saturation 
% is related to the amplitude. In 'N2Nmstat_data_ref.m', we obtain the 
% fraction of diff time that is within 0.25*m+0.125 s, the fraction of 
% unique events, after binning these events by their median amp. Here, the 
% direct waveform amp is used. 
% 
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/10/23
% Last modified date:   2023/10/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
% format short e   % Set the format to 5-digit floating point
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

cutout = 'ellipse';
ttol = 35;
ntol = 3;
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.pgc002.',cutout(1:4)));
% trange = trange(1:end-1,:);
nbst = size(trange,1);
tlen = trange(:,3)-trange(:,2);
tlensum = sum(tlen);

modname = 'timeoff_plfit_4thsta_160sps.mat';
planefit = load(strcat(rstpath, '/MAPS/',modname));
rmse = planefit.gof{4}.rmse;
offmax = round(2.0*rmse);

%%%load data
% savefile = 'deconv_stats4th_allbstsig.mat';
savefile = 'deconv_stats4th_no23_allbstsig.mat';
load(strcat(rstpath, '/MAPS/',savefile));
% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
load(strcat(rstpath, '/MAPS/',savefile));
 
% keyboard

%%
%%%param for secondary sources removed
locxyprojall = allbstsig.locxyprojall;
tarvlsplstall = allbstsig.impindepall(:,1);
nsrc = allbstsig.nsrc;
imp = allbstsig.impindepall;
off1i = allbstsig.off1i;
locxyprojalln = allbstnoi.locxyprojall;
tarvlsplstalln = allbstnoi.impindepall(:,1);
nsrcn = allbstnoi.nsrc;
impn = allbstnoi.impindepall;
off1in = allbstnoi.off1i;
supertstr = 'Secondary sources removed';
fnsuffix = [];

% %%%param for further checked at KLNB
% locxyprojall = allbstsig.locxyproj4thall;
% tarvlsplstall = allbstsig.impindep4thall(:,1);
% nsrc = allbstsig.nsrc4th;
% imp = allbstsig.impindep4thall;
% locxyprojalln = allbstnoi.locxyproj4thall;
% tarvlsplstalln = allbstnoi.impindep4thall(:,1);
% nsrcn = allbstnoi.nsrc4th;
% impn = allbstnoi.impindep4thall;
% supertstr = 'Further checked at KLNB';
% fnsuffix = '4th';

% keyboard

%% Seismic amplitude or envelope associated with each detection
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

%filtering passband for reading data, confirmed by 'spectrabursts002_4s.m'
hisig=6.3; % this will give a similar spectral shape between template and signal
losig=1.8;

% dates and ets
%%% NOTE: 'dates' is the dates that the tremor is active at the region of interest, so that the
%%% waveform at the stations on these dates you are looking at mainly comes from the sources at the
%%% region of interest. We want to see if there is a noticable change in spectra during the burst
%%% windows on these dates
dates = unique(trange(:,1));
years = unique(floor(dates/1000));
nets = length(years);

%%%moving window length in samples for running CC, envelope, etc.
%standard window length is about 0.5s, this is about the visual duration of the filtered and unfiltered
%template, although in fact to include as least one cycle of the main dipole of template
rccmwsec=0.5;
sps = 160;
rccmwlen=rccmwsec*sps;

wlensec = 6;

%%%Flag to indicate if it is necessary to recalculate everything
flagrecalc = 0;
% flagrecalc = 1;

savefile = strcat('medseisampbetweenlfes',fnsuffix,num2str(wlensec),'s.mat');

if flagrecalc

seisamp=[];
seisampn=[];

for iii = 1: nbst
  
  [iets,i,j] = indofburst(trange,iii);
  
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
            
    %bursts and 4-s detections of the same day
    rangetemp = trange(trange(:,1)==datesets(i), :);
    hfdayi = hfbnd(hfbnd(:,daycol)==datesets(i), :);  % inside bound of the day
    
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
    
    k = iii;
    disp(k);
    
    tmaxi = hfdayi(:, seccol); % starting time of max power rate of half sec inside the ellipse
    tcnti = hfdayi(:, 15);  % the center of detecting win is the 15th col
    
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
%       overshoot = rccmwlen/2;
%       overshoot = 0;
      overshoot = wlensec/2*sps;
    
      %%%%2022/09/26, obtain all information for data before decon, so that you know the threshold
      %%%%being used, although this was used mainly for noise experiment, we still use it here 
      %chop a record segment
      optseg = STAopt(max(floor(tstbuf*sps+1-overshoot-msftaddm),1): ...
        min(floor(tedbuf*sps+overshoot+msftaddm),86400*sps), :); % sta 1
      ortseg = STAort(max(floor(tstbuf*sps+1-overshoot-msftaddm),1): ...
        min(floor(tedbuf*sps+overshoot+msftaddm),86400*sps), :);
      
      %%%read in a single best alignment based on the entire win       
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
      sigstacut = detrend(sigsta(overshoot+1:end-overshoot, :));  %excluding the overshoot

      %for ort. comp
      sigstaort = zeros(size(ortdat,1), nsta);
      for ista = 1:nsta
        tmp = ortdat(:,ista+1); %best aligned, filtered
        tmp = detrend(tmp);
        sigstaort(:,ista) = tmp;
      end
      sigstaort = detrend(sigstaort(overshoot+1:end-overshoot, :));  %excluding the overshoot
      
      %envelope
      env=envelope(detrend(sigsta));
      envcut=envelope(detrend(sigstacut));
      lsig=size(sigsta,1);
      lsigcut=size(sigstacut,1);
      
      %LFEs of the same burst
      ist = sum(nsrc(1:iii-1))+1;
      ied = ist+nsrc(iii)-1;
      impi = imp(ist:ied,:);
      
      ampi=zeros(nsrc(iii),2);
      for j=1:nsrc(iii)
        tcnt=impi(j,1);
        if tcnt<1 || tcnt>lsigcut
          sprintf('detection %d in burst %d is out of range',j,iii);
        end
        
        tst=tcnt-wlensec/2*sps+1;
        ted=tcnt+wlensec/2*sps;
        
        if tst>=1 && ted<=lsigcut          
          ampi(j,1)=median(sigstacut(tst:ted,1));
          ampi(j,2)=median(envcut(tst:ted,1));
        elseif (tst<1 && tst>=1-overshoot) || (ted>lsigcut && ted<=lsigcut+overshoot)
          ampi(j,1)=median(sigsta(overshoot+tst:overshoot+ted,1));
          ampi(j,2)=median(env(overshoot+tst:overshoot+ted,1));
        else
          sprintf('window centered at detection %d in burst %d is out of range',j,iii);
        end
      end
      
      seisamp=[seisamp;ampi];
            
      %%%%%%%%%%%%%%%%%%%%%% now deal with noise %%%%%%%%%%%%%%%%%%%%%%%%
      %obtain the amp and phase spectra of records via fft
      nfft = size(optseg,1); % number of points in fft
      [xf,ft,amp,pha] = fftspectrum(optseg(:,2:end), nfft, sps,'twosided');
      
      %uniform, random phase with the same span [-pi,pi];
      mpharan = minmax(pha');
      seed = iii;
      %         seed = k;
      %         seed = 3;
      rng(seed);
      pharand = (rand(nfft,nsta)-0.5)*2*pi;  %make the phases span from -pi to pi

      %construct record with the same amplitude but random phase
      xfrand = amp.*nfft.*exp(1i.*pharand);
      optseg(:,2:end) = real(ifft(xfrand,nfft));
      
      %%%for orthogonal components, with the same phase??
      [xf,ft,amp,pha] = fftspectrum(ortseg(:,2:end), nfft, sps,'twosided');
      %         pharand = (rand(nfft,3)-0.5)*2*pi;
      xfrand = amp.*nfft.*exp(1i.*pharand);
      ortseg(:,2:end) = real(ifft(xfrand,nfft));
      
      %%%read in a single best alignment based on the entire win       
      %%%Align and compute the RCC based on the entire win, and take that as the input signal!      
      optdat = [];  % win segment of interest
      ortdat = [];
      optdat(:, 1) = optseg(1+msftaddm: end-msftaddm, 1); % time column
      ortdat(:, 1) = ortseg(1+msftaddm: end-msftaddm, 1);
      for ista = 1: nsta 
        optdat(:, ista+1) = optseg(1+msftaddm-off1in(k,ista): end-msftaddm-off1in(k,ista), ista+1);
        ortdat(:, ista+1) = ortseg(1+msftaddm-off1in(k,ista): end-msftaddm-off1in(k,ista), ista+1);
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
      sigstacut = detrend(sigsta(overshoot+1:end-overshoot, :));  %excluding the overshoot

      %for ort. comp
      sigstaort = zeros(size(ortdat,1), nsta);
      for ista = 1:nsta
        tmp = ortdat(:,ista+1); %best aligned, filtered
        tmp = detrend(tmp);
        sigstaort(:,ista) = tmp;
      end
      sigstaort = detrend(sigstaort(overshoot+1:end-overshoot, :));  %excluding the overshoot
      
      %envelope
      env=envelope(detrend(sigsta));
      envcut=envelope(detrend(sigstacut));
      lsig=size(sigsta,1);
      lsigcut=size(sigstacut,1);
      
      %LFEs of the same burst      
      ist = sum(nsrcn(1:iii-1))+1;
      ied = ist+nsrcn(iii)-1;
      impin = impn(ist:ied,:);
      
      ampin=zeros(nsrcn(iii),2);
      for j=1:nsrcn(iii)
        tcnt=impin(j,1);
        if tcnt<1 || tcnt>lsigcut
          sprintf('detection %d in burst %d is out of range',j,iii);
        end
        
        tst=tcnt-wlensec/2*sps+1;
        ted=tcnt+wlensec/2*sps;
        
        if tst>=1 && ted<=lsigcut          
          ampin(j,1)=median(sigstacut(tst:ted,1));
          ampin(j,2)=median(envcut(tst:ted,1));
        elseif (tst<1 && tst>=1-overshoot) || (ted>lsigcut && ted<=lsigcut+overshoot)
          ampin(j,1)=median(sigsta(overshoot+tst:overshoot+ted,1));
          ampin(j,2)=median(env(overshoot+tst:overshoot+ted,1));
        else
          sprintf('window centered at detection %d in burst %d is out of range',j,iii);
        end
      end
      
      seisampn=[seisampn;ampin];
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
  save(strcat(rstpath, '/MAPS/',savefile), 'seisampn','seisamp');
  
else
  load(strcat(rstpath, '/MAPS/',savefile));
end

%% fraction of 'isolated' events, eg, when m=1, evts whose minimum interevt time >0.375s
%%%2 definitions of inter-event times, one is what we have been used
%%%the other is the smaller one of the diff time to the left and right
m = 1;
dtcut = 0.25*m+0.125;

mindtinter=[];
dtinter=[];
minmamp=[];
mamp=[];
for i = 1: nbst
  if nsrc(i) == 0
    continue
  end
  ist = sum(nsrc(1:i-1))+1;
  ied = ist+nsrc(i)-1;
  impi = imp(ist:ied,:);
  if nsrc(i)>=m+1
    %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
    dtfor = diffcustom(impi(:,1), m,'forward'); %to its preceding one
    dtforpad = [zeros(m,1); dtfor];
    dtback = diffcustom(impi(:,1), m,'backward'); %to its following one
    dtbackpad = [dtback; zeros(m,1)];
    tmp1 = [dtforpad dtbackpad];  %time to N-m and N+m for each N
    %choose the min time to neighbors to find isolated ones
    tmp2 = [dtbackpad(1:m); min(tmp1(m+1: end-m, :),[],2); dtforpad(end-m+1:end)];
  else
    dtfor = [];
    tmp2 = [];
  end
  mindtinter = [mindtinter; tmp2];
  dtinter = [dtinter; dtfor];
  
  ampi = seisamp(ist:ied,:);  
  aaa=ampi(1:end-m,:);
  bbb=ampi(1+m:end,:);
  mampi=(aaa+bbb)/2;
  
  tmp3=(ampi(1:m,:) + ampi(1+m:m+m,:))/2;
  tmp4=(ampi(end-m+1:end,:) + ampi(end-m+1-m:end-m,:))/2;
  [~,ind] = min(tmp1(m+1: end-m, :),[],2);
  tmp5=[];
  for j=1:length(ind)
    if ind(j)==1
      tmp5(j,:)=(ampi(j+m,:) + ampi(j,:))/2;
    elseif ind(j)==2
      tmp5(j,:)=(ampi(j+m,:) + ampi(j+m+m,:))/2;
    end
  end
  minmampi = [tmp3; tmp5; tmp4];
  
  mamp=[mamp; mampi];
  minmamp = [minmamp; minmampi];
end

%%%%%%%%%%%% for noise
mindtintern=[];
dtintern=[];
minmampn=[];
mampn=[];
for i = 1: nbst
  if nsrcn(i) == 0
    continue
  end
  ist = sum(nsrcn(1:i-1))+1;
  ied = ist+nsrcn(i)-1;
  impin = impn(ist:ied,:);
  if nsrcn(i)>=m+1
    dtforn = diffcustom(impin(:,1), m,'forward'); %to its preceding one
    dtforpadn = [zeros(m,1); dtforn];
    dtbackn = diffcustom(impin(:,1), m,'backward'); %to its following one
    dtbackpadn = [dtbackn; zeros(m,1)]; 
    tmp1n = [dtforpadn dtbackpadn];  %time to N-m and N+m for each N     
    tmp2n = [dtbackpadn(1:m); min(tmp1n(m+1: end-m, :),[],2); dtforpadn(end-m+1:end)];
  else
    dtforn = [];
    tmp2n = [];
  end   
  mindtintern = [mindtintern; tmp2n];
  dtintern = [dtintern; dtforn];
    
  ampin = seisampn(ist:ied,:);
  if nsrcn(i)>=m+1
    aaa=ampin(1:end-m,:);
    bbb=ampin(1+m:end,:);
    mampin=(aaa+bbb)/2;  %average amplitude of each to its preceding 
    
    tmp3n=(ampin(1:m,:) + ampin(1+m:m+m,:))/2;  %the first few with no preceding
    tmp4n=(ampin(end-m+1:end,:) + ampin(end-m+1-m:end-m,:))/2;  %the last few with no following
    [~,ind] = min(tmp1n(m+1: end-m, :),[],2);
    tmp5n=[]; %for the rest, find the amp of the nearest neighbour
    for j=1:length(ind)
      if ind(j)==1
        tmp5n(j,:)=(ampin(j+m,:) + ampin(j,:))/2;
      elseif ind(j)==2
        tmp5n(j,:)=(ampin(j+m,:) + ampin(j+m+m,:))/2;
      end
    end
    minmampin = [tmp3n; tmp5n; tmp4n];    
  else
    mampin=[];
    minmampin=[];
  end
  
  mampn=[mampn; mampin];  
  minmampn = [minmampn; minmampin];

end
%if the smaller one of the time to N-m and N+m for the source N is big, it is
%isolated
nevtiso = sum(mindtinter/sps>dtcut);
fraciso = nevtiso/length(imp);
nevtison = sum(mindtintern/sps>dtcut);
fracison = nevtison/length(impn);


%% bin the interevent time by the seismic amp of at times of associated sources
nbin=5;
% binedge=0:0.1:12;
binedge=[0 0.375:0.25:ceil(max(dtinter)/sps)];  % 1st bin [0 0.375], then 0.25 increment
% color = jet(nbin);
% color = viridis(nbin);
% color = flipud(gradientblue(nbin));
color = gradientblue(nbin);
% color = ['m';'b';'c';'k';'r';];
% xran = [0 ceil(max(dtinter)/sps)];
xran = [0 7];
yran = [1e-4 1];

widin = 6.5;  % maximum width allowed is 8.5 inches
htin = 5.5;   % maximum height allowed is 11 inches
nrow = 2;
ncol = 2;
f = initfig(widin,htin,nrow,ncol); %initialize fig
pltxran = [0.08 0.98]; pltyran = [0.08 0.98]; % optimal axis location
pltxsep = 0.02; pltysep = 0.09;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); set(ax,'YScale','log');
[ampbin,indbin,n] = binxeqnum(mamp(:,2),nbin);  %bin by amp with same number
plot(ax,[0 30],[1/mode(n) 1/mode(n)],'k--');
medmampbin = zeros(nbin,1);  %median amp of each amp bin
N1dn = [];
for i = 1: nbin
  indi = indbin{i};
  medmampbin(i) = median(ampbin{i});
  dtinterbin{i} = dtinter(indi);
  N1d=histcounts(dtinterbin{i}/sps,binedge,'normalization','count');
  N1d=[N1d N1d(end)];
  N1dn(:,i) = reshape(N1d,[],1)./mode(n);
  p1d(i)=stairs(ax,binedge,N1dn(:,i),'color',color(i,:),'LineWidth',1);
  label1d{i} = sprintf('amp of %.2f',medmampbin(i));
end
legend(ax,p1d(nbin:-1:1),label1d{nbin:-1:1},'FontSize',8);
% xlabel(ax,'Time (s) from each to its preceding');
xlabel(ax,'Time delay (s)');
ylabel(ax,'Normalized count');
xlim(ax,xran);
ylim(ax,yran);
text(ax,0.02,0.40,'Data','Units','normalized','FontSize',11);
text(ax,0.02,0.25,sprintf('From each to its \npreceding'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.08,'a','FontSize',11,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
longticks(ax,2);
hold(ax,'off');

ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); set(ax,'YScale','log');
[ampbin,indbin,n] = binxeqnum(minmamp(:,2),nbin);  %bin by amp with same number
plot(ax,[0 30],[1/mode(n) 1/mode(n)],'k--');
medminmampbin = zeros(nbin,1);  %median amp of each amp bin
N2dn = [];
for i = 1: nbin
  indi = indbin{i};
  medminmampbin(i) = median(ampbin{i});
  mindtinterbin{i} = mindtinter(indi);
  N2d=histcounts(mindtinterbin{i}/sps,binedge,'normalization','count');
  N2d=[N2d N2d(end)];
  N2dn(:,i) = reshape(N2d,[],1)./mode(n);
  p2d(i)=stairs(ax,binedge,N2dn(:,i),'color',color(i,:),'LineWidth',1);
  label2d{i} = sprintf('amp of %.2f',medminmampbin(i));
end
legend(ax,p2d(nbin:-1:1),label2d{nbin:-1:1},'FontSize',8);
xlabel(ax,'Time delay (s)');
% xlabel(ax,'Time (s) from each to its nearest neighbor');
% ylabel(ax,'Normalized count');
xlim(ax,xran);
ylim(ax,yran);
text(ax,0.02,0.40,'Data','Units','normalized','FontSize',11);
text(ax,0.02,0.25,sprintf('From each to its \nnearest neighbor'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.08,'b','FontSize',11,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
longticks(ax,2);
nolabels(ax,2);
hold(ax,'off');

ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); set(ax,'YScale','log');
[ampbinn,indbinn,nn] = binxeqnum(mampn(:,2),nbin);  %bin by amp with same number
plot(ax,[0 30],[1/mode(n) 1/mode(n)],'k--');
medmampbinn = zeros(nbin,1);  %median amp of each amp bin
N1nn = [];
for i = 1: nbin
  indi = indbinn{i};
  medmampbinn(i) = median(ampbinn{i});
  dtinterbinn{i} = dtintern(indi);
  N1n=histcounts(dtinterbinn{i}/sps,binedge,'normalization','count');
  N1n=[N1n N1n(end)];
  N1nn(:,i) = reshape(N1n,[],1)./mode(n);
  p1n(i)=stairs(ax,binedge,N1nn(:,i),'color',color(i,:),'LineWidth',1);
  label1n{i} = sprintf('amp of %.2f',medmampbinn(i));
end
legend(ax,p1n(nbin:-1:1),label1n{nbin:-1:1},'FontSize',8);
% xlabel(ax,'Time (s) from each to its preceding');
ylabel(ax,'Normalized count');
xlim(ax,xran);
ylim(ax,yran);
text(ax,0.02,0.40,'Noise','Units','normalized','FontSize',11);
text(ax,0.02,0.25,sprintf('From each to its \npreceding'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.08,'c','FontSize',11,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
longticks(ax,2);
nolabels(ax,1);
hold(ax,'off');

ax=f.ax(4); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); set(ax,'YScale','log');
[ampbinn,indbinn,nn] = binxeqnum(minmampn(:,2),nbin);  %bin by amp with same number
plot(ax,[0 30],[1/mode(n) 1/mode(n)],'k--');
medminmampbinn = zeros(nbin,1);  %median amp of each amp bin
N2nn = [];
for i = 1: nbin
% for i = nbin:-1:1 
  indi = indbinn{i};
  medminmampbinn(i) = median(ampbinn{i});
  mindtinterbinn{i} = mindtintern(indi);
  N2n=histcounts(mindtinterbinn{i}/sps,binedge,'normalization','count');
  N2n=[N2n N2n(end)];
  N2nn(:,i) = reshape(N2n,[],1)./mode(n);
  p2n(i)=stairs(ax,binedge,N2nn(:,i),'color',color(i,:),'LineWidth',1);
  label2n{i} = sprintf('amp of %.2f',medminmampbinn(i));
end
legend(ax,p2n(nbin:-1:1),label2n{nbin:-1:1},'FontSize',8);
% xlabel(ax,'Time (s) from each to its nearest neighbor');
% ylabel(ax,'Normalized count');
xlim(ax,xran);
ylim(ax,yran);
text(ax,0.02,0.40,'Noise','Units','normalized','FontSize',11);
text(ax,0.02,0.25,sprintf('From each to its \nnearest neighbor'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.08,'d','FontSize',11,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
longticks(ax,2);
nolabels(ax,3);
hold(ax,'off');

fname = strcat('lfeintertime.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
keyboard


%% distribution without binning by amp 
widin = 6.5;  % maximum width allowed is 8.5 inches
htin = 3.5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f = initfig(widin,htin,nrow,ncol); %initialize fig
pltxran = [0.08 0.98]; pltyran = [0.12 0.98]; % optimal axis location
pltxsep = 0.02; pltysep = 0.09;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

xran = [0 7];
yran = [1e-4 1];

ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); set(ax,'YScale','log');
[N1d]=histcounts(dtinter/sps,binedge,'normalization','count');
[N1n]=histcounts(dtintern/sps,binedge,'normalization','count');
N1d=[N1d N1d(end)];
N1dn = reshape(N1d,[],1)./length(dtinter);
N1n=[N1n N1n(end)];
N1nn = reshape(N1n,[],1)./length(dtinter);
p1d=stairs(ax,binedge,N1dn,'b','LineWidth',1);
p1n=stairs(ax,binedge,N1nn,'r','LineWidth',1);
fitxran = [2 xran(2)]; 
ind = find(binedge>=fitxran(1) & binedge<=fitxran(2));
fitstruct=robustexpfit(binedge(ind),N1dn(ind),'log10',[1e3 -1]);
fitobj=fitstruct.fitobj;
coef=coeffvalues(fitobj); slp1d=coef(2);
yfit1d = feval(fitobj,binedge);
plot(ax,binedge,yfit1d,'b--');
fitstruct=robustexpfit(binedge(ind),N1nn(ind),'log10',[1e3 -1]);
fitobj=fitstruct.fitobj;
coef=coeffvalues(fitobj); slp1n=coef(2);
yfit1n = feval(fitobj,binedge);
plot(ax,binedge,yfit1n,'r--');
% xlabel(ax,'Time (s) from each to its preceding');
xlabel(ax,'Time delay (s)');
ylabel(ax,'Normalized count');
xlim(ax,xran);
ylim(ax,yran);
legend(ax,[p1d,p1n],'Data','Synthetic noise');
text(ax,0.02,0.35,sprintf('[%.1f, %.1f] s',fitxran(1),fitxran(2)),'Units',...
  'normalized','FontSize',9);
text(ax,0.02,0.25,sprintf('From each to its \npreceding'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.08,'a','FontSize',11,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
% [lambdahat,lambdaci] = poissfit(dtinter/sps)
% x=0:1:12;
% y=poisspdf(x,lambdahat);
% plot(ax,x,y,'r-','linew',2);
% keyboard
longticks(ax,2);
hold(ax,'off');

ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); set(ax,'YScale','log');
[N2d]=histcounts(mindtinter/sps,binedge,'normalization','count');
[N2n]=histcounts(mindtintern/sps,binedge,'normalization','count');
N2d=[N2d N2d(end)];
N2dn = reshape(N2d,[],1)./length(mindtinter);
N2n=[N2n N2n(end)];
N2nn = reshape(N2n,[],1)./length(mindtinter);
stairs(ax,binedge,N2dn,'b','LineWidth',1);
stairs(ax,binedge,N2nn,'r','LineWidth',1);
fitxran = [1.5 4]; 
ind = find(binedge>=fitxran(1) & binedge<=fitxran(2));
fitstruct=robustexpfit(binedge(ind),N2dn(ind),'log10',[1e3 -1]);
fitobj=fitstruct.fitobj;
coef=coeffvalues(fitobj); slp2d=coef(2);
yfit2d = feval(fitobj,binedge);
plot(ax,binedge,yfit2d,'b--');
fitstruct=robustexpfit(binedge(ind),N2nn(ind),'log10',[1e3 -1]);
fitobj=fitstruct.fitobj;
coef=coeffvalues(fitobj); slp2n=coef(2);
yfit2n = feval(fitobj,binedge);
plot(ax,binedge,yfit2n,'r--');
% xlabel(ax,'Time (s) from each to its preceding');
xlabel(ax,'Time delay (s)');
% ylabel(ax,'Normalized count');
xlim(ax,xran);
ylim(ax,yran);
text(ax,0.02,0.35,sprintf('[%.1f, %.1f] s',fitxran(1),fitxran(2)),'Units',...
  'normalized','FontSize',9);
text(ax,0.02,0.25,sprintf('From each to its \nnearest neighbor'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.08,'b','FontSize',11,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
longticks(ax,2);
nolabels(ax,2);
hold(ax,'off');

fname = strcat('lfeintertimelump.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
% keyboard



%%
f=initfig(12,9,2,2); %initialize fig
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); set(ax,'YScale','log');
% binedge=0:0.1:12;
[N1d]=histcounts(dtinter/sps,binedge,'normalization','count');
[N1n]=histcounts(dtintern/sps,binedge,'normalization','count');
N1d=[N1d N1d(end)];
N1n=[N1n N1n(end)];
p1d=stairs(ax,binedge,N1d,'b','LineWidth',1);
p1n=stairs(ax,binedge,N1n,'r','LineWidth',1);
% p1=histogram(ax,dtinter/sps,'binedges',binedge,'normalization','count','Facec','b');
% p2=histogram(ax,dtintern/sps,'binedges',binedge,'normalization','count','Facec','r');
fitstruct=robustexpfit(binedge(6:end),N1d(6:end),'log10',[1e3 -1]);
fitobj=fitstruct.fitobj;
coef=coeffvalues(fitobj); slp1d=coef(2);
yfit1d = feval(fitobj,binedge);
plot(ax,binedge,yfit1d,'b--');
fitstruct=robustexpfit(binedge(6:end),N1n(6:end),'log10',[1e3 -1]);
fitobj=fitstruct.fitobj;
coef=coeffvalues(fitobj); slp1n=coef(2);
yfit1n = feval(fitobj,binedge);
plot(ax,binedge,yfit1n,'r--');
% xlabel(ax,'Time (s) from each to its preceding');
% ylabel(ax,'Count');
xlim(ax,[0 6]);
legend(ax,[p1d,p1n],'Data','Synthetic noise');
% [lambdahat,lambdaci] = poissfit(dtinter/sps)
% x=0:1:12;
% y=poisspdf(x,lambdahat);
% plot(ax,x,y,'r-','linew',2);
% keyboard
ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); set(ax,'YScale','log');
N1=N1d-N1n;
p1=stairs(ax,binedge,N1,'k','LineWidth',1);
% fitstruct=robustexpfit(binedge(6:50),N1(6:50),'log10',[1e3 -1]);
% fitobj=fitstruct.fitobj;
% coef=coeffvalues(fitobj); slp1=coef(2);
% yfit1 = feval(fitobj,binedge);
% plot(ax,binedge,yfit1,'k--');
xlabel(ax,'Time (s) from each to its preceding');
ylabel(ax,'Count');
xlim(ax,[0 6]);
ylim(ax,[0 1e4]);
legend(ax,p1,'Data-Synthetic noise');
ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); set(ax,'YScale','log');
[N2d]=histcounts(mindtinter/sps,binedge,'normalization','count');
[N2n]=histcounts(mindtintern/sps,binedge,'normalization','count');
N2d=[N2d N2d(end)];
N2n=[N2n N2n(end)];
stairs(ax,binedge,N2d,'b','LineWidth',1);
stairs(ax,binedge,N2n,'r','LineWidth',1);
% histogram(ax,mindtinter/sps,'binedges',binedge,'normalization','count','Facec','b');
% histogram(ax,mindtintern/sps,'binedges',binedge,'normalization','count','Facec','r');
% fitstruct=robustexpfit(binedge(6:end),N2d(6:end),'log10',[1e3 -1]);
% fitobj=fitstruct.fitobj;
% coef=coeffvalues(fitobj); slp2d=coef(2);
% yfit2d = feval(fitobj,binedge);
% plot(ax,binedge,yfit2d,'b--');
% fitstruct=robustexpfit(binedge(6:end),N2n(6:end),'log10',[1e3 -1]);
% fitobj=fitstruct.fitobj;
% coef=coeffvalues(fitobj); slp2n=coef(2);
% yfit2n = feval(fitobj,binedge);
% plot(ax,binedge,yfit2n,'r--');
xlabel(ax,'Time (s) from each to nearest neighbor');
ylabel(ax,'Count');
xlim(ax,[0 6]);
% [lambdahat,lambdaci] = poissfit(mindtinter/sps)
% x=0:1:12;
% y=poisspdf(x,lambdahat);
% plot(ax,x,y,'r-','linew',2);
ax=f.ax(4); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); set(ax,'YScale','log');
N2=N2d-N2n;
stairs(ax,binedge,N2,'b','LineWidth',1);
% fitstruct=robustexpfit(binedge(6:50),N2(6:50),'log10',[1e3 -1]);
% fitobj=fitstruct.fitobj;
% coef=coeffvalues(fitobj); slp2=coef(2);
% yfit2 = feval(fitobj,binedge);
% plot(ax,binedge,yfit2,'k--');
xlabel(ax,'Time (s) from each to nearest neighbor');
ylabel(ax,'Count');
xlim(ax,[0 6]);
ylim(ax,[0 1e4]);

keyboard
