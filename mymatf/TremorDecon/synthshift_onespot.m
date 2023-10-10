% function synthshift_onespot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to 'synthshift_chao.m', this code is to make synthetic seismograms
% from LFE templates.
% But this code tries to put all srcs at one and only one spot, then add
% different level of noise, which is diff precent of noise that is made of
% the same amp spectrum of data and randomized phase, times the envelope of
% data (so that the noise has the same fluctuation as data)
% 1. in order to see what percentage of noise would make it look like data
% 2. now there can be 2 end-member synthetic tests: one is different saturation
%  level, diff source area, and no noise (or little); one is same spot source,
%  plus diff satur level and diff noise level
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/09/13
% Last modified date:   2023/09/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
format short e   % Set the format to 5-digit floating point
clear
clc
close all

%%%Flag to indicate if it is necessary to recalculate everything
% flagrecalc = 0;
flagrecalc = 1;

if ~flagrecalc
  load('rst_synth_onespot.mat');
else
  %% for easy testing
  defval('normflag',0); %whether to normalize templates
  defval('rccmwsec',0.5); %moving win len in sec for computing RCC
  
  set(0,'DefaultFigureVisible','on');
  %set(0,'DefaultFigureVisible','off');
  
  workpath = getenv('ALLAN');
  datapath = strcat(workpath,'/data-no-resp');
  temppath = strcat(datapath, '/templates/PGCtrio/');
  rstpath = strcat(datapath, '/PGCtrio');
  
  [scrsz, resol] = pixelperinch(1);
  
  % tempflag = 'allan';
  tempflag = 'chao';
  
  adatapath = '/home/data2/chaosong/matlab/allan/matfils/';  %path for Allan's data
  
  stas=['PGC  '
    'SSIB '
    'SILB '
    %   'LZB  '
    %   'TWKB '
    %   'MGCB '
    'KLNB '
    ]; % determine the trio and order, here the 1st sta is PGC
  nsta=size(stas,1);         %  number of stations
  
  if strcmp(tempflag,'chao')
    fam = '002';   % family number
    sps = 160;
    templensec = 60;
    ccstack = [];
    for ista = 1: nsta
      fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
        num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao_catnew');
      ccstack(:,ista) = load(fname);
    end
    STA = detrend(ccstack);
  elseif strcmp(tempflag,'allan')
    sps = 100;
    templensec = 120;
    fname = ['PGCopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh_0.5-15Hz_CC1.25-6.5Hz_le20shift';
      'SSIopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh_0.5-15Hz_CC1.25-6.5Hz_le20shift';
      'SILopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh_0.5-15Hz_CC1.25-6.5Hz_le20shift'];
    %   templensec = 60;
    %   fname = ['PGCopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh0.5-15Hz_le14sh';
    %            'SSIopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh0.5-15Hz_le14sh';
    %            'SILopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh0.5-15Hz_le14sh'];
    for ista=1:nsta
      temp=load(strcat(adatapath,fname(ista,:)));
      STA(:,ista)=detrend(temp(:,1))/350;
    end
  end
  
  % %plot the raw templates, not filtered, not best aligned
  % figure
  % ltemp = size(STA,1);
  % subplot(111)
  % hold on
  % for ista = 1: nsta
  %   plot((1:ltemp)/sps,STA(:,ista));
  % end
  
  %%%The below aligns the templates by x-correlation
  if strcmp(tempflag,'chao')
    ist = templensec*sps*4/10;  %not using whole window in case any station has very long-period energy
    ied = templensec*sps*6/10;
  elseif strcmp(tempflag,'allan')
    ist = templensec*sps*6/10;
    ied = templensec*sps*8/10;
  end
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
  end
  %normalization
  if normflag
    for ista=1:nsta
      STAtmp(:,ista)=STAtmp(:,ista)/spread(ista); % now templates are 'aligned' indeoendently by x-corr wrt. sta 1
    end
  end
  % figure
  % ltemp = size(STAtmp,1);
  % subplot(111)
  % hold on
  % for ista = 1: nsta
  %   plot((1:ltemp)/sps,STAtmp(:,ista));
  % end
  %%%The above aligns the templates by x-correlation
  
  %%%detrend, taper and bandpass templates
  tmpwlet = STAtmp; % no bandpass
  tmpwletf = STAtmp;  % bandpassed version
  fractap = sps/size(tmpwlet,1);
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
  if nsta>3
    for ista = 4: nsta
      [mcoef,offwlet1i(ista)] = xcorrmax(tmpwletf(:,1), tmpwletf(:,ista), mshiftadd, 'coeff');
    end
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
  for ista = 1: nsta
    %cut according to the zero-crossing and the time shift from the constrained CC
    green(:,ista) = tmpwlet(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
    greenf(:,ista) = tmpwletf(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
    %detrend again for caution
    green(:,ista) = detrend(green(:,ista));
    greenf(:,ista) = detrend(greenf(:,ista));
    
    %re-find the zero-crossing as the template length has changed
    [~,imin] = min(green(:,ista));
    [~,imax] = max(green(:,ista));
    [~,zcrosses(ista)] = min(abs(green(imin:imax,ista)));
    zcrosses(ista) = zcrosses(ista)+imin-1;
    ppeaks(ista) = imax;
    npeaks(ista) = imin;
    
    [~,imin] = min(greenf(:,ista));
    [~,imax] = max(greenf(:,ista));
    [~,zcrossesf(ista)] = min(abs(greenf(imin:imax,ista)));
    zcrossesf(ista) = zcrossesf(ista)+imin-1;
    ppeaksf(ista) = imax;
    npeaksf(ista) = imin;
    
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
  if nsta>3
    for ista = 4: nsta
      [mcoef,mlag] = xcorrmax(greenf(:,1), greenf(:,ista), mshiftadd, 'coeff');
      if mlag~=0   % offset in samples
        fprintf('Filtered templates are NOT best aligned at %s \n',stas(ista,:));
      end
    end
  end
  
  amprat(1,:) = minmax(greenf(:,1)')./minmax(greenf(:,2)');	% amp ratio between max at sta 3 and 2 or min
  amprat(2,:) = minmax(greenf(:,1)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min
  amprat(3,:) = minmax(greenf(:,2)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min
  spread = range(greenf);   % range of the amp of template
  
  %%%plot the unfiltered and filtered templates
  % plt_templates(green,greenf,stas,[],[],lowlet,hiwlet,sps);
  
  %just the filtered templates
  % plt_templates_bp(greenf,stas,lowlet,hiwlet,sps);
  % zcrosses
  % ppeaks
  % npeaks
  % zcrossesf
  % ppeaksf
  % npeaksf
  % keyboard
  
  
  %% make synthetic noise that has the same fluctuation as data
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
  cutout = 'ellipse';
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
  
  sps = 160;
  
  %filtering passband for reading data, confirmed by 'spectrabursts002_4s.m'
  hisig=6.3; % this will give a similar spectral shape between template and signal
  losig=1.8;
  
  %%%moving window length in samples for running CC, envelope, etc.
  %standard window length is about 0.5s, this is about the visual duration of the filtered and unfiltered
  %template, although in fact to include as least one cycle of the main dipole of template
  rccmwlen=rccmwsec*sps;
  % rccmwlen=sps/2;
  % rccmwlen=sps;
  
  idxbst = 181;
  
  [iets,i,j] = indofburst(trange,idxbst);
  
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
    % continue    % continue to the next day
  end
  
  k = idxbst;
  disp(k);
  
  tmaxi = hfdayi(:, seccol); % starting time of max power rate of half sec inside the ellipse
  tmaxo = hfdayo(:, seccol); % starting time of max power rate of half sec outside the ellipse
  tcnti = hfdayi(:, 15);  % the center of detecting win is the 15th col
  
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
  
  %%%according the linear fitting result of off12 and off13 VS. origin time, we sort of know how
  %%%much it needs to change the overall offset to change by 1 sample
  %generate overlapping windows of the same length
  indst = floor(tstbuf*sps+1);
  inded = floor(tedbuf*sps);
  
  %% making synthetic noise as the base
  %obtain the amp and phase spectra of records via fft
  nfft = size(optseg,1); % number of points in fft
  [xf,ft,amp,pha] = fftspectrum(optseg(:,2:end), nfft, sps,'twosided');
  
  %uniform, random phase with the same span [-pi,pi];
  mpharan = minmax(pha');
  seedpool = 10:10:40;
  seed = idxbst;
  %       for ise = 1: length(seedpool)
  %       seed = seedpool(ise);
  rng(seed);
  pharand = (rand(nfft,nsta)-0.5)*2*pi;  %make the phases span from -pi to pi
  
  %construct record with the same amplitude but random phase
  xfrand = amp.*nfft.*exp(1i.*pharand);
  optsegnoi = real(ifft(xfrand,nfft));
  
  xnoi1 = optsegnoi; %synthetic noise
  
  xnoi = xnoi1;
  % ind = find(impindepst(:,1)/sps >= xnoi(1) & impindepst(:,1)/sps <= xnoi(2));
  
  %%%Different scaling schemes
  %       %1. To make noise have the median amp of data
  %       sclfact = median(envelope(detrend(optseg(:,2:end))));
  %2. To make noise have the same fluctuation as data
  env = envelope(detrend(optseg(:,2:end)));
  % sclfact = env./range(env);  %normalize
  sclfact = env;
  %       %3. To make noise have the median amp of decon sources from data
  %       sclfact = 6.4461e-01*mean(spread(1:3))/2;
  %scale the noise so that it has amp fluctuation as data
  xnoi = xnoi .* sclfact;
  
  %normalize so that noi has the same median env as data
  envn = envelope(detrend(xnoi));
  xnoi = xnoi .* median(env) ./ median(envn);
  median(env) - median(envelope(detrend(xnoi)))
  
  figure
  subplot(311)
  plot((1:size(optseg,1))/sps,optseg(:,2),'r','linew',1); hold on
  plot((1:size(optseg,1))/sps,optseg(:,3),'b','linew',1);
  plot((1:size(optseg,1))/sps,optseg(:,4),'k','linew',1);
  text(0.95,0.9,'Data','Units','normalized','HorizontalAlignment','right');
  %       xlabel('Samples at 160 sps');
  xlabel('Time (s)');
  ylabel('Amplitude');
  xlim([0 size(optseg,1)]/sps);
  ylim([-1.2 1.2]);
  
  subplot(312)
  plot((1:size(optseg,1))/sps,xnoi(:,1),'r','linew',1); hold on
  plot((1:size(optseg,1))/sps,xnoi(:,2),'b','linew',1);
  plot((1:size(optseg,1))/sps,xnoi(:,3),'k','linew',1);
  text(0.95,0.9,'100% noise','Units','normalized','HorizontalAlignment','right');
  %       xlabel('Samples at 160 sps');
  xlabel('Time (s)');
  ylabel('Amplitude');
  xlim([0 size(optseg,1)]/sps);
  ylim([-1.2 1.2]);
  
  subplot(313)
  plot((1:size(optseg,1))/sps,optseg(:,2)+xnoi(:,1),'r','linew',1); hold on
  plot((1:size(optseg,1))/sps,optseg(:,3)+xnoi(:,2),'b','linew',1);
  plot((1:size(optseg,1))/sps,optseg(:,4)+xnoi(:,3),'k','linew',1);
  text(0.95,0.9,'Data+100% noise','Units','normalized','HorizontalAlignment','right');
  %       xlabel('Samples at 160 sps');
  xlabel('Time (s)');
  ylabel('Amplitude');
  xlim([0 size(optseg,1)]/sps);
  ylim([-1.2 1.2]);
  
  keyboard
  
  %% generate synthetic sources
  %%%Specify the amplitude-frequency (counts) distribution
  distr='UN'  % uniform distribution
  
  % Twin=0.5*3600+3+2*ceil(greenlen/sps); %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
  % winlen=Twin*sps+1;
  winlen = size(xnoi,1)+greenlen-1;
  Twin = winlen/sps;
  skiplen=greenlen;
  %%%for the duration of templates, there are several options
  %%%1. (ppeak-npeak)*2 of the bb template: 44;54;38
  %%%2. direct eyeballing for between zerocrossings: ~65
  %%%3. binned peak-to-peak separation for decently saturated unfiltered synthetics: ~37
  %%%4. similar to 3, but synthetics are filtered first: ~37
  tdura = 0.4;  % duration from Chao's broadband template, width is about 795-730=65 spls at 160 hz
  satn=1/tdura*Twin   % if just saturated, how many templates can be fit in? a single peak is ~20 samples wide; maybe a little less (at 100 sps).
  %Twin is window duration in seconds. Events can fall within Twin of
  %the start, but the synthetics will go to Twin*(sample rate)+Greenlen to
  %avoid checking for subscript overrrun.  When "synths" are written from "synth" in the
  %subroutine, Greenlen from the start and end will not be written.
  fracelsew=0; %0.25 %0.5; %0.6; %The ratio of elsewhere events to local events.  Zero for Discrete Ide.
  nsat=[0.4 1 2 4 10 20 40 100];  % times of saturation
  % nsat=[95 100 200];  % times of saturation
  nnsat = length(nsat);
  writes=round(nsat*satn) %how many templates to throw in, under different degrees of saturation
  
  synth=zeros(winlen+greenlen,nsta);
  
  nouts=length(writes);
  
  %%%specify which time is uniform in time
  % timetype = 'tarvl';
  timetype = 'tori';
  
  %%%specify if considering the physical size of each source
  % physicalsize = 1;
  physicalsize = 0;
  
  %%%specify regime for transformation from time offset to map location
  % ftrans = 'interpArmb';
  % ftrans = 'interpArmbreloc';
  ftrans = 'interpchao';
  
  %%%specify if forcing a min speration of arrival time for 2 events from the same spot
  % forcesep = 1;
  forcesep = 0;
  
  b=999. %>150 for uniform size distribution
  
  %which location transformation to use
  if strcmp(ftrans,'interpArmb')
    xygrid=load('xygridArmb'); %Made from /ARMMAP/MAPS/2021/interpgrid.m; format PGSS, PGSI, dx, dy
    size(xygrid) %this grid is 40 sps!
  elseif strcmp(ftrans,'interpArmbreloc') || strcmp(ftrans,'interpchao')
    [loc, indinput] = off2space002([],sps,ftrans,0);
    % loc has 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
    iup = sps/40;
    %     xygrid(:,3)=griddata(loc(:,7)/iup,loc(:,8)/iup,loc(:,1),xygrid(:,1),xygrid(:,2));
    %     xygrid(:,4)=griddata(loc(:,7)/iup,loc(:,8)/iup,loc(:,2),xygrid(:,1),xygrid(:,2));
    xygrid = [loc(:,7)/iup,loc(:,8)/iup,loc(:,1),loc(:,2)];
  end
  
  %choose only one spot, use [0,0]
  xygrid = xygrid(xygrid(:,1)==0 & xygrid(:,2)==0, :);
  
  diam=0;
  
  [synths,mommax,sources,greensts]=csplaw3c(writes,winlen,skiplen,synth,green,b,...
    xygrid,sps,fracelsew,seed,timetype,ftrans,stas);
  
  %%% a simple plot the synthetics
  figure
  nrow = length(writes)+1;
  ncol = 1;
  subplot(nrow,ncol,1)
  hold on
  tmp = synth;
  % if nsta == 3
  plot(tmp(:,1),'r');
  plot(tmp(:,2),'b');
  plot(tmp(:,3),'k');
  % else
  %   color = jet(nsta);
  %   for ista = 1: nsta
  %     plot(tmp(:,ista),'Color',color(ista,:));
  %   end
  % end
  xlim([0 40*sps]);
  axranexp(gca,6,20);
  
  for i = 1: length(writes)
    subplot(nrow,ncol,i+1)
    hold on
    tmp = synths(:,:,i);
    % if nsta == 3
    plot(tmp(:,1),'r');
    plot(tmp(:,2),'b');
    plot(tmp(:,3),'k');
    % else
    %   color = jet(nsta);
    %   for ista = 1: nsta
    %     plot(tmp(:,ista),'Color',color(ista,:));
    %   end
    % end
    text(0.95,0.9,sprintf('%.1f x saturation',nsat(i)),'Units','normalized','HorizontalAlignment',...
      'right');
    xlim([0 40*sps]);
    axranexp(gca,6,20);
  end
  xlabel(sprintf('Samples at %d Hz',sps),'FontSize',12);
  
  tmpgrid = xygrid;
  tmpgrid(:,1:2)=round(sps/40*tmpgrid(:,1:2)); % *4 to get to 160 sps from 40.
  insat = 2;  %which saturation to look at
  n=writes(insat);
  a = squeeze(sources(1:n,:,insat));
  b = a(any(a,2),:);
  source=b;
  off = tmpgrid(source(:,2),1:2); %note that 'tmpgrid' has the desired sps
  
  close all
  
  %% filter data
  for insat = 1: nnsat
    for ista = 1:nsta
      synths(:,ista,insat) = Bandpass(synths(:,ista,insat), sps, losig, hisig, 2, 2, 'butter');
    end
  end
  
  %% a glance of signal
  %add simulated noise to the current burst window
  insat = 2;
  tmp = synths(:,:,insat) + xnoi;
  
  figure
  subplot(311)
  plot((1:size(xnoi,1))/sps,synths(:,1,insat),'r','linew',1); hold on
  plot((1:size(xnoi,1))/sps,synths(:,2,insat),'b','linew',1);
  plot((1:size(xnoi,1))/sps,synths(:,3,insat),'k','linew',1);
  text(0.95,0.9,sprintf('Synthetics with sat=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','right');
  %       xlabel('Samples at 160 sps');
  xlabel('Time (s)');
  ylabel('Amplitude');
  % xlim([0 size(xnoi,1)]/sps);
  xlim([0 30]);
  ylim([-1.2 1.2]);
  
  subplot(312)
  plot((1:size(xnoi,1))/sps,xnoi(:,1),'r','linew',1); hold on
  plot((1:size(xnoi,1))/sps,xnoi(:,2),'b','linew',1);
  plot((1:size(xnoi,1))/sps,xnoi(:,3),'k','linew',1);
  text(0.95,0.9,'100% noise','Units','normalized','HorizontalAlignment','right');
  %       xlabel('Samples at 160 sps');
  xlabel('Time (s)');
  ylabel('Amplitude');
  % xlim([0 size(xnoi,1)]/sps);
  xlim([0 30]);
  ylim([-1.2 1.2]);
  
  subplot(313)
  plot((1:size(xnoi,1))/sps,tmp(:,1),'r','linew',1); hold on
  plot((1:size(xnoi,1))/sps,tmp(:,2),'b','linew',1);
  plot((1:size(xnoi,1))/sps,tmp(:,3),'k','linew',1);
  text(0.95,0.9,'Synthetics+100% noise','Units','normalized','HorizontalAlignment','right');
  %       xlabel('Samples at 160 sps');
  xlabel('Time (s)');
  ylabel('Amplitude');
  % xlim([0 size(xnoi,1)]/sps);
  xlim([0 30]);
  ylim([-1.2 1.2]);
  
  keyboard
  
  %% compose new synthetic waveform and carry out deconvolution
  %%%flag for validating if ground truth of sources can recover the record
  % testsrcflag = 1;
  testsrcflag = 0;
  
  %%%flag for validing if the spectral shapes of data and templates are similar
  % testfreqflag = 1;
  testfreqflag = 0;
  
  %%%flag for plot the data
  % pltdataflag = 1;
  pltdataflag = 0;
  
  %%%flag for plot the ground truth distribution
  % pltgtflag = 1;
  pltgtflag = 0;
  
  %%%flag for plot the decon src distribution after grouping
  % pltsrcflag1 = 1;
  pltsrcflag1 = 0;
  
  %%%flag for plot the decon src distribution after removing 2ndary src
  % pltsrcflag2 = 1;
  pltsrcflag2 = 0;
  
  %%%flag for plot the decon src distribution after checking at 4th stas
  % pltsrcflag3 = 1;
  pltsrcflag3 = 0;
  
  %different percent of noise
  perctrial = 0.1*(0:1:10)';
  ntrial = length(perctrial);
  
  impgrp = cell(nnsat,ntrial);
  imp = cell(nnsat,ntrial);
  imp4th = cell(nnsat,ntrial);
  
  %%%loop for noise level
  for iperc = 1: ntrial
    
    perc = perctrial(iperc);
    disp(perc);
    
    noiseg = xnoi .*perc;
    
    synnew = synths;
    %%%loop for saturation level
    for insat = 1: nnsat
      %   insat = 1;
      disp(nsat(insat));
      
      %add simulated noise to the current burst window
      tmp = synths(:,:,insat);
      synnew(:,:,insat) = tmp + noiseg;
      
      %%%load synthetics of certain saturation level
      STAopt = synnew(:,:,insat);
      
      %%%load sources
      n=writes(insat);
      a = squeeze(sources(1:n,:,insat));
      b = a(any(a,2),:);
      synsrc=b;
      
      tmp = xygrid(synsrc(:,2),:);
      synsrc = [synsrc(:,1) tmp(:,1:4) ones(length(tmp),1)];  %[indtarvl, off12, off13, loce, locn, amp]
      
      % keyboard
      %%%load starting indices of added sources at sta 1
      synsrcstind = squeeze(greensts{insat}{1});
      
      %% testing, extract and validate the added impulses of template
      if testsrcflag
        wlensec = 30; %how long the window to test
        bufsec = 1; %need some buffer window for later CC alignment
        buffer = bufsec*sps;
        wlensecb = wlensec+bufsec;
        sigstagt = [];  %target window of synthetic signals at all stations
        sigpnstagt = [];  %target window of synthetic signals plus the noise
        sigconvstagt = [];  %target window of reproduced synthetic signals using convolution
        noistagt = [];  %target window of noise
        synsrcgtsta = cell(nsta,1);  %sources whose zero-crossing with time range, [indtarvl rnoff12 rnoff13 amp]
        for ista = 1: nsta
          tgreen = zcrosses(ista);  % anchor time of templates, here choose as zero-crossing time
          indst = 1;  % starting index of the simulated signal to test
          inded = wlensecb*sps+indst-1; % ending index of the simulated signal to test
          source = synsrc;
          if ista == 1
            greenst = synsrcstind; % the starting index of each added template, context of full length
          elseif ista <=3
            %ind - rnoff is the arrival time in index at sta 2 and 3
            %note the sign here, if off12 >0, move 2 to the right to align with 1, meaning that 2 is
            %earlier than 1, ie., tarvl2 < tarvl1. Be consistent all the time
            greenst = synsrcstind-source(:, ista);
          else
            greenst = pred_tarvl_at4thsta(stas(ista,:),source(:,2),source(:,3),source(:,1));
          end
          %%%you don't need all impulses, only some of them contribute to the length of truncated record
          %%%cut out the green indice and sources that contribute
          %%%Note that some zero-crossings might appear later than 'indedsig', but the corresponding start
          %%%index of the green is actually smaller than 'indedsig', ie. their 'greenst' <= 'indedsig'
          induse = find(greenst>=skiplen-greenlen+indst & greenst<=inded+skiplen-1);
          greenst = greenst(induse);
          source = source(induse,:);
          %     keyboard
          
          greenzc = greenst+tgreen; % index of approximate zero-crossing
          source(:, 1) = source(:, 1)+zcrosses(1);  % now 'greenzc' should be the same as 1st col of 'source'
          impamp = zeros(max(greenzc)+20,1);  %20 is a bit arbitary
          for i = 1: length(greenzc)
            impamp(greenzc(i)) = impamp(greenzc(i))+source(i, 6);
          end
          
          %note that some length of the full simulation 'skiplen' was skipped in 'synthgen.m'
          sigpntemp = STAopt(:,ista);  % simulated signal with white noise
          noitemp = synth(skiplen:winlen,ista);   % account for the skipped part
          sigtemp = sigpntemp-noitemp; % subtract the white noise
          sigpnstagt(:,ista) = sigpntemp(indst:inded); % now focus on the part of interest
          sigstagt(:,ista) = sigtemp(indst:inded);
          noistagt(:,ista) = noitemp(indst:inded);
          
          %direct convolution
          sigconvtmp = conv(green(:,ista),impamp,'full');
          indtrc = tgreen+skiplen;  % starting index for truncation
          sigconvstagt(:,ista) = sigconvtmp(indtrc+indst-1:inded+indtrc-1);  % cut accordingly
          
          %%%Can we reproduce the synthetics with truncated convolution? --YES
          figure
          subplot(411)
          plot(green(:,ista),'r');
          xlim([0 greenlen]);
          legend('Template (unfiltered)');
          text(0.05,0.9,stas(ista,:),'HorizontalAlignment','left','Units','normalized');
          title('Reproduce synthetics with truncated convolution');
          subplot(412)
          imptemp = find(impamp>0);
          p1=stem(imptemp-skiplen, impamp(imptemp),'b','MarkerSize',4); hold on;
          % p1=stem((1:size(impamp))-skiplen, impamp,'b','MarkerSize',4); hold on;
          ax = gca;
          plot([indst indst],ax.YLim,'--','color',[.5 .5 .5]);
          plot([inded inded],ax.YLim,'--','color',[.5 .5 .5]);
          legend(p1,'Synthetic random impulses');
          xran1 = ax.XLim; axpos1 = ax.Position(1);
          subplot(413);
          plot(sigstagt(:,ista),'b'); hold on
          ax=gca; yran = ax.YLim; xran2 = ax.XLim;
          shrink(ax,xran1/xran2,1);
          ax.Position(1)=axpos1;
          plot(sigconvstagt(:,ista),'k');
          legend('Truncated synthetic signal','Truncated signal from convolution');
          subplot(414);
          plot(sigstagt(:,ista)-sigconvstagt(:,ista),'k'); hold on
          legend('Difference'); ylim(yran)
          
          %now 'source' contain impulse triplets that contribute to the trace segment of interest, but we
          %want to get ones whose zero-crossing falls into the window.
          source(:,1) = source(:,1)-skiplen;  % index after skipping
          tmp = source(source(:,1)>=indst & source(:,1)<=inded, :);
          tmp = sortrows(tmp,1,'ascend');
          synsrcgtsta{ista} = tmp; %ideally for each station this should be the same, but coincidence is possible
          %     keyboard
        end
        %notify if sources are the same at diff stations
        if ~isequaln(synsrcgtsta{1,1},synsrcgtsta{2,1}) || ~isequaln(synsrcgtsta{1,1},synsrcgtsta{3,1})
          disp('Extracted sources at different stations are not the same, re-check!');
        end
        %   keyboard
      end
      
      %% filter data (and templates)
      % %%filter data
      % optseg = [];
      % for ista = 1:nsta
      %   optseg(:,ista) = Bandpass(STAopt(:,ista), sps, losig, hisig, 2, 2, 'butter');
      % end
      optseg = STAopt;
      
      %%%important, does the data actucally have a similar spectral shape to templates
      if testfreqflag
        for i = 1: nsta
          [f] = plt_wletsig_timefreq(green(:,i),STAopt(:,i),greenf(:,i),optseg(:,i),[],[],sps,...
            [lowlet hiwlet],[losig hisig]);
          text(f.ax(end),0.95, 0.95,stas(i,:),'Units','normalized','HorizontalAlignment','right',...
            'FontSize',12);
        end
      end
      
      %% Best alignment for the whole window
      %some params
      bufsec = 0.5;
      msftaddm = bufsec*sps;  %buffer range for later CC alignment
      rccmwsec = 0.5;
      rccmwlen = rccmwsec*sps;  %window length for computing RCC
      overshoot = rccmwlen/2; %RCC points to the center of computing window, so extra room is needed
      
      %%%obtain a single best alignment based on the entire win
      optcc = detrend(optseg(1+msftaddm: end-msftaddm, :));
      msftadd = 10*sps/40;
      loffmax = 4*sps/40;
      ccmid = ceil(size(optcc,1)/2);
      ccwlen = round(size(optcc,1)-2*(msftadd+1));  % minus ensures successful shifting of records
      ccmin = 0.01;  % depending on the length of trace, cc could be very low
      iup = 1;    % times of upsampling
      [off12con,off13con,ccali,iloopoff,loopoff] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
        ccwlen,msftadd,loffmax,ccmin,iup);
      % if a better alignment cannot be achieved, use 0,0
      if off12con == msftadd+1 && off13con == msftadd+1
        off12con = 0;
        off13con = 0;
        cc12 = xcorr(optcc(:,1), optcc(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
        cc13 = xcorr(optcc(:,1), optcc(:,3),0,'normalized');
        cc23 = xcorr(optcc(:,2), optcc(:,3),0,'normalized');
        ccali = (cc12+cc13+cc23)/3;
        fprintf('Current window cannot be properly aligned, double-check needed \n');
      end
      off1i = zeros(nsta,1);
      off1i(2) = round(off12con);
      off1i(3) = round(off13con);
      
      %%%obtain a single 'observational' alignment between 4th and 1st station, note that this might
      %%%be very different from the empirical prediction from 'empioffset4thsta002'
      %%%--there are also different options here, using sig, or sig-wlet CC as in decon routine
      mcoef = zeros(nsta-3, 1);
      mlag = zeros(nsta-3, 1);
      for ista = 4: nsta
        [mcoef(ista-3),off1i(ista)] = xcorrmax(optcc(:,1), optcc(:,ista), msftadd, 'coeff');
        mlag(ista-3) = off1i(ista);
      end
      
      %%%Align and compute the RCC based on the entire win, and take that as the input signal!
      optdat = [];  % win segment of interest
      for ista = 1: nsta
        optdat(:, ista) = optseg(1+msftaddm-off1i(ista): end-msftaddm-off1i(ista), ista);
      end
      
      %location of the whole-win best alignment
      [loc0, indinput] = off2space002([off1i(2) off1i(3)],sps,ftrans,0);
      % loc has 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
      
      %%%2022/06/06, do NOT taper whatsoever!!
      sigsta = zeros(size(optdat,1), nsta);
      for ista = 1:nsta
        tmp = optdat(:,ista); %best aligned, filtered
        %detrend and taper only the data, NOT the noise
        tmp = detrend(tmp);
        sigsta(:,ista) = tmp;
      end
      %compute running CC between 3 stations
      [ircc,rcc12,rcc13,rcc23] = RunningCC3sta(sigsta,rccmwlen);
      ircc = ircc-overshoot;
      rccpair = [rcc12 rcc13 rcc23];
      %if only use the mean RCC from pair 12 and 13
      rcc = mean(rccpair(:,[1 2]), 2);
      rccsat(:,insat,iperc) = rcc;
      mrcc = median(rcc);
      
      rcc1i = zeros(length(rcc),nsta-3);
      for ista = 4:nsta
        [~,rcc1i(:,ista-3)] = RunningCC(sigsta(:,1), sigsta(:,ista), rccmwlen);
      end
      
      sigsta = detrend(sigsta(overshoot+1:end-overshoot, :));  %excluding the overshoot
      
      %zero-lag max CC
      cc12 = xcorr(sigsta(:,1), sigsta(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
      cc13 = xcorr(sigsta(:,1), sigsta(:,3),0,'normalized');
      cc23 = xcorr(sigsta(:,2), sigsta(:,3),0,'normalized');
      cc1i = zeros(nsta-3,1);
      for ista = 4:nsta
        cc1i(ista-3) = xcorr(sigsta(:,1), sigsta(:,ista),0,'normalized');
      end
      ccpair = [cc12 cc13 cc23];
      mcc = (cc12+cc13)/2;  %if only use the pair 12 and 13
      
      %%
      %if you want a qucik look of the data that feeds to the deconvolution
      if pltdataflag
        figure; hold on
        lsig = size(sigsta,1);
        if nsta == 3
          plot((1:lsig)/sps,sigsta(:,1),'r');
          plot((1:lsig)/sps,sigsta(:,2),'b');
          plot((1:lsig)/sps,sigsta(:,3),'k');
        else
          color = jet(nsta);
          for ista = 1: nsta
            p(ista)=plot((1:lsig)/sps,sigsta(:,ista),'Color',color(ista,:));
            %       label{i}=stas(i,:);
          end
        end
        ax = gca;
        axsym(ax);
        plot((1:lsig)/sps,rcc*ax.YLim(2),'o','color',[.6 .6 .6],'markersize',2);
        text(0.95,0.1,sprintf('Saturation: %.1f',nsat(insat)),'Units','normalized','HorizontalAlignment','right');
        xlim([0 20]);
        legend(p,stas);
      end
      
      %% ground truth of conrtibuting sources
      lsig = size(STAopt,1); %length of original synthetics
      synsrcgtsta = cell(nsta,1);  %sources whose zero-crossing with time range, [indtarvl rnoff12 rnoff13 amp]
      for ista = 1: nsta
        tgreen = zcrosses(ista);  % anchor time of templates, here choose as zero-crossing time
        source = synsrc;
        
        indst = 1+msftaddm+overshoot-off1i(ista);  % starting index of signal in terms of 'STAopt'
        inded = lsig-msftaddm-overshoot-off1i(ista);  % ending index of signal in terms of 'STAopt'
        
        if ista == 1
          greenst = synsrcstind; % the starting index of each added template, context of full length
        elseif ista <=3
          greenst = synsrcstind-source(:, ista)-off1i(ista);
        else
          greenst = pred_tarvl_at4thsta(stas(ista,:),source(:,2),source(:,3),source(:,1),off1i(ista));
        end
        
        %%%you don't need all impulses, only some of them contribute to the length of truncated record
        %%%cut out the green indice and sources that contribute
        %%%Note that some zero-crossings might appear later than 'indedsig', but the corresponding start
        %%%index of the green is actually smaller than 'indedsig', ie. their 'greenst' <= 'indedsig'
        induse = find(greenst>=skiplen-greenlen+indst & greenst<=inded+skiplen-1);
        greenst = greenst(induse);
        source = source(induse,:);
        
        greenzc = greenst+tgreen; % index of approximate zero-crossing
        source(:,1) = source(:,1)+zcrosses(1);  % now 'greenzc' should be the same as 1st col of 'source'
        
        %now 'source' contain impulse triplets that contribute to the trace segment of interest, but we
        %want to get ones whose zero-crossing falls into the window.
        source(:,1) = source(:,1)-skiplen;  % index after skipping
        tmp = source(source(:,1)>=indst & source(:,1)<=inded, :);
        tmp = sortrows(tmp,1,'ascend');
        synsrcgtsta{ista} = tmp; %ideally for each station this should be the same, but coincidence is possible
      end
      % %notify if sources are the same at diff stations
      % if ~isequaln(synsrcgtsta{1,1},synsrcgtsta{2,1}) || ~isequaln(synsrcgtsta{1,1},synsrcgtsta{3,1})
      %   disp('Extracted sources at different stations are not the same, re-check!');
      % else
      ista = 1;
      synsrcgt = synsrcgtsta{ista};
      %store the info of synthetic impulses, same format as in paired deconvolution,
      %9 cols, [ind1 amp1 ind2 amp2 ind3 amp3 off12 off13 off23]
      impgt = zeros(size(synsrcgt,1), 9);
      impgt(:,1) = synsrcgt(:,1)-msftaddm-overshoot;  % cut out the buffer
      impgt(:,[2 4 6]) = repmat(synsrcgt(:,6),1,3);
      %the indices of synthetic impulses need to be shifted too
      impgt(:,3) = impgt(:,1)-synsrcgt(:,2)+off1i(2); % note the sign is consistent!
      impgt(:,5) = impgt(:,1)-synsrcgt(:,3)+off1i(3);
      impgt(:,7) = synsrcgt(:,2);  % no need to shift offset, because they are 'true' offset
      impgt(:,8) = synsrcgt(:,3);
      impgt(:,9) = impgt(:,8)-impgt(:,7);  % off23, == off13 - off12 == tarvl2 - tarvl3
      impgtst = sortrows(impgt,1,'ascend');
      % end
      
      %amp and density at each pixel for ground truth sources
      ampgt = mean(impgt(:,[2 4 6]),2); %amp for all LFE catalog
      density1d = density_pixel(impgt(:,7),impgt(:,8));
      % ampgtsum = sum_pixel(impgt(:,7),impgt(:,8),ampgt);
      [impgtloc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
      density1d = [impgtloc(:,1:2) density1d(:,3)];
      % ampgtsum1d = sortrows([impgtloc(:,1:2) ampgtsum(:,3)], 3);
      
      %%%if you want to plot the ground truth
      if pltgtflag
        %plot the ground truth source offset
        spsscale = sps/40;
        loff_max = 6*spsscale;  % maximum allowed shift is 4 samples at 40 Hz, it is also true from the synthetic distribution
        xran = [-loff_max+off1i(2)-1 loff_max+off1i(2)+1];
        yran = [-loff_max+off1i(3)-1 loff_max+off1i(3)+1];
        offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
        offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
        cran = [0 lsig];
        f.fig = figure;
        f.fig.Renderer = 'painters';
        ax=gca;
        [ax,torispl,mamp] = plt_decon_imp_scatter(ax,impgt,xran,yran,cran,offxran,offyran,...
          sps,50,'mean','tarvl');
        scatter(ax,off1i(2),off1i(3),20,'ks','filled','MarkerEdgeColor','k');
        title(ax,'Ground truth');
        
        %plot the transformed ground truth source map locations
        xran = [-4 4];
        yran = [-4 4];
        cran = [0 lsig/sps];
        f.fig = figure;
        f.fig.Renderer = 'painters';
        ax=gca;
        [ax] = plt_decon_imp_scatter_space(ax,impgt,xran,yran,cran,offxran,...
          offyran,sps,50,ftrans,'mean','tarvl');
        plot(ax,xcut,ycut,'k-','linew',2);
        title(ax,'Ground truth');
        
        %%%plot the cumulative density and summed amp of detections
        [f] = plt_sum_pixel(density1d,ampgtsum1d,[-4 4],[-4 4],20,'amp','linear');
        hold(f.ax(1),'on');
        plot(f.ax(1),xcut,ycut,'k-','linew',2);
        hold(f.ax(2),'on');
        plot(f.ax(2),xcut,ycut,'k-','linew',2);
        supertit(f.ax,'Ground truth');
        
        % %plot the ground truth source map locations from grid
        % f.fig = figure;
        % f.fig.Renderer = 'painters';
        % ax=gca;
        % hold(ax,'on');
        % ax.Box='on'; grid(ax,'on');
        % plot(ax,[xran(1) xran(2)],[0 0],'--','Color',[.3 .3 .3]); hold(ax,'on')
        % plot(ax,[0 0],[yran(1) yran(2)],'--','Color',[.3 .3 .3]);
        % scatter(ax,synsrcgt(:,4),synsrcgt(:,5),synsrcgt(:,6)*30,impgt(:,1)/sps,'filled','MarkerEdgeColor',[.5 .5 .5]);
        % scatter(ax,loc0(:,1),loc0(:,2),30,'k','filled');
        % scatter(ax,shiftor(:,1),shiftor(:,2),30,'k','linew',1);
        % plot(ax,xcut,ycut,'k-','linew',2);
        % oldc = colormap(ax,'kelicol');
        % newc = flipud(oldc);
        % colormap(ax,newc);
        % c=colorbar(ax);
        % caxis(ax, cran);
        % c.Label.String = sprintf('Arrival time at PGC (s)');
        % text(ax,0.98,0.05,sprintf('%d events',size(synsrcgt,1)),'Units','normalized',...
        %   'HorizontalAlignment','right','FontSize',9);
        % ax.YAxis.FontSize = 8;
        % ax.XAxis.FontSize = 8;
        % xlabel(ax,'E (km)','FontSize',11);
        % ylabel(ax,'N (km)','FontSize',11);
        % axis(ax, 'equal');
        % xlim(ax,xran); xticks(ax,xran(1): 1 : xran(2));
        % ylim(ax,yran); yticks(ax,yran(1): 1 : yran(2));
      end
      
      %% independent deconvolution at each station
      %%%finalize the signal, noise, and template (Green's function)
      sigdecon = [];
      pred = [];
      ampit = [];
      for ista = 1:3
        wlet = greenf(:,ista);  %template here is best aligned, tapered, linear trend removed, filtered
        lwlet = length(wlet);
        sig = sigsta(:,ista); %best aligned, filtered, tapered
        lsig = length(sig);
        noi = zeros(lsig,1);
        
        dt = 1/sps;  % sampling interval
        twlet = zcrossesf(ista)*dt;
        width = 2.5;  % width for Gaussian filter
        dres_min = 0.5;  % tolerance, percentage change in residual per iteration
        mfit_min = 5e-1;  % tolerance, norm of the residual/misfit
        tdura = 0.5;  % estimate from the broadband template from fam 002
        tlen = ceil(lsig/sps);
        nit_max = round(1.5*1/tdura*tlen*nsat(insat));  % max numer of iterations
        nimp_max = round(1/tdura*tlen*nsat(insat));%a single peak is ~20 samples wide; maybe a little less (at 100 sps). ~0.4s, 1/0.4=2.5
        fpltit = 0;  % plot flag for each iteration
        fpltend = 0;  % plot flag for the final iteration
        fpltchk = 0; % plot flag for intermediate computations
        
        [sigdecon(:,ista),pred(:,ista),res,dresit,mfitit,ampit{ista},nit,fighdl] = ...
          iterdecon(sig,wlet,rcc,noi,[],dt,twlet,width,dres_min,mfit_min,nit_max,nimp_max,...
          fpltit,fpltend,fpltchk);
        
        if fpltend
          ax = fighdl{2}.ax(1);
          hold(ax,'on');
          text(ax,0.05,0.85,stas(ista,:),'unit','normalized');
          hold(ax,'off');
        end
        
        nit
        
      end
      
      %% Group nearest impulses from different stations into triplets, using moving searching range
      spsscale = sps/40;
      loff_max = 6*spsscale;  % maximum allowed shift is 4 samples at 40 Hz, it is also true from the synthetic distribution
      %note the output 'impindep' gives the arrival index of impulse at each station, after
      %alignment based upon the entire window 'off1i', and the last three cols are the arrival time
      %difference, NOT the true location yet!
      refsta = 1;
      [impindep,imppairf,indpair,sharp] = groupimptripdecon(sigdecon,ampit,rcc,loff_max,refsta);
      
      %note here 'impindepst' inherits the first 6 cols from 'impindep', but the last three cols
      %are adjusted from arrival time difference to the true location offset accounting for the best
      %alignment upon the entire window that is also used in grouping!
      impindep(:,7:8) = impindep(:,7:8)+repmat([off1i(2) off1i(3)],size(impindep,1),1); %account for prealignment
      impindepst = sortrows(impindep,1);
      impgrp{insat,iperc} = impindepst;
      
      %%%if you want to plot the deconvolved sources
      if pltsrcflag1
        %plot the scatter of offsets, accounting for prealignment offset, == true offset
        xran = [-loff_max+off1i(2)-1 loff_max+off1i(2)+1];
        yran = [-loff_max+off1i(3)-1 loff_max+off1i(3)+1];
        offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
        offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
        cran = [0 lsig];
        f.fig = figure;
        f.fig.Renderer = 'painters';
        ax=gca;
        [ax,torispl,mamp] = plt_decon_imp_scatter(ax,impindepst,xran,yran,cran,offxran,offyran,...
          sps,50,'mean','tarvl');
        scatter(ax,off1i(2),off1i(3),20,'ks','filled','MarkerEdgeColor','k');
        title(ax,'Grouped Total');
        
        %plot the scatter of sources in terms of rela locations
        xran = [-4 4];
        yran = [-4 4];
        cran = [0 lsig/sps];
        f.fig = figure;
        f.fig.Renderer = 'painters';
        ax=gca;
        [ax] = plt_decon_imp_scatter_space(ax,impindepst,xran,yran,cran,offxran,...
          offyran,sps,50,ftrans,'mean','tarvl');
        plot(ax,xcut,ycut,'k-','linew',2);
        title(ax,'Grouped Total');
        
        %plot the cumulative density and summed amp of detections
        amp = mean(impindepst(:,[2 4 6]),2); %amp for all LFE catalog
        density1d = density_pixel(impindepst(:,7),impindepst(:,8));
        ampsum = sum_pixel(impindepst(:,7),impindepst(:,8),amp);
        [imploc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
        density1d = [imploc(:,1:2) density1d(:,3)];
        ampsum1d = sortrows([imploc(:,1:2) ampsum(:,3)], 3);
        [f] = plt_sum_pixel(density1d,ampsum1d,[-4 4],[-4 4],20,'amp','linear');
        hold(f.ax(1),'on');
        plot(f.ax(1),xcut,ycut,'k-','linew',2);
        hold(f.ax(2),'on');
        plot(f.ax(2),xcut,ycut,'k-','linew',2);
        supertit(f.ax,'Grouped Total');
        
        %distribution of amp, statistically and in space
        f = initfig(10,5,1,2); %initialize fig
        ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
        plot(ax,(1:length(ampgtsum1d))/length(ampgtsum1d),ampgtsum1d(:,3),'ko-','markersize',3);
        plot(ax,(1:length(ampsum1d))/length(ampsum1d),ampsum1d(:,3),'ro-','markersize',3);
        xlabel(ax,'sorted loc index');
        ylabel(ax,'amp at pixel');
        ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
        histogram(ax,ampgtsum1d(:,3),'facecolor','k');
        histogram(ax,ampsum1d(:,3),'facecolor','r');
        plot(ax,[median(ampsum1d(:,3)) median(ampsum1d(:,3))],ax.YLim,'b--');
        ylabel(ax,'count');
        xlabel(ax,'amp at pixel');
        supertit(f.ax,'Grouped Total');
      end
      % keyboard
      
      %% distance, diff time for source pairs right after grouping
      %note the 'tsep' obtained from the deconvolved positive peaks should be identical to that if
      %obtained from the deconvolved impulses themselves, which represent the arrival indices of the
      %zero-crossing
      ista=1;
      %convert time offset to relative loc
      [imploc, indinput] = off2space002(impindepst(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
      [impindepstst, indsort] = sortrows(impindepst, (ista-1)*2+1);
      implocst = imploc(indsort, :);
      tarvlsplst = impindepstst(:,(ista-1)*2+1);
      
      %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
      m = 5;
      [dtarvl,dneloc,eucdist] = srcdistNtoNm(tarvlsplst, implocst, m);
      
      %project along min-scatter direction
      projang = 135;
      [projx,projy,projxy] = customprojection(implocst,projang);
      dprojxy = [diffcustom(projx,1,'forward') diffcustom(projy,1,'forward')];
      if pltsrcflag1
        if ~isempty(dprojxy)
          nsep = 1;
          [f] = plt_customprojdist(tarvlsplst,implocst,dtarvl{nsep},eucdist{nsep},...
            projxy,dprojxy,projang,sps,'tarvl');
          hold(f.ax(1),'on');
          plot(f.ax(1),xcut,ycut,'k-','linew',2);
        end
      end
      mprojx1nn1(insat,iperc) = median(abs(dprojxy(:,1))); %median consecutive dist along min-scatter after grouping
      
      %For each LFE source, get its distance to all other LFEs, maybe we don't care that long separation in time
      [~,dprojxy2all] = srcdistall(tarvlsplst,projxy,[0 2*sps]);
      mprojx12all(insat,iperc) = median(abs(dprojxy2all(:,1))); %median dist along min-scatter to all sources after grouping
      
      %% Remove the small-amplitude, secondary triplets from the grouped result
      %convert the sources in terms of the arrival time of zero-crossing to positve peaks' indices
      ppkindep = impindep;
      for ista = 1: 3
        ppkindep(:,(ista-1)*2+1) = ppkindep(:,(ista-1)*2+1)+ppeaksf(ista)-zcrossesf(ista);
      end
      npkindep = impindep;  %negative peaks
      for ista = 1: 3
        npkindep(:,(ista-1)*2+1) = npkindep(:,(ista-1)*2+1)+npeaksf(ista)-zcrossesf(ista);
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
      sharp(indremove, :) = [];
      impindepst = sortrows(impindep,1);
      imp{insat,iperc} = impindepst;
      
      %%%if you want to plot the deconvolved sources
      if pltsrcflag2
        %plot the scatter of offsets, accounting for prealignment offset, == true offset
        xran = [-loff_max+off1i(2)-1 loff_max+off1i(2)+1];
        yran = [-loff_max+off1i(3)-1 loff_max+off1i(3)+1];
        offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
        offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
        cran = [0 lsig];
        f.fig = figure;
        f.fig.Renderer = 'painters';
        ax=gca;
        [ax,torispl,mamp,xbnd,ybnd] = plt_decon_imp_scatter(ax,impindepst,xran,yran,cran,offxran,offyran,...
          sps,50,'mean','tarvl');
        scatter(ax,off1i(2),off1i(3),20,'ks','filled','MarkerEdgeColor','k');
        title(ax,'Secondary sources removed');
        
        %plot the scatter of sources in terms of rela locations
        xran = [-4 4];
        yran = [-4 4];
        cran = [0 lsig/sps];
        f.fig = figure;
        f.fig.Renderer = 'painters';
        ax=gca;
        [ax] = plt_decon_imp_scatter_space(ax,impindepst,xran,yran,cran,offxran,...
          offyran,sps,50,ftrans,'mean','tarvl');
        plot(ax,xcut,ycut,'k-','linew',2);
        title(ax,'Secondary sources removed');
        
        %plot the cumulative density and summed amp of detections
        amp = mean(impindepst(:,[2 4 6]),2); %amp for all LFE catalog
        density1d = density_pixel(impindepst(:,7),impindepst(:,8));
        ampsum = sum_pixel(impindepst(:,7),impindepst(:,8),amp);
        [imploc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
        density1d = [imploc(:,1:2) density1d(:,3)];
        ampsum1d = sortrows([imploc(:,1:2) ampsum(:,3)], 3);
        [f] = plt_sum_pixel(density1d,ampsum1d,[-4 4],[-4 4],20,'amp','linear');
        hold(f.ax(1),'on');
        plot(f.ax(1),xcut,ycut,'k-','linew',2);
        hold(f.ax(2),'on');
        plot(f.ax(2),xcut,ycut,'k-','linew',2);
        supertit(f.ax,'Secondary sources removed');
        
        %distribution of amp, statistically and in space
        f = initfig(10,5,1,2); %initialize fig
        ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
        plot(ax,(1:length(ampgtsum1d))/length(ampgtsum1d),ampgtsum1d(:,3),'ko-','markersize',3);
        plot(ax,(1:length(ampsum1d))/length(ampsum1d),ampsum1d(:,3),'ro-','markersize',3);
        xlabel(ax,'sorted loc index');
        ylabel(ax,'amp at pixel');
        ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
        histogram(ax,ampgtsum1d(:,3),'facecolor','k');
        histogram(ax,ampsum1d(:,3),'facecolor','r');
        plot(ax,[median(ampsum1d(:,3)) median(ampsum1d(:,3))],ax.YLim,'b--');
        ylabel(ax,'count');
        xlabel(ax,'amp at pixel');
        supertit(f.ax,'Secondary sources removed');
      end
      % keyboard
      
      %% distance, diff time for source pairs, after 2ndary sources removed
      %note the 'tsep' obtained from the deconvolved positive peaks should be identical to that if
      %obtained from the deconvolved impulses themselves, which represent the arrival indices of the
      %zero-crossing
      ista=1;
      %convert time offset to relative loc
      [imploc, indinput] = off2space002(impindepst(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
      [impindepstst, indsort] = sortrows(impindepst, (ista-1)*2+1);
      implocst = imploc(indsort, :);
      tarvlsplst = impindepstst(:,(ista-1)*2+1);
      
      %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
      m = 5;
      [dtarvl,dneloc,eucdist] = srcdistNtoNm(tarvlsplst, implocst, m);
      
      %project along min-scatter direction
      projang = 135;
      [projx,projy,projxy] = customprojection(implocst,projang);
      dprojxy = [diffcustom(projx,1,'forward') diffcustom(projy,1,'forward')];
      if pltsrcflag2
        if ~isempty(dprojxy)
          nsep = 1;
          [f] = plt_customprojdist(tarvlsplst,implocst,dtarvl{nsep},eucdist{nsep},...
            projxy,dprojxy,projang,sps,'tarvl');
          hold(f.ax(1),'on');
          plot(f.ax(1),xcut,ycut,'k-','linew',2);
        end
      end
      mprojx2nn1(insat,iperc) = median(abs(dprojxy(:,1))); %median consecutive dist along min-scatter after 2nd src removal
      
      %For each LFE source, get its distance to all other LFEs, maybe we don't care that long separation in time
      [~,dprojxy2all] = srcdistall(tarvlsplst,projxy,[0 2*sps]);
      mprojx22all(insat,iperc) = median(abs(dprojxy2all(:,1))); %median dist along min-scatter to all sources after 2nd src removal
      
      %% plot amp ratio 12 and 13, and 23, without 2ndary sources
      %%%ideally, when templates and data are not normalized, and there is no particular noise or any
      %%%other factors causing the amplitude scaling between temp and data and each station to be
      %%%vastly different, then for each deconvolved source, the direct impulse amp should be
      %%%~identical at all stations, ie., the ratio between station pairs should be ~1
      srcampr = [impindepst(:,2)./impindepst(:,4) impindepst(:,2)./impindepst(:,6) ...
        impindepst(:,4)./impindepst(:,6)];
      
      psrcamprs = [impindepst(:,2)*max(greenf(:,1))./(impindepst(:,4)*max(greenf(:,2))) ...
        impindepst(:,2)*max(greenf(:,1))./(impindepst(:,6)*max(greenf(:,3))) ...
        impindepst(:,4)*max(greenf(:,2))./(impindepst(:,6)*max(greenf(:,3)))];
      psrcamps = [impindepst(:,2)*max(greenf(:,1)) impindepst(:,4)*max(greenf(:,2)) ...
        impindepst(:,6)*max(greenf(:,3))];
      
      nsrcamprs = [impindepst(:,2)*min(greenf(:,1))./(impindepst(:,4)*min(greenf(:,2))) ...
        impindepst(:,2)*min(greenf(:,1))./(impindepst(:,6)*min(greenf(:,3))) ...
        impindepst(:,4)*min(greenf(:,2))./(impindepst(:,6)*min(greenf(:,3)))];
      nsrcamps = [impindepst(:,2)*min(greenf(:,1)) impindepst(:,4)*min(greenf(:,2))...
        impindepst(:,6)*min(greenf(:,3))];
      
      msrcampr{insat,iperc} = median(log10(srcampr), 1);
      madsrcampr{insat,iperc} = mad(log10(srcampr), 1, 1);
      mpsrcamprs{insat,iperc} = median(log10(psrcamprs), 1);
      madpsrcamprs{insat,iperc} = mad(log10(psrcamprs), 1, 1);
      mnsrcamprs{insat,iperc} = median(log10(nsrcamprs), 1);
      madnsrcamprs{insat,iperc} = mad(log10(nsrcamprs), 1, 1);
      nsrc{insat,iperc} = size(srcampr,1);
      
      %what is the deviation of amp ratio from the median for each source?
      lndevsrcampr = srcampr-median(srcampr, 1); % in linear scale
      lgdevsrcampr = log10(srcampr)-median(log10(srcampr), 1); % in log scale, note that log2+log5=log10, so this means a ratio
      
      srcamprall{insat,iperc} = srcampr;
      
      %plot histograms of source amp
      if pltsrcflag2
        f = initfig(12,6,2,3); %initialize fig
        f=plt_deconpk_rat_comb(f,srcampr,impindepst,'k');
        supertit(f.ax,'Secondary sources removed');
      end
      % keyboard
      
      %% 2ndary src removed, prediction of impulse tarvl at 4th sta given sources and empirical off14-src relation
      %%%carry out 'deconvolution' at 4th stations as well for the tarvl and amp
      modname = 'timeoff_plfit_4thsta_160sps.mat';
      planefit = load(strcat(rstpath, '/MAPS/',modname));
      
      pred4diff = [];
      predoff14 = [];
      for ista = 4:nsta
        wlet = greenf(:,ista);  %template here is best aligned, tapered, linear trend removed, filtered
        lwlet = length(wlet);
        sig = sigsta(:,ista); %best aligned, filtered, tapered
        lsig = length(sig);
        
        dt = 1/sps;  % sampling interval
        twlet = zcrossesf(ista)*dt;
        fpltit = 0;  % plot flag for each iteration
        fpltend = 0;  % plot flag for the final iteration
        fpltchk = 0; % plot flag for intermediate computations
        rmse = planefit.gof{ista}.rmse;
        %1.6*rmse seems proper from whole-win RCC, 2 from concatenated RCC
        offmax = round(2.0*rmse);
        
        [sigdecon(:,ista),pred,res,dresit,mfitit,ampit{ista},fighdl] = ...
          iterdecon_4thsta(sig,wlet,[],rcc1i(:,ista-3),[],...
          dt,twlet,impindep,stas(ista,:),off1i(ista),[],offmax,...
          fpltit,fpltend,fpltchk);
        
        ampiti = ampit{ista};
        impindep(:,9+(ista-4)*2+1) = ampiti(:,1);
        impindep(:,9+(ista-3)*2) = ampiti(:,2);
        if ista == nsta
          pred4diff(:,ista-3) = ampiti(:,end-1);  %difference between found peak and predicted arrival
          predoff14(:,ista-3) = ampiti(:,end);  %predicted off14
        end
      end
      
      %%%given zero-crossing indices, obtain corresponding positive and negative peak indices
      ppkindep = impindep;
      npkindep = impindep;
      for ista = 1: 3
        ppkindep(:,(ista-1)*2+1) = ppkindep(:,(ista-1)*2+1)+ppeaksf(ista)-zcrossesf(ista);
        npkindep(:,(ista-1)*2+1) = npkindep(:,(ista-1)*2+1)+npeaksf(ista)-zcrossesf(ista);
      end
      for ista = 4:nsta
        ppkindep(:,9+(ista-4)*2+1) = ppkindep(:,9+(ista-4)*2+1)+ppeaksf(ista)-zcrossesf(ista);
        npkindep(:,9+(ista-4)*2+1) = npkindep(:,9+(ista-4)*2+1)+npeaksf(ista)-zcrossesf(ista);
      end
      
      %       %%%plot the predicted tarvl and amp of decon impulses at 4th sta vs. waveform
      %       [f] = plt_deconpk_sigpk_comp_4thsta(sigsta(:,4:nsta),stas(4:nsta,:),...
      %         impindep(:,10:end),ppkindep4th,npkindep4th,greenf(:,4:nsta));
      
      %%% further ELIMINATE sources that fail the check at 4th stations
      trust4th = 4; % trust KLNB the most among all 4th stations
      indremove = find(impindep(:,9+(trust4th-4)*2+1)==0 & impindep(:,9+(trust4th-3)*2)==0);
      pred4difftr = pred4diff(setdiff(1:size(pred4diff,1),indremove),trust4th-3);
      predoff14tr = predoff14(setdiff(1:size(pred4diff,1),indremove),trust4th-3);
      impindep(indremove,:) = [];
      ppkindep(indremove, :) = [];
      npkindep(indremove, :) = [];
      impindepst = sortrows(impindep,1);
      imp4th{insat,iperc} = impindepst;
      
      %%%if you want to plot the deconvolved sources
      if pltsrcflag3
        %plot the scatter of offsets, accounting for prealignment offset, == true offset
        xran = [-loff_max+off1i(2)-1 loff_max+off1i(2)+1];
        yran = [-loff_max+off1i(3)-1 loff_max+off1i(3)+1];
        offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
        offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
        cran = [0 lsig];
        f.fig = figure;
        f.fig.Renderer = 'painters';
        ax=gca;
        [ax,torispl,mamp,xbnd,ybnd] = plt_decon_imp_scatter(ax,impindepst,xran,yran,cran,offxran,offyran,...
          sps,50,'mean','tarvl');
        scatter(ax,off1i(2),off1i(3),20,'ks','filled','MarkerEdgeColor','k');
        title(ax,'Checkd at 4th stas');
        
        %plot the scatter of sources in terms of rela locations
        xran = [-4 4];
        yran = [-4 4];
        cran = [0 lsig/sps];
        f.fig = figure;
        f.fig.Renderer = 'painters';
        ax=gca;
        [ax] = plt_decon_imp_scatter_space(ax,impindepst,xran,yran,cran,offxran,...
          offyran,sps,50,ftrans,'mean','tarvl');
        plot(ax,xcut,ycut,'k-','linew',2);
        title(ax,'Checkd at 4th stas');
        
        %plot the cumulative density and summed amp of detections
        amp = mean(impindepst(:,[2 4 6]),2); %amp for all LFE catalog
        density1d = density_pixel(impindepst(:,7),impindepst(:,8));
        ampsum = sum_pixel(impindepst(:,7),impindepst(:,8),amp);
        [imploc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
        density1d = [imploc(:,1:2) density1d(:,3)];
        ampsum1d = sortrows([imploc(:,1:2) ampsum(:,3)], 3);
        [f] = plt_sum_pixel(density1d,ampsum1d,[-4 4],[-4 4],20,'amp','linear');
        hold(f.ax(1),'on');
        plot(f.ax(1),xcut,ycut,'k-','linew',2);
        hold(f.ax(2),'on');
        plot(f.ax(2),xcut,ycut,'k-','linew',2);
        supertit(f.ax,'Checkd at 4th stas');
        
      end
      
      %%%
      %off14 prediction for decon srcs using plane fit model WITH alignment
      [~,off14pred] = pred_tarvl_at4thsta(stas(trust4th,:),impindepst(:,7),impindepst(:,8),...
        impindepst(:,1),0);
      %off14 computed from actually-matched arrivals at stas 1 and 4 from decon
      % +repmat(off1i(trust4th),size(impindepst,1),1)
      off14 = impindepst(:,1)-impindepst(:,9+(trust4th-4)*2+1); %after prealignment
      
      %source arrival prediction from plane fitting, including calibrating the arrival prediction in the
      %context of prealigned signal, note sign is '+'
      [~,off14gt] = pred_tarvl_at4thsta(stas(trust4th,:),impgtst(:,7),impgtst(:,8),impgtst(:,1),0);
      
      if pltsrcflag3
        % f = initfig(5,5,1,1); %initialize fig
        % ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
        % p1=scatter3(ax,impindepst(:,7),impindepst(:,8),off14pred,15,'r','filled','MarkerEdgeColor','k');
        % p2=scatter3(ax,impgtst(:,7),impgtst(:,8),off14gt,15,[.5 .5 .5],'filled','MarkerEdgeColor','k');
        % legend(ax,[p1 p2],'plane fit pred for decon srcs', 'plane fit pred for ground truth');
        % xlabel(ax,sprintf('off12 at %d sps',sps));
        % ylabel(ax,sprintf('off13 at %d sps',sps));
        % zlabel(ax,sprintf('off14 at %d sps',sps));
        % title(ax,stas(trust4th,:));
        % view(ax, 35, 27);
        % ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
        % p1=scatter3(ax,impindepst(:,7),impindepst(:,8),off14pred,15,'r','filled','MarkerEdgeColor','k');
        % p2=scatter3(ax,impindepst(:,7),impindepst(:,8),off14,15,[.5 .5 .5],'filled',...
        %   'MarkerEdgeColor','k'); %correct ground truth for record alignment
        % legend(ax,[p1 p2],'plane fit pred for decon srcs', 'actually-matched arrival diff.');
        % xlabel(ax,sprintf('off12 at %d sps',sps));
        % ylabel(ax,sprintf('off13 at %d sps',sps));
        % zlabel(ax,sprintf('off14 at %d sps',sps));
        % view(ax, 45, 10);
        
        %plot error in off14
        f = initfig(12,5,1,2); %initialize fig
        ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
        histogram(ax,pred4difftr);
        plot(ax,[median(pred4difftr) median(pred4difftr)],ax.YLim,'r--','LineWidth',1);
        plot(ax,[-offmax -offmax],ax.YLim,'k--');
        plot(ax,[offmax offmax],ax.YLim,'k--');
        xlabel(ax,'diff. in 4th arrival between pred and decon');
        ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
        histogram(ax,off14pred-off14);
        plot(ax,[median(off14pred-off14) median(off14pred-off14)],ax.YLim,'r--','LineWidth',1);
        plot(ax,[-offmax -offmax],ax.YLim,'k--');
        plot(ax,[offmax offmax],ax.YLim,'k--');
        xlabel(ax,'diff. in off14 between plane-fit and decon');
        ylabel(ax,'count');
        % ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
        % histogram(ax,off14pred-off14gt);
        % plot(ax,[-offmax -offmax],ax.YLim,'k--');
        % plot(ax,[offmax offmax],ax.YLim,'k--');
        % xlabel(ax,'diff. in plane-fit off14 between decon and ground truth');
        % ylabel(ax,'count');
      end
      
      % keyboard
      
      %%%amp ratio 12 and 13, and 23, with 4th station checked
      %%%ideally, when templates and data are not normalized, and there is no particular noise or any
      %%%other factors causing the amplitude scaling between temp and data and each station to be
      %%%vastly different, then for each deconvolved source, the direct impulse amp should be
      %%%~identical at all stations, ie., the ratio between station pairs should be ~1
      srcampr4th = [impindepst(:,2)./impindepst(:,4) impindepst(:,2)./impindepst(:,6) ...
        impindepst(:,4)./impindepst(:,6) impindepst(:,2)./impindepst(:,end)];
      
      psrcamprs4th = [impindepst(:,2)*max(greenf(:,1))./(impindepst(:,4)*max(greenf(:,2))) ...
        impindepst(:,2)*max(greenf(:,1))./(impindepst(:,6)*max(greenf(:,3))) ...
        impindepst(:,4)*max(greenf(:,2))./(impindepst(:,6)*max(greenf(:,3))) ...
        impindepst(:,2)*max(greenf(:,1))./(impindepst(:,end)*max(greenf(:,end)))];
      psrcamps4th = [impindepst(:,2)*max(greenf(:,1)) impindepst(:,4)*max(greenf(:,2)) ...
        impindepst(:,6)*max(greenf(:,3)) impindepst(:,end)*max(greenf(:,end))];
      
      nsrcamprs4th = [impindepst(:,2)*min(greenf(:,1))./(impindepst(:,4)*min(greenf(:,2))) ...
        impindepst(:,2)*min(greenf(:,1))./(impindepst(:,6)*min(greenf(:,3))) ...
        impindepst(:,4)*min(greenf(:,2))./(impindepst(:,6)*min(greenf(:,3))) ...
        impindepst(:,2)*min(greenf(:,1))./(impindepst(:,end)*min(greenf(:,end)))];
      nsrcamps4th = [impindepst(:,2)*min(greenf(:,1)) impindepst(:,4)*min(greenf(:,2))...
        impindepst(:,6)*min(greenf(:,3)) impindepst(:,end)*min(greenf(:,end))];
      
      msrcampr4th{insat,iperc} = median(log10(srcampr4th), 1);
      madsrcampr4th{insat,iperc} = mad(log10(srcampr4th), 1, 1);
      mpsrcamprs4th{insat,iperc} = median(log10(psrcamprs4th), 1);
      madpsrcamprs4th{insat,iperc} = mad(log10(psrcamprs4th), 1, 1);
      mnsrcamprs4th{insat,iperc} = median(log10(nsrcamprs4th), 1);
      madnsrcamprs4th{insat,iperc} = mad(log10(nsrcamprs4th), 1, 1);
      nsrc4th{insat,iperc} = size(srcampr4th,1);
      
      %what is the deviation of amp ratio from the median for each source?
      lndevsrcampr = srcampr-median(srcampr, 1); % in linear scale
      lgdevsrcampr = log10(srcampr)-median(log10(srcampr), 1); % in log scale, note that log2+log5=log10, so this means a ratio
      
      srcamprall4th{insat,iperc} = srcampr4th;
      
      if pltsrcflag3
        %plot histograms of source amp
        f = initfig(16,6,2,4); %initialize fig
        f = plt_deconpk_rat_comb4th(f,srcampr4th,impindepst,'k');
        supertit(f.ax,'Checkd at 4th stas');
      end
      % keyboard
      
      %% distance, diff time for source pairs, after 4th sta check
      %note the 'tsep' obtained from the deconvolved positive peaks should be identical to that if
      %obtained from the deconvolved impulses themselves, which represent the arrival indices of the
      %zero-crossing
      ista=1;
      %convert time offset to relative loc
      [imploc, indinput] = off2space002(impindepst(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
      [impindepstst, indsort] = sortrows(impindepst, (ista-1)*2+1);
      implocst = imploc(indsort, :);
      tarvlsplst = impindepstst(:,(ista-1)*2+1);
      
      %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
      m = 5;
      [dtarvl,dneloc,eucdist] = srcdistNtoNm(tarvlsplst, implocst, m);
      
      %project along min-scatter direction
      projang = 135;
      [projx,projy,projxy] = customprojection(implocst,projang);
      dprojxy = [diffcustom(projx,1,'forward') diffcustom(projy,1,'forward')];
      if pltsrcflag3
        if ~isempty(dprojxy)
          nsep = 1;
          [f] = plt_customprojdist(tarvlsplst,implocst,dtarvl{nsep},eucdist{nsep},...
            projxy,dprojxy,projang,sps,'tarvl');
          hold(f.ax(1),'on');
          plot(f.ax(1),xcut,ycut,'k-','linew',2);
        end
      end
      mprojx3nn1(insat,iperc) = median(abs(dprojxy(:,1))); %median consecutive dist along min-scatter after 2nd src removal
      
      %For each LFE source, get its distance to all other LFEs, maybe we don't care that long separation in time
      [~,dprojxy2all] = srcdistall(tarvlsplst,projxy,[0 2*sps]);
      mprojx32all(insat,iperc) = median(abs(dprojxy2all(:,1))); %median dist along min-scatter to all sources after 2nd src removal
      
      
    end %loop end for saturation level
    
    % keyboard
    
  end %loop end for percent of noise
  
  save('rst_synth_onespot.mat');
  
end %if need to recalculate
keyboard

%%
%%%loop for noise level
for iperc = 1: ntrial
  perc = perctrial(iperc);
  disp(perc);
  
  %initialize figs
  f1 = initfig(16,8,2,4,(iperc-1)*2+1); %map loc for srcs no 2ndary srcs
  supertit(f1.ax,sprintf('noise: %.1f, Secondary sources removed',perc));
  f2 = initfig(16,8,2,4,iperc*2); %map loc for srcs checked at 4th stas
  supertit(f2.ax,sprintf('noise: %.1f, Checkd at 4th stas',perc));
  
  %%%loop for saturation level
  for insat = 1: nnsat
    
    impindepst = imp{insat,iperc};
    
    %plot the scatter of sources in terms of rela locations
    % xran = [-4 4];
    % yran = [-4 4];
    % cran = [0 lsig/sps];
    % offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
    % offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
    % ax=f1.ax(insat);
    % [ax] = plt_decon_imp_scatter_space(ax,impindepst,xran,yran,cran,offxran,...
    %   offyran,sps,30,ftrans,'mean','tarvl');
    % plot(ax,xcut,ycut,'k-','linew',2);
    % text(ax,0.98,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
    %       'HorizontalAlignment','right');
    % text(ax,0.02,0.05,sprintf('%.2f',mprojx22all(insat,iperc)),'Units','normalized',...
    % 'HorizontalAlignment','left');
    
    %plot the cumulative density and summed amp of detections
    density1d = density_pixel(impindepst(:,7),impindepst(:,8));
    [imploc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
    density1d = [imploc(:,1:2) density1d(:,3)];
    xran = [-4 4];
    yran = [-4 4];
    scale = 'linear';
    ax=f1.ax(insat);
    hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
    dum = density1d;
    dum(dum(:,3)>1, :) = [];
    if strcmp(scale,'log')
      dum(:,3) = log10(dum(:,3));
    end
    scatter(ax,dum(:,1),dum(:,2),10,dum(:,3),'linew',0.2);  %, 'MarkerEdgeColor', 'w')
    dum = sortrows(density1d,3);
    dum(dum(:,3)==1, :) = [];
    if strcmp(scale,'log')
      dum(:,3) = log10(dum(:,3));
    end
    scatter(ax,dum(:,1),dum(:,2),10,dum(:,3),'filled','MarkerEdgeColor','none');
    colormap(ax,flipud(colormap(ax,'kelicol')));
    c=colorbar(ax,'SouthOutside');
    ax.CLim(2) = prctile(dum(:,3),99);
    if strcmp(scale,'log')
      c.Label.String = strcat('log_{10}(# of detections)');
    elseif strcmp(scale,'linear')
      c.Label.String = '# of detections';
    end
    axis(ax,'equal');
    axis(ax,[xran yran]);
    ax.GridLineStyle = '--';
    ax.XAxisLocation = 'top';
    text(ax,0.98,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
      'HorizontalAlignment','right');
    text(ax,0.98,0.05,sprintf('%d events',size(impindepst,1)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.02,0.05,sprintf('%.2f',mprojx22all(insat,iperc)),'Units','normalized',...
      'HorizontalAlignment','left');
    
    
    impindepst = imp4th{insat,iperc};
    
    % %plot the scatter of sources in terms of rela locations
    % xran = [-4 4];
    % yran = [-4 4];
    % cran = [0 lsig/sps];
    % offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
    % offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
    % ax=f2.ax(insat);
    % [ax] = plt_decon_imp_scatter_space(ax,impindepst,xran,yran,cran,offxran,...
    %   offyran,sps,30,ftrans,'mean','tarvl');
    % plot(ax,xcut,ycut,'k-','linew',2);
    % text(ax,0.98,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
    %       'HorizontalAlignment','right');
    
    %plot the cumulative density and summed amp of detections
    density1d = density_pixel(impindepst(:,7),impindepst(:,8));
    [imploc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
    density1d = [imploc(:,1:2) density1d(:,3)];
    xran = [-4 4];
    yran = [-4 4];
    ax=f2.ax(insat);
    hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
    dum = density1d;
    dum(dum(:,3)>1, :) = [];
    if strcmp(scale,'log')
      dum(:,3) = log10(dum(:,3));
    end
    scatter(ax,dum(:,1),dum(:,2),10,dum(:,3),'linew',0.2);  %, 'MarkerEdgeColor', 'w')
    dum = sortrows(density1d,3);
    dum(dum(:,3)==1, :) = [];
    if strcmp(scale,'log')
      dum(:,3) = log10(dum(:,3));
    end
    scatter(ax,dum(:,1),dum(:,2),10,dum(:,3),'filled','MarkerEdgeColor','none');
    colormap(ax,flipud(colormap(ax,'kelicol')));
    c=colorbar(ax,'SouthOutside');
    ax.CLim(2) = prctile(dum(:,3),99);
    if strcmp(scale,'log')
      c.Label.String = strcat('log_{10}(# of detections)');
    elseif strcmp(scale,'linear')
      c.Label.String = '# of detections';
    end
    axis(ax,'equal');
    axis(ax,[xran yran]);
    ax.GridLineStyle = '--';
    ax.XAxisLocation = 'top';
    text(ax,0.98,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
      'HorizontalAlignment','right');
    text(ax,0.98,0.05,sprintf('%d events',size(impindepst,1)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.02,0.05,sprintf('%.2f',mprojx32all(insat,iperc)),'Units','normalized',...
      'HorizontalAlignment','left');
    
  end %loop end for saturation level
end %loop end for noise level

%% Summary of amplitude ratio
%%%loop for noise level
for iperc = 1: ntrial
  perc = perctrial(iperc);
  disp(perc);
  
  %%%loop for saturation level
  for insat = 1: nnsat
    disp(nsat(insat));
    
    f3 = initfig(16,8,2,4); %plot histograms of source amp
    supertit(f3.ax,sprintf('noise: %.1f, Secondary sources removed & Checkd at 4th stas',...
      perc));
    %     f4 = initfig(16,4,1,4); %plot histograms of source amp
    %     supertit(f4.ax,'Checkd at 4th stas');
    
    impindepst = imp{insat,iperc};
    srcampr = srcamprall{insat,iperc};
    f3.ax(1:3) = plt_deconpk_rat_comb(f3.ax(1:3),srcampr,impindepst,'k','hist');
    
    impindepst = imp4th{insat,iperc};
    srcampr4th = srcamprall4th{insat,iperc};
    f3.ax(5:end) = plt_deconpk_rat_comb4th(f3.ax(5:end),srcampr4th,impindepst,'k','hist');
    
    %     keyboard
  end %loop end for saturation level
end %loop end for noise level

%% consecutive dist along min-scatter VS saturation rate & region size
f = initfig(12,4,1,3); %initialize fig
color = jet(ntrial);
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
for iperc = 1: ntrial
  p(iperc) = plot(ax,log10(nsat),mprojx12all(:,iperc),'-o','markersize',4,'color',color(iperc,:));
  label{iperc} = sprintf('noise=%.1f',perctrial(iperc));
end
legend(ax,p,label);
title(ax,'Grouped Total');
xlabel(ax,'Saturation level (log)');
ylabel(ax,'med. consec. dist. along min-scatter direc');

ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
for iperc = 1: ntrial
  plot(ax,log10(nsat),mprojx22all(:,iperc),'-o','markersize',4,'color',color(iperc,:));
end
title(ax,'Secondary sources removed');

ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
for iperc = 1: ntrial
  plot(ax,log10(nsat),mprojx32all(:,iperc),'-o','markersize',4,'color',color(iperc,:));
end
title(ax,'Checkd at 4th stas');


