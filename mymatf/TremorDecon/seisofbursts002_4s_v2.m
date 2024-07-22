% seisofbursts002_4s_v2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script now aims to compute the rcc, renv, ramp in a lower passband 
% than the high-frequency one, 1.8-6.3 Hz, which is used for filtering data
% for the purpose of deconvolution. 
% --The way of computation is similar to 25-s-win
% in 'deconv_ref_4s_exp_4thsta_fn.m'. The segmentation, alignment are the 
% same as HF. However, the length of running window is longer than 0.5 s.
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
% flagrecalc = 0;
flagrecalc = 1;

if flagrecalc

  defval('noiseflag',0);  %whether to use synthetic noises

  rccflag = 1; %1 means RCC weighting is used
  whichrcc = 0; %if rcc weighting, which pair is used, 0 is best 2 pairs; 1 is 12; 2 is 13; 3 is 23

  %Choice to make upon the actual-used alignment at 4th stations
  if noiseflag
    align14flag = 0;  %do NOT align sta 4 wrt. 1 if using noise
  else
    align14flag = 1; 
  end

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
  ntol = 3;
  trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.',...
    num2str(ntol),'.pgc002.',cutout(1:4)));
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

  ftrans = 'interpchao';

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

  %% load HF decon results for the alignment of 25-s windows
  % savefile = 'deconv_stats4th_allbstsig.mat';
  savefile = 'deconv_stats4th_no23_allbstsig.mat';
  load(strcat(rstpath, '/MAPS/',savefile));

  nsrc = allbstsig.nsrc;
  imp = allbstsig.impindepall;
  nsrc4th = allbstsig.nsrc4th;
  imp4th = allbstsig.impindep4thall;
  off1i = allbstsig.off1i;
  off1iwk = allbstsig.off1iwk;

  %% read daily data, break into windows of segments, plot
  %filtering passband for reading data
%   hisig=1.8; % this will give a similar spectral shape between template and signal
%   losig=0.5;

  hisig=6.3; % this will give a similar spectral shape between template and signal
  losig=1.8;

  sps = 160;

  %moving window length in samples for running CC, envelope, etc.
  if hisig==1.8
    rccmwsec=4;
  elseif hisig==6.3
    rccmwsec=0.5;
  end
  rccmwlen=rccmwsec*sps;
  rccmwlenhf=0.5*sps;

  idxbst = 1:size(trange,1); 
%   idxbst = 57;
  
  k = 0;  % counting the burst windows

  irccrank = cell(size(trange,1),1); %indices of rcc concatenation 
  ccwpairk = cell(size(trange,1),1); % 0-lag CC value for all subwins and all burst wins
  mrccwpairk = cell(size(trange,1),1); % median RCC value for all subwins and all burst wins
  renvcatk = cell(size(trange,1),1); %median of running envelope of same win as RCC, mean of 3 stas
  renvcatnormk = cell(size(trange,1),1); %running envelope normalized by median of the burst
  rampcatk = cell(size(trange,1),1); %median of running amplitude, similar to 'renvcatk'
  rampcatnormk = cell(size(trange,1),1); %running amplitude normalized by median of the burst
  
  rcccatsrcall = [];  %mean concat RCC among trio and RCC14 at src arrival 
  rccpairsrcall = []; %concat RCC for trio sta pairs at src arrival
  renvcatsrcall = [];  %mean concat RENV among trio at src arrival
  renvcatnormsrcall = []; %normalized
  renvpairsrcall = []; %concat RENV for each trio sta at src arrival
  rampcatsrcall = [];  %mean concat RAMP among trio at src arrival
  rampcatnormsrcall = []; %normalized
  ramppairsrcall = []; %concat RAMP for each trio sta at src arrival
  
  rcccatsrc1all = []; %concat RCC at src arrival at sta 1
  renvcatsrc1all = [];  
  renvcatnormsrc1all = [];
  diffzcsrc1all = [];
  rcccatsrczc1all = [];
  renvcatsrczc1all = [];
  renvcatnormsrczc1all = [];
        
  renvortcatsrcall = [];  %mean concat RENV among trio at src arrival
  renvortcatnormsrcall = []; %normalized
  renvpairortsrcall = []; %concat RENV for each trio sta at src arrival
  renvvertcatsrcall = [];  %mean concat RENV among trio at src arrival
  renvvertcatnormsrcall = []; %normalized
  renvpairvertsrcall = []; %concat RENV for each trio sta at src arrival

  rcccatsrc4thall = [];  %mean concat RCC among trio and RCC14 at src arrival 
  rccpairsrc4thall = []; %concat RCC for trio sta pairs at src arrival 
  renvcatsrc4thall = [];  %mean concat RENV among trio at src arrival 
  renvcatnormsrc4thall = []; %normalized
  renvpairsrc4thall = []; %concat RENV for each trio sta at src arrival 
  rampcatsrc4thall = [];  %mean concat RAMP among trio at src arrival 
  rampcatnormsrc4thall = []; %normalized
  ramppairsrc4thall = []; %concat RAMP for each trio sta at src arrival 
  rcccatsrc14thall = [];
  renvcatsrc14thall = [];  
  renvcatnormsrc14thall = [];

  for iii = 1: length(idxbst)
  
    [iets,i,j] = indofburst(trange,idxbst(iii));
  
%   for iets = 1: nets
    % dates in each ets
    year = years(iets);
    datesets = dates(floor(dates/1000)==year);
        
%     for i = 1: length(datesets)
      
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
      
      %bursts and 4-s detections of the same day
      rangetemp = trange(trange(:,1)==datesets(i), :);
      hfdayi = hfbnd(hfbnd(:,daycol)==datesets(i), :);
      hfdayo = hfout(hfout(:,daycol)==datesets(i), :);
          
      %read horizontal optimal and orthogonal components
      JDAY = num2zeropadstr(jday,3);
      MO=day2month(jday,year);     % EXTERNAL function, day2month, get the month of one particular date
      direc=[datapath, '/arch', yr,'/',MO,'/'];     % directory name
      prename=[direc,yr,'.',JDAY,'.00.00.00.0000.CN'];    %  path plus prefix of data file,
%       disp(prename);
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

  %     keyboard
%       for j = 3: size(rangetemp,1)  
%         k = k+1;  
        k = idxbst(iii);
        disp(k);

        tmaxi = hfdayi(:, seccol); % starting time of max power rate of half sec inside the ellipse
        tmaxo = hfdayo(:, seccol); % starting time of max power rate of half sec outside the ellipse  
        tcnti = hfdayi(:, 15);  % the center of detecting win is the 15th col
        tcnto = hfdayo(:, 15);  % the center of detecting win is the 15th col

        tst = rangetemp(j,2); % start and end time of bursts
        ted = rangetemp(j,3);

        %how many 4-s detections fall into the burst range 
        indtmaxi = find(tmaxi>=tst-0.1 & tmaxi<=ted+0.1);
        ninbst(k,1) = length(indtmaxi); %+-0.1 s in case of resolution while saving to file
        
        %%%%Use the start and end of the 4-s detecting window
        tstbuf = min(tcnti(indtmaxi)-2);
        tedbuf = max(tcnti(indtmaxi)+2); 
        tlenbuf = tedbuf-tstbuf;
%%              
        %max allowable shift in best alignment
        msftaddmhf = 1.5*sps+1;  %+1 for safety
        msftaddm = 2*sps+1; 
        msftaddmdiff = msftaddm-msftaddmhf;
        
        %have some overshoot, so that the resulted rcc would have the same length as the signal
%         overshoot = rccmwlen/2;
        overshoot = rccmwlenhf/2;
        
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
        
        %%%according the linear fitting result of off12 and off13 VS. origin time, we sort of know how
        %%%much it needs to change the overall offset to change by 1 sample
        %generate overlapping windows of the same length
        indst = floor(tstbuf*sps+1);
        inded = floor(tedbuf*sps);
        if tlenbuf<25
          subwsec(k) = tlenbuf;
        else
          subwsec(k) = 25;   % determined from fitting
        end
        subwlen = subwsec(k)*sps;
        %since the rcc would lose rccmwlen/2 at both ends, this results in overlapping of 'rccmwlen' in rcc 
        %across consecutive windows; if use 'ovlplen' of rccmwlen, then rcc has no overlapping at all
  %       ovlplen = rccmwlen*2;
%         ovlplen = rccmwlen;
        ovlplen = rccmwlenhf;
        windows = movingwins(indst,inded,subwlen,ovlplen,0);
        nwin =  size(windows,1);

        off1iw = off1iwk{k};  % the best alignment between sta2, sta3 wrt sta1 for each subwin

        ircccat = [];   % concatenated indices of RCC
        irccran = zeros(nwin,2);  % start and end indices (range) of RCC of all subwins
        rcccat = [];  % average concatenated RCC
        rccpaircat = [];  % concatenated RCC between each station pair, order is 12, 13, 23
        mrccwpair = []; % median of concatenated RCC between each station pair
        ccwpair = []; % 0-lag overall cc of each subwin, between each station pair, order is 12, 13, 23 
        ramppaircat = []; % concatenated running amp at sta 1, 2, 3
        renvpaircat = []; % concatenated running env at sta 1, 2, 3
        indzc1cat = [];
        rcczc1cat = []; %ONLY at the zero-crossing with pos envelo
        renvzc1cat = [];
        renvzc1catnorm = []; 
        ramppairortcat = []; % concat ramp, ort comp
        renvpairortcat = []; % concat renv, ort comp
        ramppairvertcat = []; % concat ramp, vert comp
        renvpairvertcat = []; % concat renv, vert comp
        
        rccwlendiff = rccmwlen-rccmwlenhf;
        for iwin = 1: nwin
          isubwst = windows(iwin,1)-max(floor(tstbuf*sps+1-overshoot-msftaddm),1);
          isubwed = windows(iwin,2)-max(floor(tstbuf*sps+1-overshoot-msftaddm),1);
          isubwst = isubwst -rccwlendiff/2; %+msftaddmdiff
          isubwed = isubwed +rccwlendiff/2; %+msftaddmdiff
%           isubwst = windows(iwin,1)-max(floor(tstbuf*sps+1-overshoot-msftaddm),1);
%           isubwed = windows(iwin,2)-max(floor(tstbuf*sps+1-overshoot-msftaddm),1);
          if iwin == 1
            isubwst = isubwst-overshoot;
          end
          if iwin == nwin
            isubwed = isubwed+overshoot;
          end
          
          % %align records
          % optcc = detrend(optseg(isubwst: isubwed, 2:end));
          % msftadd = (round(max(abs([off12ran off13ran])))+1)*sps/40;  %+1 for safety
          % loffmax = 4*sps/40;
          % ccmid = ceil(size(optcc,1)/2);
          % ccwlen = round(size(optcc,1)-2*(msftadd+1));  % minus ensures successful shifting of records
          % ccmin = 0.01;  % depending on the length of trace, cc could be very low
          % iup = 1;    % times of upsampling        
          % [off12con,off13con,ccaliw(iwin,1)] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
          %   ccwlen,msftadd,loffmax,ccmin,iup);        
          % % if a better alignment cannot be achieved, use 0,0
          % if off12con == msftadd+1 && off13con == msftadd+1
          %   off12con = 0;
          %   off13con = 0;
          %   cc12 = xcorr(optcc(:,1), optcc(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
          %   cc13 = xcorr(optcc(:,1), optcc(:,3),0,'normalized');
          %   cc23 = xcorr(optcc(:,2), optcc(:,3),0,'normalized');
          %   ccaliw(iwin,1) = (cc12+cc13+cc23)/3;
          %   fprintf('Tremor burst %d cannot be properly aligned, double-check needed \n',k);
          % end
          % off1iw(iwin,1) = 0;
          % off1iw(iwin,2) = round(off12con);
          % off1iw(iwin,3) = round(off13con);
          
          %Align records
          optdat = [];  % win +/-3 s, segment of interest,first 1s will be tapered 
          ortdat = [];
          vertdat = [];
          optdat(:, 1) = optseg(isubwst: isubwed, 1); % time column
          ortdat(:, 1) = ortseg(isubwst: isubwed, 1);
          vertdat(:, 1) = vertseg(isubwst: isubwed, 1);
          for ista = 1: nsta
            optdat(:, ista+1) = optseg(isubwst-off1iw(iwin,ista): isubwed-off1iw(iwin,ista), ista+1); 
            ortdat(:, ista+1) = ortseg(isubwst-off1iw(iwin,ista): isubwed-off1iw(iwin,ista), ista+1);
            vertdat(:, ista+1) = vertseg(isubwst-off1iw(iwin,ista): isubwed-off1iw(iwin,ista), ista+1);
          end
          
          subw = zeros(size(optdat,1), nsta);
          for ista = 1: nsta
            tmp = optdat(:,ista+1); %best aligned, filtered
            %detrend and taper only the data, NOT the noise
            tmp = detrend(tmp);
            subw(:,ista) = tmp;
          end
          %compute running CC between 3 stations
          [irccw,rccw12,rccw13,rccw23] = RunningCC3sta(subw,rccmwlen);
          rccw = (rccw12+rccw13+rccw23)/3;
          irccw = irccw + windows(iwin,1) - windows(1,1) - rccwlendiff/2;   %convert to global index
          if iwin == 1
            irccw = irccw - overshoot;
          end
  %         plot(irccw, rccw);

          ircccat = [ircccat; irccw];
          irccran(iwin,:) = [irccw(1) irccw(end)];
          rcccat = [rcccat; rccw];
          rccpaircat = [rccpaircat; rccw12 rccw13 rccw23]; 
          mrccwpair = [mrccwpair; median(rccw12) median(rccw13) median(rccw23)]; 
          
          ccw12 = xcorr(subw(:,1), subw(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
          ccw13 = xcorr(subw(:,1), subw(:,3),0,'normalized');
          ccw23 = xcorr(subw(:,2), subw(:,3),0,'normalized');
          ccwpair = [ccwpair; ccw12 ccw13 ccw23];

          %compute the median absolute amplitude and envelope of the same moving window
          %for the moving window at the same station, sensable to use median
%           renvmwlen = rccmwlen;
          renvmwlen = rccmwlen/2;
          renvrccmwlendif = rccmwlen-renvmwlen;
          [ir,ramp1,renv1] = Runningampenv(subw(:,1),renvmwlen,renvmwlen-1,'median');
          [~,ramp2,renv2] = Runningampenv(subw(:,2),renvmwlen,renvmwlen-1,'median');
          [~,ramp3,renv3] = Runningampenv(subw(:,3),renvmwlen,renvmwlen-1,'median');
          %looks like using the amplitude and envelope are pretty similar
          %maybe better to mean, otherwise stations PGC will be underrepresented due to a larger amp
          %variation
          ramptemp = [ramp1 ramp2 ramp3];
          ramptemp = ramptemp(1+renvrccmwlendif/2: end-renvrccmwlendif/2, :);
          ramppaircat = [ramppaircat; ramptemp]; 
          renvtemp = [renv1 renv2 renv3];
          renvtemp = renvtemp(1+renvrccmwlendif/2: end-renvrccmwlendif/2, :);
          renvpaircat = [renvpaircat; renvtemp];

          %at orthogonal component
          subwort = zeros(size(ortdat,1), nsta);
          for ista = 1: nsta
            tmp = ortdat(:,ista+1); %best aligned, filtered
            tmp = detrend(tmp);
            subwort(:,ista) = tmp;
          end
          [ir,ramp1ort,renv1ort] = Runningampenv(subwort(:,1),renvmwlen,renvmwlen-1,'median');
          [~,ramp2ort,renv2ort] = Runningampenv(subwort(:,2),renvmwlen,renvmwlen-1,'median');
          [~,ramp3ort,renv3ort] = Runningampenv(subwort(:,3),renvmwlen,renvmwlen-1,'median');
          %looks like using the amplitude and envelope are pretty similar
          %maybe better to mean, otherwise stations PGC will be underrepresented due to a larger amp
          %variation
          ramptemport = [ramp1ort ramp2ort ramp3ort];
          ramptemport = ramptemport(1+renvrccmwlendif/2: end-renvrccmwlendif/2, :);
          ramppairortcat = [ramppairortcat; ramptemport]; 
          renvtemport = [renv1ort renv2ort renv3ort];
          renvtemport = renvtemport(1+renvrccmwlendif/2: end-renvrccmwlendif/2, :);
          renvpairortcat = [renvpairortcat; renvtemport];
          
          %at vertical component
          subwvert = zeros(size(vertdat,1), nsta);
          for ista = 1: nsta
            tmp = vertdat(:,ista+1); %best aligned, filtered
            tmp = detrend(tmp);
            subwvert(:,ista) = tmp;
          end
          [ir,ramp1vert,renv1vert] = Runningampenv(subwvert(:,1),renvmwlen,renvmwlen-1,'median');
          [~,ramp2vert,renv2vert] = Runningampenv(subwvert(:,2),renvmwlen,renvmwlen-1,'median');
          [~,ramp3vert,renv3vert] = Runningampenv(subwvert(:,3),renvmwlen,renvmwlen-1,'median');
          %looks like using the amplitude and envelope are pretty similar
          %maybe better to mean, otherwise stations PGC will be underrepresented due to a larger amp
          %variation
          ramptempvert = [ramp1vert ramp2vert ramp3vert];
          ramptempvert = ramptempvert(1+renvrccmwlendif/2: end-renvrccmwlendif/2, :);
          ramppairvertcat = [ramppairvertcat; ramptempvert]; 
          renvtempvert = [renv1vert renv2vert renv3vert];
          renvtempvert = renvtempvert(1+renvrccmwlendif/2: end-renvrccmwlendif/2, :);
          renvpairvertcat = [renvpairvertcat; renvtempvert];
                    
          %
%           subwtemp = subw(1+overshoot: end-overshoot, 1); %waveform at sta 1
%           %find all positive peaks
%           [~, ppkind] = findpeaks(subwtemp);
%           ppk = [ppkind ones(length(ppkind),1)];
%           [~, npkind] = findpeaks(-subwtemp);
%           npk = [npkind -ones(length(npkind),1)];
%           pk = [ppk; npk];
%           pk = sortrows(pk,1,'ascend');
%           ind=find(diff(pk(:,2))==2); %means it is from neg to positive
%           rcczc1 = [];
%           renvzc1 = [];
%           for ipk = 1: length(ind)
%             seg = subwtemp(pk(ind(ipk)): pk(ind(ipk)+1));
%             [~,ind] = min(abs(seg));  %zero-crossing
%             ind = ind-1+pk(ind(ipk));
%             rcczc1(ipk) = (rccw12(ind)+rccw13(ind))/2;
%             renvzc1(ipk) = (renv1(ind)+renv2(ind)+renv3(ind))/3;
%           end
%           rcczc1cat = [rcczc1cat; rcczc1];
%           renvzc1cat = [renvzc1cat; renvzc1];
          
        end
        irccrank{k} = irccran;
        ccwpairk{k} = ccwpair;
        mrccwpairk{k} = mrccwpair;
        
        %if only use the mean RCC from pair 12 and 13
        rcccat = mean(rccpaircat(:,[1 2]), 2);

        rcccomb = [];
        rcccomb(:,1) = rcccat;

        %mean running envelope of all 3 stations, normlized by the median of each burst, == ratio
        renvcat = mean(renvpaircat, 2);
        renvcatk{k} = renvcat;
        renvcatnorm = renvcat./median(renvcat); %distribute around 1
        renvcatnormk{k} = renvcatnorm;
        %mean running amp of all 3 stations
        rampcat = mean(ramppaircat, 2);  
        rampcatk{k} = rampcat;
        rampcatnorm = rampcat./median(rampcat);
        rampcatnormk{k} = rampcatnorm; %distribute around 1    

        %mean running envelope of all 3 stations, ort. comp.
        renvortcat = mean(renvpairortcat, 2);
        renvortcatk{k} = renvortcat;
        renvortcatnorm = renvortcat./median(renvortcat); %distribute around 1
        renvortcatnormk{k} = renvortcatnorm;
        %mean running amp of all 3 stations, normlized by the median of each burst, == ratio, ort. comp.
        ramportcat = mean(ramppairortcat, 2);  
        ramportcatk{k} = ramportcat;
        ramportcatnorm = ramportcat./median(ramportcat);
        ramportcatnormk{k} = ramportcatnorm; %distribute around 1    

        %mean running envelope of all 3 stations, vert. comp.
        renvvertcat = mean(renvpairvertcat, 2);
        renvvertcatk{k} = renvvertcat;
        renvvertcatnorm = renvvertcat./median(renvvertcat); %distribute around 1
        renvvertcatnormk{k} = renvvertcatnorm;
        %mean running amp of all 3 stations, normlized by the median of each burst, == ratio, ort. comp.
        rampvertcat = mean(ramppairvertcat, 2);  
        rampvertcatk{k} = rampvertcat;
        rampvertcatnorm = rampvertcat./median(rampvertcat);
        rampvertcatnormk{k} = rampvertcatnorm; %distribute around 1    

        
        %%%Align and compute the RCC based on the entire win, and take that as the input signal!      
        optdat = [];  % win segment of interest
        ortdat = [];
        vertdat = [];
        optdat(:, 1) = optseg(1+msftaddm-rccwlendiff/2: end-msftaddm+rccwlendiff/2, 1); % time column
        ortdat(:, 1) = ortseg(1+msftaddm-rccwlendiff/2: end-msftaddm+rccwlendiff/2, 1);
        vertdat(:, 1) = vertseg(1+msftaddm-rccwlendiff/2: end-msftaddm+rccwlendiff/2, 1);
        for ista = 1: nsta 
          optdat(:, ista+1) = optseg(1+msftaddm-rccwlendiff/2-off1i(k,ista): end-msftaddm+rccwlendiff/2-off1i(k,ista), ista+1);
          ortdat(:, ista+1) = ortseg(1+msftaddm-rccwlendiff/2-off1i(k,ista): end-msftaddm+rccwlendiff/2-off1i(k,ista), ista+1);
          vertdat(:, ista+1) = vertseg(1+msftaddm-rccwlendiff/2-off1i(k,ista): end-msftaddm+rccwlendiff/2-off1i(k,ista), ista+1);
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
        [ircc,rcc12,rcc13,rcc23] = RunningCC3sta(sigsta,rccmwlen);
        ircc = ircc-overshoot - rccwlendiff/2;
        rcc = (rcc12+rcc13+rcc23)/3;
%         rcc1i = zeros(length(rcc),nsta-3);
%         for ista = 4:nsta
%           [~,rcc1i(:,ista-3)] = RunningCC(sigsta(:,1), sigsta(:,ista), rccmwlen);
%         end
        %if only use pairs 12 and 13
        rcc = (rcc12+rcc13)/2;
        rcccomb(:,2) = rcc;
        
        rccbst{k} = rcccomb;
        
        sigsta = detrend(sigsta(1+overshoot+rccwlendiff/2: ...
          end-overshoot-rccwlendiff/2, :));  %excluding the overshoot

        %%%find all positive peaks at sta 1
        wintemp = sigsta(:,1);
        [~, ppkind] = findpeaks(wintemp);
        ppk = [ppkind ones(length(ppkind),1)];
        [~, npkind] = findpeaks(-wintemp);
        npk = [npkind -ones(length(npkind),1)];
        pk = [ppk; npk];
        pk = sortrows(pk,1,'ascend');
        ind=find(diff(pk(:,2))==2); %means it is from neg to positive
        for ipk = 1: length(ind)
          seg = wintemp(pk(ind(ipk)): pk(ind(ipk)+1));
          [~,indzc1] = min(abs(seg));  %zero-crossing
          indzc1 = indzc1-1+pk(ind(ipk));
          indzc1cat(ipk,1) = indzc1;
          rcczc1cat(ipk,1) = rcccat(indzc1);
          renvzc1cat(ipk,1) = renvcat(indzc1);
          renvzc1catnorm(ipk,1) = renvcatnorm(indzc1);
        end
        rcczc1catk{k} = rcczc1cat;
        renvzc1catk{k} = renvzc1cat;
        renvzc1catnormk{k} = renvzc1catnorm;
        indzc1catk{k} = indzc1cat;
        
        %%%for ort. comp
        sigstaort = zeros(size(ortdat,1), nsta);
        for ista = 1:nsta
          tmp = ortdat(:,ista+1); %best aligned, filtered
          tmp = detrend(tmp);
          sigstaort(:,ista) = tmp;
        end
        
        %%%for vert. comp
        sigstavert = zeros(size(vertdat,1), nsta);
        for ista = 1:nsta
          tmp = vertdat(:,ista+1); %best aligned, filtered
          tmp = detrend(tmp);
          sigstavert(:,ista) = tmp;
        end

        %%%%%%%%%% for 3-sta srcs
        %%%what are the corresponding RCC at each source
        ist = sum(nsrc(1:k-1))+1;
        ied = ist+nsrc(k)-1;
        impindepst = imp(ist:ied,:);

        rccpairsrc = [];
        rccpairsrc(:,1) = rccpaircat(round(mean(impindepst(:,[1 3]),2)),1);
        rccpairsrc(:,2) = rccpaircat(round(mean(impindepst(:,[1 5]),2)),2);
        rccpairsrc(:,3) = rccpaircat(round(mean(impindepst(:,[3 5]),2)),3);
        rccpairsrcall = [rccpairsrcall; rccpairsrc];
        
        %use the concatenated rcc at the average arrival time of each source
        rcccatsrc = [];
        rcccatsrc(:,1) = rcccat(round(mean(impindepst(:,[1 3 5]),2)));
        rcccatsrcall = [rcccatsrcall; rcccatsrc];
        
        %%%similarly, what are the corresponding renv and ramp at each source 
        %running envelope at src arrival
        renvpairsrc = [];
        renvpairsrc(:,1) = renvpaircat(impindepst(:,1), 1);
        renvpairsrc(:,2) = renvpaircat(impindepst(:,3), 2);
        renvpairsrc(:,3) = renvpaircat(impindepst(:,5), 3);
        renvpairsrcall = [renvpairsrcall; renvpairsrc];
        renvcatsrc = [];
        renvcatsrc(:,1) = renvcat(round(mean(impindepst(:,[1 3 5]),2)));
        renvcatsrcall = [renvcatsrcall; renvcatsrc];
        renvcatnormsrc = [];
        renvcatnormsrc(:,1) = renvcatnorm(round(mean(impindepst(:,[1 3 5]),2)));
        renvcatnormsrcall = [renvcatnormsrcall; renvcatnormsrc];

        %running amplitude at src arrival
        ramppairsrc = [];
        ramppairsrc(:,1) = ramppaircat(impindepst(:,1), 1);
        ramppairsrc(:,2) = ramppaircat(impindepst(:,3), 2);
        ramppairsrc(:,3) = ramppaircat(impindepst(:,5), 3);
        ramppairsrcall = [ramppairsrcall; ramppairsrc];
        rampcatsrc = [];
        rampcatsrc(:,1) = rampcat(round(mean(impindepst(:,[1 3 5]),2)));
        rampcatsrcall = [rampcatsrcall; rampcatsrc];
        rampcatnormsrc = [];
        rampcatnormsrc(:,1) = rampcatnorm(round(mean(impindepst(:,[1 3 5]),2)));
        rampcatnormsrcall = [rampcatnormsrcall; rampcatnormsrc];

        %running envelope at src arrival, ort. comp.
        renvpairortsrc = [];
        renvpairortsrc(:,1) = renvpairortcat(impindepst(:,1), 1);
        renvpairortsrc(:,2) = renvpairortcat(impindepst(:,3), 2);
        renvpairortsrc(:,3) = renvpairortcat(impindepst(:,5), 3);
        renvpairortsrcall = [renvpairortsrcall; renvpairortsrc];
        renvortcatsrc = [];
        renvortcatsrc(:,1) = renvortcat(round(mean(impindepst(:,[1 3 5]),2)));
        renvortcatsrcall = [renvortcatsrcall; renvortcatsrc];
        renvortcatnormsrc = [];
        renvortcatnormsrc(:,1) = renvortcatnorm(round(mean(impindepst(:,[1 3 5]),2)));
        renvortcatnormsrcall = [renvortcatnormsrcall; renvortcatnormsrc];

        %running envelope at src arrival, vert. comp.
        renvpairvertsrc = [];
        renvpairvertsrc(:,1) = renvpairortcat(impindepst(:,1), 1);
        renvpairvertsrc(:,2) = renvpairortcat(impindepst(:,3), 2);
        renvpairvertsrc(:,3) = renvpairortcat(impindepst(:,5), 3);
        renvpairvertsrcall = [renvpairvertsrcall; renvpairvertsrc];
        renvvertcatsrc = [];
        renvvertcatsrc(:,1) = renvortcat(round(mean(impindepst(:,[1 3 5]),2)));
        renvvertcatsrcall = [renvvertcatsrcall; renvvertcatsrc];
        renvvertcatnormsrc = [];
        renvvertcatnormsrc(:,1) = renvvertcatnorm(round(mean(impindepst(:,[1 3 5]),2)));
        renvvertcatnormsrcall = [renvvertcatnormsrcall; renvvertcatnormsrc];
        
        %in particular at the src timing of sta 1
        rcccatsrc1 = [];
        rcccatsrc1(:,1) = rcccat(impindepst(:,1));
        rcccatsrc1all = [rcccatsrc1all; rcccatsrc1];
        renvcatsrc1 = [];
        renvcatsrc1(:,1) = renvcat(impindepst(:,1));
        renvcatsrc1all = [renvcatsrc1all; renvcatsrc1];
        renvcatnormsrc1 = [];
        renvcatnormsrc1(:,1) = renvcatnorm(impindepst(:,1));
        renvcatnormsrc1all = [renvcatnormsrc1all; renvcatnormsrc1];
        
        %at the zero-crossing at sta 1 that is closest to src arrival
        rcccatsrczc1 = [];
        renvcatsrczc1 = [];
        renvcatnormsrczc1 = [];
        diffzcsrc1 = [];
        for isrc = 1: nsrc(k)
          [diffzcsrc1(isrc,1),ind] = min(abs(impindepst(isrc,1)-indzc1cat));
          rcccatsrczc1(isrc,1) = rcccat(indzc1cat(ind));
          renvcatsrczc1(isrc,1) = renvcat(indzc1cat(ind));
          renvcatnormsrczc1(isrc,1) = renvcatnorm(indzc1cat(ind));
        end
        diffzcsrc1all = [diffzcsrc1all; diffzcsrc1];
        rcccatsrczc1all = [rcccatsrczc1all; rcccatsrczc1];
        renvcatsrczc1all = [renvcatsrczc1all; renvcatsrczc1];
        renvcatnormsrczc1all = [renvcatnormsrczc1all; renvcatnormsrczc1];

        %%%%%%%%%% for 4-sta srcs
        %%%what are the corresponding RCC at each source
        ist = sum(nsrc4th(1:k-1))+1;
        ied = ist+nsrc4th(k)-1;
        impindepst = imp4th(ist:ied,:);
        rccpairsrc4th = [];
        rccpairsrc4th(:,1) = rccpaircat(round(mean(impindepst(:,[1 3]),2)),1);
        rccpairsrc4th(:,2) = rccpaircat(round(mean(impindepst(:,[1 5]),2)),2);
        rccpairsrc4th(:,3) = rccpaircat(round(mean(impindepst(:,[3 5]),2)),3);
        rccpairsrc4thall = [rccpairsrc4thall; rccpairsrc4th];
        
        %use the concatenated rcc at the average arrival time of each source
        rcccatsrc4th = [];
        rcccatsrc4th(:,1) = rcccat(round(mean(impindepst(:,[1 3 5]),2)));
        rcccatsrc4thall = [rcccatsrc4thall; rcccatsrc4th];

        %%%similarly, what are the corresponding renv and ramp at each source 
        %running envelope at src arrival
        renvpairsrc4th = [];
        renvpairsrc4th(:,1) = renvpaircat(impindepst(:,1), 1);
        renvpairsrc4th(:,2) = renvpaircat(impindepst(:,3), 2);
        renvpairsrc4th(:,3) = renvpaircat(impindepst(:,5), 3);
        renvpairsrc4thall = [renvpairsrc4thall; renvpairsrc4th];
        renvcatsrc4th = [];
        renvcatsrc4th(:,1) = renvcat(round(mean(impindepst(:,[1 3 5]),2)));
        renvcatsrc4thall = [renvcatsrc4thall; renvcatsrc4th];
        renvcatnormsrc4th = [];
        renvcatnormsrc4th(:,1) = renvcatnorm(round(mean(impindepst(:,[1 3 5]),2)));
        renvcatnormsrc4thall = [renvcatnormsrc4thall; renvcatnormsrc4th];

        %running amplitude at src arrival
        ramppairsrc4th = [];
        ramppairsrc4th(:,1) = ramppaircat(impindepst(:,1), 1);
        ramppairsrc4th(:,2) = ramppaircat(impindepst(:,3), 2);
        ramppairsrc4th(:,3) = ramppaircat(impindepst(:,5), 3);
        ramppairsrc4thall = [ramppairsrc4thall; ramppairsrc4th];
        rampcatsrc4th = [];
        rampcatsrc4th(:,1) = rampcat(round(mean(impindepst(:,[1 3 5]),2)));
        rampcatsrc4thall = [rampcatsrc4thall; rampcatsrc4th];
        rampcatnormsrc4th = [];
        rampcatnormsrc4th(:,1) = rampcatnorm(round(mean(impindepst(:,[1 3 5]),2)));
        rampcatnormsrc4thall = [rampcatnormsrc4thall; rampcatnormsrc4th];

        %in particular at the src timing of sta 1
        rcccatsrc14th = [];
        rcccatsrc14th(:,1) = rcccat(impindepst(:,1));
        rcccatsrc14thall = [rcccatsrc14thall; rcccatsrc14th];
        renvcatsrc14th = [];
        renvcatsrc14th(:,1) = renvcat(impindepst(:,1));
        renvcatsrc14thall = [renvcatsrc14thall; renvcatsrc14th];
        renvcatnormsrc14th = [];
        renvcatnormsrc14th(:,1) = renvcatnorm(impindepst(:,1));
        renvcatnormsrc14thall = [renvcatnormsrc14thall; renvcatnormsrc14th];

%       end    
%     end   
  end
%%
  %concatenated from all 25-s wins, separated by bursts 
  allbstsig.rccbst = rccbst;
  allbstsig.rampcatk = rampcatk;
  allbstsig.rampcatnormk = rampcatnormk;
  allbstsig.renvcatk = renvcatk;
  allbstsig.renvcatnormk = renvcatnormk;
  %only at zero-crossings with positive slopes
  allbstsig.rcczc1catk = rcczc1catk;
  allbstsig.renvzc1catk = renvzc1catk;
  allbstsig.renvzc1catnormk = renvzc1catnormk;
  allbstsig.indzc1catk = indzc1catk;

  %at 3-sta src arrival time averaged at 3 stas 
  allbstsig.rcccatsrcall = rcccatsrcall;
  allbstsig.rccpairsrcall = rccpairsrcall;
  allbstsig.renvcatsrcall = renvcatsrcall;
  allbstsig.renvcatnormsrcall = renvcatnormsrcall;
  allbstsig.renvpairsrcall = renvpairsrcall;
  allbstsig.rampcatsrcall = rampcatsrcall;
  allbstsig.rampcatnormsrcall = rampcatnormsrcall;
  allbstsig.ramppairsrcall = ramppairsrcall;
  %at 3-sta src arrival time only at sta 1
  allbstsig.rcccatsrc1all = rcccatsrc1all;
  allbstsig.renvcatsrc1all = renvcatsrc1all;
  allbstsig.renvcatnormsrc1all = renvcatnormsrc1all;
  %at the zero-crossing at sta 1 that is closest to src arrival
  allbstsig.diffzcsrc1all = diffzcsrc1all;
  allbstsig.rcccatsrczc1all = rcccatsrczc1all;
  allbstsig.renvcatsrczc1all = renvcatsrczc1all;
  allbstsig.renvcatnormsrczc1all = renvcatnormsrczc1all;
  
  %at ort. comp.
  allbstsig.renvortcatk = renvortcatk;
  allbstsig.renvortcatnormk = renvortcatnormk;
  allbstsig.ramportcatk = ramportcatk;
  allbstsig.ramportcatnormk = ramportcatnormk; %distribute around 1
  allbstsig.renvortcatsrcall = renvortcatsrcall;
  allbstsig.renvortcatnormsrcall = renvortcatnormsrcall;
  allbstsig.renvpairortsrcall = renvpairortsrcall;
  
  %at vert. comp.
  allbstsig.renvvertcatk = renvvertcatk;
  allbstsig.renvvertcatnormk = renvvertcatnormk;
  allbstsig.rampvertcatk = rampvertcatk;
  allbstsig.rampvertcatnormk = rampvertcatnormk; %distribute around 1
  allbstsig.renvvertcatsrcall = renvvertcatsrcall;
  allbstsig.renvvertcatnormsrcall = renvvertcatnormsrcall;
  allbstsig.renvpairvertsrcall = renvpairvertsrcall;

  %at 4-sta src arrival time averaged at 3 stas 
  allbstsig.rcccatsrc4thall = rcccatsrc4thall;
  allbstsig.rccpairsrc4thall = rccpairsrc4thall;
  allbstsig.renvcatsrc4thall = renvcatsrc4thall;
  allbstsig.renvcatnormsrc4thall = renvcatnormsrc4thall;
  allbstsig.renvpairsrc4thall = renvpairsrc4thall;
  allbstsig.rampcatsrc4thall = rampcatsrc4thall;
  allbstsig.rampcatnormsrc4thall = rampcatnormsrc4thall;
  allbstsig.ramppairsrc4thall = ramppairsrc4thall;
  %at 4-sta src arrival time only at sta 1
  allbstsig.rcccatsrc14thall = rcccatsrc14thall;
  allbstsig.renvcatsrc14thall = renvcatsrc14thall;
  allbstsig.renvcatnormsrc14thall = renvcatnormsrc14thall;


  %save all workspace variables
  save(strcat('seisbursts_allvari',num2str(losig),'-',num2str(hisig),'.mat'),'allbstsig');

end












