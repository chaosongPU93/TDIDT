% detection002_4s.m
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this script, I want to grasp all we have learned and paid attention to
% till now to re-do the 4-s HF detection around 002 region. Since there is no
% comparison between HF and LF, this script does not have to distinguish the
% frequency band here, but we are using 1.25-6.5 Hz.
%
% Should be essentially similar to 'hf_pgc.m'; part related to R_o/o comes
% from 'hf.m'
%
% --PGC-SSIB-SILB trio in Cascadia, using data at 40 sps, then add extra
%   stations to check the validity of detections
% --There will be no filtering effect correction, because we only focus on
%   HF here, so there is no relative change in locations.
% --Upsampling to 1/4 sample at 40 sps.
% --Using the rots parameters from the new LFE catalog.
% --Using the quality control, R_o/o, (not used in any of the previous detection)
% --We will try to interpolate the spatial locations from an existing grid of
%   locations inverted from a grid of time offsets +/-13 samples at 20 sps
%   which itself is relatively regular. Previously, we directly invert for the
%   location of every source, which may be too demanding in the spatial
%   resolution, causing some unpredictable distortions. 
% --To be consistent with above, the mshift has been changed to 26 from 29 
%   samples at 40 sps, which is +/-13 at 20 sps.
% --I shall reiterate the order of the entire flow of detection here: 
%   'detection002_4s.m' --> merge all days and clean up --> DO NOT
%   do inversion, just interpolate --> plot density map, choose the 002 region
%   with high density, 1.75*1 km ellipse or ~a cicle --> find detections inside,
%   so to get the windows of bursts for future deconvolution
% 
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2022/03/10
% Last modified date:   2022/03/10
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format short e   % Set the format to 5-digit floating point
clear
clc
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
temppath = strcat(workpath, '/templates/');
%%% choose the data path here
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/PGCtrio');

FLAG = 'PGC';

fampool = [
  '002';
  ];    % family number

for ifam = 1: size(fampool,1)
  fam = fampool(ifam, :)
  
  freqflag='hf';  % flag to indicate whether to do hf or lf;
  
  % get days, 3-sta trio name, bostock's template name, zero crossing index
  % of templates
  % timoffrot= [
  %             2003 061;
  %             2003 062;
  %             2003 063;
  %             2004 196;
  %             2004 197;
  %             2004 198;
  %             2004 199;
  %             2005 254;
  %             2005 255;
  %             2005 256];
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
  
  
  %Sampling rate of the data will be used
  sps=40;     % samples per second
  
  %%% IMPORTANT, NEED change sometimes
  %Basics of the cross-correlation:  Window length, number of windows, filter parameters, etc.
  winlensec=4;     % offsec = 3 was used in first-year report
  winoffsec=1;        % window offset in sec, which is the step of a moving window
  winlen=winlensec*sps;      % length in smaples
  winoff=winoffsec*sps;      % offset in samples
  tracelen=86400*sps; %one day of data at 40 sps, overall trace length, 24*3600
  cutsec = 2;
  winbig=2*(tracelen/2-(cutsec*sps)); %ignore 2 seconds at each end of day, bigger window contains positive and negative, why tracelen/2?, -4s
  timbig=winbig/(2*sps); %half that time, in seconds, half day - 2s
  igstart=floor(tracelen/2-winbig/2)+1; %start counting seis data from here, 2*sps+1, floor(1.8)==1
  nwin=floor((winbig-winlen)/winoff)+1;    % number of windows, first one not included, ADD +1 by Chao, 2019/02/17
  
  hi=6.5;    % frequency band
  lo=1.25;
  
  npo=2;     % poles, passes of filters
  npa=2;
  if isequal(fam,'002')
    cyclskip = 0;
    mshift=26+cyclskip; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
    loopoffmax=2.1; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
    xcmaxAVEnmin=0.44; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?
    %%% xcmaxAVEnmin = x-correlation max average min ?
  end
  
  %%% 'scalefact' scales templates; 'scaleseisms' scales seisms.  Strategy changes with family.
  if isequal(fam,'002')
    whichtoplot=2;          % flag to which to plot
    scaleseisms=[1.0 0.76 0.95];       % scaleseisms scales seismograms
  elseif isequal(fam,'068')
    whichtoplot=1;
    scaleseisms=[1.0 0.6 0.6];
  else
    scaleseisms=[1.0 0.76 0.95];
  end
  
  %% START TO LOOP FOR every day
  %cycle over each day:
  for nd=1: nday      % num of rows, also num of days
    %     nd=1;
    close all
    
    year=timoffrot(nd,1);
    jday=timoffrot(nd,2);
    
    YEAR=int2str(year);
    JDAY = num2zeropadstr(jday,3);
    MO=day2month(jday,year);     % EXTERNAL function, day2month, get the month of one particular date
    IDENTIF=[YEAR,'.',JDAY,'.',fam,'.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),...
      '.',int2str(npo),int2str(npa),'.ms', int2str(mshift)]
    direc=[datapath, '/arch', YEAR,'/',MO,'/'];     % directory name
    prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];    %  path plus prefix of data file,
       
    %read horizontal optimal and orthogonal components
    [STAopt,STAort,STAnzeros,fileflag] = rd_daily_bpdata(year,jday,prename,stas,PERMSTA,POLSTA,...
      PERMROTS,POLROTS,sps,lo,hi,npo,npa,nwin,winlen,winoff,igstart);
    %rearrange the output
    timsSTA = STAopt(:,1)';
    STAopt = STAopt(:,2:end)';
    STAort = STAort(:,2:end)';
    STAnzeros = STAnzeros';
    %     %scale seismogram relatively, from historical codes, maybe not ncecessary
    %     for ista = 1: nsta
    %       STAopt(:,ista+1) = STAopt(:,ista+1)/scaleseisms(ista);
    %     end
    
    %%
    %%%%%2022/05/01, add only for presentation purpose
    [windows] = movingwins(1,20*sps,winlen,winlen-winoff,0);
    figure
    ax=gca; hold(ax,'on'); box on
    plot(1/sps:1/sps:20, STAopt(1, 5*sps+1:25*sps),'r');
    plot(1/sps/sps:1/sps:20, STAopt(2, 5*sps+1:25*sps),'b');
    plot(1/sps:1/sps:20, STAopt(3, 5*sps+1:25*sps),'k');
    text(0.05,0.95,'Velocity, optimal','Units','normalized','FontSize',12);
    text(0.95,0.95,sprintf('Bandpassed by %.1f-%.1f Hz',lo,hi),'Units','normalized','FontSize',12,...
      'HorizontalAlignment','right');
    xlabel('Time (s)');
    ylabel('Amplitude');
    %%%%%
    %%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Autocorrelation of stations.  Those that end in "sq" are the running
    %   cumulative sum, to be used later by differencing the window edpoints.
    %   (Used to be PGCauto, PGC2, SSIBauto, SSIB2, etc.)
    %   Station to itself is in a 3 x tracelen array
    %   Cross-station measurements are in their own linear array
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    STAauto=STAopt.*STAopt;
    %     STAsq=cumsum(STAauto,2);    % cumulative sum of each trace(row)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Cross-correlation between stations, with small offsets up to +/- mshift.
    %  First index is pointwise multiplication of traces; second is shifting offset.
    %  lenx is shorter than tracelen by mshift at each end (see notebook sketch)
    %  For STA12 and PGSI, SSI and SIL are shifted relative to PGC, by 1 each time through loop.
    %  For SISS, SSI is shifted relative to SILB.
    %  cumsumSTA12 etc. are the running cumulative sum of the x-correlation.
    %  PGSSx becomes STA12x, PGSI -> STA13, SISS -> STA32
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lenx=tracelen-2*mshift;     % 86400*sps-19*2
    STA12x=zeros(lenx, 2*mshift+1);    % 19*2+1
    STA13x=zeros(lenx, 2*mshift+1);    % stas: 1->PGC, 2->SSIB, 3->SILB
    STA32x=zeros(lenx, 2*mshift+1);
    %%% SEE NOTES, Chao
    for n=-mshift:mshift
      % PGC corr SSIB, 1+mshift:tracelen-mshift == 1+mshift-n:tracelen-mshift-n == lenx
      STA12x(:,n+mshift+1)=STAopt(1,1+mshift:tracelen-mshift).* ...
        STAopt(2,1+mshift-n:tracelen-mshift-n);
      % PGC corr SILB
      STA13x(:,n+mshift+1)=STAopt(1,1+mshift:tracelen-mshift).* ...
        STAopt(3,1+mshift-n:tracelen-mshift-n);
      % SILB corr SSIB
      STA32x(:,n+mshift+1)=STAopt(3,1+mshift:tracelen-mshift).* ...
        STAopt(2,1+mshift-n:tracelen-mshift-n);
    end
    %     cumsumSTA12=cumsum(STA12x);  %prev cumsumPGSS    % ==cumsum(STA12x,1),column
    %     cumsumSTA13=cumsum(STA13x);  %prev cumsumPGSI
    %     cumsumSTA32=cumsum(STA32x);  %prev cumsumSISS
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  "winbig"             is now the whole day, minus 2 sec at each end (apparently).
    %  "timbig"             is the time of half that.
    %  'igstart'            is the index of the starting sample.
    %  'winlen and winoff'  refer to the small windows.
    %  'timswin'            refers to the central times of those small windows.
    %  'sumsPGSS' (etc.)    is the cross-correlation sum over the window.  The first
    %    index refers to the window number and the second the shift over +/-mshift.
    %
    %  Normalized x-correlation:
    %    For PGSS and PGSI, for a given window PGC does not shift but SSI and
    %    SIL do.  So can compute sumsPGC2 (from the running cum. sum PGC2) just
    %    once for each window.  Same for sumsSILB2b.  But for the stations that
    %    shift, SSI and SIL (for PGC) and SSI (for SIL), must compute sumsSSIB2
    %    and sumsSILB2 upon each shift (actually, this is is easy book-keeping
    %    but not efficient).  Again, the first index refers to the window
    %    number and the second the shift over +/-mshift.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timswin=zeros(nwin,1);    %%% timswin refers to the central times of those small windows.
    sumsSTA12=zeros(nwin,2*mshift+1);   % PGC-SSIB
    sumsSTA13=zeros(nwin,2*mshift+1);   % PGC-SILB
    sumsSTA32=zeros(nwin,2*mshift+1);   % SILB-SSIB
    sumsSTA1sq=zeros(nwin,2*mshift+1);  % "sq" are the running cumulative sum, to be used later by differencing the window edpoints
    sumsSTA2sq=zeros(nwin,2*mshift+1);
    sumsSTA3sq=zeros(nwin,2*mshift+1);
    sumsSTA3Bsq=zeros(nwin,2*mshift+1); % refers the shifting SSI for SIL is not moving?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  sumsPGSS is shorter than sumsPGC2 by 2*mshift.  This is why sumsPGC2 etc
    %  is shifted by +mshift.  cumsumPGSS(1,:)=cumsum(PGSSx)(1,:) starts mshift
    %  to the right of the first data sample.  igstart is how many to the right
    %  of that.
    %  6/29/2018:  I'm pretty sure I want to subtract mshift from every iend or
    %  istart index to the right of the equal signs in the following FOR loop.
    %  The following looks like it (previously) had the proper relative shifts but the
    %  absolute time wasn't registered to the global igstart properly.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for n=1:nwin
      istart=igstart+(n-1)*winoff;     % 2*sps+1 + (n-1)*3*sps, igstart is the index of the starting sample; istart is index of each sample window
      iend=istart+winlen-1;              % + 12.5*sps, ADD -1 by Chao, 2019/02/17
      timswin(n)=timsSTA(istart+winlen/2);    % timsSTA, time serie data of STA; timswin, center of win, also == (istart+iend)/2
      sumsSTA12(n,:)=sum(STA12x(istart-mshift: iend-mshift, :));
      sumsSTA13(n,:)=sum(STA13x(istart-mshift: iend-mshift, :));
      sumsSTA32(n,:)=sum(STA32x(istart-mshift: iend-mshift, :));
      sumsSTA1sq(n,:)=sum(STAauto(1, istart: iend));
      sumsSTA3Bsq(n,:)=sum(STAauto(3, istart: iend));
      for m=-mshift:mshift
        sumsSTA2sq(n,m+mshift+1)=sum(STAauto(2, istart-m: iend-m)); %+m??? (yes).
        sumsSTA3sq(n,m+mshift+1)=sum(STAauto(3, istart-m: iend-m));
      end
      
      %         sumsSTA12(n,:)=cumsumSTA12(iend-mshift,:)-cumsumSTA12(istart-mshift-1,:);  %Yes, -mshift (6/29/18), sum of square of points between istart-mshift and iend-shift
      %         sumsSTA13(n,:)=cumsumSTA13(iend-mshift,:)-cumsumSTA13(istart-mshift-1,:);  %%% result in the sum from istart-mshift to iend-mshift
      %         sumsSTA32(n,:)=cumsumSTA32(iend-mshift,:)-cumsumSTA32(istart-mshift-1,:);
      %         sumsSTA1sq(n,:)=STAsq(1,iend)-STAsq(1,istart-1);  %PGC2 is cumsummed. Yes, +mshift.  No, no mshift (6/29/18)
      %         sumsSTA3Bsq(n,:)=STAsq(3,iend)-STAsq(3,istart-1); %Similar, for the SILB-SSIB connection, here SILB does not shift
      %         for m=-mshift:mshift
      %             sumsSTA2sq(n,m+mshift+1)=STAsq(2,iend-m)-STAsq(2,istart-1-m); %+m??? (yes).
      %             sumsSTA3sq(n,m+mshift+1)=STAsq(3,iend-m)-STAsq(3,istart-1-m);
      %         end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Denominator for the normalization.  A 2D array, nwin by 2*mshift+1.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %An attempt to bypass glitches in data.  Min value of good data typically ~10^{-2}
    glitches=1.e-7;
    sumsSTA1sq=max(sumsSTA1sq,glitches);
    sumsSTA2sq=max(sumsSTA2sq,glitches);
    sumsSTA3sq=max(sumsSTA3sq,glitches);    % return maximum between A and B
    %
    denomSTA12n=realsqrt(sumsSTA1sq.*sumsSTA2sq);    % Real square root, An error is produced if X is negative
    denomSTA13n=realsqrt(sumsSTA1sq.*sumsSTA3sq);
    denomSTA32n=realsqrt(sumsSTA3Bsq.*sumsSTA2sq);
    %
    sumsSTA12n=sumsSTA12./denomSTA12n;   % suffix 'n' means normalized
    sumsSTA13n=sumsSTA13./denomSTA13n;
    sumsSTA32n=sumsSTA32./denomSTA32n;
    [xcmaxSTA12n,imaxSTA12]=max(sumsSTA12n,[],2);   %Integer-offset max cross-correlation
    [xcmaxSTA13n,imaxSTA13]=max(sumsSTA13n,[],2);   % along row, max cc val and index in each window
    [xcmaxSTA32n,imaxSTA32]=max(sumsSTA32n,[],2);
    %Parabolic fit:
    [xmaxSTA12n,ymaxSTA12n,aSTA12]=parabol(nwin,mshift,sumsSTA12n,imaxSTA12); %Interpolated max cross-correlation
    [xmaxSTA13n,ymaxSTA13n,aSTA13]=parabol(nwin,mshift,sumsSTA13n,imaxSTA13);
    [xmaxSTA32n,ymaxSTA32n,aSTA32]=parabol(nwin,mshift,sumsSTA32n,imaxSTA32);
    
    %h=figure('Position',[0.1*wid 1 2.5*wid hite]); %center
    
    ix=sub2ind(size(denomSTA12n),(1:nwin)',imaxSTA12); %Find the linear index of the largest denominator
    ampSTA12=sqrt(denomSTA12n(ix)); %This makes amplitude linear rather than quadratic with counts.   % JUST FOR EASIER USAGE
    ampSTA1sq=sumsSTA1sq(ix); %by construction PGC2 is the same for all shifts  % sumsPGC2 becomes sumsSTA1sq
    ampSTA2sq=sumsSTA2sq(ix); % sumsSSIB2 becomes sumsSTA2sq, NOTICE: here ampSTA1sq are still sum of sqaures, i.e., quadratic
    ix=sub2ind(size(denomSTA13n),(1:nwin)',imaxSTA13);
    ampSTA13=sqrt(denomSTA13n(ix));
    ampSTA3sq=sumsSTA3sq(ix);
    ix=sub2ind(size(denomSTA32n),(1:nwin)',imaxSTA32);
    ampSTA32=sqrt(denomSTA32n(ix));
    AmpComp(1:4)=0;       % amplitude compare
    %AmpComp seems to be amplitude squared in 4s window minus amp squared in prior window,
    %divided by sum of amp squared in the two windows.  And why?
    AmpComp(5:nwin)=((ampSTA1sq(5:nwin)+ampSTA2sq(5:nwin)+ampSTA3sq(5:nwin))- ...
      (ampSTA1sq(1:nwin-4)+ampSTA2sq(1:nwin-4)+ampSTA3sq(1:nwin-4)))./ ...
      ((ampSTA1sq(5:nwin)+ampSTA2sq(5:nwin)+ampSTA3sq(5:nwin))+ ...
      (ampSTA1sq(1:nwin-4)+ampSTA2sq(1:nwin-4)+ampSTA3sq(1:nwin-4))) ;
    
    %%%%%%%%%%%%
    % add by Chao, to clear unnecessary variables to free memory in case of
    % crash
    clear denomSTA12n denomSTA13n denomSTA32n sumsSTA12 sumsSTA13 sumsSTA32
    clear sumsSTA1sq sumsSTA2sq sumsSTA3sq sumsSTA3Bsq
    %%%%%%%%%%%%
    
    %Center them
    imaxSTA12cent=imaxSTA12-mshift-1;  % "cent" is "centered"; imaxSTA12 is original 1: 2*mshift+1, corresponds to -mshift: mshift
    imaxSTA13cent=imaxSTA13-mshift-1;
    imaxSTA32cent=imaxSTA32-mshift-1;
    %%% NOTICE: the right order of a closed 3-sta pair is +13, -12, +32, where 13 means 1-->3
    iloopoff=imaxSTA13cent-imaxSTA12cent+imaxSTA32cent; %How well does the integer loop close?
    %
    xmaxSTA12n=xmaxSTA12n-mshift-1;
    xmaxSTA13n=xmaxSTA13n-mshift-1;
    xmaxSTA32n=xmaxSTA32n-mshift-1;
    loopoff=xmaxSTA13n-xmaxSTA12n+xmaxSTA32n; %How well does the interpolated loop close?
    xcmaxAVEn=(xcmaxSTA12n+xcmaxSTA13n+xcmaxSTA32n)/3;     % arithmetic average, == x-corr max average normalized
    % xcnshifts=cputime-t
    % t=cputime;
    %%% ampmax == max amplitude among all windows of all 3 2-station pairs
    ampmax=max([ampSTA12; ampSTA13; ampSTA32]);  % ';' means another row, size of the concatenate is 3*nwin, 2*mshift+1
    medxcmaxAVEn=median(xcmaxAVEn);
    xmaxSTA12ntmp=xmaxSTA12n;    % tmp == temporary
    xmaxSTA13ntmp=xmaxSTA13n;
    xmaxSTA32ntmp=xmaxSTA32n;
    
    %% find the strongest 0.5s window with main arrival
    iup=4;
    nin=0;      % subscript flag, to count the successful detections
    zerosallowed=20*winlen/160;   % 20*12.5*sps/160, BUT what is the meaning of this eqn, 5 zeros /sec == 5*12.5 /win ???
    concentration=0.5; %in seconds; how concentrated is the coherent energy within the window?
    cncntr=concentration*sps;   % in samples, 20
    offset=round(0.5*cncntr);   % +- 1/2*cncntr, 10 samples, 0.25s
    
    Roo = -1*ones(1,nwin);  % R_o/o in Peng et al. 2015
    for n=1:nwin
      %         n=1;
      %%%%%%%%%%%%%%% Detection Rejection Criteria %%%%%%%%%%%%%%%%%%%%%%
      % 1. if < min threshold == 0.3, or
      % 2. abs(loopoff) > loopoffmax == 2.1, or
      % 3. imaxSTA12/13/32 is 1 or 2*mshift+1, which is located at the edge of the range [1, 2*mshift+1], or
      % 4. too much zeros in the trace
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% xmaxSTA13n(n)-xmaxSTA12n(n)+xmaxSTA32n(n) == loopoff(n)
      if xcmaxAVEn(n)<xcmaxAVEnmin || abs(xmaxSTA13n(n)-xmaxSTA12n(n)+xmaxSTA32n(n))>loopoffmax ...
          || isequal(abs(imaxSTA12cent(n)),mshift) || isequal(abs(imaxSTA13cent(n)),mshift) ...
          || isequal(abs(imaxSTA32cent(n)),mshift) || max(STAnzeros(:,n))>zerosallowed
        xmaxSTA12ntmp(n)=mshift+1; xmaxSTA13ntmp(n)=mshift+1; xmaxSTA32ntmp(n)=mshift+1; %dummy them, if these criteria are met
      else
%%%%%%%%%%%%% START, USE THIS SECTION IF USING UPSAMPLING %%%%%%%%%%%%%%%%%%%%%%%            
            interpSTA12n=interp(sumsSTA12n(n,:),iup,3);
            interpSTA13n=interp(sumsSTA13n(n,:),iup,3);
            interpSTA32n=interp(sumsSTA32n(n,:),iup,3);
            leninterp=length(interpSTA12n);
            [xcmaxinterpSTA12n,imaxinterpSTA12]=max(interpSTA12n(1:leninterp-(iup-1)));
            [xcmaxinterpSTA13n,imaxinterpSTA13]=max(interpSTA13n(1:leninterp-(iup-1)));
            [xcmaxinterpSTA32n,imaxinterpSTA32]=max(interpSTA32n(1:leninterp-(iup-1)));
            xcmaxconprev=-99999.;  %used to be 0; not good with glitches
            for iSTA12=max(1,imaxinterpSTA12-3*iup):...
                       min(imaxinterpSTA12+3*iup,iup*(2*mshift+1)-(iup-1))
            %3 samples from peak;%intentionally wider than acceptable;%iup-1 are extrapolated points
                for iSTA13=max(1,imaxinterpSTA13-3*iup):...
                           min(imaxinterpSTA13+3*iup,iup*(2*mshift+1)-(iup-1))
                    ibangon = (iup*mshift+1)-iSTA13+iSTA12;
                    if ibangon >= 1 && ibangon<=iup*(2*mshift+1)
                        xcmaxcon=interpSTA12n(iSTA12)+interpSTA13n(iSTA13)+interpSTA32n(ibangon);
                        if xcmaxcon > xcmaxconprev
                            xcmaxconprev=xcmaxcon;
                            iSTA12bang=iSTA12;
                            iSTA13bang=iSTA13;
                        end
                    end
                end
            end
%%%%%%%%%% END, USE THIS SECTION IF USING UPSAMPLING %%%%%%%%%%%%%%%%%%%%%%% 
        %%% will result in the max xcmaxcon and corresponding iSTA12,
        %%% iSTA13, and save them into xcmaxconprev, iSTA12bang and iSTA13bang
        
        iSTA32bang=(iup*mshift+1)-iSTA13bang+iSTA12bang;
        if abs(iSTA12bang-imaxinterpSTA12) <= loopoffmax*iup && ...
            abs(iSTA13bang-imaxinterpSTA13) <= loopoffmax*iup && ...
            abs(iSTA32bang-imaxinterpSTA32) <= loopoffmax*iup && ...
            interpSTA12n(iSTA12bang)+interpSTA13n(iSTA13bang)+interpSTA32n(iSTA32bang) >= ...
            3*xcmaxAVEnmin
          %%% ALSO, sumsSTA12n(n,iSTA12bang) == sumsSTA12nn(iSTA12bang)
          
          xmaxSTA12ntmp(n)=(iSTA12bang-(iup*mshift+1))/iup;
          xmaxSTA13ntmp(n)=(iSTA13bang-(iup*mshift+1))/iup;
          xmaxSTA32ntmp(n)=(iSTA32bang-(iup*mshift+1))/iup;
                    
          %%% let us assume xmaxSTA12ntmp and so is the mean, mu of the distribution of the
          %%% real value, and now we need some estimates of standard deviation, sigma
          % 1st kind of estimate, the deviation between forced circuit and real circuit and
          % then take it as the same to both 12 and 13
          tmp12 = (imaxinterpSTA12-(iup*mshift+1));
          tmp13 = (imaxinterpSTA13-(iup*mshift+1));
          tmp32 = (imaxinterpSTA32-(iup*mshift+1));
          sigma(nin+1) = abs(tmp13+tmp32-tmp12);
          
          % 2nd kind of estimate, the deviation between CC shift and forced shift in each
          % pair, 12 and 13, so sigma is different for each pair
          sigma12(nin+1) = abs(iSTA12bang-imaxinterpSTA12);
          sigma13(nin+1) = abs(iSTA13bang-imaxinterpSTA13);
          
          %for plotting traces
          imaxSTA12wr=round(xmaxSTA12ntmp(n)); %without interpolation this is not needed.
          imaxSTA13wr=round(xmaxSTA13ntmp(n));
          %
          istart=igstart+(n-1)*winoff; %+mshift; %a better way might exist?  %ADDED mshift 10/20/12; DELETED IT 1/19/17.
          %ADDED IT BACK 10/4/2017 to fix bug.  PGC is offset from igstart by mshift before first x-correlation.
          %Not sure why mshift was added.  It changes STA12tr, STA1file etc. relative to the window that was used
          %in the original x-correlation.  This will affect the stated time of max energy (through idiff).
          %GOT RID of the mshift, yet again, 6/29/18, but only after subtracing mshift from all those istarts and
          %iends in lines 342-355.
          iend=istart+winlen-1;
          imid=round((istart+iend)/2);
          %Check power spectrum for reasonableness
          %%% pwelch is a built-in function, [Pxx F] = pwelch(X, WINDOW, NOVERLAP, NFFT, Fs)
          [STA1xx fp] = pwelch(STAopt(1,istart:iend),[],[],[],sps); %40 is sps
          STA1xx=STA1xx/max(STA1xx);    % normalization
          [STA2xx fp] = pwelch(STAopt(2,istart-imaxSTA12wr:iend-imaxSTA12wr),[],[],[],sps);  % WHY substract imaxSTA12wr ???
          STA2xx=STA2xx/max(STA2xx);
          [STA3xx fp] = pwelch(STAopt(3,istart-imaxSTA13wr:iend-imaxSTA13wr),[],[],[],sps);
          STA3xx=STA3xx/max(STA3xx);
          flo=find(fp > lo,1)-1;    % find(fp > lo,1) finds the first 1 indice that satisfies fp > lo
          fhi=find(fp > hi,1)+1;    %extra 1 for good measure
          belowcut=median([STA1xx(2:flo); STA2xx(2:flo); STA3xx(2:flo)]);
          ppeaksSTA1=findpeaks(STA1xx(flo+1:fhi));   % PKS = findpeaks(Y) finds local peaks in the data vector Y
          if length(ppeaksSTA1)>=1
            maxppeakSTA1=max(ppeaksSTA1);
          else
            maxppeakSTA1=0.;
          end
          ppeaksSTA2=findpeaks(STA2xx(flo+1:fhi));   % for STA2, use exactly the same procedure as STA1
          if length(ppeaksSTA2)>=1
            maxppeakSTA2=max(ppeaksSTA2);
          else
            maxppeakSTA2=0.;
          end
          ppeaksSTA3=findpeaks(STA3xx(flo+1:fhi));   % for STA3, still the same
          if length(ppeaksSTA3)>=1
            maxppeakSTA3=max(ppeaksSTA3);
          else
            maxppeakSTA3=0.;
          end
          abovecut=median([maxppeakSTA1 maxppeakSTA2 maxppeakSTA3]);   % relative to belowcut, remain [belowcut, abovecut]
          if abovecut > 0.9*belowcut %-1 %This checks for frequency range; make sure it's not too narrow?
            STA12tr=STAopt(1,istart:iend).*STAopt(2,istart-imaxSTA12wr:iend-imaxSTA12wr);   % see line 554, imax is already centered
            STA13tr=STAopt(1,istart:iend).*STAopt(3,istart-imaxSTA13wr:iend-imaxSTA13wr);   % see line 556
            STA32tr=STAopt(3,istart-imaxSTA13wr:iend-imaxSTA13wr).* ...
              STAopt(2,istart-imaxSTA12wr:iend-imaxSTA12wr);
            cumsumtr=cumsum(STA12tr)+cumsum(STA13tr)+cumsum(STA32tr);    % sum of the cumsum of all traces
            %%% first get the squared sum of each 0.5s window, then get the maximum and start indice
            %%% this gives the proxy of Energy of 0.5s window, i.e. the integral of \dot(E)
            %%% of eqn 1 in Rubin et al. 2013
            %%% idiff<==>win [idiff+1, idiff+20] in regional ind<==>win [istart+idiff, istart+idiff+19] in global index
            [cumsumtrdiff, idiff]=max(cumsumtr(cncntr+1:winlen)-cumsumtr(1:winlen-cncntr));
            
            %%%% Added by Chao, implemention of R_o/o in Peng et al. 2015
            imdiff = istart+idiff+offset;
            isdiff = imdiff-cncntr; %Start of strongest 0.5s, DELETE -1 by Chao, 2019/02/17, see NOTES
            iediff = imdiff+cncntr-1;  % 2*cncntr is 40sps ~= 1s in current hf
            
            %                     isdiff = istart+idiff; %Start of strongest 0.5s, DELETE -1 by Chao, 2019/02/17, see NOTES
            %                     iediff = istart+idiff-1+2*cncntr;
            
            STA11opt = STAopt(1,isdiff:iediff).* STAopt(1,isdiff:iediff);
            STA22opt = STAopt(2,isdiff-imaxSTA12wr:iediff-imaxSTA12wr).* ...
              STAopt(2,isdiff-imaxSTA12wr:iediff-imaxSTA12wr);
            STA33opt = STAopt(3,isdiff-imaxSTA13wr:iediff-imaxSTA13wr).* ...
              STAopt(3,isdiff-imaxSTA13wr:iediff-imaxSTA13wr);
            STAsumopt = STA11opt+STA22opt+STA33opt;
            STA11ort = STAort(1,isdiff:iediff).* STAort(1,isdiff:iediff);
            STA22ort = STAort(2,isdiff-imaxSTA12wr:iediff-imaxSTA12wr).* ...
              STAort(2,isdiff-imaxSTA12wr:iediff-imaxSTA12wr);
            STA33ort = STAort(3,isdiff-imaxSTA13wr:iediff-imaxSTA13wr).* ...
              STAort(3,isdiff-imaxSTA13wr:iediff-imaxSTA13wr);
            STAsumort = STA11ort+STA22ort+STA33ort;
            Eopt = sum(STAsumopt);  % in hf, it is ~ 1s interval
            Eort = sum(STAsumort);
            Roo(n) = Eopt/Eort;
            
            Roomin = 3;
            %%%% Added by Chao, implemention of R_o/o in Peng et al. 2015
            
            if Roo(n) >= Roomin
              
              %%% xcmaxAVEnbang is added by Chao, to distinguish from
              %%% xcmaxAVEn, because it is max average CC coef
              xcmaxAVEnbang(nin+1)=(sumsSTA12n(n,ceil(iSTA12bang/iup))+ ...
                sumsSTA13n(n,ceil(iSTA13bang/iup))+ ...
                sumsSTA32n(n,ceil(iSTA32bang/iup)))/3;
              
              %                         %%% 2020/08/27, i just realize that the follow way to express the average CC
              %                         %%% coeff. is actually more accurate than above, but let's not change it, i
              %                         %%% assume they should not be different too much
              %                         xcmaxAVEnbang(nin+1)=(interpSTA12n(iSTA12bang)+...
              %                                               interpSTA13n(iSTA13bang)+...
              %                                               interpSTA32n(iSTA32bang))/3;
              
              %                         % display the number of overlapping window, easy for testing
              %                         disp(n)
              %                         keyboard

            %What is amp squared in strongest coherent 1/2 sec?
            %%% amp squared is a self to self operation
            %%% NOTICE!   HERE, istart=igstart+(n-1)*winoff >= igstart,
            %%% So, isdiff >= istart >= igstart
            isdiff=istart+idiff; %Start of strongest 0.5s, DELETE -1 by Chao, 2019/02/17, see NOTES
            iediff=istart+idiff-1+cncntr;  % cncntr is 20sps == 0.5s
            dummy=STAopt(1,isdiff:iediff).^2+ ...   % point square
              STAopt(2,isdiff-imaxSTA12wr:iediff-imaxSTA12wr).^2+ ...
              STAopt(3,isdiff-imaxSTA13wr:iediff-imaxSTA13wr).^2;
            dum2=cumsum(dummy);
            %                     Ampsq(nin+1)=dum2(end);   % Ampsq(1) == amplitude square = cumsum of all dummy = squared sum
            Ampsq(nin+1)=dum2(end) / length(dummy);
            
            %%Energy in prior 2.5*cncntr seconds, with offset (assuming 0.5cncntr)
            %%% energy is a self to self operation, and is proportional to amp squared
            %%% the offset could be regarded as a buffer zone
            if isdiff > round(2.5*cncntr)+(mshift-cyclskip)+offset
              dummy=STAopt(1,isdiff-round(2.5*cncntr)-offset:isdiff-offset-1).^2+ ...
                STAopt(2,isdiff-round(2.5*cncntr)-imaxSTA12wr-offset:...
                isdiff-imaxSTA12wr-offset-1).^2+ ...
                STAopt(3,isdiff-round(2.5*cncntr)-imaxSTA13wr-offset:...
                isdiff-imaxSTA13wr-offset-1).^2;
            else
              dummy=STAopt(1,(mshift-cyclskip):isdiff-1-offset).^2+ ...
                STAopt(2,(mshift-cyclskip)-imaxSTA12wr:...
                isdiff-imaxSTA12wr-1-offset).^2+ ...
                STAopt(3,(mshift-cyclskip)-imaxSTA13wr:...
                isdiff-imaxSTA13wr-1-offset).^2;
            end
            dum2=cumsum(dummy);
            %                     Prev(nin+1)=dum2(end);    % Prev(1) == previous 1.25s window before the strongest window, length 2.5*cncntr
            Prior(nin+1)=dum2(end) / length(dummy);
            clear dummy
            
            %CC in same window (test)
            %%% cc is a cross-station operation
            dummy(1,:)=STAopt(1,istart:iend);
            dummy(2,:)=STAopt(2,istart-imaxSTA12wr:iend-imaxSTA12wr);
            dummy(3,:)=STAopt(3,istart-imaxSTA13wr:iend-imaxSTA13wr);
            denoms=dot(dummy,dummy,2);    % dot(A,B,DIM) returns the summed scalar product of A and B in the dimension DIM, 2 is row
            cc(nin+1)=(dot(dummy(1,:),dummy(2,:))/sqrt(denoms(1)*denoms(2))+...
              dot(dummy(2,:),dummy(3,:))/sqrt(denoms(2)*denoms(3))+ ...
              dot(dummy(3,:),dummy(1,:))/sqrt(denoms(3)*denoms(1)))/3;
            clear dummy
            %                     cc(nin+1)=0;
            
            %CC in prior 2.5*cncntr seconds, with offset
            if isdiff > round(2.5*cncntr)+(mshift-cyclskip)+offset
              dummy(1,:)=STAopt(1,isdiff-round(2.5*cncntr)-offset:isdiff-1-offset);   % 2.5*cnctr == 1.25s *sps
              dummy(2,:)=STAopt(2,isdiff-imaxSTA12wr-round(2.5*cncntr)-offset:...
                isdiff-imaxSTA12wr-1-offset);
              dummy(3,:)=STAopt(3,isdiff-imaxSTA13wr-round(2.5*cncntr)-offset:...
                isdiff-imaxSTA13wr-1-offset);
            else
              dummy(1,:)=STAopt(1,(mshift-cyclskip):isdiff-1-offset);
              dummy(2,:)=STAopt(2,(mshift-cyclskip)-imaxSTA12wr:...
                isdiff-imaxSTA12wr-1-offset);
              dummy(3,:)=STAopt(3,(mshift-cyclskip)-imaxSTA13wr:...
                isdiff-imaxSTA13wr-1-offset);
            end
            denoms=dot(dummy,dummy,2);   % dot(A,B,DIM) returns the summed scalar product of A and B in the dimension DIM, 2 is row
            ccprior(nin+1)=(dot(dummy(1,:),dummy(2,:))/sqrt(denoms(1)*denoms(2))+...
              dot(dummy(2,:),dummy(3,:))/sqrt(denoms(2)*denoms(3))+ ...
              dot(dummy(3,:),dummy(1,:))/sqrt(denoms(3)*denoms(1)))/3;  % ccprior125(1) means cc in prior 1.25s
            clear dummy
            %                         ccprior(nin+1)=0;
            
            
            %CC in prior 4 seconds, with offset
            if isdiff > 4*sps+(mshift-cyclskip)+offset   % >winsec*sps+19+10
              dummy(1,:)=STAopt(1,isdiff-4*sps-offset:isdiff-1-offset);   % this is 12.5s much more than 4s, maybe the previous winlen is 4s
              dummy(2,:)=STAopt(2,isdiff-imaxSTA12wr-4*sps-offset:...
                isdiff-imaxSTA12wr-1-offset);
              dummy(3,:)=STAopt(3,isdiff-imaxSTA13wr-4*sps-offset:...
                isdiff-imaxSTA13wr-1-offset);
            else
              dummy(1,:)=STAopt(1,(mshift-cyclskip):isdiff-1-offset);   % 19: isdiff-1    % MIGHT this part be a mistake? no symmetry, with offset
              dummy(2,:)=STAopt(2,(mshift-cyclskip)-imaxSTA12wr:...
                isdiff-imaxSTA12wr-1-offset);
              dummy(3,:)=STAopt(3,(mshift-cyclskip)-imaxSTA13wr:...
                isdiff-imaxSTA13wr-1-offset);
            end
            denoms=dot(dummy,dummy,2);
            ccprior4(nin+1)=(dot(dummy(1,:),dummy(2,:))/sqrt(denoms(1)*denoms(2))+...
              dot(dummy(2,:),dummy(3,:))/sqrt(denoms(2)*denoms(3))+ ...
              dot(dummy(3,:),dummy(1,:))/sqrt(denoms(3)*denoms(1)))/3;
            clear dummy
            %                         ccprior4(nin+1)=0;
            
            
            %CC in following 2.5*cncntr seconds, with offset
            if iediff+round(2.5*cncntr)+(mshift-cyclskip)+offset <= size(STAopt,2)
              dummy(1,:)=STAopt(1,iediff+1+offset:iediff+round(2.5*cncntr)+offset);  % 1.25s win after 0.5s after 2*offset
              dummy(2,:)=STAopt(2,iediff+1-imaxSTA12wr+offset:...
                iediff-imaxSTA12wr+round(2.5*cncntr)+offset);
              dummy(3,:)=STAopt(3,iediff+1-imaxSTA13wr+offset:...
                iediff-imaxSTA13wr+round(2.5*cncntr)+offset);
            else
              dummy(1,:)=STAopt(1,iediff+1+offset: size(STAopt,2)-(mshift-cyclskip));  % 1.25s win after 0.5s after 2*offset
              dummy(2,:)=STAopt(2,iediff+1-imaxSTA12wr+offset:...
                size(STAopt,2)-(mshift-cyclskip)-imaxSTA12wr);
              dummy(3,:)=STAopt(3,iediff+1-imaxSTA13wr+offset:...
                size(STAopt,2)-(mshift-cyclskip)-imaxSTA13wr);
            end
            denoms=dot(dummy,dummy,2);
            ccpost(nin+1)=(dot(dummy(1,:),dummy(2,:))/sqrt(denoms(1)*denoms(2))+...
              dot(dummy(2,:),dummy(3,:))/sqrt(denoms(2)*denoms(3))+ ...
              dot(dummy(3,:),dummy(1,:))/sqrt(denoms(3)*denoms(1)))/3;
            clear dummy
            %                         ccpost(nin+1)=0;
            
            %CC in following 4 seconds, with offset
            if iediff+4*sps+(mshift-cyclskip)+offset <= size(STAopt,2)
              dummy(1,:)=STAopt(1,iediff+1+offset:iediff+4*sps+offset);  % 1.25s win after 0.5s after 2*offset
              dummy(2,:)=STAopt(2,iediff+1-imaxSTA12wr+offset:...
                iediff-imaxSTA12wr+4*sps+offset);
              dummy(3,:)=STAopt(3,iediff+1-imaxSTA13wr+offset:...
                iediff-imaxSTA13wr+4*sps+offset);
            else
              dummy(1,:)=STAopt(1,iediff+1+offset: size(STAopt,2)-(mshift-cyclskip));  % 1.25s win after 0.5s after 2*offset
              dummy(2,:)=STAopt(2,iediff+1-imaxSTA12wr+offset:...
                size(STAopt,2)-(mshift-cyclskip)-imaxSTA12wr);
              dummy(3,:)=STAopt(3,iediff+1-imaxSTA13wr+offset:...
                size(STAopt,2)-(mshift-cyclskip)-imaxSTA13wr);
            end
            denoms=dot(dummy,dummy,2);
            ccpost4(nin+1)=(dot(dummy(1,:),dummy(2,:))/sqrt(denoms(1)*denoms(2))+...
              dot(dummy(2,:),dummy(3,:))/sqrt(denoms(2)*denoms(3))+ ...
              dot(dummy(3,:),dummy(1,:))/sqrt(denoms(3)*denoms(1)))/3;
            clear dummy
            
            %Energy in prior winlen, i.e. 16 s in current setting
            if isdiff > winlen+(mshift-cyclskip)+offset
              dummy=STAopt(1,isdiff-winlen-offset:isdiff-offset-1).^2+ ...
                STAopt(2,isdiff-winlen-imaxSTA12wr-offset:...
                isdiff-imaxSTA12wr-offset-1).^2+ ...
                STAopt(3,isdiff-winlen-imaxSTA13wr-offset:...
                isdiff-imaxSTA13wr-offset-1).^2;
            else
              dummy=STAopt(1,(mshift-cyclskip):isdiff-1-offset).^2+ ...
                STAopt(2,(mshift-cyclskip)-imaxSTA12wr:...
                isdiff-imaxSTA12wr-1-offset).^2+ ...
                STAopt(3,(mshift-cyclskip)-imaxSTA13wr:...
                isdiff-imaxSTA13wr-1-offset).^2;
            end
            dum2=cumsum(dummy);
            %                     Prior16(nin+1)=dum2(end);
            Prior16(nin+1)=dum2(end) / length(dummy);
            clear dummy
            
            %Energy in prior half of winlen, i.e. 8 s in current setting
            if isdiff > winlen/2 +(mshift-cyclskip)+offset
              dummy=STAopt(1,isdiff-winlen/2 -offset:isdiff-offset-1).^2+ ...
                STAopt(2,isdiff-winlen/2 -imaxSTA12wr-offset:...
                isdiff-imaxSTA12wr-offset-1).^2+ ...
                STAopt(3,isdiff-winlen/2 -imaxSTA13wr-offset:...
                isdiff-imaxSTA13wr-offset-1).^2;
            else
              dummy=STAopt(1,(mshift-cyclskip):isdiff-1-offset).^2+ ...
                STAopt(2,(mshift-cyclskip)-imaxSTA12wr:...
                isdiff-imaxSTA12wr-1-offset).^2+ ...
                STAopt(3,(mshift-cyclskip)-imaxSTA13wr:...
                isdiff-imaxSTA13wr-1-offset).^2;
            end
            dum2=cumsum(dummy);
            %                     Prior8(nin+1)=dum2(end);
            Prior8(nin+1)=dum2(end) / length(dummy);
            clear dummy
            
            %Energy in prior 4 s
            if isdiff > 4*sps +(mshift-cyclskip)+offset
              dummy=STAopt(1,isdiff-4*sps -offset:isdiff-offset-1).^2+ ...
                STAopt(2,isdiff-4*sps -imaxSTA12wr-offset:...
                isdiff-imaxSTA12wr-offset-1).^2+ ...
                STAopt(3,isdiff-4*sps -imaxSTA13wr-offset:...
                isdiff-imaxSTA13wr-offset-1).^2;
            else
              dummy=STAopt(1,(mshift-cyclskip):isdiff-1-offset).^2+ ...
                STAopt(2,(mshift-cyclskip)-imaxSTA12wr:...
                isdiff-imaxSTA12wr-1-offset).^2+ ...
                STAopt(3,(mshift-cyclskip)-imaxSTA13wr:...
                isdiff-imaxSTA13wr-1-offset).^2;
            end
            dum2=cumsum(dummy);
            %                     Prior4(nin+1)=dum2(end);
            Prior4(nin+1)=dum2(end) / length(dummy);
            clear dummy
            
            % energy (amp squared sum) in post 2.5*cncntr second window, with offset
            if iediff+round(2.5*cncntr)+(mshift-cyclskip)+offset <= size(STAopt,2)
              dummy=STAopt(1,iediff+1+offset:iediff+round(2.5*cncntr)+offset).^2+ ...
                STAopt(2,iediff+1-imaxSTA12wr+offset:...
                iediff-imaxSTA12wr+round(2.5*cncntr)+offset).^2+ ...
                STAopt(3,iediff+1-imaxSTA13wr+offset:...
                iediff-imaxSTA13wr+round(2.5*cncntr)+offset).^2;
            else
              dummy=STAopt(1,iediff+1+offset: size(STAopt,2)-(mshift-cyclskip)).^2+ ...
                STAopt(2,iediff+1-imaxSTA12wr+offset:...
                size(STAopt,2)-(mshift-cyclskip)-imaxSTA12wr).^2+ ...
                STAopt(3,iediff+1-imaxSTA13wr+offset:...
                size(STAopt,2)-(mshift-cyclskip)-imaxSTA13wr).^2;
            end
            dum2=cumsum(dummy);
            %                     Post(nin+1)=dum2(end);
            Post(nin+1)=dum2(end) / length(dummy);
            clear dummy
            
            % energy (amp squared sum) in post 4 s window, with offset
            if iediff+4*sps+(mshift-cyclskip)+offset <= size(STAopt,2)
              dummy=STAopt(1,iediff+1+offset:iediff+4*sps+offset).^2+ ...
                STAopt(2,iediff+1-imaxSTA12wr+offset:...
                iediff-imaxSTA12wr+4*sps+offset).^2+ ...
                STAopt(3,iediff+1-imaxSTA13wr+offset:...
                iediff-imaxSTA13wr+4*sps+offset).^2;
            else
              dummy=STAopt(1,iediff+1+offset: size(STAopt,2)-(mshift-cyclskip)).^2+ ...
                STAopt(2,iediff+1-imaxSTA12wr+offset:...
                size(STAopt,2)-(mshift-cyclskip)-imaxSTA12wr).^2+ ...
                STAopt(3,iediff+1-imaxSTA13wr+offset:...
                size(STAopt,2)-(mshift-cyclskip)-imaxSTA13wr).^2;
            end
            dum2=cumsum(dummy);
            %                     Post4(nin+1)=dum2(end);
            Post4(nin+1)=dum2(end) / length(dummy);
            clear dummy
            
            % energy (amp squared sum) in post winlen/2 window, i.e. 8 s, with offset
            if iediff+winlen/2+(mshift-cyclskip)+offset <= size(STAopt,2)
              dummy=STAopt(1,iediff+1+offset:iediff+winlen/2+offset).^2+ ...
                STAopt(2,iediff+1-imaxSTA12wr+offset:...
                iediff-imaxSTA12wr+winlen/2+offset).^2+ ...
                STAopt(3,iediff+1-imaxSTA13wr+offset:...
                iediff-imaxSTA13wr+winlen/2+offset).^2;
            else
              dummy=STAopt(1,iediff+1+offset: size(STAopt,2)-(mshift-cyclskip)).^2+ ...
                STAopt(2,iediff+1-imaxSTA12wr+offset:...
                size(STAopt,2)-(mshift-cyclskip)-imaxSTA12wr).^2+ ...
                STAopt(3,iediff+1-imaxSTA13wr+offset:...
                size(STAopt,2)-(mshift-cyclskip)-imaxSTA13wr).^2;
            end
            dum2=cumsum(dummy);
            %                     Post8(nin+1)=dum2(end);
            Post8(nin+1)=dum2(end) / length(dummy);
            clear dummy
            
            % energy (amp squared sum) in post winlen window, i.e. 16 s, with offset
            if iediff+winlen+(mshift-cyclskip)+offset <= size(STAopt,2)
              dummy=STAopt(1,iediff+1+offset:iediff+winlen+offset).^2+ ...
                STAopt(2,iediff+1-imaxSTA12wr+offset:...
                iediff-imaxSTA12wr+winlen+offset).^2+ ...
                STAopt(3,iediff+1-imaxSTA13wr+offset:...
                iediff-imaxSTA13wr+winlen+offset).^2;
            else
              dummy=STAopt(1,iediff+1+offset: size(STAopt,2)-(mshift-cyclskip)).^2+ ...
                STAopt(2,iediff+1-imaxSTA12wr+offset:...
                size(STAopt,2)-(mshift-cyclskip)-imaxSTA12wr).^2+ ...
                STAopt(3,iediff+1-imaxSTA13wr+offset:...
                size(STAopt,2)-(mshift-cyclskip)-imaxSTA13wr).^2;
            end
            dum2=cumsum(dummy);
            %                     Post16(nin+1)=dum2(end);
            Post16(nin+1)=dum2(end) / length(dummy);
            clear dummy
            
            
            %                     %%%%%%% Below is for comparison to template, + and -, start
            %                     %
            %                     %%% isdiff=istart+idiff-1, so indx == isdiff+round(cncntr/2)
            %                     indx=istart+idiff-1+round(cncntr/2); %indx should be centered ~ on zero-crossing of main arrival
            %                     traces(1,:)=STAopt(1,indx-tempbef:indx+tempaft);   % -59;+60
            %                     traces(2,:)=STAopt(2,indx-tempbef-imaxSTA12wr:indx+tempaft-imaxSTA12wr);
            %                     traces(3,:)=STAopt(3,indx-tempbef-imaxSTA13wr:indx+tempaft-imaxSTA13wr);
            %                     for ista=1:nsta
            %                         % tempxc == template x-corr, |max lag| == 0.25s == 10 samples
            %                         % cc between template and the window with the main arrival
            %                         tempxc(ista,:)=xcorr(traces(ista,:),STAtemps(ista,:),floor(cncntr/2),'coeff');
            %                     end
            %                     sumxc=sum(tempxc)/nsta;    % sum along column, average cc coeff of all 3 stats at each shift
            %                     [match(nin+1,1), ioff]=max(sumxc);  % get the max coeff and the index of the shift
            %                     %%% the follow line is to convert index to real shift value
            %                     ioff=ioff-(floor(cncntr/2)+1); %shift STAtemps by match(nin+1) (shift right for positive values)
            %                     %%% timstemp == times of template, 1: 120, 121: 240
            %                     timstemp(nin*templen+1:(nin+1)*templen)=timsSTA(indx-tempbef+ioff:indx+tempaft+ioff);
            %                     for ista=1:nsta
            %                         % tempxcneg means x-corr with negative template
            %                         tempxcneg(ista,:)=xcorr(traces(ista,:),-STAtemps(ista,:),floor(cncntr/2),'coeff');
            %                     end
            %                     sumxcneg=sum(tempxcneg)/nsta;
            %                     [match(nin+1,2), ioff]=max(sumxcneg);
            %                     ioff=ioff-(floor(cncntr/2)+1); %shift "snips" by match(nin+1) (shift right for positive values), what is snips ???
            %                     timstempneg(nin*templen+1:(nin+1)*templen)=timsSTA(indx-tempbef+ioff:indx+tempaft+ioff);
            %                     clear traces
            %                     %%%%%%% For comparison, end
            match(nin+1,1)=0;
            match(nin+1,2)=0;
            
            % STA1file, STA2file, STA3file are back-to-back windows for the day (plus sample times)
            % 1st column: time
            % 2nd column: data
            STA1file(nin*winlen+1:(nin+1)*winlen,1:2)=[timsSTA(istart:iend)' ...
              STAopt(1,istart:iend)'];
            STA2file(nin*winlen+1:(nin+1)*winlen,1:2)=[timsSTA(istart:iend)' ...
              STAopt(2,istart-imaxSTA12wr:iend-imaxSTA12wr)'];
            STA3file(nin*winlen+1:(nin+1)*winlen,1:2)=[timsSTA(istart:iend)' ...
              STAopt(3,istart-imaxSTA13wr:iend-imaxSTA13wr)'];
            STAamp(nin+1,1)=prctile(abs(STAopt(1,istart:iend)),80);  % prctile,  Percentiles of a sample, SEE NOTES
            STAamp(nin+1,2)=prctile(abs(STAopt(2,istart-imaxSTA12wr:iend-imaxSTA12wr)),80);
            STAamp(nin+1,3)=prctile(abs(STAopt(3,istart-imaxSTA13wr:iend-imaxSTA13wr)),80);
            STAamp(nin+1,:)=STAamp(nin+1,:)/STAamp(nin+1,1);   % normalization relative to 1st stat
            
            %%% imaxSTA12wr     -> cc shift of qualified win only, so remain no changed until the next qualified win;
            %%% xcmaxSTA12n(n)  -> max cc coef in all 4s win
            %%% cumsumtrdiff/cumsumtr(winlen) -> normalized cumsumtrdiff
            STA12file(nin+1,1:2)=[imaxSTA12wr xcmaxSTA12n(n)];
            STA13file(nin+1,1:2)=[imaxSTA13wr xcmaxSTA13n(n)];
            STA32file(nin+1,1:3)=[cumsumtrdiff/cumsumtr(winlen) xcmaxSTA32n(n) idiff];
            
            nin=nin+1;
            istartkeep(nin)=istart; %For adding other stations later, keep the istart of each win
            aSTA12keep(nin,:)=[timswin(n) aSTA12(n)];   % timswin == time at the center of each win
            aSTA13keep(nin,:)=[timswin(n) aSTA13(n)];
            aSTA32keep(nin,:)=[timswin(n) aSTA32(n)];
            loopoffkeep(nin,:)=[timswin(n) loopoff(n)];
            
            %%% 1. Without interpolation, xmaxSTA13ntmp == imaxSTA13wr
            %%% 2. xcmaxAVEn is arithmetic average, == x-corr coeff max average normalized in each win
            %%% 3. Ampsq -> amplitude square == squared sum in each strongest 0.5s win
            %%% 4. cumsumtrdiff -> the maximum squared sum of among all 0.5s windows
            %%% 5. timswin(n)-winlensec/2+idiff/sps -> start time in sec of the strongest 0.5s win
            %%% 6. cumsumtrdiff/cumsumtr(winlen) -> normalized cumsumtrdiff
            %%% 7. match(1) -> max cc coeff between template and the window with the main arrival, 2s win, 120 samples
            %%% 8. match(2) -> max cc coeff between negtive template and the window with the main arrival, 2s win, 120 samples
            %%% 9. Prev -> amp squared sum (energy) in previous 1.25s window before the strongest window
            %%% 10. Post -> amp squared sum (energy) in post 1.25s window before the strongest window
            %%% 11. ccprior -> cc coeff in previous 4 secs window
            %%% 12. STAamp -> normalized 80 percentile of amplitude of data
            %%% 13. xcmaxSTA12n -> max cc coeff in each 4s window
            
            %                     mapfile(nin,:)=[timswin(n) xmaxSTA13ntmp(n) xmaxSTA12ntmp(n) ...
            %                         xcmaxAVEn(n) loopoff(n) Ampsq(nin) cumsumtrdiff timswin(n)-winlensec/2+idiff/sps cumsumtrdiff/cumsumtr(winlen) ...
            %                         match(nin,1) match(nin,2) Prev(nin) Post(nin) Prev15(nin) Prev30(nin) ...
            %                         ccprior125(nin) ccprior(nin) ccpost125(nin) STAamp(nin,2) STAamp(nin,3) xcmaxSTA12n(n) xcmaxSTA13n(n) xcmaxSTA32n(n) ];
            %%% NOTICE: add by Chao, save xcmaxAVEnbang instead of
            %%% xcmaxAVEn, and the other 3 ccmax's
            
            %%% 30 cols, 2021/03/20
            mapfile(nin,:)=[timswin(n) xmaxSTA12ntmp(n) xmaxSTA13ntmp(n) ...
                            xcmaxAVEnbang(nin) loopoff(n) cumsumtrdiff ...
                            timswin(n)-winlensec/2+(idiff+1)/sps ...
                            cumsumtrdiff/cumsumtr(winlen) Ampsq(nin) ...
                            Prior(nin) Post(nin) Prior4(nin) Post4(nin) Prior8(nin) ...
                            Post8(nin) Prior16(nin) Post16(nin) cc(nin) ccprior(nin) ...
                            ccpost(nin) ccprior4(nin) ccpost4(nin)...
                            STAamp(nin,2) STAamp(nin,3) ...
                            sumsSTA12n(n,ceil(iSTA12bang/iup)) ...
                            sumsSTA13n(n,ceil(iSTA13bang/iup)) ...
                            sumsSTA32n(n,ceil(iSTA32bang/iup)) ...
                            sigma(nin) sigma12(nin) sigma13(nin)];
            else
              xmaxSTA12ntmp(n)=mshift+1; xmaxSTA13ntmp(n)=mshift+1; xmaxSTA32ntmp(n)=mshift+1;
            end
          else
            % 20 == mshift+1
            xmaxSTA12ntmp(n)=mshift+1; xmaxSTA13ntmp(n)=mshift+1; xmaxSTA32ntmp(n)=mshift+1;
          end
        else
          xmaxSTA12ntmp(n)=mshift+1; xmaxSTA13ntmp(n)=mshift+1; xmaxSTA32ntmp(n)=mshift+1;
        end
      end
    end
    
    nin
    
    if nin == 0
      fprintf('No detection is found in fam %s in day %s %s, no plots&files would be saved.\n', ...
        fam, YEAR, JDAY);
    else
      
      %%% 2nd figure, distribution of cc coefs and offsets of all windows
      figure(102)
      colormap(jet)
      %%% loopoff == xmaxSTA13n-xmaxSTA12n+xmaxSTA32n
      %%% AmpComp -> difference of (sum of squared amplitude in 4s win) between each win and its previous 4th win
      %%% scatter, 3 is the marker size, AmpComp is the color of the marker
      scatter(xmaxSTA13n-xmaxSTA12n+xmaxSTA32n,xcmaxAVEn,3,AmpComp)
      c=colorbar; 
      hold on
      plot(-50:50,xcmaxAVEnmin+zeros(101,1),'k:');
      axis([min(-5,-2.5*loopoffmax) max(5,2.5*loopoffmax) -0.2 1.0]);
      hrf = plotreflinesr(gca,-loopoffmax,'x','k');  % plots dotted lines at the locations specified
      hrf = plotreflinesr(gca,loopoffmax,'x','k');
      c.Label.String = 'Diff. in sum of amp^2 between consecutive windows';
      ylabel('Average CC');
      xlabel('Summed time offsets between station pairs (samples at 40 Hz)');
      
      box on
      %%% Tex syntax to use subscript, '_{}'
%       title([IDENTIF,'_{',num2str(lo),'-',num2str(hi),'}'])
%       print('-depsc',[rstpath,'/FIGS/up',IDENTIF,'_',num2str(winlen/sps),'s_',num2str(sps),'sps_',...
%         num2str(lo),'-',num2str(hi),'e.eps'])
%       close(102)

      figure(103)
      colormap(jet)
      %%% loopoff == xmaxSTA13n-xmaxSTA12n+xmaxSTA32n
      %%% AmpComp -> difference of (sum of squared amplitude in 4s win) between each win and its previous 4th win
      %%% scatter, 3 is the marker size, AmpComp is the color of the marker
      scatter(xmaxSTA13n-xmaxSTA12n+xmaxSTA32n,Roo,3,AmpComp);
      hold on
      text(0.9,0.9,strcat(num2str(nin),'/',num2str(sum(Roo~=-1))),'Units','normalized',...
        'HorizontalAlignment','right');
      plot(-50:50,Roomin+zeros(101,1),'k:');
      xlim([min(-5,-2.5*loopoffmax) max(5,2.5*loopoffmax)]);
      box on
      
      %% Write results into files
%       % save analytics of all detections (nin), contain double counting
%       fid = fopen([rstpath,'/MAPS/tdecmapup',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
%         num2str(winlen/sps),'s',num2str(sps),'sps'],'w+');
%       fprintf(fid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %10.3f %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f \n',...
%         mapfile(1:nin,:)');
%       fclose(fid);
%       
%       % save time traces of all detections (nin), contain double counting
%       fid = fopen(strcat(rstpath, '/MAPS/tdectraceup_',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
%         num2str(winlen/sps),'s',num2str(sps),'sps'),'w+');
%       tracefile = [STA1file(1:nin*winlen,:) STA2file(1:nin*winlen,2) STA3file(1:nin*winlen,2)];
%       fprintf(fid,'%.4f %.6f %.6f %.6f \n',tracefile');
%       fclose(fid);
%       
%       % save time traces of original rotated velocity seismogram, after bandpassing
%       fid = fopen(strcat(rstpath, '/MAPS/tdectrace1dup_',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
%         num2str(winlen/sps),'s',num2str(sps),'sps'),'w+');
%       oritracefile = [timsSTA' STAopt'];
%       fprintf(fid,'%.4f %.6f %.6f %.6f \n',oritracefile');
%       fclose(fid);
    
      
    %% METHOD 1 to add additional stations for checking 
    %%%%%%%%%%%%%%%%%%%%%%%%%% METHOD 1, Begin %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear STAopt    % clear these two variables is just to save memory space
    clear STAort
    clear STAnzeros
%     stasnew=['LZB  '
%              'TWKB '
%              'MGCB '];  % twkb lzb mgcb
    stasnew=['LZB  '
             'TWKB '
             'MGCB '
             'KLNB '];  % twkb lzb mgcb          
         
    if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
        stasnew(4,:)='KELB ';
    else
        stasnew(4,:)='KLNB ';  % remember to change it back
    end
    
    nstanew=size(stasnew,1);
    
    fileflag = ones(nstanew,1);
    STAopt=zeros(nstanew,nin*winlen);
    for istanew=1:nstanew
        found=0;
        [LIA,idx]=ismember(stasnew(istanew,:),PERMSTA,'rows');
        if LIA
            found=found+LIA;
            if strcmp(PERMSTA(idx,1:3),'PGC')
                fact=1.0e-3;
            elseif strcmp(PERMSTA(idx,1:3),'LZB')
                fact=1.8e-3;
            end
            fname = strcat(prename,'.',PERMSTA(idx,1:3),'..BHE.D.SAC');
            if isfile(fname)    % if have the data file
%                 [opt,ort,nzeros,timsSTAnew]=readperms(prename,PERMSTA,PERMROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
                [opt,ort,nzeros,timsSTAnew]=readpermsv2(prename,PERMSTA,PERMROTS,idx,sps,lo,hi,...
                                            npo,npa,fact,nwin,winlen,winoff,igstart,[]);
            else
                fileflag(istanew) = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for additional station %s in day %s / %s, will be omitted. \n',...
                        PERMSTA(idx,1:3), YEAR, JDAY);
                continue   % continue to next staion
            end         
        end
        
        [LIA,idx]=ismember(stasnew(istanew,:),POLSTA,'rows');
        if LIA
            found=found+LIA; %better be 1
            if year==2003 && jday<213
                fact=7.5e-3;
            else
                fact=1.5e-3; 
            end
            fname = strcat(prename,'.',POLSTA(idx,1:4),'..HHE.D.SAC');
            if isfile(fname)    % if have the data file
%                 [opt,ort,nzeros]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
                [opt,ort,nzeros]=readpolsv2(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,...
                                            nwin,winlen,winoff,igstart,[]);
            else
                fileflag(istanew) = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for additional station %s in day %s / %s, will be omitted. \n',...
                        POLSTA(idx,1:4), YEAR, JDAY);
                continue   % continue to next staion
            end 
        end
%         found=found
        if fileflag(istanew)
            for n=1:nin
                istart=istartkeep(n);   % saved istart of each successful window
                iend=istart+winlen-1;
                STAopt(istanew,(n-1)*winlen+1:n*winlen)=opt(istart:iend)';   % no need to divide scaleseisms here because they are all 1
            end
        end
    end
    
    stasnew=stasnew(fileflag == 1,:);
    nstanew=size(stasnew,1);
    STAopt=STAopt(fileflag == 1,:); 
    
    in=zeros(nstanew,nin);   % flag to indicate whether this win of this new sta is qualified (1)
    ioff=zeros(nstanew,nin);    % some kind of weighted offset according to the ccmax coef
    ccmaxave=zeros(nstanew,nin);
    loff=zeros(nstanew,nin);    % average sum of abs offset difference
    icc12=zeros(nstanew,nin);
    icc13=zeros(nstanew,nin);
    icc23=zeros(nstanew,nin);
    count=zeros(nstanew,1);
    %%% NOTE: for mshiftnew 's performance, 2.2 > 2.0 > 1.5 > 1.0,
    %%% meaningless to increase more
    mshiftnew=ceil(2.2*mshift); %1.8 was used in first-year report   % ceil(a) gives integer >=a, 
    ccavethres = 0.8*xcmaxAVEnmin;  %% average cc threshold
    ccindthres = 0.6*xcmaxAVEnmin;  %% individual cc threshold
    offavethres = 2;    % average differential offset threshold
    offindthres = 3;
%     offavethres = loopoffmax;    % average differential offset threshold
%     offindthres = 1.5*loopoffmax;
    
    for istanew=1:nstanew   % 1st, for each additional station
        
        for n=1:nin     % 2nd, for each trio detection 
            
            istart=(n-1)*winlen+1;  % start from 1
            iend=n*winlen;
            % recall that STA1/2/3file saves the time and shifted seismogram of each station
            %%% do xcorr between 1 new station to each of the 3-sta trio
            [cc1, lag1] = xcorr(STA1file(istart:iend,2),STAopt(istanew,istart:iend),mshiftnew,...
                                'coeff');
            [cc2, lag2] = xcorr(STA2file(istart:iend,2),STAopt(istanew,istart:iend),mshiftnew,...
                                'coeff');
            [cc3, lag3] = xcorr(STA3file(istart:iend,2),STAopt(istanew,istart:iend),mshiftnew,...
                                'coeff');
            % get the max cc coef and index
            [cc1max,icc1]=max(cc1);
            [cc2max,icc2]=max(cc2);
            [cc3max,icc3]=max(cc3);
            ccmaxave(istanew,n)=(cc1max+cc2max+cc3max)/3;
            
            % NOTE:
            % In each deteciton, since sta 1,2,3 is already aligned to
            % denote one point, thus if it can be checked by the 4th station
            % the shifts between 4th station and
            % original 3 should be close enough, i.e. all differential shifts
            % should be close to 0 theoratically. That is why it is a
            % constraint to filter suspicious detections
            icc12(istanew,n) = lag1(icc1)-lag2(icc2);
            icc13(istanew,n) = lag1(icc1)-lag3(icc3);
            icc23(istanew,n) = lag2(icc2)-lag3(icc3);
            
            % loff and lofftest are exactly the same, although different syntax, depends on user habit
            % loff the absolute sum of the offset difference            
            loff(istanew,n)=(abs(icc12(istanew,n))+abs(icc13(istanew,n))+abs(icc23(istanew,n)))/3;
            
            ioff(istanew,n)=round((cc1max * lag1(icc1) + ...
                                   cc2max * lag2(icc2) + ...
                                   cc3max * lag3(icc3)) / ...
                                  (cc1max+cc2max+cc3max));                                                    
            aaa(istanew,n)=(cc1max * lag1(icc1) + ...
                                   cc2max * lag2(icc2) + ...
                                   cc3max * lag3(icc3)) / ...
                                  (cc1max+cc2max+cc3max);              
            %%% NOTE:
            %%% in first year report, the following criteria was used, so did the hf counterpart 
            %%% if abs(icc12)<=2 && abs(icc13)<=2 && abs(icc23)<=2 && ccmaxave(istanew,n)>1*xcmaxAVEnmin
            if abs(icc12(istanew,n))<=offindthres && abs(icc13(istanew,n))<=offindthres && ...
                    abs(icc23(istanew,n))<=offindthres && ...
                    ccmaxave(istanew,n)>=ccavethres && ... 
                    loff(istanew,n)<=offavethres && ...
                    cc1max >= ccindthres && ...
                    cc2max >= ccindthres && ...
                    cc3max >= ccindthres 
                
                in(istanew,n)=1; 
                count(istanew) = count(istanew)+1;                                
                                  
                %Now shift new station (for plotting purposes)
                if ioff(istanew,n)>0    % shift to right -->
                    STAopt(istanew,istart+ioff(istanew,n):iend)= ...
                        STAopt(istanew,istart:iend-ioff(istanew,n));
                    STAopt(istanew,istart:istart+ioff(istanew,n)-1)=0;
                elseif ioff(istanew,n)<0    % shift to left <--
                    STAopt(istanew,istart:iend+ioff(istanew,n))= ...
                        STAopt(istanew,istart-ioff(istanew,n):iend);
                    STAopt(istanew,iend+ioff(istanew,n)+1:iend)=0;                    
                end
                
            end
            
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%% METHOD 1, End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Write checking result into files
    addrstfile = [in(:,1:nin); loff(:,1:nin) ; ioff(:,1:nin); ccmaxave(:,1:nin)]';      % number of columns should be 4*nstanew 
    fid = fopen([rstpath,'/MAPS/tdecmapaddup',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'add'],'w+');
    if nstanew == 3  
        fprintf(fid,'%d %d %d %.4f %.4f %.4f %d %d %d %.3f %.3f %.3f \n',addrstfile');  % 12 col
    elseif nstanew ==4
        fprintf(fid,'%d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f \n',...
            addrstfile');   % 16 col
    end
    fclose(fid);
    
    fid = fopen(strcat(rstpath, '/MAPS/tdectraceaddup_',IDENTIF,'_',num2str(lo),'-',num2str(hi),...
                '_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'add'),'w+');
    addtracefile = [STA1file(1:nin*winlen,1) STAopt(:,1:nin*winlen)'];
    if nstanew == 3
        fprintf(fid,'%.4f %.6f %.6f %.6f \n',addtracefile');    % 4 col
    elseif nstanew == 4
        fprintf(fid,'%.4f %.6f %.6f %.6f %.6f \n',addtracefile');    % 5 col
    end
    fclose(fid);
    

    %% WAY 3 to avoid double counting
    %%% The previous way to consider the duplicates, i.e. the same arrival seen
    %%% by multiple overlapping windows is to view arrivals who are separated by
    %%% less than a threshold as the same arrival, and choose the one with the 
    %%% largest CC coef.
    %%% But now we will use a different way, we choose one that are closer to the
    %%% center of the window, i.e., fairly far away from the start and the end of
    %%% the window. This issue can be seen by 'identify.m' later because some 
    %%% arrivals survived due to the largest CC are located at the end of the window
    %%% so that the col. 6 and 8 would be doubtful as the strongest arrival window
    %%% involved into the computation may not be complete.
    %%%
    %%% Back to use 'RemoveDoubleCounting2' which always keeps the largerst CC
    %%% detection
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WAY 3 start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% construct a more complete result matrix of all nin windows
    %%% structure of allrst: 30+4*nstanew cols, if 4 new stas, then will be 46 cols
    %%% UPDATED at 2021/03/20 
    %%%   1:timswin(n) 2:xmaxSTA12ntmp(n) 3:xmaxSTA13ntmp(n) 4:xcmaxAVEnbang(nin) 5:loopoff(n)
    %%%   6:cumsumtrdiff 7:timswin(n)-winlensec/2+idiff/sps 8:cumsumtrdiff/cumsumtr(winlen)
    %%%   9:Ampsq(nin) 10:Prior(nin) 11:Post(nin) 12:Prior4(nin) 13:Post4(nin) 14:Prior8(nin)
    %%%   15:Post8(nin) 16:Prior16(nin) 17:Post16(nin) 18:cc(nin) 19:ccprior(nin)
    %%%   20:ccpost(nin) 21:ccprior4(nin) 22:ccpost4(nin) 23:STAamp(nin,2) 24:STAamp(nin,3)
    %%%   25:sumsSTA12n(n,iSTA12bang) 26:sumsSTA13n(n,iSTA13bang)
    %%%   27:sumsSTA32n(n,iSTA32bang) 28:sigma(nin) 29:sigma12(nin) 30:sigma13(nin)
    %%%   in(1:nstanew) loff(1:nstanew) ioff(1:nstanew) ccmaxave(1:nstanew)  
    allrst = [mapfile(1:nin, 1:30) addrstfile];
    dtmin = concentration;      % min time during which only 1 detection is retained
    indext = ceil(winlensec/winoffsec)+4;   % index extension, the index range that needs to be checked
    if nstanew == 3
        colnum = [7 1 4 40 41 42];    % col numbers that are timing and cc coefs
    elseif nstanew == 4
        colnum = [7 1 4 43 44 45 46];
    end
    allrstsort = sortrows(allrst, [7,1]);
    [allrsttemp1,~,~] = RemoveDoubleCounting2(allrstsort,dtmin,indext,colnum,winlensec);
    % do this twice to avoid extreme case like 20.2 20.5 20.5 20.6 20.9,
    % when first round saves the 2nd 20.5 and 20.9, but 20.5 and 20.9 are
    % still duplicates
    [allrsttemp2,~,~] = RemoveDoubleCounting2(allrsttemp1,dtmin,indext,colnum,winlensec);
    bbb = allrsttemp2(2:end,7)-allrsttemp2(1:end-1,7);
    ccc = allrsttemp2(3:end,7)-allrsttemp2(1:end-2,7);
    if ~isempty(find(bbb<=dtmin,1)) || ~isempty(find(ccc<=dtmin,1))
        disp('not enough');
        [allrst_new,~,~] = RemoveDoubleCounting2(allrsttemp2,dtmin,indext,colnum,winlensec);
    else
        allrst_new = allrsttemp2;
    end
    [~, isave, ~] = intersect(allrst, allrst_new, 'row', 'stable');     % saved indexes
    idis = setdiff(1:nin,isave);    % discarded indexes
    
    %%% re-assign those variables that need to be saved
    ind_ori=find(xmaxSTA12ntmp ~= mshift+1);
    idis_ori=ind_ori(idis);
    xmaxSTA12ntmp(idis_ori) = mshift+1;
    xmaxSTA13ntmp(idis_ori) = mshift+1;
    xmaxSTA32ntmp(idis_ori) = mshift+1;
    
    nin_new = length(isave)
    for i = 1: nin_new
        STA1file_new((i-1)*winlen+1: i*winlen,:)= STA1file((isave(i)-1)*winlen+1 : ...
                                                            isave(i)*winlen, :);
        STA2file_new((i-1)*winlen+1: i*winlen,:)= STA2file((isave(i)-1)*winlen+1 : ...
                                                            isave(i)*winlen, :);
        STA3file_new((i-1)*winlen+1: i*winlen,:)= STA3file((isave(i)-1)*winlen+1 : ...
                                                            isave(i)*winlen, :);
        STAopt_new(1:nstanew,(i-1)*winlen+1: i*winlen)= STAopt(1:nstanew, (isave(i)-1)*winlen+1 : ...
                                                            isave(i)*winlen);
        
        STAamp_new(i,:) = STAamp(isave(i),:);
        
        STA12file_new(i,:) = STA12file(isave(i), :);
        STA13file_new(i,:) = STA13file(isave(i), :);
        STA32file_new(i,:) = STA32file(isave(i), :);
        
        istartkeep_new(i) = istartkeep(isave(i));
        aSTA12keep_new(i,:) = aSTA12keep(isave(i), :);
        aSTA13keep_new(i,:) = aSTA13keep(isave(i), :);
        aSTA32keep_new(i,:) = aSTA32keep(isave(i), :);
        loopoffkeep_new(i,:) = loopoffkeep(isave(i), :);
        
    end
    
    %%% Write results into files
    fid = fopen([rstpath,'/MAPS/tdecmapallup',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                 num2str(winlen/sps),'s',num2str(sps),'sps',num2str(nstanew), 'add'],'w+');
    if nstanew == 3
    fprintf(fid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %10.3f %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %.4f %.4f %.4f %d %d %d %.3f %.3f %.3f \n',...
                allrst_new(1:nin_new,:)');
%         fprintf(fid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %10.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %d %d %d %.4f %.4f %.4f %d %d %d %.3f %.3f %.3f %6.2f %6.2f %6.2f \n',allrst_new(1:nin_new,:)');
    elseif nstanew == 4
    fprintf(fid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %10.3f %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f \n',...
                allrst_new(1:nin_new,:)');
%         fprintf(fid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %10.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f %6.2f %6.2f %6.2f \n',allrst_new(1:nin_new,:)');
    end
    fclose(fid);
    
    fid = fopen(strcat(rstpath, '/MAPS/tdectraceallup_',IDENTIF,'_',num2str(lo),'-',num2str(hi),...
                '_',num2str(winlen/sps),'s',num2str(sps),'sps',num2str(nstanew), 'add'),'w+');
    tracefile = [STA1file_new(1:nin_new*winlen,:) STA2file_new(1:nin_new*winlen,2) ...
                 STA3file_new(1:nin_new*winlen,2) STAopt_new(:, 1:nin_new*winlen)'];
    if nstanew == 3
        fprintf(fid,'%.4f %.6f %.6f %.6f %.6f %.6f %.6f \n',tracefile');
    elseif nstanew == 4
        fprintf(fid,'%.4f %.6f %.6f %.6f %.6f %.6f %.6f %.6f \n',tracefile');
    end
    fclose(fid);
   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WAY 3 end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
  end
  
medlok=median(abs(loopoffkeep_new))     % median is actually effective to the 2nd column, which is loopoff
medaSTA12=median(aSTA12keep_new)        % 2nd column is aSTA12, which results from PARABOLIC FIT
medaSTA13=median(aSTA13keep_new)
medaSTA32=median(aSTA32keep_new)
  
end
  
  
  
  
  
  
  
  
  
  
  
  
