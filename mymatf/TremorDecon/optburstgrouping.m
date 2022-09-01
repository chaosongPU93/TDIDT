% optburstgrouping.m
%
% --This script tries to ask what is the appropriate threshold of the separation in time between
% consecutive detections when you group the strongest 0.5-s arrivals into tremor bursts which are
% then to be deconvolved. We have tried 35 s as the maximum, but it is rather arbitrary.
%
% --The importance of talking about it is because, if this T_sep is too big, then you definitely
% will have more longer burst windows which likely to contain more sources, both the 4-s detections
% and potentially sources that could be deconvolved. However, if there are too many sources with a
% fair scattering, it is likely that a single best alignment could be equally bad to every source.
%
% --But if you break this long window into several shorter ones in which sources occupy separatable
% regions, then you might have a better alignment for each shorter windows than for the entire
% window. Note that this alignment matters, as we compute the RCC based on the aligned waveforms.
% And if RCC could be higher in value on average than and the shape for each window is also different
% from that for the RCC computed based upon the entire window, then we might lose some real sources
% or obtain some unreal sources due to a bad RCC which is critical to the deconvolution.
%
% --So the idea is to try several different T_sep, and examine the resulted window length
% distribution; median RCC and overall CC computed upon the best alignment, averaged for all windows;
% fraction of time where there are 4-s detections, averaged for all windows. Hopefully there is an 
% optimal T_sep where the RCC reaches a maximum (even a local max) among all trial T_sep's.
%
% --But it looks like all these quantities are mototonically dereasing with the T_sep, which seems
% not helpful to determine what is the optimal T_sep for grouping the 4-s detections into tremor
% burst windows.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/05/19
% Last modified date:   2022/05/19
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

scrsz=get(0,'ScreenSize');

workpath = getenv('ALLAN');
%%% choose the data path here 
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/PGCtrio');

freqflag='hf';  % flag to indicate whether to do hf or lf;

FLAG = 'PGC'; % detector
  
fam = '002';   % family number

sps = 40;

iup = 4;  % upsample 4 times

cutout = 'ellipse';

%load detections inside the cutout boundary of interest, an output of 'locinterp002_4s.m'
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
PREFIX = strcat(fam,'.up.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),'.',...
                        int2str(npo),int2str(npa),'.ms', int2str(mshift));
hfbnd = load(strcat(rstpath, '/MAPS/eloc.pgc.tdec_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',...
            num2str(winlen/sps),'s',num2str(sps),'sps4add_',cutout(1:4)));
          
%% compute relative time in days
%sort according to day, sec of the strongest arrival of the window
daycol = 14;
seccol = 16;
hfbnd = sortrows(hfbnd, [daycol, seccol]);
ncol = size(hfbnd,2);

%obtain the relative time in days
i03 = find(hfbnd(:,daycol) < 2004*1000);
hfbnd(i03,ncol+1) = (hfbnd(i03,daycol)-2003060)+hfbnd(i03,seccol)./(3600.*24);
i04 = find(hfbnd(:,daycol) < 2005*1000 & hfbnd(:,daycol) > 2004*1000);
hfbnd(i04,ncol+1) = (hfbnd(i04,daycol)-2004194)+hfbnd(i04,seccol)./(3600.*24);
i05 = find(hfbnd(:,daycol) > 2005*1000);
hfbnd(i05,ncol+1) = (hfbnd(i05,daycol)-2005254)+hfbnd(i05,seccol)./(3600.*24);
          
%% separation in time between itself and its preceding detection
%for the detections, obtain the separation in time between itself and its preceding
%detection for the catalog of each ETS.
%1st col: occurence time
%2nd col: inter-detection time to its preceding detection
hfinter03 = interevt_time(hfbnd(i03,end));
hfinter04 = interevt_time(hfbnd(i04,end));
hfinter05 = interevt_time(hfbnd(i05,end));
hfinter = [hfinter03; hfinter04; hfinter05];
                    
%% group the tremor bursts 
% ttol = 1e-3*ones(3,1); 
% ntol = 10;
ttoltrial = 15:5:45;
for itry = 1: length(ttoltrial)
ttol = ttoltrial(itry)/86400*ones(3,1);  % this is about 4*median
ntol = 3;
[bursthf03, n03] = group_tremor_burst_indep(hfinter03,ttol(1),ntol);
[bursthf04, n04] = group_tremor_burst_indep(hfinter04,ttol(2),ntol);
[bursthf05, n05] = group_tremor_burst_indep(hfinter05,ttol(3),ntol);

% recover the index of the detections to occurence times with the same format as trange
if ~isempty(bursthf03)
  [thf03,tranhf03,ntot03] = burst_range(bursthf03,hfinter03,2003060);
  perchf03 = ntot03/length(i03)*100;
else
  thf03 = []; tranhf03 = []; ntot03 = 0; perchf03 = 0;
end
if ~isempty(bursthf04)
  [thf04,tranhf04,ntot04] = burst_range(bursthf04,hfinter04,2004194);
  perchf04 = ntot04/length(i04)*100;
else
  thf04 = []; tranhf04 = []; ntot04 = 0; perchf04 = 0; 
end
if ~isempty(bursthf05)
  [thf05,tranhf05,ntot05] = burst_range(bursthf05,hfinter05,2005254);
  perchf05 = ntot05/length(i05)*100;
else
  thf05 = []; tranhf05 = []; ntot05 = 0; perchf05 = 0;
end

%% obtain the separation time between each burst and the length in time of each burst
tsep03 = zeros(size(thf03,1),1);   %what is the separation in time between each window?
tlen03 = (thf03(:,2)-thf03(:,1))*86400;  %what is the length of each window?
for i = 2: size(thf03,1)
  tsep03(i) = (thf03(i,1)-thf03(i-1,2))*86400;
end
tsep04 = zeros(size(thf04,1),1);
tlen04 = (thf04(:,2)-thf04(:,1))*86400;
for i = 2: size(thf04,1)
  tsep04(i) = (thf04(i,1)-thf04(i-1,2))*86400;
end
tsep05 = zeros(size(thf05,1),1);
tlen05 = (thf05(:,2)-thf05(:,1))*86400;
for i = 2: size(thf05,1)
  tsep05(i) = (thf05(i,1)-thf05(i-1,2))*86400;
end
tsep = [tsep03; tsep04; tsep05];
tlen = [tlen03; tlen04; tlen05];

% figure
% subplot(221)
% scatter(tlen03,log10(tsep03),15,'r','filled'); hold on
% xlabel('Length in time (s) of each burst ');
% ylabel('log_{10}(Separation in time (s) between bursts)');
% scatter(tlen04,log10(tsep04),15,'b','filled');
% scatter(tlen05,log10(tsep05),15,'k','filled');
% ax=gca; plot(ax.XLim,[log10(ttol(1)*86400) log10(ttol(1)*86400)], 'k--');
% box on; grid on;
% legend('2003','2004','2005','T_{max}','Location','northeast');
% 
% subplot(222)
% scatter(tlen,log10(tsep),15,[.4 .4 .4],'filled'); hold on
% xlabel('Length in time (s) of each burst ');
% ylabel('log_{10}(Separation in time (s) between bursts)');
% ax=gca; plot(ax.XLim,[log10(ttol(1)*86400) log10(ttol(1)*86400)], 'k--');
% box on; grid on;
% legend('All','T_{max}','Location','northeast');
% 
% subplot(223)
% scatter(tlen03,tsep03,15,'r','filled'); hold on
% xlabel('Length in time (s) of each burst ');
% ylabel('Separation in time (s) between bursts');
% scatter(tlen04,tsep04,15,'b','filled');
% scatter(tlen05,tsep05,15,'k','filled');
% axis([-10 400 0 100]);
% ax=gca; plot(ax.XLim,[ttol(1)*86400 ttol(1)*86400], 'k--');
% box on; grid on;
% legend('2003','2004','2005','T_{max}','Location','northeast');
% 
% subplot(224)
% scatter(tlen,tsep,15,[.4 .4 .4],'filled'); hold on
% xlabel('Length in time (s) of each burst ');
% ylabel('Separation in time (s) between bursts');
% axis([-10 400 0 100]);
% ax=gca; plot(ax.XLim,[ttol(1)*86400 ttol(1)*86400], 'k--');
% box on; grid on;
% legend('All','T_{max}','Location','northeast');
% text(0.5,0.4,sprintf('<=40 s: %d/%d',sum(tsep<=40)-3,size(tsep,1)),'Units','normalized');
% text(0.5,0.5,sprintf('<=50 s: %d/%d',sum(tsep<=50)-3,size(tsep,1)),'Units','normalized');
% text(0.5,0.6,sprintf('<=60 s: %d/%d',sum(tsep<=60)-3,size(tsep,1)),'Units','normalized');
% text(0.5,0.7,sprintf('<=80 s: %d/%d',sum(tsep<=80)-3,size(tsep,1)),'Units','normalized');


%% some simple statistics of the bursts
%times of all tremor bursts
trange = [tranhf03;tranhf04;tranhf05];

tlens{itry} = tlen;
medtlen(itry) = median(tlen);

% keyboard

%% output results to file  
fid = fopen(strcat(rstpath, '/MAPS/tdec.bstran',num2str(round(ttol(1)*86400)),'s.pgc002.',cutout(1:4)),'w+');
fprintf(fid,'%d %9.2f %9.2f \n',trange');
fclose(fid);

end

%%
figure
for itry = 1: length(ttoltrial)
  subplot(2,4,itry)
  tlen = tlens{itry};
  histogram(tlen); hold on;
  text(0.1,0.9,sprintf('med: %.2f s',medtlen(itry)),'Units','normalized');
  text(0.9,0.9,sprintf('num: %d',length(tlen)),'Units','normalized','HorizontalAlignment','right');
  text(0.9,0.8,sprintf('T_{sep}: %d s',ttoltrial(itry)),'Units','normalized','HorizontalAlignment','right');
  xlabel('Window length (s)');
  ylabel('Counts');
  
end


%%
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

for itry = 1: length(ttoltrial)
ttol = ttoltrial(itry);
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.pgc002.',cutout(1:4)));
tlen = trange(:,3)-trange(:,2);


%%% read daily data, break into windows of segments, plot
% dates and ets
%%% NOTE: 'dates' is the dates that the tremor is active at the region of interest, so that the
%%% waveform at the stations on these dates you are looking at mainly comes from the sources at the
%%% region of interest. We want to see if there is a noticable change in spectra during the burst
%%% windows on these dates
dates = unique(trange(:,1));
years = unique(floor(dates/1000));
nets = length(years);

%filtering passband for reading data
hisig=6.3; % this will give a similar spectral shape between template and signal
losig=1.8;

%moving window length in samples for running CC, envelope, etc.
mwlen=sps/2;
% mwlen=sps;

k = 0;  % counting the burst windows

off1i = zeros(size(trange,1),3);  % the best alignment between sta2, sta3 wrt sta1
ccali = zeros(size(trange,1),1);  % CC using the best alignment

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
    
    %read horizontal optimal and orthogonal components
    [STAopt,STAort,~,fileflag] = rd_daily_bpdata(year,jday,datapath,stas,PERMSTA,POLSTA,...
      PERMROTS,POLROTS,sps,losig,hisig,npo,npa,[],[],[],[]);
    
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
      tbuffer = 3;   % buffer time to include some coherent precursor or coda, first 1s will be tapered 
      tstbuf = tst-tbuffer; % start and end time of bursts, buffer added
      tedbuf = ted+tbuffer;
      tlenbuf = tedbuf-tstbuf;
            
      %align records
      msftadd = (round(max(abs([off12ran off13ran])))+1)*sps/40;  %+1 for safety
%       msftadd = round(sps/8);    % maximum allowed shift between 2 traces
      optcc = STAopt(max(floor((tstbuf-1)*sps+1),1): min(floor((tedbuf+1)*sps),86400*sps), 2:end);
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
      optdat = [];  % win +/-3 s, segment of interest,first 1s will be tapered 
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
      
      %how many 4-s detections fall into the burst range 
      indtmaxi = find(tmaxi>=tst-0.1 & tmaxi<=ted+0.1);
      ninbst(k,1) = length(indtmaxi); %+-0.1 s in case of resolution while saving to file
      
      %%%taper the signal and obtain the new rcc between tapered signals
      sigsta = zeros(size(optdat,1), nsta);
      for ista = 1:nsta
        tmp = optdat(:,ista+1); %best aligned, filtered
        %detrend and taper only the data, NOT the noise
        tmp = detrend(tmp);
        ltmp = length(tmp);
        fractap = 2*sps/ltmp; % if fractap is >=1, n-point von Hann window is returned
        ptstap = fractap/2*size(tmp,1); % if fractap is >=1, n-point von Hann window is returned
        w = tukeywin(size(tmp,1),fractap);
        tmp = w.* tmp;
        tmp = detrend(tmp); %detrend again for caution
        sigsta(:,ista) = tmp;
      end
      %compute running CC between 3 stations
      [ircc,rcc12] = RunningCC(sigsta(:,1), sigsta(:,2), mwlen);
      [~,rcc13] = RunningCC(sigsta(:,1), sigsta(:,3), mwlen);
      [~,rcc23] = RunningCC(sigsta(:,2), sigsta(:,3), mwlen);
      rcc = (rcc12+rcc13+rcc23)/3;
      mrcc(k,1) = median(rcc);
      
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
            
      %overall average max CC based on the current best alignment
      mcc12 = sum(sigsta(:,1).*sigsta(:,2))./ ...
        sum(sqrt(sigsta(:,1).^2).*sqrt(sigsta(:,2).^2));
      mcc13 = sum(sigsta(:,1).*sigsta(:,3))./ ...
        sum(sqrt(sigsta(:,1).^2).*sqrt(sigsta(:,3).^2));
      mcc23 = sum(sigsta(:,2).*sigsta(:,3))./ ...
        sum(sqrt(sigsta(:,2).^2).*sqrt(sigsta(:,3).^2));
      mcc(k,1) = (mcc12+mcc13+mcc23)/3;

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
      indnodet = setdiff(setdiff(1:length(rcc), inddeti), inddeto);
      
    end
    
  end
    
end

mccs{itry} = mcc;
mrccs{itry} = mrcc;
tfracis{itry} = tfraci;

%average over all bursts windows
mmccs(itry) = mean(mcc);
mmrccs(itry) = mean(mrcc);
mtfracis(itry) = mean(tfraci);

end


%% A few more plots

fttpfree = fittype( @(a,b,x) a*x+b);

%%%Median RCC of the data win vs. fraction of time that has 4-s in-bound detections
figure
for itry = 1: length(ttoltrial)
  ax = subplot(2,4,itry);
  mrcc = mrccs{itry};
  tfraci = tfracis{itry};
  scatter(ax,tfraci,mrcc,20,1:length(tfraci),'filled'); hold(ax,'on');
  colorbar; colormap jet; box on; grid on;
  [speartmp,~] = corr(tfraci(:),mrcc,'Type','Spearman'); % Spearman coeff
  text(ax,0.95,0.9,sprintf('S: %.2f',speartmp),'unit','normalized',...
    'HorizontalAlignment','right','fontsize',10);
  [fitobj,~,~] = fit(tfraci(:),mrcc,fttpfree,'Robust','Bisquare',...
    'StartPoint',[1 1]);
  coef = coeffvalues(fitobj);
  coefci = confint(fitobj);
  xplt = min(tfraci(:)): 0.01: max(tfraci(:));
  yfit = coef(1)*xplt + coef(2);
  ylow = coefci(1,1)*xplt + coefci(1,2);
  yupp = coefci(2,1)*xplt + coefci(2,2);
  plot(ax,xplt,yfit,'k-','linew',2);
  patch(ax,[xplt fliplr(xplt)],[ylow fliplr(yupp)],'k','Facealpha',0.2,'edgecolor','none');
  text(0.9,0.8,sprintf('T_{sep}: %d s',ttoltrial(itry)),'Units','normalized',...
    'HorizontalAlignment','right');
  xlabel(ax,'Fraction of time with 4-s in-bound detections');
  ylabel(ax,'Median RCC of selected window');
  hold(ax,'off');
end

%%%Median RCC or overall CC of the data win vs. fraction of time that has 4-s in-bound detections
figure
for itry = 1: length(ttoltrial)
  ax = subplot(2,4,itry);
  mcc = mccs{itry};
  tfraci = tfracis{itry};
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
  text(0.9,0.8,sprintf('T_{sep}: %d s',ttoltrial(itry)),'Units','normalized',...
    'HorizontalAlignment','right');
  xlabel(ax,'Fraction of time with 4-s in-bound detections');
  ylabel(ax,'Overall CC selected window');
  hold(ax,'off');
end

%%
figure
plot(mmccs,'ro-','MarkerSize',8); hold on
plot(mmrccs,'bs-','MarkerSize',8);
plot(mtfracis,'kd-','MarkerSize',8);
legend('Average overall CC','Average median RCC','Average Fraction of time','Location','northeast');
grid on; box on;








