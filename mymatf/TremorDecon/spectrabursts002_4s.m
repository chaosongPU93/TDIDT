% spectrabursts002_4s.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The next script to run after 'seisbursts002_4s.m', to obtain the spectrum
% of each data window, in order to know what is the proper corner to filter 
% to bandpass the data to achieve a similar shape as the template. \
% Previously, we have been using 1.8-6.3 Hz for data, according to the 
% examination to the 9-s example in 'testiterdecon_exp.m'. To preserve the 
% most the high-freq energy in the template (as the stacking would lose
% some high-freq energy inevitably), we use a pretty high upper corner 
% passband for the template. And due to the noise at long-period in data, 1.8
% Hz is like an empirically proper lower corner. Therefore, our main goal
% is to get an ideal upper corner for the data, if it is not 6.3 hz.
%
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/06/15
% Last modified date:   2022/06/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
format short e   % Set the format to 5-digit floating point
clear
clc
close all

defval('normflag',0); %whether to normalize templates

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
  'KLNB '];     % determine the trio and order, here the 1st sta is PGC
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

%trial passband for templates
lowlet = 1.8;
hiwlet = 18;
%broader-band filtering passband for reading data
lo = 0.1;
% hi = 15;
hi = 18;
%trial passband for data wins
losig = 1.8;
hisig = 6.3;

%% prepare templates (Green's functions)
sps = 160;
templensec = 60;

ccstack = [];
for ista = 1: nsta
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
        num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao_catnew');
    ccstack(:,ista) = load(fname);
end
STA = detrend(ccstack);

% for ista=1:nsta
%     STA(:,ista)=Bandpass(STA(:,ista),sps,0.1,15,2,2,'butter');   % change 'bandpass' to 'Bandpass'
% end
%plot the raw templates, not filtered, not best aligned
% figure
% hold on
% plot(STA(:,1),'r')
% plot(STA(:,2),'b')
% plot(STA(:,3),'k')

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
% hold on
% plot(STAtmp(:,1),'r')
% plot(STAtmp(:,2),'b')
% plot(STAtmp(:,3),'k')
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
%   hiwlet=18;
%   lowlet=1.8;
  tmpwletf(:,ista) = Bandpass(tmpwlet(:,ista), sps, lowlet, hiwlet, 2, 2, 'butter');
  %detrend again for caution
  tmpwletf(:,ista)=detrend(tmpwletf(:,ista));
end

%%%constrained CC, so that only 2 offsets are independent
ccmid = round(size(tmpwletf,1)/2);
ccwlen = 10*sps;
loffmax = 5*sps/40;
ccmin = 0.01;  % depending on the length of trace, cc could be very low
iup = 1;    % times of upsampling
[off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(tmpwletf(:,1:3)',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
offwlet1i(1) = 0;
offwlet1i(2) = round(off12con);
offwlet1i(3) = round(off13con);

for ista = 4: nsta
  [mcoef,offwlet1i(ista)] = xcorrmax(tmpwletf(:,1), tmpwletf(:,ista), mshiftadd, 'coeff');
end

%%%automatically find the rough zero-crossing time, whose abs. value is closest to 0, whether + or -
[~,imin] = min(tmpwletf(:,1));
[~,imax] = max(tmpwletf(:,1));
[~,zcsta1] = min(abs(tmpwletf(imin:imax,1)));
zcsta1 = zcsta1+imin-1;
greenlen = pow2(9)*sps/40;
green = zeros(greenlen,nsta); % no bandpass
greenf = zeros(greenlen,nsta);  % bandpassed version
ppeaks = zeros(nsta,1);
npeaks = zeros(nsta,1); % negative peaks
for ista = 1: nsta
  %cut according to the zero-crossing and the time shift from the constrained CC
  green(:,ista) = tmpwlet(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenf(:,ista) = tmpwletf(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  %detrend again for caution
  green(:,ista)=detrend(green(:,ista));
  greenf(:,ista)=detrend(greenf(:,ista));
  if normflag
    %normalize by max amp
    green(:,ista) = green(:,ista)/max(abs(green(:,ista)));    % normalize
    greenf(:,ista) = greenf(:,ista)/max(abs(green(:,ista)));    % normalize
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
[off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(greenf(:,1:3)',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
if ~(off12con==0 && off13con==0)
  disp('Filtered templates are NOT best aligned');
end
for ista = 4: nsta
  [mcoef,mlag] = xcorrmax(greenf(:,1), greenf(:,ista), mshiftadd, 'coeff');
  if mlag~=0   % offset in samples
    fprintf('Filtered templates are NOT best aligned at %s \n',stas(ista,:));
  end
end

%%%plot the unfiltered and filtered templates
figure
subplot(3,1,1)
hold on
plot(green(:,1),'r')
plot(green(:,2),'b')
plot(green(:,3),'k')
text(0.95,0.9,'Raw','Units','normalized','HorizontalAlignment',...
  'right');
mx=max([abs(green(:,1)); abs(green(:,2))]);
xlim([0 greenlen])
ylim([-mx mx])
box on

subplot(3,1,2)
hold on
plot(greenf(:,1),'r')
plot(greenf(:,2),'b')
plot(greenf(:,3),'k')
text(0.95,0.9,sprintf('%.1f-%.1f Hz',lowlet,hiwlet),'Units','normalized','HorizontalAlignment',...
  'right');
mx=max([abs(greenf(:,1)); abs(greenf(:,2))]);

%%%running CC using a window length of 'cclen'
mwlen=sps/2;
[ircc,rcc12] = RunningCC(greenf(:,1), greenf(:,2), mwlen);
[~,rcc13] = RunningCC(greenf(:,1), greenf(:,3), mwlen);
[~,rcc23] = RunningCC(greenf(:,2), greenf(:,3), mwlen);
rcc = (rcc12+rcc13+rcc23)/3;
%alln(alln<0)=-10^4*yma; %just so they don't plot.
% plot(samples(cclen/2+1:greenlen-cclen/2),mx*alln,'co','markersize',2); % scale with data amp.
plot(ircc,mx*rcc,'co','markersize',2); % scale with data amp.
xlim([0 greenlen])
ylim([-mx mx])
box on
xc23=xcorr(greenf(:,2),greenf(:,3),10,'coeff');
xc13=xcorr(greenf(:,1),greenf(:,3),10,'coeff');
xc12=xcorr(greenf(:,1),greenf(:,2),10,'coeff');
[ccmax23,imax23]=max(xc23);
[ccmax13,imax13]=max(xc13);
[ccmax12,imax12]=max(xc12);

subplot(3,1,3)
%%% plot the spectrum, using 'pchave', which is capable of time-averaging a longer segment
wlen = greenlen;
polap = 75;
nfft = wlen;
Fs = sps;
%note that i request for only 1 window, so there is no robust estimate for the template, it is
%pretty much the same as 'periodogram', CHECK IT
for ista = 1: nsta
  [~,pcfw,~,~,~,~,pcwletbb(:,ista)]=pchave(green(:,ista),wlen,polap,nfft,Fs,'MAD','dpss');
  [~,~,~,~,~,~,pcwlet(:,ista)]=pchave(greenf(:,ista),wlen,polap,nfft,Fs,'MAD','dpss');
end
pcwletbb = sqrt(pcwletbb);
pcwlet = sqrt(pcwlet);
ax = gca;
loglog(ax,pcfw,pcwletbb(:,1),'r--','linew',1);
hold(ax,'on');
loglog(ax,pcfw,pcwletbb(:,2),'b--','linew',1);
loglog(ax,pcfw,pcwletbb(:,3),'k--','linew',1);

loglog(ax,pcfw,pcwlet(:,1),'color','r','linew',1.5);
loglog(ax,pcfw,pcwlet(:,2),'color','b','linew',1.5);
loglog(ax,pcfw,pcwlet(:,3),'color','k','linew',1.5);
xran = [0.1 20];
yran = [1e-4 1e1];
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,[0.1 1 10]);


%% broader-band spectra of data windows (spectral density)
% dates and ets
%%% NOTE: 'dates' is the dates that the tremor is active at the region of interest, so that the
%%% waveform at the stations on these dates you are looking at mainly comes from the sources at the
%%% region of interest. We want to see if there is a noticable change in spectra during the burst
%%% windows on these dates
dates = unique(trange(:,1));
years = unique(floor(dates/1000));
nets = length(years);

pcallbb = cell(nets, 1);  %spectral density
pcall = cell(nets, 1);  %spectral density

% keyboard
for iets = 1: nets
  % dates in each ets
  year = years(iets);
  datesets = dates(floor(dates/1000)==year);
  
  k = 0;  % counting the burst windows
  pcetsoptbb = [];
  pcetsopt = [];
  
  for i = 1: length(datesets)
    
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
  
    %read horizontal optimal and orthogonal components
    JDAY = num2zeropadstr(jday,3);
    MO=day2month(jday,year);     % EXTERNAL function, day2month, get the month of one particular date
    direc=[datapath, '/arch', yr,'/',MO,'/'];     % directory name
    prename=[direc,yr,'.',JDAY,'.00.00.00.0000.CN'];    %  path plus prefix of data file,
    disp(prename);
    [STAoptbb,~,fileflag] = rd_daily_bpdata(year,jday,prename,stas,PERMSTA,POLSTA,...
      PERMROTS,POLROTS,sps,lo,hi,npo,npa,[],[],[],[]);
    STAopt = rd_daily_bpdata(year,jday,prename,stas,PERMSTA,POLSTA,...
      PERMROTS,POLROTS,sps,losig,hisig,npo,npa,[],[],[],[]);

    if fileflag == 0    % means there are missing files
      fprintf('Day %s / %s will be omitted because of missing files. \n', yr, JDAY);
      continue    % continue to the next day
    end

    % get the windows (time ranges) of the same day, for efficiency
    rangetemp = trange(trange(:,1)==datesets(i), :);
    hfdayi = hfbnd(hfbnd(:,daycol)==datesets(i), :);  % inside bound of the day

    for j = 1: size(rangetemp,1)
      k = k+1;  
      disp(k);

      tst = rangetemp(j,2); % start and end time of bursts
      ted = rangetemp(j,3);
      
      tmaxi = hfdayi(:, seccol); % starting time of max power rate of half sec inside the ellipse
      tcnti = hfdayi(:, 15);  % the center of detecting win is the 15th col
      indtmaxi = find(tmaxi>=tst-0.1 & tmaxi<=ted+0.1);

%       %%%%Use a fixed range of time before and after the 0.5-s strongest arrival
%       tbuffer = 3;   % buffer time to include some coherent precursor or coda, first 1s will be tapered 
%       tstbuf = tst-tbuffer; % start and end time of bursts, buffer added
%       tedbuf = ted+tbuffer;
      %%%%Use the start and end of the 4-s detecting window
      tstbuf = min(tcnti(indtmaxi)-2);
      tedbuf = max(tcnti(indtmaxi)+2); 
      tlenbuf = tedbuf-tstbuf;
      wlen = greenlen;
      wlensec = wlen/sps;
      
      polap = 75;
      nfft = wlen;
      Fs = sps;
      optdatbb = STAoptbb(max(floor(tstbuf*sps),1): min(floor(tedbuf*sps-1),86400*sps), :);
      optdat = STAopt(max(floor(tstbuf*sps),1): min(floor(tedbuf*sps-1),86400*sps), :);
      
      pcoptbb = zeros(round(nfft/2)+1, nsta);
      pcopt = zeros(round(nfft/2)+1, nsta);
      %if the total data window length permits multiple overlapping windows, obtain a robust
      %estimate, however, if there is still 1 win available, non-robust estimate is returned
      if tlenbuf>=wlensec
        for ista = 1: nsta
          [pcoptbb(:,ista),pcft] = pchave(optdatbb(:,ista+1),wlen,polap,nfft,Fs,'MAD','dpss');
          pcopt(:,ista) = pchave(optdat(:,ista+1),wlen,polap,nfft,Fs,'MAD','dpss');
        end
      %if the total data window is too short, pad zeros, and there is no robust estimate for the 
      %template, use 'periodogram'
      else
        window = hann(size(optdatbb,1));
        for ista = 1: nsta
          [pcoptbb(:,ista),pcft] = periodogram(optdatbb(:,ista+1),window,nfft,Fs);
          pcopt(:,ista) = periodogram(optdat(:,ista+1),window,nfft,Fs);
%           [~,pcft,~,~,~,~,pcopt(:,ista)]=pchave(optdat(:,ista+1),wlen,polap,nfft,Fs,'MAD','dpss');
        end
      end  % if length is enough

      % stores for every rtm
      pcetsoptbb(:,:,k) = pcoptbb(:,:);   % ifreq, ista, irtm
      pcetsopt(:,:,k) = pcopt(:,:);   % ifreq, ista, irtm
      
    end % bursts in each date
    
  end % dates in each ets
  
  % store spectra of all windows of the different ETSs separately
  pcallbb{iets,1} = pcetsoptbb(:,:,:);
  pcall{iets,1} = pcetsopt(:,:,:);
end % all ets
  
% clean some temp variables to make more room
clear STAoptbb STAopt
clear optdatbb optdat 
clear pcoptbb pcopt 
clear pcetsoptbb pcetsopt  
  
%% plot the spectra of optimal components, spectral density
xran = [0.1 20];
yran = [1e-3 1e1];
[f] = plt_spectra_of_bursts(years,stas,pcft,pcallbb,xran,yran);
% orient(f.fig,'landscape');
fname = 'bb_spectra_of_bursts.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
% keyboard

yran = [5e-5 5e-1];
[f] = plt_spectra_of_bursts(years,stas,pcft,pcall,xran,yran);
% orient(f.fig,'landscape');
fname = 'bp_spectra_of_bursts.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

%% plot the spectra of normalized optimal components, spectral density
xran = [0.1 20];
yran = [5e-2 5e2];
minfnorm = 1.8; maxfnorm = 6.3;
[f] = plt_spectra_of_bursts_norm(years,stas,pcft,pcallbb,minfnorm,maxfnorm,xran,yran);  
% orient(f.fig,'landscape');
fname = 'bbnorm_spectra_of_bursts.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

yran = [1e-3 1e1];
[f] = plt_spectra_of_bursts_norm(years,stas,pcft,pcall,minfnorm,maxfnorm,xran,yran);  
% orient(f.fig,'landscape');
fname = 'bpnorm_spectra_of_bursts.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

keyboard
  
%% compare the spectral shape of templates and data windows
pcmedetsbb = zeros(round(nfft/2)+1, nsta, nets);
pcaveetsbb = zeros(round(nfft/2)+1, nsta, nets);
pcmedets = zeros(round(nfft/2)+1, nsta, nets);
pcaveets = zeros(round(nfft/2)+1, nsta, nets);
pcmedbb = zeros(round(nfft/2)+1, nsta);
pcavebb = zeros(round(nfft/2)+1, nsta);
pcmed = zeros(round(nfft/2)+1, nsta);
pcave = zeros(round(nfft/2)+1, nsta);
for ista = 1: nsta
  pcbbcomb = [];
  pccomb = [];
  for iets = 1: nets
    pcetsbb = pcallbb{iets,1};
    pcbb = squeeze(pcetsbb(:,ista,:));
    pcbb = sqrt(pcbb);  % convert power to amplitude
    pcbbcomb = [pcbbcomb pcbb];
    pcmedetsbb(:,ista,iets) = median(pcbb,2);  % median of all cols, ie, bursts
    pcaveetsbb(:,ista,iets) = mean(pcbb,2);  % mean of all cols, ie, bursts

    pcets = pcall{iets,1};
    pc = squeeze(pcets(:,ista,:));
    pc = sqrt(pc);  % convert power to amplitude
    pccomb = [pccomb pc];
    pcmedets(:,ista,iets) = median(pc,2);  % median of all cols, ie, bursts
    pcaveets(:,ista,iets) = mean(pc,2);  % mean of all cols, ie, bursts
    
  end
  pcmedbb(:,ista) = median(pcbbcomb,2);
  pcavebb(:,ista) = mean(pcbbcomb,2);
  pcmed(:,ista) = median(pccomb,2);
  pcave(:,ista) = mean(pccomb,2);
end

%%
f.fig = figure;
f.fig.Renderer = 'painters';
widin = 2.5*nsta;  % maximum width allowed is 8.5 inches
htin = 3;   % maximum height allowed is 11 inches
% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);

nrow = 1;
ncol = nsta;

for isub = 1:nrow*ncol
  f.ax(isub) = subplot(nrow,ncol,isub);
  f.ax(isub).Box = 'on';
  %     grid(f.ax(isub),'on');
end

pltxran = [0.06 0.98]; pltyran = [0.15 0.95];
pltxsep = 0.04; pltysep = 0.05; 
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

% color = ['r';'b';'k'];
for i = 1: nsta
  ax = f.ax(i);
  
  loglog(ax,pcfw,pcwletbb(:,i),'r--','linew',1);
  hold(ax,'on');
  loglog(ax,pcfw,pcwlet(:,i),'color','r','linew',1.5);
  
  loglog(ax,pcft,pcmedbb(:,i),'k--','linew',1); 
  loglog(ax,pcft,pcmed(:,i),'color','k','linew',1.5);
  
%   amprat = mean(pcmed(pcft>=losig & pcft<=hisig, i))./ ...
%     mean(pcwlet(pcfw>=losig & pcfw<=hisig, i));
%   npcmed = pcmed(:,i)/amprat;
%   loglog(ax,pcft,npcmed,'color',[.4 .4 .4],'linew',1.5);
  
  yran = [1e-4 1e0];
  xran = [0.1 20];
  % plot(ax,[1 1],yran,':','color',[0.7 0.7 0.7]);
  % plot(ax,[2 2],yran,':','color',[0.7 0.7 0.7]);
  plot(ax,[4 4],yran,':','color',[0.5 0.5 0.5]);
  % plot(ax,[8 8],yran,':','color',[0.7 0.7 0.7]);
  plot(ax,[4 8],[1e-1 2.5e-2],'-','linew',1.5,'color',[0.4 0.4 0.4]); %a line of slope -2
  plot(ax,[losig losig],yran,'-.','color',[0.5 0.5 0.5]);
  plot(ax,[hisig hisig],yran,'-.','color',[0.5 0.5 0.5]);
  plot(ax,[lowlet lowlet],yran,'-.','color',[1 0.65 0]);
  plot(ax,[hiwlet hiwlet],yran,'-.','color',[1 0.65 0]);
  xlim(ax,xran);
  ylim(ax,yran);  
  xticks(ax,[0.1 1 10]);
  text(ax,0.96,0.95,stas(i,:),'unit','normalized','HorizontalAlignment','right');
  
  xlabel(ax,'Frequency (Hz)','FontSize',10);
  longticks(ax,2);
end  
% ylabel(f.ax(1),'Spectral density (energy/Hz)','FontSize',10);
% ylabel(f.ax(1),sprintf('Spectral density (amp.^{2}/Hz)'),'FontSize',10);
ylabel(f.ax(1),sprintf('Spectral density (squared amplitude/Hz)'),'FontSize',10);
lgd = legend(f.ax(4),'BB template',sprintf('BP temp., %.1f-%.1f Hz',lowlet,hiwlet),...
  sprintf('Median of BB data'),sprintf('Med. of BP data, %.1f-%.1f Hz',losig,hisig),...
  'location','south','fontsize',7);
%make background transparent
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));


orient(f.fig,'landscape');
fname = 'spectra_summary.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));


  
  
  







