% analyze_synth.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --This script aims to do some synthetics to the synthetic seismograms other
% than deconvolution that is done particularly in 'decon_synth.m', to be 
% compared to real data. Analysis would include, tracking RCC variation with
% different saturation level and region size; tracking RCC vs. envelope (?)
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/05/26
% Last modified date:   2023/05/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e   % Set the format to 5-digit floating point
clear
clc
close all

%% for easy testing
defval('normflag',0); %whether to normalize templates
defval('rccmwsec',0.5); %moving win len in sec for computing RCC

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
% if pltflag
%     set(0,'DefaultFigureVisible','on');
% else
%     set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
% end

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

%% prepare templates (Green's functions), from 'lfetemp002_160sps.m' or Allan's templates
% tempflag = 'allan';
tempflag = 'chao';

adatapath = '/home/data2/chaosong/matlab/allan/matfils/';  %path for Allan's data

if strcmp(tempflag,'chao')
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
[off12con,off13con,cc] = constrained_cc_interp(greenf(:,1:3)',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
if ~(off12con==0 && off13con==0)
  disp('Filtered templates are NOT best aligned \n');
end
amprat(1,:) = minmax(greenf(:,1)')./minmax(greenf(:,2)');	% amp ratio between max at sta 3 and 2 or min
amprat(2,:) = minmax(greenf(:,1)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min  
amprat(3,:) = minmax(greenf(:,2)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min  
spread = range(greenf);   % range of the amp of template

%%%plot the unfiltered and filtered templates
plt_templates(green,greenf,stas,[],[],lowlet,hiwlet,sps);

%just the filtered templates
% plt_templates_bp(greenf,stas,lowlet,hiwlet,sps);

%% load synthetic seismograms
%%%specify distribution for source location  
distr='UN';  % uniform distribution

%%%specify distribution for source location  
% distrloc = 'custompdf'; %using a custom PDF function
distrloc = 'uniform'; %uniformly random in a specified region,

fracelsew=0; %0.25 %0.5; %0.6; %The ratio of elsewhere events to local events.  Zero for Discrete Ide.

%%%specify if considering the physical size of each source
physicalsize = 1;
% physicalsize = 0;

diam=0.3;% 0.5; %0.6; %

%%%specify shape of the source region
srcregion='ellipse';
% srcregion='rectangle';
% srcregion='circle';

%params of limited source region, subject to variation!
if strcmp(srcregion,'circle')
  shiftor=[0.2 0.2]; %(in km) %center of the same ellipse of my 4-s catalog
  radi=1.25; %radius
  [xcut,ycut] = circle_chao(shiftor(1),shiftor(2),radi,0.01);
elseif strcmp(srcregion,'ellipse')
  xaxis=1.75; %axis length of the same ellipse of my 4-s catalog
  yaxis=1.25;
  shiftor=[0.2 0.2]; %(in km) %center of the same ellipse of my 4-s catalog
  [xcut,ycut] = ellipse_chao(shiftor(1),shiftor(2),xaxis,yaxis,0.01,45,shiftor);
end

%%%specify regime for transformation from time offset to map location 
% ftrans = 'interpArmb';
% ftrans = 'interpArmbreloc';
ftrans = 'interpchao';

%%%specify if forcing a min speration of arrival time for 2 events from the same spot 
forcesep = 1;
% forcesep = 0;

%%%flag for validating if ground truth of sources can recover the record
% testsrcflag = 1;
testsrcflag = 0;

%%%flag for validing if the spectral shapes of data and templates are similar
% testfreqflag = 1;
testfreqflag = 0;

%times of saturation 
nsat=[0.1 0.4 1 2 4 10 100];  
% nsat=[0.4 2 10];
nnsat = length(nsat);

%%%file name prefix of synthetics
if strcmp(srcregion,'ellipse')
  fname = ['/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),'_',...
    num2str(xaxis),'-',num2str(yaxis),'.diam',num2str(diam),'.else',num2str(fracelsew,2),...
    'nsat'];
elseif strcmp(srcregion,'circle')
  fname = ['/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
    '_',num2str(radi),'.diam',num2str(diam),'.else',num2str(fracelsew,2),...
    'nsat'];
end
% fname = '/synthetics/STAS.UN.160sps.ell_2-1.25.diam0.3.else0nsat';


% for insat = 1: nnsat  
insat = 1;
disp(nsat(insat));

%%%load synthetics of certain saturation level
STAopt = load(strcat(workpath,fname,num2str(nsat(insat))));

%%%load sources
synsrc = load(strcat(workpath,fname,num2str(nsat(insat)),'_sources'));

if strcmp(distrloc, 'uniform') 
  xygrid = load([workpath,'/synthetics/synsrcloc.',srcregion(1:3),'.grd']);  
  tmp = xygrid(synsrc(:,2),:);
  synsrc = [synsrc(:,1) tmp(:,1:4) ones(length(tmp),1)];  %[indtarvl, off12, off13, loce, locn, amp]
elseif strcmp(distrloc, 'custompdf') 
  [loc, indinput] = off2space002(synsrc(:,4:5),sps,ftrans,0);
  synsrc = [synsrc(:,1) loc(:,1:4) ones(length(loc),1)];  %[indtarvl, off12, off13, loce, locn, amp]
end
  
%%%load starting indices of added sources at sta 1 
synsrcstind = load(strcat(workpath,fname,num2str(nsat(insat)),'_stind'));

skiplen=greenlen;

%% filter data (and templates)
%%filter data
hisig=6.3; % this will give a similar spectral shape between template and signal
losig=1.8;
optseg = [];
for ista = 1:nsta
  optseg(:,ista) = Bandpass(STAopt(:,ista), sps, losig, hisig, 2, 2, 'butter');
end

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
bufsec = 1;
msftaddm = bufsec*sps;  %buffer range for later CC alignment
rccmwsec = 0.5; 
rccmwlen = rccmwsec*sps;  %window length for computing RCC
overshoot = rccmwlen/2; %RCC points to the center of computing window, so extra room is needed

%%%obtain a single best alignment based on the entire win 
optcc = detrend(optseg(1+msftaddm: end-msftaddm, :));
msftadd=10*sps/40;
loffmax = 4*sps/40;
ccmid = ceil(size(optcc,1)/2);
ccwlen = round(size(optcc,1)-2*(msftadd+1));  % minus ensures successful shifting of records
ccmin = 0.01;  % depending on the length of trace, cc could be very low
iup = 1;    % times of upsampling

[off12con,off13con,ccali,iloopoff,loopoff] = constrained_cc_interp(optcc',ccmid,...
  ccwlen,msftadd,loffmax,ccmin,iup);
% if a better alignment cannot be achieved, use 0,0
if off12con == msftadd+1 && off13con == msftadd+1
  off12con = 0;
  off13con = 0;
  fprintf('Current window cannot be properly aligned, double-check needed \n');
end
off1i = zeros(nsta,1);
off1i(2) = round(off12con);
off1i(3) = round(off13con);

%%%Align and compute the RCC based on the entire win, and take that as the input signal!
optdat = [];  % win segment of interest
for ista = 1: nsta
  optdat(:, ista) = optseg(1+msftaddm-off1i(ista): end-msftaddm-off1i(ista), ista);
end

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
rcc = (rcc12+rcc13+rcc23)/3;
rccpair = [rcc12 rcc13 rcc23];
rcc1i = zeros(length(rcc),nsta-3);

sigsta = detrend(sigsta(overshoot+1:end-overshoot, :));  %excluding the overshoot

cc12 = xcorr(sigsta(:,1), sigsta(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
cc13 = xcorr(sigsta(:,1), sigsta(:,3),0,'normalized');
cc23 = xcorr(sigsta(:,2), sigsta(:,3),0,'normalized');
ccpair = [cc12 cc13 cc23];
mrcc = median(rcc);
mcc = (cc12+cc13+cc23)/3;

%if only use the mean RCC from pair 12 and 13
rcc = mean(rccpair(:,[1 2]), 2);
rccsat(:,insat) = rcc;

%a qucik look of the data
figure
lsig = size(sigsta,1);
plot((1:lsig)/sps,sigsta(:,1),'r'); hold on
plot((1:lsig)/sps,sigsta(:,2),'b');
plot((1:lsig)/sps,sigsta(:,3),'k');
ax = gca;
axsym(ax);
plot((1:lsig)/sps,rcc*ax.YLim(2),'o','color',[.6 .6 .6],'markersize',2);
text(0.95,0.1,sprintf('Saturation: %.1f',nsat(insat)),'Units','normalized','HorizontalAlignment','right');
xlim([0 20]);

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
  else
    greenst = synsrcstind-source(:, ista); %tarvl1-off is the arrival time in index at sta 2 and 3
  end
  
  %%%you don't need all impulses, only some of them contribute to the length of truncated record
  %%%cut out the green indice and sources that contribute
  %%%Note that some zero-crossings might appear later than 'indedsig', but the corresponding start
  %%%index of the green is actually smaller than 'indedsig', ie. their 'greenst' <= 'indedsig'
  induse = find(greenst>=skiplen-greenlen+indst & greenst<=inded+skiplen-1);
  greenst = greenst(induse);
  source = source(induse,:);

  greenzc = greenst+tgreen; % index of approximate zero-crossing
  source(:, 1) = source(:, 1)+tgreen;  % now 'greenzc' should be the same as 1st col of 'source'
    
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
  impgt(:,3) = impgt(:,1)-synsrcgt(:,2)-off1i(2); % note the sign is consistent!
  impgt(:,5) = impgt(:,1)-synsrcgt(:,3)-off1i(3);
  impgt(:,7) = synsrcgt(:,2);  % no need to shift offset, because they are 'true' offset
  impgt(:,8) = synsrcgt(:,3);
  impgt(:,9) = impgt(:,8)-impgt(:,7);  % off23, == off13 - off12 == tarvl2 - tarvl3
  
% end

%%%plot the ground truth source offset
spsscale = sps/40;
loff_max = 6*spsscale;  % maximum allowed shift is 4 samples at 40 Hz, it is also true from the synthetic distribution
xran = [-loff_max+off1i(2)-1 loff_max+off1i(2)+1];
yran = [-loff_max+off1i(3)-1 loff_max+off1i(3)+1];
offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
cran = [0 lsig];
f1.fig = figure;
f1.fig.Renderer = 'painters';
ax1=gca;
[ax1,torispl,mamp] = plt_decon_imp_scatter(ax1,impgt,xran,yran,cran,offxran,offyran,...
  sps,30,'mean','tarvl');
scatter(ax1,off1i(2),off1i(3),20,'ks','filled','MarkerEdgeColor','k');
title(ax1,'Ground truth');

%%%plot the transformed ground truth source map locations
xran = [-4 4];
yran = [-4 4];
cran = [0 lsig/sps];
f2.fig = figure;
f2.fig.Renderer = 'painters';
ax2=gca;
[ax2] = plt_decon_imp_scatter_space(ax2,impgt,xran,yran,cran,offxran,...
  offyran,sps,30,ftrans,'mean','tarvl');
plot(ax2,xcut,ycut,'k-','linew',2);
title(ax2,'Ground truth');

%%%plot the cumulative density of detections
[impgtloc, ~] = off2space002(impgt(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[f] = plt_cumulative_density(impgtloc,[],xran,yran,'pixel',10,10);
title(f.ax(1),'Ground truth');

%%%plot the ground truth source map locations from grid
f.fig = figure;
f.fig.Renderer = 'painters';
ax=gca;
hold(ax,'on');
ax.Box='on'; grid(ax,'on');
plot(ax,[xran(1) xran(2)],[0 0],'--','Color',[.3 .3 .3]); hold(ax,'on')
plot(ax,[0 0],[yran(1) yran(2)],'--','Color',[.3 .3 .3]);
scatter(ax,synsrcgt(:,4),synsrcgt(:,5),synsrcgt(:,6)*30,impgt(:,1)/sps,'filled','MarkerEdgeColor',[.5 .5 .5]);
scatter(ax,loc0(:,1),loc0(:,2),30,'k','filled');
scatter(ax,shiftor(:,1),shiftor(:,2),30,'k','linew',1);
plot(ax,xcut,ycut,'k-','linew',2);
oldc = colormap(ax,'kelicol');
newc = flipud(oldc);
colormap(ax,newc);
c=colorbar(ax);
caxis(ax, cran);
c.Label.String = sprintf('Arrival time at PGC (s)');
text(ax,0.98,0.05,sprintf('%d events',size(synsrcgt,1)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
ax.YAxis.FontSize = 8;
ax.XAxis.FontSize = 8;
xlabel(ax,'E (km)','FontSize',11);
ylabel(ax,'N (km)','FontSize',11);
axis(ax, 'equal');
xlim(ax,xran); xticks(ax,xran(1): 1 : xran(2));
ylim(ax,yran); yticks(ax,yran(1): 1 : yran(2));







