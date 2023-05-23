% decon_synth.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --This script is similar to but go beyond 'decon_allansyn.m'. It also
% absorbs the part of code of deconvolution previously in 'synthshift_chao.m'
% due to the purpose of separating te functionality into two main parts:
% synthetics generation and deconvolution.
% --This is the script to read and deconvolve the synthetic seismograms from
% either Allan's or my codes on using region size and temporal saturation level.
% Other than synthetic seismograms, we also know the exact information of
% the elemental LFEs, including their arrival time index (or even origin time
% depending on which algorithm used) and location (in terms of offsets 12 
% and 13, and may be map locations as well). These information will be 
% useful when comparing with the deconvolution result.
% --The deconvolution is implemented by the same algorithm in 
% 'deconv_4s_exp_4thsta_fn.m' by regarding the full window as a whole. Some 
% prealignment before decon may be necessary to optimize the RCC.
% The goal is to get the location, and distance between event 
% pairs (N and N-m, or each source to all others within some time 
% separation).
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/05/09
% Last modified date:   2023/05/09
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
% srcregion='ellipse';
% srcregion='rectangle';
srcregion='circle';

%params of limited source region
if strcmp(srcregion,'ellipse')
  xaxis=1.75; %axis length of the same ellipse of my 4-s catalog
  yaxis=1.25;
  shiftor=[0.2 0.2]; %(in km) %center of the same ellipse of my 4-s catalog
  [xcut,ycut] = ellipse_chao(shiftor(1),shiftor(2),xaxis,yaxis,0.01,45,shiftor);
elseif strcmp(srcregion,'circle')
  shiftor=[0.2 0.2]; %(in km) %center of the same ellipse of my 4-s catalog
  radi=1.25; %radius
  [xcut,ycut] = circle_chao(shiftor(1),shiftor(2),radi,0.01);
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
% nsat=[0.1 0.4 1 2 4 10 100];  
nsat=[0.4 2 10];
nnsat = length(nsat);

insat = 3;
disp(nsat(insat));

%%%load synthetics
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

% STAopt = cell(nnsat,1);
% for i = 1: nnsat
%   STAopt{i} = load(strcat(workpath,fname,num2str(nsat(i))));
% end
STAopt = load(strcat(workpath,fname,num2str(nsat(insat))));

%%%load sources
% synsrc = cell(nnsat,1);
% for i = 1: nnsat
%   synsrc{i} = load(strcat(workpath,fname,num2str(nsat(i)),'_sources'));
% end
synsrc = load(strcat(workpath,fname,num2str(nsat(insat)),'_sources'));

if strcmp(distrloc, 'uniform') 
  xygrid = load([workpath,'/synthetics/synsrcloc.grd']);  
  tmp = xygrid(synsrc(:,2),:);
  synsrc = [synsrc(:,1) tmp(:,1:4) ones(length(tmp),1)];  %[indtarvl, off12, off13, loce, locn, amp]
elseif strcmp(distrloc, 'custompdf') 
  [loc, indinput] = off2space002(synsrc(:,4:5),sps,ftrans,0);
  synsrc = [synsrc(:,1) loc(:,1:4) ones(length(loc),1)];  %[indtarvl, off12, off13, loce, locn, amp]
end
  
%%%load starting indices of added sources at sta 1 
synsrcstind = load(strcat(workpath,fname,num2str(nsat(insat)),'_stind'));

rng('default');
Twin=0.5*3600+3+2*ceil(greenlen/sps); %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
winlen=Twin*sps+1;
skiplen=greenlen;
%%%Here the noise is 'uniform' in time!
% noistd = 5e-2;
noistd = 2.e-7;
synth=noistd*(randn(winlen+greenlen+2*10,nsta)-0.5); %+2*10 a little extra, for jiggering 2nd & 3rd stations.% for ista=1:nsta

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
    else
      %ind - rnoff is the arrival time in index at sta 2 and 3
      %note the sign here, if off12 >0, move 2 to the right to align with 1, meaning that 2 is
      %earlier than 1, ie., tarvl2 < tarvl1. Be consistent all the time
      greenst = synsrcstind-source(:, ista);
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
    impamp = zeros(max(greenzc)+20,1);
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
    title('Reproduce synthetics with truncated convolution');
    subplot(412)
    imptemp = find(impamp>0);
    p1=stem(imptemp-skiplen, impamp(imptemp),'b','MarkerSize',4); hold on;
    % p1=stem((1:size(impamp))-skiplen, impamp,'b','MarkerSize',4); hold on;
    ax = gca;
    plot([indst indst],ax.YLim,'--','color',[.5 .5 .5]);
    plot([inded inded],ax.YLim,'--','color',[.5 .5 .5]);
    legend(p1,'Synthetic random impulses');
    subplot(413);
    plot(sigstagt(:,ista),'b'); hold on
    plot(sigconvstagt(:,ista),'k');
    legend('Truncated synthetic signal','Truncated signal from convolution');
    ax=gca; yran = ax.YLim;
    subplot(414);
    plot(sigstagt(:,ista)-sigconvstagt(:,ista),'k'); hold on
    legend('Difference'); ylim(yran)
    
    %now 'source' contain impulse triplets that contribute to the trace segment of interest, but we
    %want to get ones whose zero-crossing falls into the window.
    source(:,1) = source(:,1)-skiplen;  % index after skipping
    tmp = source(source(:,1)>=indst & source(:,1)<=inded, :);
    tmp = sortrows(tmp,1,'ascend');
    synsrcgtsta{ista} = tmp; %ideally for each station this should be the same, but coincidence is possible
  end
  %notify if sources are the same at diff stations
  if ~isequaln(synsrcgtsta{1,1},synsrcgtsta{2,1}) || ~isequaln(synsrcgtsta{1,1},synsrcgtsta{3,1})
    disp('Extracted sources at different stations are not the same, re-check!');
  end
end

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
  twlet = zcrosses(ista)*dt;
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

%%%plot the scatter of offsets, accounting for prealignment offset, == true offset
xran = [-loff_max+off1i(2)-1 loff_max+off1i(2)+1];
yran = [-loff_max+off1i(3)-1 loff_max+off1i(3)+1];
offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
cran = [0 lsig];
f1.fig = figure;
f1.fig.Renderer = 'painters';
ax1=gca;
[ax1,torispl,mamp] = plt_decon_imp_scatter(ax1,impindepst,xran,yran,cran,offxran,offyran,...
  sps,50,'mean','tarvl');
scatter(ax1,off1i(2),off1i(3),20,'ks','filled','MarkerEdgeColor','k');
title(ax1,'Grouped Total');

%%%plot the scatter of sources in terms of rela locations
xran = [-4 4];
yran = [-4 4];
cran = [0 lsig/sps];
f2.fig = figure;
f2.fig.Renderer = 'painters';
ax2=gca;
[ax2] = plt_decon_imp_scatter_space(ax2,impindepst,xran,yran,cran,offxran,...
  offyran,sps,50,ftrans,'mean','tarvl');
plot(ax2,xcut,ycut,'k-','linew',2);
title(ax2,'Grouped Total');


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

projang = 135;
[projx,orty,locxyproj] = customprojection(implocst,projang);
dlocxyproj = [diffcustom(projx,1,'forward') diffcustom(orty,1,'forward')];
if ~isempty(dlocxyproj)
  nsep = 1;
  ttype = 'tarvl';

  [f] = plt_customprojdist(tarvlsplst,implocst,dtarvl{nsep},eucdist{nsep},...
      locxyproj,dlocxyproj,projang,sps,ttype);
  hold(f.ax(1),'on');
  plot(f.ax(1),xcut,ycut,'k-','linew',2);
end


%% Remove the small-amplitude, secondary triplets from the grouped result
%convert the sources in terms of the arrival time of zero-crossing to positve peaks' indices
ppkindep = impindep;
for ista = 1: 3
  ppkindep(:,(ista-1)*2+1) = ppkindep(:,(ista-1)*2+1)+ppeaks(ista)-zcrosses(ista);
end
npkindep = impindep;  %negative peaks
for ista = 1: 3
  npkindep(:,(ista-1)*2+1) = npkindep(:,(ista-1)*2+1)+npeaks(ista)-zcrosses(ista);
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

xran = [-loff_max+off1i(2)-1 loff_max+off1i(2)+1];
yran = [-loff_max+off1i(3)-1 loff_max+off1i(3)+1];
offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
cran = [0 lsig];
%%%plot the scatter of offsets, accounting for prealignment offset, == true offset
f1.fig = figure;
f1.fig.Renderer = 'painters';
ax1=gca;
[ax1,torispl,mamp,xbnd,ybnd] = plt_decon_imp_scatter(ax1,impindepst,xran,yran,cran,offxran,offyran,...
  sps,50,'mean','tarvl');
scatter(ax1,off1i(2),off1i(3),20,'ks','filled','MarkerEdgeColor','k');
title(ax1,'Secondary sources removed');

%%%plot the scatter of sources in terms of rela locations
xran = [-4 4];
yran = [-4 4];
cran = [0 lsig/sps];
f2.fig = figure;
f2.fig.Renderer = 'painters';
ax2=gca;
[ax2] = plt_decon_imp_scatter_space(ax2,impindepst,xran,yran,cran,offxran,...
  offyran,sps,50,ftrans,'mean','tarvl');
plot(ax2,xcut,ycut,'k-','linew',2);
title(ax2,'Secondary sources removed');

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

%For each LFE source, get its distance to all other LFEs, maybe we don't care that long separation in time
[dt2all,dloc2all,dist2all] = srcdistall(tarvlsplst,implocst,[0 2*sps]);
% dt2allbst = [dt2allbst; dt2all];
% dloc2allbst = [dloc2allbst; dloc2all] ;
% dist2allbst = [dist2allbst; dist2all];

% %in terms of origin time?
% [imploc0, ~] = off2space002([0 0],sps,ftrans,0);  % a ref source at 0,0
% tcor = round((imploc(:,6)-imploc0(6))*sps);   % travel time difference between each source and a ref source at 0,0
% torispl = impindepst(:,1)-tcor;
% [torisplst, indsort] = sortrows(torispl,1);
% implocst = imploc(indsort, :);
% [dto2all,dloco2all,disto2all] = srcdistall(torisplst,implocst,[0 2*sps]);
% dto2allbst = [dto2allbst; dto2all];
% dloco2allbst = [dloco2allbst; dloco2all] ;
% disto2allbst = [disto2allbst; disto2all];

% %plot euclidean distance between each LFE source to all others
% f = plt_srcdistall(dt2all,dist2all,sps,40/sps,0.1,'km');
% %plot the loc diff between each LFE source to all others
% f = plt_srcdlocall(dloc2all,0.1,'km');

%between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
m = 5;
[dtarvl,dneloc,eucdist] = srcdistNtoNm(tarvlsplst, implocst, m);
% dtarvlnn1all = [dtarvlnn1all; dtarvl{1}];
% dtarvlnn2all = [dtarvlnn2all; dtarvl{2}];
% dtarvlnn3all = [dtarvlnn3all; dtarvl{3}];
% dtarvlnn4all = [dtarvlnn4all; dtarvl{4}];
% dtarvlnn5all = [dtarvlnn5all; dtarvl{5}];
% distarvlnn1all = [distarvlnn1all; eucdist{1} dneloc{1}];  % dneloc{1}(:,2)-dneloc{1}(:,1)
% distarvlnn2all = [distarvlnn2all; eucdist{2} dneloc{2}];
% distarvlnn3all = [distarvlnn3all; eucdist{3} dneloc{3}];
% distarvlnn4all = [distarvlnn4all; eucdist{4} dneloc{4}];
% distarvlnn5all = [distarvlnn5all; eucdist{5} dneloc{5}];

% %plot the loc diff between above source pairs
% f = plt_srcdlocNtoNm(dneloc,0.1,'km');
% %plot the diff time and distance between above source pairs
% f = plt_srcdistNtoNm(dtarvl,eucdist,sps,40/sps,0.1,'km');

projang = 135;
[projx,orty,locxyproj] = customprojection(implocst,projang);
dlocxyproj = [diffcustom(projx,1,'forward') diffcustom(orty,1,'forward')];
if ~isempty(dlocxyproj)
  nsep = 1;
  ttype = 'tarvl';

  [f] = plt_customprojdist(tarvlsplst,implocst,dtarvl{nsep},eucdist{nsep},...
    locxyproj,dlocxyproj,projang,sps,ttype);
  hold(f.ax(1),'on');
  plot(f.ax(1),xcut,ycut,'k-','linew',2);
end

% keyboard


%% Best alignment for the testing window
% 
% %compute the median absolute amplitude and envelope of the same moving window
% %for the moving window at the same station, sensable to use median
% [ir,ramp1,renv1] = Runningampenv(sigpnsta(:,1),mwlen,mwlen-1,'median');
% [~,ramp2,renv2] = Runningampenv(sigpnsta(:,2),mwlen,mwlen-1,'median');
% [~,ramp3,renv3] = Runningampenv(sigpnsta(:,3),mwlen,mwlen-1,'median');
% %looks like using the amplitude and envelope are pretty similar
% %maybe better to mean, otherwise stations PGC will be underrepresented due to a larger amp
% %variation
% ramp = mean([ramp1 ramp2 ramp3], 2);  % use mean or median??
% renv = mean([renv1 renv2 renv3], 2);
% 
% figure
% ax = gca; hold(ax,'on');
% scatter(ax,renv,rcc,10,ir,'filled');
% scatter(ax,median(renv),median(rcc),8,'ko','linew',1);
% [xcnt,ycnt,y1sig] = ranybinx(renv,rcc,'median',10);
% errorbar(ax,xcnt,ycnt,-y1sig,y1sig,'vertical','o','markersize',3,'color',...
%   'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor','k','CapSize',4);
% [rho1,~] = corr(renv,rcc,'Type','Spearman');
% [rho2,~] = corr(renv,rcc,'Type','Kendall');
% text(ax,0.98,0.55,sprintf('S: %.2f',rho1),'unit','normalized',...
%   'HorizontalAlignment','right','fontsize',12);
% text(ax,0.98,0.50,sprintf('K: %.2f',rho2),'unit','normalized',...
%   'HorizontalAlignment','right','fontsize',12);
% colormap(ax,'jet');
% caxis(ax,[0 size(sigpnsta,1)]);
% c=colorbar(ax,'south');
% c.Label.String = sprintf('Samples at %d Hz',sps);
% ylim(ax,[-1 1]);
% xlabel(ax,'Running envelope');
% ylabel(ax,'Running CC');
% hold(ax,'off');
% 
% %get the sources that arrived within the moving window, and estimate the source region size 
% wins = movingwins(1,size(sigpnsta,1)-1,mwlen,mwlen-1);
% nwin = size(wins,1);
% srcregsz = -10*ones(nwin,1);
% nsrci = -10*ones(nwin,1);
% dist = -10*ones(nwin,1);
% distmed = -10*ones(nwin,1);
% for i = 1: nwin
%   win = wins(i,:);  % start and end indices of the moving window
%   srci = srcpair(srcpair(:,1)>=win(1) & srcpair(:,1)<win(2), :);
%   if ~isempty(srci)
%     nsrci(i) = size(srci,1);
%     %how to estimate the 'size' of source region
%     srcregsz(i) = sqrt(polyarea(srci(:,7), srci(:,8)));
%     %find the centroid of the sources and estimate its distance to best alignment
%     srccntrd = [median(srci(:,7)) median(srci(:,8))];
%     dist(i) = sqrt((srccntrd(1)-off1i(2))^2 + (srccntrd(2)-off1i(3))^2);
%     %find the distance from each source to best alignment
%     distsrc = sqrt((srci(:,7)-off1i(2)).^2 + (srci(:,8)-off1i(3)).^2);
%     %the median of the distance
%     distmed(i) = median(distsrc);
% %     srcregsz(i) = sqrt(max(abs(srci(:,7)))^2 + max(abs(srci(:,8)))^2);
% %     srcregsz(i) = sqrt(median(abs(srci(:,7)))^2 + median(abs(srci(:,8)))^2);
%   end
% end
% 
% figure
% subplot(131)
% ax = gca; hold(ax,'on');
% ind = find(srcregsz==-10);
% scatter(ax,renv(ind),rcc(ind),6,[.5 .5 .5]);
% ind = find(srcregsz~=-10);
% scatter(ax,renv(ind),rcc(ind),10,srcregsz(ind),'filled');
% colormap(ax,'jet');
% % caxis(ax,[0 size(sigpnsta,1)]);
% c=colorbar(ax,'south');
% c.Label.String = 'Size of source region (samples)';
% % c.Label.String = 'Median size of source region (samples)';
% ylim(ax,[-1 1]);
% xlabel(ax,'Running envelope');
% ylabel(ax,'Running CC');
% hold(ax,'off');
% 
% subplot(132)
% ax = gca; hold(ax,'on');
% ind = find(dist==-10);
% scatter(ax,renv(ind),rcc(ind),6,[.5 .5 .5]);
% ind = find(dist~=-10);
% scatter(ax,renv(ind),rcc(ind),10,distmed(ind),'filled');
% colormap(ax,'jet');
% % caxis(ax,[0 size(sigpnsta,1)]);
% c=colorbar(ax,'south');
% c.Label.String = 'Median distance to best alignment (samples)';
% % c.Label.String = 'Distance of centroid to best alignment (samples)';
% % c.Label.String = 'Median size of source region (samples)';
% ylim(ax,[-1 1]);
% xlabel(ax,'Running envelope');
% ylabel(ax,'Running CC');
% hold(ax,'off');
% 
% subplot(133)
% ax = gca; hold(ax,'on');
% ind = find(dist==-10);
% scatter(ax,renv(ind),rcc(ind),6,[.5 .5 .5]);
% ind = find(dist~=-10);
% scatter(ax,renv(ind),rcc(ind),10,nsrci(ind),'filled');
% colormap(ax,'jet');
% % caxis(ax,[0 size(sigpnsta,1)]);
% c=colorbar(ax,'south');
% c.Label.String = 'Number of sources';
% ylim(ax,[-1 1]);
% xlabel(ax,'Running envelope');
% ylabel(ax,'Running CC');
% hold(ax,'off');
% 
% figure
% yyaxis('left');
% plot(sigpnsta(:,1),'r-'); hold on
% plot(sigpnsta(:,2),'b-');
% plot(sigpnsta(:,3),'k-'); %,'linew',0.5
% axranexp(gca,6,10);
% xlabel(sprintf('Samples at %d Hz',sps));
% ylabel('Amplitude');
% yyaxis('right');
% % plot(ircc,rcc,'o','color',[.7 .7 .7],'markersize',2);
% scatter(ircc,rcc,2,renv);  % scale it to the amplitude range
% colormap('jet');
% caxis([min(renv) max(renv)]);
% ylabel('Running CC','FontSize',10);
% ylim([-1.1 1.1]);
% 
% %distribution of synthetic sources, ground truth
% spsscale = sps/40;  
% loff_max = 4*spsscale;  % maximum allowed shift is 4 samples at 40 Hz, it is also true from the synthetic distribution 
% offxran = [-loff_max loff_max]; 
% offyran = [-loff_max loff_max];
% lsig = size(sigpnsta,1);
% cran = [0 lsig];
% [f] = plt_decon_imp_scatter(srcpair,offxran,offyran,cran,sps,1,'mean');
% title(gca,sprintf('Synthetic sources, using data of %d Hz',sps));
% 
% 
% %% 
% figure;
% colorn = {'r','b','k'};
% for i = 1: nsta
%   subplot(5,1,i)
%   ax = gca;
%   hold(ax,'on');
%   ax.Box = 'on';
%   grid(ax,'on');
%   indtemp = find(sigdecon(:,i)>0);  % independent impulses from other stations
%   stem(ax,indtemp,sigdecon(indtemp,i),'color',[.8 .8 .8],'MarkerSize',4); 
%   text(ax,0.9,0.9,num2str(length(indtemp)),'unit','normalized','HorizontalAlignment','right',...
%     'FontSize',10,'color',[.8 .8 .8]);
%   indtemp = find(impindep(:,(i-1)*2+1)>0);  % that can be individually paired with 3-sta among independent
%   stem(ax,impindep(indtemp,(i-1)*2+1),impindep(indtemp,(i-1)*2+2),'color',colorn{i},...
%     'MarkerSize',4);
%   text(ax,0.9,0.75,num2str(length(indtemp)),'unit','normalized','HorizontalAlignment','right',...
%     'FontSize',10,'color',colorn{i});
%   text(ax,0.05,0.9,strcat(stas(i,:)),'unit','normalized','HorizontalAlignment','left',...
%     'FontSize',12);
%   xlim(ax,[0 lsig]); hold(ax,'off');
%   if i == 1
%     title('Triplet impulses from grouped independent deconvolution');
%   end
% end
% subplot(5,1,4)
% ax = gca;
% yyaxis(ax,'left');
% hold(ax,'on');
% ax.Box = 'on';
% grid(ax,'on');
% stem(ax,impindep(:,1),impindep(:,2), 'r-','MarkerSize',4);
% stem(ax,impindep(:,3),impindep(:,4), 'b-','MarkerSize',4);
% stem(ax,impindep(:,5),impindep(:,6), 'k-','MarkerSize',4);
% text(ax,0.05,0.9,sprintf('Number of triplets: %d', length(indpair)),'fontsize',10,...
%   'Units','normalized');
% yyaxis(ax,'right');
% plot(ircc,rcc,'o','color',[.5 .5 .5],'markersize',2);
% ylim(ax,[-1 1]);
% xlim(ax,[0 lsig]);
% hold(ax,'off');
% 
% %the synthetic sources, ground truth
% subplot(5,1,5)
% ax = gca;
% yyaxis(ax,'left');
% hold(ax,'on');
% ax.Box = 'on';
% stem(ax,srcpair(:,1),srcpair(:,2), 'r-','MarkerSize',4);
% stem(ax,srcpair(:,3),srcpair(:,4), 'b-','MarkerSize',4);
% stem(ax,srcpair(:,5),srcpair(:,6), 'k-','MarkerSize',4);
% text(ax,0.05,0.9,sprintf('Number of triplets: %d',size(srcpair,1)),'fontsize',10,...
%   'Units','normalized');
% yyaxis(ax,'right');
% plot(ircc,rcc,'o','color',[.5 .5 .5],'markersize',2);
% ylim(ax,[-1 1]);
% xlim(ax,[0 lsig]);
% title(ax,'Triplet impulses from true synthetic sources');
% hold(ax,'off');
