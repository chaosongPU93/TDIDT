function allsyn = decon_synth_regandnoi_fn(srcregion,reg,perc,sat,distr,Twin,...
  tdura,distrloc,physicalsize,testsrcflag,testfreqflag,pltdataflag,pltgtflag,...
  pltsrcflag,pltsrc4thflag,normflag,tempflag,ftrans,nrun)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --The function that is called to read in synthetic seismograms generated
% from different region sizes and saturation levels, and then implement
% deconvolution to those synthetics.
% --First read and deconvolve the synthetic seismograms from
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
% --2024/05/03, for better comparison with data, we need to generate a total
% length of 3 hours of synthetics. But decon to such a long record takes to much
% time, we generate several synthetics and decon each one, then combine them. As
% randomness is embedded, this should have no negative effects.
% --See also 'decon_synth_onespot_fn.m'.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/05/03
% Last modified date:   2024/05/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% default settings
defval('srcregion','ellipse');  %specify shape of the source region
defval('reg',[1.75 1.25]); %source region size
defval('perc',0.6); %noise level
defval('sat',1); %saturation level
defval('distr','UN'); %specify the amplitude-frequency (counts) distribution as uniform
defval('Twin',0.5*3600+3+2*ceil(2048/160)); %length of each simulation
defval('tdura',0.25); %duration of templates
defval('distrloc','uniform'); %uniformly random for source location
defval('physicalsize',0); %specify if considering the physical size of each source
defval('testsrcflag',0); %flag for validating if ground truth of sources can recover the record
defval('testfreqflag',0); %flag for validing if the spectral shapes of data and templates are similar
defval('pltdataflag',0);%flag for plot the data
defval('pltgtflag',0);%flag for plot the ground truth distribution
% defval('pltsrcflag1',0);%flag for plot the decon src distribution after grouping
defval('pltsrcflag',0);%flag for plot the decon src distribution after removing 2ndary src
defval('pltsrc4thflag',0);%flag for plot the decon src distribution after checking at 4th stas
defval('normflag',0); %whether to normalize templates
defval('tempflag','chao'); %which templates to use
defval('ftrans','interpchao');  %specify regime for transformation from time offset to map location
defval('nrun',1);%flag for plot the decon src distribution after checking at 4th stas

%% Initialization
format short e   % Set the format to 5-digit floating point
% clear
% clc
% close all

%%%Flag to indicate if it is necessary to recalculate everything
% flagrecalc = 0;
% flagrecalc = 1;

% tdura = 0.25;  %must be consistent with what synthetics actually used
% tdura = 0.4;

% % if ~flagrecalc
% %
% % else
% %% for easy testing

%% Initialization
% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');

fam = '002';   % family number

stas=['PGC  '
  'SSIB '
  'SILB '
  %   'LZB  '
  %   'TWKB '
  %   'MGCB '
  'KLNB '
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

%%%load the empirical model param of time offset 14 for each addtional stations
off14mod = load(strcat(rstpath, '/MAPS/timeoff_planefitparam_4thsta_160sps'),'w+');

%% prepare templates (Green's functions), from 'lfetemp002_160sps.m' or Allan's templates
adatapath = '/home/data2/chaosong/matlab/allan/matfils/';  %path for Allan's data

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
[off12conf,off13conf,ccf] = constrained_cc_interp(tmpwletf(:,1:3)',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
offwlet1if(1) = 0;
offwlet1if(2) = round(off12conf);
offwlet1if(3) = round(off13conf);
if nsta>3
  for ista = 4: nsta
    [mcoef,offwlet1if(ista)] = xcorrmax(tmpwletf(:,1), tmpwletf(:,ista), mshiftadd, 'coeff');
  end
end
%%%for the broadband templates as well
[off12con,off13con,cc] = constrained_cc_interp(tmpwlet(:,1:3)',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
offwlet1i(1) = 0;
offwlet1i(2) = round(off12con);
offwlet1i(3) = round(off13con);
if nsta>3
  for ista = 4: nsta
    [mcoef,offwlet1i(ista)] = xcorrmax(tmpwlet(:,1), tmpwlet(:,ista), mshiftadd, 'coeff');
  end
end

%%%automatically find the rough zero-crossing time, whose abs. value is closest to 0, whether + or -
[~,imin] = min(tmpwletf(:,1));
[~,imax] = max(tmpwletf(:,1));
[~,zcsta1] = min(abs(tmpwletf(imin:imax,1)));
zcsta1 = zcsta1+imin-1;
[~,imin] = min(tmpwlet(:,1));
[~,imax] = max(tmpwlet(:,1));
[~,zcsta2] = min(abs(tmpwlet(imin:imax,1)));
zcsta2 = zcsta2+imin-1;
greenlen = pow2(9)*sps/40;
green = zeros(greenlen,nsta); % no bandpass
greenf = zeros(greenlen,nsta);  % bandpassed version
ppeaks = zeros(nsta,1); % positive peaks
npeaks = zeros(nsta,1); % negative peaks
zcrosses = zeros(nsta,1);
ppeaksf = zeros(nsta,1); % positive peaks
npeaksf = zeros(nsta,1); % negative peaks
zcrossesf = zeros(nsta,1);
for ista = 1: nsta
  %cut according to the zero-crossing and the time shift from the constrained CC
  green(:,ista) = tmpwlet(zcsta2+8*sps-greenlen+1-offwlet1i(ista): zcsta2+8*sps-offwlet1i(ista), ista);
  greenf(:,ista) = tmpwletf(zcsta1+8*sps-greenlen+1-offwlet1if(ista): zcsta1+8*sps-offwlet1if(ista), ista);
  %detrend again for caution
  green(:,ista) = detrend(green(:,ista));
  greenf(:,ista) = detrend(greenf(:,ista));
  if normflag
    %normalize by max amp
    green(:,ista) = green(:,ista)/max(abs(green(:,ista)));    % normalize
    greenf(:,ista) = greenf(:,ista)/max(abs(green(:,ista)));    % normalize
  end
  
  %re-find the zero-crossing as the template length has changed
  [~,imin] = min(green(:,ista));
  [~,imax] = max(green(:,ista));
  [~,zcrosses(ista)] = min(abs(green(imin:imax,ista)));
  zcrosses(ista) = zcrosses(ista)+imin-1;
  ppeaks(ista) = imax;
  npeaks(ista) = imin;
  
  %for bandpassed templates
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
[off12con,off13con,cc] = constrained_cc_interp(green(:,1:3)',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
if ~(off12con==0 && off13con==0)
  disp('Broadband templates are NOT best aligned \n');
end
if nsta>3
  for ista = 4: nsta
    [mcoef,mlag] = xcorrmax(green(:,1), green(:,ista), mshiftadd, 'coeff');
    if mlag~=0   % offset in samples
      fprintf('Broadband templates are NOT best aligned at %s \n',stas(ista,:));
    end
  end
end

amprat(1,:) = minmax(greenf(:,1)')./minmax(greenf(:,2)');	% amp ratio between max at sta 3 and 2 or min
amprat(2,:) = minmax(greenf(:,1)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min
amprat(3,:) = minmax(greenf(:,2)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min
spread = range(greenf);   % range of the amp of template

%% load synthetic seismograms

fracelsew=0; %0.25 %0.5; %0.6; %The ratio of elsewhere events to local events.  Zero for Discrete Ide.

%%%diameter of physical size
if physicalsize
  diam=0.15;% 0.5; %0.6; %
else
  diam=0;
end

%%%file name prefix of synthetics
if strcmp(srcregion,'ellipse')
  xaxis = reg(1); %semi-long axis
  yaxis = reg(2); %semi-short axis
  % xaxis=1.75; %axis length of the same ellipse of my 4-s catalog
  % yaxis=1.25;
  shiftor=[0.2 0.2]; %(in km) %center of the same ellipse of my 4-s catalog
  [xcut,ycut] = ellipse_chao(shiftor(1),shiftor(2),xaxis,yaxis,0.01,45,shiftor);
  fname = ['/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),'_',...
    num2str(xaxis),'-',num2str(yaxis),'.diam',num2str(diam),'.noi',num2str(perc),...
    '.else',num2str(fracelsew,2),'nsat'];
elseif strcmp(srcregion,'circle')
  shiftor=[0.2 0.2]; %(in km) %center of the same ellipse of my 4-s catalog
  radi = reg(1); %radius
  [xcut,ycut] = circle_chao(shiftor(1),shiftor(2),radi,0.01);
  fname = ['/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
    '_',num2str(radi),'.diam',num2str(diam),'.noi',num2str(perc),...
    '.else',num2str(fracelsew,2),'nsat'];
end
% fname = '/synthetics/STAS.UN.160sps.ell_2-1.25.diam0.3.else0nsat';

%%%initialize array for combined results from all runs
irccranall = [];
windowsall = [];
off1iwall = [];
ccaliwall = [];
ccwpairall = [];
mrccwpairall = [];
rcccatall = [];
off1iall = [];
rccall = [];
mrccall = [];
ccall = [];
impgtall = [];
nitall = [];
impgrpall = [];
nsrcgrpall = [];
impall = [];
nsrcall = [];
mprojxnn1all = [];
mprojx2allall = [];
srcamprall = [];
msrcamprall = [];
madsrcamprall = [];
mpsrcamprsall = [];
madpsrcamprsall = [];
mnsrcamprsall = [];
madnsrcamprsall = [];
imp4thall = [];
nsrc4thall = [];
diffoff14trall = [];
mprojxnn14thall = [];
mprojx2all4thall = [];
srcampr4thall = [];
msrcampr4thall = [];
madsrcampr4thall = [];
mpsrcamprs4thall = [];
madpsrcamprs4thall = [];
mnsrcamprs4thall = [];
madnsrcamprs4thall = [];


%%%loop for each run (each segment of synthetics to be combined later)
for irun = 1: nrun
  fprintf('part: %d/%d \n',irun,nrun);
  
  if tdura == 0.25
    %%%load synthetics of certain saturation level
    STAopt = load(strcat(workpath,fname,num2str(sat),'tdura',...
      num2str(tdura),'T',num2str(round(Twin)),'p',num2str(irun)));
    %%%load sources
    synsrc = load(strcat(workpath,fname,num2str(sat),'tdura',...
      num2str(tdura),'T',num2str(round(Twin)),'p',num2str(irun),'_sources'));
    %%%load starting indices of added sources at sta 1
    synsrcstind = load(strcat(workpath,fname,num2str(sat),'tdura',...
      num2str(tdura),'T',num2str(round(Twin)),'p',num2str(irun),'_stind'));
  end
  
  if strcmp(distrloc, 'uniform')
    if strcmp(srcregion,'ellipse')
      xygrid = load([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),...
        'sps.',srcregion(1:3),'_',num2str(xaxis),'-',num2str(yaxis),'.diam',...
        num2str(diam),'.noi',num2str(perc),'.else',num2str(fracelsew,2),'nsat',...
        num2str(sat),'tdura',num2str(tdura),...
        'T',num2str(round(Twin)),'p',num2str(irun),'_grd']);
    elseif strcmp(srcregion,'circle')
      xygrid = load([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',...
        srcregion(1:3),'_',num2str(ll(1)),num2str(ll(2)),num2str(ur(1)),...
        num2str(ur(2)),'.diam',num2str(diam),...
        '.noi',num2str(perc),'.else',num2str(fracelsew,2),'nsat',...
        num2str(sat),'tdura',num2str(tdura),...
        'T',num2str(round(Twin)),'p',num2str(irun),'_grd']);
    end
    tmp = xygrid(synsrc(:,2),:);
    synsrc = [synsrc(:,1) tmp(:,1:4) ones(length(tmp),1)];  %[indtarvl, off12, off13, loce, locn, amp]
  end
  % keyboard
  
  %%%some params related to the synthetics setting
  %       Twin=0.5*3600+3+2*ceil(greenlen/sps); %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
  winlen=Twin*sps+1;
  skiplen=greenlen;
  %%%Here the noise is 'uniform' in time!
  % noistd = 5e-2;
  %       noistd = 2.0e-4;
  noistd = 0; %NOISE-FREE
  rng('default');
  %     seed=(irun-1)*nreg+ireg;
  %     rng(seed);
  synth=noistd*(randn(winlen+greenlen+2*10,nsta)-0.5); %+2*10 a little extra, for jiggering 2nd & 3rd stations.% for ista=1:nsta
  
  %% testing, extract and validate the added impulses of template
  if testsrcflag
    srcvalidatev2(sps,zcrosses,synsrc,synsrcstind,stas,skiplen,greenlen,...
      STAopt,synth,winlen,green);
    
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
  
  %% break into short windows
  %some params
  bufsec = 1;
  msftaddm = bufsec*sps;  %buffer range for later CC alignment, +1 for safety
  rccmwsec = 0.5;
  rccmwlen = rccmwsec*sps;  %window length for computing RCC
  overshoot = rccmwlen/2; %RCC points to the center of computing window, so extra room is needed
  
  indst = 1+msftaddm;
  inded = size(optseg,1)-msftaddm;
  subwsec = 25;
  subwlen = subwsec*sps;
  %since the rcc would lose rccmwlen/2 at both ends, this results in overlapping of 'rccmwlen' in rcc
  %across consecutive windows; if use 'ovlplen' of rccmwlen, then rcc has no overlapping at all
  %       ovlplen = rccmwlen*2;
  ovlplen = rccmwlen;
  windows = movingwins(indst,inded,subwlen,ovlplen,0);
  nwin =  size(windows,1);
  
  off1iw = zeros(nwin,nsta);  % the best alignment between sta2, sta3 wrt sta1 for each subwin
  ccaliw = zeros(nwin,1+nsta-3);  % CC value using the best alignment, including 4th stas
  
  ircccat = [];   % concatenated indices of RCC
  irccran = zeros(nwin,2);  % start and end indices (range) of RCC of all subwins
  rcccat = [];  % average concatenated RCC
  rcc1icat = [];  % concatenated RCC between sta 1 and 4
  rccpaircat = [];  % concatenated RCC between each station pair, order is 12, 13, 23
  mrccwpair = []; % median of concatenated RCC between each station pair
  ccwpair = []; % 0-lag overall cc of each subwin, between each station pair, order is 12, 13, 23
  ccw1i = []; % same as above, but between sta 1 and 4
  
  for iwin = 1: nwin
    isubwst = windows(iwin,1);
    isubwed = windows(iwin,2);
    %         isubwst = windows(iwin,1)-max(floor(indst-overshoot-msftaddm),1);
    %         isubwed = windows(iwin,2)-max(floor(indst-overshoot-msftaddm),1);
    %         if iwin == 1
    %           isubwst = isubwst-overshoot;
    %         end
    %         if iwin == nwin
    %           isubwed = isubwed+overshoot;
    %         end
    
    %align records
    optcc = detrend(optseg(isubwst: isubwed,:));
    msftadd = (round(max(abs([off12ran off13ran])))+1)*sps/40;  %+1 for safety
    loffmax = 4*sps/40;
    ccmid = ceil(size(optcc,1)/2);
    ccwlen = round(size(optcc,1)-2*(msftadd+1));  % minus ensures successful shifting of records
    ccmin = 0.01;  % depending on the length of trace, cc could be very low
    iup = 1;    % times of upsampling
    [off12con,off13con,ccaliw(iwin,1)] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
      ccwlen,msftadd,loffmax,ccmin,iup);
    % if a better alignment cannot be achieved, use 0,0
    if off12con == msftadd+1 && off13con == msftadd+1
      fprintf('Short window %d cannot be properly aligned, double-check needed \n',iwin);
      off12con = 0;
      off13con = 0;
      cc12 = xcorr(optcc(:,1), optcc(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
      cc13 = xcorr(optcc(:,1), optcc(:,3),0,'normalized');
      cc23 = xcorr(optcc(:,2), optcc(:,3),0,'normalized');
      ccaliw(iwin,1) = (cc12+cc13+cc23)/3;
    end
    off1iw(iwin,1) = 0;
    off1iw(iwin,2) = round(off12con);
    off1iw(iwin,3) = round(off13con);
    
    for ista = 4: nsta
      [ccaliw(iwin,ista-2),off1iw(iwin,ista)] = xcorrmax(optcc(:,1),optcc(:,ista), 1.5*msftadd, 'coeff');
    end
    
    %Align records
    optdat = [];  % win +/-3 s, segment of interest,first 1s will be tapered
    optdat(:, 1) = optseg(isubwst: isubwed, 1); % time column
    for ista = 1: nsta
      optdat(:, ista) = optseg(isubwst-off1iw(iwin,ista): isubwed-off1iw(iwin,ista), ista);
    end
    
    subw = zeros(size(optdat,1), nsta);
    for ista = 1: nsta
      tmp = optdat(:,ista); %best aligned, filtered
      %detrend and taper only the data, NOT the noise
      tmp = detrend(tmp);
      %%%2022/06/06, do NOT taper whatsoever!!
      subw(:,ista) = tmp;
    end
    %compute running CC between 3 stations
    [irccw,rccw12,rccw13,rccw23] = RunningCC3sta(subw,rccmwlen);
    rccw = (rccw12+rccw13+rccw23)/3;
    irccw = irccw + windows(iwin,1) - windows(1,1);   %convert to global index
    %         if iwin == 1
    irccw = irccw - overshoot;
    %         end
    %         plot(irccw, rccw);
    
    ircccat = [ircccat; irccw];
    irccran(iwin,:) = [irccw(1) irccw(end)];
    rcccat = [rcccat; rccw];
    rccpaircat = [rccpaircat; rccw12 rccw13 rccw23];
    mrccwpair = [mrccwpair; median(rccw12) median(rccw13) median(rccw23)];
    
    rccw1i = [];  % between 4th stas and 1st sta
    for ista = 4:nsta
      [~,rccw1i(:,ista-3)] = RunningCC(subw(:,1), subw(:,ista), rccmwlen);
    end
    rcc1icat = [rcc1icat; rccw1i];
    
    ccw12 = xcorr(subw(:,1), subw(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
    ccw13 = xcorr(subw(:,1), subw(:,3),0,'normalized');
    ccw23 = xcorr(subw(:,2), subw(:,3),0,'normalized');
    ccwpair = [ccwpair; ccw12 ccw13 ccw23];
    tmp = zeros(1,nsta-3);
    for ista = 4:nsta
      tmp(1,ista-3) = xcorr(subw(:,1), subw(:,ista),0,'normalized');
    end
    ccw1i = [ccw1i; tmp];
  end
  lsig = Twin*sps-(greenlen+msftaddm*2+overshoot*2-2);
  
  %if only use the mean RCC from pair 12 and 13
  rcccat = mean(rccpaircat(:,[1 2]), 2);
  % rcccatall = [rcccatall; rcccat];
  
  %% Best alignment for the whole window
  %%%obtain a single best alignment based on the entire win
  optcc = detrend(optseg(1+msftaddm: end-msftaddm, :));
  %       msftadd = 10*sps/40;
  ccmid = ceil(size(optcc,1)/2);
  ccwlen = round(size(optcc,1)-2*(msftadd+1));  % minus ensures successful shifting of records
  ccmin = 0.01;  % depending on the length of trace, cc could be very low
  iup = 1;    % times of upsampling
  [off12con,off13con,ccali,iloopoff,loopoff] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
    ccwlen,msftadd,loffmax,ccmin,iup);
  % if a better alignment cannot be achieved, use 0,0
  if off12con == msftadd+1 && off13con == msftadd+1
    fprintf('Whole window cannot be properly aligned, double-check needed \n');
    off12con = 0;
    off13con = 0;
    cc12 = xcorr(optcc(:,1), optcc(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
    cc13 = xcorr(optcc(:,1), optcc(:,3),0,'normalized');
    cc23 = xcorr(optcc(:,2), optcc(:,3),0,'normalized');
    ccali = (cc12+cc13+cc23)/3;
  end
  off1i = zeros(1,nsta);
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
  
  %       %%%%%%%%%%%%%%%%%%%%%%%%%%
  %       %if you want to avoid the case when the alignment is way off the centroid
  %       %by chance while the saturation level is low, you can force it to be the
  %       %an average location, this reference value is from the abs location of the
  %       %the centroid (0.2,0.2), and prediction of off14 from plane fit model
  %       off1i(2:3) = [2 2];
  %       [~,off1i(4)] = pred_tarvl_at4thsta(stas(4,:),off1i(2),off1i(3));
  %       %%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % off1iall = [off1iall; off1i];
  
  %%%Align and compute the RCC based on the entire win, and take that as the input signal!
  optdat = [];  % win segment of interest
  for ista = 1: nsta
    optdat(:, ista) = optseg(1+msftaddm-off1i(ista): end-msftaddm-off1i(ista), ista);
  end
  
  %location of the whole-win best alignment
  loc0 = off2space002([off1i(2) off1i(3)],sps,ftrans,0);
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
  mrcc = median(rcc);
  % rccall = [rccall; rcc];
  
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
  cc = (cc12+cc13)/2;  %if only use the pair 12 and 13
  
  %%
  %if you want a qucik look of the data that feeds to the deconvolution
  if pltdataflag
    figure; hold on
    lsig = size(sigsta,1);
    %         if nsta == 3
    plot((1:lsig)/sps,sigsta(:,1),'r');
    plot((1:lsig)/sps,sigsta(:,2),'b');
    plot((1:lsig)/sps,sigsta(:,3),'k');
    %         else
    %           color = jet(nsta);
    %           for ista = 1: nsta
    %             p(ista)=plot((1:lsig)/sps,sigsta(:,ista),'Color',color(ista,:));
    %             %       label{i}=stas(i,:);
    %           end
    %         end
    ax = gca;
    axsym(ax);
    plot((1:lsig)/sps,rcc*ax.YLim(2),'o','color',[.6 .6 .6],'markersize',2);
    text(0.95,0.1,sprintf('Saturation: %.1f',sat),'Units','normalized','HorizontalAlignment','right');
    xlim([0 20]);
    legend(p,stas);
  end
  
  %% ground truth of conrtibuting sources AFTER alignment
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
  [density1d, inddup] = density_pixel(impgt(:,7),impgt(:,8));
  ampgtsum = sum_at_indices(ampgt,inddup);
  [impgtloc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
  density1d = [impgtloc(:,1:2) density1d(:,3)];
  ampgtsum1d = sortrows([impgtloc(:,1:2) ampgtsum], 3);
  
  %% if you want to plot the ground truth
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
    scatter(ax,off1i(2),off1i(3),50,'ks','filled','MarkerEdgeColor','k');
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
    scatter(ax,loc0(1),loc0(2),50,'ks','filled','MarkerEdgeColor','k');
    title(ax,'Ground truth');
    
    %%%plot the cumulative density and summed amp of detections
    cstr = {'# detections / pixel'; 'amp sum / pixel'};
    [f] = plt_sum_pixel(density1d,ampgtsum1d,[-4 4],[-4 4],20,cstr,'o','linear');
    hold(f.ax(1),'on');
    plot(f.ax(1),xcut,ycut,'k-','linew',2);
    scatter(f.ax(1),loc0(1),loc0(2),50,'ks','filled','MarkerEdgeColor','k');
    hold(f.ax(2),'on');
    plot(f.ax(2),xcut,ycut,'k-','linew',2);
    supertit(f.ax,'Ground truth');
    
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
    tlen = ceil(lsig/sps);
    nit_max = round(1.5*1/tdura*(tlen));  % max numer of iterations
    nimp_max = round(1/tdura*(tlen));%a single peak is ~20 samples wide; maybe a little less (at 100 sps). ~0.4s, 1/0.4=2.5
    fpltit = 0;  % plot flag for each iteration
    fpltend = 0;  % plot flag for the final iteration
    fpltchk = 0; % plot flag for intermediate computations
    
    [sigdecon(:,ista),pred(:,ista),res,dresit,mfitit{ista},ccchgit,ampit{ista},nit(1,ista),fighdl] = ...
      iterdecon(sig,wlet,rcccat,noi,[],dt,twlet,width,dres_min,mfit_min,nit_max,nimp_max,...
      fpltit,fpltend,fpltchk);
    
    if fpltend
      ax = fighdl{2}.ax(1);
      hold(ax,'on');
      text(ax,0.05,0.85,stas(ista,:),'unit','normalized');
      hold(ax,'off');
    end
    
    fprintf('%d-th sta; nit=%d \n',ista,nit(1,ista));
    
  end
  
  %% Group nearest impulses from different stations into triplets, using moving searching range
  spsscale = sps/40;
  loff_max = 6*spsscale;  % maximum allowed shift is 4 samples at 40 Hz, it is also true from the synthetic distribution
  %note the output 'impindep' gives the arrival index of impulse at each station, after
  %alignment based upon the entire window 'off1i', and the last three cols are the arrival time
  %difference, NOT the true location yet!
  refsta = 1;
  %     [impindep,imppairf,indpair,sharp] = groupimptripdecon(sigdecon,ampit,rcc,loff_max,refsta);
  [impindep,imppairf,indpair,sharp] = groupimptripdecon_ref_nooff23(sigdecon,ampit,irccran,rcccat,...
    off1i(1:3),off1iw(:,1:3),loff_max,refsta);
  
  %note here 'impindepst' inherits the first 6 cols from 'impindep', but the last three cols
  %are adjusted from arrival time difference to the true location offset accounting for the best
  %alignment upon the entire window that is also used in grouping!
  impindep(:,7:8) = impindep(:,7:8)+repmat([off1i(2) off1i(3)],size(impindep,1),1); %account for prealignment
  impindepst = sortrows(impindep,1);
  impgrp = impindepst;
  impgrp(:,[1 3 5]) = impgrp(:,[1 3 5])+lsig*(irun-1);
  
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
  imp = impindepst;
  imp(:,[1 3 5]) = imp(:,[1 3 5])+lsig*(irun-1);
  
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
  mprojxnn1 = median(abs(dprojxy(:,1))); %median consecutive dist along min-scatter after 2nd src removal
  
  %For each LFE source, get its distance to all other LFEs, maybe we don't care that long separation in time
  [~,dprojxy2all] = srcdistall(tarvlsplst,projxy,[0 2*sps]);
  mprojx2all = median(abs(dprojxy2all(:,1))); %median dist along min-scatter to all sources after 2nd src removal
  
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
  
  msrcampr = median(log10(srcampr), 1);
  madsrcampr = mad(log10(srcampr), 1, 1);
  mpsrcamprs = median(log10(psrcamprs), 1);
  madpsrcamprs = mad(log10(psrcamprs), 1, 1);
  mnsrcamprs = median(log10(nsrcamprs), 1);
  madnsrcamprs = mad(log10(nsrcamprs), 1, 1);
  
  %what is the deviation of amp ratio from the median for each source?
  lndevsrcampr = srcampr-median(srcampr, 1); % in linear scale
  lgdevsrcampr = log10(srcampr)-median(log10(srcampr), 1); % in log scale, note that log2+log5=log10, so this means a ratio
  
  %%%if you want to plot the deconvolved sources
  if pltsrcflag
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
    text(ax,0.98,0.95,sprintf('Satur=%.1f',sat),'Units','normalized',...
      'HorizontalAlignment','right');
    text(ax,0.02,0.05,sprintf('%.2f',mprojx2all),'Units','normalized',...
      'HorizontalAlignment','left');
    plot(ax,xcut,ycut,'k-','linew',2);
    title(ax,'Secondary sources removed');
    
    %plot the cumulative density and summed amp of detections
    [density1d, inddup] = density_pixel(impindepst(:,7),impindepst(:,8)); %count at each unique loc
    [imploc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
    density1d = [imploc(:,1:2) density1d(:,3)];
    amp = mean(impindepst(:,[2 4 6]),2); %amp for all LFE catalog
    ampsum = sum_at_indices(amp,inddup);
    ampsum1d = sortrows([imploc(:,1:2) ampsum], 3);
    cstr = {'# detections / pixel'; 'amp sum / pixel'};
    [f] = plt_sum_pixel(density1d,ampsum1d,[-4 4],[-4 4],20,cstr,'o','linear');
    hold(f.ax(1),'on');
    plot(f.ax(1),xcut,ycut,'k-','linew',2);
    hold(f.ax(2),'on');
    plot(f.ax(2),xcut,ycut,'k-','linew',2);
    supertit(f.ax,'Secondary sources removed');
    
    %plot distance between srcs N&N-m along the projected direction
    if ~isempty(dprojxy)
      nsep = 1;
      [f] = plt_customprojdist(tarvlsplst,implocst,dtarvl{nsep},eucdist{nsep},...
        projxy,dprojxy,projang,sps,'tarvl');
      hold(f.ax(1),'on');
      plot(f.ax(1),xcut,ycut,'k-','linew',2);
    end
    
    %plot histograms of source amp
    f = initfig(12,5,1,3); %initialize fig
    f.ax(1:3) = plt_deconpk_rat_comb4th(f.ax(1:3),srcampr,impindepst,'k','hist');
    supertit(f.ax,'Secondary sources removed');
  end
  % keyboard
  
  %% 2ndary src removed, prediction of impulse tarvl at 4th sta given sources and empirical off14-src relation
  %%%carry out 'deconvolution' at 4th stations as well for the tarvl and amp
  modname = 'timeoff_plfit_4thsta_160sps.mat';
  planefit = load(strcat(rstpath, '/MAPS/',modname));
  
  pred4diff = [];
  off14pred = [];
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
    
    [sigdecon(:,ista),pred,res,dresit,mfitit{ista},ampit{ista},fighdl] = ...
      iterdecon_4thsta(sig,wlet,irccran,rcc1icat(:,ista-3),[],...
      dt,twlet,impindep,stas(ista,:),off1i(ista),off1iw(:,ista),offmax,...
      fpltit,fpltend,fpltchk);
    
    ampiti = ampit{ista};
    impindep(:,9+(ista-4)*2+1) = ampiti(:,1); %impulse arrival
    impindep(:,9+(ista-3)*2) = ampiti(:,2); %impulse amp
    if ista == nsta
      pred4diff(:,ista-3) = ampiti(:,end-1);  %difference between found peak and predicted arrival
      off14pred(:,ista-3) = ampiti(:,end);  %predicted off14
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
  
  %%% further ELIMINATE sources that fail the check at 4th stations
  trust4th = 4; % trust KLNB the most among all 4th stations
  indremove = find(impindep(:,9+(trust4th-4)*2+1)==0 & impindep(:,9+(trust4th-3)*2)==0);
  pred4difftr = pred4diff(setdiff(1:size(pred4diff,1),indremove),trust4th-3);
  off14predtr = off14pred(setdiff(1:size(pred4diff,1),indremove),trust4th-3);
  impindep(indremove,:) = [];
  ppkindep(indremove, :) = [];
  npkindep(indremove, :) = [];
  [impindepst,indsort] = sortrows(impindep,1);
  pred4difftr = pred4difftr(indsort);
  imp4th = impindepst;
  imp4th(:,[1 3 5]) = imp4th(:,[1 3 5])+lsig*(irun-1);
  
  %% error in off14 between prediction and observation, expect to be Gaussian around 0
  %the actual 'off14' based on the impulses found by deconvolution
  off14 = impindep(:,1)-impindep(:,9+(trust4th-4)*2+1); %after prealignment
  diffoff14tr = off14predtr-off14;  %error in off14 between prediction and observation
  diffoff14tr = diffoff14tr(indsort);  %sorted to have the same order as saved impulses
  % diffoff14trall = [diffoff14trall; diffoff14tr];
  
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
  
  msrcampr4th = median(log10(srcampr4th), 1);
  madsrcampr4th = mad(log10(srcampr4th), 1, 1);
  mpsrcamprs4th = median(log10(psrcamprs4th), 1);
  madpsrcamprs4th = mad(log10(psrcamprs4th), 1, 1);
  mnsrcamprs4th = median(log10(nsrcamprs4th), 1);
  madnsrcamprs4th = mad(log10(nsrcamprs4th), 1, 1);
  
  %what is the deviation of amp ratio from the median for each source?
  lndevsrcampr4th = srcampr4th-median(srcampr4th, 1); % in linear scale
  lgdevsrcampr4th = log10(srcampr4th)-median(log10(srcampr4th), 1); % in log scale, note that log2+log5=log10, so this means a ratio
  
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
  mprojxnn14th = median(abs(dprojxy(:,1))); %median consecutive dist along min-scatter after 2nd src removal
  
  %For each LFE source, get its distance to all other LFEs, maybe we don't care that long separation in time
  [~,dprojxy2all] = srcdistall(tarvlsplst,projxy,[0 2*sps]);
  mprojx2all4th = median(abs(dprojxy2all(:,1))); %median dist along min-scatter to all sources after 2nd src removal
  
  %%%if you want to plot the deconvolved sources
  if pltsrc4thflag
    %       %plot the scatter of offsets, accounting for prealignment offset, == true offset
    %       xran = [-loff_max+off1i(2)-1 loff_max+off1i(2)+1];
    %       yran = [-loff_max+off1i(3)-1 loff_max+off1i(3)+1];
    %       offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
    %       offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
    %       cran = [0 lsig];
    %       f.fig = figure;
    %       f.fig.Renderer = 'painters';
    %       ax=gca;
    %       [ax,torispl,mamp,xbnd,ybnd] = plt_decon_imp_scatter(ax,impindepst,xran,yran,cran,offxran,offyran,...
    %         sps,50,'mean','tarvl');
    %       scatter(ax,off1i(2),off1i(3),20,'ks','filled','MarkerEdgeColor','k');
    %       title(ax,'Checkd at 4th stas');
    %
    %       %plot the scatter of sources in terms of rela locations
    %       xran = [-4 4];
    %       yran = [-4 4];
    %       cran = [0 lsig/sps];
    %       f.fig = figure;
    %       f.fig.Renderer = 'painters';
    %       ax=gca;
    %       [ax] = plt_decon_imp_scatter_space(ax,impindepst,xran,yran,cran,offxran,...
    %         offyran,sps,50,ftrans,'mean','tarvl');
    %       text(ax,0.98,0.95,sprintf('Satur=%.1f',sat),'Units','normalized',...
    %         'HorizontalAlignment','right');
    %       text(ax,0.02,0.05,sprintf('%.2f',mprojx2all4th),'Units','normalized',...
    %         'HorizontalAlignment','left');
    %       plot(ax,xcut,ycut,'k-','linew',2);
    %       title(ax,'Checkd at 4th stas');
    
    %plot the cumulative density and summed amp of detections
    [density1d, inddup] = density_pixel(impindepst(:,7),impindepst(:,8)); %count at each unique loc
    [imploc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
    density1d = [imploc(:,1:2) density1d(:,3)];
    amp = mean(impindepst(:,[2 4 6]),2); %amp for all LFE catalog
    ampsum = sum_at_indices(amp,inddup);
    ampsum1d = sortrows([imploc(:,1:2) ampsum], 3);
    cstr = {'# detections / pixel'; 'amp sum / pixel'};
    [f] = plt_sum_pixel(density1d,ampsum1d,[-4 4],[-4 4],20,cstr,'o','linear');
    hold(f.ax(1),'on');
    plot(f.ax(1),xcut,ycut,'k-','linew',2);
    hold(f.ax(2),'on');
    plot(f.ax(2),xcut,ycut,'k-','linew',2);
    supertit(f.ax,'Checkd at 4th stas');
    
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
    f = initfig(8,5,1,2); %initialize fig
    ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    histogram(ax,pred4difftr);
    plot(ax,[median(pred4difftr) median(pred4difftr)],ax.YLim,'r--','LineWidth',1);
    plot(ax,[-offmax -offmax],ax.YLim,'k--');
    plot(ax,[offmax offmax],ax.YLim,'k--');
    xlabel(ax,'diff. in 4th arrival between pred and decon');
    ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    histogram(ax,diffoff14tr);
    plot(ax,[median(diffoff14tr) median(diffoff14tr)],ax.YLim,'r--','LineWidth',1);
    plot(ax,[-offmax -offmax],ax.YLim,'k--');
    plot(ax,[offmax offmax],ax.YLim,'k--');
    xlabel(ax,'diff. in off14 between plane-fit and decon');
    ylabel(ax,'count');
    
    %       %plot distance between srcs N&N-m along the projected direction
    %       if ~isempty(dprojxy)
    %         nsep = 1;
    %         [f] = plt_customprojdist(tarvlsplst,implocst,dtarvl{nsep},eucdist{nsep},...
    %           projxy,dprojxy,projang,sps,'tarvl');
    %         hold(f.ax(1),'on');
    %         plot(f.ax(1),xcut,ycut,'k-','linew',2);
    %       end
    
    %       %plot histograms of source amp
    %       f = initfig(16,5,1,4); %initialize fig
    %       f.ax(1:4) = plt_deconpk_rat_comb4th(f.ax(1:4),srcampr4th,impindepst,'k','hist');
    %       supertit(f.ax,'Checkd at 4th stas');
  end
  
  %combine results from same sat and reg, but diff runs into same array
  irccranall = [irccranall; irccran+lsig*(irun-1)];
  windowsall = [windowsall; windows-windows(1,1)+1+lsig*(irun-1)];
  off1iwall = [off1iwall; off1iw];
  ccaliwall = [ccaliwall; ccaliw];
  ccwpairall = [ccwpairall; ccwpair];
  mrccwpairall = [mrccwpairall; mrccwpair];
  rcccatall = [rcccatall; rcccat];
  off1iall = [off1iall; off1i];
  rccall = [rccall; rcc];
  mrccall = [mrccall; mrcc];
  ccall = [ccall; cc];
  impgtall = [impgtall; impgt];
  nitall = [nitall; nit];
  impgrpall = [impgrpall; impgrp];
  nsrcgrpall = [nsrcgrpall; size(impgrp,1)];
  impall = [impall; imp];
  nsrcall = [nsrcall; size(imp,1)];
  mprojxnn1all = [mprojxnn1all; mprojxnn1];
  mprojx2allall = [mprojx2allall; mprojx2all];
  srcamprall = [srcamprall; srcampr];
  msrcamprall = [msrcamprall; msrcampr];
  madsrcamprall = [madsrcamprall; madsrcampr];
  mpsrcamprsall = [mpsrcamprsall; mpsrcamprs];
  madpsrcamprsall = [madpsrcamprsall; madpsrcamprs];
  mnsrcamprsall = [mnsrcamprsall; mnsrcamprs];
  madnsrcamprsall = [madnsrcamprsall; madnsrcamprs];
  imp4thall = [imp4thall; imp4th];
  nsrc4thall = [nsrc4thall; size(imp4th,1)];
  diffoff14trall = [diffoff14trall; diffoff14tr];
  mprojxnn14thall = [mprojxnn14thall; mprojxnn14th];
  mprojx2all4thall = [mprojx2all4thall; mprojx2all4th];
  srcampr4thall = [srcampr4thall; srcampr4th];
  msrcampr4thall = [msrcampr4thall; msrcampr4th];
  madsrcampr4thall = [madsrcampr4thall; madsrcampr4th];
  mpsrcamprs4thall = [mpsrcamprs4thall; mpsrcamprs4th];
  madpsrcamprs4thall = [madpsrcamprs4thall; madpsrcamprs4th];
  mnsrcamprs4thall = [mnsrcamprs4thall; mnsrcamprs4th];
  madnsrcamprs4thall = [madnsrcamprs4thall; madnsrcamprs4th];
  
  %   keyboard
  
end  %loop end for number of runs

%% Ouput everything in the form of a structure array
allsyn.irccran = irccranall;
allsyn.windows = windowsall;
allsyn.off1iw = off1iwall;
allsyn.ccaliw = ccaliwall;
allsyn.ccwpair = ccwpairall;
allsyn.mrccwpair = mrccwpairall;
allsyn.rcccat = rcccatall;
allsyn.off1i = off1iall;
allsyn.rcc = rccall;
allsyn.mrcc = mrccall;
allsyn.cc = ccall;
allsyn.impgt = impgtall;
allsyn.nit = nitall;
allsyn.nitsum = sum(nitall,1); %1 value; all runs are summed or concatenated
allsyn.impgrp = impgrpall;  %concatenation of all runs
allsyn.nsrcgrp = nsrcgrpall;
allsyn.nsrcgrpsum = sum(nsrcgrpall); %1 value; all runs are summed or concatenated
allsyn.imp = impall;
allsyn.nsrc = nsrcall;
allsyn.nsrcsum = sum(nsrcall); %1 value; all runs are summed or concatenated
allsyn.mprojxnn1 = median(mprojxnn1all);  %1 value
allsyn.mprojx2all = median(mprojx2allall);  %1 value
allsyn.srcampr = srcamprall;
allsyn.msrcampr = median(msrcamprall);  %1 value; a median of all runs, suppose each run has a similar value
allsyn.madsrcampr = median(madsrcamprall);  %1 value
allsyn.mpsrcamprs = median(mpsrcamprsall);  %1 value
allsyn.madpsrcamprs = median(madpsrcamprsall);  %1 value
allsyn.mnsrcamprs = median(mnsrcamprsall);  %1 value
allsyn.madnsrcamprs = median(madnsrcamprsall);  %1 value
allsyn.imp4th = imp4thall;
allsyn.nsrc4th = nsrc4thall;
allsyn.nsrc4thsum = sum(nsrc4thall);  %1 value; all runs are summed or concatenated
allsyn.diffoff14 = diffoff14trall;
allsyn.mprojxnn14th = median(mprojxnn14thall);  %1 value; a median of all runs, suppose each run has a similar value
allsyn.mprojx2all4th = median(mprojx2all4thall);  %1 value
allsyn.srcampr4th = srcampr4thall;
allsyn.msrcampr4th = median(msrcampr4thall);  %1 value
allsyn.madsrcampr4th = median(madsrcampr4thall);  %1 value
allsyn.mpsrcamprs4th = median(mpsrcamprs4thall);  %1 value
allsyn.madpsrcamprs4th = median(madpsrcamprs4thall);  %1 value
allsyn.mnsrcamprs4th = median(mnsrcamprs4thall);  %1 value
allsyn.madnsrcamprs4th = median(madnsrcamprs4thall);  %1 value

%% save file
savefile = strcat('rst_decon_synth_mix','_reg',num2str(xaxis),'-',num2str(yaxis),...
  '_noi',num2str(perc),'_nsat',num2str(sat),'_td',num2str(tdura),'.mat');
save(strcat(workpath,'/synthetics/',savefile), 'allsyn');







