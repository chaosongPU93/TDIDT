function synthshift_chao_fn(srcregion,xaxis,yaxis,nsat,distr,Twin,tdura,seed,distrloc,...
  timetype,physicalsize,forcesep,testsrcflag,normflag,tempflag,ftrans,pltflag,irun)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the functionalized version of 'synthshift_chao.m' so that you can
% customize the parameters in generating the synthetics from different 
% region size, noise level, simulation length, etc.
%
% It is now able to do a lot of different things by turning on/off
% different flags, for example, with uniformly and randomly 
% distributed sources with/without a physical dimension from an 
% specified size of region, with different saturation level in time.
% It could also draw samples from any custom PDF, not necessarily the
% uniform distribution. (DONE as of 2023/09/07)
% You can also choose different transformation time offset grid 
% between offset and spatial location, either Allan's grid from direct
% interpolation from JA's inverted but sparser grid, or my own 
% inversion. (DONE as of 2023/09/07)
% Rather than having sources from a limited region with a negligible
% noise level, another extreme is to put all sources at a single spot
% with variable noise level, and see how the result like, compared 
% to statistics from data. (subject to change as of 2023/09/07)
%
% Available flags are:
%   distrloc -- distribution for source location
%   timetype -- which time is uniform in time
%   physicalsize -- if considering the physical size of each source
%   srcregion -- shape of the source region
%   ftrans -- regime for location transformation
%   forcesep -- if forcing a min speration of arrival time for 2 events
%               from the same spot
%
% The algorithm of addition of ground truth sources to seismograms is
% confirmed to reproducable by using convolution between source impulses
% and templates. 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/05/02
% Last modified date:   2024/05/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% default settings
defval('srcregion','ellipse');  %specify shape of the source region
defval('xaxis',1.75); %axis length of the same ellipse of my 4-s catalog
defval('yaxis',1.25);
defval('nsat',[0.1 0.4 1 2 4 10 20 40 100]);  %saturation level
defval('distr','UN'); %specify the amplitude-frequency (counts) distribution as uniform
defval('Twin',0.5*3600+3+2*ceil(2048/160)); %length of each simulation
defval('tdura',0.25); %duration of templates
defval('seed',3); %for random number generator
defval('distrloc','uniform'); %uniformly random for source location
defval('timetype','tori');  %specify which time is used for source
defval('physicalsize',0); %specify if considering the physical size of each source
defval('forcesep', 0);  %specify if forcing a min speration of arrival time for 2 events from the same spot 
defval('testsrcflag',0);  %whether to test if synthetics can be reproduced by convolution
defval('normflag',0); %whether to normalize templates
defval('tempflag','chao'); %which templates to use
defval('ftrans','interpchao');  %specify regime for transformation from time offset to map location 
defval('pltflag',0);  %whether to show a plot of raw synthetics
defval('irun',1);  %which run out of the total # of runs of simulations

if ~isempty(irun)
  irunstr = num2str(irun);
else
  irunstr = [];
end

%% Initialization
format short e   % Set the format to 5-digit floating point
% clear
% clc
% close all

% set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');

%% prepare templates (Green's functions), from 'lfetemp002_160sps.m' or Allan's templates
% normflag = 0;

% tempflag = 'allan';
% tempflag = 'chao';

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

%%%plot the unfiltered and filtered templates
% plt_templates(green,greenf,stas,[],[],lowlet,hiwlet,sps);

%just the filtered templates
% plt_templates_bp(greenf,stas,lowlet,hiwlet,sps);

% plt_templates(green(:,1:3),greenf(:,1:3),stas(1:3,:),[],[],lowlet,hiwlet,sps);
% zcrosses
% 2*(ppeaks-npeaks)/sps
% 
% zcrossesf
% 2*(ppeaksf-npeaksf)/sps
% keyboard

%% synthetics generation
%%%Specify the amplitude-frequency (counts) distribution 
% distr='UN'  % uniform distribution

% Twin=0.5*3600+3+2*ceil(greenlen/sps); %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
winlen=Twin*sps+1;
skiplen=greenlen;
%%%for the duration of templates, there are several options
%%%1. (ppeak-npeak)*2 of the BB template: 44;54;38;38 0.275s, 0.3375s, 0.2375s, 0.2375s
%%%2. direct eyeballing for between zerocrossings: ~65, 0.4s
%%%3. binned peak-to-peak separation for decently saturated unfiltered synthetics: ~37, 0.25s
%%%4. similar to 3, but synthetics are filtered first: ~37, 0.25s
% tdura = 0.4;  % duration from Chao's broadband template, width is about 795-730=65 spls at 160 hz
% tdura = 0.25; % start to use on 2023/12/07
satn=1/tdura*Twin   % if just saturated, how many templates can be fit in? a single peak is ~20 samples wide; maybe a little less (at 100 sps).
%Twin is window duration in seconds. Events can fall within Twin of
%the start, but the synthetics will go to Twin*(sample rate)+Greenlen to
%avoid checking for subscript overrrun.  When "synths" are written from "synth" in the
%subroutine, Greenlen from the start and end will not be written.
fracelsew=0; %0.25 %0.5; %0.6; %The ratio of elsewhere events to local events.  Zero for Discrete Ide.
% nsat=[0.1 0.4 1 2 4 10 20 40 100];  % times of saturation 
% nsat=[95 100 200];  % times of saturation 
nnsat = length(nsat);
writes=round(nsat*satn) %how many templates to throw in, under different degrees of saturation

%%%Here the noise is 'uniform' in time!
% noistd = 5e-2;
% noistd = 2.e-7;
% noistd = 2.0e-4;
noistd = 0; %NOISE-FREE
rng('default');
% seed=(irun-1)*nreg+ireg;
% rng(seed);
synth=noistd*(randn(winlen+greenlen+2*10,nsta)-0.5); %+2*10 a little extra, for jiggering 2nd & 3rd stations.% for ista=1:nsta

nouts=length(writes);
% seed=round(writes(4)/5e3); %for random number generator
% seed=10;

%%%specify distribution for source location  
% distrloc = 'custompdf'; %using a custom PDF function
% distrloc = 'uniform'; %uniformly random in a specified region,

%%%specify which time is uniform in time
% timetype = 'tarvl';
% timetype = 'tori';

%%%specify if considering the physical size of each source
% physicalsize = 1;
% physicalsize = 0;

%%%specify shape of the source region
% srcregion='ellipse';
% srcregion='rectangle';
% srcregion='circle';

%%%specify regime for transformation from time offset to map location 
% ftrans = 'interpArmb';
% ftrans = 'interpArmbreloc';
% ftrans = 'interpchao';

%%%specify if forcing a min speration of arrival time for 2 events from the same spot 
% forcesep = 1;
% forcesep = 0;

if strcmp(distrloc, 'custompdf')
  
  %% load the true 4-s tremor detections, obtain its PDF that will be used to generate synthetic sources
  freqflag='hf';  % flag to indicate whether to do hf or lf;
  FLAG = 'PGC'; % detector 
  fam = '002';   % family number
  spsdect = 40;
  iup = sps/spsdect;  % upsample 4 times to 160 sps
  % load detections
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
  wlen=winlensec*spsdect;      % length in smaples
  PREFIX = strcat(fam,'.up.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),'.',...
                          int2str(npo),int2str(npa),'.ms', int2str(mshift));
  hftime = load(strcat(rstpath, '/MAPS/tdectimeori_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',...
              num2str(wlen/spsdect),'s',num2str(spsdect),'sps','4add'));
  %format as follows:
  %%% 34+4*nstanew cols, if 4 new stas, then will be 50 cols
  %%% UPDATED at 2021/06/23
  %%%   1:off12(n) 2:off13(n) 3:off12sec 4:off13sec 5:fam 6:date 7:timswin(n)
  %%%   8:timswin(n)-winlensec/2+idiff/sps 9:xcmaxAVEnbang(nin) 10:loopoff(n)
  %%%   11:cumsumtrdiff  12:cumsumtrdiff/cumsumtr(wlen)
  %%%   for the rest, n=n+4;
  %%%   9:Ampsq(nin) 10:Prior(nin) 11:Post(nin) 12:Prior4(nin) 13:Post4(nin) 14:Prior8(nin)
  %%%   15:Post8(nin) 16:Prior16(nin) 17:Post16(nin) 18:cc(nin) 19:ccprior(nin)
  %%%   20:ccpost(nin) 21:ccprior4(nin) 22:ccpost4(nin) 23:STAamp(nin,2) 24:STAamp(nin,3)
  %%%   25:sumsSTA12n(n,iSTA12bang) 26:sumsSTA13n(n,iSTA13bang)
  %%%   27:sumsSTA32n(n,iSTA32bang) 28:sigma(nin) 29:sigma12(nin) 30:sigma13(nin)
  %%%   in(1:nstanew) loff(1:nstanew) ioff(1:nstanew) ccmaxave(1:nstanew)

  %%%Interpolate from existing grids for the locations of sources
  %round to the integer sample at the 'iup' times the original sampling rate 40 sps
  offint = round(hftime(:, 1:2)*iup);

  %convert the time offset to locations
  fplt = 0;
  ind = find(offint(:,1)>=-23*iup & offint(:,1)<=25*iup & ...
    offint(:,2)>=-24*iup & offint(:,2)<=24*iup);
  [loc, indinput] = off2space002(offint(ind,:),sps,ftrans,fplt);
  % loc has 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
  hfuse = hftime(ind,:);
  hfuse = [loc hfuse];

  %outline a boundary region on top of the density map
  cutout = 'ellipse';
  x0 = 0.2;
  y0 = 0.2;
  semia = 1.75;
  semib = 1.0;
  angrot = 45;
  [x, y] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);

  %get the detections inside the cutout boundary
  bnd = [x y];
  [iin,ion] = inpolygon(hfuse(:,1),hfuse(:,2),bnd(:,1),bnd(:,2));
  isinbnd = iin | ion;
  hfbnd = hfuse(isinbnd == 1, :);
  %format as follows:
  %%% 8+34+4*nstanew cols, if 4 new stas, then will be 58 cols
  
  %%%Obtain the location density (count), bin by pixel
  % dx = 0.05;
  % dy = 0.05;
  offxran = [-4 4];
  offyran = [-4 4];
  denloc= density_pixel(hfbnd(:,1),hfbnd(:,2));
  %plot the density distribution
  figure
  scatter(denloc(:,1),denloc(:,2),6,denloc(:,3),'filled');
  axis equal
  axis([offxran offyran]);
  colormap(jet)
  c=colorbar;
  c.Label.String = '# tremor detections / pixel';
  caxis([0 max(denloc(:,3))]);
  xlabel('E (km)');
  ylabel('N (km)');
  box on; grid on;

  %%%Note that the locations in the ellipse corresponds to integer time offsets at sps*iup, so we can
  %%%avoid doing the transformation, and directly get the density and PDF in sample space, the benefit
  %%%is that, in sample sapce, it is actucally an even grid with a spacing of 1 sample 
  denoff= density_pixel(hfbnd(:,7),hfbnd(:,8));

  %%%plot the density distribution
  % figure
  % scatter(denoff(:,1),denoff(:,2),30,denoff(:,3),'filled');
  % axis equal
  % off12ran = [min(denoff(:,1))-1 max(denoff(:,1))+1];
  % off13ran = [min(denoff(:,2))-1 max(denoff(:,2))+1];
  % axis([off12ran off13ran]);
  % colormap(jet)
  % c=colorbar;
  % c.Label.String = '# tremor detections / pixel';
  % caxis([0 max(denoff(:,3))]);
  % xlabel(strcat({'Offset 12 (samples at '},num2str(sps),' sps)'));
  % ylabel(strcat({'Offset 13 (samples at '},num2str(sps),' sps)'));
  % box on; grid on;

  %Obtain the PDF
  epdfoff = denoff;
  %now the 'area' of each bin is 1*1 sample
  epdfoff(:,3) = epdfoff(:,3)/sum(epdfoff(:,3))/(1*1); % normalize to PDF, ==counts/N_total/area_of_bin

  %%%plot the PDF, note the color should be the SAME as density, as it is just a normalization
  % figure
  % scatter(epdfoff(:,1),epdfoff(:,2),30,epdfoff(:,3),'filled');
  % axis equal
  % axis([off12ran off13ran]);
  % colormap(jet)
  % c=colorbar;
  % c.Label.String = 'Probability Density';
  % caxis([0 max(epdfoff(:,3))]);
  % xlabel(strcat({'Offset 12 (samples at '},num2str(sps),' sps)'));
  % ylabel(strcat({'Offset 13 (samples at '},num2str(sps),' sps)'));
  % box on; grid on;

  %%% Note that 'pinky(x1,x2,y)' can only deal with an evenly-spaced vector x1 and x2, y is a grid of
  %%% PDF values on the grid points that are defined by x1 and x2.
  %%% Therefore, 'pinky(x1,x2,y)' can be applied when either your PDF is generated by binning upon an 
  %%% evenly-spaced grid (or ksdensity), OR, you bin by pixel in locations but return to the sample
  %%% domain, in which case the grid is even. But the latter case needs zero-padding to grid points
  %%% that is otherwise empty.
  xvec = -max(abs(epdfoff(:,1)))-1: 1: max(abs(epdfoff(:,1)))+1;
  yvec = -max(abs(epdfoff(:,2)))-1: 1: max(abs(epdfoff(:,2)))+1;
  [epdfoffpad,xgrid,ygrid,epdfoffgrid,ind1] = zeropadmat2d(epdfoff,xvec,yvec);

  %%%plot the PDF, note the color should be the SAME as density, as it is just a normalization
  figure
  dum = epdfoffpad;
  dum(dum(:,3)~=0, :) = [];
  scatter(dum(:,1),dum(:,2),6,dum(:,3),'linew',0.2);  hold on
  dum = epdfoffpad;
  dum(dum(:,3)==0, :) = [];
  scatter(dum(:,1),dum(:,2),30,dum(:,3),'filled');
  axis equal
  axis([minmax(xvec) minmax(yvec)]);
  colormap(jet)
  c=colorbar;
  c.Label.String = 'Probability Density';
  caxis([0 max(epdfoffpad(:,3))]);
  title('Zero-padded PDF for location sampling');
  xlabel(strcat({'Offset 12 (samples at '},num2str(sps),' sps)'));
  ylabel(strcat({'Offset 13 (samples at '},num2str(sps),' sps)'));
  box on; grid on;

  %save the source location grid
  fid = fopen([workpath,'/synthetics/synsrcloc.custompdf.grd'],'w');
  fprintf(fid,'%8.3f %8.3f %8.3f\n',epdfoffpad');
  fclose(fid);
  size(epdfoffpad)

  %%%generation of synthetics 
  b=999. %>150 for uniform size distribution
  %USE unfiltered templates to generate synthetics, then bandpass before deconvolution
  [synths,mommax,sources,greensts]=synthgen(writes,winlen,skiplen,synth,green,b,xgrid,...
    ygrid,epdfoffgrid,sps,fracelsew,seed,timetype,ftrans,stas);

  % %simple check if the tori or tarvl is indeed uniform in time
  % figure
  % inwrites = 5;
  % n=writes(inwrites);
  % a = sources(1:n,:);
  % b = a(any(a,2),:);
  % if strcmp(timetype,'tarvl')
  %   histogram(b(:,1),'facecolor','k');
  %   xlabel('Arrival time');
  % elseif strcmp(timetype,'tori')
  %   histogram(b(:,1),'facecolor','k'); hold on;
  %   histogram(b(:,4),'facecolor','r');
  %   xlabel('Arrival/Origin time');
  % end
  % ylabel('Count');
  % keyboard
  
elseif strcmp(distrloc,'uniform')
  
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

  %whether to consider the physical size of each source
  if physicalsize  
    %%%%%%%%%%%%% start, main revision from 'synthshift.m' %%%%%%%%%%%%%%%%%%%%%%
    offPGSS=xygrid(:,1);
    offPGSI=xygrid(:,2);
    dx=xygrid(:,3);
    dy=xygrid(:,4);
    % Make a hex grid
    llx=xygrid(1,3)+0.1; %+1 to make the grid more symmetric for this choice of diam, a, and b.
    lly=xygrid(1,4);  %lower left corner x and y
    urx=xygrid(end,3);  %upper right corner x and y    
    ury=xygrid(end,4);
    diam=0.3;  %note that if diam is <150m, in sample space you'll have too many nonunqiue sources 
    [hxq,hyq]=meshgrid(llx:diam:urx, lly:diam*sqrt(3)/2:ury); %area of a hexagon is sqrt(3)*a^2, a=diam/2
    for i=2:2:size(hxq,1)
        hxq(i,:)=hxq(i,:)+0.5*diam; %shift x every 2 rows, to create hex grid
    end
    angrot=-15*pi/180;  %rotate hex grid clockwise by 15 degs
    hexloc=complex(hxq,hyq);
    hexrot=hexloc*exp(1i*angrot);
    hxqrot=real(hexrot);
    hyqrot=imag(hexrot);
    hexPGSS=griddata(dx,dy,offPGSS,hxqrot,hyqrot);
    hexPGSI=griddata(dx,dy,offPGSI,hxqrot,hyqrot);
%     hexPGSS=griddata(dx,dy,offPGSS,hxq,hyq);
%     hexPGSI=griddata(dx,dy,offPGSI,hxq,hyq);
    len=size(hxq,1)*size(hxq,2);
    hexPGSSvect=reshape(hexPGSS,len,1);
    hexPGSIvect=reshape(hexPGSI,len,1);
    hxqrotvect=reshape(hxqrot,len,1);
    hyqrotvect=reshape(hyqrot,len,1);
    hexgrid=[hexPGSSvect, hexPGSIvect, hxqrotvect, hyqrotvect];
    %%%%%%%%%%%%% end, main revision from 'synthshift.m' %%%%%%%%%%%%%%%%%%%%%%
    
    xygrid = hexgrid;
  else
    diam=0;
  end
  
  %which shape of source region to use
  if strcmp(srcregion,'rectangle') %The following rotates 45˚ and looks for -5 < x' < 1 and 2.25 < y' < 3.5
    angrot=-45*pi/180;
    loc=complex(xygrid(:,3),xygrid(:,4));
    locrot=loc*exp(1i*angrot);
    xygrid(:,5)=real(locrot);
    xygrid(:,6)=imag(locrot);
    ll=[-5 2.25]; ur=[1 3.5];
    xygrid(xygrid(:,5)<ll(1) | xygrid(:,5)>ur(1) | xygrid(:,6)<ll(2) | xygrid(:,6)>ur(2),:)=[]; %This limits xygrid to the desired range.  PGSS, PGSI, x, y, x', y'
  elseif strcmp(srcregion,'ellipse') %The following adds 0.5 km to y, rotates 45˚, and looks for an elliptical region with major semi-axis y' < 3 km and minor semi-axis < 2 km
    angrot=-45*pi/180;
%     shiftor=[0 -0.5]; %(in km)  %center used by Allan
    shiftor=[0.2 0.2]; %(in km) %center of the same ellipse of my 4-s catalog
    loc=complex(xygrid(:,3)-shiftor(1),xygrid(:,4)-shiftor(2));
    locrot=loc*exp(1i*angrot);
    xygrid(:,5)=real(locrot);
    xygrid(:,6)=imag(locrot);
%     xaxis=0.5; %(semi-, in km)
%     yaxis=0.5; %(semi-, in km)
%     xaxis=3.0; %(semi-, in km)
%     yaxis=1.5; %(semi-, in km)
%     xaxis=2.75; %(semi-, in km)
%     yaxis=1.0; %(semi-, in km)
%     xaxis=3.0; %(semi-, in km)
%     yaxis=1.25; %(semi-, in km)
%     xaxis=2.0; %(semi-, in km)
%     yaxis=1.25; %(semi-, in km)
    %variation of source region size
%     semia = 1.75*(0.6:0.2:2.0);
%     semib = 1.25*(0.6:0.2:2.0);
%     nreg = length(semia);
%     ireg = 3;
%     xaxis = semia(ireg); %axis length of the same ellipse of my 4-s catalog
%     yaxis = semib(ireg);
    % xaxis=1.75; %axis length of the same ellipse of my 4-s catalog
    % yaxis=1.25;
    [xcut,ycut] = ellipse_chao(shiftor(1),shiftor(2),xaxis,yaxis,0.01,45,shiftor);
    xygrid(sqrt((xygrid(:,5)/xaxis).^2 + (xygrid(:,6)/yaxis).^2) > 1,:)=[]; %This limits xygrid to the desired range.  PGSS, PGSI, x, y, x', y'
  elseif strcmp(srcregion,'circle') %The following adds 0.5 km to y, rotates 45˚, and looks for an elliptical region with major semi-axis y' < 3 km and minor semi-axis < 2 km
    angrot=-45*pi/180;
%     shiftor=[0 -0.5]; %(in km)  %center used by Allan
    shiftor=[0.2 0.2]; %(in km) %center of the same ellipse of my 4-s catalog
    loc=complex(xygrid(:,3)-shiftor(1),xygrid(:,4)-shiftor(2));
    locrot=loc*exp(1i*angrot);
    xygrid(:,5)=real(locrot);
    xygrid(:,6)=imag(locrot);
    radi=1.25; %radius
    [xcut,ycut] = circle_chao(shiftor(1),shiftor(2),radi,0.01);
    xygrid(sqrt((xygrid(:,5)/radi).^2 + (xygrid(:,6)/radi).^2) > 1,:)=[]; %This limits xygrid to the desired range.  PGSS, PGSI, x, y, x', y' 
  end
  size(xygrid)

  if pltflag
    %whether to consider the physical size of each source
    figure
    subplot(1,2,1)
    axis equal
    hold on; grid on; box on
    plot(xygrid(:,1)*iup,xygrid(:,2)*iup,'.');
    xlabel(sprintf('off12 samples (%d Hz)',sps));
    ylabel(sprintf('off13 samples (%d Hz)',sps));
    tmp = ceil(max(max(abs([xygrid(:,1)*iup,xygrid(:,2)*iup]))))+1;
    xran = [-tmp tmp];
    yran = [-tmp tmp];
    xlim(xran);
    ylim(yran);
    subplot(1,2,2)
    axis equal
    hold on; grid on; box on
    scatter(shiftor(1),shiftor(2),15,'k','filled');
    plot(xcut,ycut,'k-','linew',1);
    plot(xygrid(:,3),xygrid(:,4),'.');
    text(0.95,0.05,sprintf('%d unique sources',size(xygrid,1)),'Units','normalized',...
      'HorizontalAlignment','right');
    tmp = ceil(max(max(abs([xcut,ycut]))));
    xran = [-tmp tmp];
    yran = [-tmp tmp];
    xlim(xran);
    ylim(yran);
    xlabel('E (km)');
    ylabel('N (km)');
    if physicalsize
      %SEE Chao's NOTES: research/geometry of synthetics
      rad=0.5*diam; %radius of blue circles, max. with no overlapping
      areafrac=0.5; %fraction of area of single asperity
      lod2=areafrac*2*sqrt(3)/pi; %lod2 is ("asperity diameter"/diam)^2, where diam is hex grid spacing
      asprad=0.5*sqrt(lod2)*diam; %asperity radius
      ang=0:0.1*pi:2*pi;
      xdiam=rad*cos(ang);
      ydiam=rad*sin(ang);
      xasp=asprad*cos(ang);
      yasp=asprad*sin(ang);
      for ispot=1:size(xygrid,1)
        plot(xygrid(ispot,3)+xdiam,xygrid(ispot,4)+ydiam,'b--');
        plot(xygrid(ispot,3)+xasp,xygrid(ispot,4)+yasp,'r');
      end
    end
  end

%   keyboard
  %whether to force a min speration if 2 events from the same spot separated by less than the template duration
  %USE broadband templates to generate synthetics, then bandpass before deconvolution
  if forcesep
    [synths,mommax,sources,greensts]=csplaw3d(writes,winlen,skiplen,synth,green,b,...
      xygrid,sps,fracelsew,seed,tdura,timetype,ftrans,stas); 
  else
    [synths,mommax,sources,greensts]=csplaw3c(writes,winlen,skiplen,synth,green,b,...
      xygrid,sps,fracelsew,seed,timetype,ftrans,stas);
  end
  
  % %simple check if the tori or tarvl is indeed uniform in time
  % figure
  % inwrites = 5;
  % n=writes(inwrites);
  % if forcesep
  %   a = squeeze(sources(1:n,:,inwrites));
  %   b = a(any(a,2),:);
  % else
  %   a = sources(1:n,:);
  %   b = a(any(a,2),:);
  % end
  % if strcmp(timetype,'tarvl')
  %   histogram(b(:,1),'facecolor','k');
  %   xlabel('Arrival time');
  % elseif strcmp(timetype,'tori')
  %   histogram(b(:,1),'facecolor','k'); hold on;
  %   histogram(b(:,3),'facecolor','r');
  %   xlabel('Arrival/Origin time');
  % end
  % ylabel('Count');
  % keyboard

end

%% a simple plot the synthetics
if pltflag
  % figure
  % tmp = synths(:,:,2)-synth(skiplen:winlen,:);
  % plot(tmp(:,1),'b');
  % xlim([0 31*sps]);

  figure
  nrow = length(writes)+1;
  ncol = 1;
  subplot(nrow,ncol,1)
  hold on
  tmp = synth;
  if nsta == 3
    plot(tmp(:,1),'r');
    plot(tmp(:,2),'b');
    plot(tmp(:,3),'k');
  else
    color = jet(nsta);
    for ista = 1: nsta
      plot(tmp(:,ista),'Color',color(ista,:));
    end
  end
  xlim([0 40*sps]);
  axranexp(gca,6,20);

  for i = 1: length(writes)
  subplot(nrow,ncol,i+1)
  hold on
  tmp = synths(:,:,i);
  if nsta == 3
    plot(tmp(:,1),'r');
    plot(tmp(:,2),'b');
    plot(tmp(:,3),'k');
  else
    color = jet(nsta);
    for ista = 1: nsta
      plot(tmp(:,ista),'Color',color(ista,:));
    end
  end
  text(0.95,0.9,sprintf('%.1f x saturation',nsat(i)),'Units','normalized','HorizontalAlignment',...
    'right');
  xlim([0 40*sps]);
  axranexp(gca,6,20);
  end
  xlabel(sprintf('Samples at %d Hz',sps),'FontSize',12);

end

% keyboard

%% save files
for inwrites=1:length(writes)
  n=writes(inwrites);
  %%%seismograms
  if strcmp(srcregion,'ellipse')
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(xaxis),'-',num2str(yaxis),'.diam',num2str(diam),'.else',num2str(fracelsew,2),...
      'nsat',num2str(nsat(inwrites)),'tdura',num2str(tdura),...
      'T',num2str(round(Twin)),'p',irunstr],'w');
  elseif strcmp(srcregion,'rectangle')
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(ll(1)),num2str(ll(2)),num2str(ur(1)),num2str(ur(2)),'.diam',num2str(diam),...
      '.else',num2str(fracelsew,2),'nsat',num2str(nsat(inwrites)),'tdura',...
      num2str(tdura),'T',num2str(round(Twin)),'p',irunstr],'w');
  elseif strcmp(srcregion,'circle')
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(radi),'.diam',num2str(diam),'.else',num2str(fracelsew,2),...
      'nsat',num2str(nsat(inwrites)),'tdura',num2str(tdura),...
      'T',num2str(round(Twin)),'p',irunstr],'w');  
  end
  towrite=squeeze(synths(:,:,inwrites));
  if nsta == 3
    fprintf(fid,'%13.5e %13.5e %13.5e \n', towrite');
  elseif nsta == 4
    fprintf(fid,'%13.5e %13.5e %13.5e %13.5e \n', towrite');
  end
  fclose(fid);
  
  %%%source info
  if strcmp(srcregion,'ellipse')
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(xaxis),'-',num2str(yaxis),'.diam',num2str(diam),'.else',num2str(fracelsew,2),...
      'nsat',num2str(nsat(inwrites)),'tdura',num2str(tdura),...
      'T',num2str(round(Twin)),'p',irunstr,'_sources'],'w');
  elseif strcmp(srcregion,'rectangle')
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(ll(1)),num2str(ll,2),num2str(ur(1)),num2str(ur,2),'.diam',num2str(diam),...
      '.else',num2str(fracelsew,2),'nsat',num2str(nsat(inwrites)),'tdura',...
      num2str(tdura),'T',num2str(round(Twin)),'p',irunstr,'_sources'],'w');
  elseif strcmp(srcregion,'circle')
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(radi),'.diam',num2str(diam),'.else',num2str(fracelsew,2),...
      'nsat',num2str(nsat(inwrites)),'tdura',num2str(tdura),...
      'T',num2str(round(Twin)),'p',irunstr,'_sources'],'w');  
  end
  if strcmp(distrloc,'uniform')
%     if forcesep
      a = squeeze(sources(1:n,:,inwrites));
      b = a(any(a,2),:);
%     else
%       a = sources(1:n,:);
%       b = a(any(a,2),:);
%     end
    size(b)
    truncsources=b;
    if strcmp(timetype,'tarvl')
      fprintf(fid,'%13i %5i \n', truncsources'); %arrival; loc ind of grid
    elseif strcmp(timetype,'tori')
      fprintf(fid,'%13i %5i %13i %13i \n', truncsources'); %arrival; loc ind of grid; origin; travel time
    end
    fclose(fid);
  elseif strcmp(distrloc,'custompdf')
    a = sources(1:n,:);
    b = a(any(a,2),:);
    truncsources=b;
    if strcmp(timetype,'tarvl')
      fprintf(fid,'%d %d %d %.1f \n', truncsources'); %arrival; off12; off13; moment
    elseif strcmp(timetype,'tori')
      fprintf(fid,'%d %d %d %d %d %.1f \n', truncsources'); %arrival; off12; off13; origin; travel time; moment
    end   
    fclose(fid);
  end
  
  %save the source grid location 
%   fid = fopen([workpath,'/synthetics/synsrcloc.',srcregion(1:3),'.grd'],'w');
  if strcmp(srcregion,'ellipse')
    fid = fopen([workpath,'/synthetics/synsrcloc.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(xaxis),'-',num2str(yaxis),'.diam',num2str(diam),...
      'T',num2str(round(Twin)),'p',irunstr,'_grd'],'w');
  elseif strcmp(srcregion,'rectangle')
    fid = fopen([workpath,'/synthetics/synsrcloc.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(ll(1)),num2str(ll,2),num2str(ur(1)),num2str(ur,2),'.diam',...
      num2str(diam),'T',num2str(round(Twin)),'p',irunstr,'_grd'],'w');
  elseif strcmp(srcregion,'circle')
    fid = fopen([workpath,'/synthetics/synsrcloc.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(radi),'.diam',num2str(diam),'T',num2str(round(Twin)),'p',irunstr,'_grd'],'w');  
  end
  tmpgrid = xygrid;
  tmpgrid(:,1:2)=round(sps/40*tmpgrid(:,1:2)); % *4 to get to 160 sps from 40.
  fprintf(fid,'%8.3f %8.3f %7.2f %7.2f %8.3f %8.3f\n',tmpgrid');
  fclose(fid);
  
  %%%in particular, starting indices of added greens function (i.e., templates), I need them for
  %%%forward convolution
  if strcmp(srcregion,'ellipse')
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(xaxis),'-',num2str(yaxis),'.diam',num2str(diam),'.else',num2str(fracelsew,2),...
      'nsat',num2str(nsat(inwrites)),'tdura',num2str(tdura),...
      'T',num2str(round(Twin)),'p',irunstr,'_stind'],'w');
  elseif strcmp(srcregion,'rectangle')
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(ll(1)),num2str(ll,2),num2str(ur(1)),num2str(ur,2),'.diam',num2str(diam),...
      '.else',num2str(fracelsew,2),'nsat',num2str(nsat(inwrites)),'tdura',...
      num2str(tdura),'T',num2str(round(Twin)),'p',irunstr,'_stind'],'w');
  elseif strcmp(srcregion,'circle')
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(radi),'.diam',num2str(diam),'.else',num2str(fracelsew,2),...
      'nsat',num2str(nsat(inwrites)),'tdura',num2str(tdura),...
      'T',num2str(round(Twin)),'p',irunstr,'_stind'],'w');  
  end
  tmp = squeeze(greensts{inwrites}{1});
  fprintf(fid,'%d \n', tmp');
    
end

% keyboard

%% extract and validate the added impulses of template 
if testsrcflag
  insat = 2;  %which saturation to look at
  disp(nsat(insat));
  wlensec = 30; %how long the window to test
  bufsec = 1; %need some buffer window for later CC alignment
  buffer = bufsec*sps;
  wlensecb = wlensec+bufsec;  

  sigstab = [];  %target window of synthetic signals at all stations
  sigpnstab = [];  %target window of synthetic signals plus the noise
  sigconvstab = [];  %target window of reproduced synthetic signals using convolution
  noistab = [];  %target window of noise

  srcpairstab = cell(nsta,1);  %sources whose zero-crossing with time range, [indtarvl rnoff12 rnoff13 amp]

  for ista = 1: nsta
    % ista = 3;  
    tgreen = zcrosses(ista);  % anchor time of templates, here choose as zero-crossing time
    indstsig = 1;  % starting index of the simulated signal to test
    indedsig = wlensecb*sps+indstsig-1; % ending index of the simulated signal to test
  %   source = squeeze(sources(1:size(greensts{insat}{1},1), :, insat));  % [indtori indttrvl indtarvl rnoff12 rnoff13 amp];
    n=writes(insat);
    if strcmp(distrloc,'uniform')
  %     if forcesep
        a = squeeze(sources(1:n,:,insat));
        b = a(any(a,2),:);
  %     else
  %       a = sources(1:n,:);
  %       b = a(any(a,2),:);
  %     end
      source=b;
      off = tmpgrid(source(:,2),1:2); %note that 'tmpgrid' has the desired sps
    elseif strcmp(distrloc,'custompdf')
      a = sources(1:n,:);
      b = a(any(a,2),:);
      source=b;
      off = source(:,1:2);
    end
    if ista == 1
      greenst = greensts{insat}{1}; % the starting index of each added template, context of full length
    elseif ista <=3
      %ind - rnoff is the arrival time in index at sta 2 and 3
      %note the sign here, if off12 >0, move 2 to the right to align with 1, meaning that 2 is
      %earlier than 1, ie., tarvl2 < tarvl1. Be consistent all the time
      greenst = greensts{insat}{1}-off(:, ista-1); 
    else
      greenst = pred_tarvl_at4thsta(stas(ista,:),off(:,1),off(:,2),greensts{insat}{1});
    end
    
    %%%you don't need all impulses, only some of them contribute to the length of truncated record
    %%%cut out the green indice and sources that contribute
    %%%Note that some zero-crossings might appear later than 'indedsig', but the corresponding start
    %%%index of the green is actually smaller than 'indedsig', ie. their 'greenst' <= 'indedsig'
    induse = find(greenst>=skiplen-greenlen+indstsig & greenst<=indedsig+skiplen-1);
    greenst = greenst(induse);
    source = source(induse,:);
    off = off(induse,:);
    source = [source(:,1) off];
    
    greenzc = greenst+tgreen; % index of approximate zero-crossing
    source(:, 1) = source(:, 1)+zcrosses(1);  % now 'greenzc' should be the same as 1st col of 'source'
    impamp = zeros(max(greenzc)+20,1);
    for i = 1: length(greenzc)
      impamp(greenzc(i)) = impamp(greenzc(i))+1;  %assuming every source has an amp of 1
    end
    
    %note that some length of the full simulation 'skiplen' was skipped in 'synthgen.m'
    sigpntemp = squeeze(synths(:,ista,insat));  % simulated signal with white noise
    noitemp = synth(skiplen:winlen,ista);   % account for the skipped part
    sigtemp = sigpntemp-noitemp; % subtract the white noise
    sigpnstab(:,ista) = sigpntemp(indstsig:indedsig); % now focus on the part of interest
    sigstab(:,ista) = sigtemp(indstsig:indedsig);
    noistab(:,ista) = noitemp(indstsig:indedsig);
    
    sigconvtmp = conv(green(:,ista),impamp,'full');
    indtrc = tgreen+skiplen;  % starting index for truncation
    sigconvstab(:,ista) = sigconvtmp(indtrc+indstsig-1:indedsig+indtrc-1);  % cut accordingly
    
    % figure
    % plot(sig,'k'); hold on
    % plot(sigconv+3,'b');
    % plot(sigconv-sig+1,'r-')
    
    %%%Can we reproduce the synthetics with truncated convolution? YES
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
    plot([indstsig indstsig],ax.YLim,'--','color',[.5 .5 .5]);
    plot([indedsig indedsig],ax.YLim,'--','color',[.5 .5 .5]);
    legend(p1,'Synthetic random impulses');
    xran1 = ax.XLim; axpos1 = ax.Position(1);
    subplot(413);
    plot(sigstab(:,ista),'b'); hold on
    ax=gca; xran2 = ax.XLim;
    shrink(ax,xran1/xran2,1);
    ax.Position(1)=axpos1;  
    plot(sigconvstab(:,ista),'k');
    legend('Truncated synthetic signal','Truncated signal from convolution');
    subplot(414);
    plot(sigstab(:,ista)-sigconvstab(:,ista),'k'); hold on
    legend('Difference');
    ax=gca; xran3 = ax.XLim;
    shrink(ax,xran1/xran3,1);
    ax.Position(1)=axpos1;
  %   keyboard

  %   %%%For showing purpose
  %   widin = 6;  % maximum width allowed is 8.5 inches
  %   htin = 5;   % maximum height allowed is 11 inches
  %   nrow = 4;
  %   ncol = 1;
  %   f = initfig(widin,htin,nrow,ncol);  % create a figure with subplots
  %   
  %   pltxran = [0.1 0.9]; pltyran = [0.1 0.9];
  %   pltxsep = 0.02; pltysep = 0.04;
  %   %get the locations for each axis
  %   axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
  %   
  %   ax = f.ax(1);
  %   plot(ax,green(:,ista),'r'); hold on
  %   xlim(ax,[0 greenlen]);
  %   nolabels(ax,1);
  %   legend(ax,'Template (unfiltered)');
  %   ylabel(ax,'Amplitude');
  %   shrink(ax,wlensecb*sps/greenlen,1);
  %   ax.Position(1)=axpos(1,1);
  %   longticks(ax,1.5);
  %   axsym(ax,2);
  %   axranexp(ax,6,10);
  %   text(ax,0.05,0.85,stas(ista,:),'unit','normalized');
  %   
  %   ax = f.ax(2);
  %   imptemp = find(impamp>0);
  %   p1=stem(ax,imptemp-skiplen, impamp(imptemp),'b','MarkerSize',4); hold on;
  %   % p1=stem((1:size(impamp))-skiplen, impamp,'b','MarkerSize',4); hold on;
  %   % ax = gca;
  %   % plot([indstsig indstsig],ax.YLim,'--','color',[.5 .5 .5]);
  %   % plot([indedsig indedsig],ax.YLim,'--','color',[.5 .5 .5]);
  %   nolabels(ax,1);
  %   longticks(ax,3);
  %   legend(ax,p1,'Synthetic random impulses');
  %   %We only show the impulses whose zero-cross is inside 'wlensecb', even if later ones contribute
  %   %the trace as well, but their contribution is small. To avoid talking about them, we need to taper
  %   %both ends of the signal before deconvolution
  %   xlim(ax,[0 wlensecb*sps]);  
  %   ylabel(ax,'Amplitude');
  %   axranexp(ax,2,10);
  %   
  %   ax = f.ax(3);
  %   plot(ax,sigconvstab(:,ista),'k');
  %   xlim(ax,[0 wlensecb*sps]);
  %   nolabels(ax,1);
  %   longticks(ax,3);
  %   legend(ax,'Truncated signal from direct convolution');
  %   ylabel(ax,'Amplitude');
  %   axranexp(ax,6,10);
  %   
  %   % noistd = 5.e-2;
  %   % rng(seed,'twister');
  %   % noi = noistd*randn(length(sigconv),1);
  %   % sigplnoi = sigconv+noi;
  %   ax = f.ax(4);
  %   plot(ax,sigpnstab(:,ista),'k');
  %   xlim(ax,[0 wlensecb*sps]);
  %   longticks(ax,3);
  %   legend(ax,sprintf('Gaussian noise with std=%.2f added',noistd));
  % %   legend(ax,sprintf('Uniform noise with a max amp of %.1e added',mampnoi));
  %   xlabel(ax,sprintf('Samples at %d Hz',sps));
  %   ylabel(ax,'Amplitude');
  %   axranexp(ax,6,10);
    
    
    %now 'source' contain impulse triplets that contribute to the trace segment of interest, but we
    %want to get ones whose zero-crossing falls into the window.
    source(:,1) = source(:,1)-skiplen;  % index after skipping 
    srcpairb = source(source(:,1)>=indstsig & source(:,1)<=indedsig, :);
    srcpairb = sortrows(srcpairb,1,'ascend');
    srcpairstab{ista} = srcpairb; %ideally for each station this should be the same, but coincidence is possible  
    
  end

  %notify if sources are the same at diff stations
  if ~isequaln(srcpairstab{1,1},srcpairstab{2,1}) || ~isequaln(srcpairstab{1,1},srcpairstab{3,1})
    disp('Extracted sources at different stations are not the same, re-check!');
  end

end

