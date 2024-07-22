function synthshift_onespot_fn(perctrial,nsat,distr,Twin,tdura,seed,timetype,...
  normflag,tempflag,ftrans,pltsynflag,pltnewsynflag,testsrcflag,irun)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to 'synthshift_chao.m', this code is to make synthetic seismograms
% from LFE templates.
% --This code becomes the function to make synthetics which is called by
% But this code tries to put all srcs at one and only one spot, then add
% different level of noise, which is diff precent of noise that is made of
% the same amp spectrum of synthetics and randomized phase, *[times the
% envelope of data (so that the noise has the same fluctuation as data)]
% 1. in order to see what percentage of noise would make it look like data
% 2. now there can be 2 end-member synthetic tests: one is different saturation
%  level, diff source area, and no noise (or little); one is same spot source,
%  plus diff satur level and diff noise level
% 3. It should have nothing to do with data, the amp spectrum comes from
%   the synthetics with a certain saturation level, (and a certain size)
% 4. Since technically both synthetics have a uniform amp distribution
%   (we now know synthetics have a smaller amp fluctuation compared with
%   data), you don't have to multiply the envelope of synthetics with noise.
%   You just need to match their median envelope for 100% noise
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/05/05
% Last modified date:   2024/05/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% default settings
defval('perctrial',0.1*(0:2:16)');  %ratio of noise to add into synthetics
defval('nsat',[0.1 0.4 1 2 4 10 20 40 100]);  %saturation level
defval('distr','UN'); %specify the amplitude-frequency (counts) distribution as uniform
defval('Twin',0.5*3600+3+2*ceil(2048/160)); %length of each simulation
defval('tdura',0.25); %duration of templates
defval('seed',3); %for random number generator
defval('timetype','tori');  %specify which time is used for source
defval('normflag',0); %whether to normalize templates
defval('tempflag','chao'); %which templates to use
defval('ftrans','interpchao');  %specify regime for transformation from time offset to map location
defval('pltsynflag',0);  %whether to show a plot of raw synthetics
defval('pltnewsynflag',0);  %whether to show a plot of new synthetics with added noise
defval('testsrcflag',0);  %whether to test if synthetics can be reproduced by convolution
defval('irun',[]);  %which run out of the total # of runs of simulations

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

%% generate synthetic sources
%%%Specify the amplitude-frequency (counts) distribution
%   distr='UN'  % uniform distribution

%%%specify distribution for source location
% distrloc = 'custompdf'; %using a custom PDF function
%   distrloc = 'uniform'; %uniformly random in a specified region,

%   Twin=0.5*3600+3+2*ceil(greenlen/sps); %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
winlen=Twin*sps+1;
%   winlen = size(xnoi,1)+greenlen-1;
%   Twin = winlen/sps;
skiplen=greenlen;
%%%for the duration of templates, there are several options
%%%1. (ppeak-npeak)*2 of the bb template: 44;54;38
%%%2. direct eyeballing for between zerocrossings: ~65
%%%3. binned peak-to-peak separation for decently saturated unfiltered synthetics: ~37
%%%4. similar to 3, but synthetics are filtered first: ~37
%   tdura = 0.4;  % duration from Chao's broadband template, width is about 795-730=65 spls at 160 hz
%   tdura = 0.25; % start to use on 2023/12/07
satn=1/tdura*Twin   % if just saturated, how many templates can be fit in? a single peak is ~20 samples wide; maybe a little less (at 100 sps).
%Twin is window duration in seconds. Events can fall within Twin of
%the start, but the synthetics will go to Twin*(sample rate)+Greenlen to
%avoid checking for subscript overrrun.  When "synths" are written from "synth" in the
%subroutine, Greenlen from the start and end will not be written.
fracelsew=0; %0.25 %0.5; %0.6; %The ratio of elsewhere events to local events.  Zero for Discrete Ide.
% nsat=[0.4 1 2 4 10 20 40 100];  % times of saturation
% nsat=[0.1 0.4 1 2 4 10 20 40 100];  % times of saturation
%   nsat=[0.05];  % times of saturation
nnsat = length(nsat);
writes=round(nsat*satn) %how many templates to throw in, under different degrees of saturation

synth=zeros(winlen+greenlen,nsta);

nouts=length(writes);
%   seed=round(writes(4)/5e3); %for random number generator
%   seed=2;

%%%specify which time is uniform in time
% timetype = 'tarvl';
%   timetype = 'tori';

%%%specify if considering the physical size of each source
% physicalsize = 1;
%   physicalsize = 0;

%%%specify regime for transformation from time offset to map location
% ftrans = 'interpArmb';
% ftrans = 'interpArmbreloc';
%   ftrans = 'interpchao';

%%%whether to plot to check the synthetics
%   pltsynflag = 0;
%   pltsynflag = 1;

%%%whether to plot to check the new synthetics with added noise
%   pltnewsynflag = 0;
%   pltnewsynflag = 1;

%%%specify if forcing a min speration of arrival time for 2 events from the same spot
% forcesep = 1;
%   forcesep = 0;

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

%%%choose only one spot
%   xygrid = xygrid(xygrid(:,1)==0 & xygrid(:,2)==0, :); %use [0,0]
xygrid = xygrid(xygrid(:,1)==2/iup & xygrid(:,2)==2/iup, :);  %use [2,2]

diam=0;

%USE broadband templates to generate synthetics, then bandpass before deconvolution
[synths,mommax,sources,greensts]=csplaw3c(writes,winlen,skiplen,synth,green,b,...
  xygrid,sps,fracelsew,seed,timetype,ftrans,stas);

%%
%%% a simple plot the synthetics
if pltsynflag
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
  
  %only the lowest saturation to see details
  figure
  i=1;
  tmp = synths(:,:,i);
  plot(tmp(:,1),'r'); hold on
  plot(tmp(:,2),'b');
  plot(tmp(:,3),'k');
  text(0.95,0.9,sprintf('%.1f x saturation',nsat(i)),'Units','normalized','HorizontalAlignment',...
    'right');
  xlim([0 40*sps]);
  axranexp(gca,6,20);
  xlabel(sprintf('Samples at %d Hz',sps),'FontSize',12);
  
end

%   insat = 2;  %which saturation to look at
%   n=writes(insat);
%   a = squeeze(sources(1:n,:,insat));
%   b = a(any(a,2),:);
%   source=b;
%   off = tmpgrid(source(:,2),1:2); %note that 'tmpgrid' has the desired sps
%   [~,off(:,3)] = pred_tarvl_at4thsta(stas(4,:),off(:,1),off(:,2));

%% compose new synthetic waveform with diff noise level
%different percent of noise
%   perctrial = 0.1*(0:2:16)';
ntrial = length(perctrial);

%%%loop for noise level
for iperc = 1: ntrial
  
  perc = perctrial(iperc);
  disp(perc);
  
  synnew = synths;  %time, station, sat
  
  %%%loop for saturation level
  for insat = 1: nnsat
    %   insat = 1;
    disp(nsat(insat));
    
    %%%load synthetics of certain saturation level
    %       optseg = load(strcat(workpath,fname,num2str(nsat(insat))));
    optseg = synths(:,:,insat);
    
    %       %%%load sources
    %       synsrc = load(strcat(workpath,fname,num2str(nsat(insat)),'_sources'));
    %       if strcmp(distrloc, 'uniform')
    %         xygrid = load([workpath,'/synthetics/synsrcloc.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
    %           '_',num2str(xaxis),'-',num2str(yaxis),'.diam',num2str(diam),'_grd']);
    %         tmp = xygrid(synsrc(:,2),:);
    %         synsrc = [synsrc(:,1) tmp(:,1:4) ones(length(tmp),1)];  %[indtarvl, off12, off13, loce, locn, amp]
    %       elseif strcmp(distrloc, 'custompdf')
    %         [loc, indinput] = off2space002(synsrc(:,2:3),sps,ftrans,0);
    %         synsrc = [synsrc(:,1) loc(:,1:4) ones(length(loc),1)];  %[indtarvl, off12, off13, loce, locn, amp]
    %       end
    
    %%% make noise by amp spectrum of synthetics and random phase
    %obtain the amp and phase spectra of records via fft
    nfft = size(optseg,1); % number of points in fft
    [xf,ft,amp,pha] = fftspectrum(optseg(:,1:end), nfft, sps,'twosided');
    
    %uniform, random phase with the same span [-pi,pi];
    mpharan = minmax(pha');
    
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
    env = envelope(detrend(optseg(:,1:end)));
    % sclfact = env./range(env);  %normalize
    sclfact = env;
    %       %3. To make noise have the median amp of decon sources from data
    %       sclfact = 6.4461e-01*mean(spread(1:3))/2;
    
    %scale the noise so that it has amp fluctuation as data
    %       xnoi = xnoi .* sclfact;
    
    %normalize so that noi has the same median env as data
    envn = envelope(detrend(xnoi));
    xnoi = xnoi .* median(env) ./ median(envn);
    median(env) - median(envelope(detrend(xnoi)))
    
    %%
    if pltnewsynflag
      [f] = initfig(8,4.5,3,1);
      optaxpos(f,3,1,[],[],[],0.06);
      ax = f.ax(1); hold(ax,'on'); ax.Box = 'on';
      plot(ax,(1:size(optseg,1))/sps,optseg(:,1),'r','linew',1); hold on
      plot(ax,(1:size(optseg,1))/sps,optseg(:,2),'b','linew',1);
      plot(ax,(1:size(optseg,1))/sps,optseg(:,3),'k','linew',1);
      text(ax,0.98,0.9,sprintf('syn from satur=%.1f',nsat(insat)),'Units','normalized',...
        'HorizontalAlignment','right');
      %       xlabel('Samples at 160 sps');
      %         xlabel('Time (s)');
      %         ylabel('Amplitude');
      %         xlim([0 size(optseg,1)]/sps);
      xran = [2 27];
      xlim(ax,xran);
      axranexp(ax,6,10);
      yran = ax.YLim;
      longticks(ax,4);
      %       ylim([-1.2 1.2]);
      
      ax = f.ax(2); hold(ax,'on'); ax.Box = 'on';
      plot(ax,(1:size(optseg,1))/sps,xnoi(:,1),'r','linew',1); hold on
      plot(ax,(1:size(optseg,1))/sps,xnoi(:,2),'b','linew',1);
      plot(ax,(1:size(optseg,1))/sps,xnoi(:,3),'k','linew',1);
      text(ax,0.98,0.9,'100% assembled noise','Units','normalized','HorizontalAlignment','right');
      %       xlabel('Samples at 160 sps');
      %         xlabel('Time (s)');
      %         ylabel('Amplitude');
      %         xlim([0 size(optseg,1)]/sps);
      xlim(ax,xran);
      ylim(ax,yran);
      longticks(ax,4);
      
      ax = f.ax(3); hold(ax,'on'); ax.Box = 'on';
      plot(ax,(1:size(optseg,1))/sps,optseg(:,1)+xnoi(:,1),'r','linew',1); hold on
      plot(ax,(1:size(optseg,1))/sps,optseg(:,2)+xnoi(:,2),'b','linew',1);
      plot(ax,(1:size(optseg,1))/sps,optseg(:,3)+xnoi(:,3),'k','linew',1);
      text(ax,0.98,0.9,'syn + 100% noise','Units','normalized','HorizontalAlignment','right');
      %       xlabel('Samples at 160 sps');
      xlabel(ax,'Time (s)','FontSize',10);
      ylabel(ax,'Amplitude','FontSize',10);
      %         xlim([0 30]);
      xlim(ax,xran);
      ylim(ax,yran);
      longticks(ax,4);
    end
    %       keyboard
    %%
    %some percent of assembled noise
    noiseg = xnoi .*perc;
    
    %add simulated noise to the current burst window
    tmp = synths(:,:,insat);
    synnew(:,:,insat) = tmp + noiseg;
    
    %       %%%filter data
    %       for ista = 1:nsta
    %         synnew(:,ista,insat) = Bandpass(synnew(:,ista,insat), sps, losig, hisig, 2, 2, 'butter');
    %       end
    
    %%
    %%%a glance of signal
    if pltnewsynflag
      figure
      subplot(311)
      plot((1:size(synths,1))/sps,synths(:,1,insat),'r','linew',1); hold on
      plot((1:size(synths,1))/sps,synths(:,2,insat),'b','linew',1);
      plot((1:size(synths,1))/sps,synths(:,3,insat),'k','linew',1);
      text(0.95,0.9,sprintf('Single-spot synthetics with sat=%.1f',nsat(insat)),'Units','normalized',...
        'HorizontalAlignment','right');
      %       xlabel('Samples at 160 sps');
      xlabel('Time (s)');
      ylabel('Amplitude');
      % xlim([0 size(xnoi,1)]/sps);
      xlim([0 30]);
      ax = gca; yran = ax.YLim;
      
      subplot(312)
      plot((1:size(synths,1))/sps,noiseg(:,1),'r','linew',1); hold on
      plot((1:size(synths,1))/sps,noiseg(:,2),'b','linew',1);
      plot((1:size(synths,1))/sps,noiseg(:,3),'k','linew',1);
      text(0.95,0.9,sprintf('%d%% noise',round(perc*100)),'Units','normalized','HorizontalAlignment','right');
      %       xlabel('Samples at 160 sps');
      xlabel('Time (s)');
      ylabel('Amplitude');
      % xlim([0 size(xnoi,1)]/sps);
      xlim([0 30]);
      ylim(yran);
      
      subplot(313)
      plot((1:size(synths,1))/sps,synnew(:,1,insat),'r','linew',1); hold on
      plot((1:size(synths,1))/sps,synnew(:,2,insat),'b','linew',1);
      plot((1:size(synths,1))/sps,synnew(:,3,insat),'k','linew',1);
      text(0.95,0.9,'Synthetics+noise','Units','normalized','HorizontalAlignment','right');
      %       xlabel('Samples at 160 sps');
      xlabel('Time (s)');
      ylabel('Amplitude');
      % xlim([0 size(xnoi,1)]/sps);
      xlim([0 30]);
      ylim(yran);
    end
    
    %       close all
  end
  
  %% save files
  for inwrites=1:length(writes)
    n=writes(inwrites);
    %%%seismograms
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',...
      'onespot','.noi',num2str(perc),'.else',num2str(fracelsew,2),...
      'nsat',num2str(nsat(inwrites)),'tdura',num2str(tdura),...
      'T',num2str(round(Twin)),'p',irunstr],'w');
    towrite=squeeze(synnew(:,:,inwrites));
    if nsta == 3
      fprintf(fid,'%13.5e %13.5e %13.5e \n', towrite');
    elseif nsta == 4
      fprintf(fid,'%13.5e %13.5e %13.5e %13.5e \n', towrite');
    end
    fclose(fid);
    
    %%%source info
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',...
      'onespot','.noi',num2str(perc),'.else',num2str(fracelsew,2),...
      'nsat',num2str(nsat(inwrites)),'tdura',num2str(tdura),...
      'T',num2str(round(Twin)),'p',irunstr,'_sources'],'w');
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
    
    %%%save the source grid location
    %   fid = fopen([workpath,'/synthetics/synsrcloc.',srcregion(1:3),'.grd'],'w');
    fid = fopen([workpath,'/synthetics/synsrcloc.',distr,'.',int2str(sps),'sps.',...
      'onespot','.noi',num2str(perc),...
      'T',num2str(round(Twin)),'p',irunstr,'_grd'],'w');
    tmpgrid = xygrid;
    tmpgrid(:,1:2)=round(sps/40*tmpgrid(:,1:2)); % *4 to get to 160 sps from 40.
    fprintf(fid,'%8.3f %8.3f %7.2f %7.2f %8.3f %8.3f\n',tmpgrid');
    fclose(fid);
    
    %%%in particular, starting indices of added greens function (i.e., templates), I need them for
    %%%forward convolution
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',...
      'onespot','.noi',num2str(perc),'.else',num2str(fracelsew,2),...
      'nsat',num2str(nsat(inwrites)),'tdura',num2str(tdura),...
      'T',num2str(round(Twin)),'p',irunstr,'_stind'],'w');
    tmp = squeeze(greensts{inwrites}{1});
    fprintf(fid,'%d \n', tmp');
    
  end
  
  % keyboard
  
  %% testing, extract and validate the added impulses of template
  if testsrcflag
    %%%load sources
    insat = 1;
    n=writes(insat);
    a = squeeze(sources(1:n,:,insat));
    b = a(any(a,2),:);
    synsrc=b;
    tmp = tmpgrid(synsrc(:,2),:);
    synsrc = [synsrc(:,1) tmp(:,1:4) ones(length(tmp),1)];  %[indtarvl, off12, off13, loce, locn, amp]
    
    % keyboard
    %%%load starting indices of added sources at sta 1
    synsrcstind = squeeze(greensts{insat}{1});
    
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
      %           aaa=synths(:,ista,insat);
      sigpntemp = synths(:,ista,insat);  % simulated signal with white noise
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
  
end














