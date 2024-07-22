function synthshift_regandnoi_fn(srcregion,xaxis,yaxis,perctrial,nsat,distr,...
  Twin,tdura,seed,distrloc,timetype,physicalsize,forcesep,testsrcflag,normflag,...
  tempflag,ftrans,pltsynflag,pltnewsynflag,pltflag,irun)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function called by 'synthshift_regandnoi_multibsts.m' to generate synthetics
% with diff. saturation levels made by sources from the combination of a small 
% set of region sizes and with a small set of noise levels. 
% For generating synthetics EITHER from various region sizes with no added 
% noise, refer to 'synthshift_chao_fn' and 'synthshift_chao_multibsts',
% OR from a single spot with various noise levels, refer to
% 'synthshift_onespot_fn' and 'synthshift_onespot_multibsts'.
% 
% The first part is similar to 'synthshift_chao_fn', ie., generating synthetics
% from a specific region size and a specific saturation level. For each reg size
% and sat level's synthetics, obtain the amp spectrum and use a random phase 
% spectrum, then you can assemble the 100% noise. Scale the noise, and add it to
% the original synthetics, then you get the synthetics with a certain sat level
% from a certain reg size that is contaminated by a certain noise level.
% 
%
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/07/20
% Last modified date:   2024/07/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% default settings
defval('srcregion','ellipse');  %specify shape of the source region
defval('xaxis',1.75); %axis length of the same ellipse of my 4-s catalog
defval('yaxis',1.25);
defval('perctrial',0.1*(4:2:8)');  %ratio of noise to add into synthetics
defval('nsat',[0.4 1 2 4 10 20]);  %saturation level
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
defval('pltsynflag',0);  %whether to show a plot of raw synthetics
defval('pltnewsynflag',0);  %whether to show a plot of new synthetics with added noise
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

% keyboard

%% synthetics generation
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

if strcmp(distrloc,'uniform')
  
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
  
end

%% a simple plot the synthetics
if pltsynflag
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

%% extract and validate the added impulses of template 
if testsrcflag
  insat = 2;  %which saturation to look at  
  disp(nsat(insat));
  
  srcvalidate(insat,sps,zcrosses,writes,distrloc,sources,tmpgrid,...
    greensts,stas,skiplen,greenlen,synths,synth,winlen,green);

end

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
    
    %%%load synthetics of certain saturation level for all stations
    optseg = synths(:,:,insat);
        
    %%% make noise by amp spectrum of synthetics and random phase
    %obtain the amp and phase spectra of records via fft
    nfft = size(optseg,1); % number of points in fft
    [xf,ft,amp,pha] = fftspectrum(optseg(:,1:end), nfft, sps,'twosided');
    
    %phase range of real data
    mpharan = minmax(pha');
    
    %vary the seed number, note that diff 'iperc' need to scale with the SAME
    %100% noise, so this new seed has to be irrelavant to 'iperc'
    seed2=seed*10+insat;
    rng(seed2);
    
    %uniform, random phase with the same span [-pi,pi];
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
    
    envdiff=median(env) - median(envelope(detrend(xnoi)));
    if envdiff > 1e-6
      disp('100% noise and synthetics do not have the same median env.')
    end
    
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
    
    %% add the scaled noise with assgined 'perc' to original to build new synthetics
    %some percent of assembled noise
    noiseg = xnoi .*perc;
    
    %add simulated noise to the current burst window
    tmp = synths(:,:,insat);  %original synthetics
    synnew(:,:,insat) = tmp + noiseg; %new synthetics with 'perc' noise added 
    
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
    if strcmp(srcregion,'ellipse')
      fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',...
        srcregion(1:3),'_',num2str(xaxis),'-',num2str(yaxis),'.diam',...
        num2str(diam),'.noi',num2str(perc),'.else',num2str(fracelsew,2),'nsat',...
        num2str(nsat(inwrites)),'tdura',num2str(tdura),...
        'T',num2str(round(Twin)),'p',irunstr],'w');
    elseif strcmp(srcregion,'rectangle')
      fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',...
        srcregion(1:3),'_',num2str(ll(1)),num2str(ll(2)),num2str(ur(1)),...
        num2str(ur(2)),'.diam',num2str(diam),...
        '.noi',num2str(perc),'.else',num2str(fracelsew,2),'nsat',...
        num2str(nsat(inwrites)),'tdura',num2str(tdura),...
        'T',num2str(round(Twin)),'p',irunstr],'w');
    elseif strcmp(srcregion,'circle')
      fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',...
        srcregion(1:3),'_',num2str(radi),'.diam',num2str(diam),...
        '.noi',num2str(perc),'.else',num2str(fracelsew,2),'nsat',...
        num2str(nsat(inwrites)),'tdura',num2str(tdura),...
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
      fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',...
        srcregion(1:3),'_',num2str(xaxis),'-',num2str(yaxis),'.diam',...
        num2str(diam),'.noi',num2str(perc),'.else',num2str(fracelsew,2),'nsat',...
        num2str(nsat(inwrites)),'tdura',num2str(tdura),...
        'T',num2str(round(Twin)),'p',irunstr,'_sources'],'w');
    elseif strcmp(srcregion,'rectangle')
      fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',...
        srcregion(1:3),'_',num2str(ll(1)),num2str(ll(2)),num2str(ur(1)),...
        num2str(ur(2)),'.diam',num2str(diam),...
        '.noi',num2str(perc),'.else',num2str(fracelsew,2),'nsat',...
        num2str(nsat(inwrites)),'tdura',num2str(tdura),...
        'T',num2str(round(Twin)),'p',irunstr,'_sources'],'w');
    elseif strcmp(srcregion,'circle')
      fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',...
        srcregion(1:3),'_',num2str(radi),'.diam',num2str(diam),...
        '.noi',num2str(perc),'.else',num2str(fracelsew,2),'nsat',...
        num2str(nsat(inwrites)),'tdura',num2str(tdura),...
        'T',num2str(round(Twin)),'p',irunstr,'_sources'],'w');
    end
    if strcmp(distrloc,'uniform')
      %     if forcesep
      a = squeeze(sources(1:n,:,inwrites));
      b = a(any(a,2),:);
      size(b)
      truncsources=b;
      if strcmp(timetype,'tarvl')
        fprintf(fid,'%13i %5i \n', truncsources'); %arrival; loc ind of grid
      elseif strcmp(timetype,'tori')
        fprintf(fid,'%13i %5i %13i %13i \n', truncsources'); %arrival; loc ind of grid; origin; travel time
      end
      fclose(fid);
    end
        
    %save the source grid location
    if strcmp(srcregion,'ellipse')
      fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',...
        srcregion(1:3),'_',num2str(xaxis),'-',num2str(yaxis),'.diam',...
        num2str(diam),'.noi',num2str(perc),'.else',num2str(fracelsew,2),'nsat',...
        num2str(nsat(inwrites)),'tdura',num2str(tdura),...
        'T',num2str(round(Twin)),'p',irunstr,'_grd'],'w');
    elseif strcmp(srcregion,'rectangle')
      fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',...
        srcregion(1:3),'_',num2str(ll(1)),num2str(ll(2)),num2str(ur(1)),...
        num2str(ur(2)),'.diam',num2str(diam),...
        '.noi',num2str(perc),'.else',num2str(fracelsew,2),'nsat',...
        num2str(nsat(inwrites)),'tdura',num2str(tdura),...
        'T',num2str(round(Twin)),'p',irunstr,'_grd'],'w');
    elseif strcmp(srcregion,'circle')
      fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',...
        srcregion(1:3),'_',num2str(radi),'.diam',num2str(diam),...
        '.noi',num2str(perc),'.else',num2str(fracelsew,2),'nsat',...
        num2str(nsat(inwrites)),'tdura',num2str(tdura),...
        'T',num2str(round(Twin)),'p',irunstr,'_grd'],'w');
    end
    tmpgrid = xygrid;
    tmpgrid(:,1:2)=round(sps/40*tmpgrid(:,1:2)); % *4 to get to 160 sps from 40.
    fprintf(fid,'%8.3f %8.3f %7.2f %7.2f %8.3f %8.3f\n',tmpgrid');
    fclose(fid);
    
    %%%in particular, starting indices of added greens function (i.e., templates), I need them for
    %%%forward convolution
    if strcmp(srcregion,'ellipse')
      fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',...
        srcregion(1:3),'_',num2str(xaxis),'-',num2str(yaxis),'.diam',...
        num2str(diam),'.noi',num2str(perc),'.else',num2str(fracelsew,2),'nsat',...
        num2str(nsat(inwrites)),'tdura',num2str(tdura),...
        'T',num2str(round(Twin)),'p',irunstr,'_stind'],'w');
    elseif strcmp(srcregion,'rectangle')
      fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',...
        srcregion(1:3),'_',num2str(ll(1)),num2str(ll(2)),num2str(ur(1)),...
        num2str(ur(2)),'.diam',num2str(diam),...
        '.noi',num2str(perc),'.else',num2str(fracelsew,2),'nsat',...
        num2str(nsat(inwrites)),'tdura',num2str(tdura),...
        'T',num2str(round(Twin)),'p',irunstr,'_stind'],'w');
    elseif strcmp(srcregion,'circle')
      n([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',...
        srcregion(1:3),'_',num2str(radi),'.diam',num2str(diam),...
        '.noi',num2str(perc),'.else',num2str(fracelsew,2),'nsat',...
        num2str(nsat(inwrites)),'tdura',num2str(tdura),...
        'T',num2str(round(Twin)),'p',irunstr,'_stind'],'w');
    end
    tmp = squeeze(greensts{inwrites}{1});
    fprintf(fid,'%d \n', tmp');
    
  end
  
  % keyboard
    
end

















