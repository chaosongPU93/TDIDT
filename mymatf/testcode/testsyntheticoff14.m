% testsyntheticoff14.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS code is supposed to test when the saturation level is so low that no
% arrivals at any station are too close in time. Then the deconvolution at
% each station should be able to resolve every source. The question is, in
% that case, should the off14 between the actually-matched arrivals at stas
% 1 and 4 exactly the same as synthetic and prediction from plane fit model,
% or the difference in arrivals at the 4th sta between the prediction and 
% actually-matched one. Ideally, the above statement is true. Otherwise,
% it is hard to tell if the distribution of off14 for a higher saturation
% level is real.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/07/10
% Last modified date:   2023/07/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
clc
close all

set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');

[scrsz, resol] = pixelperinch(1);

%% prepare templates (Green's functions), from 'lfetemp002_160sps.m' or Allan's templates
defval('normflag',0); %whether to normalize templates

% tempflag = 'allan';
tempflag = 'chao';

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
[off12con,off13con,cc] = constrained_cc_interp(tmpwletf(:,1:3)',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
offwlet1i(1) = 0;
offwlet1i(2) = round(off12con);
offwlet1i(3) = round(off13con);
if nsta>3 
  for ista = 4: nsta
    [mcoef,offwlet1i(ista)] = xcorrmax(tmpwletf(:,1), tmpwletf(:,ista), mshiftadd, 'coeff');
  end
end

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
  [~,imin] = min(green(:,ista));
  [~,imax] = max(green(:,ista));
  [~,zcrosses(ista)] = min(abs(green(imin:imax,ista)));
  zcrosses(ista) = zcrosses(ista)+imin-1;
  ppeaks(ista) = imax;
  npeaks(ista) = imin;
  
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

amprat(1,:) = minmax(greenf(:,1)')./minmax(greenf(:,2)');	% amp ratio between max at sta 3 and 2 or min
amprat(2,:) = minmax(greenf(:,1)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min  
amprat(3,:) = minmax(greenf(:,2)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min  
spread = range(greenf);   % range of the amp of template

%%%plot the unfiltered and filtered templates
% plt_templates(green,greenf,stas,[],[],lowlet,hiwlet,sps);

%just the filtered templates
% plt_templates_bp(greenf,stas,lowlet,hiwlet,sps);
% zcrosses
% ppeaks
% npeaks
% zcrossesf
% ppeaksf
% npeaksf
% keyboard

%% synthetics generation
%%%Specify the amplitude-frequency (counts) distribution 
distr='UN'  % uniform distribution

Twin=0.5*50+3+2*ceil(greenlen/sps); %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
winlen=Twin*sps+1;
skiplen=greenlen;
%%%for the duration of templates, there are several options
%%%1. (ppeak-npeak)*2 of the bb template: 44;54;38
%%%2. direct eyeballing for between zerocrossings: ~65
%%%3. binned peak-to-peak separation for decently saturated unfiltered synthetics: ~37
%%%4. similar to 3, but synthetics are filtered first: ~37
tdura = 0.4;  % duration from Chao's broadband template, width is about 795-730=65 spls at 160 hz
satn=1/tdura*Twin   % if just saturated, how many templates can be fit in? a single peak is ~20 samples wide; maybe a little less (at 100 sps).
    %Twin is window duration in seconds. Events can fall within Twin of
    %the start, but the synthetics will go to Twin*(sample rate)+Greenlen to
    %avoid checking for subscript overrrun.  When "synths" are written from "synth" in the
    %subroutine, Greenlen from the start and end will not be written.
fracelsew=0; %0.25 %0.5; %0.6; %The ratio of elsewhere events to local events.  Zero for Discrete Ide.
nsat=[0.1];  % times of saturation 
% nsat=[95 100 200];  % times of saturation 
writes=round(nsat*satn) %how many templates to throw in, under different degrees of saturation
rng('default');
%%%Here the noise is 'uniform' in time!
% noistd = 5e-2;
noistd = 2.e-7;
% noistd = 2.0e-4;
synth=noistd*(randn(winlen+greenlen+2*10,nsta)-0.5); %+2*10 a little extra, for jiggering 2nd & 3rd stations.% for ista=1:nsta

nouts=length(writes);
seed=round(writes(1)/5e3); %for random number generator

%%%specify which time is uniform in time
% timetype = 'tarvl';
timetype = 'tori';

%%%specify if considering the physical size of each source
physicalsize = 1;
% physicalsize = 0;

%%%specify shape of the source region
srcregion='ellipse';
% srcregion='rectangle';
% srcregion='circle';

%%%specify regime for transformation from time offset to map location 
% ftrans = 'interpArmb';
% ftrans = 'interpArmbreloc';
ftrans = 'interpchao';

%%%specify if forcing a min speration of arrival time for 2 events from the same spot 
forcesep = 1;
% forcesep = 0;

  b=999. %>150 for uniform size distribution
  xygrid=load('xygridArmb'); %Made from /ARMMAP/MAPS/2021/interpgrid.m; format PGSS, PGSI, dx, dy
  size(xygrid) %this grid is 40 sps!
  
  %which location transformation to use
  if strcmp(ftrans,'interpArmb')
    xygrid = xygrid;
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
    diam=0.15;  %note that if diam is <150m, in sample space you'll have too many nonunqiue sources 
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
    semia = 1.75*(0.6:0.2:2.0);
    semib = 1.25*(0.6:0.2:2.0);
    nreg = length(semia);
    ireg = 3;
    xaxis = semia(ireg); %axis length of the same ellipse of my 4-s catalog
    yaxis = semib(ireg);
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

%   %whether to consider the physical size of each source
%   tmp = ceil(max(max(abs([xcut,ycut]))));
%   figure
%   axis equal
%   hold on; grid on; box on
%   scatter(shiftor(1),shiftor(2),15,'k','filled');
%   plot(xcut,ycut,'k-','linew',1);
%   plot(xygrid(:,3),xygrid(:,4),'.');
%   text(0.95,0.05,sprintf('%d unique sources',size(xygrid,1)),'Units','normalized',...
%     'HorizontalAlignment','right');
%   xran = [-tmp tmp];
%   yran = [-tmp tmp];
%   xlim(xran);
%   ylim(yran);
%   xlabel('E (km)');
%   ylabel('N (km)');
%   if physicalsize
%     %SEE Chao's NOTES: research/geometry of synthetics
%     rad=0.5*diam; %radius of blue circles, max. with no overlapping
%     areafrac=0.5; %fraction of area of single asperity
%     lod2=areafrac*2*sqrt(3)/pi; %lod2 is ("asperity diameter"/diam)^2, where diam is hex grid spacing
%     asprad=0.5*sqrt(lod2)*diam; %asperity radius
%     ang=0:0.1*pi:2*pi;
%     xdiam=rad*cos(ang);
%     ydiam=rad*sin(ang);
%     xasp=asprad*cos(ang);
%     yasp=asprad*sin(ang);
%     for ispot=1:size(xygrid,1)
%       plot(xygrid(ispot,3)+xdiam,xygrid(ispot,4)+ydiam,'b--');
%       plot(xygrid(ispot,3)+xasp,xygrid(ispot,4)+yasp,'r');
%     end
%   end
  size(xygrid)

  %whether to force a min speration if 2 events from the same spot separated by less than the template duration
  %USE unfiltered templates to generate synthetics, then bandpass before deconvolution
  if forcesep
    [synths,mommax,sources,greensts]=csplaw3d(writes,winlen,skiplen,synth,green,b,...
      xygrid,sps,fracelsew,seed,tdura,timetype,ftrans,stas); 
  else
    [synths,mommax,sources,greensts]=csplaw3c(writes,winlen,skiplen,synth,green,b,...
      xygrid,sps,fracelsew,seed,timetype,ftrans,stas);
  end
  
%% a simple plot the synthetics
% figure
% tmp = synths(:,:,2)-synth(skiplen:winlen,:);
% plot(tmp(:,1),'b');
% xlim([0 31*sps]);

% figure
% nrow = length(writes)+1;
% ncol = 1;
% subplot(nrow,ncol,1)
% hold on
% tmp = synth;
% if nsta == 3
%   plot(tmp(:,1),'r');
%   plot(tmp(:,2),'b');
%   plot(tmp(:,3),'k');
% else
%   color = jet(nsta);
%   for ista = 1: nsta
%     plot(tmp(:,ista),'Color',color(ista,:));
%   end
% end
% xlim([0 40*sps]);
% axranexp(gca,6,20);
% 
% for i = 1: length(writes)
% subplot(nrow,ncol,i+1)
% hold on
% tmp = synths(:,:,i);
% if nsta == 3
%   plot(tmp(:,1),'r');
%   plot(tmp(:,2),'b');
%   plot(tmp(:,3),'k');
% else
%   color = jet(nsta);
%   for ista = 1: nsta
%     plot(tmp(:,ista),'Color',color(ista,:));
%   end
% end
% text(0.95,0.9,sprintf('%.1f x saturation',nsat(i)),'Units','normalized','HorizontalAlignment',...
%   'right');
% xlim([0 40*sps]);
% axranexp(gca,6,20);
% end
% xlabel(sprintf('Samples at %d Hz',sps),'FontSize',12);

% keyboard

%%
%%%load synthetics of certain saturation level
insat = 1;
disp(nsat(insat));

STAopt = squeeze(synths(:,:,insat));

%%%load sources
n=writes(insat);
if forcesep
  a = squeeze(sources(1:n,:,insat));
  b = a(any(a,2),:);
else
  a = sources(1:n,:);
  b = a(any(a,2),:);
end
size(b)
synsrc=b;

tmpgrid = xygrid;
tmpgrid(:,1:2)=round(sps/40*tmpgrid(:,1:2)); % *4 to get to 160 sps from 40.
tmp = tmpgrid(synsrc(:,2),:);
synsrc = [synsrc(:,1) tmp(:,1:4) ones(length(tmp),1)];  %[indtarvl, off12, off13, loce, locn, amp]

%%%load starting indices of added sources at sta 1 
synsrcstind = squeeze(greensts{insat}{1});

% keyboard

%% testing, extract and validate the added impulses of template 
testsrcflag = 0;
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
    text(0.9,0.9,stas(ista,:),'HorizontalAlignment','right','Units','normalized');
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

%% filter data (and templates)
%%filter data
hisig=6.3; % this will give a similar spectral shape between template and signal
losig=1.8;
optseg = [];
for ista = 1:nsta
  optseg(:,ista) = Bandpass(STAopt(:,ista), sps, losig, hisig, 2, 2, 'butter');
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
[off12con,off13con,ccali,iloopoff,loopoff] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
  ccwlen,msftadd,loffmax,ccmin,iup);
% if a better alignment cannot be achieved, use 0,0
if off12con == msftadd+1 && off13con == msftadd+1
  off12con = 0;
  off13con = 0;
  cc12 = xcorr(optcc(:,1), optcc(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
  cc13 = xcorr(optcc(:,1), optcc(:,3),0,'normalized');
  cc23 = xcorr(optcc(:,2), optcc(:,3),0,'normalized');
  ccali = (cc12+cc13+cc23)/3;
  fprintf('Current window cannot be properly aligned, double-check needed \n');
end
off1i = zeros(nsta,1);
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

%%%Align and compute the RCC based on the entire win, and take that as the input signal!
optdat = [];  % win segment of interest
for ista = 1: nsta
  optdat(:, ista) = optseg(1+msftaddm-off1i(ista): end-msftaddm-off1i(ista), ista);
end

%location of the whole-win best alignment
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
for ista = 4:nsta
  [~,rcc1i(:,ista-3)] = RunningCC(sigsta(:,1), sigsta(:,ista), rccmwlen);
end
sigsta = detrend(sigsta(overshoot+1:end-overshoot, :));  %excluding the overshoot

cc12 = xcorr(sigsta(:,1), sigsta(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
cc13 = xcorr(sigsta(:,1), sigsta(:,3),0,'normalized');
cc23 = xcorr(sigsta(:,2), sigsta(:,3),0,'normalized');
cc1i = zeros(nsta-3,1);
for ista = 4:nsta
  cc1i(ista-3) = xcorr(sigsta(:,1), sigsta(:,ista),0,'normalized');
end
ccpair = [cc12 cc13 cc23];
mrcc = median(rcc);
mcc = (cc12+cc13+cc23)/3;

%if only use the mean RCC from pair 12 and 13
rcc = mean(rccpair(:,[1 2]), 2);
rccsat(:,insat) = rcc;

%%
%if you want a qucik look of the data that feeds to the deconvolution
pltdataflag = 1;
if pltdataflag    
  figure; hold on
  lsig = size(sigsta,1);
  if nsta == 3
    plot((1:lsig)/sps,sigsta(:,1),'r'); 
    plot((1:lsig)/sps,sigsta(:,2),'b');
    plot((1:lsig)/sps,sigsta(:,3),'k');
  else
    color = jet(nsta);
    for ista = 1: nsta
      p(ista)=plot((1:lsig)/sps,sigsta(:,ista),'Color',color(ista,:));
%       label{i}=stas(i,:);
    end
  end
  ax = gca;
  axsym(ax);
  plot((1:lsig)/sps,rcc*ax.YLim(2),'o','color',[.6 .6 .6],'markersize',2);
  text(0.95,0.1,sprintf('Saturation: %.1f',nsat(insat)),'Units','normalized','HorizontalAlignment','right');
%   xlim([0 20]);
  legend(p,stas);
end

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

%amp and density at each pixel for ground truth sources
ampgt = mean(impgt(:,[2 4 6]),2); %amp for all LFE catalog
density1d = density_pixel(impgt(:,7),impgt(:,8));
ampgtsum = sum_pixel(impgt(:,7),impgt(:,8),ampgt);
[impgtloc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
density1d = [impgtloc(:,1:2) density1d(:,3)];
ampgtsum1d = sortrows([impgtloc(:,1:2) ampgtsum(:,3)], 3);

%%%if you want to plot the ground truth
pltgtflag= 1;
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
  scatter(ax,off1i(2),off1i(3),20,'ks','filled','MarkerEdgeColor','k');
  title(ax,'Ground truth');

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

%%%if you want to plot the deconvolved sources 
pltsrcflag1 = 1;
if pltsrcflag1
  %plot the scatter of offsets, accounting for prealignment offset, == true offset
  xran = [-loff_max+off1i(2)-1 loff_max+off1i(2)+1];
  yran = [-loff_max+off1i(3)-1 loff_max+off1i(3)+1];
  offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
  offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
  cran = [0 lsig];
  f.fig = figure;
  f.fig.Renderer = 'painters';
  ax=gca;
  [ax,torispl,mamp] = plt_decon_imp_scatter(ax,impindepst,xran,yran,cran,offxran,offyran,...
    sps,50,'mean','tarvl');
  scatter(ax,off1i(2),off1i(3),20,'ks','filled','MarkerEdgeColor','k');
  title(ax,'Grouped Total');

end
% keyboard

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

%%%if you want to plot the deconvolved sources
pltsrcflag2 = 1;
if pltsrcflag2
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

end
% keyboard

%% 2ndary src removed, prediction of impulse tarvl at 4th sta given sources and empirical off14-src relation
%%%carry out 'deconvolution' at 4th stations as well for the tarvl and amp
modname = 'timeoff_plfit_4thsta_160sps.mat';
planefit = load(strcat(rstpath, '/MAPS/',modname));

pred4diff = [];
predoff14 = [];
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
  offmax = round(1.0*rmse);
        
  [sigdecon(:,ista),pred,res,dresit,mfitit,ampit{ista},fighdl] = ...
    iterdecon_4thsta(sig,wlet,[],rcc1i(:,ista-3),[],...
    dt,twlet,impindep,stas(ista,:),off1i(ista),[],offmax,...
    fpltit,fpltend,fpltchk);
        
  ampiti = ampit{ista};
  impindep(:,9+(ista-4)*2+1) = ampiti(:,1);
  impindep(:,9+(ista-3)*2) = ampiti(:,2);
  if ista == nsta
    pred4diff(:,ista-3) = ampiti(:,end-1);  %difference between found peak and predicted arrival
    predoff14(:,ista-3) = ampiti(:,end);  %predicted off14
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
         
%       %%%plot the predicted tarvl and amp of decon impulses at 4th sta vs. waveform
%       [f] = plt_deconpk_sigpk_comp_4thsta(sigsta(:,4:nsta),stas(4:nsta,:),...
%         impindep(:,10:end),ppkindep4th,npkindep4th,greenf(:,4:nsta));
      
%%% further ELIMINATE sources that fail the check at 4th stations
trust4th = 4; % trust KLNB the most among all 4th stations
indremove = find(impindep(:,9+(trust4th-4)*2+1)==0 & impindep(:,9+(trust4th-3)*2)==0);
pred4difftr = pred4diff(setdiff(1:size(pred4diff,1),indremove),trust4th-3);
predoff14tr = predoff14(setdiff(1:size(predoff14,1),indremove),trust4th-3);
impindep(indremove,:) = [];
ppkindep(indremove, :) = [];
npkindep(indremove, :) = [];
impindepst = sortrows(impindep,1);

%%%if you want to plot the deconvolved sources
pltsrcflag3 = 1;
if pltsrcflag3
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
  title(ax,'Checkd at 4th stas');

end
% keyboard

%%
%off14 prediction for decon srcs using plane fit model WITH alignment
%'off14pred' is the same as 'predoff14tr' but in diff order, since 'predoff14tr' is computed from 
%unsorted sources 'impindep'
[~,off14pred] = pred_tarvl_at4thsta(stas(trust4th,:),impindepst(:,7),impindepst(:,8),...
  impindepst(:,1),0);

%source arrival prediction from plane fitting, including calibrating the arrival prediction in the
%context of prealigned signal, note sign is '+' 
[~,off14gt] = pred_tarvl_at4thsta(stas(trust4th,:),impgtst(:,7),impgtst(:,8),impgtst(:,1),0);

%off14 computed from actually-matched arrivals at stas 1 and 4 from decon
off14 = impindepst(:,1)-impindepst(:,9+(trust4th-4)*2+1)+repmat(off1i(trust4th),size(impindepst,1),1); %after prealignment

% f = initfig(12,5,1,1); %initialize fig
% ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');  
% p1=scatter3(ax,synsrc(:,2),synsrc(:,3),off14gt-off1i(trust4th),15,[.5 .5 .5],'filled',...
%   'MarkerEdgeColor','k'); %correct ground truth for record alignment
% p2=scatter3(ax,impindepst(:,7),impindepst(:,8),off14pred,15,[.5 .5 .5],'filled',...
%   'MarkerEdgeColor','r');
% legend(ax,[p1 p2],'ground truth', 'empirical pred for decon srcs');
% xlabel(ax,sprintf('off12 at %d sps',sps));
% ylabel(ax,sprintf('off13 at %d sps',sps));
% zlabel(ax,sprintf('off14 at %d sps',sps));
% title(ax,stas(trust4th,:));
% view(ax, 45, 10);
% ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');  
% p1=scatter3(ax,impindepst(:,7),impindepst(:,8),off14pred,15,[.5 .5 .5],'filled','MarkerEdgeColor','r');
% p2=scatter3(ax,impindepst(:,7),impindepst(:,8),off14,15,[.5 .5 .5],'filled','MarkerEdgeColor','k');
% legend(ax,[p1 p2],'empirical pred for decon srcs', 'actually matched');
% xlabel(ax,sprintf('off12 at %d sps',sps));
% ylabel(ax,sprintf('off13 at %d sps',sps));
% zlabel(ax,sprintf('off14 at %d sps',sps));
% view(ax, 45, 10);

f = initfig(12,5,1,3); %initialize fig
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');  
histogram(ax,off14pred-off14gt);
plot(ax,[-offmax -offmax],ax.YLim,'k--');
plot(ax,[offmax offmax],ax.YLim,'k--');
xlabel(ax,'diff. in plane-fit off14 between decon and ground truth');
ylabel(ax,'count');
ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); 
histogram(ax,pred4difftr);
plot(ax,[-offmax -offmax],ax.YLim,'k--');
plot(ax,[offmax offmax],ax.YLim,'k--');
xlabel(ax,'diff. in 4th arrival between pred and decon');
ylabel(ax,'count');
ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); 
histogram(ax,off14pred-off14);
plot(ax,[-offmax -offmax],ax.YLim,'k--');
plot(ax,[offmax offmax],ax.YLim,'k--');
xlabel(ax,'diff. in off14 between plane-fit and decon');
ylabel(ax,'count');


