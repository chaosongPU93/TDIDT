% testsyntheticRoo.m
%
% This script is to generate synthetic seismograms from a restricted region, both the optimal and
% orthogonal components using the corresponding templates. After the generation of a long record, we
% like to bin the record using a moving window based on the amplitude, obtain the amplitude spectrum
% of each bin for both components and obtain the amp ratio between the two. This may provide some
% insights in comparing with Allan's plot from the real data

%% Initialization
clear
clc
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');

[scrsz, res] = pixelperinch(1);

%% prepare templates (Green's functions)
ccstack = [];
sps = 160;
templensec = 60;

fam = '002';
disp(fam);

stas=['PGC ';
      'SSIB';
      'SILB';
      ];
nsta=size(stas,1);

%optimal components
for ista = 1: nsta
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
        num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao_catnew');
    ccstack(:,ista) = load(fname);
end
STA = ccstack;
%orthogonal components
for ista = 1: nsta
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
        num2str(templensec), 's_', 'BBCCS_', 'ort_Nof_Non_Chao_catnew');
    ccstackort(:,ista) = load(fname);
end
STAort = ccstackort;

for ista=1:nsta
    STA(:,ista)=Bandpass(STA(:,ista),sps,0.1,15,2,2,'butter');   % change 'bandpass' to 'Bandpass'
    STAort(:,ista)=Bandpass(STAort(:,ista),sps,0.1,15,2,2,'butter'); 
end

%%%The below aligns the templates by x-correlation
[maxses,imaxses]=max(STA,[],1);
[minses,iminses]=min(STA,[],1);
spread=maxses-minses;
% zcrosses=round(0.5*(imaxses+iminses));  % rough, assuming symmetry, Chao 2021/07/16
%automatically find the zero-crossings
zcrosses = zeros(nsta,1);
for ista = 1:nsta
    seg = STA(iminses(ista): imaxses(ista),ista);  % for zero-crossing timing, only use the main station
    [~,zcrosses(ista)] = min(abs(seg));
    zcrosses(ista) = zcrosses(ista)-1+iminses(ista);  % convert to global index
end
%now we want to cut a segment around the zero-crossing at each station
sampbef=6*sps;
sampaft=10*sps;
is=zcrosses-sampbef;
ie=zcrosses+sampaft;
for ista=1:nsta
    STAtmp(:,ista)=STA(is(ista):ie(ista),ista);  % this means templates are 'aligned' at zero-crossings
    STAtmport(:,ista)=STAort(is(ista):ie(ista),ista);
end
%x-correlation independently between each station pair 
mshiftadd=10*sps/40;
tempxc(:,1)=xcorr(STAtmp(:,2),STAtmp(:,3),mshiftadd,'coeff');
tempxc(:,2)=xcorr(STAtmp(:,1),STAtmp(:,2),mshiftadd,'coeff'); %shift STAtmp(3,:) to right for positive values
tempxc(:,3)=xcorr(STAtmp(:,1),STAtmp(:,3),mshiftadd,'coeff'); %shift STAtmp(2,:) to right for positive values
[~,imax]=max(tempxc,[],1);
imax=imax-(mshiftadd+1); %This would produce a slightly different shift, if filtered seisms were used.
imax(2)-imax(3)+imax(1)   %enclosed if it equals to 0
for ista=2:nsta
    STAtmp(mshiftadd+1:end-(mshiftadd+1),ista)=STAtmp(mshiftadd+1-imax(ista):end-(mshiftadd+1)-imax(ista),ista);
    STAtmport(mshiftadd+1:end-(mshiftadd+1),ista)=STAtmport(mshiftadd+1-imax(ista):end-(mshiftadd+1)-imax(ista),ista);
end
for ista=1:nsta
    STAtmp(:,ista)=STAtmp(:,ista)/spread(ista); % now templates are 'aligned' indeoendently by x-corr wrt. sta 1
    STAtmport(:,ista)=STAtmport(:,ista)/spread(ista); % now templates are 'aligned' indeoendently by x-corr wrt. sta 1
end
%%%The above aligns the templates by x-correlation

%%%detrend, taper and bandpass templates
tmpwlet = STAtmp; % no bandpass
tmpwletf = STAtmp;  % bandpassed version
tmpwletort = STAtmport;
tmpwletfort = STAtmport;
fractap = sps/size(tmpwlet,1);
for ista = 1: nsta
  %romve mean, linear trend of template
  tmpwlet(:,ista) = detrend(tmpwlet(:,ista));
  %and taper with tukeywin, which is actually a tapered cosine window
  w = tukeywin(size(tmpwlet(:,ista),1),fractap);
  tmpwlet(:,ista) = w.* tmpwlet(:,ista);
  %detrend again for caution
  tmpwlet(:,ista)=detrend(tmpwlet(:,ista));
  %filter the template
  hiwlet=18;
  lowlet=1.8;
  tmpwletf(:,ista) = Bandpass(tmpwlet(:,ista), sps, lowlet, hiwlet, 2, 2, 'butter');
  %detrend again for caution
  tmpwletf(:,ista)=detrend(tmpwletf(:,ista));
  
  tmpwletort(:,ista) = detrend(tmpwletort(:,ista));
  w = tukeywin(size(tmpwletort(:,ista),1),fractap);
  tmpwletort(:,ista) = w.* tmpwletort(:,ista);
  tmpwletort(:,ista)=detrend(tmpwletort(:,ista));
  hiwlet=18;
  lowlet=1.8;
  tmpwletfort(:,ista) = Bandpass(tmpwletort(:,ista), sps, lowlet, hiwlet, 2, 2, 'butter');
  tmpwletfort(:,ista)=detrend(tmpwletfort(:,ista));
end

%%%constrained CC, so that only 2 offsets are independent
ccmid = round(size(tmpwletf,1)/2);
ccwlen = 10*sps;
loffmax = 5*sps/40;
ccmin = 0.01;  % depending on the length of trace, cc could be very low
iup = 1;    % times of upsampling
[off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(tmpwletf',ccmid,...
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
greenort = zeros(greenlen,nsta); 
greenfort = zeros(greenlen,nsta); 

for ista = 1: nsta
  %cut according to the zero-crossing and the time shift from the constrained CC
  green(:,ista) = tmpwlet(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenf(:,ista) = tmpwletf(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  %detrend again for caution
  green(:,ista)=detrend(green(:,ista));
  greenf(:,ista)=detrend(greenf(:,ista));
  %normalize by max amp
  green(:,ista)=green(:,ista)/max(abs(green(:,ista)));    % normalize
  greenf(:,ista)=greenf(:,ista)/max(abs(green(:,ista)));    % normalize
  
  greenort(:,ista) = tmpwletort(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenfort(:,ista) = tmpwletfort(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenort(:,ista)=detrend(greenort(:,ista));
  greenfort(:,ista)=detrend(greenfort(:,ista));
  greenort(:,ista)=greenort(:,ista)/max(abs(green(:,ista)));    % normalize
  greenfort(:,ista)=greenfort(:,ista)/max(abs(green(:,ista)));    % normalize

  %re-find the zero-crossing as the template length has changed
  [~,imin] = min(greenf(:,ista));
  [~,imax] = max(greenf(:,ista));
  [~,zcrosses(ista)] = min(abs(greenf(imin:imax,ista)));
  zcrosses(ista) = zcrosses(ista)+imin-1;
end
%the following is just a check, because now the templates must be best aligned 
ccmid = round(size(greenf,1)/2);
ccwlen = 4*sps;
loffmax = 5*sps/40;
ccmin = 0.01;  % depending on the length of trace, cc could be very low
iup = 1;    % times of upsampling
[off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(greenf',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
if ~(off12con==0 && off13con==0)
  disp('Filtered templates are NOT best aligned');
end

%%%plot the unfiltered and filtered templates
mean(green,1)
figure
subplot(2,1,1)
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

subplot(2,1,2)
hold on
plot(greenf(:,1),'r')
plot(greenf(:,2),'b')
plot(greenf(:,3),'k')
text(0.95,0.9,sprintf('%.1f-%.1f Hz',lowlet,hiwlet),'Units','normalized','HorizontalAlignment',...
  'right');
mx=max([abs(green(:,1)); abs(green(:,2))]);

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
[ccmax23,imax23]=max(xc23)
[ccmax13,imax13]=max(xc13)
[ccmax12,imax12]=max(xc12)

%orthogonal components
figure
subplot(2,1,1)
hold on
plot(greenort(:,1),'r')
plot(greenort(:,2),'b')
plot(greenort(:,3),'k')
text(0.95,0.9,'Raw orthogonal','Units','normalized','HorizontalAlignment',...
  'right');
xlim([0 greenlen])
ylim([-mx mx])
box on

subplot(2,1,2)
hold on
plot(greenfort(:,1),'r')
plot(greenfort(:,2),'b')
plot(greenfort(:,3),'k')
text(0.95,0.9,sprintf('%.1f-%.1f Hz orthogonal',lowlet,hiwlet),'Units','normalized',...
  'HorizontalAlignment','right');
xlim([0 greenlen])
ylim([-mx mx])
box on

%% load the true 4-s tremor detections, obtain its PDF that will be used to generate synthetic sources
freqflag='hf';  % flag to indicate whether to do hf or lf;
FLAG = 'PGC'; % detector 
fam = '002';   % family number
spsdect = 40;
iup = 4;  % upsample 4 times to 160 sps
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
winlen=winlensec*spsdect;      % length in smaples
PREFIX = strcat(fam,'.up.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),'.',...
                        int2str(npo),int2str(npa),'.ms', int2str(mshift));
hftime = load(strcat(rstpath, '/MAPS/tdectimeori_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',...
            num2str(winlen/spsdect),'s',num2str(spsdect),'sps','4add'));
%format as follows:
%%% 34+4*nstanew cols, if 4 new stas, then will be 50 cols
%%% UPDATED at 2021/06/23
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

%%%Interpolate from existing grids for the locations of sources
%round to the integer sample at the 'iup' times the original sampling rate 40 sps
offint = round(hftime(:, 1:2)*iup);
          
%convert the time offset to locations
% ftrans = 'interpArmb';
ftrans = 'interpArmbreloc';
% ftrans = 'interpchao';
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
%%% UPDATED at 2021/06/23
%%%   dx,dy,lon,lat,dep,ttrvl,off12,off13 (integer samples at upsampled sps)
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

%%%Note that the locations in the ellipse corresponds to integer time offsets at sps*iup, so we can
%%%avoid doing the transformation, and directly get the density and PDF in sample space, the benefit
%%%is that, in sample sapce, it is actucally an even grid with a spacing of 1 sample 
denoff= density_pixel(hfbnd(:,7),hfbnd(:,8));

%Obtain the PDF
epdfoff = denoff;
%now the 'area' of each bin is 1*1 sample
epdfoff(:,3) = epdfoff(:,3)/sum(epdfoff(:,3))/(1*1); % normalize to PDF, ==counts/N_total/area_of_bin

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
color=colorbar;
color.Label.String = 'Probability Density';
caxis([0 max(epdfoffpad(:,3))]);
xlabel(strcat({'Offset 12 (samples at '},num2str(sps),' sps)'));
ylabel(strcat({'Offset 13 (samples at '},num2str(sps),' sps)'));
box on; grid on;


%% synthetics generation
%%%Specify the amplitude-frequency (counts) distribution
%distr='PL'
%distr='LN'
distr='UN'  % uniform distribution
%distr='DI'
%distr='EX'

%Twin=3*3600+3+2*ceil(greenlen/sps); %Early and late portions will be deleted. Twin includes only the early portion. 
Twin=2*3600+3+2*ceil(greenlen/sps); %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
%Twin=4*3600+3+2*ceil(greenlen/sps); %For DI. 
%Twin=6*3600+3+2*ceil(greenlen/sps); %For DI. 
%Twin=1*3600+3+2*ceil(greenlen/sps); %For DI. 
winlen=Twin*sps+1;
skiplen=greenlen;
tdura = 0.5;  % duration from the broadband template
% tdura = 0.75; % for bandpassed, can be ~0.3s or 0.7s depending on the definition of 'duration'
satn=1/tdura*Twin   % if just saturated, how many templates can be fit in? a single peak is ~20 samples wide; maybe a little less (at 100 sps).
    %Twin is window duration in seconds. Events can fall within Twin of
    %the start, but the synthetics will go to Twin*(sample rate)+Greenlen to
    %avoid checking for subscript overrrun.  When "synths" are written from "synth" in the
    %subroutine, Greenlen from the start and end will not be written.
fracelsew=0; %0.25 %0.5; %0.6; %The ratio of elsewhere events to local events.  Zero for Discrete Ide.
nsat=[1 5];  % times of saturation 
writes=round(nsat*satn) %how many templates to throw in, under different degrees of saturation
rng('default');
%%%Here the noise is 'uniform' in time!
% synth=2.e-7*(rand(nsta,winlen+greenlen+2*mjig)-0.5); %+2*mjig a little extra, for jiggering 2nd & 3rd stations.
mampnoi = 1e-3;
synth=mampnoi*2*(rand(winlen+greenlen+2*10,nsta)-0.5); %+2*10 a little extra, for jiggering 2nd & 3rd stations.% for ista=1:nsta
% synth(ista,:)=1.e-7*sin(1:winlen+greenlen+2*mjig); %+2*mjig a little extra, for jiggering 2nd & 3rd stations.
synthort=mampnoi*2*(rand(winlen+greenlen+2*10,nsta)-0.5); %+2*10 a little extra, for jiggering 2nd & 3rd stations.% for ista=1:nsta

nouts=length(writes);
seed=round(writes(1)/5e3); %for random number generator

if strcmp(distr,'PL') || strcmp(distr,'UN') %b>150 for uniform
    b=999. %>150 for uniform size distribution
    [synths,mommax,sources,greensts]=synthgen(writes,winlen,skiplen,synth,greenf,b,xgrid,...
      ygrid,epdfoffgrid,sps,fracelsew,seed);
    
    [synthsort,~,~,~]=synthgen(writes,winlen,skiplen,synthort,greenfort,b,xgrid,...
      ygrid,epdfoffgrid,sps,fracelsew,seed);
end

%%%a simple plot the synthetics
figure
nrow = length(writes)+1;
ncol = 1;
subplot(nrow,ncol,1)
hold on
tmp = synth;
plot(tmp(:,1),'r');
plot(tmp(:,2),'b');
plot(tmp(:,3),'k');
xlim([0 40*sps]);
axranexp(gca,6,20);

for i = 1: length(writes)
subplot(nrow,ncol,i+1)
hold on
tmp = synths(:,:,i);
plot(tmp(:,1),'r');
plot(tmp(:,2),'b');
plot(tmp(:,3),'k');
text(0.95,0.9,sprintf('%.1f x saturation',nsat(i)),'Units','normalized','HorizontalAlignment',...
  'right');
xlim([0 40*sps]);
axranexp(gca,6,20);
end
xlabel(sprintf('Samples at %d Hz',sps),'FontSize',12);

% %% extract and validate the added impulses of template 
% insat = 1;  %which saturation to look at
% wlensec = 50; %how long the window to test
% bufsec = 1; %need some buffer window for later CC alignment
% buffer = bufsec*sps;
% wlensecb = wlensec+bufsec;  
% 
% sigstab = [];  %target window of synthetic signals at all stations
% sigpnstab = [];  %target window of synthetic signals plus the noise
% sigconvstab = [];  %target window of reproduced synthetic signals using convolution
% noistab = [];  %target window of noise
% 
% srcpairstab = cell(nsta,1);  %sources whose zero-crossing with time range, [indtarvl rnoff12 rnoff13 amp]
% 
% for ista = 1: nsta
%   % ista = 3;
%   
%   tgreen = zcrosses(ista);  % anchor time of templates, here choose as zero-crossing time
%   indstsig = 1;  % starting index of the simulated signal to test
%   indedsig = wlensecb*sps+indstsig-1; % ending index of the simulated signal to test
%   source = sources(1:size(greensts{insat}{1},1), :);  % [indtori indttrvl indtarvl rnoff12 rnoff13 amp];
%   if ista == 1
%     greenst = greensts{insat}{1}; % the starting index of each added template, context of full length
%   else
%     %ind - rnoff is the arrival time in index at sta 2 and 3
%     %note the sign here, if off12 >0, move 2 to the right to align with 1, meaning that 2 is
%     %earlier than 1, ie., tarvl2 < tarvl1. Be consistent all the time
%     greenst = greensts{insat}{1}-source(:, 2+ista); 
%   end
%   
%   %%%you don't need all impulses, only some of them contribute to the length of truncated record
%   %%%cut out the green indice and sources that contribute
%   %%%Note that some zero-crossings might appear later than 'indedsig', but the corresponding start
%   %%%index of the green is actually smaller than 'indedsig', ie. their 'greenst' <= 'indedsig'
%   induse = find(greenst>=skiplen-greenlen+indstsig & greenst<=indedsig+skiplen-1);
%   greenst = greenst(induse);
%   source = source(induse,:);
%   
%   greenzc = greenst+tgreen; % index of approximate zero-crossing
%   source(:, 3) = source(:, 3)+tgreen;  % now 'greenzc' should be the same as 3rd col of 'source'
%   impamp = zeros(max(greenzc)+20,1);
%   for i = 1: length(greenzc)
%     impamp(greenzc(i)) = impamp(greenzc(i))+source(i, 6);
%   end
%   
%   %note that some length of the full simulation 'skiplen' was skipped in 'synthgen.m'
%   sigpntemp = squeeze(synths(:,ista,insat));  % simulated signal with white noise
%   noitemp = synth(skiplen:winlen,ista);   % account for the skipped part
%   sigtemp = sigpntemp-noitemp; % subtract the white noise
%   sigpnstab(:,ista) = sigpntemp(indstsig:indedsig); % now focus on the part of interest
%   sigstab(:,ista) = sigtemp(indstsig:indedsig);
%   noistab(:,ista) = noitemp(indstsig:indedsig);
%   
%   sigconvtmp = conv(greenf(:,ista),impamp,'full');
%   indtrc = tgreen+skiplen;  % starting index for truncation
%   sigconvstab(:,ista) = sigconvtmp(indtrc+indstsig-1:indedsig+indtrc-1);  % cut accordingly
%   
%   %%%Can we reproduce the synthetics with truncated convolution? YES
%   figure
%   subplot(311)
%   plot(greenf(:,ista),'r');
%   xlim([0 greenlen]);
%   legend('Template');
%   title('Reproduce synthetics with truncated convolution');
%   subplot(312)
%   imptemp = find(impamp>0);
%   p1=stem(imptemp-skiplen, impamp(imptemp),'b','MarkerSize',4); hold on;
%   % p1=stem((1:size(impamp))-skiplen, impamp,'b','MarkerSize',4); hold on;
%   ax = gca;
%   plot([indstsig indstsig],ax.YLim,'--','color',[.5 .5 .5]);
%   plot([indedsig indedsig],ax.YLim,'--','color',[.5 .5 .5]);
%   legend(p1,'Synthetic random impulses');
%   subplot(313);
%   plot(sigstab(:,ista),'b'); hold on
%   plot(sigconvstab(:,ista),'k');
%   legend('Truncated synthetic signal','Truncated signal from convolution');
%   
%   %now 'source' contain impulse triplets that contribute to the trace segment of interest, but we
%   %want to get ones whose zero-crossing falls into the window.
%   source(:,3) = source(:,3)-skiplen;  % index after skipping 
%   srcpairb = source(source(:,3)>=indstsig & source(:,3)<=indedsig, 3:6);
%   srcpairb = sortrows(srcpairb,1,'ascend');
%   srcpairstab{ista} = srcpairb; %ideally for each station this should be the same, but coincidence is possible  
%   
% end    

%% bin the synthetics upon amp, using a half-overlapping win of 18 s
bwlensec = 18;     % offsec = 3 was used in first-year report
bwoffsec = 9;        % window offset in sec, which is the step of a moving window
bwlen = bwlensec*sps;      % length in smaples
bwoff = bwoffsec*sps;      % offset in samples
insat = 2;
opt = squeeze(synths(:,:,insat));
ort = squeeze(synthsort(:,:,insat));
tracelen = size(opt,1);
%cut out time at both the start and end of simulation, useful also when you specify arrival time as
%uniform distribution
cutsec = 5;   
indst = cutsec*sps+1;
inded = tracelen-cutsec*sps;

wins = movingwins(indst,inded,bwlen,bwlen-bwoff);
nwin = size(wins,1);

%params for getting the spectrum
nfft = bwlen;
window = hann(bwlen);
Fs = sps;

ampopt = zeros(nwin,1);
psdopt = zeros(nfft/2+1,nwin);
psdort = zeros(nfft/2+1,nwin);

%%%obtain the amplitude estimate for each moving window
for i = 1: nwin
  win = wins(i,:);
  optseg = opt(win(1): win(2), :);  %(:,i)
  ortseg = ort(win(1): win(2), :);
  
  ampopt(i) = median(median(abs(optseg), 1)); % median of the segment, then median of 3 stations
  
  %compute the amplitude spectrum
  psdopts = zeros(nfft/2+1, nsta);
  for ista = 1: nsta
    [psdopts(:,ista),pft] = periodogram(optseg(:,ista),window,nfft,Fs);
  end
  psdopt(:,i) = median(sqrt(psdopts),2); %'sqrt' to get amp from power 
  
  psdorts = zeros(nfft/2+1, nsta);
  for ista = 1: nsta
    [psdorts(:,ista),~] = periodogram(ortseg(:,ista),window,nfft,Fs);
  end
  psdort(:,i) = median(sqrt(psdorts),2);
end

%bin windows based on the amplitude of the optimal component
nbin = 10;
%'ibin' indicates which bin each element of data belongs to
[ibin,binedge] = discretize(ampopt,nbin); 

figure
histogram(ampopt,nbin);
xlabel('Median abs. amplitude');
ylabel('Counts');


for i = 1: nbin
  ind = find(ibin==i);
  psdoptbin(:,i) = mean(psdopt(:,ind),2);   % mean, an average of the spectrum in the same bin
  psdortbin(:,i) = mean(psdort(:,ind),2);
  
  %ratio of optimal to orthogonal spectrum
  psdroobin(:,i) = psdoptbin(:,i)./psdortbin(:,i);  

end

%obtain the amplitude spectrum of the template
nfft = size(greenf,1);
window = hann(size(greenf,1));
Fs = sps;
psdoptgfs = zeros(nfft/2+1, nsta);
for ista = 1: nsta
    [psdoptgfs(:,ista),pftgf] = periodogram(greenf(:,ista),window,nfft,Fs);
end
psdoptgf = median(sqrt(psdoptgfs),2);

psdortgfs = zeros(nfft/2+1, nsta);
for ista = 1: nsta
    [psdortgfs(:,ista),~] = periodogram(greenfort(:,ista),window,nfft,Fs);
end
psdortgf = median(sqrt(psdortgfs),2);


%%%
figure
% subplot(131)
% loglog(pftgf,psdoptgf,'k-'); hold on
% loglog(pftgf,psdortgf,'k--');
% % axis equal
% xlim([1e-2, 1e2]);
% ylim([1e-5, 1e0]);

ax(1)=subplot(121);
color = jet(nbin);
for i = 1: nbin
  loglog(pft,psdoptbin(:,i),'-','Color',color(i,:)); hold on
end
for i = 1: nbin
  loglog(pft,psdortbin(:,i),'--','Color',color(i,:));
end
loglog(pftgf,psdoptgf,'k-','linew',1.5); hold on
loglog(pftgf,psdortgf,'k--','linew',1.5);
colormap(jet(nbin))
cbar=colorbar;
caxis([0 nbin]);
cbar.Label.String = 'Bin index';
% axis equal
xlim([1e-2, 1e2]);
ylim([1e-5, 1e0]);
% grid on
xlabel('Frequency (Hz)');
ylabel('Amplitude');

ax(2)=subplot(122);
for i = 1: nbin
  semilogx(pft,psdroobin(:,i),'-','Color',color(i,:)); hold on
end
colormap(jet(nbin))
cbar=colorbar;
caxis([0 nbin]);
cbar.Label.String = 'Bin index';
xlim([1e-2, 1e2]);
xlabel('Frequency (Hz)');
ylabel('Opt/ort');

supertit(ax([1 2]),sprintf('Saturation level: %.1f; wlen: %d s; BP templates: %.1f-%.1f Hz', ...
  nsat(insat), bwlensec, lowlet, hiwlet), 10);














