% function synthrccvssatur.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is similar in large part to 'synthshift_chao.m', especially 
% for the start. The goal is to generate synthetic seismograms from random 
% sources that are drawn from the real dat PDF from the 4-s tremor detection
% distribution, where each source is the SAME as the template. The origin 
% time of source is unifrom. The amplitude is identical of 1. 
%
% --By varying the saturation rate, we want to know the relationship 
% between rcc and saturation rate for the same window length after aligning 
% different stations. Also, for the same saturation rate, it is also possible
% to know the relation between rcc and window length for the same saturation.
%
% --Allan's recollection is the decreasing rate of rcc would be smaller wrt. 
% saturation rate. In other words, rcc would be ~constant beyond some saturation
% level.
% 
% --My feeling for rcc vs. win length is, median(rcc) would decrease with the 
% increase of window length, especially for real data where the saturation level
% might be time-variant so that a longer window might have a bigger change to 
% include a larger portion of 'noise' [noise from signal elsewhere or random 
% noise]. If this feeling is correct, then for a long window, maybe you should
% NOT use med(rcc) as the stopping criterion for iterdecon, otherwise you may 
% overfit the noise due to more iters. However, i am not sure how the synthetics
% would behave.
%
% --If you instead use a uniform random distribution for sources, then you 
% should use code 'analyze_synth'. 
%
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/05/04
% Last modified date:   2023/09/07
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
plt_templates(green,greenf,stas,[],[],lowlet,hiwlet,sps);


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

%% interpolate from existing grids for the locations of sources
%round to the integer sample at the 'iup' times the original sampling rate 
offint = round(hftime(:, 1:2)*iup);
          
%convert the time offset to locations
% ftrans = 'interpArmb';
% ftrans = 'interpArmbreloc';
ftrans = 'interpchao';
fplt = 0;
if isequal(ftrans,'interpchao')
  [loc, indinput] = off2space002(offint,sps,ftrans,fplt);
  % loc has 8 cols, format: dx,dy,lon,lat,dep,tori,off12,off13
  hfuse = hftime;
  hfuse = [loc hfuse];
elseif isequal(ftrans,'interpArmb') || isequal(ftrans,'interpArmbreloc')
  ind = find(offint(:,1)>=-23*iup & offint(:,1)<=25*iup & ...
             offint(:,2)>=-24*iup & offint(:,2)<=24*iup);
  [loc, indinput] = off2space002(offint(ind,:),sps,ftrans,fplt);
  % loc has 8 cols, format: dx,dy,lon,lat,dep,tori,off12,off13
  hfuse = hftime(ind,:);
  hfuse = [loc hfuse];
end

%outline a boundary region on top of the density map
cutout = 'ellipse';
x0 = 0.2; 
y0 = 0.2;
semia = 1.75;
semib = 1.25;
angrot = 45;
[xcut, ycut] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);

%get the detections inside the cutout boundary
bnd = [xcut ycut];
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

%%%Obtain the location density (count), bin by pixel
% dx = 0.05;
% dy = 0.05;
xran = [-4 4];
yran = [-4 4];
denloc= density_pixel(hfbnd(:,1),hfbnd(:,2));
%plot the density distribution
figure
scatter(denloc(:,1),denloc(:,2),6,denloc(:,3),'filled');
axis equal
axis([xran yran]);
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
c=colorbar;
c.Label.String = 'Probability Density';
caxis([0 max(epdfoffpad(:,3))]);
xlabel(strcat({'Offset 12 (samples at '},num2str(sps),' sps)'));
ylabel(strcat({'Offset 13 (samples at '},num2str(sps),' sps)'));
title('Data PDF in smaple space');
box on; grid on;


%% synthetics generation
%%%Specify the amplitude-frequency (counts) distribution
distr='UN'  % uniform distribution

Twin=600+3+2*ceil(greenlen/sps); %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
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
nsat=[0.1 0.2 0.4 1 2 4 10 20 40 100];  % times of saturation 
writes=round(nsat*satn) %how many templates to throw in, under different degrees of saturation
rng('default');
%%%Here the noise is 'uniform' in time!
% mampnoi = 1e-3;
% synth=mampnoi*2*(rand(winlen+greenlen+2*10,nsta)-0.5); %+2*10 a little extra, for jiggering 2nd & 3rd stations.% for ista=1:nsta
%%%Below use the Gaussian white noise
noistd = 5e-2;
synth=noistd*(randn(winlen+greenlen+2*10,nsta)); %+2*10 a little extra, for jiggering 2nd & 3rd stations.% for ista=1:nsta

nouts=length(writes);
seed=round(writes(1)/5e3); %for random number generator

if strcmp(distr,'PL') || strcmp(distr,'UN') %b>150 for uniform
    b=999. %>150 for uniform size distribution
    [synths,mommax,sources,greensts]=synthgen(writes,winlen,skiplen,synth,greenf,b,xgrid,...
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
xlim([0 45*sps]);
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
xlim([0 45*sps]);
axranexp(gca,6,20);
end
xlabel(sprintf('Samples at %d Hz',sps),'FontSize',12);


%%
%%%Is the median median. abs amplitude proportional to the square root of number of templates?
%%%The abs amp is similar to envelope
for i = 1: length(writes)
  tmp = squeeze(synths(:,:,i));
  ampmed(i) = median(median(abs(tmp),1));
  ntempsqrt(i) = sqrt(writes(i));
  ntempsqrt(i)/ampmed(i)
end
fttpfree = fittype( @(a,x) a*x);
[fitobj,~,~] = fit(ntempsqrt(:),ampmed(:),fttpfree,'Robust','Bisquare',...
  'StartPoint',[1]);
coef = coeffvalues(fitobj);
coefci = confint(fitobj);
xplt = min(ntempsqrt): 0.01: max(ntempsqrt);
yfit = coef(1)*xplt ;
ylow = coefci(1,1)*xplt;
yupp = coefci(2,1)*xplt;

figure
% ax=gca;ax,
scatter(ntempsqrt,ampmed,20,'b','filled'); hold on
plot(xplt,yfit,'k-','linew',2);
patch([xplt fliplr(xplt)],[ylow fliplr(yupp)],'k','Facealpha',0.2,'edgecolor','none');
xlabel('Square root of number of tmeplates');
ylabel('Med. of med. of abs amplitude');
box on; grid on;

    
%% median(rcc) or (cc) VS. saturation rate
wlensectry = (15:30:300)';
medrcc = zeros(length(wlensectry), length(nsat));
mcc = zeros(length(wlensectry), length(nsat));
nsrc = zeros(length(wlensectry), length(nsat));

for insat = 1: length(nsat)
  % insat = 4;  %which saturation to look at
  disp(nsat(insat));
  for iwlen = 1: length(wlensectry)
    wlensec = wlensectry(iwlen);
%     iwlen = 1;
%     wlensec = 47; %this is the median win length of from burst grouping result in 'tremorbursts002_4s.m'
    bufsec = 1; %need some buffer window for later CC alignment
    buffer = bufsec*sps;
    wlensecb = wlensec+bufsec;
    
    %generate half-overlapping windows of the same length
    indst = 2*sps/2;
    inded = winlen-skiplen-sps/2-buffer;
    [windows] = movingwins(indst,inded,wlensecb*sps,wlensecb*sps/2,0);
    nwin =  size(windows,1);
    
    rcctmp = zeros(nwin,1);  % store the result for each overlapping window, will get averaged later 
    mcctmp = zeros(nwin,1);
    nsrctmp = zeros(nwin,1);
    
    for iwin = 1: nwin
%     indstsig = 1;  % starting index of the simulated signal to test
%     indedsig = wlensecb*sps+indstsig-1; % ending index of the simulated signal to test
    indstsig = windows(iwin, 1);  % starting index of the simulated signal to test
    indedsig = windows(iwin, 2); % ending index of the simulated signal to test
    
    sigstab = [];  %target window of synthetic signals at all stations
    sigpnstab = [];  %target window of synthetic signals plus the noise
    noistab = [];  %target window of noise
    
    srcpairstab = cell(nsta,1);  %sources whose zero-crossing with time range, [indtarvl rnoff12 rnoff13 amp]
    
    for ista = 1: nsta
      % ista = 3;
      
      tgreen = zcrosses(ista);  % anchor time of templates, here choose as zero-crossing time
      source = sources(1:size(greensts{insat}{1},1), :);  % [indtori indttrvl indtarvl rnoff12 rnoff13 amp];
      if ista == 1
        greenst = greensts{insat}{1}; % the starting index of each added template, context of full length
      else
        %ind - rnoff is the arrival time in index at sta 2 and 3
        %note the sign here, if off12 >0, move 2 to the right to align with 1, meaning that 2 is
        %earlier than 1, ie., tarvl2 < tarvl1. Be consistent all the time
        greenst = greensts{insat}{1}-source(:, 2+ista);
      end
      
      %%%you don't need all impulses, only some of them contribute to the length of truncated record
      %%%cut out the green indice and sources that contribute
      %%%Note that some zero-crossings might appear later than 'indedsig', but the corresponding start
      %%%index of the green is actually smaller than 'indedsig', ie. their 'greenst' <= 'indedsig'
      induse = find(greenst>=skiplen-greenlen+indstsig & greenst<=indedsig+skiplen-1);
      greenst = greenst(induse);
      source = source(induse,:);
      
      greenzc = greenst+tgreen; % index of approximate zero-crossing
      source(:, 3) = source(:, 3)+tgreen;  % now 'greenzc' should be the same as 3rd col of 'source'
      impamp = zeros(max(greenzc)+20,1);
      for i = 1: length(greenzc)
        impamp(greenzc(i)) = impamp(greenzc(i))+source(i, 6);
      end
      
      %note that some length of the full simulation 'skiplen' was skipped in 'synthgen.m'
      sigpntemp = squeeze(synths(:,ista,insat));  % simulated signal with white noise
      noitemp = synth(skiplen:winlen,ista);   % account for the skipped part
      sigtemp = sigpntemp-noitemp; % subtract the white noise
      sigpnstab(:,ista) = sigpntemp(indstsig:indedsig); % now focus on the part of interest
      sigstab(:,ista) = sigtemp(indstsig:indedsig);
      noistab(:,ista) = noitemp(indstsig:indedsig);
      
      %now 'source' contain impulse triplets that contribute to the trace segment of interest, but we
      %want to get ones whose zero-crossing falls into the window.
      source(:,3) = source(:,3)-skiplen;  % index after skipping
      srcpairb = source(source(:,3)>=indstsig & source(:,3)<=indedsig, 3:6);
      srcpairb = sortrows(srcpairb,1,'ascend');
      srcpairstab{ista} = srcpairb; %ideally for each station this should be the same, but coincidence is possible
      
    end
    
    %notify if sources are the same at diff stations
    if ~isequaln(srcpairstab{1,1},srcpairstab{2,1}) || ~isequaln(srcpairstab{1,1},srcpairstab{3,1})
      disp('Extracted sources at different stations are not the same, re-check!');
    end
    
    
    %% Best alignment for the testing window
    ccmid = round(size(sigpnstab,1)/2);
    ccwlen = round(size(sigpnstab,1)*0.8);
    loffmax = 5*sps/40;
    ccmin = 0.01;  % depending on the length of trace, cc could be very low
    iup = 1;    % times of upsampling
    mshiftadd=10*sps/40;
    [off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(sigpnstab',ccmid,...
      ccwlen,mshiftadd,loffmax,ccmin,iup);
    if off12con == mshiftadd+1 && off13con == mshiftadd+1
        off12con = 0;
        off13con = 0;
        fprintf('%d / %d could not be better aligned \n',iwin, nwin);
    end
    off1i = zeros(nsta,1);
    off1i(2) = round(off12con);
    off1i(3) = round(off13con);
    
    %align the signals, noises, etc
    sigpnsta = zeros(wlensec*sps, nsta);
    noista = zeros(wlensec*sps, nsta);
    sigsta = zeros(wlensec*sps, nsta);
    for ista = 1: nsta
      sigpnsta(:,ista) = sigpnstab(round(buffer/2)+1-off1i(ista): end-round(buffer/2)-off1i(ista), ista);
      noista(:,ista) = noistab(round(buffer/2)+1-off1i(ista): end-round(buffer/2)-off1i(ista), ista);
      sigsta(:,ista) = sigstab(round(buffer/2)+1-off1i(ista): end-round(buffer/2)-off1i(ista), ista); % sta 2
    end
    
    srcpairb = srcpairstab{1,1};
    %store the info of synthetic impulses, same format as in paired deconvolution,
    %9 cols, [ind1 amp1 ind2 amp2 ind3 amp3 off12 off13 off23]
    srcpair = zeros(size(srcpairb,1), 9);
    srcpair(:,1) = srcpairb(:,1)-round(buffer/2);  % cut out the buffer
    srcpair(:,[2 4 6]) = repmat(srcpairb(:,4), 1,3);
    %the indices of synthetic impulses need to be shifted too
    srcpair(:,3) = srcpair(:,1)-srcpairb(:,2)-off1i(2); % note the sign is consistent!
    srcpair(:,5) = srcpair(:,1)-srcpairb(:,3)-off1i(3);
    srcpair(:,7) = srcpairb(:,2);  % no need to shift offset, because they are 'true' offset
    srcpair(:,8) = srcpairb(:,3);
    srcpair(:,9) = srcpair(:,8)-srcpair(:,7);  % off23, == off13 - off12 == tarvl2 - tarvl3
    
    %compute running CC
    mwlen=sps/2;
    % mwlen=sps;
    [ircc,rcc12] = RunningCC(sigpnsta(:,1), sigpnsta(:,2), mwlen);
    [~,rcc13] = RunningCC(sigpnsta(:,1), sigpnsta(:,3), mwlen);
    [~,rcc23] = RunningCC(sigpnsta(:,2), sigpnsta(:,3), mwlen);
    rcc = (rcc12+rcc13+rcc23)/3;  %average
    rcctmp(iwin) = median(rcc);
    
    %compute the median absolute amplitude and envelope of the same moving window
    %for the moving window at the same station, sensable to use median
    [ir,ramp1,renv1] = Runningampenv(sigpnsta(:,1),mwlen,mwlen-1,'median');
    [~,ramp2,renv2] = Runningampenv(sigpnsta(:,2),mwlen,mwlen-1,'median');
    [~,ramp3,renv3] = Runningampenv(sigpnsta(:,3),mwlen,mwlen-1,'median');
    %looks like using the amplitude and envelope are pretty similar
    %maybe better to mean, otherwise stations PGC will be underrepresented due to a larger amp
    %variation
    ramp = mean([ramp1 ramp2 ramp3], 2);  % use mean or median??
    renv = mean([renv1 renv2 renv3], 2);
    
    %overall average max CC based on the current best alignment
    mcc12 = sum(sigpnsta(:,1).*sigpnsta(:,2))./ ...
      (sqrt(sum(sigpnsta(:,1).^2)).*sqrt(sum(sigpnsta(:,2).^2)));
    mcc13 = sum(sigpnsta(:,1).*sigpnsta(:,3))./ ...
      (sqrt(sum(sigpnsta(:,1).^2)).*sqrt(sum(sigpnsta(:,3).^2)));
    mcc23 = sum(sigpnsta(:,2).*sigpnsta(:,3))./ ...
      (sqrt(sum(sigpnsta(:,2).^2)).*sqrt(sum(sigpnsta(:,3).^2)));
    mcctmp(iwin) = (mcc12+mcc13+mcc23)/3;
    
    %number of sources actually of interest
    nsrctmp(iwin) = size(srcpair,1);
    
%     %distribution of synthetic sources, ground truth
%     spsscale = sps/40;
%     loff_max = 4*spsscale;  % maximum allowed shift is 4 samples at 40 Hz, it is also true from the synthetic distribution
%     offxran = [-loff_max loff_max];
%     offyran = [-loff_max loff_max];
%     lsig = size(sigpnsta,1);
%     cran = [0 lsig];
%     [f] = plt_decon_imp_scatter(srcpair,offxran,offyran,cran,sps,50,'mean','tori');%
%     scatter(gca,off1i(2),off1i(3),20,'ks','filled','MarkerEdgeColor','k');
%     title(gca,sprintf('Synthetic sources, using data of %d Hz',sps));
    
%     %%% scatter of rela locations, and account for prealignment offset
%     xran = [-3 3];
%     yran = [-3 3];
%     cran = [0 lsig];
%     [f] = plt_decon_imp_scatter_space(srcpair,xran,yran,cran,offxran,offyran,sps,50,ftrans,'mean','tori');
%     plot(gca,xcut,ycut,'k-','linew',2);
%     title(gca,sprintf('Synthetic sources, saturation level: %.1f',nsat(insat)));
    end
    
    medrcc(iwlen, insat) = mean(rcctmp);
    mcc(iwlen, insat) = mean(mcctmp);
    nsrc(iwlen, insat) = mean(nsrctmp);
  end
end

%%
%%%For each window length, plot of median(rcc)/mcc  vs.  saturation rate
figure
color = jet(length(wlensectry));
ax(1) = subplot(221); hold on
for i = 1: length(wlensectry)
  plot(nsrc(i,:),mcc(i,:),'o-','markersize',4,'Color',color(i,:));
end
legend(num2str(wlensectry),'Location','southeast');
xlabel('Number of sources inside window');
ylabel('Overall CC');
box on; grid on

ax(2) = subplot(222); hold on
for i = 1: length(wlensectry)
  plot(nsrc(i,:),medrcc(i,:),'o-','markersize',4,'Color',color(i,:));
end
% legend(num2str(wlensectry),'Location','southeast');
xlabel('Number of sources inside window');
ylabel('Median of running CC');
box on; grid on

ax(3) = subplot(223); hold on
for i = 1: length(wlensectry)
  plot(nsat,mcc(i,:),'o-','markersize',4,'Color',color(i,:));
end
xlabel('Saturation level');
ylabel('Overall CC');
box on; grid on

ax(4) = subplot(224); hold on
for i = 1: length(wlensectry)
  plot(nsat,medrcc(i,:),'o-','markersize',4,'Color',color(i,:));
end
xlabel('Saturation level');
ylabel('Median of running CC');
box on; grid on

supertit(ax([1 2]),sprintf('CC vs. # of sources for diff. win length out of full simulation of %d s',...
  Twin),12);


%%
i = 1;
color = jet(length(wlensectry));
figure
ax(1) = subplot(221);
plot(sqrt(nsrc(i,:)),mcc(i,:),'o-','markersize',4,'Color',color(i,:));
xlabel('sqrt(Number of sources inside window)');
ylabel('Overall CC');
box on; grid on

ax(2) = subplot(222);
plot(sqrt(nsrc(i,:)),medrcc(i,:),'o-','markersize',4,'Color',color(i,:));
xlabel('sqrt(Number of sources inside window)');
ylabel('Median of running CC');
box on; grid on

ax(3) = subplot(223);
plot(sqrt(nsat),mcc(i,:),'o-','markersize',4,'Color',color(i,:));
xlabel('sqrt(Saturation level)');
ylabel('Overall CC');
box on; grid on

ax(4) = subplot(224);
plot(sqrt(nsat),medrcc(i,:),'o-','markersize',4,'Color',color(i,:));
xlabel('sqrt(Saturation level)');
ylabel('Median of running CC');
box on; grid on

supertit(ax([1 2]),sprintf('Synthetic record length: %d s / %d s',wlensec, Twin),12);








