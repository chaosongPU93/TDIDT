% function synthshift_chao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS should be the modified version of Allan's code 'synthshift',
% integrating my own codes to make synthetic seismograms from LFE
% templates.
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
% First created date:   2022/04/07
% Last modified date:   2023/09/06
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

Twin=0.5*3600+3+2*ceil(greenlen/sps); %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
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
nsat=[0.1 0.4 1 2 4 10 20 40 100];  % times of saturation 
% nsat=[95 100 200];  % times of saturation 
writes=round(nsat*satn) %how many templates to throw in, under different degrees of saturation
rng('default');
%%%Here the noise is 'uniform' in time!
% noistd = 5e-2;
% noistd = 2.e-7;
noistd = 2.0e-4;
synth=noistd*(randn(winlen+greenlen+2*10,nsta)-0.5); %+2*10 a little extra, for jiggering 2nd & 3rd stations.% for ista=1:nsta

nouts=length(writes);
seed=round(writes(1)/5e3); %for random number generator

%%%specify distribution for source location  
% distrloc = 'custompdf'; %using a custom PDF function
distrloc = 'uniform'; %uniformly random in a specified region,

%%%specify which time is uniform in time
% timetype = 'tarvl';
timetype = 'tori';

%%%specify if considering the physical size of each source
% physicalsize = 1;
physicalsize = 0;

%%%specify shape of the source region
srcregion='ellipse';
% srcregion='rectangle';
% srcregion='circle';

%%%specify regime for transformation from time offset to map location 
% ftrans = 'interpArmb';
% ftrans = 'interpArmbreloc';
ftrans = 'interpchao';

%%%specify if forcing a min speration of arrival time for 2 events from the same spot 
% forcesep = 1;
forcesep = 0;

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
    semia = 1.75*(0.6:0.2:2.0);
    semib = 1.25*(0.6:0.2:2.0);
    nreg = length(semia);
    ireg = 8;
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

%% save files
for inwrites=1:length(writes)
  n=writes(inwrites);
  %%%seismograms
  if strcmp(srcregion,'ellipse')
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(xaxis),'-',num2str(yaxis),'.diam',num2str(diam),'.else',num2str(fracelsew,2),...
      'nsat',num2str(nsat(inwrites))],'w');
  elseif strcmp(srcregion,'rectangle')
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(ll(1)),num2str(ll(2)),num2str(ur(1)),num2str(ur(2)),'.diam',num2str(diam),...
      '.else',num2str(fracelsew,2),'nsat',num2str(nsat(inwrites))],'w');
  elseif strcmp(srcregion,'circle')
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(radi),'.diam',num2str(diam),'.else',num2str(fracelsew,2),...
      'nsat',num2str(nsat(inwrites))],'w');  
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
      'nsat',num2str(nsat(inwrites)),'_sources'],'w');
  elseif strcmp(srcregion,'rectangle')
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(ll(1)),num2str(ll,2),num2str(ur(1)),num2str(ur,2),'.diam',num2str(diam),...
      '.else',num2str(fracelsew,2),'nsat',num2str(nsat(inwrites)),'_sources'],'w');
  elseif strcmp(srcregion,'circle')
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(radi),'.diam',num2str(diam),'.else',num2str(fracelsew,2),...
      'nsat',num2str(nsat(inwrites)),'_sources'],'w');  
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
      '_',num2str(xaxis),'-',num2str(yaxis),'.diam',num2str(diam),'_grd'],'w');
  elseif strcmp(srcregion,'rectangle')
    fid = fopen([workpath,'/synthetics/synsrcloc.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(ll(1)),num2str(ll,2),num2str(ur(1)),num2str(ur,2),'.diam',num2str(diam),'_grd'],'w');
  elseif strcmp(srcregion,'circle')
    fid = fopen([workpath,'/synthetics/synsrcloc.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(radi),'.diam',num2str(diam),'_grd'],'w');  
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
      'nsat',num2str(nsat(inwrites)),'_stind'],'w');
  elseif strcmp(srcregion,'rectangle')
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(ll(1)),num2str(ll,2),num2str(ur(1)),num2str(ur,2),'.diam',num2str(diam),...
      '.else',num2str(fracelsew,2),'nsat',num2str(nsat(inwrites)),'_stind'],'w');
  elseif strcmp(srcregion,'circle')
    fid = fopen([workpath,'/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
      '_',num2str(radi),'.diam',num2str(diam),'.else',num2str(fracelsew,2),...
      'nsat',num2str(nsat(inwrites)),'_stind'],'w');  
  end
  tmp = squeeze(greensts{inwrites}{1});
  fprintf(fid,'%d \n', tmp');
    
end

keyboard

%%
%what is the peak-to-peak separation in time for synthetics
figure
for i = 1: length(writes)
  tmp = synths(:,:,i);
  nsta = size(tmp,2);
  pkinds = cell(nsta,1);  % indices of peaks of the waveform
  medpksep = zeros(nsta,1); % median separation of waveform peaks
  ranpksep = zeros(nsta,2); % median separation of waveform peaks
  aa = [];
  for ista = 1: nsta
    [~, pkinds{ista,1}] = findpeaks(tmp(:,ista));
    %   medpksep(ista) = median(diff(pkinds{ista,1}));
    %   ranpksep(ista,1:2) = minmax(diff(pkinds{ista,1})');
  end
  aa = [aa; diff(pkinds{ista,1})];
  subplot(3,3,i)
  histogram(aa,'binw',2);
  xlim([0 80]);
end

figure
hisig=6.3; % this will give a similar spectral shape between template and signal
losig=1.8;
for i = 1: length(writes)
  tmp = synths(:,:,i);
  nsta = size(tmp,2);
  pkinds = cell(nsta,1);  % indices of peaks of the waveform
  medpksep = zeros(nsta,1); % median separation of waveform peaks
  ranpksep = zeros(nsta,2); % median separation of waveform peaks
  aa = [];
  tmp1 = [];
  for ista = 1: nsta
    tmp1(:,ista) = Bandpass(tmp(:,ista), sps, losig, hisig, 2, 2, 'butter');
    [~, pkinds{ista,1}] = findpeaks(tmp1(:,ista));
    %   medpksep(ista) = median(diff(pkinds{ista,1}));
    %   ranpksep(ista,1:2) = minmax(diff(pkinds{ista,1})');
  end
  aa = [aa; diff(pkinds{ista,1})];
  subplot(3,3,i)
  histogram(aa,'binw',2);
  xlim([0 80]);
end

%%%Is the median median. abs amplitude proportional to the square root of number of templates?
%%%The abs amp is similar to envelope
% figure
% hold on; box on
% for i = 1: length(writes)
%   tmp = squeeze(synths(:,:,i));
%   ampmax = median(median(abs(tmp),[],1));
%   ntempsqrt = sqrt(writes(i));
%   ntempsqrt/ampmax
%   scatter(ntempsqrt,ampmax,20,'k','filled');
% end
% xlabel('Square root of number of tmeplates');
% ylabel('Med. of med. of abs amplitude');
    
%% extract and validate the added impulses of template 
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
  
  keyboard
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


%% Best alignment for the testing window
% ccmid = round(size(sigpnstab,1)/2);
% ccwlen = round(size(sigpnstab,1)*0.8);
% loffmax = 5*sps/40;
% ccmin = 0.01;  % depending on the length of trace, cc could be very low
% iup = 1;    % times of upsampling
% mshiftadd=10*sps/40;
% [off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(sigpnstab',ccmid,...
%   ccwlen,mshiftadd,loffmax,ccmin,iup);
% off1i = zeros(nsta,1);
% off1i(2) = round(off12con);
% off1i(3) = round(off13con);
% 
% %align the signals, noises, etc
% sigpnsta = zeros(wlensec*sps, nsta);
% noista = zeros(wlensec*sps, nsta);
% sigsta = zeros(wlensec*sps, nsta);
% sigconvsta = zeros(wlensec*sps, nsta);
% for ista = 1: nsta
%   sigpnsta(:,ista) = sigpnstab(round(buffer/2)+1-off1i(ista): end-round(buffer/2)-off1i(ista), ista); 
%   noista(:,ista) = noistab(round(buffer/2)+1-off1i(ista): end-round(buffer/2)-off1i(ista), ista);
%   sigsta(:,ista) = sigstab(round(buffer/2)+1-off1i(ista): end-round(buffer/2)-off1i(ista), ista); % sta 2
%   sigconvsta(:,ista) = sigconvstab(round(buffer/2)+1-off1i(ista): end-round(buffer/2)-off1i(ista), ista); % sta 2
% end
% 
% srcpairb = srcpairstab{1,1};
% %store the info of synthetic impulses, same format as in paired deconvolution, 
% %9 cols, [ind1 amp1 ind2 amp2 ind3 amp3 off12 off13 off23]  
% srcpair = zeros(size(srcpairb,1), 9);  
% srcpair(:,1) = srcpairb(:,1)-round(buffer/2);  % cut out the buffer 
% srcpair(:,[2 4 6]) = repmat(srcpairb(:,4), 1,3);
% %the indices of synthetic impulses need to be shifted too
% srcpair(:,3) = srcpair(:,1)-srcpairb(:,2)-off1i(2); % note the sign is consistent!
% srcpair(:,5) = srcpair(:,1)-srcpairb(:,3)-off1i(3); 
% srcpair(:,7) = srcpairb(:,2);  % no need to shift offset, because they are 'true' offset
% srcpair(:,8) = srcpairb(:,3);  
% srcpair(:,9) = srcpair(:,8)-srcpair(:,7);  % off23, == off13 - off12 == tarvl2 - tarvl3
% 
% %compute running CC
% % mwlen=sps/2;
% mwlen=sps;
% [ircc,rcc12] = RunningCC(sigpnsta(:,1), sigpnsta(:,2), mwlen);
% [~,rcc13] = RunningCC(sigpnsta(:,1), sigpnsta(:,3), mwlen);
% [~,rcc23] = RunningCC(sigpnsta(:,2), sigpnsta(:,3), mwlen);
% rcc = (rcc12+rcc13+rcc23)/3;  %average
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
% %% independent iterative deconvolution
% for ista = 1:nsta
%   wlet = greenf(:,ista);
%   lwlet = length(wlet);
%   sig = sigpnsta(:,ista);
%   noi = noista(:,ista);
%   
%   %detrend and taper
%   sig = detrend(sig);
%   fractap = 0.05; % if fractap is >=1, n-point von Hann window is returned
%   ptstap = fractap/2*size(sig,1); % if fractap is >=1, n-point von Hann window is returned
%   w = tukeywin(size(sig,1),fractap);
%   % sig = w.* sig;
%   %detrend again for caution
%   sig = detrend(sig);
%   lsig = length(sig);
%   
%   dt = 1/sps;  % sampling interval
%   twlet = zcrosses(ista)*dt;
%   width = 2.5;  % width for Gaussian filter
%   dres_min = 0.5;  % tolerance, percentage change in residual per iteration
%   mfit_min = 5e-1;  % tolerance, norm of the residual/misfit
%   nit_max = 1.5*round(1/tdura*wlensec*nsat(insat));  % max numer of iterations
%   npul_max = round(1/tdura*wlensec*nsat(insat));%a single peak is ~20 samples wide; maybe a little less (at 100 sps). ~0.4s, 1/0.4=2.5
%   fpltit = 0;  % plot flag for each iteration
%   fpltend = 1;  % plot flag for the final iteration
%   fcheck = 0; % plot flag for intermediate computations
%   % rcc = [];  % running CC between diff. stations
%   rcc = ones(length(ircc) ,1); % for testing purpose
%   [sigdecon(:,ista),pred,res,dresit,mfitit,ampit,nit,fighdl] = ...
%     iterdecon(sig,wlet,rcc,noi,dt,twlet,width,dres_min,mfit_min,nit_max,npul_max,fpltit,fpltend,fcheck);
%   ax = fighdl{2}.ax(1);
%   hold(ax,'on');
%   text(ax,0.05,0.85,stas(ista,:),'unit','normalized');
%   hold(ax,'off');
%   nit
% end
% 
% 
% %% Group nearest impulses from different stations into pairs
% % %%%Way 1: find the closest impulse at each station independently
% % %different stations has different number of non-zero impulses
% % for i = 1: nsta
% %   nimp(i) = sum(sigdecon(:,i)>0);
% % end
% % [npair, ista] = max(nimp);
% % imp = [];
% % imp(:,1) = find(sigdecon(:,ista)>0);  % index of impulse
% % imp(:,2) = sigdecon(sigdecon(:,ista)>0, ista);  % amp of impulse
% % imp(:,3) = rcc(imp(:,1)-cclen/2); % rcc value at the impulse
% % imp(:,4) = imp(:,3).*imp(:,2);  % rcc*amp
% % %sort them based on the product of amplitude and rcc?
% % imp = sortrows(imp,4,'descend');
% % 
% % imppair = zeros(npair, 4*3);
% % for ip = 1: npair
% %   imppair(ip,(ista-1)*4+1:ista*4) = imp(ip,1:4);
% % end
% % 
% % spsscale = sps/40;  
% % loff_max = 4*spsscale;  % maximum allowed shift is 4 samples at 40 Hz, it is also true from the synthetic distribution 
% % for i = 1:nsta
% %   if i == ista
% %     continue
% %   end
% %   imp1 = [];
% %   imp1(:,1) = find(sigdecon(:,i)>0);
% %   imp1(:,2) = sigdecon(sigdecon(:,i)>0, i);
% %   imp1(:,3) = rcc(imp1(:,1)-cclen/2);
% %   imp1(:,4) = imp1(:,3).*imp1(:,2);
% %   for ip = 1: npair
% %     [loff,itemp] = min(abs(imp1(:,1)-imp(ip,1))); 
% %     if loff>loff_max
% %       continue
% %     end
% %     imppair(ip,(i-1)*4+1:i*4) = imp1(itemp,1:4); % pair them
% %     imp1(itemp,:) = [];   % mute the used ones
% %   end
% % end
%     
% %%%Way 2: find all pairs that are close enough, choose the one with highest weighted CC
% %different stations has different number of non-zero impulses
% npair = sum(sigdecon(:,1)>0);
% imp1 = [];
% imp1(:,1) = find(sigdecon(:,1)>0);  % index of impulse
% imp1(:,2) = sigdecon(sigdecon(:,1)>0, 1);  % amp of impulse
% imp1(:,3) = rcc(imp1(:,1)-mwlen/2); % rcc value at the impulse
% imp1(:,4) = imp1(:,3).*imp1(:,2);  % rcc*amp
% %sort them based on the product of amplitude and rcc?
% imp1 = sortrows(imp1,4,'descend');
% 
% imppair = zeros(npair, 4*3);
% imppair(:,1:4) = imp1(:,1:4);
% 
% spsscale = sps/40;  
% loff_max = 4*spsscale;  % maximum allowed shift is 4 samples at 40 Hz, it is also true from the synthetic distribution 
% 
% %info of impulses at sta 2
% imp2 = [];
% imp2(:,1) = find(sigdecon(:,2)>0);
% imp2(:,2) = sigdecon(sigdecon(:,2)>0, 2);
% imp2(:,3) = rcc(imp2(:,1)-mwlen/2);
% imp2(:,4) = imp2(:,3).*imp2(:,2);
% 
% %info of impulses at sta 3
% imp3 = [];
% imp3(:,1) = find(sigdecon(:,3)>0);
% imp3(:,2) = sigdecon(sigdecon(:,3)>0, 3);
% imp3(:,3) = rcc(imp3(:,1)-mwlen/2);
% imp3(:,4) = imp3(:,3).*imp3(:,2);
% 
% for ip = 1: npair
%   off12 = imp1(ip,1)-imp2(:,1);
%   off13 = imp1(ip,1)-imp3(:,1);
%   ind2 = find(abs(off12)<=loff_max);  
%   ind3 = find(abs(off13)<=loff_max);
%   off12c = off12(ind2);
%   off13c = off13(ind3);
%   imp2c = imp2(ind2,:);
%   imp3c = imp3(ind3,:);
%   off23c = zeros(length(ind2), length(ind3));
%   for ii = 1: length(ind2)
%     for jj = 1: length(ind3)
%       off23c(ii,jj) = off13c(jj)-off12c(ii);
%     end
%   end
%   ind23 = find(abs(off23c)<=loff_max);
%   if isempty(ind23)
%     continue
%   else
%     [sub2, sub3] = ind2sub(size(off23c),ind23); % convert linear indices to matrix subscripts
%   end
%   sumwtcoef = imp2c(sub2,4)+imp3c(sub3,4); % obtain sum of the weighted coefs at these pairs
%   [msumwtcoef,ind] = max(sumwtcoef);   % choose the pair has the highest weighted coef sum
%   bsub = [sub2(ind), sub3(ind)];
%   imppair(ip,(2-1)*4+1:2*4) = imp2(ind2(sub2(ind)),1:4);
%   imp2(ind2(sub2(ind)),:) = [];   % mute the used ones
%   imppair(ip,(3-1)*4+1:3*4) = imp3(ind3(sub3(ind)),1:4);
%   imp3(ind3(sub3(ind)),:) = [];   % mute the used ones
% 
% end
%     
% indpair = find(sum(imppair==0,2)==0); %pairs that have corresponding peaks at all stations
% imppairf = imppair(indpair,:);  % found pairs 
% impindep = imppairf(:, [1 2 5 6 9 10]); % info of impulse index and amp 
% %adding the time offset 12, 13 and 23, indicating location, off23 == off13-off12 == tarvl2-tarvl3
% impindep(:,7:9) = [impindep(:,1)-impindep(:,3) impindep(:,1)-impindep(:,5) ...
%                    impindep(:,3)-impindep(:,5)];  % again, be consistent in sign!
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
% 
% %%%scatter of offsets, accounting for prealignment offset, == true offset 
% offxran = [-loff_max loff_max]; 
% offyran = [-loff_max loff_max];
% cran = [0 lsig];
% impindepst = sortrows(impindep,1);
% impindepst(:,7:8) = impindepst(:,7:8)+repmat([off1i(2) off1i(3)],size(impindepst,1),1); %account for prealignment
% [f] = plt_decon_imp_scatter(impindepst,offxran,offyran,cran,sps,50,'mean');
% title(gca,sprintf('Independent, grouped, using data of %d Hz',sps));
% 
% [f] = plt_decon_imp_scatter(srcpair,offxran,offyran,cran,sps,50,'mean');
% title(gca,sprintf('Synthetic sources, using data of %d Hz',sps));
% 
% % %%% scatter of rela locations, and account for prealignment offset 
% % xran = [-3 3]; 
% % yran = [-3 3];
% % cran = [0 lsig];
% % [f] = plt_decon_imp_scatter_space(impindepst,xran,yran,cran,offxran,offyran,sps,50,ftrans,...
% %     'mean');
% % title(gca,sprintf('Independent, grouped, using data of %d Hz',sps));
% % 
% % [f] = plt_decon_imp_scatter_space(srcpair,xran,yran,cran,offxran,offyran,sps,50,ftrans,...
% %     'mean');
% % title(gca,sprintf('Synthetic sources, using data of %d Hz',sps));









