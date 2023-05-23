% decon_allansyn.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the script to read and deconvolve the synthetic seismograms from
% Allan's experiments on using region size and temporal saturation level.
% The deconvolution is implemented by the same algorithm in 
% 'deconv_4s_exp_4thsta_fn.m' by regarding the full window as a whole. I 
% might compare results from using different LFE templates (Allan's and
% mine). The goal is to get the location, and distance between event 
% pairs (N and N-m, or each source to all others within some time 
% separation).
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/02/02
% Last modified date:   2023/02/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e   % Set the format to 5-digit floating point
clear
clc
% close all

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

%plot the raw templates, not filtered, not best aligned
figure
ltemp = size(STA,1);
subplot(111)
hold on
for ista = 1: nsta
  plot((1:ltemp)/sps,STA(:,ista));
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
% plt_templates(green,greenf,stas,greenort,greenfort,lowlet,hiwlet,sps);

%just the filtered templates
plt_templates_bp(greenf,stas,lowlet,hiwlet,sps);

%% load data, ie., Allan's synthetics
% sps = 100;  %note I have been using 160 hz
% 
% regionflag = 'small';
% % regionflag = 'large';
% 
% if strcmp(regionflag,'small') %if source region is a small ellipse
%   fname = 'STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0.25nsat';
%   nsat = [0.1;0.4;1;2;4;10;40];
%   nnsat = length(nsat);
% elseif strcmp(regionflag,'large') %if source region is a bigger ellipse
%   fname = 'STAS.UN.100sps.ell_2.75-1.25.diam0.3.else0nsat';
%   nsat = [0.1;0.4;1;2;4;10;100];
%   nnsat = length(nsat);
% end

fname = '/synthetics/STAS.UN.160sps.ell_2-1.25.diam0.3.else0nsat';
nsat = [0.1; 0.4; 1; 2; 4; 10; 100];
nnsat = length(nsat);

STAopt = cell(nnsat,1);
for i = 1: nnsat
  STAopt{i} = load(strcat(workpath,fname,num2str(nsat(i))));
end

%%
isat = 5;

optseg = STAopt{isat};

%%%example plot of data
lsig = size(optseg,1);
figure
subplot(311)
plot((1:lsig)/sps,optseg(:,1),'r');
subplot(312)
plot((1:lsig)/sps,optseg(:,2),'b');
subplot(313)
plot((1:lsig)/sps,optseg(:,3),'k');

%filter data
hisig=6.3; % this will give a similar spectral shape between template and signal
losig=1.8;
for ista = 1:nsta
  optseg(:,ista) = Bandpass(optseg(:,ista), sps, losig, hisig, 2, 2, 'butter');
end

%here we don't further align them
optdat = optseg; 

off1i = zeros(1,nsta);

%%%2022/06/06, do NOT taper whatsoever!!
sigsta = zeros(size(optdat,1), nsta);
for ista = 1:nsta
  tmp = optdat(:,ista); %best aligned, filtered
  %detrend and taper only the data, NOT the noise
  tmp = detrend(tmp);
  sigsta(:,ista) = tmp;
end
%compute running CC between 3 stations
rccmwlen=rccmwsec*sps;
overshoot = rccmwlen/2;
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

figure
lsig = size(sigsta,1);
plot((1:lsig)/sps,sigsta(:,1),'r'); hold on
plot((1:lsig)/sps,sigsta(:,2),'b');
plot((1:lsig)/sps,sigsta(:,3),'k');
ax = gca;
axsym(ax);
plot((1:lsig)/sps,rcc*ax.YLim(2),'o','color',[.6 .6 .6],'markersize',2);
text(0.95,0.1,sprintf('Saturation: %.1f',nsat(isat)),'Units','normalized','HorizontalAlignment','right');
xlim([0 20]);


%%
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
  nit_max = round(1.5*1/tdura*(tlen));  % max numer of iterations
  nimp_max = round(1/tdura*(tlen));%a single peak is ~20 samples wide; maybe a little less (at 100 sps). ~0.4s, 1/0.4=2.5
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

%% Group nearest impulses from different stations into pairs, using moving searching range
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
  sps,50,'mean','tori');
scatter(ax1,off1i(2),off1i(3),20,'ks','filled','MarkerEdgeColor','k');
title(ax1,'Grouped Total');

%%
grid = load('tmp.grd');

%%
ind = union(impindepst(:,7:8),grid(:,1:2),'rows','stable');

%%
%note the 'tsep' obtained from the deconvolved positive peaks should be identical to that if
%obtained from the deconvolved impulses themselves, which represent the arrival indices of the
%zero-crossing
ista=1;
ftrans = 'interpchao';

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
  sps,50,'mean','tori');
scatter(ax1,off1i(2),off1i(3),20,'ks','filled','MarkerEdgeColor','k');
title(ax1,'Secondary sources removed');

% keyboard

%% what is the distance for consecutive sourcces, but in terms of arrival time?
%note the 'tsep' obtained from the deconvolved positive peaks should be identical to that if
%obtained from the deconvolved impulses themselves, which represent the arrival indices of the
%zero-crossing
ista=1;
ftrans = 'interpchao';

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
end

%%%Projected distance along specific directions, eg., propagation and its orthogonal, in terms of arrival time
[locxyproj,dlocxyproj,stats] = srcprojdistNtoNm(tarvlsplst,implocst,m,sps);
if ~isempty(dlocxyproj)
  nsep = 1;
  ttype = 'tarvl';
%   dtarvlprojall = [dtarvlprojall; dtarvl{nsep}];
%   distarvlprojall = [distarvlprojall; dlocxyproj{nsep}];
%   locxyprojall = [locxyprojall; locxyproj];
%   projangrm(k,1) = stats.angrmse;
%   projangsl(k,1) = stats.angslope;
%   projpear(k,1) = stats.pearwt;
  
  [f] = plt_srcprojdist(tarvlsplst,implocst,dtarvl{nsep},eucdist{nsep},...
    locxyproj,dlocxyproj{nsep},stats,sps,ttype);
  %         close(f.fig);
end
%       orient(f.fig,'landscape');
%       if noiseflag
%         print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',ttype,'nn1distnoi.pdf'));
%       else
%         print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',ttype,'nn1dist.pdf'));
%       end
% keyboard

%%%what are the corresponding RCC at each source
rccpairsrc = [];
rccpairsrc(:,1) = rccpair(round(mean(impindepst(:,[1 3]),2)),1);
rccpairsrc(:,2) = rccpair(round(mean(impindepst(:,[1 5]),2)),2);
rccpairsrc(:,3) = rccpair(round(mean(impindepst(:,[3 5]),2)),3);
% rccpairsrcall = [rccpairsrcall; rccpairsrc];

%use the concatenated rcc at the average arrival time of each source
rccsrc = [];
rccsrc(:,1) = rcc(round(mean(impindepst(:,[1 3 5]),2)));
% rccsrcall = [rccsrcall; rccsrc];


% keyboard



