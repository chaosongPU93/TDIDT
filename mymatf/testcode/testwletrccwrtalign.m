% testwletrccwrtalign.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is to test how CC/RCC changes wrt. to different alignments between
% templates at stations. Similarly, you can also ask how CC/RCC changes
% wrt. to different alignments between signals at stations, however, that
% would be burst dependent, as data is changing. In contrast, if using
% templates for similar analysis, result would be invariant.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/01/04
% Last modified date:   2023/01/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e   % Set the format to 5-digit floating point
clear
clc
close all

defval('normflag',0); %whether to normalize templates
defval('rccmwsec',0.5); %moving win len in sec for computing RCC

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
set(0,'DefaultFigureVisible','on');

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
ntol = 3;
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.',...
  num2str(ntol),'.pgc002.',cutout(1:4)));
tlen = trange(:,3)-trange(:,2);
nbst = size(trange,1);

%%%load the empirical model param of time offset 14 for each addtional stations
off14mod = load(strcat(rstpath, '/MAPS/timeoff_planefitparam_4thsta_160sps'),'w+');

%% prepare templates (Green's functions), from 'lfetemp002_160sps.m'
sps = 160;
templensec = 60;

ccstack = [];
for ista = 1: nsta
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
        num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao_catnew');
    ccstack(:,ista) = load(fname);
end
STA = detrend(ccstack);

ccstackort = [];
for ista = 1: nsta
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
        num2str(templensec), 's_', 'BBCCS_', 'ort_Nof_Non_Chao_catnew');
    ccstackort(:,ista) = load(fname);
end
STAort = detrend(ccstackort);

%flag of normalization
% normflg = 0;

% %plot the raw templates, not filtered, not best aligned
% figure
% subplot(211)
% hold on
% for ista = 1: nsta
%   plot(STA(:,ista));
% end
% subplot(212)
% hold on
% for ista = 1: nsta
%   plot(STAort(:,ista));
% end

%%%The below aligns the templates by x-correlation
ist = templensec*sps*4/10;  %not using whole window in case any station has very long-period energy
ied = templensec*sps*6/10;
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
    STAtmport(:,ista)=detrend(STAort(is(ista):ie(ista),ista)); 
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
    STAtmport(mshiftadd+1:end-(mshiftadd+1),ista)=detrend(STAtmport(mshiftadd+1-imax(ista):end-(mshiftadd+1)-imax(ista),ista));
end
%normalization
if normflag 
  for ista=1:nsta
      STAtmp(:,ista)=STAtmp(:,ista)/spread(ista); % now templates are 'aligned' indeoendently by x-corr wrt. sta 1
      STAtmport(:,ista)=STAtmport(:,ista)/spread(ista);
  end
end
% figure
% subplot(211)
% hold on
% for ista = 1: nsta
%   plot(STAtmp(:,ista));
% end
% subplot(212)
% hold on
% for ista = 1: nsta
%   plot(STAtmport(:,ista));
% end
%%%The above aligns the templates by x-correlation

%%%detrend, taper and bandpass templates
tmpwlet = STAtmp; % no bandpass
tmpwletf = STAtmp;  % bandpassed version
fractap = sps/size(tmpwlet,1);
tmpwletort = STAtmport; % no bandpass
tmpwletfort = STAtmport;  % bandpassed version
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
  
  %same process for orthogonal
  tmpwletort(:,ista) = detrend(tmpwletort(:,ista));
  tmpwletort(:,ista) = w.* tmpwletort(:,ista);
  tmpwletort(:,ista) = detrend(tmpwletort(:,ista));
  tmpwletfort(:,ista) = Bandpass(tmpwletort(:,ista), sps, lowlet, hiwlet, 2, 2, 'butter');
  tmpwletfort(:,ista) = detrend(tmpwletfort(:,ista));
end

%%%constrained CC, so that only 2 offsets are independent
ccmid = round(size(tmpwletf,1)/2);
ccwlen = 10*sps;
loffmax = 5*sps/40;
ccmin = 0.01;  % depending on the length of trace, cc could be very low
iup = 1;    % times of upsampling
[off12con,off13con,~] = constrained_cc_interp(tmpwletf(:,1:3)',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
offwlet1i(1) = 0;
offwlet1i(2) = round(off12con);
offwlet1i(3) = round(off13con);

% offmax = sps/10;
offmax = 28;
% offmax = 24;
off12 = -offmax: 1: offmax;
off13 = -offmax: 1: offmax;

[X1,X2] = meshgrid(off12,off13);
tmp = [X1(:) X2(:)];
npts = length(tmp);

offwlet = zeros(npts,3);
offwlet(:,2) = tmp(:,1)+offwlet1i(2);
offwlet(:,3) = tmp(:,2)+offwlet1i(3);

rccmwlen=rccmwsec*sps;

for ipt = 1: npts
  ipt
  %%%automatically find the rough zero-crossing time, whose abs. value is closest to 0, whether + or -
  [~,imin] = min(tmpwletf(:,1));
  [~,imax] = max(tmpwletf(:,1));
  [~,zcsta1] = min(abs(tmpwletf(imin:imax,1)));
  zcsta1 = zcsta1+imin-1;
  greenlen = pow2(9)*sps/40;
  greenf = zeros(greenlen,nsta);  % bandpassed version
  ppeaks = zeros(nsta,1); % positive peaks
  npeaks = zeros(nsta,1); % negative peaks
  for ista = 1: nsta
    %cut according to the zero-crossing and the time shift from the constrained CC
    greenf(:,ista) = tmpwletf(zcsta1+8*sps-greenlen+1-offwlet(ipt,ista): zcsta1+8*sps-offwlet(ipt,ista), ista);
    %detrend again for caution
    greenf(:,ista) = detrend(greenf(:,ista));
    if normflag
      %normalize by max amp
      greenf(:,ista) = greenf(:,ista)/max(abs(green(:,ista)));    % normalize
    end
    
    %re-find the zero-crossing as the template length has changed
    [~,imin] = min(greenf(:,ista));
    [~,imax] = max(greenf(:,ista));
    [~,zcrosses(ista)] = min(abs(greenf(imin:imax,ista)));
    zcrosses(ista) = zcrosses(ista)+imin-1;
    ppeaks(ista) = imax;
    npeaks(ista) = imin;
  end
  
  %compute running CC between 3 stations
  [ircc,rcc12,rcc13,rcc23] = RunningCC3sta(greenf,rccmwlen);
  rccpair = [rcc12 rcc13 rcc23];
  meanmedrcc123(ipt) = mean(median(rccpair,1));
  meanmeanrcc123(ipt) = mean(mean(rccpair,1));
  
  cc12 = xcorr(greenf(:,1), greenf(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
  cc13 = xcorr(greenf(:,1), greenf(:,3),0,'normalized');
  cc23 = xcorr(greenf(:,2), greenf(:,3),0,'normalized');
  ccpair = [cc12 cc13 cc23];
  meancc123(ipt) = mean(ccpair);
  
%   [~,ind] = min(ccpair);
%   cc(ipt) = mean(ccpair(:,setdiff(1:3,ind)), 2);
%   medrcc(ipt) = median(mean(rccpair(:,setdiff(1:3,ind)), 2));
%   maxrcc(ipt) = max(mean(rccpair(:,setdiff(1:3,ind)), 2));

  meancc12(ipt) = mean(ccpair(:,[1 2]));
  meanmedrcc12(ipt) = mean(median(rccpair(:,[1 2]),1));
  meanmeanrcc12(ipt) = mean(mean(rccpair(:,[1 2]),1));

  meancc13(ipt) = mean(ccpair(:,[1 3]));
  meanmedrcc13(ipt) = mean(median(rccpair(:,[1 3]),1));
  meanmeanrcc13(ipt) = mean(mean(rccpair(:,[1 3]),1));

  meancc23(ipt) = mean(ccpair(:,[2 3]));
  meanmedrcc23(ipt) = mean(median(rccpair(:,[2 3]),1));
  meanmeanrcc23(ipt) = mean(mean(rccpair(:,[2 3]),1));
  
end

meanmedrcc123 = reshape(meanmedrcc123,length(off13),length(off12));
meanmeanrcc123 = reshape(meanmeanrcc123,length(off13),length(off12));
meanmedrcc12 = reshape(meanmedrcc12,length(off13),length(off12));
meanmeanrcc12 = reshape(meanmeanrcc12,length(off13),length(off12));
meanmedrcc13 = reshape(meanmedrcc13,length(off13),length(off12));
meanmeanrcc13 = reshape(meanmeanrcc13,length(off13),length(off12));
meanmedrcc23 = reshape(meanmedrcc23,length(off13),length(off12));
meanmeanrcc23 = reshape(meanmeanrcc23,length(off13),length(off12));
meancc123 = reshape(meancc123,length(off13),length(off12));
meancc12 = reshape(meancc12,length(off13),length(off12));
meancc13 = reshape(meancc13,length(off13),length(off12));
meancc23 = reshape(meancc23,length(off13),length(off12));

%the cut-out boundary of 4-s detections
ftrans = 'interpchao';
x0 = 0.2; 
y0 = 0.2;
semia = 1.75;
semib = 1.25;
angrot = 45;
[xcut, ycut] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);
%what's the shape of it in sample space?
offcut = space2off002([xcut-x0, ycut-y0],sps,ftrans,0);

%%%range that encloses about 90 percent of LFE detections aftering shifting back
%%%to the origin
load('90thprcrangeoflfes.mat');

%%
%%%plot to show how does the RCC/CC change wrt. the diff alignment between sta pairs 12 and 13
%%%note that the result is gonna be burst dependent, as data is changing
%%%in contrast, if using templates for similar analysis, result would be invariant
widin = 6.5;  % maximum width allowed is 8.5 inches
htin = 6.5;   % maximum height allowed is 11 inches
nrow = 2;
ncol = 2;
f1=initfig(widin,htin,nrow,ncol);

figxran = [0.1 0.98]; figyran = [0.1 0.97];
figxsep = 0.05; figysep = 0.08;
optaxpos(f1,nrow,ncol,figxran,figyran,figxsep,figysep);

matplt = meancc12;
ax=f1.ax(1);
hold(ax,'on'); ax.Box = 'on'; 
imagesc(ax,off12/sps,off13/sps,matplt);
plot(ax,offcut(:,1)/sps,offcut(:,2)/sps,'k-','linew',1.5);
loffm = 6*4;  %at 160 Hz
detoffm = [-loffm loffm; -loffm -loffm; loffm -loffm; loffm loffm; -loffm loffm];
plot(ax,detoffm(:,1)/sps,detoffm(:,2)/sps,'-','linew',1.5,'color',[.5 .5 .5]);
%plot the range of 90 percent of detected LFEs
plot(ax,conmat10th(:,3),conmat10th(:,4),'r-','LineWidth',1.5);
[C,h] = contour(ax,X1/sps,X2/sps,matplt,'k-','ShowText','on','LineWidth',0.5);
clabel(C,h,'FontSize',8);
% scatter(ax,lfe(:,7)/sps,lfe(:,8)/sps,10,[.7 .7 .7],'filled','MarkerEdgeColor','k');
% plot(ax,minmax(off12)/sps,minmax(off13)/sps,'--','Color',[.3 .3 .3],'linew',1);
ind = find(matplt==max(matplt,[],'all'));
[sub13, sub12] = ind2sub(size(matplt),ind);
scatter(ax,off12(sub12)/sps,off13(sub13)/sps,30,'k^','linew',1.5);
% colormap(ax,'jet');
colormap(ax, flipud(colormap(ax,'kelicol')));
% caxis(ax,[-1 1]);
% c=colorbar(ax);
% c.Label.String = 'mean overall CC of pairs 12 and 13';
title(ax,'Average CC of pairs 12 and 13','FontSize',10,'FontWeight','bold');  %,'FontWeight','normal'
text(ax,0.02,0.95,'a','FontSize',11,'unit','normalized','EdgeColor','k',...
  'Margin',1,'color','k','backgroundcolor','w');
axis(ax,'equal','tight');
xlim(ax,minmax(off12)/sps);
ylim(ax,minmax(off13)/sps);
xticks(ax,-0.1: 0.1: 0.1);
yticks(ax,-0.1: 0.1: 0.1);
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
% xlabel(ax,sprintf('\\Delta{t}_{12} (s)'));
ylabel(ax,sprintf('\\Delta{t}_{13} (s)'));
longticks(ax,1);

matplt = meancc13;
ax=f1.ax(2); hold(ax,'on'); ax.Box = 'on'; 
imagesc(ax,off12/sps,off13/sps,matplt);
% plot(ax,offcut(:,1)/sps,offcut(:,2)/sps,'b-','linew',1.5);
[C,h] = contour(ax,X1/sps,X2/sps,matplt,'k-','ShowText','on','LineWidth',0.5);
clabel(C,h,'FontSize',8);
% plot(ax,minmax(off12)/sps,minmax(off13)/sps,'--','Color',[.3 .3 .3],'linew',1);
ind = find(matplt==max(matplt,[],'all'));
[sub13, sub12] = ind2sub(size(matplt),ind);
scatter(ax,off12(sub12)/sps,off13(sub13)/sps,30,'k^','linew',1.5);
% colormap(ax,'jet');
colormap(ax, flipud(colormap(ax,'kelicol')));
% caxis(ax,[-1 1]);
% c=colorbar(ax);
% c.Label.String = 'median of mean RCC of 3 pairs';
title(ax,'Average CC of pairs 12 and 23','FontSize',10,'FontWeight','normal');
text(ax,0.02,0.95,'b','FontSize',11,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w'); %,'backgroundcolor','w'
axis(ax,'equal','tight');
xlim(ax,minmax(off12)/sps);
ylim(ax,minmax(off13)/sps);
xticks(ax,-0.1: 0.1: 0.1);
yticks(ax,-0.1: 0.1: 0.1);
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
% xlabel(ax,sprintf('PGC-SSIB offset (samples at %d Hz)',sps),'FontSize',11);
% ylabel(ax,sprintf('PGC-SILB offset (samples at %d Hz)',sps),'FontSize',11);
% xlabel(ax,sprintf('\\Delta{t}_{12} (s)'));
% ylabel(ax,sprintf('\\Delta{t}_{13} (s)'));
longticks(ax,1);

matplt = meancc23;
ax=f1.ax(3); hold(ax,'on'); ax.Box = 'on'; 
imagesc(ax,off12/sps,off13/sps,matplt);
% plot(ax,offcut(:,1)/sps,offcut(:,2)/sps,'b-','linew',1.5);
[C,h] = contour(ax,X1/sps,X2/sps,matplt,'k-','ShowText','on','LineWidth',0.5);
clabel(C,h,'FontSize',8);
% scatter(ax,lfe(:,7)/sps,lfe(:,8)/sps,10,[.7 .7 .7],'filled','MarkerEdgeColor','k');
% plot(ax,minmax(off12)/sps,minmax(off13)/sps,'--','Color',[.3 .3 .3],'linew',1);
ind = find(matplt==max(matplt,[],'all'));
[sub13, sub12] = ind2sub(size(matplt),ind);
scatter(ax,off12(sub12)/sps,off13(sub13)/sps,30,'k^','linew',1.5);
% colormap(ax,'jet');
colormap(ax, flipud(colormap(ax,'kelicol')));
% caxis(ax,[-1 1]);
% c=colorbar(ax);
% c.Label.String = 'median of mean RCC of pairs 12 and 13';
title(ax,'Average CC of pairs 13 and 23','FontSize',10,'FontWeight','normal');
text(ax,0.02,0.95,'c','FontSize',11,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
axis(ax,'equal','tight');
xlim(ax,minmax(off12)/sps);
ylim(ax,minmax(off13)/sps);
xticks(ax,-0.1: 0.1: 0.1);
yticks(ax,-0.1: 0.1: 0.1);
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
xlabel(ax,sprintf('\\Delta{t}_{12} (s)'));
ylabel(ax,sprintf('\\Delta{t}_{13} (s)'));
longticks(ax,1);

matplt = meancc123;
ax=f1.ax(4);
hold(ax,'on'); ax.Box = 'on'; 
imagesc(ax,off12/sps,off13/sps,matplt);
% plot(ax,offcut(:,1)/sps,offcut(:,2)/sps,'b-','linew',1.5);
[C,h] = contour(ax,X1/sps,X2/sps,matplt,'k-','ShowText','on','LineWidth',0.5);
clabel(C,h,'FontSize',8);
% conmat = contour(ax,X1,X2,matplt,-0.15:0.01:0.3,'k-','ShowText','on'); %
% plot(ax,minmax(off12)/sps,minmax(off13)/sps,'--','Color',[.3 .3 .3],'linew',1);
ind = find(matplt==max(matplt,[],'all'));
[sub13, sub12] = ind2sub(size(matplt),ind);
scatter(ax,off12(sub12)/sps,off13(sub13)/sps,30,'k^','linew',1.5);
% colormap(ax,'jet');
colormap(ax, flipud(colormap(ax,'kelicol')));
% caxis(ax,[-1 1]);
% c=colorbar(ax);
% c.Label.String = 'mean overall CC of 3 pairs';
title(ax,'Average CC of all 3 pairs','FontSize',10,'FontWeight','normal');
text(ax,0.02,0.95,'d','FontSize',11,'unit','normalized','EdgeColor','k',...
  'Margin',1,'color','k','backgroundcolor','w');
axis(ax,'equal','tight');
xlim(ax,minmax(off12)/sps);
ylim(ax,minmax(off13)/sps);
xticks(ax,-0.1: 0.1: 0.1);
yticks(ax,-0.1: 0.1: 0.1);
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
xlabel(ax,sprintf('\\Delta{t}_{12} (s)'));
% ylabel(ax,sprintf('\\Delta{t}_{13} (s)'));
longticks(ax,1);

% orient(f1.fig,'landscape');
fname = 'wletccwrtoffset.pdf';
print(f1.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));


keyboard

%% choose a set of offset12, 13 with equal CC/RCC
ccval = 0.7;
indcol = find(conmat(1,:)==ccval);
nind = length(indcol);
for ipt = 1: nind
  erroff = conmat(:,indcol(ipt)+1: indcol(ipt)+conmat(2,indcol(ipt)))';
end
%fit an ellipse to it
ellfit = fit_ellipse(erroff(:,1),erroff(:,2));

%%%add the indication of the error ellipse
mmmax = 10;
nnmax = 10;
erroff1 = zeros((mmmax*2+1)^2,2);
kk = 0;
for mm = -mmmax:1:mmmax
  for nn = -nnmax:1:nnmax
    kk = kk+1;
    erroff1(kk,:) = [mm nn];
  end
end
ftrans = 'interpchao';
[errloc1, ~] = off2space002(erroff1,sps,ftrans,0);

F1 = scatteredInterpolant(erroff1(:,1),erroff1(:,2),errloc1(:,1),'linear');
F2 = scatteredInterpolant(erroff1(:,1),erroff1(:,2),errloc1(:,2),'linear');

errloc = [];
errloc(:,1) = F1(erroff(:,1),erroff(:,2));
errloc(:,2) = F2(erroff(:,1),erroff(:,2));


widin = 8; htin = 4;
nrow = 1; ncol = 2;
[f2] = initfig(widin,htin,nrow,ncol);
xran = [0.1 0.95]; yran = [0.1 0.95];
xsep = 0.08; ysep = 0.05;
optaxpos(f2,nrow,ncol,xran,yran,xsep,ysep);

ax=f2.ax(1);
hold(ax,'on'); ax.Box = 'on'; 
plot(ax,[-mmmax mmmax],[-nnmax nnmax],'--','Color',[.3 .3 .3],'linew',1);
plot(ax,erroff(:,1),erroff(:,2),'k','linew',1.5);
x0 = ellfit.X0_in; 
y0 = ellfit.Y0_in;
semia = ellfit.a; 
semib = ellfit.b;
rotang = -rad2deg(ellfit.phi);
rotcnt = [x0, y0];
[xe, ye] = ellipse_chao(x0,y0,semia,semib,0.01,rotang,rotcnt);
plot(ax,xe,ye,'r-');
[xa,ya] = complex_rot(x0,y0+semib,rotang,rotcnt);
plot(ax,[x0 xa],[y0 ya],'b--','linew',1);
[xb,yb] = complex_rot(x0+semia,y0,rotang,rotcnt);
plot(ax,[x0 xb],[y0 yb],'b--','linew',1);
text(ax,0.95,0.1,strcat(sprintf('%.1f; %.1f; %.1f',semia,semib,rotang),' {\circ}'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',10);
xlabel(ax,sprintf('PGC-SSIB offset (samples at %d Hz)',sps),'FontSize',11);
ylabel(ax,sprintf('PGC-SILB offset (samples at %d Hz)',sps),'FontSize',11);
axis(ax,'equal','tight');
axis(ax,[-mmmax mmmax -nnmax nnmax]);
title(ax,sprintf('Contour of CC/RCC=%.1f',ccval));

ax=f2.ax(2);
hold(ax,'on'); ax.Box = 'on'; 
plot(ax,errloc(:,1),errloc(:,2),'k','linew',1.5);
plot(ax,[-1 1],[-1 1],'--','Color',[.3 .3 .3],'linew',1);
xlabel(ax,'E (km)','FontSize',11);
ylabel(ax,'N (km)','FontSize',11);
axis(ax,'equal','tight');
axis(ax,[-1 1 -1 1]);

%% what is the shape like for sqrt(off12^2+off13^2+off23^2)==const
const = 8;

widin = 8; htin = 4;
nrow = 1; ncol = 2;
[f3] = initfig(widin,htin,nrow,ncol);
xran = [0.1 0.95]; yran = [0.1 0.95];
xsep = 0.08; ysep = 0.05;
optaxpos(f3,nrow,ncol,xran,yran,xsep,ysep);

ax=f3.ax(1);
hold(ax,'on'); ax.Box = 'on'; 
plot(ax,[-mmmax mmmax],[-nnmax nnmax],'--','Color',[.3 .3 .3],'linew',1);
syms x y
func1 = @(x,y) x.^2+y.^2-const^2;
fimplicit(ax,func1,'color',[.3 .3 .3],'linew',1);
func2 = @(x,y) x.^2+y.^2-x.*y-const^2/2;
fp=fimplicit(ax,func2,'k','linew',1.5);
%fit an ellipse to it
erroff = [reshape(fp.XData,[],1) reshape(fp.YData,[],1)];
ellfit = fit_ellipse(erroff(:,1),erroff(:,2));
x0 = ellfit.X0_in; 
y0 = ellfit.Y0_in;
semia = ellfit.a; 
semib = ellfit.b;
rotang = -rad2deg(ellfit.phi);
rotcnt = [x0, y0];
[xe, ye] = ellipse_chao(x0,y0,semia,semib,0.01,rotang,rotcnt);
plot(ax,xe,ye,'r-');
[xa,ya] = complex_rot(x0,y0+semib,rotang,rotcnt);
plot(ax,[x0 xa],[y0 ya],'b--','linew',1);
[xb,yb] = complex_rot(x0+semia,y0,rotang,rotcnt);
plot(ax,[x0 xb],[y0 yb],'b--','linew',1);
text(ax,0.95,0.1,strcat(sprintf('%.1f; %.1f; %.1f',semia,semib,rotang),' {\circ}'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',10);
xlabel(ax,sprintf('PGC-SSIB offset (samples at %d Hz)',sps),'FontSize',11);
ylabel(ax,sprintf('PGC-SILB offset (samples at %d Hz)',sps),'FontSize',11);
axis(ax,'equal','tight');
axis(ax,[-mmmax mmmax -nnmax nnmax]);
xticks(ax,-mmmax: 2: mmmax);
yticks(ax,-nnmax: 2: nnmax);
title(ax,strcat('$\sqrt(x^2+y^2-(y-x)^2)=$',sprintf(' %d',const)),'Interpreter','latex');

ax=f3.ax(2);
hold(ax,'on'); ax.Box = 'on';
errloc = [];
errloc(:,1) = F1(erroff(:,1),erroff(:,2));
errloc(:,2) = F2(erroff(:,1),erroff(:,2));
plot(ax,errloc(:,1),errloc(:,2),'k','linew',1.5);
plot(ax,[-1 1],[-1 1],'--','Color',[.3 .3 .3],'linew',1);
xlabel(ax,'E (km)','FontSize',11);
ylabel(ax,'N (km)','FontSize',11);
axis(ax,'equal','tight');
axis(ax,[-1 1 -1 1]);





