% lfetmrdloc.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We need a comparison between the locations of 4-s tremors and lfes that 
% fall inside the same 4-s window of tremor. 
%
% --Besides the direct display in map view of their locations side by side
% and look at them with eyes qualitatively like in 
% 'plt_shortwin_bst181_tmr.m', this script would compute their location
% difference.
% --One way is to plot the dloc between each tremor and all concurrent lfes,
% and regard each as a measurement
% --Another way is plot the dloc between each tremor and the mean (or median)
% of the concurrent lfes. Can you also visualize their direction locations on
% map while aalso being able tell which is which? Maybe choose a certain 
% window size?
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/05/25
% Last modified date:   2024/05/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
format short e   % Set the format to 5-digit floating point
clear
clc
close all

%% for easy testing
set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

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
  % 'LZB  '
  % 'TWKB '
  % 'MGCB '
  'KLNB ']; % determine the trio and order, here the 1st sta is PGC
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
% trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.pgc002.',cutout(1:4)));
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.',num2str(ntol),'.pgc002.',...
  cutout(1:4)));
tlen = trange(:,3)-trange(:,2);
nbst = size(trange,1);
idxbst = 1:length(trange);

dates = unique(trange(:,1));
years = unique(floor(dates/1000));
nets = length(years);
nday = length(dates);

%ppeaks-zcrosses
ppkmzc = [10;12;10;7;10;7;8];

trangeout = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.',num2str(ntol),'.pgcout002.',...
  cutout(1:4)));

%%%load the empirical model param of time offset 14 for each addtional stations
off14mod = load(strcat(rstpath, '/MAPS/timeoff_planefitparam_4thsta_160sps'),'w+');

ftrans = 'interpchao';

%the cut-out boundary of 4-s detections
cutout = 'ellipse';
x0 = 0.2; 
y0 = 0.2;
semia = 1.75;
semib = 1.25;
angrot = 45;
[xcut, ycut] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);

%%%load data, LFE catalog, 25-s-win
% savefile = 'deconv_stats4th_allbstsig.mat';
savefile = 'deconv_stats4th_no23_allbstsig.mat';
load(strcat(rstpath, '/MAPS/',savefile));
% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
load(strcat(rstpath, '/MAPS/',savefile));

%%%param for secondary sources removed
nsrc = allbstsig.nsrc;
imp = allbstsig.impindepall;
trangenew = allbstsig.trangenew;
nsrcn = allbstnoi.nsrc;
impn = allbstnoi.impindepall;
supertstr = 'Secondary sources removed';
fnsuffix = [];

%%%param for further checked at KLNB
nsrc4th = allbstsig.nsrc4th;
imp4th = allbstsig.impindep4thall;
nsrcn4th = allbstnoi.nsrc4th;
impn4th = allbstnoi.impindep4thall;
% supertstr = 'Further checked at KLNB';
% fnsuffix = '4th';

%
sps = 160;
[imploc, ~] = off2space002(imp(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

% keyboard

%% For each tremor, find lfes in the same 4-s time window
ktmr=0;
klfe=0;
for iii = 1: length(idxbst)
  
  [iets,i,j] = indofburst(trangenew,idxbst(iii));
  
  year = years(iets);
  datesets = dates(floor(dates/1000)==year);
  
  date = datesets(i);
  jday = floor(date-year*1000);
  a = jul2dat(year,jday);
  if a(1) == 9
    mo = 'Sep.';
  elseif a(1) == 7
    mo = 'Jul.';
  else
    mo = 'Mar.';
  end
  dy = num2str(a(2));
  yr = num2str(a(3));
  
  %bursts and 4-s tremor detections of the same day
  rangetemp = trangenew(trangenew(:,1)==datesets(i), :);
  hfdayi = hfbnd(hfbnd(:,daycol)==datesets(i), :);  % inside bound of the day

  tmaxi = hfdayi(:, seccol); % starting time of max power rate of half sec inside the ellipse
  tcnti = hfdayi(:, 15);  % the center of detecting win is the 15th col
    
  %%%%Use the start and end of the 4-s detecting window
  tstbuf = rangetemp(j,2); % start and end time of bursts
  tedbuf = rangetemp(j,3);
  tlenbuf = tedbuf-tstbuf;

  %4-s detections inside the burst range
  indtcnti = find(tcnti>=tstbuf & tcnti<=tedbuf); %not every tmr fall into a burst
  ntmrinbst(iii,1) = length(indtcnti);
  tmrii = hfdayi(indtcnti, :);
  tmrtime = tmrii(:, 15) - tstbuf;

  %lfes of the same burst range
  ist = sum(nsrc(1:idxbst(iii)-1))+1;
  ied = ist+nsrc(idxbst(iii))-1;
  impii = imp(ist:ied,:);  %LFEs of the same time win
  implocii = imploc(ist:ied,:);  %LFEs of the same time win
  lfetime = impii(:,1); %lfe timing, at zero-crossing
  timp = impii(:,1)+ppkmzc(1);  %LFE timing, corrected to the positive peak

  %%%for each tremor, stores the concurrent lfes, and compute the dloc
  for jj = 1: length(indtcnti)
    tmrtst = tmrtime(jj)-2;
    tmrted = tmrtime(jj)+2;
    ind = find(lfetime>=tmrtst*sps & lfetime<=tmrted*sps);
    ktmr = ktmr+1;
    nlfeintmr(ktmr,1) = length(ind);
    lfeintmr{ktmr,1} = impii(ind, :);
    lfelocintmr{ktmr,1} = implocii(ind, :);
    mlfeloc = mean(implocii(ind, 1:2), 1); %mean loc of all lfes
    mlfelocspl = mean(implocii(ind, 7:8), 1); %likely non-integer
    dlocm(ktmr,:) = tmrii(jj, 1:2) - mlfeloc;
    dlocsplm(ktmr,:) = tmrii(jj, 7:8) - mlfelocspl; %this contains non-integer
    for jjj = 1: length(ind)
      klfe = klfe+1;
      dloc(klfe,:) = tmrii(jj, 1:2) - implocii(ind(jjj), 1:2);
      dlocspl(klfe,:) = tmrii(jj, 7:8) - implocii(ind(jjj), 7:8);
    end
  end
end  
% keyboard  

%%
nrow = 2;
ncol = 2;
widin = 5.5;  % maximum width allowed is 8.5 inches
htin = 7.3;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);

% pltxran = [0.06 0.98]; pltyran = [0.05 0.98]; % optimal axis location
% pltxsep = 0.01; pltysep = 0.03;
% optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
axpos = [0.07 0.55 0.44 0.4;
         0.60 0.580 0.39 0.34;
         0.07 0.07 0.44 0.4;
         0.60 0.100 0.39 0.34];
for isub = 1:nrow*ncol
  set(f.ax(isub), 'position', axpos(isub,:));
end
  
disttype = 'km';
msize = 3;
cstr={'# events / grid'}; 
marker = 'o';
scale = 'linear';
binmethod = 'grid';
xran=[-4 4]; yran=[-4 4]; %in km
dx=0.025; dy=0.025; %in km, to distinguish dp near origin, size cannot exceed ~0.053 km 
smoothsigma=5;
ncont = 100;

%%%cumulative density for N and N-1 w/i 0.375s 
dplt = dloc;
ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
if strcmp(binmethod,'pixel')
  den1d = density_pixel(dplt(:,1),dplt(:,2));
elseif strcmp(binmethod,'grid')
  den1d = density_matrix(dplt(:,1),dplt(:,2),xran,yran,dx,dy);
end
den1d = den1d(den1d(:,3)>0, :);
den1d = sortrows(den1d,3);

dum = den1d;
dum(dum(:,3)>1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'linew',0.2);  %, 'MarkerEdgeColor', 'k'

dum = sortrows(den1d,3);
dum(dum(:,3)==1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'filled','MarkerEdgeColor','none');
text(ax,0.98,0.05,sprintf('%d pairs',size(dplt,1)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9); %,'FontName','Monospaced'
colormap(ax,'plasma');
c=colorbar(ax,'SouthOutside');
ax.CLim(2) = prctile(dum(:,3),99);
if strcmp(scale,'log10')
  c.Label.String = strcat('log_{10}(',cstr{1},')');
elseif strcmp(scale,'linear')
  c.Label.String = cstr{1};
end
%Principal component analysis
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dplt);
plot(ax,ellx,elly,'-','linew',1.5,'color','k');
% plot(ax,x0+semib*[coeff(1,2),-coeff(1,2)],y0+semib*[coeff(2,2),-coeff(2,2)],...
%   '-','linew',2,'color','b');
% plot(ax,x0+semia*[coeff(1,1),-coeff(1,1)],y0+semia*[coeff(2,1),-coeff(2,1)],...
%   ':','linew',2,'color','b');
% text(ax,0.98,0.16,strcat(num2str(round(anglegeo(2))),'$^{\,\circ}$',{'; '},...
%   num2str(round(anglegeo(1))),'$^{\,\circ}$'),'FontSize',11,...
%   'unit','normalized','interpreter','latex','HorizontalAlignment','right');
text(ax,0.98,0.17,sprintf('%d%c; %d%c',round(anglegeo(2)),char(176),...
  round(anglegeo(1)),char(176)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.98,0.11,strcat({'Asp. ratio: '},sprintf('%.2f',semia/semib)),'Units',...
  'normalized','HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,'a','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
if strcmp(disttype,'spl')
  xlabel(ax,'Diff off12 (samples)');
  ylabel(ax,'Diff off13 (samples)');
  xticks(ax,xran(1):10:xran(2));
  yticks(ax,yran(1):10:yran(2));
elseif strcmp(disttype,'km')
  xlabel(ax,'Location difference E (km)');
  ylabel(ax,'Location difference N (km)');
  xticks(ax,xran(1):1:xran(2));
  yticks(ax,yran(1):1:yran(2));
end
% plot(ax,ax.XLim,ax.YLim,'k--','linew',1);
% pos = ax.Position;
% c.Position = [pos(1), pos(2)-0.12, pos(3), 0.02];
% 0.06 0.54 0.32 0.4;
% 0.4 0.54 0.32 0.4;
c.Position = [0.07 0.55 0.44 0.02];
longticks(ax,2);
% nolabels(ax,3);

if strcmp(disttype,'spl')
  [~,xgrid,ygrid,zgrid] = ...
    zeropadmat2d(den1d,min(den1d(:,1)):1:max(den1d(:,1)),...
    min(den1d(:,2)):1:max(den1d(:,2)));
elseif strcmp(disttype,'km')
  [~,xgrid,ygrid,zgrid] = ...
    zeropadmat2d(den1d,floor(min(den1d(:,1))):dx:ceil(max(den1d(:,1))),...
    floor(min(den1d(:,2))):dy:ceil(max(den1d(:,2))));
end
if ~isempty(smoothsigma)
  zgridgf = imgaussfilt(zgrid,smoothsigma);  %smooth it a bit
else
  zgridgf = zgrid;
end
if strcmp(scale,'log10')
  zgridgf = log10(zgridgf);
end

[conmat,conobj] = contour(ax,xgrid,ygrid,zgridgf,ncont,'-','linew',1);
delete(conobj);
if ~isempty(conmat)
  contable = getContourLineCoordinates(conmat);
  conmat=table2array(contable);
end
%a line cross (0,0) in the projection direction
x = reshape(xran(1):0.001:xran(2), [], 1);
yopt = linefcn(x,tand(angle(2)),0);
plot(ax,x,yopt,'-','linew',2,'color','b');
%a line cross (0,0) in the orthogonal direction
yort = linefcn(x,tand(angle(1)),0);
plot(ax,x,yort,':','linew',2,'color','b');
hold(ax,'off');

%interpolation to obtain the intersection between PCA directions and contours
F = scatteredInterpolant(conmat(:,3),conmat(:,4),conmat(:,1),'linear','none');
normalizer=length(dloc);
ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
[ax,muopt,sigmaopt,mdistprojopt,~,countn1]=plt_dloccrssect(ax,F,x,yopt,anglegeo(2),...
  ['-';'-'],[.5 .5 .5; 0 0 1],xran,'both',normalizer);
[ax,muort,sigmaort,mdistprojort,~,countn1r]=plt_dloccrssect(ax,F,x,yort,anglegeo(1),...
  [':';':'],[.5 .5 .5; 0 0 1],xran,'both',normalizer);
text(ax,0.01,0.65,sprintf('SE \\mu=%.3f;\nSE \\sigma=%.3f;\nmed(|x|)=%.3f',...
  muopt,sigmaopt,mdistprojopt),'Units','normalized',...
  'HorizontalAlignment','left','FontSize',9,'color','b');
text(ax,0.61,0.65,sprintf('NE \\mu=%.3f;\nNE \\sigma=%.3f;\nmed(|x|)=%.3f',...
  muort,sigmaort,mdistprojort),'Units','normalized',...
  'HorizontalAlignment','left','FontSize',9,'color','b');
text(ax,0.02,0.05,'b','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
ylim(ax,[0 1.1e-3]);
% ylabel(ax,'Probability');
% keyboard

%%%cumulative density for each to all later ones in 25-s windows.
dplt = dlocm;
ax=f.ax(3); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
if strcmp(binmethod,'pixel')
  den1dm = density_pixel(dplt(:,1),dplt(:,2));
elseif strcmp(binmethod,'grid')
  den1dm = density_matrix(dplt(:,1),dplt(:,2),xran,yran,dx,dy);
end
den1dm = den1dm(den1dm(:,3)>0, :);
den1dm = sortrows(den1dm,3);

dum = den1dm;
dum(dum(:,3)>1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'linew',0.2);  %, 'MarkerEdgeColor', 'k'

dum = sortrows(den1dm,3);
dum(dum(:,3)==1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'filled','MarkerEdgeColor','none');
text(ax,0.98,0.05,sprintf('%d pairs',size(dplt,1)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
colormap(ax,'plasma');
c=colorbar(ax,'SouthOutside');
ax.CLim(2) = prctile(dum(:,3),99);
if strcmp(scale,'log10')
  c.Label.String = strcat('log_{10}(',cstr{1},')');
elseif strcmp(scale,'linear')
  c.Label.String = cstr{1};
end
%Principal component analysis
[coeff,score,anglem,anglegeom,x0,y0,semia,semib,ellx,elly]=pcaellipse(dplt);
plot(ax,ellx,elly,'-','linew',1.5,'color','k');
% plot(ax,x0+semib*[coeff(1,2),-coeff(1,2)],y0+semib*[coeff(2,2),-coeff(2,2)],...
%   '-','linew',2,'color','r');
% plot(ax,x0+semia*[coeff(1,1),-coeff(1,1)],y0+semia*[coeff(2,1),-coeff(2,1)],...
%   ':','linew',2,'color','r');
% text(ax,0.98,0.16,strcat(num2str(round(anglegeom(2))),'$^{\,\circ}$',{'; '},...
%   num2str(round(anglegeom(1))),'$^{\,\circ}$'),'FontSize',11,...
%   'unit','normalized','interpreter','latex','HorizontalAlignment','right');
text(ax,0.98,0.17,sprintf('%d%c; %d%c',round(anglegeom(2)),char(176),...
  round(anglegeom(1)),char(176)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.98,0.11,strcat({'Asp. ratio: '},sprintf('%.2f',semia/semib)),'Units',...
  'normalized','HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,'c','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
if strcmp(disttype,'spl')
%   xlabel(ax,'Diff off12 (samples)');
%   ylabel(ax,'Diff off13 (samples)');
  xticks(ax,xran(1):10:xran(2));
  yticks(ax,yran(1):10:yran(2));
elseif strcmp(disttype,'km')
  xlabel(ax,'Location difference E (km)');
  ylabel(ax,'Location difference N (km)');
  xticks(ax,xran(1):1:xran(2));
  yticks(ax,yran(1):1:yran(2));
end
% plot(ax,ax.XLim,ax.YLim,'k--','linew',1);
% pos = ax.Position;
% c.Position = [pos(1), pos(2)-0.12, pos(3), 0.02];
% 0.08 0.55 0.4 0.4
c.Position = [0.07 0.07 0.44 0.02];
longticks(ax,2);
% nolabels(ax,3);

if strcmp(disttype,'spl')
  [~,xgrid,ygrid,zgrid] = ...
    zeropadmat2d(den1dm,min(den1dm(:,1)):1:max(den1dm(:,1)),...
    min(den1dm(:,2)):1:max(den1dm(:,2)));
elseif strcmp(disttype,'km')
  [~,xgrid,ygrid,zgrid] = ...
    zeropadmat2d(den1dm,floor(min(den1dm(:,1))):dx:ceil(max(den1dm(:,1))),...
    floor(min(den1dm(:,2))):dy:ceil(max(den1dm(:,2))));
end
if ~isempty(smoothsigma)
  zgridgf = imgaussfilt(zgrid,smoothsigma);  %smooth it a bit
else
  zgridgf = zgrid;
end
if strcmp(scale,'log10')
  zgridgf = log10(zgridgf);
end

[conmatm,conobj] = contour(ax,xgrid,ygrid,zgridgf,ncont,'-','linew',1);
delete(conobj);
if ~isempty(conmatm)
  contable = getContourLineCoordinates(conmatm);
  conmatm=table2array(contable);
end
%a line cross (0,0) in the projection direction
x = reshape(xran(1):0.001:xran(2), [], 1);
yoptm = linefcn(x,tand(anglem(2)),0);
plot(ax,x,yoptm,'-','linew',2,'color','r');
%a line cross (0,0) in the orthogonal direction
yortm = linefcn(x,tand(anglem(1)),0);
plot(ax,x,yortm,':','linew',2,'color','r');
hold(ax,'off');

%interpolation to obtain the intersection between PCA directions and contours
Fm = scatteredInterpolant(conmatm(:,3),conmatm(:,4),conmatm(:,1),'linear','none');
normalizerm=length(dlocm);
ax=f.ax(4); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
[ax,muopt,sigmaopt,mdistprojopt,~,countn2]=plt_dloccrssect(ax,Fm,x,yoptm,...
  anglegeom(2),['-';'-'],[.5 .5 .5; 1 0 0],xran,'both',normalizerm);
[ax,muort,sigmaort,mdistprojort,~,countn2r]=plt_dloccrssect(ax,Fm,x,yortm,...
  anglegeom(1),[':';':'],[.5 .5 .5; 1 0 0],xran,'both',normalizerm);
text(ax,0.01,0.65,sprintf('SE \\mu=%.3f;\nSE \\sigma=%.3f;\nmed(|x|)=%.3f',...
  muopt,sigmaopt,mdistprojopt),'Units','normalized',...
  'HorizontalAlignment','left','FontSize',9,'color','r');
text(ax,0.61,0.65,sprintf('NE \\mu=%.3f;\nNE \\sigma=%.3f;\nmed(|x|)=%.3f',...
  muort,sigmaort,mdistprojort),'Units','normalized',...
  'HorizontalAlignment','left','FontSize',9,'color','r');
text(ax,0.02,0.05,'d','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
ylim(ax,[0 1.1e-3]);
% ylabel(ax,'Probability');

lgd1 = legend(f.ax(2),'SE data','SE gsfit','NE data','NE gsfit',...
  'Location','north','fontsize',7,'NumColumns',2,'orientation','horizontal');
%make background transparent
set(lgd1.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));

lgd2 = legend(f.ax(4),'SE data','SE gsfit','NE data','NE gsfit',...
  'Location','north','fontsize',7,'NumColumns',2,'orientation','horizontal');
%make background transparent
set(lgd2.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));


fname = strcat('lfetmrdloc',fnsuffix,'.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
% print(f.fig,'-dpdf',...
%   strcat('/home/chaosong/Pictures/',fname));


%% analyze the location difference between concurrent tmrs and lfes
% smoothsigma=1;  %smoothing sigma for Gaussian filtering 
% ncont=100;  %num of contour lines
% dx=1; dy=1;
% cstr={'# events / pixel'}; xran=[-40 40]; yran=[-40 40];
% [f,den1d_spl,conmat]=plt_srcdloc(dlocsplm,'spl',3,cstr,...
%   'o','linear','pixel',xran,yran,dx,dy,smoothsigma,ncont);
% fname = strcat('lfetmrdlocsplm',fnsuffix,fnsuffix2,'.pdf');
% print(f.fig,'-dpdf',...
%   strcat('/home/chaosong/Pictures/',fname));
% keyboard
 
% smoothsigma=1;  %smoothing sigma for Gaussian filtering 
% ncont=100;  %num of contour lines
% dx=1; dy=1;
% cstr={'# events / pixel'}; xran=[-40 40]; yran=[-40 40];
% [f,den1d_spl,conmat]=plt_srcdloc(dlocspl,'spl',3,cstr,...
%   'o','linear','pixel',xran,yran,dx,dy,smoothsigma,ncont);
% fname = strcat('lfetmrdlocspl',fnsuffix,'.pdf');
% print(f.fig,'-dpdf',...
%   strcat('/home/chaosong/Pictures/',fname));
% keyboard

% smoothsigma=5;
% ncont=100;  %num of contour lines
% dx=0.025; dy=0.025; % to distinguish dp near origin, size cannot exceed ~0.053 km 
% cstr={'# events / grid'}; xran=[-4 4]; yran=[-4 4];
% [f,den1dm,conmat]=plt_srcdloc(dlocm,'km',3,cstr,...
%   'o','linear','grid',xran,yran,dx,dy,smoothsigma,ncont);
% fname = strcat('lfetmrdlocm',fnsuffix,'.pdf');
% print(f.fig,'-dpdf',...
%   strcat('/home/chaosong/Pictures/',fname));
% keyboard

% smoothsigma=5;
% ncont=100;  %num of contour lines
% dx=0.025; dy=0.025; % to distinguish dp near origin, size cannot exceed ~0.053 km 
% cstr={'# events / grid'}; xran=[-4 4]; yran=[-4 4];
% [f,den1dm,conmat]=plt_srcdloc(dloc,'km',3,cstr,...
%   'o','linear','grid',xran,yran,dx,dy,smoothsigma,ncont);
% fname = strcat('lfetmrdloc',fnsuffix,'.pdf');
% print(f.fig,'-dpdf',...
%   strcat('/home/chaosong/Pictures/',fname));
% keyboard


