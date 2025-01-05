% plt_tmrlfedensity.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a specific script to plot the density of the 4-s tremor catalog
% using PGC trio, AND short-win 3-sta and 4-sta data catalogs, that is used
% in the paper for the time-saturation of tremor.
% The tremor density plot is basically the same as "plt_4stremordensit.m",
% and lfe density is basicall the same as "plt_lfecat_comparison.m". This
% script combines them into one plot, with minor edits
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/11/05
% Last modified date:   2024/11/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
format short e   % Set the format to 5-digit floating point
clear
clc
close all

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
  'LZB  '
  'TWKB '
  'MGCB '
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
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.pgc002.',cutout(1:4)));
tlen = trange(:,3)-trange(:,2);
nbst = size(trange,1);

dates = unique(trange(:,1));
years = unique(floor(dates/1000));
nets = length(years);

%%%load the empirical model param of time offset 14 for each addtional stations
off14mod = load(strcat(rstpath, '/MAPS/timeoff_planefitparam_4thsta_160sps'),'w+');

%% load other catalogs
%%% load the LFE catalog of Michael Bostock, inside and outside the rectangle in 'locinterp002_4s.m'
%obtain the location of fam 002, lon0 and lat0
ftrans = 'interpchao';
loc0 = off2space002([0 0],sps*iup,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%the cut-out boundary of 4-s detections
cutout = 'ellipse';
x0 = 0.2;
y0 = 0.2;
semia = 1.75;
semib = 1.25;
angrot = 45;
[xcut, ycut] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);

% %same rectangle in 'locinterp002_4s.m'
% EW = [-7 3];
% NS = [-3 4];
% wid = range(EW);
% hgt = range(NS);
% x0 = mean(EW);
% y0 = mean(NS);
% [x, y] = rectangle_chao(x0,y0,wid,hgt,0.01);

%%%detecting range of LFEs
loffm = 6;  %at 40 Hz
detranoff = [-loffm*ones(2*loffm+1,1) (loffm:-1:-loffm)';
  (-loffm:1:loffm)' -loffm*ones(2*loffm+1,1);
  loffm*ones(2*loffm+1,1) (-loffm:1:loffm)';
  (loffm:-1:-loffm)' loffm*ones(2*loffm+1,1)];
detranloc = off2space002(detranoff,sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%%%range that encloses about 90 percent of LFE detections aftering shifting back
%%%to the origin
load('90thprcrangeoflfes.mat');

%% load LFE catalogs
%%%load the LFE catalog, 25-s-win
% savefile = 'deconv_stats4th_allbstsig.mat';
savefile = 'deconv_stats4th_no23_allbstsig.mat';
allsig = load(strcat(rstpath, '/MAPS/',savefile));
imp = allsig.allbstsig.impindepall;
imp4th = allsig.allbstsig.impindep4thall;
off1i = allsig.allbstsig.off1i;
off1iwk = allsig.allbstsig.off1iwk;
nsrc = allsig.allbstsig.nsrc;
nsrc4th = allsig.allbstsig.nsrc4th;
irccrank = allsig.allbstsig.irccrank;
stats = allsig.allbstsig.fitstats;
%convert time offset to relative loc
sps=160;
[imploc, ~] = off2space002(imp(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[imploc4th, ~] = off2space002(imp4th(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
%shift the alignment back to the origin
imps=shiftsrctoorigin(imp,nsrc,nbst,off1i,off1iwk,irccrank);
imp4ths=shiftsrctoorigin(imp4th,nsrc4th,nbst,off1i,off1iwk,irccrank);
[implocs, ~] = off2space002(imps(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[imploc4ths, ~] = off2space002(imp4ths(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13


%% A combined plot
nrow = 1;
ncol = 3;
widin = 8.4;  % maximum width allowed is 8.5 inches
htin = 3.5;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);

% pltxran = [0.1 0.95]; pltyran = [0.15 0.9];
% pltxsep = 0.08; pltysep = 0.05;
% axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

axpos = [0.075 0.08 0.3125 0.8;
  0.42 0.08 0.25 0.8;
  0.70 0.08 0.25 0.8];
for i = 1:nrow*ncol
  set(f.ax(i), 'position', axpos(i,:));
end

% plot the cumulative density map, binning by pixel, ie., each unique detection
xran=[-6 4];
yran=[-4 4];
binmethod='pixel';
msizehf=4;
scale='linear';
hfplt=hfall;
contourflag=0;

ax = f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
ax.FontSize = 9;
xlim(ax,xran);
ylim(ax,yran);
% xticks(ax,xran(1):5:xran(2));
% yticks(ax,yran(1):5:yran(2));
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
if ~isempty(hfplt)
  %create a density matrix to store the number of detections in each small grid
  if isequal(binmethod,'grid')
    %     dxhf = 0.2;
    %     dyhf = 0.2;
    [density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
      xran,yran,dxhf,dyhf);
    marker = 'o';
    %     msizehf = 500*dxhf*dyhf;
    %     msizelf = 2.5*msizehf;
    %bin based upon pixel
  elseif isequal(binmethod,'pixel')
    density1d = density_pixel(hfplt(:,1),hfplt(:,2));
    marker = 'o';
    %     msizehf = 6;
    %     msizelf = 2.5*msizehf;
  end
  dum = density1d(density1d(:,3)>0, :);
  normalizer = max(dum(:,3));
  dum(dum(:,3)>1, :) = [];
  if strcmp(scale,'log10')
    dum(:,3) = log10(dum(:,3));
  end
  scatter(ax,dum(:,1),dum(:,2),msizehf,dum(:,3),marker,'linew',0.3);  %, 'MarkerEdgeColor', 'w')
  %   scatter(ax,dum(:,1),dum(:,2), msizehf, dum(:,3)/normalizer,marker,'linew',0.2);  %, 'MarkerEdgeColor', 'w')
  dum = sortrows(density1d(density1d(:,3)>0, :), 3);
  dum(dum(:,3)==1, :) = [];
  if strcmp(scale,'log10')
    dum(:,3) = log10(dum(:,3));
  end
  scatter(ax,dum(:,1),dum(:,2),msizehf,dum(:,3),marker,'filled','MarkerEdgeColor','none');  %,
  %   scatter(ax,dum(:,1),dum(:,2), msizehf, dum(:,3)/normalizer,marker,'filled','MarkerEdgeColor','none');  %,
  % imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
  % colormap(ax,flipud(colormap(ax,'kelicol')));
  colormap(ax,'plasma');
  c1=colorbar(ax,'SouthOutside');
  pos = ax.Position;
  c1.Position = [pos(1), pos(2)+0.05, pos(3), 0.02];
  c1pos = c1.Position;
  if strcmp(scale,'log10')
    cstr = strcat({'log_{10}(# events / '},binmethod,')');
  elseif strcmp(scale,'linear')
    cstr = strcat({'# events / '},binmethod);
  end
  c1.Label.String = cstr;
  % c1.Label.FontSize = 10;
  %   caxis(ax,[0 1.7]);
  
  %plot the high-density ellipse
  plot(ax,xcut,ycut,'k-','LineWidth',2);
  
  %plot the loc of fam 002
%   scatter(ax,0,0,12,'k^');
  %   text(ax,-0.5,-0.5,sprintf('002'),'FontSize',9,'horizontalalignment','left','fontweight','bold');
  
  %   %plot the adopted range for detecting LFEs
  %   plot(ax,detranloc(:,1)+x0,detranloc(:,2)+y0,'r-','LineWidth',1);
  
  %plot the range of 90 percent of detected LFEs
  plot(ax,conmatxy10th(:,3)+x0,conmatxy10th(:,4)+y0,'r-','LineWidth',1.5);
  
  %plot arrows to indicate the semi-long and -short axes
  [rotx, roty] = complex_rot(0,semia,-45);
  xvect = [x0 x0+rotx];
  yvect = [y0 y0+roty];
  drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1);
  [rotx, roty] = complex_rot(0,semib,45);
  xvect = [x0 x0+rotx];
  yvect = [y0 y0+roty];
  drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1);
  
  plot(ax,[x0 x0+2],[y0 y0],'k--','linewidth',1);
  
  text(ax,0.6,0.45,sprintf('\\theta'),'FontSize',11,...
    'horizontalalignment','left','fontweight','bold');
  h=text(ax,0.3,0.9,'a','FontSize',10,'horizontalalignment','left','fontweight','bold');
  set(h,'Rotation',45);
  h=text(ax,-0.5,0.5,'b','FontSize',10,'horizontalalignment','left','fontweight','bold');
  set(h,'Rotation',-45);
  text(ax,0.02,0.95,'Tremor','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'fontweight','bold');
  text(ax,0.98,0.94,sprintf('a=%.2f km',semia),'FontSize',10,'horizontalalignment',...
    'right','unit','normalized');
  text(ax,0.98,0.86,sprintf('b=%.2f km',semib),'FontSize',10,'horizontalalignment',...
    'right','unit','normalized');
  text(ax,0.98,0.78,sprintf('\\theta=%d^o',angrot),'FontSize',10,...
    'horizontalalignment','right','unit','normalized');
  % text(ax, 0.85, 0.93, '2004','FontSize',12,'unit','normalized','horizontalalignment','center',...
  %      'EdgeColor','k','Margin',2);
  % text(ax,0.85,0.8,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
  text(ax,0.02,0.06,'a','FontSize',10,'unit','normalized','EdgeColor','k',...
    'Margin',1,'backgroundcolor','w');
  text(ax,0.98,0.05,strcat(num2str(length(hfplt(:,1))),{' events'}),'FontSize',9,'unit','normalized',...
    'horizontalalignment','right');
  %%%Add contour lines if needed
  if contourflag
    [xyzgridpad,xgrid,ygrid,zgrid,ind2] = zeropadmat2d(density1d,xran(1):1:xran(2),yran(1):1:yran(2));
    zgridgf = imgaussfilt(zgrid, 1);  %smooth it a bit
    perc = 50:10:90;
    conplt = prctile(dum(:,3),perc);
    if strcmp(scale,'log10')
      zgridgf = log10(zgridgf);
    end
    conmat = contour(ax,xgrid,ygrid,zgridgf,conplt,'-','color',[.3 .3 .3]); %,'ShowText','on'
    %     conplt = prctile(dum(:,3)/normalizer,perc);
    %     conmat = contour(ax,xgrid,ygrid,zgridgf/normalizer,conplt,'-','color',[.3 .3 .3]); %,'ShowText','on'
  end
end
xlabel(ax,'E (km)','fontsize',10);
ylabel(ax,'N (km)','fontsize',10);
longticks(ax,2);
hold(ax,'off');

%%% obtain the contour and density matrix of the lfe catalogs
impplt=imp(:,7:8);
den = density_pixel(impplt(:,1),impplt(:,2));
den = den(den(:,3)>0, :);
den = sortrows(den,3);
tmploc = off2space002(den(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
den(:,1:2) = tmploc(:,1:2);

impplt=imp4th(:,7:8);
den4th = density_pixel(impplt(:,1),impplt(:,2));
den4th = den4th(den4th(:,3)>0, :);
den4th = sortrows(den4th,3);
tmploc = off2space002(den4th(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
den4th(:,1:2) = tmploc(:,1:2);

xran = [-4 4];  %in km
yran = [-4 4];
% dx = 0.025;
% dy = 0.025;
binmethod = 'pixel';
marker = 'o';
msize = 4;
disttype = 'km';
contourflag=1;
smoothsigma=[];
scale = 'linear';

%%%3-station, Short-win, data
denplt = den;
ax = f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
ax.FontSize = 9;
xlim(ax,xran);
ylim(ax,yran);
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
%scatter the density
dum = denplt;
dum(dum(:,3)>1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'linew',0.2);  %, 'MarkerEdgeColor', 'w')
dum = sortrows(denplt,3);
dum(dum(:,3)==1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'filled','MarkerEdgeColor','none');
text(ax,0.98,0.05,sprintf('%d events',sum(denplt(:,3))),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
colormap(ax,'plasma');
c2=colorbar(ax,'SouthOutside');
if strcmp(scale,'log10')
  cstr = strcat({'log_{10}(# events / '},binmethod,')');
elseif strcmp(scale,'linear')
  cstr = strcat({'# events / '},binmethod);
end
c2.Label.String = cstr;

xlabel(ax,'E (km)','FontSize',10);
% ylabel(ax,'N (km)','FontSize',10);
longticks(ax,2);
pos = ax.Position;
c2.Position = [pos(1), c1pos(2), pos(3), c1pos(4)];
c2pos = c2.Position;
% hold(ax,'off');
% caxis(ax,[0 1.7]);
plot(ax,xcut,ycut,'k-','linew',2);
text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',10);
text(ax,0.02,0.95,'LFE','FontSize',11,'unit','normalized','horizontalalignment','left',...
  'fontweight','bold');
text(ax,0.02,0.06,'b','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');

% keyboard
% %%
%%%One migration example, let's use 174
ibst = 174;
impi = findimpofburst(imp,nsrc,ibst);
imploci = findimpofburst(imploc,nsrc,ibst);
statsi = stats{ibst};

% ttol = 35;
% ntol = 3;
% tranbstbuf = load(strcat(rstpath, '/MAPS/tdec.bstranbuf',...
%   num2str(ttol),'s.pgc002.',cutout(1:4)));
tranbstbuf=[2005255,35790,36158];
tlen = tranbstbuf(3)-tranbstbuf(2);
tstbuf = tranbstbuf(2);

torispl=tarvl2tori(impi,sps,ftrans);

[iets,i,j] = indofburst(trange,ibst);

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


ax = f.ax(3); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
plot(ax,xcut,ycut,'k-','linew',2);
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
ax.FontSize = 9;
xlim(ax,xran);
ylim(ax,yran);
cran = [0 tlen];
msize = 30;
if ~isempty(impi)
  wt = mean(impi(:,[2 4 6]),2);
  wtmax = prctile(wt,95); %use percentile in case
  refscl = wt./wtmax;
  refscl(refscl>=1) = 1;  %force the larger amp to be plotted as the same size in case of saturation
  scatter(ax,imploci(:,1),imploci(:,2),msize*refscl,torispl/sps,'filled','o',...
    'MarkerEdgeColor',[.5 .5 .5]);
end
text(ax,0.99,0.05,sprintf('%d events',size(imploci,1)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
% colormap(ax,flipud(colormap(ax,'kelicol')));
colormap(ax,'viridis');
c3=colorbar(ax,'SouthOutside');
caxis(ax,cran);
c3.Label.String = sprintf('Origin (s) since %d s on %s %s %s',tstbuf,dy,mo,yr);
c3.Label.FontSize = 9;
pos = ax.Position;
c3.Position = [pos(1), c2pos(2), pos(3), c2pos(4)];
scatter(ax,xran(1)+0.1*range(xran),yran(2)-0.05*range(yran),msize,'w','filled',...
  'MarkerEdgeColor',[.5 .5 .5],'linew',1);
text(ax,0.02,0.9,strcat({'Amplitude '},'$\geq$',{' 95th prctile'}) ,'Units','normalized',...
  'HorizontalAlignment','left','FontSize',8,'interpreter','latex');
% xticks(ax,xran(1):1:xran(2));
% yticks(ax,yran(1):1:yran(2));
if ~isempty(statsi)
  angrmse = statsi.angrmse;
  [rotx, roty] = complex_rot(0,1,-angrmse);
  xarrow = [0.5-rotx 0.5+rotx];
  yarrow = [-2.5-roty -2.5+roty];
  % [xarrown,yarrown] = ds2nfu(ax,xarrow,yarrow);
  % p=annotation('arrow',xarrown,yarrown,'color','k','linestyle','-','linewidth',1.5);
  p=annotation('arrow','color','k','linestyle','-','linewidth',1.5);
  p.Parent = ax;
  p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
  % text(ax,0.6,0.2,strcat(num2str(angrmse),'$^{\,\circ}$'),'FontSize',11,...
  %     'unit','normalized','interpreter','latex','HorizontalAlignment','left');
  text(ax,0.65,0.2,sprintf('%d%c',angrmse,char(176)),...
    'Units','normalized','HorizontalAlignment','left','FontSize',10);
end
% text(ax,0.99,0.95,'Secondary removed','HorizontalAlignment','right','Units',...
%   'normalized','FontSize',9);
text(ax,0.99,0.95,'3-station','HorizontalAlignment','right','Units',...
  'normalized','FontSize',10);
text(ax,0.02,0.06,'c','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
xlabel(ax,'E (km)');
% ylabel(ax,'N (km)');
longticks(ax,2);


% %%%4-station, Short-win, data
% denplt = den4th;
% ax = f.ax(3); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
% ax.GridLineStyle = '--';
% ax.XAxisLocation = 'top';
% ax.FontSize = 9;
% xlim(ax,xran);
% ylim(ax,yran);
% plot(ax,[-100 100],[0 0],'k--');
% plot(ax,[0 0],[-100 100],'k--');
% %scatter the density
% dum = denplt;
% dum(dum(:,3)>1, :) = [];
% if strcmp(scale,'log10')
%   dum(:,3) = log10(dum(:,3));
% end
% scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'linew',0.2);  %, 'MarkerEdgeColor', 'w')
% dum = sortrows(denplt,3);
% dum(dum(:,3)==1, :) = [];
% if strcmp(scale,'log10')
%   dum(:,3) = log10(dum(:,3));
% end
% scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'filled','MarkerEdgeColor','none');
% text(ax,0.98,0.05,sprintf('%d events',sum(denplt(:,3))),'Units','normalized',...
%   'HorizontalAlignment','right','FontSize',9);
% colormap(ax,'plasma');
% c3=colorbar(ax,'SouthOutside');
% if strcmp(scale,'log10')
%   cstr = strcat({'log_{10}(# events / '},binmethod,')');
% elseif strcmp(scale,'linear')
%   cstr = strcat({'# events / '},binmethod);
% end
% c3.Label.String = cstr;

% xlabel(ax,'E (km)','FontSize',10);
% ylabel(ax,'N (km)','FontSize',10);
% longticks(ax,2);
% pos = ax.Position;
% % c3.Position = [pos(1), pos(2)-0.03, pos(3), 0.02];
% c3.Position = [pos(1), c1pos(2), pos(3), c1pos(4)];
% % hold(ax,'off');
% caxis(ax,c2.Limits);
% plot(ax,xcut,ycut,'k-','linew',2);
% text(ax,0.98,0.95,sprintf('4-station'),'Units','normalized',...
%   'HorizontalAlignment','right','FontSize',10);
% text(ax,0.02,0.95,'LFE','FontSize',11,'unit','normalized','horizontalalignment','left',...
%   'fontweight','bold');
% text(ax,0.02,0.06,'c','FontSize',10,'unit','normalized','EdgeColor','k',...
%   'Margin',1,'backgroundcolor','w');
% hold(ax,'off');
keyboard
% orient(f.fig,'landscape');
fname = 'tmrlfedensity.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

