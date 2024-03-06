% plt_4stremordensity.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a specific script to plot the density of the 4-s tremor catalog 
% using PGC trio. The high-density region is what we focus and extract the
% the burst windows for deconvolution from. Similar plot has been created 
% for agu2023 in 'compare_lfetremorcats', but i want this script to plot
% all elements, including arrows, etc.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/02/09
% Last modified date:   2024/02/09
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
% 
%% load LFE catalogs
%%%load the LFE catalog, 25-s-win
savefile = 'deconv_stats4th_allbstsig.mat';
allsig = load(strcat(rstpath, '/MAPS/',savefile));
imp = allsig.allbstsig.impindepall;
imp4th = allsig.allbstsig.impindep4thall;
%convert time offset to relative loc
sps = 160;
[imploc, indinput] = off2space002(imp(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[imploc4th, indinput] = off2space002(imp4th(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%%%load the LFE catalog, whole-win
savefile = 'deconv1win_stats4th_allbstsig.mat';
allsig1win = load(strcat(rstpath, '/MAPS/',savefile));
imp1win = allsig1win.allbstsig.impindepall;
imp1win4th = allsig1win.allbstsig.impindep4thall;
%convert time offset to relative loc
[imploc1win, ~] = off2space002(imp1win(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[imploc1win4th, ~] = off2space002(imp1win4th(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%%%load the LFE catalog, 25-s-win, noise
savefile = 'deconv_stats4th_allbstnoi.mat';
allnoi = load(strcat(rstpath, '/MAPS/',savefile));
impn = allnoi.allbstnoi.impindepall;
impn4th = allnoi.allbstnoi.impindep4thall;
%convert time offset to relative loc
[implocn, ~] = off2space002(impn(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[implocn4th, ~] = off2space002(impn4th(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%%%load the LFE catalog, whole-win, noise
savefile = 'deconv1win_stats4th_allbstnoi.mat';
allnoi1win = load(strcat(rstpath, '/MAPS/',savefile));
impn1win = allnoi1win.allbstnoi.impindepall;
impn1win4th = allnoi1win.allbstnoi.impindep4thall;
%convert time offset to relative loc
[implocn1win, ~] = off2space002(impn1win(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[implocn1win4th, ~] = off2space002(impn1win4th(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%% plot the cumulative density map, binning by pixel, ie., each unique detection         
xran=[-6 4];
yran=[-4 4];
binmethod='pixel';
msizehf=6;
scale='log10';
hfplt=hfall;
contourflag=0;


nrow = 1;
ncol = 1;
widin = 5;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.08 0.95]; pltyran = [0.15 0.9];
pltxsep = 0.08; pltysep = 0.05; 
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

ax = f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
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
  colormap(ax,'jet');
%   oldc = colormap(ax,'kelicol');
%   newc = flipud(oldc);
%   colormap(ax,newc);
  % colormap(ax, flipud(colormap(ax,'kelicol')));
  c=colorbar(ax,'SouthOutside');
  pos = ax.Position;
  c.Position = [pos(1), pos(2)-0.03, pos(3), 0.02];
  if strcmp(scale,'log10')
    cstr = strcat({'log_{10}(# detections / '},binmethod,')');
  elseif strcmp(scale,'linear')
    cstr = strcat({'# detections / '},binmethod);  
  end
  c.Label.String = cstr;
%   c.Label.String = strcat({'normalized # tremor detections / '},binmethod,')');
  c.Label.FontSize = 10;
%   caxis(ax,[0 1.7]);

  %plot the high-density ellipse
  plot(ax,xcut,ycut,'k-','LineWidth',2);

  %plot arrows to indicate the semi-long and -short axes
  [rotx, roty] = complex_rot(0,semia,-45);
  xvect = [x0 x0+rotx];
  yvect = [y0 y0+roty];
  drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5);
  [rotx, roty] = complex_rot(0,semib,45);
  xvect = [x0 x0+rotx];
  yvect = [y0 y0+roty];
  drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5);
  
  plot(ax,[x0 x0+2],[y0 y0],'k--','linewidth',1.5);
  text(ax,x0+0.4,y0+0.3,sprintf('\\theta=%d^o',angrot),'FontSize',11,...
      'horizontalalignment','left');
  h=text(ax,0.05,y0+0.2,sprintf('a=%.2fkm',semia),'FontSize',11,'horizontalalignment','left');
%   h=text(ax,0.05,y0+0.2,sprintf('a=%.2fkm',semia),'FontSize',10,'horizontalalignment','left');
  set(h,'Rotation',45);
  h=text(ax,0.2,-0.2,sprintf('b=%.2f',semib),'FontSize',11,'horizontalalignment','right');
  set(h,'Rotation',-45);
  text(ax,0.98,0.95,'4-s Tremor','FontSize',12,'unit','normalized','horizontalalignment','right',...
    'fontweight','bold');
% text(ax, 0.85, 0.93, '2004','FontSize',12,'unit','normalized','horizontalalignment','center',...
%      'EdgeColor','k','Margin',2);
% text(ax,0.85,0.8,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
%   text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
  text(ax,0.5,0.05,strcat(num2str(length(hfplt(:,1))),{' detections'}),'FontSize',10,'unit','normalized',...
      'horizontalalignment','center');
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
xlim(ax,xran);
ylim(ax,yran);
% xticks(ax,xran(1):5:xran(2));
% yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
longticks(ax,2);
hold(ax,'off');

orient(f.fig,'landscape');
fname = '4stremordensity.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

