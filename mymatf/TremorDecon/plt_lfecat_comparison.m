% plt_lfecat_comparison.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code aims to plot the hit count (density) of the LFE detected in
% different types of catalogs: whole-win 3-station data, whole-win 4-sta data, 
% short-win 3-sta data, short-win 4-sta data, whole-win 3-station noise, 
% whole-win 4-sta noise, FOR PAPER. 
% --Note that tremor has been plotted separately, so need to plot thot, but
% technically this code would use the same way of plotting.
% --For noise, we should plot the whole-win detection, instead of the short-win
% one, as using short-win might introduce some artificial change in the center
% of the short windows, thus may lead to an artificial migrating pattern.
% --This comparison is currently in sample space. the advantage is that you can
% see an elongation in the density (or contours) in data detection, but rather
% equvidimensional for noise. This indcates not only that noise is indeed
% behaving like noise cause the detections uniformly distributed, but also that
% your detection algorithm is not introducing any systematic bias.
% --See also 'compare_lfetremorcats.m', 'plt_4stremordensity.m'.  
%
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/03/22
% Last modified date:   2024/03/22
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
ntol = 3;
% trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.pgc002.',cutout(1:4)));
% trange = trange(1:end-1,:);
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.',num2str(ntol),'.pgc002.',...
  cutout(1:4)));
nbst = size(trange,1);

dates = unique(trange(:,1));
years = unique(floor(dates/1000));
nets = length(years);

trangeout = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.',num2str(ntol),'.pgcout002.',...
  cutout(1:4)));

%%%load the empirical model param of time offset 14 for each addtional stations
off14mod = load(strcat(rstpath, '/MAPS/timeoff_planefitparam_4thsta_160sps'),'w+');

%% load other catalogs
%%% load the LFE catalog of Michael Bostock, inside and outside the rectangle in 'locinterp002_4s.m'
%obtain the location of fam 002, lon0 and lat0
sps = 160;
ftrans = 'interpchao';
loc0 = off2space002([0 0],sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%the cut-out boundary of 4-s detections
cutout = 'ellipse';
x0 = 0.2; 
y0 = 0.2;
semia = 1.75;
semib = 1.25;
angrot = 45;
[xcut, ycut] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);

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
%convert time offset to relative loc
[imploc, ~] = off2space002(imp(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[imploc4th, ~] = off2space002(imp4th(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
%shift the alignment back to the origin
imps=shiftsrctoorigin(imp,nsrc,nbst,off1i,off1iwk,irccrank);
imp4ths=shiftsrctoorigin(imp4th,nsrc4th,nbst,off1i,off1iwk,irccrank);
[implocs, ~] = off2space002(imps(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[imploc4ths, ~] = off2space002(imp4ths(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%%%load the LFE catalog, whole-win
% savefile = 'deconv1win_stats4th_allbstsig.mat';
savefile = 'deconv1win_stats4th_no23_allbstsig.mat';
allsig1w = load(strcat(rstpath, '/MAPS/',savefile));
imp1w = allsig1w.allbstsig.impindepall;
imp1w4th = allsig1w.allbstsig.impindep4thall;
off1i1w = allsig1w.allbstsig.off1i;
nsrc1w = allsig1w.allbstsig.nsrc;
nsrc1w4th = allsig1w.allbstsig.nsrc4th;
%convert time offset to relative loc
[imploc1w, ~] = off2space002(imp1w(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[imploc1w4th, ~] = off2space002(imp1w4th(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
%shift the alignment back to the origin
imp1ws=shiftsrctoorigin(imp1w,nsrc1w,nbst,off1i1w);
imp1w4ths=shiftsrctoorigin(imp1w4th,nsrc1w4th,nbst,off1i1w);
[imploc1ws, ~] = off2space002(imp1ws(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[imploc1w4ths, ~] = off2space002(imp1w4ths(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%%%load the LFE catalog, 25-s-win, noise
% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
allnoi = load(strcat(rstpath, '/MAPS/',savefile));
impn = allnoi.allbstnoi.impindepall;
impn4th = allnoi.allbstnoi.impindep4thall;
off1in = allnoi.allbstnoi.off1i;
off1iwkn = allnoi.allbstnoi.off1iwk;
nsrcn = allnoi.allbstnoi.nsrc;
nsrcn4th = allnoi.allbstnoi.nsrc4th;
irccrankn = allnoi.allbstnoi.irccrank;
%convert time offset to relative loc
[implocn, ~] = off2space002(impn(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[implocn4th, ~] = off2space002(impn4th(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
%shift the alignment back to the origin
impns=shiftsrctoorigin(impn,nsrcn,nbst,off1in,off1iwkn,irccrankn);
impn4ths=shiftsrctoorigin(impn4th,nsrcn4th,nbst,off1in,off1iwkn,irccrankn);
[implocns, ~] = off2space002(impns(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[implocn4ths, ~] = off2space002(impn4ths(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%%%load the LFE catalog, whole-win, noise
% savefile = 'deconv1win_stats4th_allbstnoi.mat';
savefile = 'deconv1win_stats4th_no23_allbstnoi.mat';
allnoi1w = load(strcat(rstpath, '/MAPS/',savefile));
impn1w = allnoi1w.allbstnoi.impindepall;
impn1w4th = allnoi1w.allbstnoi.impindep4thall;
off1in1w = allnoi1w.allbstnoi.off1i;
nsrcn1w = allnoi1w.allbstnoi.nsrc;
nsrcn1w4th = allnoi1w.allbstnoi.nsrc4th;
%convert time offset to relative loc
[implocn1w, ~] = off2space002(impn1w(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[implocn1w4th, ~] = off2space002(impn1w4th(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
%shift the alignment back to the origin
impn1ws=shiftsrctoorigin(impn1w,nsrcn1w,nbst,off1in1w);
impn1w4ths=shiftsrctoorigin(impn1w4th,nsrcn1w4th,nbst,off1in1w);
[implocn1ws, ~] = off2space002(impn1ws(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[implocn1w4ths, ~] = off2space002(impn1w4ths(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
% keyboard

%% SAMPLE space, SHIFT to origin, cumulative density & cross-sections, binning by pixel, only short-win
nrow = 2;
ncol = 3;
widin = 2.5*ncol;  % maximum width allowed is 8.5 inches
htin = 7;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.06 0.98]; pltyran = [0.05 0.98]; % optimal axis location
pltxsep = 0.01; pltysep = 0.03;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
 
xran = [-40 40];  %in samples
yran = [-40 40];
dx = 1; dy = 1; 
binmethod = 'pixel';
marker = 'o';
msize = 6;
disttype = 'spl';
contourflag=1;
smoothsigma=1;
scale = 'log10';
ncont = 500;
which2plt = 'data'; %'gsfit'

%short-win data, 3-sta
ax = f.ax(1); 
impplt=imps(:,7:8);
[ax,den1ds,conmats,angle,anglegeo,c] = ...
  plt_cumulative_density_axis(ax,impplt,xran,yran,...
  [],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
% %Principal component analysis
% [coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(impplt);
hold(ax,'on');  
dxln = dx/sps/10; %in sec
loff_max = 6/40;  %in sec
x = reshape(-loff_max:dxln:loff_max, [], 1);
% x = reshape(-0.1: dxln:0.1, [], 1); %need to truncate if plotting 'gsfit'
%a line cross (0,0) in the projection direction
yopt = linefcn(x,tand(angle(2)),0);
plot(ax,x,yopt,'-','linew',2,'color','b');
%a line cross (0,0) in the orthogonal direction
yort = linefcn(x,tand(angle(1)),0);
plot(ax,x,yort,':','linew',2,'color','b');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.055, pos(3), 0.02];
caxis(ax,[0 1.9]);
% cran=c.Limits;
text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, data'),'Units','normalized',...
  'FontSize',9);
% text(ax,0.98,0.12,strcat(num2str(round(anglegeo(2))),'$^{\,\circ}$',{'; '},...
%   num2str(round(anglegeo(1))),'$^{\,\circ}$'),'FontSize',11,...
%   'unit','normalized','interpreter','latex','HorizontalAlignment','right');
text(ax,0.98,0.12,sprintf('%d%c; %d%c',round(anglegeo(2)),char(176),...
  round(anglegeo(1)),char(176)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',10);
text(ax,0.02,0.06,'a','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1);
hold(ax,'off');  

%short-win noise, 3-sta
ax = f.ax(2); 
impplt=impns(:,7:8);
[ax,den1dns,conmatns,anglen,anglegeon,c,~,conobj10th] = ...
  plt_cumulative_density_axis(ax,impplt,xran,yran,...
  [],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
hold(ax,'on');
delete(conobj10th);
plot(ax,x,yopt,'-','linew',2,'color','r');
plot(ax,x,yort,':','linew',2,'color','r');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.055, pos(3), 0.02];
caxis(ax,[0 1.9]);
xlabel(ax,'');
ylabel(ax,'');
nolabels(ax,3);
text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, noise'),'Units','normalized',...
  'FontSize',9);
% text(ax,0.98,0.12,strcat(num2str(round(anglegeo(2))),'$^{\,\circ}$',{'; '},...
%   num2str(round(anglegeo(1))),'$^{\,\circ}$'),'FontSize',11,...
%   'unit','normalized','interpreter','latex','HorizontalAlignment','right');
text(ax,0.98,0.12,sprintf('%d%c; %d%c',round(anglegeo(2)),char(176),...
  round(anglegeo(1)),char(176)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',10);
text(ax,0.02,0.06,'b','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1);
hold(ax,'off');  

%short-win cross section profile, data and noise, 3-sta 
ax = f.ax(3);
c1pos = c.Position; 
f.ax(3).Position = [c1pos(1)+c1pos(3)+0.05 c1pos(2)+0.033 0.27 0.32];
normalizer=length(imps);
hold(ax,'on');  
%contour lines are log10 scale, now revert to linear scale
F = scatteredInterpolant(conmats(:,3),conmats(:,4),10.^(conmats(:,1)),'linear','none');
Fn = scatteredInterpolant(conmatns(:,3),conmatns(:,4),10.^(conmatns(:,1)),'linear','none');
[ax,avgopt,stdopt,muopt,sigmaopt,mdistprojopt]=plt_loccrssect(ax,F,x,yopt,anglegeo(2),...
  'b','-',xran/sps,which2plt,normalizer);
[ax,avgort,stdort,muort,sigmaort,mdistprojort]=plt_loccrssect(ax,F,x,yort,anglegeo(1),...
  'b',':',xran/sps,which2plt,normalizer);
%use the SAME projection direction from data to noise
[ax,avgoptn,stdoptn,muoptn,sigmaoptn,mdistprojoptn]=plt_loccrssect(ax,Fn,x,yopt,anglegeo(2),...
  'r','-',xran/sps,which2plt,normalizer);
[ax,avgortn,stdortn,muortn,sigmaortn,mdistprojortn]=plt_loccrssect(ax,Fn,x,yort,anglegeo(1),...
  'r',':',xran/sps,which2plt,normalizer);
text(ax,0.01,0.7,sprintf('SE \\mu=%.3f;\nSE \\sigma=%.3f',muopt,sigmaopt),'Units',...
  'normalized','HorizontalAlignment','left','FontSize',8,'color','b');
text(ax,0.01,0.55,sprintf('NE \\mu=%.3f;\nNE \\sigma=%.3f',muort,sigmaort),'Units',...
  'normalized','HorizontalAlignment','left','FontSize',8,'color','b');
text(ax,0.65,0.7,sprintf('SE \\mu=%.3f;\nSE \\sigma=%.3f',muoptn,sigmaoptn),'Units',...
  'normalized','HorizontalAlignment','left','FontSize',8,'color','r');
text(ax,0.65,0.55,sprintf('NE \\mu=%.3f;\nNE \\sigma=%.3f',muortn,sigmaortn),'Units',...
  'normalized','HorizontalAlignment','left','FontSize',8,'color','r');
% text(ax,0.98,0.25,sprintf('3-station'),'Units','normalized',...
%   'HorizontalAlignment','right','FontSize',9);
ylim(ax,[0 5.5]*1e-3);
lgd = legend(ax,'SE Data','NE Data','SE Noise','NE Noise','Location','north',...
  'fontsize',7,'NumColumns',2);
%make background transparent
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));
text(ax,0.02,0.06,'c','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1);
xlabel(ax,sprintf('Projected \\Delta{t} (s)'));
hold(ax,'off');  
% keyboard  
  
%short-win data, 4-sta
ax = f.ax(4); 
impplt=imp4ths(:,7:8);
[ax,den1d4ths,conmat4ths,angle4th,anglegeo4th,c,~,conobj10th] = ...
  plt_cumulative_density_axis(ax,impplt,xran,yran,...
  [],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
hold(ax,'on');  
delete(conobj10th);
%a line cross (0,0) in the projection direction
yopt = linefcn(x,tand(angle4th(2)),0);
plot(ax,x,yopt,'-','linew',2,'color','b');
%a line cross (0,0) in the orthogonal direction
yort = linefcn(x,tand(angle4th(1)),0);
plot(ax,x,yort,':','linew',2,'color','b');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.057, pos(3), 0.02];
caxis(ax,[0 1.9]);
text(ax,0.98,0.95,sprintf('4-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, data'),'Units','normalized',...
  'FontSize',9);
% text(ax,0.98,0.12,strcat(num2str(round(anglegeo4th(2))),'$^{\,\circ}$',{'; '},...
%   num2str(round(anglegeo4th(1))),'$^{\,\circ}$'),'FontSize',11,...
%   'unit','normalized','interpreter','latex','HorizontalAlignment','right');
text(ax,0.98,0.12,sprintf('%d%c; %d%c',round(anglegeo4th(2)),char(176),...
  round(anglegeo4th(1)),char(176)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',10);
text(ax,0.02,0.06,'d','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1);
hold(ax,'off');  

%short-win noise, 4-sta
ax = f.ax(5); 
impplt=impn4ths(:,7:8);
[ax,den1dn4ths,conmatn4ths,anglen4th,anglegeon4th,c,~,conobj10th] = ...
  plt_cumulative_density_axis(ax,impplt,xran,yran,...
  [],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
hold(ax,'on');  
delete(conobj10th);
plot(ax,x,yopt,'-','linew',2,'color','r');
plot(ax,x,yort,':','linew',2,'color','r');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.057, pos(3), 0.02];
caxis(ax,[0 1.9]);
xlabel(ax,'');
ylabel(ax,'');
nolabels(ax,3);
text(ax,0.98,0.95,sprintf('4-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, noise'),'Units','normalized',...
  'FontSize',9);
% text(ax,0.98,0.12,strcat(num2str(round(anglegeo4th(2))),'$^{\,\circ}$',{'; '},...
%   num2str(round(anglegeo4th(1))),'$^{\,\circ}$'),'FontSize',11,...
%   'unit','normalized','interpreter','latex','HorizontalAlignment','right');
text(ax,0.98,0.12,sprintf('%d%c; %d%c',round(anglegeo4th(2)),char(176),...
  round(anglegeo4th(1)),char(176)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',10);
text(ax,0.02,0.06,'e','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1);
hold(ax,'off');  

%short-win cross section profile, data and noise, 4-sta 
ax = f.ax(6); 
c2pos = c.Position; 
f.ax(6).Position = [c2pos(1)+c2pos(3)+0.05 c2pos(2)+0.035 0.27 0.321];
normalizer=length(imp4ths);
hold(ax,'on');  
F = scatteredInterpolant(conmat4ths(:,3),conmat4ths(:,4),10.^(conmat4ths(:,1)),'linear','none');
Fn = scatteredInterpolant(conmatn4ths(:,3),conmatn4ths(:,4),10.^(conmatn4ths(:,1)),'linear','none');
[ax,avgopt,stdopt,muopt,sigmaopt,mdistprojopt]=plt_loccrssect(ax,F,x,yopt,anglegeo4th(2),...
  'b','-',xran/sps,which2plt,normalizer);
[ax,avgort,stdort,muort,sigmaort,mdistprojort]=plt_loccrssect(ax,F,x,yort,anglegeo4th(1),...
  'b',':',xran/sps,which2plt,normalizer);
[ax,avgoptn,stdoptn,muoptn,sigmaoptn,mdistprojoptn]=plt_loccrssect(ax,Fn,x,yopt,anglegeo4th(2),...
  'r','-',xran/sps,which2plt,normalizer);
[ax,avgortn,stdortn,muortn,sigmaortn,mdistprojortn]=plt_loccrssect(ax,Fn,x,yort,anglegeo4th(1),...
  'r',':',xran/sps,which2plt,normalizer);
text(ax,0.01,0.7,sprintf('SE \\mu=%.3f;\nSE \\sigma=%.3f',muopt,sigmaopt),'Units',...
  'normalized','HorizontalAlignment','left','FontSize',8,'color','b');
text(ax,0.01,0.55,sprintf('NE \\mu=%.3f;\nNE \\sigma=%.3f',muort,sigmaort),'Units',...
  'normalized','HorizontalAlignment','left','FontSize',8,'color','b');
text(ax,0.65,0.7,sprintf('SE \\mu=%.3f;\nSE \\sigma=%.3f',muoptn,sigmaoptn),'Units',...
  'normalized','HorizontalAlignment','left','FontSize',8,'color','r');
text(ax,0.65,0.55,sprintf('NE \\mu=%.3f;\nNE \\sigma=%.3f',muortn,sigmaortn),'Units',...
  'normalized','HorizontalAlignment','left','FontSize',8,'color','r');%,'FontName','Monospaced'
ylim(ax,[0 5.5]*1e-3);
lgd = legend(ax,'SE Data','NE Data','SE Noise','NE Noise','Location','north',...
  'fontsize',7,'NumColumns',2);
%make background transparent
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));
text(ax,0.02,0.06,'f','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1);
xlabel(ax,sprintf('Projected \\Delta{t} (s)'));
hold(ax,'off');  

fname = 'lfecatcompsftshort.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

keyboard

%% MAP space, cumulative density, directly transform density from SHIFTED sample space, short-win only
nrow = 2;
ncol = 2;
widin = 2.5*ncol;  % maximum width allowed is 8.5 inches
htin = 6.8;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.06 0.98]; pltyran = [0.05 0.98]; % optimal axis location
pltxsep = 0.02; pltysep = 0.03;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
 
xran = [-4 4];  %in km
yran = [-4 4];
% dx = 0.025;
% dy = 0.025;
binmethod = 'pixel';
marker = 'o';
msize = 5;
disttype = 'km';
contourflag=1;
smoothsigma=[];
scale = 'log10';

%%%3-station, Short-win, data
ax = f.ax(1); 
%convert contour lines from offset (s) to location (km) 
conmatsxy = contouroff2space(conmats,sps,ftrans); 
%note the density matrix is still in samples
den = den1ds;
tmploc = off2space002(den1ds(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
den(:,1:2) = tmploc(:,1:2);
ax = plt_den1d_axis(ax,den,conmatsxy,...
  xran,yran,[],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
caxis(ax,[0 1.9]);
% plot(ax,xcut,ycut,'k-','linew',2);
text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, data'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'a','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');

%%%3-station, Short-win, noise
ax = f.ax(2); 
conmatnsxy = contouroff2space(conmatns,sps,ftrans); 
den = den1dns;
tmploc = off2space002(den1dns(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
den(:,1:2) = tmploc(:,1:2);
ax = plt_den1d_axis(ax,den,conmatnsxy,...
  xran,yran,[],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
caxis(ax,[0 1.9]);
% plot(ax,xcut,ycut,'k-','linew',2);
xlabel(ax,'');
ylabel(ax,'');
nolabels(ax,3);
text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, noise'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'b','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');

%%%4-station, Short-win, data
ax = f.ax(3); 
conmat4thsxy = contouroff2space(conmat4ths,sps,ftrans);
den = den1d4ths;
tmploc = off2space002(den1d4ths(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
den(:,1:2) = tmploc(:,1:2);
ax = plt_den1d_axis(ax,den,conmat4thsxy,...
  xran,yran,[],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
caxis(ax,[0 1.9]);
% plot(ax,xcut,ycut,'k-','linew',2);
text(ax,0.98,0.95,sprintf('4-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, data'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'c','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');

%%%4-station, Short-win, noise
ax = f.ax(4); 
conmatn4thsxy = contouroff2space(conmatn4ths,sps,ftrans); 
den = den1dn4ths;
tmploc = off2space002(den1dn4ths(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
den(:,1:2) = tmploc(:,1:2);
ax = plt_den1d_axis(ax,den,conmatn4thsxy,...
  xran,yran,[],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
caxis(ax,[0 1.9]);
% plot(ax,xcut,ycut,'k-','linew',2);
xlabel(ax,'');
ylabel(ax,'');
nolabels(ax,3);
text(ax,0.98,0.95,sprintf('4-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, noise'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'d','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');


%% SAMPLE space, cumulative density, binning by pixel, ALL catalogs
nrow = 2;
ncol = 4;
widin = 2.5*ncol;  % maximum width allowed is 8.5 inches
htin = 6.8;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.06 0.98]; pltyran = [0.05 0.98]; % optimal axis location
pltxsep = 0.01; pltysep = 0.03;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
 
xran = [-40 40];  %in samples
yran = [-40 40];
dx = 1; dy = 1; 
binmethod = 'pixel';
marker = 'o';
msize = 6;
disttype = 'spl';
contourflag=1;
smoothsigma=1;
scale = 'log10';
ncont = 500;

%%%3-station, Short-win, data
ax = f.ax(1); 
impplt=imp(:,7:8);
[ax,den1d,conmat,~,~,c] = ...
  plt_cumulative_density_axis(ax,impplt,xran,yran,...
  [],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.055, pos(3), 0.02];
caxis(ax,[0 1.7]);
text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, data'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'a','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');

%%%3-station, Short-win, noise
ax = f.ax(2); 
impplt=impn(:,7:8);
[ax,den1dn,conmatn,~,~,c] = ...
  plt_cumulative_density_axis(ax,impplt,xran,yran,...
  [],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.055, pos(3), 0.02];
caxis(ax,[0 1.7]);
xlabel(ax,'');
ylabel(ax,'');
nolabels(ax,3);
text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, noise'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'b','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');

%%%3-station, Whole-win, data
ax = f.ax(3); 
impplt=imp1w(:,7:8);
[ax,den1d1w,conmat1w,~,~,c] = ...
  plt_cumulative_density_axis(ax,impplt,xran,yran,...
  [],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.055, pos(3), 0.02];
caxis(ax,[0 1.7]);
xlabel(ax,'');
ylabel(ax,'');
nolabels(ax,3);
text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Whole-win, data'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'c','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');

%%%3-station, Whole-win, noise 
ax = f.ax(4); 
impplt=impn1w(:,7:8);
[ax,den1dn1w,conmatn1w,~,~,c] = ...
  plt_cumulative_density_axis(ax,impplt,xran,yran,...
  [],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.055, pos(3), 0.02];
caxis(ax,[0 1.7]);
xlabel(ax,'');
ylabel(ax,'');
nolabels(ax,3);
text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Whole-win, noise'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'d','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');

%%%4-station, Short-win, data
ax = f.ax(5); 
impplt=imp4th(:,7:8);
[ax,den1d4th,conmat4th,~,~,c] = ...
  plt_cumulative_density_axis(ax,impplt,xran,yran,...
  [],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.055, pos(3), 0.02];
caxis(ax,[0 1.7]);
text(ax,0.98,0.95,sprintf('4-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, data'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'e','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');

%%%4-station, Short-win, noise
ax = f.ax(6); 
impplt=impn4th(:,7:8);
[ax,den1dn4th,conmatn4th,~,~,c] = ...
  plt_cumulative_density_axis(ax,impplt,xran,yran,...
  [],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.055, pos(3), 0.02];
caxis(ax,[0 1.7]);
xlabel(ax,'');
ylabel(ax,'');
nolabels(ax,3);
text(ax,0.98,0.95,sprintf('4-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, noise'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'f','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');

%%%4-station, Whole-win, data
ax = f.ax(7); 
impplt=imp1w4th(:,7:8);
[ax,den1d1w4th,conmat1w4th,~,~,c] = ...
  plt_cumulative_density_axis(ax,impplt,xran,yran,...
  [],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.055, pos(3), 0.02];
caxis(ax,[0 1.7]);
xlabel(ax,'');
ylabel(ax,'');
nolabels(ax,3);
text(ax,0.98,0.95,sprintf('4-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Whole-win, data'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'g','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');

%%%4-station, Whole-win, noise
ax = f.ax(8); 
impplt=impn1w4th(:,7:8);
[ax,den1dn1w4th,conmatn1w4th,~,~,c] = ...
  plt_cumulative_density_axis(ax,impplt,xran,yran,...
  [],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.055, pos(3), 0.02];
caxis(ax,[0 1.7]);
xlabel(ax,'');
ylabel(ax,'');
nolabels(ax,3);
text(ax,0.98,0.95,sprintf('4-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Whole-win, noise'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'h','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
% keyboard

orient(f.fig,'landscape');
fname = 'lfecatcompall.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
% keyboard


%% MAP space, cumulative density, directly transform density from sample space, ALL catalogs
nrow = 2;
ncol = 4;
widin = 2.5*ncol;  % maximum width allowed is 8.5 inches
htin = 6.8;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.06 0.98]; pltyran = [0.05 0.98]; % optimal axis location
pltxsep = 0.01; pltysep = 0.03;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
 
xran = [-4 4];  %in km
yran = [-4 4];
% dx = 0.025;
% dy = 0.025;
binmethod = 'pixel';
marker = 'o';
msize = 5;
disttype = 'km';
contourflag=1;
smoothsigma=[];
scale = 'log10';

%%%3-station, Short-win, data
ax = f.ax(1); 
%convert contour lines from offset (s) to location (km) 
conmatxy = contouroff2space(conmat,sps,ftrans); 
%note the density matrix is still in samples
den = den1d;
tmploc = off2space002(den1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
den(:,1:2) = tmploc(:,1:2);
ax = plt_den1d_axis(ax,den,conmatxy,...
  xran,yran,[],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
caxis(ax,[0 1.7]);
plot(ax,xcut,ycut,'k-','linew',2);
text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, data'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'a','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');

%%%3-station, Short-win, noise
ax = f.ax(2); 
conmatnxy = contouroff2space(conmatn,sps,ftrans); 
den = den1dn;
tmploc = off2space002(den1dn(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
den(:,1:2) = tmploc(:,1:2);
ax = plt_den1d_axis(ax,den,conmatnxy,...
  xran,yran,[],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
caxis(ax,[0 1.7]);
plot(ax,xcut,ycut,'k-','linew',2);
xlabel(ax,'');
ylabel(ax,'');
nolabels(ax,3);
text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, noise'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'b','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');

%%%3-station, Whole-win, data
ax = f.ax(3); 
conmat1wxy = contouroff2space(conmat1w,sps,ftrans);
den = den1d1w;
tmploc = off2space002(den1d1w(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
den(:,1:2) = tmploc(:,1:2);
ax = plt_den1d_axis(ax,den,conmat1wxy,...
  xran,yran,[],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
caxis(ax,[0 1.7]);
plot(ax,xcut,ycut,'k-','linew',2);
xlabel(ax,'');
ylabel(ax,'');
nolabels(ax,3);
text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Whole-win, data'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'c','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');

%%%3-station, Whole-win, noise 
ax = f.ax(4); 
conmatn1wxy = contouroff2space(conmatn1w,sps,ftrans);
den = den1dn1w;
tmploc = off2space002(den1dn1w(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
den(:,1:2) = tmploc(:,1:2);
ax = plt_den1d_axis(ax,den,conmatn1wxy,...
  xran,yran,[],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
caxis(ax,[0 1.7]);
plot(ax,xcut,ycut,'k-','linew',2);
xlabel(ax,'');
ylabel(ax,'');
nolabels(ax,3);
text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Whole-win, noise'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'d','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');

%%%4-station, Short-win, data
ax = f.ax(5); 
conmat4thxy = contouroff2space(conmat4th,sps,ftrans);
den = den1d4th;
tmploc = off2space002(den1d4th(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
den(:,1:2) = tmploc(:,1:2);
ax = plt_den1d_axis(ax,den,conmat4thxy,...
  xran,yran,[],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
caxis(ax,[0 1.7]);
plot(ax,xcut,ycut,'k-','linew',2);
text(ax,0.98,0.95,sprintf('4-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, data'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'e','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');

%%%4-station, Short-win, noise
ax = f.ax(6); 
conmatn4thxy = contouroff2space(conmatn4th,sps,ftrans); 
den = den1dn4th;
tmploc = off2space002(den1dn4th(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
den(:,1:2) = tmploc(:,1:2);
ax = plt_den1d_axis(ax,den,conmatn4thxy,...
  xran,yran,[],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
caxis(ax,[0 1.7]);
plot(ax,xcut,ycut,'k-','linew',2);
xlabel(ax,'');
ylabel(ax,'');
nolabels(ax,3);
text(ax,0.98,0.95,sprintf('4-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, noise'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'f','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');

%%%4-station, Whole-win, data
ax = f.ax(7); 
conmat1w4thxy = contouroff2space(conmat1w4th,sps,ftrans);
den = den1d1w4th;
tmploc = off2space002(den1d1w4th(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
den(:,1:2) = tmploc(:,1:2);
ax = plt_den1d_axis(ax,den,conmat1w4thxy,...
  xran,yran,[],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
caxis(ax,[0 1.7]);
plot(ax,xcut,ycut,'k-','linew',2);
xlabel(ax,'');
ylabel(ax,'');
nolabels(ax,3);
text(ax,0.98,0.95,sprintf('4-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Whole-win, data'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'g','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');

%%%4-station, Whole-win, noise
ax = f.ax(8);   
conmatn1w4thxy = contouroff2space(conmatn1w4th,sps,ftrans);
den = den1dn1w4th;
tmploc = off2space002(den1dn1w4th(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
den(:,1:2) = tmploc(:,1:2);
ax = plt_den1d_axis(ax,den,conmatn1w4thxy,...
  xran,yran,[],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
caxis(ax,[0 1.7]);
plot(ax,xcut,ycut,'k-','linew',2);
xlabel(ax,'');
ylabel(ax,'');
nolabels(ax,3);
text(ax,0.98,0.95,sprintf('4-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Whole-win, noise'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'h','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');

orient(f.fig,'landscape');
fname = 'lfecatcompall_map.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
% keyboard


%% MAP space, cumulative density, directly transform density from sample space, short-win only
nrow = 2;
ncol = 2;
widin = 2.5*ncol;  % maximum width allowed is 8.5 inches
htin = 6.8;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.06 0.98]; pltyran = [0.05 0.98]; % optimal axis location
pltxsep = 0.02; pltysep = 0.03;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
 
xran = [-4 4];  %in km
yran = [-4 4];
% dx = 0.025;
% dy = 0.025;
binmethod = 'pixel';
marker = 'o';
msize = 5;
disttype = 'km';
contourflag=1;
smoothsigma=[];
scale = 'log10';

%%%3-station, Short-win, data
ax = f.ax(1); 
%convert contour lines from offset (s) to location (km) 
conmatxy = contouroff2space(conmat,sps,ftrans); 
%note the density matrix is still in samples
den = den1d;
tmploc = off2space002(den1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
den(:,1:2) = tmploc(:,1:2);
ax = plt_den1d_axis(ax,den,conmatxy,...
  xran,yran,[],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
caxis(ax,[0 1.7]);
plot(ax,xcut,ycut,'k-','linew',2);
text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, data'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'a','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');

%%%3-station, Short-win, noise
ax = f.ax(2); 
conmatnxy = contouroff2space(conmatn,sps,ftrans); 
den = den1dn;
tmploc = off2space002(den1dn(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
den(:,1:2) = tmploc(:,1:2);
ax = plt_den1d_axis(ax,den,conmatnxy,...
  xran,yran,[],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
caxis(ax,[0 1.7]);
plot(ax,xcut,ycut,'k-','linew',2);
xlabel(ax,'');
ylabel(ax,'');
nolabels(ax,3);
text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, noise'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'b','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');

%%%4-station, Short-win, data
ax = f.ax(3); 
conmat4thxy = contouroff2space(conmat4th,sps,ftrans);
den = den1d4th;
tmploc = off2space002(den1d4th(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
den(:,1:2) = tmploc(:,1:2);
ax = plt_den1d_axis(ax,den,conmat4thxy,...
  xran,yran,[],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
caxis(ax,[0 1.7]);
plot(ax,xcut,ycut,'k-','linew',2);
text(ax,0.98,0.95,sprintf('4-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, data'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'c','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');

%%%4-station, Short-win, noise
ax = f.ax(4); 
conmatn4thxy = contouroff2space(conmatn4th,sps,ftrans); 
den = den1dn4th;
tmploc = off2space002(den1dn4th(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
den(:,1:2) = tmploc(:,1:2);
ax = plt_den1d_axis(ax,den,conmatn4thxy,...
  xran,yran,[],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
caxis(ax,[0 1.7]);
plot(ax,xcut,ycut,'k-','linew',2);
xlabel(ax,'');
ylabel(ax,'');
nolabels(ax,3);
text(ax,0.98,0.95,sprintf('4-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,sprintf('Short-win, noise'),'Units','normalized',...
  'FontSize',9);
text(ax,0.02,0.06,'d','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');

% % orient(f.fig,'landscape');
fname = 'lfecatcompshort_map.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));


% %% MAP space, cumulative density, binning by pixel 
% nrow = 2;
% ncol = 3;
% widin = 2.5*ncol;  % maximum width allowed is 8.5 inches
% htin = 7;   % maximum height allowed is 11 inches
% f = initfig(widin,htin,nrow,ncol);

% pltxran = [0.06 0.98]; pltyran = [0.05 0.98]; % optimal axis location
% pltxsep = 0.01; pltysep = 0.03;
% optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

% xran = [-4 4];  %in km
% yran = [-4 4];
% binmethod = 'pixel';
% marker = 'o';
% msize = 4;
% disttype = 'km';
% contourflag=1;
% smoothsigma=1;
% scale = 'log10';
% ncont = 500;
% which2plt = 'data'; %'gsfit'

% %short-win data, 3-sta
% ax = f.ax(1);
% %convert contour lines from offset (s) to location (km) 
% conmatxy = contouroff2space(conmat,sps,ftrans); 
% %note the density matrix is still in samples
% den = den1d;
% tmploc = off2space002(den1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
% den(:,1:2) = tmploc(:,1:2);
% ax = plt_den1d_axis(ax,den,conmatxy,xran,yran,...
%   [],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
% hold(ax,'on');  
% %Principal component analysis
% [coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(imploc(:,1:2));
% dxln = 0.01; %in sec
% x = reshape(xran(1):dxln:xran(2), [], 1);
% %a line cross (0,0) in the projection direction
% yopt = linefcn(x,tand(angle(2)),0);
% plot(ax,x,yopt,'-','linew',2,'color','b');
% %a line cross (0,0) in the orthogonal direction
% yort = linefcn(x,tand(angle(1)),0);
% plot(ax,x,yort,':','linew',2,'color','b');
% %plot the high-density 4-s tremor ellipse on top 
% plot(ax,xcut,ycut,'k-','LineWidth',2);
% % pos = ax.Position;
% % c.Position = [pos(1), pos(2)-0.055, pos(3), 0.02];
% text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
%   'HorizontalAlignment','right','FontSize',9);
% text(ax,0.02,0.95,sprintf('Short-win, data'),'Units','normalized',...
%   'FontSize',9);
% text(ax,0.98,0.12,strcat(num2str(round(anglegeo(2))),'$^{\,\circ}$',{'; '},...
%   num2str(round(anglegeo(1))),'$^{\,\circ}$'),'FontSize',11,...
%   'unit','normalized','interpreter','latex','HorizontalAlignment','right');

% hold(ax,'off');  

% %short-win noise, 3-sta
% ax = f.ax(2);
% %convert contour lines from offset (s) to location (km) 
% conmatxyn = contouroff2space(conmatn,sps,ftrans); 
% %note the density matrix is still in samples
% den = den1dn;
% tmploc = off2space002(den1dn(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
% den(:,1:2) = tmploc(:,1:2);
% ax = plt_den1d_axis(ax,den,conmatxyn,xran,yran,...
%   [],[],binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale);
% hold(ax,'on');  
% plot(ax,x,yopt,'-','linew',2,'color','r');
% plot(ax,x,yort,':','linew',2,'color','r');
% %plot the high-density 4-s tremor ellipse on top 
% plot(ax,xcut,ycut,'k-','LineWidth',2);
% % pos = ax.Position;
% % c.Position = [pos(1), pos(2)-0.055, pos(3), 0.02];
% xlabel(ax,'');
% ylabel(ax,'');
% nolabels(ax,3);
% text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
%   'HorizontalAlignment','right','FontSize',9);
% text(ax,0.02,0.95,sprintf('Short-win, noise'),'Units','normalized',...
%   'FontSize',9);
% text(ax,0.98,0.12,strcat(num2str(round(anglegeo(2))),'$^{\,\circ}$',{'; '},...
%   num2str(round(anglegeo(1))),'$^{\,\circ}$'),'FontSize',11,...
%   'unit','normalized','interpreter','latex','HorizontalAlignment','right');

% hold(ax,'off');  

% %short-win cross section profile, data and noise, 3-sta 
% ax = f.ax(3);
% ax1pos = f.ax(2).Position; 
% f.ax(3).Position = [ax1pos(1)+ax1pos(3)+0.05 ax1pos(2)+0.033 0.27 0.32];
% normalizer=length(imp);
% hold(ax,'on');  
% %contour lines are log10 scale, now revert to linear scale
% F = scatteredInterpolant(conmatxy(:,3),conmatxy(:,4),10.^(conmatxy(:,1)),'linear','none');
% Fn = scatteredInterpolant(conmatxyn(:,3),conmatxyn(:,4),10.^(conmatxyn(:,1)),'linear','none');
% [ax,avgopt,stdopt,muopt,sigmaopt,mdistprojopt]=plt_loccrssect(ax,F,x,yopt,anglegeo(2),...
%   'b','-',xran,which2plt,normalizer);
% [ax,avgort,stdort,muort,sigmaort,mdistprojort]=plt_loccrssect(ax,F,x,yort,anglegeo(1),...
%   'b',':',xran,which2plt,normalizer);
% %use the SAME projection direction from data to noise
% [ax,avgoptn,stdoptn,muoptn,sigmaoptn,mdistprojoptn]=plt_loccrssect(ax,Fn,x,yopt,anglegeo(2),...
%   'r','-',xran,which2plt,normalizer);
% [ax,avgortn,stdortn,muortn,sigmaortn,mdistprojortn]=plt_loccrssect(ax,Fn,x,yort,anglegeo(1),...
%   'r',':',xran,which2plt,normalizer);
% text(ax,0.01,0.75,sprintf('SE \\mu=%.2f;\nSE \\sigma=%.2f',muopt,sigmaopt),'Units',...
%   'normalized','HorizontalAlignment','left','FontSize',8,'color','b');
% text(ax,0.01,0.6,sprintf('NE \\mu=%.2f;\nNE \\sigma=%.2f',muort,sigmaort),'Units',...
%   'normalized','HorizontalAlignment','left','FontSize',8,'color','b');
% text(ax,0.65,0.75,sprintf('SE \\mu=%.2f;\nSE \\sigma=%.2f',muoptn,sigmaoptn),'Units',...
%   'normalized','HorizontalAlignment','left','FontSize',8,'color','r');
% text(ax,0.65,0.6,sprintf('NE \\mu=%.2f;\nNE \\sigma=%.2f',muortn,sigmaortn),'Units',...
%   'normalized','HorizontalAlignment','left','FontSize',8,'color','r');
% ylim(ax,[0 2.5]*1e-3);
% lgd = legend(ax,'SE Data','NE Data','SE Noise','NE Noise','Location','northwest',...
%   'fontsize',7,'NumColumns',2);
% %make background transparent
% set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));
% xlabel(ax,sprintf('Projected \\Delta{t} (s)'));
% hold(ax,'off');  
% keyboard  

% %% in MAP view, plot cumulative density, binning by a grid 
% nrow = 2;
% ncol = 3;
% widin = 8.3;  % maximum width allowed is 8.5 inches
% htin = 7;   % maximum height allowed is 11 inches
% f = initfig(widin,htin,nrow,ncol);

% pltxran = [0.08 0.98]; pltyran = [0.05 0.98]; % optimal axis location
% pltxsep = 0.02; pltysep = 0.06;
% optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
 
% xran = [-4 4];  %in km
% yran = [-4 4];
% dx = 0.025;
% dy = 0.025;
% binmethod = 'grid';
% marker = 'o';
% msize = 6;
% disttype = 'km';
% contourflag=0;
% smoothsigma=4;
% scale = 'log10';

% ax = f.ax(1); 
% impplt=imploc(:,1:2);
% ax = plt_cumulative_density_axis(ax,impplt,xran,yran,...
%   dx,dy,binmethod,marker,msize,disttype,contourflag,smoothsigma,scale);
% text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
%   'HorizontalAlignment','right','FontSize',9);
% text(ax,0.02,0.95,sprintf('Short-win, data'),'Units','normalized',...
%   'FontSize',9);

% ax = f.ax(4); 
% impplt=imploc4th(:,1:2);
% ax = plt_cumulative_density_axis(ax,impplt,xran,yran,...
%   dx,dy,binmethod,marker,msize,disttype,contourflag,smoothsigma,scale);
% text(ax,0.98,0.95,sprintf('4-station'),'Units','normalized',...
%   'HorizontalAlignment','right','FontSize',9);
% text(ax,0.02,0.95,sprintf('Short-win, data'),'Units','normalized',...
%   'FontSize',9);


% ax = f.ax(2); 
% impplt=imploc1w(:,1:2);
% ax = plt_cumulative_density_axis(ax,impplt,xran,yran,...
%   dx,dy,binmethod,marker,msize,disttype,contourflag,smoothsigma,scale);
% xlabel(ax,'');
% ylabel(ax,'');
% nolabels(ax,3);
% text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
%   'HorizontalAlignment','right','FontSize',9);
% text(ax,0.02,0.95,sprintf('Whole-win, data'),'Units','normalized',...
%   'FontSize',9);

% ax = f.ax(5); 
% impplt=imploc1w4th(:,1:2);
% ax = plt_cumulative_density_axis(ax,impplt,xran,yran,...
%   dx,dy,binmethod,marker,msize,disttype,contourflag,smoothsigma,scale);
% xlabel(ax,'');
% ylabel(ax,'');
% nolabels(ax,3);
% text(ax,0.98,0.95,sprintf('4-station'),'Units','normalized',...
%   'HorizontalAlignment','right','FontSize',9);
% text(ax,0.02,0.95,sprintf('Whole-win, data'),'Units','normalized',...
%   'FontSize',9);


% ax = f.ax(3); 
% impplt=implocn1w(:,1:2);
% ax = plt_cumulative_density_axis(ax,impplt,xran,yran,...
%   dx,dy,binmethod,marker,msize,disttype,contourflag,smoothsigma,scale);
% xlabel(ax,'');
% ylabel(ax,'');
% nolabels(ax,3);
% text(ax,0.98,0.95,sprintf('3-station'),'Units','normalized',...
%   'HorizontalAlignment','right','FontSize',9);
% text(ax,0.02,0.95,sprintf('Whole-win, noise'),'Units','normalized',...
%   'FontSize',9);

% ax = f.ax(6); 
% impplt=implocn1w4th(:,1:2);
% ax = plt_cumulative_density_axis(ax,impplt,xran,yran,...
%   dx,dy,binmethod,marker,msize,disttype,contourflag,smoothsigma,scale);
% xlabel(ax,'');
% ylabel(ax,'');
% nolabels(ax,3);
% text(ax,0.98,0.95,sprintf('4-station'),'Units','normalized',...
%   'HorizontalAlignment','right','FontSize',9);
% text(ax,0.02,0.95,sprintf('Whole-win, noise'),'Units','normalized',...
%   'FontSize',9);

% % orient(f.fig,'landscape');
% fname = 'lfecatcomparison_mapgrid.pdf';
% print(f.fig,'-dpdf',...
%   strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));





