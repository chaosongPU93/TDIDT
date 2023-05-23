% compare_lfetremorcats.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we have the catalog of the LFEs from every burst win of interest. We 
% also have the 4-s tremor catalog from myself, and the LFE catalog from 
% Michael Bostock (although his locations might be less useful). Then we 
% could have a comparison of both the cumulative heat map that shows the 
% active spots, and the single map of comparing LFEs and tremor locations
% and see how close they are. I also have a script of determining the 4-s
% tremor migration times using the automated algorithm, which might give
% a larger scale of tremor background. Comparison can also be time-distance
% plot.
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/01/30
% Last modified date:   2023/01/30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
format short e   % Set the format to 5-digit floating point
clear
clc
close all

%% for easy testing
defval('idxbst',1:195); %global indices of bursts to run 

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
%format: [fam yyyy mm dd sec dx dy lon lat dep magnitude number-of-stations], 12 cols
%time should point to the peak, zero-crossing?
bostcat = ReformBostock(loc0(3),loc0(4),0);

%the cut-out boundary of 4-s detections
cutout = 'ellipse';
x0 = 0.2; 
y0 = 0.2;
semia = 1.75;
semib = 1.25;
angrot = 45;
[xcut, ycut] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);

%same rectangle in 'locinterp002_4s.m'
EW = [-7 3];
NS = [-3 4];
wid = range(EW);
hgt = range(NS);
x0 = mean(EW);
y0 = mean(NS);
[x, y] = rectangle_chao(x0,y0,wid,hgt,0.01);

%%%2022/06/29, use the same ellipse to exclude fam 047
% bnd = [x y];
bnd = [xcut ycut];
[iin,ion] = inpolygon(bostcat(:,6),bostcat(:,7),bnd(:,1),bnd(:,2));
isinbnd = iin | ion;
bostcati = bostcat(isinbnd == 1, :);
bostcato = bostcat(isinbnd ~= 1, :);
clear bostcat

%%%load the tremor catalog of John Armbruster, 
%%%2022/06/29, not really very useful as the detecting window could be 128-s long (not sure)
%format: [yyyy mm dd sec dx dy lon lat dep], 9 cols;
armcat = ReformArmbrusterv2(loc0(3),loc0(4),0);
[~,ind] = min(abs(armcat(:,5))+abs(armcat(:,6)));
% armcat1 = ReformArmbruster(loc0(3),loc0(4),1);
% [~,ind1] = min(abs(armcat1(:,5))+abs(armcat1(:,6)));
[iin,ion] = inpolygon(armcat(:,5),armcat(:,6),bnd(:,1),bnd(:,2));
isinbnd = iin | ion;
armcati = armcat(isinbnd == 1, :);
armcato = armcat(isinbnd ~= 1, :);
clear armcat

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

%% Comparison of cumulative density in relative locations
%plot the cumulative density map, binning by pixel, ie., each unique detection         
xran = [-6 6];
yran = [-6 6];
[f] = plt_cumulative_density(imploc,hfall,xran,yran,'pixel',10,10);          
% supertit(f.ax, 'Density binned by pixel');
text(f.ax(1),0.95,0.95,'LFE, 25-s-win, 2ndary removed','Units','normalized',...
'HorizontalAlignment','right','FontSize',9);
hold(f.ax(2),'on');
plot(f.ax(2),xcut,ycut,'k-','LineWidth',1.5);
text(f.ax(2),0.95,0.95,'Tremor','Units','normalized',...
'HorizontalAlignment','right','FontSize',9);
% caxis(f2.ax,[0 1.4]);

[f] = plt_cumulative_density(imploc,imploc4th,xran,yran,'pixel',10,10);          
% supertit(f.ax, 'Density binned by pixel');
text(f.ax(1),0.95,0.95,'LFE, 25-s-win, 2ndary removed','Units','normalized',...
'HorizontalAlignment','right','FontSize',9);
text(f.ax(2),0.95,0.95,'LFE, 25-s-win, checked at KLNB','Units','normalized',...
'HorizontalAlignment','right','FontSize',9);
% caxis(f2.ax,[0 1.4]);

[f] = plt_cumulative_density(imploc1win,imploc1win4th,xran,yran,'pixel',10,10);          
supertit(f.ax, 'Density binned by pixel');
text(f.ax(1),0.95,0.95,'LFE, whole-win, 2ndary removed','Units','normalized',...
'HorizontalAlignment','right','FontSize',9);
text(f.ax(2),0.95,0.95,'LFE, whole-win, checked at KLNB','Units','normalized',...
'HorizontalAlignment','right','FontSize',9);
% caxis(f2.ax,[0 1.4]);

[f] = plt_cumulative_density(implocn1win,implocn1win4th,xran,yran,'pixel',10,10);          
supertit(f.ax, 'Density binned by pixel');
text(f.ax(1),0.95,0.95,'LFE, whole-win, noise, 2ndary removed','Units','normalized',...
'HorizontalAlignment','right','FontSize',9);
text(f.ax(2),0.95,0.95,'LFE, whole-win, noise, checked at KLNB','Units','normalized',...
'HorizontalAlignment','right','FontSize',9);


%% Comparison of cumulative density in sample space
[f] = plt_cumulative_density(imp(:,7:8),hfall(:,9:10)*4,[-50 50],[-50 50],'pixel',6,6,[],[],0);          
text(f.ax(1),0.95,0.95,'LFE, 25-s-win, 2ndary removed','Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(f.ax(2),0.95,0.95,'Tremor','Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
xlabel(f.ax(1),'offset 12 (samples)');
ylabel(f.ax(1),'offset 13 (samples)');
xticks(f.ax(1),-50:10:50);
yticks(f.ax(1),-50:10:50);
xlabel(f.ax(2),'offset 12 (samples)');
ylabel(f.ax(2),'offset 13 (samples)');
xticks(f.ax(2),-50:10:50);
yticks(f.ax(2),-50:10:50);

[f] = plt_cumulative_density(imp(:,7:8),imp4th(:,7:8),[-50 50],[-50 50],'pixel',6,6,[],[],1);          
text(f.ax(1),0.95,0.95,'LFE, 25-s-win, 2ndary removed','Units','normalized',...
'HorizontalAlignment','right','FontSize',9);
text(f.ax(2),0.95,0.95,'LFE, 25-s-win, checked at KLNB','Units','normalized',...
'HorizontalAlignment','right','FontSize',9);
xlabel(f.ax(1),'offset 12 (samples)');
ylabel(f.ax(1),'offset 13 (samples)');
xticks(f.ax(1),-50:10:50);
yticks(f.ax(1),-50:10:50);
xlabel(f.ax(2),'offset 12 (samples)');
ylabel(f.ax(2),'offset 13 (samples)');
xticks(f.ax(2),-50:10:50);
yticks(f.ax(2),-50:10:50);

[f] = plt_cumulative_density(imp1win(:,7:8),imp1win4th(:,7:8),[-50 50],[-50 50],'pixel',6,6,[],[],1);          
text(f.ax(1),0.95,0.95,'LFE, whole-win, 2ndary removed','Units','normalized',...
'HorizontalAlignment','right','FontSize',9);
text(f.ax(2),0.95,0.95,'LFE, whole-win, checked at KLNB','Units','normalized',...
'HorizontalAlignment','right','FontSize',9);
xlabel(f.ax(1),'offset 12 (samples)');
ylabel(f.ax(1),'offset 13 (samples)');
xticks(f.ax(1),-50:10:50);
yticks(f.ax(1),-50:10:50);
xlabel(f.ax(2),'offset 12 (samples)');
ylabel(f.ax(2),'offset 13 (samples)');
xticks(f.ax(2),-50:10:50);
yticks(f.ax(2),-50:10:50);

[f] = plt_cumulative_density(impn1win(:,7:8),impn1win4th(:,7:8),[-50 50],[-50 50],'pixel',6,6,[],[],1);          
text(f.ax(1),0.95,0.95,'LFE, whole-win, noise, 2ndary removed','Units','normalized',...
'HorizontalAlignment','right','FontSize',9);
text(f.ax(2),0.95,0.95,'LFE, whole-win, noise, checked at KLNB','Units','normalized',...
'HorizontalAlignment','right','FontSize',9);
xlabel(f.ax(1),'offset 12 (samples)');
ylabel(f.ax(1),'offset 13 (samples)');
xticks(f.ax(1),-50:10:50);
yticks(f.ax(1),-50:10:50);
xlabel(f.ax(2),'offset 12 (samples)');
ylabel(f.ax(2),'offset 13 (samples)');
xticks(f.ax(2),-50:10:50);
yticks(f.ax(2),-50:10:50);

% impden = density_pixel(imp(:,7),imp(:,8));
% [xyzgridpad,xgrid,ygrid,zgrid,ind2] = zeropadmat2d(impden,-50:1:50,-50:1:50);
% zgridgf = imgaussfilt(zgrid, 1);  %smooth it a bit
% hold(f.ax(1),'on');
% conmat = contour(f.ax(1),xgrid,ygrid,log10(zgridgf),0.1:0.2:1.5,'-','color',[.3 .3 .3],...
%   'ShowText','on'); %
% impden = density_pixel(hfall(:,9)*4,hfall(:,10)*4);
% [xyzgridpad,xgrid,ygrid,zgrid,ind2] = zeropadmat2d(impden,-50:1:50,-50:1:50);
% zgridgf = imgaussfilt(zgrid, 1);  %smooth it a bit
% hold(f.ax(2),'on');
% conmat = contour(f.ax(2),xgrid,ygrid,log10(zgridgf),0.1:0.2:1.5,'-','color',[.3 .3 .3],...
%   'ShowText','on'); %


%%
[f] = plt_cumulative_density(imp1win(:,7:8),hfbnd(:,9:10)*4,[-50 50],[-50 50],'pixel',6,6);          
supertit(f.ax, 'Density binned by pixel');
text(f.ax(1),0.95,0.95,'LFE, whole-win','Units','normalized',...
'HorizontalAlignment','right','FontSize',9);
text(f.ax(2),0.95,0.95,'Tremor','Units','normalized',...
'HorizontalAlignment','right','FontSize',9);
xlabel(f.ax(1),'offset 12');
ylabel(f.ax(1),'offset 13');
xlabel(f.ax(2),'offset 12');
ylabel(f.ax(2),'offset 13');
impden = density_pixel(imp1win(:,7),imp1win(:,8));
[xyzgridpad,xgrid,ygrid,zgrid,ind2] = zeropadmat2d(impden,-50:1:50,-50:1:50);
zgridgf = imgaussfilt(zgrid, 2);  %smooth it a bit
hold(f.ax(1),'on');
conmat = contour(f.ax(1),xgrid,ygrid,log10(zgridgf),0.1:0.2:1.5,'-','color',[.3 .3 .3],...
  'ShowText','on'); %
impden = density_pixel(hfbnd(:,9)*4,hfbnd(:,10)*4);
[xyzgridpad,xgrid,ygrid,zgrid,ind2] = zeropadmat2d(impden,-50:1:50,-50:1:50);
zgridgf = imgaussfilt(zgrid, 2);  %smooth it a bit
hold(f.ax(2),'on');
conmat = contour(f.ax(2),xgrid,ygrid,log10(zgridgf),0.1:0.2:1.5,'-','color',[.3 .3 .3],...
  'ShowText','on'); %

%%
[f] = plt_cumulative_density(impn(:,7:8),impn1win(:,7:8),[-50 50],[-50 50],'pixel',6,6);          
supertit(f.ax, 'Density binned by pixel');
text(f.ax(1),0.95,0.95,'LFE, 25-s-win, noise','Units','normalized',...
'HorizontalAlignment','right','FontSize',9);
text(f.ax(2),0.95,0.95,'LFE, whole-win, noise','Units','normalized',...
'HorizontalAlignment','right','FontSize',9);
xlabel(f.ax(1),'offset 12');
ylabel(f.ax(1),'offset 13');
xlabel(f.ax(2),'offset 12');
ylabel(f.ax(2),'offset 13');
impden = density_pixel(impn(:,7),impn(:,8));
[xyzgridpad,xgrid,ygrid,zgrid,ind2] = zeropadmat2d(impden,-50:1:50,-50:1:50);
zgridgf = imgaussfilt(zgrid, 2);  %smooth it a bit
hold(f.ax(1),'on');
conmat = contour(f.ax(1),xgrid,ygrid,log10(zgridgf),0.1:0.2:1.5,'-','color',[.3 .3 .3],...
  'ShowText','on'); %
impden = density_pixel(impn1win(:,7),impn1win(:,8));
[xyzgridpad,xgrid,ygrid,zgrid,ind2] = zeropadmat2d(impden,-50:1:50,-50:1:50);
zgridgf = imgaussfilt(zgrid, 2);  %smooth it a bit
hold(f.ax(2),'on');
conmat = contour(f.ax(2),xgrid,ygrid,log10(zgridgf),0.1:0.2:1.5,'-','color',[.3 .3 .3],...
  'ShowText','on'); %


%%
%%%plot the comparison for each single burst
nsrcraw = allsig.allbstsig.nsrcraw;
rccsrc = allsig.allbstsig.rcccatsrcall;
propang = allsig.allbstsig.propang;
proppear = allsig.allbstsig.proppear;

nsrcraw2 = allsig1win.allbstsig.nsrcraw;
rccsrc2 = allsig1win.allbstsig.rccsrcall;
propang2 = allsig1win.allbstsig.propang;
proppear2 = allsig1win.allbstsig.proppear;

nsrcrawn = allnoi.allbstnoi.nsrcraw;
rccsrcn = allnoi.allbstnoi.rcccatsrcall;
propangn = allnoi.allbstnoi.propang;
proppearn = allnoi.allbstnoi.proppear;

nsrcraw2n = allnoi1win.allbstnoi.nsrcraw;
rccsrc2n = allnoi1win.allbstnoi.rccsrcall;
propang2n = allnoi1win.allbstnoi.propang;
proppear2n = allnoi1win.allbstnoi.proppear;

clrcode = 'tarvl';
clrcode = 'rcc';

idxbst = 181;
idxbst = 1:length(trange);

dates = unique(trange(:,1));
years = unique(floor(dates/1000));
nets = length(years);
for iii = 1: length(idxbst)
  [iets,i,j] = indofburst(trange,idxbst(iii));

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

  %time and loc of LFE of that burst, from 25-win detection
  ist = sum(nsrcraw(1:(idxbst(iii)-1)))+1;
  ied = ist+nsrcraw(idxbst(iii))-1;
  lfe = imp(ist:ied, :);
  lfeloc = imploc(ist:ied, :);
  rcc = rccsrc(ist:ied, 1); %cat RCC at the average src arrival of 3 stas
  
  %time and loc of LFE of that burst, from whole-win detection
  ist = sum(nsrcraw2(1:(idxbst(iii)-1)))+1;
  ied = ist+nsrcraw2(idxbst(iii))-1;
  lfe2 = imp1win(ist:ied, :);
  lfeloc2 = imploc1win(ist:ied, :);
  rcc2 = rccsrc2(ist:ied, 1); %cat RCC at the average src arrival of 3 stas

  %time and loc of LFE of that burst, from 25-win detection, synthetic noise
  ist = sum(nsrcrawn(1:(idxbst(iii)-1)))+1;
  ied = ist+nsrcrawn(idxbst(iii))-1;
  lfen = impn(ist:ied, :);
  lfelocn = impnloc(ist:ied, :);
  rccn = rccsrcn(ist:ied, 1); %cat RCC at the average src arrival of 3 stas

  %time and loc of LFE of that burst, from whole-win detection, synthetic noise
  ist = sum(nsrcraw2n(1:(idxbst(iii)-1)))+1;
  ied = ist+nsrcraw2n(idxbst(iii))-1;
  lfe2n = impn1win(ist:ied, :);
  lfeloc2n = impnloc1win(ist:ied, :);
  rcc2n = rccsrc2n(ist:ied, 1); %cat RCC at the average src arrival of 3 stas

  %time and loc of tremor of that burst
  hfdayi = hfall(hfall(:,daycol)==trange(idxbst(iii),1), :);  % inside bound of the day 
  tmaxi = hfdayi(:, seccol); % starting time of max power rate of half sec inside the ellipse
  tcnti = hfdayi(:, 15);  % the center of detecting win is the 15th col
  %how many 4-s detections fall into the burst range 
  indtmaxi = find(tmaxi>=trange(idxbst(iii),2)-0.1 & tmaxi<=trange(idxbst(iii),3)+0.1);
  %Use the start and end of the 4-s detecting window
  tstbuf = min(tcnti(indtmaxi)-2);
  tedbuf = max(tcnti(indtmaxi)+2); 
  tlenbuf(idxbst(iii)) = tedbuf-tstbuf;
  tmrloc = hfdayi(hfdayi(:, seccol)>=tstbuf & hfdayi(:, seccol)<=tedbuf, :);

  if strcmp(clrcode, 'tarvl') 
    %%%plot
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 8;   % maximum height allowed is 11 inches
    nrow = 2;
    ncol = 2;
    f = initfig(widin,htin,nrow,ncol); %initialize fig
    xran = [0.06 0.96]; yran = [0.06 0.98];
  %   c.Label.String = clabel;
    caxis(ax,cran);
    axis(ax,'equal');
    [rotx, roty] = complex_rot(0,1,-propang(idxbst(iii)));
    xvect = [0.5-rotx 0.5+rotx];
    yvect = [-2.5-roty -2.5+roty];
    drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5);
    text(ax,0.95,0.1,strcat({'max-speed direc: '},num2str(propang(idxbst(iii))),{'{\circ}'}),...
      'FontSize',9,'unit','normalized','horizontalalignment','right');
    text(ax,0.95,0.05,sprintf('Pearson=%.2f',proppear(idxbst(iii))),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.95,'LFE, 25-s-win','Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    xlabel(ax,'E (km)');
    ylabel(ax,'N (km)');
    xticks(ax,xran(1):1:xran(2));
    yticks(ax,yran(1):1:yran(2));
    axis(ax,[xran yran]);
    longticks(ax,2);

    %LFE locations, whole-win detection
    ax=f.ax(2);
    hold(ax,'on');
    ax.Box='on'; grid(ax,'on');
    xran = [-5 5];
    yran = [-5 5];
    plot(ax,[xran(1) xran(2)],[0 0],'--','Color',[.3 .3 .3]);
    plot(ax,[0 0],[yran(1) yran(2)],'--','Color',[.3 .3 .3]);
    scatter(ax,lfeloc2(:,1),lfeloc2(:,2),15,lfe2(:,1)/sps,'filled',...
      'MarkerEdgeColor',[.5 .5 .5]);
    text(ax,0.95,0.85,sprintf('%d',size(lfeloc2,1)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    oldc = colormap(ax,'kelicol');
    newc = flipud(oldc);
    colormap(ax,newc);
    c=colorbar(ax);
    c.Label.String = clabel;
    caxis(ax,cran);
    axis(ax,'equal');
    [rotx, roty] = complex_rot(0,1,-propang2(idxbst(iii)));
    xvect = [0.5-rotx 0.5+rotx];
    yvect = [-2.5-roty -2.5+roty];
    drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5);
    text(ax,0.95,0.1,strcat({'max-speed direc: '},num2str(propang2(idxbst(iii))),{'{\circ}'}),...
      'FontSize',9,'unit','normalized','horizontalalignment','right');
    text(ax,0.95,0.05,sprintf('Pearson=%.2f',proppear2(idxbst(iii))),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.95,'LFE, whole-win','Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    xlabel(ax,'E (km)');
    ylabel(ax,'N (km)');
    xticks(ax,xran(1):1:xran(2));
    yticks(ax,yran(1):1:yran(2));
    axis(ax,[xran yran]);
    longticks(ax,2);

    %LFE synthetic noise, whole-win
    ax=f.ax(3);
    hold(ax,'on');
    ax.Box='on'; grid(ax,'on');
    xran = [-5 5];
    yran = [-5 5];
    plot(ax,[xran(1) xran(2)],[0 0],'--','Color',[.3 .3 .3]);
    plot(ax,[0 0],[yran(1) yran(2)],'--','Color',[.3 .3 .3]);
    scatter(ax,lfeloc2n(:,1),lfeloc2n(:,2),15,lfe2n(:,1)/sps,'filled',...
      'MarkerEdgeColor',[.5 .5 .5]);
    text(ax,0.95,0.85,sprintf('%d',size(lfeloc2n,1)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    oldc = colormap(ax,'kelicol');
    newc = flipud(oldc);
    colormap(ax,newc);
    c=colorbar(ax);
    c.Label.String = clabel;
    caxis(ax,cran);
    axis(ax,'equal');
    [rotx, roty] = complex_rot(0,1,-propang2n(idxbst(iii)));
    xvect = [0.5-rotx 0.5+rotx];
    yvect = [-2.5-roty -2.5+roty];
    drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5);
    text(ax,0.95,0.1,strcat({'max-speed direc: '},num2str(propang2n(idxbst(iii))),{'{\circ}'}),...
      'FontSize',9,'unit','normalized','horizontalalignment','right');
    text(ax,0.95,0.05,sprintf('Pearson=%.2f',proppear2n(idxbst(iii))),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.95,'LFE, whole-win, noise','Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    xlabel(ax,'E (km)');
    ylabel(ax,'N (km)');
    xticks(ax,xran(1):1:xran(2));
    yticks(ax,yran(1):1:yran(2));
    axis(ax,[xran yran]);
    longticks(ax,2);

    %tremor locations
    ax=f.ax(4);
    hold(ax,'on');
    ax.Box='on'; grid(ax,'on');
    plot(ax,[xran(1) xran(2)],[0 0],'--','Color',[.3 .3 .3]);
    plot(ax,[0 0],[yran(1) yran(2)],'--','Color',[.3 .3 .3]);
    scatter(ax,tmrloc(:,1),tmrloc(:,2),15,tmrloc(:,seccol)-tstbuf,'filled',...
      'MarkerEdgeColor',[.5 .5 .5]);
    text(ax,0.95,0.85,sprintf('%d',size(tmrloc,1)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.95,'4-s Tremor','Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    oldc = colormap(ax,'kelicol');
    newc = flipud(oldc);
    colormap(ax,newc);
    c=colorbar(ax);
    c.Label.String = sprintf('Tremor arrival time (s)');
    caxis(ax,cran);
    axis(ax,'equal');
    xlabel(ax,'E (km)');
    ylabel(ax,'N (km)');
    xticks(ax,xran(1):1:xran(2));
    yticks(ax,yran(1):1:yran(2));
    axis(ax,[xran yran]);
    longticks(ax,2);
  
  elseif strcmp(clrcode, 'rcc')
    %%%plot
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 8;   % maximum height allowed is 11 inches
    nrow = 2;
    ncol = 2;
    f = initfig(widin,htin,nrow,ncol); %initialize fig
    xran = [0.06 0.96]; yran = [0.06 0.98];
    xsep = 0.05; ysep = 0.07;
    optaxpos(f,nrow,ncol,xran,yran,xsep,ysep);

    %LFE locations, 25-win detection
    ax=f.ax(1);
    hold(ax,'on');
    ax.Box='on'; grid(ax,'on');
    xran = [-5 5];
    yran = [-5 5];
    plot(ax,[xran(1) xran(2)],[0 0],'--','Color',[.3 .3 .3]);
    plot(ax,[0 0],[yran(1) yran(2)],'--','Color',[.3 .3 .3]);
    scatter(ax,lfeloc(:,1),lfeloc(:,2),15,rcc,'filled',...
      'MarkerEdgeColor',[.5 .5 .5]);
    text(ax,0.95,0.85,sprintf('%d',size(lfeloc,1)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    oldc = colormap(ax,'kelicol');
    newc = flipud(oldc);
    colormap(ax,newc);
    c=colorbar(ax);
    clabel = 'RCC at average src arrival';
    cran = [0 1];
  %   c.Label.String = clabel;
    caxis(ax,cran);
    axis(ax,'equal');
    text(ax,0.95,0.05,sprintf('1st prctile: %.2f',prctile(rcc,1)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.1,sprintf('5th prctile: %.2f',prctile(rcc,5)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.15,sprintf('50th prctile: %.2f',prctile(rcc,50)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.2,sprintf('95th prctile: %.2f',prctile(rcc,95)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.95,'LFE, 25-s-win','Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    xlabel(ax,'E (km)');
    ylabel(ax,'N (km)');
    xticks(ax,xran(1):1:xran(2));
    yticks(ax,yran(1):1:yran(2));
    axis(ax,[xran yran]);
    longticks(ax,2);

    %LFE locations, whole-win detection
    ax=f.ax(2);
    hold(ax,'on');
    ax.Box='on'; grid(ax,'on');
    xran = [-5 5];
    yran = [-5 5];
    plot(ax,[xran(1) xran(2)],[0 0],'--','Color',[.3 .3 .3]);
    plot(ax,[0 0],[yran(1) yran(2)],'--','Color',[.3 .3 .3]);
    scatter(ax,lfeloc2(:,1),lfeloc2(:,2),15,rcc2,'filled',...
      'MarkerEdgeColor',[.5 .5 .5]);
    text(ax,0.95,0.85,sprintf('%d',size(lfeloc2,1)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    oldc = colormap(ax,'kelicol');
    newc = flipud(oldc);
    colormap(ax,newc);
    c=colorbar(ax);
    c.Label.String = clabel;
    caxis(ax,cran);
    axis(ax,'equal');
    text(ax,0.95,0.05,sprintf('1st prctile: %.2f',prctile(rcc2,1)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.1,sprintf('5th prctile: %.2f',prctile(rcc2,5)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.15,sprintf('50th prctile: %.2f',prctile(rcc2,50)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.2,sprintf('95th prctile: %.2f',prctile(rcc2,95)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.95,'LFE, whole-win','Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    xlabel(ax,'E (km)');
    ylabel(ax,'N (km)');
    xticks(ax,xran(1):1:xran(2));
    yticks(ax,yran(1):1:yran(2));
    axis(ax,[xran yran]);
    longticks(ax,2);

    %LFE locations, 25-win detection
    ax=f.ax(3);
    hold(ax,'on');
    ax.Box='on'; grid(ax,'on');
    xran = [-5 5];
    yran = [-5 5];
    plot(ax,[xran(1) xran(2)],[0 0],'--','Color',[.3 .3 .3]);
    plot(ax,[0 0],[yran(1) yran(2)],'--','Color',[.3 .3 .3]);
    scatter(ax,lfelocn(:,1),lfelocn(:,2),15,rccn,'filled',...
      'MarkerEdgeColor',[.5 .5 .5]);
    text(ax,0.95,0.85,sprintf('%d',size(lfelocn,1)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    oldc = colormap(ax,'kelicol');
    newc = flipud(oldc);
    colormap(ax,newc);
    c=colorbar(ax);
    clabel = 'RCC at average src arrival';
    cran = [0 1];
  %   c.Label.String = clabel;
    caxis(ax,cran);
    axis(ax,'equal');
    text(ax,0.95,0.05,sprintf('1st prctile: %.2f',prctile(rccn,1)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.1,sprintf('5th prctile: %.2f',prctile(rccn,5)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.15,sprintf('50th prctile: %.2f',prctile(rccn,50)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.20,sprintf('95th prctile: %.2f',prctile(rccn,95)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.95,'LFE, 25-s-win, noise','Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    xlabel(ax,'E (km)');
    ylabel(ax,'N (km)');
    xticks(ax,xran(1):1:xran(2));
    yticks(ax,yran(1):1:yran(2));
    axis(ax,[xran yran]);
    longticks(ax,2);

    %LFE synthetic noise, whole-win
    ax=f.ax(4);
    hold(ax,'on');
    ax.Box='on'; grid(ax,'on');
    xran = [-5 5];
    yran = [-5 5];
    plot(ax,[xran(1) xran(2)],[0 0],'--','Color',[.3 .3 .3]);
    plot(ax,[0 0],[yran(1) yran(2)],'--','Color',[.3 .3 .3]);
    scatter(ax,lfeloc2n(:,1),lfeloc2n(:,2),15,rcc2n,'filled',...
      'MarkerEdgeColor',[.5 .5 .5]);
    text(ax,0.95,0.85,sprintf('%d',size(lfeloc2n,1)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    oldc = colormap(ax,'kelicol');
    newc = flipud(oldc);
    colormap(ax,newc);
    c=colorbar(ax);
    c.Label.String = clabel;
    caxis(ax,cran);
    axis(ax,'equal');
    text(ax,0.95,0.05,sprintf('1st prctile: %.2f',prctile(rcc2n,1)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.1,sprintf('5th prctile: %.2f',prctile(rcc2n,5)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.15,sprintf('50th prctile: %.2f',prctile(rcc2n,50)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.2,sprintf('95th prctile: %.2f',prctile(rcc2n,95)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.95,0.95,'LFE, whole-win, noise','Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    xlabel(ax,'E (km)');
    ylabel(ax,'N (km)');
    xticks(ax,xran(1):1:xran(2));
    yticks(ax,yran(1):1:yran(2));
    axis(ax,[xran yran]);
    longticks(ax,2);
    
  end

  h=supertit([f.ax(1) f.ax(2)], sprintf('Burst #%s, %d, %.1f-%.1f s, sps: %d',...
    num2zeropadstr(idxbst(iii), 3),date,tstbuf,tedbuf,sps),10);
  movev(h,-0.3);
  
  %save figure
  fignm{iii,1} = sprintf('lfevstremor_bst%s.pdf',...
    num2zeropadstr(idxbst(iii),3));
  print(f.fig,'-dpdf',fullfile(rstpath, '/FIGS',fignm{iii,1}));
%   close(f.fig);

end

%% merge all figures into a single pdf file
if strcmp(clrcode, 'tarvl') 
  status = system('rm -f /home/data2/chaosong/matlab/allan/data-no-resp/PGCtrio/FIGS/lfevstremor_tarvl.pdf');
  for i = 1:size(fignm,1)
    mergenm{i} = fullfile(rstpath, '/FIGS/',fignm{i});
  end
  append_pdfs(fullfile(rstpath, '/FIGS/lfevstremor_tarvl.pdf'),mergenm);
  
elseif strcmp(clrcode, 'rcc')  
  status = system('rm -f /home/data2/chaosong/matlab/allan/data-no-resp/PGCtrio/FIGS/lfevstremor_rcc.pdf');
  for i = 1:size(fignm,1)
    mergenm{i} = fullfile(rstpath, '/FIGS/',fignm{i});
  end
  append_pdfs(fullfile(rstpath, '/FIGS/lfevstremor_rcc.pdf'),mergenm);
  
end
