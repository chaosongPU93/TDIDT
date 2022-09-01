% locinterp002_4s.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The next script to run after 'mergeday002_4s.m', is to interpolate for the
% locations of detections on the slab according to a pre-computed grid of 
% locations from an offset grid with a 1-sample spacing ranges +-13 samples
% at 20 sps. As the mshift in 'detection002_4s.m' is the same as the max range
% of the existed grid, interpolation is possible.
% After interpolation, plot density map, choose the 002 region with high 
% density, 1.75*1 km ellipse / cicle / rectangle, output detections inside to
% file.
%
% --I shall reiterate the order of the entire flow of detection here: 
%   'detection002_4s.m' --> 'mergeday002_4s.m' --> 'locinterp002_4s.m'
%   --> inter-event time, get the windows of bursts for future deconvolution
% --Here, we avoid the direct inversion using Hypoinverse, as the time
%   resolution demanded here exceeds the effective max spatial resolution
%   of hypoinverse. You can surely invert, but the location scattering would 
%   not be regular, leading to unpredictable distortion
% --Use 'off2space002' to interplate
% --I think it is fine to plot the density map in this script, but I'd rather
%   put the rest into the next script.
% --In terms of the density map, for the purpose of getting true detection 
%   distribution and construct a cutoff region, the better way of binning is
%   to base upon pixels, rather than a rectangular grid.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/03/10
% Last modified date:   2022/03/14
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

[scrsz, res] = pixelperinch(1);

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band 
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
%%% choose the data path here 
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/PGCtrio');

freqflag='hf';  % flag to indicate whether to do hf or lf;

FLAG = 'PGC'; % detector
  
fam = '002';   % family number

sps = 40;

iup = 4;  % upsample 4 times

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
winlen=winlensec*sps;      % length in smaples
PREFIX = strcat(fam,'.up.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),'.',...
                        int2str(npo),int2str(npa),'.ms', int2str(mshift));
hftime = load(strcat(rstpath, '/MAPS/tdectimeori_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',...
            num2str(winlen/sps),'s',num2str(sps),'sps','4add'));
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
  [loc, indinput] = off2space002(offint,sps*iup,ftrans,fplt);
  % loc has 8 cols, format: dx,dy,lon,lat,dep,tori,off12,off13
  hfuse = hftime;
  hfuse = [loc hfuse];
elseif isequal(ftrans,'interpArmb') || isequal(ftrans,'interpArmbreloc')
  ind = find(offint(:,1)>=-23*iup & offint(:,1)<=25*iup & ...
             offint(:,2)>=-24*iup & offint(:,2)<=24*iup);
  [loc, indinput] = off2space002(offint(ind,:),sps*iup,ftrans,fplt);
  % loc has 8 cols, format: dx,dy,lon,lat,dep,tori,off12,off13
  hfuse = hftime(ind,:);
  hfuse = [loc hfuse];
end

%%
% avemed = average_pixel(hfuse(:,1),hfuse(:,2),[hfuse(:,9:10) hfuse(:,10)-hfuse(:,9)],'median');
avemean = average_pixel(hfuse(:,1),hfuse(:,2),[hfuse(:,9:10) hfuse(:,10)-hfuse(:,9)],'mean');

%get the 1-D derivative of offset wrt. relative location
for i = 1: 3
  dzdx(:,i) = dfxdx(avemean(:,2+i), avemean(:,1));
  dzdy(:,i) = dfxdx(avemean(:,2+i), avemean(:,2));
end
%vector orientation at each point
% vecang = atan2d(dzdy,dzdx);
% vecang = -vecang; 
% vecang(vecang<0) = vecang(vecang<0)+360;
% vecang = vecang+90;
% vecang(vecang>360) = vecang(vecang>360)-90;

vecang = atan2d(dzdy,dzdx);
mvecang = median(vecang);

%%
%simply scatter the detections to check the spatial resolution
xran = [-8 3];
yran = [-4 4];
figure
for i = 1: 3
  subplot(3,1,i)
  scatter(avemean(:,1),avemean(:,2),10,avemean(:,2+i),'filled');
  colormap('jet');
  axis equal
%   axis([xran yran]);
  [rotx, roty] = complex_rot(0,1,-mvecang(i));
  xvect = [0-rotx 0+rotx];
  yvect = [0-roty 0+roty];
  drawArrow(gca,xvect,yvect,xran,yran,'linewidth',1);
  c=colorbar;
  caxis([-15 15]);
  box on
  if i==1
    title('PGC-SSIB (off12)');
  elseif i==2
    title('PGC-SILB (off13)');
  else
    title('SSIB-SILB (off23)');
  end
  longticks(gca,1);
  if i == 3
    xlabel('E (km)');
    ylabel('N (km)');
    c.Label.String = sprintf('mean offset (samples) at %d Hz',sps);
  end

end

%%
% avemed = average_pixel(hfuse(:,9),hfuse(:,10),[hfuse(:,9:10) hfuse(:,10)-hfuse(:,9)],'median');
avemean = average_pixel(hfuse(:,9),hfuse(:,10),[hfuse(:,9:10) hfuse(:,10)-hfuse(:,9)],'mean');

%simply scatter the detections to check the spatial resolution
xran = [-30 30];
yran = [-20 20];
figure
for i = 1: 3
  subplot(3,1,i)
  scatter(avemean(:,1),avemean(:,2),10,avemean(:,2+i),'filled');
  colormap('jet');
  axis equal
  axis([xran yran]);
  c=colorbar;
  caxis([-15 15]);
  box on
  if i==1
    title('PGC-SSIB (off12)');
  elseif i==2
    title('PGC-SILB (off13)');
  else
    title('SSIB-SILB (off23)');
  end
  longticks(gca,1);
  if i == 3
    xlabel('PGC-SSIB (off12) (samples)');
    ylabel('PGC-SILB (off13) (samples)');
    c.Label.String = sprintf('mean offset (samples) at %d Hz',sps);
  end

end

%%
xran = [-8 3];
yran = [-4 4];

%simply scatter the detections to check the spatial resolution
figure
scatter(loc(:,1),loc(:,2),4,'ko','filled');
colormap('jet');
axis equal
axis([xran yran]);
box on

%plot the cumulative density map, binning by a 0.05*0.05 km grid, looks similar to that is binned 
%upon pixel
[f1] = plt_cumulative_density(loc,[],[-20 20],[-20 20],'grid',6,0.05,0.05);          
f1.ax.XLim = xran;
f1.ax.YLim = yran;
title(f1.ax, 'Density binned by grid');

%% emphasize the density map by binning upon pixel, then outline the boundary
%plot the cumulative density map, binning by pixel, ie., each unique detection         
[f2] = plt_cumulative_density(loc,[],[-20 20],[-20 20],'pixel');          
f2.ax.XLim = xran;
f2.ax.YLim = yran;
title(f2.ax, 'Density binned by pixel');
caxis(f2.ax,[0 1.4]);

%outline a boundary region on top of the density map
cutout = 'ellipse';

if isequal(cutout,'circle')
  x0 = 0.2; 
  y0 = 0.2;
  radi = 1.75;
  [x, y] = circle_chao(x0,y0,radi,0.1);
elseif isequal(cutout,'ellipse')
  x0 = 0.2; 
  y0 = 0.2;
  semia = 1.75;
  semib = 1.25;
  angrot = 45;
%   [x, y] = ellipse_chao(x0,y0,semia,semib,0.01);
  [x, y] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);
elseif isequal(cutout,'rectangle')
  EW = [-7 3];
  NS = [-3 4];
  wid = range(EW);
  hgt = range(NS);
  x0 = mean(EW); 
  y0 = mean(NS);  
  [x, y] = rectangle_chao(x0,y0,wid,hgt,0.01);
%   [x, y] = rectangle_chao(x0,y0,wid,hgt,0.01,angrot,[x0,y0]);
end        
          
hold(f2.ax,'on');
plot(f2.ax,x,y,'k-','linew',2);
scatter(f2.ax,x0,y0,10,'ko','linew',1);
hold(f2.ax,'off');


%% obtain detections that fall in the cut boundary
bnd = [x y];
[iin,ion] = inpolygon(hfuse(:,1),hfuse(:,2),bnd(:,1),bnd(:,2));
isinbnd = iin | ion;
hfbnd = hfuse(isinbnd == 1, :);

%plot the cumulative density map, binning by pixel, ie., each unique detection         
[f3] = plt_cumulative_density(hfbnd,[],[-20 20],[-20 20],'pixel');          
f3.ax.XLim = xran;
f3.ax.YLim = yran;
caxis(f3.ax,[0 1.4]);

figure
scatter(hfbnd(:,1),hfbnd(:,2),4,'ko','filled');
axis equal
axis([xran yran]);
box on

keyboard

%% write into files
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

fid = fopen(strcat(rstpath, '/MAPS/eloc.pgc.tdec_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',...
            num2str(winlen/sps),'s',num2str(sps),'sps4add_',cutout(1:4)),'w+');
fprintf(fid,'%.4f %.4f %.4f %.4f %.4f %.4f %d %d %6.2f %6.2f %.4f %.4f %d %d %9.1f %10.3f %8.3f %7.2f %10.3e %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f \n',...
        hfbnd');
fclose(fid);

fid = fopen(strcat(rstpath, '/MAPS/eloc.pgc.tdec_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',...
            num2str(winlen/sps),'s',num2str(sps),'sps4add_'),'w+');
fprintf(fid,'%.4f %.4f %.4f %.4f %.4f %.4f %d %d %6.2f %6.2f %.4f %.4f %d %d %9.1f %10.3f %8.3f %7.2f %10.3e %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f \n',...
        hfuse');
fclose(fid);










