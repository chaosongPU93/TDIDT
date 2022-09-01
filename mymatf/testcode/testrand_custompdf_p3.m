% testrand_custompdf_p3.m
%
% This script serves as a test and reminder of obtaining the custom PDF/CDF from a dataset, not
% necessarily a Gaussian (normal) distribution. Given the PDF/CDF, we try to explore how to draw 
% pseudo-random numbers from it.
% 
% Part 3, this code basically is cut and pasted here from 'synthshift_chao.m' where I started there
% to test if it is doable.
%
% The purpose of this code is to try if it is possible to switch getting the PDF in location dx and
% dy to PDF in offset sample space. We try to avoid doing the transformation between locations and
% offset in generating the synthetics, as for each source, the source location provides the time
% shift of the arrival at station 2 or 3 relative to station 1. If you want to bin the detections
% inside the region of interest by a regular grid, it is easier to normalize for the PDF; but if you
% want to bin by pixel (ie., each unique source), then it is hard to get the 'size' of each bin.
% However, if you view the sources in terms of offset, they are integer samples at a specific
% sampling rate which is easier for obtaining the size (1*1 sample). Therefore, you may outline 
% a boundary in terms of dx and dy then get the corresponding off12 and off13, then get the PDF in 
% terms of offset, then draw random samples from that PDF. Alternatively, you may start with 
% outlining a boundary directly in terms of time offset. The rest is just to test if the random
% samples drawn from the data PDF indeed follow the same distribution.
% 
% 
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/04/18
% Last modified date:   2022/04/18

clear;
clc;
close all;

rstpath = strcat(datapath, '/PGCtrio');
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

%%%Interpolate from existing grids for the locations of sources
%round to the integer sample at the 'iup' times the original sampling rate 40 sps
offint = round(hftime(:, 1:2)*iup);
          
%convert the time offset to locations
% ftrans = 'interpArmb';
ftrans = 'interpArmbreloc';
% ftrans = 'interpchao';
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

%Obtain the empirical PDF binned by pixel, via normalization according to definition 
epdfloc = denloc;
%Here we don't know the exact spacing of the each 'bin' as it is not even, we can only roughly
%estimate it
%Get the pairwise Euclidean distance between all pairs of points
dist = pdist(denloc(:,1:2));
binwid = min(dist);  % choose the min as the length of a square grid
clear dist;  %release the memory
%Assume the grid is square-ish
epdfloc(:,3) = epdfloc(:,3)/sum(epdfloc(:,3))/binwid^2; % normalize to PDF, ==counts/N_total/area_of_bin

%plot the PDF, note the color should be the SAME as density, as it is just a normalization
figure
scatter(epdfloc(:,1),epdfloc(:,2),6,epdfloc(:,3),'filled');
axis equal
axis([xran yran]);
colormap(jet)
c=colorbar;
c.Label.String = 'Probability Density';
caxis([0 max(epdfloc(:,3))]);
xlabel('E (km)');
ylabel('N (km)');
box on; grid on;


%%%Note that the locations in the ellipse corresponds to integer time offsets at sps*iup, so we can
%%%avoid doing the transformation, and directly get the density and PDF in sample space, the benefit
%%%is that, in sample sapce, it is actucally an even grid with a spacing of 1 sample 
denoff= density_pixel(hfbnd(:,7),hfbnd(:,8));

%plot the density distribution
figure
scatter(denoff(:,1),denoff(:,2),30,denoff(:,3),'filled');
axis equal
off12ran = [min(denoff(:,1))-1 max(denoff(:,1))+1];
off13ran = [min(denoff(:,2))-1 max(denoff(:,2))+1];
axis([off12ran off13ran]);
colormap(jet)
c=colorbar;
c.Label.String = '# tremor detections / pixel';
caxis([0 max(denoff(:,3))]);
xlabel(strcat({'Offset 12 (samples at '},num2str(sps),' sps)'));
ylabel(strcat({'Offset 13 (samples at '},num2str(sps),' sps)'));
box on; grid on;

%Obtain the PDF
epdfoff = denoff;
%now the 'area' of each bin is 1*1 sample
epdfoff(:,3) = epdfoff(:,3)/sum(epdfoff(:,3))/(1*1); % normalize to PDF, ==counts/N_total/area_of_bin

%plot the PDF, note the color should be the SAME as density, as it is just a normalization
figure
scatter(epdfoff(:,1),epdfoff(:,2),30,epdfoff(:,3),'filled');
axis equal
axis([off12ran off13ran]);
colormap(jet)
c=colorbar;
c.Label.String = 'Probability Density';
caxis([0 max(epdfoff(:,3))]);
xlabel(strcat({'Offset 12 (samples at '},num2str(sps),' sps)'));
ylabel(strcat({'Offset 13 (samples at '},num2str(sps),' sps)'));
box on; grid on;

%% test of drawing random samples from the above PDFs
%%% Note that 'pinky(x1,x2,y)' can only deal with an evenly-spaced vector x1 and x2, y is a grid of
%%% PDF values on the grid points that are defined by x1 and x2.
%%% Therefore, 'pinky(x1,x2,y)' can be applied when either your PDF is generated by binning upon an 
%%% evenly-spaced grid (or ksdensity), OR, you bin by pixel in locations but return to the sample
%%% domain, in which case the grid is even. But the latter case needs zero-padding to grid points
%%% that is otherwise empty.
xvec = -max(abs(epdfoff(:,1)))-1: 1: max(abs(epdfoff(:,1)))+1;
yvec = -max(abs(epdfoff(:,2)))-1: 1: max(abs(epdfoff(:,2)))+1;
[epdfoffpad,xgrid,ygrid,epdfoffgrid,ind1] = zeropadmat2d(epdfoff,xvec,yvec);

%plot the PDF, note the color should be the SAME as density, as it is just a normalization
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
box on; grid on;

%plot the gridded data just to confirm the two are equivalent
figure
imagesc(xgrid(1,:)',ygrid(:,1),epdfoffgrid); hold on; % use 'imagesc' if the density is 2D matrix
dum = epdfoffpad;
dum(dum(:,3)==0, :) = [];
scatter(dum(:,1),dum(:,2),20,dum(:,3),'filled','MarkerEdgeColor','w');
axis equal
axis([minmax(xvec) minmax(yvec)]);
colormap(jet)
c=colorbar;
c.Label.String = 'Probability Density';
caxis([0 max(epdfoffpad(:,3))]);
xlabel(strcat({'Offset 12 (samples at '},num2str(sps),' sps)'));
ylabel(strcat({'Offset 13 (samples at '},num2str(sps),' sps)'));
box on; grid on;
ax=gca;
ax.YDir='normal'; %reverse

%%%Draw random samples from the pdf
rng('default')  % For reproducibility
N=5000;
rdoff=zeros(N,2);
for i=1:N
  [rdoff(i,1),rdoff(i,2)]=pinky(xgrid(1,:)',ygrid(:,1),epdfoffgrid);
end
denrdoff= density_pixel(rdoff(:,1),rdoff(:,2));

%plot the density distribution
figure
scatter(denrdoff(:,1),denrdoff(:,2),30,denrdoff(:,3),'filled');
axis equal
axis([minmax(xvec) minmax(yvec)]);
colormap(jet)
c=colorbar;
c.Label.String = '# random samples / pixel';
caxis([0 max(denrdoff(:,3))]);
xlabel(strcat({'Offset 12 (samples at '},num2str(sps),' sps)'));
ylabel(strcat({'Offset 13 (samples at '},num2str(sps),' sps)'));
box on; grid on;










