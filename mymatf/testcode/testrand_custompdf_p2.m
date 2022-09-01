% testrand_custompdf_p2.m
%
% This script serves as a test and reminder of obtaining the custom PDF/CDF from a dataset, not
% necessarily a Gaussian (normal) distribution. Given the PDF/CDF, we try to explore how to draw 
% pseudo-random numbers from it.
% 
% Part 2, a custom PDF that is not known beforehand, can not be expressed as a function of x1 and
% x2, so you need to estimate a discrete PDF from x1 and x2 first, then apply 'pinky'.
%
% --In Part 1, for 2-D case, 'pinky(x1,x2,y)' works fine to generate random numbers from any custom
%   PDF and vectors x1 and x2, if PDF is known as an analytic function or discrete values.
%   If you don't know about the PDF of a discrete
%   data set at all, ie. you only have x1 and x2 vectors (but from them you have a feeling of 
%   counts), you have a few options. 
%   Option 1 is to bin the x1 and x2 upon pixel (ie. find unique entries) or to bin upon grid point
%   with a size. This is similar to getting the density distribution from x1 and x2. Normalizing the
%   density (counts), you can obtain PDF. Then you can feed this PDF to 'pinky' to draw samples from
%   it.
%   Similarly to option 1, option 2 is to use 'histogram2' to bin for you, specify the binedges
%   carefully, and normalize it to PDF, then you can get the PDF as the height and histogram. From
%   plotting, you could see that they are almost identical.
%   Option 3 is '[f,xi] = ksdensity(x,pts)' where x is the 2-colomn vector [x1 x2], pts is the
%   location points that you want to estimate the PDF (output f) on. Then you can feed this discrete 
%   PDF to 'pinky' to draw samples from it. (help ksdensity and see the last example). But note that
%   the 'function' can NOT deal with 'icdf', meaning that if can NOT draw 2D random numbers using
%   this funtion directly. There is a big advantage of 'ksdensity' compared to binning method in
%   options 1 and 2. A histogram represents the probability distribution by establishing bins and
%   placing each data value in the appropriate bin. Because of this bin count approach, the 
%   histogram produces a discrete probability density function. This might be unsuitable for certain
%   applications, such as generating random numbers from a fitted distribution.
%   Alternatively, the kernel distribution builds the probability density function (pdf) by creating
%   an individual probability density curve for each data value, then summing the smooth curves. 
%   This approach creates one smooth, continuous probability density function for the data set.
% --All of them has its own advantage, depending on your needs, but generally I kinda like 
%   'ksdensity' more due to its smoothing nature.  
%   
% 
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/04/06
% Last modified date:   2022/04/07

clear;
clc;
close all

%% Construct a custom 2D PDF using the true density distribution of 4-s tremor detections
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

%%%Interpolate from existing grids for the locations of sources
%round to the integer sample at the 'iup' times the original sampling rate 
offint = round(hftime(:, 1:2)*iup);
          
%convert the time offset to locations
ftrans = 'interpArmb';
fplt = 0;
ind = find(offint(:,1)>=-23*iup & offint(:,1)<=25*iup & ...
  offint(:,2)>=-24*iup & offint(:,2)<=24*iup);
[loc, indinput] = off2space002(offint(ind,:),sps*iup,ftrans,fplt);
% loc has 8 cols, format: dx,dy,lon,lat,dep,tori,off12,off13
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

%% obtain the density then PDF by binning 
%Obtain the density (count), bin by each desired grid point
dx = 0.05;
dy = 0.05;
xran = [-4 4];
yran = [-4 4];

[density,xgrid,ygrid,densitygrid]= density_matrix(hfbnd(:,1),hfbnd(:,2),xran,yran,dx,dy);
emppdf = density;
emppdf(:,3) = emppdf(:,3)/sum(emppdf(:,3))/(dx*dy); % normalize to PDF, ==counts/N_total/area_of_bin
emppdfgrid = densitygrid/sum(sum(densitygrid))/(dx*dy);

%visualize the empirical PDF
figure;
hold on
% scatter(emppdf(:,1),emppdf(:,2),25,emppdf(:,3),'s','filled'); % use scatter if the density is 1D
imagesc(xgrid(1,:)',ygrid(:,1),emppdfgrid); hold on; % use 'imagesc' if the density is 2D matrix
colormap(jet)
c=colorbar;
c.Label.String = 'Probability Density';
caxis([0 max(emppdf(:,3))]);
axis equal tight
xlim(xran);
ylim(yran);
xlabel('x1');
ylabel('x2');
box on; grid on;
title('Data PDF by binning upon grids');

%%% An alternative to visualize and/or obtain the PDF is to use 'histogram2', theretically gives 
%%% the same result, 'h.Values' gives the PDF value
figure;
hold on
h = histogram2(hfbnd(:,1),hfbnd(:,2),'BinWidth',[dx dy],'FaceColor','flat',...
  'Normalization','pdf');
colormap(jet);
c=colorbar;
c.Label.String = 'Probability Density';
% view(3);
h.DisplayStyle = 'tile';
h.ShowEmptyBins = 'on';
axis equal
xlim(xran);
ylim(yran);
xlabel('x1');
ylabel('x2');
box on; grid on;
title('Data PDF by 2D histogram');

%% Draw random numbers from the designed 2D normal PDF with 'pinky'
rng('default')  % For reproducibility
N=10000;
rdnum=zeros(N,2);
for i=1:N
  [rdnum(i,1),rdnum(i,2)]=pinky(xgrid(1,:)',ygrid(:,1),emppdfgrid);
end

%visualize the empirical PDF
figure;
hold on
% scatter(emppdf(:,1),emppdf(:,2),25,emppdf(:,3),'s','filled'); % use scatter if the density is 1D
imagesc(xgrid(1,:)',ygrid(:,1),emppdfgrid); hold on; % use 'imagesc' if the density is 2D matrix
colormap(jet)
c=colorbar;
c.Label.String = 'Probability Density';
caxis([0 max(emppdf(:,3))]);
scatter(rdnum(:,1),rdnum(:,2),10,[.5 .5 .5],'+');
axis equal tight
xlim(xran);
ylim(yran);
xlabel('x1');
ylabel('x2');
box on; grid on;
title('Data PDF by binning upon grids with drawn random samples');

%% Justify if drawn numbers comply to PDF? 1. plot the density map to see if similar to the data PDF
%%% 1. plot the density map of drawn numbers to see if it is similar to the data PDF
%Obtain the density (count), bin by each pixel (intrinsic resolution)
% density = density_pixel(R(:,1), R(:,2));

%Obtain the density (count), bin by each desired grid point
[densrd,xgrid,ygrid,densgridrd] = density_matrix(rdnum(:,1),rdnum(:,2),xran,yran,dx,dy);
emppdfrd = densrd;
emppdfrd(:,3) = emppdfrd(:,3)/sum(emppdfrd(:,3))/(dx*dy); % normalize to PDF, ==counts/N_total/area_of_bin
emppdfgridrd = densgridrd/sum(sum(densgridrd))/(dx*dy);

figure;
hold on
% scatter(emppdfrd(:,1),emppdfrd(:,2),25,emppdfrd(:,3),'s','filled'); % use scatter if the density is 1D
imagesc(xgrid(1,:)',ygrid(:,1),emppdfgridrd); hold on; % use 'imagesc' if the density is 2D matrix
colormap(jet)
c=colorbar;
c.Label.String = 'Probability Density';
caxis([0 max(emppdfrd(:,3))]);
axis equal tight
xlim(xran);
ylim(yran);
xlabel('x1');
ylabel('x2');
box on
title('Random sample PDF by binning upon grids');

%%% An alternative to visualize the PDF is to use 'histogram2', theretically gives the same result 
figure;
hold on
h = histogram2(rdnum(:,1), rdnum(:,2),'BinWidth',[dx dy],'FaceColor','flat','Normalization','pdf');
colormap(jet);
c=colorbar;
c.Label.String = 'Probability Density';
% view(3);
h.DisplayStyle = 'tile';
h.ShowEmptyBins = 'on';
axis equal
xlim(xran);
ylim(yran);
xlabel('x1');
ylabel('x2');
box on; grid on;
title('Random sample PDF by 2D histogram');


%% Alternatively, use 'ksdensity' to get its smoothed PDF
%Create a two-column vector of points at which to evaluate the density.
xran = [-4 4];
yran = [-4 4];
dx = 0.05;
dy = 0.05;
x = xran(1)+dx/2: dx: xran(2)-dx/2;
y = yran(2)-dy/2: -dy: yran(1)+dy/2;
[xgrid,ygrid] = meshgrid(x, y);
ksemppdf(:,1:2) = [reshape(xgrid,[],1) reshape(ygrid,[],1)];

ksemppdf(:,3) = ksdensity(hfbnd(:,1:2),ksemppdf(:,1:2),'Function','pdf');
ksemppdfgrid = reshape(ksemppdf(:,3), length(y), length(x));

%visualize the empirical PDF
figure;
hold on
% scatter(ksemppdf(:,1),ksemppdf(:,2),25,ksemppdf(:,3),'s','filled'); % use scatter if the density is 1D
imagesc(xgrid(1,:)',ygrid(:,1),ksemppdfgrid); hold on; % use 'imagesc' if the density is 2D matrix
colormap(jet)
c=colorbar;
c.Label.String = 'Probability Density';
caxis([0 max(ksemppdf(:,3))]);
axis equal tight
xlim(xran);
ylim(yran);
xlabel('x1');
ylabel('x2');
box on; grid on;
title('Data PDF by ksdensity');

%% Draw random numbers from the designed 2D normal PDF with 'pinky'
ksrdnum=zeros(N,2);
for i=1:N
  [ksrdnum(i,1),ksrdnum(i,2)]=pinky(xgrid(1,:)',ygrid(:,1),ksemppdfgrid);
end
%visualize the empirical PDF
figure;
hold on
% scatter(ksemppdf(:,1),ksemppdf(:,2),25,ksemppdf(:,3),'s','filled'); % use scatter if the density is 1D
imagesc(xgrid(1,:)',ygrid(:,1),ksemppdfgrid); hold on; % use 'imagesc' if the density is 2D matrix
colormap(jet)
c=colorbar;
c.Label.String = 'Probability Density';
caxis([0 max(ksemppdf(:,3))]);
scatter(ksrdnum(:,1),ksrdnum(:,2),10,[.5 .5 .5],'+');
axis equal tight
xlim(xran);
ylim(yran);
xlabel('x1');
ylabel('x2');
box on; grid on;
title('Data PDF by ksdensity with drawn random samples');


%% Justify if drawn numbers comply to PDF? 1. plot the density map to see if similar to the data PDF
ksemppdfrd(:,1:2) = [reshape(xgrid,[],1) reshape(ygrid,[],1)];

ksemppdfrd(:,3) = ksdensity(ksrdnum(:,1:2),ksemppdfrd(:,1:2),'Function','pdf');
ksemppdfrdgrid = reshape(ksemppdf(:,3), length(y), length(x));

%visualize the empirical PDF
figure;
hold on
% scatter(ksemppdfrd(:,1),ksemppdfrd(:,2),25,ksemppdfrd(:,3),'s','filled'); % use scatter if the density is 1D
imagesc(xgrid(1,:)',ygrid(:,1),ksemppdfrdgrid); hold on; % use 'imagesc' if the density is 2D matrix
colormap(jet)
c=colorbar;
c.Label.String = 'Probability Density';
caxis([0 max(ksemppdf(:,3))]);
axis equal tight
xlim(xran);
ylim(yran);
xlabel('x1');
ylabel('x2');
box on; grid on;
title('Random sample PDF by ksdensity');













