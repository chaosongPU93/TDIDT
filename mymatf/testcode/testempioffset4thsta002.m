% testempioffset4thsta002.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script aims to obtain an empirical estimate of the time offset needed
% for the 4th station among LZB, TWKB, MGCB, KLNB to align with the detected
% sources. The estimate comes from the 4th stational check result to 4-s 
% tremor detections in 'detection002_4s.m'. The estimate will be used to
% determine the offset needed for a deconvolved LFE source based on its
% location, which is the output in codes like 'deconvbursts002_ref_4s_exp.m'.
% 
% --4-s tremor detections checked by diff 4th stations are diff. The same 
%   location have duplicates. Ideally, the observational (or theroratical)
%   predication of the time offset between sta 4 and 1 is smoothly varying
%   wrt the source location, which is also your eyes seeing. However, there
%   are many 'spikes' in the data sets on top of a smooth surface.
% --to smooth out these spikes, this script tries a few ways, like 2D 
%   interpolation, gradient estimate, and 2D plane fitting. The test shows
%   that interpolation does not work at all; the gradient might not be very
%   stable; plane fitting seems like the best solution so far in terms of 
%   having a 1st-order estimate for the time offset.
%
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/10/04
% Last modified date:   2022/10/04
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
  'SILB '];     % determine the trio and order, here the 1st sta is PGC
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


%%
sps = 160;

stasnew=['LZB  '
         'TWKB '
         'MGCB '
         'KLNB '];  % twkb lzb mgcb
nstasnew = size(stasnew,1); 

f1 = initfig(16,9,3,4);
dataplt = hfall;

flagsmooth = 0;

if isequaln(dataplt, hfall)
  sybsz = 5; 
  filtsigma = 5;
elseif isequaln(dataplt, hfbnd)
  sybsz = 20;
  filtsigma = 1;
end

xran = [min(dataplt(:,7))-1 max(dataplt(:,7))+1];
yran = [min(dataplt(:,8))-1 max(dataplt(:,8))+1];
xvec = xran(1): 1: xran(2);
yvec = yran(1): 1: yran(2);

for i = 1: nstasnew
  colflg = 42+i;  % column num for the flag to indicate if the 4th sta passes the check
  colloff = 46+i; % column num for the summed difference in offset needed to align 4th sta and original trio
  colioff = 50+i; % column num for the offset needed to align 4th sta and original trio
  colcc = 54+i; % column num for ave CC of aligning 4th sta and original trio
  ccavethres = 0.8*xcmaxAVEnmin;  %% average cc threshold
  offavethres = 2*sps/40;    % average differential offset threshold
  
  ind = find(dataplt(:,colflg)==1);
  loff = dataplt(ind,colloff)*sps/40;
  ioff = dataplt(ind,colioff)*sps/40;
  cc = dataplt(ind,colcc);
  off12 = dataplt(ind,7);
  off13 = dataplt(ind,8);
  mloff = median_pixel(off12,off13,loff);
  mioff = median_pixel(off12,off13,ioff);
  mcc = median_pixel(off12,off13,cc);

  ax = f1.ax(i);
  hold(ax,'on'); grid(ax,'on');
  scatter(ax,dataplt(:,7),dataplt(:,8),sybsz,[.8 .8 .8],'filled');
  scatter(ax,mloff(:,1),mloff(:,2),sybsz,mloff(:,3),'filled');
  c=colorbar(ax);
  c.Label.String = 'sum of offset diff.';
  %  colormap(ax,flipud(colormap(ax,'gray')));
  colormap(ax,'jet');
  caxis(ax,[0 offavethres]);
  axis equal
  axis(ax,[xran yran]);
  hold(ax,'off');
  
  ax = f1.ax(4+i);
  hold(ax,'on'); grid(ax,'on');
  if flagsmooth
    F = scatteredInterpolant(mioff(:,1),mioff(:,2),mioff(:,3));
    F.Method = 'natural';
    F.ExtrapolationMethod = 'none';
    [xgrid, ygrid] = meshgrid(xvec,yvec);
    mioffi = [];
    mioffi(:,1) = reshape(xgrid,[],1);
    mioffi(:,2) = reshape(ygrid,[],1);
    mioffi(:,3) = F(mioffi(:,1),mioffi(:,2));
    zgrid = reshape(mioffi(:,3),size(xgrid));
    scatter(ax,mioffi(:,1),mioffi(:,2),sybsz,mioffi(:,3),'filled');
%     zgridgf = imgaussfilt(zgrid, filtsigma);
%     scatter(ax,reshape(xgrid,[],1),reshape(ygrid,[],1),sybsz,reshape(zgridgf,[],1),'filled');
%     [xyzgridpad,xgrid,ygrid,zgrid,ind2] = zeropadmat2d(mioff,xvec,yvec,prctile(mioff(:,3),5));
    
  else
    scatter(ax,mioff(:,1),mioff(:,2),sybsz,mioff(:,3),'filled');
  end
  c=colorbar(ax);
  c.Label.String = 'off14';
  colormap(ax,'jet');
  caxis(ax,[prctile(mioff(:,3),5), prctile(mioff(:,3),95)]);
  axis equal
  axis(ax,[xran yran]);
  hold(ax,'off');
  
  ax = f1.ax(8+i);
  hold(ax,'on'); grid(ax,'on');
  scatter(ax,mcc(:,1),mcc(:,2),sybsz,mcc(:,3),'filled');
  c=colorbar(ax);
  c.Label.String = 'CC';
  colormap(ax,'jet');
  caxis(ax,[0 prctile(mcc(:,3),95)]);
  axis equal
  axis(ax,[xran yran]);
  hold(ax,'off');
  
  title(f1.ax(1),stasnew(1,:));
  title(f1.ax(2),stasnew(2,:));
  title(f1.ax(3),stasnew(3,:));
  title(f1.ax(4),stasnew(4,:));
  xlabel(f1.ax(1),'off12 (160 sps)');
  ylabel(f1.ax(1),'off13 (160 sps)');
end

%% use the gradient? 
f2 = initfig(16,9,3,4);
dataplt = hfall;

for i = 4: nstasnew
  colflg = 42+i;  % column num for the flag to indicate if the 4th sta passes the check
  colloff = 46+i; % column num for the summed difference in offset needed to align 4th sta and original trio
  colioff = 50+i; % column num for the offset needed to align 4th sta and original trio
  colcc = 54+i; % column num for ave CC of aligning 4th sta and original trio
  ccavethres = 0.8*xcmaxAVEnmin;  %% average cc threshold
  offavethres = 2;    % average differential offset threshold
  
  ind = find(dataplt(:,colflg)==1);
  loff = dataplt(ind,colloff)*sps/40;
  ioff = dataplt(ind,colioff)*sps/40;
  cc = dataplt(ind,colcc);
  off12 = dataplt(ind,7);
  off13 = dataplt(ind,8);
  mloff = median_pixel(off12,off13,loff);
  mioff = median_pixel(off12,off13,ioff);
  mcc = median_pixel(off12,off13,cc);

  ran = [-20 20 -60 60];
  ind = find(mioff(:,1)>=ran(1) & mioff(:,1)<=ran(2) & mioff(:,2)>=ran(3) & mioff(:,2)<=ran(4));
  x=mioff(ind,1); y=mioff(ind,2); z=mioff(ind,3);
  
  F = scatteredInterpolant(mioff(:,1),mioff(:,2),mioff(:,3));
  F.Method = 'natural';
  F.ExtrapolationMethod = 'none';
  [xgrid, ygrid] = meshgrid(xvec,yvec);
  mioffi = [];
  mioffi(:,1) = reshape(xgrid,[],1);
  mioffi(:,2) = reshape(ygrid,[],1);
  mioffi(:,3) = F(mioffi(:,1),mioffi(:,2));
  zgrid = reshape(mioffi(:,3),size(xgrid));
  dx = 10;
  dy = 10;
  [dzdxgrid, dzdygrid] = gradient(zgrid,dx,dy);
  dzdx = reshape(dzdxgrid,[],1);
  dzdy = reshape(dzdygrid,[],1);
  dzscalar = sqrt(dzdx.^2+dzdy.^2);
  
  ran = [-20 20 -60 60];
  ind = find(mioffi(:,1)>=ran(1) & mioffi(:,1)<=ran(2) & mioffi(:,2)>=ran(1) & mioffi(:,2)<=ran(2));
  
  meddzdx = median(dzdx(ind));
  meddzdy = median(dzdy(ind));
  meddzscalar = median(dzscalar(ind));
  gradvec = atan2d(meddzdy,meddzdx);

  
  ax = f2.ax(i);
  hold(ax,'on'); grid(ax,'on');
  scatter(ax,mioffi(:,1),mioffi(:,2),sybsz,reshape(dzdx,[],1),'filled');
  c=colorbar(ax);
  c.Label.String = 'off14';
  colormap(ax,'jet');
  caxis(ax,[prctile(reshape(dzdx,[],1),5), prctile(reshape(dzdx,[],1),95)]);
  axis equal
  axis(ax,[xran yran]);
  hold(ax,'off');
  
  ax = f2.ax(i+4);
  hold(ax,'on'); grid(ax,'on');
  scatter(ax,mioffi(:,1),mioffi(:,2),sybsz,reshape(dzdy,[],1),'filled');
  c=colorbar(ax);
  c.Label.String = 'off14';
  colormap(ax,'jet');
  caxis(ax,[prctile(reshape(dzdy,[],1),5), prctile(reshape(dzdy,[],1),95)]);
  axis equal
  axis(ax,[xran yran]);
  hold(ax,'off');

  ax = f2.ax(i+8);
  hold(ax,'on'); grid(ax,'on');
  scatter(ax,mioffi(:,1),mioffi(:,2),sybsz,reshape(dzscalar,[],1),'filled');
  c=colorbar(ax);
  c.Label.String = 'off14';
  colormap(ax,'jet');
  caxis(ax,[prctile(reshape(dzscalar,[],1),5), prctile(reshape(dzscalar,[],1),95)]);
  axis equal
  axis(ax,[xran yran]);
  hold(ax,'off');

  
end

%% 2D plane fitting
dataplt = hfall;

for i = 1: nstasnew
  colflg = 42+i;  % column num for the flag to indicate if the 4th sta passes the check
  colloff = 46+i; % column num for the summed difference in offset needed to align 4th sta and original trio
  colioff = 50+i; % column num for the offset needed to align 4th sta and original trio
  colcc = 54+i; % column num for ave CC of aligning 4th sta and original trio
  ccavethres = 0.8*xcmaxAVEnmin;  %% average cc threshold
  offavethres = 2;    % average differential offset threshold
  
  ind = find(dataplt(:,colflg)==1);
  loff = dataplt(ind,colloff)*sps/40;
  ioff = dataplt(ind,colioff)*sps/40;
  cc = dataplt(ind,colcc);
  off12 = dataplt(ind,7);
  off13 = dataplt(ind,8);
  mloff = median_pixel(off12,off13,loff);
  mioff = median_pixel(off12,off13,ioff);
  mcc = median_pixel(off12,off13,cc);

  if isequaln(dataplt, hfall)
    ran = [-20 20 -75 55];
    ind2 = find(off12>=ran(1) & off12<=ran(2) & off13>=ran(3) & off13<=ran(4));
    off12 = off12(ind2);
    off13 = off13(ind2);
    ioff = ioff(ind2);
    mloff = mloff(mloff(:,1)>=ran(1) & mloff(:,1)<=ran(2) & mloff(:,2)>=ran(3) & mloff(:,2)<=ran(4), :);
    mioff = mioff(mioff(:,1)>=ran(1) & mioff(:,1)<=ran(2) & mioff(:,2)>=ran(3) & mioff(:,2)<=ran(4), :);
    mcc = mcc(mcc(:,1)>=ran(1) & mcc(:,1)<=ran(2) & mcc(:,2)>=ran(3) & mcc(:,2)<=ran(4), :);
  end

  
  %%%plane fitting
  [xData, yData, zData] = prepareSurfaceData(mioff(:,1),mioff(:,2),mioff(:,3));
%   [xData, yData, zData] = prepareSurfaceData(off12,off13,ioff);

  % Set up fittype and options.
  ft = fittype( 'a*x+b*y+c', 'independent', {'x', 'y'}, 'dependent', 'z' );
  opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
  opts.Display = 'Off';
  opts.Robust = 'Bisquare';
  opts.StartPoint = [0.5 0.5 0.5];

  % Fit model to data.
  [fitresult, gof] = fit( [xData, yData], zData, ft, opts );
  modparam(i,:) = [fitresult.a fitresult.b fitresult.c];
  
  off14pred{i} = fitresult.a.*off12+fitresult.b.*off13+fitresult.c;


  % Plot fit with data.
  figure; hold on
  p1 = scatter3(xData, yData, zData, 15,[.5 .5 .5],'filled','MarkerEdgeColor','k');
  p2 = plot(fitresult);
  colormap jet;
  c=colorbar;
  c.Label.String = 'off14';
%   h = plot( fitresult, [xData, yData], zData );
  legend([p1 p2], 'z vs. x, y', 'plane fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
  % Label axes
  xlabel( 'x', 'Interpreter', 'none' );
  ylabel( 'y', 'Interpreter', 'none' );
  zlabel( 'z', 'Interpreter', 'none' );
  title(stasnew(i,:));
  grid on
%   axis equal
%   zlim([-20 10]);
%   view( -74.8, 5.4 );
  view(2);

end



