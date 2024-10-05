% empioffset4thsta002.m
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
% --to smooth out these spikes, there is a testing script 
%   'testempioffset4thsta002' that tries a few ways, like 2D 
%   interpolation, gradient estimate, and 2D plane fitting. The test shows
%   that interpolation does not work at all; the gradient might not be very
%   stable; plane fitting seems like the best solution so far in terms of 
%   having a 1st-order estimate for the time offset. Therefore, this script
%   is going to use plane fitting only. 
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

sps = 160;

%%
%%%obtain the mean of median of off12 (t1-t2), off13 (t1-t3), and
%%%off23 (t2-t3==off13-off12). off12+off23-off13 == 0
offtrio{1} = [hfall(:,1) hfall(:,2) hfall(:,7)]; %off12 for each loce and locn
offtrio{2} = [hfall(:,1) hfall(:,2) hfall(:,8)]; %off13 for each loce and locn
offtrio{3} = [hfall(:,1) hfall(:,2) hfall(:,8)-hfall(:,7)]; %off23 for each loce and locn
avetrio{1} = median_pixel(hfall(:,1),hfall(:,2),hfall(:,7),'median'); %median off12 at duplicate loc
avetrio{2} = median_pixel(hfall(:,1),hfall(:,2),hfall(:,8),'median'); %median off13 at duplicate loc
avetrio{3} = median_pixel(hfall(:,1),hfall(:,2),hfall(:,8)-hfall(:,7),'median'); %median off23 at duplicate loc
%%
%simply scatter the detections to check the spatial resolution
xran = [-8 3];
yran = [-4 4];
f = initfig(8.5,4,1,3); %initialize fig
for i = 1: 3
  ax=f.ax(i); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  ave = avetrio{i};
  scatter(ax,ave(:,1),ave(:,2),10,ave(:,3),'filled');
  colormap(ax,'jet');
  axis(ax,'equal');
  axis(ax,[xran yran]);
  c=colorbar(ax);
  caxis(ax,[-15 15]*sps/40);
  box on
  if i==1
    title(ax,'PGC-SSIB (off12)');
  elseif i==2
    title(ax,'PGC-SILB (off13)');
  else
    title(ax,'SSIB-SILB (off23)');
  end
  longticks(ax,1);
  if i == 1
    xlabel(ax,'E (km)');
    ylabel(ax,'N (km)');
    c.Label.String = sprintf('mean offset (samples) at %d Hz',sps);
  end

end

%%
modparammaptrio = [];
for i = 1: 3
  off = offtrio{i};
  [xData, yData, zData] = prepareSurfaceData(off(:,1),off(:,2),off(:,3)/sps);
  % Set up fittype and options.
  opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
  opts.Display = 'Off';
  opts.Robust = 'Bisquare';

  % Fit model to data.
  %%%plane fit model
  ftplane = fittype( 'a*x+b*y+c', 'independent', {'x', 'y'}, 'dependent', 'z' );
  opts.StartPoint = [0.1 0.1 0.1];
  [fitresult, gofmaptrio{i},outputmaptrio{i}] = fit( [xData, yData], zData, ftplane, opts );
  modparammaptrio(i,:) = [fitresult.a fitresult.b fitresult.c];
  offpred = fitresult.a .* xData + fitresult.b .* yData + fitresult.c;
  
  %test if the returned goodness-of-fit metrics are also weighted, answer is YES 
  sseauto = gofmaptrio{i}.sse;  %sum of squares due to error
  ssecalc = sum((zData-offpred).^2); %assuming equal weights of 1
  rmseauto = gofmaptrio{i}.rmse;  %root mean squared error
  rmsecalc = sqrt(sum((zData-offpred).^2)./(length(zData)-2));  %assuming equal weights of 1

  %%% Plot fit with data.
  widin = 8;  % maximum width allowed is 8.5 inches
  htin = 4.5;   % maximum height allowed is 11 inches
  nrow = 1;
  ncol = 2;
  f=initfig(widin,htin,nrow,ncol);

  figxran = [0.1 0.95]; figyran = [0.12 0.95];
  figxsep = 0.1; figysep = 0.07;
  axloc=optaxpos(f,nrow,ncol,figxran,figyran,figxsep,figysep);

  ax = f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  p1 = scatter3(xData, yData, zData, 15,[.3 .3 .3],'filled','MarkerEdgeColor','w');
  %     aa=feval(fitresult,xData,yData);
  %     surf(ax,xData, yData,aa);
  p2 = plot(fitresult);
  colormap('jet');
  %     c=colorbar;
  %     c.Label.String = sprintf('off14 at %d sps',sps);
  %   h = plot( fitresult, [xData, yData], zData );
  legend([p1 p2], sprintf('\\Delta{t}'), 'Plane fit', 'Location', 'northwest');
  xlabel('E (km)');
  ylabel('N (km)');
  zlabel(sprintf('\\Delta{t}_{14} (s)'));
  zlim([-1.5 1.5]);
  xlim([-10 10]);
  ylim([-10 10]);
  view(108,10);
  %     view(2);
  
  ax = f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  error = zData-offpred;
  binw=2/sps;
  histogram(ax,error,'Normalization','count','BinWidth',binw,...
    'Facec','k','edgec','none');
  p1=plot(ax,[rmseauto rmseauto],ax.YLim,'r--','linew',1.5);
  plot(ax,[-rmseauto -rmseauto],ax.YLim,'r--','linew',1.5);
  text(ax,0.98,0.95,sprintf('%s',strtrim(stas(i,:))),...
    'Units','normalized','HorizontalAlignment','right','fontsize',12);
  xran = [-0.5,0.5];
  %     [muHat,sigmaHat] = normfit(error);
  %     pdfhat = normpdf(xran(1):binw:xran(2),muHat,sigmaHat)*binw;
  %     p3=plot(ax,xran(1):binw:xran(2),pdfhat,'r','linew',2);
  xlabel(ax,sprintf('Deviation from plane fitting (s)'));
  ylabel(ax,'Count');
  xlim(ax,xran);
  xticks(ax,xran(1):0.25:xran(2));
end

%% gradient
doffdloc = zeros(nsta,2);
doffdprojloc = zeros(nsta,2);
for i = 1: nsta
  a = modparammaptrio(i,1);
  b = modparammaptrio(i,2);
  doffdloc(i,1) = a;
  doffdloc(i,2) = b;
  projang = 135;
  theta = -(projang-90);
  doffdprojloc(i,1) = a*cosd(theta)+b*sind(theta);
  doffdprojloc(i,2) = b*cosd(theta)-a*sind(theta);  
end
  
sps = 160;
doffspldloc = doffdloc*sps;
doffspldprojloc = doffdprojloc*sps;

savefile = 'timeoff_plfitmap_160sps.mat';
save(strcat(rstpath, '/MAPS/',savefile), 'modparammaptrio','gofmaptrio',...
  'outputmaptrio','doffdloc','doffdprojloc');

%% extra settings for 4th stations
stasnew=['LZB  '
         'TWKB '
         'MGCB '
         'KLNB '];  % twkb lzb mgcb
nstasnew = size(stasnew,1); 

data = hfall;

flagplt = 1;

if isequaln(data, hfall)
  sybsz = 5; 
elseif isequaln(data, hfbnd)
  sybsz = 20;
end

xran = [min(data(:,7))-1 max(data(:,7))+1];
yran = [min(data(:,8))-1 max(data(:,8))+1];
xvec = xran(1): 1: xran(2);
yvec = yran(1): 1: yran(2);


%% raw data
if flagplt

  f1 = initfig(16,6,2,4);
  f2 = initfig(16,6,2,4);
  f3 = initfig(16,6,2,4);

  for i = 4: nstasnew
    colflg = 42+i;  % column num for the flag to indicate if the 4th sta passes the check
    colloff = 46+i; % column num for the summed difference in offset needed to align 4th sta and original trio
    colioff = 50+i; % column num for the offset needed to align 4th sta and original trio
    colcc = 54+i; % column num for ave CC of aligning 4th sta and original trio
    ccavethres = 0.8*xcmaxAVEnmin;  %% average cc threshold
    offavethres = 2*sps/40;    % average differential offset threshold

    ind = find(data(:,colflg)==1);
    loff = data(ind,colloff)*sps/40;
    ioff = data(ind,colioff)*sps/40;
    cc = data(ind,colcc);
    off12 = data(ind,7);
    off13 = data(ind,8);
    mloff = median_pixel(off12,off13,loff);
    mioff = median_pixel(off12,off13,ioff);
    mcc = median_pixel(off12,off13,cc);

    ran = [-30 50 -60 40];

    %%%figure 1, averaged difference in offset between 3 pairs, 14,24,34
    %median
    ax = f1.ax(i);
    hold(ax,'on'); grid(ax,'on');
    scatter(ax,data(:,7),data(:,8),sybsz,[.8 .8 .8],'filled');
    scatter(ax,mloff(:,1),mloff(:,2),sybsz,mloff(:,3),'filled');
    text(ax,0.98,0.15,sprintf('%d/%d',length(ind),size(data,1)),'Units','normalized',...
      'HorizontalAlignment','right');
    text(ax,0.98,0.05,sprintf('=%.1f%%',length(ind)/size(data,1)*100),'Units','normalized',...
      'HorizontalAlignment','right');
    c=colorbar(ax);
    c.Label.String = 'med of averaged offset diff.';
    %  colormap(ax,flipud(colormap(ax,'gray')));
    colormap(ax,'jet');
    caxis(ax,[0 offavethres]);
    axis(ax,[xran yran],'equal');
    text(ax,0.02,0.95,stasnew(i,:),'Units','normalized','HorizontalAlignment','left');
    hold(ax,'off');

    %standard deviation
    ax = f1.ax(4+i);
    hold(ax,'on'); grid(ax,'on');
    idx = find(mloff(:,end)~=0);  %only show non-zeros, ie., has duplicates
  %   idx = 1:size(mloff,1);
    scatter(ax,mloff(idx,1),mloff(idx,2),sybsz,mloff(idx,end),'filled','MarkerEdgeColor','none');
    c=colorbar(ax);
    c.Label.String = 'std of averaged offset diff.';
  %   colormap(ax,flipud(colormap(ax,'gray')));
    colormap(ax,'jet');
    caxis(ax,[0 prctile(mloff(idx,end),95)]);
    axis(ax,'equal',ran);
    hold(ax,'off');


    %%%figure 2, observational weighted off14 (equivalent to 24 and 34)
    %median
    ax = f2.ax(i);
    hold(ax,'on'); grid(ax,'on');
    scatter(ax,mioff(:,1),mioff(:,2),sybsz,mioff(:,3),'filled');
    c=colorbar(ax);
    c.Label.String = 'med of weighted off14';
    colormap(ax,'jet');
    caxis(ax,[prctile(mioff(:,3),5), prctile(mioff(:,3),95)]);
    axis(ax,[xran yran],'equal');
    text(ax,0.02,0.95,stasnew(i,:),'Units','normalized','HorizontalAlignment','left');
    hold(ax,'off');

    %standard deviation
    ax = f2.ax(4+i);
    hold(ax,'on'); grid(ax,'on');
    idx = find(mloff(:,end)~=0);  %only show non-zeros, ie., has duplicates
  %   idx = 1:size(mloff,1);
    scatter(ax,mioff(idx,1),mioff(idx,2),sybsz,mioff(idx,end),'filled','MarkerEdgeColor','none');
    c=colorbar(ax);
    c.Label.String = 'std of weighted off14';
  %   colormap(ax,flipud(colormap(ax,'gray')));
    colormap(ax,'jet');
    caxis(ax,[0 prctile(mioff(idx,end),95)]);
    axis(ax,'equal',ran);
    hold(ax,'off');


    %%%figure 3, averaged CC
    %median
    ax = f3.ax(i);
    hold(ax,'on'); grid(ax,'on');
    scatter(ax,mcc(:,1),mcc(:,2),sybsz,mcc(:,3),'filled');
    c=colorbar(ax);
    c.Label.String = 'med of averaged CC';
    colormap(ax,'jet');
    caxis(ax,[0 prctile(mcc(:,3),95)]);
    axis(ax,[xran yran],'equal');
    text(ax,0.02,0.95,stasnew(i,:),'Units','normalized','HorizontalAlignment','left');
    hold(ax,'off');

    %standard deviation
    ax = f3.ax(4+i);
    hold(ax,'on'); grid(ax,'on');
    idx = find(mloff(:,end)~=0);  %only show non-zeros, ie., has duplicates
  %   idx = 1:size(mloff,1);
    scatter(ax,mcc(idx,1),mcc(idx,2),sybsz,mcc(idx,end),'filled','MarkerEdgeColor','none');
    c=colorbar(ax);
    c.Label.String = 'std of averaged CC';
  %   colormap(ax,flipud(colormap(ax,'gray')));
    colormap(ax,'jet');
    caxis(ax,[0 prctile(mcc(idx,end),95)]);
    axis(ax,'equal',ran);
    hold(ax,'off');
  end
  xlabel(f1.ax(1),'off12 (160 sps)');
  ylabel(f1.ax(1),'off13 (160 sps)');
  title(f1.ax(1),'averaged offset diff.');
  xlabel(f2.ax(1),'off12 (160 sps)');
  ylabel(f2.ax(1),'off13 (160 sps)');
  title(f2.ax(1),'weighted off14');
  xlabel(f3.ax(1),'off12 (160 sps)');
  ylabel(f3.ax(1),'off13 (160 sps)');
  title(f3.ax(1),'averaged CC');

end

%% 2D plane fitting
modparam = [];
modparammap = [];
off14pred = [];
for i = 4: nstasnew
  colflg = 42+i;  % column num for the flag to indicate if the 4th sta passes the check, 43-46
  colloff = 46+i; % column num for the summed difference in offset needed to align 4th sta and original trio, 47-50
  colioff = 50+i; % column num for the offset needed to align 4th sta and original trio, 51-54
  colcc = 54+i; % column num for ave CC of aligning 4th sta and original trio, 55-58
  ccavethres = 0.8*xcmaxAVEnmin;  %% average cc threshold
  offavethres = 2*sps/40;    % average differential offset threshold
  
  ind = find(data(:,colflg)==1);
  loff = data(ind,colloff)*sps/40;  %convert to the same sps 
  ioff = data(ind,colioff)*sps/40;  %convert to the same sps, essentially 'off14'
  cc = data(ind,colcc);
  off12 = data(ind,7);  %note that the location off12 and off13 are already consistent here
  off13 = data(ind,8);
  loce = data(ind,1);
  locn = data(ind,2);
  mloff = median_pixel(off12,off13,loff);
  mioff = median_pixel(off12,off13,ioff);
  mcc = median_pixel(off12,off13,cc);

%   if isequaln(data, hfall)
%     if i==1
%       ran = [-30 20 -75 35];
%     elseif i==2
%       ran = [-30 25 -70 40];
%     elseif i==3
%       ran = [-30 25 -65 40];
%     elseif i==4
%       ran = [-40 45 -70 50];
%     end
%     ind2 = find(off12>=ran(1) & off12<=ran(2) & off13>=ran(3) & off13<=ran(4));
%     off12 = off12(ind2);
%     off13 = off13(ind2);
%     ioff = ioff(ind2);
%     mloff = mloff(mloff(:,1)>=ran(1) & mloff(:,1)<=ran(2) & mloff(:,2)>=ran(3) & mloff(:,2)<=ran(4), :);
%     mioff = mioff(mioff(:,1)>=ran(1) & mioff(:,1)<=ran(2) & mioff(:,2)>=ran(3) & mioff(:,2)<=ran(4), :);
%     mcc = mcc(mcc(:,1)>=ran(1) & mcc(:,1)<=ran(2) & mcc(:,2)>=ran(3) & mcc(:,2)<=ran(4), :);
%   end

  
  %%%plane fitting
%   [xData, yData, zData] = prepareSurfaceData(mioff(:,1),mioff(:,2),mioff(:,3));
  [xData, yData, zData] = prepareSurfaceData(off12/sps,off13/sps,ioff/sps);

  % Set up fittype and options.
  opts = fitoptions( 'Method', 'NonlinearLeastSquares' ); %
  opts.Display = 'Off';
  opts.Robust = 'Bisquare';
  opts.StartPoint = [0.8 0.8 0.1];
  
  % Fit model to data.
  %%%plane fit model
  ftplane = fittype( 'a*x+b*y+c', 'independent', {'x', 'y'}, 'dependent', 'z' );
  [fitresult, gof{i},output{i}] = fit( [xData, yData], zData, ftplane, opts );
  modparam(i,:) = [fitresult.a fitresult.b fitresult.c];
  off14pred = fitresult.a .* xData + fitresult.b .* yData + fitresult.c;
  
  %test if the returned goodness-of-fit metrics are also weighted, answer is YES 
  sseauto = gof{i}.sse;  %sum of squares due to error
  ssecalc = sum((zData-off14pred).^2); %assuming equal weights of 1
  rmseauto = gof{i}.rmse;  %root mean squared error
  rmsecalc = sqrt(sum((zData-off14pred).^2)./(length(zData)-2));  %assuming equal weights of 1
  
%   %%%parabolic fit model
%   ftpara = fittype( '(x-x0)^2/a^2+(y-y0)^2/b^2+z0', 'independent', {'x', 'y'}, 'dependent', 'z' );
%   opts.StartPoint = [0.5 0.5 0.5 0.5 0.5];
%   [fitresult, gof] = fit( [xData, yData], zData, ftpara, opts );
%   modparam(i,:) = [fitresult.a fitresult.b fitresult.x0 fitresult.y0 fitresult.z0];
%   off14pred{i} = (off12-(fitresult.x0)).^2./(fitresult.a)^2+(off13-(fitresult.y0)).^2./(fitresult.b)^2+...
%     (fitresult.z0);

  %%% Plot fit with data.
  widin = 8;  % maximum width allowed is 8.5 inches
  htin = 4.5;   % maximum height allowed is 11 inches
  nrow = 1;
  ncol = 2;
  f=initfig(widin,htin,nrow,ncol);

%   figxran = [0.1 0.95]; figyran = [0.12 0.95];
%   figxsep = 0.1; figysep = 0.07;
%   axloc=optaxpos(f,nrow,ncol,figxran,figyran,figxsep,figysep);

  ax = f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  p1 = scatter3(xData, yData, zData, 15,[.3 .3 .3],'filled','MarkerEdgeColor','w');
  %     aa=feval(fitresult,xData,yData);
  %     surf(ax,xData, yData,aa);
  p2 = plot(fitresult);
  colormap('jet');
  %     c=colorbar;
  %     c.Label.String = sprintf('off14 at %d sps',sps);
  legend([p1 p2], sprintf('\\Delta{t}'), 'Plane fit', 'Location', 'northwest');
  xlabel(sprintf('\\Delta{t}_{12}^{win} (s)'));
  ylabel(sprintf('\\Delta{t}_{13}^{win} (s)'));
  zlabel(sprintf('\\Delta{t}_{14}^{win} (s)'));
  axis('equal');
  zlim([-1.5 1.5]);
  xlim([-1 1]);
  ylim([-1 1]);
  text(0.02,0.10,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1);
  view(243,5.3);
  %     view(2);

  ax = f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  error = zData-off14pred;
  binw=1/sps;
  histogram(ax,error,'Normalization','count','BinWidth',binw,...
    'Facec','k','edgec','none');
  p1=plot(ax,[mean(error) mean(error)],ax.YLim,'r--','linew',2);
  p2=plot(ax,[rmseauto rmseauto],ax.YLim,'k--','linew',1.5);
  plot(ax,[-rmseauto -rmseauto],ax.YLim,'k--','linew',1.5);
  text(ax,0.98,0.95,sprintf('%s',strtrim(stasnew(i,:))),...
    'Units','normalized','HorizontalAlignment','right','fontsize',12);
  xran = [-0.25,0.25];
  %     [muHat,sigmaHat] = normfit(error);
  %     pdfhat = normpdf(xran(1):binw:xran(2),muHat,sigmaHat)*binw;
  %     p3=plot(ax,xran(1):binw:xran(2),pdfhat,'r','linew',2);
  xlabel(ax,sprintf('Deviation from plane fit (s)'));
  ylabel(ax,'Count');
  xlim(ax,xran);
  xticks(ax,xran(1):0.25:xran(2));
  text(ax,0.02,0.10,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1);
  legend([p1 p2], 'Mean', '\pm1 RMSE', 'Location', 'northwest');
  
  fname = strcat('planefit',strtrim(stasnew(i,:)),'.pdf');
  print(f.fig,'-dpdf',...
    strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

  keyboard
  
%   %%%%%%%%%%%%%%%%%%%%%%%%%55
%   %%%What if you force 'off14' is 0 at (0,0)  
%   ftplane = fittype( 'a*x+b*y', 'independent', {'x', 'y'}, 'dependent', 'z' );
%   opts = fitoptions( 'Method', 'NonlinearLeastSquares' ); %
%   opts.Display = 'Off';
%   opts.Robust = 'Bisquare';
%   opts.StartPoint = [0.8 0.8];
%   [fitresult] = fit( [xData, yData], zData, ftplane, opts );
%   off14pred = fitresult.a .* xData + fitresult.b .* yData;
%   f=initfig;
%   ax = f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%   error = zData-off14pred;
%   binw=2/sps;
%   histogram(ax,error,'Normalization','count','BinWidth',binw,...
%     'Facec','k','edgec','none');
%   p1=plot(ax,[rmseauto rmseauto],ax.YLim,'r--','linew',1.5);
%   plot(ax,[-rmseauto -rmseauto],ax.YLim,'r--','linew',1.5);
%   text(ax,0.98,0.95,sprintf('%s',strtrim(stasnew(i,:))),...
%     'Units','normalized','HorizontalAlignment','right','fontsize',12);
%   xran = [-0.5,0.5];
%   %     [muHat,sigmaHat] = normfit(error);
%   %     pdfhat = normpdf(xran(1):binw:xran(2),muHat,sigmaHat)*binw;
%   %     p3=plot(ax,xran(1):binw:xran(2),pdfhat,'r','linew',2);
%   xlabel(ax,sprintf('Deviation from plane fitting (s)'));
%   ylabel(ax,'Count');
%   xlim(ax,xran);
%   xticks(ax,xran(1):0.25:xran(2));
%   %%%%%%%%%%%%%%%%%%%%%%%%%55
%   
% 
%   %%%%%%%%%%%%%%%%%%%%%%%%%55
%   %%%% use map locations instead
%   [xData, yData, zData] = prepareSurfaceData(loce,locn,ioff/sps);
%   % Set up fittype and options.
%   opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%   opts.Display = 'Off';
%   opts.Robust = 'Bisquare';
% 
%   % Fit model to data.
%   %%%plane fit model
%   ftplane = fittype( 'a*x+b*y+c', 'independent', {'x', 'y'}, 'dependent', 'z' );
%   opts.StartPoint = [0.1 0.1 0.1];
%   [fitresult, gofmap{i},outputmap{i}] = fit( [xData, yData], zData, ftplane, opts );
%   modparammap(i,:) = [fitresult.a fitresult.b fitresult.c];
%   off14pred = fitresult.a .* xData + fitresult.b .* yData + fitresult.c;
%   
%   %test if the returned goodness-of-fit metrics are also weighted, answer is YES 
%   sseauto = gofmap{i}.sse;  %sum of squares due to error
%   ssecalc = sum((zData-off14pred).^2); %assuming equal weights of 1
%   rmseauto = gofmap{i}.rmse;  %root mean squared error
%   rmsecalc = sqrt(sum((zData-off14pred).^2)./(length(zData)-2));  %assuming equal weights of 1
%   %%%%%%%%%%%%%%%%%%%%%%%%%55
  
%   %%% Plot fit with data.
%   widin = 8;  % maximum width allowed is 8.5 inches
%   htin = 4.5;   % maximum height allowed is 11 inches
%   nrow = 1;
%   ncol = 2;
%   f=initfig(widin,htin,nrow,ncol);
% 
%   figxran = [0.1 0.95]; figyran = [0.12 0.95];
%   figxsep = 0.1; figysep = 0.07;
%   axloc=optaxpos(f,nrow,ncol,figxran,figyran,figxsep,figysep);
% 
%   ax = f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%   p1 = scatter3(xData, yData, zData, 15,[.3 .3 .3],'filled','MarkerEdgeColor','w');
%   %     aa=feval(fitresult,xData,yData);
%   %     surf(ax,xData, yData,aa);
%   p2 = plot(fitresult);
%   colormap('jet');
%   %     c=colorbar;
%   %     c.Label.String = sprintf('off14 at %d sps',sps);
%   %   h = plot( fitresult, [xData, yData], zData );
%   legend([p1 p2], sprintf('\\Delta{t}'), 'Plane fit', 'Location', 'northwest');
%   xlabel('E (km)');
%   ylabel('N (km)');
%   zlabel(sprintf('\\Delta{t}_{14} (s)'));
%   zlim([-1.5 1.5]);
%   xlim([-10 10]);
%   ylim([-10 10]);
%   view(108,10);
%   %     view(2);
%   
%   ax = f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%   error = zData-off14pred;
%   binw=2/sps;
%   histogram(ax,error,'Normalization','count','BinWidth',binw,...
%     'Facec','k','edgec','none');
%   p1=plot(ax,[rmseauto rmseauto],ax.YLim,'r--','linew',1.5);
%   plot(ax,[-rmseauto -rmseauto],ax.YLim,'r--','linew',1.5);
%   text(ax,0.98,0.95,sprintf('%s',strtrim(stasnew(i,:))),...
%     'Units','normalized','HorizontalAlignment','right','fontsize',12);
%   xran = [-0.5,0.5];
%   %     [muHat,sigmaHat] = normfit(error);
%   %     pdfhat = normpdf(xran(1):binw:xran(2),muHat,sigmaHat)*binw;
%   %     p3=plot(ax,xran(1):binw:xran(2),pdfhat,'r','linew',2);
%   xlabel(ax,sprintf('Deviation from plane fitting (s)'));
%   ylabel(ax,'Count');
%   xlim(ax,xran);
%   xticks(ax,xran(1):0.25:xran(2));


%%%
%   figure;
%   hold on; box on; grid on
%   ax = gca;
%   yyaxis(ax,'left');
%   error = zData-off14pred;
%   histogram(ax,abs(error),'Normalization','probability','BinWidth',1*sps/40);
%   p1=plot(ax,[rmseauto rmseauto],ax.YLim,'b--');
%   plot(ax,[-rmseauto -rmseauto],ax.YLim,'b--');
% %   p2=plot(ax,[rmsecalc rmsecalc],ax.YLim,'k--');
% %   plot(ax,[-rmsecalc -rmsecalc],ax.YLim,'k--');
%   xran = [0,50];
% %   [muHat,sigmaHat] = normfit(error);
% %   pdfhat = normpdf(xran(1):xran(2),muHat,sigmaHat);
% %   p3=plot(ax,xran(1):xran(2),pdfhat,'r','linew',2);
%   xlabel('Abs diff in samples between plane fitting and observation');
% %   ylabel('PDF');
% %   legend(ax,[p1,p2,p3],'RMSE using actual weights','RMSE using weights of 1','Gaussian fitting');
%   ylabel('Probability');
%   xlim(ax,xran);
% 
%   yyaxis(ax,'right');
%   [cdfval,x] = ecdf(abs(error)); %between Nth and (N-1)th source
%   plot(ax,x,cdfval,'linew',1,'color','r');
% %   xint = xran(1):xran(2);
% %   [~,ind,~] = unique(x);
% %   cdfint = interp1(x(ind),cdfval(ind),xint,'spline');
% %   plot(ax,xint,cdfint,'linew',1,'color','r');
%   x80 = x(abs(cdfval-0.8)==min(abs(cdfval-0.8)));
%   x75 = x(abs(cdfval-0.75)==min(abs(cdfval-0.75)));
%   x50 = x(abs(cdfval-0.50)==min(abs(cdfval-0.50)));
%   frac = sum(abs(error)<=rmseauto)/length(error)*100;
%   text(ax,0.5,0.9,sprintf('<=%.1f: %.1f%%',x80,80),'Units',...
%     'normalized');
%   text(ax,0.5,0.85,sprintf('<=%.1f: %.1f%%',x75,75),'Units',...
%     'normalized');
%   text(ax,0.5,0.8,sprintf('<=%.1f: %.1f%%',rmseauto,frac),'Units',...
%     'normalized');
%   text(ax,0.5,0.75,sprintf('<=%.1f: %.1f%%',x50,50),'Units',...
%     'normalized');
%   ylabel(ax,'Empirical CDF');
%   legend(ax,[p1],'RMSE using actual weights','Location','east');
  
end

% keyboard
% fid = fopen(strcat(rstpath, '/MAPS/timeoff_planefitparam_4thsta_160sps'),'w+');
% fprintf(fid,'%.4f %.4f %.4f \n',modparam');
% fclose(fid);

% savefile = 'timeoff_plfit_4thsta_160sps.mat';
% save(strcat(rstpath, '/MAPS/',savefile), 'modparam','gof','output');

% %%
% aaa = 1.8+2*randn(1000,1);
% figure; hold on; box on; grid on;
% histogram(aaa,'Normalization','pdf');
% [muHat,sigmaHat] = normfit(aaa);
% xran = [-5,5];
% pdfhat = normpdf(xran(1):0.01:xran(2),muHat,sigmaHat);
% p3=plot(xran(1):0.01:xran(2),pdfhat,'r','linew',2);


%%
%%%if off14 = a*x+b*y+c, then d(off14)/dx = a, d(off14)/dy = b, where x and y
%%%are map locations in N and E in km. If you instead want to know the gradient
%%%along a custom coordinate frame x' and y' where x' is rotated from x by
%%%theta, so [x; y] = [cos(theta) -sin(theta); sin(theta) cos(theta)]* [x'; y'].
%%%Then, the new gradient is d(off14)/dx' = a*cos(theta)+b*sin(theta);
%%%d(off14)/dy' = b*cos(theta)-a*sin(theta);
doff14dloc = zeros(nstasnew,2);
doff14dprojloc = zeros(nstasnew,2);
for i = 1: nstasnew
  a = modparammap(i,1);
  b = modparammap(i,2);
  doff14dloc(i,1) = a;
  doff14dloc(i,2) = b;
  projang = 135;
  theta = -(projang-90);
  doff14dprojloc(i,1) = a*cosd(theta)+b*sind(theta);
  doff14dprojloc(i,2) = b*cosd(theta)-a*sind(theta);  
end
  
sps = 160;
doff14spldloc = doff14dloc*sps;
doff14spldprojloc = doff14dprojloc*sps;

savefile = 'timeoff_plfitmap_4thsta_160sps.mat';
save(strcat(rstpath, '/MAPS/',savefile), 'modparammap','gofmap','outputmap',...
  'doff14dloc','doff14dprojloc');







