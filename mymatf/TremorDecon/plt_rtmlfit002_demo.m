% plt_rtmlfit002_demo.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to plot the DEMO of analysing and linear fitting several 
% relatively large scale migrations imaged by detections around 002 region. 
% These RTMs typically migrated through the high-density ellipse around 002. The 
% purpose is to get the a sense of the migration direction, speed, and a
% the time needed to pass the ellipse given its size. 
%
% --Similar to 'miglfit002_4s.m'
% 
% --This is the specific RTM that corresponds with burst #181 used in
% deconvolution.
%
% These information might be useful for determining the appropriate max.
% threshold for the separation between detections to group them into 
% temporally-cluster tremor bursts. 
% 
% --Through the plots, looks like the direction ranges from 115;130;130;145
%   degrees, while the speed through the 002 ellipse is around 20 km/h.
%   Let's take the direction as the median, 130. The time needed to pass
%   through the ellipse is about 0.1 h = 6 m = 360 s.
% --Some useful scale: If speed is 30 km/h, then time needed for the same
%   distance 2 km is about 4 m = 240 s; Given that the location error of 
%   +-1 sample at 40 sps using PGC trio possibily is around +-0.75 km in
%   in the strike direction that is also close to the propagation direction,
%   then it means, change in 1 sample would need 1.5 m = 90 s if speed is 
%   30 km/h and 135 s if speed is 20 km/h
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/04/16
% Last modified date:   2024/04/16
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
temppath = strcat(workpath, '/templates/');
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/PGCtrio');

freqflag='hf';  % flag to indicate whether to do hf or lf;

FLAG = 'PGC'; % detector
  
fam = '002';   % family number

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

PREFIX = strcat(fam,'.up.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),'.',...
                        int2str(npo),int2str(npa),'.ms', int2str(mshift));

%load all detections
hfall = load(strcat(rstpath, '/MAPS/eloc.pgc.tdec_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',...
            num2str(winlen/sps),'s',num2str(sps),'sps4add_'));
daycol = 14;
seccol = 15;  % 2 choices here, one is the center of detecting window, other is the start of strongest .5 s 
hfall = sortrows(hfall, [daycol, seccol]);
%format as follows:
%%% 8+34+4*nstanew cols, if 4 new stas, then will be 58 cols
%%% UPDATED at 2021/06/23
%%%   dx,dy,lon,lat,dep,tori,off2,off3 (integer samples at upsampled sps)
%%%   1:off2(n) 2:off3(n) 3:off2sec 4:off3sec 5:fam 6:date 7:timswin(n)
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
% %Time ranges of RTMs
% trange = [
%   2005255   1.35*3600    1.70*3600;
%   2005255   9.55*3600   10.10*3600;
%   2005255  16.10*3600   16.55*3600;
%   2003063   1.70*3600    2.05*3600;
%   ];

%original burst win ranges
ttol = 35;
ntol = 3;
tranbst = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.',...
  num2str(ntol),'.pgc002.',cutout(1:4)));

%actually-used in decon, ranges with a little buffer 
tranbstbuf = load(strcat(rstpath, '/MAPS/tdec.bstranbuf',...
  num2str(ttol),'s.pgc002.',cutout(1:4)));
tlen = tranbstbuf(:,3)-tranbstbuf(:,2);

ttol = 1e-3*86400;
tranmig = load(strcat(rstpath, '/MAPS/migran',num2str(round(ttol)),'s.pgc002'),'w+');
% tranmig(:,2:3) = tranmig(:,2:3)/3600; 

%%%migrations in Rubin et al. 2013, correspondence is listed as to my deconvolution
%%%burst windows, and to my migration windows
tranmigrubin = [
           2003062 6.15 6.28; %s2b, <->none, <->#3
           2003062 9.02 9.16; %s2c, <->#23, <->#12
           2003062 10.65 10.81; %s2d, <->#31, <->#15
           2003062 20.85 21.25; %s2e, <->#44-47, <->#23
           2003062 21.38 21.58; %s2f, <->#48-49, <->#24
           2003063 1.71 2.04; %s2g, <->#53, <->#27
           2003063 18.03 18.21; %s2h, <->#65-66, <->#37
           2004196 11.097 11.128; %s3a, <->#73, <->none
           2004196 15.57 15.74; %s3b, <->#92-93,, <->none
           2004196 19.55 19.68; %s3c, <->#99, <->#57
           2004196 20.51 20.71; %s3d, <->#100, <->#58 
           2004197 8.5 8.8; %s3e, <->#115-116, <->#69-70
           2004197 10.38 10.72; %s3f, <->#117-118, <->#71-72
           2004197 10.88 11.00; %s3g, <->#119-120, <->#73
           2004197 11.93 12.15; %s3h, <->#121, <->#74
           2005254 7.62 7.79; %fig 6a, <->#146, <->#92
           2005254 8.43 8.80; %fig 6b, <->#147-150, <->#93
           2005254 9.58 9.71; %6c, <->#151-152, <->#94
           2005254 9.77 9.98; %6d, <->#153, <->#95
           2005254 20.43 20.49; %6e, <->#164, <->#98
           2005255 1.34 1.69; %6f, <->#170, <->#100
           2005255 9.57 10.11; %6g, <->#174-175, <->#103
           2005255 16.12 16.55; %fig 4 and 6h, <->#181, <->#109
           ];

%%%corresponding
indplt = [3; 12; 15; 23; 24; 27; 37; 57; 58; 69; 70; 71; 72; 73; 74;
          92; 93; 94; 95; 98; 100; 103; 109];

%range as 4-s tremor density
xran = [-6 4];
yran = [-4 4];

dangle = 5; %increment of trial angles
angle = 0: dangle: 360-dangle;
% slope = zeros(size(indplt,1),length(angle));
% rmse = zeros(size(indplt,1),length(angle));
% angrmse = zeros(size(indplt,1),1);
% angslope = zeros(size(indplt,1),1);
                      
% create fit object with free constraints
fttpfree = fittype( @(a,b,x) a*x+b);

%the cut-out boundary of 4-s detections
cutout = 'ellipse';
x0 = 0.2; 
y0 = 0.2;
semia = 1.75;
semib = 1.25;
angrot = 45;
[xcut, ycut] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);
        
% for i = 1: size(trange,1)
% for i = 23: size(indplt,1)
  i=23;
  % disp(trange(i,:));
  disp(tranmig(indplt(i)));
  ind = find(hfall(:,daycol)==tranmig(indplt(i),1) & ...
    hfall(:,seccol)>=tranmig(indplt(i),2) & ...
    hfall(:,seccol)<=tranmig(indplt(i),3));
  mig = hfall(ind,:);
  mig(:,seccol) = mig(:,seccol)-tranmig(indplt(i),2);
  tlenmig = tranmig(indplt(i),3)-tranmig(indplt(i),2);
  
  %%% find best angle, now is mainly to get the variation of se and slope with the trial angle
  for iang = 1: length(angle)
    %%% propagation trial of hf
    migdum = mig;
    for j = 1: size(migdum,1)
      x0 = migdum(j,1);
      y0 = migdum(j,2);
      [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),[0 0]);
      migdum(j,1) = newx;
      migdum(j,2) = newy;
    end
    % linear robust least square
    [fitobj,gof,~] = fit(migdum(:,seccol), migdum(:,1),fttpfree,'Robust','Bisquare',...
      'StartPoint',[1 1]);
    % output fit parameters
    coef = coeffvalues(fitobj);
    slope(iang) = coef(1);
    rmse(iang) = gof.rmse;    
  end
                      
  %%% best angle estimate from hf
  ind = find(slope>0);
  ind3 = find(rmse(ind)==min(rmse(ind)));     % rmse, Root Mean Squared Error, the smaller, the better
  if length(ind3) > 1
    disp(strcat({'multiple angles have same rmse: '},num2str(angle(ind3))));
  end
  angrmse = angle(ind(ind3(1)));
  
  ind6 = find(slope==max(slope)); % one with the largest slope, i.e., migrating speed
  if length(ind6) > 1
    disp(strcat({'multiple angles have same slope:'},num2str(angle(ind6))));
  end
  angslope = angle(ind6(1));
  
  %%%%%%%%%%%%%% 4-s tremors in the same window as burst 181
  indbst = find(mig(:,seccol)+tranmig(indplt(i),2)>=tranbstbuf(181,2) & ...
                mig(:,seccol)+tranmig(indplt(i),2)<=tranbstbuf(181,3));
  %%% find best angle, now is mainly to get the variation of se and slope with the trial angle
  for iang = 1: length(angle)
    %%% propagation trial of hf
    migdum = mig(indbst,:);
    for j = 1: size(migdum,1)
      x0 = migdum(j,1);
      y0 = migdum(j,2);
      [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),[0 0]);
      migdum(j,1) = newx;
      migdum(j,2) = newy;
    end
    % linear robust least square
    [fitobj,gof,~] = fit(migdum(:,seccol), migdum(:,1),fttpfree,'Robust','Bisquare',...
      'StartPoint',[1 1]);
    % output fit parameters
    coef = coeffvalues(fitobj);
    slope(iang) = coef(1);
    rmse(iang) = gof.rmse;    
  end
  %%% best angle estimate from hf
  ind = find(slope>0);
  ind3 = find(rmse(ind)==min(rmse(ind)));     % rmse, Root Mean Squared Error, the smaller, the better
  if length(ind3) > 1
    disp(strcat({'multiple angles have same rmse: '},num2str(angle(ind3))));
  end
  angrmse1 = angle(ind(ind3(1)));
  
  ind6 = find(slope==max(slope)); % one with the largest slope, i.e., migrating speed
  if length(ind6) > 1
    disp(strcat({'multiple angles have same slope:'},num2str(angle(ind6))));
  end
  angslope1 = angle(ind6(1));

                    
  %% LFE-style plot
  widin = 8.3;  % maximum width allowed is 8.5 inches
  htin = 4;   % maximum height allowed is 11 inches
  nrow = 1;
  ncol = 3;
  f = initfig(widin,htin,nrow,ncol); %initialize fig
  
  axpos = [0.08 0.15 0.42 0.7;
           0.57 0.15 0.40 0.7;
           0.60 0.48 0.16 0.3];
  for isub = 1:nrow*ncol
    set(f.ax(isub), 'position', axpos(isub,:));
  end

  % subplot 1 of figure i
  ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  hold(ax,'on');
  plot(ax,[-100 100],[0 0],'k--');
  plot(ax,[0 0],[-100 100],'k--');
  plot(ax,xcut,ycut,'k-','linew',1.5);
  migdum = sortrows(mig,seccol);
  scatter(ax,migdum(:,1),migdum(:,2), 20, migdum(:,seccol), 'filled','o',...
    'MarkerEdgeColor',[.5 .5 .5]);
  colormap(ax,'jet');
  % colormap(ax,flipud(colormap(ax,'kelicol')));
%   colormap(ax,'viridis');
  c1=colorbar(ax,'SouthOutside');
  pos = ax.Position;
  c1.Position = [0.08, 0.11, 0.42, 0.02];
  %     c.TickLabels=[];
  juldate = num2str(tranmig(indplt(i),1));
  yr = str2double(juldate(1:4));
  date = str2double(juldate(5:end));
  a = jul2dat(yr,date);
  if a(1) == 9
    mo = 'Sep.';
  elseif a(1) == 7
    mo = 'Jul.';
  else
    mo = 'Mar.';
  end
  day = num2str(a(2));
  yr = num2str(a(3));
  c1.Label.String = sprintf('Time (s) since %.1f s on %s %s %s', ...
    tranmig(indplt(i),2),day, mo, yr);
  %     c.Label.String = strcat(num2str(trange(i,1)),' of HF',' (hr)');
  caxis(ax,[0 tlenmig]);
  % text(f.ax(1),0.85,0.1,'HF','FontSize',12,'unit','normalized');
  text(ax,0.02,0.95,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1);
%   text(ax,0.04,0.1,'PGC','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
  axis(ax, 'equal');
  axis(ax,[xran yran]);
  xticks(ax,xran(1):1:xran(2));
  yticks(ax,yran(1):1:yran(2));
  ax.GridLineStyle = '--';
  ax.XAxisLocation = 'top';
  medxhf = median(migdum(:,1));
  medyhf = median(migdum(:,2));
  [rotx, roty] = complex_rot(0,1,-angrmse);
  xarrow = [medxhf-rotx medxhf+rotx];
  yarrow = [medyhf-roty medyhf+roty];
%   drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5);
  p=annotation('arrow','color',[0.5 0.5 0.5],'linestyle','-','linewidth',1.5);
  p.Parent = ax;
  p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
  xlabel(ax,'E (km)');
  ylabel(ax,'N (km)');
  % text(f.ax(1),0.5,0.1,num2str(i),'FontSize',12,'unit','normalized');
  text(ax,0.98,0.22,strcat(num2str(size(migdum,1)),{' events'}),'FontSize',8,...
    'unit','normalized','horizontalalignment','right');
  text(ax,0.98,0.17,strcat({'in '},num2str(tlenmig),{' s'}),'FontSize',8,...
    'unit','normalized','horizontalalignment','right');
  rate = sprintf('%.3f',size(migdum,1)/tlenmig);
  text(ax,0.98,0.12,strcat({'rate: '},rate),'FontSize',8,...
    'unit','normalized','horizontalalignment','right');
  % text(ax,0.98,0.05,strcat(num2str(angrmse),'$^{\,\circ}$'),'interpreter',...
  %   'latex','Units','normalized','FontSize',12,'HorizontalAlignment','right');
  text(ax,0.98,0.05,sprintf('%d%c',angrmse,char(176)),...
    'Units','normalized','HorizontalAlignment','right','FontSize',10);
  text(ax,0.98,0.95,'4-s Tremor','FontSize',11,'unit','normalized','horizontalalignment','right',...
    'fontweight','bold');
  hold(ax,'off');
  
  % subplot 2 of figure i
  ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  %%% Actually used is the pre-determined best prop direc to do the fitting
  migdum = mig;
  for j = 1: size(mig,1)
    x0 = migdum(j,1);
    y0 = migdum(j,2);
    [newx,newy] = coordinate_rot(x0,y0,-(angrmse-90),[0 0]);
    migdum(j,1) = newx;
    migdum(j,2) = newy;
  end
  indelse=setdiff((1:length(migdum))', indbst);  
  scatter(ax,migdum(indelse,seccol),migdum(indelse,1),...
    20,[.6 .6 .6],'filled','o','MarkerEdgeColor','w');
  % create fit object  
  [fitobjprop,gofprop,outprop] = fit(migdum(:,seccol), migdum(:,1),fttpfree,'Robust',...
    'Bisquare','StartPoint',[1 1]);
  %%%Some statistics
  x = migdum(:,seccol);
  y = migdum(:,1);
  stats = statsofrobustlnfit(fitobjprop,gofprop,outprop,x,y);
  %add more statistics into the structure
  stats.angrmse = angrmse;
  stats.angslope = angslope;
  fitprop = feval(fitobjprop,migdum(:,seccol));
  plot(ax,migdum(:,seccol),fitprop,'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
  
  %%%%%%%%tremors that correspond to burst win #181
  migbst = sortrows(mig(indbst,:),seccol);
  timevecbst = migbst(:,seccol);
  [~,~,statsbst] = srcprojdistNtoNm(timevecbst,migbst(:,1:2),1,1);
  migbstdum = migbst;
  for j = 1: size(migbst,1)
    x0 = migbstdum(j,1);
    y0 = migbstdum(j,2);
    [newx,newy] = coordinate_rot(x0,y0,-(statsbst.angrmse-90),[0 0]);
    migbstdum(j,1) = newx;
    migbstdum(j,2) = newy;
  end
  scatter(ax,migbstdum(:,seccol),migbstdum(:,1),20,timevecbst,'filled','o',...
    'MarkerEdgeColor',[.5 .5 .5]);  
  colormap(ax,'viridis');
  c2=colorbar(ax,'SouthOutside');
  pos = ax.Position;
  % c2.Position = [pos(1), 0.11, pos(3), 0.02];
  c2st = pos(1)+(tranbstbuf(181,2)-tranmig(indplt(i),2))/tlenmig*pos(3);
  c2len = tlen(181)/tlenmig*pos(3);
  c2.Position = [c2st, 0.11, c2len, 0.02];
  % caxis(ax,[0 tlen(181)]);
  % c2.Label.String = sprintf('Time (s) since %.1f s on %s %s %s', ...
  %   tranbstbuf(181,2),day, mo, yr);
  caxis(ax,minmax(timevecbst'));
  c2.Ticks = [1050 1250];
  c2.TickLabels = ['1050'; '1250'];
  c2.TickLength = 0.05;
  % c2.Label.String = sprintf('Time (s) since %.1f s on %s %s %s', ...
  %   tranmig(indplt(i),2),day, mo, yr);  
  % create fit object  
  [fitobjprop1] = fit(migbstdum(:,seccol), migbstdum(:,1),fttpfree,'Robust',...
    'Bisquare','StartPoint',[1 1]);
  % output fit parameters
  fitprop1 = feval(fitobjprop1,migbstdum(:,seccol));
  plot(ax,migbstdum(:,seccol),fitprop1,'-','linewidth',2,'color','k');
  xran1 = [0 tlenmig];
  yran1 = [-5 2];
  xlim(ax,xran1);
  ylim(ax,yran1);
  xlabel(ax,sprintf('Time (s) since %.1f s on %s %s %s',tranmig(indplt(i),2),day, mo, yr));
  ylabel(ax,'Projected location (km)');    
  text(ax,0.98,0.95,sprintf('Slope: %.1f m/s',statsbst.slope*1e3),'FontSize',8,...
    'unit','normalized','horizontalalignment','right');
%   text(ax,0.98,0.1,sprintf('SE: %.2f',stats.slopese),'FontSize',8,...
%     'unit','normalized','horizontalalignment','right');
  text(ax,0.98,0.9,sprintf('Pearson: %.2f',statsbst.pearwt),'FontSize',8,...
    'unit','normalized','horizontalalignment','right');
  %%%%%%%%  

  text(ax,0.98,0.1,sprintf('Slope: %.1f m/s',stats.slope*1e3),'FontSize',8,...
    'unit','normalized','horizontalalignment','right');
%   text(ax,0.98,0.1,sprintf('SE: %.2f',stats.slopese),'FontSize',8,...
%     'unit','normalized','horizontalalignment','right');
  text(ax,0.98,0.05,sprintf('Pearson: %.2f',stats.pearwt),'FontSize',8,...
    'unit','normalized','horizontalalignment','right');  
  text(ax,0.02,0.95,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1);
  % ax.GridLineStyle = '--';
  ax.XAxisLocation = 'top';
  hold(ax,'off');

  %an inset on top, map view of tremors that correspond to burst win #181
  timevecbst = migbst(:,seccol)+tranmig(indplt(i),2)-tranbstbuf(181,2);
  ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  hold(ax,'on');
  plot(ax,[-100 100],[0 0],'k--');
  plot(ax,[0 0],[-100 100],'k--');
  plot(ax,xcut,ycut,'k-','linew',1.5);
  scatter(ax,migbst(:,1),migbst(:,2), 15, timevecbst, 'filled','o',...
    'MarkerEdgeColor',[.5 .5 .5]);
  colormap(ax,'viridis');
  caxis(ax,[0 tlen(181)]);
  [rotx, roty] = complex_rot(0,1,-statsbst.angrmse);
  xarrow = [1-rotx 1+rotx];
  yarrow = [-1-roty -1+roty];
  p=annotation('arrow','color',[0.5 0.5 0.5],'linestyle','-','linewidth',1.5);
  p.Parent = ax;
  p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
  % text(ax,0.98,0.1,strcat(num2str(statsbst.angrmse),'$^{\,\circ}$'),'interpreter',...
  %   'latex','Units','normalized','FontSize',11,'HorizontalAlignment','right');
  text(ax,0.98,0.1,sprintf('%d%c',statsbst.angrmse,char(176)),...
    'Units','normalized','HorizontalAlignment','right','FontSize',10);
  text(ax,0.05,0.12,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1);
  axis(ax, 'equal');
  axis(ax,[-2.5 2.5 -2.5 2.5]);
  % axis(ax,[-4 4 -4 4]);
  xticks(ax,xran(1):1:xran(2));
  yticks(ax,yran(1):1:yran(2));
  ax.GridLineStyle = '--';
  % text(ax,0.98,0.22,strcat(num2str(size(migdum,1)),{' events'}),'FontSize',8,...
  %   'unit','normalized','horizontalalignment','right');
  nolabels(ax,3);
  longticks(ax,1.5);
  hold(ax,'off');


  % orient(f.fig,'landscape');
  fname = '4srtmdemo.pdf';
  print(f.fig,'-dpdf',...
    strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

%   keyboard
  

