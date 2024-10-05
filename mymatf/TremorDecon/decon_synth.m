% decon_synth.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --This script now becomes the driving script to run the deconvolution to
% synthetic seismograms that are generated from different region sizes and
% different saturation levels. Parameters related to which synthetics are 
% read and deconvolved, and those related to the settings of deconvolution 
% are defined here.
% --Moreover, this script analyzes the statistics of the resulting deconvolved
% catalogs. These statistics are a bit out-of-date given our new understanding
% to data. 
% --See also 'decon_synth_onespot.m' for decon to synthetics generated from  
% diff noise and saturation levels from the same spot.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/05/09
% Last modified date:   2024/05/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
format short e   % Set the format to 5-digit floating point
clear
clc
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');

fam = '002';   % family number

stas=['PGC  '
  'SSIB '
  'SILB '
  %   'LZB  '
  %   'TWKB '
  %   'MGCB '
  'KLNB '
  ]; % determine the trio and order, here the 1st sta is PGC
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

%%%load the empirical model param of time offset 14 for each addtional stations
off14mod = load(strcat(rstpath, '/MAPS/timeoff_planefitparam_4thsta_160sps'),'w+');

%%%practically used ranges of tremor bursts w/i buffer
cutout = 'ellipse';
ttol = 35;
ntol = 3;
trange = load(strcat(rstpath, '/MAPS/tdec.bstranbuf',num2str(ttol),'s.pgc002.',cutout(1:4)));
nbst = size(trange,1);
tlen = trange(:,3)-trange(:,2);
tlensum = sum(tlen);

%% implement deconvolution
tdura = 0.25;  %must be consistent with what synthetics actually used
% tdura = 0.4;

normflag = 0; %whether to normalize templates

%which templates to use
tempflag = 'chao';
  
%%%specify distribution for source location
distr='UN';  % uniform distribution
  
%%%specify distribution for source location
% distrloc = 'custompdf'; %using a custom PDF function
distrloc = 'uniform'; %uniformly random in a specified region,
  
%The ratio of elsewhere events to local events.  Zero for Discrete Ide.
fracelsew=0; %0.25 %0.5; %0.6; 
  
%%%specify if considering the physical size of each source
% physicalsize = 1;
physicalsize = 0;
  
%%%diameter of physical size
if physicalsize
  diam=0.15;% 0.5; %0.6; %
else
  diam=0;
end
  
%%%specify regime for transformation from time offset to map location
% ftrans = 'interpArmb';
% ftrans = 'interpArmbreloc';
ftrans = 'interpchao';

timetype = 'tori';
  
%%%specify if forcing a min speration of arrival time for 2 events from the same spot
% forcesep = 1;
forcesep = 0;
    
%%%flag for validating if ground truth of sources can recover the record
%   testsrcflag = 1;
testsrcflag = 0;
  
%%%flag for validing if the spectral shapes of data and templates are similar
% testfreqflag = 1;
testfreqflag = 0;
  
%%%flag for plot the data
%   pltdataflag = 1;
pltdataflag = 0;
  
%%%flag for plot the ground truth distribution
%   pltgtflag = 1;
pltgtflag = 0;
  
%%%flag for plot the decon src distribution after grouping
%   pltsrcflag1 = 1;
%   pltsrcflag1 = 0;
  
%%%flag for plot the decon src distribution after removing 2ndary src
%   pltsrcflag = 1;
pltsrcflag = 0;
  
%%%flag for plot the decon src distribution after checking at 4th stas
%   pltsrc4thflag = 1;
pltsrc4thflag = 0;
  
%%%specify shape of the source region
srcregion='ellipse';
% srcregion='rectangle';
% srcregion='circle';
  
%variation of source region size
if strcmp(srcregion,'ellipse')
  semia = 1.75*(0.6:0.2:2.0);
  semib = 1.25*(0.6:0.2:2.0);
  nreg = length(semia);
end
  
%times of saturation
% nsat=[0.1 0.4 1 2 4 10 20 40 100];
  nsat=[0.4 1 2 4 10 20 40 100];
nnsat = length(nsat);
  
%length of each simulation
sps = 160;
greenlen = pow2(9)*sps/40;
bufsec = 1;
msftaddm = bufsec*sps;  %buffer range for later CC alignment, +1 for safety
rccmwsec = 0.5;
rccmwlen = rccmwsec*sps;  %window length for computing RCC
overshoot = rccmwlen/2; %RCC points to the center of computing window, so extra room is needed

% Twin=0.5*3600+(greenlen+msftaddm*2+overshoot*2-2)/sps; %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
% nrun = 6;

Twin=3*3600+(greenlen+msftaddm*2+overshoot*2-2)/sps; %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
nrun = 1;
  
fnsuffix1 = '';

%%%Flag to indicate if it is necessary to recalculate everything
flagrecalc = 0;
% flagrecalc = 1;

if flagrecalc
  %%%loop for noise level
  for ireg = 1: nreg
    if strcmp(srcregion,'circle')
      reg = radi;
      fprintf('reg: %f*%f; %d/%d \n',radi,radi,ireg,nreg);
    elseif strcmp(srcregion,'ellipse')
      xaxis = semia(ireg); %axis length of the same ellipse of my 4-s catalog
      yaxis = semib(ireg);
      reg = [xaxis yaxis];
      fprintf('reg: %f*%f; %d/%d \n',xaxis,yaxis,ireg,nreg);
    end
      
    %%%loop for saturation level
    parfor insat = 1: nnsat
%     for insat = 2: nnsat
      sat = nsat(insat);
      fprintf('sat: %f; %d/%d \n',sat,insat,nnsat);

      decon_synth_fn(srcregion,reg,sat,distr,Twin,tdura,distrloc,...
        physicalsize,testsrcflag,testfreqflag,pltdataflag,pltgtflag,...
        pltsrcflag,pltsrc4thflag,normflag,tempflag,ftrans,nrun);
    
    end %loop end for saturation level
    
  end %loop end for region size
    
end %if need to recalculate

% keyboard

%% load decon results
%%%Flag to indicate if to discard the detections from misaligned wins 
% flagdiscard = 0;
flagdiscard = 1;

for ireg = 1: nreg
  xaxis = semia(ireg);
  yaxis = semib(ireg);
  for insat = 1: nnsat
    sat = nsat(insat);
    
    savefile = strcat('rst_decon_synth_reg',num2str(xaxis),...
      '-',num2str(yaxis),'_nsat',num2str(sat),'_td',num2str(tdura),'.mat'); %,srcregion(1:3),'_'
    load(strcat(workpath,'/synthetics/',savefile));
    
    impgt{insat,ireg} = allsyn.impgt; %ground-truth sources
    nsrcgt(insat,ireg) = size(allsyn.impgt,1);
    impgrp{insat,ireg} = allsyn.impgrp; %grouped sources
    nsrcgrp(insat,ireg) = size(allsyn.impgrp,1);
    imp{insat,ireg} = allsyn.imp;
    % nsrck{insat,ireg} = allsyn.nsrc;
    nsrc(insat,ireg) = allsyn.nsrcsum;
    % nitk{insat,ireg} = allsyn.nit;
    nit{insat,ireg} = allsyn.nitsum;
    nitmin(insat,ireg) = min(nit{insat,ireg});
    nitmean(insat,ireg) = mean(nit{insat,ireg});
    srcampr{insat,ireg} = allsyn.srcampr;
    mprojxnn1(insat,ireg) = allsyn.mprojxnn1;
    mprojx2all(insat,ireg) = allsyn.mprojx2all;
    irccran{insat,ireg} = allsyn.irccran;  %stores start and end indices of RCC
    
    imp4th{insat,ireg} = allsyn.imp4th;
    % nsrc4thk{insat,ireg} = allsyn.nsrc4th;
    nsrc4th(insat,ireg) = allsyn.nsrc4thsum;
    srcampr4th{insat,ireg} = allsyn.srcampr4th;
    mprojxnn14th(insat,ireg) = allsyn.mprojxnn14th;
    mprojx2all4th(insat,ireg) = allsyn.mprojx2all4th;

    %load misaligned windows
    savefile2 = strcat('badwins_decon_synth_reg',num2str(xaxis),...
      '-',num2str(yaxis),'_nsat',num2str(sat),'_td',num2str(tdura),'.mat');
    load(strcat(workpath,'/synthetics/',savefile2));
    badwins{insat,ireg} = allsyn.badwins;   
    
    if flagdiscard
      %grouped sources
      [impnew,nsrcnew,lsignew] = discard_badwins_decon_synth(impgrp{insat,ireg},...
        badwins{insat,ireg},irccran{insat,ireg},sps);
      impgrp{insat,ireg} = impnew;
      nsrcgrp(insat,ireg) = nsrcnew;
      lsiggrp(insat,ireg) = lsignew;

      %3-station sources
      [impnew,nsrcnew,lsignew] = discard_badwins_decon_synth(imp{insat,ireg},...
        badwins{insat,ireg},irccran{insat,ireg},sps);
      imp{insat,ireg} = impnew;
      nsrc(insat,ireg) = nsrcnew;
      lsig(insat,ireg) = lsignew;

      %4-station sources   
      [impnew,nsrcnew,lsignew] = discard_badwins_decon_synth(imp4th{insat,ireg},...
        badwins{insat,ireg},irccran{insat,ireg},sps);
      imp4th{insat,ireg} = impnew;
      nsrc4th(insat,ireg) = nsrcnew;
      lsig4th(insat,ireg) = lsignew;
    else
      lsiggrp(insat,ireg) = Twin-(greenlen+msftaddm*2+overshoot*2-2)/sps;
      lsig(insat,ireg) = Twin-(greenlen+msftaddm*2+overshoot*2-2)/sps;  %length of the signal of synthetics
      lsig4th(insat,ireg) = Twin-(greenlen+msftaddm*2+overshoot*2-2)/sps;  %length of the signal of synthetics    
    end
       
  end
end

%%%load results from data
% savefile = 'deconv_stats4th_allbstsig.mat';
savefile = 'deconv_stats4th_no23_allbstsig.mat';
allsig = load(strcat(rstpath, '/MAPS/',savefile));

%%%load results from noise
% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
allnoi = load(strcat(rstpath, '/MAPS/',savefile));


%% summarize cumulative density plot for each region size and sat level
% %%%loop for region size
% for ireg = 1: nreg
%   disp(semia(ireg));
%   disp(semib(ireg));
%   
%   %initialize figs
%   f1 = initfig(16,8,2,4,(ireg-1)*2+1); %map loc for srcs no 2ndary srcs
%   supertit(f1.ax,'Secondary sources removed');
%   f2 = initfig(16,8,2,4,ireg*2); %map loc for srcs checked at 4th stas
%   supertit(f2.ax,'Checkd at 4th stas');
% 
%   %params of limited source region, subject to variation!
%   if strcmp(srcregion,'circle')
%     shiftor=[0.2 0.2]; %(in km) %center of the same ellipse of my 4-s catalog
%     radi=1.25; %radius
%     [xcut,ycut] = circle_chao(shiftor(1),shiftor(2),radi,0.01);
%   elseif strcmp(srcregion,'ellipse')
%     xaxis = semia(ireg);
%     yaxis = semib(ireg);
%     % xaxis=1.75; %axis length of the same ellipse of my 4-s catalog
%     % yaxis=1.25;
%     shiftor=[0.2 0.2]; %(in km) %center of the same ellipse of my 4-s catalog
%     [xcut,ycut] = ellipse_chao(shiftor(1),shiftor(2),xaxis,yaxis,0.01,45,shiftor);
%   end
%   
%   %%%loop for saturation level
%   for insat = 1: nnsat
%     
%     %%%3-sta detections
%     impi = imp{insat,ireg};
%     
%     %plot the scatter of sources in terms of rela locations
%     % xran = [-4 4];
%     % yran = [-4 4];
%     % cran = [0 lsig/sps];
%     % offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
%     % offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
%     % wt = mean(impindepst(:,[2 4 6]),2);
%     % wtmax = prctile(wt,95); %use percentile in case
%     % refscl = wt./wtmax;
%     % refscl(refscl>=1) = 1;  %force the larger amp to be plotted as the same size in case of saturation
%     % ax=f1.ax(insat);
%     % [ax] = plt_decon_imp_scatter_space(ax,impindepst,xran,yran,cran,offxran,...
%     %   offyran,sps,30,ftrans,'mean','tarvl');
%     % plot(ax,xcut,ycut,'k-','linew',2);
%     % text(ax,0.98,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
%     %       'HorizontalAlignment','right');
%     
%     %plot the cumulative density and summed amp of detections
%     density1d = density_pixel(impi(:,7),impi(:,8));
%     [imploc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
%     density1d = [imploc(:,1:2) density1d(:,3)];
%     xran = [-4 4];
%     yran = [-4 4];
%     scale = 'linear';
%     ax=f1.ax(insat);
%     hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%     dum = density1d;
%     dum(dum(:,3)>1, :) = [];
%     if strcmp(scale,'log')
%       dum(:,3) = log10(dum(:,3));
%     end
%     scatter(ax,dum(:,1),dum(:,2),10,dum(:,3),'linew',0.2);  %, 'MarkerEdgeColor', 'w')
%     dum = sortrows(density1d,3);
%     dum(dum(:,3)==1, :) = [];
%     if strcmp(scale,'log')
%       dum(:,3) = log10(dum(:,3));
%     end
%     scatter(ax,dum(:,1),dum(:,2),10,dum(:,3),'filled','MarkerEdgeColor','none');
% %     colormap(ax,flipud(colormap(ax,'kelicol')));
%     colormap(ax,'plasma');
%     c=colorbar(ax,'SouthOutside');
%     ax.CLim(2) = prctile(dum(:,3),99);
%     if insat == 1
%       if strcmp(scale,'log')
%         c.Label.String = strcat('log_{10}(# of detections)');
%       elseif strcmp(scale,'linear')
%         c.Label.String = '# of detections';
%       end
%       xlabel(ax,'E (km)','FontSize',11);
%       ylabel(ax,'N (km)','FontSize',11);
%     end
%     axis(ax,'equal');
%     axis(ax,[xran yran]);
%     ax.GridLineStyle = '--';
%     ax.XAxisLocation = 'top';
%     plot(ax,xcut,ycut,'k-','linew',2);
%     text(ax,0.98,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
%       'HorizontalAlignment','right');
%     text(ax,0.98,0.05,sprintf('%d events',size(impi,1)),'Units','normalized',...
%       'HorizontalAlignment','right','FontSize',9);
%     text(ax,0.02,0.05,sprintf('%.2f',mprojx2all{insat,ireg}),'Units',...
%       'normalized','HorizontalAlignment','left');
%     longticks(ax,2);
%     
% 
%     %%%4-sta detections
%     imp4thi = imp4th{insat,ireg};
%     
%     % %plot the scatter of sources in terms of rela locations
%     % xran = [-4 4];
%     % yran = [-4 4];
%     % cran = [0 lsig/sps];
%     % offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
%     % offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
%     % ax=f2.ax(insat);
%     % [ax] = plt_decon_imp_scatter_space(ax,impindepst,xran,yran,cran,offxran,...
%     %   offyran,sps,30,ftrans,'mean','tarvl');
%     % plot(ax,xcut,ycut,'k-','linew',2);
%     % text(ax,0.98,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
%     %       'HorizontalAlignment','right');
%     
%     %plot the cumulative density and summed amp of detections
%     density1d = density_pixel(imp4thi(:,7),imp4thi(:,8));
%     [imploc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
%     density1d = [imploc(:,1:2) density1d(:,3)];
%     xran = [-4 4];
%     yran = [-4 4];
%     ax=f2.ax(insat);
%     hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%     dum = density1d;
%     dum(dum(:,3)>1, :) = [];
%     if strcmp(scale,'log')
%       dum(:,3) = log10(dum(:,3));
%     end
%     scatter(ax,dum(:,1),dum(:,2),10,dum(:,3),'linew',0.2);  %, 'MarkerEdgeColor', 'w')
%     dum = sortrows(density1d,3);
%     dum(dum(:,3)==1, :) = [];
%     if strcmp(scale,'log')
%       dum(:,3) = log10(dum(:,3));
%     end
%     scatter(ax,dum(:,1),dum(:,2),10,dum(:,3),'filled','MarkerEdgeColor','none');
% %     colormap(ax,flipud(colormap(ax,'kelicol')));
%     colormap(ax,'plasma');
%     c=colorbar(ax,'SouthOutside');
%     ax.CLim(2) = prctile(dum(:,3),99);
%     if insat == 1
%       if strcmp(scale,'log')
%         c.Label.String = strcat('log_{10}(# of detections)');
%       elseif strcmp(scale,'linear')
%         c.Label.String = '# of detections';
%       end
%       xlabel(ax,'E (km)','FontSize',11);
%       ylabel(ax,'N (km)','FontSize',11);
%     end
%     axis(ax,'equal');
%     axis(ax,[xran yran]);
%     ax.GridLineStyle = '--';
%     ax.XAxisLocation = 'top';
%     plot(ax,xcut,ycut,'k-','linew',2);
%     text(ax,0.98,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
%       'HorizontalAlignment','right');
%     text(ax,0.98,0.05,sprintf('%d events',size(imp4thi,1)),'Units','normalized',...
%       'HorizontalAlignment','right','FontSize',9);
%     text(ax,0.02,0.05,sprintf('%.2f',mprojx2all4th{insat,ireg}),'Units','normalized',...
%       'HorizontalAlignment','left');
%     longticks(ax,2);
% 
%   
%   end %loop end for saturation level  
% end %loop end for src region size

%% num of detections VS saturation rate & region size
nsrcd = allsig.allbstsig.nsrc;
nsrc4thd = allsig.allbstsig.nsrc4th;
ngrpd = allsig.allbstsig.ngrp;
nitd = allsig.allbstsig.nitk;
nitd = cat(1, nitd{:});
% nsrcn = allnoi.allbstnoi.nsrc;
% nsrc4thn = allnoi.allbstnoi.nsrc4th;
% ngrpn = allnoi.allbstnoi.ngrp;
% nitn = allnoi.allbstnoi.nitk;
% nitn = cat(1, nitn{:});

nrow = 3; ncol = 3;
widin = 8.4;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);
pltxran = [0.08 0.96]; pltyran = [0.08 0.96]; % optimal axis location
pltxsep = 0.06; pltysep = 0.08;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

color = gradientblue(nreg);

%%%ground-truth num of srcs
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log'; ax.YScale='log';
for ireg = 1: nreg
  plot(ax,nsat,nsrcgt(:,ireg),'-',...
    'linew',1,'color',color(ireg,:),'marker','o','markersize',3,...
    'markerfacec',color(ireg,:));
end
xlabel(ax,'Saturation');
ylabel(ax,'# of ground-truth sources');
% yran = [0 1];
% ylim(ax,yran);
xticks(ax,nsat);
longticks(ax,2);
% ax.YAxis.Exponent = 4;
% axranexp(ax,6,20);
hold(ax,'off');

%%%mean num of iterations
ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log'; ax.YScale='log';
for ireg = 1: nreg
  plot(ax,nsat,nitmin(:,ireg),'-',...
    'linew',1,'color',color(ireg,:),'marker','o','markersize',3,...
    'markerfacec',color(ireg,:));
end
xlabel(ax,'Saturation');
ylabel(ax,'Min # of iterations');
yran = [2e3 4e4];
ylim(ax,yran);
yticks(ax,[0.2 0.4 1 2 4]*1e4);
xticks(ax,nsat);
longticks(ax,2);
ax.YAxis.Exponent = 4;
% axranexp(ax,6,20);
% yran = ax.YLim;
hold(ax,'off');

%%%num of grouped sources
ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log'; ax.YScale='log';
for ireg = 1: nreg
  plot(ax,nsat,nsrcgrp(:,ireg),'-',...
    'linew',1,'color',color(ireg,:),'marker','o','markersize',3,...
    'markerfacec',color(ireg,:));
end
text(ax,0.98,0.95,'Grouped','HorizontalAlignment','right','Units','normalized',...
  'fontsize',10);
xlabel(ax,'Saturation');
ylabel(ax,'# of detections');
% yran = [0 1];
ylim(ax,yran);
yticks(ax,[0.2 0.4 1 2 4]*1e4);
xticks(ax,nsat);
longticks(ax,2);
ax.YAxis.Exponent = 4;
% axranexp(ax,6,20);
hold(ax,'off');

%%%3-station sources
ax=f.ax(4); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log'; ax.YScale='log';
for ireg = 1: nreg
  plot(ax,nsat,nsrc(:,ireg),'-',...
    'linew',1,'color',color(ireg,:),'marker','o','markersize',3,...
    'markerfacec',color(ireg,:));
end
text(ax,0.98,0.95,'3-station','HorizontalAlignment','right','Units','normalized',...
  'fontsize',10);
xlabel(ax,'Saturation');
ylabel(ax,'# of detections');
% yran = [4e3 1.5e4];
% yran = [0.1 1.3];
ylim(ax,yran);
yticks(ax,[0.2 0.4 1 2 4]*1e4);
xticks(ax,nsat);
longticks(ax,2);
ax.YAxis.Exponent = 4;
% yran = ax.YLim;
% yyaxis(ax,'right');
% ylabel(ax,'# of detections per sec');
% ylim(ax,yran/lsig);
hold(ax,'off');

%%%4-station sources
ax=f.ax(5); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log'; ax.YScale='log';
% yyaxis(ax,'left');
for ireg = 1: nreg
  plot(ax,nsat,nsrc4th(:,ireg),'-',...
    'linew',1,'color',color(ireg,:),'marker','o','markersize',3,...
    'markerfacec',color(ireg,:));
end
text(ax,0.98,0.95,'4-station','HorizontalAlignment','right','Units','normalized',...
  'fontsize',10);
xlabel(ax,'Saturation');
ylabel(ax,'# of detections');
% yran = [0 1];
ylim(ax,yran);
yticks(ax,[0.2 0.4 1 2 4]*1e4);
xticks(ax,nsat);
longticks(ax,2);
ax.YAxis(1).Exponent = 4;
% yran = ax.YLim;
% yyaxis(ax,'right');
% ylabel(ax,'# of detections per sec');
% ylim(ax,yran/lsig);
hold(ax,'off');

%%%4- / 3-station ratio
ax=f.ax(6); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log';
for ireg = 1: nreg
  p(ireg) = plot(ax,nsat,nsrc4th(:,ireg)./nsrc(:,ireg),'-',...
    'linew',1,'color',color(ireg,:),'marker','o','markersize',3,...
    'markerfacec',color(ireg,:));
  % p(ireg) = plot(ax,nsat,nsrc4th(:,ireg)./nsrc(:,ireg),'-','Color',...
  %   color(ireg,:),'linew',1);
  % scatter(ax,nsat,nsrc4th(:,ireg)./nsrc(:,ireg),20,color(ireg,:),'filled',...
  %   'MarkerEdgeColor','k');
  label{ireg} = sprintf('%.1fx%.1f',2*semia(ireg),2*semib(ireg));
end
p(nreg+1) = plot(ax,ax.XLim,[sum(nsrc4thd)/sum(nsrcd) sum(nsrc4thd)/sum(nsrcd)],'k--','linew',1.5);
label{nreg+1} = 'Data';
% p(nreg+2) = plot(ax,ax.XLim,[sum(nsrc4thn)/sum(nsrcn) sum(nsrc4thn)/sum(nsrcn)],'r--','linew',1.5);
% label{nreg+2} = 'Noise';
lgd=legend(f.ax(4),p,label,'NumColumns',2,'Location','south','fontsize',6);
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
lgdtit = 'Region size (km)';
title(lgd,lgdtit,'fontsize',7);
% title(ax,'Secondary sources removed');
xlabel(ax,'Saturation');
ylabel(ax,'Ratio of 4- to 3-station');
yran = [0.3 1];
ylim(ax,yran);
xticks(ax,nsat);
longticks(ax,2);
hold(ax,'off');

%%%3-station sources, detection rate 
ax=f.ax(7); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log';
for ireg = 1: nreg
  plot(ax,nsat/0.25,nsrc(:,ireg)./lsig(:,ireg),'-',...
    'linew',1,'color',color(ireg,:),'marker','o','markersize',3,...
    'markerfacec',color(ireg,:));
end
plot(ax,ax.XLim,[sum(nsrcd)/tlensum sum(nsrcd)/tlensum],'k--','linew',1.5);
% plot(ax,ax.XLim,[sum(nsrcn)/tlensum sum(nsrcn)/tlensum],'r--','linew',1.5);
text(ax,0.98,0.95,'3-station','HorizontalAlignment','right','Units','normalized',...
  'fontsize',10);
xlabel(ax,'# of ground-truth sources per sec');
ylabel(ax,'# of detections per sec');
yran = [0.1 1.3];
ylim(ax,yran);
xticks(ax,nsat/0.25);
xlim(ax,[nsat(1) nsat(end)]/0.25);
longticks(ax,2);
hold(ax,'off');

%%%4-station sources, detection rate
ax=f.ax(8); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log';
for ireg = 1: nreg
  plot(ax,nsat/0.25,nsrc4th(:,ireg)./lsig4th(:,ireg),'-',...
    'linew',1,'color',color(ireg,:),'marker','o','markersize',3,...
    'markerfacec',color(ireg,:));
end
plot(ax,ax.XLim,[sum(nsrc4thd)/tlensum sum(nsrc4thd)/tlensum],'k--','linew',1.5);
% plot(ax,ax.XLim,[sum(nsrc4thn)/tlensum sum(nsrc4thn)/tlensum],'r--','linew',1.5);
text(ax,0.98,0.95,'4-station','HorizontalAlignment','right','Units','normalized',...
  'fontsize',10);
xlabel(ax,'# of ground-truth sources per sec');
ylabel(ax,'# of detections per sec');
% yran = [0 1];
ylim(ax,yran);
xticks(ax,nsat/0.25);
xlim(ax,[nsat(1) nsat(end)]/0.25);
longticks(ax,2);
hold(ax,'off');

%%%min num of iteration, rate, rough, as containing iterations within misaligned
%%%windows
ax=f.ax(9); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log';
for ireg = 1: nreg
  plot(ax,nsat/0.25,nitmin(:,ireg)/(3*3600),'-',...
    'linew',1,'color',color(ireg,:),'marker','o','markersize',3,...
    'markerfacec',color(ireg,:));
end
plot(ax,ax.XLim,[min(sum(nitd,1))/tlensum min(sum(nitd,1))/tlensum],'k--','linew',1.5);
% plot(ax,ax.XLim,[min(sum(nitn,1))/tlensum min(sum(nitn,1))/tlensum],'r--','linew',1.5);
xlabel(ax,'# of ground-truth sources per sec');
ylabel(ax,'Min # of iterations per sec');
% yran = [0.2 4]*1e4;
% ylim(ax,yran);
% yticks(ax,(0.5:1:4)*1e4);
xticks(ax,nsat);
longticks(ax,2);
% ax.YAxis.Exponent = 4;
% axranexp(ax,6,20);
% yran = ax.YLim;
hold(ax,'off');


fname = strcat('num_det_syn',fnsuffix1,'.pdf');
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));
keyboard


%% summarize consecutive dist along min-scatter VS saturation rate & region size
% f = initfig(8,4.5,1,2); %initialize fig
% color = gradientblue(nreg);
% ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% for ireg = 1: nreg
%   % p(ireg) = plot(ax,log10(nsat),mprojx2nn1(:,ireg),'-o','markersize',4,'color',color(ireg,:));
% %   p(ireg) = plot(ax,log10(nsat),mprojx22all(:,ireg),'-o','markersize',4,...
% %     'color',color(ireg,:),'filled');
%   p(ireg) = plot(ax,log10(nsat),mprojx2all(:,ireg),'-','Color',color(ireg,:),'linew',1);
%   scatter(ax,log10(nsat),mprojx2all(:,ireg),nsrcm(:,ireg)/50,color(ireg,:),...
%     'filled','MarkerEdgeColor','k');
%   label{ireg} = sprintf('a=%.1f, b=%.1f',semia(ireg),semib(ireg));
% end
% p(nreg+1) = plot(ax,log10(nsat),0.50*ones(nnsat,1),'k--','linew',1);  %this is from data
% label{nreg+1} = 'Data';
% % legend(ax,p,label,'NumColumns',2,'Location','south');
% % title(ax,'Secondary sources removed');
% xlabel(ax,'log_{10}(Saturation)');
% ylabel(ax,'Distance (km)');
% yran = [0 1];
% ylim(ax,yran);
% longticks(ax,2);
% hold(ax,'off');
% 
% ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% for ireg = 1: nreg
%   % plot(ax,log10(nsat),mprojx3nn1(:,ireg),'-o','markersize',4,'color',color(ireg,:));
% %   plot(ax,log10(nsat),mprojx32all(:,ireg),'-o','markersize',4,...
% %     'color',color(ireg,:),'filled');
%   plot(ax,log10(nsat),mprojx2all4th(:,ireg),'-','Color',color(ireg,:),'linew',1);
%   scatter(ax,log10(nsat),mprojx2all4th(:,ireg),nsrc4thm(:,ireg)/50,color(ireg,:),...
%     'filled','MarkerEdgeColor','k');
% end
% plot(ax,log10(nsat),0.45*ones(nnsat,1),'k--','linew',1);  %this is from data
% % title(ax,'Checkd at 4th stas');
% ylim(ax,yran);
% longticks(ax,2);
% hold(ax,'off');

% ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% for ireg = 1: nreg
%   % plot(ax,log10(nsat),mprojx1nn1(:,ireg),'-o','markersize',4,'color',color(ireg,:));
%   plot(ax,log10(nsat),mprojx12all(:,ireg),'-o','markersize',4,'color',color(ireg,:));
% end
% title(ax,'Grouped Total');
% ylim(ax,yran);

% stit = supertit(f.ax,'Med. dist. along min-error direc from each to all others w/i 2 s');
% movev(stit,0.4);

% keyboard

%% Summary of amplitude ratio 
% %%%%%%%%% you can choose to plot the actual histogram for each sat and noise
%%%loop for region size
% for ireg = 1: nreg
%   disp(semia(ireg));
%   disp(semib(ireg));
%   %%%loop for saturation level
%   for insat = 1: nnsat
%     disp(nsat(insat));
%     f3 = initfig(16,8,2,4); %plot histograms of source amp
%     supertit(f3.ax,'Secondary sources removed & Checkd at 4th stas');
%     impindepst = imp{insat,ireg};
%     srcampr = srcamprall{insat,ireg};
%     f3.ax(1:3) = plt_deconpk_rat_comb(f3.ax(1:3),srcampr,impindepst,'k','hist');  
%     impindepst = imp4th{insat,ireg};
%     srcampr4th = srcamprall4th{insat,ireg};
%     f3.ax(5:end) = plt_deconpk_rat_comb4th(f3.ax(5:end),srcampr4th,impindepst,'k','hist');  
%   end %loop end for saturation level 
% end %loop end for src region size
% %%%%%%%%% you can choose to plot the actual histogram for each sat and noise

%%%%%%%%% or only plot the median & mad for each sat and noise
for ireg = 1: nreg
  label{ireg} = sprintf('%.1fx%.1f',2*semia(ireg),2*semib(ireg));
end
label{nreg+1} = 'Data';
srcampralld = allsig.allbstsig.srcamprall;  %using reference from real data
for i = 1: size(srcampralld,2)
  mamprd(i,1) = median(log10(srcampralld(:,i)));
  madamprd(i,1) =  mad(log10(srcampralld(:,i)),1);
  stdamprd(i,1) = std(log10(srcampralld(:,i)));  
end
ref{1} = mamprd;
ref{2} = stdamprd;
for ireg = 1: nreg
  for insat = 1: nnsat
    mampr{insat,ireg} = median(log10(srcampr{insat,ireg}), 1);
    madampr{insat,ireg} = mad(log10(srcampr{insat,ireg}), 1, 1);
    stdampr{insat,ireg} = std(log10(srcampr{insat,ireg}), 1);
  end
end
f = initfig(8.4,6,2,3); %initialize fig
% orient(f.fig,'landscape');
pltxran = [0.08 0.98]; pltyran = [0.08 0.98]; % optimal axis location
optaxpos(f,2,3,pltxran,pltyran,0.03,0.08);
% [f,lgd]=plt_deconpk_rat_stat(f,nsat,label,mampr,stdampr,ref,5*log10(nsrc));
[f,lgd]=plt_deconpk_rat_stat(f,nsat,label,mampr,stdampr,ref);
lgdtit = 'Region size (km)';
title(lgd(1),lgdtit,'fontsize',7);
% title(lgd(2),lgdtit,'fontsize',7);
% stit = supertit(f.ax,'Secondary sources removed');
% movev(stit,0.3);
ylim(f.ax(4:end),log10([1.2 2.4]));
% xlim(f.ax(:),[-0.5 2]);
fname = strcat('amprat_syn',fnsuffix1,'.pdf');
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));
% keyboard

%%
srcampr4thalld = allsig.allbstsig.srcampr4thall;  %using reference from real data
for i = 1: size(srcampr4thalld,2) 
  mampr4thd(i,1) = median(log10(srcampr4thalld(:,i)));
  madampr4thd(i,1) =  mad(log10(srcampr4thalld(:,i)),1);
  stdampr4thd(i,1) = std(log10(srcampr4thalld(:,i)));  
end
ref4th{1} = mampr4thd;
ref4th{2} = stdampr4thd;
for ireg = 1: nreg
  for insat = 1: nnsat
    mampr4th{insat,ireg} = median(log10(srcampr4th{insat,ireg}), 1);
    madampr4th{insat,ireg} = mad(log10(srcampr4th{insat,ireg}), 1, 1);
    stdampr4th{insat,ireg} = std(log10(srcampr4th{insat,ireg}), 1);
  end
end
f = initfig(8.4,6,2,4); %initialize fig
% orient(f.fig,'landscape');
pltxran = [0.05 0.98]; pltyran = [0.08 0.98]; % optimal axis location
optaxpos(f,2,4,pltxran,pltyran,0.02,0.08);
% [f,lgd]=plt_deconpk_rat_stat(f,nsat,label,mampr4th,stdampr4th,ref4th,5*log10(nsrc4th));
[f,lgd]=plt_deconpk_rat_stat(f,nsat,label,mampr4th,stdampr4th,ref4th);
title(lgd(1),lgdtit,'fontsize',7);
% title(lgd(2),lgdtit,'fontsize',7);
% stit = supertit(f.ax,'Checkd at 4th stas');
% movev(stit,0.3);
ylim(f.ax(5:end),log10([1.2 2.4]));
% xlim(f.ax(:),[-0.5 2]);

fname = strcat('amprat_syn4th',fnsuffix1,'.pdf');
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));
% keyboard
%%%%%%%%% or only plot the median & mad for each sat and noise

%% Best alignment for the testing window
%
% %compute the median absolute amplitude and envelope of the same moving window
% %for the moving window at the same station, sensable to use median
% [ir,ramp1,renv1] = Runningampenv(sigpnsta(:,1),mwlen,mwlen-1,'median');
% [~,ramp2,renv2] = Runningampenv(sigpnsta(:,2),mwlen,mwlen-1,'median');
% [~,ramp3,renv3] = Runningampenv(sigpnsta(:,3),mwlen,mwlen-1,'median');
% %looks like using the amplitude and envelope are pretty similar
% %maybe better to mean, otherwise stations PGC will be underrepresented due to a larger amp
% %variation
% ramp = mean([ramp1 ramp2 ramp3], 2);  % use mean or median??
% renv = mean([renv1 renv2 renv3], 2);
%
% figure
% ax = gca; hold(ax,'on');
% scatter(ax,renv,rcc,10,ir,'filled');
% scatter(ax,median(renv),median(rcc),8,'ko','linew',1);
% [xcnt,ycnt,y1sig] = ranybinx(renv,rcc,'median',10);
% errorbar(ax,xcnt,ycnt,-y1sig,y1sig,'vertical','o','markersize',3,'color',...
%   'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor','k','CapSize',4);
% [rho1,~] = corr(renv,rcc,'Type','Spearman');
% [rho2,~] = corr(renv,rcc,'Type','Kendall');
% text(ax,0.98,0.55,sprintf('S: %.2f',rho1),'unit','normalized',...
%   'HorizontalAlignment','right','fontsize',12);
% text(ax,0.98,0.50,sprintf('K: %.2f',rho2),'unit','normalized',...
%   'HorizontalAlignment','right','fontsize',12);
% colormap(ax,'jet');
% caxis(ax,[0 size(sigpnsta,1)]);
% c=colorbar(ax,'south');
% c.Label.String = sprintf('Samples at %d Hz',sps);
% ylim(ax,[-1 1]);
% xlabel(ax,'Running envelope');
% ylabel(ax,'Running CC');
% hold(ax,'off');
%
% %get the sources that arrived within the moving window, and estimate the source region size
% wins = movingwins(1,size(sigpnsta,1)-1,mwlen,mwlen-1);
% nwin = size(wins,1);
% srcregsz = -10*ones(nwin,1);
% nsrci = -10*ones(nwin,1);
% dist = -10*ones(nwin,1);
% distmed = -10*ones(nwin,1);
% for i = 1: nwin
%   win = wins(i,:);  % start and end indices of the moving window
%   srci = srcpair(srcpair(:,1)>=win(1) & srcpair(:,1)<win(2), :);
%   if ~isempty(srci)
%     nsrci(i) = size(srci,1);
%     %how to estimate the 'size' of source region
%     srcregsz(i) = sqrt(polyarea(srci(:,7), srci(:,8)));
%     %find the centroid of the sources and estimate its distance to best alignment
%     srccntrd = [median(srci(:,7)) median(srci(:,8))];
%     dist(i) = sqrt((srccntrd(1)-off1i(2))^2 + (srccntrd(2)-off1i(3))^2);
%     %find the distance from each source to best alignment
%     distsrc = sqrt((srci(:,7)-off1i(2)).^2 + (srci(:,8)-off1i(3)).^2);
%     %the median of the distance
%     distmed(i) = median(distsrc);
% %     srcregsz(i) = sqrt(max(abs(srci(:,7)))^2 + max(abs(srci(:,8)))^2);
% %     srcregsz(i) = sqrt(median(abs(srci(:,7)))^2 + median(abs(srci(:,8)))^2);
%   end
% end
%
% figure
% subplot(131)
% ax = gca; hold(ax,'on');
% ind = find(srcregsz==-10);
% scatter(ax,renv(ind),rcc(ind),6,[.5 .5 .5]);
% ind = find(srcregsz~=-10);
% scatter(ax,renv(ind),rcc(ind),10,srcregsz(ind),'filled');
% colormap(ax,'jet');
% % caxis(ax,[0 size(sigpnsta,1)]);
% c=colorbar(ax,'south');
% c.Label.String = 'Size of source region (samples)';
% % c.Label.String = 'Median size of source region (samples)';
% ylim(ax,[-1 1]);
% xlabel(ax,'Running envelope');
% ylabel(ax,'Running CC');
% hold(ax,'off');
%
% subplot(132)
% ax = gca; hold(ax,'on');
% ind = find(dist==-10);
% scatter(ax,renv(ind),rcc(ind),6,[.5 .5 .5]);
% ind = find(dist~=-10);
% scatter(ax,renv(ind),rcc(ind),10,distmed(ind),'filled');
% colormap(ax,'jet');
% % caxis(ax,[0 size(sigpnsta,1)]);
% c=colorbar(ax,'south');
% c.Label.String = 'Median distance to best alignment (samples)';
% % c.Label.String = 'Distance of centroid to best alignment (samples)';
% % c.Label.String = 'Median size of source region (samples)';
% ylim(ax,[-1 1]);
% xlabel(ax,'Running envelope');
% ylabel(ax,'Running CC');
% hold(ax,'off');
%
% subplot(133)
% ax = gca; hold(ax,'on');
% ind = find(dist==-10);
% scatter(ax,renv(ind),rcc(ind),6,[.5 .5 .5]);
% ind = find(dist~=-10);
% scatter(ax,renv(ind),rcc(ind),10,nsrci(ind),'filled');
% colormap(ax,'jet');
% % caxis(ax,[0 size(sigpnsta,1)]);
% c=colorbar(ax,'south');
% c.Label.String = 'Number of sources';
% ylim(ax,[-1 1]);
% xlabel(ax,'Running envelope');
% ylabel(ax,'Running CC');
% hold(ax,'off');
%
% figure
% yyaxis('left');
% plot(sigpnsta(:,1),'r-'); hold on
% plot(sigpnsta(:,2),'b-');
% plot(sigpnsta(:,3),'k-'); %,'linew',0.5
% axranexp(gca,6,10);
% xlabel(sprintf('Samples at %d Hz',sps));
% ylabel('Amplitude');
% yyaxis('right');
% % plot(ircc,rcc,'o','color',[.7 .7 .7],'markersize',2);
% scatter(ircc,rcc,2,renv);  % scale it to the amplitude range
% colormap('jet');
% caxis([min(renv) max(renv)]);
% ylabel('Running CC','FontSize',10);
% ylim([-1.1 1.1]);
%
% %distribution of synthetic sources, ground truth
% spsscale = sps/40;
% loff_max = 4*spsscale;  % maximum allowed shift is 4 samples at 40 Hz, it is also true from the synthetic distribution
% offxran = [-loff_max loff_max];
% offyran = [-loff_max loff_max];
% lsig = size(sigpnsta,1);
% cran = [0 lsig];
% [f] = plt_decon_imp_scatter(srcpair,offxran,offyran,cran,sps,1,'mean');
% title(gca,sprintf('Synthetic sources, using data of %d Hz',sps));
%
%
% %%
% figure;
% colorn = {'r','b','k'};
% for i = 1: nsta
%   subplot(5,1,i)
%   ax = gca;
%   hold(ax,'on');
%   ax.Box = 'on';
%   grid(ax,'on');
%   indtemp = find(sigdecon(:,i)>0);  % independent impulses from other stations
%   stem(ax,indtemp,sigdecon(indtemp,i),'color',[.8 .8 .8],'MarkerSize',4);
%   text(ax,0.9,0.9,num2str(length(indtemp)),'unit','normalized','HorizontalAlignment','right',...
%     'FontSize',10,'color',[.8 .8 .8]);
%   indtemp = find(impindep(:,(i-1)*2+1)>0);  % that can be individually paired with 3-sta among independent
%   stem(ax,impindep(indtemp,(i-1)*2+1),impindep(indtemp,(i-1)*2+2),'color',colorn{i},...
%     'MarkerSize',4);
%   text(ax,0.9,0.75,num2str(length(indtemp)),'unit','normalized','HorizontalAlignment','right',...
%     'FontSize',10,'color',colorn{i});
%   text(ax,0.05,0.9,strcat(stas(i,:)),'unit','normalized','HorizontalAlignment','left',...
%     'FontSize',12);
%   xlim(ax,[0 lsig]); hold(ax,'off');
%   if i == 1
%     title('Triplet impulses from grouped independent deconvolution');
%   end
% end
% subplot(5,1,4)
% ax = gca;
% yyaxis(ax,'left');
% hold(ax,'on');
% ax.Box = 'on';
% grid(ax,'on');
% stem(ax,impindep(:,1),impindep(:,2), 'r-','MarkerSize',4);
% stem(ax,impindep(:,3),impindep(:,4), 'b-','MarkerSize',4);
% stem(ax,impindep(:,5),impindep(:,6), 'k-','MarkerSize',4);
% text(ax,0.05,0.9,sprintf('Number of triplets: %d', length(indpair)),'fontsize',10,...
%   'Units','normalized');
% yyaxis(ax,'right');
% plot(ircc,rcc,'o','color',[.5 .5 .5],'markersize',2);
% ylim(ax,[-1 1]);
% xlim(ax,[0 lsig]);
% hold(ax,'off');
%
% %the synthetic sources, ground truth
% subplot(5,1,5)
% ax = gca;
% yyaxis(ax,'left');
% hold(ax,'on');
% ax.Box = 'on';
% stem(ax,srcpair(:,1),srcpair(:,2), 'r-','MarkerSize',4);
% stem(ax,srcpair(:,3),srcpair(:,4), 'b-','MarkerSize',4);
% stem(ax,srcpair(:,5),srcpair(:,6), 'k-','MarkerSize',4);
% text(ax,0.05,0.9,sprintf('Number of triplets: %d',size(srcpair,1)),'fontsize',10,...
%   'Units','normalized');
% yyaxis(ax,'right');
% plot(ircc,rcc,'o','color',[.5 .5 .5],'markersize',2);
% ylim(ax,[-1 1]);
% xlim(ax,[0 lsig]);
% title(ax,'Triplet impulses from true synthetic sources');
% hold(ax,'off');
