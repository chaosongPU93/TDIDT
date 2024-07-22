% decon_synth_onespot.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --This script now becomes the driving script to run the deconvolution to
% synthetic seismograms that are generated from the same spot but with
% different added noise levels and different saturation levels.
% Parameters related to which synthetics are read and deconvolved, and those
% related to the settings of deconvolution are defined here.
% --Moreover, this script analyzes the statistics of the resulting deconvolved
% catalogs. These statistics are a bit out-of-date given our new understanding
% to data. 
% --See also 'decon_synth.m' for decon to synthetics generated from different 
% region sizes and saturation levels with no noise.
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/05/05
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

%The ratio of elsewhere events to local events.  Zero for Discrete Ide.
fracelsew=0; %0.25 %0.5; %0.6;

%%%specify regime for transformation from time offset to map location
% ftrans = 'interpArmb';
% ftrans = 'interpArmbreloc';
ftrans = 'interpchao';

timetype = 'tori';

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

%different percent of noise
perctrial = 0.1*(0:2:16)';
ntrial = length(perctrial);

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

fnsuffix1 = '_onespot';

%%%Flag to indicate if it is necessary to recalculate everything
flagrecalc = 0;
% flagrecalc = 1;

if flagrecalc
  %%%loop for noise level
  for iperc = 3: ntrial
    perc = perctrial(iperc);
    fprintf('noi: %f; %d/%d \n',perc,iperc,ntrial);
    
    %%%loop for saturation level
    parfor insat = 2: nnsat
%     for insat = 2: nnsat
      sat = nsat(insat);
      fprintf('sat: %f; %d/%d \n',sat,insat,nnsat);
      
      %%%call function to carry out deconvolution to all 'nrun' synthetics with
      %%%the SAME sat and noise fraction
      decon_synth_onespot_fn(perc,sat,distr,Twin,tdura,testfreqflag,...
        pltdataflag,pltgtflag,pltsrcflag,pltsrc4thflag,normflag,tempflag,ftrans,nrun);
      
    end %loop end for saturation level
    
  end %loop end for percent of noise
    
end %if need to recalculate

% keyboard

%% load decon results
for iperc = 1: ntrial
  perc = perctrial(iperc);
  
  for insat = 1: nnsat
    sat = nsat(insat);
    
    savefile = strcat('rst_decon_synth_onespot','_noi',num2str(perc),'_nsat',...
      num2str(sat),'_td',num2str(tdura),'.mat');
    load(strcat(workpath,'/synthetics/',savefile));
    
    nsrcgt(insat,iperc) = size(allsyn.impgt,1);
    nsrcgrp(insat,iperc) = size(allsyn.impgrp,1);
    impk{insat,iperc} = allsyn.imp;
    nsrck{insat,iperc} = allsyn.nsrc;
    nsrc(insat,iperc) = allsyn.nsrcsum;
    nitk{insat,iperc} = allsyn.nit;
    nit{insat,iperc} = allsyn.nitsum;
    nitmin(insat,iperc) = min(nit{insat,iperc});
    nitmean(insat,iperc) = mean(nit{insat,iperc});
    srcamprk{insat,iperc} = allsyn.srcampr;
    mprojxnn1(insat,iperc) = allsyn.mprojxnn1;
    mprojx2all(insat,iperc) = allsyn.mprojx2all;
    
    imp4thk{insat,iperc} = allsyn.imp4th;
    nsrc4thk{insat,iperc} = allsyn.nsrc4th;
    nsrc4th(insat,iperc) = allsyn.nsrc4thsum;
    srcampr4thk{insat,iperc} = allsyn.srcampr4th;
    mprojxnn14th(insat,iperc) = allsyn.mprojxnn14th;
    mprojx2all4th(insat,iperc) = allsyn.mprojx2all4th;
    
  end
end
lsig = Twin-(greenlen+msftaddm*2+overshoot*2-2)/sps;  %length of the signal of synthetics

%%%load results from data
% savefile = 'deconv_stats4th_allbstsig.mat';
savefile = 'deconv_stats4th_no23_allbstsig.mat';
allsig = load(strcat(rstpath, '/MAPS/',savefile));

%%%load results from noise
% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
allnoi = load(strcat(rstpath, '/MAPS/',savefile));

%% summarize cumulative density plot for each noise and sat level
% %%%loop for noise level
% for iperc = 1: ntrial
%   perc = perctrial(iperc);
%   disp(perc);
  
%   %initialize figs
%   f1 = initfig(8.4,5,2,4,(iperc-1)*2+1); %map loc for srcs no 2ndary srcs
%   supertit(f1.ax,sprintf('noise: %.1f, Secondary sources removed',perc));
%   % f2 = initfig(16,8,2,4,iperc*2); %map loc for srcs checked at 4th stas
%   % supertit(f2.ax,sprintf('noise: %.1f, Checkd at 4th stas',perc));
  
%   %%%loop for saturation level
%   for insat = 1: nnsat
    
%     impi = impk{insat,iperc};
    
%     %plot the scatter of sources in terms of rela locations
%     % xran = [-4 4];
%     % yran = [-4 4];
%     % cran = [0 lsig/sps];
%     % offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
%     % offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
%     % ax=f1.ax(insat);
%     % [ax] = plt_decon_imp_scatter_space(ax,impindepst,xran,yran,cran,offxran,...
%     %   offyran,sps,30,ftrans,'mean','tarvl');
%     % plot(ax,xcut,ycut,'k-','linew',2);
%     % text(ax,0.98,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
%     %       'HorizontalAlignment','right');
%     % text(ax,0.02,0.05,sprintf('%.2f',mprojx22all(insat,iperc)),'Units','normalized',...
%     % 'HorizontalAlignment','left');
    
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
%     text(ax,0.98,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
%       'HorizontalAlignment','right');
%     text(ax,0.98,0.05,sprintf('%d events',size(impi,1)),'Units','normalized',...
%       'HorizontalAlignment','right','FontSize',9);
%     text(ax,0.02,0.05,sprintf('%.2f',mprojx2all(insat,iperc)),'Units',...
%       'normalized','HorizontalAlignment','left');
%     text(ax,0.98,0.15,sprintf('%d; %.2f',...
%       sum(impi(:,7)==2&impi(:,8)==2),...
%       sum(impi(:,7)==2&impi(:,8)==2)/size(impi,1)),...
%       'Units','normalized','HorizontalAlignment','right','FontSize',9);
%     longticks(ax,2);

    
%     %%%4-sta detections
%     imp4thi = imp4thk{insat,iperc};
    
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
    
% %     %plot the cumulative density and summed amp of detections
% %     density1d = density_pixel(imp4thi(:,7),imp4thi(:,8));
% %     [imploc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
% %     density1d = [imploc(:,1:2) density1d(:,3)];
% %     xran = [-4 4];
% %     yran = [-4 4];
% %     scale = 'linear';
% %     ax=f1.ax(insat);
% %     hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
% %     dum = density1d;
% %     dum(dum(:,3)>1, :) = [];
% %     if strcmp(scale,'log')
% %       dum(:,3) = log10(dum(:,3));
% %     end
% %     scatter(ax,dum(:,1),dum(:,2),10,dum(:,3),'linew',0.2);  %, 'MarkerEdgeColor', 'w')
% %     dum = sortrows(density1d,3);
% %     dum(dum(:,3)==1, :) = [];
% %     if strcmp(scale,'log')
% %       dum(:,3) = log10(dum(:,3));
% %     end
% %     scatter(ax,dum(:,1),dum(:,2),10,dum(:,3),'filled','MarkerEdgeColor','none');
% % %     colormap(ax,flipud(colormap(ax,'kelicol')));
% %     colormap(ax,'plasma');
% %     c=colorbar(ax,'SouthOutside');
% %     ax.CLim(2) = prctile(dum(:,3),99);
% %     if insat == 1
% %       if strcmp(scale,'log')
% %         c.Label.String = strcat('log_{10}(# of detections)');
% %       elseif strcmp(scale,'linear')
% %         c.Label.String = '# of detections';
% %       end
% %       xlabel(ax,'E (km)','FontSize',11);
% %       ylabel(ax,'N (km)','FontSize',11);
% %     end
% %     axis(ax,'equal');
% %     axis(ax,[xran yran]);
% %     ax.GridLineStyle = '--';
% %     ax.XAxisLocation = 'top';
% %     text(ax,0.98,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
% %       'HorizontalAlignment','right');
% %     text(ax,0.98,0.05,sprintf('%d events',size(imp4thi,1)),'Units','normalized',...
% %       'HorizontalAlignment','right','FontSize',9);
% %     text(ax,0.02,0.05,sprintf('%.2f',median(mprojx2all4th(insat,iperc))),'Units',...
% %       'normalized','HorizontalAlignment','left');
% %     text(ax,0.98,0.15,sprintf('%d; %.2f',...
% %       sum(imp4thi(:,7)==2&imp4thi(:,8)==2),...
% %       sum(imp4thi(:,7)==2&imp4thi(:,8)==2)/size(imp4thi,1)),...
% %       'Units','normalized','HorizontalAlignment','right','FontSize',9);
% %     longticks(ax,2);
    
%   end %loop end for saturation level

%   fname = strcat('den_map_syn','_noi',num2str(perc),'.pdf');
%   print(f1.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));

% end %loop end for noise level

% keyboard

%% num of detections VS saturation rate & region size
nsrcd = allsig.allbstsig.nsrc;
nsrc4thd = allsig.allbstsig.nsrc4th;
nsrcn = allnoi.allbstnoi.nsrc;
nsrc4thn = allnoi.allbstnoi.nsrc4th;

nrow = 2; ncol = 3;
widin = 8.4;  % maximum width allowed is 8.5 inches
htin = 6;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);
pltxran = [0.08 0.96]; pltyran = [0.08 0.96]; % optimal axis location
pltxsep = 0.06; pltysep = 0.08;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

color = gradientblue(ntrial);

ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log'; ax.YScale='log';
for iperc = 1: ntrial
  plot(ax,nsat,nsrcgt(:,iperc),'-o','Color',color(iperc,:),...
    'markersize',4,'MarkerFaceColor',color(iperc,:),...
    'MarkerEdgeColor',color(iperc,:),'LineWidth',1);
end
xlabel(ax,'Saturation');
ylabel(ax,'# of ground-truth sources');
% yran = [0 1];
% ylim(ax,yran);
xticks(ax,nsat);
longticks(ax,2);
% ax.YAxis.Exponent = 4;
hold(ax,'off');

ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log'; ax.YScale='log';
for iperc = 1: ntrial
  plot(ax,nsat,nitmean(:,iperc),'-o','Color',color(iperc,:),...
    'markersize',4,'MarkerFaceColor',color(iperc,:),...
    'MarkerEdgeColor',color(iperc,:),'LineWidth',1);
end
xlabel(ax,'Saturation');
ylabel(ax,'Average # of iterations');
% yran = [0 1];
% ylim(ax,yran);
xticks(ax,nsat);
longticks(ax,2);
ax.YAxis.Exponent = 4;
axranexp(ax,6,20);
hold(ax,'off');

ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log'; ax.YScale='log';
for iperc = 1: ntrial
  plot(ax,nsat,nsrcgrp(:,iperc),'-o','Color',color(iperc,:),...
    'markersize',4,'MarkerFaceColor',color(iperc,:),...
    'MarkerEdgeColor',color(iperc,:),'LineWidth',1);
end
xlabel(ax,'Saturation');
ylabel(ax,'# of grouped sources');
% yran = [0 1];
% ylim(ax,yran);
xticks(ax,nsat);
longticks(ax,2);
ax.YAxis.Exponent = 4;
axranexp(ax,6,20);
hold(ax,'off');

ax=f.ax(4); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log';
for iperc = 1: ntrial
  plot(ax,nsat,nsrc(:,iperc)/lsig,'-o','Color',color(iperc,:),...
    'markersize',4,'MarkerFaceColor',color(iperc,:),...
    'MarkerEdgeColor',color(iperc,:),'LineWidth',1);
end
plot(ax,ax.XLim,[sum(nsrcd)/tlensum sum(nsrcd)/tlensum],'k--','linew',1.5);
plot(ax,ax.XLim,[sum(nsrcn)/tlensum sum(nsrcn)/tlensum],'r--','linew',1.5);
text(ax,0.98,0.95,'3-station','HorizontalAlignment','right','Units','normalized',...
  'fontsize',10);
% title(ax,'Secondary sources removed');
xlabel(ax,'Saturation');
ylabel(ax,'# of detections per sec');
% yran = [5e3 1.2e4];
yran = [0.1 1.3];
ylim(ax,yran);
% yran = ax.YLim;
xticks(ax,nsat);
longticks(ax,2);
% ax.YAxis.Exponent = 4;
hold(ax,'off');

ax=f.ax(5); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log';
for iperc = 1: ntrial
  plot(ax,nsat,nsrc4th(:,iperc)/lsig,'-o','Color',color(iperc,:),...
    'markersize',4,'MarkerFaceColor',color(iperc,:),...
    'MarkerEdgeColor',color(iperc,:),'LineWidth',1);
end
plot(ax,ax.XLim,[sum(nsrc4thd)/tlensum sum(nsrc4thd)/tlensum],'k--','linew',1.5);
plot(ax,ax.XLim,[sum(nsrc4thn)/tlensum sum(nsrc4thn)/tlensum],'r--','linew',1.5);
text(ax,0.98,0.95,'4-station','HorizontalAlignment','right','Units','normalized',...
  'fontsize',10);
% title(ax,'Secondary sources removed');
xlabel(ax,'Saturation');
ylabel(ax,'# of detections per sec');
% yran = [0 1];
ylim(ax,yran);
xticks(ax,nsat);
longticks(ax,2);
% ax.YAxis.Exponent = 4;
hold(ax,'off');

ax=f.ax(6); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log';
for iperc = 1: ntrial
  p(iperc) = plot(ax,nsat,nsrc4th(:,iperc)./nsrc(:,iperc),'-o','Color',color(iperc,:),...
    'markersize',4,'MarkerFaceColor',color(iperc,:),...
    'MarkerEdgeColor',color(iperc,:),'LineWidth',1);
  label{iperc} = sprintf('noi=%.1f',perctrial(iperc));
end
p(ntrial+1) = plot(ax,ax.XLim,[sum(nsrc4thd)/sum(nsrcd) sum(nsrc4thd)/sum(nsrcd)],'k--','linew',1.5);
label{ntrial+1} = 'Data';
p(ntrial+2) = plot(ax,ax.XLim,[sum(nsrc4thn)/sum(nsrcn) sum(nsrc4thn)/sum(nsrcn)],'r--','linew',1.5);
label{ntrial+2} = 'Noise';
lgd=legend(f.ax(4),p,label,'NumColumns',2,'Location','best','fontsize',6);
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
lgdtit = 'Noise level';
title(lgd,lgdtit,'fontsize',7);
% title(ax,'Secondary sources removed');
xlabel(ax,'Saturation');
ylabel(ax,'Ratio of 4- to 3-station');
yran = [0.3 1];
ylim(ax,yran);
xticks(ax,nsat);
longticks(ax,2);
hold(ax,'off');

fname = strcat('num_det_syn',fnsuffix1,'.pdf');
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));
keyboard


%% summarize consecutive dist along min-scatter VS saturation rate & region size
% f = initfig(8,4.5,1,2); %initialize fig
% color = gradientblue(ntrial);
% ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% for iperc = 1: ntrial
% %   p(iperc) = plot(ax,log10(nsat),mprojx22all(:,iperc),'-o','markersize',4,'color',color(iperc,:));
%   p(iperc) = plot(ax,log10(nsat),mprojx2all(:,iperc),'-','Color',color(iperc,:),'linew',1);
%   scatter(ax,log10(nsat),mprojx2all(:,iperc),nsrc(:,iperc)/50,color(iperc,:),...
%     'filled','MarkerEdgeColor','k');
%   label{iperc} = sprintf('Noise=%.1f',perctrial(iperc));
% end
% p(ntrial+1) = plot(ax,log10(nsat),0.50*ones(nnsat,1),'k--','linew',1);  %this is from data
% label{ntrial+1} = 'Data';
% % legend(ax,p,label,'NumColumns',2,'Location','south');
% % title(ax,'Secondary sources removed');
% xlabel(ax,'log_{10}(Saturation)');
% ylabel(ax,'Distance (km)');
% yran = [0 1];
% ylim(ax,yran);
% longticks(ax,2);
% hold(ax,'off');

% ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% for iperc = 1: ntrial
% %   plot(ax,log10(nsat),mprojx32all(:,iperc),'-o','markersize',4,'color',color(iperc,:));
%   plot(ax,log10(nsat),mprojx2all4th(:,iperc),'-','Color',color(iperc,:),'linew',1);
%   scatter(ax,log10(nsat),mprojx2all4th(:,iperc),nsrc4th(:,iperc)/50,color(iperc,:),...
%     'filled','MarkerEdgeColor','k');
% end
% plot(ax,log10(nsat),0.45*ones(nnsat,1),'k--','linew',1);  %this is from data
% ylim(ax,yran);
% % title(ax,'Checkd at 4th stas');
% longticks(ax,2);
% hold(ax,'off');

% ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% for iperc = 1: ntrial
%   plot(ax,log10(nsat),mprojx12all(:,iperc),'-o','markersize',4,'color',color(iperc,:));
% end
% title(ax,'Grouped Total');
% ylim(ax,yran);

% stit = supertit(f.ax,'Med. dist. along min-error direc from each to all others w/i 2 s');
% movev(stit,0.4);

%% Summary of amplitude ratio
% %%%%%%%%% you can choose to plot the actual histogram for each sat and noise
% %%%loop for noise level
% for iperc = 1: ntrial
%   perc = perctrial(iperc);
%   disp(perc);
%   %%%loop for saturation level
%   for insat = 1: nnsat
%     disp(nsat(insat));    
%     f3 = initfig(16,8,2,4); %plot histograms of source amp
%     supertit(f3.ax,sprintf('noise: %.1f, Secondary sources removed & Checkd at 4th stas',...
%       perc));
%     impindepst = imp{insat,iperc};
%     srcampr = srcamprall{insat,iperc};
%     f3.ax(1:3) = plt_deconpk_rat_comb(f3.ax(1:3),srcampr,impindepst,'k','hist');
%     impindepst = imp4th{insat,iperc};
%     srcampr4th = srcamprall4th{insat,iperc};
%     f3.ax(5:end) = plt_deconpk_rat_comb4th(f3.ax(5:end),srcampr4th,impindepst,'k','hist');   
%   end %loop end for saturation level
% end %loop end for noise level
% %%%%%%%%% you can choose to plot the actual histogram for each sat and noise

%%%%%%%%% or only plot the median & mad for each sat and noise
for iperc = 1: ntrial 
  label{iperc} = sprintf('%.1f',perctrial(iperc));
end
label{ntrial+1} = 'Data';
label{ntrial+2} = 'Noise';
srcampralld = allsig.allbstsig.srcamprall;  %using reference from real data
srcampralln = allnoi.allbstnoi.srcamprall;  %using reference from real data
for i = 1: size(srcampralld,2)
  mamprd(i,1) = median(log10(srcampralld(:,i)));
  madamprd(i,1) =  mad(log10(srcampralld(:,i)),1);
  stdamprd(i,1) = std(log10(srcampralld(:,i)));  
  mamprn(i,1) = median(log10(srcampralln(:,i)));
  madamprn(i,1) =  mad(log10(srcampralln(:,i)),1);
  stdamprn(i,1) = std(log10(srcampralln(:,i)));  
end
ref{1} = mamprd;
ref{2} = stdamprd;
ref{3} = mamprn;
ref{4} = stdamprn;
for iperc = 1: ntrial
  for insat = 1: nnsat
    mampr{insat,iperc} = median(log10(srcamprk{insat,iperc}), 1);
    madampr{insat,iperc} = mad(log10(srcamprk{insat,iperc}), 1, 1);
    stdampr{insat,iperc} = std(log10(srcamprk{insat,iperc}), 1);
%     mampr{insat,iperc} = median(srcamprk{insat,iperc}, 1);
%     madampr{insat,iperc} = mad(srcamprk{insat,iperc}, 1, 1);
%     stdampr{insat,iperc} = std(srcamprk{insat,iperc}, 1);
  end
end
f = initfig(8.4,6,2,3); %initialize fig
% orient(f.fig,'landscape');
pltxran = [0.08 0.98]; pltyran = [0.08 0.98]; % optimal axis location
optaxpos(f,2,3,pltxran,pltyran,0.01,0.08);
[f,lgd]=plt_deconpk_rat_stat(f,nsat,label,mampr,stdampr,ref,5*log10(nsrc));
title(lgd(1),lgdtit,'fontsize',7);
title(lgd(2),lgdtit,'fontsize',7);
% stit = supertit(f.ax,'Secondary sources removed');
% movev(stit,0.3);
% ylim(f.ax(4:end),[0 0.4]);
% xlim(f.ax(:),[-0.5 2]);

fname = strcat('amprat_syn',fnsuffix1,'.pdf');
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));
keyboard

%%
srcampr4thalld = allsig.allbstsig.srcampr4thall;  %using reference from real data
srcampr4thalln = allnoi.allbstnoi.srcampr4thall;  %using reference from real data
for i = 1: size(srcampr4thalld,2)
  mampr4thd(i,1) = median(log10(srcampr4thalld(:,i)));
  madampr4thd(i,1) =  mad(log10(srcampr4thalld(:,i)),1);
  stdampr4thd(i,1) = std(log10(srcampr4thalld(:,i)));  
  mampr4thn(i,1) = median(log10(srcampr4thalln(:,i)));
  madampr4thn(i,1) =  mad(log10(srcampr4thalln(:,i)),1);
  stdampr4thn(i,1) = std(log10(srcampr4thalln(:,i)));  
end
ref4th{1} = mampr4thd;
ref4th{2} = stdampr4thd;
ref4th{3} = mampr4thn;
ref4th{4} = stdampr4thn;
for iperc = 1: ntrial
  for insat = 1: nnsat
    mampr4th{insat,iperc} = median(log10(srcampr4thk{insat,iperc}), 1);
    madampr4th{insat,iperc} = mad(log10(srcampr4thk{insat,iperc}), 1, 1);
    stdampr4th{insat,iperc} = std(log10(srcampr4thk{insat,iperc}), 1);
  end
end
f = initfig(8.4,6,2,4); %initialize fig
% orient(f.fig,'landscape');
pltxran = [0.08 0.98]; pltyran = [0.08 0.98]; % optimal axis location
optaxpos(f,2,4,pltxran,pltyran,0.01,0.08);
[f,lgd]=plt_deconpk_rat_stat(f,nsat,label,mampr4th,stdampr4th,ref4th,5*log10(nsrc4th));
title(lgd(1),lgdtit,'fontsize',7);
title(lgd(2),lgdtit,'fontsize',7);
% stit = supertit(f.ax,'Checkd at 4th stas');
% movev(stit,0.3);
% ylim(f.ax(5:end),[0 0.4]);
% xlim(f.ax(:),[-0.5 2]);

fname = strcat('amprat_syn4th',fnsuffix1,'.pdf');
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));

%%%%%%%%% or only plot the median & mad for each sat and noise













