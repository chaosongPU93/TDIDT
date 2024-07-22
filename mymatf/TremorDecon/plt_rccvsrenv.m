% plt_rccvsrenv.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The script to plot the scatter of running cross-correlation and running
% envelope (or amp). In two frequency bands!!! 
% --The rcc, renv, ramp in high-freq (1.8-6.3 Hz) is computed and stored by 
% 'deconv_ref_4s_exp_4thsta_fn.m'. 
% --The rcc, renv, ramp in a lower passband are computed from 
% 'seisofbursts002_4s_v2.m'.
% --the way of plotting is a bit different from 'seisofbursts002_4s.m' which
% separate out diff ETS years. This script lumps all years, and also plot 
% the corresponding values at src arrivals, and also the LF values.
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/03/21
% Last modified date:   2022/03/21
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

%%%Flag to indicate if it is necessary to recalculate everything
flagrecalc = 0;
% flagrecalc = 1;

% if flagrecalc

[scrsz, res] = pixelperinch(1);

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

ttol = 35;
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.pgc002.',cutout(1:4)));
tlen = trange(:,3)-trange(:,2);
nbst = size(trange,1);

% dates and ets
%%% NOTE: 'dates' is the dates that the tremor is active at the region of interest, so that the
%%% waveform at the stations on these dates you are looking at mainly comes from the sources at the
%%% region of interest. We want to see if there is a noticable change in spectra during the burst
%%% windows on these dates
dates = unique(trange(:,1));
years = unique(floor(dates/1000));
nets = length(years);

ftrans = 'interpchao';

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

%% read daily data, break into windows of segments, plot
%filtering passband for reading data
hisig=6.3; % this will give a similar spectral shape between template and signal
losig=1.8;

sps = 160;

%moving window length in samples for running CC, envelope, etc.
mwlen=0.5*sps;
% mwlen=sps;

%% directly load results in 'deconv_ref_4s_exp_4thsta_fn.m'
% savefile = 'deconv_stats4th_allbstsig.mat';
% savefile = 'deconv_stats4th_no23_allbstsig.mat';
% hf = load(strcat(rstpath, '/MAPS/',savefile));

savefile = strcat('seisbursts_allvari',num2str(losig),'-',num2str(hisig),'.mat');
hf = load(savefile);

rccbst = hf.allbstsig.rccbst;  %rcc of each burst
rampcatk = hf.allbstsig.rampcatk;  %ramp of each burst 
rampcatnormk = hf.allbstsig.rampcatnormk;  %normalized ramp of each burst   
renvcatk = hf.allbstsig.renvcatk;  %renv of each burst 
renvcatnormk = hf.allbstsig.renvcatnormk;  %normalized renv of each burst
rcccatsrcall = hf.allbstsig.rcccatsrcall;  %rcc at src arrival, all bsts lumped 
renvcatsrcall = hf.allbstsig.renvcatsrcall;  %renv at src arrival, all bsts lumped
renvcatnormsrcall = hf.allbstsig.renvcatnormsrcall;  %norm. renv at src, all bsts lumped
rampcatsrcall = hf.allbstsig.rampcatsrcall;  %ramp at src arrival, all bsts lumped
rampcatnormsrcall = hf.allbstsig.rampcatnormsrcall;  %norm. ramp at src, all bsts lumped

%ort. comp.
renvortcatk = hf.allbstsig.renvortcatk;  %renv of each burst 
renvortcatnormk = hf.allbstsig.renvortcatnormk;  %normalized renv of each burst
renvortcatsrcall = hf.allbstsig.renvortcatsrcall;  %renv at src arrival, all bsts lumped
renvortcatnormsrcall = hf.allbstsig.renvortcatnormsrcall;  %norm. renv at src, all bsts lumped

rcczc1catk = hf.allbstsig.rcczc1catk;
renvzc1catk = hf.allbstsig.renvzc1catk;
renvzc1catnormk = hf.allbstsig.renvzc1catnormk;
rcccatsrc1all = hf.allbstsig.rcccatsrc1all;
renvcatsrc1all = hf.allbstsig.renvcatsrc1all;
renvcatnormsrc1all = hf.allbstsig.renvcatnormsrc1all;
rcccatsrczc1all = hf.allbstsig.rcccatsrczc1all;
renvcatsrczc1all = hf.allbstsig.renvcatsrczc1all;
renvcatnormsrczc1all = hf.allbstsig.renvcatnormsrczc1all;

%lump results from all bursts
rcccomb = cat(1,rccbst{:});  %4 cols, 25-s data; whole-win data; 25-s noise; whole-win noi    
rccall = rcccomb(:,1);
rampall = cat(1,rampcatk{:});
rampnormall = cat(1,rampcatnormk{:});
renvall = cat(1,renvcatk{:});
renvnormall = cat(1,renvcatnormk{:});
renvortall = cat(1,renvortcatk{:});
renvortnormall = cat(1,renvortcatnormk{:});
renvratall = renvall./renvortall;
renvnormratall = renvnormall./renvortnormall;

runall = [log10(renvnormall) log10(rampnormall) rccall];
runsrcall = [log10(renvcatnormsrcall) log10(rampcatnormsrcall) rcccatsrcall];

rcczc1all = cat(1,rcczc1catk{:});
renvzc1all = cat(1,renvzc1catk{:});
renvzc1normall = cat(1,renvzc1catnormk{:});
runzc1all = [log10(renvzc1normall) renvzc1all rcczc1all];
runsrc1all = [log10(renvcatnormsrc1all) renvcatsrc1all rcccatsrc1all];
runsrczc1all = [log10(renvcatnormsrczc1all) renvcatsrczc1all rcccatsrczc1all];

%%
% dnpts=round(mwlen/4);
dnpts=10;
runplt = runall(1:dnpts:end,:);
runzcplt = runzc1all;
runsrcplt = runsrczc1all;

%%%put env into bins with equal num, then get the med rcc in that bin 
nbin=10;
[xbin,indbin] = binxeqnum(runplt(:,1),nbin);  %bin by env with same number
for i = 1: nbin
  indi = indbin{i};
  xcnt(i) = median(xbin{i});  %x is renv
  ycnt(i) = median(runplt(indi,3)); %y is rcc
  y1sig(i) = std(runplt(indi,3));  %1-sigma, std
end
% [spear,~] = corr(runplt(:,1),runplt(:,3),'Type','Spearman');
% [ken,~] = corr(runplt(:,1),runplt(:,3),'Type','Kendall');

%%%for the largest bin, break into 2 bins
indi = indbin{nbin};
nbin2 = 2;
runplt90 = runplt(indi,:);
[xbin90,indbin90] = binxeqnum(runplt90(:,1),nbin2);  %bin by env with same number
for i = 1: nbin2
  indi = indbin90{i};
  xcnt90(i) = median(xbin90{i});  %x is renv
  ycnt90(i) = median(runplt90(indi,3)); %y is rcc
  y1sig90(i) = std(runplt90(indi,3));  %1-sigma, std
end

%% summary plot of rcc VS renv for all bursts
nrow = 1; % rows and cols of subplots in each figure
ncol = 2; 
widin = 5.6; % size of each figure
htin = 3.8;
pltxran = [0.07 0.98]; pltyran = [0.18 0.88];
pltxsep = 0.02; pltysep = 0.03;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

scale='log10';
binmethod='grid';
xbinran=[-1 1]; ybinran=[-1 1];
% dx=0.05; dy=0.05;
% msize=50;
% marker='s';
dx=0.02; dy=0.02;
msize=8;
marker='o';

%%%all rcc and renv measurements for all bursts
ax = f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
den1d = density_matrix(runplt(:,1),runplt(:,3),xbinran,ybinran,dx,dy);
den1d = den1d(den1d(:,3)>0, :);
den1d = sortrows(den1d,3);
dum = den1d;
dum(dum(:,3)>1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'linew',0.2);  %, 'MarkerEdgeColor', 'w')
dum = sortrows(den1d,3);
dum(dum(:,3)==1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end  
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'filled','MarkerEdgeColor','none');
%   colormap(ax,flipud(colormap(ax,'kelicol')));
% colormap(ax,jet);
% colormap(ax,'plasma');
colormap(ax,flipud(colormap(ax,'plasma')));
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.03, pos(3), 0.03];
if strcmp(scale,'log10')
  cstr = strcat({'log_{10}(# points / '},binmethod,')');
elseif strcmp(scale,'linear')
  cstr = strcat({'# points / '},binmethod);  
end
c.Label.String = cstr;
%   scatter(ax,runplt(:,1),runplt(:,3),6,'MarkerFaceColor',[.2 .2 .2],'MarkerEdgeColor',...
%     'none','MarkerFaceAlpha',.2); hold on
scatter(ax,median(runplt(:,1)),median(runplt(:,3)),40,'k^','filled');
%   text(ax,median(runplt(:,1))*1.05,median(runplt(:,3))*1.05,sprintf('(%.3f, %.3f)',...
%     median(runplt(:,1)),median(runplt(:,3))),'HorizontalAlignment',...
%     'left','fontsize',10);
  % text(ax,0.98,0.05,sprintf('%d',years(iets)),'unit','normalized',...
  %   'HorizontalAlignment','right','fontsize',10);
%   text(ax,0.2,0.1,sprintf('2.5 prctile: %.3f',prctile(runtmp(:,1),2.5)),'unit','normalized',...
%     'HorizontalAlignment','left','fontsize',10);
text(ax,0.98,0.05,sprintf('%.1f-%.1f Hz',losig,hisig),'unit','normalized',...
  'HorizontalAlignment','right','fontsize',10);
errorbar(ax,xcnt,ycnt,-y1sig,y1sig,'vertical','o',...
  'markersize',4,'color','k','linewidth',0.8,'MarkerEdgeColor','k',...
  'MarkerFaceColor','none','CapSize',4);
errorbar(ax,xcnt90,ycnt90,-y1sig90,y1sig90,'vertical','o',...
  'markersize',4,'color','r','linewidth',0.8,'MarkerEdgeColor','r',...
  'MarkerFaceColor','none','CapSize',4);
% text(ax,0.02,0.1,sprintf('S: %.2f',spear),'unit','normalized',...
%   'HorizontalAlignment','left','fontsize',8);
% text(ax,0.02,0.05,sprintf('K: %.2f',ken),'unit','normalized',...
%   'HorizontalAlignment','left','fontsize',8);
%   text(ax,0.02,0.05,sprintf('%.3f; %.3f',pctlv(1,iets),pctlv(2,iets)),'unit','normalized',...
%     'HorizontalAlignment','left','fontsize',8.5,'color','k');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
axis(ax,'equal');
yran=[-1 1];
ylim(ax,yran);
yticks(ax,yran(1):0.2:yran(2));
xran=[-1 1];
xlim(ax,xran);
% xticks(ax,xran(1):0.2:xran(2));
xtks=[0.1 0.2 0.4 1 2 4 10];
xticks(ax,log10(xtks));
for i=1:length(xtks)
  xtklbls{i}=num2str(xtks(i)); 
end
xticklabels(ax,xtklbls);
longticks(ax,2);
xlabel(ax,'Normalized running envelope','FontSize',10);
ylabel(ax,'Running CC','FontSize',10);
hold(ax,'off');
%%
%%%only rcc and renv measurements at src arrivals
ax = f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% dx=0.02; dy=0.02;
% msize=2*3;
den1dzc = density_matrix(runzcplt(:,1),runzcplt(:,3),xbinran,ybinran,dx,dy);
den1dsrc = density_matrix(runsrcplt(:,1),runsrcplt(:,3),xbinran,ybinran,dx,dy);
%ratio of density
denrat = den1dzc;
denrat(:,3)=den1dsrc(:,3)./den1dzc(:,3);
% keyboard
%plot the rcc and renv at zero-crossing of positive slopes
den1dzc = den1dzc(den1dzc(:,3)>0, :);
scatter(ax,den1dzc(:,1),den1dzc(:,2),msize,[.7 .7 .7],marker,'filled','MarkerEdgeColor','none');
%plot the rcc and renv at src timings at sta 1, colorcoded by density ratio
denrat = denrat(denrat(:,3)>0, :);
denrat = sortrows(denrat,3);
dum = sortrows(denrat,3);
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end  
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'filled','MarkerEdgeColor','none');
colormap(ax,flipud(colormap(ax,'plasma')));
% colormap(ax,'plasma');
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.03, pos(3), 0.03];
if strcmp(scale,'log10')
  cstr = strcat({'log_{10}(ratio of points / '},binmethod,')');
elseif strcmp(scale,'linear')
  cstr = strcat({'ratio of points / '},binmethod);  
end
c.Label.String = cstr;
caxis(ax,[prctile(dum(:,3),1) 0]);
scatter(ax,median(runsrcplt(:,1)),median(runsrcplt(:,3)),40,'k^','filled');
% errorbar(ax,xcnt,ycnt,-y1sig,y1sig,'vertical','o',...
%   'markersize',4,'color','k','linewidth',0.8,'MarkerEdgeColor','k',...
%   'MarkerFaceColor','none','CapSize',4);
% errorbar(ax,xcnt90,ycnt90,-y1sig90,y1sig90,'vertical','o',...
%   'markersize',4,'color','k','linewidth',0.8,'MarkerEdgeColor','k',...
%   'MarkerFaceColor','none','CapSize',4);
% text(ax,0.02,0.1,sprintf('S: %.2f',spearets),'unit','normalized',...
%   'HorizontalAlignment','left','fontsize',8);
% text(ax,0.02,0.05,sprintf('K: %.2f',kenets),'unit','normalized',...
%   'HorizontalAlignment','left','fontsize',8);
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
axis(ax,'equal');
ylim(ax,yran);
yticks(ax,yran(1):0.2:yran(2));
xlim(ax,xran);
xticks(ax,log10(xtks));
xticklabels(ax,xtklbls);
longticks(ax,2);
nolabels(ax,3);
hold(ax,'off');

fname = 'rccvsrenv.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));


%% load the LF computation results
hisiglf=1.8; % this will give a similar spectral shape between template and signal
losiglf=0.5;
savefile = strcat('seisbursts_allvari',num2str(losiglf),'-',num2str(hisiglf),'.mat');
lf=load(savefile);

rccbstlf = lf.allbstsig.rccbst;  %rcc of each burst
rampcatklf = lf.allbstsig.rampcatk;  %ramp of each burst 
rampcatnormklf = lf.allbstsig.rampcatnormk;  %normalized ramp of each burst   
renvcatklf = lf.allbstsig.renvcatk;  %renv of each burst 
renvcatnormklf = lf.allbstsig.renvcatnormk;  %normalized renv of each burst
rcccatsrcalllf = lf.allbstsig.rcccatsrcall;  %rcc at src arrival, all bsts lumped 
renvcatsrcalllf = lf.allbstsig.renvcatsrcall;  %renv at src arrival, all bsts lumped
renvcatnormsrcalllf = lf.allbstsig.renvcatnormsrcall;  %norm. renv at src, all bsts lumped
rampcatsrcalllf = lf.allbstsig.rampcatsrcall;  %ramp at src arrival, all bsts lumped
rampcatnormsrcalllf = lf.allbstsig.rampcatnormsrcall;  %norm. ramp at src, all bsts lumped

%lump results from all bursts
rcccomblf = cat(1,rccbstlf{:});  %4 cols, 25-s data; whole-win data; 25-s noise; whole-win noi    
rccalllf = rcccomblf(:,1);
rampalllf = cat(1,rampcatklf{:});
rampnormalllf = cat(1,rampcatnormklf{:});
renvalllf = cat(1,renvcatklf{:});
renvnormalllf = cat(1,renvcatnormklf{:});
% runalllf = [log10(renvnormalllf) log10(rampnormalllf) rccalllf];
% runsrcalllf = [log10(renvcatnormsrcalllf) log10(rampcatnormsrcalllf) rcccatsrcalllf];
runalllf = [log10(renvalllf) rampalllf rccalllf];
runsrcalllf = [log10(renvcatsrcalllf) rampcatsrcalllf rcccatsrcalllf];

% dnpts=round(mwlen/4);
runpltlf = runalllf(1:dnpts:end,:);
runsrcpltlf = runsrcalllf;

runall = [log10(renvall) rampall rccall];
runsrcall = [log10(renvcatsrcall) rampcatsrcall rcccatsrcall];
runplt = runall(1:dnpts:end,:);
runsrcplt = runsrcall;

runratplt = renvratall(1:dnpts:end);

%%
%%%HF renv vs HF rcc, but colorcoded by LF rcc (median value at the grid point)
nbin=10;
xdata = runplt(:,1); %HF renv
ydata = runplt(:,3); %HF rcc
zdata = runpltlf(:,3);  %LF rcc
[xbin,indbin] = binxeqnum(xdata,nbin);  %bin by env with same number
for i = 1: nbin
  indi = indbin{i};
  xcnt(i) = median(xbin{i});  %x is renv
  ycnt(i) = median(ydata(indi)); %y is rcc
  y1sig(i) = std(ydata(indi));  %1-sigma, std
end
% [spear,~] = corr(xdata,ydata,'Type','Spearman');
% [ken,~] = corr(xdata,ydata,'Type','Kendall');

%%%for the largest bin, break into 2 bins
indi = indbin{nbin};
nbin2 = 2;
xdata90 = xdata(indi);
ydata90 = ydata(indi);
[xbin90,indbin90] = binxeqnum(xdata90,nbin2);  %bin by env with same number
for i = 1: nbin2
  indi = indbin90{i};
  xcnt90(i) = median(xbin90{i});  %x is renv
  ycnt90(i) = median(ydata90(indi)); %y is rcc
  y1sig90(i) = std(ydata90(indi));  %1-sigma, std
end

% dx=0.05; dy=0.05;
% msize=50;
% marker='s';
dx=0.02; dy=0.02;
msize=8;
marker='o';

xbinran=[-2.5 0]; ybinran=[-1 1]; 
[den1d,indgrid] = density_matrix(xdata,ydata,xbinran,ybinran,dx,dy);
den1d = den1d(den1d(:,3)>0, :);
indgrid = indgrid(~cellfun('isempty',indgrid));
zdatamed = median_at_indices(zdata,indgrid);  %sum of moment at each grid point
den1dmed = den1d;
den1dmed(:,3) = zdatamed;
den1dmed = sortrows(den1dmed,3);

% %%
nrow = 1; % rows and cols of subplots in each figure
ncol = 2; 
widin = 8.4; % size of each figure
htin = 4.5;
pltxran = [0.07 0.98]; pltyran = [0.18 0.88];
pltxsep = 0.02; pltysep = 0.03;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

scale='linear';
binmethod='grid';

%%%HF renv vs HF rcc, colorcoded by LF rcc
ax = f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% dum = sortrows(den1dmed,3,'descend');
dum = sortrows(den1dmed,3,'ascend');
dum(dum(:,3)==1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end  
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'filled','MarkerEdgeColor','none');
colormap(ax,flipud(colormap(ax,'plasma')));
% colormap(ax,'plasma');
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.03, pos(3), 0.03];
cstr = strcat({'median LF RCC / '},binmethod);  
c.Label.String = cstr;
% caxis(ax,minmax(dum(:,3)'));
caxis(ax,[prctile(dum(:,3),1) prctile(dum(:,3),99)]);
scatter(ax,median(xdata),median(ydata),40,'k^','filled');
% text(ax,0.98,0.05,sprintf('%.1f-%.1f',losiglf,hisiglf),'unit','normalized',...
%   'HorizontalAlignment','right','fontsize',10);
errorbar(ax,xcnt,ycnt,-y1sig,y1sig,'vertical','o',...
  'markersize',4,'color','k','linewidth',0.8,'MarkerEdgeColor','k',...
  'MarkerFaceColor','none','CapSize',4);
errorbar(ax,xcnt90,ycnt90,-y1sig90,y1sig90,'vertical','o',...
  'markersize',4,'color','r','linewidth',0.8,'MarkerEdgeColor','r',...
  'MarkerFaceColor','none','CapSize',4);
% text(ax,0.02,0.1,sprintf('S: %.2f',spear),'unit','normalized',...
%   'HorizontalAlignment','left','fontsize',8);
% text(ax,0.02,0.05,sprintf('K: %.2f',ken),'unit','normalized',...
%   'HorizontalAlignment','left','fontsize',8);
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
axis(ax,'equal');
yran=[-1 1];
ylim(ax,yran);
yticks(ax,yran(1):0.2:yran(2));
xran=xbinran;
xlim(ax,xran);
% % xticks(ax,xran(1):0.2:xran(2));
% xtks=[0.2 0.4 1 2 4];
% xticks(ax,log10(xtks));
% for i=1:length(xtks)
%   xtklbls{i}=num2str(xtks(i)); 
% end
% xticklabels(ax,xtklbls);
longticks(ax,2);
xlabel(ax,'log_{10}{HF running envelope}','FontSize',10);
ylabel(ax,'HF running CC','FontSize',10);
hold(ax,'off');

% fname = 'rccvsrenv_rcclf.pdf';
% print(f.fig,'-dpdf',...
%   strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

% %%
nbin=10;
xdata = runplt(:,1); %HF renv
ydata = runplt(:,3); %HF rcc
zdata = runratplt(:,1);  %HF opt/ort env ratio
[xbin,indbin] = binxeqnum(xdata,nbin);  %bin by env with same number
for i = 1: nbin
  indi = indbin{i};
  xcnt(i) = median(xbin{i});  %x is renv
  ycnt(i) = median(ydata(indi)); %y is rcc
  y1sig(i) = std(ydata(indi));  %1-sigma, std
end
[spear,~] = corr(xdata,ydata,'Type','Spearman');
[ken,~] = corr(xdata,ydata,'Type','Kendall');

%%%for the largest bin, break into 2 bins
indi = indbin{nbin};
nbin2 = 2;
xdata90 = xdata(indi);
ydata90 = ydata(indi);
[xbin90,indbin90] = binxeqnum(xdata90,nbin2);  %bin by env with same number
for i = 1: nbin2
  indi = indbin90{i};
  xcnt90(i) = median(xbin90{i});  %x is rcc
  ycnt90(i) = median(ydata90(indi)); %y is rcc
  y1sig90(i) = std(ydata90(indi));  %1-sigma, std
end

zdatamed = median_at_indices(zdata,indgrid);  %sum of moment at each grid point
den1dmed = den1d;
den1dmed(:,3) = zdatamed;
den1dmed = sortrows(den1dmed,3);

% %%
% nrow = 1; % rows and cols of subplots in each figure
% ncol = 1; 
% widin = 8; % size of each figure
% htin = 8;
% pltxran = [0.07 0.98]; pltyran = [0.18 0.88];
% pltxsep = 0.02; pltysep = 0.03;
% f = initfig(widin,htin,nrow,ncol);
% optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

scale = 'log10';
% scale = 'linear';
% binmethod='grid';
% xbinran=[-2.5 0]; ybinran=[-1 1];
% dx=1e-2; dy=1e-2;
% msize=4;

%%%HF (opt) renv vs HF rcc, colorcoded by HF opt/ort env ratio
ax = f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
dum = sortrows(den1dmed,3,'ascend');
dum(dum(:,3)==1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end  
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'filled','MarkerEdgeColor','none');
colormap(ax,flipud(colormap(ax,'plasma')));
% colormap(ax,'plasma');
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.03, pos(3), 0.03];
if strcmp(scale,'log10')
  cstr = strcat({'log_{10}(median opt/ort env ratio / '},binmethod,')');
elseif strcmp(scale,'linear')
  cstr = strcat({'median opt/ort env ratio / '},binmethod);  
end
c.Label.String = cstr;
% caxis(ax,minmax(dum(:,3)'));
caxis(ax,[prctile(dum(:,3),1) prctile(dum(:,3),99)]);
% scatter(ax,median(xdata),median(ydata),40,'k^','filled');
% text(ax,0.98,0.05,sprintf('%.1f-%.1f',losiglf,hisiglf),'unit','normalized',...
%   'HorizontalAlignment','right','fontsize',10);
% errorbar(ax,xcnt,ycnt,-y1sig,y1sig,'vertical','o',...
%   'markersize',4,'color','k','linewidth',0.8,'MarkerEdgeColor','k',...
%   'MarkerFaceColor','none','CapSize',4);
% errorbar(ax,xcnt90,ycnt90,-y1sig90,y1sig90,'vertical','o',...
%   'markersize',4,'color','b','linewidth',0.8,'MarkerEdgeColor','b',...
%   'MarkerFaceColor','none','CapSize',4);
% text(ax,0.02,0.1,sprintf('S: %.2f',spear),'unit','normalized',...
%   'HorizontalAlignment','left','fontsize',8);
% text(ax,0.02,0.05,sprintf('K: %.2f',ken),'unit','normalized',...
%   'HorizontalAlignment','left','fontsize',8);
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
axis(ax,'equal');
yran=[-1 1];
ylim(ax,yran);
yticks(ax,yran(1):0.2:yran(2));
xran=xbinran;
xlim(ax,xran);
% % xticks(ax,xran(1):0.2:xran(2));
% xtks=[0.2 0.4 1 2 4];
% xticks(ax,log10(xtks));
% for i=1:length(xtks)
%   xtklbls{i}=num2str(xtks(i)); 
% end
% xticklabels(ax,xtklbls);
nolabels(ax,3);
longticks(ax,2);
% xlabel(ax,'log_{10}{HF running envelope}','FontSize',10);
% ylabel(ax,'HF running CC','FontSize',10);
hold(ax,'off');

% fname = 'rccvsrenv_envrat.pdf';
% print(f.fig,'-dpdf',...
%   strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

fname = 'rccvsrenv_lf.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

