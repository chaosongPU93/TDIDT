% N2Nmstat_synth.mmax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the script to carry out a bunch of N to N-mmax statisical analysis
% to deconvoluted catalog from the synthetics. The bulk is very similar to
% what has been done to real data in 'deconv_ref_4s_exp_4thsta.mmax'. 
% --This script is now compatible with synthetics generated from either 
% different source region sizes and saturation levels with no added noise,
% or all sources from a single spot, but with different noise levels and 
% saturation levels, depending on 'singleflag'.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/10/18
% Last modified date:   2023/11/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
% format short e   % Set the format to 5-digit floating point
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

cutout = 'ellipse';
ttol = 35;
ntol = 3;
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.pgc002.',cutout(1:4)));
% trange = trange(1:end-1,:);
tlen = trange(:,3)-trange(:,2);
tlensum = sum(tlen);

modname = 'timeoff_plfit_4thsta_160sps.mat';
planefit = load(strcat(rstpath, '/MAPS/',modname));
rmse = planefit.gof{4}.rmse;
offmax = round(2.0*rmse);

%% load synthetic seismograms
%flag to decide which type of synthetics to use
singleflag = 0; 
tdura = 0.25;
if ~singleflag  %%%synthetics from different region sizes and saturation levels 
  if tdura == 0.25
    % savefile = 'rst_decon_synth.mat';
    % savefile = 'rst_decon_synthnitm8500.mat';
    savefile = 'rst_decon_synthmedwtcoef.mat';
  elseif tdura == 0.4 
    % savefile = strcat('rst_decon_synth','_td',num2str(tdura),'.mat');
    % savefile = strcat('rst_decon_synthnitm8500','_td',num2str(tdura),'.mat');
    savefile = strcat('rst_decon_synthmedwtcoef','_td',num2str(tdura),'.mat');
  end
  ttstr1 = {'Noise-free syn, '};
  load(savefile);
  nround = nreg;
else  %%%synthetics from different noise and saturation levels, sources at a single spot
  if tdura == 0.25
    % savefile = 'rst_synth_onespot.mat';
    % savefile = 'rst_synth_onespotnitm8500.mat';
    savefile = 'rst_synth_onespotmedwtcoef.mat';
  elseif tdura == 0.4
    % savefile = strcat('rst_synth_onespot','_td',num2str(tdura),'.mat');
    % savefile = strcat('rst_synth_onespotnitm8500','_td',num2str(tdura),'.mat');
    savefile = strcat('rst_synth_onespotmedwtcoef','_td',num2str(tdura),'.mat');
  end
  ttstr1 = {'Single-spot syn, '};
  load(savefile);
  nround = ntrial;
end

%% load data
% savefile = 'deconv_stats4th_allbstsig.mat';
savefile = 'deconv_stats4th_no23_allbstsig.mat';
load(strcat(rstpath, '/MAPS/',savefile));
% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
load(strcat(rstpath, '/MAPS/',savefile));
nbst = size(trange,1);
sps = 160;


%%%param for secondary sources removed
ttstr2 = '2ndary removed';
fnsuffix = [];
impplt = imp;
denom = 18275;
%below from data for reference
locxyprojalld = allbstsig.locxyprojall;
tarvlsplstalld = allbstsig.impindepall(:,1);
nsrcd = allbstsig.nsrc;
impd = allbstsig.impindepall;
nsrcn = allbstnoi.nsrc;
impn = allbstnoi.impindepall;

% %%%param for further checked at KLNB
% ttstr2 = 'Checked at KLNB';
% fnsuffix = '4th';
% impplt = imp4th;
% denom = 10547;
% %below from data for reference
% locxyprojalld = allbstsig.locxyproj4thall;
% tarvlsplstalld = allbstsig.impindep4thall(:,1);
% nsrcd = allbstsig.nsrc4th;
% impd = allbstsig.impindep4thall;
% nsrcn = allbstnoi.nsrc4th;
% impn = allbstnoi.impindep4thall;

mmax=15;
m=1;
supertstr = strcat(ttstr1,ttstr2);

keyboard

%% summarize the whole catalog, distribution of diff arrival time
[~,cntd,Nnd,~,fracd,~,mmaxnonzero,mmaxnonzeron]=...
  plt_srcdifftime_NNm_mmax([],nbst,impd,nsrcd,impn,nsrcn,mmax,sps,1,0);

f = initfig(10.5,6,2,round(nnsat/2)); %initialize fig
optaxpos(f,2,round(nnsat/2),[],[],0.04,0.05);
% tit=supertit(f.ax,supertstr);
% movev(tit,0.2);
label = [];
for iround = 1: nround
  if ~singleflag
    label{iround} = sprintf('a=%.1f, b=%.1f',2*semia(iround),2*semib(iround));
  else
    label{iround} = sprintf('Noise=%.1f',perctrial(iround));
  end
end
label{nround+1} = 'Data';
ref = Nnd;

[f,Nn,frac]=plt_difftime_NNm_syn(f,impplt,nsat,nround,label,sps,m,ref);
orient(f.fig,'landscape');
% keyboard
% print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/agu2023s3f4.pdf'));

%% summarize the whole catalog, only fractions of diff arrival time w/i dtcut
%%%summarize the whole catalog, 
f = initfig(4,4.5,1,1); %initialize fig
% tit=supertit(f.ax,supertstr);
% movev(tit,0.2);
ref = fracd;
[f,frac]=plt_frac_difftime_NNm_syn(f,impplt,nsat,nround,label,sps,m,ref);

keyboard


%% bin by amp, diff time and frac, for each satur and size/noise level
%%%first bin source by amp, then plot diff arrival time for N and N-1 for each amp bin
for insat = 1: nnsat
  disp(nsat(insat));

  for iround = 1: nround
    f = initfig(12,4,1,2); %initialize fig
    tit=supertit(f.ax,supertstr);
    movev(tit,0.3);
    m = 1;
    nbin = 5;
    [f,Nn,ampbincnt,fraci,dtarvlplt,ampplt] = ...
      plt_frac_difftime_NNm_syn_binamp(f,impplt,mmax,insat,iround,sps,m,nbin);
  keyboard

  end
end

keyboard

    
%% bin by amp, only summarize frac, for each satur and size/noise level
m = 1;
nbin = 5;
f = initfig(15,8,2,round(nnsat/2)); %initialize fig
tit=supertit(f.ax,strcat(supertstr, sprintf(', N & N-%d',m)));
movev(tit,0.2);
[f,ampbincnt,fraci,dtarvlplt,ampplt] = ...
  plt_fracdifftime_NNm_syn_binamp(f,impplt,mmax,nsat,nround,label,sps,m,nbin);

keyboard


%% fraction of event pairs w/i a diff time cut, and fraction of all catalog
mmax = 15;
[fracsrc2alld, dfracsrc2alld,mmaxzero]=frac_uniqevt_incluster(nbst,impd,nsrcd,mmax,sps,0);

f1 = initfig(15,8,2,round(nnsat/2)); %initialize fig
tit=supertit(f1.ax,strcat(supertstr,...
  {', inclusive Frac of srcs in cluster of N & N-m whose diff time w/i 0.25*m+0.125 s'}));
movev(tit,0.2);

f2 = initfig(15,8,2,round(nnsat/2)); %initialize fig
tit=supertit(f2.ax,strcat(supertstr,...
  {', exclusive Frac of srcs in cluster of N & N-m whose diff time w/i 0.25*m+0.125 s'}));
movev(tit,0.2);

label = [];
for iround = 1: nround
  if ~singleflag
    label{iround} = sprintf('a/2=%.2f,b/2=%.2f',semia(iround),semib(iround));
  else
    label{iround} = sprintf('Noise=%.1f',perctrial(iround));
  end
end
label{nround+1} = 'Data';

ref = [];
ref{1} = fracsrc2alld;
ref{2} = dfracsrc2alld;
[f1,f2,fracsrc2all,dfracsrc2all]=frac_uniqevt_incluster_syn(f1,f2,impplt,mmax,nsat,nround,label,sps,ref);

for iround = 1: nround
  ax=f1.ax(iround);
  p1(nround+1)=plot(ax,0:1:mmaxzero,fracsrc2alld,'ko-','linew',1,'markersize',4);
end
legend(f1.ax(1),p1,label);

keyboard



%% 
% figure
% subplot(311)
% h=histogram(log10(ampbfdt),'BinWidth',0.25); hold on; %,'NumBins',5
% ax=gca;
% plot(ax,ax.XLim,[200 200],'k--');
% xlabel('log_{10}{Amp}');
% ylabel('Count');
% xlim([-1.5 1]);
% text(0.95,0.9,'Earlier one of the pair','Units','normalized','HorizontalAlignment','right');
% title(sprintf('Between sources N and N-%d (s)',m));
% subplot(312)
% h=histogram(log10(ampafdt),'BinWidth',0.25); hold on;
% ax=gca;
% plot(ax,ax.XLim,[200 200],'k--');
% xlim([-1.5 1]);
% text(0.95,0.9,'Later one of the pair','Units','normalized','HorizontalAlignment','right');
% subplot(313)
% h=histogram(log10(ampdt),'BinWidth',0.25); hold on;
% ax=gca;
% plot(ax,ax.XLim,[200 200],'k--');
% xlim([-1.5 1]);
% text(0.95,0.9,'Mean of the pair','Units','normalized','HorizontalAlignment','right');

%amp distribution of source pair 
figure
subplot(121)
h=histogram(log10(ampplt),'BinWidth',0.25); hold on;
ax=gca;
plot(ax,ax.XLim,[minnum minnum],'k--');
xlim([-1.5 1]);
title('Mean of the pair for data');
ylabel(ax,'Count');
xlabel(ax,sprintf('log_{10}{amp} between sources N and N-%d (s)',m));

subplot(122)
h=histogram(log10(amppltn),'BinWidth',0.25); hold on;
ax=gca;
plot(ax,ax.XLim,[minnumn minnumn],'k--');
xlim([-1.5 1]);
title('Mean of the pair for syn noise');
ylabel(ax,'Count');
xlabel(ax,sprintf('log_{10}{amp} between sources N and N-%d (s)',m));


%%
% figure
% subplot(121)
% eucdistnn1 = allbstsig.distarvlnn1all(:,1);
% eucdistnn2 = allbstsig.distarvlnn2all(:,1);
% 
% eucdistplt = eucdistnn1;
% dlocprojplt = dlocprojnn1bst;
% 
% ax = gca;
% hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% tmp = [dtarvlplt/sps eucdistplt log10(ampplt)];
% tmp = sortrows(tmp,3,'descend');
% scatter(ax,tmp(:,1),tmp(:,2),15,tmp(:,3),'o','filled');
% oldc = colormap(ax,'kelicol');
% newc = flipud(oldc);
% colormap(ax,newc);
% c=colorbar(ax,'SouthOutside');
% c.Label.String = strcat({'log_{10}(amp)'});
% ylabel(ax,'Distance (km)');
% xlabel(ax,sprintf('Diff. arrival time between sources N and N-%d (s)',m));
% xlim(ax,[0 2]);
% ylim(ax,[0 8]);
% hold(ax,'off');
% 
% subplot(122)
% ax = gca;
% hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% tmp = [dtarvlplt/sps abs(dlocprojplt(:,1)) log10(ampplt)];
% tmp = sortrows(tmp,3,'descend');
% scatter(ax,tmp(:,1),tmp(:,2),15,tmp(:,3),'o','filled');
% oldc = colormap(ax,'kelicol');
% newc = flipud(oldc);
% colormap(ax,newc);
% c=colorbar(ax,'SouthOutside');
% c.Label.String = strcat({'log_{10}(amp)'});
% ylabel(ax,'Distance along min-rmse (km)');
% xlabel(ax,sprintf('Diff. arrival time between sources N and N-%d (s)',m));
% xlim(ax,[0 2]);
% ylim(ax,[0 4]);
% hold(ax,'off');

