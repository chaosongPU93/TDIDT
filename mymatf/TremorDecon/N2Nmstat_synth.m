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

%flag to decide which type of synthetics to use
singleflag = 0;
if ~singleflag  %%%synthetics from different region sizes and saturation levels  
  savefile = 'rst_decon_synth.mat';
  ttstr1 = {'Noise-free syn, '};
  load(savefile);
  nrounds = nreg;
else  %%%synthetics from different noise and saturation levels, sources at a single spot
  savefile = 'rst_synth_onespot.mat';
  ttstr1 = {'Single-spot syn, '};
  load(savefile);
  nrounds = ntrial;
end

%%%load data
savefile = 'deconv_stats4th_allbstsig.mat';
load(strcat(rstpath, '/MAPS/',savefile));

%%
%%%param for secondary sources removed
ttstr2 = '2ndary removed';
fnsuffix = [];
impplt = imp;
denom = 18275;

%below from data for reference
locxyprojall = allbstsig.locxyprojall;
tarvlsplstall = allbstsig.impindepall(:,1);
nsrc = allbstsig.nsrc;
imp = allbstsig.impindepall;
locxyprojalln = allbstnoi.locxyprojall;
tarvlsplstalln = allbstnoi.impindepall(:,1);
nsrcn = allbstnoi.nsrcraw;
impn = allbstnoi.impindepall;
supertstr = 'Secondary sources removed';
fnsuffix = [];
% %%%param for further checked at KLNB
% ttstr2 = 'Checked at KLNB';
% fnsuffix = '4th';
% impplt = imp4th;
% denom = 10547;

mmax=5;
m=1;
supertstr = strcat(ttstr1,ttstr2);

%% summarize the whole catalog, diff arrival time and fractions
f = initfig(15,8,2,round(nnsat/2)); %initialize fig
tit=supertit(f.ax,supertstr);
movev(tit,0.2);
label = [];
for iround = 1: nrounds
  if ~singleflag
    label{iround} = sprintf('a/2=%.2f,b/2=%.2f',semia(iround),semib(iround));
  else
    label{iround} = sprintf('noise=%.1f',perctrial(iround));
  end
end
[f,Nn,frac]=plt_difftime_NNm_syn(f,impplt,mmax,nsat,nrounds,label,sps,m);
keyboard

%%%summarize the whole catalog, only fractions
f = initfig(4,5,1,1); %initialize fig
tit=supertit(f.ax,supertstr);
movev(tit,0.2);
[f,frac]=plt_frac_difftime_NNm_syn(f,impplt,mmax,nsat,nrounds,label,sps,m);

keyboard


%% bin by amp, diff time and frac, for each satur and size/noise level
%%%first bin source by amp, then plot diff arrival time for N and N-1 for each amp bin
for insat = 1: nnsat
  disp(nsat(insat));

  for iround = 1: nrounds
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
  plt_fracdifftime_NNm_syn_binamp(f,impplt,mmax,nsat,nrounds,label,sps,m,nbin);

keyboard


%% fraction of event pairs w/i a diff time cut, and fraction of all catalog
mmax = 10;
f1 = initfig(15,8,2,round(nnsat/2)); %initialize fig
tit=supertit(f1.ax,strcat(supertstr,...
  {', inclusive Frac of srcs in cluster of N & N-m whose diff time w/i 0.25*m+0.125 s'}));
movev(tit,0.2);

f2 = initfig(15,8,2,round(nnsat/2)); %initialize fig
tit=supertit(f2.ax,strcat(supertstr,...
  {', exclusive Frac of srcs in cluster of N & N-m whose diff time w/i 0.25*m+0.125 s'}));
movev(tit,0.2);

label = [];
for iround = 1: nrounds
  if ~singleflag
    label{iround} = sprintf('a/2=%.2f,b/2=%.2f',semia(iround),semib(iround));
  else
    label{iround} = sprintf('noise=%.1f',perctrial(iround));
  end
end

[f1,f2,fracsrc2all,dfracsrc2all]=frac_uniqevt_incluster_syn(f1,f2,impplt,mmax,nsat,nrounds,label,sps);
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

