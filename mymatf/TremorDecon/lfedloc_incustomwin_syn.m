% lfedloc_incustomwin_syn.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose is the same as in 'lfedloc_incustomwin.m', but deal with synthetic
% seismograms, either from noise-free synthetics
% (diff source region sizes and saturation level) or from single-spot
% synthetics (diff noise and saturation levels).
%
% To estimate the location difference between each possible pair (the unique
% pairs are from each to all succeeding ones) in the 25-s windows that are
% used in the short-window deconvolutions. 
%
% As the synthetic sources are purely random, and indeed show no migrations, 
% while data shows migrations such that close-in-time events are spaced 
% closer than each to 12.5-s later ones. Then maybe a fair comparison between 
% data and synthetics is the median distance between each event pair in 25-s
% windows.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/05/08
% Last modified date:   2024/07/25
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
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.',...
  num2str(ntol),'.pgc002.',cutout(1:4)));
tlen = trange(:,3)-trange(:,2);
tlensum = sum(tlen);

modname = 'timeoff_plfit_4thsta_160sps.mat';
planefit = load(strcat(rstpath, '/MAPS/',modname));
rmse = planefit.gof{4}.rmse;
offmax = round(2.0*rmse);

%%%specify regime for transformation from time offset to map location
% ftrans = 'interpArmb';
% ftrans = 'interpArmbreloc';
ftrans = 'interpchao';

%% load decon results
tdura = 0.25;  %must be consistent with what synthetics actually used
% tdura = 0.4;

%times of saturation
% nsat=[0.1 0.4 1 2 4 10 20 40 100];
nsat=[0.4 1 2 4 10 20 40 100];
nnsat = length(nsat);

sps = 160;

%%%flag to decide which type of synthetics to use
singleflag = 0;
% singleflag = 1;

if ~singleflag  %%%synthetics from different region sizes and saturation levels
  %%%specify if considering the physical size of each source
  % physicalsize = 1;
  physicalsize = 0;
  
  %%%diameter of physical size
  if physicalsize
    diam=0.15;% 0.5; %0.6; %
  else
    diam=0;
  end
  
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
  nround = nreg;
  
  ttstr1 = {'Noise-free syn, '};
  fnsuffix1 = '';
  
else  %%%synthetics from different noise and saturation levels, sources at a single spot
  %different percent of noise
  perctrial = 0.1*(0:2:16)';
  ntrial = length(perctrial);
  nround = ntrial;
  
  ttstr1 = {'Single-spot syn, '};
  fnsuffix1 = '_onespot';
end

for iround = 1: nround
  for insat = 1: nnsat
    sat = nsat(insat);
    
    if ~singleflag
      xaxis = semia(iround);
      yaxis = semib(iround);
      savefile = strcat('rst_decon_synth_reg',num2str(xaxis),...
        '-',num2str(yaxis),'_nsat',num2str(sat),'_td',num2str(tdura),'.mat');
    else
      perc = perctrial(iround);
      savefile = strcat('rst_decon_synth_onespot','_noi',num2str(perc),'_nsat',...
        num2str(sat),'_td',num2str(tdura),'.mat');
    end
    load(strcat(workpath,'/synthetics/',savefile));
    
    irccran{insat,iround} = allsyn.irccran;  %stores start and end indices of RCC
    off1iw{insat,iround} = allsyn.off1iw; %stores alignment of each 25-s win
    imp{insat,iround} = allsyn.imp;
    nsrc{insat,iround} = allsyn.nsrcsum;
    imp4th{insat,iround} = allsyn.imp4th;
    nsrc4th{insat,iround} = allsyn.nsrc4thsum;
    
  end
end

%%%load real data
% savefile = 'deconv_stats4th_allbstsig.mat';
savefile = 'deconv_stats4th_no23_allbstsig.mat';
load(strcat(rstpath, '/MAPS/',savefile));
nsrcd = allbstsig.nsrc;
impd = allbstsig.impindepall;
nsrc4thd = allbstsig.nsrc4th;
imp4thd = allbstsig.impindep4thall;
irccrand = allbstsig.irccrank;  %stores start and end indices of RCC
off1iwd = allbstsig.off1iwk; %stores alignment of each 25-s win

%%%load synthetic noise
% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
load(strcat(rstpath, '/MAPS/',savefile));
nsrcn = allbstnoi.nsrc;
impn = allbstnoi.impindepall;
nsrc4thn = allbstnoi.nsrc4th;
imp4thn = allbstnoi.impindep4thall;

% keyboard

%%
%%%param for secondary sources removed
impuse = imp;
nsrcuse = nsrc;
supertstr = '3-station';
fnsuffix = [];
impduse = impd;
nsrcduse = nsrcd;
impnuse = impn;
nsrcnuse = nsrcn;

% %%%param for further checked at KLNB
% impuse = imp4th;
% nsrcuse = nsrc4th;
% supertstr = '4-station';
% fnsuffix = '4th';
% impduse = imp4thd;
% nsrcduse = nsrc4thd;
% impnuse = imp4thn;
% nsrcnuse = nsrc4thn;

%% compution of location diff between each to all later ones within a win length
% subwsectar = 25/2;  %target subwindow length in sec
subwsectar = 25;  %target subwindow length in sec
n = 1;  %between N and N-n
m = 1;  %max n to compute

%%%Flag to indicate if it is necessary to recalculate everything
flagrecalc = 0;
% flagrecalc = 1;

savefile = strcat('lfedloc_syn',fnsuffix,fnsuffix1,num2str(subwsectar),'swin.mat');

if flagrecalc
  
  for iround = 1: nround
    for insat = 1 : nnsat
      
      fprintf('round: %d/%d; sat: %d/%d',iround,nround,insat,nnsat);

      nsrci = nsrcuse{insat,iround};
      impi = impuse{insat,iround};
      [imploci, ~] = off2space002(impi(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
      
      dtnn1cat=[];  % time difference between N and N-1, always + 
      dlocnn1cat=[];  % loc difference between N and N-1, sign preserved
      dloc_splnn1cat=[];  % loc difference between N and N-1, sample space, sign preserved
      dt2allcat=[]; % time difference between each to all others, always + 
      dloc2allcat=[]; % loc difference between each to all others, sign preserved
      dloc2all_splcat=[]; % loc difference between each to all others, sample space, sign preserved

      windows = irccran{insat,iround};
      nwin =  size(windows,1);
      iwin = findwhichrange(impi(:,1),windows);
      %     keyboard
      for j = 1: nwin
        impiwin = impi(iwin==j,:);
        implociwin = imploci(iwin==j,:);
        lwinsec = round(((windows(j,2)-windows(j,1))+1)/sps); %notes the actual subwin length in sec
                
        %%%%%% diff loc within short win, able to be combined to analyse later
        %compute the diff loc between N and N-m, and each to all others in the
        %short window, no projection is applied
        [dloc,dt,dloc_spl,dloc2all,dt2all,dloc2all_spl]=...
          dloc_evtcustom(impiwin,implociwin,sps,ftrans,m,'tarvl');
        dtnn1cat=[dtnn1cat; dt{n}];
        dlocnn1cat=[dlocnn1cat; dloc{n}];
        dloc_splnn1cat=[dloc_splnn1cat; dloc_spl{n}];
        dloc2allcat=[dloc2allcat; dloc2all];
        dt2allcat=[dt2allcat; dt2all];
        dloc2all_splcat=[dloc2all_splcat; dloc2all_spl];
        %%%%%%%%%%
      end
      
      dtk{insat,iround} = dtnn1cat;
      dlock{insat,iround} = dlocnn1cat;
      dlocsplk{insat,iround} = dloc_splnn1cat;
      dt2allk{insat,iround} = dt2allcat;
      dloc2allk{insat,iround} = dloc2allcat;
      dloc2all_splk{insat,iround} = dloc2all_splcat;
      
      %%%plot the density, contours and cross-sections
      smoothsigma=5;
      ncont=100;  %num of contour lines
      dx=0.025; dy=0.025; % to distinguish dp near origin, size cannot exceed ~0.053 km 
      cstr={'# events / grid'}; xran=[-4 4]; yran=[-4 4];
      [f,den1d,conmat,muopt,sigmaopt,~,dprojxopt,countnopt,...
        muort,sigmaort,~,dprojxort,countnort] = ...
        plt_srcdloc(dloc2allcat,'km',3,cstr,...
        'o','linear','grid',xran,yran,dx,dy,smoothsigma,ncont);
      close(f.fig);
      %Principal component analysis
      [coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc2allcat);
      [~,~,projxy] = customprojection(dloc2allcat(:,1:2),anglegeo(2));
      mdistprojopt = median(abs(projxy(:,1)));
      mdistprojort = median(abs(projxy(:,2)));
      medprojopt = median(projxy(:,1));
      medprojort = median(projxy(:,2)); 
    
      pcaaxes{insat,iround} = anglegeo;
      pcaellp{insat,iround} = [semia; semib];
      den1dk{insat,iround} = den1d;
      conmatk{insat,iround} = conmat;

      muoptk{insat,iround} = muopt;
      sigmaoptk{insat,iround} = sigmaopt;
      mdistprojoptk{insat,iround} = mdistprojopt;
      medprojoptk{insat,iround} = medprojopt;
      dprojxoptk{insat,iround} = dprojxopt;
      countnoptk{insat,iround} = countnopt;

      muortk{insat,iround} = muort;
      sigmaortk{insat,iround} = sigmaort;
      mdistprojortk{insat,iround} = mdistprojort;
      medprojortk{insat,iround} = medprojort;
      dprojxortk{insat,iround} = dprojxort;
      countnortk{insat,iround} = countnort;      

    end
  end
  
  save(strcat(workpath,'/synthetics/',savefile), 'dtk','dlocsplk','dlock',...
    'dt2allk','dloc2allk','dloc2all_splk','pcaaxes','pcaellp','den1dk','conmatk',...
    'muoptk','sigmaoptk','mdistprojoptk','medprojoptk','dprojxoptk','countnoptk',...
    'muortk','sigmaortk','mdistprojortk','medprojortk','dprojxortk','countnortk');

else
  load(strcat(workpath,'/synthetics/',savefile));
end

%results from data
muoptd = 0.022;
sigmaoptd = 1.09;
mdistprojoptd = 0.58;
medprojoptd = 0.018;
muortd = -0.022;
sigmaortd = 2.05;
mdistprojortd = 1.08;
medprojortd = -0.012;
%results from pure noise
muoptn = -0.008;
sigmaoptn = 1.44;
mdistprojoptn = 0.71;
medprojoptn = 0.0;
muortn = -0.028;
sigmaortn = 2.50;
mdistprojortn = 1.19;
medprojortn = -0.018;

%% summary plot for the distance
nrow = 2; % rows and cols of subplots in each figure
ncol = 3;
widin = 8.4; % size of each figure
htin = 6;
pltxran = [0.07 0.98]; pltyran = [0.18 0.88];
pltxsep = 0.07; pltysep = 0.05;
f = initfig(widin,htin,nrow,ncol);
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

color = gradientblue(nround);

%%%plot the median of loc diff of the whole set in short PCA directions 
ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); ax.XScale='log';
for iround = 1: nround
  for insat = 1: nnsat
    medprojopt(insat,iround) = medprojoptk{insat,iround};  %median of rcc
  end
  % p(iround) = plot(ax,nsat,medprojopt(:,iround),'-','Color',color(iround,:),'linew',1);
  % scatter(ax,nsat,medprojopt(:,iround),20,color(iround,:),'filled','MarkerEdgeColor','k');
  p(iround) = plot(ax,nsat,medprojopt(:,iround),'-',...
    'linew',1,'color',color(iround,:),'marker','o','markersize',3,...
    'markerfacec',color(iround,:));
  if ~singleflag
    label{iround} = sprintf('%.1fx%.1f',2*semia(iround),2*semib(iround));
  else
    label{iround} = sprintf('%.1f',perctrial(iround));
  end
end
p(nround+1) = plot(ax,ax.XLim,[medprojoptd medprojoptd],'k--','linew',1.5);
label{nround+1} = 'Data';
if singleflag
  p(nround+2) = plot(ax,ax.XLim,[medprojoptn medprojoptn],'r--','linew',1.5);
  label{nround+2} = 'Noise';
end
lgd=legend(ax,p,label,'NumColumns',2,'Location','south','fontsize',6);
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
if ~singleflag
  lgdtit = 'Region size (km)';
else
  lgdtit = 'Noise level';
end
title(lgd,strcat(lgdtit),'fontsize',7); %,'; short PCA'
text(ax,0.99,0.93,'Along short PCA','HorizontalAlignment','right','Units','normalized',...
  'fontsize',8);
text(ax,0.02,0.93,'a','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
xlabel(ax,'Saturation');
ylabel(ax,sprintf('med(x) (km)'));
longticks(ax,2);
ylim(ax,[-0.3 0.3]);
xticks(ax,nsat);

%%%plot the median of abs loc diff of the whole set in short PCA directions
ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); ax.XScale='log';
for iround = 1: nround
  for insat = 1: nnsat
    mdistprojopt(insat,iround) = mdistprojoptk{insat,iround};  %median of rcc
  end
  % plot(ax,nsat,mdistprojopt(:,iround),'-','Color',color(iround,:),'linew',1);
  % scatter(ax,nsat,mdistprojopt(:,iround),20,color(iround,:),'filled','MarkerEdgeColor','k');
  plot(ax,nsat,mdistprojopt(:,iround),'-',...
    'linew',1,'color',color(iround,:),'marker','o','markersize',3,...
    'markerfacec',color(iround,:));
end
plot(ax,ax.XLim,[mdistprojoptd mdistprojoptd],'k--','linew',1.5);
if singleflag
  plot(ax,ax.XLim,[mdistprojoptn mdistprojoptn],'r--','linew',1.5);
end
text(ax,0.99,0.93,'Along short PCA','HorizontalAlignment','right','Units','normalized',...
  'fontsize',8);
text(ax,0.02,0.93,'b','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
xlabel(ax,'Saturation');
ylabel(ax,sprintf('med(|x|) (km)'));
longticks(ax,2);
ylim(ax,[0 1.5]);
xticks(ax,nsat);

%%%plot the std of gaussian fit to the cross-sections in short PCA directions
ax=f.ax(3); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); ax.XScale='log';
for iround = 1: nround
  for insat = 1: nnsat
    sigmaopt(insat,iround) = sigmaoptk{insat,iround};  %median of rcc
  end
  % plot(ax,nsat,sigmaopt(:,iround),'-','Color',color(iround,:),'linew',1);
  % scatter(ax,nsat,sigmaopt(:,iround),20,color(iround,:),'filled','MarkerEdgeColor','k');
  plot(ax,nsat,sigmaopt(:,iround),'-',...
    'linew',1,'color',color(iround,:),'marker','o','markersize',3,...
    'markerfacec',color(iround,:));
end
plot(ax,ax.XLim,[sigmaoptd sigmaoptd],'k--','linew',1.5);
if singleflag
  plot(ax,ax.XLim,[sigmaoptn sigmaoptn],'r--','linew',1.5);
end
text(ax,0.99,0.93,'Along short PCA','HorizontalAlignment','right','Units','normalized',...
  'fontsize',8);
text(ax,0.02,0.93,'c','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
xlabel(ax,'Saturation');
ylabel(ax,sprintf('\\sigma of GS fit (km)'));
longticks(ax,2);
ylim(ax,[0 3.5]);
xticks(ax,nsat);

%%
%%%plot the median of loc diff of the whole set in long PCA directions 
ax=f.ax(1+ncol); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); ax.XScale='log';
p=[]; label=[];
for iround = 1: nround
  for insat = 1: nnsat
    medprojort(insat,iround) = medprojortk{insat,iround};  %median of rcc
  end
  % p(iround) = plot(ax,nsat,medprojort(:,iround),'-','Color',color(iround,:),'linew',1);
  % scatter(ax,nsat,medprojort(:,iround),20,color(iround,:),'filled','MarkerEdgeColor','k');
  p(iround) = plot(ax,nsat,medprojort(:,iround),'-',...
    'linew',1,'color',color(iround,:),'marker','o','markersize',3,...
    'markerfacec',color(iround,:));
  if ~singleflag
    label{iround} = sprintf('%.1fx%.1f',2*semia(iround),2*semib(iround));
  else
    label{iround} = sprintf('%.1f',perctrial(iround));
  end
end
p(nround+1) = plot(ax,ax.XLim,[medprojortd medprojortd],'k--','linew',1.5);
label{nround+1} = 'Data';
if singleflag
  p(nround+2) = plot(ax,ax.XLim,[medprojortn medprojortn],'r--','linew',1.5);
  label{nround+2} = 'Noise';
end
% lgd=legend(ax,p,label,'NumColumns',2,'Location','south','fontsize',6);
% set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
% title(lgd,strcat(lgdtit),'fontsize',7); %,'; long PCA'
text(ax,0.99,0.93,'Along long PCA','HorizontalAlignment','right','Units','normalized',...
  'fontsize',8);
text(ax,0.02,0.93,'d','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
xlabel(ax,'Saturation');
ylabel(ax,sprintf('med(x) (km)'));
longticks(ax,2);
ylim(ax,[-0.3 0.3]);
xticks(ax,nsat);

%%%plot the median of abs loc diff of the cross-sections in short PCA directions
ax=f.ax(2+ncol); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); ax.XScale='log';
for iround = 1: nround
  for insat = 1: nnsat
    mdistprojort(insat,iround) = mdistprojortk{insat,iround};  %median of rcc
  end
  % plot(ax,nsat,mdistprojort(:,iround),'-','Color',color(iround,:),'linew',1);
  % scatter(ax,nsat,mdistprojort(:,iround),20,color(iround,:),'filled','MarkerEdgeColor','k');
  plot(ax,nsat,mdistprojort(:,iround),'-',...
    'linew',1,'color',color(iround,:),'marker','o','markersize',3,...
    'markerfacec',color(iround,:));
end
plot(ax,ax.XLim,[mdistprojortd mdistprojortd],'k--','linew',1.5);
if singleflag
  plot(ax,ax.XLim,[mdistprojortn mdistprojortn],'r--','linew',1.5);
end
text(ax,0.99,0.93,'Along long PCA','HorizontalAlignment','right','Units','normalized',...
  'fontsize',8);
text(ax,0.02,0.93,'e','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
xlabel(ax,'Saturation');
ylabel(ax,sprintf('med(|x|) (km)'));
longticks(ax,2);
ylim(ax,[0 1.5]);
xticks(ax,nsat);

%%%plot the std of gaussian fit to the cross-sections in long PCA directions
ax=f.ax(3+ncol); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); ax.XScale='log';
for iround = 1: nround
  for insat = 1: nnsat
    sigmaort(insat,iround) = sigmaortk{insat,iround};  %median of rcc
  end
  % plot(ax,nsat,sigmaort(:,iround),'-','Color',color(iround,:),'linew',1);
  % scatter(ax,nsat,sigmaort(:,iround),20,color(iround,:),'filled','MarkerEdgeColor','k');
  plot(ax,nsat,sigmaort(:,iround),'-',...
    'linew',1,'color',color(iround,:),'marker','o','markersize',3,...
    'markerfacec',color(iround,:));
end
plot(ax,ax.XLim,[sigmaortd sigmaortd],'k--','linew',1.5);
if singleflag
  plot(ax,ax.XLim,[sigmaortn sigmaortn],'r--','linew',1.5);
end
text(ax,0.99,0.93,'Along long PCA','HorizontalAlignment','right','Units','normalized',...
  'fontsize',8);
text(ax,0.02,0.93,'f','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
xlabel(ax,'Saturation');
ylabel(ax,sprintf('\\sigma of GS fit (km)'));
longticks(ax,2);
ylim(ax,[0 3.5]);
xticks(ax,nsat);

fname = strcat('lfedloc2all',num2str(subwsectar),'swin_syn',fnsuffix,fnsuffix1,'.pdf');
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));
% keyboard
