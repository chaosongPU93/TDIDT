% plt_lfedloc_syn.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose is the same as in 'plt_lfedloc.m', but deal with synthetic
% seismograms, either from noise-free synthetics
% (diff source region sizes and saturation level) or from single-spot
% synthetics (diff noise and saturation levels).
%
% Create a plot combining the lfe diff loc from consecutive events separated
% by less than 0.375 s and that from each to all later ones in 25-s windows.
% Basically this script takes in the interested parts from 
% 'lfedloc_incluster_syn.m' and 'lfedloc_incustomwin_syn.m', and parts for 
% plotting in 'plt_srcdloc.m'. 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/04/25
% Last modified date:   2024/04/25
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

%% load computed diff location results
fnsuffix = [];

%times of saturation
nsat=[0.1 0.4 1 2 4 10 20 40 100];
% nsat=[0.4 1 2 4 10 20 40 100];
% nsat=[0.4 2 10];
nnsat = length(nsat);

savefile = strcat('lfedlocNN1_syn',fnsuffix,'.mat');
rst = load(strcat(workpath,'/synthetics/',savefile));
dloc = rst.dlock;
muopt = rst.muoptk;
sigmaopt = rst.sigmaoptk;
mdistprojopt = rst.mdistprojoptk;
dprojxopt = rst.dprojxoptk;
countnopt = rst.countnoptk;
muoprt = rst.muoptk;
sigmaort = rst.sigmaoptk;
mdistprojort = rst.mdistprojoptk;
dprojxort = rst.dprojxoptk;
countnort = rst.countnoptk;


subwsectar = 25;  %target subwindow length in sec
savefile = strcat('lfedloc2all_syn',fnsuffix,num2str(subwsectar),'swin.mat');
rst25swin = load(strcat(workpath,'/synthetics/',savefile));
dloc2all = rst25swin.dloc2allk;
muopt2all = rst25swin.muoptk;
sigmaopt2all = rst25swin.sigmaoptk;
mdistprojopt2all = rst25swin.mdistprojoptk;
dprojxopt2all = rst25swin.dprojxoptk;
countnopt2all = rst25swin.countnoptk;
muoprt2all = rst25swin.muoptk;
sigmaort2all = rst25swin.sigmaoptk;
mdistprojort2all = rst25swin.mdistprojoptk;
dprojxort2all = rst25swin.dprojxoptk;
countnort2all = rst25swin.countnoptk;

%%%flag to decide which type of synthetics to use
singleflag = 0;

if ~singleflag  %%%synthetics from different region sizes and saturation levels 
  %variation of source region size
  semia = 1.75*(0.6:0.2:2.0);
  semib = 1.25*(0.6:0.2:2.0);
  nreg = length(semia);
  nround = nreg;
else
  %different percent of noise
  perctrial = 0.1*(0:2:16)';
  ntrial = length(perctrial);
  nround = ntrial;
end

%%% between N and N-1 < 0.375 s, SE direction
nrow = 2; ncol = 4;
f=initfig(8.4,6,nrow,ncol); %initialize fig
color = gradientblue(nround);
for insat = 1 : nnsat
  ax=f.ax(insat); hold(ax,'on'); ax.Box='on'; grid(ax,'on');

  for iround = 1: nround    
    p(iround)=plot(ax,dprojxopt{insat,iround},countnopt{insat,iround},'-',...
      'linew',1,'color',color(iround,:));
    if ~singleflag
      label{iround} = sprintf('%.2f, %.2f',semia(iround),semib(iround));
    else
      label{iround} = sprintf('Noi=%.1f',perctrial(iround));
    end
  end
  text(ax,0.98,0.05,sprintf('Sat=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','right');

  if insat == (nrow-1)*ncol+1 
    legend(ax,p,label,'Orientation','horizontal','NumColumns',2,'Location','south');
    xlabel(ax,'Projected location (km)');
    ylabel(ax,'Normalized count');  
  end
  xlim(ax,[-4 4]);  
end

%%% between N and N-1 < 0.375 s, NE direction
nrow = 2; ncol = 4;
f=initfig(8.4,6,nrow,ncol); %initialize fig
color = gradientblue(nround);
for insat = 1 : nnsat
  ax=f.ax(insat); hold(ax,'on'); ax.Box='on'; grid(ax,'on');

  for iround = 1: nround    
    p(iround)=plot(ax,dprojxort{insat,iround},countnort{insat,iround},'-',...
      'linew',1,'color',color(iround,:));
    if ~singleflag
      label{iround} = sprintf('%.2f, %.2f',semia(iround),semib(iround));
    else
      label{iround} = sprintf('Noi=%.1f',perctrial(iround));
    end
  end
  text(ax,0.98,0.05,sprintf('Sat=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','right');

  if insat == (nrow-1)*ncol+1 
    legend(ax,p,label,'Orientation','horizontal','NumColumns',2,'Location','south');
    xlabel(ax,'Projected location (km)');
    ylabel(ax,'Normalized count');  
  end
  xlim(ax,[-4 4]);
end

%%% between each and all succeeding ones, SE direction
nrow = 2; ncol = 4;
f=initfig(8.4,6,nrow,ncol); %initialize fig
color = gradientblue(nround);
for insat = 1 : nnsat
  ax=f.ax(insat); hold(ax,'on'); ax.Box='on'; grid(ax,'on');

  for iround = 1: nround    
    p(iround)=plot(ax,dprojxopt2all{insat,iround},countnopt2all{insat,iround},'-',...
      'linew',1,'color',color(iround,:));
    if ~singleflag
      label{iround} = sprintf('%.1fx%.1f',2*semia(iround),2*semib(iround));
    else
      label{iround} = sprintf('%.1f',perctrial(iround));
    end
  end
  text(ax,0.98,0.05,sprintf('Sat=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','right');

  if insat == (nrow-1)*ncol+1 
    legend(ax,p,label,'Orientation','horizontal','NumColumns',2,'Location','south');
    xlabel(ax,'Projected location (km)');
    ylabel(ax,'Normalized count');  
  end
  xlim(ax,[-4 4]);
end

%%% between each and all succeeding ones, NE direction
nrow = 2; ncol = 4;
f=initfig(8.4,6,nrow,ncol); %initialize fig
color = gradientblue(nround);
for insat = 1 : nnsat
  ax=f.ax(insat); hold(ax,'on'); ax.Box='on'; grid(ax,'on');

  for iround = 1: nround    
    p(iround)=plot(ax,dprojxort2all{insat,iround},countnort2all{insat,iround},'-',...
      'linew',1,'color',color(iround,:));
    if ~singleflag
      label{iround} = sprintf('%.2f, %.2f',semia(iround),semib(iround));
    else
      label{iround} = sprintf('Noi=%.1f',perctrial(iround));
    end
  end
  text(ax,0.98,0.05,sprintf('Sat=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','right');

  if insat == (nrow-1)*ncol+1 
    legend(ax,p,label,'Orientation','horizontal','NumColumns',2,'Location','south');
    xlabel(ax,'Projected location (km)');
    ylabel(ax,'Normalized count');  
  end
  xlim(ax,[-4 4]);
end







