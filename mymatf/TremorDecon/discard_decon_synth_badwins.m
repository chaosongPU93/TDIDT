% discard_decon_synth_badwins.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In 'find_decon_synth_badwins.m', we identify and note down the 25-s windows
% that cannot be well aligned, now this code is to read synthetic catalog 
% detections, and associate them with 25-s windows that they belong to. 
% Then we discard detections that belong to those misaligned windows.
%
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/07/25
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
      savefile1 = strcat('rst_decon_synth_reg',num2str(xaxis),...
        '-',num2str(yaxis),'_nsat',num2str(sat),'_td',num2str(tdura),'.mat');
      savefile2 = strcat('badwins_decon_synth_reg',num2str(xaxis),...
        '-',num2str(yaxis),'_nsat',num2str(sat),'_td',num2str(tdura),'.mat');
    else
      perc = perctrial(iround);
      savefile1 = strcat('rst_decon_synth_onespot','_noi',num2str(perc),'_nsat',...
        num2str(sat),'_td',num2str(tdura),'.mat');
      savefile2 = strcat('badwins_decon_synth_onespot','_noi',num2str(perc),'_nsat',...
        num2str(sat),'_td',num2str(tdura),'.mat');
    end

    %load synthetic catalogs
    load(strcat(workpath,'/synthetics/',savefile1));
    irccran{insat,iround} = allsyn.irccran;  %stores start and end indices of RCC
    imp{insat,iround} = allsyn.imp;
    nsrc{insat,iround} = allsyn.nsrcsum;
    imp4th{insat,iround} = allsyn.imp4th;
    nsrc4th{insat,iround} = allsyn.nsrc4thsum;

    %load misaligned windows
    load(strcat(workpath,'/synthetics/',savefile2));
    badwins{insat,iround} = allsyn.badwins;  
    
  end
end

%%%param for secondary sources removed
impuse = imp;
nsrcuse = nsrc;

% %%%param for further checked at KLNB
% impuse = imp4th;
% nsrcuse = nsrc4th;

%%
%%%Flag to indicate if it is necessary to recalculate everything
% flagrecalc = 0;
flagrecalc = 1;

if flagrecalc
  tic
  for iround = 1: nround
    for insat = 1 : nnsat
      
      nsrci = nsrcuse{insat,iround};
      impi = impuse{insat,iround};

      badwinsi = badwins{insat,iround};

      windows = irccran{insat,iround};
      nwin =  size(windows,1);
      goodwinsi = setdiff((1:nwin)', badwinsi);
      ngoodw = length(goodwinsi);

      iwin = findwhichrange(impi(:,1),windows);

      %find all detections that belong to the good wins
      impilump = [];
      for j = 1: ngoodw
        goodwin = goodwinsi(j);
        impiwin = impi(iwin==goodwin, :);
        impilump = [impilump; impiwin];
      end
      impilump = sortrows(impilump,[1]);

    end
  end
  toc
  
end










