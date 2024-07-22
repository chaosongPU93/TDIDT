% lfedloc_incluster_syn.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose is the same as in 'lfedloc_incluster.m', but deal with synthetic
% seismograms, either from noise-free synthetics
% (diff source region sizes and saturation level) or from single-spot 
% synthetics (diff noise and saturation levels).
%
% Separated from 'circsrcmodel.m', now this becomes the script that specifically
% compute the distance or location difference between source pair N and N-n
% in the context of a event cluster contains m+1 events. For example, you can
% ask what is distance between consecutive events in all types of clusters,
% ie., doublets, triplets, quaduaplets, to the unpper limit (till the point
% there is a nonzero number of clusters in the catalog).  You can also look 
% at individual clusters, but that would be the main purpose of 'circsrcmodel.m'.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/05/08
% Last modified date:   2024/05/08
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
% singleflag = 0;
singleflag = 1;

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
  
else  %%%synthetics from different noise and saturation levels, sources at a single spot
  %different percent of noise
  perctrial = 0.1*(0:2:16)';
  ntrial = length(perctrial);
  nround = ntrial;
  
  ttstr1 = {'Single-spot syn, '};
end

for iround = 1: nround
  for insat = 1: nnsat
    sat = nsat(insat);
    
    if ~singleflag
      xaxis = semia(iround);
      yaxis = semib(iround);
      savefile = strcat('rst_decon_synth_',srcregion(1:3),'_',num2str(xaxis),...
        '-',num2str(yaxis),'_nsat',num2str(sat),'_td',num2str(tdura),'.mat');
    else
      perc = perctrial(iround);
      savefile = strcat('rst_decon_synth_onespot','_noi',num2str(perc),'_nsat',...
        num2str(sat),'_td',num2str(tdura),'.mat');
    end
    load(strcat(workpath,'/synthetics/',savefile));
    
    imp{insat,iround} = allsyn.imp;
    nsrc{insat,iround} = allsyn.nsrcsum;
    imp4th{insat,iround} = allsyn.imp4th;
    nsrc4th{insat,iround} = allsyn.nsrc4thsum;
    
  end
end
% keyboard

%%
%%%param for secondary sources removed
impuse = imp;
nsrcuse = nsrc;
supertstr = '3-station';
fnsuffix = [];

% %%%param for further checked at KLNB
% impuse = imp4th;
% nsrcuse = nsrc4th;
% supertstr = '4-station';
% fnsuffix = '4th';

%% compution of location diff between event pair N and N-1 that are separated by <0.375s
%%%Flag to indicate if it is necessary to recalculate everything
% flagrecalc = 0;
flagrecalc = 1;

savefile = strcat('lfedlocNN1_syn',fnsuffix,'.mat');

if flagrecalc
  
%   mmaxk = cell(nnsat,nround);
%   dlock = cell(nnsat,nround);
%   dlocsplk = cell(nnsat,nround);
%   den1dk = cell(nnsat,nround); 
%   conmatk = cell(nnsat,nround); 
  for iround = 1: nround
    for insat = 1 : nnsat

      nsrci = nsrcuse{insat,iround};
      impi = impuse{insat,iround};
      %%%% FOR 3-station catalog
      %%%determine the 'mmax' for which the resulting number of clusters is nonzero
      timetype = 'tarvl';
      nbst = 1;
      mmax=getmmaxcluster(nbst,impi,nsrci,sps,timetype);

      %%%for a cluster, not only the time separation between N and N-m needs to
      %%%be smaller than 'dtcut', but also the max time separation between each
      %%%consecutive events needs to smaller than 0.25+0.125 s.
      %%%ie, doublet means a cluster of 2 events ONLY occur as doublets
      [catclus,catclusbst,catimp,catuimp,catmedamp,catdtnnm]=...
        evtcluster_ex(nbst,impi,nsrci,mmax,sps,timetype);
      mmaxk{insat,iround} = mmax;

      %%%get all eligible, unique N and N-1 event pairs from all types of exclusive clusters
      n=1;  %to decide what event pair in the cluster to look at, N and N-n
      lumpst = n; %to decide which cluster to start with lumping 
      lumped = mmax; %to decide which cluster to end with lumping 
      [imppairst,imppaired,mampcont,amppair,dt,dloc,dlocspl,imppairuevt]=...
        lump_evtpairs_fromcluster(catclus,n,lumpst,lumped,sps,ftrans);
      dtk{insat,iround} = dt;
      dlock{insat,iround} = dloc;
      dlocsplk{insat,iround} = dlocspl;

      %%%plot the density, contours and cross-sections
      smoothsigma=5;
      ncont=100;  %num of contour lines
      dx=0.025; dy=0.025; % to distinguish dp near origin, size cannot exceed ~0.053 km 
      cstr={'# events / grid'}; xran=[-4 4]; yran=[-4 4];
      [f,den1d,conmat,muopt,sigmaopt,mdistprojopt,dprojxopt,countnopt,...
        muort,sigmaort,mdistprojort,dprojxort,countnort] = ...
        plt_srcdloc(dloc,'km',3,cstr,...
        'o','linear','grid',xran,yran,dx,dy,smoothsigma,ncont);
      den1dk{insat,iround} = den1d;
      conmatk{insat,iround} = conmat;
      muoptk{insat,iround} = muopt;
      sigmaoptk{insat,iround} = sigmaopt;
      mdistprojoptk{insat,iround} = mdistprojopt;
      dprojxoptk{insat,iround} = dprojxopt;
      countnoptk{insat,iround} = countnopt;
      muortk{insat,iround} = muort;
      sigmaortk{insat,iround} = sigmaort;
      mdistprojortk{insat,iround} = mdistprojort;
      dprojxortk{insat,iround} = dprojxort;
      countnortk{insat,iround} = countnort;      
    end
  end

  save(strcat(workpath,'/synthetics/',savefile), 'mmaxk','dtk','dlocsplk','dlock',...
    'den1dk','conmatk','muoptk','sigmaoptk','mdistprojoptk','dprojxoptk',...
    'countnoptk','muortk','sigmaortk','mdistprojortk','dprojxortk','countnortk');

end











