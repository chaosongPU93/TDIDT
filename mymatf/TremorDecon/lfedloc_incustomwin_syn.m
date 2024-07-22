% lfedloc_incustomwin_syn.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose is the same as in 'lfedloc_incustomwin.m', but deal with synthetic
% seismograms, either from noise-free synthetics
% (diff source region sizes and saturation level) or from single-spot
% synthetics (diff noise and saturation levels).
%
% Separated from 'circsrcmodel2.m', now this becomes the script that specifically
% compute the distance or location difference between source pair N and N-n
% in the context of a custom window with a certain length. For example, you can
% ask what is distance between consecutive events among all events included in
% the window of a certain length. Results from each window are also lumped
% together to get a reliable statistics. You can also look at individual
% windows, but that would be the main purpose of 'circsrcmodel2.m'.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/05/08
% Last modified date:   2024/05/08
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
nsat=[0.1 0.4 1 2 4 10 20 40 100];
% nsat=[0.4 1 2 4 10 20 40 100];
nnsat = length(nsat);

sps = 160;

%%%flag to decide which type of synthetics to use
singleflag = 0;

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
    
    irccran{insat,iround} = allsyn.irccran;  %stores start and end indices of RCC
    off1iw{insat,iround} = allsyn.off1iw; %stores alignment of each 25-s win
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

%% compution of location diff between each to all later ones within a win length
% subwsectar = 25/2;  %target subwindow length in sec
subwsectar = 25;  %target subwindow length in sec
n = 1;  %between N and N-n
m = 1;  %max n to compute

%%%Flag to indicate if it is necessary to recalculate everything
% flagrecalc = 0;
flagrecalc = 1;

savefile = strcat('lfedloc2all_syn',fnsuffix,num2str(subwsectar),'swin.mat');

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
        %       mdtnn1(k,1)=median(dt{n});
        %       mdlocnn1(k,:)=median(dloc{n});
        %       mdloc_splnn1(k,:)=median(dloc_spl{n});
        %       mdloc2all(k,:)=median(dloc2all);
        %       mdt2all(k,1)=median(dt2all);
        %       mdloc2all_spl(k,:)=median(dloc2all_spl);
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
      [f,den1d,conmat,muopt,sigmaopt,mdistprojopt,dprojxopt,countnopt,...
        muort,sigmaort,mdistprojort,dprojxort,countnort] = ...
        plt_srcdloc(dloc2allcat,'km',3,cstr,...
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
  
  save(strcat(workpath,'/synthetics/',savefile), 'dtk','dlocsplk','dlock',...
    'dt2allk','dloc2allk','dloc2all_splk','den1dk','conmatk',...
    'muoptk','sigmaoptk','mdistprojoptk','dprojxoptk','countnoptk',...
    'muortk','sigmaortk','mdistprojortk','dprojxortk','countnortk');

end







