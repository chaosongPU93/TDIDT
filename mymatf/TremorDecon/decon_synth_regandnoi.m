% decon_synth_regandnoi.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --This script now becomes the driving script to run the deconvolution to
% synthetic seismograms that are generated from different region sizes and
% different saturation levels, plus some level of noise.
% --Synthetics are generated by 'synthshift_regandnoi_multibsts.m' and
% 'synthshift_regandnoi_fn.m'.
%
% --Parameters related to which synthetics are
% read and deconvolved, and those related to the settings of deconvolution
% are defined here.
% --Moreover, this script analyzes the statistics of the resulting deconvolved
% catalogs. These statistics are a bit out-of-date given our new understanding
% to data.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/07/20
% Last modified date:   2024/07/20
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
  semia = 1.75*(0.8:0.2:1.2);
  semib = 1.25*(0.8:0.2:1.2);
  nreg = length(semia);
end

%times of saturation
% nsat=[0.1 0.4 1 2 4 10 20 40 100];
% nsat=[0.4 1 2 4 10 20 40 100];
nsat=[0.4 1 2 4 10 20];  % times of saturation, # of arrivals per duration
nnsat = length(nsat);

%different percent of noise
perctrial = 0.1*(4:2:8)';
ntrial = length(perctrial);

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

fnsuffix1 = '_mix';

%%%Flag to indicate if it is necessary to recalculate everything
% flagrecalc = 0;
flagrecalc = 1;

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
    
    %%%loop for noise level
    for iperc = 1: ntrial
      perc = perctrial(iperc);
      fprintf('noi: %f; %d/%d \n',perc,iperc,ntrial);
      
      %%%loop for saturation level
      %     parfor insat = 1: nnsat
      for insat = 1: nnsat
        sat = nsat(insat);
        fprintf('sat: %f; %d/%d \n',sat,insat,nnsat);
        
        decon_synth_regandnoi_fn(srcregion,reg,perc,sat,distr,Twin,tdura,...
          distrloc,physicalsize,testsrcflag,testfreqflag,pltdataflag,pltgtflag,...
          pltsrcflag,pltsrc4thflag,normflag,tempflag,ftrans,nrun);
        
      end %loop end for saturation level
      
    end %loop end for percent of noise
    
    keyboard
  end %loop end for region size
  
end %if need to recalculate

% keyboard



















