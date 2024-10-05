% find_decon_synth_badwins.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --We have noticed that in the deconvolution of synthetics, for some parameter
% combinations, some windows could not be well aligned in the same way as 
% data (there are only few cases for data that windows cannot be aligned).
% This failture in alignment and then the forced aligning at 0,0 can lead to
% artifacts to the synthetic detections. Ideally, detections need be discarded
% from the final catalog. 
% However, it is impossible to just cut off the related 25-s windows from the
% full simulation. An alternative may be to identify and note down the misaligned
% windows, and discard the detections that belong to the misaligned windows.
%
% Thus, this code aims to rerun the segmentation, and identify the misaligned
% windows. Then later codes need to correlate detections to 25-s windows, and
% get rid of them and save the reminder to a new catalog. 
%
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/07/25
% Last modified date:   2024/07/25
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
  semia = 1.75*(0.6:0.2:2.0);
  semib = 1.25*(0.6:0.2:2.0);
  nreg = length(semia);
end
  
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
  
fnsuffix1 = '';

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
      
    %%%loop for saturation level
    % parfor insat = 1: nnsat
    for insat = 1: nnsat
      sat = nsat(insat);
      fprintf('sat: %f; %d/%d \n',sat,insat,nnsat);

      fracelsew=0; %0.25 %0.5; %0.6; %The ratio of elsewhere events to local events.  Zero for Discrete Ide.

      %%%diameter of physical size
      if physicalsize
        diam=0.15;% 0.5; %0.6; %
      else
        diam=0;
      end
        
      %%%file name prefix of synthetics
      if strcmp(srcregion,'ellipse')
        xaxis = reg(1); %semi-long axis
        yaxis = reg(2); %semi-short axis 
        % xaxis=1.75; %axis length of the same ellipse of my 4-s catalog
        % yaxis=1.25;
        shiftor=[0.2 0.2]; %(in km) %center of the same ellipse of my 4-s catalog
        [xcut,ycut] = ellipse_chao(shiftor(1),shiftor(2),xaxis,yaxis,0.01,45,shiftor);
        fname = ['/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),'_',...
          num2str(xaxis),'-',num2str(yaxis),'.diam',num2str(diam),'.else',num2str(fracelsew,2),...
          'nsat'];
      elseif strcmp(srcregion,'circle')
        shiftor=[0.2 0.2]; %(in km) %center of the same ellipse of my 4-s catalog
        radi = reg(1); %radius
        [xcut,ycut] = circle_chao(shiftor(1),shiftor(2),radi,0.01);    
        fname = ['/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
          '_',num2str(radi),'.diam',num2str(diam),'.else',num2str(fracelsew,2),...
          'nsat'];
      end

      %%%initialize array for combined results from all runs
      badwinsall = [];

      tic
      %%%loop for each run (each segment of synthetics to be combined later)
      for irun = 1: nrun
        fprintf('part: %d/%d \n',irun,nrun);
        
        if tdura == 0.25
          %%%load synthetics of certain saturation level
          STAopt = load(strcat(workpath,fname,num2str(sat),'tdura',...
            num2str(tdura),'T',num2str(round(Twin)),'p',num2str(irun)));
          %%%load sources
          synsrc = load(strcat(workpath,fname,num2str(sat),'tdura',...
            num2str(tdura),'T',num2str(round(Twin)),'p',num2str(irun),'_sources'));
          %%%load starting indices of added sources at sta 1
          synsrcstind = load(strcat(workpath,fname,num2str(sat),'tdura',...
            num2str(tdura),'T',num2str(round(Twin)),'p',num2str(irun),'_stind'));
        elseif tdura == 0.4
          STAopt = load(strcat(workpath,fname,num2str(sat)));
          synsrc = load(strcat(workpath,fname,num2str(sat),'_sources'));
          synsrcstind = load(strcat(workpath,fname,num2str(sat),'_stind'));
        end
          
        if strcmp(distrloc, 'uniform') 
          if strcmp(srcregion,'ellipse')
            xygrid = load([workpath,'/synthetics/synsrcloc.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
              '_',num2str(xaxis),'-',num2str(yaxis),'.diam',num2str(diam),...
              'T',num2str(round(Twin)),'p',num2str(irun),'_grd']);
          elseif strcmp(srcregion,'circle')
            xygrid = load([workpath,'/synthetics/synsrcloc.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
              '_',num2str(radi),'.diam',num2str(diam),'T',num2str(round(Twin)),'p',num2str(irun),'_grd']);  
          end
          tmp = xygrid(synsrc(:,2),:);
          synsrc = [synsrc(:,1) tmp(:,1:4) ones(length(tmp),1)];  %[indtarvl, off12, off13, loce, locn, amp]
        elseif strcmp(distrloc, 'custompdf')
          loc = off2space002(synsrc(:,2:3),sps,ftrans,0);
          synsrc = [synsrc(:,1) loc(:,1:4) ones(length(loc),1)];  %[indtarvl, off12, off13, loce, locn, amp]
        end
        % keyboard
  
        %%%some params related to the synthetics setting
        %       Twin=0.5*3600+3+2*ceil(greenlen/sps); %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
        winlen=Twin*sps+1;
        skiplen=greenlen;
        %%%Here the noise is 'uniform' in time!
        % noistd = 5e-2;
        %       noistd = 2.0e-4;
        noistd = 0; %NOISE-FREE
        rng('default');
        %     seed=(irun-1)*nreg+ireg; 
        %     rng(seed);
        synth=noistd*(randn(winlen+greenlen+2*10,nsta)-0.5); %+2*10 a little extra, for jiggering 2nd & 3rd stations.% for ista=1:nsta

        %% filter data (and templates)
        %%filter data
        hisig=6.3; % this will give a similar spectral shape between template and signal
        losig=1.8;
        optseg = [];
        for ista = 1:nsta
          optseg(:,ista) = Bandpass(STAopt(:,ista), sps, losig, hisig, 2, 2, 'butter');
        end

        %% break into short windows
        %some params
        bufsec = 1;
        msftaddm = bufsec*sps;  %buffer range for later CC alignment, +1 for safety
        rccmwsec = 0.5;
        rccmwlen = rccmwsec*sps;  %window length for computing RCC
        overshoot = rccmwlen/2; %RCC points to the center of computing window, so extra room is needed
        
        indst = 1+msftaddm;
        inded = size(optseg,1)-msftaddm;
        subwsec = 25;
        subwlen = subwsec*sps;
        %since the rcc would lose rccmwlen/2 at both ends, this results in overlapping of 'rccmwlen' in rcc
        %across consecutive windows; if use 'ovlplen' of rccmwlen, then rcc has no overlapping at all
        %       ovlplen = rccmwlen*2;
        ovlplen = rccmwlen;
        windows = movingwins(indst,inded,subwlen,ovlplen,0);
        nwin =  size(windows,1);
    
        off1iw = zeros(nwin,nsta);  % the best alignment between sta2, sta3 wrt sta1 for each subwin
        ccaliw = zeros(nwin,1+nsta-3);  % CC value using the best alignment, including 4th stas
        
        ircccat = [];   % concatenated indices of RCC
        irccran = zeros(nwin,2);  % start and end indices (range) of RCC of all subwins
        rcccat = [];  % average concatenated RCC
        rcc1icat = [];  % concatenated RCC between sta 1 and 4
        rccpaircat = [];  % concatenated RCC between each station pair, order is 12, 13, 23
        mrccwpair = []; % median of concatenated RCC between each station pair
        ccwpair = []; % 0-lag overall cc of each subwin, between each station pair, order is 12, 13, 23
        ccw1i = []; % same as above, but between sta 1 and 4
        
        k = 0;
        badwins = [];
        for iwin = 1: nwin
%           iwin
          isubwst = windows(iwin,1);
          isubwed = windows(iwin,2);
          
          %align records
          optcc = detrend(optseg(isubwst: isubwed,:));
          msftadd = (round(max(abs([off12ran off13ran])))+1)*sps/40;  %+1 for safety
          loffmax = 4*sps/40;
          ccmid = ceil(size(optcc,1)/2);
          ccwlen = round(size(optcc,1)-2*(msftadd+1));  % minus ensures successful shifting of records
          ccmin = 0.01;  % depending on the length of trace, cc could be very low
          iup = 1;    % times of upsampling
          [off12con,off13con,ccaliw(iwin,1)] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
            ccwlen,msftadd,loffmax,ccmin,iup);
          % if a better alignment cannot be achieved, use 0,0
          if off12con == msftadd+1 && off13con == msftadd+1
            fprintf('Short window %d cannot be properly aligned, double-check needed \n',iwin);
            k=k+1;
            badwins(k,1) = iwin;   
          end
          
        end

        %combine results from same sat and reg, but diff runs into same array 
        badwinsall = [badwinsall; badwins];

      end
      toc
      %% Ouput everything in the form of a structure array
      allsyn.badwins = badwinsall;

      %% save file      
      savefile = strcat('badwins_decon_synth','_reg',num2str(xaxis),'-',num2str(yaxis),...
        '_nsat',num2str(sat),'_td',num2str(tdura),'.mat');
      save(strcat(workpath,'/synthetics/',savefile), 'allsyn');

    end %loop end for saturation level
    
  end %loop end for region size
    
end %if need to recalculate




