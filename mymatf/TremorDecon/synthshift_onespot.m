% function synthshift_onespot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to 'synthshift_chao.m', this code is to make synthetic seismograms
% from LFE templates.
% But this code tries to put all srcs at one and only one spot, then add
% different level of noise, which is diff precent of noise that is made of
% the same amp spectrum of synthetics and randomized phase, *[times the 
% envelope of data (so that the noise has the same fluctuation as data)]
% 1. in order to see what percentage of noise would make it look like data
% 2. now there can be 2 end-member synthetic tests: one is different saturation
%  level, diff source area, and no noise (or little); one is same spot source,
%  plus diff satur level and diff noise level
% 3. It should have nothing to do with data, the amp spectrum comes from
%   the synthetics with a certain saturation level, (and a certain size)
% 4. Since technically both synthetics have a uniform amp distribution 
%   (we now know synthetics have a smaller amp fluctuation compared with 
%   data), you don't have to multiply the envelope of synthetics with noise. 
%   You just need to match their median envelope for 100% noise 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/09/13
% Last modified date:   2023/09/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
format short e   % Set the format to 5-digit floating point
clear
clc
close all

%%%Flag to indicate if it is necessary to recalculate everything
% flagrecalc = 0;
flagrecalc = 1;

tdura = 0.25;  %must be consistent with what synthetics actually used, start to use on 2023/12/07
% tdura = 0.4;

if ~flagrecalc
  if tdura == 0.25
    %   load(strcat('rst_synth_onespot','_td',num2str(tdura),'.mat'));
    %   load(strcat('rst_synth_onespotnitm8500','_td',num2str(tdura),'.mat'));
    load(strcat('rst_synth_onespotmedwtcoef','_td',num2str(tdura),'.mat'));
  elseif tdura == 0.4
    load('rst_synth_onespot.mat');
%   load('rst_synth_onespotnitm8500.mat');
%   load('rst_synth_onespotmedwtcoef.mat');
  end
  
else
  %% for easy testing
  defval('normflag',0); %whether to normalize templates
  defval('rccmwsec',0.5); %moving win len in sec for computing RCC
  
  set(0,'DefaultFigureVisible','on');
  %set(0,'DefaultFigureVisible','off');
  
  workpath = getenv('ALLAN');
  datapath = strcat(workpath,'/data-no-resp');
  temppath = strcat(datapath, '/templates/PGCtrio/');
  rstpath = strcat(datapath, '/PGCtrio');
  
  [scrsz, resol] = pixelperinch(1);
  
  % tempflag = 'allan';
  tempflag = 'chao';
  
  adatapath = '/home/data2/chaosong/matlab/allan/matfils/';  %path for Allan's data
  
  stas=['PGC  '
    'SSIB '
    'SILB '
    %   'LZB  '
    %   'TWKB '
    %   'MGCB '
    'KLNB '
    ]; % determine the trio and order, here the 1st sta is PGC
  nsta=size(stas,1);         %  number of stations
  
  if strcmp(tempflag,'chao')
    fam = '002';   % family number
    sps = 160;
    templensec = 60;
    ccstack = [];
    for ista = 1: nsta
      fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
        num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao_catnew');
      ccstack(:,ista) = load(fname);
    end
    STA = detrend(ccstack);
  elseif strcmp(tempflag,'allan')
    sps = 100;
    templensec = 120;
    fname = ['PGCopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh_0.5-15Hz_CC1.25-6.5Hz_le20shift';
      'SSIopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh_0.5-15Hz_CC1.25-6.5Hz_le20shift';
      'SILopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh_0.5-15Hz_CC1.25-6.5Hz_le20shift'];
    %   templensec = 60;
    %   fname = ['PGCopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh0.5-15Hz_le14sh';
    %            'SSIopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh0.5-15Hz_le14sh';
    %            'SILopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh0.5-15Hz_le14sh'];
    for ista=1:nsta
      temp=load(strcat(adatapath,fname(ista,:)));
      STA(:,ista)=detrend(temp(:,1))/350;
    end
  end
  
  % %plot the raw templates, not filtered, not best aligned
  % figure
  % ltemp = size(STA,1);
  % subplot(111)
  % hold on
  % for ista = 1: nsta
  %   plot((1:ltemp)/sps,STA(:,ista));
  % end
  
  %%%The below aligns the templates by x-correlation
  if strcmp(tempflag,'chao')
    ist = templensec*sps*4/10;  %not using whole window in case any station has very long-period energy
    ied = templensec*sps*6/10;
  elseif strcmp(tempflag,'allan')
    ist = templensec*sps*6/10;
    ied = templensec*sps*8/10;
  end
  [maxses,imaxses]=max(STA(ist:ied,:),[],1);
  [minses,iminses]=min(STA(ist:ied,:),[],1);
  spread=maxses-minses;
  imaxses = imaxses+ist-1;  %convert to global indices
  iminses = iminses+ist-1;
  % zcrosses=round(0.5*(imaxses+iminses));  % rough, assuming symmetry, Chao 2021/07/16
  %automatically find the zero-crossings
  zcrosses = zeros(nsta,1);
  for ista = 1:nsta
    seg = detrend(STA(iminses(ista): imaxses(ista),ista));  % for zero-crossing timing, only use the main station
    [~,zcrosses(ista)] = min(abs(seg));
    zcrosses(ista) = zcrosses(ista)-1+iminses(ista);  % convert to global index
  end
  %now we want to cut a segment around the zero-crossing at each station
  sampbef=6*sps;
  sampaft=10*sps;
  is=zcrosses-sampbef;
  ie=zcrosses+sampaft;
  for ista=1:nsta
    STAtmp(:,ista)=detrend(STA(is(ista):ie(ista),ista));  % this means templates are 'aligned' at zero-crossings
  end
  %x-correlation independently between each station pair
  mshiftadd=10*sps/40;
  tempxc(:,1)=xcorr(STAtmp(:,2),STAtmp(:,3),mshiftadd,'coeff');
  tempxc(:,2)=xcorr(STAtmp(:,1),STAtmp(:,2),mshiftadd,'coeff'); %shift STAtmp(3,:) to right for positive values
  tempxc(:,3)=xcorr(STAtmp(:,1),STAtmp(:,3),mshiftadd,'coeff'); %shift STAtmp(2,:) to right for positive values
  for ista = 4: nsta
    tempxc(:,ista)=xcorr(STAtmp(:,1),STAtmp(:,ista),mshiftadd,'coeff'); %shift STAtmp(2,:) to right for positive values
  end
  [~,imax]=max(tempxc,[],1);
  imax=imax-(mshiftadd+1); %This would produce a slightly different shift, if filtered seisms were used.
  imax(2)-imax(3)+imax(1);   %enclosed if it equals to 0
  for ista=2:nsta
    STAtmp(mshiftadd+1:end-(mshiftadd+1),ista)=detrend(STAtmp(mshiftadd+1-imax(ista):end-(mshiftadd+1)-imax(ista),ista));
  end
  %normalization
  if normflag
    for ista=1:nsta
      STAtmp(:,ista)=STAtmp(:,ista)/spread(ista); % now templates are 'aligned' indeoendently by x-corr wrt. sta 1
    end
  end
  % figure
  % ltemp = size(STAtmp,1);
  % subplot(111)
  % hold on
  % for ista = 1: nsta
  %   plot((1:ltemp)/sps,STAtmp(:,ista));
  % end
  %%%The above aligns the templates by x-correlation
  
  %%%detrend, taper and bandpass templates
  tmpwlet = STAtmp; % no bandpass
  tmpwletf = STAtmp;  % bandpassed version
  fractap = sps/size(tmpwlet,1);
  for ista = 1: nsta
    %romve mean, linear trend of template
    tmpwlet(:,ista) = detrend(tmpwlet(:,ista));
    %and taper with tukeywin, which is actually a tapered cosine window
    w = tukeywin(size(tmpwlet(:,ista),1),fractap);
    tmpwlet(:,ista) = w.* tmpwlet(:,ista);
    %detrend again for caution
    tmpwlet(:,ista) = detrend(tmpwlet(:,ista));
    %filter the template
    hiwlet=18;
    lowlet=1.8;
    tmpwletf(:,ista) = Bandpass(tmpwlet(:,ista), sps, lowlet, hiwlet, 2, 2, 'butter');
    %detrend again for caution
    tmpwletf(:,ista) = detrend(tmpwletf(:,ista));
  end
  
  %%%constrained CC, so that only 2 offsets are independent
  ccmid = round(size(tmpwletf,1)/2);
  ccwlen = 10*sps;
  loffmax = 5*sps/40;
  ccmin = 0.01;  % depending on the length of trace, cc could be very low
  iup = 1;    % times of upsampling
  [off12con,off13con,cc] = constrained_cc_interp(tmpwletf(:,1:3)',ccmid,...
    ccwlen,mshiftadd,loffmax,ccmin,iup);
  offwlet1i(1) = 0;
  offwlet1i(2) = round(off12con);
  offwlet1i(3) = round(off13con);
  if nsta>3
    for ista = 4: nsta
      [mcoef,offwlet1i(ista)] = xcorrmax(tmpwletf(:,1), tmpwletf(:,ista), mshiftadd, 'coeff');
    end
  end
  
  %%%automatically find the rough zero-crossing time, whose abs. value is closest to 0, whether + or -
  [~,imin] = min(tmpwletf(:,1));
  [~,imax] = max(tmpwletf(:,1));
  [~,zcsta1] = min(abs(tmpwletf(imin:imax,1)));
  zcsta1 = zcsta1+imin-1;
  greenlen = pow2(9)*sps/40;
  green = zeros(greenlen,nsta); % no bandpass
  greenf = zeros(greenlen,nsta);  % bandpassed version
  ppeaks = zeros(nsta,1); % positive peaks
  npeaks = zeros(nsta,1); % negative peaks
  for ista = 1: nsta
    %cut according to the zero-crossing and the time shift from the constrained CC
    green(:,ista) = tmpwlet(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
    greenf(:,ista) = tmpwletf(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
    %detrend again for caution
    green(:,ista) = detrend(green(:,ista));
    greenf(:,ista) = detrend(greenf(:,ista));
    
    %re-find the zero-crossing as the template length has changed
    [~,imin] = min(green(:,ista));
    [~,imax] = max(green(:,ista));
    [~,zcrosses(ista)] = min(abs(green(imin:imax,ista)));
    zcrosses(ista) = zcrosses(ista)+imin-1;
    ppeaks(ista) = imax;
    npeaks(ista) = imin;
    
    [~,imin] = min(greenf(:,ista));
    [~,imax] = max(greenf(:,ista));
    [~,zcrossesf(ista)] = min(abs(greenf(imin:imax,ista)));
    zcrossesf(ista) = zcrossesf(ista)+imin-1;
    ppeaksf(ista) = imax;
    npeaksf(ista) = imin;
    
  end
  %the following is just a check, because now the templates must be best aligned
  ccmid = round(size(greenf,1)/2);
  ccwlen = 4*sps;
  loffmax = 5*sps/40;
  ccmin = 0.01;  % depending on the length of trace, cc could be very low
  iup = 1;    % times of upsampling
  [off12con,off13con,cc] = constrained_cc_interp(greenf(:,1:3)',ccmid,...
    ccwlen,mshiftadd,loffmax,ccmin,iup);
  if ~(off12con==0 && off13con==0)
    disp('Filtered templates are NOT best aligned \n');
  end
  if nsta>3
    for ista = 4: nsta
      [mcoef,mlag] = xcorrmax(greenf(:,1), greenf(:,ista), mshiftadd, 'coeff');
      if mlag~=0   % offset in samples
        fprintf('Filtered templates are NOT best aligned at %s \n',stas(ista,:));
      end
    end
  end
  
  amprat(1,:) = minmax(greenf(:,1)')./minmax(greenf(:,2)');	% amp ratio between max at sta 3 and 2 or min
  amprat(2,:) = minmax(greenf(:,1)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min
  amprat(3,:) = minmax(greenf(:,2)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min
  spread = range(greenf);   % range of the amp of template
  
  %%%plot the unfiltered and filtered templates
%   plt_templates(green,greenf,stas,[],[],lowlet,hiwlet,sps);
  plt_templates(green(:,1:3),greenf(:,1:3),stas(1:3,:),[],[],lowlet,hiwlet,sps);
  
  %just the filtered templates
  % plt_templates_bp(greenf,stas,lowlet,hiwlet,sps);
  % zcrosses
  % ppeaks
  % npeaks
  % zcrossesf
  % ppeaksf
  % npeaksf
  % keyboard
  
   
  %% generate synthetic sources
  %%%Specify the amplitude-frequency (counts) distribution
  distr='UN'  % uniform distribution
  
  %%%specify distribution for source location
  % distrloc = 'custompdf'; %using a custom PDF function
  distrloc = 'uniform'; %uniformly random in a specified region,
  
  Twin=0.5*3600+3+2*ceil(greenlen/sps); %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
  winlen=Twin*sps+1;
%   winlen = size(xnoi,1)+greenlen-1;
%   Twin = winlen/sps;
  skiplen=greenlen;
  %%%for the duration of templates, there are several options
  %%%1. (ppeak-npeak)*2 of the bb template: 44;54;38
  %%%2. direct eyeballing for between zerocrossings: ~65
  %%%3. binned peak-to-peak separation for decently saturated unfiltered synthetics: ~37
  %%%4. similar to 3, but synthetics are filtered first: ~37
%   tdura = 0.4;  % duration from Chao's broadband template, width is about 795-730=65 spls at 160 hz
%   tdura = 0.25; % start to use on 2023/12/07
  satn=1/tdura*Twin   % if just saturated, how many templates can be fit in? a single peak is ~20 samples wide; maybe a little less (at 100 sps).
  %Twin is window duration in seconds. Events can fall within Twin of
  %the start, but the synthetics will go to Twin*(sample rate)+Greenlen to
  %avoid checking for subscript overrrun.  When "synths" are written from "synth" in the
  %subroutine, Greenlen from the start and end will not be written.
  fracelsew=0; %0.25 %0.5; %0.6; %The ratio of elsewhere events to local events.  Zero for Discrete Ide.
  nsat=[0.4 1 2 4 10 20 40 100];  % times of saturation
%   nsat=[0.05];  % times of saturation
  nnsat = length(nsat);
  writes=round(nsat*satn) %how many templates to throw in, under different degrees of saturation
  
  synth=zeros(winlen+greenlen,nsta);
  
  nouts=length(writes);
  seed=round(writes(4)/5e3); %for random number generator
%   seed=2;
  
  %%%specify which time is uniform in time
  % timetype = 'tarvl';
  timetype = 'tori';
  
  %%%specify if considering the physical size of each source
  % physicalsize = 1;
  physicalsize = 0;
  
  %%%specify regime for transformation from time offset to map location
  % ftrans = 'interpArmb';
  % ftrans = 'interpArmbreloc';
  ftrans = 'interpchao';
  
  %%%whether to plot to check the synthetics
  pltsynflag = 0;
%   pltsynflag = 1;
  
  %%%whether to plot to check the new synthetics with added noise
  pltnewsynflag = 0;
%   pltnewsynflag = 1;
  
  %%%specify if forcing a min speration of arrival time for 2 events from the same spot
  % forcesep = 1;
  forcesep = 0;
  
  b=999. %>150 for uniform size distribution
  
  %which location transformation to use
  if strcmp(ftrans,'interpArmb')
    xygrid=load('xygridArmb'); %Made from /ARMMAP/MAPS/2021/interpgrid.m; format PGSS, PGSI, dx, dy
    size(xygrid) %this grid is 40 sps!
  elseif strcmp(ftrans,'interpArmbreloc') || strcmp(ftrans,'interpchao')
    [loc, indinput] = off2space002([],sps,ftrans,0);
    % loc has 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
    iup = sps/40;
    %     xygrid(:,3)=griddata(loc(:,7)/iup,loc(:,8)/iup,loc(:,1),xygrid(:,1),xygrid(:,2));
    %     xygrid(:,4)=griddata(loc(:,7)/iup,loc(:,8)/iup,loc(:,2),xygrid(:,1),xygrid(:,2));
    xygrid = [loc(:,7)/iup,loc(:,8)/iup,loc(:,1),loc(:,2)];
  end
  
  %choose only one spot
%   xygrid = xygrid(xygrid(:,1)==0 & xygrid(:,2)==0, :); %use [0,0]
  xygrid = xygrid(xygrid(:,1)==2/iup & xygrid(:,2)==2/iup, :);  %use [2,2]
  
  diam=0;
  
  [synths,mommax,sources,greensts]=csplaw3c(writes,winlen,skiplen,synth,green,b,...
    xygrid,sps,fracelsew,seed,timetype,ftrans,stas);
 
  %%
  %%% a simple plot the synthetics
  if pltsynflag
    figure
    nrow = length(writes)+1;
    ncol = 1;
    subplot(nrow,ncol,1)
    hold on
    tmp = synth;
    % if nsta == 3
    plot(tmp(:,1),'r');
    plot(tmp(:,2),'b');
    plot(tmp(:,3),'k');
    % else
    %   color = jet(nsta);
    %   for ista = 1: nsta
    %     plot(tmp(:,ista),'Color',color(ista,:));
    %   end
    % end
    xlim([0 40*sps]);
    axranexp(gca,6,20);

    for i = 1: length(writes)
      subplot(nrow,ncol,i+1)
      hold on
      tmp = synths(:,:,i);
      % if nsta == 3
      plot(tmp(:,1),'r');
      plot(tmp(:,2),'b');
      plot(tmp(:,3),'k');
      % else
      %   color = jet(nsta);
      %   for ista = 1: nsta
      %     plot(tmp(:,ista),'Color',color(ista,:));
      %   end
      % end
      text(0.95,0.9,sprintf('%.1f x saturation',nsat(i)),'Units','normalized','HorizontalAlignment',...
        'right');
      xlim([0 40*sps]);
      axranexp(gca,6,20);
    end
    xlabel(sprintf('Samples at %d Hz',sps),'FontSize',12);
    
    %only the lowest saturation to see details
    figure
    i=1;
    tmp = synths(:,:,i); 
    plot(tmp(:,1),'r'); hold on
    plot(tmp(:,2),'b');
    plot(tmp(:,3),'k');
    text(0.95,0.9,sprintf('%.1f x saturation',nsat(i)),'Units','normalized','HorizontalAlignment',...
      'right');
    xlim([0 40*sps]);
    axranexp(gca,6,20);
    xlabel(sprintf('Samples at %d Hz',sps),'FontSize',12);

  end
  
  tmpgrid = xygrid;
  tmpgrid(:,1:2)=round(sps/40*tmpgrid(:,1:2)); % *4 to get to 160 sps from 40.
  
%   insat = 2;  %which saturation to look at
%   n=writes(insat);
%   a = squeeze(sources(1:n,:,insat));
%   b = a(any(a,2),:);
%   source=b;
%   off = tmpgrid(source(:,2),1:2); %note that 'tmpgrid' has the desired sps
%   [~,off(:,3)] = pred_tarvl_at4thsta(stas(4,:),off(:,1),off(:,2));
  
  %% compose new synthetic waveform and carry out deconvolution
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
  pltsrcflag1 = 0;
  
  %%%flag for plot the decon src distribution after removing 2ndary src
%   pltsrcflag2 = 1;
  pltsrcflag2 = 0;
  
  %%%flag for plot the decon src distribution after checking at 4th stas
%   pltsrcflag3 = 1;
  pltsrcflag3 = 0;
  
%   %%%%%%%%%%%%% params of synthetics from diff sat levels and region sizes %%%%%%%%%%%%
%   %%%specify shape of the source region
%   srcregion='ellipse';
%   % srcregion='rectangle';
%   % srcregion='circle';
%   
%   %variation of source region size
%   if strcmp(srcregion,'ellipse')
%     semia = 1.75*(0.6:0.2:2.0);
%     semib = 1.25*(0.6:0.2:2.0);
%     nreg = length(semia);
%   end
%   
%   %'ireg' is chosen based on CC and median sep from noise-free synthetic experiment comparison with data
%   ireg = 5;
%   disp(semia(ireg));
%   disp(semib(ireg));
%   
%   %params of limited source region, subject to variation!
%   if strcmp(srcregion,'circle')
%     shiftor=[0.2 0.2]; %(in km) %center of the same ellipse of my 4-s catalog
%     radi=1.25; %radius
%     [xcut,ycut] = circle_chao(shiftor(1),shiftor(2),radi,0.01);
%   elseif strcmp(srcregion,'ellipse')
%     xaxis = semia(ireg);
%     yaxis = semib(ireg);
%     % xaxis=1.75; %axis length of the same ellipse of my 4-s catalog
%     % yaxis=1.25;
%     shiftor=[0.2 0.2]; %(in km) %center of the same ellipse of my 4-s catalog
%     [xcut,ycut] = ellipse_chao(shiftor(1),shiftor(2),xaxis,yaxis,0.01,45,shiftor);
%   end
%   
%   %%%file name prefix of synthetics
%   if strcmp(srcregion,'ellipse')
%     fname = ['/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),'_',...
%       num2str(xaxis),'-',num2str(yaxis),'.diam',num2str(diam),'.else',num2str(fracelsew,2),...
%       'nsat'];
%   elseif strcmp(srcregion,'circle')
%     fname = ['/synthetics/STAS.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
%       '_',num2str(radi),'.diam',num2str(diam),'.else',num2str(fracelsew,2),...
%       'nsat'];
%   end
%   %%%%%%%%%%%%% params of synthetics from diff sat levels and region sizes %%%%%%%%%%%%
  
  %different percent of noise
  perctrial = 0.1*(0:2:16)';
  ntrial = length(perctrial);
  
  impgrp = cell(nnsat,ntrial);
  imp = cell(nnsat,ntrial);
  imp4th = cell(nnsat,ntrial);
  
  %%%loop for noise level
  for iperc = 1: ntrial
    
    perc = perctrial(iperc);
    disp(perc);
    
    synnew = synths;  %time, station, sat
    
    %%%loop for saturation level
    for insat = 1: nnsat
      %   insat = 1;
      disp(nsat(insat));
      
      %%%load synthetics of certain saturation level
%       optseg = load(strcat(workpath,fname,num2str(nsat(insat))));
      optseg = synths(:,:,insat);
      
%       %%%load sources
%       synsrc = load(strcat(workpath,fname,num2str(nsat(insat)),'_sources'));
%       if strcmp(distrloc, 'uniform')
%         xygrid = load([workpath,'/synthetics/synsrcloc.',distr,'.',int2str(sps),'sps.',srcregion(1:3),...
%           '_',num2str(xaxis),'-',num2str(yaxis),'.diam',num2str(diam),'_grd']);
%         tmp = xygrid(synsrc(:,2),:);
%         synsrc = [synsrc(:,1) tmp(:,1:4) ones(length(tmp),1)];  %[indtarvl, off12, off13, loce, locn, amp]
%       elseif strcmp(distrloc, 'custompdf')
%         [loc, indinput] = off2space002(synsrc(:,2:3),sps,ftrans,0);
%         synsrc = [synsrc(:,1) loc(:,1:4) ones(length(loc),1)];  %[indtarvl, off12, off13, loce, locn, amp]
%       end
      
      %%% make noise by amp spectrum of synthetics and random phase
      %obtain the amp and phase spectra of records via fft
      nfft = size(optseg,1); % number of points in fft
      [xf,ft,amp,pha] = fftspectrum(optseg(:,1:end), nfft, sps,'twosided');
      
      %uniform, random phase with the same span [-pi,pi];
      mpharan = minmax(pha');
      
      rng(seed);
      pharand = (rand(nfft,nsta)-0.5)*2*pi;  %make the phases span from -pi to pi
      
      %construct record with the same amplitude but random phase
      xfrand = amp.*nfft.*exp(1i.*pharand);
      optsegnoi = real(ifft(xfrand,nfft));
      
      xnoi1 = optsegnoi; %synthetic noise
      
      xnoi = xnoi1;
      % ind = find(impindepst(:,1)/sps >= xnoi(1) & impindepst(:,1)/sps <= xnoi(2));
      
      %%%Different scaling schemes
%       %1. To make noise have the median amp of data
%       sclfact = median(envelope(detrend(optseg(:,2:end))));
      %2. To make noise have the same fluctuation as data
      env = envelope(detrend(optseg(:,1:end)));
      % sclfact = env./range(env);  %normalize
      sclfact = env;
%       %3. To make noise have the median amp of decon sources from data
%       sclfact = 6.4461e-01*mean(spread(1:3))/2;

      %scale the noise so that it has amp fluctuation as data
%       xnoi = xnoi .* sclfact;
      
      %normalize so that noi has the same median env as data
      envn = envelope(detrend(xnoi));
      xnoi = xnoi .* median(env) ./ median(envn);
      median(env) - median(envelope(detrend(xnoi)))
      
      %%
      if pltnewsynflag
        [f] = initfig(8,4.5,3,1);
        optaxpos(f,3,1,[],[],[],0.06);
        ax = f.ax(1); hold(ax,'on'); ax.Box = 'on';
        plot(ax,(1:size(optseg,1))/sps,optseg(:,1),'r','linew',1); hold on
        plot(ax,(1:size(optseg,1))/sps,optseg(:,2),'b','linew',1);
        plot(ax,(1:size(optseg,1))/sps,optseg(:,3),'k','linew',1);
        text(ax,0.98,0.9,sprintf('syn from satur=%.1f',nsat(insat)),'Units','normalized',...
          'HorizontalAlignment','right');
        %       xlabel('Samples at 160 sps');
%         xlabel('Time (s)');
%         ylabel('Amplitude');
%         xlim([0 size(optseg,1)]/sps);
        xran = [2 27];
        xlim(ax,xran);
        axranexp(ax,6,10);
        yran = ax.YLim; 
        longticks(ax,4);
  %       ylim([-1.2 1.2]);

        ax = f.ax(2); hold(ax,'on'); ax.Box = 'on';
        plot(ax,(1:size(optseg,1))/sps,xnoi(:,1),'r','linew',1); hold on
        plot(ax,(1:size(optseg,1))/sps,xnoi(:,2),'b','linew',1);
        plot(ax,(1:size(optseg,1))/sps,xnoi(:,3),'k','linew',1);
        text(ax,0.98,0.9,'100% assembled noise','Units','normalized','HorizontalAlignment','right');
        %       xlabel('Samples at 160 sps');
%         xlabel('Time (s)');
%         ylabel('Amplitude');
%         xlim([0 size(optseg,1)]/sps);
        xlim(ax,xran);
        ylim(ax,yran);
        longticks(ax,4);
        
        ax = f.ax(3); hold(ax,'on'); ax.Box = 'on';
        plot(ax,(1:size(optseg,1))/sps,optseg(:,1)+xnoi(:,1),'r','linew',1); hold on
        plot(ax,(1:size(optseg,1))/sps,optseg(:,2)+xnoi(:,2),'b','linew',1);
        plot(ax,(1:size(optseg,1))/sps,optseg(:,3)+xnoi(:,3),'k','linew',1);
        text(ax,0.98,0.9,'syn + 100% noise','Units','normalized','HorizontalAlignment','right');
        %       xlabel('Samples at 160 sps');
        xlabel(ax,'Time (s)','FontSize',10);
        ylabel(ax,'Amplitude','FontSize',10);
%         xlim([0 30]);
        xlim(ax,xran);
        ylim(ax,yran);
        longticks(ax,4);
      end
%       keyboard
      %%
      %some percent of assembled noise
      noiseg = xnoi .*perc;
    
      %add simulated noise to the current burst window
      tmp = synths(:,:,insat);
      synnew(:,:,insat) = tmp + noiseg;
      
%       %%%filter data
%       for ista = 1:nsta
%         synnew(:,ista,insat) = Bandpass(synnew(:,ista,insat), sps, losig, hisig, 2, 2, 'butter');
%       end
      
      %%
      %%%a glance of signal
      if pltnewsynflag
        figure
        subplot(311)
        plot((1:size(synths,1))/sps,synths(:,1,insat),'r','linew',1); hold on
        plot((1:size(synths,1))/sps,synths(:,2,insat),'b','linew',1);
        plot((1:size(synths,1))/sps,synths(:,3,insat),'k','linew',1);
        text(0.95,0.9,sprintf('Single-spot synthetics with sat=%.1f',nsat(insat)),'Units','normalized',...
          'HorizontalAlignment','right');
        %       xlabel('Samples at 160 sps');
        xlabel('Time (s)');
        ylabel('Amplitude');
        % xlim([0 size(xnoi,1)]/sps);
        xlim([0 30]);
        ax = gca; yran = ax.YLim;
        
        subplot(312)
        plot((1:size(synths,1))/sps,noiseg(:,1),'r','linew',1); hold on
        plot((1:size(synths,1))/sps,noiseg(:,2),'b','linew',1);
        plot((1:size(synths,1))/sps,noiseg(:,3),'k','linew',1);
        text(0.95,0.9,sprintf('%d%% noise',round(perc*100)),'Units','normalized','HorizontalAlignment','right');
        %       xlabel('Samples at 160 sps');
        xlabel('Time (s)');
        ylabel('Amplitude');
        % xlim([0 size(xnoi,1)]/sps);
        xlim([0 30]);
        ylim(yran);
        
        subplot(313)
        plot((1:size(synths,1))/sps,synnew(:,1,insat),'r','linew',1); hold on
        plot((1:size(synths,1))/sps,synnew(:,2,insat),'b','linew',1);
        plot((1:size(synths,1))/sps,synnew(:,3,insat),'k','linew',1);
        text(0.95,0.9,'Synthetics+noise','Units','normalized','HorizontalAlignment','right');
        %       xlabel('Samples at 160 sps');
        xlabel('Time (s)');
        ylabel('Amplitude');
        % xlim([0 size(xnoi,1)]/sps);
        xlim([0 30]);
        ylim(yran);
      end
      
%       close all
  
      %% load synthetics of certain saturation level
      STAopt = synnew(:,:,insat);
      
      %%%load sources
      n=writes(insat);
      a = squeeze(sources(1:n,:,insat));
      b = a(any(a,2),:);
      synsrc=b;
      tmp = tmpgrid(synsrc(:,2),:);
      synsrc = [synsrc(:,1) tmp(:,1:4) ones(length(tmp),1)];  %[indtarvl, off12, off13, loce, locn, amp]

      % keyboard
      %%%load starting indices of added sources at sta 1
      synsrcstind = squeeze(greensts{insat}{1});
      
      %% testing, extract and validate the added impulses of template
      if testsrcflag
        wlensec = 30; %how long the window to test
        bufsec = 1; %need some buffer window for later CC alignment
        buffer = bufsec*sps;
        wlensecb = wlensec+bufsec;
        sigstagt = [];  %target window of synthetic signals at all stations
        sigpnstagt = [];  %target window of synthetic signals plus the noise
        sigconvstagt = [];  %target window of reproduced synthetic signals using convolution
        noistagt = [];  %target window of noise
        synsrcgtsta = cell(nsta,1);  %sources whose zero-crossing with time range, [indtarvl rnoff12 rnoff13 amp]
        for ista = 1: nsta
          tgreen = zcrosses(ista);  % anchor time of templates, here choose as zero-crossing time
          indst = 1;  % starting index of the simulated signal to test
          inded = wlensecb*sps+indst-1; % ending index of the simulated signal to test
          source = synsrc;
          if ista == 1
            greenst = synsrcstind; % the starting index of each added template, context of full length
          elseif ista <=3
            %ind - rnoff is the arrival time in index at sta 2 and 3
            %note the sign here, if off12 >0, move 2 to the right to align with 1, meaning that 2 is
            %earlier than 1, ie., tarvl2 < tarvl1. Be consistent all the time
            greenst = synsrcstind-source(:, ista);
          else
            greenst = pred_tarvl_at4thsta(stas(ista,:),source(:,2),source(:,3),source(:,1));
          end
          %%%you don't need all impulses, only some of them contribute to the length of truncated record
          %%%cut out the green indice and sources that contribute
          %%%Note that some zero-crossings might appear later than 'indedsig', but the corresponding start
          %%%index of the green is actually smaller than 'indedsig', ie. their 'greenst' <= 'indedsig'
          induse = find(greenst>=skiplen-greenlen+indst & greenst<=inded+skiplen-1);
          greenst = greenst(induse);
          source = source(induse,:);
          %     keyboard
          
          greenzc = greenst+tgreen; % index of approximate zero-crossing
          source(:, 1) = source(:, 1)+zcrosses(1);  % now 'greenzc' should be the same as 1st col of 'source'
          impamp = zeros(max(greenzc)+20,1);  %20 is a bit arbitary
          for i = 1: length(greenzc)
            impamp(greenzc(i)) = impamp(greenzc(i))+source(i, 6);
          end
          
          %note that some length of the full simulation 'skiplen' was skipped in 'synthgen.m'
%           aaa=synths(:,ista,insat);
          sigpntemp = synths(:,ista,insat);  % simulated signal with white noise
          noitemp = synth(skiplen:winlen,ista);   % account for the skipped part
          sigtemp = sigpntemp-noitemp; % subtract the white noise
          sigpnstagt(:,ista) = sigpntemp(indst:inded); % now focus on the part of interest
          sigstagt(:,ista) = sigtemp(indst:inded);
          noistagt(:,ista) = noitemp(indst:inded);
          
          %direct convolution
          sigconvtmp = conv(green(:,ista),impamp,'full');
          indtrc = tgreen+skiplen;  % starting index for truncation
          sigconvstagt(:,ista) = sigconvtmp(indtrc+indst-1:inded+indtrc-1);  % cut accordingly
          
          %%%Can we reproduce the synthetics with truncated convolution? --YES
          figure
          subplot(411)
          plot(green(:,ista),'r');
          xlim([0 greenlen]);
          legend('Template (unfiltered)');
          text(0.05,0.9,stas(ista,:),'HorizontalAlignment','left','Units','normalized');
          title('Reproduce synthetics with truncated convolution');
          subplot(412)
          imptemp = find(impamp>0);
          p1=stem(imptemp-skiplen, impamp(imptemp),'b','MarkerSize',4); hold on;
          % p1=stem((1:size(impamp))-skiplen, impamp,'b','MarkerSize',4); hold on;
          ax = gca;
          plot([indst indst],ax.YLim,'--','color',[.5 .5 .5]);
          plot([inded inded],ax.YLim,'--','color',[.5 .5 .5]);
          legend(p1,'Synthetic random impulses');
          xran1 = ax.XLim; axpos1 = ax.Position(1);
          subplot(413);
          plot(sigstagt(:,ista),'b'); hold on
          ax=gca; yran = ax.YLim; xran2 = ax.XLim;
          shrink(ax,xran1/xran2,1);
          ax.Position(1)=axpos1;
          plot(sigconvstagt(:,ista),'k');
          legend('Truncated synthetic signal','Truncated signal from convolution');
          subplot(414);
          plot(sigstagt(:,ista)-sigconvstagt(:,ista),'k'); hold on
          legend('Difference'); ylim(yran)
          
          %now 'source' contain impulse triplets that contribute to the trace segment of interest, but we
          %want to get ones whose zero-crossing falls into the window.
          source(:,1) = source(:,1)-skiplen;  % index after skipping
          tmp = source(source(:,1)>=indst & source(:,1)<=inded, :);
          tmp = sortrows(tmp,1,'ascend');
          synsrcgtsta{ista} = tmp; %ideally for each station this should be the same, but coincidence is possible
          %     keyboard
        end
        %notify if sources are the same at diff stations
        if ~isequaln(synsrcgtsta{1,1},synsrcgtsta{2,1}) || ~isequaln(synsrcgtsta{1,1},synsrcgtsta{3,1})
          disp('Extracted sources at different stations are not the same, re-check!');
        end
        %   keyboard
      end
      
      %% filter data (and templates)
      %%filter data
      hisig=6.3; % this will give a similar spectral shape between template and signal
      losig=1.8;
      optseg = [];
      for ista = 1:nsta
        optseg(:,ista) = Bandpass(STAopt(:,ista), sps, losig, hisig, 2, 2, 'butter');
      end
      
      %%%important, does the data actucally have a similar spectral shape to templates
      if testfreqflag
        for i = 1: nsta
          [f] = plt_wletsig_timefreq(green(:,i),STAopt(:,i),greenf(:,i),optseg(:,i),[],[],sps,...
            [lowlet hiwlet],[losig hisig]);
          text(f.ax(end),0.95, 0.95,stas(i,:),'Units','normalized','HorizontalAlignment','right',...
            'FontSize',12);
        end
      end
      
      %% Best alignment for the whole window
      %some params
      bufsec = 1;
      msftaddm = bufsec*sps;  %buffer range for later CC alignment
      rccmwsec = 0.5;
      rccmwlen = rccmwsec*sps;  %window length for computing RCC
      overshoot = rccmwlen/2; %RCC points to the center of computing window, so extra room is needed
      
      %%%obtain a single best alignment based on the entire win
      optcc = detrend(optseg(1+msftaddm: end-msftaddm, :));
      msftadd = 10*sps/40;
      loffmax = 4*sps/40;
      ccmid = ceil(size(optcc,1)/2);
      ccwlen = round(size(optcc,1)-2*(msftadd+1));  % minus ensures successful shifting of records
      ccmin = 0.01;  % depending on the length of trace, cc could be very low
      iup = 1;    % times of upsampling
      [off12con,off13con,ccali,iloopoff,loopoff] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
        ccwlen,msftadd,loffmax,ccmin,iup);
      % if a better alignment cannot be achieved, use 0,0
      if off12con == msftadd+1 && off13con == msftadd+1
        off12con = 0;
        off13con = 0;
        cc12 = xcorr(optcc(:,1), optcc(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
        cc13 = xcorr(optcc(:,1), optcc(:,3),0,'normalized');
        cc23 = xcorr(optcc(:,2), optcc(:,3),0,'normalized');
        ccali = (cc12+cc13+cc23)/3;
        fprintf('Current window cannot be properly aligned, double-check needed \n');
      end
      off1i = zeros(nsta,1);
      off1i(2) = round(off12con);
      off1i(3) = round(off13con);
      
      %%%obtain a single 'observational' alignment between 4th and 1st station, note that this might
      %%%be very different from the empirical prediction from 'empioffset4thsta002'
      %%%--there are also different options here, using sig, or sig-wlet CC as in decon routine
      mcoef = zeros(nsta-3, 1);
      mlag = zeros(nsta-3, 1);
      for ista = 4: nsta
        [mcoef(ista-3),off1i(ista)] = xcorrmax(optcc(:,1), optcc(:,ista), msftadd, 'coeff');
        mlag(ista-3) = off1i(ista);
      end
%       keyboard

      %if you want to avoid the case when the alignment is way off the centroid 
      %by chance while the saturation level is low, you can force it to be the 
      %an average location, this reference value is from the abs location of the
      %the centroid (0.2,0.2), and prediction of off14 from plane fit model
      off1i(2:3) = [2 2];
      [~,off1i(4)] = pred_tarvl_at4thsta(stas(4,:),off1i(2),off1i(3));
      
      %%%Align and compute the RCC based on the entire win, and take that as the input signal!
      optdat = [];  % win segment of interest
      for ista = 1: nsta
        optdat(:, ista) = optseg(1+msftaddm-off1i(ista): end-msftaddm-off1i(ista), ista);
      end
      
      %location of the whole-win best alignment
      [loc0, indinput] = off2space002([off1i(2) off1i(3)],sps,ftrans,0);
      % loc has 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
      
      %%%2022/06/06, do NOT taper whatsoever!!
      sigsta = zeros(size(optdat,1), nsta);
      for ista = 1:nsta
        tmp = optdat(:,ista); %best aligned, filtered
        %detrend and taper only the data, NOT the noise
        tmp = detrend(tmp);
        sigsta(:,ista) = tmp;
      end
      %compute running CC between 3 stations
      [ircc,rcc12,rcc13,rcc23] = RunningCC3sta(sigsta,rccmwlen);
      ircc = ircc-overshoot;
      rccpair = [rcc12 rcc13 rcc23];
      %if only use the mean RCC from pair 12 and 13
      rcc = mean(rccpair(:,[1 2]), 2);
      rccsat(:,insat,iperc) = rcc;
      mrcc = median(rcc);
      
      rcc1i = zeros(length(rcc),nsta-3);
      for ista = 4:nsta
        [~,rcc1i(:,ista-3)] = RunningCC(sigsta(:,1), sigsta(:,ista), rccmwlen);
      end
      
      sigsta = detrend(sigsta(overshoot+1:end-overshoot, :));  %excluding the overshoot
      
      %zero-lag max CC
      cc12 = xcorr(sigsta(:,1), sigsta(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
      cc13 = xcorr(sigsta(:,1), sigsta(:,3),0,'normalized');
      cc23 = xcorr(sigsta(:,2), sigsta(:,3),0,'normalized');
      cc1i = zeros(nsta-3,1);
      for ista = 4:nsta
        cc1i(ista-3) = xcorr(sigsta(:,1), sigsta(:,ista),0,'normalized');
      end
      ccpair = [cc12 cc13 cc23];
      mcc = (cc12+cc13)/2;  %if only use the pair 12 and 13
      
      %%
      %if you want a qucik look of the data that feeds to the deconvolution
      if pltdataflag
        figure; hold on
        lsig = size(sigsta,1);
        if nsta == 3
          plot((1:lsig)/sps,sigsta(:,1),'r');
          plot((1:lsig)/sps,sigsta(:,2),'b');
          plot((1:lsig)/sps,sigsta(:,3),'k');
        else
          color = jet(nsta);
          for ista = 1: nsta
            p(ista)=plot((1:lsig)/sps,sigsta(:,ista),'Color',color(ista,:));
            %       label{i}=stas(i,:);
          end
        end
        ax = gca;
        axsym(ax);
        plot((1:lsig)/sps,rcc*ax.YLim(2),'o','color',[.6 .6 .6],'markersize',2);
        text(0.95,0.1,sprintf('Saturation: %.1f',nsat(insat)),'Units','normalized','HorizontalAlignment','right');
        xlim([0 20]);
        legend(p,stas);
      end
      
      %% ground truth of conrtibuting sources AFTER alignment
      lsig = size(STAopt,1); %length of original synthetics
      synsrcgtsta = cell(nsta,1);  %sources whose zero-crossing with time range, [indtarvl rnoff12 rnoff13 amp]
      for ista = 1: nsta
        tgreen = zcrosses(ista);  % anchor time of templates, here choose as zero-crossing time
        source = synsrc;
        
        indst = 1+msftaddm+overshoot-off1i(ista);  % starting index of signal in terms of 'STAopt'
        inded = lsig-msftaddm-overshoot-off1i(ista);  % ending index of signal in terms of 'STAopt'
        
        if ista == 1
          greenst = synsrcstind; % the starting index of each added template, context of full length
        elseif ista <=3
          greenst = synsrcstind-source(:, ista)-off1i(ista);
        else
          greenst = pred_tarvl_at4thsta(stas(ista,:),source(:,2),source(:,3),source(:,1),off1i(ista));
        end
        
        %%%you don't need all impulses, only some of them contribute to the length of truncated record
        %%%cut out the green indice and sources that contribute
        %%%Note that some zero-crossings might appear later than 'indedsig', but the corresponding start
        %%%index of the green is actually smaller than 'indedsig', ie. their 'greenst' <= 'indedsig'
        induse = find(greenst>=skiplen-greenlen+indst & greenst<=inded+skiplen-1);
        greenst = greenst(induse);
        source = source(induse,:);
        
        greenzc = greenst+tgreen; % index of approximate zero-crossing
        source(:,1) = source(:,1)+zcrosses(1);  % now 'greenzc' should be the same as 1st col of 'source'
        
        %now 'source' contain impulse triplets that contribute to the trace segment of interest, but we
        %want to get ones whose zero-crossing falls into the window.
        source(:,1) = source(:,1)-skiplen;  % index after skipping
        tmp = source(source(:,1)>=indst & source(:,1)<=inded, :);
        tmp = sortrows(tmp,1,'ascend');
        synsrcgtsta{ista} = tmp; %ideally for each station this should be the same, but coincidence is possible
      end
      % %notify if sources are the same at diff stations
      % if ~isequaln(synsrcgtsta{1,1},synsrcgtsta{2,1}) || ~isequaln(synsrcgtsta{1,1},synsrcgtsta{3,1})
      %   disp('Extracted sources at different stations are not the same, re-check!');
      % else
      ista = 1;
      synsrcgt = synsrcgtsta{ista};
      %store the info of synthetic impulses, same format as in paired deconvolution,
      %9 cols, [ind1 amp1 ind2 amp2 ind3 amp3 off12 off13 off23]
      impgt = zeros(size(synsrcgt,1), 9);
      impgt(:,1) = synsrcgt(:,1)-msftaddm-overshoot;  % cut out the buffer
      impgt(:,[2 4 6]) = repmat(synsrcgt(:,6),1,3);
      %the indices of synthetic impulses need to be shifted too
      impgt(:,3) = impgt(:,1)-synsrcgt(:,2)+off1i(2); % note the sign is consistent!
      impgt(:,5) = impgt(:,1)-synsrcgt(:,3)+off1i(3);
      impgt(:,7) = synsrcgt(:,2);  % no need to shift offset, because they are 'true' offset
      impgt(:,8) = synsrcgt(:,3);
      impgt(:,9) = impgt(:,8)-impgt(:,7);  % off23, == off13 - off12 == tarvl2 - tarvl3
      impgtst = sortrows(impgt,1,'ascend');
      % end
      
      %amp and density at each pixel for ground truth sources
      ampgt = mean(impgt(:,[2 4 6]),2); %amp for all LFE catalog
      [density1d, inddup] = density_pixel(impgt(:,7),impgt(:,8));
      ampgtsum = sum_at_indices(ampgt,inddup);
      [impgtloc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
      density1d = [impgtloc(:,1:2) density1d(:,3)];
      ampgtsum1d = sortrows([impgtloc(:,1:2) ampgtsum], 3);
%       keyboard
      %%%if you want to plot the ground truth
      if pltgtflag
        %plot the ground truth source offset
        spsscale = sps/40;
        loff_max = 6*spsscale;  % maximum allowed shift is 4 samples at 40 Hz, it is also true from the synthetic distribution
        xran = [-loff_max+off1i(2)-1 loff_max+off1i(2)+1];
        yran = [-loff_max+off1i(3)-1 loff_max+off1i(3)+1];
        offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
        offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
        cran = [0 lsig];
        f.fig = figure;
        f.fig.Renderer = 'painters';
        ax=gca;
        [ax,torispl,mamp] = plt_decon_imp_scatter(ax,impgt,xran,yran,cran,offxran,offyran,...
          sps,50,'mean','tarvl');
        scatter(ax,off1i(2),off1i(3),20,'ks','filled','MarkerEdgeColor','k');
        title(ax,'Ground truth');
        
        %plot the transformed ground truth source map locations
        xran = [-4 4];
        yran = [-4 4];
        cran = [0 lsig/sps];
        f.fig = figure;
        f.fig.Renderer = 'painters';
        ax=gca;
        [ax] = plt_decon_imp_scatter_space(ax,impgt,xran,yran,cran,offxran,...
          offyran,sps,50,ftrans,'mean','tarvl');
        title(ax,'Ground truth');

        % %%%plot the cumulative density and summed amp of detections
        % cstr = {'# detections / pixel'; 'amp sum / pixel'};
        % [f] = plt_sum_pixel(density1d,ampgtsum1d,[-4 4],[-4 4],20,cstr,'o','linear');
        % supertit(f.ax,'Ground truth');
        
        % %plot the ground truth source map locations from grid
        % f.fig = figure;
        % f.fig.Renderer = 'painters';
        % ax=gca;
        % hold(ax,'on');
        % ax.Box='on'; grid(ax,'on');
        % plot(ax,[xran(1) xran(2)],[0 0],'--','Color',[.3 .3 .3]); hold(ax,'on')
        % plot(ax,[0 0],[yran(1) yran(2)],'--','Color',[.3 .3 .3]);
        % scatter(ax,synsrcgt(:,4),synsrcgt(:,5),synsrcgt(:,6)*30,impgt(:,1)/sps,'filled','MarkerEdgeColor',[.5 .5 .5]);
        % scatter(ax,loc0(:,1),loc0(:,2),30,'k','filled');
        % scatter(ax,shiftor(:,1),shiftor(:,2),30,'k','linew',1);
        % plot(ax,xcut,ycut,'k-','linew',2);
        % oldc = colormap(ax,'kelicol');
        % newc = flipud(oldc);
        % colormap(ax,newc);
        % c=colorbar(ax);
        % caxis(ax, cran);
        % c.Label.String = sprintf('Arrival time at PGC (s)');
        % text(ax,0.98,0.05,sprintf('%d events',size(synsrcgt,1)),'Units','normalized',...
        %   'HorizontalAlignment','right','FontSize',9);
        % ax.YAxis.FontSize = 8;
        % ax.XAxis.FontSize = 8;
        % xlabel(ax,'E (km)','FontSize',11);
        % ylabel(ax,'N (km)','FontSize',11);
        % axis(ax, 'equal');
        % xlim(ax,xran); xticks(ax,xran(1): 1 : xran(2));
        % ylim(ax,yran); yticks(ax,yran(1): 1 : yran(2));
      end
%       keyboard
      
      %% independent deconvolution at each station
      %%%finalize the signal, noise, and template (Green's function)
      sigdecon = [];
      pred = [];
      ampit = [];
      for ista = 1:3
        wlet = greenf(:,ista);  %template here is best aligned, tapered, linear trend removed, filtered
        lwlet = length(wlet);
        sig = sigsta(:,ista); %best aligned, filtered, tapered
        lsig = length(sig);
        noi = zeros(lsig,1);
        
        dt = 1/sps;  % sampling interval
        twlet = zcrossesf(ista)*dt;
        width = 2.5;  % width for Gaussian filter
        dres_min = 0.5;  % tolerance, percentage change in residual per iteration
        mfit_min = 5e-1;  % tolerance, norm of the residual/misfit
        tlen = ceil(lsig/sps);
%         tdura = 0.4;  % estimate from the broadband template from fam 002
%         nit_max = round(1.5*1/tdura*tlen*nsat(insat));  % max numer of iterations
%         nimp_max = round(1/tdura*tlen*nsat(insat));%a single peak is ~20 samples wide; maybe a little less (at 100 sps). ~0.4s, 1/0.4=2.5
%         tdura = 0.25; 
        nit_max = round(1.5*1/tdura*(tlen));  % max numer of iterations
%         nit_max = 8500;
        nimp_max = round(1/tdura*(tlen));%a single peak is ~20 samples wide; maybe a little less (at 100 sps). ~0.4s, 1/0.4=2.5
        fpltit = 0;  % plot flag for each iteration
        fpltend = 0;  % plot flag for the final iteration
        fpltchk = 0; % plot flag for intermediate computations
        
        [sigdecon(:,ista),pred(:,ista),res,dresit,mfitit,ampit{ista},nit(ista),fighdl] = ...
          iterdecon(sig,wlet,rcc,noi,[],dt,twlet,width,dres_min,mfit_min,nit_max,nimp_max,...
          fpltit,fpltend,fpltchk);
        
        if fpltend
          ax = fighdl{2}.ax(1);
          hold(ax,'on');
          text(ax,0.05,0.85,stas(ista,:),'unit','normalized');
          hold(ax,'off');
        end
        
        nitsyn{insat,iperc} = nit;
        
      end
      
      %% Group nearest impulses from different stations into triplets, using moving searching range
      spsscale = sps/40;
      loff_max = 6*spsscale;  % maximum allowed shift is 4 samples at 40 Hz, it is also true from the synthetic distribution
      %note the output 'impindep' gives the arrival index of impulse at each station, after
      %alignment based upon the entire window 'off1i', and the last three cols are the arrival time
      %difference, NOT the true location yet!
      refsta = 1;
      [impindep,imppairf,indpair,sharp] = groupimptripdecon(sigdecon,ampit,rcc,loff_max,refsta);
      
      %note here 'impindepst' inherits the first 6 cols from 'impindep', but the last three cols
      %are adjusted from arrival time difference to the true location offset accounting for the best
      %alignment upon the entire window that is also used in grouping!
      impindep(:,7:8) = impindep(:,7:8)+repmat([off1i(2) off1i(3)],size(impindep,1),1); %account for prealignment
      impindepst = sortrows(impindep,1);
      impgrp{insat,iperc} = impindepst;

      
      %% distance, diff time for source pairs right after grouping
      %note the 'tsep' obtained from the deconvolved positive peaks should be identical to that if
      %obtained from the deconvolved impulses themselves, which represent the arrival indices of the
      %zero-crossing
      ista=1;
      %convert time offset to relative loc
      [imploc, indinput] = off2space002(impindepst(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
      [impindepstst, indsort] = sortrows(impindepst, (ista-1)*2+1);
      implocst = imploc(indsort, :);
      tarvlsplst = impindepstst(:,(ista-1)*2+1);
      
      %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
      m = 5;
      [dtarvl,dneloc,eucdist] = srcdistNtoNm(tarvlsplst, implocst, m);
      
      %project along min-scatter direction
      projang = 135;
      [projx,projy,projxy] = customprojection(implocst,projang);
      dprojxy = [diffcustom(projx,1,'forward') diffcustom(projy,1,'forward')];
      mprojx1nn1(insat,iperc) = median(abs(dprojxy(:,1))); %median consecutive dist along min-scatter after grouping
      
      %For each LFE source, get its distance to all other LFEs, maybe we don't care that long separation in time
      [~,dprojxy2all] = srcdistall(tarvlsplst,projxy,[0 2*sps]);
      mprojx12all(insat,iperc) = median(abs(dprojxy2all(:,1))); %median dist along min-scatter to all sources after grouping

      %%%if you want to plot the deconvolved sources
      if pltsrcflag1
        %plot the scatter of offsets, accounting for prealignment offset, == true offset
        xran = [-loff_max+off1i(2)-1 loff_max+off1i(2)+1];
        yran = [-loff_max+off1i(3)-1 loff_max+off1i(3)+1];
        offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
        offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
        cran = [0 lsig];
        f.fig = figure;
        f.fig.Renderer = 'painters';
        ax=gca;
        [ax,torispl,mamp] = plt_decon_imp_scatter(ax,impindepst,xran,yran,cran,offxran,offyran,...
          sps,50,'mean','tarvl');
        scatter(ax,off1i(2),off1i(3),20,'ks','filled','MarkerEdgeColor','k');
        title(ax,'Grouped Total');
        
        %plot the scatter of sources in terms of rela locations
        xran = [-4 4];
        yran = [-4 4];
        cran = [0 lsig/sps];
        f.fig = figure;
        f.fig.Renderer = 'painters';
        ax=gca;
        [ax] = plt_decon_imp_scatter_space(ax,impindepst,xran,yran,cran,offxran,...
          offyran,sps,50,ftrans,'mean','tarvl');
        plot(ax,xcut,ycut,'k-','linew',2);
        title(ax,'Grouped Total');
        
        %plot the cumulative density and summed amp of detections
        [density1d, inddup] = density_pixel(impindepst(:,7),impindepst(:,8)); %count at each unique loc
        [imploc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
        density1d = [imploc(:,1:2) density1d(:,3)];
        amp = mean(impindepst(:,[2 4 6]),2); %amp for all LFE catalog
        ampsum = sum_at_indices(amp,inddup);
        ampsum1d = sortrows([imploc(:,1:2) ampsum], 3);
        cstr = {'# detections / pixel'; 'amp sum / pixel'};
        [f] = plt_sum_pixel(density1d,ampsum1d,[-4 4],[-4 4],20,cstr,'o','linear');
        hold(f.ax(1),'on');
        plot(f.ax(1),xcut,ycut,'k-','linew',2);
        hold(f.ax(2),'on');
        plot(f.ax(2),xcut,ycut,'k-','linew',2);
        supertit(f.ax,'Grouped Total');
        
        %distribution of amp, statistically and in space
        f = initfig(10,5,1,2); %initialize fig
        ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
        plot(ax,(1:length(ampgtsum1d))/length(ampgtsum1d),ampgtsum1d(:,3),'ko-','markersize',3);
        plot(ax,(1:length(ampsum1d))/length(ampsum1d),ampsum1d(:,3),'ro-','markersize',3);
        xlabel(ax,'sorted loc index');
        ylabel(ax,'amp at pixel');
        ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
        histogram(ax,ampgtsum1d(:,3),'facecolor','k');
        histogram(ax,ampsum1d(:,3),'facecolor','r');
        plot(ax,[median(ampsum1d(:,3)) median(ampsum1d(:,3))],ax.YLim,'b--');
        ylabel(ax,'count');
        xlabel(ax,'amp at pixel');
        supertit(f.ax,'Grouped Total');

        %plot distance between srcs N&N-m along the projected direction
        if ~isempty(dprojxy)
          nsep = 1;
          [f] = plt_customprojdist(tarvlsplst,implocst,dtarvl{nsep},eucdist{nsep},...
            projxy,dprojxy,projang,sps,'tarvl');
          hold(f.ax(1),'on');
          plot(f.ax(1),xcut,ycut,'k-','linew',2);
        end
      end
      % keyboard
      
      %% Remove the small-amplitude, secondary triplets from the grouped result
      %convert the sources in terms of the arrival time of zero-crossing to positve peaks' indices
      ppkindep = impindep;
      for ista = 1: 3
        ppkindep(:,(ista-1)*2+1) = ppkindep(:,(ista-1)*2+1)+ppeaksf(ista)-zcrossesf(ista);
      end
      npkindep = impindep;  %negative peaks
      for ista = 1: 3
        npkindep(:,(ista-1)*2+1) = npkindep(:,(ista-1)*2+1)+npeaksf(ista)-zcrossesf(ista);
      end
      
      %use function 'removesecondarysrc' to remove the smaller amplitude, secondary triplets that
      %are too close in time to large triplets
      minsrcsep = 20;   % min allowed separation between deconvolved peaks
      [ppkindepsave,indremove] = removesecondarysrc(ppkindep,sigsta(:,1:3));
      %       [pkindepsave,indremove] = removesecondarysrc(pkindep,sigsta,minsrcsep);
      
      %REMOVE the secondary sources from the grouped result
      impindep(indremove, :) = [];
      ppkindep(indremove, :) = [];
      npkindep(indremove, :) = [];
      sharp(indremove, :) = [];
      impindepst = sortrows(impindep,1);
      imp{insat,iperc} = impindepst;
      
      %% distance, diff time for source pairs, after 2ndary sources removed
      %note the 'tsep' obtained from the deconvolved positive peaks should be identical to that if
      %obtained from the deconvolved impulses themselves, which represent the arrival indices of the
      %zero-crossing
      ista=1;
      %convert time offset to relative loc
      [imploc, indinput] = off2space002(impindepst(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
      [impindepstst, indsort] = sortrows(impindepst, (ista-1)*2+1);
      implocst = imploc(indsort, :);
      tarvlsplst = impindepstst(:,(ista-1)*2+1);
      
      %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
      m = 5;
      [dtarvl,dneloc,eucdist] = srcdistNtoNm(tarvlsplst, implocst, m);
      
      %project along min-scatter direction
      projang = 135;
      [projx,projy,projxy] = customprojection(implocst,projang);
      dprojxy = [diffcustom(projx,1,'forward') diffcustom(projy,1,'forward')];
      mprojx2nn1(insat,iperc) = median(abs(dprojxy(:,1))); %median consecutive dist along min-scatter after 2nd src removal
      
      %For each LFE source, get its distance to all other LFEs, maybe we don't care that long separation in time
      [~,dprojxy2all] = srcdistall(tarvlsplst,projxy,[0 2*sps]);
      mprojx22all(insat,iperc) = median(abs(dprojxy2all(:,1))); %median dist along min-scatter to all sources after 2nd src removal
      
      %% plot amp ratio 12 and 13, and 23, without 2ndary sources
      %%%ideally, when templates and data are not normalized, and there is no particular noise or any
      %%%other factors causing the amplitude scaling between temp and data and each station to be
      %%%vastly different, then for each deconvolved source, the direct impulse amp should be
      %%%~identical at all stations, ie., the ratio between station pairs should be ~1
      srcampr = [impindepst(:,2)./impindepst(:,4) impindepst(:,2)./impindepst(:,6) ...
        impindepst(:,4)./impindepst(:,6)];
      
      psrcamprs = [impindepst(:,2)*max(greenf(:,1))./(impindepst(:,4)*max(greenf(:,2))) ...
        impindepst(:,2)*max(greenf(:,1))./(impindepst(:,6)*max(greenf(:,3))) ...
        impindepst(:,4)*max(greenf(:,2))./(impindepst(:,6)*max(greenf(:,3)))];
      psrcamps = [impindepst(:,2)*max(greenf(:,1)) impindepst(:,4)*max(greenf(:,2)) ...
        impindepst(:,6)*max(greenf(:,3))];
      
      nsrcamprs = [impindepst(:,2)*min(greenf(:,1))./(impindepst(:,4)*min(greenf(:,2))) ...
        impindepst(:,2)*min(greenf(:,1))./(impindepst(:,6)*min(greenf(:,3))) ...
        impindepst(:,4)*min(greenf(:,2))./(impindepst(:,6)*min(greenf(:,3)))];
      nsrcamps = [impindepst(:,2)*min(greenf(:,1)) impindepst(:,4)*min(greenf(:,2))...
        impindepst(:,6)*min(greenf(:,3))];
      
      msrcampr{insat,iperc} = median(log10(srcampr), 1);
      madsrcampr{insat,iperc} = mad(log10(srcampr), 1, 1);
      mpsrcamprs{insat,iperc} = median(log10(psrcamprs), 1);
      madpsrcamprs{insat,iperc} = mad(log10(psrcamprs), 1, 1);
      mnsrcamprs{insat,iperc} = median(log10(nsrcamprs), 1);
      madnsrcamprs{insat,iperc} = mad(log10(nsrcamprs), 1, 1);
      nsrc{insat,iperc} = size(srcampr,1);
      
      %what is the deviation of amp ratio from the median for each source?
      lndevsrcampr = srcampr-median(srcampr, 1); % in linear scale
      lgdevsrcampr = log10(srcampr)-median(log10(srcampr), 1); % in log scale, note that log2+log5=log10, so this means a ratio
      
      srcamprall{insat,iperc} = srcampr;
      
      %%
      %%%if you want to plot the deconvolved sources
      if pltsrcflag2
        %plot the scatter of offsets, accounting for prealignment offset, == true offset
        xran = [-loff_max+off1i(2)-1 loff_max+off1i(2)+1];
        yran = [-loff_max+off1i(3)-1 loff_max+off1i(3)+1];
        offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
        offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
        cran = [0 lsig];
        f.fig = figure;
        f.fig.Renderer = 'painters';
        ax=gca;
        [ax,torispl,mamp,xbnd,ybnd] = plt_decon_imp_scatter(ax,impindepst,xran,yran,cran,offxran,offyran,...
          sps,50,'mean','tarvl');
        scatter(ax,off1i(2),off1i(3),20,'ks','filled','MarkerEdgeColor','k');
        title(ax,'Secondary sources removed');
        
        %plot the scatter of sources in terms of rela locations
        xran = [-4 4];
        yran = [-4 4];
        cran = [0 lsig/sps];
        f.fig = figure;
        f.fig.Renderer = 'painters';
        ax=gca;
        [ax] = plt_decon_imp_scatter_space(ax,impindepst,xran,yran,cran,offxran,...
          offyran,sps,50,ftrans,'mean','tarvl');
        text(ax,0.98,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
          'HorizontalAlignment','right');
        text(ax,0.02,0.05,sprintf('%.2f',mprojx22all(insat,iperc)),'Units','normalized',...
          'HorizontalAlignment','left');
        title(ax,'Secondary sources removed');
        
        %plot the cumulative density and summed amp of detections
        [density1d, inddup] = density_pixel(impindepst(:,7),impindepst(:,8)); %count at each unique loc
        [imploc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
        density1d = [imploc(:,1:2) density1d(:,3)];
        amp = mean(impindepst(:,[2 4 6]),2); %amp for all LFE catalog
        ampsum = sum_at_indices(amp,inddup);
        ampsum1d = sortrows([imploc(:,1:2) ampsum], 3);
        cstr = {'# detections / pixel'; 'amp sum / pixel'};
        [f] = plt_sum_pixel(density1d,ampsum1d,[-4 4],[-4 4],20,cstr,'o','linear');
        text(f.ax(1),0.98,0.15,sprintf('%d; %.2f located at the ground truth',...
          sum(impindepst(:,7)==2&impindepst(:,8)==2),...
          sum(impindepst(:,7)==2&impindepst(:,8)==2)/size(impindepst,1)),...
          'Units','normalized','HorizontalAlignment','right','FontSize',9);
        supertit(f.ax,'Secondary sources removed');
        
%         %distribution of amp, statistically and in space
%         f = initfig(10,5,1,2); %initialize fig
%         ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%         plot(ax,(1:length(ampgtsum1d))/length(ampgtsum1d),ampgtsum1d(:,3),'ko-','markersize',3);
%         plot(ax,(1:length(ampsum1d))/length(ampsum1d),ampsum1d(:,3),'ro-','markersize',3);
%         xlabel(ax,'sorted loc index');
%         ylabel(ax,'amp at pixel');
%         ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%         histogram(ax,ampgtsum1d(:,3),'facecolor','k');
%         histogram(ax,ampsum1d(:,3),'facecolor','r');
%         plot(ax,[median(ampsum1d(:,3)) median(ampsum1d(:,3))],ax.YLim,'b--');
%         ylabel(ax,'count');
%         xlabel(ax,'amp at pixel');
%         supertit(f.ax,'Secondary sources removed');
        
        if ~isempty(dprojxy)
          nsep = 1;
          [f] = plt_customprojdist(tarvlsplst,implocst,dtarvl{nsep},eucdist{nsep},...
            projxy,dprojxy,projang,sps,'tarvl');
        end
        
        %plot histograms of source amp
        f = initfig(12,5,1,3); %initialize fig
        f.ax(1:3) = plt_deconpk_rat_comb4th(f.ax(1:3),srcampr,impindepst,'k','hist');
        supertit(f.ax,'Secondary sources removed');
      end
      % keyboard
            
      %% 2ndary src removed, prediction of impulse tarvl at 4th sta given sources and empirical off14-src relation
      %%%carry out 'deconvolution' at 4th stations as well for the tarvl and amp
      modname = 'timeoff_plfit_4thsta_160sps.mat';
      planefit = load(strcat(rstpath, '/MAPS/',modname));
      
      pred4diff = [];
      predoff14 = [];
      for ista = 4:nsta
        wlet = greenf(:,ista);  %template here is best aligned, tapered, linear trend removed, filtered
        lwlet = length(wlet);
        sig = sigsta(:,ista); %best aligned, filtered, tapered
        lsig = length(sig);
        
        dt = 1/sps;  % sampling interval
        twlet = zcrossesf(ista)*dt;
        fpltit = 0;  % plot flag for each iteration
        fpltend = 0;  % plot flag for the final iteration
        fpltchk = 0; % plot flag for intermediate computations
        rmse = planefit.gof{ista}.rmse;
        %1.6*rmse seems proper from whole-win RCC, 2 from concatenated RCC
        offmax = round(2.0*rmse);
        
        [sigdecon(:,ista),pred,res,dresit,mfitit,ampit{ista},fighdl] = ...
          iterdecon_4thsta(sig,wlet,[],rcc1i(:,ista-3),[],...
          dt,twlet,impindep,stas(ista,:),off1i(ista),[],offmax,...
          fpltit,fpltend,fpltchk);
        
        ampiti = ampit{ista};
        impindep(:,9+(ista-4)*2+1) = ampiti(:,1);
        impindep(:,9+(ista-3)*2) = ampiti(:,2);
        if ista == nsta
          pred4diff(:,ista-3) = ampiti(:,end-1);  %difference between found peak and predicted arrival
          predoff14(:,ista-3) = ampiti(:,end);  %predicted off14
        end
      end
      
      %%%given zero-crossing indices, obtain corresponding positive and negative peak indices
      ppkindep = impindep;
      npkindep = impindep;
      for ista = 1: 3
        ppkindep(:,(ista-1)*2+1) = ppkindep(:,(ista-1)*2+1)+ppeaksf(ista)-zcrossesf(ista);
        npkindep(:,(ista-1)*2+1) = npkindep(:,(ista-1)*2+1)+npeaksf(ista)-zcrossesf(ista);
      end
      for ista = 4:nsta
        ppkindep(:,9+(ista-4)*2+1) = ppkindep(:,9+(ista-4)*2+1)+ppeaksf(ista)-zcrossesf(ista);
        npkindep(:,9+(ista-4)*2+1) = npkindep(:,9+(ista-4)*2+1)+npeaksf(ista)-zcrossesf(ista);
      end
      
      %       %%%plot the predicted tarvl and amp of decon impulses at 4th sta vs. waveform
      %       [f] = plt_deconpk_sigpk_comp_4thsta(sigsta(:,4:nsta),stas(4:nsta,:),...
      %         impindep(:,10:end),ppkindep4th,npkindep4th,greenf(:,4:nsta));
      
      %%% further ELIMINATE sources that fail the check at 4th stations
      trust4th = 4; % trust KLNB the most among all 4th stations
      indremove = find(impindep(:,9+(trust4th-4)*2+1)==0 & impindep(:,9+(trust4th-3)*2)==0);
      pred4difftr = pred4diff(setdiff(1:size(pred4diff,1),indremove),trust4th-3);
      predoff14tr = predoff14(setdiff(1:size(pred4diff,1),indremove),trust4th-3);
      impindep(indremove,:) = [];
      ppkindep(indremove, :) = [];
      npkindep(indremove, :) = [];
      impindepst = sortrows(impindep,1);
      imp4th{insat,iperc} = impindepst;
      
      %%%
      %off14 prediction for decon srcs using plane fit model WITH alignment
      [~,off14pred] = pred_tarvl_at4thsta(stas(trust4th,:),impindepst(:,7),impindepst(:,8),...
        impindepst(:,1),0);
      %off14 computed from actually-matched arrivals at stas 1 and 4 from decon
      % +repmat(off1i(trust4th),size(impindepst,1),1)
      off14 = impindepst(:,1)-impindepst(:,9+(trust4th-4)*2+1); %after prealignment
      
      %source arrival prediction from plane fitting, including calibrating the arrival prediction in the
      %context of prealigned signal, note sign is '+'
      [~,off14gt] = pred_tarvl_at4thsta(stas(trust4th,:),impgtst(:,7),impgtst(:,8),impgtst(:,1),0);
      
      %%%amp ratio 12 and 13, and 23, with 4th station checked
      %%%ideally, when templates and data are not normalized, and there is no particular noise or any
      %%%other factors causing the amplitude scaling between temp and data and each station to be
      %%%vastly different, then for each deconvolved source, the direct impulse amp should be
      %%%~identical at all stations, ie., the ratio between station pairs should be ~1
      srcampr4th = [impindepst(:,2)./impindepst(:,4) impindepst(:,2)./impindepst(:,6) ...
        impindepst(:,4)./impindepst(:,6) impindepst(:,2)./impindepst(:,end)];
      
      psrcamprs4th = [impindepst(:,2)*max(greenf(:,1))./(impindepst(:,4)*max(greenf(:,2))) ...
        impindepst(:,2)*max(greenf(:,1))./(impindepst(:,6)*max(greenf(:,3))) ...
        impindepst(:,4)*max(greenf(:,2))./(impindepst(:,6)*max(greenf(:,3))) ...
        impindepst(:,2)*max(greenf(:,1))./(impindepst(:,end)*max(greenf(:,end)))];
      psrcamps4th = [impindepst(:,2)*max(greenf(:,1)) impindepst(:,4)*max(greenf(:,2)) ...
        impindepst(:,6)*max(greenf(:,3)) impindepst(:,end)*max(greenf(:,end))];
      
      nsrcamprs4th = [impindepst(:,2)*min(greenf(:,1))./(impindepst(:,4)*min(greenf(:,2))) ...
        impindepst(:,2)*min(greenf(:,1))./(impindepst(:,6)*min(greenf(:,3))) ...
        impindepst(:,4)*min(greenf(:,2))./(impindepst(:,6)*min(greenf(:,3))) ...
        impindepst(:,2)*min(greenf(:,1))./(impindepst(:,end)*min(greenf(:,end)))];
      nsrcamps4th = [impindepst(:,2)*min(greenf(:,1)) impindepst(:,4)*min(greenf(:,2))...
        impindepst(:,6)*min(greenf(:,3)) impindepst(:,end)*min(greenf(:,end))];
      
      msrcampr4th{insat,iperc} = median(log10(srcampr4th), 1);
      madsrcampr4th{insat,iperc} = mad(log10(srcampr4th), 1, 1);
      mpsrcamprs4th{insat,iperc} = median(log10(psrcamprs4th), 1);
      madpsrcamprs4th{insat,iperc} = mad(log10(psrcamprs4th), 1, 1);
      mnsrcamprs4th{insat,iperc} = median(log10(nsrcamprs4th), 1);
      madnsrcamprs4th{insat,iperc} = mad(log10(nsrcamprs4th), 1, 1);
      nsrc4th{insat,iperc} = size(srcampr4th,1);
      
      %what is the deviation of amp ratio from the median for each source?
      lndevsrcampr = srcampr-median(srcampr, 1); % in linear scale
      lgdevsrcampr = log10(srcampr)-median(log10(srcampr), 1); % in log scale, note that log2+log5=log10, so this means a ratio
      
      srcamprall4th{insat,iperc} = srcampr4th;
            
      %% distance, diff time for source pairs, after 4th sta check
      %note the 'tsep' obtained from the deconvolved positive peaks should be identical to that if
      %obtained from the deconvolved impulses themselves, which represent the arrival indices of the
      %zero-crossing
      ista=1;
      %convert time offset to relative loc
      [imploc, indinput] = off2space002(impindepst(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
      [impindepstst, indsort] = sortrows(impindepst, (ista-1)*2+1);
      implocst = imploc(indsort, :);
      tarvlsplst = impindepstst(:,(ista-1)*2+1);
      
      %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
      m = 5;
      [dtarvl,dneloc,eucdist] = srcdistNtoNm(tarvlsplst, implocst, m);
      
      %project along min-scatter direction
      projang = 135;
      [projx,projy,projxy] = customprojection(implocst,projang);
      dprojxy = [diffcustom(projx,1,'forward') diffcustom(projy,1,'forward')];
      mprojx3nn1(insat,iperc) = median(abs(dprojxy(:,1))); %median consecutive dist along min-scatter after 2nd src removal
      
      %For each LFE source, get its distance to all other LFEs, maybe we don't care that long separation in time
      [~,dprojxy2all] = srcdistall(tarvlsplst,projxy,[0 2*sps]);
      mprojx32all(insat,iperc) = median(abs(dprojxy2all(:,1))); %median dist along min-scatter to all sources after 2nd src removal
      
      %%
      %%%if you want to plot the deconvolved sources
      if pltsrcflag3
        %plot the scatter of offsets, accounting for prealignment offset, == true offset
        xran = [-loff_max+off1i(2)-1 loff_max+off1i(2)+1];
        yran = [-loff_max+off1i(3)-1 loff_max+off1i(3)+1];
        offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
        offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
        cran = [0 lsig];
        f.fig = figure;
        f.fig.Renderer = 'painters';
        ax=gca;
        [ax,torispl,mamp,xbnd,ybnd] = plt_decon_imp_scatter(ax,impindepst,xran,yran,cran,offxran,offyran,...
          sps,50,'mean','tarvl');
        scatter(ax,off1i(2),off1i(3),20,'ks','filled','MarkerEdgeColor','k');
        title(ax,'Checkd at 4th stas');
        
        %plot the scatter of sources in terms of rela locations
        xran = [-4 4];
        yran = [-4 4];
        cran = [0 lsig/sps];
        f.fig = figure;
        f.fig.Renderer = 'painters';
        ax=gca;
        [ax] = plt_decon_imp_scatter_space(ax,impindepst,xran,yran,cran,offxran,...
          offyran,sps,50,ftrans,'mean','tarvl');
        text(ax,0.98,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
          'HorizontalAlignment','right');
        text(ax,0.02,0.05,sprintf('%.2f',mprojx32all(insat,iperc)),'Units','normalized',...
          'HorizontalAlignment','left');
        title(ax,'Checkd at 4th stas');
        
        %plot the cumulative density and summed amp of detections
        [density1d, inddup] = density_pixel(impindepst(:,7),impindepst(:,8)); %count at each unique loc
        [imploc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
        density1d = [imploc(:,1:2) density1d(:,3)];
        amp = mean(impindepst(:,[2 4 6]),2); %amp for all LFE catalog
        ampsum = sum_at_indices(amp,inddup);
        ampsum1d = sortrows([imploc(:,1:2) ampsum], 3);
        cstr = {'# detections / pixel'; 'amp sum / pixel'};
        [f] = plt_sum_pixel(density1d,ampsum1d,[-4 4],[-4 4],20,cstr,'o','linear');
        text(f.ax(1),0.98,0.15,sprintf('%d; %.2f located at the ground truth',...
          sum(impindepst(:,7)==2&impindepst(:,8)==2),...
          sum(impindepst(:,7)==2&impindepst(:,8)==2)/size(impindepst,1)),...
          'Units','normalized','HorizontalAlignment','right','FontSize',9);
        supertit(f.ax,'Checkd at 4th stas');
              
        % f = initfig(5,5,1,1); %initialize fig
        % ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
        % p1=scatter3(ax,impindepst(:,7),impindepst(:,8),off14pred,15,'r','filled','MarkerEdgeColor','k');
        % p2=scatter3(ax,impgtst(:,7),impgtst(:,8),off14gt,15,[.5 .5 .5],'filled','MarkerEdgeColor','k');
        % legend(ax,[p1 p2],'plane fit pred for decon srcs', 'plane fit pred for ground truth');
        % xlabel(ax,sprintf('off12 at %d sps',sps));
        % ylabel(ax,sprintf('off13 at %d sps',sps));
        % zlabel(ax,sprintf('off14 at %d sps',sps));
        % title(ax,stas(trust4th,:));
        % view(ax, 35, 27);
        % ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
        % p1=scatter3(ax,impindepst(:,7),impindepst(:,8),off14pred,15,'r','filled','MarkerEdgeColor','k');
        % p2=scatter3(ax,impindepst(:,7),impindepst(:,8),off14,15,[.5 .5 .5],'filled',...
        %   'MarkerEdgeColor','k'); %correct ground truth for record alignment
        % legend(ax,[p1 p2],'plane fit pred for decon srcs', 'actually-matched arrival diff.');
        % xlabel(ax,sprintf('off12 at %d sps',sps));
        % ylabel(ax,sprintf('off13 at %d sps',sps));
        % zlabel(ax,sprintf('off14 at %d sps',sps));
        % view(ax, 45, 10);
        
        %plot error in off14
        f = initfig(12,5,1,2); %initialize fig
        ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
        histogram(ax,pred4difftr);
        plot(ax,[median(pred4difftr) median(pred4difftr)],ax.YLim,'r--','LineWidth',1);
        plot(ax,[-offmax -offmax],ax.YLim,'k--');
        plot(ax,[offmax offmax],ax.YLim,'k--');
        xlabel(ax,'diff. in 4th arrival between pred and decon');
        ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
        histogram(ax,off14pred-off14);
        plot(ax,[median(off14pred-off14) median(off14pred-off14)],ax.YLim,'r--','LineWidth',1);
        plot(ax,[-offmax -offmax],ax.YLim,'k--');
        plot(ax,[offmax offmax],ax.YLim,'k--');
        xlabel(ax,'diff. in off14 between plane-fit and decon');
        ylabel(ax,'count');
        % ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
        % histogram(ax,off14pred-off14gt);
        % plot(ax,[-offmax -offmax],ax.YLim,'k--');
        % plot(ax,[offmax offmax],ax.YLim,'k--');
        % xlabel(ax,'diff. in plane-fit off14 between decon and ground truth');
        % ylabel(ax,'count');
        
        %plot distance between srcs N&N-m along the projected direction
        if ~isempty(dprojxy)
          nsep = 1;
          [f] = plt_customprojdist(tarvlsplst,implocst,dtarvl{nsep},eucdist{nsep},...
            projxy,dprojxy,projang,sps,'tarvl');
        end
        
        %plot histograms of source amp
        f = initfig(16,5,1,4); %initialize fig
        f.ax(1:4) = plt_deconpk_rat_comb4th(f.ax(1:4),srcampr4th,impindepst,'k','hist');
        supertit(f.ax,'Checkd at 4th stas');
      end
      
      % keyboard
      
    end %loop end for saturation level
    
    % keyboard
    
  end %loop end for percent of noise
  
  if tdura == 0.25
%     save(strcat('rst_synth_onespot','_td',num2str(tdura),'.mat'));
%     save(strcat('rst_synth_onespotnitm8500','_td',num2str(tdura),'.mat'));
    save(strcat('rst_synth_onespotmedwtcoef','_td',num2str(tdura),'.mat'));
  elseif tdura == 0.4 
%     save('rst_synth_onespot.mat');
%     save('rst_synth_onespotnitm8500.mat');
    save('rst_synth_onespotmedwtcoef.mat');
  end
  
end %if need to recalculate

%%%load data
savefile = 'deconv_stats4th_allbstsig.mat';
load(strcat(rstpath, '/MAPS/',savefile));

keyboard

%%
%%%loop for noise level
for iperc = 1: ntrial
  perc = perctrial(iperc);
  disp(perc);
  
  %initialize figs
  f1 = initfig(16,8,2,4,(iperc-1)*2+1); %map loc for srcs no 2ndary srcs
  supertit(f1.ax,sprintf('noise: %.1f, Secondary sources removed',perc));
  f2 = initfig(16,8,2,4,iperc*2); %map loc for srcs checked at 4th stas
  supertit(f2.ax,sprintf('noise: %.1f, Checkd at 4th stas',perc));
  
  %%%loop for saturation level
  for insat = 1: nnsat
    
    impindepst = imp{insat,iperc};
    
    %plot the scatter of sources in terms of rela locations
    % xran = [-4 4];
    % yran = [-4 4];
    % cran = [0 lsig/sps];
    % offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
    % offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
    % ax=f1.ax(insat);
    % [ax] = plt_decon_imp_scatter_space(ax,impindepst,xran,yran,cran,offxran,...
    %   offyran,sps,30,ftrans,'mean','tarvl');
    % plot(ax,xcut,ycut,'k-','linew',2);
    % text(ax,0.98,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
    %       'HorizontalAlignment','right');
    % text(ax,0.02,0.05,sprintf('%.2f',mprojx22all(insat,iperc)),'Units','normalized',...
    % 'HorizontalAlignment','left');
    
    %plot the cumulative density and summed amp of detections
    density1d = density_pixel(impindepst(:,7),impindepst(:,8));
    [imploc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
    density1d = [imploc(:,1:2) density1d(:,3)];
    xran = [-4 4];
    yran = [-4 4];
    scale = 'linear';
    ax=f1.ax(insat);
    hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
    dum = density1d;
    dum(dum(:,3)>1, :) = [];
    if strcmp(scale,'log')
      dum(:,3) = log10(dum(:,3));
    end
    scatter(ax,dum(:,1),dum(:,2),10,dum(:,3),'linew',0.2);  %, 'MarkerEdgeColor', 'w')
    dum = sortrows(density1d,3);
    dum(dum(:,3)==1, :) = [];
    if strcmp(scale,'log')
      dum(:,3) = log10(dum(:,3));
    end
    scatter(ax,dum(:,1),dum(:,2),10,dum(:,3),'filled','MarkerEdgeColor','none');
    colormap(ax,flipud(colormap(ax,'kelicol')));
    c=colorbar(ax,'SouthOutside');
    ax.CLim(2) = prctile(dum(:,3),99);
    if strcmp(scale,'log')
      c.Label.String = strcat('log_{10}(# of detections)');
    elseif strcmp(scale,'linear')
      c.Label.String = '# of detections';
    end
    axis(ax,'equal');
    axis(ax,[xran yran]);
    ax.GridLineStyle = '--';
    ax.XAxisLocation = 'top';
    text(ax,0.98,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
      'HorizontalAlignment','right');
    text(ax,0.98,0.05,sprintf('%d events',size(impindepst,1)),'Units','normalized',...
      'HorizontalAlignment','right','FontSize',9);
    text(ax,0.02,0.05,sprintf('%.2f',mprojx22all(insat,iperc)),'Units','normalized',...
      'HorizontalAlignment','left');
    text(ax,0.98,0.15,sprintf('%d; %.2f',...
      sum(impindepst(:,7)==2&impindepst(:,8)==2),...
      sum(impindepst(:,7)==2&impindepst(:,8)==2)/size(impindepst,1)),...
      'Units','normalized','HorizontalAlignment','right','FontSize',9);

    
%     impindepst = imp4th{insat,iperc};
    
    % %plot the scatter of sources in terms of rela locations
    % xran = [-4 4];
    % yran = [-4 4];
    % cran = [0 lsig/sps];
    % offxran = [-loff_max+off1i(2) loff_max+off1i(2)];
    % offyran = [-loff_max+off1i(3) loff_max+off1i(3)];
    % ax=f2.ax(insat);
    % [ax] = plt_decon_imp_scatter_space(ax,impindepst,xran,yran,cran,offxran,...
    %   offyran,sps,30,ftrans,'mean','tarvl');
    % plot(ax,xcut,ycut,'k-','linew',2);
    % text(ax,0.98,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
    %       'HorizontalAlignment','right');
    
%     %plot the cumulative density and summed amp of detections
%     density1d = density_pixel(impindepst(:,7),impindepst(:,8));
%     [imploc, ~] = off2space002(density1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
%     density1d = [imploc(:,1:2) density1d(:,3)];
%     xran = [-4 4];
%     yran = [-4 4];
%     ax=f2.ax(insat);
%     hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%     dum = density1d;
%     dum(dum(:,3)>1, :) = [];
%     if strcmp(scale,'log')
%       dum(:,3) = log10(dum(:,3));
%     end
%     scatter(ax,dum(:,1),dum(:,2),10,dum(:,3),'linew',0.2);  %, 'MarkerEdgeColor', 'w')
%     dum = sortrows(density1d,3);
%     dum(dum(:,3)==1, :) = [];
%     if strcmp(scale,'log')
%       dum(:,3) = log10(dum(:,3));
%     end
%     scatter(ax,dum(:,1),dum(:,2),10,dum(:,3),'filled','MarkerEdgeColor','none');
%     colormap(ax,flipud(colormap(ax,'kelicol')));
%     c=colorbar(ax,'SouthOutside');
%     ax.CLim(2) = prctile(dum(:,3),99);
%     if strcmp(scale,'log')
%       c.Label.String = strcat('log_{10}(# of detections)');
%     elseif strcmp(scale,'linear')
%       c.Label.String = '# of detections';
%     end
%     axis(ax,'equal');
%     axis(ax,[xran yran]);
%     ax.GridLineStyle = '--';
%     ax.XAxisLocation = 'top';
%     text(ax,0.98,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
%       'HorizontalAlignment','right');
%     text(ax,0.98,0.05,sprintf('%d events',size(impindepst,1)),'Units','normalized',...
%       'HorizontalAlignment','right','FontSize',9);
%     text(ax,0.02,0.05,sprintf('%.2f',mprojx32all(insat,iperc)),'Units','normalized',...
%       'HorizontalAlignment','left');
%     text(ax,0.98,0.15,sprintf('%d; %.2f',...
%       sum(impindepst(:,7)==2&impindepst(:,8)==2),...
%       sum(impindepst(:,7)==2&impindepst(:,8)==2)/size(impindepst,1)),...
%       'Units','normalized','HorizontalAlignment','right','FontSize',9);
    
  end %loop end for saturation level
end %loop end for noise level

%% num of detections VS saturation rate & region size
for iperc = 1: ntrial
  for insat = 1: nnsat    
    nsrcm(insat,iperc) = nsrc{insat,iperc};
    nsrc4thm(insat,iperc) = nsrc4th{insat,iperc};
    nit = nitsyn{insat,iperc};
    nitmin(insat,iperc) = min(nit);
  end
end

f = initfig(10.5,4,1,3); %initialize fig
optaxpos(f,1,3,[],[],0.06);
color = jet(ntrial);
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
for iperc = 1: ntrial
  plot(ax,log10(nsat),nitmin(:,iperc),'-o','Color',color(iperc,:),...
    'markersize',4,'MarkerFaceColor',color(iperc,:),...
    'MarkerEdgeColor',color(iperc,:),'LineWidth',1);
end
xlabel(ax,'log_{10}(Saturation)');
ylabel(ax,'Number of iterations');
% yran = [0 1];
% ylim(ax,yran);
longticks(ax,2);
hold(ax,'off');

ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
for iperc = 1: ntrial
  p(iperc) = plot(ax,log10(nsat),nsrcm(:,iperc),'-o','Color',color(iperc,:),...
    'markersize',4,'MarkerFaceColor',color(iperc,:),...
    'MarkerEdgeColor',color(iperc,:),'LineWidth',1);
  label{iperc} = sprintf('Noise=%.1f',perctrial(iperc));
end
legend(ax,p,label,'NumColumns',2,'Location','south');
% title(ax,'Secondary sources removed');
xlabel(ax,'log_{10}(Saturation)');
ylabel(ax,'Number of detections');
yran = [6e2 2.2e3];
ylim(ax,yran);
% yran = ax.YLim;
longticks(ax,2);
ax.YAxis.Exponent = 3;
hold(ax,'off');

ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
for iperc = 1: ntrial
  plot(ax,log10(nsat),nsrc4thm(:,iperc),'-o','Color',color(iperc,:),...
    'markersize',4,'MarkerFaceColor',color(iperc,:),...
    'MarkerEdgeColor',color(iperc,:),'LineWidth',1);
end
% title(ax,'Secondary sources removed');
xlabel(ax,'log_{10}(Saturation)');
ylabel(ax,'Number of detections');
% yran = [0 1];
ylim(ax,yran);
longticks(ax,2);
ax.YAxis.Exponent = 3;
hold(ax,'off');


%% summarize consecutive dist along min-scatter VS saturation rate & region size
f = initfig(8,4.5,1,2); %initialize fig
color = jet(ntrial);
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
for iperc = 1: ntrial
%   p(iperc) = plot(ax,log10(nsat),mprojx22all(:,iperc),'-o','markersize',4,'color',color(iperc,:));
  p(iperc) = plot(ax,log10(nsat),mprojx22all(:,iperc),'-','Color',color(iperc,:),'linew',1);
  scatter(ax,log10(nsat),mprojx22all(:,iperc),nsrcm(:,iperc)/50,color(iperc,:),...
    'filled','MarkerEdgeColor','k');
  label{iperc} = sprintf('Noise=%.1f',perctrial(iperc));
end
p(ntrial+1) = plot(ax,log10(nsat),0.50*ones(nnsat,1),'k--','linew',1);  %this is from data
label{ntrial+1} = 'Data';
legend(ax,p,label,'NumColumns',2,'Location','south');
% title(ax,'Secondary sources removed');
xlabel(ax,'log_{10}(Saturation)');
ylabel(ax,'Distance (km)');
yran = [0 1];
ylim(ax,yran);
longticks(ax,2);
hold(ax,'off');

ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
for iperc = 1: ntrial
%   plot(ax,log10(nsat),mprojx32all(:,iperc),'-o','markersize',4,'color',color(iperc,:));
  plot(ax,log10(nsat),mprojx32all(:,iperc),'-','Color',color(iperc,:),'linew',1);
  scatter(ax,log10(nsat),mprojx32all(:,iperc),nsrc4thm(:,iperc)/50,color(iperc,:),...
    'filled','MarkerEdgeColor','k');
end
plot(ax,log10(nsat),0.45*ones(nnsat,1),'k--','linew',1);  %this is from data
ylim(ax,yran);
% title(ax,'Checkd at 4th stas');
longticks(ax,2);
hold(ax,'off');

% ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% for iperc = 1: ntrial
%   plot(ax,log10(nsat),mprojx12all(:,iperc),'-o','markersize',4,'color',color(iperc,:));
% end
% title(ax,'Grouped Total');
% ylim(ax,yran);

stit = supertit(f.ax,'Med. dist. along min-error direc from each to all others w/i 2 s');
movev(stit,0.4);

%% Summary of amplitude ratio
% %%%%%%%%% you can choose to plot the actual histogram for each sat and noise
% %%%loop for noise level
% for iperc = 1: ntrial
%   perc = perctrial(iperc);
%   disp(perc);
%   %%%loop for saturation level
%   for insat = 1: nnsat
%     disp(nsat(insat));    
%     f3 = initfig(16,8,2,4); %plot histograms of source amp
%     supertit(f3.ax,sprintf('noise: %.1f, Secondary sources removed & Checkd at 4th stas',...
%       perc));
%     impindepst = imp{insat,iperc};
%     srcampr = srcamprall{insat,iperc};
%     f3.ax(1:3) = plt_deconpk_rat_comb(f3.ax(1:3),srcampr,impindepst,'k','hist');
%     impindepst = imp4th{insat,iperc};
%     srcampr4th = srcamprall4th{insat,iperc};
%     f3.ax(5:end) = plt_deconpk_rat_comb4th(f3.ax(5:end),srcampr4th,impindepst,'k','hist');   
%   end %loop end for saturation level
% end %loop end for noise level
% %%%%%%%%% you can choose to plot the actual histogram for each sat and noise

%%%%%%%%% or only plot the median & mad for each sat and noise
f = initfig(12,8,2,3); %initialize fig
orient(f.fig,'landscape');
optaxpos(f,2,3,[],[],0.05,0.08);
for iperc = 1: ntrial 
  label{iperc} = sprintf('Noise=%.1f',perctrial(iperc));
end
label{ntrial+1} = 'Data';
srcampralld = allbstsig.srcamprall;  %using reference from real data
for i = 1: size(srcampralld,2)
  mamprd(i,1) = median(log10(srcampralld(:,i)));
  madamprd(i,1) =  mad(log10(srcampralld(:,i)),1);
end
ref{1} = mamprd;
ref{2} = madamprd;
f=plt_deconpk_rat_stat(f,nsat,label,msrcampr,madsrcampr,ref);
% stit = supertit(f.ax,'Secondary sources removed');
% movev(stit,0.3);
ylim(f.ax(4:end),[0 0.3]);
xlim(f.ax(:),[-0.5 2]);
print(f.fig,'-dpdf','-bestfit','/home/chaosong/Pictures/agu2023s2f7.pdf');

f = initfig(15,8,2,4); %initialize fig
orient(f.fig,'landscape');
optaxpos(f,2,3,[],[],0.05,0.08);
srcampralld = allbstsig.srcampr4thall;  %using reference from real data
for i = 1: size(srcampralld,2)
  mamprd(i,1) = median(log10(srcampralld(:,i)));
  madamprd(i,1) =  mad(log10(srcampralld(:,i)),1);
end
ref{1} = mamprd;
ref{2} = madamprd;
f=plt_deconpk_rat_stat(f,nsat,label,msrcampr4th,madsrcampr4th,ref);
stit = supertit(f.ax,'Checkd at 4th stas');
movev(stit,0.3);
% ylim(f.ax(5:end),[0 0.3]);
xlim(f.ax(:),[-0.5 2]);
%%%%%%%%% or only plot the median & mad for each sat and noise





