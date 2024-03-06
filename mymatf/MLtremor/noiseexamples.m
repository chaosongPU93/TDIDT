% noiseexamples.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given we already have the 4-s tremor detections and the associated tremor
% bursts, target the periods where no tremors were found, then extract the 
% seismograms as noise. Show a bunch of examples, and plot them.
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/02/28
% Last modified date:   2024/02/28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
format short e   % Set the format to 5-digit floating point
clear
clc
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');

freqflag='hf';  % flag to indicate whether to do hf or lf;

FLAG = 'PGC'; % detector
  
fam = '002';   % family number

[timoffrot,~] = GetDays4Stack(fam);
nday = size(timoffrot, 1);

%Use new LFE catalog
CATA = 'new';

% get permanent and polaris station rotation parameters, based on 40-sps data
sft2=0;     % centroid shift of station 2
sft3=0;     % centroid shift of station 3
[PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,sft2,sft3);
if ~isequal(CATA, 'fixed')
  reftime = PERMROTS(1,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
  PERMROTS(:,4) = PERMROTS(:,4)-reftime;    % to make sure that 1st station is 0
  POLROTS(:,4) = POLROTS(:,4)-reftime;
end
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;     % convert 2-3 columns to rad
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
  'LZB'];
POLSTA=['SSIB '           % polaris station names
  'SILB '
  'KLNB '
  'MGCB '
  'TWKB '];

stas=['PGC  '
  'SSIB '
  'SILB '
  ];     % determine the trio and order, here the 1st sta is PGC
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
cntcol = 15;
maxcol = 16;  % 2 choices here, one is the center of detecting window, other is the start of strongest .5 s 
hfbnd = sortrows(hfbnd, [daycol, maxcol]);
off12ran = minmax(hfbnd(:,9)');
off13ran = minmax(hfbnd(:,10)');

%load all detections, an output of 'locinterp002_4s.m'
hfall = load(strcat(rstpath, '/MAPS/eloc.pgc.tdec_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',...
            num2str(winlen/sps),'s',num2str(sps),'sps4add_'));
hfall = sortrows(hfall, [daycol, maxcol]);
%get ones outside your boundary
hfout = setdiff(hfall,hfbnd,'rows');           
hfout = sortrows(hfout, [daycol, maxcol]);
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

% ttol = 35;
% trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.pgc002.',cutout(1:4)));

savefile = 'deconv_stats4th_allbstsig.mat';
load(strcat(rstpath, '/MAPS/',savefile));
trange = allbstsig.trangenew;

tlen = trange(:,3)-trange(:,2);
nbst = size(trange,1);

%% get the range of noise time windows
dates = unique(trange(:,1));
trangenoi = [];
for i=1:length(dates)
  tmp=trange(trange(:,1)==dates(i), :);
  aa=[];
  aa(:,1)=dates(i)*ones(length(tmp)+1, 1);
  aa(:,2:3)=[1 tmp(1,2); tmp(1:end-1,3) tmp(2:end,2); tmp(end,3) 86400];
  trangenoi = [trangenoi; aa];
end

%% generate random windows with fixed length
wlen = 4*sps;
ovlplen=0;
windows=[];
for i = 1:length(trangenoi)
  indst = trangenoi(i,2)*sps+1;
  inded = trangenoi(i,3)*sps;
  win = movingwins(indst,inded,wlen,ovlplen,0);
  win = [trangenoi(i,1)*ones(length(win), 1) win];
  windows = [windows; win];
end

%%
ntrace = 100;
rng('default');
ind = randi(length(windows),ntrace,1);

tnoiuse = windows(ind,:);

%% read daily data, break into windows of segments, plot
%filtering passband for reading data
hisig=8; % this will give a similar spectral shape between template and signal
losig=1;

% dates and ets
%%% NOTE: 'dates' is the dates that the tremor is active at the region of interest, so that the
%%% waveform at the stations on these dates you are looking at mainly comes from the sources at the
%%% region of interest. We want to see if there is a noticable change in spectra during the burst
%%% windows on these dates
dates = unique(tnoiuse(:,1));
years = unique(floor(dates/1000));
nets = length(years);

trace = zeros(ntrace,wlen,3,nsta);  % 'ntrace' with each being 'wlen' long, 3 components, at 'nsta' stations
traceopt = zeros(ntrace,wlen,nsta); % 'ntrace' with each being 'wlen' long, 1 opt comp., at 'nsta' stations
k = 0;
for iets = 1: nets
  % dates in each ets
  year = years(iets);
  datesets = dates(floor(dates/1000)==year);
  
  for i = 1: length(datesets)
    
    date = datesets(i);
    jday = floor(date-year*1000);
    a = jul2dat(year,jday);
    mo = a(1);
    if mo == 9
      mo = 'Sep.';
    elseif mo == 7
      mo = 'Jul.';
    else
      mo = 'Mar.';
    end
    dy = num2str(a(2));
    yr = num2str(a(3));
    
    %read horizontal optimal and orthogonal components
    JDAY = num2zeropadstr(jday,3);
    MO=day2month(jday,year);     % EXTERNAL function, day2month, get the month of one particular date
    direc=[datapath, '/arch', yr,'/',MO,'/'];     % directory name
    datafnm=[direc,yr,'.',JDAY,'.00.00.00.0000.CN'];    %  path plus prefix of data file,
    disp(datafnm);
    [STAopt,~,~,fileflag] = rd_daily_bpdata(year,jday,datafnm,stas,PERMSTA,POLSTA,...
      PERMROTS,POLROTS,sps,losig,hisig,npo,npa,[],[],[],[]);
    if fileflag == 0    % means there are missing files
        fprintf('Day %s / %s will be omitted because of missing files. \n', YEAR, JDAY);
        continue    % continue to the next day
    end

    
    %read the raw N, E, Z data
    if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
        POLSTA(3,:)='KELB ';
    else
        POLSTA(3,:)='KLNB ';  % remember to change it back
    end
    
    fileflag = 1;   % 1 means all files exist, the file status is normal
    
    %%% loop for each station for reading data
    for ista = 1: nsta
        found = 0;
        [LIA, idx] = ismember(stas(ista, :), PERMSTA, 'rows');
        
        %%% if station is a permanent one
        if LIA
            found = found+ LIA;
            if strcmp(PERMSTA(idx, 1:3),'PGC')
                fact = 1.0e-3;
            elseif strcmp(PERMSTA(idx, 1:3),'LZB')
                fact = 1.8e-3;
            end
            fname = strcat(datafnm,'.',PERMSTA(idx,1:3),'..BHE.D.SAC');
            if isfile(fname)    % if have the data file
                % 1. this is for data without removing station response
                %   [opt,ort,timeperm]=readperm_1daynofilterv2(datafnm, PERMSTA, PERMROTS, idx, sps, fact);
                % 2. this is for data with no response
                [dataE, dataN, ~] = readperms_norotsv2(datafnm, PERMSTA, idx, sps, fact);
                [dataZ] = readperms_norotsv2(datafnm, PERMSTA, idx, sps, fact, 'Z');
            else
                fileflag = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for station %s in day %s %s, this day will be omitted. \n',...
                        PERMSTA(idx,1:3), YEAR, JDAY);
                break   % break the entire station loop
            end
            
        end
        
        %%% if station is a polaris one
        [LIA, idx] = ismember(stas(ista, :), POLSTA, 'rows');
        if LIA
            found = found+ LIA;
            if year==2003 && jday<213 && ~strcmp(POLSTA(idx,1:4),'KELB')      % should result from some criteria
              %     if year==2003 && jday<213     % should result from some criteria
              fact=7.5e-3;
            elseif year==2003 && jday<213 && strcmp(POLSTA(idx,1:4),'KELB')      % should result from some criteria
              fact=5.0e-3;
            else
              fact=1.5e-3;
            end

            fname = strcat(datafnm,'.',POLSTA(idx,1:4),'..HHE.D.SAC');
            if isfile(fname)    % if have the data file
                %                 % 1. this is for data without removing station response
                %                 [opt,ort,nzeros]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
                % 2. this is for data with no response
                [dataE, dataN, ~] = readpols_norotsv2(datafnm, POLSTA, idx, sps, fact);
                [dataZ] = readpols_norotsv2(datafnm, POLSTA, idx, sps, fact, 'Z');
            else
                fileflag = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for station %s in day %s / %s. \n',POLSTA(idx,1:4), YEAR, JDAY);
                break   % break the entire station loop
            end
            
        end
        
        STANEZ(:, 1, ista) = Bandpass(dataN,sps,losig,hisig,npo,npa,'butter');
        STANEZ(:, 2, ista) = Bandpass(dataE,sps,losig,hisig,npo,npa,'butter');
        STANEZ(:, 3, ista) = Bandpass(dataZ,sps,losig,hisig,npo,npa,'butter');
%         STAE(:, ista) = Bandpass(dataE,sps,lo,hi,npo,npa,'butter');
%         STAZ(:, ista) = Bandpass(dataZ,sps,lo,hi,npo,npa,'butter');
    end
    
    if fileflag == 0    % means there are missing files
        fprintf('Day %s / %s will be omitted because of missing files. \n', YEAR, JDAY);
        continue    % continue to the next day
    end
    
    %4-s detections of the same day
    tnoiusei = tnoiuse(tnoiuse(:,1)==datesets(i), :);
    
    for j = 1: size(tnoiusei,1)
      tst = tnoiusei(j,2);
      ted = tnoiusei(j,3);
      
      k = k+1;
      trace(k,:,:,:) = STANEZ(tst:ted, :, :);
      traceopt(k,:,:) = STAopt(tst:ted, 2:end);
      
%         trace1(k,:,:,ista) = [STAN(tst:ted, 1) STAE(tst:ted, 1) STAZ(tst:ted, 1)];
%         trace2(k,:,:) = [STAN(tst:ted, 2) STAE(tst:ted, 2) STAZ(tst:ted, 2)];
%         trace3(k,:,:) = [STAN(tst:ted, 3) STAE(tst:ted, 3) STAZ(tst:ted, 3)];
%       end
      
    end
      
  end  

end

%%
htin = 8.5;   % maximum height allowed is 11 inches
widin = 10.5;  % maximum width allowed is 8.5 inches
nrow = 1;
ncol = 3;

for ista = 1: nsta

  f = initfig(widin,htin,nrow,ncol); %initialize fig
  figxran = [0.06 0.96]; figyran = [0.06 0.96];
  figxsep = 0.05; figysep = 0.06;
  optaxpos(f,nrow,ncol,figxran,figyran,figxsep,figysep);
  comp = ['N';'E';'Z'];
  for i = 1: 3
    ax=f.ax(i); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
    for j = 1: ntrace
      plot(ax,(1:wlen)/sps,trace(j,:,i,ista)+j*0.1,'k-');
    end
    text(ax,0.02,0.98,sprintf('%.1f--%.1f Hz',losig,hisig),'Units','normalized',...
      'HorizontalAlignment','left');
    text(ax,0.98,0.98,'Noise','Units','normalized',...
      'HorizontalAlignment','right');
    xlabel(ax,'Time (s)');
    title(ax,sprintf('%s velocity at %s',comp(i,:),stas(ista,:)));
    ylim(ax,[0 1.1*ntrace*0.1]);
    xlim(ax,[0 ceil(wlen/sps)]);
    longticks(ax,4);
  end
  ylabel(f.ax(1),'Amplitude');

  orient(f.fig,'landscape');
  fname = strcat('noiseexpnez',strtrim(stas(ista,:)),'.pdf');
  print(f.fig,'-dpdf',...
  strcat('/home/chaosong/Pictures/',fname));

end


%%
% htin = 8.5;   % maximum height allowed is 11 inches
% widin = 10.5;  % maximum width allowed is 8.5 inches
% nrow = 1;
% ncol = 3;
% 
% f = initfig(widin,htin,nrow,ncol); %initialize fig
% figxran = [0.06 0.96]; figyran = [0.06 0.96];
% figxsep = 0.05; figysep = 0.06;
% optaxpos(f,nrow,ncol,figxran,figyran,figxsep,figysep);
% comp = ['optimal'];
% 
% for i = 1: nsta
%   ax=f.ax(i); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
%   for j = 1: ntrace
%     plot(ax,(1:wlen)/sps,traceopt(j,:,i)+j*0.1,'k-');
%   end
%   text(ax,0.02,0.98,sprintf('%.1f--%.1f Hz',losig,hisig),'Units','normalized',...
%     'HorizontalAlignment','left');
%   text(ax,0.98,0.98,'Noise','Units','normalized',...
%     'HorizontalAlignment','right');
%   xlabel(ax,'Time (s)');
%   title(ax,sprintf('%s velocity at %s',comp(1,:),stas(i,:)));
%   ylim(ax,[0 1.1*ntrace*0.1]);
%   xlim(ax,[0 ceil(wlen/sps)]);
%   longticks(ax,4);
% end
% ylabel(f.ax(1),'Amplitude');
% 
% fname = strcat('noiseexpopt.pdf');
% print(f.fig,'-dpdf',...
%   strcat('/home/chaosong/Pictures/',fname));


%%
htin = 8.5;   % maximum height allowed is 11 inches
widin = 7;  % maximum width allowed is 8.5 inches
nrow = 1;
ncol = 2;

f = initfig(widin,htin,nrow,ncol); %initialize fig
figxran = [0.06 0.96]; figyran = [0.06 0.96];
figxsep = 0.05; figysep = 0.06;
optaxpos(f,nrow,ncol,figxran,figyran,figxsep,figysep);
comp = ['optimal'];
vsep = 0.4;
ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
  for j = 1: ntrace/2
    plot(ax,(1:wlen)/sps,traceopt(j,:,1)+j*vsep,'r-');
    plot(ax,(1:wlen)/sps,traceopt(j,:,2)+j*vsep,'b-');
    plot(ax,(1:wlen)/sps,traceopt(j,:,3)+j*vsep,'k-');
  end
  text(ax,0.02,0.98,sprintf('%.1f--%.1f Hz',losig,hisig),'Units','normalized',...
    'HorizontalAlignment','left');
  text(ax,0.98,0.98,'Noise','Units','normalized',...
    'HorizontalAlignment','right');
  xlabel(ax,'Time (s)');
  title(ax,sprintf('%s velocity',comp(1,:)));
  ylim(ax,[0 1.1*ntrace/2*vsep]);
  xlim(ax,[0 ceil(wlen/sps)]);
  longticks(ax,4);
ylabel(f.ax(1),'Amplitude');

ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
  for j = 1: ntrace/2
    plot(ax,(1:wlen)/sps,traceopt(j+ntrace/2,:,1)+j*vsep,'r-');
    plot(ax,(1:wlen)/sps,traceopt(j+ntrace/2,:,2)+j*vsep,'b-');
    plot(ax,(1:wlen)/sps,traceopt(j+ntrace/2,:,3)+j*vsep,'k-');
  end
  text(ax,0.02,0.98,sprintf('%.1f--%.1f Hz',losig,hisig),'Units','normalized',...
    'HorizontalAlignment','left');
  text(ax,0.98,0.98,'Noise','Units','normalized',...
    'HorizontalAlignment','right');
  xlabel(ax,'Time (s)');
  title(ax,sprintf('%s velocity',comp(1,:)));
  ylim(ax,[0 1.1*ntrace/2*vsep]);
  xlim(ax,[0 ceil(wlen/sps)]);
  longticks(ax,4);
ylabel(f.ax(1),'Amplitude');

orient(f.fig,'landscape');
fname = strcat('noiseexpopt.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/chaosong/Pictures/',fname));


