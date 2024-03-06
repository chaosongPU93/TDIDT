% tremorexamples.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given we already have the 4-s tremor detections, extract the seismograms
% of a bunch of examples, and plot them.
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/02/27
% Last modified date:   2024/02/27
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

ttol = 35;
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.pgc002.',cutout(1:4)));
tlen = trange(:,3)-trange(:,2);
nbst = size(trange,1);

% dates and ets
%%% NOTE: 'dates' is the dates that the tremor is active at the region of interest, so that the
%%% waveform at the stations on these dates you are looking at mainly comes from the sources at the
%%% region of interest. We want to see if there is a noticable change in spectra during the burst
%%% windows on these dates
dates = unique(trange(:,1));
years = unique(floor(dates/1000));
nets = length(years);

%% load LFE catalogs
%%% load the LFE catalog of Michael Bostock, inside and outside the rectangle in 'locinterp002_4s.m'
%obtain the location of fam 002, lon0 and lat0
ftrans = 'interpchao';
% loc0 = off2space002([0 0],sps*iup,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
%format: [fam yyyy mm dd sec dx dy lon lat dep magnitude number-of-stations], 12 cols
%time should point to the peak, zero-crossing?
% bostcat = ReformBostock(loc0(3),loc0(4),0);

bosdir = ('/home/data2/chaosong/matlab/allan/BOSTOCK');
lfefnm = ('newlfeloc');
lfeloc = load(fullfile(bosdir, lfefnm));
loc0 = lfeloc(lfeloc(:,1)==2,:);
% bostcat = ReformBostock(loc0(3),loc0(2),1);
%%%if use the lumped catalog that combine unique events from 002 and 246
bostcat = load(fullfile(bosdir, '002-246_lumped.2003-2005_cull_NEW_chao'));

%all MB's LFEs of the SAME dates that I have been using
bostdayi = cell(length(dates),1);
for i = 1: length(dates)
  date = dates(i);
  year = floor(date/1000);
  jday = floor(date-year*1000);
  a = jul2dat(year,jday);

  temp = bostcat(bostcat(:,2)==year & bostcat(:,3)==a(1) & bostcat(:,4)==a(2),:);
  bostdayi{i} = sortrows(temp, 5);
  nbostdayi(i) = size(temp,1);
  bomsumday(i) = sum(temp(:,end));
  
end
bostdayia = cat(1,bostdayi{:});
bostcat = bostdayia;

%the cut-out boundary of 4-s detections
cutout = 'ellipse';
x0 = 0.2; 
y0 = 0.2;
semia = 1.75;
semib = 1.25;
angrot = 45;
[xcut, ycut] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);

%%%2022/06/29, use the same ellipse to exclude fam 047
% bnd = [x y];
bnd = [xcut ycut];
[iin,ion] = inpolygon(bostcat(:,6),bostcat(:,7),bnd(:,1),bnd(:,2));
isinbnd = iin | ion;
bostcati = bostcat(isinbnd == 1, :);
bostcato = bostcat(isinbnd ~= 1, :);
clear bostcat

%% let's just choose the random 100
ntrace = 100;
wlen = 4*sps;
rng('default');
ind = randi(size(bostcati,1),ntrace,1);
bostuse = bostcati(ind,:);

%% read daily data, break into windows of segments, plot
%filtering passband for reading data
hisig=8; % this will give a similar spectral shape between template and signal
losig=1;

% dates and ets
%%% NOTE: 'dates' is the dates that the tremor is active at the region of interest, so that the
%%% waveform at the stations on these dates you are looking at mainly comes from the sources at the
%%% region of interest. We want to see if there is a noticable change in spectra during the burst
%%% windows on these dates
bb = bostuse(:,2)*10000+bostuse(:,3)*100+bostuse(:,4);
dates = unique(bb);
years = unique(floor(dates/10000));
nets = length(years);

trace = zeros(ntrace,wlen,3,nsta);  % 'ntrace' with each being 'wlen' long, 3 components, at 'nsta' stations
traceopt = zeros(ntrace,wlen,nsta); % 'ntrace' with each being 'wlen' long, 1 opt comp., at 'nsta' stations
k = 0;
for iets = 1: nets
  % dates in each ets
  year = years(iets);
  datesets = dates(floor(dates/10000)==year);
  
  for i = 1: length(datesets)
    
    date = datesets(i);
    mo = floor((date-year*10000)/100);
    dy = date-year*10000-mo*100;
    jday = dat2jul(mo,dy,year);
    
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
    yr = num2str(year);
    
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
    bostdayi = bostuse(bostuse(:,2)==year & bostuse(:,3)==a(1) & bostuse(:,4)==a(2),:);
    bostdayi = sortrows(bostdayi, 5);
        
    for j = 1: size(bostdayi,1)
      tcnt = bostdayi(j,5);
      tst = round(tcnt*sps)-wlen/2+1;
      ted = round(tcnt*sps)+wlen/2;
      
      k = k+1;
      trace(k,:,:,:) = STANEZ(tst:ted, :, :);
      traceopt(k,:,:) = STAopt(tst:ted, 2:end);
            
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
    text(ax,0.98,0.98,'LFE','Units','normalized',...
      'HorizontalAlignment','right');
    xlabel(ax,'Time (s)');
    title(ax,sprintf('%s velocity at %s',comp(i,:),stas(ista,:)));
    ylim(ax,[0 1.1*ntrace*0.1]);
    xlim(ax,[0 ceil(wlen/sps)]);
    longticks(ax,4);
  end
  ylabel(f.ax(1),'Amplitude');
  
  orient(f.fig,'landscape');
  fname = strcat('lfeexpnez',strtrim(stas(ista,:)),'.pdf');
  print(f.fig,'-dpdf',...
  strcat('/home/chaosong/Pictures/',fname));

end


%%
% htin = 7;   % maximum height allowed is 11 inches
% widin = 8;  % maximum width allowed is 8.5 inches
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
%   text(ax,0.98,0.98,'LFE','Units','normalized',...
%     'HorizontalAlignment','right');
%   xlabel(ax,'Time (s)');
%   title(ax,sprintf('%s velocity at %s',comp(1,:),stas(i,:)));
%   ylim(ax,[0 1.1*ntrace*0.1]);
%   xlim(ax,[0 ceil(wlen/sps)]);
%   longticks(ax,4);
% end
% ylabel(f.ax(1),'Amplitude');
% 
% fname = strcat('lfeexpopt.pdf');
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
  text(ax,0.98,0.98,'LFE','Units','normalized',...
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
  text(ax,0.98,0.98,'LFE','Units','normalized',...
    'HorizontalAlignment','right');
  xlabel(ax,'Time (s)');
  title(ax,sprintf('%s velocity',comp(1,:)));
  ylim(ax,[0 1.1*ntrace/2*vsep]);
  xlim(ax,[0 ceil(wlen/sps)]);
  longticks(ax,4);
ylabel(f.ax(1),'Amplitude');

orient(f.fig,'landscape');
fname = strcat('lfeexpopt.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/chaosong/Pictures/',fname));






