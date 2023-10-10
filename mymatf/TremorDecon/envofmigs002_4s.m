% envofmigs002_4s.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is to obtain and plot the envelope of several stations other
% than PGC trio, a few hours around the migrations. To see if the 4th or 
% more stations are suitable for deconvolution as well. Through the waxing
% waning, you would see if the envelope of other station would have a similar
% shape to that of the trio stations.
% --2023/09/27, revise the code, as we now know KLNB is the only promising
%   4th sta, so env of other stations is omitted. Also I want to add the 
%   temporal density of LFEs from other fams.
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/07/25
% Last modified date:   2023/09/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
format short e   % Set the format to 5-digit floating point
clear
clc
close all

% set(0,'DefaultFigureVisible','on');
set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, resol] = pixelperinch(1);

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
  'KLNB '
  % 'LZB  '
  % 'TWKB '
  % 'MGCB '
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

ttol = 35;
ntol = 3;
tranbst = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.pgc002.',cutout(1:4)));
tlen = tranbst(:,3)-tranbst(:,2);

ttol = 1e-3*86400;
tranmig = load(strcat(rstpath, '/MAPS/migran',num2str(round(ttol)),'s.pgc002'),'w+');
tranmig(:,2:3) = tranmig(:,2:3)/3600;

% dates and ets
%%% NOTE: 'dates' is the dates that the tremor is active at the region of interest, so that the
%%% waveform at the stations on these dates you are looking at mainly comes from the sources at the
%%% region of interest. We want to see if there is a noticable change in spectra during the burst
%%% windows on these dates
dates = unique(tranbst(:,1));
years = unique(floor(dates/1000));
nets = length(years);

% figure
% histogram(tranmig(:,3)-tranmig(:,2));

%% load other catalogs
%%%load the LFE catalog of Michael Bostock, inside and outside the rectangle in 'locinterp002_4s.m'
%%%his time is approximately corrected to the postive waveform peak 
%obtain the location of fam 002, lon0 and lat0
ftrans = 'interpchao';
% loc0 = off2space002([0 0],sps*iup,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
% %format: [fam yyyy mm dd sec dx dy lon lat dep magnitude number-of-stations], 12 cols
% %time should point to the peak, zero-crossing?
% bostcat = ReformBostock(loc0(3),loc0(4),0);

bosdir = ('/home/data2/chaosong/matlab/allan/BOSTOCK');
lfefnm = ('newlfeloc');
lfeloc = load(fullfile(bosdir, lfefnm));
loc0 = lfeloc(lfeloc(:,1)==2,:);
%%%if still use the catalog separated by each fam
%this specifies which MB LFE catalog to use, 'AMR' is updated within 2-4 Hz
%while 'NEW' is within 0.5-1 Hz
flag = 'NEW';
% flag = 'AMR';
bostname = strcat('/BOSTOCK/update20230916/total_mag_detect_0000_cull_',flag,'.txt');
bostcatof = ReformBostock(loc0(3),loc0(2),1,bostname);
bostcatof = bostcatof(bostcatof(:,2)==2003 | bostcatof(:,2)==2004 | bostcatof(:,2)==2005, :);
bostcatof = bostcatof(bostcatof(:,1)~=2 & bostcatof(:,1)~=47 & bostcatof(:,1)~=246, :);
bostcatof(:,13) = mag2moment(bostcatof(:,11));

%%%if use the lumped catalog that combine unique events from 002 and 246
bostcat = load(strcat(bosdir,'/002-246_lumped.2003-2005_cull_',flag,'_chao'));

%the cut-out boundary of 4-s detections
cutout = 'ellipse';
x0 = 0.2; 
y0 = 0.2;
semia = 1.75;
semib = 1.25;
angrot = 45;
[xcut, ycut] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);

%same rectangle in 'locinterp002_4s.m'
EW = [-7 3];
NS = [-3 4];
wid = range(EW);
hgt = range(NS);
x0 = mean(EW);
y0 = mean(NS);
[x, y] = rectangle_chao(x0,y0,wid,hgt,0.01);

%%%2022/06/29, use the same ellipse to exclude fam 047
% bnd = [x y];
bnd = [xcut ycut];
[iin,ion] = inpolygon(bostcat(:,6),bostcat(:,7),bnd(:,1),bnd(:,2));
isinbnd = iin | ion;
bostcati = bostcat(isinbnd == 1, :);
% bostcato = bostcat(isinbnd ~= 1, :);
clear bostcat

%%%load the tremor catalog of John Armbruster
%%%2022/06/29, not really very useful as the detecting window could be 128-s long (not sure)
%format: [yyyy mm dd sec dx dy lon lat dep], 9 cols;
armcat = ReformArmbrusterv2(loc0(3),loc0(4),0);
[~,ind] = min(abs(armcat(:,5))+abs(armcat(:,6)));
% armcat1 = ReformArmbruster(loc0(3),loc0(4),1);
% [~,ind1] = min(abs(armcat1(:,5))+abs(armcat1(:,6)));
[iin,ion] = inpolygon(armcat(:,5),armcat(:,6),bnd(:,1),bnd(:,2));
isinbnd = iin | ion;
armcati = armcat(isinbnd == 1, :);
armcato = armcat(isinbnd ~= 1, :);
clear armcat

%% How did i not include many of MB's LFEs?
%all MB's LFEs of the SAME dates that I have been using
bostdayi = cell(length(dates),1);
for i = 1: length(dates)
  date = dates(i);
  year = floor(date/1000);
  jday = floor(date-year*1000);
  a = jul2dat(year,jday);

  temp = bostcati(bostcati(:,2)==year & bostcati(:,3)==a(1) & bostcati(:,4)==a(2),:);
  bostdayi{i} = sortrows(temp, 5);
  nbostdayi(i) = size(temp,1);
  
  temp = bostcatof(bostcatof(:,2)==year & bostcatof(:,3)==a(1) & bostcatof(:,4)==a(2),:);
  bostdayof{i} = sortrows(temp, 5);
end
bostdayia = cat(1,bostdayi{:});   %MB 002/246 lfes of the same dates I case
bostdayofa = cat(1,bostdayof{:});   %MB other fams' lfes of the same dates I case
clear bostdayof

%%

%filtering passband for reading data, confirmed by 'spectrabursts002_4s.m'
hisig=6.3; % this will give a similar spectral shape between template and signal
losig=1.8;

k = 0;  %deconvolution burst win count
n = 0;  %auto-determined migration win count

widin = 8.5;  % maximum width allowed is 8.5 inches
htin = 10.5;   % maximum height allowed is 11 inches
nrow = 4;
ncol = 1;
pltxran = [0.06 0.96]; pltyran = [0.06 0.96]; % optimal axis location
pltxsep = 0.03; pltysep = 0.04;
ifig = 0;

inbst = []; % indices of bursts fall in the local night time 10pm--5am 
idbst = []; % indices of bursts fall in the local day time
inmig = []; % indices of migrations fall in the local night time 10pm--5am 
idmig = []; % indices of migrations fall in the local day time

for i = 1: length(dates)  % dates in each ets
  date = dates(i);
  year = floor(date/1000);
  jday = floor(date-year*1000);
  a = jul2dat(year,jday);
  if a(1) == 9
    mo = 'Sep.';
  elseif a(1) == 7
    mo = 'Jul.';
  else
    mo = 'Mar.';
  end
  dy = num2str(a(2));
  yr = num2str(a(3));
  
  %%%check the meeting notes as of 2022/08/09 for details
  if year == 2003  %in 2003, target dates are outside daylight saving time, UTC-8
    night(1) = 6;  % corresponding 10 pm local time (current day-1)
    night(2) = 13; % corresponding 5 am local time (current day)
  else  %in 2004 and 2005, target dates are inside daylight saving time, UTC-7
    night(1) = 5;  
    night(2) = 12;
  end

  %Bostock's LFE catalog on the same date
  bostdayi = bostdayia(bostdayia(:,2)==year & bostdayia(:,3)==a(1) & bostdayia(:,4)==a(2),:);
  bostdayi = sortrows(bostdayi, 5);
  % bostdayo = bostcato(bostcato(:,2)==year & bostcato(:,3)==a(1) & bostcato(:,4)==a(2),:);
  % bostdayo = sortrows(bostdayo, 5);
  
  % %Armbruster's tremor catalog on the same date
  % armdayi = armcati(armcati(:,1)==year & armcati(:,2)==a(1) & armcati(:,3)==a(2),:);
  % armdayi = sortrows(armdayi, 4);
  % armdayo = armcato(armcato(:,1)==year & armcato(:,2)==a(1) & armcato(:,3)==a(2),:);
  % armdayo = sortrows(armdayo, 4);
  
  %bursts and 4-s detections of the same day
  rangetemp = tranbst(tranbst(:,1)==date, :);
  hfdayi = hfbnd(hfbnd(:,daycol)==date, :);  % inside bound of the day
  hfdayo = hfout(hfout(:,daycol)==date, :);  % outside bound of the day
  k = k+size(rangetemp,1);
  
  %migrations of the same day
  mig = tranmig(tranmig(:,1)==date, :);
  n = n+size(mig,1);
  
  tmaxi = hfdayi(:, seccol); % starting time of max power rate of half sec inside the ellipse
  tmaxo = hfdayo(:, seccol); % starting time of max power rate of half sec outside the ellipse
  tcnti = hfdayi(:, 15);  % the center of detecting win is the 15th col
  tcnto = hfdayo(:, 15);  % the center of detecting win is the 15th col
  tbosti = bostdayi(:,5); % (arrival) time of Bostock's LFE catalog inside the rectangle
  % tbosto = bostdayo(:,5); % (arrival) time of Bostock's LFE catalog outside the rectangle
  % tarmi = armdayi(:,4); % (arrival) time of Armbruster's tremor catalog inside the rectangle
  % tarmo = armdayo(:,4); % (arrival) time of Armbruster's tremor catalog outside the rectangle

  %Bostock's LFEs of other fams 
  bostdayofi = bostdayofa(bostdayofa(:,2)==year & bostdayofa(:,3)==a(1) & bostdayofa(:,4)==a(2),:);
  bostdayofi = sortrows(bostdayofi, 5);
  tbostof = bostdayofi(:,5);
  %compute the temporal density for 1-min interval
  subwlen = 60*sps;
  ovlplen = 0;
  twins = movingwins(1,86400*sps,subwlen,ovlplen,0);
  twincnt = floor(mean(twins,2));
  %find which subwin this lfe's time belongs to
  iwin = findwhichrange(tbostof*sps,twins);
  nwin = length(twins);
  tempden = zeros(nwin,1);
  for j = 1:nwin
    tempden(j) = sum(iwin == j);
  end
  tempden = [twincnt tempden];
  tempden = tempden(tempden(:,2)>0, :); 

  %read horizontal optimal and orthogonal components
  JDAY = num2zeropadstr(jday,3);
  MO=day2month(jday,year);     % EXTERNAL function, day2month, get the month of one particular date
  direc=[datapath, '/arch', yr,'/',MO,'/'];     % directory name
  prename=[direc,yr,'.',JDAY,'.00.00.00.0000.CN'];    %  path plus prefix of data file,
  %     disp(prename);
  [STAopt,STAort,~,fileflag] = rd_daily_bpdata(year,jday,prename,stas,PERMSTA,POLSTA,...
    PERMROTS,POLROTS,sps,losig,hisig,npo,npa,[],[],[],[]);
  
  if fileflag == 0    % means there are missing files
    fprintf('Day %s / %s will be omitted because of missing files. \n', yr, JDAY);
    continue    % continue to the next day
  end
 
  %%%obtain the mean envelope
  tmp = detrend(STAopt(:,2:end));
%   fractap = 20*sps/size(STAopt,1); % if fractap is >=1, n-point von Hann window is returned
%   ptstap = fractap/2*size(STAopt,1); % if fractap is >=1, n-point von Hann window is returned
%   w = tukeywin(size(STAopt,1),fractap);
%   tmp = w.* tmp;
  tmp = detrend(tmp); %detrend again for caution
  envup = envelope(tmp);

  % %%%obtain the CC for each 1-min
  % msftaddm = 0.5*sps;
  % msftadd = 20;
  % ccmin = 0.01;  % depending on the length of trace, cc could be very low
  % for j = 1:nwin
  %   envupi = envup(twins(j,1): twins(j,2), 1:3);
  %   envcc = envupi(1+msftaddm: end-msftaddm, :);
  %   ccmid = ceil(size(envcc,1)/2);
  %   ccwlen = round(size(envcc,1)-2*(msftadd+1));
  %   iup = 1;    % times of upsampling
  %   %for the whole win, ask how well can you align the noise, does not matter if the location is
  %   %the most reliable, but the resulting CC is the highest possible
  %   off1i(1) = 0;
  %   [off1i(2),off1i(3),cc] = constrained_cc_loose(envcc(:,1:3)',ccmid,...
  %     ccwlen,msftadd,ccmin,iup);
  %   %if a better alignment cannot be achieved, use 0,0
  %   if off1i(2) == msftadd+1 && off1i(3) == msftadd+1
  %     off1i(2) = 0;
  %     off1i(3) = 0;
  %   end
  %   %Align records
  %   envali = [];  % win +/-3 s, segment of interest,first 1s will be tapered 
  %   for ista = 1: 3
  %     envali(:, ista) = envupi(1+msftaddm-off1i(ista): end-msftaddm-off1i(ista), ista);
  %   end
  %   ccw12 = xcorr(envali(:,1), envali(:,2),0,'normalized');  %0-lag maximum cc based on current alignment
  %   ccw13 = xcorr(envali(:,1), envali(:,3),0,'normalized');
  %   % ccw23 = xcorr(envali(:,2), envali(:,3),0,'normalized');
  %   ccaliw(j,1) = (ccw12+ccw13)/2;
  % end
  % ccaliw = [twincnt ccaliw];
  % ccaliw = ccaliw(ccaliw(:,2)>0, :); 

  %%
  %%%start the plot, and plot the mean envelope
  ifig = ifig+1;
  f = initfig(widin,htin,nrow,ncol,ifig);
  axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

  sublen = 1;  %length of each window is 1 hr
    
  nwin = round(24/sublen);
  
  nsubtot = 0;
  isub = 0;
  
  for jj = 1: nwin
    xran = [(jj-1)*sublen jj*sublen];
    indmig = find(mig(:,2)>=xran(1) & mig(:,3)<=xran(2));
    
    indbst = find(rangetemp(:,2)/3600>=xran(1) & rangetemp(:,3)/3600<=xran(2));
    if isempty(indmig) && isempty(indbst)
      continue
    else
      isub = isub+1;
      if isub > nrow*ncol
        fname{ifig,1} = sprintf('envelope002_%dhFig%s_%d_%.1f-%.1f.pdf',...
          sublen,num2zeropadstr(ifig,2),date,losig,hisig);
        print(f.fig,'-dpdf','-fillpage',fullfile(rstpath, '/FIGS',fname{ifig,1}));
        %move on to next figure f1
        ifig = ifig+1;  %initialize another figure
        f = initfig(widin,htin,nrow,ncol,ifig);
        axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
        
        isub = isub-nrow*ncol;  %refresh

      end
        
      ax = f.ax(isub);
      hold(ax,'on');
      ax.Box = 'on';
      for ista = 1: nsta
        plot(ax,STAopt(:,1)/3600,envup(:,ista)+0.25*(ista-1));
      end
      ylim(ax,[-0.05 2]);
      xlim(ax,xran);
      longticks(ax,3);
      
      %%%plot the ranges of 4-s migrations from automated grouping
      for j = 1: size(mig,1)
        barst = mig(j,2);
        bared = mig(j,3);
        yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.99;
        plot(ax,[barst bared],[yloc yloc],'-','linew',2,'color',[0 0 .7]);
        if barst>=xran(1) && barst<=xran(2)
          yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.96;
          text(ax,barst,yloc,num2str(j+n-size(mig,1)),'FontSize',6,'HorizontalAlignment',...
            'left','color',[0 0 .7]);
        end
        
        %distinguish and note migrations in local night or day time
        if barst>=night(1) && bared<=night(2)
          inmig = [inmig, j+n-size(mig,1)];
        else
          idmig = [idmig, j+n-size(mig,1)];
        end
      end

      %%%plot all burst windows in deconvolution
      for j = 1: size(rangetemp,1)
        
        tst = rangetemp(j,2); % start and end time of bursts
        ted = rangetemp(j,3);
        
        %how many 4-s detections fall into the burst range
        indtmaxi = find(tmaxi>=tst-0.1 & tmaxi<=ted+0.1);
        
        %       %%%%Use a fixed range of time before and after the 0.5-s strongest arrival
        %       tbuffer = 3;   % buffer time to include some coherent precursor or coda, first 1s will be tapered
        %       tstbuf = tst-tbuffer; % start and end time of bursts, buffer added
        %       tedbuf = ted+tbuffer;
        %%%%Use the start and end of the 4-s detecting window
        tstbuf = min(tcnti(indtmaxi)-2);
        tedbuf = max(tcnti(indtmaxi)+2);
        tlenbuf = tedbuf-tstbuf;
        yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.93;
        plot(ax,[tstbuf tedbuf]/3600,[yloc yloc],'k-','linew',2.5);
        if tstbuf/3600>=xran(1) && tstbuf/3600<=xran(2)
          yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.9;
          text(ax,tstbuf/3600,yloc,num2str(j+k-size(rangetemp,1)),'FontSize',6,'HorizontalAlignment',...
            'left');
        end
        
        %distinguish and note bursts in local night or day time
        if tstbuf/3600>=night(1) && tedbuf/3600<=night(2)
          inbst = [inbst, j+k-size(rangetemp,1)];
        else
          idbst = [idbst, j+k-size(rangetemp,1)];
        end
        
      end

      %plot MB 002/246 LFEs of the same date
      yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.87;
      scatter(ax,tbosti/3600, yloc*ones(size(tbosti)),6,'r','filled');

      %plot the temporal density of MB other fams' LFEs of the same date, in terms of density in 1-win windows
      yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.82;
      scatter(ax,tempden(:,1)/sps/3600, yloc*ones(size(tempden(:,1))),20,tempden(:,2),'filled');%,'MarkerEdgeColor','k'
      colormap(ax,flipud(colormap(ax,'gray')));

      % %plot the mean CC in 1-win windows
      % yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.77;
      % scatter(ax,ccaliw(:,1)/sps/3600, yloc*ones(size(ccaliw(:,1))),40,ccaliw(:,2)*5,...
      %   's','filled');%,'MarkerEdgeColor','k'
      % colormap(ax,flipud(colormap(ax,'gray')));

      %label if the time segment is during the local day time or night time 6pm--6am
      if xran(1)>=night(1) && xran(2)<=night(2)
        text(ax,0.5,0.5,'N','Units','normalized','FontSize',10);
      else
        text(ax,0.5,0.5,'D','Units','normalized','FontSize',10);
      end
      
      if isub == 1
        for ista = 1: nsta
          text(ax,xran(2)-sublen*0.02,0.25*(ista-1)+0.1,strtrim(stas(ista,:)),...
          'HorizontalAlignment','right');
          if year == 2003
            title(ax,'N: local night time 10pm--5am (UTC-8); D: day');
          else
            title(ax,'N: local night time 10pm--5am (UTC-7); D: day');
          end  
        end
      end
      
      if isub == nrow
        xlabel(ax,sprintf('Time (hr) on %s %s, %s',mo,dy,yr),'FontSize',10);
        ylabel(ax,'Envelope','FontSize',10);
      end
      
    end
%      keyboard  
 
  end  

  %if the same date is finished, print it
  %print the current figure f1
  xlabel(ax,sprintf('Time (hr) on %s %s, %s',mo,dy,yr),'FontSize',10);
  ylabel(ax,'Envelope','FontSize',10);
%   keyboard  
  fname{ifig,1} = sprintf('envelope002_%dhFig%s_%d_%.1f-%.1f.pdf',...
    sublen,num2zeropadstr(ifig,2),date,losig,hisig);
  print(f.fig,'-dpdf','-fillpage',fullfile(rstpath, '/FIGS',fname{ifig,1}));

  close all

end

%% merge all figures into a single pdf file
status = system('rm -f /home/data2/chaosong/matlab/allan/data-no-resp/PGCtrio/FIGS/envelope002v2.pdf');
for i = 1:size(fname,1)
  mergenm{i} = fullfile(rstpath, '/FIGS/',fname{i});
end
append_pdfs(fullfile(rstpath, '/FIGS/envelope002v2.pdf'),mergenm);

%delete separated files
status = system('rm -f /home/data2/chaosong/matlab/allan/data-no-resp/PGCtrio/FIGS/envelope002_*.pdf');










