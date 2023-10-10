% migrelatedbursts.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is to document the correspondance between the 4-s migrations
% & the burst windows used in deconvolution. Ideally, they should show 
% similar location and temporal variation, if there is some migrating pattern
% for both. Reference migrations come from Rubin et al. 2013.
% The product is the mean envelope with 4-s tremors, automatic migration 
% windows, automic burst windows, Allan's migrations in Rubin et al. 2013.
%
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/06/30
% Last modified date:   2023/09/30
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
  'SILB '];     % determine the trio and order, here the 1st sta is PGC
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
tranbst = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(round(ttol)),'s.pgc002.',cutout(1:4)));
tlen = tranbst(:,3)-tranbst(:,2);

ttol = 1e-3*86400;
tranmig = load(strcat(rstpath, '/MAPS/migran',num2str(round(ttol)),'s.pgc002'),'w+');
tranmig(:,2:3) = tranmig(:,2:3)/3600; 

%%%migrations in Rubin et al. 2013, correspondence is listed as to my deconvolution
%%%burst windows, and to my migration windows
tranmigrubin = [
           2003062 6.15 6.28; %s2b, <->none, <->#3
           2003062 9.02 9.16; %s2c, <->#23, <->#12
           2003062 10.65 10.81; %s2d, <->#31, <->#15
           2003062 20.85 21.25; %s2e, <->#44-47, <->#23
           2003062 21.38 21.58; %s2f, <->#48-49, <->#24
           2003063 1.71 2.04; %s2g, <->#53, <->#27
           2003063 18.03 18.21; %s2h, <->#65-66, <->#37
           2004196 11.097 11.128; %s3a, <->#73, <->none
           2004196 15.57 15.74; %s3b, <->#92-93,, <->none
           2004196 19.55 19.68; %s3c, <->#99, <->#57
           2004196 20.51 20.71; %s3d, <->#100, <->#58 
           2004197 8.5 8.8; %s3e, <->#115-116, <->#69-70
           2004197 10.38 10.72; %s3f, <->#117-118, <->#71-72
           2004197 10.88 11.00; %s3g, <->#119-120, <->#73
           2004197 11.93 12.15; %s3h, <->#121, <->#74
           2005254 7.62 7.79; %fig 6a, <->#146, <->#92
           2005254 8.43 8.80; %fig 6b, <->#147-150, <->#93
           2005254 9.58 9.71; %6c, <->#151-152, <->#94
           2005254 9.77 9.98; %6d, <->#153, <->#95
           2005254 20.43 20.49; %6e, <->#164, <->#98
           2005255 1.34 1.69; %6f, <->#170, <->#100
           2005255 9.57 10.11; %6g, <->#174-175, <->#103
           2005255 16.12 16.55; %fig 4 and 6h, <->#181, <->#109
           ];

%% load other catalogs
%%% load the LFE catalog of Michael Bostock, inside and outside the rectangle in 'locinterp002_4s.m'
%obtain the location of fam 002, lon0 and lat0
loc0 = off2space002([0 0],sps*iup,'interpchao',0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
%format: [fam yyyy mm dd sec dx dy lon lat dep magnitude number-of-stations], 12 cols
%time should point to the peak, zero-crossing?
bostcat = ReformBostock(loc0(3),loc0(4),0);

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
bostcato = bostcat(isinbnd ~= 1, :);
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


%%
% dates and ets
%%% NOTE: 'dates' is the dates that the tremor is active at the region of interest, so that the
%%% waveform at the stations on these dates you are looking at mainly comes from the sources at the
%%% region of interest. We want to see if there is a noticable change in spectra during the burst
%%% windows on these dates
dates = unique(tranbst(:,1));
years = unique(floor(dates/1000));
nets = length(years);

%filtering passband for reading data, confirmed by 'spectrabursts002_4s.m'
hisig=6.3; % this will give a similar spectral shape between template and signal
losig=1.8;

k = 0;

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
  
  %Bostock's LFE catalog on the same date
  bostdayi = bostcati(bostcati(:,2)==year & bostcati(:,3)==a(1) & bostcati(:,4)==a(2),:);
  bostdayi = sortrows(bostdayi, 5);
  bostdayo = bostcato(bostcato(:,2)==year & bostcato(:,3)==a(1) & bostcato(:,4)==a(2),:);
  bostdayo = sortrows(bostdayo, 5);
  
  %Armbruster's tremor catalog on the same date
  armdayi = armcati(armcati(:,1)==year & armcati(:,2)==a(1) & armcati(:,3)==a(2),:);
  armdayi = sortrows(armdayi, 4);
  armdayo = armcato(armcato(:,1)==year & armcato(:,2)==a(1) & armcato(:,3)==a(2),:);
  armdayo = sortrows(armdayo, 4);
  
  %bursts and 4-s detections of the same day
  rangetemp = tranbst(tranbst(:,1)==date, :);
  hfdayi = hfbnd(hfbnd(:,daycol)==date, :);  % inside bound of the day
  hfdayo = hfout(hfout(:,daycol)==date, :);  % outside bound of the day
  k = k+size(rangetemp,1);
  
  %migrations of the same day
  migrubin = tranmigrubin(tranmigrubin(:,1)==date, :);
  mig = tranmig(tranmig(:,1)==date, :);
  
  tmaxi = hfdayi(:, seccol); % starting time of max power rate of half sec inside the ellipse
  tmaxo = hfdayo(:, seccol); % starting time of max power rate of half sec outside the ellipse
  tcnti = hfdayi(:, 15);  % the center of detecting win is the 15th col
  tcnto = hfdayo(:, 15);  % the center of detecting win is the 15th col
  tbosti = bostdayi(:,5); % (arrival) time of Bostock's LFE catalog inside the rectangle
  tbosto = bostdayo(:,5); % (arrival) time of Bostock's LFE catalog outside the rectangle
  tarmi = armdayi(:,4); % (arrival) time of Armbruster's tremor catalog inside the rectangle
  tarmo = armdayo(:,4); % (arrival) time of Armbruster's tremor catalog outside the rectangle
    
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
  tmp = detrend(STAopt(:,2:4));
%   fractap = 20*sps/size(STAopt,1); % if fractap is >=1, n-point von Hann window is returned
%   ptstap = fractap/2*size(STAopt,1); % if fractap is >=1, n-point von Hann window is returned
%   w = tukeywin(size(STAopt,1),fractap);
%   tmp = w.* tmp;
  tmp = detrend(tmp); %detrend again for caution
  [envup,~] = envelope(tmp);
  menv = mean(envup,2);
  
  %%%start the plot, and plot the mean envelope
  widin = 10.5;  % maximum width allowed is 8.5 inches
  htin = 6;   % maximum height allowed is 11 inches
  nrow = 6;
  ncol = 1;
  f = initfig(widin,htin,nrow,ncol);

  pltxran = [0.06 0.96]; pltyran = [0.06 0.96]; % optimal axis location
  pltxsep = 0.03; pltysep = 0.04;  
  axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);


  % %get the locations for each axis
  % axpos = [0.08 0.87 0.85 0.1;
  %          0.08 0.74 0.85 0.1;
  %          0.08 0.585 0.85 0.1;
  %          0.08 0.1 0.8 0.23;
  %          ];
  % for isub = 1:nrow*ncol
  %   set(f.ax(isub), 'position', axpos(isub,:));
  % end

  for isub = 1: nrow
    ax = f.ax(isub);
    hold(ax,'on');
    ax.Box = 'on';
    plot(ax,STAopt(:,1)/3600,menv,'g-');
    ylim(ax,[0 1]);
    xran = [(isub-1)*24/nrow isub*24/nrow];
    xlim(ax,xran);
  
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
      yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.65;
      plot(ax,[tstbuf tedbuf]/3600,[yloc yloc],'k-','linew',2.5);
      if tstbuf/3600>=xran(1) && tedbuf/3600<=xran(2)
        yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.5;
        text(ax,tstbuf/3600,yloc,num2str(j+k-size(rangetemp,1)),'FontSize',6,'HorizontalAlignment',...
          'left');
      end
    end

    %%%plot all 4-s detection windows in-bound
    ind = find(tcnti/3600>=xran(1) & tcnti/3600<=xran(2));
    for j = 1: length(ind)
      barst = tcnti(ind(j))-2;
      bared = tcnti(ind(j))+2;
      yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.85;
      plot(ax,[barst bared]/3600,[yloc yloc],'-','linew',2.5,'color','b');
    end
    
    %%%plot all 4-s detection windows out-bound
    ind = find(tcnto/3600>=xran(1) & tcnto/3600<=xran(2));
    for j = 1: length(ind)
      barst = tcnto(ind(j))-2;
      bared = tcnto(ind(j))+2;
      yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.85;
      plot(ax,[barst bared]/3600,[yloc yloc],'-','linew',2.5,'color',[.4 .4 .4]);
    end
  
    %%%plot the ranges of confirmed 4-s migrations from Rubin et al. 2013
    for j = 1: size(migrubin,1)
      barst = migrubin(j,2);
      bared = migrubin(j,3);
      yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.75;
      plot(ax,[barst bared],[yloc yloc],'-','linew',2,'color','r');
    end
    
    %%%plot the ranges of 4-s migrations from automated grouping 
    for j = 1: size(mig,1)
      barst = mig(j,2);
      bared = mig(j,3);
      yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.7;
      plot(ax,[barst bared],[yloc yloc],'-','linew',1.5,'color',[0 0 .7]);
    end
    
    %%%plot the bostock's LFE catalog outside rectangle
    yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.95;
    ind = find(tbosto/3600>=xran(1) & tbosto/3600<=xran(2));
    scatter(ax,tbosto(ind)/3600, yloc*ones(size(tbosto(ind))),1,'m'); % bostock
    %%%plot the bostock's LFE catalog inside rectangle
    ind = find(tbosti/3600>=xran(1) & tbosti/3600<=xran(2));
    scatter(ax,tbosti(ind)/3600, yloc*ones(size(tbosti(ind))),1,'m','filled');

    if isub == 1
      text(ax,0.02,0.9,sprintf('%d',size(hfdayi,1)),'Units','normalized',...
        'HorizontalAlignment','left','color','b','FontSize',8);
      text(ax,0.02,0.75,sprintf('%d',size(hfdayo,1)),'Units','normalized',...
        'HorizontalAlignment','left','color',[.4 .4 .4],'FontSize',8);
      text(ax,0.02,0.6,sprintf('%d',size(migrubin,1)),'Units','normalized',...
        'HorizontalAlignment','left','color','r','FontSize',8);
      text(ax,0.02,0.45,sprintf('%d',size(rangetemp,1)),'Units','normalized',...
        'HorizontalAlignment','left','FontSize',8);
    end
    
    if isub == nrow
      xlabel(ax,sprintf('Time (hr) on %s %s %s',dy,mo,yr),'FontSize',12);
      ylabel(ax,'Mean envelope','FontSize',12);
    end
  end
  
  % keyboard
  %save figure
  fignm{i,1} = sprintf('meanenv%d.pdf',date);
  orient(f.fig,'landscape');
  print(f.fig,'-dpdf','-bestfit',fullfile(rstpath, '/FIGS',fignm{i,1}));
end  


%% merge all figures into a single pdf file
status = system('rm -f /home/data2/chaosong/matlab/allan/data-no-resp/PGCtrio/FIGS/meanenv.pdf');
for i = 1:size(fignm,1)
  mergenm{i} = fullfile(rstpath, '/FIGS/',fignm{i});
end
append_pdfs(fullfile(rstpath, '/FIGS/meanenv.pdf'),mergenm);

%delete separated files
status = system('rm -f /home/data2/chaosong/matlab/allan/data-no-resp/PGCtrio/FIGS/meanenv?.pdf');

keyboard

%% which burst windows have reliable/less ambiguious pre-event noise
gp{1} = 1;
gp{2} = 9:10;
gp{3} = 11;
gp{4} = 14:15;
gp{5} = 16;
gp{6} = 17;
gp{7} = 18:19;
gp{8} = 20:22;
gp{9} = 20:22;
gp{10} = 23;
gp{11} = 25:26;
gp{12} = 27:28;
gp{13} = 29:30;
gp{14} = 31;
gp{15} = 32;
gp{16} = 33;
gp{17} = 34;
gp{18} = 35;
gp{19} = 36;
gp{20} = 37;
gp{21} = 38;
gp{22} = 45:47;
gp{23} = 48:49;
gp{24} = 50;
gp{25} = 51:52;
gp{26} = 53;
gp{27} = 54;
gp{28} = 56:57;
gp{29} = 63:64;
gp{30} = 65:67;
gp{31} = 68:69;
gp{32} = 75:78;
gp{33} = 80;
gp{34} = 81:82;
gp{35} = 83:84;
gp{36} = 85;
gp{37} = 90:91;
gp{38} = 92:93;
gp{39} = 94;
gp{40} = 95:96;
gp{41} = 97;
gp{42} = 98;
gp{43} = 99;
gp{44} = 100;
gp{45} = 101;
gp{46} = 103;
gp{47} = 104:105;
gp{48} = 104:105;
gp{49} = 106:107;
gp{50} = 108:109;
gp{51} = 110:111;
gp{52} = 112:113;
gp{53} = 114;
gp{54} = 115:116;
gp{55} = 117:118;
gp{56} = 119:120;
gp{57} = 123;
gp{58} = 124:127;
gp{59} = 130;
gp{60} = 138;
gp{61} = 142;
gp{62} = 143:144;
gp{63} = 145;
gp{64} = 146;
gp{65} = 147:150;
gp{66} = 151:152;
gp{67} = 153;
gp{68} = 154;
gp{69} = 155:159;
gp{70} = 160:162;
gp{71} = 164;
gp{72} = 165:167; %?
gp{73} = 170;
gp{74} = 171;
gp{75} = 173;
gp{76} = 174:175;
gp{77} = 176:178;
gp{78} = 179:180;
gp{79} = 181;
gp{80} = 182:183;
gp{81} = 184:185;
gp{82} = 186:187;
gp{83} = 188;
gp{84} = 189:191;
gp{85} = 194:195;


























