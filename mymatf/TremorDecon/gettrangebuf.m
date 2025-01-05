% gettrangebuf.m
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % obtain the actually used "trange" that includes the date, starting and ending
  % times of the 196 burts.
  %
  %
  %
  % Chao Song, chaosong@princeton.edu
  % First created date:   2024/12/18
  % Last modified date:   2024/12/18
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  format short e   % Set the format to 5-digit floating point
  clear
  clc
  close all
  
  %% for easy testing
  defval('idxbst',1:196); %global indices of bursts to run 
  defval('normflag',0); %whether to normalize templates
  defval('noiseflag',0);  %whether to use synthetic noises
  defval('pltflag',0);  %whether to plot figs for each burst
  defval('rccmwsec',0.5); %moving win len in sec for computing RCC
  defval('alignflag',1); %align the short wins between trio stations 
  
  rccflag = 1; %1 means RCC weighting is used
  whichrcc = 0; %if rcc weighting, which pair is used, 0 is best 2 pairs; 1 is 12; 2 is 13; 3 is 23
  
  %Choice to make upon the actual-used alignment at 4th stations
  if noiseflag
    align14flag = 0;  %do NOT align sta 4 wrt. 1 if using noise
  else
    align14flag = 1; 
  end
  
  %% Initialization
  %%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
  %%% AND if using the same family, same station trio
  % set(0,'DefaultFigureVisible','on');
  % set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
  
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
    'LZB  '
    'TWKB '
    'MGCB '
    'KLNB ']; % determine the trio and order, here the 1st sta is PGC
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
  trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.',...
    num2str(ntol),'.pgc002.',cutout(1:4)));
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


%%%load the empirical model param of time offset 14 for each addtional stations
off14mod = load(strcat(rstpath, '/MAPS/timeoff_planefitparam_4thsta_160sps'),'w+');

%% load other catalogs
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

% %if you choose only fam 002 regardless
% bostcati = bostcati(bostcati(:,1)==2,:);
%convert moment mag to moment, Mw = (2/3)*log_10(M0)-10.7, where M0 has the unit of dyne.cm
%(10^-7 N.m), so Mw = (2/3)*log_10(M0)-6 if M0 has the unit of N.m
bostcati(:,13) = 10.^(1.5*(6+bostcati(:,11)));

% bostcato = bostcato(bostcato(:,1)~=2 & bostcato(:,1)~=47 & bostcato(:,1)~=246,:);
bostcato(:,13) = 10.^(1.5*(6+bostcato(:,11)));


%%%load the tremor catalog of John Armbruster, 
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
  

trangenew = trange; 
for iii = 1: length(idxbst)
  
  [iets,i,j] = indofburst(trange,idxbst(iii));
  
% for iets = 3: nets
  % dates in each ets
  year = years(iets);
  datesets = dates(floor(dates/1000)==year);
    
%   for i = 2: length(datesets)
    
    date = datesets(i);
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
    rangetemp = trange(trange(:,1)==datesets(i), :);
    hfdayi = hfbnd(hfbnd(:,daycol)==datesets(i), :);  % inside bound of the day
    hfdayo = hfout(hfout(:,daycol)==datesets(i), :);  % outside bound of the day

    k = idxbst(iii);
    disp(k);

    tmaxi = hfdayi(:, seccol); % starting time of max power rate of half sec inside the ellipse
    tmaxo = hfdayo(:, seccol); % starting time of max power rate of half sec outside the ellipse  
    tcnti = hfdayi(:, 15);  % the center of detecting win is the 15th col
    tcnto = hfdayo(:, 15);  % the center of detecting win is the 15th col
    tbosti = bostdayi(:,5); % (arrival) time of Bostock's LFE catalog inside the rectangle
    tbosto = bostdayo(:,5); % (arrival) time of Bostock's LFE catalog outside the rectangle
    tarmi = armdayi(:,4); % (arrival) time of Armbruster's tremor catalog inside the rectangle
    tarmo = armdayo(:,4); % (arrival) time of Armbruster's tremor catalog outside the rectangle
    
    tst = rangetemp(j,2); % start and end time of bursts
    ted = rangetemp(j,3);

    %how many 4-s detections fall into the burst range 
    indtmaxi = find(tmaxi>=tst-0.1 & tmaxi<=ted+0.1);
    ninbst(k,1) = length(indtmaxi); %+-0.1 s in case of resolution while saving to file
    
%       %%%%Use a fixed range of time before and after the 0.5-s strongest arrival
%       tbuffer = 3;   % buffer time to include some coherent precursor or coda, first 1s will be tapered 
%       tstbuf = tst-tbuffer; % start and end time of bursts, buffer added
%       tedbuf = ted+tbuffer;
    %%%%Use the start and end of the 4-s detecting window
    tstbuf = min(tcnti(indtmaxi)-2);
    tedbuf = max(tcnti(indtmaxi)+2); 
    tlenbuf = tedbuf-tstbuf;
    trangenew(iii,2)=tstbuf;
    trangenew(iii,3)=tedbuf;

end

fid = fopen(strcat(rstpath, '/MAPS/tdec.bstranbuf',...
  num2str(ttol),'s.pgc002.',cutout(1:4)),'w+');
fprintf(fid,'%d %9.2f %9.2f \n',trangenew');
fclose(fid);
