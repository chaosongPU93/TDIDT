% lfecatmoment.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we have the catalog of the LFEs from every burst win of interest. 
% What remains is to add the other constraints we have - basically an estimate
% of individual event moments from Bostock, and estimate of the total slip
% in each SSE from geodesy.  Some things that would be useful:  (a) What 
% fraction of Bostock's detections, in terms of moment and total number, are
% within your burst time windows?  (b) What is the moment of your detections
% that are not in his catalog, compared to the moment (as determined by 
% Bostock) of those that are?  It would be relatively straightforward to look
% at your relative deconvolved amplitudes, but I found the best correlation
% to be between the Bostock moment and my amplitude-times-duration-squared
% (the expected correlation, since moment is two integrations of moment 
% acceleration), with duration determined by stretching the dipole.  Although
% even that correlation wasn't perfect.  
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/03/22
% Last modified date:   2023/03/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
format short e   % Set the format to 5-digit floating point
clear
clc
close all

%% for easy testing
defval('idxbst',1:195); %global indices of bursts to run 

set(0,'DefaultFigureVisible','on');
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
  % 'LZB  '
  % 'TWKB '
  % 'MGCB '
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
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.pgc002.',cutout(1:4)));
% trange = trange(1:end-1,:);
% trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.',num2str(ntol),'.pgc002.',...
%   cutout(1:4)));
tlen = trange(:,3)-trange(:,2);
nbst = size(trange,1);

dates = unique(trange(:,1));
years = unique(floor(dates/1000));
nets = length(years);

trangeout = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.',num2str(ntol),'.pgcout002.',...
  cutout(1:4)));


%%%load the empirical model param of time offset 14 for each addtional stations
off14mod = load(strcat(rstpath, '/MAPS/timeoff_planefitparam_4thsta_160sps'),'w+');

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
bostcati = bostcat(isinbnd == 1, :);  %inside boundary
% bostcato = bostcat(isinbnd ~= 1, :);  %outside boundary
% clear bostcat

%convert moment mag to moment, Mw = (2/3)*log_10(M0)-10.7, where M0 has the unit of dyne.cm
%(10^-7 N.m), so Mw = (2/3)*log_10(M0)-6 if M0 has the unit of N.m
bostcati(:,13) = mag2moment(bostcati(:,11));
% bostcato(:,13) = mag2moment(bostcato(:,11));

%if you choose only fam 002 regardless
% bostcati = bostcati(bostcati(:,1)==2,:);
% bostcato = bostcato(bostcato(:,1)~=2 & bostcato(:,1)~=47 & bostcato(:,1)~=246,:);


%% plot to see if 002 and 246 are still strongly overlapping
% aaa = unique(bostcati(:,1));
% bbb = [2003; 2004; 2005];
% cat = bostcati(bostcati(:,1)==2,:);
% cat2 = bostcati(bostcati(:,1)==47,:);
% cat3 = bostcati(bostcati(:,1)==246,:);
% 
% color=['r';'k'];
% figure
% for i = 1: length(bbb)
%   subplot(3,1,i);
%   hold on; box on; 
%   cat = bostcati(bostcati(:,2)==bbb(i),:);
%   for j = 1: length(aaa)
%     ind=find(cat(:,1)==aaa(j));
%     tmp=cat(ind,:);
%     time=[];
%     for k = 1:size(tmp,1)
%       time(k)=dat2jul(tmp(k,3),tmp(k,4),tmp(k,2))+tmp(k,5)/86400;
%     end
%     scatter(time,tmp(:,11),15,color(j,:),'filled');
%   end
%   xlabel(sprintf('Julian date on %d',bbb(i)));
%   ylabel('Magnitude');
% end


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
  bomsumday(i) = sum(temp(:,end));
  
  temp = bostcatof(bostcatof(:,2)==year & bostcatof(:,3)==a(1) & bostcatof(:,4)==a(2),:);
  bostdayof{i} = sortrows(temp, 5);
end
bostdayia = cat(1,bostdayi{:});   %MB 002/246 lfes of the same dates I case
bomsumdaya = sum(bomsumday);
factsamed = length(bostdayia)/length(bostcati);
% factsamedi = nbostdayi/length(bostcati);
bostdayofa = cat(1,bostdayof{:});   %MB other fams' lfes of the same dates I case
clear bostdayof

%so what happened to days I never focused, do I miss a lot of MB's LFEs?
bostotherd=setdiff(bostcati,bostdayia,'rows','stable');
factotherd = length(bostotherd)/length(bostcati);
aa = bostotherd(:,2)*10000+bostotherd(:,3)*100+bostotherd(:,4);
odate = unique(aa);
odatenum = zeros(length(odate),1);
for i = 1:length(odate)
  odatenum(i) = sum(aa==odate(i));
end

bb = bostdayia(:,2)*10000+bostdayia(:,3)*100+bostdayia(:,4);
idate = unique(bb);

%%
%%%daily # LFEs for dates analysed and days not 
figure
subplot(211)
box on; grid on; hold on;
scatter(1:length(idate),nbostdayi);
tklbl = [];
minnum = 40;
for i = 1: length(idate)
  tklbl{i} = num2str(idate(i));
end
ax=gca;
plot(ax.XLim,[minnum minnum],'k--');
text(1,minnum,'Min # 002 LFEs to select','VerticalAlignment','top');
text(0.98,0.95,sprintf('Total: %d/%d=%.1f%%',sum(nbostdayi),length(bostcati),factsamed*100),...
  'Units','normalized','HorizontalAlignment','right');
xticks(1:length(idate));
xticklabels(tklbl);
xtickangle(45);  
ylabel('# MB 002 LFEs');
title('Dates in 2003/2004/2005 analysed');

subplot(212)
box on; grid on; hold on;
scatter(1:length(odate),odatenum);
tklbl = [];
for i = 1:length(odate)
  tklbl{i} = num2str(odate(i));
end
ax=gca;
plot(ax.XLim,[minnum minnum],'k--');
text(1,minnum,'Min # 002 LFEs to select','VerticalAlignment','top');
text(1,odatenum(3),'Not chosen because no data at LZB','VerticalAlignment','top','color','r');
text(0.98,0.95,sprintf('Total: %d/%d',sum(odatenum),length(bostcati)),'Units','normalized',...
  'HorizontalAlignment','right');
xticks(1:length(odate));
xticklabels(tklbl);
xtickangle(45);  
ylabel('# MB 002 LFEs');
title('Dates in 2003/2004/2005 NOT analysed');

%%
%%%mag distribution of LFEs for dates analysed and days not 
figure
box on; grid on; hold on;
histogram(log10(bostdayia(:,end)),'Normalization','count','FaceColor','b','BinWidth',0.05);
histogram(log10(bostotherd(:,end)),'Normalization','count','FaceColor','r','BinWidth',0.05);
ax=gca;
plot([median(log10(bostdayia(:,end))) median(log10(bostdayia(:,end)))],ax.YLim,'b--','linew',2);
plot([median(log10(bostotherd(:,end))) median(log10(bostotherd(:,end)))],ax.YLim,'r--','linew',2);
legend('MB LFEs in my analysed days','Those in my ignored days');
xlabel('log_{10}{moment of MB LFEs}');
ylabel('Probability');

      
%%
% % bostname = ('/BOSTOCK/002_culled.mags.2003-2005');
% bostname = ('/BOSTOCK/002-246_culled.mags.2003-2005no2.6');
% allanver = ReformBostock(loc0(3),loc0(2),1,bostname);
% allanver(:,13) = 10.^(1.5*(6+allanver(:,11)));
% allanveri = cell(length(dates),1);
% for i = 1: length(dates)
%   date = dates(i);
%   year = floor(date/1000);
%   jday = floor(date-year*1000);
%   a = jul2dat(year,jday);
% 
%   temp = allanver(allanver(:,2)==year & allanver(:,3)==a(1) & allanver(:,4)==a(2),:);
%   allanveri{i} = sortrows(temp, 5);
%   allanverimsum(i) = sum(temp(:,end));
% end
% allanveria = cat(1,allanveri{:});
% allanverimsuma = sum(allanverimsum);

% keyboard

%%
%%%load the tremor catalog of John Armbruster, 
%%%2022/06/29, not really very useful as the detecting window could be 128-s long (not sure)
%format: [yyyy mm dd sec dx dy lon lat dep], 9 cols;
armcat = ReformArmbrusterv2(loc0(3),loc0(4),0);
[~,ind] = min(abs(armcat(:,5))+abs(armcat(:,6)));
% armcat1 = ReformArmbruster(loc0(3),loc0(4),1);
% [~,ind1] = min(abs(armcat1(:,5))+abs(armcat1(:,6)));
[iin,ion] = inpolygon(armcat(:,5),armcat(:,6),bnd(:,1),bnd(:,2));
isinbnd = iin | ion;
armcati = armcat(isinbnd == 1, :);  %inside boundary
armcato = armcat(isinbnd ~= 1, :);  %outside boundary
clear armcat

%% load LFE catalogs
%%%load the LFE catalog, 25-s-win
savefile = 'deconv_stats4th_allbstsig.mat';
allsig = load(strcat(rstpath, '/MAPS/',savefile));
imp = allsig.allbstsig.impindepall;
imp4th = allsig.allbstsig.impindep4thall;
nsrc = allsig.allbstsig.nsrcraw;
nsrc4th = allsig.allbstsig.nsrc;
rccbst = allsig.allbstsig.rccbst; %cat rcc
off1iwk = allsig.allbstsig.off1iwk; %alignment center of each subwin
off1ik = allsig.allbstsig.off1i;  %alignment of whole win
ccwpairk = allsig.allbstsig.ccwpairk; %0-lag max cc of each subwin
mrccwpairk = allsig.allbstsig.mrccwpairk; %median rcc if each subwin
ndetwin4thk = allsig.allbstsig.ndetwin4thk; %num of detections checked at 4th sta in each subwin
%convert time offset to relative loc
sps = 160;
[imploc, indinput] = off2space002(imp(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[imploc4th, indinput] = off2space002(imp4th(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%%%load the LFE catalog, whole-win
savefile = 'deconv1win_stats4th_allbstsig.mat';
allsig1win = load(strcat(rstpath, '/MAPS/',savefile));
imp1win = allsig1win.allbstsig.impindepall;
imp1win4th = allsig1win.allbstsig.impindep4thall;
nsrc1win = allsig1win.allbstsig.nsrcraw;
nsrc1win4th = allsig1win.allbstsig.nsrc;
%convert time offset to relative loc
[imploc1win, ~] = off2space002(imp1win(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[imploc1win4th, ~] = off2space002(imp1win4th(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%%%load the LFE catalog, 25-s-win, noise
savefile = 'deconv_stats4th_allbstnoi.mat';
allnoi = load(strcat(rstpath, '/MAPS/',savefile));
impn = allnoi.allbstnoi.impindepall;
impn4th = allnoi.allbstnoi.impindep4thall;
nsrcn = allnoi.allbstnoi.nsrcraw;
nsrcn4th = allnoi.allbstnoi.nsrc;
rccbstn = allnoi.allbstnoi.rccbst;
off1iwkn = allnoi.allbstnoi.off1iwk;
off1ikn = allnoi.allbstnoi.off1i;
%convert time offset to relative loc
[implocn, ~] = off2space002(impn(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[implocn4th, ~] = off2space002(impn4th(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%%%load the LFE catalog, whole-win, noise
savefile = 'deconv1win_stats4th_allbstnoi.mat';
allnoi1win = load(strcat(rstpath, '/MAPS/',savefile));
impn1win = allnoi1win.allbstnoi.impindepall;
impn1win4th = allnoi1win.allbstnoi.impindep4thall;
nsrcn1win = allnoi1win.allbstnoi.nsrcraw;
nsrcn1win4th = allnoi1win.allbstnoi.nsrc;
%convert time offset to relative loc
[implocn1win, ~] = off2space002(impn1win(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[implocn1win4th, ~] = off2space002(impn1win4th(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%%
impuse = imp;
nsrcuse = nsrc;
% impuse = imp4th;
% nsrcuse = nsrc4th;
rccuse = rccbst;
rcccol = 1;

%%
%ppeaks-zcrosses
ppkmzc = [10;12;10;7;10;7;8];

idxbst = 1:length(trange);
% idxbst = 181;

pltflag = 0;

%filtering passband for reading data, confirmed by 'spectrabursts002_4s.m'
hisig=6.3; % this will give a similar spectral shape between template and signal
losig=1.8;

rccmwsec = 0.5;

%moving window length in samples for running CC, envelope, etc.
%standard window length is about 0.5s, this is about the visual duration of the filtered and unfiltered
%template, although in fact to include as least one cycle of the main dipole of template
rccmwlen=rccmwsec*sps;
% rccmwlen=sps/2;
% rccmwlen=sps;

off1ic = zeros(size(trange,1),nsta);  % single best alignment 'computed' between ALL stas wrt 1 for entire win
off1i = zeros(size(trange,1),nsta);  % single best alignment 'actually used' between ALL stas wrt 1 for entire win
off14pred = zeros(size(trange,1),nsta-3); %empirical pred of off14 from plane fit given single best alignment
ccali = zeros(size(trange,1),1);  % CC value using the best alignment

nbosti = zeros(length(idxbst),1);
ratbosti = zeros(length(idxbst),1);
bost = cell(length(idxbst),1);  %all bostocks's lfes inside burst wins (and inside region)
bostmsum = zeros(length(idxbst),1);
bostday = [];

impi = cell(length(idxbst),1); %all my own lfes 
impasum = zeros(length(idxbst),1);
bostcomm = cell(length(idxbst),1);  %bostock's lfes that correspond to one of my lfes, inside burst wins (and inside region)
impcomm = cell(length(idxbst),1);   %my corresponding lfes
impuniq = cell(length(idxbst),1);   %my corresponding lfes
ampcomm = cell(length(idxbst),1);   %max and std of amp of my LFEs close to MB's in time
bostmiss = cell(length(idxbst),1);	%bostock's lfes that i missed, inside burst wins (and inside region)   
ncomm = zeros(length(idxbst),1); %num of commons
ratcomm = zeros(length(idxbst),1);  %ratio of commons

for iii = 1: length(idxbst)
  
  [iets,i,j] = indofburst(trange,idxbst(iii));
  
  year = years(iets);
  datesets = dates(floor(dates/1000)==year);
  
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
%   bostdayo = bostcato(bostcato(:,2)==year & bostcato(:,3)==a(1) & bostcato(:,4)==a(2),:);
%   bostdayo = sortrows(bostdayo, 5);
  
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
%   disp(k);
  
  tmaxi = hfdayi(:, seccol); % starting time of max power rate of half sec inside the ellipse
  tmaxo = hfdayo(:, seccol); % starting time of max power rate of half sec outside the ellipse
  tcnti = hfdayi(:, 15);  % the center of detecting win is the 15th col
  tcnto = hfdayo(:, 15);  % the center of detecting win is the 15th col
  tbosti = bostdayi(:,5); % (peak arrival) time of Bostock's LFE catalog inside the rectangle
%   tbosto = bostdayo(:,5); % (peak arrival) time of Bostock's LFE catalog outside the rectangle
  tarmi = armdayi(:,4); % (peak arrival?) time of Armbruster's tremor catalog inside the rectangle
  tarmo = armdayo(:,4); % (peak arrival?) time of Armbruster's tremor catalog outside the rectangle
  
  tst = rangetemp(j,2); % start and end time of bursts
  ted = rangetemp(j,3);
  
  %how many 4-s detections fall into the burst range
  indtmaxi = find(tmaxi>=tst-0.1 & tmaxi<=ted+0.1);
  ninbst(k,1) = length(indtmaxi); %+-0.1 s in case of resolution while saving to file
  
  %%%%Use the start and end of the 4-s detecting window
  tstbuf = min(tcnti(indtmaxi)-2);
  tedbuf = max(tcnti(indtmaxi)+2);
  tlenbuf = tedbuf-tstbuf;
      
  %plot the strongest 0.5-s arrival outside ellipse
  indto = find(tmaxo>=tstbuf & tmaxo<=tedbuf);
  %plot the strongest 0.5-s arrival inside ellipse
  indti = find(tmaxi>=tstbuf & tmaxi<=tedbuf);
%   %plot the bostock's LFE catalog outside rectangle
%   indbo = find(tbosto>=tstbuf & tbosto<=tedbuf);
  %plot the bostock's LFE catalog inside rectangle
  indbi = find(tbosti>=tstbuf & tbosti<=tedbuf);

  %my lfes
  ist = sum(nsrcuse(1:idxbst(iii)-1))+1;
  ied = ist+nsrcuse(idxbst(iii))-1;
  impii = impuse(ist:ied,:);  %LFEs of the same time win
  timp = impii(:,1)+ppkmzc(1);
  
  rcci = rccuse{k}; %rcc, 1st col is 25-s-win, 2nd col is whole-win
  
  if pltflag
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
    
    %max allowable shift in best alignment
    msftaddm = 1.5*sps+1;  %+1 for safety
    
    %have some overshoot, so that the resulted rcc would have the same length as the signal
    overshoot = rccmwlen/2;
    
    %%%%2022/09/26, obtain all information for data before decon, so that you know the threshold
    %%%%being used, although this was used mainly for noise experiment, we still use it here
    %chop a record segment
    optseg = STAopt(max(floor(tstbuf*sps+1-overshoot-msftaddm),1): ...
      min(floor(tedbuf*sps+overshoot+msftaddm),86400*sps), :); % sta 1
    ortseg = STAort(max(floor(tstbuf*sps+1-overshoot-msftaddm),1): ...
      min(floor(tedbuf*sps+overshoot+msftaddm),86400*sps), :);
    
    %%%obtain a single best alignment based on the entire win
    %       optcc = optseg(:, 2:end);
    optcc = detrend(optseg(1+msftaddm: end-msftaddm, 2:end));
    msftadd = (round(max(abs([off12ran off13ran])))+1)*sps/40;  %+1 for safety
    loffmax = 5*sps/40;
    ccmid = ceil(size(optcc,1)/2);
    ccwlen = round(size(optcc,1)-2*(msftadd+1));
    ccmin = 0.01;  % depending on the length of trace, cc could be very low
    iup = 1;    % times of upsampling
    [off12con,off13con,ccali(k),iloopoff,loopoff] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
      ccwlen,msftadd,loffmax,ccmin,iup);
    % if a better alignment cannot be achieved, use 0,0
    if off12con == msftadd+1 && off13con == msftadd+1
      off12con = 0;
      off13con = 0;
      fprintf('Tremor burst %d cannot be properly aligned, double-check needed \n',k);
    end
    off1ic(k,1) = 0;
    off1ic(k,2) = round(off12con);
    off1ic(k,3) = round(off13con);
    
    %%%obtain a single 'observational' alignment between 4th and 1st station, note that this might
    %%%be very different from the empirical prediction from 'empioffset4thsta002'
    %%%--there are also different options here, using sig, or sig-wlet CC as in decon routine
    mcoef = zeros(nsta-3, 3);
    mlag = zeros(nsta-3, 3);
    envrat = zeros(nsta-3, 3);
    for ista = 4: nsta
      [mcoef(ista-3, 1),off1ic(k,ista)] = xcorrmax(optcc(:,1), optcc(:,ista), msftadd, 'coeff');
      mlag(ista-3, 1) = off1ic(k,ista);
      envrat(ista-3, 1) = median(envelope(optcc(:,ista)))./median(envelope(optcc(:,1)));
      %do an overall CC between 4th and 2nd/3rd stas, to see which one they are most coehrent with
      for jjj = 2:3
        [mcoef(ista-3, jjj),mlag(ista-3, jjj)] = xcorrmax(optcc(:,1), optcc(:,ista), msftadd, 'coeff');
        envrat(ista-3, jjj) = median(envelope(optcc(:,ista)))./median(envelope(optcc(:,jjj)));
      end
      %empirical prediction from plane fitting in 'empioffset4thsta002'
      off14pred(k,ista-3) = round(off14mod(ista-3,1).*off1ic(k,2) + off14mod(ista-3,2).*off1ic(k,3) + ...
        off14mod(ista-3,3));
    end
    %       off1ia(k,4:end) = off14pred(k,:);
    
    %for real data, USE the best whole-win alignment before decon
    off1i(k,1:3) = off1ic(k,1:3);
    
    off1i(k,4:end) = off1ic(k,4:end); %if you trust overall alignment at 4th stations
    
    %%%Align and compute the RCC based on the entire win, and take that as the input signal!
    optdat = [];  % win segment of interest
    ortdat = [];
    optdat(:, 1) = optseg(1+msftaddm: end-msftaddm, 1); % time column
    ortdat(:, 1) = ortseg(1+msftaddm: end-msftaddm, 1);
    for ista = 1: nsta
      optdat(:, ista+1) = optseg(1+msftaddm-off1i(k,ista): end-msftaddm-off1i(k,ista), ista+1);
      ortdat(:, ista+1) = ortseg(1+msftaddm-off1i(k,ista): end-msftaddm-off1i(k,ista), ista+1);
    end
    
    %%%taper the signal and obtain the new rcc between tapered signals
    %%%2022/06/06, do NOT taper whatsoever!!
    sigsta = zeros(size(optdat,1), nsta);
    for ista = 1:nsta
      tmp = optdat(:,ista+1); %best aligned, filtered
      tmp = detrend(tmp);
      sigsta(:,ista) = tmp;
    end
    %compute running CC between 3 stations
    [ircc,rcc12,rcc13,rcc23] = RunningCC3sta(sigsta,rccmwlen);
    ircc = ircc-overshoot;
    rcc = (rcc12+rcc13)/2;
    rccpair = [rcc12 rcc13 rcc23];
    
    %compute the median absolute amplitude and envelope of the same moving window
    %for the moving window at the same station, sensable to use median
    [ir,ramp1,renv1] = Runningampenv(sigsta(:,1),rccmwlen,rccmwlen-1,'median');
    [~,ramp2,renv2] = Runningampenv(sigsta(:,2),rccmwlen,rccmwlen-1,'median');
    [~,ramp3,renv3] = Runningampenv(sigsta(:,3),rccmwlen,rccmwlen-1,'median');
    %looks like using the amplitude and envelope are pretty similar
    %maybe better to mean, otherwise stations PGC will be underrepresented due to a larger amp
    %variation
    ramp = mean([ramp1 ramp2 ramp3], 2);  % use mean or median??
    renv = mean([renv1 renv2 renv3], 2);
    %       ramp = median([ramp1 ramp2 ramp3], 2);  % use mean or median??
    %       renv = median([renv1 renv2 renv3], 2);
    sigsta = detrend(sigsta(overshoot+1:end-overshoot, :));  %excluding the overshoot
    
    %for ort. comp
    sigstaort = zeros(size(ortdat,1), nsta);
    for ista = 1:nsta
      tmp = ortdat(:,ista+1); %best aligned, filtered
      tmp = detrend(tmp);
      sigstaort(:,ista) = tmp;
    end
    [irccort,rcc12,rcc13,rcc23] = RunningCC3sta(sigstaort,rccmwlen);
    irccort = irccort-overshoot;
    rccort = (rcc12+rcc13)/2;
    
    %compute the median absolute amplitude and envelope of the same moving window
    %for the moving window at the same station, sensable to use median
    [ir,ramp1ort,renv1ort] = Runningampenv(sigstaort(:,1),rccmwlen,rccmwlen-1,'median');
    [~,ramp2ort,renv2ort] = Runningampenv(sigstaort(:,2),rccmwlen,rccmwlen-1,'median');
    [~,ramp3ort,renv3ort] = Runningampenv(sigstaort(:,3),rccmwlen,rccmwlen-1,'median');
    %looks like using the amplitude and envelope are pretty similar
    %maybe better to mean, otherwise stations PGC will be underrepresented due to a larger amp
    %variation
    ramport = mean([ramp1ort ramp2ort ramp3ort], 2);  % use mean or median??
    renvort = mean([renv1ort renv2ort renv3ort], 2);
%     ramport = median([ramp1ort ramp2ort ramp3ort], 2);  % use mean or median??
%     renvort = median([renv1ort renv2ort renv3ort], 2);
    
    sigstaort = detrend(sigstaort(overshoot+1:end-overshoot, :));  %excluding the overshoot
    
    renvrat = renv./renvort;
    
    time = optdat(overshoot+1:end-overshoot, 1);

    %%
    %plot the overview of both the signal and noise together
    widin = 12;  % maximum width allowed is 8.5 inches
    htin = 4;   % maximum height allowed is 11 inches
    nrow = 1;
    ncol = 1;
    f = initfig(widin,htin,nrow,ncol); %initialize fig
    
    xran = [0.08 0.92]; yran = [0.15 0.90];
    xsep = 0.02; ysep = 0.04;
    optaxpos(f,nrow,ncol,xran,yran,xsep,ysep);
    
    ym = max(abs(sigsta(:)));
    yran = [-1.1*ym 1.5*ym];
    
    ax = f.ax(1); hold(ax,'on');
    yyaxis(ax,'left');
    plot(ax,time-tstbuf,sigsta(:,1),'r-','linew',0.5);
    plot(ax,time-tstbuf,sigsta(:,2),'b-','linew',0.5);
    plot(ax,time-tstbuf,sigsta(:,3),'k-','linew',0.5); %,'linew',0.5
    xran = [time(1)-tstbuf time(end)-tstbuf]; % x axis range in sec
    xlim(ax,xran);
    ylim(ax,yran);
    longticks(ax,5);
    %plot the strongest 0.5-s arrival outside ellipse
    for ii = 1: length(indto)
      barst = tmaxo(indto(ii)); % arrival of strongest 0.5-s, in sec
      bared = tmaxo(indto(ii))+0.5; % last for 0.5 s
      yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.03;
      plot(ax,[barst bared]-tstbuf,[yloc yloc],'-','linew',2.5,'color',[.6 .6 .6]);
    end
    %plot the strongest 0.5-s arrival inside ellipse
    for ii = 1: length(indti)
      barst = tmaxi(indti(ii)); % arrival of strongest 0.5-s, in sec
      bared = tmaxi(indti(ii))+0.5; % last for 0.5 s
      plot(ax,[barst bared]-tstbuf,[yloc yloc],'k-','linew',2.5);
    end
    %       %plot the armbruster's tremor catalog outside rectangle
    %       ind = find(tarmo>=xran(1) & tarmo<=xran(2));
    %       yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.97;
    %       scatter(ax,tarmo(ind), yloc*ones(size(tarmo(ind))),10,'m','linew',1); % armbruster
    %       %plot the armbruster's tremor catalog inside rectangle
    %       ind = find(tarmi>=xran(1) & tarmi<=xran(2));
    %       scatter(ax,tarmi(ind), yloc*ones(size(tarmi(ind))),10,'m','filled');
    %plot the bostock's LFE catalog outside rectangle
    yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.06;
%     scatter(ax,tbosto(indbo)-tstbuf, yloc*ones(size(tbosto(indbo))),10,'g','linew',1); % bostock
%     plot(ax,[tbosto(indbo)-tstbuf tbosto(indbo)-tstbuf],ax.YLim,'k:');%[yran(1) 0]
    %plot the bostock's LFE catalog inside rectangle
    scatter(ax,tbosti(indbi)-tstbuf, yloc*ones(size(tbosti(indbi))),10,'g','filled');
    plot(ax,[tbosti(indbi)-tstbuf tbosti(indbi)-tstbuf],ax.YLim,'k--');
    %   plot(ax,[tstbuf tstbuf],ax.YLim,'--','Color',[.5 .5 .5],'linew',1);
    %   text(ax,tstbuf,yloc,'start','fontsize',8);
    %   plot(ax,[tedbuf tedbuf],ax.YLim,'--','Color',[.5 .5 .5],'linew',1);
    %   text(ax,tedbuf,yloc,'end','fontsize',8,'HorizontalAlignment','right');
    %plot my 3-sta LFE catalog's time at PGC, corrected from zero-crossing to postive peak
    yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.09;
    scatter(ax,timp/sps, yloc*ones(size(timp)),10,'r','filled'); % bostock
    %plot my 4-sta LFE catalog's time at PGC, corrected from zero-crossing to postive peak
    yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.12;
    scatter(ax,timp4th/sps, yloc*ones(size(timp4th)),10,'r','filled');
    text(ax,0.01,0.2,num2str(off1i(k,2)),'unit','normalized',...
      'HorizontalAlignment','left','fontsize',10,'color','b'); % offset 12
    text(ax,0.03,0.2,num2str(off1i(k,3)),'unit','normalized',...
      'HorizontalAlignment','left','fontsize',10,'color','k'); % offset 13
    text(ax,0.06,0.2,sprintf('%.3f',ccali(k)),'unit','normalized',...
      'HorizontalAlignment','left','fontsize',10,'color','k'); % mean CC of 3 pairs
%     text(ax,0.01,0.95,sprintf('W%d-T%d-B%d-S%d',k,length(indti),length(indbi),nsrcuse(idxbst(iii))),...
%       'unit','normalized','HorizontalAlignment','left','fontsize',10); % burst number
    title(ax,sprintf('W%d-T%d-B%d-S%d',k,length(indti),length(indbi),nsrcuse(idxbst(iii)))); % burst number
    xlabel(sprintf('Time (s) on %s %s %s',dy,mo,yr),'FontSize',12);
    ylabel(ax,'Amplitude','FontSize',12);
    yyaxis(ax,'right');
%     plot(ax,(1:size(renvrat,1))/sps,log10(renvrat),'o','markersize',0.5,'Color','c'); % scale with data amp.    
    plot(ax,(1:size(renvrat,1))/sps,log10(renvrat),'-','Color','c','linew',0.5); % scale with data amp.    
    axsym(ax,2,10);
%     plot(ax,(1:size(rcci,1))/sps,yran(2)/2*rcci(:,rcccol)+yran(2)/2,'o','markersize',0.5,'Color',[.5 .5 .5]); % scale with data amp.
    plot(ax,(1:size(rcci,1))/sps,ax.YLim(2)/4*rcci(:,rcccol)+ax.YLim(2)/4*3,'-','Color',[.6 .6 .6],'linew',0.5); % scale with data amp.
    %       plot(ax,ir/sps+tstbuf,renv*1.2+0.3*yran(2),'Color',[.6 .6 .6],'linew',0.5);  % scale it to the amplitude range
%     ylim(ax,[-1 1]*yran(2));
    cc = xcorr(rcci(:,rcccol), log10(renvrat),0,'normalized');
    text(ax,0.06,0.8,sprintf('%.3f',cc),'unit','normalized',...
      'HorizontalAlignment','left','fontsize',10,'color','k'); % offset 13
    ylabel(ax,'log_{10}{Running env rat opt/ort}','FontSize',12);
    hold(ax,'off');
    
  end
  
  %which bostock's LFEs are within my time burst window? Fraction? Moment?
  bost{iii} = bostdayi(indbi,:); %all bostocks's lfes inside burst wins (and inside region)
  impi{iii} = impuse(ist:ied,:); %all my own lfes 

  temp1 = [];
  temp2 = [];
  temp3 = [];
  temp4 = [];
  for j = 1: length(indbi)
    % [toff,ind] = min(abs(timp/sps-(tbosti(indbi(j))-tstbuf))); %which my LFE is closest to Bostock's?
    % if toff<= 0.25/2  %if the difference is small enough, then there is a correspondance
    %   temp1 = [temp1; bostdayi(indbi(j),:)];
    %   temp2 = [temp2; impii(ind,:)];
    % else
    %   temp3 = [temp3; bostdayi(indbi(j),:)];
    % end

    maxtoff = 0.5;
    ind = find(abs(timp/sps-(tbosti(indbi(j))-tstbuf)) <= maxtoff); %find all close enough
    if ~isempty(ind)
      temp1 = [temp1; bostdayi(indbi(j),:)];
      if length(ind) == 1
        temp2 = [temp2; ind];  %if just one, use it
        temp4 = [temp4; mean(impii(ind,[2 4 6]),2) 1 0];
      else
        amptemp = mean(impii(ind,[2 4 6]),2);
        [ampmax,ind2] = max(amptemp);
        temp2 = [temp2; ind(ind2)];  %if multiple, choose the one with max amp
        temp4 = [temp4; ampmax length(ind) std(amptemp)];
      end
    else
      temp3 = [temp3; bostdayi(indbi(j),:)];  %otherwise it is a 'miss'
    end

  end
  bostcomm{iii} = temp1;  %bostock's lfes that correspond to one of my lfes, inside burst wins (and inside region)
  impcomm{iii} = impii(temp2,:);   %my corresponding lfes
  if isempty(temp2)
    impuniq{iii} = impii;
  else
    if ~isequaln(unique(temp2),temp2)
      fprintf('%d DUPLICATES! ',k);
    end
    tmp = setdiff(1:size(impii,1),temp2);
    impuniq{iii} = impii(tmp,:);
  end
  ampcomm{iii} = temp4;
  bostmiss{iii} = temp3;	%bostock's lfes that i missed, inside burst wins (and inside region) 
  ncomm(iii) = size(temp1,1); %num of commons
  ratcomm(iii) = size(temp1,1)/length(indbi);  %ratio of commons  
  
%   if ~isempty(temp1)
%     %diff time between consecutive common MB's LFEs
%     temp1 = sortrows(temp1,5,'ascend');
%     bostcommdt{iii} = diffcustom(temp1(:,5),1,'forward');
% 
%     %diff time and amp between each common MB's LFEs and all others within say, 20 s
%     targetind = target2all(temp1(:,5),[0 20]);
%     dtime = [];
%     pairmom = [];
%     for kk = 1: size(targetind,1)
%       ind = targetind{kk};
%       aa = abs(temp1(ind,5)-temp1(kk,5));  %diff time
%       dtime = [dtime; aa];
%       bb = [temp1(kk,end)*ones(length(ind),1) temp1(ind,end)]; %moment of the pair
%       pairmom = [pairmom; bb];
%     end
%     bostcommdtall{iii} = dtime; 
%     bostcommpairmom{iii} = pairmom; 
%   else
%     bostcommdt{iii} = [];
%     bostcommdtall{iii} = []; 
%     bostcommpairmom{iii} = []; 
%   end
  
end
% keyboard

%lump all bursts
bostcomma = cat(1,bostcomm{:});
impcomma = cat(1,impcomm{:});
impuniqa = cat(1,impuniq{:});
ampcomma = cat(1,ampcomm{:});
bostmissa = cat(1,bostmiss{:});
bosta = cat(1,bost{:});
% bostcommdta = cat(1,bostcommdt{:});
% bostcommdtalla = cat(1,bostcommdtall{:});
% bostcommpairmoma = cat(1,bostcommpairmom{:});

%%
f = initfig(10,5,1,2); %initialize fig
ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
ind = find(ampcomma(:,2)>1);
histogram(ax,ampcomma(ind,2));
xlabel(ax,sprintf("# of my lfes w/in %.2fs of MB's",2*maxtoff));
ylabel(ax,'Count');

ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
histogram(ax,ampcomma(ind,3));
xlabel(ax,'std of their amp');
ylabel(ax,'Count');


%% distribution in diff time between common MB's sources
% f = initfig(10,8,2,2); %initialize fig
% ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% histogram(ax,bostcommdta,'binw',0.25,'facec','k');
% xlabel(ax,'diff time between consecutive sources');
% ylabel(ax,'Count');
% ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% histogram(ax,bostcommdtalla,'binw',0.25,'facec','k');
% xlabel(ax,'diff time between sources and all others');
% ylabel(ax,'Count');
% ax=f.ax(3); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% scatter(ax,bostcommdtalla,log10(bostcommpairmoma(:,2)),20,'ko');
% xlabel(ax,'diff time between sources and all others');
% ylabel(ax,'Moment of the pair (latter)');
% ax=f.ax(4); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% p(1)=histogram(ax,log10(bostcommpairmoma(:,1)),'facec','r');
% p(2)=histogram(ax,log10(bostcommpairmoma(:,2)),'facec','k');
% legend(ax,p,'former','latter');
% xlabel(ax,'Moment of the pair');
% ylabel(ax,'Count');
% 
% keyboard

%% plot of correlation between MB LFE moments and my LFE amp
widin = 8;  % maximum width allowed is 8.5 inches
htin = 7;   % maximum height allowed is 11 inches
nrow = 2;
ncol = 2;
f = initfig(widin,htin,nrow,ncol); %initialize fig
  
fttpfree = fittype( @(a,b,x) a*x+b);

ax=f.ax(1);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
lfeamp = mean(impcomma(:,[2 4 6]),2);
% lfeamp = median(impcomma(:,[2 4 6]),2);
y = log10(bostcomma(:,end));
x = log10(lfeamp);
scatter(ax,x,y,15,'filled');
%linear robust least square
[fitobj,gof,output] = fit(x,y,fttpfree,'Robust','Bisquare','StartPoint',[1 1]);
% [fitobj,gof,output] = fit(x,y,'poly1');
stats = statsofrobustlnfit(fitobj,gof,output,x,y);
fitx = linspace(-1.2,1.2,100);
fity = feval(fitobj,fitx);
plot(ax,fitx,fity,'--','linewidth',2,'color','r');
text(ax,0.95,0.95,sprintf('y = %.2f x + %.2f',stats.slope,stats.intcpt),'HorizontalAlignment',...
  'right','Units','normalized');
text(ax,0.95,0.90,sprintf('weighted pearson = %.2f',stats.pearwt),'HorizontalAlignment',...
  'right','Units','normalized');
% text(ax,0.95,0.25,sprintf('In burst: %d; %.2f; %.2e',nbostia,rbostia,bostiams),'HorizontalAlignment',...
%   'right','Units','normalized');
% text(ax,0.95,0.2,sprintf('Outside: %d; %.2f; %.2e',nbostimis,rbostimis,bostimisms),'HorizontalAlignment',...
%   'right','Units','normalized');
text(ax,0.95,0.05,sprintf('%d/%d in common',length(x),size(bosta,1)),'HorizontalAlignment',...
  'right','Units','normalized');
ylabel(ax,"Y-->log_{10}{Michael's moments (NM)}");
xlabel(ax,"X-->log_{10}{My amps of common LFEs}");

ax=f.ax(2);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
histogram(ax,ratcomm,'BinWidth',0.1);
xlabel(ax,'Ratio of common LFEs');
ylabel(ax,'Count');

ax=f.ax(3);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
histogram(ax,bostcomma(:,end)./lfeamp);
xlabel(ax,'Ratio');
ylabel(ax,'Count');

ax=f.ax(4);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
histogram(ax,y-x);
xlabel(ax,'Y-X [log_{10}{ratio}]');
ylabel(ax,'Count');

% supertit(f.ax(1:2),'data, 25-s wins, 2ndary sources removed');
% supertit(f.ax(1:2),'data, 25-s wins, 4th station checked');
% supertit(f.ax(1:2),'data, whole-win, 2ndary sources removed');
% supertit(f.ax(1:2),'data, 25-s wins, 4th station checked');

% keyboard

%% infer my total moment according to relation above
%%%Estimated my moment for my whole catalog of LFEs
lfeamp = mean(impuse(:,[2 4 6]),2); %amp for all LFE catalog
% lfeamp = median(impuse(:,[2 4 6]),2);
fitx = log10(lfeamp); 
fity = feval(fitobj,fitx);  %use correlation above to infer the log10 of moment
lfemoall = 10.^(fity);
mosumallb = sum(lfemoall)  %total moment for my catalog

lfeamp = mean(impcomma(:,[2 4 6]),2); %amp for common LFE catalog
fitx = log10(lfeamp); 
fity = feval(fitobj,fitx);  %use correlation above to infer the log10 of moment
lfemocomm = 10.^(fity);
mosumcomm = sum(lfemocomm)  %total moment for common lfes

%%%Bostock's moment sum 
bmosummiss = sum(bostmissa(:,end)) %moment sum of missed
bmosumcomm = sum(bostcomma(:,end)) %moment sum of commons
bmosum = sum(bosta(:,end)) %moment sum of all

%% moment correction
%the above moment could be deviated from truth in several ways, see NOTE for details:
%1. My burst wins missed half of MB's LFEs on same dates same years, if my catalog can be projected
%to those missing times, then the moment will increase
%2. My catalog is smaller after checking at KLNB, moment will decrease
%3. MB's LFEs' moments could be overestimate as he looks at 1-2 Hz, moment will decrease
%4. if amp is proportional to square root of moments, then moment will increase

%%%first of all, is there any systematic bias in amp/mag of different set of LFEs?
%%%distribution of mag of common and missed LFEs
widin = 15;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 3;
f = initfig(widin,htin,nrow,ncol); %initialize fig

ax=f.ax(1);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
p1=histogram(ax,log10(bostdayia(:,end)),'FaceColor','k','Normalization','count','BinWidth',0.05);
p2=histogram(ax,log10(bosta(:,end)),'FaceColor','b','Normalization','count','BinWidth',0.05);
fracinbst = length(bosta)/length(bostdayia)
%for days we share, why did I miss a lot of MB's LFEs? 
bostoutbst = setdiff(bostdayia,bosta,'rows','stable');  %same dates, but outside bursts
p3=histogram(ax,log10(bostoutbst(:,end)),'FaceColor','r','Normalization','count','BinWidth',0.05);
xlabel(ax,'log_{10}{moment of MB LFEs}');
% ylabel(ax,'Probability');
ylabel(ax,'Count');
legend(ax,[p1,p2,p3],'MB LFEs in my analysed days','Those in my burst wins',...
  'Those not in my burst wins');

ax=f.ax(2);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
p1=histogram(ax,log10(bosta(:,end)),'FaceColor','k','Normalization','count','BinWidth',0.05);
p2=histogram(ax,log10(bostcomma(:,end)),'FaceColor','b','Normalization','count','BinWidth',0.05);
fraccomm = length(bostcomma)/length(bosta)
p3=histogram(ax,log10(bostmissa(:,end)),'FaceColor','r','Normalization','count','BinWidth',0.05);
xlabel(ax,'log_{10}{moment of MB LFEs}');
% ylabel(ax,'Probability');
ylabel(ax,'Count');
legend(ax,[p1,p2,p3],'MB LFEs in my burst wins','Those I share','Those I miss');

ax=f.ax(3);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
p1=histogram(ax,log10(mean(impuse(:,[2 4 6]),2)),'FaceColor','k','Normalization','count',...
  'BinWidth',0.1);
p3=histogram(ax,log10(mean(impuniqa(:,[2 4 6]),2)),'FaceColor','r','Normalization','count','BinWidth',0.1);
p2=histogram(ax,log10(mean(impcomma(:,[2 4 6]),2)),'FaceColor','b','Normalization','count',...
  'BinWidth',0.1);
xlabel(ax,'log_{10}{Amp of My LFEs}');
% ylabel(ax,'Probability');
ylabel(ax,'Count');
legend(ax,[p1,p2,p3],'My LFEs in my burst wins','Those MB shares','Those unique');

%%%correction #2, after checking at KLNB

%%%correction #4
f = initfig(5,5,1,1);
ax=f.ax(1);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
p1=histogram(ax,log10(lfemoall),'FaceColor','b','Normalization','count',...
  'BinWidth',0.1);
text(ax,0.95,0.7,sprintf('before: %.1e Nm',sum(lfemoall)),'HorizontalAlignment',...
  'right','Units','normalized');
ratio = lfemoall./prctile(lfemoall,5);
ratio(ratio<1) = 1;
lfemoall = lfemoall.*ratio;

p2=histogram(ax,log10(lfemoall),'FaceColor','r','Normalization','count',...
  'BinWidth',0.1);
text(ax,0.95,0.8,sprintf('after: %.1e Nm',sum(lfemoall)),'HorizontalAlignment',...
  'right','Units','normalized');
xlabel(ax,'log_{10}{moment (Nm)}');
ylabel(ax,'Count');
legend(ax,[p1,p2],'before correction','after correction');

%%%correction #1, project to those missing times, given the bias is negligible
lfemoall = lfemoall./fracinbst; %direct projection as the bias is low

mosumalla = sum(lfemoall)  %total moment for my catalog
mocor = mosumalla/mosumallb

% keyboard

%% plot implied cumulative moment and slip
%%%moment per pixel in sample space
[denuniq1d, inddup] = density_pixel(impuse(:,7),impuse(:,8)); %count at each unique loc
mosumuniq = sum_at_indices(lfemoall,inddup);
mosumuniq1d = [denuniq1d(:,1:2) mosumuniq];  %sum of all at each unique loc
cstr = {'# detections / pixel'; 'moment (Nm) / pixel'};
[f] = plt_sum_pixel(denuniq1d,mosumuniq1d,[-40 40],[-40 40],20,cstr);
%%%moment per grid in map view
dx = 0.4; dy = 0.4;
loc = off2space002(impuse(:,7:8),sps,ftrans,0);
[den1d,indices] = density_matrix(loc(:,1),loc(:,2),[-4 4],[-4 4],dx,dy);
den1d = den1d(den1d(:,3)>0, :);
indices = indices(~cellfun('isempty',indices));
mosum = sum_at_indices(lfemoall,indices);
mosum1d = [den1d(:,1:2) mosum];
cstr = {'# detections / grid'; 'moment (Nm) / grid'};
[f] = plt_sum_pixel(den1d,mosum1d,[-4 4],[-4 4],300,cstr,'s');
hold(f.ax(1),'on');
plot(f.ax(1),xcut,ycut,'k-','linew',2);
hold(f.ax(1),'off');

%%%moment per detection per pixel in sample space
moaveuniq1d = [mosumuniq1d(:,1:2) mosumuniq1d(:,3)./denuniq1d(:,3)];
cstr = {'# detections / pixel'; 'moment per detection / pixel'};
[f] = plt_sum_pixel(denuniq1d,moaveuniq1d,[-40 40],[-40 40],20,cstr);

%%%moment per detection per grid in map view
moave1d = [den1d(:,1:2) mosum./den1d(:,3)];
cstr = {'# detections / grid'; 'moment per detection / grid'};
[f] = plt_sum_pixel(den1d,moave1d,[-4 4],[-4 4],300,cstr,'s');
hold(f.ax(1),'on');
plot(f.ax(1),xcut,ycut,'k-','linew',2);
hold(f.ax(1),'off');
errgrid1 = [0 0;
           1 0;
           -1 0;
           0 1;
           0 -1
           ];
errloc1 = off2space002(errgrid1,sps,ftrans,0);
sft1 = [1.5 -3];
errgrid2 = 4*errgrid1;
errloc2 = off2space002(errgrid2,sps,ftrans,0);
sft2 = [2.5 -3];
hold(f.ax(2),'on');
plot(f.ax(2),errloc1([1 2],1)+sft1(1), errloc1([1 2],2)+sft1(2),'r-','linew',1);
plot(f.ax(2),errloc1([1 3],1)+sft1(1), errloc1([1 3],2)+sft1(2),'r-','linew',1);
plot(f.ax(2),errloc1([1 4],1)+sft1(1), errloc1([1 4],2)+sft1(2),'r-','linew',1);
plot(f.ax(2),errloc1([1 5],1)+sft1(1), errloc1([1 5],2)+sft1(2),'r-','linew',1);
plot(f.ax(2),errloc2([1 2],1)+sft2(1), errloc2([1 2],2)+sft2(2),'b-','linew',1);
plot(f.ax(2),errloc2([1 3],1)+sft2(1), errloc2([1 3],2)+sft2(2),'b-','linew',1);
plot(f.ax(2),errloc2([1 4],1)+sft2(1), errloc2([1 4],2)+sft2(2),'b-','linew',1);
plot(f.ax(2),errloc2([1 5],1)+sft2(1), errloc2([1 5],2)+sft2(2),'b-','linew',1);
hold(f.ax(2),'off');

%%
%%%convert moment to average slip per grid in map view
locuniq = off2space002(mosumuniq1d(:,1:2),sps,ftrans,0);
[dengrid1d,indgrid] = density_matrix(locuniq(:,1),locuniq(:,2),[-4 4],[-4 4],dx,dy);
dengrid1d = dengrid1d(dengrid1d(:,3)>0, :);
indgrid = indgrid(~cellfun('isempty',indgrid));
mosumgrid = sum_at_indices(mosumuniq,indgrid);  % sum of moment in each grid cell
mosumgrid1d = [dengrid1d(:,1:2) mosumgrid];
mu = 3e10;
slipavegrid = mosumgrid / (mu*dx*dy)/ 1e3; %has the unit of mm .* dengrid1d(:,3)/ sum(dengrid1d(:,3))
slipavegrid1d = [dengrid1d(:,1:2) slipavegrid];

%burst #181
iii = 181;
ist = sum(nsrcuse(1:idxbst(iii)-1))+1;
ied = ist+nsrcuse(idxbst(iii))-1;
impi = impuse(ist:ied,:);  %LFEs of the same time win
[denuniq1di, inddupi] = density_pixel(impi(:,7),impi(:,8)); %count at each unique loc
mosumuniqi = sum_at_indices(lfemoall(ist:ied),inddupi);
mosumuniq1di = [denuniq1di(:,1:2) mosumuniqi];  %sum of all at each unique loc
moaveuniq1di = [mosumuniq1di(:,1:2) mosumuniq1di(:,3)./denuniq1di(:,3)];

locuniqi = off2space002(mosumuniq1di(:,1:2),sps,ftrans,0);
[dengrid1di,indgridi] = density_matrix(locuniqi(:,1),locuniqi(:,2),[-4 4],[-4 4],dx,dy);
dengrid1di = dengrid1di(dengrid1di(:,3)>0, :);
indgridi = indgridi(~cellfun('isempty',indgridi));
mosumgridi = sum_at_indices(mosumuniqi,indgridi);  % sum of moment in each grid cell
slipavegridi = mosumgridi / (mu*dx*dy)/ 1e3; %has the unit of mm .* dengrid1d(:,3)/ sum(dengrid1d(:,3))
slipavegrid1di = [dengrid1di(:,1:2) slipavegridi];

loci = loc(ist:ied,:);
[den1di,indicesi] = density_matrix(loci(:,1),loci(:,2),[-4 4],[-4 4],dx,dy);
den1di = den1di(den1di(:,3)>0, :);
indicesi = indicesi(~cellfun('isempty',indicesi));


%%%density and slip for #181
f = initfig(9,4.5,1,2); %initialize fig
pltxran = [0.05 0.95]; pltyran = [0.05 0.95];
pltxsep = 0.06; pltysep = 0.06; 
optaxpos(f,1,2,pltxran,pltyran,pltxsep,pltysep);

%cumulative density
ax=f.ax(1);
hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
dum = den1di;
dum(dum(:,3)>1, :) = [];
dum(:,3) = log10(dum(:,3));
scatter(ax,dum(:,1),dum(:,2),220,dum(:,3),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dum = sortrows(den1di,3);
dum(dum(:,3)==1, :) = [];
dum(:,3) = log10(dum(:,3));
scatter(ax,dum(:,1),dum(:,2),220,dum(:,3),'s','filled','MarkerEdgeColor','none');
% oldc = colormap(ax,'kelicol');
% newc = flipud(oldc);
% colormap(ax,newc);
colormap(ax,flipud(colormap(ax,'kelicol')));
c=colorbar(ax,'SouthOutside');
ax.CLim(2) = prctile(dum(:,3),99);
c.Label.String = strcat('log_{10}(# detections / grid)');
axis(ax,'equal');
axis(ax,[-4 4 -4 4]);
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
%the cut-out boundary of 4-s detections
plot(ax,xcut,ycut,'k-','linew',2);
scatter(ax,loci(:,1),loci(:,2),15,[.5 .5 .5],'o','filled',...
  'MarkerEdgeColor','none');
longticks(ax,1.5);

%average slip
ax=f.ax(2);
hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
slipavegrid1di = sortrows(slipavegrid1di,3);
scatter(ax,slipavegrid1di(:,1),slipavegrid1di(:,2),220,slipavegrid1di(:,3),'s','filled',...
  'MarkerEdgeColor','none');
colormap(ax,flipud(colormap(ax,'kelicol')));
c=colorbar(ax,'SouthOutside'); %
caxis(ax,[0 prctile(slipavegrid1di(:,3),99)]);
% c.Label.String = 'log_{10}(average slip (mm)/ grid)';
c.Label.String = 'average slip (mm)/ grid';
axis(ax,'equal');
axis(ax,[-4 4 -4 4]);
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
%the cut-out boundary of 4-s detections
plot(ax,xcut,ycut,'k-','linew',2);
longticks(ax,1.5);

%%%what's average slip over the entire area?
indbnd = convhull(slipavegrid1d(:,1),slipavegrid1d(:,2));
area = 0.5*polyarea(slipavegrid1d(indbnd,1),slipavegrid1d(indbnd,2));
slipave = sum(lfemoall) ./ mu ./ area ./ 1e3; %has the unit of mm 

%%%what's average slip over each concentric ellipses?
semia = 1.75*(0.6:0.2:2.0);
semib = 1.25*(0.6:0.2:2.0);
nreg = length(semia);
for ireg = 1: nreg
  xaxis = semia(ireg);
  yaxis = semib(ireg);
  % xaxis=1.75; %axis length of the same ellipse of my 4-s catalog
  % yaxis=1.25;
  shiftor=[0.2 0.2]; %(in km) %center of the same ellipse of my 4-s catalog
  [xcut1,ycut1] = ellipse_chao(shiftor(1),shiftor(2),xaxis,yaxis,0.01,45,shiftor);
  bnd = [xcut1 ycut1];
  [iin,ion] = inpolygon(loc(:,1),loc(:,2),bnd(:,1),bnd(:,2));
  isinbnd = iin | ion;
  mosumell = sum(lfemoall(isinbnd == 1));  %total moment for my catalog
  area = pi*xaxis*yaxis;
  slipaveell(ireg) = mosumell ./ mu ./ area ./ 1e3; %has the unit of mm 
end

%%
%%%average slip for the whole catalog
f = initfig(9,9,2,2); %initialize fig
pltxran = [0.05 0.95]; pltyran = [0.05 0.95];
pltxsep = 0.06; pltysep = 0.06; 
optaxpos(f,2,2,pltxran,pltyran,pltxsep,pltysep);

%slip
ax=f.ax(1);
hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
slipavegrid1d = sortrows(slipavegrid1d,3);
% scatter(ax,slipavegrid1d(:,1),slipavegrid1d(:,2),300,log10(slipavegrid1d(:,3)),'s','filled',...
%   'MarkerEdgeColor','none');
scatter(ax,slipavegrid1d(:,1),slipavegrid1d(:,2),220,slipavegrid1d(:,3),'s','filled',...
  'MarkerEdgeColor','none');
colormap(ax,flipud(colormap(ax,'kelicol')));
c=colorbar(ax,'SouthOutside'); %
caxis(ax,[0 prctile(slipavegrid1d(:,3),99)]);
% c.Label.String = 'log_{10}(average slip (mm)/ grid)';
c.Label.String = 'average slip (mm)/ grid';
axis(ax,'equal');
axis(ax,[-4 4 -4 4]);
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
%the cut-out boundary of 4-s detections
plot(ax,xcut,ycut,'k-','linew',2);
loci = loc(ist:ied,:);
scatter(ax,loci(:,1),loci(:,2),15,[.5 .5 .5],'o','filled',...
  'MarkerEdgeColor','none');
longticks(ax,1.5);

%slip contour
ax=f.ax(2);
hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
plot(ax,xcut,ycut,'k-','linew',2);
off1iwa = [];
plta = [];
for iii = 1: length(idxbst)
  off1iw = off1iwk{idxbst(iii)};
  off1iwa = [off1iwa; off1iw(:,2:3)];

  ccwpair = ccwpairk{idxbst(iii)};
  ccw = mean(ccwpair(:,1:2),2); %0-lag max cc of each subwin
  mrccwpair =  mrccwpairk{idxbst(iii)};
  mrccw = mean(mrccwpair(:,1:2),2); %median rcc of each subwin
  ndetw4th = ndetwin4thk{idxbst(iii)}; %num of detections checked at 4th sta in each subwin
  plta = [plta; ccw];
end
[denoff1iwa, ind] = density_pixel(off1iwa(:,1),off1iwa(:,2)); %count at each unique loc
medplta = median_at_indices(plta,ind);
locoff1iw = off2space002(denoff1iwa(:,1:2),sps,ftrans,0);
% keyboard
scatter(ax,locoff1iw(:,1),locoff1iw(:,2),10,medplta,'^','filled',...
  'MarkerEdgeColor','none');
%   scatter(ax,locoff1iw(:,1),locoff1iw(:,2),10,'c','^','filled',...
%     'MarkerEdgeColor','none');
colormap(ax,flipud(colormap(ax,'kelicol')));
c=colorbar(ax,'SouthOutside'); %
caxis(ax,[prctile(plta,1) prctile(plta,99)]);
c.Label.String = 'max CC of 25-s win';  %max CC, 'median RCC', 'num of det checked at KLNB'
% [xyzgridpad,xgrid,ygrid,zgrid,ind2] = zeropadmat2d(slipavegrid1d,-4+dx/2: dx: 4-dx/2,...
%   -4+dy/2: dy: 4-dy/2);
% zgridgf = imgaussfilt(zgrid, 1);  %smooth it a bit
% perc = 10:20:90;
% conplt = round(prctile(slipavegrid1d(:,3),perc));
% % conplt = 5: 20: 160;
% conmat = contour(ax,xgrid,ygrid,zgridgf,conplt,'-','color','k','ShowText','on'); %
% contour(ax,xgrid,ygrid,zgridgf,round([median(slipavegrid1d(:,3)),median(slipavegrid1d(:,3))]),...
%   '-','color','r','linew',1.5);
for ireg = 1: 2: nreg
  xaxis = semia(ireg);
  yaxis = semib(ireg);
  shiftor=[0.2 0.2]; %(in km) %center of the same ellipse of my 4-s catalog
  [xcut1,ycut1] = ellipse_chao(shiftor(1),shiftor(2),xaxis,yaxis,0.01,45,shiftor);
  plot(ax,xcut1,ycut1,'-','linew',1,'color',[.5 .5 .5]);
  text(ax,0.2,0.98*max(ycut1),sprintf('%d',round(slipaveell(ireg))),'fontsize',9,...
    'VerticalAlignment','top');
end
axis(ax,'equal');
axis(ax,[-4 4 -4 4]);
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
longticks(ax,1.5);
% keyboard

%slip + 4-s tremor loc
ax=f.ax(3);
hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
scatter(ax,slipavegrid1d(:,1),slipavegrid1d(:,2),100,slipavegrid1d(:,3),'s','filled',...
  'MarkerEdgeColor','none');
colormap(ax,flipud(colormap(ax,'kelicol')));
c=colorbar(ax,'SouthOutside'); %
caxis(ax,[0 prctile(slipavegrid1d(:,3),99)]);
% c.Label.String = 'log_{10}(average slip (mm)/ grid)';
c.Label.String = 'average slip (mm)/ grid';
axis(ax,'equal');
axis(ax,[-7 5 -6 6]);
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
scatter(ax,hfall(:,1),hfall(:,2),3,[.3 .3 .3],'o','filled',...
  'MarkerEdgeColor','none');
longticks(ax,1.5);

%slip contour + 4-s tremor density
ax=f.ax(4);
hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
dentmr = density_matrix(hfall(:,1),hfall(:,2),[-7 5],[-6 6],dx,dy);
dentmr = dentmr(dentmr(:,3)>0, :);
scatter(ax,dentmr(:,1),dentmr(:,2),100,dentmr(:,3),'s','filled',...
  'MarkerEdgeColor','none');
colormap(ax,flipud(colormap(ax,'kelicol')));
c=colorbar(ax,'SouthOutside'); %
caxis(ax,[0 prctile(dentmr(:,3),99)]);
c.Label.String = '# tremor / grid';
% conmat = contour(ax,xgrid,ygrid,zgridgf,conplt,'-','color',[.5 .5 .5],'ShowText','on','linew',1); %
% contour(ax,xgrid,ygrid,zgridgf,round([median(slipavegrid1d(:,3)),median(slipavegrid1d(:,3))]),...
%   '-','color','r','linew',1.5); %
plot(ax,xcut,ycut,'k-','linew',2);
axis(ax,'equal');
axis(ax,[-7 5 -6 6]);
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
longticks(ax,1.5);


keyboard


%%
% %%%get the moment sum for each single burst
% lfemoi = cell(length(idxbst),1);
% lfemocommi = cell(length(idxbst),1);
% mosumi = zeros(length(idxbst),1);
% mosumcommi = zeros(length(idxbst),1);
% bmosumi = zeros(length(idxbst),1);
% bmosummissi = zeros(length(idxbst),1);
% bmosumcommi = zeros(length(idxbst),1);
% for iii = 1: length(idxbst)
%   temp = impi{iii};
%   if ~isempty(temp)
%     lfeamp = mean(temp(:,[2 4 6]),2); 
%     fitx = log10(lfeamp); 
%     fity = feval(fitobj,fitx);  %use correlation above to infer the log10 of moment
%     lfemoi{iii} = 10.^(fity);
%     mosumi(iii) = sum(lfemoi{iii});  %total moment for my catalog
%   end
%   
%   temp = impcomm{iii};
%   if ~isempty(temp)
%     lfeamp = mean(temp(:,[2 4 6]),2); 
%     fitx = log10(lfeamp); 
%     fity = feval(fitobj,fitx);  %use correlation above to infer the log10 of moment
%     lfemocommi{iii} = 10.^(fity);
%     mosumcommi(iii) = sum(lfemocommi{iii}); 
%   end
% 
%   temp = bost{iii};
%   if ~isempty(temp)
%     bmosumi(iii) = sum(temp(:,end));
%   end
%   temp = bostmiss{iii};
%   if ~isempty(temp)
%     bmosummissi(iii) = sum(temp(:,end));
%   end
%   temp = bostcomm{iii};
%   if ~isempty(temp)
%     bmosumcommi(iii) = sum(temp(:,end));
%   end
% 
% end

keyboard


%% how many MB's LFEs fall in burst wins from 4-s tremors outside ellipse?
for iii = 1: length(trangeout)
  
  [iets,i,j] = indofburst(trangeout,iii);
  
  year = years(iets);
  datesets = dates(floor(dates/1000)==year);
  
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
  
  %Bostock's LFEs NOT in my burst wins inside ellipse
  bostoutbsti = bostoutbst(bostoutbst(:,2)==year & bostoutbst(:,3)==a(1) & bostoutbst(:,4)==a(2),:);
  bostoutbsti = sortrows(bostoutbsti, 5);  
  tbostobi = bostoutbsti(:,5);    
    
  %burst wins outside ellipse
  rangeouttemp = trangeout(trangeout(:,1)==datesets(i), :);
  tst = rangeouttemp(j,2); % start and end time of bursts
  ted = rangeouttemp(j,3);
  
  %which bostock's LFEs are within burst wins outside ellipse?
  indobi = find(tbostobi>=tst & tbostobi<=ted);
  bostob{iii} = bostoutbsti(indobi,:); %all bostocks's lfes inside burst wins (and inside region)

end
bostoba = cat(1,bostob{:});
fracoba = length(bostoba)/length(bostoutbst)
fracoba2 = length(bostoba)/length(bostdayia)

%% how many MB's LFEs fall in ANY 4-s tremor detecting window
bostobiit = [];
bostobiot = [];
for i = 1: nday
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
  
  %Bostock's LFEs NOT in my burst wins inside ellipse
  bostoutbsti = bostoutbst(bostoutbst(:,2)==year & bostoutbst(:,3)==a(1) & bostoutbst(:,4)==a(2),:);
  bostoutbsti = sortrows(bostoutbsti, 5);  
  tbostobi = bostoutbsti(:,5);
  
  hfdayi = hfbnd(hfbnd(:,daycol)==date, :);  % inside bound of the day
  hfdayo = hfout(hfout(:,daycol)==date, :);  % outside bound of the day
  tcnti = hfdayi(:, 15);  % the center of detecting win is the 15th col
  tcnto = hfdayo(:, 15);  % the center of detecting win is the 15th col

  indobiit = [];
  indobiot = [];
  for j = 1: length(tbostobi)
    if ~isempty(find(abs(tbostobi(j)-tcnti)<=2, 1))
      indobiit = [indobiit; j];
    end      
    if ~isempty(find(abs(tbostobi(j)-tcnto)<=2, 1))
      indobiot = [indobiot; j];
    end  
  end
  bostobiit = [bostobiit; bostoutbsti(unique(indobiit), :)];
  bostobiot = [bostobiot; bostoutbsti(unique(indobiot), :)];  
end
fracobiit = length(bostobiit)/length(bostoutbst)
fracobiot = length(bostobiot)/length(bostoutbst)
fracobit = fracobiit+fracobiot


%% plot of comparing times of LFEs from different catalog
for i = 1: length(dates)
  
  date = dates(i);
  year = floor(date/1000);
  jday = floor(date-year*1000);
  a = jul2dat(year,jday);
  dy = num2str(a(2));
  yr = num2str(a(3));

  %Bostock's LFE catalog on the same date
  bostdayi = bostdayia(bostdayia(:,2)==year & bostdayia(:,3)==a(1) & bostdayia(:,4)==a(2),:);
  bostdayi = sortrows(bostdayi, 5); %MB 002/246 LFEs of a certain date
  bostoutbsti = bostoutbst(bostoutbst(:,2)==year & bostoutbst(:,3)==a(1) & bostoutbst(:,4)==a(2),:);
  bostoutbsti = sortrows(bostoutbsti, 5); %MB 002/246 LFEs, same date, but outside bursts
  bost = bosta(bosta(:,2)==year & bosta(:,3)==a(1) & bosta(:,4)==a(2),:);
  bost = sortrows(bost, 5); %MB 002/246 LFEs, same date, inside bursts
  bostmiss = bostmissa(bostmissa(:,2)==year & bostmissa(:,3)==a(1) & bostmissa(:,4)==a(2),:);
  bostmiss = sortrows(bostmiss, 5);
  
  %bursts and 4-s detections of the same day
  rangetemp = trange(trange(:,1)==date, :);
  hfdayi = hfbnd(hfbnd(:,daycol)==date, :);  % inside bound of the day
  hfdayo = hfout(hfout(:,daycol)==date, :);  % outside bound of the day
  
  %bursts and 4-s detections of the same day outside the ellipse
  rangeouttemp = trangeout(trangeout(:,1)==date, :);
  
  tmaxi = hfdayi(:, seccol); % starting time of max power rate of half sec inside the ellipse
  tmaxo = hfdayo(:, seccol); % starting time of max power rate of half sec outside the ellipse
  tcnti = hfdayi(:, 15);  % the center of detecting win is the 15th col
  tcnto = hfdayo(:, 15);  % the center of detecting win is the 15th col
  tbosti = bostdayi(:,5); % (peak arrival) time of Bostock's LFE catalog inside the rectangle
  tbostobi = bostoutbsti(:,5);    
  tbost = bost(:,5);
  tbostmiss = bostmiss(:,5);
  
  %Bostock's LFEs of other fams 
  bostdayofi = bostdayofa(bostdayofa(:,2)==year & bostdayofa(:,3)==a(1) & bostdayofa(:,4)==a(2),:);
  bostdayofi = sortrows(bostdayofi, 5);
  tbostof = bostdayofi(:,5);
  %compute the temporal density
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
  
  widin = 15;  % maximum width allowed is 8.5 inches
  htin = 9;   % maximum height allowed is 11 inches
  nrow = 6;
  ncol = 1;
  f = initfig(widin,htin,nrow,ncol); %initialize fig
  xran = [0.06 0.96]; yran = [0.06 0.98];
  xsep = 0.05; ysep = 0.04;
  optaxpos(f,nrow,ncol,xran,yran,xsep,ysep);

  for j = 1: nrow
    ax=f.ax(j);
    hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    ylim(ax,[0 1]);
    
    %%%plot 4-s tremor detections outside ellipse
    barst = tcnto-2;
    bared = tcnto+2;
    yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.1;
    for ii = 1: length(tcnto)
      p1=plot(ax,[barst(ii) bared(ii)]/3600,[yloc yloc],'-','linew',2.5,'color',[.6 .6 .6]);
    end
    %%%plot auto-determined tremor burst windows from 4-s tremors outside ellipse
    yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.2;
    for ii = 1: size(rangeouttemp,1)
      tst = rangeouttemp(ii,2); % start and end time of bursts
      ted = rangeouttemp(ii,3);
      p2=plot(ax,[tst ted]/3600,[yloc yloc],'c-','linew',2.5);
    end
    %%%plot 4-s tremor detections inside ellipse
    barst = tcnti-2;
    bared = tcnti+2;
    yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.3;
    for ii = 1: length(tcnti)
      p3=plot(ax,[barst(ii) bared(ii)]/3600,[yloc yloc],'k-','linew',2.5);
    end
    %%%plot auto-determined tremor burst windows, same windows for deconvolution
    yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.4;
    for ii = 1: size(rangetemp,1)
      tst = rangetemp(ii,2); % start and end time of bursts
      ted = rangetemp(ii,3);
      indtmaxi = find(tmaxi>=tst-0.1 & tmaxi<=ted+0.1);
      tstbuf = min(tcnti(indtmaxi)-2);
      tedbuf = max(tcnti(indtmaxi)+2);
      p4=plot(ax,[tstbuf tedbuf]/3600,[yloc yloc],'b-','linew',2.5);
    end
    
    %%%plot MB's LFEs of the same dates, but outside my bursts
    yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.55;
    p5=scatter(ax,tbostobi/3600, yloc*ones(size(tbostobi)),5,'r','linew',0.5);
    %%%plot MB's LFEs fall into my burst windows
    yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.65;
    p6=scatter(ax,tbost/3600, yloc*ones(size(tbost)),6,'r','filled');
%     %%%plot MB's LFEs common with mine
%     bostcomm = bostcomma(bostcomma(:,2)==year & bostcomma(:,3)==a(1) & bostcomma(:,4)==a(2),:);
%     bostcomm = sortrows(bostcomm, 5);
%     tbostcomm = bostcomm(:,5);
%     yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.85;
%     p6=scatter(ax,tbostcomm/3600, yloc*ones(size(tbostcomm)),12,'g','filled');
    %%%plot MB's LFEs inside bursts but I missed
    yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.75;
    p7=scatter(ax,tbostmiss/3600, yloc*ones(size(tbostmiss)),6,'g','filled');
    %%%plot MB's LFEs of other fams, in terms of density in 1-win windows
    yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.92;
    p8=scatter(ax,tempden(:,1)/sps/3600, yloc*ones(size(tempden(:,1))),20,tempden(:,2),'filled');%,'MarkerEdgeColor','k'
    colormap(ax,flipud(colormap(ax,'gray')));
    xlim(ax,[(j-1)*24/nrow j*24/nrow]);
    longticks(ax,5);
    nolabels(ax,2);
    if j==1
      legend(ax,[p1 p2 p3 p4 p5 p6 p7 p8],'4-s tremors out ell','burst wins out ell',...
        '4-s tremors in ell','burst wins','MB LFEs out bursts','MB LFEs in bursts',...
        'MB in-burst LFEs I missed','# MB LFEs of other fams per min',...
        'location','best');
    end
    if j==nrow
      xlabel(sprintf('Time (hr) on %d %d %d, %d',a(3),a(1),a(2),date),'FontSize',12);
    end
    
  end
  orient(f.fig,'landscape');
  f1name{i,1} = strcat('mblfes',num2zeropadstr(i,2),'.pdf');
  print(f.fig,'-dpdf','-fillpage',strcat('/home/chaosong/Pictures/',f1name{i,1}));
  
end  

%% merge all figures into a single pdf file
status = system('rm -f /home/chaosong/Pictures/mblfessum.pdf');
for i = 1:size(f1name,1)
  fname{i} = strcat('/home/chaosong/Pictures/',f1name{i});
end
append_pdfs('/home/chaosong/Pictures/mblfessum.pdf',fname);





