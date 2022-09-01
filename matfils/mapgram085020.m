%Sits at one spot (rotations; delays) and looks for the cross-correlation.
%Plots seismograms under locations.
clear all
close all
format short e

%timoffrot=[2005 254 05 19 48 +085 +020  80 115  55];   
%timoffrot=[2005 255 05 19 48 +085 +020  80 115  55];   
%timoffrot=[2003 063 05 19 48 +085 +020  80 115  55];   %2003 needs larger amplitudes, different passband.
timoffrot=[2004 196 05 19 48 +085 +020  80 115  55];   
%timoffrot=[2004 197 05 19 48 +085 +020  80 115  55];   
WAVES=[int2str(timoffrot(1)),'.',int2str(timoffrot(2)),'.',int2str(timoffrot(6)),'.', ...
    int2str(timoffrot(7)),'.',int2str(timoffrot(8)),'.',int2str(timoffrot(9)),'.',int2str(timoffrot(10))]
hi=6.;
lo=1.5;
% hi=8.;
% lo=2.;
year=timoffrot(1);
YEAR=int2str(year);
jday=timoffrot(2);
if jday <= 9
    JDAY=['00',int2str(jday)];
elseif jday<= 99
    JDAY=['0',int2str(jday)];
else
    JDAY=int2str(jday)
end
timlook=3600*timoffrot(3)+60*timoffrot(4)+timoffrot(5)
rotPGC=pi*(timoffrot(8)-90)/180;
rotSSIB=pi*(timoffrot(9)-90)/180;
rotSILB=pi*(timoffrot(10)-90)/180;
SSIBsoff=timoffrot(6);
SILBsoff=timoffrot(7);
SSIBtoff=SSIBsoff/40;
SILBtoff=SILBsoff/40;
IDENTIF=[YEAR,'.',JDAY,'.',int2str(timoffrot(6)),'.',int2str(timoffrot(7)), ...
    '.',int2str(timoffrot(8)),'.',int2str(timoffrot(9)),'.',int2str(timoffrot(10))]

%Read Armbruster's detections:
ArmCat=load(['/data2/arubin/CNDC/',YEAR,'/list.',YEAR,'.pgsssiMAJ']);
%ArmCat=load('/data2/arubin/CNDC/',YEAR,'/list.',YEAR,'.pgsssiMIN');
detects=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
vectoff=[ArmCat(:,5)-SSIBsoff ArmCat(:,6)-SILBsoff (ArmCat(:,6)-SILBsoff)-(ArmCat(:,5)-SSIBsoff)];
distoff=vectoff.*vectoff;
distoff2=sum(distoff,2);
gridoff=max(abs(vectoff),[],2); %Try this instead.  Should change w/mshift (defined later, though)
m0=0;
Detects=0;
for n=1:length(detects)
    %if distoff2(n)<=9
    if gridoff(n)<=10
%      if isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff+1]) 
         %ArmCat(n,5:6)
         m0=m0+1;
         Detects(m0)=detects(n);
     end
end
%abs(xmaxPGSIn(n)-xmaxPGSSn(n)+xmaxSISSn(n))>loopoffma
ArmCat=load(['/data2/arubin/CNDC/',YEAR,'/list.',YEAR,'.pgsssiMAJ_new']);
detects2_8=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
vectoff=[ArmCat(:,5)-SSIBsoff ArmCat(:,6)-SILBsoff];
distoff=vectoff.*vectoff;
distoff2=sum(distoff,2);
gridoff=max(abs(vectoff),[],2); %Try this instead.  Should change w/mshift (defined later, though)
m1=0;
Detects2_8=0;
for n=1:length(detects2_8)
    %if distoff2(n)<=9
    if gridoff(n)<=10
%      if isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff+1]) 
         m1=m1+1;
         Detects2_8(m1)=detects2_8(n);
     end
end

%Read data:
%direc=[YEAR,'/SEP/'];
direc=[YEAR,'/JUL/'];
%direc=[YEAR,'/MAR/'];
prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];
PGCEdat=[prename,'.PGC..BHE.D.SAC'];
PGCNdat=[prename,'.PGC..BHN.D.SAC'];
SSIBEdat=[prename,'.SSIB..HHE.D.SAC'];
SSIBNdat=[prename,'.SSIB..HHN.D.SAC'];
SILBEdat=[prename,'.SILB..HHE.D.SAC'];
SILBNdat=[prename,'.SILB..HHN.D.SAC'];

[PGCE,HdrDataPGC,tnuPGC,pobjPGC,timsPGC]=readsac(PGCEdat,0,'l');
[PGCN,~,~,~,~]=readsac(PGCNdat,0,'l');
[SSIBE,HdrDataSSIB,tnuSSIB,pobjSSIB,timsSSIB]=readsac(SSIBEdat,0,'l');
[SSIBN,~,~,~,~]=readsac(SSIBNdat,0,'l');
[SILBE,HdrDataSILB,tnuSILB,pobjSILB,timsSILB]=readsac(SILBEdat,0,'l');
[SILBN,~,~,~,~]=readsac(SILBNdat,0,'l');

tracelenPGC=length(PGCE);
tracelenSSIB=length(SSIBE);
tracelenSILB=length(SILBE);
    timPGCfirst=timsPGC(1);
    timSSIBfirst=timsSSIB(1); %0.01s earlier than PGC
    timSILBfirst=timsSILB(1); %0.01s earlier than PGC
    timPGClast=timsPGC(tracelenPGC);
    timSSIBlast=timsSSIB(tracelenSSIB);
    timSILBlast=timsSILB(tracelenSILB);
    timstart=timPGCfirst; timend=timPGClast; 
    timwin=(timPGClast-timPGCfirst)/2;
    startdiff=timPGCfirst-timSSIBfirst
    enddiff=timSSIBlast-timPGClast

%cosine taper before filtering:
x=(0:pi/80:pi/2-pi/80)';
% %Seems to be necessary at the start of each day for PGCE:
    PGCE(1:80)=0.;
    PGCE(81:120)=sin(x).*PGCE(81:120); %Only at start of day!
%PGCE(1:40)=sin(x).*PGCE(1:40);
PGCN(1:40)=sin(x).*PGCN(1:40);
x=flipud(x);
PGCE(tracelenPGC-39:tracelenPGC)=sin(x).*PGCE(tracelenPGC-39:tracelenPGC);
PGCN(tracelenPGC-39:tracelenPGC)=sin(x).*PGCN(tracelenPGC-39:tracelenPGC);
x=(0:pi/200:pi/2-pi/200)';
SSIBE(1:100)=sin(x).*SSIBE(1:100);
SSIBN(1:100)=sin(x).*SSIBN(1:100);
SILBE(1:100)=sin(x).*SILBE(1:100);
SILBN(1:100)=sin(x).*SILBN(1:100);
x=flipud(x);
SSIBE(tracelenSSIB-99:tracelenSSIB)=sin(x).*SSIBE(tracelenSSIB-99:tracelenSSIB);
SSIBN(tracelenSSIB-99:tracelenSSIB)=sin(x).*SSIBN(tracelenSSIB-99:tracelenSSIB);
SILBE(tracelenSILB-99:tracelenSSIB)=sin(x).*SILBE(tracelenSILB-99:tracelenSILB);
SILBN(tracelenSILB-99:tracelenSSIB)=sin(x).*SILBN(tracelenSILB-99:tracelenSILB);

%Filter data:
npo=2;
npa=1;
[PGCEf]=1.6e-4*bandpass(PGCE,40,lo,hi,npo,npa,'butter');
[PGCNf]=1.6e-4*bandpass(PGCN,40,lo,hi,npo,npa,'butter');
%[PGCZf]=1.6e-4*bandpass(PGCZ,40,lo,hi,npo,npa,'butter');
[SSIBEf]=4.0e-3*bandpass(SSIBE,100,lo,hi,npo,npa,'butter'); %times 5 for early 2003
[SSIBNf]=4.0e-3*bandpass(SSIBN,100,lo,hi,npo,npa,'butter'); %times 5 for early 2003
[SILBEf]=4.0e-3*bandpass(SILBE,100,lo,hi,npo,npa,'butter'); %times 5 for early 2003
[SILBNf]=4.0e-3*bandpass(SILBN,100,lo,hi,npo,npa,'butter'); %times 5 for early 2003

SSIBEfd = resample(SSIBEf,2,5);
SSIBNfd = resample(SSIBNf,2,5);
SILBEfd = resample(SILBEf,2,5);
SILBNfd = resample(SILBNf,2,5);

PGC=PGCEf+1i*PGCNf;
SSIB=SSIBEfd+1i*SSIBNfd;
SILB=SILBEfd+1i*SILBNfd;
PGCrot=PGC*exp(1i*rotPGC);
SSIBrot=SSIB*exp(1i*rotSSIB);
SILBrot=SILB*exp(1i*rotSILB);

if SSIBtoff > 0
    SSIBrot(1:tracelenPGC-SSIBsoff)=SSIBrot(SSIBsoff+1:tracelenPGC);
    SSIBrot(tracelenPGC-SSIBsoff+1:tracelenPGC)=0;
else
    SSIBrot(-SSIBsoff+1:tracelenPGC)=SSIBrot(1:tracelenPGC+SSIBsoff);
    SSIBrot(1:-SSIBsoff)=0;
end
if SILBtoff > 0
    SILBrot(1:tracelenPGC-SILBsoff)=SILBrot(SILBsoff+1:tracelenPGC);
    SILBrot(tracelenPGC-SILBsoff+1:tracelenPGC)=0;
else
    SILBrot(-SILBsoff+1:tracelenPGC)=SILBrot(1:tracelenPGC+SILBsoff);
    SILBrot(1:-SILBsoff)=0;
end

realPGC=real(PGCrot);
realSSIB=real(SSIBrot);
realSILB=real(SILBrot);

plotstart=1; %(timlook-timstart)*40
plotend=tracelenPGC; %(timlook-timstart)*40
ltmax=max([real(PGCrot(plotstart:plotend)); ...
          real(SSIBrot(plotstart:plotend)); ...
          real(SILBrot(plotstart:plotend))]);
ltmax=min(ltmax,1.);
ltmin=min([real(PGCrot(plotstart:plotend)); ...
          real(SSIBrot(plotstart:plotend)); ...
          real(SILBrot(plotstart:plotend))]);
ltmin=max(ltmin,-1.);

%%%%%%%%%
rotmap=(-28./180.)*pi; %this is for map view
rotprop=(22./180.)*pi; %this is for projections

% IDENTIF4='map2003.061.85.20.80.115.55_2-8-4s';
% IDENTIF128='map2003.061.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER=0;
% tees(1,:)=[0 12000];
% tees(2,:)=[12000 21000];
% tees(3,:)=[21000 25000];
% tees(4,:)=[25000 32000];
% tees(5,:)=[32000 53000];
% tees(6,:)=[53000 61800];
% tees(7,:)=[61800 72000];
% tees(8,:)=[72000 82000];
% tees(9,:)=[82000 86400];
% tees(10,:)=[21000 86400];
% ntees=10;
% % IDENTIF4='map2003.062.85.20.80.115.55_2-8-4s';
% % IDENTIF128='map2003.062.85.20.80.115.55_2-6-128s';
% % ALL128='all.85.20.128s.fort.13.sor';
% % OLDER128='map2003.061.85.20.80.115.55_2-6-128s';
% % OLDER=1;
% % tees(1,:)=[12000 16000];
% % tees(2,:)=[16000 19900];
% % tees(3,:)=[16099 16492];
% % tees(4,:)=[16419 16883];
% % tees(5,:)=[17763 17947];
% % tees(6,:)=[18787 19195];
% % tees(7,:)=[20800 25000];
% % tees(8,:)=[22147 22595];
% % tees(9,:)=[22746 23450];
% % tees(10,:)=[25091 25788];
% % tees(11,:)=[25000 30800];
% % tees(12,:)=[26049 26769];
% % tees(13,:)=[29023 29042];
% % tees(14,:)=[30000 42000];
% % tees(15,:)=[31331 31545];
% % tees(16,:)=[32491 33130];
% % tees(17,:)=[34131 34459];
% % tees(18,:)=[35691 36227];
% % tees(19,:)=[37571 37739];
% % tees(20,:)=[38339 38891];
% % tees(21,:)=[40597 41243];
% % tees(22,:)=[41243 41523];
% % tees(23,:)=[43129 43540];
% % tees(24,:)=[43540 46347];
% % tees(25,:)=[46347 46763];
% % tees(26,:)=[46763 52195];
% % tees(27,:)=[52195 52591];
% % tees(28,:)=[52591 54891];
% % tees(29,:)=[54891 55115];
% % tees(30,:)=[55115 57563];
% % tees(31,:)=[57563 58019];
% % tees(32,:)=[58019 61595];
% % tees(33,:)=[61595 61753];
% % tees(34,:)=[61753 63113];
% % tees(35,:)=[63113 63457];
% % tees(36,:)=[63457 68403];
% % tees(37,:)=[68403 69011];
% % tees(38,:)=[69011 72147];
% % tees(39,:)=[72147 72443];
% % tees(40,:)=[72443 75059];
% % tees(41,:)=[75059 76469];
% % tees(42,:)=[76971 77683];
% % tees(43,:)=[77683 81943];
% % tees(44,:)=[81943 82776];
% % tees(45,:)=[84200 85200];
% % tees(46,:)=[0 12000];
% % ntees=46;
% IDENTIF4='map2003.063.85.20.80.115.55_2-8-4s';
% IDENTIF128='map2003.063.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER128='map2003.061-062.85.20.80.115.55_2-6-128s';
% OLDER=1;
% tees(1,:)=[6152 7355];
% tees(2,:)=[8671 8788];
% tees(3,:)=[8901 9091];
% tees(4,:)=[10397 10475];
% tees(5,:)=[12323 12627];
% tees(6,:)=[13981 14476];
% tees(7,:)=[21403 22576];
% tees(8,:)=[25398 25920];
% tees(9,:)=[31853 32411];
% tees(10,:)=[32025 33741];
% tees(11,:)=[37908 38378];
% tees(12,:)=[40177 41000];
% tees(13,:)=[40702 42166];
% tees(14,:)=[42707 43154];
% tees(15,:)=[63683 64095];
% tees(16,:)=[64935 65534];
% tees(17,:)=[65534 65885];
% tees(18,:)=[65789 66487];
% tees(19,:)=[67896 68091];
% tees(20,:)=[69456 69892];
% tees(21,:)=[78919 79308];
% tees(22,:)=[79987 81280];
% tees(23,:)=[84456 84900];
% tees(24,:)=[0 86400];
% ntees=24;

% IDENTIF4='map2004.196.85.20.80.115.55_2-6-4s';
% IDENTIF128='map2004.196.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER=0;
% tees(1,:)=[16000 27877];
% tees(2,:)=[27877 27971];
% tees(3,:)=[35087 35277]; 
% tees(4,:)=[36403 37614];
% tees(5,:)=[37788 38971];
% tees(6,:)=[38571 38755];
% tees(7,:)=[38891 39132]; % tees(5,:)=[38971 40179];
% tees(8,:)=[39115 39287];
% tees(9,:)=[39467 40179];
% tees(10,:)=[39955 40197];
% tees(11,:)=[40373 40477];
% tees(12,:)=[40487 40503];
% tees(13,:)=[40899 41395];
% tees(14,:)=[41645 42441];
% tees(15,:)=[42715 43187];
% tees(16,:)=[43621 43751];
% tees(17,:)=[44227 44539];
% tees(18,:)=[44715 44827];
% tees(19,:)=[45275 45403];
% tees(20,:)=[46264 46396];
% tees(21,:)=[46304 47023];
% tees(22,:)=[47043 47188];
% tees(23,:)=[48375 49034];
% tees(24,:)=[49619 49763];
% tees(25,:)=[51035 52363];
% tees(26,:)=[51035 51700];
% tees(27,:)=[51035 51442];
% tees(28,:)=[51442 51555];
% tees(29,:)=[54771 54907];
% tees(30,:)=[56075 56643];
% tees(31,:)=[56819 57603];
% tees(32,:)=[62331 62576];
% tees(33,:)=[62576 63683];
% tees(34,:)=[64807 65280];
% tees(35,:)=[69098 69459];
% tees(36,:)=[69891 70115];
% tees(37,:)=[70395 70854];
% tees(38,:)=[73855 74566];
% tees(39,:)=[75915 76027];
% tees(40,:)=[78827 79122];
% tees(41,:)=[81539 81860];
% tees(42,:)=[85841 86379];
% tees(43,:)=[16000 86400];
% ntees=43;
IDENTIF4='map2004.197.85.20.80.115.55_2-6-4s';
IDENTIF128='map2004.197.85.20.80.115.55_2-6-128s';
ALL128='all.85.20.128s.fort.13.sor';
OLDER128='map2004.196.85.20.80.115.55_2-6-128s';
OLDER=1;
tees(1,:)=[5384 5742];
tees(2,:)=[8265 8797];   %tees(2,:)=[8265 8500];
tees(3,:)=[13055 14218];
tees(4,:)=[13867 14218];
tees(5,:)=[14787 14931];
tees(6,:)=[15883 16475];
tees(7,:)=[16083 16475];
tees(8,:)=[16450 17000];
tees(9,:)=[17363 17775];
tees(10,:)=[18359 18765];
tees(11,:)=[19678 20382];
tees(12,:)=[21703 22003];
tees(13,:)=[23862 24575];
tees(14,:)=[24043 24427];
tees(15,:)=[24939 25251];
tees(16,:)=[25227 26091];
tees(17,:)=[27431 28673];
tees(18,:)=[28531 28673];
tees(19,:)=[30232 31013];
tees(20,:)=[30621 31644]; % tees(13,:)=[30621 31013];
tees(21,:)=[36790 38587]; % tees(14,:)=[36790 38071];
tees(22,:)=[37373 38587];
tees(23,:)=[39195 39600];
tees(24,:)=[42959 43720];
tees(25,:)=[45211 45323];
tees(26,:)=[45879 46059];
tees(27,:)=[46107 46528];
tees(28,:)=[55874 57049];
tees(29,:)=[69603 71128];
tees(30,:)=[69603 70700];
tees(31,:)=[71312 71424];
tees(32,:)=[77128 78097]; %possible split?
tees(33,:)=[78295 78850]; %79330];
tees(34,:)=[79347 79387];
tees(35,:)=[83411 83917];
tees(36,:)=[84379 85059]; %possible split?
tees(37,:)=[0 86400];
ntees=37;
% IDENTIF4='map2004.198.85.20.80.115.55_2-6-4s';
% IDENTIF128='map2004.198.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER128='map2004.196-197.85.20.80.115.55_2-6-128s';
% OLDER=1;
% tees(1,:)=[9200 11000];
% tees(2,:)=[11000 12000];
% tees(3,:)=[0 86400];   
% ntees=3;
% % IDENTIF4='map2004.199.85.20.80.115.55_2-6-4s';
% % IDENTIF128='map2004.198.85.20.80.115.55_2-6-128s';
% % ALL128='all.85.20.128s.fort.13.sor';
% % OLDER128='map2004.196-197.85.20.80.115.55_2-6-128s';
% % OLDER=1;
% % tees(1,:)=[0 2000];
% % tees(2,:)=[2000 4000];
% % tees(3,:)=[4200 4300];   
% % tees(4,:)=[4750 5000];   
% % tees(5,:)=[6150 6300];   
% % tees(6,:)=[6300 8000];   
% % tees(7,:)=[8550 10000];   
% % tees(8,:)=[0 10000];   
% % ntees=8;

% % IDENTIF4='map2005.253.85.20.80.115.55_2-6-4s';
% % IDENTIF128='map2005.253.85.20.80.115.55_2-6-128s';
% % ALL128='all.85.20.128s.fort.13.sor';
% % OLDER=0;
% % tees(1,:)=[0 86400];
% % ntees=1;
% IDENTIF4='map2005.254.85.20.80.115.55_2-6-4s';
% IDENTIF128='map2005.254.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER=0;
% tees(1,:)=[0 14000];
% tees(2,:)=[14000 16500];
% tees(3,:)=[14000 18000];
% tees(4,:)=[14100 14400];
% tees(5,:)=[14600 14800];
% tees(6,:)=[14800 15000];
% tees(7,:)=[15000 15300];
% tees(8,:)=[18090 18261];
% tees(9,:)=[21106 21348];
% tees(10,:)=[24330 24596];
% tees(11,:)=[27426 28044];
% tees(12,:)=[30342 31676];
% tees(13,:)=[30342 31140];
% tees(14,:)=[31434 31676];
% tees(15,:)=[34498 49900];
% tees(16,:)=[34498 34942];
% tees(17,:)=[35202 35916];
% tees(18,:)=[48754 49932];
% tees(19,:)=[54226 54972];
% tees(20,:)=[68850 69116];
% tees(21,:)=[73372 73788];
% tees(22,:)=[79970 81020];
% tees(23,:)=[82654 83014];
% tees(24,:)=[0 86400];
% ntees=24;
% % IDENTIF4='map2005.255.85.20.80.115.55_2-6-4s';
% % IDENTIF128='map2005.255.85.20.80.115.55_2-6-128s';
% % ALL128='all.85.20.128s.fort.13.sor';
% % OLDER128='map2005.254.85.20.80.115.55_2-6-128s';
% % OLDER=1;
% % tees(1,:)=[803 1083];
% % tees(2,:)=[4835 6060];
% % tees(3,:)=[8069 8146];
% % tees(4,:)=[8069 8322];
% % tees(5,:)=[18321 19291];
% % tees(6,:)=[22739 23395];
% % tees(7,:)=[26291 26539];
% % tees(8,:)=[34451 36408];
% % tees(9,:)=[38955 39251];
% % tees(10,:)=[43013 43867];
% % tees(11,:)=[45971 46649];
% % tees(12,:)=[48573 48755];
% % tees(13,:)=[50523 51796];
% % tees(14,:)=[53447 53925];
% % tees(15,:)=[54109 54979];
% % tees(16,:)=[58019 59571];
% % tees(17,:)=[61859 62587];
% % tees(18,:)=[63669 63979];
% % tees(19,:)=[67633 68595];
% % tees(20,:)=[75059 75533];
% % tees(21,:)=[0 86400];
% % ntees=21;
% IDENTIF4='map2005.256.85.20.80.115.55_2-6-4s';
% IDENTIF128='map2005.256.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER128='map2005.255.85.20.80.115.55_2-6-128s';
% OLDER=1;
% tees(1,:)=[3000 4500];
% tees(2,:)=[3000 5000];
% tees(3,:)=[12000 16000];
% tees(4,:)=[16000 17200];
% tees(5,:)=[79000 83000];
% tees(6,:)=[80100 80200];
% tees(7,:)=[84000 85900];
% tees(8,:)=[85000 85200];
% tees(9,:)=[0 86400];
% ntees=9;
% % IDENTIF4='map2005.257.85.20.80.115.55_2-6-4s';
% % IDENTIF128='map2005.257.85.20.80.115.55_2-6-128s';
% % ALL128='all.85.20.128s.fort.13.sor';
% % OLDER128='map2005.255.85.20.80.115.55_2-6-128s';
% % OLDER=1;
% % tees(1,:)=[2650 2850];
% % tees(2,:)=[82000 86000];
% % tees(3,:)=[0 86400];
% % ntees=3;

% locs4=load(IDENTIF4);
% locs128=load(IDENTIF128);
% locsAll=load(ALL128);
% if OLDER==1
%     olderlocs=load(OLDER128);
% end
% nin4=length(locs4);
% nin128=length(locs128);
orig4=load(IDENTIF4);
orig128=load(IDENTIF128);
origAll=load(ALL128);
if OLDER==1
    origolder=load(OLDER128);
end
nin4=length(orig4);
nin128=length(orig128);

locs4=0.33*orig4(:,2)+0.6i*orig4(:,3);
locs4=locs4*exp(1i*rotmap);
locs128=0.33*orig128(:,2)+0.6i*orig128(:,3);
locs128=locs128*exp(1i*rotmap);
locsAll=0.33*origAll(:,1)+0.6i*origAll(:,2);
locsAll=locsAll*exp(1i*rotmap);
if OLDER==1
    olderlocs=0.33*origolder(:,2)+0.6i*origolder(:,3);
    olderlocs=olderlocs*exp(1i*rotmap);
end
prop4=0.33*orig4(:,2)+0.6i*orig4(:,3);
prop4=prop4*exp(1i*rotprop);
prop128=0.33*orig128(:,2)+0.6i*orig128(:,3);
prop128=prop128*exp(1i*rotprop);
propAll=0.33*origAll(:,1)+0.6i*origAll(:,2);
propAll=propAll*exp(1i*rotprop);
if OLDER==1
    olderprop=0.33*origolder(:,2)+0.6i*origolder(:,3);
    olderprop=olderprop*exp(1i*rotprop);
end

clo=-5;
chi=2;
dotsz=8;
scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite;
for k=ntees:-1:1
    if orig4(1,1) >= tees(k,1)-1
        istart4=1;
    end
    for i=2:nin4-1
        if orig4(i-1,1) < tees(k,1)-1 && orig4(i,1) >= tees(k,1)-1
            istart4=i;
        end
        if orig4(i+1,1) > tees(k,2)+1 && orig4(i,1) <= tees(k,2)+1
            iend4=i;
        end
    end
    if orig4(nin4,1) <= tees(k,2)-1
        iend4=nin4;
    end
    if orig128(1,1) >= tees(k,1)-1
        istart128=1;
    end
    for i=2:nin128-1
        if orig128(i-1,1) < tees(k,1)-1 && orig128(i,1) >= tees(k,1)-1
            istart128=i;
        end
        if orig128(i+1,1) > tees(k,2)+1 && orig128(i,1) <= tees(k,2)+1
            iend128=i;
        end
    end
    if orig128(nin128,1) <= tees(k,2)-1
        iend128=nin128;
    end
    if orig128(1,1) >= tees(k,2)-1
        iend128=istart128;
    end
    if orig128(nin128,1) <= tees(k,1)-1
        istart128=iend128;
    end

    h=figure('Position',[wid 1 wid hite]);
    %subplot(2,1,1,'align')
    axes('Units','normalized','position',[0.1 0.55 0.5 0.5])
    colormap(jet) 
    plot(real(locsAll),imag(locsAll),'ko','MarkerSize',2,'linewidth',0.2)
    hold on
    plot(real(locs128(1:istart128)),imag(locs128(1:istart128)),'ro','MarkerSize',2,'linewidth',0.2)
    if OLDER==1
        plot(real(olderlocs),imag(olderlocs),'ro','MarkerSize',2,'linewidth',0.2)
    end
    scatter(real(locs128(istart128:iend128)),imag(locs128(istart128:iend128)),90,orig128(istart128:iend128,1)...
        ,'linewidth',1)
    % First on top
    scatter3(real(locs4(istart4:iend4)),imag(locs4(istart4:iend4)),86400-orig4(istart4:iend4,1),25,orig4(istart4:iend4,1),'filled')
    % Last on top
    %scatter3(real(locs4(istart4:iend4)),imag(locs4(istart4:iend4)),orig4(istart4:iend4,1),40,orig4(istart4:iend4,1),'filled')
    axis equal
    axis([-6 5 -6 5])
    xlabel('km E')
    ylabel('km N')
    caxis([tees(k,1) tees(k,2)])
    colorbar
    title(IDENTIF4)
    box on
    
    %subplot(2,1,2,'align')
    axes('position',[0.1 0.2 0.5 0.5])
    colormap(jet) 
    plot(real(locsAll),imag(locsAll),'ko','MarkerSize',2,'linewidth',0.2)
    hold on
    plot(real(locs128(1:istart128)),imag(locs128(1:istart128)),'ro','MarkerSize',2,'linewidth',0.2)
    if OLDER==1
        plot(real(olderlocs),imag(olderlocs),'ro','MarkerSize',2,'linewidth',0.2)
    end
    scatter3(real(locs4(istart4:iend4)),imag(locs4(istart4:iend4)),86400-orig4(istart4:iend4,1),25,log(orig4(istart4:iend4,7)),'filled')
    axis equal
    axis([-6 5 -6 5])
    xlabel('km E')
    ylabel('km N')
    caxis([clo chi])
    %caxis([-4.5 2.5])
    %caxis([min(log(orig4(istart4:iend4,7))) max(log(orig4(istart4:iend4,7)))])
    colorbar
    box on
    
    axes('position',[0.1 0.15 0.8 0.08])
    hold on
    hrf = plotreflinesr(gca,orig4(:,1),'x','r');
    plot(timsPGC,real(PGCrot),'g');
    plot(timsPGC,real(SSIBrot),'b');
    plot(timsPGC,real(SILBrot),'k');
    axis([tees(k,1) tees(k,2) ltmin ltmax]);
    %xlabel(WAVES)
    box on
    
    axes('position',[0.1 0.05 0.8 0.08])
    colormap(jet) 
    hold on
    scatter3(orig4(:,1),real(prop4),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
    axis([tees(k,1) tees(k,2) -5 5])
    ylabel('km along strike')
    caxis([clo chi])
    %colorbar
    box on

    axes('position',[0.7 0.3 0.08/(8/10.5) 0.8*(8/10.5)])
    colormap(jet) 
    hold on
    scatter3(imag(prop4),orig4(:,1),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
    axis([-6 4 tees(k,1) tees(k,2)])
    ylabel('km along dip')
    caxis([clo chi])
    %colorbar
    box on

    set(h,'PaperPosition',[0.25 0.25 8 10.5])
    print(h,'-depsc',['MAPGRAMS/tmp',int2str(k),'.eps'])
      
end

