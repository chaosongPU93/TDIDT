%Sits at one spot (rotations; delays) and looks for the cross-correlation.
clear all
close all
format short e

% % %timoffrot=[2005 254 05 19 48 +085 +020  80 115  55];   
% % %timoffrot=[2005 255 05 19 48 +085 +020  80 115  55];   
% % timoffrot=[2003 062 05 19 48 +085 +020  80 115  55];   %2003 needs larger amplitudes, different passband.
% % %timoffrot=[2003 063 05 19 48 +085 +020  80 115  55];   %2003 needs larger amplitudes, different passband.
% % %timoffrot=[2004 196 05 19 48 +085 +020  80 115  55];   
% % %timoffrot=[2004 197 05 19 48 +085 +020  80 115  55];   
% % WAVES=[int2str(timoffrot(1)),'.',int2str(timoffrot(2)),'.',int2str(timoffrot(6)),'.', ...
% %     int2str(timoffrot(7)),'.',int2str(timoffrot(8)),'.',int2str(timoffrot(9)),'.',int2str(timoffrot(10))]
% % hi=6.;
% % lo=1.5;
% % hi=8.;
% % lo=2.;
% % year=timoffrot(1);
% % YEAR=int2str(year);
% % jday=timoffrot(2);
% % if jday <= 9
% %     JDAY=['00',int2str(jday)];
% % elseif jday<= 99
% %     JDAY=['0',int2str(jday)];
% % else
% %     JDAY=int2str(jday)
% % end
% % timlook=3600*timoffrot(3)+60*timoffrot(4)+timoffrot(5)
% % rotPGC=pi*(timoffrot(8)-90)/180;
% % rotSSIB=pi*(timoffrot(9)-90)/180;
% % rotSILB=pi*(timoffrot(10)-90)/180;
% % SSIBsoff=timoffrot(6);
% % SILBsoff=timoffrot(7);
% % SSIBtoff=SSIBsoff/40;
% % SILBtoff=SILBsoff/40;
% % IDENTIF=[YEAR,'.',JDAY,'.',int2str(timoffrot(6)),'.',int2str(timoffrot(7)), ...
% %     '.',int2str(timoffrot(8)),'.',int2str(timoffrot(9)),'.',int2str(timoffrot(10))]
% % 
% % %Read Armbruster's detections:
% % ArmCat=load(['/data2/arubin/CNDC/',YEAR,'/list.',YEAR,'.pgsssiMAJ']);
% % %ArmCat=load('/data2/arubin/CNDC/',YEAR,'/list.',YEAR,'.pgsssiMIN');
% % detects=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
% % vectoff=[ArmCat(:,5)-SSIBsoff ArmCat(:,6)-SILBsoff (ArmCat(:,6)-SILBsoff)-(ArmCat(:,5)-SSIBsoff)];
% % distoff=vectoff.*vectoff;
% % distoff2=sum(distoff,2);
% % gridoff=max(abs(vectoff),[],2); %Try this instead.  Should change w/mshift (defined later, though)
% % m0=0;
% % Detects=0;
% % for n=1:length(detects)
% %     %if distoff2(n)<=9
% %     if gridoff(n)<=10
% % %      if isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff]) ||...
% % %         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff]) ||...
% % %         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff]) ||...
% % %         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff-1]) ||...
% % %         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff+1]) ||...
% % %         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff-1]) ||...
% % %         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff+1]) ||...
% % %         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff-1]) ||...
% % %         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff+1]) 
% %          %ArmCat(n,5:6)
% %          m0=m0+1;
% %          Detects(m0)=detects(n);
% %      end
% % end
% % %abs(xmaxPGSIn(n)-xmaxPGSSn(n)+xmaxSISSn(n))>loopoffma
% % ArmCat=load(['/data2/arubin/CNDC/',YEAR,'/list.',YEAR,'.pgsssiMAJ_new']);
% % detects2_8=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
% % vectoff=[ArmCat(:,5)-SSIBsoff ArmCat(:,6)-SILBsoff];
% % distoff=vectoff.*vectoff;
% % distoff2=sum(distoff,2);
% % gridoff=max(abs(vectoff),[],2); %Try this instead.  Should change w/mshift (defined later, though)
% % m1=0;
% % Detects2_8=0;
% % for n=1:length(detects2_8)
% %     %if distoff2(n)<=9
% %     if gridoff(n)<=10
% % %      if isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff]) ||...
% % %         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff]) ||...
% % %         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff]) ||...
% % %         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff-1]) ||...
% % %         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff+1]) ||...
% % %         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff-1]) ||...
% % %         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff+1]) ||...
% % %         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff-1]) ||...
% % %         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff+1]) 
% %          m1=m1+1;
% %          Detects2_8(m1)=detects2_8(n);
% %      end
% % end
% % 
% % %Read data:
% % %direc=[YEAR,'/SEPT/'];
% % direc=[YEAR,'/MAR/'];
% % %direc=[YEAR,'/MAR/'];
% % prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];
% % PGCEdat=[prename,'.PGC..BHE.D.SAC'];
% % PGCNdat=[prename,'.PGC..BHN.D.SAC'];
% % SSIBEdat=[prename,'.SSIB..HHE.D.SAC'];
% % SSIBNdat=[prename,'.SSIB..HHN.D.SAC'];
% % SILBEdat=[prename,'.SILB..HHE.D.SAC'];
% % SILBNdat=[prename,'.SILB..HHN.D.SAC'];
% % 
% % [PGCE,HdrDataPGC,tnuPGC,pobjPGC,timsPGC]=readsac(PGCEdat,0,'l');
% % [PGCN,~,~,~,~]=readsac(PGCNdat,0,'l');
% % [SSIBE,HdrDataSSIB,tnuSSIB,pobjSSIB,timsSSIB]=readsac(SSIBEdat,0,'l');
% % [SSIBN,~,~,~,~]=readsac(SSIBNdat,0,'l');
% % [SILBE,HdrDataSILB,tnuSILB,pobjSILB,timsSILB]=readsac(SILBEdat,0,'l');
% % [SILBN,~,~,~,~]=readsac(SILBNdat,0,'l');
% % 
% % tracelenPGC=length(PGCE);
% % tracelenSSIB=length(SSIBE);
% % tracelenSILB=length(SILBE);
% %     timPGCfirst=timsPGC(1);
% %     timSSIBfirst=timsSSIB(1); %0.01s earlier than PGC
% %     timSILBfirst=timsSILB(1); %0.01s earlier than PGC
% %     timPGClast=timsPGC(tracelenPGC);
% %     timSSIBlast=timsSSIB(tracelenSSIB);
% %     timSILBlast=timsSILB(tracelenSILB);
% %     timstart=timPGCfirst; timend=timPGClast; 
% %     timwin=(timPGClast-timPGCfirst)/2;
% %     startdiff=timPGCfirst-timSSIBfirst
% %     enddiff=timSSIBlast-timPGClast
% % 
% % %cosine taper before filtering:
% % x=(0:pi/80:pi/2-pi/80)';
% % % %Seems to be necessary at the start of each day for PGCE:
% %     PGCE(1:80)=0.;
% %     PGCE(81:120)=sin(x).*PGCE(81:120); %Only at start of day!
% % %PGCE(1:40)=sin(x).*PGCE(1:40);
% % PGCN(1:40)=sin(x).*PGCN(1:40);
% % x=flipud(x);
% % PGCE(tracelenPGC-39:tracelenPGC)=sin(x).*PGCE(tracelenPGC-39:tracelenPGC);
% % PGCN(tracelenPGC-39:tracelenPGC)=sin(x).*PGCN(tracelenPGC-39:tracelenPGC);
% % x=(0:pi/200:pi/2-pi/200)';
% % SSIBE(1:100)=sin(x).*SSIBE(1:100);
% % SSIBN(1:100)=sin(x).*SSIBN(1:100);
% % SILBE(1:100)=sin(x).*SILBE(1:100);
% % SILBN(1:100)=sin(x).*SILBN(1:100);
% % x=flipud(x);
% % SSIBE(tracelenSSIB-99:tracelenSSIB)=sin(x).*SSIBE(tracelenSSIB-99:tracelenSSIB);
% % SSIBN(tracelenSSIB-99:tracelenSSIB)=sin(x).*SSIBN(tracelenSSIB-99:tracelenSSIB);
% % SILBE(tracelenSILB-99:tracelenSSIB)=sin(x).*SILBE(tracelenSILB-99:tracelenSILB);
% % SILBN(tracelenSILB-99:tracelenSSIB)=sin(x).*SILBN(tracelenSILB-99:tracelenSILB);
% % 
% % %Filter data:
% % npo=2;
% % npa=1;
% % [PGCEf]=1.6e-4*bandpass(PGCE,40,lo,hi,npo,npa,'butter');
% % [PGCNf]=1.6e-4*bandpass(PGCN,40,lo,hi,npo,npa,'butter');
% % %[PGCZf]=1.6e-4*bandpass(PGCZ,40,lo,hi,npo,npa,'butter');
% % [SSIBEf]=5*4.0e-3*bandpass(SSIBE,100,lo,hi,npo,npa,'butter'); %times 5 for early 2003
% % [SSIBNf]=5*4.0e-3*bandpass(SSIBN,100,lo,hi,npo,npa,'butter'); %times 5 for early 2003
% % [SILBEf]=5*4.0e-3*bandpass(SILBE,100,lo,hi,npo,npa,'butter'); %times 5 for early 2003
% % [SILBNf]=5*4.0e-3*bandpass(SILBN,100,lo,hi,npo,npa,'butter'); %times 5 for early 2003
% % 
% % SSIBEfd = resample(SSIBEf,2,5);
% % SSIBNfd = resample(SSIBNf,2,5);
% % SILBEfd = resample(SILBEf,2,5);
% % SILBNfd = resample(SILBNf,2,5);
% % 
% % PGC=PGCEf+1i*PGCNf;
% % SSIB=SSIBEfd+1i*SSIBNfd;
% % SILB=SILBEfd+1i*SILBNfd;
% % PGCrot=PGC*exp(1i*rotPGC);
% % SSIBrot=SSIB*exp(1i*rotSSIB);
% % SILBrot=SILB*exp(1i*rotSILB);
% % 
% % if SSIBtoff > 0
% %     SSIBrot(1:tracelenPGC-SSIBsoff)=SSIBrot(SSIBsoff+1:tracelenPGC);
% %     SSIBrot(tracelenPGC-SSIBsoff+1:tracelenPGC)=0;
% % else
% %     SSIBrot(-SSIBsoff+1:tracelenPGC)=SSIBrot(1:tracelenPGC+SSIBsoff);
% %     SSIBrot(1:-SSIBsoff)=0;
% % end
% % if SILBtoff > 0
% %     SILBrot(1:tracelenPGC-SILBsoff)=SILBrot(SILBsoff+1:tracelenPGC);
% %     SILBrot(tracelenPGC-SILBsoff+1:tracelenPGC)=0;
% % else
% %     SILBrot(-SILBsoff+1:tracelenPGC)=SILBrot(1:tracelenPGC+SILBsoff);
% %     SILBrot(1:-SILBsoff)=0;
% % end
% % 
% % realPGC=real(PGCrot);
% % realSSIB=real(SSIBrot);
% % realSILB=real(SILBrot);
% % 
% % plotstart=1; %(timlook-timstart)*40
% % plotend=tracelenPGC; %(timlook-timstart)*40
% % ltmax=max([real(PGCrot(plotstart:plotend)); ...
% %           real(SSIBrot(plotstart:plotend)); ...
% %           real(SILBrot(plotstart:plotend))]);
% % ltmax=min(ltmax,1.);
% % ltmin=min([real(PGCrot(plotstart:plotend)); ...
% %           real(SSIBrot(plotstart:plotend)); ...
% %           real(SILBrot(plotstart:plotend))]);
% % ltmin=max(ltmin,-1.);

% IDENTIF4='map2005.259.36.75.80.365.385_2-6-4s';
% IDENTIF128='map2005.259.36.75.80.365.385_2-6-128s';
% tees(1,:)=[45191 47637];
% tees(2,:)=[48995 49285];   
% tees(3,:)=[47840 48500];   
% tees(4,:)=[50145 50485];
% tees(5,:)=[51003 51504];
% tees(6,:)=[54263 54489];
% tees(7,:)=[54473 54554];
% tees(8,:)=[54566 55088];
% tees(9,:)=[56124 56611];
% tees(10,:)=[56124 57070];
% tees(11,:)=[57479 58071];
% tees(12,:)=[60008 61141];
% tees(13,:)=[60601 61141];
% tees(14,:)=[68068 68443];
% tees(15,:)=[70001 70875];
% tees(16,:)=[79017 79688];
% tees(17,:)=[45191 46400];
% tees(18,:)=[46400 47637];
% ntees=18;
% % IDENTIF4='map2005.259.45.72.65.355.370_2-6-4s';
% % IDENTIF128='map2005.259.45.72.65.355.370_2-6-128s';
% % tees(1,:)=[10353 10473];
% % tees(2,:)=[10550 10720];   
% % tees(3,:)=[11027 11905];
% % tees(4,:)=[12509 13063];
% % tees(5,:)=[13807 13920];
% % tees(6,:)=[13993 14441];
% % tees(7,:)=[15545 15802];
% % tees(8,:)=[16907 17275];
% % tees(9,:)=[20065 20773];
% % tees(10,:)=[23057 23475];
% % tees(11,:)=[24435 24635];
% % tees(12,:)=[26164 26377];
% % tees(13,:)=[27099 27403];
% % tees(14,:)=[28081 28315];
% % tees(15,:)=[29731 29883];
% % tees(16,:)=[30206 30875];
% % tees(17,:)=[31579 31841];
% % tees(18,:)=[31834 32318]; 
% % tees(19,:)=[34963 35682];
% % tees(20,:)=[35539 35950];
% % tees(21,:)=[36721 38145];
% % tees(22,:)=[41235 41509]; 
% % tees(23,:)=[43917 44855];
% % tees(24,:)=[44855 45155];
% % tees(25,:)=[46952 47915];
% % tees(26,:)=[50047 53536];
% % tees(27,:)=[58263 58624];
% % tees(28,:)=[62373 63533];
% % tees(29,:)=[63787 64050];
% % tees(30,:)=[64077 64112];
% % tees(31,:)=[64221 64837];
% % tees(32,:)=[66755 67061]; 
% % tees(33,:)=[68620 71250];
% % tees(34,:)=[70755 71075];
% % tees(35,:)=[72731 72859];
% % ntees=35;
% IDENTIF4='map2004.197.86.27.75.120.45_2-6-4s';
% IDENTIF128='map2004.197.86.27.75.120.45_2-6-128s';
% tees(1,:)=[8478 8659];
% tees(2,:)=[13230 14218];  
% tees(3,:)=[14787 14931];
% tees(4,:)=[15883 16475];
% tees(5,:)=[17363 17775];
% tees(6,:)=[18359 18765];
% tees(7,:)=[19699 20382];
% tees(8,:)=[21703 22003];
% tees(9,:)=[23862 24575];
% tees(10,:)=[24939 25251];
% tees(11,:)=[28531 28673];
% tees(12,:)=[30232 31013];
% tees(13,:)=[30621 31013];
% tees(14,:)=[36790 38071];
% tees(15,:)=[42959 43544];
% tees(16,:)=[46107 46414];
% tees(17,:)=[55874 57049];
% tees(18,:)=[69639 70700];
% tees(19,:)=[71312 71424];
% tees(20,:)=[78355 78850]; %79330];
% tees(21,:)=[84443 85059];
% tees(22,:)=[85059 85450]; 
% tees(23,:)=[36790 38587];
% ntees=23;
% % IDENTIF4='map2004.196.86.13.80.115.55_2-6-4s';
% % IDENTIF128='map2004.196.86.13.80.115.55_2-6-128s';
% % tees(1,:)=[27877 27971];
% % tees(2,:)=[35087 35277]; 
% % tees(3,:)=[36403 37614];
% % tees(4,:)=[37788 38971];
% % tees(5,:)=[38971 40179];
% % tees(6,:)=[44227 44539];
% % tees(7,:)=[46304 47023];
% % tees(8,:)=[51035 52363];
% % tees(12,:)=[51035 51600];
% % tees(10,:)=[51035 51700];
% % tees(9,:)=[56819 57603];
% % tees(11,:)=[56819 57350];
% % ntees=12;
% IDENTIF4='MAPS/map2004.197.85.20.80.115.55_2-6-4s';
% IDENTIF128='MAPS/map2004.197.85.20.80.115.55_2-6-128s';
% OLDER128='MAPS/map2004.196.85.20.80.115.55_2-6-128s';
% OLDER=1;
% tees(1,:)=[5384 5742];
% tees(2,:)=[8265 8797];   %tees(2,:)=[8265 8500];
% tees(23,:)=[13055 14175];
% tees(3,:)=[13867 14115];
% tees(4,:)=[16083 16475];
% tees(5,:)=[19678 20382];
% tees(6,:)=[24043 24427];
% tees(7,:)=[24939 25251];
% tees(8,:)=[25227 26091];
% tees(9,:)=[27431 28673];
% tees(10,:)=[30867 31644];
% tees(11,:)=[37373 38587];
% tees(12,:)=[39195 39600];
% tees(13,:)=[42987 43720];
% tees(14,:)=[45211 45323];
% tees(15,:)=[45879 46059];
% tees(16,:)=[46157 46528];
% tees(17,:)=[69603 71128];
% tees(18,:)=[77128 78097]; %possible split?
% tees(19,:)=[78295 78843];
% tees(20,:)=[79347 79387];
% tees(21,:)=[83411 83917];
% tees(22,:)=[84379 84975]; %possible split?
% tees(24,:)=[36790 38587];
% ntees=24;
% % IDENTIF4='MAPS/map2004.196.85.20.80.115.55_2-6-4s';
% % IDENTIF128='MAPS/map2004.196.85.20.80.115.55_2-6-128s';
% % YOUNGER128='MAPS/map2004.197.85.20.80.115.55_2-6-128s';
% % OLDER=0;
% % tees(1,:)=[38571 38755];
% % tees(2,:)=[38891 39132];
% % tees(3,:)=[39115 39287];
% % tees(4,:)=[39467 40179];
% % tees(5,:)=[39955 40197];
% % tees(6,:)=[40373 40477];
% % tees(7,:)=[40487 40503];
% % tees(8,:)=[40899 41395];
% % tees(9,:)=[41645 42441];
% % tees(10,:)=[42715 43187];
% % tees(11,:)=[43621 43751];
% % tees(12,:)=[44227 44539];
% % tees(13,:)=[44715 44827];
% % tees(14,:)=[45275 45403];
% % tees(15,:)=[46264 46396];
% % tees(16,:)=[47043 47188];
% % tees(17,:)=[48375 49034];
% % tees(18,:)=[49619 49763];
% % tees(19,:)=[51035 51442];
% % tees(20,:)=[51442 51555];
% % tees(21,:)=[54771 54907];
% % tees(22,:)=[56075 56643];
% % tees(23,:)=[56819 57572];
% % tees(24,:)=[62331 62576];
% % tees(25,:)=[62576 63683];
% % tees(26,:)=[64807 65280];
% % tees(27,:)=[69098 69459];
% % tees(28,:)=[69891 70115];
% % tees(29,:)=[70395 70854];
% % tees(30,:)=[73855 74566];
% % tees(31,:)=[75915 76027];
% % tees(32,:)=[78827 79122];
% % tees(33,:)=[81539 81860];
% % tees(34,:)=[85841 86379];
% % ntees=34;
%IDENTIF4='map2003.063.86.27.75.125.45_2-8-4s';
%IDENTIF128='map2003.063.86.27.75.125.45_2-8-128s';
%tees(13,:)=[6215 7355];
%tees(14,:)=[40177 41000];
%tees(15,:)=[41250 41403];
%tees(1,:)=[6152 6997];
%tees(2,:)=[8671 8788];
%tees(3,:)=[21403 21953];
%tees(4,:)=[25398 25920];
%tees(5,:)=[31853 32411];
%tees(6,:)=[32806 33239];
%tees(7,:)=[37908 38378];
%tees(8,:)=[40177 41403];
%tees(9,:)=[63683 64095];
%tees(10,:)=[67896 68091];
%tees(11,:)=[69456 69892];
%tees(12,:)=[40597 41271];
%ntees=15;
% % IDENTIF4='map2003.062.86.13.80.110.55_2-8-4s';
% % IDENTIF128='map2003.062.86.13.80.110.55_2-8-128s';
% % tees(1,:)=[16100 16492];
% % tees(2,:)=[16683 16738];
% % tees(3,:)=[17826 17905];
% % tees(4,:)=[18979 19173];
% % tees(5,:)=[22153 22595];
% % tees(6,:)=[22746 23450];
% % tees(7,:)=[25105 25788];
% % tees(8,:)=[26049 26769];
% % tees(9,:)=[29023 29042];
% % tees(10,:)=[32493 33130];
% % tees(11,:)=[35691 36227];
% % tees(12,:)=[40597 41271];
% % tees(13,:)=[41271 41505];
% % tees(14,:)=[46351 46730];
% % tees(15,:)=[63113 63449];
% % tees(16,:)=[75059 76401];
% % tees(17,:)=[76981 77606];
% % ntees=17;
% IDENTIF4='MAPS/map2003.063.85.20.80.115.55_2-8-4s';
% IDENTIF128='MAPS/map2003.063.85.20.80.115.55_2-8-128s';
% OLDER128='MAPS/map2003.062.85.20.80.115.55_2-8-128s';
% OLDER=1;
% tees(1,:)=[6215 7355];
% tees(2,:)=[8901 9091];
% tees(3,:)=[10397 10475];
% tees(4,:)=[12323 12627];
% tees(5,:)=[13981 14476];
% tees(6,:)=[21458 22576];
% tees(7,:)=[25641 25918];
% tees(8,:)=[32025 33741];
% tees(9,:)=[40179 40911];
% tees(10,:)=[40702 42166];
% tees(11,:)=[42707 43154];
% tees(12,:)=[64935 65534];
% tees(13,:)=[65534 65885];
% tees(14,:)=[65789 66487];
% tees(15,:)=[78919 79308];
% tees(16,:)=[79987 81280];
% tees(17,:)=[84456 84900];
% ntees=17;
% % IDENTIF4='MAPS/map2003.062.85.20.80.115.55_2-8-4s';
% % IDENTIF128='MAPS/map2003.062.85.20.80.115.55_2-8-128s';
% % YOUNGER128='MAPS/map2003.063.85.20.80.115.55_2-8-128s';
% % OLDER=0;
% % tees(1,:)=[16099 16387];
% % tees(2,:)=[16419 16883];
% % tees(3,:)=[17763 17947];
% % tees(4,:)=[18787 19195];
% % tees(5,:)=[22147 22536];
% % tees(6,:)=[22923 23115];
% % tees(7,:)=[25091 25353];
% % tees(8,:)=[31331 31545];
% % tees(9,:)=[32491 32946];
% % tees(10,:)=[34131 34459];
% % tees(11,:)=[35787 36131];
% % tees(12,:)=[37571 37739];
% % tees(13,:)=[38339 38891];
% % tees(14,:)=[40597 41070];
% % tees(15,:)=[41243 41523];
% % tees(16,:)=[43129 43540];
% % tees(17,:)=[46347 46763];
% % tees(18,:)=[52195 52591];
% % tees(19,:)=[54891 55115];
% % tees(20,:)=[57563 58019];
% % tees(21,:)=[61595 61753];
% % tees(22,:)=[63113 63457];
% % tees(23,:)=[68403 69011];
% % tees(24,:)=[72147 72443];
% % tees(25,:)=[75067 76469];
% % tees(26,:)=[76971 77683];
% % tees(27,:)=[81943 82776];
% % tees(28,:)=[25600 30800];
% % ntees=28;
% IDENTIF4='map2005.255.86.27.75.125.45_2-6-4s';
% IDENTIF128='map2005.255.86.27.75.125.45_2-6-128s';
% tees(1,:)=[803 1083];
% tees(2,:)=[4835 5749];
% tees(3,:)=[8069 8146];
% tees(4,:)=[8069 8322];
% tees(5,:)=[18321 19291];
% tees(6,:)=[26291 26539];
% tees(7,:)=[34451 36283];
% tees(8,:)=[38955 39251];
% tees(9,:)=[43013 43867];
% tees(10,:)=[45971 46649];
% tees(11,:)=[48573 48755];
% tees(12,:)=[50523 51796];
% tees(13,:)=[53447 53925];
% tees(14,:)=[58019 59420];
% tees(15,:)=[61859 62587];
% tees(16,:)=[67633 68595];
% tees(17,:)=[75059 75533];
% ntees=17;
IDENTIF4='MAPS/map2005.254.85.20.80.115.55_2-6-4s';
IDENTIF128='MAPS/map2005.254.85.20.80.115.55_2-6-128s';
YOUNGER128='MAPS/map2005.255.85.20.80.115.55_2-6-128s';
OLDER=0;
% tees(1,:)=[18090 18261];
% tees(2,:)=[21106 21348];
% tees(3,:)=[24330 24596];
tees(1,:)=[27426 28044];
% tees(5,:)=[30342 31140];
%tees(2,:)=[31434 31676];
tees(2,:)=[34498 34942];
tees(3,:)=[35202 35916];
%tees(4,:)=[48754 49900];
% tees(10,:)=[54226 54972];
% tees(11,:)=[68850 69116];
% tees(12,:)=[73372 73788];
% tees(13,:)=[79970 81020];
% tees(14,:)=[82654 83014];
% tees(16,:)=[14000 16500];
% tees(17,:)=[14000 18000];
ntees=3;
% % % IDENTIF4='map2005.254.85.20.80.115.55_2-6-3s'; %3s not as impressive(?)
% % % IDENTIF128='map2005.254.85.20.80.115.55_2-6-128s';
% % % tees(1,:)=[18090 18261];
% % % tees(2,:)=[21106 21348];
% % % tees(3,:)=[24330 24596];
% % % tees(4,:)=[27426 28044];
% % % tees(5,:)=[30342 31140];
% % % tees(6,:)=[31434 31676];
% % % tees(7,:)=[34498 34942];
% % % tees(8,:)=[35202 35916];
% % % tees(9,:)=[48754 49932];
% % % tees(10,:)=[54226 54972];
% % % tees(11,:)=[68850 69116];
% % % tees(12,:)=[73372 73788];
% % % tees(13,:)=[79970 81020];
% % % tees(14,:)=[82654 83014];
% % % tees(15,:)=[48754 49900];
% % % ntees=15;
locs4=load(IDENTIF4);
locs128=load(IDENTIF128);
if OLDER==1
    olderlocs=load(OLDER128);
else
    youngerlocs=load(YOUNGER128);
end
nin4=length(locs4);
nin128=length(locs128);
%figure
for k=1:ntees
    for i=2:nin4-1
        if locs4(i-1,1) < tees(k,1)-1 && locs4(i,1) >= tees(k,1)-1
            istart4=i;
        end
        if locs4(i+1,1) > tees(k,2)+1 && locs4(i,1) <= tees(k,2)+1
            iend4=i;
        end
    end
    for i=2:nin128-1
        if locs128(i-1,1) < tees(k,1)-1 && locs128(i,1) >= tees(k,1)-1
            istart128=i;
        end
        if locs128(i+1,1) > tees(k,2)+1 && locs128(i,1) <= tees(k,2)+1
            iend128=i;
        end
    end

    colormap(jet) 
    %hold on
    if k==1
        axes('position',[0.03 0.62 0.36 0.36],'FontSize',6,'DataAspectRatio',[1 1 1]) %,'XLim',[-3.4 3.4],'YLim',[-3.4 3.4])
        axis equal
        axis([-3.4 3.4 -3.4 3.4])
        text(2.8, -3, 'a','FontSize',8);
    elseif k==2
        axes('position',[0.34 0.62 0.36 0.36],'FontSize',6,'DataAspectRatio',[1 1 1]) %,'XLim',[-3.4 3.4],'YLim',[-3.4 3.4])
        %position
        axis equal
        axis([-3.4 3.4 -3.4 3.4])
        text(2.8, -3, 'b','FontSize',8);
   elseif k==3
        axes('position',[0.05 0.35 0.21 0.21],'FontSize',6,'DataAspectRatio',[1 1 1]) %,'XLim',[-4 4],'YLim',[-4 4])
        axis equal
        axis([-3.4 3.4 -3.4 3.4])
        text(2.8, -3, 'c','FontSize',8);
%     elseif k==4
%         subplot('position',[0.25 0.45 0.2 0.2],'FontSize',6,'DataAspectRatio',[1 1 1],'XLim',[-4 4],'YLim',[-4 4])
%     elseif k==5
%         subplot('position',[0.25 0.55 0.2 0.2])
    end    
    hold on
    if OLDER==0
        plot(0.33*youngerlocs(:,2),0.6*youngerlocs(:,3),'ro','MarkerSize',1,'linewidth',0.2)
    end
    hold on 
    plot(0.33*locs128(iend128:length(locs128),2),0.6*locs128(iend128:length(locs128),3),'ro','MarkerSize',1,'linewidth',0.2)
    hold on 
    plot(0.33*locs128(1:istart128,2),0.6*locs128(1:istart128,3),'ko','MarkerSize',1,'linewidth',0.2)
    hold on 
    if OLDER==1
        plot(0.33*olderlocs(:,2),0.6*olderlocs(:,3),'ko','MarkerSize',1,'linewidth',0.2)
    end
    if k==3
        scatter(0.33*locs128(istart128:iend128,2),0.6*locs128(istart128:iend128,3),60,locs128(istart128:iend128,1)...
        ,'linewidth',0.6)
    else
        scatter(0.33*locs128(istart128:iend128,2),0.6*locs128(istart128:iend128,3),80,locs128(istart128:iend128,1)...
        ,'linewidth',0.6)
    end
    if k==3
        scatter(0.33*locs4(istart4:iend4,2),0.6*locs4(istart4:iend4,3),10,locs4(istart4:iend4,1),'filled')
    else
        scatter(0.33*locs4(istart4:iend4,2),0.6*locs4(istart4:iend4,3),15,locs4(istart4:iend4,1),'filled')
    end
    timestamp=[int2str(tees(k,1)),'-',int2str(tees(k,2)),' s'];
    text(-3.1, 3, timestamp,'FontSize',6);
    text(-3.1, 2.5, 'Sept 11 2005','FontSize',6);
%     axis equal
%     axis([-4 4 -4 4])
    if k==1 || k==3
        %xlabel('PGC-SILB: KM ESE','FontSize',7)
        ylabel('PGC-SSIB: KM NNE','FontSize',7)
    end
    %xlabel('PGSI - KM ESE')
    %ylabel('PGSS - KM NNE')
    caxis([tees(k,1) tees(k,2)])
    %colorbar
    %title(IDENTIF4)
    box on
    %subplot(2,1,2)
%     axes('position',[0.1 0.07 0.6 0.1])
%     hold on
%     hrf = plotreflinesr(gca,locs4(:,1),'x','r');
%     plot(timsPGC,real(PGCrot),'g');
%     plot(timsPGC,real(SSIBrot),'b');
%     plot(timsPGC,real(SILBrot),'k');
%     axis([tees(k,1) tees(k,2) ltmin ltmax]);
%     xlabel(WAVES)
%     box on
    end
IDENTIF4='MAPS/map2005.255.85.20.80.115.55_2-6-4s';
IDENTIF128='MAPS/map2005.255.85.20.80.115.55_2-6-128s';
OLDER128='MAPS/map2005.254.85.20.80.115.55_2-6-128s';
OLDER=1;
tees(1,:)=[5010 6060]; 
%tees(2,:)=[18760 19291];
%tees(3,:)=[22739 23395];
tees(2,:)=[35046 36408];
%tees(5,:)=[50565 51692];
%tees(6,:)=[54109 54979];
%tees(3,:)=[58474 59571];
%tees(8,:)=[63669 63979];
%tees(9,:)=[54200 54450];
%tees(10,:)=[54450 54979];
%tees(11,:)=[54200 54350];
ntees=2;
locs4=load(IDENTIF4);
locs128=load(IDENTIF128);
if OLDER==1
    olderlocs=load(OLDER128);
else
    youngerlocs=load(YOUNGER128);
end
nin4=length(locs4);
nin128=length(locs128);
%figure
%pbaspect('manual')
%pbaspect([4 1 1])
for k=1:ntees
    for i=2:nin4-1
        if locs4(i-1,1) < tees(k,1)-1 && locs4(i,1) >= tees(k,1)-1
            istart4=i;
        end
        if locs4(i+1,1) > tees(k,2)+1 && locs4(i,1) <= tees(k,2)+1
            iend4=i;
        end
    end
    for i=2:nin128-1
        if locs128(i-1,1) < tees(k,1)-1 && locs128(i,1) >= tees(k,1)-1
            istart128=i;
        end
        if locs128(i+1,1) > tees(k,2)+1 && locs128(i,1) <= tees(k,2)+1
            iend128=i;
        end
    end

    %subplot(2,1,1)
    %axes('position',[0.1 0.25 0.7 0.7])
    colormap(jet) 
    %hold on
    if k==1
        axes('position',[0.26 0.35 0.21 0.21],'FontSize',6,'DataAspectRatio',[1 1 1],'XLim',[-4 4],'YLim',[-4 4])
        text(2.8, -3, 'd','FontSize',8);
    elseif k==2
        axes('position',[0.47 0.35 0.21 0.21],'FontSize',6,'DataAspectRatio',[1 1 1],'XLim',[-4 4],'YLim',[-4 4])
        text(2.8, -3, 'e','FontSize',8);
%     elseif k==3
%         axes('position',[0.05 0.1 0.21 0.21],'FontSize',6,'DataAspectRatio',[1 1 1],'XLim',[-4 4],'YLim',[-4 4])
    end   
    hold on
    if OLDER==0
        plot(0.33*youngerlocs(:,2),0.6*youngerlocs(:,3),'ro','MarkerSize',1,'linewidth',0.2)
    end
    hold on 
    plot(0.33*locs128(iend128:length(locs128),2),0.6*locs128(iend128:length(locs128),3),'ro','MarkerSize',1,'linewidth',0.2)
    plot(0.33*locs128(1:istart128,2),0.6*locs128(1:istart128,3),'ko','MarkerSize',1,'linewidth',0.2)
    if OLDER==1
        plot(0.33*olderlocs(:,2),0.6*olderlocs(:,3),'ko','MarkerSize',1,'linewidth',0.2)
    end
    scatter(0.33*locs128(istart128:iend128,2),0.6*locs128(istart128:iend128,3),60,locs128(istart128:iend128,1)...
        ,'linewidth',0.6)
    scatter(0.33*locs4(istart4:iend4,2),0.6*locs4(istart4:iend4,3),10,locs4(istart4:iend4,1),'filled')
    axis equal
    axis([-3.4 3.4 -3.4 3.4])
    timestamp=[int2str(tees(k,1)),'-',int2str(tees(k,2)),' s'];
    text(-3.1, 3, timestamp,'FontSize',6);
    text(-3.1, 2.5, 'Sept 12 2005','FontSize',6);
%     if k==1
%         %xlabel('PGC-SILB: KM ESE','FontSize',7)
%         ylabel('PGC-SSIB: KM NNE','FontSize',7)
%     end
    %xlabel('PGSI - KM ESE')
    %ylabel('PGSS - KM NNE')
    caxis([tees(k,1) tees(k,2)])
    %colorbar
    %title(IDENTIF4)
    box on
    %subplot(2,1,2)
end

IDENTIF4='MAPS/map2004.196.85.20.80.115.55_2-6-4s';
IDENTIF128='MAPS/map2004.196.85.20.80.115.55_2-6-128s';
YOUNGER128='MAPS/map2004.197.85.20.80.115.55_2-6-128s';
OLDER=0;
% tees(1,:)=[38571 38755];
% tees(2,:)=[38891 39132];
% tees(3,:)=[39115 39287];
% tees(4,:)=[39467 40179];
% tees(5,:)=[39955 40197];
% tees(6,:)=[40373 40477];
% tees(7,:)=[40487 40503];
% tees(8,:)=[40899 41395];
% tees(9,:)=[41645 42441];
% tees(10,:)=[42715 43187];
% tees(11,:)=[43621 43751];
% tees(12,:)=[44227 44539];
% tees(13,:)=[44715 44827];
% tees(14,:)=[45275 45403];
% tees(15,:)=[46264 46396];
% tees(16,:)=[47043 47188];
% tees(17,:)=[48375 49034];
% tees(18,:)=[49619 49763];
% tees(19,:)=[51035 51442];
% tees(20,:)=[51442 51555];
% tees(21,:)=[54771 54907];
tees(1,:)=[56075 56643];
% tees(23,:)=[56819 57572];
% tees(24,:)=[62331 62576];
% tees(25,:)=[62576 63683];
% tees(26,:)=[64807 65280];
% tees(27,:)=[69098 69459];
% tees(28,:)=[69891 70115];
% tees(29,:)=[70395 70854];
% tees(30,:)=[73855 74566];
% tees(31,:)=[75915 76027];
% tees(32,:)=[78827 79122];
% tees(33,:)=[81539 81860];
% tees(34,:)=[85841 86379];
ntees=1;
locs4=load(IDENTIF4);
locs128=load(IDENTIF128);
if OLDER==1
    olderlocs=load(OLDER128);
else
    youngerlocs=load(YOUNGER128);
end
nin4=length(locs4);
nin128=length(locs128);
for k=1:ntees
    for i=2:nin4-1
        if locs4(i-1,1) < tees(k,1)-1 && locs4(i,1) >= tees(k,1)-1
            istart4=i;
        end
        if locs4(i+1,1) > tees(k,2)+1 && locs4(i,1) <= tees(k,2)+1
            iend4=i;
        end
    end
    for i=2:nin128-1
        if locs128(i-1,1) < tees(k,1)-1 && locs128(i,1) >= tees(k,1)-1
            istart128=i;
        end
        if locs128(i+1,1) > tees(k,2)+1 && locs128(i,1) <= tees(k,2)+1
            iend128=i;
        end
    end

% %     colormap(jet) 
% %     hold on
    if k==1
        axes('position',[0.05 0.08 0.21 0.21],'FontSize',6,'DataAspectRatio',[1 1 1],'XLim',[-3.4 3.4],'YLim',[-3.4 3.4])
        text(2.8, -3, 'f','FontSize',8);
%     elseif k==2
%         axes('position',[0.05 0.15 0.21 0.21],'FontSize',6,'DataAspectRatio',[1 1 1],'XLim',[-4 4],'YLim',[-4 4])
%     elseif k==3
%         axes('position',[0.05 0.15 0.21 0.21],'FontSize',6,'DataAspectRatio',[1 1 1],'XLim',[-4 4],'YLim',[-4 4])
    end    
    hold on
    if OLDER==0
        plot(0.33*youngerlocs(:,2),0.6*youngerlocs(:,3),'ro','MarkerSize',1,'linewidth',0.2)
    end
    hold on 
    plot(0.33*locs128(iend128:length(locs128),2),0.6*locs128(iend128:length(locs128),3),'ro','MarkerSize',1,'linewidth',0.2)
    plot(0.33*locs128(1:istart128,2),0.6*locs128(1:istart128,3),'ko','MarkerSize',1,'linewidth',0.2)
    if OLDER==1
        plot(0.33*olderlocs(:,2),0.6*olderlocs(:,3),'ko','MarkerSize',1,'linewidth',0.2)
    end
    scatter(0.33*locs128(istart128:iend128,2),0.6*locs128(istart128:iend128,3),60,locs128(istart128:iend128,1)...
        ,'linewidth',0.6)
    scatter(0.33*locs4(istart4:iend4,2),0.6*locs4(istart4:iend4,3),10,locs4(istart4:iend4,1),'filled')
%     axis equal
%     axis([-4 4 -4 4])
    timestamp=[int2str(tees(k,1)),'-',int2str(tees(k,2)),' s'];
    text(-3.1, 3, timestamp,'FontSize',6);
    text(-3.1, 2.5, 'July 14 2004','FontSize',6);
    if k==1
        xlabel('PGC-SILB: KM ESE','FontSize',7)
        ylabel('PGC-SSIB: KM NNE','FontSize',7)
    end
    %xlabel('PGSI - KM ESE')
    %ylabel('PGSS - KM NNE')
    caxis([tees(k,1) tees(k,2)])
    %colorbar
    %title(IDENTIF4)
    box on
end

IDENTIF4='MAPS/map2003.062.85.20.80.115.55_2-8-4s';
IDENTIF128='MAPS/map2003.062.85.20.80.115.55_2-8-128s';
YOUNGER128='MAPS/map2003.063.85.20.80.115.55_2-8-128s';
OLDER=0;
% tees(1,:)=[16099 16387];
% tees(2,:)=[16419 16883];
% tees(3,:)=[17763 17947];
% tees(4,:)=[18787 19195];
% tees(5,:)=[22147 22536];
% tees(6,:)=[22923 23115];
% tees(7,:)=[25091 25353];
% tees(8,:)=[31331 31545];
% tees(9,:)=[32491 32946];
% tees(10,:)=[34131 34459];
% tees(11,:)=[35787 36131];
% tees(12,:)=[37571 37739];
% tees(13,:)=[38339 38891];
% tees(14,:)=[40597 41070];
% tees(15,:)=[41243 41523];
% tees(16,:)=[43129 43540];
% tees(17,:)=[46347 46763];
% tees(18,:)=[52195 52591];
% tees(19,:)=[54891 55115];
% tees(20,:)=[57563 58019];
% tees(21,:)=[61595 61753];
% tees(22,:)=[63113 63457];
% tees(23,:)=[68403 69011];
% tees(24,:)=[72147 72443];
tees(1,:)=[75067 76469];
tees(2,:)=[76971 77683];
% tees(27,:)=[81943 82776];
% tees(28,:)=[25600 30800];
ntees=2;
locs4=load(IDENTIF4);
locs128=load(IDENTIF128);
if OLDER==1
    olderlocs=load(OLDER128);
else
    youngerlocs=load(YOUNGER128);
end
nin4=length(locs4);
nin128=length(locs128);
%figure
%pbaspect('manual')
%pbaspect([4 1 1])
for k=1:ntees
    for i=2:nin4-1
        if locs4(i-1,1) < tees(k,1)-1 && locs4(i,1) >= tees(k,1)-1
            istart4=i;
        end
        if locs4(i+1,1) > tees(k,2)+1 && locs4(i,1) <= tees(k,2)+1
            iend4=i;
        end
    end
    for i=2:nin128-1
        if locs128(i-1,1) < tees(k,1)-1 && locs128(i,1) >= tees(k,1)-1
            istart128=i;
        end
        if locs128(i+1,1) > tees(k,2)+1 && locs128(i,1) <= tees(k,2)+1
            iend128=i;
        end
    end

    %subplot(2,1,1)
    %axes('position',[0.1 0.25 0.7 0.7])
%     colormap(jet) 
%     hold on
    if k==1
        subplot('position',[0.26 0.08 0.21 0.21],'FontSize',6,'DataAspectRatio',[1 1 1],'XLim',[-3.4 3.4],'YLim',[-3.4 3.4])
        text(2.8, -3, 'g','FontSize',8);
        k
    elseif k==2
        subplot('position',[0.47 0.08 0.21 0.21],'FontSize',6,'DataAspectRatio',[1 1 1],'XLim',[-3.4 3.4],'YLim',[-3.4 3.4])
        text(2.8, -3, 'h','FontSize',8);
        k
%     elseif k==3
%         subplot('position',[0.45 0.15 0.21 0.2],'FontSize',6,'DataAspectRatio',[1 1 1],'XLim',[-4 4],'YLim',[-4 4])
    end    
    if OLDER==0
        plot(0.33*youngerlocs(:,2),0.6*youngerlocs(:,3),'ro','MarkerSize',1,'linewidth',0.2)
    end
    hold on 
    plot(0.33*locs128(iend128:length(locs128),2),0.6*locs128(iend128:length(locs128),3),'ro','MarkerSize',1,'linewidth',0.2)
    plot(0.33*locs128(1:istart128,2),0.6*locs128(1:istart128,3),'ko','MarkerSize',1,'linewidth',0.2)
    if OLDER==1
        plot(0.33*olderlocs(:,2),0.6*olderlocs(:,3),'ko','MarkerSize',1,'linewidth',0.2)
    end
    scatter(0.33*locs128(istart128:iend128,2),0.6*locs128(istart128:iend128,3),60,locs128(istart128:iend128,1)...
        ,'linewidth',0.6)
    scatter(0.33*locs4(istart4:iend4,2),0.6*locs4(istart4:iend4,3),10,locs4(istart4:iend4,1),'filled')
    axis equal
    axis([-3.4 3.4 -3.4 3.4])
    if k==1
        text(2.8, -3, 'g','FontSize',8);
    elseif k==2
        text(2.8, -3, 'h','FontSize',8);
    end
    timestamp=[int2str(tees(k,1)),'-',int2str(tees(k,2)),' s'];
    text(-3.1, 3, timestamp,'FontSize',6);
    text(-3.1, 2.5, 'Mar 03 2003','FontSize',6);
    %if k==1
        xlabel('PGC-SILB: KM ESE','FontSize',7)
        %ylabel('PGSS - KM NNE','FontSize',8)
    %end
    caxis([tees(k,1) tees(k,2)])
    %colorbar
    %title(IDENTIF4)
    box on
    %subplot(2,1,2)
end

print('-depsc',['tmp',int2str(k),'.eps'])

