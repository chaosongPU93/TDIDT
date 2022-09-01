%Sits at one spot (rotations; delays) and looks for the cross-correlation.
close all
clear all
t=cputime; tt=cputime;
format short e
%Rotations & offsets:
%timoffrot=[227 14 18 28 +058 -074  90 85 50]; %Quiet day
%timoffrot=[227 20 08 20 +028 +004 100 380 150];
%timoffrot=[246 00 61 00 +060 -085  85  65  45];  %Near beginning; %isolated; 100%
%timoffrot=[246 04 21 16 +062 -084  90  70  45];
%timoffrot=[246 04 21 56 +061 -084  90  60  50]; %Also quite isolated
%timoffrot=[246 04 22 20 +060 -085  85  60  40];
%timoffrot=[246 04 26 36 +059 -085  90  60  45];
%timoffrot=[246 04 25 16 +060 -085  90  60  45]; %GOOD FOR BROWN
%timoffrot=[246 04 26 52 +059 -085  90  60  45];
%timoffrot=[246 09 45 56 +061 -083  95  65  50];
%timoffrot=[246 11 08 52 +064 -083  90  60  50];
%timoffrot=[246 12 53 24 +062 -083  95  60  50];
%timoffrot=[246 09 47 40 +060 -084  90  60  50];
%timoffrot=[246 09 56 12 +062 -085  95  55  50]; %Also small
%timoffrot=[246 10 06 36 +061 -085  85  60  45]; %36415-36420, e.g.
%timoffrot=[246 11 11 16 +063 -084  90  60  50];
%timoffrot=[246 12 56 12 +062 -085  90  55  50];
%timoffrot=[246 17 17 48 +062 -087  80  60  40];
%timoffrot=[247 04 59 24 +033 -100  85 370  80];
%timoffrot=[247 10 15 24 +051 -078  85  60  55 ];
%timoffrot=[247 11 14 36 +056 -077  90  55  65];
%timoffrot=[247 14 18 28 +058 -074  90 85 50]; %(the original; catalog has 85  75  55)
%timoffrot=[247 14 42 28 +059 -073  85  65  65];  %Small signal!
%timoffrot=[250 05 09 24 +065 -068  80  75  45]; %Near start of a small migration.  V. messy.
%timoffrot=[250 05 23 56 +060 -072  85  65  50];
%timoffrot=[250 05 44 44 +056 -075  85  70  55]; %Possible EGFs here. FIRST LOOK!
     %timoffrot=[250 05 44 44 -104 +085  85  70  55]; %FAKE!
     %timoffrot=[250 05 44 44 +016 -035  85  70  55]; %FAKE!
%timoffrot=[250 09 13 40 +105 -031  85  50  80];
%timoffrot=[250 11 20 36 +105 -030  90  50  80]; %Big; good for looking for "missed" energy
%timoffrot=[250 11 25 00 +093 -035  85  70  90];
%timoffrot=[250 11 28 36 +090 -037  80  75  65];
%timoffrot=[250 11 30 44 +090 -036  75  60  70]; %early in a migration episode
%timoffrot=[250 11 34 20 +086 -032  80 105  65];
     %timoffrot=[250 11 34 20 -074 +128  80 105  65]; %FAKE!
%timoffrot=[250 11 34 20 +086 -032 +024  80 105  65  80];
%timoffrot=[250 11 47 00 +080 -042  75  70  65];
%timoffrot=[250 15 17 24 +104 -031  90  55  80]; %Big!
%timoffrot=[250 06 21 00 +100 -036  85  50  70]; %Also big, not quite as.jjj
%timoffrot=[253 01 50 20 +093 -024  85  75  85];
%timoffrot=[254 06 47 08 +086 +019  80 105  55]; %Dominated by last half of last panel
%timoffrot=[254 07 40 52 +086 +019  80 110  55]; %Consistent with lxtwmg %(whole 2nd panel, e.g.)
%timoffrot=[254 09 57 08 +086 +019  80 115  50]; %2nd panel quite impressive, again.
%timoffrot=[254 13 48 28 +086 +019  80 115  50];
%timoffrot=[254 13 50 04 +086 +019  80 115  50]; %BIG SPIKE, panel 2.  Coherence starts slightly earlier?
%timoffrot=[254 15 11 56 +086 +019  80 120  55]; %Messy; dominated by panel 2
%timoffrot=[254 20 29 00 +086 +019  85 115  55]; %Messy; mix of in-phase/out-of-phase end of panel 2
%timoffrot=[254 22 17 56 +086 +019  75 105  45];
%timoffrot=[254 05 54 12 +086 +018  80 115  50 ];
timoffrot=[254 04 05 00 +087 +018  85 110  45]; %Arrival of front; not much to see
timoffrot=[254 09 38 40 +085 +020 -002 -034 -079 -182 -144  85 110  45  85 340]; %particular front
%timoffrot=[254 09 38 40 +085 +020 -002 -034 -079 +020 +085  85 110  45  85 340];
timoffrot=[254 07 40 44 +085 +020 -002 -034 -079 -182 -144  85 110  45  85 340]; %SNB, SILBz spot
timoffrot=[255 05 19 40 +085 +020 -002 -034 -079 -182 -144  85 110  45  85 340]; %4STNs
timoffrot=[255 09 58 20 +085 +020 -002 -034 -079 -182 -144  85 110  45  85 340]; %4STNs
timoffrot=[255 17 22 04 +085 +020 -002 -034 -079 -182 -144  85 110  45  85 340]; %4STNs
timoffrot=[255 09 60 20 +085 +020 -002 -034 -079 -182 -144  85 110  45  85 340]; %4STNs
timoffrot=[255 02 14 44 +085 +020 -002 -034 -079 -182 -144  85 110  45  85 340]; %4STNs
timoffrot=[255 16 26 04 +085 +020 -002 -034 -079 -182 -144  85 110  45  85 340]; %4STNs
timoffrot=[256 01 04 04 +085 +020 -002 -034 -079 -182 -144  85 110  45  85 340]; %4STNs
timoffrot=[256 22 56 36 +085 +020 -002 -034 -079 -182 -144  85 110  45  85 340]; %4STNs
%timoffrot=[254 07 40 44 +085 +020 -002 -034 -079 +020 +085  85 110  45  85 340]; %SNB, SILBz spot
%timoffrot=[254 07 42 12 +086 +016  80 110  60];  %Quite amazing migration
%timoffrot=[254 07 41 00 +085 +020  80 115  55];  %Quite amazing 
%***********************
jday=timoffrot(1);
JDAY=int2str(jday)
timlook=3600*timoffrot(2)+60*timoffrot(3)+timoffrot(4);
rotPGC=pi*(timoffrot(12)-90)/180;
rotSSIB=pi*(timoffrot(13)-90)/180;
rotSILB=pi*(timoffrot(14)-90)/180;
rotKLNB=pi*(timoffrot(15)-90)/180;
rotVGZ=pi*(timoffrot(16)-90)/180;

SSIBsoff=timoffrot(5);
SILBsoff=timoffrot(6);
KLNBsoff=timoffrot(7);
VGZsoff=timoffrot(8);
SNBpsoff=timoffrot(9);
SILpsoff=timoffrot(10);
SSIpsoff=timoffrot(11);

SSIBtoff=SSIBsoff/40;
SILBtoff=SILBsoff/40;
KLNBtoff=KLNBsoff/40;
VGZtoff=VGZsoff/40;
SNBptoff=SNBpsoff/40;
SILptoff=SILpsoff/40;
SSIptoff=SSIpsoff/40;

%Read data:
direc='2005/SEPT/';
prename=[direc,'2005.',JDAY,'.00.00.00.0000.CN'];
PGCEdat=[prename,'.PGC..BHE.D.SAC'];
PGCNdat=[prename,'.PGC..BHN.D.SAC'];
SSIBEdat=[prename,'.SSIB..HHE.D.SAC'];
SSIBNdat=[prename,'.SSIB..HHN.D.SAC'];
SILBEdat=[prename,'.SILB..HHE.D.SAC'];
SILBNdat=[prename,'.SILB..HHN.D.SAC'];
KLNBEdat=[prename,'.KLNB..HHE.D.SAC'];
KLNBNdat=[prename,'.KLNB..HHN.D.SAC'];
VGZEdat=[prename,'.VGZ..BHE.D.SAC'];
VGZNdat=[prename,'.VGZ..BHN.D.SAC'];
SNBZdat=[prename,'.SNB..BHE.D.SAC'];
SILBZdat=[prename,'.SILB..HHZ.D.SAC'];
SSIBZdat=[prename,'.SSIB..HHZ.D.SAC'];

[PGCE,HdrDataPGC,tnuPGC,pobjPGC,timsPGC]=readsac(PGCEdat,0,'l');
[PGCN,~,~,~,~]=readsac(PGCNdat,0,'l');
[SSIBE,HdrDataSSIB,tnuSSIB,pobjSSIB,timsSSIB]=readsac(SSIBEdat,0,'l');
[SSIBN,~,~,~,~]=readsac(SSIBNdat,0,'l');
[SSIBZ,~,~,~,~]=readsac(SSIBZdat,0,'l');
[SILBE,HdrDataSILB,tnuSILB,pobjSILB,timsSILB]=readsac(SILBEdat,0,'l');
[SILBN,~,~,~,~]=readsac(SILBNdat,0,'l');
[SILBZ,~,~,~,~]=readsac(SILBZdat,0,'l');
[KLNBE,HdrDataKLNB,tnuKLNB,pobjKLNB,timsKLNB]=readsac(KLNBEdat,0,'l');
[KLNBN,~,~,~,~]=readsac(KLNBNdat,0,'l');
[VGZE,HdrDataVGZ,tnuVGZ,pobjVGZ,timsVGZ]=readsac(VGZEdat,0,'l');
[VGZN,~,~,~,~]=readsac(VGZNdat,0,'l');
[SNBZ,HdrDataSNB,tnuSNB,pobjSNB,timsSNB]=readsac(SNBZdat,0,'l');
%
rdsac=cputime-t
t=cputime;
    % Truncate for quicker calculation.  100sps stations are longer for
    % later interpolation
    timwin=1*60+1; %4*
    timstart=max(0, timlook-timwin); timend=min(86400, timlook+timwin);
    PGCE=PGCE(timstart*40+1:timend*40)-mean(PGCE(timstart*40+1:timend*40));
    PGCN=PGCN(timstart*40+1:timend*40)-mean(PGCN(timstart*40+1:timend*40)); 
    VGZE=VGZE(timstart*40+1:timend*40)-mean(VGZE(timstart*40+1:timend*40));
    VGZN=VGZN(timstart*40+1:timend*40)-mean(VGZN(timstart*40+1:timend*40)); 
    SNBZ=SNBZ(timstart*40+1:timend*40)-mean(SNBZ(timstart*40+1:timend*40)); 
    SSIBE=SSIBE(timstart*100:timend*100+1)-mean(SSIBE(timstart*100:timend*100+1));
    SSIBN=SSIBN(timstart*100:timend*100+1)-mean(SSIBN(timstart*100:timend*100+1));
    SSIBZ=SSIBZ(timstart*100:timend*100+1)-mean(SSIBZ(timstart*100:timend*100+1));
    SILBE=SILBE(timstart*100:timend*100+1)-mean(SILBE(timstart*100:timend*100+1));
    SILBN=SILBN(timstart*100:timend*100+1)-mean(SILBN(timstart*100:timend*100+1));
    SILBZ=SILBZ(timstart*100:timend*100+1)-mean(SILBZ(timstart*100:timend*100+1));
    KLNBE=KLNBE(timstart*100:timend*100+1)-mean(KLNBE(timstart*100:timend*100+1));
    KLNBN=KLNBN(timstart*100:timend*100+1)-mean(KLNBN(timstart*100:timend*100+1));
    timsPGC=timsPGC(timstart*40+1:timend*40);
    timsVGZ=timsVGZ(timstart*40+1:timend*40);
    timsSNB=timsSNB(timstart*40+1:timend*40);
    timsSSIB=timsSSIB(timstart*100:timend*100+1);
    timsSILB=timsSILB(timstart*100:timend*100+1);
    timsKLNB=timsKLNB(timstart*100:timend*100+1);
%shrten=cputime-t
%t=cputime;
%
tracelenPGC=length(PGCE);
tracelenVGZ=length(VGZE);
tracelenSNB=length(SNBZ);
tracelenSSIB=length(SSIBE);
tracelenSILB=length(SILBE);
tracelenKLNB=length(KLNBE);
    timPGCfirst=timsPGC(1);
    timSSIBfirst=timsSSIB(1); %0.01s earlier than PGC
    timPGClast=timsPGC(tracelenPGC);
    timSSIBlast=timsSSIB(tracelenSSIB);
    startdiff=timPGCfirst-timSSIBfirst
    enddiff=timSSIBlast-timPGClast

%cosine taper before filtering:
x=(0:pi/80:pi/2-pi/80)';
% %Seems to be necessary at the start of each day for PGCE:
%     PGCE(1:80)=0.;
%     PGCE(81:120)=sin(x).*PGCE(81:120); %Only at start of day!
%
PGCE(1:40)=sin(x).*PGCE(1:40);
PGCN(1:40)=sin(x).*PGCN(1:40);
VGZE(1:40)=sin(x).*VGZE(1:40);
VGZN(1:40)=sin(x).*VGZN(1:40);
SNBZ(1:40)=sin(x).*SNBZ(1:40);
x=flipud(x);
PGCE(tracelenPGC-39:tracelenPGC)=sin(x).*PGCE(tracelenPGC-39:tracelenPGC);
PGCN(tracelenPGC-39:tracelenPGC)=sin(x).*PGCN(tracelenPGC-39:tracelenPGC);
VGZE(tracelenVGZ-39:tracelenVGZ)=sin(x).*VGZE(tracelenVGZ-39:tracelenVGZ);
VGZN(tracelenVGZ-39:tracelenVGZ)=sin(x).*VGZN(tracelenVGZ-39:tracelenVGZ);
SNBZ(tracelenSNB-39:tracelenSNB)=sin(x).*SNBZ(tracelenSNB-39:tracelenSNB);
x=(0:pi/200:pi/2-pi/200)';
SSIBE(1:100)=sin(x).*SSIBE(1:100);
SSIBN(1:100)=sin(x).*SSIBN(1:100);
SSIBZ(1:100)=sin(x).*SSIBZ(1:100);
SILBE(1:100)=sin(x).*SILBE(1:100);
SILBN(1:100)=sin(x).*SILBN(1:100);
SILBZ(1:100)=sin(x).*SILBZ(1:100);
KLNBE(1:100)=sin(x).*KLNBE(1:100);
KLNBN(1:100)=sin(x).*KLNBN(1:100);
x=flipud(x);
SSIBE(tracelenSSIB-99:tracelenSSIB)=sin(x).*SSIBE(tracelenSSIB-99:tracelenSSIB);
SSIBN(tracelenSSIB-99:tracelenSSIB)=sin(x).*SSIBN(tracelenSSIB-99:tracelenSSIB);
SSIBZ(tracelenSSIB-99:tracelenSSIB)=sin(x).*SSIBZ(tracelenSSIB-99:tracelenSSIB);
SILBE(tracelenSILB-99:tracelenSILB)=sin(x).*SILBE(tracelenSILB-99:tracelenSILB);
SILBN(tracelenSILB-99:tracelenSILB)=sin(x).*SILBN(tracelenSILB-99:tracelenSILB);
SILBZ(tracelenSILB-99:tracelenSILB)=sin(x).*SILBZ(tracelenSILB-99:tracelenSILB);
KLNBE(tracelenKLNB-99:tracelenKLNB)=sin(x).*KLNBE(tracelenKLNB-99:tracelenKLNB);
KLNBN(tracelenKLNB-99:tracelenKLNB)=sin(x).*KLNBN(tracelenKLNB-99:tracelenKLNB);

% if(rem(tracelen,2)==0)
%     SeisData(tracelen)=[];
%     tims(tracelen)=[];
%     tracelen=tracelen-1;
% end
%Looks like filtering the 40-Hz data followed by interpolation introduces 
%less high-frequency stuff that does not appear in the original filtered
%40-Hz data:
%Filter data:
hi=6.;
lo=1.5;
% hi=8.;
% lo=2;
npo=2;
npa=1;
[PGCEf]=1.6e-4*bandpass(PGCE,40,lo,hi,npo,npa,'butter');
[PGCNf]=1.6e-4*bandpass(PGCN,40,lo,hi,npo,npa,'butter');
[VGZEf]=1.2e-3*bandpass(VGZE,40,lo,hi,npo,npa,'butter');
[VGZNf]=1.2e-3*bandpass(VGZN,40,lo,hi,npo,npa,'butter');
[SNBZf]=5*3.2e-4*bandpass(SNBZ,40,lo,hi,npo,npa,'butter');
[SSIBEf]=1.4*4.0e-3*bandpass(SSIBE,100,lo,hi,npo,npa,'butter');
[SSIBNf]=1.4*4.0e-3*bandpass(SSIBN,100,lo,hi,npo,npa,'butter');
[SSIBZf]=2*1.4*4.0e-3*bandpass(SSIBZ,100,lo,hi,npo,npa,'butter');
[SILBEf]=1.2*4.0e-3*bandpass(SILBE,100,lo,hi,npo,npa,'butter');
[SILBNf]=1.2*4.0e-3*bandpass(SILBN,100,lo,hi,npo,npa,'butter');
[SILBZf]=2*1.2*4.0e-3*bandpass(SILBZ,100,lo,hi,npo,npa,'butter');
[KLNBEf]=0.8*4.0e-3*bandpass(KLNBE,100,lo,hi,npo,npa,'butter');
[KLNBNf]=0.8*4.0e-3*bandpass(KLNBN,100,lo,hi,npo,npa,'butter');
%fltr=cputime-t
%t=cputime;
%Decimate the 100 sps data:
SSIBEfd = interp1(timsSSIB,SSIBEf,timsPGC,'linear')';
SSIBNfd = interp1(timsSSIB,SSIBNf,timsPGC,'linear')';
SSIBZfd = interp1(timsSSIB,SSIBZf,timsPGC,'linear')';
SILBEfd = interp1(timsSILB,SILBEf,timsPGC,'linear')';
SILBNfd = interp1(timsSILB,SILBNf,timsPGC,'linear')';
SILBZfd = interp1(timsSILB,SILBZf,timsPGC,'linear')';
KLNBEfd = interp1(timsKLNB,KLNBEf,timsPGC,'linear')';
KLNBNfd = interp1(timsKLNB,KLNBNf,timsPGC,'linear')';
% SSIBEfd = resample(SSIBEf,2,5);
% SSIBNfd = resample(SSIBNf,2,5);
% SILBEfd = resample(SILBEf,2,5);
% SILBNfd = resample(SILBNf,2,5);
%dmate=cputime-t
%t=cputime;

PGC=PGCEf+1i*PGCNf;
VGZ=VGZEf+1i*VGZNf;
SNBz=SNBZf;
SSIB=SSIBEfd+1i*SSIBNfd;
SSIBz=SSIBZfd;
SILB=SILBEfd+1i*SILBNfd;
SILBz=SILBZfd;
KLNB=KLNBEfd+1i*KLNBNfd;
PGCrot=PGC*exp(1i*rotPGC);
VGZrot=VGZ*exp(1i*rotVGZ);
SSIBrot=SSIB*exp(1i*rotSSIB);
SILBrot=SILB*exp(1i*rotSILB);
KLNBrot=KLNB*exp(1i*rotKLNB);
    %rtate=cputime-t
    %t=cputime;

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

if KLNBtoff > 0
    KLNBrot(1:tracelenPGC-KLNBsoff)=KLNBrot(KLNBsoff+1:tracelenPGC);
    KLNBrot(tracelenPGC-KLNBsoff+1:tracelenPGC)=0;
else
    KLNBrot(-KLNBsoff+1:tracelenPGC)=KLNBrot(1:tracelenPGC+KLNBsoff);
    KLNBrot(1:-KLNBsoff)=0;
end

if VGZtoff > 0
    VGZrot(1:tracelenPGC-VGZsoff)=VGZrot(VGZsoff+1:tracelenPGC);
    VGZrot(tracelenPGC-VGZsoff+1:tracelenPGC)=0;
else
    VGZrot(-VGZsoff+1:tracelenPGC)=VGZrot(1:tracelenPGC+VGZsoff);
    VGZrot(1:-VGZsoff)=0;
end

if SNBptoff > 0
    SNBz(1:tracelenPGC-SNBpsoff)=SNBz(SNBpsoff+1:tracelenPGC);
    SNBz(tracelenPGC-SNBpsoff+1:tracelenPGC)=0;
else
    SNBz(-SNBpsoff+1:tracelenPGC)=SNBz(1:tracelenPGC+SNBpsoff);
    SNBz(1:-SNBpsoff)=0;
end

if SSIptoff > 0
    SSIBz(1:tracelenPGC-SSIpsoff)=SSIBz(SSIpsoff+1:tracelenPGC);
    SSIBz(tracelenPGC-SSIpsoff+1:tracelenPGC)=0;
    SSIsrot(1:tracelenPGC-SSIpsoff)=SSIBrot(SSIpsoff+1:tracelenPGC);
    SSIsrot(tracelenPGC-SSIpsoff+1:tracelenPGC)=0;
    SSIsrot=SSIsrot'; %Not sure why needed
else
    SSIBz(-SSIpsoff+1:tracelenPGC)=SSIBz(1:tracelenPGC+SSIpsoff);
    SSIBz(1:-SSIpsoff)=0;
    SSIsrot(-SSIpsoff+1:tracelenPGC)=SSIBrot(1:tracelenPGC+SSIpsoff);
    SSIsrot(1:-SSIpsoff)=0;
    SSIsrot=SSIsrot'; %Not sure why needed
end

if SILptoff > 0
    SILBz(1:tracelenPGC-SILpsoff)=SILBz(SILpsoff+1:tracelenPGC);
    SILBz(tracelenPGC-SILpsoff+1:tracelenPGC)=0;
    SILsrot(1:tracelenPGC-SILpsoff)=SILBrot(SILpsoff+1:tracelenPGC);
    SILsrot(tracelenPGC-SILpsoff+1:tracelenPGC)=0;
    SILsrot=SILsrot'; %not sure why needed
else
    SILBz(-SILpsoff+1:tracelenPGC)=SILBz(1:tracelenPGC+SILpsoff);
    SILBz(1:-SILpsoff)=0;
    SILsrot(-SILpsoff+1:tracelenPGC)=SILBrot(1:tracelenPGC+SILpsoff);
    SILsrot(1:-SILpsoff)=0;
    SILsrot=SILsrot'; %not sure why needed
end
    %tracesht=cputime-t
    %t=cputime;

plotstart=1; %(timlook-timstart)*40
plotend=(2*timwin)*40; %(timlook-timstart)*40
ltmax=max([real(PGCrot(plotstart:plotend)); ...
          real(SSIBrot(plotstart:plotend)); ...
          real(SILBrot(plotstart:plotend))]);
ltmin=min([real(PGCrot(plotstart:plotend)); ...
          real(SSIBrot(plotstart:plotend)); ...
          real(SILBrot(plotstart:plotend))]);

realPGC=real(PGCrot);
realSSIB=real(SSIBrot);
realSILB=real(SILBrot);
realSSIs=real(SSIsrot); %with P offsets
realSILs=real(SILsrot); %with P offsets
PGC2=dot(realPGC,realPGC);
SSIB2=dot(realSSIB,realSSIB);
SILB2=dot(realSILB,realSILB);
%SSIs2=dot(realSSIs,realSSIs);
%SILs2=dot(realSILs,realSILs);
% 
%Rotations
%Rotating the vector by -alpha is the same as rotating the axes by alpha.
SILBnext=realSILs+1i*SILBz;
SSIBnext=realSSIs+1i*SSIBz;
angs=0:-pi/36:-2*pi+pi/36;
nang=length(angs);
for n=1:nang
    SILBnrot(:,n)=SILBnext*exp(1i*angs(n));
    realSILBz(:,n)=real(SILBnrot(:,n));
    SILBz2=dot(realSILBz(:,n),realSILBz(:,n));
    SILBSILB(n)=dot(realSILBz(:,n),realSILB)/sqrt(SILBz2*SILB2);
    SILBSSIB(n)=dot(realSILBz(:,n),realSSIB)/sqrt(SILBz2*SSIB2);
    SILBPGC(n)=dot(realSILBz(:,n),realPGC)/sqrt(SILBz2*PGC2);

    SSIBnrot(:,n)=SSIBnext*exp(1i*angs(n));
    realSSIBz(:,n)=real(SSIBnrot(:,n));
    SSIBz2=dot(realSSIBz(:,n),realSSIBz(:,n));
    SSIBSILB(n)=dot(realSSIBz(:,n),realSILB)/sqrt(SSIBz2*SILB2);
    SSIBSSIB(n)=dot(realSSIBz(:,n),realSSIB)/sqrt(SSIBz2*SSIB2);
    SSIBPGC(n)=dot(realSSIBz(:,n),realPGC)/sqrt(SSIBz2*PGC2);
end
[maxxc,imaxSILBSILB]=max(SILBSILB)
[maxxc,imaxSILBSSIB]=max(SILBSSIB)
[maxxc,imaxSILBPGC]=max(SILBPGC)
[maxxc,imaxSSIBSILB]=max(SSIBSILB)
[maxxc,imaxSSIBSSIB]=max(SSIBSSIB)
[maxxc,imaxSSIBPGC]=max(SSIBPGC)
%[maxkeep(keeptot),ikeep]=max(keeptmpx); %keep track of the kept maxima

figure %Superposed traces, 150s window:
subplot(4,1,1); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
plot(timsPGC,realSILBz(:,imaxSILBSILB),'c','LineWidth',2);
% plot(timsPGC,realSSIBz(:,imaxSSIBSILB),'c','LineWidth',2);
axis([timlook-timwin timlook-timwin/2 ltmin ltmax]);
subplot(4,1,2); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
plot(timsPGC,realSILBz(:,imaxSILBSILB),'c','LineWidth',2);
%plot(timsPGC,realSSIBz(:,imaxSSIBSILB),'c','LineWidth',2);
axis([timlook-timwin/2 timlook ltmin ltmax]);
subplot(4,1,3); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
plot(timsPGC,realSILBz(:,imaxSILBSILB),'c','LineWidth',2);
%plot(timsPGC,realSSIBz(:,imaxSSIBSILB),'c','LineWidth',2);
axis([timlook timlook+timwin/2 ltmin ltmax]);
subplot(4,1,4); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
plot(timsPGC,realSILBz(:,imaxSILBSILB),'c','LineWidth',2);
%plot(timsPGC,realSSIBz(:,imaxSSIBSILB),'c','LineWidth',2);
axis([timlook+timwin/2 timlook+timwin ltmin ltmax]);


% mshift=0;
% lenx=tracelenPGC-2*mshift;
% 
% templen=3*40;
% tempoff=40/2;
% %tempoff=10*40;
% ntemp=floor((lenx-templen)/tempoff);
% %ntemp=floor((lenx-templen-10*40)/tempoff);
% timstemp=zeros(ntemp,1);
% PGCtx=zeros(tracelenPGC-templen+1,ntemp); %no +1 so you can subtract ___(j-1)
% PGCtn=zeros(tracelenPGC-templen+1,ntemp);
% SSIBtx=zeros(tracelenPGC-templen+1,ntemp); %no +1 so you can subtract ___(j-1)
% SSIBtn=zeros(tracelenPGC-templen+1,ntemp);
% SILBtx=zeros(tracelenPGC-templen+1,ntemp); %no +1 so you can subtract ___(j-1)
% SILBtn=zeros(tracelenPGC-templen+1,ntemp);
% alltx=zeros(tracelenPGC-templen+1,ntemp);
% alltn=zeros(tracelenPGC-templen+1,ntemp);
% tottn=zeros(tracelenPGC-templen+1,ntemp);
% n05=0; j05=0; f05=0;
% for n=1:ntemp;
%     n
%     istart=(n-1)*tempoff+mshift+1; %+mshift to coincide with ...?
%     %istart=(n-1)*tempoff+mshift+1+5*40; %+mshift to coincide with ...?
%     iend=istart+templen-1;
%     PGCtempl=realPGC(istart:iend);
%     SILBtempl=realSILB(istart:iend);
%     SSIBtempl=realSSIB(istart:iend);
%     PGCt2=dot(PGCtempl,PGCtempl);
%     SSIBt2=dot(SSIBtempl,SSIBtempl);
%     SILBt2=dot(SILBtempl,SILBtempl);
%         %     for j=2:istart-40 %so you can subtract ___(j-1) %!!!JUST LOOK AFTER TEMPLATE, FOR NOW!!!
%         %         PGCtx(j,n)=dot(PGCtempl,realPGC(j:j+templen-1));
%         %         PGCtarg2=PGC2(j+templen-1)-PGC2(j-1);
%         %         PGCtn(j,n)=PGCtx(j,n)/sqrt(PGCt2*PGCtarg2);
%         %         SSIBtx(j,n)=dot(SSIBtempl,realSSIB(j:j+templen-1));
%         %         SSIBtarg2=SSIB2(j+templen-1)-SSIB2(j-1);
%         %         SSIBtn(j,n)=SSIBtx(j,n)/sqrt(SSIBt2*SSIBtarg2);
%         %         SILBtx(j,n)=dot(SILBtempl,realSILB(j:j+templen-1));
%         %         SILBtarg2=SILB2(j+templen-1)-SILB2(j-1);
%         %         SILBtn(j,n)=SILBtx(j,n)/sqrt(SILBt2*SILBtarg2);
%         %         alltx(j,n)=PGCtx(j,n)+SSIBtx(j,n)+SILBtx(j,n);
%         %         alltn(j,n)=alltx(j,n)/sqrt((PGCt2+SSIBt2+SILBt2)*(PGCtarg2+SSIBtarg2+SILBtarg2));
%         %         tottn(j,n)=(PGCtn(j,n)+SSIBtn(j,n)+SILBtn(j,n))/3.0;
%         %         if alltn(j,n) > 0.45
%         %         %if tottn(j,n) > 0.52
%         %             alltn(j,n)
%         %             f05=f05+1;
%         %             n05=n05+1; j05=j05+1;
%         %             tempfindt(f05,:)=timsPGC(istart:iend);
%         %             targfindt(f05,:)=timsPGC(j:j+templen-1);
%         %             tempfindPGC(:,f05)=PGCtempl;
%         %             tempfindSSIB(:,f05)=SSIBtempl;
%         %             tempfindSILB(:,f05)=SILBtempl;
%         %             targfindPGC(:,f05)=realPGC(j:j+templen-1);
%         %             targfindSSIB(:,f05)=realSSIB(j:j+templen-1);
%         %             targfindSILB(:,f05)=realSILB(j:j+templen-1);
%         %             find05(f05,:)=[n,j]; %Why is this n,j and alltn j,n?
%         %         end
%         %     end
%     for j=istart+40:tracelenPGC-templen+1 
%     %for j=istart+1:tracelenPGC-templen+1 
%         PGCtx(j,n)=dot(PGCtempl,realPGC(j:j+templen-1));
%         PGCtarg2=PGC2(j+templen-1)-PGC2(j-1);
%         PGCtn(j,n)=PGCtx(j,n)/sqrt(PGCt2*PGCtarg2);
%         SSIBtx(j,n)=dot(SSIBtempl,realSSIB(j:j+templen-1));
%         SSIBtarg2=SSIB2(j+templen-1)-SSIB2(j-1);
%         SSIBtn(j,n)=SSIBtx(j,n)/sqrt(SSIBt2*SSIBtarg2);
%         SILBtx(j,n)=dot(SILBtempl,realSILB(j:j+templen-1));
%         SILBtarg2=SILB2(j+templen-1)-SILB2(j-1);
%         SILBtn(j,n)=SILBtx(j,n)/sqrt(SILBt2*SILBtarg2);
%         alltx(j,n)=PGCtx(j,n)+SSIBtx(j,n)+SILBtx(j,n);
%         alltn(j,n)=alltx(j,n)/sqrt((PGCt2+SSIBt2+SILBt2)*(PGCtarg2+SSIBtarg2+SILBtarg2));
%         tottn(j,n)=(PGCtn(j,n)+SSIBtn(j,n)+SILBtn(j,n))/3.0;
%         if alltn(j,n) > 0.4
%         %if tottn(j,n) > 0.42
%             alltn(j,n)
%             f05=f05+1;
%             n05=n05+1; j05=j05+1;
%             tempfindt(f05,:)=timsPGC(istart:iend);
%             targfindt(f05,:)=timsPGC(j:j+templen-1);
%             tempfindPGC(:,f05)=PGCtempl;
%             tempfindSSIB(:,f05)=SSIBtempl;
%             tempfindSILB(:,f05)=SILBtempl;
%             targfindPGC(:,f05)=realPGC(j:j+templen-1);
%             targfindSSIB(:,f05)=realSSIB(j:j+templen-1);
%             targfindSILB(:,f05)=realSILB(j:j+templen-1);
%             find05(f05,:)=[n,j];
%         end
%     end
% end
% %lump things
% find05'
% out=zeros(f05,1); 
% keeptot=0;
% for i=1:f05
%     if out(i) < 0.5 %if not previously kicked out (out=0 instead of 1)
%         clear keeptmp
%         clear keeptmpx
%         keeptot=keeptot+1; %how many will be kept in the end
%         inset=1; %how many in the set to compare; at least this one
%         keeptmp(inset)=i; 
%         keeptmpx(inset)=alltn(find05(i,2),find05(i,1)); %alltn(n,j) of the i_th member of f05
%         keeptmpPGCtn(inset)=PGCtn(find05(i,2),find05(i,1));
%         keeptmpSSIBtn(inset)=SSIBtn(find05(i,2),find05(i,1));
%         keeptmpSILBtn(inset)=SILBtn(find05(i,2),find05(i,1));
%         keeptmpnj(inset,:)=[find05(i,1), find05(i,2)];
%         for k=i+1:f05
%             if isequal(find05(i,1),find05(k,1)) && find05(k,2)-find05(i,2)<=80 %Same template, targets closer than 80 samples
%                 inset=inset+1;
%                 keeptmp(inset)=k; out(k)=1;
%                 keeptmpx(inset)=alltn(find05(k,2),find05(k,1));
%                 keeptmpPGCtn(inset)=PGCtn(find05(k,2),find05(k,1));
%                 keeptmpSSIBtn(inset)=SSIBtn(find05(k,2),find05(k,1));
%                 keeptmpSILBtn(inset)=SILBtn(find05(k,2),find05(k,1));
%                 keeptmpnj(inset,:)=[find05(k,1), find05(k,2)];
%             elseif find05(k,1)-find05(i,1)<=85.0/tempoff && abs(find05(k,2)-find05(i,2))<=85 %Templates within a few (4); targets?
%                 inset=inset+1;
%                 keeptmp(inset)=k; out(k)=1;
%                 keeptmpx(inset)=alltn(find05(k,2),find05(k,1));
%                 keeptmpPGCtn(inset)=PGCtn(find05(k,2),find05(k,1));
%                 keeptmpSSIBtn(inset)=SSIBtn(find05(k,2),find05(k,1));
%                 keeptmpSILBtn(inset)=SILBtn(find05(k,2),find05(k,1));
%                 keeptmpnj(inset,:)=[find05(k,1), find05(k,2)];
%            end
%         end
%         keeptmp
%         keeptmpx
%         [maxkeep(keeptot),ikeep]=max(keeptmpx); %keep track of the kept maxima
%         ikeep %ikeep varies only from 1 to inset, the number in the set being compared
%         keeptmp(ikeep) %keeptmp(ikeep) varies between 1 and f05
%         tempkeept(keeptot,:)=tempfindt(keeptmp(ikeep),:);
%         targkeept(keeptot,:)=targfindt(keeptmp(ikeep),:);
%         tempkeepPGC(:,keeptot)=tempfindPGC(:,keeptmp(ikeep));
%         tempkeepSSIB(:,keeptot)=tempfindSSIB(:,keeptmp(ikeep));
%         tempkeepSILB(:,keeptot)=tempfindSILB(:,keeptmp(ikeep));
%         targkeepPGC(:,keeptot)=targfindPGC(:,keeptmp(ikeep));
%         targkeepSSIB(:,keeptot)=targfindSSIB(:,keeptmp(ikeep));
%         targkeepSILB(:,keeptot)=targfindSILB(:,keeptmp(ikeep));
%         keepPGCtn(keeptot)=keeptmpPGCtn(ikeep);
%         keepSSIBtn(keeptot)=keeptmpSSIBtn(ikeep);
%         keepSILBtn(keeptot)=keeptmpSILBtn(ikeep);
%         keepnj(keeptot,:)=keeptmpnj(ikeep,:);
%     end
% end
% keeptot
% [maxkeepsort,mksi] = sort(maxkeep,'descend');
% ncheck=floor((lenx-templen)/templen); %looks for the max in a string of pt-wise cc m'ments
% histmat=zeros(ncheck,ntemp);
% for n=1:ncheck
%     istart=(n-1)*templen+1;
%     iend=n*templen;
%     histmat(n,:)=max(alltn(istart:iend,:));
% end
% nlump=floor(ntemp/4);
% hist2mat=zeros(ncheck,nlump);
% for n=1:nlump
%     istart=(n-1)*4+1;
%     iend=n*4;
%     hist2mat(:,n)=max(histmat(:,istart:iend),[],2);
% end
% [m,n]=size(hist2mat);
% hist2=reshape(hist2mat,1,m*n);
% [look,x]=hist(hist2,50);
% look2=[x;look];
% fid = fopen('look.txt','w');
% fprintf(fid,'%8.4f  %8i\n',look2);
% fclose(fid);
%         
%         % orient landscape
%         % print('-dpdf',[JDAY,'_',int2str(timlook),'_4.pdf'])tempkeept
% slope=(ltmax-ltmin)/160;
% for ifig=1:floor(keeptot/8)
%     figure
%     for isub=1:8
%         icount=(ifig-1)*8+isub;
%         crossoff=(2*ltmax+(0:159));
%         for ij=1:icount-1  %Get rid of "duplicates" (same small offsets for template and target)
%             if abs((tempkeept(mksi(icount),1)-tempkeept(mksi(ij),1)) ...
%                   -(targkeept(mksi(icount),1)-targkeept(mksi(ij),1))) <=0.1 && ...
%                abs((tempkeept(mksi(icount),1)-tempkeept(mksi(ij),1))) <= 4.0
% %                     targkeepPGC(:,mksi(icount))=0.;
% %                     targkeepSSIB(:,mksi(icount))=0.;
% %                     targkeepSILB(:,mksi(icount))=0.;
%                     crossoff=(ltmin+slope*(0:159));
%             end
%         end
%         subplot(4,6,(isub-1)*3+1,'align'); 
%         hold on
%         plot(targkeepPGC(:,mksi(icount)),'r');
%         plot(tempkeepPGC(:,mksi(icount)),'g');
%         plot(crossoff,'k--')
%         text(10, 0.8*ltmax, int2str(icount));
%         text(100, 0.8*ltmax, num2str(maxkeep(mksi(icount))),'fontsize',7);
%         text(100, 0.8*ltmin, num2str(keepPGCtn(mksi(icount))),'fontsize',7);
%         text(10, 0.8*ltmin, 'PGC','fontsize',7);
%         text(60, 0.6*ltmax, int2str(tempkeept(mksi(icount),floor(templen/2))),'fontsize',7);
%         text(110, 0.6*ltmax, int2str(targkeept(mksi(icount),floor(templen/2))),'fontsize',7);
%         axis([0 160 ltmin ltmax]);
%         set(gca,'fontsize',8);
%         box on
%         subplot(4,6,(isub-1)*3+2,'align'); 
%         hold on
%         plot(targkeepSSIB(:,mksi(icount)),'b');
%         plot(tempkeepSSIB(:,mksi(icount)),'g');
%         plot(crossoff,'k--')
%         text(10, 0.8*ltmax, int2str(icount));
%         text(100, 0.8*ltmax, num2str(maxkeep(mksi(icount))),'fontsize',7);
%         text(100, 0.8*ltmin, num2str(keepSSIBtn(mksi(icount))),'fontsize',7);
%         text(10, 0.8*ltmin, 'SSIB','fontsize',7);
%         text(60, 0.6*ltmax, int2str(tempkeept(mksi(icount),floor(templen/2))),'fontsize',7);
%         text(110, 0.6*ltmax, int2str(targkeept(mksi(icount),floor(templen/2))),'fontsize',7);
%         %text(100, 0.6*ltmax, int2str(keepnj(mksi(icount),:)));
%         axis([0 160 ltmin ltmax]);
%         set(gca,'fontsize',8);
%         box on
%         subplot(4,6,(isub-1)*3+3,'align'); 
%         hold on
%         plot(targkeepSILB(:,mksi(icount)),'k');
%         plot(tempkeepSILB(:,mksi(icount)),'g');
%         plot(crossoff,'k--')
%         text(10, 0.8*ltmax, int2str(icount));
%         text(100, 0.8*ltmax, num2str(maxkeep(mksi(icount))),'fontsize',7);
%         text(100, 0.8*ltmin, num2str(keepSILBtn(mksi(icount))),'fontsize',7);
%         text(10, 0.8*ltmin, 'SILB','fontsize',7);
%         text(60, 0.6*ltmax, int2str(tempkeept(mksi(icount),floor(templen/2))),'fontsize',7);
%         text(110, 0.6*ltmax, int2str(targkeept(mksi(icount),floor(templen/2))),'fontsize',7);
%         %text(100, 0.6*ltmax, int2str(keepnj(mksi(icount),:)));
%         axis([0 160 ltmin ltmax]);
%         set(gca,'fontsize',8);
%         box on
%     end
%     orient landscape
%     print('-dpdf',[JDAY,'_',int2str(timlook),'_Brown',int2str(ifig),'.pdf'])
% end

%cputime-t;
tot=cputime-tt
