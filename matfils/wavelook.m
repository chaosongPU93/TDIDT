%Sits at one spot (rotations; delays) and looks for the cross-correlation.
close all
clear all
t=cputime; tt=cputime;
format short e
%Rotations & offsets:
timoffrot=[192 03 25 16 +059 -085 +055 -174  000 -318  000  90  60  45  140 30]; %good for KLN?
timoffrot=[121 23 50 17 +048 -071 +051 -174  000 -071 +048  95  90  50  140 30]; %Impressive signal
%timoffrot=[121 01 23 36 +058 -085 +055 -174  000 -085 +058  95  90  50  140 30]; %
timoffrot=[254 09 38 40 +085 +020 -002 -034 -079 -182 -144  85 110  45  85 340]; %particular front
timoffrot=[254 09 55 16 +085 +020 -002 -034 -079 -182 -144  85 110  45  85 340]; %particular front
%timoffrot=[254 09 55 16 +085 +020 -002 -034 -079 +020 +085  85 110  45  85 340]; %particular front
timoffrot=[254 13 47 34 +088 +020 -002 -034 -079 -182 -144  85 110  45  85 340]; %particular front
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
    timwin=3*60+1; %4*
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
hi=8.;
lo=2;
% hi=10.;
% lo=3;
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
else
    SSIBz(-SSIpsoff+1:tracelenPGC)=SSIBz(1:tracelenPGC+SSIpsoff);
    SSIBz(1:-SSIpsoff)=0;
end

if SILptoff > 0
    SILBz(1:tracelenPGC-SILpsoff)=SILBz(SILpsoff+1:tracelenPGC);
    SILBz(tracelenPGC-SILpsoff+1:tracelenPGC)=0;
else
    SILBz(-SILpsoff+1:tracelenPGC)=SILBz(1:tracelenPGC+SILpsoff);
    SILBz(1:-SILpsoff)=0;
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

timstart=timlook-8; timend=timlook+8;
figure %Superposed traces, 150s window:
subplot(6,1,1); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([timstart timend ltmin ltmax]);
subplot(6,1,2); 
plot(timsPGC,real(KLNBrot),'c');
axis([timstart timend ltmin ltmax]);
subplot(6,1,3); 
plot(timsPGC,real(VGZrot),'c');
xlim([timstart timend]);
subplot(6,1,4); 
plot(timsPGC,SNBz,'c');
xlim([timstart timend]);
subplot(6,1,5); 
plot(timsPGC,SILBz,'c');
xlim([timstart timend]);
subplot(6,1,6); 
plot(timsPGC,SSIBz,'c');
xlim([timstart timend]);

