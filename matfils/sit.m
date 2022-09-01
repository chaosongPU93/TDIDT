%Sits at one spot (rotations; delays) and looks for the cross-correlation.
t=cputime; tt=cputime;
close all
format long e
%Rotations & offsets:
%timoffrot=[227 06 23 48 +017 +067  90 220 165]; %A single 0.2 detection
%timoffrot=[227 14 18 28 +058 -074  90 85 50]; %Quiet day
%timoffrot=[227 20 08 20 +028 +004 100 380 150];
%timoffrot=[246 00 61 00 +060 -085  85  65  45];  %Near beginning;
%isolated; 100%
%timoffrot=[246 04 21 16 +062 -084  90  70  45];
%timoffrot=[246 04 19 16 +062 -084  90  70  45]; %Added
%timoffrot=[246 04 21 56 +061 -084  90  60  50]; %Also quite isolated
%timoffrot=[246 04 22 20 +060 -085  85  60  40];
 %timoffrot=[246 04 25 16 +060 -085  90  60  45];
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
 timoffrot=[247 14 18 28 +058 -074  90 85 50]; %(the original; catalog has 85  75  55)
%timoffrot=[247 14 42 28 +059 -073  85  65  65];  %Small signal!
%timoffrot=[250 05 09 24 +065 -068  80  75  45]; %Near start of a small migration.  V. messy.
%timoffrot=[250 05 23 56 +060 -072  85  65  50];
 %timoffrot=[250 05 44 44 +056 -075  85  70  55]; %Possible EGFs here
%timoffrot=[250 09 13 40 +105 -031  85  50  80];
%timoffrot=[250 11 20 36 +105 -030  90  50  80]; %Big; good for looking for "missed" energy
%timoffrot=[250 11 25 00 +093 -035  85  70  90];
%timoffrot=[250 11 28 36 +090 -037  80  75  65];
%timoffrot=[250 11 30 44 +090 -036  75  60  70]; %early in a migration episode
%timoffrot=[250 11 34 20 +086 -032  80 105  65];
%timoffrot=[250 11 34 20 +086 -032 +024  80 105  65  80];
%timoffrot=[250 11 47 00 +080 -042  75  70  65];
%timoffrot=[250 15 17 24 +104 -031  90  55  80]; %Big!
%timoffrot=[250 06 21 00 +100 -036  85  50  70]; %Also big, not quite as.
%timoffrot=[253 01 50 20 +093 -024  85  75  85];
%timoffrot=[254 06 47 08 +086 +019  80 105  55]; %Dominated by last half of last panel
%timoffrot=[254 07 40 52 +086 +019  80 110  55]; %Consistent with lxtwmg (whole 2nd panel, e.g.)
%timoffrot=[254 09 57 08 +086 +019  80 115  50];
%timoffrot=[254 13 48 28 +086 +019  80 115  50];
%timoffrot=[254 13 50 04 +086 +019  80 115  50]; %BIG SPIKE, panel 1
%timoffrot=[254 15 11 56 +086 +019  80 120  55]; %Messy; dominated by panel 2
%timoffrot=[254 20 29 00 +086 +019  85 115  55]; %Messy; mix of in-phase/out-of-phase end of panel 2
%timoffrot=[254 22 17 56 +086 +019  75 105  45];
jday=timoffrot(1);
JDAY=int2str(jday)
timlook=3600*timoffrot(2)+60*timoffrot(3)+timoffrot(4)
rotPGC=pi*(timoffrot(7)-90)/180;
rotSSIB=pi*(timoffrot(8)-90)/180;
rotSILB=pi*(timoffrot(9)-90)/180;
SSIBsoff=timoffrot(5);
SILBsoff=timoffrot(6);
SSIBtoff=SSIBsoff/40;
SILBtoff=SILBsoff/40;
%
% SSIBsoff=51;
% SILBsoff=-78;
% SSIBsoff=93;
% SILBsoff=-24;

%Read Armbruster's detections:
ArmCat=load('/data2/arubin/ARMB/BC2005/list.2005.pgsssiMAJ');
detects=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
m=0;
Detects=0;
for n=1:length(detects)
     %if ArmCat(n,5)==SSIBsoff && ArmCat(n,6)==SILBsoff
     if isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff]) ||...
        isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff]) ||...
        isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff]) ||...
        isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff-1]) ||...
        isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff+1]) ||...
        isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff-1]) ||...
        isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff+1]) ||...
        isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff-1]) ||...
        isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff+1]) 
         m=m+1;
         Detects(m)=detects(n);
     end
end
ArmCat=load('/data2/arubin/ARMB/BC2005/list.2005.pgsssiMAJ_new');
detects2_8=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
m=0;
Detects2_8=0;
for n=1:length(detects2_8)
     if isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff]) ||...
        isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff]) ||...
        isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff]) ||...
        isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff-1]) ||...
        isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff+1]) ||...
        isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff-1]) ||...
        isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff+1]) ||...
        isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff-1]) ||...
        isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff+1]) 
         m=m+1;
         Detects2_8(m)=detects2_8(n);
     end
end

%Read data:
prename=['2005.',JDAY,'.00.00.00.0000.CN'];
PGCEdat=[prename,'.PGC..BHE.D.SAC'];
PGCNdat=[prename,'.PGC..BHN.D.SAC'];
PGCZdat=[prename,'.PGC..BHZ.D.SAC'];
SSIBEdat=[prename,'.SSIB..HHE.D.SAC'];
SSIBNdat=[prename,'.SSIB..HHN.D.SAC'];
SILBEdat=[prename,'.SILB..HHE.D.SAC'];
SILBNdat=[prename,'.SILB..HHN.D.SAC'];

[PGCE,HdrDataPGC,tnuPGC,pobjPGC,timsPGC]=readsac(PGCEdat,0,'l');
[PGCN,~,~,~,~]=readsac(PGCNdat,0,'l');
[PGCZ,~,~,~,~]=readsac(PGCZdat,0,'l');
[SSIBE,HdrDataSSIB,tnuSSIB,pobjSSIB,timsSSIB]=readsac(SSIBEdat,0,'l');
[SSIBN,~,~,~,~]=readsac(SSIBNdat,0,'l');
[SILBE,HdrDataSILB,tnuSILB,pobjSILB,timsSILB]=readsac(SILBEdat,0,'l');
[SILBN,~,~,~,~]=readsac(SILBNdat,0,'l');
%
rdsac=cputime-t
t=cputime;
    % Truncate for quicker calculation.  100sps stations are longer for
    % later interpolation
    timwin=12*60;
    timstart=max(0, timlook-timwin); timend=min(86400, timlook+timwin);
    PGCE=PGCE(timstart*40+1:timend*40);
    PGCN=PGCN(timstart*40+1:timend*40);
    PGCZ=PGCZ(timstart*40+1:timend*40);
    SSIBE=SSIBE(timstart*100:timend*100+1);
    SSIBN=SSIBN(timstart*100:timend*100+1);
    SILBE=SILBE(timstart*100:timend*100+1);
    SILBN=SILBN(timstart*100:timend*100+1);
    timsPGC=timsPGC(timstart*40+1:timend*40);
    timsSSIB=timsSSIB(timstart*100:timend*100+1);
    timsSILB=timsSILB(timstart*100:timend*100+1);
shrten=cputime-t
t=cputime;
%
tracelenPGC=length(PGCE);
tracelenSSIB=length(SSIBE);
tracelenSILB=length(SILBE);
    timPGCfirst=timsPGC(1);
    timSSIBfirst=timsSSIB(1);
    timSILBfirst=timsSILB(1);
    timPGClast=timsPGC(tracelenPGC);
    timSSIBlast=timsSSIB(tracelenSSIB);
    timSILBlast=timsSILB(tracelenSILB);
    startdiff=timPGCfirst-timSSIBfirst
    enddiff=timSSIBlast-timPGClast
%pause

%cosine taper before filtering:
x=(0:pi/80:pi/2-pi/80)';
% %Seems to be necessary at the start of each day for PGCE:
%     PGCE(1:80)=0.;
%     PGCE(81:120)=sin(x).*PGCE(81:120); %Only at start of day!
%
PGCE(1:40)=sin(x).*PGCE(1:40);
PGCN(1:40)=sin(x).*PGCN(1:40);
PGCZ(1:40)=sin(x).*PGCZ(1:40);
x=flipud(x);
PGCE(tracelenPGC-39:tracelenPGC)=sin(x).*PGCE(tracelenPGC-39:tracelenPGC);
PGCN(tracelenPGC-39:tracelenPGC)=sin(x).*PGCN(tracelenPGC-39:tracelenPGC);
PGCZ(tracelenPGC-39:tracelenPGC)=sin(x).*PGCZ(tracelenPGC-39:tracelenPGC);
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
% %Check up on tapering:
% figure
% plot(PGCE);
% xlim([0 140]);
% figure
% plot(PGCN);
% xlim([tracelenPGC-79 tracelenPGC]);
% figure
% plot(SSIBE);
% xlim([tracelenSSIB-199 tracelenSSIB]);
% figure
% plot(SSIBN);
% xlim([0 200]);
% figure
% plot(SILBE);
% xlim([tracelenSILB-199 tracelenSILB]);
% figure
% plot(SILBN);
% xlim([0 200]);

% if(rem(tracelen,2)==0)
%     SeisData(tracelen)=[];
%     tims(tracelen)=[];
%     tracelen=tracelen-1;
% end
%Looks like filtering the 40-Hz data followed by interpolation introduces 
%less high-frequency stuff that does not appear in the original filtered
%40-Hz data:
%PGCE=spline(tims,SeisData,timsx);
%[PGCE]=bandpass(PGCE,100,1.0,5.0,2,1,'butter');

%Filter data:
hi=6.;
lo=1.5;
% hi=1.;
% lo=0.1;
npo=2;
npa=1;
[PGCEf]=1.6e-4*bandpass(PGCE,40,lo,hi,npo,npa,'butter');
[PGCNf]=1.6e-4*bandpass(PGCN,40,lo,hi,npo,npa,'butter');
[PGCZf]=1.6e-4*bandpass(PGCZ,40,lo,hi,npo,npa,'butter');
[SSIBEf]=1.4*4.0e-3*bandpass(SSIBE,100,lo,hi,npo,npa,'butter');
[SSIBNf]=1.4*4.0e-3*bandpass(SSIBN,100,lo,hi,npo,npa,'butter');
[SILBEf]=1.2*4.0e-3*bandpass(SILBE,100,lo,hi,npo,npa,'butter');
[SILBNf]=1.2*4.0e-3*bandpass(SILBN,100,lo,hi,npo,npa,'butter');
fltr=cputime-t
t=cputime;
%Decimate the 100 sps data:
SSIBEfd = interp1(timsSSIB,SSIBEf,timsPGC,'linear')';
SSIBNfd = interp1(timsSSIB,SSIBNf,timsPGC,'linear')';
SILBEfd = interp1(timsSILB,SILBEf,timsPGC,'linear')';
SILBNfd = interp1(timsSILB,SILBNf,timsPGC,'linear')';
dmate=cputime-t
t=cputime;
% SSIBEfd(1:6)
% SSIBNfd(1:6)
% SILBEfd(1:6)
% SILBNfd(1:6)

% %Check up on interpolated data:
% figure
% plot(timsPGC,PGCEf);
% xlim([51545.5 51583]);
% figure
% plot(timsPGC,PGCNf);
% xlim([51545.5 51583]);
% figure
% plot(timsSSIB,SSIBEf);
% hold on
% plot(timsPGC,SSIBEfd,'r');
% xlim([51545.5 51583]);
% figure
% plot(timsSSIB,SSIBNf);
% hold on
% plot(timsPGC,SSIBNfd,'r');
% xlim([51545.5 51583]);
% figure
% plot(timsSILB,SILBEf);
% hold on
% plot(timsPGC,SILBEfd,'r');
% xlim([51545.5 51583]);
% figure
% plot(timsSILB,SILBNf);
% hold on
% plot(timsPGC,SILBNfd,'r');
% xlim([51545.5 51583]);

PGC=PGCEf+1i*PGCNf;
SSIB=SSIBEfd+1i*SSIBNfd;
SILB=SILBEfd+1i*SILBNfd;
PGCrot=PGC*exp(1i*rotPGC);
SSIBrot=SSIB*exp(1i*rotSSIB);
SILBrot=SILB*exp(1i*rotSILB);
    rtate=cputime-t
    t=cputime;

timsSSIB=timsPGC-SSIBtoff;
timsSILB=timsPGC-SILBtoff;
    tshft=cputime-t
    t=cputime;

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
    tracesht=cputime-t
    t=cputime;

plotstart=(timlook-timstart-75)*40;
plotend=(timlook-timstart+75)*40;
ltmax=max([real(PGCrot(plotstart:plotend)); ...
          real(SSIBrot(plotstart:plotend)); ...
          real(SILBrot(plotstart:plotend))])
ltmin=min([real(PGCrot(plotstart:plotend)); ...
          real(SSIBrot(plotstart:plotend)); ...
          real(SILBrot(plotstart:plotend))])
%Pointwise x-correlation:
PGSSx=real(PGCrot).*real(SSIBrot);
PGSIx=real(PGCrot).*real(SILBrot);
SSSIx=real(SSIBrot).*real(SILBrot);
    ptwise=cputime-t
    t=cputime;
lxmax=max([PGSSx(plotstart:plotend); PGSIx(plotstart:plotend); SSSIx(plotstart:plotend)]);
lxmin=min([PGSSx(plotstart:plotend); PGSIx(plotstart:plotend); SSSIx(plotstart:plotend)]);
figure %Pointwise x-correlation, 150s window:
subplot(4,1,1); 
hold on
plot(timsPGC,real(PGCrot),'r');
%plot(timsPGC,5*PGCZf,'g');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([timlook-75 timlook-37.5 ltmin ltmax]);
%xlim([timlook-75 timlook-37.5]);
subplot(4,1,2); 
hold on
plot(timsPGC,PGSIx,'b');
plot(timsPGC,PGSSx,'r');
plot(timsPGC,SSSIx,'k'); 
axis([timlook-75 timlook-37.5 lxmin lxmax]);
subplot(4,1,3); 
hold on
plot(timsPGC,real(PGCrot),'r');
%plot(timsPGC,5*PGCZf,'g');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([timlook-37.5 timlook ltmin ltmax]);
%xlim([timlook-37.5 timlook]);
subplot(4,1,4); 
hold on
plot(timsPGC,PGSIx,'b');
plot(timsPGC,PGSSx,'r');
plot(timsPGC,SSSIx,'k');
axis([timlook-37.5 timlook lxmin lxmax]);

figure
subplot(4,1,1); 
hold on
plot(timsPGC,real(PGCrot),'r');
%plot(timsPGC,5*PGCZf,'g');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
% text(40, 140.,'M2.4 earthquake, 45 km depth, 1-5 Hz','FontSize',12)
% text(40, 90.,'Rotated & shifted to maximize coherence','FontSize',12)
text(timlook+1, 0.8*ltmax,'PGC','Color',[1 0 0],'FontSize',7)
text(timlook+1, 0.5*ltmax,'SSIB','Color',[0 0 1],'FontSize',7)
text(timlook+1, -0.6*ltmax,'SILB','Color',[0 0 0],'FontSize',7)
XLABEL('sec')
YLABEL('Velocity counts')
axis([timlook timlook+37.5 ltmin ltmax]);
box on
%xlim([timlook timlook+37.5]);
subplot(4,1,2); 
hold on
plot(timsPGC,PGSIx,'b');
plot(timsPGC,PGSSx,'r');
plot(timsPGC,SSSIx,'k');
text(timlook+1, 0.8*lxmax,'PGC dot SSIB','Color',[1 0 0],'FontSize',6)
text(timlook+1, 0.6*lxmax,'PGC dot SILB','Color',[0 0 1],'FontSize',6)
text(timlook+1, 0.4*lxmax,'SSIB dot SILB','Color',[0 0 0],'FontSize',6)
XLABEL('sec')
axis([timlook timlook+37.5 lxmin lxmax]);
box on
subplot(4,1,3); 
hold on
plot(timsPGC,real(PGCrot),'r');
%plot(timsPGC,5*PGCZf,'g');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
text(timlook+38.5, 0.8*ltmax,'PGC','Color',[1 0 0],'FontSize',7)
text(timlook+38.5, 0.5*ltmax,'SSIB','Color',[0 0 1],'FontSize',7)
text(timlook+38.5, -0.6*ltmax,'SILB','Color',[0 0 0],'FontSize',7)
XLABEL('sec')
YLABEL('Velocity counts')
axis([timlook+37.5 timlook+75 ltmin ltmax]);
box on
%xlim([timlook+37.5 timlook+75]);
subplot(4,1,4); 
hold on
plot(timsPGC,PGSIx,'b');
plot(timsPGC,PGSSx,'r');
plot(timsPGC,SSSIx,'k');
text(timlook+38.5, 0.8*lxmax,'PGC dot SSIB','Color',[1 0 0],'FontSize',6)
text(timlook+38.5, 0.6*lxmax,'PGC dot SILB','Color',[0 0 1],'FontSize',6)
text(timlook+38.5, 0.4*lxmax,'SSIB dot SILB','Color',[0 0 0],'FontSize',6)
XLABEL('sec')
axis([timlook+37.5 timlook+75 lxmin lxmax]);
box on
%orient landscape
print('-depsc',[JDAY,'_',int2str(timlook),'_sit2.eps'])

PGSSxcum=cumsum(PGSSx);
PGSIxcum=cumsum(PGSIx);
SSSIxcum=cumsum(SSSIx);
    cmsm=cputime-t
    t=cputime;
figure %Running sum x-correlation pairs, with detections, all day
% plot(timsPGC,PGSIxcum-PGSIxcum(2057320),'r');
% hold on
% plot(timsPGC,PGSSxcum-PGSSxcum(2057320),'k');
% plot(timsPGC,SSSIxcum-SSSIxcum(2057320),'b');
plot(timsPGC,PGSIxcum,'b');
hold on
plot(timsPGC,PGSSxcum,'r');
plot(timsPGC,SSSIxcum,'k');
hrf = plotreflinesr(gca,detects,'x','k'); 
hrf = plotreflinesr(gca,Detects,'x','r'); 
hrf = plotreflinesr(gca,detects2_8,'x','y'); 
hrf = plotreflinesr(gca,Detects2_8,'x','g'); 
xlim([timlook-75 timlook+75]);
%xlim([timsPGC(1) timsPGC(tracelenPGC)]);
%xlim([0 86400]);
%xlim ([51433 51583]);
    plotcmsm=cputime-t
    t=cputime;

%Autocorrelation of stations:
PGCauto=real(PGCrot).*real(PGCrot);
PGC2=cumsum(PGCauto);
SSIBauto=real(SSIBrot).*real(SSIBrot);
SSIB2=cumsum(SSIBauto);
SILBauto=real(SILBrot).*real(SILBrot);
SILB2=cumsum(SILBauto);
PGCautoi=imag(PGCrot).*imag(PGCrot);
PGC2i=cumsum(PGCautoi);
SSIBautoi=imag(SSIBrot).*imag(SSIBrot);
SSIB2i=cumsum(SSIBautoi);
SILBautoi=imag(SILBrot).*imag(SILBrot);
SILB2i=cumsum(SILBautoi);
    autocorrandsumboth=cputime-t
    t=cputime;

%Get 3-sec sums:
tmp=reshape(PGSSx,120,tracelenPGC/120);
PGSSx3=sum(tmp);
tmp=reshape(PGSIx,120,tracelenPGC/120);
PGSIx3=sum(tmp);
tmp=reshape(SSSIx,120,tracelenPGC/120);
SSSIx3=sum(tmp);
tmp=reshape(PGCauto,120,tracelenPGC/120);
PGCauto3=sum(tmp);
tmp=reshape(SSIBauto,120,tracelenPGC/120);
SSIBauto3=sum(tmp);
tmp=reshape(SILBauto,120,tracelenPGC/120);
SILBauto3=sum(tmp);
timsPGC3=timsPGC(1)+1.5:3:timsPGC(tracelenPGC)-1.;
    threesecs=cputime-t
    t=cputime;

    figure %3-sec window ptwise x-correlation pairs and individual station strengths, +/-1000 sec 
    lmax=max([PGSIx3 PGSSx3 SSSIx3]);
    lmin=min([PGSIx3 PGSSx3 SSSIx3]);
    subplot(2,1,1);
    plot(timsPGC3,PGSIx3,'b');
    hold on
    plot(timsPGC3,PGSSx3,'r');
    plot(timsPGC3,SSSIx3,'k');
        pltthreesecs=cputime-t
        t=cputime;
    hrf = plotreflinesr(gca,detects,'x','k'); 
    hrf = plotreflinesr(gca,Detects,'x','r'); 
    hrf = plotreflinesr(gca,detects2_8,'x','y'); 
    hrf = plotreflinesr(gca,Detects2_8,'x','g'); 
        pltreflines=cputime-t
        t=cputime;
    text(timstart+60, 0.8*lmax,'(smoothed) PGC dot SSIB (3 sec)','Color',[1 0 0],'FontSize',9)
    text(timstart+60, 0.63*lmax,'(smoothed) PGC dot SILB','Color',[0 0 1],'FontSize',9)
    text(timstart+60, 0.46*lmax,'(smoothed) SSIB dot SILB','Color',[0 0 0],'FontSize',9)
    XLABEL('sec')
    axis([timstart timend lmin lmax])
    box on
    subplot(2,1,2);
    lmax=max([PGCauto3 SSIBauto3 SILBauto3]);
    lmin=min([PGCauto3 SSIBauto3 SILBauto3]);
    plot(timsPGC3,PGCauto3,'r');
    hold on
    plot(timsPGC3,SSIBauto3,'b');
    plot(timsPGC3,SILBauto3,'k');
        pltthreesecs=cputime-t
        t=cputime;
    hrf = plotreflinesr(gca,detects,'x','k'); 
    hrf = plotreflinesr(gca,Detects,'x','r'); 
    hrf = plotreflinesr(gca,detects2_8,'x','y'); 
    hrf = plotreflinesr(gca,Detects2_8,'x','g'); 
    text(timstart+60, 0.8*lmax,'(smoothed) PGC^2 (3 sec)','Color',[1 0 0],'FontSize',10)
    text(timstart+60, 0.65*lmax,'(smoothed) SSIB^2','Color',[0 0 1],'FontSize',10)
    text(timstart+60, 0.5*lmax,'(smoothed) SILB^2','Color',[0 0 0],'FontSize',10)
    XLABEL('sec')
    axis([timstart timend lmin lmax]);
    %orient landscape
    print('-depsc',[JDAY,'_',int2str(timlook),'_sit4.eps'])

pltreflines=cputime-t
        t=cputime;
    %xlim([0 86400]);
    %axis([timlook-1000 timlook+1000 -0.5 2]);
    %axis([10000 40000 -0.1 1]);
    %axis([34500 37500 -0.1 1]);

    figure %Running sum signal strength at each station, all day
    plot(timsPGC,PGC2,'r');
    hold on
    plot(timsPGC,PGC2i,'r--');
    plot(timsPGC,SSIB2,'b');
    plot(timsPGC,SSIB2i,'b--');
    plot(timsPGC,SILB2,'k');
    plot(timsPGC,SILB2i,'k--');
    hrf = plotreflinesr(gca,detects,'x','k');
    hrf = plotreflinesr(gca,Detects,'x','r');
    hrf = plotreflinesr(gca,detects2_8,'x','y');
    hrf = plotreflinesr(gca,Detects2_8,'x','g');
    xlim([timsPGC(1) timsPGC(tracelenPGC)]);

% %%%ORTHOGONAL COMPONENTS%%%
% figure %Individual orth. station traces, 150s window
% %lmax=max([real(PGCrot) real(SSIBrot) real(SILBrot)])
% subplot(4,1,1); 
% hold on
% plot(timsPGC,imag(PGCrot),'r');
% plot(timsSSIB,imag(SSIBrot),'b');
% plot(timsSILB,imag(SILBrot),'k');
% axis([timlook-75 timlook-37.5 lmin lmax]);
% %axis([timlook-75 timlook -0.2 0.2]);
% subplot(4,1,2); 
% hold on
% plot(timsPGC,imag(PGCrot),'r');
% plot(timsSSIB,imag(SSIBrot),'b');
% plot(timsSILB,imag(SILBrot),'k');
% axis([timlook-37.5 timlook lmin lmax]);
% %axis([timlook timlook+75 -0.2 0.2]);
% subplot(4,1,3); 
% hold on
% plot(timsPGC,imag(PGCrot),'r');
% plot(timsSSIB,imag(SSIBrot),'b');
% plot(timsSILB,imag(SILBrot),'k');
% axis([timlook timlook+37.5 lmin lmax]);
% %axis([timlook timlook+75 -0.2 0.2]);
% subplot(4,1,4); 
% hold on
% plot(timsPGC,imag(PGCrot),'r');
% plot(timsSSIB,imag(SSIBrot),'b');
% plot(timsSILB,imag(SILBrot),'k');
% axis([timlook+37.5 timlook+75 lmin lmax]);
% %axis([timlook timlook+75 -0.2 0.2]);
% pltrace=cputime-t
%     t=cputime;
% 
% %Pointwise x-correlation:
% PGSSxi=imag(PGCrot).*imag(SSIBrot);
% PGSIxi=imag(PGCrot).*imag(SILBrot);
% SSSIxi=imag(SSIBrot).*imag(SILBrot);
%     ptwise=cputime-t
%     t=cputime;
% figure %Pointwise x-correlation, 150s window:
% subplot(4,1,1); 
% hold on
% plot(timsPGC,PGSIxi,'b');
% plot(timsPGC,PGSSxi,'r');
% plot(timsPGC,SSSIxi,'k'); 
% axis([timlook-75 timlook-37.5 lmin lmax]);
% subplot(4,1,2); 
% hold on
% plot(timsPGC,PGSIxi,'b');
% plot(timsPGC,PGSSxi,'r');
% plot(timsPGC,SSSIxi,'k');
% axis([timlook-37.5 timlook lmin lmax]);
% subplot(4,1,3); 
% hold on
% plot(timsPGC,PGSIxi,'b');
% plot(timsPGC,PGSSxi,'r');
% plot(timsPGC,SSSIxi,'k');
% axis([timlook timlook+37.5 lmin lmax]);
% subplot(4,1,4); 
% hold on
% plot(timsPGC,PGSIxi,'b');
% plot(timsPGC,PGSSxi,'r');
% plot(timsPGC,SSSIxi,'k');
% axis([timlook+37.5 timlook+75 lmin lmax]);
%     pltptwise=cputime-t
%     t=cputime;
% 
% % PGSSxi=imag(PGCrot).*imag(SSIBrot);
% % PGSIxi=imag(PGCrot).*imag(SILBrot);
% % SSSIxi=imag(SSIBrot).*imag(SILBrot);
% PGSSxicum=cumsum(PGSSxi);
% PGSIxicum=cumsum(PGSIxi);
% SSSIxicum=cumsum(SSSIxi);
%     orthxcorrandsum=cputime-t
%     t=cputime;
% figure %Running sum orthogonal x-correlation pairs, all day
% plot(timsPGC,PGSIxicum,'b');
% hold on
% plot(timsPGC,PGSSxicum,'r');
% plot(timsPGC,SSSIxicum,'k');
% xlim([timsPGC(1) timsPGC(tracelenPGC)]);
%     plotcmsmorth=cputime-t
%     t=cputime;

cputime-t
cputime-tt
