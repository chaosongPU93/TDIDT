%Sits at one spot (rotations; delays) and looks for the cross-correlation.
t=cputime; tt=cputime;
close all
format short e
%Rotations & offsets:
%timoffrot=[227 14 18 28 +058 -074  90 85 50]; %Quiet day
%timoffrot=[227 20 08 20 +028 +004 100 380 150];
%timoffrot=[246 00 61 00 +060 -085  85  65  45];  %Near beginning; %isolated; 100%
%timoffrot=[246 04 21 16 +062 -084  90  70  45];
%timoffrot=[246 04 21 56 +061 -084  90  60  50]; %Also quite isolated
%timoffrot=[246 04 22 20 +060 -085  85  60  40];
timoffrot=[246 04 25 16 +060 -085  90  60  45];
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
%timoffrot=[254 07 40 52 +086 +019  80 110  55]; %Consistent with lxtwmg %(whole 2nd panel, e.g.)
%timoffrot=[254 09 57 08 +086 +019  80 115  50]; %2nd panel quite impressive, again.
%timoffrot=[254 13 48 28 +086 +019  80 115  50];
%timoffrot=[254 13 50 04 +086 +019  80 115  50]; %BIG SPIKE, panel 2.  Coherence starts slightly earlier?
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

% %Read Armbruster's detections:
% ArmCat=load('/data2/arubin/ARMB/BC2005/list.2005.pgsssiMAJ');
% detects=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
% m=0;
% Detects=0;
% for n=1:length(detects)
%      %if ArmCat(n,5)==SSIBsoff && ArmCat(n,6)==SILBsoff
%      if isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff+1]) 
%          m=m+1;
%          Detects(m)=detects(n);
%      end
% end
% ArmCat=load('/data2/arubin/ARMB/BC2005/list.2005.pgsssiMAJ_new');
% detects2_8=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
% m=0;
% Detects2_8=0;
% for n=1:length(detects2_8)
%      if isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff+1]) 
%          m=m+1;
%          Detects2_8(m)=detects2_8(n);
%      end
% end

%Read data:
prename=['2005.',JDAY,'.00.00.00.0000.CN'];
PGCEdat=[prename,'.PGC..BHE.D.SAC'];
PGCNdat=[prename,'.PGC..BHN.D.SAC'];
%PGCZdat=[prename,'.PGC..BHZ.D.SAC'];
SSIBEdat=[prename,'.SSIB..HHE.D.SAC'];
SSIBNdat=[prename,'.SSIB..HHN.D.SAC'];
SILBEdat=[prename,'.SILB..HHE.D.SAC'];
SILBNdat=[prename,'.SILB..HHN.D.SAC'];

[PGCE,HdrDataPGC,tnuPGC,pobjPGC,timsPGC]=readsac(PGCEdat,0,'l');
[PGCN,~,~,~,~]=readsac(PGCNdat,0,'l');
%[PGCZ,~,~,~,~]=readsac(PGCZdat,0,'l');
[SSIBE,HdrDataSSIB,tnuSSIB,pobjSSIB,timsSSIB]=readsac(SSIBEdat,0,'l');
[SSIBN,~,~,~,~]=readsac(SSIBNdat,0,'l');
[SILBE,HdrDataSILB,tnuSILB,pobjSILB,timsSILB]=readsac(SILBEdat,0,'l');
[SILBN,~,~,~,~]=readsac(SILBNdat,0,'l');
% timsPGC(1)
% timsSSIB(1)
% timsSILB(1)
%
rdsac=cputime-t
t=cputime;
    % Truncate for quicker calculation.  100sps stations are longer for
    % later interpolation
    timwin=4*60+1;
    timstart=max(0, timlook-timwin); timend=min(86400, timlook+timwin);
    PGCE=PGCE(timstart*40+1:timend*40)-mean(PGCE(timstart*40+1:timend*40));
    PGCN=PGCN(timstart*40+1:timend*40)-mean(PGCN(timstart*40+1:timend*40)); 
    %PGCZ=PGCZ(timstart*40+1:timend*40);
    SSIBE=SSIBE(timstart*100:timend*100+1)-mean(SSIBE(timstart*100:timend*100+1));
    SSIBN=SSIBN(timstart*100:timend*100+1)-mean(SSIBN(timstart*100:timend*100+1));
    SILBE=SILBE(timstart*100:timend*100+1)-mean(SILBE(timstart*100:timend*100+1));
    SILBN=SILBN(timstart*100:timend*100+1)-mean(SILBN(timstart*100:timend*100+1));
    timsPGC=timsPGC(timstart*40+1:timend*40);
    timsSSIB=timsSSIB(timstart*100:timend*100+1);
    timsSILB=timsSILB(timstart*100:timend*100+1);
%shrten=cputime-t
%t=cputime;
%
tracelenPGC=length(PGCE);
tracelenSSIB=length(SSIBE);
tracelenSILB=length(SILBE);
    timPGCfirst=timsPGC(1);
    timSSIBfirst=timsSSIB(1); %0.01s earlier than PGC
    timSILBfirst=timsSILB(1); %0.01s earlier than PGC
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
%PGCZ(1:40)=sin(x).*PGCZ(1:40);
x=flipud(x);
PGCE(tracelenPGC-39:tracelenPGC)=sin(x).*PGCE(tracelenPGC-39:tracelenPGC);
PGCN(tracelenPGC-39:tracelenPGC)=sin(x).*PGCN(tracelenPGC-39:tracelenPGC);
%PGCZ(tracelenPGC-39:tracelenPGC)=sin(x).*PGCZ(tracelenPGC-39:tracelenPGC);
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
npo=2;
npa=1;
[PGCEf]=1.6e-4*bandpass(PGCE,40,lo,hi,npo,npa,'butter');
[PGCNf]=1.6e-4*bandpass(PGCN,40,lo,hi,npo,npa,'butter');
%[PGCZf]=1.6e-4*bandpass(PGCZ,40,lo,hi,npo,npa,'butter');
[SSIBEf]=4.0e-3*bandpass(SSIBE,100,lo,hi,npo,npa,'butter');
[SSIBNf]=4.0e-3*bandpass(SSIBN,100,lo,hi,npo,npa,'butter');
[SILBEf]=4.0e-3*bandpass(SILBE,100,lo,hi,npo,npa,'butter');
[SILBNf]=4.0e-3*bandpass(SILBN,100,lo,hi,npo,npa,'butter');
%fltr=cputime-t
%t=cputime;
%Decimate the 100 sps data:
SSIBEfd = interp1(timsSSIB,SSIBEf,timsPGC,'linear')';
SSIBNfd = interp1(timsSSIB,SSIBNf,timsPGC,'linear')';
SILBEfd = interp1(timsSILB,SILBEf,timsPGC,'linear')';
SILBNfd = interp1(timsSILB,SILBNf,timsPGC,'linear')';
%dmate=cputime-t
%t=cputime;
% SSIBEfd(1:6)
% SSIBNfd(1:6)
% SILBEfd(1:6)
% SILBNfd(1:6)

PGC=PGCEf+1i*PGCNf;
SSIB=SSIBEfd+1i*SSIBNfd;
SILB=SILBEfd+1i*SILBNfd;
PGCrot=PGC*exp(1i*rotPGC);
SSIBrot=SSIB*exp(1i*rotSSIB);
SILBrot=SILB*exp(1i*rotSILB);
    %rtate=cputime-t
    %t=cputime;

timsSSIB=timsPGC-SSIBtoff; %No longer used, I think
timsSILB=timsPGC-SILBtoff; %No longer used, I think
    %tshft=cputime-t
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
    %tracesht=cputime-t
    %t=cputime;

realPGC=real(PGCrot);
realSSIB=real(SSIBrot);
realSILB=real(SILBrot);
aimagPGC=imag(PGCrot); %A quick/dirty way to try the same on the orthogonal components
aimagSSIB=imag(SSIBrot);
aimagSILB=imag(SILBrot);

plotstart=1; %(timlook-timstart)*40
plotend=(2*timwin)*40; %(timlook-timstart)*40
ltmax=max([real(PGCrot(plotstart:plotend)); ...
          real(SSIBrot(plotstart:plotend)); ...
          real(SILBrot(plotstart:plotend))]);
ltmin=min([real(PGCrot(plotstart:plotend)); ...
          real(SSIBrot(plotstart:plotend)); ...
          real(SILBrot(plotstart:plotend))]);

figure %Superposed traces, 150s window:
subplot(4,1,1); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([timlook-timwin timlook-timwin/2 ltmin ltmax]);
%xlim([timlook-75 timlook-37.5]);
subplot(4,1,2); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([timlook-timwin/2 timlook ltmin ltmax]);
%xlim([timlook-37.5 timlook]);
subplot(4,1,3); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([timlook timlook+timwin/2 ltmin ltmax]);
%xlim([timlook timlook+37.5]);
subplot(4,1,4); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([timlook+timwin/2 timlook+timwin ltmin ltmax]);
%xlim([timlook+37.5 timlook+75]);

%Autocorrelation of stations:
PGCauto=realPGC.*realPGC;
PGC2=cumsum(PGCauto);
SSIBauto=realSSIB.*realSSIB;
SSIB2=cumsum(SSIBauto);
SILBauto=realSILB.*realSILB;
SILB2=cumsum(SILBauto);
%
PGCiauto=aimagPGC.*aimagPGC;
PGCi2=cumsum(PGCiauto);
SSIBiauto=aimagSSIB.*aimagSSIB;
SSIBi2=cumsum(SSIBiauto);
SILBiauto=aimagSILB.*aimagSILB;
SILBi2=cumsum(SILBiauto);
%

mshift=11;
lenx=tracelenPGC-2*mshift;
PGSSx=zeros(lenx, 2*mshift+1);
PGSIx=zeros(lenx, 2*mshift+1);
SISSx=zeros(lenx, 2*mshift+1);
PGSSix=zeros(lenx, 2*mshift+1);
PGSIix=zeros(lenx, 2*mshift+1);
SISSix=zeros(lenx, 2*mshift+1);
% First index is pointwise multiplication of traces; second is shifting offset.
for n=-mshift:mshift;
    PGSSx(:,n+mshift+1)=realPGC(1+mshift:tracelenPGC-mshift).* ...
        realSSIB(1+mshift-n:tracelenPGC-mshift-n);
    PGSIx(:,n+mshift+1)=realPGC(1+mshift:tracelenPGC-mshift).* ...
        realSILB(1+mshift-n:tracelenPGC-mshift-n);
    SISSx(:,n+mshift+1)=realSILB(1+mshift:tracelenPGC-mshift).* ...
        realSSIB(1+mshift-n:tracelenPGC-mshift-n);
    PGSSix(:,n+mshift+1)=aimagPGC(1+mshift:tracelenPGC-mshift).* ...
        aimagSSIB(1+mshift-n:tracelenPGC-mshift-n);
    PGSIix(:,n+mshift+1)=aimagPGC(1+mshift:tracelenPGC-mshift).* ...
        aimagSILB(1+mshift-n:tracelenPGC-mshift-n);
    SISSix(:,n+mshift+1)=aimagSILB(1+mshift:tracelenPGC-mshift).* ...
        aimagSSIB(1+mshift-n:tracelenPGC-mshift-n);
end
cumsumPGSS=cumsum(PGSSx);
cumsumPGSI=cumsum(PGSIx);
cumsumSISS=cumsum(SISSx);
cumsumPGSSi=cumsum(PGSSix);
cumsumPGSIi=cumsum(PGSIix);
cumsumSISSi=cumsum(SISSix);

timbig=3*60;
winbig=2*timbig*40;
igstart=floor(tracelenPGC/2-winbig/2)+1;
winlen=4*40;
winoff=40/2;
nwin=floor((winbig-winlen)/winoff);
timswin=zeros(nwin,1);
sumsPGSS=zeros(nwin,2*mshift+1);
sumsPGSI=zeros(nwin,2*mshift+1);
sumsSISS=zeros(nwin,2*mshift+1);
sumsPGC2=zeros(nwin,2*mshift+1);
sumsSSIB2=zeros(nwin,2*mshift+1);
sumsSILB2=zeros(nwin,2*mshift+1);
sumsSILB2b=zeros(nwin,2*mshift+1);
sumsPGSSi=zeros(nwin,2*mshift+1);
sumsPGSIi=zeros(nwin,2*mshift+1);
sumsSISSi=zeros(nwin,2*mshift+1);
sumsPGCi2=zeros(nwin,2*mshift+1);
sumsSSIBi2=zeros(nwin,2*mshift+1);
sumsSILBi2=zeros(nwin,2*mshift+1);
sumsSILBi2b=zeros(nwin,2*mshift+1);
for n=1:nwin;
    istart=igstart+(n-1)*winoff;
    iend=istart+winlen;
    timswin(n)=timsPGC(istart+winlen/2); %Gets rid of "mshifts" by referencing to the global igstart.
    %First index is window #; second is shifts over +-mshift.
    sumsPGSS(n,:)=cumsumPGSS(iend,:)-cumsumPGSS(istart-1,:); %Summed x-correlations
    sumsPGSI(n,:)=cumsumPGSI(iend,:)-cumsumPGSI(istart-1,:);
    sumsSISS(n,:)=cumsumSISS(iend,:)-cumsumSISS(istart-1,:);
    sumsPGC2(n,:)=PGC2(iend+mshift)-PGC2(istart+mshift-1);  %PGC doesn't shift, so do this here.  PGC2 is cumsummed. Yes, +mshift.
    sumsSILB2b(n,:)=SILB2(iend+mshift)-SILB2(istart+mshift-1); %Similar, for the SILB-SSIB connection.
    for m=-mshift:mshift;
        sumsSSIB2(n,m+mshift+1)=SSIB2(iend+mshift-m)-SSIB2(istart+mshift-1-m); %+m??? (yes).
        sumsSILB2(n,m+mshift+1)=SILB2(iend+mshift-m)-SILB2(istart+mshift-1-m);
    end

    sumsPGSSi(n,:)=cumsumPGSSi(iend,:)-cumsumPGSSi(istart-1,:); %Summed x-correlations
    sumsPGSIi(n,:)=cumsumPGSIi(iend,:)-cumsumPGSIi(istart-1,:);
    sumsSISSi(n,:)=cumsumSISSi(iend,:)-cumsumSISSi(istart-1,:);
    sumsPGCi2(n,:)=PGCi2(iend+mshift)-PGCi2(istart+mshift-1);  %PGC doesn't shift, so do this here.  PGC2 is cumsummed. Yes, +mshift.
    sumsSILBi2b(n,:)=SILBi2(iend+mshift)-SILBi2(istart+mshift-1); %Similar, for the SILB-SSIB connection.
    for m=-mshift:mshift;
        sumsSSIBi2(n,m+mshift+1)=SSIBi2(iend+mshift-m)-SSIBi2(istart+mshift-1-m); %+m??? (yes).
        sumsSILBi2(n,m+mshift+1)=SILBi2(iend+mshift-m)-SILBi2(istart+mshift-1-m);
    end
end
denomPGSSn=realsqrt(sumsPGC2.*sumsSSIB2);
denomPGSIn=realsqrt(sumsPGC2.*sumsSILB2);
denomSISSn=realsqrt(sumsSILB2b.*sumsSSIB2);
denomPGSSin=realsqrt(sumsPGCi2.*sumsSSIBi2);
denomPGSIin=realsqrt(sumsPGCi2.*sumsSILBi2);
denomSISSin=realsqrt(sumsSILBi2b.*sumsSSIBi2);
sumsPGSSn=sumsPGSS./denomPGSSn;
sumsPGSIn=sumsPGSI./denomPGSIn;
sumsSISSn=sumsSISS./denomSISSn;
sumsPGSSin=sumsPGSSi./denomPGSSin;
sumsPGSIin=sumsPGSIi./denomPGSIin;
sumsSISSin=sumsSISSi./denomSISSin;
%sumsPGSSn'
[xcmaxPGSSn,imaxPGSS]=max(sumsPGSSn,[],2);
[xcmaxPGSIn,imaxPGSI]=max(sumsPGSIn,[],2);
[xcmaxSISSn,imaxSISS]=max(sumsSISSn,[],2);
[xcmaxPGSSin,imaxPGSSi]=max(sumsPGSSin,[],2);
[xcmaxPGSIin,imaxPGSIi]=max(sumsPGSIin,[],2);
[xcmaxSISSin,imaxSISSi]=max(sumsSISSin,[],2);
%Parabolic fit:
%%% ADD by Chao, 2019/02/21, parabol has been changed to output only 3 para
%%% function [xmaxSTAS,ymaxSTAS,a]=parabol(nwin,mshift,sumsSTAS,imaxSTAS)
[xmaxPGSSn,ymaxPGSSn,xloPGSSn,xhiPGSSn]=parabol(nwin,mshift,sumsPGSSn,imaxPGSS);
[xmaxPGSIn,ymaxPGSIn,xloPGSIn,xhiPGSIn]=parabol(nwin,mshift,sumsPGSIn,imaxPGSI);
[xmaxSISSn,ymaxSISSn,xloSISSn,xhiSISSn]=parabol(nwin,mshift,sumsSISSn,imaxSISS);
[xmaxPGSSin,ymaxPGSSin,xloPGSSin,xhiPGSSin]=parabol(nwin,mshift,sumsPGSSin,imaxPGSSi);
[xmaxPGSIin,ymaxPGSIin,xloPGSIin,xhiPGSIin]=parabol(nwin,mshift,sumsPGSIin,imaxPGSIi);
[xmaxSISSin,ymaxSISSin,xloSISSin,xhiSISSin]=parabol(nwin,mshift,sumsSISSin,imaxSISSi);
ix=sub2ind(size(denomPGSSn),(1:nwin)',imaxPGSS);
ixi=sub2ind(size(denomPGSSin),(1:nwin)',imaxPGSSi);
ampPGSS=sqrt(denomPGSSn(ix)); %This makes amplitude linear rather than quadratic with counts.
ampPGSSi=sqrt(denomPGSSin(ixi)); %This makes amplitude linear rather than quadratic with counts.
xcmaxPGSSin=sumsPGSSin(ix); %This evaluates orthogonal x-correlation value at main maximum.

ix=sub2ind(size(denomPGSIn),(1:nwin)',imaxPGSI);
ixi=sub2ind(size(denomPGSIin),(1:nwin)',imaxPGSIi);
ampPGSI=sqrt(denomPGSIn(ix));
ampPGSIi=sqrt(denomPGSIin(ixi));
xcmaxPGSIin=sumsPGSIin(ix);

ix=sub2ind(size(denomSISSn),(1:nwin)',imaxSISS);
ixi=sub2ind(size(denomSISSin),(1:nwin)',imaxSISSi);
ampSISS=sqrt(denomSISSn(ix));
ampSISSi=sqrt(denomSISSin(ixi));
xcmaxSISSin=sumsSISSin(ix);
%Center them
imaxPGSS=imaxPGSS-mshift-1;
imaxPGSI=imaxPGSI-mshift-1;
imaxSISS=imaxSISS-mshift-1;
xmaxPGSSn=xmaxPGSSn-mshift-1;
xmaxPGSIn=xmaxPGSIn-mshift-1;
xmaxSISSn=xmaxSISSn-mshift-1;
imaxPGSSi=imaxPGSSi-mshift-1;
imaxPGSIi=imaxPGSIi-mshift-1;
imaxSISSi=imaxSISSi-mshift-1;
xmaxPGSSin=xmaxPGSSin-mshift-1;
xmaxPGSIin=xmaxPGSIin-mshift-1;
xmaxSISSin=xmaxSISSin-mshift-1;

ampmax=max([ampPGSS; ampPGSI; ampSISS]);
figure 
subplot(4,1,1); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,2.5+zeros(nwin,1),'k:');
plot(timswin,ampSISS*5/ampmax,'k+','MarkerSize',2);
plot(timswin,ampPGSS*5/ampmax,'b+','MarkerSize',2);
plot(timswin,ampPGSI*5/ampmax,'r+','MarkerSize',2);
plot(timswin,xcmaxSISSn*5,'k');
plot(timswin,xcmaxPGSSn*5,'b');
plot(timswin,xcmaxPGSIn*5,'r');
plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',4);
plot(timswin,xmaxPGSSn,'bs','MarkerSize',4);
plot(timswin,xmaxPGSIn,'ro','MarkerSize',4);
axis([timlook-timbig timlook-timbig/2 -5 5]);
xlabel('sec')
ylabel('samples')
box on
subplot(4,1,2); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,2.5+zeros(nwin,1),'k:');
plot(timswin,ampSISS*5/ampmax,'k+','MarkerSize',2);
plot(timswin,ampPGSS*5/ampmax,'b+','MarkerSize',2);
plot(timswin,ampPGSI*5/ampmax,'r+','MarkerSize',2);
plot(timswin,xcmaxSISSn*5,'k');
plot(timswin,xcmaxPGSSn*5,'b');
plot(timswin,xcmaxPGSIn*5,'r');
plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',4);
plot(timswin,xmaxPGSSn,'bs','MarkerSize',4);
plot(timswin,xmaxPGSIn,'ro','MarkerSize',4);
axis([timlook-timbig/2 timlook -5 5]);
xlabel('sec')
ylabel('samples')
box on
subplot(4,1,3); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,2.5+zeros(nwin,1),'k:');
plot(timswin,ampSISS*5/ampmax,'k+','MarkerSize',2);
plot(timswin,ampPGSS*5/ampmax,'b+','MarkerSize',2);
plot(timswin,ampPGSI*5/ampmax,'r+','MarkerSize',2);
plot(timswin,xcmaxSISSn*5,'k');
plot(timswin,xcmaxPGSSn*5,'b');
plot(timswin,xcmaxPGSIn*5,'r');
plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',4);
plot(timswin,xmaxPGSSn,'bs','MarkerSize',4);
plot(timswin,xmaxPGSIn,'ro','MarkerSize',4);
axis([timlook timlook+timbig/2 -5 5]);
xlabel('sec')
ylabel('samples')
box on
subplot(4,1,4); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,2.5+zeros(nwin,1),'k:');
plot(timswin,ampPGSS*5/ampmax,'k+','MarkerSize',2);
plot(timswin,ampPGSS*5/ampmax,'b+','MarkerSize',2);
plot(timswin,ampPGSI*5/ampmax,'r+','MarkerSize',2);
plot(timswin,xcmaxSISSn*5,'k');
plot(timswin,xcmaxPGSSn*5,'b');
plot(timswin,xcmaxPGSIn*5,'r');
plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',4);
plot(timswin,xmaxPGSSn,'bs','MarkerSize',4);
plot(timswin,xmaxPGSIn,'ro','MarkerSize',4);
axis([timlook+timbig/2 timlook+timbig -5 5]);
xlabel('sec')
ylabel('samples')
box on
orient landscape
print('-depsc',[JDAY,'_',int2str(timlook),'_wiggli2.eps'])

figure 
subplot(4,1,1); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,2.5+zeros(nwin,1),'k:');
plot(timswin,ampSISSi*5/ampmax,'k+','MarkerSize',2);
plot(timswin,ampPGSSi*5/ampmax,'b+','MarkerSize',2);
plot(timswin,ampPGSIi*5/ampmax,'r+','MarkerSize',2);
plot(timswin,xcmaxSISSin*5,'k');
plot(timswin,xcmaxPGSSin*5,'b');
plot(timswin,xcmaxPGSIin*5,'r');
plot(timswin,(xmaxPGSIin-xmaxPGSSin+xmaxSISSin),'k*','MarkerSize',4);
plot(timswin,xmaxPGSSin,'bs','MarkerSize',4);
plot(timswin,xmaxPGSIin,'ro','MarkerSize',4);
axis([timlook-timbig timlook-timbig/2 -5 5]);
xlabel('sec')
ylabel('samples')
box on
subplot(4,1,2); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,2.5+zeros(nwin,1),'k:');
plot(timswin,ampSISSi*5/ampmax,'k+','MarkerSize',2);
plot(timswin,ampPGSSi*5/ampmax,'b+','MarkerSize',2);
plot(timswin,ampPGSIi*5/ampmax,'r+','MarkerSize',2);
plot(timswin,xcmaxSISSin*5,'k');
plot(timswin,xcmaxPGSSin*5,'b');
plot(timswin,xcmaxPGSIin*5,'r');
plot(timswin,(xmaxPGSIin-xmaxPGSSin+xmaxSISSin),'k*','MarkerSize',4);
plot(timswin,xmaxPGSSin,'bs','MarkerSize',4);
plot(timswin,xmaxPGSIin,'ro','MarkerSize',4);
axis([timlook-timbig/2 timlook -5 5]);
xlabel('sec')
ylabel('samples')
box on
subplot(4,1,3); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,2.5+zeros(nwin,1),'k:');
plot(timswin,ampSISSi*5/ampmax,'k+','MarkerSize',2);
plot(timswin,ampPGSSi*5/ampmax,'b+','MarkerSize',2);
plot(timswin,ampPGSIi*5/ampmax,'r+','MarkerSize',2);
plot(timswin,xcmaxSISSin*5,'k');
plot(timswin,xcmaxPGSSin*5,'b');
plot(timswin,xcmaxPGSIin*5,'r');
plot(timswin,(xmaxPGSIin-xmaxPGSSin+xmaxSISSin),'k*','MarkerSize',4);
plot(timswin,xmaxPGSSin,'bs','MarkerSize',4);
plot(timswin,xmaxPGSIin,'ro','MarkerSize',4);
axis([timlook timlook+timbig/2 -5 5]);
xlabel('sec')
ylabel('samples')
box on
subplot(4,1,4); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,2.5+zeros(nwin,1),'k:');
plot(timswin,ampPGSSi*5/ampmax,'k+','MarkerSize',2);
plot(timswin,ampPGSSi*5/ampmax,'b+','MarkerSize',2);
plot(timswin,ampPGSIi*5/ampmax,'r+','MarkerSize',2);
plot(timswin,xcmaxSISSin*5,'k');
plot(timswin,xcmaxPGSSin*5,'b');
plot(timswin,xcmaxPGSIin*5,'r');
plot(timswin,(xmaxPGSIin-xmaxPGSSin+xmaxSISSin),'k*','MarkerSize',4);
plot(timswin,xmaxPGSSin,'bs','MarkerSize',4);
plot(timswin,xmaxPGSIin,'ro','MarkerSize',4);
axis([timlook+timbig/2 timlook+timbig -5 5]);
xlabel('sec')
ylabel('samples')
box on
orient landscape
print('-depsc',[JDAY,'_',int2str(timlook),'_wiggli3.eps'])

figure %Superposed traces, 150s window:
subplot(4,1,1); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([timlook-timbig timlook-timbig/2 ltmin ltmax]);
subplot(4,1,2); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,2.5+zeros(nwin,1),'k:');
% plot(timswin,ampSISS*5/ampmax,'k+','MarkerSize',5);
% plot(timswin,ampPGSS*5/ampmax,'b+','MarkerSize',5);
% plot(timswin,ampPGSI*5/ampmax,'r+','MarkerSize',5);
plot(timswin,xcmaxSISSn*5,'k');
plot(timswin,xcmaxPGSSn*5,'b');
plot(timswin,xcmaxPGSIn*5,'r');
plot(timswin,xmaxPGSSn,'bs','MarkerSize',5);
plot(timswin,xmaxPGSIn,'ro','MarkerSize',5);
plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',5);
axis([timlook-timbig timlook-timbig/2 -5 5]);
subplot(4,1,3); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([timlook-timbig/2 timlook ltmin ltmax]);
subplot(4,1,4); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,2.5+zeros(nwin,1),'k:');
% plot(timswin,ampSISS*5/ampmax,'k+','MarkerSize',5);
% plot(timswin,ampPGSS*5/ampmax,'b+','MarkerSize',5);
% plot(timswin,ampPGSI*5/ampmax,'r+','MarkerSize',5);
plot(timswin,xcmaxSISSn*5,'k');
plot(timswin,xcmaxPGSSn*5,'b');
plot(timswin,xcmaxPGSIn*5,'r');
plot(timswin,xmaxPGSSn,'bs','MarkerSize',5);
plot(timswin,xmaxPGSIn,'ro','MarkerSize',5);
plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',5);
axis([timlook-timbig/2 timlook -5 5]);
orient landscape
print('-dpdf',[JDAY,'_',int2str(timlook),'wiggli4.pdf'])

figure
subplot(4,1,1); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([timlook timlook+timbig/2 ltmin ltmax]);
subplot(4,1,2); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,2.5+zeros(nwin,1),'k:');
% plot(timswin,ampSISS*5/ampmax,'k+','MarkerSize',5);
% plot(timswin,ampPGSS*5/ampmax,'b+','MarkerSize',5);
% plot(timswin,ampPGSI*5/ampmax,'r+','MarkerSize',5);
plot(timswin,xcmaxSISSn*5,'k--');
plot(timswin,xcmaxPGSSn*5,'b--');
plot(timswin,xcmaxPGSIn*5,'r--');
plot(timswin,xmaxPGSSn,'bs','MarkerSize',5);
plot(timswin,xmaxPGSIn,'ro','MarkerSize',5);
plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',5);
axis([timlook timlook+timbig/2 -5 5]);
subplot(4,1,3); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([timlook+timbig/2 timlook+timbig ltmin ltmax]);
subplot(4,1,4); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,2.5+zeros(nwin,1),'k:');
% plot(timswin,ampSISS*5/ampmax,'k+','MarkerSize',5);
% plot(timswin,ampPGSS*5/ampmax,'b+','MarkerSize',5);
% plot(timswin,ampPGSI*5/ampmax,'r+','MarkerSize',5);
plot(timswin,xcmaxSISSn*5,'k--');
plot(timswin,xcmaxPGSSn*5,'b--');
plot(timswin,xcmaxPGSIn*5,'r--');
plot(timswin,xmaxPGSSn,'bs','MarkerSize',5);
plot(timswin,xmaxPGSIn,'ro','MarkerSize',5);
plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',5);
axis([timlook+timbig/2 timlook+timbig -5 5]);
orient landscape
print('-dpdf',[JDAY,'_',int2str(timlook),'wiggli5.pdf'])


%cputime-t;
tot=cputime-tt
