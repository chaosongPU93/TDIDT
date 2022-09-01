%Sits at one spot (rotations; delays) and looks for the cross-correlation.
t=cputime; tt=cputime;
close all
format short e
%Rotations & offsets:
   timoffrot=[247 15 35 08 -073 -146  90  50 385]; %(the original; catalog has 85  75  55)
%   timoffrot=[250 11 34 20 +086 -032  80 105  65]; %migration
 %  timoffrot=[250 11 34 28 -032 -142  80  60 385];
   %timoffrot=[259 12 28 44 +044 +071  70 350 375];
   %timoffrot=[259 13 29 24 +036 +075  75 370 390];
jday=timoffrot(1);
JDAY=int2str(jday)
timlook=3600*timoffrot(2)+60*timoffrot(3)+timoffrot(4)
rotPGC=pi*(timoffrot(7)-90)/180;
rotSILB=pi*(timoffrot(8)-90)/180;
rotVGZ=pi*(timoffrot(9)-90)/180;
SILBsoff=timoffrot(5);
VGZsoff=timoffrot(6);
SILBtoff=SILBsoff/40;
VGZtoff=VGZsoff/40;
%
% SILBsoff=51;
% VGZsoff=-78;
% SILBsoff=93;
% VGZsoff=-24;

%Read Armbruster's detections:
ArmCat=load('/data2/arubin/ARMB/BC2005/list.2005.pgsssiMAJ');
%ArmCat=load('/data2/arubin/ARMB/BC2005/list.2005.pgsssiMIN');
detects=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
vectoff=[ArmCat(:,5)-SILBsoff ArmCat(:,6)-VGZsoff];
distoff=vectoff.*vectoff;
distoff2=sum(distoff,2);
m=0;
Detects=0;
for n=1:length(detects)
     if distoff2(n)<=9
%      if isequal(ArmCat(n,5:6),[SILBsoff VGZsoff]) ||...
%         isequal(ArmCat(n,5:6),[SILBsoff-1 VGZsoff]) ||...
%         isequal(ArmCat(n,5:6),[SILBsoff+1 VGZsoff]) ||...
%         isequal(ArmCat(n,5:6),[SILBsoff VGZsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SILBsoff VGZsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SILBsoff-1 VGZsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SILBsoff-1 VGZsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SILBsoff+1 VGZsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SILBsoff+1 VGZsoff+1]) 
         m=m+1;
         Detects(m)=detects(n);
     end
end
ArmCat=load('/data2/arubin/ARMB/BC2005/list.2005.pgsssiMAJ_new');
detects2_8=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
vectoff=[ArmCat(:,5)-SILBsoff ArmCat(:,6)-VGZsoff];
distoff=vectoff.*vectoff;
distoff2=sum(distoff,2);
m=0;
Detects2_8=0;
for n=1:length(detects2_8)
    if distoff2(n)<=9
%      if isequal(ArmCat(n,5:6),[SILBsoff VGZsoff]) ||...
%         isequal(ArmCat(n,5:6),[SILBsoff-1 VGZsoff]) ||...
%         isequal(ArmCat(n,5:6),[SILBsoff+1 VGZsoff]) ||...
%         isequal(ArmCat(n,5:6),[SILBsoff VGZsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SILBsoff VGZsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SILBsoff-1 VGZsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SILBsoff-1 VGZsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SILBsoff+1 VGZsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SILBsoff+1 VGZsoff+1]) 
         m=m+1;
         Detects2_8(m)=detects2_8(n);
     end
end

%Read data:
direc='2005/SEPT/';
prename=[direc,'2005.',JDAY,'.00.00.00.0000.CN'];
PGCEdat=[prename,'.PGC..BHE.D.SAC'];
PGCNdat=[prename,'.PGC..BHN.D.SAC'];
%PGCZdat=[prename,'.PGC..BHZ.D.SAC'];
SILBEdat=[prename,'.SILB..HHE.D.SAC'];
SILBNdat=[prename,'.SILB..HHN.D.SAC'];
VGZEdat=[prename,'.VGZ..BHE.D.SAC'];
VGZNdat=[prename,'.VGZ..BHN.D.SAC'];

[PGCE,HdrDataPGC,tnuPGC,pobjPGC,timsPGC]=readsac(PGCEdat,0,'l');
[PGCN,~,~,~,~]=readsac(PGCNdat,0,'l');
%[PGCZ,~,~,~,~]=readsac(PGCZdat,0,'l');
[SILBE,HdrDataSILB,tnuSILB,pobjSILB,timsSILB]=readsac(SILBEdat,0,'l');
[SILBN,~,~,~,~]=readsac(SILBNdat,0,'l');
[VGZE,HdrDataVGZ,tnuVGZ,pobjVGZ,timsVGZ]=readsac(VGZEdat,0,'l');
[VGZN,~,~,~,~]=readsac(VGZNdat,0,'l');
% timsPGC(1)
% timsSILB(1)
% timsVGZ(1)
%
rdsac=cputime-t
t=cputime;
%     % Truncate for quicker calculation.  100sps stations are longer for
%     % later interpolation
%     timwin=4*60+1;
%     timstart=max(0, timlook-timwin); timend=min(86400, timlook+timwin);
%     PGCE=PGCE(timstart*40+1:timend*40)-mean(PGCE(timstart*40+1:timend*40));
%     PGCN=PGCN(timstart*40+1:timend*40)-mean(PGCN(timstart*40+1:timend*40)); 
%     %PGCZ=PGCZ(timstart*40+1:timend*40);
%     SILBE=SILBE(timstart*100:timend*100+1)-mean(SILBE(timstart*100:timend*100+1));
%     SILBN=SILBN(timstart*100:timend*100+1)-mean(SILBN(timstart*100:timend*100+1));
%     VGZE=VGZE(timstart*100:timend*100+1)-mean(VGZE(timstart*100:timend*100+1));
%     VGZN=VGZN(timstart*100:timend*100+1)-mean(VGZN(timstart*100:timend*100+1));
%     timsPGC=timsPGC(timstart*40+1:timend*40);
%     timsSILB=timsSILB(timstart*100:timend*100+1);
%     timsVGZ=timsVGZ(timstart*100:timend*100+1);
% %shrten=cputime-t
% %t=cputime;
% %
tracelenPGC=length(PGCE);
tracelenSILB=length(SILBE);
tracelenVGZ=length(VGZE);
    timPGCfirst=timsPGC(1);
    timSILBfirst=timsSILB(1); %0.01s earlier than PGC
    timVGZfirst=timsVGZ(1); %0.01s earlier than PGC
    timPGClast=timsPGC(tracelenPGC);
    timSILBlast=timsSILB(tracelenSILB);
    timVGZlast=timsVGZ(tracelenVGZ);
    timstart=timPGCfirst; timend=timPGClast; 
    timwin=(timPGClast-timPGCfirst)/2;
    startdiff=timPGCfirst-timSILBfirst
    enddiff=timSILBlast-timPGClast

%cosine taper before filtering:
x=(0:pi/80:pi/2-pi/80)';
% %Seems to be necessary at the start of each day for PGCE:
    PGCE(1:80)=0.;
    PGCE(81:120)=sin(x).*PGCE(81:120); %Only at start of day!
    VGZE(1:80)=0.;
    VGZE(81:120)=sin(x).*VGZE(81:120); %Only at start of day!
%PGCE(1:40)=sin(x).*PGCE(1:40);
PGCN(1:40)=sin(x).*PGCN(1:40);
VGZN(1:40)=sin(x).*VGZN(1:40);
x=flipud(x);
PGCE(tracelenPGC-39:tracelenPGC)=sin(x).*PGCE(tracelenPGC-39:tracelenPGC);
PGCN(tracelenPGC-39:tracelenPGC)=sin(x).*PGCN(tracelenPGC-39:tracelenPGC);
VGZE(tracelenVGZ-39:tracelenVGZ)=sin(x).*VGZE(tracelenVGZ-39:tracelenVGZ);
VGZN(tracelenVGZ-39:tracelenVGZ)=sin(x).*VGZN(tracelenVGZ-39:tracelenVGZ);
x=(0:pi/200:pi/2-pi/200)';
SILBE(1:100)=sin(x).*SILBE(1:100);
SILBN(1:100)=sin(x).*SILBN(1:100);
x=flipud(x);
SILBE(tracelenSILB-99:tracelenSILB)=sin(x).*SILBE(tracelenSILB-99:tracelenSILB);
SILBN(tracelenSILB-99:tracelenSILB)=sin(x).*SILBN(tracelenSILB-99:tracelenSILB);

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
[SILBEf]=4.0e-3*bandpass(SILBE,100,lo,hi,npo,npa,'butter');
[SILBNf]=4.0e-3*bandpass(SILBN,100,lo,hi,npo,npa,'butter');
[VGZEf]=8*1.6e-4*bandpass(VGZE,40,lo,hi,npo,npa,'butter');
[VGZNf]=8*1.6e-4*bandpass(VGZN,40,lo,hi,npo,npa,'butter');
%fltr=cputime-t
%t=cputime;
%Decimate the 100 sps data:
SILBEfd = interp1(timsSILB,SILBEf,timsPGC,'linear')';
SILBNfd = interp1(timsSILB,SILBNf,timsPGC,'linear')';
%dmate=cputime-t
%t=cputime;
% SILBEfd(1:6)
% SILBNfd(1:6)
% VGZEfd(1:6)
% VGZNfd(1:6)

PGC=PGCEf+1i*PGCNf;
SILB=SILBEfd+1i*SILBNfd;
VGZ=VGZEf+1i*VGZNf;
PGCrot=PGC*exp(1i*rotPGC);
SILBrot=SILB*exp(1i*rotSILB);
VGZrot=VGZ*exp(1i*rotVGZ);
    %rtate=cputime-t
    %t=cputime;

timsSILB=timsPGC-SILBtoff; %No longer used, I think
timsVGZ=timsPGC-VGZtoff; %No longer used, I think
    %tshft=cputime-t
    %t=cputime;

if SILBtoff > 0
    SILBrot(1:tracelenPGC-SILBsoff)=SILBrot(SILBsoff+1:tracelenPGC);
    SILBrot(tracelenPGC-SILBsoff+1:tracelenPGC)=0;
else
    SILBrot(-SILBsoff+1:tracelenPGC)=SILBrot(1:tracelenPGC+SILBsoff);
    SILBrot(1:-SILBsoff)=0;
end
if VGZtoff > 0
    VGZrot(1:tracelenPGC-VGZsoff)=VGZrot(VGZsoff+1:tracelenPGC);
    VGZrot(tracelenPGC-VGZsoff+1:tracelenPGC)=0;
else
    VGZrot(-VGZsoff+1:tracelenPGC)=VGZrot(1:tracelenPGC+VGZsoff);
    VGZrot(1:-VGZsoff)=0;
end
    %tracesht=cputime-t
    %t=cputime;

realPGC=real(PGCrot);
realSILB=real(SILBrot);
realVGZ=real(VGZrot);
% realPGC=imag(PGCrot); %A quick/dirty way to try the same on the orthogonal components
% realSILB=imag(SILBrot);
% realVGZ=imag(VGZrot);

plotstart=1; %(timlook-timstart)*40
plotend=tracelenPGC; %(timlook-timstart)*40
ltmax=max([real(PGCrot(plotstart:plotend)); ...
          real(SILBrot(plotstart:plotend)); ...
          real(VGZrot(plotstart:plotend))]);
ltmax=min(ltmax,0.75);
ltmin=min([real(PGCrot(plotstart:plotend)); ...
          real(SILBrot(plotstart:plotend)); ...
          real(VGZrot(plotstart:plotend))]);
ltmin=max(ltmin,-0.75);

figure %Superposed traces, 150s window:
subplot(4,1,1); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SILBrot),'b');
plot(timsPGC,real(VGZrot),'k');
axis([0 timwin/2 ltmin ltmax]);
%xlim([timlook-75 timlook-37.5]);
subplot(4,1,2); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SILBrot),'b');
plot(timsPGC,real(VGZrot),'k');
axis([timwin/2 timwin ltmin ltmax]);
%xlim([timlook-37.5 timlook]);
subplot(4,1,3); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SILBrot),'b');
plot(timsPGC,real(VGZrot),'k');
axis([timwin 3*timwin/2 ltmin ltmax]);
%xlim([timlook timlook+37.5]);
subplot(4,1,4); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SILBrot),'b');
plot(timsPGC,real(VGZrot),'k');
axis([3*timwin/2 2*timwin ltmin ltmax]);
%xlim([timlook+37.5 timlook+75]);
traceplt=cputime-t
t=cputime;

%Autocorrelation of stations:
PGCauto=realPGC.*realPGC;
PGC2=cumsum(PGCauto);
SILBauto=realSILB.*realSILB;
SILB2=cumsum(SILBauto);
VGZauto=realVGZ.*realVGZ;
VGZ2=cumsum(VGZauto);

mshift=11;
lenx=tracelenPGC-2*mshift;
PGSIx=zeros(lenx, 2*mshift+1);
PGVGx=zeros(lenx, 2*mshift+1);
VGSIx=zeros(lenx, 2*mshift+1);
% First index is pointwise multiplication of traces; second is shifting offset.
for n=-mshift:mshift;
    PGSIx(:,n+mshift+1)=realPGC(1+mshift:tracelenPGC-mshift).* ...
        realSILB(1+mshift-n:tracelenPGC-mshift-n);
    PGVGx(:,n+mshift+1)=realPGC(1+mshift:tracelenPGC-mshift).* ...
        realVGZ(1+mshift-n:tracelenPGC-mshift-n);
    VGSIx(:,n+mshift+1)=realVGZ(1+mshift:tracelenPGC-mshift).* ...
        realSILB(1+mshift-n:tracelenPGC-mshift-n);
end
cumsumPGSI=cumsum(PGSIx);
cumsumPGVG=cumsum(PGVGx);
cumsumVGSI=cumsum(VGSIx);
xcshifts=cputime-t
t=cputime;

winbig=2*(tracelenPGC/2-120);
timbig=winbig/(2*40);
igstart=floor(tracelenPGC/2-winbig/2)+1;
winlen=4*40;
winoff=40/1;
nwin=floor((winbig-winlen)/winoff);
timswin=zeros(nwin,1);
sumsPGSI=zeros(nwin,2*mshift+1);
sumsPGVG=zeros(nwin,2*mshift+1);
sumsVGSI=zeros(nwin,2*mshift+1);
sumsPGC2=zeros(nwin,2*mshift+1);
sumsSILB2=zeros(nwin,2*mshift+1);
sumsVGZ2=zeros(nwin,2*mshift+1);
sumsVGZ2b=zeros(nwin,2*mshift+1);
for n=1:nwin;
    istart=igstart+(n-1)*winoff;
    iend=istart+winlen;
    timswin(n)=timsPGC(istart+winlen/2); %Gets rid of "mshifts" by referencing to the global igstart.
    %First index is window #; second is shifts over +-mshift.
    sumsPGSI(n,:)=cumsumPGSI(iend,:)-cumsumPGSI(istart-1,:); %Summed x-correlations
    sumsPGVG(n,:)=cumsumPGVG(iend,:)-cumsumPGVG(istart-1,:);
    sumsVGSI(n,:)=cumsumVGSI(iend,:)-cumsumVGSI(istart-1,:);
    sumsPGC2(n,:)=PGC2(iend+mshift)-PGC2(istart+mshift-1);  %PGC doesn't shift, so do this here.  PGC2 is cumsummed. Yes, +mshift.
    sumsVGZ2b(n,:)=VGZ2(iend+mshift)-VGZ2(istart+mshift-1); %Similar, for the VGZ-SILB connection.
    for m=-mshift:mshift;
        sumsSILB2(n,m+mshift+1)=SILB2(iend+mshift-m)-SILB2(istart+mshift-1-m); %+m??? (yes).
        sumsVGZ2(n,m+mshift+1)=VGZ2(iend+mshift-m)-VGZ2(istart+mshift-1-m);
    end
end
denomPGSIn=realsqrt(sumsPGC2.*sumsSILB2);
denomPGVGn=realsqrt(sumsPGC2.*sumsVGZ2);
denomVGSIn=realsqrt(sumsVGZ2b.*sumsSILB2);
sumsPGSIn=sumsPGSI./denomPGSIn;
%sumsPGSIn'
sumsPGVGn=sumsPGVG./denomPGVGn;
sumsVGSIn=sumsVGSI./denomVGSIn;
[xcmaxPGSIn,imaxPGSI]=max(sumsPGSIn,[],2);
[xcmaxPGVGn,imaxPGVG]=max(sumsPGVGn,[],2);
[xcmaxVGSIn,imaxVGSI]=max(sumsVGSIn,[],2);
%Parabolic fit:
% [xmaxPGSIn,ymaxPGSIn,xloPGSIn,xhiPGSIn]=parabol(nwin,mshift,sumsPGSIn,imaxPGSI);
% [xmaxPGVGn,ymaxPGVGn,xloPGVGn,xhiPGVGn]=parabol(nwin,mshift,sumsPGVGn,imaxPGVG);
% [xmaxVGSIn,ymaxVGSIn,xloVGSIn,xhiVGSIn]=parabol(nwin,mshift,sumsVGSIn,imaxVGSI);
[xmaxPGSIn,ymaxPGSIn]=parabol(nwin,mshift,sumsPGSIn,imaxPGSI);
[xmaxPGVGn,ymaxPGVGn]=parabol(nwin,mshift,sumsPGVGn,imaxPGVG);
[xmaxVGSIn,ymaxVGSIn]=parabol(nwin,mshift,sumsVGSIn,imaxVGSI);
ix=sub2ind(size(denomPGSIn),(1:nwin)',imaxPGSI);
ampPGSI=sqrt(denomPGSIn(ix)); %This makes amplitude linear rather than quadratic with counts.
ampPGC2=sumsPGC2(ix); %by construction PGC2 is the same for all shifts
ampSILB2=sumsSILB2(ix);
ix=sub2ind(size(denomPGVGn),(1:nwin)',imaxPGVG);
ampPGVG=sqrt(denomPGVGn(ix));
ampVGZ2=sumsVGZ2(ix);
ix=sub2ind(size(denomVGSIn),(1:nwin)',imaxVGSI);
ampVGSI=sqrt(denomVGSIn(ix));
AmpComp(1:4)=0;
AmpComp(5:nwin)=((ampPGC2(5:nwin)+ampSILB2(5:nwin)+ampVGZ2(5:nwin))- ...
                (ampPGC2(1:nwin-4)+ampSILB2(1:nwin-4)+ampVGZ2(1:nwin-4)))./ ...
                ((ampPGC2(5:nwin)+ampSILB2(5:nwin)+ampVGZ2(5:nwin))+ ...
                (ampPGC2(1:nwin-4)+ampSILB2(1:nwin-4)+ampVGZ2(1:nwin-4))) ;
%Center them
imaxPGSI=imaxPGSI-mshift-1;
imaxPGVG=imaxPGVG-mshift-1;
imaxVGSI=imaxVGSI-mshift-1;
xmaxPGSIn=xmaxPGSIn-mshift-1;
xmaxPGVGn=xmaxPGVGn-mshift-1;
xmaxVGSIn=xmaxVGSIn-mshift-1;
xcmaxAVEn=(xcmaxPGSIn+xcmaxPGVGn+xcmaxVGSIn)/3;
xcnshifts=cputime-t
t=cputime;

ampmax=max([ampPGSI; ampPGVG; ampVGSI]);
figure 
subplot(4,1,1,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
%plot(timswin,3.5+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*7,'g');
plot(timswin,xmaxPGSIn,'bs','MarkerSize',2);
plot(timswin,xmaxPGVGn,'ro','MarkerSize',2);
plot(timswin,(xmaxPGVGn-xmaxPGSIn+xmaxVGSIn),'k*','MarkerSize',2);
%plot(timswin,xmaxVGSIn,'k*','MarkerSize',2);
axis([0 timbig/2 -7 7]);
xlabel('sec')
ylabel('samples')
box on
subplot(4,1,2,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
%plot(timswin,3.5+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*7,'g');
plot(timswin,xmaxPGSIn,'bs','MarkerSize',2);
plot(timswin,xmaxPGVGn,'ro','MarkerSize',2);
plot(timswin,(xmaxPGVGn-xmaxPGSIn+xmaxVGSIn),'k*','MarkerSize',2);
%plot(timswin,xmaxVGSIn,'k*','MarkerSize',2);
axis([timbig/2 timbig -7 7]);
xlabel('sec')
ylabel('samples')
box on
subplot(4,1,3,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
%plot(timswin,3.5+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*7,'g');
plot(timswin,xmaxPGSIn,'bs','MarkerSize',2);
plot(timswin,xmaxPGVGn,'ro','MarkerSize',2);
plot(timswin,(xmaxPGVGn-xmaxPGSIn+xmaxVGSIn),'k*','MarkerSize',2);
%plot(timswin,xmaxVGSIn,'k*','MarkerSize',2);
axis([timbig 3*timbig/2 -7 7]);
xlabel('sec')
ylabel('samples')
box on
subplot(4,1,4,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
%plot(timswin,3.5+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*7,'g');
plot(timswin,xmaxPGSIn,'bs','MarkerSize',2);
plot(timswin,xmaxPGVGn,'ro','MarkerSize',2);
plot(timswin,(xmaxPGVGn-xmaxPGSIn+xmaxVGSIn),'k*','MarkerSize',2);
%plot(timswin,xmaxVGSIn,'k*','MarkerSize',2);
axis([3*timbig/2 2*timbig -7 7]);
xlabel('sec')
ylabel('samples')
box on
% orient landscape
% print('-depsc',[JDAY,'.eps'])

xmaxPGSIntmp=xmaxPGSIn;
xmaxPGVGntmp=xmaxPGVGn;
xmaxVGSIntmp=xmaxVGSIn;
for n=1:nwin;
    if xcmaxAVEn(n)<0.3 || abs(xmaxPGVGn(n)-xmaxPGSIn(n)+xmaxVGSIn(n))>2
        xmaxPGSIntmp(n)=20; xmaxPGVGntmp(n)=20; xmaxVGSIntmp(n)=20; 
    end
end
figure 
subplot(4,1,1,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*7,'g');
plot(timswin,xmaxPGSIntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGVGntmp,'ro','MarkerSize',2);
plot(timswin,xmaxVGSIntmp,'k*','MarkerSize',2);
axis([0 timbig/2 -7 7]);
xlabel('sec')
ylabel('samples')
box on
subplot(4,1,2,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*7,'g');
plot(timswin,xmaxPGSIntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGVGntmp,'ro','MarkerSize',2);
plot(timswin,xmaxVGSIntmp,'k*','MarkerSize',2);
axis([timbig/2 timbig -7 7]);
xlabel('sec')
ylabel('samples')
box on
subplot(4,1,3,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*7,'g');
plot(timswin,xmaxPGSIntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGVGntmp,'ro','MarkerSize',2);
plot(timswin,xmaxVGSIntmp,'k*','MarkerSize',2);
axis([timbig 3*timbig/2 -7 7]);
xlabel('sec')
ylabel('samples')
box on
subplot(4,1,4,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*7,'g');
plot(timswin,xmaxPGSIntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGVGntmp,'ro','MarkerSize',2);
plot(timswin,xmaxVGSIntmp,'k*','MarkerSize',2);
axis([3*timbig/2 2*timbig -7 7]);
xlabel('sec')
ylabel('samples')
box on
%orient landscape
%print('-depsc',[JDAY,'.eps'])

loopoff=xmaxPGVGn-xmaxPGSIn+xmaxVGSIn;
for n=1:nwin;
    if abs(loopoff(n))>2.0
        loopoff(n)=20; 
    end
end
figure 
subplot(4,1,1,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xmaxPGSIntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGVGntmp,'ro','MarkerSize',2);
plot(timswin,xmaxVGSIntmp,'k*','MarkerSize',2);
axis([0 timbig/2 -6 6]);
subplot(4,1,2,'align'); 
hold on
hrf = plotreflinesr(gca,detects,'x','k'); 
hrf = plotreflinesr(gca,Detects,'x','r'); 
hrf = plotreflinesr(gca,detects2_8,'x','g'); 
hrf = plotreflinesr(gca,Detects2_8,'x','b'); 
plot(timswin,0.5+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn,'k');
plot(timswin,0.1*abs(loopoff),'go','MarkerSize',2);
plot(timsPGC,min(0,real(PGCrot)),'r');
plot(timsPGC,min(0,real(SILBrot)),'b');
plot(timsPGC,min(0,real(VGZrot)),'k');
axis([0 timbig/2 -1 1]);
subplot(4,1,3,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xmaxPGSIntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGVGntmp,'ro','MarkerSize',2);
plot(timswin,xmaxVGSIntmp,'k*','MarkerSize',2);
axis([timbig/2 timbig -6 6]);
subplot(4,1,4,'align'); 
hold on
hrf = plotreflinesr(gca,detects,'x','k'); 
hrf = plotreflinesr(gca,Detects,'x','r'); 
hrf = plotreflinesr(gca,detects2_8,'x','g'); 
hrf = plotreflinesr(gca,Detects2_8,'x','b'); 
plot(timswin,0.5+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn,'k');
plot(timswin,0.1*abs(loopoff),'go','MarkerSize',2);
plot(timsPGC,min(0,real(PGCrot)),'r');
plot(timsPGC,min(0,real(SILBrot)),'b');
plot(timsPGC,min(0,real(VGZrot)),'k');
axis([timbig/2 timbig -1 1]);
%orient landscape
%print('-dpdf',[JDAY,'_',int2str(timlook),'_1.pdf'])

figure
subplot(4,1,1,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xmaxPGSIntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGVGntmp,'ro','MarkerSize',2);
plot(timswin,xmaxVGSIntmp,'k*','MarkerSize',2);
axis([timbig 3*timbig/2 -6 6]);
subplot(4,1,2,'align'); 
hold on
hrf = plotreflinesr(gca,detects,'x','k'); 
hrf = plotreflinesr(gca,Detects,'x','r'); 
hrf = plotreflinesr(gca,detects2_8,'x','g'); 
hrf = plotreflinesr(gca,Detects2_8,'x','b'); 
plot(timswin,0.5+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn,'k');
plot(timswin,0.1*abs(loopoff),'go','MarkerSize',2);
plot(timsPGC,min(0,real(PGCrot)),'r');
plot(timsPGC,min(0,real(SILBrot)),'b');
plot(timsPGC,min(0,real(VGZrot)),'k');
axis([timbig 3*timbig/2 -1 1]);
subplot(4,1,3,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xmaxPGSIntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGVGntmp,'ro','MarkerSize',2);
plot(timswin,xmaxVGSIntmp,'k*','MarkerSize',2);
axis([3*timbig/2 2*timbig -5 5]);
subplot(4,1,4,'align'); 
hold on
hrf = plotreflinesr(gca,detects,'x','k'); 
hrf = plotreflinesr(gca,Detects,'x','r'); 
hrf = plotreflinesr(gca,detects2_8,'x','g'); 
hrf = plotreflinesr(gca,Detects2_8,'x','b'); 
plot(timswin,0.5+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn,'k');
plot(timswin,0.1*abs(loopoff),'go','MarkerSize',2);
plot(timsPGC,min(0,real(PGCrot)),'r');
plot(timsPGC,min(0,real(SILBrot)),'b');
plot(timsPGC,min(0,real(VGZrot)),'k');
axis([3*timbig/2 2*timbig -1 1]);
% orient landscape
% print('-dpdf',[JDAY,'_',int2str(timlook),'_2.pdf'])

figure
colormap(jet)
scatter(xmaxPGVGn-xmaxPGSIn+xmaxVGSIn,xcmaxAVEn,3,AmpComp)
axis([-5 5 -0.2 1.0])
colorbar

%cputime-t;
tot=cputime-tt
