%Sits at one spot (rotations; delays) and looks for the cross-correlation.
clear all
t=cputime; tt=cputime;
close all
format short e
%Rotations & offsets:
    timoffrot=[2005 255 05 19 48 +085 -003  080  115  090]; %Starts today.  Slow migration?
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
vectoff=[ArmCat(:,5)-SSIBsoff ArmCat(:,6)-SILBsoff];
distoff=vectoff.*vectoff;
distoff2=sum(distoff,2);
m=0;
Detects=0;
for n=1:length(detects)
     if distoff2(n)<=9
%      if isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff+1]) 
         m=m+1;
         Detects(m)=detects(n);
     end
end
ArmCat=load(['/data2/arubin/CNDC/',YEAR,'/list.',YEAR,'.pgsssiMAJ_new']);
detects2_8=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
vectoff=[ArmCat(:,5)-SSIBsoff ArmCat(:,6)-SILBsoff];
distoff=vectoff.*vectoff;
distoff2=sum(distoff,2);
m=0;
Detects2_8=0;
for n=1:length(detects2_8)
    if distoff2(n)<=9
%      if isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff+1]) 
         m=m+1;
         Detects2_8(m)=detects2_8(n);
     end
end

%Read data:
direc=[YEAR,'/SEPT/'];
prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];
PGCEdat=[prename,'.PGC..BHE.D.SAC'];
PGCNdat=[prename,'.PGC..BHN.D.SAC'];
%PGCZdat=[prename,'.PGC..BHZ.D.SAC'];
SSIBEdat=[prename,'.SSIB..HHE.D.SAC'];
SSIBNdat=[prename,'.SSIB..HHN.D.SAC'];
SILBEdat=[prename,'.KLNB..HHE.D.SAC'];
SILBNdat=[prename,'.KLNB..HHN.D.SAC'];

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
%     % Truncate for quicker calculation.  100sps stations are longer for
%     % later interpolation
%     timwin=4*60+1;
%     timstart=max(0, timlook-timwin); timend=min(86400, timlook+timwin);
%     PGCE=PGCE(timstart*40+1:timend*40)-mean(PGCE(timstart*40+1:timend*40));
%     PGCN=PGCN(timstart*40+1:timend*40)-mean(PGCN(timstart*40+1:timend*40)); 
%     %PGCZ=PGCZ(timstart*40+1:timend*40);
%     SSIBE=SSIBE(timstart*100:timend*100+1)-mean(SSIBE(timstart*100:timend*100+1));
%     SSIBN=SSIBN(timstart*100:timend*100+1)-mean(SSIBN(timstart*100:timend*100+1));
%     SILBE=SILBE(timstart*100:timend*100+1)-mean(SILBE(timstart*100:timend*100+1));
%     SILBN=SILBN(timstart*100:timend*100+1)-mean(SILBN(timstart*100:timend*100+1));
%     timsPGC=timsPGC(timstart*40+1:timend*40);
%     timsSSIB=timsSSIB(timstart*100:timend*100+1);
%     timsSILB=timsSILB(timstart*100:timend*100+1);
% %shrten=cputime-t
% %t=cputime;
% %
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
npo=2;
npa=1;
[PGCEf]=1.6e-4*bandpass(PGCE,40,lo,hi,npo,npa,'butter');
[PGCNf]=1.6e-4*bandpass(PGCN,40,lo,hi,npo,npa,'butter');
%[PGCZf]=1.6e-4*bandpass(PGCZ,40,lo,hi,npo,npa,'butter');
[SSIBEf]=4.0e-3*bandpass(SSIBE,100,lo,hi,npo,npa,'butter'); %times 5 for early 2003
[SSIBNf]=4.0e-3*bandpass(SSIBN,100,lo,hi,npo,npa,'butter'); 
[SILBEf]=4.0e-3*bandpass(SILBE,100,lo,hi,npo,npa,'butter');
[SILBNf]=4.0e-3*bandpass(SILBN,100,lo,hi,npo,npa,'butter');
fltr=cputime-t
t=cputime;

%Decimate the 100 sps data:
% SSIBEfd = interp1(timsSSIB,SSIBEf,timsPGC,'linear')';
% SSIBNfd = interp1(timsSSIB,SSIBNf,timsPGC,'linear')';
% SILBEfd = interp1(timsSILB,SILBEf,timsPGC,'linear')';
% SILBNfd = interp1(timsSILB,SILBNf,timsPGC,'linear')';
% dmate=cputime-t
% t=cputime;
SSIBEfd = resample(SSIBEf,2,5);
SSIBNfd = resample(SSIBNf,2,5);
SILBEfd = resample(SILBEf,2,5);
SILBNfd = resample(SILBNf,2,5);
dmate=cputime-t
t=cputime;

PGC=PGCEf+1i*PGCNf;
SSIB=SSIBEfd+1i*SSIBNfd;
SILB=SILBEfd+1i*SILBNfd;
PGCrot=PGC*exp(1i*rotPGC);
SSIBrot=SSIB*exp(1i*rotSSIB);
SILBrot=SILB*exp(1i*rotSILB);
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
    %tracesht=cputime-t
    %t=cputime;

realPGC=real(PGCrot);
realSSIB=real(SSIBrot);
realSILB=real(SILBrot);
% realPGC=imag(PGCrot); %A quick/dirty way to try the same on the orthogonal components
% realSSIB=imag(SSIBrot);
% realSILB=imag(SILBrot);

plotstart=1; %(timlook-timstart)*40
plotend=tracelenPGC; %(timlook-timstart)*40
ltmax=max([real(PGCrot(plotstart:plotend)); ...
          real(SSIBrot(plotstart:plotend)); ...
          real(SILBrot(plotstart:plotend))]);
ltmax=min(ltmax,0.75);
ltmin=min([real(PGCrot(plotstart:plotend)); ...
          real(SSIBrot(plotstart:plotend)); ...
          real(SILBrot(plotstart:plotend))]);
ltmin=max(ltmin,-0.75);

figure %Superposed traces, 150s window:
subplot(4,1,1); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([0 timwin/2 ltmin ltmax]);
%xlim([timlook-75 timlook-37.5]);
subplot(4,1,2); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([timwin/2 timwin ltmin ltmax]);
%xlim([timlook-37.5 timlook]);
subplot(4,1,3); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([timwin 3*timwin/2 ltmin ltmax]);
%xlim([timlook timlook+37.5]);
subplot(4,1,4); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([3*timwin/2 2*timwin ltmin ltmax]);
%xlim([timlook+37.5 timlook+75]);
traceplt=cputime-t
t=cputime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Autocorrelation of stations.  Those that end in "2" are the running
%   cumulative sum, to be used later by differncing the window edpoints.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PGCauto=realPGC.*realPGC;
PGC2=cumsum(PGCauto);
SSIBauto=realSSIB.*realSSIB;
SSIB2=cumsum(SSIBauto);
SILBauto=realSILB.*realSILB;
SILB2=cumsum(SILBauto);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Cross-correlation between stations, with small offsets up to +/- mshift.
%  First index is pointwise multiplication of traces; second is shifting offset.
%  lenx is shorter than tracelenPGC by mshift at each end (?)
%  For PGSS and PGSI, SSI and SIL are shifted relative to PGC, by 1 each time through loop.
%  For SISS, SSI is shifted relative to SILB.
%  cumsumPGSS etc. are the running cumulative sum of the x-correlation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mshift=11;
lenx=tracelenPGC-2*mshift;
PGSSx=zeros(lenx, 2*mshift+1);
PGSIx=zeros(lenx, 2*mshift+1);
SISSx=zeros(lenx, 2*mshift+1);
for n=-mshift:mshift;
    PGSSx(:,n+mshift+1)=realPGC(1+mshift:tracelenPGC-mshift).* ...
        realSSIB(1+mshift-n:tracelenPGC-mshift-n);
    PGSIx(:,n+mshift+1)=realPGC(1+mshift:tracelenPGC-mshift).* ...
        realSILB(1+mshift-n:tracelenPGC-mshift-n);
    SISSx(:,n+mshift+1)=realSILB(1+mshift:tracelenPGC-mshift).* ...
        realSSIB(1+mshift-n:tracelenPGC-mshift-n);
end
cumsumPGSS=cumsum(PGSSx);
cumsumPGSI=cumsum(PGSIx);
cumsumSISS=cumsum(SISSx);
xcshifts=cputime-t
t=cputime;

winbig=2*(tracelenPGC/2-120);
timbig=winbig/(2*40);
igstart=floor(tracelenPGC/2-winbig/2)+1;
winlen=4*40;
winoff=40/1;
%winoff=40/2;
nwin=floor((winbig-winlen)/winoff);
timswin=zeros(nwin,1);
sumsPGSS=zeros(nwin,2*mshift+1);
sumsPGSI=zeros(nwin,2*mshift+1);
sumsSISS=zeros(nwin,2*mshift+1);
sumsPGC2=zeros(nwin,2*mshift+1);
sumsSSIB2=zeros(nwin,2*mshift+1);
sumsSILB2=zeros(nwin,2*mshift+1);
sumsSILB2b=zeros(nwin,2*mshift+1);
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
end
denomPGSSn=realsqrt(sumsPGC2.*sumsSSIB2);
denomPGSIn=realsqrt(sumsPGC2.*sumsSILB2);
denomSISSn=realsqrt(sumsSILB2b.*sumsSSIB2);
sumsPGSSn=sumsPGSS./denomPGSSn;
%sumsPGSSn'
sumsPGSIn=sumsPGSI./denomPGSIn;
sumsSISSn=sumsSISS./denomSISSn;
[xcmaxPGSSn,imaxPGSS]=max(sumsPGSSn,[],2); %Integer-offset max cross-correlation
[xcmaxPGSIn,imaxPGSI]=max(sumsPGSIn,[],2);
[xcmaxSISSn,imaxSISS]=max(sumsSISSn,[],2);
%Parabolic fit:
% [xmaxPGSSn,ymaxPGSSn,xloPGSSn,xhiPGSSn]=parabol(nwin,mshift,sumsPGSSn,imaxPGSS); %These return some poor measure of "confidence".
% [xmaxPGSIn,ymaxPGSIn,xloPGSIn,xhiPGSIn]=parabol(nwin,mshift,sumsPGSIn,imaxPGSI);
% [xmaxSISSn,ymaxSISSn,xloSISSn,xhiSISSn]=parabol(nwin,mshift,sumsSISSn,imaxSISS);
[xmaxPGSSn,ymaxPGSSn]=parabol(nwin,mshift,sumsPGSSn,imaxPGSS); %Interpolated max cross-correlation
[xmaxPGSIn,ymaxPGSIn]=parabol(nwin,mshift,sumsPGSIn,imaxPGSI);
[xmaxSISSn,ymaxSISSn]=parabol(nwin,mshift,sumsSISSn,imaxSISS);

ix=sub2ind(size(denomPGSSn),(1:nwin)',imaxPGSS);
ampPGSS=sqrt(denomPGSSn(ix)); %This makes amplitude linear rather than quadratic with counts.
ampPGC2=sumsPGC2(ix); %by construction PGC2 is the same for all shifts
ampSSIB2=sumsSSIB2(ix);
ix=sub2ind(size(denomPGSIn),(1:nwin)',imaxPGSI);
ampPGSI=sqrt(denomPGSIn(ix));
ampSILB2=sumsSILB2(ix);
ix=sub2ind(size(denomSISSn),(1:nwin)',imaxSISS);
ampSISS=sqrt(denomSISSn(ix));
AmpComp(1:4)=0;
AmpComp(5:nwin)=((ampPGC2(5:nwin)+ampSSIB2(5:nwin)+ampSILB2(5:nwin))- ...
                (ampPGC2(1:nwin-4)+ampSSIB2(1:nwin-4)+ampSILB2(1:nwin-4)))./ ...
                ((ampPGC2(5:nwin)+ampSSIB2(5:nwin)+ampSILB2(5:nwin))+ ...
                (ampPGC2(1:nwin-4)+ampSSIB2(1:nwin-4)+ampSILB2(1:nwin-4))) ;
%Center them
imaxPGSS=imaxPGSS-mshift-1;
imaxPGSI=imaxPGSI-mshift-1;
imaxSISS=imaxSISS-mshift-1;
iloopoff=imaxPGSI-imaxPGSS+imaxSISS; %How well does the integer loop close?
xmaxPGSSn=xmaxPGSSn-mshift-1;
xmaxPGSIn=xmaxPGSIn-mshift-1;
xmaxSISSn=xmaxSISSn-mshift-1;
loopoff=xmaxPGSIn-xmaxPGSSn+xmaxSISSn; %How well does the interpolated loop close?
xcmaxAVEn=(xcmaxPGSSn+xcmaxPGSIn+xcmaxSISSn)/3;
loopoffmax=1.5
xcmaxAVEnmin=0.4
xcnshifts=cputime-t
t=cputime;

ampmax=max([ampPGSS; ampPGSI; ampSISS]);
figure 
subplot(4,1,1,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
%plot(timswin,3.5+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*8,'g');
plot(timswin,xmaxPGSSn,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIn,'ro','MarkerSize',2);
plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',2);
%plot(timswin,xmaxSISSn,'k*','MarkerSize',2);
axis([0 timbig/2 -8 8]);
ylabel('samples')
box on
subplot(4,1,2,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
%plot(timswin,3.5+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*8,'g');
plot(timswin,xmaxPGSSn,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIn,'ro','MarkerSize',2);
plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',2);
%plot(timswin,xmaxSISSn,'k*','MarkerSize',2);
axis([timbig/2 timbig -8 8]);
ylabel('samples')
box on
subplot(4,1,3,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
%plot(timswin,3.5+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*8,'g');
plot(timswin,xmaxPGSSn,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIn,'ro','MarkerSize',2);
plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',2);
%plot(timswin,xmaxSISSn,'k*','MarkerSize',2);
axis([timbig 3*timbig/2 -8 8]);
ylabel('samples')
box on
subplot(4,1,4,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
%plot(timswin,3.5+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*8,'g');
plot(timswin,xmaxPGSSn,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIn,'ro','MarkerSize',2);
plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',2);
%plot(timswin,xmaxSISSn,'k*','MarkerSize',2);
axis([3*timbig/2 2*timbig -8 8]);
xlabel('sec')
ylabel('samples')
box on
title([IDENTIF,'_{',int2str(lo),'-',int2str(hi),'}'])
orient landscape
print('-depsc',['FIGS/',IDENTIF,'_',int2str(lo),'-',int2str(hi),'a.eps'])

medxcmaxAVEn=median(xcmaxAVEn)
xmaxPGSSntmp=xmaxPGSSn;
xmaxPGSIntmp=xmaxPGSIn;
xmaxSISSntmp=xmaxSISSn;
%loopoff=xmaxPGSIn-xmaxPGSSn+xmaxSISSn; 
iup=4;
nin=0;
for n=1:nwin;
    if xcmaxAVEn(n)<xcmaxAVEnmin || abs(xmaxPGSIn(n)-xmaxPGSSn(n)+xmaxSISSn(n))>loopoffmax
        xmaxPGSSntmp(n)=20; xmaxPGSIntmp(n)=20; xmaxSISSntmp(n)=20; 
    else
        interpPGSSn=interp(sumsPGSSn(n,:),iup,3);
        interpPGSIn=interp(sumsPGSIn(n,:),iup,3);
        interpSISSn=interp(sumsSISSn(n,:),iup,3);
        [xcmaxinterpPGSSn,imaxinterpPGSS]=max(interpPGSSn);
        [xcmaxinterpPGSIn,imaxinterpPGSI]=max(interpPGSIn);
        [xcmaxinterpSISSn,imaxinterpSISS]=max(interpSISSn);
        xcmaxconprev=0.;
        for iPGSS=max(1,imaxinterpPGSS-3*iup):min(imaxinterpPGSS+3*iup,iup*(2*mshift+1)-(iup-1)) %3 samples from peak; 
                                                                                            %intentionally wider than acceptable;
                                                                                            %iup-1 are extrapolated points
            for iPGSI=max(1,imaxinterpPGSI-3*iup):min(imaxinterpPGSI+3*iup,iup*(2*mshift+1)-(iup-1))
                ibangon = (iup*mshift+1)-iPGSI+iPGSS;
                if ibangon >= 1 && ibangon<=iup*(2*mshift+1)
                    xcmaxcon=interpPGSSn(iPGSS)+interpPGSIn(iPGSI)+interpSISSn(ibangon);
                    if xcmaxcon > xcmaxconprev
                        xcmaxconprev=xcmaxcon;
                        iPGSSbang=iPGSS;
                        iPGSIbang=iPGSI;
                    end
                end
            end
        end
        iSISSbang=(iup*mshift+1)-iPGSIbang+iPGSSbang;
        if abs(iPGSSbang-imaxinterpPGSS) <= loopoffmax*iup && ...
           abs(iPGSIbang-imaxinterpPGSI) <= loopoffmax*iup && ...
           abs(iSISSbang-imaxinterpSISS) <= loopoffmax*iup && ...
           interpPGSSn(iPGSSbang)+interpPGSIn(iPGSIbang)+interpSISSn(iSISSbang) >= 3*xcmaxAVEnmin
            xmaxPGSSntmp(n)=(iPGSSbang-(iup*mshift+1))/iup;
            xmaxPGSIntmp(n)=(iPGSIbang-(iup*mshift+1))/iup;
            xmaxSISSntmp(n)=(iSISSbang-(iup*mshift+1))/iup;
%             interpPGSSn(iPGSSbang)
%             interpPGSIn(iPGSIbang)
%             interpSISSn(iSISSbang)
%             xcmaxcon
%             xcmaxPGSSn(n)
%             xcmaxPGSIn(n)
%             xcmaxSISSn(n)
            
            xmaxPGSSntmp(n)
            xmaxPGSIntmp(n)
            xmaxSISSntmp(n)
            xmaxPGSSn(n)
            xmaxPGSIn(n)
            xmaxSISSn(n)
            xcmaxAVEn(n)
            nin=nin+1;
        else
            xmaxPGSSntmp(n)=20; xmaxPGSIntmp(n)=20; xmaxSISSntmp(n)=20; 
        end
%         if nin>=3500
%             figure
%             plot(1:23,sumsPGSSn(n,:),'b*')
%             hold on
%             plot(0.3333*(1:iup*23)+0.6667,interpPGSSn,'bo')
%             plot(1:23,sumsPGSIn(n,:),'r*')
%             plot(0.3333*(1:iup*23)+0.6667,interpPGSIn,'ro')
%             plot(1:23,sumsSISSn(n,:),'g*')
%             plot(0.3333*(1:iup*23)+0.6667,interpSISSn,'go')
%         end
%         pause
    end
end
figure 
subplot(4,1,1,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*8,'g');
plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
axis([0 timbig/2 -8 8]);
ylabel('samples')
box on
subplot(4,1,2,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*8,'g');
plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
axis([timbig/2 timbig -8 8]);
ylabel('samples')
box on
subplot(4,1,3,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*8,'g');
plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
axis([timbig 3*timbig/2 -8 8]);
ylabel('samples')
box on
subplot(4,1,4,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*8,'g');
plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
axis([3*timbig/2 2*timbig -8 8]);
xlabel('sec')
ylabel('samples')
box on
% title([IDENTIF,'_{',int2str(lo),'-',int2str(hi),'}'])
% orient landscape
% print('-depsc',['FIGS/',IDENTIF,'_',int2str(lo),'-',int2str(hi),'b.eps'])

for n=1:nwin;
    if abs(loopoff(n))>2.0
        loopoff(n)=30; 
    end
end
figure 
subplot(4,1,1,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
axis([0 timbig/2 -8 8]);
box on
subplot(4,1,2,'align'); 
hold on
hrf = plotreflinesr(gca,detects,'x','k'); 
hrf = plotreflinesr(gca,Detects,'x','r'); 
hrf = plotreflinesr(gca,detects2_8,'x','g'); 
hrf = plotreflinesr(gca,Detects2_8,'x','b'); 
plot(timswin,0.5+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn,'k');
plot(timswin,-1+0.1*abs(loopoff),'g.','MarkerSize',1);
plot(timsPGC,min(0,real(PGCrot)),'r');
plot(timsPGC,min(0,real(SSIBrot)),'b');
plot(timsPGC,min(0,real(SILBrot)),'k');
axis([0 timbig/2 -1 1]);
box on
subplot(4,1,3,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
axis([timbig/2 timbig -8 8]);
box on
subplot(4,1,4,'align'); 
hold on
hrf = plotreflinesr(gca,detects,'x','k'); 
hrf = plotreflinesr(gca,Detects,'x','r'); 
hrf = plotreflinesr(gca,detects2_8,'x','g'); 
hrf = plotreflinesr(gca,Detects2_8,'x','b'); 
plot(timswin,0.5+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn,'k');
plot(timswin,-1+0.1*abs(loopoff),'g.','MarkerSize',1);
plot(timsPGC,min(0,real(PGCrot)),'r');
plot(timsPGC,min(0,real(SSIBrot)),'b');
plot(timsPGC,min(0,real(SILBrot)),'k');
axis([timbig/2 timbig -1 1]);
box on
% title([IDENTIF,'_{',int2str(lo),'-',int2str(hi),'}'])
% orient landscape
% print('-depsc',['FIGS/',IDENTIF,'_',int2str(lo),'-',int2str(hi),'c.eps'])

figure
subplot(4,1,1,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
axis([timbig 3*timbig/2 -8 8]);
box on
subplot(4,1,2,'align'); 
hold on
hrf = plotreflinesr(gca,detects,'x','k'); 
hrf = plotreflinesr(gca,Detects,'x','r'); 
hrf = plotreflinesr(gca,detects2_8,'x','g'); 
hrf = plotreflinesr(gca,Detects2_8,'x','b'); 
plot(timswin,0.5+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn,'k');
plot(timswin,-1+0.1*abs(loopoff),'g.','MarkerSize',1);
plot(timsPGC,min(0,real(PGCrot)),'r');
plot(timsPGC,min(0,real(SSIBrot)),'b');
plot(timsPGC,min(0,real(SILBrot)),'k');
axis([timbig 3*timbig/2 -1 1]);
box on
subplot(4,1,3,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
axis([3*timbig/2 2*timbig -8 8]);
box on
subplot(4,1,4,'align'); 
hold on
hrf = plotreflinesr(gca,detects,'x','k'); 
hrf = plotreflinesr(gca,Detects,'x','r'); 
hrf = plotreflinesr(gca,detects2_8,'x','g'); 
hrf = plotreflinesr(gca,Detects2_8,'x','b'); 
plot(timswin,0.5+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn,'k');
plot(timswin,-1+0.1*abs(loopoff),'g.','MarkerSize',1);
plot(timsPGC,min(0,real(PGCrot)),'r');
plot(timsPGC,min(0,real(SSIBrot)),'b');
plot(timsPGC,min(0,real(SILBrot)),'k');
axis([3*timbig/2 2*timbig -1 1]);
box on
% title([IDENTIF,'_{',int2str(lo),'-',int2str(hi),'}'])
% orient landscape
% print('-depsc',['FIGS/',IDENTIF,'_',int2str(lo),'-',int2str(hi),'d.eps'])

figure
colormap(jet)
scatter(xmaxPGSIn-xmaxPGSSn+xmaxSISSn,xcmaxAVEn,3,AmpComp)
axis([-5 5 -0.2 1.0])
colorbar

%cputime-t;
tot=cputime-tt
