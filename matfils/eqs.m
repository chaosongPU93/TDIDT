close all
prename='2003.285.03.42.12'
%prename='2005.020.01.04.30';
%prename='2005.152.23.25.54'
%prename='2005.334.22.57.48'
hi=5.;
lo=1.;
datname=[prename,'.0000.CN.PGC..BHE.D.SAC'];
[SeisData,HdrData,tnu,pobj,tims]=readsac(datname,0,'l');
tracelen=length(SeisData);
if(rem(tracelen,2)==0)
    SeisData(tracelen)=[];
    tims(tracelen)=[];
    tracelen=tracelen-1;
end
%Looks like filtering the 40-Hz data followed by interpolation introduces 
%less high-frequency stuff that does not appear in the original filtered
%40-Hz data:
%PGCE=spline(tims,SeisData,timsx);
%[PGCE]=bandpass(PGCE,100,1.0,5.0,2,1,'butter');
[PGCE]=1.6e-4*bandpass(SeisData,40,lo,hi,2,1,'butter');
timsx=0:0.01:0.025*(tracelen-1);
PGCE=spline(tims,PGCE,timsx)'; %transpose to make it match SILB?

datname=[prename,'.0000.CN.PGC..BHN.D.SAC'];
[SeisData,HdrData,tnu,pobj,tims]=readsac(datname,0,'l');
tracelen=length(SeisData);
if(rem(tracelen,2)==0)
    SeisData(tracelen)=[];
    tims(tracelen)=[];
    tracelen=tracelen-1;
end
[PGCN]=1.6e-4*bandpass(SeisData,40,lo,hi,2,1,'butter');
timsx=0:0.01:0.025*(tracelen-1);
PGCN=spline(tims,PGCN,timsx)'; %transpose to make it match SILB
PGC=PGCE+1i*PGCN;
plot(timsx,real(PGC));
hold on
plot(timsx,imag(PGC),'r-');
xlim([37 43]);

%scalef=10.;
datname=[prename,'.0000.CN.SILB..HHE.D.SAC'];
[SeisData,HdrData,tnu,pobj,tims]=readsac(datname,0,'l');
[SILBE]=1.4*4.0e-3*bandpass(SeisData,100,lo,hi,2,1,'butter');
datname=[prename,'.0000.CN.SILB..HHN.D.SAC'];
[SeisData,HdrData,tnu,pobj,tims]=readsac(datname,0,'l');
[SILBN]=1.4*4.0e-3*bandpass(SeisData,100,lo,hi,2,1,'butter');
SILB=SILBE+1i*SILBN;
figure
plot(tims,real(SILB));
hold on
plot(tims,imag(SILB),'r-');
xlim([37 43]);

%SILBE(1:tracelen-110)=PGCE(111:tracelen);  %to test

mlag=250;
PGSIx=xcorr(real(PGC),real(SILB),mlag);
figure
plot(PGSIx,'r-x');
[~,imaxx]=max(PGSIx)
imaxx=imaxx-(mlag+1) %This makes a logical lag of zero equal to zero
[~,imaxPGC]=max(abs(real(PGC)));

%set longer and shorter windows:
ipre=imaxPGC-300; ipost=imaxPGC+600; ishorter=100;
PGCshort=PGC(ipre:ipost);
timsPGCshort=timsx(ipre:ipost);
lenPGs=length(PGCshort);

%taper SILB
SILBshort=SILB(ipre-imaxx:ipost-imaxx);
timsSILBshort=tims(ipre-imaxx:ipost-imaxx)+0.01*imaxx;
x=(pi/200:pi/200:pi/2)';
SILBshort(1:ishorter)=0;
SILBshort(ishorter+1:200)=sin(x).*SILBshort(ishorter+1:200);
lenSIs=length(SILBshort);
x=flipud(x);
SILBshort(lenSIs-(ishorter-1):lenSIs)=0;
SILBshort(lenSIs-(ishorter-1)-100:lenSIs-ishorter)=...
    sin(x(1:100)).*SILBshort(lenSIs-(ishorter-1)-100:lenSIs-ishorter);

%Rotations
%Rotating the vector by -alpha is the same as rotating the axes by alpha.
angs=0:-pi/36:-2*pi+pi/36;
rots=exp(1i*angs);
numrots=length(rots);
PGCrot=PGCshort*rots(1:numrots/2);
SILBrot=SILBshort*rots;
% figure
% plot(real(PGCrot(:,numrots/4+1)),imag(PGCshort),'r-')
mlag=40;
maxxc=zeros(numrots,numrots/2);
for Prot=1:numrots/2
    for Srot=1:numrots
        xc=xcorr(real(PGCrot(:,Prot)),real(SILBrot(:,Srot)),mlag);
        maxxc(Srot,Prot)=max(xc);
        minxc(Srot,Prot)=min(xc);
    end
end
[SImax,iSImax]=max(maxxc);
[PGmax,iPGmax]=max(SImax);
[SImin,iSImin]=min(maxxc);
[PGmin,iPGmin]=min(SImax);
PGCrotfin=real(PGCrot(:,iPGmax));
SILBrotfin=real(SILBrot(:,iSImax(iPGmax)));
PGCrotfini=imag(PGCrot(:,iPGmax));
SILBrotfini=imag(SILBrot(:,iSImax(iPGmax)));
[iPGmax iSImax(iPGmax)]

PGSIx=xcorr(PGCrotfin,SILBrotfin,mlag);
[~,imaxxupdat]=max(PGSIx)
imaxxupdat=imaxxupdat-(mlag+1) %This sets a logical lag of 0 to 0 (new win)
timsSILBshort=timsSILBshort+0.01*imaxxupdat;
figure
plot(PGSIx,'r-x')

figure
plot(timsPGCshort,PGCrotfin,'b-')
hold on
plot(timsSILBshort,SILBrotfin,'r-')
figure
plot(timsPGCshort,PGCrotfini,'k-')
hold on
plot(timsSILBshort,SILBrotfini,'r-')
PGSI_x=PGCrotfin(ishorter:lenPGs-ishorter).*...
    SILBrotfin(ishorter-imaxxupdat:lenPGs-ishorter-imaxxupdat);
PGSI_xi=PGCrotfini(ishorter:lenPGs-ishorter).*...
    SILBrotfini(ishorter-imaxxupdat:lenPGs-ishorter-imaxxupdat);
figure
plot(timsPGCshort(ishorter:lenPGs-ishorter),PGSI_x)
hold on 
plot(timsPGCshort(ishorter:lenPGs-ishorter),PGSI_xi,'r-')
PGSI_xcum=cumsum(PGSI_x);
figure
plot(timsPGCshort(ishorter:lenPGs-ishorter),PGSI_xcum);
hold on
PGSI_xicum=cumsum(PGSI_xi);
plot(timsPGCshort(ishorter:lenPGs-ishorter),PGSI_xicum,'r-');

figure
subplot(2,1,1)
plot(timsSILBshort,SILBrotfin,'r-')
hold on
plot(timsPGCshort,PGCrotfin,'b-')
text(40, 140.,'M2.4 earthquake, 45 km depth, 1-5 Hz','FontSize',12)
text(40, 90.,'Rotated & shifted to maximize coherence','FontSize',12)
text(37.4, 40.,'SILB','Color',[1 0 0],'FontSize',12)
text(37.4, -60.,'PGC','Color',[0 0 1],'FontSize',12)
XLABEL('sec')
YLABEL('Velocity counts')
axis([37 44 -200 200]);
subplot(2,1,2)
plot(timsPGCshort(ishorter:lenPGs-ishorter),PGSI_x)
hold on 
plot(timsPGCshort(ishorter:lenPGs-ishorter),PGSI_xi,'r-')
text(40, 7700.,'Pointwise cross-correlation, PGC & SILB','FontSize',12)
text(40.8, 1500.,'Maximally coherent component','Color',[0 0 1],'FontSize',12)
text(40.8, -2300.,'Orthogonal component','Color',[1 0 0],'FontSize',12)
%PGSI_xcum=cumsum(PGSI_x);
axis([37 44 -5000 10000]);
XLABEL('sec')
%orient landscape
print('-depsc',[prename,'.eps'])

%plot(timsPGCshort(ishorter:lenPGs-ishorter),PGCrotfin,'b-')
% hold on
% plot(SILBrotfin,imag(SILBshort),'r-')

% mlag=20;
% PGSIx=xcorr(real(PGCshort),real(SILBshort),mlag);
% figure
% plot(PGSIx,'r-x');
% [~,imaxxupdat]=max(PGSIx)
% imaxxupdat=imaxxupdat-(mlag+1) %This sets a logical lag of 0 to 0 (new win)
% PGSI_Ex=real(PGCshort(ishorter:lenPGs-ishorter)).*...
%     real(SILBshort(ishorter-imaxxupdat:lenPGs-ishorter-imaxxupdat));
% figure
% plot(timsPGCshort(ishorter:lenPGs-ishorter),PGSI_Ex);
% PGSI_Excum=cumsum(PGSI_Ex);
% figure
% plot(timsPGCshort(ishorter:lenPGs-ishorter),PGSI_Excum);
% 
