%Sits at one spot (rotations; delays) for one day and looks for the cross-correlation.
%Now also plots aligned 4-s windows.  As of July 2012.
format short e
clear all
set(0,'DefaultFigureVisible','on');

   %timoffrot(1,:)=[2003 062 07 42 12 +086 +020  80 115  50];
   timoffrot(1,:)=[2003 063 07 42 12 +086 +020  80 115  50];
   %timoffrot(3,:)=[2003 064 07 42 12 +086 +020  80 115  50];
   %timoffrot(4,:)=[2003 065 07 42 12 +086 +020  80 115  50];
%    
%    timoffrot(1,:)=[2005 254 07 42 12 +086 +020  80 115  50];
%    timoffrot(2,:)=[2005 255 07 42 12 +086 +020  80 115  50];
%    timoffrot(3,:)=[2005 256 07 42 12 +086 +020  80 115  50];
%    timoffrot(4,:)=[2005 257 07 42 12 +086 +020  80 115  50];
%   
%    timoffrot(1,:)=[2004 196 07 42 12 +086 +020  80 115  50];
%    timoffrot(2,:)=[2004 197 07 42 12 +086 +020  80 115  50];
%    timoffrot(3,:)=[2004 198 07 42 12 +086 +020  80 115  50];
%    timoffrot(4,:)=[2004 199 07 42 12 +086 +020  80 115  50];

nd=length(timoffrot(:,1))
close all

%mshift=12; %16 %11 %14 %12_for_+085+020 %15_for_+059-073  %19 for 128-s?
hi=6.;
lo=1.5;
% hi=12.;
% lo=1;
% hi=8.;
% lo=2.;
% hi=6.;
% lo=0.75;

year=timoffrot(nd,1);
YEAR=int2str(year);
jday=timoffrot(nd,2);
if jday <= 9
    JDAY=['00',int2str(jday)];
elseif jday<= 99
    JDAY=['0',int2str(jday)];
else
    JDAY=int2str(jday)
end
MO=day2month(jday,year)
timlook=3600*timoffrot(nd,3)+60*timoffrot(nd,4)+timoffrot(nd,5)
rotPGC=pi*(timoffrot(nd,8)-90)/180;
rotSSIB=pi*(timoffrot(nd,9)-90)/180;
rotSILB=pi*(timoffrot(nd,10)-90)/180;
SSIBsoff=timoffrot(nd,6);
SILBsoff=timoffrot(nd,7);
SSIBtoff=SSIBsoff/40;
SILBtoff=SILBsoff/40;
IDENTIF=[YEAR,'.',JDAY,'.',int2str(timoffrot(nd,6)),'.',int2str(timoffrot(nd,7)), ...
    '.',int2str(timoffrot(nd,8)),'.',int2str(timoffrot(nd,9)),'.',int2str(timoffrot(nd,10))]

% %Read Armbruster's detections:
% ArmCat=load(['/data2/arubin/CNDC/',YEAR,'/list.',YEAR,'.pgsssiMAJ']);
% detects=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
% vectoff=[ArmCat(:,5)-SSIBsoff ArmCat(:,6)-SILBsoff (ArmCat(:,6)-SILBsoff)-(ArmCat(:,5)-SSIBsoff)];
% distoff=vectoff.*vectoff;
% distoff2=sum(distoff,2);
% gridoff=max(abs(vectoff),[],2); %Try this instead.  Should change w/mshift (defined later, though)
% m0=0;
% Detects=0;
% for n=1:length(detects)
%     %if distoff2(n)<=9
%     if gridoff(n)<=mshift-1
%          m0=m0+1;
%          Detects(m0)=detects(n);
%      end
% end
% ArmCat=load(['/data2/arubin/CNDC/',YEAR,'/list.',YEAR,'.pgsssiMAJ_new']);
% detects2_8=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
% vectoff=[ArmCat(:,5)-SSIBsoff ArmCat(:,6)-SILBsoff];
% distoff=vectoff.*vectoff;
% distoff2=sum(distoff,2);
% gridoff=max(abs(vectoff),[],2); %Try this instead.  Should change w/mshift (defined later, though)
% m1=0;
% Detects2_8=0;
% for n=1:length(detects2_8)
%     %if distoff2(n)<=9
%     if gridoff(n)<=mshift-1
%          m1=m1+1;
%          Detects2_8(m1)=detects2_8(n);
%      end
% end

%Read data:
direc=[YEAR,'/',MO,'/'];
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

winbig=2*(tracelenPGC/2-120);
timbig=winbig/(2*40);
igstart=floor(tracelenPGC/2-winbig/2)+1;

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
if year==2003 && jday<213
    [SSIBEf]=20.0e-3*bandpass(SSIBE,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
    [SSIBNf]=20.0e-3*bandpass(SSIBN,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
    [SILBEf]=20.0e-3*bandpass(SILBE,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0* 
    [SILBNf]=20.0e-3*bandpass(SILBN,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
else
    [SSIBEf]=4.0e-3*bandpass(SSIBE,100,lo,hi,npo,npa,'butter'); 
    [SSIBNf]=4.0e-3*bandpass(SSIBN,100,lo,hi,npo,npa,'butter'); 
    [SILBEf]=4.0e-3*bandpass(SILBE,100,lo,hi,npo,npa,'butter'); 
    [SILBNf]=4.0e-3*bandpass(SILBN,100,lo,hi,npo,npa,'butter'); 
end


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

% % Not efficient but easy book-keeping:
% SSIBauto=realSSIB.*realSSIB;
% SSIB2=cumsum(SSIBauto);
% SILBauto=realSILB.*realSILB;
% SILB2=cumsum(SILBauto);

shortlen=30*40;
scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite;

mshift=5; %max shift of SSIB relative to SILB
wins=[21800 22580];
for n=1:size(wins,1)
    plotstart=wins(n,1)*40;
    plotend=wins(n,2)*40;
    longlen=plotend-plotstart+1;
    SILBauto=realSILB(plotstart:plotend).*realSILB(plotstart:plotend);
    SILB2=cumsum(SILBauto);
    SSIBauto=realSSIB(plotstart-mshift:plotend+mshift).*realSSIB(plotstart-mshift:plotend+mshift); %Longer b.c. of mshift.
    SSIB2=cumsum(SSIBauto);
    SISSx=zeros(longlen,2*mshift+1);
    for k=-mshift:mshift;
        SISSx(:,k+mshift+1)=realSILB(plotstart:plotend).*realSSIB(plotstart+k:plotend+k);
    end
    cumsumSISS=cumsum(SISSx);
    SISSnum=cumsumSISS(21:longlen,:)-cumsumSISS(1:longlen-20,:); %half-sec window, think of as centered at sample 10.
    SILBden=SILB2(21:longlen)-SILB2(1:longlen-20); %A long vector, 20 shorter than the long window
    SSIBden=SSIB2(21:longlen+2*mshift)-SSIB2(1:longlen+2*mshift-20); %A long vector, 20 shorter than the long window+2*mshift
    SISSn=zeros(longlen,2*mshift+1);  %Make as long as seismograms, but will remain zeros for first and last 10 samples
    for k=-mshift:mshift;
        SISSn(11:longlen-10,mshift+k+1)=SISSnum(:,mshift+k+1)./realsqrt(SILBden.*SSIBden(mshift+k+1:longlen+mshift+k-20));
    end
    SISSn(SISSn<0)=-2;

    ltmax=max([realSSIB(plotstart:plotend); realSILB(plotstart:plotend)]);
    ltmin=min([realSSIB(plotstart:plotend); realSILB(plotstart:plotend)]);
    npages=ceil((wins(n,2)-wins(n,1))/30.);
    for j=1:npages
        h=figure('Position',[wid 1 wid hite]); %center
        istart=plotstart+(j-1)*shortlen;
        iend=istart+shortlen;
        inormstart=istart-plotstart+1;
        inormend=inormstart+shortlen;

%                    SISSnum=cumsumSISStr(21:winlen)-cumsumSISStr(1:winlen-20);
%                    SSIBden=SSIB2(21:winlen)-SSIB2(1:winlen-20);
%                    subplot(3*nrow,mcol,3*(n-1)*mcol+2*mcol+m,'align');
%                    PGSSn=PGSSnum./realsqrt(PGCden.*SSIBden);
%                    PGSIn=PGSInum./realsqrt(PGCden.*SILBden);
%                    SISSn=SISSnum./realsqrt(SSIBden.*SILBden);
%                    alln=(PGSSn+PGSIn+SISSn)/3;
%                    idiff=SISSfile(nt,3);
%                    maxxc=max(alln(idiff:idiff+20)); %endpoints of alln window should be endpoints of 1 sec cumsum interval
%                    plot(PGCfile(winlen*(nt-1)+11:winlen*nt-10,1),alln,'c','linewidth',3)
%                    hold on
%                    plot(PGCfile(winlen*(nt-1)+11:winlen*nt-10,1),PGSSn,'b')
%                    plot(PGCfile(winlen*(nt-1)+11:winlen*nt-10,1),PGSIn,'r')
%                    plot(PGCfile(winlen*(nt-1)+11:winlen*nt-10,1),SISSn,'k')
%                    plot(PGCfile(winlen*(nt-1)+1:winlen*nt,1),0.75*maxxc*ones(winlen,1),'k:');
%                    plot(PGCfile(winlen*(nt-1)+1:winlen*nt,1),0.65*maxxc*ones(winlen,1),'k:');

        for k=-mshift:mshift
            subplot(11,1,k+mshift+1,'align'); 
            plot(timsPGC(istart:iend),SISSn(inormstart:inormend,k+mshift+1)*(ltmax-ltmin)*0.5+0.5*(ltmax+ltmin),'co','markersize',1)
            hold on
            plot(timsPGC(istart:iend),realSILB(istart:iend),'k');
            plot(timsPGC(istart:iend),realSSIB(istart+k:iend+k),'b');
            axis([timsPGC(istart) timsPGC(iend) ltmin ltmax]);
            box on
        end
       set(h,'PaperPosition',[0.25 0.25 8 10.5])
       orient landscape
       if j <= 9
           print(h,'-depsc',['ARMMAP/WIGS/',YEAR,'.',JDAY,'.',int2str(0),int2str(j),'.eps'])
       else
           print(h,'-depsc',['ARMMAP/WIGS/',YEAR,'.',JDAY,'.',int2str(j),'.eps'])
       end
    end
end

% subplot(4,1,2); 
% hold on
% plot(timsPGC,real(PGCrot),'r');
% plot(timsPGC,real(SSIBrot),'b');
% plot(timsPGC,real(SILBrot),'k');
% axis([timwin/2 timwin ltmin ltmax]);
% %xlim([timlook-37.5 timlook]);
% subplot(4,1,3); 
% hold on
% plot(timsPGC,real(PGCrot),'r');
% plot(timsPGC,real(SSIBrot),'b');
% plot(timsPGC,real(SILBrot),'k');
% axis([timwin 3*timwin/2 ltmin ltmax]);
% %xlim([timlook timlook+37.5]);
% subplot(4,1,4); 
% hold on
% plot(timsPGC,real(PGCrot),'r');
% plot(timsPGC,real(SSIBrot),'b');
% plot(timsPGC,real(SILBrot),'k');
% axis([3*timwin/2 2*timwin ltmin ltmax]);
% %xlim([timlook+37.5 timlook+75]);
% traceplt=cputime-t
% t=cputime;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   Autocorrelation of stations.  Those that end in "2" are the running
% %   cumulative sum, to be used later by differncing the window edpoints.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PGCauto=realPGC.*realPGC;
% PGC2=cumsum(PGCauto);
% SSIBauto=realSSIB.*realSSIB;
% SSIB2=cumsum(SSIBauto);
% SILBauto=realSILB.*realSILB;
% SILB2=cumsum(SILBauto);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  Cross-correlation between stations, with small offsets up to +/- mshift.
% %  First index is pointwise multiplication of traces; second is shifting offset.
% %  lenx is shorter than tracelenPGC by mshift at each end (see notebook sketch)
% %  For PGSS and PGSI, SSI and SIL are shifted relative to PGC, by 1 each time through loop.
% %  For SISS, SSI is shifted relative to SILB.
% %  cumsumPGSS etc. are the running cumulative sum of the x-correlation.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %mshift=16;
% lenx=tracelenPGC-2*mshift;
% PGSSx=zeros(lenx, 2*mshift+1);
% PGSIx=zeros(lenx, 2*mshift+1);
% SISSx=zeros(lenx, 2*mshift+1);
% for n=-mshift:mshift;
%     PGSSx(:,n+mshift+1)=realPGC(1+mshift:tracelenPGC-mshift).* ...
%         realSSIB(1+mshift-n:tracelenPGC-mshift-n);
%     PGSIx(:,n+mshift+1)=realPGC(1+mshift:tracelenPGC-mshift).* ...
%         realSILB(1+mshift-n:tracelenPGC-mshift-n);
%     SISSx(:,n+mshift+1)=realSILB(1+mshift:tracelenPGC-mshift).* ...
%         realSSIB(1+mshift-n:tracelenPGC-mshift-n);
% end
% cumsumPGSS=cumsum(PGSSx);
% cumsumPGSI=cumsum(PGSIx);
% cumsumSISS=cumsum(SISSx);
% xcshifts=cputime-t
% t=cputime;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  "winbig" is now the whole day, minus 3 sec at each end (apparently).
% %  "timbig" is the time of half that.
% %  igstart is the index of the starting sample.
% %  winlen and winoff refer to the small windows.
% %  timswin refers to the central times of those small windows.
% %  sumsPGSS (etc.) is the cross-correlation sum over the window.  The first
% %    index refers to the window number and the second the shift over +/-mshift.
% %  Normalized x-correlation:
% %    For PGSS and PGSI, for a given window PGC does not shift but SSI and 
% %    SIL do.  So can compute sumsPGC2 (from the running cum. sum PGC2) just
% %    once for each window.  Same for sumsSILB2b.  But for the stations that
% %    shift, SSI and SIL (for PGC) and SSI (for SIL), must compute sumsSSIB2 
% %    and sumsSILB2 upon each shift (actually, this is is easy book-keeping
% %    but not efficient).  Again, the first index refers to the window
% %    number and the second the shift over +/-mshift.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % winbig=2*(tracelenPGC/2-120);
% % timbig=winbig/(2*40);
% % igstart=floor(tracelenPGC/2-winbig/2)+1;
% % nwin=floor((winbig-winlen)/winoff);
% timswin=zeros(nwin,1);
% sumsPGSS=zeros(nwin,2*mshift+1);
% sumsPGSI=zeros(nwin,2*mshift+1);
% sumsSISS=zeros(nwin,2*mshift+1);
% sumsPGC2=zeros(nwin,2*mshift+1);
% sumsSSIB2=zeros(nwin,2*mshift+1);
% sumsSILB2=zeros(nwin,2*mshift+1);
% sumsSILB2b=zeros(nwin,2*mshift+1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  sumsPGSS is shorter than sumsPGC2 by 2*mshift.  This is why sumsPGC2 etc
% %  is shifted by +mshift.  cumsumPGSS(1,:)=cumsum(PGSSx)(1,:) starts mshift
% %  to the right of the first data sample.  igstart is how many to the right
% %  of that.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for n=1:nwin;
%     istart=igstart+(n-1)*winoff;
%     iend=istart+winlen;
%     timswin(n)=timsPGC(istart+winlen/2); 
%     sumsPGSS(n,:)=cumsumPGSS(iend,:)-cumsumPGSS(istart-1,:); 
%     sumsPGSI(n,:)=cumsumPGSI(iend,:)-cumsumPGSI(istart-1,:);
%     sumsSISS(n,:)=cumsumSISS(iend,:)-cumsumSISS(istart-1,:);
%     sumsPGC2(n,:)=PGC2(iend+mshift)-PGC2(istart+mshift-1);  %PGC2 is cumsummed. Yes, +mshift.
%     sumsSILB2b(n,:)=SILB2(iend+mshift)-SILB2(istart+mshift-1); %Similar, for the SILB-SSIB connection.
%     for m=-mshift:mshift;
%         sumsSSIB2(n,m+mshift+1)=SSIB2(iend+mshift-m)-SSIB2(istart+mshift-1-m); %+m??? (yes).
%         sumsSILB2(n,m+mshift+1)=SILB2(iend+mshift-m)-SILB2(istart+mshift-1-m);
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  Denominator for the normalization.  A 2D array, nwin by 2*mshift+1.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %An attempt to bypass glitches in data.  Min value of good data typically ~10^{-2}
% glitches=1.e-7;
% sumsPGC2=max(sumsPGC2,glitches);
% sumsSSIB2=max(sumsSSIB2,glitches);
% sumsSILB2=max(sumsSILB2,glitches);
% %
% denomPGSSn=realsqrt(sumsPGC2.*sumsSSIB2);
% denomPGSIn=realsqrt(sumsPGC2.*sumsSILB2);
% denomSISSn=realsqrt(sumsSILB2b.*sumsSSIB2);
% %
% % denomPGSSn=max(denomPGSSn,eps);
% % denomPGSIn=max(denomPGSIn,eps);
% % denomSISSn=max(denomSISSn,eps);
% %
% sumsPGSSn=sumsPGSS./denomPGSSn;
% sumsPGSIn=sumsPGSI./denomPGSIn;
% sumsSISSn=sumsSISS./denomSISSn;
% [xcmaxPGSSn,imaxPGSS]=max(sumsPGSSn,[],2); %Integer-offset max cross-correlation
% [xcmaxPGSIn,imaxPGSI]=max(sumsPGSIn,[],2);
% [xcmaxSISSn,imaxSISS]=max(sumsSISSn,[],2);
% %Parabolic fit:
% % [xmaxPGSSn,ymaxPGSSn,xloPGSSn,xhiPGSSn]=parabol(nwin,mshift,sumsPGSSn,imaxPGSS); %These return some poor measure of "confidence".
% % [xmaxPGSIn,ymaxPGSIn,xloPGSIn,xhiPGSIn]=parabol(nwin,mshift,sumsPGSIn,imaxPGSI);
% % [xmaxSISSn,ymaxSISSn,xloSISSn,xhiSISSn]=parabol(nwin,mshift,sumsSISSn,imaxSISS);
% [xmaxPGSSn,ymaxPGSSn,aPGSS]=parabol(nwin,mshift,sumsPGSSn,imaxPGSS); %Interpolated max cross-correlation
% [xmaxPGSIn,ymaxPGSIn,aPGSI]=parabol(nwin,mshift,sumsPGSIn,imaxPGSI);
% [xmaxSISSn,ymaxSISSn,aSISS]=parabol(nwin,mshift,sumsSISSn,imaxSISS);
% 
% ix=sub2ind(size(denomPGSSn),(1:nwin)',imaxPGSS);
% ampPGSS=sqrt(denomPGSSn(ix)); %This makes amplitude linear rather than quadratic with counts.
% ampPGC2=sumsPGC2(ix); %by construction PGC2 is the same for all shifts
% ampSSIB2=sumsSSIB2(ix);
% ix=sub2ind(size(denomPGSIn),(1:nwin)',imaxPGSI);
% ampPGSI=sqrt(denomPGSIn(ix));
% ampSILB2=sumsSILB2(ix);
% ix=sub2ind(size(denomSISSn),(1:nwin)',imaxSISS);
% ampSISS=sqrt(denomSISSn(ix));
% AmpComp(1:4)=0;
% AmpComp(5:nwin)=((ampPGC2(5:nwin)+ampSSIB2(5:nwin)+ampSILB2(5:nwin))- ...
%                 (ampPGC2(1:nwin-4)+ampSSIB2(1:nwin-4)+ampSILB2(1:nwin-4)))./ ...
%                 ((ampPGC2(5:nwin)+ampSSIB2(5:nwin)+ampSILB2(5:nwin))+ ...
%                 (ampPGC2(1:nwin-4)+ampSSIB2(1:nwin-4)+ampSILB2(1:nwin-4))) ;
% %Center them
% imaxPGSS=imaxPGSS-mshift-1;
% imaxPGSI=imaxPGSI-mshift-1;
% imaxSISS=imaxSISS-mshift-1;
% iloopoff=imaxPGSI-imaxPGSS+imaxSISS; %How well does the integer loop close?
% xmaxPGSSn=xmaxPGSSn-mshift-1;
% xmaxPGSIn=xmaxPGSIn-mshift-1;
% xmaxSISSn=xmaxSISSn-mshift-1;
% loopoff=xmaxPGSIn-xmaxPGSSn+xmaxSISSn; %How well does the interpolated loop close?
% xcmaxAVEn=(xcmaxPGSSn+xcmaxPGSIn+xcmaxSISSn)/3;
% xcnshifts=cputime-t
% t=cputime;
% 
% ampmax=max([ampPGSS; ampPGSI; ampSISS]);
% % figure 
% % subplot(4,1,1,'align'); 
% % hold on
% % plot(timswin,zeros(nwin,1),'k:');
% % %plot(timswin,3.5+zeros(nwin,1),'k:');
% % plot(timswin,xcmaxAVEn*mshift,'g');
% % plot(timswin,xmaxPGSSn,'bs','MarkerSize',2);
% % plot(timswin,xmaxPGSIn,'ro','MarkerSize',2);
% % plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',2);
% % %plot(timswin,xmaxSISSn,'k*','MarkerSize',2);
% % axis([0 timbig/2 -mshift mshift]);
% % ylabel('samples')
% % box on
% % subplot(4,1,2,'align'); 
% % hold on
% % plot(timswin,zeros(nwin,1),'k:');
% % %plot(timswin,3.5+zeros(nwin,1),'k:');
% % plot(timswin,xcmaxAVEn*mshift,'g');
% % plot(timswin,xmaxPGSSn,'bs','MarkerSize',2);
% % plot(timswin,xmaxPGSIn,'ro','MarkerSize',2);
% % plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',2);
% % %plot(timswin,xmaxSISSn,'k*','MarkerSize',2);
% % axis([timbig/2 timbig -mshift mshift]);
% % ylabel('samples')
% % box on
% % subplot(4,1,3,'align'); 
% % hold on
% % plot(timswin,zeros(nwin,1),'k:');
% % %plot(timswin,3.5+zeros(nwin,1),'k:');
% % plot(timswin,xcmaxAVEn*mshift,'g');
% % plot(timswin,xmaxPGSSn,'bs','MarkerSize',2);
% % plot(timswin,xmaxPGSIn,'ro','MarkerSize',2);
% % plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',2);
% % %plot(timswin,xmaxSISSn,'k*','MarkerSize',2);
% % axis([timbig 3*timbig/2 -mshift mshift]);
% % ylabel('samples')
% % box on
% % subplot(4,1,4,'align'); 
% % hold on
% % plot(timswin,zeros(nwin,1),'k:');
% % %plot(timswin,3.5+zeros(nwin,1),'k:');
% % plot(timswin,xcmaxAVEn*mshift,'g');
% % plot(timswin,xmaxPGSSn,'bs','MarkerSize',2);
% % plot(timswin,xmaxPGSIn,'ro','MarkerSize',2);
% % plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',2);
% % %plot(timswin,xmaxSISSn,'k*','MarkerSize',2);
% % axis([3*timbig/2 2*timbig -mshift mshift]);
% % xlabel('sec')
% % ylabel('samples')
% % box on
% % title([IDENTIF,'_{',int2str(lo),'-',int2str(hi),'}'])
% % orient landscape
% % print('-depsc',['FIGS/',IDENTIF,'-',int2str(winlen/40),'s_',int2str(lo),'-',int2str(hi),'a.eps'])
% 
% medxcmaxAVEn=median(xcmaxAVEn)
% xmaxPGSSntmp=xmaxPGSSn;
% xmaxPGSIntmp=xmaxPGSIn;
% xmaxSISSntmp=xmaxSISSn;
% %loopoff=xmaxPGSIn-xmaxPGSSn+xmaxSISSn; 
% iup=4;
% nin=0;
% zerosallowed=20*winlen/160;
% for n=1:nwin
% %     if timswin(n) > 35204
% %         [0 n timswin(n) xcmaxAVEn(n) xmaxPGSIn(n)-xmaxPGSSn(n)+xmaxSISSn(n)]
% %         pause
% %     end
%     if xcmaxAVEn(n)<xcmaxAVEnmin || abs(xmaxPGSIn(n)-xmaxPGSSn(n)+xmaxSISSn(n))>loopoffmax ...
%             || isequal(abs(imaxPGSS(n)),mshift) || isequal(abs(imaxPGSI(n)),mshift) || isequal(abs(imaxSISS(n)),mshift) ...
%             || nzerosPGC(n)>zerosallowed || nzerosSSIB(n)>zerosallowed || nzerosSILB(n)>zerosallowed
%         xmaxPGSSntmp(n)=20; xmaxPGSIntmp(n)=20; xmaxSISSntmp(n)=20; 
% %         if timswin(n) > 35204
% %             [1 n timswin(n) xcmaxAVEn(n) xmaxPGSIn(n)-xmaxPGSSn(n)+xmaxSISSn(n)]
% %             pause
% %         end
% %         if abs(timswin(n)-55158) <= 0.4
% %             [1 timswin(n) xcmaxAVEn(n) xmaxPGSIn(n)-xmaxPGSSn(n)+xmaxSISSn(n)]
% %             pause
% %         end
%     else
%         interpPGSSn=interp(sumsPGSSn(n,:),iup,3);
%         interpPGSIn=interp(sumsPGSIn(n,:),iup,3);
%         interpSISSn=interp(sumsSISSn(n,:),iup,3);
%         leninterp=length(interpPGSSn);
%         [xcmaxinterpPGSSn,imaxinterpPGSS]=max(interpPGSSn(1:leninterp-(iup-1)));
%         [xcmaxinterpPGSIn,imaxinterpPGSI]=max(interpPGSIn(1:leninterp-(iup-1)));
%         [xcmaxinterpSISSn,imaxinterpSISS]=max(interpSISSn(1:leninterp-(iup-1)));
%         if abs((imaxinterpPGSS-mshift*iup-1)/iup-imaxPGSS(n))>0.5 || ...
%            abs((imaxinterpPGSI-mshift*iup-1)/iup-imaxPGSI(n))>0.5 || ...
%            abs((imaxinterpSISS-mshift*iup-1)/iup-imaxSISS(n))>0.5
%                 [n (imaxinterpPGSS-mshift*iup-1)/iup imaxPGSS(n) ...
%                    (imaxinterpPGSI-mshift*iup-1)/iup imaxPGSI(n) ... 
%                    (imaxinterpSISS-mshift*iup-1)/iup imaxSISS(n)]
%         end
%         xcmaxconprev=-99999.;  %used to be 0; not good with glitches
%         for iPGSS=max(1,imaxinterpPGSS-3*iup):min(imaxinterpPGSS+3*iup,iup*(2*mshift+1)-(iup-1)) %3 samples from peak; 
%                                                                                  %intentionally wider than acceptable;
%                                                                                  %iup-1 are extrapolated points
%             for iPGSI=max(1,imaxinterpPGSI-3*iup):min(imaxinterpPGSI+3*iup,iup*(2*mshift+1)-(iup-1))
%                 ibangon = (iup*mshift+1)-iPGSI+iPGSS;
%                 if ibangon >= 1 && ibangon<=iup*(2*mshift+1)
%                     xcmaxcon=interpPGSSn(iPGSS)+interpPGSIn(iPGSI)+interpSISSn(ibangon);
%                     if xcmaxcon > xcmaxconprev
%                         xcmaxconprev=xcmaxcon;
%                         iPGSSbang=iPGSS;
%                         iPGSIbang=iPGSI;
%                     end
%                 end
%             end
%         end
%         iSISSbang=(iup*mshift+1)-iPGSIbang+iPGSSbang;
% %         if abs(timswin(n)-55158) <= 0.4
% %             [2 timswin(n) xcmaxAVEn(n) xmaxPGSIn(n)-xmaxPGSSn(n)+xmaxSISSn(n) ...
% %             iPGSSbang imaxinterpPGSS xcmaxinterpPGSSn iPGSIbang imaxinterpPGSI xcmaxinterpPGSIn ...
% %             iSISSbang imaxinterpSISS xcmaxinterpSISSn]
% %             pause
% %         end
%         if abs(iPGSSbang-imaxinterpPGSS) <= loopoffmax*iup && ...
%            abs(iPGSIbang-imaxinterpPGSI) <= loopoffmax*iup && ...
%            abs(iSISSbang-imaxinterpSISS) <= loopoffmax*iup && ...
%            interpPGSSn(iPGSSbang)+interpPGSIn(iPGSIbang)+interpSISSn(iSISSbang) >= 3*xcmaxAVEnmin
%             xmaxPGSSntmp(n)=(iPGSSbang-(iup*mshift+1))/iup;
%             xmaxPGSIntmp(n)=(iPGSIbang-(iup*mshift+1))/iup;
%             xmaxSISSntmp(n)=(iSISSbang-(iup*mshift+1))/iup;
%             
%             %for plotting traces
%             imaxPGSSwr=round(xmaxPGSSntmp(n));
%             imaxPGSIwr=round(xmaxPGSIntmp(n));
% 
%             istart=igstart+(n-1)*winoff+mshift; %a better way might exist?  %ADDED mshift 10/20/12
%             iend=istart+winlen-1;
%             imid=round((istart+iend)/2);
%             %Check power spectrum
%             [PGCxx fp] = pwelch(realPGC(istart:iend),[],[],[],40); %40 is sps
%             PGCxx=PGCxx/max(PGCxx);
%             [SSIBxx fp] = pwelch(realSSIB(istart-imaxPGSSwr:iend-imaxPGSSwr),[],[],[],40);
%             SSIBxx=SSIBxx/max(SSIBxx);
%             [SILBxx fp] = pwelch(realSILB(istart-imaxPGSIwr:iend-imaxPGSIwr),[],[],[],40);
%             SILBxx=SILBxx/max(SILBxx);
%             flo=find(fp > lo,1)-1;
%             fhi=find(fp > hi,1)+1; %extra 1 for good measure
%             belowcut=median([PGCxx(2:flo); SSIBxx(2:flo); SILBxx(2:flo)]);
%             ppeaksPGC=findpeaks(PGCxx(flo+1:fhi));
%             if length(ppeaksPGC)>=1
%                 maxppeakPGC=max(ppeaksPGC);
%             else
%                 maxppeakPGC=0.;
%             end
%             ppeaksSSIB=findpeaks(SSIBxx(flo+1:fhi));
%             if length(ppeaksSSIB)>=1
%                 maxppeakSSIB=max(ppeaksSSIB);
%             else
%                 maxppeakSSIB=0.;
%             end
%             ppeaksSILB=findpeaks(SILBxx(flo+1:fhi));
%             if length(ppeaksSILB)>=1
%                 maxppeakSILB=max(ppeaksSILB);
%             else
%                 maxppeakSILB=0.;
%             end
%             abovecut=median([maxppeakPGC maxppeakSSIB maxppeakSILB]);
%             if abovecut > 0.9*belowcut %-1
%                 PGSStr=realPGC(istart:iend).*realSSIB(istart-imaxPGSSwr:iend-imaxPGSSwr);
%                 PGSItr=realPGC(istart:iend).*realSILB(istart-imaxPGSIwr:iend-imaxPGSIwr);
%                 SISStr=realSILB(istart-imaxPGSIwr:iend-imaxPGSIwr).*realSSIB(istart-imaxPGSSwr:iend-imaxPGSSwr);
%                 cumsumtr=cumsum(PGSStr)+cumsum(PGSItr)+cumsum(SISStr);
%                 [cumsumtrdiff idiff]=max(cumsumtr(41:winlen)-cumsumtr(1:winlen-40));
% 
%                 PGCfile(nin*winlen+1:(nin+1)*winlen,1:2)=[timsPGC(istart:iend)' realPGC(istart:iend)];
%                 SSIBfile(nin*winlen+1:(nin+1)*winlen,1:2)=[timsPGC(istart:iend)' realSSIB(istart-imaxPGSSwr:iend-imaxPGSSwr)];
%                 SILBfile(nin*winlen+1:(nin+1)*winlen,1:2)=[timsPGC(istart:iend)' realSILB(istart-imaxPGSIwr:iend-imaxPGSIwr)];
%                 PGSSfile(nin+1,1:2)=[imaxPGSSwr xcmaxPGSSn(n)];
%                 PGSIfile(nin+1,1:2)=[imaxPGSIwr xcmaxPGSIn(n)];
%                 SISSfile(nin+1,1:3)=[cumsumtrdiff/cumsumtr(winlen) xcmaxSISSn(n) idiff];
% 
%                 nin=nin+1;
%                 aPGSSkeep(nin,:)=[timswin(n) aPGSS(n)];
%                 aPGSIkeep(nin,:)=[timswin(n) aPGSI(n)];
%                 aSISSkeep(nin,:)=[timswin(n) aSISS(n)];
%                 loopoffkeep(nin,:)=[timswin(n) loopoff(n)];
%                 mapfile(nin,:)=[timswin(n) xmaxPGSIntmp(n) xmaxPGSSntmp(n) ...
%                     xcmaxAVEn(n) loopoff(n) AmpComp(n) cumsumtrdiff timswin(n)-2+idiff/40. cumsumtrdiff/cumsumtr(winlen)];
% %                 if timswin(n) > 35204
% %                     [2 n timswin(n) xcmaxAVEn(n)]
% %                     pause
% %                 end
%             else
%                 xmaxPGSSntmp(n)=20; xmaxPGSIntmp(n)=20; xmaxSISSntmp(n)=20;
% %                 if timswin(n) > 35204
% %                     [3 n timswin(n) xcmaxAVEn(n)]
% %                     pause
% %                 end
%            end
%         else
%             xmaxPGSSntmp(n)=20; xmaxPGSIntmp(n)=20; xmaxSISSntmp(n)=20; 
% %             if timswin(n) > 35204
% %                 [4 n timswin(n) xcmaxAVEn(n)]
% %                 pause
% %             end
%         end
%     end
% end
% %fid = fopen(['ARMMAP/MAPS/map',IDENTIF,'_',int2str(lo),'-',int2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
% fid = fopen(['ARMMAP/map',IDENTIF,'_',int2str(lo),'-',int2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
% fprintf(fid,'%9.1f %9.5f %9.5f %9.4f %9.2f %9.3f %10.4f %10.4f %8.3f\n',mapfile(1:nin,:)');
% fclose(fid);
% figure 
% subplot(4,1,1,'align'); 
% hold on
% plot(timswin,xcmaxAVEnmin*mshift+zeros(nwin,1),'k:');
% plot(timsPGC(winlen:2*winlen),7+zeros(winlen+1,1),'k','linewidth',2);
% plot(timswin,zeros(nwin,1),'k:');
% plot(timswin,xcmaxAVEn*mshift,'g');
% plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
% plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
% plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
% axis([0 timbig/2 -mshift mshift]);
% ylabel('samples')
% box on
% subplot(4,1,2,'align'); 
% hold on
% plot(timswin,xcmaxAVEnmin*mshift+zeros(nwin,1),'k:');
% plot(timswin,zeros(nwin,1),'k:');
% plot(timswin,xcmaxAVEn*mshift,'g');
% plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
% plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
% plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
% axis([timbig/2 timbig -mshift mshift]);
% ylabel('samples')
% box on
% subplot(4,1,3,'align'); 
% hold on
% plot(timswin,xcmaxAVEnmin*mshift+zeros(nwin,1),'k:');
% plot(timswin,zeros(nwin,1),'k:');
% plot(timswin,xcmaxAVEn*mshift,'g');
% plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
% plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
% plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
% axis([timbig 3*timbig/2 -mshift mshift]);
% ylabel('samples')
% box on
% subplot(4,1,4,'align'); 
% hold on
% plot(timswin,xcmaxAVEnmin*mshift+zeros(nwin,1),'k:');
% plot(timswin,zeros(nwin,1),'k:');
% plot(timswin,xcmaxAVEn*mshift,'g');
% plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
% plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
% plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
% axis([3*timbig/2 2*timbig -mshift mshift]);
% xlabel('sec')
% ylabel('samples')
% box on
% title([IDENTIF,'_{',int2str(lo),'-',int2str(hi),'}'])
% orient landscape
% print('-depsc',['FIGS/',IDENTIF,'-',int2str(winlen/40),'s_',int2str(lo),'-',int2str(hi),'b.eps'])
% fid = fopen(['HILBERTS/xcmax',IDENTIF,'_',int2str(lo),'-',int2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
% fprintf(fid,'%9.3f %9.5f\n',[timswin xcmaxAVEn]');
% fclose(fid);
% 
% for n=1:nwin;
%     if abs(loopoff(n))>2.0
%         loopoff(n)=30; 
%     end
% end
% figure 
% subplot(4,1,1,'align'); 
% hold on
% plot(timswin,zeros(nwin,1),'k:');
% plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
% plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
% plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
% plot(timsPGC(winlen:2*winlen),7+zeros(winlen+1,1),'k','linewidth',2);
% axis([0 timbig/2 -mshift mshift]);
% title('blue = PGSS ;   red = PGSI ;  black = SISS')
% box on
% subplot(4,1,2,'align'); 
% hold on
% hrf = plotreflinesr(gca,detects,'x','k'); 
% hrf = plotreflinesr(gca,Detects,'x','r'); 
% hrf = plotreflinesr(gca,detects2_8,'x','g'); 
% hrf = plotreflinesr(gca,Detects2_8,'x','b'); 
% plot(timswin,xcmaxAVEnmin+zeros(nwin,1),'k:');
% plot(timswin,xcmaxAVEn,'k');
% plot(timswin,-1+0.1*abs(loopoff),'g.','MarkerSize',1);
% plot(timsPGC,min(0,real(PGCrot)),'r');
% plot(timsPGC,min(0,real(SSIBrot)),'b');
% plot(timsPGC,min(0,real(SILBrot)),'k');
% axis([0 timbig/2 -1 1]);
% box on
% subplot(4,1,3,'align'); 
% hold on
% plot(timswin,zeros(nwin,1),'k:');
% plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
% plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
% plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
% axis([timbig/2 timbig -mshift mshift]);
% box on
% subplot(4,1,4,'align'); 
% hold on
% hrf = plotreflinesr(gca,detects,'x','k'); 
% hrf = plotreflinesr(gca,Detects,'x','r'); 
% hrf = plotreflinesr(gca,detects2_8,'x','g'); 
% hrf = plotreflinesr(gca,Detects2_8,'x','b'); 
% plot(timswin,xcmaxAVEnmin+zeros(nwin,1),'k:');
% plot(timswin,xcmaxAVEn,'k');
% plot(timswin,-1+0.1*abs(loopoff),'g.','MarkerSize',1);
% plot(timsPGC,min(0,real(PGCrot)),'r');
% plot(timsPGC,min(0,real(SSIBrot)),'b');
% plot(timsPGC,min(0,real(SILBrot)),'k');
% axis([timbig/2 timbig -1 1]);
% box on
% title([IDENTIF,'_{',int2str(lo),'-',int2str(hi),'}'])
% orient landscape
% print('-depsc',['FIGS/',IDENTIF,'-',int2str(winlen/40),'s_',int2str(lo),'-',int2str(hi),'c.eps'])
% 
% figure
% subplot(4,1,1,'align'); 
% hold on
% plot(timswin,zeros(nwin,1),'k:');
% plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
% plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
% plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
% plot(timsPGC(tracelenPGC/2+winlen:tracelenPGC/2+2*winlen),7+zeros(winlen+1,1),'k','linewidth',2);
% axis([timbig 3*timbig/2 -mshift mshift]);
% box on
% subplot(4,1,2,'align'); 
% hold on
% hrf = plotreflinesr(gca,detects,'x','k'); 
% hrf = plotreflinesr(gca,Detects,'x','r'); 
% hrf = plotreflinesr(gca,detects2_8,'x','g'); 
% hrf = plotreflinesr(gca,Detects2_8,'x','b'); 
% plot(timswin,xcmaxAVEnmin+zeros(nwin,1),'k:');
% plot(timswin,xcmaxAVEn,'k');
% plot(timswin,-1+0.1*abs(loopoff),'g.','MarkerSize',1);
% plot(timsPGC,min(0,real(PGCrot)),'r');
% plot(timsPGC,min(0,real(SSIBrot)),'b');
% plot(timsPGC,min(0,real(SILBrot)),'k');
% axis([timbig 3*timbig/2 -1 1]);
% box on
% subplot(4,1,3,'align'); 
% hold on
% plot(timswin,zeros(nwin,1),'k:');
% plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
% plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
% plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
% axis([3*timbig/2 2*timbig -mshift mshift]);
% box on
% subplot(4,1,4,'align'); 
% hold on
% hrf = plotreflinesr(gca,detects,'x','k'); 
% hrf = plotreflinesr(gca,Detects,'x','r'); 
% hrf = plotreflinesr(gca,detects2_8,'x','g'); 
% hrf = plotreflinesr(gca,Detects2_8,'x','b'); 
% plot(timswin,xcmaxAVEnmin+zeros(nwin,1),'k:');
% plot(timswin,xcmaxAVEn,'k');
% plot(timswin,-1+0.1*abs(loopoff),'g.','MarkerSize',1);
% plot(timsPGC,min(0,real(PGCrot)),'r');
% plot(timsPGC,min(0,real(SSIBrot)),'b');
% plot(timsPGC,min(0,real(SILBrot)),'k');
% axis([3*timbig/2 2*timbig -1 1]);
% box on
% title([IDENTIF,'_{',int2str(lo),'-',int2str(hi),'}'])
% orient landscape
% print('-depsc',['FIGS/',IDENTIF,'-',int2str(winlen/40),'s_',int2str(lo),'-',int2str(hi),'d.eps'])
% 
% figure
% colormap(jet)
% scatter(xmaxPGSIn-xmaxPGSSn+xmaxSISSn,xcmaxAVEn,3,AmpComp)
% hold on 
% plot(-5:5,xcmaxAVEnmin+zeros(11,1),'k:');
% axis([-5 5 -0.2 1.0])
% hrf = plotreflinesr(gca,-1.5,'x','k');colorbar
% hrf = plotreflinesr(gca,1.5,'x','k');colorbar
% box on
% title([IDENTIF,'_{',int2str(lo),'-',int2str(hi),'}'])
% print('-depsc',['FIGS/',IDENTIF,'-',int2str(winlen/40),'s_',int2str(lo),'-',int2str(hi),'e.eps'])
% 
% if winlen<=200
% scrsz=get(0,'ScreenSize');
% nt=0;
% nrow=4;
% mcol=6;
% for ifig=1:floor(nin/(nrow*mcol))+1
%     figure('Position',[scrsz(3)/10 scrsz(4)/10 4*scrsz(3)/5 9*scrsz(4)/10]);
%     for n = 1:nrow
%         for m = 1:mcol
%             nt=nt+1;
%             if nt <= nin
%                  %if PGSSfile(nt,1) >= 10 && PGSSfile(nt,1) <= 16 && PGSIfile(nt,1) >= 2 && PGSIfile(nt,1) <= 8
%                     subplot(3*nrow,mcol,3*(n-1)*mcol+m,'align');
%                     plot(PGCfile(winlen*(nt-1)+1:winlen*nt,1),PGCfile(winlen*(nt-1)+1:winlen*nt,2),'r')
%                     hold on
%                     plot(SSIBfile(winlen*(nt-1)+1:winlen*nt,1),SSIBfile(winlen*(nt-1)+1:winlen*nt,2),'b')
%                     plot(SILBfile(winlen*(nt-1)+1:winlen*nt,1),SILBfile(winlen*(nt-1)+1:winlen*nt,2),'k')
%                     is = PGCfile(winlen*(nt-1)+1,1);
%                     ien= PGCfile(winlen*nt,1);
%                     %yma=0.4;
%                     yma=max(max([PGCfile(winlen*(nt-1)+1:winlen*nt,2) SSIBfile(winlen*(nt-1)+1:winlen*nt,2) ...
%                         SILBfile(winlen*(nt-1)+1:winlen*nt,2)]));
%                     ymi=min(min([PGCfile(winlen*(nt-1)+1:winlen*nt,2) SSIBfile(winlen*(nt-1)+1:winlen*nt,2) ...
%                         SILBfile(winlen*(nt-1)+1:winlen*nt,2)]));
%                     xvect=[is is+2*(yma-ymi)*(winlen/160.)]; %amplitude bar originally scaled for 4-s window
%                     yma=2.4*max(yma,-ymi);
%                     yvect=[-0.9*yma -0.9*yma];
%                     plot(xvect,yvect,'r','linewidth',3)
%                     plot([is+1/hi is+1/lo],[-0.8*yma -0.8*yma],'k','linewidth',3)
%                     text(is+0.2, 0.66*yma, int2str(PGSIfile(nt,1)),'fontsize',6);
%                     text(ien-0.6, 0.66*yma, int2str(PGSSfile(nt,1)),'fontsize',6);
%                     box on
%                     axis([is ien -yma yma])
%                     set(gca,'XTick',[is is+2],'fontsize',6);
% 
%                     subplot(3*nrow,mcol,3*(n-1)*mcol+mcol+m,'align');
% % The following for spectra
% %                     [PGCxx f] = pwelch(PGCfile(winlen*(nt-1)+1:winlen*nt,2)-mean(PGCfile(winlen*(nt-1)+1:winlen*nt,2)),[],[],[],40);
% %                     [SSIBxx f] = pwelch(SSIBfile(winlen*(nt-1)+1:winlen*nt,2)-mean(SSIBfile(winlen*(nt-1)+1:winlen*nt,2)),[],[],[],40);
% %                     [SILBxx f] = pwelch(SILBfile(winlen*(nt-1)+1:winlen*nt,2)-mean(SILBfile(winlen*(nt-1)+1:winlen*nt,2)),[],[],[],40);
% % %                     [PGCxx f] = pwelch(PGCfile(winlen*(nt-1)+1:winlen*nt,2),[],[],[],40);
% % %                     [SSIBxx f] = pwelch(SSIBfile(winlen*(nt-1)+1:winlen*nt,2),[],[],[],40);
% % %                     [SILBxx f] = pwelch(SILBfile(winlen*(nt-1)+1:winlen*nt,2),[],[],[],40);
% % %                     [Cpgss,f] = mscohere(PGCfile(winlen*(nt-1)+1:winlen*nt,2),SSIBfile(winlen*(nt-1)+1:winlen*nt,2), ...
% % %                         [],[],[],40);
% % %                     [Cpgsi,f] = mscohere(PGCfile(winlen*(nt-1)+1:winlen*nt,2),SILBfile(winlen*(nt-1)+1:winlen*nt,2), ...
% % %                         [],[],[],40);
% % %                     [Csiss,f] = mscohere(SILBfile(winlen*(nt-1)+1:winlen*nt,2),SSIBfile(winlen*(nt-1)+1:winlen*nt,2), ...
% % %                         [],[],[],40);
% % %                     plot(f,Cpgss,'b')
% % %                     hold on
% % %                     plot(f,Cpgsi,'r')
% % %                     plot(f,Csiss,'k')
% %                     plot(f,PGCxx/max(PGCxx),'r+')
% %                     hold on
% %                     plot(f,SSIBxx/max(SSIBxx),'b+')
% %                     plot(f,SILBxx/max(SILBxx),'k+')
% %                     xlim([0 hi+2])
% %                     %xlim([0.1 hi+1])
% %                     %ylim([0.01 1])
% %                     set(gca,'XTick',(1:10),'fontsize',6);
% % The above for spectra
%                     PGSStr=PGCfile(winlen*(nt-1)+1:winlen*nt,2).*SSIBfile(winlen*(nt-1)+1:winlen*nt,2);            
%                     PGSItr=PGCfile(winlen*(nt-1)+1:winlen*nt,2).*SILBfile(winlen*(nt-1)+1:winlen*nt,2);
%                     SISStr=SSIBfile(winlen*(nt-1)+1:winlen*nt,2).*SILBfile(winlen*(nt-1)+1:winlen*nt,2);
%                     % find the peaks etc.
%                     avedots=(PGSStr+PGSItr+SISStr)/3.;
%                     [peaks, locs]=findpeaks(avedots,'minpeakdistance',3);
%                     npeaks=length(peaks);
%                     [maxpk, imaxpk]=max(peaks);
%                     pks=zeros(npeaks,2);
%                     pks(:,2)=peaks;
%                     pks(:,1)=PGCfile(winlen*(nt-1)+locs);
%                     pksort=sortrows(pks,2);
%                     rat14=maxpk/pksort(npeaks-4,2); %ratio of max to 4th largest anywhere in window
%                     if imaxpk==1 
%                         maxpkp=peaks(2);
%                         maxpkm=-9e9;
%                         pkwid=2*(pks(2,1)-pks(1,1));
%                         pksid12=maxpk/maxpkp;
%                         pksid13=maxpk/peaks(3);
%                     elseif imaxpk==npeaks
%                         maxpkp=-9e9;
%                         maxpkm=peaks(npeaks-1);
%                         pkwid=2*(pks(npeaks,1)-pks(npeaks-1,1));
%                         pksid12=maxpk/maxpkm;
%                         pksid13=maxpk/peaks(npeaks-2);
%                     else
%                         maxpkp=peaks(imaxpk+1);
%                         maxpkm=peaks(imaxpk-1);
%                         pkwid=pks(imaxpk+1,1)-pks(imaxpk-1,1);
%                         pksid12=maxpk/max(maxpkp,maxpkm);
%                         pksid13=maxpk/min(maxpkp,maxpkm);
%                     end
%                     cumsumPGSStr=cumsum(PGSStr);
%                     cumsumPGSItr=cumsum(PGSItr);
%                     cumsumSISStr=cumsum(SISStr);
%                      yma=1.1; %yma=max(avedots);
%                      ymi=-0.1; %ymi=min(avedots);
%                     %hold on
% %                     plot(PGCfile(winlen*(nt-1)+1:winlen*nt,1),PGSStr,'b')
% %                     plot(PGCfile(winlen*(nt-1)+1:winlen*nt,1),PGSItr,'r')
% %                     plot(PGCfile(winlen*(nt-1)+1:winlen*nt,1),SISStr,'k')
% %                     plot(PGCfile(winlen*(nt-1)+1:winlen*nt,1),avedots,'c') %,'linewidth',2)
% %                     plot(PGCfile(winlen*(nt-1)+1:winlen*nt,1),zeros(winlen,1),'k')
% %                     plot(pks(:,1),pks(:,2),'r+')
% % The following for running-sum dot product
%                     plot(PGCfile(winlen*(nt-1)+1:winlen*nt,1),cumsumPGSStr,'b')
%                     hold on
%                     plot(PGCfile(winlen*(nt-1)+1:winlen*nt,1),cumsumPGSItr,'r')
%                     plot(PGCfile(winlen*(nt-1)+1:winlen*nt,1),cumsumSISStr,'k')
%                     yma=max(max([cumsumPGSStr cumsumPGSItr cumsumSISStr]));
%                     ymi=min(min([cumsumPGSStr cumsumPGSItr cumsumSISStr]));
%                     axis([is ien ymi yma])
%                     set(gca,'XTick',(0:20),'fontsize',6);
% % The above for running-sum dot product
%                     box on
%                     text(is+0.1, ymi+0.82*(yma-ymi), num2str(PGSIfile(nt,2),3),'fontsize',6);
%                     text(is+0.1, ymi+0.64*(yma-ymi), num2str(PGSSfile(nt,2),3),'fontsize',6);
%                     text(is+0.1, ymi+0.46*(yma-ymi), num2str(SISSfile(nt,2),3),'fontsize',6);
%                     text(ien-0.8, ymi+0.1*(yma-ymi), num2str(SISSfile(nt,1),3),'fontsize',6);
%                     %axis([lo-1 hi+1 ymi yma])
%                     %axis tight
%                     %set(gca,'XTick',[is is+2],'fontsize',6);
%                     pkfile(nt,:)=[PGCfile(winlen*(nt-1)+1+(winlen/2)) pks(imaxpk,1) maxpk pksid12 pksid13 pkwid rat14];
% 
%                     PGCauto=PGCfile(winlen*(nt-1)+1:winlen*nt,2).*PGCfile(winlen*(nt-1)+1:winlen*nt,2);
%                     PGC2=cumsum(PGCauto);
%                     SSIBauto=SSIBfile(winlen*(nt-1)+1:winlen*nt,2).*SSIBfile(winlen*(nt-1)+1:winlen*nt,2);
%                     SSIB2=cumsum(SSIBauto);
%                     SILBauto=SILBfile(winlen*(nt-1)+1:winlen*nt,2).*SILBfile(winlen*(nt-1)+1:winlen*nt,2);
%                     SILB2=cumsum(SILBauto);
%                     PGSSnum=cumsumPGSStr(21:winlen)-cumsumPGSStr(1:winlen-20);
%                     PGSInum=cumsumPGSItr(21:winlen)-cumsumPGSItr(1:winlen-20);
%                     SISSnum=cumsumSISStr(21:winlen)-cumsumSISStr(1:winlen-20);
%                     PGCden=PGC2(21:winlen)-PGC2(1:winlen-20);
%                     SSIBden=SSIB2(21:winlen)-SSIB2(1:winlen-20);
%                     SILBden=SILB2(21:winlen)-SILB2(1:winlen-20);
%                     subplot(3*nrow,mcol,3*(n-1)*mcol+2*mcol+m,'align');
%                     PGSSn=PGSSnum./realsqrt(PGCden.*SSIBden);
%                     PGSIn=PGSInum./realsqrt(PGCden.*SILBden);
%                     SISSn=SISSnum./realsqrt(SSIBden.*SILBden);
%                     alln=(PGSSn+PGSIn+SISSn)/3;
%                     idiff=SISSfile(nt,3);
%                     maxxc=max(alln(idiff:idiff+20)); %endpoints of alln window should be endpoints of 1 sec cumsum interval
%                     plot(PGCfile(winlen*(nt-1)+11:winlen*nt-10,1),alln,'c','linewidth',3)
%                     hold on
%                     plot(PGCfile(winlen*(nt-1)+11:winlen*nt-10,1),PGSSn,'b')
%                     plot(PGCfile(winlen*(nt-1)+11:winlen*nt-10,1),PGSIn,'r')
%                     plot(PGCfile(winlen*(nt-1)+11:winlen*nt-10,1),SISSn,'k')
%                     plot(PGCfile(winlen*(nt-1)+1:winlen*nt,1),0.75*maxxc*ones(winlen,1),'k:');
%                     plot(PGCfile(winlen*(nt-1)+1:winlen*nt,1),0.65*maxxc*ones(winlen,1),'k:');
%                     axis([is ien -0.5 1])
%                     set(gca,'XTick',(0:2),'fontsize',6);
%                  %end
%             end
%         end
%     end
%     orient landscape
%     if ifig <= 9
%         print('-depsc',['ARMMAP/WIGS/',IDENTIF,'-',int2str(lo),'-',int2str(hi),'.',int2str(0),int2str(0),int2str(ifig),'.eps'])
%     elseif ifig <= 99
%         print('-depsc',['ARMMAP/WIGS/',IDENTIF,'-',int2str(lo),'-',int2str(hi),'.',int2str(0),int2str(ifig),'.eps'])
%     else
%         print('-depsc',['ARMMAP/WIGS/',IDENTIF,'-',int2str(lo),'-',int2str(hi),'.',int2str(ifig),'.eps'])
%     end
% 
% end
% %fid = fopen(['ARMMAP/MAPS/pks',IDENTIF,'_',int2str(lo),'-',int2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
% fid = fopen(['ARMMAP/pks',IDENTIF,'_',int2str(lo),'-',int2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
% fprintf(fid,'%9.1f %9.3f %10.6f %9.3f %9.3f %9.3f %9.3f \n',pkfile(1:nin,:)');
% fclose(fid);
% end
% 
% medlok=median(abs(loopoffkeep))
% medaPGSS=median(aPGSSkeep)
% medaPGSI=median(aPGSIkeep)
% medaSISS=median(aSISSkeep)
% 
% 
% %cputime-t;
% tot=cputime-tt
% end
