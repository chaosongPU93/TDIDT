%Reads in seismograms and Bostock detections. Filters and rotates seisms, etc.  Stacks based on 
%times of B. detections, using shifts of up to 5 samples to maximize cross-correlation 
%of the central half of the window (at SSIB and SILB only). Puts SILB arrival at
%center of 5-sec window.  Isolates that time +/-12 samples.  Plots energy in that 
%0.6-sec window (both coherent and amp. squared) compared to the total in
%window.  Writes to .En files in BOSTOCK/NEW, with columns 7,8,9 being
%coherent energy (1st 2.2 sec, central 0.6 sec, last 2.2 sec), and 10,11,12
%being the same for amp squared.  

format short e
clear all
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');
scrsz=get(0,'ScreenSize');
% wid=scrsz(3)/3.1;
% hite=scrsz(4);
% scrat=wid/hite;
% h=figure('Position',[wid 1 wid hite]); %center

    timoffrot(1,:)=[2003 062 07 42 12 +086 +020  80 120  50]; %120 for 2-8 Hz
    timoffrot(2,:)=[2003 063 07 42 12 +086 +020  80 115  50];
    timoffrot(3,:)=[2003 064 07 42 12 +086 +020  80 115  50];
    
    timoffrot(4,:)=[2004 196 07 42 12 +086 +020  80 115  50];
    timoffrot(5,:)=[2004 197 07 42 12 +086 +020  80 115  50];
    timoffrot(6,:)=[2004 198 07 42 12 +086 +020  80 115  50];
    timoffrot(7,:)=[2004 199 07 42 12 +086 +020  80 115  50];
   
    timoffrot(8,:)=[2005 254 07 42 12 +086 +020  80 115  50];
    timoffrot(9,:)=[2005 255 07 42 12 +086 +020  80 115  50];
    timoffrot(10,:)=[2005 256 07 42 12 +086 +020  80 115  50];
    timoffrot(11,:)=[2005 257 07 42 12 +086 +020  80 115  50];

%246  30304  1 3284.125 1.632 10
%  2  30304  2 3108.825 1.681 13
%  2  30304  2 3115.600 1.727 15
%  2  30304  2 3127.925 1.697 17
     bostname(1,:)='BOSTOCK/NEW/002_2003.062';
     bostname(2,:)='BOSTOCK/NEW/002_2003.063';
     bostname(3,:)='BOSTOCK/NEW/002_2003.064';

     bostname(4,:)='BOSTOCK/NEW/002_2004.196';
     bostname(5,:)='BOSTOCK/NEW/002_2004.197';
     bostname(6,:)='BOSTOCK/NEW/002_2004.198';
     bostname(7,:)='BOSTOCK/NEW/002_2004.199';

     bostname(8,:)='BOSTOCK/NEW/002_2005.254';
     bostname(9,:)='BOSTOCK/NEW/002_2005.255';
     bostname(10,:)='BOSTOCK/NEW/002_2005.256';
     bostname(11,:)='BOSTOCK/NEW/002_2005.257';

hi=6.;
lo=1.5;
hihi=8.;
lolo=0.5;
winlen=4.*40;  %Full window
%stack=zeros(size(winlen));
%stackSI=stack';
%stackSS=stack';
%stackPG=stack';

for nd=1:length(timoffrot(:,1))

    bostocks=load(bostname(nd,:));

    bostsec=3600*(bostocks(:,3)-1)+bostocks(:,4)+7*0.025; %(the latter to correct for the diff between Bostock and Rubin filtering.)
    boststart=round(bostsec*40)-winlen/2;
    bostend=boststart+winlen-1;

    year=timoffrot(nd,1);
    YEAR=int2str(year);
    jday=timoffrot(nd,2);
    if jday <= 9
        JDAY=['00',int2str(jday)];
    elseif jday<= 99
        JDAY=['0',int2str(jday)];
    else
        JDAY=int2str(jday);
    end
    MO=day2month(jday,year);
    rotPGC=pi*(timoffrot(nd,8)-90)/180;
    rotSSIB=pi*(timoffrot(nd,9)-90)/180;
    rotSILB=pi*(timoffrot(nd,10)-90)/180;
    SSIBsoff=timoffrot(nd,6);
    SILBsoff=timoffrot(nd,7);
    SSIBtoff=SSIBsoff/40;
    SILBtoff=SILBsoff/40;
    IDENTIF=[YEAR,'.',JDAY,'.',int2str(timoffrot(nd,6)),'.',int2str(timoffrot(nd,7)), ...
        '.',int2str(timoffrot(nd,8)),'.',int2str(timoffrot(nd,9)),'.',int2str(timoffrot(nd,10))]

    %Read data:
    direc=[YEAR,'/',MO,'/'];
    prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];
    PGCEdat=[prename,'.PGC..BHE.D.SAC'];
    PGCNdat=[prename,'.PGC..BHN.D.SAC'];
%   PGCZdat=[prename,'.PGC..BHZ.D.SAC'];
    SSIBEdat=[prename,'.SSIB..HHE.D.SAC'];
    SSIBNdat=[prename,'.SSIB..HHN.D.SAC'];
%   SSIBZdat=[prename,'.SSIB..HHZ.D.SAC'];
    SILBEdat=[prename,'.SILB..HHE.D.SAC'];
    SILBNdat=[prename,'.SILB..HHN.D.SAC'];
%   SILBZdat=[prename,'.SILB..HHZ.D.SAC'];

    [PGCE,HdrDataPGC,tnuPGC,pobjPGC,timsPGC]=readsac(PGCEdat,0,'l');
    [PGCN,~,~,~,~]=readsac(PGCNdat,0,'l');
%   [PGCZ,~,~,~,~]=readsac(PGCZdat,0,'l');
    [SSIBE,HdrDataSSIB,tnuSSIB,pobjSSIB,timsSSIB]=readsac(SSIBEdat,0,'l');
    [SSIBN,~,~,~,~]=readsac(SSIBNdat,0,'l');
%   [SSIBZ,~,~,~,~]=readsac(SSIBZdat,0,'l');
    [SILBE,HdrDataSILB,tnuSILB,pobjSILB,timsSILB]=readsac(SILBEdat,0,'l');
    [SILBN,~,~,~,~]=readsac(SILBNdat,0,'l');
%   [SILBZ,~,~,~,~]=readsac(SILBZdat,0,'l');

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
%   PGCZ(1:40)=sin(x).*PGCZ(1:40);
    x=flipud(x);
    PGCE(tracelenPGC-39:tracelenPGC)=sin(x).*PGCE(tracelenPGC-39:tracelenPGC);
    PGCN(tracelenPGC-39:tracelenPGC)=sin(x).*PGCN(tracelenPGC-39:tracelenPGC);
%   PGCZ(tracelenPGC-39:tracelenPGC)=sin(x).*PGCZ(tracelenPGC-39:tracelenPGC);
    x=(0:pi/200:pi/2-pi/200)';
    SSIBE(1:100)=sin(x).*SSIBE(1:100);
    SSIBN(1:100)=sin(x).*SSIBN(1:100);
%   SSIBZ(1:100)=sin(x).*SSIBZ(1:100);
    SILBE(1:100)=sin(x).*SILBE(1:100);
    SILBN(1:100)=sin(x).*SILBN(1:100);
%   SILBZ(1:100)=sin(x).*SILBZ(1:100);
    x=flipud(x);
    SSIBE(tracelenSSIB-99:tracelenSSIB)=sin(x).*SSIBE(tracelenSSIB-99:tracelenSSIB);
    SSIBN(tracelenSSIB-99:tracelenSSIB)=sin(x).*SSIBN(tracelenSSIB-99:tracelenSSIB);
%   SSIBZ(tracelenSSIB-99:tracelenSSIB)=sin(x).*SSIBZ(tracelenSSIB-99:tracelenSSIB);
    SILBE(tracelenSILB-99:tracelenSSIB)=sin(x).*SILBE(tracelenSILB-99:tracelenSILB);
    SILBN(tracelenSILB-99:tracelenSSIB)=sin(x).*SILBN(tracelenSILB-99:tracelenSILB);
%   SILBZ(tracelenSILB-99:tracelenSSIB)=sin(x).*SILBZ(tracelenSILB-99:tracelenSILB);

    %Filter data:
    npo=2;
    npa=1;
    [PGCEf]=1.6e-4*bandpass(PGCE,40,lo,hi,npo,npa,'butter');
    [PGCNf]=1.6e-4*bandpass(PGCN,40,lo,hi,npo,npa,'butter');
%   [PGCZf]=1.6e-4*bandpass(PGCZ,40,lo,hi,npo,npa,'butter');
    clear PGCE
    clear PGCN
%   clear PGCZ
    if year==2003 && jday<213
        [SSIBEf]=20.0e-3*bandpass(SSIBE,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
        [SSIBNf]=20.0e-3*bandpass(SSIBN,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
%       [SSIBZf]=20.0e-3*bandpass(SSIBZ,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
        [SILBEf]=20.0e-3*bandpass(SILBE,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0* 
        [SILBNf]=20.0e-3*bandpass(SILBN,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
%       [SILBZf]=20.0e-3*bandpass(SILBZ,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
        [SSIBEbb]=20.0e-3*bandpass(SSIBE,100,lolo,hihi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
        [SSIBNbb]=20.0e-3*bandpass(SSIBN,100,lolo,hihi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
        [SILBEbb]=20.0e-3*bandpass(SILBE,100,lolo,hihi,npo,npa,'butter'); %times 5 until Aug 2003 5.0* 
        [SILBNbb]=20.0e-3*bandpass(SILBN,100,lolo,hihi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
    else
        [SSIBEf]=4.0e-3*bandpass(SSIBE,100,lo,hi,npo,npa,'butter'); 
        [SSIBNf]=4.0e-3*bandpass(SSIBN,100,lo,hi,npo,npa,'butter'); 
%       [SSIBZf]=4.0e-3*bandpass(SSIBZ,100,lo,hi,npo,npa,'butter'); 
        [SILBEf]=4.0e-3*bandpass(SILBE,100,lo,hi,npo,npa,'butter'); 
        [SILBNf]=4.0e-3*bandpass(SILBN,100,lo,hi,npo,npa,'butter'); 
%       [SILBZf]=4.0e-3*bandpass(SILBZ,100,lo,hi,npo,npa,'butter'); 
        [SSIBEbb]=4.0e-3*bandpass(SSIBE,100,lolo,hihi,npo,npa,'butter'); 
        [SSIBNbb]=4.0e-3*bandpass(SSIBN,100,lolo,hihi,npo,npa,'butter'); 
        [SILBEbb]=4.0e-3*bandpass(SILBE,100,lolo,hihi,npo,npa,'butter'); 
        [SILBNbb]=4.0e-3*bandpass(SILBN,100,lolo,hihi,npo,npa,'butter'); 
    end
    clear SSIBE
    clear SSIBN
%   clear SSIBZ
    clear SILBE
    clear SILBN
%   clear SILBZ

    SSIBEfd = resample(SSIBEf,2,5);
    SSIBNfd = resample(SSIBNf,2,5);
%   SSIBZfd = resample(SSIBZf,2,5);
    SILBEfd = resample(SILBEf,2,5);
    SILBNfd = resample(SILBNf,2,5);
%   SILBZfd = resample(SILBZf,2,5);
    SSIBEbbd = resample(SSIBEbb,2,5);
    SSIBNbbd = resample(SSIBNbb,2,5);
    SILBEbbd = resample(SILBEbb,2,5);
    SILBNbbd = resample(SILBNbb,2,5);
    clear SSIBEf
    clear SSIBNf
%   clear SSIBZf
    clear SILBEf
    clear SILBNf
%   clear SILBZf
    clear SSIBEbb
    clear SSIBNbb
    clear SILBEbb
    clear SILBNbb
 
    PGC=PGCEf+1i*PGCNf;
    SSIB=SSIBEfd+1i*SSIBNfd;
    SILB=SILBEfd+1i*SILBNfd;
    SSIBbb=SSIBEbbd+1i*SSIBNbbd;
    SILBbb=SILBEbbd+1i*SILBNbbd;
    PGCrot=PGC*exp(1i*rotPGC);
    SSIBrot=SSIB*exp(1i*rotSSIB);
    SILBrot=SILB*exp(1i*rotSILB);
    SSIBbbrot=SSIBbb*exp(1i*rotSSIB);
    SILBbbrot=SILBbb*exp(1i*rotSILB);
    clear PGCEf
    clear PGCEf
    clear SSIBEfd
    clear SSIBNfd
    clear SILBEfd
    clear SILBNfd
    clear SSIBEbbd
    clear SSIBNbbd
    clear SILBEbbd
    clear SILBNbbd
    clear PGC
    clear SSIB
    clear SILB

    if SSIBtoff > 0
        SSIBrot(1:tracelenPGC-SSIBsoff)=SSIBrot(SSIBsoff+1:tracelenPGC);
        SSIBrot(tracelenPGC-SSIBsoff+1:tracelenPGC)=0;
        SSIBbbrot(1:tracelenPGC-SSIBsoff)=SSIBbbrot(SSIBsoff+1:tracelenPGC);
        SSIBbbrot(tracelenPGC-SSIBsoff+1:tracelenPGC)=0;
    else
        SSIBrot(-SSIBsoff+1:tracelenPGC)=SSIBrot(1:tracelenPGC+SSIBsoff);
        SSIBrot(1:-SSIBsoff)=0;
        SSIBbbrot(-SSIBsoff+1:tracelenPGC)=SSIBbbrot(1:tracelenPGC+SSIBsoff);
        SSIBbbrot(1:-SSIBsoff)=0;
    end
    if SILBtoff > 0
        SILBrot(1:tracelenPGC-SILBsoff)=SILBrot(SILBsoff+1:tracelenPGC);
        SILBrot(tracelenPGC-SILBsoff+1:tracelenPGC)=0;
        SILBbbrot(1:tracelenPGC-SILBsoff)=SILBbbrot(SILBsoff+1:tracelenPGC);
        SILBbbrot(tracelenPGC-SILBsoff+1:tracelenPGC)=0;
    else
        SILBrot(-SILBsoff+1:tracelenPGC)=SILBrot(1:tracelenPGC+SILBsoff);
        SILBrot(1:-SILBsoff)=0;
        SILBbbrot(-SILBsoff+1:tracelenPGC)=SILBbbrot(1:tracelenPGC+SILBsoff);
        SILBbbrot(1:-SILBsoff)=0;
    end

    realPGC=real(PGCrot);
    realSSIB=real(SSIBrot);
    realSILB=real(SILBrot);
    realSSIBbb=real(SSIBbbrot);
    realSILBbb=real(SILBbbrot);

    % % Not efficient but easy book-keeping:
    % SSIBauto=realSSIB.*realSSIB;
    % SSIB2=cumsum(SSIBauto);
    % SILBauto=realSILB.*realSILB;
    % SILB2=cumsum(SILBauto);

    mlag=5;
    lenb=length(bostsec);
    imaxx=zeros(1,lenb);
    for n=1:lenb
        if bostend(n) <= length(realSILB)
            SILBtr=realSILB(boststart(n):bostend(n)); %do this over the narrower bandpass
            SSIBtr=realSSIB(boststart(n):bostend(n)); %do this over the narrower bandpass
	    SISScc=xcorr(SILBtr(winlen/2-20:winlen/2+20),SSIBtr(winlen/2-20:winlen/2+20),mlag);  %correlate over a 1-sec sub-window
	    [~,imaxx(n)]=max(SISScc);
	    imaxx(n)=imaxx(n)-(mlag+1); %This makes a logical lag of zero equal to zero
%             figure
%             hold on
% 	          plot(realSILB(boststart(n):bostend(n)),'b')
% 	          plot(realSSIB(boststart(n)-imaxx:bostend(n)-imaxx),'k')
%             title(int2str(bostsec(n)))
            SISSx=realSILB(boststart(n):bostend(n)).*realSSIB(boststart(n)-imaxx(n):bostend(n)-imaxx(n));
            SISSx(SISSx<0)=0;
            cumSISSx=cumsum(SISSx);
            SILBauto=realSILB(boststart(n):bostend(n)).*realSILB(boststart(n):bostend(n));
            SSIBauto=realSSIB(boststart(n)-imaxx(n):bostend(n)-imaxx(n)).*realSSIB(boststart(n)-imaxx(n):bostend(n)-imaxx(n));
            SILB2=cumsum(SILBauto);
            SSIB2=cumsum(SSIBauto);
            bostocks(n,7:12)=[cumSISSx(winlen/2-12) cumSISSx(winlen/2+18)-cumSISSx(winlen/2-12) cumSISSx(end)-cumSISSx(winlen/2+18) ...
                              SILB2(winlen/2-12)+SSIB2(winlen/2-12) SILB2(winlen/2+18)-SILB2(winlen/2-12)+SSIB2(winlen/2+18)-SSIB2(winlen/2-12) SILB2(end)-SILB2(winlen/2+18)+SSIB2(end)-SSIB2(winlen/2+18)];
            %The central window measuring strength of the arrival extends from -12 to +18; 0.75 sec
            %the above over the narrower bandpass
        end
    end
    percohere=bostocks(:,8)./(bostocks(:,7)+bostocks(:,8)+bostocks(:,9));
    percamp=bostocks(:,11)./(bostocks(:,10)+bostocks(:,11)+bostocks(:,12));
    fid = fopen([bostname(nd,:),num2str(lo),'-',num2str(hi),'.En'],'w'); %again, over the narrower bandpass
    fprintf(fid,'%3i %7i %3i %9.3f %6.3f %3i %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n',bostocks(:,:)');
    fclose(fid);
%     subplot(2,1,1,'align')
%     plot(bostocks(:,5),bostocks(:,8)./(bostocks(:,7)+bostocks(:,8)+bostocks(:,9)),'bo','markersize',2)
%     hold on
%     subplot(2,1,2,'align')
%     plot(bostocks(:,5),bostocks(:,11)./(bostocks(:,10)+bostocks(:,11)+bostocks(:,12)),'ro','markersize',2)
%     hold on
    
    nt=0;
    nrow=4;
    mcol=6;
    for ifig=1:floor(lenb/(nrow*mcol))+1
        figure('Position',[scrsz(3)/10 scrsz(4)/10 4*scrsz(3)/5 9*scrsz(4)/10]);
        for n = 1:nrow
            for m = 1:mcol
                nt=nt+1;
                if nt <= lenb
                    subplot(3*nrow,mcol,3*(n-1)*mcol+m,'align');
                    hold on
                    plot(realSILB(boststart(nt):bostend(nt)),'b')
                    plot(realSSIB(boststart(nt)-imaxx(nt):bostend(nt)-imaxx(nt)),'r')
                    yma=max(max([realSILB(boststart(nt):bostend(nt)) realSSIB(boststart(nt)-imaxx(nt):bostend(nt)-imaxx(nt))]));
                    ymi=min(min([realSILB(boststart(nt):bostend(nt)) realSSIB(boststart(nt)-imaxx(nt):bostend(nt)-imaxx(nt))]));
                    plot([winlen/2-10 winlen/2+10],[1.2*yma 1.2*yma],'k','linewidth',3)
                    xvect=[0 winlen*(yma-ymi)/2]; %amplitude bar, all the way across for 2
                    yma=2.*max(yma,-ymi);
                    yvect=[-0.9*yma -0.9*yma];
                    plot(xvect,yvect,'r','linewidth',3)
                    plot([40/hi 40/lo],[-0.8*yma -0.8*yma],'k','linewidth',3)
                    if n == 1 && m == 1
                        title(['002-',YEAR,'.',JDAY])
                    end
                    text(3*winlen/4, -0.8*yma, int2str(bostsec(nt)),'fontsize',6);
                    text(0.05*winlen, 0.75*yma, num2str(percohere(nt)),'fontsize',6);
                    text(3.*winlen/4, 0.75*yma, num2str(percamp(nt)),'fontsize',6);
                    box on
                    axis([0 winlen -yma yma])
                    set(gca,'XTick',[0 40 80 120 160],'fontsize',6);
                    
                    subplot(3*nrow,mcol,3*(n-1)*mcol+mcol+m,'align');
                    hold on
                    plot(realSILBbb(boststart(nt):bostend(nt)),'b')
                    plot(realSSIBbb(boststart(nt)-imaxx(nt):bostend(nt)-imaxx(nt)),'r')
                    yma=max(max([realSILBbb(boststart(nt):bostend(nt)) realSSIBbb(boststart(nt)-imaxx(nt):bostend(nt)-imaxx(nt))]));
                    ymi=min(min([realSILBbb(boststart(nt):bostend(nt)) realSSIBbb(boststart(nt)-imaxx(nt):bostend(nt)-imaxx(nt))]));
                    plot([winlen/2-10 winlen/2+10],[1.2*yma 1.2*yma],'k','linewidth',3)
                    xvect=[0 winlen*(yma-ymi)/2]; %amplitude bar, all the way across for 2
                    yma=2.*max(yma,-ymi);
                    yvect=[-0.9*yma -0.9*yma];
                    plot(xvect,yvect,'r','linewidth',3)
                    plot([40/hihi 40/lolo],[-0.8*yma -0.8*yma],'k','linewidth',3)
                    text(3*winlen/4, -0.8*yma, int2str(bostsec(nt)),'fontsize',6);
                    box on
                    axis([0 winlen -yma yma])
                    set(gca,'XTick',[0 40 80 120 160],'fontsize',6);

                    subplot(3*nrow,mcol,3*(n-1)*mcol+2*mcol+m,'align');
%                     [SILBxx f] = pwelch(realSILBbb(boststart(nt):bostend(nt))-mean(realSILBbb(boststart(nt):bostend(nt))),[],[],[],40);
%                     [SSIBxx f] = pwelch(realSSIBbb(boststart(nt)-imaxx(nt):bostend(nt)-imaxx(nt))-mean(realSSIBbb(boststart(nt)-imaxx(nt):bostend(nt)-imaxx(nt))),[],[],[],40);
                    [SILBxx f] = pwelch(detrend(realSILBbb(boststart(nt):bostend(nt))),[],[],[],40);
                    [SSIBxx f] = pwelch(detrend(realSSIBbb(boststart(nt)-imaxx(nt):bostend(nt)-imaxx(nt))),[],[],[],40);
                    loglog(f,SILBxx,'b+')
                    hold on
                    loglog(f,SSIBxx,'r+')
                    yma=max(max([SILBxx SSIBxx]));
                    ymi=0.00316*yma;
                    loglog(f,(SILBxx./SSIBxx)*ymi,'k')
                    [Csiss,f] = mscohere(detrend(realSILBbb(boststart(nt):bostend(nt))),detrend(realSSIBbb(boststart(nt)-imaxx(nt):bostend(nt)-imaxx(nt))),[],[],[],40);
                    semilogx(f,(yma-ymi)*Csiss+ymi,'g+')
                    axis([lolo hihi ymi yma])
                    box on
                    set(gca,'XTick',([0.5 1 2 3 4 5 6 7 8]),'fontsize',6);
                    set(gca,'YTick',([0.01*yma yma]),'fontsize',6);
                end
            end
        end
        orient landscape
        if ifig <= 9
            print('-depsc',['BWIGS/','002-',YEAR,'.',JDAY,'-',num2str(lolo),'-',int2str(hihi),'.',int2str(0),int2str(0),int2str(ifig),'.eps'])
        elseif ifig <= 99
            print('-depsc',['BWIGS/','002-',YEAR,'.',JDAY,'-',num2str(lolo),'-',int2str(hihi),'.',int2str(0),int2str(ifig),'.eps'])
        else
            print('-depsc',['BWIGS/','002-',YEAR,'.',JDAY,'-',num2str(lolo),'-',int2str(hihi),'.',int2str(ifig),'.eps'])
        end
    end
end
% subplot(2,1,1,'align')
% ylim([-0.5 1.5])
% ylabel('coherent energy 0.6s/4s')
% title('2004 - 2005')
% subplot(2,1,2,'align')
% ylim([0. 1.])
% ylabel('amp. squared 0.6s/4s')
% xlabel('Bostock magnitude')
% set(h,'PaperPosition',[0.25 0.25 8 10.5])
% print('-depsc','2004-2005BE.eps')
% %print(h,'-depsc',['temps_',int2str(timoffrot(1,1)),'-',int2str(timoffrot(size(timoffrot,1),1)),'_',num2str(lo),'-',int2str(hi),'Hz.eps'])
