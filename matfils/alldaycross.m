%Reads in seismograms; filters and rotates them, etc.  Plots SILB and SSIB on top of one another for specified 
%windows during one day.  For timeshifts up to +/- 5 samples gets most of central "hot spot". 
%Plots running normalized x-correlation for each of these offsets.  One printed page gets 30 seconds of 
%activity and all +/-5 samples of traces (11).
format short e
clear all
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

%    timoffrot(1,:)=[2003 062 07 42 12 +086 +020  80 115  50]; %SSIB 120 for 2-8 Hz
%    timoffrot(1,:)=[2003 063 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2003 064 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2003 065 07 42 12 +086 +020  80 115  50];
%   
%    timoffrot(1,:)=[2004 196 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2004 197 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2004 198 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2004 199 07 42 12 +086 +020  80 115  50];
%    
%    timoffrot(1,:)=[2005 254 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2005 255 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2005 256 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2005 257 07 42 12 +086 +020  80 115  50];

  
   %bostocks=load('BOSTOCK/NEW/002-246_2003.062');
   %bostocks=load('BOSTOCK/NEW/002-246_2003.063');
   %bostocks=load('BOSTOCK/NEW/002-246_2003.064');
   %bostocks=load('BOSTOCK/NEW/002-246_2004.196');
   %bostocks=load('BOSTOCK/NEW/002-246_2004.197');
   %bostocks=load('BOSTOCK/NEW/002-246_2004.198');
   %bostocks=load('BOSTOCK/NEW/002-246_2004.199');
   %bostocks=load('BOSTOCK/NEW/002-246_2005.254');
   %bostocks=load('BOSTOCK/NEW/002-246_2005.255');
   %bostocks=load('BOSTOCK/NEW/002-246_2005.256');
   %bostocks=load('BOSTOCK/NEW/002-246_2005.257');
   %bostsec=3600*(bostocks(:,3)-1)+bostocks(:,4)+4*0.025; %(the latter to correct for the diff between Bostock and Rubin filtering.  But why isn't it 8 samples?)

days=[2003.062;
      2003.063;
      2003.064;
      2004.196;
      2004.197;
      2004.198;
      2004.199;
      2005.254;
      2005.255;
      2005.256;
      2005.257];
ndays=length(days);
  
bads=[62200 62450; 
      0.025 0.025;
      0.025 0.025;
      2650 2675;
      43650 44100;
      0.025 0.025;
      0.025 0.025;
      0.025 0.025;
      0.025 0.025;
      18875 18925;
      70620 70645];
bads=bads*40;

infiles=['ARMMAP/MAPS/map2003.062.86.20.80.120.50_2-8-ms12-4s.dx.onepeak0.5.inbox';
         'ARMMAP/MAPS/map2003.063.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox';
         'ARMMAP/MAPS/map2003.064.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox';
         'ARMMAP/MAPS/map2004.196.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox'; 
         'ARMMAP/MAPS/map2004.197.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox';
         'ARMMAP/MAPS/map2004.198.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox';
         'ARMMAP/MAPS/map2004.199.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox';
         'ARMMAP/MAPS/map2005.254.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox';
         'ARMMAP/MAPS/map2005.255.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox';
         'ARMMAP/MAPS/map2005.256.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox';
         'ARMMAP/MAPS/map2005.257.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox'];

outfiles=['ARMMAP/MAPS/map2003.062.86.20.80.120.50_2-8-ms12-4s.dx.onepeak0.5';
         'ARMMAP/MAPS/map2003.063.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5';
         'ARMMAP/MAPS/map2003.064.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5';
         'ARMMAP/MAPS/map2004.196.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5'; 
         'ARMMAP/MAPS/map2004.197.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5';
         'ARMMAP/MAPS/map2004.198.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5';
         'ARMMAP/MAPS/map2004.199.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5';
         'ARMMAP/MAPS/map2005.254.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5';
         'ARMMAP/MAPS/map2005.255.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5';
         'ARMMAP/MAPS/map2005.256.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5';
         'ARMMAP/MAPS/map2005.257.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5'];

hi=6.;
lo=1.5;
npo=2;
npa=1;

rotSSIB=pi*(115-90)/180;
rotSILB=pi*(50-90)/180;
SSIBsoff=86;
SILBsoff=20;
SSIBtoff=SSIBsoff/40;
SILBtoff=SILBsoff/40;

cushion=0; %20; %5; %where, if the max cc coefficient is for that offset, you don't add the dot product.
nshift=0; %6+cushion; %max (-) shift of SSIB relative to SILB (NW corner)
mshift=0; %5+cushion; %max shift of SSIB relative to SILB (SE edge)
%cclen=floor(40/lo)+mod(floor(40/lo),2) %length of running cc window

scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite;

tottotcross=0.;
for iday=1:ndays
    year=floor(days(iday));
    YEAR=int2str(year);
    jday=round(1000*(days(iday)-year))
    if jday <= 9
        JDAY=['00',int2str(jday)];
    elseif jday<= 99
        JDAY=['0',int2str(jday)];
    else
        JDAY=int2str(jday);
    end
    MO=day2month(jday,year);
    %IDENTIF=[YEAR,'.',JDAY,'.',int2str(timoffrot(nd,6)),'.',int2str(timoffrot(nd,7)), ...
    %   '.',int2str(timoffrot(nd,8)),'.',int2str(timoffrot(nd,9)),'.',int2str(timoffrot(nd,10))]
    %Read data:
    direc=[YEAR,'/',MO,'/'];
    prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];
    SSIBEdat=[prename,'.SSIB..HHE.D.SAC'];
    SSIBNdat=[prename,'.SSIB..HHN.D.SAC'];
    SILBEdat=[prename,'.SILB..HHE.D.SAC'];
    SILBNdat=[prename,'.SILB..HHN.D.SAC'];

    [SSIBE,HdrDataSSIB,tnuSSIB,pobjSSIB,timsSSIB]=readsac(SSIBEdat,0,'l');
    [SSIBN,~,~,~,~]=readsac(SSIBNdat,0,'l');
    [SILBE,HdrDataSILB,tnuSILB,pobjSILB,timsSILB]=readsac(SILBEdat,0,'l');
    [SILBN,~,~,~,~]=readsac(SILBNdat,0,'l');

    tracelenSSIB=length(SSIBE);
    tracelenSILB=length(SILBE);

    %cosine taper before filtering:
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

    tracelenSSIB=length(SSIBEfd);
    tracelenSILB=length(SILBEfd);

    SSIB=SSIBEfd+1i*SSIBNfd;
    SILB=SILBEfd+1i*SILBNfd;
    SSIBrot=SSIB*exp(1i*rotSSIB);
    SILBrot=SILB*exp(1i*rotSILB);

    if SSIBtoff > 0
        SSIBrot(1:tracelenSSIB-SSIBsoff)=SSIBrot(SSIBsoff+1:tracelenSSIB);
        SSIBrot(tracelenSSIB-SSIBsoff+1:tracelenSSIB)=0;
    else
        SSIBrot(-SSIBsoff+1:tracelenSSIB)=SSIBrot(1:tracelenSSIB+SSIBsoff);
        SSIBrot(1:-SSIBsoff)=0;
    end
    if SILBtoff > 0
        SILBrot(1:tracelenSILB-SILBsoff)=SILBrot(SILBsoff+1:tracelenSILB);
        SILBrot(tracelenSILB-SILBsoff+1:tracelenSILB)=0;
    else
        SILBrot(-SILBsoff+1:tracelenSILB)=SILBrot(1:tracelenSILB+SILBsoff);
        SILBrot(1:-SILBsoff)=0;
    end

    realSSIB=real(SSIBrot);
    realSILB=real(SILBrot);
    small=1.e-8;
    realSSIB(bads(iday,1):bads(iday,2))=small;
    realSILB(bads(iday,1):bads(iday,2))=small;
    
    rubins=load(infiles(iday,:));
    rubouts=load(outfiles(iday,:));
    LIA=ismember(rubouts(:,1),rubins(:,1));
    rubouts(LIA,:)=[];

    offsetstart=100;
    detects=windows(days(iday));
    %wins=outsidewindows(days(iday));
    wins=[offsetstart 86400-offsetstart];
    nwins=size(wins,1);
    shortwin=30*40;
    totlen=0;
    totcross=0;
    ends=0*(1:nwins);
    for iwin=1:nwins 
        istart=wins(iwin,1)*40;
        iend=wins(iwin,2)*40;
        winlen=iend-istart+1;
        SILBauto=realSILB(istart:iend).*realSILB(istart:iend);
        SILB2=cumsum(SILBauto);
        %SSIBauto=realSSIB(istart-nshift:iend+mshift).*realSSIB(istart-nshift:iend+mshift); %Longer b.c. of mshift.
        SSIBauto=realSSIB(istart:iend).*realSSIB(istart:iend); %Not as precise, but should be OK
        SSIB2=cumsum(SSIBauto);
        SISSx=zeros(winlen,nshift+mshift+1);
        for k=-nshift:mshift;
            SISSx(:,k+nshift+1)=realSILB(istart:iend).*realSSIB(istart+k:iend+k);
        end
        if nshift+mshift==0
            cumsumSISS=cumsum(SISSx,1); 
        else
            cumsumSISS=cumsum(SISSx); 
        end
        ends(iwin)=totlen+winlen;
        cumsquared(totlen+1:totlen+winlen)=0.5*(SILB2+SSIB2);
        SISSxamalg=0*(1:winlen);
        ie=0; %need to initialize for when ie doesn't change (winlen<shortlen)
        for j=1:floor(winlen/shortwin)
            is=1+(j-1)*shortwin;
            ie=is+shortwin-1;
            [~,kmax]=max(cumsumSISS(ie,:)-cumsumSISS(is,:)); %find the offset with the max cc value
            if kmax<=cushion || kmax>nshift+mshift-cushion+1
                SISSxamalg(is:ie)=0;
            else
                SISSxamalg(is:ie)=SISSx(is:ie,kmax);
            end
        end %now add the remainder
        is=ie+1;
        ie=winlen;
        [~,kmax]=max(cumsumSISS(ie,:)-cumsumSISS(is,:)); %find the offset with the max cc value
            if kmax<=cushion || kmax>nshift+mshift-cushion+1
                SISSxamalg(is:ie)=0;
            else
                SISSxamalg(is:ie)=SISSx(is:ie,kmax);
            end
        cumsumSISS=cumsum(SISSxamalg);
        cumcross(totlen+1:totlen+winlen)=cumsumSISS;
        totlen=totlen+winlen;
        totcross=totcross+cumcross(end);
    end
    tottotcross=tottotcross+totcross;
    h=figure('Position',[wid/3 1 2.5*wid hite]); %center
    xlims=[0 21600;
           21600 43200;
           43200 64800;
           64800 86400];
    for ipan=1:4
        subplot(4,1,ipan,'align')
        hold on
        plot(offsetstart:0.025:86400-offsetstart,cumsquared/cumsquared(end),'k')
        hrf = plotreflinesr(gca,rubouts,'x','k');
        hrf = plotreflinesr(gca,rubins,'x','r');
        plot(offsetstart:0.025:86400-offsetstart,cumcross/cumcross(end),'r')
        hrf = plotreflinesr(gca,detects(:,1),'x','g');
        hrf = plotreflinesr(gca,detects(:,2),'x','b');
        if ipan==1
            title([YEAR,'.',JDAY,'  ',num2str(lo),'-',num2str(hi),' Hz, jiggering ', ...
                   int2str(-nshift+cushion),':',int2str(mshift-cushion),'  cushion = ',int2str(cushion)])
        end
        if ipan==4
            text(xlims(ipan,1)+2000,-0.06,num2str(cumsquared(end)))
            text(xlims(ipan,1)+5000,-0.06,num2str(cumcross(end)))
        end
        xlim(xlims(ipan,:))
        ylim([-0.11 1.1])
        box on
    end
    clear cumsquared
    clear cumcross
    clear ends
    text(83000,-0.06,num2str(tottotcross))
    set(h,'PaperPosition',[0.25 0.25 8 10.5])
    print(h,'-depsc',['BFIGS/',YEAR,'.',JDAY,'_',num2str(lo),'-',num2str(hi),'Hz', ...
        int2str(-nshift+cushion),':',int2str(mshift-cushion),'cush',int2str(cushion),'.eps'])
end
tottotcross
    
