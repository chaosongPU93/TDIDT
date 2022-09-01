%Reads in seismograms; filters and rotates them, etc.  Plots SILB and SSIB on top of one another for specified 
%windows during one day.  For timeshifts up to +/- 5 samples gets most of central "hot spot". 
%Plots running normalized x-correlation for each of these offsets.  One printed page gets 30 seconds of 
%activity and all +/-5 samples of traces (11).
format short e
clear all %UNCOMMENT
close all %UNCOMMENT
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

%    timoffrot(1,:)=[2003 062 07 42 12 +086 +020  80 120  50]; %120 for 2-8 Hz
     timoffrot(1,:)=[2003 063 07 42 12 +086 +020  80 115  50];
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

%  426.0    5.25    0.75  0.405  1.20     0.163     0.174   424.350     0.619    -0.911     1.467
% 1724.0    8.00   11.25  0.423 -0.16     0.371     1.861  1724.950     0.772     4.274     2.398
% 3320.0    8.75    1.75  0.400  1.32     0.167     0.466  3320.400     0.477    -0.995     2.679
% 3340.0   -0.25    1.50  0.422  0.41     0.550     0.647  3338.675     0.545     1.133    -2.131
% .dx files; output from wiggledaywig3 etc.  Columns 2 and 3 are xmaxPGSIntmp, xmaxPGSSntmp.
   %rubins=load('ARMMAP/MAPS/map2003.062.86.20.80.120.50_2-8-ms12-4s.dx.onepk.5');
    rubins=load('ARMMAP/MAPS/map2003.063.86.20.80.115.50_2-6-ms12-4s.dx.onepk.5');
   %rubins=load('ARMMAP/MAPS/map2003.064.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5');
   %rubins=load('ARMMAP/MAPS/map2004.196.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5'); %previously .inbox; missing red bars top/bottom.
   %rubins=load('ARMMAP/MAPS/map2004.199.86.20.80.115.50_2-6-ms19-4s.dx.onepk.5');
   %rubins=load('ARMMAP/MAPS/map2005.254.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5');
   %rubins=load('ARMMAP/MAPS/map2005.255.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5');
   %rubins=load('ARMMAP/MAPS/map2005.256.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5');
   %rubins=load('ARMMAP/MAPS/map2005.257.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5');
    
%246  30304  1 3284.125 1.632 10
%  2  30304  2 3108.825 1.681 13
%  2  30304  2 3115.600 1.727 15
%  2  30304  2 3127.925 1.697 17
   %bostocks=load('BOSTOCK/NEW/002-246_2003.062');
    bostocks=load('BOSTOCK/NEW/002-246_2003.063');
   %bostocks=load('BOSTOCK/NEW/002-246_2003.064');
   %bostocks=load('BOSTOCK/NEW/002-246_2004.196');
   %bostocks=load('BOSTOCK/NEW/002-246_2004.199');
   %bostocks=load('BOSTOCK/NEW/002-246_2005.254');
   %bostocks=load('BOSTOCK/NEW/002-246_2005.255');
   %bostocks=load('BOSTOCK/NEW/002-246_2005.256');
   %bostocks=load('BOSTOCK/NEW/002-246_2005.257');

   bostsec=3600*(bostocks(:,3)-1)+bostocks(:,4)+4*0.025; %(the latter to correct for the diff between Bostock and Rubin filtering.  But why isn't it 8 samples?)

%  wins=[25475 25581; %2003.062 (postdates 30-s)
%        25680 25788;
%        26035 27065;
%        27170 27297;
%        27905 28005;
%        28150 28179;
%        28250 28279;
%        28450 28507;
%        28830 29051;
%        29405 29454;
%        29750 29922;
%        30120 30149;
%        30285 30473;
%        30820 31101;
%        31190 31243;
%        31700 31856;
%        32505 32739;
%        34170 34408;
%        37580 37712; 
%        38415 38789;
%        40585 40860;
%        41315 41509;
%        43290 43543;
%        46340 46618;
%        49325 49490;
%        52165 52413;
%        54915 55075;
%        61530 61757;
%        68420 68638;
%        68720 68972;
%        72285 72397;
%        75135 75792;
%        75890 75991;
%        77195 77690;
%        77970 77995;
%        82495 82714;
%        84840 84924;
%        85045 85094];
     wins=[6605 7073; %2003.063
          8830 9034;
          10340 10456;
          12440 12521;
          21795 22579;
          32350 33058;
          41375 41692;
          50200 50431;
          65000 65489;
          80425 80995];
%    wins=[51060 51554; %2003.064
%          51730 51814;
%          53485 53660];

%   wins=[39460 40060; %2004.196 (postdates 30-s)
%         40360 40540;
%         40940 41364;
%         41955 42275;
%         42910 43032;
%         43610 43755;
%         44220 44355;
%         46260 46575;
%         47080 47142;
%         48360 48773;];
%     wins=[51080 51500; %2004.196 (postdates 30-s)
%           56230 56624;
%           56845 57315;
%           62500 63148;
%           64935 65270;
%           70420 70684;
%           73940 74277;
%           78885 79124;
%           81520 81768;
%           85830 86199];
%     wins=[3040 3412; %2004.199 (postdates 30-s)
%           4210 4277;
%           4860 4916;
%           6175 6255;
%           8930 9201];

%    wins=[14640 15270; %2005.254 (postdates 30-s)
%          15740 15859; 
%          18090 18170; 
%          20035 20115; 
%          20560 20631; 
%          21155 21331; 
%          24335 24525; 
%          27430 27694;
%          30330 31065;
%          34560 34945;
%          35260 35903;
%          49140 50144;
%          54290 54719;
%          68770 69082;
%          73360 73765; 
%          79970 80340
% 	     82915 83015];
%    wins=[1000 1029; %2005.255 (postdates 30-s)
% 	     5305 5800;
% 	     8060 8149;
% 	     19000 19245;
% 	     22980 23030;
%          35680 36362;
%          51030 51544;
%          54205 54359;
%          58930 59423;
%          62305 62626;
%          68200 68544;
%          75540 75638];
%    wins=[3655 4115; %2005.256 (postdates 30-s)
%          81410 81468;
%          81760 81903;
%          82535 82597;
%          84655 84772]; 
%    wins=[2705 2848]; %2005.257 (postdates 30-s)

nd=size(timoffrot,1);

hi=6.;
lo=1.5;
hi=4.5;
lo=1.;
%lo=0.75;
% hi=8.;
% lo=0.5;
% hi=1.5;
% lo=0.5;

% hi=3.;
% lo=0.75;
% hi=6.;
% lo=0.75;
% hi=6.;
% lo=0.5;
% hi=2.;
% lo=0.5;
% hi=1.;
% lo=0.25;
% hi=1.5;
% lo=0.375;

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
imagSILB=imag(SILBrot);

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

nshift=5; %max (-) shift of SSIB relative to SILB (NW corner)
mshift=5; %max shift of SSIB relative to SILB (SE edge)
cclen=floor(40/lo)+mod(floor(40/lo),2) %length of running cc window
ndata=0; %UNCOMMENT
for n=1:size(wins,1)
    plotstart=wins(n,1)*40;
    seclen=(wins(n,2)-wins(n,1))+30-mod((wins(n,2)-wins(n,1)),30);
    plotend=plotstart+seclen*40;
    longlen=plotend-plotstart+1;
    SILBauto=realSILB(plotstart:plotend).*realSILB(plotstart:plotend);
    SILB2=cumsum(SILBauto);
    SSIBauto=realSSIB(plotstart-nshift:plotend+mshift).*realSSIB(plotstart-nshift:plotend+mshift); %Longer b.c. of mshift.
    SSIB2=cumsum(SSIBauto);
    SISSx=zeros(longlen,nshift+mshift+1);
    for k=-nshift:mshift;
        SISSx(:,k+nshift+1)=realSILB(plotstart:plotend).*realSSIB(plotstart+k:plotend+k);
    end
    cumsumSISS=cumsum(SISSx);
    SISSnum=cumsumSISS(cclen+1:longlen,:)-cumsumSISS(1:longlen-cclen,:); %has length longlen minus cclen (in samples); individually, think of as centered at cclen/2.
    SILBden=SILB2(cclen+1:longlen)-SILB2(1:longlen-cclen); %A long vector, cclen shorter than the long window
    SSIBden=SSIB2(cclen+1:longlen+nshift+mshift)-SSIB2(1:longlen+nshift+mshift-cclen); %A long vector, cclen shorter than the long window+2*mshift
    SISSn=zeros(longlen,nshift+mshift+1);  %Make as long as seismograms, but will remain zeros for first and last cclen/2 samples
    for k=-nshift:mshift;
        SISSn(floor(cclen/2)+1:longlen-floor(cclen/2),nshift+k+1)=SISSnum(:,nshift+k+1)./realsqrt(SILBden.*SSIBden(nshift+k+1:longlen+nshift+k-cclen));
    end
    runningcc=SISSn(:,nshift+1);  %a placeholder, for now
    runningSSamp=runningcc;  %a placeholder, for now
    runningSSfreq=runningcc;  %a placeholder, for now
    SISSn(SISSn<0)=-2; %So they don't plot.
    cumSISSn=cumsum(SISSn);
%     ltmax=max([realSSIB(plotstart:plotend); realSILB(plotstart:plotend)]);
%     ltmin=min([realSSIB(plotstart:plotend); realSILB(plotstart:plotend)]);

    inrub=1;
    while rubins(inrub,1)<wins(n,1) %If the center of the detection window is within the desired range
        inrub=inrub+1;
    end
    inrubstart=inrub;
    while rubins(inrub,1)<wins(n,2) && inrub < length(rubins)
        inrub=inrub+1;
    end
    inrubend=inrub;

    inbos=1;
    while bostsec(inbos)<wins(n,1) %If the center of the detection window is within the desired range
        inbos=inbos+1;
    end
    inbosstart=inbos+1;
    while bostsec(inbos)<wins(n,2) && inbos < length(bostsec)
        inbos=inbos+1;
    end
    inbosend=inbos;

    hSI=hilbert(realSILB(plotstart:plotend));
    SIamp=abs(hSI);
    hSIort=hilbert(imagSILB(plotstart:plotend));
    SIortamp=abs(hSIort);
%     SIfreq = 40*diff(unwrap(angle(hSI)))/(2*pi);
%     SIfreq(end+1)=SIfreq(end);
    xprime=0*(1:longlen);
    yprime=0*(1:longlen);
    xprime(2:longlen-1)=real(hSI(3:end))-real(hSI(1:end-2));
    xprime(1)=real(hSI(2))-real(hSI(1));
    xprime(longlen)=real(hSI(longlen))-real(hSI(longlen-1));
    yprime(2:longlen-1)=imag(hSI(3:end))-imag(hSI(1:end-2));
    yprime(1)=imag(hSI(2))-imag(hSI(1));
    yprime(longlen)=imag(hSI(longlen))-imag(hSI(longlen-1));
    SIfreq=((40/2)/(2*pi))*(realSILB(plotstart:plotend).*yprime'-xprime'.*imag(hSI))./abs(hSI).^2; %divided by 2 b.c. of centered f.d.
    %%%%%SISSx(:,k+nshift+1)=realSILB(plotstart:plotend).*realSSIB(plotstart+k:plotend+k); %Here for reference
    plotstartSS=plotstart-nshift; 
    plotendSS=plotend+mshift; 
    hSS=hilbert(realSSIB(plotstartSS:plotendSS)); %longer because of potential shifts
    SSlen=plotendSS-plotstartSS+1;
    SSamp=abs(hSS);
    xprime=0*(1:SSlen);
    yprime=0*(1:SSlen);
    xprime(2:SSlen-1)=real(hSS(3:end))-real(hSS(1:end-2));
    xprime(1)=real(hSS(2))-real(hSS(1));
    xprime(SSlen)=real(hSS(SSlen))-real(hSS(SSlen-1));
    yprime(2:SSlen-1)=imag(hSS(3:end))-imag(hSS(1:end-2));
    yprime(1)=imag(hSS(2))-imag(hSS(1));
    yprime(SSlen)=imag(hSS(SSlen))-imag(hSS(SSlen-1));
    SSfreq=((40/2)/(2*pi))*(realSSIB(plotstartSS:plotendSS).*yprime'-xprime'.*imag(hSS))./abs(hSS).^2; %divided by 2 b.c. of centered f.d.
    
    npages=seclen/30;
    for j=1:npages
        %h=figure('Position',[wid/3 1 2.5*wid hite]); %UNCOMMENT
        istart=plotstart+(j-1)*shortlen+1; %added +1 late
        iend=istart+shortlen;
        inormstart=istart-plotstart;
        inormend=inormstart+shortlen;
        iscmplx=istart-plotstart;
        iecmplx=iscmplx+shortlen;
        ltmax=max([realSSIB(istart-nshift:iend+mshift); realSILB(istart:iend)]);
        ltmin=min([realSSIB(istart-nshift:iend+mshift); realSILB(istart:iend)]);
        ltmax=1.5*max(ltmax,-ltmin);
        [dummy,kmax]=max(cumSISSn(iecmplx,:)-cumSISSn(iscmplx,:));
        runningcc(iscmplx:iecmplx)=SISSn(iscmplx:iecmplx,kmax); %make the running cc equal to that of the most coherent 30 sec.
        kmax=kmax-1-nshift %to center
        runningSSamp(iscmplx:iecmplx)=SSamp(iscmplx+nshift+kmax:iecmplx+nshift+kmax); %to account for shift of SSIB w.r.t SILB
        runningSSfreq(iscmplx:iecmplx)=SSfreq(iscmplx+nshift+kmax:iecmplx+nshift+kmax); %same as above
%         for k=-nshift:mshift
%             subplot(nshift+mshift+2,1,k+nshift+2,'align'); 
%             hold on
%             for inr=inrubstart:inrubend  %inefficient
%                 if k == -nshift
%                     if round(rubins(inr,2)-rubins(inr,3)) < k
%                         irstart=(rubins(inr,1)-2)*40-round(rubins(inr,2));
%                         irend=irstart+4*40-1;
%                         plot(timsPGC(irstart:irend),0.9*ltmax,'r','linewidth',3)
%                     end
%                 end
%                 if (k == mshift) && (round(rubins(inr,2)-rubins(inr,3)) > k)
%                     irstart=(rubins(inr,1)-2)*40-round(rubins(inr,2));
%                     irend=irstart+4*40-1;
%                     plot(timsPGC(irstart:irend),0.9*ltmax,'r','linewidth',3)
%                 end
%                 if round(rubins(inr,2)-rubins(inr,3))==k 
%                     irstart=(rubins(inr,1)-2)*40-round(rubins(inr,2));
%                     irend=irstart+4*40-1;
%                     plot(timsPGC(irstart:irend),realPGC(irstart+round(rubins(inr,2)):irend+round(rubins(inr,2))),'r')
%                     text(timsPGC(irstart),-0.9*ltmax,int2str(round(rubins(inr,2))),'fontsize',6);
%                 end
%             end
%             if k==kmax
%                 plot(timsPGC(istart:iend),SIamp(iscmplx:iecmplx),'g');
%             end
%             plot(timsPGC(istart:iend),SISSn(inormstart:inormend,k+nshift+1)*ltmax,'co','markersize',1)
%             if k==0
%                 for inb=inbosstart:inbosend
%                     plot(bostsec,0.85*ltmax,'ko','MarkerSize',4,'MarkerFaceColor','k')
%                 end
%             end
%             plot(timsPGC(istart:iend),realSILB(istart:iend),'k');
%             plot(timsPGC(istart:iend),realSSIB(istart+k:iend+k),'b');
%             axis([timsPGC(istart) timsPGC(iend) -ltmax ltmax]);
%             box on
%        end
%        subplot(nshift+mshift+2,1,1,'align'); 
%        hold on
%        %semilogy(timsPGC(istart:iend),SIfreq(iscmplx:iecmplx),'b'); %DON'T KNOW WHY THIS DOESN'T WORK! 
%        plot(timsPGC(istart:iend),SIfreq(iscmplx:iecmplx),'b');
%        plot(timsPGC(istart:iend),SIamp(iscmplx:iecmplx)./SIortamp(iscmplx:iecmplx),'r');
%        axis([timsPGC(istart) timsPGC(iend) 0.5 10]);
%        set(gca,'ytick',[1 3 5 10]);
%        title([IDENTIF,'  ',num2str(lo),'-',num2str(hi),' Hz'])
%        box on
%        
%        set(h,'PaperPosition',[0.25 0.25 8 10.5])
%        orient landscape
%        %if j <= 9
%        print(h,'-depsc',['ARMMAP/CMPWIGS/',YEAR,'.',JDAY,'.',int2str(timsPGC(istart)),'_',num2str(lo),'-',num2str(hi),'eye.eps'])
%        %else
%        %    print(h,'-depsc',['ARMMAP/WIGS/',YEAR,'.',JDAY,'.',int2str(wins(n,1)),'_',num2str(lo),'-',num2str(hi),'_',int2str(j),'eye.eps'])
%        %end
    end
    cmplxseisms(ndata+1:ndata+longlen,1)=SIamp;
    cmplxseisms(ndata+1:ndata+longlen,2)=SIfreq;
    cmplxseisms(ndata+1:ndata+longlen,3)=SIamp./SIortamp;
    cmplxseisms(ndata+1:ndata+longlen,4)=runningcc;
    cmplxseisms(ndata+1:ndata+longlen,5)=runningSSamp;
    cmplxseisms(ndata+1:ndata+longlen,6)=runningSSfreq;
    ndata=ndata+longlen;
end
aa=cmplxseisms;
tf1=cmplxseisms(:,1)<1.e-3;
tf2=cmplxseisms(:,1)>10;
tf3=cmplxseisms(:,2)<0.1;
tf4=cmplxseisms(:,2)>10;
tf5=cmplxseisms(:,3)<1.5;
TFall=tf1 | tf2 | tf3 | tf4 | tf5 ;
aa(TFall,:)=[];
aa2=aa;
tf6=aa(:,4)<0.5;
aa2(tf6,:)=[];
aa2=sortrows(aa2,-1);
%values=hist3(aa(:,1:2),[201 201]);
%xb = linspace(min(aa(:,1)),max(aa(:,1)),size(values,1)+1);
%yb = linspace(min(aa(:,2)),max(aa(:,2)),size(values,1)+1);
%imagesc(xb,yb,values')
figure
loglog(cmplxseisms(:,1),cmplxseisms(:,2),'bo','markersize',1)
hold on
loglog(aa(:,1),aa(:,2),'ro','markersize',1)
loglog(aa2(:,1),aa2(:,2),'go','markersize',1)
in=200;
navs=floor(size(aa2,1)/in);
x=0*(1:navs);
y=0*(1:navs);
for i=1:navs
    x(i)=median(aa2((i-1)*in+1:i*in,1));
    y(i)=median(aa2((i-1)*in+1:i*in,2));
end
loglog(x,y,'k','linewidth',2)
xlim([0.001 2])
ylim([0.1 10])
xlabel('log envelope amplitude')
ylabel('log frequency (Hz)')
title(['SILB  ',YEAR,'.',JDAY,':',num2str(lo),'-',num2str(hi),' Hz'])
print('-depsc',['Amp_v_FreqSILB_',YEAR,'.',JDAY,'_',num2str(lo),'-',num2str(hi),'.eps'])

aa2=sortrows(aa2,-5);
figure
loglog(cmplxseisms(:,5),cmplxseisms(:,6),'bo','markersize',1)
hold on
loglog(aa(:,5),aa(:,6),'ro','markersize',1)
loglog(aa2(:,5),aa2(:,6),'go','markersize',1)
in=200;
navs=floor(size(aa2,1)/in);
x=0*(1:navs);
y=0*(1:navs);
for i=1:navs
    x(i)=median(aa2((i-1)*in+1:i*in,5));
    y(i)=median(aa2((i-1)*in+1:i*in,6));
end
loglog(x,y,'k','linewidth',2)
xlim([0.001 2])
ylim([0.1 10])
xlabel('log envelope amplitude')
ylabel('log frequency (Hz)')
title(['SSIB  ',YEAR,'.',JDAY,':',num2str(lo),'-',num2str(hi),' Hz'])
print('-depsc',['Amp_v_FreqSILB_',YEAR,'.',JDAY,'_',num2str(lo),'-',num2str(hi),'.eps'])
