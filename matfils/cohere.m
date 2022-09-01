%Reads in seismograms; filters and rotates them, etc.  Can plot coherence 
%on optimal component, or amplitude spectrum on both optimal and orthogonal components,
%and their ratio, for SILB and SSIB on top of one another for specified 
%windows during one day.  
format short e
clear all
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

   %timoffrot(1,:)=[2003 062 07 42 12 +086 +020  80 120  50]; %120 for 2-8 Hz
   %timoffrot(1,:)=[2003 063 07 42 12 +086 +020  80 115  50];
   %timoffrot(1,:)=[2003 064 07 42 12 +086 +020  80 115  50];
   %timoffrot(1,:)=[2003 065 07 42 12 +086 +020  80 115  50];
%    
%    timoffrot(1,:)=[2005 254 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2005 255 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2005 256 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2005 257 07 42 12 +086 +020  80 115  50];
%   
     timoffrot(1,:)=[2004 196];
%    timoffrot(1,:)=[2004 197 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2004 198 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2004 199];

%    wins=[26045 27297; %2003.062
%          27916 28005;
%          28839 29051;
%          29757 29922;
%          38424 38789];   %wins=[26185 26275];
%    wins=[21800 22640; %2003.063
%          32300 32900;
%          32900 33080;
%          41382 41652;
%          41652 41832];
%    wins=[41473 41528]; %41503 41528]; %2003.063

%   wins=[39460 40060; %2004.196
%         40360 40540;
%         40940 41364;
%         41955 42275;
%         42910 43032;
%         43610 43755;
%         44220 44355;
%         46260 46575;
%         47080 47142;
%         48360 48773;
%         51080 51500; 
%         56230 56624;
%         56845 57315;
%         62500 63148;
%         64935 65270;
%         70420 70684;
%         73940 74277;
%         78885 79124;
%         81520 81768;
%         85830 86199];
%   wins=[42960 42990]; %2004.196
    wins=[51380 51399]; %2004.196
%     wins=[4245 4260; %2004.199
%           8968 8985;
%           9000 9015;
%           9055 9065];

%    wins=[14600 15320; %2005.254
%          21130 21400; 
%          27400 27580;
%          27580 27670;
%          27670 27730;
%          30300 30630;
%          30630 31080;
%          34500 34980;
%          35200 35710;
%          35710 35950;
%          49140 49920;
%          49920 50160;
%          54260 54770;
%          68760 69090;
%          73350 73800; 
%          79950 80370];
%      wins=[5310 5800; %2005.255
%            8064 8149;
%            19006 19245;
%            22984 23030;
%            35688 36362;
%            51033 51548;
%            54210 54359;
%            58935 59423;
%            62310 62626;
%            68205 68544;
%            75545 75638];

nd=size(timoffrot,1);

% hi=6.;
% lo=1.5;
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

PERMSTA=['PGC'
         'LZB'];
POLSTA=['SSIB'
        'SILB'
        'KLNB'
        'TWKB'
        'MGCB'];
PERMROTS=[0 90 32 00;  %PGC , Yajun's "0" changed to 90.
          0 90 54 00]; %LZB
 POLROTS=[6 85 33 86;  %SSIB
          0 90 39 20;  %SILB
          0 90  7 00;  %KLNB
          4 70 48 00;  %MGCB
          4 75 38 00]; %TWKB
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;
POLROTS(:,1)=round(POLROTS(:,1)*(100/40));

yr=timoffrot(nd,1);
YEAR=int2str(yr);
jday=timoffrot(nd,2);
if jday <= 9
    JDAY=['00',int2str(jday)];
elseif jday<= 99
    JDAY=['0',int2str(jday)];
else
    JDAY=int2str(jday);
end
MO=day2month(jday,yr);
IDENTIF=[YEAR,'.',JDAY]

direc=[YEAR,'/',MO,'/'];
prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];
scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite;
xero=[1.e-5 1;
      1.e5 1];
  
for ista=1:size(POLSTA,1)
    hi=30;
    lo=0.01;
    sps=100;
    STAEdat=[prename,'.',POLSTA(ista,:),'..HHE.D.SAC']; %BHE for permstas.
    STANdat=[prename,'.',POLSTA(ista,:),'..HHN.D.SAC'];
    [STAE,HdrDataSTA,tnuSTA,pobjSTA,timsSTA]=readsac(STAEdat,0,'l');
    [STAN,~,~,~,~]=readsac(STANdat,0,'l');
    tracelen=length(STAE);
    %cosine taper before filtering:
    x=(0:pi/200:pi/2-pi/200)';
    STAE(1:100)=sin(x).*STAE(1:100);
    STAN(1:100)=sin(x).*STAN(1:100);
    x=flipud(x);
    STAE(tracelen-99:tracelen)=sin(x).*STAE(tracelen-99:tracelen);
    STAN(tracelen-99:tracelen)=sin(x).*STAN(tracelen-99:tracelen);
    %Filter data:
    npo=2;
    npa=1;
    if yr==2003 && jday<213
        [STAEf]=20.0e-3*bandpass(STAE,sps,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
        [STANf]=20.0e-3*bandpass(STAN,sps,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
    else
        [STAEf]=4.0e-3*bandpass(STAE,sps,lo,hi,npo,npa,'butter'); 
        [STANf]=4.0e-3*bandpass(STAN,sps,lo,hi,npo,npa,'butter'); 
    end
    STA=STAEf+1i*STANf;
    STAfs=STA*exp(1i*POLROTS(ista,2));
    STAslow=-real(STAfs);
    STAfast=imag(STAfs);
    len=length(STA);
    STAslow(20:len-20)=STAslow(20+POLROTS(ista,1):len-20+POLROTS(ista,1));
    
    STAsc=STAslow+1i*STAfast;
    STAscrot=STAsc*exp(1i*POLROTS(ista,3));

%     SSIBEfd = resample(SSIBEf,2,5);
%     SSIBNfd = resample(SSIBNf,2,5);
%     SILBEfd = resample(SILBEf,2,5);
%     SILBNfd = resample(SILBNf,2,5);

    STAsoff=POLROTS(ista,4)
    STAtoff=round(STAsoff/40); %really 40.  That's how it's listed.

    if STAtoff > -1
        STAscrot(1:tracelen-STAsoff)=STAscrot(STAsoff+1:tracelen);
        STAscrot(tracelen-STAsoff+1:tracelen)=0;
    else
        STAscrot(-STAsoff+1:tracelen)=STAscrot(1:tracelen+STAsoff);
        STAscrot(1:-STAsoff)=0;
    end
    STAort=real(STAscrot);
    STAopt=imag(STAscrot);

    for n=1:size(wins,1)
        istart=wins(n,1)*sps;
        iend=wins(n,2)*sps;
        %[Csiss,f] = mscohere(realSILB(istart:iend),realSSIB(istart:iend),[],[],[],sps);   
        [STAoptxx f1] = pwelch(STAopt(istart:iend),[],[],[],sps);
        [STAortxx f1] = pwelch(STAort(istart:iend),[],[],[],sps);
        h=figure('Position',[wid/10 hite/3 2.5*wid hite/1.7]);
        loglog(f1,STAoptxx,'b')
        hold on
        loglog(f1,STAortxx,'r')
        loglog(f1,STAoptxx./STAortxx,'g')
        %loglog(f,10.^(4*(Csiss-1)),'c')
        loglog(xero(:,1),xero(:,2),'k')
        loglog(xero(:,1),10*xero(:,2),'k')
        loglog(xero(:,1),1.e-4*xero(:,2),'k')
        xlim([0.005 30])
        title([POLSTA(ista,:),'  ',YEAR,'.',JDAY,'  ',int2str(wins(n,1)),'-',int2str(wins(n,2)),' sec'])
        set(h,'PaperPosition',[0.25 5 8 5])
        print(h,'-depsc',['spectr',POLSTA(ista,:),YEAR,'.',JDAY,'_',int2str(wins(n,1)),'-',int2str(wins(n,2)),'.eps'])
    end
end

for ista=1:size(PERMSTA,1)
    hi=18;
    lo=0.01;
    sps=40;
    STAEdat=[prename,'.',PERMSTA(ista,:),'..BHE.D.SAC']; %BHE for permstas.
    STANdat=[prename,'.',PERMSTA(ista,:),'..BHN.D.SAC'];
    [STAE,HdrDataSTA,tnuSTA,pobjSTA,timsSTA]=readsac(STAEdat,0,'l');
    [STAN,~,~,~,~]=readsac(STANdat,0,'l');
    tracelen=length(STAE);
    %cosine taper before filtering:
    x=(0:pi/80:pi/2-pi/80)';
    %Seems to be necessary at the start of each day for PGCE:
    STAE(1:80)=0.;
    STAE(81:120)=sin(x).*STAE(81:120); %Only at start of day!
    STAN(1:40)=sin(x).*STAN(1:40);
    x=flipud(x);
    STAE(tracelen-39:tracelen)=sin(x).*STAE(tracelen-39:tracelen);
    STAN(tracelen-39:tracelen)=sin(x).*STAN(tracelen-39:tracelen);
    %Filter data:
    npo=2;
    npa=1;
    [STAEf]=1.6e-4*bandpass(STAE,sps,lo,hi,npo,npa,'butter');
    [STANf]=1.6e-4*bandpass(STAN,sps,lo,hi,npo,npa,'butter');

    STA=STAEf+1i*STANf;
    if strcmp(PERMSTA(ista,:),'LZB')
        STA=10*STA;
    end
    STAfs=STA*exp(1i*PERMROTS(ista,2));
    STAslow=-real(STAfs);
    STAfast=imag(STAfs);
    len=length(STA);
    STAslow(20:len-20)=STAslow(20+PERMROTS(ista,1):len-20+PERMROTS(ista,1));
    
    STAsc=STAslow+1i*STAfast;
    STAscrot=STAsc*exp(1i*PERMROTS(ista,3));

    STAsoff=PERMROTS(ista,4)
    STAtoff=round(STAsoff/40); 

    if STAtoff > -1
        STAscrot(1:tracelen-STAsoff)=STAscrot(STAsoff+1:tracelen);
        STAscrot(tracelen-STAsoff+1:tracelen)=0;
    else
        STAscrot(-STAsoff+1:tracelen)=STAscrot(1:tracelen+STAsoff);
        STAscrot(1:-STAsoff)=0;
    end
    STAort=real(STAscrot);
    STAopt=imag(STAscrot);

    for n=1:size(wins,1)
        istart=wins(n,1)*sps;
        iend=wins(n,2)*sps;
        %[Csiss,f] = mscohere(realSILB(istart:iend),realSSIB(istart:iend),[],[],[],sps);   
        [STAoptxx f1] = pwelch(STAopt(istart:iend),[],[],[],sps);
        [STAortxx f1] = pwelch(STAort(istart:iend),[],[],[],sps);
        h=figure('Position',[wid/10 hite/3 2.5*wid hite/1.7]);
        loglog(f1,STAoptxx,'b')
        hold on
        loglog(f1,STAortxx,'r')
        loglog(f1,STAoptxx./STAortxx,'g')
        %loglog(f,10.^(4*(Csiss-1)),'c')
        loglog(xero(:,1),xero(:,2),'k')
        loglog(xero(:,1),10*xero(:,2),'k')
        loglog(xero(:,1),1.e-4*xero(:,2),'k')
        xlim([0.005 18])
        title([PERMSTA(ista,:),'  ',YEAR,'.',JDAY,'  ',int2str(wins(n,1)),'-',int2str(wins(n,2)),' sec'])
        set(h,'PaperPosition',[0.25 5 8 5])
        print(h,'-depsc',['spectr',PERMSTA(ista,:),YEAR,'.',JDAY,'_',int2str(wins(n,1)),'-',int2str(wins(n,2)),'.eps'])
    end
end

    %Filter data:
%     npo=2;
%     npa=1;
