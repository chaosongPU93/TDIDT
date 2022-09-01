%Reads in seismograms and Bostock detections. Filters and rotates seisms, etc.  Stacks based on times of B. detections.
%At SILB, SSIB, PGC in user-specified passband.
format short e
clear all
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

%GET RID OF -22.675 for families other than 002!!!
% 002
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
% 047 (Bad stuff on day 197???)
%     timoffrot(1,:)=[2004 196 07 42 12 +086 +020  80 115  50];
%     timoffrot(2,:)=[2004 197 07 42 12 +086 +020  80 115  50];
%     timoffrot(3,:)=[2004 198 07 42 12 +086 +020  80 115  50];
%     timoffrot(4,:)=[2004 199 07 42 12 +086 +020  80 115  50];
%     timoffrot(5,:)=[2004 200 07 42 12 +086 +020  80 115  50];
%     timoffrot(6,:)=[2005 254 07 42 12 +086 +020  80 115  50];
%     timoffrot(7,:)=[2005 255 07 42 12 +086 +020  80 115  50];
%     timoffrot(8,:)=[2005 256 07 42 12 +086 +020  80 115  50];
%     timoffrot(9,:)=[2005 260 07 42 12 +086 +020  80 115  50];
% 068
%     timoffrot(1,:)=[2004 198 07 42 12 +086 +020  80 115  50];
%     timoffrot(2,:)=[2004 199 07 42 12 +086 +020  80 115  50];
%     timoffrot(3,:)=[2004 200 07 42 12 +086 +020  80 115  50];
%     timoffrot(4,:)=[2004 201 07 42 12 +086 +020  80 115  50];
%     timoffrot(5,:)=[2005 256 07 42 12 +086 +020  80 115  50];
%     timoffrot(6,:)=[2005 257 07 42 12 +086 +020  80 115  50];
%     timoffrot(7,:)=[2005 258 07 42 12 +086 +020  80 115  50];
%     timoffrot(8,:)=[2005 259 07 42 12 +086 +020  80 115  50];
%     timoffrot(9,:)=[2005 260 07 42 12 +086 +020  80 115  50];
%     timoffrot(10,:)=[2005 261 07 42 12 +086 +020  80 115  50];
    
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
%047 040714 7 1689.350 1.448 9
%      bostname(1,:)='BOSTOCK/NEW/047_2004.196';
%      bostname(2,:)='BOSTOCK/NEW/047_2004.197';
%      bostname(3,:)='BOSTOCK/NEW/047_2004.198';
%      bostname(4,:)='BOSTOCK/NEW/047_2004.199';
%      bostname(5,:)='BOSTOCK/NEW/047_2004.200';
%      bostname(6,:)='BOSTOCK/NEW/047_2005.254';
%      bostname(7,:)='BOSTOCK/NEW/047_2005.255';
%      bostname(8,:)='BOSTOCK/NEW/047_2005.256';
%      bostname(9,:)='BOSTOCK/NEW/047_2005.260';
%047 040714 8 1375.075 1.603 9
%      bostname(1,:)='BOSTOCK/NEW/068_2004.198';
%      bostname(2,:)='BOSTOCK/NEW/068_2004.199';
%      bostname(3,:)='BOSTOCK/NEW/068_2004.200';
%      bostname(4,:)='BOSTOCK/NEW/068_2004.201';
%      bostname(5,:)='BOSTOCK/NEW/068_2005.256';
%      bostname(6,:)='BOSTOCK/NEW/068_2005.257';
%      bostname(7,:)='BOSTOCK/NEW/068_2005.258';
%      bostname(8,:)='BOSTOCK/NEW/068_2005.259';
%      bostname(9,:)='BOSTOCK/NEW/068_2005.260';
%      bostname(10,:)='BOSTOCK/NEW/068_2005.261';
     
PERMSTA=['PGC'
         'LZB'];
POLSTA=['SSIB'
        'SILB'
        %'KLNB'
        'TWKB'
        'MGCB'];
PERMROTS=[0 90 32 00;  %PGC , Yajun's "0" changed to 90.
          0 90 54 00]; %LZB
% 002
 POLROTS=[6 85 33 86;  %SSIB from Yajun
          0 90 39 20;  %SILB
          %0 90  7 00;  %KLNB
          4 75 38 00;  %MGCB
          4 70 48 00]; %TWKB
 %POLROTS=[6 85 33 86;  %SSIB from Yajun
 %POLROTS=[0 90 -25 86;  %SSIB from John
% 047
%  POLROTS=[4 80 29 86;  %SSIB
%           0 90 41 20;  %SILB
%           %0 90 00 00;  %KLNB
%           2 85 49 00;  %MGCB
%           4 65 24 00]; %TWKB
% 068
%  POLROTS=[2 50 25 00; %KLNB
%           4 60 14 00; %TWKB
%           5 65 41 00]; %MGCB
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;
POLROTS(:,1)=round(POLROTS(:,1)*(100/40)); %Polaris stations at 100 sps; table has split time at 40 sps.
nsta=size(POLSTA,1);
sps=100;
%hi=15.;
%lo=0.14;
%lo=1;
hi=6.;
lo=1.5;
hi=7;
lo=1.5;
% hi=16.;
% lo=0.75;
npo=2;
npa=2;
winlen=60*100;
stack=zeros(nsta,winlen);
stackort=zeros(nsta,winlen);
STAopt=zeros(nsta,sps*24*3600);
STAort=zeros(nsta,sps*24*3600);
scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite;
xero=[0 0;
      6000 0];
h=figure('Position',[0.1*wid 1 2.5*wid hite]); %center

for nd=1:length(timoffrot(:,1))

    bostocks=load(bostname(nd,:));

    bostsec=3600*(bostocks(:,3)-1)+bostocks(:,4)-22.675; %22.675 b.c. 002_ catalogs already "corrected" for PGC start; from tempseps.f
    bostsamp=round(bostsec*100);

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
    IDENTIF=[YEAR,'.',JDAY,'.',int2str(timoffrot(nd,6)),'.',int2str(timoffrot(nd,7)), ...
        '.',int2str(timoffrot(nd,8)),'.',int2str(timoffrot(nd,9)),'.',int2str(timoffrot(nd,10))]
    direc=[YEAR,'/',MO,'/'];
    prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];

    %Read data:
    for ista=1:nsta
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
        if year==2003 && jday<213
            [STAEf]=20.0e-3*bandpass(STAE,sps,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
            [STANf]=20.0e-3*bandpass(STAN,sps,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
        else
            [STAEf]=4.0e-3*bandpass(STAE,sps,lo,hi,npo,npa,'butter'); 
            [STANf]=4.0e-3*bandpass(STAN,sps,lo,hi,npo,npa,'butter'); 
        end
        STA=STAEf+1i*STANf;
        STAfastslow=STA*exp(-1i*POLROTS(ista,2));
        STAslow=real(STAfastslow);
        STAfast=imag(STAfastslow);
        len=length(STA);
        STAslow(20:len-20)=STAslow(20+POLROTS(ista,1):len-20+POLROTS(ista,1));
        STAsplitcorrected=(STAslow+1i*STAfast)*exp(1i*POLROTS(ista,2));
        STAsplitcorrectedrot=STAsplitcorrected*exp(-1i*POLROTS(ista,3));
        STAscrot=STAsplitcorrectedrot;

        STAsoff=POLROTS(ista,4);
        STAtoff=round(STAsoff/40); %really 40.  That's how it's listed.

        if STAtoff > -1
            STAscrot(1:tracelen-STAsoff)=STAscrot(STAsoff+1:tracelen);
            STAscrot(tracelen-STAsoff+1:tracelen)=0;
        else
            STAscrot(-STAsoff+1:tracelen)=STAscrot(1:tracelen+STAsoff);
            STAscrot(1:-STAsoff)=0;
        end
        STAopt(ista,:)=real(STAscrot);
        STAort(ista,:)=imag(STAscrot);
    end

    for n=1:length(bostsec)
        if bostsamp(n)+winlen-1 <= len
            stack=stack+STAopt(:,bostsamp(n):bostsamp(n)+winlen-1);
            stackort=stackort+STAort(:,bostsamp(n):bostsamp(n)+winlen-1);
        end
        for ista=1:nsta
            subplot(3,2,ista,'align')
            plot(stack(ista,:),'b')
            %xlim([825 1045]) %xlim([895 975]) SILB
            %xlim([891 1111]) %xlim([961 1041]) SSIB
            %xlim([804 1024]) %xlim([874 954]) PGC
            xmin=1500; xmax=3500;
            xmin=1800; xmax=3000;
            %xlim([xmin xmax]) 
            ylim([min(stack(ista,xmin:xmax)) max(stack(ista,xmin:xmax))])
            title([bostname(1,13:15),'  ',POLSTA(ista,:),'  ',int2str(timoffrot(1,1)),'-',int2str(timoffrot(size(timoffrot,1),1)), ...
                '  ',num2str(lo),'-',int2str(hi),' Hz  ',int2str(npa),' pass'])
        end
        drawnow
    end
end
disp=zeros(nsta,xmax-xmin+1);
ymin=0*(1:nsta);
ymax=0*(1:nsta);
for ista=1:nsta
    disp(ista,:)=cumsum(stack(ista,xmin:xmax));
    ymin(ista)=min(stack(ista,xmin:xmax)); 
    ymax(ista)=max(stack(ista,xmin:xmax));
    mindisp=min(disp(ista,:));
    maxdisp=max(disp(ista,:));
    disp(ista,:)=disp(ista,:)*ymin(ista)/mindisp;
    subplot(3,2,ista,'align')
    hold on
%     if lo < 0.5
%         plot(xmin:xmax,disp(ista,:),'r')
%     end
    plot(xero(:,1),xero(:,2),'g--')
    box on
end
[stackmins minlocs]=min(stack(:,xmin:xmax),[],2);
minlocs=minlocs+xmin-1;
meanloc=round(mean(minlocs));
tostack=zeros(nsta,xmax-xmin+1+400);
subplot(3,2,ista+1,'align')
hold on
for ista=1:nsta
    tostack(ista,:)=stack(ista,xmin+(minlocs(ista)-meanloc)-200:xmax+(minlocs(ista)-meanloc)+200);
    if ista==1
        plot(0:0.01:(size(tostack,2)-1)/100,tostack(ista,:)/(-ymin(ista)),'b')
    else
        plot(0:0.01:(size(tostack,2)-1)/100,tostack(ista,:)/(-ymin(ista)),'r')
    end
end
stacksum=sum(tostack);
% plot(0:0.01:(size(tostack,2)-1)/100,-stacksum/min(stacksum),'r')
% % dispsum=cumsum(stacksum);
% % plot(0:0.01:(size(tostack,2)-1)/100,-dispsum/min(dispsum),'r')
% plot(xero(:,1),xero(:,2),'g--')
axis([0 (size(tostack,2)-1)/100 -1 1])
xlabel('seconds')
title('Stack')
box on
set(h,'PaperPosition',[0.25 0.25 8 10.5])
print(h,'-depsc',['temps_',bostname(1,13:15),'_',int2str(timoffrot(1,1)),'-',int2str(timoffrot(size(timoffrot,1),1)),'_', ...
    num2str(lo),'-',num2str(hi),'Hz_',int2str(npa),'pass.eps'])

g=figure('Position',[0.1*wid 1 2.5*wid hite]); %center
for ista=1:nsta
    subplot(3,2,ista,'align')
    hold on
    plot(stackort(ista,:),'b')
    plot(xero(:,1),xero(:,2),'g--')
    %xmin=1500; xmax=3500;
    xlim([xmin xmax]) 
    ylim([min(stack(ista,xmin:xmax)) max(stack(ista,xmin:xmax))])
    title([bostname(1,13:15),'   ',POLSTA(ista,:),' orth. ',int2str(timoffrot(1,1)),'-',int2str(timoffrot(size(timoffrot,1),1)), ...
        '  ',num2str(lo),'-',int2str(hi),' Hz'])
    box on
end

