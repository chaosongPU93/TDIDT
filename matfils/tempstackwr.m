%Reads in seismograms and Bostock detections. Filters and rotates seisms, etc.  Stacks based on times of B. detections.
%At SILB, SSIB, PGC in user-specified passband.
%This one writes out E, S, and N for later manipulation.
format short e
clear all
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

%GET RID OF -22.675 for families other than 002!!!
% 002
timoffrot=[2003 062;
           2003 063;
%          2003 064; %At least one bad window at SSIB
%          2003 068;
           2004 196;
           2004 197;
           2004 198;
           2004 199;
%          2004 200;
           2005 254;
           2005 255;
           2005 256;
           2005 257];
%          2005 260];
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
     bostname=['BOSTOCK/NEW/002_2003.062';
               'BOSTOCK/NEW/002_2003.063';
              %'BOSTOCK/NEW/002_2003.064';
               'BOSTOCK/NEW/002_2004.196';
               'BOSTOCK/NEW/002_2004.197';
               'BOSTOCK/NEW/002_2004.198';
               'BOSTOCK/NEW/002_2004.199';
               'BOSTOCK/NEW/002_2005.254';
               'BOSTOCK/NEW/002_2005.255';
               'BOSTOCK/NEW/002_2005.256';
               'BOSTOCK/NEW/002_2005.257'];
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
     
PERMSTA=['PGC  '
         'LZB  '
         'LZBz '];
POLSTA =['SSIB '
         'SILB '
         'KLNB '
         'MGCB '
         'TWKB '
         'SILBz'
         'SSIBz'];

sps=100;
stas=['PGC  '
      'SSIB '
      'SILB '];
nsta=size(stas,1);
hi=6.5;
lo=1.25;
% %
% hi=16.;
% lo=0.5;
% %
% hi=12.;
% lo=0.75;
% %
% hi=12.;
% lo=0.25;
npo=2;
npa=2;
%npa=1;

winlen=60*sps;
stackE=zeros(nsta,winlen);
stackN=zeros(nsta,winlen);
stackZ=zeros(nsta,winlen);
staE=zeros(nsta,sps*24*3600);
staN=zeros(nsta,sps*24*3600);
staZ=zeros(nsta,sps*24*3600);
scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite;
xero=[0 0;
      6000*sps/100 0];
h=figure('Position',[0.1*wid 1 2.5*wid hite]); %center
%I think the following used only in Glitches; these are dummy parameters.
nwin=1; winoff=1; igstart=1;

for nd=1:length(timoffrot(:,1))
    timoffrot(nd,:)
    bostocks=load(bostname(nd,:));

    bostsec=3600*(bostocks(:,3)-1)+bostocks(:,4)-22.675; %22.675 b.c. 002_ catalogs already "corrected" for PGC start; from tempseps.f
    bostsamp=round(bostsec*sps);

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
    direc=[YEAR,'/',MO,'/'];
    prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];

    for ista=1:nsta
        found=0;
        [LIA,idx]=ismember(stas(ista,:),PERMSTA,'rows');
        if LIA
            found=found+LIA;
            if strcmp(PERMSTA(idx,1:3),'PGC')
                fact=1.6e-4;
            elseif strcmp(PERMSTA(idx,1:3),'LZB')
                fact=4.e-3;
            end
            [datE,datN,datZ,timsSTAperm]=readpermsnorots(prename,PERMSTA,idx,sps,lo,hi,npo,npa,fact);
        end
        [LIA,idx]=ismember(stas(ista,:),POLSTA,'rows');
        if LIA
            found=found+LIA; %better be 1
            if year==2003 && jday<213
                fact=20.0e-3;
            else
                fact=4.0e-3; 
            end
            [datE,datN,datZ,timsSTApol]=readpolsnorots(prename,POLSTA,idx,sps,lo,hi,npo,npa,fact);
        end
        found=found
        staE(ista,:)=datE;
        staN(ista,:)=datN;
        staZ(ista,:)=datZ;
    end
    len=size(staE,2);

    for n=1:length(bostsec)
        if bostsamp(n)+winlen-1 <= len
            stackE=stackE+staE(:,bostsamp(n):bostsamp(n)+winlen-1);
            stackN=stackN+staN(:,bostsamp(n):bostsamp(n)+winlen-1);
            stackZ=stackZ+staZ(:,bostsamp(n):bostsamp(n)+winlen-1);
        end
        for ista=1:nsta
            subplot(3,nsta,ista,'align')
            plot(stackE(ista,:),'b')
            xmin=round(1800*sps/100); xmax=round(3000*sps/100);
            ylim([min(stackE(ista,xmin:xmax)) max(stackE(ista,xmin:xmax))])
            title([bostname(1,13:15),'  ',stas(ista,:),'  ',int2str(timoffrot(1,1)),'-',int2str(timoffrot(size(timoffrot,1),1)), ...
                '  ',num2str(lo),'-',num2str(hi),' Hz  ',int2str(npa),' pass  E'])
            subplot(3,nsta,nsta+ista,'align')
            plot(stackN(ista,:),'r')
            xmin=round(1800*sps/100); xmax=round(3000*sps/100);
            ylim([min(stackN(ista,xmin:xmax)) max(stackN(ista,xmin:xmax))])
            title('N')
            subplot(3,nsta,2*nsta+ista,'align')
            plot(stackZ(ista,:),'k')
            xmin=round(1800*sps/100); xmax=round(3000*sps/100);
            ylim([min(stackZ(ista,xmin:xmax)) max(stackZ(ista,xmin:xmax))])
            title('Z')
        end
        drawnow
    end
end

for ista=1:nsta
    savdat(:,1)=stackE(ista,:);
    savdat(:,2)=stackN(ista,:);
    savdat(:,3)=stackZ(ista,:);
    if stas(ista,4) == ' '
        fid = fopen([stas(ista,1:3),'_',bostname(1,13:15),'_',num2str(lo),'-',num2str(hi),'Hz_',int2str(npa),'pass'],'w');
    else
        fid = fopen([stas(ista,1:4),'_',bostname(1,13:15),'_',num2str(lo),'-',num2str(hi),'Hz_',int2str(npa),'pass'],'w');
    end
    fprintf(fid,'%9.5f %9.5f %9.5f\n',savdat');
    fclose(fid);
end
         
