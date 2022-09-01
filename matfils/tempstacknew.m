%Reads in seismograms and Bostock detections. Filters and rotates seisms, etc.  Stacks based on times of B. detections.
%At SILB, SSIB, PGC in user-specified passband.
format short e
clear all
%close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

fam='065'
if isequal(fam,'002')
%GET RID OF -22.675 for families other than 002!!!
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
% 002 %PG SS SI: 80 115  50
PERMROTS=[0 90 32 0;  %PGC , Yajun's "0" changed to 90. 4th column here is shift w.r.t. PGC
          0 90 54 9;  %LZB
          0  0  0 9]; %LZBz
POLROTS=[6 85 33  0; %86;  %SSIB from Yajun
%POLROTS=[0 90 -25  0; %86;  %SSIB from JArmb.
         0 90 39  0; %20;  %SILB
         0 90  7 -4;  %KLNB
         4 70 48 -26; %MGCB
         4 75 38 -5;  %TWKB
         0  0  0 20;  %SILBz
         0  0  0 86]; %SSIBz
stas=['PGC  '
      'SSIB '
      'SILB '];
correction=-22.675;
elseif isequal(fam,'068')
% 068 %GET RID OF -22.675 for families other than 002!!!
timoffrot=[2004 198;
           2004 199;
           2004 200;
           2004 201;
           2005 256;
           2005 257;
           2005 258;
           2005 259;
           2005 260;
           2005 261];
bostname=['BOSTOCK/NEW/068_2004.198';
          'BOSTOCK/NEW/068_2004.199';
          'BOSTOCK/NEW/068_2004.200';
          'BOSTOCK/NEW/068_2004.201';
          'BOSTOCK/NEW/068_2005.256';
          'BOSTOCK/NEW/068_2005.257';
          'BOSTOCK/NEW/068_2005.258';
          'BOSTOCK/NEW/068_2005.259';
          'BOSTOCK/NEW/068_2005.260';
          'BOSTOCK/NEW/068_2005.261'];
PERMROTS=[0  0  0  0;  %PGC  4th column here is shift w.r.t. TWKB based on 2004/2005 catalogs.  Better set to zero?
          8 65  6  0;  %LZB (offset is -6)
          0  0  0  0]; %LZBz
POLROTS=[0  0  0  0;  %SSIB
         0  0  0  0;  %SILB
         2 50 25  0;  %KLNB (delay w.r.t. TWKB is fictitious.
         5 65 41  0;  %MGCB (offset is 1)
         4 60 14  0;  %TWKB
         0  0  0  0;  %SILBz
         0  0  0  0]; %SSIBz
stas=['TWKB '
      'LZB  '
      'MGCB '];
correction=0;
elseif isequal(fam,'065')
% timoffrot=[2003 066;
%            2003 067;
%            2003 068;
timoffrot=[2004 200;
           2004 201;
           2004 204;
           2005 259;
           2005 260];
% bostname=['BOSTOCK/NEW/065_2003.066';
%           'BOSTOCK/NEW/065_2003.067';
%           'BOSTOCK/NEW/065_2003.068';
bostname=['BOSTOCK/NEW/065_2004.200';
          'BOSTOCK/NEW/065_2004.201';
          'BOSTOCK/NEW/065_2004.204';
          'BOSTOCK/NEW/065_2005.259';
          'BOSTOCK/NEW/065_2005.260'];
PERMROTS=[0  0  0  0;  %PGC  4th column here is shift w.r.t. TWKB based on 2004/2005 catalogs.  Better set to zero?
          0  0  0  0;  %LZB (offset is -6)
          0  0  0  0]; %LZBz
POLROTS=[0   0   0  0;  %SSIB
%        1  20 155  0;  %SILB
         1 -70  65  0;  %SILB
%        3 160 120  0;  %KLNB (delay w.r.t. TWKB is fictitious.
         3  70  30  0;  %KLNB (delay w.r.t. TWKB is fictitious.
         0   0   0  0;  %MGCB (offset is 1)
%        9 155 110  0;  %TWKB
         9  65  20  0;  %TWKB
         0   0   0  0;  %SILBz
         0   0   0  0]; %SSIBz
stas=['TWKB '
      'KLNB '
      'SILB '];
correction=0;
end
nsta=size(stas,1);
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
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;
sps=100; %40; 
PERMROTS(:,1)=round(PERMROTS(:,1)*(sps/40)); %Polaris stations at 100 sps; table has split time at 40 sps.
POLROTS(:,1)=round(POLROTS(:,1)*(sps/40)); %Polaris stations at 100 sps; table has split time at 40 sps.
hi=6.5;
lo=1.25;
% hi=15;
% lo=1;
npo=2;
npa=2;

winlen=60*sps;
stack=zeros(nsta,winlen);
stackort=zeros(nsta,winlen);
STAopt=zeros(nsta,sps*24*3600);
STAort=zeros(nsta,sps*24*3600);
scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite;
xero=[0 0;
      6000*sps/100 0];
h=figure('Position',[0.1*wid 1 2.5*wid hite]); %center
%I think the following used only in Glitches; these are dummy parameters.
nwin=1; winoff=1; igstart=1;

tot=0;
for nd=1:length(timoffrot(:,1))
    timoffrot(nd,:)
    bostocks=load(bostname(nd,:));
    tot=tot+size(bostocks,1)
    bostsec=3600*(bostocks(:,3)-1)+bostocks(:,4)-correction; %-22.675 b.c. 002_ catalogs already "corrected" for PGC start; from tempseps.f
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
            [opt,ort,nzeros,timsSTAperm]=readperms(prename,PERMSTA,PERMROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
        end
        [LIA,idx]=ismember(stas(ista,:),POLSTA,'rows');
        if LIA
            found=found+LIA; %better be 1
            if year==2003 && jday<213
                fact=20.0e-3;
            else
                fact=4.0e-3; 
            end
            [opt,ort,nzeros,timsSTApol]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
        end
        found=found
        STAopt(ista,:)=opt;
        STAort(ista,:)=ort;
    end
    len=size(STAopt,2);

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
            xmin=round(1800*sps/100); xmax=round(3000*sps/100);
            %xlim([xmin xmax]) 
            ylim([min(stack(ista,xmin:xmax)) max(stack(ista,xmin:xmax))])
            title([bostname(1,13:15),'  ',stas(ista,:),'  ',int2str(timoffrot(1,1)),'-',int2str(timoffrot(size(timoffrot,1),1)), ...
                '  ',num2str(lo),'-',num2str(hi),' Hz  ',int2str(npa),' pass'])
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
for ista=1:nsta
    subplot(3,2,ista,'align')
    hold on
    plot(stackort(ista,:),'g--')
    box on
end
[stackmins minlocs]=min(stack(:,xmin:xmax),[],2);
minlocs=minlocs+xmin-1;
[stackmaxs maxlocs]=max(stack(:,xmin:xmax),[],2);
maxlocs=maxlocs+xmin-1;
zerolocs=round(0.5*(minlocs+maxlocs));
%meanloc=round(mean(minlocs));
subplot(3,2,ista+1,'align')
plotszm=250;
plotszp=950;
%tostack=zeros(nsta,xmax-xmin+1+2*plotsz);
tostack=zeros(nsta,1+plotszm+plotszp);
hold on
for ista=1:nsta
    %tostack(ista,:)=stack(ista,xmin+(minlocs(ista)-meanloc)-plotsz:xmax+(minlocs(ista)-meanloc)+plotsz);
    tostack(ista,:)=stack(ista,zerolocs(ista)-plotszm:zerolocs(ista)+plotszp);
    ym=max(ymax(ista),-ymin(ista));
    if ista==1
        plot(tostack(ista,:)/ym,'b')
    elseif ista==2
        plot(tostack(ista,:)/ym,'r')
    else
        plot(tostack(ista,:)/ym,'k')
    end
end
%SSSI=xcorr(tostack(2,:),tostack(3,:),20,'coeff');
stacksum=sum(tostack);
% plot(0:0.01:(size(tostack,2)-1)/100,-stacksum/min(stacksum),'r')
% % dispsum=cumsum(stacksum);
% % plot(0:0.01:(size(tostack,2)-1)/100,-dispsum/min(dispsum),'r')
% plot(xero(:,1),xero(:,2),'g--')
%axis([0 (size(tostack,2)-1)/sps -1 1])
axis([0 size(tostack,2)-1 -1 1])
xlabel('samples')
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
    title([bostname(1,13:15),'   ',stas(ista,:),' orth. ',int2str(timoffrot(1,1)),'-',int2str(timoffrot(size(timoffrot,1),1)), ...
        '  ',num2str(lo),'-',num2str(hi),' Hz'])
    box on
end

for ista=1:nsta
    savdat(:,1)=stack(ista,:);
    if stas(ista,4) == ' '
        fid = fopen([stas(ista,1:3),'opt_',bostname(1,13:15),'_',num2str(lo),'-',num2str(hi),'Hz_',int2str(npa),'pass_',num2str(sps),'sps'],'w');
    else
        fid = fopen([stas(ista,1:4),'opt_',bostname(1,13:15),'_',num2str(lo),'-',num2str(hi),'Hz_',int2str(npa),'pass_',num2str(sps),'sps'],'w');
    end
    fprintf(fid,'%9.5f\n',savdat');
    fclose(fid);
end
