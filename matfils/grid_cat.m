%locate the detections.  Scatter-plot different attributes (max amplitude;
%max cc; etc.)

%%%%
% Hi, Chao.  You'll only need the first few sections of grid_cat.m to see
% how p085p020 is used (don't worry about the name; I used to use offsets
% of 85,20 but now use 86,20; you can see how that's used in centPGSS and
% cent PGSI).
% Allan
%%%%

clear
close all
rads=pi/180.;
erad=6372.028;
lat0=48.0+26.32/60.;    % lat of the event inverted from (0,0), off12,off13
lon0=123.0+35.07/60.;   % lon of the event inverted from (0,0), off12,off13
srad=erad*cos(lat0*rads)

fam='002';
A=load('p085p020');
mm=25; nn=25;
dlat=reshape(A(:,3),mm,nn); %degrees
mlat=reshape(A(:,4),mm,nn); %minutes
lat=dlat+mlat/60.;
%%% Vq = interp2(V,K) returns the interpolated values on a refined grid 
%%% formed by repeatedly halving the intervals K times in each dimension.
%%% This results in 2^K-1 interpolated points between sample values.
lat2=interp2(lat,3,'spline');    % 24*(2^3-1)+25=193
dlon=reshape(A(:,5),mm,nn); %degrees
mlon=reshape(A(:,6),mm,nn); %minutes
lon=dlon+mlon/60.;
lon2=interp2(lon,3,'spline');
dep=reshape(A(:,7),mm,nn);

%%%%%%%%% algorithm %%%%%%%%%%%%%%%%
% arr2 = arr1 - off12 + 86(centPGSS)
% arr3 = arr1 - off13 + 20(centPGSI)
%  ||
% arr2 - arr1 = 86(centPGSS) - off12 == 1st col in A
% arr3 - arr1 = 20(centPGSI) - off13 == 2nd col in A
%  ||
% off12 = 86(centPGSS) - (arr2 - arr1)
% off13 = 20(centPGSI) - (arr3 - arr1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
centPGSS=86
centPGSI=20
offPGSS=reshape(centPGSS-A(:,1),mm,nn); %PGSS is uniform in columns; changes across rows.
offPGSI=reshape(centPGSI-A(:,2),mm,nn); %PGSI is uniform in rows; changes down columns.
offPGSS2=interp2(offPGSS,3,'spline'); %scalar is number of times grid is halved
offPGSI2=interp2(offPGSI,3,'spline'); %scalar is number of times grid is halved
offPGSS1=offPGSS2(1,:); %a 1-d array of offsets at 1/4 sample
offPGSI1=offPGSI2(:,1);
mPGSS=length(offPGSS1); %
nPGSI=length(offPGSI1);
% nPGSI=length(offPGSS2); %This and the previous line look iffy!  Should be offPGSS2 and offPGSI2?  Probably works only b.c. grid is square.
% In p085p020, order is PGSS PISI.  PGSI changes more rapidly.
% in output of wiggledaywig2, order is [timswin(n) xmaxPGSIntmp(n) xmaxPGSSntmp(n) ...

%%%%%% 
% fname= 'mapall_002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
% %fname= 'mapall_002.loff1.5.ccmin0.44.nponpa22_1.5-6-ms26-4s';
% for n=1:size(fname,1)
% n=1;
% len=0;
% for i=1:size(fname,2)
%     len=len+(1-strcmp(fname(n,i),' ')); %strcmp returns 1 if equal, 0 if not
% end
% locfile=load(fname(n,1:len));
%%%%%

%%%%%
workpath = getenv('ALLAN');
datapath = strcat(workpath, '/PGCtrio/MAPS');
cd(datapath);
system('rm -f mapall.002.loff2.1.ccmin0.44.nponpa22.ms19_1.25-6.5_4s40sps');
system('cat map*.002.loff2.1.ccmin0.44.nponpa22.ms19_1.25-6.5_4s40sps > mapall.002.loff2.1.ccmin0.44.nponpa22.ms19_1.25-6.5_4s40sps');
locfile = load('mapall.002.loff2.1.ccmin0.44.nponpa22.ms19_1.25-6.5_4s40sps');
%%%%%


scrsz=get(0,'ScreenSize');
wid=scrsz(3);
hite=scrsz(4);
scrat=wid/hite;
%h=figure('Position',[0.3*wid hite/2 wid/2 hite/2]); %center
%h=figure('Position',[1 1 wid hite]); %center


ncol=size(locfile,2);
locfile(:,16:19)=0;
PGSIs=locfile(:,2);
PGSSs=locfile(:,3);
nevs=length(PGSIs);
%fid = fopen([fname(n,1:len),'.dx'],'w');
nhits=zeros(nPGSI,mPGSS);
maxcc=zeros(nPGSI,mPGSS);
allmaxcc=zeros(nPGSI,mPGSS,40);
frac=zeros(nPGSI,mPGSS);
allfrac=zeros(nPGSI,mPGSS,40);
Ampsq=zeros(nPGSI,mPGSS);
allAmpsq=zeros(nPGSI,mPGSS,40);
UArrivoPrior=zeros(nPGSI,mPGSS);
allUArrivoPrior=zeros(nPGSI,mPGSS,40);
for i=1:nevs
    if PGSIs(i)>=min(offPGSI1) && PGSIs(i)<=max(offPGSI1) && PGSSs(i)>=min(offPGSS1) && PGSSs(i)<=max(offPGSS1) %if w/in Armbruster's grid
        for j=1:nPGSI
            if (abs(PGSIs(i)-offPGSI1(j))) <= 1.e-6
                PGSImatch=j;                            %This just assigns PGSI offset
            end
        end
        for j=1:mPGSS
            if (abs(PGSSs(i)-offPGSS1(j))) <= 1.e-6
                PGSSmatch=j;                            %This just assigns PGSS offset
            end
        end
        nhits(PGSImatch,PGSSmatch)=nhits(PGSImatch,PGSSmatch)+1;
        maxcc(PGSImatch,PGSSmatch)=max(maxcc(PGSImatch,PGSSmatch),locfile(i,4)); %max CC
        allmaxcc(PGSImatch,PGSSmatch,nhits(PGSImatch,PGSSmatch))=locfile(i,4); %This stores them all.  Makes maxcc kinda obsolete.
        frac(PGSImatch,PGSSmatch)=max(frac(PGSImatch,PGSSmatch),locfile(i,9));   %dot prod (arrival/window)
        allfrac(PGSImatch,PGSSmatch,nhits(PGSImatch,PGSSmatch))=locfile(i,9); %This stores them all.  Makes frac kinda obsolete.
        Ampsq(PGSImatch,PGSSmatch)=max(Ampsq(PGSImatch,PGSSmatch),locfile(i,6)); %Amp^2 arrival
        allAmpsq(PGSImatch,PGSSmatch,nhits(PGSImatch,PGSSmatch))=locfile(i,6); %This stores them all.  Makes Ampsq kinda obsolete.
        UArrivoPrior(PGSImatch,PGSSmatch)=max(UArrivoPrior(PGSImatch,PGSSmatch),2.5*locfile(i,6)/locfile(i,12)); %Amp^2 arrival
        allUArrivoPrior(PGSImatch,PGSSmatch,nhits(PGSImatch,PGSSmatch))=2.5*locfile(i,6)/locfile(i,12); %This stores all; makes prev obsolete.
        %%%% IMPORTANT %%%%
        dy=rads*(lat2(PGSImatch,PGSSmatch)-lat0)*erad;
        dx=-rads*(lon2(PGSImatch,PGSSmatch)-lon0)*srad; %(minus bcause 123 is -123)
        %%%%%%%%%%%%%%%%%%%
        locfile(i,14)=dy;                               %This assigns north coord
        locfile(i,15)=dx;                               %This assigns east coord
        angrot=-45*pi/180;
        spot=complex(dx,dy);
        spotrot=spot*exp(1i*angrot);    % *exp(1i*angrot) means rotate angrot counter-clockwise
        locfile(i,16)=real(spotrot);                    %This assigns pseudo-east coord
        locfile(i,17)=imag(spotrot);                    %This assigns pseudo-north coord
        locfile(i,18)=PGSImatch;
        locfile(i,19)=PGSSmatch;
    else
        locfile(i,14)=9.9e6;                               %This assigns north coord
        locfile(i,15)=9.9e6;                               %This assigns east coord
        locfile(i,16)=9.9e6;                    %This assigns pseudo-east coord
        locfile(i,17)=9.9e6;                    %This assigns pseudo-north coord
        locfile(i,18)=1;
        locfile(i,19)=1;
    end
end
medmaxcc=zeros(nPGSI,mPGSS);
medfrac=zeros(nPGSI,mPGSS);
medAmpsq=zeros(nPGSI,mPGSS);
medUArrivoPrior=zeros(nPGSI,mPGSS);
for j=1:nPGSI
    for k=1:mPGSS
        medmaxcc(j,k)=median(allmaxcc(j,k,1:nhits(j,k)));
        medfrac(j,k)=median(allfrac(j,k,1:nhits(j,k)));
        medAmpsq(j,k)=median(allAmpsq(j,k,1:nhits(j,k)));
        medUArrivoPrior(j,k)=median(allUArrivoPrior(j,k,1:nhits(j,k)));
    end
end
h=figure('Position',[wid/3 1 wid/2.8 hite]);  %Splits.  These variable names make no logical sense; just inherited from a code that works.
SEbord=-1.35; NWbord=1.6; %SEbord=1.3; NWbord=1.5;
SE=[-10 SEbord;
    10 SEbord];
NW=[-10 NWbord;
    10 NWbord];
angrot=45*pi/180;
SEline=complex(SE(:,1),SE(:,2));
NWline=complex(NW(:,1),NW(:,2));
SEline=SEline*exp(1i*angrot);
NWline=NWline*exp(1i*angrot);
subplot(2,1,1,'align')
scatter(locfile(:,15),locfile(:,14),14,locfile(:,2),'filled')
hold on
caxis([-20 20])
colormap(jet)
colorbar
axis equal
%axis([-8 3 -4 4])
axis([-10 6 -7 7])
axis([-8 5 -6 6])
%axis([-12 8 -9 9])
title('PGSI')
plot(SEline,'r--','linewidth',2)
plot(NWline,'r--','linewidth',2)
box on
%
subplot(2,1,2,'align')
scatter(locfile(:,15),locfile(:,14),14,locfile(:,3),'filled')
hold on
caxis([-20 20])
colormap(jet)
colorbar
axis equal
%axis([-8 3 -4 4])
axis([-10 6 -7 7])
axis([-8 5 -6 6])
%axis([-12 8 -9 9])
title('PGSS')
plot(SEline,'r--','linewidth',2)
plot(NWline,'r--','linewidth',2)
box on
set(h,'PaperPosition',[0.25 0.25 8 10])
print(h,'-depsc',['grid',fam,'offsets.eps'])

h=figure('Position',[wid/3 1 wid/2 hite]);  %Splits.  These variable names make no logical sense; just inherited from a code that works.
angs=0:pi/50:2*pi;
circ=0.8*exp(1i*angs);
SEbord=-1.35; NWbord=1.6; %SEbord=1.3; NWbord=1.5;
SE=[-10 SEbord;
    10 SEbord];
NW=[-10 NWbord;
    10 NWbord];
angrot=45*pi/180;
SEline=complex(SE(:,1),SE(:,2));
NWline=complex(NW(:,1),NW(:,2));
SEline=SEline*exp(1i*angrot);
NWline=NWline*exp(1i*angrot);
scatter(locfile(:,15),locfile(:,14),14,locfile(:,2)-locfile(:,3),'filled')
hold on
northoffset=-0.5; %(in km)
eastoffset=-0.1; %(in km)
plot(real(circ)+eastoffset,imag(circ)+northoffset,'k','linewidth',2)
caxis([-22 22])
colormap(jet)
colorbar
axis equal
%axis([-7 3 -4 4])
axis([-10 6 -7 7])
axis([-8 5 -6 6])
%axis([-12 8 -9 9])
title('SSSI')
plot(SEline,'r--','linewidth',2)
plot(NWline,'r--','linewidth',2)
box on
set(h,'PaperPosition',[0.25 0.25 8 10])
print(h,'-depsc',['grid',fam,'SSSIoffsets.eps'])

for i=1:nevs
    locfile(i,20)=nhits(locfile(i,18),locfile(i,19));
    locfile(i,21)=maxcc(locfile(i,18),locfile(i,19));
    locfile(i,22)=frac(locfile(i,18),locfile(i,19));
    locfile(i,23)=Ampsq(locfile(i,18),locfile(i,19));
    locfile(i,24)=UArrivoPrior(locfile(i,18),locfile(i,19));
    locfile(i,25)=medmaxcc(locfile(i,18),locfile(i,19));
    locfile(i,26)=medfrac(locfile(i,18),locfile(i,19));
    locfile(i,27)=medAmpsq(locfile(i,18),locfile(i,19));
    locfile(i,28)=medUArrivoPrior(locfile(i,18),locfile(i,19));
end

h=figure('Position',[wid/3 1 wid/2 hite]);
colormap(jet)
dumlocfile=locfile;
A=find(locfile(:,20)>1); dumlocfile(A,:)=[];
scatter(dumlocfile(:,15),dumlocfile(:,14),10,log10(dumlocfile(:,20)))
hold on
dumlocfile=locfile;
A=find(locfile(:,20)==1); dumlocfile(A,:)=[];
scatter(dumlocfile(:,15),dumlocfile(:,14),10,log10(dumlocfile(:,20)),'filled')
SEbord=-1.35; NWbord=1.6; %SEbord=1.3; NWbord=1.5;
SE=[-10 SEbord;
    10 SEbord];
NW=[-10 NWbord;
    10 NWbord];
angrot=45*pi/180;
SEline=complex(SE(:,1),SE(:,2));
NWline=complex(NW(:,1),NW(:,2));
SEline=SEline*exp(1i*angrot);
NWline=NWline*exp(1i*angrot);
plot(SEline,'r--','linewidth',2)
plot(NWline,'r--','linewidth',2)
angs=0:pi/20:2*pi;
circ=0.8*exp(1i*angs);
%plot(real(circ),imag(circ)-0.5,'k')
%plot(real(circ)-0.1,imag(circ)-0.5,'k','linewidth',2)
plot([-6 -5],[3 3],'k','linewidth',4)
caxis([0 1.3])
xlabel('km East','fontsize',13);
ylabel('km North','fontsize',13)
axis equal
box on
colorbar
%axis([-6 3.2 -3.5 4.5])
%axis([-7 3 -4 4])
axis([-10 6 -7 7])
axis([-8 5 -6 6])
title('Log_{10} (# hits)')
set(h,'PaperPosition',[0.25 0.25 8 10])
print(h,'-depsc',['nhits',fam,'.eps'])

h=figure('Position',[wid/3 1 wid/2 hite]);
colormap(jet)
scatter(locfile(:,15),locfile(:,14),10,locfile(:,21),'filled') %maxcc
hold on
plot(SEline,'r--','linewidth',2)
plot(NWline,'r--','linewidth',2)
plot([-6 -5],[3 3],'k','linewidth',4)
%caxis([0 1.3])
xlabel('km East','fontsize',13);
ylabel('km North','fontsize',13)
axis equal
box on
colorbar
%axis([-6 3.2 -3.5 4.5])
%axis([-7 3 -4 4])
axis([-10 6 -7 7])
axis([-8 5 -6 6])
title('max(Max CC)')
set(h,'PaperPosition',[0.25 0.25 8 10])
print(h,'-depsc',['maxCC',fam,'.eps'])
%
h=figure('Position',[wid/3 1 wid/2 hite]);
scatter(locfile(:,15),locfile(:,14),10,locfile(:,25),'filled') %medmaxcc
hold on
plot(SEline,'r--','linewidth',2)
plot(NWline,'r--','linewidth',2)
plot([-6 -5],[3 3],'k','linewidth',4)
caxis([0.44 0.67])
xlabel('km East','fontsize',13);
ylabel('km North','fontsize',13)
axis equal
box on
colorbar
%axis([-6 3.2 -3.5 4.5])
%axis([-7 3 -4 4])
axis([-10 6 -7 7])
axis([-8 5 -6 6])
title('med(Max CC)')
set(h,'PaperPosition',[0.25 0.25 8 10])
print(h,'-depsc',['maxCCmed',fam,'.eps'])

h=figure('Position',[wid/3 1 wid/2 hite]);
colormap(jet)
scatter(locfile(:,15),locfile(:,14),10,locfile(:,22),'filled') %maxcc
hold on
plot(SEline,'r--','linewidth',2)
plot(NWline,'r--','linewidth',2)
plot([-6 -5],[3 3],'k','linewidth',4)
%caxis([0 1.3])
xlabel('km East','fontsize',13);
ylabel('km North','fontsize',13)
axis equal
box on
colorbar
%axis([-6 3.2 -3.5 4.5])
%axis([-7 3 -4 4])
axis([-10 6 -7 7])
axis([-8 5 -6 6])
title('Max[Dot Prod (arrival/window)]')
set(h,'PaperPosition',[0.25 0.25 8 10])
print(h,'-depsc',['dotprod',fam,'.eps'])
%
h=figure('Position',[wid/3 1 wid/2 hite]);
colormap(jet)
scatter(locfile(:,15),locfile(:,14),10,locfile(:,26),'filled') %maxcc
hold on
plot(SEline,'r--','linewidth',2)
plot(NWline,'r--','linewidth',2)
plot([-6 -5],[3 3],'k','linewidth',4)
%caxis([0 1.3])
xlabel('km East','fontsize',13);
ylabel('km North','fontsize',13)
axis equal
box on
colorbar
%axis([-6 3.2 -3.5 4.5])
%axis([-7 3 -4 4])
axis([-10 6 -7 7])
axis([-8 5 -6 6])
title('Med[Dot Prod (arrival/window)]')
set(h,'PaperPosition',[0.25 0.25 8 10])
print(h,'-depsc',['dotprodmed',fam,'.eps'])

h=figure('Position',[wid/3 1 wid/2 hite]);
colormap(jet)
scatter(locfile(:,15),locfile(:,14),10,log10(locfile(:,23)),'filled') %maxcc
hold on
plot(SEline,'r--','linewidth',2)
plot(NWline,'r--','linewidth',2)
plot([-6 -5],[3 3],'k','linewidth',4)
caxis([min(log10(locfile(:,23))) 1.5])
xlabel('km East','fontsize',13);
ylabel('km North','fontsize',13)
axis equal
box on
colorbar
%axis([-6 3.2 -3.5 4.5])
%axis([-7 3 -4 4])
axis([-10 6 -7 7])
axis([-8 5 -6 6])
title('Max Log_{10} (Amp^2 arrival)')
set(h,'PaperPosition',[0.25 0.25 8 10])
print(h,'-depsc',['Ampsq',fam,'.eps'])
%
h=figure('Position',[wid/3 1 wid/2 hite]);
colormap(jet)
scatter(locfile(:,15),locfile(:,14),10,log10(locfile(:,27)),'filled') %maxcc
hold on
plot(SEline,'r--','linewidth',2)
plot(NWline,'r--','linewidth',2)
plot([-6 -5],[3 3],'k','linewidth',4)
caxis([min(log10(locfile(:,23))) 1.5])
xlabel('km East','fontsize',13);
ylabel('km North','fontsize',13)
axis equal
box on
colorbar
%axis([-6 3.2 -3.5 4.5])
%axis([-7 3 -4 4])
axis([-10 6 -7 7])
axis([-8 5 -6 6])
title('Med Log_{10} (Amp^2 arrival)')
set(h,'PaperPosition',[0.25 0.25 8 10])
print(h,'-depsc',['Ampsqmed',fam,'.eps'])

h=figure('Position',[wid/3 1 wid/2 hite]);
colormap(jet)
scatter(locfile(:,15),locfile(:,14),10,log10(locfile(:,24)),'filled') %energy arrival/prior
hold on
plot(SEline,'r--','linewidth',2)
plot(NWline,'r--','linewidth',2)
plot([-6 -5],[3 3],'k','linewidth',4)
%caxis([min(log10(locfile(:,23))) 1.5])
xlabel('km East','fontsize',13);
ylabel('km North','fontsize',13)
axis equal
box on
colorbar
%axis([-6 3.2 -3.5 4.5])
%axis([-7 3 -4 4])
axis([-10 6 -7 7])
axis([-8 5 -6 6])
title('Max Log_{10} (Amp^2 arrival/prior)')
set(h,'PaperPosition',[0.25 0.25 8 10])
print(h,'-depsc',['ArrivoPrior',fam,'.eps'])
%
h=figure('Position',[wid/3 1 wid/2 hite]);
colormap(jet)
scatter(locfile(:,15),locfile(:,14),10,log10(locfile(:,28)),'filled') %energy arrival/prior
hold on
plot(SEline,'r--','linewidth',2)
plot(NWline,'r--','linewidth',2)
plot([-6 -5],[3 3],'k','linewidth',4)
%caxis([min(log10(locfile(:,23))) 1.5])
xlabel('km East','fontsize',13);
ylabel('km North','fontsize',13)
axis equal
box on
colorbar
%axis([-6 3.2 -3.5 4.5])
%axis([-7 3 -4 4])
axis([-10 6 -7 7])
axis([-8 5 -6 6])
title('Med Log_{10} (Amp^2 arrival/prior)')
set(h,'PaperPosition',[0.25 0.25 8 10])
print(h,'-depsc',['ArrivoPriormed',fam,'.eps'])
