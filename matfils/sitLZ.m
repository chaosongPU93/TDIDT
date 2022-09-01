%Sits at one spot (rotations; delays) and looks for the cross-correlation.
%close all
%Rotations & offsets:
% timoffrot=[253 14 52 52 +000 -021  40  55  55];
% timoffrot=[247 10 52 44 +069 +104 140 155 140];
timoffrot=[254 07 40 52 -011 -032 390 385 385];
jday=timoffrot(1);
JDAY=int2str(jday)
timlook=3600*timoffrot(2)+60*timoffrot(3)+timoffrot(4)
rotLZB=pi*(timoffrot(7)-90)/180;
rotTWKB=pi*(timoffrot(8)-90)/180;
rotMGCB=pi*(timoffrot(9)-90)/180;
TWKBsoff=timoffrot(5);
MGCBsoff=timoffrot(6);
TWKBtoff=TWKBsoff/40;
MGCBtoff=MGCBsoff/40;

%Read Armbruster's detections:
ArmCat=load('/data2/arubin/ARMB/LZTWMG/list.2005.lztwmgMAJ');
detects=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
m=0;
for n=1:length(detects)
     %if ArmCat(n,5)==TWKBsoff && ArmCat(n,6)==MGCBsoff
     if isequal(ArmCat(n,5:6),[TWKBsoff MGCBsoff]) ||...
        isequal(ArmCat(n,5:6),[TWKBsoff-1 MGCBsoff]) ||...
        isequal(ArmCat(n,5:6),[TWKBsoff+1 MGCBsoff]) ||...
        isequal(ArmCat(n,5:6),[TWKBsoff MGCBsoff-1]) ||...
        isequal(ArmCat(n,5:6),[TWKBsoff MGCBsoff+1]) ||...
        isequal(ArmCat(n,5:6),[TWKBsoff-1 MGCBsoff-1]) ||...
        isequal(ArmCat(n,5:6),[TWKBsoff-1 MGCBsoff+1]) ||...
        isequal(ArmCat(n,5:6),[TWKBsoff+1 MGCBsoff-1]) ||...
        isequal(ArmCat(n,5:6),[TWKBsoff+1 MGCBsoff+1]) 
         m=m+1;
         Detects(m)=detects(n);
     end
end
ArmCat=load('/data2/arubin/ARMB/LZTWMG/list.2005.lztwmg.64');
detects2_8=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
m=0;
for n=1:length(detects2_8)
     if isequal(ArmCat(n,5:6),[TWKBsoff MGCBsoff]) ||...
        isequal(ArmCat(n,5:6),[TWKBsoff-1 MGCBsoff]) ||...
        isequal(ArmCat(n,5:6),[TWKBsoff+1 MGCBsoff]) ||...
        isequal(ArmCat(n,5:6),[TWKBsoff MGCBsoff-1]) ||...
        isequal(ArmCat(n,5:6),[TWKBsoff MGCBsoff+1]) ||...
        isequal(ArmCat(n,5:6),[TWKBsoff-1 MGCBsoff-1]) ||...
        isequal(ArmCat(n,5:6),[TWKBsoff-1 MGCBsoff+1]) ||...
        isequal(ArmCat(n,5:6),[TWKBsoff+1 MGCBsoff-1]) ||...
        isequal(ArmCat(n,5:6),[TWKBsoff+1 MGCBsoff+1]) 
         m=m+1;
         Detects2_8(m)=detects2_8(n);
     end
end

%Read data:
prename=['2005.',JDAY,'.00.00.00.0000.CN'];
LZBEdat=[prename,'.LZB..BHE.D.SAC'];
LZBNdat=[prename,'.LZB..BHN.D.SAC'];
TWKBEdat=[prename,'.TWKB..HHE.D.SAC'];
TWKBNdat=[prename,'.TWKB..HHN.D.SAC'];
MGCBEdat=[prename,'.MGCB..HHE.D.SAC'];
MGCBNdat=[prename,'.MGCB..HHN.D.SAC'];

[LZBE,HdrDataLZB,tnuLZB,pobjLZB,timsLZB]=readsac(LZBEdat,0,'l');
[LZBN,~,~,~,~]=readsac(LZBNdat,0,'l');
[TWKBE,HdrDataTWKB,tnuTWKB,pobjTWKB,timsTWKB]=readsac(TWKBEdat,0,'l');
[TWKBN,~,~,~,~]=readsac(TWKBNdat,0,'l');
[MGCBE,HdrDataMGCB,tnuMGCB,pobjMGCB,timsMGCB]=readsac(MGCBEdat,0,'l');
[MGCBN,~,~,~,~]=readsac(MGCBNdat,0,'l');
tracelenLZB=length(LZBE);
tracelenTWKB=length(TWKBE);
tracelenMGCB=length(MGCBE);

%cosine taper before filtering:
x=(0:pi/80:pi/2-pi/80)';
LZBE(1:40)=sin(x).*LZBE(1:40);
LZBN(1:40)=sin(x).*LZBN(1:40);
x=flipud(x);
LZBE(tracelenLZB-39:tracelenLZB)=sin(x).*LZBE(tracelenLZB-39:tracelenLZB);
LZBN(tracelenLZB-39:tracelenLZB)=sin(x).*LZBN(tracelenLZB-39:tracelenLZB);
% For day 254, TWK:
TWKBE(571000:589000)=0;
TWKBE(1981000:1999000)=0;
TWKBN(571000:589000)=0;
TWKBN(1981000:1999000)=0;
x=(0:pi/2000:pi/2-pi/2000)';
TWKBE(570001:571000)=sin(x).*TWKBE(570001:571000);
TWKBE(1980001:1981000)=sin(x).*TWKBE(1980001:1981000);
TWKBN(570001:571000)=sin(x).*TWKBN(570001:571000);
TWKBN(1980001:1981000)=sin(x).*TWKBN(1980001:1981000);
x=flipud(x);
TWKBE(589001:590000)=sin(x).*TWKBE(589001:590000);
TWKBE(1999001:2000000)=sin(x).*TWKBE(1999001:2000000);
TWKBN(589001:590000)=sin(x).*TWKBN(589001:590000);
TWKBN(1999001:2000000)=sin(x).*TWKBN(1999001:2000000);
%
x=(0:pi/200:pi/2-pi/200)';
TWKBE(1:100)=sin(x).*TWKBE(1:100);
TWKBN(1:100)=sin(x).*TWKBN(1:100);
MGCBE(1:100)=sin(x).*MGCBE(1:100);
MGCBN(1:100)=sin(x).*MGCBN(1:100);
x=flipud(x);
TWKBE(tracelenTWKB-99:tracelenTWKB)=sin(x).*TWKBE(tracelenTWKB-99:tracelenTWKB);
TWKBN(tracelenTWKB-99:tracelenTWKB)=sin(x).*TWKBN(tracelenTWKB-99:tracelenTWKB);
MGCBE(tracelenMGCB-99:tracelenTWKB)=sin(x).*MGCBE(tracelenMGCB-99:tracelenMGCB);
MGCBN(tracelenMGCB-99:tracelenTWKB)=sin(x).*MGCBN(tracelenMGCB-99:tracelenMGCB);
% %Check up on tapering:
% figure 
% hold on
% plot(timsTWKB,TWKBE,'b');
% plot(timsTWKB,TWKBN,'k');
% xlim([5000 21000]);
% figure
% plot(LZBE);
% xlim([0 140]);
% figure
% plot(LZBN);
% xlim([tracelenLZB-79 tracelenLZB]);
% figure
% plot(TWKBE);
% xlim([tracelenTWKB-199 tracelenTWKB]);
% figure
% plot(TWKBN);
% xlim([0 200]);
% figure
% plot(MGCBE);
% xlim([tracelenMGCB-199 tracelenMGCB]);
% figure
% plot(MGCBN);
% xlim([0 200]);

% if(rem(tracelen,2)==0)
%     SeisData(tracelen)=[];
%     tims(tracelen)=[];
%     tracelen=tracelen-1;
% end
%Looks like filtering the 40-Hz data followed by interpolation introduces 
%less high-frequency stuff that does not appear in the original filtered
%40-Hz data:
%LZBE=spline(tims,SeisData,timsx);
%[LZBE]=bandpass(LZBE,100,1.0,5.0,2,1,'butter');

%Filter data:
hi=6.;
lo=1.5;
% hi=1.;
% lo=0.1;
npo=2;
npa=1;
[LZBEf]=4.0e-3*bandpass(LZBE,40,lo,hi,npo,npa,'butter');
[LZBNf]=4.0e-3*bandpass(LZBN,40,lo,hi,npo,npa,'butter');
[TWKBEf]=4.0e-3*bandpass(TWKBE,100,lo,hi,npo,npa,'butter');
[TWKBNf]=4.0e-3*bandpass(TWKBN,100,lo,hi,npo,npa,'butter');
[MGCBEf]=4.0e-3*bandpass(MGCBE,100,lo,hi,npo,npa,'butter');
[MGCBNf]=4.0e-3*bandpass(MGCBN,100,lo,hi,npo,npa,'butter');
%Decimate the 100 sps data:
TWKBEfd = interp1(timsTWKB,TWKBEf,timsLZB,'linear')';
TWKBNfd = interp1(timsTWKB,TWKBNf,timsLZB,'linear')';
MGCBEfd = interp1(timsMGCB,MGCBEf,timsLZB,'linear')';
MGCBNfd = interp1(timsMGCB,MGCBNf,timsLZB,'linear')';

% %Check up on interpolated data:
% figure
% plot(timsLZB,LZBEf);
% xlim([51545.5 51583]);
% figure
% plot(timsLZB,LZBNf);
% xlim([51545.5 51583]);
% figure
% plot(timsTWKB,TWKBEf);
% hold on
% plot(timsLZB,TWKBEfd,'r');
% xlim([51545.5 51583]);
% figure
% plot(timsTWKB,TWKBNf);
% hold on
% plot(timsLZB,TWKBNfd,'r');
% xlim([51545.5 51583]);
% figure
% plot(timsMGCB,MGCBEf);
% hold on
% plot(timsLZB,MGCBEfd,'r');
% xlim([51545.5 51583]);
% figure
% plot(timsMGCB,MGCBNf);
% hold on
% plot(timsLZB,MGCBNfd,'r');
% xlim([51545.5 51583]);

LZB=LZBEf+1i*LZBNf;
TWKB=TWKBEfd+1i*TWKBNfd;
MGCB=MGCBEfd+1i*MGCBNfd;
LZBrot=LZB*exp(1i*rotLZB);
TWKBrot=TWKB*exp(1i*rotTWKB);
MGCBrot=MGCB*exp(1i*rotMGCB);

timsTWKB=timsLZB-TWKBtoff;
timsMGCB=timsLZB-MGCBtoff;
%
% figure 
% plot(timsTWKB,real(TWKBrot),'k');
% xlim([5000 21000]);
%
%Check traces w/Armbruster:
figure %Individual station traces, 150s window
subplot(4,1,1); 
hold on
plot(timsLZB,real(LZBrot),'r');
plot(timsTWKB,real(TWKBrot),'b');
plot(timsMGCB,real(MGCBrot),'k');
xlim([timlook-75 timlook-37.5]);
%axis([timlook-75 timlook -0.2 0.2]);
subplot(4,1,2); 
hold on
plot(timsLZB,real(LZBrot),'r');
plot(timsTWKB,real(TWKBrot),'b');
plot(timsMGCB,real(MGCBrot),'k');
xlim([timlook-37.5 timlook]);
%axis([timlook timlook+75 -0.2 0.2]);
subplot(4,1,3); 
hold on
plot(timsLZB,real(LZBrot),'r');
plot(timsTWKB,real(TWKBrot),'b');
plot(timsMGCB,real(MGCBrot),'k');
xlim([timlook timlook+37.5]);
%axis([timlook timlook+75 -0.2 0.2]);
subplot(4,1,4); 
hold on
plot(timsLZB,real(LZBrot),'r');
plot(timsTWKB,real(TWKBrot),'b');
plot(timsMGCB,real(MGCBrot),'k');
xlim([timlook+37.5 timlook+75]);
%axis([timlook timlook+75 -0.2 0.2]);

if TWKBtoff > 0
    TWKBrot(1:tracelenLZB-TWKBsoff)=TWKBrot(TWKBsoff+1:tracelenLZB);
    TWKBrot(tracelenLZB-TWKBsoff+1:tracelenLZB)=0;
else
    TWKBrot(-TWKBsoff+1:tracelenLZB)=TWKBrot(1:tracelenLZB+TWKBsoff);
    TWKBrot(1:-TWKBsoff)=0;
end
if MGCBtoff > 0
    MGCBrot(1:tracelenLZB-MGCBsoff)=MGCBrot(MGCBsoff+1:tracelenLZB);
    MGCBrot(tracelenLZB-MGCBsoff+1:tracelenLZB)=0;
else
    MGCBrot(-MGCBsoff+1:tracelenLZB)=MGCBrot(1:tracelenLZB+MGCBsoff);
    MGCBrot(1:-MGCBsoff)=0;
end

% %Check up w/Armbruster:
% figure
% subplot(2,1,1); 
% plot(timsLZB,real(LZBrot),'r');
% hold on
% plot(timsLZB,real(TWKBrot),'b');
% plot(timsLZB,real(MGCBrot),'k');
% axis([51508 51545.5 -0.2 0.2]);
% subplot(2,1,2); 
% plot(timsLZB,real(LZBrot),'r');
% hold on
% plot(timsLZB,real(TWKBrot),'b');
% plot(timsLZB,real(MGCBrot),'k');
% axis([51560 51570 -0.2 0.2]);
% %axis([51545.5 51583 -0.2 0.2]);

%Pointwise x-correlation:
LZTWx=real(LZBrot).*real(TWKBrot);
LZMGx=real(LZBrot).*real(MGCBrot);
TWMGx=real(TWKBrot).*real(MGCBrot);
figure %Pointwise x-correlation, 150s window:
subplot(4,1,1); 
hold on
plot(timsLZB,LZMGx,'b');
plot(timsLZB,LZTWx,'r');
plot(timsLZB,TWMGx,'k'); 
xlim([timlook-75 timlook-37.5]);
subplot(4,1,2); 
hold on
plot(timsLZB,LZMGx,'b');
plot(timsLZB,LZTWx,'r');
plot(timsLZB,TWMGx,'k');
xlim([timlook-37.5 timlook]);
subplot(4,1,3); 
hold on
plot(timsLZB,LZMGx,'b');
plot(timsLZB,LZTWx,'r');
plot(timsLZB,TWMGx,'k');
xlim([timlook timlook+37.5]);
subplot(4,1,4); 
hold on
plot(timsLZB,LZMGx,'b');
plot(timsLZB,LZTWx,'r');
plot(timsLZB,TWMGx,'k');
xlim([timlook+37.5 timlook+75]);

LZTWxcum=cumsum(LZTWx);
LZMGxcum=cumsum(LZMGx);
TWMGxcum=cumsum(TWMGx);
figure %Running sum x-correlation pairs, with detections, all day
% plot(timsLZB,LZMGxcum-LZMGxcum(2057320),'r');
% hold on
% plot(timsLZB,LZTWxcum-LZTWxcum(2057320),'k');
% plot(timsLZB,TWMGxcum-TWMGxcum(2057320),'b');
plot(timsLZB,LZMGxcum,'b');
hold on
plot(timsLZB,LZTWxcum,'r');
plot(timsLZB,TWMGxcum,'k');
hrf = plotreflinesr(gca,detects,'x','k'); 
hrf = plotreflinesr(gca,Detects,'x','r'); 
hrf = plotreflinesr(gca,detects2_8,'x','y'); 
hrf = plotreflinesr(gca,Detects2_8,'x','g'); 
xlim([0 86400]);
%xlim ([51433 51583]);

LZTWxi=imag(LZBrot).*imag(TWKBrot);
LZMGxi=imag(LZBrot).*imag(MGCBrot);
TWMGxi=imag(TWKBrot).*imag(MGCBrot);
LZTWxicum=cumsum(LZTWxi);
LZMGxicum=cumsum(LZMGxi);
TWMGxicum=cumsum(TWMGxi);
figure %Running sum orthogonal x-correlation pairs, all day
plot(timsLZB,LZMGxicum,'b');
hold on
plot(timsLZB,LZTWxicum,'r');
plot(timsLZB,TWMGxicum,'k');

%Autocorrelation of ststaions:
LZBauto=real(LZBrot).*real(LZBrot);
LZB2=cumsum(LZBauto);
TWKBauto=real(TWKBrot).*real(TWKBrot);
TWKB2=cumsum(TWKBauto);
MGCBauto=real(MGCBrot).*real(MGCBrot);
MGCB2=cumsum(MGCBauto);
LZBautoi=imag(LZBrot).*imag(LZBrot);
LZB2i=cumsum(LZBautoi);
TWKBautoi=imag(TWKBrot).*imag(TWKBrot);
TWKB2i=cumsum(TWKBautoi);
MGCBautoi=imag(MGCBrot).*imag(MGCBrot);
MGCB2i=cumsum(MGCBautoi);

%Get 3-sec sums:
tmp=reshape(LZTWx,120,28800);
LZTWx3=sum(tmp);
tmp=reshape(LZMGx,120,28800);
LZMGx3=sum(tmp);
tmp=reshape(TWMGx,120,28800);
TWMGx3=sum(tmp);
tmp=reshape(LZBauto,120,28800);
LZBauto3=sum(tmp);
tmp=reshape(TWKBauto,120,28800);
TWKBauto3=sum(tmp);
tmp=reshape(MGCBauto,120,28800);
MGCBauto3=sum(tmp);
timsLZB3=1.5:3:86400-1.5;
figure %3-sec window ptwise x-correlation pairs and individual station strengths, +/-1000 sec 
subplot(2,1,1);
plot(timsLZB3,LZMGx3,'b');
hold on
plot(timsLZB3,LZTWx3,'r');
plot(timsLZB3,TWMGx3,'k');
hrf = plotreflinesr(gca,detects,'x','k'); 
hrf = plotreflinesr(gca,Detects,'x','r'); 
hrf = plotreflinesr(gca,detects2_8,'x','y'); 
hrf = plotreflinesr(gca,Detects2_8,'x','g'); 
%xlim([0 86400]);
xlim([timlook-1000 timlook+1000]);
%axis([10000 40000 -0.3 0.3]);
%axis([34500 37500 -0.3 0.3]);
subplot(2,1,2);
plot(timsLZB3,LZBauto3,'r');
hold on
plot(timsLZB3,TWKBauto3,'b');
plot(timsLZB3,MGCBauto3,'k');
hrf = plotreflinesr(gca,detects,'x','k'); 
hrf = plotreflinesr(gca,Detects,'x','r'); 
hrf = plotreflinesr(gca,detects2_8,'x','y'); 
hrf = plotreflinesr(gca,Detects2_8,'x','g'); 
%xlim([0 86400]);
axis([timlook-1000 timlook+1000 -0.5 2]);
%axis([10000 40000 -0.1 1]);
%axis([34500 37500 -0.1 1]);

figure %Running sum signal strength at each station, all day
plot(timsLZB,LZB2,'r');
hold on
plot(timsLZB,LZB2i,'r--');
plot(timsLZB,TWKB2,'b');
plot(timsLZB,TWKB2i,'b--');
plot(timsLZB,MGCB2,'k');
plot(timsLZB,MGCB2i,'k--');
hrf = plotreflinesr(gca,detects,'x','k'); 
hrf = plotreflinesr(gca,Detects,'x','r'); 
hrf = plotreflinesr(gca,detects2_8,'x','y'); 
hrf = plotreflinesr(gca,Detects2_8,'x','g'); 
xlim([0 86400]);



