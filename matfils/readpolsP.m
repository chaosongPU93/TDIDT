function [optP]=readpolsP(prename,POLSTA,POLPROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

STAEdat=[prename,'.',POLSTA(idx,1:4),'..HHE.D.SAC']; %BHE for permstas.
STANdat=[prename,'.',POLSTA(idx,1:4),'..HHN.D.SAC'];
STAZdat=[prename,'.',POLSTA(idx,1:4),'..HHZ.D.SAC'];
[STAE,HdrDataSTA,tnuSTA,pobjSTA,timsSTApol]=readsac(STAEdat,0,'l');
[STAN,~,~,~,~]=readsac(STANdat,0,'l');
[STAZ,~,~,~,~]=readsac(STAZdat,0,'l');
tracelen=length(STAE); %temporary; before decimation

%cosine taper over 1 second (at 100 sps) before filtering:
x=(0:pi/200:pi/2-pi/200)';
STAE(1:100)=sin(x).*STAE(1:100);
STAN(1:100)=sin(x).*STAN(1:100);
STAZ(1:100)=sin(x).*STAZ(1:100);
x=flipud(x);
STAE(tracelen-99:tracelen)=sin(x).*STAE(tracelen-99:tracelen);
STAN(tracelen-99:tracelen)=sin(x).*STAN(tracelen-99:tracelen);
STAZ(tracelen-99:tracelen)=sin(x).*STAZ(tracelen-99:tracelen);

%Filter data:
[STAEf]=fact*bandpass(STAE,100,lo,hi,npo,npa,'butter');
[STANf]=fact*bandpass(STAN,100,lo,hi,npo,npa,'butter'); 
[STAZf]=fact*bandpass(STAZ,100,lo,hi,npo,npa,'butter'); 
if sps == 100
    %no decimation
    STAEfd=STAEf;
    STANfd=STANf;
    STAZfd=STAZf;
else
    %decimate the 100 sps data to 40 (assume one or the other)
    STAEfd=resample(STAEf,2,5);
    STANfd=resample(STANf,2,5);
    STAZfd=resample(STAZf,2,5);
end
tracelen=length(STAEfd); %after decimation

%Rotate 
incP=POLPROTS(:,1)*pi/180; %incidence angle
bazP=POLPROTS(:,2)*pi/180; %back-azimuth
%The following works ONLY if there is ONLY ONE staP:
m1=[cos(incP) -sin(incP)*sin(bazP) -sin(incP)*cos(bazP)];
% m2=[sin(incP) cos(incP)*sin(bazP) cos(incP)*cos(bazP)];
% m3=[0 -cos(bazP) sin(bazP)];
optP=m1(1)*STAZfd+m1(2)*STAEfd+m1(3)*STANfd;
%optP=STANfd; %Just a test; N works well for 068 SILB ...
%Q=m2(1)*STAz+m2(2)*STAe+m2(3)*STAn;
%T=m3(1)*STAz+m3(2)*STAe+m3(3)*STAn;

%STA=STAEfd+1i*STANfd;
%STAfastslow=STA*exp(-1i*POLROTS(idx,2));
%STAslow=real(STAfastslow);
%STAfast=imag(STAfastslow);
%len=length(STA);
%off=round(10*sps/40);
%STAslow(off:len-off)=STAslow(off+POLROTS(idx,1):len-off+POLROTS(idx,1));
%STAsplitcorrected=(STAslow+1i*STAfast)*exp(1i*POLROTS(idx,2));
%STAscrot=STAsplitcorrected*exp(-1i*POLROTS(idx,3));

%Timeshift
%The folloiwng works ONLY if there is ONLY ONE staP:
STAsoff=round(POLPROTS(:,3)*sps/40);  %input offsets given at 40 sps
if STAsoff > -1
    optP(1:tracelen-STAsoff)=optP(STAsoff+1:tracelen);
    optP(tracelen-STAsoff+1:tracelen)=0;
else
    optP(-STAsoff+1:tracelen)=optP(1:tracelen+STAsoff);
    optP(1:-STAsoff)=0;
end
