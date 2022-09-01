function [datE,datN,datZ,timsSTApol]=readpolsnorots(prename,POLSTA,idx,sps,lo,hi,npo,npa,fact)
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
STAE(1:100)=(1.-cos(x)).*STAE(1:100);
STAN(1:100)=(1.-cos(x)).*STAN(1:100);
STAZ(1:100)=(1.-cos(x)).*STAZ(1:100);
x=flipud(x);
STAE(tracelen-99:tracelen)=(1.-cos(x)).*STAE(tracelen-99:tracelen);
STAN(tracelen-99:tracelen)=(1.-cos(x)).*STAN(tracelen-99:tracelen);
STAZ(tracelen-99:tracelen)=(1.-cos(x)).*STAZ(tracelen-99:tracelen);

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
tracelen=length(STAEfd); %after decimation/upsampling

datE=STAEfd;
datN=STANfd;
datZ=STAZfd;
