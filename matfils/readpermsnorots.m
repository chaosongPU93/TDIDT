function [datE,datN,datZ,timsSTAperm]=readpermsnorots(prename,PERMSTA,idx,sps,lo,hi,npo,npa,fact)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
STAEdat=[prename,'.',PERMSTA(idx,1:3),'..BHE.D.SAC']; %HHE for polstas.
STANdat=[prename,'.',PERMSTA(idx,1:3),'..BHN.D.SAC'];
STAZdat=[prename,'.',PERMSTA(idx,1:3),'..BHZ.D.SAC'];
[STAE,HdrDataSTA,tnuSTA,pobjSTA,timsSTAperm]=readsac(STAEdat,0,'l');
[STAN,~,~,~,~]=readsac(STANdat,0,'l');
[STAZ,~,~,~,~]=readsac(STAZdat,0,'l');
tracelen=length(STAE);

%cosine taper over 1 second (at 40 sps) before filtering:
x=(0:pi/80:pi/2-pi/80)';
STAE(1:80)=0.;
STAE(81:120)=(1.-cos(x)).*STAE(81:120); %Only at start of day on E component!
STAN(1:40)=(1.-cos(x)).*STAN(1:40);
STAZ(1:40)=(1.-cos(x)).*STAZ(1:40);
x=flipud(x);
STAE(tracelen-39:tracelen)=(1.-cos(x)).*STAE(tracelen-39:tracelen);
STAN(tracelen-39:tracelen)=(1.-cos(x)).*STAN(tracelen-39:tracelen);
STAZ(tracelen-39:tracelen)=(1.-cos(x)).*STAZ(tracelen-39:tracelen);

%Filter data:
[STAEf]=fact*bandpass(STAE,40,lo,hi,npo,npa,'butter'); 
[STANf]=fact*bandpass(STAN,40,lo,hi,npo,npa,'butter'); 
[STAZf]=fact*bandpass(STAZ,40,lo,hi,npo,npa,'butter'); 
if sps == 40
    %no decimation
    STAEfd=STAEf;
    STANfd=STANf;
    STAZfd=STAZf;
else
    %upsample the 40 sps data to 100 (assume one or the other)
    STAEfd=resample(STAEf,5,2);
    STANfd=resample(STANf,5,2);
    STAZfd=resample(STAZf,5,2);
end
tracelen=length(STAEfd); %after decimation/upsampling

datE=STAEfd;
datN=STANfd;
datZ=STAZfd;
