function [opt,ort,nzeros,timsSTAperm]=readperms(prename,PERMSTA,PERMROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
STAEdat=[prename,'.',PERMSTA(idx,1:3),'..BHE.D.SAC']; %HHE for polstas.
STANdat=[prename,'.',PERMSTA(idx,1:3),'..BHN.D.SAC'];
[STAE,HdrDataSTA,tnuSTA,pobjSTA,timsSTAperm]=readsac(STAEdat,0,'l');
[STAN,~,~,~,~]=readsac(STANdat,0,'l');
tracelen=length(STAE);

STAsoffspsorig=round(PERMROTS(idx,4))  %input offsets given at 40 sps, same as PERM stations.
%find glitches
nzerosE=glitchesPERM(STAE,nwin,winlen,winoff,igstart,STAsoffspsorig);
nzerosN=glitchesPERM(STAN,nwin,winlen,winoff,igstart,STAsoffspsorig);
nzeros=max(nzerosE,nzerosN);

%cosine taper over 1 second (at 40 sps) before filtering:
x=(0:pi/80:pi/2-pi/80)';
STAE(1:80)=0.;
STAE(81:120)=sin(x).*STAE(81:120); %Only at start of day on E component!
STAN(1:40)=sin(x).*STAN(1:40);
x=flipud(x);
STAE(tracelen-39:tracelen)=sin(x).*STAE(tracelen-39:tracelen);
STAN(tracelen-39:tracelen)=sin(x).*STAN(tracelen-39:tracelen);

%Rather broadband filter before integrating:
[STAEf]=bandpass(STAE,40,0.2,15,npo,npa,'butter');
[STANf]=bandpass(STAN,40,0.2,15,npo,npa,'butter'); 
%integrate
STAE=inteFD(STAEf,0.025);
STAN=inteFD(STANf,0.025);
    
%Filter data:
[STAEf]=fact*bandpass(STAE,40,lo,hi,npo,npa,'butter'); 
[STANf]=fact*bandpass(STAN,40,lo,hi,npo,npa,'butter'); 
if sps == 40
    %no decimation
    STAEfd=STAEf;
    STANfd=STANf;
else
    %upsample the 40 sps data to 100 (assume one or the other)
    STAEfd=resample(STAEf,5,2);
    STANfd=resample(STANf,5,2);
end
tracelen=length(STAEfd); %after decimation/upsampling

%Rotate & split-correct
STA=STAEfd+1i*STANfd;
STAfastslow=STA*exp(-1i*PERMROTS(idx,2));
STAslow=real(STAfastslow);
STAfast=imag(STAfastslow);
len=length(STA);
off=round(10*sps/40);
STAslow(off:len-off)=STAslow(off+PERMROTS(idx,1):len-off+PERMROTS(idx,1));
STAsplitcorrected=(STAslow+1i*STAfast)*exp(1i*PERMROTS(idx,2));
STAscrot=STAsplitcorrected*exp(-1i*PERMROTS(idx,3));

%Timeshift
STAsoff=round(PERMROTS(idx,4)*sps/40);  %input offsets given at 40 sps. 
if STAsoff > -1
    STAscrot(1:tracelen-STAsoff)=STAscrot(STAsoff+1:tracelen);
    STAscrot(tracelen-STAsoff+1:tracelen)=0;
else
    STAscrot(-STAsoff+1:tracelen)=STAscrot(1:tracelen+STAsoff);
    STAscrot(1:-STAsoff)=0;
end

opt=real(STAscrot);
ort=imag(STAscrot);

