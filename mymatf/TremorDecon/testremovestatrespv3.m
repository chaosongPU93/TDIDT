%testremovestatrespv3.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a follow-up script to 'testremovestatrespv2.m', in order to compare 
% the processed data at the same station that used different instrument 
% response files instrument response removal (IRR).
% 
% From 'testremovestatrespv2.m', IRR22 using RESP at LZB in 2003 looks very 
% different from IRR19 using SAC_PZ, I want to know why.
% Here i compare the IRR22 using RESP, IRR22 using PZ, and IRR19, to see if
% it is due to the intrinsic difference between SAC_PZ and RESP file at this
% station.
%
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/06/16
% Last modified date:   2022/06/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
clear

%%
workpath = getenv('ALLAN');
pz19fn = strcat(workpath,'/data-no-resp/arch2003/MAR/2003.061.00.00.00.0000.CN.LZB..BHZ.D.SAC');
pz22fn = strcat(workpath,'/data-no-resp/testdiffresp/lzbrespvspz/2003.061.00.00.00.0000.CN.LZB..BHZ.D.SAC');
resp22fn = strcat(workpath,'/data-no-resp/2003/MAR/2003.061.00.00.00.0000.CN.LZB..BHZ.D.SAC');

sps=40;     % samples per second
hi=6.5;    % frequency band
lo=1.25;
npo=2;     % poles, passes of filters
npa=2;

[pz19,Hdr2,~,~,t2]=readsac(pz19fn,0,'l');
[pz19]=Bandpass(pz19,40,lo,hi,npo,npa,'butter');
[pz22,Hdr3,~,~,t3]=readsac(pz22fn,0,'l');
[pz22]=Bandpass(pz22,40,lo,hi,npo,npa,'butter');
[resp22,Hdr4,~,~,t4]=readsac(resp22fn,0,'l');
[resp22]=Bandpass(resp22,40,lo,hi,npo,npa,'butter');

figure
subplot(411)
plot(t3,pz22,'r-'); hold on;
plot(t2,pz19,'k-'); 
text(0.95,0.8,'IRR pz22 vs. IRR pz19','Units','normalized','HorizontalAlignment','right');
% title('LZB, E, 0.04-20 Hz');
% title('LZB, N, 0.04-20 Hz');
title('LZB, Z, 0.04-20 Hz');
xlabel('Time (s)');
ylabel('Amplitude (nm/s)');

subplot(412)
plot(t4,resp22,'r-'); hold on;
plot(t2,pz19,'k-'); 
text(0.95,0.8,'IRR resp22 vs. IRR pz19','Units','normalized','HorizontalAlignment','right');

subplot(413)
plot(t4,resp22,'r-'); hold on;
plot(t3,pz22,'k-');
text(0.95,0.8,'IRR resp22 vs. IRR pz22','Units','normalized','HorizontalAlignment','right');

subplot(414)
plot(t3,resp22-pz22,'r-'); hold on;
plot(t2,pz19-pz22,'k-');
ylabel('difference');


