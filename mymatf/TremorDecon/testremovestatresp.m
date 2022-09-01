%testremovestatresp.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is to compare the data at the same station that unders 
% different process of instrument response removal (IRR). We will compare
% the raw data, IRR using 2019 SAC_PZs, IRR using 2019 RESPs, IRR using 
% 2021/2017/2016 RESPs. Check '/home/data2/chaosong/Notes/20220608.md' for
% how the RESPs are obtained, and '/home/data2/chaosong/matlab/allan/
% data-no-resp/testdiffresp.sh' for the IRR process. 
% 
%
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/06/09
% Last modified date:   2022/06/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
clear

%%
workpath = getenv('ALLAN');
rawfn = strcat(workpath,'/2003/MAR/2003.062.00.00.00.0000.CN.SSIB..HHE.D.SAC');
pz19fn = strcat(workpath,'/data-no-resp/2003/MAR/2003.062.00.00.00.0000.CN.SSIB..HHE.D.SAC');
resp19fn = strcat(workpath,'/data-no-resp/testdiffresp/ver2019/2003.062.00.00.00.0000.CN.SSIB..HHE.D.SAC');
resp21fn = strcat(workpath,'/data-no-resp/testdiffresp/ver2021/2003.062.00.00.00.0000.CN.SSIB..HHE.D.SAC');

[raw,Hdr1,~,~,t1]=readsac(rawfn,0,'l');
[pz19,Hdr2,~,~,t2]=readsac(pz19fn,0,'l');
[resp19,Hdr3,~,~,t3]=readsac(resp19fn,0,'l');
[resp21,Hdr4,~,~,t4]=readsac(resp21fn,0,'l');

figure
subplot(411)
plot(t2,pz19,'k-'); hold on;
plot(t1,raw,'r-');
text(0.95,0.8,'IRR pz19 vs. raw','Units','normalized','HorizontalAlignment','right');
title('SSIB, E, 0.04-20 Hz');
% title('SSIB, N, 0.04-20 Hz');
% title('SSIB, Z, 0.04-20 Hz');
xlabel('Time (s)');
ylabel('Amplitude (nm/s)');

subplot(412)
plot(t2,pz19,'k-'); hold on;
plot(t3,resp19,'r-');
text(0.95,0.8,'IRR pz19 vs. IRR resp19','Units','normalized','HorizontalAlignment','right');

subplot(413)
plot(t4,resp21,'k-'); hold on;
plot(t3,resp19,'r-');
text(0.95,0.8,'IRR resp21 vs. IRR resp19','Units','normalized','HorizontalAlignment','right');

subplot(414)
plot(t2,pz19-resp19,'k-'); hold on;
plot(t3,resp21./resp19,'r-');
ylabel('Ratio or difference');




