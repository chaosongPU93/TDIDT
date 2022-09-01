% test the equation 3-4 in Ide2019JGR

clc
clear
close all

Fs = 100;
t2 = 0: 1/Fs: 3;
t2 = reshape(t2, [], 1);
taus = 0.5;
macc = (1-2*pi/taus .*t2) .* exp(-2*pi/taus .*t2);
macc = macc/max(macc);
macc = reshape(macc, [], 1);

t1 = -3: 1/Fs: 0-1/Fs;
t1 = reshape(t1, [], 1);
% y = dirac(t2);
% idx = y == inf;
% y(idx) = 1;
% plot(t2,y);

t = [t1; t2];
macc = [zeros(length(t1),1); macc];

lo = 1;
hi = 8;
npo = 2;
npa = 1;
maccbp2o1a = Bandpass(macc, Fs, lo, hi, npo, npa, 'butter');
figure
plot(t,macc,'k-','linew',2); hold on
plot(t,maccbp2o1a,'r','linew',2);
title(sprintf('%d-%d Hz, npo: %d, npa: %d', lo,hi,npo,npa));
ylim([-0.2 1.2]);


lo = 1;
hi = 8;
npo = 1;
npa = 2;
maccbp1o2a = Bandpass(macc, Fs, lo, hi, npo, npa, 'butter');
figure
plot(t,macc,'k-','linew',2); hold on
plot(t,maccbp1o2a,'r','linew',2);
title(sprintf('%d-%d Hz, npo: %d, npa: %d', lo,hi,npo,npa));
ylim([-0.2 1.2]);


lo = 1;
hi = 8;
npo = 2;
npa = 2;
maccbp2o2a = Bandpass(macc, Fs, lo, hi, npo, npa, 'butter');
figure
plot(t,macc,'k-','linew',2); hold on
plot(t,maccbp2o2a,'r','linew',2);
title(sprintf('%d-%d Hz, npo: %d, npa: %d', lo,hi,npo,npa));
ylim([-0.2 1.2]);


lo = 2;
hi = 8;
npo = 2;
npa = 2;
maccbp2o2a = Bandpass(macc, Fs, lo, hi, npo, npa, 'butter');
figure
plot(t,macc,'k-','linew',2); hold on
plot(t,maccbp2o2a,'r','linew',2);
title(sprintf('%d-%d Hz, npo: %d, npa: %d', lo,hi,npo,npa));
ylim([-0.2 1.2]);



