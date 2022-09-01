%%% test different methods doing PSD estimate in matlab
close all
clear


%% warm-up, use periodogram
fs = 1000;  % sampling rate
t = 0:1/fs:1-1/fs;  % time, 1s
x = cos(2*pi*100*t) + sin(2*pi*150*t) + randn(size(t)); % composed signal

% Obtain the periodogram PSD estimate with 95%-confidence bounds. Plot the periodogram along with
% the confidence interval and zoom in on the frequency region of interest near 100 and 150 Hz.
[pxx,f,pxxc] = periodogram(x,rectwin(length(x)),length(x),fs,...
    'ConfidenceLevel',0.95);

plot(f,10*log10(pxx))
hold on
plot(f,10*log10(pxxc),'-.')

xlim([85 175])
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
title('Periodogram with 95%-Confidence Bounds')

%% power spectral density estimates using fft