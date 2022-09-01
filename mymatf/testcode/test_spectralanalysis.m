%%% This is the test for comparing power density estimates using fft or
%%% periodogram.
%%% SUMMARY:
%%%     Periodogram can obtain PSD much easier than fft!
%%%

%% Even-Length Input with Sample Rate
% clear;
% close all;

% an even-length signal sampled at 1 kHz
% Create a signal consisting of a 100 Hz sine wave in N(0,1) additive noise.
% The signal length is 1000 samples
Fs = 1000;  
t = 0:1/Fs:1-1/Fs;  % can regard the time as only 1 s long
x = cos(2*pi*100*t) + randn(size(t));   

% Obtain the periodogram using fft. The signal is real-valued and has even
% length. Because the signal is real-valued, you only need power estimates
% for the positive or negative frequencies.
% In order to conserve the total power, multiply all frequencies that occur
% in both sets -- the positive and negative frequencies -- by a factor of 2.
% Zero frequency (DC) and the Nyquist frequency do not occur twice. Plot
% the result.
N = length(x);
xdft = fft(x);  % dft of x
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;   % psd
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/N:Fs/2;

figure
plot(freq,10*log10(psdx))
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')

% comparison
figure
periodogram(x,rectwin(N),N,Fs, 'psd')

mxerr = max(psdx'-periodogram(x,rectwin(N),N,Fs, 'psd'))

%% input with Normalized Frequency
% Use fft to produce a periodogram for an input using normalized frequency.
% Create a signal consisting of a sine wave in N(0,1) additive noise. The
% sine wave has an angular frequency of pi/4 rad/sample. 
clear;
close all;

n = 0:999;
x = cos(pi/4*n) + randn(size(n));
N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(2*pi*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:(2*pi)/N:pi;

figure
plot(freq/pi,10*log10(psdx))
grid on
title('Periodogram Using FFT')
xlabel('Normalized Frequency (\times\pi rad/sample)') 
ylabel('Power/Frequency (dB/rad/sample)')

figure
periodogram(x,rectwin(length(x)),length(x))

mxerr = max(psdx'-periodogram(x,rectwin(length(x)),length(x)))

%% Complex-Valued Input with Normalized Frequency
% Use fft to produce a periodogram for a complex-valued input with
% normalized frequency. The signal is a complex exponential with an angular
% frequency of pi/4 rad/sample in complex-valued N(0,1) noise.

clear;
close all;
n = 0:999;
x = exp(1j*pi/4*n) + [1 1j]*randn(2,length(n))/sqrt(2);

% Use fft to obtain the periodogram. Because the input is complex-valued,
% obtain the periodogram from (-pi, pi] rad/sample. Plot the result.
N = length(x);
xdft = fft(x);
psdx = (1/(2*pi*N)) * abs(xdft).^2;
freq = 0:(2*pi)/N:2*pi-(2*pi)/N;

figure
plot(freq/pi,10*log10(psdx))
grid on
title('Periodogram Using FFT')
xlabel('Normalized Frequency (\times\pi rad/sample)') 
ylabel('Power/Frequency (dB/rad/sample)')

figure
% use two-sided option!
periodogram(x,rectwin(length(x)),length(x),'twosided')

mxerr = max(psdx'-periodogram(x,rectwin(length(x)),length(x),'twosided'))


 





