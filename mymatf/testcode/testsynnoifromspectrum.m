%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test if you can recover the time-domain signal from its amplitude and phase
% sepctra through the usage of 'fft' and 'ifft'.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/01/13
% Last modified date:   2022/01/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% close all
clc

fs = 100;                                % sample frequency (Hz)
t = 0:1/fs:10-1/fs;                      % 10 second span time vector
x = (1.3)*sin(2*pi*15*t) ...             % 15 Hz component
  + (1.7)*sin(2*pi*40*(t-2)) ...         % 40 Hz component
  ;

nfft = length(x);          % number of samples

[xf,f,amp,pha,power,psd] = fftspectrum(x, nfft, fs,'twosided');

%direct ifft
xrec = ifft(xf,nfft);

%synthesize from amp and phase spectra
XF = amp*nfft.*exp(1i*pha);
xrec1 = ifft(XF,nfft);

%uniform, random phase with the same span [-pi,pi];
pharand = (rand(1,nfft)-0.5)*2*pi;
% pharand = (randn(1,nfft)-0.5)*2*pi;

%construct record with the same amplitude but random phase
XFrand = amp*nfft.*exp(1i*pharand);
xrecrand = real(ifft(XFrand,nfft));


%%
% figure
% plot(pha/pi,'b'); hold on
% plot(pharand/pi),'r';
% ylim([-1.5 1.5]);

figure
plot(x,'b'); hold on
plot(xrecrand,'r');

% %%
% figure
% plot(x); hold on
% plot(xrec);
% 
% %%
% figure
% plot(x); hold on
% plot(xrec1);

%%

