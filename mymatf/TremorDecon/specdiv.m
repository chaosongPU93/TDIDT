function stdec = specdiv(st,wt,dt,tzcwt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function sigdecon = specdiv(sig, wlet)
%
% This function implements the most basic deconvolution via the plain spectral
% division in frequency domain, without any regularization for stability. It
% comes naturally with the risk of instability via amplifying noise (near-zero)
% because of a small value on the denominator.
%
% INPUT:
%   st: signal in time domain
%   wt: wavelet in time domain
%   dt: sampling interval in [sec]
%   tzcwt: time in [sec] of the zero-crossing of the wavelet in order to shift
% OUTPUT:
%   stdec: deconvolved signal, ie., impulses
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/01/13
% Last modified date:   2022/01/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of FFT
lst = length(st);
lwt = length(wt);
nfft = pow2(nextpow2(max(lst, lwt)));

%transform the signal and wavelet from time domain to frequency domain
%express the time-domain of signal and wavelet in frequency domain
[Sf,~] = time2freq(st,nfft);
[Wf,~] = time2freq(wt,nfft);

%direct division
Sfdec = Sf./Wf;

%rototation in complex domain to accomondate the time shift in [sec]
w=(0:nfft/2)*2*pi/(nfft*dt);  % angular frequencies, note that the wrapping nature of FFT
w=[w,-fliplr(w(2:end-1))];
w=reshape(w,[],1);
Sfdec=Sfdec.*exp(-1i*w*(tzcwt-dt));  % the actual time to shift is minus 1 sample

%transform back to time domain
stdec = real(ifft(Sfdec));
stdec = stdec(1:lst); % deconvoluted signal





