function [stdec,wtdec]=specdiv_damp(st,wt,nt,dt,tzcwt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [stdec,wtdec]=specdiv_damp(st,wt,dt,tzcwt,nt)
%
% This function implements the deconvolution via spectral division in frequency
% domain, with regularization for stability using constant damping factor. 
% From Wiener deconvolution, the damping factor effectively prewhitens the power
% spectrum in the denominator and prevents the instability that would be caused
% by small values. Equivalent to adding white noise to the signal and effectively
% shifts the power spectrum of the denominator (auto-CC of the wavelet) vertically.
%
% According to the cross-correlation theorem : the cross-correlation 
% between two signals is equal to the product of fourier transform of one
% signal multiplied by complex conjugate of fourier transform of another
% signal. After doing this, when we take the ifft of the product signal, 
% we get a peak which indicates the shift between two signals.
%
%
% INPUT:
%   st: signal in time domain
%   wt: wavelet in time domain
%   dt: sampling interval in [sec]
%   tzcwt: time in [sec] of the zero-crossing of the wavelet in order to shift
%   nt: simply a constant damping factor, or compute from a noise in time domain
% OUTPUT:
%   stdec: deconvolved signal, ie., impulses
%   wtdec: deconvolved wavelet itself, ie., a single impulse of amplitude 1 at the
%          zero-crossing of the wavelet
%
% Reference literature:
%   Wiener N (1964) Extrapolation, interpolation, and smoothing of stationary time 
%     series: with engineering applications. MIT Press, Cambridge, MA
%   Gurrola H, Baker GE, Minster JB (1995) Simultaneous time-domain deconvolution
%     with application to the computation of receiver functions. Geophys J Int 
%     120:537-543
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
[Sf, ~] = time2freq(st,nfft);
[Wf, Wfconj] = time2freq(wt,nfft);

% [Wgtf, ~] = time2freq(wgt,nfft);
% Sf = conv(Sf,Wgtf,'same');
% Sf = Sf.*Wgtf;

%damping factor, could be a constant value or derived directly from noise 'nt', 
%depending on the input size
if length(nt) >1   % meaning that this is a noise trace
%   nt = nt(1:nfft);
%   [Nf, Nfconj] = time2freq(nt,length(nt));
  [Nf, Nfconj] = time2freq(nt,nfft);
  noise = Nf.* Nfconj;
  damp = max(abs(noise));
else  % meaning that this is just a constant value
  damp = nt;
end

%deconvolution in terms of frequency-domain division with regularization
Sfdec = Sf.* Wfconj ./ (Wf.* Wfconj + damp); % deconvolution of signal
Wfdec = Wf.* Wfconj ./ (Wf.* Wfconj + damp); % deconvolution of wavelet itself

%rototation in complex domain to accomondate the time shift in [sec]
w=(0:nfft/2)*2*pi/(nfft*dt);  % angular frequencies, note that the wrapping nature of FFT
w=[w,-fliplr(w(2:end-1))];
w=reshape(w,[],1);
Sfdec=Sfdec.*exp(-1i*w*(tzcwt-dt));  % the actual time to shift is minus 1 sample
Wfdec=Wfdec.*exp(-sqrt(-1)*w*(tzcwt-dt));

%transform back to time domain
stdec=real(ifft(Sfdec,nfft));
stdec=stdec(1:lst); % deconvoluted signal
wtdec=real(ifft(Wfdec,nfft));
wtdec=wtdec(1:lst); % deconvoluted wavelet, should be just a single impulse of amplitude 1




