function [stdec,wtdec] = specdiv_arraycond(st,wt,nt,dt,tzcwt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [stdec,wtdec] = specdiv_arraycond(st,wt,dt,tzcwt,nt)
%
% This function implements the deconvolution via spectral division in frequency
% domain, with regularization for stability, using the array-conditioned 
% method. It creates an optimal inverse filter that is automatic and independently
% determined for each frequency present in the Fourier transform (Chen et al. 2010)
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
%   nt: simply a constant value, or compute from a noise in time domain
% OUTPUT:
%   stdec: deconvolved signal, ie., impulses
%   wtdec: deconvolved wavelet itself, ie., a single impulse of amplitude 1 at the
%          zero-crossing of the wavelet
%
% Reference Literature:
%   Chen CW, Miller DE, Djikpesse HA, Haldorsen JBU, Rondenay S (2010) Array-
%     conditioned deconvolution of multiple-component teleseismic recordings.
%     Geophys J Int 182:967-976
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

%denorminator
denorm = Wf.* Wfconj;

%frequency-dependent damping, derived directly from noise 'nt', 
% nt = nt(1:nfft);
% [Nf, Nfconj] = time2freq(nt,length(nt));
[Nf, Nfconj] = time2freq(nt,nfft);
noise = Nf.* Nfconj;
denorm = denorm+noise;

%deconvolution in terms of frequency-domain division with regularization  
Sfdec = Sf.* Wfconj ./ denorm; % deconvolution of signal
Wfdec = Wf.* Wfconj ./ denorm; % deconvolution of wavelet itself

%I don't quite understand this part, Chao, 2022/01/13
%SR: this is the cos^2 filter discussed in Park and Levin, 2000,
%on page 1509 (top of right column), with fc set to 1 Hz - probably
%comes from Helffrich's code.
N2=floor(numel(Sfdec)/2)+1;
for j=1:N2
  fac=cos(pi/2*(j-1)/N2)^2;
  Sfdec(j,:)=Sfdec(j, :)*fac;
end

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