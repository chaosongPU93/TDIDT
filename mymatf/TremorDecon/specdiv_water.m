function [stdec,wtdec] = specdiv_water(st,wt,nt,dt,tzcwt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [stdec,wtdec] = specdiv_water(st,wt,dt,tzcwt,nt)
%
% This function implements the deconvolution via spectral division in frequency
% domain, with regularization for stability, using the water level method 
% (Menke, 1984). It is very similar in practice to damping factor deconvolution.
% Instead of adding a regularization parameter to every entry in the denominator,
% values below a specified threshold are replaced with a chosen water level value.
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
% Reference literature:
%   Menke, W (1984) Geophysical data analysis: discrete inverse theory. Academic Press, Inc., New York
%     City, N
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

%water level, could be a constant value or derived directly from noise 'nt', 
%depending on the input size
if length(nt) >1   % meaning that this is a noise trace
%   nt = nt(1:nfft);
%   [Nf, Nfconj] = time2freq(nt,length(nt));
  [Nf, Nfconj] = time2freq(nt,nfft);
  noise = Nf.* Nfconj;
  water = max(abs(noise));
else  % meaning that this is just a constant value
  water = nt;
end

%modify the denorminator according to the water level
denorm = Wf.* Wfconj;
denorm(real(denorm)<water)=water;

%deconvolution in terms of frequency-domain division with regularization  
Sfdec = Sf.* Wfconj ./ denorm; % deconvolution of signal
Wfdec = Wf.* Wfconj ./ denorm; % deconvolution of wavelet itself

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
