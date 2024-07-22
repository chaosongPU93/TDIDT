function mfit = msfitvsiter(ampit,sigsta,greenf,zcrosses,ista,sps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mfit = msfitvsiter(ampit,sigsta,greenf,zcrosses,ista,sps)
%
% According to the arrival times and amplitudes of deconvolved sources 'ampit' 
% at each iteration, and templates of particular components (orthogonal or 
% vertical, optimal has been computed inside decon function), we can compute the
% predicted waveform for each iteration. Compared with data of particular 
% components, residual can be obtained, so is the variance. This code only
% returns the misfit as the variance of the residual at each iteration.
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/06/21
% Last modified date:   2024/06/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('ista',1);
defval('sps',160);

wlet = greenf(:,ista);  %template here is best aligned, tapered, linear trend removed, filtered
sig = sigsta(:,ista); %best aligned, filtered, tapered
lsig = length(sig);
nfft = lsig;
dt = 1/sps;  % sampling interval
ampitista=ampit{ista};
nit = size(ampitista,1);
twlet = zcrosses(ista)*dt;
itwlet = round(twlet/dt); % time shift of main arrival of wavelet in samples

pred = zeros(nfft,1);  % predefine the prediction array
mfit = zeros(nit,1);  %misfit, ie, variance of residual
for it = 1:nit
  impchg = zeros(nfft,1);   % the change of impulse resulting from the current iteration
  idximp = ampitista(it, 1);   %zero-crossing index and amp of decon impulse at this iter
  amp = ampitista(it, 2); %amp of impulse
  impchg(idximp) = amp;
  %change in predicted signal by convolving the change of impulse with wavelet
  predchg = conv(wlet, impchg, 'full');
  predchg = predchg(itwlet:nfft+itwlet-1);  % cut accordingly
  pred = pred + predchg;  % update the prediction
  res = sig-pred;   % residual
  mfit(it) = var(res);
end




