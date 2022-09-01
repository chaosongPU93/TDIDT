function [ab]=convol(a,b,nf,ndt)
%convolution

afr=fft(a,nf);
bfr=fft(b,nf);

afr=afr.*bfr*ndt;

ab=real(ifft(afr,nf));

return
