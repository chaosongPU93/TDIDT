function [ab]=gfilter(a,b,nf,ndt)

afr=fft(a,nf);
afr=afr.*b*ndt;

ab=real(ifft(afr,nf));

return
