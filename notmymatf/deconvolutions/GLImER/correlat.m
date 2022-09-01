function [x]=correlat(u,w,nf)

x=real(ifft(fft(u,nf).*conj(fft(w,nf)),nf));

return
