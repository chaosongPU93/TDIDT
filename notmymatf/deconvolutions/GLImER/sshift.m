function [b]=sshift(a,nf,ndt,shift)
%

afr=fft(a,nf);

s=round(shift/ndt);
p=2*pi*(1:nf).*s/nf;
afr=afr.*(cos(p)-1i.*sin(p));

%zurueck in time domain
b=real(ifft(afr,nf))/cos(2*pi*s/nf);

return