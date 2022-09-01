function [lrf,qrf]=spectraldivision(w,u,ndt,tshift,regul)

% Function spectraldivision(w,u,ndt,tshift,regul) is a standard spectral
% division frequency domain deconvolution. The regularization can be chosen
% by the variable "regul", this can be 'con', 'wat', or 'fqd' for constant
% damping factor, waterlevel, or frequency-dependent damping, respectively.
% "ndt" is the sampling rate in seconds and "tshift" is the time before
% onset of the direct wave.
% w=lcomp;
% u=qcomp;
%
% output is the receiver function "qrf" and the direct wave deconvolved by
% itself

N=length(w);

% pre-event noise needed for regularisation parameter
wn=zeros(1,N);
wn(1:(tshift-3)/ndt)=w(1:(tshift-3)/ndt);

%number of points in fft
nf=2^nextpow2(N);
nft=nf/2+1;

%fourier transform
uf=fft(u,nf);
wf=fft(w,nf);
wnf=fft(wn,nf);

%denominator and regularisation
den=wf.*conj(wf);
noise=wnf.*conj(wnf);

den0=den;

% which regularization do you want?
freqdep=strcmp(regul,'fqd');
const=strcmp(regul,'con');
water=strcmp(regul,'wat');

if (freqdep==0) && (const==0) && (water==0)
   error(['Regularization not defined (your input: regul=',regul ...
       ').' ...
       ' Use either "fqd" for frequency-dependent' ...
       'or "con" for constant value regularization' ...
       'or "wat" for water-level.'])
end

% constant damping factor regularization
if const==1
   eps=max(real(noise));
   den=den+eps; 
end

% frequency-dependent damping 
if freqdep==1
   den=den+noise;
end

% waterlevel regularization
if water==1
eps=(max(real(noise)));
den(real(den)<eps)=eps;
end

% numerator
num=uf.*conj(wf);
numl=wf.*conj(wf);

% deconvolution
rfq=num./den;
rfl=numl./den;

if freqdep==1
     N2=floor(numel(rfq)/2)+1;
for i=1:N2
fac=cos(pi/2*(i-1)/N2)^2;
rfq(:,i)=rfq(:,i)*fac;
end

end

% back transformation
w=(0:nf/2)*2*pi/(nf*ndt);
w=[w,-fliplr(w(2:end-1))];
rfq=rfq.*exp(-sqrt(-1)*w*tshift);
rfl=rfl.*exp(-sqrt(-1)*w*tshift);

qrft=ifft(rfq,nf);
qrf=real(qrft);

qrf=qrf(1:N);
lrft=ifft(rfl,nf);
lrf=real(lrft(1:N));
lrf=lrf(1:N);
