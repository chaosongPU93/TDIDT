function [lrf,qrf,eps]=timedeco_noise(w,u,ndt,tshift)

% Implementation of the time deconvolution described by Gurrola 1995
% "Simultaneous time-domain deconvolution with application  to  the
% computation of receiver functions".
% w: denominator, u: nominator, eps: weighting/regularisaton parameter, 
% ndt: sampling interval, tshift: time before direct wave onset
% The following function is needed: glev.m

N=length(w);
M=length(u);
N2=2*N-1;
%number of points in fft
nf=2^nextpow2(N);

w=w';

%noise; regularisation parameter eps is set to the maximum of the pre-event noise
noise=w(1:(tshift-2)/ndt);
noiseM=convmtx(noise,length(noise));
noised=noiseM'*noiseM;
eps=max(max(noised));

u01=zeros(tshift/ndt,1);
u02=zeros(N-tshift/ndt-1,1);
u2=[u01; u'; u02];
w2=[u01; w ; u02];
%creating convolution matrix
CM=convmtx(w,N); 
dumm1=(eye(N)+1/eps*(CM'*CM));

%Calculate rf
dummq=1/eps*(CM'*u2);
qrff=glev(dumm1,dummq);
qrf=qrff(1:N);
dumml=1/eps*(CM'*w2);
lrff=glev(dumm1,dumml);
lrf=lrff(1:N);


