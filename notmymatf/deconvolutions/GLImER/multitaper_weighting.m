function [rf,rfss]=multitaper_weighting(rf_f,var_rf,tshift,dt)

% using output of matlab_multitaper which is the variance and the rf in the
% frequency domain, to do inverse-variance weighting and stacking
% (note: up to now without moveout correction)

[M,N]=size(rf_f);
rf_stack(1:N)=0;

for i=1:M
    for j=1:N
if var_rf(i,j)==0
   var_rf(i,j)=10^(-30); 
end
    end
end

%Continue with RF
for i=1:M
rf_tmp1(i,:)=rf_f(i,:).*((1./(var_rf(i,:))));
rf_tmp(i,:)=rf_f(i,:).*((1./(var_rf(i,:))))./sum(1./var_rf);
end

rf_stack=sum(rf_tmp1);
rf_stack=rf_stack./(sum(1./var_rf));

% frequency-domain taper
N2=floor(N/2)+1;
for j=1:M
for i=1:N2
fac=cos(pi/2*(i-1)/N2)^2;
rf_ff(i)=rf_stack(i)*fac;
rf_tmpp(j,i)=rf_tmp(j,i)*fac;
end
end

rfs=real(ifft([rf_ff,conj(fliplr(rf_ff(2:N2-1)))]));

for i=1:M
rf_tmp_new(i,:)=real(ifft([rf_tmpp(i,:),conj(fliplr(rf_tmpp(i,2:N2-1)))]));
end

for i=1:M
rf(i,1:N-(N-tshift/dt)) =rf_tmp_new(i,N-tshift/dt+1:N); 
rf(i,tshift/dt+1:N)=rf_tmp_new(i,1:N-tshift/dt);
end

%g=0;
rfss=rfs(N-tshift/dt+1:N); 
rfss(tshift/dt+1:N)=rfs(1:N-tshift/dt);
    
