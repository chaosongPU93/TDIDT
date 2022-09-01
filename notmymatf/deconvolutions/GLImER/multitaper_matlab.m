function [rf,lrf,var_rf,tmpp] = multitaper_matlab(P,D,dt,tshift,regul)

% Findings: when P arrival is not in the center of the window, the
% amplitudes are not unity at the beginning and decreasing from there on.
% Instead they peak at the time shift which corresponds to the middle index
% in the P time window.

% TB  = time bandwidth product (usually between 2 and 4)
% NT  = number of tapers to use, has to be <= 2*TB-1
% TB = 4; NT = 7; %choice of TB = 4, NT = 3 is supposed to be optimal
% t0 = -5; t1 = max(time); % This defines the beginning and end of the lag
% times for the receiver function

%%% changed by KS June 2016, added coherence and variance estimate; output
%%% RF in time domain (rf), RF in freq. domain (tmpp), and variance of RF
%%% (var_rf); can be used as input for multitaper_weighting.m
%%% the input "regul" defines which type of regularization is used,
%%% regul='fqd' defines frequency dependent regularisation from pre-event
%%% noise, regul='con' defines adding a constant value (here maximum of
%%% pre-event noise) as regularisation

%wavelet always in the center of the window
win_len = tshift*2;

Nwin = round(win_len/dt);

% Fraction of overlap overlap between moving time windows. As your TB
% increases, the frequency smearing gets worse, which means that the RFs
% degrate at shorter and shorter lag times. Therefore, as you increase TB,
% you should also increase Poverlap.
Poverlap = 0.75; 

% Length of waveforms;
nh = length(P);

% Create moving time windowed slepians
starts = 1:round((1-Poverlap)*Nwin):nh-Nwin+1;

%tapernumber, bandwith
%%% is in general put to NT=3, and bandwidth to 2.5
TB=2.5;
%NT=2*TB-1; %4 tapers
NT=3;

% Construct Slepians
[Etmp,lambdas] = dpss(Nwin,TB,NT);
E = zeros(length(starts)*NT,nh);
n = 0;

NUM = zeros(size(P));  DEN = zeros(size(D));
DUM = zeros(size(D));

%% finding frequency dependent regularisation parameter DEN_noise
% added: KS 26.06.2016
Pn=zeros(size(P));

Pn(3/dt:(tshift-5)/dt)=P(3/dt:(tshift-5)/dt); 
% pre-event noise: starting 3s after trace start 
% stop 10s before theoretical start of P
% wave to aviod including it

DEN_noise = zeros(size(P));

%% Multitaper
for j = 1:length(starts)
    for k = 1:NT
        n = n + 1;
        E(n,starts(j):starts(j)+Nwin-1) = transpose(Etmp(:,k));
        
        tmp1 = fft(E(k,:).*P); 
        tmp2 = fft(E(n,:).*D); % allow moving time windows

       
        NUM = NUM + lambdas(k)*conj(tmp1).*tmp2;
        DEN = DEN + lambdas(k)*conj(tmp1).*tmp1;
        
        % DUM only from D trace (converted wave component) used in
        % coherence estimate
        DUM = DUM + lambdas(k)*conj(tmp2).*tmp2;
        
        % pre-event noise        
        tmp1n = fft(E(n,:).*Pn); % always stick to first time window
        DEN_noise = DEN_noise + lambdas(k)*conj(tmp1n).*tmp1n;
        
    end
end


%% Calculate optimal RF with frequency-dependend regularisation
freqdep=strcmp(regul,'fqd');
const=strcmp(regul,'con');

if freqdep==1
 tmpp=(NUM)./((DEN+DEN_noise));
 tmpp_l=DEN./(DEN+DEN_noise);

 N2=floor(nh/2)+1;
for i=1:N2
fac=cos(pi/2*(i-1)/N2)^2;
tmpp(:,i)=tmpp(:,i)*fac;
end

end

if const==1
% ordinary regularisation with adding only a constant value
 eps=(max(real(DEN_noise)));
 tmpp=NUM./(DEN + eps);
 tmpp_l=DEN./(DEN + eps);
end

if (freqdep==0) && (const==0)
   error(['Regularization not defined (your input: regul=',regul ...
       ').' ...
       ' Use either "fqd" for frequency-dependent' ...
       'or "con" for constant value regularization.'])
end
%%% RF without variance weighting
tmp1 = real(ifft(tmpp));
tmp1_l = real(ifft(tmpp_l));

% Interpolate to desired
N=length(P);
rf=tmp1(N-tshift/dt+1:N); 
rf(tshift/dt+1:N)=tmp1(1:N-tshift/dt);

lrf=tmp1_l(N-tshift/dt+1:N); 
lrf(tshift/dt+1:N)=tmp1_l(1:N-tshift/dt);

%%%%
%%% Coherence and Variance of RF
% added: KS 26.06.2016
C_rf=((NUM))./sqrt((DUM.*DEN));

for ii=1:length(C_rf)
var_rf(ii) = ((1-abs(C_rf(ii))^2)/((NT-1)*abs(C_rf(ii))^2))*(abs(tmpp(ii))^2);
end

end
