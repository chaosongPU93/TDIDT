%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test directivity effect
% to see if it is possible to fit the corner frequencies at all 4
% stations at KLNB, PGC, SSIB, SILB, from the result of 
% 'plt_tremor_spectra_ETS', with just the directivity effect.
%
% --This is just a simple inversion for rupture length, direction and 
% speed, so there is no fancy algorithm other than a plain iterative
% method to increase the stability or speed.
% --According to the result, it seems hard to ascribe the apprent
% 'corner' frequency difference at different stations to the directivity
% effect.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/11/10
% Last modified date:   2021/11/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear;
close all;
clc;

Cs = 4.6; % shear wave speed
k = 0.5; % velocity ratio between the rupture speed and shear wave speed
Vr = k*Cs; % rupture speed

azi = [2.4; 22.7; 22.3; 50.6];  % the azimuth wrt center of fam 002 at KLNB, PGC, SSIB, SILB
eta = [143; 143; 131; 140];  % take-off angle from center of fam 002
alpha = 28;
L = 0.95;  % rupture length, or length of the fault

fc_obs = [3.5; 3.5; 3.0; 3.8]; % observation data

term1 = (1-k.*cosd(alpha-azi).*sind(eta));  % just to simplify the expression
term2 = -k.*sind(alpha-azi).*sind(eta);  % just to simplify the expression

fc_pred = k.*Cs./L./term1; % predicted data
dfc = fc_obs-fc_pred; % difference between observation and prediction
dfc = reshape(dfc,[],1);

niter = 0;  % number of iteration

%store results of each iteration, this should be called iteration 0
Lit(niter+1) = L;
kit(niter+1) = k;
alphait(niter+1) = alpha;
dfcit(:, niter+1) = dfc;

while mean(abs(dfc)) > 1e-3 && niter <=20
  %update iteration number 
  niter = niter+1;
  
  %1st deriatives at previous model parameters 
  dfcdL = -k.*Cs./L.^2./term1;
  dfcdL = reshape(dfcdL, [],1);
  dfcdk = Cs./L./term1.^2;
  dfcdk = reshape(dfcdk, [],1);
  dfcdalpha = k.*Cs.*term2./L./term1.^2;
  dfcdalpha = reshape(dfcdalpha, [],1);  

%   dfc = [dfcdvrat dfcdalpha dfcdL]*[dvrat; dalpha; dL];
  
  %linearized matrix G
  G = [dfcdL dfcdk dfcdalpha];
  `
  %generalized inverse of G
%   Gmg = inv(G'*G)*G';
  Gmg = (G'*G)\G';
  
  %model change
  dm = Gmg*dfc;
  
  %update model
  L = L+dm(1);
  k = k+dm(2);
  alpha = alpha+dm(3);

  %update prediction and difference
  term1 = (1-k.*cosd(alpha-azi).*sind(eta));  % just to simplify the expression
  term2 = -k.*sind(alpha-azi).*sind(eta);  % just to simplify the expression
  fc_pred = k.*Cs./L./term1; % predicted data
  dfc = fc_obs-fc_pred; % difference between observation and prediction
  dfc = reshape(dfc,[],1);
  
  %update stored results of each iteration
  Lit(niter+1) = L;
  kit(niter+1) = k;
  alphait(niter+1) = alpha;
  dfcit(:, niter+1) = dfc;
  
end

%% plot
iter=6;
figure
subplot(4,1,1)
hold on;
box on;
grid on;
plot(0:iter, Lit(1:iter+1),'linew',2);
legend('rupture length L (km)','location','best');

subplot(4,1,2)
hold on;
box on;
grid on;
plot(0:iter, kit(1:iter+1),'linew',2);
legend('velocity ratio','location','best');

subplot(4,1,3)
hold on;
box on;
grid on;
plot(0:iter, alphait(1:iter+1),'linew',2);
legend('rupture direction (deg)','location','best');

subplot(4,1,4)
hold on;
box on;
grid on;
plot(0:iter, dfcit(1,1:iter+1),'linew',1);
plot(0:iter, dfcit(2,1:iter+1),'linew',1);
plot(0:iter, dfcit(3,1:iter+1),'linew',1);
plot(0:iter, dfcit(4,1:iter+1),'linew',1);
legend('model misfit at KLNB','model misfit at PGC','model misfit at SSIB','model misfit at SILB',...
  'location','best');













