function [m0dot,time]=BSE1d_Ide(C0,dt,alpha,sigma,Vslip,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [m0dot,time]=BSE1d_Ide(C0,dt,alpha,sigma,Vslip,mu)
%
% This function is to generate moment rate function wrt time, for the 
% Brownian Slow Earthquake model proposed by Satoshi Ide in 2008 and further 
% developed. It assumes a Brown motion in the radius of the fault area, such
% that a stochastic differential equation is used to describe r. To prevent
% r goes to infinity, a damping term proportional to r is needed in the SDE.
% --We start with eq1 of ide2008, assuming the initial state of radius is 0.
% Then the fault area S and moment rate M0dot is directly obtain from radius.
% --Note that if you start with eq2 in ide&yabe2019, then dS can be complex
% numbers, as S can be negative. r can also be negative, but if you r^2, 
% there is no practical problem.
% --default values for a few quantities are chosen from ide2008 and 
% ide&yabe2019
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/11/10
% Last modified date:   2023/11/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% clear
% clc
% close all

defval('C0',pi); %geomerty constant, from ide&yabe2019
defval('dt',0.1); %time step [s], from ide&yabe2019
defval('alpha',1/100); %damping coefficient [s^-1], 1/alpha is characteristic damping time, from ide&yabe2019
defval('sigma',400); %diffusion coefficient, [m/s^0.5], from ide2008
defval('Vslip',5e-6); %shear modulus, [Pa], from ide2008
defval('mu',4e10); %slip velocity, [m/s], from ide2008

% rng('default');
meandB = 0;  %mean of Gaussian white noise
vardB = dt;  %variance of Gaussian white noise
stddB = sqrt(vardB);
npts = round(10000/dt);
%Gaussian white noise with zero mean and a variance of time increment dt
dBt = stddB.*randn(npts,1) + meandB; %create normal random numbers following N(ave, vari)

%%%%%%%%%%%% if start with eq1 of ide&yabe2019
rt = zeros(npts+1,1); %radius of fault area, and take initial state as 0
for i = 1: npts
  %Stochastic Differential Equation of radius
  drt = -alpha*rt(i)*dt + sigma*dBt(i); %chang of rt, from eq1 of ide2008
  rt(i+1) = rt(i)+drt;  
end
S = C0*rt.^2; %fault area
%%%%%%%%%%%% if start with eq1 of ide&yabe2019

% %%%%%%%%%%%% if start with eq2 of ide&yabe2019
% S = zeros(npts+1,1);  %slip area as a function of time
% % r0 = 4e4; %from ide2008, initial radius of a circular fault at time 0
% % S(1) = C0*r0^2;  %initial slip area at time 0
% for i = 1: npts
%   dS = (C0*sigma^2 - 2*alpha*S(i))*dt + 2*sigma*sqrt(C0*S(i))*dBt(i);
% %   S(i+1) = S(i)+dS;
%   S(i+1) = S(1)+dS;
% %   if S(i+1)<0
% %     break
% %   end
% end
% S = real(S);
% %%%%%%%%%%%% if use eq2 of ide&yabe2019

m0dot = mu*Vslip*S; %moment rate
time = (0:npts)*dt;

%%
f = initfig(8,4,1,1);
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
plot(ax,time,m0dot);
text(ax,0.98,0.9,...
  sprintf('\\alpha=%.1e; \\sigma=%.1e m/s^{0.5}; V_{slip}=%.1e m/s \n',...
  alpha,sigma,Vslip),'HorizontalAlignment','right',...
  'Units','normalized','FontSize',10);
ylabel(ax,'Moment rate (Nm/s)');
xlabel(ax,'Time (s)');




