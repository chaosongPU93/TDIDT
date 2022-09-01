% testfminsearchbnd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to test the usage of 'fminsearchbnd' that implements the downhill
% simplex method with the constraint of boundary to the model params,
% and apply this method to our case of directivity effect
% to invert for the rupture speed, direction and length, in order to fit 
% the observed corner frequency.
%
% This function is written by  John D'Errico (2021). 
% fminsearchbnd, fminsearchcon (https://www.mathworks.com/matlabcentral/
%   fileexchange/8277-fminsearchbnd-fminsearchcon), MATLAB Central File 
%   Exchange. Retrieved November 18, 2021.  
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/11/18
% Last modified date:   2021/11/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear
close all

Cs = 4.6; % shear wave speed
azi = [2.4; 22.7; 22.3; 50.6];  % the azimuth wrt center of fam 002 at KLNB, PGC, SSIB, SILB
eta = [143; 143; 131; 140];  % take-off angle from center of fam 002
fc_obs = [3.5; 3.5; 3.0; 3.8]; % observation data
m0 = [0.8, 28, 1.8];
objtype = 'L1norm';

mlbnd = [0.001, 0, 0.001];  % lower boundary of model params
mubnd = [1.5, 360, 10];   % upper boundary of model params
[m,fval,history] = runfminsearchbnd(objtype,m0,mlbnd,mubnd,Cs,azi,eta,fc_obs);

fc_pred = m(1).*Cs./m(3)./(1-m(1).*cosd(m(2)-azi).*sind(eta)); % model-predicted data


%%
figure
subplot(411)
ax = gca;
hold(ax,'on');
ax.Box = 'on';
grid(ax,'on'); 
scatter(ax,0:size(history.fval,1)-1, history.fval,30,[0.7 0.7 0.7],'filled','o',...
            'MarkerEdgeColor','k');
title(ax, sprintf('Final model: Vr/Cs=%.2f, alpha=%.1f^o, L=%.2fkm',m(1),m(2),m(3)));
text(ax, 0.9, 0.8, sprintf('Final function value: %.5f', fval),'unit','normalized',...
  'horizontalalignment','right','fontsize',12);
xlabel(ax,'Iteration number');
ylabel(ax,'Objective function value');

subplot(412)
ax = gca;
hold(ax,'on');
ax.Box = 'on';
grid(ax,'on'); 
k = history.x(:,1);   % velocity ratio between the rupture speed and shear wave speed
scatter(ax,0:size(k,1)-1, k,30,[0.7 0.7 0.7],'filled','o',...
            'MarkerEdgeColor','k');
xlabel(ax,'Iteration number');
ylabel(ax,'Vr/Cs');

subplot(413)
ax = gca;
hold(ax,'on');
ax.Box = 'on';
grid(ax,'on'); 
alpha = history.x(:,2);  % rupture direction
scatter(ax,0:size(alpha,1)-1, alpha,30,[0.7 0.7 0.7],'filled','o',...
            'MarkerEdgeColor','k');
xlabel(ax,'Iteration number');
ylabel(ax,'Rupture direction alpha (^o)');

subplot(414)
ax = gca;
hold(ax,'on');
ax.Box = 'on';
grid(ax,'on'); 
L = history.x(:,3);  % rupture length
scatter(ax,0:size(L,1)-1, L,30,[0.7 0.7 0.7],'filled','o',...
            'MarkerEdgeColor','k');
xlabel(ax,'Iteration number');
ylabel(ax,'Rupture length (km)');

%%
function [m,fval,history] = runfminsearchbnd(objtype,m0,mlbnd,mubnd,Cs,azi,eta,fc_obs)
% clear
% clc
% close all
history.x = [];
history.fval = [];


% different form of objective function
if strcmp(objtype, 'L1norm')
  fc_misfit = @(m,Cs,azi,eta,fc_obs)...
    norm(m(1).*Cs./m(3)./(1-m(1).*cosd(m(2)-azi).*sind(eta))-fc_obs, 1);  % minimize L-1 norm
elseif strcmp(objtype, 'L2norm')
  fc_misfit = @(m,Cs,azi,eta,fc_obs)...
    norm(m(1).*Cs./m(3)./(1-m(1).*cosd(m(2)-azi).*sind(eta))-fc_obs); % minimize L-2 norm
end

% Cs = 4.6; % shear wave speed
% azi = [2.4; 22.7; 22.3; 50.6];  % the azimuth wrt center of fam 002 at KLNB, PGC, SSIB, SILB
% eta = [143; 143; 131; 140];  % take-off angle from center of fam 002
% fc_obs = [3.5; 3.5; 3.0; 3.8]; % observation data
fobj = @(m)fc_misfit(m,Cs,azi,eta,fc_obs);

% m0 = [0.8, 28, 1.8];
options = optimset('Display','iter','PlotFcns',@optimplotfval,'MaxFunEvals',1000,...
  'OutputFcn',@outfun); 
options.TolFun = 1.e-7;
options.TolX = 1.e-7;
[m,fval] = fminsearchbnd(fobj,m0,mlbnd,mubnd,options);

k = m(1);   % velocity ratio between the rupture speed and shear wave speed
alpha = m(2);  % rupture direction
L = m(3);  % rupture length


  function stop = outfun(x,optimValues,state)
    stop = false;
    
    switch state
      case 'init'
%         hold on
      case 'iter'
        % Concatenate current point and objective function
        % value with history. x must be a row vector.
        history.fval = [history.fval; optimValues.fval];
        history.x = [history.x; x];
      case 'done'
%         hold off
      otherwise
        
    end
  end

end
