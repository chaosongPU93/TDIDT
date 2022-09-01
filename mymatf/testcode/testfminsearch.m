% testfminsearch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to test the usage of built-in 'fminsearch' that implements the downhill
% simplex method, and apply this method to our case of directivity effect
% to invert for the rupture speed, direction and length, in order to fit 
% the observed corner frequency
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/11/18
% Last modified date:   2021/11/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear
% close all

Vs = 4.6; % shear wave speed
alpha = [2.4; 22.7; 22.3; 50.6];  % the azimuth wrt center of fam 002 at KLNB, PGC, SSIB, SILB
eta = [143; 143; 131; 140];  % take-off angle from center of fam 002
% fc_obs = [3.5; 3.5; 3.0; 3.8]; % observation data
% fc_obs = [3.6; 3.4; 3.0; 2.9]; % observation data
fc_obs = [5.27; 4.54; 4.0; 4.0]; % observation data, chosen from median spectra
% fc_obs = [3.4277e+00; 3.7797e+00; 4.0988e+00; 4.2099e+00];  % forward prediction of m0 = [0.5, 60, 0.8] 

m0 = [0.8, 28, 1.8];
% m0 = [0.2, 150, 0.2];
% m0 = [0.5, 60, 0.8];

%%% type of objective function,
%   L1norm: sum(abs(X))
%   L2norm: sum(abs(X).^2)^(1/2)
%   Median: median(abs(X))
objtype = 'Median';   

[m,fval,history] = runfminsearch(objtype,m0,Vs,alpha,eta,fc_obs);

fc_pred = m(1).*Vs./m(3)./(1-m(1).*cosd(m(2)-alpha).*sind(eta)); % model-predicted data
% fc_pred_m0 = m0(1).*Vs./m0(3)./(1-m0(1).*cosd(m0(2)-alpha).*sind(eta)); % model-predicted data


%%
f.fig = figure;
f.fig.Renderer = 'painters';
widin = 5;  % maximum width allowed is 8.5 inches
htin = 7;   % maximum height allowed is 11 inches
% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);

nrow = length(m)+1;
ncol = 1;

for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
%     grid(f.ax(isub),'on');
end

pltxran = [0.1 0.95]; pltyran = [0.1 0.95];
pltxsep = 0.05; pltysep = 0.05; 
axpos = optaxposition(nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

for isub = 1:nrow*ncol
  set(f.ax(isub), 'position', axpos(isub,:));
end

ax = f.ax(1);
hold(ax,'on');
ax.Box = 'on';
grid(ax,'on'); 
scatter(ax,0:size(history.fval,1)-1, history.fval,20,[0.8 0.8 0.8],'filled','o',...
            'MarkerEdgeColor','k');
title(ax, [sprintf('Final model: Vr/Vs=%.2f, ',m(1)), '\theta=', ...
           sprintf('%.1f^o, L=%.2fkm',m(2),m(3))]); %, 'interpreter','latex'
text(ax, 0.9, 0.8, sprintf('Final function value: %.2e', fval),'unit','normalized',...
  'horizontalalignment','right','fontsize',12);
xlabel(ax,'Iteration number');
ylabel(ax,'Objective function value');

ax = f.ax(2);
hold(ax,'on');
ax.Box = 'on';
grid(ax,'on'); 
k = history.x(:,1);   % velocity ratio between the rupture speed and shear wave speed
scatter(ax,0:size(k,1)-1, k,20,[0.8 0.8 0.8],'filled','o',...
            'MarkerEdgeColor','k');
xlabel(ax,'Iteration number');
ylabel(ax,'Vr/Vs');

ax = f.ax(3);
hold(ax,'on');
ax.Box = 'on';
grid(ax,'on'); 
theta = history.x(:,2);  % rupture direction
scatter(ax,0:size(theta,1)-1, theta,20,[0.8 0.8 0.8],'filled','o',...
            'MarkerEdgeColor','k');
xlabel(ax,'Iteration number');
ylabel(ax,'Rupture direction \theta (^o)');

ax = f.ax(4);
hold(ax,'on');
ax.Box = 'on';
grid(ax,'on'); 
L = history.x(:,3);  % rupture length
scatter(ax,0:size(L,1)-1, L,20,[0.8 0.8 0.8],'filled','o',...
            'MarkerEdgeColor','k');
xlabel(ax,'Iteration number');
ylabel(ax,'Rupture length L (km)');

%%
function [m,fval,history] = runfminsearch(objtype,m0,Vs,alpha,eta,fc_obs)
% clear
% clc
% close all
history.x = [];
history.fval = [];


% different form of objective function
if strcmp(objtype, 'L1norm')
  fc_misfit = @(m,Vs,alpha,eta,fc_obs)...
    norm(m(1).*Vs./m(3)./(1-m(1).*cosd(m(2)-alpha).*sind(eta))-fc_obs, 1);  % minimize L-1 norm
elseif strcmp(objtype, 'L2norm')
  fc_misfit = @(m,Vs,alpha,eta,fc_obs)...
    norm(m(1).*Vs./m(3)./(1-m(1).*cosd(m(2)-alpha).*sind(eta))-fc_obs); % minimize L-2 norm
elseif strcmp(objtype, 'Median')
  fc_misfit = @(m,Vs,alpha,eta,fc_obs)...
    median(abs(m(1).*Vs./m(3)./(1-m(1).*cosd(m(2)-alpha).*sind(eta))-fc_obs)); % variant that i call 'Median' 
end

% Cs = 4.6; % shear wave speed
% alpha = [2.4; 22.7; 22.3; 50.6];  % the azimuth wrt center of fam 002 at KLNB, PGC, SSIB, SILB
% eta = [143; 143; 131; 140];  % take-off angle from center of fam 002
% fc_obs = [3.5; 3.5; 3.0; 3.8]; % observation data
fobj = @(m)fc_misfit(m,Vs,alpha,eta,fc_obs);

% m0 = [0.8, 28, 1.8];
options = optimset('Display','iter','PlotFcns',@optimplotfval,'MaxFunEvals',2000,...
  'MaxIter',1000,'OutputFcn',@outfun); %

[m,fval] = fminsearch(fobj,m0,options);

k = m(1);   % velocity ratio between the rupture speed and shear wave speed
theta = m(2);  % rupture direction
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
