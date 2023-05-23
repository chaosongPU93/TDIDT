function [f,predgrp,resgrp,predgrpl,resgrpl,l2normred]=predsig_conv_imptemp(sigsta,optdat,impindepst,...
  greenf,zcrosses,overshoot,stas,flagplt,sps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script obtains the predicted (modeled) waveform by carring out a 
% convolution between the arrival time of deconvolved impulses (usually 
% the 2ndary sources have been removed) and template at the corresponding
% stations. It is capable of dealing with both the optimal and orthogonal
% components, as the tarvl is the same at both comp., so you just change
% the 'sigsta' and 'greenf' and other related inputs to the comp. you want.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/10/12
% Last modified date:   2022/10/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('flagplt',1)
defval('sps',160);

[lsig,nsta] = size(sigsta);
predgrp = zeros(lsig,nsta);  % predefine the prediction array
impfull = zeros(lsig,nsta);  % a full array of impulses containing zeros
predgrpl = zeros(lsig+overshoot*2, nsta);   % longer including an overshoot on both sides for alignment
impfulll = zeros(lsig+overshoot*2, nsta);
mcoef = zeros(nsta,1);
mlag = zeros(nsta,1);
for ista = 1: nsta
  if ista <= 3
    impfull(impindepst(:,(ista-1)*2+1), ista) = impindepst(:,ista*2); % the non-zero index has the amp
  else
    ind = find(impindepst(:,9+(ista-4)*2+1)>0); %in case there is any source' arrival is less than 0
    impfull(impindepst(ind,9+(ista-4)*2+1), ista) = impindepst(ind,9+(ista-3)*2);
  end
  predtmp = conv(greenf(:,ista), impfull(:,ista), 'full');
  itwlet = zcrosses(ista);
  predgrp(:,ista) = predtmp(itwlet:lsig+itwlet-1);  % cut accordingly
  
  %obtain a longer predication and residual for alignment in each subwin
  if ista <= 3
    impfulll(impindepst(:,(ista-1)*2+1)+overshoot, ista) = impindepst(:,ista*2);
  else
    impfulll(impindepst(:,9+(ista-4)*2+1)+overshoot, ista) = impindepst(:,9+(ista-3)*2);
  end
  predtmpl = conv(greenf(:,ista), impfulll(:,ista), 'full');
  predgrpl(:,ista) = predtmpl(itwlet:lsig+overshoot*2+itwlet-1);  % cut accordingly
  
  %obtain a CC to see if pred and sig are also well aligned
  msftadd = 20;
  [coef,lag] = xcorr(sigsta(:,ista), predgrp(:,ista), msftadd, 'coeff');
  [mcoef(ista), idx] = max(coef);   % max of master raw cc
  mlag(ista) = lag(idx);
  
  %shift to the best alignment
  predgrpa(:,ista) = predgrpl(1+overshoot-mlag(ista): end-overshoot-mlag(ista),ista);
  
end

resgrp = sigsta - predgrp;  %residual, difference between signal and prediction
resgrpa = sigsta - predgrpa;  %residual, difference between aligned signal and prediction
resgrpl = detrend(optdat(:,2:end)) - predgrpl;  %longer residual including an overshoot on both sides for alignment

l2normred = zeros(nsta, 3);  %residual reduction
for i = 1: nsta
  l2normred(i, :) = [norm(sigsta(:,i)) norm(resgrp(:,i)) ...
    (norm(sigsta(:,i))-norm(resgrp(:,i)))/norm(sigsta(:,i))*100];  
end

if flagplt
%   %%%comparison between the prediction and longer prediction, they'd have a phase shift of overshoot
%   figure;
%   subplot(211)
%   plot(predtmp,'b');hold on;
%   plot(predtmpl,'r');
%   subplot(212)
%   plot(predgrp(:,end),'b');hold on;
%   plot(predgrpa(:,end),'k');hold on;
%   plot(predgrpl(:,end),'r');

  %%%signal vs. residual, the absolute and relative res should be small
  f.fig = figure;
  subplot(3,2,1)
  hold on
  box on; grid on
  plot((1:lsig)/sps,sigsta(:,1),'r');
  plot((1:lsig)/sps,sigsta(:,2),'b');
  plot((1:lsig)/sps,sigsta(:,3),'k');
  ym = max(abs(sigsta(:)));
  yran=1.2*[-ym ym];
  ylim(yran);
  text(0.95,0.9,sprintf('%.2f; %.2f; %.2f',l2normred(1, 1),l2normred(2, 1),l2normred(3, 1)),...
    'Units','normalized','HorizontalAlignment','right');
  text(0.95,0.1,sprintf('%s; %s; %s',strtrim(stas(1,:)),strtrim(stas(2,:)),...
    strtrim(stas(3,:))),'Units','normalized','HorizontalAlignment','right');
  
  subplot(3,2,2)
  hold on
  box on; grid on
  plot((1:lsig)/sps,resgrp(:,1),'r');
  plot((1:lsig)/sps,resgrp(:,2),'b');
  plot((1:lsig)/sps,resgrp(:,3),'k');
  ylim(yran);
  text(0.95,0.9,sprintf('%.2f; %.2f; %.2f',l2normred(1, 2),l2normred(2, 2),l2normred(3, 2)),...
    'Units','normalized','HorizontalAlignment','right');
  text(0.95,0.8,sprintf('%.1f%%; %.1f%%; %.1f%%',l2normred(1, 3),l2normred(2, 3),l2normred(3, 3)),...
    'Units','normalized','HorizontalAlignment','right');
  
  
  if nsta >3
    for ista = 4:nsta
      subplot(3,2,2+ista-3)
      hold on
      box on; grid on
      plot((1:lsig)/sps,sigsta(:,ista),'k');
      plot((1:lsig)/sps,resgrp(:,ista),'r');
      ylim(yran);
      text(0.95,0.9,sprintf('%.2f; %.2f; %.1f%%',l2normred(ista, 1),l2normred(ista, 2),l2normred(ista, 3)),...
        'Units','normalized','HorizontalAlignment','right');
      text(0.95,0.1,sprintf('%s',strtrim(stas(ista,:))),...
        'Units','normalized','HorizontalAlignment','right');

    end
  end
  
else
  f = [];
end

% keyboard








