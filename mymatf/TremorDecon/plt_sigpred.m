function [f] = plt_sigpred(sigsta,predgrp,impindepst,sps,tmaxi,tmaxo,tbosti,tbosto,...
  ircc,rcc,irccp,rccp,mwlen,tstbuf,tedbuf,dy,mo,yr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_sigpred(sigsta,predgrp,impindepst,sps,tstbuf,tedbuf,tmaxi,tmaxo,tbosti,tbosto,...
%   ircc,rcc,dy,mo,yr)
%
% Plot the comparision between the signal waveform fed into deconvoution
% and the prediction computed by the convolution of templates (Green's 
% function) and the preserved impulse sources that are grouped, and the 
% secondary sources are removed afterwards.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/06/22
% Last modified date:   2022/06/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f.fig = figure;
f.fig.Renderer = 'painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 4;   % maximum height allowed is 11 inches
[scrsz, resol] = pixelperinch(1);
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);
nrow = 2;
ncol = 1;
for isub = 1:nrow*ncol
  f.ax(isub) = subplot(nrow,ncol,isub);
end
pltxran = [0.1 0.9]; pltyran = [0.15 0.9];
pltxsep = 0.02; pltysep = 0.05;
%get the locations for each axis
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

lsig = size(sigsta,1);

ax=f.ax(1);
yyaxis(ax,'left');
plot(ax,(1:lsig)/sps, sigsta(:,1), 'r-'); xlim(ax,[0,lsig/sps]); hold(ax,'on');
plot(ax,(1:lsig)/sps, sigsta(:,2), 'b-'); plot(ax,(1:lsig)/sps, sigsta(:,3), 'k-');
ym = max(abs(sigsta(:)));
yran=1.2*[-ym ym];
ylim(ax,yran); ax.Box='on'; grid(ax,'on');
%plot the strongest 0.5-s arrival outside ellipse
ind = find(tmaxo>=tstbuf & tmaxo<=tedbuf);
yloc = (yran(1)+range(yran)*0.05);
for ii = 1: length(ind)
  barst = tmaxo(ind(ii))-tstbuf; % arrival of strongest 0.5-s, in sec
  bared = tmaxo(ind(ii))-tstbuf+0.5; % last for 0.5 s
  plot(ax,[barst bared],[yloc yloc],'-','linew',2.5,'color',[.6 .6 .6]);
end
%plot the strongest 0.5-s arrival inside ellipse
ind = find(tmaxi>=tstbuf & tmaxi<=tedbuf);
for ii = 1: length(ind)
  barst = tmaxi(ind(ii))-tstbuf; % arrival of strongest 0.5-s, in sec
  bared = tmaxi(ind(ii))-tstbuf+0.5; % last for 0.5 s
  plot(ax,[barst bared],[yloc yloc],'k-','linew',2.5);
end
%plot the bostock's LFE catalog outside rectangle
ind = find(tbosto>=tstbuf & tbosto<=tedbuf);
scatter(ax,tbosto(ind)-tstbuf, yloc*ones(size(tbosto(ind))),10,'m','linew',1); % bostock
%plot the bostock's LFE catalog inside rectangle
ind = find(tbosti>=tstbuf & tbosti<=tedbuf);
scatter(ax,tbosti(ind)-tstbuf, yloc*ones(size(tbosti(ind))),10,'m','filled');
yyaxis(ax,'right');
plot(ax,ircc/sps,rcc,'o','color',[.6 .6 .6],'markersize',2);
ylim(ax,[-1 1]);
text(ax,0.98,0.9,'Signal','Units','normalized','HorizontalAlignment','right','FontSize',10);
longticks(ax,2); nolabels(ax,1); hold(ax,'off');

ax=f.ax(2);
yyaxis(ax,'left');
plot(ax,(1:lsig)/sps, predgrp(:,1), 'r-'); xlim(ax,[0,lsig/sps]); hold(ax,'on');
plot(ax,(1:lsig)/sps, predgrp(:,2), 'b-'); plot(ax,(1:lsig)/sps, predgrp(:,3), 'k-');
ylim(ax,yran); ax.Box='on'; grid(ax,'on');
scatter(ax,impindepst(:,1)/sps,yloc*ones(size(impindepst(:,1))),8,'r','filled');
ylabel(ax,'Amplitude','FontSize',12);
yyaxis(ax,'right');
% [irccp,rccp12] = RunningCC(predgrp(:,1), predgrp(:,2), mwlen);
% [~,rccp13] = RunningCC(predgrp(:,1), predgrp(:,3), mwlen);
% [~,rccp23] = RunningCC(predgrp(:,2), predgrp(:,3), mwlen);
% rccp = (rccp12+rccp13+rccp23)/3;
plot(ax,irccp/sps,rccp,'o','color',[.6 .6 .6],'markersize',2);
ylim(ax,[-1 1]);
text(ax,0.98,0.9,'Prediction','Units','normalized','HorizontalAlignment','right','FontSize',10);
%       xlabel(sprintf('Samples at %d Hz',round(1/dt)),'FontSize',12);
xlabel(sprintf('Time (s) since %.1f s on %s %s %s',tstbuf,dy,mo,yr),'FontSize',12);
ylabel('Running CC','FontSize',12);
longticks(ax,2); hold(ax,'off');
