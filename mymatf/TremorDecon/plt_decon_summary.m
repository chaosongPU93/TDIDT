function [f] = plt_decon_summary(sigsta,resgrp,pkindepsave,impindepst,sps,tmaxi,tmaxo,tbosti,tbosto,...
  ircc,rcc,mwlen,off1i,loff_max,tstbuf,tedbuf,dy,mo,yr,k,ftrans,xcut,ycut)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_decon_summary(sigsta,resgrp,pkindep,impindep,sps,tmaxi,tmaxo,tbosti,tbosto,...
%   ircc,rcc,off1i,loff_max,tstbuf,tedbuf,dy,mo,yr)
%
% This function would create a summary figure for each data win, so that we 
% will have a feeling for all migrations.
%
% --It contains the comparison between
%   the deconvolved positive peaks and signal waveform peaks at each station;
%   the signal waveform with the 4-s tremor detections and M. Bostock's LFE
%   detections indicated; residual waveform with the arrival time of deconvolved
%   impulses at sta 1 indicated; scatter of sources in terms of offsets, 
%   accounting for prealignment offset; the scatter of sources in terms of rela 
%   locations.
%
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/06/29
% Last modified date:   2022/06/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%compose a summary figure for each data win and save it, so that we will have a feeling for
%%%all migrations
f.fig = figure;
f.fig.Renderer = 'painters';
widin = 8.5;  % maximum width allowed is 8.5 inches
htin = 10;   % maximum height allowed is 11 inches
% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, resol] = pixelperinch(1);
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);
nrow = 7;
ncol = 1;
for isub = 1:nrow*ncol
  f.ax(isub) = subplot(nrow,ncol,isub);
end

%get the locations for each axis
axpos = [0.08 0.87 0.85 0.1;
         0.08 0.755 0.85 0.1;
         0.08 0.64 0.85 0.1;
         0.08 0.525 0.85 0.1;
         0.08 0.41 0.85 0.1;
         0.08 0.02 0.4 0.36;
         0.56 0.02 0.4 0.36;
         ];
for isub = 1:nrow*ncol
  set(f.ax(isub), 'position', axpos(isub,:));
end

color = ['r';'b';'k'];

ym = max(abs(sigsta(:)));
yran=1.2*[-ym ym];

lsig = size(sigsta,1); 
nsta = size(sigsta,2); 

for i = 1: nsta
  ax=f.ax(i);
  hold(ax,'on');
  plot(ax, (1:lsig)/sps, sigsta(:,i), '-','Color',color(i,:)); 
  xlim(ax,[0,lsig/sps]);  ylim(ax,yran);
  [pkhgt, pk] = findpeaks(sigsta(:,i));
  scatter(ax,pk/sps,pkhgt,8,color(i,:));
  for ii = 1: size(pkindepsave,1)
    plot(ax,[pkindepsave(ii,1) pkindepsave(ii,1)]/sps, ax.YLim, '--','Color',color(i,:));
  end
  text(ax,0.98,0.9,sprintf('Signal %d',i),'Units','normalized','HorizontalAlignment','right',...
    'FontSize',10);
  text(ax,0.02,0.9,num2str(size(pkindepsave,1)),'fontsize',10,...
    'Units','normalized');
  ax.Box='on'; grid(ax,'on'); 
  longticks(ax,3); nolabels(ax,1);  
%   if i==1
%     h=supertit(ax, sprintf('Window %s, %d, %.1f-%.1f s, Passband: %.1f-%.1f Hz, sps: %d',...
%       num2zeropadstr(k, 3),date,tstbuf,tedbuf,losig,hisig,sps),10);
%     movev(h,-0.05);
%   end
  hold(ax,'off');
end

ax=f.ax(4);
hold(ax,'on');
yyaxis(ax,'right');
% plot(ax,ircc/sps,rcc,'o','color',[.7 .7 .7],'markersize',1);
plot(ax,ircc/sps,rcc,'color',[.7 .7 .7]);
ylim(ax,[-1 1]);
yyaxis(ax,'left');
for i = 1: nsta
  plot(ax,(1:lsig)/sps, sigsta(:,i), '-','Color',color(i,:)); 
end
xlim(ax,[0,lsig/sps]);  ylim(ax,yran); 
ax.Box='on'; grid(ax,'on');
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
text(ax,0.98,0.9,'Signal','Units','normalized','HorizontalAlignment','right','FontSize',10);
longticks(ax,3); nolabels(ax,1); hold(ax,'off');

ax=f.ax(5);
hold(ax,'on');
yyaxis(ax,'right');
[irccr,rccr12] = RunningCC(resgrp(:,1), resgrp(:,2), mwlen);
[~,rccr13] = RunningCC(resgrp(:,1), resgrp(:,3), mwlen);
[~,rccr23] = RunningCC(resgrp(:,2), resgrp(:,3), mwlen);
rccr = (rccr12+rccr13+rccr23)/3;
% plot(ax,irccr/sps,rccr,'o','color',[.7 .7 .7],'markersize',1);
plot(ax,irccr/sps,rccr,'color',[.7 .7 .7]);
ylim(ax,[-1 1]);
ylabel(ax,'Running CC','FontSize',11);
yyaxis(ax,'left');
for i = 1: nsta
  plot(ax,(1:lsig)/sps, resgrp(:,i), '-','Color',color(i,:)); 
end
xlim(ax,[0,lsig/sps]);  ylim(ax,yran); 
ax.Box='on'; grid(ax,'on');
scatter(ax,impindepst(:,1)/sps,yloc*ones(size(impindepst(:,1))),8,'r','filled');
text(ax,0.02,0.9,num2str(size(impindepst,1)),'fontsize',10,...
  'Units','normalized');
ylabel(ax,'Amplitude','FontSize',11);
text(ax,0.98,0.9,'Residual','Units','normalized','HorizontalAlignment','right','FontSize',10);
%       xlabel(sprintf('Samples at %d Hz',round(1/dt)),'FontSize',12);
xlabel(ax,sprintf('Time (s) since %.1f s on %s %s %s',tstbuf,dy,mo,yr),'FontSize',11);
longticks(ax,3); hold(ax,'off');

%%%plot the scatter of offsets, accounting for prealignment offset, == true offset
ax=f.ax(6);
xran = [-loff_max+off1i(k,2)-1 loff_max+off1i(k,2)+1];
yran = [-loff_max+off1i(k,3)-1 loff_max+off1i(k,3)+1];
offxran = [-loff_max+off1i(k,2) loff_max+off1i(k,2)];
offyran = [-loff_max+off1i(k,3) loff_max+off1i(k,3)];
cran = [0 lsig];
[ax] = plt_decon_imp_scatter(ax,impindepst,xran,yran,cran,offxran,offyran,...
  sps,50,'mean','tori');
scatter(ax,off1i(k,2),off1i(k,3),20,'ks','filled','MarkerEdgeColor','k');
hold(ax,'off');

ax=f.ax(7);
xran = [-3 3];
yran = [-3 3];
cran = [0 lsig/sps];
[ax] = plt_decon_imp_scatter_space(ax,impindepst,xran,yran,cran,offxran,offyran,...
  sps,50,ftrans,'mean','tori');
plot(ax,xcut,ycut,'k-','linew',2);
hold(ax,'off');

