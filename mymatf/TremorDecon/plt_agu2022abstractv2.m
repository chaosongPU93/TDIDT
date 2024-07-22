function [f] = plt_agu2022abstractv2(greenf,sigsta,impindepst,sps,xzoom,off1iw,loff_max,tstbuf,dy,mo,yr,ftrans)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_agu2022abstractv2(sigsta,resgrp,pkindepsave,impindepst,sps,tmaxi,tmaxo,...
%   tbosti,tbosto,ircccat,rcccat,off1iw,off1i,loff_max,tstbuf,tedbuf,dy,mo,yr)
%
% This function would create a simplified version of the plot for the 2022
% AGU abstract. It will contain template waveforms; signal waveforms and zoom-in
% at a short period so that you can see the wax and wane, maybe no running
% cc, with dots indicating the decon sources; map-view of sources distribution
% only of that short window for demonstration, no boundaries.
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
widin = 8.4;  % maximum width allowed is 8.5 inches
htin = 10;   % maximum height allowed is 11 inches
% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, resol] = pixelperinch(1);
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);
nrow = 5;
ncol = 1;
for isub = 1:nrow*ncol
  f.ax(isub) = subplot(nrow,ncol,isub);
end

%get the locations for each axis
axpos = [0.08 0.87 0.85 0.1;
         0.08 0.74 0.85 0.1;
         0.08 0.59 0.85 0.1;
         0.08 0.215 0.38 0.38;
         0.54 0.215 0.38 0.38;
         ];
for isub = 1:nrow*ncol
  set(f.ax(isub), 'position', axpos(isub,:));
end

color = ['r';'b';'k'];

ym = max(abs(sigsta(:)));
yran=1.5*[-ym ym];

lsig = size(sigsta,1); 
nsta = size(sigsta,2); 
lwlet = size(greenf,1);

%%%seismograms of signal
ax=f.ax(1);
hold(ax,'on');
% yyaxis(ax,'right');
% % plot(ax,ircccat/sps,rcccat,'o','color',[.7 .7 .7],'markersize',1);
% plot(ax,ircccat/sps,rcccat,'color',[.7 .7 .7]);
% ylim(ax,[-1 1]);
% yyaxis(ax,'left');
for i = 1: nsta
  plot(ax,(1:lsig)/sps, sigsta(:,i), '-','Color',color(i,:));
end
xlim(ax,[0,lsig/sps]);  ylim(ax,yran); 
ax.Box='on'; grid(ax,'on');
%plot the deconvolved sources
yloc = (yran(1)+range(yran)*0.05);
scatter(ax,impindepst(:,1)/sps,yloc*ones(size(impindepst(:,1))),6,[.7 .7 .7],'filled');
% text(ax,0.5,0.9,num2str(size(impindepst,1)),'fontsize',10,...
%   'Units','normalized');
%plot the deconvolved sources
yloc = (yran(1)+range(yran)*0.05);
% xzoom = [20 35];
% xzoom = [185 200];
% xzoom = [225 245];
ind = find(impindepst(:,1)/sps >= xzoom(1) & impindepst(:,1)/sps <= xzoom(2));
scatter(ax,impindepst(ind,1)/sps,yloc*ones(size(impindepst(ind,1))),8,[.3 .3 .3],'filled');
plot(ax,[xzoom(1) xzoom(1)],yran,'--','color',[.3 .3 .3],'linew',1);
plot(ax,[xzoom(2) xzoom(2)],yran,'--','color',[.3 .3 .3],'linew',1);
text(ax,0.98,0.9,'Signal','Units','normalized','HorizontalAlignment','right','FontSize',10);
% text(ax,0.01,0.9,'a','FontSize',11,'unit','normalized');
% legend(ax,p,'PGC','SSIB','SILB','Location','north','NumColumns',3,'fontsize',8);
longticks(ax,3); 
hold(ax,'off');

%%%plot the scatter of sources in terms of rela locations
ax=f.ax(5);
xran = [-4 4];
yran = [-4 4];
% cran = [0 lsig/sps];
cran = xzoom;
[ax,tori] = plt_decon_imp_scatter_space_ref(ax,impindepst(ind,:),xran,yran,cran,off1iw,loff_max,...
  sps,35,ftrans,'mean','tori','none');
% text(ax,0.02,0.96,'d','FontSize',11,'unit','normalized','HorizontalAlignment','left');
% pos = ax.Position;
% c.Position = [pos(1)+pos(3)+0.085, pos(2), 0.02, pos(3)];
% plot(ax,xcut,ycut,'k-','linew',2);
xticks(ax,xran(1): 1 : xran(2));
yticks(ax,yran(1): 1 : yran(2));
hold(ax,'off');

%%%plot the scatter of sources in terms of rela locations
ax=f.ax(4);
cran = [0 lsig/sps];
ax = plt_decon_imp_scatter_space_ref(ax,impindepst,xran,yran,cran,off1iw,loff_max,...
  sps,35,ftrans,'mean','tori','none');
% text(ax,0.02,0.96,'e','FontSize',11,'unit','normalized','HorizontalAlignment','left');
% pos = ax.Position;
% c.Position = [pos(1)+pos(3)+0.085, pos(2), 0.02, pos(3)];
% plot(ax,xcut,ycut,'k-','linew',2);
xticks(ax,xran(1): 1 : xran(2));
yticks(ax,yran(1): 1 : yran(2));
angrmse = 140;
[rotx, roty] = complex_rot(0,1,-angrmse);
xvect = [0.5-rotx 0.5+rotx];
yvect = [-2.5-roty -2.5+roty];
drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1);
% text(ax,0.63,0.2,strcat(num2str(angrmse),'^{\circ}'),'FontSize',9,...
%   'unit','normalized','horizontalalignment','left');
text(ax,0.62,0.18,strcat(num2str(angrmse),'$^{\circ}$'),'FontSize',9,...
  'unit','normalized','interpreter','latex');
hold(ax,'off');

%%%seismograms of signal, zoom-in
ax=f.ax(2);
hold(ax,'on');
%%%obtain a single best alignment based on the zoom-in segment
msftadd = 20; 
optcc = sigsta(xzoom(1)*sps+1: xzoom(2)*sps, :);
ccmid = ceil(size(optcc,1)/2);
ccwlen = round(size(optcc,1)-2*(msftadd+1));
loffmax = 4*sps/40;
ccmin = 0.01;  % depending on the length of trace, cc could be very low
iup = 1;    % times of upsampling
[off12con,off13con,ccali,iloopoff,loopoff] = constrained_cc_interp(optcc',ccmid,...
  ccwlen,msftadd,loffmax,ccmin,iup);
% if a better alignment cannot be achieved, use 0,0
if off12con == msftadd+1 && off13con == msftadd+1
  off12con = 0;
  off13con = 0;
  fprintf('This segment cannot be properly aligned, double-check needed \n');
end
sigseg(:, 1) = sigsta(xzoom(1)*sps+1: xzoom(2)*sps, 1);
sigseg(:, 2) = sigsta(xzoom(1)*sps+1-off12con: xzoom(2)*sps-off12con, 2); % sta 2
sigseg(:, 3) = sigsta(xzoom(1)*sps+1-off13con: xzoom(2)*sps-off13con, 3); % sta 3

mwlen = sps/2;
[ircc,rcc12] = RunningCC(sigseg(:,1), sigseg(:,2), mwlen);
[~,rcc13] = RunningCC(sigseg(:,1), sigseg(:,3), mwlen);
[~,rcc23] = RunningCC(sigseg(:,2), sigseg(:,3), mwlen);
% ircc = ircc+*sps;
rcc = (rcc12+rcc13+rcc23)/3;

lseg = size(sigseg,1);
% yyaxis(ax,'left');
for i = 1: nsta
  p(i)=plot(ax,(1:lseg)/sps, sigseg(:,i), '-','Color',color(i,:),'linew',1); 
end
xlim(ax,[0 lseg/sps]);  
xticks(ax,0: 2: lseg/sps);
tkval = xzoom(1): 2: xzoom(2);
tklbl = cell(length(tkval),1);
for jj = 1: length(tkval)
  tklbl{jj} = num2str(tkval(jj));
end
xticklabels(ax,tklbl);
yran=1.5*[-ym ym];
ylim(ax,yran); 
% yyaxis(ax,'left');
plot(ax,ircc/sps,yran(2)*rcc,'o','markersize',1,'Color',[.5 .5 .5]); % scale with data amp.
ax.Box='on'; grid(ax,'on');
%plot the deconvolved sources
yloc = (yran(1)+range(yran)*0.05);
ind = find(impindepst(:,1)/sps >= xzoom(1) & impindepst(:,1)/sps <= xzoom(2));
% scatter(ax,impindepst(ind,1)/sps-xzoom(1),yloc*ones(size(impindepst(ind,1))),8,'r','filled');
scatter(ax,impindepst(ind,1)/sps-xzoom(1),yloc*ones(size(impindepst(ind,1))),10,tori,'filled',...
  'MarkerEdgeColor',[.5 .5 .5]);
oldc = colormap(ax,'kelicol');
newc = flipud(oldc);
colormap(ax,newc);
caxis(ax,xzoom);
text(ax,0.98,0.9,'Signal zoom-in','Units','normalized','HorizontalAlignment','right','FontSize',10);
% text(ax,0.01,0.9,'b','FontSize',11,'unit','normalized');
ylabel(ax,'Amplitude','FontSize',10);
xlabel(ax,sprintf('Time (s) since %.1f s on %s %s %s',tstbuf,dy,mo,yr),'FontSize',10);
% legend(ax,p,'PGC','SSIB','SILB','Location','north','NumColumns',3,'fontsize',8);
longticks(ax,3); 
hold(ax,'off');

% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(2),0,yran(2));
[xe,ye] = ds2nfu(f.ax(1),xzoom(1),yran(1));
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color',[.3 .3 .3],'linestyle','-');
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(2),lseg/sps,yran(2));
[xe,ye] = ds2nfu(f.ax(1),xzoom(2),yran(1));
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color',[.3 .3 .3],'linestyle','-');

%%%seismograms of templates
ax=f.ax(3);
hold(ax,'on');
for i = 1: nsta
  p(i)=plot(ax,(1:lwlet)/sps, greenf(:,i), '-','Color',color(i,:),'linew',1); 
end
ax.Box='on'; grid(ax,'on');
xlim(ax,[0,lwlet/sps]);
% ylim(ax,yran/2.5); 
ylim(ax,[-0.5 0.5]); 
xticks(ax,0: 2: lwlet/sps);
text(ax,0.98,0.9,'Template','Units','normalized','HorizontalAlignment','right','FontSize',10);
% text(ax,0.01,0.9,'c','FontSize',11,'unit','normalized');
ylabel(ax,'Amplitude','FontSize',10);
xlabel(ax,'Time (s)','FontSize',10);
legend(ax,p,'PGC','SSIB','SILB','Location','southeast','NumColumns',3,'fontsize',8);
longticks(ax,3); 
if lwlet <= lseg
  loc = getpos(f.ax(3));
  shrink(f.ax(3),lseg/lwlet,1);
  f.ax(3).Position(1)=loc(1);
else
  loc = getpos(f.ax(2));
  shrink(f.ax(2),lwlet/lseg,1);
  f.ax(2).Position(1)=loc(1);
end











