function f=plt_seiszoomedin(sigsta,xzoom,sps,tstbuf,tbosto,...
  tbosti,timp,timp4th,ircc,rcc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_seiszoomedin(sigseg,xzoom,sps,tstbuf,tbosto,...
%   tbosti,timp,timp4th,ircc,rcc)
%
% Function just to ease the plotting of the seismograms at trio stations
% plus a 4th station, along with the deconvolved LFE catalogs, Bostock's
% detections, rcc, etc.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/05/01
% Last modified date:   2023/05/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%obtain a single best alignment based on the zoom-in segment
msftadd = 20;
optcc = sigsta(xzoom(1)*sps+1: xzoom(2)*sps, :);
ccmid = ceil(size(optcc,1)/2);
ccwlen = round(size(optcc,1)-2*(msftadd+1));
loffmax = 4*sps/40;
ccmin = 0.01;  % depending on the length of trace, cc could be very low
iup = 1;    % times of upsampling
[off12con,off13con,ccali] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
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

[mcoef,off17] = xcorrmax(optcc(:,1), optcc(:,7), msftadd, 'coeff');
sigseg(:, 4) = sigsta(xzoom(1)*sps+1-off17: xzoom(2)*sps-off17, nsta); % sta 3

[ircc,rcc12] = RunningCC(sigseg(:,1), sigseg(:,2), rccmwlen);
[~,rcc13] = RunningCC(sigseg(:,1), sigseg(:,3), rccmwlen);
[~,rcc23] = RunningCC(sigseg(:,2), sigseg(:,3), rccmwlen);
% ircc = ircc+*sps;
rcc = (rcc12+rcc13+rcc23)/3;

widin = 15;  % maximum width allowed is 8.5 inches
htin = 3;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 1;
f = initfig(widin,htin,nrow,ncol); %initialize fig

xran = [0.08 0.92]; yran = [0.15 0.90];
xsep = 0.02; ysep = 0.04;
optaxpos(f,nrow,ncol,xran,yran,xsep,ysep);

color = ['r';'b';'k';'c'];
stas = ['PGC ';'SSIB';'SILB';'KLNB'];

nsta = size(sigseg,2);
lseg = size(sigseg,1);

ym = max(abs(sigseg(:)));
yran=1.2*[-ym ym];

tkval = xzoom(1): 2: xzoom(2);
tklbl = cell(length(tkval),1);
for jj = 1: length(tkval)
  tklbl{jj} = num2str(tkval(jj));
end

ax=f.ax(1);
hold(ax,'on'); ax.Box='on'; grid(ax,'on');
for i = 1: nsta
  p(i)=plot(ax,(1:lseg)/sps, sigseg(:,i), '-','Color',color(i,:),'linew',0.5);
end
plot(ax,ircc/sps,yran(2)*rcc,'o','markersize',1,'Color',[.5 .5 .5]); % scale with data amp.
ind = find(timp/sps >= xzoom(1) & timp/sps <= xzoom(2));
if ~isempty(ind)
  yloc = (yran(1)+range(yran)*0.03);
  scatter(ax,timp(ind)/sps-xzoom(1),...
    yloc*ones(size(timp(ind))),14,[.7 .7 .7],'filled');
end
ind = find(timp4th/sps >= xzoom(1) & timp4th/sps <= xzoom(2));
if ~isempty(ind)
  yloc = (yran(1)+range(yran)*0.07);
  scatter(ax,timp4th(ind)/sps-xzoom(1),...
    yloc*ones(size(timp4th(ind))),14,[.3 .3 .3],'filled');
end
%plot the bostock's LFE catalog inside ellipse
indbo = find(tbosto-tstbuf >= xzoom(1) & tbosto-tstbuf <= xzoom(2) );
if ~isempty(indbo)
  yloc = (yran(1)+range(yran)*0.11);
  scatter(ax,tbosto(indbo)-tstbuf-xzoom(1), yloc*ones(size(tbosto(indbo))),14,'g','linew',1); % bostock
end
%plot the bostock's LFE catalog inside rectangle
indbi = find(tbosti-tstbuf >= xzoom(1) & tbosti-tstbuf <= xzoom(2));
if ~isempty(indbi)
  yloc = (yran(1)+range(yran)*0.11);
  scatter(ax,tbosti(indbi)-tstbuf-xzoom(1), yloc*ones(size(tbosti(indbi))),14,'g','filled');
end
xlim(ax,[0 lseg/sps]);
xticks(ax,0: 2: lseg/sps);
xticklabels(ax,tklbl);
ylim(ax,yran);
text(ax,0.98,0.9,'Signal zoom-in','Units','normalized','HorizontalAlignment','right','FontSize',10);
ylabel(ax,'Amplitude','FontSize',10);
xlabel(ax,sprintf('Time (s) since %.1f s',tstbuf),'FontSize',10);
legend(ax,p,stas(1:nsta,:),'Location','north','NumColumns',nsta,'fontsize',8);



