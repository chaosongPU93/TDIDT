function f=plt_002lfe_signalzoom(greenf,sigsta,impindepst,sps,xzoom,tstbuf,dy,mo,yr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
% This is the function to plot the seismograms of the stacked LFE templates
% 'greenf' that are bandpassed. And a zoom-in window of the signal seismograms
% 'signal', this is supposed to be a figure shown in the paper. 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/02/08
% Last modified date:   2024/02/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%compose a summary figure for each data win and save it, so that we will have a feeling for
%%%all migrations
widin = 9;  % maximum width allowed is 8.5 inches
htin = 3;   % maximum height allowed is 11 inches
nrow = 2;
ncol = 1;
f=initfig(widin,htin,nrow,ncol);

figxran = [0.08 0.94]; figyran = [0.15 0.9];
figxsep = 0.05; figysep = 0.05;
axpos=optaxpos(f,nrow,ncol,figxran,figyran,figxsep,figysep);


color = ['r';'b';'k'];

ym = max(abs(sigsta(:)));
yran=1.2*[-ym ym];

lsig = size(sigsta,1); 
nsta = 3; 
lwlet = size(greenf,1);
% lwlet = 12*sps;

%%%LFE templates
ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on'); 
for i = 1: nsta
  p(i)=plot(ax,(1:lwlet)/sps, greenf(1:lwlet,i), '-','Color',color(i,:),'linew',1); 
end
legend(ax,p,{'PGC','SSIB','SILB'},'Location','southeast','NumColumns',3,...
  'fontsize',8);
text(ax,0.02,0.2,'1.8-18 Hz','Units','normalized','HorizontalAlignment','left',...
  'FontSize',10);
text(ax,0.98,0.9,'LFE templates at family 002','Units','normalized',...
  'HorizontalAlignment','right','FontSize',10);
% text(ax,0.02,0.9,'b','FontSize',11,'unit','normalized');
xlim(ax,[0 lwlet/sps]);  
ylim(ax,yran/2);
shrink(ax,range(xzoom)*sps/lwlet,1);
ax.Position(1)=axpos(1,1);
xticks(ax,0: 2: lwlet/sps);
longticks(ax,1.5); 
nolabels(ax,1);

%%%seismograms of signal, zoom-in
ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on'); 
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
yyaxis(ax,'left');
for i = 1: nsta
  p(i)=plot(ax,(1:lseg)/sps, sigseg(:,i), '-','Color',color(i,:),'linew',1); 
end
xlim(ax,[0 lseg/sps]);  
xticks(ax,0: 2: lseg/sps);
% tkval = xzoom(1): 2: xzoom(2);
% tklbl = cell(length(tkval),1);
% for jj = 1: length(tkval)
%   tklbl{jj} = num2str(tkval(jj));
% end
% xticklabels(ax,tklbl);
ylim(ax,yran); 
text(ax,0.99,0.9,'Tremor velocity signal','Units','normalized',...
  'HorizontalAlignment','right','FontSize',10);
text(ax,0.01,0.2,'1.8-6.3 Hz','Units','normalized','HorizontalAlignment','left',...
  'FontSize',10);
% text(ax,0.01,0.9,'c','FontSize',11,'unit','normalized');
ylabel(ax,'Amplitude','FontSize',10);
xlabel(ax,sprintf('Time (s) since %.1f s on %s %s %s',tstbuf,dy,mo,yr),'FontSize',10);
longticks(ax,3); 
%plot the deconvolved sources
yloc = (yran(1)+range(yran)*0.05);
ind = find(impindepst(:,1)/sps >= xzoom(1) & impindepst(:,1)/sps <= xzoom(2));
% scatter(ax,impindepst(ind,1)/sps-xzoom(1),yloc*ones(size(impindepst(ind,1))),8,'r','filled');
scatter(ax,impindepst(ind,1)/sps-xzoom(1),yloc*ones(size(impindepst(ind,1))),10,'k');
oldc = colormap(ax,'kelicol');
newc = flipud(oldc);
colormap(ax,newc);
caxis(ax,xzoom);

yyaxis(ax,'right');
plot(ax,ircc/sps,yran(2)*rcc,'o','markersize',1,'Color',[.5 .5 .5]); % scale with data amp.
ylabel(ax,'Running CC','FontSize',10);
hold(ax,'off');

orient(f.fig,'landscape');
fname = 'lfetempvstremor.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
