% migoflfes.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to identify the potential migrations of LFEs among
% the burst windows we have analyzed for deconvolution. 
%
% --To recognize which ones are potential migrations, we will rely on the
% whole-win catalog. 1st step is to look by eyes, hopefully it can be quantified
% to some extent by Pearson correlation?
% --To tell if the discarded sources from 4-sta check are consistent with
% the migrating pattern of the preserved, we will rely on the
% whole-win catalog. 
% --to analyse the distribution of the resduals from linear fit along the 
% propagation direction, we will rely on the short-win catalog. 
% --For potential migrations, how much does the prop. direction vary between 
% preserved and discarded srcs? And the prop. speed?
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/04/29
% Last modified date:   2024/04/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
% format short e   % Set the format to 5-digit floating point
clear
clc
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, resol] = pixelperinch(1);

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');

%load the tremor bursts with buffered ranges used in deconvolution
cutout = 'ellipse';
ttol = 35;
ntol = 3;
trange = load(strcat(rstpath, '/MAPS/tdec.bstranbuf',num2str(ttol),'s.pgc002.',cutout(1:4)));
nbst = size(trange,1);
tlen = trange(:,3)-trange(:,2);
tlensum = sum(tlen);

modname = 'timeoff_plfit_4thsta_160sps.mat';
planefit = load(strcat(rstpath, '/MAPS/',modname));
rmse = planefit.gof{4}.rmse;
offmax = round(2.0*rmse);

sps = 160;

ftrans = 'interpchao';

%%%load data
% savefile = 'deconv_stats4th_allbstsig.mat';
savefile = 'deconv_stats4th_no23_allbstsig.mat';
allsig = load(strcat(rstpath, '/MAPS/',savefile));
% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
allnoi = load(strcat(rstpath, '/MAPS/',savefile));

% savefile = 'deconv1win_stats4th_allbstsig.mat';
savefile = 'deconv1win_stats4th_no23_allbstsig.mat';
allsig1w = load(strcat(rstpath, '/MAPS/',savefile));
% savefile = 'deconv1win_stats4th_allbstnoi.mat';
savefile = 'deconv1win_stats4th_no23_allbstnoi.mat';
allnoi1w = load(strcat(rstpath, '/MAPS/',savefile));

% keyboard

%%
fitstats=allsig1w.allbstsig.fitstats;
fitstats4th=allsig1w.allbstsig.fitstats4th;
fitstatsrem=allsig1w.allbstsig.fitstatsrem;
nsrc = allsig.allbstsig.nsrc;
nsrc4th = allsig.allbstsig.nsrc4th;

for i=1:nbst
  stats=fitstats{i};
  if ~isempty(stats)
    slopese(i,1)=stats.slopese;
    pear(i,1)=stats.pearwt;
    rsquare(i,1)=stats.gof.rsquare;
    adjrsquare(i,1)=stats.gof.adjrsquare;
    slope(i,1)=stats.slope;
  else
    slopese(i,1)=nan;
    pear(i,1)=nan;
    rsquare(i,1)=nan;
    adjrsquare(i,1)=nan;
    slope(i,1)=nan;
  end
end
  
% %%
f=initfig(8,8,2,2);
ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
histogram(ax,slope*1e3,'binw',1);
xlabel(ax,'Slope (m/s)');
ylabel(ax,'Count');
xlim(ax,[0 50]);
text(ax,0.98,0.95,'ALL','HorizontalAlignment','right','Units','normalized');

ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
histogram(ax,slopese*1e3,'binw',1);
xlabel(ax,'SE of slope (m/s)');
ylabel(ax,'Count');
xlim(ax,[0 50]);
text(ax,0.98,0.95,'ALL','HorizontalAlignment','right','Units','normalized');

ax=f.ax(3); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
histogram(ax,pear,'binw',0.05);
xlabel(ax,'Pearson coeff');
ylabel(ax,'Count');
text(ax,0.98,0.95,'ALL','HorizontalAlignment','right','Units','normalized');

ax=f.ax(4); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
histogram(ax,rsquare,'binw',0.05,'facecolor','b');
histogram(ax,adjrsquare,'binw',0.05,'facecolor','r');
xlabel(ax,'R^2');
ylabel(ax,'Count');
text(ax,0.98,0.95,'ALL','HorizontalAlignment','right','Units','normalized');
legend(ax,'Original','DOF-adjusted');

%%
indmig = [14; 23; 25; 31; 33; 35; 42; 43; 44; 48; 49; 50; 53; 57; 60; 62; 66;
          67; 80; 82; 85; 91; 96; 99; 100; 104; 106; 110; 112; 115; 117; 122;
          138; 146; 151; 153; 164; 170; 172; 174; 181; 184; 194; 195];

indhi = [23; 35; 53; 60; 99; 100; 115; 164; 170; 174; 181; 184];

f=initfig(8,8,2,2);
ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
scatter(ax,slopese*1e3,adjrsquare,15,'k','filled');
scatter(ax,slopese(indmig)*1e3,adjrsquare(indmig),15,'b','filled');

xlabel(ax,'SE of slope (m/s)');
ylabel(ax,'Adjusted R^2');
xlim(ax,[0 50]);

ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
scatter(ax,slopese*1e3,rsquare,15,'k','filled');
scatter(ax,slopese(indmig)*1e3,rsquare(indmig),15,'b','filled');
xlabel(ax,'SE of slope (m/s)');
ylabel(ax,'R^2');
xlim(ax,[0 50]);

ax=f.ax(3); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
scatter(ax,slopese*1e3,pear,15,'k','filled');
scatter(ax,slopese(indmig)*1e3,pear(indmig),15,'b','filled');
xlabel(ax,'SE of slope (m/s)');
ylabel(ax,'Pearson coeff');
xlim(ax,[0 50]);

%%
f=initfig(8.4,3.5,1,3);
pltxran = [0.08 0.98]; pltyran = [0.12 0.98];
pltxsep = 0.08; pltysep = 0.03;
optaxpos(f,1,3,pltxran,pltyran,pltxsep,pltysep);

ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
scatter(ax,slope*1e3,adjrsquare,15,slopese*1e3,'o','filled','MarkerEdgeColor','none');
colormap(ax,'plasma');
c=colorbar(ax,'northoutside');
c.Label.String = 'SE of slope (m/s)';
caxis(ax,[0 5]);
xlim(ax,[0 50]);
scatter(ax,slope(indmig)*1e3,adjrsquare(indmig),30,slopese(indmig)*1e3,'^',...
  'filled','MarkerEdgeColor','k');
xlabel(ax,'Slope (m/s)');
ylabel(ax,'Adjusted R^2');

ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
scatter(ax,slope*1e3,rsquare,15,slopese*1e3,'o','filled','MarkerEdgeColor','none');
colormap(ax,'plasma');
c=colorbar(ax,'northoutside');
c.Label.String = 'SE of slope (m/s)';
caxis(ax,[0 5]);
xlim(ax,[0 50]);
scatter(ax,slope(indmig)*1e3,rsquare(indmig),30,slopese(indmig)*1e3,'^',...
  'filled','MarkerEdgeColor','k');
xlabel(ax,'Slope (m/s)');
ylabel(ax,'R^2');

ax=f.ax(3); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
scatter(ax,slope*1e3,pear,15,slopese*1e3,'o','filled','MarkerEdgeColor','none');
colormap(ax,'plasma');
c=colorbar(ax,'northoutside');
c.Label.String = 'SE of slope (m/s)';
caxis(ax,[0 5]);
xlim(ax,[0 50]);
scatter(ax,slope(indmig)*1e3,pear(indmig),30,slopese(indmig)*1e3,'^',...
  'filled','MarkerEdgeColor','k');
xlabel(ax,'Slope (m/s)');
ylabel(ax,'Pearson coeff');

%%
f=initfig(8.4,3.5,1,3);
pltxran = [0.08 0.98]; pltyran = [0.12 0.98];
pltxsep = 0.08; pltysep = 0.03;
optaxpos(f,1,3,pltxran,pltyran,pltxsep,pltysep);

ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
scatter(ax,slope*1e3,adjrsquare,15,slopese*1e3,'o','filled','MarkerEdgeColor','none');
colormap(ax,'plasma');
c=colorbar(ax,'northoutside');
c.Label.String = 'SE of slope (m/s)';
caxis(ax,[0 5]);
xlim(ax,[0 50]);
scatter(ax,slope(indhi)*1e3,adjrsquare(indhi),30,slopese(indhi)*1e3,'^',...
  'filled','MarkerEdgeColor','k');
xlabel(ax,'Slope (m/s)');
ylabel(ax,'Adjusted R^2');

ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
scatter(ax,slope*1e3,rsquare,15,slopese*1e3,'o','filled','MarkerEdgeColor','none');
colormap(ax,'plasma');
c=colorbar(ax,'northoutside');
c.Label.String = 'SE of slope (m/s)';
caxis(ax,[0 5]);
xlim(ax,[0 50]);
scatter(ax,slope(indhi)*1e3,rsquare(indhi),30,slopese(indhi)*1e3,'^',...
  'filled','MarkerEdgeColor','k');
xlabel(ax,'Slope (m/s)');
ylabel(ax,'R^2');

ax=f.ax(3); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
scatter(ax,slope*1e3,pear,15,slopese*1e3,'o','filled','MarkerEdgeColor','none');
colormap(ax,'plasma');
c=colorbar(ax,'northoutside');
c.Label.String = 'SE of slope (m/s)';
caxis(ax,[0 5]);
xlim(ax,[0 50]);
scatter(ax,slope(indhi)*1e3,pear(indhi),30,slopese(indhi)*1e3,'^',...
  'filled','MarkerEdgeColor','k');
xlabel(ax,'Slope (m/s)');
ylabel(ax,'Pearson coeff');

%%
for i=1:nbst
  stats=fitstats{i};
  if ~isempty(stats)
    angrmse(i,1)=stats.angrmse;
    pear(i,1)=stats.pearwt;
    slope(i,1)=stats.slope;
  else
    angrmse(i,1)=nan;
    pear(i,1)=nan;
    slope(i,1)=nan;
  end
  stats4th=fitstats4th{i};
  if ~isempty(stats4th)
    angrmse4th(i,1)=stats4th.angrmse;
    pear4th(i,1)=stats4th.pearwt;
    slope4th(i,1)=stats4th.slope;
  else
    angrmse4th(i,1)=nan;
    pear4th(i,1)=nan;
    slope4th(i,1)=nan;
  end
  statsrem=fitstatsrem{i};
  if ~isempty(statsrem)
    angrmserem(i,1)=statsrem.angrmse;
    pearrem(i,1)=statsrem.pearwt;
    sloperem(i,1)=statsrem.slope;
  else
    angrmserem(i,1)=nan;
    pearrem(i,1)=nan;
    sloperem(i,1)=nan;
  end
end

%direction for ALL bursts
f=initfig(4,4,1,1);
ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
histogram(ax,angrmse,'binw',5);
histogram(ax,angrmse4th,'binw',5,'facecolor','r');
xlabel(ax,'Direction (^o)');
ylabel(ax,'Count');
text(ax,0.98,0.95,'ALL','HorizontalAlignment','right','Units','normalized');

%speed for ALL bursts
f=initfig(8,4,1,2);
ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
histogram(ax,slope*1e3,'binw',1);
histogram(ax,slope4th*1e3,'binw',1,'facecolor','r');
xlabel(ax,'Speed (m/s)');
ylabel(ax,'Count');
text(ax,0.98,0.95,'ALL','HorizontalAlignment','right','Units','normalized');
ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
histogram(ax,slope*1e3,'binw',1);
histogram(ax,slope4th*1e3,'binw',1,'facecolor','r');
xlabel(ax,'Speed (m/s)');
ylabel(ax,'Count');
xlim(ax,[0 50]);
text(ax,0.98,0.95,'ALL','HorizontalAlignment','right','Units','normalized');

% fname = strcat('aaa.pdf');
% print(f.fig,'-dpdf',...
%   strcat('/home/chaosong/Pictures/',fname));
% keyboard


%% indices of a few representative migrations
indmig = [14; 23; 25; 31; 33; 35; 42; 43; 44; 48; 49; 50; 53; 57; 60; 62; 66;
          67; 80; 82; 85; 91; 96; 99; 100; 104; 106; 110; 112; 115; 117; 122;
          138; 146; 151; 153; 164; 170; 172; 174; 181; 184; 194; 195];

indhi = [23; 35; 53; 60; 99; 100; 115; 164; 170; 174; 181; 184];

%direction for SELECTED bursts
f=initfig(4,4,1,1);
ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
histogram(ax,angrmse(indmig),'binw',5);
histogram(ax,angrmse4th(indmig),'binw',5,'facecolor','r');
xlabel(ax,'Direction (^o)');
ylabel(ax,'Count');
text(ax,0.98,0.95,'Decent','HorizontalAlignment','right','Units','normalized');
legend(ax,'3-sta','4-sta','Location','south');

%speed for SELECTED bursts
f=initfig(8,4,1,2);
ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
histogram(ax,slope(indmig)*1e3,'binw',1);
histogram(ax,slope4th(indmig)*1e3,'binw',1,'facecolor','r');
xlabel(ax,'Speed (m/s)');
ylabel(ax,'Count');
text(ax,0.98,0.95,'Decent','HorizontalAlignment','right','Units','normalized');
ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
histogram(ax,slope(indhi)*1e3,'binw',1);
histogram(ax,slope4th(indhi)*1e3,'binw',1,'facecolor','r');
xlabel(ax,'Speed (m/s)');
ylabel(ax,'Count');
% xlim(ax,[0 50]);
text(ax,0.98,0.95,'High-quality','HorizontalAlignment','right','Units','normalized');
% fname = strcat('aaa.pdf');
% print(f.fig,'-dpdf',...
%   strcat('/home/chaosong/Pictures/',fname));
% keyboard

indplt = indmig;
nmig = length(indplt);

angdiff1=[]; angdiff2=[]; angdiff3=[]; slopediff1=[]; slopediff2=[]; slopediff3=[];
for i = 1:nmig
  angdiff1(i,1) = angrmse(indplt(i))-angrmse4th(indplt(i));
  angdiff2(i,1) = angrmse(indplt(i))-angrmserem(indplt(i));
  angdiff3(i,1) = angrmse4th(indplt(i))-angrmserem(indplt(i));
  slopediff1(i,1) = slope(indplt(i))/slope4th(indplt(i));
  slopediff2(i,1) = slope(indplt(i))/sloperem(indplt(i));
  slopediff3(i,1) = slope4th(indplt(i))/sloperem(indplt(i));
end

%%
nrow = 1; % rows and cols of subplots in each figure
ncol = 2;
widin = 6.5; % size of each figure
htin = 3.2;
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.08 0.98]; pltyran = [0.15 0.98];
pltxsep = 0.1; pltysep = 0.03;
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
scatter(ax,1:nmig, angdiff1,15,'k','filled'); hold on
% scatter(ax,1:nmig, angdiff2,15,'b','filled');
scatter(ax,1:nmig, angdiff3,15,'b','filled');
ylim(ax,[-30 30]);
xlabel(ax,'Burst #');
ylabel(ax,sprintf('Difference in prop. direction (%c)',char(176)));
legend(ax,'3-sta vs. 4-sta (preserved)','preserved vs. discarded');
text(ax,0.98,0.12,sprintf('%d%% w/i \x00B1 5%c, %d%% w/i \x00B1 10%c',...
  round(sum(abs(angdiff1)<=5)/nmig*100),char(176),...
  round(sum(abs(angdiff1)<=10)/nmig*100),char(176)),...
  'Units','normalized','HorizontalAlignment','right','color','k'); %\x2264
text(ax,0.98,0.05,sprintf('%d%% w/i \x00B1 5%c, %d%% w/i \x00B1 10%c',...
  round(sum(abs(angdiff3)<=5)/nmig*100),char(176),...
  round(sum(abs(angdiff3)<=10)/nmig*100),char(176)),...
  'Units','normalized','HorizontalAlignment','right','color','b');
% text(ax,0.98,0.95,'Decent','HorizontalAlignment','right','Units','normalized');
text(ax,0.02,0.95,'a','FontSize',11,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');

ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
ax.YScale = 'log';
scatter(ax,1:nmig, slopediff1,15,'k','filled'); hold on
% scatter(1:nmig, angdiff2,15,'b','filled');
scatter(ax,1:nmig, slopediff3,15,'b','filled');
% ylim(ax,[-30 30]);
xlabel(ax,'Burst #');
ylabel(ax,'Ratio of prop. speed');
% legend(ax,'3-sta / 4-sta (preserved)','preserved / discarded','location','northwest');
text(ax,0.98,0.12,sprintf('%d%% w/i [0.8,1.2], %d%% w/i [0.5,2]',...
  round(sum(slopediff1>=0.8 & slopediff1<=1.2)/nmig*100),...
  round(sum(slopediff1>=0.5 & slopediff1<=2)/nmig*100)),...
  'Units','normalized','HorizontalAlignment','right','color','k');
text(ax,0.98,0.05,sprintf('%d%% w/i [0.8,1.2], %d%% w/i [0.5,2]',...
  round(sum(slopediff3>=0.8 & slopediff3<=1.2)/nmig*100),...
  round(sum(slopediff3>=0.5 & slopediff3<=2)/nmig*100)),...
  'Units','normalized','HorizontalAlignment','right','color','b');
% text(ax,0.98,0.25,'Decent','HorizontalAlignment','right','Units','normalized');
% text(ax,0.98,0.25,'High-quality','HorizontalAlignment','right','Units','normalized');
ylim(ax,[0.25 4]);
yticks(ax,[0.25 0.5 0.8 1 1.2 2 4]);
text(ax,0.02,0.95,'b','FontSize',11,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');

fname = strcat('lfemigstats.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
% keyboard

%%
fitstats=allsig.allbstsig.fitstats;
fitstats4th=allsig.allbstsig.fitstats4th;
fitstatsrem=allsig.allbstsig.fitstatsrem;
nsrc = allsig.allbstsig.nsrc;
nsrc4th = allsig.allbstsig.nsrc4th;

pear=zeros(nbst,1);
pear4th=zeros(nbst,1);
resi=cell(nbst,1);
resi4th=cell(nbst,1);
for i=1:nbst
  stats=fitstats{i};
  pear(i)=stats.pearwt;
  resi{i}=stats.output.residuals;
  stats4th=fitstats4th{i};
  pear4th(i)=stats4th.pearwt;
  resi4th{i}=stats4th.output.residuals;
end

pearhi = prctile(pear,90);
indhi = find(pear>=pearhi);

%%
%%%let's select 12 of them to plot the linear fit residual
indplt = [23; 35; 53; 60; 99; 100; 115; 164; 170; 174; 181; 184];
nrow = 3; % rows and cols of subplots in each figure
ncol = 4;
widin = 8.4; % size of each figure
htin = 6;
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.06 0.98]; pltyran = [0.06 0.98];
pltxsep = 0.04; pltysep = 0.04;
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

xran=[-3 3];
yran=[0 0.16];
binw=0.1;
X=xran(1)+binw/2:binw:xran(2)-binw/2;
edges=X;
for i = 1: length(indplt)
  resplt = resi{indplt(i)};
  resplt4th = resi4th{indplt(i)};
  ax=f.ax(i); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
  h1=histogram(ax,resplt,'BinEdges',edges,'normalization','probability',...
    'Facec','b','edgec','none');
  h2=histogram(ax,resplt4th,'BinEdges',edges,'normalization','probability',...
    'Facec','r','edgec','none');
  [mu1,sigma1] = normfit(resplt);
  Y1=normpdf(X,mu1,sigma1)*binw;
  plot(ax,X,Y1,'Color','b','linew',2);
  [mu2,sigma2] = normfit(resplt4th);
  Y2=normpdf(X,mu2,sigma2)*binw;
  plot(ax,X,Y2,'Color','r','linew',2);
  text(ax,0.02,0.80,sprintf('\\mu=%.2f;\n\\sigma=%.2f km',mu1,sigma1),'Units',...
    'normalized','HorizontalAlignment','left','FontSize',8,'color','b');
  text(ax,0.65,0.80,sprintf('\\mu=%.2f;\n\\sigma=%.2f km',mu2,sigma2),'Units',...
    'normalized','HorizontalAlignment','left','FontSize',8,'color','r');
  text(ax,0.98,0.95,sprintf('burst #%d, %d/%d evts',indplt(i),nsrc(indplt(i)),...
    nsrc4th(indplt(i))),'Units','normalized','HorizontalAlignment','right',...
    'FontSize',8);
  xlim(ax,xran);
  ylim(ax,yran);
  if i == (nrow-1)*ncol+1
    xlabel(ax,'Residual (km)');
    ylabel(ax,'Probability');
  end
end








