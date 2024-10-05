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

cutout = 'ellipse';
ttol = 35;
ntol = 3;
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.',...
  num2str(ntol),'.pgc002.',cutout(1:4)));
tlen = trange(:,3)-trange(:,2);
nbst = size(trange,1);
tlensum = sum(tlen);

modname = 'timeoff_plfit_4thsta_160sps.mat';
planefit = load(strcat(rstpath, '/MAPS/',modname));
rmse = planefit.gof{4}.rmse;
offmax = round(2.0*rmse);

sps = 160;
ftrans = 'interpchao';
[imploc0, ~] = off2space002([0 0],sps,ftrans,0);  % a ref source at 0,0

%%%load data
% savefile = 'deconv_stats4th_allbstsig.mat';
savefile = 'deconv_stats4th_no23_allbstsig.mat';
load(strcat(rstpath, '/MAPS/',savefile));
% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
load(strcat(rstpath, '/MAPS/',savefile));
 
% keyboard

%%
%%%param for secondary sources removed
locxyprojall = allbstsig.locxyprojall;
tarvlsplstall = allbstsig.impindepall(:,1);
nsrc = allbstsig.nsrc;
imp = allbstsig.impindepall;
locxyprojalln = allbstnoi.locxyprojall;
tarvlsplstalln = allbstnoi.impindepall(:,1);
nsrcn = allbstnoi.nsrc;
impn = allbstnoi.impindepall;
supertstr = 'Secondary sources removed';
fnsuffix = [];

% %%%param for further checked at KLNB
% locxyprojall = allbstsig.locxyproj4thall;
% tarvlsplstall = allbstsig.impindep4thall(:,1);
% nsrc = allbstsig.nsrc4th;
% imp = allbstsig.impindep4thall;
% locxyprojalln = allbstnoi.locxyproj4thall;
% tarvlsplstalln = allbstnoi.impindep4thall(:,1);
% nsrcn = allbstnoi.nsrc4th;
% impn = allbstnoi.impindep4thall;
% supertstr = 'Further checked at KLNB';
% fnsuffix = '4th';

impuse = imp;
nsrcuse = nsrc;
fnsuffix2 = [];
% impuse = impn;
% nsrcuse = nsrcn;
% fnsuffix2 = 'noi';

%% using new ways to generate EXCLUSIVE clusters from each other
%%%determine the 'mmax' for which the resulting number of clusters is nonzero
timetype = 'tarvl';
mmax=getmmaxcluster(nbst,impuse,nsrcuse,sps,timetype);
%%%for a cluster, not only the time separation between N and N-m needs to
%%%be smaller than 'dtcut', but also the max time separation between each
%%%consecutive events needs to smaller than 0.25+0.125 s. 
%%%ie, doublet means a cluster of 2 events ONLY occur as doublets
[catclus,catclusbst,catimp,catuimp,catmedamp,catdtnnm]=...
  evtcluster_ex(nbst,impuse,nsrcuse,mmax,sps,timetype);

%% various types of distance between events in the clusters
mdistnn1l=[];
mdist2alll=[];
mdprojx1nn1l=[];
mdprojx12alll=[];
mdprojx2nn1l=[];
mdprojx22alll=[];
distnn1catl=[];
dist2allcatl=[];
dprojxy1nn1catl=[];
dprojxy12allcatl=[];
dprojxy2nn1catl=[];
dprojxy22allcatl=[];

minm = 6;
for m = minm: mmax
  % m=6;  % a cluster of m+1 events
  catclusm = catclus{m};
  ncluster = size(catclusm,1); %num of clusters, consecu. clusters may share events!

  %%%compute the distance w/i each cluster
  %this would take a while
  [mdistnn1,mdist2all,pca2vec,projang1,mdprojx1nn1,mdprojx12all,...
    projang2,mdprojx2nn1,mdprojx22all,distnn1cat,dist2allcat,...
    dprojxy1nn1cat,dprojxy12allcat,dprojxy2nn1cat,dprojxy22allcat]=...
    dist_evtcluster(catclusm,sps,ftrans,timetype);
  mdistnn1l = [mdistnn1l; mdistnn1];
  mdist2alll = [mdist2alll; mdist2all];
  mdprojx1nn1l = [mdprojx1nn1l; mdprojx1nn1];
  mdprojx12alll = [mdprojx12alll; mdprojx12all];
  mdprojx2nn1l = [mdprojx2nn1l; mdprojx2nn1];
  mdprojx22alll = [mdprojx22alll; mdprojx22all];
  distnn1catl = [distnn1catl; distnn1cat];
  dist2allcatl = [dist2allcatl; dist2allcat];
  dprojxy1nn1catl = [dprojxy1nn1catl; dprojxy1nn1cat];
  dprojxy12allcatl = [dprojxy12allcatl; dprojxy12allcat];
  dprojxy2nn1catl = [dprojxy2nn1catl; dprojxy2nn1cat];
  dprojxy22allcatl = [dprojxy22allcatl; dprojxy22allcat];

end

%%
widin = 7.5;  % maximum width allowed is 8.5 inches
htin = 3;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 3;
f = initfig(widin,htin,nrow,ncol);
pltxran = [0.06 0.98]; pltyran = [0.12 0.98]; % optimal axis location
pltxsep = 0.015; pltysep = 0.05;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

xran = [0.05 4];
yran = [0.05 4];
xtks = [0.1 0.5 1 2 4];
ytks = [0.05 0.1 0.5 1 2 4];

ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); ax.XScale='log'; ax.YScale='log';
plot(ax,xran,yran,'--','color',[.4 .4 .4],'linew',1);
scatter(ax,mdistnn1l,mdist2alll,20,'k','filled','markeredgec','w');
axis(ax,'equal');
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xtks);
yticks(ax,ytks);
xlabel(ax,'Med. dist. (km) of consecutive pairs');
ylabel(ax,'Med. dist. (km) of all pairs');
text(ax,0.99,0.05,'Absolute','Units','normalized','HorizontalAlignment',...
  'right','FontSize',9);
text(ax,0.02,0.94,'a','FontSize',10,'unit','normalized','HorizontalAlignment',...
  'left','EdgeColor','k','Margin',1,'backgroundcolor','w');

ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); ax.XScale='log'; ax.YScale='log';
plot(ax,xran,yran,'--','color',[.4 .4 .4],'linew',1);
scatter(ax,mdprojx1nn1l,mdprojx12alll,20,'k','filled','markeredgec','w');
axis(ax,'equal');
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xtks);
yticks(ax,ytks);
nolabels(ax,2);
xlabel(ax,'Med. dist. (km) of consecutive pairs');
text(ax,0.99,0.05,'Along short PCA axis','Units','normalized','HorizontalAlignment',...
  'right','FontSize',9);
text(ax,0.02,0.94,'b','FontSize',10,'unit','normalized','HorizontalAlignment',...
  'left','EdgeColor','k','Margin',1,'backgroundcolor','w');

ax=f.ax(3); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); ax.XScale='log'; ax.YScale='log';
plot(ax,xran,yran,'--','color',[.4 .4 .4],'linew',1);
scatter(ax,mdprojx2nn1l,mdprojx22alll,20,'k','filled','markeredgec','w');
axis(ax,'equal');
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xtks);
yticks(ax,ytks);
nolabels(ax,2);
xlabel(ax,'Med. dist. (km) of consecutive pairs');
text(ax,0.99,0.05,'Along SW-NE','Units','normalized','HorizontalAlignment',...
  'right','FontSize',9);
text(ax,0.02,0.94,'c','FontSize',10,'unit','normalized','HorizontalAlignment',...
  'left','EdgeColor','k','Margin',1,'backgroundcolor','w');

fname = strcat('dist_evtclus_sum_mgeq',num2str(minm+1),'.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));









