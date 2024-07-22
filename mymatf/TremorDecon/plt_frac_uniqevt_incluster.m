% plt_frac_uniqevt_incluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the script in particular to plot the comparison of the fraction of 
% the catalog in terms of unique events in different type of clusters between
% 3-sta data, 3-sta noise, 4-sta data and 4-sta noise.
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/04/04
% Last modified date:   2024/04/04
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

cutout = 'ellipse';
ttol = 35;
ntol = 3;
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.pgc002.',cutout(1:4)));
% trange = trange(1:end-1,:);
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
load(strcat(rstpath, '/MAPS/',savefile));
% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
load(strcat(rstpath, '/MAPS/',savefile));
  
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
supertstr = '3-station';

%%%param for further checked at KLNB
locxyprojall4th = allbstsig.locxyproj4thall;
tarvlsplstall4th = allbstsig.impindep4thall(:,1);
nsrc4th = allbstsig.nsrc4th;
imp4th = allbstsig.impindep4thall;
locxyprojalln4th = allbstnoi.locxyproj4thall;
tarvlsplstalln4th = allbstnoi.impindep4thall(:,1);
nsrcn4th = allbstnoi.nsrc4th;
impn4th = allbstnoi.impindep4thall;
supertstr4th = '4-station';

%% compution of fractions for diff catalogs
%%%determine the 'mmax' for which the resulting number of clusters is nonzero
timetype = 'tarvl';
mmax=getmmaxcluster(nbst,imp,nsrc,sps,timetype);
%%%for a cluster, not only the time separation between N and N-m needs to
%%%be smaller than 'dtcut', but also the max time separation between each
%%%consecutive events needs to smaller than 0.25+0.125 s. 
%%%ie, doublet means a cluster of 2 events ONLY occur as doublets
[catclus,catclusbst,catimp,catuimp,catmedamp,catdtnnm]=...
  evtcluster_ex(nbst,imp,nsrc,mmax,sps,timetype);

%%%fraction of unique events ONLY occurring as certain clusters
[fracuimp,nuimp,fracuimpsum]=frac_uniqevt_incluster2(catuimp,catclus,imp,nsrc,mmax);

%%%% FOR NOISE, 3-station catalog
mmaxn=getmmaxcluster(nbst,impn,nsrcn,sps,timetype);
[catclusn,catclusbstn,catimpn,catuimpn,catmedampn,catdtnnmn]=...
  evtcluster_ex(nbst,impn,nsrcn,mmaxn,sps,timetype);
[fracuimpn,~,fracuimpnsum]=frac_uniqevt_incluster2(catuimpn,catclusn,impn,nsrcn,mmaxn);

%%%% FOR DATA, 4-station catalog
mmax4th=getmmaxcluster(nbst,imp4th,nsrc4th,sps,timetype);
[catclus4th,catclusbst4th,catimp4th,catuimp4th,catmedamp4th,catdtnnm4th]=...
  evtcluster_ex(nbst,imp4th,nsrc4th,mmax4th,sps,timetype);
[fracuimp4th,~,fracuimp4thsum]=frac_uniqevt_incluster2(catuimp4th,catclus4th,imp4th,nsrc4th,mmax4th);

%%%% FOR NOISE, 4-station catalog
mmaxn4th=getmmaxcluster(nbst,impn4th,nsrcn4th,sps,timetype);
[catclusn4th,catclusbstn4th,catimpn4th,catuimpn4th,catmedampn4th,catdtnnmn4th]=...
  evtcluster_ex(nbst,impn4th,nsrcn4th,mmaxn4th,sps,timetype);
[fracuimpn4th,~,fracuimpn4thsum]=frac_uniqevt_incluster2(catuimpn4th,catclusn4th,impn4th,nsrcn4th,mmaxn4th);


%%
%%%PLOT
nrow = 1; % rows and cols of subplots in each figure
ncol = 2; 
widin = 7; % size of each figure
htin = 3.5;
pltxran = [0.1 0.98]; pltyran = [0.15 0.95];
pltxsep = 0.07; pltysep = 0.03;
f = initfig(widin,htin,nrow,ncol);
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

ax = f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
p(1)=plot(ax,1:mmax+1,fracuimpsum,'k-','Linew',1.5,'marker','o','markersize',3,...
  'markerfacec','k');
p(2)=plot(ax,1:mmaxn+1,fracuimpnsum,'r-','Linew',1.5,'marker','o','markersize',3,...
  'markerfacec','r');
legend(ax,p,'Data','Noise','Location','east');
ax.YAxis.Scale = 'log'; %make y axis log10 scale
xlim(ax,[1 mmax+2]);
ylim(ax,[5e-2 1e2]);
xticks(ax,0:2:mmax+2);
yticks(ax,[0.1 0.2 0.5 1 2 5 10 20 50 100]);
% xlabel(ax,'# of events in cluster','FontSize',10);
% ylabel(ax,'% of catalog in such clusters','FontSize',10);
xlabel(ax,'m (# of events in cluster)','FontSize',10);
ylabel(ax,'% of catalog in clusters of >=m events','FontSize',10);
text(ax,0.98,0.95,supertstr,'HorizontalAlignment','right','Units','normalized',...
  'fontsize',10);
text(ax,0.02,0.05,'a','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);

ax = f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
p(1)=plot(ax,1:mmax4th+1,fracuimp4thsum,'k-','Linew',1.5,'marker','o',...
  'markersize',3,'markerfacec','k');
p(2)=plot(ax,1:mmaxn4th+1,fracuimpn4thsum,'r-','Linew',1.5,'marker','o',...
  'markersize',3,'markerfacec','r');
ax.YAxis.Scale = 'log'; %make y axis log10 scale
xlim(ax,[1 mmax+2]);
ylim(ax,[5e-2 1e2]);
xticks(ax,0:2:mmax+2);
yticks(ax,[0.1 0.2 0.5 1 2 5 10 20 50 100]);
text(ax,0.98,0.95,supertstr4th,'HorizontalAlignment','right','Units','normalized',...
  'fontsize',10);
text(ax,0.02,0.05,'b','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);

fname = strcat('fracunievtclus.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));





