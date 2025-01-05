% plt_vrvssnr_bst.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the script in particular to plot a summary figure for the scattering
% of variance reduction vs the signal to noise ratio, or the histogram of
% variance reduction for several selected bursts or all bursts.
%
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/11/23
% Last modified date:   2024/11/23
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
trange = load(strcat(rstpath, '/MAPS/tdec.bstranbuf',num2str(ttol),'s.pgc002.',cutout(1:4)));
nbst = size(trange,1);
tlen = trange(:,3)-trange(:,2);
tlensum = sum(tlen);

modname = 'timeoff_plfit_4thsta_160sps.mat';
planefit = load(strcat(rstpath, '/MAPS/',modname));
rmse = planefit.gof{4}.rmse;
offmax = round(2.0*rmse);

stas=['PGC  '
  'SSIB '
  'SILB '
  'LZB  '
  'TWKB '
  'MGCB '
  'KLNB ']; % determine the trio and order, here the 1st sta is PGC
nsta=size(stas,1);         %  number of stations

sps = 160;

ftrans = 'interpchao';

detecttype = '1win';
% detecttype = [];

%%%load data
%choose the window length in sec for computing RCC
% rccwin = 0.25;
rccwin = 0.5;

if rccwin == 0.5
  savefile = 'deconv1win_stats4th_no23_allbstsig_indhi.mat';  %whole-win detection
%   savefile = 'deconv1win_stats4th_no23_allbstsig.mat';  %whole-win detection
elseif rccwin == 0.25
  savefile = 'deconv_stats4th_no23_allbstsig0.25s.mat';
end
load(strcat(rstpath, '/MAPS/',savefile));

%%%whether to save the figure
savefig = 1;
% savefig = 0;

%%
%%%param for secondary sources removed
vr = allbstsig.varred;
vrort = allbstsig.varredort;
vrvert = allbstsig.varredvert;
menv = allbstsig.menv;
menvort = allbstsig.menvort;
menvvert = allbstsig.menvvert;
% mamp = allbstsig.mamp;
snrf = allbstsig.snrf;
snrfort = allbstsig.snrfort;
snrfvert = allbstsig.snrfvert;
supertstr = '3-station';

%%%param for further checked at KLNB
vr4th = allbstsig.varred4th;
vrort4th = allbstsig.varredort4th;
vrvert4th = allbstsig.varredvert4th;
supertstr4th = '4-station';

%% uncomment if only looking at the 6 other bursts
indhi = [99; 100; 115; 164; 170; 174; 181];
vr = vr(indhi, :, :);
vrort = vrort(indhi, :, :);
vrvert = vrvert(indhi, :, :);
vr4th = vr4th(indhi, :, :);
vrort4th = vrort4th(indhi, :, :);
vrvert4th = vrvert4th(indhi, :, :);
menv = menv(indhi, :);
menvort = menvort(indhi, :);
menvvert = menvvert(indhi, :);

%%
% savefile = strcat('deconrst047.mat');
% if isfile(strcat(rstpath, '/MAPS/',savefile))
%   deconrst047 = load(strcat(rstpath, '/MAPS/',savefile));
% end
% snrf = deconrst047.snrf;
% snrfort = deconrst047.snrfort;
% snrfvert = deconrst047.snrfvert;

%% PLOT, SCATTER OF VR OF SEISMOGRAMS VS SNR OF TEMPLATES
widin = 5.5;  % maximum width allowed is 8.5 inches
htin = 6;   % maximum height allowed is 11 inches
nrow = 2;
ncol = 2;
f = initfig(widin,htin,nrow,ncol);
pltxran = [0.1 0.98]; pltyran = [0.1 0.98]; % optimal axis location
pltxsep = 0.1; pltysep = 0.08;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
pltsta=[1 2 3];
p=[];
%%%for others bursts, use open symbols 
for i = 1: size(vr,1)-1
  varred = squeeze(vr(i,:,:));
  varredort = squeeze(vrort(i,:,:));
  varredvert = squeeze(vrvert(i,:,:));
  p(1)=scatter(ax,snrf(pltsta),varred(pltsta,3),20,'ro');
  p(2)=scatter(ax,snrfort(pltsta),varredort(pltsta,3),25,'bs');
  p(3)=scatter(ax,snrfvert(pltsta),varredvert(pltsta,3),20,'k^');
end
%%%for burst 181, use solid symbols 
for i = size(vr,1):size(vr,1)
  varred = squeeze(vr(i,:,:));
  varredort = squeeze(vrort(i,:,:));
  varredvert = squeeze(vrvert(i,:,:));
  scatter(ax,snrf(pltsta),varred(pltsta,3),20*2.5,'ro','filled');
  scatter(ax,snrfort(pltsta),varredort(pltsta,3),25*3,'bs','filled');
  scatter(ax,snrfvert(pltsta),varredvert(pltsta,3),20*3,'k^','filled');
end
meanvr181 = mean(varred(pltsta,3))
meanvrort181 = mean(varredort(pltsta,3))
meanvrvert181 = mean(varredvert(pltsta,3))
xlabel(ax,'LFE template SNR');
ylabel(ax,'Variance reduction (%)');
%       title(ax,'Right after grouping');
text(ax,0.01,0.94,'Grouped sources',...
  'Units','normalized','HorizontalAlignment','left','fontsize',10);
text(ax,0.01,0.85,'7 selected bursts',...
  'Units','normalized','HorizontalAlignment','left','fontsize',10);
text(ax,0.98,0.94,'a','FontSize',10,'unit','normalized','HorizontalAlignment',...
  'right','EdgeColor','k','Margin',1,'backgroundcolor','w');
xlim(ax,[0 80]);
% ylim(ax,[0 100]);
ylim(ax,[-20 100]);

ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
pltsta=[1 2 3 7];
%%%for others bursts, use open symbols 
for i = 1: size(vr,1)-1
  varred4th = squeeze(vr4th(i,:,:));
  varredort4th = squeeze(vrort4th(i,:,:));
  varredvert4th = squeeze(vrvert4th(i,:,:));
  scatter(ax,snrf(pltsta),varred4th(pltsta,3),20,'ro');
  scatter(ax,snrfort(pltsta),varredort4th(pltsta,3),25,'bs');
  scatter(ax,snrfvert(pltsta),varredvert4th(pltsta,3),20,'k^');
end
%%%for burst 181, use solid symbols 
for i = size(vr,1):size(vr,1)
  varred4th = squeeze(vr4th(i,:,:));
  varredort4th = squeeze(vrort4th(i,:,:));
  varredvert4th = squeeze(vrvert4th(i,:,:));
  scatter(ax,snrf(pltsta),varred4th(pltsta,3),20*2.5,'ro','filled');
  scatter(ax,snrfort(pltsta),varredort4th(pltsta,3),25*3,'bs','filled');
  scatter(ax,snrfvert(pltsta),varredvert4th(pltsta,3),20*3,'k^','filled');
end
pltsta=[1 2 3];
meanvr4th181 = mean(varred4th(pltsta,3))
meanvrort4th181 = mean(varredort4th(pltsta,3))
meanvrvert4th181 = mean(varredvert4th(pltsta,3))
lgd=legend(ax,p,'Optimal','Orthogonal','Vertical','location','southeast');
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
xlabel(ax,'LFE template SNR');
ylabel(ax,'Variance reduction (%)');
%       title(ax,'4-station srcs');
text(ax,0.01,0.94,'4-station sources',...
  'Units','normalized','HorizontalAlignment','left','fontsize',10);
text(ax,0.01,0.85,'7 selected bursts',...
  'Units','normalized','HorizontalAlignment','left','fontsize',10);
text(ax,0.98,0.94,'b','FontSize',10,'unit','normalized','HorizontalAlignment',...
  'right','EdgeColor','k','Margin',1,'backgroundcolor','w');
xlim(ax,[0 80]);
% ylim(ax,[0 100]);
ylim(ax,[-20 100]);

%% FOR all bursts
if rccwin == 0.5
  savefile = 'deconv1win_stats4th_no23_allbstsig.mat';  %whole-win detection
elseif rccwin == 0.25
  savefile = 'deconv_stats4th_no23_allbstsig0.25s.mat';
end
load(strcat(rstpath, '/MAPS/',savefile));

% %%
%%%param for secondary sources removed
vr = allbstsig.varred;
vrort = allbstsig.varredort;
vrvert = allbstsig.varredvert;
menv = allbstsig.menv;
menvort = allbstsig.menvort;
menvvert = allbstsig.menvvert;
% mamp = allbstsig.mamp;
snrf = allbstsig.snrf;
snrfort = allbstsig.snrfort;
snrfvert = allbstsig.snrfvert;
supertstr = '3-station';

%%%param for further checked at KLNB
vr4th = allbstsig.varred4th;
vrort4th = allbstsig.varredort4th;
vrvert4th = allbstsig.varredvert4th;
supertstr4th = '4-station';

%% PLOT, MEAN VR OF ALL STATIONS, HISTOGRAM
% widin = 5.5;  % maximum width allowed is 8.5 inches
% htin = 3;   % maximum height allowed is 11 inches
% nrow = 1;
% ncol = 2;
% f = initfig(widin,htin,nrow,ncol);
% pltxran = [0.1 0.98]; pltyran = [0.15 0.98]; % optimal axis location
% pltxsep = 0.1; pltysep = 0.05;
% optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

ax=f.ax(3); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
pltsta=[1 2 3];
p=[];
varred = mean(squeeze(vr(:,pltsta,3)), 2);  %mean of a few stations
varredort = mean(squeeze(vrort(:,pltsta,3)), 2);
varredvert = mean(squeeze(vrvert(:,pltsta,3)), 2);
p(1)=histogram(ax,varred,'facecolor','r','BinWidth',2.5);  %total num is the num of bursts
p(2)=histogram(ax,varredort,'facecolor','b','BinWidth',2.5);
p(3)=histogram(ax,varredvert,'facecolor','k','BinWidth',2.5);
medvarred = median(varred)
medvarredort = median(varredort)
medvarredvert = median(varredvert)
xlabel(ax,'Mean variance reduction (%)');
ylabel(ax,'Count');
text(ax,0.01,0.94,'Grouped sources',...
  'Units','normalized','HorizontalAlignment','left','fontsize',10);
text(ax,0.01,0.85,'All bursts',...
  'Units','normalized','HorizontalAlignment','left','fontsize',10);
text(ax,0.98,0.94,'c','FontSize',10,'unit','normalized','HorizontalAlignment',...
  'right','EdgeColor','k','Margin',1,'backgroundcolor','w');
ylim(ax,[0 50]);
xlim(ax,[-20 100]);
xticks(ax,-20:20:100);

ax=f.ax(4); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
% pltsta=[1 2 3 7]; %if using the mean of all 4 stas
pltsta=[1 2 3]; %if using the mean of the trio stas
varred4th = mean(squeeze(vr4th(:,pltsta,3)), 2);
varredort4th = mean(squeeze(vrort4th(:,pltsta,3)), 2);
varredvert4th = mean(squeeze(vrvert4th(:,pltsta,3)), 2);
histogram(ax,varred4th,'facecolor','r','BinWidth',2.5);
histogram(ax,varredort4th,'facecolor','b','BinWidth',2.5);
histogram(ax,varredvert4th,'facecolor','k','BinWidth',2.5);
medvarred4th = median(varred4th)
medvarredort4th = median(varredort4th)
medvarredvert4th = median(varredvert4th)
lgd=legend(ax,p,'Optimal','Orthogonal','Vertical','location','east');
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
xlabel(ax,'Mean variance reduction (%)');
ylabel(ax,'Count');
text(ax,0.01,0.94,'4-station sources',...
  'Units','normalized','HorizontalAlignment','left','fontsize',10);
text(ax,0.01,0.85,'All bursts',...
  'Units','normalized','HorizontalAlignment','left','fontsize',10);
text(ax,0.98,0.94,'d','FontSize',10,'unit','normalized','HorizontalAlignment',...
  'right','EdgeColor','k','Margin',1,'backgroundcolor','w');
ylim(ax,[0 50]);
xlim(ax,[-20 100]);
xticks(ax,-20:20:100);

fname = strcat('vrvssnr_bsts',detecttype,'.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

keyboard



%% PLOT, MEAN VR OF ALL STATIONS VS MEAN ENVELOPE
widin = 5.5;  % maximum width allowed is 8.5 inches
htin = 3;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f = initfig(widin,htin,nrow,ncol);
pltxran = [0.1 0.98]; pltyran = [0.15 0.98]; % optimal axis location
pltxsep = 0.1; pltysep = 0.05;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
pltsta=[1 2 3];
mmenv = mean(menv(:,pltsta),2);
mmenvort = mean(menvort(:,pltsta),2);
mmenvvert = mean(menvvert(:,pltsta),2);
p(1)=scatter(ax,mmenv,varred,20,'ro');
p(2)=scatter(ax,mmenvort,varredort,24,'bs');
p(3)=scatter(ax,mmenvvert,varredvert,20,'k^');
lgd=legend(ax,p,'Optimal','Orthogonal','Vertical','location','southeast');
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
xlabel(ax,'Mean median envelope');
ylabel(ax,'Mean variance reduction (%)');
text(ax,0.01,0.94,'Grouped sources',...
  'Units','normalized','HorizontalAlignment','left','fontsize',10);
text(ax,0.98,0.94,'a','FontSize',10,'unit','normalized','HorizontalAlignment',...
  'right','EdgeColor','k','Margin',1,'backgroundcolor','w');
% ylim(ax,[0 100]);
ylim(ax,[-20 100]);

ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
pltsta=[1 2 3 7];
mmenv4th = mean(menv(:,pltsta),2);
mmenvort4th = mean(menvort(:,pltsta),2);
mmenvvert4th = mean(menvvert(:,pltsta),2);
p(1)=scatter(ax,mmenv4th,varred4th,20,'ro');
p(2)=scatter(ax,mmenvort4th,varredort4th,24,'bs');
p(3)=scatter(ax,mmenvvert4th,varredvert4th,20,'k^');
xlabel(ax,'Mean median envelope');
ylabel(ax,'Mean variance reduction (%)');
text(ax,0.01,0.94,'4-station sources',...
  'Units','normalized','HorizontalAlignment','left','fontsize',10);
text(ax,0.98,0.94,'b','FontSize',10,'unit','normalized','HorizontalAlignment',...
  'right','EdgeColor','k','Margin',1,'backgroundcolor','w');
% ylim(ax,[0 100]);
ylim(ax,[-20 100]);





