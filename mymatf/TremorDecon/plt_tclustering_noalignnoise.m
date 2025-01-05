% plt_tclustering_noalignnoise.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the script in particular to plot a comparison figure of time clustering
% between the noise catalogs that were obtained from noise aligned to as good
% as we can over 25-s wins, and from noise not aligned at all.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/11/12
% Last modified date:   2024/11/12
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

sps = 160;

ftrans = 'interpchao';

%%%load data
%choose the window length in sec for computing RCC 
% rccwin = 0.25;
rccwin = 0.5;

if rccwin == 0.5
%   savefile = 'deconv1win_stats4th_no23_allbstnoi.mat';  %whole-win detection
  % savefile = 'deconv_stats4th_allbstnoi.mat';
  savefile = 'deconv_stats4th_no23_allbstnoi.mat';
elseif rccwin == 0.25
  savefile = 'deconv_stats4th_no23_allbstnoi0.25s.mat';
end
load(strcat(rstpath, '/MAPS/',savefile));

%%%whether to save the figure
% savefig = 1;
savefig = 0;

%%
%%%param for secondary sources removed
nsrcn = allbstnoi.nsrc;
impn = allbstnoi.impindepall;
supertstr = '3-station';

%%%param for further checked at KLNB
nsrcn4th = allbstnoi.nsrc4th;
impn4th = allbstnoi.impindep4thall;
supertstr4th = '4-station';

%% inter-event time, ie. time from each to its preceeding event, N to N-1
m = 1;
nbst = size(trange,1);
dtcut = 0.25*m+0.125;
[ampn{1},dtintern{1}]=med_amp_incluster(nbst,impn,nsrcn,m);
[ampn4th{1},dtintern4th{1}]=med_amp_incluster(nbst,impn4th,nsrcn4th,m);


%% load syn noise catalogs from noise not aligned at all
if rccwin == 0.5
  savefile = 'deconv_stats4th_no23_allbstnoi_noalign.mat';
end
load(strcat(rstpath, '/MAPS/',savefile));

%%%param for secondary sources removed
nsrcn = allbstnoi.nsrc;
impn = allbstnoi.impindepall;
supertstr = '3-station';

%%%param for further checked at KLNB
nsrcn4th = allbstnoi.nsrc4th;
impn4th = allbstnoi.impindep4thall;
supertstr4th = '4-station';

%%%inter-event time, ie. time from each to its preceeding event, N to N-1
m = 1;
nbst = size(trange,1);
dtcut = 0.25*m+0.125;
[ampn{2},dtintern{2}]=med_amp_incluster(nbst,impn,nsrcn,m);
[ampn4th{2},dtintern4th{2}]=med_amp_incluster(nbst,impn4th,nsrcn4th,m);


%% plot the inter-event time distribution separately
widin = 5.5;  % maximum width allowed is 8.5 inches
htin = 3;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f = initfig(widin,htin,nrow,ncol); %initialize fig
pltxran = [0.1 0.98]; pltyran = [0.15 0.95]; % optimal axis location
pltxsep = 0.08; pltysep = 0.05;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

%%% histogram of inter-event time, linear scale
xran = [0 2];
yran = [0 0.15];

binwdt = 0.05;
nx = round(xran(2)/binwdt)+1;
binedge = (xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2)';
bincnt = (xran(1): binwdt: xran(2))';

ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
patdtcut = [0 yran(2);
          dtcut yran(2);
          dtcut yran(1);
          0 yran(1);
          0 yran(2)];
patch(ax,patdtcut(:,1),patdtcut(:,2),'k','Facealpha',0.15,'edgecolor','none');
dtinterplt = dtintern{1};
Nn=histcounts(dtinterplt(:,1)/sps,binedge,'normalization','count');
Nnn = reshape(Nn,[],1)./length(dtintern{1});
fracn = sum(dtinterplt(:,1)/sps<=dtcut)/length(dtinterplt);
p(1)=plot(ax,bincnt,Nnn,'color','r','LineWidth',1.5);
label{1}='Synthetic noise';
dtinterplt = dtintern{2};
Nn=histcounts(dtinterplt(:,1)/sps,binedge,'normalization','count');
Nnn = reshape(Nn,[],1)./length(dtintern{1});
fracn = sum(dtinterplt(:,1)/sps<=dtcut)/length(dtinterplt);
p(2)=plot(ax,bincnt,Nnn,'color','b','LineWidth',1.5);
label{2}='Noise not aligned';
% fracdif = (sum(dtarvl/sps<=dtcut)-sum(dtarvln/sps<=dtcut)) / ...
%   (length(dtarvl)-length(dtarvln));
% p(3)=plot(ax,cnt,Nn-Nnn,'color','k','LineWidth',1.5);
% label{3}='Data - Synthetic noise';
% text(ax,0.25,0.65,sprintf('%.2f',fracn),'Color','r','HorizontalAlignment',...
%   'left','Units','normalized');
% text(ax,0.7,0.56,sprintf('%.2f',fracdif),'Color','k','HorizontalAlignment','left','Units','normalized');
text(ax,0.5,0.94,supertstr,'HorizontalAlignment','left','Units','normalized',...
  'fontsize',10);
text(ax,0.98,0.94,'a','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w','HorizontalAlignment','right');
legend(ax,p,label,'Location','east');
ylabel(ax,'Normalized count');
xlabel(ax,'Inter-event time (s)','FontSize',10);
xlim(ax,xran);
% if nsep == 1
%   yran = [0 0.2];
% elseif nsep == 2
%   yran = [0 0.06];
% else
%   yran = [0 0.04];
% end
ylim(ax,yran);
longticks(ax,2);
hold(ax,'off');
% keyboard

ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
patdtcut = [0 yran(2);
          dtcut yran(2);
          dtcut yran(1);
          0 yran(1);
          0 yran(2)];
patch(ax,patdtcut(:,1),patdtcut(:,2),'k','Facealpha',0.15,'edgecolor','none');
dtinterplt = dtintern4th{1};
Nn4th=histcounts(dtinterplt(:,1)/sps,binedge,'normalization','count');
Nn4thn = reshape(Nn4th,[],1)./length(dtintern{1});
fracn4th = sum(dtinterplt(:,1)/sps<=dtcut)/length(dtinterplt);
p(1)=plot(ax,bincnt,Nn4thn,'color','r','LineWidth',1.5);
label{1}='Synthetic noise';
dtinterplt = dtintern4th{2};
Nn4th=histcounts(dtinterplt(:,1)/sps,binedge,'normalization','count');
Nn4thn = reshape(Nn4th,[],1)./length(dtintern{1});
fracn4th = sum(dtinterplt(:,1)/sps<=dtcut)/length(dtinterplt);
p(2)=plot(ax,bincnt,Nn4thn,'color','b','LineWidth',1.5);
label{2}='Noise not aligned';
% fracdif = (sum(dtarvl/sps<=dtcut)-sum(dtarvln/sps<=dtcut)) / ...
%   (length(dtarvl)-length(dtarvln));
% p(3)=plot(ax,cnt,Nn-Nnn,'color','k','LineWidth',1.5);
% label{3}='Data - Synthetic noise';
% text(ax,0.25,0.65,sprintf('%.2f',fracn4th),'Color','r','HorizontalAlignment',...
%   'left','Units','normalized');
% text(ax,0.7,0.56,sprintf('%.2f',fracdif),'Color','k','HorizontalAlignment','left','Units','normalized');
text(ax,0.5,0.94,supertstr4th,'HorizontalAlignment','left','Units','normalized',...
  'fontsize',10);
text(ax,0.98,0.94,'b','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w','HorizontalAlignment','right');
% legend(ax,p,label,'Location','east');
% ylabel(ax,'Normalized count','FontSize',10);
xlabel(ax,'Inter-event time (s)','FontSize',10);
% xlabel(ax,sprintf('Diff. arrival between sources N and N-%d (s)',m));
xlim(ax,xran);
% if nsep == 1
%   yran = [0 0.2];
% elseif nsep == 2
%   yran = [0 0.06];
% else
%   yran = [0 0.04];
% end
ylim(ax,yran);
longticks(ax,2);
hold(ax,'off');

% if savefig
%   if rccwin == 0.25
%     fname = sprintf('interevttime%.2fs.pdf',rccwin);
%   else
%     fname = sprintf('interevttime.pdf');
%   end
%   print(f.fig,'-dpdf',...
%     strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
% end

% keyboard


%% plot the probability of inter-event time > t separately
widin = 6.5;  % maximum width allowed is 8.5 inches
htin = 3.5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f = initfig(widin,htin,nrow,ncol); %initialize fig
pltxran = [0.08 0.98]; pltyran = [0.12 0.98]; % optimal axis location
pltxsep = 0.07; pltysep = 0.05;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

xran = [0 50];
yran = [1e-4 1];
% xran = [0 7];
% yran = [5e-4 1];

ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); set(ax,'YScale','log');
xlim(ax,xran);
ylim(ax,yran);

dtinterplt = dtintern{1};
aa=sort(dtinterplt(:,1)/sps);
bb = prob_geq(aa);  %probability of inter-event time > t0
p(1)=plot(ax,bb(:,1),bb(:,2),'r-','LineWidth',1.5);
label{1}='Synthetic noise';
fitxran = [2 xran(2)]; 
ind = find(bb(:,1)>=fitxran(1) & bb(:,1)<=fitxran(2));
fitstruct=robustexpfit(bb(ind,1),bb(ind,2),'log10',[1e3 -1]);
fitobj=fitstruct.fitobj;
coef=coeffvalues(fitobj); slpn=coef(2); intcptn=coef(1);
xfit = xran(1):1e-2:xran(2);
yfitn = feval(fitobj,xfit);
plot(ax,xfit,yfitn,'r--','LineWidth',1);
text(ax,0.58,0.62,sprintf('slope: %.3f;\nintercept: %.3f',slpn,intcptn),'Units','normalized','FontSize',9,...
  'Color','r');

dtinterplt = dtintern{2};
aa=sort(dtinterplt(:,1)/sps);
bb = prob_geq(aa);  %probability of inter-event time > t0
p(2)=plot(ax,bb(:,1),bb(:,2),'b-','LineWidth',1.5);
label{2}='Noise not aligned';
ind = find(bb(:,1)>=fitxran(1) & bb(:,1)<=fitxran(2));
fitstruct=robustexpfit(bb(ind,1),bb(ind,2),'log10',[1e3 -1]);
fitobj=fitstruct.fitobj;
coef=coeffvalues(fitobj); slpn=coef(2); intcptn=coef(1);
xfit = xran(1):1e-2:xran(2);
yfitn = feval(fitobj,xfit);
plot(ax,xfit,yfitn,'b--','LineWidth',1);
text(ax,0.58,0.8,sprintf('slope: %.3f;\nintercept: %.3f',slpn,intcptn),'Units','normalized','FontSize',9,...
  'Color','b');

text(ax,0.58,0.5,sprintf('fit w/i [%.1f, %.1f] s',fitxran(1),fitxran(2)),'Units',...
  'normalized','FontSize',9);

plot(ax,ax.XLim,[1 1]/length(dtintern{1}),'k:');

text(ax,0.58,0.94,supertstr,'HorizontalAlignment','left','Units','normalized',...
  'fontsize',10);
text(ax,0.98,0.94,'a','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w','HorizontalAlignment','right');
xlabel(ax,'Time t (s)','FontSize',10);
ylabel(ax,sprintf('Prob(inter-event time \x2265 t)'),'FontSize',10);

ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); set(ax,'YScale','log');
xlim(ax,xran);
ylim(ax,yran);

dtinterplt = dtintern4th{1};
aa=sort(dtinterplt(:,1)/sps);
bb = prob_geq(aa);  %probability of inter-event time > t0
plot(ax,bb(:,1),bb(:,2),'r-','LineWidth',1.5);
ind = find(bb(:,1)>=fitxran(1) & bb(:,1)<=fitxran(2));
fitstruct=robustexpfit(bb(ind,1),bb(ind,2),'log10',[1e3 -1]);
fitobj=fitstruct.fitobj;
coef=coeffvalues(fitobj); slpn4th=coef(2); intcptn4th=coef(1);
yfitn = feval(fitobj,xfit);
plot(ax,xfit,yfitn,'r--','LineWidth',1);
text(ax,0.02,0.17,sprintf('slope: %.3f;\nintercept: %.3f',slpn4th,intcptn4th),'Units','normalized','FontSize',9,...
  'Color','r');

dtinterplt = dtintern4th{2};
aa=sort(dtinterplt(:,1)/sps);
bb = prob_geq(aa);  %probability of inter-event time > t0
plot(ax,bb(:,1),bb(:,2),'b-','LineWidth',1.5);
ind = find(bb(:,1)>=fitxran(1) & bb(:,1)<=fitxran(2));
fitstruct=robustexpfit(bb(ind,1),bb(ind,2),'log10',[1e3 -1]);
fitobj=fitstruct.fitobj;
coef=coeffvalues(fitobj); slpn4th=coef(2); intcptn4th=coef(1);
yfitn = feval(fitobj,xfit);
plot(ax,xfit,yfitn,'b--','LineWidth',1);
text(ax,0.02,0.3,sprintf('slope: %.3f;\nintercept: %.3f',slpn4th,intcptn4th),'Units','normalized','FontSize',9,...
  'Color','b');

text(ax,0.02,0.05,sprintf('fit w/i [%.1f, %.1f] s',fitxran(1),fitxran(2)),'Units',...
  'normalized','FontSize',9);

plot(ax,ax.XLim,[1 1]/length(dtintern{1}),'k:');

text(ax,0.58,0.94,supertstr4th,'HorizontalAlignment','left','Units','normalized',...
  'fontsize',10);
text(ax,0.98,0.94,'b','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w','HorizontalAlignment','right');
xlabel(ax,'Time t (s)','FontSize',10);
% ylabel(ax,sprintf('Prob(inter-event time \x2265 t)'),'FontSize',10);

% if savefig
%   if rccwin == 0.25
%     fname = sprintf('probintert%ds%.2fs.pdf',xran(2),rccwin);
%   else
%     fname = sprintf('probintert%ds.pdf',xran(2));
%   end
%   print(f.fig,'-dpdf',...
%     strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
% end

% keyboard