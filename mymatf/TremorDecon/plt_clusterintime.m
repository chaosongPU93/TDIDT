% plt_clusterintime.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the script in particular to plot a summary figure for the clustering
% of LFEs in time, including the inter-event time (N to N-1) in linear scale;
% in log scale for data and noise, 3-sta and 4-sta cats, and finally the
% fraction of events in clusters. In total, 6 panels.
% A combination of 'plt_difftime_NNm', '' 
% and 'plt_frac_uniqevt_incluster'
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
% savefile = 'deconv_stats4th_no23_allbstsig0.25s.mat';
load(strcat(rstpath, '/MAPS/',savefile));
% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
% savefile = 'deconv_stats4th_no23_allbstnoi0.25s.mat';
load(strcat(rstpath, '/MAPS/',savefile));
 

%%
%%%param for secondary sources removed
nsrc = allbstsig.nsrc;
imp = allbstsig.impindepall;
nsrcn = allbstnoi.nsrc;
impn = allbstnoi.impindepall;
supertstr = '3-station';

%%%param for further checked at KLNB
nsrc4th = allbstsig.nsrc4th;
imp4th = allbstsig.impindep4thall;
nsrcn4th = allbstnoi.nsrc4th;
impn4th = allbstnoi.impindep4thall;
supertstr4th = '4-station';


%%
widin = 8.4;  % maximum width allowed is 8.5 inches
htin = 5.5;   % maximum height allowed is 11 inches
nrow = 2;
ncol = 3;
f = initfig(widin,htin,nrow,ncol); %initialize fig
pltxran = [0.07 0.99]; pltyran = [0.08 0.98]; % optimal axis location
pltxsep = 0.07; pltysep = 0.05;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);


%%%inter-event time, ie. time from each to its preceeding event, N to N-1
m = 1;
nbst = size(trange,1);
dtcut = 0.25*m+0.125;
[amp,dtinter]=med_amp_incluster(nbst,imp,nsrc,m);
[ampn,dtintern]=med_amp_incluster(nbst,impn,nsrcn,m);
[amp4th,dtinter4th]=med_amp_incluster(nbst,imp4th,nsrc4th,m);
[ampn4th,dtintern4th]=med_amp_incluster(nbst,impn4th,nsrcn4th,m);

%% histogram of inter-event time, linear scale
xran = [0 2];
yran = [0 0.15];

binwdt = 0.05;
nx = round(xran(2)/binwdt)+1;
binedge = (xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2)';
bincnt = (xran(1): binwdt: xran(2))';

ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
patnoi = [0 yran(2);
          dtcut yran(2);
          dtcut yran(1);
          0 yran(1);
          0 yran(2)];
patch(ax,patnoi(:,1),patnoi(:,2),'k','Facealpha',0.15,'edgecolor','none');
Nd=histcounts(dtinter(:,1)/sps,binedge,'normalization','count');
Ndn = reshape(Nd,[],1)./length(dtinter);
frac = sum(dtinter(:,1)/sps<=dtcut)/length(dtinter);
Nn=histcounts(dtintern(:,1)/sps,binedge,'normalization','count');
Nnn = reshape(Nn,[],1)./length(dtinter);
fracn = sum(dtintern(:,1)/sps<=dtcut)/length(dtintern);
p(1)=plot(ax,bincnt,Ndn,'color','k','LineWidth',1.5);
label{1}='Data';
p(2)=plot(ax,bincnt,Nnn,'color','r','LineWidth',1.5);
label{2}='Synthetic noise';
% fracdif = (sum(dtarvl/sps<=dtcut)-sum(dtarvln/sps<=dtcut)) / ...
%   (length(dtarvl)-length(dtarvln));
% p(3)=plot(ax,cnt,Nn-Nnn,'color','k','LineWidth',1.5);
% label{3}='Data - Synthetic noise';
text(ax,0.25,0.75,sprintf('%.2f',frac),'Color','k','HorizontalAlignment',...
  'left','Units','normalized');
text(ax,0.25,0.65,sprintf('%.2f',fracn),'Color','r','HorizontalAlignment',...
  'left','Units','normalized');
% text(ax,0.7,0.56,sprintf('%.2f',fracdif),'Color','k','HorizontalAlignment','left','Units','normalized');
text(ax,0.5,0.94,supertstr,'HorizontalAlignment','left','Units','normalized',...
  'fontsize',10);
text(ax,0.98,0.94,'a','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w','HorizontalAlignment','right');
legend(ax,p,label,'Location','east');
% ylabel(ax,'Normalized count');
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

ax=f.ax(4); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
patnoi = [0 yran(2);
          dtcut yran(2);
          dtcut yran(1);
          0 yran(1);
          0 yran(2)];
patch(ax,patnoi(:,1),patnoi(:,2),'k','Facealpha',0.15,'edgecolor','none');
Nd4th=histcounts(dtinter4th(:,1)/sps,binedge,'normalization','count');
Nd4thn = reshape(Nd4th,[],1)./length(dtinter4th);
frac4th = sum(dtinter4th(:,1)/sps<=dtcut)/length(dtinter4th);
Nn4th=histcounts(dtintern4th(:,1)/sps,binedge,'normalization','count');
Nn4thn = reshape(Nn4th,[],1)./length(dtinter4th);
fracn4th = sum(dtintern4th(:,1)/sps<=dtcut)/length(dtintern4th);
p(1)=plot(ax,bincnt,Nd4thn,'color','k','LineWidth',1.5);
label{1}='Data';
p(2)=plot(ax,bincnt,Nn4thn,'color','r','LineWidth',1.5);
label{2}='Synthetic noise';
% fracdif = (sum(dtarvl/sps<=dtcut)-sum(dtarvln/sps<=dtcut)) / ...
%   (length(dtarvl)-length(dtarvln));
% p(3)=plot(ax,cnt,Nn-Nnn,'color','k','LineWidth',1.5);
% label{3}='Data - Synthetic noise';
text(ax,0.25,0.75,sprintf('%.2f',frac4th),'Color','k','HorizontalAlignment',...
  'left','Units','normalized');
text(ax,0.25,0.65,sprintf('%.2f',fracn4th),'Color','r','HorizontalAlignment',...
  'left','Units','normalized');
% text(ax,0.7,0.56,sprintf('%.2f',fracdif),'Color','k','HorizontalAlignment','left','Units','normalized');
text(ax,0.5,0.94,supertstr4th,'HorizontalAlignment','left','Units','normalized',...
  'fontsize',10);
text(ax,0.98,0.94,'d','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w','HorizontalAlignment','right');
legend(ax,p,label,'Location','east');
ylabel(ax,'Normalized count');
xlabel(ax,'Inter-event time dt (s)');
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

%% histogram of inter-event time, specific bins, log scale
% xran = [0 7];
% yran = [1e-4 1];
% binedge=[0 0.375:0.25:ceil(max(dtinter(:,1))/sps)]';  % 1st bin [0 0.375], then 0.25 increment
% bincnt=mean([binedge(1:end-1) binedge(2:end)],2);
% 
% ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); set(ax,'YScale','log');
% [Nd]=histcounts(dtinter(:,1)/sps,binedge,'normalization','count');
% Ndn = reshape(Nd,[],1)./length(dtinter);
% [Nn]=histcounts(dtintern(:,1)/sps,binedge,'normalization','count');
% Nnn = reshape(Nn,[],1)./length(dtinter);
% p(1)=plot(ax,bincnt,Ndn,'k-','LineWidth',1.5,'marker','o','markersize',3,...
%   'markerfacec','k');
% label{1}='Data';
% p(2)=plot(ax,bincnt,Nnn,'r-','LineWidth',1.5);%,'marker','o','markersize',3,'markerfacec','r'
% label{2}='Synthetic noise';
% fitxran = [2 xran(2)]; 
% ind = find(bincnt>=fitxran(1) & bincnt<=fitxran(2));
% fitstruct=robustexpfit(bincnt(ind),Ndn(ind),'log10',[1e3 -1]);
% fitobj=fitstruct.fitobj;
% coef=coeffvalues(fitobj); slpd=coef(2);
% yfitd = feval(fitobj,bincnt);
% plot(ax,bincnt,yfitd,'k--','LineWidth',1);
% 
% fitstruct=robustexpfit(bincnt(ind),Nnn(ind),'log10',[1e3 -1]);
% fitobj=fitstruct.fitobj;
% coef=coeffvalues(fitobj); slpn=coef(2);
% yfitn = feval(fitobj,bincnt);
% plot(ax,bincnt,yfitn,'r--','LineWidth',1);
% % xlabel(ax,'Time (s) from each to its preceding');
% % xlabel(ax,'Time delay (s)');
% % ylabel(ax,'Normalized count');
% xlim(ax,xran);
% ylim(ax,yran);
% text(ax,0.02,0.1,sprintf('fit w/i [%.1f, %.1f] s',fitxran(1),fitxran(2)),'Units',...
%   'normalized','FontSize',9);
% % text(ax,0.02,0.25,sprintf('From each to its \npreceding'),'Units','normalized',...
% %   'FontSize',9);
% text(ax,0.5,0.94,supertstr,'HorizontalAlignment','left','Units','normalized',...
%   'fontsize',10);
% text(ax,0.98,0.94,'b','FontSize',10,'unit','normalized','EdgeColor','k',...
%   'Margin',1,'backgroundcolor','w','HorizontalAlignment','right');
% % [lambdahat,lambdaci] = poissfit(dtinter/sps)
% % x=0:1:12;
% % y=poisspdf(x,lambdahat);
% % plot(ax,x,y,'r-','linew',2);
% % keyboard
% longticks(ax,2);
% hold(ax,'off');
% 
% 
% ax=f.ax(5); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); set(ax,'YScale','log');
% [Nd4th]=histcounts(dtinter4th(:,1)/sps,binedge,'normalization','count');
% Nd4thn = reshape(Nd4th,[],1)./length(dtinter4th);
% [Nn4th]=histcounts(dtintern4th(:,1)/sps,binedge,'normalization','count');
% Nn4thn = reshape(Nn4th,[],1)./length(dtinter4th);
% p(1)=plot(ax,bincnt,Nd4thn,'k-','LineWidth',1.5,'marker','o','markersize',3,...
%   'markerfacec','k');
% label{1}='Data';
% p(2)=plot(ax,bincnt,Nn4thn,'r-','LineWidth',1.5);%,'marker','o','markersize',3,'markerfacec','r'
% label{2}='Synthetic noise';
% fitxran = [2 xran(2)]; 
% ind = find(bincnt>=fitxran(1) & bincnt<=fitxran(2));
% fitstruct=robustexpfit(bincnt(ind),Nd4thn(ind),'log10',[1e3 -1]);
% fitobj=fitstruct.fitobj;
% coef=coeffvalues(fitobj); slpd=coef(2);
% yfitd4th = feval(fitobj,bincnt);
% plot(ax,bincnt,yfitd4th,'k--','LineWidth',1);
% 
% fitstruct=robustexpfit(bincnt(ind),Nn4thn(ind),'log10',[1e3 -1]);
% fitobj=fitstruct.fitobj;
% coef=coeffvalues(fitobj); slpn=coef(2);
% yfitn4th = feval(fitobj,bincnt);
% plot(ax,bincnt,yfitn4th,'r--','LineWidth',1);
% % xlabel(ax,'Time (s) from each to its preceding');
% % xlabel(ax,'Time delay (s)');
% xlabel(ax,'Inter-event time (s)');
% ylabel(ax,'Normalized count');
% xlim(ax,xran);
% ylim(ax,yran);
% % legend(ax,p,'Data','Synthetic noise');
% text(ax,0.02,0.1,sprintf('fit w/i [%.1f, %.1f] s',fitxran(1),fitxran(2)),'Units',...
%   'normalized','FontSize',9);
% % text(ax,0.02,0.25,sprintf('From each to its \npreceding'),'Units','normalized',...
% %   'FontSize',9);
% text(ax,0.5,0.94,supertstr4th,'HorizontalAlignment','left','Units','normalized',...
%   'fontsize',10);
% text(ax,0.98,0.94,'e','FontSize',10,'unit','normalized','EdgeColor','k',...
%   'Margin',1,'backgroundcolor','w','HorizontalAlignment','right');
% % [lambdahat,lambdaci] = poissfit(dtinter/sps)
% % x=0:1:12;
% % y=poisspdf(x,lambdahat);
% % plot(ax,x,y,'r-','linew',2);
% % keyboard
% longticks(ax,2);
% hold(ax,'off');

%% probability of inter-event time > t0, log scale
xran = [0 20];
yran = [1e-4 1];

ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); set(ax,'YScale','log');
xlim(ax,xran);
ylim(ax,yran);
aa=sort(dtinter(:,1)/sps);
% bb=(length(aa):-1:1)/length(aa);
[unival, ~, ind] = unique(aa);
counts = accumarray(ind, 1);
countsleq=zeros(length(counts),1);
for i=1:length(counts)
  countsleq(i,1)=sum(counts(i:end));
end
% Combine the unique values and their counts into a two-column matrix
bb = [unival countsleq/length(aa)];

p(1)=plot(ax,bb(:,1),bb(:,2),'k-','LineWidth',1.5);
label{1}='Data';
fitxran = [2 xran(2)]; 
ind = find(bb(:,1)>=fitxran(1) & bb(:,1)<=fitxran(2));
fitstruct=robustexpfit(bb(ind,1),bb(ind,2),'log10',[1e3 -1]);
fitobj=fitstruct.fitobj;
coef=coeffvalues(fitobj); slpd=coef(2);
xfit = xran(1):1e-2:xran(2);
yfitd = feval(fitobj,xfit);
plot(ax,xfit,yfitd,'k--','LineWidth',1);

aa=sort(dtintern(:,1)/sps);
% bb=(length(aa):-1:1)/length(aa);
[unival, ~, ind] = unique(aa);
counts = accumarray(ind, 1);
countsleq=zeros(length(counts),1);
for i=1:length(counts)
  countsleq(i,1)=sum(counts(i:end));
end
% Combine the unique values and their counts into a two-column matrix
bb = [unival countsleq/length(aa)];

p(2)=plot(ax,bb(:,1),bb(:,2),'r-','LineWidth',1.5);
label{2}='Synthetic noise';
ind = find(bb(:,1)>=fitxran(1) & bb(:,1)<=fitxran(2));
fitstruct=robustexpfit(bb(ind,1),bb(ind,2),'log10',[1e3 -1]);
fitobj=fitstruct.fitobj;
coef=coeffvalues(fitobj); slpn=coef(2);
yfitn = feval(fitobj,xfit);
plot(ax,xfit,yfitn,'r--','LineWidth',1);

text(ax,0.45,0.75,sprintf('slope: %.3f',slpd),'Units','normalized','FontSize',9,...
  'Color','k');
text(ax,0.45,0.65,sprintf('slope: %.3f',slpn),'Units','normalized','FontSize',9,...
  'Color','r');
text(ax,0.45,0.55,sprintf('fit w/i [%.1f, %.1f] s',fitxran(1),fitxran(2)),'Units',...
  'normalized','FontSize',9);

plot(ax,ax.XLim,[1 1]/length(dtinter(:,1)),'k:');

text(ax,0.5,0.94,supertstr,'HorizontalAlignment','left','Units','normalized',...
  'fontsize',10);
text(ax,0.98,0.94,'b','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w','HorizontalAlignment','right');
% xlabel(ax,'Inter-event time dt (s)');
% ylabel(ax,'Prob(dt>=t)');

ax=f.ax(5); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); set(ax,'YScale','log');
xlim(ax,xran);
ylim(ax,yran);
aa=sort(dtinter4th(:,1)/sps);
% bb=(length(aa):-1:1)/length(aa);
[unival, ~, ind] = unique(aa);
counts = accumarray(ind, 1);
countsleq=zeros(length(counts),1);
for i=1:length(counts)
  countsleq(i,1)=sum(counts(i:end));
end
% Combine the unique values and their counts into a two-column matrix
bb = [unival countsleq/length(aa)];
plot(ax,bb(:,1),bb(:,2),'k-','LineWidth',1.5);
ind = find(bb(:,1)>=fitxran(1) & bb(:,1)<=fitxran(2));
fitstruct=robustexpfit(bb(ind,1),bb(ind,2),'log10',[1e3 -1]);
fitobj=fitstruct.fitobj;
coef=coeffvalues(fitobj); slpd=coef(2);
yfitd = feval(fitobj,xfit);
plot(ax,xfit,yfitd,'k--','LineWidth',1);

aa=sort(dtintern4th(:,1)/sps);
% bb=(length(aa):-1:1)/length(aa);
[unival, ~, ind] = unique(aa);
counts = accumarray(ind, 1);
countsleq=zeros(length(counts),1);
for i=1:length(counts)
  countsleq(i,1)=sum(counts(i:end));
end
% Combine the unique values and their counts into a two-column matrix
bb = [unival countsleq/length(aa)];
plot(ax,bb(:,1),bb(:,2),'r-','LineWidth',1.5);
ind = find(bb(:,1)>=fitxran(1) & bb(:,1)<=fitxran(2));
fitstruct=robustexpfit(bb(ind,1),bb(ind,2),'log10',[1e3 -1]);
fitobj=fitstruct.fitobj;
coef=coeffvalues(fitobj); slpn=coef(2);
yfitn = feval(fitobj,xfit);
plot(ax,xfit,yfitn,'r--','LineWidth',1);

text(ax,0.02,0.3,sprintf('slope: %.3f',slpd),'Units','normalized','FontSize',9,...
  'Color','k');
text(ax,0.02,0.2,sprintf('slope: %.3f',slpn),'Units','normalized','FontSize',9,...
  'Color','r');
text(ax,0.02,0.1,sprintf('fit w/i [%.1f, %.1f] s',fitxran(1),fitxran(2)),'Units',...
  'normalized','FontSize',9);

plot(ax,ax.XLim,[1 1]/length(dtinter(:,1)),'k:');
text(ax,0.5,0.94,supertstr4th,'HorizontalAlignment','left','Units','normalized',...
  'fontsize',10);
text(ax,0.98,0.94,'e','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w','HorizontalAlignment','right');
xlabel(ax,'Inter-event time dt (s)');
% ylabel(ax,sprintf('Prob(dt%ct)',char(242)));
ylabel(ax,sprintf('Prob(dt \x2265 t)'));
% ylabel(ax,'Prob(dt \x2265 t)');


%% fraction of catalog of events in clusters vs. # of events in clusters
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


ax = f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
p(1)=plot(ax,1:mmax+1,fracuimpsum,'k-','Linew',1.5,'marker','o','markersize',3,...
  'markerfacec','k');
p(2)=plot(ax,1:mmaxn+1,fracuimpnsum,'r-','Linew',1.5,'marker','o','markersize',3,...
  'markerfacec','r');
% legend(ax,p,'Data','Noise','Location','east');
ax.YAxis.Scale = 'log'; %make y axis log10 scale
xlim(ax,[1 mmax+2]);
ylim(ax,[5e-2 1e2]);
xticks(ax,0:2:mmax+2);
yticks(ax,[0.1 0.2 0.5 1 2 5 10 20 50 100]);
% xlabel(ax,'# of events in cluster','FontSize',10);
% ylabel(ax,'% of catalog in such clusters','FontSize',10);
text(ax,0.5,0.94,supertstr,'HorizontalAlignment','left','Units','normalized',...
  'fontsize',10);
text(ax,0.98,0.94,'c','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w','HorizontalAlignment','right');

ax = f.ax(6); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
p(1)=plot(ax,1:mmax4th+1,fracuimp4thsum,'k-','Linew',1.5,'marker','o',...
  'markersize',3,'markerfacec','k');
p(2)=plot(ax,1:mmaxn4th+1,fracuimpn4thsum,'r-','Linew',1.5,'marker','o',...
  'markersize',3,'markerfacec','r');
ax.YAxis.Scale = 'log'; %make y axis log10 scale
xlim(ax,[1 mmax+2]);
ylim(ax,[5e-2 1e2]);
xticks(ax,0:2:mmax+2);
yticks(ax,[0.1 0.2 0.5 1 2 5 10 20 50 100]);
xlabel(ax,'m (# of events in cluster)','FontSize',10);
ylabel(ax,sprintf('%% of catalog in clusters of \x2265 m events'),'FontSize',10);
% ylabel(ax,'% of catalog in clusters of >=m events','FontSize',10);
text(ax,0.5,0.94,supertstr4th,'HorizontalAlignment','left','Units','normalized',...
  'fontsize',10);
text(ax,0.98,0.94,'f','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w','HorizontalAlignment','right');

fname = strcat('tclusteringnn1.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));


%%
% widin = 8.4;  % maximum width allowed is 8.5 inches
% htin = 4.5;   % maximum height allowed is 11 inches
% nrow = 1;
% ncol = 2;
% f = initfig(widin,htin,nrow,ncol); %initialize fig
% pltxran = [0.08 0.99]; pltyran = [0.1 0.98]; % optimal axis location
% pltxsep = 0.08; pltysep = 0.05;
% optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
% 
% xran = [0 20];
% yran = [1e-4 1];
% % xran = [0 7];
% % yran = [1e-3 1];
% 
% ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); set(ax,'YScale','log');
% xlim(ax,xran);
% ylim(ax,yran);
% aa=sort(dtinter(:,1)/sps);
% % bb=(length(aa):-1:1)/length(aa);
% [unival, ~, ind] = unique(aa);
% counts = accumarray(ind, 1);
% countsleq=zeros(length(counts),1);
% for i=1:length(counts)
%   countsleq(i,1)=sum(counts(i:end));
% end
% % Combine the unique values and their counts into a two-column matrix
% bb = [unival countsleq/length(aa)];
% 
% p(1)=plot(ax,bb(:,1),bb(:,2),'k-','LineWidth',1.5);
% label{1}='Data';
% fitxran = [2 xran(2)]; 
% ind = find(bb(:,1)>=fitxran(1) & bb(:,1)<=fitxran(2));
% fitstruct=robustexpfit(bb(ind,1),bb(ind,2),'log10',[1e3 -1]);
% fitobj=fitstruct.fitobj;
% coef=coeffvalues(fitobj); slpd=coef(2);
% xfit = xran(1):1e-2:xran(2);
% yfitd = feval(fitobj,xfit);
% plot(ax,xfit,yfitd,'k--','LineWidth',1);
% 
% aa=sort(dtintern(:,1)/sps);
% % bb=(length(aa):-1:1)/length(aa);
% [unival, ~, ind] = unique(aa);
% counts = accumarray(ind, 1);
% countsleq=zeros(length(counts),1);
% for i=1:length(counts)
%   countsleq(i,1)=sum(counts(i:end));
% end
% % Combine the unique values and their counts into a two-column matrix
% bb = [unival countsleq/length(aa)];
% 
% p(2)=plot(ax,bb(:,1),bb(:,2),'r-','LineWidth',1.5);
% label{2}='Synthetic noise';
% ind = find(bb(:,1)>=fitxran(1) & bb(:,1)<=fitxran(2));
% fitstruct=robustexpfit(bb(ind,1),bb(ind,2),'log10',[1e3 -1]);
% fitobj=fitstruct.fitobj;
% coef=coeffvalues(fitobj); slpn=coef(2);
% yfitn = feval(fitobj,xfit);
% plot(ax,xfit,yfitn,'r--','LineWidth',1);
% 
% text(ax,0.02,0.3,sprintf('slope: %.3f',slpd),'Units','normalized','FontSize',9,...
%   'Color','k');
% text(ax,0.02,0.2,sprintf('slope: %.3f',slpn),'Units','normalized','FontSize',9,...
%   'Color','r');
% text(ax,0.02,0.1,sprintf('fit w/i [%.1f, %.1f] s',fitxran(1),fitxran(2)),'Units',...
%   'normalized','FontSize',9);
% 
% plot(ax,ax.XLim,[1 1]/length(dtinter(:,1)),'k:');
% 
% text(ax,0.5,0.94,supertstr,'HorizontalAlignment','left','Units','normalized',...
%   'fontsize',10);
% xlabel(ax,'Inter-event time dt (s)');
% ylabel(ax,'P(dt>=t)');
% 
% ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); set(ax,'YScale','log');
% xlim(ax,xran);
% ylim(ax,yran);
% aa=sort(dtinter4th(:,1)/sps);
% % bb=(length(aa):-1:1)/length(aa);
% [unival, ~, ind] = unique(aa);
% counts = accumarray(ind, 1);
% countsleq=zeros(length(counts),1);
% for i=1:length(counts)
%   countsleq(i,1)=sum(counts(i:end));
% end
% % Combine the unique values and their counts into a two-column matrix
% bb = [unival countsleq/length(aa)];
% plot(ax,bb(:,1),bb(:,2),'k-','LineWidth',1.5);
% ind = find(bb(:,1)>=fitxran(1) & bb(:,1)<=fitxran(2));
% fitstruct=robustexpfit(bb(ind,1),bb(ind,2),'log10',[1e3 -1]);
% fitobj=fitstruct.fitobj;
% coef=coeffvalues(fitobj); slpd=coef(2);
% yfitd = feval(fitobj,xfit);
% plot(ax,xfit,yfitd,'k--','LineWidth',1);
% 
% aa=sort(dtintern4th(:,1)/sps);
% % bb=(length(aa):-1:1)/length(aa);
% [unival, ~, ind] = unique(aa);
% counts = accumarray(ind, 1);
% countsleq=zeros(length(counts),1);
% for i=1:length(counts)
%   countsleq(i,1)=sum(counts(i:end));
% end
% % Combine the unique values and their counts into a two-column matrix
% bb = [unival countsleq/length(aa)];
% plot(ax,bb(:,1),bb(:,2),'r-','LineWidth',1.5);
% ind = find(bb(:,1)>=fitxran(1) & bb(:,1)<=fitxran(2));
% fitstruct=robustexpfit(bb(ind,1),bb(ind,2),'log10',[1e3 -1]);
% fitobj=fitstruct.fitobj;
% coef=coeffvalues(fitobj); slpn=coef(2);
% yfitn = feval(fitobj,xfit);
% plot(ax,xfit,yfitn,'r--','LineWidth',1);
% 
% text(ax,0.02,0.3,sprintf('slope: %.3f',slpd),'Units','normalized','FontSize',9,...
%   'Color','k');
% text(ax,0.02,0.2,sprintf('slope: %.3f',slpn),'Units','normalized','FontSize',9,...
%   'Color','r');
% text(ax,0.02,0.1,sprintf('fit w/i [%.1f, %.1f] s',fitxran(1),fitxran(2)),'Units',...
%   'normalized','FontSize',9);
% 
% plot(ax,ax.XLim,[1 1]/length(dtinter(:,1)),'k:');
% text(ax,0.5,0.94,supertstr4th,'HorizontalAlignment','left','Units','normalized',...
%   'fontsize',10);
% xlabel(ax,'Inter-event time dt (s)');
% ylabel(ax,'P(dt>=t)');
% 
% fname = sprintf('probintert%ds.pdf',xran(2));
% print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));
% 
% keyboard













