% circsrcmodel.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the script to assume a circular source model. Now we have locations
% of all deconvolved LFE sources where each source is represented by a circle
% with a series of diameter. This script aims to show several things:
% 1. for any cluster of m+1 consecutive events, what is the distribution
% 2. To summarize all these images, consider the ratio of overlapping area.
% with a certain diameter, sources may or may not overlap. in total, they have
% a covered area. Meanwhile, there is a nominal area which is the sum of the
% area of each circle. 1 minus this ratio of area indicates the overlapping 
% extent.
% 
% 2024/02/06, mapping from samples to relation location already contains error,
% so when you try to obtain the location difference, or distance, start with 
% the difference in samples first. For example, if you ask the difference in 
% certain directions, first obtain the difference in off12 and off13, then map
% the diff to location, then you can project along any direction to get the
% loc difference along that direction, or the distance (abs of loc diff) in
% that direction, etc.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/01/04
% Last modified date:   2024/01/04
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
[imploc0, ~] = off2space002([0 0],sps,ftrans,0);  % a ref source at 0,0

%%%load data
savefile = 'deconv_stats4th_allbstsig.mat';
load(strcat(rstpath, '/MAPS/',savefile));
savefile = 'deconv_stats4th_allbstnoi.mat';
load(strcat(rstpath, '/MAPS/',savefile));
 
% keyboard

%%
% %%%param for secondary sources removed
% locxyprojall = allbstsig.locxyprojall;
% tarvlsplstall = allbstsig.impindepall(:,1);
% nsrc = allbstsig.nsrc;
% imp = allbstsig.impindepall;
% locxyprojalln = allbstnoi.locxyprojall;
% tarvlsplstalln = allbstnoi.impindepall(:,1);
% nsrcn = allbstnoi.nsrc;
% impn = allbstnoi.impindepall;
% supertstr = 'Secondary sources removed';
% fnsuffix = [];

%%%param for further checked at KLNB
locxyprojall = allbstsig.locxyproj4thall;
tarvlsplstall = allbstsig.impindep4thall(:,1);
nsrc = allbstsig.nsrc4th;
imp = allbstsig.impindep4thall;
locxyprojalln = allbstnoi.locxyproj4thall;
tarvlsplstalln = allbstnoi.impindep4thall(:,1);
nsrcn = allbstnoi.nsrc4th;
impn = allbstnoi.impindep4thall;
supertstr = 'Further checked at KLNB';
fnsuffix = '4th';

%choice of diameter of circular sources
diams = [0.1 0.3 0.5 0.7];

% impuse = imp;
% nsrcuse = nsrc;
impuse = impn;
nsrcuse = nsrcn;

%to the 'mmax' for which the resulting number of clusters is nonzero
mmax=1;
timetype = 'tarvl';
dtcut = 0.25*mmax+0.125;
[~,impuni]=evtcluster(nbst,impuse,nsrcuse,mmax,dtcut,sps,timetype);
while ~isempty(impuni)
  mmax=mmax+1;
  dtcut = 0.25*mmax+0.125;
  [~,impuni]=evtcluster(nbst,impuse,nsrcuse,mmax,dtcut,sps,timetype);
end
mmax=mmax-1;

%% loc difference between events pairs N and N-n in clusters of m+1 events (N and N-m w/i dtcut)
% mmax=14; %to define cluster of events within N and N-m
n=1;  %to decide which event pair inthe cluster to look at
dtlump=cell(mmax,1);
dloclump=cell(mmax,1);
dlocspllump=cell(mmax,1);
dprojloclump=cell(mmax,1);
dprojloccslump=cell(mmax,1);
% dprojlocspllump=cell(mmax,1);
% dprojlocsplcslump=cell(mmax,1);
for m=n:mmax
  % m=5;  %to define cluster of events within N and N-m
  timetype = 'tarvl';
  % timetype = 'tori';
  %obtain clusters of events defined by N and N-m w/i dtcut
  dtcut = 0.25*m+0.125;
  [impcluster,impuni]=evtcluster(nbst,impuse,nsrcuse,m,dtcut,sps,timetype);
  implocclus = off2space002(impcluster(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

  %for each cluster, obtain their diff location between events N and N-n
  % timetype = 'tori';
  [~,dt,k,dloc_spl]=dloc_evtcluster(impcluster,implocclus,sps,ftrans,m,n,timetype);
  temp = off2space002(dloc_spl(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
  dloc=temp(:,1:2);
  dtlump{m}=dt;
  dloclump{m}=dloc;
  dlocspllump{m}=dloc_spl;
  
%   %plot the diff location in map
%   [coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc_spl);
%   cstr={'# events / pixel'}; xran=[-50 50]; yran=[-50 50]; dx=1; dy=1;
%   [f,den1d_spl,conmat]=plt_srcdlocinmap(dt/sps,dloc_spl,[],'spl',timetype,...
%     6,cstr,'o','linear','pixel',xran,yran,dx,dy,1,2);
%   ax=f.ax(2);
%   hold(ax,'on');
%   plot(ax,ellx,elly,'k-','linew',2);
%   plot(ax,x0+10*[coeff(1,2),-coeff(1,2)],y0+10*[coeff(2,2),-coeff(2,2)],'k--','linew',2);
%   plot(ax,x0+10*[coeff(1,1),-coeff(1,1)],y0+10*[coeff(2,1),-coeff(2,1)],'k--','linew',2);
%   text(ax,0.98,0.15,sprintf('%d ^o',round(anglegeo(2))),'Units','normalized',...
%     'HorizontalAlignment','right','FontSize',9);
%   text(ax,0.98,0.1,sprintf('asprat: %.1f',semia/semib),'Units','normalized',...
%     'HorizontalAlignment','right','FontSize',9);
%   hold(ax,'off');
  
%   %%
%   [coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc); 
%   cstr={'# events / pixel'}; xran=[-4 4]; yran=[-4 4]; dx=0.1; dy=0.1;
%   [f,den1d]=plt_srcdlocinmap(dt/sps,dloc,[],'km',timetype,...
%     10,cstr,'o','linear','pixel',xran,yran,dx,dy,0);
%   ax=f.ax(2);
%   hold(ax,'on');
%   plot(ax,ellx,elly,'k-','linew',2);
%   plot(ax,x0+[coeff(1,2),-coeff(1,2)],y0+[coeff(2,2),-coeff(2,2)],'k--','linew',2);
%   plot(ax,x0+[coeff(1,1),-coeff(1,1)],y0+[coeff(2,1),-coeff(2,1)],'k--','linew',2);
%   text(ax,0.98,0.15,sprintf('%d ^o',round(anglegeo(2))),'Units','normalized',...
%     'HorizontalAlignment','right','FontSize',9);
%   text(ax,0.98,0.1,sprintf('asprat: %.1f',semia/semib),'Units','normalized',...
%     'HorizontalAlignment','right','FontSize',9);
%   hold(ax,'off');
%   keyboard
    
  %there is some elongation along N45E, then we want to project along that direction
  projang=135;
  % projang=ang2; %note that 'ang2', 2nd pc axis from PCA is very close to the eyeball value 135
  [~,~,dprojloc] = customprojection(dloc(:,1:2),projang);
%   f=plt_srcdlocinmap(dt/sps,dprojloc,[],'km',timetype,...
%     10,cstr,'o','linear','pixel',xran,yran,dx,dy,0);
  ndloc(m)=size(dprojloc,1);
  dprojloclump{m}=dprojloc;
  
  % %histogram of diff location in various ways
  % f = initfig(12,4.5,1,3);
  % xran=[-4 4]; binw=0.1; legendstr1={'E','N'}; normopt='countdensity';
  % f.ax(1)=plt_dlochist(f.ax(1),dloc,xran,binw,legendstr1,normopt);
  % legendstr2={'S45E','N45E'};
  % f.ax(2)=plt_dlochist(f.ax(2),dprojloc,xran,binw,legendstr2,normopt);
  % % [~,mu1(m),sigma1(m),mu2(m),sigma2(m)]=plt_dlochist([],dprojloc,xran,binw,legendstr2,normopt);
  %focus on a narrow strip along the projection direction
  dprojloccs=dprojloc(abs(dprojloc(:,2))<=0.2*sqrt(2), :);
  dprojloccslump{m}=dprojloccs;
  % legendstr3={'S45E-cs'};
  % f.ax(3)=plt_dlochist(f.ax(3),dprojloccs(:,1),xran,binw,legendstr3,normopt);

  % keyboard
end

%%
st=n;
dt=cat(1,dtlump{st:end});
dloc_spl=cat(1,dlocspllump{st:end});
dloc=cat(1,dloclump{st:end});
% dprojloc=cat(1,dprojloclump{st:end});
% dprojloccs=cat(1,dprojloccslump{st:end});

%%
%plot the diff location in map
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc_spl);
cstr={'# events / pixel'}; xran=[-50 50]; yran=[-50 50]; dx=1; dy=1;
[f,den1d_spl,conmat]=plt_srcdlocinmap(dt/sps,dloc_spl,[],'spl',timetype,...
  6,cstr,'o','linear','pixel',xran,yran,dx,dy,1,1);
ax=f.ax(2);
hold(ax,'on');
plot(ax,ellx,elly,'k-','linew',2);
plot(ax,x0+10*[coeff(1,2),-coeff(1,2)],y0+10*[coeff(2,2),-coeff(2,2)],'k--','linew',2);
plot(ax,x0+10*[coeff(1,1),-coeff(1,1)],y0+10*[coeff(2,1),-coeff(2,1)],'k--','linew',2);
text(ax,0.98,0.15,sprintf('%d ^o',round(anglegeo(2))),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.98,0.1,sprintf('asprat: %.1f',semia/semib),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
hold(ax,'off');
% keyboard

%%
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc);
cstr={'# events / pixel'}; xran=[-4 4]; yran=[-4 4]; dx=0.1; dy=0.1;
[f,den1d]=plt_srcdlocinmap(dt/sps,dloc,[],'km',timetype,...
  10,cstr,'o','linear','pixel',xran,yran,dx,dy,0);
ax=f.ax(2);
hold(ax,'on');
plot(ax,ellx,elly,'k-','linew',2);
plot(ax,x0+[coeff(1,2),-coeff(1,2)],y0+[coeff(2,2),-coeff(2,2)],'k--','linew',2);
plot(ax,x0+[coeff(1,1),-coeff(1,1)],y0+[coeff(2,1),-coeff(2,1)],'k--','linew',2);
text(ax,0.98,0.15,sprintf('%d ^o',round(anglegeo(2))),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.98,0.1,sprintf('asprat: %.1f',semia/semib),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
hold(ax,'off');

% %% if choose a strip with a finite width in off12 or off13.
% widspl=2;
% ind1=find(abs(dloc_spl(:,1))<=widspl);
% ind2=find(abs(dloc_spl(:,2))<=widspl);
% f = initfig(8,4.5,1,2);
% ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on'); axis(ax, 'equal');
% scatter(ax,dloc_spl(ind2,1),dloc_spl(ind2,2),3,[.8 .8 .8],'filled');
% scatter(ax,dloc_spl(ind1,1),dloc_spl(ind1,2),3,[.2 .2 .2],'filled');
% xlim(ax,[-50 50]); ylim(ax,[-50 50]);
% xlabel(ax,'Diff off12 (samples)'); ylabel(ax,'Diff off13 (samples)');
% ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on'); axis(ax, 'equal');
% scatter(ax,dloc(ind2,1),dloc(ind2,2),3,[.8 .8 .8],'filled');
% scatter(ax,dloc(ind1,1),dloc(ind1,2),3,[.2 .2 .2],'filled');
% xlim(ax,[-4 4]); ylim(ax,[-4 4]);
% xlabel(ax,'Diff E loc (km)'); ylabel(ax,'Diff N loc (km)');
% [coeff1,score,angle,anglegeo1,x0,y0]=pcaellipse(dloc(ind1,1:2));
% plot(ax,x0+10*[coeff1(1,1),-coeff1(1,1)],y0+10*[coeff1(2,1),-coeff1(2,1)],...
%   '-','linew',2,'Color',[.2 .2 .2]);
% text(ax,0.98,0.15,sprintf('%d ^o',round(anglegeo1(1))),'Units','normalized',...
%   'HorizontalAlignment','right','FontSize',9);
% [coeff2,score,angle,anglegeo2,x0,y0]=pcaellipse(dloc(ind2,1:2));
% plot(ax,x0+10*[coeff2(1,1),-coeff2(1,1)],y0+10*[coeff2(2,1),-coeff2(2,1)],...
%   '-','linew',2,'Color',[.8 .8 .8]);
% text(ax,0.98,0.1,sprintf('%d ^o',round(anglegeo2(1))),'Units','normalized',...
%   'HorizontalAlignment','right','FontSize',9);
% % keyboard
% 
% dlocsplplt=dloc_spl(ind1,:);
% dlocplt=dloc(ind1,:);
% [~,~,dprojlocplt] = customprojection(dlocplt(:,1:2),anglegeo1(1));
% f = initfig(8,4.5,1,2);
% supertit(f.ax,'abs(off12)<=2 samples',10);
% xran=[-40 40]; binw=1; legendstr1={'off13'}; normopt='count';
% f.ax(1)=plt_dlochist(f.ax(1),dlocsplplt(:,2),xran,binw,legendstr1,normopt);
% xlabel(f.ax(1),'Location difference (samples)');
% 
% f = initfig(8,4.5,1,2);
% supertit(f.ax,'abs(off12)<=2 samples',10);
% xran=[-4 4]; binw=0.1; legendstr1={'E','N'}; normopt='countdensity';
% f.ax(1)=plt_dlochist(f.ax(1),dlocplt,xran,binw,legendstr1,normopt);
% legendstr2={sprintf('%d ^o',round(anglegeo1(1)))};
% f.ax(2)=plt_dlochist(f.ax(2),dprojlocplt(:,1),xran,binw,legendstr2,normopt);  
% distprojlocplt=median(abs(dprojlocplt(:,1)));

%% if choose a strip with a finite width in the PCA direction
[~,~,dprojloc] = customprojection(dloc(:,1:2),anglegeo(2));

wid1=0.25;
cstr={'# events / pixel'}; xran=[-4 4]; yran=[-4 4]; dx=0.1; dy=0.1;
[f,den1d]=plt_srcdlocinmap(dt/sps,dprojloc,[],'km',timetype,...
  10,cstr,'o','linear','pixel',xran,yran,dx,dy,0);
ax=f.ax(2);
hold(ax,'on');
plot(ax,ax.XLim,[wid1 wid1],'-','linew',2,'color',[.7 .7 .7]);
plot(ax,ax.XLim,[-wid1 -wid1],'-','linew',2,'color',[.7 .7 .7]);
hold(ax,'off');

%%
f = initfig(8,8,2,2);
supertit(f.ax(1:2),sprintf('all, %d measurements',size(dloc,1)),10);
xran=[-4 4]; binw=0.1; legendstr1={'E','N'}; normopt='countdensity';
f.ax(1)=plt_dlochist(f.ax(1),dloc,xran,binw,legendstr1,normopt);
legendstr2={sprintf('%d ^o',round(anglegeo(2))), ...
  sprintf('%d ^o',round(anglegeo(1)))};
f.ax(2)=plt_dlochist(f.ax(2),dprojloc,xran,binw,legendstr2,normopt);  
distprojloc=median(abs(dprojloc(:,1)));

%strip along the short direction
ind1=find(abs(dprojloc(:,2))<=wid1);
dlocplt1=dloc(ind1,:);
dprojlocplt1=dprojloc(ind1,:);
title(f.ax(3),sprintf('A strip along %d^o of %.1f km wide, %d measurements',...
  round(anglegeo(2)),2*wid1,size(dlocplt1,1)));
% legendstr1={'E','N'};
% f.ax(3)=plt_dlochist(f.ax(3),dlocplt,xran,binw,legendstr1,normopt);
legendstr2={sprintf('%d ^o',round(anglegeo(2)))};
f.ax(3)=plt_dlochist(f.ax(3),dprojlocplt1(:,1),xran,binw,legendstr2,normopt);  
distprojlocplt=median(abs(dprojlocplt1(:,1)));

%strip along the short direction
wid2=0.15;
ind2=find(abs(dprojloc(:,1))<=wid2);
dlocplt2=dloc(ind2,:);
dprojlocplt2=dprojloc(ind2,:);
title(f.ax(4),sprintf('A strip along %d^o of %.1f km wide, %d measurements',...
  round(anglegeo(1)),2*wid2,size(dlocplt2,1)));
% xran=[-4 4]; binw=0.1; legendstr1={'E','N'}; normopt='countdensity';
% f.ax(3)=plt_dlochist(f.ax(3),dlocplt,xran,binw,legendstr1,normopt);
legendstr2={sprintf('%d ^o',round(anglegeo(1)))};
f.ax(4)=plt_dlochist(f.ax(4),dprojlocplt2(:,2),xran,binw,legendstr2,normopt);  
distprojlocplt=median(abs(dprojlocplt2(:,1)));

%%
for dmu=5:5:30
x=(1:100)';
mu0=45;
sigma=10;
g1=normpdf(x,mu0,sigma);
g2=normpdf(x,mu0+dmu,sigma);
figure; plot(x,g1,'k'); hold on; plot(x,g2,'k');
g3=g1+g2;
plot(x,g3,'r');
g4=g3/sum(g3);
plot(x,g4,'r','linew',2);
% y=randpdf(g4,x,[120,1]);
y=randpdf(x,g4,[10000,1]);
% histogram(y,'binw',2,'Normalization','pdf');
[MUHAT1,SIGMAHAT1] = normfit(y)
g5=normpdf(x,MUHAT1,SIGMAHAT1);
plot(x,g5,'b','linew',2);
text(0.1,0.9,sprintf('\\Delta\\mu=%d, \\sigma=%d, %.1f',dmu,sigma,SIGMAHAT1),'Units','normalized',...
  'HorizontalAlignment','left');
end

%%
%histogram of diff location in various ways
f = initfig(12,4.5,1,3);
xran=[-4 4]; binw=0.1; legendstr1={'E','N'}; normopt='countdensity';
f.ax(1)=plt_dlochist(f.ax(1),dloc,xran,binw,legendstr1,normopt);
legendstr2={'S45E','N45E'};
f.ax(2)=plt_dlochist(f.ax(2),dprojloc,xran,binw,legendstr2,normopt);
legendstr3={'S45E-cs'};
f.ax(3)=plt_dlochist(f.ax(3),dprojloccs(:,1),xran,binw,legendstr3,normopt);  

distproj=median(abs(dprojloc(:,1)))

keyboard

%%%summarize the diff location between N and N-n for different clusters with m+1 events
f = initfig(8,4.5,1,2);
color=jet(mmax);
X=xran(1):binw:xran(2);
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
for m=n:mmax
  Y1=normpdf(X,mu1(m),sigma1(m));
  p(m) = plot(ax,X,Y1,'Color',color(m,:),'linew',1);
  label{m} = sprintf('m=%d, \\sigma=%.2f',m,sigma1(m));
end
ylabel(ax,'PDF of Gaussian fit');
xlabel(ax,'Location difference along S45E (km)');
legend(ax,p(n:mmax),label{n:mmax},'NumColumns',1,'Location','east');
ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
for m=1:mmax
  Y2=normpdf(X,mu2(m),sigma2(m));
  p(m) = plot(ax,X,Y2,'Color',color(m,:),'linew',1);
  label{m} = sprintf('m=%d, \\sigma=%.2f',m,sigma2(m));
end
ylabel(ax,'PDF of Gaussian fit');
xlabel(ax,'Location difference along N45E (km)');
legend(ax,p(n:mmax),label{n:mmax},'NumColumns',1,'Location','east');

keyboard


%% various types of distance between events in the clusters w/i a diff time cut
% for m=1:mmax
m=5;
timetype = 'tarvl';
% timetype = 'tori';
dtcut = 0.25*m+0.125;
[ampplt,dtplt,indimpdtcut] = med_amp_incluster(nbst,imp,nsrc,m,dtcut,sps,timetype);
amp=mean(imp(:,[2 4 6]),2);
ampch=prctile(amp,5);

[impcluster,impuni,ncluster]=evtcluster(nbst,imp,nsrc,m,dtcut,sps,timetype);
imploc = off2space002(impuni(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
implocclus = off2space002(impcluster(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
imploc0 = off2space002([0 0],sps,ftrans,0);  % a ref source at 0,0

%%%compute the distance 
%this would take a while
tic
[mdistnn1,mdist2all,pca2vec,projang1,mprojx1nn1,mprojx12all,...
  projang2,mprojx2nn1,mprojx22all]=...
  dist_evtcluster(impcluster,implocclus,sps,ftrans,m,timetype);
% [mdistnn1,mdist2all,pca2vec,projang1,mprojx1nn1,mprojx12all]=...
%   dist_evtcluster(impcluster,implocclus,sps,ftrans,m,'tori');
toc
% keyboard

%%%summarize the median distance comparison between consecutive ones, and each to all others
f = initfig(12,4.5,1,3);
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
p1=histogram(ax,mdistnn1,'binw',0.05,'normalization','count','Facec','b');
p2=histogram(ax,mdist2all,'binw',0.05,'normalization','count','Facec','r');
plot(ax,[median(mdistnn1) median(mdistnn1)],ax.YLim,'b--','LineWidth',1);
plot(ax,[median(mdist2all) median(mdist2all)],ax.YLim,'r--','LineWidth',1);
xlabel(ax,'Absolute distance (km)');
ylabel(ax,'Count');
legend(ax,[p1,p2],'N and N-1','each to others');
ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
histogram(ax,mprojx1nn1,'binw',0.05,'normalization','count','Facec','b');
histogram(ax,mprojx12all,'binw',0.05,'normalization','count','Facec','r');
plot(ax,[median(mprojx1nn1) median(mprojx1nn1)],ax.YLim,'b--','LineWidth',1);
plot(ax,[median(mprojx12all) median(mprojx12all)],ax.YLim,'r--','LineWidth',1);
xlabel(ax,'Distance along the 2nd PC (km)');
ylabel(ax,'Count');
ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
histogram(ax,mprojx2nn1,'binw',0.05,'normalization','count','Facec','b');
histogram(ax,mprojx22all,'binw',0.05,'normalization','count','Facec','r');
plot(ax,[median(mprojx2nn1) median(mprojx2nn1)],ax.YLim,'b--','LineWidth',1);
plot(ax,[median(mprojx22all) median(mprojx22all)],ax.YLim,'r--','LineWidth',1);
xlabel(ax,'Distance along the prop. direc. (km)');
ylabel(ax,'Count');
keyboard

%% normalized total covered area for the clustered events assuming a diameter
for ii = 1: length(diams)
  tic
  diam=diams(ii); %dimension of source, would be diamter if circle
  radi=0.5*diam; %radius of blue circles, max. with no overlapping
  for i = 1: ncluster
    ist=(i-1)*(m+1)+1;
    ied=i*(m+1);
    arearat(i,ii) = area_of_overlap_circs_grid(implocclus(ist:ied,1:2),...
      radi,0.01,0);
    % [arearat(i,ii)] = area_of_overlap_circs_MC(implocclus(ist:ied,1:2),...
    %   radi,5e4,10,1);
    % keyboard
    % close all  
  end
  toc
end

%%
%plot the distribution of total covered area normalized by norminal area
f = initfig(8,8,2,2); %initialize fig
binw=0.02;
binedge=0:binw:1;
for ii = 1:length(diams)
  diam=diams(ii);
  ax=f.ax(ii); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  histogram(ax,arearat(:,ii),'BinEdges',binedge);
  text(ax,0.98,0.95,sprintf('diameter: %.1f km',diam),'Units','normalized',...
      'HorizontalAlignment','right');
  xlabel(ax,'Normalized area');
  ylabel(ax,'Count');
  xlim(ax,[0 1]);
  %there is a lower limit of the ratio of covered area
  plot(ax,[1/(m+1) 1/(m+1)],ax.YLim,'k--');
end
supertit(f.ax(1:2),sprintf('m=%d, ncluster=%d',m,ncluster),10);
% keyboard
figname = sprintf('arearat_m%d_nclus%d.pdf',m,ncluster);
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',figname));

%% plot locations of srcs with finite size in each clusters assuming a diameter
for ii = 1: length(diams)
  diam=diams(ii);
  radi=0.5*diam;

  nrow1 = 4; % rows and cols of subplots in each figure
  ncol1 = 5; 
  widin1 = 11; % size of each figure
  htin1 = 8.5;
  pltxran1 = [0.05 0.96]; pltyran1 = [0.05 0.96]; % optimal axis location
  pltxsep1 = 0.03; pltysep1 = 0.03;
  ifig1 = 0; % figure number

  ifig1 = ifig1+1;
  f1 = initfig(widin1,htin1,nrow1,ncol1,ifig1);
  axpos1 = optaxpos(f1,nrow1,ncol1,pltxran1,pltyran1,pltxsep1,pltysep1);
  axtit = f1.ax(1: ncol1);
  supertit(axtit, sprintf('Fig %s, m=%d, nclus=%d, diam=%.1f km',...
    num2zeropadstr(ifig1, 3),m,ncluster,diam),10);
  xlabel(f1.ax((nrow1-1)*ncol1+1),'E (km)','fontsize',10);
  ylabel(f1.ax((nrow1-1)*ncol1+1),'N (km)','fontsize',10);
  orient(f1.fig,'landscape');
  isub = 0;  % count of the total subplots on the current figure f1

  for i = 1: ncluster
    %events in each cluster
    ist=(i-1)*(m+1)+1;
    ied=i*(m+1);
    implocplt=implocclus(ist:ied,:);
    impplt=impcluster(ist:ied,:);
    
    isub = isub+1; 
    if isub > nrow1*ncol1
      %print the current figure f1
      f1name{ifig1,1} = sprintf('clusterloc_Fig%s_m%d_%.1fkm.pdf',...
        num2zeropadstr(ifig1, 3),m,diam);
      print(f1.fig,'-dpdf','-fillpage',strcat('/home/chaosong/Pictures/',f1name{ifig1,1}));
      %move on to next figure f1
      ifig1 = ifig1+1;  %initialize another figure
      f1 = initfig(widin1,htin1,nrow1,ncol1,ifig1);
      axpos1 = optaxpos(f1,nrow1,ncol1,pltxran1,pltyran1,pltxsep1,pltysep1);
      axtit = f1.ax(1: ncol1);
      supertit(axtit, sprintf('Fig %s, m=%d, nclus=%d, diam=%.1f km',...
        num2zeropadstr(ifig1, 3),m,ncluster,diam),10);
      xlabel(f1.ax((nrow1-1)*ncol1+1),'E (km)','fontsize',10);
      ylabel(f1.ax((nrow1-1)*ncol1+1),'N (km)','fontsize',10);
      orient(f1.fig,'landscape');
      
      isub = isub-nrow1*ncol1;  %refresh
    end

    %in terms of origin time
    tcor = round((implocplt(:,6)-imploc0(6))*sps);   % travel time difference between each source and a ref source at 0,0
    torispl = impplt(:,1)-tcor;    
    [torisplst, indsort] = sortrows(torispl,1);
    implocpltst = implocplt(indsort, :);
    
    %plot locations of events in each cluster
    color=jet(m+1);
    ax=f1.ax(isub); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    axis(ax, 'equal');
    for j=1: m+1
      [xcut,ycut] = circle_chao(implocpltst(j,1),implocpltst(j,2),radi,0.1);
      plot(ax,xcut,ycut,'-','color',color(j,:),'linew',1);
    end
    x0=mean(implocpltst(:,1));
    y0=mean(implocpltst(:,2));
    plot(ax,x0+[pca2vec(i,1), -pca2vec(i,1)],y0+[pca2vec(i,2), -pca2vec(i,2)],...
      'k-','linew',1);
    text(ax,0.02,0.9,sprintf('amp:%.2f; %.1fx',ampplt(i),ampplt(i)/ampch),...
      'Units','normalized','HorizontalAlignment','left','fontsize',7);   
    text(ax,0.02,0.85,sprintf('arearat:%.2f',arearat(i,ii)),...
      'Units','normalized','HorizontalAlignment','left','fontsize',7);
    text(ax,0.98,0.15,sprintf('%.2f; %.2f km',mdistnn1(i,1),mdist2all(i,1)),...
      'Units','normalized','HorizontalAlignment','right','fontsize',7);
    text(ax,0.98,0.1,sprintf('%d / %d ^o; %.2f; %.2f km',round(projang1(i,1)),...
      round(projang1(i,1))+180,mprojx1nn1(i,1),mprojx12all(i,1)),...
      'Units','normalized','HorizontalAlignment','right','fontsize',7);
    text(ax,0.98,0.05,sprintf('%d ^o; %.2f; %.2f km',round(projang2(i,1)),...
      mprojx2nn1(i,1),mprojx22all(i,1)),...
      'Units','normalized','HorizontalAlignment','right','fontsize',7);
    xlim(ax,[-4 4]);
    ylim(ax,[-4 4]);
    longticks(ax,2);
    hold(ax,'off');
  end
  %if all clusters have been plotted, print the current figure f1
  f1name{ifig1,1} = sprintf('clusterloc_Fig%s_m%d_%.1fkm.pdf',...
    num2zeropadstr(ifig1, 3),m,diam);
  print(f1.fig,'-dpdf','-fillpage',strcat('/home/chaosong/Pictures/',f1name{ifig1,1}));

  %merge all figures into a single pdf file
  mergef1name=sprintf('clusterloc_m%d_%.1fkm.pdf',m,diam);
  str=strcat('rm -f /home/chaosong/Pictures/',mergef1name);
  status = system(str);
  % status = system('rm -f /home/chaosong/Pictures/clusterloc.pdf');
  for i = 1:size(f1name,1)
    fname{i} = strcat('/home/chaosong/Pictures/',f1name{i});
  end
  append_pdfs(strcat('/home/chaosong/Pictures/',mergef1name),fname);
  status = system('rm -f /home/chaosong/Pictures/clusterloc_Fig*.pdf');
  % keyboard
  close all
end
  




