% lfedloc_incluster.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to 'circsrcmodel.m', this is another script to assume a circular
% source model. Now we have locations
% of all deconvolved LFE sources where each source is represented by a circle
% with a series of diameter. This script aims to show several things:
% 1. in any certain period of time, eg, 20 s, during which sources won't
% migrate for too far, what are distribution of the sources?
% --This code can be considered as the supplementary information for getting
% the distance from each source to all other sources whose ABS arrival time
% difference is within 10 s, +/- 10s == 20 s.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/01/08
% Last modified date:   2024/01/08
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
  [dloc,dt,k,dloc_spl]=dloc_evtcluster(impcluster,implocclus,sps,ftrans,m,n,timetype);
  %note that 'dloc' is NOT the same as the direct inversion of 'dloc_spl'
  dtlump{m}=dt;
  dloclump{m}=dloc;
  dlocspllump{m}=dloc_spl;
  
%   %plot the diff location in samples, bin by pixel
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
%  
%   %plot the diff location in map, bin by grid
%   [coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc); 
%   cstr={'# events / grid'}; xran=[-4 4]; yran=[-4 4]; dx=0.1; dy=0.1;
%   [f,den1d,conmat]=plt_srcdlocinmap(dt/sps,dloc,[],'km',timetype,...
%     10,cstr,'o','linear','grid',xran,yran,dx,dy,1,2);
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
%plot the diff location in samples, bin by pixel
smoothsigma=1;
contintvl=1;
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc_spl);
cstr={'# events / pixel'}; xran=[-50 50]; yran=[-50 50]; dx=1; dy=1;
[f,den1d_spl,conmat]=plt_srcdlocinmap(dt/sps,dloc_spl,[],'spl',timetype,...
  8,cstr,'o','linear','pixel',xran,yran,dx,dy,1,smoothsigma,contintvl);
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
ax=f.ax(3);
hold(ax,'on');
text(ax,0.75,0.1,sprintf('smooth sigma: %.1f',smoothsigma),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.75,0.05,sprintf('interval: %d',contintvl),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
hold(ax,'off');
% keyboard

%%
%plot the diff location in map, bin by grid
smoothsigma=1.5;
contintvl=1;
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc);
cstr={'# events / grid'}; xran=[-4 4]; yran=[-4 4]; dx=0.1; dy=0.1;
[f,den1d,conmat]=plt_srcdlocinmap(dt/sps,dloc,[],'km',timetype,...
  8,cstr,'o','linear','grid',xran,yran,dx,dy,1,smoothsigma,contintvl);
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
ax=f.ax(3);
hold(ax,'on');
F = scatteredInterpolant(conmat(:,3),conmat(:,4),conmat(:,1),'linear','none');
%a line cross (0,0) in the projection direction
x = reshape(xran(1):0.01:xran(2), [], 1);
yopt = linefcn(x,tand(angle(2)),0);
yopt1 = linefcn(x,tand(angle(2)),0.25/cosd(angle(2)));
yopt2 = linefcn(x,tand(angle(2)),-0.25/cosd(angle(2)));
plot(ax,x,yopt,'k-','linew',2);
plot(ax,x,yopt1,'k--','linew',1);
plot(ax,x,yopt2,'k-.','linew',1);
%a line cross (0,0) in the orthogonal direction
yort = linefcn(x,tand(angle(1)),0);
yort1 = linefcn(x,tand(angle(1)),-0.25/cosd(angle(1)));
yort2 = linefcn(x,tand(angle(1)),0.25/cosd(angle(1)));
plot(ax,x,yort,'-','linew',2,'color',[.5 .5 .5]);
plot(ax,x,yort1,'--','linew',1,'color',[.5 .5 .5]);
plot(ax,x,yort2,'-.','linew',1,'color',[.5 .5 .5]);
text(ax,0.75,0.1,sprintf('smooth sigma: %.1f',smoothsigma),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.75,0.05,sprintf('interval: %d',contintvl),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);

%% overall diff loc distribution and the projection
f = initfig(8,4.5,1,2);
optaxpos(f,1,2);
supertit(f.ax(1:2),sprintf('all, %d dpts',size(dloc,1)),10);
xran=[-4 4]; binw=0.1; legendstr1={'E','N'};
normopt='custom'; normalizer=length(dloc);
f.ax(1)=plt_dlochist(f.ax(1),dloc,xran,binw,legendstr1,normopt,normalizer);

legendstr2={sprintf('%d ^o',round(anglegeo(2))), ...
  sprintf('%d ^o',round(anglegeo(1)))};
[~,~,dprojloc] = customprojection(dloc(:,1:2),anglegeo(2));
f.ax(2)=plt_dlochist(f.ax(2),dprojloc,xran,binw,legendstr2,normopt,normalizer);
distprojloc=median(abs(dprojloc(:,1)));

%% if choosing a line intersecting with the contour
f = initfig(8,4.5,1,2);
optaxpos(f,1,2);
[f.ax(1),muopt,sigmaopt,mdistprojopt]=plt_dloccrssect(f.ax(1),F,x,yopt,anglegeo(2),...
  [0 0 0; 0 0 1],xran,normalizer,yopt1,yopt2);

[f.ax(2),muort,sigmaort,mdistprojort]=plt_dloccrssect(f.ax(2),F,x,yort,anglegeo(1),...
  [.5 .5 .5; 1 0 0],xran,normalizer,yort1,yort2);

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
cstr={'# events / grid'}; xran=[-4 4]; yran=[-4 4]; dx=0.1; dy=0.1;
[f,den1d]=plt_srcdlocinmap(dt/sps,dprojloc,[],'km',timetype,...
  10,cstr,'o','linear','grid',xran,yran,dx,dy,0);
ax=f.ax(2);
hold(ax,'on');
plot(ax,ax.XLim,[wid1 wid1],'-','linew',2,'color',[.7 .7 .7]);
plot(ax,ax.XLim,[-wid1 -wid1],'-','linew',2,'color',[.7 .7 .7]);
hold(ax,'off');
ax=f.ax(1);
hold(ax,'on');
plot(ax,ax.XLim,[wid1 wid1],'-','linew',2,'color',[.7 .7 .7]);
plot(ax,ax.XLim,[-wid1 -wid1],'-','linew',2,'color',[.7 .7 .7]);
hold(ax,'off');

%%
widpool=[0.3 0.5 0.7 0.9]/2;
for i=4:length(widpool)
  wid=widpool(i);
  f = initfig(8,8,2,2);
  supertit(f.ax(1:2),sprintf('all, %d dpts',size(dloc,1)),10);
  xran=[-4 4]; binw=0.2; legendstr1={'E','N'}; 
  normopt='custom'; normalizer=length(dloc);
  f.ax(1)=plt_dlochist(f.ax(1),dloc,xran,binw,legendstr1,normopt,normalizer);
  legendstr2={sprintf('%d ^o',round(anglegeo(2))), ...
    sprintf('%d ^o',round(anglegeo(1)))};
  f.ax(2)=plt_dlochist(f.ax(2),dprojloc,xran,binw,legendstr2,normopt,normalizer);  
  distprojloc=median(abs(dprojloc(:,1)));

  %strip along the short direction
  ind1=find(abs(dprojloc(:,2))<=wid);
  dlocplt1=dloc(ind1,:);
  dprojlocplt1=dprojloc(ind1,:);
  dlocsplplt1=dloc_spl(ind1,:);
  title(f.ax(3),sprintf('Along %d^o, %.1f km wide, %d dpts',...
    round(anglegeo(2)),2*wid,size(dlocplt1,1)));
  % legendstr1={'E','N'};
  % f.ax(3)=plt_dlochist(f.ax(3),dlocplt,xran,binw,legendstr1,normopt);
  normopt='custom'; normalizer=length(dlocplt1);
  legendstr2={sprintf('%d ^o',round(anglegeo(2)))};
  f.ax(3)=plt_dlochist(f.ax(3),dprojlocplt1(:,1),xran,binw,legendstr2,...
    normopt,normalizer);  
  distprojlocplt=median(abs(dprojlocplt1(:,1)));
  
  xedge=xran(1)+binw/2:binw:xran(2)-binw/2;
  xcnt=0.5*(xedge(1:end-1)+xedge(2:end));
  nunispl=zeros(length(xedge)-1,1);
  nuni=zeros(length(xedge)-1,1);
  Nnspl=zeros(length(xedge)-1,1);
  Nn=zeros(length(xedge)-1,1);
  for j=1:length(xedge)-1
    ind=find(dprojlocplt1(:,1)>=xedge(j)&dprojlocplt1(:,1)<xedge(j+1));
    if ~isempty(ind)      
      nunispl(j)=size(unique(dlocsplplt1(ind,:),'row'),1);
      nuni(j)=size(unique(dlocplt1(ind,:),'row'),1);
      Nnspl(j)=length(ind)/nunispl(j);
      Nn(j)=length(ind)/nuni(j);
    end
  end
  f1=initfig(8,4,1,2);
  ax=f1.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
  bar(ax,xcnt,nuni,1,'stacked','b','facea',0.6);
  bar(ax,xcnt,nunispl,1,'stacked','r','facea',0.6);
  legend(ax,'relative loc.','sample space');
  ylabel(ax,'Unique points');
  xlabel(ax,'Location difference (km)');
  ax=f1.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
  bar(ax,xcnt,Nn,1,'stacked','b','facea',0.6);
  bar(ax,xcnt,Nnspl,1,'stacked','r','facea',0.6);
  ylabel(ax,'Normalized count');
  xlabel(ax,'Location difference (km)');
  
  %strip along the short direction
  % wid2=0.15;
  ind2=find(abs(dprojloc(:,1))<=wid);
  dlocplt2=dloc(ind2,:);
  dprojlocplt2=dprojloc(ind2,:);
  title(f.ax(4),sprintf('Along %d^o, %.1f km wide, %d dpts',...
    round(anglegeo(1)),2*wid,size(dlocplt2,1)));
  % xran=[-4 4]; binw=0.1; legendstr1={'E','N'}; normopt='countdensity';
  % f.ax(3)=plt_dlochist(f.ax(3),dlocplt,xran,binw,legendstr1,normopt);
  legendstr2={sprintf('%d ^o',round(anglegeo(1)))};
  f.ax(4)=plt_dlochist(f.ax(4),dprojlocplt2(:,2),xran,binw,legendstr2,...
    normopt,normalizer);    
  distprojlocplt=median(abs(dprojlocplt2(:,1)));
end

%%
for dmu=5:5:30
x=(1:100)';
mu0=45;
sigma=2;
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
dmu=0.4;
x=(-4:0.01:4)';
mu0=-0.2;
sigma=1;
g1=normpdf(x,mu0,sigma);
g2=normpdf(x,mu0+dmu,sigma);
figure; 
subplot(121)
plot(x,g1,'k'); hold on; plot(x,g2,'k');
g3=g1+g2;
plot(x,g3,'r');
text(0.1,0.9,sprintf('\\Delta\\mu=%.2f, \\sigma=%.2f',dmu,sigma),'Units','normalized',...
  'HorizontalAlignment','left');
fitobj = fit(x,g3,'gauss2','StartPoint',[1 mu0 1 1 mu0+dmu 1]);
coef = coeffvalues(fitobj);
subplot(122)
g4=feval(fitobj,x);
plot(x,g4,'b','linew',1); hold on;
a1 = coef(1);
mu1 = coef(2);
sigma1 = coef(3);
a2 = coef(4);
mu2 = coef(5);
sigma2 = coef(6);
g5 = a1.*exp(-((x-mu1)./sigma1).^2);
g6 = a2.*exp(-((x-mu2)./sigma2).^2);
plot(x,g5,'k'); hold on; plot(x,g6,'k');
text(0.1,0.9,sprintf('\\mu=%.2f; %.2f, \n \\sigma=%.2f; %.2f',mu1,sigma1,...
  mu2,sigma2),'Units','normalized',...
  'HorizontalAlignment','left');


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
  Y1=normpdf(X,muopt(m),sigmaopt(m));
  p(m) = plot(ax,X,Y1,'Color',color(m,:),'linew',1);
  label{m} = sprintf('m=%d, \\sigma=%.2f',m,sigmaopt(m));
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


