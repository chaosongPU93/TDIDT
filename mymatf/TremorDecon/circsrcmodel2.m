% circsrcmodel2.m
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
savefile = 'deconv_stats4th_allbstsig.mat';
load(strcat(rstpath, '/MAPS/',savefile));
savefile = 'deconv_stats4th_allbstnoi.mat';
load(strcat(rstpath, '/MAPS/',savefile));

% keyboard

%%
%%%param for secondary sources removed
locxyprojall = allbstsig.locxyprojall;
tarvlsplstall = allbstsig.impindepall(:,1);
nsrc = allbstsig.nsrc;
imp = allbstsig.impindepall;
off1i = allbstsig.off1i;
trangenew = allbstsig.trangenew;
locxyprojalln = allbstnoi.locxyprojall;
tarvlsplstalln = allbstnoi.impindepall(:,1);
nsrcn = allbstnoi.nsrc;
impn = allbstnoi.impindepall;
off1in = allbstnoi.off1i;
supertstr = 'Secondary sources removed';
fnsuffix = [];

% %%%param for further checked at KLNB
% locxyprojall = allbstsig.locxyproj4thall;
% tarvlsplstall = allbstsig.impindep4thall(:,1);
% nsrc = allbstsig.nsrc4th;
% imp = allbstsig.impindep4thall;
% off1i = allbstsig.off1i;
% trangenew = allbstsig.trangenew;
% locxyprojalln = allbstnoi.locxyproj4thall;
% tarvlsplstalln = allbstnoi.impindep4thall(:,1);
% nsrcn = allbstnoi.nsrc4th;
% impn = allbstnoi.impindep4thall;
% off1in = allbstnoi.off1i;
% supertstr = 'Further checked at KLNB';
% fnsuffix = '4th';

% keyboard

% impuse = imp;
% nsrcuse = nsrc;
impuse = impn;
nsrcuse = nsrcn;

[imploc, ~] = off2space002(impuse(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[imploc0, ~] = off2space002([0 0],sps,ftrans,0);  % a ref source at 0,0
amp=mean(impuse(:,[2 4 6]),2);
ampch=prctile(amp,5);

tlennew = trangenew(:,3)-trangenew(:,2);
tlensumnew = sum(tlennew);

timetype = 'tarvl';

%% distance between events, including PCA
subwsectar = 25/2;  %target subwindow length in sec
k = 0;  % count of the total subplots
n = 1;  %between N and N-n
m = 1;  %max n to compute
dtnn1cat=[];
dlocnn1cat=[];
dloc_splnn1cat=[];
dloc2allcat=[];
dt2allcat=[];
dloc2all_splcat=[];
distnn1cat=[];
dist2allcat=[];
dprojxy1nn1cat=[];
dprojxy12allcat=[];
dprojxy2nn1cat=[];
dprojxy22allcat=[];
%%
for i = 1: nbst
  disp(i)
  ist = sum(nsrcuse(1:i-1))+1;
  ied = ist+nsrcuse(i)-1;
  impi = impuse(ist:ied,:);
  imploci = imploc(ist:ied,:);
  %     %Principal component analysis
  %     [coeff,score,latent,tsquared,explained] = pca(imploci(:,1:2));
  %     %each column in coeff represent the unit vector of each principal component
  %     angle=atan2d(coeff(2,2),coeff(1,2));
  %     [anggeo,anggeooppo]=angatan2d2geo(angle);
  %     ang=min([anggeo,anggeooppo]);
  
  %bursts and 4-s detections of the same day
  indst=1;
  inded=tlennew(i)*sps;
  if tlennew(i)<subwsectar
    subwsec=tlennew(i);   %actually used subwindow length in sec  
  else
    subwsec=subwsectar;
  end
  subwlen=subwsec*sps;
  ovlplen=0;
  windows = movingwins(indst,inded,subwlen,ovlplen,0);
  nwin =  size(windows,1);
  iwin = findwhichrange(impi(:,1),windows);
  for j = 1: nwin
    impiwin = impi(iwin==j,:);
    implociwin = imploci(iwin==j,:);
    lwinsec = round(((windows(j,2)-windows(j,1))+1)/sps); %notes the actual subwin length in sec
    
    if size(impiwin,1) >1
      k=k+1;
      ampi(k,1) = median(mean(impiwin(:,[2 4 6]),2));

      %%%%%% diff loc within short win, able to be combined to analyse later
      %compute the diff loc between N and N-m, and each to all others in the
      %short window, no projection is applied
      [dloc,dt,dloc_spl,dloc2all,dt2all,dloc2all_spl]=dloc_evtcustom(impiwin,implociwin,sps,ftrans,m,'tarvl');
      mdtnn1(k,1)=median(dt{n});
      mdlocnn1(k,:)=median(dloc{n});
      mdloc_splnn1(k,:)=median(dloc_spl{n});
      mdloc2all(k,:)=median(dloc2all);
      mdt2all(k,1)=median(dt2all);
      mdloc2all_spl(k,:)=median(dloc2all_spl);
      dtnn1cat=[dtnn1cat; dt{n}];
      dlocnn1cat=[dlocnn1cat; dloc{n}];
      dloc_splnn1cat=[dloc_splnn1cat; dloc_spl{n}];
      dloc2allcat=[dloc2allcat; dloc2all];
      dt2allcat=[dt2allcat; dt2all];
      dloc2all_splcat=[dloc2all_splcat; dloc2all_spl];
      %%%%%%%%%%

      %%%%%% distance, etc for each short win
      %compute the abs distance between N and N-m, and each to all others in the
      %short window. For all srcs in each win, apply the PCA, and project to 
      %this direction. Note that this direction is not the same as to all srcs
      %combined
      [distnn1,dist2all,pcavec,ang1,dprojxy1nn1,dprojxy12all,...
        ang2,dprojxy2nn1,dprojxy22all]=dist_evtcustom(impiwin,implociwin,sps,ftrans,'tarvl'); %tori
      pca2vec(k,:)=pcavec;
      projang1(k,1)=ang1;
      projang2(k,1)=ang2;
      mdistnn1(k,1)=median(distnn1);
      mdist2all(k,1)=median(dist2all);
      mdprojx1nn1(k,1)=median(abs(dprojxy1nn1(:,1)));
      mdprojx12all(k,1)=median(abs(dprojxy12all(:,1)));
      mdprojx2nn1(k,1)=median(abs(dprojxy2nn1(:,1)));
      mdprojx22all(k,1)=median(abs(dprojxy22all(:,1)));
      distnn1cat=[distnn1cat; distnn1];
      dist2allcat=[dist2allcat; dist2all];
      dprojxy1nn1cat=[dprojxy1nn1cat; dprojxy1nn1];
      dprojxy12allcat=[dprojxy12allcat; dprojxy12all];
      dprojxy2nn1cat=[dprojxy2nn1cat; dprojxy2nn1];
      dprojxy22allcat=[dprojxy22allcat; dprojxy22all];
      %%%%%%%%%%
      
    end
  end
end
% keyboard

%% which type of location difference or distance to look at
% %%%%%%%% between consecutive ones
% dloc_spl = dloc_splnn1cat;
% dloc = dlocnn1cat;
% dt = dtnn1cat;
% dist = distnn1cat;
% dprojxy1 = dprojxy1nn1cat;
% dprojxy2 = dprojxy2nn1cat;
% str=sprintf('between N and N-%d',n);
%%%%%%%% between each to all others
dloc_spl = dloc2all_splcat;
dloc = dloc2allcat;
dt = dt2allcat;
dist = dist2allcat;
dprojxy1 = dprojxy12allcat;
dprojxy2 = dprojxy22allcat;
str=sprintf('between each to others');

%% summariz the location difference between consecutive ones, and each to all others
%plot the diff location in samples, bin by pixel
smoothsigma=1;
contintvl=1;
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc_spl);
cstr={'# events / pixel'}; xran=[-40 40]; yran=[-40 40]; dx=1; dy=1;
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
supertit(f.ax, strcat({'location difference '},str));
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
smoothsigma=1;
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
hold(ax,'off');

supertit(f.ax, strcat({'Map-view, location difference '},str));


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


%% if choose a strip with a finite width in the PCA direction
[~,~,dprojloc] = customprojection(dloc(:,1:2),anglegeo(2));

wid=0.45;
cstr={'# events / grid'}; xran=[-4 4]; yran=[-4 4]; dx=0.1; dy=0.1;
[f,den1d]=plt_srcdlocinmap(dt/sps,dprojloc,[],'km',timetype,...
  8,cstr,'o','linear','grid',xran,yran,dx,dy,0);
ax=f.ax(2);
hold(ax,'on');
plot(ax,ax.XLim,[wid wid],'-','linew',2,'color',[.7 .7 .7]);
plot(ax,ax.XLim,[-wid -wid],'-','linew',2,'color',[.7 .7 .7]);
hold(ax,'off');
ax=f.ax(1);
hold(ax,'on');
plot(ax,ax.XLim,[wid wid],'-','linew',2,'color',[.7 .7 .7]);
plot(ax,ax.XLim,[-wid -wid],'-','linew',2,'color',[.7 .7 .7]);
hold(ax,'off');

f = initfig(8,8,2,2);
supertit(f.ax(1:2),sprintf('all, %d dpts',size(dloc,1)),10);
xran=[-4 4]; binw=0.2; legendstr1={'E','N'}; 
normopt='custom'; normalizer=length(dloc);
f.ax(1)=plt_dlochist(f.ax(1),dloc,xran,binw,legendstr1,normopt,normalizer);
legendstr2={sprintf('%d ^o',round(anglegeo(2))), ...
  sprintf('%d ^o',round(anglegeo(1)))};
f.ax(2)=plt_dlochist(f.ax(2),dprojloc,xran,binw,legendstr2,normopt,normalizer);  
distprojlocnn1=median(abs(dprojloc(:,1)));

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
distprojlocplt1=median(abs(dprojlocplt1(:,1)));

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
dlocsplplt2=dloc_spl(ind2,:);
title(f.ax(4),sprintf('Along %d^o, %.1f km wide, %d dpts',...
  round(anglegeo(1)),2*wid,size(dlocplt2,1)));
% xran=[-4 4]; binw=0.1; legendstr1={'E','N'}; normopt='countdensity';
% f.ax(3)=plt_dlochist(f.ax(3),dlocplt,xran,binw,legendstr1,normopt);
legendstr2={sprintf('%d ^o',round(anglegeo(1)))};
f.ax(4)=plt_dlochist(f.ax(4),dprojlocplt2(:,2),xran,binw,legendstr2,...
  normopt,normalizer);    
distprojlocplt2=median(abs(dprojlocplt2(:,1)));


%% summarize the distance comparison between consecutive ones, and each to all others
% f = initfig(12,4.5,1,3);
% ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% % p1=histogram(ax,distnn1cat,'binw',0.05,'normalization','count','Facec','b');
% p2=histogram(ax,dist2allcat,'binw',0.05,'normalization','count','Facec','r');
% % plot(ax,[median(distnn1cat) median(distnn1cat)],ax.YLim,'b--','LineWidth',1);
% plot(ax,[median(dist2allcat) median(dist2allcat)],ax.YLim,'r--','LineWidth',1);
% xlabel(ax,'Absolute distance (km)');
% ylabel(ax,'Count');
% legend(ax,[p2],'each to others');
% ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% % histogram(ax,abs(dprojxy1nn1cat(:,1)),'binw',0.05,'normalization','count','Facec','b');
% histogram(ax,abs(dprojxy12allcat(:,1)),'binw',0.05,'normalization','count','Facec','r');
% % plot(ax,[median(abs(dprojxy1nn1cat(:,1))) median(abs(dprojxy1nn1cat(:,1)))],ax.YLim,'b--','LineWidth',1);
% plot(ax,[median(abs(dprojxy12allcat(:,1))) median(abs(dprojxy12allcat(:,1)))],ax.YLim,'r--','LineWidth',1);
% xlabel(ax,'Distance along the 2nd PC (km)');
% ylabel(ax,'Count');
% ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% % histogram(ax,abs(dprojxy2nn1cat(:,1)),'binw',0.05,'normalization','count','Facec','b');
% histogram(ax,abs(dprojxy22allcat(:,1)),'binw',0.05,'normalization','count','Facec','r');
% % plot(ax,[median(abs(dprojxy2nn1cat(:,1))) median(abs(dprojxy2nn1cat(:,1)))],ax.YLim,'b--','LineWidth',1);
% plot(ax,[median(abs(dprojxy22allcat(:,1))) median(abs(dprojxy22allcat(:,1)))],ax.YLim,'r--','LineWidth',1);
% xlabel(ax,'Distance along the prop. direc. (km)');
% ylabel(ax,'Count');
% keyboard

%% summarize the MEDIAN distance comparison between consecutive ones, and each to all others
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
histogram(ax,mdprojx1nn1,'binw',0.05,'normalization','count','Facec','b');
histogram(ax,mdprojx12all,'binw',0.05,'normalization','count','Facec','r');
plot(ax,[median(mdprojx1nn1) median(mdprojx1nn1)],ax.YLim,'b--','LineWidth',1);
plot(ax,[median(mdprojx12all) median(mdprojx12all)],ax.YLim,'r--','LineWidth',1);
xlabel(ax,'Distance along the 2nd PC (km)');
ylabel(ax,'Count');
ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
histogram(ax,mdprojx2nn1,'binw',0.05,'normalization','count','Facec','b');
histogram(ax,mdprojx22all,'binw',0.05,'normalization','count','Facec','r');
plot(ax,[median(mdprojx2nn1) median(mdprojx2nn1)],ax.YLim,'b--','LineWidth',1);
plot(ax,[median(mdprojx22all) median(mdprojx22all)],ax.YLim,'r--','LineWidth',1);
xlabel(ax,'Distance along the prop. direc. (km)');
ylabel(ax,'Count');
figname = sprintf('distcomp_wlen%ds_nwin%d.pdf',subwsec,k);
orient(f.fig,'landscape');
print(f.fig,'-dpdf','-bestfit',strcat('/home/chaosong/Pictures/',figname));
keyboard

%% plot locations of srcs with finite size in short wins assuming a diameter
%choice of diameter of circular sources
diams = [0.1 0.3 0.5 0.7];
% arearat = cell(length(diams),1);
for ii = 4: length(diams)
  diam=diams(ii);
  radi=0.5*diam;
  
  subwsec=20;
  
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
  supertit(axtit, sprintf('Fig %s, wlen=%d s, diam=%.1f km',...
    num2zeropadstr(ifig1, 3),subwsec,diam),10);
  xlabel(f1.ax((nrow1-1)*ncol1+1),'E (km)','fontsize',10);
  ylabel(f1.ax((nrow1-1)*ncol1+1),'N (km)','fontsize',10);
  orient(f1.fig,'landscape');
  isub = 0;  % count of the total subplots on the current figure f1
  k = 0;  % count of the total subplots
  
  for i = 1: nbst
    %   i=181;
    ist = sum(nsrcuse(1:i-1))+1;
    ied = ist+nsrcuse(i)-1;
    impi = impuse(ist:ied,:);
    imploci = imploc(ist:ied,:);
    
    %bursts and 4-s detections of the same day
    indst=1;
    inded=tlen(i)*sps;
    if tlen(i)<20
      subwsec=tlen(i);
    else
      subwsec=20;
    end
    subwlen=subwsec*sps;
    ovlplen=0;
    windows = movingwins(indst,inded,subwlen,ovlplen,0);
    nwin =  size(windows,1);
    iwin = findwhichrange(impi(:,1),windows);
    % impiwin = cells(nwin,1);
    % implociwin = cells(nwin,1);
    for j = 1: nwin
      impiwin = impi(iwin==j,:);
      implociwin = imploci(iwin==j,:);
      lwinsec = round(((windows(j,2)-windows(j,1))+1)/sps);
      
      k=k+1;
      isub = isub+1;
      if isub > nrow1*ncol1
        %print the current figure f1
        f1name{ifig1,1} = sprintf('circloc_Fig%s_wlen%ds_%.1fkm.pdf',...
          num2zeropadstr(ifig1, 3),subwsec,diam);
        print(f1.fig,'-dpdf','-fillpage',strcat('/home/chaosong/Pictures/',f1name{ifig1,1}));
        %move on to next figure f1
        ifig1 = ifig1+1;  %initialize another figure
        f1 = initfig(widin1,htin1,nrow1,ncol1,ifig1);
        axpos1 = optaxpos(f1,nrow1,ncol1,pltxran1,pltyran1,pltxsep1,pltysep1);
        axtit = f1.ax(1: ncol1);
        supertit(axtit, sprintf('Fig %s, wlen=%d s, diam=%.1f km',...
          num2zeropadstr(ifig1, 3),subwsec,diam),10);
        xlabel(f1.ax((nrow1-1)*ncol1+1),'E (km)','fontsize',10);
        ylabel(f1.ax((nrow1-1)*ncol1+1),'N (km)','fontsize',10);
        orient(f1.fig,'landscape');
        
        isub = isub-nrow1*ncol1;  %refresh
      end
      
      %in terms of origin time
      % if strcmp(timetype,'tori')
        [tevt,impiwinst,indsort]=tarvl2tori(impiwin,sps,ftrans,1);  %return the origin time 
        implociwinst = implociwin(indsort, :);
      % end
      
      %plot locations of events in each cluster
      nevtsiwin = size(impiwin,1);
      color=jet(nevtsiwin);
      ax=f1.ax(isub); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
      axis(ax, 'equal');
      % scatter(ax,implociwinst(:,1),implociwinst(:,2),5,'k','filled');
      for jj=1: nevtsiwin
        [xcut,ycut] = circle_chao(implociwinst(jj,1),implociwinst(jj,2),radi,0.1);
        plot(ax,xcut,ycut,'-','color',color(jj,:),'linew',1);
      end
      x0=mean(implociwinst(:,1));
      y0=mean(implociwinst(:,2));
      plot(ax,x0+[0,pca2vec(k,1)],y0+[0,pca2vec(k,2)],'k-','linew',1);
      plot(ax,x0-[0,pca2vec(k,1)],y0-[0,pca2vec(k,2)],'k-','linew',1);
      arearat(k,ii) = area_of_overlap_circs_grid(implociwin(:,1:2),radi,0.01,0);
      text(ax,0.02,0.95,sprintf('%d events in %d s',nevtsiwin,lwinsec),...
        'Units','normalized','HorizontalAlignment','left','fontsize',7);
      text(ax,0.02,0.9,sprintf('amp:%.2f; %.1fx',ampi(k,1),ampi(k,1)/ampch),...
        'Units','normalized','HorizontalAlignment','left','fontsize',7);
      text(ax,0.02,0.85,sprintf('arearat:%.2f',arearat(k,ii)),...
        'Units','normalized','HorizontalAlignment','left','fontsize',7);
      text(ax,0.98,0.15,sprintf('%.2f; %.2f km',mdistnn1(k,1),mdist2all(k,1)),...
        'Units','normalized','HorizontalAlignment','right','fontsize',7);
      text(ax,0.98,0.1,sprintf('%d / %d ^o; %.2f; %.2f km',round(projang1(k,1)),...
        round(projang1(k,1))+180,mprojx1nn1(k,1),mprojx12all(k,1)),...
        'Units','normalized','HorizontalAlignment','right','fontsize',7);
      text(ax,0.98,0.05,sprintf('%d ^o; %.2f; %.2f km',round(projang2(k,1)),...
        mprojx2nn1(k,1),mprojx22all(k,1)),...
        'Units','normalized','HorizontalAlignment','right','fontsize',7);
      xlim(ax,[-4 4]);
      ylim(ax,[-4 4]);
      longticks(ax,2);
      hold(ax,'off');
      %       keyboard
    end
  end
  % arearat{ii}=arearati;
  
  %if all clusters have been plotted, print the current figure f1
  f1name{ifig1,1} = sprintf('circloc_Fig%s_wlen%ds_%.1fkm.pdf',...
    num2zeropadstr(ifig1, 3),subwsec,diam);
  print(f1.fig,'-dpdf','-fillpage',strcat('/home/chaosong/Pictures/',f1name{ifig1,1}));
  
  %merge all figures into a single pdf file
  mergef1name=sprintf('circloc_wlen%ds_%.1fkm.pdf',subwsec,diam);
  str=strcat('rm -f /home/chaosong/Pictures/',mergef1name);
  status = system(str);
  % status = system('rm -f /home/chaosong/Pictures/circloc.pdf');
  for i = 1:size(f1name,1)
    fname{i} = strcat('/home/chaosong/Pictures/',f1name{i});
  end
  append_pdfs(strcat('/home/chaosong/Pictures/',mergef1name),fname);
  status = system('rm -f /home/chaosong/Pictures/circloc_Fig*.pdf');
  % keyboard
  close all
  
end

%%
%plot the distribution of total covered area normalized by norminal area
f = initfig(8,8,2,2); %initialize fig
binw=0.02;
binedge=0:binw:1;
for ii = 1:length(diams)
  diam=diams(ii);
  % arearati=arearat{ii};
  ax=f.ax(ii); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  histogram(ax,arearat(:,ii),'BinEdges',binedge);
  text(ax,0.98,0.95,sprintf('diameter: %.1f km',diam),'Units','normalized',...
    'HorizontalAlignment','right');
  xlabel(ax,'Normalized area');
  ylabel(ax,'Count');
  xlim(ax,[0 1]);
end
supertit(f.ax(1:2),sprintf('wlen=%d s, nwin=%d',subwsec,k),10);
% keyboard
figname = sprintf('arearat_wlen%ds_nwin%d.pdf',subwsec,k);
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',figname));

