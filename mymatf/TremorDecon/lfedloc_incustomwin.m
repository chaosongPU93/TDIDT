% lfedloc_incustomwin.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Separated from 'circsrcmodel2.m', now this becomes the script that specifically
% compute the distance or location difference between source pair N and N-n
% in the context of a custom window with a certain length. For example, you can
% ask what is distance between consecutive events among all events included in
% the window of a certain length. Results from each window are also lumped
% together to get a reliable statistics. You can also look at individual
% windows, but that would be the main purpose of 'circsrcmodel2.m'.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/03/06
% Last modified date:   2024/03/06
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
off1i = allbstsig.off1i;  %stores alignment of the whole win
trangenew = allbstsig.trangenew;
% windowsk = allbstsig.windowsk;  %stores start and end indices of 25-s OVERLAPPING windows
irccrank = allbstsig.irccrank;  %stores start and end indices of RCC
off1iwk = allbstsig.off1iwk; %stores alignment of each 25-s win
locxyprojalln = allbstnoi.locxyprojall;
tarvlsplstalln = allbstnoi.impindepall(:,1);
nsrcn = allbstnoi.nsrc;
impn = allbstnoi.impindepall;
off1in = allbstnoi.off1i;
off1iwkn = allbstnoi.off1iwk;
supertstr = 'Secondary sources removed';
fnsuffix = [];

% %%%param for further checked at KLNB
% locxyprojall = allbstsig.locxyproj4thall;
% tarvlsplstall = allbstsig.impindep4thall(:,1);
% nsrc = allbstsig.nsrc4th;
% imp = allbstsig.impindep4thall;
% off1i = allbstsig.off1i;
% trangenew = allbstsig.trangenew;
% irccrank = allbstsig.irccrank;  %stores start and end indices of RCC
% off1iwk = allbstsig.off1iwk; %stores alignment of each 25-s win
% locxyprojalln = allbstnoi.locxyproj4thall;
% tarvlsplstalln = allbstnoi.impindep4thall(:,1);
% nsrcn = allbstnoi.nsrc4th;
% impn = allbstnoi.impindep4thall;
% off1in = allbstnoi.off1i;
% supertstr = 'Further checked at KLNB';
% fnsuffix = '4sta';

% keyboard

% impuse = imp;
% nsrcuse = nsrc;
% fnsuffix2 = [];
impuse = impn;
nsrcuse = nsrcn;
fnsuffix2 = 'noi';

[imploc, ~] = off2space002(impuse(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[imploc0, ~] = off2space002([0 0],sps,ftrans,0);  % a ref source at 0,0
amp=mean(impuse(:,[2 4 6]),2);
ampch=prctile(amp,5);

tlennew = trangenew(:,3)-trangenew(:,2);
tlensumnew = sum(tlennew);

timetype = 'tarvl';

%% distance between events, including PCA
% subwsectar = 25/2;  %target subwindow length in sec
subwsectar = 25;  %target subwindow length in sec
n = 1;  %between N and N-n
m = 1;  %max n to compute

%%%Flag to indicate if it is necessary to recalculate everything
flagrecalc = 0;
% flagrecalc = 1;

if isequaln(impuse,imp)
  savefile = strcat('lfedloc',fnsuffix,num2str(subwsectar),'swin.mat');
elseif isequaln(impuse,impn)
  savefile = strcat('lfedloc',fnsuffix,num2str(subwsectar),'swin_noi.mat');
end

if flagrecalc
  
  k = 0;  % count of the total subplots
  dtnn1cat=[];  % time difference between N and N-1, always + 
  dlocnn1cat=[];  % loc difference between N and N-1, sign preserved
  dloc_splnn1cat=[];  % loc difference between N and N-1, sample space, sign preserved
  dt2allcat=[]; % time difference between each to all others, always + 
  dloc2allcat=[]; % loc difference between each to all others, sign preserved
  dloc2all_splcat=[]; % loc difference between each to all others, sample space, sign preserved
  distnn1cat=[];  % abs distance between N and N-1
  dist2allcat=[]; % abs distance between each to all others
  dprojxy1nn1cat=[];  % loc diff between N and N-1 projected along PCA, sign preserved
  dprojxy12allcat=[]; % loc diff between each to others projected along PCA, sign preserved
  dprojxy2nn1cat=[];  % loc diff between N and N-1 projected along prop direct, sign preserved
  dprojxy22allcat=[]; % loc diff between each to others projected along prop direct, sign preserved
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
    if subwsectar == 25
      windows = irccrank{i};
    else
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
    end
    
    nwin =  size(windows,1);
    iwin = findwhichrange(impi(:,1),windows);
%     keyboard
    for j = 1: nwin
      impiwin = impi(iwin==j,:);
      implociwin = imploci(iwin==j,:);
      lwinsec = round(((windows(j,2)-windows(j,1))+1)/sps); %notes the actual subwin length in sec

      if size(impiwin,1) >1
        k=k+1;

        %%%%%% diff loc within short win, able to be combined to analyse later
        %compute the diff loc between N and N-m, and each to all others in the
        %short window, no projection is applied
        [dloc,dt,dloc_spl,dloc2all,dt2all,dloc2all_spl]=...
          dloc_evtcustom(impiwin,implociwin,sps,ftrans,m,'tarvl');
  %       mdtnn1(k,1)=median(dt{n});
  %       mdlocnn1(k,:)=median(dloc{n});
  %       mdloc_splnn1(k,:)=median(dloc_spl{n});
  %       mdloc2all(k,:)=median(dloc2all);
  %       mdt2all(k,1)=median(dt2all);
  %       mdloc2all_spl(k,:)=median(dloc2all_spl);
        dtnn1cat=[dtnn1cat; dt{n}];
        dlocnn1cat=[dlocnn1cat; dloc{n}];
        dloc_splnn1cat=[dloc_splnn1cat; dloc_spl{n}];
        dloc2allcat=[dloc2allcat; dloc2all];
        dt2allcat=[dt2allcat; dt2all];
        dloc2all_splcat=[dloc2all_splcat; dloc2all_spl];
        %%%%%%%%%%

        %%%%%% distance, etc for each short win
        %compute the abs distance between N and N-m, and each to all others in the
        %short window. The dist is the sqrt of 'dloc'
        %For all srcs in each win, apply the PCA, and project to 
        %this direction. Note that this direction is not the same as to all srcs
        %combined
        [distnn1,dist2all,pcavec,projang1,dprojxy1nn1,dprojxy12all,...
          projang2,dprojxy2nn1,dprojxy22all]=dist_evtcustom(impiwin,implociwin,sps,ftrans,'tarvl'); %tori
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
  
  save(strcat(rstpath, '/MAPS/',savefile), 'dtnn1cat','dlocnn1cat','dloc_splnn1cat','dloc2allcat',...
    'dt2allcat','dloc2all_splcat','distnn1cat','dist2allcat','dprojxy1nn1cat','dprojxy12allcat',...
    'dprojxy2nn1cat','dprojxy22allcat');
  
else
%   savefile = strcat('lfedloc',fnsuffix,num2str(subwsectar),'swin.mat');
  load(strcat(rstpath, '/MAPS/',savefile));
end

% keyboard

%% type of location difference or distance to look at, between consec. ones, or each to others
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

%%
smoothsigma=5;
ncont=100;  %num of contour lines
dx=0.025; dy=0.025; % to distinguish dp near origin, size cannot exceed ~0.053 km 
cstr={'# events / grid'}; xran=[-4 4]; yran=[-4 4];
[f,den1d,conmat]=plt_srcdloc(dloc,'km',3,cstr,...
  'o','linear','grid',xran,yran,dx,dy,smoothsigma,ncont);
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc);
[~,~,projxy] = customprojection(dloc(:,1:2),anglegeo(2));
mdistprojopt = median(abs(projxy(:,1)))
medprojopt = median(projxy(:,1))
mdistprojort = median(abs(projxy(:,2)))
medprojort = median(projxy(:,2)) 
  
fname = strcat('lfedloc2all',num2str(subwsectar),'swin',fnsuffix,fnsuffix2,'.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

keyboard

%% summarize the location diff in sample space 
%plot the diff location in samples, bin by pixel
smoothsigma=1;  %smoothing sigma for Gaussian filtering 
ncont=100;  %num of contour lines
% dx=1; dy=1;
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc_spl);
cstr={'# events / pixel'}; xran=[-50 50]; yran=[-50 50];
[f,den1d_spl,conmat]=plt_srcdlocinmap(dt/sps,dloc_spl,[],'spl',timetype,...
  8,cstr,'o','linear','pixel',xran,yran,[],[],smoothsigma,ncont);
ax=f.ax(2);
hold(ax,'on');
plot(ax,ellx,elly,'k-','linew',1);
% plot(ax,x0+10*[coeff(1,2),-coeff(1,2)],y0+10*[coeff(2,2),-coeff(2,2)],'k--','linew',2);
% plot(ax,x0+10*[coeff(1,1),-coeff(1,1)],y0+10*[coeff(2,1),-coeff(2,1)],'k--','linew',2);
text(ax,0.98,0.15,sprintf('%d ^o',round(anglegeo(2))),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.98,0.1,sprintf('asprat: %.1f',semia/semib),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
% supertit(f.ax, strcat({'location difference '},str));
hold(ax,'off');
ax=f.ax(3);
hold(ax,'on');
%interpolation to obtain the intersection between PCA directions and contours
F = scatteredInterpolant(conmat(:,3),conmat(:,4),conmat(:,1),'linear','none');
%a line cross (0,0) in the projection direction
x = reshape(xran(1):0.05:xran(2), [], 1);
yopt = linefcn(x,tand(angle(2)),0);
plot(ax,x,yopt,'k-','linew',2);
%a line cross (0,0) in the orthogonal direction
yort = linefcn(x,tand(angle(1)),0);
plot(ax,x,yort,'-','linew',2,'color',[.5 .5 .5]);
text(ax,0.95,0.1,sprintf('smooth sigma: %.1f',smoothsigma),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.95,0.05,sprintf('%d contours',ncont),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);

%if choosing a line intersecting with the contour
normalizer=length(dloc);
f = initfig(8,4.5,1,2);
optaxpos(f,1,2);
[f.ax(1),muopt,sigmaopt,mdistprojopt]=plt_dloccrssect(f.ax(1),F,x,yopt,anglegeo(2),...
  [0 0 0; 0 0 1],xran,normalizer); %,yopt1,yopt2
[f.ax(2),muort,sigmaort,mdistprojort]=plt_dloccrssect(f.ax(2),F,x,yort,anglegeo(1),...
  [.5 .5 .5; 1 0 0],xran,normalizer); %,yort1,yort2

% [locall, indinput] = off2space002([],sps,ftrans); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
% Fx = scatteredInterpolant(locall(:,7),locall(:,8),locall(:,1),'linear','none');
% Fy = scatteredInterpolant(locall(:,7),locall(:,8),locall(:,2),'linear','none');
% 
% grp = unique(conmat(:,2));
% ngrp = max(grp);
% conmatxy = conmat;
% for igrp = 1:ngrp
%   ind = find(conmat(:,2)==igrp);  
%   conmatxy(ind,3) = Fx(conmat(ind,3),conmat(ind,4));
%   conmatxy(ind,4) = Fy(conmat(ind,3),conmat(ind,4));
% end
% % maxcont = max(conmatxy(:,1));
% % contintvl = conmatxy(1,1);
% contval=unique(conmat(:,1));
% % ncont = round(conmatxy(end,1)/conmatxy(1,1));
% color = jet(ncont);
% xran=[-4 4]; yran=[-4 4];
% f=initfig;
% ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
% for icont = 1:ncont
%   for igrp = 1:ngrp
%     ind = find(conmat(:,1)==contval(icont) & conmat(:,2)==igrp);
%     plot(ax,conmatxy(ind,3),conmatxy(ind,4),'-','Color',color(icont,:));
%   end
% end
% colormap(ax,'jet');
% c=colorbar(ax,'SouthOutside');
% ax.CLim(2) = conmatxy(end,1);
% c.Label.String = cstr{1}; 
% F = scatteredInterpolant(conmatxy(:,3),conmatxy(:,4),conmatxy(:,1),'linear','none');
% %using the direction from absolute locations to project
% [~,~,angle,anglegeo]=pcaellipse(dloc);
% %a line cross (0,0) in the projection direction
% x = reshape(xran(1):0.01:xran(2), [], 1);
% yopt = linefcn(x,tand(angle(2)),0);
% % yopt1 = linefcn(x,tand(angle(2)),0.25/cosd(angle(2)));
% % yopt2 = linefcn(x,tand(angle(2)),-0.25/cosd(angle(2)));
% plot(ax,x,yopt,'k-','linew',2);
% % plot(ax,x,yopt1,'k--','linew',1);
% % plot(ax,x,yopt2,'k-.','linew',1);
% %a line cross (0,0) in the orthogonal direction
% yort = linefcn(x,tand(angle(1)),0);
% % yort1 = linefcn(x,tand(angle(1)),-0.25/cosd(angle(1)));
% % yort2 = linefcn(x,tand(angle(1)),0.25/cosd(angle(1)));
% plot(ax,x,yort,'-','linew',2,'color',[.5 .5 .5]);
% % plot(ax,x,yort1,'--','linew',1,'color',[.5 .5 .5]);
% % plot(ax,x,yort2,'-.','linew',1,'color',[.5 .5 .5]);
% ax.GridLineStyle = '--';
% ax.XAxisLocation = 'top';
% xlim(ax,xran);
% ylim(ax,yran);
% xlabel(ax,'Diff E loc (km)');
% ylabel(ax,'Diff N loc (km)');
% longticks(ax,2);
% hold(ax,'off');
% 
% normalizer=length(dloc);
% f = initfig(8,4.5,1,2);
% optaxpos(f,1,2);
% [f.ax(1),muopt,sigmaopt,mdistprojopt]=plt_dloccrssect(f.ax(1),F,x,yopt,anglegeo(2),...
%   [0 0 0; 0 0 1],xran,normalizer); %,yopt1,yopt2
% 
% [f.ax(2),muort,sigmaort,mdistprojort]=plt_dloccrssect(f.ax(2),F,x,yort,anglegeo(1),...
%   [.5 .5 .5; 1 0 0],xran,normalizer); %,yort1,yort2


%% summarize the location diff in map space
% smoothsigmap = 4:0.5:5;
% for i = 1: length(smoothsigmap)
%     smoothsigma=smoothsigmap(i);
%plot the diff location in map, bin by grid
smoothsigma=5;  %sigma of Gaussian filter for smoothing 
ncont=100;  %contour interval
dx=0.025; dy=0.025;
%PCA analysis
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc);
cstr={'# events / grid'}; xran=[-4 4]; yran=[-4 4];
[f,den1d,conmat]=plt_srcdlocinmap(dt/sps,dloc,[],'km',timetype,...
  3,cstr,'o','linear','grid',xran,yran,dx,dy,smoothsigma,ncont);
fname = strcat('lfedloc2all',num2str(subwsectar),'swin',fnsuffix,fnsuffix2,'.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

%%
%%%histogram of overall diff loc distribution and the projection
nrow = 1; % rows and cols of subplots in each figure
ncol = 2; 
widin = 6; % size of each figure
htin = 3.5;
pltxran = [0.1 0.96]; pltyran = [0.12 0.9];
pltxsep = 0.07; pltysep = 0.03;
f = initfig(widin,htin,nrow,ncol);
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
supertit(f.ax(1:2),sprintf('All, %d dpts',size(dloc,1)),10);
xran=[-4 4]; binw=0.1; legendstr1={'E','N'};
normopt='custom'; normalizer=length(dloc);
f.ax(1)=plt_dlochist(f.ax(1),dloc,xran,binw,legendstr1,normopt,normalizer);

legendstr2={sprintf('%d ^o',round(anglegeo(2))), ...
  sprintf('%d ^o',round(anglegeo(1)))};
[~,~,dprojloc] = customprojection(dloc(:,1:2),anglegeo(2));
f.ax(2)=plt_dlochist(f.ax(2),dprojloc,xran,binw,legendstr2,normopt,normalizer);
xlabel(f.ax(2),'');
ylabel(f.ax(2),'');
distprojloc=median(abs(dprojloc(:,1)));

%%
%%%if choosing a line intersecting with the contour
F = scatteredInterpolant(conmat(:,3),conmat(:,4),conmat(:,1),'linear','none');
x = reshape(xran(1):0.001:xran(2), [], 1);
yopt = linefcn(x,tand(angle(2)),0); %a line cross (0,0) in the projection direction
yort = linefcn(x,tand(angle(1)),0);%a line cross (0,0) in the orthogonal direction
normalizer=length(dloc);
nrow = 1; % rows and cols of subplots in each figure
ncol = 2; 
widin = 6; % size of each figure
htin = 3.5;
pltxran = [0.1 0.96]; pltyran = [0.12 0.9];
pltxsep = 0.07; pltysep = 0.03;
f = initfig(widin,htin,nrow,ncol);
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

[f.ax(1),muopt,sigmaopt,mdistprojopt]=plt_dloccrssect(f.ax(1),F,x,yopt,anglegeo(2),...
  [0 0 0; 0 0 1],xran,normalizer);%,yopt1,yopt2

[f.ax(2),muort,sigmaort,mdistprojort]=plt_dloccrssect(f.ax(2),F,x,yort,anglegeo(1),...
  [.5 .5 .5; 1 0 0],xran,normalizer);%,yort1,yort2
xlabel(f.ax(2),'');
ylabel(f.ax(2),'');

fname = strcat('lfedlocproj2all',num2str(subwsectar),'swin',fnsuffix,fnsuffix2,'.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

keyboard

%% demo of generating random samples from the 2D Gaussian from fitting contours 
%plot the diff location in samples, bin by pixel
smoothsigma=1;
ncont=100;
dx=1; dy=1;
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc_spl);
cstr={'# events / pixel'}; xran=[-50 50]; yran=[-50 50];
[f,den1d_spl,conmat,pks2d]=plt_srcdlocinmap(dt/sps,dloc_spl,[],'spl',timetype,...
  8,cstr,'o','linear','pixel',xran,yran,dx,dy,1,smoothsigma,ncont);
ax=f.ax(2);
hold(ax,'on');
plot(ax,ellx,elly,'k-','linew',2);
% plot(ax,x0+10*[coeff(1,2),-coeff(1,2)],y0+10*[coeff(2,2),-coeff(2,2)],'k--','linew',2);
% plot(ax,x0+10*[coeff(1,1),-coeff(1,1)],y0+10*[coeff(2,1),-coeff(2,1)],'k--','linew',2);
text(ax,0.98,0.15,sprintf('%d ^o',round(anglegeo(2))),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.98,0.1,sprintf('asprat: %.1f',semia/semib),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
supertit(f.ax, strcat({'location difference '},str));
hold(ax,'off');
ax=f.ax(3);
hold(ax,'on');
F = scatteredInterpolant(conmat(:,3),conmat(:,4),conmat(:,1),'linear','none');
%a line cross (0,0) in the projection direction
x = reshape(xran(1):0.05:xran(2), [], 1);
yopt = linefcn(x,tand(angle(2)),0);
plot(ax,x,yopt,'k-','linew',2);
%a line cross (0,0) in the orthogonal direction
yort = linefcn(x,tand(angle(1)),0);
plot(ax,x,yort,'-','linew',2,'color',[.5 .5 .5]);
text(ax,0.75,0.1,sprintf('smooth sigma: %.1f',smoothsigma),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.75,0.05,sprintf('%d contours',ncont),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
hold(ax,'off');

npks2d = size(pks2d,1); %must be 2
signx = pks2d(1,1)*pks2d(2,1);  %must be negative
signy = pks2d(1,2)*pks2d(2,2);  %must be negative
dx=diff(pks2d(:,1));
dy=diff(pks2d(:,2));
theta=atan2d(dy,dx);  
if theta < 0 
  theta=180+theta;   %must be very close to angle(1)
end
distpk = sqrt(dx^2+dy^2); %must be >=4

% %% if choosing a line intersecting with the contour
normalizer=length(dloc);
f = initfig(8,4.5,1,2);
optaxpos(f,1,2);
[f.ax(1),muopt,sigmaopt,mdistprojopt,pksopt]=plt_dloccrssect(f.ax(1),F,x,yopt,anglegeo(2),...
  [0 0 0; 0 0 1],xran,normalizer); %,yopt1,yopt2

[f.ax(2),muort,sigmaort,mdistprojort,pksort]=plt_dloccrssect(f.ax(2),F,x,yort,anglegeo(1),...
  [.5 .5 .5; 1 0 0],xran,normalizer); %,yort1,yort2
%only focus on peaks near the origin
pksort=pksort(abs(pksort(:,1))<=5, :);
%to be similar to the case like data, there are 3 peaks, one close to 0,
%other 2 roughly symmetric to origin
npks = size(pksort,1);

%%%%%%%%%%%%% random sampling from 2D Gaussian
% xvec = xran(1):dx:xran(2);
% yvec = yran(1):dy:yran(2);
% [X1,X2] = meshgrid(xvec,yvec);
% X = [X1(:) X2(:)];
% mu = [0 0];
% sigmax = sigmaort;
% sigmay = sigmaopt;
% rotang = angle(1);
% covar = mvncov(sigmax,sigmay,rotang);
% [eigvec,eigval] = eig(covar);
% sigmaeig = [sqrt(eigval(1,1)) sqrt(eigval(2,2))]
% pdf2 = mvnpdf(X,mu,covar);
% figure
% imagesc(xvec,yvec,reshape(pdf2,length(yvec),length(xvec))); hold on;% substitute for surf, but need to make 'y' axis direction as normal
% ax = gca; ax.YDir = 'normal';
% axis equal tight
% colormap('jet');
% xlabel('x1');
% ylabel('x2');
% c=colorbar(ax,'SouthOutside');
% c.Label.String = 'Probability Density';
% R = mvnrnd(mu,covar,size(dloc_spl,1));
% scatter(R(:,1),R(:,2),10,[.5 .5 .5],'+');
% dloc_spluse = R;
%%%%%%%%%%%%% random sampling from 2D Gaussian

%%%%%%%%%%%%% Bootstrapping 
[~,ind] = bootstrp(1,[],dloc_spl);
dloc_spluse = dloc_spl(ind,:);
%%%%%%%%%%%%% Bootstrapping 

smoothsigma=1;
ncont=100;
dx=1; dy=1;
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc_spluse);
cstr={'# events / grid'}; xran=[-50 50]; yran=[-50 50];
[f,den1d_spl,conmat,pks2ddemo]=...
  plt_srcdlocinmap(ones(size(dloc_spluse,1),1)/sps,dloc_spluse,[],'spl',timetype,...
  8,cstr,'o','linear','grid',xran,yran,dx,dy,1,smoothsigma,ncont);
ax=f.ax(2);
hold(ax,'on');
plot(ax,ellx,elly,'k-','linew',2);
% plot(ax,x0+10*[coeff(1,2),-coeff(1,2)],y0+10*[coeff(2,2),-coeff(2,2)],'k--','linew',2);
% plot(ax,x0+10*[coeff(1,1),-coeff(1,1)],y0+10*[coeff(2,1),-coeff(2,1)],'k--','linew',2);
text(ax,0.98,0.15,sprintf('%d ^o',round(anglegeo(2))),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.98,0.1,sprintf('asprat: %.1f',semia/semib),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
supertit(f.ax, strcat({'location difference '},str));
hold(ax,'off');
ax=f.ax(3);
hold(ax,'on');
F = scatteredInterpolant(conmat(:,3),conmat(:,4),conmat(:,1),'linear','none');
%a line cross (0,0) in the projection direction
x = reshape(xran(1):0.05:xran(2), [], 1);
yopt = linefcn(x,tand(angle(2)),0);
plot(ax,x,yopt,'k-','linew',2);
%a line cross (0,0) in the orthogonal direction
yort = linefcn(x,tand(angle(1)),0);
plot(ax,x,yort,'-','linew',2,'color',[.5 .5 .5]);
text(ax,0.75,0.1,sprintf('smooth sigma: %.1f',smoothsigma),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.75,0.05,sprintf('%d contours',ncont),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);

% npks2ddemo = size(pks2ddemo,1); %must be 2
% signxdemo = pks2ddemo(1,1)*pks2ddemo(2,1);  %must be negative
% signydemo = pks2ddemo(1,2)*pks2ddemo(2,2);  %must be negative
% dx=diff(pks2ddemo(:,1));
% dy=diff(pks2ddemo(:,2));
% thetademo=atan2d(dy,dx);  
% if thetademo < 0 
%   thetademo=180+thetademo;   %must be very close to angle(1)
% end
% distpkdemo = sqrt(dx^2+dy^2); %must be >=4

% %% if choosing a line intersecting with the contour
normalizer=length(dloc);
f = initfig(8,4.5,1,2);
optaxpos(f,1,2);
[f.ax(1),muoptrnddemo,sigmaoptrnddemo,~,pksoptdemo]=plt_dloccrssect(f.ax(1),F,x,yopt,anglegeo(2),...
  [0 0 0; 0 0 1],xran,normalizer); %,yopt1,yopt2

[f.ax(2),muortrnddemo,sigmaortrnddemo,~,pksortdemo]=plt_dloccrssect(f.ax(2),F,x,yort,anglegeo(1),...
  [.5 .5 .5; 1 0 0],xran,normalizer); %,yort1,yort2


%% generating random samples from the 2D Gaussian from fitted contours
% xran=[-50 50]; yran=[-50 50];
% dx=1; dy=1;
% xvec = xran(1):dx:xran(2);
% yvec = yran(1):dy:yran(2);
% [X1,X2] = meshgrid(xvec,yvec);
% X = [X1(:) X2(:)];
% mu = [0 0];
% sigmax = sigmaort;
% sigmay = sigmaopt;
% rotang = angle(1);
% covar = mvncov(sigmax,sigmay,rotang);
% [eigvec,eigval] = eig(covar);
% sigmaeig = [sqrt(eigval(1,1)) sqrt(eigval(2,2))]
% 
% ntrial = 10000;
% 
% flagrecalc = 0;
% % flagrecalc = 1;
% if flagrecalc
% 
%   rng('default');
% 
%   lnoptdloc=cell(ntrial,1);
%   lnoptcount=cell(ntrial,1);
%   lnortdloc=cell(ntrial,1);
%   lnortcount=cell(ntrial,1);
%   pks2drnd=cell(ntrial,1);
%   pksoptrnd=cell(ntrial,1);
%   pksortrnd=cell(ntrial,1);
%   anglernd=cell(ntrial,1);
%   anglegeornd=cell(ntrial,1);
%   
%   % f = initfig(8,4.5,1,2);
%   % optaxpos(f,1,2);
%   for i = 1: ntrial
%     R = mvnrnd(mu,covar,size(dloc_spl,1));
%     [coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(R);
%     den1d = density_matrix(R(:,1),R(:,2),xran,yran,dx,dy);
%     den1d = den1d(den1d(:,3)>0, :);
%     den1d = sortrows(den1d,3);
%     [~,xgrid,ygrid,zgrid] = ...
%       zeropadmat2d(den1d,min(den1d(:,1)):1:max(den1d(:,1)),...
%       min(den1d(:,2)):1:max(den1d(:,2)));
%     zgridgf = imgaussfilt(zgrid,smoothsigma);  %smooth it a bit
%     anglernd{i} = angle;
%     anglegeornd{i} = anglegeo;
%     
%     %find local maxima near the origin of the 2D data 
% %     cent=FastPeakFind(zgrid, 50, 1);
%     [pkshgt,locs_y,locs_x]=peaks2(zgridgf);
%     pks=[reshape(xgrid(1,locs_x),[],1) ...
%          reshape(ygrid(locs_y,1),[],1) ...
%          reshape(pkshgt,[],1)];
% %     ind=find(abs(pksxy(:,1))<=5 & abs(pksxy(:,2))<=5);
%     %only focus on peaks near the origin
%     pks = pks(abs(pks(:,1))<=5 & abs(pks(:,2))<=5, :);
%     pks2drnd{i} = pks;
% 
%     f2=figure; ax=gca; axis(ax, 'equal');
%     [conmat,conobj] = contour(ax,xgrid,ygrid,zgridgf,100,'-');
%     close(f2);
%     contable = getContourLineCoordinates(conmat);
%     conmat=table2array(contable);
% 
%     F = scatteredInterpolant(conmat(:,3),conmat(:,4),conmat(:,1),'linear','none');
%     %a line cross (0,0) in the projection direction
%     x = reshape(xran(1):0.05:xran(2), [], 1);
%     yopt = linefcn(x,tand(angle(2)),0);
%     %a line cross (0,0) in the orthogonal direction
%     yort = linefcn(x,tand(angle(1)),0);
% 
%   %   ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%     crssecx = x; crssecy = yopt; angle = anglegeo(2); color=[0 0 0 0.3; 0 0 1 0.15];
%     count = F(crssecx,crssecy);
%     dprojx = customprojection([crssecx crssecy],angle);
%     dprojx = dprojx(~isnan(count));
%     count = count(~isnan(count));
%     lnoptdloc{i} = dprojx;
%     lnoptcount{i} = count;
%     fitobj = fit(dprojx,count,'gauss1','StartPoint',[10 0 1]);
%     Y1=feval(fitobj,dprojx);
%     coef = coeffvalues(fitobj);
%     muoptrnd(i,1) = coef(2);
%     sigmaoptrnd(i,1) = coef(3);
%     [pkhgt,locs] = findpeaks(count,dprojx);
%     pks=[reshape(locs,[],1) reshape(pkhgt,[],1)];
%     %only focus on peaks near the origin
%     pks = pks(abs(pks(:,1))<=5, :);
%     pksoptrnd{i} = pks;
%   %   plot(ax,dprojx,count / normalizer,'-','Color',color(1,:),'linew',1);
%   %   plot(ax,dprojx,Y1 / normalizer,'-','Color',color(2,:),'linew',1);
%   %   xlim(ax,xran);
%   %   hold(ax,'off');
% 
%   %   ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%     crssecx = x; crssecy = yort; angle = anglegeo(1); color=[1 0 0 0.3; 1 0 0 0.15];
%     count = F(crssecx,crssecy);
%     dprojx = customprojection([crssecx crssecy],angle);
%     dprojx = dprojx(~isnan(count));
%     count = count(~isnan(count));
%     lnortdloc{i} = dprojx;
%     lnortcount{i} = count;
%     fitobj = fit(dprojx,count,'gauss1','StartPoint',[10 0 1]);
%     Y1=feval(fitobj,dprojx);
%     coef = coeffvalues(fitobj);
%     muortrnd(i,1) = coef(2);
%     sigmaortrnd(i,1) = coef(3);
%     [pkhgt,locs] = findpeaks(count,dprojx);
%     pks=[reshape(locs,[],1) reshape(pkhgt,[],1)];
%     %only focus on peaks near the origin
%     pks = pks(abs(pks(:,1))<=5, :);
%     pksortrnd{i} = pks;
%   %   plot(ax,dprojx,count / normalizer,'-','Color',color(1,:),'linew',1);
%   %   plot(ax,dprojx,Y1 / normalizer,'-','Color',color(2,:),'linew',1);
%   %   xlim(ax,xran);
%   %   hold(ax,'off');
% 
%   end
%   savefile = strcat('randsamplfedloc',num2str(subwsectar),'swin.mat');
%   save(strcat(rstpath, '/MAPS/',savefile), 'lnoptdloc','lnoptcount','lnortdloc',...
%     'lnortcount','muoptrnd','sigmaoptrnd','muortrnd','sigmaortrnd','pks2drnd',...
%     'pksoptrnd','pksortrnd','anglernd','anglegeornd');
% 
% else
%   savefile = strcat('randsamplfedloc',num2str(subwsectar),'swin.mat');
%   load(strcat(rstpath, '/MAPS/',savefile));
% end


%% Bootstrapping test
ntrial = 10000;

flagrecalc = 0;
% flagrecalc = 1;
if flagrecalc

  %%%bootstrap random sampling with replacement, will return the SAME size
  [~,indboot] = bootstrp(ntrial,[],dloc_spl);
  %   %%%random sampling with no replaceent, will return desired size
  %   for i = 1: ntrial
  %     indboot(:,i) = randsample(size(dloc_spl,1), round(0.8*size(dloc_spl,1)));
  %   end
  
  lnoptdlocboot=cell(ntrial,1);
  lnoptcountboot=cell(ntrial,1);
  lnortdlocboot=cell(ntrial,1);
  lnortcountboot=cell(ntrial,1);
  pks2dboot=cell(ntrial,1);
  pksoptboot=cell(ntrial,1);
  pksortboot=cell(ntrial,1);
  angleboot=cell(ntrial,1);
  anglegeoboot=cell(ntrial,1);

  for i = 1: ntrial
    
    dloc_spluse = dloc_spl(indboot(:,i),:);
%     dtuse = dt(indboot,:);
    
    [coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc_spluse);    
    den1d = density_pixel(dloc_spluse(:,1),dloc_spluse(:,2));
    den1d = den1d(den1d(:,3)>0, :);
    den1d = sortrows(den1d,3);    
    [~,xgrid,ygrid,zgrid] = ...
      zeropadmat2d(den1d,min(den1d(:,1)):1:max(den1d(:,1)),...
      min(den1d(:,2)):1:max(den1d(:,2)));
    smoothsigma=1;
    zgridgf = imgaussfilt(zgrid,smoothsigma);  %smooth it a bit
    angleboot{i} = angle;
    anglegeoboot{i} = anglegeo;
    
    %find local maxima near the origin of the 2D data 
%     cent=FastPeakFind(zgrid, 50, 1);
    [pkshgt,locs_y,locs_x]=peaks2(zgridgf);
    pks=[reshape(xgrid(1,locs_x),[],1) ...
         reshape(ygrid(locs_y,1),[],1) ...
         reshape(pkshgt,[],1)];
%     ind=find(abs(pksxy(:,1))<=5 & abs(pksxy(:,2))<=5);
    %only focus on peaks near the origin
    pks = pks(abs(pks(:,1))<=5 & abs(pks(:,2))<=5, :);
    pks2dboot{i} = pks;
    
    f2=figure; ax=gca; axis(ax, 'equal');
    [conmat,conobj] = contour(ax,xgrid,ygrid,zgridgf,100,'-');
    close(f2);
    contable = getContourLineCoordinates(conmat);
    conmat=table2array(contable);
    
    F = scatteredInterpolant(conmat(:,3),conmat(:,4),conmat(:,1),'linear','none');
    %a line cross (0,0) in the projection direction
    x = reshape(xran(1):0.05:xran(2), [], 1);
    yopt = linefcn(x,tand(angle(2)),0);
    %a line cross (0,0) in the orthogonal direction
    yort = linefcn(x,tand(angle(1)),0);
    
    %   ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    crssecx = x; crssecy = yopt; angle = anglegeo(2); color=[0 0 0 0.3; 0 0 1 0.15];
    count = F(crssecx,crssecy);
    dprojx = customprojection([crssecx crssecy],angle);
    dprojx = dprojx(~isnan(count));
    count = count(~isnan(count));
    lnoptdlocboot{i} = dprojx;
    lnoptcountboot{i} = count;
    fitobj = fit(dprojx,count,'gauss1','StartPoint',[10 0 1]);
    Y1=feval(fitobj,dprojx);
    coef = coeffvalues(fitobj);
    muoptboot(i) = coef(2);
    sigmaoptboot(i) = coef(3);
    [pkhgt,locs] = findpeaks(count,dprojx);
    pks=[reshape(locs,[],1) reshape(pkhgt,[],1)];
    %only focus on peaks near the origin
    pks = pks(abs(pks(:,1))<=5, :);
    pksoptboot{i} = pks;
    %   plot(ax,dprojx,count / normalizer,'-','Color',color(1,:),'linew',1);
    %   plot(ax,dprojx,Y1 / normalizer,'-','Color',color(2,:),'linew',1);
    %   xlim(ax,xran);
    %   hold(ax,'off');
    
    %   ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    crssecx = x; crssecy = yort; angle = anglegeo(1); color=[1 0 0 0.3; 1 0 0 0.15];
    count = F(crssecx,crssecy);
    dprojx = customprojection([crssecx crssecy],angle);
    dprojx = dprojx(~isnan(count));
    count = count(~isnan(count));
    lnortdlocboot{i} = dprojx;
    lnortcountboot{i} = count;
    fitobj = fit(dprojx,count,'gauss1','StartPoint',[10 0 1]);
    Y1=feval(fitobj,dprojx);
    coef = coeffvalues(fitobj);
    muortboot(i) = coef(2);
    sigmaortboot(i) = coef(3);
    [pkhgt,locs] = findpeaks(count,dprojx);
    pks=[reshape(locs,[],1) reshape(pkhgt,[],1)];
    %only focus on peaks near the origin
    pks = pks(abs(pks(:,1))<=5, :);
    pksortboot{i} = pks;
    %   plot(ax,dprojx,count / normalizer,'-','Color',color(1,:),'linew',1);
    %   plot(ax,dprojx,Y1 / normalizer,'-','Color',color(2,:),'linew',1);
    %   xlim(ax,xran);
    %   hold(ax,'off');
    
  end
  savefile = strcat('bootstrplfedloc',num2str(subwsectar),'swin.mat');
  save(strcat(rstpath, '/MAPS/',savefile), 'lnoptdlocboot','lnoptcountboot','lnortdlocboot',...
    'lnortcountboot','muoptboot','sigmaoptboot','muortboot','sigmaortboot','pks2dboot',...
    'pksoptboot','pksortboot','angleboot','anglegeoboot');

else
  savefile = strcat('bootstrplfedloc',num2str(subwsectar),'swin.mat');
  load(strcat(rstpath, '/MAPS/',savefile));
end
   
keyboard

%% analyzing 2D peaks from random sampling tests
% %%%if using random sampling from a 2D Gaussian
% savefile = strcat('randsamplfedloc',num2str(subwsectar),'swin.mat');
% load(strcat(rstpath, '/MAPS/',savefile));
% angleuse = anglernd;
% pks2duse = pks2drnd;

%%%if using boootstrapping 
savefile = strcat('bootstrplfedloc',num2str(subwsectar),'swin.mat');
load(strcat(rstpath, '/MAPS/',savefile));
angleuse = angleboot;
pks2duse = pks2dboot;

nplt = 2000;
nrow = round(ntrial/nplt);
% f = initfig(8,8,nrow,1);
% figxran = [0.06 0.96]; figyran = [0.06 0.96];
% figxsep = 0.05; figysep = 0.03;
% optaxpos(f,nrow,1,figxran,figyran,figxsep,figysep);
ncount = 0;
for isub = 1: nrow
%   ax=f.ax(isub); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%   ylim(ax,[-5 5]);
  for i = 2000*(isub-1)+1: 2000*isub
%     %%%%%% analyze the peaks on the projected curve
%     dprojx = lnortdloc{i};
%     count = lnortcount{i};
%     [pkhgt,locs,w,p] = findpeaks(count,dprojx);
%     pks=[reshape(locs,[],1) reshape(pkhgt,[],1)];
%     %only focus on peaks near the origin
%     pks=pks(abs(locs)<=5, :);
%     %to be similar to the case like data, there are 3 peaks, one close to 0,
%     %other 2 roughly symmetric to origin
%     npks = size(pks,1);
%     iflag = 0;
%     %if 2 peaks, and symmetric
%     if npks < 2
%       continue
%     elseif npks == 2  
%       loc = pks(:,1);
% %       loc = abs(pks(:,1));  %absolute distance to 0
%       diffloc = diff(abs(loc));  %difference in distance
%       if loc(1)*loc(2)<0 && abs(loc(1))>=2 && abs(loc(2))>=2 && ...
%           abs(diffloc)<=0.5
%         iflag = 1;
%       end 
%     %if 3 peaks, 1 at the center and 2 symmetric   
%     elseif npks == 3
%       loc0 = pks(2,1); %peak that should be closest to 0
%       loc = pks([1 3],1);  %absolute distance to 0
%       diffloc = diff(abs(loc));  %difference in distance
%       if loc(1)*loc(2)<0 && abs(loc(1))>=2 && abs(loc(2))>=2 && ...
%           abs(diffloc)<=0.5 && abs(loc0)<=0.5
%         iflag = 1;
%       end 
%     else
%       loc = pks(:,1);
%       [abslocst,ind] = sort(abs(loc));
%       if loc(ind(2))*loc(ind(3))<0 && abslocst(2)>=2 && abslocst(3)>=2 && ...
%           abs(abslocst(2)-abslocst(3))<=0.5 && abslocst(1)<=0.5
%         iflag = 1;
%       end
%     end
%     %%%%%% analyze the peaks on the projected curve

    %%%%%% analyze the 2d peaks on the smoothed hit count
    angle = angleuse{i};
    pks2d = pks2duse{i};
    npks2d = size(pks2d,1); %must be 2
    iflag = 0;
    if npks2d < 2 || npks2d > 3
      continue
    elseif npks2d == 2  
      signx = pks2d(1,1)*pks2d(2,1);  %must be negative
      signy = pks2d(1,2)*pks2d(2,2);  %must be negative
      dx=diff(pks2d(:,1));
      dy=diff(pks2d(:,2));
      theta=atan2d(dy,dx);
      if theta < 0
        theta=180+theta;   %must be very close to angle(1)
      end
      distpk = sqrt(dx^2+dy^2); %must be >=4
      if signx<0 && signy<0 && abs(theta-angle(1))<=5 && distpk>=4
        iflag = 1;
      end
    elseif npks2d == 3
      distpk0 = sqrt(pks2d(:,1).^2+pks2d(:,2).^2);
      [distpk0st,ind] = sort(distpk0);
      if distpk0st > sqrt(2)
        continue
      else
        pks2d = pks2d(ind(2:end), :);
        signx = pks2d(1,1)*pks2d(2,1);  %must be negative
        signy = pks2d(1,2)*pks2d(2,2);  %must be negative
        dx=diff(pks2d(:,1));
        dy=diff(pks2d(:,2));
        theta=atan2d(dy,dx);
        if theta < 0
          theta=180+theta;   %must be very close to angle(1)
        end
        distpk = sqrt(dx^2+dy^2); %must be >=4
        if signx<0 && signy<0 && abs(theta-angle(1))<=5 && distpk>=4
          iflag = 1;
        end
      end
    end
    ncount = ncount+iflag;
%     if iflag
%       plot(ax,ncount*ones(size(pks,1),1),pks(:,1),'ko-','markersize',2);
%     end
  end
%   longticks(ax,6);
end
% title(f.ax(1),sprintf('%d / %d = %.1f%%',ncount,ntrial,ncount/ntrial*100),...
%   'Units','normalized','HorizontalAlignment','right','FontSize',9);

%% if choose a strip with a finite width in the PCA direction
[~,~,dprojloc] = customprojection(dloc(:,1:2),anglegeo(2));

wid=0.45; % half width of the strip
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

% %%%%%%%%%%% test how many unique points in sample or map space got sampled
% xedge=xran(1)+binw/2:binw:xran(2)-binw/2;
% xcnt=0.5*(xedge(1:end-1)+xedge(2:end));
% nunispl=zeros(length(xedge)-1,1);
% nuni=zeros(length(xedge)-1,1);
% Nnspl=zeros(length(xedge)-1,1);
% Nn=zeros(length(xedge)-1,1);
% for j=1:length(xedge)-1
%   ind=find(dprojlocplt1(:,1)>=xedge(j)&dprojlocplt1(:,1)<xedge(j+1));
%   if ~isempty(ind)      
%     nunispl(j)=size(unique(dlocsplplt1(ind,:),'row'),1);
%     nuni(j)=size(unique(dlocplt1(ind,:),'row'),1);
%     Nnspl(j)=length(ind)/nunispl(j);
%     Nn(j)=length(ind)/nuni(j);
%   end
% end
% f1=initfig(8,4,1,2);
% ax=f1.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% bar(ax,xcnt,nuni,1,'stacked','b','facea',0.6);
% bar(ax,xcnt,nunispl,1,'stacked','r','facea',0.6);
% legend(ax,'relative loc.','sample space');
% ylabel(ax,'Unique points');
% xlabel(ax,'Location difference (km)');
% ax=f1.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% bar(ax,xcnt,Nn,1,'stacked','b','facea',0.6);
% bar(ax,xcnt,Nnspl,1,'stacked','r','facea',0.6);
% ylabel(ax,'Normalized count');
% xlabel(ax,'Location difference (km)');
% %%%%%%%%%%% test how many unique points in sample or map space got sampled

%% summarize the lumped dist comparison between consecutive ones, and each to all others
f = initfig(12,4.5,1,3);
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% p1=histogram(ax,distnn1cat,'binw',0.05,'normalization','count','Facec','b');
p2=histogram(ax,dist,'binw',0.05,'normalization','count');
% plot(ax,[median(distnn1cat) median(distnn1cat)],ax.YLim,'b--','LineWidth',1);
plot(ax,[median(dist) median(dist)],ax.YLim,'k--','LineWidth',1);
xlabel(ax,'Absolute distance (km)');
ylabel(ax,'Count');
legend(ax,[p2],str);
ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% histogram(ax,abs(dprojxy1nn1cat(:,1)),'binw',0.05,'normalization','count','Facec','b');
histogram(ax,abs(dprojxy1(:,1)),'binw',0.05,'normalization','count');
% plot(ax,[median(abs(dprojxy1nn1cat(:,1))) median(abs(dprojxy1nn1cat(:,1)))],ax.YLim,'b--','LineWidth',1);
plot(ax,[median(abs(dprojxy1(:,1))) median(abs(dprojxy1(:,1)))],ax.YLim,'k--','LineWidth',1);
xlabel(ax,'Distance along the 2nd PC (km)');
ylabel(ax,'Count');
ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% histogram(ax,abs(dprojxy2nn1cat(:,1)),'binw',0.05,'normalization','count','Facec','b');
histogram(ax,abs(dprojxy2(:,1)),'binw',0.05,'normalization','count');
% plot(ax,[median(abs(dprojxy2nn1cat(:,1))) median(abs(dprojxy2nn1cat(:,1)))],ax.YLim,'b--','LineWidth',1);
plot(ax,[median(abs(dprojxy2(:,1))) median(abs(dprojxy2(:,1)))],ax.YLim,'k--','LineWidth',1);
xlabel(ax,'Distance along the prop. direc. (km)');
ylabel(ax,'Count');
% keyboard

