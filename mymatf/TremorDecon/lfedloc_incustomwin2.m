% lfedloc_incustomwin2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Very similar to 'lfedloc_incustomwin.m', focusing on all sources
% within a window of length of 25 s,
% Instead of computing the distance between each to all later ones
% within the same window, here we compute the distance between each src and
% all others that are at least 12.5 s later in time.
% --The reason behind this script is because we see some difference between
% distance w/i close-in-time srcs and that between each to all in 25-s windows,
% and we wonder if this difference owes to the propagation of srcs. That's why
% we want to see the distance between srcs that are separated by at least 
% 12.5 s. 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/04/23
% Last modified date:   2024/04/23
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
% windowsk = allbstsig.windowsk;  %stores start and end indices of 25-s windows
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

impuse = imp;
nsrcuse = nsrc;
fnsuffix2 = [];
% impuse = impn;
% nsrcuse = nsrcn;
% fnsuffix2 = 'noi';

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

%%%Flag to indicate if it is necessary to recalculate everything
flagrecalc = 0;
% flagrecalc = 1;

if isequaln(impuse,imp)
  savefile = strcat('lfedloc',fnsuffix,num2str(subwsectar),'shalfwin.mat');
elseif isequaln(impuse,impn)
  savefile = strcat('lfedloc',fnsuffix,num2str(subwsectar),'shalfwin_noi.mat');
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

%     windows = windowsk{i};
%     irccran = irccrank{i};

    %if the burst window length is too short, skip it 
    if tlennew(i)<subwsectar/2
      continue
    end

    %bursts and 4-s detections of the same day
%     indst=1;
%     inded=tlennew(i)*sps;
%     if tlennew(i)<subwsectar
%       subwsec=tlennew(i);   %actually used subwindow length in sec  
%     else
%       subwsec=subwsectar;
%     end
%     subwlen=subwsec*sps;
%     ovlplen=0;
%     windows = movingwins(indst,inded,subwlen,ovlplen,0);

    windows = irccrank{i};
    
    nwin =  size(windows,1);
    iwin = findwhichrange(impi(:,1),windows);
%     keyboard
    for j = 1: nwin
      impiwin = impi(iwin==j,:);
      implociwin = imploci(iwin==j,:);
      lwinsec = round(((windows(j,2)-windows(j,1))+1)/sps); %notes the actual subwin length in sec
      
      n1 = size(impiwin,1);
      for j1 = 1: n1
        ind = find(impiwin(:,1)-impiwin(j1,1) >= subwsectar/2*sps);
        n2 = length(ind);
        if isempty(n2)
          continue
        else
          dt2all = impiwin(ind,1) - impiwin(j1,1);
          dloc2all_spl = impiwin(ind,7:8) - impiwin(j1,7:8);
          dloc2all = implociwin(ind,1:2) - implociwin(j1,1:2);
          dt2allcat = [dt2allcat; dt2all];
          dloc2all_splcat = [dloc2all_splcat; dloc2all_spl];
          dloc2allcat = [dloc2allcat; dloc2all];
        end
      end
      
%       %distinguish the 1st and 2nd halves
%       ind1 = find(impiwin(:,1)>=windows(j,1) & impiwin(:,1)<=round(windows(j,2)/2));
%       ind2 = find(impiwin(:,1)>=round(windows(j,2)/2)+1 & impiwin(:,1)<=windows(j,2));
%       impiwin1 = impiwin(ind1,:);
%       implociwin1 = implociwin(ind1,:);
%       n1 = size(impiwin1,1);
%       impiwin2 = impiwin(ind2,:);
%       implociwin2 = implociwin(ind2,:);
%       n2 = size(impiwin2,1);
% 
%       if n1>=1 && n2>=1
%         % k=k+1;
%         %compute the diff loc between each in the 1st half-win and each in 2nd half-win
%         for j1 = 1: n1
%           for j2 = 1: n2
%             dt2all = impiwin2(j2,1) - impiwin1(j1,1);
%             dloc2all_spl = impiwin2(j2,7:8) - impiwin1(j1,7:8);
%             dloc2all = implociwin2(j2,1:2) - implociwin1(j1,1:2);
%             dt2allcat = [dt2allcat; dt2all];
%             dloc2all_splcat = [dloc2all_splcat; dloc2all_spl];
%             dloc2allcat = [dloc2allcat; dloc2all];
%           end
%         end
%       end
      
    end
  end
  
  save(strcat(rstpath, '/MAPS/',savefile), 'dloc2allcat','dt2allcat','dloc2all_splcat');
  
else
  load(strcat(rstpath, '/MAPS/',savefile));
end

% keyboard

%% type of location difference or distance to look at, between consec. ones, or each to others
%%%%%%%% between each to all others
dloc_spl = dloc2all_splcat;
dloc = dloc2allcat;
dt = dt2allcat;
str=sprintf('between each in 1st half and each in 2nd half');

[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc_spl);
axlena=2*semia
axlenb=2*semib

[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc);
axlena=2*semia
axlenb=2*semib

%% plot the diff location in samples, bin by pixel
smoothsigma=1;  %smoothing sigma for Gaussian filtering 
ncont=100;  %num of contour lines
dx=[]; dy=[]; 
cstr={'# events / pixel'}; xran=[-40 40]; yran=[-40 40];
[f,den1d,conmat]=plt_srcdloc(dloc_spl,'spl',3,cstr,...
  'o','linear','pixel',xran,yran,dx,dy,smoothsigma,ncont);


%% plot the diff location in map km, bin by grid
smoothsigma=5;
ncont=100;  %num of contour lines
dx=0.025; dy=0.025; % to distinguish dp near origin, size cannot exceed ~0.053 km 
cstr={'# events / grid'}; xran=[-4 4]; yran=[-4 4];
[f,den1d,conmat,conobj]=plt_srcdloc(dloc,'km',3,cstr,...
  'o','linear','grid',xran,yran,dx,dy,smoothsigma,ncont);
fname = strcat('lfedloc2all',num2str(subwsectar),'shalfwin',fnsuffix,fnsuffix2,'.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

% keyboard

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
  [0 0 0; 0 0 1],xran,'both',normalizer); %,yopt1,yopt2
[f.ax(2),muort,sigmaort,mdistprojort]=plt_dloccrssect(f.ax(2),F,x,yort,anglegeo(1),...
  [.5 .5 .5; 1 0 0],xran,'both',normalizer); %,yort1,yort2

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
keyboard
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

% fname = strcat('lfedlocproj2all',num2str(subwsectar),'swin',fnsuffix,fnsuffix2,'.pdf');
% print(f.fig,'-dpdf',...
%   strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

keyboard

