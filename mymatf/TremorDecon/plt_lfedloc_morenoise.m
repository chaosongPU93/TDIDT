% plt_lfedloc.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a plot combining the lfe diff loc from consecutive events separated
% by less than 0.375 s and that from each to all later ones in 25-s windows.
% Basically this script takes in the interested parts from
% 'lfedloc_incluster.m' and 'lfedloc_incustomwin.m', and parts for plotting in
% 'plt_srcdloc.m'.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/04/25
% Last modified date:   2024/04/25
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
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.',...
  num2str(ntol),'.pgc002.',cutout(1:4)));
tlen = trange(:,3)-trange(:,2);
nbst = size(trange,1);
tlensum = sum(tlen);

sps = 160;
ftrans = 'interpchao';

%%%load noise
% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
load(strcat(rstpath, '/MAPS/',savefile));

%%%param for secondary sources removed
locxyprojalln = allbstnoi.locxyprojall;
tarvlsplstalln = allbstnoi.impindepall(:,1);
nsrcn = allbstnoi.nsrc;
impn = allbstnoi.impindepall;
off1in = allbstnoi.off1i;
off1iwkn = allbstnoi.off1iwk;
supertstr = 'Secondary sources removed';
fnsuffix = [];

impuse = impn;
nsrcuse = nsrcn;
fnsuffix2 = 'noi';

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

%% get all eligible, unique N and N-1 event pairs from all types of exclusive clusters
n=1;  %to decide what event pair in the cluster to look at, N and N-n
lumpst = n; %to decide which cluster to start with lumping
lumped = mmax; %to decide which cluster to end with lumping
[imppairst,imppaired,mampcont,amppair,dt,dloc,dlocspl,imppairuevt]=...
  lump_evtpairs_fromcluster(catclus,n,lumpst,lumped,sps,ftrans);

dloclump = [];
dloclump = [dloclump; dloc];

%%% data for 25-s wins
subwsectar = 25;  %target subwindow length in sec
% wintype = 'full';
wintype = 'half';

%%%asking each src to all following srcs in 25-s win
if strcmp(wintype,'full')
  if isequaln(impuse,impn)
    savefile = strcat('lfedloc',fnsuffix,num2str(subwsectar),'swin_noi.mat');
  end
  %%%asking only each src to others at least 12.5 s later in 25-s win
elseif strcmp(wintype,'half')
  if isequaln(impuse,impn)
    savefile = strcat('lfedloc',fnsuffix,num2str(subwsectar),'shalfwin_noi.mat');
  end
end
load(strcat(rstpath, '/MAPS/',savefile));

dloc2all = dloc2allcat;
dloc2alllump = [];
dloc2alllump = [dloc2alllump; dloc2all];


%% load more noise from diff. seed numbers
seedsclpool = 3: 10: 63;
npool=length(seedsclpool);

for iseed = 1: npool
  seedscl=seedsclpool(iseed);
  sclmethod = 'add';
  rccwin = 0.5;
  if rccwin == 0.5
    savefile = sprintf('deconv_stats4th_no23_allbstnoi_sd%s%d.mat',sclmethod,seedscl);
  elseif rccwin == 0.25
    savefile = 'deconv_stats4th_no23_allbstnoi0.25s.mat';
  end
  load(strcat(rstpath, '/MAPS/',savefile));
  
  %%%param for secondary sources removed
  locxyprojalln = allbstnoi.locxyprojall;
  tarvlsplstalln = allbstnoi.impindepall(:,1);
  nsrcn = allbstnoi.nsrc;
  impn = allbstnoi.impindepall;
  off1in = allbstnoi.off1i;
  trangenew = allbstnoi.trangenew;
  % windowsk = allbstsig.windowsk;  %stores start and end indices of 25-s windows
  irccrank = allbstnoi.irccrank;  %stores start and end indices of RCC
  off1iwkn = allbstnoi.off1iwk;
  supertstr = 'Secondary sources removed';
  fnsuffix = [];
  
  impuse = impn;
  nsrcuse = nsrcn;
  fnsuffix2 = 'noi';
  
  %%%%% consecutive events <=0.375 s
  %%%using new ways to generate EXCLUSIVE clusters from each other
  %%%determine the 'mmax' for which the resulting number of clusters is nonzero
  timetype = 'tarvl';
  mmax=getmmaxcluster(nbst,impuse,nsrcuse,sps,timetype);
  %%%for a cluster, not only the time separation between N and N-m needs to
  %%%be smaller than 'dtcut', but also the max time separation between each
  %%%consecutive events needs to smaller than 0.25+0.125 s.
  %%%ie, doublet means a cluster of 2 events ONLY occur as doublets
  [catclus,catclusbst,catimp,catuimp,catmedamp,catdtnnm]=...
    evtcluster_ex(nbst,impuse,nsrcuse,mmax,sps,timetype);
  
  %% get all eligible, unique N and N-1 event pairs from all types of exclusive clusters
  n=1;  %to decide what event pair in the cluster to look at, N and N-n
  lumpst = n; %to decide which cluster to start with lumping
  lumped = mmax; %to decide which cluster to end with lumping
  [imppairst,imppaired,mampcont,amppair,dt,dloc,dlocspl,imppairuevt]=...
    lump_evtpairs_fromcluster(catclus,n,lumpst,lumped,sps,ftrans);
  
  dloclump = [dloclump; dloc];


  %%%%% each event to all later ones >=12.5 s
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
    
    %if the burst window length is too short, skip it
    if tlennew(i)<subwsectar/2
      continue
    end
    
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
      
    end
  end
  
  dloc2all = dloc2allcat;
  dloc2alllump = [dloc2alllump; dloc2all];

end

%%
nrow = 2;
ncol = 2;
widin = 5.6;  % maximum width allowed is 8.5 inches
htin = 7.4;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);

% pltxran = [0.06 0.98]; pltyran = [0.05 0.98]; % optimal axis location
% pltxsep = 0.01; pltysep = 0.03;
% optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
%   axpos = [0.05 0.55 0.45 0.4;
%           0.60 0.582 0.39 0.337;
%           0.05 0.06 0.45 0.4;
%           0.60 0.090 0.39 0.337];
axpos = [0.07 0.55 0.44 0.4;
         0.60 0.580 0.39 0.34;
         0.07 0.09 0.44 0.4;
         0.60 0.120 0.39 0.34];
for isub = 1:nrow*ncol
  set(f.ax(isub), 'position', axpos(isub,:));
end

disttype = 'km';
msize = 3;
cstr={'# events / grid'};
marker = 'o';
scale = 'linear';
binmethod = 'grid';
xran=[-4 4]; yran=[-4 4]; %in km
dx=0.025; dy=0.025; %in km, to distinguish dp near origin, size cannot exceed ~0.053 km
smoothsigma=5;
ncont = 100;

%%%cumulative density for N and N-1 w/i 0.375s
dplt = dloclump;
ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
if strcmp(binmethod,'pixel')
  den1d = density_pixel(dplt(:,1),dplt(:,2));
elseif strcmp(binmethod,'grid')
  den1d = density_matrix(dplt(:,1),dplt(:,2),xran,yran,dx,dy);
end
den1d = den1d(den1d(:,3)>0, :);
den1d = sortrows(den1d,3);

dum = den1d;
dum(dum(:,3)>1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'linew',0.2);  %, 'MarkerEdgeColor', 'k'

dum = sortrows(den1d,3);
dum(dum(:,3)==1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'filled','MarkerEdgeColor','none');
text(ax,0.98,0.05,sprintf('%d pairs',size(dplt,1)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9); %,'FontName','Monospaced'
colormap(ax,'plasma');
c=colorbar(ax,'SouthOutside');
ax.CLim(2) = prctile(dum(:,3),99);
if strcmp(scale,'log10')
  c.Label.String = strcat('log_{10}(',cstr{1},')');
elseif strcmp(scale,'linear')
  c.Label.String = cstr{1};
end
%Principal component analysis
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dplt);
plot(ax,ellx,elly,'-','linew',1.5,'color','k');
% plot(ax,x0+semib*[coeff(1,2),-coeff(1,2)],y0+semib*[coeff(2,2),-coeff(2,2)],...
%   '-','linew',2,'color','b');
% plot(ax,x0+semia*[coeff(1,1),-coeff(1,1)],y0+semia*[coeff(2,1),-coeff(2,1)],...
%   ':','linew',2,'color','b');
% text(ax,0.98,0.16,strcat(num2str(round(anglegeo(2))),'$^{\,\circ}$',{'; '},...
%   num2str(round(anglegeo(1))),'$^{\,\circ}$'),'FontSize',11,...
%   'unit','normalized','interpreter','latex','HorizontalAlignment','right');
text(ax,0.98,0.17,sprintf('%d%c; %d%c',round(anglegeo(2)),char(176),...
  round(anglegeo(1)),char(176)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.98,0.11,strcat({'Asp. ratio: '},sprintf('%.2f',semia/semib)),'Units',...
  'normalized','HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,'a','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
if strcmp(disttype,'spl')
  xlabel(ax,'Diff off12 (samples)');
  ylabel(ax,'Diff off13 (samples)');
  xticks(ax,xran(1):10:xran(2));
  yticks(ax,yran(1):10:yran(2));
elseif strcmp(disttype,'km')
  xlabel(ax,'Location difference E (km)');
  ylabel(ax,'Location difference N (km)');
  xticks(ax,xran(1):1:xran(2));
  yticks(ax,yran(1):1:yran(2));
end
% plot(ax,ax.XLim,ax.YLim,'k--','linew',1);
% pos = ax.Position;
% c.Position = [pos(1), pos(2)-0.12, pos(3), 0.02];
% 0.06 0.54 0.32 0.4;
% 0.4 0.54 0.32 0.4;
c.Position = [0.07 0.55 0.44 0.02];
longticks(ax,2);
% nolabels(ax,3);

if strcmp(disttype,'spl')
  [~,xgrid,ygrid,zgrid] = ...
    zeropadmat2d(den1d,min(den1d(:,1)):1:max(den1d(:,1)),...
    min(den1d(:,2)):1:max(den1d(:,2)));
elseif strcmp(disttype,'km')
  [~,xgrid,ygrid,zgrid] = ...
    zeropadmat2d(den1d,floor(min(den1d(:,1))):dx:ceil(max(den1d(:,1))),...
    floor(min(den1d(:,2))):dy:ceil(max(den1d(:,2))));
end
if ~isempty(smoothsigma)
  zgridgf = imgaussfilt(zgrid,smoothsigma);  %smooth it a bit
else
  zgridgf = zgrid;
end
if strcmp(scale,'log10')
  zgridgf = log10(zgridgf);
end

[conmat,conobj] = contour(ax,xgrid,ygrid,zgridgf,ncont,'-','linew',1);
delete(conobj);
if ~isempty(conmat)
  contable = getContourLineCoordinates(conmat);
  conmat=table2array(contable);
end
%a line cross (0,0) in the projection direction
x = reshape(xran(1):0.001:xran(2), [], 1);
yopt = linefcn(x,tand(angle(2)),0);
plot(ax,x,yopt,'-','linew',2,'color','b');
%a line cross (0,0) in the orthogonal direction
yort = linefcn(x,tand(angle(1)),0);
plot(ax,x,yort,':','linew',2,'color','b');
hold(ax,'off');

%interpolation to obtain the intersection between PCA directions and contours
F = scatteredInterpolant(conmat(:,3),conmat(:,4),conmat(:,1),'linear','none');
normalizer=length(dloc);
ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
[~,muopt,sigmaopt,mdistprojopt,~,countn1,gsfit1]=plt_dloccrssect([],F,x,yopt,...
  anglegeo(2),['-';'-'],[.5 .5 .5; 0 0 1],xran,'both',normalizer);
[~,muort,sigmaort,mdistprojort,~,countn1r,gsfit1r]=plt_dloccrssect([],F,x,yort,...
  anglegeo(1),[':';':'],[.5 .5 .5; 0 0 1],xran,'both',normalizer);
%rather than the median dist of the cross-section, get the median of the whole set
x0=mean(dplt(:,1));
y0=mean(dplt(:,2));
[~,~,projxy] = customprojection(dplt(:,1:2),anglegeo(2));
mdistprojopt = median(abs(projxy(:,1)));
mdistprojort = median(abs(projxy(:,2)));
medprojopt = median(projxy(:,1));
medprojort = median(projxy(:,2));
%   meanprojopt = mean(projxy(:,1));
%   meanprojort = mean(projxy(:,2));
%   text(ax,0.01,0.65,sprintf('SE \\mu=%.3f;\nSE \\sigma=%.2f;\nmed(|x|)=%.2f',...
%     muopt,sigmaopt,mdistprojopt),'Units','normalized',...
%     'HorizontalAlignment','left','FontSize',9,'color','b');
%   text(ax,0.62,0.65,sprintf('NE \\mu=%.3f;\nNE \\sigma=%.2f;\nmed(|x|)=%.2f',...
%     muort,sigmaort,mdistprojort),'Units','normalized',...
%     'HorizontalAlignment','left','FontSize',9,'color','b');
text(ax,0.01,0.65,sprintf('SE: \\sigma=%.2f\nmed(x)=%.3f\nmed(|x|)=%.2f',...
  sigmaopt,medprojopt,mdistprojopt),'Units','normalized',...
  'HorizontalAlignment','left','FontSize',9,'color','b');
text(ax,0.61,0.65,sprintf('NE: \\sigma=%.2f\nmed(x)=%.3f\nmed(|x|)=%.2f',...
  sigmaort,medprojort,mdistprojort),'Units','normalized',...
  'HorizontalAlignment','left','FontSize',9,'color','b');
text(ax,0.02,0.05,'b','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
ylim(ax,[0 7.0e-4]);
yticks(ax,0: 1e-4: 5.0e-4);
%   ylabel(ax,'Probability');

%%%cumulative density for each to all later ones in 25-s windows.
dplt = dloc2alllump;
ax=f.ax(3); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
if strcmp(binmethod,'pixel')
  den1d = density_pixel(dplt(:,1),dplt(:,2));
elseif strcmp(binmethod,'grid')
  den1d = density_matrix(dplt(:,1),dplt(:,2),xran,yran,dx,dy);
end
den1d = den1d(den1d(:,3)>0, :);
den1d = sortrows(den1d,3);

dum = den1d;
dum(dum(:,3)>1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'linew',0.2);  %, 'MarkerEdgeColor', 'k'

dum = sortrows(den1d,3);
dum(dum(:,3)==1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'filled','MarkerEdgeColor','none');
text(ax,0.98,0.05,sprintf('%d pairs',size(dplt,1)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
colormap(ax,'plasma');
c=colorbar(ax,'SouthOutside');
ax.CLim(2) = prctile(dum(:,3),99);
if strcmp(scale,'log10')
  c.Label.String = strcat('log_{10}(',cstr{1},')');
elseif strcmp(scale,'linear')
  c.Label.String = cstr{1};
end
%Principal component analysis
[coeff,score,angle2all,anglegeo2all,x0,y0,semia,semib,ellx,elly]=pcaellipse(dplt);
plot(ax,ellx,elly,'-','linew',1.5,'color','k');
% plot(ax,x0+semib*[coeff(1,2),-coeff(1,2)],y0+semib*[coeff(2,2),-coeff(2,2)],...
%   '-','linew',2,'color','r');
% plot(ax,x0+semia*[coeff(1,1),-coeff(1,1)],y0+semia*[coeff(2,1),-coeff(2,1)],...
%   ':','linew',2,'color','r');
% text(ax,0.98,0.16,strcat(num2str(round(anglegeo2all(2))),'$^{\,\circ}$',{'; '},...
%   num2str(round(anglegeo2all(1))),'$^{\,\circ}$'),'FontSize',11,...
%   'unit','normalized','interpreter','latex','HorizontalAlignment','right');
text(ax,0.98,0.17,sprintf('%d%c; %d%c',round(anglegeo2all(2)),char(176),...
  round(anglegeo2all(1)),char(176)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.98,0.11,strcat({'Asp. ratio: '},sprintf('%.2f',semia/semib)),'Units',...
  'normalized','HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,'c','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
if strcmp(disttype,'spl')
  %   xlabel(ax,'Diff off12 (samples)');
  ylabel(ax,'Diff off13 (samples)');
  xticks(ax,xran(1):10:xran(2));
  yticks(ax,yran(1):10:yran(2));
elseif strcmp(disttype,'km')
  % xlabel(ax,'Location difference E (km)');
  ylabel(ax,'Location difference N (km)');
  xticks(ax,xran(1):1:xran(2));
  yticks(ax,yran(1):1:yran(2));
end
% plot(ax,ax.XLim,ax.YLim,'k--','linew',1);
% pos = ax.Position;
% c.Position = [pos(1), pos(2)-0.12, pos(3), 0.02];
% 0.08 0.55 0.4 0.4
c.Position = [0.07 0.09 0.44 0.02];
longticks(ax,2);
% nolabels(ax,3);

if strcmp(disttype,'spl')
  [~,xgrid,ygrid,zgrid] = ...
    zeropadmat2d(den1d,min(den1d(:,1)):1:max(den1d(:,1)),...
    min(den1d(:,2)):1:max(den1d(:,2)));
elseif strcmp(disttype,'km')
  [~,xgrid,ygrid,zgrid] = ...
    zeropadmat2d(den1d,floor(min(den1d(:,1))):dx:ceil(max(den1d(:,1))),...
    floor(min(den1d(:,2))):dy:ceil(max(den1d(:,2))));
end
if ~isempty(smoothsigma)
  zgridgf = imgaussfilt(zgrid,smoothsigma);  %smooth it a bit
else
  zgridgf = zgrid;
end
if strcmp(scale,'log10')
  zgridgf = log10(zgridgf);
end

[conmat2all,conobj] = contour(ax,xgrid,ygrid,zgridgf,ncont,'-','linew',1);
delete(conobj);
if ~isempty(conmat2all)
  contable = getContourLineCoordinates(conmat2all);
  conmat2all=table2array(contable);
end
%a line cross (0,0) in the projection direction
x = reshape(xran(1):0.001:xran(2), [], 1);
yopt2all = linefcn(x,tand(angle2all(2)),0);
plot(ax,x,yopt2all,'-','linew',2,'color','r');
%a line cross (0,0) in the orthogonal direction
yort2all = linefcn(x,tand(angle2all(1)),0);
plot(ax,x,yort2all,':','linew',2,'color','r');
hold(ax,'off');

%interpolation to obtain the intersection between PCA directions and contours
F2all = scatteredInterpolant(conmat2all(:,3),conmat2all(:,4),conmat2all(:,1),'linear','none');
normalizer2all=length(dloc2all);
ax=f.ax(4); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
[~,muopt,sigmaopt,mdistprojopt,~,countn2,gsfit2]=plt_dloccrssect([],F2all,x,...
  yopt2all,anglegeo2all(2),['-';'-'],[.5 .5 .5; 1 0 0],xran,'both',normalizer2all);
[~,muort,sigmaort,mdistprojort,~,countn2r,gsfit2r]=plt_dloccrssect([],F2all,x,...
  yort2all,anglegeo2all(1),[':';':'],[.5 .5 .5; 1 0 0],xran,'both',normalizer2all);
%rather than the median dist of the cross-section, get the median of the whole set
x0=mean(dplt(:,1));
y0=mean(dplt(:,2));
[~,~,projxy] = customprojection(dplt(:,1:2),anglegeo2all(2));
mdistprojopt = median(abs(projxy(:,1)));
mdistprojort = median(abs(projxy(:,2)));
medprojopt = median(projxy(:,1));
medprojort = median(projxy(:,2));
%   meanprojopt = mean(projxy(:,1))
%   meanprojort = mean(projxy(:,2))
%   text(ax,0.01,0.65,sprintf('SE \\mu=%.3f;\nSE \\sigma=%.2f;\nmed(|x|)=%.2f',...
%     muopt,sigmaopt,mdistprojopt),'Units','normalized',...
%     'HorizontalAlignment','left','FontSize',9,'color','r');
%   text(ax,0.63,0.65,sprintf('NE \\mu=%.3f;\nNE \\sigma=%.2f;\nmed(|x|)=%.2f',...
%     muort,sigmaort,mdistprojort),'Units','normalized',...
%     'HorizontalAlignment','left','FontSize',9,'color','r');
text(ax,0.01,0.65,sprintf('SE: \\sigma=%.2f\nmed(x)=%.3f\nmed(|x|)=%.2f',...
  sigmaopt,medprojopt,mdistprojopt),'Units','normalized',...
  'HorizontalAlignment','left','FontSize',9,'color','r');
text(ax,0.61,0.65,sprintf('NE: \\sigma=%.2f\nmed(x)=%.3f\nmed(|x|)=%.2f',...
  sigmaort,medprojort,mdistprojort),'Units','normalized',...
  'HorizontalAlignment','left','FontSize',9,'color','r');
text(ax,0.02,0.05,'d','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
ylim(ax,[0 7.0e-4]);
yticks(ax,0: 1e-4: 5.0e-4);
%   ylabel(ax,'Probability');

%add the Gaussian fit along the SE direction
ax=f.ax(2); hold(ax,'on');
normalizer1 = max(gsfit2)/max(gsfit1)*normalizer;
[ax,~,~,~,~,~,~,~,pother]=plt_dloccrssect(ax,F2all,x,yopt2all,...
  anglegeo2all(2),['-';'-'],[.5 .5 .5; 1 0 0],xran,'gsfit',normalizer1,1,1);
[ax,~,~,~,~,~,~,~,popt]=plt_dloccrssect(ax,F,x,yopt,...
  anglegeo(2),['-';'-'],[.5 .5 .5; 0 0 1],xran,'both',normalizer,1,2);
[ax,~,~,~,~,~,~,~,port]=plt_dloccrssect(ax,F,x,yort,...
  anglegeo(1),[':';':'],[.5 .5 .5; 0 0 1],xran,'both',normalizer,1,2);
lgd1 = legend(ax,[popt port pother],'SE gsfit','SE profile','NE gsfit','NE profile','SE gsfit in (f)',...
  'Location','north','fontsize',7,'NumColumns',2,'orientation','horizontal');
%make background transparent
set(lgd1.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));

%add the Gaussian fit along the NE direction
ax=f.ax(4); hold(ax,'on');
normalizer2 = max(gsfit1r)/max(gsfit2r)*normalizer2all;
[ax,~,~,~,~,~,~,~,pother]=plt_dloccrssect(ax,F,x,yort,anglegeo(1),...
  [':';':'],[.5 .5 .5; 0 0 1],xran,'gsfit',normalizer2,1,1);
[ax,~,~,~,~,~,~,~,popt]=plt_dloccrssect(ax,F2all,x,...
  yopt2all,anglegeo2all(2),['-';'-'],[.5 .5 .5; 1 0 0],xran,'both',normalizer2all,1,2);
[ax,~,~,~,~,~,~,~,port]=plt_dloccrssect(ax,F2all,x,...
  yort2all,anglegeo2all(1),[':';':'],[.5 .5 .5; 1 0 0],xran,'both',normalizer2all,1,2);
lgd2 = legend(ax,[popt port pother],'SE gsfit','SE profile','NE gsfit','NE profile','NE gsfit in (c)',...
  'Location','north','fontsize',7,'NumColumns',2,'orientation','horizontal');
%make background transparent
set(lgd2.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));

keyboard
fname = strcat('lfedloc',fnsuffix,fnsuffix2,'_more.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));




