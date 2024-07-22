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
% locxyprojalln = allbstnoi.locxyproj4thall;
% tarvlsplstalln = allbstnoi.impindep4thall(:,1);
% nsrcn = allbstnoi.nsrc4th;
% impn = allbstnoi.impindep4thall;
% off1in = allbstnoi.off1i;
% supertstr = 'Further checked at KLNB';
% fnsuffix = '4sta';

impuse = imp;
nsrcuse = nsrc;
fnsuffix2 = [];
% impuse = impn;
% nsrcuse = nsrcn;
% fnsuffix2 = 'noi';

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

% %%
% smoothsigma=5;
% ncont=100;  %num of contour lines
% dx=0.025; dy=0.025; % to distinguish dp near origin, size cannot exceed ~0.053 km 
% cstr={'# events / grid'}; xran=[-4 4]; yran=[-4 4];
% [f,den1d,conmat,conobj]=plt_srcdloc(dloc,'km',3,cstr,...
%   'o','linear','grid',xran,yran,dx,dy,smoothsigma,ncont);

%% data for 25-s wins
subwsectar = 25;  %target subwindow length in sec
% wintype = 'full';
wintype = 'half';

%%%asking each src to all following srcs in 25-s win
if strcmp(wintype,'full')
  if isequaln(impuse,imp)
    savefile = strcat('lfedloc',fnsuffix,num2str(subwsectar),'swin.mat');
  elseif isequaln(impuse,impn)
    savefile = strcat('lfedloc',fnsuffix,num2str(subwsectar),'swin_noi.mat');
  end
%%%asking only each src to others at least 12.5 s later in 25-s win
elseif strcmp(wintype,'half')
  if isequaln(impuse,imp)
    savefile = strcat('lfedloc',fnsuffix,num2str(subwsectar),'shalfwin.mat');
  elseif isequaln(impuse,impn)
    savefile = strcat('lfedloc',fnsuffix,num2str(subwsectar),'shalfwin_noi.mat');
  end
end

load(strcat(rstpath, '/MAPS/',savefile));

dlocspl2all = dloc2all_splcat;
dloc2all = dloc2allcat;
dt2all = dt2allcat;

%%
if isequaln(impuse,imp)
  nrow = 2;
  ncol = 3;
  widin = 8.4;  % maximum width allowed is 8.5 inches
  htin = 7.5;   % maximum height allowed is 11 inches
  f = initfig(widin,htin,nrow,ncol);

  % pltxran = [0.06 0.98]; pltyran = [0.05 0.98]; % optimal axis location
  % pltxsep = 0.01; pltysep = 0.03;
  % optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
  axpos = [0.05 0.54 0.3 0.4;
          0.36 0.54 0.3 0.4;
          0.72 0.570 0.27 0.340;
          0.05 0.06 0.3 0.4;
          0.36 0.06 0.3 0.4;
          0.72 0.090 0.27 0.340];
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
  dplt = dloc;
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
  plot(ax,x0+semib*[coeff(1,2),-coeff(1,2)],y0+semib*[coeff(2,2),-coeff(2,2)],...
    '-','linew',2,'color','b');
  plot(ax,x0+semia*[coeff(1,1),-coeff(1,1)],y0+semia*[coeff(2,1),-coeff(2,1)],...
    ':','linew',2,'color','b');
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
  c.Position = [0.05 0.54 0.3 0.02];
  longticks(ax,2);
  % nolabels(ax,3);
  hold(ax,'off');

  %%%contours of density for N and N-1 w/i 0.375s 
  ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
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
  colormap(ax,'plasma');
  c=colorbar(ax,'SouthOutside');
  %     ax.CLim(2) = prctile(dum(:,3),99);
  ax.CLim(2) = max(conmat(:,1));
  if strcmp(scale,'log10')
    c.Label.String = strcat('log_{10}(',cstr{1},')');
  elseif strcmp(scale,'linear')
    c.Label.String = cstr{1};
  end
  ax.GridLineStyle = '--';
  ax.XAxisLocation = 'top';
  xlim(ax,xran);
  ylim(ax,yran);
  if strcmp(disttype,'spl')
    %     xlabel(ax,'Diff off12 (samples)');
    %     ylabel(ax,'Diff off13 (samples)');
    xticks(ax,xran(1):5:xran(2));
    yticks(ax,yran(1):5:yran(2));
  elseif strcmp(disttype,'km')
    %     xlabel(ax,'Diff E loc (km)');
    %     ylabel(ax,'Diff N loc (km)');
    xticks(ax,xran(1):1:xran(2));
    yticks(ax,yran(1):1:yran(2));
  end
  % pos = ax.Position;
  % c.Position = [pos(1), pos(2)-0.12, pos(3), 0.02];
  c.Position = [0.36 0.54 0.3 0.02];
  longticks(ax,2);
  nolabels(ax,3);
  text(ax,0.02,0.95,'b','FontSize',10,'unit','normalized','EdgeColor','k',...
    'Margin',1,'backgroundcolor','w');
  hold(ax,'off');

  %interpolation to obtain the intersection between PCA directions and contours
  F = scatteredInterpolant(conmat(:,3),conmat(:,4),conmat(:,1),'linear','none');
  normalizer=length(dloc);
  ax=f.ax(3); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
  [ax,muopt,sigmaopt,mdistprojopt,~,countn1,gsfit1]=plt_dloccrssect(ax,F,x,yopt,...
    anglegeo(2),['-';'-'],[.5 .5 .5; 0 0 1],xran,'both',normalizer);
  [ax,muort,sigmaort,mdistprojort,~,countn1r,gsfit1r]=plt_dloccrssect(ax,F,x,yort,...
    anglegeo(1),[':';':'],[.5 .5 .5; 0 0 1],xran,'both',normalizer);
  text(ax,0.01,0.65,sprintf('SE \\mu=%.3f;\nSE \\sigma=%.3f;\nmed(|x|)=%.3f',...
    muopt,sigmaopt,mdistprojopt),'Units','normalized',...
    'HorizontalAlignment','left','FontSize',9,'color','b');
  text(ax,0.62,0.65,sprintf('NE \\mu=%.3f;\nNE \\sigma=%.3f;\nmed(|x|)=%.3f',...
    muort,sigmaort,mdistprojort),'Units','normalized',...
    'HorizontalAlignment','left','FontSize',9,'color','b');
  text(ax,0.02,0.05,'c','FontSize',10,'unit','normalized','EdgeColor','k',...
    'Margin',1,'backgroundcolor','w');
  ylim(ax,[0 1.6e-4]);
  yticks(ax,0: 0.5e-4: 1.5e-4);
  % ylabel(ax,'Probability');

  % supertit(f.ax(1:3),'Between n and n-1 closer than 0.375 s',9);

  %%
  %%%cumulative density for each to all later ones in 25-s windows.
  dplt = dloc2all;
  ax=f.ax(4); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
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
  plot(ax,x0+semib*[coeff(1,2),-coeff(1,2)],y0+semib*[coeff(2,2),-coeff(2,2)],...
    '-','linew',2,'color','r');
  plot(ax,x0+semia*[coeff(1,1),-coeff(1,1)],y0+semia*[coeff(2,1),-coeff(2,1)],...
    ':','linew',2,'color','r');
  % text(ax,0.98,0.16,strcat(num2str(round(anglegeo2all(2))),'$^{\,\circ}$',{'; '},...
  %   num2str(round(anglegeo2all(1))),'$^{\,\circ}$'),'FontSize',11,...
  %   'unit','normalized','interpreter','latex','HorizontalAlignment','right');
  text(ax,0.98,0.17,sprintf('%d%c; %d%c',round(anglegeo2all(2)),char(176),...
    round(anglegeo2all(1)),char(176)),'Units','normalized',...
    'HorizontalAlignment','right','FontSize',9);
  text(ax,0.98,0.11,strcat({'Asp. ratio: '},sprintf('%.2f',semia/semib)),'Units',...
    'normalized','HorizontalAlignment','right','FontSize',9);
  text(ax,0.02,0.95,'d','FontSize',10,'unit','normalized','EdgeColor','k',...
    'Margin',1,'backgroundcolor','w');
  ax.GridLineStyle = '--';
  ax.XAxisLocation = 'top';
  xlim(ax,xran);
  ylim(ax,yran);
  if strcmp(disttype,'spl')
  %   xlabel(ax,'Diff off12 (samples)');
  %   ylabel(ax,'Diff off13 (samples)');
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
  % 0.08 0.55 0.4 0.4
  c.Position = [0.05 0.06 0.3 0.02];
  longticks(ax,2);
  % nolabels(ax,3);
  hold(ax,'off');

  %%%contours of density for each to all later ones in 25-s windows.
  ax=f.ax(5); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
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
  colormap(ax,'plasma');
  c=colorbar(ax,'SouthOutside');
  %     ax.CLim(2) = prctile(dum(:,3),99);
  ax.CLim(2) = max(conmat2all(:,1));
  if strcmp(scale,'log10')
    c.Label.String = strcat('log_{10}(',cstr{1},')');
  elseif strcmp(scale,'linear')
    c.Label.String = cstr{1};
  end
  ax.GridLineStyle = '--';
  ax.XAxisLocation = 'top';
  xlim(ax,xran);
  ylim(ax,yran);
  if strcmp(disttype,'spl')
  %   xlabel(ax,'Diff off12 (samples)');
  %   ylabel(ax,'Diff off13 (samples)');
    xticks(ax,xran(1):5:xran(2));
    yticks(ax,yran(1):5:yran(2));
  elseif strcmp(disttype,'km')
  %   xlabel(ax,'Diff E loc (km)');
  %   ylabel(ax,'Diff N loc (km)');
    xticks(ax,xran(1):1:xran(2));
    yticks(ax,yran(1):1:yran(2));
  end
  % pos = ax.Position;
  % c.Position = [pos(1), pos(2)-0.12, pos(3), 0.02];
  c.Position = [0.36 0.06 0.3 0.02];
  longticks(ax,2);
  nolabels(ax,3);
  text(ax,0.02,0.95,'e','FontSize',10,'unit','normalized','EdgeColor','k',...
    'Margin',1,'backgroundcolor','w');
  hold(ax,'off');

  %interpolation to obtain the intersection between PCA directions and contours
  F2all = scatteredInterpolant(conmat2all(:,3),conmat2all(:,4),conmat2all(:,1),'linear','none');
  normalizer2all=length(dloc2all);
  ax=f.ax(6); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
  [ax,muopt,sigmaopt,mdistprojopt,~,countn2,gsfit2]=plt_dloccrssect(ax,F2all,x,yopt2all,...
    anglegeo2all(2),['-';'-'],[.5 .5 .5; 1 0 0],xran,'both',normalizer2all);
  [ax,muort,sigmaort,mdistprojort,~,countn2r,gsfit2r]=plt_dloccrssect(ax,F2all,x,yort2all,...
    anglegeo2all(1),[':';':'],[.5 .5 .5; 1 0 0],xran,'both',normalizer2all);
  text(ax,0.01,0.65,sprintf('SE \\mu=%.3f;\nSE \\sigma=%.3f;\nmed(|x|)=%.3f',...
    muopt,sigmaopt,mdistprojopt),'Units','normalized',...
    'HorizontalAlignment','left','FontSize',9,'color','r');
  text(ax,0.62,0.65,sprintf('NE \\mu=%.3f;\nNE \\sigma=%.3f;\nmed(|x|)=%.3f',...
    muort,sigmaort,mdistprojort),'Units','normalized',...
    'HorizontalAlignment','left','FontSize',9,'color','r');
  text(ax,0.02,0.05,'f','FontSize',10,'unit','normalized','EdgeColor','k',...
    'Margin',1,'backgroundcolor','w');
  ylim(ax,[0 1.6e-4]);
  yticks(ax,0: 0.5e-4: 1.5e-4);
%   ylabel(ax,'Probability');

  %add the Gaussian fit along the SE direction 
  normalizer1 = max(gsfit2)/max(gsfit1)*normalizer;
  [f.ax(3)]=plt_dloccrssect(f.ax(3),F2all,x,yopt2all,...
    anglegeo2all(2),['-';'-'],[.5 .5 .5; 1 0 0],xran,'gsfit',normalizer1);

  %add the Gaussian fit along the NE direction 
  normalizer2 = max(gsfit1r)/max(gsfit2r)*normalizer2all;
  [f.ax(6)]=plt_dloccrssect(f.ax(6),F,x,yort,anglegeo(1),...
    [':';':'],[.5 .5 .5; 0 0 1],xran,'gsfit',normalizer2);

  lgd1 = legend(f.ax(3),'SE data','SE gsfit','NE data','NE gsfit','SE gsfit in (f)',...
    'Location','north','fontsize',7,'NumColumns',2,'orientation','horizontal');
  %make background transparent
  set(lgd1.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));

  lgd2 = legend(f.ax(6),'SE data','SE gsfit','NE data','NE gsfit','NE gsfit in (c)',...
    'Location','north','fontsize',7,'NumColumns',2,'orientation','horizontal');
  %make background transparent
  set(lgd2.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));

  
elseif isequaln(impuse,impn)
  nrow = 2;
  ncol = 2;
  widin = 5.64;  % maximum width allowed is 8.5 inches
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
         0.07 0.07 0.44 0.4;
         0.60 0.100 0.39 0.34];
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
  dplt = dloc;
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
  [ax,muopt,sigmaopt,mdistprojopt,~,countn1,gsfit1]=plt_dloccrssect(ax,F,x,yopt,...
    anglegeo(2),['-';'-'],[.5 .5 .5; 0 0 1],xran,'both',normalizer);
  [ax,muort,sigmaort,mdistprojort,~,countn1r,gsfit1r]=plt_dloccrssect(ax,F,x,yort,...
    anglegeo(1),[':';':'],[.5 .5 .5; 0 0 1],xran,'both',normalizer);
  text(ax,0.01,0.65,sprintf('SE \\mu=%.3f;\nSE \\sigma=%.3f;\nmed(|x|)=%.3f',...
    muopt,sigmaopt,mdistprojopt),'Units','normalized',...
    'HorizontalAlignment','left','FontSize',9,'color','b');
  text(ax,0.62,0.65,sprintf('NE \\mu=%.3f;\nNE \\sigma=%.3f;\nmed(|x|)=%.3f',...
    muort,sigmaort,mdistprojort),'Units','normalized',...
    'HorizontalAlignment','left','FontSize',9,'color','b');
  text(ax,0.02,0.05,'b','FontSize',10,'unit','normalized','EdgeColor','k',...
    'Margin',1,'backgroundcolor','w');
  ylim(ax,[0 1.6e-4]);
  yticks(ax,0: 0.5e-4: 1.5e-4);
%   ylabel(ax,'Probability');

  %%%cumulative density for each to all later ones in 25-s windows.
  dplt = dloc2all;
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
  %   ylabel(ax,'Diff off13 (samples)');
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
  % 0.08 0.55 0.4 0.4
  c.Position = [0.07 0.07 0.44 0.02];
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
  [ax,muopt,sigmaopt,mdistprojopt,~,countn2,gsfit2]=plt_dloccrssect(ax,F2all,x,...
    yopt2all,anglegeo2all(2),['-';'-'],[.5 .5 .5; 1 0 0],xran,'both',normalizer2all);
  [ax,muort,sigmaort,mdistprojort,~,countn2r,gsfit2r]=plt_dloccrssect(ax,F2all,x,...
    yort2all,anglegeo2all(1),[':';':'],[.5 .5 .5; 1 0 0],xran,'both',normalizer2all);
  text(ax,0.01,0.65,sprintf('SE \\mu=%.3f;\nSE \\sigma=%.3f;\nmed(|x|)=%.3f',...
    muopt,sigmaopt,mdistprojopt),'Units','normalized',...
    'HorizontalAlignment','left','FontSize',9,'color','r');
  text(ax,0.63,0.65,sprintf('NE \\mu=%.3f;\nNE \\sigma=%.3f;\nmed(|x|)=%.3f',...
    muort,sigmaort,mdistprojort),'Units','normalized',...
    'HorizontalAlignment','left','FontSize',9,'color','r');
  text(ax,0.02,0.05,'d','FontSize',10,'unit','normalized','EdgeColor','k',...
    'Margin',1,'backgroundcolor','w');
  ylim(ax,[0 1.6e-4]);
  yticks(ax,0: 0.5e-4: 1.5e-4);
%   ylabel(ax,'Probability');

  %add the Gaussian fit along the SE direction 
  normalizer1 = max(gsfit2)/max(gsfit1)*normalizer;
  [f.ax(2)]=plt_dloccrssect(f.ax(2),F2all,x,yopt2all,...
    anglegeo2all(2),['-';'-'],[.5 .5 .5; 1 0 0],xran,'gsfit',normalizer1);

  %add the Gaussian fit along the NE direction 
  normalizer2 = max(gsfit1r)/max(gsfit2r)*normalizer2all;
  [f.ax(4)]=plt_dloccrssect(f.ax(4),F,x,yort,anglegeo(1),...
    [':';':'],[.5 .5 .5; 0 0 1],xran,'gsfit',normalizer2);

  lgd1 = legend(f.ax(2),'SE data','SE gsfit','NE data','NE gsfit','SE gsfit in (f)',...
    'Location','north','fontsize',7,'NumColumns',2,'orientation','horizontal');
  %make background transparent
  set(lgd1.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));

  lgd2 = legend(f.ax(4),'SE data','SE gsfit','NE data','NE gsfit','NE gsfit in (c)',...
    'Location','north','fontsize',7,'NumColumns',2,'orientation','horizontal');
  %make background transparent
  set(lgd2.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));

end

fname = strcat('lfedloc',fnsuffix,fnsuffix2,'.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));














