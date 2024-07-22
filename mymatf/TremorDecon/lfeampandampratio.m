% lfeampandampratio.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to specifically analyze the amplitude of deconvolved
% impulses. 
% --For example, plotting the map-vew of amp at each station. The idea came 
% from the small gradient in the distribution of average moment over grids
% which might suggest the SE side is lower in amplitude than the NE region.
% Plotting individual amplitude instead the average of 3 stas may be helpful
% --Since amplitude information is already extracted, we can also plotted
% out the amp ratios between stations in terms of histograms. This is also 
% something to be compared with synthetics.  
%
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/05/24
% Last modified date:   2024/05/24
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

stas=['PGC  '
  'SSIB '
  'SILB '
  % 'LZB  '
  % 'TWKB '
  % 'MGCB '
  'KLNB ']; % determine the trio and order, here the 1st sta is PGC
nsta=size(stas,1);         %  number of stations

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
nsrc = allbstsig.nsrc;
imp = allbstsig.impindepall;
off1i = allbstsig.off1i;
off1iwk = allbstsig.off1iwk;
irccrank = allbstsig.irccrank;
ampr = allbstsig.srcamprall;
mamprk = allbstsig.msrcampr;
madamprk = allbstsig.madsrcampr;
nsrcn = allbstnoi.nsrc;
impn = allbstnoi.impindepall;
amprn = allbstnoi.srcamprall;
mamprnk = allbstnoi.msrcampr;
madamprnk = allbstnoi.madsrcampr;
supertstr = 'Secondary sources removed';
fnsuffix = [];
%shift the alignment back to the origin
impsft=shiftsrctoorigin(imp,nsrc,nbst,off1i,off1iwk,irccrank);
[implocsft, ~] = off2space002(impsft(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%%%param for further checked at KLNB
nsrc4th = allbstsig.nsrc4th;
imp4th = allbstsig.impindep4thall;
ampr4th = allbstsig.srcampr4thall;
mamprk4th = allbstsig.msrcampr4th;
madamprk4th = allbstsig.madsrcampr4th;
nsrcn4th = allbstnoi.nsrc4th;
impn4th = allbstnoi.impindep4thall;
amprn4th = allbstnoi.srcampr4thall;
mamprnk4th = allbstnoi.msrcampr4th;
madamprnk4th = allbstnoi.madsrcampr4th;
% supertstr = 'Further checked at KLNB';
% fnsuffix = '4th';
imp4ths=shiftsrctoorigin(imp4th,nsrc4th,nbst,off1i,off1iwk,irccrank);
[imploc4ths, ~] = off2space002(imp4ths(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%% load decon results
tdura = 0.25;  %must be consistent with what synthetics actually used
% tdura = 0.4;

%times of saturation
% nsat=[0.1 0.4 1 2 4 10 20 40 100];
nsat=[0.4 1 2 4 10 20 40 100];
nnsat = length(nsat);

sps = 160;

%%%flag to decide which type of synthetics to use
singleflag = 0;
% singleflag = 1;

if ~singleflag  %%%synthetics from different region sizes and saturation levels
  %%%specify if considering the physical size of each source
  % physicalsize = 1;
  physicalsize = 0;
  
  %%%diameter of physical size
  if physicalsize
    diam=0.15;% 0.5; %0.6; %
  else
    diam=0;
  end
  
  %%%specify shape of the source region
  srcregion='ellipse';
  % srcregion='rectangle';
  % srcregion='circle';
  
  %variation of source region size
  if strcmp(srcregion,'ellipse')
    semia = 1.75*(0.6:0.2:2.0);
    semib = 1.25*(0.6:0.2:2.0);
    nreg = length(semia);
  end
  nround = nreg;
  
  ttstr1 = {'Noise-free syn, '};
  fnsuffix1 = '';
  
else  %%%synthetics from different noise and saturation levels, sources at a single spot
  %different percent of noise
  perctrial = 0.1*(0:2:16)';
  ntrial = length(perctrial);
  nround = ntrial;
  
  ttstr1 = {'Single-spot syn, '};
  fnsuffix1 = '_onespot';
end

for iround = 1: nround
  for insat = 1: nnsat
    sat = nsat(insat);
    
    if ~singleflag
      xaxis = semia(iround);
      yaxis = semib(iround);
      savefile = strcat('rst_decon_synth_reg',num2str(xaxis),...
        '-',num2str(yaxis),'_nsat',num2str(sat),'_td',num2str(tdura),'.mat'); %,srcregion(1:3),'_'
    else
      perc = perctrial(iround);
      savefile = strcat('rst_decon_synth_onespot','_noi',num2str(perc),'_nsat',...
        num2str(sat),'_td',num2str(tdura),'.mat');
    end
    load(strcat(workpath,'/synthetics/',savefile));
      
    impsyn{insat,iround} = allsyn.imp;
    nsrcsyn(insat,iround) = allsyn.nsrcsum;
    srcamprsyn{insat,iround} = allsyn.srcampr;
    imp4thsyn{insat,iround} = allsyn.imp4th;
    nsrc4thsyn(insat,iround) = allsyn.nsrc4thsum;
    srcampr4thsyn{insat,iround} = allsyn.srcampr4th;

  end
end

%%
impuse = imp;
nsrcuse = nsrc;

keyboard
%% amplitude 
%amp at each sta 
for i = 1:3
  lfeamp(:,i) = impuse(:, i*2);
end
%mean amp of 3 stations for the whole catalog
lfeamp(:,4) = mean(impuse(:,[2 4 6]),2);

[denuniq1d, inddup] = density_pixel(impuse(:,7),impuse(:,8)); %count at each unique loc
for i = 1:4
  %%%amp summed at each unique pixel in sample space
  ampsumuniq = sum_at_indices(lfeamp(:,i),inddup);  %sum of all at each unique pixel
  ampsumuniq1d{i} = [denuniq1d(:,1:2) ampsumuniq]; %matrix, moment sum at unique pixels
  
  %%%amp per detection per pixel in sample space
  ampaveuniq1d{i} = [denuniq1d(:,1:2) ampsumuniq./denuniq1d(:,3)];

  %%%amp summed at each grid in map view
  dx = 0.2; dy = 0.2;
  %convert uniq pixel loc to map loc
  tmp = ampsumuniq1d{i};
  locuniq = off2space002(tmp(:,1:2),sps,ftrans,0);
  %bin uniq loc by grid
  [denuniqgrid1d,indgrid] = density_matrix(locuniq(:,1),locuniq(:,2),[-5 5],[-5 5],dx,dy);
  denuniqgrid1d = denuniqgrid1d(denuniqgrid1d(:,3)>0, :);
  indgrid = indgrid(~cellfun('isempty',indgrid));
  ampsumgrid = sum_at_indices(ampsumuniq,indgrid);  %sum of amp at each grid point
  ampsumgrid1d{i} = [denuniqgrid1d(:,1:2) ampsumgrid]; %matrix, amp sum at grid points
  dengrid = sum_at_indices(denuniq1d(:,3),indgrid); %num of events at each grid point
  dengrid1d{i} = [denuniqgrid1d(:,1:2) dengrid]; %matrix, num of events at grid points
  
  %%%amp per detection per grid in map view
  tmp2 = ampsumgrid1d{i};
  ampavegrid1d{i} = [tmp2(:,1:2) ampsumgrid./dengrid];

end

% %% plot of amp averaged over pixel
% nrow = 2;
% ncol = 2;
% widin = 5.3;  % maximum width allowed is 8.5 inches
% htin = 7;   % maximum height allowed is 11 inches
% f = initfig(widin,htin,nrow,ncol);

% pltxran = [0.08 0.98]; pltyran = [0.03 0.98]; % optimal axis location
% pltxsep = 0.02; pltysep = 0.03;
% optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

% % if isequaln(impuse,imps)
% %   xran = [-32 32];
% %   yran = [-32 32];
% % elseif isequaln(impuse,imp) 
%   xran = [-40 40];
%   yran = [-40 40];
% % end
% dx = 1; dy = 1;
% msize = 5;
% scale = 'log10';
% symbol = 'o';
% panels = ['a'; 'b'; 'c'; 'd'];

% % if isequaln(impuse,imps)
% %   cran = [-1.5 0.3];
% % elseif isequaln(impuse,imp) 
%   cran = [-1.5 0.6];
% % end

% for i = 1: 4

%   ax=f.ax(i);
%   hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%   cstr = 'Average moment (Nm) / pixel';
%   sumz1d = ampaveuniq1d{i};
%   sumz1d = sortrows(sumz1d,3);
%   if strcmp(scale,'log10')
%     sumz1d(:,3) = log10(sumz1d(:,3));
%   end
%   scatter(ax,sumz1d(:,1)/sps,sumz1d(:,2)/sps,msize,sumz1d(:,3),symbol,'filled',...
%     'MarkerEdgeColor','none');
%   colormap(ax,'jet');
%   % colormap(ax,'plasma');
%   % colormap(ax,'parula');
%   % colormap(ax,flipud(colormap(ax,'gray')));
%   c=colorbar(ax,'SouthOutside');
%   %   ax.CLim(2) = prctile(sumz1d(:,3),99);
%   caxis(ax,cran);
%   if strcmp(scale,'log10')
%     c.Label.String = strcat('log_{10}(',cstr,')');
%   elseif strcmp(scale,'linear')
%     c.Label.String = cstr;
%   end
%   pos = ax.Position;
%   c.Position = [pos(1), pos(2)-0.09, pos(3), 0.02];
%   text(ax,0.02,0.93,panels(i,:),'FontSize',11,'unit','normalized','EdgeColor','k',...
%     'Margin',1);
%   if i<4
%     text(ax,0.98,0.05,stas(i,:),'FontSize',11,'unit','normalized','HorizontalAlignment',...
%       'right');
%   else
%     text(ax,0.98,0.05,'Mean','FontSize',11,'unit','normalized','HorizontalAlignment',...
%       'right');
%   end
%   axis(ax, 'equal');
%   ax.GridLineStyle = '--';
%   ax.XAxisLocation = 'top';
%   xlim(ax,xran/sps);
%   ylim(ax,yran/sps);
%   longticks(ax,2);
%   if rem(i,2)==0
%     nolabels(ax,3);
%   else
%     xlabel(ax,sprintf('\\Delta{t}_{12} (s)'));
%     ylabel(ax,sprintf('\\Delta{t}_{13} (s)'));
%   end
%   hold(ax,'off');
% end

% if isequaln(impuse,imps)
%   fname = strcat('lfeampsft.pdf');
% elseif isequaln(impuse,imp) 
%   fname = strcat('lfeamp.pdf');
% end
% print(f.fig,'-dpdf',...
%   strcat('/home/chaosong/Pictures/',fname));
% keyboard

%% plot of amp averaged over grids, ie, amp per detecton per grid
nrow = 2;
ncol = 2;
widin = 5.3;  % maximum width allowed is 8.5 inches
htin = 7;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.08 0.98]; pltyran = [0.03 0.98]; % optimal axis location
pltxsep = 0.02; pltysep = 0.03;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

xran = [-4 4];
yran = [-4 4];
dx = 0.2; dy = 0.2;
msize = 30;
scale = 'log10';
symbol = 's';
panels = ['a'; 'b'; 'c'; 'd'];

% if isequaln(impuse,imps)
%   cran = [-1.5 0.3];
% elseif isequaln(impuse,imp) 
  cran = [-1.2 0.2];
% end

%Principal component analysis
locuse = off2space002(impuse(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(locuse(:,1:2));
% plot(ax,ellx,elly,'-','linew',1.5,'color','k');
% plot(ax,x0+semib*[coeff(1,2),-coeff(1,2)],y0+semib*[coeff(2,2),-coeff(2,2)],...
%   '-','linew',2,'color','b');
% plot(ax,x0+semia*[coeff(1,1),-coeff(1,1)],y0+semia*[coeff(2,1),-coeff(2,1)],...
%   ':','linew',2,'color','b');
%a line cross (0,0) in the projection direction
x = reshape(xran(1):0.001:xran(2), [], 1);
yopt = linefcn(x,tand(angle(2)),0);
%a line cross (0,0) in the orthogonal direction
yort = linefcn(x,tand(angle(1)),0);

for i = 1: 4

  ax=f.ax(i);
  hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
  cstr = 'Average amplitude / grid';
  sumz1d = ampavegrid1d{i};
  sumz1d = sortrows(sumz1d,3);
  if strcmp(scale,'log10')
    sumz1d(:,3) = log10(sumz1d(:,3));
  end
  scatter(ax,sumz1d(:,1),sumz1d(:,2),msize,sumz1d(:,3),symbol,'filled',...
    'MarkerEdgeColor','none');
  colormap(ax,'jet');
  % colormap(ax,'plasma');
  % colormap(ax,'parula');
  % colormap(ax,flipud(colormap(ax,'gray')));
  c=colorbar(ax,'SouthOutside');
  %   ax.CLim(2) = prctile(sumz1d(:,3),99);
  caxis(ax,cran);
  if strcmp(scale,'log10')
    c.Label.String = strcat('log_{10}(',cstr,')');
  elseif strcmp(scale,'linear')
    c.Label.String = cstr;
  end
  pos = ax.Position;
  c.Position = [pos(1), pos(2)-0.09, pos(3), 0.02];
  plot(ax,x,yopt,'-','linew',2,'color','b');
  plot(ax,x,yort,':','linew',2,'color','b');
  text(ax,0.02,0.93,panels(i,:),'FontSize',11,'unit','normalized','EdgeColor','k',...
    'Margin',1,'backgroundcolor','w');
  if i<4
    text(ax,0.98,0.05,stas(i,:),'FontSize',11,'unit','normalized','HorizontalAlignment',...
      'right');
  else
    text(ax,0.98,0.05,'Mean','FontSize',11,'unit','normalized','HorizontalAlignment',...
      'right');
  end
  axis(ax, 'equal');
  ax.GridLineStyle = '--';
  ax.XAxisLocation = 'top';
  xlim(ax,xran);
  ylim(ax,yran);
  xticks(ax,xran(1):1:xran(2));
  yticks(ax,yran(1):1:yran(2));
  longticks(ax,2);
  if rem(i,2)==0
    nolabels(ax,3);
  else
    xlabel(ax,'E (km)');
    ylabel(ax,'N (km)');  
  end
  
  hold(ax,'off');

end

if isequaln(impuse,impsft)
  fname = strcat('lfeampsft_map.pdf');
elseif isequaln(impuse,imp) 
  fname = strcat('lfeamp_map.pdf');
end
print(f.fig,'-dpdf',...
  strcat('/home/chaosong/Pictures/',fname));
% keyboard

%% plot of cross-sections
nrow = 2;
ncol = 2;
widin = 5.3;  % maximum width allowed is 8.5 inches
htin = 7;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.08 0.98]; pltyran = [0.03 0.98]; % optimal axis location
pltxsep = 0.02; pltysep = 0.03;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

xran = [-4 4];
dx = 0.2; dy = 0.2;
msize = 30;
scale = 'linear';
if strcmp(scale,'log10') 
  yran = [-1 0.4];
elseif strcmp(scale,'linear') 
  yran = [-0.4 1];
end
symbol = 's';
panels = ['a'; 'b'; 'c'; 'd'];
ncont = 100;

for i = 1: 4

  ax=f.ax(i);
  hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');

  sumz1d = ampavegrid1d{i};
  sumz1d = sortrows(sumz1d,3);
  [~,xgrid,ygrid,zgrid] = ...
    zeropadmat2d(sumz1d,floor(min(sumz1d(:,1))):dx:ceil(max(sumz1d(:,1))),...
    floor(min(sumz1d(:,2))):dy:ceil(max(sumz1d(:,2))));

  zgridgf = zgrid;
  if strcmp(scale,'log10')
    zgridgf = log10(zgridgf);
  end
  [conmat,conobj] = contour(ax,xgrid,ygrid,zgridgf,ncont,'-','linew',1);
  delete(conobj);
  if ~isempty(conmat)
    contable = getContourLineCoordinates(conmat);
    conmat=table2array(contable);
  end
  F = scatteredInterpolant(conmat(:,3),conmat(:,4),conmat(:,1),'linear','none');
  count = F(x,yopt);
  dprojx = customprojection([x yopt],anglegeo(2));
  dprojx = dprojx(~isnan(count));
  count = count(~isnan(count));
  p(1)=plot(ax,dprojx,count,'-','Color','b','linew',2);%[.5 .5 .5]
  % count = F(x,yort);
  % dprojx = customprojection([x yort],anglegeo(1));
  % dprojx = dprojx(~isnan(count));
  % count = count(~isnan(count));
  % p(2)=plot(ax,dprojx,count,':','Color',[.5 .5 .5],'linew',2);
  p(2)=plot(ax,dprojx,flipud(count),'-','Color','r','linew',2);%[.5 .5 .5]
  ampdif = count-flipud(count);
  % ampdif(dprojx>0) = -ampdif(dprojx>0);
  p(3)=plot(ax,dprojx,ampdif,'-','Color','k','linew',2);%[.5 .5 .5]
  xlim(ax,xran);
  ylim(ax,yran);
  if i~=3
    nolabels(ax,3);
  else
    xlabel(ax,'Projected location (km)');
    ylblstr = 'Average amplitude / grid';
    if strcmp(scale,'log10')
      ylabel(ax,strcat('log_{10}(',ylblstr,')'));
    elseif strcmp(scale,'linear')
      ylabel(ax,ylblstr);
    end
    lgd = legend(ax,p,'SE Data','SE Data flipped','diff.','Location','northwest',...
      'fontsize',7,'NumColumns',2);
    %make background transparent
    set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));

  end
  hold(ax,'off');

end

if isequaln(impuse,impsft)
  fname = strcat('lfeampsft_crssec_',scale,'.pdf');
elseif isequaln(impuse,imp) 
  fname = strcat('lfeamp_crssec_',scale,'.pdf');
end
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));
keyboard


%% src amplitude ratio, data vs noise
nrow = 1;
ncol = 3;
widin = 7;  % maximum width allowed is 8.5 inches
htin = 3;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.08 0.98]; pltyran = [0.08 0.98]; % optimal axis location
pltxsep = 0.02; pltysep = 0.03;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

%%%which synthetics to load
insat = 4;
iround = 6;
nsat(insat);
xaxis = semia(iround);
yaxis = semib(iround);
imps=impsyn{insat,iround};
nsrcs=nsrcsyn(insat,iround);
amprs = srcamprsyn{insat,iround};
imp4ths = imp4thsyn{insat,iround};
nsrc4ths = nsrc4thsyn(insat,iround);
ampr4ths = srcampr4thsyn{insat,iround};

%3-station events
[f.ax(1:3),mampr,madampr,stdampr,p(1)] = plt_deconpk_rat_comb4th(f.ax(1:3),ampr,imp,'k','hist');
[f.ax(1:3),mamprs,madamprs,stdamprs,p(2)] = plt_deconpk_rat_comb4th(f.ax(1:3),amprs,imps,'r','hist');
% [f.ax(1:3),mamprn,madamprn] = plt_deconpk_rat_comb4th(f.ax(1:3),amprn,impn,'r','hist');
text(f.ax(1),0.99,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(f.ax(1),0.01,0.7,sprintf('med=%.2f;\nstd=%.2f',mampr(1),stdampr(1)),...
  'Units','normalized','HorizontalAlignment','left','FontSize',9,'color','k');
text(f.ax(1),0.99,0.7,sprintf('med=%.2f;\nstd=%.2f',mamprs(1),stdamprs(1)),...
  'Units','normalized','HorizontalAlignment','right','FontSize',9,'color','r');
text(f.ax(2),0.99,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(f.ax(2),0.01,0.7,sprintf('med=%.2f;\nstd=%.2f',mampr(2),stdampr(2)),...
  'Units','normalized','HorizontalAlignment','left','FontSize',9,'color','k');
text(f.ax(2),0.99,0.7,sprintf('med=%.2f;\nstd=%.2f',mamprs(2),stdamprs(2)),...
  'Units','normalized','HorizontalAlignment','right','FontSize',9,'color','r');
text(f.ax(3),0.99,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(f.ax(3),0.01,0.7,sprintf('med=%.2f;\nstd=%.2f',mampr(3),stdampr(3)),...
  'Units','normalized','HorizontalAlignment','left','FontSize',9,'color','k');
text(f.ax(3),0.99,0.7,sprintf('med=%.2f;\nstd=%.2f',mamprs(3),stdamprs(3)),...
  'Units','normalized','HorizontalAlignment','right','FontSize',9,'color','r');

% ylim(f.ax(:),[0 1.5e3]);
ylim(f.ax(:),[0 0.12]);  
xlim(f.ax(:),[-1 1]);
xticks(f.ax(:),-1: 0.5: 1);

for i = 1:nrow*ncol
  ax=f.ax(i);
  if i == 1
%     xlabel(ax,'');
%     nolabels(ax,1);
  elseif i < 5
    nolabels(ax,3);
    xlabel(ax,'');
    ylabel(ax,'');
  elseif i > 5
    nolabels(ax,2);
    ylabel(ax,'');
  end
end

% legend(f.ax(1),p,'Data, short-win',sprintf('Syn, %.1fx%.1f km, %d',...
%   2*semia(iround),2*semib(iround),nsat(insat)),'location','northwest','fontsize',8);
legend(f.ax(1),p,'Data','Synthetics','location','northwest','fontsize',8);

fname = strcat('lfeamprhist3sta.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));


%% src amplitude ratio, data vs noise
nrow = 2;
ncol = 4;
widin = 8.3;  % maximum width allowed is 8.5 inches
htin = 5.5;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.08 0.98]; pltyran = [0.08 0.98]; % optimal axis location
pltxsep = 0.02; pltysep = 0.03;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

%%%which synthetics to load
insat = 4;
iround = 6;
nsat(insat);
xaxis = semia(iround);
yaxis = semib(iround);
imps=impsyn{insat,iround};
nsrcs=nsrcsyn(insat,iround);
amprs = srcamprsyn{insat,iround};
imp4ths = imp4thsyn{insat,iround};
nsrc4ths = nsrc4thsyn(insat,iround);
ampr4ths = srcampr4thsyn{insat,iround};

%3-station events
[f.ax(1:3),mampr,madampr,stdampr] = plt_deconpk_rat_comb4th(f.ax(1:3),ampr,imp,'k','hist');
[f.ax(1:3),mamprs,madamprs,stdamprs] = plt_deconpk_rat_comb4th(f.ax(1:3),amprs,imps,'r','hist');
% [f.ax(1:3),mamprn,madamprn] = plt_deconpk_rat_comb4th(f.ax(1:3),amprn,impn,'r','hist');
text(f.ax(1),0.99,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(f.ax(1),0.01,0.65,sprintf('med=%.2f;\nstd=%.2f',mampr(1),stdampr(1)),...
  'Units','normalized','HorizontalAlignment','left','FontSize',9,'color','k');
text(f.ax(1),0.99,0.65,sprintf('med=%.2f;\nstd=%.2f',mamprs(1),stdamprs(1)),...
  'Units','normalized','HorizontalAlignment','right','FontSize',9,'color','r');
text(f.ax(2),0.99,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(f.ax(2),0.01,0.65,sprintf('med=%.2f;\nstd=%.2f',mampr(2),stdampr(2)),...
  'Units','normalized','HorizontalAlignment','left','FontSize',9,'color','k');
text(f.ax(2),0.99,0.65,sprintf('med=%.2f;\nstd=%.2f',mamprs(2),stdamprs(2)),...
  'Units','normalized','HorizontalAlignment','right','FontSize',9,'color','r');
text(f.ax(3),0.99,0.95,sprintf('3-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(f.ax(3),0.01,0.65,sprintf('med=%.2f;\nstd=%.2f',mampr(3),stdampr(3)),...
  'Units','normalized','HorizontalAlignment','left','FontSize',9,'color','k');
text(f.ax(3),0.99,0.65,sprintf('med=%.2f;\nstd=%.2f',mamprs(3),stdamprs(3)),...
  'Units','normalized','HorizontalAlignment','right','FontSize',9,'color','r');

%4-station events
p=[];
[f.ax(5:8),mampr4th,madampr4th,stdampr4th,p(1)] = plt_deconpk_rat_comb4th(f.ax(5:8),ampr4th,imp4th,'k','hist');
[f.ax(5:8),mampr4ths,madampr4ths,stdampr4ths,p(2)] = plt_deconpk_rat_comb4th(f.ax(5:8),ampr4ths,imp4ths,'r','hist');
% [f.ax(5:8),mamprn4th,madamprn4th,p(2)] = plt_deconpk_rat_comb4th(f.ax(5:8),amprn4th,impn4th,'r','hist');
text(f.ax(5),0.99,0.95,sprintf('4-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(f.ax(5),0.01,0.65,sprintf('med=%.2f;\nstd=%.2f',mampr4th(1),stdampr4th(1)),...
  'Units','normalized','HorizontalAlignment','left','FontSize',9,'color','k');
text(f.ax(5),0.99,0.65,sprintf('med=%.2f;\nstd=%.2f',mampr4ths(1),stdampr4ths(1)),...
  'Units','normalized','HorizontalAlignment','right','FontSize',9,'color','r');
text(f.ax(6),0.99,0.95,sprintf('4-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(f.ax(6),0.01,0.65,sprintf('med=%.2f;\nstd=%.2f',mampr4th(2),stdampr4th(2)),...
  'Units','normalized','HorizontalAlignment','left','FontSize',9,'color','k');
text(f.ax(6),0.99,0.65,sprintf('med=%.2f;\nstd=%.2f',mampr4ths(2),stdampr4ths(2)),...
  'Units','normalized','HorizontalAlignment','right','FontSize',9,'color','r');
text(f.ax(7),0.99,0.95,sprintf('4-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(f.ax(7),0.01,0.65,sprintf('med=%.2f;\nstd=%.2f',mampr4th(3),stdampr4th(3)),...
  'Units','normalized','HorizontalAlignment','left','FontSize',9,'color','k');
text(f.ax(7),0.99,0.65,sprintf('med=%.2f;\nstd=%.2f',mampr4ths(3),stdampr4ths(3)),...
  'Units','normalized','HorizontalAlignment','right','FontSize',9,'color','r');
text(f.ax(8),0.99,0.95,sprintf('4-station'),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(f.ax(8),0.01,0.65,sprintf('med=%.2f;\nstd=%.2f',mampr4th(4),stdampr4th(4)),...
  'Units','normalized','HorizontalAlignment','left','FontSize',9,'color','k');
text(f.ax(8),0.99,0.65,sprintf('med=%.2f;\nstd=%.2f',mampr4ths(4),stdampr4ths(4)),...
  'Units','normalized','HorizontalAlignment','right','FontSize',9,'color','r');

% ylim(f.ax(:),[0 1.5e3]);
ylim(f.ax(:),[0 0.1]);  
xlim(f.ax(:),[-1 1]);
xticks(f.ax(:),-1: 0.5: 1);

for i = 1:nrow*ncol
  ax=f.ax(i);
  if i == 1
    xlabel(ax,'');
    nolabels(ax,1);
  elseif i < 5
    nolabels(ax,3);
    xlabel(ax,'');
    ylabel(ax,'');
  elseif i > 5
    nolabels(ax,2);
    ylabel(ax,'');
  end
end

legend(f.ax(5),p,'Data, short-win',sprintf('Syn, %.1fx%.1f km, %d',...
  2*semia(iround),2*semib(iround),nsat(insat)),'location','northwest','fontsize',8);

delete(f.ax(4));

fname = strcat('lfeamprhist.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/chaosong/Pictures/',fname));

