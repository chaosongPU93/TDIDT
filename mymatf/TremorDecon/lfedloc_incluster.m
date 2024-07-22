% lfedloc_incluster.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Separated from 'circsrcmodel.m', now this becomes the script that specifically
% compute the distance or location difference between source pair N and N-n
% in the context of a event cluster contains m+1 events. For example, you can
% ask what is distance between consecutive events in all types of clusters,
% ie., doublets, triplets, quaduaplets, to the unpper limit (till the point
% there is a nonzero number of clusters in the catalog).  You can also look 
% at individual clusters, but that would be the main purpose of 'circsrcmodel.m'.
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
[imploc0, ~] = off2space002([0 0],sps,ftrans,0);  % a ref source at 0,0

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
locxyprojalln = allbstnoi.locxyprojall;
tarvlsplstalln = allbstnoi.impindepall(:,1);
nsrcn = allbstnoi.nsrc;
impn = allbstnoi.impindepall;
supertstr = 'Secondary sources removed';
fnsuffix = [];

% %%%param for further checked at KLNB
% locxyprojall = allbstsig.locxyproj4thall;
% tarvlsplstall = allbstsig.impindep4thall(:,1);
% nsrc = allbstsig.nsrc4th;
% imp = allbstsig.impindep4thall;
% locxyprojalln = allbstnoi.locxyproj4thall;
% tarvlsplstalln = allbstnoi.impindep4thall(:,1);
% nsrcn = allbstnoi.nsrc4th;
% impn = allbstnoi.impindep4thall;
% supertstr = 'Further checked at KLNB';
% fnsuffix = '4th';

% impuse = imp;
% nsrcuse = nsrc;
% fnsuffix2 = [];
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

%% during example burst 181, are there clusters?
idxbst = 181;
ist = sum(nsrcuse(1:idxbst-1))+1;
ied = ist+nsrcuse(idxbst)-1;
impi = impuse(ist:ied,:);
timp = impi(:,1)/sps;

catclusibst = [];
nclusibst = zeros(mmax,1);
k = 0;
for m=mmax:-1:1
  catclusbstm = catclusbst{m};
  catclusibstm = catclusbstm{idxbst};
  if isempty(catclusibstm)
    nclusibst(m) = 0;
    continue
  else
    nclusibst(m) = length(catclusibstm);
  end
  
  for j = 1: length(catclusibstm)
    k = k+1;
    impclusibstmj = catclusibstm{j};
    catclusibst{k,1} = catclusibstm{j};
  end
end

for i=1:k
  aaa = catclusibst{i};
  timpi = aaa(:,1)/sps;
end

% keyboard

%% get all eligible, unique event pairs from exclusive clusters
n=1;  %to decide what event pair in the cluster to look at, N and N-n
catimppairst = cell(mmax,1);  %for a certain category, start event of the pair
catimppaired = cell(mmax,1);  %for a certain category, end event of the pair
catdt=cell(mmax,1); %for a certain category, diff time between each event pair
catdloc=cell(mmax,1); %for a certain category, diff loc between each event pair
catdlocspl=cell(mmax,1); %for a certain category, diff loc in spls between each event pair

%%%'m' decides what type of cluster start to look at, eg., if look at all N &
%%%N-1 in doublets, triplets and above, then n=1, and m starts from 1 to mmax
for m=n:mmax
  catclususe = catclus{m};
  nclus = size(catclususe,1); %num of clusters, consecu. clusters may share events!
  imppairst = []; %start event of the pair for ALL clusters, but same category
  imppaired = []; %end event of the pair for ALL clusters, but same category
  for i = 1: nclus  %loop over all clusters
    impi = catclususe{i};
    impipairst = impi(1:end-n, :);  %list of starting src of a PAIR
    impipaired = impi(1+n:end, :);  %list of starting src of a PAIR
    imppairst = [imppairst; impipairst];
    imppaired = [imppaired; impipaired];    
  end
  if size(imppairst,1) ~= nclus*(m+1-n)
    error('number does not match');
  end
  [~,iunist] = unique(imppairst,'rows','stable');
  idupst = setdiff((1:size(imppairst,1))', iunist); %indices for duplicates of start event
  [~,iunied] = unique(imppaired,'rows','stable');
  iduped = setdiff((1:size(imppaired,1))', iunied); %indices for duplicates of end event
  %%%if a index shows in both 'idupst' and 'iduped', then the src PAIR is a duplicate
  idup = intersect(idupst,iduped);
  %%%make sure for certain category of cluster, no duplicate 
  imppairst(idup,:) = []; %remove the duplicates from both lists
  imppaired(idup,:) = [];
  
  %store src pairs for a certain category of cluster
  catimppairst{m} = imppairst;
  catimppaired{m} = imppaired;
  
  %for found UNIQUE event PAIRS
  dt = imppaired(:,1) - imppairst(:,1);
  dloc_spl = imppaired(:,7:8) - imppairst(:,7:8);
  %for diff in loc in map view, do the mapping first, then diff
  implocpairst = off2space002(imppairst(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
  implocpaired = off2space002(imppaired(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
  dloc = implocpaired(:,1:2) - implocpairst(:,1:2);
  catdt{m}=dt;
  catdloc{m}=dloc;
  catdlocspl{m}=dloc_spl;  
end
% keyboard

%%%'lumpst' determine from which category of clusters to lump together, must be >= n
lumpst = n; %e.g., n=1,lumpst=n, means lump all N and N-1 from doublets and above
imppairstlump=cat(1,catimppairst{lumpst:end});
imppairedlump=cat(1,catimppaired{lumpst:end});

%%%in case of duplicate event PAIRS between categories, find them and remove
[~,iunist] = unique(imppairstlump,'rows','stable');
idupst = setdiff((1:size(imppairstlump,1))', iunist);
[~,iunied] = unique(imppairedlump,'rows','stable');
iduped = setdiff((1:size(imppairedlump,1))', iunied);
%if a index shows in both 'idupst' and 'iduped', then the src PAIR is a duplicate
idup = intersect(idupst,iduped);
%make sure for certain category of cluster, no duplicate
imppairstlump(idup,:) = []; %remove the duplicates from both lists
imppairedlump(idup,:) = [];

%obtain the diff location, etc. of the LUMPPED src PAIRS
dt = imppairedlump(:,1) - imppairstlump(:,1);
dloc_spl = imppairedlump(:,7:8) - imppairstlump(:,7:8);
implocpairstlump = off2space002(imppairstlump(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
implocpairedlump = off2space002(imppairedlump(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
dloc = implocpairedlump(:,1:2) - implocpairstlump(:,1:2);

%if there is perfectly no duplicates, for example, you can compute num of event
%pairs for a particular lump starting from 'lumpst' and N-n
npairideal = 0;
for i = lumpst:mmax
  inpair = size(catclus{i}, 1)*(i+1-n);
  npairideal = npairideal+inpair;
end
npair = size(imppairstlump,1);  %but in fact, 'npair' < 'npairideal'

%%
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
[f,den1d,conmat]=plt_srcdloc(dloc,'km',3,cstr,...
  'o','linear','grid',xran,yran,dx,dy,smoothsigma,ncont);
fname = strcat('lfedlocNN1clus',fnsuffix,fnsuffix2,'.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

keyboard


%% plot the diff location in samples, bin by pixel
smoothsigma=1;  %smoothing sigma for Gaussian filtering 
ncont=100;  %num of contour lines
% dx=1; dy=1;
%PCA analysis
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc_spl);
cstr={'# events / pixel'}; xran=[-40 40]; yran=[-40 40];
[f,den1d_spl,conmat,pks]=plt_srcdlocinmap(dt/sps,dloc_spl,[],'spl',timetype,...
  8,cstr,'o','linear','pixel',xran,yran,[],[],smoothsigma,ncont);
ax=f.ax(2);
hold(ax,'on');
plot(ax,ellx,elly,'k-','linew',2);
% plot(ax,x0+10*[coeff(1,2),-coeff(1,2)],y0+10*[coeff(2,2),-coeff(2,2)],'k--','linew',2);
% plot(ax,x0+10*[coeff(1,1),-coeff(1,1)],y0+10*[coeff(2,1),-coeff(2,1)],'k--','linew',2);
text(ax,0.98,0.15,sprintf('%d ^o',round(anglegeo(2))),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.98,0.1,sprintf('asprat: %.1f',semia/semib),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
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
% end

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
% %   contval=contintvl*icont;
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
% x = reshape(xran(1):0.001:xran(2), [], 1);
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

%% summarize the location diff in map space, bin by grid
%%%plot the diff location in map, bin by grid
smoothsigma=5;
ncont=100;  %num of contour lines
dx=0.025; dy=0.025; % to distinguish dp near origin, size cannot exceed ~0.053 km 
%PCA analysis
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc);
cstr={'# events / grid'}; xran=[-4 4]; yran=[-4 4];
[f,den1d,conmat]=plt_srcdlocinmap(dt/sps,dloc,[],'km',timetype,...
  3,cstr,'o','linear','grid',xran,yran,dx,dy,smoothsigma,ncont);
fname = strcat('lfedlocNN1clus',fnsuffix,fnsuffix2,'.pdf');
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
xran=[-4 4]; binw=dx*smoothsigma; legendstr1={'E','N'};
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
  [0 0 0; 0 0 1],xran,normalizer); %,yopt1,yopt2

[f.ax(2),muort,sigmaort,mdistprojort]=plt_dloccrssect(f.ax(2),F,x,yort,anglegeo(1),...
  [.5 .5 .5; 1 0 0],xran,normalizer); %,yort1,yort2
xlabel(f.ax(2),'');
ylabel(f.ax(2),'');

fname = strcat('lfedlocprojNN1clus',fnsuffix,fnsuffix2,'.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

keyboard

% end

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


% keyboard
% %% loc difference between events pairs N and N-n in clusters of m+1 events (N and N-m w/i dtcut)
% n=1;  %to decide which event pair inthe cluster to look at
% catdt=cell(mmax,1);
% catdloc=cell(mmax,1);
% catdlocspl=cell(mmax,1);
% % dprojloclump=cell(mmax,1);
% % dprojloccslump=cell(mmax,1);
% % dprojlocspllump=cell(mmax,1);
% % dprojlocsplcslump=cell(mmax,1);
% for m=n:mmax
%   % m=5;  %to define cluster of events within N and N-m
%   timetype = 'tarvl';
%   % timetype = 'tori';
%   %obtain clusters of events defined by N and N-m w/i dtcut
%   dtcut = 0.25*m+0.125;
%   [impcluster,impuni]=evtcluster(nbst,impuse,nsrcuse,m,dtcut,sps,timetype);
%   implocclus = off2space002(impcluster(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
% 
%   %for each cluster, obtain their diff location between events N and N-n
%   % timetype = 'tori';
%   [dloc,dt,k,dloc_spl]=dloc_evtcluster(impcluster,implocclus,sps,ftrans,m,n,timetype);
%   %note that 'dloc' is NOT the same as the direct inversion of 'dloc_spl'
%   catdt{m}=dt;
%   catdloc{m}=dloc;
%   catdlocspl{m}=dloc_spl;
%   
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
    
%   %there is some elongation along N45E, then we want to project along that direction
%   projang=135;
%   projang=ang2; %note that 'ang2', 2nd pc axis from PCA is very close to the eyeball value 135
%   [~,~,dprojloc] = customprojection(dloc(:,1:2),projang);
%   f=plt_srcdlocinmap(dt/sps,dprojloc,[],'km',timetype,...
%     10,cstr,'o','linear','pixel',xran,yran,dx,dy,0);
%   ndloc(m)=size(dprojloc,1);
%   dprojloclump{m}=dprojloc;
  
%   %histogram of diff location in various ways
%   f = initfig(12,4.5,1,3);
%   xran=[-4 4]; binw=0.1; legendstr1={'E','N'}; normopt='countdensity';
%   f.ax(1)=plt_dlochist(f.ax(1),dloc,xran,binw,legendstr1,normopt);
%   legendstr2={'S45E','N45E'};
%   f.ax(2)=plt_dlochist(f.ax(2),dprojloc,xran,binw,legendstr2,normopt);
%   % [~,mu1(m),sigma1(m),mu2(m),sigma2(m)]=plt_dlochist([],dprojloc,xran,binw,legendstr2,normopt);
%   %focus on a narrow strip along the projection direction
%   dprojloccs=dprojloc(abs(dprojloc(:,2))<=0.2*sqrt(2), :);
%   dprojloccslump{m}=dprojloccs;
%   legendstr3={'S45E-cs'};
%   f.ax(3)=plt_dlochist(f.ax(3),dprojloccs(:,1),xran,binw,legendstr3,normopt);
% 
%  keyboard
% end
% 
% lumpst=n;
% dt=cat(1,catdt{lumpst:end});
% dloc_spl=cat(1,catdlocspl{lumpst:end});
% dloc=cat(1,catdloc{lumpst:end});
% dprojloc=cat(1,dprojloclump{st:end});
% dprojloccs=cat(1,dprojloccslump{st:end});

%% demo of generating random samples from the 2D Gaussian from fitting contours 
%plot the diff location in samples, bin by pixel
smoothsigma=1;
ncont=100;
dx=1; dy=1;
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc_spl);
cstr={'# events / pixel'}; xran=[-40 40]; yran=[-40 40];
[f,den1d_spl,conmat,pks2d]=plt_srcdlocinmap(dt/sps,dloc_spl,[],'spl',timetype,...
  8,cstr,'o','linear','pixel',xran,yran,dx,dy,smoothsigma,ncont);
ax=f.ax(2);
hold(ax,'on');
plot(ax,ellx,elly,'k-','linew',2);
% plot(ax,x0+10*[coeff(1,2),-coeff(1,2)],y0+10*[coeff(2,2),-coeff(2,2)],'k--','linew',2);
% plot(ax,x0+10*[coeff(1,1),-coeff(1,1)],y0+10*[coeff(2,1),-coeff(2,1)],'k--','linew',2);
text(ax,0.98,0.15,sprintf('%d ^o',round(anglegeo(2))),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.98,0.1,sprintf('asprat: %.1f',semia/semib),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
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

distpk0 = sqrt(pks2d(:,1).^2+pks2d(:,2).^2);
[distpk0st,ind] = sort(distpk0);
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

rng('default');
xvec = xran(1):dx:xran(2);
yvec = yran(1):dy:yran(2);
[X1,X2] = meshgrid(xvec,yvec);
X = [X1(:) X2(:)];
mu = [0 0];
sigmax = sigmaort;
sigmay = sigmaopt;
rotang = angle(1);
covar = mvncov(sigmax,sigmay,rotang);
[eigvec,eigval] = eig(covar);
sigmaeig = [sqrt(eigval(1,1)) sqrt(eigval(2,2))]
pdf2 = mvnpdf(X,mu,covar);
figure
imagesc(xvec,yvec,reshape(pdf2,length(yvec),length(xvec))); hold on;% substitute for surf, but need to make 'y' axis direction as normal
ax = gca; ax.YDir = 'normal';
axis equal tight
colormap('jet');
xlabel('x1');
ylabel('x2');
c=colorbar(ax,'SouthOutside');
c.Label.String = 'Probability Density';
R = mvnrnd(mu,covar,size(dloc_spl,1));
scatter(R(:,1),R(:,2),10,[.5 .5 .5],'+');

smoothsigma=1;
ncont=100;
dx=1; dy=1;
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(R);
cstr={'# events / grid'}; xran=[-40 40]; yran=[-40 40];
[f,den1d_spl,conmat,pks2ddemo]=plt_srcdlocinmap(ones(size(dloc_spl,1),1)/sps,R,[],'spl',timetype,...
  8,cstr,'o','linear','grid',xran,yran,dx,dy,smoothsigma,ncont);
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

% %% if choosing a line intersecting with the contour
normalizer=length(dloc);
f = initfig(8,4.5,1,2);
optaxpos(f,1,2);
[f.ax(1),muoptrnddemo,sigmaoptrnddemo,~,pksoptdemo]=plt_dloccrssect(f.ax(1),F,x,yopt,anglegeo(2),...
  [0 0 0; 0 0 1],xran,normalizer); %,yopt1,yopt2

[f.ax(2),muortrnddemo,sigmaortrnddemo,~,pksortdemo]=plt_dloccrssect(f.ax(2),F,x,yort,anglegeo(1),...
  [.5 .5 .5; 1 0 0],xran,normalizer); %,yort1,yort2

%% generating random samples from the 2D Gaussian from fitting contours
xran=[-40 40]; yran=[-40 40];
dx=1; dy=1;
xvec = xran(1):dx:xran(2);
yvec = yran(1):dy:yran(2);
[X1,X2] = meshgrid(xvec,yvec);
X = [X1(:) X2(:)];
mu = [0 0];
sigmax = sigmaort;
sigmay = sigmaopt;
rotang = angle(1);
covar = mvncov(sigmax,sigmay,rotang);
[eigvec,eigval] = eig(covar);
sigmaeig = [sqrt(eigval(1,1)) sqrt(eigval(2,2))]

ntrial = 10000;

% flagrecalc = 0;
flagrecalc = 1;
if flagrecalc

  rng('default');
  
  lnoptdloc=cell(ntrial,1);
  lnoptcount=cell(ntrial,1);
  lnortdloc=cell(ntrial,1);
  lnortcount=cell(ntrial,1);
  pks2drnd=cell(ntrial,1);
  pksoptrnd=cell(ntrial,1);
  pksortrnd=cell(ntrial,1);
  anglernd=cell(ntrial,1);
  anglegeornd=cell(ntrial,1);
  
  % f = initfig(8,4.5,1,2);
  % optaxpos(f,1,2);
  for i = 1: ntrial
    R = mvnrnd(mu,covar,size(dloc_spl,1));
    [coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(R);
    den1d = density_matrix(R(:,1),R(:,2),xran,yran,dx,dy);
    den1d = den1d(den1d(:,3)>0, :);
    den1d = sortrows(den1d,3);
    [~,xgrid,ygrid,zgrid] = ...
      zeropadmat2d(den1d,min(den1d(:,1)):1:max(den1d(:,1)),...
      min(den1d(:,2)):1:max(den1d(:,2)));
    smoothsigma = 1;
    zgridgf = imgaussfilt(zgrid,smoothsigma);  %smooth it a bit
    anglernd{i} = angle;
    anglegeornd{i} = anglegeo;
    
    %find local maxima near the origin of the 2D data 
%     cent=FastPeakFind(zgrid, 50, 1);
    [pkshgt,locs_y,locs_x]=peaks2(zgridgf);
    pks=[reshape(xgrid(1,locs_x),[],1) ...
         reshape(ygrid(locs_y,1),[],1) ...
         reshape(pkshgt,[],1)];
%     ind=find(abs(pksxy(:,1))<=5 & abs(pksxy(:,2))<=5);
    %only focus on peaks near the origin
    pks = pks(abs(pks(:,1))<=5 & abs(pks(:,2))<=5, :);
    pks2drnd{i} = pks;

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
    lnoptdloc{i} = dprojx;
    lnoptcount{i} = count;
    fitobj = fit(dprojx,count,'gauss1','StartPoint',[10 0 1]);
    Y1=feval(fitobj,dprojx);
    coef = coeffvalues(fitobj);
    muoptrnd(i) = coef(2);
    sigmaoptrnd(i) = coef(3);
    [pkhgt,locs] = findpeaks(count,dprojx);
    pks=[reshape(locs,[],1) reshape(pkhgt,[],1)];
    %only focus on peaks near the origin
    pks = pks(abs(pks(:,1))<=5, :);
    pksoptrnd{i} = pks;
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
    lnortdloc{i} = dprojx;
    lnortcount{i} = count;
    fitobj = fit(dprojx,count,'gauss1','StartPoint',[10 0 1]);
    Y1=feval(fitobj,dprojx);
    coef = coeffvalues(fitobj);
    muortrnd(i) = coef(2);
    sigmaortrnd(i) = coef(3);
    [pkhgt,locs] = findpeaks(count,dprojx);
    pks=[reshape(locs,[],1) reshape(pkhgt,[],1)];
    %only focus on peaks near the origin
    pks = pks(abs(pks(:,1))<=5, :);
    pksortrnd{i} = pks;
    %   plot(ax,dprojx,count / normalizer,'-','Color',color(1,:),'linew',1);
    %   plot(ax,dprojx,Y1 / normalizer,'-','Color',color(2,:),'linew',1);
    %   xlim(ax,xran);
    %   hold(ax,'off');
    
  end
  savefile = strcat('randsamplfedloc.mat');
  save(strcat(rstpath, '/MAPS/',savefile), 'lnoptdloc','lnoptcount','lnortdloc',...
    'lnortcount','muoptrnd','sigmaoptrnd','muortrnd','sigmaortrnd','pks2drnd',...
    'pksoptrnd','pksortrnd','anglernd','anglegeornd');

else
  savefile = strcat('randsamplfedloc.mat');
  load(strcat(rstpath, '/MAPS/',savefile));
end


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
  savefile = strcat('bootstrplfedloc.mat');
  save(strcat(rstpath, '/MAPS/',savefile), 'lnoptdlocboot','lnoptcountboot','lnortdlocboot',...
    'lnortcountboot','muoptboot','sigmaoptboot','muortboot','sigmaortboot','pks2dboot',...
    'pksoptboot','pksortboot','angleboot','anglegeoboot');

else
  savefile = strcat('bootstrplfedloc.mat');
  load(strcat(rstpath, '/MAPS/',savefile));
end
   
 
%% analyzing 2D peaks from random sampling tests
% %%%if using random sampling from a 2D Gaussian
% savefile = strcat('randsamplfedloc.mat');
% load(strcat(rstpath, '/MAPS/',savefile));
% angleuse = anglernd;
% pks2duse = pks2drnd;

%%%if using boootstrapping 
savefile = strcat('bootstrplfedloc.mat');
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
  %find peaks in the along-orthogonal profiles
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


%% some test of addition of 2 Guassians
for dmu=5:5:30
x=(1:100)';
mu0=45;
sigmasq=2;
g1=normpdf(x,mu0,sigmasq);
g2=normpdf(x,mu0+dmu,sigmasq);
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
text(0.1,0.9,sprintf('\\Delta\\mu=%d, \\sigma=%d, %.1f',dmu,sigmasq,SIGMAHAT1),'Units','normalized',...
  'HorizontalAlignment','left');
end

%% some test of addition of 2 Guassians
dmu=0.4;
x=(-4:0.01:4)';
mu0=-0.2;
sigmasq=1;
g1=normpdf(x,mu0,sigmasq);
g2=normpdf(x,mu0+dmu,sigmasq);
figure; 
subplot(121)
plot(x,g1,'k'); hold on; plot(x,g2,'k');
g3=g1+g2;
plot(x,g3,'r');
text(0.1,0.9,sprintf('\\Delta\\mu=%.2f, \\sigma=%.2f',dmu,sigmasq),'Units','normalized',...
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


