% scaled_clusterintime.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The data and noise catalogs have difference sizes but spans on a same length
% in time, therefore the clustering in time for data can be larger partially due
% to that bias. 
% So we need to know how much clustering you would get for a random-in-time
% catalog, like the noise catalog but with the (larger) size of the data 
% catalog.
% --Here's an empirical approach:  Take the number of data detections in each 
% burst, randomly distribute those detections at discrete times separated by
% 0.25 s (with exclusion, meaning only 1 event per 0.25 s), and do the 
% clustering analysis on that synthetic catalog.
% --The next level of sophistication, which I think would be necessary, would
% be to estimate the average excess (percentage) of delays at 0.25 s from 
% your dashed lines in the middle panels for the NOISE, subtract that number 
% of detections from the random scheme outlined above, and then add that number
% back in, forcing them to (randomly) occur at previously-unoccupied times 
% 0.25 s away from the times of existing detections.  This would quantify how  
% much clustering you're seeing in the data, compared to the expectation for
% a random-in-time catalog of the same size (random except for the excess
% at 0.25 s)
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/07/07
% Last modified date:   2024/07/07
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
load(strcat(rstpath, '/MAPS/',savefile));
% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
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


%% compution of fractions for diff catalogs
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

%% Scheme 1
imprand = [];
for i = 1: nbst
  if nsrc(i) == 0
    continue
  end
  ist = sum(nsrc(1:i-1))+1;
  ied = ist+nsrc(i)-1;
  impi = imp(ist:ied,:);
  tevt = impi(:,1);
  tran = trange(i,:);
%   tlen(i)
  
  %divide the time ranges into 0.25-s bins
  binw = 0.25;
  binedges = (tran(2):0.25:tran(3))-tran(2);
  bincnts = (binedges(1:end-1)+binedges(2:end))/2;
  nbins = length(binedges)-1;
  
  if nbins < nsrc(i)
    disp('# of bins are not enough to fit all sources');
  end
  
  %generate random numbers with permutation
%   rng('default');
  seed=i;
  rng(seed);
  indbin = randperm(nbins,nsrc(i));
  indbin = sort(indbin);
  tevtrand = bincnts(indbin);
  impirand = impi;
  impirand(:,1) = reshape(tevtrand*sps,[],1);
  imprand = [imprand; impirand];
end

%%%% FOR RANDOM CATALOGS, 3-station catalog
mmaxr=getmmaxcluster(nbst,imprand,nsrc,sps,timetype);
[catclusr,catclusbstr,catimpr,catuimpr,catmedampr,catdtnnmr]=...
  evtcluster_ex(nbst,imprand,nsrc,mmaxr,sps,timetype);
[fracuimpr,~,fracuimprsum]=frac_uniqevt_incluster2(catuimpr,catclusr,imprand,nsrc,mmaxr);

%%
p=[];
f=initfig(5,5,1,1);
ax = f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
p(1)=plot(ax,1:mmax+1,fracuimpsum,'k-','Linew',1.5,'marker','o','markersize',3,...
  'markerfacec','k');
p(2)=plot(ax,1:mmaxn+1,fracuimpnsum,'r-','Linew',1.5,'marker','o','markersize',3,...
  'markerfacec','r');
p(3)=plot(ax,1:mmaxr+1,fracuimprsum,'b-','Linew',1.5,'marker','o','markersize',3,...
  'markerfacec','b');
legend(ax,p,'Data','Noise','Random','Location','east');
ax.YAxis.Scale = 'log'; %make y axis log10 scale
xlim(ax,[1 mmax+2]);
ylim(ax,[2e-2 1e2]);
xticks(ax,0:2:mmax+2);
yticks(ax,[0.05 0.1 0.2 0.5 1 2 5 10 20 50 100]);
% xlabel(ax,'# of events in cluster','FontSize',10);
% ylabel(ax,'% of catalog in such clusters','FontSize',10);
text(ax,0.5,0.94,supertstr,'HorizontalAlignment','left','Units','normalized',...
  'fontsize',10);
% text(ax,0.98,0.94,'c','FontSize',10,'unit','normalized','EdgeColor','k',...
%   'Margin',1,'backgroundcolor','w','HorizontalAlignment','right');


%% scheme 2, put aside the excess first, then add them back 
xran = [0 7];

m = 1;
nbst = size(trange,1);
dtcut = 0.25*m+0.125;
[amp,dtinter]=med_amp_incluster(nbst,imp,nsrc,m);
[ampn,dtintern]=med_amp_incluster(nbst,impn,nsrcn,m);

binedge=[0 0.375:0.25:ceil(max(dtinter(:,1))/sps)]';  % 1st bin [0 0.375], then 0.25 increment
bincnt=mean([binedge(1:end-1) binedge(2:end)],2);

[Nd]=histcounts(dtinter(:,1)/sps,binedge,'normalization','count');
Ndn = reshape(Nd,[],1)./length(dtinter);
[Nn]=histcounts(dtintern(:,1)/sps,binedge,'normalization','count');
Nnn = reshape(Nn,[],1)./length(dtinter);
fitxran = [2 xran(2)]; 
ind = find(bincnt>=fitxran(1) & bincnt<=fitxran(2));
fitstruct=robustexpfit(bincnt(ind),Ndn(ind),'log10',[1e3 -1]);
fitobj=fitstruct.fitobj;
coef=coeffvalues(fitobj); slpd=coef(2);
yfitd = feval(fitobj,bincnt);

fitstruct=robustexpfit(bincnt(ind),Nnn(ind),'log10',[1e3 -1]);
fitobj=fitstruct.fitobj;
coef=coeffvalues(fitobj); slpn=coef(2);
yfitn = feval(fitobj,bincnt);

%%
usefrac = 0;
if usefrac
  fracsusp = (Nnn(1)-yfitn(1))*length(dtinter)/length(dtintern);
  nexcess = round(fracsusp*length(imp));
  nrest = length(imp)-nexcess;
else
  nexcess = round((Nnn(1)-yfitn(1))*length(dtinter));
  nrest = length(imp)-nexcess;
end

%%
%divide the time ranges into 0.25-s bins
binw = 0.25;
binedges = 0:0.25:tlensum;
bincnts = (binedges(1:end-1)+binedges(2:end))/2;
nbins = length(binedges)-1;

if nbins < nrest
  disp('# of bins are not enough to fit all sources');
end

% seed=0;
% rng(seed);
rng('default');
indbin = randperm(nbins,nrest);
indbin = sort(indbin);
for i = 1:nexcess
  adjbin = findadjacentinteger(nbins,indbin);
  seed=i;
  rng(seed);
  indnew = randi(length(adjbin),1);
  indbin = [indbin adjbin(indnew)];
  indbin = sort(indbin);
end

tevtrand = bincnts(indbin);
imprand2 = imp;
imprand2(:,1) = reshape(tevtrand*sps,[],1);

%%
%%%% FOR RANDOM CATALOGS, 3-station catalog
timetype = 'tarvl';
mmaxr2=getmmaxcluster(1,imprand2,length(imp),sps,timetype);
[catclusr2,catclusbstr2,catimpr2,catuimpr2,catmedampr2,catdtnnmr2]=...
  evtcluster_ex(1,imprand2,length(imp),mmaxr2,sps,timetype);
[fracuimpr2,~,fracuimpr2sum]=frac_uniqevt_incluster2(catuimpr2,catclusr2,imprand2,1,mmaxr2);

%%
p=[];
f=initfig(5,5,1,1);
ax = f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
p(1)=plot(ax,1:mmax+1,fracuimpsum,'k-','Linew',1.5,'marker','o','markersize',3,...
  'markerfacec','k');
p(2)=plot(ax,1:mmaxn+1,fracuimpnsum,'r-','Linew',1.5,'marker','o','markersize',3,...
  'markerfacec','r');
p(3)=plot(ax,1:mmaxr2+1,fracuimpr2sum,'b-','Linew',1.5,'marker','o','markersize',3,...
  'markerfacec','b');
legend(ax,p,'Data','Noise','Random','Location','east');
ax.YAxis.Scale = 'log'; %make y axis log10 scale
xlim(ax,[1 mmax+2]);
ylim(ax,[2e-2 1e2]);
xticks(ax,0:2:mmax+2);
yticks(ax,[0.05 0.1 0.2 0.5 1 2 5 10 20 50 100]);
% xlabel(ax,'# of events in cluster','FontSize',10);
% ylabel(ax,'% of catalog in such clusters','FontSize',10);
text(ax,0.5,0.94,supertstr,'HorizontalAlignment','left','Units','normalized',...
  'fontsize',10);
% text(ax,0.98,0.94,'c','FontSize',10,'unit','normalized','EdgeColor','k',...
%   'Margin',1,'backgroundcolor','w','HorizontalAlignment','right');










