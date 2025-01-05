% plt_tclustering_byamp.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the script to create a specific plot of the clustering of LFE 
% detections in time, varying with different amplitude, 
% FOR PAPER Song&Rubin2024 part 1.
%
% This plot contains 2 panels. 
% The 2nd panel shows the fraction of arrival time difference
% between events N and N-m that is within the dtcut == 0.25*m+0.125
% s, varying with amplitude, for different m until 5. Here the amp is the 
% median impulse amp of ALL events in between N and N-m. 
% The 1st panel shows the arrival time difference between events N and N-1,
% histogram in log scale, and binned by amplitude. In this case, amp is the 
% average of the median seismic amp of a 6-s window centered around the event
% arrival of the 2 events, different from the 1st panel. 
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/08/10
% Last modified date:   2024/08/10
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

%%%load data
%choose the window length in sec for computing RCC 
% rccwin = 0.25;
rccwin = 0.5;

if rccwin == 0.5
  % savefile = 'deconv_stats4th_allbstsig.mat';
  savefile = 'deconv_stats4th_no23_allbstsig.mat';
elseif rccwin == 0.25
  savefile = 'deconv_stats4th_no23_allbstsig0.25s.mat';
end
load(strcat(rstpath, '/MAPS/',savefile));

if rccwin == 0.5
  % savefile = 'deconv_stats4th_allbstnoi.mat';
  savefile = 'deconv_stats4th_no23_allbstnoi.mat';
elseif rccwin == 0.25
  savefile = 'deconv_stats4th_no23_allbstnoi0.25s.mat';
end
load(strcat(rstpath, '/MAPS/',savefile));


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

impuse = imp;
nsrcuse = nsrc;
% impuse = impn;
% nsrcuse = nsrcn;

typepltnoi = 1; %plot noise
% typepltnoi = 2; %plot data -noise

%% arrival diff. for N and N-1, binned by seismic envelope amp
%%%load the median seismic amp of a 6-s window centered at the detection
%%%2 cols, 1st is the direct amp, 2nd is the envelope
wlensec = 6;

if rccwin==0.25
  savefile = sprintf('medseisampbetweenlfes%s%ds%.2fs.mat',fnsuffix,wlensec,rccwin);
else
  savefile = sprintf('medseisampbetweenlfes%s%ds.mat',fnsuffix,wlensec);
  % savefile = strcat('medseisampbetweenlfes',fnsuffix,num2str(wlensec),'s.mat');
end
load(strcat(rstpath, '/MAPS/',savefile));

%link the seismic amp to specific event pairs N and N-m, use the average of two
m = 1;
dtinter=[];
mamp=[];
for i = 1: nbst
  if nsrcuse(i) == 0
    continue
  end
  ist = sum(nsrcuse(1:i-1))+1;
  ied = ist+nsrcuse(i)-1;
  impi = impuse(ist:ied,:);
  if nsrcuse(i)>=m+1
    %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
    dtfor = diffcustom(impi(:,1), m,'forward'); %to its preceding one
  else
    dtfor = [];
  end
  dtinter = [dtinter; dtfor];
  
  ampi = seisamp(ist:ied,:);  
  aaa=ampi(1:end-m,:);
  bbb=ampi(1+m:end,:);
  mampi=(aaa+bbb)/2;
  
  mamp=[mamp; mampi];
end

% nbin=5;
nbin=10;
binedge=[0 0.375:0.25:ceil(max(dtinter)/sps)];  % 1st bin [0 0.375], then 0.25 increment
bincnt=(binedge(1:end-1)+binedge(2:end))/2;
[ampbin,indbin,n] = binxeqnum(mamp(:,2),nbin);  %bin by amp with same number
medmampbin = zeros(nbin,1);  %median amp of each amp bin
prob1dgeq = cell(nbin,1); %probability of inter-event time > t0, NO time binning 
N1dn = [];
N1dgeq = [];
N1dgeqn = [];
for i = 1: nbin
  indi = indbin{i};
  medmampbin(i) = median(ampbin{i});
  dtinterbin{i} = dtinter(indi);
  N1d=histcounts(dtinterbin{i}/sps,binedge,'normalization','count');
  % N1d=[N1d N1d(end)];
  N1dn(:,i) = reshape(N1d,[],1)./mode(n);
  for j = 1: length(N1d)
    N1dgeq(j,i) = sum(N1d(j:end));
  end
  N1dgeqn(:,i) = N1dgeq(:,i)./length(dtinterbin{i});

  aa=sort(dtinterbin{i}/sps);
  bb = prob_geq(aa);  %probability of inter-event time > t0
  prob1dgeq{i} = bb;
end



%% Count of clusters binned by amp, for each m, similar to Fig. 13c
%fraction of catalog of events in clusters vs. # of events in clusters
%%%determine the 'mmax' for which the resulting number of clusters is nonzero
timetype = 'tarvl';
mmax=getmmaxcluster(nbst,impuse,nsrcuse,sps,timetype);
%%%for a cluster, not only the time separation between N and N-m needs to
%%%be smaller than 'dtcut', but also the max time separation between each
%%%consecutive events needs to smaller than 0.25+0.125 s. 
%%%ie, doublet means a cluster of 2 events ONLY occur as doublets
[catclus,catclusbst,catimp,catuimp,catmedamp,catdtnnm]=...
  evtcluster_ex(nbst,impuse,nsrcuse,mmax,sps,timetype);
%%%fraction of unique events ONLY occurring as certain clusters
[fracuimp,nuimp,fracuimpsum]=frac_uniqevt_incluster2(catuimp,catclus,impuse,nsrcuse,mmax);

%%%get isolated events,  singletons
%all unique events in clusters of >=2 events
catuimplump = cat(1, catuimp{:});

%there should be no duplicates between clusters of different m
[dupinds,induni,inddup]=finddupind(catuimplump);

%deal with singletons, give them amp too
[~, clusind] = ismember(catuimplump,imp,'rows');
isoind = setdiff((1:size(imp,1))', clusind);
isoimp = imp(isoind, :);
ibst = findwhichburst(isoimp,imp,nsrc);
nisoimp = size(isoimp,1);
isoamp = mean(isoimp(:,[2 4 6]),2);
isoimp = [isoimp isoamp (1:nisoimp)' ibst ones(nisoimp,1)];
  
%for other clusters, give each event in the cluster sama amp, as the median amp of the cluster
for m=1:mmax
  catimpm = catimp{m};
  catuimpm = catuimp{m};
  catibstm = catclus{m,2};  %idendifier of which burst this cluster belongs to

  catmedampm = catmedamp{m};
  temp = catimpm;
  ncol = size(temp,2);
  for j = 1: size(catmedampm,1)
    temp((j-1)*(m+1)+1: j*(m+1), ncol+1) = catmedampm(j,1); %median amp of cluster j
    temp((j-1)*(m+1)+1: j*(m+1), ncol+2) = j; %keep track of cluster j for each m
    temp((j-1)*(m+1)+1: j*(m+1), ncol+3) = catibstm(j);  %keep track of which burst
  end
  
  %find duplicates within the clusters of m, expected to exist
  dupinds=finddupind(catimpm);
%   isequal(size(dupinds,1), size(catimpm,1)-size(catuimpm,1))

  %duplicates have diff amp from diff clusters, use the median to reconcile
  for j = 1:size(dupinds,1)
    dupindsj = dupinds{j};
    temp(dupindsj, ncol+1) = median(temp(dupindsj, ncol+1));
  end

  %aafter reconciling, only keep the first unique one
  [~,ind] = unique(temp(:,1:ncol+1),'rows','stable');
  tempuni = temp(ind,:);

  tempuni(:,ncol+4) = m+1;  %keep track of m

  clusimpuni{m,1} = tempuni;

end
%lump in m
clusimp = cat(1,clusimpuni{:});

%combine singletons and other clusters
allimp = [isoimp; clusimp];

%sort by amp first, 'm', then cluster #, finally burst # 
allimp = sortrows(allimp,[ncol+1 ncol+4 ncol+2 ncol+3]);
%%
%%%to avoid breaking the same cluster into diff amp bins, manually choose # in each bin
% nbin2 = 5;
nbin2 = 10;

if ~strcmp(fnsuffix,'4th')
  n2 = [3713; 3714; 3716; 3719; 3714];
  bindiv = [1 3713; 3714 7427; 7428 11143; 11144 14862; 14863 18576];
else
  n2 = [2179; 2179; 2181; 2178; 2177];
  bindiv = [1 2179; 2180 4358; 4359 6539; 6540 8717; 8718 10894];
end

if nbin2 ~=5
  [~,~,n2] = binxeqnum(allimp(:,ncol+1),nbin2);  %bin by amp with same number
  tmp = cumsum(n2);
  bindiv = [[1; tmp(1:end-1)+1] tmp];
end

N2d = [];
N2dn = [];
N2dgeq = [];
N2dgeqn = [];
medmampbin2 = zeros(nbin2,1);  %median amp of each amp bin
for i = 1: nbin2 
  indi = bindiv(i,1): bindiv(i,2);  
  impi = allimp(indi,:);
  medmampbin2(i) = median(impi(:,ncol+1));
  for m=1:mmax+1
    N2d(m,i) = sum(impi(:,ncol+4)==m);
    N2dgeq(m,i) = sum(impi(:,ncol+4)>=m);
  end
  N2dn(:,i) = N2d(:,i) ./ n2(i) *100;
  N2dgeqn(:,i) = N2dgeq(:,i) ./ n2(i) *100;
end

%% choose which one to plot
%%%if a uses cumu prob, and b uses the cumu counts (>= greater equal to)
N1duse = prob1dgeq;
N2duse = N2dgeqn;
yaxflag = 'cumu';

% %%%if use the cumulative binned counts (>= greater equal to)
% N1duse = N1dgeqn;
% N2duse = N2dgeqn;
% yaxflag = 'cumubin';

% %%%if use the binned counts
% N1duse = N1dn;
% N2duse = N2dn;
% yaxflag = 'binned';

%% PLOT
widin = 6;  % maximum width allowed is 8.5 inches
htin = 3.5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f = initfig(widin,htin,nrow,ncol); %initialize fig
pltxran = [0.10 0.98]; pltyran = [0.15 0.96];
pltxsep = 0.10; pltysep = 0.05;
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
% axpos = [0.07 0.1 0.275 0.8;
%          0.36 0.1 0.275 0.8;  
%          0.715 0.1 0.275 0.8];
% for isub = 1:nrow*ncol
%   set(f.ax(isub), 'position', axpos(isub,:));
% end

%%%arrival diff. for N and N-1, binned by seismic envelope amp; binned or cumulative
color = gradientblue(nbin);
xran = [0 7];
% yran = [1e-4 1];
yran = [5e-4 1];
m = 1;
dtcut = 0.25*m+0.125;
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); set(ax,'YScale','log');
patdtcut = [0 yran(2);
          dtcut yran(2);
          dtcut yran(1);
          0 yran(1);
          0 yran(2)];
patch(ax,patdtcut(:,1),patdtcut(:,2),'k','Facealpha',0.15,'edgecolor','none');
plot(ax,[0 30],[1/mode(n) 1/mode(n)],'k--');
for i = 1: nbin
  if ~strcmp(yaxflag, 'cumu')
    % p1d(i)=stairs(ax,binedge,N1duse(:,i),'color',color(i,:),'LineWidth',1);
    p1d(i)=plot(ax,bincnt,N1duse(:,i),'-','Color',color(i,:),'linew',1,...
      'marker','o','markersize',2.5,'markerfacec',color(i,:));
  else
    bb=N1duse{i};
    p1d(i)=plot(ax,bb(:,1),bb(:,2),'-','Color',color(i,:),'linew',1);
  end
  label1d{i} = sprintf('%.2f',medmampbin(i)); %amp of 
end
lgd=legend(ax,p1d(nbin:-1:1),label1d{nbin:-1:1},'FontSize',8,'NumColumns',2,...
  'location','northeast');
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
lgdtit = 'Median envelope amp.';
title(lgd,strcat(lgdtit),'fontsize',8); %,'; short PCA'
% xlabel(ax,'Time (s) from each to its preceding');
% xlabel(ax,'Time delay (s)');
% xlabel(ax,'Inter-event time dt (s)','FontSize',10);
xlabel(ax,'Time t (s)','FontSize',10);
% ylabel(ax,'Normalized count');
if strcmp(yaxflag, 'binned')
  ylabel(ax,'Number in bin (normalized)','FontSize',10);
elseif strcmp(yaxflag, 'cumubin')
  ylabel(ax,sprintf('Prob(inter-event time \x2265 t)'),'FontSize',10);
elseif strcmp(yaxflag, 'cumu') 
  ylabel(ax,sprintf('Prob(inter-event time \x2265 t)'),'FontSize',10);
end
xlim(ax,xran);
ylim(ax,yran);
text(ax,0.02,0.05,'a','FontSize',11,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
longticks(ax,2);
hold(ax,'off');

%%%Count of clusters binned by amp, for each m, similar to Fig. 13c; binned or cumulative
color2 = gradientblue(nbin2);
ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');  
for i = 1: nbin2
  p2d(i)=plot(ax,1:mmax+1,N2duse(:,i),'-','Color',color2(i,:),'linew',1,...
    'marker','o','markersize',2.5,'markerfacec',color2(i,:));    
  label2d{i} = sprintf('%.2f',medmampbin2(i)); %amp of 
end
lgd2=legend(ax,p2d(nbin2:-1:1),label2d{nbin2:-1:1},'FontSize',8,'NumColumns',2,...
  'location','northeast');
set(lgd2.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
lgdtit2 = 'Median impulse amp.';
title(lgd2,strcat(lgdtit2),'fontsize',8);
xlabel(ax,'m (# of events in cluster)','FontSize',10);
if strcmp(yaxflag, 'binned')
  ylabel(ax,sprintf('%% of catalog in clusters of m events'),'FontSize',10);
else
  ylabel(ax,sprintf('%% of catalog in clusters of \x2265 m events'),'FontSize',10);
end
ax.YAxis.Scale = 'log'; %make y axis log10 scale
xlim(ax,[1 mmax+2]);
ylim(ax,[5e-2 1e2]);
xticks(ax,0:2:mmax+2);
yticks(ax,[0.1 0.2 0.5 1 2 5 10 20 50 100]);
text(ax,0.02,0.05,'b','FontSize',11,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
longticks(ax,2);
hold(ax,'off');


% %%%Frac. of arrival diff. w/i (m-1/2)*0.25 s, binned by amp
% mmax = 5;
% % scale = 'linear';
% scale = 'log10';
% ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% [ax,ampbincnt,Nbnall,fracdtb,Nbnnall,fracdtbn]=...
%   plt_fracdt_NNm_mmax(ax,nbst,imp,nsrc,impn,nsrcn,mmax,sps,typepltnoi,scale);
% text(ax,0.02,0.05,'b','FontSize',11,'unit','normalized','EdgeColor','k',...
%   'Margin',1,'backgroundcolor','w');

if ~strcmp(yaxflag, 'cumu')
  if rccwin == 0.25
    fname = sprintf('tclustering_%dampbins%s%.2fs.pdf',nbin,fnsuffix,rccwin);
  else
    fname = sprintf('tclustering_%dampbins%s.pdf',nbin,fnsuffix);
  end  
else
  if rccwin == 0.25
    fname = sprintf('tclustering_%dampbins_cumu%s%.2fs.pdf',nbin,fnsuffix,rccwin);
  else
    fname = sprintf('tclustering_%dampbins_cumu%s.pdf',nbin,fnsuffix);
  end    
end
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

    