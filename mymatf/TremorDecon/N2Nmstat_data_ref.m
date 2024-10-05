% N2Nmstat_data_ref.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the script to carry out a bunch of N to N-m statisical analysis
% to deconvoluted catalog from the real data. The codes are cut from 
% 'deconv_ref_4s_exp_4thsta.m' because they are getting more and more 
% complicated and different from the original scope of 
% 'deconv_ref_4s_exp_4thsta.m'
% --Instead of asking the distance from each src to all others within 2s,
% ie, close-in-time sources, we now move on to ask within 10s (20-s range),
% not only the distance, but also the number of srcs matter. 20-s range 
% is able to constrain the migration to be less than 100m, 20% of the 
% distance, if there is any migration. --- 2024/01/17
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/10/23
% Last modified date:   2023/10/23
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

% dtarvlnn1 = allbstsig.dtarvlnn1all;
% dtarvlnn2 = allbstsig.dtarvlnn2all;
% dtarvlnn3 = allbstsig.dtarvlnn3all;
% dtarvlnn4 = allbstsig.dtarvlnn4all;
% dtarvlnn5 = allbstsig.dtarvlnn5all;
% dtarvlnn6 = allbstsig.dtarvlnn6all;
% dtarvlnn7 = allbstsig.dtarvlnn7all;
% dtarvlnn8 = allbstsig.dtarvlnn8all;
% dtarvlnn1n = allbstnoi.dtarvlnn1all;
% dtarvlnn2n = allbstnoi.dtarvlnn2all;
% dtarvlnn3n = allbstnoi.dtarvlnn3all;
% dtarvlnn4n = allbstnoi.dtarvlnn4all;
% dtarvlnn5n = allbstnoi.dtarvlnn5all;
% dtarvlnn6n = allbstnoi.dtarvlnn6all;
% dtarvlnn7n = allbstnoi.dtarvlnn7all;
% dtarvlnn8n = allbstnoi.dtarvlnn8all;
% dtarvlplt = dtarvlnn1;
% dtarvlpltn = dtarvlnn1n;
% m = 5;

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

% keyboard

%% one plot showing dt_{N, N-1} for data and noise, then summarize fractions for a few m
%%%Median amplitude for srcs w/i the cluster
m = 1;
nbst = size(trange,1);
dtcut = 0.25*m+0.125;
[ampplt,dtarvlplt]=med_amp_incluster(nbst,imp,nsrc,m);  %median amp and dt between N and N-m
[amppltn,dtarvlpltn]=med_amp_incluster(nbst,impn,nsrcn,m);
% keyboard
%%
widin = 8.4;  % maximum width allowed is 8.5 inches
htin = 3.5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 3;
f = initfig(widin,htin,nrow,ncol); %initialize fig

axpos = [0.07 0.1 0.275 0.8;
         0.36 0.1 0.275 0.8;  
         0.715 0.1 0.275 0.8];
for isub = 1:nrow*ncol
  set(f.ax(isub), 'position', axpos(isub,:));
end

xran = [0 2];
nbin = 8;

ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%%%%%%%%%% if bin by amp with a equal number
[ampbin,indbin,n] = binxeqnum(ampplt(:,1),nbin);
% color = flipud(gradientblue(nbin));
color = gradientblue(nbin);
binwdt = 0.05;
nx = round(xran(2)/binwdt)+1;
Nn = zeros(nx, nbin);
edges = xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2;
cnt = xran(1): binwdt: xran(2);
for i = 1: nbin
  ind = indbin{i};
  dtarvli = dtarvlplt(ind,1);
  ampbincntm(i) = median(ampbin{i});
  N=histcounts(dtarvli/sps,edges,'normalization','count');
  Nn(:,i) = N./mode(n);
  fracm(i) = sum(dtarvli/sps<=dtcut)/mode(n);
  p(i)=plot(ax,cnt,Nn(:,i),'color',color(i,:),'LineWidth',1);
  label{i} = sprintf('amp of %.1f',ampbincntm(i)); %
end
%%%%%%%%%% if bin by amp with a equal number
% Nnm = mean(Nn,2);
% p(nbin+1)=plot(ax,cnt,Nnm,'k-','LineWidth',1.5);
% label{nbin+1} = sprintf('mean');  %median
legend(ax,p(nbin:-1:1),label{nbin:-1:1});
ylabel(ax,'Normalized count');
xlabel(ax,sprintf('Time delay (s) between N and N-%d',m));
xlim(ax,xran);
yran = [0 0.2];
ylim(ax,yran);
longticks(ax,2);
text(ax,0.98,0.2,'Data','Units','normalized','HorizontalAlignment','right','FontSize',12);
yticks(ax,xran(1):0.05:xran(2));
hold(ax,'off');

ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%%%%%%%%%% if bin by amp with a equal number
[ampbinn,indbinn,nn] = binxeqnum(amppltn(:,1),nbin);
Nnn = zeros(nx, nbin);
for i = 1: nbin
  ind = indbinn{i};
  dtarvlin = dtarvlpltn(ind,1);
  ampbinncntm(i) = median(ampbinn{i});
  [N,edges]=histcounts(dtarvlin/sps,edges,'normalization','count');
  Nnn(:,i) = N./mode(n);
  if typepltnoi == 1 
    p(i)=plot(ax,cnt,Nnn(:,i),'color',color(i,:),'LineWidth',1);
    fracnm(i) = sum(dtarvlin/sps<=dtcut)/mode(nn);
  elseif typepltnoi == 2 
    p(i)=plot(ax,cnt,Nn(:,i)-Nnn(:,i),'color',color(i,:),'LineWidth',1);
    fracnm(i) = (fracm(i)*mode(n)-sum(dtarvlin/sps<=dtcut)) / (mode(n)-mode(nn));  
  end
  label{i} = sprintf('amp of %.1f',ampbinncntm(i));
end
%%%%%%%%%% if bin by amp with a equal number
% Nnnm = mean(Nnn,2);
% label{nbin+1} = sprintf('mean');  %median
if typepltnoi == 1 
  text(ax,0.98,0.2,'Noise','Units','normalized',...
    'HorizontalAlignment','right','FontSize',12);
  % p(nbin+1)=plot(ax,cnt,Nnnm,'k-','LineWidth',1.5);
  legend(ax,p(nbin:-1:1),label{nbin:-1:1});
elseif typepltnoi == 2 
  title(ax,'Data - Synthetic noise');
  % p(nbin+1)=plot(ax,cnt,Nnm-Nnnm,'k-','LineWidth',1.5);
end
% ylabel(ax,'Normalized count');
% xlabel(ax,sprintf('Arrival time difference N to N-%d (s)',m));
xlim(ax,xran);
ylim(ax,yran);
longticks(ax,2);
yticks(ax,xran(1):0.05:xran(2));
nolabels(ax,3);
hold(ax,'off');

mmax = 5;
% scale = 'linear';
scale = 'log10';
ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
[ax,ampbincnt,Nbnall,fracdtb,Nbnnall,fracdtbn]=...
  plt_fracdt_NNm_mmax(ax,nbst,imp,nsrc,impn,nsrcn,mmax,sps,typepltnoi,scale);

keyboard

fname = strcat('dtnn1binbyamp.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
%%
%%%summarize the whole catalog, diff arrival time and fractions, data and noise together
f1 = initfig(5.5,4,1,1); %initialize fig
tit=supertit(f1.ax,supertstr);
movev(tit,0.2);

f2 = initfig(12,4,1,2); %initialize fig
tit=supertit(f2.ax,supertstr);
movev(tit,0.2);

[f1,f2,Nn,frac,Nnn,fracn,fracdif]=...
  plt_difftime_NNm(f1,f2,dtarvlplt(:,1),dtarvlpltn(:,1),sps,typepltnoi,m);
keyboard

%%
%%%Bin by amp, then plot diff time distribution, and frac w/i some 'dtcut'
f = initfig(10,8,2,2); %initialize fig
% tit=supertit(f.ax,supertstr);
% movev(tit,0.2);
m = 1;
nbin = 8;
[f,Nn,fracm,ampbincntm,Nnn,fracnm,ampbinncntm]=...
  plt_fracdifftime_NNm(f,ampplt(:,1),dtarvlplt(:,1),amppltn(:,1),...
  dtarvlpltn(:,1),sps,typepltnoi,m,nbin);
orient(f.fig,'landscape');
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/agu2023s2f2.pdf'));
keyboard


%% for a few m, bin by median amp for all events within N&N-m, fraction of diff time measurements w/i dtcut  
mmax = 5;
nbst = size(trange,1);
scale = 'linear';
% scale = 'log10';
if strcmp(scale,'linear')
  widin=5; htin=5; nrow=1; ncol=1;
elseif strcmp(scale,'log10')
  widin=5; htin=5; nrow=1; ncol=1;
end
f = initfig(widin,htin,nrow,ncol); %initialize fig
xran = [0.15 0.9]; yran = [0.15 0.9];
xsep = 0.05; ysep = 0.05;
optaxpos(f,nrow,ncol,xran,yran,xsep,ysep);

[f.ax(1),ampbincnt,Nbnall,fracdtb,Nbnnall,fracdtbn]=...
  plt_fracdt_NNm_mmax(f.ax(1),nbst,imp,nsrc,impn,nsrcn,mmax,sps,typepltnoi,scale);

orient(f.fig,'landscape');
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/agu2023s2f3.pdf'));
keyboard




% %% fraction of 'isolated' events, eg, when m=1, evts whose minimum interevt time >0.375s
% %%%2 definitions of inter-event times, one is what we have been used
% %%%the other is the smaller one of the diff time to the left and right
% nbst = size(trange,1);
% m = 1;
% dtcut = 0.25*m+0.125;
% mindtinter=[];
% dtinter=[];
% mindtintern=[];
% dtintern=[];
% for i = 1: nbst
%   if nsrc(i) == 0
%     continue
%   end
%   ist = sum(nsrc(1:i-1))+1;
%   ied = ist+nsrc(i)-1;
%   impi = imp(ist:ied,:);
%   ist = sum(nsrcn(1:i-1))+1;
%   ied = ist+nsrcn(i)-1;
%   impin = impn(ist:ied,:);
%   %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
%   dtfor = diffcustom(impi(:,1), m,'forward'); %to its preceding one
%   dtforpad = [zeros(m,1); dtfor];
%   dtback = diffcustom(impi(:,1), m,'backward'); %to its following one
%   dtbackpad = [dtback; zeros(m,1)]; 
%   tmp1 = [dtforpad dtbackpad];  %time to N-m and N+m for each N
%   dtforn = diffcustom(impin(:,1), m,'forward'); %to its preceding one
%   dtforpadn = [zeros(m,1); dtforn];
%   dtbackn = diffcustom(impin(:,1), m,'backward'); %to its following one
%   dtbackpadn = [dtbackn; zeros(m,1)]; 
%   tmp1n = [dtforpadn dtbackpadn];  %time to N-m and N+m for each N     
%   %choose the min time to neighbors to find isolated ones
%   tmp2 = [dtbackpad(1:m); min(tmp1(m+1: end-m, :),[],2); dtforpad(1:m)];
%   mindtinter = [mindtinter; tmp2];
%   dtinter = [dtinter; dtfor];
%   tmp2n = [dtbackpadn(1:m); min(tmp1n(m+1: end-m, :),[],2); dtforpadn(1:m)];
%   mindtintern = [mindtintern; tmp2n];
%   dtintern = [dtintern; dtforn];
% end
% %if the smaller one of the time to N-m and N+m for the source N is big, it is
% %isolated
% nevtiso = sum(mindtinter/sps>dtcut);
% fraciso = nevtiso/length(imp);
% nevtison = sum(mindtintern/sps>dtcut);
% fracison = nevtison/length(impn);

%% using new ways to generate EXCLUSIVE clusters from each other
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
[fracuimp,nuimp]=frac_uniqevt_incluster2(catuimp,catclus,imp,nsrc,mmax);


%%%% FOR NOISE, 3-station catalog
mmaxn=getmmaxcluster(nbst,impn,nsrcn,sps,timetype);
[catclusn,catclusbstn,catimpn,catuimpn,catmedampn,catdtnnmn]=...
  evtcluster_ex(nbst,impn,nsrcn,mmaxn,sps,timetype);
fracuimpn=frac_uniqevt_incluster2(catuimpn,catclusn,impn,nsrcn,mmaxn);

%%
%%%PLOT
nrow = 1; % rows and cols of subplots in each figure
ncol = 1; 
widin = 4; % size of each figure
htin = 3;
pltxran = [0.15 0.95]; pltyran = [0.15 0.95];
pltxsep = 0.02; pltysep = 0.03;
f = initfig(widin,htin,nrow,ncol);
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

ax = f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
p(1)=plot(ax,1:mmax+1,fracuimp,'k-','Linew',1.5,'marker','o','markersize',4,...
  'markerfacec','k');
p(2)=plot(ax,1:mmaxn+1,fracuimpn,'r-','Linew',1.5,'marker','o','markersize',4,...
  'markerfacec','r');
legend(ax,p,'Data','Noise','Location','northeast');
ax.YAxis.Scale = 'log'; %make y axis log10 scale
yticks(ax,[0.1 0.2 0.5 1 2 5 10 20 50 100]);
xlabel(ax,'# of events in cluster','FontSize',10);
ylabel(ax,'% of catalog in such clusters','FontSize',10);

%% # of unique event pairs from ALL clusters? eg., N and N-1 separated by<=0.375s
%%%Then you need to combine unique event pairs from different exclusive clusters
imppairlump = cell(mmax,2); %src info for unique evt pairs combined from diff clusters
amppairlump = cell(mmax,1); %mean amp for such evet pairs
ampcontlump = cell(mmax,1); %median amp of all continous evts between the pair, not just the pair itself  
dtpairlump = cell(mmax,1);  %diff time between such evet pairs
for n = 1:mmax
  % n=1;  %to decide what event pair in the cluster to look at, N and N-n
  catimppairst = cell(mmax,1);  %for a certain category, start event of the pair
  catimppaired = cell(mmax,1);  %for a certain category, end event of the pair
  catmampcont = cell(mmax,1);  %for a certain category, median amp of ALL evts between event pair
  
  %%%'m' decides what type of cluster start to look at, eg., if look at all N &
  %%%N-1 in doublets, triplets and above, then n=1, and m starts from 1 to mmax
  for m=n:mmax
    catclususe = catclus{m};
    nclus = size(catclususe,1); %num of clusters, consecu. clusters may share events!
    imppairst = []; %start event of the pair for ALL clusters, but same category
    imppaired = []; %end event of the pair for ALL clusters, but same category
    mampcont = []; %median amp of all continous evts between the start and end of the pair for ALL clusters, but same category
    for i = 1: nclus  %loop over all clusters
      impi = catclususe{i};
      impipairst = impi(1:end-n, :);  %list of starting src of a PAIR
      impipaired = impi(1+n:end, :);  %list of starting src of a PAIR
      imppairst = [imppairst; impipairst];
      imppaired = [imppaired; impipaired]; 
      %median amp of all continous events between 'imppairst' and 'imppaired' 
      mampicont = zeros(size(impipairst,1), 1);
      for j = 1: size(impipairst,1)
        impcont = impi(j: j+n,:);
        mampicont(j,1) = median(mean(impcont(:,[2 4 6]),2));
      end
      mampcont = [mampcont; mampicont]; 
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
    %%%make sure for certain category of cluster, no duplicates 
    imppairst(idup,:) = []; %remove the duplicates from both lists
    imppaired(idup,:) = [];
    mampcont(idup) = [];

    %store src pairs for a certain category of cluster
    catimppairst{m} = imppairst;
    catimppaired{m} = imppaired;
    catmampcont{m} = mampcont;
  end
  % keyboard

  %%%'lumpst' determine from which category of clusters to lump together, must be >= n
  lumpst = n; %e.g., n=1,lumpst=n, means lump all N and N-1 from doublets and above
  imppairstlump=cat(1,catimppairst{lumpst:end});
  imppairedlump=cat(1,catimppaired{lumpst:end});
  mampcontlump=cat(1,catmampcont{lumpst:end});

  %%%in case of duplicate event PAIRS between categories, find them and remove
  [~,iunist] = unique(imppairstlump,'rows','stable');
  idupst = setdiff((1:size(imppairstlump,1))', iunist);
  [~,iunied] = unique(imppairedlump,'rows','stable');
  iduped = setdiff((1:size(imppairedlump,1))', iunied);
  %if a index shows in both 'idupst' and 'iduped', then the src PAIR is a duplicate
  idup = intersect(idupst,iduped);
  %make sure for certain category of cluster, no duplicates
  imppairstlump(idup,:) = []; %remove the duplicates from both lists
  imppairedlump(idup,:) = [];
  mampcontlump(idup,:) = [];
  
  %store for each n
  imppairlump{n,1} = imppairstlump;
  imppairlump{n,2} = imppairedlump;
  amppairlump{n} = mean([mean(imppairstlump(:,[2 4 6]),2) ...
                      mean(imppairedlump(:,[2 4 6]),2)], 2); %mean amp of the pair
  ampcontlump{n} = mampcontlump;  %median amp of all evts between the pair                  
  dtpairlump{n} = imppairedlump(:,1)-imppairstlump(:,1);              
  nuimppair(n,1) = size(imppairstlump,1); %num of unique evt pairs

  %%%unlike in 'lfedloc_incluster', essentially, we want to know the UNIQUE EVENTS
  %%%associated with these event pairs
  imppairevtlump = [imppairstlump; imppairedlump];
  imppairuevtlump = unique(imppairevtlump,'rows','stable');
  nuimpNNn(n,1) = size(imppairuevtlump,1);  %num of unique evts associated with pairs
end

%%
nrow = 1; % rows and cols of subplots in each figure
ncol = 1; 
widin = 5; % size of each figure
htin = 4;
pltxran = [0.15 0.95]; pltyran = [0.15 0.95];
pltxsep = 0.02; pltysep = 0.03;
f = initfig(widin,htin,nrow,ncol);
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

ax = f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% %%%note there can be duplicates in 'nuimpNNn' for different n, so the sum is NOT
% %%%the whol catalog
% plot(ax,1:mmax,nuimpNNn/length(impuse)*100,'k-','Linew',1.5,'marker','o','markersize',6);
% ylabel(ax,'% of catalog of such events','FontSize',10);
%%%if plotting 'nuimppair', not sure how to normalize it
plot(ax,1:mmax,nuimppair,'k-','Linew',1.5,'marker','o','markersize',6);
ylabel(ax,'# of such event pairs','FontSize',10);
xlabel(ax,'Event pair N and N-m w/i 0.25*m+0.125 s','FontSize',10);

% %% for a few m, fraction of unique events inside clusters w/i a diff time cut, wrt. all catalog for data
% mmax = 15;
% nbst = size(trange,1);
% 
% [fracsrc2all, dfracsrc2all,mmaxzero, f]=frac_uniqevt_incluster(nbst,imp,nsrc,mmax,sps);
% [fracsrc2alln, dfracsrc2alln,mmaxzeron, f]=frac_uniqevt_incluster(nbst,impn,nsrcn,mmax,sps);
% aa = (fracsrc2all*size(imp,1) - fracsrc2alln*size(impn,1))/100;
% fracsrc2alldmn = aa/(size(imp,1)-size(impn,1))*100;
% bb = (dfracsrc2all*size(imp,1) - dfracsrc2alln*size(impn,1))/100;
% dfracsrc2alldmn = bb/(size(imp,1)-size(impn,1))*100;
% % f=initfig;
% % ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% % plot(ax,1:mmaxzero,dfracsrc2alldmn(1:mmaxzero),'ro-','linew',1,'markersize',4);
% % xlabel(ax,'m');
% % ylabel(ax,'Frac of srcs');
% % title(ax,'Frac of srcs in cluster of N & N-m whose diff time w/i 0.25*m+0.125 s');
% % legend(ax,'exclusive, data - synthetic noise');
% % xlim(ax,[0 mmax]);
% 
% keyboard

%% 
if isequaln(impuse,imp)
  mplt = 4;  
else
  mplt = 2; 
end

% for m = 1:mmax
%   ndpts = size(amppairlump{m}, 1);
%   if ndpts < 100
%     break
%   end
% end
% mplt = m-1;

catmedampbin = cell(mplt,1);  %median amp of the events in a cluster of m evts
catmeddtnnm = cell(mplt,1); %diff time between evts N and N-m for clusters of m evts
catmedamppair = cell(mplt,1); %mean amp of evt pairs, lumped all cats of clusters
catmeddtpair = cell(mplt,1);  %diff time between evt pair N and N-n
catmedampcont = cell(mplt,1); %median amp of all continous evts between the evt pair N and N-n, lumped all cats of clusters
catmeddtcont = cell(mplt,1);  %diff time between evt pair N and N-n

for m=1:mplt
  if isequaln(impuse,imp)
    if m<2
      nbin = 8;
    elseif m<4
      nbin = 5;
    else
      nbin = 3;
    end
  else
    if m<2
      nbin = 5;
    elseif m<4
      nbin = 3;
    end
  end

  %%%bin by the median amp of the events in a cluster of m evts
  aa = catmedamp{m};  %median amp of the events in a cluster of m evts
  medamp = aa(:,1);
  bb = catdtnnm{m}; %diff time between evts N and N-m for clusters of m evts
  dtnnm = bb(:,1);
  [ampbin,indbin] = binxeqnum(log10(medamp),nbin);  %bin by amp with same number
  medampbin = zeros(nbin,1);  %median amp of each amp bin
  meddtnnm = zeros(nbin,1);  %median diff time of each amp bin
  for i = 1: nbin
    indi = indbin{i};
    medampbin(i) = median(ampbin{i});
    meddtnnm(i) = median(dtnnm(indi))/sps;
  end
  catmedampbin{m} = medampbin;
  catmeddtnnm{m} = meddtnnm;
  ndpts(m) = size(ampbin{1},1);
    
  %%%bin by the mean amp of evt pair N and N-n, lumped all cats of clusters
  amppair = amppairlump{m}; %mean amp of evt pairs, lumped all cats of clusters
  dtpair = dtpairlump{m}; %diff time between evt pair N and N-n
  [ampbin,indbin] = binxeqnum(log10(amppair),nbin);  %bin by amp with same number
  medamppair = zeros(nbin,1);  %median amp of each amp bin
  meddtpair = zeros(nbin,1);  %median diff time of each amp bin
  for i = 1: nbin
    indi = indbin{i};
    medamppair(i) = median(ampbin{i});
    meddtpair(i) = median(dtpair(indi))/sps;
  end
  catmedamppair{m} = medamppair;
  catmeddtpair{m} = meddtpair;
  ndptspair(m) = size(ampbin{1},1);
  
  %%%bin by the median amp of all continous evts between the evt pair N and N-n, lumped all cats of clusters
  ampcont = ampcontlump{m}; %median amp of all continous evts between the evt pair N and N-n, lumped all cats of clusters
  dtpair = dtpairlump{m}; %diff time between evt pair N and N-m
  [ampbin,indbin] = binxeqnum(log10(ampcont),nbin);  %bin by amp with same number
  medampcont = zeros(nbin,1);  %median amp of each amp bin
  meddtpair = zeros(nbin,1);  %median diff time of each amp bin
  for i = 1: nbin
    indi = indbin{i};
    medampcont(i) = median(ampbin{i});
    meddtpair(i) = median(dtpair(indi))/sps;
  end
  catmedampcont{m} = medampcont;
  catmeddtcont{m} = meddtpair;
  ndptscont(m) = size(ampbin{1},1);
  
end
%%
nrow = 1; % rows and cols of subplots in each figure
ncol = 3; 
widin = 12; % size of each figure
htin = 5;
pltxran = [0.1 0.95]; pltyran = [0.15 0.95];
pltxsep = 0.1; pltysep = 0.03;
f = initfig(widin,htin,nrow,ncol);
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

color = plasma(mplt);
xran = [-1.2 0.3];
yran = [0.15 0.35];
ax = f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
for m=1:mplt
%   if m==1
%     aa = catmedamp{m};
%     medamp = aa(:,1);
%     bb = catdtnnm{m};
%     dtnnm = bb(:,1);
%     scatter(ax,log10(medamp),dtnnm/sps,15,[.7 .7 .7],'filled');
%   end
%   histogram(ax,log10(medamp),nbin,'Normalization','probability');
  medampbin = catmedampbin{m};
  meddtnnm = catmeddtnnm{m};
  plot(ax,medampbin,meddtnnm-0.25*(m-1),'-','Color',color(m,:),'linew',1);
  scatter(ax,medampbin,meddtnnm-0.25*(m-1),30*sqrt(ndpts(m)/ndpts(1)),color(m,:),...
    'filled','MarkerEdgeColor','k');
end
xlabel(ax,'log_{10}{Med of med amp of all events of the cluster}');
ylabel(ax,'Med of time diff. (s) between N and N-m of the cluster');
ylim(ax,yran);
xlim(ax,xran);

ax = f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
p = [];
label = [];
for m=1:mplt
  medamppair = catmedamppair{m};
  meddtpair = catmeddtpair{m};  
  p(m)=plot(ax,medamppair,meddtpair-0.25*(m-1),'-','Color',color(m,:),'linew',1);
  scatter(ax,medamppair,meddtpair-0.25*(m-1),30*sqrt(ndptspair(m)/ndptspair(1)),color(m,:),...
    'filled','MarkerEdgeColor','k');
  label{m} = sprintf('m=%d',m);
end
xlabel(ax,'log_{10}{Med of mean amp of event pair N and N-m}');
ylabel(ax,'Med of time diff. (s) between event pair N and N-m');
ylim(ax,yran);
xlim(ax,xran);
legend(ax,p,label,'Location','south');

ax = f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
for m=1:mplt
  medampcont = catmedampcont{m};
  meddtpair = catmeddtcont{m};  
  plot(ax,medampcont,meddtpair-0.25*(m-1),'-','Color',color(m,:),'linew',1);
  scatter(ax,medampcont,meddtpair-0.25*(m-1),30*sqrt(ndptscont(m)/ndptscont(1)),color(m,:),...
    'filled','MarkerEdgeColor','k');
end
xlabel(ax,'log_{10}{Med of med amp of all evts between pair N and N-m}');
ylabel(ax,'Med of time diff. (s) between event pair N and N-m');
ylim(ax,yran);
xlim(ax,xran);

keyboard

% %% for a few m, diff time distribution between N&N-m, and fraction w/i dtcut  
% mmax = 15;
% nbst = size(trange,1);
% f1 = initfig(10,4.5,1,2); %initialize fig
% tit=supertit(f1.ax,strcat({'Diff. arrival between N & N-m, '},supertstr));
% movev(tit,0.3);
% [f1,cnt,Nn,Nnn,frac,fracn,mmaxnonzero,mmaxnonzeron]=...
%   plt_srcdifftime_NNm_mmax(f1,nbst,imp,nsrc,impn,nsrcn,mmax,sps,typepltnoi);

% %% for a few m, bin by median amp for all events within N&N-m, fraction of unique events w/i dtcut  
% mmax = 15;
% nbst = size(trange,1);
% scale = 'linear';
% % scale = 'log';
% if strcmp(scale,'linear')
%   widin=10; htin=5; nrow=1; ncol=2;
% elseif strcmp(scale,'log')
%   widin=5; htin=5; nrow=1; ncol=1;
% end
% f2 = initfig(widin,htin,nrow,ncol); %initialize fig
% xran = [0.15 0.9]; yran = [0.15 0.9];
% xsep = 0.05; ysep = 0.05;
% optaxpos(f2,nrow,ncol,xran,yran,xsep,ysep);
% 
% tit=supertit(f2.ax,strcat({'Diff. arrival between N & N-m, binned by amp, '},supertstr));
% movev(tit,0.2);
% [f2,ampbincnt,fracevtb,fracevtbn]=...
%   plt_fracevt_NNm_mmax(f2,nbst,imp,nsrc,impn,nsrcn,mmax,sps,typepltnoi,scale);
% keyboard



%%
% %%%abs distance along min-rmse direction between each source and all others whose arrival 
% %%%separation is <=2 s
% widin = 11;  % maximum width allowed is 8.5 inches
% htin = 5;   % maximum height allowed is 11 inches
% nrow = 1;
% ncol = 2;
% f = initfig(widin,htin,nrow,ncol); %initialize fig
% 
% ax=f.ax(1);
% hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% binwdist = 0.1;
% binedge = (0: binwdist: 50*binwdist)';
% %if looking at distance in map view, use plain count or normalize by area
% %%%projection along the min-rmse direction
% %%%for data
% nsrcprop = nsrc;
% nsrcprop(nsrcprop<=2)=0;
% dlocproj2allbst = [];
% dlocprojnn1bst = [];
% dlocprojnn2bst = [];
% for i = 1: size(trange,1)
%   ist = sum(nsrc(1:i-1))+1;
%   ied = ist+nsrc(i)-1;
%   tarvlsplst = tarvlsplstall(ist:ied);
%   ist = sum(nsrcprop(1:i-1))+1;
%   ied = ist+nsrcprop(i)-1;
%   locxyproj = locxyprojall(ist:ied, :);
%   if ~isempty(tarvlsplst) && ~isempty(locxyproj) 
%     [~,dlocproj2all] = srcdistall(tarvlsplst,locxyproj,[0 2*sps]);
%     dlocproj2allbst = [dlocproj2allbst; dlocproj2all];
%     [~,dlocproj] = srcdistNtoNm(tarvlsplst,locxyproj,2);
%     dlocprojnn1bst = [dlocprojnn1bst; dlocproj{1}];
%     dlocprojnn2bst = [dlocprojnn2bst; dlocproj{2}];
%   end
% end
% dist2all = abs(dlocproj2allbst(:,1));
% [bincnt,binhgt,count,normalizer] = histbinbyarea(dist2all,binedge,'countdensity');
% % binhgt = binhgt./length(dist2all);
% binhgt = count./length(dist2all);
% p1=bar(ax,bincnt,binhgt,1,'stacked','b','facea',0.6);
% % stairs(ax,binedge,[binhgt; binhgt(end)],'k','LineWidth',1);
% plot(ax,[median(dist2all) median(dist2all)],ax.YLim,'b--','linew',1.5);
% text(ax,median(dist2all)+0.1,ax.YLim(2)-0.1*range(ax.YLim),...
%   sprintf('%.2f',median(dist2all)),'HorizontalAlignment','left');
% %%%for noise
% nsrcpropn = nsrcn;
% nsrcpropn(nsrcpropn<=2)=0;
% dlocproj2allbstn = [];
% for i = 1: size(trange,1)
% %   i
%   ist = sum(nsrcn(1:i-1))+1;
%   ied = ist+nsrcn(i)-1;
%   tarvlsplst = tarvlsplstalln(ist:ied);
%   ist = sum(nsrcpropn(1:i-1))+1;
%   ied = ist+nsrcpropn(i)-1;
%   locxyproj = locxyprojalln(ist:ied, :);
%   if ~isempty(tarvlsplst) && ~isempty(locxyproj) 
%     [~,dlocproj2all] = srcdistall(tarvlsplst,locxyproj,[0 2*sps]);
%     dlocproj2allbstn = [dlocproj2allbstn; dlocproj2all];
%   end
% end
% dist2alln = abs(dlocproj2allbstn(:,1));
% [bincnt,binhgt,count,normalizer] = histbinbyarea(dist2alln,binedge,'countdensity');
% % binhgt = binhgt./length(dist2all);
% binhgt = count./length(dist2all);
% p2=bar(ax,bincnt,binhgt,1,'stacked','r','facea',0.6);
% % stairs(ax,binedge,[binhgt; binhgt(end)],'k','LineWidth',1);
% plot(ax,[median(dist2alln) median(dist2alln)],ax.YLim,'r--','linew',1.5);
% text(ax,median(dist2alln)+0.1,ax.YLim(2)-0.2*range(ax.YLim),...
%   sprintf('%.2f',median(dist2alln)),'HorizontalAlignment','left');
% xlabel(ax,'Dist. (km) along min-scatter direc. between each source and all others with diff. arrival time \leq 2 s');
% ylabel(ax,'Normalized count');
% legend(ax,[p1,p2],'Data','Synthetic noise','location','east');
% xlim(ax,[0 4]);
% longticks(ax,2);
% hold(ax,'off');
% 
% ax=f.ax(2);
% hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% binwdt = 0.05;
% [N,edges]=histcounts(dtarvlplt/sps,'binwidth',binwdt,'normalization','count');
% N = N./length(dtarvlplt);
% bincnt = (edges(1:end-1)+edges(2:end))/2;
% p1=bar(ax,bincnt,N,1,'stacked','b','facea',0.6,'edgecolor','k');
% % p1=histogram(ax,dtarvlnn1/sps,'binwidth',binwdt,'normalization','count',...
% %   'facecolor','b','EdgeColor','k','facea',0.6);
% % p1.Values = p1.Values ./ length(dtarvlnn1);
% plot(ax,[median(dtarvlplt/sps) median(dtarvlplt/sps)],ax.YLim, '--', ...
%   'Color', 'b', 'linew', 1.5);
% text(ax,median(dtarvlplt/sps)+0.02,ax.YLim(2)-0.1*range(ax.YLim),...
%   sprintf('%.2f',median(dtarvlplt/sps)),'HorizontalAlignment','left');
% 
% [N,edges]=histcounts(dtarvlpltn/sps,'binwidth',binwdt,'normalization','count');
% N = N./length(dtarvlplt);
% bincnt = (edges(1:end-1)+edges(2:end))/2;
% p2=bar(ax,bincnt,N,1,'stacked','r','facea',0.6,'edgecolor','k');
% plot(ax,[median(dtarvlpltn/sps) median(dtarvlpltn/sps)],ax.YLim, '--', ...
%   'Color', 'r', 'linew', 1.5);
% text(ax,median(dtarvlpltn/sps)+0.02,ax.YLim(2)-0.2*range(ax.YLim),...
%   sprintf('%.2f',median(dtarvlpltn/sps)),'HorizontalAlignment','left');
% ylabel(ax,'Normalized count');
% xlabel(ax,sprintf('Diff. arrival time between sources N and N-%d (s)',nsep));
% legend(ax,[p1,p2],'Data','Synthetic noise','location','east');
% xlim(ax,[0 1.5]);
% % ylim(ax,[0 3.5]);
% longticks(ax,2);
% hold(ax,'off');
% tit=supertit(f.ax,supertstr);
% movev(tit,0.2);
% 
% orient(f.fig,'landscape');
% fname = strcat(sprintf('nn%ddiff',nsep),fnsuffix,'.pdf');
% print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));

% keyboard

%% amp, but only for N & N-m pair, ignoring sources in between
% ampafdt = [];
% ampbfdt = [];
% for i = 1: size(trange,1)
%   ist = sum(nsrc(1:i-1))+1;
%   ied = ist+nsrc(i)-1;
%   impi = imp(ist:ied,:);
%   impaf = impi(1+nsep:end,:);
%   ampaf = mean(impaf(:,[2 4 6]),2);
%   ampafdt = [ampafdt; ampaf];
%   impbf = impi(1:end-nsep,:);
%   ampbf = mean(impbf(:,[2 4 6]),2);
%   ampbfdt = [ampbfdt; ampbf];
% end
% ampdt = mean([ampbfdt ampafdt],2);
% minnum = 500;
% 
% ampafdtn = [];
% ampbfdtn = [];
% for i = 1: size(trange,1)
%   ist = sum(nsrcn(1:i-1))+1;
%   ied = ist+nsrcn(i)-1;
%   impi = impn(ist:ied,:);
%   impaf = impi(1+nsep:end,:);
%   ampaf = mean(impaf(:,[2 4 6]),2);
%   ampafdtn = [ampafdtn; ampaf];
%   impbf = impi(1:end-nsep,:);
%   ampbf = mean(impbf(:,[2 4 6]),2);
%   ampbfdtn = [ampbfdtn; ampbf];
% end
% ampdtn = mean([ampbfdtn ampafdtn],2);
% minnumn = 250;
% 
% ampplt = ampdt;
% amppltn = ampdtn;

% figure
% subplot(311)
% h=histogram(log10(ampbfdt),'BinWidth',0.25); hold on; %,'NumBins',5
% ax=gca;
% plot(ax,ax.XLim,[200 200],'k--');
% xlabel('log_{10}{Amp}');
% ylabel('Count');
% xlim([-1.5 1]);
% text(0.95,0.9,'Earlier one of the pair','Units','normalized','HorizontalAlignment','right');
% title(sprintf('Between sources N and N-%d (s)',nsep));
% subplot(312)
% h=histogram(log10(ampafdt),'BinWidth',0.25); hold on;
% ax=gca;
% plot(ax,ax.XLim,[200 200],'k--');
% xlim([-1.5 1]);
% text(0.95,0.9,'Later one of the pair','Units','normalized','HorizontalAlignment','right');
% subplot(313)
% h=histogram(log10(ampdt),'BinWidth',0.25); hold on;
% ax=gca;
% plot(ax,ax.XLim,[200 200],'k--');
% xlim([-1.5 1]);
% text(0.95,0.9,'Mean of the pair','Units','normalized','HorizontalAlignment','right');


%%
% figure
% subplot(121)
% eucdistnn1 = allbstsig.distarvlnn1all(:,1);
% eucdistnn2 = allbstsig.distarvlnn2all(:,1);
% 
% eucdistplt = eucdistnn1;
% dlocprojplt = dlocprojnn1bst;
% 
% ax = gca;
% hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% tmp = [dtarvlplt/sps eucdistplt log10(ampplt)];
% tmp = sortrows(tmp,3,'descend');
% scatter(ax,tmp(:,1),tmp(:,2),15,tmp(:,3),'o','filled');
% oldc = colormap(ax,'kelicol');
% newc = flipud(oldc);
% colormap(ax,newc);
% c=colorbar(ax,'SouthOutside');
% c.Label.String = strcat({'log_{10}(amp)'});
% ylabel(ax,'Distance (km)');
% xlabel(ax,sprintf('Diff. arrival time between sources N and N-%d (s)',nsep));
% xlim(ax,[0 2]);
% ylim(ax,[0 8]);
% hold(ax,'off');
% 
% subplot(122)
% ax = gca;
% hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% tmp = [dtarvlplt/sps abs(dlocprojplt(:,1)) log10(ampplt)];
% tmp = sortrows(tmp,3,'descend');
% scatter(ax,tmp(:,1),tmp(:,2),15,tmp(:,3),'o','filled');
% oldc = colormap(ax,'kelicol');
% newc = flipud(oldc);
% colormap(ax,newc);
% c=colorbar(ax,'SouthOutside');
% c.Label.String = strcat({'log_{10}(amp)'});
% ylabel(ax,'Distance along min-rmse (km)');
% xlabel(ax,sprintf('Diff. arrival time between sources N and N-%d (s)',nsep));
% xlim(ax,[0 2]);
% ylim(ax,[0 4]);
% hold(ax,'off');

