function [clusibst,nclusibst,tclusibst]=clusters_in_burst(idxbst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is made to associate all types of event clusters in the
% deconvolution catalog to each tremor burst. The existing code 'evtcluster_ex.m'
% finds clusters of 'm' and then categorize them according to 'm'. We want to 
% look at in a different way. For each burst, we want to find all types of
% clusters.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/03/06
% Last modified date:   2024/03/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('idxbst',[]);

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
% format short e   % Set the format to 5-digit floating point
% clear
% clc
% close all

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

%% during example burst 181, are there clusters?
% clusbst = cell(nbst,1); %clusters of different m in each burst
% nclusbst = cell(nbst,1);  %num of clusters of different m in each burst
% tclusbst = cell(nbst,1);  %start and end time of each cluster in each burst

% for idxbst = 181:nbst
%   idxbst = 181;
  
  ist = sum(nsrcuse(1:idxbst-1))+1;
  ied = ist+nsrcuse(idxbst)-1;
  impi = impuse(ist:ied,:);
  timp = impi(:,1)/sps;
  
  clusibst = [];
  nclusibst = zeros(mmax,1);
  tclusibst = [];
  k = 0;
  for m=1:mmax
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
      clusibst{k,1} = catclusibstm{j};
    end
  end
  
  for i=1:k
    aaa = clusibst{i};
    timpi = aaa(:,1)/sps;
    tclusibst(i,1) = timpi(1);
    tclusibst(i,2) = timpi(end);
    tclusibst(i,3) = timpi(end)-timpi(1);
    tclusibst(i,4) = length(timpi);
  end
  
  [tmp,ind] = sort(tclusibst(:,1));
  tclusibst = tclusibst(ind,:);
  
  tmp=[];
  for i=1:k
    tmp{i,1} = clusibst{ind(i)};
  end
  clusibst = tmp;
  
%   clusbst{idxbst} = clusibst;
%   nclusbst{idxbst} = nclusibst;
%   tclusbst{idxbst} = tclusibst;
  
  for i=1:k-1
    if tclusibst(i,2)>tclusibst(i+1,1)
      disp('two clusters overlaps in time!');
    end
  end
  % keyboard
  
  %which clusters are within a zoom-in window
  xzoom = [5 30];
  ind = find(tclusibst(:,1)>=xzoom(1) & tclusibst(:,2)<=xzoom(2));
  clusibstxzoom = clusibst(ind);
  tclusibstxzoom = tclusibst(ind);
  
% end





