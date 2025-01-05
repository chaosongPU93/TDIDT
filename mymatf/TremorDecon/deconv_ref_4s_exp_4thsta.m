% deconv_ref_4s_exp_4thsta.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the 'driver' script to call 'deconv_ref_4s_exp_4thsta_fn' to do the
% analysis based on the desired 'indofburst'. 
% so that outputs from either from quieter burst windows (presumbly night
% times) or from noisier windows presumbly day times), could be directly
% compared in one single plot or more. Here i am more interested into how
% the signal waveform coherency between 4th stations and trio stations changes
% wrt to the background noise level (night or day). Maybe the prediction at
% a 4th station would be better if the burst was during night times given a 
% set of deconvolved sources using trio stations.
% --Instead of asking the distance from each src to all others within 2s,
% ie, close-in-time sources, we now move on to ask within 10s (20-s range),
% not only the distance, but also the number of srcs matter. 20-s range 
% is able to constrain the migration to be less than 100m, 20% of the 
% distance, if there is any migration. --- 2024/01/17
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/10/13
% Last modified date:   2022/10/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
% format short e   % Set the format to 5-digit floating point
clear
clc
close all

% set(0,'DefaultFigureVisible','on');
set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

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
% trange = trange(1:end-1,:);
tlen = trange(:,3)-trange(:,2);
tlensum = sum(tlen);

modname = 'timeoff_plfit_4thsta_160sps.mat';
planefit = load(strcat(rstpath, '/MAPS/',modname));
rmse = planefit.gof{4}.rmse;
offmax = round(2.0*rmse);

%% some options for burst windows depending on the assumed noise level 
%empirically determined indices of bursts fall into local day times (noisier)
%or night times (quieter)
inbst = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,...
  34,35,56,57,58,59,60,61,62,70,71,72,73,74,75,76,77,78,79,80,81,82,110,111,112,113,114,115,116,117,...
  118,119,120,142,143,144,145,146,147,148,149,150,151,152,153,154,172,173,174,175]';
idbst = [36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,63,64,65,66,67,68,69,83,84,85,...
  86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,121,122,123,124,...
  125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,155,156,157,158,159,160,161,...
  162,163,164,165,166,167,168,169,170,171,176,177,178,179,180,181,182,183,184,185,186,187,188,189,...
  190,191,192,193,194,195]';

%indices of bursts whose CC of sig-wlet cc between sta 1 and 4 are high, above 80th percentile
ind41 = [1,4,6,10,12,15,21,22,23,24,29,46,52,56,74,77,78,82,83,102,111,113,114,121,124,125,...
  129,145,151,152,155,166,167,169,174,175,180,183,187]';

%indices of bursts whose CC of sig-wlet cc between sta 14, 24, 34 are all high, above 80th percentile
ind4123 = [1,12,15,46,56,77,102,114,125,129,155,180]';

%indices of bursts whose CC of envelope between sta 1,2,3 are high, above 75th percentile
indenv123 = [1,3,8,10,31,67,78,81,82,91,102,108,114,116,121,129,146,153,167,175];

%indices of bursts whose CC of sig-wlet cc between sta 1,2,3 are high, above 75th percentile
indswcc123 = [1,2,3,6,8,18,21,46,56,77,78,82,83,84,102,107,114,125,127,145,169];

sps = 160;

%% call function 'deconv_ref_4s_exp_rand_fn', high-correlation VS low-correlation bursts, DATA VS NOISE
% flagrecalc = 0;
% % flagrecalc = 1;
% 
% if flagrecalc
%   normflag = 0; %do not normalize the templates
%   pltflag = 0;  %do not create summary plots for each choice of inputs
%   % rccmwsec = 0.25; %use 0.5s or 0.25s 
%   rccmwsec = 0.5; %use 0.5s or 0.25s
%   
%   %%%high-correlation bursts using real data
%   noiseflag = 0;
%   hibstsig = deconv_ref_4s_exp_4thsta_fn(indenv123,normflag,noiseflag,pltflag,rccmwsec);
%   savefile = 'deconv_stats4th_hienvccbstsig.mat';
%   save(strcat(rstpath, '/MAPS/',savefile), 'hibstsig');
%   
%   %%%high-correlation bursts using synthetic noise
%   noiseflag = 1;
%   hibstnoi = deconv_ref_4s_exp_4thsta_fn(indenv123,normflag,noiseflag,pltflag,rccmwsec);
%   savefile = 'deconv_stats4th_hienvccbstnoi.mat';
%   save(strcat(rstpath, '/MAPS/',savefile), 'hibstnoi');
% 
% else
%   savefile = 'deconv_stats4th_hienvccbstsig.mat';
%   load(strcat(rstpath, '/MAPS/',savefile));
%   savefile = 'deconv_stats4th_hienvccbstnoi.mat';
%   load(strcat(rstpath, '/MAPS/',savefile));
% end
% 
% bstsig = hibstsig;
% bstnoi = hibstnoi;

% keyboard

%% call function 'deconv_ref_4s_exp_rand_fn', all bursts, DATA VS NOISE
%%%Flag to indicate if it is necessary to recalculate everything
flagrecalc = 0;
% flagrecalc = 1;

if flagrecalc
  normflag = 0; %do not normalize the templates
  pltflag = 0;  %do not create summary plots for each burst
%   rccmwsec = 0.25; %use 0.5s or 0.25s 
  rccmwsec = 0.5; %use 0.5s or 0.25s
  
  %%all bursts using real data
  noiseflag = 0;
  allbstsig = deconv_ref_4s_exp_4thsta_fn(1:size(trange,1),normflag,noiseflag,pltflag,rccmwsec); %
  % savefile = 'deconv_stats4th_allbstsig.mat'; 
  savefile = 'deconv_stats4th_no23_allbstsig.mat';
%   savefile = 'deconv_stats4th_no23_allbstsig0.25s.mat';
  save(strcat(rstpath, '/MAPS/',savefile), 'allbstsig');

  %%%all bursts using synthetic noise
  noiseflag = 1;
  allbstnoi = deconv_ref_4s_exp_4thsta_fn(1:size(trange,1),normflag,noiseflag,pltflag,rccmwsec);  %
  % savefile = 'deconv_stats4th_allbstnoi.mat';
  savefile = 'deconv_stats4th_no23_allbstnoi.mat';
%   savefile = 'deconv_stats4th_no23_allbstnoi0.25s.mat';
  save(strcat(rstpath, '/MAPS/',savefile), 'allbstnoi');
  
else
%   savefile = 'deconv_stats4th_allbstsig.mat';
  savefile = 'deconv_stats4th_no23_allbstsig.mat';
  load(strcat(rstpath, '/MAPS/',savefile));
%   savefile = 'deconv_stats4th_allbstnoi.mat';
  savefile = 'deconv_stats4th_no23_allbstnoi.mat';
  load(strcat(rstpath, '/MAPS/',savefile));
end

keyboard

% bstsig = allbstsig;
% bstnoi = allbstnoi;

% %%
% idxbst = 1:size(trange,1);
% for iii = 1: length(idxbst)  
%   f1name{idxbst(iii),1} = sprintf('shortwinmapprojboth_bst%s.pdf',...
%     num2zeropadstr(idxbst(iii), 3));
% end
% for i = 1:size(f1name,1)
%   fname{i} = fullfile(rstpath, '/FIGS/',f1name{i});
% end
% 
% status = system('rm -f /home/data2/chaosong/matlab/allan/data-no-resp/PGCtrio/FIGS/shortwinmapprojboth.pdf');
% append_pdfs('/home/data2/chaosong/matlab/allan/data-no-resp/PGCtrio/FIGS/shortwinmapprojboth.pdf',fname{:});
% status = system('rm -f /home/data2/chaosong/matlab/allan/data-no-resp/PGCtrio/FIGS/shortwinmapprojboth_bst*.pdf');
        
%% call function 'deconv_ref_4s_exp_rand_fn', all bursts, DATA VS NOISE, use 0.25-s win for RCC
%%%Flag to indicate if it is necessary to recalculate everything
% flagrecalc = 0;
% % flagrecalc = 1;
% 
% if flagrecalc
%   normflag = 0; %do not normalize the templates
%   pltflag = 0;  %do not create summary plots for each choice of inputs
%   rccmwsec = 0.25; %use 0.5s or 0.25s 
% %   rccmwsec = 0.5; %use 0.5s or 0.25s
%   
%   %%%all bursts using real data
%   noiseflag = 0;
%   allbstsig = deconv_ref_4s_exp_4thsta_fn(1:size(trange,1),normflag,noiseflag,pltflag,rccmwsec); %
%   savefile = 'deconv_stats4th_allbstsig0.25s.mat';
%   save(strcat(rstpath, '/MAPS/',savefile), 'allbstsig');
%   
%   %%%all bursts using synthetic noise
%   noiseflag = 1;
%   allbstnoi = deconv_ref_4s_exp_4thsta_fn(1:size(trange,1),normflag,noiseflag,pltflag,rccmwsec);  %
%   savefile = 'deconv_stats4th_allbstnoi0.25s.mat';
%   save(strcat(rstpath, '/MAPS/',savefile), 'allbstnoi');
% else
%   savefile = 'deconv_stats4th_allbstsig0.25s.mat';
%   load(strcat(rstpath, '/MAPS/',savefile));
%   savefile = 'deconv_stats4th_allbstnoi0.25s.mat';
%   load(strcat(rstpath, '/MAPS/',savefile));
% end
% 
% bstsig = allbstsig;
% bstnoi = allbstnoi;
% keyboard
 
% params for loading data
%%%param for secondary sources removed
%data
locxyprojall = allbstsig.locxyprojall;
tarvlsplstall = allbstsig.impindepall(:,1);
nsrc = allbstsig.nsrc;
dtarvlnn1 = allbstsig.dtarvlnn1all;
dtarvlnn2 = allbstsig.dtarvlnn2all;
dtarvlnn3 = allbstsig.dtarvlnn3all;
dtarvlnn4 = allbstsig.dtarvlnn4all;
dtarvlnn5 = allbstsig.dtarvlnn5all;
imp = allbstsig.impindepall;
dt2all = allbstsig.dt2allbst;
dist2all = allbstsig.dist2allbst;
dloc2all = allbstsig.dloc2allbst;
dneloc{1,1} = allbstsig.distarvlnn1all(:,2:3);
dneloc{2,1} = allbstsig.distarvlnn2all(:,2:3);
dneloc{3,1} = allbstsig.distarvlnn3all(:,2:3);
dneloc{4,1} = allbstsig.distarvlnn4all(:,2:3);
dneloc{5,1} = allbstsig.distarvlnn5all(:,2:3);
eucdist{1,1} = allbstsig.distarvlnn1all(:,1);
eucdist{2,1} = allbstsig.distarvlnn2all(:,1);
eucdist{3,1} = allbstsig.distarvlnn3all(:,1);
eucdist{4,1} = allbstsig.distarvlnn4all(:,1);
eucdist{5,1} = allbstsig.distarvlnn5all(:,1);
dtarvl{1,1} = allbstsig.dtarvlnn1all;
dtarvl{2,1} = allbstsig.dtarvlnn2all;
dtarvl{3,1} = allbstsig.dtarvlnn3all;
dtarvl{4,1} = allbstsig.dtarvlnn4all;
dtarvl{5,1} = allbstsig.dtarvlnn5all;
srcamprall = allbstsig.srcamprall;
%noise
locxyprojalln = allbstnoi.locxyprojall;
tarvlsplstalln = allbstnoi.impindepall(:,1);
nsrcn = allbstnoi.nsrc;
dtarvlnn1n = allbstnoi.dtarvlnn1all;
dtarvlnn2n = allbstnoi.dtarvlnn2all;
dtarvlnn3n = allbstnoi.dtarvlnn3all;
dtarvlnn4n = allbstnoi.dtarvlnn4all;
dtarvlnn5n = allbstnoi.dtarvlnn5all;
impn = allbstnoi.impindepall;
dt2alln = allbstnoi.dt2allbst;
dist2alln = allbstnoi.dist2allbst;
dloc2alln = allbstnoi.dloc2allbst;
dnelocn{1,1} = allbstnoi.distarvlnn1all(:,2:3);
dnelocn{2,1} = allbstnoi.distarvlnn2all(:,2:3);
dnelocn{3,1} = allbstnoi.distarvlnn3all(:,2:3);
dnelocn{4,1} = allbstnoi.distarvlnn4all(:,2:3);
dnelocn{5,1} = allbstnoi.distarvlnn5all(:,2:3);
eucdistn{1,1} = allbstnoi.distarvlnn1all(:,1);
eucdistn{2,1} = allbstnoi.distarvlnn2all(:,1);
eucdistn{3,1} = allbstnoi.distarvlnn3all(:,1);
eucdistn{4,1} = allbstnoi.distarvlnn4all(:,1);
eucdistn{5,1} = allbstnoi.distarvlnn5all(:,1);
dtarvln{1,1} = allbstnoi.dtarvlnn1all;
dtarvln{2,1} = allbstnoi.dtarvlnn2all;
dtarvln{3,1} = allbstnoi.dtarvlnn3all;
dtarvln{4,1} = allbstnoi.dtarvlnn4all;
dtarvln{5,1} = allbstnoi.dtarvlnn5all;
srcampralln = allbstnoi.srcamprall;
supertstr = 'Secondary sources removed';
fnsuffix = [];
nsta = 3;

% %%%param for further checked at KLNB
% %data
% locxyprojall = allbstsig.locxyproj4thall;
% tarvlsplstall = allbstsig.impindep4thall(:,1);
% nsrc = allbstsig.nsrc4th;
% dtarvlnn1 = allbstsig.dtarvlnn14thall;
% dtarvlnn2 = allbstsig.dtarvlnn24thall;
% dtarvlnn3 = allbstsig.dtarvlnn34thall;
% imp = allbstsig.impindep4thall;
% dt2all = allbstsig.dt2all4thbst;
% dist2all = allbstsig.dist2all4thbst;
% dloc2all = allbstsig.dloc2all4thbst;
% dneloc{1,1} = allbstsig.distarvlnn14thall(:,2:3);
% dneloc{2,1} = allbstsig.distarvlnn24thall(:,2:3);
% dneloc{3,1} = allbstsig.distarvlnn34thall(:,2:3);
% dneloc{4,1} = allbstsig.distarvlnn44thall(:,2:3);
% dneloc{5,1} = allbstsig.distarvlnn54thall(:,2:3);
% eucdist{1,1} = allbstsig.distarvlnn14thall(:,1);
% eucdist{2,1} = allbstsig.distarvlnn24thall(:,1);
% eucdist{3,1} = allbstsig.distarvlnn34thall(:,1);
% eucdist{4,1} = allbstsig.distarvlnn44thall(:,1);
% eucdist{5,1} = allbstsig.distarvlnn54thall(:,1);
% dtarvl{1,1} = allbstsig.dtarvlnn14thall;
% dtarvl{2,1} = allbstsig.dtarvlnn24thall;
% dtarvl{3,1} = allbstsig.dtarvlnn34thall;
% dtarvl{4,1} = allbstsig.dtarvlnn44thall;
% dtarvl{5,1} = allbstsig.dtarvlnn54thall;
% srcamprall = allbstsig.srcampr4thall;
% %noise
% locxyprojalln = allbstnoi.locxyproj4thall;
% tarvlsplstalln = allbstnoi.impindep4thall(:,1);
% nsrcn = allbstnoi.nsrc4th;
% dtarvlnn1n = allbstnoi.dtarvlnn14thall;
% dtarvlnn2n = allbstnoi.dtarvlnn24thall;
% dtarvlnn3n = allbstnoi.dtarvlnn34thall;
% impn = allbstnoi.impindep4thall;
% dt2alln = allbstnoi.dt2all4thbst;
% dist2alln = allbstnoi.dist2all4thbst;
% dloc2alln = allbstnoi.dloc2all4thbst;
% dnelocn{1,1} = allbstnoi.distarvlnn14thall(:,2:3);
% dnelocn{2,1} = allbstnoi.distarvlnn24thall(:,2:3);
% dnelocn{3,1} = allbstnoi.distarvlnn34thall(:,2:3);
% dnelocn{4,1} = allbstnoi.distarvlnn44thall(:,2:3);
% dnelocn{5,1} = allbstnoi.distarvlnn54thall(:,2:3);
% eucdistn{1,1} = allbstnoi.distarvlnn14thall(:,1);
% eucdistn{2,1} = allbstnoi.distarvlnn24thall(:,1);
% eucdistn{3,1} = allbstnoi.distarvlnn34thall(:,1);
% eucdistn{4,1} = allbstnoi.distarvlnn44thall(:,1);
% eucdistn{5,1} = allbstnoi.distarvlnn54thall(:,1);
% dtarvln{1,1} = allbstnoi.dtarvlnn14thall;
% dtarvln{2,1} = allbstnoi.dtarvlnn24thall;
% dtarvln{3,1} = allbstnoi.dtarvlnn34thall;
% dtarvln{4,1} = allbstnoi.dtarvlnn44thall;
% dtarvln{5,1} = allbstnoi.dtarvlnn54thall;
% srcampralln = allbstnoi.srcampr4thall;
% supertstr = 'Further checked at KLNB';
% fnsuffix = '4th';
% nsta = 4;

dtarvlplt = dtarvlnn1;
dtarvlpltn = dtarvlnn1n;
mmax = 5;
m = 1;
tmax2all = 10;  %max time in sec from each src to all others to check distance
% dtarvlplt = dtarvlnn2;
% dtarvlpltn = dtarvlnn2n;
% nsep = 2;
% dtarvlplt = dtarvlnn3;
% dtarvlpltn = dtarvlnn3n;
% nsep = 3;

% keyboard

%%
%%%abs distance along min-rmse direction between each source and all others whose arrival 
%%%separation is <=2 s
widin = 11;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f = initfig(widin,htin,nrow,ncol); %initialize fig

ax=f.ax(1);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
binwdist = 0.1;
binedge = (0: binwdist: 50*binwdist)';
%if looking at distance in map view, use plain count or normalize by area
%%%projection along the min-rmse direction
%%%for data
nsrcprop = nsrc;
nsrcprop(nsrcprop<=2)=0;
dlocproj2allbst = [];
dlocprojnn1bst = [];
dlocprojnn2bst = [];
for i = 1: size(trange,1)
  ist = sum(nsrc(1:i-1))+1;
  ied = ist+nsrc(i)-1;
  tarvlsplst = tarvlsplstall(ist:ied);
  ist = sum(nsrcprop(1:i-1))+1;
  ied = ist+nsrcprop(i)-1;
  locxyproj = locxyprojall(ist:ied, :);
  if ~isempty(tarvlsplst) && ~isempty(locxyproj) 
    [~,dlocproj2all] = srcdistall(tarvlsplst,locxyproj,[0 tmax2all*sps]);
    dlocproj2allbst = [dlocproj2allbst; dlocproj2all];
    [~,dlocproj] = srcdistNtoNm(tarvlsplst,locxyproj,2);
    dlocprojnn1bst = [dlocprojnn1bst; dlocproj{1}];
    dlocprojnn2bst = [dlocprojnn2bst; dlocproj{2}];
  end
end
distpjnn1 = abs(dlocprojnn1bst(:,1));
distpj2all = abs(dlocproj2allbst(:,1));
[bincnt,binhgt,count,normalizer] = histbinbyarea(distpj2all,binedge,'countdensity');
% binhgt = binhgt./length(dist2all);
binhgt = count./length(distpj2all);
p1=bar(ax,bincnt,binhgt,1,'stacked','b','facea',0.6);
% stairs(ax,binedge,[binhgt; binhgt(end)],'k','LineWidth',1);
plot(ax,[median(distpj2all) median(distpj2all)],ax.YLim,'b--','linew',1.5);
text(ax,median(distpj2all)+0.1,ax.YLim(2)-0.1*range(ax.YLim),...
  sprintf('%.2f',median(distpj2all)),'HorizontalAlignment','left');
%%%for noise
nsrcpropn = nsrcn;
nsrcpropn(nsrcpropn<=2)=0;
dlocproj2allbstn = [];
for i = 1: size(trange,1)
%   i
  ist = sum(nsrcn(1:i-1))+1;
  ied = ist+nsrcn(i)-1;
  tarvlsplst = tarvlsplstalln(ist:ied);
  ist = sum(nsrcpropn(1:i-1))+1;
  ied = ist+nsrcpropn(i)-1;
  locxyproj = locxyprojalln(ist:ied, :);
  if ~isempty(tarvlsplst) && ~isempty(locxyproj) 
    [~,dlocproj2all] = srcdistall(tarvlsplst,locxyproj,[0 tmax2all*sps]);
    dlocproj2allbstn = [dlocproj2allbstn; dlocproj2all];
  end
end
distpj2alln = abs(dlocproj2allbstn(:,1));
[bincnt,binhgt,count,normalizer] = histbinbyarea(distpj2alln,binedge,'countdensity');
% binhgt = binhgt./length(distpj2all);
binhgt = count./length(distpj2all);
p2=bar(ax,bincnt,binhgt,1,'stacked','r','facea',0.6);
% stairs(ax,binedge,[binhgt; binhgt(end)],'k','LineWidth',1);
plot(ax,[median(distpj2alln) median(distpj2alln)],ax.YLim,'r--','linew',1.5);
text(ax,median(distpj2alln)+0.1,ax.YLim(2)-0.2*range(ax.YLim),...
  sprintf('%.2f',median(distpj2alln)),'HorizontalAlignment','left');
xlabel(ax,'Dist. (km) along min-scatter direc. between each source and all others with diff. arrival time \leq 2 s');
ylabel(ax,'Normalized count');
legend(ax,[p1,p2],'Data','Synthetic noise','location','east');
xlim(ax,[0 4]);
longticks(ax,2);
hold(ax,'off');

ax=f.ax(2);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
binwdt = 0.05;
[N,edges]=histcounts(dtarvlplt/sps,'binwidth',binwdt,'normalization','count');
N = N./length(dtarvlplt);
bincnt = (edges(1:end-1)+edges(2:end))/2;
p1=bar(ax,bincnt,N,1,'stacked','b','facea',0.6,'edgecolor','k');
% p1=histogram(ax,dtarvlnn1/sps,'binwidth',binwdt,'normalization','count',...
%   'facecolor','b','EdgeColor','k','facea',0.6);
% p1.Values = p1.Values ./ length(dtarvlnn1);
plot(ax,[median(dtarvlplt/sps) median(dtarvlplt/sps)],ax.YLim, '--', ...
  'Color', 'b', 'linew', 1.5);
text(ax,median(dtarvlplt/sps)+0.02,ax.YLim(2)-0.1*range(ax.YLim),...
  sprintf('%.2f',median(dtarvlplt/sps)),'HorizontalAlignment','left');

[N,edges]=histcounts(dtarvlpltn/sps,'binwidth',binwdt,'normalization','count');
N = N./length(dtarvlplt);
bincnt = (edges(1:end-1)+edges(2:end))/2;
p2=bar(ax,bincnt,N,1,'stacked','r','facea',0.6,'edgecolor','k');
plot(ax,[median(dtarvlpltn/sps) median(dtarvlpltn/sps)],ax.YLim, '--', ...
  'Color', 'r', 'linew', 1.5);
text(ax,median(dtarvlpltn/sps)+0.02,ax.YLim(2)-0.2*range(ax.YLim),...
  sprintf('%.2f',median(dtarvlpltn/sps)),'HorizontalAlignment','left');
ylabel(ax,'Normalized count');
xlabel(ax,sprintf('Diff. arrival time between sources N and N-%d (s)',m));
legend(ax,[p1,p2],'Data','Synthetic noise','location','east');
xlim(ax,[0 1.5]);
% ylim(ax,[0 3.5]);
longticks(ax,2);
hold(ax,'off');
tit=supertit(f.ax,supertstr);
movev(tit,0.2);

% orient(f.fig,'landscape');
% fname = strcat(sprintf('nn%ddiff',nsep),fnsuffix,'.pdf');
% print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));

% keyboard

%% in map view, 2ndary sources removed
% f=plt_srcdlocinmap(dt2all,dloc2all,sps,'km');

dto2allbst=[];  %diff. origin time between event pairs
dloco2allbst=[];  
disto2allbst=[];
impstall = [];
tcorall = [];
for i = 1: size(trange,1)
  ist = sum(nsrc(1:i-1))+1;
  ied = ist+nsrc(i)-1;
  impi = imp(ist:ied,:);
  %convert time offset to relative loc
  ftrans = 'interpchao';
  [imploc, indinput] = off2space002(impi(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
  [imploc0, ~] = off2space002([0 0],sps,ftrans,0);  % a ref source at 0,0
  tcor = round((imploc(:,6)-imploc0(6))*sps);   % travel time difference between each source and a ref source at 0,0
  torispl = impi(:,1)-tcor;    
  [torisplst, indsort] = sortrows(torispl,1); %note, sort by origin time now
  implocst = imploc(indsort, :);
  impst = impi(indsort, :);
  %I suspect the diff loc and distance between pairs should be the same as using tarvl, but correspondance with time 
  %should be different
  [dto2all,dloco2all,disto2all] = srcdistall(torisplst,implocst,[0 tmax2all*sps]);
  dto2allbst = [dto2allbst; dto2all];
  dloco2allbst = [dloco2allbst; dloco2all] ;
  disto2allbst = [disto2allbst; disto2all];
  impstall = [impstall; impst];
  tcorall = [tcorall; tcor(indsort)];
end
cstr = {'# events / grid'};
f=plt_srcdlocinmap(dto2allbst/sps,dloco2allbst,[],'km','tori',...
  10,cstr,'s','log10','grid',[-4 4],[-4 4],0.1,0.1);

%plot euclidean distance between each LFE source to all others
tlensec = 14981;
f = plt_srcdistall(dto2allbst,disto2allbst,sps,40/sps,tlensec,0.1,'km');
%plot the loc diff between each LFE source to all others
f = plt_srcdlocall(dloco2allbst,0.1,'km');

%plot the loc diff between above source pairs
f = plt_srcdlocNtoNm(dneloc,0.1,'km');
%plot the diff time and distance between above source pairs
f = plt_srcdistNtoNm(dtarvl,eucdist,sps,40/sps,tlensec,0.1,'km');

%% projection along the min-rmse direction
dlocproj2allbst = [];
dto2allbst = [];
tmp1=[];
tmp2=[];
tmp3=[];
tmp4=[];
tmp5=[];
for i = 1: size(trange,1)
  ist = sum(nsrc(1:i-1))+1;
  ied = ist+nsrc(i)-1;
  tarvlsplst = tarvlsplstall(ist:ied);
  locxyproj = locxyprojall(ist:ied, :);
  [dto2all,dlocproj2all] = srcdistall(tarvlsplst,locxyproj,[0 tmax2all*sps]);
  dto2allbst = [dto2allbst; dto2all];
  dlocproj2allbst = [dlocproj2allbst; dlocproj2all];
  
  [~,dlocproj] = srcdistNtoNm(tarvlsplst,locxyproj,mmax);
  tmp1 = [tmp1; dlocproj{1}];
  tmp2 = [tmp2; dlocproj{2}];
  tmp3 = [tmp3; dlocproj{3}];
  tmp4 = [tmp4; dlocproj{4}];
  tmp5 = [tmp5; dlocproj{5}];
end
dlocprojbst{1,1} = tmp1;
dlocprojbst{2,1} = tmp2;
dlocprojbst{3,1} = tmp3;
dlocprojbst{4,1} = tmp4;
dlocprojbst{5,1} = tmp5;

cstr = {'# events / grid'};
f=plt_srcdlocinmap(dto2allbst/sps,dlocproj2allbst,[],'km','tarvl',...
  10,cstr,'s','log10','grid',[-4 4],[-4 4],0.1,0.1);
keyboard
%%
%plot the projected loc diff between each LFE source to all others
[f,mpjdist2all] = plt_srcdlocall(dlocproj2allbst,0.1,'km');
xlabel(f.ax(1),'Diff loc along min-rmse direc. (km)');
xlabel(f.ax(2),'Abs dist. along min-rmse direc. (km)');
xlabel(f.ax(3),'Diff loc along ort direc. (km)');
xlabel(f.ax(4),'Abs dist. along ort direc. (km)');
xlim(f.ax([1 3]),[-5 5]);
xlim(f.ax([2 4]),[0 5]);

%plot the projected loc diff between above source pairs
[f,mpjdist] = plt_srcdlocNtoNm(dlocprojbst,0.1,'km');
xlabel(f.ax(1),'Diff loc along min-rmse direc. (km)');
xlabel(f.ax(2),'Abs dist. along min-rmse direc. (km)');
xlabel(f.ax(3),'Diff loc along ort direc. (km)');
xlabel(f.ax(4),'Abs dist. along ort direc. (km)');
xlim(f.ax([1 3]),[-5 5]);
xlim(f.ax([2 4]),[0 5]);

keyboard

%% src amp ratio, data vs noise
%%%combine the direct deconvolved amp ratio between all station pairs of
%%%all burst wins, and summarize into one histogram
f3 = initfig(3*nsta,3,1,nsta); %plot histograms of source amp
orient(f3.fig,'landscape');
% tit=supertit(f3.ax,supertstr);
% movev(tit,0.2);
[f3.ax(1:end),mampr,madampr] = plt_deconpk_rat_comb4th(f3.ax(1:end),srcamprall,imp,'k','hist');
[f3.ax(1:end),mamprn,madamprn] = plt_deconpk_rat_comb4th(f3.ax(1:end),srcampralln,impn,'r','hist');
keyboard



%%
%%%the direct/scaled deconvolved pos/neg source peak ratio between all station pairs, for each
%%%burst win separately
nsrc = allbstsig.nsrc;
msrcampr = allbstsig.msrcampr;
madsrcampr = allbstsig.madsrcampr;
f1 = initfig(16,5,1,size(msrcampr,2));
f1 = plt_deconpk_rat4th(f1,msrcampr,madsrcampr,nsrc,'b');

nsrc = allbstnoi.nsrc;
msrcampr = allbstnoi.msrcampr;
madsrcampr = allbstnoi.madsrcampr;
f1 = plt_deconpk_rat4th(f1,msrcampr,madsrcampr,nsrc,'r');

%%
%%%deviation of source amp ratio from some median vs. RCC
lgdevsrcamprall = allbstsig.lgdevsrcamprall;
rccpairsrcall = allbstsig.rccpairsrcall;
rcccatsrcall = allbstsig.rcccatsrcall;
f3 = initfig(15,7,2,size(lgdevsrcamprall,2)); %initialize fig
f3 = plt_deconpk_ratdevvsrcc4th(f3,lgdevsrcamprall,rccpairsrcall,rcccatsrcall,'b');

lgdevsrcamprall = allbstnoi.lgdevsrcamprall;
rccpairsrcall = allbstnoi.rccpairsrcall;
rcccatsrcall = allbstnoi.rcccatsrcall;
f3 = plt_deconpk_ratdevvsrcc4th(f3,lgdevsrcamprall,rccpairsrcall,rcccatsrcall,'r');

%%
%%%diff between predicted arrival and selected peak at 4th sta
pred4offtrall = allbstsig.pred4offtrall;
f4 = initfig(5,5,1,size(pred4offtrall,2)); %initialize fig
f4 = plt_errorof4thtarvlpred(f4,pred4offtrall,offmax,'b');

pred4offtrall = allbstnoi.pred4offtrall;
f4 = plt_errorof4thtarvlpred(f4,pred4offtrall,offmax,'r');

%%
%%%preserved sources' amp ratio between 4th and 1st stas
srcamprall = allbstsig.srcamprall;
impindep4thall = allbstsig.impindep4thall;
f5 = initfig(12,5,1,3); %initialize fig
f5 = plt_deconpk_rat14(f5,impindep4thall,srcamprall,'b');

srcamprall = allbstnoi.srcamprall;
impindep4thall = allbstnoi.impindep4thall;
f5 = plt_deconpk_rat14(f5,impindep4thall,srcamprall,'r');

%%
%%%cloest waveform peak ratio between sta pairs VS. decon src amp ratio between same pairs
clppkhtwfall = allbstsig.clppkhtwfall;
psrcampsall = allbstsig.psrcampsall;
clnpkhtwfall = allbstsig.clnpkhtwfall;
nsrcampsall = allbstsig.nsrcampsall;
f6 = initfig(15,7,2,size(clppkhtwfall,2)); %initialize fig
f6 = plt_deconpkratvswfpkrat4th(f6,clppkhtwfall,psrcampsall,clnpkhtwfall,nsrcampsall);
% f6 = initfig(15,7,2,4); %initialize fig
% clppkhtwfall = [clppkhtwfall(:,1:3) clppkhtwfall(:,end)];  %I only want KLNB
% clnpkhtwfall = [clnpkhtwfall(:,1:3) clnpkhtwfall(:,end)];  %I only want KLNB
% f6 = plt_deconpkratvswfpkrat4th(f6,clppkhtwfall,psrcampsall,clnpkhtwfall,nsrcampsall);

clppkhtwfall = allbstnoi.clppkhtwfall;
psrcampsall = allbstnoi.psrcampsall;
clnpkhtwfall = allbstnoi.clnpkhtwfall;
nsrcampsall = allbstnoi.nsrcampsall;
f7 = initfig(15,7,2,size(clppkhtwfall,2)); %initialize fig
f7 = plt_deconpkratvswfpkrat4th(f7,clppkhtwfall,psrcampsall,clnpkhtwfall,nsrcampsall);

%%
%%%cloest waveform peak ratio between sta pairs VS. decon src amp ratio between same pairs
clppkhtwfall = allbstsig.clppkhtwfall;
clnpkhtwfall = allbstsig.clnpkhtwfall;
clpkspanwfall = clppkhtwfall-clnpkhtwfall;
psrcampsall = allbstsig.psrcampsall;
nsrcampsall = allbstsig.nsrcampsall;
srcampsspanall = psrcampsall-nsrcampsall;
f7 = initfig(15,4,1,size(clppkhtwfall,2)); %initialize fig
f7 = plt_deconpkratvswfpkrat4th(f7,clpkspanwfall,srcampsspanall);


%% how many sources from noise compared to data, after checked at 4th stations
nsrcs = allbstsig.nsrc;  %number of sources from signal
nsrcn = allbstnoi.nsrc;  %number of sources from synthetic noise

frac = nsrcn./nsrcs; %fraction during day (noisy) times

f8 = initfig(5,4,1,1); %initialize fig
ax = f8.ax(1);
hold(ax,'on');
grid(ax,'on');
histogram(ax,frac,'BinWidth',0.05);
plot(ax,[median(frac) median(frac)],ax.YLim,'--','color','r','linew',1);
ylabel(ax,'Count of burst wins');
xlabel(ax,'Fraction');
% title(ax,'Burst wins with high-correlation between 14');
title(ax,'All burst wins');

keyboard

%% what is the propagation direction and how good is the estimate
propflag = 'map';
% propflag = 'spl';

checkflag = 0;  %whether using the result after 4th-sta check
% checkflag = 1;

if strcmp(propflag,'map')
  if checkflag == 0
    propang = allbstsig.propang;
    proppear = allbstsig.proppear;
  else
    propang = allbstsig.propang4th;
    proppear = allbstsig.proppear4th;
  end
  pearmin = 0.5;
elseif strcmp(propflag,'spl')
  if checkflag == 0
    propang = allbstsig.propspang;
    proppear = allbstsig.propsppear;
  else
    propang = allbstsig.propspang4th;
    proppear = allbstsig.propsppear4th;
  end
  pearmin = 0.5;
end

widin = 12; htin = 5;
nrow = 1; ncol = 3;
f9 = initfig(widin,htin,nrow,ncol);
xran = [0.1 0.95]; yran = [0.1 0.95];
xsep = 0.05; ysep = 0.05;
optaxpos(f9,nrow,ncol,xran,yran,xsep,ysep);

ax=f9.ax(1); 
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
scatter(ax,propang,proppear,20,'k','filled');
plot(ax,ax.XLim,[pearmin pearmin],'r--');
xlabel(ax,'Propagation direction');
ylabel(ax,'Weighted Pearson coeff');

ax=f9.ax(2);
polaraxes(ax);
theta = propang;
theta = deg2rad(theta);
binedge = deg2rad(2.5: 5: (270+92.5));
polarhistogram(theta,'binedges',binedge,'normalization','count',...
               'facec','k','facea',0.8); 
ax = gca; hold(ax,'on');
ax.RTickLabelRotation = 45;
ax.ThetaMinorTick = 'on';
ax.TickLength = [0.02 0];
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.FontSize = 10;
% ax.RTick = 0:1:5;
% ax.LineWidth = 1.5;
ax.Box = 'on';
ax.GridAlpha = 0.3;
text(ax,0,0.9*ax.RLim(2),'N','HorizontalAlignment',"center",'FontSize',14);
text(ax,deg2rad(90),0.9*ax.RLim(2),'E','HorizontalAlignment',"center",'FontSize',14);
text(ax,deg2rad(180),0.9*ax.RLim(2),'S','HorizontalAlignment',"center",'FontSize',14);
text(ax,deg2rad(270),0.9*ax.RLim(2),'W','HorizontalAlignment',"center",'FontSize',14);


ax=f9.ax(3);
polaraxes(ax);
propang = propang(proppear>=pearmin);
proppear = proppear(proppear>=pearmin);
theta = propang;
theta = deg2rad(theta);
binedge = deg2rad(2.5: 5: 92.5);
polarhistogram(theta(propang>0 & propang<=90),'binedges',binedge,'normalization','count',...
               'facec',[0 1 1],'facea',0.8);
ax = gca; hold(ax,'on');
binedge = deg2rad(2.5: 5: 92.5 + 90);           
polarhistogram(theta(propang>90 & propang<=180),'binedges',binedge,'normalization','count',...
               'facec','k','facea',0.8);
binedge = deg2rad(2.5: 5: 92.5 + 180);           
polarhistogram(theta(propang>180 & propang<=270),'binedges',binedge,'normalization','count',...
               'facec',[1 96/255 0],'facea',0.8);
binedge = deg2rad(2.5: 5: 92.5 + 270);           
polarhistogram(theta(propang>270 & propang<=360),'binedges',binedge,'normalization','count',...
               'facec','k','facea',0.2);
ax.RTickLabelRotation = 45;
ax.ThetaMinorTick = 'on';
ax.TickLength = [0.02 0];
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.FontSize = 10;
% ax.RTick = 0:1:5;
% ax.LineWidth = 1.5;
ax.Box = 'on';
ax.GridAlpha = 0.3;
text(ax,deg2rad(60),0.8*ax.RLim(2),num2str(sum(propang>0 & propang<=90)),'HorizontalAlignment','center',...
     'FontSize',11);           
text(ax,deg2rad(130),0.8*ax.RLim(2),num2str(sum(propang>90 & propang<=180)),'HorizontalAlignment','center',...
     'FontSize',11);           
text(ax,deg2rad(225),0.8*ax.RLim(2),num2str(sum(propang>180 & propang<=270)),'HorizontalAlignment','center',...
     'FontSize',11);
text(ax,deg2rad(315),0.8*ax.RLim(2),num2str(sum(propang>270 & propang<=360)),'HorizontalAlignment','center',...
     'FontSize',11); 

%%
%%%the separation in arrival time vs. distance between sources N--N-1, N--N-2, and N--N-3,
%%%for data and noise
widin = 7;  % maximum width allowed is 8.5 inches
htin = 3.5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f1 = initfig(widin,htin,nrow,ncol); %initialize fig
xran = [0.1 0.94]; yran = [0.12 0.96];
xsep = 0.08; ysep = 0.05;
optaxpos(f1,nrow,ncol,xran,yran,xsep,ysep);

dtarvlnn1all = allbstsig.dtarvlnn1all;
distarvlnn1all = allbstsig.distarvlnn1all;
f1 = plt_absdistvsdt(f1,[dtarvlnn1all,distarvlnn1all],1,sps,'tarvl','b');
% median(distarvlnn1all(dtarvlnn1all/sps<=1))
% median(distarvlnn1all(dtarvlnn1all/sps<=0.5))
% median(distarvlnn1all)
% median(dtarvlnn1all/sps)
% legend(f1.ax(3),'Data','location','east');

dtarvlnn1all = allbstnoi.dtarvlnn1all;
distarvlnn1all = allbstnoi.distarvlnn1all;
f1 = plt_absdistvsdt(f1,[dtarvlnn1all,distarvlnn1all],1,sps,'tarvl','r');
ylim(f1.ax(1),[0 0.6]);
ylim(f1.ax(2),[0 3]);


%% in sample sapce, 2ndary sources removed
% % allbstsig = allbstnoi;
% dt2all = allbstsig.dt2allbst;
% dist2allsp = allbstsig.dist2allspbst;
% dloc2allsp = allbstsig.dloc2allspbst;
% imp = allbstsig.impindepall;
% f=plt_srcdlocinmap(dt2all,dloc2allsp,sps,'spl');
% 
% nsrc = allbstsig.nsrcraw;
% dto2allbst=[];  %diff. origin time between event pairs
% dloco2allbst=[];  
% disto2allbst=[];
% impstall = [];
% tcorall = [];
% for i = 1: size(trange,1)
%   ist = sum(nsrc(1:i-1))+1;
%   ied = ist+nsrc(i)-1;
%   impi = imp(ist:ied,:);
%   if ~isempty(impi)
%     %convert time offset to relative loc
%     ftrans = 'interpchao';
%     [imploc, indinput] = off2space002(impi(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
%     [imploc0, ~] = off2space002([0 0],sps,ftrans,0);  % a ref source at 0,0
%     tcor = round((imploc(:,6)-imploc0(6))*sps);   % travel time difference between each source and a ref source at 0,0
%     torispl = impi(:,1)-tcor;
%     [torisplst, indsort] = sortrows(torispl,1); %note, sort by origin time now
%     implocst = imploc(indsort, :);
%     impst = impi(indsort, :);
%     %I suspect the diff loc and distance between pairs should be the same as using tarvl, but correspondance with time
%     %should be different
%     [dto2all,dloco2all,disto2all] = srcdistall(torisplst,impst(:,7:8),[0 tmax2all*sps]);
%     
%   else
%     dto2all = [];
%     dloco2all = [];
%     disto2all = [];
%     impst = [];
%     dto2all = [];
%     indsort = [];
%     tcor = [];
%   end
%   dto2allbst = [dto2allbst; dto2all];
%   dloco2allbst = [dloco2allbst; dloco2all] ;
%   disto2allbst = [disto2allbst; disto2all];
%   impstall = [impstall; impst];
%   tcorall = [tcorall; tcor(indsort)];
% end
% f=plt_srcdlocinmap(dto2allbst,dloco2allbst,sps,'spl','tori');
% 
% %plot the diff time and euclidean distance between each LFE source to all others
% tlensec = 14981;  %total length of time in secs for all bursts combined
% f = plt_srcdistall(dt2all,dist2allsp,sps,40/sps,tlensec,1,'spl');
% f = plt_srcdistall(dto2allbst,disto2allbst,sps,40/sps,tlensec,1,'spl','tori');
% 
% %plot the loc diff between each LFE source to all others
% f = plt_srcdlocall(dloc2allsp,1,'spl');
% 
% doffset{1,1} = allbstsig.distarvlspnn1all(:,2:3);
% doffset{2,1} = allbstsig.distarvlspnn2all(:,2:3);
% doffset{3,1} = allbstsig.distarvlspnn3all(:,2:3);
% doffset{4,1} = allbstsig.distarvlspnn4all(:,2:3);
% doffset{5,1} = allbstsig.distarvlspnn5all(:,2:3);
% eucdistsp{1,1} = allbstsig.distarvlspnn1all(:,1);
% eucdistsp{2,1} = allbstsig.distarvlspnn2all(:,1);
% eucdistsp{3,1} = allbstsig.distarvlspnn3all(:,1);
% eucdistsp{4,1} = allbstsig.distarvlspnn4all(:,1);
% eucdistsp{5,1} = allbstsig.distarvlspnn5all(:,1);
% dtarvl{1,1} = allbstsig.dtarvlnn1all;
% dtarvl{2,1} = allbstsig.dtarvlnn2all;
% dtarvl{3,1} = allbstsig.dtarvlnn3all;
% dtarvl{4,1} = allbstsig.dtarvlnn4all;
% dtarvl{5,1} = allbstsig.dtarvlnn5all;
% %plot the diff time and euclidean distance between above source pairs
% f = plt_srcdistNtoNm(dtarvl,eucdistsp,sps,40/sps,tlensec,1,'spl');
% %plot the loc diff between above source pairs
% f = plt_srcdlocNtoNm(doffset,1,'spl');
% 
% %%%if you want to copy several figs to a new summary to put everything into 
% % f1 = initfig(10,9,2,3);
% % for i = 1:3
% %   axchild=f.ax(i).Children;
% %   copyobj(axchild, f1.ax(i));
% % end
