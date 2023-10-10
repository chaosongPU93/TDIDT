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
% 
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
  pltflag = 0;  %do not create summary plots for each choice of inputs
  % rccmwsec = 0.25; %use 0.5s or 0.25s 
  rccmwsec = 0.5; %use 0.5s or 0.25s
  
  %%all bursts using real data
  noiseflag = 0;
  allbstsig = deconv_ref_4s_exp_4thsta_fn(1:size(trange,1),normflag,noiseflag,pltflag,rccmwsec); %
  savefile = 'deconv_stats4th_allbstsig.mat';
  save(strcat(rstpath, '/MAPS/',savefile), 'allbstsig');
  
  %%%all bursts using synthetic noise
  noiseflag = 1;
  allbstnoi = deconv_ref_4s_exp_4thsta_fn(1:size(trange,1),normflag,noiseflag,pltflag,rccmwsec);  %
  savefile = 'deconv_stats4th_allbstnoi.mat';
  save(strcat(rstpath, '/MAPS/',savefile), 'allbstnoi');
else
  savefile = 'deconv_stats4th_allbstsig.mat';
  load(strcat(rstpath, '/MAPS/',savefile));
  savefile = 'deconv_stats4th_allbstnoi.mat';
  load(strcat(rstpath, '/MAPS/',savefile));
end

% keyboard

% bstsig = allbstsig;
% bstnoi = allbstnoi;

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
 
%% a specific plot for Allan's presentation, similar to last fig of AGU2022 poster
% %%%param for secondary sources removed
% locxyprojall = allbstsig.locxyprojall;
% tarvlsplstall = allbstsig.impindepall(:,1);
% nsrc = allbstsig.nsrcraw;
% dtarvlnn1 = allbstsig.dtarvlnn1all;
% dtarvlnn2 = allbstsig.dtarvlnn2all;
% dtarvlnn3 = allbstsig.dtarvlnn3all;
% dtarvlnn4 = allbstsig.dtarvlnn4all;
% dtarvlnn5 = allbstsig.dtarvlnn5all;
% dtarvlnn6 = allbstsig.dtarvlnn6all;
% dtarvlnn7 = allbstsig.dtarvlnn7all;
% dtarvlnn8 = allbstsig.dtarvlnn8all;
% imp = allbstsig.impindepall;
% locxyprojalln = allbstnoi.locxyprojall;
% tarvlsplstalln = allbstnoi.impindepall(:,1);
% nsrcn = allbstnoi.nsrcraw;
% dtarvlnn1n = allbstnoi.dtarvlnn1all;
% dtarvlnn2n = allbstnoi.dtarvlnn2all;
% dtarvlnn3n = allbstnoi.dtarvlnn3all;
% dtarvlnn4n = allbstnoi.dtarvlnn4all;
% dtarvlnn5n = allbstnoi.dtarvlnn5all;
% dtarvlnn6n = allbstnoi.dtarvlnn6all;
% dtarvlnn7n = allbstnoi.dtarvlnn7all;
% dtarvlnn8n = allbstnoi.dtarvlnn8all;
% impn = allbstnoi.impindepall;
% supertstr = 'Secondary sources removed';
% fnsuffix = [];

%%%param for further checked at KLNB
locxyprojall = allbstsig.locxyproj4thall;
tarvlsplstall = allbstsig.impindep4thall(:,1);
nsrc = allbstsig.nsrc;
dtarvlnn1 = allbstsig.dtarvlnn14thall;
dtarvlnn2 = allbstsig.dtarvlnn24thall;
dtarvlnn3 = allbstsig.dtarvlnn34thall;
dtarvlnn4 = allbstsig.dtarvlnn44thall;
dtarvlnn5 = allbstsig.dtarvlnn54thall;
dtarvlnn6 = allbstsig.dtarvlnn64thall;
dtarvlnn7 = allbstsig.dtarvlnn74thall;
dtarvlnn8 = allbstsig.dtarvlnn84thall;
imp = allbstsig.impindep4thall;
locxyprojalln = allbstnoi.locxyproj4thall;
tarvlsplstalln = allbstnoi.impindep4thall(:,1);
nsrcn = allbstnoi.nsrc;
dtarvlnn1n = allbstnoi.dtarvlnn14thall;
dtarvlnn2n = allbstnoi.dtarvlnn24thall;
dtarvlnn3n = allbstnoi.dtarvlnn34thall;
dtarvlnn4n = allbstnoi.dtarvlnn44thall;
dtarvlnn5n = allbstnoi.dtarvlnn54thall;
dtarvlnn6n = allbstnoi.dtarvlnn64thall;
dtarvlnn7n = allbstnoi.dtarvlnn74thall;
dtarvlnn8n = allbstnoi.dtarvlnn84thall;
impn = allbstnoi.impindep4thall;
supertstr = 'Further checked at KLNB';
fnsuffix = '4th';

% typepltnoi = 1; %plot noise
typepltnoi = 2; %plot data -noise

dtarvlplt = dtarvlnn1;
dtarvlpltn = dtarvlnn1n;
nsep = 1;
% dtarvlplt = dtarvlnn2;
% dtarvlpltn = dtarvlnn2n;
% nsep = 2;
% dtarvlplt = dtarvlnn3;
% dtarvlpltn = dtarvlnn3n;
% nsep = 3;

% %% refer to NOTE 20230707 
% %%%For quantifying noise, one possible way is to obtain the average of 2 
% %%%inter-event time every detection (prior and posterior), VS its amp, 
% %%%for real data
% dtmean = [];
% amp = [];
% dtmeann = [];
% ampn = [];
% for i = 1: size(trange,1)
%   ist = sum(nsrc(1:i-1))+1;
%   ied = ist+nsrc(i)-1;
%   impi = imp(ist:ied,:);
%   if nsrc(i) < 3
%     continue
%   end
%   dtbfi = diff(impi(1:end-1,1)); 
%   dtafi = diff(impi(2:end,1));
%   dtmeani = mean([dtbfi dtafi], 2);
%   ampi = mean(impi(2:end-1,[2 4 6]),2);
%   dtmean = [dtmean; dtmeani];
%   amp = [amp; ampi];

%   ist = sum(nsrcn(1:i-1))+1;
%   ied = ist+nsrcn(i)-1;
%   impi = impn(ist:ied,:);
%   if nsrcn(i) < 3
%     continue
%   end
%   dtbfi = diff(impi(1:end-1,1)); 
%   dtafi = diff(impi(2:end,1));
%   dtmeani = mean([dtbfi dtafi], 2);
%   ampi = mean(impi(2:end-1,[2 4 6]),2);
%   dtmeann = [dtmeann; dtmeani];
%   ampn = [ampn; ampi];
% end
% nbin = 25;
% [dtbin,indbin,n] = binxeqnum(dtmean,nbin);
% [dtbinn,indbinn,nn] = binxeqnum(dtmeann,nbin);
% for i = 1: nbin
%   ampbini = amp(indbin{i});
%   ampbincnt(i) = median(ampbini);
%   dtbincnt(i) = median(dtbin{i});
%   ampbinin = amp(indbinn{i});
%   ampbincntn(i) = median(ampbinin);
%   dtbincntn(i) = median(dtbinn{i});
% end
% widin = 6;  % maximum width allowed is 8.5 inches
% htin = 9;   % maximum height allowed is 11 inches
% nrow = 3;
% ncol = 1;
% f = initfig(widin,htin,nrow,ncol); %initialize fig
% ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% scatter(ax,dtmean/sps,amp,20);
% plot(ax,[0.25 0.25],[0 10],'r--','linew',2);
% scatter(ax,dtbincnt/sps,ampbincnt,20,'k','filled');
% xlabel(ax,'mean tsep (s)');
% ylabel(ax,'mean amp');
% xlim(ax,[0 4]);
% ylim(ax,[0 4]);
% title(ax,'Data');

% ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% scatter(ax,dtmeann/sps,ampn,20);
% plot(ax,[0.25 0.25],[0 10],'r--','linew',2);
% scatter(ax,dtbincntn/sps,ampbincntn,20,'k','filled');
% xlabel(ax,'mean tsep (s)');
% ylabel(ax,'mean amp');
% xlim(ax,[0 4]);
% ylim(ax,[0 4]);
% title(ax,'Synthetic noise');

% ax=f.ax(3); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% plot(ax,[0.25 0.25],[0 10],'k--','linew',1);
% scatter(ax,dtbincnt/sps,ampbincnt,20,'b','filled');
% scatter(ax,dtbincntn/sps,ampbincntn,20,'r','filled');
% xlabel(ax,'mean tsep (s)');
% ylabel(ax,'mean amp');
% xlim(ax,[0 4]);
% ylim(ax,[0 1]);
% keyboard


% %% fraction of event pairs w/i a diff time cut, and fraction of all catalog
% % nsep = 14;
% for nsep = 1:15
% 
% dcut = 0.25*nsep+0.125;
% nsrcsep = nsrc-nsep;
% nsrcsep(nsrcsep<0) = 0; 
% ndcutpair = zeros(size(trange,1), 1);
% ndcutsrc = zeros(size(trange,1), 1);
% ndcutsrc2 = zeros(size(trange,1), 1);
% imppair = cell(size(trange,1), 1);
% imppairuni = cell(size(trange,1), 1);
% impcont = cell(size(trange,1), 1);
% impcontuni = cell(size(trange,1), 1);
% indcont = cell(size(trange,1), 1);
% indcontuni = cell(size(trange,1), 1);
% m = 15;
% 
% for i = 1: size(trange,1)
%   if nsrc(i) == 0
%     continue
%   end
%   ist = sum(nsrc(1:i-1))+1;
%   ied = ist+nsrc(i)-1;
%   impi = imp(ist:ied,:);
%   %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
%   dtarvl = srcdistNtoNm(impi(:,1),impi(:,7:8),m);
% %   dtarvlnnsepall = [dtarvlnnsepall; dtarvl{nsep}];
%   if isempty(dtarvl{nsep})
%     continue
%   end
%   impbf = impi(1:end-nsep,:);
%   impaf = impi(1+nsep:end,:);
%   if ~isequal(size(impbf,1),length(dtarvl{nsep}))
%     disp('Check');
%   end
%   ind = find(dtarvl{nsep}/sps <= dcut);
%   ndcutpair(i,1) = length(ind);
%   tmp = [];
%   for j = 1: length(ind)
%     tmp = [tmp; impbf(ind(j),:); impaf(ind(j),:)];
%   end
%   imppair{i,1} = tmp;
%   % impsort = sortrows(imppair,1);
%   imppairuni{i,1} = unique(tmp,'rows','stable');
%   % impall = [impall; impuni];
%   ndcutsrc(i,1) = size(imppairuni{i,1},1);
% 
%   tmp = [];
%   for j = 1: length(ind)
%     tmp = [tmp; impi(ind(j):ind(j)+nsep,:)];
%   end
%   impcont{i,1} = tmp;
%   impcontuni{i,1} = unique(tmp,'rows','stable');
%   % impall = [impall; impuni];
%   ndcutsrc2(i,1) = size(impcontuni{i,1},1);
% 
%   tmp = [];
%   for j = 1: length(ind)
%     tmp = [tmp; reshape(ind(j):ind(j)+nsep, [], 1)];
%   end
%   indcont{i,1} = tmp;
%   indcontuni{i,1} = unique(tmp,'rows','stable');
% end
% 
% imppairunia = cat(1,imppairuni{:});
% ndcutsrcall = sum(ndcutsrc);
% fracsrcall = ndcutsrcall/sum(nsrc)*100;
% 
% impcontunia = cat(1,impcontuni{:});
% ndcutsrc2all(nsep,1) = sum(ndcutsrc2);
% fracsrc2all(nsep,1) = ndcutsrc2all(nsep,1)/sum(nsrc)*100;
% 
% ndcutpairall(nsep,1) = sum(ndcutpair);
% fracpairall = ndcutpairall(nsep,1)/sum(nsrcsep)*100;
% 
% end
% 
% for nsep = 1:15
%   dcut = 0.25*nsep+0.125;
%   fprintf('%d clusters of %d consecutive events (%d/%d unique) w/i %.3f s \n',...
%     ndcutpairall(nsep,1), nsep+1, ndcutsrc2all(nsep,1), sum(nsrc), dcut);
% end
% 
% fracsrc2all = [100; fracsrc2all];
% fprintf('%.3f \n',fracsrc2all);
% dfracsrc2all = fracsrc2all(1:end-1) - fracsrc2all(2:end);
% fprintf('%.3f \n',dfracsrc2all);
% 
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
    [~,dlocproj2all] = srcdistall(tarvlsplst,locxyproj,[0 2*sps]);
    dlocproj2allbst = [dlocproj2allbst; dlocproj2all];
    [~,dlocproj] = srcdistNtoNm(tarvlsplst,locxyproj,2);
    dlocprojnn1bst = [dlocprojnn1bst; dlocproj{1}];
    dlocprojnn2bst = [dlocprojnn2bst; dlocproj{2}];
  end
end
dist2all = abs(dlocproj2allbst(:,1));
[bincnt,binhgt,count,normalizer] = histbinbyarea(dist2all,binedge,'countdensity');
% binhgt = binhgt./length(dist2all);
binhgt = count./length(dist2all);
p1=bar(ax,bincnt,binhgt,1,'stacked','b','facea',0.6);
% stairs(ax,binedge,[binhgt; binhgt(end)],'k','LineWidth',1);
plot(ax,[median(dist2all) median(dist2all)],ax.YLim,'b--','linew',1.5);
text(ax,median(dist2all)+0.1,ax.YLim(2)-0.1*range(ax.YLim),...
  sprintf('%.2f',median(dist2all)),'HorizontalAlignment','left');
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
    [~,dlocproj2all] = srcdistall(tarvlsplst,locxyproj,[0 2*sps]);
    dlocproj2allbstn = [dlocproj2allbstn; dlocproj2all];
  end
end
dist2alln = abs(dlocproj2allbstn(:,1));
[bincnt,binhgt,count,normalizer] = histbinbyarea(dist2alln,binedge,'countdensity');
% binhgt = binhgt./length(dist2all);
binhgt = count./length(dist2all);
p2=bar(ax,bincnt,binhgt,1,'stacked','r','facea',0.6);
% stairs(ax,binedge,[binhgt; binhgt(end)],'k','LineWidth',1);
plot(ax,[median(dist2alln) median(dist2alln)],ax.YLim,'r--','linew',1.5);
text(ax,median(dist2alln)+0.1,ax.YLim(2)-0.2*range(ax.YLim),...
  sprintf('%.2f',median(dist2alln)),'HorizontalAlignment','left');
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
xlabel(ax,sprintf('Diff. arrival time between sources N and N-%d (s)',nsep));
legend(ax,[p1,p2],'Data','Synthetic noise','location','east');
xlim(ax,[0 1.5]);
% ylim(ax,[0 3.5]);
longticks(ax,2);
hold(ax,'off');
tit=supertit(f.ax,supertstr);
movev(tit,0.2);

orient(f.fig,'landscape');
fname = strcat(sprintf('nn%ddiff',nsep),fnsuffix,'.pdf');
% print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));

% keyboard

% %% CDF 
% [cdfval,x] = ecdf(dtarvlnn1/sps);
% figure;
% plot(x,cdfval,'r','linew',1); hold on ;
% aa = sum(dtarvlnn1/sps<=0.375)/length(dtarvlnn1)
% [cdfval,x] = ecdf(dtarvlnn2/sps);
% plot(x,cdfval,'b','linew',1); hold on ;
% bb = sum(dtarvlnn2/sps<=0.625)/length(dtarvlnn2)

%%
%%%summarize the whole catalog, diff arrival time and fractions
f = initfig(5.5,4,1,1); %initialize fig
tit=supertit(f.ax,supertstr);
movev(tit,0.2);
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
binwdt = 0.05;
dtcut = 0.25*nsep+0.125;
xran = [0 2];
nx = round(xran(2)/binwdt)+1;
edges = xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2;
cnt = xran(1): binwdt: xran(2);
N=histcounts(dtarvlplt/sps,edges,'normalization','count');
Nn = N./length(dtarvlplt);
frac = sum(dtarvlplt/sps<=dtcut)/length(dtarvlplt);
p(1)=plot(ax,cnt,Nn,'color','b','LineWidth',1);
label{1}='Data';
N=histcounts(dtarvlpltn/sps,edges,'normalization','count');
Nnn = N./length(dtarvlplt);
fracn = sum(dtarvlpltn/sps<=dtcut)/length(dtarvlpltn);
p(2)=plot(ax,cnt,Nnn,'color','r','LineWidth',1);
label{2}='Synthetic noise';
fracdif = (sum(dtarvlplt/sps<=dtcut)-sum(dtarvlpltn/sps<=dtcut)) / ...
  (length(dtarvlplt)-length(dtarvlpltn));
p(3)=plot(ax,cnt,Nn-Nnn,'color','k','LineWidth',1.5);
label{3}='Data - Synthetic noise';
text(ax,0.7,0.7,sprintf('%.2f',frac),'Color','b','HorizontalAlignment','left','Units','normalized');
text(ax,0.7,0.63,sprintf('%.2f',fracn),'Color','r','HorizontalAlignment','left','Units','normalized');
text(ax,0.7,0.56,sprintf('%.2f',fracdif),'Color','k','HorizontalAlignment','left','Units','normalized');
legend(ax,p,label);
ylabel(ax,'Normalized count');
xlabel(ax,sprintf('Diff. arrival between sources N and N-%d (s)',nsep));
xlim(ax,xran);
if nsep == 1
  yran = [0 0.2];
elseif nsep == 2
  yran = [0 0.06];
end
ylim(ax,yran);
longticks(ax,2);
hold(ax,'off');
% keyboard

%%
f = initfig(12,4.5,1,2); %initialize fig
tit=supertit(f.ax,supertstr);
movev(tit,0.3);
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
plot(ax,cnt,Nn,'color','b','LineWidth',1);
text(ax,0.95,0.85,sprintf('Fraction w/i %.3f s: %.2f',dtcut,frac),'HorizontalAlignment','right',...
  'Units','normalized','FontSize',12);
ylabel(ax,'Normalized count');
xlabel(ax,sprintf('Diff. arrival between sources N and N-%d (s)',nsep));
xlim(ax,xran);
ylim(ax,yran);
longticks(ax,2);
title(ax,'Data');
hold(ax,'off');

ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
plot(ax,cnt,Nnn,'color','b','LineWidth',1);
text(ax,0.95,0.85,sprintf('Fraction w/i %.3f s: %.2f',dtcut,fracn),'HorizontalAlignment','right',...
  'Units','normalized','FontSize',12);
ylabel(ax,'Normalized count');
xlabel(ax,sprintf('Diff. arrival between sources N and N-%d (s)',nsep));
xlim(ax,xran);
ylim(ax,yran);
longticks(ax,2);
title(ax,'Synthetic noise');
hold(ax,'off');

% keyboard



%% 
%%%first bin source by amp, then plot diff arrival time for N and N-1 for each amp bin
ampafdt = [];
ampbfdt = [];

for i = 1: size(trange,1)
  ist = sum(nsrc(1:i-1))+1;
  ied = ist+nsrc(i)-1;
  impi = imp(ist:ied,:);
  impaf = impi(1+nsep:end,:);
  ampaf = mean(impaf(:,[2 4 6]),2);
  ampafdt = [ampafdt; ampaf];
  impbf = impi(1:end-nsep,:);
  ampbf = mean(impbf(:,[2 4 6]),2);
  ampbfdt = [ampbfdt; ampbf];
end
ampdt = mean([ampbfdt ampafdt],2);
minnum = 500;

ampafdtn = [];
ampbfdtn = [];
for i = 1: size(trange,1)
  ist = sum(nsrcn(1:i-1))+1;
  ied = ist+nsrcn(i)-1;
  impi = impn(ist:ied,:);
  impaf = impi(1+nsep:end,:);
  ampaf = mean(impaf(:,[2 4 6]),2);
  ampafdtn = [ampafdtn; ampaf];
  impbf = impi(1:end-nsep,:);
  ampbf = mean(impbf(:,[2 4 6]),2);
  ampbfdtn = [ampbfdtn; ampbf];
end
ampdtn = mean([ampbfdtn ampafdtn],2);
minnumn = 250;

ampplt = ampdt;
amppltn = ampdtn;

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

%amp distribution of source pair 
figure
subplot(121)
h=histogram(log10(ampplt),'BinWidth',0.25); hold on;
ax=gca;
plot(ax,ax.XLim,[minnum minnum],'k--');
xlim([-1.5 1]);
title('Mean of the pair for data');
ylabel(ax,'Count');
xlabel(ax,sprintf('log_{10}{amp} between sources N and N-%d (s)',nsep));

subplot(122)
h=histogram(log10(amppltn),'BinWidth',0.25); hold on;
ax=gca;
plot(ax,ax.XLim,[minnumn minnumn],'k--');
xlim([-1.5 1]);
title('Mean of the pair for syn noise');
ylabel(ax,'Count');
xlabel(ax,sprintf('log_{10}{amp} between sources N and N-%d (s)',nsep));


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

%%
f = initfig(12,8,2,2); %initialize fig
tit=supertit(f.ax,supertstr);
movev(tit,0.2);
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% %%%%%%%%%% if bin by amp with a equal width
% nbin = h.NumBins;
% color = jet(nbin-1);
% binwdt = 0.1;
% iplt = 0;
% for i = 1: nbin-1 
%   ind = find(log10(ampplt)>=h.BinEdges(i) & log10(ampplt)<h.BinEdges(i+1));
%   dtarvli = dtarvlplt(ind);
%   if length(dtarvli)>=minnum
%     iplt = iplt+1;
%     dtarvlimed(iplt) = median(dtarvli)/sps;
%     ampbincnt(iplt) = (h.BinEdges(i)+h.BinEdges(i+1))/2;
%     [N,edges]=histcounts(dtarvli/sps,'binwidth',binwdt,'normalization','count');
%     edges = edges-binwdt/2;
%   %   edges = -binwdt/2: binwdt: 2+binwdt/2;
%     [N,edges]=histcounts(dtarvli/sps,edges,'normalization','count');
%     Nn = N./length(dtarvli);
% %     p(iplt)=stairs(ax,edges,[Nn Nn(end)],'color',color(iplt,:),'LineWidth',1);
%     cnt = (edges(1:end-1)+edges(2:end))/2;
%     p(iplt)=plot(ax,cnt,Nn,'color',color(iplt,:),'LineWidth',1);
%     label{iplt} = sprintf('%.2f',ampbincnt(iplt));
%   end  
% end
% %%%%%%%%%% if bin by amp with a equal width

%%%%%%%%%% if bin by amp with a equal number
nbin = 5;
[ampbin,indbin,n] = binxeqnum(ampplt,nbin);
color = jet(nbin);
% color = gray(nbin+1);
% color = flipud(color(1:end-1,:));
% color = flipud(kelicmap(nbin));
binwdt = 0.05;
% dtcut1 = 0.5*nsep*sps;
% dtcut2 = dtcut1+0.125*sps;
dtcut = 0.25*nsep+0.125;
if nsep == 3
  xran = [0 5];
else
  xran = [0 2];
end
nx = round(xran(2)/binwdt)+1;
Nn = zeros(nx, nbin);
edges = xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2;
cnt = xran(1): binwdt: xran(2);
for i = 1: nbin
% for i = 1: 1
  ind = indbin{i};
  dtarvli = dtarvlplt(ind);
  % dtarvlimed1(i) = median(dtarvli(dtarvli<=dtcut1))/sps;
  % dtarvlimed2(i) = median(dtarvli(dtarvli<=dtcut2))/sps;
  dtarvlimed(i) = median(dtarvli)/sps;
  dtarvlimode(i) = mode(dtarvli)/sps;
  ampbincnt(i) = median(ampbin{i});
  % [N,edges]=histcounts(dtarvli/sps,'binwidth',binwdt,'normalization','count');
  % edges = edges-binwdt/2;
  N=histcounts(dtarvli/sps,edges,'normalization','count');
  Nn(:,i) = N./mode(n);
  frac(i) = sum(dtarvli/sps<=dtcut)/mode(n);
  % cnt = (edges(1:end-1)+edges(2:end))/2;
  p(i)=plot(ax,cnt,Nn(:,i),'color',color(i,:),'LineWidth',1);
  label{i} = sprintf('amp of %.1f',ampbincnt(i));
  % label{i} = sprintf('%d/%dth amp',i,nbin);
  % keyboard
end
%%%%%%%%%% if bin by amp with a equal number
Nnm = mean(Nn,2);
% Nnm = median(Nn,2);
% for i = 1: nbin
%   p(i)=plot(ax,cnt,Nn(:,i)-Nnm,'color',color(i,:),'LineWidth',1);
%   label{i} = sprintf('amp of %.1f - mean',ampbincnt(i));
% end
p(nbin+1)=plot(ax,cnt,Nnm,'k-','LineWidth',1.5);
label{nbin+1} = sprintf('mean');  %median
legend(ax,p,label);
ylabel(ax,'Normalized count');
xlabel(ax,sprintf('Diff. arrival between sources N and N-%d (s)',nsep));
xlim(ax,xran);
if nsep == 1
  yran = [0 0.2];
elseif nsep == 2
  yran = [0 0.06];
else
  yran = ax.YLim;
end
ylim(ax,yran);
longticks(ax,2);
hold(ax,'off');
title(ax,'Data');
% keyboard

ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% %%%%%%%%%% if bin by amp with a equal width
% iplt = 0;
% for i = 1: h.NumBins-1 
%   ind = find(log10(amppltn)>=h.BinEdges(i) & log10(amppltn)<h.BinEdges(i+1));
%   dtarvlin = dtarvlpltn(ind);  
%   if length(dtarvlin)>=minnumn
%     iplt = iplt+1;
%     dtarvlinmed(iplt) = median(dtarvlin)/sps;
%     ampbinncnt(iplt) = (h.BinEdges(i)+h.BinEdges(i+1))/2;
%     [N,edges]=histcounts(dtarvlin/sps,'binwidth',binwdt,'normalization','count');
%     edges = edges-binwdt/2;
%   %   edges = -binwdt/2: binwdt: 2+binwdt/2;
%     [N,edges]=histcounts(dtarvlin/sps,edges,'normalization','count');
%     Nn = N./length(dtarvlin);
% %     p(iplt)=stairs(ax,edges,[Nn Nn(end)],'color',color(iplt,:),'LineWidth',1);
%     cnt = (edges(1:end-1)+edges(2:end))/2;
%     p(iplt)=plot(ax,cnt,Nn,'color',color(iplt,:),'LineWidth',1);
%     label{iplt} = sprintf('%.2f',ampbinncnt(iplt));
%   end  
% end
% %%%%%%%%%% if bin by amp with a equal width

%%%%%%%%%% if bin by amp with a equal number
[ampbinn,indbinn,nn] = binxeqnum(amppltn,nbin);
Nnn = zeros(nx, nbin);
for i = 1: nbin
% for i = 1: 1
  ind = indbinn{i};
  dtarvlin = dtarvlpltn(ind);
  % dtarvlinmed1(i) = median(dtarvlin(dtarvlin<=dtcut1))/sps;
  % dtarvlinmed2(i) = median(dtarvlin(dtarvlin<=dtcut2))/sps;
  dtarvlinmed(i) = median(dtarvlin)/sps;
  dtarvlinmode(i) = mode(dtarvlin)/sps;
  ampbinncnt(i) = median(ampbinn{i});
  % [N,edges]=histcounts(dtarvlin/sps,'binwidth',binwdt,'normalization','count');
  % edges = edges-binwdt/2;
  [N,edges]=histcounts(dtarvlin/sps,edges,'normalization','count');
  Nnn(:,i) = N./mode(n);
  % cnt = (edges(1:end-1)+edges(2:end))/2;
  if typepltnoi == 1 
    p(i)=plot(ax,cnt,Nnn(:,i),'color',color(i,:),'LineWidth',1);
    fracn(i) = sum(dtarvlin/sps<=dtcut)/mode(nn);
  elseif typepltnoi == 2 
    p(i)=plot(ax,cnt,Nn(:,i)-Nnn(:,i),'color',color(i,:),'LineWidth',1);
    fracn(i) = (frac(i)*mode(n)-sum(dtarvlin/sps<=dtcut)) / (mode(n)-mode(nn));  
    % fracn(i) = (frac(i)*mode(n)-sum(dtarvlin/sps<=dtcut)) / mode(n);  
  end
  label{i} = sprintf('amp of %.1f',ampbinncnt(i));
end
%%%%%%%%%% if bin by amp with a equal number
Nnnm = mean(Nnn,2);
% Nnm = median(Nn,2);
% for i = 1: nbin
%   p(i)=plot(ax,cnt,Nnn(:,i)-Nnnm,'color',color(i,:),'LineWidth',1);
%   label{i} = sprintf('amp of %.1f - mean',ampbinncnt(i));
% end
label{nbin+1} = sprintf('mean');  %median
if typepltnoi == 1 
  title(ax,'Synthetic noise');
  p(nbin+1)=plot(ax,cnt,Nnnm,'k-','LineWidth',1.5);
  legend(ax,p,label);
elseif typepltnoi == 2 
  title(ax,'Data - Synthetic noise');
  p(nbin+1)=plot(ax,cnt,Nnm-Nnnm,'k-','LineWidth',1.5);
end
ylabel(ax,'Normalized count');
xlabel(ax,sprintf('Diff. arrival between sources N and N-%d (s)',nsep));
xlim(ax,xran);
ylim(ax,yran);
longticks(ax,2);
hold(ax,'off');
% keyboard

ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% plot(ax,ampbincnt,dtarvlimed,'k-'); %if bin by amp with a equal width
% p1=plot(ax,log(ampbincnt),dtarvlimed1,'k-');  %if bin by amp with a equal number
% p2=plot(ax,log(ampbincnt),dtarvlimed2,'k--'); 
% p3=plot(ax,log(ampbincnt),dtarvlimed,'k-.'); 
% if nsep==1
%   yran=[0.25 0.75];
%   % yran=[0.1 0.5];
% %   legend(ax,[p1 p2 p3],'w/i 0.5 s','w/i 0.625 s','all');
% elseif nsep==2
%   yran=[0.1 0.6]; 
%   % yran=[0 0.4];
% %   legend(ax,[p1 p2 p3],'w/i 1 s','w/i 1.125 s','all');
% elseif nsep==3
%   yran=[0.5 4.5];
% %   legend(ax,[p1 p2 p3],'w/i 1.5 s','w/i 1.625 s','all');
% end
yran=[0 1];
plot(ax,log(ampbincnt),frac,'k-','linew',1,'marker','o','markersize',4,'markerfacec','k');
xlabel(ax,'Median log_{10}{amp}');
% ylabel(ax,'Median diff. arrival (s)');
ylabel(ax,sprintf('Frac. of diff. arrival w/i %.3f s',dtcut));
ylim(ax,yran);
title(ax,'Data');

ax=f.ax(4); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% plot(ax,ampbinncnt,dtarvlinmed,'k-');
% plot(ax,log(ampbinncnt),dtarvlinmed1,'k-');
% plot(ax,log(ampbinncnt),dtarvlinmed2,'k--'); 
% plot(ax,log(ampbinncnt),dtarvlinmed,'k-.');
if typepltnoi == 1 
  plot(ax,log(ampbinncnt),fracn,'k-','linew',1,'marker','o','markersize',4,'markerfacec','k');
  title(ax,'Synthetic noise');
elseif typepltnoi == 2 
  plot(ax,log(ampbincnt),fracn,'k-','linew',1,'marker','o','markersize',4,'markerfacec','k');
  title(ax,'Data - Synthetic noise');
end
xlabel(ax,'Median log_{10}{amp}');
% ylabel(ax,'Median diff. arrival (s)');
ylabel(ax,sprintf('Frac. of diff. arrival w/i %.3f s',dtcut));
ylim(ax,yran);

keyboard

orient(f.fig,'landscape');
fname = strcat(sprintf('nn%d%dbinsdifftime',nsep,nbin),fnsuffix,'.pdf');
print(f.fig,'-dpdf','-fillpage',strcat('/home/chaosong/Pictures/',fname));



%% in sample sapce, 2ndary sources removed
% allbstsig = allbstnoi;
dt2all = allbstsig.dt2allbst;
dist2allsp = allbstsig.dist2allspbst;
dloc2allsp = allbstsig.dloc2allspbst;
imp = allbstsig.impindepall;
f=plt_srcdlocinmap(dt2all,dloc2allsp,sps,'spl');

nsrc = allbstsig.nsrcraw;
dto2allbst=[];  %diff. origin time between event pairs
dloco2allbst=[];  
disto2allbst=[];
impstall = [];
tcorall = [];
for i = 1: size(trange,1)
  ist = sum(nsrc(1:i-1))+1;
  ied = ist+nsrc(i)-1;
  impi = imp(ist:ied,:);
  if ~isempty(impi)
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
    [dto2all,dloco2all,disto2all] = srcdistall(torisplst,impst(:,7:8),[0 2*sps]);
    
  else
    dto2all = [];
    dloco2all = [];
    disto2all = [];
    impst = [];
    dto2all = [];
    indsort = [];
    tcor = [];
  end
  dto2allbst = [dto2allbst; dto2all];
  dloco2allbst = [dloco2allbst; dloco2all] ;
  disto2allbst = [disto2allbst; disto2all];
  impstall = [impstall; impst];
  tcorall = [tcorall; tcor(indsort)];
end
f=plt_srcdlocinmap(dto2allbst,dloco2allbst,sps,'spl','tori');

%plot the diff time and euclidean distance between each LFE source to all others
tlensec = 14981;  %total length of time in secs for all bursts combined
f = plt_srcdistall(dt2all,dist2allsp,sps,40/sps,tlensec,1,'spl');
f = plt_srcdistall(dto2allbst,disto2allbst,sps,40/sps,tlensec,1,'spl','tori');

%plot the loc diff between each LFE source to all others
f = plt_srcdlocall(dloc2allsp,1,'spl');

doffset{1,1} = allbstsig.distarvlspnn1all(:,2:3);
doffset{2,1} = allbstsig.distarvlspnn2all(:,2:3);
doffset{3,1} = allbstsig.distarvlspnn3all(:,2:3);
doffset{4,1} = allbstsig.distarvlspnn4all(:,2:3);
doffset{5,1} = allbstsig.distarvlspnn5all(:,2:3);
eucdistsp{1,1} = allbstsig.distarvlspnn1all(:,1);
eucdistsp{2,1} = allbstsig.distarvlspnn2all(:,1);
eucdistsp{3,1} = allbstsig.distarvlspnn3all(:,1);
eucdistsp{4,1} = allbstsig.distarvlspnn4all(:,1);
eucdistsp{5,1} = allbstsig.distarvlspnn5all(:,1);
dtarvl{1,1} = allbstsig.dtarvlnn1all;
dtarvl{2,1} = allbstsig.dtarvlnn2all;
dtarvl{3,1} = allbstsig.dtarvlnn3all;
dtarvl{4,1} = allbstsig.dtarvlnn4all;
dtarvl{5,1} = allbstsig.dtarvlnn5all;
%plot the diff time and euclidean distance between above source pairs
f = plt_srcdistNtoNm(dtarvl,eucdistsp,sps,40/sps,tlensec,1,'spl');
%plot the loc diff between above source pairs
f = plt_srcdlocNtoNm(doffset,1,'spl');

%%%if you want to copy several figs to a new summary to put everything into 
% f1 = initfig(10,9,2,3);
% for i = 1:3
%   axchild=f.ax(i).Children;
%   copyobj(axchild, f1.ax(i));
% end

%% in sample sapce, further checked by the 4th station KLNB
dt2all = allbstsig.dt2all4thbst;
dist2allsp = allbstsig.dist2allsp4thbst;
dloc2allsp = allbstsig.dloc2allsp4thbst;
imp = allbstsig.impindep4thall;
f=plt_srcdlocinmap(dt2all,dloc2allsp,sps,'spl');

nsrc = allbstsig.nsrc;
dto2allbst=[];  %diff. origin time between event pairs
dloco2allbst=[];  
disto2allbst=[];
impstall = [];
tcorall = [];
for i = 1: size(trange,1)
  ist = sum(nsrc(1:i-1))+1;
  ied = ist+nsrc(i)-1;
  impi = imp(ist:ied,:);
  if ~isempty(impi)
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
    [dto2all,dloco2all,disto2all] = srcdistall(torisplst,impst(:,7:8),[0 2*sps]);
    
  else
    dto2all = [];
    dloco2all = [];
    disto2all = [];
    impst = [];
    dto2all = [];
    indsort = [];
    tcor = [];
  end
  dto2allbst = [dto2allbst; dto2all];
  dloco2allbst = [dloco2allbst; dloco2all] ;
  disto2allbst = [disto2allbst; disto2all];
  impstall = [impstall; impst];
  tcorall = [tcorall; tcor(indsort)];
end
f=plt_srcdlocinmap(dto2allbst,dloco2allbst,sps,'spl','tori');

%plot euclidean distance between each LFE source to all others
f = plt_srcdistall(dt2all,dist2allsp,sps,40/sps,tlensec,1,'spl');
f = plt_srcdistall(dto2allbst,disto2allbst,sps,40/sps,tlensec,1,'spl','tori');

%plot the loc diff between each LFE source to all others
f = plt_srcdlocall(dloc2allsp,1,'spl');

doffset{1,1} = allbstsig.distarvlspnn14thall(:,2:3);
doffset{2,1} = allbstsig.distarvlspnn24thall(:,2:3);
doffset{3,1} = allbstsig.distarvlspnn34thall(:,2:3);
doffset{4,1} = allbstsig.distarvlspnn44thall(:,2:3);
doffset{5,1} = allbstsig.distarvlspnn54thall(:,2:3);
eucdistsp{1,1} = allbstsig.distarvlspnn14thall(:,1);
eucdistsp{2,1} = allbstsig.distarvlspnn24thall(:,1);
eucdistsp{3,1} = allbstsig.distarvlspnn34thall(:,1);
eucdistsp{4,1} = allbstsig.distarvlspnn44thall(:,1);
eucdistsp{5,1} = allbstsig.distarvlspnn54thall(:,1);
dtarvl{1,1} = allbstsig.dtarvlnn14thall;
dtarvl{2,1} = allbstsig.dtarvlnn24thall;
dtarvl{3,1} = allbstsig.dtarvlnn34thall;
dtarvl{4,1} = allbstsig.dtarvlnn44thall;
dtarvl{5,1} = allbstsig.dtarvlnn54thall;
%plot the loc diff between above source pairs
f = plt_srcdlocNtoNm(doffset,1,'spl');
%plot the diff time and distance between above source pairs
f = plt_srcdistNtoNm(dtarvl,eucdistsp,sps,40/sps,tlensec,1,'spl');


%% in map view, 2ndary sources removed
dt2all = allbstsig.dt2allbst;
dist2all = allbstsig.dist2allbst;
dloc2all = allbstsig.dloc2allbst;
imp = allbstsig.impindepall;
f=plt_srcdlocinmap(dt2all,dloc2all,sps,'km');

nsrc = allbstsig.nsrcraw;
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
  [dto2all,dloco2all,disto2all] = srcdistall(torisplst,implocst,[0 2*sps]);
  dto2allbst = [dto2allbst; dto2all];
  dloco2allbst = [dloco2allbst; dloco2all] ;
  disto2allbst = [disto2allbst; disto2all];
  impstall = [impstall; impst];
  tcorall = [tcorall; tcor(indsort)];
end
f=plt_srcdlocinmap(dto2allbst,dloco2allbst,sps,'km','tori');

%plot euclidean distance between each LFE source to all others
tlensec = 14981;
f = plt_srcdistall(dt2all,dist2all,sps,40/sps,tlensec,0.1,'km');
%plot the loc diff between each LFE source to all others
f = plt_srcdlocall(dloc2all,0.1,'km');

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
%plot the loc diff between above source pairs
f = plt_srcdlocNtoNm(dneloc,0.1,'km');
%plot the diff time and distance between above source pairs
f = plt_srcdistNtoNm(dtarvl,eucdist,sps,40/sps,tlensec,0.1,'km');

%%
%%%projection along the min-rmse direction
locxyprojall = allbstsig.locxyprojall;
tarvlsplstall = allbstsig.impindepall(:,1);
nsrc = allbstsig.nsrcraw;
dlocproj2allbst = [];
% dt2allchkbst = [];
m = 5;
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
  [~,dlocproj2all] = srcdistall(tarvlsplst,locxyproj,[0 2*sps]);
%   dt2allchkbst = [dt2allchkbst; dt2allchk];
  dlocproj2allbst = [dlocproj2allbst; dlocproj2all];
  
  [~,dlocproj] = srcdistNtoNm(tarvlsplst,locxyproj,m);
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

% f=plt_srcdlocinmap(dt2all,dlocproj2allbst,sps,'km');

%plot the projected loc diff between each LFE source to all others
f = plt_srcdlocall(dlocproj2allbst,0.1,'km');
xlabel(f.ax(1),'Diff loc along min-rmse direc.');
xlabel(f.ax(2),'Abs diff loc along min-rmse direc.');
xlabel(f.ax(3),'Diff loc along ort direc. to min-rmse');
xlabel(f.ax(4),'Abs diff loc along ort direc. to min-rmse');
xlim(f.ax([1 3]),[-5 5]);
xlim(f.ax([2 4]),[0 5]);

%plot the projected loc diff between above source pairs
f = plt_srcdlocNtoNm(dlocprojbst,0.1,'km');
xlabel(f.ax(1),'Diff loc along min-rmse direc.');
xlabel(f.ax(2),'Abs diff loc along min-rmse direc.');
xlabel(f.ax(3),'Diff loc along ort direc. to min-rmse');
xlabel(f.ax(4),'Abs diff loc along ort direc. to min-rmse');
xlim(f.ax([1 3]),[-5 5]);
xlim(f.ax([2 4]),[0 5]);

%% in map view, further checked by the 4th station KLNB
dt2all = allbstsig.dt2all4thbst;
dist2all = allbstsig.dist2all4thbst;
dloc2all = allbstsig.dloc2all4thbst;
%plot euclidean distance between each LFE source to all others
f = plt_srcdistall(dt2all,dist2all,sps,40/sps,tlensec,0.1,'km');
%plot the loc diff between each LFE source to all others
f = plt_srcdlocall(dloc2all,0.1,'km');

dneloc{1,1} = allbstsig.distarvlnn14thall(:,2:3);
dneloc{2,1} = allbstsig.distarvlnn24thall(:,2:3);
dneloc{3,1} = allbstsig.distarvlnn34thall(:,2:3);
dneloc{4,1} = allbstsig.distarvlnn44thall(:,2:3);
dneloc{5,1} = allbstsig.distarvlnn54thall(:,2:3);
eucdist{1,1} = allbstsig.distarvlnn14thall(:,1);
eucdist{2,1} = allbstsig.distarvlnn24thall(:,1);
eucdist{3,1} = allbstsig.distarvlnn34thall(:,1);
eucdist{4,1} = allbstsig.distarvlnn44thall(:,1);
eucdist{5,1} = allbstsig.distarvlnn54thall(:,1);
dtarvl{1,1} = allbstsig.dtarvlnn14thall;
dtarvl{2,1} = allbstsig.dtarvlnn24thall;
dtarvl{3,1} = allbstsig.dtarvlnn34thall;
dtarvl{4,1} = allbstsig.dtarvlnn44thall;
dtarvl{5,1} = allbstsig.dtarvlnn54thall;
%plot the loc diff between above source pairs
f = plt_srcdlocNtoNm(dneloc,0.1,'km');
%plot the diff time and distance between above source pairs
f = plt_srcdistNtoNm(dtarvl,eucdist,sps,40/sps,tlensec,0.1,'km');

%%%projection along the min-rmse direction
locxyprojall = allbstsig.locxyproj4thall;
tarvlsplstall = allbstsig.impindep4thall(:,1);
nsrc = allbstsig.nsrc;
dlocproj2allbst = [];
% dt2allchkbst = [];
m = 5;
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
  [~,dlocproj2all] = srcdistall(tarvlsplst,locxyproj,[0 2*sps]);
%   dt2allchkbst = [dt2allchkbst; dt2allchk];
  dlocproj2allbst = [dlocproj2allbst; dlocproj2all];
  
  [~,dlocproj] = srcdistNtoNm(tarvlsplst,locxyproj,m);
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

%plot the projected loc diff between each LFE source to all others
f = plt_srcdlocall(dlocproj2allbst,0.1,'km');
xlabel(f.ax(1),'Diff loc along min-rmse direc.');
xlabel(f.ax(2),'Abs diff loc along min-rmse direc.');
xlabel(f.ax(3),'Diff loc along ort direc. to min-rmse');
xlabel(f.ax(4),'Abs diff loc along ort direc. to min-rmse');
xlim(f.ax([1 3]),[-5 5]);
xlim(f.ax([2 4]),[0 5]);

%plot the projected loc diff between above source pairs
f = plt_srcdlocNtoNm(dlocprojbst,0.1,'km');
xlabel(f.ax(1),'Diff loc along min-rmse direc.');
xlabel(f.ax(2),'Abs diff loc along min-rmse direc.');
xlabel(f.ax(3),'Diff loc along ort direc. to min-rmse');
xlabel(f.ax(4),'Abs diff loc along ort direc. to min-rmse');
xlim(f.ax([1 3]),[-5 5]);
xlim(f.ax([2 4]),[0 5]);




%%
%%%the direct/scaled deconvolved pos/neg source peak ratio between all station pairs, for each
%%%burst win separately
nsrc = allbstsig.nsrc;
msrcampr = allbstsig.msrcampr;
madsrcampr = allbstsig.madsrcampr;
f1 = initfig(16,5,1,size(msrcampr,2));
f1 = plt_deconpk_rat4th(f1,msrcampr,madsrcampr,nsrc,'b');

nsrc = allallallallbstnoi.nsrc;
msrcampr = allallallallbstnoi.msrcampr;
madsrcampr = allallallallbstnoi.madsrcampr;
f1 = plt_deconpk_rat4th(f1,msrcampr,madsrcampr,nsrc,'r');

%%
%%%combine the direct/scaled deconvolved pos/neg source peak ratio between all station pairs of
%%%all burst wins, and summarize into one histogram
srcamprall = allbstsig.srcamprall;
impindep4thall = allbstsig.impindep4thall;
f2 = initfig(16,9,2,size(srcamprall,2));
f2 = plt_deconpk_rat_comb4th(f2,srcamprall,impindep4thall,'b');

srcamprall = allallallallbstnoi.srcamprall;
impindep4thall = allallallallbstnoi.impindep4thall;
f2 = plt_deconpk_rat_comb4th(f2,srcamprall,impindep4thall,'r');

%%
%%%deviation of source amp ratio from some median vs. RCC
lgdevsrcamprall = allbstsig.lgdevsrcamprall;
rccpairsrcall = allbstsig.rccpairsrcall;
rcccatsrcall = allbstsig.rcccatsrcall;
f3 = initfig(15,7,2,size(lgdevsrcamprall,2)); %initialize fig
f3 = plt_deconpk_ratdevvsrcc4th(f3,lgdevsrcamprall,rccpairsrcall,rcccatsrcall,'b');

lgdevsrcamprall = allallallallbstnoi.lgdevsrcamprall;
rccpairsrcall = allallallallbstnoi.rccpairsrcall;
rcccatsrcall = allallallallbstnoi.rcccatsrcall;
f3 = plt_deconpk_ratdevvsrcc4th(f3,lgdevsrcamprall,rccpairsrcall,rcccatsrcall,'r');

%%
%%%diff between predicted arrival and selected peak at 4th sta
pred4offtrall = allbstsig.pred4offtrall;
f4 = initfig(5,5,1,size(pred4offtrall,2)); %initialize fig
f4 = plt_errorof4thtarvlpred(f4,pred4offtrall,offmax,'b');

pred4offtrall = allallallallbstnoi.pred4offtrall;
f4 = plt_errorof4thtarvlpred(f4,pred4offtrall,offmax,'r');

%%
%%%preserved sources' amp ratio between 4th and 1st stas
srcamprall = allbstsig.srcamprall;
impindep4thall = allbstsig.impindep4thall;
f5 = initfig(12,5,1,3); %initialize fig
f5 = plt_deconpk_rat14(f5,impindep4thall,srcamprall,'b');

srcamprall = allallallallbstnoi.srcamprall;
impindep4thall = allallallallbstnoi.impindep4thall;
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

clppkhtwfall = allallallallbstnoi.clppkhtwfall;
psrcampsall = allallallallbstnoi.psrcampsall;
clnpkhtwfall = allallallallbstnoi.clnpkhtwfall;
nsrcampsall = allallallallbstnoi.nsrcampsall;
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
nsrcn = allallallallbstnoi.nsrc;  %number of sources from synthetic noise

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

dtarvlnn1all = allallallallbstnoi.dtarvlnn1all;
distarvlnn1all = allallallallbstnoi.distarvlnn1all;
f1 = plt_absdistvsdt(f1,[dtarvlnn1all,distarvlnn1all],1,sps,'tarvl','r');
ylim(f1.ax(1),[0 0.6]);
ylim(f1.ax(2),[0 3]);









