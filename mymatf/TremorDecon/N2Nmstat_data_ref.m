% N2Nmstat_data_ref.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the script to carry out a bunch of N to N-m statisical analysis
% to deconvoluted catalog from the real data. The codes are cut from 
% 'deconv_ref_4s_exp_4thsta.m' because they are getting more and more 
% complicated and different from the original scope of 
% 'deconv_ref_4s_exp_4thsta.m'
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
tlen = trange(:,3)-trange(:,2);
tlensum = sum(tlen);

modname = 'timeoff_plfit_4thsta_160sps.mat';
planefit = load(strcat(rstpath, '/MAPS/',modname));
rmse = planefit.gof{4}.rmse;
offmax = round(2.0*rmse);

sps = 160;

%%%load data
savefile = 'deconv_stats4th_allbstsig.mat';
load(strcat(rstpath, '/MAPS/',savefile));
savefile = 'deconv_stats4th_allbstnoi.mat';
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

typepltnoi = 1; %plot noise
% typepltnoi = 2; %plot data -noise

keyboard

%% for a few m, fraction of unique events inside clusters w/i a diff time cut, wrt. all catalog for data
mmax = 15;
nbst = size(trange,1);

[fracsrc2all, dfracsrc2all, f]=frac_uniqevt_incluster(nbst,imp,nsrc,mmax,sps);
[fracsrc2alln, dfracsrc2alln, f]=frac_uniqevt_incluster(nbst,impn,nsrcn,mmax,sps);

% keyboard

%% for a few m, diff time distribution between N&N-m, and fraction w/i dtcut  
f1 = initfig(10,4.5,1,2); %initialize fig
tit=supertit(f1.ax,strcat({'Diff. arrival between N & N-m, '},supertstr));
movev(tit,0.3);
mmax = 15;
nbst = size(trange,1);
[f,cnt,Nn,Nnn,frac,fracn,dtarvl,dtarvln,mmaxnonzero,mmaxnonzeron]=...
  plt_srcdifftime_NNm_mmax(f1,nbst,imp,nsrc,impn,nsrcn,mmax,sps,typepltnoi);

%% for a few m, bin by median amp for all events within N&N-m, fraction of diff time measurements w/i dtcut  
f2 = initfig(10,9,2,2); %initialize fig
tit=supertit(f2.ax,strcat({'Diff. arrival between N & N-m, binned by amp, '},supertstr));
movev(tit,0.2);
[f2,ampbincnt,Nbn,fracb,Nbnn,fracbn]=...
  plt_fracdifftime_NNm_mmax(f2,nbst,imp,nsrc,impn,nsrcn,mmax,sps,typepltnoi);

keyboard


%% for a certain m, similar plots for N & N-m pairs 
%%%Median amplitude for srcs w/i the cluster
m = 1;
nbst = size(trange,1);
dtcut = 0.25*m+0.125;
[ampplt,dtarvlplt]=med_amp_incluster(nbst,imp,nsrc,m);
[amppltn,dtarvlpltn]=med_amp_incluster(nbst,impn,nsrcn,m);

%%
%%%summarize the whole catalog, diff arrival time and fractions, data and noise together
f1 = initfig(5.5,4,1,1); %initialize fig
tit=supertit(f1.ax,supertstr);
movev(tit,0.2);

f2 = initfig(12,4,1,2); %initialize fig
tit=supertit(f2.ax,supertstr);
movev(tit,0.2);

[f1,f2,Nn,frac,Nnn,fracn,fracdif]=...
  plt_difftime_NNm(f1,f2,dtarvlplt,dtarvlpltn,dtcut,sps,typepltnoi,m);

%%
%%%Bin by amp, then plot diff time distribution, and frac w/i some 'dtcut'
f = initfig(12,8,2,2); %initialize fig
tit=supertit(f.ax,supertstr);
movev(tit,0.2);
[f,Nn,fracm,ampbincntm,Nnn,fracnm,ampbinncntm]=...
  plt_fracdifftime_NNm(f,ampplt,dtarvlplt,amppltn,dtarvlpltn,dtcut,sps,typepltnoi,m);

keyboard

orient(f.fig,'landscape');
fname = strcat(sprintf('nn%d%dbinsdifftime',m,nbin),fnsuffix,'.pdf');
print(f.fig,'-dpdf','-fillpage',strcat('/home/chaosong/Pictures/',fname));

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

