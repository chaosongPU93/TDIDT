% plt_rccdatavsnoise.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This piece of code is to plot the comparison of RCC between data and noise,
% potentially it also provides a template for comparing with synthetics. 
% From decon routines, we computed the RCC of both. Here we try to plot the
% empirical PDF and CDF.
% The PDF is estimated by 'ksdensity.m', it can be done to each burst, and a 
% median of all bursts can be extracted. Or, lump all bursts then compute.
% CDF is estimated by 'ecdf', or a plain way (see code). Similarly, it can be 
% done to each burst, and a median of all bursts can be extracted. Or, lump all 
% bursts then compute.
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/06/24
% Last modified date:   2024/06/24
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

%%%load the LFE catalog, 25-s-win
% savefile = 'deconv_stats4th_allbstsig.mat';
savefile = 'deconv_stats4th_no23_allbstsig.mat';
allsig = load(strcat(rstpath, '/MAPS/',savefile));
rccbst = allsig.allbstsig.rccbst; %cat rcc; whole-win rcc; same for noise
ccwpairk = allsig.allbstsig.ccwpairk; %zero-lag CC of 25-s-window
mrccwpairk = allsig.allbstsig.mrccwpairk; %median RCC of 25-s-windows

%%%for RCC, noise catalog contains data and noise, 25-s-win, noise
% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
allnoi = load(strcat(rstpath, '/MAPS/',savefile));
rccbstn = allnoi.allbstnoi.rccbst; %cat rcc; whole-win rcc; same for noise
ccwpairkn = allnoi.allbstnoi.ccwpairk; %zero-lag CC of 25-s-window
mrccwpairkn = allnoi.allbstnoi.mrccwpairk; %median RCC of 25-s-windows

%%%load the LFE catalog, whole-win
% savefile = 'deconv1win_stats4th_allbstsig.mat';
savefile = 'deconv1win_stats4th_no23_allbstsig.mat';
allsig1win = load(strcat(rstpath, '/MAPS/',savefile));
rccbst1win = allsig1win.allbstsig.rccbst; %whole-win rcc; same for noise
ccpair1win = allsig1win.allbstsig.ccpair; %whole-win 0-lag cc paris
mcc1win = allsig1win.allbstsig.mcc; %mean whole-win 0-lag cc

%%
%%%load the RCC and CC value from real data, all bursts
% load(strcat('rcc',num2str(rccmwsec),'.mat'));
templump = cell(nbst,1);  % lumped RCC from all bursts, from 25-s windows
templump1 = cell(nbst,1);  % lumped RCC from all bursts, from whole-windows
templump2 = cell(nbst,1); % lumped zero-lag CC of whole-window from all bursts
templump3 = cell(nbst,1); % lumped zero-lag CC of 25-s-window from all bursts
templump4 = cell(nbst,1);  % lumped RCC from all bursts, from 25-s windows, NOISE
templump5 = cell(nbst,1);  % lumped zero-lag CC from all bursts, from 25-s windows, NOISE
for i = 1: nbst
  temp = rccbst{i};
  %%%cat rcc from 25-s-wins
  temp1 = sort(temp(:,1),'ascend');
  temp1(:,2) = (1:size(temp1,1))/size(temp1,1);
  templump{i} = temp1;  
  [cdfval,x] = ecdf(temp(:,1));
  templumpc{i} = [x cdfval];
  [pdfval,x] = ksdensity(temp(:,1),(-1:1e-2:1)');
  templumpp{i} = [x pdfval];
  
  %%%0-lag cc from 25-s-wins
  temp3 = ccwpairk{i};  %all 25-s subwins, pairs 12, 13, 23, as col
  temp3 = (temp3(:,1)+temp3(:,2))./2; %retain the mean of 12 and 13
  templump3{i} = temp3;
  
  %%%rcc from the whole-win
  temp1 = sort(temp(:,2),'ascend');
  temp1(:,2) = (1:size(temp1,1))/size(temp1,1);
  templump1{i} = temp1;
  [cdfval,x] = ecdf(temp(:,2));
  templump1c{i} = [x cdfval];
  [pdfval,x] = ksdensity(temp(:,2),(-1:1e-2:1)');
  templump1p{i} = [x pdfval];
  
  %%%0-lag cc from the whole-win
  temp2 = mcc1win(i); % the mean of 12 and 13 of the whole win,
  templump2{i} = temp2;

  %%%%%%%% for noise
  %%%cat rcc from 25-s-wins
  temp = rccbstn{i};
  temp4 = sort(temp(:,3),'ascend');
  temp4(:,2) = (1:size(temp4,1))/size(temp4,1);
  templump4{i} = temp4;  
  [cdfval,x] = ecdf(temp(:,3));
  templump4c{i} = [x cdfval];
  [pdfval,x] = ksdensity(temp(:,3),(-1:1e-2:1)');
  templump4p{i} = [x pdfval];
  
  %%%0-lag cc from 25-s-wins
  temp5 = ccwpairkn{i};  %all 25-s subwins, pairs 12, 13, 23, as col
  temp5 = (temp5(:,1)+temp5(:,2))./2; %retain the mean of 12 and 13
  templump5{i} = temp5; %lump them
  %%%%%%%% for noise
end
%%
%%%bin based on the lumped RCC for some 'median' value, 25-s-win version
rccplt = cat(1,templump{:}); %data
%%%2. bin x by equal number
nbin = 200;
[xbin,indbin,n] = binxeqnum(rccplt(:,1),nbin);
for i = 1: nbin
  ind = indbin{i};
  rccreal(i,1) = median(xbin{i});
  rccreal(i,2) = median(rccplt(ind,2));
end

rcclumpcdf = sort(rccplt(:,1),'ascend');
rcclumpcdf(:,2) = (1:size(rcclumpcdf,1))/size(rcclumpcdf,1);
[pdfval,x] = ksdensity(rccplt(:,1),(-1:1e-3:1)');
rcclumppdf = [x pdfval];


rccpltc = cat(1,templumpc{:}); %data
[xbin,indbin,n] = binxeqnum(rccpltc(:,1),nbin);
for i = 1: nbin
  ind = indbin{i};
  rccrealc(i,1) = median(xbin{i});
  rccrealc(i,2) = median(rccpltc(ind,2));
end

for i = 1:nbst
  aa=templumpp{i};
  rccpltp(:,i) = aa(:,2); %data
end
rccrealp(:,1) = aa(:,1);
rccrealp(:,2) = median(rccpltp,2);

%%%same as above, but for noise
rccnplt = cat(1,templump4{:}); %data
%%%2. bin x by equal number
[xbin,indbin,n] = binxeqnum(rccnplt(:,1),nbin);
for i = 1: nbin
  ind = indbin{i};
  rccnoi(i,1) = median(xbin{i});
  rccnoi(i,2) = median(rccnplt(ind,2));
end

rccnlumpcdf = sort(rccnplt(:,1),'ascend');
rccnlumpcdf(:,2) = (1:size(rccnlumpcdf,1))/size(rccnlumpcdf,1);
[pdfval,x] = ksdensity(rccnplt(:,1),(-1:1e-3:1)');
rccnlumppdf = [x pdfval];


rccnpltc = cat(1,templump4c{:}); %data
[xbin,indbin,n] = binxeqnum(rccnpltc(:,1),nbin);
for i = 1: nbin
  ind = indbin{i};
  rccnoic(i,1) = median(xbin{i});
  rccnoic(i,2) = median(rccnpltc(ind,2));
end

for i = 1:nbst
  aa=templump4p{i};
  rccnpltp(:,i) = aa(:,2); %data
end
rccnoip(:,1) = aa(:,1);
rccnoip(:,2) = median(rccnpltp,2);

% %median lumped rcc from all data bursts, 25-s-win version
% mrccreal = median(templump(:,1));
% %median lumped rcc from all data bursts, whole-win version
% mrccreal1win = median(templump1(:,1));
% %median lumped 0-lag CC from all data bursts, 25-s-win version
% mccreal = median(templump3);
% %median lumped 0-lag CC from all data bursts, whole-win version
% mccreal1win = median(templump2);
% %median lumped rcc from all data bursts, 25-s-win version, NOISE
% mrccnoi = median(templump4(:,1));
% %median lumped 0-lag CC from all data bursts, 25-s-win version
% mccnoi = median(templump5);

%%
widin = 7;  % maximum width allowed is 8.5 inches
htin = 3.5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f = initfig(widin,htin,nrow,ncol); %initialize fig
xran = [0.08 0.98]; yran = [0.15 0.95];
xsep = 0.08; ysep = 0.04;
optaxpos(f,nrow,ncol,xran,yran,xsep,ysep);

%empirical PDF
ax=f.ax(1);
hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% plot(ax,rccrealp(:,1),rccrealp(:,2),'k-','LineWidth',1.5);
% plot(ax,rccnoip(:,1),rccnoip(:,2),'r-','LineWidth',1.5);
p(1)=plot(ax,rcclumppdf(:,1),rcclumppdf(:,2),'k-','LineWidth',1.5);
p(2)=plot(ax,rccnlumppdf(:,1),rccnlumppdf(:,2),'r-','LineWidth',1.5);
label{1}='median of data';
label{2}='median of noise';
lgd=legend(ax,p,label,'Location','south');
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));
xlabel(ax,'RCC');
ylabel(ax,'PDF');
xticks(ax,-1: 0.5: 1);

%empirical CDF
ax=f.ax(2);
hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% p(1)=plot(ax,rccreal(:,1),rccreal(:,2),'k-','LineWidth',1.5);
% plot(ax,rccrealc(:,1),rccrealc(:,2),'b--','LineWidth',1);
% [~,ind]=min(abs(rccreal(:,2)-0.5));
% scatter(ax,rccreal(ind,1),0.5,30,'kv');
% p(2)=plot(ax,rccnoic(:,1),rccnoic(:,2),'r-','LineWidth',1.5);
% plot(ax,rccnoic(:,1),rccnoic(:,2),'c--','LineWidth',1);
% [~,ind]=min(abs(rccnoi(:,2)-0.5));
% scatter(ax,rccnoi(ind,1),0.5,30,'rv');
p(1)=plot(ax,rcclumpcdf(:,1),rcclumpcdf(:,2),'k-','LineWidth',1.5);
[~,ind]=min(abs(rcclumpcdf(:,2)-0.5));
scatter(ax,rcclumpcdf(ind,1),0.5,50,'k^','markerfacecolor','w');
p(2)=plot(ax,rccnlumpcdf(:,1),rccnlumpcdf(:,2),'r-','LineWidth',1.5);
[~,ind]=min(abs(rccnlumpcdf(:,2)-0.5));
scatter(ax,rccnlumpcdf(ind,1),0.5,50,'r^','markerfacecolor','w');
xlabel(ax,'RCC');
ylabel(ax,'CDF');
xticks(ax,-1: 0.5: 1);

fname = 'rcccdf.pdf';
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

% ax=f.ax(1);
% hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% for i=1:nbst
%   rcci=templump{i};
%   plot(ax,rcci(:,1),rcci(:,2),'-','color',[0 0 0  0.15],'LineWidth',1);
% end
% for i=1:nbst
%   rccni=templump4{i};
%   plot(ax,rccni(:,1),rccni(:,2),'-','color',[1 0 0  0.15],'LineWidth',1);
% end
% p(1)=plot(ax,rccreal(:,1),rccreal(:,2),'k-','LineWidth',1);
% p(2)=plot(ax,rccnoi(:,1),rccnoi(:,2),'r-','LineWidth',1);
% % label{1}='Med of data';
% % label{2}='Med of noise';
% % legend(ax,p,label,'Location','best');
% xlabel(ax,'RCC');
% ylabel(ax,'CDF');
% xticks(ax,-1: 0.2: 1);

% 
% ax=f.ax(1);
% hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% scatter(ax,rccplt(:,1),rccplt(:,2),5,'ko','filled');
% % scatter(ax,rccnplt(:,1),rccnplt(:,2),5,'ko','filled');





