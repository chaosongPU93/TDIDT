% circsrcmodel.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the script to assume a circular source model. Now we have locations
% of all deconvolved LFE sources where each source is represented by a circle
% with a series of diameter. This script aims to show several things:
% 1. for any cluster of m+1 consecutive events, what is the distribution
% 2. To summarize all these images, consider the ratio of overlapping area.
% with a certain diameter, sources may or may not overlap. in total, they have
% a covered area. Meanwhile, there is a nominal area which is the sum of the
% area of each circle. 1 minus this ratio of area indicates the overlapping 
% extent.
% 
% --2024/02/06, mapping from samples to relation location already contains error,
% so when you try to obtain the location difference, or distance, start with 
% the difference in samples first. For example, if you ask the difference in 
% certain directions, first obtain the difference in off12 and off13, then map
% the diff to location, then you can project along any direction to get the
% loc difference along that direction, or the distance (abs of loc diff) in
% that direction, etc.
%
% --2024/02/15, looks like if you need the loc difference in map-view, you have
% to first invert time offsets to map locations, then obtain the difference.
% Despite the possible location error propagation, you need it. Since if you 
% start from difference in samples, then transform them to map locations,
% you are essentially assuming one of the two sources is located at the 
% origin, which is not the case in reality, even if the actual difference in
% the two ways may be small.
%
% --Lots of information is shown in the resulting plot of each cluster,
% including the distance between N and N-1, and each to others, projection
% direction, etc.
% --Other stuff, like the location difference (sign preserved) are now moved
% to 'lfedloc_incuster.m' for lumped statistics.

%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/01/04
% Last modified date:   2024/02/15
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
trange = load(strcat(rstpath, '/MAPS/tdec.bstran',num2str(ttol),'s.',...
  num2str(ntol),'.pgc002.',cutout(1:4)));
tlen = trange(:,3)-trange(:,2);
nbst = size(trange,1);
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

%choice of diameter of circular sources
% diams = [0.1 0.3 0.5 0.7];
diams = 0.75;

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

%% various types of distance between events in the clusters
m=6;  % a cluster of m+1 events
catclusm = catclus{m};
ncluster = size(catclusm,1); %num of clusters, consecu. clusters may share events!

%%%compute the distance w/i each cluster
%this would take a while
tic
[mdistnn1,mdist2all,pca2vec,projang1,mdprojx1nn1,mdprojx12all,...
  projang2,mdprojx2nn1,mdprojx22all,distnn1cat,dist2allcat,...
  dprojxy1nn1cat,dprojxy12allcat,dprojxy2nn1cat,dprojxy22allcat]=...
  dist_evtcluster(catclusm,sps,ftrans,timetype);
toc


%% various types of distance between events in the clusters
% % for m=1:mmax
% m=5;  % a cluster of m+1 events
% timetype = 'tarvl';
% % timetype = 'tori';
% dtcut = 0.25*m+0.125;
% [ampplt,dtplt,indimpdtcut] = med_amp_incluster(nbst,imp,nsrc,m,dtcut,sps,timetype);
% 
% %%%compute the distance w/i each cluster
% %this would take a while
% tic
% [mdistnn1,mdist2all,pca2vec,projang1,mprojx1nn1,mprojx12all,...
%   projang2,mprojx2nn1,mprojx22all]=...
%   dist_evtcluster(impcluster,implocclus,sps,ftrans,m,timetype);
% % [mdistnn1,mdist2all,pca2vec,projang1,mprojx1nn1,mprojx12all]=...
% %   dist_evtcluster(impcluster,implocclus,sps,ftrans,m,'tori');
% toc
% % keyboard

%% summarize the median distance comparison between consecutive ones, and each to all others
% f = initfig(12,4.5,1,3);
% binw=0.1;
% ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% p2=histogram(ax,dist2allcat,'binw',binw,'normalization','count','Facec','r','edgec','none');
% p1=histogram(ax,distnn1cat,'binw',binw,'normalization','count','Facec','b','edgec','none');
% plot(ax,[median(distnn1cat) median(distnn1cat)],ax.YLim,'b--','LineWidth',1);
% plot(ax,[median(dist2allcat) median(dist2allcat)],ax.YLim,'r--','LineWidth',1);
% xlabel(ax,'Absolute distance (km)');
% ylabel(ax,'Count');
% legend(ax,[p1,p2],'N and N-1','each to others');
% ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% plt2 = abs(dprojxy12allcat(:,1));
% histogram(ax,plt2,'binw',binw,'normalization','count','Facec','r','edgec','none');
% plt1 = abs(dprojxy1nn1cat(:,1));
% histogram(ax,plt1,'binw',binw,'normalization','count','Facec','b','edgec','none');
% plot(ax,[median(plt1) median(plt1)],ax.YLim,'b--','LineWidth',1);
% plot(ax,[median(plt2) median(plt2)],ax.YLim,'r--','LineWidth',1);
% xlabel(ax,'Distance along the 2nd PC (km)');
% ylabel(ax,'Count');
% ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% plt2 = abs(dprojxy22allcat(:,1));
% histogram(ax,plt2,'binw',binw,'normalization','count','Facec','r','edgec','none');
% plt1 = abs(dprojxy2nn1cat(:,1));
% histogram(ax,plt1,'binw',binw,'normalization','count','Facec','b','edgec','none');
% plot(ax,[median(plt1) median(plt1)],ax.YLim,'b--','LineWidth',1);
% plot(ax,[median(plt2) median(plt2)],ax.YLim,'r--','LineWidth',1);
% xlabel(ax,'Distance along the prop. direc. (km)');
% ylabel(ax,'Count');
% % keyboard

%% normalized total covered area for the clustered events assuming a diameter
% for ii = 1: length(diams)
%   tic
%   diam=diams(ii); %dimension of source, would be diamter if circle
%   radi=0.5*diam; %radius of blue circles, max. with no overlapping
%   for i = 1: ncluster
%     imp = catclusm{i};
%     imploc = off2space002(imp(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
%     arearat(i,ii) = area_of_overlap_circs_grid(imploc(:,1:2),radi,0.01,0);
%   end
%   toc
% end

%% plot the distribution of total covered area normalized by norminal area
% f = initfig(8,8,2,2); %initialize fig
% binw=0.02;
% binedge=0:binw:1;
% for ii = 1:length(diams)
%   diam=diams(ii);
%   ax=f.ax(ii); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%   histogram(ax,arearat(:,ii),'BinEdges',binedge);
%   text(ax,0.98,0.95,sprintf('diameter: %.1f km',diam),'Units','normalized',...
%       'HorizontalAlignment','right');
%   xlabel(ax,'Normalized area');
%   ylabel(ax,'Count');
%   xlim(ax,[0 1]);
%   %there is a lower limit of the ratio of covered area
%   plot(ax,[1/(m+1) 1/(m+1)],ax.YLim,'k--');
% end
% supertit(f.ax(1:2),sprintf('m=%d, ncluster=%d',m,ncluster),10);
% % keyboard
% figname = sprintf('arearat_m%d_nclus%d.pdf',m,ncluster);
% print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',figname));

%% plot locations of srcs with finite size in each clusters assuming a diameter
close all

ampm = catmedamp{m};
amp=mean(imp(:,[2 4 6]),2);
ampch=prctile(amp,5);
nevt = m+1;

for ii = 1: length(diams)
  diam=diams(ii);
  radi=0.5*diam;

  nrow1 = 5; % rows and cols of subplots in each figure
  ncol1 = 6; 
  widin1 = 10; % size of each figure
  htin1 = 8.4;
  pltxran1 = [0.04 0.98]; pltyran1 = [0.04 0.98]; % optimal axis location
  pltxsep1 = 0.008; pltysep1 = 0.01;
  ifig1 = 0; % figure number

  ifig1 = ifig1+1;
  f1 = initfig(widin1,htin1,nrow1,ncol1,ifig1);
  axpos1 = optaxpos(f1,nrow1,ncol1,pltxran1,pltyran1,pltxsep1,pltysep1);
  axtit = f1.ax(1: ncol1);
%   supertit(axtit, sprintf('Fig %s, nevt=%d, nclus=%d, diam=%.1f km',...
%     num2zeropadstr(ifig1, 3),nevt,ncluster,diam),10);
  xlabel(f1.ax((nrow1-1)*ncol1+1),'E (km)','fontsize',10);
  ylabel(f1.ax((nrow1-1)*ncol1+1),'N (km)','fontsize',10);
  orient(f1.fig,'landscape');
  isub = 0;  % count of the total subplots on the current figure f1

  for i = 1: ncluster
    %events in each cluster
    impplt = catclusm{i};
    tevt=impplt(:,1);
    implocplt = off2space002(impplt(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
    ampplt = ampm(:,1);
    
    isub = isub+1; 
    if isub > nrow1*ncol1
      %print the current figure f1
      ind = setdiff(1:nrow1*ncol1, (nrow1-1)*ncol1+1);
      nolabels(f1.ax(ind),3);
      f1name{ifig1,1} = sprintf('clusterloc_Fig%s_m%d_%.2fkm.pdf',...
        num2zeropadstr(ifig1, 3),m,diam);
      print(f1.fig,'-dpdf',strcat('/home/chaosong/Pictures/',f1name{ifig1,1}));
      %move on to next figure f1
      ifig1 = ifig1+1;  %initialize another figure
      f1 = initfig(widin1,htin1,nrow1,ncol1,ifig1);
      axpos1 = optaxpos(f1,nrow1,ncol1,pltxran1,pltyran1,pltxsep1,pltysep1);
      axtit = f1.ax(1: ncol1);
%       supertit(axtit, sprintf('Fig %s, nevt=%d, nclus=%d, diam=%.1f km',...
%         num2zeropadstr(ifig1, 3),nevt,ncluster,diam),10);
      xlabel(f1.ax((nrow1-1)*ncol1+1),'E (km)','fontsize',10);
      ylabel(f1.ax((nrow1-1)*ncol1+1),'N (km)','fontsize',10);
      orient(f1.fig,'landscape');
      
      isub = isub-nrow1*ncol1;  %refresh
    end

    %in terms of origin time
    if strcmp(timetype,'tori')
      [torisplst,~,indsort]=tarvl2tori(impplt,sps,ftrans,1);  %return the origin time 
      implocplt = implocplt(indsort, :);
      tevt = torisplst;
    end
    
    %plot locations of events in each cluster
%     color=jet(nevt);
%     color=gradientblue(nevt);
    color=plasma(nevt);
    ax=f1.ax(isub); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    axis(ax, 'equal');
    for j=1: nevt
      [xcut,ycut] = circle_chao(implocplt(j,1),implocplt(j,2),radi,0.1);
      plot(ax,xcut,ycut,'-','color',color(j,:),'linew',1);
    end
    x0=mean(implocplt(:,1));
    y0=mean(implocplt(:,2));
    plot(ax,x0+[pca2vec(i,1), -pca2vec(i,1)],y0+[pca2vec(i,2), -pca2vec(i,2)],...
      'k-','linew',1);
    text(ax,0.02,0.95,sprintf('%d/%d, m=%d',i,ncluster,nevt),...
      'Units','normalized','HorizontalAlignment','left','fontsize',7);   
%     text(ax,0.02,0.95,sprintf('amp:%.2f; %.1fx',ampplt(i),ampplt(i)/ampch),...
%       'Units','normalized','HorizontalAlignment','left','fontsize',7);   
    % text(ax,0.02,0.95,sprintf('arearat: %.2f',arearat(i,ii)),...
    %   'Units','normalized','HorizontalAlignment','left','fontsize',7);
    text(ax,0.99,0.19,sprintf('Abs: %.2f; %.2f km',mdistnn1(i,1),mdist2all(i,1)),...
      'Units','normalized','HorizontalAlignment','right','fontsize',7);
    text(ax,0.99,0.12,sprintf('%d%c: %.2f; %.2f km',round(projang1(i,1)),...
      char(176),mdprojx1nn1(i,1),mdprojx12all(i,1)),...
      'Units','normalized','HorizontalAlignment','right','fontsize',7);
    text(ax,0.99,0.05,sprintf('%d%c: %.2f; %.2f km',round(projang2(i,1)),...
      char(176),mdprojx2nn1(i,1),mdprojx22all(i,1)),...
      'Units','normalized','HorizontalAlignment','right','fontsize',7);
    xlim(ax,[-4 4]);
    ylim(ax,[-4 4]);
    longticks(ax,2);
%     if isub ~= (nrow1-1)*ncol1+1
%       nolabels(ax,3);
%     end
    hold(ax,'off');
  end
  
  if isub<nrow1*ncol1
    delete(f1.ax(isub+1:nrow1*ncol1));
    if rem(isub,ncol1) ==0
      indax=(floor(isub/ncol1)-1)*ncol1+1;
    else
      indax=floor(isub/ncol1)*ncol1+1;
    end
    xlabel(f1.ax(indax),'E (km)','fontsize',10);
    ylabel(f1.ax(indax),'N (km)','fontsize',10);
    xticks(f1.ax(indax),-4:2:4);
    yticks(f1.ax(indax),-4:2:4);
  end
  
  for iii=1:isub
    if iii ~= floor(isub/ncol1)*ncol1+1
      nolabels(f1.ax(iii),3);
    end
  end
  
  %if all clusters have been plotted, print the current figure f1
  f1name{ifig1,1} = sprintf('clusterloc_Fig%s_m%d_%.2fkm.pdf',...
    num2zeropadstr(ifig1, 3),m,diam);
  print(f1.fig,'-dpdf',strcat('/home/chaosong/Pictures/',f1name{ifig1,1}));

  %%
  %merge all figures into a single pdf file
  mergef1name=sprintf('clusterloc_m%d_%.2fkm.pdf',m,diam);
  str=strcat('rm -f /home/chaosong/Pictures/',mergef1name);
  status = system(str);
  % status = system('rm -f /home/chaosong/Pictures/clusterloc.pdf');
  for jj = 1:size(f1name,1)
    fname{jj} = strcat('/home/chaosong/Pictures/',f1name{jj});
  end
  append_pdfs(strcat('/home/chaosong/Pictures/',mergef1name),fname);
%   status = system('rm -f /home/chaosong/Pictures/clusterloc_Fig*.pdf');
  % keyboard
%   close all
end
  




