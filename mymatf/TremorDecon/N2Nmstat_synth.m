% N2Nmstat_synth.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the script to carry out a bunch of N to N-m statisical analysis
% to deconvoluted catalog from the synthetics generated from different 
% source region sizes and saturation levels. The bulk is very similar to
% what has been done to real data in 'deconv_ref_4s_exp_4thsta.m'.
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/10/18
% Last modified date:   2023/10/18
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

singleflag = 0;

if ~singleflag
  %%%synthetics from different region sizes and saturation levels
  savefile = 'rst_decon_synth.mat';
  ttstr1 = {'Noise-free synthetics, '};
  load(savefile);
  nrounds = nreg;
else
  %%%synthetics from different noise and saturation levels, sources at a single spot
  savefile = 'rst_synth_onespot.mat';
  ttstr1 = {'Single-spot synthetics, '};
  load(savefile);
  nrounds = ntrial;
end


%%
%%%param for secondary sources removed
ttstr2 = 'Secondary sources removed';
fnsuffix = [];
impplt = imp;
denom = 18275;
% %%%param for further checked at KLNB
% ttstr2 = 'Further checked at KLNB';
% fnsuffix = '4th';
% impplt = imp4th;
% denom = 10547;

m=5;
nsep=1;
supertstr = strcat(ttstr1,ttstr2);

%% summarize the whole catalog, diff arrival time and fractions
f = initfig(15,8,2,round(nnsat/2)); %initialize fig
tit=supertit(f.ax,supertstr);
movev(tit,0.2);

color = jet(nrounds);
  
%%%loop for saturation level
for insat = 1: nnsat
  disp(nsat(insat));
  
  ax=f.ax(insat); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  binwdt = 0.05;
  dtcut = 0.25*nsep+0.125;
  xran = [0 2];
  if nsep == 1
    yran = [0 0.2];
  elseif nsep == 2
    yran = [0 0.06];
  end
  nx = round(xran(2)/binwdt)+1;
  edges = xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2;
  cnt = xran(1): binwdt: xran(2);
  
  %%%loop for region size or noise level
  for iround = 1: nrounds
    impi = impplt{insat,iround};
    nsrc = size(impi,1);
    
    %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
    dtarvl = srcdistNtoNm(impi(:,1),impi(:,7:8),m);
    %   dtarvlnnsepall = [dtarvlnnsepall; dtarvl{nsep}];
    if isempty(dtarvl{nsep})
      continue
    end
    dtarvlplt = dtarvl{nsep};
    
    N=histcounts(dtarvlplt/sps,edges,'normalization','count');
    Nn = N./length(dtarvlplt);
%     Nn = N./denom;
    frac(insat,iround) = sum(dtarvlplt/sps<=dtcut)/length(dtarvlplt);
    p(iround) = plot(ax,cnt,Nn,'color',color(iround,:),'LineWidth',1);
    if ~singleflag
      label{iround} = sprintf('a/2=%.2f,b/2=%.2f',semia(iround),semib(iround));
    else
      label{iround} = sprintf('noise=%.1f',perctrial(iround));
    end

  end
  %   text(ax,0.95,0.85,sprintf('Fraction w/i %.3f s: %.2f',dtcut,frac),'HorizontalAlignment','right',...
  %     'Units','normalized','FontSize',12);
  if insat == 1
    ylabel(ax,'Normalized count');
    xlabel(ax,sprintf('Diff. arrival between sources N and N-%d (s)',nsep));
    legend(ax,p,label);
  end
  text(ax,0.02,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','left');  
  xlim(ax,xran);
  ylim(ax,yran);
  longticks(ax,2);
  hold(ax,'off');
end

%%%summarize the whole catalog, diff arrival time and fractions
f = initfig(4,5,1,1); %initialize fig
tit=supertit(f.ax,supertstr);
movev(tit,0.2);
color = jet(nrounds);
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%%%loop for region size
for iround = 1: nrounds
  p(iround) = plot(ax,log10(nsat),frac(:,iround),'-o','markersize',4,'color',color(iround,:));
  if ~singleflag
    label{iround} = sprintf('a/2=%.2f,b/2=%.2f',semia(iround),semib(iround));
  else
    label{iround} = sprintf('noise=%.1f',perctrial(iround));
  end
end
legend(ax,p,label);
xlabel(ax,'log_{10}(Saturation)');
ylabel(ax,sprintf('Frac. of diff. arrival w/i %.3f s',dtcut));
yran = [0 1];
ylim(ax,yran);
xlim(ax,[-0.5 2]);

keyboard


%% bin by amp, diff time and frac, for each satur and size/noise level
%%%first bin source by amp, then plot diff arrival time for N and N-1 for each amp bin
%%%loop for saturation level
for insat = 1: nnsat
  disp(nsat(insat));
  %%%loop for region size or noise level
  for iround = 1: nrounds
    impi = impplt{insat,iround};
    nsrc = size(impi,1);
    %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
    dtarvl = srcdistNtoNm(impi(:,1),impi(:,7:8),m);
    %   dtarvlnnsepall = [dtarvlnnsepall; dtarvl{nsep}];
    if isempty(dtarvl{nsep})
      continue
    end
    dtarvlplt = dtarvl{nsep};
    
    ampcont = [];
    for i = 1: nsrc-nsep
      impcont = impi(i: i+nsep,:);
      ampcont(i,1) = median(mean(impcont(:,[2 4 6]),2));
    end
    ampplt = ampcont;
    
    %%%
    f = initfig(12,4,1,2); %initialize fig
    tit=supertit(f.ax,supertstr);
    movev(tit,0.2);
    ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    %%%%%%%%%% if bin by amp with a equal number
    nbin = 5;
    [ampbin,indbin,n] = binxeqnum(ampplt,nbin);
    color = jet(nbin);
    binwdt = 0.05;
    dtcut = 0.25*nsep+0.125;
    if nsep < 3
      xran = [0 2];
    else
      xran = [0 2*ceil(dtcut)];
    end
    nx = round(xran(2)/binwdt)+1;
    Nn = zeros(nx, nbin);
    edges = xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2;
    cnt = xran(1): binwdt: xran(2);
    label = []; p = []; fraci = [];
    for i = 1: nbin
      ind = indbin{i};
      dtarvli = dtarvlplt(ind);
      ampbincnt(i) = median(ampbin{i});
      N=histcounts(dtarvli/sps,edges,'normalization','count');
      Nn(:,i) = N./mode(n);
      fraci(i) = sum(dtarvli/sps<=dtcut)/mode(n);
      p(i)=plot(ax,cnt,Nn(:,i),'color',color(i,:),'LineWidth',1);
      label{i} = sprintf('amp of %.1f',ampbincnt(i));
      % label{i} = sprintf('%d/%dth amp',i,nbin);
      % keyboard
    end
    %%%%%%%%%% if bin by amp with a equal number
    Nnm = mean(Nn,2);
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
    
    ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    yran=[0 1];
    plot(ax,log(ampbincnt),fraci,'k-','linew',1,'marker','o','markersize',4,'markerfacec','k');
    xlabel(ax,'Median log_{10}{amp}');
    % ylabel(ax,'Median diff. arrival (s)');
    ylabel(ax,sprintf('Frac. of diff. arrival w/i %.3f s',dtcut));
    ylim(ax,yran);
    title(ax,'Data');
    
  end
  
end
keyboard
    
%% bin by amp, only summarize frac, for each satur and size/noise level
f = initfig(15,8,2,round(nnsat/2)); %initialize fig
tit=supertit(f.ax,strcat(supertstr, sprintf(', N & N-%d',nsep)));
movev(tit,0.2);

color = jet(nrounds);

%%%loop for saturation level
for insat = 1: nnsat
  disp(nsat(insat));
  
  ax=f.ax(insat); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  
  %%%loop for region size or noise level
  for iround = 1: nrounds
    impi = impplt{insat,iround};
    nsrc = size(impi,1);
    %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
    dtarvl = srcdistNtoNm(impi(:,1),impi(:,7:8),m);
    %   dtarvlnnsepall = [dtarvlnnsepall; dtarvl{nsep}];
    if isempty(dtarvl{nsep})
      continue
    end
    dtarvlplt = dtarvl{nsep};
    
    ampcont = [];
    for i = 1: nsrc-nsep
      impcont = impi(i: i+nsep,:);
      ampcont(i,1) = median(mean(impcont(:,[2 4 6]),2));
    end
    ampplt = ampcont;
    
    %%%
    %%%%%%%%%% if bin by amp with a equal number
    nbin = 5;
    [ampbin,indbin,n] = binxeqnum(ampplt,nbin);
    binwdt = 0.05;
    dtcut = 0.25*nsep+0.125;
    if nsep < 3
      xran = [0 2];
    else
      xran = [0 2*ceil(dtcut)];
    end
    nx = round(xran(2)/binwdt)+1;
    Nn = zeros(nx, nbin);
    edges = xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2;
    cnt = xran(1): binwdt: xran(2);
    fraci = [];
    for i = 1: nbin
      ind = indbin{i};
      dtarvli = dtarvlplt(ind);
      ampbincnt(i) = median(ampbin{i});
      N=histcounts(dtarvli/sps,edges,'normalization','count');
      Nn(:,i) = N./mode(n);
      fraci(i,iround) = sum(dtarvli/sps<=dtcut)/mode(n);
    end
    p(iround) = plot(ax,log(ampbincnt),fraci(:,iround),'-','Color',color(iround,:),'linew',1,...
      'marker','o','markersize',4,'markerfacec',color(iround,:));
    if ~singleflag
      label{iround} = sprintf('a/2=%.2f,b/2=%.2f',semia(iround),semib(iround));
    else
      label{iround} = sprintf('noise=%.1f',perctrial(iround));
    end
  end
  if insat == 1
    xlabel(ax,'Median log_{10}{amp}');
    % ylabel(ax,'Median diff. arrival (s)');
    ylabel(ax,sprintf('Frac. of diff. arrival w/i %.3f s',dtcut));
    legend(ax,p,label);
  end
  text(ax,0.02,0.05,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','left');
  yran=[0 1];
  ylim(ax,yran);
  longticks(ax,2);
  hold(ax,'off');
    
end

%% fraction of event pairs w/i a diff time cut, and fraction of all catalog
% nsep = 14;
for nsep = 1:m

dcut = 0.25*nsep+0.125;
nsrcsep = nsrc-nsep;
nsrcsep(nsrcsep<0) = 0; 
ndcutpair = zeros(size(trange,1), 1);
ndcutsrc = zeros(size(trange,1), 1);
ndcutsrc2 = zeros(size(trange,1), 1);
imppair = cell(size(trange,1), 1);
imppairuni = cell(size(trange,1), 1);
impcont = cell(size(trange,1), 1);
impcontuni = cell(size(trange,1), 1);
indcont = cell(size(trange,1), 1);
indcontuni = cell(size(trange,1), 1);
m = 15;

for i = 1: size(trange,1)
  if nsrc(i) == 0
    continue
  end
  ist = sum(nsrc(1:i-1))+1;
  ied = ist+nsrc(i)-1;
  impi = imp(ist:ied,:);
  %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
  dtarvl = srcdistNtoNm(impi(:,1),impi(:,7:8),m);
%   dtarvlnnsepall = [dtarvlnnsepall; dtarvl{nsep}];
  if isempty(dtarvl{nsep})
    continue
  end
  impbf = impi(1:end-nsep,:);
  impaf = impi(1+nsep:end,:);
  if ~isequal(size(impbf,1),length(dtarvl{nsep}))
    disp('Check');
  end
  ind = find(dtarvl{nsep}/sps <= dcut);
  ndcutpair(i,1) = length(ind);
  tmp = [];
  for j = 1: length(ind)
    tmp = [tmp; impbf(ind(j),:); impaf(ind(j),:)];
  end
  imppair{i,1} = tmp;
  % impsort = sortrows(imppair,1);
  imppairuni{i,1} = unique(tmp,'rows','stable');
  % impall = [impall; impuni];
  ndcutsrc(i,1) = size(imppairuni{i,1},1);

  tmp = [];
  for j = 1: length(ind)
    tmp = [tmp; impi(ind(j):ind(j)+nsep,:)];
  end
  impcont{i,1} = tmp;
  impcontuni{i,1} = unique(tmp,'rows','stable');
  % impall = [impall; impuni];
  ndcutsrc2(i,1) = size(impcontuni{i,1},1);

  tmp = [];
  for j = 1: length(ind)
    tmp = [tmp; reshape(ind(j):ind(j)+nsep, [], 1)];
  end
  indcont{i,1} = tmp;
  indcontuni{i,1} = unique(tmp,'rows','stable');
end

imppairunia = cat(1,imppairuni{:});
ndcutsrcall = sum(ndcutsrc);
fracsrcall = ndcutsrcall/sum(nsrc)*100;

impcontunia = cat(1,impcontuni{:});
ndcutsrc2all(nsep,1) = sum(ndcutsrc2);
fracsrc2all(nsep,1) = ndcutsrc2all(nsep,1)/sum(nsrc)*100;

ndcutpairall(nsep,1) = sum(ndcutpair);
fracpairall = ndcutpairall(nsep,1)/sum(nsrcsep)*100;

end

for nsep = 1:15
  dcut = 0.25*nsep+0.125;
  fprintf('%d clusters of %d consecutive events (%d/%d unique) w/i %.3f s \n',...
    ndcutpairall(nsep,1), nsep+1, ndcutsrc2all(nsep,1), sum(nsrc), dcut);
end

fracsrc2all = [100; fracsrc2all];
fprintf('%.3f \n',fracsrc2all);
dfracsrc2all = fracsrc2all(1:end-1) - fracsrc2all(2:end);
fprintf('%.3f \n',dfracsrc2all);

keyboard


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



% keyboard



%% 

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

