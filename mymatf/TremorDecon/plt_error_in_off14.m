% plt_error_in_off14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the script in particular to plot the comparison of the fraction of 
% the catalog in terms of unique events in different type of clusters between
% 3-sta data, 3-sta noise, 4-sta data and 4-sta noise.
% --2024/05/07, add the result from synthetics to this code as well.
% 
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/04/04
% Last modified date:   2024/05/07
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

%%
%%%load data
% savefile = 'deconv_stats4th_allbstsig.mat';
savefile = 'deconv_stats4th_no23_allbstsig.mat';
allsig = load(strcat(rstpath, '/MAPS/',savefile));
% savefile = 'deconv1win_stats4th_allbstsig.mat';
savefile = 'deconv1win_stats4th_no23_allbstsig.mat';
allsig1win = load(strcat(rstpath, '/MAPS/',savefile));

%first is the difference (error) in the found arrival peak from decon and
%predicted arrival peak from plane fit
pred4diff = allsig.allbstsig.pred4difftrall;
diffoff14 = allsig.allbstsig.diffoff14trall;

%second is the difference (error) in the off14 from arrival peaks from decon and
%predicted off14 from plane fit (note that off14 is the direct product of plane
%fit)
pred4diff1win = allsig1win.allbstsig.pred4difftrall;
diffoff141win = allsig1win.allbstsig.diffoff14trall;

%%%load results from noise
% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
allnoi = load(strcat(rstpath, '/MAPS/',savefile));
diffoff14n = allnoi.allbstnoi.diffoff14trall;


%% load decon results
tdura = 0.25;  %must be consistent with what synthetics actually used
% tdura = 0.4;

%times of saturation
% nsat=[0.1 0.4 1 2 4 10 20 40 100];
nsat=[0.4 1 2 4 10 20 40 100];
nnsat = length(nsat);

sps = 160;

%%%flag to decide which type of synthetics to use
singleflag = 0;
% singleflag = 1;

if ~singleflag  %%%synthetics from different region sizes and saturation levels
  %%%specify if considering the physical size of each source
  % physicalsize = 1;
  physicalsize = 0;
  
  %%%diameter of physical size
  if physicalsize
    diam=0.15;% 0.5; %0.6; %
  else
    diam=0;
  end
  
  %%%specify shape of the source region
  srcregion='ellipse';
  % srcregion='rectangle';
  % srcregion='circle';
  
  %variation of source region size
  if strcmp(srcregion,'ellipse')
    semia = 1.75*(0.6:0.2:2.0);
    semib = 1.25*(0.6:0.2:2.0);
    nreg = length(semia);
  end
  nround = nreg;
  
  ttstr1 = {'Noise-free syn, '};
  fnsuffix1 = '';
  
else  %%%synthetics from different noise and saturation levels, sources at a single spot
  %different percent of noise
  perctrial = 0.1*(0:2:16)';
  ntrial = length(perctrial);
  nround = ntrial;
  
  ttstr1 = {'Single-spot syn, '};
  fnsuffix1 = '_onespot';
end

for iround = 1: nround
  for insat = 1: nnsat
    sat = nsat(insat);
    
    if ~singleflag
      xaxis = semia(iround);
      yaxis = semib(iround);
      savefile = strcat('rst_decon_synth_reg',num2str(xaxis),...
        '-',num2str(yaxis),'_nsat',num2str(sat),'_td',num2str(tdura),'.mat'); %,srcregion(1:3),'_'
    else
      perc = perctrial(iround);
      savefile = strcat('rst_decon_synth_onespot','_noi',num2str(perc),'_nsat',...
        num2str(sat),'_td',num2str(tdura),'.mat');
    end
    load(strcat(workpath,'/synthetics/',savefile));
      
    diffoff14syn{insat,iround} = allsyn.diffoff14;

  end
end

% savefile ='rst_decon_synthmedwtcoef_td0.25.mat';
% load(strcat(workpath,'/synthetics/',savefile));
% diffoff14syn=allsyn.diffoff14k;

%%
%%%plot the difference (error) in the off14 from arrival peaks from decon 
%%%andpredicted off14 from plane fit (note that off14 is the direct product of 
%%%planefit), MORE IMPORTANT!
f = initfig(4,4.5,1,1); %initialize fig
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
binw=1;
binedge=(-offmax-0.5: binw: offmax+0.5)/sps;
% histogram(ax,diffoff14/sps,'Normalization','probability','BinEdges',binedge,...
%     'Facec','k','edgec','none');
[N]=histcounts(diffoff14/sps,'BinEdges',binedge,'normalization','probability');
N=[N N(end)];
p(1)=stairs(ax,binedge,N,'-','linew',1.5,'color','k');
label{1} = 'Data';
% plot(ax,[-offmax -offmax]/sps,ax.YLim,'k--');
% plot(ax,[offmax offmax]/sps,ax.YLim,'k--');
% aa=diffoff14;
% bb=aa(aa<=8);
% mean(bb)
% median(bb)
% histogram(ax,bb/sps,'Normalization','pdf','BinEdges',binedge,...
%     'Facec','b','edgec','none');
[mu,sigma]=normfit(diffoff14/sps);
differror = round(mu*sps);
x=(-offmax:0.2:offmax)/sps; yfit=normpdf(x,mu,sigma)*binw/sps;
% plot(ax,x,yfit,'k-','LineWidth',2);
% text(ax,0.99,0.95,'Short-win','HorizontalAlignment','right','Units',...
%   'normalized','FontSize',10);  
% text(ax,0.99,0.95,'Data','HorizontalAlignment','right','Units',...
%   'normalized','FontSize',10);  

%%%load the result from one combo of syntheics
%let's say choose the same size of black ellipse, a saturation of 4
%%%possible reg size could be 4th or 5th; possible noi level 7th or 8th
insat = 4;
iround = 6;
nsat(insat)
xaxis = semia(iround)
yaxis = semib(iround)
diffoff14trsyn = diffoff14syn{insat,iround};
% histogram(ax,diffoff14tr1win);
% histogram(ax,diffoff14trsyn/sps,'Normalization','probability','BinEdges',binedge,...
%     'Facec','k','edgec','none');
[Ns]=histcounts(diffoff14trsyn/sps,'BinEdges',binedge,'normalization','probability');
Ns=[Ns Ns(end)];
p(2)=stairs(ax,binedge,Ns,'-','linew',1.5,'color','r');
% label{2} = sprintf('Syn, %.1fx%.1f km, %d',2*semia(iround),...
%   2*semib(iround),nsat(insat));
label{2} = 'Synthetics';
plot(ax,[median(diffoff14) median(diffoff14)]/sps,ax.YLim,'k--','LineWidth',1);
plot(ax,[median(diffoff14trsyn) median(diffoff14trsyn)]/sps,ax.YLim,'r--','LineWidth',1);
% plot(ax,[-offmax -offmax]/sps,ax.YLim,'k--');
% plot(ax,[offmax offmax]/sps,ax.YLim,'k--');
% aa=diffoff141win;
% bb=aa(aa<=9);
% mean(bb)
% median(bb)
% histogram(ax,bb/sps,'Normalization','pdf','BinEdges',binedge,...
%     'Facec','b','edgec','none');
[mus,sigmas]=normfit(diffoff14trsyn/sps);
differrors = round(mus*sps);
yfit=normpdf(x,mus,sigmas)*binw/sps;
% plot(ax,x,yfit,'k-','LineWidth',2);
% xlabel(ax,'diff. in off14 between plane-fit and decon');
% ylabel(ax,'count');
% text(ax,0.99,0.95,'Synthetics','HorizontalAlignment','right','Units',...
%   'normalized','FontSize',10); 
legend(ax,p,label,'Location','south');
xlabel(ax,sprintf('\\Delta{t}_{12} (s)'));
xlabel(ax,sprintf('Diff. in \\Delta{t}_{14}^{abs} between deconvolution and plane-fit'));
ylabel(ax,'Probability');

  
fname = strcat('diffoff14_demo',fnsuffix1,'.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

% keyboard

% ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% % histogram(ax,diffoff14tr1win);
% histogram(ax,diffoff141win/sps,'Normalization','probability','BinEdges',binedge,...
%     'Facec','k','edgec','none');
% plot(ax,[median(diffoff141win) median(diffoff141win)]/sps,ax.YLim,'r--','LineWidth',1.5);
% plot(ax,[-offmax -offmax]/sps,ax.YLim,'k--');
% plot(ax,[offmax offmax]/sps,ax.YLim,'k--');
% % aa=diffoff141win;
% % bb=aa(aa<=9);
% % mean(bb)
% % median(bb)
% % histogram(ax,bb/sps,'Normalization','pdf','BinEdges',binedge,...
% %     'Facec','b','edgec','none');
% [mu1win,sigma1win]=normfit(diffoff141win/sps);
% differror1win = round(mu1win*sps);
% yfit=normpdf(x,mu1win,sigma1win)*binw/sps;
% plot(ax,x,yfit,'k-','LineWidth',2);
% % xlabel(ax,'diff. in off14 between plane-fit and decon');
% % ylabel(ax,'count');
% text(ax,0.99,0.95,'Whole-win','HorizontalAlignment','right','Units',...
%   'normalized','FontSize',10);  
% keyboard

%%
% %%%plot the difference (error) in the found arrival peak from decon and
% %%%predicted arrival peak from plane fit
% f = initfig(12,5,1,2); %initialize fig
% ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% histogram(ax,pred4difftr/sps,'Normalization','count','BinEdges',binedge,...
%     'Facec','k','edgec','none');
% plot(ax,[median(pred4difftr) median(pred4difftr)]/sps,ax.YLim,'r--','LineWidth',1);
% plot(ax,[-offmax -offmax]/sps,ax.YLim,'k--');
% plot(ax,[offmax offmax]/sps,ax.YLim,'k--');
% xlabel(ax,'diff. in 4th arrival between pred and decon');
% ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
% binw=1;
% binedge=(-offmax-0.5: binw: offmax+0.5)/sps;
% histogram(ax,pred4difftr1win/sps,'Normalization','count','BinEdges',binedge,...
%     'Facec','k','edgec','none');
% plot(ax,[median(pred4difftr1win) median(pred4difftr1win)]/sps,ax.YLim,'r--','LineWidth',1);
% plot(ax,[-offmax -offmax]/sps,ax.YLim,'k--');
% plot(ax,[offmax offmax]/sps,ax.YLim,'k--');
% xlabel(ax,'diff. in 4th arrival between pred and decon');


%%
%%%plot the difference (error) in the off14 from arrival peaks from decon 
%%%andpredicted off14 from plane fit (note that off14 is the direct product of 
%%%planefit), MORE IMPORTANT!
nrow = 2; ncol = 4;
widin = 8.4;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);
pltxran = [0.08 0.98]; pltyran = [0.08 0.98]; % optimal axis location
pltxsep = 0.01; pltysep = 0.02;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

color = gradientblue(nround);
for insat = 1 : nnsat
  ax=f.ax(insat); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  xlim(ax,[-0.1 0.1]);  
  ylim(ax,[0 0.17]);
  
  binw=1;
  binedge=(-offmax-0.5: binw: offmax+0.5)/sps;
  
%   histogram(ax,diffoff14syn/sps,'Normalization','probability','BinEdges',binedge,...
%       'Facec','k','edgec','none');
  for iround = 1: nround
    aa = diffoff14syn{insat,iround}/sps;
    [N]=histcounts(aa,'BinEdges',binedge,'normalization','probability');
    N=[N N(end)];
    p(iround)=stairs(ax,binedge,N,'-','linew',1,'color',color(iround,:));
    if ~singleflag
      label{iround} = sprintf('%.1fx%.1f',2*semia(iround),2*semib(iround));
    else
      label{iround} = sprintf('%.1f',perctrial(iround));
    end
  end
  aa = diffoff14/sps;
  [N]=histcounts(aa,'BinEdges',binedge,'normalization','probability');
  N=[N N(end)];
  p(nround+1)=stairs(ax,binedge,N,'-','linew',1.5,'color','k');
  label{nround+1} = 'Data';  
  if singleflag
    aa = diffoff14n/sps;
    [N]=histcounts(aa,'BinEdges',binedge,'normalization','probability');
    N=[N N(end)];
    p(nround+2)=stairs(ax,binedge,N,'-','linew',1.5,'color','r');
    label{nround+2} = 'Noise';  
  end
%   plot(ax,[median(aa) median(aa)],ax.YLim,'r--','LineWidth',1.5);
%   plot(ax,[-offmax -offmax]/sps,ax.YLim,'k--');
%   plot(ax,[offmax offmax]/sps,ax.YLim,'k--');
  [mu,sigma]=normfit(aa);
  differror = round(mu*sps);
%   x=(-offmax:0.2:offmax)/sps; yfit=normpdf(x,mu,sigma)*binw/sps;
%   plot(ax,x,yfit,'k-','LineWidth',2);
%   text(ax,0.99,0.95,'Noise-free synthetics','HorizontalAlignment','right','Units',...
%     'normalized','FontSize',10);  
  text(ax,0.99,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','right');

  if insat == (nrow-1)*ncol+1 
    lgd=legend(ax,p,label,'NumColumns',2,'Location','best','fontsize',6);  %'Orientation','vertical'
    set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
    if ~singleflag
      lgdtit = 'Region size (km)';
    else
      lgdtit = 'Noise level';
    end
    title(lgd,lgdtit,'fontsize',7);
    xlabel(ax,'Diff. in off14 between plane-fit and decon');
    ylabel(ax,'Probability');
  else
    nolabels(ax,3);
  end

end

fname = strcat('diffoff14_syn',fnsuffix1,'.pdf');
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));


