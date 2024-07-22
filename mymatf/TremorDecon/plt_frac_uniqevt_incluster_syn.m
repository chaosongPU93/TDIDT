% plt_frac_uniqevt_incluster_syn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to 'plt_frac_uniqevt_incluster', this is the script in particular to 
% plot the comparison of the fraction of 
% the catalog in terms of unique events in different type of clusters between
% 3-sta and 4-sta catalogs, from synthetics, either from noise-free synthetics
% (diff source region sizes and saturation level) or from single-spot 
% synthetics (diff noise and saturation levels).
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/04/04
% Last modified date:   2024/04/04
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
tlensum = sum(tlen);

modname = 'timeoff_plfit_4thsta_160sps.mat';
planefit = load(strcat(rstpath, '/MAPS/',modname));
rmse = planefit.gof{4}.rmse;
offmax = round(2.0*rmse);

%%%specify regime for transformation from time offset to map location
% ftrans = 'interpArmb';
% ftrans = 'interpArmbreloc';
ftrans = 'interpchao';

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
    
    imp{insat,iround} = allsyn.imp;
    nsrc{insat,iround} = allsyn.nsrcsum;
    imp4th{insat,iround} = allsyn.imp4th;
    nsrc4th{insat,iround} = allsyn.nsrc4thsum;
    
  end
end

% savefile ='rst_decon_synthmedwtcoef_td0.25.mat';
% load(strcat(workpath,'/synthetics/',savefile));
% imp=allsyn.impk;
% nsrc=allsyn.nsrcsum;
% imp4th=allsyn.imp4thk;
% nsrc4th=allsyn.nsrc4thsum;


% keyboard
%%
%%%param for secondary sources removed
supertstr = '3-station';
fnsuffix = [];

%%%param for further checked at KLNB
supertstr4th = '4-station';
fnsuffix4th = '4th';

%% compution of fractions for diff catalogs
for iround = 1: nround
  for insat = 1 : nnsat
    
    nsrci = nsrc{insat,iround};
    impi = imp{insat,iround};
    %%%% FOR 3-station catalog
    %%%determine the 'mmax' for which the resulting number of clusters is nonzero
    timetype = 'tarvl';
    nbst = 1;
    mmax=getmmaxcluster(nbst,impi,nsrci,sps,timetype);
    %%%for a cluster, not only the time separation between N and N-m needs to
    %%%be smaller than 'dtcut', but also the max time separation between each
    %%%consecutive events needs to smaller than 0.25+0.125 s.
    %%%ie, doublet means a cluster of 2 events ONLY occur as doublets
    [catclus,catclusbst,catimp,catuimp,catmedamp,catdtnnm]=...
      evtcluster_ex(nbst,impi,nsrci,mmax,sps,timetype);
    mmaxk{insat,iround} = mmax;
    
    %%%fraction of unique events ONLY occurring as certain clusters
    [fracuimp,nuimp,fracuimpsum]=frac_uniqevt_incluster2(catuimp,catclus,...
      impi,nsrci,mmax);
    fracuimpk{insat,iround} = fracuimp;
    fracuimpsumk{insat,iround} = fracuimpsum;
    
    %%%% FOR 4-station catalog
    nsrc4thi = nsrc4th{insat,iround};
    imp4thi = imp4th{insat,iround};
    mmax4th=getmmaxcluster(nbst,imp4thi,nsrc4thi,sps,timetype);
    [catclus4th,catclusbst4th,catimp4th,catuimp4th,catmedamp4th,catdtnnm4th]=...
      evtcluster_ex(nbst,imp4thi,nsrc4thi,mmax4th,sps,timetype);
    [fracuimp4th,~,fracuimp4thsum]=frac_uniqevt_incluster2(catuimp4th,catclus4th,...
      imp4thi,nsrc4thi,mmax4th);
    mmax4thk{insat,iround} = mmax4th;
    fracuimp4thk{insat,iround} = fracuimp4th;
    fracuimp4thsumk{insat,iround} = fracuimp4thsum;

  end
end

%% compution of fractions for data
%%%load data
% savefile = 'deconv_stats4th_allbstsig.mat';
savefile = 'deconv_stats4th_no23_allbstsig.mat';
load(strcat(rstpath, '/MAPS/',savefile));
% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
load(strcat(rstpath, '/MAPS/',savefile));
  
nsrcd = allbstsig.nsrc;
impd = allbstsig.impindepall;
nsrcn = allbstnoi.nsrc;
impn = allbstnoi.impindepall;
nsrc4thd = allbstsig.nsrc4th;
imp4thd = allbstsig.impindep4thall;
nsrc4thn = allbstnoi.nsrc4th;
imp4thn = allbstnoi.impindep4thall;

cutout = 'ellipse';
ttol = 35;
ntol = 3;
trange = load(strcat(rstpath,'/MAPS/tdec.bstran',num2str(ttol),'s.pgc002.',cutout(1:4)));
% trange = trange(1:end-1,:);
nbst = size(trange,1);

%%%determine the 'mmax' for which the resulting number of clusters is nonzero
timetype = 'tarvl';
mmaxd=getmmaxcluster(nbst,impd,nsrcd,sps,timetype);
%%%for a cluster, not only the time separation between N and N-m needs to
%%%be smaller than 'dtcut', but also the max time separation between each
%%%consecutive events needs to smaller than 0.25+0.125 s. 
%%%ie, doublet means a cluster of 2 events ONLY occur as doublets
[catclus,catclusbst,catimp,catuimp,catmedamp,catdtnnm]=...
  evtcluster_ex(nbst,impd,nsrcd,mmaxd,sps,timetype);

%%%fraction of unique events ONLY occurring as certain clusters
[fracuimpd,nuimpd,fracuimpsumd]=frac_uniqevt_incluster2(catuimp,catclus,impd,nsrcd,mmaxd);

%%%% FOR DATA, 4-station catalog
mmax4thd=getmmaxcluster(nbst,imp4thd,nsrc4thd,sps,timetype);
[catclus4th,catclusbst4th,catimp4th,catuimp4th,catmedamp4th,catdtnnm4th]=...
  evtcluster_ex(nbst,imp4thd,nsrc4thd,mmax4thd,sps,timetype);
[fracuimp4thd,~,fracuimp4thsumd]=frac_uniqevt_incluster2(catuimp4th,catclus4th,imp4thd,nsrc4thd,mmax4thd);

if singleflag
  %%%% FOR NOISE, 3-station catalog
  mmaxn=getmmaxcluster(nbst,impn,nsrcn,sps,timetype);
  [catclusn,catclusbstn,catimpn,catuimpn,catmedampn,catdtnnmn]=...
    evtcluster_ex(nbst,impn,nsrcn,mmaxn,sps,timetype);
  [fracuimpn,~,fracuimpsumn]=frac_uniqevt_incluster2(catuimpn,catclusn,impn,nsrcn,mmaxn);

  %%%% FOR NOISE, 4-station catalog
  mmax4thn=getmmaxcluster(nbst,imp4thn,nsrc4thn,sps,timetype);
  [catclus4thn,catclusbst4thn,catimp4thn,catuimp4thn,catmedamp4thn,catdtnnm4thn]=...
    evtcluster_ex(nbst,imp4thn,nsrc4thn,mmax4thn,sps,timetype);
  [fracuimp4thn,~,fracuimp4thsumn]=frac_uniqevt_incluster2(catuimp4thn,catclus4thn,imp4thn,nsrc4thn,mmax4thn);
end

%%
%%%PLOT
% nrow = 1; % rows and cols of subplots in each figure
% ncol = 2; 
% widin = 7; % size of each figure
% htin = 3.5;
% pltxran = [0.1 0.98]; pltyran = [0.15 0.95];
% pltxsep = 0.07; pltysep = 0.03;
% f = initfig(widin,htin,nrow,ncol);
% axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

nrow = 2; ncol = 4;
widin = 8.4;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);
pltxran = [0.08 0.98]; pltyran = [0.08 0.98]; % optimal axis location
pltxsep = 0.01; pltysep = 0.02;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

color = gradientblue(nround);
p=[]; label=[];
for insat = 1: nnsat
  
  ax=f.ax(insat); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.YAxis.Scale = 'log';
  
  for iround = 1: nround
    p(iround)=plot(ax,1:mmaxk{insat,iround}+1,fracuimpsumk{insat,iround},'-',...
      'linew',1,'color',color(iround,:),'marker','o','markersize',3,...
      'markerfacec',color(iround,:));
    if ~singleflag
      label{iround} = sprintf('%.1fx%.1f',2*semia(iround),2*semib(iround));
    else
      label{iround} = sprintf('%.1f',perctrial(iround));
    end
  end
  p(nround+1) = plot(ax,1:mmaxd+1,fracuimpsumd,'k--','Linew',1.5,'marker','o',...
    'markersize',3,'markerfacec','k');
  label{nround+1} = 'Data';
  if singleflag
    p(nround+2) = plot(ax,1:mmaxn+1,fracuimpsumn,'r--','Linew',1.5,'marker','o',...
      'markersize',3,'markerfacec','r');
    label{nround+2} = 'Noise';
  end
  text(ax,0.98,0.05,sprintf('Sat=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','right');
  if insat == (nrow-1)*ncol+1 
    lgd=legend(ax,p,label,'NumColumns',2,'Location','best','fontsize',6); %,'Orientation','horizontal'
    set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
    if ~singleflag
      lgdtit = 'Region size (km)';
    else
      lgdtit = 'Noise level';
    end
    title(lgd,lgdtit,'fontsize',7);
    % xlabel(ax,'# of events in cluster','FontSize',10);
    % ylabel(ax,'% of catalog in such clusters','FontSize',10);
    xlabel(ax,'m (# of events in cluster)','FontSize',10);
    ylabel(ax,'% of catalog in clusters of >=m events','FontSize',10);
  else
    nolabels(ax,3);
  end
%   xlim(ax,[1 max(cat(1,mmaxk{:}))+2]);
%   xticks(ax,0:2:max(cat(1,mmaxk{:}))+2);
  xlim(ax,[1 mmaxd+2]);
  ylim(ax,[5e-2 1e2]);
  yticks(ax,[0.1 0.2 0.5 1 2 5 10 20 50 100]);
  text(ax,0.98,0.95,supertstr,'HorizontalAlignment','right','Units','normalized',...
    'fontsize',10);
%   text(ax,0.02,0.05,'a','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
end
% supertit(f.ax,ttstr1);
fname = strcat('fracunievtclus_syn',fnsuffix1,'.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/chaosong/Pictures/',fname));

%%
nrow = 2; ncol = 4;
widin = 8.4;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);
pltxran = [0.08 0.98]; pltyran = [0.08 0.98]; % optimal axis location
pltxsep = 0.01; pltysep = 0.02;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

color = gradientblue(nround);
p=[]; label=[];
for insat = 1: nnsat
  
  ax=f.ax(insat); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.YAxis.Scale = 'log';
  
  for iround = 1: nround
    p(iround)=plot(ax,1:mmax4thk{insat,iround}+1,fracuimp4thsumk{insat,iround},'-',...
      'linew',1,'color',color(iround,:),'marker','o','markersize',3,...
      'markerfacec',color(iround,:));
    if ~singleflag
      label{iround} = sprintf('%.1fx%.1f',2*semia(iround),2*semib(iround));
    else
      label{iround} = sprintf('%.1f',perctrial(iround));
    end
  end
  p(nround+1) = plot(ax,1:mmax4thd+1,fracuimp4thsumd,'k--','Linew',1.5,'marker',...
    'o','markersize',3,'markerfacec','k');
  label{nround+1} = 'Data';
  if singleflag
    p(nround+2) = plot(ax,1:mmax4thn+1,fracuimp4thsumn,'r--','Linew',1.5,'marker',...
      'o','markersize',3,'markerfacec','r');
    label{nround+2} = 'Noise';
  end
  text(ax,0.98,0.05,sprintf('Sat=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','right');
  if insat == (nrow-1)*ncol+1 
    lgd=legend(ax,p,label,'NumColumns',2,'Location','best','fontsize',6); %,'Orientation','horizontal'
    set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
    if ~singleflag
      lgdtit = 'Region size (km)';
    else
      lgdtit = 'Noise level';
    end
    title(lgd,lgdtit,'fontsize',7);
    % xlabel(ax,'# of events in cluster','FontSize',10);
    % ylabel(ax,'% of catalog in such clusters','FontSize',10);
    xlabel(ax,'m (# of events in cluster)','FontSize',10);
    ylabel(ax,'% of catalog in clusters of >=m events','FontSize',10);
  else
    nolabels(ax,3);
  end
%   xlim(ax,[1 max(cat(1,mmaxk{:}))+2]);
%   xticks(ax,0:2:max(cat(1,mmaxk{:}))+2);
  xlim(ax,[1 mmaxd+2]);
  ylim(ax,[5e-2 1e2]);
  yticks(ax,[0.1 0.2 0.5 1 2 5 10 20 50 100]);
  text(ax,0.98,0.95,supertstr4th,'HorizontalAlignment','right','Units','normalized',...
    'fontsize',10);
%   text(ax,0.02,0.05,'a','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
end
% supertit(f.ax,ttstr1);

fname = strcat('fracunievtclus_syn4th',fnsuffix1,'.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/chaosong/Pictures/',fname));



