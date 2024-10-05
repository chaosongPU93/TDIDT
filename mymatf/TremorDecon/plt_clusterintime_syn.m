% plt_clusterintime_syn.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is similar to 'plt_clusterintime.m' which is for the data 
% catalog, to obtain the clusting in time of the deconvolved events.
% --synthetics have a smaller amp variation than data, we know that. So rather
% than getting the results after the binning by amp (either envelope amp, or
% impulse amp), we directly plot the results with amp binning using different
% synthetics.
% --The results are going to be compared with the similar version for data, ie.,
% Figure 13 in the paper. It can also be compared with fig. 16, the amp-binned
% version. 
% --for an easier comparison, data is plotted on top.
%  
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/09/17
% Last modified date:   2024/09/17
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

sps = 160;

ftrans = 'interpchao';

%%%load data
% savefile = 'deconv_stats4th_allbstsig.mat';
savefile = 'deconv_stats4th_no23_allbstsig.mat';
% savefile = 'deconv_stats4th_no23_allbstsig0.25s.mat';
load(strcat(rstpath, '/MAPS/',savefile));

% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
% savefile = 'deconv_stats4th_no23_allbstnoi0.25s.mat';
load(strcat(rstpath, '/MAPS/',savefile));

%% load decon results
tdura = 0.25;  %must be consistent with what synthetics actually used
% tdura = 0.4;

%times of saturation
% nsat=[0.1 0.4 1 2 4 10 20 40 100];
nsat=[0.4 1 2 4 10 20 40 100];
nnsat = length(nsat);

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
  
else  %%%synthetics from different noise and saturation levels, sources at a single spot
  %different percent of noise
  perctrial = 0.1*(0:2:16)';
  ntrial = length(perctrial);
  nround = ntrial;
  
  ttstr1 = {'Single-spot syn, '};
end

for iround = 1: nround
  for insat = 1: nnsat
    sat = nsat(insat);
    
    if ~singleflag
      xaxis = semia(iround);
      yaxis = semib(iround);
      savefile = strcat('rst_decon_synth_reg',num2str(xaxis),...
        '-',num2str(yaxis),'_nsat',num2str(sat),'_td',num2str(tdura),'.mat');
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
% keyboard

%%
%%%param for secondary sources removed
impuse = imp;
nsrcuse = nsrc;
supertstr = '3-station';
fnsuffix = [];
nsrcd = allbstsig.nsrc;
impd = allbstsig.impindepall;
nsrcn = allbstnoi.nsrc;
impn = allbstnoi.impindepall;

% %%%param for further checked at KLNB
% impuse = imp4th;
% nsrcuse = nsrc4th;
% supertstr = '4-station';
% fnsuffix = '4th';
% nsrcd = allbstsig.nsrc4th;
% impd = allbstsig.impindep4thall;
% nsrcn = allbstnoi.nsrc4th;
% impn = allbstnoi.impindep4thall;

%% computation, the probability of inter-event time > t separately
for iround = 1: nround
  for insat = 1 : nnsat
    nsrcs = nsrcuse{insat,iround};
    imps = impuse{insat,iround};
    
    %inter-event time, ie. time from each to its preceeding event, N to N-1
    m = 1;
    dtcut = 0.25*m+0.125;
    nbst = 1;
    [amp,dtinter]=med_amp_incluster(nbst,imps,nsrcs,m);
    %     [amp4th,dtinter4th]=med_amp_incluster(nbst,impuse4th,nsrcuse4th,m);
    
    dt=sort(dtinter(:,1)/sps);
    prob = prob_geq(dt);  %probability of inter-event time > t0
    
    probk{insat,iround} = prob;
  end
end

%%%% FOR DATA
nbst = size(trange,1);
[ampd,dtinterd]=med_amp_incluster(nbst,impd,nsrcd,m);
dt=sort(dtinterd(:,1)/sps);
probd = prob_geq(dt);  %probability of inter-event time > t0

if singleflag
  %%%% FOR NOISE
  [ampn,dtintern]=med_amp_incluster(nbst,impn,nsrcn,m);
  dt=sort(dtintern(:,1)/sps);
  probn = prob_geq(dt);  %probability of inter-event time > t0
end

%% plot,
nrow = 2; ncol = 4;
widin = 8.4;  % maximum width allowed is 8.5 inches
htin = 5 ;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);
pltxran = [0.08 0.98]; pltyran = [0.08 0.98]; % optimal axis location
pltxsep = 0.01; pltysep = 0.02;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

color = gradientblue(nround);

p=[]; label=[];
for insat = 1 : nnsat
  
  ax=f.ax(insat); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); set(ax,'YScale','log');
  
  for iround = 1: nround
    prob = probk{insat,iround};
    p(iround)=plot(ax,prob(:,1),prob(:,2),'-','linew',1,'color',color(iround,:));
%     fitxran = [2 xran(2)]; 
%     ind = find(bb(:,1)>=fitxran(1) & bb(:,1)<=fitxran(2));
%     fitstruct=robustexpfit(bb(ind,1),bb(ind,2),'log10',[1e3 -1]);
%     fitobj=fitstruct.fitobj;
%     coef=coeffvalues(fitobj); slpd=coef(2); intcptd=coef(1);
%     xfit = xran(1):1e-2:xran(2);
%     yfitd = feval(fitobj,xfit);
%     plot(ax,xfit,yfitd,'k--','LineWidth',1);

    if ~singleflag
      label{iround} = sprintf('%.2f, %.2f',semia(iround),semib(iround));
    else
      label{iround} = sprintf('Noi=%.1f',perctrial(iround));
    end
    
  end
  p(nround+1)=plot(ax,probd(:,1),probd(:,2),'-','linew',1.5,'color','k');
  label{nround+1} = 'Data';
  if singleflag
    p(nround+2)=plot(ax,probn(:,1),probn(:,2),'-','linew',1.5,'color','r');
    label{nround+2} = 'Noise';
  end
  text(ax,0.99,0.95,supertstr,'HorizontalAlignment','right','Units','normalized',...
    'fontsize',10);
  text(ax,0.99,0.85,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','right');
  if insat == (nrow-1)*ncol+1
    lgd=legend(ax,p,label,'NumColumns',2,'Location','southwest',...
      'fontsize',6); %,'Orientation','horizontal'
    set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
    if ~singleflag
      lgdtit = 'Region size (km)';
    else
      lgdtit = 'Noise level';
    end
    title(lgd,lgdtit,'fontsize',7);
    xlabel(ax,'Time (s) from each to its preceding');
    ylabel(ax,'Probability');
  elseif rem(insat,ncol)==1
    ylabel(ax,'Probability');
    nolabels(ax,1);
  elseif insat > (nrow-1)*ncol
    nolabels(ax,2);
  else
    nolabels(ax,3);
  end  
  
  %   xran = [0 20];
  %   yran = [1e-4 1];
  xran = [0 7];
  yran = [5e-4 1];
  xlim(ax,xran);
  ylim(ax,yran);
  longticks(ax,2);
end

%% compution, fraction of catalog of events in clusters vs. # of events in clusters
for iround = 1: nround
  for insat = 1 : nnsat
    nsrcs = nsrcuse{insat,iround};
    imps = impuse{insat,iround};
    %%%% FOR 3-station catalog
    %%%determine the 'mmax' for which the resulting number of clusters is nonzero
    timetype = 'tarvl';
    nbst = 1;
    mmax=getmmaxcluster(nbst,imps,nsrcs,sps,timetype);
    %%%for a cluster, not only the time separation between N and N-m needs to
    %%%be smaller than 'dtcut', but also the max time separation between each
    %%%consecutive events needs to smaller than 0.25+0.125 s.
    %%%ie, doublet means a cluster of 2 events ONLY occur as doublets
    [catclus,catclusbst,catimp,catuimp,catmedamp,catdtnnm]=...
      evtcluster_ex(1,imps,nsrcs,mmax,sps,timetype);
    mmaxk{insat,iround} = mmax;
    
    %%%fraction of unique events ONLY occurring as certain clusters
    [fracuimp,nuimp,fracuimpsum]=frac_uniqevt_incluster2(catuimp,catclus,imps,...
      nsrcs,mmax);
    fracuimpk{insat,iround} = fracuimp;
    fracuimpsumk{insat,iround} = fracuimpsum;
    
  end
end

%%%% FOR DATA
nbst = size(trange,1);
mmaxd=getmmaxcluster(nbst,impd,nsrcd,sps,timetype);
[catclusd,catclusbstd,catimpd,catuimpd,catmedampd,catdtnnmd]=...
  evtcluster_ex(nbst,impd,nsrcd,mmaxd,sps,timetype);
[fracuimpd,nuimpd,fracuimpsumd]=frac_uniqevt_incluster2(catuimpd,catclusd,impd,nsrcd,mmaxd);

if singleflag
  %%%% FOR NOISE
  mmaxn=getmmaxcluster(nbst,impn,nsrcn,sps,timetype);
  [catclusn,catclusbstn,catimpn,catuimpn,catmedampn,catdtnnmn]=...
    evtcluster_ex(nbst,impn,nsrcn,mmaxn,sps,timetype);
  [fracuimpn,~,fracuimpsumn]=frac_uniqevt_incluster2(catuimpn,catclusn,impn,nsrcn,mmaxn);
end

%% plot,
nrow = 2; ncol = 4;
widin = 8.4;  % maximum width allowed is 8.5 inches
htin = 5 ;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);
pltxran = [0.08 0.98]; pltyran = [0.08 0.98]; % optimal axis location
pltxsep = 0.01; pltysep = 0.02;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

color = gradientblue(nround);

p=[]; label=[];
for insat = 1 : nnsat
  
  ax=f.ax(insat); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); set(ax,'YScale','log');
  
  for iround = 1: nround
    
    mmax = mmaxk{insat,iround};
    fracuimpsum = fracuimpsumk{insat,iround};
    
    p(iround)=plot(ax,1:mmax+1,fracuimpsum,'-',...
      'linew',1,'color',color(iround,:),'marker','o','markersize',3,...
      'markerfacec',color(iround,:));
    if ~singleflag
      label{iround} = sprintf('%.2f, %.2f',semia(iround),semib(iround));
    else
      label{iround} = sprintf('Noi=%.1f',perctrial(iround));
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
  text(ax,0.99,0.95,supertstr,'HorizontalAlignment','right','Units','normalized',...
    'fontsize',10);
  text(ax,0.99,0.85,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','right');
  if insat == (nrow-1)*ncol+1 
    lgd=legend(ax,p,label,'NumColumns',2,'Location','southwest',...
      'fontsize',6); %,'Orientation','horizontal'
    set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
    if ~singleflag
      lgdtit = 'Region size (km)';
    else
      lgdtit = 'Noise level';
    end
    title(lgd,lgdtit,'fontsize',7);
    xlabel(ax,'m (# of events in cluster)','FontSize',10);
    ylabel(ax,'% of catalog in clusters of >=m events','FontSize',10);
  elseif rem(insat,ncol)==1
%     ylabel(ax,'% of catalog in clusters of >=m events','FontSize',10);
    nolabels(ax,1);
  elseif insat > (nrow-1)*ncol
    nolabels(ax,2);
  else
    nolabels(ax,3);
  end

  xlim(ax,[1 mmaxd+2]);
  ylim(ax,[5e-2 1e2]);
  yticks(ax,[0.1 0.2 0.5 1 2 5 10 20 50 100]);

end








