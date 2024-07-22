% lfeinterevttime_syn.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to 'lfeinterevttime', this cript to obtain the inter-event time 
% of the deconvolved LFE catalog, either 3-station or 4-station catalog, from
% synthetic seismograms. Right now, we use single-spot synthetics with different
% nosie level tp mimic the different waveform amplitude in the case of real
% data.
%
% --The whole analysis is another way to address whether is the saturation 
% is related to the amplitude. In 'N2Nmstat_data_ref.m', we obtain the 
% fraction of diff time that is within 0.25*m+0.125 s, the fraction of 
% unique events, after binning these events by their median amp. Here, the 
% direct waveform amp is used. 

% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/05/07
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

%%%load real data
% savefile = 'deconv_stats4th_allbstsig.mat';
savefile = 'deconv_stats4th_no23_allbstsig.mat';
load(strcat(rstpath, '/MAPS/',savefile));
nsrcd = allbstsig.nsrc;
impd = allbstsig.impindepall;
nsrc4thd = allbstsig.nsrc4th;
imp4thd = allbstsig.impindep4thall;

%%%load synthetic noise
% savefile = 'deconv_stats4th_allbstnoi.mat';
savefile = 'deconv_stats4th_no23_allbstnoi.mat';
load(strcat(rstpath, '/MAPS/',savefile));
nsrcn = allbstnoi.nsrc;
impn = allbstnoi.impindepall;
nsrc4thn = allbstnoi.nsrc4th;
imp4thn = allbstnoi.impindep4thall;

% keyboard

%%
% %%%param for secondary sources removed
% impuse = imp;
% nsrcuse = nsrc;
% supertstr = '3-station';
% fnsuffix = [];
% impduse = impd;
% nsrcduse = nsrcd;
% impnuse = impn;
% nsrcnuse = nsrcn;

%%%param for further checked at KLNB
impuse = imp4th;
nsrcuse = nsrc4th;
supertstr = '4-station';
fnsuffix = '4th';
impduse = imp4thd;
nsrcduse = nsrc4thd;
impnuse = imp4thn;
nsrcnuse = nsrc4thn;

%% fraction of 'isolated' events, eg, when m=1, evts whose minimum interevt time >0.375s
%%%2 definitions of inter-event times, one is what we have been used
%%%the other is the smaller one of the diff time to the left and right
m = 1;
dtcut = 0.25*m+0.125;
mindtinter = cell(nnsat,nround);
dtinter = cell(nnsat,nround);
for iround = 1: nround
  for insat = 1 : nnsat
    
    nsrci = nsrcuse{insat,iround};
    impi = impuse{insat,iround};
    
    %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
    dtfor = diffcustom(impi(:,1), m,'forward'); %to its preceding one
    dtforpad = [zeros(m,1); dtfor];
    dtback = diffcustom(impi(:,1), m,'backward'); %to its following one
    dtbackpad = [dtback; zeros(m,1)];
    tmp1 = [dtforpad dtbackpad];  %time to N-m and N+m for each N
    %choose the min time to neighbors to find isolated ones
    tmp2 = [dtbackpad(1:m); min(tmp1(m+1: end-m, :),[],2); dtforpad(1:m)];

    mindtinter{insat,iround} = tmp2;
    dtinter{insat,iround} = dtfor;
  end
end

%%%for data
mindtinterd=[];
dtinterd=[];
for i = 1: nbst
  if nsrcduse(i) == 0
    continue
  end
  ist = sum(nsrcduse(1:i-1))+1;
  ied = ist+nsrcduse(i)-1;
  impi = impduse(ist:ied,:);
  if nsrcduse(i)>=m+1
    %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
    dtfor = diffcustom(impi(:,1), m,'forward'); %to its preceding one
    dtforpad = [zeros(m,1); dtfor];
    dtback = diffcustom(impi(:,1), m,'backward'); %to its following one
    dtbackpad = [dtback; zeros(m,1)];
    tmp1 = [dtforpad dtbackpad];  %time to N-m and N+m for each N
    %choose the min time to neighbors to find isolated ones
    tmp2 = [dtbackpad(1:m); min(tmp1(m+1: end-m, :),[],2); dtforpad(1:m)];
  else
    tmp2 = [];
  end
  mindtinterd = [mindtinterd; tmp2];
  dtinterd = [dtinterd; dtfor];
end

%%%for pure noise
mindtintern=[];
dtintern=[];
for i = 1: nbst
  if nsrcnuse(i) == 0
    continue
  end
  ist = sum(nsrcnuse(1:i-1))+1;
  ied = ist+nsrcnuse(i)-1;
  impi = impnuse(ist:ied,:);
  if nsrcnuse(i)>=m+1
    %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
    dtfor = diffcustom(impi(:,1), m,'forward'); %to its preceding one
    dtforpad = [zeros(m,1); dtfor];
    dtback = diffcustom(impi(:,1), m,'backward'); %to its following one
    dtbackpad = [dtback; zeros(m,1)];
    tmp1 = [dtforpad dtbackpad];  %time to N-m and N+m for each N
    %choose the min time to neighbors to find isolated ones
    tmp2 = [dtbackpad(1:m); min(tmp1(m+1: end-m, :),[],2); dtforpad(1:m)];
  else
    tmp2 = [];
  end
  mindtintern = [mindtintern; tmp2];
  dtintern = [dtintern; dtfor];
end


%% inter-event time, no amp binning
% binedge=0:0.1:12;
binedge=[0 0.375:0.25:12];  % 1st bin [0 0.375], then 0.25 increment
color = gradientblue(nround);

%%%Time from each to its preceding
nrow = 2; ncol = 4;
widin = 8.4;  % maximum width allowed is 8.5 inches
htin = 5 ;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);
pltxran = [0.08 0.98]; pltyran = [0.08 0.98]; % optimal axis location
pltxsep = 0.01; pltysep = 0.02;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

p=[]; label=[];
for insat = 1 : nnsat
  
  ax=f.ax(insat); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); set(ax,'YScale','log');
  for iround = 1: nround
    [N]=histcounts(dtinter{insat,iround}/sps,binedge,'normalization','probability');
    N=[N N(end)];
    p(iround)=stairs(ax,binedge,N,'-','linew',1,'color',color(iround,:));
    % p1=histogram(ax,dtinter/sps,'binedges',binedge,'normalization','count','Facec','b');
%     fitstruct=robustexpfit(binedge(6:end),N1d(6:end),'log10',[1e3 -1]);
%     fitobj=fitstruct.fitobj;
%     coef=coeffvalues(fitobj); slp1d=coef(2);
%     yfit1d = feval(fitobj,binedge);
%     plot(ax,binedge,yfit1d,'b--');
    if ~singleflag
      label{iround} = sprintf('%.1fx%.1f',2*semia(iround),2*semib(iround));
    else
      label{iround} = sprintf('%.1f',perctrial(iround));
    end
  end
  aa = dtinterd/sps;
  [N]=histcounts(aa,binedge,'normalization','Probability');
  N=[N N(end)];
  p(nround+1)=stairs(ax,binedge,N,'-','linew',1.5,'color','k');
  label{nround+1} = 'Data';  
  if singleflag
    aa = dtintern/sps;
    [N]=histcounts(aa,binedge,'normalization','Probability');
    N=[N N(end)];
    p(nround+2)=stairs(ax,binedge,N,'-','linew',1.5,'color','r');
    label{nround+2} = 'Noise';
  end
  text(ax,0.99,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','right');
  text(ax,0.01,0.95,supertstr,'HorizontalAlignment','left','Units','normalized',...
    'fontsize',10);
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
    % ylabel(ax,'Count');
    ylabel(ax,'Probability');
  else
    nolabels(ax,3);
  end
  xlim(ax,[0 7]);
  ylim(ax,[1e-4 1e0]);

end

fname = strcat('dt2preceding_syn',fnsuffix,fnsuffix1,'.pdf');
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));
% keyboard

%%
%%%Time from each to nearest neighbor
nrow = 2; ncol = 4;
widin = 8.4;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);
pltxran = [0.08 0.98]; pltyran = [0.08 0.98]; % optimal axis location
pltxsep = 0.01; pltysep = 0.02;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

for insat = 1 : nnsat
  
  ax=f.ax(insat); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); set(ax,'YScale','log');
  for iround = 1: nround 
    [N]=histcounts(mindtinter{insat,iround}/sps,binedge,'normalization','Probability');
    N=[N N(end)];
    p(iround)=stairs(ax,binedge,N,'-','linew',1,'color',color(iround,:));
    % p1=histogram(ax,dtinter/sps,'binedges',binedge,'normalization','count','Facec','b');
%     fitstruct=robustexpfit(binedge(6:end),N1d(6:end),'log10',[1e3 -1]);
%     fitobj=fitstruct.fitobj;
%     coef=coeffvalues(fitobj); slp1d=coef(2);
%     yfit1d = feval(fitobj,binedge);
%     plot(ax,binedge,yfit1d,'b--');
  end
  aa = mindtinterd/sps;
  [N]=histcounts(aa,binedge,'normalization','Probability');
  N=[N N(end)];
  p(nround+1)=stairs(ax,binedge,N,'-','linew',1.5,'color','k');
  if singleflag
    aa = mindtintern/sps;
    [N]=histcounts(aa,binedge,'normalization','Probability');
    N=[N N(end)];
    p(nround+2)=stairs(ax,binedge,N,'-','linew',1.5,'color','r');
  end
  text(ax,0.99,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','right');
  text(ax,0.01,0.95,supertstr,'HorizontalAlignment','left','Units','normalized',...
    'fontsize',10);
  if insat == (nrow-1)*ncol+1 
    lgd=legend(ax,p,label,'NumColumns',2,'Location','southwest',...
      'fontsize',6);
    set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
    if ~singleflag
      lgdtit = 'Region size (km)';
    else
      lgdtit = 'Noise level';
    end
    title(lgd,lgdtit,'fontsize',7);
    xlabel(ax,'Time (s) from each to nearest neighbor');
    % ylabel(ax,'Count');
    ylabel(ax,'Probability');
  else
    nolabels(ax,3);
  end
  xlim(ax,[0 7]);
  ylim(ax,[1e-4 1e0]);
end

fname = strcat('dt2nearest_syn',fnsuffix,fnsuffix1,'.pdf');
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));

