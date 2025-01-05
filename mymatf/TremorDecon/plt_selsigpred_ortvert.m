function f = plt_selsigpred_ortvert(sigstaort,predgrport,varredort,...
  sigstavert,predgrpvert,varredvert,greenfort,greenfvert,...
  stas,pltsta,sps,xzoom,patchwins,detecttype,saveflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f = plt_selsigpredres(sigsta,predgrp,resgrp,varred1,stas,pltsta,sps,xzoom,detecttype,saveflag)
%
% This is a simple function to plot signal, prediction from deconvolved 
% sources, and residual at selected station indices defined in 'pltsta'. 
% Text out the relative change in variance, or misfit.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/06/24
% Last modified date:   2024/06/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('sps',160);
defval('xzoom',[]);
defval('detecttype',[]);  %default is short-win detections
defval('saveflag',1);

nsta = length(pltsta);
lwlet = 12*sps;
lsigzoom = range(xzoom)*sps;

% ym = max(abs(sigstavert(:)));
% yran=1.5*[-ym ym];
ym = 0.5;
yran = [-ym ym];

lsig = size(sigstavert,1);

varredort = squeeze(varredort);
varredvert = squeeze(varredvert);

nrow = nsta+2;
ncol = 1;
widin = 8.3;  % maximum width allowed is 8.5 inches
htin = 1.4*(nsta+1);   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.07 0.99]; pltyran = [0.08 0.98]; % optimal axis location
pltxsep = 0.05; pltysep = 0.015;
ywid = (range(pltyran)-nsta*pltysep)/(nsta+1)-0.01;
xwid = range(pltxran);
if nsta ==3
  color = ['r';'b';'k'];
  lgdlbl = {'PGC';'SSIB';'SILB'};
  panelnm = ['a';'b';'c';'d';'e'];
  axpos = [pltxran(1) pltyran(1)+(ywid+pltysep)*3+0.03 xwid*lwlet/lsigzoom ywid;
           pltxran(1)+xwid*(1-lwlet/lsigzoom) pltyran(1)+(ywid+pltysep)*3+0.03 xwid*lwlet/lsigzoom ywid;
           pltxran(1) pltyran(1)+(ywid+pltysep)*2 xwid ywid;
           pltxran(1) pltyran(1)+ywid+pltysep xwid ywid;
           pltxran(1) pltyran(1) xwid ywid;
           ];
elseif nsta ==4
  color = ['r';'b';'k';'c'];
  lgdlbl = {'PGC';'SSIB';'SILB';'KLNB'};
  panelnm = ['a';'b';'c';'d';'e';'f'];
  axpos = [pltxran(1) pltyran(1)+(ywid+pltysep)*4+0.03 xwid*lwlet/lsigzoom ywid;
           pltxran(1)+xwid*(1-lwlet/lsigzoom) pltyran(1)+(ywid+pltysep)*4+0.03 xwid*lwlet/lsigzoom ywid;
           pltxran(1) pltyran(1)+(ywid+pltysep)*3 xwid ywid;
           pltxran(1) pltyran(1)+(ywid+pltysep)*2 xwid ywid;
           pltxran(1) pltyran(1)+ywid+pltysep xwid ywid;
           pltxran(1) pltyran(1) xwid ywid;
           ];
end
for i = 1:nrow*ncol
  set(f.ax(i), 'position', axpos(i,:));
end

%%%LFE templates, bandpased, orthogonal
ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on'); 
for i = 1: nsta
  p(i)=plot(ax,(1:lwlet)/sps, greenfort(1:lwlet,i), '-','Color',color(i,:),'linew',1); 
end
xlim(ax,[0 lwlet/sps]);  
% ylim(ax,yran/2);
ylim(ax,[-0.1 0.1]);
% shrink(ax,range(xzoom)*sps/lwlet,1);
% ax.Position(1)=axpos(1,1);
xticks(ax,0: 2: lwlet/sps);
if ~isempty(patchwins)
  winnoi=patchwins(1,:);
  winsig=patchwins(2,:);
  %patch a gray shaded part for zoomed-in in other plots
  patnoi = [winnoi(1) ax.YLim(2);
    winnoi(2) ax.YLim(2);
    winnoi(2) ax.YLim(1);
    winnoi(1) ax.YLim(1);
    winnoi(1) ax.YLim(2)];
  patch(ax,patnoi(:,1),patnoi(:,2),'k','Facealpha',0.15,'edgecolor','none');
  patsig = [winsig(1) ax.YLim(2);
    winsig(2) ax.YLim(2);
    winsig(2) ax.YLim(1);
    winsig(1) ax.YLim(1);
    winsig(1) ax.YLim(2)];
  patch(ax,patsig(:,1),patsig(:,2),'c','Facealpha',0.15,'edgecolor','none');
  
  staplt=1:nsta;
  envf = envelope(detrend(greenfort(:,staplt)));
  snrf = range(envf(winsig(1)*sps: winsig(2)*sps, :),1)./ ...
    median(envf(winnoi(1)*sps: winnoi(2)*sps, :),1);
  for ista = 1: nsta
    text(ax,1+(ista-1-nsta)*0.1,0.7,sprintf('%.1f;',snrf(:,ista)),'Units',...
      'normalized','FontSize',9,'color',color(ista,:));
  end  
end
legend(ax,p,lgdlbl,'Location','south','Orientation','horizontal',...
  'fontsize',8);
text(ax,0.98,0.9,'Orthogonal templates','Units','normalized','HorizontalAlignment','right',...
  'FontSize',10);
% text(ax,0.98,0.9,'LFE templates at family 002','Units','normalized',...
%   'HorizontalAlignment','right','FontSize',10);
text(ax,0.1,0.9,'1.8-18 Hz','Units','normalized','HorizontalAlignment','left',...
  'FontSize',10);
text(ax,0.02,0.88,panelnm(1,:),'FontSize',10,'unit','normalized',...
    'EdgeColor','k','Margin',1,'backgroundcolor','w'); 
longticks(ax,2); 
% nolabels(ax,1);

%%%LFE templates, bandpased, orthogonal
ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on'); 
for i = 1: nsta
  p(i)=plot(ax,(1:lwlet)/sps, greenfvert(1:lwlet,i), '-','Color',color(i,:),'linew',1); 
end
xlim(ax,[0 lwlet/sps]);  
% ylim(ax,yran/2);
ylim(ax,[-0.1 0.1]);
% shrink(ax,range(xzoom)*sps/lwlet,1);
% ax.Position(1)=axpos(1,1);
xticks(ax,0: 2: lwlet/sps);
if ~isempty(patchwins)
  winnoi=patchwins(1,:);
  winsig=patchwins(2,:);
  %patch a gray shaded part for zoomed-in in other plots
  patnoi = [winnoi(1) ax.YLim(2);
    winnoi(2) ax.YLim(2);
    winnoi(2) ax.YLim(1);
    winnoi(1) ax.YLim(1);
    winnoi(1) ax.YLim(2)];
  patch(ax,patnoi(:,1),patnoi(:,2),'k','Facealpha',0.15,'edgecolor','none');
  patsig = [winsig(1) ax.YLim(2);
    winsig(2) ax.YLim(2);
    winsig(2) ax.YLim(1);
    winsig(1) ax.YLim(1);
    winsig(1) ax.YLim(2)];
  patch(ax,patsig(:,1),patsig(:,2),'c','Facealpha',0.15,'edgecolor','none');
  
  staplt=1:nsta;
  envf = envelope(detrend(greenfvert(:,staplt)));
  snrf = range(envf(winsig(1)*sps: winsig(2)*sps, :),1)./ ...
    median(envf(winnoi(1)*sps: winnoi(2)*sps, :),1);
  for ista = 1: nsta
    text(ax,1+(ista-1-nsta)*0.1,0.7,sprintf('%.1f;',snrf(:,ista)),'Units',...
      'normalized','FontSize',9,'color',color(ista,:));
  end  
end
% legend(ax,p,{'PGC','SSIB','SILB'},'Location','southeast','Orientation','horizontal',...
%   'fontsize',8);
text(ax,0.1,0.9,'1.8-18 Hz','Units','normalized','HorizontalAlignment','left',...
  'FontSize',10);
text(ax,0.98,0.9,'Vertical templates','Units','normalized','HorizontalAlignment','right',...
  'FontSize',10);
% text(ax,0.98,0.9,'LFE templates at family 002','Units','normalized',...
%   'HorizontalAlignment','right','FontSize',10);
text(ax,0.02,0.88,panelnm(2,:),'FontSize',10,'unit','normalized',...
    'EdgeColor','k','Margin',1,'backgroundcolor','w'); 
longticks(ax,2); 
nolabels(ax,2);

p=[];
for jj = 3:nsta+2
  ista = pltsta(jj-2);
  ax = f.ax(jj);
  hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
  p(1)=plot(ax,(1:lsig)/sps,sigstaort(:,ista)+0.5*ym,'k','linew',1);
  p(2)=plot(ax,(1:lsig)/sps,predgrport(:,ista)+0.5*ym,'r','linew',1);
  plot(ax,(1:lsig)/sps,sigstavert(:,ista)-0.5*ym,'k','linew',1);
  plot(ax,(1:lsig)/sps,predgrpvert(:,ista)-0.5*ym,'r','linew',1);
%   text(ax,0.98,0.9,sprintf('%.2f; %.2f; %.1f%%',l2normred(ista,1),l2normred(ista,2),l2normred(ista,3)),...
%     'Units','normalized','HorizontalAlignment','right','fontsize',10);
  text(ax,0.05,0.9,sprintf('Orthogonal, VR: %.1f%%',varredort(ista,3)),...
    'Units','normalized','HorizontalAlignment','left','fontsize',10);
  text(ax,0.05,0.1,sprintf('Vertical, VR: %.1f%%',varredvert(ista,3)),...
    'Units','normalized','HorizontalAlignment','left','fontsize',10);
  text(ax,0.99,0.92,sprintf('%s',strtrim(stas(ista,:))),...
    'Units','normalized','HorizontalAlignment','right','fontsize',12);
  ylim(ax,yran);
  yticks(ax,-0.25:0.25:0.25);
  if isempty(xzoom)
    xlim(ax,[0 lsig]/sps);
  else
    xlim(ax,xzoom);
    xticks(ax,xzoom(1): 2: xzoom(2));
  end
  longticks(ax,6);
  if jj~=nrow
    nolabels(ax,1);
  else
    xlabel(ax,'Time (s)','fontsize',10);
    ylabel(ax,'Amplitude','fontsize',10);
  end
  if jj==3
    legend(ax,p,'Signal','Prediction','Location','west',...
      'Orientation','horizontal');
  end
  text(ax,0.01,0.88,panelnm(jj,:),'FontSize',10,'unit','normalized',...
    'EdgeColor','k','Margin',1,'backgroundcolor','w'); 
end

if saveflag
  if nsta == 3
    fname = strcat('vrgrp',detecttype,'_ort_vert.pdf');
    print(f.fig,'-dpdf',...
      strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
  elseif nsta == 4 && strcmp(strtrim(stas(pltsta(end),:)),'KLNB')
    fname = strcat('vr4th',detecttype,'_ort_vert.pdf');
    print(f.fig,'-dpdf',...
      strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));    
  end
end

