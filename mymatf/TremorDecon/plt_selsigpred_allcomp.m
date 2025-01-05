function f = plt_selsigpred_allcomp(sigsta,predgrp,varred,sigstaort,predgrport,...
  varredort,sigstavert,predgrpvert,varredvert,greenfort,greenfvert,...
  stas,pltsta,sps,xzoom,impindepst,patchwins,normflag,detecttype,saveflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f = plt_selsigpred_allcomp(sigsta,predgrp,varred,sigstaort,predgrport,...
%   varredort,sigstavert,predgrpvert,varredvert,greenfort,greenfvert,...
%   stas,pltsta,sps,xzoom,impindepst,patchwins,normflag,detecttype,saveflag)
%
% This is a function to plot signal, prediction from deconvolved 
% sources at selected station indices defined in 'pltsta' for ALL components
% opt, ort, vert. Then templates at ort and vert comp. are plotted.
% Text out the relative change in variance, or misfit.
% --If the source 'impindepst' is not empty, then plot the sources in terms
% of a shade area from the start of the early shoulder and end of the dipole
% to reveal the effect. 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/07/17
% Last modified date:   2024/07/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('sps',160);
defval('xzoom',[]);
defval('patchwins',[]);
defval('normflag',0);  %default is do not normalize the seismograms
defval('detecttype',[]);  %default is short-win detections
defval('saveflag',1);

nsta = length(pltsta);
lwlet = 12*sps;
lsigzoom = range(xzoom)*sps;

% ym = max(abs(sigstavert(:)));
% yran=1.5*[-ym ym];
ym = 2;
yran = [-ym ym];

sigsta = sigsta(:,pltsta);
sigstaort = sigstaort(:,pltsta);
sigstavert = sigstavert(:,pltsta);
predgrp = predgrp(:,pltsta);
predgrport = predgrport(:,pltsta);
predgrpvert = predgrpvert(:,pltsta);

%if need to scale the seismograms
if normflag
  maxval = range(sigsta(xzoom(1)*sps+1:xzoom(2)*sps, :));
  sigsta = sigsta ./ maxval;
  predgrp = predgrp ./ maxval;
  maxvalort = range(sigstaort(xzoom(1)*sps+1:xzoom(2)*sps, :));
  sigstaort = sigstaort ./ maxvalort;
  predgrport = predgrport ./ maxvalort;
  maxvalvert = range(sigstavert(xzoom(1)*sps+1:xzoom(2)*sps, :));
  sigstavert = sigstavert ./ maxvalvert;
  predgrpvert = predgrpvert ./ maxvalvert;
  % ppkindep(:,[2 4 6]) = ppkindep(:,[2 4 6]) ./ normalizer;
  % ppkindepsave(:,[2 4 6]) = ppkindepsave(:,[2 4 6]) ./ normalizer;
  sclort = maxval./maxvalort;
  sclvert = maxval./maxvalvert;
end

scale = 1;  % to scale the ort and vert. comps as they have smaller amplitude

lsig = size(sigstavert,1);

varredort = squeeze(varredort);
varredvert = squeeze(varredvert);

% keyboard

nrow = nsta+2;
ncol = 1;
widin = 8.3;  % maximum width allowed is 8.5 inches
htin = 1.6*(nsta+1);   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.07 0.99]; pltyran = [0.07 0.98]; % optimal axis location
pltxsep = 0.05; pltysep = 0.015;
ywid = (range(pltyran)-nsta*pltysep)/(nsta+1)-0.01;
xwid = range(pltxran);
if nsta ==3
  color = ['r';'b';'k'];
  lgdlbl = {'PGC';'SSIB';'SILB'};
  panelnm = ['a';'b';'c';'d';'e'];
  axpos = [
           pltxran(1) pltyran(1)+(ywid+pltysep)*3+0.05 xwid ywid;
           pltxran(1) pltyran(1)+(ywid+pltysep)*2+0.05 xwid ywid;
           pltxran(1) pltyran(1)+(ywid+pltysep)*1+0.05 xwid ywid;
           pltxran(1) pltyran(1) xwid*lwlet/lsigzoom ywid;
           pltxran(1)+xwid*(1-lwlet/lsigzoom) pltyran(1) xwid*lwlet/lsigzoom ywid;
           ];
elseif nsta ==4
  color = ['r';'b';'k';'c'];
  lgdlbl = {'PGC';'SSIB';'SILB';'KLNB'};
  panelnm = ['a';'b';'c';'d';'e';'f'];
  axpos = [
           pltxran(1) pltyran(1)+(ywid+pltysep)*4+0.05 xwid ywid;
           pltxran(1) pltyran(1)+(ywid+pltysep)*3+0.05 xwid ywid;
           pltxran(1) pltyran(1)+(ywid+pltysep)*2+0.05 xwid ywid;
           pltxran(1) pltyran(1)+(ywid+pltysep)*1+0.05 xwid ywid;
           pltxran(1) pltyran(1) xwid*lwlet/lsigzoom ywid;
           pltxran(1)+xwid*(1-lwlet/lsigzoom) pltyran(1) xwid*lwlet/lsigzoom ywid;           ];
end
for i = 1:nrow*ncol
  set(f.ax(i), 'position', axpos(i,:));
end

%%%signals, predictions, at 3 components, top, mid, bottom
for jj = 1:nsta
  ista = jj;
  ax = f.ax(jj);
  hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
  ylim(ax,yran);
%   yticks(ax,-0.25:0.25:0.25);
  if isempty(xzoom)
    xlim(ax,[0 lsig]/sps);
  else
    xlim(ax,xzoom);
    xticks(ax,xzoom(1): 2: xzoom(2));
  end
  %decon sources
  ind = find(impindepst(:,1)/sps >= xzoom(1) & impindepst(:,1)/sps <= xzoom(2));
  %a region around the zero-crossing supposed to affect most the decon of sources
  % zc2lt1ppk = 64; %from dipole zero-crossing to the left, till the zc start of 1st pos peak (shoulder)
  zc2lt1ppk = 60; %from dipole zero-crossing to the left, till the zc start of 1st pos peak (shoulder)
  zc2rt2npk = 56; %from dipole zero-crossing to the right, till the zc end of 2nd neg peak
  zc2lt1npkst = 24; %from dipole zero-crossing to the left, till the zc start of 1st neg peak(dipole)
  zc2rt2ppked = 22; %from dipole zero-crossing to the right, till the zc end of 2nd pos peak(dipole)
  
  %patch a shaded region around each source
  for isrc = 1: length(ind)
    if ista <=3 
      tcol = (pltsta(ista)-1)*2+1;
    else
      tcol = (pltsta(ista)-4)*2+10;
    end
    %option 1, from shoulder to end of pos dipole peak 
    patsrc = [(impindepst(ind(isrc),tcol)-zc2lt1ppk)/sps ax.YLim(2);
              (impindepst(ind(isrc),tcol)+zc2rt2ppked)/sps ax.YLim(2);
              (impindepst(ind(isrc),tcol)+zc2rt2ppked)/sps ax.YLim(1);
              (impindepst(ind(isrc),tcol)-zc2lt1ppk)/sps ax.YLim(1);
              (impindepst(ind(isrc),tcol)-zc2lt1ppk)/sps ax.YLim(2)];
%     %option 2, only include the neg dipole peak 
%     patsrc = [(impindepst(ind(isrc),tcol)-zc2lt1npkst)/sps ax.YLim(2);
%               (impindepst(ind(isrc),tcol)+0)/sps ax.YLim(2);
%               (impindepst(ind(isrc),tcol)+0)/sps ax.YLim(1);
%               (impindepst(ind(isrc),tcol)-zc2lt1npkst)/sps ax.YLim(1);
%               (impindepst(ind(isrc),tcol)-zc2lt1npkst)/sps ax.YLim(2)];              
    patch(ax,patsrc(:,1),patsrc(:,2),[0 100 0]/255,'Facealpha',0.15,'edgecolor','none');%192
  end
  %olive 128,128,0; dark green 0,100,0
  p(1)=plot(ax,(1:lsig)/sps,sigsta(:,ista)+0.5*ym,'k','linew',1);
  p(2)=plot(ax,(1:lsig)/sps,predgrp(:,ista)+0.5*ym,'r','linew',1);
  plot(ax,(1:lsig)/sps,sigstaort(:,ista)*scale,'k','linew',1);
  plot(ax,(1:lsig)/sps,predgrport(:,ista)*scale,'r','linew',1);
  plot(ax,(1:lsig)/sps,sigstavert(:,ista)*scale-0.5*ym,'k','linew',1);
  plot(ax,(1:lsig)/sps,predgrpvert(:,ista)*scale-0.5*ym,'r','linew',1);

  %plot decon sources
  yloc = (yran(1)+range(yran)*0.05);
  scatter(ax,impindepst(ind,tcol)/sps,yloc*ones(size(impindepst(ind,1))),10,'k');
% keyboard
%   text(ax,0.98,0.9,sprintf('%.2f; %.2f; %.1f%%',l2normred(ista,1),l2normred(ista,2),l2normred(ista,3)),...
%     'Units','normalized','HorizontalAlignment','right','fontsize',10);
  text(ax,0.05,0.85,sprintf('Optimal; VR: %.1f%%',varred(pltsta(ista),3)),...
    'Units','normalized','HorizontalAlignment','left','fontsize',8);
  text(ax,0.05,0.6,sprintf('Orthogonal, x%.1f; VR: %.1f%%',sclort(ista),...
    varredort(pltsta(ista),3)),'Units','normalized','HorizontalAlignment','left','fontsize',8);
  text(ax,0.05,0.35,sprintf('Vertical, x%.1f; VR: %.1f%%',sclvert(ista),...
    varredvert(pltsta(ista),3)),'Units','normalized','HorizontalAlignment','left','fontsize',8);
  text(ax,0.99,0.92,sprintf('%s',strtrim(stas(pltsta(ista),:))),...
    'Units','normalized','HorizontalAlignment','right','fontsize',12);
  ylabel(ax,'Norm. amplitude','fontsize',10);
  if ista~=nsta
    nolabels(ax,1);
  else
    xlabel(ax,'Time (s)','fontsize',10);
  end
  if ista==1
    legend(ax,p,'Signal','Prediction','Location','north',...
      'Orientation','horizontal','fontsize',8);
  end
  text(ax,0.01,0.88,panelnm(ista,:),'FontSize',10,'unit','normalized',...
    'EdgeColor','k','Margin',1,'backgroundcolor','w'); 
  longticks(ax,6);
end
% keyboard

%%%LFE templates, bandpased, orthogonal
ax=f.ax(nsta+1); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on'); 
for i = 1: nsta
  p(i)=plot(ax,(1:lwlet)/sps, greenfort(1:lwlet,i), '-','Color',color(i,:),'linew',1); 
end
xlim(ax,[0 lwlet/sps]);  
% ylim(ax,yran/2);
ylim(ax,[-0.1 0.1]);
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
locmaxamp = 4.8625;
plot(ax,[locmaxamp locmaxamp],[0.05 0.1],'k-','linewidth',1);
%an arrow pointing at the max amp of the opt. templates
[rotx, roty] = complex_rot(0,0.05,180);
xarrow = [locmaxamp locmaxamp+rotx];
yarrow = [0.1 0.1+roty];
a=annotation('arrow','color','k','linestyle','-','linewidth',1,'HeadLength',6,'HeadWidth',6);%
a.Parent = ax;
% a.Position = [locmaxamp, 0.1, 0, -0.05];
a.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];

lgd=legend(ax,p,lgdlbl,'Location','south','Orientation','horizontal',...
  'fontsize',8);
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
text(ax,0.98,0.9,'Orthogonal templates','Units','normalized','HorizontalAlignment','right',...
  'FontSize',10);
% text(ax,0.98,0.9,'LFE templates at family 002','Units','normalized',...
%   'HorizontalAlignment','right','FontSize',10);
text(ax,0.1,0.9,'1.8-18 Hz','Units','normalized','HorizontalAlignment','left',...
  'FontSize',10);
text(ax,0.02,0.88,panelnm(nsta+1,:),'FontSize',10,'unit','normalized',...
    'EdgeColor','k','Margin',1,'backgroundcolor','w'); 
longticks(ax,2);
xlabel(ax,'Time (s)','fontsize',10);
ylabel(ax,'Amplitude','fontsize',10);
% shrink(ax,range(xzoom)*sps/lwlet,1);
% ax.Position(1)=axpos(1,1);
xticks(ax,0: 2: lwlet/sps);
% nolabels(ax,1);

%%%LFE templates, bandpased, vertical
p=[];
ax=f.ax(nsta+2); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on'); 
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
locmaxamp = 4.8625;
plot(ax,[locmaxamp locmaxamp],[0.05 0.1],'k-','linewidth',1);
%an arrow pointing at the max amp of the opt. templates
[rotx, roty] = complex_rot(0,0.05,180);
xarrow = [locmaxamp locmaxamp+rotx];
yarrow = [0.1 0.1+roty];
a=annotation('arrow','color','k','linestyle',':','linewidth',1,'HeadLength',6,...
  'HeadWidth',6);
a.Parent = ax;
a.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];

% legend(ax,p,{'PGC','SSIB','SILB'},'Location','southeast','Orientation','horizontal',...
%   'fontsize',8);
text(ax,0.1,0.9,'1.8-18 Hz','Units','normalized','HorizontalAlignment','left',...
  'FontSize',10);
text(ax,0.98,0.9,'Vertical templates','Units','normalized','HorizontalAlignment','right',...
  'FontSize',10);
% text(ax,0.98,0.9,'LFE templates at family 002','Units','normalized',...
%   'HorizontalAlignment','right','FontSize',10);
text(ax,0.02,0.88,panelnm(nsta+2,:),'FontSize',10,'unit','normalized',...
    'EdgeColor','k','Margin',1,'backgroundcolor','w'); 
xlabel(ax,'Time (s)','fontsize',10);
longticks(ax,2); 
nolabels(ax,2);

if saveflag
  if nsta == 3
  %   orient(f.fig,'landscape');
    fname = strcat('vrgrp',detecttype,'.pdf');
    print(f.fig,'-dpdf',...
      strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
  elseif nsta == 4 && strcmp(strtrim(stas(pltsta(end),:)),'KLNB')
  %   orient(f.fig,'landscape');
    fname = strcat('vr4th',detecttype,'.pdf');
    print(f.fig,'-dpdf',...
      strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
    
  end
end

