function f = plt_selsigpred(sigsta,predgrp,varred,stas,pltsta,sps,xzoom,detecttype,saveflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f = plt_selsigpred(sigsta,predgrp,varred1,stas,pltsta,sps,xzoom,detecttype,saveflag)
%
% This is a simple function to plot signal, prediction from deconvolved 
% sources, at selected station indices defined in 'pltsta'. 
% Text out the relative change in variance, or misfit.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/11/17
% Last modified date:   2022/11/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('sps',160);
defval('xzoom',[]);
defval('detecttype',[]);  %default is short-win detections
defval('saveflag',1);

nsta = length(pltsta);
if nsta ==3
  panelnm = ['a';'b';'c'];
elseif nsta ==4
  panelnm = ['a';'b';'c';'d'];
end
widin = 8.3;  % maximum width allowed is 8.5 inches
htin = 1.0*nsta;   % maximum height allowed is 11 inches
nrow = nsta;
ncol = 1;
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.05 0.98]; pltyran = [0.1 0.98]; % optimal axis location
pltxsep = 0.05; pltysep = 0.03;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

% ym = max(abs(sigsta(:)));
% yran=1.2*[-ym ym];
ym = 1;
yran= [-ym ym];

lsig = size(sigsta,1);

varred = squeeze(varred);

for jj = 1:length(pltsta)
  ista = pltsta(jj);
  ax = f.ax(jj);
  hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
  p(1)=plot(ax,(1:lsig)/sps,sigsta(:,ista),'k','linew',1);
  p(2)=plot(ax,(1:lsig)/sps,predgrp(:,ista),'r','linew',1);
%   text(ax,0.98,0.9,sprintf('%.2f; %.2f; %.1f%%',l2normred(ista,1),l2normred(ista,2),l2normred(ista,3)),...
%     'Units','normalized','HorizontalAlignment','right','fontsize',10);
  text(ax,0.99,0.1,sprintf('VR: %.1f%%',varred(ista,3)),...
    'Units','normalized','HorizontalAlignment','right','fontsize',10);
  text(ax,0.99,0.9,sprintf('%s',strtrim(stas(ista,:))),...
    'Units','normalized','HorizontalAlignment','right','fontsize',12);
  ylim(ax,yran);
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
  end
  ylabel(ax,'Amplitude','fontsize',10);
  if jj==1
    legend(ax,p,'Signal','Prediction','Location','south',...
      'Orientation','horizontal');%,'Residual'
  end
  text(ax,0.01,0.82,panelnm(jj,:),'FontSize',10,'unit','normalized',...
    'EdgeColor','k','Margin',1,'backgroundcolor','w');
end

if saveflag
  if nsta == 3
  %   orient(f.fig,'landscape');
    fname = strcat('varredaftgrp',detecttype,'.pdf');
    print(f.fig,'-dpdf',...
      strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
  elseif nsta == 4 && strcmp(strtrim(stas(pltsta(end),:)),'KLNB')
  %   orient(f.fig,'landscape');
    fname = strcat('varred4th',detecttype,'.pdf');
    print(f.fig,'-dpdf',...
      strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
    
  end
end

