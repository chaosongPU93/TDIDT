function f1 = plt_selsigpredres(sigsta,predgrp,resgrp,l2normred,stas,pltsta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f1 = plt_selsigpredres(sigsta,predgrp,resgrp,resred,pltsta)
%
% This is a simple function to plot signal, prediction from deconvolved 
% sources, and residual at selected station indices defined in 'pltsta'. 
% Text out the relative change in L2-norm, or misfit.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/11/17
% Last modified date:   2022/11/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

widin = 12;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
nrow = length(pltsta);
ncol = 1;
f1 = initfig(widin,htin,nrow,ncol);

pltxran = [0.06 0.96]; pltyran = [0.06 0.96];
pltxsep = 0.02; pltysep = 0.02;
optaxpos(f1,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

ym = max(abs(sigsta(:)));
yran=1.3*[-ym ym];

l2normred = squeeze(l2normred);

for jj = 1:length(pltsta)
  ista = pltsta(jj);
  ax = f1.ax(jj);
  hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
  p1=plot(ax,sigsta(:,ista)+0.15,'k');
  p2=plot(ax,predgrp(:,ista)+0.15,'r');
  p3=plot(ax,resgrp(:,ista)-0.15,'Color',[.7 .7 .7]);
  text(ax,0.98,0.9,sprintf('%.2f; %.2f; %.1f%%',l2normred(ista,1),l2normred(ista,2),l2normred(ista,3)),...
    'Units','normalized','HorizontalAlignment','right','fontsize',10);
  text(ax,0.98,0.1,sprintf('%s',strtrim(stas(ista,:))),...
    'Units','normalized','HorizontalAlignment','right','fontsize',11);
  ylim(ax,yran);
  xlim(ax,[0 size(sigsta,1)]);
  longticks(ax,6);
  if jj~=nrow
    nolabels(ax,1);
  else
    xlabel(ax,'Samples','fontsize',11);
    ylabel(ax,'Amplitude','fontsize',11);
  end
  if jj==1
    legend(ax,[p1,p2,p3],'Signal','Prediction','Residual','Location','west',...
      'NumColumns',3)
  end
end


