function f = plt_impsigpredres_1sta(sigdecon,sigsta,pred,rcc,stas,ista,sps,...
  xzoom,detecttype,saveflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f = plt_selsigpredres(sigsta,predgrp,resgrp,varred1,stas,pltsta,sps)
%
% At a user specified station, plot the impulses, signal, prediction, 
% and residual. Additionally, the running CC between sig and pred is 
% computed and plotted. Note that this script is designed to be called
% RIGHT AFTER deconvolution at one station to easily see how the 
% the signal is explained (quantified by the variance reduction) by the
% deconvolution. A similar and slightly more informative plot is 'f2' in
% 'iterdecon.m'. 
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/03/08
% Last modified date:   2024/03/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('ista',1);
defval('sps',160);
defval('xzoom',[]);
defval('detecttype',[]);  %default is short-win detections
defval('saveflag',1);  %default is short-win detections

lsig = size(sigsta,1);

widin = 8.3;  % maximum width allowed is 8.5 inches
htin = 4;   % maximum height allowed is 11 inches
nrow = 3;
ncol = 1;
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.06 0.94]; pltyran = [0.1 0.98]; % optimal axis location
pltxsep = 0.05; pltysep = 0.025;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'off');
yyaxis(ax,'left');
imptemp = find(sigdecon(:,ista)>0);
p1=stem(ax,imptemp/sps,sigdecon(imptemp,ista),'Color',[1 .3 0],'MarkerSize',2.5);
% text(ax,0.9,0.9,num2str(length(imptemp)),'unit','normalized',...
%   'HorizontalAlignment','right','FontSize',10);
text(ax,0.98,0.9,strcat(stas(ista,:)),'unit','normalized',...
  'HorizontalAlignment','right','FontSize',12);
text(ax,0.01,0.85,'a','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
ylim(ax,[0 4.5]);
if isempty(xzoom)
  xlim(ax,[0 lsig]/sps);
else
  xlim(ax,xzoom);
  xticks(ax,xzoom(1): 2: xzoom(2));
end
ax.YColor=ax.XColor;
ylabel(ax,'Amplitude','fontsize',10);
yyaxis(ax,'right');
ircc = 1:length(rcc);
% ind = find(ircc/sps>=xroom())
p2=plot(ax,ircc/sps,rcc,'o','color',[.8 .8 .8],'markersize',1);
legend(ax,[p1,p2],'Impulses','Original RCC','location','north','Orientation','horizontal');
ylabel(ax,'Running CC','FontSize',10);
ylim(ax,[-1 1]);
longticks(ax,6);
nolabels(ax,1);
ax.YColor=[.6 .6 .6];
hold(ax,'off');

ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
ym = max(abs([sigsta(:,ista); pred(:,ista)]));
yran=1.2*[-ym ym];
yyaxis(ax,'left');
p1=plot(ax,(1:lsig)/sps,sigsta(:,ista), 'k-','linew',1);
p2=plot(ax,(1:lsig)/sps,pred(:,ista), 'r-','linew',1);
text(ax,0.01,0.15,'b','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
ylim(ax,yran);
if isempty(xzoom)
  xlim(ax,[0 lsig]/sps);
else
  xlim(ax,xzoom);
  xticks(ax,xzoom(1): 2: xzoom(2));
end
ax.YColor=ax.XColor;
ylabel(ax,'Amplitude','fontsize',10);
yyaxis(ax,'right');
rccmwlen = sps/2;
[irccsp, rccsp] = RunningCC(sigsta(:,ista), pred(:,ista), rccmwlen);
p3=plot(ax,irccsp/sps,rccsp,'o','color',[.8 .8 .8],'markersize',1);
ylabel(ax,'Running CC','FontSize',10);
ylim(ax,[-1 1]);
legend(ax,[p1,p2,p3],'Signal','Prediction','Sig-Pred RCC','location',...
  'south','Orientation','horizontal');
longticks(ax,6);
nolabels(ax,1);
ax.YColor=[.6 .6 .6];
hold(ax,'off');

ax=f.ax(3); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
res=sigsta(:,ista)-pred(:,ista);
% yyaxis(ax,'left');
p1=plot(ax,(1:lsig)/sps,res,'Color',[.5 .5 .5],'linew',1);
varred=(var(sigsta(:,ista))-var(res))/var(sigsta(:,ista))*100;
text(ax,0.98,0.1,sprintf('VR: %.1f%%',varred),...
  'Units','normalized','HorizontalAlignment','right','fontsize',10);
text(ax,0.01,0.85,'c','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
ylim(ax,yran);
if isempty(xzoom)
  xlim(ax,[0 lsig]/sps);
else
  xlim(ax,xzoom);
  xticks(ax,xzoom(1): 2: xzoom(2));
end
xlabel(ax,'Time (s)','fontsize',10);
ylabel(ax,'Amplitude','fontsize',10);
% ax.YColor=ax.XColor;
% yyaxis(ax,'right');
% [irccsr, rccsr] = RunningCC(sigsta(:,ista), res, rccmwlen);
% p2=plot(ax,irccsr/sps,rccsr,'o','color',[.8 .8 .8],'markersize',1);
% ylabel(ax,'Running CC','FontSize',10);
% ylim(ax,[-1 1]);
legend(ax,[p1],'Residual','location','south');
% legend(ax,[p1,p2],'Residual','Sig-Res RCC','location','south');
longticks(ax,6);
hold(ax,'off');

if saveflag
  % orient(f.fig,'landscape');
  fname = strcat('sigvspredaftdecon_',strtrim(stas(ista,:)),detecttype,'.pdf');
  print(f.fig,'-dpdf',...
    strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
end







