function [f,snrf,snrfort,snrfvert]=plt_templates_bp_bysta(greenf,greenfort,greenfvert,stas,staplt,...
  lowlet,hiwlet,sps,patchwins)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plt_templates_bp_bysta(greenf,greenfort,greenfvert,stas,staplt,...
%   lowlet,hiwlet,sps,patchwins)
%
% this function simply plot the bandpassed filtered templates (or green's 
% functions), on all input components (optimal, orthogonal, vertical). All in
% one panel for each station.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/09/12
% Last modified date:   2022/09/12 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('patchwins',[]);

compname = {'Opt.','Ort.','Vert.'};

nsta = length(staplt);
ncol = 2;
nrow = ceil(nsta/ncol);
widin = 8.4;  % maximum width allowed is 8.5 inches
htin = 1.5*nrow;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.06 0.98]; pltyran = [0.1 0.98]; % optimal axis location
pltxsep = 0.04; pltysep = 0.015;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

greenf=greenf(:,staplt);
greenfort=greenfort(:,staplt);
greenfvert=greenfvert(:,staplt);

if ~isempty(patchwins)
  winnoi=patchwins(1,:);
  winsig=patchwins(2,:);
  envf = envelope(detrend(greenf));
  snrf = range(envf(winsig(1)*sps: winsig(2)*sps, :),1)./ ...
    median(envf(winnoi(1)*sps: winnoi(2)*sps, :),1);
  envfort = envelope(detrend(greenfort));
  snrfort = range(envfort(winsig(1)*sps: winsig(2)*sps, :),1)./ ...
    median(envfort(winnoi(1)*sps: winnoi(2)*sps, :),1);
  envfvert = envelope(detrend(greenfvert));
  snrfvert = range(envfvert(winsig(1)*sps: winsig(2)*sps, :),1)./ ...
    median(envfvert(winnoi(1)*sps: winnoi(2)*sps, :),1);
else
  snrf = []; snrfort = []; snrfort = []; 
end

mx=1.2*max(max(abs(greenf(:,:))));

for i = 1: nsta
  ax = f.ax(i);
  hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%   yyaxis(ax,'left');
  lwlet = size(greenf,1);
%   if nsta == 3
    color=['r';'b';'k'];
%     p(1)=plot(ax,(1:lwlet)/sps,greenplt(:,1),'r-','linew',1);
%     p(2)=plot(ax,(1:lwlet)/sps,greenplt(:,2),'b-','linew',1);
%     p(3)=plot(ax,(1:lwlet)/sps,greenplt(:,3),'k-','linew',1);
%   elseif nsta == 4
%     color=['r';'b';'k';'c'];
%     p(1)=plot(ax,(1:lwlet)/sps,greenplt(:,1),'r-','linew',1);
%     p(2)=plot(ax,(1:lwlet)/sps,greenplt(:,2),'b-','linew',1);
%     p(3)=plot(ax,(1:lwlet)/sps,greenplt(:,3),'k-','linew',1);
%     p(4)=plot(ax,(1:lwlet)/sps,greenplt(:,4),'c-','linew',1);
%   else
%     color=jet(nsta);
%   end
%   for ista = 1: nsta
%     p(ista)=plot(ax,(1:lwlet)/sps,greenplt(:,ista),'-','color',color(ista,:),'linew',1);
%   end
  p(1)=plot(ax,(1:lwlet)/sps,greenf(:,i),'-','color',color(1,:),'linew',1);
  p(2)=plot(ax,(1:lwlet)/sps,greenfort(:,i),'-','color',color(2,:),'linew',1);
  p(3)=plot(ax,(1:lwlet)/sps,greenfvert(:,i),'-','color',color(3,:),'linew',1);
  text(ax,0.98,0.9,sprintf('%.1f-%.1f Hz, %s',lowlet,hiwlet,stas(staplt(i),:)),...
    'Units','normalized','HorizontalAlignment','right','fontsize',9);
%   text(ax,0.01,0.85,panelname{i},'FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
%   xlim(ax,[0 lwlet]/sps);
  xlim(ax,[0 12]);
%   if i == 1
    ylim(ax,[-mx mx]);
%   else
%     ylim(ax,[-mx mx]/4);
%   end
  if i<nsta
    nolabels(ax,1);
  end
  
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
    
    text(ax,0.01+(1-1)*0.1,0.9,sprintf('%.1f;',snrf(:,i)),'Units',...
      'normalized','FontSize',9,'color',color(1,:));
    text(ax,0.01+(2-1)*0.1,0.9,sprintf('%.1f;',snrfort(:,i)),'Units',...
      'normalized','FontSize',9,'color',color(2,:));
    text(ax,0.01+(3-1)*0.1,0.9,sprintf('%.1f;',snrfvert(:,i)),'Units',...
      'normalized','FontSize',9,'color',color(3,:));
    
    
  end
  
%   ax.YColor=ax.XColor;
%   yyaxis(ax,'right');
%   %%%running CC using a window length of 'cclen'
%   mwlen=sps/2;
%   [ircc,rcc12] = RunningCC(greenplt(:,1), greenplt(:,2), mwlen);
%   [~,rcc13] = RunningCC(greenplt(:,1), greenplt(:,3), mwlen);
%   [~,rcc23] = RunningCC(greenplt(:,2), greenplt(:,3), mwlen);
%   rcc = (rcc12+rcc13)/3;
%   %alln(alln<0)=-10^4*yma; %just so they don't plot.
%   % plot(samples(cclen/2+1:greenlen-cclen/2),mx*alln,'co','markersize',2); % scale with data amp.
% %   plot(ircc/sps,mx*rcc,'co','markersize',2); % scale with data amp.
%   plot(ax,ircc/sps,rcc,'o','color',[.8 .8 .8],'markersize',1);
%   ylabel(ax,'Running CC','FontSize',10);
%   ylim(ax,[-1 1]);
% %   nolabels(ax,1);
%   ax.YColor=[.6 .6 .6];
  
  longticks(ax,6);
  hold(ax,'off');
  
end

legend(f.ax(1),p,compname,'Location','southeast','Orientation','horizontal');
xlabel(f.ax((nrow-1)*ncol+1),'Time (s)');
ylabel(f.ax((nrow-1)*ncol+1),'Amplitude');

delete(f.ax(nsta+1:end));






