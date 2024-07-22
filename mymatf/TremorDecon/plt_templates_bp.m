function f=plt_templates_bp(greenf,greenfort,greenfvert,stas,staplt,...
  lowlet,hiwlet,sps,patchwins)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_templates_bp(greenf,greenfort,greenfvert,stas,staplt,...
%   lowlet,hiwlet,sps,patchwins)
%
% this function simply plot the bandpassed filtered templates (or green's 
% functions), either optimal component or orthogonal comp.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/09/12
% Last modified date:   2022/09/12 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('patchwins',[]);

greenall{1} = greenf;
compname{1} = 'Optimal';
panelname{1} = 'a';
if isempty(greenfort) && isempty(greenfvert)
  nrow = 1;
elseif ~isempty(greenfort) && isempty(greenfvert)
  nrow = 2;
  greenall{2} = greenfort;
  compname{2} = 'Orthogonal';
  panelname{2} = 'b';
elseif isempty(greenfort) && ~isempty(greenfvert)
  nrow = 2;
  greenall{2} = greenfvert;
  compname{2} = 'Vertical';
  panelname{2} = 'b';
elseif ~isempty(greenfort) && ~isempty(greenfvert)
  nrow = 3;
  greenall{2} = greenfort;
  greenall{3} = greenfvert;
  compname{2} = 'Orthogonal';
  compname{3} = 'Vertical';
  panelname{2} = 'b';
  panelname{3} = 'c';
end

widin = 6;  % maximum width allowed is 8.5 inches
htin = 1.5*nrow;   % maximum height allowed is 11 inches
ncol = 1;
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.1 0.95]; pltyran = [0.1 0.98]; % optimal axis location
pltxsep = 0.05; pltysep = 0.015;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

greenplt = greenall{1};
mx=1.2*max(max(abs(greenplt(:,:))));

for i = 1: nrow
  greenplt = greenall{i};
  greenplt = greenplt(:,staplt);
  ax = f.ax(i);
  hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%   yyaxis(ax,'left');
  lwlet = size(greenplt,1);
  nsta = size(greenplt,2);
  if nsta == 3
    color=['r';'b';'k'];
%     p(1)=plot(ax,(1:lwlet)/sps,greenplt(:,1),'r-','linew',1);
%     p(2)=plot(ax,(1:lwlet)/sps,greenplt(:,2),'b-','linew',1);
%     p(3)=plot(ax,(1:lwlet)/sps,greenplt(:,3),'k-','linew',1);
  elseif nsta == 4
    color=['r';'b';'k';'c'];
%     p(1)=plot(ax,(1:lwlet)/sps,greenplt(:,1),'r-','linew',1);
%     p(2)=plot(ax,(1:lwlet)/sps,greenplt(:,2),'b-','linew',1);
%     p(3)=plot(ax,(1:lwlet)/sps,greenplt(:,3),'k-','linew',1);
%     p(4)=plot(ax,(1:lwlet)/sps,greenplt(:,4),'c-','linew',1);
  else
    color=jet(nsta);
  end
  for ista = 1: nsta
    p(ista)=plot(ax,(1:lwlet)/sps,greenplt(:,ista),'-','color',color(ista,:),'linew',1);
  end
  text(ax,0.98,0.9,sprintf('%.1f-%.1f Hz, %s',lowlet,hiwlet,compname{i}),'Units','normalized',...
    'HorizontalAlignment','right','fontsize',9);
  text(ax,0.01,0.85,panelname{i},'FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
  greenlen = size(greenplt,1);
  xlim(ax,[0 greenlen]/sps);
  if i == 1
    ylim(ax,[-mx mx]);
  else
    ylim(ax,[-mx mx]/4);
  end
  if i<nrow
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
    
    envf = envelope(detrend(greenplt));
    snrf = range(envf(winsig(1)*sps: winsig(2)*sps, :),1)./ ...
      median(envf(winnoi(1)*sps: winnoi(2)*sps, :),1);
    for ista = 1: size(greenplt,2)
      text(ax,0.01+(ista-1)*0.07,0.1,sprintf('%.1f;',snrf(:,ista)),'Units',...
        'normalized','FontSize',9,'color',color(ista,:));
    end    
    
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

legend(f.ax(1),p,stas(staplt,:),'Location','southeast','Orientation','horizontal');
xlabel(f.ax(nrow),'Time (s)');
ylabel(f.ax(nrow),'Amplitude');








