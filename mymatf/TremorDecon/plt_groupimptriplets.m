function [f] = plt_groupimptriplets(sigdecon,impindep,stas,ircc,rcc,renv,sps,xzoom,detecttype,saveflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_groupimptriplets(sigdecon,impindep,stas,ircc,rcc)
%
% plot the individually deconvolved impulses at each station in gray, and 
% the grouped impulse triplets in differernt colors. A summary of grouped
% result only comes afterwards.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/06/21
% Last modified date:   2022/06/21 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('sps',160);
defval('xzoom',[]);
defval('detecttype',[]);  %default is short-win detections
defval('saveflag',1);

% widin = 7.5;
widin = 8.3;
htin = 5;
nrow = 4;
ncol = 1;
pltxran = [0.04 0.94]; pltyran = [0.08 0.98]; % optimal axis location
pltxsep = 0.05; pltysep = 0.02;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

color = ['r';'b';'k'];

nsta = size(sigdecon,2);
lsig = size(sigdecon,1);
nsrc = size(impindep,1);
color2 = jet(nsrc);
% color2 = kelicmap(nsrc);
% color2 = brewermap(nsrc,'reds');

for ista = 1: nsta
  ax = f.ax(ista); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'off');
  indtemp = find(sigdecon(:,ista)>0);  % independent impulses from other stations
  stem(ax,indtemp/sps,sigdecon(indtemp,ista),'color',[.7 .7 .7],'MarkerSize',3);
  indzoom = find(indtemp/sps>=xzoom(1) & indtemp/sps<=xzoom(2));
  text(ax,0.05,0.9,num2str(length(indzoom)),'unit','normalized','HorizontalAlignment','left',...
    'FontSize',10,'color',[.7 .7 .7]);
  indtemp = find(impindep(:,(ista-1)*2+1)>0);  % that can be individually paired with 3-sta among independent
  for i = 1: nsrc
    stem(ax,impindep(i,(ista-1)*2+1)/sps,impindep(i,(ista-1)*2+2),...
      'color',color2(i,:),'MarkerFaceColor',color2(i,:),'MarkerSize',4);
  end
%   stem(ax,impindep(indtemp,(ista-1)*2+1),impindep(indtemp,(ista-1)*2+2),'color',color{ista},...
%     'MarkerSize',4);
  indzoom = find(impindep(:,(ista-1)*2+1)/sps>=xzoom(1) & impindep(:,(ista-1)*2+1)/sps<=xzoom(2));
  text(ax,0.05,0.75,num2str(length(indzoom)),'unit','normalized','HorizontalAlignment','left',...
    'FontSize',10,'color',color(ista,:));
  text(ax,0.98,0.9,strcat(stas(ista,:)),'unit','normalized','HorizontalAlignment','right',...
    'FontSize',12);
  if isempty(xzoom)
    xlim(ax,[0 lsig]/sps);
  else
    xlim(ax,xzoom);
  end
  ylim(ax,[0 4.5]);
%   if ista == 1
%     title(ax,'Grouped impulse triplets from independent deconvolution');
%   end
  ylabel(ax,'Amplitude','FontSize',10);
  longticks(ax,6);
  nolabels(ax,1);
  hold(ax,'off');
end
text(f.ax(1),0.01,0.85,'a','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
text(f.ax(2),0.01,0.85,'b','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
text(f.ax(3),0.01,0.85,'c','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);

ax = f.ax(4); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'off');
yyaxis(ax,'right');
% plot(ax,ircc/sps,rcc,'o','color',[.8 .8 .8],'markersize',1);
scatter(ax,ircc/sps,rcc,2,renv,'o','filled');
ind=find(ircc>=xzoom(1)*sps & ircc<=xzoom(2)*sps);
colormap(ax,'jet');
cran = [prctile(renv(ind),1) prctile(renv(ind),99)];
caxis(ax,cran);
ylabel(ax,'Running CC','FontSize',10);
ylim(ax,[-1 1]);
if isempty(xzoom)
  xlim(ax,[0 lsig]/sps);
else
  xlim(ax,xzoom);
end
ax.YColor=[.6 .6 .6];

yyaxis(ax,'left');
p(1)=stem(ax,impindep(:,1)/sps,impindep(:,2), 'r-','MarkerSize',3);
p(2)=stem(ax,impindep(:,3)/sps,impindep(:,4), 'b-','MarkerSize',3);
p(3)=stem(ax,impindep(:,5)/sps,impindep(:,6), 'k-','MarkerSize',3);
% text(ax,0.05,0.9,sprintf('Number of triplets: %d', size(impindep,1)),'fontsize',10,...
%   'Units','normalized');
text(f.ax(4),0.01,0.85,'d','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
ylabel(ax,'Amplitude','FontSize',10);
xlabel(ax,'Time (s)','FontSize',10);
ylim(ax,[0 4.5]);
ax.YColor=ax.XColor;
legend(ax,p,stas(1:nsta,:),'Orientation','horizontal','Location','north');
longticks(ax,6);
hold(ax,'off');

if saveflag
  % orient(f.fig,'landscape');
  fname = strcat('groupedimpulses',detecttype,'.pdf');
  print(f.fig,'-dpdf',...
    strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
end

% keyboard
