function [f,tsep,pkist,indpk,indreverse] = plt_tsep_deconpk(pkindepsave,sps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = comp_deconpk_sigpk(sigsta,pkindep,indremove,fpltremv)
%
% Plot the separation in arrival time between preserved positive peaks after 
% removing the secondary ones, to see if they can be too close to each other
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/06/21
% Last modified date:   2022/06/21 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pkist = [];
indpk = [];
tsep = [];
nsta = 3;
for ista = 1:nsta
  pki = pkindepsave(:,(ista-1)*2+1);
  [pkist(:,ista),indpk(:,ista)] = sort(pki,'ascend');
  tsep(:,ista) = diff(pkist(:,ista))/sps;
end
mintsep = min(tsep,[],1);

indreverse = find(indpk(:,1)~=indpk(:,2) | indpk(:,1)~=indpk(:,3) | indpk(:,2)~=indpk(:,3));

f.fig = figure;
f.fig.Renderer = 'painters';
axall(1) = subplot(1,2,1);
ax = gca;
hold(ax,'on');
ax.Box = 'on';
grid(ax,'on');

h1=histogram(ax,tsep(:,1),'FaceColor','r','BinWidth',0.1);
h2=histogram(ax,tsep(:,2),'FaceColor','b','BinWidth',0.1);
h3=histogram(ax,tsep(:,3),'FaceColor','k','BinWidth',0.1);
% text(ax,0.95,0.8,sprintf('%d; %d; %d',mintsep(1),mintsep(2),mintsep(3)),'Units','normalized',...
%   'HorizontalAlignment','right');
mbincnt1 = h1.BinEdges(h1.BinCounts == max(h1.BinCounts))+(h1.BinWidth)/2;
mbincnt2 = h2.BinEdges(h2.BinCounts == max(h2.BinCounts))+(h2.BinWidth)/2;
mbincnt3 = h3.BinEdges(h3.BinCounts == max(h3.BinCounts))+(h3.BinWidth)/2;
text(ax,0.95,0.6,sprintf('%.2f',mbincnt1),'Units','normalized',...
  'HorizontalAlignment','right');
text(ax,0.95,0.5,sprintf('%.2f',mbincnt2),'Units','normalized',...
  'HorizontalAlignment','right');
text(ax,0.95,0.4,sprintf('%.2f',mbincnt3),'Units','normalized',...
  'HorizontalAlignment','right');
text(ax,0.95,0.9,sprintf('%d',size(pkindepsave,1)),'Units','normalized',...
  'HorizontalAlignment','right');
xlim(ax,[0 2]);
legend(ax,'PGC','SSIB','SILB','Location','north');
xlabel(ax,'Separation time (s)');
ylabel(ax,'Counts');
hold(ax,'off');

axall(2) = subplot(1,2,2);
ax = gca;
hold(ax,'on');
ax.Box = 'on';
grid(ax,'on');
[cdfval,x] = ecdf(tsep(:,1));
plot(ax,x,cdfval,'r','linew',1);
[cdfval,x] = ecdf(tsep(:,2));
plot(ax,x,cdfval,'b','linew',1);
[cdfval,x] = ecdf(tsep(:,3));
plot(ax,x,cdfval,'k','linew',1);
text(ax,0.95,0.25,sprintf('25 perc: %.2f; %.2f; %.2f',prctile(tsep(:,1),25),prctile(tsep(:,2),25),...
  prctile(tsep(:,3),25)),'Units','normalized','HorizontalAlignment','right');
text(ax,0.95,0.15,sprintf('10 perc: %.2f; %.2f; %.2f',prctile(tsep(:,1),10),prctile(tsep(:,2),10),...
  prctile(tsep(:,3),10)),'Units','normalized','HorizontalAlignment','right');
text(ax,0.95,0.05,sprintf('5 perc: %.2f; %.2f; %.2f',prctile(tsep(:,1),5),prctile(tsep(:,2),5),...
  prctile(tsep(:,3),5)),'Units','normalized','HorizontalAlignment','right');
xlim(ax,[0 2]);
legend(ax,'PGC','SSIB','SILB','Location','northwest');
xlabel(ax,'Separation time (s)');
ylabel(ax,'Empirical CDF');
hold(ax,'off');
      
h=supertit(axall,'Separation time between consecutive positive peaks in arrival',10);
movev(h,0.2);
%       title('imp. must point to diff. waveform peaks for >=1 sta');
% title(ax,'imp. must point to diff. waveform peaks for all 3 stas');
%       title('imp separated >= 20 samples for 3 sta');
      
      
