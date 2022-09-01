function [f] = plt_groupimptriplets(sigdecon,impindep,stas,ircc,rcc)
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

f.fig=figure;
f.fig.Renderer = 'painters';
color = {'r','b','k'};

nsta = size(sigdecon,2);
lsig = size(sigdecon,1);
for ista = 1: nsta
  subplot(4,1,ista)
  ax = gca;
  hold(ax,'on');
  ax.Box = 'on';
  grid(ax,'on');
  indtemp = find(sigdecon(:,ista)>0);  % independent impulses from other stations
  stem(ax,indtemp,sigdecon(indtemp,ista),'color',[.7 .7 .7],'MarkerSize',4);
  text(ax,0.9,0.9,num2str(length(indtemp)),'unit','normalized','HorizontalAlignment','right',...
    'FontSize',10,'color',[.8 .8 .8]);
  indtemp = find(impindep(:,(ista-1)*2+1)>0);  % that can be individually paired with 3-sta among independent
  stem(ax,impindep(indtemp,(ista-1)*2+1),impindep(indtemp,(ista-1)*2+2),'color',color{ista},...
    'MarkerSize',4);
  text(ax,0.9,0.75,num2str(length(indtemp)),'unit','normalized','HorizontalAlignment','right',...
    'FontSize',10,'color',color{ista});
  text(ax,0.05,0.9,strcat(stas(ista,:)),'unit','normalized','HorizontalAlignment','left',...
    'FontSize',12);
  xlim(ax,[0 lsig]); hold(ax,'off');
  if ista == 1
    title('Grouped impulse triplets from independent deconvolution');
  end
end
subplot(4,1,4)
ax = gca;
yyaxis(ax,'left');
hold(ax,'on');
ax.Box = 'on';
grid(ax,'on');
stem(ax,impindep(:,1),impindep(:,2), 'r-','MarkerSize',4);
stem(ax,impindep(:,3),impindep(:,4), 'b-','MarkerSize',4);
stem(ax,impindep(:,5),impindep(:,6), 'k-','MarkerSize',4);
text(ax,0.05,0.9,sprintf('Number of triplets: %d', size(impindep,1)),'fontsize',10,...
  'Units','normalized');
ylabel('Amplitude');
xlabel('Samples');
yyaxis(ax,'right');
plot(ircc,rcc,'o','color',[.5 .5 .5],'markersize',2);
ylabel('Running CC');
ylim(ax,[-1 1]);
xlim(ax,[0 lsig]);
hold(ax,'off');




