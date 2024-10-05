function [f] = plt_secondarysrc(ppkindep,sigsta,ppkindepsave,indremove,stas,...
  xzoom,sps,normflag,detecttype,saveflag)
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
defval('xzoom',[]);
defval('sps',160);
defval('normflag',0);  %default is do not normalize the seismograms
defval('detecttype',[]);  %default is short-win detections
defval('saveflag',1);

nrow = 3;
ncol = 1;
widin = 8.3;  % maximum width allowed is 8.5 inches
htin = 4;   % maximum height allowed is 11 inches
pltxran = [0.05 0.98]; pltyran = [0.1 0.98]; % optimal axis location
pltxsep = 0.05; pltysep = 0.025;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

color = ['r';'b';'k'];
nsta = size(sigsta,2);
lsig = size(sigsta,1);

%if need to scale the seismograms
if normflag
%   maxval = max(sigsta);
  maxval = range(sigsta(xzoom(1)*sps+1:xzoom(2)*sps, :));
  sigsta = sigsta ./ maxval;
  % ppkindep(:,[2 4 6]) = ppkindep(:,[2 4 6]) ./ normalizer;
  % ppkindepsave(:,[2 4 6]) = ppkindepsave(:,[2 4 6]) ./ normalizer;
end
scale = 2;

for ista = 1: nsta
  ax = f.ax(ista); hold(ax,'on'); ax.Box = 'on';
  stem(ax,ppkindep(indremove,(ista-1)*2+1)/sps,ppkindep(indremove,ista*2),...
    'Color',[.6 .6 .6],'MarkerSize',2.5,'linew',0.8);
  stem(ax,ppkindepsave(:,(ista-1)*2+1)/sps,ppkindepsave(:,ista*2),...
    'Color',color(ista,:),'MarkerSize',3,'linew',0.8);
  plot(ax,(1:lsig)/sps,sigsta(:,ista)*scale+3.0, 'k-','linew',1);
  text(ax,0.99,0.9,strcat(stas(ista,:)),'unit','normalized',...
    'HorizontalAlignment','right','FontSize',12);
  ylim(ax,[0 4.5]);
  %       xlim(ax,[0 lsig]/sps);
  xlim(ax,xzoom);
  longticks(ax,6);
  noticks(ax,2);
  if ista <3
    nolabels(ax,3);
  else
    xlabel(ax,'Time (s)','FontSize',10);
  end
  ylabel(ax,'Norm. amplitude','FontSize',10);
  hold(ax,'off');
end

if saveflag
  % orient(f.fig,'landscape');
  fname = strcat('remove2ndarysrc',detecttype,'.pdf');
  print(f.fig,'-dpdf',...
    strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
end
