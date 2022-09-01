function [f] = plt_cumulative_density(hfplt,lfplt,xran,yran,binmethod,msizehf,dxhf,dyhf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [f] = plt_cumulative_density(hfplt,lfplt,xran,yran)
% This function is to plot the cumulative density map of the hf and lf 
% catalog to see the activated region of the catalogs
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/09/28
% Last modified date:   2021/09/28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('hfplt',[]);
defval('lfplt',[]);
defval('xran',[-20 20]);
defval('yran',[-20 20]);
defval('binmethod','grid');
defval('dxhf',0.2);
defval('dyhf',0.2);
if isequal(binmethod,'grid')
  defval('msizehf',10);
elseif isequal(binmethod,'pixel')
  defval('msizehf',6);
end

[scrsz, res] = pixelperinch(1);

f.fig=figure;
f.fig.Renderer='Painters';
nrow = 1;
ncol = (~isempty(hfplt))+(~isempty(lfplt));
widin = 4*ncol;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

pltxran = [0.08 0.95]; pltyran = [0.15 0.9];
pltxsep = 0.08; pltysep = 0.05; 
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

for isub = 1:nrow*ncol
  set(f.ax(isub), 'position', axpos(isub,:));
end    

% xran = [-20 25];
% yran = [-20 20];

% subplot 1 of figure i
% hfplt = hftime(hftime(:,13) > 2004*1000 & hftime(:,13) < 2005*1000, :);
% lfplt = lftime(lftime(:,13) > 2004*1000 & lftime(:,13) < 2005*1000, :);

ax = f.ax(1);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
if ~isempty(hfplt)
  %create a density matrix to store the number of detections in each small grid
  if isequal(binmethod,'grid') 
%     dxhf = 0.2;
%     dyhf = 0.2;
    [density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
      xran,yran,dxhf,dyhf);
    marker = 'o';
%     msizehf = 500*dxhf*dyhf;
    msizelf = 2.5*msizehf;
  %bin based upon pixel   
  elseif isequal(binmethod,'pixel') 
    density1d = density_pixel(hfplt(:,1),hfplt(:,2));
    marker = 'o';
%     msizehf = 6;
    msizelf = 2.5*msizehf;
  end
  dumhf = density1d(density1d(:,3)>0, :);
  dumhf(dumhf(:,3)>1, :) = [];
  scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),marker,'linew',0.2);  %, 'MarkerEdgeColor', 'w')
  dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
  dumhf(dumhf(:,3)==1, :) = [];
  scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),marker,'filled','MarkerEdgeColor','none');  %, 
  % imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
  oldcmap = colormap(ax,'jet');
  % colormap(ax, flipud(oldcmap) );
  c=colorbar(ax,'SouthOutside');
  pos = ax.Position;
  c.Position = [pos(1), pos(2)-0.03, pos(3), 0.02];
  c.Label.String = strcat({'log_{10}(# tremor detections / '},binmethod,')');
  c.Label.FontSize = 10;
% text(ax, 0.85, 0.93, '2004','FontSize',12,'unit','normalized','horizontalalignment','center',...
%      'EdgeColor','k','Margin',2);
% text(ax,0.85,0.8,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
  text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
  text(ax,0.5,0.93,strcat(num2str(length(hfplt(:,1))),{' detections'}),'FontSize',10,'unit','normalized',...
      'horizontalalignment','center');
end
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
% xticks(ax,xran(1):5:xran(2));
% yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

% subplot 2 of figure i
if ~isempty(lfplt)
  ax = f.ax(2);
  hold(ax,'on');
  plot(ax,[-100 100],[0 0],'k--');
  plot(ax,[0 0],[-100 100],'k--');
  ax.FontSize = 9;
  %create a density matrix to store the number of detections in each small grid
  if isequal(binmethod,'grid') 
    dxlf = 2.5*dxhf;
    dylf = 2.5*dyhf;
    [density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
      xran,yran,dxlf,dylf);
  %bin based upon pixel   
  elseif isequal(binmethod,'pixel')
    density1d = density_pixel(lfplt(:,1),lfplt(:,2));
  end
  
  dumhf = density1d(density1d(:,3)>0, :);
  dumhf(dumhf(:,3)>1, :) = [];
  scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)),marker,'linew',0.2);  %, 'MarkerEdgeColor', 'w')
  dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
  dumhf(dumhf(:,3)==1, :) = [];
  scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)),marker,'filled');  %, 'MarkerEdgeColor', 'w')
  % imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
  oldcmap = colormap(ax,'jet');
  % colormap(ax, flipud(oldcmap) );
  c=colorbar(ax,'SouthOutside');
  pos = ax.Position;
  c.Position = [pos(1), pos(2)-0.03, pos(3), 0.02];
%   if isequal(binmethod,'grid') 
    c.Label.String = strcat({'log_{10}(# tremor detections / '},binmethod,')');
%   elseif isequal(binmethod,'pixel')
%     
%   end
  c.Label.FontSize = 10;
% text(ax, 0.85, 0.93, '2004','FontSize',12,'unit','normalized','horizontalalignment','center',...
%      'EdgeColor','k','Margin',2);
% text(ax,0.85,0.8,'LF','FontSize',12,'unit','normalized','horizontalalignment','center');
  text(ax,0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
  text(ax,0.5,0.93,strcat(num2str(length(lfplt(:,1))),{' detections'}),'FontSize',10,'unit','normalized',...
      'horizontalalignment','center');
  ax.Box = 'on';
  grid(ax, 'on');
  axis(ax, 'equal');
  ax.GridLineStyle = '--';
  ax.XAxisLocation = 'top';
  xlim(ax,xran);
  ylim(ax,yran);
%   xticks(ax,xran(1):5:xran(2));
%   yticks(ax,yran(1):5:yran(2));
  xlabel(ax,'E (km)','fontsize',11);
  ylabel(ax,'N (km)','fontsize',11);
  hold(ax,'off');
end
