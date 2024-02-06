function [f] = plt_sum_pixel(density1d,sumz1d,xran,yran,msize,cstr,symbol,scale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_sum_pixel(density1d,sumz1d,xran,yran,msize,cstr,symbol,scale)
%
% Regardless of binning method, plot the cumulative density and related 
% quantity that is summed at indices with multiple detections at single 
% bins. The 'scale' of data can be either 'linear' or 'log10'.
% 
%  
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2023/05/27
% Last modified date:   2023/05/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('symbol','o');
defval('scale','log10');

htin = 5;   % maximum height allowed is 11 inches
nrow = 1;
if ~isempty(sumz1d)
  widin = 10;  % maximum width allowed is 8.5 inches
  ncol = 2;
else
  widin = 5;  % maximum width allowed is 8.5 inches
  ncol = 1;
end
f = initfig(widin,htin,nrow,ncol); %initialize fig

%cumulative density
ax=f.ax(1);
hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
dum = density1d;
dum(dum(:,3)>1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),symbol,'linew',0.2);  %, 'MarkerEdgeColor', 'w')
dum = sortrows(density1d,3);
dum(dum(:,3)==1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),symbol,'filled','MarkerEdgeColor','none');
text(ax,0.98,0.05,sprintf('%d events',sum(density1d(:,3))),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
oldc = colormap(ax,'kelicol');
newc = flipud(oldc);
colormap(ax,newc);
c=colorbar(ax,'SouthOutside');
ax.CLim(2) = prctile(dum(:,3),99);
if strcmp(scale,'log10')
  c.Label.String = strcat('log_{10}(',cstr{1},')');
elseif strcmp(scale,'linear')
  c.Label.String = cstr{1};  
end
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

if ~isempty(sumz1d)
  %sum of some quantity
  ax=f.ax(2);
  hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
  sumz1d = sortrows(sumz1d,3);
  if strcmp(scale,'log10')
    sumz1d(:,3) = log10(sumz1d(:,3));
  end
  scatter(ax,sumz1d(:,1),sumz1d(:,2),msize,sumz1d(:,3),symbol,'filled','MarkerEdgeColor','none');
  oldc = colormap(ax,'kelicol');
  newc = flipud(oldc);
  colormap(ax,newc);
  c=colorbar(ax,'SouthOutside');
%   ax.CLim(2) = prctile(sumz1d(:,3),99);
  caxis(ax,[prctile(sumz1d(:,3),1) prctile(sumz1d(:,3),99)]);
  if strcmp(scale,'log10')
    c.Label.String = strcat('log_{10}(',cstr{2},')');
  elseif strcmp(scale,'linear')
    c.Label.String = cstr{2};  
  end
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
end

