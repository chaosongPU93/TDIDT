function [f] = plt_sum_pixel(density1d,sumz1d,xran,yran,msize,cstr,scale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [sumz1d,indices] = sum_pixel(x,y,z)
%
% Beyond 'density_pixel', sometimes you not only want the cumulative count
% of data (x,y) at each pixel, but also you want the sum of z at those 
% unique points. The sum of z and indices of unique pixels are returned.
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2023/05/27
% Last modified date:   2023/05/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('scale','log');

widin = 12;  % maximum width allowed is 8.5 inches
htin = 6;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f = initfig(widin,htin,nrow,ncol); %initialize fig

%cumulative density
ax=f.ax(1);
hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
dum = density1d;
dum(dum(:,3)>1, :) = [];
if strcmp(scale,'log')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),'o','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dum = sortrows(density1d,3);
dum(dum(:,3)==1, :) = [];
if strcmp(scale,'log')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),'o','filled','MarkerEdgeColor','none');
oldc = colormap(ax,'kelicol');
newc = flipud(oldc);
colormap(ax,newc);
c=colorbar(ax,'SouthOutside');
if strcmp(scale,'log')
  c.Label.String = strcat('log_{10}(# detections / pixel)');
elseif strcmp(scale,'linear')
  c.Label.String = strcat('# detections / pixel');  
end
axis(ax,[xran yran],'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';

%sum of some quantity
ax=f.ax(2);
hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
sumz1d = sortrows(sumz1d,3);
if strcmp(scale,'log')
  sumz1d(:,3) = log10(sumz1d(:,3));
end
scatter(ax,sumz1d(:,1),sumz1d(:,2),msize,sumz1d(:,3),'o','filled','MarkerEdgeColor','none');
oldc = colormap(ax,'kelicol');
newc = flipud(oldc);
colormap(ax,newc);
c=colorbar(ax,'SouthOutside');
if strcmp(scale,'log')
  c.Label.String = strcat({'log_{10}('},cstr,' / pixel)');
elseif strcmp(scale,'linear')
  c.Label.String = strcat(cstr,' / pixel'); 
end
axis(ax,[xran yran],'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';


