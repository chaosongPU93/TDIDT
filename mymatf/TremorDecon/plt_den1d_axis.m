function [ax,den1d,conmat,conobj] = plt_den1d_axis(ax,den1d,conmat,xran,yran,...
  dx,dy,binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ax,den1d,conmat,conobj] = plt_den1d_axis(ax,den1d,conmat,...
%   xran,yran,dx,dy,binmethod,marker,msize,disttype,contourflag,smoothsigma,scale)
% 
% Given an axis of a figure, this function is to plot the cumulative density 
% matrix 'den1d' that contains 3 cols as locations in x and y directions and
% count at that location. It is similar to 'plt_cumulative_density_axis', but
% this function directly deals with a density matrix rather than raw locations.
% See also 'plt_cumulative_density', 'plt_srcdlocinmap', 'plt_4stremordensity'.
%
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/03/22
% Last modified date:   2024/03/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('xran',[-4 4]);
defval('yran',[-4 4]);
defval('binmethod','pixel');
defval('marker','o');
defval('msize',6);
defval('dx',[]);
defval('dy',[]);
defval('disttype','km');
defval('contourflag',0);
defval('smoothsigma',1);
defval('scale','log10');

hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
ax.FontSize = 9;
xlim(ax,xran);
ylim(ax,yran);
% xticks(ax,xran(1):5:xran(2));
% yticks(ax,yran(1):5:yran(2));
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
%scatter the density
dum = den1d;
dum(dum(:,3)>1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'linew',0.2);  %, 'MarkerEdgeColor', 'w')
dum = sortrows(den1d,3);
dum(dum(:,3)==1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'filled','MarkerEdgeColor','none');
text(ax,0.98,0.05,sprintf('%d events',sum(den1d(:,3))),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
colormap(ax,'plasma');
% colormap(ax,'jet');
% colormap(ax,flipud(colormap(ax,'kelicol')));
c=colorbar(ax,'SouthOutside');
if strcmp(scale,'log10')
  cstr = strcat({'log_{10}(# events / '},binmethod,')');
elseif strcmp(scale,'linear')
  cstr = strcat({'# detections / '},binmethod);
end
c.Label.String = cstr;
%plot contour lines
if contourflag && isempty(conmat)
  %compute contour lines, make it dense
  if strcmp(disttype,'spl')
    xvec = min(den1d(:,1)):1:max(den1d(:,1));
    yvec = min(den1d(:,2)):1:max(den1d(:,2));
  elseif strcmp(disttype,'km')
    xvec = floor(min(den1d(:,1))):dx:ceil(max(den1d(:,1)));
    yvec = floor(min(den1d(:,2))):dy:ceil(max(den1d(:,2)));
  % [~,xgrid,ygrid,zgrid] = zeropadmat2d(den1d,xvec,yvec);
  end
  [~,xgrid,ygrid,zgrid] = zeropadmat2d(den1d,xvec,yvec);

  %smoothing with a Gaussian filter
  if ~isempty(smoothsigma)
    zgridgf = imgaussfilt(zgrid,smoothsigma);
  else
    zgridgf = zgrid;
  end
  if strcmp(scale,'log10')
    zgridgf = log10(zgridgf);
  end

%   perc = 50:10:90;
%   conplt = prctile(dum(:,3),perc);  
  if strcmp(disttype,'spl')
    % [conmat,conobj] = contour(ax,xgrid/sps,ygrid/sps,zgridgf,conplt,'-','color',[.2 .2 .2]);
    [conmat,conobj] = contour(ax,xgrid/sps,ygrid/sps,zgridgf,ncont,'-','color',[.2 .2 .2]);
  %   conmat = contourc(zgridgf,ncont);
  elseif strcmp(disttype,'km')
    % [conmat,conobj] = contour(ax,xgrid,ygrid,zgridgf,conplt,'-','color',[.2 .2 .2]);
    [conmat,conobj] = contour(ax,xgrid,ygrid,zgridgf,ncont,'-','color',[.2 .2 .2]);
  %   conmat = contourc(xgrid,ygrid,zgridgf,ncont);
  end
  if ~isempty(conmat)
    contable = getContourLineCoordinates(conmat);
    conmat=table2array(contable);
  end
  delete(conobj);
  
% elseif contourflag && ~isempty(conmat)
%   grp = unique(conmat(:,2));
%   ngrp = max(grp);
%   contval = unique(conmat(:,1));
%   ncont = length(contval);
%   for icont = 1:ncont
%   %   contval=contintvl*icont;
%     for igrp = 1:ngrp
%       ind = find(conmat(:,1)==contval(icont) & conmat(:,2)==igrp);
%       plot(ax,conmat(ind,3),conmat(ind,4),'-','color',[.2 .2 .2],'linew',0.5);
%     end
%   end
  
end

%plot contour lines
if contourflag && ~isempty(conmat)
  ncontplt = 5;
  contst = ncont*0.75;
  contintvl = floor((ncont-contst)/ncontplt);
  contval=unique(conmat(:,1),'stable');
  for icont = contst: contintvl: ncont-contintvl
      tmp = conmat(conmat(:,1)==contval(icont),:);
      grp = unique(tmp(:,2));
      ngrp = max(grp);
      for igrp = 1:ngrp
        ind = find(tmp(:,2)==igrp);
        plot(ax,tmp(ind,3),tmp(ind,4),'-','color',[.2 .2 .2],'linew',0.5);
      end
  end
end


if strcmp(disttype,'spl')
  xlabel(ax,sprintf('\\Delta{t}_{12} (s)'));
  ylabel(ax,sprintf('\\Delta{t}_{13} (s)'));
elseif strcmp(disttype,'km')
  xlabel(ax,'E (km)');
  ylabel(ax,'N (km)');
end
longticks(ax,2);
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.07, pos(3), 0.02];
% hold(ax,'off');

