function [ax,den1d,conmat,angle,anglegeo,c,conmat10th,conobj10th] = plt_cumulative_density_axis(ax,locxy,xran,yran,...
  dx,dy,binmethod,marker,msize,disttype,smoothsigma,contourflag,ncont,scale,sps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax = plt_cumulative_density_axis(ax,locxy,xran,yran,binmethod,msize,...
%   dx,dy,disttype,contourflag,smoothsigma,scale)
% 
% Given an axis of a figure, this function is to plot the cumulative density 
% of a data array that contains 2 cols as locations in x and y directions.
% You can also choose to plot the contour lines at specified values (the default
% is to plot 5 contour lines at 50th to 90th percentiles of the density matrix).
% See also 'plt_cumulative_density', 'plt_srcdlocinmap', 'plt_4stremordensity'.
%
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/03/22
% Last modified date:   2024/03/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('xran',[-40 40]);
defval('yran',[-40 40]);
defval('binmethod','pixel');
defval('marker','o');
defval('msize',6);
defval('dx',[]);
defval('dy',[]);
defval('disttype','spl');
defval('contourflag',0);  %whether to PLOT a few contour lines
defval('smoothsigma',1);
defval('scale','log10');
defval('sps',160);

if strcmp(binmethod,'pixel')
  dx = 1; dy = 1; 
end

hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
ax.FontSize = 9;
xlim(ax,xran/sps);
ylim(ax,yran/sps);
% xticks(ax,xran(1):5:xran(2));
% yticks(ax,yran(1):5:yran(2));
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
if strcmp(binmethod,'pixel')
  den1d = density_pixel(locxy(:,1),locxy(:,2));
  cstr = '# events / pixel';
elseif strcmp(binmethod,'grid')
  den1d = density_matrix(locxy(:,1),locxy(:,2),xran,yran,dx,dy);
  cstr = '# events / grid';
end
den1d = den1d(den1d(:,3)>0, :);
den1d = sortrows(den1d,3);
%scatter the density
dum = den1d;
dum(dum(:,3)>1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
if strcmp(disttype,'spl')
  scatter(ax,dum(:,1)/sps,dum(:,2)/sps,msize,dum(:,3),marker,'linew',0.2);  %, 'MarkerEdgeColor', 'w')
elseif strcmp(disttype,'km')  
  scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'linew',0.2);  %, 'MarkerEdgeColor', 'w')
end
dum = sortrows(den1d,3);
dum(dum(:,3)==1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
if strcmp(disttype,'spl')
  scatter(ax,dum(:,1)/sps,dum(:,2)/sps,msize,dum(:,3),marker,'filled','MarkerEdgeColor','none');
elseif strcmp(disttype,'km')  
  scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'filled','MarkerEdgeColor','none');
end  
text(ax,0.98,0.05,sprintf('%d events',size(locxy,1)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9); %,'FontName','Monospaced'
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

%Principal component analysis
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(locxy);
% plot(ax,ellx/sps,elly/sps,'-','linew',1.5,'color','k');
% plot(ax,(x0+semib*[coeff(1,2),-coeff(1,2)])/sps, ...
%   (y0+semib*[coeff(2,2),-coeff(2,2)])/sps,'--','linew',1.5,'color','k');
% plot(ax,(x0+semia*[coeff(1,1),-coeff(1,1)])/sps, ...
%   (y0+semia*[coeff(2,1),-coeff(2,1)])/sps,'--','linew',1.5,'color','k');
% dxln = dx/sps/10; %in sec
% loff_max = 6/40;  %in sec
% x = reshape(-loff_max:dxln:loff_max, [], 1);
% % x = reshape(-0.1: dxln:0.1, [], 1);
% %a line cross (0,0) in the projection direction
% yopt = linefcn(x,tand(angle(2)),0);
% plot(ax,x,yopt,'-','linew',1.5,'color','k');
% %a line cross (0,0) in the orthogonal direction
% yort = linefcn(x,tand(angle(1)),0);
% plot(ax,x,yort,'--','linew',1.5,'color','k');

% text(ax,0.98,0.16,strcat(num2str(round(anglegeo(2))),'$^{\,\circ}$',{'; '},...
%   num2str(round(anglegeo(1))),'$^{\,\circ}$'),'FontSize',...
%   11,'unit','normalized','interpreter','latex','HorizontalAlignment','right');
%   text(ax,0.98,0.10,sprintf('%d^o',round(anglegeo(2))),'Units','normalized',...
%     'HorizontalAlignment','right','FontSize',9);
% text(ax,0.98,0.11,strcat({'Asp. ratio: '},sprintf('%.1f',semia/semib)),'Units',...
%   'normalized','HorizontalAlignment','right','FontSize',9);

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
zgridgf = imgaussfilt(zgrid,smoothsigma);
if strcmp(scale,'log10')
  zgridgf = log10(zgridgf);
end
% xvec = reshape(xgrid,[],1);
% yvec = reshape(ygrid,[],1);
% zvect = reshape(zgridgf,[],1);
% 
% perc = 50:10:90;
% conplt = prctile(dum(:,3),perc); 
dum = den1d;
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
%find the grid density that include a certain percentile of data
bb=cumsum(den1d(:,3))/size(locxy,1);
conplt5th = dum(abs(bb-0.05)==min(abs(bb-0.05)),3);
conplt10th = dum(abs(bb-0.1)==min(abs(bb-0.1)),3);
conpltprc = conplt10th;

if strcmp(disttype,'spl')
  % [conmat,conobj] = contour(ax,xgrid/sps,ygrid/sps,zgridgf,conplt,'-','color',[.2 .2 .2]);
  [conmat,conobj] = contour(ax,xgrid/sps,ygrid/sps,zgridgf,ncont,'-','color',[.2 .2 .2]);
%   conmat = contourc(zgridgf,ncont);
  [conmat10th,conobj10th] = contour(ax,xgrid/sps,ygrid/sps,zgridgf,...
    [conpltprc,conpltprc],'-','color','k','linew',1);
elseif strcmp(disttype,'km')
  % [conmat,conobj] = contour(ax,xgrid,ygrid,zgridgf,conplt,'-','color',[.2 .2 .2]);
  [conmat,conobj] = contour(ax,xgrid,ygrid,zgridgf,ncont,'-','color',[.2 .2 .2]);
%   conmat = contourc(xgrid,ygrid,zgridgf,ncont);
  [conmat10th,conobj10th] = contour(ax,xgrid,ygrid,zgridgf,...
    [conpltprc,conpltprc],'-','color','k','linew',1);
end
if ~isempty(conmat)
  contable = getContourLineCoordinates(conmat);
  conmat=table2array(contable);
  contable = getContourLineCoordinates(conmat10th);
  conmat10th=table2array(contable);
%   conmat10th = conmat10th(conmat10th(:,2)==3,:);
%   conmatxy10th = contouroff2space(conmat10th);
%   conmatxy10th = conmatxy10th(conmatxy10th(:,2)==3,:);
%   save('90thprcrangeoflfes.mat','conmat10th','conmatxy10th');
end
  
%plot contour lines
if contourflag
  delete(conobj);
  % delete(conobj10th);
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

% keyboard