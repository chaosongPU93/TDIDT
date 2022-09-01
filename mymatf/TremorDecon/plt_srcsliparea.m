function ax = plt_src_slip_area(ax,imp,xran,yran,cran,radi,sps,ftrans,fsty)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax = plt_src_slip_area(ax,imp,xran,yran,cran,radi,sps,ftrans,fsty)
%
% This function is plot the aggregate source region of each source, which is
% represented by a circle centered at the source location with a diameter
% approximated by Vs * Tdura. The area with the most overlapping is assumed 
% to be the region that is slipping during that period of time.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/06/06
% Last modified date:   2022/06/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('loff_max',16);  % allowed searching range for offset12 and 13
defval('maxsybsz',75);  % symbol size of the impulse with the largest amp
defval('ftrans','interpchao');  % default flag of off2space transformation 
defval('fwt','mean');   % as default, use mean amplitude as the weight
defval('ftime','tarvl');   % as default, use arrival time index of impulses at sta 1 as the color scheme

%convert time offset to relative loc
[imploc, indinput] = off2space002(imp(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13 

%plot
plot(ax,[xran(1) xran(2)],[0 0],'--','Color',[.3 .3 .3]); hold(ax,'on');
plot(ax,[0 0],[yran(1) yran(2)],'--','Color',[.3 .3 .3]);

% %outline the boundary of the allowed offset 
% xbndfull = [];
% ybndfull = [];
% for iwin = 1: size(off1iw,1)
% %   iwin=1;
%   off12m = loff_max;
%   off13m = loff_max;
%   offxran = [off1iw(iwin,2)-off12m off1iw(iwin,2)+off12m];
%   offyran = [off1iw(iwin,3)-off13m off1iw(iwin,3)+off13m];
%   xbndlow = [offxran(1):offxran(2), offxran(2)*ones(1,off13m)];
%   xbndupp = [offxran(2)-1:-1:offxran(1), offxran(1)*ones(1,off13m)];
%   ybndlow = [offyran(1)*ones(1, off12m), offyran(1):offyran(2)];
%   ybndupp = [offyran(2)*ones(1, off12m-1), offyran(2):-1:offyran(1)];
%   xbnd = [xbndlow xbndupp]';
%   ybnd = [ybndlow ybndupp]';
%   %convert time offset to relative loc
%   [bndloc, ~] = off2space002([xbnd,ybnd],sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
%   plot(ax,bndloc(:,1),bndloc(:,2),'-','Color',[.7 .7 .7],'linew',2);
%   xbndfull = [xbndfull; xbnd];
%   ybndfull = [ybndfull; ybnd];
% end
% indbnd = boundary(xbndfull,ybndfull,1);
% %convert time offset to relative loc
% [bndloc, ~] = off2space002([xbndfull(indbnd),ybndfull(indbnd)],sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
% plot(ax,bndloc(:,1),bndloc(:,2),'-','Color',[.4 .4 .4],'linew',2);

if isequal(fsty,'shade')
  %%%option 1 to define 'relative origin time', time 0 is the origin time of a source at 0,0 and
  %%%arrive at the start of your window at PGC
  [imploc0, ~] = off2space002([0 0],sps,ftrans,0);  % a ref source at 0,0
  tcor = imploc(:,6)-imploc0(6);
  imptime = imp(indinput,1)/sps-tcor;
%   cran = [0 max(imptime)];
  color = colormatch(imptime,cran,'kelicol');

  for i = 1: size(imploc,1)
    x0 = imploc(i,1);
    y0 = imploc(i,2);
    [x, y] = circle_chao(x0,y0,radi,0.1);
%     patch(x,y,'k','Facealpha',0.1,'edgecolor','none');
    patch(x,y,color(i,:),'Facealpha',0.15,'edgecolor','none'); 
    
  end
  
elseif isequal(fsty,'line')
  for i = 1: size(imploc,1)
    x0 = imploc(i,1);
    y0 = imploc(i,2);
    [x, y] = circle_chao(x0,y0,radi,0.1);
    scatter(ax,x0,y0,10,'ko','filled');
    plot(ax,x,y,'-','Color',[.5 .5 .5],'linew',1);
  end
end
oldc = colormap(ax,'kelicol');
newc = flipud(oldc);
colormap(ax,newc);
c=colorbar(ax);
c.Label.String = 'Relative origin time (s)';
c.Label.FontSize = 9;
caxis(ax,cran);
ax.Box = 'on'; grid(ax, 'on');
text(ax,0.98,0.95,num2str(size(imploc,1)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',10);
ax.YAxis.FontSize = 8;
ax.XAxis.FontSize = 8;
xlabel(ax,'E (km)','FontSize',11);
ylabel(ax,'N (km)','FontSize',11);
axis(ax, 'equal');
xlim(ax,xran); xticks(ax,xran(1): 0.5 : xran(2));
ylim(ax,yran); yticks(ax,yran(1): 0.5 : yran(2));

% keyboard
