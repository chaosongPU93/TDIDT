function [ax,imptime,wt,bndloccvhl,c] = plt_decon_imp_scatter_space_ref(ax,imp,xran,yran,cran,off1iw,loff_max,...
  sps,maxsybsz,ftrans,fwt,ftime,fbnd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ax,imptime,wt,bndloccvhl,c] = plt_decon_imp_scatter_space(ax,imp,xran,yran,cran,offxran,offyran,...
%   sps,maxsybsz,ftrans,fwt,ftime,fbnd)
%
% This function is to plot the scatter of the offsets in time between
% the deconvolved impulses at trio stations in terms of relative spatial 
% locations x and y which are transformed from the off12 and off13, 
% colorcoded by the index of impulse at station 1. Sympol size is 
% proportionally scaled by the amplitude of impulse at station 1. The 
% largest ampltiude impulse will be plotted with a size of 'maxsybsz', 
% where the amplitude is either the 'median' or 'mean' of the same 
% group of the impulses. 
% The colorcoding now has another option other than using the 'tarvl' at
% station 1, ie., the arrival index of impulse at sta 1, you can also use
% the option 'tori' to use the relative origin time index corrected by 
% the actual variation in source locations. 
% 
% --2022/06/06, adding the option to plot either the convex composite 
%   boundary without the individual ones, or, vice versa.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/02/12
% Last modified date:   2022/05/03 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('loff_max',16);  % allowed searching range for offset12 and 13
defval('maxsybsz',75);  % symbol size of the impulse with the largest amp
defval('ftrans','interpchao');  % default flag of off2space transformation 
defval('fwt','mean');   % as default, use mean amplitude as the weight
defval('ftime','tarvl');   % as default, use arrival time index of impulses at sta 1 as the color scheme
defval('fbnd','comb');   % as default, plot the combined boundary from all subwins

%convert time offset to relative loc
[imploc, indinput] = off2space002(imp(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%plot
% f.fig = figure;
% f.fig.Renderer = 'painters';
plot(ax,[xran(1) xran(2)],[0 0],'--','Color',[.3 .3 .3]); hold(ax,'on');
plot(ax,[0 0],[yran(1) yran(2)],'--','Color',[.3 .3 .3]);

%outline the boundary of the allowed offset 
xbndfull = [];
ybndfull = [];
for iwin = 1: size(off1iw,1)
%   iwin=1;
  off12m = loff_max;
  off13m = loff_max;
  offxran = [off1iw(iwin,2)-off12m off1iw(iwin,2)+off12m];
  offyran = [off1iw(iwin,3)-off13m off1iw(iwin,3)+off13m];
  xbndlow = [offxran(1):offxran(2), offxran(2)*ones(1,off13m)];
  xbndupp = [offxran(2)-1:-1:offxran(1), offxran(1)*ones(1,off13m)];
  ybndlow = [offyran(1)*ones(1, off12m), offyran(1):offyran(2)];
  ybndupp = [offyran(2)*ones(1, off12m-1), offyran(2):-1:offyran(1)];
  xbnd = [xbndlow xbndupp]';
  ybnd = [ybndlow ybndupp]';
  if isequal(fbnd,'each') || isequal(fbnd,'both')
    %convert time offset to relative loc
    [bndloc, ~] = off2space002([xbnd,ybnd],sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
    plot(ax,bndloc(:,1),bndloc(:,2),'-','Color',[.7 .7 .7],'linew',2);
  end
  xbndfull = [xbndfull; xbnd];
  ybndfull = [ybndfull; ybnd];
end
% indbnd = boundary(xbndfull,ybndfull,1);
indbnd = convhull(xbndfull,ybndfull); % convex hull, recommended this way
if isequal(fbnd,'comb') || isequal(fbnd,'both')
  %convert time offset to relative loc
  [bndloccvhl, ~] = off2space002([xbndfull(indbnd),ybndfull(indbnd)],sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
  plot(ax,bndloccvhl(:,1),bndloccvhl(:,2),'-','Color',[.4 .4 .4],'linew',2);
else
  bndloccvhl = [];
end

if isequal(fwt,'median')
  wt = median(imp(indinput,[2 4 6]),2);
elseif isequal(fwt,'mean')
  wt = mean(imp(indinput,[2 4 6]),2);
end
%make the impulse with largest weighted amplitude to have the size of 'maxsybsz', then scale the 
%other relatively
% wtmax = max(wt);
wtmax = prctile(wt,95); %use percentile in case 
refscl = wt./wtmax; 
refscl(refscl>=1) = 1;  %force the larger amp to be plotted as the same size in case of saturation

%how to colorcode the impulses? use arrival time or origin time?
if isequal(ftime,'tarvl')
  imptime = imp(indinput,1)/sps;
elseif isequal(ftime,'tori')
  %%%option 1 to define 'relative origin time', time 0 is the origin time of a source at 0,0 and
  %%%arrive at the start of your window at PGC
  [imploc0, ~] = off2space002([0 0],sps,ftrans,0);  % a ref source at 0,0
  tcor = imploc(:,6)-imploc0(6);
  imptime = imp(indinput,1)/sps-tcor;
%   %%%option 2 to define 'relative origin time', time 0 is the true origin time of the first source 
%   %%%of your window
%   imptime = imp(indinput,1)/sps-imploc(:,6);
%   imptime = imptime - min(imptime);
%   cran = [0 ceil(max(imptime))];
end

[imptimest, indsort] = sortrows(imptime,1);  % sort based on the time (no matter origin or arrival)
scatter(ax,imploc(indsort,1),imploc(indsort,2),maxsybsz*refscl,imptimest,'filled',...
  'MarkerEdgeColor',[.5 .5 .5]);
% scatter(ax,imploc(:,1),imploc(:,2),maxsybsz*refscl,imptime,'filled','MarkerEdgeColor',[.5 .5 .5]);
% scatter(ax,imploc(:,1),imploc(:,2),maxsybsz*refscl,imptime,'filled');
oldc = colormap(ax,'kelicol');
newc = flipud(oldc);
colormap(ax,newc);
c=colorbar(ax);
if isequal(ftime,'tarvl')
%   c.Label.String = sprintf('Index of impulse at PGC at %d Hz', sps);
  c.Label.String = 'Relative arrival time (s) at PGC';
elseif isequal(ftime,'tori')
%   c.Label.String = sprintf('Traveltime corrected index of impulse at PGC at %d Hz', sps);
%   c.Label.String = 'Traveltime corrected time (s) of impulse at PGC since the start';
%   c.Label.String = 'Relative origin time (s) of the impulse';
  c.Label.String = 'Relative origin time (s)';
end
c.Label.FontSize = 9;
% pos = ax.Position;
% c.Position = [pos(1)+pos(3)+0.01, pos(2), 0.02, pos(4)];
caxis(ax,cran);
ax.Box = 'on'; grid(ax, 'on');
% scatter(ax,xran(1)+0.05*range(xran),yran(2)-0.05*range(yran),maxsybsz,'w','filled',...
%   'MarkerEdgeColor',[.5 .5 .5]);
% text(ax,0.02,0.9,strcat({'Amp. \geq '},sprintf('%.1f',wtmax)),'Units','normalized',...
%   'HorizontalAlignment','left','FontSize',8);
% text(ax,0.98,0.15,sprintf('med. of abs.: %.2f km, %.2f km',median(abs(imploc(:,1))),...
%   median(abs(imploc(:,2)))),'Units','normalized','HorizontalAlignment','right');
% text(ax,0.98,0.1,sprintf('wgt. med. of abs.: %.2f km, %.2f km',...
%   wt_median(abs(imploc(:,1)),wt),...
%   wt_median(abs(imploc(:,2)),wt)),'Units','normalized',...
%   'HorizontalAlignment','right');
% text(ax,0.98,0.05,sprintf('wgt. mean of abs.: %.2f km, %.2f km',...
%   wt_mean(abs(imploc(:,1)),wt),...
%   wt_mean(abs(imploc(:,2)),wt)),'Units','normalized',...
%   'HorizontalAlignment','right');
text(ax,0.98,0.05,sprintf('%d events',size(imploc,1)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
ax.YAxis.FontSize = 8;
ax.XAxis.FontSize = 8;
xlabel(ax,'E (km)','FontSize',10);
ylabel(ax,'N (km)','FontSize',10);
axis(ax, 'equal');
xlim(ax,xran); xticks(ax,xran(1): 1 : xran(2));
ylim(ax,yran); yticks(ax,yran(1): 1 : yran(2));
% pos = ax.Position;
% c.Position = [pos(1)+pos(3)+0.02, pos(2), 0.03, pos(4)];

% polyarea(ax,imploc(:,1),imploc(:,2))

% keyboard



