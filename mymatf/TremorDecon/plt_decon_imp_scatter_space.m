function [ax,imptime,wt,bndloc] = plt_decon_imp_scatter_space(ax,imp,xran,yran,cran,offxran,offyran,...
  sps,maxsybsz,ftrans,fwt,ftime,fbnd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ax,imptime,wt] = plt_decon_imp_scatter_space(ax,imp,xran,yran,cran,offxran,offyran,...
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
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/02/12
% Last modified date:   2022/05/03 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('offxran',[-16 16]);  % offset12  boundary range
defval('offyran',[-16 16]);  % offset13  boundary range
defval('maxsybsz',75);  % symbol size of the impulse with the largest amp
defval('ftrans','interpchao');  % default flag of off2space transformation 
defval('fwt','mean');   % as default, use mean amplitude as the weight
defval('ftime','tarvl');   % as default, use arrival time index of impulses at sta 1 as the color scheme
defval('fbnd','each');   % as default, plot the search boundary

%convert time offset to relative loc
[imploc, indinput] = off2space002(imp(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%plot
% f.fig = figure;
% f.fig.Renderer = 'painters';
plot(ax,[xran(1) xran(2)],[0 0],'--','Color',[.3 .3 .3]); hold(ax,'on');
plot(ax,[0 0],[yran(1) yran(2)],'--','Color',[.3 .3 .3]);
%outline the boundary of the allowed offset 
off12m = round(range(offxran)/2);
off13m = round(range(offyran)/2);
xbndlow = [offxran(1):offxran(2), offxran(2)*ones(1,off13m)];
xbndupp = [offxran(2)-1:-1:offxran(1), offxran(1)*ones(1,off13m)];
ybndlow = [offyran(1)*ones(1, off12m), offyran(1):offyran(2)];
ybndupp = [offyran(2)*ones(1, off12m-1), offyran(2):-1:offyran(1)];
xbnd = [xbndlow xbndupp]';
ybnd = [ybndlow ybndupp]';
%convert time offset to relative loc
[bndloc, ~] = off2space002([xbnd,ybnd],sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
if isequal(fbnd,'each')
  plot(ax,bndloc(:,1),bndloc(:,2),'-','Color',[.4 .4 .4],'linew',2);
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
  c.Label.String = 'Arrival time at PGC (s)';
elseif isequal(ftime,'tori')
%   c.Label.String = sprintf('Traveltime corrected index of impulse at PGC at %d Hz', sps);
%   c.Label.String = 'Traveltime corrected time (s) of impulse at PGC since the start';
%   c.Label.String = 'Relative origin time (s) of the impulse';
  c.Label.String = 'Relative origin time (s)';
end
c.Label.FontSize = 9;
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
xlabel(ax,'E (km)','FontSize',11);
ylabel(ax,'N (km)','FontSize',11);
axis(ax, 'equal');
xlim(ax,xran); xticks(ax,xran(1): 1 : xran(2));
ylim(ax,yran); yticks(ax,yran(1): 1 : yran(2));
% pos = ax.Position;
% c.Position = [pos(1)+pos(3)+0.02, pos(2), 0.03, pos(4)];

% polyarea(ax,imploc(:,1),imploc(:,2))

% keyboard