function [ax,imptime,wt,xbndcvhl,ybndcvhl,c] = plt_decon_imp_scatter_ref(ax,imp,xran,yran,cran,off1iw,loff_max,...
  sps,maxsybsz,fwt,ftime,fbnd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ax,imptime,wt,xbndcvhl,ybndcvhl,c] = plt_decon_imp_scatter_ref(ax,imp,xran,yran,cran,off1iw,loff_max,...
%     sps,maxsybsz,fwt,ftime,fbnd)
%
% Similar to 'plt_decon_imp_scatter.m', This function is to plot the 
% scatter of the offsets in time between
% the deconvolved impulses at trio stations, off12 and off13, 
% colorcoded by the index of impulse at station 1. 
% The difference is, this is designed for 'refinement' where triplets
% were grouped using a variable allowed searching range relative the 
% best alignement upon each subwin of the entire data win. Therefore,
% part of the code that plots the combined boundary needs to be changed
% mainly.
%
% Sympol size is 
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
%   boundary without the individual ones, or, vice versa, or both.
%   Meanwhile, provide another map of individual boundaries shifted
%   to the same origin 0,0 in 'plt_decon_imp_scatter_ref_sft.m'
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/05/25
% Last modified date:   2022/05/25 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('maxsybsz',75);  % symbol size of the impulse with the largest amp
defval('fwt','mean');   % as default, use 'mean' amplitude as the weight
defval('ftime','tarvl');   % as default, use arrival time index of impulses at sta 1 as the color scheme
defval('fbnd','comb');   % as default, plot the combined boundary from all subwins

% f.fig = figure;
% f.fig.Renderer = 'painters';
plot(ax, [xran(1) xran(2)],[0 0],'--','Color',[.3 .3 .3]); hold(ax,'on');
plot(ax, [0 0],[yran(1) yran(2)],'--','Color',[.3 .3 .3]); 

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
    plot(ax, xbnd,ybnd,'-','Color',[.7 .7 .7],'linew',2);
  end
  xbndfull = [xbndfull; xbnd];
  ybndfull = [ybndfull; ybnd];  
end
% indbnd = boundary(xbndfull,ybndfull,1);
indbnd = convhull(xbndfull,ybndfull); % convex hull, recommended this way
xbndcvhl = xbndfull(indbnd);
ybndcvhl = ybndfull(indbnd);
if isequal(fbnd,'comb') || isequal(fbnd,'both')
  plot(ax, xbndcvhl,ybndcvhl,'-','Color',[.4 .4 .4],'linew',2);
end

if isequal(fwt,'median')
  wt = median(imp(:,[2 4 6]),2);
elseif isequal(fwt,'mean')
  wt = mean(imp(:,[2 4 6]),2);
end
%make the impulse with largest weighted amplitude to have the size of 'maxsybsz', then scale the 
%other relatively
% wtmax = max(wt);
wtmax = prctile(wt,95); %use percentile in case 
refscl = wt./wtmax; 
refscl(refscl>=1) = 1;  %force the larger amp to be plotted as the same size in case of saturation

%how to colorcode the impulses? use arrival time or origin time?
if isequal(ftime,'tarvl')
  imptime = imp(:,1);
elseif isequal(ftime,'tori')
  %convert time offset to relative loc, but we mainly need the travel time estimate from sources
  ftrans = 'interpchao';
  [imploc, ~] = off2space002(imp(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
  %%%option 1 to define 'relative origin time', time 0 is the origin time of a source at 0,0 and
  %%%arrive at the start of your window at PGC
  [imploc0, ~] = off2space002([0 0],sps,ftrans,0);  % a ref source at 0,0
  tcor = round((imploc(:,6)-imploc0(6))*sps);   % travel time difference between each source and a ref source at 0,0
  imptime = imp(:,1)-tcor;
%   %%%option 2 to define 'relative origin time', time 0 is the true origin time of the first source 
%   %%%of your window
%   imptime = imp(:,1)-round(imploc(:,6)*sps);
%   imptime = imptime - min(imptime);
%   cran = [0 ceil(max(imptime))];
end

[imptimest, indsort] = sortrows(imptime,1);  % sort based on the time (no matter origin or arrival)
scatter(ax,imp(indsort,7),imp(indsort,8),maxsybsz*refscl,imptimest,'filled',...
  'MarkerEdgeColor',[.5 .5 .5]);
% scatter(ax, imp(:,7),imp(:,8),maxsybsz*refscl,imptime,'filled','MarkerEdgeColor',[.5 .5 .5]);
% scatter(ax, imp(:,7),imp(:,8),maxsybsz*refscl,imptime,'filled');
oldc = colormap(ax, 'kelicol');
newc = flipud(oldc);
colormap(ax, newc);
c=colorbar(ax);
caxis(ax, cran);
if isequal(ftime,'tarvl')
  c.Label.String = sprintf('Index of impulse at PGC at %d Hz', sps);
elseif isequal(ftime,'tori')
%   c.Label.String = sprintf('Traveltime corrected index of impulse at PGC at %d Hz', sps);
%   c.Label.String = sprintf('Relative origin time (samples) of the impulse at %d Hz', sps);
  c.Label.String = sprintf('Relative origin time (samples)');
end
c.Label.FontSize = 9;
% pos = ax.Position;
% c.Position = [pos(1)+pos(3)*16/15, pos(2), pos(3)*1/20, pos(4)];
ax.Box = 'on'; grid(ax, 'on');
scatter(ax,xran(1)+0.05*range(xran),yran(2)-0.05*range(yran),maxsybsz,'w','filled',...
  'MarkerEdgeColor',[.5 .5 .5]);
% scatter(ax,xran(1)+0.05*range(xran),yran(2)-0.05*range(yran),maxsybsz,[.5 .5 .5],'filled');
text(ax,0.02,0.9,strcat({'Amp. \geq '},sprintf('%.1f',wtmax)),'Units','normalized',...
  'HorizontalAlignment','left','FontSize',8);
% text(ax,0.98,0.15,sprintf('med. of abs.: %.2f, %.2f, %.2f',median(abs(imp(:,7))),...
%   median(abs(imp(:,8))),median(abs(imp(:,9)))),'Units','normalized',...
%   'HorizontalAlignment','right');
% text(ax,0.98,0.1,sprintf('wgt. med. of abs.: %.2f, %.2f, %.2f',wt_median(abs(imp(:,7)),wt),...
%   wt_median(abs(imp(:,8)),wt),wt_median(abs(imp(:,9)),wt)),'Units','normalized',...
%   'HorizontalAlignment','right');
% text(ax,0.98,0.05,sprintf('wgt. mean of abs.: %.2f, %.2f, %.2f',wt_mean(abs(imp(:,7)),wt),...
%   wt_mean(abs(imp(:,8)),wt),wt_mean(abs(imp(:,9)),wt)),'Units','normalized',...
%   'HorizontalAlignment','right');
text(ax,0.98,0.95,num2str(size(imp,1)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',10);
ax.YAxis.FontSize = 8;
ax.XAxis.FontSize = 8;
xlabel(ax,sprintf('PGC-SSIB offset (samples at %d Hz)',sps),'FontSize',10);
ylabel(ax,sprintf('PGC-SILB offset (samples at %d Hz)',sps),'FontSize',10);
axis(ax, 'equal');
xlim(ax,xran); xticks(ax,xran(1): max(0.5,2) : xran(2));
ylim(ax,yran); yticks(ax,yran(1): max(0.5,2) : yran(2));
% pos = ax.Position;
% c.Position = [pos(1)+pos(3)+0.02, pos(2), 0.03, pos(4)];
% c.Position(3) = 0.03;
% c.Position(1) = c.Position(1) - 0.01;

% keyboard

