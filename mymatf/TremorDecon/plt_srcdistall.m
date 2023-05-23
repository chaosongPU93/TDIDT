function f=plt_srcdistall(dt2all,dist2all,sps,binwdt,tlensec,binwdist,disttype,timetype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_srcdistall(dt2all,dist2all,sps,binwdt,binwdist,disttype)
%
% Function just to ease the plotting of euclidean distance between each LFE 
% source to all other LFEs.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/03/08
% Last modified date:   2023/03/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('disttype','spl');
defval('timetype','tarvl');

widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
nrow = 2;
ncol = 2;
f = initfig(widin,htin,nrow,ncol); %initialize fig
xran = [0.10 0.95]; yran = [0.10 0.95];
xsep = 0.09; ysep = 0.08;
optaxpos(f,nrow,ncol,xran,yran,xsep,ysep);

h=supertit(f.ax(1:ncol),'Between each source and all others');
movev(h,0.1);

ax = f.ax(1);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
scatter(ax,dt2all/sps,dist2all,15,'ko');
if strcmp(disttype,'spl')
  ylabel(ax,'Distance (samples)');
elseif strcmp(disttype,'km')
  ylabel(ax,'Distance (km)');
end
if strcmp(timetype,'tarvl')
  xlabel(ax,'Differential arrival time (s)');
elseif strcmp(timetype,'tori')
  xlabel(ax,'Differential origin time (s)');
end
hold(ax,'off');

ax = f.ax(2);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% hist=histogram(ax,dt2all/sps,'binwidth',binwdt,'normalization','count','facecolor','k');
% hist.BinEdges = hist.BinEdges-binwdt/2;
% [N,edges]=histcounts(dt2all/sps,'binwidth',binwdt,'normalization','count');
% edges = edges-binwdt/2;
edges = -binwdt/2: binwdt: max(dt2all/sps)+binwdt/2;
[N,edges]=histcounts(dt2all/sps,edges,'normalization','count');
Nn = N / (tlensec/binwdt);  
% bar(ax,cnt(1:end-1),Nn,1,'stacked','k','facea',0.6);
stairs(ax,edges,[Nn Nn(end)],'k','LineWidth',1);
% keyboard
plot(ax,[median(dt2all/sps) median(dt2all/sps)],ax.YLim,'r--','linew',2);
text(ax,0.95,0.95,sprintf('med=%.2f',median(dt2all/sps)),'Units','normalized',...
  'HorizontalAlignment','right');
ylabel(ax,'Norm. count');
if strcmp(timetype,'tarvl')
  xlabel(ax,'Differential arrival time (s)');
elseif strcmp(timetype,'tori')
  xlabel(ax,'Differential origin time (s)');
end
xlim(ax,[-binwdt f.ax(1).XLim(2)+binwdt]);
hold(ax,'off');

ax = f.ax(3);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
% hist=histogram(ax,dist2all,'binwidth',binwdist,'normalization','pdf','facecolor','k');
% hist.BinEdges = hist.BinEdges-binwdist/2;
binedge = (0: binwdist: 50*binwdist)';
if strcmp(disttype,'spl')
  %if looking at distance in sample space, normalize by integer samples circumvented
  [bincnt,binhgt,count,normalizer] = histbinbypixelinrad(dist2all,binedge,'countdensity');
elseif strcmp(disttype,'km')
  %if looking at distance in map view, normalize by area
  % binedge = [0 binwdist/2: binwdist: 50*binwdist]';
  [bincnt,binhgt,count,normalizer] = histbinbyarea(dist2all,binedge,'countdensity');
end
% bar(ax,bincnt,binhgt,1,'stacked','k','facea',0.6);
stairs(ax,binedge,[binhgt; binhgt(end)],'k','LineWidth',1);
plot(ax,[median(dist2all) median(dist2all)],ax.YLim,'r--','linew',2);
text(ax,0.95,0.95,sprintf('med=%.2f',median(dist2all)),'Units','normalized',...
  'HorizontalAlignment','right');
if strcmp(disttype,'spl')
  xlabel(ax,'Distance (samples)');
  ylabel(ax,'Count / # of pixels in bin');
elseif strcmp(disttype,'km')
  xlabel(ax,'Distance (km)');
  ylabel(ax,'Count / area (km) of bin');
end
ax.XLim(1) = -binwdist;
hold(ax,'off');

ax = f.ax(4);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
nbin = length(edges)-1;
color = jet(nbin);
binedge = (0: binwdist: 50*binwdist)';
for i = 1: nbin
  distplt = dist2all((dt2all/sps)>=edges(i) & (dt2all/sps)<edges(i+1));
  if strcmp(disttype,'spl')
    %if looking at distance in sample space, normalize by integer samples circumvented
    [bincnt,binhgt,count,normalizer] = histbinbypixelinrad(distplt,binedge,'countdensity');
  elseif strcmp(disttype,'km')
    %if looking at distance in map view, normalize by area
    % binedge = [0 binwdist/2: binwdist: 50*binwdist]';
    [bincnt,binhgt,count,normalizer] = histbinbyarea(distplt,binedge,'countdensity');
  end
  p(i)=stairs(ax,binedge,[binhgt; binhgt(end)],'color',color(i,:),'LineWidth',1);
  label{i} = sprintf('[%.3f, %.3f] s',edges(i),edges(i+1));
end
legend(ax,p,label);
if strcmp(disttype,'spl')
  xlabel(ax,'Distance (samples)');
  ylabel(ax,'Count / # of pixels in bin');
elseif strcmp(disttype,'km')
  xlabel(ax,'Distance (km)');
  ylabel(ax,'Count / area (km) of bin');
end
ax.XLim(1) = -binwdist;
hold(ax,'off');


% keyboard

