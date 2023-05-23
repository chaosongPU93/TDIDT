function f=plt_srcdistNtoNm(dtarvl,eucdist,sps,binwdt,tlensec,binwdist,disttype,timetype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_srcdistNtoNm(dtarvlnn1,distnn1,dtarvlnn2,distnn2,dtarvlnn3,distnn3,dt2all,dist2all,sps)%
%
% Function just to ease the plotting of differential time and distance between
% sources N and N-1 until N and N-4, and each source to all others.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/02/20
% Last modified date:   2023/02/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('disttype','spl');
defval('timetype','tarvl');

widin = 10;  % maximum width allowed is 8.5 inches
htin = 4;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 3;
f = initfig(widin,htin,nrow,ncol); %initialize fig
xran = [0.06 0.96]; yran = [0.12 0.9];
xsep = 0.07; ysep = 0.08;
optaxpos(f,nrow,ncol,xran,yran,xsep,ysep);

h=supertit(f.ax(1:ncol),'Between different source pairs');
movev(h,0.3);

% binwdt = 8/sps;
% if strcmp(disttype,'spl')
%   binwdist = 1;
% elseif strcmp(disttype,'km')
%   binwdist = 0.1;
% end

nsep = size(dtarvl,1);
color = jet(nsep);
symb = ['o';'^';'v';'s';'d';'p';'h';'>';'<';'+';'x';'*';];
ax = f.ax(1);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
for i = 1: nsep
  scatter(ax,dtarvl{i}/sps,eucdist{i},20,color(i,:),symb(i,:));
  plot(ax,ax.XLim,[median(eucdist{i}) median(eucdist{i})],'--','Color',color(i,:));
end
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
xlim(ax,[0 2]);
hold(ax,'off');

ax = f.ax(2);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
for i = 1: nsep
%   hist=histogram(ax,dtarvl{i}/sps,'binwidth',binwdt,'normalization','pdf',...
%     'EdgeColor',color(i,:),'facecolor','none','LineWidth',1);
%   hist.BinEdges = hist.BinEdges-binwdt/2;
  [N,edges]=histcounts(dtarvl{i}/sps,'binwidth',binwdt,'normalization','count');
  edges = edges-binwdt/2;
%   edges = -binwdt/2: binwdt: 2+binwdt/2;
  [N,edges]=histcounts(dtarvl{i}/sps,edges,'normalization','count');
  Nn = N / (tlensec/binwdt);  
  p(i)=stairs(ax,edges,[Nn Nn(end)],'color',color(i,:),'LineWidth',1);
  label{i} = sprintf('N and N-%d',i);
end
legend(ax,p,label);
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
for i = 1: nsep
%   [N,edges]=histcounts(eucdist{i},'binwidth',binwdist,'normalization','pdf');
%   edges = edges-binwdist/2;
%   [N,edges]=histcounts(eucdist{i},edges,'normalization','pdf');
  binedge = (0: binwdist: 50*binwdist)';
  if strcmp(disttype,'spl')
    %if looking at distance in sample space, normalize by integer samples circumvented
    [bincnt,binhgt,count,normalizer] = histbinbypixelinrad(eucdist{i},binedge,'countdensity');
  elseif strcmp(disttype,'km')
    %if looking at distance in map view, normalize by area
    % binedge = [0 binwdist/2: binwdist: 50*binwdist]';
    [bincnt,binhgt,count,normalizer] = histbinbyarea(eucdist{i},binedge,'countdensity');
  end
  stairs(ax,binedge,[binhgt; binhgt(end)],'color',color(i,:),'LineWidth',1);
  plot(ax,[median(eucdist{i}) median(eucdist{i})],ax.YLim,'--','Color',color(i,:));
  text(ax,0.95,0.5-i*0.05,sprintf('med=%.2f',median(eucdist{i})),'Units','normalized',...
    'HorizontalAlignment','right');
end
if strcmp(disttype,'spl')  
  xlabel(ax,'Distance (samples)');
  ylabel(ax,'Count / # of pixels in bin');
elseif strcmp(disttype,'km')
  xlabel(ax,'Distance (km)');
  ylabel(ax,'Count / area (km) of bin');
end
ax.XLim(1) = -binwdist;
hold(ax,'off');



