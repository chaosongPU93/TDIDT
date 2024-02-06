function [f,mdist2all]=plt_srcdlocall(dloc2all,binw,disttype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,mdist2all]=plt_srcdlocall(dloc2all,binw,disttype)
%
% Function just to ease the plotting of differential location between
% between each LFE source to all others
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/03/08
% Last modified date:   2023/03/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('binw',1);
defval('disttype','spl');

widin = 6;  % maximum width allowed is 8.5 inches
if strcmp(disttype,'spl')
  nrow = 3;
  htin = 9;   % maximum height allowed is 11 inches
elseif strcmp(disttype,'km')
  nrow = 2;
  htin = 6.5;   % maximum height allowed is 11 inches
end
ncol = 2;
f = initfig(widin,htin,nrow,ncol); %initialize fig

xran = [0.1 0.96]; yran = [0.1 0.96];
xsep = 0.09; ysep = 0.1;
optaxpos(f,nrow,ncol,xran,yran,xsep,ysep);

if strcmp(disttype,'spl')
  X = -50:binw:50;
elseif strcmp(disttype,'km')
  X = -5:binw:5;
end
% supertit(f.ax(1:ncol),'Between each source and all others');

ax=f.ax(1);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
hist=histogram(ax,dloc2all(:,1),'binwidth',binw,'normalization','countdensity','facecolor','k');
% hist.BinEdges = hist.BinEdges-binw/2;
[MUHAT,SIGMAHAT] = normfit(dloc2all(:,1));
Y=normpdf(X,MUHAT,SIGMAHAT)*size(dloc2all,1);
plot(ax,X,Y,'Color','r','linew',2);
text(ax,0.02,0.95,sprintf('\\mu=%.2f, \\sigma=%.2f',MUHAT,SIGMAHAT),'Units','normalized',...
  'HorizontalAlignment','left');
if strcmp(disttype,'spl')
  xlabel(ax,'Diff off12 (samples)');
  ylabel(ax,'Count per sample');
elseif strcmp(disttype,'km')
  xlabel(ax,'Diff E loc (km)');
  ylabel(ax,'Count per km');
end
axsym(ax,1);
ax.YAxis.Exponent = 4;
longticks(ax,2);
hold(ax,'off');

ax=f.ax(2);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
[N,edges]=histcounts(abs(dloc2all(:,1)),'binwidth',binw,'normalization','countdensity');
% edges = edges-binw/2;
% [N,edges]=histcounts(abs(dloc2all(:,1)),edges,'normalization','countdensity');
Nn = N;
if strcmp(disttype,'spl')
  Nn(1) = N(1)*2;
  Nn = Nn*(sum(N)/sum(Nn));
end
% bar(ax,cnt(1:end-1),Nn,1,'stacked','k','facea',0.6);
stairs(ax,edges,[Nn Nn(end)],'k','LineWidth',1);
plot(ax,[median(abs(dloc2all(:,1))) median(abs(dloc2all(:,1)))],ax.YLim,'r--','linew',2);
text(ax,0.98,0.95,sprintf('med=%.2f km',median(abs(dloc2all(:,1)))),'Units','normalized',...
  'HorizontalAlignment','right');
if strcmp(disttype,'spl')
  xlabel(ax,'Abs diff off12 (samples)');
  ylabel(ax,'Count per sample');
elseif strcmp(disttype,'km')
  xlabel(ax,'Abs diff E loc (km)');
  ylabel(ax,'Count per km');
end
ax.YAxis.Exponent = 4;
longticks(ax,2);
hold(ax,'off');

ax=f.ax(3);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
hist=histogram(ax,dloc2all(:,2),'binwidth',binw,'normalization','countdensity','facecolor','k');
% hist.BinEdges = hist.BinEdges-binw/2;
[MUHAT,SIGMAHAT] = normfit(dloc2all(:,2));
Y=normpdf(X,MUHAT,SIGMAHAT)*size(dloc2all,1);
plot(ax,X,Y,'Color','r','linew',2);
text(ax,0.02,0.95,sprintf('\\mu=%.2f, \\sigma=%.2f',MUHAT,SIGMAHAT),'Units','normalized',...
  'HorizontalAlignment','left');
if strcmp(disttype,'spl')
  xlabel(ax,'Diff off13 (samples)');
  ylabel(ax,'Count per sample');
elseif strcmp(disttype,'km')
  xlabel(ax,'Diff N loc (km)');
  ylabel(ax,'Count per km');
end
axsym(ax,1);
ax.YAxis.Exponent = 3;
longticks(ax,2);
hold(ax,'off');

ax=f.ax(4);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
[N,edges]=histcounts(abs(dloc2all(:,2)),'binwidth',binw,'normalization','countdensity');
% edges = edges-binw/2;
% [N,edges]=histcounts(abs(dloc2all(:,2)),edges,'normalization','countdensity');
Nn = N;
if strcmp(disttype,'spl')
  Nn(1) = N(1)*2;
  Nn = Nn*(sum(N)/sum(Nn));
end
% bar(ax,cnt(1:end-1),Nn,1,'stacked','k','facea',0.6);
stairs(ax,edges,[Nn Nn(end)],'k','LineWidth',1);
plot(ax,[median(abs(dloc2all(:,2))) median(abs(dloc2all(:,2)))],ax.YLim,'r--','linew',2);
text(ax,0.98,0.95,sprintf('med=%.2f km',median(abs(dloc2all(:,2)))),'Units','normalized',...
  'HorizontalAlignment','right');
if strcmp(disttype,'spl')
  xlabel(ax,'Abs diff off13 (samples)');
  ylabel(ax,'Count per sample');
elseif strcmp(disttype,'km')
  xlabel(ax,'Abs diff N loc (km)');
  ylabel(ax,'Count per km');
end
ax.YAxis.Exponent = 4;
longticks(ax,2);
hold(ax,'off');

if strcmp(disttype,'spl')
  ax=f.ax(5);
  hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');  
  hist=histogram(ax,dloc2all(:,2)-dloc2all(:,1),'binwidth',binw,'normalization',...
    'countdensity','facecolor','k');
%   hist.BinEdges = hist.BinEdges-binw/2;
  [MUHAT,SIGMAHAT] = normfit(dloc2all(:,2)-dloc2all(:,1));
  Y=normpdf(X,MUHAT,SIGMAHAT)*size(dloc2all,1);
  plot(ax,X,Y,'Color','r','linew',2);
  text(ax,0.02,0.95,sprintf('\\mu=%.2f, \\sigma=%.2f',MUHAT,SIGMAHAT),'Units','normalized',...
    'HorizontalAlignment','left');
  xlabel(ax,'Diff off23 (samples)');
  ylabel(ax,'Count per sample');
  axsym(ax,1);
  hold(ax,'off');

  ax=f.ax(6);
  hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
  [N,edges]=histcounts(abs(dloc2all(:,2)-dloc2all(:,1)),'binwidth',binw,'normalization','countdensity');
  % edges = edges-binw/2;
  % [N,edges]=histcounts(abs(dloc2all(:,2)-dloc2all(:,1)),edges,'normalization','countdensity');
  Nn = N;
  if strcmp(disttype,'spl')
    Nn(1) = N(1)*2;
    Nn = Nn*(sum(N)/sum(Nn));
  end
  % bar(ax,cnt(1:end-1),Nn,1,'stacked','k','facea',0.6);
  stairs(ax,edges,[Nn Nn(end)],'k','LineWidth',1);
  plot(ax,[median(abs(dloc2all(:,2)-dloc2all(:,1))) median(abs(dloc2all(:,2)-dloc2all(:,1)))],...
    ax.YLim,'r--','linew',2);
  text(ax,0.98,0.95,sprintf('med=%.2f',median(abs(dloc2all(:,2)-dloc2all(:,1)))),'Units','normalized',...
    'HorizontalAlignment','right');
  xlabel(ax,'Abs diff off23 (samples)');
  ylabel(ax,'Count per sample');
  hold(ax,'off');
end

mdist2all = [median(abs(dloc2all(:,1))) median(abs(dloc2all(:,2)))];




