function f=plt_srcdlocNtoNm(dlocxy,binw,disttype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_srcdlocNtoNm(dlocxy,disttype)
%
% Function just to ease the plotting of differential location between
% sources N and N-1 until N and N-m, where m is determined by input data size.
% Plot the direct difference.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/02/20
% Last modified date:   2023/02/20
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

m = size(dlocxy,1);
color = jet(m);

if strcmp(disttype,'spl')
  X = -50:binw:50;
elseif strcmp(disttype,'km')
  X = -5:binw:5;
end
% supertit(f.ax(1:ncol),'Between different source pairs');

ax=f.ax(1);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
for i = m:-1:1
  aa=dlocxy{i};
  [N,edges]=histcounts(aa(:,1),'binwidth',binw,'normalization','countdensity');
%   edges = edges-binw/2;
%   [N,edges]=histcounts(aa(:,1),edges,'normalization','countdensity');
  stairs(ax,edges,[N N(end)],'color',color(i,:),'LineWidth',1);
  [MUHAT,SIGMAHAT] = normfit(aa(:,1));
  Y=normpdf(X,MUHAT,SIGMAHAT)*size(aa,1);
  plot(ax,X,Y,'Color',color(i,:),'linew',1);
  text(ax,0.02,0.5-i*0.05,sprintf('\\mu=%.2f, \\sigma=%.2f',MUHAT,SIGMAHAT),'Units','normalized',...
    'HorizontalAlignment','left'); %,'Color',color(i,:)
end
if strcmp(disttype,'spl')
  xlabel(ax,'Diff off12 (samples)');
  ylabel(ax,'Count per sample');
elseif strcmp(disttype,'km')
  xlabel(ax,'Diff E loc (km)');
  ylabel(ax,'Count per km');
end
axsym(ax,1);
% xlim(ax,[-50 50]);
ax.YAxis.Exponent = 3;
longticks(ax,2);
hold(ax,'off');

ax=f.ax(2);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
for i = m:-1:1
  aa=dlocxy{i};
  [N,edges]=histcounts(abs(aa(:,1)),'binwidth',binw,'normalization','countdensity');
%   edges = edges-binw/2;
%   [N,edges]=histcounts(abs(aa(:,1)),edges,'normalization','countdensity');
  Nn = N;
  if strcmp(disttype,'spl')
    Nn(1) = N(1)*2;
    Nn = Nn*(sum(N)/sum(Nn));
  end
  p(i)=stairs(ax,edges,[Nn Nn(end)],'color',color(i,:),'LineWidth',1);
  plot(ax,[median(abs(aa(:,1))) median(abs(aa(:,1)))],ax.YLim,'--','Color',color(i,:));
  text(ax,0.98,0.5-i*0.05,sprintf('med=%.2f km',median(abs(aa(:,1)))),'Units','normalized',...
    'HorizontalAlignment','right'); %,'Color',color(i,:)
  label{i} = sprintf('N and N-%d',i);
end
if strcmp(disttype,'spl')
  xlabel(ax,'Abs diff off12 (samples)');
  ylabel(ax,'Count per sample');
elseif strcmp(disttype,'km')
  xlabel(ax,'Abs diff E loc (km)');
  ylabel(ax,'Count per km');
end
legend(ax,p,label);
ax.XLim(1) = -binw;
ax.YAxis.Exponent = 4;
longticks(ax,2);
hold(ax,'off');

ax=f.ax(3);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
for i = m:-1:1
  aa=dlocxy{i};
  [N,edges]=histcounts(aa(:,2),'binwidth',binw,'normalization','countdensity');
%   edges = edges-binw/2;
%   [N,edges]=histcounts(aa(:,2),edges,'normalization','countdensity');
  stairs(ax,edges,[N N(end)],'color',color(i,:),'LineWidth',1);
  [MUHAT,SIGMAHAT] = normfit(aa(:,2));
  Y=normpdf(X,MUHAT,SIGMAHAT)*size(aa,1);
  plot(ax,X,Y,'Color',color(i,:),'linew',1);
  text(ax,0.02,0.5-i*0.05,sprintf('\\mu=%.2f, \\sigma=%.2f',MUHAT,SIGMAHAT),'Units','normalized',...
    'HorizontalAlignment','left'); %,'Color',color(i,:)
end
if strcmp(disttype,'spl')
  xlabel(ax,'Diff off13 (samples)');
  ylabel(ax,'Count per sample');
elseif strcmp(disttype,'km')
  xlabel(ax,'Diff N loc (km)');
  ylabel(ax,'Count per km');
end
axsym(ax,1);
% xlim(ax,[-50 50]);
ax.YAxis.Exponent = 3;
longticks(ax,2);
hold(ax,'off');

ax=f.ax(4);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
for i = m:-1:1
  aa=dlocxy{i};
  [N,edges]=histcounts(abs(aa(:,2)),'binwidth',binw,'normalization','countdensity');
%   edges = edges-binw/2;
%   [N,edges]=histcounts(abs(aa(:,2)),edges,'normalization','countdensity');
  Nn = N;
  if strcmp(disttype,'spl')
    Nn(1) = N(1)*2;
    Nn = Nn*(sum(N)/sum(Nn));
  end
  stairs(ax,edges,[Nn Nn(end)],'color',color(i,:),'LineWidth',1);
  plot(ax,[median(abs(aa(:,2))) median(abs(aa(:,2)))],ax.YLim,'--','Color',color(i,:));
  text(ax,0.98,0.5-i*0.05,sprintf('med=%.2f km',median(abs(aa(:,2)))),'Units','normalized',...
    'HorizontalAlignment','right'); %,'Color',color(i,:)
end
if strcmp(disttype,'spl')
  xlabel(ax,'Abs diff off13 (samples)');
  ylabel(ax,'Count per sample');
elseif strcmp(disttype,'km')
  xlabel(ax,'Abs diff N loc (km)');
  ylabel(ax,'Count per km');
end
ax.XLim(1) = -binw;
ax.YAxis.Exponent = 3;
longticks(ax,2);
hold(ax,'off');

if strcmp(disttype,'spl')
  ax=f.ax(5);
  hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
  for i = m:-1:1
    aa=dlocxy{i};
    [N,edges]=histcounts(aa(:,2)-aa(:,1),'binwidth',binw,'normalization','countdensity');
%     edges = edges-binw/2;
%     [N,edges]=histcounts(aa(:,2)-aa(:,1),edges,'normalization','countdensity');
    stairs(ax,edges,[N N(end)],'color',color(i,:),'LineWidth',1);
    [MUHAT,SIGMAHAT] = normfit(aa(:,2)-aa(:,1));
    Y=normpdf(X,MUHAT,SIGMAHAT)*size(aa,1);
    plot(ax,X,Y,'Color',color(i,:),'linew',1);
    text(ax,0.02,0.5-i*0.05,sprintf('\\mu=%.2f, \\sigma=%.2f',MUHAT,SIGMAHAT),'Units','normalized',...
      'HorizontalAlignment','left'); %,'Color',color(i,:)
  end
  xlabel(ax,'Diff off23 (samples)');
  ylabel(ax,'Count per sample');
  axsym(ax,1);
  longticks(ax,2);
  hold(ax,'off');

  ax=f.ax(6);
  hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
  for i = m:-1:1
    aa=dlocxy{i};
    [N,edges]=histcounts(abs(aa(:,2)-aa(:,1)),'binwidth',binw,'normalization','countdensity');
%     edges = edges-binw/2;
%     [N,edges]=histcounts(abs(aa(:,2)-aa(:,1)),edges,'normalization','countdensity');
    Nn = N;
    if strcmp(disttype,'spl')
      Nn(1) = N(1)*2;
      Nn = Nn*(sum(N)/sum(Nn));
    end
    stairs(ax,edges,[Nn Nn(end)],'color',color(i,:),'LineWidth',1);
    plot(ax,[median(abs(aa(:,2)-aa(:,1))) median(abs(aa(:,2)-aa(:,1)))],ax.YLim,'--','Color',color(i,:));
    text(ax,0.98,0.5-i*0.05,sprintf('med=%.2f',median(abs(aa(:,2)-aa(:,1)))),'Units','normalized',...
      'HorizontalAlignment','right'); %,'Color',color(i,:)
  end
  xlabel(ax,'Abs diff off23 (samples)');
  ylabel(ax,'Count per sample');
  ax.XLim(1) = -binw;
  longticks(ax,2);
  hold(ax,'off');
end


