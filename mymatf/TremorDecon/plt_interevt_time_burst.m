function [f] = plt_interevt_time_burst(burstt,intert,refjuldate,ymax,ttol1,ttol2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_interevt_time_burst(ttol1,ttol2,ymax,hfinter03,hfinter04,hfinter05)
%
% Different from 'plt_interevt_time_diffets.m' which plots all ETSs,
% this function is adaptive to plot the following for a single ETS in detail:
%
% 1.
%   if 'burstt' is not empty, plot the found tremor bursts as shaded area.
% 2.
%   if 'intert' is to plot the separation in time between itself and its 
%   preceding detection, i.e., inter-event time.
% 
% Therefore, you can plot either 1 or 2, or both, depending on if the input
% entry is empty of not
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/03/15
% Last modified date:   2022/03/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

defval('ymax',2e-3);
defval('ttol1',[]);
defval('ttol2',[]);

% evaluate the max range of time to determine rows of plots suitable
if ~isempty(intert)
  xran = [floor(min(intert(:,1))) ceil(max(intert(:,1)))];
else
  xran = [floor(min(burstt(:))) ceil(max(burstt(:)))];
end

ixlen = 0.5;  % each row means 1 day.
nrow = ceil(range(xran)/ixlen);
if nrow>4
  ixlen = 1;  % each row means 2 days
  nrow = ceil(range(xran)/ixlen);
end
ncol = 1;
widin = 8;
htin = 1.5*nrow;

%intialize figure
[f] = initfig(widin,htin,nrow,ncol);

%set optimal axis postions
pltxran = [0.1 0.9]; pltyran = [0.15 0.95];
pltxsep = 0.02; pltysep = 0.03; 
xpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

msizehf = 3;

for i = 1: nrow
  ax = f.ax(i);
  hold(ax,'on');
  if ~isempty(burstt)
    for j = 1: size(burstt,1)
      patarea = [burstt(j,1) -ymax;
        burstt(j,2) -ymax;
        burstt(j,2) ymax;
        burstt(j,1) ymax;
        burstt(j,1) -ymax];
      patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.3,'edgecolor','none');
    end
  end
  if ~isempty(intert)
    scatter(ax, intert(:,1), intert(:,2)*86400, msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w')
  end
  if ~isempty(ttol1)
    plot(ax, ax.XLim, [ttol1(1) ttol1(1)], 'k--');
  end
  if ~isempty(ttol2)
    plot(ax, ax.XLim, [ttol2(1) ttol2(1)], 'b--');
  end
  ylim(ax,[0, ymax]);
  xlim(ax,[xran(1)+(i-1)*ixlen xran(1)+i*ixlen]);
  hold(ax,'off');

end

if isequal(refjuldate, 2003060)
  stday = {'Mar. 1, 2003 (2003060)'};
elseif isequal(refjuldate, 2004194)
  stday = {'Jul. 12, 2004 (2004194)'};
elseif isequal(refjuldate, 2005254)
  stday = {'Sep. 11, 2005 (2005254)'};
end
ETS = num2str(floor(refjuldate/1000));
text(f.ax(1),0.92,0.9,ETS,'FontSize',10,'unit','normalized','horizontalalignment','left',...
  'EdgeColor','k','Margin',2);
ylabel(f.ax(nrow), 'Time (s) to the preceding detection','FontSize',10);
xlabel(f.ax(nrow), strcat({'Time (day) since '},stday),'FontSize',10);



