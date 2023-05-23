function f=plt_srcdlocinmap(dt2all,dloc2all,sps,disttype,timetype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_srcdlocinmap(dt2all,dloc2all,sps,disttype,timetype)
%
% Function to create a similar plot to figure 4 of Rubin & Gillard (JGRSE,
% 2000) and figure 4 of Rubin (JGRSE, 2002). For each detection, it gets 
% placed at the origin, and then all the other detections in the adopted
% differential time bin (or overall within some range) get plotted in their
% relative position on the fault plane. Here we talk about the event pairs
% composed by each detection and all the other (that occur within some 
% time range).
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/03/08
% Last modified date:   2023/03/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('disttype','spl');
defval('timetype','tarvl');

% %convert time offset to relative loc, but we mainly need the travel time estimate from sources
% ftrans = 'interpchao';
% [imploc, ~] = off2space002(imp(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
% %%%option 1 to define 'relative origin time', time 0 is the origin time of a source at 0,0 and
% %%%arrive at the start of your window at PGC
% [imploc0, ~] = off2space002([0 0],sps,ftrans,0);  % a ref source at 0,0
% tcor = round((imploc(:,6)-imploc0(6))*sps);   % travel time difference between each source and a ref source at 0,0
% imptime = imp(:,1)-tcor;

widin = 12;  % maximum width allowed is 8.5 inches
htin = 6;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f = initfig(widin,htin,nrow,ncol); %initialize fig

xran = [0.1 0.96]; yran = [0.06 0.96];
xsep = 0.1; ysep = 0.08;
optaxpos(f,nrow,ncol,xran,yran,xsep,ysep);

% nbin = length(edges)-1;
% color = jet(nbin);
% binedge = (0: binwdist: 50*binwdist)';
binwdt = 0.25;
edges = -binwdt/2: binwdt: max(dt2all/sps)+binwdt/2;
nbin = length(edges)-1;
dlocplt = dloc2all;
dtplt = dt2all;
for i = 1: nbin
  dlocplt = dloc2all((dt2all/sps)>=edges(i) & (dt2all/sps)<edges(i+1),:);
  dtplt = dt2all((dt2all/sps)>=edges(i) & (dt2all/sps)<edges(i+1),:);
  
  dplt = [dtplt dlocplt];
  dpltst = dplt;
  dpltst = sortrows(dplt,1,'descend');

  ax=f.ax(1);
  hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
  scatter(ax,dpltst(:,2),dpltst(:,3),5,dpltst(:,1)/sps,'filled');
  colormap(ax,'jet');
  c=colorbar(ax,'east');
  if strcmp(timetype,'tarvl')
    c.Label.String = strcat({'Differential arrival time (s)'});
  elseif strcmp(timetype,'tori')
    c.Label.String = strcat({'Differential origin time (s)'});
  end
  caxis(ax,[min(dtplt/sps) max(dtplt/sps)]);
  axis(ax, 'equal');
  if strcmp(disttype,'spl')
    xlabel(ax,'Diff off12 (samples)');
    ylabel(ax,'Diff off13 (samples)');
    axis(ax,[-50 50 -50 50]);
  elseif strcmp(disttype,'km')
    xlabel(ax,'Diff E loc (km)');
    ylabel(ax,'Diff N loc (km)');
    axis(ax,[-4 4 -4 4]);
  end
  plot(ax,ax.XLim,ax.YLim,'k--','linew',1);
  hold(ax,'off');

  ax=f.ax(2);
  hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
  density1d = density_pixel(dpltst(:,2),dpltst(:,3));
  dum = density1d(density1d(:,3)>0, :);
  dum(dum(:,3)>1, :) = [];
  scatter(ax,dum(:,1),dum(:,2),5,log(dum(:,3)),'o','linew',0.2);  %, 'MarkerEdgeColor', 'w')
  dum = sortrows(density1d(density1d(:,3)>0, :), 3);
  dum(dum(:,3)==1, :) = [];
  scatter(ax,dum(:,1),dum(:,2),5,log(dum(:,3)),'o','filled','MarkerEdgeColor','none');  %,
  colormap(ax,'jet');
  c=colorbar(ax,'east');
  c.Label.String = strcat({'log_{10}(# of detections / pixel)'});
  axis(ax, 'equal');
  if strcmp(disttype,'spl')
    xlabel(ax,'Diff off12 (samples)');
    ylabel(ax,'Diff off13 (samples)');
    axis(ax,[-50 50 -50 50]);
  elseif strcmp(disttype,'km')
    xlabel(ax,'Diff E loc (km)');
    ylabel(ax,'Diff N loc (km)');
    axis(ax,[-4 4 -4 4]);
  end
  plot(ax,ax.XLim,ax.YLim,'k--','linew',1);
  hold(ax,'off');

end


% keyboard





