function [f,den1d,conmat,conobj]=plt_srcdlocinmap(dt,dloc,dtran,disttype,timetype,...
  msize,cstr,symbol,scale,bintype,xran,yran,dx,dy,contourflag,smoothsigma,intvl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_srcdlocinmap(dt,dloc,disttype,timetype,bintype)
%
% Function to create a similar plot to figure 4 of Rubin & Gillard (JGRSE,
% 2000) and figure 4 of Rubin (JGRSE, 2002). For each detection, it gets 
% placed at the origin, and then all the other detections in the adopted
% differential time bin (or overall within some range) get plotted in their
% relative position on the fault plane. Here we talk about the event pairs
% composed by either each detection and all the other (that occur within  
% some time range), or consecutive events 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/03/08
% Last modified date:   2024/01/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('dtran',[]);
defval('disttype','km');
defval('timetype','tarvl');
defval('scale','log10');
defval('bintype','grid');
defval('bintype','grid');
defval('dx',[]);
defval('dy',[]);
defval('contourflag',1);
defval('smoothsigma',[]);
defval('intvl',2);

% %convert time offset to relative loc, but we mainly need the travel time estimate from sources
% ftrans = 'interpchao';
% [imploc, ~] = off2space002(imp(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
% %%%option 1 to define 'relative origin time', time 0 is the origin time of a source at 0,0 and
% %%%arrive at the start of your window at PGC
% [imploc0, ~] = off2space002([0 0],sps,ftrans,0);  % a ref source at 0,0
% tcor = round((imploc(:,6)-imploc0(6))*sps);   % travel time difference between each source and a ref source at 0,0
% imptime = imp(:,1)-tcor;

htin = 5.5;   % maximum height allowed is 11 inches
nrow = 1;
if contourflag
  widin = 12;  % maximum width allowed is 8.5 inches
  ncol = 3;
else
  widin = 8;  % maximum width allowed is 8.5 inches
  ncol = 2;  
end
f = initfig(widin,htin,nrow,ncol); %initialize fig

figxran = [0.06 0.96]; figyran = [0.06 0.96];
figxsep = 0.05; figysep = 0.06;
optaxpos(f,nrow,ncol,figxran,figyran,figxsep,figysep);

% nbin = length(edges)-1;
% color = jet(nbin);
% binedge = (0: binwdist: 50*binwdist)';
% binwdt = 0.25;
% edges = -binwdt/2: binwdt: max(dt)+binwdt/2;
% nbin = length(edges)-1;
% dlocplt = dloc;
% dtplt = dt;
% for i = 1: nbin
  % dlocplt = dloc((dt)>=edges(i) & (dt)<edges(i+1),:);
  % dtplt = dt((dt)>=edges(i) & (dt)<edges(i+1),:);
  if ~isempty(dtran)
    dtplt = dt(dt>=dtran(1) & dt<=dtran(2),:);
    dlocplt = dloc(dt>=dtran(1) & dt<=dtran(2),:);
  else
    dtplt = dt;
    dlocplt = dloc;
  end
  
  dplt = [dtplt dlocplt];
  dpltst = sortrows(dplt,1,'descend');

  %scatter plot of all data points
  ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on'); axis(ax, 'equal');
  scatter(ax,dpltst(:,2),dpltst(:,3),msize,dpltst(:,1),'filled');
%   keyboard
  colormap(ax,'jet');
%   colormap(ax,flipud(colormap(ax,'kelicol')));
  c=colorbar(ax,'SouthOutside');
%   caxis(ax,[min(dpltst(:,1)) max(dpltst(:,1))]);
  ax.CLim(2) = prctile(dpltst(:,1),99);
  if strcmp(timetype,'tarvl')
    c.Label.String = strcat({'Differential arrival time (s)'});
  elseif strcmp(timetype,'tori')
    c.Label.String = strcat({'Differential origin time (s)'});
  end
  ax.GridLineStyle = '--';
  ax.XAxisLocation = 'top';
  xlim(ax,xran);
  ylim(ax,yran);
  if strcmp(disttype,'spl')
    xlabel(ax,'Diff off12 (samples)');
    ylabel(ax,'Diff off13 (samples)');
  elseif strcmp(disttype,'km')
    xlabel(ax,'Diff E loc (km)');
    ylabel(ax,'Diff N loc (km)');
  end
  % plot(ax,ax.XLim,ax.YLim,'k--','linew',1);
  longticks(ax,2);
  hold(ax,'off');

  %cumulative density, bin type determined by 'bintype'
  ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
  if strcmp(bintype,'pixel')
    den1d = density_pixel(dlocplt(:,1),dlocplt(:,2));
  elseif strcmp(bintype,'grid')
    den1d = density_matrix(dlocplt(:,1),dlocplt(:,2),xran,yran,dx,dy);
  end
  den1d = den1d(den1d(:,3)>0, :);
  den1d = sortrows(den1d,3);
  
  dum = den1d;
  dum(dum(:,3)>1, :) = [];
  if strcmp(scale,'log10')
    dum(:,3) = log10(dum(:,3));
  end
  scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),symbol,'linew',0.2);  %, 'MarkerEdgeColor', 'w')
  
  dum = sortrows(den1d,3);
  dum(dum(:,3)==1, :) = [];
  if strcmp(scale,'log10')
    dum(:,3) = log10(dum(:,3));
  end  
  scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),symbol,'filled','MarkerEdgeColor','none');
  text(ax,0.98,0.05,sprintf('%d events',size(dlocplt,1)),'Units','normalized',...
    'HorizontalAlignment','right','FontSize',9);
  colormap(ax,'jet');
%   colormap(ax,flipud(colormap(ax,'kelicol')));
  c=colorbar(ax,'SouthOutside');
  ax.CLim(2) = prctile(dum(:,3),99);
  if strcmp(scale,'log10')
    c.Label.String = strcat('log_{10}(',cstr{1},')');
  elseif strcmp(scale,'linear')
    c.Label.String = cstr{1};  
  end
  ax.GridLineStyle = '--';
  ax.XAxisLocation = 'top';
  xlim(ax,xran);
  ylim(ax,yran);
  % xticks(ax,xran(1):5:xran(2));
  % yticks(ax,yran(1):5:yran(2));
  if strcmp(disttype,'spl')
    xlabel(ax,'Diff off12 (samples)');
    ylabel(ax,'Diff off13 (samples)');
  elseif strcmp(disttype,'km')
    xlabel(ax,'Diff E loc (km)');
    ylabel(ax,'Diff N loc (km)');
  end
  % plot(ax,ax.XLim,ax.YLim,'k--','linew',1);
  longticks(ax,2);
  hold(ax,'off');

  if contourflag
    %contours of cumulative density
    ax=f.ax(3); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
    if strcmp(disttype,'spl')
      [~,xgrid,ygrid,zgrid] = ...
        zeropadmat2d(den1d,min(den1d(:,1)):1:max(den1d(:,1)),...
        min(den1d(:,2)):1:max(den1d(:,2)));
    elseif strcmp(disttype,'km')
      [~,xgrid,ygrid,zgrid] = ...
        zeropadmat2d(den1d,floor(min(den1d(:,1))):dx:ceil(max(den1d(:,1))),...
        floor(min(den1d(:,2))):dx:ceil(max(den1d(:,2))));
    end
    if ~isempty(smoothsigma)
      zgridgf = imgaussfilt(zgrid,smoothsigma);  %smooth it a bit
    else
      zgridgf = zgrid;
    end

%     %%%%%% test the effect of diff smothing sigma
%     f2 = initfig(12,9,3,4); %initialize fig
%     figxran = [0.06 0.96]; figyran = [0.06 0.96];
%     figxsep = 0.03; figysep = 0.03;
%     optaxpos(f2,3,4,figxran,figyran,figxsep,figysep);    
%     conplt=1:intvl:prctile(den1d(:,3),99);
%     ax=f2.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
%     [conmat,conobj] = contour(ax,xgrid,ygrid,zgrid,conplt,'-'); %,'ShowText','on','color',[.3 .3 .3]
% %     imagesc(ax,[xgrid(1,1) xgrid(1,end)], [ygrid(1,1) ygrid(end,1)],zgrid);
%     text(ax,0.9,0.9,'raw','Units','normalized','HorizontalAlignment','right');
%     colorbar(ax);
%     colormap(ax,'jet');    
%     smoothsigma = 0.5:0.2:2.5;
%     for i = 1: length(smoothsigma)
%       zgridgf = imgaussfilt(zgrid,smoothsigma(i));
%       ax=f2.ax(i+1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
%       [conmat,conobj] = contour(ax,xgrid,ygrid,zgridgf,conplt,'-'); %,'ShowText','on','color',[.3 .3 .3]
% %       imagesc(ax,[xgrid(1,1) xgrid(1,end)], [ygrid(1,1) ygrid(end,1)],zgridgf);
%       text(ax,0.9,0.9,sprintf('%.1f',smoothsigma(i)),'Units','normalized','HorizontalAlignment','right');
%       colorbar(ax);
%       colormap(ax,'jet');
%     end
%     %%%%%% test the effect of diff smothing sigma  
    
    % perc = 50:10:90;
    % conplt = prctile(dum(:,3),perc);
%     conplt=1:intvl:prctile(den1d(:,3),99);
    conplt=1:intvl:max(den1d(:,3));
    if strcmp(scale,'log10')
      zgridgf = log10(zgridgf);
    end
    [conmat,conobj] = contour(ax,xgrid,ygrid,zgridgf,conplt,'-'); %,'ShowText','on','color',[.3 .3 .3]
    if ~isempty(conmat) 
      contable = getContourLineCoordinates(conmat);
      conmat=table2array(contable);
    end
    colormap(ax,'jet');
%     colormap(ax,flipud(colormap(ax,'kelicol')));
    c=colorbar(ax,'SouthOutside');
    ax.CLim(2) = prctile(dum(:,3),99);
    if strcmp(scale,'log10')
      c.Label.String = strcat('log_{10}(',cstr{1},')');
    elseif strcmp(scale,'linear')
      c.Label.String = cstr{1};
    end
    ax.GridLineStyle = '--';
    ax.XAxisLocation = 'top';
    xlim(ax,xran);
    ylim(ax,yran);
    % xticks(ax,xran(1):5:xran(2));
    % yticks(ax,yran(1):5:yran(2));
    if strcmp(disttype,'spl')
      xlabel(ax,'Diff off12 (samples)');
      ylabel(ax,'Diff off13 (samples)');
    elseif strcmp(disttype,'km')
      xlabel(ax,'Diff E loc (km)');
      ylabel(ax,'Diff N loc (km)');
    end
    longticks(ax,2);
    hold(ax,'off');
  end

% keyboard





