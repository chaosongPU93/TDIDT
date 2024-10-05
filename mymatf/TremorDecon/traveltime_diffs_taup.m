% traveltime_diffs_taup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to use the customized velocity model Puget Sound P3 
% velocity model created by 'taupcreate' to compute the theoretical S-wave 
% travel time difference at the Cascadia stations, given a set of source 
% locations. 
% --We would like to know, for example, if srcs changes in loc by 1 km in diff
%   directions, what is their arrival time / travel time difference of S-wave
%   between different stations, e.g., change in off12, off13, off14 for change
%   in src loc.
% --Basically, if you have a set of src, then you can obtain the travel time
%   using 'taup', then you will have all time differences.
% --Using these time offsets, you can ask the question slightly differently,
%   how much change in src loc is needed to change the time offset by a quarter 
%   period of seismogram (or characteristic duration of LFE templates), which is
%   about 0.25 s. A quarter is 0.0625 s. A quarter of period change in offset
%   roughly would shift from max. CC to 0 (ie, sources would have deconstructive
%   interference), and half would shift to min CC. 
% --See also 'traveltime_diffps_taup.m' and 'traveltime_diff_taupvshypo' for 
%   different purposes.
% 
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2024/04/13
% Last modified date:   2024/04/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clc;
clear;
close all;

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

sps = 160;
ftrans = 'interpchao';


%% the currently used location grid, inversion results from hypoinverse
%location of family 002, and the origin
loc0 = off2space002([0 0],sps,ftrans,0);  % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

%location of all grid points
allgrid = off2space002([],sps,ftrans,0);

%% a bunch of change in relative locations to the family 002
%in E direction
dx1 = (-3:0.5:3)'; %in E
dy1 = zeros(size(dx1));
% dxy1 = [dx1 dy1];

%in N direction
dy2 = (-3:0.5:3)'; %in N
dx2 = zeros(size(dy2));
% dxy2 = [dx2 dy2];

%in SE direction, but expressed in N and E coordinate frame
[dx3,dy3] = coordinate_rot(dx1,dy1,45);
% dxy3 = [dx3 dy3];

%in NE direction, but expressed in N and E coordinate frame
[dx4,dy4] = coordinate_rot(dx2,dy2,45);
% dxy4 = [dx4 dy4];

% %%%plot of source location
% f=initfig;
% ax = f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
% axis(ax, 'equal');
% scatter(ax,dx1,dy1,40,'r','filled');
% scatter(ax,dx2,dy2,40,'b','filled');
% scatter(ax,dx3,dy3,40,'k','filled');
% scatter(ax,dx4,dy4,40,'c','filled');
% xlim(ax,xran);
% ylim(ax,yran);
% xlabel(ax,'E (km)');
% ylabel(ax,'N (km)');


%% absolute locations including depth
%obtain the lon and lat
[lon1,lat1] = relaloc2absloc(dx1,dy1,loc0(3),loc0(4));
[lon2,lat2] = relaloc2absloc(dx2,dy2,loc0(3),loc0(4));
[lon3,lat3] = relaloc2absloc(dx3,dy3,loc0(3),loc0(4));
[lon4,lat4] = relaloc2absloc(dx4,dy4,loc0(3),loc0(4));

%obtain the lon and lat, replying on scatter interpolation
F = scatteredInterpolant(allgrid(:,3),allgrid(:,4),allgrid(:,5),'linear','none');
dep1 = F(lon1,lat1);
dep2 = F(lon2,lat2);
dep3 = F(lon3,lat3);
dep4 = F(lon4,lat4);

% %%%plot of source location colorcoded by depth
% f=initfig;
% ax = f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
% % axis(ax, 'equal');
% scatter(ax,lon1,lat1,40,dep1,'filled');
% scatter(ax,lon2,lat2,40,dep2,'filled');
% scatter(ax,lon3,lat3,40,dep3,'filled');
% scatter(ax,lon4,lat4,40,dep4,'filled');
% colormap(ax,'viridis'); 
% c=colorbar(ax);
% c.Label.String = 'Depth';
% daspect(ax,[1 cos(loc0(4)*pi/180) 1]);
% % daspect(ax,[cos(loc0(4)*pi/180) 1 1]);
% xlabel(ax,'Lon');
% ylabel(ax,'Lat');


%% compute the travel time difference using Taup at all stations
stas=['PGC  '
      'SSIB '
      'SILB '
      'KLNB '
      'LZB  '
      'TWKB '
      'MGCB '
      ];
nsta = size(stas,1);

fid = fopen(fullfile('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary','staloc.txt'), 'r');
datalist = textscan(fid, '%s %s %s \n');  % 3 cell arrays, station name, lat, lon 
stanm = char(length(datalist{1}), 5);  % station name array
% stalist = [];
for i = 1: length(datalist{1})
  str = datalist{1}{i};
  if length(str)==3
    str = strcat(str,{'  '});
  elseif length(str)==4
    str = strcat(str,{' '});
  end
  stanm(i,1:5) = char(str);
end

staloc = [];  % station loc, lon, lat
for ista = 1: size(stas,1)
  [~,idx] = ismember(stas(ista,:),stanm,'rows');
  staloc(ista,:) = [str2double(datalist{3}{idx}) str2double(datalist{2}{idx})];
end

pugetp3 = '/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/pugetp3.taup';  % pre-created vel model for Taup

tts1 = zeros(length(lon1), nsta);  % S-wave travel time
tts2 = zeros(length(lon2), nsta);  % S-wave travel time
tts3 = zeros(length(lon3), nsta);  % S-wave travel time
tts4 = zeros(length(lon4), nsta);  % S-wave travel time
for j = 1: nsta
  %for src locs varying in E direction
  for i = 1: length(lon1)
    tt=tauptime('mod',pugetp3,'dep',dep1(i,1),'ph','s,S','evt',[lat1(i,1) lon1(i,1)],...
      'sta',[staloc(j,2) staloc(j,1)]);
    if ~strcmp(tt(1).phase,'s')
      fprintf('At i=%d, j=%d, s take-off angle might be >90 degrees, ie., downgoing wave \n.',...
        i,j);
      tts1(i,j) = 0;   % choose the first S arrival
    else
      tts1(i,j) = tt(1).time;   % choose the first S arrival
    end
  end
  
  %for src locs varying in N direction
  for i = 1: length(lon2)
    tt=tauptime('mod',pugetp3,'dep',dep2(i,1),'ph','s,S','evt',[lat2(i,1) lon2(i,1)],...
      'sta',[staloc(j,2) staloc(j,1)]);
    if ~strcmp(tt(1).phase,'s')
      fprintf('At i=%d, j=%d, s take-off angle might be >90 degrees, ie., downgoing wave \n.',...
        i,j);
      tts2(i,j) = 0;   % choose the first S arrival
    else
      tts2(i,j) = tt(1).time;   % choose the first S arrival
    end
  end
  
  %for src locs varying in SE direction
  for i = 1: length(lon3)
    tt=tauptime('mod',pugetp3,'dep',dep3(i,1),'ph','s,S','evt',[lat3(i,1) lon3(i,1)],...
      'sta',[staloc(j,2) staloc(j,1)]);
    if ~strcmp(tt(1).phase,'s')
      fprintf('At i=%d, j=%d, s take-off angle might be >90 degrees, ie., downgoing wave \n.',...
        i,j);
      tts3(i,j) = 0;   % choose the first S arrival
    else
      tts3(i,j) = tt(1).time;   % choose the first S arrival
    end
  end
  
  %for src locs varying in NE direction
  for i = 1: length(lon4)
    tt=tauptime('mod',pugetp3,'dep',dep4(i,1),'ph','s,S','evt',[lat4(i,1) lon4(i,1)],...
      'sta',[staloc(j,2) staloc(j,1)]);
    if ~strcmp(tt(1).phase,'s')
      fprintf('At i=%d, j=%d, s take-off angle might be >90 degrees, ie., downgoing wave \n.',...
        i,j);
      tts4(i,j) = 0;   % choose the first S arrival
    else
      tts4(i,j) = tt(1).time;   % choose the first S arrival
    end
  end  
end

%% Travel time from each src to different stations 
% nrow = 3;
% ncol = 3;
% widin = 8.4;  % maximum width allowed is 8.5 inches
% htin = 8;   % maximum height allowed is 11 inches
% f = initfig(widin,htin,nrow,ncol);
% pltxran = [0.1 0.95]; pltyran = [0.1 0.95];
% pltxsep = 0.03; pltysep = 0.03; 
% axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
% delete(f.ax(8:9));

% for i = 1: nsta
%   ax = f.ax(i); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%   axis(ax, 'equal');
%   scatter(ax,dx1,dy1,50,tts1(:,i),'filled');
%   scatter(ax,dx2,dy2,50,tts2(:,i),'filled');
%   scatter(ax,dx3,dy3,50,tts3(:,i),'filled');
%   scatter(ax,dx4,dy4,50,tts4(:,i),'filled');
%   colormap(ax,'viridis'); 
%   colorbar(ax);
%   title(ax,stas(i,:));
%   xran = [-4 4];
%   yran = [-4 4];
%   xlim(ax,xran);
%   ylim(ax,yran);
% end
% supertit(f.ax(1:nsta),'Travel time');

%% Travel time gradient at the same station
dtts1 = tts1-tts1(floor(length(lon1)/2)+1, :);
dtts2 = tts2-tts2(floor(length(lon1)/2)+1, :);
dtts3 = tts3-tts3(floor(length(lon1)/2)+1, :);
dtts4 = tts4-tts4(floor(length(lon1)/2)+1, :);
%gradient, doff/dloc
dtts1dloc = dtts1./dx1;
dtts2dloc = dtts2./dx1;
dtts3dloc = dtts3./dx1;
dtts4dloc = dtts4./dx1;
indnan=find(isnan(dtts1dloc(:,1)));
ind=setdiff(1:length(dx1), indnan);

%median of gradient
dttsdloc = [median(dtts1dloc(ind,:),1); median(dtts2dloc(ind,:),1); ...
               median(dtts3dloc(ind,:),1); median(dtts4dloc(ind,:),1)]; 

% nrow = 3;
% ncol = 3;
% widin = 8.4;  % maximum width allowed is 8.5 inches
% htin = 8;   % maximum height allowed is 11 inches
% f = initfig(widin,htin,nrow,ncol);
% pltxran = [0.1 0.95]; pltyran = [0.1 0.95];
% pltxsep = 0.03; pltysep = 0.03; 
% axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
% delete(f.ax(8:9));
% 
% for i = 1: nsta
%   ax = f.ax(i); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%   axis(ax, 'equal');
%   scatter(ax,dx1,dy1,50,dtts1(:,i),'filled');
%   scatter(ax,dx2,dy2,50,dtts2(:,i),'filled');
%   scatter(ax,dx3,dy3,50,dtts3(:,i),'filled');
%   scatter(ax,dx4,dy4,50,dtts4(:,i),'filled');
% %   colormap(ax,'viridis'); 
%   colormap('bluewhitered'); 
%   colorbar(ax);
%   title(ax,stas(i,:));
%   xran = [-4 4];
%   yran = [-4 4];
%   xlim(ax,xran);
%   ylim(ax,yran);
% end
% supertit(f.ax(1:nsta),'Travel time change wrt location');
% 
% nrow = 3;
% ncol = 3;
% widin = 8.4;  % maximum width allowed is 8.5 inches
% htin = 8;   % maximum height allowed is 11 inches
% f = initfig(widin,htin,nrow,ncol);
% pltxran = [0.1 0.95]; pltyran = [0.1 0.95];
% pltxsep = 0.03; pltysep = 0.03; 
% axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
% delete(f.ax(8:9));
% 
% for i = 1: nsta
%   ax = f.ax(i); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%   axis(ax, 'equal');
%   scatter(ax,dx1(ind),dy1(ind),50,dtts1dloc(ind,i),'filled');
%   scatter(ax,dx2(ind),dy2(ind),50,dtts2dloc(ind,i),'filled');
%   scatter(ax,dx3(ind),dy3(ind),50,dtts3dloc(ind,i),'filled');
%   scatter(ax,dx4(ind),dy4(ind),50,dtts4dloc(ind,i),'filled');
%   text(ax,3.0,0.8,sprintf('%.3f',dttsdloc(1,i)),'fontsize',9,...
%     'HorizontalAlignment','center');
%   text(ax,0,3.6,sprintf('%.3f s/km',dttsdloc(2,i)),'fontsize',9,...
%     'HorizontalAlignment','center');
%   text(ax,3.0,-3.0,sprintf('%.3f',dttsdloc(3,i)),'fontsize',9,...
%     'HorizontalAlignment','center');
%   text(ax,3.0,3.0,sprintf('%.3f',dttsdloc(4,i)),'fontsize',9,...
%     'HorizontalAlignment','center');  
%   colormap(ax,'viridis'); 
%   colorbar(ax);
%   title(ax,stas(i,:));
%   xran = [-4 4];
%   yran = [-4 4];
%   xlim(ax,xran);
%   ylim(ax,yran);
% end
% supertit(f.ax(1:nsta),'Travel time gradient');

%% Change in travel time difference (doff1i) wrt origin
%%%Travel time difference between sta 1 and stat i from each src 
%%%this is what we called 'off1i'

%for each src loc, obtain the travel time diff between sta 1 and sta i, PGC
%following the convention: off1i = ttrvl1 - ttrvli
off1i1 = tts1(:,1)-tts1;
off1i2 = tts2(:,1)-tts2;
off1i3 = tts3(:,1)-tts3;
off1i4 = tts4(:,1)-tts4;

% nrow = 2;
% ncol = 3;
% widin = 8.4;  % maximum width allowed is 8.5 inches
% htin = 6;   % maximum height allowed is 11 inches
% f = initfig(widin,htin,nrow,ncol);
% pltxran = [0.1 0.95]; pltyran = [0.1 0.95];
% pltxsep = 0.03; pltysep = 0.03; 
% optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
% 
% for i = 1: nsta-1
%   ax = f.ax(i); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%   axis(ax, 'equal');
%   scatter(ax,dx1,dy1,50,off1i1(:,i+1),'filled');
%   scatter(ax,dx2,dy2,50,off1i2(:,i+1),'filled');
%   scatter(ax,dx3,dy3,50,off1i3(:,i+1),'filled');
%   scatter(ax,dx4,dy4,50,off1i4(:,i+1),'filled');
%   colormap(ax,'viridis'); 
% %   colormap(ax,flipud(colormap(ax,'kelicol')));
%   colorbar(ax);
%   title(ax,strcat(stas(1,:),'--',stas(i+1,:)));
%   xlim(ax,xran);
%   ylim(ax,yran);
% end
% supertit(f.ax,'Travel time difference (t1-ti, off1i)');


%%%Change in Travel time difference between sta 1 and stat i wrt the change in
%%%src (src j - origin), this is what we called 'doff1i' 
%between same sta pair, diff in offset wrt diff in src loc
doff1i1 = off1i1-off1i1(floor(length(lon1)/2)+1, :);
doff1i2 = off1i2-off1i2(floor(length(lon2)/2)+1, :);
doff1i3 = off1i3-off1i3(floor(length(lon3)/2)+1, :);
doff1i4 = off1i4-off1i4(floor(length(lon4)/2)+1, :);

% nrow = 2;
% ncol = 3;
% widin = 8.4;  % maximum width allowed is 8.5 inches
% htin = 6;   % maximum height allowed is 11 inches
% f = initfig(widin,htin,nrow,ncol);
% pltxran = [0.1 0.95]; pltyran = [0.1 0.95];
% pltxsep = 0.03; pltysep = 0.03; 
% optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
% for i = 1: nsta-1
%   ax = f.ax(i); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%   axis(ax, 'equal');
%   scatter(ax,dx1,dy1,50,doff1i1(:,i+1),'filled');
%   scatter(ax,dx2,dy2,50,doff1i2(:,i+1),'filled');
%   scatter(ax,dx3,dy3,50,doff1i3(:,i+1),'filled');
%   scatter(ax,dx4,dy4,50,doff1i4(:,i+1),'filled');
% %   colormap(ax,'viridis'); 
% %   colormap(ax,flipud(colormap(ax,'kelicol')));
% %   colormap(ax,'redblue'); 
%   colormap('bluewhitered'); 
%   c=colorbar(ax);
%   title(ax,strcat(stas(1,:),'--',stas(i+1,:)));
%   cm = max(abs(c.Limits));
%   caxis(ax,[-cm cm]);
%   xran = [-4 4];
%   yran = [-4 4];
%   xlim(ax,xran);
%   ylim(ax,yran);
% end
% supertit(f.ax,'Change in travel time difference (doff1i) wrt origin');

%% Gradient of travel time difference (doff1i/dloc)
%%%Gradient of travel time difference between sta 1 and stat i. ie, this is what
%%%we called 'doff1i/dloc' 
%gradient, doff/dloc
doff1i1dloc = doff1i1./dx1;
doff1i2dloc = doff1i2./dx1;
doff1i3dloc = doff1i3./dx1;
doff1i4dloc = doff1i4./dx1;
indnan=find(isnan(doff1i1dloc(:,1)));
ind=setdiff(1:length(dx1), indnan);

% nrow = 2;
% ncol = 3;
% widin = 8.4;  % maximum width allowed is 8.5 inches
% htin = 6;   % maximum height allowed is 11 inches
% f = initfig(widin,htin,nrow,ncol);
% pltxran = [0.1 0.95]; pltyran = [0.1 0.95];
% pltxsep = 0.03; pltysep = 0.03; 
% optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
% for i = 1: nsta-1
%   ax = f.ax(i); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%   axis(ax, 'equal');
%   xlim(ax,xran);
%   ylim(ax,yran);
%   [rotx, roty] = complex_rot(0,yran(2),-0);
%   xarrow = [0-rotx 0+rotx];
%   yarrow = [0-roty 0+roty];
%   p=annotation('arrow','color','k','linestyle','-','linewidth',1);
%   p.Parent = ax;
%   p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
%   [rotx, roty] = complex_rot(0,yran(2),-90);
%   xarrow = [0-rotx 0+rotx];
%   yarrow = [0-roty 0+roty];
%   p=annotation('arrow','color','k','linestyle','-','linewidth',1);
%   p.Parent = ax;
%   p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
%   [rotx, roty] = complex_rot(0,yran(2),-45);
%   xarrow = [0-rotx 0+rotx];
%   yarrow = [0-roty 0+roty];
%   p=annotation('arrow','color','k','linestyle','--','linewidth',1);
%   p.Parent = ax;
%   p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
%   [rotx, roty] = complex_rot(0,yran(2),-135);
%   xarrow = [0-rotx 0+rotx];
%   yarrow = [0-roty 0+roty];
%   p=annotation('arrow','color','k','linestyle','--','linewidth',1);
%   p.Parent = ax;
%   p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
%   scatter(ax,dx1(ind),dy1(ind),50,doff1i1dloc(ind,i+1),'filled');
%   scatter(ax,dx2(ind),dy2(ind),50,doff1i2dloc(ind,i+1),'filled');
%   scatter(ax,dx3(ind),dy3(ind),50,doff1i3dloc(ind,i+1),'filled');
%   scatter(ax,dx4(ind),dy4(ind),50,doff1i4dloc(ind,i+1),'filled');
%   
%   colormap(ax,'viridis'); 
% %   colormap(ax,flipud(colormap(ax,'kelicol')));
% %   colormap(ax,'redblue'); 
% %   colormap('bluewhitered'); 
%   c=colorbar(ax);
%   title(ax,strcat(stas(1,:),'--',stas(i+1,:)));
% %   cm = max(abs(c.Limits));
% %   caxis(ax,[-cm cm]);
% end
% supertit(f.ax,'Gradient of travel time difference (doff1i/dloc)');

%%%instead of plotting in terms color, just use the numbers, how many quarters
%%%of period of tremor seismogram, 0.25/4 s
doff1i1dlocp = doff1i1./dx1/(0.25/4);
doff1i2dlocp = doff1i2./dx1/(0.25/4);
doff1i3dlocp = doff1i3./dx1/(0.25/4);
doff1i4dlocp = doff1i4./dx1/(0.25/4);


%% re-arrange like we do in the envelope CC
%%%between trio stations, 1/2/3
off1231 = [tts1(:,1)-tts1(:,2) tts1(:,1)-tts1(:,3) tts1(:,2)-tts1(:,3)];
off1232 = [tts2(:,1)-tts2(:,2) tts2(:,1)-tts2(:,3) tts2(:,2)-tts2(:,3)];
off1233 = [tts3(:,1)-tts3(:,2) tts3(:,1)-tts3(:,3) tts3(:,2)-tts3(:,3)];
off1234 = [tts4(:,1)-tts4(:,2) tts4(:,1)-tts4(:,3) tts4(:,2)-tts4(:,3)];

%%%between each of the 4th station, 4/5/6/7 and each of trio, 1/2/3
k = 0;
for i = 1:3
  for j = 4:nsta
    k=k+1;
    off141(:,k) = tts1(:,i)-tts1(:,j);
    off142(:,k) = tts2(:,i)-tts2(:,j);
    off143(:,k) = tts3(:,i)-tts3(:,j);
    off144(:,k) = tts4(:,i)-tts4(:,j);
  end
end

%between same sta pair, diff in offset wrt diff in src loc
doff1231 = off1231-off1231(floor(length(lon1)/2)+1, :);
doff1232 = off1232-off1232(floor(length(lon1)/2)+1, :);
doff1233 = off1233-off1233(floor(length(lon1)/2)+1, :);
doff1234 = off1234-off1234(floor(length(lon1)/2)+1, :);
doff141 = off141-off141(floor(length(lon1)/2)+1, :);
doff142 = off142-off142(floor(length(lon1)/2)+1, :);
doff143 = off143-off143(floor(length(lon1)/2)+1, :);
doff144 = off144-off144(floor(length(lon1)/2)+1, :);

%between same sta pair, GRADIENT in offset wrt diff in src loc
doff1231dloc = doff1231./dx1;
doff1232dloc = doff1232./dx1;
doff1233dloc = doff1233./dx1;
doff1234dloc = doff1234./dx1;
doff141dloc = doff141./dx1;
doff142dloc = doff142./dx1;
doff143dloc = doff143./dx1;
doff144dloc = doff144./dx1;
indnan=find(isnan(doff1i1dloc(:,1)));
ind=setdiff(1:length(dx1), indnan);

%%%looks like the gradient is pretty stable, can be represented by a median
doff123dloc = [median(doff1231dloc(ind,:),1); median(doff1232dloc(ind,:),1); ...
               median(doff1233dloc(ind,:),1); median(doff1234dloc(ind,:),1)];

doff14dloc = [median(doff141dloc(ind,:),1); median(doff142dloc(ind,:),1); ...
              median(doff143dloc(ind,:),1); median(doff144dloc(ind,:),1)];

%%%the max. gradient in terms of magnitude and direction
%results using E-N and SE-NE should be close
doff123dlocm = [sqrt(doff123dloc(1,:).^2+doff123dloc(2,:).^2); ...
                sqrt(doff123dloc(3,:).^2+doff123dloc(4,:).^2)]; %magnitude
doff123dloca = [90-atan2d(doff123dloc(2,:), doff123dloc(1,:)); ...
                90+45-atan2d(doff123dloc(4,:), doff123dloc(3,:))]; %direction

doff14dlocm = [sqrt(doff14dloc(1,:).^2+doff14dloc(2,:).^2); ...
               sqrt(doff14dloc(3,:).^2+doff14dloc(4,:).^2)]; %magnitude

doff14dloca = [90-atan2d(doff14dloc(2,:), doff14dloc(1,:)); ...
               90+45-atan2d(doff14dloc(4,:), doff14dloc(3,:))]; %direction


%%%instead of plotting in terms color, just use the numbers, how many quarters
%%%of period of tremor seismogram, 0.25/4 s
doff1231dlocp = doff1231dloc/(0.25/4);
doff1232dlocp = doff1232dloc/(0.25/4);
doff1233dlocp = doff1233dloc/(0.25/4);
doff1234dlocp = doff1234dloc/(0.25/4);
doff141dlocp = doff141dloc/(0.25/4);
doff142dlocp = doff142dloc/(0.25/4);
doff143dlocp = doff143dloc/(0.25/4);
doff144dlocp = doff144dloc/(0.25/4);

%%%in other words, how much location change is needed for a quarter of period 
%%%change in the travel time difference?
dlocdoff123p1d = 1./doff1231dlocp;
dlocdoff123p2d = 1./doff1232dlocp;
dlocdoff123p3d = 1./doff1233dlocp;
dlocdoff123p4d = 1./doff1234dlocp;
dlocdoff14p1d = 1./doff141dlocp;
dlocdoff14p2d = 1./doff142dlocp;
dlocdoff14p3d = 1./doff143dlocp;
dlocdoff14p4d = 1./doff144dlocp;

%%%looks like the gradient is pretty stable, can be represented by a median
doff123dlocp = [median(doff1231dlocp(ind,:),1); median(doff1232dlocp(ind,:),1); ...
                median(doff1233dlocp(ind,:),1); median(doff1234dlocp(ind,:),1)]; 

doff14dlocp = [median(doff141dlocp(ind,:),1); median(doff142dlocp(ind,:),1); ...
               median(doff143dlocp(ind,:),1); median(doff144dlocp(ind,:),1)];
             
%%%in other words, how much location change is needed for a quarter of period 
%%%change in the travel time difference?
dlocdoff123p = 1./doff123dlocp;
dlocdoff14p = 1./doff14dlocp;

%%%the min change in loc that results in a quarter period change in offset
%results using E-N and SE-NE should be close
%do NOT re-compute
dlocdoff123pm = (0.25/4)./doff123dlocm; %magnitude
dlocdoff123pa = doff123dloca; %direction

dlocdoff14pm = (0.25/4)./doff14dlocm; %magnitude
dlocdoff14pa = doff14dloca; %direction

%% max gradient of travel time difference (doff1i/dloc), mag and direcction 
nrow = 1;
ncol = 2;
widin = ncol*2.2;
htin = nrow*2.4;
pltxran = [0.10 0.96]; pltyran = [0.15 0.96]; % optimal axis location
pltxsep = 0.025; pltysep = 0.02;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

xran = [-3 3];
yran = [-3 3];

%for trio stations
ax = f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
axis(ax, 'equal');
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,-4:1:4);
yticks(ax,-4:1:4);
color = ['r';'b';'k';'c'];
refscl=1;
for i = 1:3
% i=1
  mag = doff123dlocm(2,i)/(0.25/4);
  ang = doff123dloca(2,i);
  [rotx, roty] = complex_rot(0,refscl*mag,-ang);
  xarrow = [0 0+rotx];
  yarrow = [0 0+roty];
  if i==3
    a=annotation('arrow','color',color(i,:),'linestyle',':','linewidth',1,'HeadLength',6,...
      'HeadWidth',6);
  else
    a=annotation('arrow','color',color(i,:),'linestyle','-','linewidth',1,'HeadLength',6,...
      'HeadWidth',6);
  end
  a.Parent = ax;
  a.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
  if i==3
    label{i} = sprintf('%s-%s',strtrim(stas(2,:)),strtrim(stas(3,:)));
  else
    label{i} = sprintf('%s-%s',strtrim(stas(1,:)),strtrim(stas(i+1,:)));
  end
  text(ax,0.02,0.25-(i-1)*0.1,label{i},'color',color(i,:),'Units','normalized',...
    'HorizontalAlignment','left');
end
%reference arrow, a quarter period
[rotx, roty] = complex_rot(0,refscl,-90);
xarrow = [1 1+rotx];
yarrow = [2 2+roty];
a=annotation('arrow','color',[.5 .5 .5],'linestyle','-','linewidth',1,'HeadLength',6,...
  'HeadWidth',6);
a.Parent = ax;
a.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
text(ax,0.5,2.5,'0.0625 s/km','fontsize',9,'HorizontalAlignment','left');
text(ax,0.02,0.92,'a','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
% legend(ax,'PGC-SSIB','PGC-SILB','SSIB-SILB','location','west');
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');

%for 4th station & 1st station
ax = f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
axis(ax, 'equal');
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,-4:1:4);
yticks(ax,-4:1:4);
label=[];
for i = 1:4
  mag = doff14dlocm(2,i)/(0.25/4);
  ang = doff14dloca(2,i);
  [rotx, roty] = complex_rot(0,refscl*mag,-ang);
  xarrow = [0 0+rotx];
  yarrow = [0 0+roty];
  a=annotation('arrow','color',color(i,:),'linestyle','-','linewidth',1,'HeadLength',6,...
    'HeadWidth',6);
  a.Parent = ax;
  a.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
  label{i} = sprintf('%s-%s',strtrim(stas(1,:)),strtrim(stas(i+3,:)));
  text(ax,0.55,0.35-(i-1)*0.1,label{i},'color',color(i,:),'Units','normalized',...
    'HorizontalAlignment','left');
end  
%reference arrow, a quarter period
[rotx, roty] = complex_rot(0,refscl,-90);
xarrow = [1 1+rotx];
yarrow = [2 2+roty];
a=annotation('arrow','color',[.5 .5 .5],'linestyle','-','linewidth',1,'HeadLength',6,...
  'HeadWidth',6);
a.Parent = ax;
a.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
text(ax,0.5,2.5,'0.0625 s/km','fontsize',9,'HorizontalAlignment','left');
text(ax,0.02,0.92,'b','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
% legend(ax,label,'Location','best');
xlabel(ax,'E (km)');
% ylabel(ax,'N (km)');
nolabels(ax,2);

fname = strcat('maxdoffdloc.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
% keyboard

%% min km change in location would lead to a quarter of a period change in 'doff'
nrow = 1;
ncol = 2;
widin = ncol*2.2;
htin = nrow*2.4;
pltxran = [0.10 0.96]; pltyran = [0.15 0.96]; % optimal axis location
pltxsep = 0.025; pltysep = 0.02;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

xran = [-2 2];
yran = [-2 2];

%for trio stations
ax = f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
axis(ax, 'equal');
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,-4:1:4);
yticks(ax,-4:1:4);
color = ['r';'b';'k';'c'];
refscl = 1;
for i = 1:3
  mag = dlocdoff123pm(2,i);
  ang = dlocdoff123pa(2,i);
  [rotx, roty] = complex_rot(0,refscl*mag,-ang);
  xarrow = [0 0+rotx];
  yarrow = [0 0+roty];
  if i==3
    a=annotation('arrow','color',color(i,:),'linestyle',':','linewidth',1,'HeadLength',6,...
      'HeadWidth',6);
  else
    a=annotation('arrow','color',color(i,:),'linestyle','-','linewidth',1,'HeadLength',6,...
      'HeadWidth',6);
  end
  a.Parent = ax;
  a.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
  if i==3
    label{i} = sprintf('%s-%s',strtrim(stas(2,:)),strtrim(stas(3,:)));
  else
    label{i} = sprintf('%s-%s',strtrim(stas(1,:)),strtrim(stas(i+1,:)));
  end
  text(ax,0.02,0.25-(i-1)*0.1,label{i},'color',color(i,:),'Units','normalized',...
    'HorizontalAlignment','left');
end
%reference arrow, length of 1 km/0.0625 s
[rotx, roty] = complex_rot(0,refscl,-90);
xarrow = [0.5 0.5+rotx];
yarrow = [1.5 1.5+roty];
a=annotation('arrow','color',[.5 .5 .5],'linestyle','-','linewidth',1,'HeadLength',6,...
  'HeadWidth',6);
a.Parent = ax;
a.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
text(ax,0.1,1.7,'1 km/0.0625 s','fontsize',9,'HorizontalAlignment','left');
text(ax,0.02,0.92,'a','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
% legend(ax,p,label,'Location','best');
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');

%for 4th station & 1st station
ax = f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
axis(ax, 'equal');
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,-4:1:4);
yticks(ax,-4:1:4);
label=[];
for i = 1:4
  mag = dlocdoff14pm(2,i);
  ang = dlocdoff14pa(2,i);
  [rotx, roty] = complex_rot(0,refscl*mag,-ang);
  xarrow = [0 0+rotx];
  yarrow = [0 0+roty];
  a=annotation('arrow','color',color(i,:),'linestyle','-','linewidth',1,'HeadLength',6,...
    'HeadWidth',6);
  a.Parent = ax;
  a.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
  label{i} = sprintf('%s-%s',strtrim(stas(1,:)),strtrim(stas(i+3,:)));
  text(ax,0.55,0.35-(i-1)*0.1,label{i},'color',color(i,:),'Units','normalized',...
    'HorizontalAlignment','left');
end  
%reference arrow, length of 1 km/0.0625 s
[rotx, roty] = complex_rot(0,refscl,-90);
xarrow = [0.5 0.5+rotx];
yarrow = [1.5 1.5+roty];
a=annotation('arrow','color',[.5 .5 .5],'linestyle','-','linewidth',1,'HeadLength',6,...
  'HeadWidth',6);
a.Parent = ax;
a.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
text(ax,0.1,1.7,'1 km/0.0625 s','fontsize',9,'HorizontalAlignment','left');
text(ax,0.02,0.92,'b','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
% legend(ax,p,label,'Location','best');
xlabel(ax,'E (km)');
% ylabel(ax,'N (km)');
nolabels(ax,2);

fname = strcat('mindlocdoff.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

keyboard




%% Gradient of travel time difference (doff1i/dloc)
nrow = 1;
ncol = 3;
widin = ncol*2.2;
htin = nrow*2.4;
pltxran = [0.10 0.96]; pltyran = [0.15 0.96]; % optimal axis location
pltxsep = 0.025; pltysep = 0.02;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

xran = [-5 5];
yran = [-5 5];

for i = 1:3
  ax = f.ax(i); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
  axis(ax, 'equal');
  xlim(ax,xran);
  ylim(ax,yran);
  xticks(ax,-4:2:4);
  yticks(ax,-4:2:4);
  [rotx, roty] = complex_rot(0,4,-0);
  xarrow = [0-rotx 0+rotx];
  yarrow = [0-roty 0+roty];
  p=annotation('arrow','color','k','linestyle','-','linewidth',1,'HeadLength',6,...
    'HeadWidth',6);
  p.Parent = ax;
  p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
  [rotx, roty] = complex_rot(0,4,-90);
  xarrow = [0-rotx 0+rotx];
  yarrow = [0-roty 0+roty];
  p=annotation('arrow','color','k','linestyle','-','linewidth',1,'HeadLength',6,...
    'HeadWidth',6);
  p.Parent = ax;
  p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
  [rotx, roty] = complex_rot(0,4,-45);
  xarrow = [0-rotx 0+rotx];
  yarrow = [0-roty 0+roty];
  p=annotation('arrow','color','k','linestyle','-','linewidth',1,'HeadLength',6,...
    'HeadWidth',6);
  p.Parent = ax;
  p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
  [rotx, roty] = complex_rot(0,4,-135);
  xarrow = [0-rotx 0+rotx];
  yarrow = [0-roty 0+roty];
  p=annotation('arrow','color','k','linestyle','-','linewidth',1,'HeadLength',6,...
    'HeadWidth',6);
  p.Parent = ax;
  p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
  %gradients at each src location along each direction
  scatter(ax,dx1,dy1,40,doff1231dloc(:,i),'filled');
  scatter(ax,dx2,dy2,40,doff1232dloc(:,i),'filled');
  scatter(ax,dx3,dy3,40,doff1233dloc(:,i),'filled');
  scatter(ax,dx4,dy4,40,doff1234dloc(:,i),'filled');
  colormap(ax,'viridis'); 
%   colormap(ax,flipud(colormap(ax,'kelicol')));
%   colormap(ax,'redblue'); 
%   colormap('bluewhitered'); 
%   c=colorbar(ax);
  %median of gradients along each direction
  [~,indmax] = max(abs(doff123dloc(:,i)));
  pltloc = [3.5 0.8;
            0 4.4;
            3.5 -3.5;
            3.5 3.5];
  text(ax,pltloc(indmax,1),pltloc(indmax,2),sprintf('%.3f',doff123dloc(indmax,i)),...
    'fontsize',9,'HorizontalAlignment','center','FontWeight','bold');
  ind = setdiff(1:4, indmax);
  for j=1:length(ind)
    text(ax,pltloc(ind(j),1),pltloc(ind(j),2),sprintf('%.3f',doff123dloc(ind(j),i)),...
      'fontsize',9,'HorizontalAlignment','center');
  end
  text(ax,4.9,4.4,'(s/km)','fontsize',9,'HorizontalAlignment','right');
%   text(ax,3.5,0.8,sprintf('%.3f',doff123dloc(1,i)),'fontsize',9,...
%     'HorizontalAlignment','center');
%   text(ax,0,4.4,sprintf('%.3f s/km',doff123dloc(2,i)),'fontsize',9,...
%     'HorizontalAlignment','center');
%   text(ax,3.5,-3.5,sprintf('%.3f',doff123dloc(3,i)),'fontsize',9,...
%     'HorizontalAlignment','center');
%   text(ax,3.5,3.5,sprintf('%.3f',doff123dloc(4,i)),'fontsize',9,...
%     'HorizontalAlignment','center');
  if i==3
    text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(2,:)),strtrim(stas(3,:))),...
      'Units','normalized','HorizontalAlignment','left');
%     title(ax,strcat(stas(2,:),'--',stas(3,:)));
  else
    text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(1,:)),strtrim(stas(i+1,:))),...
      'Units','normalized','HorizontalAlignment','left');
%     title(ax,strcat(stas(1,:),'--',stas(i+1,:)));
  end
%   cm = max(abs(c.Limits));
%   caxis(ax,[-cm cm]);
  if i==1
    xlabel(ax,'E (km)');
    ylabel(ax,'N (km)');
  else
    nolabels(ax,3);
  end
end
% supertit(f.ax,'Gradient of travel time difference (doff1i/dloc) w/i trio stations');
fname = strcat('doffdloctrio.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

%%
nrow = 3;
ncol = nsta-3;
widin = ncol*2.1;
htin = nrow*2.1;
pltxran = [0.06 0.96]; pltyran = [0.06 0.96]; % optimal axis location
pltxsep = 0.02; pltysep = 0.02;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

for i = 1:nrow
  for j = 1:ncol
    k = (i-1)*ncol+j;
    ax = f.ax(k); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
    axis(ax, 'equal');
    xlim(ax,xran);
    ylim(ax,yran);
    xticks(ax,-4:2:4);
    yticks(ax,-4:2:4);
    [rotx, roty] = complex_rot(0,4,-0);
    xarrow = [0-rotx 0+rotx];
    yarrow = [0-roty 0+roty];
    p=annotation('arrow','color','k','linestyle','-','linewidth',1,'HeadLength',6,...
    'HeadWidth',6);
    p.Parent = ax;
    p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
    [rotx, roty] = complex_rot(0,4,-90);
    xarrow = [0-rotx 0+rotx];
    yarrow = [0-roty 0+roty];
    p=annotation('arrow','color','k','linestyle','-','linewidth',1,'HeadLength',6,...
    'HeadWidth',6);
    p.Parent = ax;
    p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
    [rotx, roty] = complex_rot(0,4,-45);
    xarrow = [0-rotx 0+rotx];
    yarrow = [0-roty 0+roty];
    p=annotation('arrow','color','k','linestyle','-','linewidth',1,'HeadLength',6,...
    'HeadWidth',6);
    p.Parent = ax;
    p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
    [rotx, roty] = complex_rot(0,4,-135);
    xarrow = [0-rotx 0+rotx];
    yarrow = [0-roty 0+roty];
    p=annotation('arrow','color','k','linestyle','-','linewidth',1,'HeadLength',6,...
    'HeadWidth',6);
    p.Parent = ax;
    p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
    scatter(ax,dx1,dy1,40,doff141dloc(:,k),'filled');
    scatter(ax,dx2,dy2,40,doff142dloc(:,k),'filled');
    scatter(ax,dx3,dy3,40,doff143dloc(:,k),'filled');
    scatter(ax,dx4,dy4,40,doff144dloc(:,k),'filled');
    colormap(ax,'viridis');
    %   colormap(ax,flipud(colormap(ax,'kelicol')));
    %   colormap(ax,'redblue');
    %     colormap('bluewhitered');
    %     c=colorbar(ax);
    [~,indmax] = max(abs(doff14dloc(:,k)));
    pltloc = [3.5 0.8;
              0 4.4;
              3.5 -3.5;
              3.5 3.5];
    text(ax,pltloc(indmax,1),pltloc(indmax,2),sprintf('%.3f',doff14dloc(indmax,k)),...
      'fontsize',9,'HorizontalAlignment','center','FontWeight','bold');
    ind = setdiff(1:4, indmax);
    for jj=1:length(ind)
      text(ax,pltloc(ind(jj),1),pltloc(ind(jj),2),sprintf('%.3f',doff14dloc(ind(jj),k)),...
        'fontsize',9,'HorizontalAlignment','center');
    end
    text(ax,4.9,4.4,'(s/km)','fontsize',9,'HorizontalAlignment','right');
%     text(ax,3.5,0.8,sprintf('%.3f',doff14dloc(1,k)),'fontsize',9,...
%       'HorizontalAlignment','center');
%     text(ax,0,4.4,sprintf('%.3f s/km',doff14dloc(2,k)),'fontsize',9,...
%       'HorizontalAlignment','center');
%     text(ax,3.5,-3.5,sprintf('%.3f',doff14dloc(3,k)),'fontsize',9,...
%       'HorizontalAlignment','center');
%     text(ax,3.5,3.5,sprintf('%.3f',doff14dloc(4,k)),'fontsize',9,...
%       'HorizontalAlignment','center');
    text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(i,:)),strtrim(stas(j+3,:))),'Units',...
      'normalized','HorizontalAlignment','left');
%     title(ax,strcat(stas(i,:),'--',stas(j+3,:)));
%     cm = max(abs(c.Limits));
%     caxis(ax,[-cm cm]);
    if i == nrow && j==1
      xlabel(ax,'E (km)');
      ylabel(ax,'N (km)');
    else
      nolabels(ax,3);
    end
    longticks(ax,1.5);
  end
end
% supertit(f.ax,'Gradient of travel time difference (doff1i/dloc) w/i 4th and trio stas');
fname = strcat('doffdloc4thvstrio.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

%% how many km change in location would lead to a quarter of a period change in 'doff'
nrow = 1;
ncol = 3;
widin = ncol*2.2;
htin = nrow*2.4;
pltxran = [0.10 0.96]; pltyran = [0.15 0.96]; % optimal axis location
pltxsep = 0.025; pltysep = 0.02;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

xran = [-5 5];
yran = [-5 5];

for i = 1:3
  ax = f.ax(i); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
  axis(ax, 'equal');
  xlim(ax,xran);
  ylim(ax,yran);
  xticks(ax,-4:2:4);
  yticks(ax,-4:2:4);
  [rotx, roty] = complex_rot(0,4,-0);
  xarrow = [0-rotx 0+rotx];
  yarrow = [0-roty 0+roty];
  p=annotation('arrow','color','k','linestyle','-','linewidth',1,'HeadLength',6,...
    'HeadWidth',6);
  p.Parent = ax;
  p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
  [rotx, roty] = complex_rot(0,4,-90);
  xarrow = [0-rotx 0+rotx];
  yarrow = [0-roty 0+roty];
  p=annotation('arrow','color','k','linestyle','-','linewidth',1,'HeadLength',6,...
    'HeadWidth',6);
  p.Parent = ax;
  p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
  [rotx, roty] = complex_rot(0,4,-45);
  xarrow = [0-rotx 0+rotx];
  yarrow = [0-roty 0+roty];
  p=annotation('arrow','color','k','linestyle','-','linewidth',1,'HeadLength',6,...
    'HeadWidth',6);
  p.Parent = ax;
  p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)];
  [rotx, roty] = complex_rot(0,4,-135);
  xarrow = [0-rotx 0+rotx];
  yarrow = [0-roty 0+roty];
  p=annotation('arrow','color','k','linestyle','-','linewidth',1,'HeadLength',6,...
    'HeadWidth',6);
  p.Parent = ax;
  p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
  %how many km needed for a quarter period change in 'doff' at each src location along each direction
  scatter(ax,dx1,dy1,40,dlocdoff123p1d(:,i),'filled');
  scatter(ax,dx2,dy2,40,dlocdoff123p2d(:,i),'filled');
  scatter(ax,dx3,dy3,40,dlocdoff123p3d(:,i),'filled');
  scatter(ax,dx4,dy4,40,dlocdoff123p4d(:,i),'filled');
  colormap(ax,'viridis'); 
%   colormap(ax,flipud(colormap(ax,'kelicol')));
%   colormap(ax,'redblue'); 
%   colormap('bluewhitered'); 
%   c=colorbar(ax);
  %median along each direction
  [~,indmin] = min(abs(dlocdoff123p(:,i)));
  pltloc = [3.5 0.8;
            0 4.4;
            3.5 -3.5;
            3.5 3.5];
  text(ax,pltloc(indmin,1),pltloc(indmin,2),sprintf('%.1f',dlocdoff123p(indmin,i)),...
    'fontsize',9,'HorizontalAlignment','center','FontWeight','bold');
  ind = setdiff(1:4, indmin);
  for j=1:length(ind)
    text(ax,pltloc(ind(j),1),pltloc(ind(j),2),sprintf('%.1f',dlocdoff123p(ind(j),i)),...
      'fontsize',9,'HorizontalAlignment','center');
  end
  text(ax,4.9,4.4,'(km)','fontsize',9,'HorizontalAlignment','right');
%   text(ax,3.5,0.8,sprintf('%.1f',dlocdoff123p(1,i)),'fontsize',9,...
%     'HorizontalAlignment','center');
%   text(ax,0,4.4,sprintf('%.1f km',dlocdoff123p(2,i)),'fontsize',9,...
%     'HorizontalAlignment','center');
%   text(ax,3.5,-3.5,sprintf('%.1f',dlocdoff123p(3,i)),'fontsize',9,...
%     'HorizontalAlignment','center');
%   text(ax,3.5,3.5,sprintf('%.1f',dlocdoff123p(4,i)),'fontsize',9,...
%     'HorizontalAlignment','center');
  if i==3
    text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(2,:)),strtrim(stas(3,:))),...
      'Units','normalized','HorizontalAlignment','left');
%     title(ax,strcat(stas(2,:),'--',stas(3,:)));
  else
    text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(1,:)),strtrim(stas(i+1,:))),...
      'Units','normalized','HorizontalAlignment','left');
%     title(ax,strcat(stas(1,:),'--',stas(i+1,:)));
  end
%   cm = max(abs(c.Limits));
%   caxis(ax,[-cm cm]);
  if i==1
    xlabel(ax,'E (km)');
    ylabel(ax,'N (km)');
  else
    nolabels(ax,3);
  end
end
% supertit(f.ax,'Dloc needed for a 0.25/4-s change in travel time difference w/i trio stations');
fname = strcat('dlocdofftrio.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

%%
nrow = 3;
ncol = nsta-3;
widin = ncol*2.1;
htin = nrow*2.1;
pltxran = [0.06 0.96]; pltyran = [0.06 0.96]; % optimal axis location
pltxsep = 0.02; pltysep = 0.02;
f = initfig(widin,htin,nrow,ncol);
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

for i = 1:nrow
  for j = 1:ncol
    k = (i-1)*ncol+j;
    ax = f.ax(k); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
    axis(ax, 'equal');
    xlim(ax,xran);
    ylim(ax,yran);
    xticks(ax,-4:2:4);
    yticks(ax,-4:2:4);
    [rotx, roty] = complex_rot(0,4,-0);
    xarrow = [0-rotx 0+rotx];
    yarrow = [0-roty 0+roty];
    p=annotation('arrow','color','k','linestyle','-','linewidth',1,'HeadLength',6,...
    'HeadWidth',6);
    p.Parent = ax;
    p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
    [rotx, roty] = complex_rot(0,4,-90);
    xarrow = [0-rotx 0+rotx];
    yarrow = [0-roty 0+roty];
    p=annotation('arrow','color','k','linestyle','-','linewidth',1,'HeadLength',6,...
    'HeadWidth',6);
    p.Parent = ax;
    p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
    [rotx, roty] = complex_rot(0,4,-45);
    xarrow = [0-rotx 0+rotx];
    yarrow = [0-roty 0+roty];
    p=annotation('arrow','color','k','linestyle','-','linewidth',1,'HeadLength',6,...
    'HeadWidth',6);
    p.Parent = ax;
    p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
    [rotx, roty] = complex_rot(0,4,-135);
    xarrow = [0-rotx 0+rotx];
    yarrow = [0-roty 0+roty];
    p=annotation('arrow','color','k','linestyle','-','linewidth',1,'HeadLength',6,...
    'HeadWidth',6);
    p.Parent = ax;
    p.Position = [xarrow(1), yarrow(1), xarrow(2)-xarrow(1), yarrow(2)-yarrow(1)] ;
    scatter(ax,dx1,dy1,40,dlocdoff14p1d(:,k),'filled');
    scatter(ax,dx2,dy2,40,dlocdoff14p2d(:,k),'filled');
    scatter(ax,dx3,dy3,40,dlocdoff14p3d(:,k),'filled');
    scatter(ax,dx4,dy4,40,dlocdoff14p4d(:,k),'filled');
      colormap(ax,'viridis');
    %   colormap(ax,flipud(colormap(ax,'kelicol')));
    %   colormap(ax,'redblue');
%     colormap('bluewhitered');
%     c=colorbar(ax);
    [~,indmin] = min(abs(dlocdoff14p(:,k)));
    pltloc = [3.5 0.8;
              0 4.4;
              3.5 -3.5;
              3.5 3.5];
    text(ax,pltloc(indmin,1),pltloc(indmin,2),sprintf('%.1f',dlocdoff14p(indmin,k)),...
      'fontsize',9,'HorizontalAlignment','center','FontWeight','bold');
    ind = setdiff(1:4, indmin);
    for jj=1:length(ind)
      text(ax,pltloc(ind(jj),1),pltloc(ind(jj),2),sprintf('%.1f',dlocdoff14p(ind(jj),k)),...
        'fontsize',9,'HorizontalAlignment','center');
    end
    text(ax,4.9,4.4,'(km)','fontsize',9,'HorizontalAlignment','right');
%     text(ax,3.5,0.8,sprintf('%.1f',dlocdoff14p(1,k)),'fontsize',9,...
%       'HorizontalAlignment','center');
%     text(ax,0,4.4,sprintf('%.1f km',dlocdoff14p(2,k)),'fontsize',9,...
%       'HorizontalAlignment','center');
%     text(ax,3.5,-3.5,sprintf('%.1f',dlocdoff14p(3,k)),'fontsize',9,...
%       'HorizontalAlignment','center');
%     text(ax,3.5,3.5,sprintf('%.1f',dlocdoff14p(4,k)),'fontsize',9,...
%       'HorizontalAlignment','center');
    text(ax,0.02,0.05,sprintf('%s-%s',strtrim(stas(i,:)),strtrim(stas(j+3,:))),'Units',...
      'normalized','HorizontalAlignment','left');
%     title(ax,strcat(stas(i,:),'--',stas(j+3,:)));
%     cm = max(abs(c.Limits));
%     caxis(ax,[-cm cm]);
    if i == nrow && j==1
      xlabel(ax,'E (km)');
      ylabel(ax,'N (km)');
    else
      nolabels(ax,3);
    end
    longticks(ax,1.5);
  end
end
% supertit(f.ax,'Dloc needed for a 0.25/4-s change in travel time difference w/i 4th and trio stas');
fname = strcat('dlocdoff4thvstrio.pdf');
print(f.fig,'-dpdf',...
  strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
 





