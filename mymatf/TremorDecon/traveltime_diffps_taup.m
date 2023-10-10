% function traveltime_diffps_taup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to use the customized velocity model Puget Sound P3 
% velocity model created by 'taupcreate' to compute the theoretical travel
% time difference at the Cascadia stations, given a set of source locations
% between phases P and S. We would like to know, for example, if srcs are
% 1 km apart, what is their arrival time / travel time difference, OR, to'
% make the P-S trave time difference to be different by half a cycle at one 
% station, so that the constructive stacking at S-wave timings would be 
% deconstructive for P-wave, how apart the srcs need to be?
% The idea is to get a map of S-P travel time difference at all stations 
% from a grid of locations (which come from the inversion using hypoinverse
% of a grid of differential travel times between the trio stations, off12 
% and off13). 
% 
% 
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2023/07/26
% Last modified date:   2023/07/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clc;
clear;
close all;

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, resol] = pixelperinch(1);

%% generate a grid of offsets around fam 002, then transform the offset grid to relative locations
sps = 160;
offmax = round(4*sps/40); % make it compatible to different sps
%grid of time offsets
i = 1;
for off12 = -offmax: 1: offmax
  for off13 = -offmax: 1: offmax
    offset(i,1:2) = [off12 off13];  % offset in samples at requested sampling rate.
    i = i+1;
  end
end

%transform to locations
ftrans = 'interpchao';
[evtloc, indinput] = off2space002(offset,sps,ftrans,0); 
% loc has 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
nevt = size(evtloc,1);

figure
subplot(1,2,1)
axis equal
hold on; grid on; box on
plot(evtloc(:,7),evtloc(:,8),'.');
xlabel(sprintf('off12 samples (%d Hz)',sps));
ylabel(sprintf('off13 samples (%d Hz)',sps));
tmp = ceil(max(max(abs([evtloc(:,7),evtloc(:,8)]))))+1;
xran = [-tmp tmp];
yran = [-tmp tmp];
xlim(xran);
ylim(yran);
subplot(1,2,2)
axis equal
hold on; grid on; box on
plot(evtloc(:,1),evtloc(:,2),'.');
text(0.95,0.05,sprintf('%d unique sources',size(evtloc,1)),'Units','normalized',...
  'HorizontalAlignment','right');
tmp = ceil(max(max(abs([evtloc(:,1),evtloc(:,2)]))));
xran = [-tmp tmp];
yran = [-tmp tmp];
xlim(xran);
ylim(yran);
xlabel('E (km)');
ylabel('N (km)');


%% compute the travel time difference using Taup at all stations
stas=['PGC  '
      'SSIB '
      'SILB '
      'LZB  '
      'TWKB '
      'MGCB '
      'KLNB '
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
ttp = zeros(nevt, nsta);  % S-wave travel time
tts = zeros(nevt, nsta);  % S-wave travel time
tic

%%%
[loc0, indinput] = off2space002([0,0],sps,ftrans,0); 
for i = 1: size(loc0,1)
  for j = 1:3
    tt=tauptime('mod',pugetp3,'dep',loc0(i,5),'ph','s,S','evt',[loc0(i,4) loc0(i,3)],...
        'sta',[staloc(j,2) staloc(j,1)]); 
  end
end

for i = 1: nevt
  for j = 1: nsta
    tt=tauptime('mod',pugetp3,'dep',evtloc(i,5),'ph','p,P','evt',[evtloc(i,4) evtloc(i,3)],...
      'sta',[staloc(j,2) staloc(j,1)]); 
    if ~strcmp(tt(1).phase,'p')
      fprintf('At i=%d, j=%d, p take-off angle might be >90 degrees, ie., downgoing wave \n.',...
        i,j);
      ttp(i,j) = 0;   % choose the first P arrival
    else
      ttp(i,j) = tt(1).time;   % choose the first P arrival
    end
    tt=tauptime('mod',pugetp3,'dep',evtloc(i,5),'ph','s,S','evt',[evtloc(i,4) evtloc(i,3)],...
      'sta',[staloc(j,2) staloc(j,1)]);
    if ~strcmp(tt(1).phase,'s')
      fprintf('At i=%d, j=%d, s take-off angle might be >90 degrees, ie., downgoing wave \n.',...
        i,j);
      tts(i,j) = 0;   % choose the first S arrival
    else
      tts(i,j) = tt(1).time;   % choose the first P arrival
    end
  end
end
toc
          
dttps = tts-ttp;  % differential traveltime realative to main sta PGC, in sec
dttpsspl = dttps*sps; % differential traveltime in samples from taup

%%
f = initfig(13,9,3,3); %initialize fig
xran = [0.06 0.96]; yran = [0.06 0.98];
xsep = 0.05; ysep = 0.02;
optaxpos(f,3,3,xran,yran,xsep,ysep);
for ista = 1: nsta
  ax=f.ax(ista);
  hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
  scatter(ax,evtloc(:,7),evtloc(:,8),80,dttpsspl(:,ista),'s','filled','MarkerEdgeColor','none');
  tmp = ceil(max(max(abs([evtloc(:,7),evtloc(:,8)]))))+1;
  xran = [-tmp tmp];
  yran = [-tmp tmp];
  axis(ax,[xran yran]);
  axis(ax,'equal');
  oldc = colormap(ax,'kelicol');
  newc = flipud(oldc);
  colormap(ax,newc);
  c=colorbar(ax);
  title(ax,stas(ista,:))
  if ista==1
    xlabel(ax,sprintf('off12 samples (%d Hz)',sps));
    ylabel(ax,sprintf('off13 samples (%d Hz)',sps));
    c.Label.String = sprintf('S-P trvl time samples (%d Hz)',sps);
  end
end
delete(f.ax(8:9));

%%
f = initfig(13,9,3,3); %initialize fig
xran = [0.06 0.96]; yran = [0.06 0.98];
xsep = 0.05; ysep = 0.02;
optaxpos(f,3,3,xran,yran,xsep,ysep);
for ista = 1: nsta
  ax=f.ax(ista);
  hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
  scatter(ax,evtloc(:,1),evtloc(:,2),40,dttpsspl(:,ista),'o','filled','MarkerEdgeColor','none');
  tmp = ceil(max(max(abs([evtloc(:,1),evtloc(:,2)]))));
  xran = [-tmp tmp];
  yran = [-tmp tmp];
  axis(ax,[xran yran],'equal');
  oldc = colormap(ax,'kelicol');
  newc = flipud(oldc);
  colormap(ax,newc);
  c=colorbar(ax);
  title(ax,stas(ista,:))
  if ista==1
    xlabel(ax,'E (km)');
    ylabel(ax,'N (km)');
    c.Label.String = sprintf('S-P trvl time samples (%d Hz)',sps);
  end
end
delete(f.ax(8:9));



