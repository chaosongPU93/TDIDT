% testgradtrvlpgc.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the traveltime from a grid of sources to PGC, so to estimate the 
% gradient of traveltime wrt to location in E-W and N-S direction.
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2023/09/12
% Last modified date:   2023/09/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
format short e   % Set the format to 5-digit floating point
clear
clc
close all

set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

sps = 40;
ftrans = 'interpchao';
[loc, indinput] = off2space002([],sps,ftrans,0);
% loc has 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13

f = initfig(10,5,1,2); %initialize fig
ax = f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
scatter(ax,loc(:,7),loc(:,8),20,loc(:,6),'filled');
oldc = colormap(ax,'kelicol');
newc = flipud(oldc);
colormap(ax,newc);
c=colorbar(ax);
c.Label.String = 'traveltime to PGC (s)';
axis(ax,'equal');
axis(ax,[-30 30 -30 30]);
xlabel(ax,sprintf('PGC-SSIB offset in samples at %d sps', sps),'fontsize',11);
ylabel(ax,sprintf('PGC-SILB offset in samples at %d sps', sps),'fontsize',11);

ax = f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
scatter(ax,loc(:,1),loc(:,2),20,loc(:,6),'filled');
oldc = colormap(ax,'kelicol');
newc = flipud(oldc);
colormap(ax,newc);
c=colorbar(ax);
c.Label.String = 'traveltime to PGC (s)';
axis(ax,'equal');
axis(ax,[-30 20 -30 20]);
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);



