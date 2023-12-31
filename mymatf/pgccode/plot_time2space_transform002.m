% function plot_time2space_transform002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to plot the location of points from a square of 9 points
% around the center of LFE fam 002, to see how the time transforms into the space
% coordinates.
%
%
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2020/08/19
% Last modified date:   2020/08/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all
clc

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/forsummary');

offcont = load(fullfile(rstpath,'offset_transfromtest'));
loccont = load(fullfile(rstpath,'eventloc.002.transform'));

%%% convert lfe location to km relative to 043
relacont = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),loccont(end,1),loccont(end,2));
relacont(:,1) = dx;
relacont(:,2) = dy;

f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 4;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);

nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

ax = f.ax(1);
hold(ax,'on');

sps = 40;
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
plot(ax,offcont(:,1)*sps,offcont(:,2)*sps,'b-','linew',1);
scatter(ax,offcont(:,1)*sps,offcont(:,2)*sps, 20, [0.5 0.5 0.5], 'filled','o');
strlbl = ['A';'B';'C';'D';'E';'F';'G';'H';'I']; 
text(ax,offcont(:,1)*sps,offcont(:,2)*sps+2,strlbl(:,:),'fontsize',12);

ax.Box='on';
grid(ax,'on');
ax.GridLineStyle = '--';
axis(ax,'equal');
axis(ax,[-30 30 -30 30]);
xlabel(ax,'off12 (sample) 40 sps');
ylabel(ax,'off13 (sample) 40 sps');
hold(ax,'off');


ax = f.ax(2);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
plot(ax,relacont(:,1),relacont(:,2),'b-','linew',1);
scatter(ax,relacont(:,1),relacont(:,2), 20, [0.5 0.5 0.5], 'filled','o');
strlbl = ['A';'B';'C';'D';'E';'F';'G';'H';'I']; 
text(ax,relacont(:,1),relacont(:,2)+1,strlbl(:,:),'fontsize',12);

ax.Box='on';
grid(ax,'on');
ax.GridLineStyle = '--';
axis(ax,'equal');
axis(ax,[-20 20 -20 20]);
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
hold(ax,'off');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/time2space_transform002.pdf'));


















