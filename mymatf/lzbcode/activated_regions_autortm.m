% function activated_regions_autortm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is to use the density plot to lighten up the activated regions
% during the ETS periods in hf and lf
% Use the new automatic time ranges of migrations that come from 'identify_RTMs.m'
% 
%
%    2020/09/10, i am adding a fam 001 back again for last try, 13 fam in total now
%                this is for adding new fams 006, 001
%
% Chao Song, chaosong@princeton.edu
% First created date:   2020/12/29
% Last modified date:   2020/12/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all
clc

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/lzbrst');

nfampool = ['002';
            '043';
            '141';
            '047';
            '010';
            '144';
            '099';
            '068';
            '125';
            '147';
            '017';
            '006';
            '001';
%             '158';      % 158, 20200916,testing purpose
            ];     % 2020/09/07, i am adding a new family 006 to the pool, final version
        
nfam = size(nfampool,1);
disp(nfam); 

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
         'LZB'];
POLSTA=['SSIB '           % polaris station names
        'SILB '
        'KLNB '
        'MGCB '
        'TWKB '];
    
stas=['TWKB '
      'LZB  '
      'MGCB '];     % determine the trio and order         

  
% load files
winlenhf = 4;

winlenlf = 16;
lofflf = 4;
ccminlf = 0.35;

SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.dcutnodou.',SUFFIXhf);
hftime = load(fname);
% 25 cols, format is:
%   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date
%   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum
%


SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.dcutnodou.',SUFFIXlf);
lftime = load(fname);

% this is inverted from (0,0) of all fams, same order, location of control points
loccont = [-123.492667 48.451500 38.1400; 
           -123.772167 48.493000 35.5900; 
           -123.863167 48.528167 35.2100;
           -123.603333 48.440167 36.7100;
           -123.800167 48.408833 34.5200;
           -123.893333 48.536500 35.0700;
           -123.864500 48.498667 34.8800;
           -123.753333 48.525667 36.2000;
           -123.703667 48.502667 36.4100;
           -123.814333 48.538667 35.7900;
           -123.838500 48.544833 35.6600;
           -123.908000 48.494167 34.5100;       % 006
           -123.879667 48.446167 34.2600;       % 001
%            -123.967000 48.458667 33.7800;       % 158, 20200916,testing purpose
           ];
       

relacont = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),loccont(2,1),loccont(2,2));
relacont(:,1) = dx;
relacont(:,2) = dy;


% sort according to day, sec
hftime = sortrows(hftime, [13, 15]);
lftime = sortrows(lftime, [13, 15]);

xran = [-20 25];
yran = [-20 20];

%%

f2.fig=figure;
f2.fig.Renderer='Painters';
widin = 12;  % maximum width allowed is 8.5 inches
htin = 10;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 3;
ncol = 2;
for isub = 1:nrow*ncol
    f2.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f2.ax(1), 'position', [ 0.1, 0.71, 0.28, 0.28]);
set(f2.ax(2), 'position', [ 0.5, 0.71, 0.28, 0.28]);
set(f2.ax(3), 'position', [ 0.1, 0.38, 0.28, 0.28]);
set(f2.ax(4), 'position', [ 0.5, 0.38, 0.28, 0.28]);
set(f2.ax(5), 'position', [ 0.1, 0.05, 0.28, 0.28]);
set(f2.ax(6), 'position', [ 0.5, 0.05, 0.28, 0.28]);

% xran = [-14 9];
% yran = [-9 13];

msizehf = 4;
msizelf = 10;

% subplot 1 of figure i
hfplt = hftime(hftime(:,13) < 2004*1000, :);
lfplt = lftime(lftime(:,13) < 2004*1000, :);

ax = f2.ax(1);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
%     ind = find(density1d>=0);
%     scatter(ax,xloc(ind),yloc(ind), 4, log(density1d(ind)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
%     aaa = density2d';
%     bbb = aaa(aaa>=1);
%     ccc = log(bbb);
%     imagesc(ax,xran+dx, yran+dy, ccc);
% create a density matrix to store the number of detections in each small grid
dxhf = 0.2;
dyhf = 0.2;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'EastOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, '2003','unit','normalized');
text(ax, 0.7, 0.8, 'HF','unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');

% subplot 2 of figure i
ax = f2.ax(2);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxlf = 0.5;
dylf = 0.5;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
    xran,yran,dxlf,dylf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'EastOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, '2003','unit','normalized');
text(ax, 0.7, 0.8, 'LF','unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');


% subplot 3 of figure i
hfplt = hftime(hftime(:,13) > 2004*1000 & hftime(:,13) < 2005*1000, :);
lfplt = lftime(lftime(:,13) > 2004*1000 & lftime(:,13) < 2005*1000, :);

ax = f2.ax(3);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxhf = 0.2;
dyhf = 0.2;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'EastOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, '2004','unit','normalized');
text(ax, 0.7, 0.8, 'HF','unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');

% subplot 4 of figure i
ax = f2.ax(4);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxlf = 0.5;
dylf = 0.5;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
    xran,yran,dxlf,dylf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'EastOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, '2004','unit','normalized');
text(ax, 0.7, 0.8, 'LF','unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');


% subplot 5 of figure i
hfplt = hftime(hftime(:,13) > 2005*1000, :);
lfplt = lftime(lftime(:,13) > 2005*1000, :);

ax = f2.ax(5);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxhf = 0.2;
dyhf = 0.2;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'EastOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, '2005','unit','normalized');
text(ax, 0.7, 0.8, 'HF','unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');

% subplot 6 of figure i
ax = f2.ax(6);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;

% create a density matrix to store the number of detections in each small grid
dxlf = 0.5;
dylf = 0.5;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
    xran,yran,dxlf,dylf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'EastOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, '2005','unit','normalized');
text(ax, 0.7, 0.8, 'LF','unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');

print(f2.fig,'-depsc2',strcat(rstpath,'/density_eachETS.eps'));


%%
trangeauto = [
    2004195   4.1784e+04   4.2666e+04 % 2004195,4.15e+4,4.60e+4;
    2004197   3.3831e+04   3.5500e+04
    2004197   4.0944e+04   4.2168e+04
    2004197   5.4704e+04   5.5897e+04
    2004197   5.6240e+04   5.7257e+04
    2004197   6.9920e+04   7.1142e+04
    2004197   7.7484e+04   7.9841e+04
    2004197   7.9919e+04   8.1483e+04
    2004197   8.3766e+04   8.6638e+04
    2004198   5.9230e+03   9.7670e+03
    2004198   1.9998e+04   2.1646e+04
    2004198   4.2602e+04   4.3751e+04
    2004198   5.3476e+04   5.8080e+04
    2004198   5.9624e+04   6.0467e+04
    2004198   6.2936e+04   6.3915e+04
    2004198   7.3510e+04   7.6059e+04
    2004198   8.4192e+04   8.5952e+04
    2004199   4.1370e+03   6.1680e+03
    2004199   4.2815e+04   4.3333e+04
    2004199   4.7777e+04   4.9156e+04
    2004199   4.9714e+04   5.0644e+04
    2004199   8.0831e+04   8.3462e+04
    2004200   1.1676e+04   1.5829e+04
    2004200   1.5884e+04   1.7155e+04
    2004200   1.9075e+04   2.0140e+04
    2004200   4.8079e+04   4.8696e+04
    2004201   8.3400e+02   1.6590e+03
    2004203   1.8492e+04   2.0806e+04
    2005255   4.6069e+04   4.6979e+04
    2005255   5.8020e+04   5.9341e+04
    2005255   6.0351e+04   6.1587e+04
    2005255   6.7238e+04   6.8565e+04
    2005256   8.0030e+03   1.0384e+04
    2005256   1.3382e+04   1.4376e+04
    2005256   2.0407e+04   2.2324e+04
    2005256   2.7689e+04   2.9202e+04
    2005256   2.9685e+04   3.1233e+04
    2005256   3.4786e+04   3.5268e+04
    2005256   5.5237e+04   5.7052e+04
    2005256   6.0693e+04   6.4846e+04
    2005256   7.1030e+04   7.3944e+04
    2005256   7.6351e+04   7.6924e+04
    2005256   7.6984e+04   7.8694e+04
    2005257   1.5943e+04   1.8685e+04
    2005257   2.1040e+04   2.5395e+04
    2005257   2.5465e+04   2.8085e+04
    2005257   3.5320e+04   3.7628e+04
    2005257   3.8483e+04   4.5322e+04
    2005257   6.1419e+04   6.3758e+04
    2005257   7.3563e+04   8.2073e+04
    2005257   8.2129e+04   8.3410e+04
    2005258   2.4787e+04   2.6247e+04
    2005258   3.4970e+04   3.9993e+04
    2005259   1.9280e+03   4.0760e+03
    2005259   5.2020e+03   7.7440e+03
    2005259   3.7267e+04   4.1249e+04
    2005259   7.4225e+04   7.4753e+04
    2005260   2.6210e+03   3.4680e+03
    2005260   3.5740e+03   6.7080e+03
    2005260   9.3790e+03   1.0387e+04
    ];

angbestauto = [
   260
   155
   165
   220
   185
   195
   195
   195
   260
   235
   205
   230
    80
    85
    65
   235
   250
   310
   250
   135
   195
   140
   170
    65
   215
   240
   260
   120
   345
   165
   110
   155
   215
    80
   170
   140
   125
   260
   240
   175
   185
   260
   195
   170
   250
   165
   335
   165
   245
   230
    10
   250
   130
   140
   210
   115
   245
   120
   160
    20   
   ];


%% 
% this is new time ranges that come from 'identify_RTMs.m'
trangevis = [
    2004195   4.1784e+04   4.2666e+04 % 2004195,4.15e+4,4.60e+4;
    2004197   3.3831e+04   3.5500e+04
    2004197   4.0944e+04   4.2168e+04
    2004197   5.4704e+04   5.5897e+04
    2004197   5.6240e+04   5.7257e+04
    2004197   6.9920e+04   7.1142e+04
    2004197   7.7484e+04   7.9841e+04
    2004197   7.9919e+04   8.1483e+04
    2004197   8.3766e+04   8.4888e+04  % old 6, new 34, ACCEPTED
    2004197   8.5320e+04   8.6638e+04  % old 8, new 34, ACCEPTED
    2004198   5.9230e+03   9.7670e+03
    2004198   1.9998e+04   2.1646e+04
    2004198   4.2602e+04   4.3751e+04
    2004198   5.3476e+04   5.8080e+04
    2004198   6.2936e+04   6.3915e+04
    2004198   7.3510e+04   7.6059e+04
    2004198   8.4192e+04   8.5952e+04
    2004199   4.1370e+03   6.1680e+03
    2004199   4.2815e+04   4.3333e+04
    2004199   4.7777e+04   4.9156e+04
    2004199   4.9714e+04   5.0644e+04
    2004199   8.0831e+04   8.3462e+04
    2004200   1.1676e+04   1.5829e+04
    2004200   1.5884e+04   1.7155e+04
    2004200   1.9075e+04   2.0140e+04
    2004200   4.8079e+04   4.8696e+04
    2004201   8.3400e+02   1.6590e+03
    2004203   1.8492e+04   2.0806e+04
    2005255   4.6069e+04   4.6979e+04
    2005255   5.8020e+04   5.9341e+04
    2005255   6.0351e+04   6.1587e+04
    2005255   6.7238e+04   6.8565e+04   % new 10, old 22, questionable direction, 155-->65
    2005256   3590         4987         % new 16, old 25, questionable direction, 175-->245
    2005256   8.0030e+03   9.8000e+03   % new 17, old 26, change ending time as old, ACCEPTED
    2005256   1.3382e+04   1.4376e+04
    2005256   2.0407e+04   2.2324e+04
    2005256   2.8260e+04   2.9202e+04   % new 24, old 28, change start time as 7.85, ACCEPTED
    2005256   2.9685e+04   3.1233e+04
    2005256   3.4786e+04   3.5268e+04
    2005256   5.2020e+04   5.3208e+04   % new 32, old 30, try a new range, ACCEPTED
    2005256   5.5548e+04   5.7052e+04   % new 33, change start time as 15.43, ACCEPTED
    2005256   6.0693e+04   6.4846e+04
    2005256   7.1496e+04   7.3944e+04   % new 35, old 31, change start time as 19.86, ACCEPTED
    2005256   7.6351e+04   7.6924e+04
    2005256   7.6984e+04   7.8694e+04
    2005257   6120         7776         % new 41, old 33, change start time as 1.7, ACCEPTED 
    2005257   1.5943e+04   1.8685e+04
    2005257   2.1040e+04   2.5395e+04
    2005257   2.5465e+04   2.8085e+04
    2005257   3.5320e+04   3.7628e+04
    2005257   3.8483e+04   4.5322e+04
    2005257   6.1419e+04   6.3758e+04
    2005257   7.3563e+04   8.2073e+04
    2005257   8.2129e+04   8.3410e+04
    2005258   2.4787e+04   2.6247e+04
    2005258   3.4970e+04   3.9993e+04
    2005259   1.9280e+03   4.0760e+03
    2005259   5.2020e+03   7.7440e+03
    2005259   3.7267e+04   4.1249e+04
    2005259   7.4225e+04   7.4753e+04
    2005260   2.6210e+03   3.4680e+03
    2005260   3.5740e+03   6.7080e+03
    2005260   9.3790e+03   1.0387e+04
];

angbestvis = [
   260
   155
   165
   220
   185
   195
   195
   195
   230
   270
   235
   205
   230
    80
    65
   235
   250
   310
   250
   135
   195
   140
   170
    65
   215
   240
   260
   120
   345
   165
   110
   155  %65
   175  %245
   215
    80
   170
   105
   125
   260
   175
   255
   175
   180
   260
   195
   195
   170
   250
   165
   335
   165
   245
   230
    10
   250
   130
   140
   210
   115
   245
   120
   160
    20      
   ];


angbestvis(32) = 65; %67.5 112.5
angbestvis(33) = 245; %247.5 270


%%
% all detections in RTMs, lf and hf, that are located at main lzb region in each propagation group
mighfa = [];
miglfa = [];

trange = trangeauto;
angbest = angbestauto;

for i = 1: length(trange)

%     disp(trange(i,:));
    indhf = find(hftime(:,13)==trange(i,1) & hftime(:,15)>=trange(i,2) & ...
                 hftime(:,15)<=trange(i,3));
    mighf = hftime(indhf,:);
    
    mighfdum = mighf;
    for j = 1: size(mighf,1)
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbestvis(i)-90),0,0);
        mighfdum(j,1) = newx;
        mighfdum(j,2) = newy;
    end
    
    indlf = find(lftime(:,13)==trange(i,1) & lftime(:,15)>=trange(i,2) & ...
                 lftime(:,15)<=trange(i,3));
    miglf = lftime(indlf,:); 
    
    miglfdum = miglf;
    for j = 1: size(miglf,1)
        x0 = miglfdum(j,1);
        y0 = miglfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbestvis(i)-90),0,0);
        miglfdum(j,1) = newx;
        miglfdum(j,2) = newy;
    end
    
    
    mintime = max(min(mighfdum(:,15)), min(miglfdum(:,15)));
    maxtime = min(max(mighfdum(:,15)), max(miglfdum(:,15)));
    mighfdum2 = mighfdum(mighfdum(:,15)>=mintime & mighfdum(:,15)<=maxtime, :);
    miglfdum2 = miglfdum(miglfdum(:,15)>=mintime & miglfdum(:,15)<=maxtime, :);

        
    % rotate back the projected locations to get the original coordinates
    tmp = mighfdum2;
    [newx,newy] = coordinate_rot(tmp(:,1),tmp(:,2),(angbestvis(i)-90),0,0);
    tmp(:,1) = newx;
    tmp(:,2) = newy;
    mighfa = [mighfa; tmp];
    tmp = miglfdum2;
    [newx,newy] = coordinate_rot(tmp(:,1),tmp(:,2),(angbestvis(i)-90),0,0);
    tmp(:,1) = newx;
    tmp(:,2) = newy;
    miglfa = [miglfa; tmp];
   
    
end

mighfavis = [];
miglfavis = [];

trange = trangevis;
angbest = angbestvis;

for i = 1: length(trange)

%     disp(trange(i,:));
    indhf = find(hftime(:,13)==trange(i,1) & hftime(:,15)>=trange(i,2) & ...
                 hftime(:,15)<=trange(i,3));
    mighf = hftime(indhf,:);
    
    mighfdum = mighf;
    for j = 1: size(mighf,1)
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbestvis(i)-90),0,0);
        mighfdum(j,1) = newx;
        mighfdum(j,2) = newy;
    end
    
    indlf = find(lftime(:,13)==trange(i,1) & lftime(:,15)>=trange(i,2) & ...
                 lftime(:,15)<=trange(i,3));
    miglf = lftime(indlf,:); 
    
    miglfdum = miglf;
    for j = 1: size(miglf,1)
        x0 = miglfdum(j,1);
        y0 = miglfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbestvis(i)-90),0,0);
        miglfdum(j,1) = newx;
        miglfdum(j,2) = newy;
    end
    
    mintime = max(min(mighfdum(:,15)), min(miglfdum(:,15)));
    maxtime = min(max(mighfdum(:,15)), max(miglfdum(:,15)));
    mighfdum2 = mighfdum(mighfdum(:,15)>=mintime & mighfdum(:,15)<=maxtime, :);
    miglfdum2 = miglfdum(miglfdum(:,15)>=mintime & miglfdum(:,15)<=maxtime, :);
        
    % rotate back the projected locations to get the original coordinates
    tmp = mighfdum2;
    [newx,newy] = coordinate_rot(tmp(:,1),tmp(:,2),(angbestvis(i)-90),0,0);
    tmp(:,1) = newx;
    tmp(:,2) = newy;
    mighfavis = [mighfavis; tmp];
    tmp = miglfdum2;
    [newx,newy] = coordinate_rot(tmp(:,1),tmp(:,2),(angbestvis(i)-90),0,0);
    tmp(:,1) = newx;
    tmp(:,2) = newy;
    miglfavis = [miglfavis; tmp];  
    
end

%%

f2.fig=figure;
f2.fig.Renderer='Painters';
widin = 12;  % maximum width allowed is 8.5 inches
htin = 10;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 3;
ncol = 2;
for isub = 1:nrow*ncol
    f2.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f2.ax(1), 'position', [ 0.1, 0.71, 0.28, 0.28]);
set(f2.ax(2), 'position', [ 0.5, 0.71, 0.28, 0.28]);
set(f2.ax(3), 'position', [ 0.1, 0.38, 0.28, 0.28]);
set(f2.ax(4), 'position', [ 0.5, 0.38, 0.28, 0.28]);
set(f2.ax(5), 'position', [ 0.1, 0.05, 0.28, 0.28]);
set(f2.ax(6), 'position', [ 0.5, 0.05, 0.28, 0.28]);

% xran = [-14 9];
% yran = [-9 13];

msizehf = 4;
msizelf = 10;

% subplot 1 of figure i
hfplt = hftime(hftime(:,13) > 2004*1000, :);
lfplt = lftime(lftime(:,13) > 2004*1000, :);

ax = f2.ax(1);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
%     ind = find(density1d>=0);
%     scatter(ax,xloc(ind),yloc(ind), 4, log(density1d(ind)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
%     aaa = density2d';
%     bbb = aaa(aaa>=1);
%     ccc = log(bbb);
%     imagesc(ax,xran+dx, yran+dy, ccc);
% create a density matrix to store the number of detections in each small grid
dxhf = 0.2;
dyhf = 0.2;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'EastOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, '2004+2005','unit','normalized');
text(ax, 0.7, 0.8, 'HF','unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');

% subplot 2 of figure i
ax = f2.ax(2);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxlf = 0.5;
dylf = 0.5;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
    xran,yran,dxlf,dylf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'EastOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, '2004+2005','unit','normalized');
text(ax, 0.7, 0.8, 'LF','unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');


% subplot 3 of figure i
hfplt = mighfa;
lfplt = miglfa;

ax = f2.ax(3);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxhf = 0.2;
dyhf = 0.2;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'EastOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, 'auto','unit','normalized');
text(ax, 0.7, 0.8, 'HF','unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');

% subplot 4 of figure i
ax = f2.ax(4);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxlf = 0.5;
dylf = 0.5;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
    xran,yran,dxlf,dylf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'EastOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, 'auto','unit','normalized');
text(ax, 0.7, 0.8, 'LF','unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');


% subplot 5 of figure i
hfplt = mighfavis;
lfplt = miglfavis;

ax = f2.ax(5);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxhf = 0.2;
dyhf = 0.2;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'EastOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, 'auto+visual','unit','normalized');
text(ax, 0.7, 0.8, 'HF','unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');

% subplot 6 of figure i
ax = f2.ax(6);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;

% create a density matrix to store the number of detections in each small grid
dxlf = 0.5;
dylf = 0.5;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
    xran,yran,dxlf,dylf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'EastOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, 'auto+visual','unit','normalized');
text(ax, 0.7, 0.8, 'LF','unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');

print(f2.fig,'-depsc2',strcat(rstpath,'/density_rtm.eps'));


%%
trange = trangevis;
angbest = angbestvis;

mighfne = [];
miglfne = [];
mighfse = [];
miglfse = [];
mighfsw = [];
miglfsw = [];
mighfnw = [];
miglfnw = [];

for i = 1: length(trange)

%     disp(trange(i,:));
    indhf = find(hftime(:,13)==trange(i,1) & hftime(:,15)>=trange(i,2) & ...
                 hftime(:,15)<=trange(i,3));
    mighf = hftime(indhf,:);
    
    mighfdum = mighf;
    for j = 1: size(mighf,1)
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbestvis(i)-90),0,0);
        mighfdum(j,1) = newx;
        mighfdum(j,2) = newy;
    end
    
    indlf = find(lftime(:,13)==trange(i,1) & lftime(:,15)>=trange(i,2) & ...
                 lftime(:,15)<=trange(i,3));
    miglf = lftime(indlf,:); 
    
    miglfdum = miglf;
    for j = 1: size(miglf,1)
        x0 = miglfdum(j,1);
        y0 = miglfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbestvis(i)-90),0,0);
        miglfdum(j,1) = newx;
        miglfdum(j,2) = newy;
    end
    
    mintime = max(min(mighfdum(:,15)), min(miglfdum(:,15)));
    maxtime = min(max(mighfdum(:,15)), max(miglfdum(:,15)));
    mighfdum2 = mighfdum(mighfdum(:,15)>=mintime & mighfdum(:,15)<=maxtime, :);
    miglfdum2 = miglfdum(miglfdum(:,15)>=mintime & miglfdum(:,15)<=maxtime, :);
    
    % rotate back the projected locations to get the original coordinates
    tmphf = mighfdum2;
    [newx,newy] = coordinate_rot(tmphf(:,1),tmphf(:,2),(angbestvis(i)-90),0,0);
    tmphf(:,1) = newx;
    tmphf(:,2) = newy;
    tmplf = miglfdum2;
    [newx,newy] = coordinate_rot(tmplf(:,1),tmplf(:,2),(angbestvis(i)-90),0,0);
    tmplf(:,1) = newx;
    tmplf(:,2) = newy;
        
    if angbestvis(i)>=0 && angbestvis(i)<90
        mighfne = [mighfne; tmphf];
        miglfne = [miglfne; tmplf];    
    elseif angbestvis(i)>=90 && angbestvis(i)<180
        mighfse = [mighfse; tmphf];
        miglfse = [miglfse; tmplf];
    elseif angbestvis(i)>=180 && angbestvis(i)<270
        mighfsw = [mighfsw; tmphf];
        miglfsw = [miglfsw; tmplf];
    else
        mighfnw = [mighfnw; tmphf];
        miglfnw = [miglfnw; tmplf];
    end
end

%% RTM density, NE vs. SW
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.1, 0.55, 0.4, 0.4]);
set(f.ax(2), 'position', [ 0.55, 0.55, 0.4, 0.4]);
set(f.ax(3), 'position', [ 0.1, 0.1, 0.4, 0.4]);
set(f.ax(4), 'position', [ 0.55, 0.1, 0.4, 0.4]);

msizehf = 4;
msizelf = 10;

% subplot 1 of figure i
hfplt = mighfne;
lfplt = miglfne;

ax = f.ax(1);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
%     ind = find(density1d>=0);
%     scatter(ax,xloc(ind),yloc(ind), 4, log(density1d(ind)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
%     aaa = density2d';
%     bbb = aaa(aaa>=1);
%     ccc = log(bbb);
%     imagesc(ax,xran+dx, yran+dy, ccc);
% create a density matrix to store the number of detections in each small grid
dxhf = 0.2;
dyhf = 0.2;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, 'NE','unit','normalized');
text(ax, 0.7, 0.8, 'HF','unit','normalized');
text(ax, 0.7, 0.7, strcat(num2str(size(hfplt,1)),{' detections'}),'unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');

% subplot 2 of figure i
ax = f.ax(2);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxlf = 0.5;
dylf = 0.5;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
    xran,yran,dxlf,dylf);
dumlf = density1d(density1d(:,3)>0, :);
dumlf(dumlf(:,3)>1, :) = [];
scatter(ax,dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumlf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumlf(dumlf(:,3)==1, :) = [];
scatter(ax,dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, 'NE','unit','normalized');
text(ax, 0.7, 0.8, 'LF','unit','normalized');
text(ax, 0.7, 0.7, strcat(num2str(size(lfplt,1)),{' detections'}),'unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');

% subplot 3 of figure i
hfplt = mighfsw;
lfplt = miglfsw;

ax = f.ax(3);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
%     ind = find(density1d>=0);
%     scatter(ax,xloc(ind),yloc(ind), 4, log(density1d(ind)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
%     aaa = density2d';
%     bbb = aaa(aaa>=1);
%     ccc = log(bbb);
%     imagesc(ax,xran+dx, yran+dy, ccc);
% create a density matrix to store the number of detections in each small grid
dxhf = 0.2;
dyhf = 0.2;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, 'SW','unit','normalized');
text(ax, 0.7, 0.8, 'HF','unit','normalized');
text(ax, 0.7, 0.7, strcat(num2str(size(hfplt,1)),{' detections'}),'unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');

% subplot 4 of figure i
ax = f.ax(4);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxlf = 0.5;
dylf = 0.5;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
    xran,yran,dxlf,dylf);
dumlf = density1d(density1d(:,3)>0, :);
dumlf(dumlf(:,3)>1, :) = [];
scatter(ax,dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumlf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumlf(dumlf(:,3)==1, :) = [];
scatter(ax,dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, 'SW','unit','normalized');
text(ax, 0.7, 0.8, 'LF','unit','normalized');
text(ax, 0.7, 0.7, strcat(num2str(size(lfplt,1)),{' detections'}),'unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');
print(f.fig,'-depsc2',strcat(rstpath,'/density_rtm_NE+SW.eps'));


%% RTM density, NE vs. SW
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.1, 0.55, 0.4, 0.4]);
set(f.ax(2), 'position', [ 0.55, 0.55, 0.4, 0.4]);
set(f.ax(3), 'position', [ 0.1, 0.1, 0.4, 0.4]);
set(f.ax(4), 'position', [ 0.55, 0.1, 0.4, 0.4]);

msizehf = 4;
msizelf = 10;

% subplot 1 of figure i
hfplt = mighfse;
lfplt = miglfse;

ax = f.ax(1);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
%     ind = find(density1d>=0);
%     scatter(ax,xloc(ind),yloc(ind), 4, log(density1d(ind)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
%     aaa = density2d';
%     bbb = aaa(aaa>=1);
%     ccc = log(bbb);
%     imagesc(ax,xran+dx, yran+dy, ccc);
% create a density matrix to store the number of detections in each small grid
dxhf = 0.2;
dyhf = 0.2;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, 'SE','unit','normalized');
text(ax, 0.7, 0.8, 'HF','unit','normalized');
text(ax, 0.7, 0.7, strcat(num2str(size(hfplt,1)),{' detections'}),'unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');

% subplot 2 of figure i
ax = f.ax(2);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxlf = 0.5;
dylf = 0.5;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
    xran,yran,dxlf,dylf);
dumlf = density1d(density1d(:,3)>0, :);
dumlf(dumlf(:,3)>1, :) = [];
scatter(ax,dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumlf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumlf(dumlf(:,3)==1, :) = [];
scatter(ax,dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, 'SE','unit','normalized');
text(ax, 0.7, 0.8, 'LF','unit','normalized');
text(ax, 0.7, 0.7, strcat(num2str(size(lfplt,1)),{' detections'}),'unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');

% subplot 3 of figure i
hfplt = mighfnw;
lfplt = miglfnw;

ax = f.ax(3);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
%     ind = find(density1d>=0);
%     scatter(ax,xloc(ind),yloc(ind), 4, log(density1d(ind)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
%     aaa = density2d';
%     bbb = aaa(aaa>=1);
%     ccc = log(bbb);
%     imagesc(ax,xran+dx, yran+dy, ccc);
% create a density matrix to store the number of detections in each small grid
dxhf = 0.2;
dyhf = 0.2;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, 'NW','unit','normalized');
text(ax, 0.7, 0.8, 'HF','unit','normalized');
text(ax, 0.7, 0.7, strcat(num2str(size(hfplt,1)),{' detections'}),'unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');

% subplot 4 of figure i
ax = f.ax(4);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxlf = 0.5;
dylf = 0.5;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
    xran,yran,dxlf,dylf);
dumlf = density1d(density1d(:,3)>0, :);
dumlf(dumlf(:,3)>1, :) = [];
scatter(ax,dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumlf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumlf(dumlf(:,3)==1, :) = [];
scatter(ax,dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
c.Label.String = strcat('log_{10}(N) of detections');
%     c.Label.FontSize = 11;
%     caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
% plot(ax,xdiv,ydiv,'b--','linew',2);
text(ax, 0.7, 0.9, 'NW','unit','normalized');
text(ax, 0.7, 0.8, 'LF','unit','normalized');
text(ax, 0.7, 0.7, strcat(num2str(size(lfplt,1)),{' detections'}),'unit','normalized');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
hold(ax,'off');
print(f.fig,'-depsc2',strcat(rstpath,'/density_rtm_SE+NW.eps'));















