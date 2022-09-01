% function loc_resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to obtain the location resolution of the detection in
% HF and LF. In time domain, HF detection is precise to 1/4 sample in 40 sps,
% on the contrary, LF is precise to 1 sample in 20 sps. Thus the resolution in
% time is 1:8. Choose several control points, (off12, off13) = (+-1, 0) and
% (0, +-1) to obtain the resulting location could illustrate the location
% resolution.
% It seems that at least in fam 002 with PGC detector, the (-1/4,0) and (0, -1/4)
% in 40 sps are indistinguishable. But in LZB detector, it should n't have this
% issue. The reason comes from the grid size of the input slab. However, the
% slab itself is interpolated to a denser grid. Doesm't really make sense to
% use a denser interpolation.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2020/08/13
% Last modified date:   2021/01/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
% close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);

% set path
workpath = getenv('MHOME');
pgcpath = strcat(workpath, '/Seisbasics/hypoinverse/forsummary');
lzbpath = strcat(workpath, '/Seisbasics/hypoinverse/lzbrst');

%% PGC trio
%%% this is inverted from (0,0) relative to fam 002, location of control points
contpgc = [-123.585000 48.436667 36.8800];

%%% load the location resolution points
respgchf = load(fullfile(pgcpath,'eventloc.002.locres.hf'));
respgclf = load(fullfile(pgcpath,'eventloc.002.locres.lf'));

%%% convert to km
[dx, dy] = absloc2relaloc(respgchf(1),respgchf(2),contpgc(1),contpgc(2));
vecpgchf = [dx dy];
magpgchf = sqrt(dx.^2+dy.^2);

[dx, dy] = absloc2relaloc(respgclf(1),respgclf(2),contpgc(1),contpgc(2));
vecpgclf = [dx dy];
magpgclf = sqrt(dx.^2+dy.^2);

%% LZB trio
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
    '017'];
nfam = length(nfampool);

if nfam == 11
    %%% this is inverted from (0,0) of all fams, same order, location of control points
    %%% inverted from /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/chat00p.reffam
    % this is inverted from (0,0) of all fams, same order, location of control points
    contlzb = [-123.492667 48.451500 38.1400;
        -123.772167 48.493000 35.5900;
        -123.863167 48.528167 35.2100;
        -123.603333 48.440167 36.7100;
        -123.800167 48.408833 34.5200;
        -123.893333 48.536500 35.0700;
        -123.864500 48.498667 34.8800;
        -123.753333 48.525667 36.2000;
        -123.703667 48.502667 36.4100;
        -123.814333 48.538667 35.7900;
        -123.838500 48.544833 35.6600];
    
elseif nfam == 12
    % this is inverted from (0,0) of all fams, same order, location of control points
    contlzb = [-123.492667 48.451500 38.1400;
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
        -123.879667 48.446167 34.2600];
end

%%% convert lfe location to km relative to 043
loclfe = zeros(size(contlzb,1),2);
for i = 1: size(contlzb,1)
    [dx, dy] = absloc2relaloc(contlzb(i,1),contlzb(i,2),contlzb(2,1),contlzb(2,2));
    loclfe(i,:) = [dx dy];
end

%%% load the location resolution points
reslzbhf = load(fullfile(lzbpath,'LZBall12famlocres.hf40'));
reslzblf = load(fullfile(lzbpath,'LZBall12famlocres.lf'));
reslzbhf = reslzbhf(1: 4*nfam, :);
reslzblf = reslzblf(1: 4*nfam, :);

%%% convert each location resolution to km relative to its own fam
veclzbhf = zeros(size(reslzbhf,1),2);
maglzbhf = zeros(size(reslzbhf));
azilzbhf = zeros(size(reslzbhf));
for i = 1: size(reslzbhf,1)
    [dx, dy] = absloc2relaloc(reslzbhf(i,1),reslzbhf(i,2),contlzb(ceil(i/4),1),contlzb(ceil(i/4),2));
    veclzbhf(i,:) = [dx dy];
    maglzbhf(i,1) = sqrt(dx.^2+dy.^2);
    tmp = rad2deg(atan2(dy,dx));
    if tmp < 0
        tmp = 360+tmp;
    end
    tmp = 360 - tmp;
    tmp = tmp+90;
    if tmp > 360
        tmp = tmp-360;
    end
    azilzbhf(i,1) = tmp;
end

veclzblf = zeros(size(reslzblf,1),2);
maglzblf = zeros(size(reslzblf));
azilzblf = zeros(size(reslzblf));
for i = 1: size(reslzblf,1)
    [dx, dy] = absloc2relaloc(reslzblf(i,1),reslzblf(i,2),contlzb(ceil(i/4),1),contlzb(ceil(i/4),2));
    veclzblf(i,:) = [dx dy];
    maglzblf(i,1) = sqrt(dx.^2+dy.^2);
    tmp = rad2deg(atan2(dy,dx));
    if tmp < 0
        tmp = 360+tmp;
    end
    tmp = 360 - tmp;
    tmp = tmp+90;
    if tmp > 360
        tmp = tmp-360;
    end
    azilzblf(i,1) = tmp;
end


%% plot
xran = [-15 25];
yran = [-20 20];

f.fig=figure;
f.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.1, 0.64, 0.32, 0.32]);
set(f.ax(2), 'position', [ 0.5, 0.64, 0.32, 0.32]);

ax = f.ax(1);
hold(ax,'on');

for i = 1: nfam
    x = [loclfe(i,1)+ veclzbhf((i-1)*4+1,1), loclfe(i,1)+ veclzbhf((i-1)*4+2,1)];
    y = [loclfe(i,2)+ veclzbhf((i-1)*4+1,2), loclfe(i,2)+ veclzbhf((i-1)*4+2,2)];
    plot(ax,x,y,'color',[0.6 0.6 0.6],'linewidth',1.5);
    
    x = [loclfe(i,1)+ veclzbhf((i-1)*4+3,1), loclfe(i,1)+ veclzbhf((i-1)*4+4,1)];
    y = [loclfe(i,2)+ veclzbhf((i-1)*4+3,2), loclfe(i,2)+ veclzbhf((i-1)*4+4,2)];
    plot(ax,x,y,'color',[0.6 0.6 0.6],'linewidth',1.5);
    
end

ax.Box='on';
grid(ax,'on');
ax.GridLineStyle = '--';
axis(ax,'equal');
xlim(ax,xran);
ylim(ax,yran);

xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
hold(ax,'off');


ax = f.ax(2);
hold(ax,'on');

for i = 1: nfam
    x = [loclfe(i,1)+ veclzblf((i-1)*4+1,1), loclfe(i,1)+ veclzblf((i-1)*4+2,1)];
    y = [loclfe(i,2)+ veclzblf((i-1)*4+1,2), loclfe(i,2)+ veclzblf((i-1)*4+2,2)];
    plot(ax,x,y,'color',[0.6 0.6 0.6],'linewidth',1.5);
    
    x = [loclfe(i,1)+ veclzblf((i-1)*4+3,1), loclfe(i,1)+ veclzblf((i-1)*4+4,1)];
    y = [loclfe(i,2)+ veclzblf((i-1)*4+3,2), loclfe(i,2)+ veclzblf((i-1)*4+4,2)];
    plot(ax,x,y,'color',[0.6 0.6 0.6],'linewidth',1.5);
    
end

ax.Box='on';
grid(ax,'on');
ax.GridLineStyle = '--';
axis(ax,'equal');
xlim(ax,xran);
ylim(ax,yran);


xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
hold(ax,'off');

%%% save figure
print(f.fig,'-dpdf',strcat(lzbpath,'/LZB.loc_resolution40.pdf'));







