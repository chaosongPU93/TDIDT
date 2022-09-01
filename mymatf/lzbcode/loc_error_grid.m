% function loc_error_grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to plot the grid of location error in map that is 
% transformed from the grid of error in +/-2 samples (gaussian or uniform
% distributed), using the reference fam 043. 
% The purpose is to see the variation of error with the azimuth, as previously
% we cared more about the error 'ellipse' rather than a grid of error (eg. in
% 'filtering_effect_mag.m', or 'loc_resolution_v2.m')
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/08/24
% Last modified date:   2021/08/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
format short e   % Set the format to 5-digit floating point
clc
clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);

% set path
workpath = getenv('MHOME');
lzbpath = strcat(workpath, '/Seisbasics/hypoinverse/lzbrst');


%% LZB trio
% fams used
nfampool = [
            '002';
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
            ];
nfam = size(nfampool, 1);

% relative arrival time to main station for each fam
contlzboff = [
%               14 -22;
              -14 -1;
%               -23 7;
%               2 -15;
%               -25 -3;
%               -26 10;
%               -26 6;
%               -8 -1;
%               -4 -6;
%               -15 4;
%               -18 6;
%               -32 9;
%               -32 5;
              ];

          %%% fc is inverted with (lf-hf_12, lf-hf_13) relative to all fams
%%% the same content as /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/LZBallfamfc, which is
%%% inverted by hypoinverse from /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/chat00p.fc
fclzb = [
    -123.509500 48.457000 38.0200;
    -123.801000 48.561333 36.2500;
    -123.896333 48.594000 35.7800;
    -123.638833 48.474000 36.7200;
    -123.797167 48.440333 34.8100;
    -123.925000 48.599333 35.5600;
    -123.898667 48.555833 35.2500;
    -123.772167 48.575333 36.7300;
    -123.734833 48.562667 36.8900;
    -123.837500 48.587500 36.3200;
    -123.867833 48.590000 36.0100;
    -123.930667 48.545167 34.8600;       % 006
    -123.892500 48.467333 34.3600;       % 001
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VERSION 1 for the locations of lfe families, by inversion from hypoinverse
%%% this is inverted from (0,0) of all fams, same order, location of control points
%%% inverted from /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/chat00p.reffam
contlzb = [
    -123.492667 48.451500 38.1400;
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
    ];

%%% load the location resolution points, +-2 samples
reslzb = load(fullfile(lzbpath,'evtloc.13fam_locres_hf2spl'));
ifam = 2;  
reslzb = reslzb((ifam-1)*4+1: ifam*4, :);
[dx, dy] = absloc2relaloc(reslzb(:,1),reslzb(:,2),contlzb(ifam,1),contlzb(ifam,2));
vecreslzb = [dx dy];

%%% load the location error grid, +-2 samples
errlzb = load(fullfile(lzbpath,'evtloc.offset_043_locerr2_40_gaus'));
% errlzb = load(fullfile(lzbpath,'evtloc.offset_043_locerr2_40_unif'));
[dx, dy] = absloc2relaloc(errlzb(:,1),errlzb(:,2),contlzb(ifam,1),contlzb(ifam,2));
vecerrlzb = [dx dy];

%%% load the offset error grid in samples
errlzbspl = load(fullfile(lzbpath, 'offset_043_locerr2_40_gaus'));

%%% convert each filterring correction to km relative to its own fam
ifam = 2;   % choose fam 043 as the representative
% ifam = 10;   % choose fam 147 as the representative
[dx, dy] = absloc2relaloc(fclzb(ifam,1),fclzb(ifam,2),contlzb(ifam,1),contlzb(ifam,2));
fcvecrep = [dx dy];
fcmagrep = sqrt(dx.^2+dy.^2);
tmp = rad2deg(atan2(dx, dy));
if tmp < 0
    fcazirep = 360+tmp;
else
    fcazirep = tmp;
end

vecerrlzb = vecerrlzb+fcvecrep;
vecreslzb = vecreslzb+fcvecrep;
fixpt = [0 0];
% fixpt = fixpt + fcvecrep;

%%
f.fig=figure;
f.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

nrow = 1;
ncol = 3;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.4, 0.88]);
set(f.ax(2), 'position', [ 0.58, 0.1, 0.4, 0.88]);

ax = f.ax(1);
hold(ax,'on');
% sample error indicated by +-2 sample grid
scatter(ax,errlzbspl(:,3),errlzbspl(:,4),10,'k','filled');
[k1,v1] = boundary(errlzbspl(:,3),errlzbspl(:,4),0);
plot(ax, errlzbspl(k1,3),errlzbspl(k1,4),'b','linew',1.5);
text(ax, 0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax, 0.95,0.25, 'Gaussian-','fontsize',11,'unit','normalized',...
    'horizontalalignment','right');
text(ax, 0.95,0.18, 'distributed','fontsize',11,'unit','normalized',...
    'horizontalalignment','right');

ax.Box='on';
grid(ax,'on');
ax.GridLineStyle = '--';
axis(ax,'equal');
axis(ax,[-3 3 -3 3]);
xticks(ax, -3: 1: 3);
yticks(ax, -3: 1: 3);
xlabel(ax,'Offset between TWKB-LZB (samples)','fontsize',12)
ylabel(ax,'Offset between TWKB-MGCB (samples)','fontsize',12)

    
ax = f.ax(2);
hold(ax,'on');    
% loc error indicated by +-2 sample grid
scatter(ax,vecerrlzb(:,1),vecerrlzb(:,2),6,'k','filled');
% loc resolution indicated by +-2 sample cross
x = [vecreslzb(1,1), vecreslzb(2,1)];
y = [vecreslzb(1,2), vecreslzb(2,2)];
plot(ax,x,y,'k--','linewidth',1.5);

x = [vecreslzb(3,1), vecreslzb(4,1)];
y = [vecreslzb(3,2), vecreslzb(4,2)];
plot(ax,x,y,'k--','linewidth',1.5);

[k1,v1] = boundary(vecerrlzb(:,1),vecerrlzb(:,2),0);
% scatter(ax, vecerrlzb(k1,1),vecerrlzb(k1,2),20,'b','filled','o');
plot(ax, vecerrlzb(k1,1),vecerrlzb(k1,2),'b','linew',1.5);

x = [0; fcvecrep(1)];
y = [0; fcvecrep(2)];
plot(ax,x,y,'color',[0.6 0.6 0.6],'linewidth',2);
    
enex = [fixpt(1)-2; fixpt(1)+2];
aziene = 90-45/2;
enea = tand(90-aziene);
eneb = fixpt(2)-enea*(fixpt(1));
eney = linefcn(enex,enea,eneb);
plot(ax, enex,eney,'r','linew',1);
text(ax, 2.1, 1, 'ENE-WSW', 'fontsize',10);

azie = 90;
ea = tand(90-azie);
eb = fixpt(2)-ea*(fixpt(1));
ey = linefcn(enex,ea,eb);
plot(ax, enex,ey,'r','linew',1);
text(ax, 2.1, 0, 'E-W', 'fontsize',10);

text(ax, 0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax, -0.5,-1.2,nfampool(ifam,:),'FontSize',11,'EdgeColor','k','Margin',2);
ax.Box='on';
grid(ax,'on');
ax.GridLineStyle = '--';
axis(ax,'equal');
axis(ax,[-5 9 -2 12]);
xticks(ax, -5: 1: 9);
yticks(ax, -2: 1: 12);
xlabel(ax,'E (km)','fontsize',12)
ylabel(ax,'N (km)','fontsize',12)

patarea = [0 3;
           8.7 3;
           8.7 11.7;
           0 11.7;
           0 3;];
patch(ax,patarea(:,1),patarea(:,2),'w','edgecolor',[0.3 0.3 0.3]);
hold(ax,'off');

%%%% a small inset to show the difference in projection of fe correction
% projection of center of the ellipse onto E and ENE direction
fcmagenerep = fcmagrep*cos(deg2rad(360-fcazirep+aziene));
fcmageneortrep = fcmagrep*sin(deg2rad(360-fcazirep+aziene));
fcmagerep = fcmagrep*cos(deg2rad(360-fcazirep+azie));
fcmageortrep = fcmagrep*sin(deg2rad(360-fcazirep+azie));
fcdiffrep = fcmagerep-fcmagenerep;
fcdiffrepa = abs(fcmagerep)-abs(fcmagenerep);
for i = 1: size(vecerrlzb,1)
    % decompose the absolute mag on its own azimuth to 2 components on ENE azi and its ortho
    x = vecerrlzb(i,1);
    y = vecerrlzb(i,2);
    fcmag(i) = sqrt(x.^2+y.^2);
    tmp = rad2deg(atan2(x, y));
    if tmp < 0
        fcazi(i) = 360+tmp;
    else
        fcazi(i) = tmp;
    end
    fcmagene(i) = fcmag(i)*cos(deg2rad(360-fcazi(i)+aziene));
    fcmageneort(i) = fcmag(i)*sin(deg2rad(360-fcazi(i)+aziene));
    fcmage(i) = fcmag(i)*cos(deg2rad(360-fcazi(i)+azie));
    fcmageort(i) = fcmag(i)*sin(deg2rad(360-fcazi(i)+azie));
    fcdiff(i) = fcmage(i) - fcmagene(i);
    fcdiffa(i) = abs(fcmage(i))-abs(fcmagene(i));
%     x0 = loclfe(i,1);
%     y0 = loclfe(i,2);
%     [loclfenew(i,1),loclfenew(i,2)] = coordinate_rot(x0,y0,-(aziene-90),0,0);
end

set(f.ax(3), 'position', [ 0.775, 0.505, 0.175, 0.32]);
ax = f.ax(3);
hold(ax,'on');
ax.Box='on';
% grid(ax,'on');
% ax.GridLineStyle = '--';
ax.FontSize = 7;
binw = 0.2;
fcdiff = fcdiff-fcdiffrep;
h=histogram(ax,fcdiff,'binw',0.2,'facec',[0.6 0.6 0.6],'facea',0.8);  %,'normalization','pdf'
% plot(ax, [fcdiffrep fcdiffrep],ax.YLim,'b--','linewidth',2);
[mu,sigma,muCI,sigmaCI] = normfit(fcdiff', 0.05);
xx = mu-3:0.01:mu+3;
synpdf = normpdf(xx,mu,sigma);
synct = synpdf*binw*length(fcdiff);      % synthetic counts
plot(ax,xx, synct, 'k-', 'linew', 1.5);
xlabel(ax,'Differential error E-ENE (km)','fontsize',8)
ylabel(ax,'Count','fontsize',8)
% if ifam == 2 
%     xlim(ax, [-6 0]);
% elseif ifam == 10
%     xlim(ax, [-5 1]);
% end
xlim(ax, [-3 3]);
xticks(ax, -2:1:2);
ylim(ax, [0 700]);

print(f.fig,'-dpdf',...
    strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2020/Figure/loc_error_grid_gaussian',...
    nfampool(ifam,:),'.pdf'));

figure
ax = gca;
hold(ax, 'on');
binw = 0.2;
h=histogram(ax,fcdiff,'binw',0.2,'facec',[0.6 0.6 0.6],'facea',0.8);  %,'normalization','pdf'
% plot(ax, [fcdiffrep fcdiffrep],ax.YLim,'b--','linewidth',2);
[mu,sigma,muCI,sigmaCI] = normfit(fcdiff', 0.05);
xx = mu-3:0.01:mu+3;
synpdf = normpdf(xx,mu,sigma);
synct = synpdf*binw*length(fcdiff);      % synthetic counts
plot(ax,xx, synct, 'k-', 'linew', 1.5);
% if ifam == 2 
%     xlim(ax, [-6 0]);
% elseif ifam == 10
%     xlim(ax, [-5 1]);
% end
xlim(ax, [-3 3]);
xticks(ax, -2:1:2);
ylim(ax, [0 700]);
xlabel(ax,'Differential error E-ENE (km)','fontsize',12)
ylabel(ax,'Count','fontsize',12)


%%
%%%% try the N-S direction as well, get the difference between E and N
azin = 0;
fcmagnrep = fcmagrep*cos(deg2rad(360-fcazirep+azin));
fcdiffrepen = fcmagerep-fcmagnrep;
for i = 1: size(vecerrlzb,1)
    % decompose the absolute mag on its own azimuth to 2 components on ENE azi and its ortho
    x = vecerrlzb(i,1);
    y = vecerrlzb(i,2);
    fcmag(i) = sqrt(x.^2+y.^2);
    tmp = rad2deg(atan2(x, y));
    if tmp < 0
        fcazi(i) = 360+tmp;
    else
        fcazi(i) = tmp;
    end
    fcmagn(i) = fcmag(i)*cos(deg2rad(360-fcazi(i)+azin));
    fcdiffen(i) = fcmage(i) - fcmagn(i);
end
fcdiffen = fcdiffen - fcdiffrepen;
figure
ax = gca;
hold(ax, 'on');
binw = 0.2;
h=histogram(ax,fcdiffen,'binw',0.2,'facec',[0.6 0.6 0.6],'facea',0.8);  %,'normalization','pdf'
% plot(ax, [fcdiffrep fcdiffrep],ax.YLim,'b--','linewidth',2);
[mu2,sigma2,muCI,sigmaCI] = normfit(fcdiffen', 0.05);
xx = mu2-5:0.01:mu2+5;
synpdf = normpdf(xx,mu2,sigma2);
synct = synpdf*binw*length(fcdiffen);      % synthetic counts
plot(ax,xx, synct, 'k-', 'linew', 1.5);
% if ifam == 2 
%     xlim(ax, [-6 0]);
% elseif ifam == 10
%     xlim(ax, [-5 1]);
% end
xlim(ax, [-5 5]);
xlabel(ax,'Differential error E-N (km)','fontsize',12)
ylabel(ax,'Count','fontsize',12)

%%%% try the NNW direction as well, get the difference between NNW and ENE
azinnw = 360-45/2;
fcmagnnwrep = fcmagrep*cos(deg2rad(360-fcazirep+azinnw));
fcdiffrepnnw = fcmagnnwrep-fcmagenerep;
for i = 1: size(vecerrlzb,1)
    % decompose the absolute mag on its own azimuth to 2 components on ENE azi and its ortho
    x = vecerrlzb(i,1);
    y = vecerrlzb(i,2);
    fcmag(i) = sqrt(x.^2+y.^2);
    tmp = rad2deg(atan2(x, y));
    if tmp < 0
        fcazi(i) = 360+tmp;
    else
        fcazi(i) = tmp;
    end
    fcmagnnw(i) = fcmag(i)*cos(deg2rad(360-fcazi(i)+azinnw));
    fcdiffnnw(i) = fcmagnnw(i) - fcmagene(i);
end
fcdiffnnw = fcdiffnnw - fcdiffrepnnw;
figure
ax = gca;
hold(ax, 'on');
binw = 0.2;
h=histogram(ax,fcdiffnnw,'binw',0.2,'facec',[0.6 0.6 0.6],'facea',0.8);  %,'normalization','pdf'
% plot(ax, [fcdiffrep fcdiffrep],ax.YLim,'b--','linewidth',2);
[mu3,sigma3,muCI,sigmaCI] = normfit(fcdiffnnw', 0.05);
xx = mu3-5:0.01:mu3+5;
synpdf = normpdf(xx,mu3,sigma3);
synct = synpdf*binw*length(fcdiffnnw);      % synthetic counts
plot(ax,xx, synct, 'k-', 'linew', 1.5);
% if ifam == 2 
%     xlim(ax, [-6 0]);
% elseif ifam == 10
%     xlim(ax, [-5 1]);
% end
xlim(ax, [-5 5]);
xlabel(ax,'Differential error NNW-ENE (km)','fontsize',12)
ylabel(ax,'Count','fontsize',12)









