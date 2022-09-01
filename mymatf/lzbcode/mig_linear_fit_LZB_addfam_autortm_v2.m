% function mig_linear_fit_LZB_addfam_autortm_v2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to use the poropagation direction estimate fit the hf
% the migration linearly to get the propagation slope and error distribution.
% Use the new automatic time ranges of migrations that come from 'identify_RTMs_v2.m'
%
%    2020/09/10, i am adding a fam 001 back again for last try, 13 fam in total now
%                this is for adding new fams 006, 001
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/01/07
% Last modified date:   2021/01/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all
clc

% set(0,'DefaultFigureVisible','on');
set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

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

% this is new time ranges that come from 'identify_RTMs_v2.m'
trange = [
    2004195   4.1773e+04   4.3362e+04   % combine, y
    2004197   4.0974e+04   4.2205e+04   % speed direct 90, y, larger pear
    2004197   5.4453e+04   5.5867e+04   % speed direct 205, y, but use rmse 220
    2004197   5.6270e+04   5.7227e+04   % speed direct 175, y, but use rmse 185
    2004197   7.9949e+04   8.1453e+04   % speed direct 260, y, but use rmse 195
    2004197   8.3402e+04   8.4888e+04   % divided into 3, 1 discarded, y
    2004197   8.5320e+04   8.6659e+04   % divided into 3, 1 discarded, y
    2004198   5.9530e+03   9.7940e+03   % speed direct 270, y, but use rmse 235
    2004198   1.9151e+04   2.1616e+04   % acceptted, y
    2004198   4.2632e+04   4.3795e+04   % acceptted, y
    2004198   5.3506e+04   5.8050e+04   % acceptted, y
    2004198   6.2966e+04   6.3885e+04   % acceptted, y
    2004198   7.3540e+04   7.6029e+04   % acceptted, y 
    2004198   8.4222e+04   8.7248e+04   % combine, y
    2004199   4.1670e+03   6.2660e+03   % speed direct 280, y, larger pear
    2004199   4.2845e+04   4.3359e+04   % speed direct 270, y, but use rmse 250
    2004199   4.6340e+04   4.9126e+04   % combine, care speed direct, y, larger pear
    2004199   4.9744e+04   5.0614e+04   % speed direct 180, y, but use rmse 195
    2004199   8.0861e+04   8.2008e+04   % divided into 2, y, care speed direct, y, but use rmse 90
    2004199   8.2008e+04   8.3515e+04   % divided into 2, y, care speed direct, RECHECK PEAR, slightly larger pear    
    2004200   1.2600e+04   1.5799e+04   % modified time, y
    2004200   1.5914e+04   1.7125e+04   % acceptted, y
    2004200   1.9104e+04   1.9900e+04   % check OLD param. to determine, y, use OLD end time
    2004200   4.8109e+04   4.8666e+04   % acceptted, y
    2004201   4.3700e+02   1.6290e+03   % speed direct 225, y, but use rmse 260
    2004203   1.6586e+04   2.0776e+04   % combine them, y
    2005255   3.4400e+04   3.5612e+04   % speed direct 80, y, larger pear
    2005255   6.0381e+04   6.1557e+04   % speed direct 150, y, but use rmse 110
    2005255   6.7268e+04   6.8535e+04   % speed direct 80, y, larger pear
    2005255   8.4296e+04   8.5634e+04   % combine them, y
    2005256   3.6200e+03   4.9570e+03   % speed direct 265, y, larger pear
    2005256   8.0330e+03   9.8280e+03   % change end time, y
    2005256   1.3412e+04   1.4346e+04   % speed direct 65, y, but use rmse 80   
    2005256   1.9069e+04   2.0203e+04   % speed direct 195, y, larger pear
    2005256   2.1492e+04   2.2294e+04   % change start time, y
    2005256   2.8080e+04   2.9172e+04   % change start time, y
    2005256   3.4776e+04   3.5238e+04   % change start time, y
    2005256   3.6900e+04   3.8520e+04   % change times, y
    2005256   4.6440e+04   4.8421e+04   % change start time, y, check speed direc, y, a bit awful
    2005256   5.2020e+04   5.2596e+04   % change times, y
    2005256   6.2640e+04   6.4865e+04   % change start time, y, check speed direc, y, but use rmse 130
    2005256   7.2720e+04   7.3914e+04   % divided into 2, y, but check speed direc, y
    2005256   7.6381e+04   7.7940e+04   % change end time, y
    2005257   6.1200e+03   7.7460e+03   % change start time, y, but check speed direc, y
    2005257   1.2240e+04   1.3248e+04   % CHANGED time, st-3.68, y, use speed direc, larger pear, SUSPICOUS
    2005257   1.6920e+04   1.8655e+04   % change start time, y
    2005257   2.1070e+04   2.5365e+04   % accepted, y 
    2005257   3.9000e+04   4.5418e+04   % change start time, y
    2005257   6.1449e+04   6.4012e+04   % accepted, y
    2005257   7.3440e+04   7.8840e+04   % change end time, y, but check speed direc, y, but use rmse 235
    2005257   8.2159e+04   8.3380e+04   % speed direct 350, y, larger pear
    2005258   3.5000e+04   3.6540e+04   % divided into 3, y, but check speed direc, RECHECK PEAR, slightly larger pear, SUSPICOUS
    2005258   3.6540e+04   3.8160e+04   % divided into 3, y
    2005258   3.8160e+04   4.0021e+04   % divided into 3, y, but check speed direc y, slightly larger pear 
    2005259   1.5400e+02   1.6340e+03   % speed direct 75, y, larger pear
    2005259   1.9560e+03   4.0460e+03   % speed direct 185, y, but use rmse 140
    2005259   5.2320e+03   7.7140e+03   % check OLD param. to determine, y, but check speed direc, larger pear, use old direc 255
    2005259   3.7000e+04   4.1293e+04   % speed direct 85, RECHECK PEAR, y, but use rmse 115
    2005260   2.6510e+03   3.4380e+03   % speed direct 80, y, but use rmse 115
    2005260   3.6040e+03   4.5000e+03   % divided into 2, y
    2005260   4.6080e+03   6.6780e+03   % divided into 2, y, but check speed direc, y, larger pear
    2005260   5.6300e+04   5.8200e+04   % use old times, y
    2005260   9.4090e+03   1.0357e+04   % speed direct 55, y, but use rmse 115
    ];

% hardcopy of angbest from either min. rmse or max. slope of accepted ones ntranacc
angbest = [
    250   % rmse
    90   % speed
    220   % rmse
    185   % rmse
    195   % rmse
    225   % rmse
    270   % rmse
    235   % rmse
    200   % rmse
    230   % rmse
    80   % rmse
    65   % rmse
    235   % rmse
    245   % rmse
    280   % speed
    250   % rmse
    90   % speed
    195   % rmse
    90   % rmse
    195   % speed
    130   % rmse
    70   % rmse
    255   % rmse
    250   % rmse
    260   % rmse
    120   % rmse
    80    % speed
    110    % rmse
    80    % speed
    240    % rmse
    265    % speed
    220    % rmse
    80    % rmse
    195    % speed
    115    % rmse
    105    % rmse
    260    % rmse
    160    % rmse
    65    % awful, speed
    260    % rmse
    130    % rmse
    230    % speed
    250    % rmse
    235    % speed
    70    % SUSPICOUS, speed
    115    % rmse
    250    % rmse
    165    % rmse
    245    % rmse
    235    % rmse
    350    % speed
    60    % SUSPICOUS, speed
    155    % rmse
    125    % speed
    75    % speed
    140    % rmse
    255    % old direc 255, neither
    115    % rmse
    115    % rmse
    250    % rmse
    220    % speed
    235    % rmse
    25    % rmse
    ];

indspeed = [
    2
    15
    17
    20
    27
    29
    31
    34
    39
    42
    44
    45
    51
    52
    54
    55
    61
    ];

% angbest(32) = 65; %67.5 112.5
% angbest(33) = 245; %247.5 270


%% linear fitting and median distance 
resprophf = nan(length(trange),1200);
resproplf = nan(length(trange),500);

medallhf = nan(length(trange),4);
medalllf = nan(length(trange),4);
ranvechf = nan(length(trange),2);
ranveclf = nan(length(trange),2);
ranvechf95 = nan(length(trange),2);
ranveclf95 = nan(length(trange),2);
ranvechf98 = nan(length(trange),2);
ranveclf98 = nan(length(trange),2);
ranvechf99 = nan(length(trange),2);
ranveclf99 = nan(length(trange),2);

% array of vertical distance from LF from fitted HF line 
vertdistraw = nan(300 , length(trange));    % raw vertical distance
vertdistwt = nan(300 , length(trange));     % weighted vertical distance
weight = nan(300 , length(trange));     % weights from regression
indlfindpen = nan(300 , length(trange));    % store the index of indepedent LF detections
indltdiv = nan(300 , length(trange));    % store the index that is lower than the division
indgtdiv = nan(300 , length(trange));    % store the index that is greater than the division

% percentage of detections from specific fam relative to all in each RTM
famperchf = zeros(length(trange), nfam);
famperclf = zeros(length(trange), nfam);

% choose the end of RTM 42 as the division point, get its location x y
x0 = -3.5474e+00;
y0 = 2.0923e+00;
% the division line passes through this point and orientates sse
xdiv = -12:0.01:0;
a1 = 1/tan(deg2rad(180-22.5));
b1 = y0-a1*x0;
ydiv = linefcn(xdiv,a1,b1);

xran = [-16 24];
yran = [-20 20];

angle = 0:5:355;

slopehf = zeros(size(trange,1),length(angle));
rmsehf = zeros(size(trange,1),length(angle));
angrmsehf = zeros(size(trange,1),1);
angslopehf = zeros(size(trange,1),1);

% create fit object with free constraints
fttpfree = fittype( @(a,b,x) a*x+b);
        
for i = 1: length(trange)
% for i = 1: 211
%     i=32;
    disp(trange(i,:));
    indhf = find(hftime(:,13)==trange(i,1) & hftime(:,15)>=trange(i,2) & ...
                 hftime(:,15)<=trange(i,3));
    mighf = hftime(indhf,:);
    
    indlf = find(lftime(:,13)==trange(i,1) & lftime(:,15)>=trange(i,2) & ...
                 lftime(:,15)<=trange(i,3));
    miglf = lftime(indlf,:); 
    
    %%% find best angle, now is mainly to get the variation of se and slope with the trial angle
    for iang = 1: length(angle)
        %%% propagation trial of hf
        mighfdum = mighf;
        for j = 1: size(mighf,1)
            x0 = mighfdum(j,1);
            y0 = mighfdum(j,2);
            [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),0,0);
            mighfdum(j,1) = newx;
            mighfdum(j,2) = newy;
        end        
        % linear robust least square
        [fitobj,gof,~] = fit(mighfdum(:,15)/3600, mighfdum(:,1),fttpfree,'Robust','Bisquare',...
                                  'StartPoint',[1 1]);
        % output fit parameters
        coef = coeffvalues(fitobj);
        slopehf(i,iang) = coef(1);
        rmsehf(i,iang) = gof.rmse;
            
    end

    %%% best angle estimate from hf
    ind = find(slopehf(i,:)>0);
    ind3 = find(rmsehf(i,ind)==min(rmsehf(i,ind)));     % rmse, Root Mean Squared Error, the smaller, the better
    if length(ind3) > 1
        disp(strcat({'multiple angles have same rmse: '},num2str(angle(ind3))));
    end
    angrmsehf(i) = angle(ind(ind3(1)));

    ind6 = find(slopehf(i,:)==max(slopehf(i,:))); % one with the largest slope, i.e., migrating speed
    if length(ind6) > 1
        disp(strcat({'multiple angles have same slope:'},num2str(angle(ind6))));
    end
    angslopehf(i) = angle(ind6(1));
    
    %%% Actually used is the pre-determined best prop direc to do the fitting
    mighfdum = mighf;
    for j = 1: size(mighf,1)
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        mighfdum(j,1) = newx;
        mighfdum(j,2) = newy;
    end

    miglfdum = miglf;
    for j = 1: size(miglf,1)
        x0 = miglfdum(j,1);
        y0 = miglfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        miglfdum(j,1) = newx;
        miglfdum(j,2) = newy;
    end
    
    %%% define and position the figure frame and axes of each plot
    f.fig=figure;
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 9;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    nrow = 2;
    ncol = 2;    
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
    end

    %%% reposition
    set(f.ax(1), 'position', [ 0.1, 0.64, 0.32, 0.32]);
    set(f.ax(2), 'position', [ 0.5, 0.64, 0.32, 0.32]);
    set(f.ax(3), 'position', [ 0.1, 0.35, 0.32, 0.22]);
    set(f.ax(4), 'position', [ 0.5, 0.35, 0.32, 0.22]);
%     set(f.ax(5), 'position', [ 0.1, 0.078, 0.32, 0.22]);
%     set(f.ax(6), 'position', [ 0.5, 0.078, 0.32, 0.22]);
    
    % subplot 1 of figure i
    hold(f.ax(1),'on');
    plot(f.ax(1),[-100 100],[0 0],'k--');
    plot(f.ax(1),[0 0],[-100 100],'k--');
    f.ax(1).FontSize = 9;
    mighf = sortrows(mighf,-15);
    scatter(f.ax(1),mighf(:,1),mighf(:,2), 20, mighf(:,15)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(1),'jet');
    c=colorbar(f.ax(1),'SouthOutside');
    pos = f.ax(1).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    juldate = num2str(trange(i,1));
    yr = str2double(juldate(1:4));
    date = str2double(juldate(5:end));
    a = jul2dat(yr,date);
    mo = a(1);
    if mo == 9 
        mo = {' Sep. '};
    elseif mo == 7
        mo = {' Jul. '};
    else
        mo = {' Mar. '};
    end
    day = num2str(a(2));
    yr = num2str(a(3));
    c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
%     c.Label.String = strcat(num2str(trange(i,1)),' of HF',' (hr)');
    c.Label.FontSize = 11;
    caxis(f.ax(1),[trange(i,2)/3600 trange(i,3)/3600])
    text(f.ax(1),0.85,0.1,'HF','FontSize',12,'unit','normalized');
    text(f.ax(1),0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(1),0.04,0.1,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(1).Box = 'on';
    grid(f.ax(1), 'on');
    axis(f.ax(1), 'equal');
    f.ax(1).GridLineStyle = '--';
    f.ax(1).XAxisLocation = 'top';
    medxhf = median(mighf(:,1));
    medyhf = median(mighf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medxhf-rotx medxhf+rotx];
    yvect = [medyhf-roty medyhf+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medxhf-rotx medxhf+rotx];
    yvect = [medyhf-roty medyhf+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
    xticks(f.ax(1),xran(1):5:xran(2));
    yticks(f.ax(1),yran(1):5:yran(2));    
    xlabel(f.ax(1),'E (km)','fontsize',11);
    ylabel(f.ax(1),'N (km)','fontsize',11);
    text(f.ax(1),0.5,0.1,num2str(i),'FontSize',12,'unit','normalized');
    text(f.ax(1),0.83,0.95,strcat(num2str(size(mighf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(1),0.83,0.90,strcat({'in '},num2str(trange(i,3)-trange(i,2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(mighf,1)/(trange(i,3)-trange(i,2)));
    text(f.ax(1),0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(1),0.83,0.80,strcat(num2str(angbest(i)),{'{\circ}'}),'FontSize',10,...
         'unit','normalized','horizontalalignment','center'); 
    hold(f.ax(1),'off');

    % subplot 2 of figure i
    hold(f.ax(2),'on');
    plot(f.ax(2),[-100 100],[0 0],'k--');
    plot(f.ax(2),[0 0],[-100 100],'k--');
    f.ax(2).FontSize = 9;
    miglf = sortrows(miglf,-15);
    scatter(f.ax(2),miglf(:,1),miglf(:,2), 20, miglf(:,15)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(2),'jet');
    c=colorbar(f.ax(2),'SouthOutside');
    pos = f.ax(2).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
%     c.Label.String = strcat(num2str(trange(i,1)),' of LF',' (hr)');
    c.Label.FontSize = 11;
    caxis(f.ax(2),[trange(i,2)/3600 trange(i,3)/3600])
    text(f.ax(2),0.85,0.1,'LF','FontSize',12,'unit','normalized');
    text(f.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(2),0.04,0.1,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(2).Box = 'on';
    grid(f.ax(2), 'on');
    axis(f.ax(2), 'equal');
    f.ax(2).GridLineStyle = '--';
    f.ax(2).XAxisLocation = 'top';
    medxlf = median(miglf(:,1));
    medylf = median(miglf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
    xticks(f.ax(2),xran(1):5:xran(2));
    yticks(f.ax(2),yran(1):5:yran(2));    
    xlabel(f.ax(2),'E (km)','fontsize',11);
    ylabel(f.ax(2),'N (km)','fontsize',11);
    text(f.ax(2),0.5,0.1,num2str(i),'FontSize',12,'unit','normalized');
    text(f.ax(2),0.83,0.95,strcat(num2str(size(miglf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(2),0.83,0.9,strcat({'in '},num2str(trange(i,3)-trange(i,2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(miglf,1)/(trange(i,3)-trange(i,2)));
    text(f.ax(2),0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(2),0.83,0.80,strcat(num2str(angbest(i)),{'{\circ}'}),'FontSize',10,...
         'unit','normalized','horizontalalignment','center');
    hold(f.ax(2),'off');
    
    % subplot 3 of figure i
    hold(f.ax(3),'on');
    f.ax(3).FontSize = 9;
    scatter(f.ax(3),mighfdum(:,15)/3600,mighfdum(:,1),20,[0.6 0.6 0.6],'filled','o',...
            'MarkerEdgeColor','k');    
    scatter(f.ax(3),miglfdum(:,15)/3600,miglfdum(:,1),20,[0.6 1 1],'filled','o',...
            'MarkerEdgeColor','k');
    % create fit object
    
%     indpool = [12;13;16;19;21;22;29;32;36;40];
%     indpool = [indpool; 3;5;7;17;20;23;25;26;30];
%     if ismember(i,indpool)
        mintime = max(min(mighfdum(:,15)), min(miglfdum(:,15)));
        maxtime = min(max(mighfdum(:,15)), max(miglfdum(:,15)));
        mighfdum2 = mighfdum(mighfdum(:,15)>=mintime & mighfdum(:,15)<=maxtime, :);
        miglfdum2 = miglfdum(miglfdum(:,15)>=mintime & miglfdum(:,15)<=maxtime, :);
%     else
%         mighfdum2 = mighfdum;
%         miglfdum2 = miglfdum;
%     end

    [fitobjhfprop,gofhfprop,outphfprop] = fit(mighfdum2(:,15)/3600, mighfdum2(:,1),fttpfree,'Robust',...
                                              'Bisquare','StartPoint',[1 1]);
    % output fit parameters
    coefprop = coeffvalues(fitobjhfprop);
    slopeprophf(i) = coefprop(1);
    intcptprophf(i) = coefprop(2);
    fitprophf = feval(fitobjhfprop,mighfdum2(:,15)/3600);
    plot(f.ax(3),mighfdum2(:,15)/3600,fitprophf,'-','linewidth',2,'color',[0.6 0.6 0.6]);
        
%     intcptproplf = ones(size(miglfdum(:,end)))\(miglfdum(:,1)-slopeprophf*miglfdum(:,end)/3600); % direct LS fit
%     fitproplf = slopeprophf*miglfdum(:,end)/3600+intcptproplf;
    a=slopeprophf(i);
    fttpfix = fittype( @(b,x) a*x+b);
    [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum2(:,15)/3600, miglfdum2(:,1),fttpfix,'Robust',...
                                              'Bisquare','StartPoint',1);
%     [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum(:,10)/3600, miglfdum(:,1),fttpfix,'Robust',...
%                                               'Bisquare','Weights',miglfdum(:,11));
    fitproplf = feval(fitobjlfprop,miglfdum2(:,15)/3600);
    intcptproplf(i) = coeffvalues(fitobjlfprop);
    offset(i) = intcptproplf(i)-intcptprophf(i);
    plot(f.ax(3),miglfdum2(:,15)/3600,fitproplf,'-','linewidth',2,'color',[0.6 1 1]);
    f.ax(3).Box = 'on';
    grid(f.ax(3), 'on');
    f.ax(3).GridLineStyle = '--';
    xran1 = [trange(i,2)/3600 trange(i,3)/3600];
%     yran1 = [round(min([mighfdum(:,1);miglfdum(:,1)]))-1 ...
%              round(max([mighfdum(:,1);miglfdum(:,1)]))+1];
    aa = round(prctile([mighfdum(:,1);miglfdum(:,1)], 98));
    bb = round(prctile([mighfdum(:,1);miglfdum(:,1)], 2));
    aap = aa + ceil((aa-bb)/6);
    bbp = bb - ceil((aa-bb)/4);
    yran1 = [bbp aap];
    xlim(f.ax(3),xran1);
    ylim(f.ax(3),yran1);
    text(f.ax(3),0.04,0.91,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    xlabel(f.ax(3),'Time (hr)','fontsize',11);
    ylabel(f.ax(3),'Dist. along prop. (km)','fontsize',11); 
    
    % compute the HF weights in robust linear regression, see NOTES above
    rhf = outphfprop.residuals;   % usual residuals
    x = [mighfdum2(:,15)/3600];
    hatmat = x*inv(x'*x)*x';
    h = zeros(size(hatmat,1),1);    % leverage of least square
    for jj = 1 : size(hatmat,1)
        h(jj) = hatmat(jj,jj);
    end    
    radj = rhf./sqrt(1-h);      % adjusted residuals
    K = 4.685;
    s = mad(rhf,1)/0.6745;
    u = radj/(K*s);
    wthf = zeros(length(u),1);    % rubust weight of next iteration
    for jj = 1 : length(u)
        if abs(u(jj)) < 1
            wthf(jj) = (1-(u(jj))^2)^2;
        else
            wthf(jj) = 0;
        end
    end
    
    % get the standard error of the estimated parameters, may indicate the compare the quality
    % of fitting cross different RTMs, NOTE that rmse, mse are definitely not a good measure
    % the obtained CI is comfirmed to be correct by comparing it with the 'fitobjhfprop'
    slopesehf = gofhfprop.rmse./sqrt(sum((x-mean(x)).^2));
    est = slopeprophf(i);
    slopeprophfCI(i,:) = confidence_interval_general(est,slopesehf,length(x)-2,95);
    interceptsehf = slopesehf.*sqrt(sum(x.^2)./length(x));
    est = intcptprophf(i);
    intcptprophfCI(i,:) = confidence_interval_general(est,interceptsehf,length(x)-2,95);
    sehf(i, :) = [gofhfprop.rmse slopesehf interceptsehf];     
    
    x = miglfdum2(:,15)/3600;
    slopeself = goflfprop.rmse./sqrt(sum((x-mean(x)).^2));
    interceptself = slopeself.*sqrt(sum(x.^2)./length(x));
    self(i, :) = [goflfprop.rmse slopeself interceptself];
    
    x = mighfdum2(:,15)/3600;
    y = mighfdum2(:,1);
    x_bar = wt_mean(x,wthf);
    y_bar = wt_mean(y,wthf);
    x_var = sum(wthf.*(x-x_bar).^2) / sum(wthf);
    y_var = sum(wthf.*(y-y_bar).^2) / sum(wthf);
    xy_cov = sum(wthf.*(x-x_bar).*(y-y_bar)) / sum(wthf);
    pearwthf(i) = xy_cov / sqrt(x_var*y_var);
    
    % compute the LF weights in robust linear regression, see NOTES above
    rlf = outplfprop.residuals;   % usual residuals
    x = [miglfdum2(:,15)/3600];
    hatmat = x*inv(x'*x)*x';
    h = zeros(size(hatmat,1),1);    % leverage of least square
    for jj = 1 : size(hatmat,1)
        h(jj) = hatmat(jj,jj);
    end    
    radj = rlf./sqrt(1-h);      % adjusted residuals
    K = 4.685;
    s = mad(rlf,1)/0.6745;
    u = radj/(K*s);
    wtlf = zeros(length(u),1);    % rubust weight of next iteration
    for jj = 1 : length(u)
        if abs(u(jj)) < 1
            wtlf(jj) = (1-(u(jj))^2)^2;
        else
            wtlf(jj) = 0;
        end
    end
    
    % get the vertical distance from LF detections to the fitted HF line
    fittemp = feval(fitobjhfprop,miglfdum2(:,15)/3600);
    verttemp = miglfdum2(:,1) - fittemp;     % can be + or -
    vertdistraw(1: length(verttemp), i) = verttemp;     % raw vertical distance
    verttempwt = wtlf.*verttemp;      % weighted vertical distance by multiplying the weight
    vertdistwt(1: length(verttempwt), i) = verttempwt;
    weight(1: length(verttempwt), i) =  wtlf;

%     % for mig 9, there are a few parts to throw away
%     if i == 9
%         ind9 = find(miglfdum2(:,1)>-10);
%     % for mig 4, only the early portion that passes the 002 region is used    
%     elseif i == 4
%         ind4 = find(miglfdum2(:,15)<5836);
%     end
        
    % get the information of propagation vector indicating range etc.
    [tmpx,tmpy] = coordinate_rot(medxhf,medyhf,-(angbest(i)-90),0,0);
    medallhf(i,1:4) = [medxhf medyhf tmpx tmpy];
    [tmpx,tmpy] = coordinate_rot(medxlf,medylf,-(angbest(i)-90),0,0);
    medalllf(i,1:4) = [medxlf medylf tmpx tmpy];
%     if i == 4
%         aaa = median(mighf(mighf(:,15)<5836,1));
%         bbb = median(mighf(mighf(:,15)<5836,2));
%         [tmpx,tmpy] = coordinate_rot(aaa,bbb,-(angbest(i)-90),0,0);
%         medallhf(i,1:4) = [aaa bbb tmpx tmpy];
%         aaa = median(miglf(miglf(:,15)<5836,1));
%         bbb = median(miglf(miglf(:,15)<5836,2));
%         [tmpx,tmpy] = coordinate_rot(aaa,bbb,-(angbest(i)-90),0,0);
%         medalllf(i,1:4) = [aaa bbb tmpx tmpy];
%     end
    
    ranvechf(i,1:2) = [min(mighfdum(:,1)) max(mighfdum(:,1))];
    ranveclf(i,1:2) = [min(miglfdum(:,1)) max(miglfdum(:,1))];
    
    ranvechf95(i,1:2) = [prctile(mighfdum(:,1),5) prctile(mighfdum(:,1),95)];
    ranveclf95(i,1:2) = [prctile(miglfdum(:,1),5) prctile(miglfdum(:,1),95)];
    
    ranvechf98(i,1:2) = [prctile(mighfdum(:,1),2) prctile(mighfdum(:,1),98)];
    ranveclf98(i,1:2) = [prctile(miglfdum(:,1),2) prctile(miglfdum(:,1),98)];
    
    nlt = sum(mighfdum(:,1)< ranvechf98(i,1));      % number that lower than the current range
    ngt = sum(mighfdum(:,1)> ranvechf98(i,2));      % number that greater than the current range
    
    mighfsort = sort(mighfdum(:,1));    % sort in descend order
    if nlt < 2
        ranvechf98(i,1) = mighfsort(3);     % so at least throw away 2 abnormals
    end
    
    if ngt < 2
        ranvechf98(i,2) = mighfsort(end-2);     % so at least throw away 2 abnormals
    end
    
%     if i == 4
%         ranvechf98(i,1:2) = [-18.3 -3.5];
%         ranveclf98(i,1:2) = [-21.4 -3.5];
%     end

    ranvechf99(i,1:2) = [prctile(mighfdum(:,1),1) prctile(mighfdum(:,1),99)];
    ranveclf99(i,1:2) = [prctile(miglfdum(:,1),1) prctile(miglfdum(:,1),99)];
                
    % use the median location and angbest to express the line denoting the propagation
    y2 = medallhf(i,2);
    x2 = medallhf(i,1);
    a2 = 1/tan(deg2rad(angbest(i)));
    b2 = y2-a2*x2;
    y = linefcn(x,a2,b2);
    % get the crossing point of this RTM and the division line
    [x0,y0] = linecrossing(a1,b1,a2,b2);
    
    % rotate this point to get the division along the propagation
    [propdiv,~] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
    tmp = find(fittemp < propdiv);
    indltdiv(1: length(tmp), i) = tmp;  
    tmp = find(fittemp >= propdiv);
    indgtdiv(1: length(tmp), i) = tmp;  
    
    fittemp = feval(fitobjhfprop,miglfdum2(:,15)/3600);
    verttemp = miglfdum2(:,1) - fittemp;     % can be + or -
    vertdistraw(1: length(verttemp), i) = verttemp;     % raw vertical distance
    verttempwt = wtlf.*verttemp;      % weighted vertical distance by multiplying the weight
    vertdistwt(1: length(verttempwt), i) = verttempwt;
    weight(1: length(verttempwt), i) =  wtlf;
    
    %%% find the independent LF detections, i.e. non-overlapping time window
    [miglfindpen,indindpen] = find_independent_detection(miglfdum2, 0);
    inumlf(i) = length(indindpen);
    indlfindpen(1: inumlf(i), i) =  indindpen;
    
    text(f.ax(3),0.4,0.2,sprintf('Slope: %.1f km/h',slopeprophf(i)),'FontSize',8,...
         'unit','normalized','horizontalalignment','left');
    text(f.ax(3),0.94,0.2,sprintf('SE: %.2f',slopesehf),'FontSize',8,...
         'unit','normalized','horizontalalignment','right'); 
    text(f.ax(3),0.4,0.13,sprintf('Pearson: %.3f',pearwthf(i)),'FontSize',8,...
         'unit','normalized','horizontalalignment','left');
    text(f.ax(3),0.94,0.12,sprintf('SE_{LF}: %.2f',slopeself),'FontSize',8,...
         'unit','normalized','horizontalalignment','right');  
    text(f.ax(3),0.4,0.05,sprintf('LF-HF: %.2f km',offset(i)),'FontSize',10,'unit','normalized'); 
    hold(f.ax(3),'off');    
    
    % subplot 4 of figure i
    hold(f.ax(4),'on');
    f.ax(4).FontSize = 9;
    f.ax(4).Box = 'on';
    yyaxis(f.ax(4),'left');
    plot(f.ax(4),angle,rmsehf(i,:),'o-','linew',1,'color','k','markers',4);
    if angbest(i) == angrmsehf(i)
        plot(f.ax(4),[angrmsehf(i) angrmsehf(i)],f.ax(4).YLim,'--','linew',1.5,'color','k');
    else
        plot(f.ax(4),[angrmsehf(i) angrmsehf(i)],f.ax(4).YLim,'--','linew',0.5,'color','k');
    end
    f.ax(4).YColor = 'k';
    ylabel(f.ax(4),'RMSE of HF','fontsize',11);
    yyaxis(f.ax(4),'right');
    plot(f.ax(4),angle,slopehf(i,:),'^:','linew',1,'color',[0.12 0.56 1],'markers',4);
    if angbest(i) == angslopehf(i)
        plot(f.ax(4),[angslopehf(i) angslopehf(i)],f.ax(4).YLim,'--','linew',1.5,...
             'color',[0.12 0.56 1]);    %[0 0 0.5];
    else
        plot(f.ax(4),[angslopehf(i) angslopehf(i)],f.ax(4).YLim,'--','linew',0.5,...
             'color',[0.12 0.56 1]);
    end    
    f.ax(4).YColor = [0.12 0.56 1]; %[0 0 0.5];
    ylabel(f.ax(4),'Slope of HF (km/h)','fontsize',11);
    
    xlim(f.ax(4),[0 360]);
%     ymax = f.ax(4).YLim(2)+0.1;
%     ylim(f.ax(4),[0 ymax]);
%     ylim(f.ax(4),[0 1]);
    xlabel(f.ax(4),'Trial propagation direction (deg)','fontsize',11);
    hold(f.ax(4),'off');
    
%     % subplot 4 of figure i
%     hold(f.ax(4),'on');
%     f.ax(4).FontSize = 9;
    numhf(i) = length(mighfdum2(:,1));
    numlf(i) = length(miglfdum2(:,1));
    resprophf(i,1:numhf(i)) = outphfprop.residuals;
    resproplf(i,1:numlf(i)) = outplfprop.residuals+(intcptproplf(i)-intcptprophf(i));
%     [xlochf, a, b, pdfvalhf] = weighted_bar_pdf(resprophf(i,1:numhf(i)), wthf, 0.5, 'int');
%     hfHdl = bar(f.ax(4),xlochf, pdfvalhf,1,'stacked');
%     hfHdl(1).FaceColor = [0.6 0.6 0.6];
%     [xloclf, ~, ~, pdfvallf] = weighted_bar_pdf(resproplf(i,1:numlf(i)), wtlf, 0.5, 'int');
%     lfHdl = bar(f.ax(4),xloclf, pdfvallf,1,'stacked','facea',0.6);
%     lfHdl(1).FaceColor = [0.6 1 1];


% returns estimates of normal distribution parameters, mu and sigma
% under 95% confidence interval
    [muhf,sigmahf,muhfCI,sigmahfCI] = normfit((resprophf(i,1:numhf(i)))', 0.05, [], wthf);
    [mulf,sigmalf,mulfCI,sigmalfCI] = normfit((resproplf(i,1:numlf(i)))', 0.05, [], wtlf);
%     text(f.ax(4),0.04,0.91,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
%     text(f.ax(4),0.4,0.92,sprintf('LF-HF = %.2f km',offset(i)),'FontSize',11,'unit','normalized');
%     f.ax(4).Box = 'on';
%     grid(f.ax(4), 'on');
%     f.ax(4).GridLineStyle = '--';
%     ymax = f.ax(4).YLim(2)+0.1;
%     ylim(f.ax(4),[0 ymax]);
% %     ylim(f.ax(4),[0 1]);
%     xlabel(f.ax(4),'Residual in prop. (km)','fontsize',11);
%     ylabel(f.ax(4),'PDF estimate','fontsize',11);
%     hold(f.ax(4),'off');
    
    
    %%% save figure
    print(f.fig,'-dpdf',strcat(rstpath,'/LZB.',num2str(nfam),'autortmv2.mig.proj.lfit',...
          num2str(trange(i,1)),'_',num2str(trange(i,2)),'-',num2str(trange(i,3)),...
          '.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),'.',...
          num2str(ccminlf),'.pdf'));
%     print(f.fig,'-dpdf',strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2020/Figures/LZB',...
%           num2str(i),'.pdf'));
    
    % get the percentage of detections from specific fam
    for j = 1: nfam
        famperchf(i,j) = sum(mighfdum2(:,end) == str2double(nfampool(j,:))) / size(mighfdum2,1);
        famperclf(i,j) = sum(miglfdum2(:,end) == str2double(nfampool(j,:))) / size(miglfdum2,1);
    end
        
end
[tmp,tmpind] = sort(sehf(:,2),'descend');
sortsehf = [tmpind tmp];
[tmp,tmpind] = sort(pearwthf(:),'ascend');
sortpearwthf = [tmpind tmp];


%% plot the migration direction in map
%%% define and position the figure frame and axes of each plot
f123.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f123.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 1;
for isub = 1:nrow*ncol
    f123.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
% set(f123.ax(1), 'position', [ 0.1, 0.1, 0.35, 0.75]);
% set(f123.ax(2), 'position', [ 0.55, 0.1, 0.35, 0.75]);

xran = [-16 24];
yran = [-20 20];

% subplot 1
ax = f123.ax(1);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';

for i = 1: size(trange,1)
    [rotx1, roty1] = complex_rot(0, ranvechf98(i,2)-medallhf(i,3), -angbest(i));
    [rotx2, roty2] = complex_rot(0, medallhf(i,3)-ranvechf98(i,1), -angbest(i));
    xvect = [medallhf(i,1)-rotx2 medallhf(i,1)+rotx1];
    yvect = [medallhf(i,2)-roty2 medallhf(i,2)+roty1];
    drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1,'linestyle','-','color',[0.7 0.7 0.7]);
    
end

for i = 1: size(trange,1)
    scatter(ax,medallhf(i,1),medallhf(i,2), 30, i, 'filled','o');  %, 'MarkerEdgeColor', 'w')
end
colormap(ax,'jet');
c=colorbar(ax,'SouthOutside');
% pos = ax.Position;
% c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
%     juldate = num2str(trange(i,1));
%     yr = str2double(juldate(1:4));
%     date = str2double(juldate(5:end));
%     a = jul2dat(yr,date);
%     mo = a(1);
%     if mo == 9
%         mo = {' Sep. '};
%     elseif mo == 7
%         mo = {' Jul. '};
%     else
%         mo = {' Mar. '};
%     end
%     day = num2str(a(2));
%     yr = num2str(a(3));
%     c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
%     c.Label.FontSize = 11;
caxis(ax,[1 size(trange,1)]);
text(ax,0.85,0.1,'HF','FontSize',12,'unit','normalized');
text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.85,0.93,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

% % subplot 2
% ax = f123.ax(2);
% hold(ax,'on');
% plot(ax,[-100 100],[0 0],'k--');
% plot(ax,[0 0],[-100 100],'k--');
% ax.FontSize = 9;
% 
% for i = 1: size(trange,1)
%     [rotx1, roty1] = complex_rot(0, ranveclf98(i,2)-medalllf(i,3), -angbest(i));
%     [rotx2, roty2] = complex_rot(0, medalllf(i,3)-ranveclf98(i,1), -angbest(i));
%     xvect = [medalllf(i,1)-rotx2 medalllf(i,1)+rotx1];
%     yvect = [medalllf(i,2)-roty2 medalllf(i,2)+roty1];
%     drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1,'linestyle','-','color',[0.7 0.7 0.7]);
%     
% end
% 
% for i = 1: size(trange,1)
%     scatter(ax,medalllf(i,1),medalllf(i,2), 30, i, 'filled','o');  %, 'MarkerEdgeColor', 'w')
% end
% colormap(ax,'jet');
% c=colorbar(ax,'SouthOutside');
% % pos = ax.Position;
% % c.Position = [pos(1), pos(2)-0.03, pos(3), 0.02];
% %     c.TickLabels=[];
% %     juldate = num2str(trange(i,1));
% %     yr = str2double(juldate(1:4));
% %     date = str2double(juldate(5:end));
% %     a = jul2dat(yr,date);
% %     mo = a(1);
% %     if mo == 9
% %         mo = {' Sep. '};
% %     elseif mo == 7
% %         mo = {' Jul. '};
% %     else
% %         mo = {' Mar. '};
% %     end
% %     day = num2str(a(2));
% %     yr = num2str(a(3));
% %     c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
% %     c.Label.FontSize = 11;
% caxis(ax,[1 size(trange,1)]);
% text(ax,0.85,0.1,'LF','FontSize',12,'unit','normalized');
% text(ax,0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% text(ax,0.85,0.93,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% ax.Box = 'on';
% grid(ax, 'on');
% axis(ax, 'equal');
% ax.GridLineStyle = '--';
% ax.XAxisLocation = 'top';
% xticks(ax,xran(1):5:xran(2));
% yticks(ax,yran(1):5:yran(2));
% xlabel(ax,'E (km)','fontsize',11);
% ylabel(ax,'N (km)','fontsize',11);
% hold(ax,'off');

% for i = 1: size(trange,1)
%     if angbest(i)>0 && angbest(i)<=90
%         disp(i)
%     end
% end

%%% save figure
print(f123.fig,'-dpdf',strcat(rstpath,'/autortmv2.prop_direction_map',num2str(nfam),'.pdf'));
% print(f.fig,'-dpdf',strcat(rstpath,'/unnamed.',num2str(winlenhf),'_',...
%     num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));
    
%% plot the migration direction in 4 quadrants in map
%%% define and position the figure frame and axes of each plot
f123.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f123.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f123.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f123.ax(1), 'position', [ 0.1, 0.55, 0.4, 0.4]);
set(f123.ax(2), 'position', [ 0.55, 0.55, 0.4, 0.4]);
set(f123.ax(3), 'position', [ 0.1, 0.1, 0.4, 0.4]);
set(f123.ax(4), 'position', [ 0.55, 0.1, 0.4, 0.4]);

xran = [-16 24];
yran = [-20 20];

% subplot 1
ax = f123.ax(1);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';

% indplt = find(angbest>0 & angbest<=90);
% indplt = find(angbest>=40 & angbest<=100);
% same as setdiff(find(angbest>=40 & angbest<=100),[27,29])
indplt = [2,11,12,17,19,22,33,39,45,52,55];
largediffind = [2,3,12,15,25,27,32,34,36,38,41,46,49,50,55,56,57,62];   % large offset diff
poorlfind = setdiff(largediffind, [15,38,41,49,50,57,62]);  % offset diff is caused by poor LF SE
indplt = setdiff(indplt, poorlfind);    %
for i = 1: length(indplt)
    [rotx1, roty1] = complex_rot(0, ranvechf98(indplt(i),2)-medallhf(indplt(i),3), -angbest(indplt(i)));
    [rotx2, roty2] = complex_rot(0, medallhf(indplt(i),3)-ranvechf98(indplt(i),1), -angbest(indplt(i)));
    xvect = [medallhf(indplt(i),1)-rotx2 medallhf(indplt(i),1)+rotx1];
    yvect = [medallhf(indplt(i),2)-roty2 medallhf(indplt(i),2)+roty1];
    drawArrow(ax,xvect,yvect,xran,yran,'linewidth',0.5,'linestyle','-','color',[0.7 0.7 0.7],...
              'HeadLength', 5, 'HeadWidth', 3);
    
end

for i = 1: length(indplt)
    scatter(ax,medallhf(indplt(i),1),medallhf(indplt(i),2), 20, indplt(i), 'filled','o');  %, 'MarkerEdgeColor', 'w')
end
colormap(ax,'jet');
c=colorbar(ax,'SouthOutside');
c.Label.String = strcat('RTM number in chronological order');
caxis(ax,[1 size(trange,1)]);
text(ax,0.85,0.1,'NE','FontSize',12,'unit','normalized');
text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.85,0.93,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

% subplot 2
ax = f123.ax(2);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';

indplt = find(angbest>90 & angbest<=180);
for i = 1: length(indplt)
    [rotx1, roty1] = complex_rot(0, ranvechf98(indplt(i),2)-medallhf(indplt(i),3), -angbest(indplt(i)));
    [rotx2, roty2] = complex_rot(0, medallhf(indplt(i),3)-ranvechf98(indplt(i),1), -angbest(indplt(i)));
    xvect = [medallhf(indplt(i),1)-rotx2 medallhf(indplt(i),1)+rotx1];
    yvect = [medallhf(indplt(i),2)-roty2 medallhf(indplt(i),2)+roty1];
    drawArrow(ax,xvect,yvect,xran,yran,'linewidth',0.5,'linestyle','-','color',[0.7 0.7 0.7],...
              'HeadLength', 5, 'HeadWidth', 3);
    
end

for i = 1: length(indplt)
    scatter(ax,medallhf(indplt(i),1),medallhf(indplt(i),2), 20, indplt(i), 'filled','o');  %, 'MarkerEdgeColor', 'w')
end
colormap(ax,'jet');
c=colorbar(ax,'SouthOutside');
c.Label.String = strcat('RTM number in chronological order');
caxis(ax,[1 size(trange,1)]);
text(ax,0.85,0.1,'SE','FontSize',12,'unit','normalized');
text(ax,0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.85,0.93,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

% subplot 3
ax = f123.ax(3);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';

% indplt = find(angbest>180 & angbest<=270);
% indplt = find(angbest>=220 & angbest<=280);
% same as setdiff(find(angbest>=220 & angbest<=280),[1,31])
indplt = [3,6,7,8,10,13,14,15,16,23,24,25,30,32,37,40,42,43,44,47,49,50,57,60,61,62];
indplt = setdiff(indplt, poorlfind);
for i = 1: length(indplt)
    [rotx1, roty1] = complex_rot(0, ranvechf98(indplt(i),2)-medallhf(indplt(i),3), -angbest(indplt(i)));
    [rotx2, roty2] = complex_rot(0, medallhf(indplt(i),3)-ranvechf98(indplt(i),1), -angbest(indplt(i)));
    xvect = [medallhf(indplt(i),1)-rotx2 medallhf(indplt(i),1)+rotx1];
    yvect = [medallhf(indplt(i),2)-roty2 medallhf(indplt(i),2)+roty1];
    drawArrow(ax,xvect,yvect,xran,yran,'linewidth',0.5,'linestyle','-','color',[0.7 0.7 0.7],...
              'HeadLength', 5, 'HeadWidth', 3);
    
end

for i = 1: length(indplt)
    scatter(ax,medallhf(indplt(i),1),medallhf(indplt(i),2), 20, indplt(i), 'filled','o');  %, 'MarkerEdgeColor', 'w')
end
colormap(ax,'jet');
c=colorbar(ax,'SouthOutside');
c.Label.String = strcat('RTM number in chronological order');
caxis(ax,[1 size(trange,1)]);
text(ax,0.85,0.1,'SW','FontSize',12,'unit','normalized');
text(ax,0.04,0.93,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.85,0.93,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

% subplot 4
ax = f123.ax(4);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';

indplt = find(angbest>270 & angbest<=360);
for i = 1: length(indplt)
    [rotx1, roty1] = complex_rot(0, ranvechf98(indplt(i),2)-medallhf(indplt(i),3), -angbest(indplt(i)));
    [rotx2, roty2] = complex_rot(0, medallhf(indplt(i),3)-ranvechf98(indplt(i),1), -angbest(indplt(i)));
    xvect = [medallhf(indplt(i),1)-rotx2 medallhf(indplt(i),1)+rotx1];
    yvect = [medallhf(indplt(i),2)-roty2 medallhf(indplt(i),2)+roty1];
    drawArrow(ax,xvect,yvect,xran,yran,'linewidth',0.5,'linestyle','-','color',[0.7 0.7 0.7],...
              'HeadLength', 5, 'HeadWidth', 3);
    
end

for i = 1: length(indplt)
    scatter(ax,medallhf(indplt(i),1),medallhf(indplt(i),2), 20, indplt(i), 'filled','o');  %, 'MarkerEdgeColor', 'w')
end
colormap(ax,'jet');
c=colorbar(ax,'SouthOutside');
c.Label.String = strcat('RTM number in chronological order');
caxis(ax,[1 size(trange,1)]);
text(ax,0.85,0.1,'NW','FontSize',12,'unit','normalized');
text(ax,0.04,0.93,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.85,0.93,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

%%% save figure
print(f123.fig,'-dpdf',strcat(rstpath,'/autortmv2.prop_direction_4q_map',num2str(nfam),'.pdf'));
% print(f.fig,'-dpdf',strcat(rstpath,'/unnamed.',num2str(winlenhf),'_',...
%     num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));


%% obtain the vertical distance from lf points to hf line in certain regions

% NOTE here, 2020/09/22, u can choose which vertical distance to use
% vertdist = vertdistwt;
vertdist = vertdistraw;

% SW migs
indsw = [3,6,7,8,10,13,14,15,16,23,24,25,30,32,37,40,42,43,44,47,49,50,57,60,61,62];    % index of SW migrations
% indsw = setdiff(indsw, poorlfind);
indsw = setdiff(indsw, largediffind);
aaa = find(self(:,2)>=3);
indsw = setdiff(indsw, aaa);
% indsw = setdiff(indsw,largediffind);

vdistsw = [];
wtsw = [];
ivdistsw = [];  % dist and weight of independent detections
iwtsw = [];
% for detections that are greater than division
vdistswgt = [];
wtswgt = [];
ivdistswgt = [];  
iwtswgt = [];
% for detections that are lower than division
vdistswlt = [];
wtswlt = [];
ivdistswlt = [];  
iwtswlt = [];

for i = 1: length(indsw)
    iind = indlfindpen(1:inumlf(indsw(i)), indsw(i));  % index of independent detection in that migration indsw(i)
    indlt = indltdiv(~isnan(indltdiv(:,indsw(i))), indsw(i));   
    indgt = indgtdiv(~isnan(indgtdiv(:,indsw(i))), indsw(i));
%     if indsw(i) == 9
%         temp1 = vertdist(1:numlf(indsw(i)), indsw(i));
%         temp2 = temp1(ind9);
%         vdistsw = [vdistsw; temp2];
%         temp2 = temp1(intersect(indgt, ind9));
%         vdistswgt = [vdistswgt; temp2];
%         temp2 = temp1(intersect(indlt, ind9));
%         vdistswlt = [vdistswlt; temp2];
%         temp2 = temp1(intersect(iind, ind9));
%         ivdistsw = [ivdistsw; temp2];
%         temp2 = temp1(intersect(intersect(iind, ind9), indgt));
%         ivdistswgt = [ivdistswgt; temp2]; 
%         temp2 = temp1(intersect(intersect(iind, ind9), indlt));
%         ivdistswlt = [ivdistswlt; temp2];
%         
%         temp1 = weight(1:numlf(indsw(i)), indsw(i));
%         temp2 = temp1(ind9);
%         wtsw = [wtsw; temp2];
%         temp2 = temp1(intersect(indgt, ind9));
%         wtswgt = [wtswgt; temp2];
%         temp2 = temp1(intersect(indlt, ind9));
%         wtswlt = [wtswlt; temp2];
%         temp2 = temp1(intersect(iind, ind9));
%         iwtsw = [iwtsw; temp2];
%         temp2 = temp1(intersect(intersect(iind, ind9), indgt));
%         iwtswgt = [iwtswgt; temp2]; 
%         temp2 = temp1(intersect(intersect(iind, ind9), indlt));
%         iwtswlt = [iwtswlt; temp2];
%     else
        temp1 = vertdist(1:numlf(indsw(i)), indsw(i));
        vdistsw = [vdistsw; temp1];
        vdistswgt = [vdistswgt; temp1(indgt)];
        vdistswlt = [vdistswlt; temp1(indlt)];
        ivdistsw = [ivdistsw; temp1(iind)];
        ivdistswgt = [ivdistswgt; temp1(intersect(iind, indgt))];
        ivdistswlt = [ivdistswlt; temp1(intersect(iind, indlt))];
        
        temp1 = weight(1:numlf(indsw(i)), indsw(i));
        wtsw = [wtsw; temp1];
        wtswgt = [wtswgt; temp1(indgt)];
        wtswlt = [wtswlt; temp1(indlt)];
        iwtsw = [iwtsw; temp1(iind)];
        iwtswgt = [iwtswgt; temp1(intersect(iind, indgt))];
        iwtswlt = [iwtswlt; temp1(intersect(iind, indlt))];
%     end    
end


% NE migs
indne = [2,11,12,17,19,22,33,39,45,52,55];    % index of NE migrations
% indne = setdiff(indne, poorlfind);
indne = setdiff(indne, largediffind);
aaa = find(self(:,2)>=3);
indne = setdiff(indne, aaa);
% indne = setdiff(indne,largediffind);

vdistne = [];
wtne = [];
ivdistne = [];  % dist and weight of independent detections
iwtne = [];
% for detections that are greater than division
vdistnegt = [];
wtnegt = [];
ivdistnegt = [];  
iwtnegt = [];
% for detections that are lower than division
vdistnelt = [];
wtnelt = [];
ivdistnelt = [];  
iwtnelt = [];
for i = 1: length(indne)
    iind = indlfindpen(1:inumlf(indne(i)), indne(i));  % index of independent detection in that migration indsw(i)
    indlt = indltdiv(~isnan(indltdiv(:,indne(i))), indne(i));   
    indgt = indgtdiv(~isnan(indgtdiv(:,indne(i))), indne(i));
    temp1 = vertdist(1:numlf(indne(i)), indne(i));
    vdistne = [vdistne; temp1];
    vdistnegt = [vdistnegt; temp1(indgt)];
    vdistnelt = [vdistnelt; temp1(indlt)];
    ivdistne = [ivdistne; temp1(iind)];
    ivdistnegt = [ivdistnegt; temp1(intersect(iind, indgt))];
    ivdistnelt = [ivdistnelt; temp1(intersect(iind, indlt))];
    
    temp1 = weight(1:numlf(indne(i)), indne(i));
    wtne = [wtne; temp1];
    wtnegt = [wtnegt; temp1(indgt)];
    wtnelt = [wtnelt; temp1(indlt)];
    iwtne = [iwtne; temp1(iind)];
    iwtnegt = [iwtnegt; temp1(intersect(iind, indgt))];
    iwtnelt = [iwtnelt; temp1(intersect(iind, indlt))];

end


% NE mig that passes through fam 002 region
% ind002ne = [22];
% ind002ne = [22,21];
ind002ne = [27,29];
vdist002ne = [];
wt002ne = [];
ivdist002ne = [];  % dist and weight of independent detections
iwt002ne = [];
% for detections that are greater than division
vdist002negt = [];
wt002negt = [];
ivdist002negt = [];  
iwt002negt = [];
% for detections that are lower than division
vdist002nelt = [];
wt002nelt = [];
ivdist002nelt = [];  
iwt002nelt = [];
for i = 1: length(ind002ne)
    iind = indlfindpen(1:inumlf(ind002ne(i)), ind002ne(i));  % index of independent detection in that migration indsw(i)
    indlt = indltdiv(~isnan(indltdiv(:,ind002ne(i))), ind002ne(i));
    indgt = indgtdiv(~isnan(indgtdiv(:,ind002ne(i))), ind002ne(i));
    temp1 = vertdist(1:numlf(ind002ne(i)), ind002ne(i));
    vdist002ne = [vdist002ne; temp1];
    vdist002negt = [vdist002negt; temp1(indgt)];
    vdist002nelt = [vdist002nelt; temp1(indlt)];
    ivdist002ne = [ivdist002ne; temp1(iind)];
    ivdist002negt = [ivdist002negt; temp1(intersect(iind, indgt))];
    ivdist002nelt = [ivdist002nelt; temp1(intersect(iind, indlt))];

    temp1 = weight(1:numlf(ind002ne(i)), ind002ne(i));
    wt002ne = [wt002ne; temp1];
    wt002negt = [wt002negt; temp1(indgt)];
    wt002nelt = [wt002nelt; temp1(indlt)];
    iwt002ne = [iwt002ne; temp1(iind)];
    iwt002negt = [iwt002negt; temp1(intersect(iind, indgt))];
    iwt002nelt = [iwt002nelt; temp1(intersect(iind, indlt))];
end


% SW mig that passes through fam 002 region
% ind002sw = [25];
% ind002sw = [25,4];
ind002sw = [31];
vdist002sw = [];
wt002sw = [];
ivdist002sw = [];  % dist and weight of independent detections
iwt002sw = [];
% for detections that are greater than division
vdist002swgt = [];
wt002swgt = [];
ivdist002swgt = [];  
iwt002swgt = [];
% for detections that are lower than division
vdist002swlt = [];
wt002swlt = [];
ivdist002swlt = [];  
iwt002swlt = [];
for i = 1: length(ind002sw)
    iind = indlfindpen(1:inumlf(ind002sw(i)), ind002sw(i));  % index of independent detection in that migration indsw(i)
    indlt = indltdiv(~isnan(indltdiv(:,ind002sw(i))), ind002sw(i));
    indgt = indgtdiv(~isnan(indgtdiv(:,ind002sw(i))), ind002sw(i));
%     if ind002sw(i) == 4
%         temp1 = vertdist(1:numlf(ind002sw(i)), ind002sw(i));
%         temp2 = temp1(ind4);
%         vdist002sw = [vdist002sw; temp2];
%         temp2 = temp1(intersect(indgt, ind4));
%         vdist002swgt = [vdist002swgt; temp2];
%         temp2 = temp1(intersect(indlt, ind4));
%         vdist002swlt = [vdist002swlt; temp2];
%         temp2 = temp1(intersect(iind, ind4));
%         ivdist002sw = [ivdist002sw; temp2];
%         temp2 = temp1(intersect(intersect(iind, ind4), indgt));
%         ivdist002swgt = [ivdist002swgt; temp2]; 
%         temp2 = temp1(intersect(intersect(iind, ind4), indlt));
%         ivdist002swlt = [ivdist002swlt; temp2];
%         
%         temp1 = weight(1:numlf(ind002sw(i)), ind002sw(i));
%         temp2 = temp1(ind4);
%         wt002sw = [wt002sw; temp2];
%         temp2 = temp1(intersect(indgt, ind4));
%         wt002swgt = [wt002swgt; temp2];
%         temp2 = temp1(intersect(indlt, ind4));
%         wt002swlt = [wt002swlt; temp2];
%         temp2 = temp1(intersect(iind, ind4));
%         iwt002sw = [iwt002sw; temp2];
%         temp2 = temp1(intersect(intersect(iind, ind4), indgt));
%         iwt002swgt = [iwt002swgt; temp2]; 
%         temp2 = temp1(intersect(intersect(iind, ind4), indlt));
%         iwt002swlt = [iwt002swlt; temp2];              
%         
%     else
        temp1 = vertdist(1:numlf(ind002sw(i)), ind002sw(i));
        vdist002sw = [vdist002sw; temp1];
        vdist002swgt = [vdist002swgt; temp1(indgt)];
        vdist002swlt = [vdist002swlt; temp1(indlt)];
        ivdist002sw = [ivdist002sw; temp1(iind)];
        ivdist002swgt = [ivdist002swgt; temp1(intersect(iind, indgt))];
        ivdist002swlt = [ivdist002swlt; temp1(intersect(iind, indlt))];
        
        temp1 = weight(1:numlf(ind002sw(i)), ind002sw(i));
        wt002sw = [wt002sw; temp1];
        wt002swgt = [wt002swgt; temp1(indgt)];
        wt002swlt = [wt002swlt; temp1(indlt)];
        iwt002sw = [iwt002sw; temp1(iind)];
        iwt002swgt = [iwt002swgt; temp1(intersect(iind, indgt))];
        iwt002swlt = [iwt002swlt; temp1(intersect(iind, indlt))];                
%     end
end

%% plot the percentage of detections from specific fams for selected RTMs
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 6;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 4;
ncol = 1;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

famcheck = ['099'; '017'; '043'; '068'];
for i = 1: size(famcheck,1)
    [~,ind] = ismember(famcheck(i,:),nfampool,'rows');
    indplt = [indsw indne]; 
    perc = famperclf(indplt, ind);
    
    ax = f.ax(i);
    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');
    plot(ax, perc, 'bo-','markersize',4);
    text(ax, 0.02, 0.9, famcheck(i,:),'unit','normalized');
    for j = 1: length(indplt) 
        text(ax, j, perc(j)+0.05, num2str(indplt(j)));
    end
    ax.YLim = [0 0.5];
    ax.YTick = 0: 0.1: 0.5;
    plot(ax,[12.5 12.5],ax.YLim,'k--','linew',1.5);
    hold(ax, 'off');
end
print(f.fig,'-dpdf',strcat(rstpath,'/autortmv2.fam_perc_dist.pdf'));


f.fig=figure;
widin = 5;  % maximum width allowed is 8.5 inches
htin = 10;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 4;
ncol = 1;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

famcheck = ['099'; '017'; '043'; '068'];
for i = 1: size(famcheck,1)
    [~,ind] = ismember(famcheck(i,:),nfampool,'rows');
    indplt = [indsw indne]; 
    perc = famperclf(indplt, ind);
    offs = offset(indplt);
    offs = offs';
    
    [rho1,pval] = corr(perc, offs,'Type','Pearson');
    [rho2,pval] = corr(perc, offs,'Type','Spearman');
    [rho3,pval] = corr(perc, offs,'Type','Kendall');
    
    ax = f.ax(i);
    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');
%     scatter(ax, perc, offs,12,'b');
    scatter(ax, famperclf(indsw, ind), -offset(indsw),12,'r','filled');
    scatter(ax, famperclf(indne, ind), offset(indne),12,'b','filled');
    text(ax, 0.02, 0.1, famcheck(i,:),'unit','normalized');
    for j = 1: length(indplt) 
        text(ax, perc(j), offs(j)-0.2, num2str(indplt(j)), 'fontsize',8);
    end
    text(ax,0.7, 0.3, strcat({'Pearson: '},sprintf('%.2f',rho1)),'unit','normalized');
    text(ax,0.7, 0.2, strcat({'Spearman: '},sprintf('%.2f',rho2)),'unit','normalized');
    text(ax,0.7, 0.1, strcat({'Kendall: '},sprintf('%.2f',rho3)),'unit','normalized');
    ax.XLim = [0 0.5];
    ax.XTick = 0: 0.1: 0.5;
%     plot(ax,[12.5 12.5],ax.YLim,'k--','linew',1.5);
    xlabel(ax,'Percentage of detections');
    ylabel(ax,'LF lags or leads HF (km)');
    hold(ax, 'off');
end
print(f.fig,'-dpdf',strcat(rstpath,'/autortmv2.fam_perc_dist_offset.pdf'));


%% Change sign of distance to geographical direction
%%% NOTE, 2020/10/27, recall that this 'vertical' distance is LF-HF in the propagation direction,
%%% in which positive means LF is leading, negative is lagging. But this is limited to the
%%% propagation. So if change to absolute geographycal direction, say WSW-ENW, if look at wsw
%%% propagating migrations, if lf leads, then distance should be +, but you need to flip the sign to
%%% convert to geo direction, since lf is in fact at the wsw of hf, - in geo direction. Same if lf
%%% lags. On the other hand, if look at ene group, if lf lags, the raw distance should be -, in
%%% fact, it means lf is at the wsw of hf, still - in the chosen geo coordinates, thus sign of this
%%% group remains unchanged.
%%%
%%% even divided spatially into two parts, for the same propagationn, sign of the offset is the same 

%%% currently focus on main LZB region and 002 region only

% main LZB region 
vdistsw = -vdistsw;     % flip the sign for SW group, check note above
ivdistsw = -ivdistsw;
vdistswgt = -vdistswgt;     % flip the sign for SW group, check note above
ivdistswgt = -ivdistswgt;
vdistswlt = -vdistswlt;     % flip the sign for SW group, check note above
ivdistswlt = -ivdistswlt;
vdistne = vdistne;     % keep the sign for NE group, check note above
ivdistne = ivdistne;
vdistnegt = vdistnegt;     % keep the sign for NE group, check note above
ivdistnegt = ivdistnegt;
vdistnelt = vdistnelt;     % keep the sign for NE group, check note above
ivdistnelt = ivdistnelt;

% 002 region
vdist002sw = -vdist002sw;     % flip the sign for SW group, check note above
ivdist002sw = -ivdist002sw;
vdist002swgt = -vdist002swgt;     % flip the sign for SW group, check note above
ivdist002swgt = -ivdist002swgt;
vdist002swlt = -vdist002swlt;     % flip the sign for SW group, check note above
ivdist002swlt = -ivdist002swlt;
vdist002ne = vdist002ne;     % keep the sign for NE group, check note above
ivdist002ne = ivdist002ne;
vdist002negt = vdist002negt;     % keep the sign for NE group, check note above
ivdist002negt = ivdist002negt;
vdist002nelt = vdist002nelt;     % keep the sign for NE group, check note above
ivdist002nelt = ivdist002nelt;


%% test the NULL hypothesis that SW and NE group has the same mean, without assuming same variance
%%% NOTE:
%%% test if the possibility of this null hypothesis is higher than a significance level (alpha). If
%%% so, it means the sample mean difference is NOT significant, then fail to reject this null 
%%% hypothesis. Otherwise, if p-value is smaller than alpha, then reject it and regard the
%%% difference is significant.

%%% 1. for SW and NE group at main LZB region, i.e. the NW part of the study region
% a. use significance level 0.05; frequency weights (2), and unequal variance (0)
[testmaina,probmaina,statsmaina] = ttest2_wt(vdistsw, wtsw, vdistne, wtne, 0.05, 2, 0)

% b. use significance level 0.05; reliability weights (1), and unequal variance (0)
[testmainb,probmainb,statsmainb] = ttest2_wt(vdistsw, wtsw, vdistne, wtne, 0.05, 1, 0)

% c. try equal variance (1) as well, since the variances in fact are very close 
[testmainc,probmainc,statsmainc] = ttest2_wt(vdistsw, wtsw, vdistne, wtne, 0.05, 2, 1)


%%% 2. for SW and NE group at fam 002 region
% a. use significance level 0.05; frequency weights (2), and unequal variance (0)
[test002a,prob002a,stats002a] = ttest2_wt(vdist002sw, wt002sw, vdist002ne, wt002ne, 0.05, 2, 0)

% b. use significance level 0.05; reliability weights (1), and unequal variance (0)
[test002b,prob002b,stats002b] = ttest2_wt(vdist002sw, wt002sw, vdist002ne, wt002ne, 0.05, 1, 0)

% c. try equal variance (1) as well, since the variances in fact are very close 
[test002c,prob002c,stats002c] = ttest2_wt(vdist002sw, wt002sw, vdist002ne, wt002ne, 0.05, 2, 1)


%% plot the combined migrations, unshifted all
%%% 1. for SW and NE group at main LZB region, i.e. the NW part of the study region
% indplt = union(indsw,indne);    % index to plot
xran = [-15 5];
yran = [-15 15];
binw = 0.5;
conf = 99;
[f,barsw,pdfxlocsw,pdfvalsw,barne,pdfxlocne,pdfvalne,lgd] = ...
    plt_hist_combined_RTMs(indne,indsw,ranvechf98,medallhf,angbest,xran,yran,vdistsw,wtsw,...
                           vdistne,wtne,binw,conf,'int');
text(f.ax(1),0.5,0.06,'All detections','FontSize',12,'unit','normalized',...
     'horizontalalignment','center');
text(f.ax(1),0.5,0.12,'Main LZB region,','FontSize',12,'unit','normalized',...
     'horizontalalignment','center');
xlim(f.ax(2), [-12.5 12.5]); 
xvect = [-7 -12.5+0.5];
yvect = [0.015 0.015];
drawArrow(f.ax(2),xvect,yvect,f.ax(2).XLim,f.ax(2).YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
% arrow1 = annotation('arrow','linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);
% arrow1.Parent = f.ax(2);
% arrow1.Position = [-7, 0.015, -5, 0] ;
xvect = [7 12.5-0.5];
yvect = [0.015 0.015];
drawArrow(f.ax(2),xvect,yvect,f.ax(2).XLim,f.ax(2).YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
text(f.ax(2),0.05,0.1,strcat({'WSW'}),'fontsize',11,'unit','normalized');
text(f.ax(2),0.95,0.1,strcat({'ENE'}),'fontsize',11,'unit','normalized',...
     'horizontalalignment','right');
lgd.String = [{'ENE propagation'},{'WSW propagation'}];
text(f.ax(3),0.95,0.86,strcat({'10X'}),'fontsize',10,'unit','normalized','EdgeColor','k',...
     'Margin',2,'horizontalalignment','right');
text(f.ax(3),0.05,0.15,strcat({'99% CI'}),'fontsize',12,'unit','normalized','horizontalalignment','left'); 

xlimit = f.ax(2).XLim/10;
xlim(f.ax(3),xlimit); 
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(3),xlimit(1),0.2);
[xe,ye] = ds2nfu(f.ax(2),xlimit(1),0.3);
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(3),xlimit(2),0.2);
[xe,ye] = ds2nfu(f.ax(2),xlimit(2),0.3);
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/autortmv2.vdist_mainSW+NE_all',num2str(nfam),'.pdf'));


%%% 2. for SW and NE group at fam 002 region
% indplt = union(ind002sw,ind002ne);    % index to plot
xran = [2 22];
yran = [-20 10];
[f,barsw002,pdfxloc002sw,pdfval002sw,barne002,pdfxloc002ne,pdfval002ne,lgd] = ...
    plt_hist_combined_RTMs(ind002ne,ind002sw,ranvechf98,medallhf,angbest,xran,yran,...
                           vdist002sw,wt002sw,vdist002ne,wt002ne,binw,conf,'int');
text(f.ax(1),0.5,0.06,'All detections','FontSize',12,'unit','normalized',...
     'horizontalalignment','center');
text(f.ax(1),0.5,0.12,'Family 002 region,','FontSize',12,'unit','normalized',...
     'horizontalalignment','center');
text(f.ax(2),0.05,0.1,strcat({'SW'}),'fontsize',11,'unit','normalized');
text(f.ax(2),0.95,0.1,strcat({'NE'}),'fontsize',11,'unit','normalized',...
     'horizontalalignment','right');
lgd.String = [{'NE propagation'},{'SW propagation'}];
% legend(f.ax(2),[barne002,barsw002],{'NE propagation','SW propagation'},'fontsize',7,...
%        'numcolumns',2,...
%        'Position',[0.55+0.35*0.048  c.Position(2)+ywid*2/3*0.93  0.35*0.9  ywid*2/3*0.05]); 
text(f.ax(3),0.95,0.86,strcat({'5X'}),'fontsize',10,'unit','normalized','EdgeColor','k',...
     'Margin',2,'horizontalalignment','right');
text(f.ax(3),0.95,0.15,strcat({'99% CI'}),'fontsize',12,'unit','normalized','horizontalalignment',...
     'right'); 
  
xlimit = f.ax(2).XLim/5;
xlim(f.ax(3),xlimit); 
[xs,ys] = ds2nfu(f.ax(3),xlimit(1),0.2);
[xe,ye] = ds2nfu(f.ax(2),xlimit(1),0.3);
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');
[xs,ys] = ds2nfu(f.ax(3),xlimit(2),0.2);
[xe,ye] = ds2nfu(f.ax(2),xlimit(2),0.3);
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/autortmv2.vdist_002SW+NE_all',num2str(nfam),'.pdf'));


% %% divide each group spatially into 2 parts, regroup and put NE or SW parts together
% indplt = union(indsw,indne);    % index to plot
% xran = [-15 5];
% yran = [-15 15];
% binw = 1;
% conf = 99;
% [f,barsw1,barne1,barsw2,barne2,lgd1,lgd2] = ...
%     plt_hist_combined_RTMs_v2(indne,indsw,ranvechf98,medallhf,angbest,xran,yran,binw,vdistswlt,wtswlt,...
%                            vdistnegt,wtnegt,vdistswgt,wtswgt,vdistnelt,wtnelt,conf);
% hold(f.ax(1),'on');
% xdiv = -7:0.01:0;
% a1 = 1/tan(deg2rad(180-22.5));
% b1 = y0-a1*x0;
% ydiv = linefcn(xdiv,a1,b1);
% plot(f.ax(1),xdiv,ydiv,'k--','linew',2);
% text(f.ax(1),0.55,0.06,'All detections','FontSize',12,'unit','normalized',...
%      'horizontalalignment','center');
% text(f.ax(1),0.55,0.12,'Main LZB region,','FontSize',12,'unit','normalized',...
%      'horizontalalignment','center');
%  
% text(f.ax(2),0.05,0.1,strcat({'WSW'}),'fontsize',10,'unit','normalized');
% text(f.ax(2),0.95,0.1,strcat({'ENE'}),'fontsize',10,'unit','normalized',...
%      'horizontalalignment','right');
% lgd1.String = [{'ENE propagation'},{'WSW propagation'}];
% 
% text(f.ax(3),0.95,0.86,strcat({'5X'}),'fontsize',10,'unit','normalized','EdgeColor','k',...
%      'Margin',2,'horizontalalignment','right');
% text(f.ax(3),0.03,0.1,strcat({'NE portion'}),'fontsize',11,'unit','normalized','horizontalalignment',...
%      'left'); 
% text(f.ax(3),0.03,0.3,strcat({'99% CI'}),'fontsize',11,'unit','normalized','horizontalalignment',...
%      'left');
% 
% text(f.ax(4),0.05,0.1,strcat({'WSW'}),'fontsize',10,'unit','normalized');
% text(f.ax(4),0.95,0.1,strcat({'ENE'}),'fontsize',10,'unit','normalized',...
%      'horizontalalignment','right');
% % lgd2.String = [{'ENE propagation'},{'WSW propagation'}];
% 
% text(f.ax(5),0.95,0.86,strcat({'5X'}),'fontsize',10,'unit','normalized','EdgeColor','k',...
%      'Margin',2,'horizontalalignment','right');
% text(f.ax(5),0.03,0.1,strcat({'SW portion'}),'fontsize',11,'unit','normalized','horizontalalignment',...
%      'left');  
% text(f.ax(5),0.03,0.3,strcat({'99% CI'}),'fontsize',11,'unit','normalized','horizontalalignment',...
%      'left'); 
% 
%  %%% save figure
% print(f.fig,'-dpdf',strcat(rstpath,'/vdist_unshifted_mainSW+NE_2portions_all',num2str(nfam),'.pdf'));
% 
% 
% %%% Also do the test to the 2 portions
% % 1. the NE portion
% % a. use significance level 0.05; frequency weights (2), and unequal variance (0)
% [testnea,probnea,statsnea] = ttest2_wt(vdistswlt,wtswlt,vdistnegt,wtnegt, 0.05, 2, 0)
% 
% % 2. the SW portion
% % a. use significance level 0.05; frequency weights (2), and unequal variance (0)
% [testswa,probswa,statsswa] = ttest2_wt(vdistswgt,wtswgt,vdistnelt,wtnelt, 0.05, 2, 0)


%% test the NULL hypothesis to those independent detections
%%% NOTE:
%%% test if the possibility of this null hypothesis is higher than a significance level (alpha). If
%%% so, it means the sample mean difference is NOT significant, then fail to reject this null 
%%% hypothesis. Otherwise, if p-value is smaller than alpha, then reject it and regard the
%%% difference is significant.

%%% 1. for SW and NE group at main LZB region, i.e. the NW part of the study region
% a. use significance level 0.05; frequency weights (2), and unequal variance (0)
[testmaina,probmaina,statsmaina] = ttest2_wt(ivdistsw, iwtsw, ivdistne, iwtne, 0.05, 2, 0)

% b. use significance level 0.05; reliability weights (1), and unequal variance (0)
[testmainb,probmainb,statsmainb] = ttest2_wt(ivdistsw, iwtsw, ivdistne, iwtne, 0.05, 1, 0)

% c. try equal variance (1) as well, since the variances in fact are very close
[testmainc,probmainc,statsmainc] = ttest2_wt(ivdistsw, iwtsw, ivdistne, iwtne, 0.05, 2, 1)

%%% 2. for SW and NE group at fam 002 region
% a. use significance level 0.05; frequency weights (2), and unequal variance (0)
[test002a,prob002a,stats002a] = ttest2_wt(ivdist002sw, iwt002sw, ivdist002ne, iwt002ne, 0.05, 2, 0)

% b. use significance level 0.05; reliability weights (1), and unequal variance (0)
[test002b,prob002b,stats002b] = ttest2_wt(ivdist002sw, iwt002sw, ivdist002ne, iwt002ne, 0.05, 1, 0)

% c. try equal variance (1) as well, since the variances in fact are very close
[test002c,prob002c,stats002c] = ttest2_wt(ivdist002sw, iwt002sw, ivdist002ne, iwt002ne, 0.05, 2, 1)


%% plot the combined migrations, unshifted independent detections
%%% 1. for SW and NE group at main LZB region, i.e. the NW part of the study region
% indplt = union(indsw,indne);    % index to plot
xran = [-15 5];
yran = [-15 15];
binw = 0.5;
conf = 99;
[f,barisw,pdfxlocisw,pdfvalisw,barine,pdfxlocine,pdfvaline,lgd] = ...
    plt_hist_combined_RTMs(indne,indsw,ranvechf98,medallhf,angbest,xran,yran,ivdistsw,iwtsw,...
                           ivdistne,iwtne,binw,conf,'int');
text(f.ax(1),0.5,0.06,'Independent detections','FontSize',12,'unit','normalized',...
     'horizontalalignment','center');
text(f.ax(1),0.5,0.12,'Main LZB region,','FontSize',12,'unit','normalized',...
     'horizontalalignment','center'); 
xlim(f.ax(2), [-12.5 12.5]); 
xvect = [-7 -12.5+0.5];
yvect = [0.015 0.015];
drawArrow(f.ax(2),xvect,yvect,f.ax(2).XLim,f.ax(2).YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
% arrow1 = annotation('arrow','linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);
% arrow1.Parent = f.ax(2);
% arrow1.Position = [-7, 0.015, -5, 0] ;
xvect = [7 12.5-0.5];
yvect = [0.015 0.015];
drawArrow(f.ax(2),xvect,yvect,f.ax(2).XLim,f.ax(2).YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
text(f.ax(2),0.05,0.1,strcat({'WSW'}),'fontsize',11,'unit','normalized');
text(f.ax(2),0.95,0.1,strcat({'ENE'}),'fontsize',11,'unit','normalized',...
     'horizontalalignment','right');
lgd.String = [{'ENE propagation'},{'WSW propagation'}];
text(f.ax(3),0.95,0.86,strcat({'10X'}),'fontsize',10,'unit','normalized','EdgeColor','k',...
     'Margin',2,'horizontalalignment','right');
text(f.ax(3),0.05,0.15,strcat({'99% CI'}),'fontsize',12,'unit','normalized','horizontalalignment','left'); 

xlimit = f.ax(2).XLim/10;
xlim(f.ax(3),xlimit); 
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(3),xlimit(1),0.2);
[xe,ye] = ds2nfu(f.ax(2),xlimit(1),0.3);
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(3),xlimit(2),0.2);
[xe,ye] = ds2nfu(f.ax(2),xlimit(2),0.3);
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');

print(f.fig,'-dpdf',strcat(rstpath,'/autortmv2.vdist_mainSW+NE_ind',num2str(nfam),...
      '.pdf'));
  
  
%%% 2. for SW and NE group at fam 002 region
% indplt = union(ind002sw,ind002ne);    % index to plot
xran = [2 22];
yran = [-20 10];
[f,barisw002,pdfxloc002isw,pdfval002isw,barine002,pdfxloc002ine,pdfval002ine,lgd] = ...
    plt_hist_combined_RTMs(ind002ne,ind002sw,ranvechf98,medallhf,angbest,xran,yran,...
                           ivdist002sw,iwt002sw,ivdist002ne,iwt002ne,binw,conf,'int');
text(f.ax(1),0.5,0.06,'Independent detections','FontSize',12,'unit','normalized',...
     'horizontalalignment','center');
text(f.ax(1),0.5,0.12,'Family 002 region,','FontSize',12,'unit','normalized',...
     'horizontalalignment','center');
text(f.ax(2),0.05,0.1,strcat({'SW'}),'fontsize',11,'unit','normalized');
text(f.ax(2),0.95,0.1,strcat({'NE'}),'fontsize',11,'unit','normalized',...
     'horizontalalignment','right');
lgd.String = [{'NE propagation'},{'SW propagation'}];
% legend(f.ax(2),[barne002,barsw002],{'NE propagation','SW propagation'},'fontsize',7,...
%        'numcolumns',2,...
%        'Position',[0.55+0.35*0.048  c.Position(2)+ywid*2/3*0.93  0.35*0.9  ywid*2/3*0.05]); 
text(f.ax(3),0.95,0.86,strcat({'5X'}),'fontsize',10,'unit','normalized','EdgeColor','k',...
     'Margin',2,'horizontalalignment','right');
text(f.ax(3),0.95,0.15,strcat({'99% CI'}),'fontsize',12,'unit','normalized','horizontalalignment',...
     'right'); 
  
xlimit = f.ax(2).XLim/5;
xlim(f.ax(3),xlimit); 
[xs,ys] = ds2nfu(f.ax(3),xlimit(1),0.2);
[xe,ye] = ds2nfu(f.ax(2),xlimit(1),0.3);
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');
[xs,ys] = ds2nfu(f.ax(3),xlimit(2),0.2);
[xe,ye] = ds2nfu(f.ax(2),xlimit(2),0.3);
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');

print(f.fig,'-dpdf',strcat(rstpath,'/autortmv2.vdist_002SW+NE_ind',num2str(nfam),...
      '.pdf'));


% %% independent detections, regroup them, put NE or SW parts from 2 propagation groups together
% indplt = union(indsw,indne);    % index to plot
% xran = [-15 5];
% yran = [-15 15];
% binw = 1;
% conf = 99;
% [f,barsw1,barne1,barsw2,barne2,lgd1,lgd2] = ...
%     plt_hist_combined_RTMs_v2(indne,indsw,ranvechf98,medallhf,angbest,xran,yran,binw,ivdistswlt,iwtswlt,...
%                            ivdistnegt,iwtnegt,ivdistswgt,iwtswgt,ivdistnelt,iwtnelt,conf);
% hold(f.ax(1),'on');
% xdiv = -7:0.01:0;
% a1 = 1/tan(deg2rad(180-22.5));
% b1 = y0-a1*x0;
% ydiv = linefcn(xdiv,a1,b1);
% plot(f.ax(1),xdiv,ydiv,'k--','linew',2);
% text(f.ax(1),0.55,0.06,'Independent detections','FontSize',12,'unit','normalized',...
%      'horizontalalignment','center');
% text(f.ax(1),0.55,0.12,'Main LZB region,','FontSize',12,'unit','normalized',...
%      'horizontalalignment','center');
%  
% text(f.ax(2),0.05,0.1,strcat({'WSW'}),'fontsize',10,'unit','normalized');
% text(f.ax(2),0.95,0.1,strcat({'ENE'}),'fontsize',10,'unit','normalized',...
%      'horizontalalignment','right');
% lgd1.String = [{'ENE propagation'},{'WSW propagation'}];
% 
% text(f.ax(3),0.95,0.86,strcat({'5X'}),'fontsize',10,'unit','normalized','EdgeColor','k',...
%      'Margin',2,'horizontalalignment','right');
% text(f.ax(3),0.03,0.1,strcat({'NE portion'}),'fontsize',11,'unit','normalized','horizontalalignment',...
%      'left'); 
% text(f.ax(3),0.03,0.3,strcat({'99% CI'}),'fontsize',11,'unit','normalized','horizontalalignment',...
%      'left');
% 
% text(f.ax(4),0.05,0.1,strcat({'WSW'}),'fontsize',10,'unit','normalized');
% text(f.ax(4),0.95,0.1,strcat({'ENE'}),'fontsize',10,'unit','normalized',...
%      'horizontalalignment','right');
% % lgd2.String = [{'ENE propagation'},{'WSW propagation'}];
% 
% text(f.ax(5),0.95,0.86,strcat({'5X'}),'fontsize',10,'unit','normalized','EdgeColor','k',...
%      'Margin',2,'horizontalalignment','right');
% text(f.ax(5),0.03,0.1,strcat({'SW portion'}),'fontsize',11,'unit','normalized','horizontalalignment',...
%      'left');  
% text(f.ax(5),0.03,0.3,strcat({'99% CI'}),'fontsize',11,'unit','normalized','horizontalalignment',...
%      'left'); 
% 
%  %%% save figure
% print(f.fig,'-dpdf',strcat(rstpath,'/vdist_unshifted_mainSW+NE_2portions_ind',num2str(nfam),'.pdf'));
% 
% %%% Also do the test to the 2 portions
% % 1. the NE portion
% % a. use significance level 0.05; frequency weights (2), and unequal variance (0)
% [testnea,probnea,statsnea] = ttest2_wt(ivdistswlt,iwtswlt,ivdistnegt,iwtnegt, 0.05, 2, 0)
% 
% % 2. the SW portion
% % a. use significance level 0.05; frequency weights (2), and unequal variance (0)
% [testswa,probswa,statsswa] = ttest2_wt(ivdistswgt,iwtswgt,ivdistnelt,iwtnelt, 0.05, 2, 0)


% %% set and save random seed state for reproduction 
% N = 1000;   % sampling times
% seed1 = cell(1,N);
% seed2 = cell(1,N);
% seed3 = cell(1,N);
% seed4 = cell(1,N);
% for i = 1: N
% %     seed{i} = rng(i,'twister');    % fix and perserve the random seed state fro reproductivity
%     seed1{i} = rng('shuffle','twister');
% end


% %% Try the random sampling from data PDF, all detections
% %%% for testing code segment, see 'mig_linear_fit_LZB_addfam_Miscellaneous.m'
% 
% %%% 1. for SW and NE group at main LZB region, i.e. the NW part of the study region
% % compute the mean of the randoms generated and the misfits of PDFs using the above 3 ways
% % in order: normpdf, ksdensity, normfit
% % NE group
% % NOTE: 'meanrandne' tracks the weighted data mean e.g. meanne = wt_mean(vdistne,wtne) in all cases 
% N = 1000;   % sampling times
% [probdiff,mdiff,refval] = random_sampling_test(N,vdistne,wtne,pdfxlocne,pdfvalne,...
%                                     vdistsw,wtsw,pdfxlocsw,pdfvalsw);
% % plot the distribution of generated random numbers
% binw = 0.05;
% f1 = plt_random_mean_combo(mdiff,refval,binw);
% 
% %%% 2. for SW and NE group at fam 002 region
% [probdiff002,mdiff002,refval002] = random_sampling_test(N,vdist002ne,wt002ne,...
%                                              pdfxloc002ne,pdfval002ne,vdist002sw,wt002sw,...
%                                              pdfxloc002sw,pdfval002sw);
% % plot the distribution of generated random numbers
% binw = 0.1;
% f2 = plt_random_mean_combo(mdiff002,refval002,binw);
% 
% %%%%%%%%%%%%%%%%%%%%%% ABONDONED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % stack a student's t distribution on the normfitting generated random number distribution
% % hold(f11.ax(6), 'on');
% % xx = -0.4: 0.001: 0.4;
% % stdmean = std(meandiff(:,3), 1);
% % synpdf = tpdf(xx, N-1) * stdmean / sqrt(N);
% % synt = tpdf(xx, statsmaina.df) * N * 0.02;
% 
% % xa = meanrandne(:,3)-meanne;
% % xb = meanrandsw(:,3)-meansw;
% % synt = xx;  %/ sqrt(var(xa)/N + var(xb)/N);
% % syndf = (var(xa)/N + var(xb)/N).^2 / ( (var(xa)/N)^2/(N-1) + (var(xb)/N)^2/(N-1) );
% % synpdf = tpdf(synt, syndf);
% % plot(xx, synpdf, 'k-', 'linew', 1.5);
% 
% 
% % Try the random sampling from data PDF, independent detections only
% 
% %%% 1. for SW and NE group at main LZB region, i.e. the NW part of the study region
% % compute the mean of the randoms generated and the misfits of PDFs using the above 3 ways
% % in order: normpdf, ksdensity, normfit
% % NE group
% % NOTE: 'meanrandne' tracks the weighted data mean e.g. meanne = wt_mean(vdistne,wtne) in all cases 
% 
% N = 1000;   % sampling times
% [probdiffi,mdiffi,refvali] = random_sampling_test(N,ivdistne,iwtne,pdfxlocine,pdfvaline,...
%                                                   ivdistsw,iwtsw,pdfxlocisw,pdfvalisw);
% % plot the distribution of generated random numbers
% binw = 0.05;
% f3 = plt_random_mean_combo(mdiffi,refvali,binw);
% 
% %%% 2. for SW and NE group at fam 002 region
% [probdiff002i,mdiff002i,refval002i] = random_sampling_test(N,ivdist002ne,iwt002ne,...
%                                              pdfxloc002ine,pdfval002ine,ivdist002sw,iwt002sw,...
%                                              pdfxloc002isw,pdfval002isw);
% % plot the distribution of generated random numbers
% binw = 0.1;
% f4 = plt_random_mean_combo(mdiff002i,refval002i,binw);                                          


%% plot the migration direction distribution
theta = angbest;
theta = deg2rad(theta);
figure
% p1 = polarhistogram(theta(angbest>0 & angbest<=90),'binw',pi/36,'normalization','count','facec',...
%                [0 1 1],'facea',0.8);
% binedge = p1.BinEdges + pi/72;
% delete(p1);
binedge = deg2rad(2.5: 5: 92.5);
polarhistogram(theta(angbest>0 & angbest<=90),'binedges',binedge,'normalization','count',...
               'facec',[0 1 1],'facea',0.8); hold on
binedge = deg2rad(2.5: 5: 92.5 + 90);           
polarhistogram(theta(angbest>90 & angbest<=180),'binedges',binedge,'normalization','count',...
               'facec','k','facea',0.8);
binedge = deg2rad(2.5: 5: 92.5 + 180);           
polarhistogram(theta(angbest>180 & angbest<=270),'binedges',binedge,'normalization','count',...
               'facec',[1 96/255 0],'facea',0.8);
binedge = deg2rad(2.5: 5: 92.5 + 270);           
polarhistogram(theta(angbest>270 & angbest<=360),'binedges',binedge,'normalization','count',...
               'facec','k','facea',0.2);
ax = gca;
ax.RTickLabelRotation = 45;
ax.ThetaMinorTick = 'on';
ax.TickLength = [0.02 0];
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.FontSize = 12;
ax.RTick = 0:1:5;
% ax.LineWidth = 1.5;
ax.Box = 'on';
ax.GridAlpha = 0.3;
text(deg2rad(60),3.5,num2str(sum(angbest>0 & angbest<=90)),'HorizontalAlignment','center',...
     'FontSize',11);           
text(deg2rad(130),3.5,num2str(sum(angbest>90 & angbest<=180)),'HorizontalAlignment','center',...
     'FontSize',11);           
text(deg2rad(225),3.5,num2str(sum(angbest>180 & angbest<=270)),'HorizontalAlignment','center',...
     'FontSize',11);
text(deg2rad(315),3.5,num2str(sum(angbest>270 & angbest<=360)),'HorizontalAlignment','center',...
     'FontSize',11); 
text(0,4.6,'N','HorizontalAlignment',"center",'FontSize',14);
text(deg2rad(90),4.6,'E','HorizontalAlignment',"center",'FontSize',14);
text(deg2rad(180),4.6,'S','HorizontalAlignment',"center",'FontSize',14);
text(deg2rad(270),4.6,'W','HorizontalAlignment',"center",'FontSize',14);
text(deg2rad(315),7.4,'b','HorizontalAlignment',"center",'FontSize',12,'EdgeColor','k','Margin',2);
text(deg2rad(45),7,'LZB','HorizontalAlignment',"center",'FontSize',14,'EdgeColor','k','Margin',2);

%%% save figure
print('-dpdf',strcat(rstpath,'/LZB.',num2str(nfam),'autortmv2.mig.direc.distri.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));
% clear ax


%% summary of results with error bars from lf fitting
CIFcn = @(x,p)std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * tinv(abs([0,1]-(1-p/100)/2),...
              sum(~isnan(x(:)))-1) + mean(x(:),'omitnan');
f5.fig=figure;
f5.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 6;   % maximum height allowed is 11 inches
set(f5.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f5.ax(isub) = subplot(nrow,ncol,isub);
end

set(f5.ax(1), 'position', [ 0.1, 0.1, 0.4, 0.85]);
set(f5.ax(2), 'position', [ 0.55, 0.1, 0.4, 0.85]);

% index of main LZB region RTMs
indne = [11, 38, 42, 45, 47];
indsw = [8,9,12,14,19,29,30,34,36,37,41,46];
indlzb = [fliplr(indne), fliplr(indsw)];

% index of 002 region RTMs
ind002ne = [21,22];
ind002sw = [4,25];
ind002 = [fliplr(ind002ne), fliplr(ind002sw)];


% index of all other RTMs
indall = 1: 1: size(trange,1);
tmp = union(indlzb,ind002);
dump = [2,3,7,23,40];
tmp2 = union(tmp, dump);
indelse = setdiff(indall, tmp2);
indelse = fliplr(indelse);

ax = f5.ax(1);
hold(ax,'on');
patarea = [0 -100;
           -100 -100;
           -100 100;
            0 100;
            0 -100];
% patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
plot(ax,[0 0],[-100 100],'--','color',[0.6 0.6 0.6],'linew',0.7);

% main LZB region RTMs 
ii = length(indlzb)+ length(ind002);
indplt = indlzb;
for i = 1: length(indplt)
    hfN = numhf(indplt(i));
    
    lfN = numlf(indplt(i));
    lf=resproplf(indplt(i),1:lfN);
    lfm = mean(lf);
    lfstd = std(lf);
    lfsem = lfstd/sqrt(lfN);
    ci95 = tinv([0.025 0.975], lfN-1);
    offcor = bsxfun(@times, lfsem, ci95(:));
    lfci95 = lfm+offcor;
    lfCI95 = CIFcn(lf,95);
    
    ypos = ii-i+1;
    if i > 5
        ypos = ii-i;
    end
    if angbest(indplt(i))>=0 && angbest(indplt(i))<90       % NE
        xpos = offset(indplt(i));   % sign unchanged
        e1 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[0 1 1],'CapSize',5);
    elseif angbest(indplt(i))>=180 && angbest(indplt(i))<270    % SW
        xpos = -offset(indplt(i));   % sign flipped
        e3 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[1 96/255 0],'CapSize',5);
    end
    text(ax,xpos+0.1,ypos+0.4,num2str(indplt(i)),'fontsize',7,...
        'HorizontalAlignment','left');
    text(ax,4.9,ypos,strcat(num2str(hfN),'; ',num2str(lfN)),'fontsize',7,...
        'HorizontalAlignment','right');
end
plot(ax,[-100 100],[ii-length(indne) ii-length(indne)],'--','color',[0 0 0],'linew',0.5);
plot(ax,[-100 100],[ypos-1 ypos-1],'-','color',[0 0 0],'linew',1.5);
text(ax,0.04,0.96,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.75,0.96,'Main LZB region','FontSize',10,'unit','normalized',...
     'horizontalalignment','center','EdgeColor','k','Margin',2);

    
% 002 region RTMs
ii = ypos-4;
indplt = ind002;
for i = 1: length(indplt)
    hfN = numhf(indplt(i));
    
    lfN = numlf(indplt(i));
    lf=resproplf(indplt(i),1:lfN);
    lfm = mean(lf);
    lfstd = std(lf);
    lfsem = lfstd/sqrt(lfN);
    ci95 = tinv([0.025 0.975], lfN-1);
    offcor = bsxfun(@times, lfsem, ci95(:));
    lfci95 = lfm+offcor;
    lfCI95 = CIFcn(lf,95);
    
    ypos = ii-i+1;
    if i == 2
        ypos = -1.5;
    elseif i == 3
        ypos = -3.5;
    elseif i == 4
        ypos = -5;
    end
    if angbest(indplt(i))>=0 && angbest(indplt(i))<90       % NE
        xpos = offset(indplt(i));   % sign unchanged
        e1 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[0 1 1],'CapSize',5);
    elseif angbest(indplt(i))>=180 && angbest(indplt(i))<270    % SW
        xpos = -offset(indplt(i));   % sign flipped
        e3 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[1 96/255 0],'CapSize',5);
    elseif angbest(indplt(i))>=90 && angbest(indplt(i))<180     % SE
        xpos = offset(indplt(i));   % sign unchanged
        e2 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor','k','CapSize',5);
    end
    text(ax,xpos+0.1,ypos+0.4,num2str(indplt(i)),'fontsize',7,...
        'HorizontalAlignment','left');
    text(ax,4.9,ypos,strcat(num2str(hfN),'; ',num2str(lfN)),'fontsize',7,...
        'HorizontalAlignment','right');
end
scatter(ax,-0.66, -0.5, 12,'MarkerEdgeColor','k','linewidth',1);
scatter(ax,-1.27, -0.5, 12,'MarkerEdgeColor','k','linewidth',1);
plot(ax,[-0.66 -1.13], [-0.5 -0.5], 'k:','linewidth',1.5);
scatter(ax, 0.22, -4, 12,'MarkerEdgeColor','k','linewidth',1);
scatter(ax,-0.19, -4, 12,'MarkerEdgeColor','k','linewidth',1);
plot(ax,[0.22 -0.19], [-4 -4], 'k:','linewidth',1.5);
plot(ax,[-100 100],[ii-length(ind002ne)-0.5 ii-length(ind002ne)-0.5],'--','color',[0 0 0],'linew',0.5);
text(ax,0.75,0.3,'Fam 002 region','FontSize',10,'unit','normalized',...
     'horizontalalignment','center','EdgeColor','k','Margin',2);
yran = [-8,24];
xran = [-5,5];
% axis(ax,[xran yran]);

xvect = [-1.5 -4.5];
yvect = [-7 -7];
drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
xvect = [1.5 4.5];
yvect = [-7 -7];
drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
text(ax,0.2,0.06,strcat({'SW'}),'fontsize',10,'fontw','bold','unit','normalized',...
     'horizontalalignment','center');
text(ax,0.8,0.06,strcat({'NE'}),'fontsize',10,'fontw','bold','unit','normalized',...
     'horizontalalignment','center');
 
drawArrow(ax,[-4.3 -1 ],[17 17],xran,yran,'linewidth',1,'linestyle','-','color',[0.65 0.65 0.65]);
text(ax,-4.8,17,'I','FontSize',13,'VerticalAlignment','middle'); 
drawArrow(ax,[-4.3 -1 ],[6 6],xran,yran,'linewidth',1,'linestyle','-','color',[0.65 0.65 0.65]);
text(ax,-4.8,6,'II','FontSize',13,'VerticalAlignment','middle');
drawArrow(ax,[-4.3 -3.6 ],[-1.5 -1.5],xran,yran,'linewidth',1,'linestyle','-','color',...
          [0.65 0.65 0.65]);
text(ax,-4.8,-1.5,'III','FontSize',13,'VerticalAlignment','middle');
 
xlabel(ax,'LF-HF (km)','fontsize',11);
% ylabel(ax,'RTM number in chronological order','fontsize',11);
legend(ax,[e1, e3],{'NE propagation','SW propagation'},'FontSize',7,...
       'Position',[0.149 0.62 0.065 0.05]);       % ,'edgecolor',[0.5 0.5 0.5]
xticks(ax,-5:1:5);
ax.YTick = [];
ax.Box='on';
grid(ax,'on');
set(ax,'GridLineStyle','--');
ax.FontSize = 9;
hold(ax,'off');


% Other RTMs
ax = f5.ax(2);
hold(ax,'on');
patarea = [0 -100;
           -100 -100;
           -100 100;
            0 100;
            0 -100];
patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
ii = length(indlzb)+ length(ind002);
indplt = indelse;
for i = 1: length(indplt)
    hfN = numhf(indplt(i));
    
    lfN = numlf(indplt(i));
    lf=resproplf(indplt(i),1:lfN);
    lfm = mean(lf);
    lfstd = std(lf);
    lfsem = lfstd/sqrt(lfN);
    ci95 = tinv([0.025 0.975], lfN-1);
    offcor = bsxfun(@times, lfsem, ci95(:));
    lfci95 = lfm+offcor;
    lfCI95 = CIFcn(lf,95);
    
    ypos = ii-i+1;
    if angbest(indplt(i))>=0 && angbest(indplt(i))<90       % NE
        xpos = offset(indplt(i));   % sign unchanged
        e1 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[0 1 1],'CapSize',5);
    elseif angbest(indplt(i))>=90 && angbest(indplt(i))<180     % SE
        xpos = offset(indplt(i));   % sign unchanged
        e2 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor','k','CapSize',5);
    elseif angbest(indplt(i))>=180 && angbest(indplt(i))<270    % SW
        xpos = offset(indplt(i));   % sign unchanged
        e3 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[1 96/255 0],'CapSize',5);
    end
    text(ax,xpos+0.1,ypos+0.4,num2str(indplt(i)),'fontsize',7,...
        'HorizontalAlignment','left');
    text(ax,4.9,ypos,strcat(num2str(hfN),'; ',num2str(lfN)),'fontsize',7,...
        'HorizontalAlignment','right');
end

text(ax,0.04,0.96,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.75,0.96,'Other RTMs','FontSize',10,'unit','normalized',...
     'horizontalalignment','center','EdgeColor','k','Margin',2);
% text(ax,0.7,0.97,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.2,0.07,'LF','FontSize',10,'unit','normalized','HorizontalAlignment','center','fontw',...
    'bold');
text(ax,0.2,0.04,'lags','FontSize',10,'unit','normalized','HorizontalAlignment','center','fontw',...
    'bold');
text(ax,0.8,0.07,'LF','FontSize',10,'unit','normalized','HorizontalAlignment','center','fontw',...
    'bold');
text(ax,0.8,0.04,'leads','FontSize',10,'unit','normalized','HorizontalAlignment','center','fontw',...
    'bold');
yran = [-3,24];
xran = [-5,5];
drawArrow(ax,[-4.2 -2.5],[5 5],xran,yran,'linewidth',1.5,'linestyle','-','color',[0.65 0.65 0.65]);
text(ax,-4.8,5,'IV','FontSize',13,'VerticalAlignment','middle','backgroundcolor','w',...
     'Margin',0.5);
% axis(ax,[xran yran]);
xlabel(ax,'LF-HF (km)','fontsize',11);
% ylabel(ax,'RTM number in chronological order','fontsize',11);
legend(ax,[e1, e2, e3],{'NE propagation','SE propagation','SW propagation'},'FontSize',7,...
       'Position',[0.599 0.6 0.065 0.05]);  
xticks(ax,-5:1:5);
ax.YTick = [];
% yticks(0:1:8);
ax.Box='on';
grid(ax,'on');
set(ax,'GridLineStyle','--');
ax.FontSize = 9;
hold(ax,'off');

%%% save figure
print(f5.fig,'-dpdf',strcat(rstpath,'/LZB',num2str(nfam),'autortmv2.mig.lfit.sum.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));

% the RTM correspondence with PGC trio
indcor = [21,22,23,25];
offcor = zeros(length(indcor),3);
for i = 1: length(indcor)
    hfN = numhf(indcor(i));
    
    lfN = numlf(indcor(i));
    lf=resproplf(indcor(i),1:lfN);
    lfm = mean(lf);
    lfstd = std(lf);
    lfsem = lfstd/sqrt(lfN);
    ci95 = tinv([0.025 0.975], lfN-1);
    offcor(i,2:3) = bsxfun(@times, lfsem, ci95(:));
    offcor(i,1) = offset(indcor(i));
end


%% summary plots in histogram
f6.fig=figure;
set(f6.fig,'Position',[scrsz(3)/10 2*scrsz(4)/10 0.2*scrsz(3) 0.3*scrsz(4)]);
f6.ax=gca;
hold on
% histogram(f6.ax,offset(offset<=0),'binwidth',0.5,'normalization','probability','facecolor',...
%           [0.4 0.4 0.4]);
% histogram(f6.ax,offset(offset>0),'binwidth',0.5,'normalization','probability','facecolor',...
%           [1 1 1]);  
histogram(f6.ax,offset,'binwidth',0.5,'normalization','probability','facecolor',...
          [0.3 0.3 0.3]);
plot(f6.ax,[0 0],[-100 100],'b--','linew',1);
text(f6.ax,0.82,0.92,'LZB','FontSize',14,'unit','normalized','EdgeColor','k','Margin',2);
fracneg = vpa(length(find(offset<=0))/length(offset),2)*100;
text(f6.ax,0.76,0.72,sprintf('%d%%',fracneg),'FontSize',14,'unit','normalized','HorizontalAlignment','center');
text(f6.ax,0.76,0.65,'negative','FontSize',14,'unit','normalized','HorizontalAlignment','center');
xlabel(f6.ax,'Offset between HF and LF (LF-HF) (km)','fontsize',12);
ylabel(f6.ax,'Fraction','fontsize',12);
ylim(f6.ax,[0,0.3]);
xlim(f6.ax,[-4,4]);
% yticks(0:1:8);
box on
grid on
set(f6.ax,'GridLineStyle','--');

%%% save figure
print(f6.fig,'-dpdf',strcat(rstpath,'/LZB',num2str(nfam),'autortmv2.mig.lfit.hist.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));



%% suumary plots in barplot
f7.fig=figure;
set(f7.fig,'Position',[scrsz(3)/10 2*scrsz(4)/10 0.2*scrsz(3) 0.3*scrsz(4)]);
f7.ax=gca;
hold on

xloc = [-1.75 -1.25 -0.75 -0.25 0.25 0.75 1.25 3.25];
h = [];
for i = 1: length(xloc)
    ibar = find(offset>xloc(i)-0.25 & offset<=xloc(i)+0.25);
    ib1 = find(angbest(ibar)>=0 & angbest(ibar)<90);
    h1 = length(ib1)/length(offset);
    ib2 = find(angbest(ibar)>=90 & angbest(ibar)<180);
    h2 = length(ib2)/length(offset);
    ib3 = find(angbest(ibar)>=180 & angbest(ibar)<270);
    h3 = length(ib3)/length(offset);
    h = [h;h1 h2 h3];
end
H = bar(f7.ax,xloc,h,1,'stacked','facea',0.8);
H(1).FaceColor = [0 1 1];
H(2).FaceColor = 'k';
H(3).FaceColor = [1 96/255 0];
plot(f7.ax,[0 0],[-100 100],'b--','linew',1);
text(f7.ax,0.82,0.92,'LZB','FontSize',14,'unit','normalized','EdgeColor','k','Margin',2);
fracneg = vpa(length(find(offset<=0))/length(offset),2)*100;
fracpos = vpa(length(find(offset>0))/length(offset),2)*100;
text(f7.ax,0.76,0.72,sprintf('%d%%',fracpos),'FontSize',14,'unit','normalized','HorizontalAlignment','center');
text(f7.ax,0.76,0.65,'positve','FontSize',14,'unit','normalized','HorizontalAlignment','center');
text(f7.ax,0.18,0.72,sprintf('%d%%',fracneg),'FontSize',14,'unit','normalized','HorizontalAlignment','center');
text(f7.ax,0.18,0.65,'negative','FontSize',14,'unit','normalized','HorizontalAlignment','center');
text(0.04,0.93,'c','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
xlabel(f7.ax,'Offset between HF and LF (LF-HF) (km)','fontsize',12);
ylabel(f7.ax,'Fraction','fontsize',12);
ylim(f7.ax,[0,0.3]);
xlim(f7.ax,[-4,4]);
xticks(-4:1:4);
box on
grid on
set(f7.ax,'GridLineStyle','--');

%%% save figure
print(f7.fig,'-dpdf',strcat(rstpath,'/LZB',num2str(nfam),'autortmv2.mig.lfit.barplot.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% combine the above three figures into one for the paper
CIFcn = @(x,p)std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * tinv(abs([0,1]-(1-p/100)/2),...
              sum(~isnan(x(:)))-1) + mean(x(:),'omitnan');
          
[scrsz, res] = pixelperinch(1);
f111.fig=figure;
f111.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 7;   % maximum height allowed is 11 inches
set(f111.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol-1
    f111.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
% ax = gca;
set(f111.ax(1), 'position', [ 0.1, 0.1, 0.45, 0.85]);
set(f111.ax(2), 'position', [ 0.62, 0.1, 0.35, 0.4]);
set(f111.ax(3), 'position', [ 0.62, 0.55, 0.35, 0.4]);


ax = f111.ax(1);
hold(ax,'on');
patarea = [0 0;
           -5 0;
           -5 50;
            0 50;
            0 0];
patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
% plot([0 0],[-100 100],'--','color',[0.5 0.5 0.5]);
for i = 1: length(trange)
    hfN = numhf(i);
    hf=resprophf(i,1:hfN);
    hfm = mean(hf);
    hfstd = std(hf);
    hfsem = hfstd/sqrt(hfN);
    ci95 = tinv([0.025 0.975], hfN-1);
    hfci95n = bsxfun(@times, hfsem, ci95(:));
    hfci95 = hfm+hfci95n;
    hfCI95 = CIFcn(hf,95);
    
    lfN = numlf(i);
    lf=resproplf(i,1:lfN);
    lfm = mean(lf);
    lfstd = std(lf);
    lfsem = lfstd/sqrt(lfN);
    ci95 = tinv([0.025 0.975], lfN-1);
    lfci95n = bsxfun(@times, lfsem, ci95(:));
    lfci95 = lfm+lfci95n;
    lfCI95 = CIFcn(lf,95);
    
    offset(i) = (intcptproplf(i)-intcptprophf(i));
    ypos = i;
    if angbest(i)>=0 && angbest(i)<90
        e1 = errorbar(ax,offset(i),ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[0 1 1]);
    elseif angbest(i)>=90 && angbest(i)<180
        e2 = errorbar(ax,offset(i),ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor','k');
    elseif angbest(i)>=180 && angbest(i)<270
        e3 = errorbar(ax,offset(i),ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[1 96/255 0]);
    end
    text(ax,4.9,ypos,strcat(num2str(hfN),'; ',num2str(lfN)),'fontsize',7,...
        'HorizontalAlignment','right');
    %%% next time use 'units', 'normalized' instead of defaults; also 
end

text(ax,0.02,0.97,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.7,0.97,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.2,0.13,'LF','FontSize',11,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(ax,0.2,0.1,'lags','FontSize',11,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(ax,0.8,0.13,'LF','FontSize',11,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(ax,0.8,0.1,'leads','FontSize',11,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(ax,-1.6,7.4,'I','FontSize',12);
yran = [0,50];
xran = [-5,5];
drawArrow(ax,[-3.5 -1.4 ],[11 11],xran,yran,'linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);
drawArrow(ax,[-3.5 -2.8 ],[13 13],xran,yran,'linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);
drawArrow(ax,[-3.5 -1 ],[35 35],xran,yran,'linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);
text(ax,-4.2,11,'II','FontSize',12,'VerticalAlignment','middle','backgroundcolor','w',...
     'Margin',0.5);
text(ax,-4.2,13,'III','FontSize',12,'VerticalAlignment','middle','backgroundcolor','w',...
     'Margin',0.5);
text(ax,-4.2,35,'IV','FontSize',12,'VerticalAlignment','middle','backgroundcolor','w',...
     'Margin',0.5);

xlabel(ax,'Offset between HF and LF (LF-HF) (km)','fontsize',11);
ylabel(ax,'RTM number in chronological order','fontsize',11);
legend(ax,[e1, e2, e3],{'NE propagation','SE propagation','SW propagation'},'FontSize',8,...
        'Position',[0.113 0.52 0.16 0.08],'edgecolor',[0.5 0.5 0.5]);
% ylim(ax,[0,50]);
% xlim(ax,[-5,5]);
xticks(ax,-5:1:5);
% yticks(0:1:8);
ax.Box='on';
grid(ax,'on');
set(ax,'GridLineStyle','--');
ax.FontSize = 9;
hold(ax,'off');


ax = f111.ax(3);
theta = angbest;
theta = deg2rad(theta);
polarhistogram(theta(angbest>=0 & angbest<90),'binw',pi/36,'normalization','count','facec',...
               [0 1 1],'facea',0.8); hold on
polarhistogram(theta(angbest>=90 & angbest<180),'binw',pi/36,'normalization','count','facec',...
               'k','facea',0.8);
polarhistogram(theta(angbest>=180 & angbest<270),'binw',pi/36,'normalization','count','facec',...
               [1 96/255 0],'facea',0.8);
polarhistogram(theta(angbest>=270 & angbest<360),'binw',pi/36,'normalization','count','facec',...
               'w','facea',0.8);
ax = gca;
ax.RTickLabelRotation = 45;
ax.ThetaMinorTick = 'on';
ax.TickLength = [0.02 0];
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.FontSize = 9;
ax.RTick = 0:1:5;
% ax.LineWidth = 1.5;
ax.Box = 'on';
ax.GridAlpha = 0.3;
text(ax,deg2rad(60),3.5,num2str(sum(angbest>=0 & angbest<90)),'HorizontalAlignment','center',...
     'FontSize',11);           
text(ax,deg2rad(130),3.5,num2str(sum(angbest>=90 & angbest<180)),'HorizontalAlignment','center',...
     'FontSize',11);           
text(ax,deg2rad(225),3.5,num2str(sum(angbest>=180 & angbest<270)),'HorizontalAlignment','center',...
     'FontSize',11); 
text(ax,0,4.6,'N','HorizontalAlignment',"center",'FontSize',12);
text(ax,deg2rad(90),4.6,'E','HorizontalAlignment',"center",'FontSize',12);
text(ax,deg2rad(180),4.6,'S','HorizontalAlignment',"center",'FontSize',12);
text(ax,deg2rad(270),4.6,'W','HorizontalAlignment',"center",'FontSize',12);
text(ax,deg2rad(308),7.3,'b','HorizontalAlignment',"center",'FontSize',11,'EdgeColor','k','Margin',2);
% text(ax,deg2rad(45),7,'LZB','HorizontalAlignment',"center",'FontSize',12,'EdgeColor','k','Margin',2);
hold(ax,'off');

ax = f111.ax(2);
hold(ax,'on');
% xloc = [-1.75 -1.25 -0.75 -0.25 0.25 0.75 1.25 3.25];
xloc = -3.75:0.5:3.75;
h = [];
for i = 1: length(xloc)
    ibar = find(offset>xloc(i)-0.25 & offset<=xloc(i)+0.25);
    ib1 = find(angbest(ibar)>=0 & angbest(ibar)<90);
    h1 = length(ib1)/length(offset);
    ib2 = find(angbest(ibar)>=90 & angbest(ibar)<180);
    h2 = length(ib2)/length(offset);
    ib3 = find(angbest(ibar)>=180 & angbest(ibar)<270);
    h3 = length(ib3)/length(offset);
    h = [h;h1 h2 h3];
end
H = bar(ax,xloc,h,1,'stacked','facea',0.8);
H(1).FaceColor = [0 1 1];
H(2).FaceColor = 'k';
H(3).FaceColor = [1 96/255 0];
patarea = [0 0;
           -4 0;
           -4 0.3;
            0 0.3;
            0 0];
patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
fracneg = vpa(length(find(offset<=0))/length(offset),2)*100;
fracpos = vpa(length(find(offset>0))/length(offset),2)*100;
med = median(offset);
plot(ax,[med med],[-100 100],'b--','linew',1.5);
text(ax,0.76,0.72,sprintf('%d%%',fracpos),'FontSize',12,'unit','normalized','HorizontalAlignment','center');
text(ax,0.76,0.65,'positve','FontSize',12,'unit','normalized','HorizontalAlignment','center');
text(ax,0.18,0.72,sprintf('%d%%',fracneg),'FontSize',12,'unit','normalized','HorizontalAlignment','center');
text(ax,0.18,0.65,'negative','FontSize',12,'unit','normalized','HorizontalAlignment','center');
text(ax,0.03,0.93,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
xlabel(ax,'Offset between HF and LF (LF-HF) (km)','fontsize',11);
ylabel(ax,'Fraction','fontsize',11);
ylim(ax,[0,0.3]);
xlim(ax,[-4,4]);
xticks(ax,-4:1:4);
ax.Box='on';
grid(ax,'on');
set(ax,'GridLineStyle','--');
ax.FontSize = 9;
hold(ax,'off');


print(f111.fig,'-dpdf',strcat(rstpath,'/LZB',num2str(nfam),'autortmv2.sum.mig.lfit.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));


%% way 2 to combine the above three figures into one, plot separate histograms instead of stack bars 
CIFcn = @(x,p)std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * tinv(abs([0,1]-(1-p/100)/2),...
              sum(~isnan(x(:)))-1) + mean(x(:),'omitnan');
          
f22.fig=figure;
f22.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 7;   % maximum height allowed is 11 inches
set(f22.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

nrow = 2;
ncol = 3;
for isub = 1:nrow*ncol-1
    f22.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
set(f22.ax(1), 'position', [ 0.1, 0.1, 0.45, 0.85]);  % summary
set(f22.ax(2), 'position', [ 0.62, 0.1, 0.35, 0.14]);     % hist red
set(f22.ax(3), 'position', [ 0.62, 0.24, 0.35, 0.14]);    % hist black
set(f22.ax(4), 'position', [ 0.62, 0.38, 0.35, 0.14]);    % hist blue
set(f22.ax(5), 'position', [ 0.62, 0.55, 0.35, 0.4]);     % azimuth

ax = f22.ax(1);
hold(ax,'on');
patarea = [0 0;
           -5 0;
           -5 50;
            0 50;
            0 0];
patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
% plot([0 0],[-100 100],'--','color',[0.5 0.5 0.5]);
for i = 1: length(trange)
    hfN = numhf(i);
    hf=resprophf(i,1:hfN);
    hfm = mean(hf);
    hfstd = std(hf);
    hfsem = hfstd/sqrt(hfN);
    ci95 = tinv([0.025 0.975], hfN-1);
    hfci95n = bsxfun(@times, hfsem, ci95(:));
    hfci95 = hfm+hfci95n;
    hfCI95 = CIFcn(hf,95);
    
    lfN = numlf(i);
    lf=resproplf(i,1:lfN);
    lfm = mean(lf);
    lfstd = std(lf);
    lfsem = lfstd/sqrt(lfN);
    ci95 = tinv([0.025 0.975], lfN-1);
    lfci95n = bsxfun(@times, lfsem, ci95(:));
    lfci95 = lfm+lfci95n;
    lfCI95 = CIFcn(lf,95);
    
    offset(i) = (intcptproplf(i)-intcptprophf(i));
    ypos = i;
    if angbest(i)>=0 && angbest(i)<90
        e1 = errorbar(ax,offset(i),ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[0 1 1]);
    elseif angbest(i)>=90 && angbest(i)<180
        e2 = errorbar(ax,offset(i),ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor','k');
    elseif angbest(i)>=180 && angbest(i)<270
        e3 = errorbar(ax,offset(i),ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[1 96/255 0]);
    end
    text(ax,4.9,ypos,strcat(num2str(hfN),'; ',num2str(lfN)),'fontsize',7,...
        'HorizontalAlignment','right');
    %%% next time use 'units', 'normalized' instead of defaults; also 
end

text(ax,0.02,0.97,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.7,0.97,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.2,0.13,'LF','FontSize',11,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(ax,0.2,0.1,'lags','FontSize',11,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(ax,0.8,0.13,'LF','FontSize',11,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(ax,0.8,0.1,'leads','FontSize',11,'unit','normalized','HorizontalAlignment','center','fontw','bold');
% text(ax,-1.6,7.4,'I','FontSize',12);
yran = [0,50];
xran = [-5,5];
drawArrow(ax,[-3.5 -1.4 ],[11 11],xran,yran,'linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);
drawArrow(ax,[-3.5 -2.8 ],[13 13],xran,yran,'linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);
drawArrow(ax,[-3.5 -1 ],[35 35],xran,yran,'linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);
text(ax,-4.2,11,'II','FontSize',12,'VerticalAlignment','middle','backgroundcolor','w',...
     'Margin',0.5);
text(ax,-4.2,13,'III','FontSize',12,'VerticalAlignment','middle','backgroundcolor','w',...
     'Margin',0.5);
text(ax,-4.2,35,'IV','FontSize',12,'VerticalAlignment','middle','backgroundcolor','w',...
     'Margin',0.5);

xlabel(ax,'Offset between HF and LF (LF-HF) (km)','fontsize',11);
ylabel(ax,'RTM number in chronological order','fontsize',11);
legend(ax,[e1, e2, e3],{'NE propagation','SE propagation','SW propagation'},'FontSize',8,...
        'Position',[0.113 0.52 0.16 0.08],'edgecolor',[0.5 0.5 0.5]);
% ylim(ax,[0,50]);
% xlim(ax,[-5,5]);
xticks(ax,-5:1:5);
% yticks(0:1:8);
ax.Box='on';
grid(ax,'on');
set(ax,'GridLineStyle','--');
ax.FontSize = 9;
hold(ax,'off');

ax = f22.ax(5);
theta = angbest;
theta = deg2rad(theta);
polarhistogram(theta(angbest>=0 & angbest<90),'binw',pi/36,'normalization','count','facec',...
               [0 1 1],'facea',0.8); hold on
polarhistogram(theta(angbest>=90 & angbest<180),'binw',pi/36,'normalization','count','facec',...
               'k','facea',0.8);
polarhistogram(theta(angbest>=180 & angbest<270),'binw',pi/36,'normalization','count','facec',...
               [1 96/255 0],'facea',0.8);
polarhistogram(theta(angbest>=270 & angbest<360),'binw',pi/36,'normalization','count','facec',...
               'w','facea',0.8);
ax = gca;
ax.RTickLabelRotation = 45;
ax.ThetaMinorTick = 'on';
ax.TickLength = [0.02 0];
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.FontSize = 9;
ax.RTick = 0:1:5;
% ax.LineWidth = 1.5;
ax.Box = 'on';
ax.GridAlpha = 0.3;
text(ax,deg2rad(60),3.5,num2str(sum(angbest>=0 & angbest<90)),'HorizontalAlignment','center',...
     'FontSize',11);           
text(ax,deg2rad(130),3.5,num2str(sum(angbest>=90 & angbest<180)),'HorizontalAlignment','center',...
     'FontSize',11);           
text(ax,deg2rad(225),3.5,num2str(sum(angbest>=180 & angbest<270)),'HorizontalAlignment','center',...
     'FontSize',11); 
text(ax,0,4.6,'N','HorizontalAlignment',"center",'FontSize',12);
text(ax,deg2rad(90),4.6,'E','HorizontalAlignment',"center",'FontSize',12);
text(ax,deg2rad(180),4.6,'S','HorizontalAlignment',"center",'FontSize',12);
text(ax,deg2rad(270),4.6,'W','HorizontalAlignment',"center",'FontSize',12);
text(ax,deg2rad(308),7.3,'b','HorizontalAlignment',"center",'FontSize',11,'EdgeColor','k','Margin',2);
% text(ax,deg2rad(45),7,'LZB','HorizontalAlignment',"center",'FontSize',12,'EdgeColor','k','Margin',2);
hold(ax,'off');


ax = f22.ax(4);
hold(ax,'on');
% xloc = [-1.75 -1.25 -0.75 -0.25 0.25 0.75 1.25 3.25];
xloc = -3.75:0.5:3.75;
h = [];
for i = 1: length(xloc)
    ibar = find(offset>xloc(i)-0.25 & offset<=xloc(i)+0.25);
    ib1 = find(angbest(ibar)>=0 & angbest(ibar)<90);
    h1 = length(ib1)/length(offset);
    ib2 = find(angbest(ibar)>=90 & angbest(ibar)<180);
    h2 = length(ib2)/length(offset);
    ib3 = find(angbest(ibar)>=180 & angbest(ibar)<270);
    h3 = length(ib3)/length(offset);
    h = [h;h1 h2 h3];
end
H = bar(ax,xloc,h(:,1),1,'stacked','facea',0.8);
H.FaceColor = [0 1 1];
patarea = [0 0;
           -4 0;
           -4 0.3;
            0 0.3;
            0 0];
patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
offset1 = offset(angbest>=0 & angbest<90);
fracneg = vpa(length(find(offset1<=0))/length(offset1),2)*100;
fracpos = vpa(length(find(offset1>0))/length(offset1),2)*100;
med = median(offset1);
plot(ax,[med med],[-100 100],'b--','linew',1.5);
text(ax,0.8,0.7,sprintf('%d%%',fracpos),'FontSize',10,'unit','normalized','HorizontalAlignment','center');
text(ax,0.8,0.52,'positve','FontSize',10,'unit','normalized','HorizontalAlignment','center');
text(ax,0.18,0.7,sprintf('%d%%',fracneg),'FontSize',10,'unit','normalized','HorizontalAlignment','center');
text(ax,0.18,0.52,'negative','FontSize',10,'unit','normalized','HorizontalAlignment','center');
text(ax,0.03,0.82,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
ylim(ax,[0,0.2]);
xlim(ax,[-4,4]);
xticks(ax,-4:1:4);
yticks(ax,[0.05 0.1 0.15]);
ax.XTickLabel = [];
ax.Box='on';
grid(ax,'on');
set(ax,'GridLineStyle','--');
ax.FontSize = 9;
hold(ax,'off');


ax = f22.ax(3);
hold(ax,'on');
H = bar(ax,xloc,h(:,2),1,'stacked','facea',0.8);
H.FaceColor = 'k';
patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
offset2 = offset(angbest>=90 & angbest<180);
fracneg = vpa(length(find(offset2<=0))/length(offset2),2)*100;
fracpos = vpa(length(find(offset2>0))/length(offset2),2)*100;
med = median(offset2);
plot(ax,[med med],[-100 100],'b--','linew',1.5);
text(ax,0.8,0.7,sprintf('%d%%',fracpos),'FontSize',10,'unit','normalized','HorizontalAlignment','center');
text(ax,0.8,0.52,'positve','FontSize',10,'unit','normalized','HorizontalAlignment','center');
text(ax,0.18,0.7,sprintf('%d%%',fracneg),'FontSize',10,'unit','normalized','HorizontalAlignment','center');
text(ax,0.18,0.52,'negative','FontSize',10,'unit','normalized','HorizontalAlignment','center');
ylim(ax,[0,0.2]);
xlim(ax,[-4,4]);
xticks(ax,-4:1:4);
yticks(ax,[0.05 0.1 0.15]);
ax.XTickLabel = [];
ax.Box='on';
grid(ax,'on');
set(ax,'GridLineStyle','--');
ax.FontSize = 9;
ylabel(ax,'Fraction','fontsize',10);
hold(ax,'off');

ax = f22.ax(2);
hold(ax,'on');
H = bar(ax,xloc,h(:,3),1,'stacked','facea',0.8);
H.FaceColor = [1 96/255 0];
patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
offset3 = offset(angbest>=180 & angbest<270);
fracneg = vpa(length(find(offset3<=0))/length(offset3),2)*100;
fracpos = vpa(length(find(offset3>0))/length(offset3),2)*100;
med = median(offset3);
plot(ax,[med med],[-100 100],'b--','linew',1.5);
text(ax,0.8,0.7,sprintf('%d%%',fracpos),'FontSize',10,'unit','normalized','HorizontalAlignment','center');
text(ax,0.8,0.52,'positve','FontSize',10,'unit','normalized','HorizontalAlignment','center');
text(ax,0.18,0.7,sprintf('%d%%',fracneg),'FontSize',10,'unit','normalized','HorizontalAlignment','center');
text(ax,0.18,0.52,'negative','FontSize',10,'unit','normalized','HorizontalAlignment','center');
ylim(ax,[0,0.2]);
xlim(ax,[-4,4]);
xticks(ax,-4:1:4);
yticks(ax,[0.05 0.1 0.15]);
ax.Box='on';
grid(ax,'on');
set(ax,'GridLineStyle','--');
ax.FontSize = 9;
xlabel(ax,'Offset between HF and LF (LF-HF) (km)','fontsize',10);
hold(ax,'off');

print(f22.fig,'-dpdf',strcat(rstpath,'/LZB',num2str(nfam),'autortmv2.sum.v2.mig.lfit.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));


%% plot selected migrations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind = [11,13,35];
xran = [-15 25];
yran = [-20 20];
for i = 1: length(ind)
% for i = 1: 21
%     i=3;
    disp(trange(ind(i),:));
    indhf = find(hftime(:,13)==trange(ind(i),1) & hftime(:,15)>=trange(ind(i),2) & ...
                 hftime(:,15)<=trange(ind(i),3));
    mighf = hftime(indhf,:);
    
    mighfdum = mighf;
    for j = 1: size(mighf,1)
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(ind(i))-90),0,0);
        mighfdum(j,1) = newx;
        mighfdum(j,2) = newy;
    end
    
    indlf = find(lftime(:,13)==trange(ind(i),1) & lftime(:,15)>=trange(ind(i),2) & ...
                 lftime(:,15)<=trange(ind(i),3));
    miglf = lftime(indlf,:); 
    
    miglfdum = miglf;
    for j = 1: size(miglf,1)
        x0 = miglfdum(j,1);
        y0 = miglfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(ind(i))-90),0,0);
        miglfdum(j,1) = newx;
        miglfdum(j,2) = newy;
    end
    
    %%% define and position the figure frame and axes of each plot
    f8.fig = figure;
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 9;   % maximum height allowed is 11 inches
    set(f8.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

    nrow = 3;
    ncol = 2;    
    for isub = 1:nrow*ncol
        f8.ax(isub) = subplot(nrow,ncol,isub);
    end

    %%% reposition
    set(f8.ax(1), 'position', [ 0.1, 0.64, 0.32, 0.32]);
    set(f8.ax(2), 'position', [ 0.5, 0.64, 0.32, 0.32]);
    set(f8.ax(3), 'position', [ 0.1, 0.35, 0.32, 0.22]);
    set(f8.ax(4), 'position', [ 0.5, 0.35, 0.32, 0.22]);
    set(f8.ax(5), 'position', [ 0.1, 0.078, 0.32, 0.22]);
    set(f8.ax(6), 'position', [ 0.5, 0.078, 0.32, 0.22]);
    
    % subplot 1 of figure ind(i)
    hold(f8.ax(1),'on');
    plot(f8.ax(1),[-100 100],[0 0],'k--');
    plot(f8.ax(1),[0 0],[-100 100],'k--');
    f8.ax(1).FontSize = 9;
    mighf = sortrows(mighf,-15);
    scatter(f8.ax(1),mighf(:,1),mighf(:,2), 20, mighf(:,15)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f8.ax(1),'jet');
    c=colorbar(f8.ax(1),'SouthOutside');
    pos = f8.ax(1).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    juldate = num2str(trange(ind(i),1));
    yr = str2double(juldate(1:4));
    date = str2double(juldate(5:end));
    a = jul2dat(yr,date);
    mo = a(1);
    if mo == 9 
        mo = {' Sep. '};
    elseif mo == 7
        mo = {' Jul. '};
    else
        mo = {' Mar. '};
    end
    day = num2str(a(2));
    yr = num2str(a(3));
    c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
    c.Label.FontSize = 11;
    caxis(f8.ax(1),[trange(ind(i),2)/3600 trange(ind(i),3)/3600])
    text(f8.ax(1),0.85,0.1,'HF','FontSize',12,'unit','normalized');
    text(f8.ax(1),0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f8.ax(1),0.04,0.1,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f8.ax(1).Box = 'on';
    grid(f8.ax(1), 'on');
    axis(f8.ax(1), 'equal');
    f8.ax(1).GridLineStyle = '--';
    f8.ax(1).XAxisLocation = 'top';
    medxlf = median(mighf(:,1));
    medylf = median(mighf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(ind(i)));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f8.ax(1),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(ind(i)));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f8.ax(1),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
    xticks(f8.ax(1),xran(1):5:xran(2));
    yticks(f8.ax(1),yran(1):5:yran(2));    
    xlabel(f8.ax(1),'E (km)','fontsize',11);
    ylabel(f8.ax(1),'N (km)','fontsize',11);
    if i == 1
        text(f8.ax(1),0.5,0.1,'II','FontSize',12,'unit','normalized');
    elseif i == 2
        text(f8.ax(1),0.5,0.1,'III','FontSize',12,'unit','normalized');
    else
        text(f8.ax(1),0.5,0.1,'IV','FontSize',12,'unit','normalized');
    end
    text(f8.ax(1),0.83,0.95,strcat(num2str(size(mighf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f8.ax(1),0.83,0.90,strcat({'in '},num2str(trange(ind(i),3)-trange(ind(i),2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(mighf,1)/(trange(ind(i),3)-trange(ind(i),2)));
    text(f8.ax(1),0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    hold(f8.ax(1),'off');

    % subplot 2 of figure ind(i)
    hold(f8.ax(2),'on');
    plot(f8.ax(2),[-100 100],[0 0],'k--');
    plot(f8.ax(2),[0 0],[-100 100],'k--');
    f8.ax(2).FontSize = 9;
    miglf = sortrows(miglf,-15);
    scatter(f8.ax(2),miglf(:,1),miglf(:,2), 20, miglf(:,15)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f8.ax(2),'jet');
    c=colorbar(f8.ax(2),'SouthOutside');
    pos = f8.ax(2).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
    c.Label.FontSize = 11;
    caxis(f8.ax(2),[trange(ind(i),2)/3600 trange(ind(i),3)/3600])
    text(f8.ax(2),0.85,0.1,'LF','FontSize',12,'unit','normalized');
    text(f8.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f8.ax(2),0.04,0.1,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f8.ax(2).Box = 'on';
    grid(f8.ax(2), 'on');
    axis(f8.ax(2), 'equal');
    f8.ax(2).GridLineStyle = '--';
    f8.ax(2).XAxisLocation = 'top';
    medxlf = median(miglf(:,1));
    medylf = median(miglf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(ind(i)));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f8.ax(2),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(ind(i)));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f8.ax(2),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
    xticks(f8.ax(2),xran(1):5:xran(2));
    yticks(f8.ax(2),yran(1):5:yran(2));    
    xlabel(f8.ax(2),'E (km)','fontsize',11);
    ylabel(f8.ax(2),'N (km)','fontsize',11);
    if i == 1
        text(f8.ax(2),0.5,0.1,'II','FontSize',12,'unit','normalized');
    elseif i == 2
        text(f8.ax(2),0.5,0.1,'III','FontSize',12,'unit','normalized');
    else
        text(f8.ax(2),0.5,0.1,'IV','FontSize',12,'unit','normalized');
    end
    text(f8.ax(2),0.83,0.95,strcat(num2str(size(miglf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f8.ax(2),0.83,0.9,strcat({'in '},num2str(trange(ind(i),3)-trange(ind(i),2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(miglf,1)/(trange(ind(i),3)-trange(ind(i),2)));
    text(f8.ax(2),0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    hold(f8.ax(2),'off');
    
    % subplot 3 of figure ind(i)
    hold(f8.ax(3),'on');
    f8.ax(3).FontSize = 9;
    scatter(f8.ax(3),mighfdum(:,15)/3600,mighfdum(:,1),20,[0.6 0.6 0.6],'filled','o',...
            'MarkerEdgeColor','k');   
    scatter(f8.ax(3),miglfdum(:,15)/3600,miglfdum(:,1),20,[0.6 1 1],'filled','o',...
            'MarkerEdgeColor','k');
%     scatter(f8.ax(3),mighfdum(:,10)/3600,mighfdum(:,1),30,[105/255 105/255 105/255],'filled','o');  %, 'MarkerEdgeColor', 'w')    
%     scatter(f8.ax(3),miglfdum(:,10)/3600,miglfdum(:,1),30,miglfdum(:,11),'filled','o');  %, 'MarkerEdgeColor', 'w')
    
    % create fit object
    %indpool = [12;13;16;19;21;22;29;36;40];
    if i == 2   % 13 in original order 
        mintime = max(min(mighfdum(:,15)), min(miglfdum(:,15)));
        maxtime = min(max(mighfdum(:,15)), max(miglfdum(:,15)));
        mighfdum2 = mighfdum(mighfdum(:,15)>=mintime & mighfdum(:,15)<=maxtime, :);
        miglfdum2 = miglfdum(miglfdum(:,15)>=mintime & miglfdum(:,15)<=maxtime, :);
    else
        mighfdum2 = mighfdum;
        miglfdum2 = miglfdum;
    end
    
    [fitobjhfprop,gofhfprop,outphfprop] = fit(mighfdum2(:,15)/3600, mighfdum2(:,1),fttpfree,'Robust',...
                                              'Bisquare','StartPoint',[1 1]);
    % output fit parameters
    coefprop = coeffvalues(fitobjhfprop);
    slopeprophf(ind(i)) = coefprop(1);
    intcptprophf(ind(i)) = coefprop(2);
    fitprophf = feval(fitobjhfprop,mighfdum2(:,15)/3600);
    plot(f8.ax(3),mighfdum2(:,15)/3600,fitprophf,'-','linewidth',2,'color',[0.6 0.6 0.6]);
        
%     intcptproplf = ones(size(miglfdum(:,end)))\(miglfdum(:,1)-slopeprophf*miglfdum(:,end)/3600); % direct LS fit
%     fitproplf = slopeprophf*miglfdum(:,end)/3600+intcptproplf;
    a=slopeprophf(ind(i));
    fttpfix = fittype( @(b,x) a*x+b);
    [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum2(:,15)/3600, miglfdum2(:,1),fttpfix,'Robust',...
                                              'Bisquare','StartPoint',1);
%     [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum(:,10)/3600, miglfdum(:,1),fttpfix,'Robust',...
%                                               'Bisquare','Weights',miglfdum(:,11));
    fitproplf = feval(fitobjlfprop,miglfdum2(:,15)/3600);
    intcptproplf(ind(i)) = coeffvalues(fitobjlfprop);
    offset(ind(i)) = intcptproplf(ind(i))-intcptprophf(ind(i));
    plot(f8.ax(3),miglfdum2(:,15)/3600,fitproplf,'-','linewidth',2,'color',[0.6 1 1]);
    f8.ax(3).Box = 'on';
    grid(f8.ax(3), 'on');
    f8.ax(3).GridLineStyle = '--';
    xran1 = [trange(ind(i),2)/3600 trange(ind(i),3)/3600];
    yran1 = [round(min([mighfdum(:,1);miglfdum(:,1)]))-1 ...
             round(max([mighfdum(:,1);miglfdum(:,1)]))+1];
    xlim(f8.ax(3),xran1);
    ylim(f8.ax(3),yran1);
    text(f8.ax(3),0.04,0.91,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    xlabel(f8.ax(3),'Time (hr)','fontsize',11);
    ylabel(f8.ax(3),'Dist. along prop. (km)','fontsize',11); 
    hold(f8.ax(3),'off');
    
    % compute the HF weights in robust linear regression, see NOTES above
    rhf = outphfprop.residuals;   % usual residuals
    x = [mighfdum2(:,15)/3600];
    hatmat = x*inv(x'*x)*x';
    h = zeros(size(hatmat,1),1);    % leverage of least square
    for jj = 1 : size(hatmat,1)
        h(jj) = hatmat(jj,jj);
    end    
    radj = rhf./sqrt(1-h);      % adjusted residuals
    K = 4.685;
    s = mad(rhf,1)/0.6745;
    u = radj/(K*s);
    wthf = zeros(length(u),1);    % rubust weight of next iteration
    for jj = 1 : length(u)
        if abs(u(jj)) < 1
            wthf(jj) = (1-(u(jj))^2)^2;
        else
            wthf(jj) = 0;
        end
    end
    
    % compute the LF weights in robust linear regression, see NOTES above
    rlf = outplfprop.residuals;   % usual residuals
    x = [miglfdum2(:,15)/3600];
    hatmat = x*inv(x'*x)*x';
    h = zeros(size(hatmat,1),1);    % leverage of least square
    for jj = 1 : size(hatmat,1)
        h(jj) = hatmat(jj,jj);
    end       
    radj = rlf./sqrt(1-h);      % adjusted residuals
    K = 4.685;
    s = mad(rlf,1)/0.6745;
    u = radj/(K*s);
    wtlf = zeros(length(u),1);    % rubust weight of next iteration
    for jj = 1 : length(u)
        if abs(u(jj)) < 1
            wtlf(jj) = (1-(u(jj))^2)^2;
        else
            wtlf(jj) = 0;
        end
    end
    
    % subplot 4 of figure ind(i)
    hold(f8.ax(4),'on');
    f8.ax(4).FontSize = 9;
    numhf(ind(i)) = length(mighfdum2(:,1));
    numlf(ind(i)) = length(miglfdum2(:,1));
    resprophf(ind(i),1:numhf(ind(i))) = outphfprop.residuals;
    resproplf(ind(i),1:numlf(ind(i))) = outplfprop.residuals+(intcptproplf(ind(i))-intcptprophf(ind(i)));
    [xlochf, a, b, pdfvalhf] = weighted_bar_pdf(resprophf(ind(i),1:numhf(ind(i))), wthf, 0.5, 'int');
    hfHdl = bar(f8.ax(4),xlochf, pdfvalhf,1,'stacked');
    hfHdl(1).FaceColor = [0.6 0.6 0.6];
    [xloclf, ~, ~, pdfvallf] = weighted_bar_pdf(resproplf(ind(i),1:numlf(ind(i))), wtlf, 0.5, 'int');
    lfHdl = bar(f8.ax(4),xloclf, pdfvallf,1,'stacked','facea',0.6);
    lfHdl(1).FaceColor = [0.6 1 1];
%     histogram(f8.ax(4),resprophf(ind(i),1:numhf(ind(i))),'binwidth',0.5,'normalization','pdf','facecolor',...
%               [1 0 0]);
%     histogram(f8.ax(4),resproplf(ind(i),1:numlf(ind(i))),'binwidth',0.5,'normalization','pdf','facecolor',...
%               [0.6 1 1],'facealpha',0.6);
%     pdhf = fitdist(resprophf(ind(i),1:numhf(ind(i)))','Normal');    % fit a distribution of residual, assuming it is normal distributed
%     pdlf = fitdist(resproplf(ind(i),1:numlf(ind(i)))','Normal');
%     muhf = pdhf.mu;    % fitted paramaters
%     mulf = pdlf.mu;
%     sigmahf = pdhf.sigma;
%     sigmalf = pdlf.sigma;
%     cihf = paramci(pdhf,'Alpha',0.05);     % give the 95% confidence interval of distribution param
%     cilf = paramci(pdlf,'Alpha',0.05);    
%     pdffithf = pdf(pdhf,min(resprophf(ind(i),1:numhf(ind(i)))):0.05:max(resprophf(ind(i),1:numhf(ind(i)))));
%     pdffitlf = pdf(pdlf,min(resproplf(ind(i),1:numlf(ind(i)))):0.05:max(resproplf(ind(i),1:numlf(ind(i)))));
%     plot(f8.ax(4),min(resprophf(ind(i),1:numhf(ind(i)))):0.05:max(resprophf(ind(i),1:numhf(ind(i)))),pdffithf,'-',...
%          'color',[1 0 0],'linewidth',2);
%     plot(f8.ax(4),min(resproplf(ind(i),1:numlf(ind(i)))):0.05:max(resproplf(ind(i),1:numlf(ind(i)))),pdffitlf,'-',...
%          'color',[0.6 1 1],'linewidth',2);
    text(f8.ax(4),0.04,0.91,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f8.ax(4),0.4,0.92,sprintf('LF-HF = %.2f km',offset(ind(i))),'FontSize',11,'unit','normalized');
    f8.ax(4).Box = 'on';
    grid(f8.ax(4), 'on');
    f8.ax(4).GridLineStyle = '--';
    ymax = f8.ax(4).YLim(2)+0.1;
    ylim(f8.ax(4),[0 ymax]);
    xlabel(f8.ax(4),'Residual in prop. (km)','fontsize',11);
    ylabel(f8.ax(4),'PDF estimate','fontsize',11);
    hold(f8.ax(4),'off');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% It seems that median of the previous detections is a promising way to show offsets as well,
    %%% so try to deal with it the same way
    medhf = [];
    for j = 1: size(mighfdum2,1)
        medhf(j,1) = median(mighfdum2(1:j,1));
        medhf(j,2) = median(mighfdum2(1:j,2));
        medhf(j,3) = mighfdum2(j,15);
    end
    
    medlf = [];
    for j = 1: size(miglfdum2,1)
        medlf(j,1) = median(miglfdum2(1:j,1));
        medlf(j,2) = median(miglfdum2(1:j,2));
        medlf(j,3) = miglfdum2(j,15);
    end
       
    
    % subplot 5 of figure ind(i)
    hold(f8.ax(5),'on');
    f8.ax(5).FontSize = 9;
    scatter(f8.ax(5),medhf(:,3)/3600,medhf(:,1),20,[0.6 0.6 0.6],'filled','o','MarkerEdgeColor',...
            'k');
    scatter(f8.ax(5),medlf(:,3)/3600,medlf(:,1),20,[0.6 1 1],'filled','o','MarkerEdgeColor',...
            'k');            
    f8.ax(5).Box = 'on';
    grid(f8.ax(5), 'on');
    f8.ax(5).GridLineStyle = '--';
    xran1 = [trange(ind(i),2)/3600 trange(ind(i),3)/3600];
    yran1 = [round(min([medhf(:,1);medlf(:,1)]))-1 ...
             round(max([medhf(:,1);medlf(:,1)]))+1];
    xlim(f8.ax(5),xran1);
    ylim(f8.ax(5),yran1);
    text(f8.ax(5),0.97,0.1,'e','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2,...
         'HorizontalAlignment','right');
     xlabel(f8.ax(5),'Time (hr)','fontsize',11);
    ylabel(f8.ax(5),'Med. dist. along prop. (km)','fontsize',11);
    hold(f8.ax(5),'off');
    
    % subplot 6 of figure ind(i)
    hold(f8.ax(6),'on');
    f8.ax(6).FontSize = 9;
    scatter(f8.ax(6),medhf(:,3)/3600,medhf(:,2),20,[0.6 0.6 0.6],'filled','o','MarkerEdgeColor',...
            'k');
    scatter(f8.ax(6),medlf(:,3)/3600,medlf(:,2),20,[0.6 1 1],'filled','o','MarkerEdgeColor',...
            'k');          
    f8.ax(6).Box = 'on';
    grid(f8.ax(6), 'on');
    f8.ax(6).GridLineStyle = '--';
    xran1 = [trange(ind(i),2)/3600 trange(ind(i),3)/3600];
    yran1 = [round(min([medhf(:,2);medlf(:,2)]))-1 ...
             round(max([medhf(:,2);medlf(:,2)]))+1];
    xlim(f8.ax(6),xran1);
    ylim(f8.ax(6),yran1);
    text(f8.ax(6),0.97,0.1,'f','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2,...
         'HorizontalAlignment','right');
    xlabel(f8.ax(6),'Time (hr)','fontsize',11);
    ylabel(f8.ax(6),'Med. dist. along ort. (km)','fontsize',11);
    hold(f8.ax(6),'off');
    
    %%% save figure
    print(f8.fig,'-dpdf',strcat(rstpath,'/LZB',num2str(nfam),'autortmv2.sel.mig.proj.lfit.meddist',...
          num2str(trange(ind(i),1)),'_',num2str(trange(ind(i),2)),'-',num2str(trange(ind(i),3)),...
          '.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),...
          '.pdf'));
end











    
    
    
    