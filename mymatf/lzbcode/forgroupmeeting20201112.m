% function mig_linear_fit_LZB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to use the poropagation direction estimate fit the hf
% the migration linearly to get the propagation slope and error distribution.
%
% this is for adding new fams 006, 001
%
% Chao Song, chaosong@princeton.edu
% First created date:   2020/09/10
% Last modified date:   2020/09/10
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


trange = [
           2004195,4.15e+4,4.60e+4;
           2004196,7.55e+4,7.70e+4;
           2004197,0.03e+4,0.35e+4;
           2004197,0.50e+4,0.75e+4;
           2004197,4.60e+4,4.72e+4;
           2004197,8.35e+4,8.46e+4;
           2004197,8.488e+4,8.52e+4;
           2004197,8.55e+4,8.64e+4;
           2004198,0.58e+4,0.98e+4;
           2004198,1.90e+4,2.18e+4;
           2004198,5.40e+4,5.80e+4;
           2004198,7.35e+4,7.60e+4;
           2004198,8.35e+4,8.64e+4;
           2004199,0.50e+4,0.63e+4;
           2004199,4.61e+4,4.91e+4;
           2004199,4.98e+4,5.08e+4;
           2004199,8.10e+4,8.20e+4;
           2004200,1.38e+4,1.80e+4;
           2004200,1.85e+4,1.99e+4;
           2004203,1.60e+4,2.20e+4;
           2005255,3.42e+4,3.58e+4;
           2005255,6.70e+4,6.85e+4;
           2005255,7.50e+4,7.57e+4;
           2005255,8.42e+4,8.56e+4;
           2005256,0.35e+4,0.50e+4;
           2005256,0.80e+4,0.98e+4;
           2005256,2.15e+4,2.23e+4;
           2005256,2.82e+4,2.95e+4;
           2005256,3.45e+4,3.53e+4;
           2005256,5.20e+4,5.26e+4;
           2005256,7.275e+4,7.40e+4;
           2005256,7.60e+4,7.80e+4;
           2005257,0.60e+4,0.70e+4;
           2005257,2.10e+4,2.50e+4;
           2005257,3.90e+4,4.40e+4;
           2005257,6.10e+4,6.40e+4;
           2005257,7.36e+4,7.80e+4;
           2005258,3.52e+4,3.66e+4;
           2005258,3.80e+4,4.00e+4;
           2005259,0.20e+4,0.40e+4;
           2005259,0.50e+4,0.77e+4;
           2005259,3.70e+4,4.26e+4;
           2005259,7.20e+4,7.58e+4;
           2005260,0.36e+4,0.43e+4;
           2005260,0.43e+4,0.57e+4;
           2005260,5.63e+4,5.82e+4;
           2005261,0.82e+4,1.10e+4;];

xran = [-15 25];
yran = [-20 20];

angbest = [265;
215;
250;
255;
235;
220;
245;
250;
235;
230;
75;
235;
240;
250;
120;
195;
85;
115;
255;
115;
60;
35;  % weird, check timing, lf is also bad, use neither
70;
240;
215;
215;
115;
110;
250;
260;
270;
250;
195;
245;
165;
245;
240;
90;
185;
135;
255;
80;
240;
245;
70;
235;
70;];


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


% create fit object with free constraints
fttpfree = fittype( @(a,b,x) a*x+b);
        
for i = 1: length(trange)
% for i = 1: 211
%     i=3;
    disp(trange(i,:));
    indhf = find(hftime(:,13)==trange(i,1) & hftime(:,15)>=trange(i,2) & ...
                 hftime(:,15)<=trange(i,3));
    mighf = hftime(indhf,:);
    
    mighfdum = mighf;
    for j = 1: size(mighf,1)
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
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
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        miglfdum(j,1) = newx;
        miglfdum(j,2) = newy;
    end
    
    %%% define and position the figure frame and axes of each plot
    f.fig=figure;
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 9;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    nrow = 3;
    ncol = 2;    
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
    end

    %%% reposition
    set(f.ax(1), 'position', [ 0.1, 0.64, 0.32, 0.32]);
    set(f.ax(2), 'position', [ 0.5, 0.64, 0.32, 0.32]);
    set(f.ax(3), 'position', [ 0.1, 0.35, 0.32, 0.22]);
    set(f.ax(4), 'position', [ 0.5, 0.35, 0.32, 0.22]);
    set(f.ax(5), 'position', [ 0.1, 0.078, 0.32, 0.22]);
    set(f.ax(6), 'position', [ 0.5, 0.078, 0.32, 0.22]);
    
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
    hold(f.ax(2),'off');
    
    % subplot 3 of figure i
    hold(f.ax(3),'on');
    f.ax(3).FontSize = 9;
    scatter(f.ax(3),mighfdum(:,15)/3600,mighfdum(:,1),20,[1 0 0],'filled','o',...
            'MarkerEdgeColor','k');   
    scatter(f.ax(3),miglfdum(:,15)/3600,miglfdum(:,1),20,[0.6 1 1],'filled','o',...
            'MarkerEdgeColor','k');
%     scatter(f.ax(3),mighfdum(:,10)/3600,mighfdum(:,1),30,[105/255 105/255 105/255],'filled','o');  %, 'MarkerEdgeColor', 'w')    
%     scatter(f.ax(3),miglfdum(:,10)/3600,miglfdum(:,1),30,miglfdum(:,11),'filled','o');  %, 'MarkerEdgeColor', 'w')
    % create fit object
    
    indpool = [12;13;16;19;21;22;29;32;36;40];
    indpool = [indpool; 3;5;7;17;20;23;25;26;30];
    if ismember(i,indpool)
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
    slopeprophf(i) = coefprop(1);
    intcptprophf(i) = coefprop(2);
    fitprophf = feval(fitobjhfprop,mighfdum2(:,15)/3600);
    plot(f.ax(3),mighfdum2(:,15)/3600,fitprophf,'-','linewidth',2,'color',[1 0 0]);
        
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
    yran1 = [round(min([mighfdum(:,1);miglfdum(:,1)]))-1 ...
             round(max([mighfdum(:,1);miglfdum(:,1)]))+1];
    xlim(f.ax(3),xran1);
    ylim(f.ax(3),yran1);
    text(f.ax(3),0.04,0.91,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    xlabel(f.ax(3),'Time (hr)','fontsize',11);
    ylabel(f.ax(3),'Dist. along prop. (km)','fontsize',11); 
    hold(f.ax(3),'off');
    
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
    
    % get the vertical distance from LF detections to the fitted HF line
    fittemp = feval(fitobjhfprop,miglfdum2(:,15)/3600);
    verttemp = miglfdum2(:,1) - fittemp;     % can be + or -
    vertdistraw(1: length(verttemp), i) = verttemp;     % raw vertical distance
    verttempwt = wtlf.*verttemp;      % weighted vertical distance by multiplying the weight
    vertdistwt(1: length(verttempwt), i) = verttempwt;
    weight(1: length(verttempwt), i) =  wtlf;
    
    % for mig 9, there are a few parts to throw away
    if i == 9
        ind9 = find(miglfdum2(:,1)>-10);
    % for mig 4, only the early portion that passes the 002 region is used    
    elseif i == 4
        ind4 = find(miglfdum2(:,15)<5836);
    end
        
    % get the information of propagation vector indicating range etc.
    [tmpx,tmpy] = coordinate_rot(medxhf,medyhf,-(angbest(i)-90),0,0);
    medallhf(i,1:4) = [medxhf medyhf tmpx tmpy];
    [tmpx,tmpy] = coordinate_rot(medxlf,medylf,-(angbest(i)-90),0,0);
    medalllf(i,1:4) = [medxlf medylf tmpx tmpy];
    if i == 4
        aaa = median(mighf(mighf(:,15)<5836,1));
        bbb = median(mighf(mighf(:,15)<5836,2));
        [tmpx,tmpy] = coordinate_rot(aaa,bbb,-(angbest(i)-90),0,0);
        medallhf(i,1:4) = [aaa bbb tmpx tmpy];
        aaa = median(miglf(miglf(:,15)<5836,1));
        bbb = median(miglf(miglf(:,15)<5836,2));
        [tmpx,tmpy] = coordinate_rot(aaa,bbb,-(angbest(i)-90),0,0);
        medalllf(i,1:4) = [aaa bbb tmpx tmpy];
    end
    
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
    
    if i == 4
        ranvechf98(i,1:2) = [-18.3 -3.5];
        ranveclf98(i,1:2) = [-21.4 -3.5];
    end

    ranvechf99(i,1:2) = [prctile(mighfdum(:,1),1) prctile(mighfdum(:,1),99)];
    ranveclf99(i,1:2) = [prctile(miglfdum(:,1),1) prctile(miglfdum(:,1),99)];
            
    %%% find the independent LF detections, i.e. non-overlapping time window
    [miglfindpen,indindpen] = find_independent_detection(miglfdum2, 0);
    inumlf(i) = length(indindpen);
    indlfindpen(1: inumlf(i), i) =  indindpen; 
    
    % subplot 4 of figure i
    hold(f.ax(4),'on');
    f.ax(4).FontSize = 9;
    numhf(i) = length(mighfdum2(:,1));
    numlf(i) = length(miglfdum2(:,1));
    resprophf(i,1:numhf(i)) = outphfprop.residuals;
    resproplf(i,1:numlf(i)) = outplfprop.residuals+(intcptproplf(i)-intcptprophf(i));
    [xlochf, a, b, pdfvalhf] = weighted_bar_pdf(resprophf(i,1:numhf(i)), wthf, 0.5, 'int');
    hfHdl = bar(f.ax(4),xlochf, pdfvalhf,1,'stacked');
    hfHdl(1).FaceColor = [1 0 0];
    [xloclf, ~, ~, pdfvallf] = weighted_bar_pdf(resproplf(i,1:numlf(i)), wtlf, 0.5, 'int');
    lfHdl = bar(f.ax(4),xloclf, pdfvallf,1,'stacked','facea',0.6);
    lfHdl(1).FaceColor = [0.6 1 1];
%     histogram(f.ax(4),resprophf(i,1:numhf(i)),'binwidth',0.5,'normalization','pdf','facecolor',...
%               [1 0 0]);
%     histogram(f.ax(4),resproplf(i,1:numlf(i)),'binwidth',0.5,'normalization','pdf','facecolor',...
%               [0.6 1 1],'facealpha',0.6);

%     pdhf = fitdist(resprophf(i,1:numhf(i))','Normal');    % fit a distribution of residual, assuming it is normal distributed
%     pdlf = fitdist(resproplf(i,1:numlf(i))','Normal');
%     muhf = pdhf.mu;    % fitted paramaters
%     mulf = pdlf.mu;
%     sigmahf = pdhf.sigma;
%     sigmalf = pdlf.sigma;
%     cihf = paramci(pdhf,'Alpha',0.05);     % give the 95% confidence interval of distribution param
%     cilf = paramci(pdlf,'Alpha',0.05);    
%     pdffithf = pdf(pdhf,min(resprophf(i,1:numhf(i))):0.05:max(resprophf(i,1:numhf(i))));
%     pdffitlf = pdf(pdlf,min(resproplf(i,1:numlf(i))):0.05:max(resproplf(i,1:numlf(i))));
%     plot(f.ax(4),min(resprophf(i,1:numhf(i))):0.05:max(resprophf(i,1:numhf(i))),pdffithf,'-',...
%          'color',[1 0 0],'linewidth',2);
%     plot(f.ax(4),min(resproplf(i,1:numlf(i))):0.05:max(resproplf(i,1:numlf(i))),pdffitlf,'-',...
%          'color',[0.6 1 1],'linewidth',2);
    text(f.ax(4),0.04,0.91,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(4),0.4,0.92,sprintf('LF-HF = %.2f km',offset(i)),'FontSize',11,'unit','normalized');
    f.ax(4).Box = 'on';
    grid(f.ax(4), 'on');
    f.ax(4).GridLineStyle = '--';
    ymax = f.ax(4).YLim(2)+0.1;
    ylim(f.ax(4),[0 ymax]);
%     ylim(f.ax(4),[0 1]);
    xlabel(f.ax(4),'Residual in prop. (km)','fontsize',11);
    ylabel(f.ax(4),'PDF estimate','fontsize',11);
    hold(f.ax(4),'off');
    
    
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
       
    
    % subplot 5 of figure i
    hold(f.ax(5),'on');
    f.ax(5).FontSize = 9;
    scatter(f.ax(5),medhf(:,3)/3600,medhf(:,1),20,[1 0 0],'filled','o','MarkerEdgeColor',...
            'k');
    scatter(f.ax(5),medlf(:,3)/3600,medlf(:,1),20,[0.6 1 1],'filled','o','MarkerEdgeColor',...
            'k');            
    f.ax(5).Box = 'on';
    grid(f.ax(5), 'on');
    f.ax(5).GridLineStyle = '--';
    xran1 = [trange(i,2)/3600 trange(i,3)/3600];
    yran1 = [round(min([medhf(:,1);medlf(:,1)]))-1 ...
             round(max([medhf(:,1);medlf(:,1)]))+1];
    xlim(f.ax(5),xran1);
    ylim(f.ax(5),yran1);
    text(f.ax(5),0.97,0.1,'e','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2,...
         'HorizontalAlignment','right');
    xlabel(f.ax(5),'Time (hr)','fontsize',11);
    ylabel(f.ax(5),'Med. dist. along prop. (km)','fontsize',11);
    hold(f.ax(5),'off');
    
    % subplot 6 of figure i
    hold(f.ax(6),'on');
    f.ax(6).FontSize = 9;
    scatter(f.ax(6),medhf(:,3)/3600,medhf(:,2),20,[1 0 0],'filled','o','MarkerEdgeColor',...
            'k');
    scatter(f.ax(6),medlf(:,3)/3600,medlf(:,2),20,[0.6 1 1],'filled','o','MarkerEdgeColor',...
            'k');          
    f.ax(6).Box = 'on';
    grid(f.ax(6), 'on');
    f.ax(6).GridLineStyle = '--';
    xran1 = [trange(i,2)/3600 trange(i,3)/3600];
    yran1 = [round(min([medhf(:,2);medlf(:,2)]))-1 ...
             round(max([medhf(:,2);medlf(:,2)]))+1];
    xlim(f.ax(6),xran1);
    ylim(f.ax(6),yran1);
    text(f.ax(6),0.97,0.1,'f','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2,...
         'HorizontalAlignment','right');
    xlabel(f.ax(6),'Time (hr)','fontsize',11);
    ylabel(f.ax(6),'Med. dist. along ort. (km)','fontsize',11);
    hold(f.ax(6),'off');
    
%     %%% save figure
%     print(f.fig,'-dpdf',strcat(rstpath,'/LZB.',num2str(nfam),'.mig.proj.lfit.meddist',...
%           num2str(trange(i,1)),'_',num2str(trange(i,2)),'-',num2str(trange(i,3)),...
%           '.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),'.',...
%           num2str(ccminlf),'.pdf'));
        
end

    

%% obtain the vertical distance from lf points to hf line in certain regions

% NOTE here, 2020/09/22, u can choose which vertical distance to use
% vertdist = vertdistwt;
vertdist = vertdistraw;

% SW migs
indsw = [8,9,12,14,19,29,30,34,36,37,41,46];    % index of SW migrations
vdistsw = [];
wtsw = [];
ivdistsw = [];  % dist and weight of independent detections
iwtsw = [];
for i = 1: length(indsw)
    iind = indlfindpen(1:inumlf(indsw(i)), indsw(i));  % index of independent detection in that migration indsw(i)
    if indsw(i) == 9
        temp1 = vertdist(1:numlf(indsw(i)), indsw(i));
        temp2 = temp1(ind9);
        vdistsw = [vdistsw; temp2];
        temp2 = temp1(intersect(iind, ind9));
        ivdistsw = [ivdistsw; temp2]; 
        temp1 = weight(1:numlf(indsw(i)), indsw(i));
        temp2 = temp1(ind9);
        wtsw = [wtsw; temp2];
        temp2 = temp1(intersect(iind, ind9));
        iwtsw = [iwtsw; temp2];
    else
        temp1 = vertdist(1:numlf(indsw(i)), indsw(i));
        vdistsw = [vdistsw; temp1];
        ivdistsw = [ivdistsw; temp1(iind)];
        temp1 = weight(1:numlf(indsw(i)), indsw(i));
        wtsw = [wtsw; temp1];
        iwtsw = [iwtsw; temp1(iind)];
    end
end
pdsw = fitdist(vdistsw,'Normal');  
musw = pdsw.mu;    % fitted paramaters
sigmasw = pdsw.sigma;
pdffitsw = pdf(pdsw,-50:0.05:50);

% NE migs
indne = [11, 38, 42, 45, 47];    % index of NE migrations
vdistne = [];
wtne = [];
ivdistne = [];  % dist and weight of independent detections
iwtne = [];
for i = 1: length(indne)
    iind = indlfindpen(1:inumlf(indne(i)), indne(i));  % index of independent detection in that migration indsw(i)
    temp1 = vertdist(1:numlf(indne(i)), indne(i));
    vdistne = [vdistne; temp1];
    ivdistne = [ivdistne; temp1(iind)];
    temp1 = weight(1:numlf(indne(i)), indne(i));
    wtne = [wtne; temp1];
    iwtne = [iwtne; temp1(iind)];
end
pdne = fitdist(vdistne,'Normal');  
mune = pdne.mu;    % fitted paramaters
sigmane = pdne.sigma;
pdffitne = pdf(pdne,-50:0.05:50);

% NE+SW migs that pass through similar region
temp = union(indsw,indne);
indswne = union(temp,[22,25]);
vdistswne = [];
wtswne = [];
ivdistswne = [];  % dist and weight of independent detections
iwtswne = [];
for i = 1: length(indswne)
    iind = indlfindpen(1:inumlf(indswne(i)), indswne(i));  % index of independent detection in that migration indsw(i)
    if indswne(i) == 9
        temp1 = vertdist(1:numlf(indswne(i)), indswne(i));
        temp2 = temp1(ind9);
        vdistswne = [vdistswne; temp2];
        temp2 = temp1(intersect(iind, ind9));
        ivdistswne = [ivdistswne; temp2];
        temp1 = weight(1:numlf(indswne(i)), indswne(i));
        temp2 = temp1(ind9);
        wtswne = [wtswne; temp2];
        temp2 = temp1(intersect(iind, ind9));
        iwtswne = [iwtswne; temp2];
    else
        temp1 = vertdist(1:numlf(indswne(i)), indswne(i));
        vdistswne = [vdistswne; temp1];
        ivdistswne = [ivdistswne; temp1(iind)];
        temp1 = weight(1:numlf(indswne(i)), indswne(i));
        wtswne = [wtswne; temp1];
        iwtswne = [iwtswne; temp1(iind)];
    end
end
pdswne = fitdist(vdistswne,'Normal');  
muswne = pdswne.mu;    % fitted paramaters
sigmaswne = pdswne.sigma;
pdffitswne = pdf(pdswne,-50:0.05:50);

% 2 NE and SW migs that pass through fam 002 region
ind002 = [22,25];
vdist002 = [];
wt002 = [];
ivdist002 = [];  % dist and weight of independent detections
iwt002 = [];
for i = 1: length(ind002)
    iind = indlfindpen(1:inumlf(ind002(i)), ind002(i));  % index of independent detection in that migration indsw(i)
    temp1 = vertdist(1:numlf(ind002(i)), ind002(i));
    vdist002 = [vdist002; temp1];
    ivdist002 = [ivdist002; temp1(iind)];
    temp1 = weight(1:numlf(ind002(i)), ind002(i));
    wt002 = [wt002; temp1];
    iwt002 = [iwt002; temp1(iind)];

end
pd002 = fitdist(vdist002,'Normal');  
mu002 = pd002.mu;    % fitted paramaters
sig002 = pd002.sigma;
pdffit002 = pdf(pd002,-50:0.05:50);

% NE mig that passes through fam 002 region
ind002ne = [22];
ind002ne = [22,21,23];
vdist002ne = [];
wt002ne = [];
ivdist002ne = [];  % dist and weight of independent detections
iwt002ne = [];
for i = 1: length(ind002ne)
    iind = indlfindpen(1:inumlf(ind002ne(i)), ind002ne(i));  % index of independent detection in that migration indsw(i)
    temp1 = vertdist(1:numlf(ind002ne(i)), ind002ne(i));
    vdist002ne = [vdist002ne; temp1];
    ivdist002ne = [ivdist002ne; temp1(iind)];
    temp1 = weight(1:numlf(ind002ne(i)), ind002ne(i));
    wt002ne = [wt002ne; temp1];
    iwt002ne = [iwt002ne; temp1(iind)];
end
pd002ne = fitdist(vdist002ne,'Normal');  
mu002ne = pd002ne.mu;    % fitted paramaters
sig002ne = pd002ne.sigma;
pdffit002ne = pdf(pd002ne,-50:0.05:50);

% SW mig that passes through fam 002 region
ind002sw = [25];
ind002sw = [25,4];
vdist002sw = [];
wt002sw = [];
ivdist002sw = [];  % dist and weight of independent detections
iwt002sw = [];
for i = 1: length(ind002sw)
    iind = indlfindpen(1:inumlf(ind002sw(i)), ind002sw(i));  % index of independent detection in that migration indsw(i)
    if ind002sw(i) == 4
        temp1 = vertdist(1:numlf(ind002sw(i)), ind002sw(i));
        temp2 = temp1(ind4);
        vdist002sw = [vdist002sw; temp2];
        temp2 = temp1(intersect(iind, ind4));
        ivdist002sw = [ivdist002sw; temp2];
        temp1 = weight(1:numlf(ind002sw(i)), ind002sw(i));
        temp2 = temp1(ind4);
        wt002sw = [wt002sw; temp2];
        temp2 = temp1(intersect(iind, ind4));
        iwt002sw = [iwt002sw; temp2];
    else
        temp1 = vertdist(1:numlf(ind002sw(i)), ind002sw(i));
        vdist002sw = [vdist002sw; temp1];
        ivdist002sw = [ivdist002sw; temp1(iind)];
        temp1 = weight(1:numlf(ind002sw(i)), ind002sw(i));
        wt002sw = [wt002sw; temp1];
        iwt002sw = [iwt002sw; temp1(iind)];
    end
end
pd002sw = fitdist(vdist002sw,'Normal');  
mu002sw = pd002sw.mu;    % fitted paramaters
sig002sw = pd002sw.sigma;
pdffit002sw = pdf(pd002sw,-50:0.05:50);

% all migs excluding 
indall = 1: 1: size(trange,1);
indelse = setdiff(indall, indswne);
vdistelse = [];
wtelse = [];
ivdistelse = [];  % dist and weight of independent detections
iwtelse = [];
for i = 1: length(indelse)
    iind = indlfindpen(1:inumlf(indelse(i)), indelse(i));  % index of independent detection in that migration indsw(i)
    if indelse(i) == 4
        temp1 = vertdist(1:numlf(indelse(i)), indelse(i));
        temp2 = temp1(ind4);
        vdistelse = [vdistelse; temp2];
        temp2 = temp1(intersect(iind, ind4));
        ivdistelse = [ivdistelse; temp2];
        temp1 = weight(1:numlf(indelse(i)), indelse(i));
        temp2 = temp1(ind4);
        wtelse = [wtelse; temp2];
        temp2 = temp1(intersect(iind, ind4));
        iwtelse = [iwtelse; temp2];
    elseif indelse(i) == 9
        temp1 = vertdist(1:numlf(indelse(i)), indelse(i));
        temp2 = temp1(ind9);
        vdistelse = [vdistelse; temp2];
        temp2 = temp1(intersect(iind, ind9));
        ivdistelse = [ivdistelse; temp2];
        temp1 = weight(1:numlf(indelse(i)), indelse(i));
        temp2 = temp1(ind9);
        wtelse = [wtelse; temp2];
        temp2 = temp1(intersect(iind, ind9));
        iwtelse = [iwtelse; temp2];
    else
        temp1 = vertdist(1:numlf(indelse(i)), indelse(i));
        vdistelse = [vdistelse; temp1];
        ivdistelse = [ivdistelse; temp1(iind)];
        temp1 = weight(1:numlf(indelse(i)), indelse(i));
        wtelse = [wtelse; temp1];
        iwtelse = [iwtelse; temp1(iind)];
    end
end
pdelse = fitdist(vdistelse,'Normal');  
muelse = pdelse.mu;    % fitted paramaters
sigmaelse = pdelse.sigma;
pdffitelse = pdf(pdelse,-50:0.05:50);


% Change sign of distance to geographical direction
%%% NOTE, 2020/10/27, recall that this 'vertical' distance is LF-HF in the propagation direction,
%%% in which positive means LF is leading, negative is lagging. But this is limited to the
%%% propagation. So if change to absolute geographycal direction, say WSW-ENW, if look at wsw
%%% propagating migrations, if lf leads, then distance should be +, but you need to flip the sign to
%%% convert to geo direction, since lf is in fact at the wsw of hf, - in geo direction. Same if lf
%%% lags. On the other hand, if look at ene group, if lf lags, the raw distance should be -, in
%%% fact, it means lf is at the wsw of hf, still - in the chosen geo coordinates, thus sign of this
%%% group remains unchanged.

%%% currently focus on main LZB region and 002 region only

% main LZB region 
vdistsw = -vdistsw-0.5;     % flip the sign for SW group, check note above
ivdistsw = -ivdistsw-0.5;
vdistne = vdistne+0.15;     % keep the sign for NE group, check note above
ivdistne = ivdistne+0.15;


%% plot the combined migrations, unshifted all
%%% 1. for SW and NE group at main LZB region, i.e. the NW part of the study region
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 4;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.35, 0.85]);

% subplot 1
hold(f.ax(1),'on');
xran = [-15 10];
yran = [-10 10];
indplt = union(indsw,indne);    % index to plot
[f.ax(1),c] = plt_mig_prop_onmap(f.ax(1), indplt, ranvechf98, medallhf, angbest, xran,yran);
% text(f.ax(1),0.95,0.06,'main LZB region','FontSize',11,'unit','normalized','EdgeColor','k',...
%     'Margin',2,'horizontalalignment','right');
% text(f.ax(1),0.5,0.06,'All detections','FontSize',12,'unit','normalized',...
%      'horizontalalignment','center');
text(f.ax(1),0.5,0.12,'Main LZB region','FontSize',12,'unit','normalized',...
     'horizontalalignment','center'); 
text(f.ax(1),0.95,0.92,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2,...
    'horizontalalignment','right');
text(f.ax(1),0.04,0.92,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
hold(f.ax(1),'off');
f.ax(1).Position(2) = f.ax(1).Position(2)-0.05;
% c.Position(2) = c.Position(2)-0.1;

% subplot 2 
binw = 1;

ywid = f.ax(1).Position(4)+f.ax(1).Position(2)-(c.Position(2))-0.022;
% ywid = f.ax(1).Position(4);
set(f.ax(2), 'position', [ 0.55, c.Position(2), 0.35, ywid ]);

meansw = wt_mean(vdistsw,wtsw);
meanne = wt_mean(vdistne,wtne);
% sigmasw = sqrt(wt_var(vdistsw,wtsw,2));
% neff = sum(wtsw);
% conf = 95;
% CIsw = confidence_interval(meansw,sigmasw,neff,conf);
% sigmane = sqrt(wt_var(vdistne,wtne,2));
% neff = sum(wtne);
% % conf = 99;
% CIne = confidence_interval(meanne,sigmane,neff,conf);

hold(f.ax(2),'on');
[f.ax(2), barsw, pdfxlocsw, ~, ~, pdfvalsw] = plt_weighted_dist(f.ax(2), vdistsw, wtsw, binw,'dec');
barsw(1).FaceAlpha = 1;
[f.ax(2), barne, pdfxlocne, ~, ~, pdfvalne] = plt_weighted_dist(f.ax(2), vdistne, wtne, binw,'dec');
barne(1).FaceColor = [0 1 1];
plot(f.ax(2),[meansw meansw],[-100 100],'r--','linew',1.5);    % median of wt dist
plot(f.ax(2),[meanne meanne],[-100 100],'b--','linew',1.5);    % median of wt dist
% errorbar(f.ax(2),meansw,0.26,CIsw(1)-meansw,CIsw(2)-meansw,'horizontal','o','markersize',3,'color',...
%          'r','linew',1,'MarkerEdgeColor','r','MarkerFaceColor','r');
% errorbar(f.ax(2),meanne,0.26,CIne(1)-meanne,CIne(2)-meanne,'horizontal','o','markersize',3,'color',...
%          'b','linew',1,'MarkerEdgeColor','b','MarkerFaceColor','b');
xvect = [-7 -12];
yvect = [0.015 0.015];
drawArrow(f.ax(2),xvect,yvect,f.ax(2).XLim,f.ax(2).YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
% arrow1 = annotation('arrow','linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);
% arrow1.Parent = f.ax(2);
% arrow1.Position = [-7, 0.015, -5, 0] ;
xvect = [7 12];
yvect = [0.015 0.015];
drawArrow(f.ax(2),xvect,yvect,f.ax(2).XLim,f.ax(2).YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
% arrow2 = annotation('arrow','linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);
% arrow2.Parent = f.ax(2);
% arrow2.Position = [7, 0.015, 5, 0] ;
text(f.ax(2),0.05,0.1,strcat({'WSW'}),'fontsize',11,'unit','normalized');
text(f.ax(2),0.95,0.1,strcat({'ENE'}),'fontsize',11,'unit','normalized',...
     'horizontalalignment','right');
text(f.ax(2),0.2,0.6,num2str(length(vdistne)),'fontsize',12,'unit','normalized',...
     'horizontalalignment','center','color','k');
text(f.ax(2),0.2,0.54,strcat({'detections'}),'fontsize',12,'unit','normalized',...
     'horizontalalignment','center','color','k');
text(f.ax(2),0.85,0.6,num2str(length(vdistsw)),'fontsize',12,'unit','normalized',...
     'horizontalalignment','center','color','k');
text(f.ax(2),0.85,0.54,strcat({'detections'}),'fontsize',12,'unit','normalized',...
     'horizontalalignment','center','color','k');
text(f.ax(2),0.04,0.86,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
legend(f.ax(2),[barne,barsw],{'ENE propagation','WSW propagation'},'fontsize',7,...
       'numcolumns',2,...
       'Position',[0.55+0.35*0.048  c.Position(2)+ywid*0.93  0.35*0.9  ywid*0.05]);
hold(f.ax(2),'off');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/temp1',num2str(nfam),'.pdf'));


%% plot the combined migrations, unshifted all
%%% 1. for SW and NE group at main LZB region, i.e. the NW part of the study region
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 6;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol-1
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.35, 0.85]);

% subplot 1
hold(f.ax(1),'on');
xran = [-15 5];
yran = [-15 15];
indplt = union(indsw,indne);    % index to plot
[f.ax(1),c] = plt_mig_prop_onmap(f.ax(1), indplt, ranvechf98, medallhf, angbest, xran,yran);
% text(f.ax(1),0.95,0.06,'main LZB region','FontSize',11,'unit','normalized','EdgeColor','k',...
%     'Margin',2,'horizontalalignment','right');
text(f.ax(1),0.5,0.06,'All detections','FontSize',12,'unit','normalized',...
     'horizontalalignment','center');
text(f.ax(1),0.5,0.12,'Main LZB region,','FontSize',12,'unit','normalized',...
     'horizontalalignment','center'); 
text(f.ax(1),0.04,0.95,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
hold(f.ax(1),'off');
f.ax(1).Position(2) = f.ax(1).Position(2)-0.05;
% c.Position(2) = c.Position(2)-0.1;

% subplot 2 
binw = 1;

ywid = f.ax(1).Position(4)+f.ax(1).Position(2)-(c.Position(2)) - 0.05;
% ywid = f.ax(1).Position(4)-f.ax(1).Position(2)+(c.Position(2)) - 0.06;
set(f.ax(2), 'position', [ 0.55, c.Position(2), 0.35, ywid*2/3 ]);

meansw = wt_mean(vdistsw,wtsw);
meanne = wt_mean(vdistne,wtne);
sigmasw = sqrt(wt_var(vdistsw,wtsw,2));
neff = sum(wtsw);
conf = 95;
CIsw = confidence_interval(meansw,sigmasw,neff,conf);
sigmane = sqrt(wt_var(vdistne,wtne,2));
neff = sum(wtne);
% conf = 99;
CIne = confidence_interval(meanne,sigmane,neff,conf);

hold(f.ax(2),'on');
[f.ax(2), barsw, pdfxlocsw, ~, ~, pdfvalsw] = plt_weighted_dist(f.ax(2), vdistsw, wtsw, binw,'dec');
barsw(1).FaceAlpha = 1;
[f.ax(2), barne, pdfxlocne, ~, ~, pdfvalne] = plt_weighted_dist(f.ax(2), vdistne, wtne, binw,'dec');
barne(1).FaceColor = [0 1 1];
plot(f.ax(2),[meansw meansw],[-100 100],'r--','linew',1.5);    % median of wt dist
plot(f.ax(2),[meanne meanne],[-100 100],'b--','linew',1.5);    % median of wt dist
errorbar(f.ax(2),meansw,0.26,CIsw(1)-meansw,CIsw(2)-meansw,'horizontal','o','markersize',3,'color',...
         'r','linew',1,'MarkerEdgeColor','r','MarkerFaceColor','r');
errorbar(f.ax(2),meanne,0.26,CIne(1)-meanne,CIne(2)-meanne,'horizontal','o','markersize',3,'color',...
         'b','linew',1,'MarkerEdgeColor','b','MarkerFaceColor','b');
xvect = [-7 -12];
yvect = [0.015 0.015];
drawArrow(f.ax(2),xvect,yvect,f.ax(2).XLim,f.ax(2).YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
% arrow1 = annotation('arrow','linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);
% arrow1.Parent = f.ax(2);
% arrow1.Position = [-7, 0.015, -5, 0] ;
xvect = [7 12];
yvect = [0.015 0.015];
drawArrow(f.ax(2),xvect,yvect,f.ax(2).XLim,f.ax(2).YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
% arrow2 = annotation('arrow','linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);
% arrow2.Parent = f.ax(2);
% arrow2.Position = [7, 0.015, 5, 0] ;
text(f.ax(2),0.05,0.1,strcat({'WSW'}),'fontsize',11,'unit','normalized');
text(f.ax(2),0.95,0.1,strcat({'ENE'}),'fontsize',11,'unit','normalized',...
     'horizontalalignment','right');
text(f.ax(2),0.2,0.6,num2str(length(vdistne)),'fontsize',12,'unit','normalized',...
     'horizontalalignment','center','color','k');
text(f.ax(2),0.2,0.54,strcat({'detections'}),'fontsize',12,'unit','normalized',...
     'horizontalalignment','center','color','k');
text(f.ax(2),0.85,0.6,num2str(length(vdistsw)),'fontsize',12,'unit','normalized',...
     'horizontalalignment','center','color','k');
text(f.ax(2),0.85,0.54,strcat({'detections'}),'fontsize',12,'unit','normalized',...
     'horizontalalignment','center','color','k');
% text(f.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
legend(f.ax(2),[barne,barsw],{'ENE propagation','WSW propagation'},'fontsize',7,...
       'numcolumns',2,...
       'Position',[0.55+0.35*0.048  c.Position(2)+ywid*2/3*0.93  0.35*0.9  ywid*2/3*0.05]);
hold(f.ax(2),'off');

set(f.ax(3), 'position', [ 0.55, c.Position(2)+ywid*2/3+0.05, 0.35, ywid*1/3 ]);
ax = f.ax(3);
hold(ax,'on');
patarea = [0 0;
           -100 0;
           -100 1;
            0 1;
            0 0];
patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.1,'edgecolor','none');
plot(ax,[meansw meansw],[-100 100],'r--','linew',2);    % median of wt dist
plot(ax,[meanne meanne],[-100 100],'b--','linew',2);    % median of wt dist
errorbar(ax,meansw,0.26,CIsw(1)-meansw,CIsw(2)-meansw,'horizontal','o','markersize',6,'color',...
         'r','linewidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','CapSize',10);
errorbar(ax,meanne,0.26,CIne(1)-meanne,CIne(2)-meanne,'horizontal','o','markersize',6,'color',...
         'b','linewidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','CapSize',10);
text(ax,0.95,0.85,strcat({'10X'}),'fontsize',10,'unit','normalized','EdgeColor','k',...
     'Margin',2,'horizontalalignment','right');
% text(ax,0.05,0.15,strcat({'99% CI'}),'fontsize',11,'unit','normalized','EdgeColor','k',...
%      'Margin',2,'horizontalalignment','left'); 
text(ax,0.05,0.15,strcat(num2str(conf),{'% CI'}),'fontsize',12,'unit','normalized','horizontalalignment','left'); 
text(ax,0.04,0.85,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
xlimit = f.ax(2).XLim/10;
xlim(ax,xlimit); 
ylim(ax,[0.2 0.3]);
% xticks(xlimit(1): 0.2: xlimit(2));
yticks(ax,[0.2 0.25 0.3]);
% xlabel(ax,'LF-HF (km)','fontsize',11);
ylabel(ax,'PDF estimate','fontsize',11);
ax.Box = 'on';
grid(ax, 'on');
ax.GridLineStyle = '--';
hold(ax,'off');

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
print(f.fig,'-dpdf',strcat(rstpath,'/temp2',num2str(nfam),'.pdf'));


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

f111.fig=figure;
widin = 4;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
set(f111.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 3;
ncol = 1;
for isub = 1:nrow*ncol
    f111.ax(isub) = subplot(nrow,ncol,isub);
end

ax = f111.ax(1);
hold(ax,'on');
box(ax,'on');
mm = statsmaina.meandiff;
ci = statsmaina.ci;
% plot(ax,[meansw meansw],[-100 100],'r--','linew',2);    % median of wt dist
plot(ax,[0 0],[0 2],'k--','linew',1);
errorbar(ax,mm,1,ci(1)-mm,ci(2)-mm,'horizontal','o','markersize',6,'color',...
         'b','linewidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','CapSize',10);
text(ax,0.2,0.85,strcat(num2str(conf),{'% CI'}),'fontsize',12,'unit','normalized','horizontalalignment','left'); 
text(ax,0.99,0.85,strcat({'equal-variance'}),'fontsize',12,'unit','normalized','horizontalalignment','right'); 
yticks(ax,[0 1 2]);
xlim(ax,[-0.1 0.5]);
hold(ax,'off');

ax = f111.ax(2);
hold(ax,'on');
box(ax,'on');
mm = statsmainb.meandiff;
ci = statsmainb.ci;
% plot(ax,[meansw meansw],[-100 100],'r--','linew',2);    % median of wt dist
plot(ax,[0 0],[0 2],'k--','linew',1);
errorbar(ax,mm,1,ci(1)-mm,ci(2)-mm,'horizontal','o','markersize',6,'color',...
         'b','linewidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','CapSize',10);
text(ax,0.2,0.85,strcat(num2str(conf),{'% CI'}),'fontsize',12,'unit','normalized','horizontalalignment','left'); 
text(ax,0.99,0.85,strcat({'equal-variance'}),'fontsize',12,'unit','normalized','horizontalalignment','right'); 
yticks(ax,[0 1 2]);
xlim(ax,[-0.1 0.5]);
hold(ax,'off');

ax = f111.ax(3);
hold(ax,'on');
box(ax,'on');
mm = statsmainc.meandiff;
ci = statsmainc.ci;
% plot(ax,[meansw meansw],[-100 100],'r--','linew',2);    % median of wt dist
plot(ax,[0 0],[0 2],'k--','linew',1);
errorbar(ax,mm,1,ci(1)-mm,ci(2)-mm,'horizontal','o','markersize',6,'color',...
         'b','linewidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','CapSize',10);
text(ax,0.2,0.85,strcat(num2str(conf),{'% CI'}),'fontsize',12,'unit','normalized','horizontalalignment','left'); 
text(ax,0.99,0.85,strcat({'unequal-variance'}),'fontsize',12,'unit','normalized','horizontalalignment','right'); 
yticks(ax,[0 1 2]);
xlim(ax,[-0.1 0.5]);
hold(ax,'off');

%%% save figure
print(f111.fig,'-dpdf',strcat(rstpath,'/temp4',num2str(nfam),'.pdf'));


%% Try the random sampling from data PDF, all detections
%%% for testing code segment, see 'mig_linear_fit_LZB_addfam_Miscellaneous.m'

%%% 1. for SW and NE group at main LZB region, i.e. the NW part of the study region
% compute the mean of the randoms generated and the misfits of PDFs using the above 3 ways
% in order: normpdf, ksdensity, normfit
% NE group
% NOTE: 'meanrandne' tracks the weighted data mean e.g. meanne = wt_mean(vdistne,wtne) in all cases 
N = 1000;   % sampling times
[probdiff,mdiff,refval] = random_sampling_test(N,vdistne,wtne,pdfxlocne,pdfvalne,...
                                    vdistsw,wtsw,pdfxlocsw,pdfvalsw);
                                
                                
%% plot the distribution of generated random numbers
binw = 0.02;

f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 7;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

[f.ax(1),f.ax(2)] = plt_random_mean_simple(f.ax(1),f.ax(2),mdiff(:,1),refval,binw);

[f.ax(3),f.ax(4)] = plt_random_mean_simple(f.ax(3),f.ax(4),mdiff(:,3),refval,binw);

title(f.ax(1),'rand gen from normal fitting', 'fontsize',12);
title(f.ax(3),'rand gen from custom pdf', 'fontsize',12);

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/temp3',num2str(nfam),'.pdf'));

% %%% 2. for SW and NE group at fam 002 region
% [probdiff002,mdiff002,refval002] = random_sampling_test(N,vdist002ne,wt002ne,...
%                                              pdfxloc002ne,pdfval002ne,vdist002sw,wt002sw,...
%                                              pdfxloc002sw,pdfval002sw);
% % plot the distribution of generated random numbers
% binw = 0.1;
% f2 = plt_random_mean_simple(mdiff002,refval002,binw);

%%%%%%%%%%%%%%%%%%%%%% ABONDONED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stack a student's t distribution on the normfitting generated random number distribution
% hold(f11.ax(6), 'on');
% xx = -0.4: 0.001: 0.4;
% stdmean = std(meandiff(:,3), 1);
% synpdf = tpdf(xx, N-1) * stdmean / sqrt(N);
% synt = tpdf(xx, statsmaina.df) * N * 0.02;

% xa = meanrandne(:,3)-meanne;
% xb = meanrandsw(:,3)-meansw;
% synt = xx;  %/ sqrt(var(xa)/N + var(xb)/N);
% syndf = (var(xa)/N + var(xb)/N).^2 / ( (var(xa)/N)^2/(N-1) + (var(xb)/N)^2/(N-1) );
% synpdf = tpdf(synt, syndf);
% plot(xx, synpdf, 'k-', 'linew', 1.5);                            


%% plot the migration direction distribution
theta = angbest;
theta = deg2rad(theta);
figure2
polarhistogram(theta(angbest>0 & angbest<=90),'binw',pi/36,'normalization','count','facec',...
               [0 1 1],'facea',0.8); hold on
polarhistogram(theta(angbest>90 & angbest<=180),'binw',pi/36,'normalization','count','facec',...
               'k','facea',0.8);
polarhistogram(theta(angbest>180 & angbest<=270),'binw',pi/36,'normalization','count','facec',...
               [1 96/255 0],'facea',0.8);
polarhistogram(theta(angbest>270 & angbest<=360),'binw',pi/36,'normalization','count','facec',...
               'w','facea',0.8);
ax = gca;
ax.RTickLabelRotation = 45;
ax.ThetaMinorTick = 'on';
ax.TickLength = [0.02 0];
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.FontSize = 12;
ax.RTick = 0:1:5;
ax.LineWidth = 1.5;
ax.Box = 'on';
text(0,4.6,'N','HorizontalAlignment',"center",'FontSize',14);
text(deg2rad(90),4.6,'E','HorizontalAlignment',"center",'FontSize',14);
text(deg2rad(180),4.6,'S','HorizontalAlignment',"center",'FontSize',14);
text(deg2rad(270),4.6,'W','HorizontalAlignment',"center",'FontSize',14);
text(deg2rad(315),7.4,'b','HorizontalAlignment',"center",'FontSize',12,'EdgeColor','k','Margin',2);
text(deg2rad(45),7,'LZB','HorizontalAlignment',"center",'FontSize',14,'EdgeColor','k','Margin',2);

%%% save figure
print('-dpdf',strcat(rstpath,'/LZB.',num2str(nfam),'mig.direc.distri.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));
clear ax


%% summary of results with error bars from lf fitting
CIFcn = @(x,p)std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * tinv(abs([0,1]-(1-p/100)/2),...
              sum(~isnan(x(:)))-1) + mean(x(:),'omitnan');
f5.fig=figure;
set(f5.fig,'Position',[scrsz(3)/10 2*scrsz(4)/10 2*scrsz(3)/6 0.7*scrsz(4)]);
hold on
plot([0 0],[-100 100],'--','color',[0.5 0.5 0.5]);
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
    if angbest(i)>0 && angbest(i)<=90
        e1 = errorbar(offset(i),ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',6,'color',...
         'k','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[0 1 1]);
    elseif angbest(i)>90 && angbest(i)<=180
        e2 = errorbar(offset(i),ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',6,'color',...
         'k','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k');
    elseif angbest(i)>180 && angbest(i)<=270
        e3 = errorbar(offset(i),ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',6,'color',...
         'k','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[1 96/255 0]);
    end
    text(4.8,ypos,strcat(num2str(hfN),'; ',num2str(lfN)),'fontsize',9,...
        'HorizontalAlignment','right');
    %%% next time use 'units', 'normalized' instead of defaults; also 
end

text(0.02,0.97,'a','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
text(0.9,0.97,'LZB','FontSize',14,'unit','normalized','EdgeColor','k','Margin',2);
text(0.2,0.13,'LF','FontSize',14,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(0.2,0.1,'lags','FontSize',14,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(0.8,0.13,'LF','FontSize',14,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(0.8,0.1,'leads','FontSize',14,'unit','normalized','HorizontalAlignment','center','fontw','bold');
xlabel('Offset between HF and LF (LF-HF) (km)','fontsize',12);
ylabel('Migration number from earlier to later','fontsize',12);
legend([e1, e2, e3],{'NE propagation','SE propagation','SW propagation'},'FontSize',10,...
        'Position',[0.17 0.52 0.2 0.08]);
ylim([0,50]);
xlim([-5,5]);
% yticks(0:1:8);
box on
grid on
set(gca,'GridLineStyle','--')

%%% save figure
print(f5.fig,'-dpdf',strcat(rstpath,'/LZB',num2str(nfam),'.mig.lfit.sum.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));


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
print(f6.fig,'-dpdf',strcat(rstpath,'/LZB',num2str(nfam),'.mig.lfit.hist.',num2str(winlenhf),'_',...
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
    ib1 = find(angbest(ibar)>0 & angbest(ibar)<=90);
    h1 = length(ib1)/length(offset);
    ib2 = find(angbest(ibar)>90 & angbest(ibar)<=180);
    h2 = length(ib2)/length(offset);
    ib3 = find(angbest(ibar)>180 & angbest(ibar)<=270);
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
print(f7.fig,'-dpdf',strcat(rstpath,'/LZB',num2str(nfam),'.mig.lfit.barplot.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));









