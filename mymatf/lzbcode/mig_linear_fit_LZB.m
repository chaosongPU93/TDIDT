% function mig_linear_fit_LZB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to use the poropagation direction estimate fit the hf
% the migration linearly to get the propagation slope and error distribution.
%
% Chao Song, chaosong@princeton.edu
% First created date:   2019/11/29
% Last modified date:   2019/11/29
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
            '017'];

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
fname = strcat(rstpath, '/evtloc.allfam.dcutnodou.',SUFFIXhf);
hftime = load(fname);
% 25 cols, format is:
%   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date
%   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum
%


SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.allfam.dcutnodou.',SUFFIXlf);
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
           -123.838500 48.544833 35.6600];
       

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
           2005256,7.24e+4,7.40e+4;
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
resprophf = nan(length(trange)+1,200);
resproplf = nan(length(trange)+1,50);

medallhf = nan(length(trange),4);
medalllf = nan(length(trange),4);
ranvechf = nan(length(trange),2);
ranveclf = nan(length(trange),2);

% create fit object with free constraints
fttpfree = fittype( @(a,b,x) a*x+b);
        
for i = 1: length(trange)
% for i = 1: 21
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
    text(f.ax(1),0.85,0.15,'HF','FontSize',12,'unit','normalized');
    text(f.ax(1),0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
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
    text(f.ax(2),0.85,0.15,'LF','FontSize',12,'unit','normalized');
    text(f.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
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
    % get the information of propagation vector indicating range etc.
    [tmpx,tmpy] = coordinate_rot(medxhf,medyhf,-(angbest(i)-90),0,0);
    medallhf(i,1:4) = [medxhf medyhf tmpx tmpy];
    [tmpx,tmpy] = coordinate_rot(medxlf,medylf,-(angbest(i)-90),0,0);
    medalllf(i,1:4) = [medxlf medylf tmpx tmpy];
    
    ranvechf(i,1:2) = [min(mighfdum(:,1)) max(mighfdum(:,1))];
    ranveclf(i,1:2) = [min(miglfdum(:,1)) max(miglfdum(:,1))];
    
    
    % subplot 4 of figure i
    hold(f.ax(4),'on');
    numhf(i) = length(mighfdum2(:,1));
    numlf(i) = length(miglfdum2(:,1));
    resprophf(i,1:numhf(i)) = outphfprop.residuals;
    resproplf(i,1:numlf(i)) = outplfprop.residuals+(intcptproplf(i)-intcptprophf(i));
    pdhf = fitdist(resprophf(i,1:numhf(i))','Normal');    % fit a distribution of residual, assuming it is normal distributed
    pdlf = fitdist(resproplf(i,1:numlf(i))','Normal');
    muhf = pdhf.mu;    % fitted paramaters
    mulf = pdlf.mu;
    sigmahf = pdhf.sigma;
    sigmalf = pdlf.sigma;
    cihf = paramci(pdhf,'Alpha',0.05);     % give the 95% confidence interval of distribution param
    cilf = paramci(pdlf,'Alpha',0.05);    
    pdffithf = pdf(pdhf,min(resprophf(i,1:numhf(i))):0.05:max(resprophf(i,1:numhf(i))));
    pdffitlf = pdf(pdlf,min(resproplf(i,1:numlf(i))):0.05:max(resproplf(i,1:numlf(i))));
    histogram(f.ax(4),resprophf(i,1:numhf(i)),'binwidth',0.5,'normalization','pdf','facecolor',...
              [1 0 0]);
    histogram(f.ax(4),resproplf(i,1:numlf(i)),'binwidth',0.5,'normalization','pdf','facecolor',...
              [0.6 1 1],'facealpha',0.6);
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
    for j = 1: size(mighfdum2,1)
        medhf(j,1) = median(mighfdum2(1:j,1));
        medhf(j,2) = median(mighfdum2(1:j,2));
        medhf(j,3) = mighfdum2(j,15);
    end
    
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
    
    %%% save figure
    print(f.fig,'-dpdf',strcat(rstpath,'/LZB.mig.proj.lfit.meddist',num2str(trange(i,1)),...
        '_',num2str(trange(i,2)),'-',num2str(trange(i,3)),'.',num2str(winlenhf),'_',...
        num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));
end


%% plot the migration direction in map
%%% define and position the figure frame and axes of each plot
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.35, 0.75]);
set(f.ax(2), 'position', [ 0.55, 0.1, 0.35, 0.75]);


% subplot 1
ax = f.ax(1);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;

for i = 1: size(trange,1)
    [rotx1, roty1] = complex_rot(0, ranvechf(i,2)-medallhf(i,3), -angbest(i));
    [rotx2, roty2] = complex_rot(0, medallhf(i,3)-ranvechf(i,1), -angbest(i));
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
text(ax,0.85,0.15,'HF','FontSize',12,'unit','normalized');
text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

% subplot 2
ax = f.ax(2);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;

for i = 1: size(trange,1)
    [rotx1, roty1] = complex_rot(0, ranveclf(i,2)-medalllf(i,3), -angbest(i));
    [rotx2, roty2] = complex_rot(0, medalllf(i,3)-ranveclf(i,1), -angbest(i));
    xvect = [medalllf(i,1)-rotx2 medalllf(i,1)+rotx1];
    yvect = [medalllf(i,2)-roty2 medalllf(i,2)+roty1];
    drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1,'linestyle','-','color',[0.7 0.7 0.7]);
    
end

for i = 1: size(trange,1)
    scatter(ax,medalllf(i,1),medalllf(i,2), 30, i, 'filled','o');  %, 'MarkerEdgeColor', 'w')
end
colormap(ax,'jet');
c=colorbar(ax,'SouthOutside');
% pos = ax.Position;
% c.Position = [pos(1), pos(2)-0.03, pos(3), 0.02];
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
text(ax,0.85,0.15,'LF','FontSize',12,'unit','normalized');
text(ax,0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

for i = 1: size(trange,1)
    if angbest(i)>0 && angbest(i)<=90
        disp(i)
    end
end

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/prop_direction_map.pdf'));
% print(f.fig,'-dpdf',strcat(rstpath,'/unnamed.',num2str(winlenhf),'_',...
%     num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));
    

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
print('-dpdf',strcat(rstpath,'/LZB.mig.direc.distri.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));



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
print(f5.fig,'-dpdf',strcat(rstpath,'/LZB.mig.lfit.sum.',num2str(winlenhf),'_',...
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
print(f6.fig,'-dpdf',strcat(rstpath,'/LZB.mig.lfit.hist.',num2str(winlenhf),'_',...
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
print(f7.fig,'-dpdf',strcat(rstpath,'/LZB.mig.lfit.barplot.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% combine the above three figures into one for the paper
CIFcn = @(x,p)std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * tinv(abs([0,1]-(1-p/100)/2),...
              sum(~isnan(x(:)))-1) + mean(x(:),'omitnan');
          
[scrsz, res] = pixelperinch(1);
f.fig=figure;
f.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 7;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol-1
    f.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.45, 0.85]);
set(f.ax(2), 'position', [ 0.62, 0.1, 0.35, 0.4]);
set(f.ax(3), 'position', [ 0.62, 0.55, 0.35, 0.4]);


ax = f.ax(1);
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
    if angbest(i)>0 && angbest(i)<=90
        e1 = errorbar(ax,offset(i),ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[0 1 1]);
    elseif angbest(i)>90 && angbest(i)<=180
        e2 = errorbar(ax,offset(i),ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor','k');
    elseif angbest(i)>180 && angbest(i)<=270
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


ax = f.ax(3);
theta = angbest;
theta = deg2rad(theta);
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
ax.FontSize = 9;
ax.RTick = 0:1:5;
% ax.LineWidth = 1.5;
ax.Box = 'on';
ax.GridAlpha = 0.3;
text(ax,deg2rad(60),3.5,num2str(sum(angbest>0 & angbest<=90)),'HorizontalAlignment','center',...
     'FontSize',11);           
text(ax,deg2rad(130),3.5,num2str(sum(angbest>90 & angbest<=180)),'HorizontalAlignment','center',...
     'FontSize',11);           
text(ax,deg2rad(225),3.5,num2str(sum(angbest>180 & angbest<=270)),'HorizontalAlignment','center',...
     'FontSize',11); 
text(ax,0,4.6,'N','HorizontalAlignment',"center",'FontSize',12);
text(ax,deg2rad(90),4.6,'E','HorizontalAlignment',"center",'FontSize',12);
text(ax,deg2rad(180),4.6,'S','HorizontalAlignment',"center",'FontSize',12);
text(ax,deg2rad(270),4.6,'W','HorizontalAlignment',"center",'FontSize',12);
text(ax,deg2rad(308),7.3,'b','HorizontalAlignment',"center",'FontSize',11,'EdgeColor','k','Margin',2);
% text(ax,deg2rad(45),7,'LZB','HorizontalAlignment',"center",'FontSize',12,'EdgeColor','k','Margin',2);
hold(ax,'off');

ax = f.ax(2);
hold(ax,'on');
% xloc = [-1.75 -1.25 -0.75 -0.25 0.25 0.75 1.25 3.25];
xloc = -3.75:0.5:3.75;
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


print(f.fig,'-dpdf',strcat(rstpath,'/LZB.sum.mig.lfit.',num2str(winlenhf),'_',...
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
    if angbest(i)>0 && angbest(i)<=90
        e1 = errorbar(ax,offset(i),ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[0 1 1]);
    elseif angbest(i)>90 && angbest(i)<=180
        e2 = errorbar(ax,offset(i),ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor','k');
    elseif angbest(i)>180 && angbest(i)<=270
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

ax = f22.ax(5);
theta = angbest;
theta = deg2rad(theta);
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
ax.FontSize = 9;
ax.RTick = 0:1:5;
% ax.LineWidth = 1.5;
ax.Box = 'on';
ax.GridAlpha = 0.3;
text(ax,deg2rad(60),3.5,num2str(sum(angbest>0 & angbest<=90)),'HorizontalAlignment','center',...
     'FontSize',11);           
text(ax,deg2rad(130),3.5,num2str(sum(angbest>90 & angbest<=180)),'HorizontalAlignment','center',...
     'FontSize',11);           
text(ax,deg2rad(225),3.5,num2str(sum(angbest>180 & angbest<=270)),'HorizontalAlignment','center',...
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
    ib1 = find(angbest(ibar)>0 & angbest(ibar)<=90);
    h1 = length(ib1)/length(offset);
    ib2 = find(angbest(ibar)>90 & angbest(ibar)<=180);
    h2 = length(ib2)/length(offset);
    ib3 = find(angbest(ibar)>180 & angbest(ibar)<=270);
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
offset1 = offset(angbest>0 & angbest<=90);
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
offset2 = offset(angbest>90 & angbest<=180);
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
offset3 = offset(angbest>180 & angbest<=270);
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

print(f22.fig,'-dpdf',strcat(rstpath,'/LZB.sum.v2.mig.lfit.',num2str(winlenhf),'_',...
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
    
%     %%% reposition
%     set(f8.ax(1), 'position', [ 0.1, 0.64, 0.35, 0.35]);
%     set(f8.ax(2), 'position', [ 0.55, 0.64, 0.35, 0.35]);
%     set(f8.ax(3), 'position', [ 0.1, 0.39, 0.35, 0.20]);
%     set(f8.ax(4), 'position', [ 0.55, 0.39, 0.35, 0.20]);
%     set(f8.ax(5), 'position', [ 0.1, 0.14, 0.35, 0.20]);
%     set(f8.ax(6), 'position', [ 0.55, 0.14, 0.35, 0.20]);
    
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
    if i == 3
        c.Label.String = strcat({'Time (hr) on 14 Sep. 2005'});
    else
        c.Label.String = strcat({'Time (hr) on 16 Jul. 2004'});
    end
    c.Label.FontSize = 11;
    caxis(f8.ax(1),[trange(ind(i),2)/3600 trange(ind(i),3)/3600])
    text(f8.ax(1),0.85,0.15,'HF','FontSize',12,'unit','normalized');
    text(f8.ax(1),0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
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
    if i == 3
        c.Label.String = strcat({'Time (hr) on 14 Sep. 2005'});
    else
        c.Label.String = strcat({'Time (hr) on 16 Jul. 2004'});
    end
    c.Label.FontSize = 11;
    caxis(f8.ax(2),[trange(ind(i),2)/3600 trange(ind(i),3)/3600])
    text(f8.ax(2),0.85,0.15,'LF','FontSize',12,'unit','normalized');
    text(f8.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
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
    hold(f8.ax(2),'off');
    
    % subplot 3 of figure ind(i)
    hold(f8.ax(3),'on');
    f8.ax(3).FontSize = 9;
    scatter(f8.ax(3),mighfdum(:,15)/3600,mighfdum(:,1),20,[1 0 0],'filled','o',...
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
    
%     mighfdum2 = mighfdum;
%     miglfdum2 = miglfdum;

    [fitobjhfprop,gofhfprop,outphfprop] = fit(mighfdum2(:,15)/3600, mighfdum2(:,1),fttpfree,'Robust',...
                                              'Bisquare','StartPoint',[1 1]);
    % output fit parameters
    coefprop = coeffvalues(fitobjhfprop);
    slopeprophf(ind(i)) = coefprop(1);
    intcptprophf(ind(i)) = coefprop(2);
    fitprophf = feval(fitobjhfprop,mighfdum2(:,15)/3600);
    plot(f8.ax(3),mighfdum2(:,15)/3600,fitprophf,'-','linewidth',2,'color',[1 0 0]);
        
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
    
    % subplot 4 of figure ind(i)
    hold(f8.ax(4),'on');
    f8.ax(4).FontSize = 9;
    numhf(ind(i)) = length(mighfdum2(:,1));
    numlf(ind(i)) = length(miglfdum2(:,1));
    resprophf(ind(i),1:numhf(ind(i))) = outphfprop.residuals;
    resproplf(ind(i),1:numlf(ind(i))) = outplfprop.residuals+(intcptproplf(ind(i))-intcptprophf(ind(i)));
    pdhf = fitdist(resprophf(ind(i),1:numhf(ind(i)))','Normal');    % fit a distribution of residual, assuming it is normal distributed
    pdlf = fitdist(resproplf(ind(i),1:numlf(ind(i)))','Normal');
    muhf = pdhf.mu;    % fitted paramaters
    mulf = pdlf.mu;
    sigmahf = pdhf.sigma;
    sigmalf = pdlf.sigma;
    cihf = paramci(pdhf,'Alpha',0.05);     % give the 95% confidence interval of distribution param
    cilf = paramci(pdlf,'Alpha',0.05);    
    pdffithf = pdf(pdhf,min(resprophf(ind(i),1:numhf(ind(i)))):0.05:max(resprophf(ind(i),1:numhf(ind(i)))));
    pdffitlf = pdf(pdlf,min(resproplf(ind(i),1:numlf(ind(i)))):0.05:max(resproplf(ind(i),1:numlf(ind(i)))));
    histogram(f8.ax(4),resprophf(ind(i),1:numhf(ind(i))),'binwidth',0.5,'normalization','pdf','facecolor',...
              [1 0 0]);
    histogram(f8.ax(4),resproplf(ind(i),1:numlf(ind(i))),'binwidth',0.5,'normalization','pdf','facecolor',...
              [0.6 1 1],'facealpha',0.6);
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
    for j = 1: size(mighfdum2,1)
        medallhf(j,1) = median(mighfdum2(1:j,1));
        medallhf(j,2) = median(mighfdum2(1:j,2));
        medallhf(j,3) = mighfdum2(j,15);
    end
    
    for j = 1: size(miglfdum2,1)
        medalllf(j,1) = median(miglfdum2(1:j,1));
        medalllf(j,2) = median(miglfdum2(1:j,2));
        medalllf(j,3) = miglfdum2(j,15);
    end
       
    
    % subplot 5 of figure ind(i)
    hold(f8.ax(5),'on');
    f8.ax(5).FontSize = 9;
    scatter(f8.ax(5),medallhf(:,3)/3600,medallhf(:,1),20,[1 0 0],'filled','o','MarkerEdgeColor',...
            'k');
    scatter(f8.ax(5),medalllf(:,3)/3600,medalllf(:,1),20,[0.6 1 1],'filled','o','MarkerEdgeColor',...
            'k');            
    f8.ax(5).Box = 'on';
    grid(f8.ax(5), 'on');
    f8.ax(5).GridLineStyle = '--';
    xran1 = [trange(ind(i),2)/3600 trange(ind(i),3)/3600];
    yran1 = [round(min([medallhf(:,1);medalllf(:,1)]))-1 ...
             round(max([medallhf(:,1);medalllf(:,1)]))+1];
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
    scatter(f8.ax(6),medallhf(:,3)/3600,medallhf(:,2),20,[1 0 0],'filled','o','MarkerEdgeColor',...
            'k');
    scatter(f8.ax(6),medalllf(:,3)/3600,medalllf(:,2),20,[0.6 1 1],'filled','o','MarkerEdgeColor',...
            'k');          
    f8.ax(6).Box = 'on';
    grid(f8.ax(6), 'on');
    f8.ax(6).GridLineStyle = '--';
    xran1 = [trange(ind(i),2)/3600 trange(ind(i),3)/3600];
    yran1 = [round(min([medallhf(:,2);medalllf(:,2)]))-1 ...
             round(max([medallhf(:,2);medalllf(:,2)]))+1];
    xlim(f8.ax(6),xran1);
    ylim(f8.ax(6),yran1);
    text(f8.ax(6),0.97,0.1,'f','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2,...
         'HorizontalAlignment','right');
    xlabel(f8.ax(6),'Time (hr)','fontsize',11);
    ylabel(f8.ax(6),'Med. dist. along ort. (km)','fontsize',11);
    hold(f8.ax(6),'off');
    
    %%% save figure
    print(f8.fig,'-dpdf',strcat(rstpath,'/LZB.sel.mig.proj.lfit.meddist',num2str(trange(ind(i),1)),...
        '_',num2str(trange(ind(i),2)),'-',num2str(trange(ind(i),3)),'.',num2str(winlenhf),'_',...
        num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));
end







