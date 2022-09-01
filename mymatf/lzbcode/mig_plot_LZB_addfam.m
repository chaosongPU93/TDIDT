% function mig_plot_LZB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to find the best propagation direction for migration, and then fit the hf
% the migration linearly to get the propagation slope and error distribution.
%
% NOTES:
%   2020/09/08, i am adding a new fam 006, so that i want to check if original
%               migrations are still reasonable, ABORTED!!
%   2020/09/10, i am adding a fam 001 back again for last try, 13 fam in total now
%
% Chao Song, chaosong@princeton.edu
% First created date:   2019/11/21
% Last modified date:   2020/09/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
% close all
% clc

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

scrsz=get(0,'ScreenSize');

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
            '158';      % 158, 20200916,testing purpose
            ];     % 2020/09/07, i am adding a new family 006 to the pool, final version

nfam = size(nfampool,1);
disp(nfam);

winlenhf = 4;

winlenlf = 16;
lofflf = 4;
ccminlf = 0.35;

SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.dcutnodou.',SUFFIXhf);
disp(fname);
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
           -123.967000 48.458667 33.7800;       % 158, 20200916,testing purpose
           ];
       

relacont = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),loccont(2,1),loccont(2,2));
relacont(:,1) = dx;
relacont(:,2) = dy;


% sort according to day, sec
hftime = sortrows(hftime, [13, 15]);
lftime = sortrows(lftime, [13, 15]);

% %%%%%%%%%%%%%%%%%%% FOR TESTING %%%%%%%%%%%%%%%%%%%%%%%
% %%% 3. hf and lf migrations worthwhile to check the projection in propagation direction, maybe some
% %%%     of them don't have many contemporaneous detections
% dateall = unique(hftime(:,12));
% 
% trange = [2005258,3.40e+4,4.10e+4];
% 
% indhf = find(hftime(:,12)==trange(1) & hftime(:,14)>=trange(2) & ...
%                  hftime(:,14)<=trange(3));
% mighf = hftime(indhf,:);
% 
% indlf = find(lftime(:,12)==trange(1) & lftime(:,14)>=trange(2) & ...
%                  lftime(:,14)<=trange(3));
% miglf = lftime(indlf,:);
% 
% %%% define and position the figure frame and axes of each plot
%     f.fig=figure;
%     set(f.fig,'Position',[scrsz(3)/10 0*scrsz(4)/10 2*scrsz(3)/5 4*scrsz(4)/5]);
%     nrow = 2;
%     ncol = 2;    
%     for isub = 1:nrow*ncol
%         f.ax(isub) = subplot(nrow,ncol,isub);
%     end
% 
% %     %%% reposition
% %     xlen = xran(i,2)-xran(i,1);
% %     ylen = yran(i,2)-yran(i,1);
% %     subfxlen = 0.35;    
% %     if xlen >= ylen
% %         subfylen = ylen/xlen*0.4;
% %     else
% %         subfylen = ylen/xlen*0.35;
% %     end
% %     set(f.ax(1), 'position', [ 0.1, 0.5, subfxlen, subfylen]);
% %     set(f.ax(2), 'position', [ 0.55, 0.5, subfxlen, subfylen]);
% %     set(f.ax(3), 'position', [ 0.1, 0.2, 0.35, 0.25]);
% %     set(f.ax(4), 'position', [ 0.55, 0.2, 0.35, 0.25]);
%     
%     % subplot 1 of figure i
%     hold(f.ax(1),'on');
%     dumhf = sortrows(mighf,-14);
%     scatter(f.ax(1),dumhf(:,1),dumhf(:,2), 30, dumhf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
%     colormap(f.ax(1),'jet');
%     c=colorbar(f.ax(1),'SouthOutside');
%     pos = f.ax(1).Position;
%     c.Position = [pos(1), pos(2), pos(3), 0.02];
% %     c.TickLabels=[];
%     c.Label.String = strcat(num2str(trange(1)),' of HF',' (hr)');
%     c.Label.FontSize = 12;
%     caxis(f.ax(1),[trange(2)/3600 trange(3)/3600])
%     plot(f.ax(1),[-100 100],[0 0],'k--');
%     plot(f.ax(1),[0 0],[-100 100],'k--');
%     f.ax(1).Box = 'on';
%     grid(f.ax(1), 'on');
%     axis(f.ax(1), 'equal');
%     xlim(f.ax(1),[-15 15]);
%     ylim(f.ax(1),[-10 10]);
%     f.ax(1).GridLineStyle = '--';
%     f.ax(1).XAxisLocation = 'top';
%     xlabel(f.ax(1),'E (km)','fontsize',12);
%     ylabel(f.ax(1),'N (km)','fontsize',12); 
%     hold(f.ax(1),'off');
% 
%     % subplot 2 of figure i
%     hold(f.ax(2),'on');
%     dumlf = sortrows(miglf,-14);
%     scatter(f.ax(2),dumlf(:,1),dumlf(:,2), 30, dumlf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
%     colormap(f.ax(2),'jet');
%     c=colorbar(f.ax(2),'SouthOutside');
%     pos = f.ax(2).Position;
%     c.Position = [pos(1), pos(2), pos(3), 0.02];
% %     c.TickLabels=[];
%     c.Label.String = strcat(num2str(trange(1)),' of LF',' (hr)');
%     c.Label.FontSize = 12;
%     caxis(f.ax(2),[trange(2)/3600 trange(3)/3600])
%     plot(f.ax(2),[-100 100],[0 0],'k--');
%     plot(f.ax(2),[0 0],[-100 100],'k--');
%     f.ax(2).Box = 'on';
%     grid(f.ax(2), 'on');
%     axis(f.ax(2), 'equal');
%     xlim(f.ax(2),[-15 15]);
%     ylim(f.ax(2),[-10 10]);
%     f.ax(2).GridLineStyle = '--';
%     f.ax(2).XAxisLocation = 'top';
%     xlabel(f.ax(2),'E (km)','fontsize',12);
%     ylabel(f.ax(2),'N (km)','fontsize',12); 
%     hold(f.ax(2),'off');


%% 
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

% trange(11,:) =  [2004198       53000       59000];
% trange(38,:) =  [2005258       35000       36800];
% trange(42,:) =  [2005259       36500       43100];
% trange(40,:) =  [2005259       2500       3450];

% xran = [-30 50];
% yran = [-30 40];
xran = [-20 20];
yran = [-20 20];

indsw = [8,9,12,14,19,29,30,34,36,37,41,46];    % index of SW migrations
indne = [11, 38, 42, 45, 47];    % index of NE migrations
temp = union(indsw,indne);
indswne = union(temp,[22,25]);
indall = 1: 1: size(trange,1);

ind158 = [8,9,11,14,15,18,20,28,31,34,38,40,42,44,45,47];

indplt = ind158;

% for i = 4: length(indplt)
% for i = 9:9
i=16;
    disp(trange(indplt(i),:));
    indhf = find(hftime(:,13)==trange(indplt(i),1) & hftime(:,15)>=trange(indplt(i),2) & ...
                 hftime(:,15)<=trange(indplt(i),3));
    mighf = hftime(indhf,:);

    indlf = find(lftime(:,13)==trange(indplt(i),1) & lftime(:,15)>=trange(indplt(i),2) & ...
                 lftime(:,15)<=trange(indplt(i),3));
    miglf = lftime(indlf,:);

    %%% define and position the figure frame and axes of each plot
    f.fig=figure;
    set(f.fig,'Position',[scrsz(3)/10 0*scrsz(4)/10 2*scrsz(3)/5 4*scrsz(4)/5]);
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
    
    % subplot 1 of figure indplt(i)
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
    juldate = num2str(trange(indplt(i),1));
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
%     c.Label.String = strcat(num2str(trange(indplt(i),1)),' of HF',' (hr)');
    c.Label.FontSize = 11;
    caxis(f.ax(1),[trange(indplt(i),2)/3600 trange(indplt(i),3)/3600])
    text(f.ax(1),0.85,0.15,'HF','FontSize',12,'unit','normalized');
    text(f.ax(1),0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(1).Box = 'on';
    grid(f.ax(1), 'on');
    axis(f.ax(1), 'equal');
    f.ax(1).GridLineStyle = '--';
    f.ax(1).XAxisLocation = 'top';
    medxhf = median(mighf(:,1));
    medyhf = median(mighf(:,2));
%     [rotx, roty] = complex_rot(0,5,-angbest(indplt(i)));
%     xvect = [medxhf-rotx medxhf+rotx];
%     yvect = [medyhf-roty medyhf+roty];
%     drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5);
%     [rotx, roty] = complex_rot(-5,0,-angbest(indplt(i)));
%     xvect = [medxhf-rotx medxhf+rotx];
%     yvect = [medyhf-roty medyhf+roty];
%     drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
    xticks(f.ax(1),xran(1):5:xran(2));
    yticks(f.ax(1),yran(1):5:yran(2));    
    xlabel(f.ax(1),'E (km)','fontsize',11);
    ylabel(f.ax(1),'N (km)','fontsize',11);
    text(f.ax(1),0.5,0.1,num2str(indplt(i)),'FontSize',12,'unit','normalized');
    text(f.ax(1),0.83,0.95,strcat(num2str(size(mighf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(1),0.83,0.90,strcat({'in '},num2str(trange(indplt(i),3)-trange(indplt(i),2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(mighf,1)/(trange(indplt(i),3)-trange(indplt(i),2)));
    text(f.ax(1),0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    xlim(f.ax(1),xran);
    ylim(f.ax(1),yran);
    hold(f.ax(1),'off');

    % subplot 2 of figure indplt(i)
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
%     c.Label.String = strcat(num2str(trange(indplt(i),1)),' of LF',' (hr)');
    c.Label.FontSize = 11;
    caxis(f.ax(2),[trange(indplt(i),2)/3600 trange(indplt(i),3)/3600])
    text(f.ax(2),0.85,0.15,'LF','FontSize',12,'unit','normalized');
    text(f.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(2).Box = 'on';
    grid(f.ax(2), 'on');
    axis(f.ax(2), 'equal');
    f.ax(2).GridLineStyle = '--';
    f.ax(2).XAxisLocation = 'top';
    medxlf = median(miglf(:,1));
    medylf = median(miglf(:,2));
%     [rotx, roty] = complex_rot(0,5,-angbest(indplt(i)));
%     xvect = [medxlf-rotx medxlf+rotx];
%     yvect = [medylf-roty medylf+roty];
%     drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5);
%     [rotx, roty] = complex_rot(-5,0,-angbest(indplt(i)));
%     xvect = [medxlf-rotx medxlf+rotx];
%     yvect = [medylf-roty medylf+roty];
%     drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
    xticks(f.ax(2),xran(1):5:xran(2));
    yticks(f.ax(2),yran(1):5:yran(2));    
    xlabel(f.ax(2),'E (km)','fontsize',11);
    ylabel(f.ax(2),'N (km)','fontsize',11);
    text(f.ax(2),0.5,0.1,num2str(indplt(i)),'FontSize',12,'unit','normalized');
    text(f.ax(2),0.83,0.95,strcat(num2str(size(miglf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(2),0.83,0.9,strcat({'in '},num2str(trange(indplt(i),3)-trange(indplt(i),2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(miglf,1)/(trange(indplt(i),3)-trange(indplt(i),2)));
    text(f.ax(2),0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    xlim(f.ax(2),xran);
    ylim(f.ax(2),yran);
    hold(f.ax(2),'off');

% end








