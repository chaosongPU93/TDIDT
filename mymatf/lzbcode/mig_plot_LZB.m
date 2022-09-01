% function mig_plot_LZB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to find the best propagation direction for migration, and then fit the hf
% the migration linearly to get the propagation slope and error distribution.
%
% NOTES:
%   2020/09/08, i am adding a new fam 006, so that i want to check if original
%               migrations are still reasonable, ABORTED!!
%
% Chao Song, chaosong@princeton.edu
% First created date:   2019/11/21
% Last modified date:   2020/09/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

scrsz=get(0,'ScreenSize');

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/lzbrst');

winlenhf = 4;

winlenlf = 16;
lofflf = 4;
ccminlf = 0.35;

SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.allfam.nodcutnodou.',SUFFIXhf);
hftime = load(fname);
% 16 cols, format is:
%   E(043) N(043) E(own) N(own) lon lat dep off12 off13 off12sec off13sec date main_arrival_time
%   cent_of_win avecc famnum


SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.allfam.nodcutnodou.',SUFFIXlf);
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
hftime = sortrows(hftime, [12, 14]);
lftime = sortrows(lftime, [12, 14]);

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
trange = [2003060,7.34e+4,7.42e+4;
          2003060,8.48e+4,8.64e+4;
          2003061,4.85e+4,4.99e+4;
          2003061,5.65e+4,6.00e+4;
          2003061,6.15e+4,6.30e+4;
          2003062,1.05e+4,1.20e+4;
          2003063,8.52e+4,8.64e+4;
          2004192,7.30e+4,7.60e+4;
          2004193,7.30e+4,7.45e+4;
          2004193,8.25e+4,8.35e+4;
          2004195,2.85e+4,3.40e+4;
          2004195,4.15e+4,4.60e+4;
          2004196,7.55e+4,7.70e+4;
          2004197,0.05e+4,0.35e+4;
          2004197,0.50e+4,0.75e+4;
          2004197,0.80e+4,1.00e+4;
          2004197,3.09e+4,3.18e+4;
          2004197,4.60e+4,4.72e+4;
          2004197,5.05e+4,5.28e+4;
          2004197,8.35e+4,8.64e+4;
          2004198,0.98e+4,1.17e+4;
          2004198,2.00e+4,2.18e+4;
          2004198,4.25e+4,4.37e+4;
          2004198,5.40e+4,5.80e+4;  % good
          2004198,7.40e+4,7.58e+4;
          2004198,8.37e+4,8.58e+4;
          2004199,0.42e+4,0.62e+4;
          2004199,3.20e+4,3.38e+4;
          2004199,3.52e+4,3.75e+4;
          2004199,4.55e+4,5.06e+4;
          2004199,8.10e+4,8.35e+4;
          2004200,1.10e+4,1.30e+4;
          2004200,1.30e+4,1.90e+4;
          2004200,1.90e+4,2.00e+4;
          2004200,2.05e+4,2.20e+4;
          2004200,3.80e+4,4.10e+4;
          2004200,4.70e+4,4.90e+4;
          2004200,8.42e+4,8.60e+4;
          2004201,0.20e+4,0.90e+4;
          2004202,4.45e+4,4.60e+4;
          2004202,4.62e+4,4.80e+4;
          2004203,1.60e+4,1.80e+4;
          2004203,1.85e+4,2.50e+4;
          2004203,4.87e+4,5.58e+4;
          2004203,6.20e+4,6.35e+4;
          2004203,8.25e+4,8.40e+4;
          2004206,1.35e+4,1.50e+4;
          2004206,1.57e+4,1.65e+4;
          2005254,6.80e+4,7.00e+4;
          2005255,2.28e+4,2.40e+4;
          2005255,3.40e+4,3.57e+4;
          2005255,5.12e+4,5.20e+4;
          2005255,5.96e+4,6.30e+4;
          2005255,6.70e+4,6.90e+4;
          2005255,7.40e+4,7.60e+4;
          2005255,8.40e+4,8.50e+4;
          2005256,0.35e+4,0.50e+4;
          2005256,1.18e+4,1.41e+4;
          2005256,1.41e+4,1.65e+4;
          2005256,2.15e+4,2.26e+4;
          2005256,2.93e+4,3.14e+4;
          2005256,6.08e+4,6.38e+4;
          2005256,7.18e+4,7.40e+4;
          2005256,7.60e+4,7.80e+4;
          2005258,3.50e+4,3.66e+4;
          2005258,3.75e+4,4.00e+4;
          2005259,0.00e+4,0.18e+4;
          2005259,0.50e+4,0.77e+4;
          2005259,3.70e+4,3.85e+4;
          2005259,4.00e+4,4.20e+4;
          2005259,7.00e+4,7.20e+4;
          2005259,7.28e+4,7.60e+4;
          2005260,0.10e+4,0.22e+4;
          2005260,0.25e+4,0.35e+4;
          2005260,0.35e+4,0.45e+4;
          2005260,0.45e+4,0.60e+4;
          2005260,1.45e+4,1.65e+4;
          2005260,5.62e+4,5.82e+4;
          2005260,6.80e+4,6.93e+4;
          2005261,0.20e+4,0.40e+4;
          2005261,0.90e+4,1.15e+4;];
            
trange = [
           2004195,2.80e+4,3.20e+4;
           2004195,3.30e+4,3.55e+4;
           2004195,4.15e+4,4.60e+4;
           2004196,7.55e+4,7.70e+4;
           2004197,0.03e+4,0.35e+4;
           2004197,0.50e+4,0.75e+4;
           2004197,4.43e+4,4.60e+4;
           2004197,4.60e+4,4.72e+4;
           2004197,8.35e+4,8.46e+4;
           2004197,8.488e+4,8.52e+4;
           2004197,8.55e+4,8.64e+4;
           2004198,0.58e+4,0.98e+4;
           2004198,0.98e+4,1.20e+4;
           2004198,1.90e+4,2.18e+4;
           2004198,5.40e+4,5.80e+4;
           2004198,7.35e+4,7.60e+4;
           2004198,8.35e+4,8.64e+4;
           2004199,0.50e+4,0.63e+4;
           2004199,0.63e+4,0.80e+4;
           2004199,3.22e+4,3.40e+4;
           2004199,3.55e+4,3.65e+4;
           2004199,4.61e+4,4.91e+4;
           2004199,4.98e+4,5.08e+4;
           2004199,8.10e+4,8.20e+4;
           2004200,1.10e+4,1.28e+4;
           2004200,1.38e+4,1.80e+4;
           2004200,1.85e+4,1.99e+4;
           2004203,1.60e+4,2.20e+4;
           2005254,5.61e+4,5.68e+4;
           2005255,3.42e+4,3.58e+4;
           2005255,4.60e+4,4.66e+4;
           2005255,5.13e+4,5.30e+4;
           2005255,6.12e+4,6.25e+4;
           2005255,6.70e+4,6.85e+4;
           2005255,7.50e+4,7.57e+4;
           2005255,8.42e+4,8.56e+4;
           2005256,0.35e+4,0.50e+4;
           2005256,0.80e+4,0.98e+4;
           2005256,1.90e+4,2.00e+4;
           2005256,2.15e+4,2.23e+4;
           2005256,2.82e+4,2.95e+4;
           2005256,3.25e+4,3.34e+4;
           2005256,3.45e+4,3.53e+4;
           2005256,5.20e+4,5.26e+4;
           2005256,7.24e+4,7.40e+4;
           2005256,7.60e+4,7.80e+4;
           2005257,0.60e+4,0.70e+4;
           2005257,2.10e+4,2.50e+4;
           2005257,3.90e+4,4.40e+4;
           2005257,6.10e+4,6.40e+4;
           2005257,7.36e+4,7.80e+4;
           2005257,7.90e+4,8.30e+4;
           2005258,3.52e+4,3.66e+4;
           2005258,3.80e+4,4.00e+4;
           2005259,0.20e+4,0.40e+4;
           2005259,0.50e+4,0.77e+4;
           2005259,3.70e+4,4.26e+4;
           2005259,7.20e+4,7.58e+4;
           2005260,0.23e+4,0.36e+4;
           2005260,0.36e+4,0.43e+4;
           2005260,0.43e+4,0.57e+4;
           2005260,5.63e+4,5.82e+4;
           2005261,0.82e+4,1.10e+4;];

trange = [ 
           2004197,8.55e+4,8.64e+4;
           2004199,0.63e+4,0.80e+4;
           2005255,3.42e+4,3.58e+4;
           2005256,3.25e+4,3.34e+4;
           2005257,7.90e+4,8.30e+4;
           2005258,3.80e+4,4.00e+4;
           2005259,0.20e+4,0.40e+4;
           2005260,0.23e+4,0.36e+4;       % lf and hf totally distinct
           2005261,0.82e+4,1.10e+4;
            ];
        
trange = [
           2003061,36000,43200;
           2003063,34200,47800;
           2003064,16200,21600];
        
% %%% weird detections
% trange = [2003061,3.90e+4,4.20e+4;
%           2003064,1.80e+4,2.00e+4;
%           2004196,6.90e+4,7.20e+4;
%            2004197,0.80e+4,1.00e+4;
%            2004197,5.10e+4,5.30e+4;
%            ];        
       
xran = [-30 50];
yran = [-30 40];

for i = 1: size(trange)
% for i = 9:9
% i=3;
trange(i,:)
    indhf = find(hftime(:,12)==trange(i,1) & hftime(:,14)>=trange(i,2) & ...
                 hftime(:,14)<=trange(i,3));
    mighf = hftime(indhf,:);

    indlf = find(lftime(:,12)==trange(i,1) & lftime(:,14)>=trange(i,2) & ...
                 lftime(:,14)<=trange(i,3));
    miglf = lftime(indlf,:);

    %%% define and position the figure frame and axes of each plot
    f.fig=figure;
    set(f.fig,'Position',[scrsz(3)/10 0*scrsz(4)/10 2*scrsz(3)/5 4*scrsz(4)/5]);
    nrow = 2;
    ncol = 2;    
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
    end

    %%% reposition
    xlen = xran(2)-xran(1);
    ylen = yran(2)-yran(1);
    subfxlen = 0.35;    
    if xlen >= ylen
        subfylen = ylen/xlen*0.4;
    else
        subfylen = ylen/xlen*0.35;
    end
    set(f.ax(1), 'position', [ 0.1, 0.5, subfxlen, subfylen]);
    set(f.ax(2), 'position', [ 0.55, 0.5, subfxlen, subfylen]);
    set(f.ax(3), 'position', [ 0.1, 0.2, 0.35, 0.25]);
    set(f.ax(4), 'position', [ 0.55, 0.2, 0.35, 0.25]);
    
    % subplot 1 of figure i
    hold(f.ax(1),'on');
    dumhf = sortrows(mighf,-14);
    scatter(f.ax(1),dumhf(:,1),dumhf(:,2), 30, dumhf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(1),'jet');
    c=colorbar(f.ax(1),'SouthOutside');
    pos = f.ax(1).Position;
    c.Position = [pos(1), pos(2), pos(3), 0.02];
%     c.TickLabels=[];
    c.Label.String = strcat(num2str(trange(i,1)),' of HF',' (hr)');
    c.Label.FontSize = 12;
    caxis(f.ax(1),[trange(i,2)/3600 trange(i,3)/3600])
    plot(f.ax(1),[-100 100],[0 0],'k--');
    plot(f.ax(1),[0 0],[-100 100],'k--');
    f.ax(1).Box = 'on';
    grid(f.ax(1), 'on');
    axis(f.ax(1), 'equal');
    xlim(f.ax(1),xran);
    ylim(f.ax(1),yran);
    f.ax(1).GridLineStyle = '--';
    f.ax(1).XAxisLocation = 'top';
    xlabel(f.ax(1),'E (km)','fontsize',12);
    ylabel(f.ax(1),'N (km)','fontsize',12); 
    hold(f.ax(1),'off');

    % subplot 2 of figure i
    hold(f.ax(2),'on');
    dumlf = sortrows(miglf,-14);
    scatter(f.ax(2),dumlf(:,1),dumlf(:,2), 30, dumlf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(2),'jet');
    c=colorbar(f.ax(2),'SouthOutside');
    pos = f.ax(2).Position;
    c.Position = [pos(1), pos(2), pos(3), 0.02];
%     c.TickLabels=[];
    c.Label.String = strcat(num2str(trange(i,1)),' of LF',' (hr)');
    c.Label.FontSize = 12;
    caxis(f.ax(2),[trange(i,2)/3600 trange(i,3)/3600])
    plot(f.ax(2),[-100 100],[0 0],'k--');
    plot(f.ax(2),[0 0],[-100 100],'k--');
    f.ax(2).Box = 'on';
    grid(f.ax(2), 'on');
    axis(f.ax(2), 'equal');
    xlim(f.ax(2),xran);
    ylim(f.ax(2),yran);
    f.ax(2).GridLineStyle = '--';
    f.ax(2).XAxisLocation = 'top';
    xlabel(f.ax(2),'E (km)','fontsize',12);
    ylabel(f.ax(2),'N (km)','fontsize',12); 
    hold(f.ax(2),'off');

end








