% function project_mig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to project all migrations to a specific direction (could be strike, dip, 
% propagation direction or any other directions. And plot them of course 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2019/10/03
% Last modified date:   2019/10/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/forsummary');
datapath = strcat(getenv('ALLAN'),'/data-no-resp/PGCtrio');

% load files
fam = '002';
winlenhf = 4;

winlenlf = 8;
lofflf = 4;
ccminlf = 0.15;

SUFFIXhf = strcat('hf.alldays.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/eventloc.',fam,'.',SUFFIXhf);
hfmaplocall = load(fname);
% 10 cols, format is:
%   lon lat dep off12 off13 off12sec off13sec date strongest_arr_time center_of_win

SUFFIXlf = strcat('lf.alldays.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/eventloc.',fam,'.',SUFFIXlf);
lfmaplocall = load(fname);

% convert absolute loc to relative loc
[ind,~] = ind2sub(size(hfmaplocall), find(hfmaplocall(:,4)==0 & hfmaplocall(:,5)==0, 1, 'first'));
lon0 = hfmaplocall(ind,1);
lat0 = hfmaplocall(ind,2);
if isempty(lon0)
    lon0 = -123.5850;
    lat0 = 48.4367;
end

hfrelalocall = hfmaplocall;
[dx, dy] = absloc2relaloc(hfmaplocall(:,1),hfmaplocall(:,2),lon0,lat0);
hfrelalocall(:,1) = dx;
hfrelalocall(:,2) = dy;

lfrelalocall = lfmaplocall;
[dx, dy] = absloc2relaloc(lfmaplocall(:,1),lfmaplocall(:,2),lon0,lat0);
lfrelalocall(:,1) = dx;
lfrelalocall(:,2) = dy;

%%
%%% 3. hf and lf migrations worthwhile to check the projection in propagation direction, maybe some
%%%    of them don't have many contemporaneous detections, also not all of them need a linear
%%%    fitting.
trange3 = [2004197,8.45e+4,8.55e+4;
           2004198,8.35e+4,8.62e+4;
           2004199,0.20e+4,0.37e+4;
           2005255,3.42e+4,3.65e+4;
           2005255,5.80e+4,5.96e+4;
           2005255,6.15e+4,6.26e+4;
           2005255,6.70e+4,6.90e+4;
           2005255,7.50e+4,7.60e+4;
           2005256,0.36e+4,0.55e+4;
           2005256,7.62e+4,8.30e+4;
           2005256,8.39e+4,8.46e+4;
           2005256,8.47e+4,8.55e+4;];
angle = [280;
         45;
         225;
         100;
         135;
         100;
         135;
         90;
         240;
         120;
         30;
         250];
     
usefulind = [1,4,5,6,7,8,9,10];     
       
for i = 1: length(trange3)
% i=7
    indhf = find(hfrelalocall(:,8)==trange3(i,1) & hfrelalocall(:,10)>=trange3(i,2) & ...
                 hfrelalocall(:,10)<=trange3(i,3));
    mighf = hfrelalocall(indhf,:);
    mighfdum = mighf;
    for j = 1: size(mighf,1)
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angle(i)-90),0,0);
        mighfdum(j,1) = newx;
        mighfdum(j,2) = newy;
    end
    
    indlf = find(lfrelalocall(:,8)==trange3(i,1) & lfrelalocall(:,10)>=trange3(i,2) & ...
                 lfrelalocall(:,10)<=trange3(i,3));
    miglf = lfrelalocall(indlf,:);
    miglfdum = miglf;
    for j = 1: size(miglf,1)
        x0 = miglfdum(j,1);
        y0 = miglfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angle(i)-90),0,0);
        miglfdum(j,1) = newx;
        miglfdum(j,2) = newy;
    end
    
    %%% define and position the figure frame and axes of each plot 
    f.fig=figure(i);
    set(f.fig,'Position',[scrsz(3)/10 scrsz(4)/10 2*scrsz(3)/5 4*scrsz(4)/5]);
    nrow = 2;
    ncol = 2;    
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
    end
    
    %%% reposition
    set(f.ax(1), 'position', [ 0.1, 0.55, 0.35, 0.35]);
    set(f.ax(2), 'position', [ 0.55, 0.55, 0.35, 0.35]);
    set(f.ax(3), 'position', [ 0.1, 0.1, 0.35, 0.35]);
    set(f.ax(4), 'position', [ 0.55, 0.1, 0.35, 0.35]);
    
    % subplot 1 of figure i
    hold(f.ax(1),'on');
    scatter(f.ax(1),mighf(:,1),mighf(:,2), 40, mighf(:,10)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(1),'jet');
    c=colorbar(f.ax(1),'SouthOutside');
    c.Position = [0.1, 0.55, 0.35, 0.02];
    c.TickLabels=[];
    caxis(f.ax(1),[trange3(i,2)/3600 trange3(i,3)/3600])
    plot(f.ax(1),[-100 100],[0 0],'k--');
    plot(f.ax(1),[0 0],[-100 100],'k--');
    
    f.ax(1).Box = 'on';
    grid(f.ax(1), 'on');
    axis(f.ax(1), 'equal');
    f.ax(1).GridLineStyle = '--';
    f.ax(1).XAxisLocation = 'top';
    xran = [round(min(mighf(:,1)))-1 round(max(mighf(:,1)))+1];
    yran = [round(min(mighf(:,2)))-1 round(max(mighf(:,2)))+1];
    [rotx, roty] = complex_rot(0,5,-angle(i));
    xvect = [-rotx rotx];
    yvect = [-roty roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5);    
%     xlim(f.ax(1),xran);
%     ylim(f.ax(1),yran);
    xlabel(f.ax(1),'E (km)','fontsize',12);
    ylabel(f.ax(1),'N (km)','fontsize',12); 
    hold(f.ax(1),'off');

    % subplot 2 of figure i
    hold(f.ax(2),'on');
    scatter(f.ax(2),miglf(:,1),miglf(:,2), 40, miglf(:,10)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(2),'jet');
    c=colorbar(f.ax(2),'SouthOutside');
    c.Position = [0.55, 0.55, 0.35, 0.02];
    c.TickLabels=[];
    caxis(f.ax(2),[trange3(i,2)/3600 trange3(i,3)/3600])
    plot(f.ax(2),[-100 100],[0 0],'k--');
    plot(f.ax(2),[0 0],[-100 100],'k--');
    f.ax(2).Box = 'on';
    grid(f.ax(2), 'on');
    axis(f.ax(2), 'equal');
    f.ax(2).GridLineStyle = '--';
    f.ax(2).XAxisLocation = 'top';
    xran = [round(min(mighf(:,1)))-1 round(max(mighf(:,1)))+1];
    yran = [round(min(mighf(:,2)))-1 round(max(mighf(:,2)))+1];
    xlim(f.ax(2),xran);
    ylim(f.ax(2),yran);
    xlabel(f.ax(2),'E (km)','fontsize',12);
    ylabel(f.ax(2),'N (km)','fontsize',12); 
    hold(f.ax(2),'off');
    
    % subplot 3 of figure i
    hold(f.ax(3),'on');
    alpha(f.ax(3),0.1);
    scatter(f.ax(3),mighfdum(:,10)/3600,mighfdum(:,1),40,[105/255 105/255 105/255],'filled','o');  %, 'MarkerEdgeColor', 'w')    
    scatter(f.ax(3),miglfdum(:,10)/3600,miglfdum(:,1),40,miglfdum(:,11),'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(3),'jet');
    colorbar(f.ax(3));
    f.ax(3).Box = 'on';
    grid(f.ax(3), 'on');
    f.ax(3).GridLineStyle = '--';
    xran = [trange3(i,2)/3600 trange3(i,3)/3600];
    yran = [round(min(mighfdum(:,1)))-1 round(max(mighfdum(:,1)))+1];
    xlim(f.ax(3),xran);
    ylim(f.ax(3),yran);
    xlabel(f.ax(3),'Time (hr)','fontsize',12);
    ylabel(f.ax(3),'Distance along propagation (km)','fontsize',12); 
    hold(f.ax(3),'off');
    
    % subplot 4 of figure i
    hold(f.ax(4),'on');
    alpha(f.ax(4),0.1);
    scatter(f.ax(4),mighfdum(:,10)/3600,mighfdum(:,2),40,[105/255 105/255 105/255],'filled','o');  %, 'MarkerEdgeColor', 'w')
    scatter(f.ax(4),miglfdum(:,10)/3600,miglfdum(:,2),40,miglfdum(:,11),'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(4),'jet');
    colorbar(f.ax(4));
    f.ax(4).Box = 'on';
    grid(f.ax(4), 'on');
    f.ax(4).GridLineStyle = '--';
    xran = [trange3(i,2)/3600 trange3(i,3)/3600];
    yran = [round(min(mighfdum(:,1)))-1 round(max(mighfdum(:,1)))+1];
    xlim(f.ax(4),xran);
    ylim(f.ax(4),yran);
    xlabel(f.ax(4),'Time (hr)','fontsize',12);
    ylabel(f.ax(4),'Distance along orthogonal (km)','fontsize',12);
    hold(f.ax(4),'off');
    
    title(gca,strcat(num2str(trange3(i,1)),'\_',num2str(trange3(i,2)),'-',num2str(trange3(i,3))));
    
    %%% save figure i
    print(f.fig,'-depsc',strcat(rstpath,'/',fam,'projection.contemp_migration.',num2str(trange3(i,1)),...
          '_',num2str(trange3(i,2)),'-',num2str(trange3(i,3)),'.',num2str(winlenhf),'_',...
          num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.timerela.eps'));
    close(f.fig)  
end



















