% function plot_contemp_mig_rela
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to plot the locations of contemporal hf and lf
% migrations inverted by hypoinverse in km relative to offset (0,0)
% 
% Looks like figure 4b in Rubin&Armbruster (2013)
% 
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/08/21
% Last modified date:   2019/08/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

scrsz=get(0,'ScreenSize');

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/forsummary');
datapath = strcat(getenv('ALLAN'),'/data-no-resp/PGCtrio');

% load files
fam = '002';
winlenhf = 4;

winlenlf = 16;
lofflf = 4;
ccminlf = 0.5;


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

SUFFIXhf = strcat('hf.alldays.count.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/eventloc.',fam,'.',SUFFIXhf);
hfmapcount = load(fname);
% 8 cols, format is:
%   lon lat dep off12 off13 off12sec off13sec count

hfrelacount = hfmapcount;
[dx, dy] = absloc2relaloc(hfmapcount(:,1),hfmapcount(:,2),lon0,lat0);

hfrelacount(:,1) = dx;
hfrelacount(:,2) = dy;


      
%% plot the selected contemporal migrations 
trange = [2004199,0.20e+4,0.37e+4;
          2005255,3.42e+4,3.65e+4;
          2005255,5.80e+4,5.96e+4;
          2005255,6.15e+4,6.26e+4;
          2005255,6.70e+4,6.90e+4;
          2005255,7.50e+4,7.60e+4;
          2005256,0.36e+4,0.55e+4;
          2005256,7.70e+4,8.55e+4;];          
meddif = zeros(length(trange), 4);    
    
for i = 1: length(trange)
% i=3;
    %%% load contemporal migration results and convert abs map loc to
    %%% rela loc
    fname = strcat(rstpath,'/eventloc.contemp_mig_hf.',fam,'_',num2str(trange(i,1)),...
                   '_',num2str(trange(i,2)),'-',num2str(trange(i,3)),'.',num2str(winlenhf),'_',...
                   num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf));
    conmigmaphf = load(fname);
    conmigrelahf = conmigmaphf;
    [dx, dy] = absloc2relaloc(conmigmaphf(:,1),conmigmaphf(:,2),lon0,lat0);
    conmigrelahf(:,1) = dx;
    conmigrelahf(:,2) = dy;
    conmigrelahf(:,end) = conmigrelahf(:,end)/3600;
       
    fname = strcat(rstpath,'/eventloc.contemp_mig_lf.',fam,'_',num2str(trange(i,1)),...
                   '_',num2str(trange(i,2)),'-',num2str(trange(i,3)),'.',num2str(winlenhf),'_',...
                   num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf));
    conmigmaplf = load(fname);
    conmigrelalf = conmigmaplf;
    [dx, dy] = absloc2relaloc(conmigmaplf(:,1),conmigmaplf(:,2),lon0,lat0);
    conmigrelalf(:,1) = dx;
    conmigrelalf(:,2) = dy;    
    conmigrelalf(:,end) = conmigrelalf(:,end)/3600;
    
    
    %%% calculate the magnitude (length) of difference
    mag = sqrt( (conmigrelahf(:,1)-conmigrelalf(:,1)).^2 + (conmigrelahf(:,2)-conmigrelalf(:,2)).^2 );
    meanmag = mean(mag);
    medmag = median(mag);
    
    %%% remember to examine the results, compared to hf_lf_temopral
    hflfdif = [];
    hflfdif(:,1:2) = conmigrelahf(:,1:2) - conmigrelalf(:,1:2);
    hflfdif(:,3:4) = conmigrelahf(:,4:5) - conmigrelalf(:,4:5)/20*40;   % offset in lf is based on 20 sps
    meddif(i,:) = median(hflfdif);
    
    if(min(conmigmaplf(:,10))< trange(i,2) || min(conmigmaphf(:,10)) <trange(i,2) ...
            ||  max(conmigmaplf(:,10)) >trange(i,3) || max(conmigmaphf(:,10)) >trange(i,3))
        disp('wrong');
    else
        disp('right');
    end

    
    %%% load hf seismic trace file in that day.
    datestr = num2str(trange(i,1));
    yr = datestr(1:4);
    day = datestr(5:end);
    sps=40;
    winlen=winlenhf*sps;      % length in smaples
    hi=6.5;
    lo=1.25;
    loff = 2.1;
    ccmin = 0.44;
    npo = 2;
    npa = 2;
    mshift = 29;
    nstanew = 4;
    IDENTIF = strcat(yr,'.',day,'.',fam,'.loff',num2str(loff),'.ccmin',num2str(ccmin),'.nponpa', ...
                     num2str(npo),num2str(npa),'.ms',num2str(mshift));
    oritracehf = load(strcat(datapath, '/MAPS/seistraceoneday_',IDENTIF,'_',num2str(lo),'-', ...
                    num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps'));   % original trace
    
    dettracehf = load(strcat(datapath,'/MAPS/seistraceall_',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                    num2str(winlen/sps),'s',num2str(sps),'sps',num2str(nstanew), 'newsta'));    % detected windows
    
    %%% load lf seismic trace file in that day as well
    sps=20;
    winlen=winlenlf*sps;      % length in smaples
    hi=1.25;
    lo=0.5;
    loff = lofflf;
    ccmin = ccminlf;
    npo = 2;
    npa = 2;
    mshift = 39;
    nstanew = 4;
    IDENTIF = strcat(yr,'.',day,'.',fam,'.loff',num2str(loff),'.ccmin',num2str(ccmin),'.nponpa', ...
                     num2str(npo),num2str(npa),'.ms',num2str(mshift));
    oritracelf = load(strcat(datapath, '/MAPS/seistraceoneday_',IDENTIF,'_',num2str(lo),'-', ...
                    num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps'));   % original trace
    
    dettracelf = load(strcat(datapath,'/MAPS/seistraceall_',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                    num2str(winlen/sps),'s',num2str(sps),'sps',num2str(nstanew), 'newsta'));
                
    %%% define and position the figure frame and axes of each plot 
    f.fig=figure(i);
    set(f.fig,'Position',[scrsz(3)/10 scrsz(4)/10 2*scrsz(3)/4 4*scrsz(4)/5]);
    nrow = 3;
    ncol = 2;    
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
    end
    
    %%% reposition
    fixht = 0.06;
    set(f.ax(1), 'position', [ (1-0.9)/ncol, 0.5, 1/ncol*0.8 0.4]);
    set(f.ax(2), 'position', [ (2-0.9)/ncol, 0.5, 1/ncol*0.8 0.4]);
    set(f.ax(3), 'position', [ (1-0.9)/ncol, 0.5-fixht, 1/ncol*0.8 fixht]);
    set(f.ax(4), 'position', [ (2-0.9)/ncol, 0.5-fixht, 1/ncol*0.8 fixht]);
    set(f.ax(5), 'position', [ (1-0.9)/ncol, 0.5-2*fixht, 1/ncol*0.8 fixht]);
    set(f.ax(6), 'position', [ (2-0.9)/ncol, 0.5-2*fixht, 1/ncol*0.8 fixht]);

    % subplot 1 of figure i
    hold(f.ax(1),'on');
%     if i = 8
%         ind = find(conmigrelahf(:,end)<=83000/3600);
%         conmigrelahf = conmigrelahf(ind,:);
%         conmigrelahf = 
    scatter(f.ax(1),conmigrelahf(:,1),conmigrelahf(:,2), 40, conmigrelahf(:,end), 'filled','o', 'MarkerEdgeColor', 'k');   %, 'MarkerEdgeColor', 'w')
    scatter(f.ax(1),conmigrelalf(:,1),conmigrelalf(:,2), 60, conmigrelalf(:,end), 'filled','s', 'MarkerEdgeColor', 'k');   %, 'MarkerEdgeColor', 'w')       
    colormap(f.ax(1),'jet');
    c=colorbar(f.ax(1),'SouthOutside');
    c.Position = [(1-0.9)/ncol+1/ncol*0.8/6, 0.5, 1/ncol*0.8*4/6, 0.02];
    c.TickLabels=[];
    caxis(f.ax(1),[trange(i,2)/3600 trange(i,3)/3600]);
    plot(f.ax(1),[-100 100],[0 0],'k--'); 
    plot(f.ax(1),[0 0],[-100 100],'k--'); 
    f.ax(1).Box = 'on';
    grid(f.ax(1), 'on');
    axis(f.ax(1), 'equal');
    f.ax(1).GridLineStyle = '--';
    xlabel(f.ax(1),'E (km)','fontsize',12);
    ylabel(f.ax(1),'N (km)','fontsize',12);
    f.ax(1).XAxisLocation = 'top';
    xran = [-8 4];
    yran = [-4 6];
    xlen = xran(2)-xran(1);
    ylen = yran(2)-yran(1);
    for j = 1: length(conmigrelahf(:,1))
        x = [conmigrelahf(j,1) conmigrelalf(j,1)];
        y = [conmigrelahf(j,2) conmigrelalf(j,2)];        
        drawArrow(f.ax(1),x,y,xran,yran,'color',[128/255 128/255 128/255]);    
    end
    x = [0 -meddif(i,1)];
    y = [0 -meddif(i,2)];
    text(f.ax(1),xran(1)+0.15*xlen, yran(1)+0.8*ylen, strcat({'median= '},num2str(-meddif(i,1)),{', '},num2str(-meddif(i,2))),'fontsize',12);
    drawArrow(f.ax(1),x,y,xran,yran,'linewidth',1.5);   % draw the median difference from hf to lf 
    hold(f.ax(1),'off');

%    %%%%%%%%%%%%% could use arrow, but not very good %%%%%%%%%%%%%%%%%%%%%%
%    %  ax=gca;
%    %  axis(ax);       % must do this before call arrow and after changing axis    
%    %  arrow(conmigrelahf(:,1:2),
%    %  conmigrelalf(:,1:2),7,'BaseAngle',25,'TipAngle',30); hold on  % way 1

%    %  arrow(conmigrelahf(:,1:2),conmigrelalf(:,1:2),7,'Type','line','Color','k'); hold on     % way 2
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    %%%%%%%%%%%%% could use quiver, but options is limiteed %%%%%%%%%%%%%%%
%    % axis([-8 4 -4 6])
%    % df = conmigrelalf(:,1:2) - conmigrelahf(:,1:2);
%    % quiver(ax,conmigrelahf(:,1), conmigrelahf(:,2), df(:,1), df(:,2), 0, 'k-');     % way 3
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    % subplot 2 of figure i 
    hold(f.ax(2),'on');
    dumhf = hfrelacount;
    dumhf(dumhf(:,8)==1, :) = [];
    scatter(f.ax(2),dumhf(:,1),dumhf(:,2), 2,'ko');
    ind = find(hfrelalocall(:,8)==trange(i,1) & hfrelalocall(:,10)>=trange(i,2) & hfrelalocall(:,10)<=trange(i,3));
    mighf = sortrows(hfrelalocall(ind,:),-10);
    scatter(f.ax(2),mighf(:,1),mighf(:,2), 40, mighf(:,end)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(2),'jet');
    c=colorbar(f.ax(2),'SouthOutside');
    c.Position = [(2-0.9)/ncol+1/ncol*0.8/6, 0.5, 1/ncol*0.8*4/6, 0.02];
    c.TickLabels=[];
    caxis(f.ax(2),[trange(i,2)/3600 trange(i,3)/3600])
    plot(f.ax(2),[-100 100],[0 0],'k--');
    plot(f.ax(2),[0 0],[-100 100],'k--');
    f.ax(2).Box = 'on';
    grid(f.ax(2), 'on');
    axis(f.ax(2), 'equal');
    f.ax(2).GridLineStyle = '--';
    f.ax(2).XAxisLocation = 'top';
    xran = [-8 4];
    yran = [-4 6];
    xlim(f.ax(2),xran);
    ylim(f.ax(2),yran);
    xlabel(f.ax(2),'E (km)','fontsize',12);
    ylabel(f.ax(2),'N (km)','fontsize',12);   
    hold(f.ax(2),'off');
    
    % subplot 3 of figure i
    hold(f.ax(3),'on');
    f.ax(3).Box = 'on';
    orithf = oritracehf(:,1:2);    % only use time and 1st station
    orithf(:,1) = orithf(:,1)/3600;     % convert time in sec to hr  
    scalehf = max([orithf(100:end-100,2); -orithf(100:end-100,2)]);
    orithf(:,2) = orithf(:,2)/scalehf;     % normalize to -1,1  
    conthf = [];
    for ii = 1: size(conmigrelahf,1)
        indtemp = find(orithf(:,1)>=conmigrelahf(ii,end)-winlenhf/2/3600 & orithf(:,1)<=conmigrelahf(ii,end)+winlenhf/2/3600);
        conthf = [conthf; orithf(indtemp,:)];   % contemporaneous hf trace
    end        
    plot(f.ax(3),orithf(:,1),orithf(:,2),'-','color',[160/255 160/255 160/255],'linewidth',0.5);
    plot(f.ax(3),conthf(:,1),conthf(:,2),'r.');    
    yran = [-1.5 1.5];
    ylim(f.ax(3),yran);
    tlen = (trange(i,3)-trange(i,2))/3600;
    xran = [trange(i,2)/3600-0.25*tlen trange(i,3)/3600+0.25*tlen];
    xlim(f.ax(3),xran);
    plot(f.ax(3),[trange(i,2)/3600 trange(i,2)/3600], yran, 'k--');    
    plot(f.ax(3),[trange(i,3)/3600 trange(i,3)/3600], yran, 'k--');    %,'linewidth',0.5
    xlabel(f.ax(3),strcat({'Hours on '},yr,{' '},day),'fontsize',12);
    f.ax(3).XTickLabel = {};
    hold(f.ax(3),'off');
    
    % subplot 4 of figure i
    hold(f.ax(4),'on');
    f.ax(4).Box = 'on';
    detthf = dettracehf(:,1:2);
    detthf(:,1) = detthf(:,1)/3600;     % convert time in sec to hr  
    detthf(:,2) = detthf(:,2)/scalehf;    % use the same scale to norm
    plot(f.ax(4),orithf(:,1),orithf(:,2),'-','color',[160/255 160/255 160/255],'linewidth',0.5);
    plot(f.ax(4),detthf(:,1),detthf(:,2),'r.');
    ylim(f.ax(4),yran);
    xlim(f.ax(4),xran);
    plot(f.ax(4),[trange(i,2)/3600 trange(i,2)/3600], yran, 'k--');    
    plot(f.ax(4),[trange(i,3)/3600 trange(i,3)/3600], yran, 'k--');    %,'linewidth',0.5
    f.ax(4).XTickLabel={};
    hold(f.ax(4),'off');
    
    % subplot 5 of figure i
    hold(f.ax(5),'on');
    f.ax(5).Box = 'on';
    oritlf = oritracelf(:,1:2);    % only use time and 1st station
    oritlf(:,1) = oritlf(:,1)/3600;     % convert time in sec to hr  
%     scalelf = max([oritlf(100:end-100,2); -oritlf(100:end-100,2)]);
    oritlf(:,2) = oritlf(:,2)/scalehf;     % normalize to -1,1
    contlf = [];
    for ii = 1: size(conmigrelalf,1)
        indtemp = find(oritlf(:,1)>=conmigrelalf(ii,end)-winlenlf/2/3600 & oritlf(:,1) <=conmigrelalf(ii,end)+winlenlf/2/3600);
        contlf = [contlf; oritlf(indtemp,:)];   % contemporaneousLf trace
    end    
    plot(f.ax(5),oritlf(:,1),oritlf(:,2),'-','color',[160/255 160/255 160/255],'linewidth',1);
    plot(f.ax(5),contlf(:,1),contlf(:,2),'r.'); 
    ylim(f.ax(5),yran);
    xlim(f.ax(5),xran);
    plot(f.ax(5),[trange(i,2)/3600 trange(i,2)/3600], yran, 'k--');    
    plot(f.ax(5),[trange(i,3)/3600 trange(i,3)/3600], yran, 'k--');    %,'linewidth',0.5
    xlabel(f.ax(5),strcat({'Hours on '},yr,{' '},day),'fontsize',12);
    hold(f.ax(5),'off');
    
    % subplot 6 of figure i
    hold(f.ax(6),'on');
    f.ax(6).Box = 'on';
    dettlf = dettracelf(:,1:2);
    dettlf(:,1) = dettlf(:,1)/3600;     % convert time in sec to hr  
    dettlf(:,2) = dettlf(:,2)/scalehf;    % use the same scale to norm
    plot(f.ax(6),oritlf(:,1),oritlf(:,2),'-','color',[160/255 160/255 160/255],'linewidth',1);
    plot(f.ax(6),dettlf(:,1),dettlf(:,2),'r.');
    ylim(f.ax(6),yran);
    xlim(f.ax(6),xran);
    plot(f.ax(6),[trange(i,2)/3600 trange(i,2)/3600], yran, 'k--');    
    plot(f.ax(6),[trange(i,3)/3600 trange(i,3)/3600], yran, 'k--');    %,'linewidth',0.5
    xlabel(f.ax(6),strcat({'Hours on '},yr,{' '},day),'fontsize',12);
    hold(f.ax(6),'off');
    
    %%% save figure i
    print(f.fig,'-depsc',strcat(rstpath,'/',fam,'.contemp_migration.',num2str(trange(i,1)),...
          '_',num2str(trange(i,2)),'-',num2str(trange(i,3)),'.',num2str(winlenhf),'_',...
          num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.timerela.eps'));
    close(f.fig)  
    
end    
    
    
    
    
    
    