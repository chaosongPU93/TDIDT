% function mig_linear_fit_v2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to find the best propagation direction for migration, and then fit the hf
% the migration linearly to get the propagation slope and error distribution.
%
% distinguish from mig_linear_fit, mig_linear_fit_v2 changes the origin from 
% fam 002 in PGC trio to fam 043 in LZB trio to better comparison
%
% Chao Song, chaosong@princeton.edu
% First created date:   2020/09/24
% Last modified date:   2020/11/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
clc
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);

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

% relative coordinates to fam 043
loc043 = [-123.772167 48.493000 35.5900;];
[dx, dy] = absloc2relaloc(hfmaplocall(:,1),hfmaplocall(:,2),loc043(1),loc043(2));
hfrelalocall = [dx dy hfmaplocall];     % now increases to 12 cols

[dx, dy] = absloc2relaloc(lfmaplocall(:,1),lfmaplocall(:,2),loc043(1),loc043(2));
lfrelalocall = [dx dy lfmaplocall];     % now increases to 12 cols

% hfrelalocall = hfmaplocall;
[dx, dy] = absloc2relaloc(hfmaplocall(:,1),hfmaplocall(:,2),lon0,lat0);
% hfrelalocall(:,1) = dx;
% hfrelalocall(:,2) = dy;
hfrelalocall = [dx dy hfrelalocall];    % now increases to 14 cols

% lfrelalocall = lfmaplocall;
[dx, dy] = absloc2relaloc(lfmaplocall(:,1),lfmaplocall(:,2),lon0,lat0);
% lfrelalocall(:,1) = dx;
% lfrelalocall(:,2) = dy;
lfrelalocall = [dx dy lfrelalocall];    % now increases to 14 cols    

% so now the first 2 cols are  locations relative to fam 002, next 2 cols are relative to 043


%%
%%% 3. hf and lf migrations worthwhile to check the projection in propagation direction, maybe some
%%%     of them don't have many contemporaneous detections

% trange1 can be fitted by one slope
trange1 = [
          2005255,3.42e+4,3.65e+4;
          2005255,6.15e+4,6.26e+4;
          2005255,6.70e+4,6.90e+4;
          2005255,7.50e+4,7.60e+4;
          2005256,0.36e+4,0.55e+4;
          2005256,7.62e+4,8.30e+4;];     %2004197,8.45e+4,8.55e+4;
% angref = [290;
%          135;
%          135;
%          125;
%          135;
%          90;
%          225;
%          120;];
angbestl1 = zeros(length(trange1),1);
angbestl2 = zeros(length(trange1),1);
angbestl3 = zeros(length(trange1),1);
angbestl4 = zeros(length(trange1),1);
angbestl5 = zeros(length(trange1),1);

xran = [-10 30];
yran = [-30 10];

resprophf = nan(length(trange1)+1,300);
resproplf = nan(length(trange1)+1,100);

medallhf = nan(length(trange1),4);
medalllf = nan(length(trange1),4);
ranvechf98 = nan(length(trange1),2);

% array of vertical distance from LF from fitted HF line 
vertdistraw = nan(300 , length(trange1));    % raw vertical distance
vertdistwt = nan(300 , length(trange1));     % weighted vertical distance
weight = nan(300 , length(trange1));     % weights from regression
indlfindpen = nan(300 , length(trange1));    % store the index of indepedent LF detections

% create fit object
fttpfree = fittype( @(a,b,x) a*x+b);
        
for i = 1: length(trange1)
%     i=3;
    disp(trange1(i,:));
    indhf = find(hfrelalocall(:,12)==trange1(i,1) & hfrelalocall(:,14)>=trange1(i,2) & ...
                 hfrelalocall(:,14)<=trange1(i,3));
    mighf = hfrelalocall(indhf,:);
    
    angle = 0:5:360;
    
    l1norm = zeros(length(angle),1);
    l2norm = zeros(length(angle),1);
    slope = zeros(length(angle),1);
    sse = zeros(length(angle),1);
    rmse = zeros(length(angle),1);
    rsquare = zeros(length(angle),1);
    for iang = 1: length(angle)
        mighfdum = mighf;
        for j = 1: size(mighf,1)
            x0 = mighfdum(j,3);
            y0 = mighfdum(j,4);
            [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),0,0);
            mighfdum(j,3) = newx;
            mighfdum(j,4) = newy;
        end
        
        %%% linear robust least square
        [fitobj,gof,output] = fit(mighfdum(:,14)/3600, mighfdum(:,3),fttpfree,'Robust','Bisquare',...
                                  'StartPoint',[1 1]);
        % output fit parameters
        coef = coeffvalues(fitobj);
        slope(iang) = coef(1);
        fitprophf = feval(fitobj,mighfdum(:,14)/3600);
        l1norm(iang) = sum(abs(mighfdum(:,3)-fitprophf))/(length(mighfdum(:,3)));
        l2norm(iang) = sum((mighfdum(:,3)-fitprophf).^2)/(length(mighfdum(:,3)));
        sse(iang) = gof.sse;
        rmse(iang) = gof.rmse;
        rsquare(iang) = gof.rsquare;    % r-square, defined as ratio of the sum of squares of the regression and the total sum of squares
    
    end
    ind = find(slope>0);
    ind1 = find(l1norm(ind)==min(l1norm(ind)));
    angbestl1(i) = angle(ind(ind1(1)));
    
    ind2 = find(l2norm(ind)==min(l2norm(ind)));
    angbestl2(i) = angle(ind(ind2(1)));
    
    ind3 = find(rmse(ind)==min(rmse(ind)));     % rmse, Root Mean Squared Error, the smaller, the better
    angbestl3(i) = angle(ind(ind3(1)));
    
    ind4 = find(rsquare(ind)==max(rsquare(ind)));   % R-square, 0-1, the larger the better
    angbestl4(i) = angle(ind(ind4(1)));
    
    ind5 = find(sse(ind)==min(sse(ind)));   % sum of square due to error, the smaller, the better
    angbestl5(i) = angle(ind(ind5(1)));
      
    angbest(i) = median([angbestl3(i),angbestl4(i),angbestl5(i)]); 
    
    mighfdum = mighf;
    for j = 1: size(mighf,1)
        x0 = mighfdum(j,3);
        y0 = mighfdum(j,4);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        mighfdum(j,3) = newx;
        mighfdum(j,4) = newy;
    end
    
    indlf = find(lfrelalocall(:,12)==trange1(i,1) & lfrelalocall(:,14)>=trange1(i,2) & ...
                 lfrelalocall(:,14)<=trange1(i,3));
    miglf = lfrelalocall(indlf,:); 
    miglfdum = miglf;
    for j = 1: size(miglf,1)
        x0 = miglfdum(j,3);
        y0 = miglfdum(j,4);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        miglfdum(j,3) = newx;
        miglfdum(j,4) = newy;
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
    mighf = sortrows(mighf,-14);
    scatter(f.ax(1),mighf(:,3),mighf(:,4), 20, mighf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(1),'jet');
    c=colorbar(f.ax(1),'SouthOutside');
    pos = f.ax(1).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    juldate = num2str(trange1(i,1));
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
%     c.Label.String = strcat(num2str(trange1(i,1)),' of HF',' (hr)');
    c.Label.FontSize = 11;
    caxis(f.ax(1),[trange1(i,2)/3600 trange1(i,3)/3600])
    text(f.ax(1),0.85,0.1,'HF','FontSize',12,'unit','normalized');
    text(f.ax(1),0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(1),0.04,0.1,'PGC','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(1).Box = 'on';
    grid(f.ax(1), 'on');
    axis(f.ax(1), 'equal');
    f.ax(1).GridLineStyle = '--';
    f.ax(1).XAxisLocation = 'top';
    medxhf = median(mighf(:,3));
    medyhf = median(mighf(:,4));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medxhf-rotx medxhf+rotx];
    yvect = [medyhf-roty medyhf+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medxhf-rotx medxhf+rotx];
    yvect = [medyhf-roty medyhf+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',...
              [0.4 0.4 0.4]);
    xticks(f.ax(1),xran(1):5:xran(2));
    yticks(f.ax(1),yran(1):5:yran(2)); 
    xlabel(f.ax(1),'E (km)','fontsize',11);
    ylabel(f.ax(1),'N (km)','fontsize',11);
    if i == 1
        text(f.ax(1),0.5,0.1,num2str(i),'FontSize',12,'unit','normalized');
    else
        text(f.ax(1),0.5,0.1,num2str(i+1),'FontSize',12,'unit','normalized');
    end
    text(f.ax(1),0.83,0.95,strcat(num2str(size(mighf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(1),0.83,0.90,strcat({'in '},num2str(trange1(i,3)-trange1(i,2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(mighf,1)/(trange1(i,3)-trange1(i,2)));
    text(f.ax(1),0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    hold(f.ax(1),'off');

    % subplot 2 of figure i
    hold(f.ax(2),'on');
    plot(f.ax(2),[-100 100],[0 0],'k--');
    plot(f.ax(2),[0 0],[-100 100],'k--');
    f.ax(2).FontSize = 9;
    miglf = sortrows(miglf,-14);
    scatter(f.ax(2),miglf(:,3),miglf(:,4), 20, miglf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(2),'jet');
    c=colorbar(f.ax(2),'SouthOutside');
    pos = f.ax(2).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
%     c.Label.String = strcat(num2str(trange1(i,1)),' of LF',' (hr)');
    c.Label.FontSize = 11;
    caxis(f.ax(2),[trange1(i,2)/3600 trange1(i,3)/3600])
    text(f.ax(2),0.85,0.1,'LF','FontSize',12,'unit','normalized');
    text(f.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(2),0.04,0.1,'PGC','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(2).Box = 'on';
    grid(f.ax(2), 'on');
    axis(f.ax(2), 'equal');
    f.ax(2).GridLineStyle = '--';
    f.ax(2).XAxisLocation = 'top';
    medxlf = median(miglf(:,3));
    medylf = median(miglf(:,4));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',...
              [0.4 0.4 0.4]);
    xticks(f.ax(2),xran(1):5:xran(2));
    yticks(f.ax(2),yran(1):5:yran(2));    
    xlabel(f.ax(2),'E (km)','fontsize',11);
    ylabel(f.ax(2),'N (km)','fontsize',11);
    if i == 1
        text(f.ax(2),0.5,0.1,num2str(i),'FontSize',12,'unit','normalized');
    else
        text(f.ax(2),0.5,0.1,num2str(i+1),'FontSize',12,'unit','normalized');
    end
    text(f.ax(2),0.83,0.95,strcat(num2str(size(miglf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(2),0.83,0.9,strcat({'in '},num2str(trange1(i,3)-trange1(i,2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(miglf,1)/(trange1(i,3)-trange1(i,2)));
    text(f.ax(2),0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    hold(f.ax(2),'off');
    
    % subplot 3 of figure i
    hold(f.ax(3),'on');
    f.ax(3).FontSize = 9;
    scatter(f.ax(3),mighfdum(:,14)/3600,mighfdum(:,3),20,[0.6 0.6 0.6],'filled','o',...
            'MarkerEdgeColor','k');      
    scatter(f.ax(3),miglfdum(:,14)/3600,miglfdum(:,3),20,[0.6 1 1],'filled','o',...
            'MarkerEdgeColor','k');
%     scatter(f.ax(3),mighfdum(:,14)/3600,mighfdum(:,3),30,[105/255 105/255 105/255],'filled','o');
%     scatter(f.ax(3),miglfdum(:,14)/3600,miglfdum(:,3),30,miglfdum(:,11),'filled','o');  

    [fitobjhfprop,gofhfprop,outphfprop] = fit(mighfdum(:,14)/3600, mighfdum(:,3),fttpfree,'Robust',...
                                              'Bisquare','StartPoint',[1 1]);
    % output fit parameters
    coefprop = coeffvalues(fitobjhfprop);
    slopeprophf(i) = coefprop(1);
    intcptprophf(i) = coefprop(2);
    fitprophf = feval(fitobjhfprop,mighfdum(:,14)/3600);
    plot(f.ax(3),mighfdum(:,14)/3600,fitprophf,'-','linewidth',2,'color',[0.6 0.6 0.6]);
        
%     intcptproplf = ones(size(miglfdum(:,end)))\(miglfdum(:,3)-slopeprophf*miglfdum(:,end)/3600); % direct LS fit
%     fitproplf = slopeprophf*miglfdum(:,end)/3600+intcptproplf;
    a=slopeprophf(i);
    fttpfix = fittype( @(b,x) a*x+b);
    [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum(:,14)/3600, miglfdum(:,3),fttpfix,'Robust',...
                                              'Bisquare','StartPoint',1);
%     [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum(:,14)/3600, miglfdum(:,3),fttpfix,'Robust',...
%                                               'Bisquare','Weights',miglfdum(:,11));
    fitproplf = feval(fitobjlfprop,miglfdum(:,14)/3600);
    intcptproplf(i) = coeffvalues(fitobjlfprop);
    offset(i) = intcptproplf(i)-intcptprophf(i);
    plot(f.ax(3),miglfdum(:,14)/3600,fitproplf,'-','linewidth',2,'color',[0.6 1 1]);
    f.ax(3).Box = 'on';
    grid(f.ax(3), 'on');
    f.ax(3).GridLineStyle = '--';
    xran1 = [trange1(i,2)/3600 trange1(i,3)/3600];
    yran1 = [round(min([mighfdum(:,3);miglfdum(:,3)]))-1 ...
             round(max([mighfdum(:,3);miglfdum(:,3)]))+1];
    xlim(f.ax(3),xran1);
    ylim(f.ax(3),yran1);
    text(f.ax(3),0.04,0.91,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    xlabel(f.ax(3),'Time (hr)','fontsize',11);
    ylabel(f.ax(3),'Dist. along prop. (km)','fontsize',11); 
    hold(f.ax(3),'off');
    
    % compute the HF weights in robust linear regression, see NOTES above
    rhf = outphfprop.residuals;   % usual residuals
    x = [mighfdum(:,14)/3600];
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
    x = [miglfdum(:,14)/3600];
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
    fittemp = feval(fitobjhfprop,miglfdum(:,14)/3600);
    verttemp = miglfdum(:,3) - fittemp;     % can be + or -
    vertdistraw(1: length(verttemp), i) = verttemp;     % raw vertical distance
    verttempwt = wtlf.*verttemp;      % weighted vertical distance by multiplying the weight
    vertdistwt(1: length(verttempwt), i) = verttempwt;
    weight(1: length(verttempwt), i) =  wtlf;
        
    % get the information of propagation vector indicating range etc.
    [tmpx,tmpy] = coordinate_rot(medxhf,medyhf,-(angbest(i)-90),0,0);
    medallhf(i,1:4) = [medxhf medyhf tmpx tmpy];
    [tmpx,tmpy] = coordinate_rot(medxlf,medylf,-(angbest(i)-90),0,0);
    medalllf(i,1:4) = [medxlf medylf tmpx tmpy];
    
    ranvechf98(i,1:2) = [prctile(mighfdum(:,3),2) prctile(mighfdum(:,3),98)];
        
    nlt = sum(mighfdum(:,3)< ranvechf98(i,1));      % number that lower than the current range
    ngt = sum(mighfdum(:,3)> ranvechf98(i,2));      % number that greater than the current range
    
    mighfsort = sort(mighfdum(:,3));    % sort in descend order
    if nlt < 2
        ranvechf98(i,1) = mighfsort(3);     % so at least throw away 2 abnormals
    end
    
    if ngt < 2
        ranvechf98(i,2) = mighfsort(end-2);     % so at least throw away 2 abnormals
    end
    
%     %%% find the independent LF detections, i.e. non-overlapping time window
%     [miglfindpen,indindpen] = find_independent_detection(miglfdum, 0);
%     inumlf(i) = length(indindpen);
%     indlfindpen(1: inumlf(i), i) =  indindpen;
    
    % subplot 4 of figure i
    hold(f.ax(4),'on');
    f.ax(4).FontSize = 9;
    numhf(i) = length(mighfdum(:,3));
    numlf(i) = length(miglfdum(:,3));
    resprophf(i,1:numhf(i)) = outphfprop.residuals;
    resproplf(i,1:numlf(i)) = outplfprop.residuals+(intcptproplf(i)-intcptprophf(i));
    [xlochf, ~, ~, pdfvalhf] = weighted_bar_pdf(resprophf(i,1:numhf(i)), wthf, 0.5, 'int');
    hfHdl = bar(f.ax(4),xlochf, pdfvalhf,1,'stacked');
    hfHdl(1).FaceColor = [0.6 0.6 0.6];
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
%     plot(min(resprophf(i,1:numhf(i))):0.05:max(resprophf(i,1:numhf(i))), pdffithf, 'r-');
%     plot(min(resproplf(i,1:numlf(i))):0.05:max(resproplf(i,1:numlf(i))), pdffitlf, '-','color',...
%          [0.6 1 1]);
    text(f.ax(4),0.04,0.91,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(4),0.4,0.92,sprintf('LF-HF = %.2f km',offset(i)),'FontSize',11,'unit','normalized');
    f.ax(4).Box = 'on';
    grid(f.ax(4), 'on');
    f.ax(4).GridLineStyle = '--';
    ymax = f.ax(4).YLim(2)+0.1;
    ylim(f.ax(4),[0 ymax]);
    xlabel(f.ax(4),'Residual in prop. (km)','fontsize',11);
    ylabel(f.ax(4),'PDF estimate','fontsize',11);
    hold(f.ax(4),'off');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% It seems that median of the previous detections is a promising way to show offsets as well,
    %%% so try to deal with it the same way
    medhf = [];
    for j = 1: size(mighfdum,1)
        medhf(j,1) = median(mighfdum(1:j,3));
        medhf(j,2) = median(mighfdum(1:j,4));
        medhf(j,3) = mighfdum(j,14);
    end
    
    medlf = [];
    for j = 1: size(miglfdum,1)
        medlf(j,1) = median(miglfdum(1:j,3));
        medlf(j,2) = median(miglfdum(1:j,4));
        medlf(j,3) = miglfdum(j,14);
    end
       
    
    % subplot 5 of figure i
    hold(f.ax(5),'on');
    f.ax(5).FontSize = 9;
    scatter(f.ax(5),medhf(:,3)/3600,medhf(:,1),20,[0.6 0.6 0.6],'filled','o','MarkerEdgeColor',...
            'k');
    scatter(f.ax(5),medlf(:,3)/3600,medlf(:,1),20,[0.6 1 1],'filled','o','MarkerEdgeColor',...
            'k');            
    f.ax(5).Box = 'on';
    grid(f.ax(5), 'on');
    f.ax(5).GridLineStyle = '--';
    xran1 = [trange1(i,2)/3600 trange1(i,3)/3600];
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
    scatter(f.ax(6),medhf(:,3)/3600,medhf(:,2),20,[0.6 0.6 0.6],'filled','o','MarkerEdgeColor',...
            'k');
    scatter(f.ax(6),medlf(:,3)/3600,medlf(:,2),20,[0.6 1 1],'filled','o','MarkerEdgeColor',...
            'k');          
    f.ax(6).Box = 'on';
    grid(f.ax(6), 'on');
    f.ax(6).GridLineStyle = '--';
    xran1 = [trange1(i,2)/3600 trange1(i,3)/3600];
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
    print(f.fig,'-dpdf',strcat(rstpath,'/',fam,'mig.proj.lfit.meddist',num2str(trange1(i,1)),...
          '_',num2str(trange1(i,2)),'-',num2str(trange1(i,3)),'.',num2str(winlenhf),'_',...
          num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));

end   


%%
% trange2 needs to 2 slopes, visually seperate the time      
trange2 = [2005255,5.80e+4,5.96e+4];
divtime = 58830;
xran = [-8 13];
yran = [-13 10];

xran = [-10 30];
yran = [-30 10];

indhf = find(hfrelalocall(:,12)==trange2(1,1) & hfrelalocall(:,14)>=trange2(1,2) & ...
             hfrelalocall(:,14)<=trange2(1,3));
mighf = hfrelalocall(indhf,:);

indlf = find(lfrelalocall(:,12)==trange2(1,1) & lfrelalocall(:,14)>=trange2(1,2) & ...
                 lfrelalocall(:,14)<=trange2(1,3));
miglf = lfrelalocall(indlf,:);

% miglfp1(4,:)=[];

angle = 0:5:360;

slope = zeros(length(angle),1);
sse = zeros(length(angle),1);
rmse = zeros(length(angle),1);
rsquare = zeros(length(angle),1);

%%% linear robust least square
fttpfree = fittype( @(a,b,x) a*x+b);
    
for iang = 1: length(angle)
    
    mighfdum = mighf;
    for j = 1: size(mighf,1)
        x0 = mighfdum(j,3);
        y0 = mighfdum(j,4);
        [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),0,0);
        mighfdum(j,3) = newx;
        mighfdum(j,4) = newy;
    end    

    % create fit object
    [fitobj,gof,output] = fit(mighfdum(:,14)/3600, mighfdum(:,3),fttpfree,'Robust','Bisquare',...
                              'StartPoint',[1 1]);
    % output fit parameters
    coef = coeffvalues(fitobj);
    slope(iang) = coef(1);
    fithf = feval(fitobj,mighfdum(:,14)/3600);
    sse(iang) = gof.sse;
    rmse(iang) = gof.rmse;
    rsquare(iang) = gof.rsquare;    % r-square, defined as ratio of the sum of squares of the regression and the total sum of squares
    
    
end

ind3 = find(rmse(ind)==min(rmse(ind)));     % rmse, Root Mean Squared Error, the smaller, the better
angbestl3 = angle(ind(ind3(1)));

ind4 = find(rsquare(ind)==max(rsquare(ind)));   % R-square, 0-1, the larger the better
angbestl4 = angle(ind(ind4(1)));

ind5 = find(sse(ind)==min(sse(ind)));   % sum of square due to error, the smaller, the better
angbestl5 = angle(ind(ind5(1)));

angbest(7) = median([angbestl3,angbestl4,angbestl5]);

% propagation on the propagation direction
mighfdum = mighf;
for j = 1: size(mighf,1)
    x0 = mighfdum(j,3);
    y0 = mighfdum(j,4);
    [newx,newy] = coordinate_rot(x0,y0,-(angbest(7)-90),0,0);
    mighfdum(j,3) = newx;
    mighfdum(j,4) = newy;
end
miglfdum = miglf;
for j = 1: size(miglf,1)
    x0 = miglfdum(j,3);
    y0 = miglfdum(j,4);
    [newx,newy] = coordinate_rot(x0,y0,-(angbest(7)-90),0,0);
    miglfdum(j,3) = newx;
    miglfdum(j,4) = newy;
end

% break to 2 segments
mighfp1dum = mighfdum(mighfdum(:,14)<divtime,:);
mighfp2dum = mighfdum(mighfdum(:,14)>=divtime,:);
miglfp1dum = miglfdum(miglfdum(:,14)<divtime,:);
miglfp2dum = miglfdum(miglfdum(:,14)>=divtime,:);


%%% define and position the figure frame and axes of each plot
f3.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f3.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 3;
ncol = 2;
for isub = 1:nrow*ncol
    f3.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f3.ax(1), 'position', [ 0.1, 0.64, 0.32, 0.32]);
set(f3.ax(2), 'position', [ 0.5, 0.64, 0.32, 0.32]);
set(f3.ax(3), 'position', [ 0.1, 0.35, 0.32, 0.22]);
set(f3.ax(4), 'position', [ 0.5, 0.35, 0.32, 0.22]);
set(f3.ax(5), 'position', [ 0.1, 0.078, 0.32, 0.22]);
set(f3.ax(6), 'position', [ 0.5, 0.078, 0.32, 0.22]);

% subplot 1 of figure i
hold(f3.ax(1),'on');
plot(f3.ax(1),[-100 100],[0 0],'k--');
plot(f3.ax(1),[0 0],[-100 100],'k--');
f3.ax(1).FontSize = 9;
mighf = sortrows(mighf,-14);
scatter(f3.ax(1),mighf(:,3),mighf(:,4), 20, mighf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
colormap(f3.ax(1),'jet');
c=colorbar(f3.ax(1),'SouthOutside');
pos = f3.ax(1).Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
% c.TickLabels=[];
juldate = num2str(trange2(1,1));
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
% c.Label.String = strcat(num2str(trange2(1,1)),' of HF',' (hr)');
c.Label.FontSize = 11;
caxis(f3.ax(1),[trange2(1,2)/3600 trange2(1,3)/3600])
text(f3.ax(1),0.85,0.1,'HF','FontSize',12,'unit','normalized');
text(f3.ax(1),0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(f3.ax(1),0.04,0.1,'PGC','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
f3.ax(1).Box = 'on';
grid(f3.ax(1), 'on');
axis(f3.ax(1), 'equal');
f3.ax(1).GridLineStyle = '--';
f3.ax(1).XAxisLocation = 'top';
medxhf = median(mighf(:,3));
medyhf = median(mighf(:,4));
[rotx, roty] = complex_rot(0,5,-angbest(7));
xvect = [medxhf-rotx medxhf+rotx];
yvect = [medyhf-roty medyhf+roty];
drawArrow(f3.ax(1),xvect,yvect,xran,yran,'linewidth',1.5);
[rotx, roty] = complex_rot(-5,0,-angbest(7));
xvect = [medxhf-rotx medxhf+rotx];
yvect = [medyhf-roty medyhf+roty];
drawArrow(f3.ax(1),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
xticks(f3.ax(1),xran(1):5:xran(2));
yticks(f3.ax(1),yran(1):5:yran(2));
xlabel(f3.ax(1),'E (km)','fontsize',11);
ylabel(f3.ax(1),'N (km)','fontsize',11);
text(f3.ax(1),0.5,0.1,num2str(2),'FontSize',12,'unit','normalized');
text(f3.ax(1),0.83,0.95,strcat(num2str(size(mighf,1)),{' detections'}),'FontSize',8,...
    'unit','normalized','horizontalalignment','center');
text(f3.ax(1),0.83,0.90,strcat({'in '},num2str(trange2(1,3)-trange2(1,2)),{' s'}),'FontSize',8,...
    'unit','normalized','horizontalalignment','center');
rate = sprintf('%.3f',size(mighf,1)/(trange2(1,3)-trange2(1,2)));
text(f3.ax(1),0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
    'unit','normalized','horizontalalignment','center');
hold(f3.ax(1),'off');

% subplot 2 of figure i
hold(f3.ax(2),'on');
plot(f3.ax(2),[-100 100],[0 0],'k--');
plot(f3.ax(2),[0 0],[-100 100],'k--');
f3.ax(2).FontSize = 9;
miglf = sortrows(miglf,-14);
scatter(f3.ax(2),miglf(:,3),miglf(:,4), 20, miglf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
colormap(f3.ax(2),'jet');
c=colorbar(f3.ax(2),'SouthOutside');
pos = f3.ax(2).Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
% c.TickLabels=[];
c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
% c.Label.String = strcat(num2str(trange2(1,1)),' of LF',' (hr)');
c.Label.FontSize = 11;
caxis(f3.ax(2),[trange2(1,2)/3600 trange2(1,3)/3600])
text(f3.ax(2),0.85,0.1,'LF','FontSize',12,'unit','normalized');
text(f3.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(f3.ax(2),0.04,0.1,'PGC','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
f3.ax(2).Box = 'on';
grid(f3.ax(2), 'on');
axis(f3.ax(2), 'equal');
f3.ax(2).GridLineStyle = '--';
f3.ax(2).XAxisLocation = 'top';
medxlf = median(miglf(:,3));
medylf = median(miglf(:,4));
[rotx, roty] = complex_rot(0,5,-angbest(7));
xvect = [medxlf-rotx medxlf+rotx];
yvect = [medylf-roty medylf+roty];
drawArrow(f3.ax(2),xvect,yvect,xran,yran,'linewidth',1.5);
[rotx, roty] = complex_rot(-5,0,-angbest(7));
xvect = [medxlf-rotx medxlf+rotx];
yvect = [medylf-roty medylf+roty];
drawArrow(f3.ax(2),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
xticks(f3.ax(2),xran(1):5:xran(2));
yticks(f3.ax(2),yran(1):5:yran(2));
xlabel(f3.ax(2),'E (km)','fontsize',11);
ylabel(f3.ax(2),'N (km)','fontsize',11);
text(f3.ax(2),0.5,0.1,num2str(2),'FontSize',12,'unit','normalized');
text(f3.ax(2),0.83,0.95,strcat(num2str(size(miglf,1)),{' detections'}),'FontSize',8,...
    'unit','normalized','horizontalalignment','center');
text(f3.ax(2),0.83,0.9,strcat({'in '},num2str(trange2(1,3)-trange2(1,2)),{' s'}),'FontSize',8,...
    'unit','normalized','horizontalalignment','center');
rate = sprintf('%.3f',size(miglf,1)/(trange2(1,3)-trange2(1,2)));
text(f3.ax(2),0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
    'unit','normalized','horizontalalignment','center');
hold(f3.ax(2),'off');

% subplot 3 of figure i
hold(f3.ax(3),'on');
f3.ax(3).FontSize = 9;
scatter(f3.ax(3),mighfp1dum(:,14)/3600,mighfp1dum(:,3),20,[0.6 0.6 0.6],'filled','o',...
            'MarkerEdgeColor','k');   
scatter(f3.ax(3),mighfp2dum(:,14)/3600,mighfp2dum(:,3),20,[0.6 0.6 0.6],'filled','o',...
            'MarkerEdgeColor','k');   
scatter(f3.ax(3),miglfp1dum(:,14)/3600,miglfp1dum(:,3),20,[0.6 1 1],'filled','o',...
            'MarkerEdgeColor','k');
scatter(f3.ax(3),miglfp2dum(:,14)/3600,miglfp2dum(:,3),20,[0.6 1 1],'filled','o',...
            'MarkerEdgeColor','k');
%%% segment 1
% create fit object
[fitobjhfp1,gofhfp1,outphfp1] = fit(mighfp1dum(:,14)/3600, mighfp1dum(:,3),fttpfree,'Robust',...
                                    'Bisquare','StartPoint',[1 1]);
% output fit parameters
coefp1 = coeffvalues(fitobjhfp1);
slopep1hf = coefp1(1);
intcptp1hf = coefp1(2);
fitp1hf = feval(fitobjhfp1,mighfp1dum(:,14)/3600);
plot(f3.ax(3),mighfp1dum(:,14)/3600,fitp1hf,'-','linewidth',2,'color',[0.6 0.6 0.6]);

% intcptp1lf = ones(size(miglfp1dum(:,end)))\(miglfp1dum(:,1)-slopep1hf*miglfp1dum(:,end)/3600);
a=slopep1hf;
fttpfix = fittype( @(b,x) a*x+b);
[fitobjlfp1,goflfp1,outplfp1] = fit(miglfp1dum(:,14)/3600, miglfp1dum(:,3),fttpfix,'Robust',...
                                    'Bisquare','StartPoint',1);
% [fitobjlfp1,goflfp1,outplfp1] = fit(miglfp1dum(:,14)/3600, miglfp1dum(:,1),fttpfix,'Robust',...
%                                     'Bisquare','Weights',miglfp1dum(:,11));
intcptp1lf = coeffvalues(fitobjlfp1);
fitp1lf = feval(fitobjlfp1,miglfp1dum(:,14)/3600);
plot(f3.ax(3),miglfp1dum(:,14)/3600,fitp1lf,'-','linewidth',2,'color',[0.6 1 1]);

%%% segment 2
% create fit object
[fitobjhfp2,gofhfp2,outphfp2] = fit(mighfp2dum(:,14)/3600, mighfp2dum(:,3),fttpfree,'Robust',...
                                    'Bisquare','StartPoint',[1 1]);
% output fit parameters
coefp2 = coeffvalues(fitobjhfp2);
slopep2hf = coefp2(1);
intcptp2hf = coefp2(2);
fitp2hf = feval(fitobjhfp2,mighfp2dum(:,14)/3600);
plot(f3.ax(3),mighfp2dum(:,14)/3600,fitp2hf,'-','linewidth',2,'color',[0.6 0.6 0.6]);
% intcptp1lf = ones(size(miglfp1dum(:,end)))\(miglfp1dum(:,1)-slopep1hf*miglfp1dum(:,end)/3600);
a=slopep2hf;
fttpfix = fittype( @(b,x) a*x+b);
[fitobjlfp2,goflfp2,outplfp2] = fit(miglfp2dum(:,14)/3600, miglfp2dum(:,3),fttpfix,'Robust',...
                                    'Bisquare','StartPoint',1);
% [fitobjlfp2,goflfp2,outplfp2] = fit(miglfp2dum(:,14)/3600, miglfp2dum(:,1),fttpfix,'Robust',...
%                                     'Bisquare','Weights',miglfp2dum(:,11));
intcptp2lf = coeffvalues(fitobjlfp2);
fitp2lf = feval(fitobjlfp2,miglfp2dum(:,14)/3600);
offset = (intcptp1lf-intcptp1hf+intcptp2lf-intcptp2hf)/2;
plot(f3.ax(3),miglfp2dum(:,14)/3600,fitp2lf,'-','linewidth',2,'color',[0.6 1 1]);

f3.ax(3).Box = 'on';
grid(f3.ax(3), 'on');
f3.ax(3).GridLineStyle = '--';
xran1 = [trange2(1,2)/3600 trange2(1,3)/3600];
yran1 = [round(min([mighfp1dum(:,3);mighfp2dum(:,3);miglfp1dum(:,3);miglfp2dum(:,3)]))-1 ...
         round(max([mighfp1dum(:,3);mighfp2dum(:,3);miglfp1dum(:,3);miglfp2dum(:,3)]))+1];
xlim(f3.ax(3),xran1);
ylim(f3.ax(3),yran1);
text(f3.ax(3),0.04,0.91,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
xlabel(f3.ax(3),'Time (hr)','fontsize',11);
ylabel(f3.ax(3),'Dist. along prop. (km)','fontsize',11);
hold(f3.ax(3),'off');

% compute the HF weights in robust linear regression, see NOTES above
rhf = outphfp1.residuals;   % usual residuals
x = [mighfp1dum(:,14)/3600];
hatmat = x*inv(x'*x)*x';
h = zeros(size(hatmat,1),1);    % leverage of least square
for jj = 1 : size(hatmat,1)
    h(jj) = hatmat(jj,jj);
end
radj = rhf./sqrt(1-h);      % adjusted residuals
K = 4.685;
s = mad(rhf,1)/0.6745;
u = radj/(K*s);
wthfp1 = zeros(length(u),1);    % rubust weight of next iteration
for jj = 1 : length(u)
    if abs(u(jj)) < 1
        wthfp1(jj) = (1-(u(jj))^2)^2;
    else
        wthfp1(jj) = 0;
    end
end
    
rhf = outphfp2.residuals;   % usual residuals
x = [mighfp2dum(:,14)/3600];
hatmat = x*inv(x'*x)*x';
h = zeros(size(hatmat,1),1);    % leverage of least square
for jj = 1 : size(hatmat,1)
    h(jj) = hatmat(jj,jj);
end
radj = rhf./sqrt(1-h);      % adjusted residuals
K = 4.685;
s = mad(rhf,1)/0.6745;
u = radj/(K*s);
wthfp2 = zeros(length(u),1);    % rubust weight of next iteration
for jj = 1 : length(u)
    if abs(u(jj)) < 1
        wthfp2(jj) = (1-(u(jj))^2)^2;
    else
        wthfp2(jj) = 0;
    end
end
wthf = [wthfp1; wthfp2];

% compute the LF weights in robust linear regression, see NOTES above
rlf = outplfp1.residuals;   % usual residuals
x = [miglfp1dum(:,14)/3600];
hatmat = x*inv(x'*x)*x';
h = zeros(size(hatmat,1),1);    % leverage of least square
for jj = 1 : size(hatmat,1)
    h(jj) = hatmat(jj,jj);
end
radj = rlf./sqrt(1-h);      % adjusted residuals
K = 4.685;
s = mad(rlf,1)/0.6745;
u = radj/(K*s);
wtlfp1 = zeros(length(u),1);    % rubust weight of next iteration
for jj = 1 : length(u)
    if abs(u(jj)) < 1
        wtlfp1(jj) = (1-(u(jj))^2)^2;
    else
        wtlfp1(jj) = 0;
    end
end
    
rlf = outplfp2.residuals;   % usual residuals
x = [miglfp2dum(:,14)/3600];
hatmat = x*inv(x'*x)*x';
h = zeros(size(hatmat,1),1);    % leverage of least square
for jj = 1 : size(hatmat,1)
    h(jj) = hatmat(jj,jj);
end
radj = rlf./sqrt(1-h);      % adjusted residuals
K = 4.685;
s = mad(rlf,1)/0.6745;
u = radj/(K*s);
wtlfp2 = zeros(length(u),1);    % rubust weight of next iteration
for jj = 1 : length(u)
    if abs(u(jj)) < 1
        wtlfp2(jj) = (1-(u(jj))^2)^2;
    else
        wtlfp2(jj) = 0;
    end
end
wtlf = [wtlfp1; wtlfp2];

% get the vertical distance from LF detections to the fitted HF line
fittemp1 = feval(fitobjhfp1,miglfp1dum(:,14)/3600);
verttemp1 = miglfp1dum(:,3) - fittemp1;     % can be + or -
fittemp2 = feval(fitobjhfp2,miglfp2dum(:,14)/3600);
verttemp2 = miglfp2dum(:,3) - fittemp2;     % can be + or -
verttemp = [verttemp1; verttemp2];
vertdistraw(1: length(verttemp), 7) = verttemp;     % raw vertical distance
verttempwt = wtlf.*verttemp;      % weighted vertical distance by multiplying the weight
vertdistwt(1: length(verttempwt), 7) = verttempwt;
weight(1: length(verttempwt), 7) =  wtlf;

% %%% find the independent LF detections, i.e. non-overlapping time window
% [miglfindpen,indindpen] = find_independent_detection(miglfdum, 0);
% inumlf(i) = length(indindpen);
% indlfindpen(1: inumlf(i), i) =  indindpen;
    
% subplot 4 of figure i
hold(f3.ax(4),'on');
f3.ax(4).FontSize = 9;
numhf(7) = length(mighfdum(:,3));
numlf(7) = length(miglfdum(:,3));
resprophf(7,1:numhf(7)) = [outphfp1.residuals; outphfp2.residuals];
resproplf(7,1:numlf(7)) = [outplfp1.residuals+(intcptp1lf-intcptp1hf); outplfp2.residuals+(intcptp2lf-intcptp2hf)];
[xlochf, ~, ~, pdfvalhf] = weighted_bar_pdf(resprophf(7,1:numhf(7)), wthf, 0.5, 'int');
hfHdl = bar(f3.ax(4),xlochf, pdfvalhf,1,'stacked');
hfHdl(1).FaceColor = [0.6 0.6 0.6];
[xloclf, ~, ~, pdfvallf] = weighted_bar_pdf(resproplf(7,1:numlf(7)), wtlf, 0.5, 'int');
lfHdl = bar(f3.ax(4),xloclf, pdfvallf,1,'stacked','facea',0.6);
lfHdl(1).FaceColor = [0.6 1 1];
% histogram(f3.ax(4),resprophf(7,1:numhf(7)),'binwidth',0.5,'normalization','pdf','facecolor',...
%           [1 0 0]);
% histogram(f3.ax(4),resproplf(7,1:numlf(7)),'binwidth',0.5,'normalization','pdf','facecolor',...
%           [0.6 1 1],'facealpha',0.6);
% pdhf = fitdist(resprophf(7,1:numhf(7))','Normal');    % fit a distribution of residual, assuming it is normal distributed
% pdlf = fitdist(resproplf(7,1:numlf(7))','Normal');
% muhf = pdhf.mu;    % fitted paramaters
% mulf = pdlf.mu;
% sigmahf = pdhf.sigma;
% sigmalf = pdlf.sigma;
% cihf = paramci(pdhf,'Alpha',0.05);     % give the 95% confidence interval of distribution param
% cilf = paramci(pdlf,'Alpha',0.05);
% pdffithf = pdf(pdhf,min(resprophf(7,1:numhf(7))):0.05:max(resprophf(7,1:numhf(7))));
% pdffitlf = pdf(pdlf,min(resproplf(7,1:numlf(7))):0.05:max(resproplf(7,1:numlf(7))));
% plot(min(resprophf(7,1:numhf(7))):0.05:max(resprophf(7,1:numhf(7))), pdffithf, 'r-');
% plot(min(resproplf(7,1:numlf(7))):0.05:max(resproplf(7,1:numlf(7))), pdffitlf, '-','color',...
%          [153/255 255/255 255/255]);
text(f3.ax(4),0.04,0.91,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(f3.ax(4),0.40,0.92,sprintf('LF-HF = %.2f km',offset),'FontSize',11,'unit','normalized');
f3.ax(4).Box = 'on';
grid(f3.ax(4), 'on');
f3.ax(4).GridLineStyle = '--';
ymax = f3.ax(4).YLim(2)+0.1;
ylim(f3.ax(4),[0 ymax]);
xlabel(f3.ax(4),'Residual in prop. (km)','fontsize',11);
ylabel(f3.ax(4),'PDF estimate','fontsize',11);
hold(f3.ax(4),'off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% It seems that median of the previous detections is a promising way to show offsets as well,
%%% so try to deal with it the same way
medhf = [];
for j = 1: size(mighfdum,1)
    medhf(j,1) = median(mighfdum(1:j,3));
    medhf(j,2) = median(mighfdum(1:j,4));
    medhf(j,3) = mighfdum(j,14);
end

medlf = [];
for j = 1: size(miglfdum,1)
    medlf(j,1) = median(miglfdum(1:j,3));
    medlf(j,2) = median(miglfdum(1:j,4));
    medlf(j,3) = miglfdum(j,14);
end


% subplot 5 of figure i
hold(f3.ax(5),'on');
f3.ax(5).FontSize = 9;
scatter(f3.ax(5),medhf(:,3)/3600,medhf(:,1),20,[0.6 0.6 0.6],'filled','o','MarkerEdgeColor',...
        'k');
scatter(f3.ax(5),medlf(:,3)/3600,medlf(:,1),20,[0.6 1 1],'filled','o','MarkerEdgeColor',...
        'k');            
f3.ax(5).Box = 'on';
grid(f3.ax(5), 'on');
f3.ax(5).GridLineStyle = '--';
xran1 = [trange2(1,2)/3600 trange2(1,3)/3600];
yran1 = [round(min([medhf(:,1);medlf(:,1)]))-1 ...
         round(max([medhf(:,1);medlf(:,1)]))+1];
xlim(f3.ax(5),xran1);
ylim(f3.ax(5),yran1);
text(f3.ax(5),0.97,0.1,'e','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2,...
         'HorizontalAlignment','right');
xlabel(f3.ax(5),'Time (hr)','fontsize',11);
ylabel(f3.ax(5),'Med. dist. along prop. (km)','fontsize',11);
hold(f3.ax(5),'off');

% subplot 4 of figure i
hold(f3.ax(6),'on');
f3.ax(6).FontSize = 9;
scatter(f3.ax(6),medhf(:,3)/3600,medhf(:,2),20,[0.6 0.6 0.6],'filled','o','MarkerEdgeColor',...
        'k');
scatter(f3.ax(6),medlf(:,3)/3600,medlf(:,2),20,[0.6 1 1],'filled','o','MarkerEdgeColor',...
        'k');          
f3.ax(6).Box = 'on';
grid(f3.ax(6), 'on');
f3.ax(6).GridLineStyle = '--';
xran1 = [trange2(1,2)/3600 trange2(1,3)/3600];
yran1 = [round(min([medhf(:,2);medlf(:,2)]))-1 ...
         round(max([medhf(:,2);medlf(:,2)]))+1];
xlim(f3.ax(6),xran1);
ylim(f3.ax(6),yran1);
text(f3.ax(6),0.97,0.1,'f','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2,...
         'HorizontalAlignment','right');
xlabel(f3.ax(6),'Time (hr)','fontsize',11);
ylabel(f3.ax(6),'Med. dist. along ort. (km)','fontsize',11);
hold(f3.ax(6),'off');

%%% save figure
print(f3.fig,'-dpdf',strcat(rstpath,'/',fam,'mig.proj.lfit.meddist',num2str(trange2(1,1)),...
    '_',num2str(trange2(1,2)),'-',num2str(trange2(1,3)),'.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));

    
%% plot the vertical distance from lf points to hf line in certain regions

% NOTE here, 2020/09/22, u can choose which vertical distance to use
vertdist = vertdistwt;
% vertdist = vertdistraw;

% SW+NE migs
ind002 = [3,5];    % index of SW migrations
vdist002 = [];
for i = 1: length(ind002)
    vdist002 = [vdist002; vertdist(1:numlf(ind002(i)), ind002(i))];
end
pd002 = fitdist(vdist002,'Normal');  
mu002 = pd002.mu;    % fitted paramaters
sig002 = pd002.sigma;
pdffit002 = pdf(pd002,-50:0.05:50);

%%% plot 3, NW regions + region of 002, in single plot
f11.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 4.5;   % maximum height allowed is 11 inches
set(f11.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f11.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f11.ax(1), 'position', [ 0.1, 0.1, 0.35, 0.75]);
% set(f.ax(2), 'position', [ 0.55, 0.1, 0.35, 0.75]);

% subplot 1
ax = f11.ax(1);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
xran = [-15 20];
yran = [-25 5];
indplt = ind002;    % index to plot
for i = 1: length(indplt)
    [rotx1, roty1] = complex_rot(0, ranvechf98(indplt(i),2)-medallhf(indplt(i),3), ...
                                 -angbest(indplt(i)));
    [rotx2, roty2] = complex_rot(0, medallhf(indplt(i),3)-ranvechf98(indplt(i),1), ...
                                 -angbest(indplt(i)));
    xvect = [medallhf(indplt(i),1)-rotx2 medallhf(indplt(i),1)+rotx1];
    yvect = [medallhf(indplt(i),2)-roty2 medallhf(indplt(i),2)+roty1];
    drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1,'linestyle','-','color',[0.7 0.7 0.7]);
   
%     text(ax,0.5,0.1,num2str(indplt(i)),'FontSize',12,'unit','normalized');
end

for i = 1: length(indplt)
    if indplt(i) == 1
        scatter(ax,medallhf(indplt(i),1),medallhf(indplt(i),2), 30, indplt(i), 'filled','o');  %, 'MarkerEdgeColor', 'w')
    else
        scatter(ax,medallhf(indplt(i),1),medallhf(indplt(i),2), 30, indplt(i)+1, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    end
end
colormap(ax,'jet');
c=colorbar(ax,'SouthOutside');
caxis(ax,[1 size(trange1,1)+1]);
c.Label.String = strcat('RTM number in chronological order');
c.Label.FontSize = 11;
text(ax,0.85,0.1,'HF','FontSize',12,'unit','normalized');
text(ax,0.04,0.91,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.85,0.91,'PGC','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
% ax.XAxisLocation = 'top';
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
axis(ax,[xran yran]);
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

% subplot 2 
ywid = f11.ax(1).Position(4)+f11.ax(1).Position(2)-(c.Position(2)-0.00);
set(f11.ax(2), 'position', [ 0.55, c.Position(2)-0.00, 0.35, ywid]);

ax = f11.ax(2);
hold(ax,'on');
ax.FontSize = 9;
patarea = [0 0;
           -100 0;
           -100 1;
            0 1;
            0 0];
patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
histogram(ax,vdist002,'binwidth',0.5,'normalization','probability','facec',...
          [0.6 0.6 0.6],'facea',0.8);
% plot(ax,-50:0.05:50,pdffitswne,'r-','linew',1);
plot(ax,[median(vdist002) median(vdist002)],[-100 100],'b--','linew',1.5);
% plot(ax,[muswne muswne],[-100 100],'b-','linew',1.5);

ax.FontSize = 9;
text(ax,0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
ax.Box = 'on';
grid(ax, 'on');
ax.GridLineStyle = '--';
ylim(ax,[0 0.24]);
xlim(ax,[-6 6]);
% xticks(ax,-4:1:4);
% yticks(ax,[0.05 0.1 0.15]);
xlabel(ax,'LF lags or leads relative to HF (km)','fontsize',11);
ylabel(ax,'Fraction','fontsize',11);
hold(ax,'off');

%%% save figure
print(f11.fig,'-dpdf',strcat(rstpath,'/',fam,'vertical_dist_SWNE002_migs.pdf'));


%% summary of results
CIFcn = @(x,p)std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * tinv(abs([0,1]-(1-p/100)/2),...
              sum(~isnan(x(:)))-1) + mean(x(:),'omitnan');
f5.fig=figure;
set(f5.fig,'Position',[scrsz(3)/10 2*scrsz(4)/10 2*scrsz(3)/6 2*scrsz(4)/4]);

ax = gca;
hold(ax,'on');
patarea = [0 0;
           -5 0;
           -5 8;
            0 8;
            0 0];
patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
% plot(ax,[0 0],[-100 100],'--','color',[0.5 0.5 0.5]);
for i = 1: 6
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
    
    offset = (intcptproplf(i)-intcptprophf(i));
    if i==1
        ypos = i;
        e1 = errorbar(ax,offset,ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor','k');
    elseif i>1 && i~=5
        ypos = i+1;
        errorbar(ax,offset,ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor','k');
    elseif i==5
        ypos = i+1;
        e3 = errorbar(ax,offset,ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[1 96/255 0]);
    end
%     text(ax,offset,ypos+0.3,strcat('HFs: ',num2str(hfN),'; LFs: ',num2str(lfN)),'fontsize',10,...
%         'HorizontalAlignment','center');
    text(ax,4.8,ypos,strcat(num2str(hfN),'; ',num2str(lfN)),'fontsize',9,...
        'HorizontalAlignment','right');
    %%% next time use 'units', 'normalized' instead of defaults; also 
end

i=7;
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

offset = 1/2*((intcptp1lf-intcptp1hf)+(intcptp2lf-intcptp2hf));
e2 = errorbar(ax,offset,2,lfci95n(1),lfci95n(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k');
% text(ax,offset,2+0.3,strcat('HFs: ',num2str(hfN),'; LFs: ',num2str(lfN)),'fontsize',8,...
%         'HorizontalAlignment','center');
text(ax,4.8,2,strcat(num2str(hfN),'; ',num2str(lfN)),'fontsize',9,...
        'HorizontalAlignment','right');

text(ax,0.02,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);    
text(ax,0.7,0.93,'PGC','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% text(ax,0.25,0.25,'LF','FontSize',11,'unit','normalized','HorizontalAlignment','center','fontw','bold');
% text(ax,0.25,0.18,'lags','FontSize',11,'unit','normalized','HorizontalAlignment','center','fontw','bold');
% text(ax,0.75,0.25,'LF','FontSize',11,'unit','normalized','HorizontalAlignment','center','fontw','bold');
% text(ax,0.75,0.18,'leads','FontSize',11,'unit','normalized','HorizontalAlignment','center','fontw','bold');
yran = [0,8];
xran = [-5,5];
drawArrow(ax,[0.6 -0.4 ],[7 7],xran,yran,'linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);

indcor = [21,22,23,25];    
offcor = [-1.2733e+00  -2.1203e+00   2.1203e+00;
        -6.6332e-01  -1.0538e+00   1.0538e+00;
        -1.9470e+00  -1.1471e+00   1.1471e+00;
        -2.2395e-01  -4.3540e-01   4.3540e-01];

% offcor = [-1.2407e+00  -1.6066e+00   1.6066e+00;
%   -9.5821e-01  -1.1124e+00   1.1124e+00;
%   -9.2154e-01  -1.1331e+00   1.1331e+00;
%   -2.7611e-01  -4.4252e-01   4.4252e-01];

num = [166 36;
       123 61;
       51 42;
       125 84];    
   
errorbar(ax,offcor(1,1),1.3,offcor(1,2),offcor(1,3),'horizontal','o','markersize',4.5,'color',...
         'b','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6 1 1]);
errorbar(ax,offcor(2,1),4.3,offcor(2,2),offcor(2,3),'horizontal','o','markersize',4.5,'color',...
         'b','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6 1 1]);
errorbar(ax,offcor(3,1),5.3,offcor(3,2),offcor(3,3),'horizontal','o','markersize',4.5,'color',...
         'b','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[0.6 1 1]);   
errorbar(ax,offcor(4,1),6.3,offcor(4,2),offcor(4,3),'horizontal','o','markersize',4.5,'color',...
         'b','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[1 96/255 0]);

text(ax,offcor(1,1),1.5,strcat(num2str(indcor(1))),'fontsize',9,...
        'HorizontalAlignment','left');     
text(ax,offcor(2,1),4.5,strcat(num2str(indcor(2))),'fontsize',9,...
        'HorizontalAlignment','left');
text(ax,offcor(3,1),5.5,strcat(num2str(indcor(3))),'fontsize',9,...
        'HorizontalAlignment','left');    
text(ax,offcor(4,1),6.5,strcat(num2str(indcor(4))),'fontsize',9,...
        'HorizontalAlignment','left');
     
text(ax,offcor(1,1)+offcor(1,3)+0.1,1.3,strcat(num2str(num(1,1)),'; ',num2str(num(1,2))),'fontsize',9,...
        'HorizontalAlignment','left');     
text(ax,offcor(2,1)+offcor(2,3)+0.1,4.3,strcat(num2str(num(2,1)),'; ',num2str(num(2,2))),'fontsize',9,...
        'HorizontalAlignment','left');
text(ax,offcor(3,1)+offcor(3,3)+0.1,5.3,strcat(num2str(num(3,1)),'; ',num2str(num(3,2))),'fontsize',9,...
        'HorizontalAlignment','left');    
text(ax,offcor(4,1)+offcor(4,3)+0.1,6.3,strcat(num2str(num(4,1)),'; ',num2str(num(4,2))),'fontsize',9,...
        'HorizontalAlignment','left');
    
% text(ax,0.8,7,'I','FontSize',13,'VerticalAlignment','middle','backgroundcolor','w',...
%      'Margin',1);
xlabel(ax,'Offset between HF and LF (LF-HF) (km)','fontsize',11);
ylabel(ax,'RTM number in chronological order','fontsize',11);
legend(ax,[e1, e2, e3],{'SE propagation','SE (2 portions)','SW propagation'},'FontSize',8,...
        'Position',[0.2 0.70 0.15 0.08],'edgecolor',[0.5 0.5 0.5]);
% ylim(ax,[0,8]);
% xlim(ax,[-5,5]);
xticks(ax,-5:1:5);
yticks(ax,0:1:8);
ax.Box='on';
grid(ax,'on');
set(ax,'GridLineStyle','--');
ax.FontSize = 9;     
     
%%% save figure
print(f5.fig,'-dpdf',strcat(rstpath,'/',fam,'mig.lfit.sum.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));


%% plot the migration direction distribution
theta = angbest;
theta = deg2rad(theta);
figure
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
ax.RTick = 0:1:3;
ax.LineWidth = 1.5;
ax.Box = 'on';
text(0,1.8,'N','HorizontalAlignment',"center",'FontSize',14);
text(deg2rad(90),1.8,'E','HorizontalAlignment',"center",'FontSize',14);
text(deg2rad(180),1.8,'S','HorizontalAlignment',"center",'FontSize',14);
text(deg2rad(270),1.8,'W','HorizontalAlignment',"center",'FontSize',14);
text(deg2rad(315),3,'b','HorizontalAlignment',"center",'FontSize',12,'EdgeColor','k','Margin',2);
text(deg2rad(45),2.8,'PGC','HorizontalAlignment',"center",'FontSize',14,'EdgeColor','k','Margin',2);

%%% save figure
print('-dpdf',strcat(rstpath,'/',fam,'mig.direc.distri.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% combine the above two figures into one for the paper
CIFcn = @(x,p)std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * tinv(abs([0,1]-(1-p/100)/2),...
              sum(~isnan(x(:)))-1) + mean(x(:),'omitnan');
[scrsz, res] = pixelperinch(1);
f.fig=figure;
f.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 3;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.45, 0.85]);
set(f.ax(2), 'position', [ 0.6, 0.1, 0.35, 0.85]);

ax = f.ax(1);
hold(ax,'on');
patarea = [0 0;
           -5 0;
           -5 8;
            0 8;
            0 0];
patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
% plot(ax,[0 0],[-100 100],'--','color',[0.5 0.5 0.5]);
for i = 1: 6
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
    
    offset = (intcptproplf(i)-intcptprophf(i));
    if i==1
        ypos = i;
        e1 = errorbar(ax,offset,ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor','k');
    elseif i>1 && i~=5
        ypos = i+1;
        errorbar(ax,offset,ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor','k');
    elseif i==5
        ypos = i+1;
        e3 = errorbar(ax,offset,ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[1 96/255 0]);
    end
%     text(ax,offset,ypos+0.3,strcat('HFs: ',num2str(hfN),'; LFs: ',num2str(lfN)),'fontsize',10,...
%         'HorizontalAlignment','center');
    text(ax,4.8,ypos,strcat(num2str(hfN),'; ',num2str(lfN)),'fontsize',9,...
        'HorizontalAlignment','right');
    %%% next time use 'units', 'normalized' instead of defaults; also 
end

i=7;
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

offset = 1/2*((intcptp1lf-intcptp1hf)+(intcptp2lf-intcptp2hf));
e2 = errorbar(ax,offset,2,lfci95n(1),lfci95n(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k');
% text(ax,offset,2+0.3,strcat('HFs: ',num2str(hfN),'; LFs: ',num2str(lfN)),'fontsize',8,...
%         'HorizontalAlignment','center');
text(ax,4.8,2,strcat(num2str(hfN),'; ',num2str(lfN)),'fontsize',9,...
        'HorizontalAlignment','right');

text(ax,0.02,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);    
text(ax,0.7,0.93,'PGC','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.25,0.25,'LF','FontSize',11,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(ax,0.25,0.18,'lags','FontSize',11,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(ax,0.75,0.25,'LF','FontSize',11,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(ax,0.75,0.18,'leads','FontSize',11,'unit','normalized','HorizontalAlignment','center','fontw','bold');
yran = [0,8];
xran = [-5,5];
drawArrow(ax,[0.6 -0.4 ],[7 7],xran,yran,'linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);

text(ax,0.8,7,'I','FontSize',13,'VerticalAlignment','middle','backgroundcolor','w',...
     'Margin',1);
xlabel(ax,'Offset between HF and LF (LF-HF) (km)','fontsize',11);
ylabel(ax,'RTM number in chronological order','fontsize',11);
legend(ax,[e1, e2, e3],{'SE propagation','SE (2 portions)','SW propagation'},'FontSize',8,...
        'Position',[0.12 0.70 0.15 0.08],'edgecolor',[0.5 0.5 0.5]);
% ylim(ax,[0,8]);
% xlim(ax,[-5,5]);
xticks(ax,-5:1:5);
yticks(ax,0:1:8);
ax.Box='on';
grid(ax,'on');
set(ax,'GridLineStyle','--');
ax.FontSize = 9;


ax = f.ax(2);
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
ax=gca;
ax.RTickLabelRotation = 45;
ax.ThetaMinorTick = 'on';
ax.TickLength = [0.02 0];
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.FontSize = 9;
ax.RTick = 0:1:3;
% ax.LineWidth = 1;
ax.Box = 'on';
ax.GridAlpha = 0.3;
text(ax,deg2rad(140),1.3,num2str(sum(angbest>90 & angbest<=180)),'HorizontalAlignment','center',...
     'FontSize',11);           
text(ax,deg2rad(205),1.3,num2str(sum(angbest>180 & angbest<=270)),'HorizontalAlignment','center',...
     'FontSize',11);           
text(ax,0,1.8,'N','HorizontalAlignment',"center",'FontSize',12);
text(ax,deg2rad(90),1.8,'E','HorizontalAlignment',"center",'FontSize',12);
text(ax,deg2rad(180),1.8,'S','HorizontalAlignment',"center",'FontSize',12);
text(ax,deg2rad(270),1.8,'W','HorizontalAlignment',"center",'FontSize',12);
text(ax,deg2rad(307),2.9,'b','HorizontalAlignment',"center",'FontSize',11,'EdgeColor','k','Margin',2);
% text(ax,deg2rad(45),2.8,'PGC','HorizontalAlignment',"center",'FontSize',12,'EdgeColor','k','Margin',2);


%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/',fam,'sum.mig.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% choose one representative RTM and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trange1 can be fitted by one slope
trange1 = [
          2005256,7.62e+4,8.30e+4;];     %2004197,8.45e+4,8.55e+4;

angbestl1 = zeros(length(trange1),1);
angbestl2 = zeros(length(trange1),1);
angbestl3 = zeros(length(trange1),1);
angbestl4 = zeros(length(trange1),1);
angbestl5 = zeros(length(trange1),1);

% xran = [ -20 12];
% yran = [ -20 12];
xran = [-10 30];
yran = [-30 10];

resprophf = nan(length(trange1)+1,200);
resproplf = nan(length(trange1)+1,50);
for i = 1: size(trange1,1)
%     i=3;
    disp(trange1(i,:));
    indhf = find(hfrelalocall(:,12)==trange1(i,1) & hfrelalocall(:,14)>=trange1(i,2) & ...
                 hfrelalocall(:,14)<=trange1(i,3));
    mighf = hfrelalocall(indhf,:);
    
    angle = 0:5:360;
    
    l1norm = zeros(length(angle),1);
    l2norm = zeros(length(angle),1);
    slope = zeros(length(angle),1);
    sse = zeros(length(angle),1);
    rmse = zeros(length(angle),1);
    rsquare = zeros(length(angle),1);
    for iang = 1: length(angle)
        mighfdum = mighf;
        for j = 1: size(mighf,1)
            x0 = mighfdum(j,3);
            y0 = mighfdum(j,4);
            [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),0,0);
            mighfdum(j,3) = newx;
            mighfdum(j,4) = newy;
        end
        
        %%% linear robust least square
        [fitobj,gof,output] = fit(mighfdum(:,14)/3600, mighfdum(:,3),fttpfree,'Robust','Bisquare',...
                                  'StartPoint',[1 1]);
        % output fit parameters
        coef = coeffvalues(fitobj);
        slope(iang) = coef(1);
        fitprophf = feval(fitobj,mighfdum(:,14)/3600);
        l1norm(iang) = sum(abs(mighfdum(:,3)-fitprophf))/(length(mighfdum(:,3)));
        l2norm(iang) = sum((mighfdum(:,3)-fitprophf).^2)/(length(mighfdum(:,3)));
        sse(iang) = gof.sse;
        rmse(iang) = gof.rmse;
        rsquare(iang) = gof.rsquare;    % r-square, defined as ratio of the sum of squares of the regression and the total sum of squares
    
    end
    ind = find(slope>0);
    ind1 = find(l1norm(ind)==min(l1norm(ind)));
    angbestl1(i) = angle(ind(ind1(1)));
    
    ind2 = find(l2norm(ind)==min(l2norm(ind)));
    angbestl2(i) = angle(ind(ind2(1)));
    
    ind3 = find(rmse(ind)==min(rmse(ind)));     % rmse, Root Mean Squared Error, the smaller, the better
    angbestl3(i) = angle(ind(ind3(1)));
    
    ind4 = find(rsquare(ind)==max(rsquare(ind)));   % R-square, 0-1, the larger the better
    angbestl4(i) = angle(ind(ind4(1)));
    
    ind5 = find(sse(ind)==min(sse(ind)));   % sum of square due to error, the smaller, the better
    angbestl5(i) = angle(ind(ind5(1)));
      
    angbest(i) = median([angbestl3(i),angbestl4(i),angbestl5(i)]); 
    
    mighfdum = mighf;
    for j = 1: size(mighf,1)
        x0 = mighfdum(j,3);
        y0 = mighfdum(j,4);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        mighfdum(j,3) = newx;
        mighfdum(j,4) = newy;
    end
    
    indlf = find(lfrelalocall(:,12)==trange1(i,1) & lfrelalocall(:,14)>=trange1(i,2) & ...
                 lfrelalocall(:,14)<=trange1(i,3));
    miglf = lfrelalocall(indlf,:); 
    miglfdum = miglf;
    for j = 1: size(miglf,1)
        x0 = miglfdum(j,3);
        y0 = miglfdum(j,4);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        miglfdum(j,3) = newx;
        miglfdum(j,4) = newy;
    end
    
    %%% define and position the figure frame and axes of each plot
    f.fig = figure;
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
    mighf = sortrows(mighf,-14);
    scatter(f.ax(1),mighf(:,3),mighf(:,4), 20, mighf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(1),'jet');
    c=colorbar(f.ax(1),'SouthOutside');
    pos = f.ax(1).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    juldate = num2str(trange1(i,1));
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
    caxis(f.ax(1),[trange1(i,2)/3600 trange1(i,3)/3600])
    text(f.ax(1),0.85,0.1,'HF','FontSize',12,'unit','normalized');
    text(f.ax(1),0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(1),0.04,0.1,'PGC','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(1).Box = 'on';
    grid(f.ax(1), 'on');
    axis(f.ax(1), 'equal');
    f.ax(1).GridLineStyle = '--';
    f.ax(1).XAxisLocation = 'top';
    medxhf = median(mighf(:,3));
    medyhf = median(mighf(:,4));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medxhf-rotx medxhf+rotx];
    yvect = [medyhf-roty medyhf+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medxhf-rotx medxhf+rotx];
    yvect = [medyhf-roty medyhf+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',...
              [0.4 0.4 0.4]);
    xticks(f.ax(1),xran(1):5:xran(2));
    yticks(f.ax(1),yran(1):5:yran(2)); 
    xlabel(f.ax(1),'E (km)','fontsize',11);
    ylabel(f.ax(1),'N (km)','fontsize',11);
    text(f.ax(1),0.5,0.1,'I','FontSize',12,'unit','normalized');
    text(f.ax(1),0.83,0.95,strcat(num2str(size(mighf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(1),0.83,0.90,strcat({'in '},num2str(trange1(i,3)-trange1(i,2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(mighf,1)/(trange1(i,3)-trange1(i,2)));
    text(f.ax(1),0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    hold(f.ax(1),'off');

    % subplot 2 of figure i
    hold(f.ax(2),'on');
    plot(f.ax(2),[-100 100],[0 0],'k--');
    plot(f.ax(2),[0 0],[-100 100],'k--');
    f.ax(2).FontSize = 9;
    miglf = sortrows(miglf,-14);
    scatter(f.ax(2),miglf(:,3),miglf(:,4), 20, miglf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(2),'jet');
    c=colorbar(f.ax(2),'SouthOutside');
    pos = f.ax(2).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
%     c.Label.String = strcat({'Time (hr) on 13 Sep. 2005'});
    c.Label.FontSize = 11;
    caxis(f.ax(2),[trange1(i,2)/3600 trange1(i,3)/3600])
    text(f.ax(2),0.85,0.1,'LF','FontSize',12,'unit','normalized');
    text(f.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(2),0.04,0.1,'PGC','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(2).Box = 'on';
    grid(f.ax(2), 'on');
    axis(f.ax(2), 'equal');
    f.ax(2).GridLineStyle = '--';
    f.ax(2).XAxisLocation = 'top';
    medxlf = median(miglf(:,3));
    medylf = median(miglf(:,4));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',...
              [0.4 0.4 0.4]);
    xticks(f.ax(2),xran(1):5:xran(2));
    yticks(f.ax(2),yran(1):5:yran(2)); 
    xlabel(f.ax(2),'E (km)','fontsize',11);
    ylabel(f.ax(2),'N (km)','fontsize',11); 
    text(f.ax(2),0.5,0.1,'I','FontSize',12,'unit','normalized');
    text(f.ax(2),0.83,0.95,strcat(num2str(size(miglf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(2),0.83,0.9,strcat({'in '},num2str(trange1(i,3)-trange1(i,2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(miglf,1)/(trange1(i,3)-trange1(i,2)));
    text(f.ax(2),0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    hold(f.ax(2),'off');
    
    % subplot 3 of figure i
    hold(f.ax(3),'on');
    f.ax(3).FontSize = 9;
    scatter(f.ax(3),mighfdum(:,14)/3600,mighfdum(:,3),20,[0.6 0.6 0.6],'filled','o',...
            'MarkerEdgeColor','k');      
    scatter(f.ax(3),miglfdum(:,14)/3600,miglfdum(:,3),20,[0.6 1 1],'filled','o',...
            'MarkerEdgeColor','k');
%     scatter(f.ax(3),mighfdum(:,14)/3600,mighfdum(:,3),30,[105/255 105/255 105/255],'filled','o');
%     scatter(f.ax(3),miglfdum(:,14)/3600,miglfdum(:,3),30,miglfdum(:,11),'filled','o');  

    [fitobjhfprop,gofhfprop,outphfprop] = fit(mighfdum(:,14)/3600, mighfdum(:,3),fttpfree,'Robust',...
                                              'Bisquare','StartPoint',[1 1]);
    % output fit parameters
    coefprop = coeffvalues(fitobjhfprop);
    slopeprophf(i) = coefprop(1);
    intcptprophf(i) = coefprop(2);
    fitprophf = feval(fitobjhfprop,mighfdum(:,14)/3600);
    plot(f.ax(3),mighfdum(:,14)/3600,fitprophf,'-','linewidth',2,'color',[0.6 0.6 0.6]);
        
%     intcptproplf = ones(size(miglfdum(:,end)))\(miglfdum(:,3)-slopeprophf*miglfdum(:,end)/3600); % direct LS fit
%     fitproplf = slopeprophf*miglfdum(:,end)/3600+intcptproplf;
    a=slopeprophf(i);
    fttpfix = fittype( @(b,x) a*x+b);
    [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum(:,14)/3600, miglfdum(:,3),fttpfix,'Robust',...
                                              'Bisquare','StartPoint',1);
%     [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum(:,14)/3600, miglfdum(:,3),fttpfix,'Robust',...
%                                               'Bisquare','Weights',miglfdum(:,11));
    fitproplf = feval(fitobjlfprop,miglfdum(:,14)/3600);
    intcptproplf(i) = coeffvalues(fitobjlfprop);
    offset(i) = intcptproplf(i)-intcptprophf(i);
    plot(f.ax(3),miglfdum(:,14)/3600,fitproplf,'-','linewidth',2,'color',[0.6 1 1]);
    f.ax(3).Box = 'on';
    grid(f.ax(3), 'on');
    f.ax(3).GridLineStyle = '--';
    xran1 = [trange1(i,2)/3600 trange1(i,3)/3600];
    yran1 = [round(min([mighfdum(:,3);miglfdum(:,3)]))-1 ...
             round(max([mighfdum(:,3);miglfdum(:,3)]))+1];
    xlim(f.ax(3),xran1);
    ylim(f.ax(3),yran1);
    text(f.ax(3),0.04,0.91,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    xlabel(f.ax(3),'Time (hr)','fontsize',11);
    ylabel(f.ax(3),'Dist. along prop. (km)','fontsize',11); 
    hold(f.ax(3),'off');
    
    % compute the HF weights in robust linear regression, see NOTES above
    rhf = outphfprop.residuals;   % usual residuals
    x = [mighfdum(:,14)/3600];
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
    x = [miglfdum(:,14)/3600];
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
    
    
    % subplot 4 of figure i
    hold(f.ax(4),'on');
    f.ax(4).FontSize = 9;
    numhf(i) = length(mighfdum(:,3));
    numlf(i) = length(miglfdum(:,3));
    resprophf(i,1:numhf(i)) = outphfprop.residuals;
    resproplf(i,1:numlf(i)) = outplfprop.residuals+(intcptproplf(i)-intcptprophf(i));
    [xlochf, ~, ~, pdfvalhf] = weighted_bar_pdf(resprophf(i,1:numhf(i)), wthf, 0.5, 'int');
    hfHdl = bar(f.ax(4),xlochf, pdfvalhf,1,'stacked');
    hfHdl(1).FaceColor = [0.6 0.6 0.6];
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
%     plot(min(resprophf(i,1:numhf(i))):0.05:max(resprophf(i,1:numhf(i))), pdffithf, 'r-');
%     plot(min(resproplf(i,1:numlf(i))):0.05:max(resproplf(i,1:numlf(i))), pdffitlf, '-','color',...
%          [0.6 1 1]);
    text(f.ax(4),0.04,0.91,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(4),0.4,0.92,sprintf('LF-HF = %.2f km',offset(i)),'FontSize',11,'unit','normalized');
    f.ax(4).Box = 'on';
    grid(f.ax(4), 'on');
    f.ax(4).GridLineStyle = '--';
    ymax = f.ax(4).YLim(2)+0.1;
    ylim(f.ax(4),[0 ymax]);
    xlabel(f.ax(4),'Residual in prop. (km)','fontsize',11);
    ylabel(f.ax(4),'PDF estimate','fontsize',11);
    hold(f.ax(4),'off');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% It seems that median of the previous detections is a promising way to show offsets as well,
    %%% so try to deal with it the same way
    medhf = [];
    for j = 1: size(mighfdum,1)
        medhf(j,1) = median(mighfdum(1:j,3));
        medhf(j,2) = median(mighfdum(1:j,4));
        medhf(j,3) = mighfdum(j,14);
    end
    
    medlf = [];
    for j = 1: size(miglfdum,1)
        medlf(j,1) = median(miglfdum(1:j,3));
        medlf(j,2) = median(miglfdum(1:j,4));
        medlf(j,3) = miglfdum(j,14);
    end
       
    
    % subplot 5 of figure i
    hold(f.ax(5),'on');
    f.ax(5).FontSize = 9;
    scatter(f.ax(5),medhf(:,3)/3600,medhf(:,1),20,[0.6 0.6 0.6],'filled','o','MarkerEdgeColor',...
            'k');
    scatter(f.ax(5),medlf(:,3)/3600,medlf(:,1),20,[0.6 1 1],'filled','o','MarkerEdgeColor',...
            'k');            
    f.ax(5).Box = 'on';
    grid(f.ax(5), 'on');
    f.ax(5).GridLineStyle = '--';
    xran1 = [trange1(i,2)/3600 trange1(i,3)/3600];
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
    scatter(f.ax(6),medhf(:,3)/3600,medhf(:,2),20,[0.6 0.6 0.6],'filled','o','MarkerEdgeColor',...
            'k');
    scatter(f.ax(6),medlf(:,3)/3600,medlf(:,2),20,[0.6 1 1],'filled','o','MarkerEdgeColor',...
            'k');          
    f.ax(6).Box = 'on';
    grid(f.ax(6), 'on');
    f.ax(6).GridLineStyle = '--';
    xran1 = [trange1(i,2)/3600 trange1(i,3)/3600];
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
    print(f.fig,'-dpdf',strcat(rstpath,'/',fam,'sel.mig.proj.lfit.meddist',num2str(trange1(i,1)),...
          '_',num2str(trange1(i,2)),'-',num2str(trange1(i,3)),'.',num2str(winlenhf),'_',...
          num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));
end    

 
















