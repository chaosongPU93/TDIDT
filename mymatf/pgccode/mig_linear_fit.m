% function mig_linear_fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to find the best propagation direction for migration, and then fit the hf
% the migration linearly to get the propagation slope and error distribution.
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

xran = [-16 12; -10 10; -21 7; -7 11; -21 13; -21 11];
yran = [-17 13; -2 10; -17 12; -10 11; -23 9; -20 12];


% create fit object
fttpfree = fittype( @(a,b,x) a*x+b);
        
for i = 1: length(trange1)
%     i=3;
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
            x0 = mighfdum(j,1);
            y0 = mighfdum(j,2);
            [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),0,0);
            mighfdum(j,1) = newx;
            mighfdum(j,2) = newy;
        end
        
        %%% linear robust least square
        [fitobj,gof,output] = fit(mighfdum(:,14)/3600, mighfdum(:,1),fttpfree,'Robust','Bisquare');
        % output fit parameters
        coef = coeffvalues(fitobj);
        slope(iang) = coef(1);
        fitprophf = feval(fitobj,mighfdum(:,14)/3600);
        l1norm(iang) = sum(abs(mighfdum(:,1)-fitprophf))/(length(mighfdum(:,1)));
        l2norm(iang) = sum((mighfdum(:,1)-fitprophf).^2)/(length(mighfdum(:,1)));
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
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        mighfdum(j,1) = newx;
        mighfdum(j,2) = newy;
    end
    
    indlf = find(lfrelalocall(:,8)==trange1(i,1) & lfrelalocall(:,10)>=trange1(i,2) & ...
                 lfrelalocall(:,10)<=trange1(i,3));
    miglf = lfrelalocall(indlf,:); 
    miglfdum = miglf;
    for j = 1: size(miglf,1)
        x0 = miglfdum(j,1);
        y0 = miglfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        miglfdum(j,1) = newx;
        miglfdum(j,2) = newy;
    end
    
    %%% define and position the figure frame and axes of each plot
    f4.fig=figure(2*i-1);
    set(f4.fig,'Position',[scrsz(3)/10 0*scrsz(4)/10 2*scrsz(3)/5 4*scrsz(4)/5]);
    nrow = 2;
    ncol = 2;    
    for isub = 1:nrow*ncol
        f4.ax(isub) = subplot(nrow,ncol,isub);
    end

    %%% reposition
    xlen = xran(i,2)-xran(i,1);
    ylen = yran(i,2)-yran(i,1);
    subfxlen = 0.35;    
    if xlen >= ylen
        subfylen = ylen/xlen*0.4;
    else
        subfylen = ylen/xlen*0.35;
    end
    set(f4.ax(1), 'position', [ 0.1, 0.5, subfxlen, subfylen]);
    set(f4.ax(2), 'position', [ 0.55, 0.5, subfxlen, subfylen]);
    set(f4.ax(3), 'position', [ 0.1, 0.2, 0.35, 0.25]);
    set(f4.ax(4), 'position', [ 0.55, 0.2, 0.35, 0.25]);
    
    % subplot 1 of figure i
    hold(f4.ax(1),'on');
    plot(f4.ax(1),[-100 100],[0 0],'k--');
    plot(f4.ax(1),[0 0],[-100 100],'k--');
    mighf = sortrows(mighf,-14);
    scatter(f4.ax(1),mighf(:,1),mighf(:,2), 30, mighf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f4.ax(1),'jet');
    c=colorbar(f4.ax(1),'SouthOutside');
    pos = f4.ax(1).Position;
    c.Position = [pos(1), pos(2), pos(3), 0.02];
%     c.TickLabels=[];
    c.Label.String = strcat(num2str(trange1(i,1)),' of HF',' (hr)');
    c.Label.FontSize = 12;
    caxis(f4.ax(1),[trange1(i,2)/3600 trange1(i,3)/3600])
    text(f4.ax(1),0.85,0.93,'HF','FontSize',15,'unit','normalized');
    text(f4.ax(1),0.04,0.93,'a','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
    f4.ax(1).Box = 'on';
    grid(f4.ax(1), 'on');
    axis(f4.ax(1), 'equal');
    f4.ax(1).GridLineStyle = '--';
    f4.ax(1).XAxisLocation = 'top';
    medx = median(mighf(:,1));
    medy = median(mighf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
    drawArrow(f4.ax(1),xvect,yvect,xran(i,:),yran(i,:),'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
    drawArrow(f4.ax(1),xvect,yvect,xran(i,:),yran(i,:),'linewidth',1.5,'linestyle','--','color',...
              [0.4 0.4 0.4]);
    xlabel(f4.ax(1),'E (km)','fontsize',12);
    ylabel(f4.ax(1),'N (km)','fontsize',12); 
    hold(f4.ax(1),'off');

    % subplot 2 of figure i
    hold(f4.ax(2),'on');
    plot(f4.ax(2),[-100 100],[0 0],'k--');
    plot(f4.ax(2),[0 0],[-100 100],'k--');
    miglf = sortrows(miglf,-14);
    scatter(f4.ax(2),miglf(:,1),miglf(:,2), 30, miglf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f4.ax(2),'jet');
    c=colorbar(f4.ax(2),'SouthOutside');
    pos = f4.ax(2).Position;
    c.Position = [pos(1), pos(2), pos(3), 0.02];
%     c.TickLabels=[];
    c.Label.String = strcat(num2str(trange1(i,1)),' of LF',' (hr)');
    c.Label.FontSize = 12;
    caxis(f4.ax(2),[trange1(i,2)/3600 trange1(i,3)/3600])
    text(f4.ax(2),0.85,0.93,'LF','FontSize',15,'unit','normalized');
    text(f4.ax(2),0.04,0.93,'b','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
    f4.ax(2).Box = 'on';
    grid(f4.ax(2), 'on');
    axis(f4.ax(2), 'equal');
    f4.ax(2).GridLineStyle = '--';
    f4.ax(2).XAxisLocation = 'top';
    medx = median(miglf(:,1));
    medy = median(miglf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
    drawArrow(f4.ax(2),xvect,yvect,xran(i,:),yran(i,:),'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
    drawArrow(f4.ax(2),xvect,yvect,xran(i,:),yran(i,:),'linewidth',1.5,'linestyle','--','color',...
              [0.4 0.4 0.4]);
    xlabel(f4.ax(2),'E (km)','fontsize',12);
    ylabel(f4.ax(2),'N (km)','fontsize',12); 
    hold(f4.ax(2),'off');
    
    % subplot 3 of figure i
    hold(f4.ax(3),'on');
    scatter(f4.ax(3),mighfdum(:,14)/3600,mighfdum(:,1),30,'filled','ro','MarkerEdgeColor', 'k');   
    scatter(f4.ax(3),miglfdum(:,14)/3600,miglfdum(:,1),30,[0.6 1 1],'filled','o',...
            'MarkerEdgeColor','k');  %, 'MarkerEdgeColor', 'w')
%     scatter(f.ax(3),mighfdum(:,14)/3600,mighfdum(:,1),30,[105/255 105/255 105/255],'filled','o');
%     scatter(f.ax(3),miglfdum(:,14)/3600,miglfdum(:,1),30,miglfdum(:,11),'filled','o');  
    % linear robust least square
    fttpfree = fittype( @(a,b,x) a*x+b);
    % create fit object
    [fitobjhfprop,gofhfprop,outphfprop] = fit(mighfdum(:,14)/3600, mighfdum(:,1),fttpfree,'Robust',...
                                              'Bisquare');
    % output fit parameters
    coefprop = coeffvalues(fitobjhfprop);
    slopeprophf(i) = coefprop(1);
    intcptprophf(i) = coefprop(2);
    fitprophf = feval(fitobjhfprop,mighfdum(:,14)/3600);
    plot(f4.ax(3),mighfdum(:,14)/3600,fitprophf,'r-','linewidth',2);
        
%     intcptproplf = ones(size(miglfdum(:,end)))\(miglfdum(:,1)-slopeprophf*miglfdum(:,end)/3600); % direct LS fit
%     fitproplf = slopeprophf*miglfdum(:,end)/3600+intcptproplf;
    a=slopeprophf(i);
    fttpfix = fittype( @(b,x) a*x+b);
    [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum(:,14)/3600, miglfdum(:,1),fttpfix,'Robust',...
                                              'Bisquare');
%     [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum(:,14)/3600, miglfdum(:,1),fttpfix,'Robust',...
%                                               'Bisquare','Weights',miglfdum(:,11));
    fitproplf = feval(fitobjlfprop,miglfdum(:,14)/3600);
    intcptproplf(i) = coeffvalues(fitobjlfprop);
    offset(i) = intcptproplf(i)-intcptprophf(i);
    plot(f4.ax(3),miglfdum(:,14)/3600,fitproplf,'-','linewidth',2,'color',[0.6 1 1]);
    f4.ax(3).Box = 'on';
    grid(f4.ax(3), 'on');
    f4.ax(3).GridLineStyle = '--';
    xran1 = [trange1(i,2)/3600 trange1(i,3)/3600];
    yran1 = [round(min([mighfdum(:,1);miglfdum(:,1)]))-1 ...
             round(max([mighfdum(:,1);miglfdum(:,1)]))+1];
    xlim(f4.ax(3),xran1);
    ylim(f4.ax(3),yran1);
    text(f4.ax(3),0.04,0.92,'c','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
    xlabel(f4.ax(3),'Time (hr)','fontsize',12);
    ylabel(f4.ax(3),'Dist. along propagation (km)','fontsize',12); 
    hold(f4.ax(3),'off');
    
    % subplot 4 of figure i
    hold(f4.ax(4),'on');
    numhf(i) = length(mighfdum(:,1));
    numlf(i) = length(miglfdum(:,1));
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
    histogram(f4.ax(4),resprophf(i,1:numhf(i)),'binwidth',0.5,'normalization','pdf','facecolor','r');
    histogram(f4.ax(4),resproplf(i,1:numlf(i)),'binwidth',0.5,'normalization','pdf','facecolor',...
              [0.6 1 1],'facealpha',0.6);
%     plot(min(resprophf(i,1:numhf(i))):0.05:max(resprophf(i,1:numhf(i))), pdffithf, 'r-');
%     plot(min(resproplf(i,1:numlf(i))):0.05:max(resproplf(i,1:numlf(i))), pdffitlf, '-','color',...
%          [0.6 1 1]);
    text(f4.ax(4),0.04,0.92,'d','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
    text(f4.ax(4),0.50,0.92,sprintf('lf-hf = %.2f km',offset(i)),'FontSize',12,'unit','normalized');
    f4.ax(4).Box = 'on';
    grid(f4.ax(4), 'on');
    f4.ax(4).GridLineStyle = '--';
    ylim(f4.ax(4),[0 1]);
    xlabel(f4.ax(4),'Residual in propagation (km)','fontsize',12);
    ylabel(f4.ax(4),'PDF estimate','fontsize',12);
    hold(f4.ax(4),'off');
    
%%%%%%%%%%% for orthogonal direction, ignored now %%%%%%%%%%%    
%     % subplot 5 of figure i
%     hold(f.ax(4),'on');
%     alpha(f.ax(4),0.1);
%     scatter(f.ax(4),mighfdum(:,end)/3600,mighfdum(:,2),30,'filled','ko');  %, 'MarkerEdgeColor', 'w')
%     scatter(f.ax(4),miglfdum(:,end)/3600,miglfdum(:,2),30,[180/255 180/255 180/255],'filled','o');  %, 'MarkerEdgeColor', 'w')
%     % linear robust least square
%     % create fit object
%     [fitobjort,gofort,~] = fit(mighfdum(:,end)/3600, mighfdum(:,2), 'poly1',fitopt);
%     % output fit parameters
%     coefort = coeffvalues(fitobjort);
%     slopeorthf = coefort(1);
%     intcptorthf = coefort(2);
%     fitorthf = feval(fitobjort,mighfdum(:,end)/3600);
%     plot(f.ax(4),mighfdum(:,end)/3600,fitorthf,'r-','linewidth',2);
%     intcptortlf = ones(size(miglfdum(:,end)))\(miglfdum(:,2)-slopeorthf*miglfdum(:,end)/3600);
%     fitortlf = slopeorthf*miglfdum(:,end)/3600+intcptortlf;
%     plot(f.ax(4),miglfdum(:,end)/3600,fitortlf,'b-','linewidth',2);
%     f.ax(4).Box = 'on';
%     grid(f.ax(4), 'on');
%     f.ax(4).GridLineStyle = '--';
%     xran = [trange1(i,2)/3600 trange1(i,3)/3600];
%     yran = [round(min(mighfdum(:,1)))-1 round(max(mighfdum(:,1)))+1];
%     xlim(f.ax(4),xran);
%     ylim(f.ax(4),yran);
%     xlabel(f.ax(4),'Time (hr)','fontsize',12);
%     ylabel(f.ax(4),'Distance along orthogonal (km)','fontsize',12);
%     hold(f.ax(4),'off');
%     
%     % subplot 6 of figure i
%     hold(f.ax(6),'on');
%     resorthf = mighfdum(:,2)-fitorthf;
%     resortlf = miglfdum(:,2)-fitortlf;
%     histogram(f.ax(6),resorthf,'binwidth',1,'normalization','pdf','facecolor','k');
%     histogram(f.ax(6),resortlf,'binwidth',1,'normalization','pdf','facecolor',[180/255 180/255 180/255],'facealpha',0.8);
%     f.ax(6).Box = 'on';
%     grid(f.ax(6), 'on');
%     f.ax(6).GridLineStyle = '--';
%     xlabel(f.ax(6),'Residual in orthogonal (km)','fontsize',12);
%     ylabel(f.ax(6),'PDF','fontsize',12);
%     hold(f.ax(6),'off');
%%%%%%%%%%% for orthogonal direction, ignored now %%%%%%%%%%%


    %%% save figure
    print(f4.fig,'-dpdf',strcat(rstpath,'/',fam,'mig.proj.lfit',num2str(trange1(i,1)),...
          '_',num2str(trange1(i,2)),'-',num2str(trange1(i,3)),'.',num2str(winlenhf),'_',...
          num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% It seems that median of the previous detections is a promising way to show offsets as well,
    %%% so try to deal with it the same way
    for j = 1: size(mighf,1)
        medhf(j,1) = median(mighfdum(1:j,1));
        medhf(j,2) = median(mighfdum(1:j,2));
        medhf(j,3) = mighfdum(j,14);
    end
    
    for j = 1: size(miglf,1)
        medlf(j,1) = median(miglfdum(1:j,1));
        medlf(j,2) = median(miglfdum(1:j,2));
        medlf(j,3) = miglfdum(j,14);
    end
    
    %%% define and position the figure frame and axes of each plot
    f2.fig=figure;
    set(f2.fig,'Position',[scrsz(3)/10 0*scrsz(4)/10 2*scrsz(3)/5 4*scrsz(4)/5]);
    nrow = 2;
    ncol = 2;    
    for isub = 1:nrow*ncol
        f2.ax(isub) = subplot(nrow,ncol,isub);
    end

    %%% reposition
    xlen = xran(i,2)-xran(i,1);
    ylen = yran(i,2)-yran(i,1);
    subfxlen = 0.35;    
    if xlen >= ylen
        subfylen = ylen/xlen*0.4;
    else
        subfylen = ylen/xlen*0.35;
    end
    set(f2.ax(1), 'position', [ 0.1, 0.5, subfxlen, subfylen]);
    set(f2.ax(2), 'position', [ 0.55, 0.5, subfxlen, subfylen]);
    set(f2.ax(3), 'position', [ 0.1, 0.2, 0.35, 0.25]);
    set(f2.ax(4), 'position', [ 0.55, 0.2, 0.35, 0.25]);
    
    % subplot 1 of figure i
    hold(f2.ax(1),'on');
    plot(f2.ax(1),[-100 100],[0 0],'k--');
    plot(f2.ax(1),[0 0],[-100 100],'k--');
    mighf = sortrows(mighf,-14);
    scatter(f2.ax(1),mighf(:,1),mighf(:,2), 30, mighf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f2.ax(1),'jet');
    c=colorbar(f2.ax(1),'SouthOutside');
    pos = f2.ax(1).Position;
    c.Position = [pos(1), pos(2), pos(3), 0.02];
%     c.TickLabels=[];
    c.Label.String = strcat(num2str(trange1(i,1)),' of HF',' (hr)');
    c.Label.FontSize = 12;
    caxis(f2.ax(1),[trange1(i,2)/3600 trange1(i,3)/3600])
    text(f2.ax(1),0.85,0.93,'HF','FontSize',15,'unit','normalized');
    text(f2.ax(1),0.04,0.93,'a','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
    f2.ax(1).Box = 'on';
    grid(f2.ax(1), 'on');
    axis(f2.ax(1), 'equal');
    f2.ax(1).GridLineStyle = '--';
    f2.ax(1).XAxisLocation = 'top';
    medx = median(mighf(:,1));
    medy = median(mighf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
    drawArrow(f2.ax(1),xvect,yvect,xran(i,:),yran(i,:),'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
    drawArrow(f2.ax(1),xvect,yvect,xran(i,:),yran(i,:),'linewidth',1.5,'linestyle','--','color',...
              [0.4 0.4 0.4]);
    xlabel(f2.ax(1),'E (km)','fontsize',12);
    ylabel(f2.ax(1),'N (km)','fontsize',12); 
    hold(f2.ax(1),'off');

    % subplot 2 of figure i
    hold(f2.ax(2),'on');
    plot(f2.ax(2),[-100 100],[0 0],'k--');
    plot(f2.ax(2),[0 0],[-100 100],'k--');
    miglf = sortrows(miglf,-14);
    scatter(f2.ax(2),miglf(:,1),miglf(:,2), 30, miglf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f2.ax(2),'jet');
    c=colorbar(f2.ax(2),'SouthOutside');
    pos = f2.ax(2).Position;
    c.Position = [pos(1), pos(2), pos(3), 0.02];
%     c.TickLabels=[];
    c.Label.String = strcat(num2str(trange1(i,1)),' of LF',' (hr)');
    c.Label.FontSize = 12;
    caxis(f2.ax(2),[trange1(i,2)/3600 trange1(i,3)/3600])
    text(f2.ax(2),0.85,0.93,'LF','FontSize',15,'unit','normalized');
    text(f2.ax(2),0.04,0.93,'b','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
    f2.ax(2).Box = 'on';
    grid(f2.ax(2), 'on');
    axis(f2.ax(2), 'equal');
    f2.ax(2).GridLineStyle = '--';
    f2.ax(2).XAxisLocation = 'top';
    medx = median(miglf(:,1));
    medy = median(miglf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
    drawArrow(f2.ax(2),xvect,yvect,xran(i,:),yran(i,:),'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
    drawArrow(f2.ax(2),xvect,yvect,xran(i,:),yran(i,:),'linewidth',1.5,'linestyle','--',...
              'color',[0.4 0.4 0.4]);
    xlabel(f2.ax(2),'E (km)','fontsize',12);
    ylabel(f2.ax(2),'N (km)','fontsize',12); 
    hold(f2.ax(2),'off');
    
    
    % subplot 3 of figure i
    hold(f2.ax(3),'on');
    scatter(f2.ax(3),medhf(:,3)/3600,medhf(:,1),30,'filled','ro','MarkerEdgeColor', 'k');
    scatter(f2.ax(3),medlf(:,3)/3600,medlf(:,1),30,[0.6 1 1],'filled','o','MarkerEdgeColor', 'k');            
    f2.ax(3).Box = 'on';
    grid(f2.ax(3), 'on');
    f2.ax(3).GridLineStyle = '--';
    xran1 = [trange1(i,2)/3600 trange1(i,3)/3600];
    yran1 = [round(min([medhf(:,1);medlf(:,1)]))-1 ...
             round(max([medhf(:,1);medlf(:,1)]))+1];
    xlim(f2.ax(3),xran1);
    ylim(f2.ax(3),yran1);
    text(f2.ax(3),0.04,0.92,'e','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
    xlabel(f2.ax(3),'Time (hr)','fontsize',12);
    ylabel(f2.ax(3),'Med. dist. along propagation (km)','fontsize',12);
    hold(f2.ax(3),'off');
    
    % subplot 4 of figure i
    hold(f2.ax(4),'on');
    scatter(f2.ax(4),medhf(:,3)/3600,medhf(:,2),30,'filled','ro','MarkerEdgeColor', 'k');
    scatter(f2.ax(4),medlf(:,3)/3600,medlf(:,2),30,[0.6 1 1],'filled','o','MarkerEdgeColor', 'k');          
    f2.ax(4).Box = 'on';
    grid(f2.ax(4), 'on');
    f2.ax(4).GridLineStyle = '--';
    xran1 = [trange1(i,2)/3600 trange1(i,3)/3600];
    yran1 = [round(min([medhf(:,2);medlf(:,2)]))-1 ...
             round(max([medhf(:,2);medlf(:,2)]))+1];
    xlim(f2.ax(4),xran1);
    ylim(f2.ax(4),yran1);
    text(f2.ax(4),0.04,0.92,'f','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
    xlabel(f2.ax(4),'Time (hr)','fontsize',12);
    ylabel(f2.ax(4),'Med. dist. along orthogonal (km)','fontsize',12);
    hold(f2.ax(4),'off');
    
    %%% save figure
    print(f2.fig,'-dpdf',strcat(rstpath,'/',fam,'mig.proj.meddist',num2str(trange1(i,1)),...
          '_',num2str(trange1(i,2)),'-',num2str(trange1(i,3)),'.',num2str(winlenhf),'_',...
          num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));
    
end   

%%
% trange2 needs to 2 slopes, visually seperate the time      
trange2 = [2005255,5.80e+4,5.96e+4];
divtime = 58830;
xran = [-8 13];
yran = [-13 10];

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
for iang = 1: length(angle)
    
    mighfdum = mighf;
    for j = 1: size(mighf,1)
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),0,0);
        mighfdum(j,1) = newx;
        mighfdum(j,2) = newy;
    end    
    %%% linear robust least square
    fttpfree = fittype( @(a,b,x) a*x+b);
    % create fit object
    [fitobj,gof,output] = fit(mighfdum(:,14)/3600, mighfdum(:,1),fttpfree,'Robust','Bisquare');
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
    x0 = mighfdum(j,1);
    y0 = mighfdum(j,2);
    [newx,newy] = coordinate_rot(x0,y0,-(angbest(7)-90),0,0);
    mighfdum(j,1) = newx;
    mighfdum(j,2) = newy;
end
miglfdum = miglf;
for j = 1: size(miglf,1)
    x0 = miglfdum(j,1);
    y0 = miglfdum(j,2);
    [newx,newy] = coordinate_rot(x0,y0,-(angbest(7)-90),0,0);
    miglfdum(j,1) = newx;
    miglfdum(j,2) = newy;
end

% break to 2 segments
mighfp1dum = mighfdum(mighfdum(:,14)<divtime,:);
mighfp2dum = mighfdum(mighfdum(:,14)>=divtime,:);
miglfp1dum = miglfdum(miglfdum(:,14)<divtime,:);
miglfp2dum = miglfdum(miglfdum(:,14)>=divtime,:);


%%% define and position the figure frame and axes of each plot
f3.fig=figure;
set(f3.fig,'Position',[scrsz(3)/10 0*scrsz(4)/10 2*scrsz(3)/5 4*scrsz(4)/5]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f3.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
xlen = xran(1,2)-xran(1,1);
ylen = yran(1,2)-yran(1,1);
subfxlen = 0.35;
if xlen >= ylen
    subfylen = ylen/xlen*0.4;
else
    subfylen = ylen/xlen*0.35;
end
set(f3.ax(1), 'position', [ 0.1, 0.5, subfxlen, subfylen]);
set(f3.ax(2), 'position', [ 0.55, 0.5, subfxlen, subfylen]);
set(f3.ax(3), 'position', [ 0.1, 0.2, 0.35, 0.25]);
set(f3.ax(4), 'position', [ 0.55, 0.2, 0.35, 0.25]);

% subplot 1 of figure i
hold(f3.ax(1),'on');
plot(f3.ax(1),[-100 100],[0 0],'k--');
plot(f3.ax(1),[0 0],[-100 100],'k--');
scatter(f3.ax(1),mighf(:,1),mighf(:,2), 30, mighf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
colormap(f3.ax(1),'jet');
c=colorbar(f3.ax(1),'SouthOutside');
pos = f3.ax(1).Position;
c.Position = [pos(1), pos(2), pos(3), 0.02];
% c.TickLabels=[];
c.Label.String = strcat(num2str(trange2(1,1)),' of HF',' (hr)');
c.Label.FontSize = 12;
caxis(f3.ax(1),[trange2(1,2)/3600 trange2(1,3)/3600])
text(f3.ax(1),0.85,0.93,'HF','FontSize',15,'unit','normalized');
text(f3.ax(1),0.04,0.93,'a','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
f3.ax(1).Box = 'on';
grid(f3.ax(1), 'on');
axis(f3.ax(1), 'equal');
f3.ax(1).GridLineStyle = '--';
f3.ax(1).XAxisLocation = 'top';
medx = median(mighf(:,1));
medy = median(mighf(:,2));
[rotx, roty] = complex_rot(0,5,-angbest(7));
xvect = [medx-rotx medx+rotx];
yvect = [medy-roty medy+roty];
drawArrow(f3.ax(1),xvect,yvect,xran,yran,'linewidth',1.5);
[rotx, roty] = complex_rot(-5,0,-angbest(7));
xvect = [medx-rotx medx+rotx];
yvect = [medy-roty medy+roty];
drawArrow(f3.ax(1),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
xlabel(f3.ax(1),'E (km)','fontsize',12);
ylabel(f3.ax(1),'N (km)','fontsize',12);
hold(f3.ax(1),'off');

% subplot 2 of figure i
hold(f3.ax(2),'on');
plot(f3.ax(2),[-100 100],[0 0],'k--');
plot(f3.ax(2),[0 0],[-100 100],'k--');
scatter(f3.ax(2),miglf(:,1),miglf(:,2), 30, miglf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
colormap(f3.ax(2),'jet');
c=colorbar(f3.ax(2),'SouthOutside');
pos = f3.ax(2).Position;
c.Position = [pos(1), pos(2), pos(3), 0.02];
% c.TickLabels=[];
c.Label.String = strcat(num2str(trange2(1,1)),' of LF',' (hr)');
c.Label.FontSize = 12;
caxis(f3.ax(2),[trange2(1,2)/3600 trange2(1,3)/3600])
text(f3.ax(2),0.85,0.93,'LF','FontSize',15,'unit','normalized');
text(f3.ax(2),0.04,0.93,'b','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
f3.ax(2).Box = 'on';
grid(f3.ax(2), 'on');
axis(f3.ax(2), 'equal');
f3.ax(2).GridLineStyle = '--';
f3.ax(2).XAxisLocation = 'top';
medx = median(miglf(:,1));
medy = median(miglf(:,2));
[rotx, roty] = complex_rot(0,5,-angbest(7));
xvect = [medx-rotx medx+rotx];
yvect = [medy-roty medy+roty];
drawArrow(f3.ax(2),xvect,yvect,xran,yran,'linewidth',1.5);
[rotx, roty] = complex_rot(-5,0,-angbest(7));
xvect = [medx-rotx medx+rotx];
yvect = [medy-roty medy+roty];
drawArrow(f3.ax(2),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
xlabel(f3.ax(2),'E (km)','fontsize',12);
ylabel(f3.ax(2),'N (km)','fontsize',12);
hold(f3.ax(2),'off');

% subplot 3 of figure i
hold(f3.ax(3),'on');
scatter(f3.ax(3),mighfp1dum(:,14)/3600,mighfp1dum(:,1),30,'filled','ro','MarkerEdgeColor', 'k');
scatter(f3.ax(3),mighfp2dum(:,14)/3600,mighfp2dum(:,1),30,'filled','ro','MarkerEdgeColor', 'k');
scatter(f3.ax(3),miglfp1dum(:,14)/3600,miglfp1dum(:,1),30,[153/255 255/255 255/255],'filled','o','MarkerEdgeColor', 'k'); 
scatter(f3.ax(3),miglfp2dum(:,14)/3600,miglfp2dum(:,1),30,[153/255 255/255 255/255],'filled','o','MarkerEdgeColor', 'k');
%%% segment 1
% linear robust least square
fttpfree = fittype( @(a,b,x) a*x+b);
% create fit object
[fitobjhfp1,gofhfp1,outphfp1] = fit(mighfp1dum(:,14)/3600, mighfp1dum(:,1),fttpfree,'Robust','Bisquare');
% output fit parameters
coefp1 = coeffvalues(fitobjhfp1);
slopep1hf = coefp1(1);
intcptp1hf = coefp1(2);
fitp1hf = feval(fitobjhfp1,mighfp1dum(:,14)/3600);
plot(f3.ax(3),mighfp1dum(:,14)/3600,fitp1hf,'r-','linewidth',2);

% intcptp1lf = ones(size(miglfp1dum(:,end)))\(miglfp1dum(:,1)-slopep1hf*miglfp1dum(:,end)/3600);
a=slopep1hf;
fttpfix = fittype( @(b,x) a*x+b);
[fitobjlfp1,goflfp1,outplfp1] = fit(miglfp1dum(:,14)/3600, miglfp1dum(:,1),fttpfix,'Robust','Bisquare');
% [fitobjlfp1,goflfp1,outplfp1] = fit(miglfp1dum(:,14)/3600, miglfp1dum(:,1),fttpfix,'Robust',...
%                                     'Bisquare','Weights',miglfp1dum(:,11));
intcptp1lf = coeffvalues(fitobjlfp1);
fitp1lf = feval(fitobjlfp1,miglfp1dum(:,14)/3600);
plot(f3.ax(3),miglfp1dum(:,14)/3600,fitp1lf,'-','linewidth',2,'color',[153/255 255/255 255/255]);

%%% segment 2
% linear robust least square
fttpfree = fittype( @(a,b,x) a*x+b);
% create fit object
[fitobjhfp2,gofhfp2,outphfp2] = fit(mighfp2dum(:,14)/3600, mighfp2dum(:,1),fttpfree,'Robust','Bisquare');
% output fit parameters
coefp2 = coeffvalues(fitobjhfp2);
slopep2hf = coefp2(1);
intcptp2hf = coefp2(2);
fitp2hf = feval(fitobjhfp2,mighfp2dum(:,14)/3600);
plot(f3.ax(3),mighfp2dum(:,14)/3600,fitp2hf,'r-','linewidth',2);
% intcptp1lf = ones(size(miglfp1dum(:,end)))\(miglfp1dum(:,1)-slopep1hf*miglfp1dum(:,end)/3600);
a=slopep2hf;
fttpfix = fittype( @(b,x) a*x+b);
[fitobjlfp2,goflfp2,outplfp2] = fit(miglfp2dum(:,14)/3600, miglfp2dum(:,1),fttpfix,'Robust','Bisquare');
% [fitobjlfp2,goflfp2,outplfp2] = fit(miglfp2dum(:,14)/3600, miglfp2dum(:,1),fttpfix,'Robust',...
%                                     'Bisquare','Weights',miglfp2dum(:,11));
intcptp2lf = coeffvalues(fitobjlfp2);
fitp2lf = feval(fitobjlfp2,miglfp2dum(:,14)/3600);
offset = (intcptp1lf-intcptp1hf+intcptp2lf-intcptp2hf)/2;
plot(f3.ax(3),miglfp2dum(:,14)/3600,fitp2lf,'-','linewidth',2,'color',[153/255 255/255 255/255]);

f3.ax(3).Box = 'on';
grid(f3.ax(3), 'on');
f3.ax(3).GridLineStyle = '--';
xran1 = [trange2(1,2)/3600 trange2(1,3)/3600];
yran1 = [round(min([mighfp1dum(:,1);mighfp2dum(:,1)]))-1 ...
         round(max([mighfp1dum(:,1);mighfp2dum(:,1)]))+1];
xlim(f3.ax(3),xran1);
ylim(f3.ax(3),yran1);
text(f3.ax(3),0.04,0.92,'c','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
xlabel(f3.ax(3),'Time (hr)','fontsize',12);
ylabel(f3.ax(3),'Dist. along propagation (km)','fontsize',12);
hold(f3.ax(3),'off');

% subplot 4 of figure i
hold(f3.ax(4),'on');
numhf(7) = length(mighfdum(:,1));
numlf(7) = length(miglfdum(:,1));
resprophf(7,1:numhf(7)) = [outphfp1.residuals; outphfp2.residuals];
resproplf(7,1:numlf(7)) = [outplfp1.residuals+(intcptp1lf-intcptp1hf); outplfp2.residuals+(intcptp2lf-intcptp2hf)];
pdhf = fitdist(resprophf(7,1:numhf(7))','Normal');    % fit a distribution of residual, assuming it is normal distributed
pdlf = fitdist(resproplf(7,1:numlf(7))','Normal');
muhf = pdhf.mu;    % fitted paramaters
mulf = pdlf.mu;
sigmahf = pdhf.sigma;
sigmalf = pdlf.sigma;
cihf = paramci(pdhf,'Alpha',0.05);     % give the 95% confidence interval of distribution param
cilf = paramci(pdlf,'Alpha',0.05);
pdffithf = pdf(pdhf,min(resprophf(7,1:numhf(7))):0.05:max(resprophf(7,1:numhf(7))));
pdffitlf = pdf(pdlf,min(resproplf(7,1:numlf(7))):0.05:max(resproplf(7,1:numlf(7))));
histogram(f3.ax(4),resprophf(7,1:numhf(7)),'binwidth',0.5,'normalization','pdf','facecolor','r');
histogram(f3.ax(4),resproplf(7,1:numlf(7)),'binwidth',0.5,'normalization','pdf','facecolor',...
          [153/255 255/255 255/255],'facealpha',0.6);
% plot(min(resprophf(7,1:numhf(7))):0.05:max(resprophf(7,1:numhf(7))), pdffithf, 'r-');
% plot(min(resproplf(7,1:numlf(7))):0.05:max(resproplf(7,1:numlf(7))), pdffitlf, '-','color',...
%          [153/255 255/255 255/255]);
text(f3.ax(4),0.04,0.92,'d','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
text(f3.ax(4),0.50,0.92,sprintf('lf-hf = %.2f km',offset),'FontSize',12,'unit','normalized');
f3.ax(4).Box = 'on';
grid(f3.ax(4), 'on');
f3.ax(4).GridLineStyle = '--';
ylim(f3.ax(4),[0 1]);
xlabel(f3.ax(4),'Residual in propagation (km)','fontsize',12);
ylabel(f3.ax(4),'PDF','fontsize',12);
hold(f3.ax(4),'off');

%%% save figure
print(f3.fig,'-dpdf',strcat(rstpath,'/',fam,'mig.proj.lfit',num2str(trange2(1,1)),...
    '_',num2str(trange2(1,2)),'-',num2str(trange2(1,3)),'.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% It seems that median of the previous detections is a promising way to show offsets as well,
%%% so try to deal with it the same way
for j = 1: size(mighf,1)
    medhf(j,1) = median(mighfdum(1:j,1));
    medhf(j,2) = median(mighfdum(1:j,2));
    medhf(j,3) = mighfdum(j,14);
end

for j = 1: size(miglf,1)
    medlf(j,1) = median(miglfdum(1:j,1));
    medlf(j,2) = median(miglfdum(1:j,2));
    medlf(j,3) = miglfdum(j,14);
end

f4.fig=figure;
set(f4.fig,'Position',[0.05*scrsz(3) 0.1*scrsz(4) 0.4*scrsz(3) 0.9*scrsz(4)]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f4.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
xlen = xran(1,2)-xran(1,1);
ylen = yran(1,2)-yran(1,1);
subfxlen = 0.35;
if xlen >= ylen
    subfylen = ylen/xlen*0.4;
else
    subfylen = ylen/xlen*0.35;
end
set(f4.ax(1), 'position', [ 0.1, 0.5, subfxlen, subfylen]);
set(f4.ax(2), 'position', [ 0.55, 0.5, subfxlen, subfylen]);
set(f4.ax(3), 'position', [ 0.1, 0.2, 0.35, 0.25]);
set(f4.ax(4), 'position', [ 0.55, 0.2, 0.35, 0.25]);

% subplot 1 of figure i
hold(f4.ax(1),'on');
plot(f4.ax(1),[-100 100],[0 0],'k--');
plot(f4.ax(1),[0 0],[-100 100],'k--');
scatter(f4.ax(1),mighf(:,1),mighf(:,2), 30, mighf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
colormap(f4.ax(1),'jet');
c=colorbar(f4.ax(1),'SouthOutside');
pos = f4.ax(1).Position;
c.Position = [pos(1), pos(2), pos(3), 0.02];
% c.TickLabels=[];
c.Label.String = strcat(num2str(trange2(1,1)),' of HF',' (hr)');
c.Label.FontSize = 12;
caxis(f4.ax(1),[trange2(1,2)/3600 trange2(1,3)/3600])
text(f4.ax(1),0.85,0.93,'HF','FontSize',15,'unit','normalized');
text(f4.ax(1),0.04,0.93,'a','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
f4.ax(1).Box = 'on';
grid(f4.ax(1), 'on');
axis(f4.ax(1), 'equal');
f4.ax(1).GridLineStyle = '--';
f4.ax(1).XAxisLocation = 'top';
medx = median(mighf(:,1));
medy = median(mighf(:,2));
[rotx, roty] = complex_rot(0,5,-angbest(7));
xvect = [medx-rotx medx+rotx];
yvect = [medy-roty medy+roty];
drawArrow(f4.ax(1),xvect,yvect,xran,yran,'linewidth',1.5);
[rotx, roty] = complex_rot(-5,0,-angbest(7));
xvect = [medx-rotx medx+rotx];
yvect = [medy-roty medy+roty];
drawArrow(f4.ax(1),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
%     xlim(f4.ax(1),xran);
%     ylim(f4.ax(1),yran);
xlabel(f4.ax(1),'E (km)','fontsize',12);
ylabel(f4.ax(1),'N (km)','fontsize',12);
hold(f4.ax(1),'off');

% subplot 2 of figure i
hold(f4.ax(2),'on');
plot(f4.ax(2),[-100 100],[0 0],'k--');
plot(f4.ax(2),[0 0],[-100 100],'k--');
scatter(f4.ax(2),miglf(:,1),miglf(:,2), 30, miglf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
colormap(f4.ax(2),'jet');
c=colorbar(f4.ax(2),'SouthOutside');
pos = f4.ax(2).Position;
c.Position = [pos(1), pos(2), pos(3), 0.02];
% c.TickLabels=[];
c.Label.String = strcat(num2str(trange2(1,1)),' of LF',' (hr)');
c.Label.FontSize = 12;
caxis(f4.ax(2),[trange2(1,2)/3600 trange2(1,3)/3600])
text(f4.ax(2),0.85,0.93,'LF','FontSize',15,'unit','normalized');
text(f4.ax(2),0.04,0.93,'b','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
f4.ax(2).Box = 'on';
grid(f4.ax(2), 'on');
axis(f4.ax(2), 'equal');
f4.ax(2).GridLineStyle = '--';
f4.ax(2).XAxisLocation = 'top';
medx = median(miglf(:,1));
medy = median(miglf(:,2));
[rotx, roty] = complex_rot(0,5,-angbest(7));
xvect = [medx-rotx medx+rotx];
yvect = [medy-roty medy+roty];
drawArrow(f4.ax(2),xvect,yvect,xran,yran,'linewidth',1.5);
[rotx, roty] = complex_rot(-5,0,-angbest(7));
xvect = [medx-rotx medx+rotx];
yvect = [medy-roty medy+roty];
drawArrow(f4.ax(2),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
xlabel(f4.ax(2),'E (km)','fontsize',12);
ylabel(f4.ax(2),'N (km)','fontsize',12);
hold(f4.ax(2),'off');
    

% subplot 3 of figure i
hold(f4.ax(3),'on');
scatter(f4.ax(3),medhf(:,3)/3600,medhf(:,1),30,'filled','ro','MarkerEdgeColor', 'k');
scatter(f4.ax(3),medlf(:,3)/3600,medlf(:,1),30,[153/255 255/255 255/255],'filled','o',...
    'MarkerEdgeColor', 'k');
f4.ax(3).Box = 'on';
grid(f4.ax(3), 'on');
f4.ax(3).GridLineStyle = '--';
xran1 = [trange2(1,2)/3600 trange2(1,3)/3600];
yran1 = [round(min([medhf(:,1);medlf(:,1)]))-1 ...
         round(max([medhf(:,1);medlf(:,1)]))+1];
xlim(f4.ax(3),xran1);
ylim(f4.ax(3),yran1);
text(f4.ax(3),0.04,0.92,'e','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
xlabel(f4.ax(3),'Time (hr)','fontsize',12);
ylabel(f4.ax(3),'Med. dist. along propagation (km)','fontsize',12);
hold(f4.ax(3),'off');

% subplot 4 of figure i
hold(f4.ax(4),'on');
scatter(f4.ax(4),medhf(:,3)/3600,medhf(:,2),30,'filled','ro','MarkerEdgeColor', 'k');
scatter(f4.ax(4),medlf(:,3)/3600,medlf(:,2),30,[153/255 255/255 255/255],'filled','o',...
    'MarkerEdgeColor', 'k');
f4.ax(4).Box = 'on';
grid(f4.ax(4), 'on');
f4.ax(4).GridLineStyle = '--';
xran1 = [trange2(1,2)/3600 trange2(1,3)/3600];
yran1 = [round(min([medhf(:,2);medlf(:,2)]))-1 ...
         round(max([medhf(:,2);medlf(:,2)]))+1];
xlim(f4.ax(4),xran1);
ylim(f4.ax(4),yran1);
text(f4.ax(4),0.04,0.92,'f','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
xlabel(f4.ax(4),'Time (hr)','fontsize',12);
ylabel(f4.ax(4),'Med. dist. along orthogonal (km)','fontsize',12);
hold(f4.ax(4),'off');

%%% save figure
print(f4.fig,'-dpdf',strcat(rstpath,'/',fam,'mig.proj.meddist',num2str(trange2(1,1)),...
    '_',num2str(trange2(1,2)),'-',num2str(trange2(1,3)),'.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));

    
    
%% summary of results
CIFcn = @(x,p)std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * tinv(abs([0,1]-(1-p/100)/2),...
              sum(~isnan(x(:)))-1) + mean(x(:),'omitnan');
f5.fig=figure;
set(f5.fig,'Position',[scrsz(3)/10 2*scrsz(4)/10 2*scrsz(3)/6 2*scrsz(4)/4]);
hold on
plot([0 0],[-100 100],'--','color',[0.5 0.5 0.5]);
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
        e1 = errorbar(offset,ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',8,'color',...
         'k','linewidth',1.2,'MarkerEdgeColor','k','MarkerFaceColor','k');
    elseif i>1 && i~=5
        ypos = i+1;
        errorbar(offset,ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',8,'color',...
         'k','linewidth',1.2,'MarkerEdgeColor','k','MarkerFaceColor','k');
    elseif i==5
        ypos = i+1;
        e3 = errorbar(offset,ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',8,'color',...
         'k','linewidth',1.2,'MarkerEdgeColor','k','MarkerFaceColor',[1 96/255 0]);
    end
    text(offset,ypos+0.3,strcat('HFs: ',num2str(hfN),'; LFs: ',num2str(lfN)),'fontsize',10,...
        'HorizontalAlignment','center');
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
e2 = errorbar(offset,2,lfci95n(1),lfci95n(2),'horizontal','o','markersize',8,'color',...
         'k','linewidth',1.2,'MarkerEdgeColor','k');
text(offset,2+0.3,strcat('HFs: ',num2str(hfN),'; LFs: ',num2str(lfN)),'fontsize',10,...
        'HorizontalAlignment','center');

text(0.02,0.96,'a','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);    
text(0.9,0.95,'PGC','FontSize',14,'unit','normalized','EdgeColor','k','Margin',2);
text(0.2,0.15,'LF','FontSize',14,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(0.2,0.1,'lags','FontSize',14,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(0.8,0.15,'LF','FontSize',14,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(0.8,0.1,'leads','FontSize',14,'unit','normalized','HorizontalAlignment','center','fontw','bold');
xlabel('Offset between HF and LF (LF-HF) (km)','fontsize',12);
ylabel('Migration number from earlier to later','fontsize',12);
legend([e1, e2, e3],{'SE propagation','SE (2 portions)','SW propagation'},'FontSize',11,...
        'Position',[0.185 0.69 0.16 0.08]);
ylim([0,8]);
xlim([-4,4]);
xticks(-4:1:4);
yticks(0:1:8);
box on
grid on
set(gca,'GridLineStyle','--')

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

xran = [ -20 12];
yran = [ -20 12];

resprophf = nan(length(trange1)+1,200);
resproplf = nan(length(trange1)+1,50);
for i = 1: size(trange1,1)
%     i=3;
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
            x0 = mighfdum(j,1);
            y0 = mighfdum(j,2);
            [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),0,0);
            mighfdum(j,1) = newx;
            mighfdum(j,2) = newy;
        end
        
        %%% linear robust least square
        % create fit object
        fttpfree = fittype( @(a,b,x) a*x+b);
        [fitobj,gof,output] = fit(mighfdum(:,14)/3600, mighfdum(:,1),fttpfree,'Robust','Bisquare');
        % output fit parameters
        coef = coeffvalues(fitobj);
        slope(iang) = coef(1);
        fitprophf = feval(fitobj,mighfdum(:,14)/3600);
        l1norm(iang) = sum(abs(mighfdum(:,1)-fitprophf))/(length(mighfdum(:,1)));
        l2norm(iang) = sum((mighfdum(:,1)-fitprophf).^2)/(length(mighfdum(:,1)));
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
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        mighfdum(j,1) = newx;
        mighfdum(j,2) = newy;
    end
    
    indlf = find(lfrelalocall(:,12)==trange1(i,1) & lfrelalocall(:,14)>=trange1(i,2) & ...
                 lfrelalocall(:,14)<=trange1(i,3));
    miglf = lfrelalocall(indlf,:); 
    miglfdum = miglf;
    for j = 1: size(miglf,1)
        x0 = miglfdum(j,1);
        y0 = miglfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        miglfdum(j,1) = newx;
        miglfdum(j,2) = newy;
    end
    
    %%% define and position the figure frame and axes of each plot
    f4.fig = figure;
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 9;   % maximum height allowed is 11 inches
    set(f4.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

    nrow = 3;
    ncol = 2;    
    for isub = 1:nrow*ncol
        f4.ax(isub) = subplot(nrow,ncol,isub);
    end


    %%% reposition
    set(f4.ax(1), 'position', [ 0.1, 0.64, 0.32, 0.32]);
    set(f4.ax(2), 'position', [ 0.5, 0.64, 0.32, 0.32]);
    set(f4.ax(3), 'position', [ 0.1, 0.35, 0.32, 0.22]);
    set(f4.ax(4), 'position', [ 0.5, 0.35, 0.32, 0.22]);
    set(f4.ax(5), 'position', [ 0.1, 0.078, 0.32, 0.22]);
    set(f4.ax(6), 'position', [ 0.5, 0.078, 0.32, 0.22]);
    
    % subplot 1 of figure i
    hold(f4.ax(1),'on');
    plot(f4.ax(1),[-100 100],[0 0],'k--');
    plot(f4.ax(1),[0 0],[-100 100],'k--');
    f4.ax(1).FontSize = 9;
    mighf = sortrows(mighf,-14);
    scatter(f4.ax(1),mighf(:,1),mighf(:,2), 20, mighf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f4.ax(1),'jet');
    c=colorbar(f4.ax(1),'SouthOutside');
    pos = f4.ax(1).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    c.Label.String = strcat({'Time (hr) on 13 Sep. 2005'});
    c.Label.FontSize = 11;
    caxis(f4.ax(1),[trange1(i,2)/3600 trange1(i,3)/3600])
    text(f4.ax(1),0.85,0.15,'HF','FontSize',12,'unit','normalized');
    text(f4.ax(1),0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f4.ax(1).Box = 'on';
    grid(f4.ax(1), 'on');
    axis(f4.ax(1), 'equal');
    f4.ax(1).GridLineStyle = '--';
    f4.ax(1).XAxisLocation = 'top';
    medx = median(mighf(:,1));
    medy = median(mighf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
    drawArrow(f4.ax(1),xvect,yvect,xran(i,:),yran(i,:),'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
    drawArrow(f4.ax(1),xvect,yvect,xran(i,:),yran(i,:),'linewidth',1.5,'linestyle','--','color',...
              [0.4 0.4 0.4]);
    xlabel(f4.ax(1),'E (km)','fontsize',11);
    ylabel(f4.ax(1),'N (km)','fontsize',11);
    text(f4.ax(1),0.5,0.1,'I','FontSize',12,'unit','normalized');
    hold(f4.ax(1),'off');

    % subplot 2 of figure i
    hold(f4.ax(2),'on');
    plot(f4.ax(2),[-100 100],[0 0],'k--');
    plot(f4.ax(2),[0 0],[-100 100],'k--');
    f4.ax(2).FontSize = 9;
    miglf = sortrows(miglf,-14);
    scatter(f4.ax(2),miglf(:,1),miglf(:,2), 20, miglf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f4.ax(2),'jet');
    c=colorbar(f4.ax(2),'SouthOutside');
    pos = f4.ax(2).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    c.Label.String = strcat({'Time (hr) on 13 Sep. 2005'});
    c.Label.FontSize = 11;
    caxis(f4.ax(2),[trange1(i,2)/3600 trange1(i,3)/3600])
    text(f4.ax(2),0.85,0.15,'LF','FontSize',12,'unit','normalized');
    text(f4.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f4.ax(2).Box = 'on';
    grid(f4.ax(2), 'on');
    axis(f4.ax(2), 'equal');
    f4.ax(2).GridLineStyle = '--';
    f4.ax(2).XAxisLocation = 'top';
    medx = median(miglf(:,1));
    medy = median(miglf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
    drawArrow(f4.ax(2),xvect,yvect,xran(i,:),yran(i,:),'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
    drawArrow(f4.ax(2),xvect,yvect,xran(i,:),yran(i,:),'linewidth',1.5,'linestyle','--','color',...
              [0.4 0.4 0.4]);
    xlabel(f4.ax(2),'E (km)','fontsize',11);
    ylabel(f4.ax(2),'N (km)','fontsize',11); 
    text(f4.ax(2),0.5,0.1,'I','FontSize',12,'unit','normalized');
    hold(f4.ax(2),'off');
    
    % subplot 3 of figure i
    hold(f4.ax(3),'on');
    f4.ax(3).FontSize = 9;
    scatter(f4.ax(3),mighfdum(:,14)/3600,mighfdum(:,1),20,'filled','ro','MarkerEdgeColor', 'k');   
    scatter(f4.ax(3),miglfdum(:,14)/3600,miglfdum(:,1),20,[0.6 1 1],'filled','o',...
            'MarkerEdgeColor','k');  %, 'MarkerEdgeColor', 'w')
%     scatter(f.ax(3),mighfdum(:,14)/3600,mighfdum(:,1),30,[105/255 105/255 105/255],'filled','o');
%     scatter(f.ax(3),miglfdum(:,14)/3600,miglfdum(:,1),30,miglfdum(:,11),'filled','o');  
    % linear robust least square
    fttpfree = fittype( @(a,b,x) a*x+b);
    % create fit object
    [fitobjhfprop,gofhfprop,outphfprop] = fit(mighfdum(:,14)/3600, mighfdum(:,1),fttpfree,'Robust',...
                                              'Bisquare');
    % output fit parameters
    coefprop = coeffvalues(fitobjhfprop);
    slopeprophf(i) = coefprop(1);
    intcptprophf(i) = coefprop(2);
    fitprophf = feval(fitobjhfprop,mighfdum(:,14)/3600);
    plot(f4.ax(3),mighfdum(:,14)/3600,fitprophf,'r-','linewidth',2);
        
%     intcptproplf = ones(size(miglfdum(:,end)))\(miglfdum(:,1)-slopeprophf*miglfdum(:,end)/3600); % direct LS fit
%     fitproplf = slopeprophf*miglfdum(:,end)/3600+intcptproplf;
    a=slopeprophf(i);
    fttpfix = fittype( @(b,x) a*x+b);
    [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum(:,14)/3600, miglfdum(:,1),fttpfix,'Robust',...
                                              'Bisquare');
%     [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum(:,14)/3600, miglfdum(:,1),fttpfix,'Robust',...
%                                               'Bisquare','Weights',miglfdum(:,11));
    fitproplf = feval(fitobjlfprop,miglfdum(:,14)/3600);
    intcptproplf(i) = coeffvalues(fitobjlfprop);
    offset(i) = intcptproplf(i)-intcptprophf(i);
    plot(f4.ax(3),miglfdum(:,14)/3600,fitproplf,'-','linewidth',2,'color',[0.6 1 1]);
    f4.ax(3).Box = 'on';
    grid(f4.ax(3), 'on');
    f4.ax(3).GridLineStyle = '--';
    xran1 = [trange1(i,2)/3600 trange1(i,3)/3600];
    yran1 = [round(min([mighfdum(:,1);miglfdum(:,1)]))-1 ...
             round(max([mighfdum(:,1);miglfdum(:,1)]))+1];
    xlim(f4.ax(3),xran1);
    ylim(f4.ax(3),yran1);
    text(f4.ax(3),0.04,0.91,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    xlabel(f4.ax(3),'Time (hr)','fontsize',11);
    ylabel(f4.ax(3),'Dist. along prop. (km)','fontsize',11); 
    hold(f4.ax(3),'off');
    
    % subplot 4 of figure i
    hold(f4.ax(4),'on');
    f4.ax(4).FontSize = 9;
    numhf(i) = length(mighfdum(:,1));
    numlf(i) = length(miglfdum(:,1));
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
    histogram(f4.ax(4),resprophf(i,1:numhf(i)),'binwidth',0.5,'normalization','pdf','facecolor','r');
    histogram(f4.ax(4),resproplf(i,1:numlf(i)),'binwidth',0.5,'normalization','pdf','facecolor',...
              [0.6 1 1],'facealpha',0.6);
%     plot(min(resprophf(i,1:numhf(i))):0.05:max(resprophf(i,1:numhf(i))), pdffithf, 'r-');
%     plot(min(resproplf(i,1:numlf(i))):0.05:max(resproplf(i,1:numlf(i))), pdffitlf, '-','color',...
%          [0.6 1 1]);
    text(f4.ax(4),0.04,0.91,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f4.ax(4),0.4,0.92,sprintf('LF-HF = %.2f km',offset(i)),'FontSize',12,'unit','normalized');
    f4.ax(4).Box = 'on';
    grid(f4.ax(4), 'on');
    f4.ax(4).GridLineStyle = '--';
    ylim(f4.ax(4),[0 0.5]);
    xlabel(f4.ax(4),'Residual in prop. (km)','fontsize',11);
    ylabel(f4.ax(4),'PDF estimate','fontsize',11);
    hold(f4.ax(4),'off');
    

    %%% It seems that median of the previous detections is a promising way to show offsets as well,
    %%% so try to deal with it the same way
    for j = 1: size(mighf,1)
        medhf(j,1) = median(mighfdum(1:j,1));
        medhf(j,2) = median(mighfdum(1:j,2));
        medhf(j,3) = mighfdum(j,14);
    end
    
    for j = 1: size(miglf,1)
        medlf(j,1) = median(miglfdum(1:j,1));
        medlf(j,2) = median(miglfdum(1:j,2));
        medlf(j,3) = miglfdum(j,14);
    end
   
    
    % subplot 5 of figure i
    hold(f4.ax(5),'on');
    f4.ax(5).FontSize = 9;
    scatter(f4.ax(5),medhf(:,3)/3600,medhf(:,1),20,'filled','ro','MarkerEdgeColor', 'k');
    scatter(f4.ax(5),medlf(:,3)/3600,medlf(:,1),20,[0.6 1 1],'filled','o','MarkerEdgeColor', 'k');            
    f4.ax(5).Box = 'on';
    grid(f4.ax(5), 'on');
    f4.ax(5).GridLineStyle = '--';
    xran1 = [trange1(i,2)/3600 trange1(i,3)/3600];
    yran1 = [round(min([medhf(:,1);medlf(:,1)]))-1 ...
             round(max([medhf(:,1);medlf(:,1)]))+1];
    xlim(f4.ax(5),xran1);
    ylim(f4.ax(5),yran1);
    text(f4.ax(5),0.97,0.1,'e','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2,...
         'HorizontalAlignment','right');
    xlabel(f4.ax(5),'Time (hr)','fontsize',11);
    ylabel(f4.ax(5),'Med. dist. along prop. (km)','fontsize',11);
    hold(f4.ax(5),'off');
    
    % subplot 6 of figure i
    hold(f4.ax(6),'on');
    f4.ax(6).FontSize = 9;
    scatter(f4.ax(6),medhf(:,3)/3600,medhf(:,2),20,'filled','ro','MarkerEdgeColor', 'k');
    scatter(f4.ax(6),medlf(:,3)/3600,medlf(:,2),20,[0.6 1 1],'filled','o','MarkerEdgeColor', 'k');          
    f4.ax(6).Box = 'on';
    grid(f4.ax(6), 'on');
    f4.ax(6).GridLineStyle = '--';
    xran1 = [trange1(i,2)/3600 trange1(i,3)/3600];
    yran1 = [round(min([medhf(:,2);medlf(:,2)]))-1 ...
             round(max([medhf(:,2);medlf(:,2)]))+1];
    xlim(f4.ax(6),xran1);
    ylim(f4.ax(6),yran1);
    text(f4.ax(6),0.97,0.1,'f','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2,...
         'HorizontalAlignment','right');
    xlabel(f4.ax(6),'Time (hr)','fontsize',11);
    ylabel(f4.ax(6),'Med. dist. along ort. (km)','fontsize',11);
    hold(f4.ax(6),'off');
    
    %%% save figure
    print(f4.fig,'-dpdf',strcat(rstpath,'/',fam,'sel.mig.proj.meddist',num2str(trange1(i,1)),...
          '_',num2str(trange1(i,2)),'-',num2str(trange1(i,3)),'.',num2str(winlenhf),'_',...
          num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));
end    


% %% try to tackle the migration with the distint direction
% trangex = [2005256,0.36e+4,0.55e+4;];
% 
% xran = [-21 13];
% yran = [-23 9];
% 
% resprophfx = nan(length(trangex)+1,200);
% resproplfx = nan(length(trangex)+1,50);
% 
% angbestx = 125;
% 
% indhf = find(hfrelalocall(:,8)==trangex(1) & hfrelalocall(:,10)>=trangex(2) & ...
%     hfrelalocall(:,10)<=trangex(3));
% mighf = hfrelalocall(indhf,:);
% mighfdum = mighf;
% for j = 1: size(mighf,1)
%     x0 = mighfdum(j,1);
%     y0 = mighfdum(j,2);
%     [newx,newy] = coordinate_rot(x0,y0,-(angbestx-90),0,0);
%     mighfdum(j,1) = newx;
%     mighfdum(j,2) = newy;
% end
% 
% indlf = find(lfrelalocall(:,8)==trangex(1) & lfrelalocall(:,10)>=trangex(2) & ...
%     lfrelalocall(:,10)<=trangex(3));
% miglf = lfrelalocall(indlf,:);
% miglfdum = miglf;
% for j = 1: size(miglf,1)
%     x0 = miglfdum(j,1);
%     y0 = miglfdum(j,2);
%     [newx,newy] = coordinate_rot(x0,y0,-(angbestx-90),0,0);
%     miglfdum(j,1) = newx;
%     miglfdum(j,2) = newy;
% end
% 
% for j = 1: size(mighf,1)
%     medhf(j,1) = median(mighfdum(1:j,1));
%     medhf(j,2) = median(mighfdum(1:j,2));
%     medhf(j,3) = mighfdum(j,10);        
% end
% 
% for j = 1: size(miglf,1)
%     medlf(j,1) = median(miglfdum(1:j,1));
%     medlf(j,2) = median(miglfdum(1:j,2));
%     medlf(j,3) = miglfdum(j,10);        
% end
% 
% %%% define and position the figure frame and axes of each plot
% f4.fig=figure;
% set(f4.fig,'Position',[scrsz(3)/10 0*scrsz(4)/10 2*scrsz(3)/5 4*scrsz(4)/5]);
% nrow = 2;
% ncol = 2;
% for isub = 1:nrow*ncol
%     f4.ax(isub) = subplot(nrow,ncol,isub);
% end
% 
% %%% reposition
% xlen = xran(2)-xran(1);
% ylen = yran(2)-yran(1);
% subfxlen = 0.35;
% if xlen >= ylen
%     subfylen = ylen/xlen*0.4;
% else
%     subfylen = ylen/xlen*0.35;
% end
% set(f4.ax(1), 'position', [ 0.1, 0.5, subfxlen, subfylen]);
% set(f4.ax(2), 'position', [ 0.55, 0.5, subfxlen, subfylen]);
% set(f4.ax(3), 'position', [ 0.1, 0.2, 0.35, 0.25]);
% set(f4.ax(4), 'position', [ 0.55, 0.2, 0.35, 0.25]);
% 
% % subplot 1 of figure i
% hold(f4.ax(1),'on');
% scatter(f4.ax(1),mighf(:,1),mighf(:,2), 30, mighf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
% colormap(f4.ax(1),'jet');
% c=colorbar(f4.ax(1),'SouthOutside');
% pos = f4.ax(1).Position;
% c.Position = [pos(1), pos(2), pos(3), 0.02];
% %     c.TickLabels=[];
% c.Label.String = strcat(num2str(trangex(1)),' of HF',' (hr)');
% c.Label.FontSize = 12;
% caxis(f4.ax(1),[trangex(2)/3600 trangex(3)/3600])
% plot(f4.ax(1),[-100 100],[0 0],'k--');
% plot(f4.ax(1),[0 0],[-100 100],'k--');
% text(f4.ax(1),(xran(1)+xran(2))/2,yran(2)-5,'HF','FontSize',15);
% text(f4.ax(1),xran(1)+2,yran(2)-2,'(a)','FontSize',13);
% f4.ax(1).Box = 'on';
% grid(f4.ax(1), 'on');
% axis(f4.ax(1), 'equal');
% f4.ax(1).GridLineStyle = '--';
% f4.ax(1).XAxisLocation = 'top';
% [rotx, roty] = complex_rot(0,5,-angbestx);
% xvect = [-rotx rotx];
% yvect = [-roty roty];
% drawArrow(f4.ax(1),xvect,yvect,xran,yran,'linewidth',1.5);
% xlabel(f4.ax(1),'E (km)','fontsize',12);
% ylabel(f4.ax(1),'N (km)','fontsize',12);
% hold(f4.ax(1),'off');
% 
% % subplot 2 of figure i
% hold(f4.ax(2),'on');
% scatter(f4.ax(2),miglf(:,1),miglf(:,2), 30, miglf(:,14)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
% colormap(f4.ax(2),'jet');
% c=colorbar(f4.ax(2),'SouthOutside');
% pos = f4.ax(2).Position;
% c.Position = [pos(1), pos(2), pos(3), 0.02];
% %     c.TickLabels=[];
% c.Label.String = strcat(num2str(trangex(1)),' of LF',' (hr)');
% c.Label.FontSize = 12;
% caxis(f4.ax(2),[trangex(2)/3600 trangex(3)/3600])
% plot(f4.ax(2),[-100 100],[0 0],'k--');
% plot(f4.ax(2),[0 0],[-100 100],'k--');
% text(f4.ax(2),(xran(1)+xran(2))/2,yran(2)-5,'LF','FontSize',15);
% text(f4.ax(2),xran(1)+2,yran(2)-2,'(b)','FontSize',13);
% f4.ax(2).Box = 'on';
% grid(f4.ax(2), 'on');
% axis(f4.ax(2), 'equal');
% f4.ax(2).GridLineStyle = '--';
% f4.ax(2).XAxisLocation = 'top';
% xlim(f4.ax(2),xran);
% ylim(f4.ax(2),yran);
% xlabel(f4.ax(2),'E (km)','fontsize',12);
% ylabel(f4.ax(2),'N (km)','fontsize',12);
% hold(f4.ax(2),'off');
% 
% % subplot 3 of figure i
% hold(f4.ax(3),'on');
% scatter(f4.ax(3),medhf(:,3)/3600,medhf(:,1),30,'filled','ro','MarkerEdgeColor', 'k');
% scatter(f4.ax(3),medlf(:,3)/3600,medlf(:,1),30,[153/255 255/255 255/255],'filled','o','MarkerEdgeColor', 'k');
% f4.ax(3).Box = 'on';
% grid(f4.ax(3), 'on');
% f4.ax(3).GridLineStyle = '--';
% xran1 = [trangex(2)/3600 trangex(3)/3600];
% yran1 = [round(min(medhf(:,1)))-1 round(max(medhf(:,1)))+1];
% xlim(f4.ax(3),xran1);
% ylim(f4.ax(3),yran1);
% text(f4.ax(3),xran1(1)+0.04,yran1(2)-2,'(c)','FontSize',13);
% xlabel(f4.ax(3),'Time (hr)','fontsize',12);
% ylabel(f4.ax(3),'Median distance along propagation (km)','fontsize',12);
% hold(f4.ax(3),'off');
% 
% % subplot 4 of figure i
% hold(f4.ax(4),'on');
% scatter(f4.ax(4),medhf(:,3)/3600,medhf(:,2),30,'filled','ro','MarkerEdgeColor', 'k');
% scatter(f4.ax(4),medlf(:,3)/3600,medlf(:,2),30,[153/255 255/255 255/255],'filled','o','MarkerEdgeColor', 'k');
% f4.ax(4).Box = 'on';
% grid(f4.ax(4), 'on');
% f4.ax(4).GridLineStyle = '--';
% xran1 = [trangex(2)/3600 trangex(3)/3600];
% yran1 = [round(min(medhf(:,2)))-1 round(max(medhf(:,2)))+1];
% xlim(f4.ax(4),xran1);
% ylim(f4.ax(4),yran1);
% text(f4.ax(4),xran1(1)+0.04,yran1(2)-2,'(d)','FontSize',13);
% xlabel(f4.ax(4),'Time (hr)','fontsize',12);
% ylabel(f4.ax(4),'Median distance along orthogonal (km)','fontsize',12);
% hold(f4.ax(4),'off');


% % subplot 3 of figure i
% hold(f4.ax(3),'on');
% scatter(f4.ax(3),mighfdum(:,14)/3600,mighfdum(:,1),30,'filled','ro','MarkerEdgeColor', 'k');
% scatter(f4.ax(3),miglfdum(:,14)/3600,miglfdum(:,1),30,[153/255 255/255 255/255],'filled','o','MarkerEdgeColor', 'k');  %, 'MarkerEdgeColor', 'w')
% %     scatter(f4.ax(3),mighfdum(:,14)/3600,mighfdum(:,1),30,[105/255 105/255 105/255],'filled','o');  %, 'MarkerEdgeColor', 'w')
% %     scatter(f4.ax(3),miglfdum(:,14)/3600,miglfdum(:,1),30,miglfdum(:,11),'filled','o');  %, 'MarkerEdgeColor', 'w')
% % linear robust least square
% fttpfree = fittype( @(a,b,x) a*x+b);
% % create fit object
% [fitobjhfprop,gofhfprop,outphfprop] = fit(mighfdum(:,14)/3600, mighfdum(:,1),fttpfree,'Robust','Bisquare');
% % output fit parameters
% coefprop = coeffvalues(fitobjhfprop);
% slopeprophfx = coefprop(1);
% intcptprophfx = coefprop(2);
% fitprophf = feval(fitobjhfprop,mighfdum(:,14)/3600);
% plot(f4.ax(3),mighfdum(:,14)/3600,fitprophf,'r-','linewidth',2);
% 
% %     intcptproplf = ones(size(miglfdum(:,end)))\(miglfdum(:,1)-slopeprophf*miglfdum(:,end)/3600); % direct LS fit
% %     fitproplf = slopeprophf*miglfdum(:,end)/3600+intcptproplf;
% a=slopeprophfx;
% fttpfix = fittype( @(b,x) a*x+b);
% [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum(:,14)/3600, miglfdum(:,1),fttpfix,'Robust','Bisquare');
% %     [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum(:,14)/3600, miglfdum(:,1),fttpfix,'Robust',...
% %                                               'Bisquare','Weights',miglfdum(:,11));
% fitproplf = feval(fitobjlfprop,miglfdum(:,14)/3600);
% intcptproplfx = coeffvalues(fitobjlfprop);
% plot(f4.ax(3),miglfdum(:,14)/3600,fitproplf,'-','linewidth',2,'color',[153/255 255/255 255/255]);
% f4.ax(3).Box = 'on';
% grid(f4.ax(3), 'on');
% f4.ax(3).GridLineStyle = '--';
% xran1 = [trangex(2)/3600 trangex(3)/3600];
% yran1 = [round(min(mighfdum(:,1)))-1 round(max(mighfdum(:,1)))+1];
% xlim(f4.ax(3),xran1);
% ylim(f4.ax(3),yran1);
% text(f4.ax(3),xran1(1)+0.04,yran1(2)-2,'(c)','FontSize',13);
% xlabel(f4.ax(3),'Time (hr)','fontsize',12);
% ylabel(f4.ax(3),'Distance along propagation (km)','fontsize',12);
% hold(f4.ax(3),'off');
% 
% % subplot 4 of figure i
% hold(f4.ax(4),'on');
% numhf = length(mighfdum(:,1));
% numlf = length(miglfdum(:,1));
% resprophf(1:numhf) = outphfprop.residuals;
% resproplf(1:numlf) = outplfprop.residuals+(intcptproplfx-intcptprophfx);
% pdhf = fitdist(resprophf(1:numhf)','Normal');    % fit a distribution of residual, assuming it is normal distributed
% pdlf = fitdist(resproplf(1:numlf)','Normal');
% muhf = pdhf.mu;    % fitted paramaters
% mulf = pdlf.mu;
% sigmahf = pdhf.sigma;
% sigmalf = pdlf.sigma;
% cihf = paramci(pdhf,'Alpha',0.05);     % give the 95% confidence interval of distribution param
% cilf = paramci(pdlf,'Alpha',0.05);
% pdffithf = pdf(pdhf,min(resprophf(1:numhf)):0.05:max(resprophf(1:numhf)));
% pdffitlf = pdf(pdlf,min(resproplf(1:numlf)):0.05:max(resproplf(1:numlf)));
% histogram(f4.ax(4),resprophf(1:numhf),'binwidth',0.5,'normalization','pdf','facecolor','r');
% histogram(f4.ax(4),resproplf(1:numlf),'binwidth',0.5,'normalization','pdf','facecolor',[153/255 255/255 255/255],'facealpha',0.6);
% plot(min(resprophf(1:numhf)):0.05:max(resprophf(1:numhf)), pdffithf, 'r-');
% plot(min(resproplf(1:numlf)):0.05:max(resproplf(1:numlf)), pdffitlf, '-','color',[153/255 255/255 255/255]);
% text(f4.ax(4),min([resproplf(1:numlf) resprophf(1:numhf)]),0.9,'(d)','FontSize',13);
% f4.ax(4).Box = 'on';
% grid(f4.ax(4), 'on');
% f4.ax(4).GridLineStyle = '--';
% ylim(f4.ax(4),[0 1]);
% xlabel(f4.ax(4),'Residual in propagation (km)','fontsize',12);
% ylabel(f4.ax(4),'Normalization by PDF estimate','fontsize',12);
% hold(f4.ax(4),'off');
    
%     %%% save figure
%     print(f4.fig,'-depsc2',strcat(rstpath,'/',fam,'proj_statistics.contemp_migration.',num2str(trange1(i,1)),...
%           '_',num2str(trange1(i,2)),'-',num2str(trange1(i,3)),'.',num2str(winlenhf),'_',...
%           num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.timerela.eps'));
% end   
















