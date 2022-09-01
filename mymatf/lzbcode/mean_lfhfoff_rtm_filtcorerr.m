% function mean_lfhfoff_rtm_filtcorerr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to summarize every piece of information and result that
% are related to estimate the offset between the lf and hf inside each rtm 
% groups considering the measurement error and sysmatic error in the filtering
% effect correction. This leads to the figure D1 and D2 in the appendix D in 
% the paper Song&Rubin2021JGR
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/08/25
% Last modified date:   2021/09/08
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

%%% load the location resolution points, +-1 & 2 samples
reslzb1 = load(fullfile(rstpath,'evtloc.13fam_locres_hf1spl'));
reslzb2 = load(fullfile(rstpath,'evtloc.13fam_locres_hf2spl'));

% sort according to day, sec
hftime = sortrows(hftime, [13, 15]);
lftime = sortrows(lftime, [13, 15]);

% this is new time ranges that come from 'identify_RTMs_v3.m'
trange = [
    2004195   4.1773e+04   4.3362e+04   % combine, y
    2004197   4.0974e+04   4.2205e+04   % speed direct 90, y, larger pear
    2004197   5.4453e+04   5.5867e+04   % speed direct 205, y, but use rmse 220
    2004197   5.6270e+04   5.7227e+04   % speed direct 175, y, but use rmse 185
%     2004197   7.9949e+04   8.1453e+04   % speed direct 260, y, but use rmse 195
    2004197   8.3402e+04   8.4888e+04   % divided into 3, 1 discarded, y
    2004197   8.5320e+04   8.6659e+04   % divided into 3, 1 discarded, y
    2004198   5.9530e+03   9.7940e+03   % speed direct 270, y, but use rmse 235
    2004198   1.9151e+04   2.1616e+04   % acceptted, y
    2004198   4.2632e+04   4.3795e+04   % acceptted, y
    2004198   5.3506e+04   5.8050e+04   % acceptted, y
    2004198   6.2966e+04   6.3885e+04   % acceptted, y
    2004198   7.3540e+04   7.6029e+04   % acceptted, y 
    2004198   8.4222e+04   8.6400e+04   % combine, y
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
%     2005255   3.4400e+04   3.5612e+04   % speed direct 80, y, larger pear
    2005255   6.0381e+04   6.1557e+04   % speed direct 150, y, but use rmse 110
    2005255   6.7268e+04   6.8535e+04   % speed direct 80, y, larger pear
    2005255   8.4296e+04   8.5634e+04   % combine them, y
%     2005256   3.6200e+03   4.9570e+03   % speed direct 265, y, larger pear
    2005256   8.0330e+03   9.8280e+03   % change end time, y
    2005256   1.3412e+04   1.4346e+04   % speed direct 65, y, but use rmse 80   
%     2005256   1.9069e+04   2.0203e+04   % speed direct 195, y, larger pear
    2005256   2.1492e+04   2.2294e+04   % change start time, y
    2005256   2.8080e+04   2.9172e+04   % change start time, y
    2005256   3.4776e+04   3.5238e+04   % change start time, y
    2005256   3.6900e+04   3.8520e+04   % change times, y
%     2005256   4.6440e+04   4.8421e+04   % change start time, y, check speed direc, y, a bit awful
    2005256   5.2020e+04   5.2596e+04   % change times, y
    2005256   6.2640e+04   6.4865e+04   % change start time, y, check speed direc, y, but use rmse 130
    2005256   7.2720e+04   7.3914e+04   % divided into 2, y, but check speed direc, y
    2005256   7.6381e+04   7.7940e+04   % change end time, y
    2005257   6.1200e+03   7.7460e+03   % change start time, y, but check speed direc, y
%     2005257   1.2240e+04   1.3248e+04   % CHANGED time, st-3.68, y, use speed direc, larger pear, SUSPICOUS
    2005257   1.6920e+04   1.8655e+04   % change start time, y
    2005257   2.1070e+04   2.5365e+04   % accepted, y 
    2005257   3.9000e+04   4.5418e+04   % change start time, y
    2005257   6.1449e+04   6.4012e+04   % accepted, y
    2005257   7.3440e+04   7.8840e+04   % change end time, y, but check speed direc, y, but use rmse 235
    2005257   8.2159e+04   8.3380e+04   % speed direct 350, y, larger pear
    2005258   3.5000e+04   3.6540e+04   % divided into 3, y, but check speed direc, RECHECK PEAR, slightly larger pear, SUSPICOUS
    2005258   3.6540e+04   3.8160e+04   % divided into 3, y
    2005258   3.8160e+04   4.0021e+04   % divided into 3, y, but check speed direc y, slightly larger pear 
%     2005259   1.5400e+02   1.6340e+03   % speed direct 75, y, larger pear
%     2005259   1.9560e+03   4.0460e+03   % speed direct 185, y, but use rmse 140
    2005259   5.2320e+03   7.7140e+03   % check OLD param. to determine, y, but check speed direc, larger pear, use old direc 255
    2005259   3.7000e+04   4.1293e+04   % speed direct 85, RECHECK PEAR, y, but use rmse 115
%     2005260   2.6510e+03   3.4380e+03   % speed direct 80, y, but use rmse 115
    2005260   3.6040e+03   4.5000e+03   % divided into 2, y
%     2005260   4.6080e+03   6.6780e+03   % divided into 2, y, but check speed direc, y, larger pear
    2005260   5.6300e+04   5.8200e+04   % use old times, y
%     2005260   9.4090e+03   1.0357e+04   % speed direct 55, y, but use rmse 115
    ]; 
    
% hardcopy of angbest from either min. rmse or max. slope of accepted ones ntranacc
angbest = [
   250
    90
   220
   185
   225
   270
   235
   200
   230
    80
    65
   235
   245
   310
   250
   120
   195
    90
   140
   130
    70
   255
   250
   260
   120
   110
    80
   240
   220
    80
   115
   105
   260
    95
   260
   130
   230
   250
   195
   115
   250
   165
   245
   235
    10
    95
    85
   185
   255
   115
   250
   235
   ];

indspeed = [
    2
    27
    34
    37
    47
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

xran = [-20 25];
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
%     i=46;
    disp(trange(i,:));
    indhf = find(hftime(:,13)==trange(i,1) & hftime(:,15)>=trange(i,2) & ...
                 hftime(:,15)<=trange(i,3));
    mighf = hftime(indhf,:);
    
    indlf = find(lftime(:,13)==trange(i,1) & lftime(:,15)>=trange(i,2) & ...
                 lftime(:,15)<=trange(i,3));
    miglf = lftime(indlf,:); 
    
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
    
    % detections involved into fitting 
    mintime = max(min(mighfdum(:,15)), min(miglfdum(:,15)));
    maxtime = min(max(mighfdum(:,15)), max(miglfdum(:,15)));
    mighfdum2 = mighfdum(mighfdum(:,15)>=mintime & mighfdum(:,15)<=maxtime, :);
    miglfdum2 = miglfdum(miglfdum(:,15)>=mintime & miglfdum(:,15)<=maxtime, :);

    mighf = sortrows(mighf,-15);
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
    medxhf = median(mighf(:,1));
    medyhf = median(mighf(:,2));

    miglf = sortrows(miglf,-15);
    medxlf = median(miglf(:,1));
    medylf = median(miglf(:,2));
    
    [fitobjhfprop,gofhfprop,outphfprop] = fit(mighfdum2(:,15)/3600, mighfdum2(:,1),fttpfree,'Robust',...
                                              'Bisquare','StartPoint',[1 1]);
    % output fit parameters
    coefprop = coeffvalues(fitobjhfprop);
    slopeprophf(i) = coefprop(1);
    intcptprophf(i) = coefprop(2);
    fitprophf = feval(fitobjhfprop,mighfdum2(:,15)/3600);
        
    a=slopeprophf(i);
    fttpfix = fittype( @(b,x) a*x+b);
    [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum2(:,15)/3600, miglfdum2(:,1),fttpfix,'Robust',...
                                              'Bisquare','StartPoint',1);
    fitproplf = feval(fitobjlfprop,miglfdum2(:,15)/3600);
    intcptproplf(i) = coeffvalues(fitobjlfprop);
    offset(i) = intcptproplf(i)-intcptprophf(i);
    
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
      
    %%% find the independent LF detections, i.e. non-overlapping time window
    [miglfindpen,indindpen] = find_independent_detection(miglfdum2, 0);
    inumlf(i) = length(indindpen);
    indlfindpen(1: inumlf(i), i) =  indindpen;
    
    
        
    numhf(i) = length(mighfdum2(:,1));
    numlf(i) = length(miglfdum2(:,1));
    resprophf(i,1:numhf(i)) = outphfprop.residuals;
    resproplf(i,1:numlf(i)) = outplfprop.residuals+(intcptproplf(i)-intcptprophf(i));

    % returns estimates of normal distribution parameters, mu and sigma
    % under 95% confidence interval
    [muhf,sigmahf,muhfCI,sigmahfCI] = normfit((resprophf(i,1:numhf(i)))', 0.05, [], wthf);
    [mulf,sigmalf,mulfCI,sigmalfCI] = normfit((resproplf(i,1:numlf(i)))', 0.05, [], wtlf);

end

[tmp,tmpind] = sort(sehf(:,2),'descend');
sortsehf = [tmpind tmp];
[tmp,tmpind] = sort(pearwthf(:),'ascend');
sortpearwthf = [tmpind tmp];
[tmp,tmpind] = sort(self(:,2),'descend');
sortself = [tmpind tmp];

% close all

%% obtain the vertical distance from lf points to hf line in certain regions

% NOTE here, 2020/09/22, u can choose which vertical distance to use
% vertdist = vertdistwt;
vertdist = vertdistraw;

% SW migs
% ones within angle range and excluding those at distant fams, eg. 002, 010
indsw = setdiff(find(angbest>=220 & angbest<=280),[1]); % [3,5,6,7,9,12,13,15,22,23,24,28,29,33,35,37,38,41,43,44,49,51,52]
% ones have poor lf constraints
indpoorlf = sortself(sortself(:,2)>=3, 1); %[40,4,5,33,3,15,11,23,35,21,45,17,2,31,30,24,32,39,6,9,27,37,51,19,48,26,46,22,8,29,18,47]
% ones have large difference in offsets due to a different selection of propagation direction
indlargediff = [2,3,11,14,24,29,32,34,36,40,43,44,49,52];  %
% so you can have several options, comment others:  USE OPTION 4!
% % 1. all within angles
% indsw = indsw;  % [3,5,6,7,9,12,13,15,22,23,24,28,29,33,35,37,38,41,43,44,49,51,52]
% % 2. all excluding EITHER poor LF constraints OR there will be very different offsets
% tmp = union(indpoorlf, indlargediff);
% indsw = setdiff(indsw, tmp); % [7,12,13,28,38,41]
% % 3. all excluding ones have large difference in offsets
% indsw = setdiff(indsw, indlargediff);  % [5,6,7,9,12,13,15,22,23,28,33,35,37,38,41,51]
% 4. all excluding ones have large difference in offsets due to LF constraint
tmp = intersect(indlargediff, indpoorlf);
indsw = setdiff(indsw, tmp);   % [5,6,7,9,12,13,15,22,23,28,33,35,37,38,41,43,44,49,51,52]

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

        temp1 = vertdist(1:numlf(indsw(i)), indsw(i));
%         temp1 = vertdist(1:numlf(indsw(i)), indsw(i))*cos(deg2rad(angbest(indsw(i)) - 250));
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
end


% NE migs
% ones within angle range and excluding those at distant fams, eg. 002, 010
indne = setdiff(find(angbest>=40 & angbest<=100),[27]); % [2,10,11,18,21,30,34,46,47]
% ones have poor lf constraints
indpoorlf = sortself(sortself(:,2)>=3, 1); %[40,4,5,33,3,15,11,23,35,21,45,17,2,31,30,24,32,39,6,9,27,37,51,19,48,26,46,22,8,29,18,47]
% ones have large difference in offsets due to a different selection of propagation direction
indlargediff = [2,3,11,14,24,29,32,34,36,40,43,44,49,52];  %
% so you can have several options, comment out others:  USE OPTION 4!
% % 1. all within angles
% indne = indne;  % [2,10,11,18,21,30,34,46,47]
% % 2. all excluding EITHER poor LF constraints OR there will be very different offsets
% tmp = union(indpoorlf, indlargediff);
% indne = setdiff(indne, tmp); % [10]
% % 3. all excluding ones have large difference in offsets
% indne = setdiff(indne, indlargediff);  % [10,18,21,30,46,47]
% 4. all excluding ones have large difference in offsets due to LF constraint
tmp = intersect(indlargediff, indpoorlf);
indne = setdiff(indne, tmp);   % [10,18,21,30,34,46,47]

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
%     temp1 = vertdist(1:numlf(indne(i)), indne(i))*cos(deg2rad(angbest(indne(i)) - 70));
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

%%% if we are even more strict, only use ones whose diff in offset is smaller than 0.5 km from the
%%% previous NE and SW migration sets
indsmalldiff = [1,5,6,8,9,10,13,15,16,17,18,19,20,21,22,23,25,31,33,38,41,42,45,46,47,48,50,51];
indne2 = intersect(indne, indsmalldiff);
indsw2 = intersect(indsw, indsmalldiff);

vdistsw2 = [];
wtsw2 = [];
ivdistsw2 = [];  % dist and weight of independent detections
iwtsw2 = [];
for i = 1: length(indsw2)
    iind = indlfindpen(1:inumlf(indsw2(i)), indsw2(i));  % index of independent detection in that migration indsw(i)
        temp1 = vertdist(1:numlf(indsw2(i)), indsw2(i));
%         temp1 = vertdist(1:numlf(indsw(i)), indsw(i))*cos(deg2rad(angbest(indsw(i)) - 250));
        vdistsw2 = [vdistsw2; temp1];
        ivdistsw2 = [ivdistsw2; temp1(iind)];
        
        temp1 = weight(1:numlf(indsw2(i)), indsw2(i));
        wtsw2 = [wtsw2; temp1];
        iwtsw2 = [iwtsw2; temp1(iind)];
end

vdistne2 = [];
wtne2 = [];
ivdistne2 = [];  % dist and weight of independent detections
iwtne2 = [];
for i = 1: length(indne2)
    iind = indlfindpen(1:inumlf(indne2(i)), indne2(i));  % index of independent detection in that migration indsw(i)
    temp1 = vertdist(1:numlf(indne2(i)), indne2(i));
%     temp1 = vertdist(1:numlf(indne(i)), indne(i))*cos(deg2rad(angbest(indne(i)) - 70));
    vdistne2 = [vdistne2; temp1];
    ivdistne2 = [ivdistne2; temp1(iind)];
    
    temp1 = weight(1:numlf(indne2(i)), indne2(i));
    wtne2 = [wtne2; temp1];
    iwtne2 = [iwtne2; temp1(iind)];

end

%%% combine the other unused migrations
indlzb = [indne; indsw];
% index of all other RTMs
indall = 1: 1: size(trange,1);
indelse = setdiff(indall, indlzb);

%%% EXCLUDE the 3 migrations at 002, 010, etc.
indexc = [1; 26; 27];
% indexc = [1; 27];
% indexc = [1; 26];
% indexc = [27];
indelse = setdiff(indelse, indexc);

vdistelse = [];
wtelse = [];
ivdistelse = [];  % dist and weight of independent detections
iwtelse = [];
for i = 1: length(indelse)
    iind = indlfindpen(1:inumlf(indelse(i)), indelse(i));  % index of independent detection in that migration indsw(i)
    temp1 = vertdist(1:numlf(indelse(i)), indelse(i));
%     temp1 = vertdist(1:numlf(indne(i)), indne(i))*cos(deg2rad(angbest(indne(i)) - 70));
    vdistelse = [vdistelse; temp1];
    ivdistelse = [ivdistelse; temp1(iind)];
    
    temp1 = weight(1:numlf(indelse(i)), indelse(i));
    wtelse = [wtelse; temp1];
    iwtelse = [iwtelse; temp1(iind)];

end

%%% in this particular case, we combine the oppoing migrations into one group
vdistlzb = [];
wtlzb = [];
ivdistlzb = [];  % dist and weight of independent detections
iwtlzb = [];
for i = 1: length(indlzb)
    iind = indlfindpen(1:inumlf(indlzb(i)), indlzb(i));  % index of independent detection in that migration indsw(i)
    temp1 = vertdist(1:numlf(indlzb(i)), indlzb(i));
%     temp1 = vertdist(1:numlf(indne(i)), indne(i))*cos(deg2rad(angbest(indne(i)) - 70));
    vdistlzb = [vdistlzb; temp1];
    ivdistlzb = [ivdistlzb; temp1(iind)];
    
    temp1 = weight(1:numlf(indlzb(i)), indlzb(i));
    wtlzb = [wtlzb; temp1];
    iwtlzb = [iwtlzb; temp1(iind)];

end

%% rose diagram for other rtms using different weighting
% f66.fig = figure;
% widin = 5;  % maximum width allowed is 8.5 inches
% htin = 10;   % maximum height allowed is 11 inches
% set(f66.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
% nrow = 2;
% ncol = 1;
% % for isub = 1:nrow*ncol
% %     f66.ax(isub) = subplot(nrow,ncol,isub);
% % end
% 
% subplot(nrow,ncol,1)
% theta = angbest;
% theta = deg2rad(theta);
% thetaelse = theta(indelse);
% 
% binw = 5;
% binedge = deg2rad(2.5: binw: (92.5 + 270));
% bincnt = deg2rad(5:binw:360);
% nbin = round(360/5);
% counts = zeros(nbin, 1);
% for i = 1: length(indelse)
%     num = numlf(indelse(i));
%     for j = 1: nbin
%         if thetaelse(i)>bincnt(j)-deg2rad(binw/2) && thetaelse(i)<=bincnt(j)+deg2rad(binw/2)
%             break
%         else
%             continue
%         end
%     end
%     counts(j) = counts(j) + num;
% end
% sum(counts)
% sum(numlf(indelse))
% 
% polarhistogram('binedges',binedge,'BinCounts',counts,...
%                'facec',[0.6 0.6 0.6],'facea',0.8); hold on
% ax = gca;
% ax.RTickLabelRotation = 45;
% ax.ThetaMinorTick = 'on';
% ax.TickLength = [0.02 0];
% ax.ThetaZeroLocation = 'top';
% ax.ThetaDir = 'clockwise';
% ax.FontSize = 9;
% % ax.RTick = 0:1:5;
% % ax.LineWidth = 1.5;
% ax.Box = 'on';
% ax.GridAlpha = 0.3;
% text(ax,0,(ax.RLim(2))*9/10,'N','HorizontalAlignment',"center",'FontSize',12);
% text(ax,deg2rad(90),(ax.RLim(2))*9/10,'E','HorizontalAlignment',"center",'FontSize',12);
% text(ax,deg2rad(180),(ax.RLim(2))*9/10,'S','HorizontalAlignment',"center",'FontSize',12);
% text(ax,deg2rad(270),(ax.RLim(2))*9/10,'W','HorizontalAlignment',"center",'FontSize',12);
% % text(ax,deg2rad(315),7.4,'b','HorizontalAlignment',"center",'FontSize',11,'EdgeColor','k','Margin',2);
% % text(ax,deg2rad(45),7.3,'TWKB','HorizontalAlignment',"center",'FontSize',11,'EdgeColor','k','Margin',2);
% 
% 
% subplot(nrow,ncol,2)
% polarhistogram(theta,'binedges',binedge,...
%                'facec',[0.6 0.6 0.6],'facea',0.8); hold on
% ax = gca;
% ax.RTickLabelRotation = 45;
% ax.ThetaMinorTick = 'on';
% ax.TickLength = [0.02 0];
% ax.ThetaZeroLocation = 'top';
% ax.ThetaDir = 'clockwise';
% ax.FontSize = 9;
% % ax.RTick = 0:1:5;
% % ax.LineWidth = 1.5;
% ax.Box = 'on';
% ax.GridAlpha = 0.3;
% text(ax,0,(ax.RLim(2))*9/10,'N','HorizontalAlignment',"center",'FontSize',12);
% text(ax,deg2rad(90),(ax.RLim(2))*9/10,'E','HorizontalAlignment',"center",'FontSize',12);
% text(ax,deg2rad(180),(ax.RLim(2))*9/10,'S','HorizontalAlignment',"center",'FontSize',12);
% text(ax,deg2rad(270),(ax.RLim(2))*9/10,'W','HorizontalAlignment',"center",'FontSize',12);

%%% save figure
% print(f66.fig,'-dpdf',strcat(rstpath,'/autortmv3.other_direc_wtnum',num2str(nfam),'.pdf'));


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
vdistne = vdistne;     % keep the sign for NE group, check note above
ivdistne = ivdistne;

vdistsw2 = -vdistsw2;     % flip the sign for SW group, check note above
ivdistsw2 = -ivdistsw2;
vdistne2 = vdistne2;     % keep the sign for NE group, check note above
ivdistne2 = ivdistne2;

vdistelse = vdistelse;     % keep the sign for other RTMs, so negative always mean LF lags
ivdistelse = ivdistelse;


%% plot one histogram for the opposing rtm wrt. migrating front, like fig c5
%%% All detections
xran = [-15 10];
yran = [-15 15];
binw = 0.5;
conf = 99;
[f,barsw,pdfxlocsw,pdfvalsw,~,~,~,~] = ...
    plt_hist_combined_RTMs([],indsw,ranvechf98,medallhf,angbest,xran,yran,vdistsw,wtsw,...
                           [],[],binw,conf,'int');
text(f.ax(1),0.5,0.06,'All LF detections','FontSize',12,'unit','normalized',...
     'horizontalalignment','center');
% text(f.ax(1),0.5,0.12,'Main LZB region,','FontSize',12,'unit','normalized',...
%      'horizontalalignment','center');
xlim(f.ax(2), [-12.5 12.5]); 
xvect = [-7 -12.5+0.5];
yvect = [0.015 0.015];
drawArrow(f.ax(2),xvect,yvect,f.ax(2).XLim,f.ax(2).YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
xvect = [7 12.5-0.5];
yvect = [0.015 0.015];
drawArrow(f.ax(2),xvect,yvect,f.ax(2).XLim,f.ax(2).YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
text(f.ax(2),0.05,0.1,strcat({'WSW'}),'fontsize',11,'unit','normalized');
text(f.ax(2),0.95,0.1,strcat({'ENE'}),'fontsize',11,'unit','normalized',...
     'horizontalalignment','right');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/autortmv3.vdist_mainSW_all',num2str(nfam),'.pdf'));

keyboard

%%
%%% fc is inverted with (lf-hf_12, lf-hf_13) relative to all fams
%%% the same content as /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/LZBallfamfc, which is
%%% inverted by hypoinverse from /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/chat00p.fc
fclzb = [-123.509500 48.457000 38.0200;
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
contlzb = [-123.492667 48.451500 38.1400;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% convert each filterring correction to km relative to its own fam
fcveclzb = zeros(size(fclzb,1),2);
fcmaglzb = zeros(size(fclzb,1),1);
for i = 1: size(fclzb,1)
    [dx, dy] = absloc2relaloc(fclzb(i,1),fclzb(i,2),contlzb(i,1),contlzb(i,2));
    fcveclzb(i,:) = [dx dy];
    fcmaglzb(i) = sqrt(dx.^2+dy.^2);
end
fcazilzb = zeros(size(fclzb,1),1);
for i = 1: nfam
    tmp = rad2deg(atan2(fcveclzb(i,1),fcveclzb(i,2)));
    if tmp < 0
        fcazilzb(i) = 360+tmp;
    else
        fcazilzb(i) = tmp;
    end
end


%%% convert lfe location to km relative to 043
loclfe = zeros(size(contlzb,1),2);
for i = 1: size(contlzb,1)
    [dx, dy] = absloc2relaloc(contlzb(i,1),contlzb(i,2),contlzb(2,1),contlzb(2,2));
    loclfe(i,:) = [dx dy];
end

dtdir = ('/home/data2/chaosong/matlab/allan/2004/JUL/');
flist= ('datalist');
outf = ('staloc.txt');
stainfo = getstainfo(dtdir, flist, outf);
ind=[3,11,4,5,8,6];
stainfo = stainfo(ind,:);

stalat = str2double(stainfo(:,2));
stalon = str2double(stainfo(:,3));
[staloc(:,1),staloc(:,2)] = absloc2relaloc(stalon,stalat,contlzb(2,1),contlzb(2,2));

%%% load the location resolution points, +-2 samples
reslzb = load(fullfile(rstpath,'evtloc.13fam_locres_hf2spl'));
ifam = 2;   % choose fam 043 as the representative
reslzb = reslzb((ifam-1)*4+1: ifam*4, :);

% vecreslzb = zeros(size(reslzb,1),2);
[dx, dy] = absloc2relaloc(reslzb(:,1),reslzb(:,2),contlzb(ifam,1),contlzb(ifam,2));
vecreslzb = [dx dy];

% keyboard
%%
%%% load the location error grid, +-2 samples
errlzb = load(fullfile(rstpath,'evtloc.offset_043_locerr2_40_gaus_sq'));
% errlzb = load(fullfile(lzbpath,'evtloc.offset_043_locerr2_40_unif'));
[dx, dy] = absloc2relaloc(errlzb(:,1),errlzb(:,2),contlzb(ifam,1),contlzb(ifam,2));
vecerrlzb = [dx dy];

% keyboard
%%% load the offset error grid in samples
errlzbspl = load(fullfile(rstpath, 'offset_043_locerr2_40_gaus_sq'));

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

meansw = wt_mean(vdistsw,wtsw); % weighted mean of single group
meanne = wt_mean(vdistne,wtne); % weighted mean of single group
meanlzb = wt_mean(vdistlzb,wtlzb);  % weighted mean of combined groups
mswnedif = meansw - meanne;    % note that it is still negative for empirical FEC, from fig7
mswneave = (meansw + meanne)/2;     % direct average, ignore the different weights
% if talking about the separation in geo coord. as in fig 7, then we can view WSW as the dominant
% direction of WSW group, and E rather than ENE as the dominant direction of ENE group, then the
% offset of E group has to projected onto the ENE direction. In fig 7 we didn't do that as we
% thought they are ideally opposite
mswnedif2 = wt_mean(-vdistsw,wtsw) - wt_mean(vdistne,wtne); % separation in map view as in fig 7
mswnedif2proj = wt_mean(-vdistsw,wtsw) - wt_mean(vdistne,wtne)*cosd(30); % separation in map view projected

% keyboard
%%% for each black dot, we will create a new version of the histogram like fig c5, by using the 
%%% difference between the new and old FEC along migration propagation, to correct each LF-HF 
%%% measurement in each of the migrations, then remake the fig c5, in the reference frame of the
% vdistelseAll = [];
% wtelseAll = [];
% vdistlzbAll = [];
% wtlzbAll = [];
for i = 1: size(vecerrlzb,1)
    x = vecerrlzb(i,1);
    y = vecerrlzb(i,2);
    fcmag(i) = sqrt(x.^2+y.^2);
    tmp = rad2deg(atan2(x, y));
    if tmp < 0
        fcazi(i) = 360+tmp;
    else
        fcazi(i) = tmp;
    end
    
    %%% for other unused migrations in fig 8b
    vdistelsen = [];
    wtelsen = wtelse;   % weight is the same, as it is a bulk shift to the regression
    for j = 1: length(indelse)
        % empirical FEC along the mig j
        azi1 = angbest(indelse(j));
        fcmag1rep = fcmagrep*cos(deg2rad(360-fcazirep+azi1));
        % possible FEC i along the mig j
        fcmag1 = fcmag(i)*cos(deg2rad(360-fcazi(i)+azi1));
        % difference between FEC
        fcdiff = fcmag1rep - fcmag1;
        % new LF-HF measurement in the same migration, new = old+fcdiff
        temp1 = vertdist(1:numlf(indelse(j)), indelse(j));
        temp1 = temp1 + fcdiff;
        vdistelsen = [vdistelsen; temp1];
        
    end
    % get the new mean and CI like fig c5 for other migs in fig 8b, using updated LF-HF measurements  
    meanelsen(i) = wt_mean(vdistelsen,wtelsen);
    sigmaelsen(i) = sqrt(wt_var(vdistelsen,wtelsen,2));
    neffn = sum(wtelsen);
    conf = 99;
    CIelsen(i, :) = confidence_interval(meanelsen(i),sigmaelsen(i),neffn,conf);
%     % combine lf-hf measurements from each new FEC
%     vdistelseAll = [vdistelseAll; vdistelsen];
%     wtelseAll = [wtelseAll; wtelsen];
    
    
    %%% First combine them for the quasi-opposite migrations in fig 8a
    vdistlzbn = [];
    wtlzbn = wtlzb;   % weight is the same, as it is a bulk shift to the regression
    for j = 1: length(indlzb)
        % empirical FEC along the mig j
        azi1 = angbest(indlzb(j));
        fcmag1rep = fcmagrep*cos(deg2rad(360-fcazirep+azi1));
        % possible FEC i along the mig j
        fcmag1 = fcmag(i)*cos(deg2rad(360-fcazi(i)+azi1));
        % difference between FEC
        fcdiff = fcmag1rep - fcmag1;
        % new LF-HF measurement in the same migration, new = old+fcdiff
        temp1 = vertdist(1:numlf(indlzb(j)), indlzb(j));
        temp1 = temp1 + fcdiff;
        vdistlzbn = [vdistlzbn; temp1];
    end
    
    % get the new mean and CI like fig c5 for other migs in fig 8b, using updated LF-HF measurements  
    meanlzbn(i) = wt_mean(vdistlzbn,wtlzbn);
    sigmalzbn(i) = sqrt(wt_var(vdistlzbn,wtlzbn,2));
    neffn = sum(wtlzbn);
    conf = 99;
    CIlzbn(i, :) = confidence_interval(meanlzbn(i),sigmalzbn(i),neffn,conf);
%     % combine lf-hf measurements from each new FEC
%     vdistlzbAll = [vdistlzbAll; vdistlzbn];
%     wtlzbAll = [wtlzbAll; wtlzbn];
    
    %%% Second separate them for the quasi-opposite migrations in fig 8a
    vdistswn = [];
    wtswn = wtsw;   % weight is the same, as it is a bulk shift to the regression
    for j = 1: length(indsw)
        % empirical FEC along the mig j
        azi1 = angbest(indsw(j));
        fcmag1rep = fcmagrep*cos(deg2rad(360-fcazirep+azi1));
        % possible FEC i along the mig j
        fcmag1 = fcmag(i)*cos(deg2rad(360-fcazi(i)+azi1));
        % difference between FEC
        fcdiff = fcmag1rep - fcmag1;
        % new LF-HF measurement in the same migration, new = old+fcdiff
        temp1 = vertdist(1:numlf(indsw(j)), indsw(j));
        temp1 = temp1 + fcdiff;
        vdistswn = [vdistswn; temp1];
    end
    vdistnen = [];
    wtnen = wtne;   % weight is the same, as it is a bulk shift to the regression
    for j = 1: length(indne)
        % empirical FEC along the mig j
        azi1 = angbest(indne(j));
        fcmag1rep = fcmagrep*cos(deg2rad(360-fcazirep+azi1));
        % possible FEC i along the mig j
        fcmag1 = fcmag(i)*cos(deg2rad(360-fcazi(i)+azi1));
        % difference between FEC
        fcdiff = fcmag1rep - fcmag1;
        % new LF-HF measurement in the same migration, new = old+fcdiff
        temp1 = vertdist(1:numlf(indne(j)), indne(j));
        temp1 = temp1 + fcdiff;
        vdistnen = [vdistnen; temp1];
    end
    % get the new mean and CI like fig c5 for other migs in fig 8b, using updated LF-HF measurements  
    meanswn(i) = wt_mean(vdistswn,wtswn); % weighted mean of single group
    meannen(i) = wt_mean(vdistnen,wtnen); % weighted mean of single group
    mswnendif(i) = meanswn(i)-meannen(i); % note that it is still negative for empirical FEC, from fig7
    mswnenave(i) = (meanswn(i)+meannen(i))/2; % direct average, ignore the different weights
    mswnendif2(i) = wt_mean(-vdistswn,wtswn) - wt_mean(vdistnen,wtnen); % separation in map view as in fig 7
    mswnendif2proj(i) = wt_mean(-vdistswn,wtswn) - wt_mean(vdistnen,wtnen)*cosd(30); % separation in map view projected

    % obtain the max. difference between the three groups
    meansort = sort([meanswn(i) meannen(i) meanelsen(i)]);
    maxdif(i) = meansort(end)-meansort(1);
end
% keyboard

%% plot the components of 9 corrections projected along the azimuth of rep. fam, 043
famcheck = ['043';
            '068';
            '147';
            '141';
            '099';
            '006';
            '125';
            '017';
            '144';];

fcmaglzbrep = zeros(size(famcheck,1),1);
for i = 1: size(famcheck,1)
    [~,ind] = ismember(famcheck(i,:),nfampool,'rows');
    % FEC of fam i of interest along the FEC 0f fam 043
    fcmaglzbrep(i) = fcmaglzb(ind)*cos(deg2rad(360-fcazilzb(ind)+fcazirep));
end
% figure
% ax = gca;
% hold(ax, 'on');
mfcmaglzbrep = mean(fcmaglzbrep);
sigfcmaglzbrep = std(fcmaglzbrep);
neffn = length(fcmaglzbrep);
conf = 95;
CIfcmaglzbrep = confidence_interval(mfcmaglzbrep,sigfcmaglzbrep,neffn,conf);
binw = 0.5;
% histogram(ax,fcmaglzbrep,'binw',binw,'facec',[0.6 0.6 0.6],'facea',0.6);  %,'normalization','pdf'[0 1 1]
[mu1,sigma1,muCI,sigmaCI] = normfit(fcmaglzbrep, 0.01);
% xx = mu1+ax.XLim(1): 0.01: mu1+ax.XLim(2);
% synpdf = normpdf(xx,mu1,sigma1);
% synct = synpdf*binw*neffn;      % synthetic counts
% plot(ax,xx, synct, 'b-', 'linew', 1.5);

% plot(ax,[mfcmaglzbrep mfcmaglzbrep],ax.YLim,'k--','linew',1.5);    % median of wt dist
% errorbar(ax,mfcmaglzbrep,ax.YLim(2)*8/9,CIfcmaglzbrep(1)-mfcmaglzbrep,CIfcmaglzbrep(2)-mfcmaglzbrep,...
%         'horizontal','o','markersize',3,'color',...
%          'k','linew',1,'MarkerEdgeColor','k','MarkerFaceColor','k');
% hold(ax, 'off');

figure
ax = gca;
hold(ax, 'on');
scatter(ax,fcmaglzbrep,ones(neffn,1),20,'filled');
plot(ax,[mfcmaglzbrep mfcmaglzbrep],ax.YLim,'k--','linew',1.5);    % median of wt dist
errorbar(ax,mfcmaglzbrep,1,CIfcmaglzbrep(1)-mfcmaglzbrep,CIfcmaglzbrep(2)-mfcmaglzbrep,...
        'horizontal','o','markersize',3,'color',...
         'k','linew',1,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold(ax, 'off');

close all


%% summary of results with error bars from lf fitting
CIFcn = @(x,p)std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * tinv(abs([0,1]-(1-p/100)/2),...
              sum(~isnan(x(:)))-1) + mean(x(:),'omitnan');
f5.fig=figure;
f5.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 6.5;   % maximum height allowed is 11 inches
set(f5.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol-1
    f5.ax(isub) = subplot(nrow,ncol,isub);
end

set(f5.ax(1), 'position', [ 0.1, 0.1, 0.4, 0.85]);
set(f5.ax(2), 'position', [ 0.53, 0.1, 0.4, 0.85]);

% % index of main LZB region RTMs
indneplt = [10,18,21,30,34,46,47];
indswplt = [5,6,7,9,12,13,15,22,23,28,33,35,37,38,41,43,44,49,51,52];
indlzbplt = [fliplr(indneplt), fliplr(indswplt)];

% % index of 002 region RTMs
% ind002ne = [21,22];
% ind002sw = [4,25];
% ind002 = [fliplr(ind002ne), fliplr(ind002sw)];
% ind002 = [];
% 
% % index of all other RTMs
% indall = 1: 1: size(trange,1);
% tmp = union(indlzb,ind002);
% indelse = setdiff(indall, tmp);
indelseplt = fliplr(indelse);

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
ii = length(indlzbplt);
indplt = indlzbplt;
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
    if i > 7
        ypos = ii-i;
    end
    if ~isempty(intersect(indne,indplt(i)))   % NE
        xpos = offset(indplt(i));   % sign unchanged
        e1 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[0 1 1],'CapSize',5);
    elseif ~isempty(intersect(indsw,indplt(i))) % SW
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
text(ax,0.04,0.96,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.96,0.96,'Oppositely-propagating RTMs','FontSize',10,'unit','normalized',...
     'horizontalalignment','right','EdgeColor','k','Margin',2);

yran = [-3,30];
xran = [-5,5];
% axis(ax,[xran yran]);
xvect = [-1.5 -4.5];
yvect = [-2 -2];
drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
xvect = [1.5 4.5];
yvect = [-2 -2];
drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
text(ax,0.2,0.05,strcat({'WSW'}),'fontsize',10,'fontw','bold','unit','normalized',...
     'horizontalalignment','center');
text(ax,0.8,0.05,strcat({'ENE'}),'fontsize',10,'fontw','bold','unit','normalized',...
     'horizontalalignment','center');
 
drawArrow(ax,[-4.2 -1],[5 5],xran,yran,'linewidth',1,'linestyle','-','color',[0.65 0.65 0.65]);
text(ax,-4.8,5,'I','FontSize',13,'VerticalAlignment','middle'); 
drawArrow(ax,[-4.2 -1],[13 13],xran,yran,'linewidth',1,'linestyle','-','color',[0.65 0.65 0.65]);
text(ax,-4.8,13,'II','FontSize',13,'VerticalAlignment','middle');
drawArrow(ax,[-4.2 -1],[16 16],xran,yran,'linewidth',1,'linestyle','-','color',[0.65 0.65 0.65]);
text(ax,-4.8,16,'III','FontSize',13,'VerticalAlignment','middle');
drawArrow(ax,[-4.2 -1],[21 21],xran,yran,'linewidth',1,'linestyle','-','color',[0.65 0.65 0.65]);
text(ax,-4.8,21,'IV','FontSize',13,'VerticalAlignment','middle');

 
xlabel(ax,'LF-HF (km)','fontsize',11);
% ylabel(ax,'RTM number in chronological order','fontsize',11);
legend(ax,[e1, e3],{'ENE propagation','WSW propagation'},'FontSize',7,...
       'Position',[0.156 0.63 0.065 0.05]);       % ,'edgecolor',[0.5 0.5 0.5]
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
ii = length(indlzbplt);
indplt = indelseplt;
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
    if angbest(indplt(i))>0 && angbest(indplt(i))<=90       % NE
        xpos = offset(indplt(i));   % sign unchanged
        e1 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[0 1 1],'CapSize',5);
    elseif angbest(indplt(i))>90 && angbest(indplt(i))<=180     % SE
        xpos = offset(indplt(i));   % sign unchanged
        e2 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor','k','CapSize',5);
    elseif angbest(indplt(i))>180 && angbest(indplt(i))<=270    % SW
        xpos = offset(indplt(i));   % sign unchanged
        e3 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[1 96/255 0],'CapSize',5);
    else
        xpos = offset(indplt(i));   % sign unchanged
        e4 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5],'CapSize',5);        
    end
    text(ax,xpos+0.1,ypos+0.4,num2str(indplt(i)),'fontsize',7,...
        'HorizontalAlignment','left');
    text(ax,4.9,ypos,strcat(num2str(hfN),'; ',num2str(lfN)),'fontsize',7,...
        'HorizontalAlignment','right');
end

text(ax,0.04,0.96,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.96,0.96,'Other RTMs','FontSize',10,'unit','normalized',...
     'horizontalalignment','right','EdgeColor','k','Margin',2);
% text(ax,0.7,0.97,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.2,0.07,'LF','FontSize',10,'unit','normalized','HorizontalAlignment','center','fontw',...
    'bold');
text(ax,0.2,0.04,'lags','FontSize',10,'unit','normalized','HorizontalAlignment','center','fontw',...
    'bold');
text(ax,0.8,0.07,'LF','FontSize',10,'unit','normalized','HorizontalAlignment','center','fontw',...
    'bold');
text(ax,0.8,0.04,'leads','FontSize',10,'unit','normalized','HorizontalAlignment','center','fontw',...
    'bold');
text(ax,0.55,0.15,'(RTMs #1, #26 and #27 are excluded)','FontSize',9,'unit','normalized',...
    'HorizontalAlignment','center');
yran = [0 30];
xran = [-5 5];
% drawArrow(ax,[-4.2 -2.5],[5 5],xran,yran,'linewidth',1.5,'linestyle','-','color',[0.65 0.65 0.65]);
% text(ax,-4.8,5,'IV','FontSize',13,'VerticalAlignment','middle','backgroundcolor','w',...
%      'Margin',0.5);
axis(ax,[xran yran]);
xlabel(ax,'LF-HF (km)','fontsize',11);
% ylabel(ax,'RTM number in chronological order','fontsize',11);
legend(ax,[e1, e2, e3, e4],{'NE propagation','SE propagation','SW propagation', 'NW propagation'},...
       'FontSize',7,'Position',[0.58 0.65 0.065 0.05]);  
xticks(ax,-5:1:5);
ax.YTick = [];
% yticks(0:1:8);
ax.Box='on';
grid(ax,'on');
set(ax,'GridLineStyle','--');
ax.FontSize = 9;
hold(ax,'off');


% the migration 40, LF is too bad, nonsense offset and error, in order to show it as well  
set(f5.ax(3), 'position', [ 0.81, 0.748, 0.06, 0.02]);
ax = f5.ax(3);
hfN = numhf(40);
lfN = numlf(40);
lf=resproplf(40,1:lfN);
lfm = mean(lf);
lfstd = std(lf);
lfsem = lfstd/sqrt(lfN);
ci95 = tinv([0.025 0.975], lfN-1);
offcor = bsxfun(@times, lfsem, ci95(:));
xpos = offset(40);   % sign unchanged
errorbar(ax,xpos,0.25,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor','k','CapSize',5);
text(ax,xpos+0.1,0.5,num2str(40),'fontsize',7,...
        'HorizontalAlignment','left');
% xlabel(ax,'LF-HF (km)','fontsize',6);    
yran = [0 0.5];
xran = [5 25];
axis(ax,[xran yran]);
xticks(ax,[5 15 25]);
ax.YTick = [];
% yticks(0:1:8);
ax.Box='off';
ax.FontSize = 6;
hold(ax,'off');

%%% save figure
print(f5.fig,'-dpdf',strcat(rstpath,'/LZB',num2str(nfam),'autortmv3.mig.lfit.sum.exc3mig',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));

keyboard

%% plot the variation of the other rtm
% f.fig=figure;
% f.fig.Renderer='Painters';
% widin = 8;  % maximum width allowed is 8.5 inches
% htin = 5;   % maximum height allowed is 11 inches
% set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
% 
% nrow = 1;
% ncol = 2;
% for isub = 1:nrow*ncol
%     f.ax(isub) = subplot(nrow,ncol,isub);
%     f.ax(isub).Box='on';
%     grid(f.ax(isub),'on');
%     f.ax(isub).GridLineStyle = '--';
% end
% 
% % % reposition
% % set(f.ax(1), 'position', [ 0.1, 0.56, 0.35, 0.35]);
% % set(f.ax(2), 'position', [ 0.55, 0.56, 0.35, 0.35]);
% % set(f.ax(3), 'position', [ 0.1, 0.1, 0.34, 0.41]);
% % set(f.ax(4), 'position', [ 0.55, 0.1, 0.35, 0.41]);
% 
% ax = f.ax(1);
% hold(ax,'on');
% xlim(ax, [-3 3]);
% xticks(ax, -3: 1: 3);
% ylim(ax, [0 600]);
% 
% binw = 0.1;
% 
% h1=histogram(ax,meanelsen,'binw',binw,'facec',[0.6 0.6 0.6],'facea',0.6);  %,'normalization','pdf'[0 1 1]
% % plot(ax, [fcdiffrep fcdiffrep],ax.YLim,'b--','linewidth',2);
% [mu1,sigma1,muCI,sigmaCI] = normfit(meanelsen', 0.01);
% xx = mu1+ax.XLim(1): 0.01: mu1++ax.XLim(2);
% synpdf = normpdf(xx,mu1,sigma1);
% synct = synpdf*binw*length(meanelsen);      % synthetic counts
% plot(ax,xx, synct, 'b-', 'linew', 1.5);
% 
% %%% every meanelsen here has the same weight of 1
% mmean = mean(meanelsen);
% sigmean = std(meanelsen);
% neffnmean = length(meanelsen);
% conf = 95;
% CImean = confidence_interval(mmean,sigmean,neffnmean,conf);
% 
% plot(ax,[mmean mmean],ax.YLim,'k--','linew',1.5);    % median of wt dist
% errorbar(ax,mmean,ax.YLim(2)*8/9,CImean(1)-mmean,CImean(2)-mmean,'horizontal','o','markersize',3,'color',...
%          'k','linew',1,'MarkerEdgeColor','k','MarkerFaceColor','k');
% % lgd = legend(ax,[h1,h2],{'E-WSW','SE-SSE'},'fontsize',8,...
% %        'numcolumns',1,...
% %        'location','northeast');
% text(ax, 0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% text(ax, 0.96,0.7,strcat({'$\mu$='},sprintf('%.2f',mmean),' km'),'FontSize',13,...
%      'unit','normalized','Interpreter','latex','color','b','horizontalalignment','right');
% text(ax, 0.96,0.6,strcat({'$\sigma$='},sprintf('%.2f',sigma1),' km'),'FontSize',13,...
%      'unit','normalized','Interpreter','latex','color','b','horizontalalignment','right');
% text(ax, 0.96,0.93,'Other RTMs','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1,...
%      'horizontalalignment','right'); 
% xlabel(ax,'mean LF-HF (km)','fontsize',12)
% ylabel(ax,'Count','fontsize',12)
% % if ifam == 2 
% %     xlim(ax, [-6 0]);
% % elseif ifam == 10
% %     xlim(ax, [-5 1]);
% % end
% hold(ax,'off');
% 
% 
% ax = f.ax(2);
% hold(ax,'on');
% xlim(ax, [0 1]);
% xticks(ax, 0: 0.2: 1);
% % ylim(ax, [0 2000]);
% 
% binw = 0.02;
% 
% CIwid = CIelsen(:,2)-CIelsen(:,1);
% h1=histogram(ax,CIwid,'binw',binw,'facec',[0.6 0.6 0.6],'facea',0.6);  %,'normalization','pdf'[0 1 1]
% % plot(ax, [fcdiffrep fcdiffrep],ax.YLim,'b--','linewidth',2);
% [mu2,sigma2,muCI,sigmaCI] = normfit(CIwid', 0.01);
% % xx = mu2+ax.XLim(1): 0.01: mu2++ax.XLim(2);
% xx = -6:0.01:6;
% synpdf = normpdf(xx,mu2,sigma2);
% synct = synpdf*binw*length(CIwid);      % synthetic counts
% % plot(ax,xx, synct, 'b-', 'linew', 1.5);
% 
% %%% every meanelsen here has the same weight of 1
% mCIwid = mean(CIwid);
% sigCIwid = std(CIwid);
% neffnCIwid = length(CIwid);
% conf = 95;
% CICIwid = confidence_interval(mCIwid,sigCIwid,neffnCIwid,conf);
% 
% % plot(ax,[mCIwid mCIwid],ax.YLim,'k--','linew',1.5);    % median of wt dist
% % errorbar(ax,mCIwid,ax.YLim(2)*8/9,CICIwid(1)-mCIwid,CICIwid(2)-mCIwid,'horizontal','o','markersize',3,'color',...
% %          'k','linew',1,'MarkerEdgeColor','k','MarkerFaceColor','k');
% text(ax, 0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% % text(ax, 0.96,0.7,strcat({'$\mu$='},sprintf('%.2f',mCIwid),' km'),'FontSize',13,...
% %      'unit','normalized','Interpreter','latex','color','b','horizontalalignment','right');
% % text(ax, 0.96,0.6,strcat({'$\sigma$='},sprintf('%.2f',sigma2),' km'),'FontSize',13,...
% %      'unit','normalized','Interpreter','latex','color','b','horizontalalignment','right');
% text(ax, 0.96,0.93,'Other RTMs','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1,...
%      'horizontalalignment','right'); 
% xlabel(ax,'Width of 99% CI of mean LF-HF (km)','fontsize',12)
% ylabel(ax,'Count','fontsize',12)
% hold(ax,'off');
%                                  
% %%% save figure
% print(f.fig,'-dpdf',strcat(rstpath,'/autortmv3.mean_LF-HF_FECerror_otherrtm.pdf'));
% % keyboard

%% plot one histogram for the opposing rtm wrt. migrating front, like fig c5
% %%% All detections
% xran = [-15 10];
% yran = [-15 15];
% binw = 0.5;
% conf = 99;
% [f,barlzb,pdfxloclzb,pdfvallzb] = ...
%     plt_hist_combined_RTMs_otherrtm(indlzb,ranvechf98,medallhf,angbest,xran,yran,vdistlzb,wtlzb,...
%                            binw,conf,'int');
% text(f.ax(1),0.5,0.06,'All detections','FontSize',12,'unit','normalized',...
%      'horizontalalignment','center');
% ylim(f.ax(2),[0 0.25]);
% delete(f.ax(4))  
% 
% %%% save figure
% print(f.fig,'-dpdf',strcat(rstpath,'/autortmv3.vdist_main_all',num2str(nfam),'.pdf'));
% close all

%% replot fig c5, add inset of rose diagram
%%% All detections
xran = [-15 25];
yran = [-20 15];
binw = 0.5;
conf = 99;
[f88,barelse,pdfxlocelse,pdfvalelse] = ...
    plt_hist_combined_RTMs_otherrtm(indelse,ranvechf98,medallhf,angbest,xran,yran,vdistelse,wtelse,...
                           binw,conf,'int');
hold(f88.ax(1), 'on');
for i = 1: length(indexc)
    [rotx1, roty1] = complex_rot(0, ranvechf98(indexc(i),2)-medallhf(indexc(i),3), ...
                                 -angbest(indexc(i)));
    [rotx2, roty2] = complex_rot(0, medallhf(indexc(i),3)-ranvechf98(indexc(i),1), ...
                                 -angbest(indexc(i)));
    xvect = [medallhf(indexc(i),1)-rotx2 medallhf(indexc(i),1)+rotx1];
    yvect = [medallhf(indexc(i),2)-roty2 medallhf(indexc(i),2)+roty1];
    drawArrow(f88.ax(1),xvect,yvect,xran,yran,'linewidth',1,'linestyle','--','color','k');
end
for i = 1: length(indexc)
    scatter(f88.ax(1),medallhf(indexc(i),1),medallhf(indexc(i),2), 30, indexc(i), 'filled','o');  %, 'MarkerEdgeColor', 'w')
end                       
% text(f88.ax(1),0.43,0.06,'#1, #26 and #27 will be excluded','FontSize',11,'unit','normalized',...
%      'horizontalalignment','center');
% text(f88.ax(1),0.43,0.12,'Other RTMs,','FontSize',11,'unit','normalized',...
%      'horizontalalignment','center');
 
% keyboard
% add an inset of the rose diagram weighted by the LF number
patarea = [6 -3;
           24 -3;
           24 14;
           6 14;
           6 -3;];
patch(f88.ax(1),patarea(:,1),patarea(:,2),'w','edgecolor','none');
hold(f88.ax(1),'off');

ax1loc = f88.ax(1).Position;
set(f88.ax(4), 'position', [ax1loc(1)+ax1loc(3)/2+0.03, ax1loc(2)+ax1loc(4)/2+0.01, ...
                          ax1loc(3)*2/5-0.01, ax1loc(4)*2/5-0.01]);
ax = f88.ax(4);
theta = angbest;
theta = deg2rad(theta);
thetaelse = theta(indelse);
binw = 5;
binedge = deg2rad(2.5: binw: (92.5 + 270));
bincnt = deg2rad(5:binw:360);
nbin = round(360/5);
counts = zeros(nbin, 1);
for i = 1: length(indelse)
    num = numlf(indelse(i));
    for j = 1: nbin
        if thetaelse(i)>bincnt(j)-deg2rad(binw/2) && thetaelse(i)<=bincnt(j)+deg2rad(binw/2)
            break
        else
            continue
        end
    end
    counts(j) = counts(j) + num;
end
polarhistogram('binedges',binedge,'BinCounts',counts,...
               'facec',[0.7 0.7 0.7],'facea',0.6); hold on
ax = gca;
ax.RTickLabelRotation = 45;
ax.ThetaMinorTick = 'on';
ax.TickLength = [0.02 0];
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.FontSize = 7;
% ax.RTick = 0:1:5;
% ax.LineWidth = 1.5;
ax.Box = 'on';
ax.GridAlpha = 0.3;
hold(ax,'off');

hold(f88.ax(2),'on');
text(f88.ax(2),0.85,0.47,'(166 excluded)','fontsize',10,'unit','normalized',...
     'horizontalalignment','center','color','k');
hold(f88.ax(2),'off');

 
%%% save figure
print(f88.fig,'-dpdf',strcat(rstpath,'/autortmv3.vdist_otherrtm_exc3mig_',num2str(nfam),'.pdf'));
                                 
                       
%% plot the variation of the opposing rtm
% f.fig=figure;
% f.fig.Renderer='Painters';
% widin = 8;  % maximum width allowed is 8.5 inches
% htin = 5;   % maximum height allowed is 11 inches
% set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
% 
% nrow = 1;
% ncol = 2;
% for isub = 1:nrow*ncol
%     f.ax(isub) = subplot(nrow,ncol,isub);
%     f.ax(isub).Box='on';
%     grid(f.ax(isub),'on');
%     f.ax(isub).GridLineStyle = '--';
% end
% 
% ax = f.ax(1);
% hold(ax,'on');
% xlim(ax, [-3 3]);
% xticks(ax, -3: 1: 3);
% % ylim(ax, [0 600]);
% 
% binw = 0.1;
% h1=histogram(ax,meanlzbn,'binw',binw,'facec',[0.6 0.6 0.6],'facea',0.6);  %,'normalization','pdf'[0 1 1]
% % plot(ax, [fcdiffrep fcdiffrep],ax.YLim,'b--','linewidth',2);
% [mu1,sigma1,muCI,sigmaCI] = normfit(meanlzbn', 0.01);
% xx = mu1+ax.XLim(1): 0.01: mu1++ax.XLim(2);
% synpdf = normpdf(xx,mu1,sigma1);
% synct = synpdf*binw*length(meanlzbn);      % synthetic counts
% plot(ax,xx, synct, 'b-', 'linew', 1.5);
% 
% %%% every meanelsen here has the same weight of 1
% mmeanlzb = mean(meanlzbn);
% sigmeanlzb = std(meanlzbn);
% neffnmeanlzb = length(meanlzbn);
% conf = 95;
% CImeanlzb = confidence_interval(mmeanlzb,sigmeanlzb,neffnmeanlzb,conf);
% 
% plot(ax,[mmeanlzb mmeanlzb],ax.YLim,'k--','linew',1.5);    % median of wt dist
% errorbar(ax,mmeanlzb,ax.YLim(2)*8/9,CImeanlzb(1)-mmeanlzb,CImeanlzb(2)-mmeanlzb,'horizontal','o','markersize',3,'color',...
%          'k','linew',1,'MarkerEdgeColor','k','MarkerFaceColor','k');
% % lgd = legend(ax,[h1,h2],{'E-WSW','SE-SSE'},'fontsize',8,...
% %        'numcolumns',1,...
% %        'location','northeast');
% text(ax, 0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% text(ax, 0.96,0.7,strcat({'$\mu$='},sprintf('%.2f',mmeanlzb),' km'),'FontSize',13,...
%      'unit','normalized','Interpreter','latex','color','b','horizontalalignment','right');
% text(ax, 0.96,0.6,strcat({'$\sigma$='},sprintf('%.2f',sigma1),' km'),'FontSize',13,...
%      'unit','normalized','Interpreter','latex','color','b','horizontalalignment','right');
% text(ax, 0.96,0.93,'Opposing RTMs','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1,...
%      'horizontalalignment','right'); 
% xlabel(ax,'mean LF-HF (km)','fontsize',12)
% ylabel(ax,'Count','fontsize',12)
% % if ifam == 2 
% %     xlim(ax, [-6 0]);
% % elseif ifam == 10
% %     xlim(ax, [-5 1]);
% % end
% hold(ax,'off');
% 
% 
% ax = f.ax(2);
% hold(ax,'on');
% xlim(ax, [0 1]);
% xticks(ax, 0: 0.2: 1);
% % ylim(ax, [0 2000]);
% 
% binw = 0.02;
% 
% CIwid = CIlzbn(:,2)-CIlzbn(:,1);
% h1=histogram(ax,CIwid,'binw',binw,'facec',[0.6 0.6 0.6],'facea',0.6);  %,'normalization','pdf'[0 1 1]
% % plot(ax, [fcdiffrep fcdiffrep],ax.YLim,'b--','linewidth',2);
% [mu2,sigma2,muCI,sigmaCI] = normfit(CIwid', 0.01);
% % xx = mu2+ax.XLim(1): 0.01: mu2++ax.XLim(2);
% xx = -6:0.01:6;
% synpdf = normpdf(xx,mu2,sigma2);
% synct = synpdf*binw*length(CIwid);      % synthetic counts
% % plot(ax,xx, synct, 'b-', 'linew', 1.5);
% 
% %%% every meanelsen here has the same weight of 1
% mCIwid = mean(CIwid);
% sigCIwid = std(CIwid);
% neffnCIwid = length(CIwid);
% conf = 95;
% CICIwid = confidence_interval(mCIwid,sigCIwid,neffnCIwid,conf);
% 
% % plot(ax,[mCIwid mCIwid],ax.YLim,'k--','linew',1.5);    % median of wt dist
% % errorbar(ax,mCIwid,ax.YLim(2)*8/9,CICIwid(1)-mCIwid,CICIwid(2)-mCIwid,'horizontal','o','markersize',3,'color',...
% %          'k','linew',1,'MarkerEdgeColor','k','MarkerFaceColor','k');
% text(ax, 0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% % text(ax, 0.96,0.7,strcat({'$\mu$='},sprintf('%.2f',mCIwid),' km'),'FontSize',13,...
% %      'unit','normalized','Interpreter','latex','color','b','horizontalalignment','right');
% % text(ax, 0.96,0.6,strcat({'$\sigma$='},sprintf('%.2f',sigma2),' km'),'FontSize',13,...
% %      'unit','normalized','Interpreter','latex','color','b','horizontalalignment','right');
% text(ax, 0.96,0.93,'Opposing RTMs','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1,...
%      'horizontalalignment','right'); 
% xlabel(ax,'Width of 99% CI of mean LF-HF (km)','fontsize',12)
% ylabel(ax,'Count','fontsize',12)
% hold(ax,'off');
%                                  
% %%% save figure
% print(f.fig,'-dpdf',strcat(rstpath,'/autortmv3.mean_LF-HF_FECerror_main.pdf'));    
% % keyboard

%% What if we combine ALL lf-hf from ALL FEC into one histogram?
% f.fig=figure;
% f.fig.Renderer='Painters';
% widin = 8;  % maximum width allowed is 8.5 inches
% htin = 5;   % maximum height allowed is 11 inches
% set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
% 
% nrow = 1;
% ncol = 2;
% for isub = 1:nrow*ncol
%     f.ax(isub) = subplot(nrow,ncol,isub);
%     f.ax(isub).Box='on';
%     grid(f.ax(isub),'on');
%     f.ax(isub).GridLineStyle = '--';
% end
% 
% ax = f.ax(1);
% hold(ax,'on');
% ylim(ax, [0 0.3]);
% 
% meanlzbAll = wt_mean(vdistlzbAll,wtlzbAll);
% sigmalzbAll = sqrt(wt_var(vdistlzbAll,wtlzbAll,2));
% neffnAll = sum(wtlzbAll);
% conf = 99;
% CIlzbAll = confidence_interval(meanlzbAll,sigmalzbAll,neffnAll,conf);
% 
% patarea = [0 0;
%            -100 0;
%            -100 1;
%             0 1;
%             0 0];
% patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
% binw = 0.5;
% [ax, barlzbAll, pdfxloclzbAll, ~, ~, pdfvallzbAll] = plt_weighted_dist(ax, vdistlzbAll, wtlzbAll, ...
%     binw,'int');
% barlzbAll(1).FaceColor = [0.7 0.7 0.7];
% plot(ax,[meanlzbAll meanlzbAll],ax.YLim,'k--','linew',1.5);    % median of wt dist
% errorbar(ax,meanlzbAll,ax.YLim(2)*8/9,CIlzbAll(1)-meanlzbAll,CIlzbAll(2)-meanlzbAll,...
%          'horizontal','o','markersize',3,'color',...
%          'k','linew',1,'MarkerEdgeColor','k','MarkerFaceColor','k');
% text(ax, 0.96,0.7,strcat({'$\mu$='},sprintf('%.2f',meanlzbAll),' km'),'FontSize',13,...
%      'unit','normalized','Interpreter','latex','color','b','horizontalalignment','right');
% text(ax, 0.96,0.6,strcat({'$\sigma$='},sprintf('%.2f',sigmalzbAll),' km'),'FontSize',13,...
%      'unit','normalized','Interpreter','latex','color','b','horizontalalignment','right');
%  
% [mu1,sigma1,mu1CI,sigma1CI] = normfit(vdistlzbAll, 0.01,[], wtlzbAll);
% % xx = mu2+ax.XLim(1): 0.01: mu2++ax.XLim(2);
% xx = -20:0.01:20;
% synpdf = normpdf(xx,mu1,sigma1);
% synct = synpdf;      % synthetic counts
% plot(ax,xx, synct, 'k-', 'linew', 1.5);
% text(ax, 0.96,0.5,strcat({'$\mu$='},sprintf('%.2f',mu1),' km'),'FontSize',13,...
%      'unit','normalized','Interpreter','latex','color','k','horizontalalignment','right');
% text(ax, 0.96,0.4,strcat({'$\sigma$='},sprintf('%.2f',sigma1),' km'),'FontSize',13,...
%      'unit','normalized','Interpreter','latex','color','k','horizontalalignment','right');
% 
% % ylim(ax,[0 0.2]);
% xlim(ax, [-12.5 12.5]); 
% text(ax,0.05,0.15,strcat({'LF lags'}),'fontsize',11,'unit','normalized');
% text(ax,0.96,0.15,strcat({'LF leads'}),'fontsize',11,'unit','normalized',...
%      'horizontalalignment','right');
% text(ax, 0.96,0.93,'Opposing RTMs','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1,...
%      'horizontalalignment','right');
% hold(ax,'off');
% 
% 
% ax = f.ax(2);
% hold(ax,'on');
% ylim(ax, [0 0.3]);
% 
% meanelseAll = wt_mean(vdistelseAll,wtelseAll);
% sigmaelseAll = sqrt(wt_var(vdistelseAll,wtelseAll,2));
% neffnAll = sum(wtelseAll);
% conf = 99;
% CIelseAll = confidence_interval(meanelseAll,sigmaelseAll,neffnAll,conf);
% 
% patarea = [0 0;
%            -100 0;
%            -100 1;
%             0 1;
%             0 0];
% patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
% binw = 0.5;
% [ax, barelseAll, pdfxlocelseAll, ~, ~, pdfvalelseAll] = plt_weighted_dist(ax, vdistelseAll, wtelseAll, ...
%     binw,'int');
% barelseAll(1).FaceColor = [0.7 0.7 0.7];
% plot(ax,[meanelseAll meanelseAll],ax.YLim,'k--','linew',1.5);    % median of wt dist
% errorbar(ax,meanelseAll,ax.YLim(2)*8/9,CIelseAll(1)-meanelseAll,CIelseAll(2)-meanelseAll,...
%          'horizontal','o','markersize',3,'color',...
%          'k','linew',1,'MarkerEdgeColor','k','MarkerFaceColor','k');
% text(ax, 0.96,0.7,strcat({'$\mu$='},sprintf('%.2f',meanelseAll),' km'),'FontSize',13,...
%      'unit','normalized','Interpreter','latex','color','b','horizontalalignment','right');
% text(ax, 0.96,0.6,strcat({'$\sigma$='},sprintf('%.2f',sigmaelseAll),' km'),'FontSize',13,...
%      'unit','normalized','Interpreter','latex','color','b','horizontalalignment','right');
%  
% [mu2,sigma2,mu2CI,sigma2CI] = normfit(vdistelseAll, 0.01,[], wtelseAll);
% % xx = mu2+ax.XLim(1): 0.01: mu2++ax.XLim(2);
% xx = -20:0.01:20;
% synpdf = normpdf(xx,mu2,sigma2);
% synct = synpdf;      % synthetic counts
% plot(ax,xx, synct, 'k-', 'linew', 1.5);
% text(ax, 0.96,0.5,strcat({'$\mu$='},sprintf('%.2f',mu2),' km'),'FontSize',13,...
%      'unit','normalized','Interpreter','latex','color','k','horizontalalignment','right');
% text(ax, 0.96,0.4,strcat({'$\sigma$='},sprintf('%.2f',sigma2),' km'),'FontSize',13,...
%      'unit','normalized','Interpreter','latex','color','k','horizontalalignment','right');
% 
%  % ylim(ax,[0 0.2]);
% xlim(ax, [-12.5 12.5]); 
% % xvect = [-7 -12.5+0.5];
% % yvect = [0.015 0.015];
% % drawArrow(ax,xvect,yvect,ax.XLim,ax.YLim,'linewidth',1.5,'linestyle','-',...
% %           'color',[0.5 0.5 0.5]);
% % arrow1 = annotation('arrow','linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);
% % arrow1.Parent = f.ax(2);
% % arrow1.Position = [-7, 0.015, -5, 0] ;
% % xvect = [7 12.5-0.5];
% % yvect = [0.015 0.015];
% % drawArrow(ax,xvect,yvect,ax.XLim,ax.YLim,'linewidth',1.5,'linestyle','-',...
% %           'color',[0.5 0.5 0.5]);
% text(ax,0.05,0.15,strcat({'LF lags'}),'fontsize',11,'unit','normalized');
% text(ax,0.96,0.15,strcat({'LF leads'}),'fontsize',11,'unit','normalized',...
%      'horizontalalignment','right');
% text(ax, 0.96,0.93,'Other RTMs','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1,...
%      'horizontalalignment','right');
% % text(ax,0.85,0.62,num2str(length(vdistelse)),'fontsize',12,'unit','normalized',...
% %      'horizontalalignment','center','color','k');
% % text(ax,0.85,0.55,strcat({'detections'}),'fontsize',12,'unit','normalized',...
% %      'horizontalalignment','center','color','k');
% % text(f.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% % lgd = legend(ax,[barne,barsw],{'ENE propagation','WSW propagation'},'fontsize',7,...
% %        'numcolumns',2,...
% %        'Position',[0.62+0.35*0.048  c.Position(2)+ywid*2/3*0.93  0.35*0.9  ywid*2/3*0.05]);
% hold(ax,'off');
% 
% %%% save figure
% print(f.fig,'-dpdf',strcat(rstpath,'/autortmv3.LF-HF_AllFECerror_main_other.pdf'));    
% 
% keyboard
% close all

%% Create another summary figure to fit in these subplots
f.fig=figure;
f.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box='on';
    grid(f.ax(isub),'on');
    f.ax(isub).GridLineStyle = '--';
end

% reposition
set(f.ax(1), 'position', [ 0.1, 0.535, 0.35, 0.35]);
set(f.ax(2), 'position', [ 0.55, 0.535, 0.35, 0.35]);
set(f.ax(3), 'position', [ 0.1, 0.1, 0.35, 0.3885]);
set(f.ax(4), 'position', [ 0.55, 0.1, 0.35, 0.3885]);


%%% subplot 1, scatter of FEC from a common origin, excluding some different ones, i.e. 010, 047, 
%%% 002, 001
ax = f.ax(1);
hold(ax,'on');
xran = [-3.5 -0.5];
yran = [ 5 8];
axis(ax,'equal');

x = [0; fcvecrep(1)];
y = [0; fcvecrep(2)];
plot(ax,x,y,'color',[0.6 0.6 0.6],'linewidth',2.5);
text(ax,x(2)+0.1,y(2)+0.1,nfampool(ifam,:),'FontSize',11,'horizontalalignment','left');

famcheck = ['043';
            '068';
            '147';
            '141';
            '099';
            '006';
            '125';
            '017';
            '144';];

xx = [];
yy = [];
for i = 1: size(famcheck,1)
    [~,ind] = ismember(famcheck(i,:),nfampool,'rows');
    xx(i,1) = [fcveclzb(ind,1)];
    yy(i,1) = [fcveclzb(ind,2)];
    scatter(ax,xx(i),yy(i),20,'k','filled');
end

[k1,v1] = boundary(xx,yy,0);
plot(ax, xx(k1),yy(k1),'k','linew',1.5);

azi1 = 90;  % E
azix1 = [-2; -1];
azia1 = tand(90-azi1);
azib1 = 6.2-azia1*(-2);
aziy1 = linefcn(azix1,azia1,azib1);
% plot(ax, azix1,aziy1,'r','linew',1);
drawArrow(ax,azix1,aziy1,xran,yran,'linewidth',1,'linestyle','-','color','b');
text(ax, -1.1, 6.0, 'E', 'fontsize',11);

azi2 = 270-45/2;  % WSW
azix2 = [-2; -2.8];
azia2 = tand(90-azi2);
azib2 = 6.2-azia2*(-2);
aziy2 = linefcn(azix2,azia2,azib2);
% plot(ax, azix2,aziy2,'r','linew',1);
drawArrow(ax,azix2,aziy2,xran,yran,'linewidth',1,'linestyle','-','color','b');
text(ax, -3, 5.7, 'WSW', 'fontsize',11);

azi3 = 90-45/2;  % ENE
azix3 = [-2; -1.2];
azia3 = tand(90-azi3);
azib3 = 6.2-azia3*(-2);
aziy3 = linefcn(azix3,azia3,azib3);
plot(ax, azix3,aziy3,'b--','linew',1);
% drawArrow(ax,azix3,aziy3,xran,yran,'linewidth',1,'linestyle','--','color','r');
text(ax, -1.5, 6.7, '(ENE)', 'fontsize',11);

text(ax, 0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
xticks(ax, xran(1): 0.5: xran(2));
yticks(ax, yran(1): 0.5: yran(2));
xlabel(ax,'E (km)','fontsize',12)
ylabel(ax,'N (km)','fontsize',12)

hold(ax,'off');


%%% subplot 2, the FEC error ellipse in map view
ax = f.ax(2);
hold(ax,'on');
xran = [-7 5];
yran = [-1 11];
axis(ax,'equal');

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
plot(ax, vecerrlzb(k1,1),vecerrlzb(k1,2),'k','linew',1.5);

x = [0; fcvecrep(1)];
y = [0; fcvecrep(2)];
plot(ax,x,y,'color',[0.5 0.5 0.5],'linewidth',2.5);

azi1 = 90;  % E
azix1 = [fixpt(1); fixpt(1)+2];
azia1 = tand(90-azi1);
azib1 = fixpt(2)-azia1*(fixpt(1));
aziy1 = linefcn(azix1,azia1,azib1);
% plot(ax, azix1,aziy1,'r','linew',1);
drawArrow(ax,azix1,aziy1,xran,yran,'linewidth',1,'linestyle','-','color','b');
text(ax, 2.2, 0.1, 'E', 'fontsize',11);

azi2 = 270-45/2;  % WSW
azix2 = [fixpt(1); fixpt(1)-1.5];
azia2 = tand(90-azi2);
azib2 = fixpt(2)-azia2*(fixpt(1));
aziy2 = linefcn(azix2,azia2,azib2);
% plot(ax, azix2,aziy2,'r','linew',1);
drawArrow(ax,azix2,aziy2,xran,yran,'linewidth',1,'linestyle','-','color','b');
text(ax, -2.5, 0.1, 'WSW', 'fontsize',11);

azi3 = 90-45/2;  % ENE
azix3 = [fixpt(1); fixpt(1)+1.5];
azia3 = tand(90-azi3);
azib3 = fixpt(2)-azia3*(fixpt(1));
aziy3 = linefcn(azix3,azia3,azib3);
plot(ax, azix3,aziy3,'b--','linew',1);
% drawArrow(ax,azix3,aziy3,xran,yran,'linewidth',1,'linestyle','--','color','r');
text(ax, 1, 1.2, '(ENE)', 'fontsize',11);
    
text(ax, 0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax, 0.96,0.93,nfampool(ifam,:),'FontSize',11,'unit','normalized','EdgeColor','k','Margin',1,...
     'horizontalalignment','right');
xticks(ax, xran(1): 1: xran(2));
yticks(ax, yran(1): 1: yran(2));
xlabel(ax,'E (km)','fontsize',12)
ylabel(ax,'N (km)','fontsize',12)

hold(ax,'off');


%%% subplot 3, generating a new histogram like fig c5 for each FEC, plot the histogram of its
%%% mean, for both opposing rtm group and other unused rtm group
ax = f.ax(3);
hold(ax,'on');
xlim(ax, [-3 3]);
xticks(ax, -3: 1: 3);
ylim(ax, [0 820]);
patarea = [0 0;
           -100 0;
           -100 1000;
            0 1000;
            0 0];
patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');

binw = 0.1;
h2=histogram(ax,meanelsen,'binw',binw,'facec',[0.5 0.5 0.5],'facea',0.6);  %,'normalization','pdf'[0 1 1]
[mu2,sigma2,muCI,sigmaCI] = normfit(meanelsen', 0.01);
xx = ax.XLim(1): 0.01: ax.XLim(2);
synpdf2 = normpdf(xx,mu2,sigma2);
synct2 = synpdf2*binw*length(meanelsen);      % synthetic counts
plot(ax,xx, synct2, 'k-', 'linew', 1.5);
%%% every meanelsen here has the same weight of 1
mmean = mean(meanelsen);
sigmean = std(meanelsen);
plot(ax,[mmean mmean],ax.YLim,'k--','linew',1.5);    % median of wt dist

h1=histogram(ax,meanlzbn,'binw',binw,'facec',[0 1 1],'facea',0.7);  %,'normalization','pdf'[0 1 1]
[mu1,sigma1,muCI,sigmaCI] = normfit(meanlzbn', 0.01);
synpdf1 = normpdf(xx,mu1,sigma1);
synct1 = synpdf1*binw*length(meanlzbn);      % synthetic counts
plot(ax,xx, synct1, 'b-', 'linew', 1.5);
%%% every meanelsen here has the same weight of 1
mmeanlzb = mean(meanlzbn);
sigmeanlzb = std(meanlzbn);
plot(ax,[mmeanlzb mmeanlzb],ax.YLim,'b--','linew',1.5);    % median of wt dist

text(ax, 0.04,0.93,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax, 0.98,0.62,strcat({'$\mu$= '},sprintf('%.2f',mmean),' km'),'FontSize',13,...
     'unit','normalized','Interpreter','latex','color','k','horizontalalignment','right');
text(ax, 0.98,0.55,strcat({'$\sigma$= '},sprintf('%.2f',sigma2),' km'),'FontSize',13,...
     'unit','normalized','Interpreter','latex','color','k','horizontalalignment','right');
text(ax, 0.02,0.62,strcat({'$\mu$= '},sprintf('%.2f',mmeanlzb),' km'),'FontSize',13,...
     'unit','normalized','Interpreter','latex','color','b','horizontalalignment','left');
text(ax, 0.02,0.55,strcat({'$\sigma$= '},sprintf('%.2f',sigma1),' km'),'FontSize',13,...
     'unit','normalized','Interpreter','latex','color','b','horizontalalignment','left');
% text(ax, 0.96,0.93,'Other RTMs','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1,...
%      'horizontalalignment','right'); 
% text(ax, 0.96,0.93,'Opposing RTMs','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1,...
%      'horizontalalignment','right'); 

xlabel(ax,'Mean LF-HF (km)','fontsize',12)
ylabel(ax,'Count','fontsize',12)

hold(ax,'off');


%%% subplot 4, scatter plot between the mean LF-HF from each FEC and the separation between the two
%%% groups in fig 8a, in the ref. frame of mig. front
ax = f.ax(4);
hold(ax,'on');
xran = [-2 2];
yran = [-2.5 2.5];
% yran = [-3 3];
xlim(ax, xran);
ylim(ax, yran);
xticks(ax, xran(1): 0.5: xran(2));
yticks(ax, yran(1): 0.5: yran(2));
axis(ax,'equal');

patarea = [0 -100;
           -100 -100;
           -100 100;
            0 100;
            0 -100];
patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
plot(ax,[-100 100],[0 0],'--','color',[0.6 0.6 0.6],'linew',1);
plot(ax,[0 0],[-100 100],'--','color',[0.6 0.6 0.6],'linew',1);

% option 1, weighted average vs. direct sep in ref frame of moving front
scatter(ax, meanlzbn, mswnendif,15,[0 1 1],'filled','markeredgecolor','k');
scatter(ax, meanlzb, mswnedif,25,'b','filled','markeredgecolor','k');
[minval, ind] = min(abs(meanlzbn)+abs(mswnendif))
scatter(ax, meanlzbn(ind), mswnendif(ind),40,'r','filled','markeredgecolor','k');

% % option 2, weighted average vs. sep in ref frame of geo coord.
% scatter(ax, meanlzbn, mswnendif2,15,[0 1 1],'filled','markeredgecolor','k');
% scatter(ax, meanlzb, mswnedif2,25,'b','filled','markeredgecolor','k');
% [minval, ind] = min(abs(meanlzbn)+abs(mswnendif2))
% scatter(ax, meanlzbn(ind), mswnendif2(ind),25,'k','filled','markeredgecolor','k');

% % option 3, direct average vs. direct sep in ref frame of moving front
% scatter(ax, mswnenave, mswnendif,15,[0 1 1],'filled','markeredgecolor','k');
% scatter(ax, mswneave, mswnedif,25,'b','filled','markeredgecolor','k');
% [minval, ind] = min(abs(mswnenave)+abs(mswnendif))
% scatter(ax, mswnenave(ind), mswnendif(ind),25,'k','filled','markeredgecolor','k');

% % option 4, direct average vs. sep in ref frame of geo coord. 
% scatter(ax, mswnenave, mswnendif2,15,[0 1 1],'filled','markeredgecolor','k');
% scatter(ax, mswneave, mswnedif2,25,'b','filled','markeredgecolor','k');
% [minval, ind] = min(abs(mswnenave)+abs(mswnendif2))
% scatter(ax, mswnenave(ind), mswnendif2(ind),25,'k','filled','markeredgecolor','k');

% % option 5, weighted average vs. projected sep in ref frame of geo coord.
% scatter(ax, meanlzbn, mswnendif2proj,15,[0 1 1],'filled','markeredgecolor','k');
% scatter(ax, meanlzb, mswnedif2proj,25,'b','filled','markeredgecolor','k');
% [minval, ind] = min(abs(meanlzbn)+abs(mswnendif2proj))
% scatter(ax, meanlzbn(ind), mswnendif2proj(ind),25,'k','filled','markeredgecolor','k');

% % option 6, direct average vs. projected sep in ref frame of geo coord.
% scatter(ax, mswnenave, mswnendif2proj,15,[0 1 1],'filled','markeredgecolor','k');
% scatter(ax, mswneave, mswnedif2proj,25,'b','filled','markeredgecolor','k');
% [minval, ind] = min(abs(mswnenave)+abs(mswnendif2proj))
% scatter(ax, mswnenave(ind), mswnendif2proj(ind),25,'k','filled','markeredgecolor','k');

text(ax, 0.04,0.93,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);

xlabel(ax,'Mean LF-HF (km)','fontsize',12)
ylabel(ax,'Difference in mean LF-HF [WSW-E] (km)','fontsize',12)
hold(ax,'off');

%%%%%%%%%%%%%%%%%%%%%%
% meanelsen(ind)
% meanlzbn(ind)
% indneg = find(meanelsen<0 & meanlzbn <0);
% hold(f.ax(2), 'on');
% scatter(f.ax(2),vecerrlzb(indneg,1),vecerrlzb(indneg,2),6,[0.7 0.7 0.7 ],'filled');
% hold(f.ax(2), 'off');

% figure
% aa = meanelsen(indneg)-meanlzbn(indneg);
% histogram(aa);
% bb = abs(aa);
% [bbsort, indsort] = sort(bb);
% indlook = indsort(1:round(length(indsort)/10));
% indlook = find(bb<=0.05);

% aa(indneg(indlook))
% ind111 = find(meanelsen<0 & meanlzbn <0 & abs(meanelsen-meanlzbn)<=0.05);
% length(ind111)
% [~,ind222] = min(abs(fcazi(ind111)-fcazirep));
% [~,ind222] = min(abs(fcmag(ind111)-fcmagrep));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% find the black dot that makes the max. spreading between offset of se, ne, else groups the
%%%%%% smallest, i.e. consistent.
diftol = 0.4;
ind111 = find(maxdif <= diftol);
length(ind111)
disp('point minimize max. distance');
[minmdif, ind222] = min(maxdif)
meanlzbn(ind222)
meanelsen(ind222)
meanswn(ind222)
meannen(ind222)
errlzbspl(ind222,3)
errlzbspl(ind222,4)

disp('point minimize systematic error')
[maxfcmag, indtmp] = max(fcmag(ind111));
ind333 = ind111(indtmp)
maxdif(ind333)
syserr = fcmagrep-maxfcmag
meanlzbn(ind333)
meanelsen(ind333)
meanswn(ind333)
meannen(ind333)
errlzbspl(ind333,3)
errlzbspl(ind333,4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


hold(f.ax(2), 'on');
scatter(f.ax(2), vecerrlzb(ind,1), vecerrlzb(ind,2),40,'r','filled','markeredgecolor','k');
% scatter(f.ax(2),vecerrlzb(indneg(indlook),1),vecerrlzb(indneg(indlook),2),6,'b','filled');
scatter(f.ax(2),vecerrlzb(ind111,1),vecerrlzb(ind111,2),20,[0.75 0.75 0.75],'filled',...
        'markeredgecolor','k');
scatter(f.ax(2),vecerrlzb(ind222,1),vecerrlzb(ind222,2),40,'w','^','filled',...
        'markeredgecolor','k');
scatter(f.ax(2),vecerrlzb(ind333,1),vecerrlzb(ind333,2),40,'w','filled',...
        'markeredgecolor','k');
hold(f.ax(2), 'off');

hold(f.ax(3), 'on');
scatter(f.ax(3),meanlzbn(ind333),700,40,[0 1 1],'filled',...
        'markeredgecolor','k');
scatter(f.ax(3),meanelsen(ind333),700,40,[0.5 0.5 0.5],'filled',...
        'markeredgecolor','k');
lgd = legend(f.ax(3),[h1,h2],{'Opposing RTMs','Other RTMs'},'fontsize',8,...
       'numcolumns',1,...
       'location','northeast');
hold(f.ax(3), 'off');

hold(f.ax(4), 'on');
scatter(f.ax(4), meanlzbn(ind111), mswnendif(ind111),20,[0.7 0.7 0.7],'filled',...
        'markeredgecolor','k');
scatter(f.ax(4), meanlzbn(ind222), mswnendif(ind222),40,'w','^','filled',...
        'markeredgecolor','k');
scatter(f.ax(4), meanlzbn(ind333), mswnendif(ind333),40,'w','filled','markeredgecolor','k');
hold(f.ax(4), 'off');


% figure
% plot(meanelsen(ind111),'b'); hold on
% plot(meanlzbn(ind111),'r');
% plot(abs(meanelsen(ind111)-meanlzbn(ind111)), 'k');
% 
% figure
% scatter(meanlzbn(ind111), meanelsen(ind111), 'b'); hold on
% plot(-1:0.01:1,-1:0.01:1,'r-');
% scatter(meanlzbn(ind111), meanelsen(ind111)-meanlzbn(ind111),'k');
% scatter(meanlzbn(ind111(ind222)), meanelsen(ind111(ind222)),'b','filled');
% scatter(meanlzbn(ind111(ind222)), meanelsen(ind111(ind222))-meanlzbn(ind111(ind222)),'k','filled');


% keyboard
print(f.fig,'-dpdf',...
    strcat(rstpath,'/filtercor_error_summary_v3_exc3mig.pdf'));














    
    
    
    