% function mig_prop_direct_LZB_addfam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to find the best propagation direction for migration
% different from 'mig_prop_direct_LZB' since it contains a new fam 006
%   
% NOTES:
%   2020/09/08, i am adding a new fam 006, so that i want to check if original
%               migrations are still reasonable
%   2020/09/10, i am adding a fam 001 back again for last try, 13 fam in total now
% 
% Chao Song, chaosong@princeton.edu
% First created date:   2020/09/08
% Last modified date:   2020/09/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
% format short e   % Set the format to 5-digit floating point
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
           ];


relacont = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),loccont(2,1),loccont(2,2));
relacont(:,1) = dx;
relacont(:,2) = dy;


% sort according to day, sec
hftime = sortrows(hftime, [13, 15]);
lftime = sortrows(lftime, [13, 15]);


%%
trange = [
           2004195,4.15e+4,4.60e+4;
           2004196,7.55e+4,7.70e+4;
           2004197,0.03e+4,0.35e+4;
           2004197,0.50e+4,0.75e+4;
           2004197,4.60e+4,4.72e+4;
           2004197,8.388e+4,8.46e+4;
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
           2005256,0.37e+4,0.50e+4;
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

% trange(31,:) =  [2005256       72750       74000];
% trange(40,:) =  [2005259       2500       3450];
trange(40,2) = 0.65*3600;
trange(40,3) = 0.92*3600;

angbestl1 = zeros(length(trange),1);
angbestl2 = zeros(length(trange),1);
angbestl3 = zeros(length(trange),1);
angbestl4 = zeros(length(trange),1);
angbestl5 = zeros(length(trange),1);

xran = [-15 25];
yran = [-20 20];

resprophf = nan(length(trange)+1,200);
resproplf = nan(length(trange)+1,50);

% create fit object with free constraints
fttpfree = fittype( @(a,b,x) a*x+b);

indcheck = [40];
indplt = indcheck;

for i = 1: length(indplt)
% for i = 41: 63
    disp(trange(indplt(i),:));
    indhf = find(hftime(:,13)==trange(indplt(i),1) & hftime(:,15)>=trange(indplt(i),2) & ...
                 hftime(:,15)<=trange(indplt(i),3));
    mighf = hftime(indhf,:);

    indlf = find(lftime(:,13)==trange(indplt(i),1) & lftime(:,15)>=trange(indplt(i),2) & ...
                 lftime(:,15)<=trange(indplt(i),3));
    miglf = lftime(indlf,:);
    
    angle = 0:5:360;
    
    l1normhf = zeros(length(angle),1);
    l2normhf = zeros(length(angle),1);
    slopehf = zeros(length(angle),1);
    ssehf = zeros(length(angle),1);
    rmsehf = zeros(length(angle),1);
    rsquarehf = zeros(length(angle),1);
    
    l1normlf = zeros(length(angle),1);
    l2normlf = zeros(length(angle),1);
    slopelf = zeros(length(angle),1);
    sself = zeros(length(angle),1);
    rmself = zeros(length(angle),1);
    rsquarelf = zeros(length(angle),1);
    
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
        slopehf(iang) = coef(1);
        fitprophf = feval(fitobj,mighfdum(:,15)/3600);
        l1normhf(iang) = sum(abs(mighfdum(:,1)-fitprophf))/(length(mighfdum(:,1)));
        l2normhf(iang) = sum((mighfdum(:,1)-fitprophf).^2)/(length(mighfdum(:,1)));
        ssehf(iang) = gof.sse;
        rmsehf(iang) = gof.rmse;
        rsquarehf(iang) = gof.rsquare;    % r-square, defined as ratio of the sum of squares of the regression and the total sum of squares
    
        %%% propagation trial of lf
        miglfdum = miglf;
        for j = 1: size(miglf,1)
            x0 = miglfdum(j,1);
            y0 = miglfdum(j,2);
            [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),0,0);
            miglfdum(j,1) = newx;
            miglfdum(j,2) = newy;
        end
        % linear robust least square
        [fitobj,gof,~] = fit(miglfdum(:,15)/3600, miglfdum(:,1),fttpfree,'Robust','Bisquare',...
                                  'StartPoint',[0.5 5]);
        % output fit parameters
        coef = coeffvalues(fitobj);
        slopelf(iang) = coef(1);
        fitproplf = feval(fitobj,miglfdum(:,15)/3600);
        l1normlf(iang) = sum(abs(miglfdum(:,1)-fitproplf))/(length(miglfdum(:,1)));
        l2normlf(iang) = sum((miglfdum(:,1)-fitproplf).^2)/(length(miglfdum(:,1)));
        sself(iang) = gof.sse;
        rmself(iang) = gof.rmse;
        rsquarelf(iang) = gof.rsquare;    % r-square, defined as ratio of the sum of squares of the regression and the total sum of squares
        
    end

    %%% best angle estimate from hf
    ind = find(slopehf>0);
    ind1 = find(l1normhf(ind)==min(l1normhf(ind)));
    angbestl1(indplt(i)) = angle(ind(ind1(1)));
    ind2 = find(l2normhf(ind)==min(l2normhf(ind)));
    angbestl2(indplt(i)) = angle(ind(ind2(1)));
    ind3 = find(rmsehf(ind)==min(rmsehf(ind)));     % rmse, Root Mean Squared Error, the smaller, the better
    angbestl3(indplt(i)) = angle(ind(ind3(1)));
    ind4 = find(rsquarehf(ind)==max(rsquarehf(ind)));   % R-square, 0-1, the larger the better
    angbestl4(indplt(i)) = angle(ind(ind4(1)));
    ind5 = find(ssehf(ind)==min(ssehf(ind)));   % sum of square due to error, the smaller, the better
    angbestl5(indplt(i)) = angle(ind(ind5(1)));
    disp([angbestl3(indplt(i)),angbestl4(indplt(i)),angbestl5(indplt(i))])
    angbesthf(indplt(i)) = median([angbestl3(indplt(i)),angbestl4(indplt(i)),angbestl5(indplt(i))]);
    
    %%% best angle estimate from hf
    ind = find(slopelf>0);
    ind1 = find(l1normlf(ind)==min(l1normlf(ind)));
    angbestl1(indplt(i)) = angle(ind(ind1(1)));
    ind2 = find(l2normlf(ind)==min(l2normlf(ind)));
    angbestl2(indplt(i)) = angle(ind(ind2(1)));
    ind3 = find(rmself(ind)==min(rmself(ind)));     % rmse, Root Mean Squared Error, the smaller, the better
    angbestl3(indplt(i)) = angle(ind(ind3(1)));
    ind4 = find(rsquarelf(ind)==max(rsquarelf(ind)));   % R-square, 0-1, the larger the better
    angbestl4(indplt(i)) = angle(ind(ind4(1)));
    ind5 = find(sself(ind)==min(sself(ind)));   % sum of square due to error, the smaller, the better
    angbestl5(indplt(i)) = angle(ind(ind5(1)));
    disp([angbestl3(indplt(i)),angbestl4(indplt(i)),angbestl5(indplt(i))])
    angbestlf(indplt(i)) = median([angbestl3(indplt(i)),angbestl4(indplt(i)),angbestl5(indplt(i))]);
    
    % determine if more inspection is needed
    if abs(angbesthf(indplt(i))-angbestlf(indplt(i))) > 20
        disp('Difference in propagation direction between hf & lf is too large!');
    end
    
end

angbesthf = angbesthf';
angbestlf = angbestlf';

angbest = [angbesthf angbestlf];

angbest2 = [265	245;
215	230;
180	250;
325	255;
235	235;
220	255;
165	245;
270	250;
235	265;
200	230;
75	85;
235	240;
240	245;
250	260;
120	85;
195	255;
90	85;
170	115;
255	270;
115	75;
160	35;
155	30;
345	70;
240	250;
175	215;
215	230;
115	105;
110	85;
250	215;
260	250;
180	270;
250	255;
195	225;
245	260;
165	85;
245	230;
240	255;
90	75;
185	60;
135	280;
210	255;
80	90;
240	255;
245	265;
70	90;
235	250;
70	80; 
            ];
