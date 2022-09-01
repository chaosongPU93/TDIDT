%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is to merge the different fams or centroids, as long as they 
% are doing detection towards the sam data, especially for PGC trio. 
%
% The reason for doing this is because results in different fams or
% centroids might have duplicates or double counting the same thing. The
% same thing means the same tremor window with same offsets theoratically,
% but sometimes there can be a little numerical error. The only criteria
% used here is still based on the main arrival timing should >= the average
% duration time, otherwise only the detection with highest CC would be
% regarded as an independent one. 
%  
%
% NOTEs:
%   1. remove the double counting first from the original detections of all
%   fams
%   2.
%
% NOTES:
%   NO set 'countmin' to the original detections (set to 1), 2019/08/23 
% 
%
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/09/10
% Last modified date:   2019/09/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
% close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

workpath = getenv('ALLAN');
%%% choose the data path here 
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/PGCtrio');

fam='002';     % family number
freqflag='hf';  % flag to indicate whether to do hf or lf;

% get days, 3-sta trio name, bostock's template name, zero crossing index
% of templates
[timoffrot,stas,bostname,tempoffs] = GetDays(fam,freqflag);


%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
         'LZB'];
POLSTA=['SSIB '           % polaris station names
        'SILB '
        'KLNB '
        'MGCB '
        'TWKB '];
    
%% Important parameters same as that during detections    
nsta=size(stas,1);         %  number of stations
sps=40;     % samples per second

%%% IMPORTANT, NEED to change sometimes
%Basics of the cross-correlation:  Window length, number of windows, filter parameters, etc.
winlensec=4;
winoffsec=1;        % window offset in sec, which is the step of a moving window
winlen=winlensec*sps;      % length in smaples
hi=6.5;
lo=1.25;
npo=2;     % poles, passes of filters
npa=2;
cyclskip = 0;
mshift=29+cyclskip; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
loopoffmax=2.1; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
xcmaxAVEnmin=0.44; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?


mshiftnew=ceil(2.2*mshift);
scrsz=get(0,'ScreenSize');

sft2=-20;     % centroid shift of station 2
sft3=20;     % centroid shift of station 3
lolim = min(-mshift-sft2, -mshift-sft3);
hilim = max(mshift-sft2, mshift-sft3);

matrange = -lolim+hilim+1; 

hitallow = 220;     % max number of hit allowed to iniate the matrix
offmat14 = NaN(matrange, matrange, hitallow);     % +mshiftnew+1 is just for easy visuallization
countmat14 = zeros(matrange, matrange);     % count of hit at the each spot
meanmat14 = NaN(matrange, matrange);  % mean of offset14 at each spot, 
stdmat14 = NaN(matrange, matrange);   % standard deviation of offset14 at each spot
tbgmat14 = zeros(matrange, matrange, hitallow);
tctmat14 = zeros(matrange, matrange, hitallow);
tdatemat14 = zeros(matrange, matrange, hitallow);

offmat15 = NaN(matrange, matrange, hitallow);
countmat15 = zeros(matrange, matrange);
meanmat15 = NaN(matrange, matrange, hitallow);
stdmat15 = NaN(matrange, matrange, hitallow);
tbgmat15 = zeros(matrange, matrange, hitallow);
tctmat15 = zeros(matrange, matrange, hitallow);
tdatemat15 = zeros(matrange, matrange, hitallow);

offmat16 = NaN(matrange, matrange, hitallow);
countmat16 = zeros(matrange, matrange);
meanmat16 = NaN(matrange, matrange, hitallow);
stdmat16 = NaN(matrange, matrange, hitallow);
tbgmat16 = zeros(matrange, matrange, hitallow);
tctmat16 = zeros(matrange, matrange, hitallow);
tdatemat16 = zeros(matrange, matrange, hitallow);

offmat17 = NaN(matrange, matrange, hitallow);
countmat17 = zeros(matrange, matrange);
meanmat17 = NaN(matrange, matrange, hitallow);
stdmat17 = NaN(matrange, matrange, hitallow);
tbgmat17 = zeros(matrange, matrange, hitallow);
tctmat17 = zeros(matrange, matrange, hitallow);
tdatemat17 = zeros(matrange, matrange, hitallow);

% offmat12 = NaN(2*mshift+1, 2*mshift+1, hitallow);     % +mshiftnew+1 is just for easy visuallization
% offmat13 = NaN(2*mshift+1, 2*mshift+1, hitallow);     % +mshiftnew+1 is just for easy visuallization
countmat23 = zeros(matrange, matrange);
tbgmat23 = zeros(matrange, matrange, hitallow);
tctmat23 = zeros(matrange, matrange, hitallow);
tdatemat23 = zeros(matrange, matrange, hitallow);

%%% load the template filtering effect result
lolf = 0.5;
hilf = 1.25;
lohf = 1.25;
hihf = 6.5;
winsechf = winlensec;
winseclf = 20;
PREFIX = strcat(fam,'.lf',num2str(lolf),'-',num2str(hilf),'.hf',num2str(lohf),'-',num2str(hihf),...
                    '.hfwin',num2str(winsechf),'.lfwin',num2str(winseclf),'.resp-corrected');
fname = strcat(rstpath, '/MAPS/template_filtering_effect',  '_',PREFIX,'_','40sps');
filtcor = load(fname);
spsratio = sps/40;
filthf = -filtcor(:,1)*spsratio;   % lf template shift due to filterig effect, but the sign convention is opposite to detection here  

%% START TO LOOP FOR every day
%cycle over each day:
for nd=1:length(timoffrot(:,1))      % num of rows, also num of days
    
    year=timoffrot(nd,1);
    YEAR=int2str(year);
    jday=timoffrot(nd,2);
    if jday <= 9
        JDAY=['00',int2str(jday)];
    elseif jday<= 99
        JDAY=['0',int2str(jday)];
    else
        JDAY=int2str(jday);
    end
    
    % load one centroid or one fam
    IDENTIF=[YEAR,'.',JDAY,'.',fam,'.loff',num2str(loopoffmax),'.ccmin',num2str(xcmaxAVEnmin),'.nponpa',int2str(npo),int2str(npa),'.ms', int2str(mshift)]    
    
    fname1 = strcat(rstpath,'/MAPS/mapallrst',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', '4newsta');   
    fname2 = strcat(rstpath,'/MAPS/mapallrst',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', '3newsta');   
    
    if isfile(fname1)
        rst1= load(fname1);
        %%% additional stations for checking detections
        stasnew=['LZB  '
                 'TWKB '
                 'MGCB '
                 'KLNB '];
        if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
            POLSTA(3,:)='KELB ';
            stasnew(4,:)='KELB ';
        else
            POLSTA(3,:)='KLNB ';  % remember to change it back
            stasnew(4,:)='KLNB ';
        end
        nstanew=size(stasnew,1);
    else
        rst1 = load(fname2);
        stasnew=['LZB  '
                 'TWKB '
                 'MGCB '];
        nstanew=size(stasnew,1);
    end
    
    
    % load another
    fname1 = strcat(rstpath,'/MAPS/mapallrst_newcentroid_',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', '4newsta');
    fname2 = strcat(rstpath,'/MAPS/mapallrst_newcentroid_',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', '3newsta');
    
    if isfile(fname1)
        rst2 = load(fname1);
        %%% additional stations for checking detections
        stasnew=['LZB  '
                 'TWKB '
                 'MGCB '
                 'KLNB '];
        if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
            POLSTA(3,:)='KELB ';
            stasnew(4,:)='KELB ';
        else
            POLSTA(3,:)='KLNB ';  % remember to change it back
            stasnew(4,:)='KLNB ';
        end
        nstanew=size(stasnew,1);
    else
        rst2 = load(fname2);
        stasnew=['LZB  '
                 'TWKB '
                 'MGCB '];
        nstanew=size(stasnew,1);
    end

    rst2(:,2) = rst2(:,2) -sft3;    % convert to a uniform offset system
    rst2(:,3) = rst2(:,3) -sft2;
    
    % construct a new large matrix
    rst = [rst1; rst2];
    rstsort = sortrows(rst, [7,2,3]);   % sort according to main arrival time, could add offset as well
    dtmin = 0.5;      % min time during which only 1 detection is retained
    colnum = [7 4 26 27 28 29];
    indext = 4;    % the extreme case is that every fam contains a duplicate
    [rstnew,indsave,inddis] = RemoveDoubleCounting(rstsort,dtmin,indext,colnum);
    
    %%%%% this part is only for original 3 station trio before adding check
    % load those offsets
    ndSTA12off = rstnew(:, 3);    % 3rd col of mapfile
    ndSTA13off = rstnew(:, 2);    % 2nd col of mapfile
    ndtct23 = rstnew(:, 1);    % timswin == time at the center of each win
    ndtbg23 = rstnew(:, 7);  % begin time in sec of the strongest arrival diapole
    
    date = floor(str2double(strcat(YEAR,JDAY)));
    
    for i = 1: length(ndSTA12off)
    
        countmat23(ndSTA13off(i)-lolim+1, ndSTA12off(i)-lolim+1) = ...
            countmat23(ndSTA13off(i)-lolim+1, ndSTA12off(i)-lolim+1)+ 1;
        % temp is used to indicate the times one spot has been hit
        temp = countmat23(ndSTA13off(i)-lolim+1, ndSTA12off(i)-lolim+1);
        % temp is also the index of the 3rd dimension
        tctmat23(ndSTA13off(i)-lolim+1, ndSTA12off(i)-lolim+1, temp) = ...
            ndtct23(i);
        
        tbgmat23(ndSTA13off(i)-lolim+1, ndSTA12off(i)-lolim+1, temp) = ...
            ndtbg23(i);
        
        tdatemat23(ndSTA13off(i)-lolim+1, ndSTA12off(i)-lolim+1, temp) = ...
            date;
        
    end
    
    %%%%%
    
    
    % get the 3D offset and count matrix
    [temp_count14,temp_off14,temp_tct14,temp_tbg14,temp_tdate14] = ...
        Addcheck_offcount(1,nstanew,date,-lolim,rstnew,countmat14,offmat14,tctmat14,tbgmat14,tdatemat14);
    countmat14 = temp_count14;
    offmat14 = temp_off14;
    tctmat14 = temp_tct14;
    tbgmat14 = temp_tbg14;
    tdatemat14 = temp_tdate14;
    
    [temp_count15,temp_off15,temp_tct15,temp_tbg15,temp_tdate15] = ...
        Addcheck_offcount(2,nstanew,date,-lolim,rstnew,countmat15,offmat15,tctmat15,tbgmat15,tdatemat15);
    countmat15 = temp_count15;
    offmat15 = temp_off15;
    tctmat15 = temp_tct15;
    tbgmat15 = temp_tbg15;
    tdatemat15 = temp_tdate15;
    
    if nstanew > 2
        
        [temp_count16,temp_off16,temp_tct16,temp_tbg16,temp_tdate16] = ...
            Addcheck_offcount(3,nstanew,date,-lolim,rstnew,countmat16,offmat16,tctmat16,tbgmat16,tdatemat16);
        countmat16 = temp_count16;
        offmat16 = temp_off16;
        tctmat16 = temp_tct16;
        tbgmat16 = temp_tbg16;
        tdatemat16 = temp_tdate16;
        
        if nstanew > 3
            
            [temp_count17,temp_off17,temp_tct17,temp_tbg17,temp_tdate17] = ...
                Addcheck_offcount(4,nstanew,date,-lolim,rstnew,countmat17,offmat17,tctmat17,tbgmat17,tdatemat17);
            countmat17 = temp_count17;
            offmat17 = temp_off17;
            tctmat17 = temp_tct17;
            tbgmat17 = temp_tbg17;
            tdatemat17 = temp_tdate17;
                       
        end
        
    end
end

% re-define stasnew and nstanew
stasnew=['LZB  '
         'TWKB '
         'MGCB '
         'KLNB '];   % in 2003, KLNB was named KELB
nstanew=size(stasnew,1);

%%
maxvect = [max(max(countmat14)), max(max(countmat15)), max(max(countmat16)),max(max(countmat17)), max(max(countmat23))];
hitmax = max(maxvect);
if hitmax > hitallow
    fprintf(strcat("ERROR, need to increase the 'hitallow' at least to ", num2str(hitmax), '\n'))
    return;
end    

%%%%% this part is only for original 3 station trio before adding check
% get the index that has offsets, which is also the spot that is activated
% for at least once

% get the subscripts of those activated spots
[ioff13, ioff12] = ind2sub(size(countmat23), find(countmat23 ~= 0));

%%% ADD 1 more constraint to check here.
%%% 1. total hit count in the vicinity (3*3-sample rectangle) is reasonable
countmin = 10;   % min counts allowed in the vicinity
countmin = 1;
newflag = zeros(matrange, matrange);
for i = 1: length(ioff13)
    
    totcount = sum(reshape(countmat23(max(ioff13(i)-1,1): min(ioff13(i)+1,matrange), ...
                                      max(ioff12(i)-1,1): min(ioff12(i)+1,matrange)), [], 1));
    if totcount >= countmin
       
        newflag(ioff13(i), ioff12(i)) = 1;

    end
end    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% still use ioff13, 12, in order to comment out the following line to
% compare the results before and after applying the new constraints easily

[ioff13, ioff12] = ind2sub(size(newflag), find(newflag ~= 0));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the mean & std matrix and vector
STA23count = zeros(length(ioff13), 1);
STA23tbg = [];
STA23tct = [];
STA23tdate = [];

for i = 1: length(ioff13)
%     i=13;
    count = countmat23(ioff13(i), ioff12(i));
    STA23count(i) = count;
    
    tmp1(1:count,1) = ioff12(i);
    tmp1(1:count,2) = ioff13(i);
    
    tmp2 = tbgmat23(ioff13(i), ioff12(i), 1:count);
    tmp2 = reshape(tmp2, count,1);
    STA23tbg = [STA23tbg; tmp2];
    clear tmp2 ;
    
    tmp2 = tctmat23(ioff13(i), ioff12(i), 1:count);
    tmp2 = reshape(tmp2, count,1);
    STA23tct = [STA23tct; tmp2];
    clear tmp2 ;
    
    tmp2 = tdatemat23(ioff13(i), ioff12(i), 1:count);
    tmp2 = reshape(tmp2, count,1);
    tmp3 = [tmp1 tmp2];
    STA23tdate = [STA23tdate; tmp3];
    clear tmp1 tmp2 tmp3;
    
end

% get the offset 12 & 13 vector
STA13off = ioff13 +lolim-1;
STA12off = ioff12 +lolim-1;

%%% filtering correction
STA12off = STA12off - filthf(2);  % since the total offset 12 contains filtering effect, so substract to correct
STA13off = STA13off - filthf(3);

STA13offsec = STA13off/sps;
STA12offsec = STA12off/sps;

STA23tdate(:, 5) = STA23tdate(:, 3);
STA23tdate(:,1) = STA23tdate(:,1)+lolim-1-filthf(2);     % add filtering correction as well
STA23tdate(:,2) = STA23tdate(:,2)+lolim-1-filthf(3);
STA23tdate(:,3:4)= STA23tdate(:,1:2)/sps;

%%% STA23time structure
%%% off12 off13 off12sec off13sec year-day timing-of-strongest-arrival
%%% timing-of-center-of-detection-window
STA23time = [STA23tdate STA23tbg STA23tct];
%%%%%

% get the mean & std matrix and vector of the additional stations
stdtol = 0.1;   % tolerence of std
countmin = 5;   % min counts allowed in the vicinity
[STA14meansec,STA14stdsec,STA134offsec,STA124offsec,STA134off,STA124off,STA14count,STA14tbg,STA14tct,STA14tdate] = ...
    Addcheck_meanstd(1,filthf,countmat14,offmat14,tbgmat14,tctmat14,tdatemat14,-lolim,matrange,sps,stdtol,countmin);

%%% STA14time has the same structure as STA23time
STA14time = [STA14tdate STA14tbg STA14tct];


[STA15meansec,STA15stdsec,STA135offsec,STA125offsec,STA135off,STA125off,STA15count,STA15tbg,STA15tct,STA15tdate] = ...
    Addcheck_meanstd(2,filthf,countmat15,offmat15,tbgmat15,tctmat15,tdatemat15,-lolim,matrange,sps,stdtol,countmin);
STA15time = [STA15tdate STA15tbg STA15tct];

       
[STA16meansec,STA16stdsec,STA136offsec,STA126offsec,STA136off,STA126off,STA16count,STA16tbg,STA16tct,STA16tdate] = ...
    Addcheck_meanstd(3,filthf,countmat16,offmat16,tbgmat16,tctmat16,tdatemat16,-lolim,matrange,sps,stdtol,countmin);
STA16time = [STA16tdate STA16tbg STA16tct];


[STA17meansec,STA17stdsec,STA137offsec,STA127offsec,STA137off,STA127off,STA17count,STA17tbg,STA17tct,STA17tdate] = ...
    Addcheck_meanstd(4,filthf,countmat17,offmat17,tbgmat17,tctmat17,tdatemat17,-lolim,matrange,sps,stdtol,countmin);
STA17time = [STA17tdate STA17tbg STA17tct];


%% PLOT
%%% the final integrated map mean & standard deviation of offset for all days
figure('Position',[scrsz(3)/10 scrsz(4)/10 3*scrsz(3)/5 4.5*scrsz(4)/5]);

nrow = nstanew+1;
ncol = 3;
msize = 3;

% additional sta 1, mean offset
ax1=subplot(nrow, ncol, 1,'align');
scatter(STA124offsec, STA134offsec, msize, STA14meansec, 'filled','s');    %, 'MarkerEdgeColor', 'w')
colormap(ax1,jet)
colorbar
caxis([-1.5, 1.5]);
box on
axis equal
xlabel('STA12 offset (s)')
ylabel('STA13 offset (s)')
axis([-1 1.5 -1.5 1]);
title(strcat(num2str(lo),'-',num2str(hi),'\_mean\_', strtrim(stasnew(1, :))));

% additional sta 1, offset std
ax2=subplot(nrow, ncol, 2,'align');
scatter(STA124offsec, STA134offsec, msize, STA14stdsec, 'filled','s');    %, 'MarkerEdgeColor', 'w')
colormap(ax2,jet)
colorbar
caxis([-1, 1]);
box on
axis equal
xlabel('STA12 offset (s)')
ylabel('STA13 offset (s)')
axis([-1 1.5 -1.5 1]);
title(strcat(num2str(lo),'-',num2str(hi),'\_std\_', strtrim(stasnew(1, :))));    

% additional sta 1, hit count
ax3=subplot(nrow, ncol, 3,'align');
scatter(STA124offsec, STA134offsec, msize, STA14count, 'filled','s');    %, 'MarkerEdgeColor', 'w')
oldcmap = colormap(ax3, hot);
colormap(ax3, flipud(oldcmap) );
colorbar
caxis([0, 5]);
box on
axis equal
xlabel('STA12 offset (s)')
ylabel('STA13 offset (s)')
axis([-1 1.5 -1.5 1]);
title(strcat(num2str(lo),'-',num2str(hi),'\_counts\_', strtrim(stasnew(1, :))));

% additional sta 2, mean offset
ax4=subplot(nrow, ncol, 4,'align');
scatter(STA125offsec, STA135offsec, msize, STA15meansec, 'filled','s');    %, 'MarkerEdgeColor', 'w')
colormap(ax4,jet)
colorbar
caxis([-1.5, 1.5]);
box on
axis equal
xlabel('STA12 offset (s)')
ylabel('STA13 offset (s)')
axis([-1 1.5 -1.5 1]);
title(strcat(num2str(lo),'-',num2str(hi),'\_mean\_', strtrim(stasnew(2, :))));

% additional sta 2, offset std
ax5=subplot(nrow, ncol, 5,'align');
scatter(STA125offsec, STA135offsec, msize, STA15stdsec, 'filled','s');    %, 'MarkerEdgeColor', 'w')
colormap(ax5,jet)
colorbar
caxis([-1, 1]);
box on
axis equal
xlabel('STA12 offset (s)')
ylabel('STA13 offset (s)')
axis([-1 1.5 -1.5 1]);
title(strcat(num2str(lo),'-',num2str(hi),'\_std\_', strtrim(stasnew(2, :))));

% additional sta 2, hit count
ax6=subplot(nrow, ncol, 6,'align');
scatter(STA125offsec, STA135offsec, msize, STA15count, 'filled','s');    %, 'MarkerEdgeColor', 'w')
oldcmap = colormap(ax6, hot);
colormap(ax6, flipud(oldcmap) );
colorbar
caxis([0, 5]);
box on
axis equal
xlabel('STA12 offset (s)')
ylabel('STA13 offset (s)')
axis([-1 1.5 -1.5 1]);
title(strcat(num2str(lo),'-',num2str(hi),'\_counts\_', strtrim(stasnew(2, :))));

% additional sta 3, mean offset
ax7=subplot(nrow, ncol, 7,'align');
scatter(STA126offsec, STA136offsec, msize, STA16meansec, 'filled','s');    %, 'MarkerEdgeColor', 'w')
colormap(ax7,jet)
colorbar
caxis([-1.5, 1.5]);
box on
axis equal
xlabel('STA12 offset (s)')
ylabel('STA13 offset (s)')
axis([-1 1.5 -1.5 1]);
title(strcat(num2str(lo),'-',num2str(hi),'\_mean\_', strtrim(stasnew(3, :))));

% additional sta 3, offset std
ax8=subplot(nrow, ncol, 8,'align');
scatter(STA126offsec, STA136offsec, msize, STA16stdsec, 'filled','s');    %, 'MarkerEdgeColor', 'w')
colormap(ax8,jet)
colorbar
caxis([-1, 1]);
box on
axis equal
xlabel('STA12 offset (s)')
ylabel('STA13 offset (s)')
axis([-1 1.5 -1.5 1]);
title(strcat(num2str(lo),'-',num2str(hi),'\_std\_', strtrim(stasnew(3, :))));

% additional sta 3, hit count
ax9=subplot(nrow, ncol, 9,'align');
scatter(STA126offsec, STA136offsec, msize, STA16count, 'filled','s');    %, 'MarkerEdgeColor', 'w')
oldcmap = colormap(ax9, hot);
colormap(ax9, flipud(oldcmap) );
colorbar
caxis([0, 5]);
box on
axis equal
xlabel('STA12 offset (s)')
ylabel('STA13 offset (s)')
axis([-1 1.5 -1.5 1]);
title(strcat(num2str(lo),'-',num2str(hi),'\_counts\_', strtrim(stasnew(3, :))));

% additional sta 4, mean offset
ax10=subplot(nrow, ncol, 10,'align');
scatter(STA127offsec, STA137offsec, msize, STA17meansec, 'filled','s');    %, 'MarkerEdgeColor', 'w')
colormap(ax10,jet)
colorbar
caxis([-1.5, 1.5]);
box on
axis equal
xlabel('STA12 offset (s)')
ylabel('STA13 offset (s)')
axis([-1 1.5 -1.5 1]);
title(strcat(num2str(lo),'-',num2str(hi),'\_mean\_', strtrim(stasnew(4, :))));

% additional sta 4, offset std
ax11=subplot(nrow, ncol, 11,'align');
scatter(STA127offsec, STA137offsec, msize, STA17stdsec, 'filled','s');    %, 'MarkerEdgeColor', 'w')
colormap(ax11,jet)
colorbar
caxis([-1, 1]);
box on
axis equal
xlabel('STA12 offset (s)')
ylabel('STA13 offset (s)')
axis([-1 1.5 -1.5 1]);
title(strcat(num2str(lo),'-',num2str(hi),'\_std\_', strtrim(stasnew(4, :))));

% additional sta 4, hit count
ax12=subplot(nrow, ncol, 12,'align');
scatter(STA127offsec, STA137offsec, msize, STA17count, 'filled','s');    %, 'MarkerEdgeColor', 'w')
oldcmap = colormap(ax12, hot);
colormap(ax12, flipud(oldcmap) );
colorbar
caxis([0, 5]);
box on
axis equal
xlabel('STA12 offset (s)')
ylabel('STA13 offset (s)')
axis([-1 1.5 -1.5 1]);
title(strcat(num2str(lo),'-',num2str(hi),'\_counts\_', strtrim(stasnew(4, :))));

% original trio hit count
ax13=subplot(nrow, ncol, 13,'align');
scatter(STA12offsec, STA13offsec, msize, STA23count, 'filled','s');    %, 'MarkerEdgeColor', 'w')
oldcmap = colormap(ax13, hot);
colormap(ax13, flipud(oldcmap) );
colorbar
box on
axis equal
caxis([0, 20]);
xlabel('STA12 offset (s)')
ylabel('STA13 offset (s)')
axis([-1 1.5 -1.5 1]);
title(strcat(num2str(lo),'-',num2str(hi),'\_trio-counts\_'));

%% save the image & the matfile
PREFIX = strcat(fam,'.loff',num2str(loopoffmax),'.ccmin',num2str(xcmaxAVEnmin),'.nponpa',int2str(npo),int2str(npa),'.ms', int2str(mshift));

savefile23 = [STA12off STA13off STA12offsec STA13offsec STA23count];
fid = fopen(strcat(rstpath, '/MAPS/maporioffcount_combined_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'newsta'),'w+');
fprintf(fid,'%d %d %.4f %.4f %d \n',savefile23');
fclose(fid);
fid = fopen(strcat(rstpath, '/MAPS/timeori_combined_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'newsta'),'w+');
fprintf(fid,'%d %d %.4f %.4f %d %.4f %.1f \n',STA23time');
fclose(fid);

    
print('-depsc',strcat(rstpath,'/FIGS/',PREFIX,'_','AddStaOffsetCheck_combined', '_',num2str(lo),'-',num2str(hi),'_',strtrim(stasnew(1, :)),'&',strtrim(stasnew(2, :)),'&',strtrim(stasnew(3, :)),'_more.eps'));
%save(strcat(workpath,'/PGCtrio/Matrix/',PREFIX,'_','AddStaOffsetCheck', '_',num2str(lo),'-',num2str(hi),'_',strtrim(stasnew(1, :)),'&',strtrim(stasnew(2, :)),'&',strtrim(stasnew(3, :)),'.mat'));

savefile14 = [STA124off STA134off STA124offsec STA134offsec STA14meansec STA14stdsec STA14count];
fid = fopen(strcat(rstpath, '/MAPS/mapaddcheck_', stasnew(1,:), '_combined_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'newsta'),'w+');
fprintf(fid,'%d %d %.4f %.4f %.4f %.4f %d \n',savefile14');
fclose(fid);
fid = fopen(strcat(rstpath, '/MAPS/timeadd_', stasnew(1,:), '_combined_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'newsta'),'w+');
fprintf(fid,'%d %d %.4f %.4f %d %.4f %.1f \n',STA14time');
fclose(fid);


savefile15 = [STA125off STA135off STA125offsec STA135offsec STA15meansec STA15stdsec STA15count];
fid = fopen(strcat(rstpath, '/MAPS/mapaddcheck_', stasnew(2,:), '_combined_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'newsta'),'w+');
fprintf(fid,'%d %d %.4f %.4f %.4f %.4f %d \n',savefile15');
fclose(fid);
fid = fopen(strcat(rstpath, '/MAPS/timeadd_', stasnew(2,:), '_combined_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'newsta'),'w+');
fprintf(fid,'%d %d %.4f %.4f %d %.4f %.1f \n',STA15time');
fclose(fid);


savefile16 = [STA126off STA136off STA126offsec STA136offsec STA16meansec STA16stdsec STA16count];
fid = fopen(strcat(rstpath, '/MAPS/mapaddcheck_', stasnew(3,:), '_combined_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'newsta'),'w+');
fprintf(fid,'%d %d %.4f %.4f %.4f %.4f %d \n',savefile16');
fclose(fid);
fid = fopen(strcat(rstpath, '/MAPS/timeadd_', stasnew(3,:), '_combined_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'newsta'),'w+');
fprintf(fid,'%d %d %.4f %.4f %d %.4f %.1f \n',STA16time');
fclose(fid);

savefile17 = [STA127off STA137off STA127offsec STA137offsec STA17meansec STA17stdsec STA17count];
fid = fopen(strcat(rstpath, '/MAPS/mapaddcheck_', stasnew(4,:), '_combined_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'newsta'),'w+');
fprintf(fid,'%d %d %.4f %.4f %.4f %.4f %d \n',savefile17');
fclose(fid);
fid = fopen(strcat(rstpath, '/MAPS/timeadd_', stasnew(4,:), '_combined_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'newsta'),'w+');
fprintf(fid,'%d %d %.4f %.4f %d %.4f %.1f \n',STA17time');
fclose(fid);




    