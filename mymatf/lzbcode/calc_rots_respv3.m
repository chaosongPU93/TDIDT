%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to get the shear-splitting parameter and rotation
% angles from scratch, following the method illustrated in NOTES. First to
% deal with the family 002 as a test, and will add other families in the
% future. version 3
%
% USE the newest version, old versions are preserved for future
%
% USE original data (with station response)
%
% Steps:
%   1.
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/07/24
% Last modified date:   2019/07/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]
wid=scrsz(3);       % width
hite=scrsz(4);       % height

% set path
workpath = getenv('ALLAN');
temppath = strcat(workpath, '/templates/');
if ~exist(temppath, 'dir')
    mkdir(temppath);
end
datapath = workpath;
% datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/PGCtrio');

ofampool = ['002';
    '013';
    '025';
    '028';
    '056';
    '084';
    '099';
    '115';
    '125';
    '147';
    '149'];

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


ifam=1;
close all

fam = nfampool(ifam, :)

%%% get days to be stacked
timoffrot = Readbostock(fam);
ind=[];
timoffrot(ind, :)=[];

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
    'LZB'];
POLSTA=['SSIB '           % polaris station names
    'SILB '
    'KLNB '
    'MGCB '
    'TWKB '];

%  station trio used, 1st station is the main station
%     stas=['LZB  '
%         'PGC  '
%         'SSIB '
%         'SILB '
%         'TWKB '
%         'MGCB '
%         'KLNB '];
stas=['PGC  '
    'SSIB '
    'SILB '
    'LZB  '
    'TWKB '
    'MGCB '
    'KLNB '];

tempoffs=[1211 1297 1231 1221 1207 1186 1207]; %these are the zero crossings

% number of used stations
nsta=size(stas,1);         %  number of stations



%% stack at bostock's detection timings on N and E components

%%% data parameters
sps=40;     % samples per second
lenofday = 24*3600;
STAE=zeros(nsta,sps * lenofday);
STAN=zeros(nsta,sps * lenofday);

%%% desired template parameters
templensec = 60;
templen = templensec * sps;
% before = templen/2-1;
% after = templen/2;
before = templen/8;
after = 7*templen/8-1;
extra = 6*sps;    % use extra more samples than the original window.
stackexE = zeros(nsta, templen+ 2*extra);
stackexN = zeros(nsta, templen+ 2*extra);


nday = size(timoffrot, 1)
% timoffrot=timoffrot(5:9, :);
nday = size(timoffrot, 1);

if nday == 0
    disp("No day with detections found in this family of BOSTOCK's catalog");
    return
end


% get the all catalog LFEs in that family
bostname = ('/BOSTOCK/total_mag_detect_0000_cull_NEW.txt');
catalog = load(strcat(workpath, bostname));
famnum = str2double(fam);
dateall = catalog(famnum == catalog(:, 1), :);

%%% loop for day
nstack = 0;
nLFE = 0;
fprintf('Part 1: Raw Template \n');

for id = 1: nday
    %     nd = 1;
    close all
    year = timoffrot(id,1);
    YEAR = int2str(year);
    jday = timoffrot(id,2);
    if year == 2003
        date = jday-62+30303;
    elseif year == 2004
        date = jday-196+40714;
    elseif year == 2005
        date = jday-254+50911;
    end
    bostocks = dateall(date == dateall(:, 2), :);
    bostsec = 3600*(bostocks(:,3)-1)+bostocks(:,4);
    bostsamp = round(bostsec*40);
    nLFE = nLFE + size(bostsamp, 1);
    
    if jday <= 9
        JDAY=['00',int2str(jday)];
    elseif jday<= 99
        JDAY=['0',int2str(jday)];
    else
        JDAY=int2str(jday);
    end
    %     MO = 'SEP';
    MO=day2month(jday,year);     % EXTERNAL function, day2month, get the month of one particular date
    direc=[datapath, '/', YEAR,'/',MO,'/'];     % directory name
    fprintf('%s \n',direc);
    datafnm = [direc, YEAR,'.',JDAY,'.00.00.00.0000.CN'];
    fprintf('Processing date: %s / %s \n',YEAR, JDAY);
    nlfe1day = size(bostocks, 1)
    
    if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
        POLSTA(3,:)='KELB ';
    else
        POLSTA(3,:)='KLNB ';  % remember to change it back
    end
    
    fileflag = 1;   % 1 means all files exist, the file status is normal
    
    %%% loop for each station for reading data
    for ista = 1: nsta
        found = 0;
        [LIA, idx] = ismember(stas(ista, :), PERMSTA, 'rows');
        
        %%% if station is a permanent one
        if LIA
            found = found+ LIA;
            if strcmp(PERMSTA(idx, 1:3),'PGC')
                fact = 1.0e-3;
            elseif strcmp(PERMSTA(idx, 1:3),'LZB')
                fact = 1.8e-3;
            end
            fname = strcat(datafnm,'.',PERMSTA(idx,1:3),'..BHE.D.SAC');
            if isfile(fname)    % if have the data file
                % 1. this is for data without removing station response
                %   [opt,ort,timeperm]=readperm_1daynofilterv2(datafnm, PERMSTA, PERMROTS, idx, sps, fact);
                % 2. this is for data with no response
                [dataE, dataN, timeperm] = readperms_norots(datafnm, PERMSTA, idx, sps, fact);
            else
                fileflag = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for station %s in day %s %s, this day will be omitted. \n',PERMSTA(idx,1:3), YEAR, JDAY);
                break   % break the entire station loop
            end
            
        end
        
        %%% if station is a polaris one
        [LIA, idx] = ismember(stas(ista, :), POLSTA, 'rows');
        if LIA
            found = found+ LIA;
            if year == 2003 && jday < 213
                fact = 7.5e-3;
            else
                fact = 1.5e-3;
            end
            fname = strcat(datafnm,'.',POLSTA(idx,1:4),'..HHE.D.SAC');
            if isfile(fname)    % if have the data file
                %                 % 1. this is for data without removing station response
                %                 [opt,ort,nzeros]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
                % 2. this is for data with no response
                [dataE, dataN, ~] = readpols_norots(datafnm, POLSTA, idx, sps, fact);
            else
                fileflag = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for station %s in day %s / %s. \n',POLSTA(idx,1:4), YEAR, JDAY);
                break   % break the entire station loop
            end
            
        end
        
        STAE(ista, :) = dataE;
        STAN(ista, :) = dataN;
    end
    
    if fileflag == 0    % means there are missing files
        fprintf('Day %s / %s will be omitted because of missing files. \n', YEAR, JDAY);
        continue    % continue to the next day
    end
    
    len = size(STAE,2);
    for n=1: size(bostsamp, 1)
        % set the timing be at the center of the window
        if (bostsamp(n)+ after+ extra <= len) && ...
                (bostsamp(n)- before- extra >= 1)
            
            %%% it is better to use a longer window than templen to make
            %%% sure the shifting in step2 would preserve the frequency
            %%% content as much as possible
            windata = STAE(:, bostsamp(n)- before- extra: ...
                bostsamp(n)+ after+ extra);
            % remove mean and linear trend
            tmp = detrend(windata');
            stackexE = stackexE+ tmp';
            
            windata2 = STAN(:, bostsamp(n)- before- extra: ...
                bostsamp(n)+ after+ extra);
            tmp2 = detrend(windata2');
            stackexN = stackexN+ tmp2';
            
            nstack = nstack + 1;
            
            % datamat is 3D, for storing all qualified wins for stacking
            % size is (nsta, nstack, templen)
            dmatE(:, nstack, :) = STAE(:, bostsamp(n)- before- extra: ...
                bostsamp(n)+ after+ extra);
            dmatN(:, nstack, :) = STAN(:, bostsamp(n)- before- extra: ...
                bostsamp(n)+ after+ extra);
            
        end
    end
    
end     % loop for days

stackE = stackexE(:, 1+extra: end-extra)/nstack;
stackN = stackexN(:, 1+extra: end-extra)/nstack;



%% rotate to get PERMROTS and POLAROTS parameters

%%%%%%%% Meanning of each column in PERMROTS/POLROTS %%%%%%%%
%%% 1. offset in samples between fast/slow direction, unit in 40 sps, same station
%%% 2. rotation angle to get fast/slow direction, same station
%%% 3. rotation angle to maximize the energy/particle motion/polarization, same station
%%% 4. offset in samples between the arrival times at different stations, unit in 40 sps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

splitrot = 0:5:85;     % shear splitting rotation angle array
polarrot = 1:1:89;     % polarization rotation angle array

splitrotrad=splitrot*pi/180;     % convert angle to rad
polarrotrad=polarrot*pi/180;

ROTS = zeros(nsta, 4);  % rots parameter,


for ista = 1:nsta
    
    % ista=1;
    bef = 40;
    aft = 40;
    seglen = bef+aft;
    
    % STA will be the same length
    STA = stackE(ista,tempoffs(ista)-bef+1: tempoffs(ista)+aft) + 1i * stackN(ista,tempoffs(ista)-bef+1: tempoffs(ista)+aft);
    
    % remove trend and mean
    STA = detrend(STA(:));
    
    for irot = 1: length(splitrotrad)
        STAfastslow = STA * exp(-1i * splitrotrad(irot));
        STAslow = real(STAfastslow);
        STAfast = imag(STAfastslow);
        
        % remove trend and mean
        STAslow = detrend(STAslow(:));
        STAfast = detrend(STAfast(:));
        
        [coef, lag] = xcorr(STAslow, STAfast, 'coeff', sps/4); % 5*sps,
        [maxcoef, idx] = max(coef);
        lagsamp = lag(idx);
        coefmat(irot) = maxcoef;
        lagmat(irot) = lagsamp;
    end
    
    maxlag = lagmat(coefmat==max(coefmat));
    lagslowfast(ista) = maxlag;
    maxsplit = splitrot(coefmat==max(coefmat));
    
    if maxlag > 0   % if staslow is indeed the slower one
        ROTS(ista, 1) = maxlag;
        ROTS(ista, 2) = maxsplit;
    else
        ROTS(ista, 1) = -maxlag;
        %%% make sure after rotation the slow one is at real axis direction,
        %%% although plus 90 would cause the difference in sign, the
        %%% following reverse operation would take care of that
        ROTS(ista, 2) = maxsplit+90;
    end
    
    figure   % plot individual correction results
    
    STAfastslow = STA * exp(-1i * ROTS(ista, 2)*pi/180);
    STAslow = real(STAfastslow);
    STAfast = imag(STAfastslow);
    
    % remove trend and mean
    STAslow = detrend(STAslow(:));
    STAfast = detrend(STAfast(:));
    
    % plot
    plot(STAslow/max(STAslow) + 1, 'linewidth', 2); hold on
    plot(STAfast/max(STAfast) + 2, 'linewidth', 2); hold on
    
    off = round(10*sps/40);
    %%% 1st column of PERMROTS == offset in samples between fast/slow direction
    STAslow(off: seglen-off) = STAslow(off+ ROTS(ista, 1): seglen- off+ ROTS(ista, 1));
    STAsplitcorrected = (STAslow+ 1i * STAfast) * exp(1i* ROTS(ista, 2)*pi/180);
    
    % plot
    plot(STAslow/max(STAslow) + 3, 'linewidth', 2); hold on
    plot(STAfast/max(STAfast) + 4, 'linewidth', 2); hold on
    
    
    for irot = 1: length(polarrotrad)
        STAscrot = STAsplitcorrected * exp(-1i * polarrotrad(irot));
        STAopt = real(STAscrot);
        STAort = imag(STAscrot);
        
        % remove trend and mean
        STAopt = detrend(STAopt(:));
        STAort = detrend(STAort(:));
        
        %             optseg = STAopt(tempoffs(ista)-bef+1: tempoffs(ista)+aft);
        %             ortseg = STAort(tempoffs(ista)-bef+1: tempoffs(ista)+aft);
        engopt = sum(STAopt.^2);
        engort = sum(STAort.^2);
        engratio(irot) = engopt/engort;
    end
    
    maxratio(ista)=max(engratio);
    maxpolar = polarrot(engratio==max(engratio));
    
    %%% usually this maxpolar angle can make sure real axis is optimal
    %%% direction
    ROTS(ista, 3) = maxpolar;
    
    STAscrot = STAsplitcorrected* exp(-1i* maxpolar*pi/180);
    STAopt = real(STAscrot);
    STAort = imag(STAscrot);
    
    % remove trend and mean
    STAopt = detrend(STAopt(:));
    STAort = detrend(STAort(:));
    
    % plot
    plot(STAopt + 5, 'linewidth', 2); hold on
    plot(STAort + 6, 'linewidth', 2); hold on
    
    STAoptmat(ista,:) = STAopt;
    STAortmat(ista,:) = STAort;
    
    %%% this is for saving the result in a longer segment
    STAlong = stackE(ista,:) + 1i * stackN(ista,:);
    STAfastslowlong = STAlong * exp(-1i * ROTS(ista, 2)*pi/180);
    STAslowlong = real(STAfastslowlong);
    STAfastlong = imag(STAfastslowlong);
    STAslowlong(off: templen-off) = STAslowlong(off+ ROTS(ista, 1): templen- off+ ROTS(ista, 1));
    STAsplitcorrectedlong = (STAslowlong+ 1i * STAfastlong) * exp(1i* ROTS(ista, 2)*pi/180);
    STAscrotlong = STAsplitcorrectedlong* exp(-1i* ROTS(ista, 3)*pi/180);
    STAscrotlongmat(ista,:) = STAscrotlong;
    
end

%%% get the offset between different stations
lo = 0.5;
hi = 6.5;
npo = 2;
npa = 2;

for ista = 1: nsta
    ampmax = max(STAoptmat(ista, :));
    ampmin = min(STAoptmat(ista, :));
    norm = max(ampmax, -ampmin);
    dummy = STAoptmat(ista, :)/norm;
    %%% bandpass function contains remove trend and mean
    STAoptnorm(ista,:) = Bandpass(dummy, sps, lo, hi, npo, npa, 'butter');
end

lagmat = zeros(nsta,1);
coefmat = zeros(nsta,1);
for ista = 1: nsta
    [coef, lag] = xcorr(STAoptnorm(ista,:), STAoptnorm(1,:), 'coeff'); % 5*sps,
    [maxcoef, idx] = max(coef);
    lagsamp = lag(idx);
    coefmat(ista) = maxcoef;
    lagmat(ista) = lagsamp;
    ROTS(ista,4) = lagmat(ista)+tempoffs(ista)-tempoffs(1);
end

figure      % plot normalized final optimal results
for ista = 1: nsta
    off=4;
    plot(STAoptnorm(ista, off+lagmat(ista): seglen-off+lagmat(ista)) + 1* ista, 'linewidth', 2); hold on
    text(50,1*ista,stas(ista,:));
end

%%% check if the orthogonal component is indeed very small
figure('Position',[scrsz(3)/10 scrsz(4)/10 4*scrsz(3)/5 2*scrsz(4)/5]);  % plot individual optimal and orthogonal components
for ista = 1: nsta
    off=4;
    subplot(2,round(nsta/2),ista)
    plot(STAoptmat(ista, off+lagmat(ista): seglen-off+lagmat(ista)), 'k', 'linewidth', 2); hold on
    plot(STAortmat(ista, off+lagmat(ista): seglen-off+lagmat(ista)), 'r', 'linewidth', 2); hold on
    text(50,0,stas(ista,:));
end

%%% check the overall one-day data, to see if the final operations are
%%% correct
for ista = 1: nsta
    if ROTS(ista,4) > -1
        STAscrotlongmat(ista, 1:templen-ROTS(ista,4))=STAscrotlongmat(ista, ROTS(ista,4)+1:templen);
        STAscrotlongmat(ista, templen-ROTS(ista,4)+1:templen)=0;
    else
        STAscrotlongmat(ista, -ROTS(ista,4)+1:templen)=STAscrotlongmat(ista, 1:templen+ROTS(ista,4));
        STAscrotlongmat(ista, 1:-ROTS(ista,4))=0;
    end
    STAoptlongmat(ista,:)=real(STAscrotlongmat(ista,:));
    STAortlongmat(ista,:)=imag(STAscrotlongmat(ista,:));
end

figure('Position',[scrsz(3)/10 scrsz(4)/10 4*scrsz(3)/5 2*scrsz(4)/5]);  % plot individual optimal and orthogonal components
for ista = 1: nsta
    subplot(2,round(nsta/2),ista)
    plot(STAoptlongmat(ista, :), 'k', 'linewidth', 2); hold on
    plot(STAortlongmat(ista, :), 'r', 'linewidth', 2); hold on
    text(50,0,stas(ista,:));
end

%%% the offset lagmat11 between the final optimal templates
%%% is expected to zeros
for ista = 1: nsta
    [coef, lag] = xcorr(STAoptlongmat(ista,:), STAoptlongmat(1,:), 'coeff',sps/4); % 5*sps,
    [maxcoef, idx] = max(coef);
    lagsamp = lag(idx);
    coefmat11(ista) = maxcoef;
    lagmat11(ista) = lagsamp;
end
ROTS(:,4) = ROTS(:,4)+ lagmat11';

%     % convert the offset in 40 sps unit to the sps specified.
%     %%% 1st column of PERMROTS == offset in samples between fast and slow components
%     fastslowoff = round(POLROTS(idx,1)*sps/40);  % this offset is based on 40sps data, even for polaris stations
%     %%% 4th column of PERMROTS == offset in samples between the arrival times at different stations
%     STAsoff = round(POLROTS(idx,4)*sps/40);  %input offsets given at 40 sps.
%
%     %Rotate & split-correct
%     STA = STAEd+ 1i* STANd;
%     STAfastslow = STA* exp(-1i* POLROTS(idx,2));   %%% 2th column of PERMROTS, == rotation angle to get fast/slow direction
%     STAslow = real(STAfastslow);
%     STAfast = imag(STAfastslow);
%     len = length(STA);
%     off = round(10*sps/40);
%     %%% 1st column of PERMROTS == offset in samples between fast/slow direction
%     STAslow(off: len-off) = STAslow(off+ fastslowoff: len- off+ fastslowoff);
%     STAsplitcorrected = (STAslow+ 1i* STAfast)* exp(1i* POLROTS(idx,2));
%     STAscrot = STAsplitcorrected* exp(-1i* POLROTS(idx,3));
%
%     %Timeshift
%     if STAsoff > -1
%         STAscrot(1: tracelen- STAsoff) = STAscrot(STAsoff+ 1: tracelen);
%         STAscrot(tracelen- STAsoff+ 1: tracelen) = 0;
%     else
%         STAscrot(-STAsoff+ 1: tracelen) = STAscrot(1: tracelen+ STAsoff);
%         STAscrot(1: -STAsoff) = 0;
%     end
%
%     dataE = real(STAscrot);
%     dataN = imag(STAscrot);

