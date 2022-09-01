% function get_filtering_effect_PGCv2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the function to get the filtering effect under different frequency
% band for PGC trio
%
% Mainly to compare the difference if use the original seperated Bostock's
% catalog (BOSTOCK/NEW/002-246_2003.062') instead of the catalog in one file
% '/BOSTOCK/total_mag_detect_0000_cull_NEW.txt'
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/07/27
% Last modified date:   2019/08/16
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
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
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

fam = nfampool(ifam, :)

% permanent stations rotation angles, NAMES can be found in Rubin &
% Armbruster 2013, Armbruster et al., 2014

%%%%%% Meanning of each column in PERMROTS/POLROTS %%%%%%%%
% 1. offset in samples between fast/slow direction for SAME station,
%    given at 40 sps
% 2. rotation angle to get fast/slow direction for SAME station
% 3. rotation angle to maximize the energy/particle
%    motion/polarization, for SAME station
% 4. offset in samples between the arrival times at different
%    stations, given at 40 sps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get days, 3-sta trio name, bostock's template name, zero crossing index
% of templates
freqflag='hf';
[timoffrot,~,bostname,~] = GetDays(fam,freqflag);

% get rots
FLAG = 'PGC';
CATA = 'new';
[PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,0,0);

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
    'LZB'];
POLSTA=['SSIB '           % polaris station names
    'SILB '
    'KLNB '
    'MGCB '
    'TWKB '];

stas=['PGC  '
    'SSIB '
    'SILB '
    'LZB  '
    'TWKB '
    'MGCB '
    'KLNB '];

% number of used stations
nsta=size(stas,1);         %  number of stations

% convert angles to rads
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;     % convert 2-3 columns to rad
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;

%%% data parameters
sps=40;     % samples per second
lenofday = 24*3600;
STAopt=zeros(nsta,sps * lenofday);
STAort=zeros(nsta,sps * lenofday);

%%% desired template parameters
templensec = 60;
templen = templensec * sps;
before = templen/2-1;
after = templen/2;
% before = templen/8;
% after = 7*templen/8-1;
extra = 6*sps;    % use extra more samples than the original window.
stackex = zeros(nsta, templen+ 2*extra);
stackexort = zeros(nsta, templen+ 2*extra);


nday = size(timoffrot, 1)
% timoffrot=timoffrot(5:9, :);
nday = size(timoffrot, 1);

if nday == 0
    disp("No day with detections found in this family of BOSTOCK's catalog");
end
%%% loop for day
nstack = 0;
nLFE = 0;
fprintf('Part 1: Raw Template \n');

for id = 1: nday
    %     nd = 1;
    %Bostock's Detections:
    bostocks=load(bostname(id,:));     % load
    % Speculation: col3-->hour  col4-->sec
    bostsec=3600*(bostocks(:,3)-1)+bostocks(:,4)+(22.83-22.675); %002; 22.675 b.c. 002_ catalogs already "corrected" for PGC start; from tempseps.f
    %22.83 a refinement.
    %bostsec=3600*(bostocks(:,3)-1)+bostocks(:,4)+22.625; %068; TWKB comes in at 905.
    bostsamp=round(bostsec*sps);    % round to nearest integer, round(4.4)=4, round(4.5)=5
    %Which days of data to read?
    year = timoffrot(id,1);
    YEAR = int2str(year);
    jday = timoffrot(id,2);
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
    fprintf('%s \n', direc);
    datafnm = [direc, YEAR,'.',JDAY,'.00.00.00.0000.CN'];
    fprintf('Processing date: %s / %s \n',YEAR, JDAY);
    nlfe1day = size(bostsamp, 1)
    
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
                [opt, ort, timeperm] = readperm_nofilterv2(datafnm, ...
                    PERMSTA, PERMROTS, idx, sps, fact);
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
                [opt, ort, timepola] = readpola_nofilterv2(datafnm, ...
                    POLSTA, POLROTS, idx, sps, fact);
            else
                fileflag = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for station %s in day %s / %s. \n',POLSTA(idx,1:4), YEAR, JDAY);
                break   % break the entire station loop
            end
            
        end
        
        STAopt(ista, :) = opt;
        STAort(ista, :) = ort;
    end
    
    if fileflag == 0    % means there are missing files
        fprintf('Day %s / %s will be omitted because of missing files. \n', YEAR, JDAY);
        continue    % continue to the next day
    end
    
    len = size(STAopt,2);
    for n=1: size(bostsamp, 1)
        % set the timing be at the center of the window
        if (bostsamp(n)+ after+ extra <= len) && ...
                (bostsamp(n)- before- extra >= 1)
            
            %%% it is better to use a longer window than templen to make
            %%% sure the shifting in step2 would preserve the frequency
            %%% content as much as possible
            windata = STAopt(:, bostsamp(n)- before- extra: ...
                bostsamp(n)+ after+ extra);
            % remove mean and linear trend
            tmp = detrend(windata');
            stackex = stackex+ tmp';
            
            windata2 = STAort(:, bostsamp(n)- before- extra: ...
                bostsamp(n)+ after+ extra);
            tmp2 = detrend(windata2');
            stackexort = stackexort+ tmp2';
            
            nstack = nstack + 1;
            
            % datamat is 3D, for storing all qualified wins for stacking
            % size is (nsta, nstack, templen)
            datamat(:, nstack, :) = STAopt(:, bostsamp(n)- before- extra: ...
                bostsamp(n)+ after+ extra);
        end
    end
    
end     % loop for days

%% Averaging and and plotting

%%% cut off extra points, stack should be templen long
stack = stackex(:, 1+extra: end-extra);
stackort = stackexort(:, 1+extra: end-extra);

%%% averaging
stack = stack/ nstack;
stackort = stackort/nstack;

%%% detrend
tmp = detrend(stack');
stack = tmp';
tmp2 = detrend(stackort');
stackort = tmp2';

%%% plot the 1-step stacked template & write into file
%%% plot

% %%% figure 1
% figure('Position',[0.1*wid 0.1*hite 0.3*wid 0.3*hite]);
% for ista=1: nsta
%     plot(stack(ista, :) + 0.1* ista, 'linewidth', 2); hold on
%     title('direct stack');
% end
% box on


%% normalization
stackno = zeros(nsta, templen);
stackortno = zeros(nsta, templen);
mid = templen/2;

for ista = 1: nsta
    
    ampmax = max(stack(ista, mid-sps/2+1: mid+sps/2));
    ampmin = min(stack(ista, mid-sps/2+1: mid+sps/2));
    if ampmax >= -ampmin
        fprintf('At station %s positive is bigger \n', strtrim(stas(ista, :)));
    else
        fprintf('At station %s negative is bigger \n', strtrim(stas(ista, :)));
    end
    norm = max(ampmax, -ampmin);
    stackno(ista, :) = stack(ista, :)/ norm;
    stackortno(ista, :) = stackort(ista, :)/ norm;
end

%%% plot the normalized 1-step template
%%% figure 2
figure %('Position',[0.1*wid 0.2*hite 0.3*wid 0.3*hite]);
for ista = 1: nsta
    plot(stackno(ista, mid-10*sps: mid+10*sps) + 1* ista, 'linewidth', 1); hold on
    text(40,1*ista,stas(ista,:),'fontsize',12);
    %title('normalized direct stack');
end
box on

%%% figure 3
figure('Position',[0.1*wid 0.3*hite 0.3*wid 0.3*hite]);
for ista = 1: nsta
    plot(stack(ista, mid-5*sps: mid+5*sps) + 0.1* ista, 'linewidth', 2); hold on
    title('direct stack-short window');
end
box on



%% bandpass and compare the cc results under different freq band

zcf =  @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);

%%%% HF
lohf = 1.25;
hihf = 6.5;
npo = 2;
npa = 2;
spshf = 40;
%             winsechf = 4:4:60;
winsechf = 4;    % window length in sec for CC
winlenhf = winsechf*spshf;
nwin = length(winsechf);
[num, denom] = rat(spshf/40);
coefmathf1 = zeros(nsta, nwin);     % cc coef and lag relative to PGC, sta 1
lagmathf1 = zeros(nsta, nwin);
coefmathf2 = zeros(nsta, nwin);     % cc coef and lag relative to SILB, sta 3
lagmathf2 = zeros(nsta, nwin);
for iwin = 1: nwin
    for ista = 1: nsta
        dummy = Bandpass(stack(ista, :), 40, lohf, hihf, npo, npa, 'butter');
        temphf(ista,:) = resample(dummy, num, denom);
    end
    
    figure
    for ista = 1: nsta
        ampmax = max(temphf(ista, mid-sps/2+1: mid+sps/2));
        ampmin = min(temphf(ista, mid-sps/2+1: mid+sps/2));
        norm = max(ampmax, -ampmin);
        dummy = temphf(ista, :)/norm;
        plot(dummy(mid-10*sps: mid+10*sps) + 1* ista, 'linewidth', 1); hold on
        text(40,1*ista,stas(ista,:),'fontsize',12);
    end
    box on
    %title('Normalized Direct stack-bandpassed 1.25-6.5 hz');
    
    %
    %         figure
    figure
    for ista = 1: nsta
        [coefhf, laghf] = xcorr(temphf(ista,mid-winlenhf(iwin)/2+1:mid+winlenhf(iwin)/2), ...
            temphf(1,mid-winlenhf(iwin)/2+1:mid+winlenhf(iwin)/2), 'coeff'); % cc relative to PGC
        [maxcoef, idx] = max(coefhf);
        lagsamp = laghf(idx);
        coefmathf1(ista,iwin) = maxcoef;
        lagmathf1(ista,iwin) = lagsamp;
        
        
        plot(laghf,coefhf + 1* ista, 'linewidth', 1); hold on
        text(-140,1*ista,stas(ista,:),'fontsize',12);
        vertical_cursors
        
        if maxcoef == coefhf(idx)
            dec1sphf1(ista,iwin) = (2*coefhf(idx)-coefhf(idx-1)-coefhf(idx+1))/2*coefhf(idx);
            dec3sphf1(ista,iwin) = (2*coefhf(idx)-coefhf(idx-3)-coefhf(idx+3))/2*coefhf(idx);
            dec5sphf1(ista,iwin) = (2*coefhf(idx)-coefhf(idx-5)-coefhf(idx+5))/2*coefhf(idx);
            zc = zcf(coefhf);
            ntotal = length(zc);
            npos = length(find(coefhf(zc)>=0));
            nneg = length(find(coefhf(zc)<0));
            if ntotal == npos
                disp('all positive')
            elseif ntotal == nneg
                disp('all negative')
            else
                disp('mixed')
            end   
            zccent = zc-idx;
            izcpos = find(zccent>0,1,'first');
            zcpos = zccent(izcpos)+idx;
            izcneg = find(zccent<0,1,'last');
            zcneg = zccent(izcneg)+idx;
            deczchf1(ista,iwin) = (2*coefhf(idx)-coefhf(zcpos)-coefhf(zcneg))/2*coefhf(idx);
            widzchf1(ista,iwin) = abs(zcpos-zcneg);
        end
        %
        %             % if lagsamp > 0, shift this station to left relative to the
        %             % reference station to align them, i.e. the arrival of this sta
        %             % is late than the reference, so: arr2 == arr1 +lagsamp
        %             % and vise versa.
        %             plot(temphf(ista, mid+lagsamp-5*spshf+1: mid+lagsamp+5*spshf) + 1* ista, 'linewidth', 2); hold on
        %             text(50,1*ista,stas(ista,:));
    end
    %         title('Aligned direct stack-short window-bandpass 1.25-6.5 hz');
    %         vertical_cursors
    
    for ista = 1: nsta
        [coefhf, laghf] = xcorr(temphf(ista,mid-winlenhf(iwin)/2+1:mid+winlenhf(iwin)/2), ...
            temphf(3,mid-winlenhf(iwin)/2+1:mid+winlenhf(iwin)/2), 'coeff'); % cc relative to SILB
        [maxcoef, idx] = max(coefhf);
        lagsamp = laghf(idx);
        coefmathf2(ista,iwin) = maxcoef;
        lagmathf2(ista,iwin) = lagsamp;
        
        if maxcoef == coefhf(idx)
            dec1sphf2(ista,iwin) = (2*coefhf(idx)-coefhf(idx-1)-coefhf(idx+1))/2*coefhf(idx);
            dec3sphf2(ista,iwin) = (2*coefhf(idx)-coefhf(idx-3)-coefhf(idx+3))/2*coefhf(idx);
            dec5sphf2(ista,iwin) = (2*coefhf(idx)-coefhf(idx-5)-coefhf(idx+5))/2*coefhf(idx);
            zc = zcf(coefhf);
            ntotal = length(zc);
            npos = length(find(coefhf(zc)>=0));
            nneg = length(find(coefhf(zc)<0));
            if ntotal == npos
                disp('all positive')
            elseif ntotal == nneg
                disp('all negative')
            else
                disp('mixed')
            end   
            zccent = zc-idx;
            izcpos = find(zccent>0,1,'first');
            zcpos = zccent(izcpos)+idx;
            izcneg = find(zccent<0,1,'last');
            zcneg = zccent(izcneg)+idx;
            deczchf2(ista,iwin) = (2*coefhf(idx)-coefhf(zcpos)-coefhf(zcneg))/2*coefhf(idx);
            widzchf2(ista,iwin) = abs(zcpos-zcneg);
        end
        
    end
    
end


%%%% LF
lolf = 0.5;
hilf = 1.25;
npo = 2;
npa = 2;
spslf = 40;
%             winseclf = 4:4:60;
winseclf = 16;    % window length in sec for CC
winlenlf = winseclf*spslf;
nwin = length(winseclf);
[num, denom] = rat(spslf/40);
coefmatlf1 = zeros(nsta, nwin);
lagmatlf1 = zeros(nsta, nwin);
coefmatlf2 = zeros(nsta, nwin);
lagmatlf2 = zeros(nsta, nwin);
for iwin = 1: nwin
    for ista = 1: nsta
        dummy = Bandpass(stack(ista, :), 40, lolf, hilf, npo, npa, 'butter');
        templf(ista,:)=resample(dummy, num, denom);
        
    end
    
    figure
    for ista = 1: nsta
        ampmax = max(templf(ista, mid-sps+1: mid+sps));
        ampmin = min(templf(ista, mid-sps+1: mid+sps));
        norm = max(ampmax, -ampmin);
        dummy = templf(ista, :)/norm;
        plot(dummy(mid-10*spslf: mid+10*spslf) + 1* ista, 'linewidth', 1); hold on
        text(40,1*ista,stas(ista,:),'fontsize',12);
    end
    box on
%     title('Direct stack-short window-bandpass 0.5-1.25 hz');
    %
    %         figure
    figure
    for ista = 1: nsta
        [coeflf, laglf] = xcorr(templf(ista,mid-winlenlf(iwin)/2+1:mid+winlenlf(iwin)/2), ...
            templf(1,mid-winlenlf(iwin)/2+1:mid+winlenlf(iwin)/2), 'coeff'); % cc relative to PGC
        [maxcoef, idx] = max(coeflf);
        lagsamp = laglf(idx);
        coefmatlf1(ista, iwin) = maxcoef;
        lagmatlf1(ista, iwin) = lagsamp;
        
        plot(laglf,coeflf + 1* ista, 'linewidth', 1); hold on
        text(-580,1*ista,stas(ista,:),'fontsize',12);
        vertical_cursors
        
        if maxcoef == coeflf(idx)
            dec1splf1(ista,iwin) = (2*coeflf(idx)-coeflf(idx-1)-coeflf(idx+1))/2*coeflf(idx);
            dec3splf1(ista,iwin) = (2*coeflf(idx)-coeflf(idx-3)-coeflf(idx+3))/2*coeflf(idx);
            dec5splf1(ista,iwin) = (2*coeflf(idx)-coeflf(idx-5)-coeflf(idx+5))/2*coeflf(idx);
            zc = zcf(coeflf);
            ntotal = length(zc);
            npos = length(find(coeflf(zc)>=0));
            nneg = length(find(coeflf(zc)<0));
            if ntotal == npos
                disp('all positive')
            elseif ntotal == nneg
                disp('all negative')
            else
                disp('mixed')
            end   
            zccent = zc-idx;
            izcpos = find(zccent>0,1,'first');
            zcpos = zccent(izcpos)+idx;
            izcneg = find(zccent<0,1,'last');
            zcneg = zccent(izcneg)+idx;
            deczclf1(ista,iwin) = (2*coeflf(idx)-coeflf(zcpos)-coeflf(zcneg))/2*coeflf(idx);
            widzclf1(ista,iwin) = abs(zcpos-zcneg);
        end
        %
        %             % if lagsamp > 0, shift this station to left relative to the
        %             % reference station to align them, i.e. the arrival of this sta
        %             % is late than the reference, so: arr2 == arr1 +lagsamp
        %             % and vise versa.
        %             plot(templf(ista, mid+lagsamp-5*spslf+1: mid+lagsamp+5*spslf) + 1* ista, 'linewidth', 2); hold on
        %             text(50,1*ista,stas(ista,:));
    end
    %         title('Aligned direct stack-short window-bandpass 0.5-1.25 hz');
    %         vertical_cursors
    
    for ista = 1: nsta
        [coeflf, laglf] = xcorr(templf(ista,mid-winlenlf(iwin)/2+1:mid+winlenlf(iwin)/2), ...
            templf(3,mid-winlenlf(iwin)/2+1:mid+winlenlf(iwin)/2), 'coeff'); % cc relative to SILB
        [maxcoef, idx] = max(coeflf);
        lagsamp = laglf(idx);
        coefmatlf2(ista, iwin) = maxcoef;
        lagmatlf2(ista, iwin) = lagsamp;
        
        if maxcoef == coeflf(idx)
            dec1splf2(ista,iwin) = (2*coeflf(idx)-coeflf(idx-1)-coeflf(idx+1))/2*coeflf(idx);
            dec3splf2(ista,iwin) = (2*coeflf(idx)-coeflf(idx-3)-coeflf(idx+3))/2*coeflf(idx);
            dec5splf2(ista,iwin) = (2*coeflf(idx)-coeflf(idx-5)-coeflf(idx+5))/2*coeflf(idx);
            zc = zcf(coeflf);
            ntotal = length(zc);
            npos = length(find(coeflf(zc)>=0));
            nneg = length(find(coeflf(zc)<0));
            if ntotal == npos
                disp('all positive')
            elseif ntotal == nneg
                disp('all negative')
            else
                disp('mixed')
            end   
            zccent = zc-idx;
            izcpos = find(zccent>0,1,'first');
            zcpos = zccent(izcpos)+idx;
            izcneg = find(zccent<0,1,'last');
            zcneg = zccent(izcneg)+idx;
            deczclf2(ista,iwin) = (2*coeflf(idx)-coeflf(zcpos)-coeflf(zcneg))/2*coeflf(idx);
            widzclf2(ista,iwin) = abs(zcpos-zcneg);
        end
        
    end
    
end

difflag = lagmathf1-lagmatlf1;

figure
for i = 1: 3
    plot(difflag(i,:),'o-','linewidth',2); hold on
    
end
legend(stas(:,:),'location','best');


%%% save the filtering effect result
% NOTE: all filtering effects are offsets based on 40 sps!!!
savemat = [lagmathf1 lagmatlf1 lagmathf1-lagmatlf1];
PREFIX = strcat(fam,'.lf',num2str(lolf),'-',num2str(hilf),'.hf',num2str(lohf),'-',num2str(hihf),...
    '.hfwin',num2str(winsechf),'.lfwin',num2str(winseclf),'.resp-corrected');
fid = fopen(strcat(rstpath, '/MAPS/template_filtering_effect',  '_',PREFIX,'_',num2str(sps),'sps'),'w+');
fprintf(fid,'%d %d %d \n',savemat');
fclose(fid);



% keyboard





















