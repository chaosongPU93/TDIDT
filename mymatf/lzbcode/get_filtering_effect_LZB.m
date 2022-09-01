% function get_filtering_effect_LZB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the function to get the filtering effect under different frequency
% band for LZB trio, should be compatible to all fams
% Using the Chao's rotation parameters to read and rotate the data, we will
% get the direct stacks. They should be already aligned. But when we use
% different frequency band to filter the template, the alignment would
% change, so will the data do. Therefore, we want to get the HF-LF
% difference from the template to correct for data.
%
% Features:
%   1. Based on the all corrections by Chao
%   2. Obtain the templates using corrections at each station per family
%   3. Apply different passbands
%
%  ABANDONED!!!! NOT USED anymore
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/09/03
% Last modified date:   2019/09/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% There are actually 2 ways to do this.
%%% WAY 1:
%%%     Since Make_templatesv2.m already generates the unfiltered and
%%%     unnormed LFE templates saved in files using the rots parameters
%%%     from Chao's new calculations, so the easiest way to get the
%%%     filterring effect is to read these files and alppy different
%%%     passbands.
%%% 
%%% WAY 2:
%%%     Repeat the operations in making templates and stack from scratch to
%%%     check if they are indeed the same. Here use the way 2 only.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% WAY 2: 
%% Initialization
format short e   % Set the format to 5-digit floating point
clear
% close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]
wid=scrsz(3);       % width
hite=scrsz(4);       % height

% set path
workpath = getenv('ALLAN');
%%% choose the data path here 
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/');
rstpath = strcat(datapath, '/LZBtrio');

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

% for ifam = 1: length(nfampool)
% fam = nfampool(ifam, :)
fam = '068'

% get days, 3-sta trio name, bostock's template name, zero crossing index
% of templates
[timoffrot,~] = GetDays4Stack(fam);

% get permanent and polaris station rotation parameters
sft2=0;     % centroid shift of station 2
sft3=0;     % centroid shift of station 3
FLAG = 'LZB';
CATA = 'new';
[PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,sft2,sft3);
reftime = POLROTS(5,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
PERMROTS(:,4) = PERMROTS(:,4)-reftime;    % to make sure that 1st station is 0
POLROTS(:,4) = POLROTS(:,4)-reftime;

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
% before = templen/2-1;
% after = templen/2;
before = templen/8;
after = 7*templen/8-1;
extra = 6*sps;    % use extra more samples than the original window.
stackex = zeros(nsta, templen+ 2*extra);
stackexort = zeros(nsta, templen+ 2*extra);

nday = size(timoffrot, 1)

if nday == 0
    disp("No day with detections found in this family of BOSTOCK's catalog");
else
    
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
        fprintf('%s \n', direc);
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
    
%     %%% plot the 1-step stacked template & write into file
%     %%% figure 1
%     figure('Position',[0.1*wid 0.1*hite 0.3*wid 0.3*hite]);
%     for ista=1: nsta
%         plot(stack(ista, :) + 0.1* ista, 'linewidth', 2); hold on
%         title('direct stack');
%     end
%     box on

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
    
%     %%% plot the normalized 1-step template
%     %%% figure 2
%     figure('Position',[0.1*wid 0.2*hite 0.3*wid 0.3*hite]);
%     for ista = 1: nsta
%         plot(stackno(ista, :) + 2* ista, 'linewidth', 2); hold on
%         title('normalized direct stack');
%     end
%     box on
%     
%     %%% figure 3
%     figure('Position',[0.1*wid 0.3*hite 0.3*wid 0.3*hite]);
%     for ista = 1: nsta
%         plot(stack(ista, mid-5*sps: mid+5*sps) + 0.1* ista, 'linewidth', 2); hold on
%         title('direct stack-short window');
%     end
%     box on

%%% blocks above is basically the same as that in Make_templatesv2.m  %%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% bandpass and compare the cc results under different freq band
    %%%% HF
    lohf = 1.25;
    hihf = 6.5;
    npo = 2;
    npa = 2;
    spshf = 40;
            winsechf = 4:4:60;
%     winsechf = 4;    % window length in sec for CC
    winlenhf = winsechf*spshf;
    nwin = length(winsechf);
    [num, denom] = rat(spshf/40);
    coefmathf = zeros(nsta, nwin);
    lagmathf = zeros(nsta, nwin);
    for iwin = 1: nwin
        for ista = 1: nsta
            ampmax = max(stack(ista, mid-sps/2+1: mid+sps/2));
            ampmin = min(stack(ista, mid-sps/2+1: mid+sps/2));
            norm = max(ampmax, -ampmin);
            dummy = stack(ista, :)/norm;
            dummy1 = Bandpass(dummy, 40, lohf, hihf, npo, npa, 'butter');
            temphf(ista,:) = resample(dummy1, num, denom);
        end
        
%         figure
%         for ista = 1: nsta
%             plot(temphf(ista, mid-5*spshf: mid+5*spshf) + 1* ista, 'linewidth', 2); hold on
%             text(50,1*ista,stas(ista,:));
%         end
%         box on
%         title('Direct stack-short window-bandpass 1.25-6.5 hz');
        
%         figure
        for ista = 1: nsta
            [coefhf, laghf] = xcorr(temphf(ista,mid-winlenhf(iwin)/2+1:mid+winlenhf(iwin)/2), ...
                temphf(5,mid-winlenhf(iwin)/2+1:mid+winlenhf(iwin)/2), 'coeff',spshf/8);    % cc relative to TWKB
            [maxcoef, idx] = max(coefhf);
            lagsamp = laghf(idx);
            coefmathf(ista,iwin) = maxcoef;
            lagmathf(ista,iwin) = lagsamp;
            
%             % if lagsamp > 0, shift this station to left relative to the
%             % reference station to align them, i.e. the arrival of this sta
%             % is late than the reference, so: arr2 == arr1 +lagsamp
%             % and vise versa.
%             plot(temphf(ista, mid+lagsamp-5*spshf+1: mid+lagsamp+5*spshf) + 1* ista, 'linewidth', 2); hold on
%             text(50,1*ista,stas(ista,:));
        end
%         title('Aligned direct stack-short window-bandpass 1.25-6.5 hz');
%         vertical_cursors
        
    end
    
    
    %%%% LF
    lolf = 0.5;
    hilf = 1.25;
    npo = 2;
    npa = 2;
    spslf = 40;
            winseclf = 4:4:60;
%     winseclf = 20;    % window length in sec for CC
    winlenlf = winseclf*spslf;
    nwin = length(winseclf);
    [num, denom] = rat(spslf/40);
    coefmatlf = zeros(nsta, nwin);
    lagmatlf = zeros(nsta, nwin);
    for iwin = 1: nwin
        for ista = 1: nsta
            ampmax = max(stack(ista, mid-sps/2+1: mid+sps/2));
            ampmin = min(stack(ista, mid-sps/2+1: mid+sps/2));
            norm = max(ampmax, -ampmin);
            dummy = stack(ista, :)/norm;
            dummy1 = Bandpass(dummy, 40, lolf, hilf, npo, npa, 'butter');
            templf(ista,:)=resample(dummy1, num, denom);
            
        end
        
%         figure
%         for ista = 1: nsta
%             plot(templf(ista, mid-5*spslf: mid+5*spslf) + 1* ista, 'linewidth', 2); hold on
%             text(50,1*ista,stas(ista,:));
%         end
%         box on
%         title('Direct stack-short window-bandpass 0.5-1.25 hz');
        
%         figure
        for ista = 1: nsta
            [coeflf, laglf] = xcorr(templf(ista,mid-winlenlf(iwin)/2+1:mid+winlenlf(iwin)/2), ...
                templf(5,mid-winlenlf(iwin)/2+1:mid+winlenlf(iwin)/2), 'coeff',spslf/4); % cc relative to TWKB
            [maxcoef, idx] = max(coeflf);
            lagsamp = laglf(idx);
            coefmatlf(ista, iwin) = maxcoef;
            lagmatlf(ista, iwin) = lagsamp;
            
%             % if lagsamp > 0, shift this station to left relative to the
%             % reference station to align them, i.e. the arrival of this sta
%             % is late than the reference, so: arr2 == arr1 +lagsamp
%             % and vise versa.
%             plot(templf(ista, mid+lagsamp-5*spslf+1: mid+lagsamp+5*spslf) + 1* ista, 'linewidth', 2); hold on
%             text(50,1*ista,stas(ista,:));
        end
%         title('Aligned direct stack-short window-bandpass 0.5-1.25 hz');
%         vertical_cursors
        
        
 
        
    end
    
    difflag = lagmathf-lagmatlf;
    
    figure
    for i = 1: nsta
        plot(difflag(i,:),'o-','linewidth',2); hold on
        
    end
    legend(stas(:,:),'location','best');
    title(strcat(fam,',hf-lf,4-60s win,'));
    
    
%     %%% save the filtering effect result
%     % NOTE: all filtering effects are offsets based on 40 sps!!!
%     savemat = [lagmathf lagmatlf lagmathf-lagmatlf];
%     PREFIX = strcat(fam,'.lf',num2str(lolf),'-',num2str(hilf),'.hf',num2str(lohf),'-',num2str(hihf),...
%         '.hfwin',num2str(winsechf),'.lfwin',num2str(winseclf),'.resp-corrected');
%     fid = fopen(strcat(rstpath, '/MAPS/template_filtering_effect',  '_',PREFIX,'_',num2str(sps),'sps'),'w+');
%     fprintf(fid,'%d %d %d \n',savemat');
%     fclose(fid);
    
end     % end if for nday ==0

% end     % end loop for fam


% keyboard
























