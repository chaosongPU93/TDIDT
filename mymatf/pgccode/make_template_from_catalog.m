%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is rewritten from Allan's temstackshift.m
% 
% Assume now you already have a catalog by hand, from others' LFE
% detections, or from your own rough detections, this script trys to make a
% template based on the catalog.
% TWO-STEP Precedure:
%   1. stack directly according to the timing of Bostock's catalog with 
%      cutted window. (stack at each station in each LFE family)
%   2. use the stacked template in 1 to cc with every cutted window to
%      make a time shift, then stack again to get a finer template
%
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/02/28
% Last modified date:   2019/03/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e
clear
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');
scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]
wid=scrsz(3);       % width
hite=scrsz(4);       % height
scrat=wid/hite;      % scratch, ratio of width/height
abspath= ('/home/data2/chaosong/matlab/allan/');
rstpath = ('/home/data2/chaosong/matlab/allan/templates/');   
if ~exist(rstpath, 'dir')
  mkdir(rstpath);
end

%% Set rotation parameters, valid in the same family
% WHEN CHANGING FAMILIES CHANGE: (1)timoffrot(dates) (2)Bostnames (3)bostsec
% (4)PERMROTS and POLROTS (5)correction (6)?

% Bostock's LFE family ID 
fam = '002';
fprintf('Family ID is %s \n', fam);

%%%% Different dates, bostnames, PERMROTS and POLROTS and correction  
%%% CASE 002
if isequal(fam, '002')
    % now this is consistent with chaoSTA.m
    timoffrot = [2003 062;
                 2003 063;
                 2004 196;
                 2004 197;
                 2004 198;
                 2004 199;
                 2005 254;
                 2005 255;
                 2005 256];
    % OR, for only one day, use the following 1 line
%     timoffrot= [2005 255;
%                 2005 256];

    %%% bostname is the directory of template from Michael Bostock
    bostname=['BOSTOCK/NEW/002-246_2003.062';   % directory name, path
              'BOSTOCK/NEW/002-246_2003.063';
              'BOSTOCK/NEW/002-246_2004.196';
              'BOSTOCK/NEW/002-246_2004.197';
              'BOSTOCK/NEW/002-246_2004.198';
              'BOSTOCK/NEW/002-246_2004.199';
              'BOSTOCK/NEW/002-246_2005.254';
              'BOSTOCK/NEW/002-246_2005.255';
              'BOSTOCK/NEW/002-246_2005.256'];
    %%% IF only have data of one day, also uncomment the following 1 line
%     bostname=['BOSTOCK/NEW/002-246_2005.255';
%               'BOSTOCK/NEW/002-246_2005.256'];

    %%% permanent stations rotation angles, NAMES can be found in Rubin &
    %%% Armbruster 2013, Armbruster et al., 2014
    
    %%%%%%%% Meanning of each column in PERMROTS/POLROTS %%%%%%%%
    %%% 1. offset in samples between fast/slow direction
    %%% 2. rotation angle to get fast/slow direction
    %%% 3. rotation angle to maximize the energy/particle motion/polarization
    %%% 4. offset in samples between the arrival times at different stations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PERMROTS=[0 90 32 0;  %PGC , Yajun's "0" changed to 90.  Don't think that matters, if no splitting.
              0 90 54 9]; %LZB
    %%% POLARIS stations rotation angles
    POLROTS=[6 85 33 86;  %SSIB from Yajun          
             0 90 39 20;  %SILB
             0 90  7 -4;  %KLNB
             4 70 48 -26; %MGCB
             4 75 38 -5]; %TWKB
    %%% stas = used names of 3-station pair
%     %%%%%%%% for PGC trio %%%%%%%%%%
%     stas=['PGC  '         % NOTICE the blank spaces, stas have 5 columns
%          'SSIB '
%          'SILB '];
%     correction = 22.83-22.675;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%% for LZB trio %%%%%%%%%%
    stas=['TWKB '
          'LZB  '
          'MGCB '];
    correction = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%% CASE 068
elseif isequal(fam,'068')
    timoffrot=[2004 198; 
               2004 199;
               2004 200;
               2004 201;
               2005 256;
               2005 257;
               2005 258;
               2005 259;
               2005 260;
               2005 261];
    bostname=['BOSTOCK/NEW/068_2004.198';       % what is the data structure of these file??
              'BOSTOCK/NEW/068_2004.199';
              'BOSTOCK/NEW/068_2004.200';
              'BOSTOCK/NEW/068_2004.201';
              'BOSTOCK/NEW/068_2005.256';
              'BOSTOCK/NEW/068_2005.257';
              'BOSTOCK/NEW/068_2005.258';
              'BOSTOCK/NEW/068_2005.259';
              'BOSTOCK/NEW/068_2005.260';
              'BOSTOCK/NEW/068_2005.261'];
    PERMROTS=[0  0  0  0;  %PGC
              8 65  6 -8]; %LZB
    POLROTS =[0  0  0  0;  %SSIB 
              0  0  0  0;  %SILB
              2 50 25  0;  %KLNB (erroneous delay w.r.t. TWKB last column)     % w.r.t. with respect to
              5 65 41 -1;  %MGCB
              4 60 14  0]; %TWKB
    stas=['TWKB '
          'LZB  '
          'MGCB '];
    correction = 0;
    
%%% CASE 065
elseif isequal(fam,'065')
    timoffrot=[2004 200;
               2004 201;
               2004 204;
               2005 259;
               2005 260];
    bostname=['BOSTOCK/NEW/065_2004.200';
              'BOSTOCK/NEW/065_2004.201';
              'BOSTOCK/NEW/065_2004.204';
              'BOSTOCK/NEW/065_2005.259';
              'BOSTOCK/NEW/065_2005.260'];
    PERMROTS=[0  0  0  0;  %PGC  4th column here is shift w.r.t. TWKB based on 2004/2005 catalogs.  Better set to zero?
              0  0  0  0];  %LZB (offset is -6)
    POLROTS=[0   0   0  0;  %SSIB
             1 -70  65  0;  %SILB
             3  70  30  0;  %KLNB (delay w.r.t. TWKB is fictitious.
             0   0   0  0;  %MGCB (offset is 1)
             9  65  20  0];  %TWKB
    stas=['TWKB '
          'KLNB '
          'SILB '];
    correction=0;

end
%%%%

%%% names of all available stations
nsta = size(stas, 1);
nday = size(timoffrot, 1);
PERMSTA=['PGC'      % permanent station names
         'LZB'];
POLSTA=['SSIB '     % polaris station names
        'SILB '
        'KLNB '
        'MGCB '
        'TWKB '];

%%% convert angles to rads
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;

%%% data parameters
sps = 40;   % sps is what we want in the end, not intrinsic sps in data 
lenofday = 24*3600;
STAopt=zeros(nsta,sps * lenofday);
STAort=zeros(nsta,sps * lenofday);

%%% desired template parameters
templensec = 60;
templen = templensec * sps;
before = templen/2-1;
after = templen/2;
extra = sps;    % use extra more samples than the original window.
stackex = zeros(nsta, templen+ 2*extra);
stackexort = zeros(nsta, templen+ 2*extra);


%% STEP 1: In a selected family, loop for every day, read data & catalog, stack 
%%% loop for day
nstack = 0;
nLFE = 0;
fprintf('Part 1: Raw Template \n');
for id = 1: nday
%     id = 1;
    bostock = load(strcat(abspath, bostname(id,:)));
    bostsec = 3600*(bostock(:, 3)-1)+ bostock(:, 4)+ correction;
    bostsamp = round(bostsec*sps);
    nLFE = nLFE + size(bostsamp, 1);
    year=timoffrot(id,1);
    YEAR=int2str(year);
    jday=timoffrot(id,2);
    if jday <= 9
        JDAY=['00',int2str(jday)];
    elseif jday<= 99
        JDAY=['0',int2str(jday)];
    else
        JDAY=int2str(jday);
    end
    MO=day2month(jday,year);
    direc=[abspath, YEAR,'/',MO,'/'];
    datafnm = [direc, YEAR,'.',JDAY,'.00.00.00.0000.CN'];
    fprintf('Processing date: %s / %s \n',YEAR, JDAY);
    
    %%% loop for each station for reading data
    for ista = 1: nsta
        found = 0;
        [LIA, idx] = ismember(stas(ista, :), PERMSTA, 'rows');
        
        %%% if station is a permanent one 
        if LIA
            found = found+ LIA;
            if strcmp(PERMSTA(idx, 1:3),'PGC')
                fact = 1.6e-4;
            elseif strcmp(PERMSTA(idx, 1:3),'LZB')
                fact = 4.e-3;
            end
            [opt, ort, timeperm] = readperm_1daynofilter(datafnm, ...
                PERMSTA, PERMROTS, idx, sps, fact);
        end
        
        %%% if station is a polaris one
        [LIA, idx] = ismember(stas(ista, :), POLSTA, 'rows');
        if LIA
            found = found+ LIA;
            if year == 2003 && jday < 213
                fact = 20.0e-3;
            else
                fact = 4.0e-3; 
            end
            [opt, ort, timepola] = readpola_1daynofilter(datafnm, ...
                POLSTA, POLROTS, idx, sps, fact);
        end
        
        STAopt(ista, :) = opt;
        STAort(ista, :) = ort;
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

%%% Averaging and normalization 
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

%%% normalization
stackno = zeros(nsta, templen);
stackortno = zeros(nsta, templen);
for ista = 1: nsta
    ampmax = max(stack(ista, templen/2-sps+1:templen/2+sps));
    ampmin = min(stack(ista, templen/2-sps+1:templen/2+sps));
    if ampmax >= -ampmin
        fprintf('At station %s positive is bigger \n', strtrim(stas(ista, :)));
    else
        fprintf('At station %s negative is bigger \n', strtrim(stas(ista, :)));
    end
    norm = max(ampmax, -ampmin);
    stackno(ista, :) = stack(ista, :)/ norm;
    stackortno(ista, :) = stackort(ista, :)/ norm; 
end

%%% plot the 1-step stacked template & write into file
%%% plot
figure('Position',[0.1*wid 0.1*hite 0.3*wid 0.3*hite]);
for ista = 1: nsta
    plot(stackno(ista, templen/2-2*sps: templen/2+2*sps) + 2* ista, 'linewidth', 2); hold on
    title('normalized direct stack-4s');
end
box on

figure('Position',[0.1*wid 0.2*hite 0.3*wid 0.3*hite]);
for ista = 1: nsta
    plot(stackno(ista, :) + 2* ista, 'linewidth', 2); hold on
    title('normalized direct stack');
end
box on

figure('Position',[0.1*wid 0.3*hite 0.3*wid 0.3*hite]);
for ista=1: nsta
    plot(stack(ista, :) + 0.1* ista, 'linewidth', 2); hold on
    title('direct stack');
end
box on

%%% write into files
for ista = 1: nsta
    fid = fopen(strcat(rstpath, fam, '_', strtrim(stas(ista, :)), '_', ...
        num2str(sps), 'sps_', num2str(templensec), 'sec_', num2str(nstack), ...
        'DirectStacks_', 'opt_Nofilter_Nonorm.txt'), 'w+');
    fprintf(fid, '%f \n', stack(ista, :)');
    fclose(fid);
end

for ista = 1: nsta
    fid = fopen(strcat(rstpath, fam, '_', strtrim(stas(ista, :)), '_', ...
        num2str(sps), 'sps_', num2str(templensec), 'sec_', num2str(nstack), ...
        'DirectStacks_', 'opt_Nofilter_Norm.txt'), 'w+');
    fprintf(fid, '%f \n', stackno(ista, :)');
    fclose(fid);
end


%% PART 2: Multi-CC the 1-step template and each window to get a finer template
%%% TO avoid adding zeros to the array due to shifting, 

nstack2 = zeros(nsta, 1);
% according to template, the duration of dipole seems to be ~0.5 sec (10
% samples), could even as small as 5 samples
offsetmax = sps/4;  
% 0.25/ 0.3 could be a reasonable cc coef threshold
coefmin = 0.25;    
stack2 = zeros(nsta, templen);
%%% SEE NOTES for proper freq. band choices
lo = 2;
hi = 8;
npo = 2;
npa = 2;
coefmat = zeros(nsta, nstack);
lagmat = zeros(nsta, nstack);
fprintf('Part 2: Precise Template \n');
% figure
for ista = 1: nsta    % loop for station
%     ista = 1; 
    for istack = 1: nstack
%         istack = 1;
        windata = zeros(1, templen+ 2*extra);
        windata(:) = datamat(ista, istack, :);
        template = stackex(ista, :)/ nstack;
        
        fwindata = Bandpass(windata, sps, lo, hi, npo, npa, 'butter');
        ftemplate = Bandpass(template, sps, lo, hi, npo, npa, 'butter');
        %%% to SEE how bad is our windowed data before/after filtering
%         subplot(4,1,1)
%         plot(windata(templen/2-4*sps+1:templen/2+4*sps), 'b');
%         subplot(4,1,2)
%         plot(fwindata(templen/2-4*sps+1:templen/2+4*sps), 'r');
%         subplot(4,1,3)
%         plot(template(templen/2-4*sps+1:templen/2+4*sps), 'b');
%         subplot(4,1,4)
%         plot(ftemplate(templen/2-4*sps+1:templen/2+4*sps), 'r');
%         pause(2)
        %%%
        
        %%% e.g. if lag = argmax(xcorr(A, B)) = 4, A needs to shift 4 to
        %%% the left to align with B, (i.e. A-B=4)
        win = fwindata(templen/2+extra-2*sps+1:templen/2+extra+2*sps);
        tem = ftemplate(templen/2+extra-2*sps+1:templen/2+extra+2*sps);
        [coef, lag] = xcorr(win, tem, offsetmax, 'coeff');
        [maxcoef, idx] = max(coef);
        lagsamp = lag(idx);
        coefmat(ista, istack) = maxcoef;
        lagmat(ista, istack) = lagsamp;
        
        if (abs(lagsamp) <= offsetmax)  && (maxcoef >= coefmin)
            newdata(:) = windata(1+lagsamp+extra: lagsamp+extra+templen);
            %%% to TEST if they are indeed aligned
%             subplot(2,1,1)
%             plot(fwindata(templen/2+extra+lagsamp-2*sps+1:templen/2+extra+lagsamp+2*sps), 'b');
%             subplot(2,1,2)
%             plot(ftemplate(templen/2+extra-2*sps+1:templen/2+extra+2*sps), 'r');
%             pause(1)
            %%%%
            
            nstack2(ista) = nstack2(ista)+ 1;
            winlag(ista, nstack2(ista)) = lagsamp;
            tmp = detrend(newdata');
            stack2(ista, :) = stack2(ista, :)+ tmp';        
        end
        
        %%% This was the test for long/short cc comparison ,now in a
        %%% separate test file
%         wincut = detrend(windata(templen+extra-offsetmax+1: templen+extra+offsetmax));
%         temcut = detrend(template(templen+extra-offsetmax+1: templen+extra+offsetmax));
%         [coef, lag] = xcorr(wincut, temcut, 'coeff');
%         [maxcoef, idx] = max(coef);
%         lagsamp = lag(idx);
%         
%         if (abs(lagsamp) < offsetmax)  && (maxcoef >= coefmin)
%             newdata(:) = windata(1+lagsamp+extra: lagsamp+extra+templen);
%             nstack2(ista) = nstack2(ista)+ 1;
%             winlag(ista, nstack2(ista)) = lagsamp;
%             tmp = detrend(newdata');
%             stack2(ista, :) = stack2(ista, :)+ tmp';        
%         end
        %%%
        
    end
    
end     % loop for station

%%% Averaging and normalization 
%%% averaging
for ista = 1: nsta
    stack2(ista, :) = stack2(ista, :)/ nstack2(ista);
end

%%% detrend
tmp = detrend(stack2');
stack2 = tmp';

%%% normalization
stack2no = zeros(nsta, templen);
for ista = 1: nsta
    ampmax = max(stack2(ista, templen/2-2*sps+1:templen/2+2*sps));
    ampmin = min(stack2(ista, templen/2-2*sps+1:templen/2+2*sps));
    if ampmax >= -ampmin
        fprintf('At station %s positive is bigger \n', strtrim(stas(ista, :)));
    else
        fprintf('At station %s negative is bigger \n', strtrim(stas(ista, :)));
    end
    norm = max(ampmax, -ampmin);
    stack2no(ista, :) = stack2(ista, :)/ norm;
end

%%% plot the 2-step stacked template & write into file
%%% plot
figure('Position',[0.3*wid 0.1*hite 0.3*wid 0.3*hite]);
for ista = 1: nsta
    plot(stack2no(ista, templen/2-2*sps: templen/2+2*sps) + 2* ista, 'linewidth', 2); hold on
    title('normalized CC stack');
end
box on

figure('Position',[0.3*wid 0.2*hite 0.3*wid 0.3*hite]);
for ista=1: nsta
    plot(stack2(ista, :) + 0.1* ista, 'linewidth', 2); hold on
    title('CC stack');
end
box on


figure
title('max cc coef vs sample shift');
for ista = 1: nsta
    subplot(3,1,ista)
    scatter(lagmat(ista, :), coefmat(ista, :)); hold on
    plot([-offsetmax, offsetmax], [coefmin, coefmin], 'r-', 'linewidth', 2);
end
box on

figure
title('remaining sample shift');
for ista = 1: nsta
    subplot(3,1,ista)
    scatter(1: nstack2(ista), winlag(ista, 1:nstack2(ista)));
end
box on



%%% write into files
% num2str(nstack2(ista)), 
for ista = 1: nsta
    fid = fopen(strcat(rstpath, fam, '_', strtrim(stas(ista, :)), '_', ...
        num2str(sps), 'sps_', num2str(templensec), 'sec_', ...
        'CCStacks_', 'opt_Nofilter_Nonorm.txt'), ...
        'w+');
    fprintf(fid, '%f \n', stack2(ista, :)');
    fclose(fid);
end

% num2str(nstack2(ista)), 
for ista = 1: nsta
    fid = fopen(strcat(rstpath, fam, '_', strtrim(stas(ista, :)), '_', ...
        num2str(sps), 'sps_', num2str(templensec), 'sec_', ...
        'CCStacks_', 'opt_Nofilter_Norm.txt'), ...
        'w+');
    fprintf(fid, '%f \n', stack2no(ista, :)');
    fclose(fid);
end

%%% calculate the SNR
for ista = 1: nsta
    zoom = stack(ista, templen/2-2*sps: templen/2+2*sps);
    signalE1(ista) = sum(zoom(90-1*sps+1: 90+1*sps).* ...
        zoom(90-1*sps+1: 90+1*sps));
    noiseE1(ista) = sum(zoom(1: 90-1*sps).* ...
        zoom(1: 90-1*sps))+ ...
        sum(zoom(90+1*sps+1: end).* ...
        zoom(90+1*sps+1: end));
    snr1(ista) = signalE1(ista)/noiseE1(ista);
end

for ista = 1: nsta
    zoom2 = stack2(ista, templen/2-2*sps: templen/2+2*sps);
    signalE2(ista) = sum(zoom2(90-1*sps+1: 90+1*sps).* ...
        zoom2(90-1*sps+1: 90+1*sps));
    noiseE2(ista) = sum(zoom2(1: 90-1*sps).* ...
        zoom2(1: 90-1*sps))+ ...
        sum(zoom2(90+1*sps+1: end).* ...
        zoom2(90+1*sps+1: end));
    snr2(ista) = signalE2(ista)/noiseE2(ista);
end

%%% plot the 2 kinds of template together
figure
for ista = 1: nsta
    subplot(3,1,ista);
    plot(stackno(ista, templen/2-2*sps: templen/2+2*sps), 'b', 'linewidth', 2); hold on
    plot(stack2no(ista, templen/2-2*sps: templen/2+2*sps), 'r', 'linewidth', 2); 
    title(strcat('Normalized stacked template with/without CC at Station', ' ', strtrim(stas(ista, :))));
    legend('Direct stack', 'CC stack');
    xlabel('Sample number');
    ylabel('Normalized amplitude');
end

figure
for ista = 1: nsta
    subplot(3,1,ista);
    plot(stack(ista, templen/2-2*sps: templen/2+2*sps), 'b', 'linewidth', 2); hold on
    plot(stack2(ista, templen/2-2*sps: templen/2+2*sps), 'r', 'linewidth', 2); 
    title(strcat('Stacked template with/without CC at Station', ' ', strtrim(stas(ista, :))));
    legend('Direct stack', 'CC stack');
    xlabel('Sample number');
    ylabel('Amplitude');
end

















