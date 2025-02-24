function [dstack,ccstack,dstackort,ccstackort,dstackvert,ccstackvert] = ...
  mk_bbtemp_PGC(fam,CATA,sps,templensec,ccmethod,ccbp,plflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the function to generate the broadband LFE family template from Bostock's
% catalog, use the shear splitting and polarization parameters
%?scaleseisms?  (theoratically) either from Allan's (Yajun) calculations for PGC trio.
% Saved for future use. 
%
% --Note the difference between 'mk_bbtemp_PGC' and 'mk_bbtemp_LZB' is only until 
%   the first few lines in 'Step 1'. The rest are the same!
% --2022/02/09, Add flag 'CATA' to indicate the lfe catalog type for fam 002; turn
%   the passband of CC the raw templates into inputs
% --USE data removed from station response.
% --Assume now you already have a catalog by hand, from others' LFE
%   detections, or from your own rough detections, this script trys to make a
%   template based on the catalog.
% --2022/11/30, I changed data direc to 'arch????' (necessary); changed 'fact' 
%   at KLNB to be 5e-3 in 2003, and 1.5e-3 after 2003 so that KLNB has a similar
%   amp level as PGC (?); changed 'scaleseisms' to 1 for all stas (?)
%   
%   
%
% INPUT:
%   fam: family string
%   sps: desired sampling rate (original is uniformed to 40)
%   ccmethod: choice number to indicate which criteria to use in CC, 1 with threshold, 2 without
%   plflag: plotting flag whether to plot the templates (1) or not (0)
%
% OUTPUT:
%   dstack: direct stacks
%   ccstack: cross-coorelating corrected stacks
%
% TWO-STEP Precedure:
%   1. stack directly according to the timing of Bostock's catalog with 
%      cutted window. (stack at each station in each LFE family)
%   2. use the stacked template in 1 to cc with every cutted window to
%      make a time shift, then stack again to get a finer template
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/11/02
% Last modified date:   2022/11/30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% default value for easy debugging
defval('fam', '002');
defval('CATA', 'new');
defval('sps', 160);
defval('templensec', 60);
defval('ccmethod', 2);
defval('ccbp',[2 8]);
defval('plflag', 0);

%% Initialization
format short e   % Set the format to 5-digit floating point
% clear
% close all

% set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]
wid=scrsz(3);       % width
hite=scrsz(4);       % height

% set path
workpath = getenv('ALLAN');
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');

% old LFE family pool, from Yajun's selections
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

% the newest families from Bostock are different, so select the closest
% one as new family pool, also selected family have the most LFEs
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

% ifam=1;
% 
% fam = nfampool(ifam, :);
fprintf('Family ID is %s \n', fam);

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
if isequal(fam,'002') && (isequal(CATA,'fixed') || isequal(CATA,'old'))
  [timoffrot,~,bostname,~] = GetDays(fam,freqflag); % from old LFE catalog, return days whorthwhile to analyze 
else
  [timoffrot,~] = GetDays4Stack(fam); % automatically return days whorthwhile to analyze based on number of LFEs in that fam
end

% get rots
sft2=0;     % centroid shift of station 2
sft3=0;     % centroid shift of station 3
FLAG = 'PGC';
% CATA = 'fixed';
[PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,sft2,sft3);
reftime = PERMROTS(1,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
PERMROTS(:,4) = PERMROTS(:,4)-reftime;    % to make sure that 1st station is 0
POLROTS(:,4) = POLROTS(:,4)-reftime;


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

% stas=['PGC  '
%       'SSIB '
%       'SILB '];

% number of used stations
nsta=size(stas,1);         %  number of stations

% convert angles to rads
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;     % convert 2-3 columns to rad
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;

%%%2022/11/30, do we really need this to scale seismograms again after after 'rddata'??
% if isequal(fam,'002') && (isequal(CATA,'fixed') || isequal(CATA,'old')) && nsta==3
%     scaleseisms=[1.0 0.76 0.95];       % scaleseisms scales seismograms
% elseif isequal(fam,'068') && nsta==3
%     scaleseisms=[1.0 0.6 0.6];
% elseif nsta >3
%     scaleseisms=[1.0 0.76 0.95 ones(1,nsta-3)];    
% end
scaleseisms = ones(1,nsta);


%% stack to templates
%%% data parameters
% sps=40;     % samples per second
lenofday = 24*3600;
STAopt=zeros(nsta,sps * lenofday);
STAort=zeros(nsta,sps * lenofday);

%%% desired template parameters
% templensec = 60;
templen = templensec * sps;
% before = 1*templen-1;
% after = 1*templen;
before = 1/2*templen-1;
after = 1/2*templen;
switch fam
    case '001'
        offsec = 22.6;
    case '002'
        offsec = 22.7;
    case '006'
        offsec = 22.6;
    case '010'
        offsec = 22.625;
    case '017'
        offsec = 22.6;
    case '043'
        offsec = 22.75;
    case '047'
        offsec = 22.775;
    case '068'
        offsec = 22.55;
    case '099'
        offsec = 22.625;
    case '125'
        offsec = 22.6;
    case '141'
        offsec = 22.625;
    case '144'
        offsec = 23;
    case '147'
        offsec = 22.65;
    case '243'
        offsec = 24.1;
    case '240'    
        offsec = 25.1;
    otherwise
        offsec = 22.6;
end
if isequal(fam,'002') && (isequal(CATA,'fixed') || isequal(CATA,'old'))
  offsec = 0;
end

before = before - offsec*sps;
after = after + offsec*sps;
extra = 6*sps;    % use extra more samples than the original window.
% stacktmp = zeros(nsta, 1*templen+ 2*extra);
% stacktmport = zeros(nsta, 1*templen+ 2*extra);

nday = size(timoffrot, 1)

if nday == 0
    disp("No day with detections found in this family of BOSTOCK's catalog");
    return
end

% get the all catalog LFEs in that family
if ~(isequal(fam,'002') && (isequal(CATA,'fixed') || isequal(CATA,'old')))  
  bostname = ('/BOSTOCK/total_mag_detect_0000_cull_NEW.txt');
  catalog = load(strcat(workpath, bostname));
  famnum = str2double(fam);
  dateall = catalog(famnum == catalog(:, 1), :);
end

%% STEP 1: In a selected family, loop for every day, read data & catalog, stack
%%% loop for day
nstack = 0;     % nstack should be the same for all stations in the same family, differs w/ fam
nLFE = 0;
fprintf('Part 1: Raw Template \n');

for id = 1: nday
    %     id = 1;
    %Which days of data to read?
    year = timoffrot(id,1);
    YEAR = int2str(year);
    jday = timoffrot(id,2);
      
    if isequal(fam,'002') && (isequal(CATA,'fixed') || isequal(CATA,'old'))
      %Bostock's Detections:
      bostocks=load(bostname(id,:));     % load
      % Speculation: col3-->hour  col4-->sec
      bostsec=3600*(bostocks(:,3)-1)+bostocks(:,4)+(22.83-22.675); %002; 22.675 b.c. 002_ catalogs already "corrected" for PGC start; from tempseps.f
      %22.83 a refinement.
      %bostsec=3600*(bostocks(:,3)-1)+bostocks(:,4)+22.625; %068; TWKB comes in at 905.
      bostsamp=round(bostsec*sps);    % round to nearest integer, round(4.4)=4, round(4.5)=5
    else
      if year == 2003
        date = jday-62+30303;
      elseif year == 2004
        date = jday-196+40714;
      elseif year == 2005
        date = jday-254+50911;
      end
      bostocks = dateall(date == dateall(:, 2), :);
      bostsec = 3600*(bostocks(:,3)-1)+bostocks(:,4);
      bostsamp = round(bostsec*sps);
    end
    
    nlfe1day = size(bostsamp, 1)
    nLFE = nLFE + size(bostsamp, 1);
    
    if jday <= 9
        JDAY=['00',int2str(jday)];
    elseif jday<= 99
        JDAY=['0',int2str(jday)];
    else
        JDAY=int2str(jday);
    end
    MO=day2month(jday,year);     % EXTERNAL function, day2month, get the month of one particular date
    direc=[datapath, '/arch', YEAR,'/',MO,'/'];     % directory name
    fprintf('%s \n', direc);
    prename = [direc, YEAR,'.',JDAY,'.00.00.00.0000.CN'];
    fprintf('Processing date: %s-%s; %d/%d \n',YEAR, JDAY,id,nday);
    
    if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
        POLSTA(3,:)='KELB ';
    else
        POLSTA(3,:)='KLNB ';  % remember to change it back
    end
    
    %in case you mix-use the 'KLNB' or 'KELB' in the station name, this would correct it for you
    if ismember('KLNB ',stas,'rows') || ismember('KELB ',stas,'rows')
      [~,idx1]=ismember('KELB ',stas,'rows');
      [~,idx2]=ismember('KLNB ',stas,'rows');
      if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
        stas(idx1+idx2,:)='KELB ';
      else
        stas(idx1+idx2,:)='KLNB ';  % remember to change it back
      end
    end
    
    fileflag = 1;   % 1 means all files exist, the file status is normal
    
    %%% loop for each station for reading data
    for ista = 1: nsta
        found = 0;
        fname = [];
        
        %%% if station is a permanent one
        [LIA, idx] = ismember(stas(ista, :), PERMSTA, 'rows');
        if LIA
            found = found+ LIA;
            if strcmp(PERMSTA(idx, 1:3),'PGC')
                fact = 1.0e-3;
            elseif strcmp(PERMSTA(idx, 1:3),'LZB')
                fact = 1.8e-3;
            end
            fname = strcat(prename,'.',PERMSTA(idx,1:3),'..BHE.D.SAC');
            if isfile(fname)    % if have the data file
                % 1. this is for data without removing station response
                %   [opt,ort,timeperm]=readperm_1daynofilterv2(datafnm, PERMSTA, PERMROTS, idx, sps, fact);
                % 2. this is for data with no response
                [opt, ort, timeperm] = readperm_nofilterv2(prename, ...
                    PERMSTA, PERMROTS, idx, sps, fact);
                vert = readperm_nofilterv2(prename, ...
                    PERMSTA, PERMROTS, idx, sps, fact, 'Z');
            else
                fileflag = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for station %s in day %s / %s, this day will be omitted. \n',...
                        PERMSTA(idx,1:3), YEAR, JDAY);
                break   % break the entire station loop
            end
            
        end
        
        %%% if station is a polaris one
        [LIA, idx] = ismember(stas(ista, :), POLSTA, 'rows');
        if LIA
            found = found+ LIA;
%             if year == 2003 && jday < 213   % should result from some criteria
            if year==2003 && jday<213 && ~strcmp(POLSTA(idx,1:4),'KELB')    
                fact = 7.5e-3;
            elseif year==2003 && jday<213 && strcmp(POLSTA(idx,1:4),'KELB')      % should result from some criteria
                fact = 5.0e-3;
            else
                fact = 1.5e-3;
            end
            fname = strcat(prename,'.',POLSTA(idx,1:4),'..HHE.D.SAC');
            if isfile(fname)    % if have the data file
                %                 % 1. this is for data without removing station response
                %                 [opt,ort,nzeros]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
                % 2. this is for data with no response
                [opt, ort, ~] = readpola_nofilterv2(prename, ...
                    POLSTA, POLROTS, idx, sps, fact);
                vert = readpola_nofilterv2(prename, ...
                    POLSTA, POLROTS, idx, sps, fact, 'Z');
            else
                fileflag = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for station %s in day %s / %s,this day will be omitted. \n',...
                        POLSTA(idx,1:4), YEAR, JDAY);
                break   % break the entire station loop
            end
            
        end
        
        STAopt(ista, :) = opt/scaleseisms(ista);
        STAort(ista, :) = ort/scaleseisms(ista);
        STAvert(ista, :) = vert/scaleseisms(ista);
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
%             windata = STAopt(:, bostsamp(n)- before- extra: ...
%                 bostsamp(n)+ after+ extra);
%             % remove mean and linear trend
%             tmp = detrend(windata');
%             stacktmp = stacktmp+ tmp';
%             
%             windata2 = STAort(:, bostsamp(n)- before- extra: ...
%                 bostsamp(n)+ after+ extra);
%             tmp2 = detrend(windata2');
%             stacktmport = stacktmport+ tmp2';
            
            nstack = nstack + 1;
            
            % datamat is 3D, for storing all qualified wins for stacking
            % size is (nsta, nstack, templen)
            datamat(:, nstack, :) = STAopt(:, bostsamp(n)- before- extra: ...
                bostsamp(n)+ after+ extra);
              
            dataortmat(:, nstack, :) = STAort(:, bostsamp(n)- before- extra: ...
                bostsamp(n)+ after+ extra);  
              
            datavertmat(:, nstack, :) = STAvert(:, bostsamp(n)- before- extra: ...
                bostsamp(n)+ after+ extra);  
              
        end
    end
    
end     % loop for days

%% find the location of zero-crossing, and adjust the window to be centered at zero-crossing
%this step is replaced by direct correction in time, original step can be found in 'testmk_bbtemp'
%now the indzero is basically the center of current segment, no need to find it any more
% [~,indmax] = max(stacktmp(:, :),[],2);
% indmaxmed = median(indmax);
% [~,indmin] = min(stacktmp(:, :),[],2);
% indminmed = median(indmin);
% indzero = round((indminmed+indmaxmed)/2);
% indzero = size(stacktmp, 2)/2;
% indst = indzero - templen/2 -extra;
% inded = indzero + templen/2 +extra-1;

windata = [];
stackex = zeros(nsta, templen+ 2*extra);
stackortex = zeros(nsta, templen+ 2*extra);
stackvertex = zeros(nsta, templen+ 2*extra);
for i = 1: nstack
%     windata(:,:) = datamat(:, i, indst: inded);
    windata(:,:) = datamat(:, i, :);
    % remove mean and linear trend
    tmp = detrend(windata');
    stackex = stackex+ tmp';
    
    windata(:,:) = dataortmat(:, i, :);
    tmp = detrend(windata');
    stackortex = stackortex+ tmp';
    
    windata(:,:) = datavertmat(:, i, :);
    tmp = detrend(windata');
    stackvertex = stackvertex+ tmp';    
end

%% Averaging & normalization
%%% cut off extra points, stack should be templen long
dstack = stackex(:, 1+extra: end-extra);
dstackort = stackortex(:, 1+extra: end-extra);
dstackvert = stackvertex(:, 1+extra: end-extra);

%%% averaging
dstack = dstack/ nstack;
dstackort = dstackort/nstack;
dstackvert = dstackvert/nstack;

%%% detrend
tmp = detrend(dstack');
dstack = tmp';
tmp2 = detrend(dstackort');
dstackort = tmp2';
tmp3 = detrend(dstackvert');
dstackvert = tmp3';

%%% normalization
mid = templen/2;
ampmax = max(dstack(:, mid-sps+1: mid+sps),[],2);
ampmin = min(dstack(:, mid-sps+1: mid+sps),[],2);   %/2
for ista = 1: nsta
    if ampmax(ista) >= -ampmin(ista)
        fprintf('At station %s positive is bigger \n', strtrim(stas(ista, :)));
    else
        fprintf('At station %s negative is bigger \n', strtrim(stas(ista, :)));
    end
end
norm = max(ampmax, -ampmin);
dstackno = dstack./ norm;


%% plotting if desired
if plflag
    %%% plot the 1-step stacked template
    %%% figure 1
    figure('Position',[0.1*wid 0.1*hite 0.3*wid 0.3*hite]);
    for ista=1: nsta
        plot(dstack(ista, :) + 0.1* ista, 'linewidth', 2); hold on
        title('direct stack');
    end
    box on
    grid on
    
    %%% plot the normalized 1-step template
    %%% figure 2
    figure('Position',[0.1*wid 0.2*hite 0.3*wid 0.3*hite]);
    color=['r','b','k'];
    for ista = 1: nsta
        plot((1:10*sps+1)/sps,dstackno(ista, mid-5*sps: mid+5*sps) + 1* ista, 'linewidth', 1,'color',...
            color(ista)); hold on
        %     plot(stackno(ista, mid-5*sps: mid+5*sps) + 1* ista, 'linewidth', 1,'color',...
        %          color(ista)); hold on
        text(1,1*ista,stas(ista,:),'fontsize',12);
    end
    title('normalized direct stack');
    xlabel('Time (s)');
    % xlabel('Samples (on 40 sps)');
    ylabel('Normalized Amp.');
    box on
    grid on
    
    %%% figure 3
    figure('Position',[0.1*wid 0.3*hite 0.3*wid 0.3*hite]);
    for ista = 1: nsta
        plot(dstack(ista, mid-5*sps: mid+5*sps) + 0.1* ista, 'linewidth', 2); hold on
        title('direct stack-short window');
    end
    box on
    grid on
    
end


%% PART 2: Multi-CC the 1-step template and each window to get a finer template
% ccmethod = 1;   % which method is using, 1 with threshold; 2 without

if ccmethod == 1    % start for CCMETHOD choice
    %%%%%%%%%%%%%%%% WAY 1, stack only the ones that pass the cc and offset thresholds %%%%%%%%%%%%%%%%%
    %%% the previous 'extra' is tO avoid adding zeros to the array due to shifting,
    nstack2 = 0;   % here nstack2 can vary w/ different stations
    % according to template, the duration of dipole seems to be ~0.5 sec (10
    % samples), could even as small as 5 samples
    offsetmax = sps/4;
    % 0.25/ 0.3 could be a reasonable cc coef threshold
    coefmin = 0.25;
    ccstack = zeros(nsta, templen);
    ccstackort = zeros(nsta, templen);
    ccstackvert = zeros(nsta, templen);
    %%% SEE NOTES for proper freq. band choices
%     lo = 2;     % empirical 
%     hi = 8;
    npo = 2;
    npa = 2;
    coefmat = zeros(nsta, nstack);
    lagmat = zeros(nsta, nstack);
    fprintf('Part 2: Precise Template \n');
    % figure
    for ista = 1: nsta    % loop for station
        %     ista = 1;
        dummy(:,:) = datamat(ista,:,:);
%         dummy(:,:) = datamat(ista,:, indst: inded);
        windata = dummy';
        dummy2 = stackex(ista, :)/ nstack;
        template = dummy2';
        for istack = 1: nstack
            fwindata(:,istack) = Bandpass(windata(:,istack),sps,ccbp(1),ccbp(2),npo,npa,'butter');  % also contains detrend
        end
        ftemplate = Bandpass(template,sps,ccbp(1),ccbp(2),npo,npa,'butter');
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
        %%% the left to align with B, (i.e. A-B=4), B is the reference station
        win = fwindata(templen/2+extra-2*sps+1:templen/2+extra+2*sps,:);  % the length is 4*sps
        tem = ftemplate(templen/2+extra-2*sps+1:templen/2+extra+2*sps,:);
        win = detrend(win);
        tem = detrend(tem);
        for istack = 1: nstack
            [coef, lag] = xcorr(win(:,istack), tem, offsetmax, 'coeff');
            [maxcoef, idx] = max(coef);
            lagsamp = lag(idx);
            coefmat(ista, istack) = maxcoef;
            lagmat(ista, istack) = lagsamp;
        end
        
    end   % loop for station
            
    windata = [];
    avecoef = sum(coefmat,1);    
    for istack = 1: nstack
        
        if (abs(lagmat(:, istack)) <= offsetmax) && (avecoef(istack) >= coefmin) ...
                && (coefmat(:,istack) >= 2/3*coefmin)
            windata = datamat(:,istack,:);
%             windata = datamat(:,istack, indst: inded);
            newdata = windata(:,1+lagmat(ista, istack)+extra: ...
                                lagmat(ista, istack)+extra+templen);
            nstack2 = nstack2+ 1;
            winlag(:, nstack2) = lagmat(:, istack);
            tmp = detrend(newdata');
            ccstack = ccstack+ tmp';
            
            windata = dataortmat(:,istack, :);
            newdata = windata(:,1+lagmat(ista, istack)+extra: ...
                                lagmat(ista, istack)+extra+templen);
            tmp = detrend(newdata');
            ccstackort = ccstackort+ tmp';
            
            windata = datavertmat(:,istack, :);
            newdata = windata(:,1+lagmat(ista, istack)+extra: ...
                                lagmat(ista, istack)+extra+templen);
            tmp = detrend(newdata');
            ccstackvert = ccstackvert+ tmp';
            
        end
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% WAY 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ccmethod == 2
    %%%%%%%%%%%%%%%% WAY 2, stack all, only use cc to correct time %%%%%%%%%%%%%%%%%
    %%% the previous 'extra' is tO avoid adding zeros to the array due to shifting,
    nstack2 = nstack;
    % according to template, the duration of dipole seems to be ~0.5 sec (10
    % samples), could even as small as 5 samples
    offsetmax = sps/4;
    ccstack = zeros(nsta, templen);
    ccstackort = zeros(nsta, templen);
    ccstackvert = zeros(nsta, templen);
    %%% SEE NOTES for proper freq. band choices
%     lo = 2;     % empirical 
%     hi = 8;
    npo = 2;
    npa = 2;
    coefmat = zeros(nsta, nstack);
    lagmat = zeros(nsta, nstack);
    fprintf('Part 2: Precise Template \n');
    % figure
    for ista = 1: nsta    % loop for station
        %     ista = 1;
        dummy(:,:) = datamat(ista,:,:);
%         dummy(:,:) = datamat(ista,:, indst: inded);
        windata = dummy';
        dummy2 = stackex(ista, :)/ nstack;
        template = dummy2';
        %             tic
        for istack = 1: nstack
            fwindata(:,istack) = Bandpass(windata(:,istack),sps,ccbp(1),ccbp(2),npo,npa,'butter');  % also contains detrend
        end
        %             toc
        ftemplate = Bandpass(template,sps,ccbp(1),ccbp(2),npo,npa,'butter');
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
        %%% the left to align with B, (i.e. A-B=4), B is the reference station
        win = fwindata(templen/2+extra-2*sps+1:templen/2+extra+2*sps,:);  % the length is 4*sps
        tem = ftemplate(templen/2+extra-2*sps+1:templen/2+extra+2*sps,:);
        win = detrend(win);
        tem = detrend(tem);
        dummyort(:,:) = dataortmat(ista,:, :);
        windataort = dummyort';
        dummyvert(:,:) = datavertmat(ista,:, :);
        windatavert = dummyvert';
        
        f=initfig(5,2.5,1,1);
        ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
        for istack = 1: nstack
            [coef, lag] = xcorr(win(:,istack), tem, offsetmax, 'coeff');
            [maxcoef, idx] = max(coef);
            lagsamp = lag(idx);
            coefmat(ista, istack) = maxcoef;
            lagmat(ista, istack) = lagsamp;
            
            newdata(:) = windata(1+lagsamp+extra: lagsamp+extra+templen,istack);
            winlag(ista, istack) = lagsamp;
            tmp = detrend(newdata');
            ccstack(ista, :) = ccstack(ista, :)+ tmp';
            
            tmp=Bandpass(tmp,sps,ccbp(1),ccbp(2),npo,npa,'butter');
            plot(ax,(1:length(tmp'))/sps,tmp','-','color',[0 0 0  0.15],'linew',1);
            
            newdata(:) = windataort(1+lagsamp+extra: lagsamp+extra+templen,istack);
            tmp = detrend(newdata');
            ccstackort(ista, :) = ccstackort(ista, :)+ tmp';
            
            newdata(:) = windatavert(1+lagsamp+extra: lagsamp+extra+templen,istack);
            tmp = detrend(newdata');
            ccstackvert(ista, :) = ccstackvert(ista, :)+ tmp';
        end
        ccstack(ista, :)=Bandpass(ccstack(ista, :),sps,ccbp(1),ccbp(2),npo,npa,'butter');
        plot(ax,(1:length(tmp'))/sps,detrend(ccstack(ista, :)'/nstack),'-',...
          'color','r','linew',2);
                
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
        
        
    end     % loop for station
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% WAY 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end     % end for CCMETHOD choice

%% Averaging and normalization 
%%% averaging
for ista = 1: nsta
    ccstack(ista, :) = ccstack(ista, :)/ nstack2;
    ccstackort(ista, :) = ccstackort(ista, :)/ nstack2;
    ccstackvert(ista, :) = ccstackvert(ista, :)/ nstack2;
end

%%% detrend
tmp = detrend(ccstack');
ccstack = tmp';
tmp = detrend(ccstackort');
ccstackort = tmp';
tmp = detrend(ccstackvert');
ccstackvert = tmp';

%%% normalization
ampmax = max(ccstack(:, mid-sps/2+1: mid+sps/2),[],2);
ampmin = min(ccstack(:, mid-sps/2+1: mid+sps/2),[],2);
for ista = 1: nsta
    if ampmax(ista) >= -ampmin(ista)
        fprintf('At station %s positive is bigger \n', strtrim(stas(ista, :)));
    else
        fprintf('At station %s negative is bigger \n', strtrim(stas(ista, :)));
    end
end
norm = max(ampmax, -ampmin);
ccstackno = ccstack./ norm;


%% plotting if desired
if plflag
    %%% plot the 2-step stacked template 
    %%% figure 1
    figure('Position',[0.1*wid 0.1*hite 0.3*wid 0.3*hite]);
    for ista=1: nsta
        plot(ccstack(ista, :) + 0.1* ista, 'linewidth', 2); hold on
        title('CC stack');
    end
    box on
    grid on
    
    %%% plot the normalized 2-step template
    %%% figure 2
    figure('Position',[0.1*wid 0.2*hite 0.3*wid 0.3*hite]);
    color=['r','b','k'];
    for ista = 1: nsta
        plot((1:10*sps+1)/sps,ccstackno(ista, mid-5*sps: mid+5*sps) + 1* ista, 'linewidth', 1,...
            'color',color(ista)); hold on
        %     plot(stack2no(ista, mid-5*sps: mid+5*sps) + 1* ista, 'linewidth', 1,'color',...
        %          color(ista)); hold on
        text(1,1*ista,stas(ista,:),'fontsize',12);
    end
    title('normalized CC stack');
    xlabel('Time (s)');
    % xlabel('Samples (on 40 sps)');
    ylabel('Normalized Amp.');
    box on
    grid on
    
    if ccmethod == 1
        %%% figure 3
        figure
        supertit(gca,'max cc coef vs sample shift of all traces');
        for ista = 1: nsta
            subplot(3,1,ista)
            scatter(lagmat(ista, :), coefmat(ista, :)); hold on
            plot([-offsetmax, offsetmax], [coefmin, coefmin], 'r-', 'linewidth', 2);
        end
        box on
        grid on
    end
    
    %%% figure 4
    figure
    supertit(gca,'sample shift of preserved traces passing threshold');
    for ista = 1: nsta
        subplot(3,1,ista)
        scatter(1: nstack2, winlag(ista, 1:nstack2));
    end
    box on
    grid on

end


%% calculate the SNR
for ista = 1: nsta
    s = dstack(ista, mid-0.5*sps+1: mid+0.5*sps);
    Es(ista) = sum(s.^2);
    np1 = dstack(ista, mid-1.5*sps+1: mid-0.5*sps);
    np2 = dstack(ista, mid+0.5*sps+1: mid+1.5*sps);
    Enp1(ista) = sum(np1.^2);
    Enp2(ista) = sum(np2.^2);
    En(ista) = (Enp1(ista)+Enp2(ista))/2;
    snr(ista) = Es(ista)/En(ista);
end

for ista = 1: nsta
    s = ccstack(ista, mid-0.5*sps+1: mid+0.5*sps);
    Es2(ista) = sum(s.^2);
    np1 = ccstack(ista, mid-1.5*sps+1: mid-0.5*sps);
    np2 = ccstack(ista, mid+0.5*sps+1: mid+1.5*sps);
    Enp1(ista) = sum(np1.^2);
    Enp2(ista) = sum(np2.^2);
    En2(ista) = (Enp1(ista)+Enp2(ista))/2;
    snr2(ista) = Es2(ista)/En2(ista);
end

%% plot if desired
if plflag
    %%% plot the 2 kinds of template together
    figure('Position',[scrsz(3)/10 scrsz(4)/10 1*scrsz(3)/2 1*scrsz(4)/2]);    
    for ista = 1: nsta
        subplot(nsta,2,2*ista-1);
        plot(dstack(ista, mid-5*sps: mid+5*sps), 'k', 'linewidth', 1); hold on
        plot(ccstack(ista, mid-5*sps: mid+5*sps), 'r', 'linewidth', 1);
        title(strcat('Stacked template with/without CC at station', {' '}, strtrim(stas(ista, :))));
        legend('Direct stack', 'CC stack');
        xlabel('Samples');
        ylabel('Amplitude');
        box on
        grid on
    end
    
    for ista = 1: nsta
        subplot(nsta,2,2*ista);
        plot(dstackno(ista, mid-5*sps: mid+5*sps), 'k', 'linewidth', 1); hold on
        plot(ccstackno(ista, mid-5*sps: mid+5*sps), 'r', 'linewidth', 1);
        title(strcat('Normalized stacked template with/without CC at station', {' '}, ...
              strtrim(stas(ista, :))));
        legend('Direct stack', 'CC stack');
        xlabel('Samples');
        ylabel('Normalized amplitude');
        box on
        grid on
    end

    print(f.fig,'-depsc2',strcat(temppath,'/',fam,'.broadbandstack.eps'));
    
    close(f.fig);
end

% keyboard




















