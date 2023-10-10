%allcomptemp002.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the script to obtain the vertical (Z) and horizontal (N and E) 
% templates of fam 002, to see if we can see the P-wave arrival other than
% S-wave. Many segments of codes are similar to 'calc_rots_noresp.m' and 
% 'mk_bbtemp_pgc.m'.
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2023/07/18
% Last modified date:   2023/07/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e   % Set the format to 5-digit floating point
clear
clc
close all

%% default value for easy debugging
defval('fam', '002');
defval('sps', 160);
defval('templensec', 60);
defval('ccmethod', 2);
defval('ccbp',[2 8]);
defval('plflag', 1);
defval('lo', 0.5);
defval('hi', 6.5);
defval('sps', 40);
defval('bef', 25);
defval('aft', 35);

%% Initialization
% set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]
wid=scrsz(3);       % width
hite=scrsz(4);       % height

% set path
workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
    'LZB'];
POLSTA=['SSIB '           % polaris station names
    'SILB '
    'KLNB '
    'MGCB '
    'TWKB '];

% station trio used, 1st station is the main station
stas=['PGC  '
    'SSIB '
    'SILB '
    'LZB  '
    'TWKB '
    'MGCB '
    'KLNB '];
% number of used stations
nsta=size(stas,1);         %  number of stations

fprintf('Family ID is %s \n', fam);
  
[timoffrot,tempoffs] = GetDays4Stack(fam);
tempoffs = round(tempoffs*sps/40);

freqflag='hf';  % flag to indicate whether to do hf or lf;
FLAG = 'PGC'; % detector
CATA = 'new'; %Use new LFE catalog

%%%%%%%%%%%% VERSION 3, FROM AUTOMATIC DAYS AND NEW CATALOG %%%%%%%%%%%%%%%%%
spsrot = 40;
PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',...
  num2str(aft),'PGC');
ROTS = load(strcat(datapath, '/split_chao/pgctrio/',...
  PREFIX,'_',num2str(spsrot),'sps',num2str(nsta),'stas_newcat'));


%% stack at bostock's detection timings on N and E components
%%% data parameters
% sps=40;     % samples per second
lenofday = 24*3600;
STAU=zeros(nsta,sps * lenofday);
STAN=zeros(nsta,sps * lenofday);
STAE=zeros(nsta,sps * lenofday);

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
% stackex = zeros(nsta, templen+ 2*extra);

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
nstack = 0;
nLFE = 0;
fprintf('Part 1: Raw Template \n');

dmatU = [];
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
    
    nlfe1day = size(bostocks, 1)
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
    direc=[datapath, '/arch', YEAR,'/',MO,'/'];     % directory name
    fprintf('%s \n',direc);
    datafnm = [direc, YEAR,'.',JDAY,'.00.00.00.0000.CN'];
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
            fname = strcat(datafnm,'.',PERMSTA(idx,1:3),'..BHZ.D.SAC');
            if isfile(fname)    % if have the data file
                % 1. this is for data without removing station response
                %   [opt,ort,timeperm]=readperm_1daynofilterv2(datafnm, PERMSTA, PERMROTS, idx, sps, fact);
                % 2. this is for data with no response
                dataU = readperms_norotsv2(datafnm, PERMSTA, idx, sps, fact, 'Z');
                [dataE, dataN, ~] = readperms_norotsv2(datafnm, PERMSTA, idx, sps, fact);
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
            fname = strcat(datafnm,'.',POLSTA(idx,1:4),'..HHZ.D.SAC');
            if isfile(fname)    % if have the data file
                %                 % 1. this is for data without removing station response
                %                 [opt,ort,nzeros]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
                % 2. this is for data with no response
                dataU = readpols_norotsv2(datafnm, POLSTA, idx, sps, fact, 'Z');
                [dataE, dataN, ~] = readpols_norotsv2(datafnm, POLSTA, idx, sps, fact);
            else
                fileflag = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for station %s in day %s / %s,this day will be omitted. \n',...
                        POLSTA(idx,1:4), YEAR, JDAY);
                break   % break the entire station loop
            end
            
        end
        STAU(ista, :) = dataU;
        STAN(ista, :) = dataN;
        STAE(ista, :) = dataE;
    end
    
    if fileflag == 0    % means there are missing files
        fprintf('Day %s / %s will be omitted because of missing files. \n', YEAR, JDAY);
        continue    % continue to the next day
    end
    
    len = size(STAU,2);
    for n=1: size(bostsamp, 1)
        % set the timing be at the center of the window
        if (bostsamp(n)+ after+ extra <= len) && ...
                (bostsamp(n)- before- extra >= 1)
            
            %%% it is better to use a longer window than templen to make
            %%% sure the shifting in step2 would preserve the frequency
            %%% content as much as possible
%             windata = STAU(:, bostsamp(n)- before- extra: ...
%                 bostsamp(n)+ after+ extra);
%             % remove mean and linear trend
%             tmp = detrend(windata');
%             stackexU = stackexU+ tmp';
%             
            nstack = nstack + 1;
            
            % datamat is 3D, for storing all qualified wins for stacking
            % size is (nsta, nstack, templen)
            dmatU(:, nstack, :) = STAU(:, bostsamp(n)- before- extra: ...
              bostsamp(n)+ after+ extra);
            dmatN(:, nstack, :) = STAN(:, bostsamp(n)- before- extra: ...
              bostsamp(n)+ after+ extra);
            dmatE(:, nstack, :) = STAE(:, bostsamp(n)- before- extra: ...
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
stackexU = zeros(nsta, templen+ 2*extra);
stackexN = zeros(nsta, templen+ 2*extra);
stackexE = zeros(nsta, templen+ 2*extra);
for i = 1: nstack
    % remove mean and linear trend
    windata(:,:) = dmatU(:, i, :);
    tmp = detrend(windata');
    stackexU = stackexU+ tmp';  
    windata(:,:) = dmatN(:, i, :);
    tmp = detrend(windata');
    stackexN = stackexN+ tmp';  
    windata(:,:) = dmatE(:, i, :);
    tmp = detrend(windata');
    stackexE = stackexE+ tmp';  
end

%% Averaging & normalization
%cut off extra points after aligning the templates at diff stations due to station loc difference
for ista = 1: nsta
  dstackU(ista,:) = stackexU(ista, 1+extra+ROTS(ista,4)*sps/spsrot: end-extra+ROTS(ista,4)*sps/spsrot);
  dstackN(ista,:) = stackexN(ista, 1+extra+ROTS(ista,4)*sps/spsrot: end-extra+ROTS(ista,4)*sps/spsrot);
  dstackE(ista,:) = stackexE(ista, 1+extra+ROTS(ista,4)*sps/spsrot: end-extra+ROTS(ista,4)*sps/spsrot);
end

%%% averaging
dstackU = dstackU/nstack;
dstackN = dstackN/nstack;
dstackE = dstackE/nstack;

%%% detrend
% tmp = detrend(dstackU');
dstackU = detrend(dstackU')';
dstackN = detrend(dstackN')';
dstackE = detrend(dstackE')';

mid = templen/2;

%% plotting if desired
if plflag
    %%% plot the 1-step stacked template
    %%% figure 1
    figure('Position',[0.1*wid 0.1*hite 0.3*wid 0.3*hite]);
    subplot(1,3,1); hold on; box on; grid on
    for ista=1: nsta
        plot(dstackU(ista, :) + 0.1* ista, 'linewidth', 2);
        title('direct stack');
    end
    subplot(1,3,2); hold on; box on; grid on
    for ista=1: nsta
        plot(dstackN(ista, :) + 0.1* ista, 'linewidth', 2);
        title('direct stack');
    end
    subplot(1,3,3); hold on; box on; grid on
    for ista=1: nsta
        plot(dstackE(ista, :) + 0.1* ista, 'linewidth', 2);
        title('direct stack');
    end
end

%% filter before all operations
npo = 2;
npa = 2;
hiwlet=18;
lowlet=1.8;
for ista =1: nsta
  dstackU(ista,:) = Bandpass(dstackU(ista,:), sps, lowlet, hiwlet, npo, npa, 'butter');
  dstackN(ista,:) = Bandpass(dstackN(ista,:), sps, lowlet, hiwlet, npo, npa, 'butter');
  dstackE(ista,:) = Bandpass(dstackE(ista,:), sps, lowlet, hiwlet, npo, npa, 'butter');
end

%% roughly cut a window around P and S-wave arrival
pwin = 4600+(-sps*6: -1);
swin = 4600+(0: sps*6-1);
% pseis = dstackU(:,pwin);
% sseis = dstackU(:,swin);
for ista = 1:nsta
  envU(ista,:) = envelope(dstackU(ista,:));
  envN(ista,:) = envelope(dstackN(ista,:));
  envE(ista,:) = envelope(dstackE(ista,:));
end
%Z COMP
[maxpenvu,maxipu] = max(envU(:,pwin),[],2);
maxipu = maxipu+pwin(1)-1;
[maxsenvu,maxisu] = max(envU(:,swin),[],2);
maxisu = maxisu+swin(1)-1;
spenvurat = maxsenvu./maxpenvu;
%N COMP
[maxpenvn,maxipn] = max(envN(:,pwin),[],2);
maxipn = maxipn+pwin(1)-1;
[maxsenvn,maxisn] = max(envN(:,swin),[],2);
maxisn = maxisn+swin(1)-1;
spenvnrat = maxsenvn./maxpenvn;
%E COMP
[maxpenve,maxipe] = max(envE(:,pwin),[],2);
maxipe = maxipe+pwin(1)-1;
[maxsenve,maxise] = max(envE(:,swin),[],2);
maxise = maxise+swin(1)-1;
spenverat = maxsenve./maxpenve;

%% plot after filtering 
figure
subplot(1,3,1)  % plot Z component, direct stack 
for ista = 1:nsta
  plot((1:templen)/sps, dstackU(ista,:)*10 + 1*ista); hold on
  plot((1:templen)/sps, envU(ista,:)*10 + 1*ista, 'k-');
  text(20.5,1*ista,stas(ista,:));
  text(25,1*ista+0.15,sprintf('%.2f',spenvurat(ista)));
end
ax=gca;
plot([pwin(1) pwin(1)]/sps,ax.YLim,'k--','linew',1);
plot([pwin(end) pwin(end)]/sps,ax.YLim,'k--','linew',1);
title(sprintf('direct stack, Z comp, %.2f-%.2f Hz',lowlet, hiwlet));
xlabel('Time (s)');
axis([20 40 0 nsta+1]);

subplot(1,3,2)  % plot N component, direct stack 
for ista = 1:nsta
  plot((1:templen)/sps, dstackN(ista,:)*10 + 1*ista); hold on
  plot((1:templen)/sps, envN(ista,:)*10 + 1*ista, 'k-');
  text(20.5,1*ista,stas(ista,:));
  text(25,1*ista+0.15,sprintf('%.2f',spenvnrat(ista)));
end
ax=gca;
plot([pwin(1) pwin(1)]/sps,ax.YLim,'k--','linew',1);
plot([pwin(end) pwin(end)]/sps,ax.YLim,'k--','linew',1);
title(sprintf('direct stack, N comp, %.2f-%.2f Hz',lowlet, hiwlet));
xlabel('Time (s)');
axis([20 40 0 nsta+1]);

subplot(1,3,3)  % plot Z component, direct stack 
for ista = 1:nsta
  plot((1:templen)/sps, dstackE(ista,:)*10 + 1*ista); hold on
  plot((1:templen)/sps, envE(ista,:)*10 + 1*ista, 'k-');
  text(20.5,1*ista,stas(ista,:));
  text(25,1*ista+0.15,sprintf('%.2f',spenverat(ista)));
end
ax=gca;
plot([pwin(1) pwin(1)]/sps,ax.YLim,'k--','linew',1);
plot([pwin(end) pwin(end)]/sps,ax.YLim,'k--','linew',1);
title(sprintf('direct stack, E comp, %.2f-%.2f Hz',lowlet, hiwlet));
xlabel('Time (s)');
axis([20 40 0 nsta+1]);

vertical_cursors;
keyboard

%% PART 2: Multi-CC the 1-step template and each window to get a finer template
%%%%%%%%%%%%%%%% WAY 2, stack all, only use cc to correct time %%%%%%%%%%%%%%%%%
%%% the previous 'extra' is tO avoid adding zeros to the array due to shifting,
nstack2 = nstack;
% according to template, the duration of dipole seems to be ~0.5 sec (10
% samples), could even as small as 5 samples
offsetmax = sps/4;
ccstackU = zeros(nsta, templen);
%%% SEE NOTES for proper freq. band choices
npo = 2;
npa = 2;
% pwin = 4600+[-sps*6 -1];
% swin = 4600+[0 sps*6-1];
pwin = 4600+(-sps*6: -1);
swin = 4600+(0: sps*6-1);
coefmat = zeros(nsta, nstack);
lagmat = zeros(nsta, nstack);
fprintf('Part 2: Precise Template \n');
% figure
for ista = 1: nsta    % loop for station
  %     ista = 1;
  dummy=[];
%   dummy(:,:) = dmatU(ista,:,:);
  dummy(:,:) = dmatU(ista,:,1+extra+ROTS(ista,4)*sps/spsrot: end-extra+ROTS(ista,4)*sps/spsrot);
  windata = detrend(dummy');

  dummy2=[];
%   dummy2 = stackex(ista, :)/ nstack;
%   dummy2 = stackex(ista, 1+extra+ROTS(ista,4)*sps/spsrot: end-extra+ROTS(ista,4)*sps/spsrot)/ nstack;
  dummy2 = dstackU(ista, :);
  template = detrend(dummy2');
  %             tic
  fwindata = [];
  for istack = 1: nstack
    fwindata(:,istack) = Bandpass(windata(:,istack),sps,ccbp(1),ccbp(2),npo,npa,'butter');  % also contains detrend
  end
  %             toc
  ftemplate = Bandpass(template,sps,ccbp(1),ccbp(2),npo,npa,'butter');
%   %%% to SEE how bad is our windowed data before/after filtering
%   figure
%   subplot(4,1,1)
%   plot(windata(:, istack), 'b');
%   subplot(4,1,2)
%   plot(fwindata(:, istack), 'r'); hold on; ax=gca;
%   plot([pwin(1) pwin(1)],ax.YLim,'k--');
%   plot([pwin(end) pwin(end)],ax.YLim,'k--');  
%   subplot(4,1,3)
%   plot(template, 'b');
%   subplot(4,1,4)
%   plot(ftemplate, 'r'); hold on; ax=gca;
%   plot([pwin(1) pwin(1)],ax.YLim,'k--');
%   plot([pwin(end) pwin(end)],ax.YLim,'k--');
% %           pause(2)
%   %%%
  
  %%% e.g. if lag = argmax(xcorr(A, B)) = 4, A needs to shift 4 to
  %%% the left to align with B, (i.e. A-B=4), B is the reference station
  win = detrend(fwindata(pwin,:));
  tem = detrend(ftemplate(pwin,:));

  for istack = 1: nstack
    [coef, lag] = xcorr(win(:,istack), tem, offsetmax, 'coeff');
    [maxcoef, idx] = max(coef);
    lagsamp = lag(idx);
    coefmat(ista, istack) = maxcoef;
    lagmat(ista, istack) = lagsamp;
    
%     newdata(:) = windata(1+lagsamp+extra: lagsamp+extra+templen,istack);
    newdata(:) = dmatU(ista,istack,1+extra+ROTS(ista,4)*sps/spsrot+lagsamp: ...
      end-extra+ROTS(ista,4)*sps/spsrot+lagsamp);
    winlag(ista, istack) = lagsamp;
    tmp = detrend(newdata');
    ccstackU(ista, :) = ccstackU(ista, :)+ tmp';
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
  
  
end     % loop for station
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WAY 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Averaging and normalization 
%%% averaging
for ista = 1: nsta
    ccstackU(ista, :) = ccstackU(ista, :)/ nstack2;
end

%%% detrend
ccstackU = detrend(ccstackU')';

%% plotting if desired
if plflag
  %%% plot the 2-step stacked template
  %%% figure 1
  figure('Position',[0.1*wid 0.1*hite 0.3*wid 0.3*hite]);
  for ista=1: nsta
    plot(ccstackU(ista, :) + 0.1* ista, 'linewidth', 2); hold on
    title('CC stack');
  end
  box on
  grid on
  
  %%% figure 2
  figure('Position',[0.1*wid 0.3*hite 0.3*wid 0.3*hite]);
  for ista = 1: nsta
    plot(ccstackU(ista, mid-5*sps: mid+5*sps) + 0.1* ista, 'linewidth', 2); hold on
    title('CC stack-short window');
  end
  box on
  grid on
  
  %%% figure 3
  figure
  supertit(gca,'sample shift of preserved traces passing threshold');
  for ista = 1: nsta
    subplot(nsta,1,ista)
    scatter(1: nstack2, winlag(ista, 1:nstack2));
  end
  box on
  grid on
end

%% filterring before all operations
for ista =1: nsta
    ccstackU(ista,:) = Bandpass(ccstackU(ista,:), sps, lowlet, hiwlet, npo, npa, 'butter');
end

%% 
figure  % plot original stack, fig 1, in case your 'tempoffs' is not pretty accurate
subplot(1,2,2)  % plot E component, cc stack  
for ista = 1:nsta
  plot((1:templen)/sps, ccstackU(ista,:)*10 + 1*ista); hold on
  tmpenvU(ista,:) = envelope(ccstackU(ista,:));
  plot((1:templen)/sps, tmpenvU(ista,:)*10 + 1*ista, 'k-');
  text(0.5,1*ista,stas(ista,:));
end
title(sprintf('cc stack w/in pwin, %.2f-%.2f Hz',lowlet, hiwlet));
xlabel('Time (s)');

vertical_cursors
keyboard


%     %%% fam 002
%     fam = '002';
%     splitoff = 12;  % the twkb trio are poor, pgc trio and klnb are fine
%     staoff = 8;
%     mainsta = 1;
% %     lo = 1.25;
% %     hi = 6.5;
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
