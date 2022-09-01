% function get_filtering_effect_PGCv3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the function to get the filtering effect under different frequency
% band for PGC trio. Use the original seperated Bostock's
% catalog (BOSTOCK/NEW/002-246_2003.062') instead of the catalog in one file
% '/BOSTOCK/total_mag_detect_0000_cull_NEW.txt'. 
% 
% Use the similar way to do
% the cross correlatation in Allan's detection code with all 3 stations,
% under the constraint that the enclosed sum of the offset need to be zero
% so that it necessarily would reach the individual max cc.
%
%   filter with different bands when reading individual data, then stack.
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/09/22
% Last modified date:   2019/09/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all

% set(0,'DefaultFigureVisible','on');
set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]
wid=scrsz(3);       % width
hite=scrsz(4);       % height

% set path
workpath = getenv('ALLAN');
temppath = strcat(workpath, '/templates/');
if ~exist(temppath, 'dir')
    mkdir(temppath);
end
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
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

% stas=['PGC  '
%     'SSIB '
%     'SILB '
%     'LZB  '
%     'TWKB '
%     'MGCB '
%     'KLNB '];

stas=['PGC  '
      'SSIB '
      'SILB '];

% number of used stations
nsta=size(stas,1);         %  number of stations

% convert angles to rads
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;     % convert 2-3 columns to rad
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;

if isequal(fam,'002')
    scaleseisms=[1.0 0.76 0.95];       % scaleseisms scales seismograms
elseif isequal(fam,'068')
    scaleseisms=[1.0 1.0 1.0];
end

%% HF 
%%% data parameters
sps=80;     % samples per second
lenofday = 24*3600;
STAopt=zeros(nsta,sps * lenofday);
STAort=zeros(nsta,sps * lenofday);

% upsps = 80;
% [num, denom] = rat(80/sps);
wlensechf = 4;
wlenhf = wlensechf*sps;
lohf = 1.25;
hihf = 6.5;
npo = 2;
npa = 2;
mshift = 29;
loopoffmax = 2.1;
xcmaxAVEnmin = 0.44;

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
    prename = [direc, YEAR,'.',JDAY,'.00.00.00.0000.CN'];
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
            fname = strcat(prename,'.',PERMSTA(idx,1:3),'..BHE.D.SAC');
            if isfile(fname)    % if have the data file
                % 1. this is for data without removing station response
                %   [opt,ort,timeperm]=readperm_1daynofilterv2(datafnm, PERMSTA, PERMROTS, idx, sps, fact);
                % 2. this is for data with no response
%                 [opt, ort, timeperm] = readperm_nofilterv2(datafnm, ...
%                     PERMSTA, PERMROTS, idx, sps, fact);
                [opt,ort,~,~]=readpermsv2(prename,PERMSTA,PERMROTS,idx,sps,lohf,hihf,npo,npa,fact,1,wlenhf,1,1);
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
            fname = strcat(prename,'.',POLSTA(idx,1:4),'..HHE.D.SAC');
            if isfile(fname)    % if have the data file
                %                 % 1. this is for data without removing station response
                %                 [opt,ort,nzeros]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
                % 2. this is for data with no response
%                 [opt, ort, timepola] = readpola_nofilterv2(datafnm, ...
%                     POLSTA, POLROTS, idx, sps, fact);
                [opt,ort,~]=readpolsv2(prename,POLSTA,POLROTS,idx,sps,lohf,hihf,npo,npa,fact,1,wlenhf,1,1);
            else
                fileflag = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for station %s in day %s / %s. \n',POLSTA(idx,1:4), YEAR, JDAY);
                break   % break the entire station loop
            end
            
        end
        
        STAopt(ista, :) = opt/scaleseisms(ista);
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

%%% figure 1
figure('Position',[0.1*wid 0.1*hite 0.3*wid 0.3*hite]);
for ista=1: nsta
    plot(stack(ista, :) + 0.1* ista, 'linewidth', 2); hold on
    title('direct stack');
end
box on


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
figure('Position',[0.1*wid 0.2*hite 0.3*wid 0.3*hite]);
for ista = 1: nsta
    plot(stackno(ista, :) + 2* ista, 'linewidth', 2); hold on
    title('normalized direct stack');
end
box on

%%% figure 3
figure('Position',[0.1*wid 0.3*hite 0.3*wid 0.3*hite]);
for ista = 1: nsta
    plot(stack(ista, mid-5*sps: mid+5*sps) + 0.1* ista, 'linewidth', 2); hold on
    title('direct stack-short window');
end
box on

stackuse = stack(1:3,:);
stackauto = stackuse.*stackuse;
lenx = templen-2*mshift;     
stack12x = zeros(lenx, 2*mshift+1);    
stack13x = zeros(lenx, 2*mshift+1);    % stas: 1->PGC, 2->SSIB, 3->SILB
stack32x = zeros(lenx, 2*mshift+1);

for n=-mshift:mshift
    % PGC corr SSIB, 1+mshift:tracelen-mshift == 1+mshift-n:tracelen-mshift-n == lenx
    stack12x(:,n+mshift+1)=stackuse(1,1+mshift:templen-mshift).* ...
        stackuse(2,1+mshift-n:templen-mshift-n);
    % PGC corr SILB
    stack13x(:,n+mshift+1)=stackuse(1,1+mshift:templen-mshift).* ...
        stackuse(3,1+mshift-n:templen-mshift-n);
    % SILB corr SSIB
    stack32x(:,n+mshift+1)=stackuse(3,1+mshift:templen-mshift).* ...
        stackuse(2,1+mshift-n:templen-mshift-n);
end

sumstack12=zeros(1,2*mshift+1);   % PGC-SSIB
sumstack13=zeros(1,2*mshift+1);   % PGC-SILB
sumstack32=zeros(1,2*mshift+1);   % SILB-SSIB
sumstack1sq=zeros(1,2*mshift+1);  % "sq" are the running cumulative sum, to be used later by differencing the window edpoints
sumstack2sq=zeros(1,2*mshift+1);
sumstack3sq=zeros(1,2*mshift+1);
sumstack3Bsq=zeros(1,2*mshift+1);

istart = mid-wlenhf/2;
iend = istart+wlenhf-1;
sumstack12(1,:) = sum(stack12x(istart-mshift: iend-mshift, :));
sumstack13(1,:) = sum(stack13x(istart-mshift: iend-mshift, :));
sumstack32(1,:) = sum(stack32x(istart-mshift: iend-mshift, :));
sumstack1sq(1,:) = sum(stackauto(1, istart: iend));
sumstack3Bsq(1,:) = sum(stackauto(3, istart: iend));
for m = -mshift:mshift
    sumstack2sq(1,m+mshift+1)=sum(stackauto(2, istart-m: iend-m)); %+m??? (yes).
    sumstack3sq(1,m+mshift+1)=sum(stackauto(3, istart-m: iend-m));
end

%An attempt to bypass glitches in data.  Min value of good data typically ~10^{-2}
glitches=1.e-7;
sumstack1sq=max(sumstack1sq,glitches);
sumstack2sq=max(sumstack2sq,glitches);
sumstack3sq=max(sumstack3sq,glitches);    % return maximum between A and B

denomstack12n=realsqrt(sumstack1sq.*sumstack2sq);    % Real square root, An error is produced if X is negative
denomstack13n=realsqrt(sumstack1sq.*sumstack3sq);
denomstack32n=realsqrt(sumstack3Bsq.*sumstack2sq);

sumstack12n=sumstack12./denomstack12n;   % suffix 'n' means normalized
sumstack13n=sumstack13./denomstack13n;
sumstack32n=sumstack32./denomstack32n;

[xcmaxstack12n,imaxstack12]=max(sumstack12n,[],2);   %Integer-offset max cross-correlation
[xcmaxstack13n,imaxstack13]=max(sumstack13n,[],2);   % along row, max cc val and index in each window
[xcmaxstack32n,imaxstack32]=max(sumstack32n,[],2);

%Parabolic fit:
[xmaxstack12n,ymaxstack12n,aSTA12]=parabol(1,mshift,sumstack12n,imaxstack12); %Interpolated max cross-correlation
[xmaxstack13n,ymaxstack13n,aSTA13]=parabol(1,mshift,sumstack13n,imaxstack13);
[xmaxstack32n,ymaxstack32n,aSTA32]=parabol(1,mshift,sumstack32n,imaxstack32);

%Center them
imaxstack12cent=imaxstack12-mshift-1;  % "cent" is "centered"; imaxSTA12 is original 1: 2*mshift+1, corresponds to -mshift: mshift
imaxstack13cent=imaxstack13-mshift-1;
imaxstack32cent=imaxstack32-mshift-1;
%%% NOTICE: the right order of a closed 3-sta pair is +13, -12, +32, where 13 means 1-->3
iloopoff=imaxstack13cent-imaxstack12cent+imaxstack32cent; %How well does the integer loop close?
%
xmaxstack12n=xmaxstack12n-mshift-1;
xmaxstack13n=xmaxstack13n-mshift-1;
xmaxstack32n=xmaxstack32n-mshift-1;
loopoff=xmaxstack13n-xmaxstack12n+xmaxstack32n; %How well does the interpolated loop close?
xcmaxAVEn=(xcmaxstack12n+xcmaxstack13n+xcmaxstack32n)/3;
medxcmaxAVEn=median(xcmaxAVEn);
xmaxstack12ntmphf=xmaxstack12n;    % tmp == temporary
xmaxstack13ntmphf=xmaxstack13n;
xmaxstack32ntmphf=xmaxstack32n;

if xcmaxAVEn<xcmaxAVEnmin || abs(loopoff)>loopoffmax || isequal(abs(imaxstack12cent),mshift) || isequal(abs(imaxstack13cent),mshift) || isequal(abs(imaxstack32cent),mshift)
    xmaxstack12ntmphf=mshift+1; xmaxstack13ntmphf=mshift+1; xmaxstack32ntmphf=mshift+1; %dummy them, if these criteria are met
    disp('WRONG! The basic criteria is not met');
else
    xcmaxconprev=-99999.;  %used to be 0; not good with glitches
    %%% NOTE here are using offset imaxSTA12/13/32 without centering, i.e. [1, 2shift+1]
    % Usually, the loop is happened between imaxSTA12n +- floor(loopoffmax+1)
    % width of which is 2*floor(loopoffmax+1)
    %%% floor(2.5)=2; floor(-2.6)=-3
    for iSTA12 =     max(1,imaxstack12-floor(loopoffmax+1)): min(imaxstack12+floor(loopoffmax+1),2*mshift+1)
        for iSTA13 = max(1,imaxstack13-floor(loopoffmax+1)): min(imaxstack13+floor(loopoffmax+1),2*mshift+1)
            ibangon = (mshift+1)-iSTA13+iSTA12;     %%% SEE NOTES #2019/03/17# page 66 to understand better
            %%% i.e., -mshift <= -iSTA13+iSTA12 <= mshift
            if ibangon >= 1 && ibangon <= 2*mshift+1
                xcmaxcon=sumstack12n(iSTA12)+sumstack13n(iSTA13)+sumstack32n(ibangon);
                if xcmaxcon > xcmaxconprev
                    xcmaxconprev=xcmaxcon;
                    iSTA12bang=iSTA12;
                    iSTA13bang=iSTA13;
                end
            end
        end
    end
    %%% will result in the max xcmaxcon and corresponding iSTA12,
    %%% iSTA13, and save them into xcmaxconprev, iSTA12bang and iSTA13bang
    iSTA32bang=(mshift+1)-iSTA13bang+iSTA12bang;
    if abs(iSTA12bang-imaxstack12) <= loopoffmax && ...  %not sure if these 3 lines are satisfied automatically ...
            abs(iSTA13bang-imaxstack13) <= loopoffmax && ...  % SHOULD not be, i think, could be floor(loopoffmax+1) > loopoffmax
            abs(iSTA32bang-imaxstack32) <= loopoffmax && ...
            sumstack12n(iSTA12bang)+sumstack13n(iSTA13bang)+sumstack32n(iSTA32bang) >= 3*xcmaxAVEnmin   % xcmaxAVEnmin, predetermined
        %%% ALSO, sumsSTA12n(n,iSTA12bang) == sumsSTA12nn(iSTA12bang)
        
        xmaxstack12ntmphf=iSTA12bang-(mshift+1); %without interpolation this is just centering.
        xmaxstack13ntmphf=iSTA13bang-(mshift+1);
        xmaxstack32ntmphf=iSTA32bang-(mshift+1);
        
        %%% xcmaxAVEnbang is added by Chao, to distinguish from
        %%% xcmaxAVEn, because it is max average CC coef
        xcmaxAVEnbang=(sumstack12n(iSTA12bang)+sumstack13n(iSTA13bang)+sumstack32n(iSTA32bang))/3;
    end
end

xmaxstack12ntmphf
xmaxstack13ntmphf
xmaxstack32ntmphf
if xmaxstack13ntmphf-xmaxstack12ntmphf+xmaxstack32ntmphf ~=0
    disp('WRONG! Loopoff is not enclosed');
else
    disp('Loopoff is 0');
end
filthf = [xmaxstack12ntmphf;  xmaxstack13ntmphf];

%%%% end for hf 

%% LF
%%% data parameters
sps=80;     % samples per second
lenofday = 24*3600;
STAopt=zeros(nsta,sps * lenofday);
STAort=zeros(nsta,sps * lenofday);

% upsps = 80;
% [num, denom] = rat(80/sps);
wlenseclf = 16;
wlenlf = wlenseclf*sps;
lolf = 0.5;
hilf = 1.25;
npo = 2;
npa = 2;
cyclskip = 20;
mshift=19+cyclskip;
loopoffmax = 4;
xcmaxAVEnmin = 0.5;

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
    prename = [direc, YEAR,'.',JDAY,'.00.00.00.0000.CN'];
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
            fname = strcat(prename,'.',PERMSTA(idx,1:3),'..BHE.D.SAC');
            if isfile(fname)    % if have the data file
                % 1. this is for data without removing station response
                %   [opt,ort,timeperm]=readperm_1daynofilterv2(datafnm, PERMSTA, PERMROTS, idx, sps, fact);
                % 2. this is for data with no response
%                 [opt, ort, timeperm] = readperm_nofilterv2(datafnm, ...
%                     PERMSTA, PERMROTS, idx, sps, fact);
                [opt,ort,~,~]=readpermsv2(prename,PERMSTA,PERMROTS,idx,sps,lolf,hilf,npo,npa,fact,1,wlenlf,1,1);
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
            fname = strcat(prename,'.',POLSTA(idx,1:4),'..HHE.D.SAC');
            if isfile(fname)    % if have the data file
                %                 % 1. this is for data without removing station response
                %                 [opt,ort,nzeros]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
                % 2. this is for data with no response
%                 [opt, ort, timepola] = readpola_nofilterv2(datafnm, ...
%                     POLSTA, POLROTS, idx, sps, fact);
                [opt,ort,~]=readpolsv2(prename,POLSTA,POLROTS,idx,sps,lolf,hilf,npo,npa,fact,1,wlenlf,1,1);
            else
                fileflag = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for station %s in day %s / %s. \n',POLSTA(idx,1:4), YEAR, JDAY);
                break   % break the entire station loop
            end
            
        end
        
        STAopt(ista, :) = opt/scaleseisms(ista);
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

%%% figure 1
figure('Position',[0.1*wid 0.1*hite 0.3*wid 0.3*hite]);
for ista=1: nsta
    plot(stack(ista, :) + 0.1* ista, 'linewidth', 2); hold on
    title('direct stack');
end
box on


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
figure('Position',[0.1*wid 0.2*hite 0.3*wid 0.3*hite]);
for ista = 1: nsta
    plot(stackno(ista, :) + 2* ista, 'linewidth', 2); hold on
    title('normalized direct stack');
end
box on

%%% figure 3
figure('Position',[0.1*wid 0.3*hite 0.3*wid 0.3*hite]);
for ista = 1: nsta
    plot(stack(ista, mid-5*sps: mid+5*sps) + 0.1* ista, 'linewidth', 2); hold on
    title('direct stack-short window');
end
box on

stackuse = stack(1:3,:);
stackauto = stackuse.*stackuse;
lenx = templen-2*mshift;     
stack12x = zeros(lenx, 2*mshift+1);    
stack13x = zeros(lenx, 2*mshift+1);    % stas: 1->PGC, 2->SSIB, 3->SILB
stack32x = zeros(lenx, 2*mshift+1);

for n=-mshift:mshift
    % PGC corr SSIB, 1+mshift:tracelen-mshift == 1+mshift-n:tracelen-mshift-n == lenx
    stack12x(:,n+mshift+1)=stackuse(1,1+mshift:templen-mshift).* ...
        stackuse(2,1+mshift-n:templen-mshift-n);
    % PGC corr SILB
    stack13x(:,n+mshift+1)=stackuse(1,1+mshift:templen-mshift).* ...
        stackuse(3,1+mshift-n:templen-mshift-n);
    % SILB corr SSIB
    stack32x(:,n+mshift+1)=stackuse(3,1+mshift:templen-mshift).* ...
        stackuse(2,1+mshift-n:templen-mshift-n);
end

sumstack12=zeros(1,2*mshift+1);   % PGC-SSIB
sumstack13=zeros(1,2*mshift+1);   % PGC-SILB
sumstack32=zeros(1,2*mshift+1);   % SILB-SSIB
sumstack1sq=zeros(1,2*mshift+1);  % "sq" are the running cumulative sum, to be used later by differencing the window edpoints
sumstack2sq=zeros(1,2*mshift+1);
sumstack3sq=zeros(1,2*mshift+1);
sumstack3Bsq=zeros(1,2*mshift+1);

istart = mid-wlenlf/2;
iend = istart+wlenlf-1;
sumstack12(1,:) = sum(stack12x(istart-mshift: iend-mshift, :));
sumstack13(1,:) = sum(stack13x(istart-mshift: iend-mshift, :));
sumstack32(1,:) = sum(stack32x(istart-mshift: iend-mshift, :));
sumstack1sq(1,:) = sum(stackauto(1, istart: iend));
sumstack3Bsq(1,:) = sum(stackauto(3, istart: iend));
for m = -mshift:mshift
    sumstack2sq(1,m+mshift+1)=sum(stackauto(2, istart-m: iend-m)); %+m??? (yes).
    sumstack3sq(1,m+mshift+1)=sum(stackauto(3, istart-m: iend-m));
end

%An attempt to bypass glitches in data.  Min value of good data typically ~10^{-2}
glitches=1.e-7;
sumstack1sq=max(sumstack1sq,glitches);
sumstack2sq=max(sumstack2sq,glitches);
sumstack3sq=max(sumstack3sq,glitches);    % return maximum between A and B

denomstack12n=realsqrt(sumstack1sq.*sumstack2sq);    % Real square root, An error is produced if X is negative
denomstack13n=realsqrt(sumstack1sq.*sumstack3sq);
denomstack32n=realsqrt(sumstack3Bsq.*sumstack2sq);

sumstack12n=sumstack12./denomstack12n;   % suffix 'n' means normalized
sumstack13n=sumstack13./denomstack13n;
sumstack32n=sumstack32./denomstack32n;

[xcmaxstack12n,imaxstack12]=max(sumstack12n,[],2);   %Integer-offset max cross-correlation
[xcmaxstack13n,imaxstack13]=max(sumstack13n,[],2);   % along row, max cc val and index in each window
[xcmaxstack32n,imaxstack32]=max(sumstack32n,[],2);

%Parabolic fit:
[xmaxstack12n,ymaxstack12n,aSTA12]=parabol(1,mshift,sumstack12n,imaxstack12); %Interpolated max cross-correlation
[xmaxstack13n,ymaxstack13n,aSTA13]=parabol(1,mshift,sumstack13n,imaxstack13);
[xmaxstack32n,ymaxstack32n,aSTA32]=parabol(1,mshift,sumstack32n,imaxstack32);

%Center them
imaxstack12cent=imaxstack12-mshift-1;  % "cent" is "centered"; imaxSTA12 is original 1: 2*mshift+1, corresponds to -mshift: mshift
imaxstack13cent=imaxstack13-mshift-1;
imaxstack32cent=imaxstack32-mshift-1;
%%% NOTICE: the right order of a closed 3-sta pair is +13, -12, +32, where 13 means 1-->3
iloopofflf=imaxstack13cent-imaxstack12cent+imaxstack32cent; %How well does the integer loop close?
%
xmaxstack12n=xmaxstack12n-mshift-1;
xmaxstack13n=xmaxstack13n-mshift-1;
xmaxstack32n=xmaxstack32n-mshift-1;
loopofflf=xmaxstack13n-xmaxstack12n+xmaxstack32n; %How well does the interpolated loop close?
xcmaxAVEn=(xcmaxstack12n+xcmaxstack13n+xcmaxstack32n)/3;
medxcmaxAVEn=median(xcmaxAVEn);
xmaxstack12ntmplf=xmaxstack12n;    % tmp == temporary
xmaxstack13ntmplf=xmaxstack13n;
xmaxstack32ntmplf=xmaxstack32n;

if xcmaxAVEn<xcmaxAVEnmin || abs(loopofflf)>loopoffmax || isequal(abs(imaxstack12cent),mshift) || isequal(abs(imaxstack13cent),mshift) || isequal(abs(imaxstack32cent),mshift)
    xmaxstack12ntmplf=mshift+1; xmaxstack13ntmplf=mshift+1; xmaxstack32ntmplf=mshift+1; %dummy them, if these criteria are met
    disp('WRONG! The basic criteria is not met');
else
    xcmaxconprev=-99999.;  %used to be 0; not good with glitches
    %%% NOTE here are using offset imaxSTA12/13/32 without centering, i.e. [1, 2shift+1]
    % Usually, the loop is happened between imaxSTA12n +- floor(loopoffmax+1)
    % width of which is 2*floor(loopoffmax+1)
    %%% floor(2.5)=2; floor(-2.6)=-3
    for iSTA12 =     max(1,imaxstack12-floor(loopoffmax+1)): min(imaxstack12+floor(loopoffmax+1),2*mshift+1)
        for iSTA13 = max(1,imaxstack13-floor(loopoffmax+1)): min(imaxstack13+floor(loopoffmax+1),2*mshift+1)
            ibangon = (mshift+1)-iSTA13+iSTA12;     %%% SEE NOTES #2019/03/17# page 66 to understand better
            %%% i.e., -mshift <= -iSTA13+iSTA12 <= mshift
            if ibangon >= 1 && ibangon <= 2*mshift+1
                xcmaxcon=sumstack12n(iSTA12)+sumstack13n(iSTA13)+sumstack32n(ibangon);
                if xcmaxcon > xcmaxconprev
                    xcmaxconprev=xcmaxcon;
                    iSTA12bang=iSTA12;
                    iSTA13bang=iSTA13;
                end
            end
        end
    end
    %%% will result in the max xcmaxcon and corresponding iSTA12,
    %%% iSTA13, and save them into xcmaxconprev, iSTA12bang and iSTA13bang
    iSTA32bang=(mshift+1)-iSTA13bang+iSTA12bang;
    if abs(iSTA12bang-imaxstack12) <= loopoffmax && ...  %not sure if these 3 lines are satisfied automatically ...
            abs(iSTA13bang-imaxstack13) <= loopoffmax && ...  % SHOULD not be, i think, could be floor(loopoffmax+1) > loopoffmax
            abs(iSTA32bang-imaxstack32) <= loopoffmax && ...
            sumstack12n(iSTA12bang)+sumstack13n(iSTA13bang)+sumstack32n(iSTA32bang) >= 3*xcmaxAVEnmin   % xcmaxAVEnmin, predetermined
        %%% ALSO, sumsSTA12n(n,iSTA12bang) == sumsSTA12nn(iSTA12bang)
        
        xmaxstack12ntmplf=iSTA12bang-(mshift+1); %without interpolation this is just centering.
        xmaxstack13ntmplf=iSTA13bang-(mshift+1);
        xmaxstack32ntmplf=iSTA32bang-(mshift+1);
        
        %%% xcmaxAVEnbang is added by Chao, to distinguish from
        %%% xcmaxAVEn, because it is max average CC coef
        xcmaxAVEnbang=(sumstack12n(iSTA12bang)+sumstack13n(iSTA13bang)+sumstack32n(iSTA32bang))/3;
    end
end

xmaxstack12ntmplf
xmaxstack13ntmplf
xmaxstack32ntmplf
if xmaxstack13ntmplf-xmaxstack12ntmplf+xmaxstack32ntmplf ~=0
    disp('WRONG! Loopoff is not enclosed');
else
    disp('Loopoff in lf is 0');
end
filtlf = [xmaxstack12ntmplf;  xmaxstack13ntmplf];

%%%% end for lf

%%% save the filtering effect result
savemat = [filthf filtlf filthf-filtlf];
PREFIX = strcat(fam,'.lf',num2str(lolf),'-',num2str(hilf),'.hf',num2str(lohf),'-',num2str(hihf),...
    '.hfwin',num2str(wlensechf),'.lfwin',num2str(wlenseclf),'.',num2str(sps),'sps','.noresp');
fid = fopen(strcat(rstpath, '/MAPS/template_filtering_effect_PGC_enclosed','_',PREFIX),'w+');
fprintf(fid,'%d %d %d \n',savemat');
fclose(fid);
