%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the script to merge the detections from all families if you want 
% to remove duplicates base on time and CC coef first (BEFORE hypoinverse),
% then set a distance cutoff after inversion for lat/lon or km;
% The input is already been removed off the double counting. Since for each
% fam, it wiped put the same data, so the detections should still have a lot of
% duplicates, preserve the ones within the one duration that have the highest 
% CC value. In the mean time, detections in each family have a different
% reference (0,0) point, so the merging can also convert all to the same
% frame (similar to changing centroid).
%
% NOTE :
%   1*. This should be ran right after hf&lfdetection and before 
%       plotaddcheck_LZB_allfam 
%       IF want to set a distance cutoff first, then merge fams to 
%       remove duplicates base on time and CC coef, SKIP this script and go ahead
%       to run plotaddcheck_LZB_allfam, while the fcflag should be 1
%        
%   2. run with fcflag=0 (DEFAULT) would save the no-duplicate results for 
%       plotaddcheck_LZB_allfam and figures; run with fcflag=1 would correct 
%       the filtering effect for no-duplicate results to generate and save new
%        figures, but not save files
%   3. run with convflag=0 (DEFAULT) would keep the results in its own fam frame;
%       run with convflag=1 would convert results of all fams to the same frame
%       of the reference fam, which is now 043, can be changed easily.
%   4. IN DEFAULT, in order to write rst files for next-step plotaddcheck scripts, set
%       fcflag=0 and convflag=0; other options ONLY save plots to eyeball analysis
%   
%        
% Conversion algorithm:
%       (off12')_fam - (off12')_ref = (off12)_fam - (off12)_ref
%   where ' means the detection in that family, no ' means the original
%   reference point in each family, i.e. the location of that family 
%   (considered as 0,0 in each fam). Since we already have the travel time
%   difference between sta 12 & 13 in each family (4th col in rots para), 
%   what we need to do here is to convert them to a uniform frame.
%   Thus, relative to a selected fam, all detections in other fam can be
%   written in the new frame as:
%       (off12')_ref = (off12')_fam - [(off12)_fam - (off12)_ref]
%
% Features:
%   1. read detection offsets and trace files
%   2. option flag to correct filtering effect (1/0)
%   3. option flag to convert all fams to the same frame (1/0)
%   4. remove double counting, duplicates, NECESSARY
%   5. re-write those files (only when fcflag=0 and convflag=0)
%   6. re-plot the offsets variation w/ time after removing duplicates (always save plots)
%   
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/09/02
% Last modified date:   2019/11/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
scrsz=get(0,'ScreenSize');

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band 
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
%%% choose the data path here 
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/LZBtrio');

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
    
stas=['TWKB '
      'LZB  '
      'MGCB '];     % determine the trio and order 
  
fcflag = 1;     % flag to do the filtering effect correction

convflag = 1;   % flag to convert all offsets to same ref fam

FLAG = 'LZB';

%% Store detection of all days for all fam, & do the filtering effect correction if want
%%% for hf
hfallold = [];
spshf=40;     % samples per second
wlensechf=4;
woffsechf=1;        % window offset in sec, which is the step of a moving window
wlenhf=wlensechf*spshf;      % length in smaples
hihf=6.5;
lohf=1.25;
npohf=2;     % poles, passes of filters
npahf=2;
cskiphf = 0;
msfthf=14+cskiphf; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
loffhf=1.5; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
xchf=0.4;

%%% for lf
lfallold = [];
spslf=20;     % samples per second
wlenseclf=16;       % CHECK in detection script for what were used in first-year report
woffseclf=8;        % window offset in sec, which is the step of a moving window
wlenlf=wlenseclf*spslf;      % length in smaples
hilf=1.25;
lolf=0.5;
npolf=2;     % poles, passes of filters
npalf=2;
cskiplf = 20;
msftlf=14+cskiplf; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
lofflf=4; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
xclf=0.35;

timoffrot = [2003 060;
            2003 061;
            2003 062; % dates of the year, two columns: year Julian_day
            2003 063;
            2003 064;
            %             2003 065;     % LZB only has data from 2003 060-064
            %             2003 066;
            %             2003 067;
            2004 194;
            2004 195;
            2004 196;
            2004 197;
            2004 198;
            2004 199;
            2004 200;
            2004 201;
            2004 202;
            2004 203;
            2005 254;
            2005 255;
            2005 256;
            2005 257;
            2005 258;
            2005 259;
            2005 260;
            2005 261];
nday = size(timoffrot, 1);

%% read and proceed HF 
disp('Start reading hf files ...');

for ifam = 1: length(nfampool)
    fam = nfampool(ifam, :);    
%     [timoffrot,~] = GetDays4Stack(fam);
    for id=1:nday    % num of rows, also num of days       
        %Which days of data to read?
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
        
        %%% load hf detection result
        IDENTIF=[YEAR,'.',JDAY,'.',fam,'.lo',num2str(loffhf),'.cc',num2str(xchf),'.',...
                 int2str(npohf),int2str(npahf),'.ms', int2str(msfthf)]
        iuphf = 4;
        if iuphf == 1
            fhf1 = strcat(rstpath,'/MAPS/mapallrst',IDENTIF,'_',num2str(lohf),'-',num2str(hihf),'_',...
                            num2str(wlenhf/spshf),'s',num2str(spshf),'sps', '4nsta');
            fhf2 = strcat(rstpath,'/MAPS/mapallrst',IDENTIF,'_',num2str(lohf),'-',num2str(hihf),'_',...
                            num2str(wlenhf/spshf),'s',num2str(spshf),'sps', '3nsta');
        else
            fhf1 = strcat(rstpath,'/MAPS/mapallrst_up',IDENTIF,'_',num2str(lohf),'-',num2str(hihf),'_',...
                            num2str(wlenhf/spshf),'s',num2str(spshf),'sps', '4nsta');
            fhf2 = strcat(rstpath,'/MAPS/mapallrst_up',IDENTIF,'_',num2str(lohf),'-',num2str(hihf),'_',...
                            num2str(wlenhf/spshf),'s',num2str(spshf),'sps', '3nsta');
        end
        %%%% IT is important to remember the strcture of your data array
        %%% structure of allrst:
        %%%   timswin(n) xmaxSTA13ntmp(n) xmaxSTA12ntmp(n) xcmaxAVEnbang(nin) loopoff(n)
        %%%   cumsumtrdiff timswin(n)-winlensec/2+idiff/sps cumsumtrdiff/cumsumtr(winlen)
        %%%   match(nin,1) match(nin,2) sumsSTA12n(n,iSTA12bang) sumsSTA13n(n,iSTA13bang)
        %%%   sumsSTA32n(n,iSTA32bang) in(1:nstanew) loff(1:nstanew) ioff(1:nstanew) ccmaxave(1:nstanew)
        hftemp = [];
        if isfile(fhf1)
            hftemp = load(fhf1);
        elseif isfile(fhf2)
            hftemp = load(fhf2);
        else
            fprintf('Day %s %s of fam %s is skipped because of no hf detections.\n',YEAR, JDAY, fam);
            continue
        end
                
        %%% whether or not to do filtering effect correction
        if fcflag
            disp('Start correcting the filtering effect ...');
            % load filtering effect result file, NOTE that the sample rate used in obtaining the
            % filtering effect is 80, might need to change accordingly
             PREFIX = strcat(fam,'.lf',num2str(lolf),'-',num2str(hilf),'.hf',num2str(lohf),'-', ...
                             num2str(hihf),'.hw',num2str(wlensechf),'.lw',num2str(wlenseclf),...
                             '.lo',num2str(lofflf),'.ccm',num2str(xclf),'.','80sps');
            fname = strcat(rstpath, '/MAPS/tempfeff_LZB_enc','_',PREFIX);
            filtcor = load(fname);
        
            spsrathf = spshf/80;    
            filthf = filtcor(:,1)*spsrathf;   % hf template shift due to filterig effect, sign convention is same
            hftemp(:,2) = hftemp(:,2)-filthf(2);    % 2nd col is offset 13. i.e. the 2nd filthf
            hftemp(:,3) = hftemp(:,3)-filthf(1);    % 3rd col is offset 12
        end
        
        %%% whether to convert offsets from different fams to a uniform frame
        %%% ALGORITHM is confirmed to be true!
        if convflag
            disp('Start convert offsets from all fams to a uniform frame ...');
            CATA = 'new';
            [PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,0,0);
            reftfam = POLROTS(5,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
            PERMROTS(:,4) = PERMROTS(:,4)-reftfam;    % to make sure that 1st station is 0
            POLROTS(:,4) = POLROTS(:,4)-reftfam;
            
            reffam = '043';     %043
            [refperm, refpol] = GetRotsCommon(FLAG,reffam,CATA,datapath,0,0);
            reftref = refpol(5,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
            refperm(:,4) = refperm(:,4)-reftref;    % to make sure that 1st station is 0
            refpol(:,4) = refpol(:,4)-reftref;
            
            % conversion time constants, NOTE that the slow/fast time offset and offset between stas
            % is all calculated based on 40 sps
            conv12 = PERMROTS(2,4) - refperm(2,4);   % station 2, LZB
            conv13 = POLROTS(4,4) - refpol(4,4);    % station 3, MGCB
            
            spsrathf = spshf/40;    
            hftemp(:,2) = hftemp(:,2)-conv13*spsrathf;  % get the new off based on the ref fam frame
            hftemp(:,3) = hftemp(:,3)-conv12*spsrathf;
            
        end
        
        %%% store result
        famcol = ones(size(hftemp,1),1)*str2double(fam);
        datecol = ones(size(hftemp,1),1)*str2double(strcat(YEAR,JDAY));
        hftemp = [famcol datecol hftemp];
        hfallold = [hfallold; hftemp];
                
    end
end

%% read and proceed LF 
disp('Start reading lf files ...');
for ifam = 1: length(nfampool)
    fam = nfampool(ifam, :);    
%     [timoffrot,~] = GetDays4Stack(fam);
    for id=1:nday    % num of rows, also num of days       
        %Which days of data to read?
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
        
        %%% load hf detection result
        IDENTIF=[YEAR,'.',JDAY,'.',fam,'.lo',num2str(lofflf),'.cc',num2str(xclf),'.',...
                 int2str(npolf),int2str(npalf),'.ms', int2str(msftlf)]
        iuplf = 1;
        if iuplf == 1
            flf1 = strcat(rstpath,'/MAPS/mapallrst',IDENTIF,'_',num2str(lolf),'-',num2str(hilf),'_',...
                          num2str(wlenlf/spslf),'s',num2str(spslf),'sps', '4nsta');
            flf2 = strcat(rstpath,'/MAPS/mapallrst',IDENTIF,'_',num2str(lolf),'-',num2str(hilf),'_',...
                          num2str(wlenlf/spslf),'s',num2str(spslf),'sps', '3nsta');
        else
            flf1 = strcat(rstpath,'/MAPS/mapallrst_up',IDENTIF,'_',num2str(lolf),'-',num2str(hilf),'_',...
                          num2str(wlenlf/spslf),'s',num2str(spslf),'sps', '4nsta');
            flf2 = strcat(rstpath,'/MAPS/mapallrst_up',IDENTIF,'_',num2str(lolf),'-',num2str(hilf),'_',...
                          num2str(wlenlf/spslf),'s',num2str(spslf),'sps', '3nsta');
        end
        lftemp = [];
        if isfile(flf1)
            lftemp = load(flf1);
        elseif isfile(flf2)
            lftemp = load(flf2);
        else
            fprintf('Day %s %s of fam %s is skipped because of no lf detections.\n',YEAR, JDAY, fam);
            continue
        end
        
        %%% whether or not to do filtering effect correction
        if fcflag
            disp('Start correcting the filtering effect ...');
            % load filtering effect result file, NOTE that the sample rate used in obtaining the
            % filtering effect is 80, might need to change accordingly
             PREFIX = strcat(fam,'.lf',num2str(lolf),'-',num2str(hilf),'.hf',num2str(lohf),'-', ...
                             num2str(hihf),'.hw',num2str(wlensechf),'.lw',num2str(wlenseclf),...
                             '.lo',num2str(lofflf),'.ccm',num2str(xclf),'.','80sps');
            fname = strcat(rstpath, '/MAPS/tempfeff_LZB_enc','_',PREFIX);
            filtcor = load(fname);
        
            spsratlf = spslf/80;
            filtlf = filtcor(:,2)*spsratlf;   % hf template shift due to filterig effect, sign convention is same
            lftemp(:,2) = lftemp(:,2)-filtlf(2);    % 2nd col is offset 13. i.e. the 2nd filthf
            lftemp(:,3) = lftemp(:,3)-filtlf(1);    % 3rd col is offset 12
            
        end
                    
        %%% whether to convert offsets from different fams to a uniform frame
        %%% ALGORITHM is confirmed to be true!
        if convflag
            disp('Start convert offsets from all fams to a uniform frame ...');
            CATA = 'new';
            [PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,0,0);
            reftfam = POLROTS(5,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
            PERMROTS(:,4) = PERMROTS(:,4)-reftfam;    % to make sure that 1st station is 0
            POLROTS(:,4) = POLROTS(:,4)-reftfam;
            
            reffam = '043';     %043
            [refperm, refpol] = GetRotsCommon(FLAG,reffam,CATA,datapath,0,0);
            reftref = refpol(5,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
            refperm(:,4) = refperm(:,4)-reftref;    % to make sure that 1st station is 0
            refpol(:,4) = refpol(:,4)-reftref;
            
            % conversion time constants, NOTE that the slow/fast time offset and offset between stas
            % is all calculated based on 40 sps
            conv12 = PERMROTS(2,4) - refperm(2,4);   % station 2, LZB
            conv13 = POLROTS(4,4) - refpol(4,4);    % station 3, MGCB    
                        % lf
            spsratlf = spslf/40;
            lftemp(:,2) = lftemp(:,2)-conv13*spsratlf;  % note the diff sps in lf, 
            lftemp(:,3) = lftemp(:,3)-conv12*spsratlf;
        end
        
        %%% store result
        famcol = ones(size(lftemp,1),1)*str2double(fam);
        datecol = ones(size(lftemp,1),1)*str2double(strcat(YEAR,JDAY));
        lftemp = [famcol datecol lftemp];
        lfallold = [lfallold; lftemp];

    end
end


%% sort according to time, remove double counting
disp('Start removing double counts ...');

%%% for hf
%%% sort according to date and main arrival time, could add offset
hfallsort = sortrows(hfallold, [2,9,4,5]);
dateall = unique(hfallsort(:,2));
hfallnew = [];
for id = 1: length(dateall)
    hfdaysort = hfallsort(hfallsort(:,2) == dateall(id), :);
    dtminhf = 0.5;      % min time during which only 1 detection is retained
    colnum = [9 6 28 29 30 31];
    indexthf = length(nfampool);    % the extreme case is that every fam contains a duplicate
    [hfdaynew,indsave,inddis] = RemoveDoubleCounting(hfdaysort,dtminhf,indexthf,colnum);    
    hfallnew = [hfallnew; hfdaynew];    % hfdaynew is the file of all fam for one day  
    
end
        
%%% for lf
%%% sort according to date and main arrival time, could add offset
lfallsort = sortrows(lfallold, [2,9,4,5]);
dateall = unique(lfallsort(:,2));
lfallnew = [];
for id = 1: length(dateall)
    lfdaysort = lfallsort(lfallsort(:,2) == dateall(id), :);
    dtminlf = 1.1;      % min time during which only 1 detection is retained
    colnum = [9 6 28 29 30 31];
    indextlf = length(nfampool);
    [lfdaynew,~,~] = RemoveDoubleCounting(lfdaysort,dtminlf,indextlf,colnum);
    lfallnew = [lfallnew; lfdaynew];
    
end



%% rewrite all results and plot if needed

%%% if convflag == 0, means arrays need to be sorted that results in each day of each fam would be
%%% saved separately (when fcflag == 0), so do the plots
if convflag == 0
    if fcflag == 0
        disp('Start re-writing files ...');
    end
    % sort according to fam number,date,time
    hfallnew = sortrows(hfallnew, [1,2,3]);
    lfallnew = sortrows(lfallnew, [1,2,3]);
    for ifam = 1: length(nfampool)
        fam = nfampool(ifam, :)
        nday = size(timoffrot, 1);
        for id = 1: nday
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
            date = str2double(strcat(YEAR,JDAY));
            
            %%% for hf, re-assign the results to each day in each fam
            rstmathf = hfallnew(hfallnew(:,1)==str2double(fam)&hfallnew(:,2)==date, 3:end);
            nin_new = size(rstmathf,1);
            if isempty(rstmathf)
                fprintf('Date %s in fam %s is empty! \n',date, fam);
            else
                rstoldhf = hfallold(hfallold(:,1)==str2double(fam)&hfallold(:,2)==date, 3:end);
                %%% obtain the indexes of remaining detections, in order to
                %%% preserve the corresponding seismic traces only
                [~, iold, ~] = intersect(rstoldhf, rstmathf, 'row', 'stable');
                
                %%% Write results into files
                IDENTIF=[YEAR,'.',JDAY,'.',fam,'.lo',num2str(loffhf),'.cc',num2str(xchf),'.',...
                    int2str(npohf),int2str(npahf),'.ms', int2str(msfthf)]
                if isfile(fhf1) && fcflag==0  % save files when only removing double counting is used
                    %%% save detections
                    if iuphf ==1
                        fid = fopen(strcat(rstpath,'/MAPS/mapallrstNoDou',IDENTIF,'_',num2str(lohf),...
                              '-',num2str(hihf),'_',num2str(wlenhf/spshf),'s',num2str(spshf),'sps',...
                              '4nsta'),'w+');
                    else
                        fid = fopen(strcat(rstpath,'/MAPS/mapallrstUpNoDou',IDENTIF,'_',num2str(lohf),...
                              '-',num2str(hihf),'_',num2str(wlenhf/spshf),'s',num2str(spshf),'sps',...
                              '4nsta'),'w+');
                    end
                    fprintf(fid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %10.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f \n',rstmathf');
                    fclose(fid);
                    %%% read old seismic traces and save new ones
                    traceoldhf = load(strcat(rstpath, '/MAPS/traceall_',IDENTIF,'_',num2str(lohf),...
                        '-',num2str(hihf),'_',num2str(wlenhf/spshf),'s',num2str(spshf),...
                        'sps', '4nsta'));
                    tracenewhf = zeros(nin_new*wlenhf, size(traceoldhf,2));
                    for i = 1: nin_new
                        tracenewhf((i-1)*wlenhf+1: i*wlenhf,:)= traceoldhf((iold(i)-1)*wlenhf+1 : ...
                            iold(i)*wlenhf, :);
                    end
                    if iuphf ==1
                        fid = fopen(strcat(rstpath, '/MAPS/traceallNoDou_',IDENTIF,'_',...
                                    num2str(lohf),'-',num2str(hihf),'_',num2str(wlenhf/spshf),'s',...
                                    num2str(spshf),'sps', '4nsta'),'w+');
                    else
                        fid = fopen(strcat(rstpath, '/MAPS/traceallUpNoDou_',IDENTIF,'_',...
                                    num2str(lohf),'-',num2str(hihf),'_',num2str(wlenhf/spshf),'s',...
                                    num2str(spshf),'sps', '4nsta'),'w+');
                    end
                    fprintf(fid,'%.4f %.6f %.6f %.6f %.6f %.6f %.6f %.6f \n',tracenewhf');
                    fclose(fid);
                    
                elseif isfile(fhf2) && fcflag==0  % save files when only removing double counting is used
                    %%% save detections
                    if iuphf ==1
                        fid = fopen(strcat(rstpath,'/MAPS/mapallrstNoDou',IDENTIF,'_',num2str(lohf),...
                              '-',num2str(hihf),'_',num2str(wlenhf/spshf),'s',num2str(spshf),'sps',...
                              '3nsta'),'w+');
                    else
                        fid = fopen(strcat(rstpath,'/MAPS/mapallrstUpNoDou',IDENTIF,'_',num2str(lohf),...
                              '-',num2str(hihf),'_',num2str(wlenhf/spshf),'s',num2str(spshf),'sps',...
                              '3nsta'),'w+');
                    end
                    fprintf(fid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %10.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %d %d %d %.4f %.4f %.4f %d %d %d %.3f %.3f %.3f \n',rstmathf');
                    fclose(fid);
                    %%% read old seismic traces and save new ones
                    traceoldhf = load(strcat(rstpath, '/MAPS/traceall_',IDENTIF,'_',num2str(lohf),...
                        '-',num2str(hihf),'_',num2str(wlenhf/spshf),'s',num2str(spshf),...
                        'sps','3nsta'));
                    tracenewhf = zeros(nin_new*wlenhf, size(traceoldhf,2));
                    for i = 1: nin_new
                        tracenewhf((i-1)*wlenhf+1: i*wlenhf,:)= traceoldhf((iold(i)-1)*wlenhf+1 : ...
                            iold(i)*wlenhf, :);
                    end
                    if iuphf ==1
                        fid = fopen(strcat(rstpath, '/MAPS/traceallNoDou_',IDENTIF,'_',...
                                    num2str(lohf),'-',num2str(hihf),'_',num2str(wlenhf/spshf),'s',...
                                    num2str(spshf),'sps', '3nsta'),'w+');
                    else
                        fid = fopen(strcat(rstpath, '/MAPS/traceallUpNoDou_',IDENTIF,'_',...
                                    num2str(lohf),'-',num2str(hihf),'_',num2str(wlenhf/spshf),'s',...
                                    num2str(spshf),'sps', '3nsta'),'w+');
                    end
                    fprintf(fid,'%.4f %.6f %.6f %.6f %.6f %.6f %.6f \n',tracenewhf');
                    fclose(fid);
                end
                
                %%% plot the hf offset variation against time
                %%% if convflag = 0, plots have to be separated by fam and day!!
                disp('Start ploting and saving figures...')
                f.fig = figure(id);
                set(f.fig,'Position',[scrsz(3)/10 scrsz(4)/10 2*scrsz(3)/5 1*scrsz(4)/2]);
                ymax = 40;
                msize = 1;
                f.ax(1) = subplot(4,1,1,'align');
                hold(f.ax(1),'on');
                plot(f.ax(1),rstmathf(:,1),rstmathf(:,3)*40/spshf,'bs','MarkerSize',msize);
                plot(f.ax(1),rstmathf(:,1),rstmathf(:,2)*40/spshf,'ro','MarkerSize',msize);
                axis(f.ax(1),[0 86400/4 -ymax ymax]);
                yticks(f.ax(1),-ymax:20:ymax);
                %             ylabel(f.ax(1),'Samples converted to 40-sps base')
                title(f.ax(1),strcat(fam,'-',YEAR,JDAY,'-',num2str(lohf),'-',num2str(hihf)));
                grid(f.ax(1),'on');
                box(f.ax(1),'on');
                hold(f.ax(1),'off');
                
                f.ax(2) = subplot(4,1,2,'align');
                hold(f.ax(2),'on');
                plot(f.ax(2),rstmathf(:,1),rstmathf(:,3)*40/spshf,'bs','MarkerSize',msize);
                plot(f.ax(2),rstmathf(:,1),rstmathf(:,2)*40/spshf,'ro','MarkerSize',msize);
                axis(f.ax(2),[86400/4 86400/2 -ymax ymax]);
                yticks(f.ax(2),-ymax:20:ymax);
                ylabel(f.ax(2),'Samples converted to 40-sps base')
                grid(f.ax(2),'on');
                box(f.ax(2),'on');
                hold(f.ax(2),'off');
                
                f.ax(3) = subplot(4,1,3,'align');
                hold(f.ax(3),'on');
                plot(f.ax(3),rstmathf(:,1),rstmathf(:,3)*40/spshf,'bs','MarkerSize',msize);
                plot(f.ax(3),rstmathf(:,1),rstmathf(:,2)*40/spshf,'ro','MarkerSize',msize);
                axis(f.ax(3),[86400/2 3*86400/4 -ymax ymax]);
                yticks(f.ax(3),-ymax:20:ymax);
                %             ylabel(f.ax(3),'Samples converted to 40-sps base')
                grid(f.ax(3),'on');
                box(f.ax(3),'on');
                hold(f.ax(3),'off');
                
                f.ax(4) = subplot(4,1,4,'align');
                hold(f.ax(4),'on');
                plot(f.ax(4),rstmathf(:,1),rstmathf(:,3)*40/spshf,'bs','MarkerSize',msize);
                plot(f.ax(4),rstmathf(:,1),rstmathf(:,2)*40/spshf,'ro','MarkerSize',msize);
                axis(f.ax(4),[3*86400/4 86400 -ymax ymax]);
                yticks(f.ax(4),-ymax:20:ymax);
                %             ylabel(f.ax(4),'Samples converted to 40-sps base')
                grid(f.ax(4),'on');
                box(f.ax(4),'on');
                xlabel(f.ax(4), 'Sec');
                hold(f.ax(4),'off');
                
                orient(f.fig,'landscape');    % is used to set up the paper orientation of a Figure or Model window for printing
                if fcflag==0    % no filtering correction & no conversion
                    print(f.fig,'-depsc',[rstpath,'/ANA/',fam,'_hf_',YEAR,JDAY,'.eps']);
                else    % filtering correction & no conversion
                    print(f.fig,'-depsc',[rstpath,'/ANA/',fam,'_hf_fc',YEAR,JDAY,'.eps']);  % filtering corrected
                end
                close(f.fig);
            end
            
            %%% for lf, re-assign the results to each day in each fam
            rstmatlf = lfallnew(lfallnew(:,1)==str2double(fam)&lfallnew(:,2)==date, 3:end);
            nin_new = size(rstmatlf,1);
            if isempty(rstmatlf)
                fprintf('Date %s in fam %s is empty! \n',date, fam);
            else 
                rstoldlf = lfallold(lfallold(:,1)==str2double(fam)&lfallold(:,2)==date, 3:end);
                %%% obtain the indexes of remaining detections, in order to
                %%% preserve the corresponding seismic traces only
                [~, iold, ~] = intersect(rstoldlf, rstmatlf, 'row', 'stable');
                
                %%% Write results into files
                IDENTIF=[YEAR,'.',JDAY,'.',fam,'.lo',num2str(lofflf),'.cc',num2str(xclf),'.',...
                    int2str(npolf),int2str(npalf),'.ms', int2str(msftlf)]
                if isfile(flf1) && fcflag==0 % save files when only removing double counting is used
                    %%% save detections
                    if iuplf ==1
                        fid = fopen(strcat(rstpath,'/MAPS/mapallrstNoDou',IDENTIF,'_',num2str(lolf),...
                              '-',num2str(hilf),'_',num2str(wlenlf/spslf),'s',num2str(spslf),'sps',...
                              '4nsta'),'w+');
                    else
                        fid = fopen(strcat(rstpath,'/MAPS/mapallrstUpNoDou',IDENTIF,'_',num2str(lolf),...
                              '-',num2str(hilf),'_',num2str(wlenlf/spslf),'s',num2str(spslf),'sps',...
                              '4nsta'),'w+');
                    end
                    fprintf(fid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %10.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f \n',rstmatlf');
                    fclose(fid);
                    %%% read old seismic traces and save new ones
                    traceoldlf = load(strcat(rstpath, '/MAPS/traceall_',IDENTIF,'_',num2str(lolf),'-',...
                        num2str(hilf),'_',num2str(wlenlf/spslf),'s',num2str(spslf),...
                        'sps','4nsta'));
                    tracenewlf = zeros(nin_new*wlenlf, size(traceoldlf,2));
                    for i = 1: nin_new
                        tracenewlf((i-1)*wlenlf+1: i*wlenlf,:)= traceoldlf((iold(i)-1)*wlenlf+1 : ...
                            iold(i)*wlenlf, :);
                    end
                    if iuplf ==1
                        fid = fopen(strcat(rstpath, '/MAPS/traceallNoDou_',IDENTIF,'_',num2str(lolf),...
                              '-',num2str(hilf),'_',num2str(wlenlf/spslf),'s',num2str(spslf),'sps',...
                              '4nsta'),'w+');
                    else
                        fid = fopen(strcat(rstpath, '/MAPS/traceallUpNoDou_',IDENTIF,'_',num2str(lolf),...
                              '-',num2str(hilf),'_',num2str(wlenlf/spslf),'s',num2str(spslf),'sps',...
                              '4nsta'),'w+');
                    end
                    fprintf(fid,'%.4f %.6f %.6f %.6f %.6f %.6f %.6f %.6f \n',tracenewlf');
                    fclose(fid);
                    
                elseif isfile(flf2) && fcflag==0 % save files when only removing double counting is used
                    if iuplf ==1
                        fid = fopen(strcat(rstpath,'/MAPS/mapallrstNoDou',IDENTIF,'_',num2str(lolf),...
                              '-',num2str(hilf),'_',num2str(wlenlf/spslf),'s',num2str(spslf),'sps',...
                              '3nsta'),'w+');
                    else
                        fid = fopen(strcat(rstpath,'/MAPS/mapallrstUpNoDou',IDENTIF,'_',num2str(lolf),...
                              '-',num2str(hilf),'_',num2str(wlenlf/spslf),'s',num2str(spslf),'sps',...
                              '3nsta'),'w+');
                    end
                    fprintf(fid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %10.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %d %d %d %.4f %.4f %.4f %d %d %d %.3f %.3f %.3f \n',rstmatlf');
                    fclose(fid);
                    %%% read old seismic traces and save new ones
                    traceoldlf = load(strcat(rstpath, '/MAPS/traceall_',IDENTIF,'_',num2str(lolf),...
                        '-',num2str(hilf),'_',num2str(wlenlf/spslf),'s',num2str(spslf),...
                        'sps','3nsta'));
                    tracenewlf = zeros(nin_new*wlenlf, size(traceoldlf,2));
                    for i = 1: nin_new
                        tracenewlf((i-1)*wlenlf+1: i*wlenlf,:)= traceoldlf((iold(i)-1)*wlenlf+1 : ...
                            iold(i)*wlenlf, :);
                    end
                    if iuplf ==1
                        fid = fopen(strcat(rstpath, '/MAPS/traceallNoDou_',IDENTIF,'_',num2str(lolf),...
                              '-',num2str(hilf),'_',num2str(wlenlf/spslf),'s',num2str(spslf),'sps',...
                              '3nsta'),'w+');
                    else
                        fid = fopen(strcat(rstpath, '/MAPS/traceallUpNoDou_',IDENTIF,'_',num2str(lolf),...
                              '-',num2str(hilf),'_',num2str(wlenlf/spslf),'s',num2str(spslf),'sps',...
                              '3nsta'),'w+');
                    end
                    fprintf(fid,'%.4f %.6f %.6f %.6f %.6f %.6f %.6f \n',tracenewlf');
                    fclose(fid);
                end
                
                %%% plot the lf offset variation against time
                f2.fig = figure(id);
                set(f2.fig,'Position',[scrsz(3)/10 scrsz(4)/10 2*scrsz(3)/5 1*scrsz(4)/2]);
                ymax = 40;
                msize = 1;
                f2.ax(1) = subplot(4,1,1,'align');
                hold(f2.ax(1),'on');
                plot(f2.ax(1),rstmatlf(:,1),rstmatlf(:,3)*40/spslf,'bs','MarkerSize',msize);
                plot(f2.ax(1),rstmatlf(:,1),rstmatlf(:,2)*40/spslf,'ro','MarkerSize',msize);
                axis(f2.ax(1),[0 86400/4 -ymax ymax]);
                yticks(f2.ax(1),-ymax:20:ymax);
                %             ylabel(f2.ax(1),'Samples converted to 40-sps base')
                title(f2.ax(1),strcat(fam,'-',YEAR,JDAY,'-',num2str(lolf),'-',num2str(hilf)));
                grid(f2.ax(1),'on');
                box(f2.ax(1),'on');
                hold(f2.ax(1),'off');
                
                f2.ax(2) = subplot(4,1,2,'align');
                hold(f2.ax(2),'on');
                plot(f2.ax(2),rstmatlf(:,1),rstmatlf(:,3)*40/spslf,'bs','MarkerSize',msize);
                plot(f2.ax(2),rstmatlf(:,1),rstmatlf(:,2)*40/spslf,'ro','MarkerSize',msize);
                axis(f2.ax(2),[86400/4 86400/2 -ymax ymax]);
                yticks(f2.ax(2),-ymax:20:ymax);
                ylabel(f2.ax(2),'Samples converted to 40-sps base')
                grid(f2.ax(2),'on');
                box(f2.ax(2),'on');
                hold(f2.ax(2),'off');
                
                f2.ax(3) = subplot(4,1,3,'align');
                hold(f2.ax(3),'on');
                plot(f2.ax(3),rstmatlf(:,1),rstmatlf(:,3)*40/spslf,'bs','MarkerSize',msize);
                plot(f2.ax(3),rstmatlf(:,1),rstmatlf(:,2)*40/spslf,'ro','MarkerSize',msize);
                axis(f2.ax(3),[86400/2 3*86400/4 -ymax ymax]);
                yticks(f2.ax(3),-ymax:20:ymax);
                %             ylabel(f2.ax(3),'Samples converted to 40-sps base')
                grid(f2.ax(3),'on');
                box(f2.ax(3),'on');
                hold(f2.ax(3),'off');
                
                f2.ax(4) = subplot(4,1,4,'align');
                hold(f2.ax(4),'on');
                plot(f2.ax(4),rstmatlf(:,1),rstmatlf(:,3)*40/spslf,'bs','MarkerSize',msize);
                plot(f2.ax(4),rstmatlf(:,1),rstmatlf(:,2)*40/spslf,'ro','MarkerSize',msize);
                axis(f2.ax(4),[3*86400/4 86400 -ymax ymax]);
                yticks(f2.ax(4),-ymax:20:ymax);
                %             ylabel(f2.ax(4),'Samples converted to 40-sps base');
                grid(f2.ax(4),'on');
                box(f2.ax(4),'on');
                xlabel(f2.ax(4), 'Sec');
                hold(f2.ax(4),'off');
                
                orient(f2.fig,'landscape');    % is used to set up the paper orientation of a Figure or Model window for printing
                if fcflag==0    % no filtering correction & no conversion
                    print(f2.fig,'-depsc',[rstpath,'/ANA/',fam,'_lf_',YEAR,JDAY,'.eps']);
                else    % filtering correction & no conversion
                    print(f2.fig,'-depsc',[rstpath,'/ANA/',fam,'_lf_fc',YEAR,JDAY,'.eps']);  % filtering corrected
                end
                close(f2.fig);
            end
            
        end     % end for day
    end % end for fam
    
%%% IF the convflag = 1, then all fams can be plotted together, so that plots are separated only by
%%% the day
else
    disp('Start ploting and saving figures...')
    % sort according to date,time
    hfallnew = sortrows(hfallnew, [2,3]);
    lfallnew = sortrows(lfallnew, [2,3]);
    
    %%% for HF
    dateall = unique(hfallnew(:,2));
    for id = 1: length(dateall)
        hfday = hfallnew(hfallnew(:,2) == dateall(id), :);
        date = num2str(dateall(id));
        YEAR = date(1:4);
        JDAY = date(5:7);
        %%% plot the hf offset variation against time
        f3.fig = figure(id);
        set(f3.fig,'Position',[scrsz(3)/10 scrsz(4)/10 2*scrsz(3)/5 1*scrsz(4)/2]);
        ymax = 60;
        msize = 1;
        f3.ax(1) = subplot(4,1,1,'align');
        hold(f3.ax(1),'on');
        plot(f3.ax(1),hfday(:,3),hfday(:,5)*40/spshf,'bs','MarkerSize',msize);
        plot(f3.ax(1),hfday(:,3),hfday(:,4)*40/spshf,'ro','MarkerSize',msize);
        axis(f3.ax(1),[0 86400/4 -ymax ymax]);
        yticks(f3.ax(1),-ymax:20:ymax);
        %             ylabel(f2.ax(1),'Samples converted to 40-sps base')
        title(f3.ax(1),strcat('Origin converted to',{' '},reffam,'-',YEAR,JDAY,'-',num2str(lohf),...
              '-',num2str(hihf)));
        grid(f3.ax(1),'on');
        box(f3.ax(1),'on');
        hold(f3.ax(1),'off');
        
        f3.ax(2) = subplot(4,1,2,'align');
        hold(f3.ax(2),'on');
        plot(f3.ax(2),hfday(:,3),hfday(:,5)*40/spshf,'bs','MarkerSize',msize);
        plot(f3.ax(2),hfday(:,3),hfday(:,4)*40/spshf,'ro','MarkerSize',msize);
        axis(f3.ax(2),[86400/4 86400/2 -ymax ymax]);
        yticks(f3.ax(2),-ymax:20:ymax);
        ylabel(f3.ax(2),'Samples converted to 40-sps base')
        grid(f3.ax(2),'on');
        box(f3.ax(2),'on');
        hold(f3.ax(2),'off');
        
        f3.ax(3) = subplot(4,1,3,'align');
        hold(f3.ax(3),'on');
        plot(f3.ax(3),hfday(:,3),hfday(:,5)*40/spshf,'bs','MarkerSize',msize);
        plot(f3.ax(3),hfday(:,3),hfday(:,4)*40/spshf,'ro','MarkerSize',msize);
        axis(f3.ax(3),[86400/2 3*86400/4 -ymax ymax]);
        yticks(f3.ax(3),-ymax:20:ymax);
        %             ylabel(f2.ax(3),'Samples converted to 40-sps base')
        grid(f3.ax(3),'on');
        box(f3.ax(3),'on');
        hold(f3.ax(3),'off');
        
        f3.ax(4) = subplot(4,1,4,'align');
        hold(f3.ax(4),'on');
        plot(f3.ax(4),hfday(:,3),hfday(:,5)*40/spshf,'bs','MarkerSize',msize);
        plot(f3.ax(4),hfday(:,3),hfday(:,4)*40/spshf,'ro','MarkerSize',msize);
        axis(f3.ax(4),[3*86400/4 86400 -ymax ymax]);
        yticks(f3.ax(4),-ymax:20:ymax);
        %             ylabel(f2.ax(4),'Samples converted to 40-sps base');
        grid(f3.ax(4),'on');
        box(f3.ax(4),'on');
        xlabel(f3.ax(4), 'Sec');
        hold(f3.ax(4),'off');
        
        orient(f3.fig,'landscape');    % is used to set up the paper orientation of a Figure or Model window for printing
        if fcflag==0    % no filtering correction & no conversion
            print(f3.fig,'-depsc',[rstpath,'/ANA/','allfam_hf_conv_',YEAR,JDAY,'.eps']);
        else    % filtering correction & no conversion
            print(f3.fig,'-depsc',[rstpath,'/ANA/','allfam_hf_fc_conv_',YEAR,JDAY,'.eps']);  % filtering corrected
        end
        close(f3.fig);
    end
    
    %%% for LF
    dateall = unique(lfallnew(:,2));
    for id = 1: length(dateall)
        lfday = lfallnew(lfallnew(:,2) == dateall(id), :);
        date = num2str(dateall(id));
        YEAR = date(1:4);
        JDAY = date(5:7);
        %%% plot the lf offset variation against time
        f4.fig = figure(id);
        set(f4.fig,'Position',[scrsz(3)/10 scrsz(4)/10 2*scrsz(3)/5 1*scrsz(4)/2]);
        ymax = 60;
        msize = 1;
        f4.ax(1) = subplot(4,1,1,'align');
        hold(f4.ax(1),'on');
        plot(f4.ax(1),lfday(:,3),lfday(:,5)*40/spslf,'bs','MarkerSize',msize);
        plot(f4.ax(1),lfday(:,3),lfday(:,4)*40/spslf,'ro','MarkerSize',msize);
        axis(f4.ax(1),[0 86400/4 -ymax ymax]);
        yticks(f4.ax(1),-ymax:20:ymax);
        %             ylabel(f2.ax(1),'Samples converted to 40-sps base')
        title(f4.ax(1),strcat('Origin converted to',{' '},reffam,'-',YEAR,JDAY,'-',num2str(lolf),...
              '-',num2str(hilf)));
        grid(f4.ax(1),'on');
        box(f4.ax(1),'on');
        hold(f4.ax(1),'off');
        
        f4.ax(2) = subplot(4,1,2,'align');
        hold(f4.ax(2),'on');
        plot(f4.ax(2),lfday(:,3),lfday(:,5)*40/spslf,'bs','MarkerSize',msize);
        plot(f4.ax(2),lfday(:,3),lfday(:,4)*40/spslf,'ro','MarkerSize',msize);
        axis(f4.ax(2),[86400/4 86400/2 -ymax ymax]);
        yticks(f4.ax(2),-ymax:20:ymax);
        ylabel(f4.ax(2),'Samples converted to 40-sps base')
        grid(f4.ax(2),'on');
        box(f4.ax(2),'on');
        hold(f4.ax(2),'off');
        
        f4.ax(3) = subplot(4,1,3,'align');
        hold(f4.ax(3),'on');
        plot(f4.ax(3),lfday(:,3),lfday(:,5)*40/spslf,'bs','MarkerSize',msize);
        plot(f4.ax(3),lfday(:,3),lfday(:,4)*40/spslf,'ro','MarkerSize',msize);
        axis(f4.ax(3),[86400/2 3*86400/4 -ymax ymax]);
        yticks(f4.ax(3),-ymax:20:ymax);
        %             ylabel(f2.ax(3),'Samples converted to 40-sps base')
        grid(f4.ax(3),'on');
        box(f4.ax(3),'on');
        hold(f4.ax(3),'off');
        
        f4.ax(4) = subplot(4,1,4,'align');
        hold(f4.ax(4),'on');
        plot(f4.ax(4),lfday(:,3),lfday(:,5)*40/spslf,'bs','MarkerSize',msize);
        plot(f4.ax(4),lfday(:,3),lfday(:,4)*40/spslf,'ro','MarkerSize',msize);
        axis(f4.ax(4),[3*86400/4 86400 -ymax ymax]);
        yticks(f4.ax(4),-ymax:20:ymax);
        %             ylabel(f2.ax(4),'Samples converted to 40-sps base');
        grid(f4.ax(4),'on');
        box(f4.ax(4),'on');
        xlabel(f4.ax(4), 'Sec');
        hold(f4.ax(4),'off');
        
        orient(f4.fig,'landscape');    % is used to set up the paper orientation of a Figure or Model window for printing
        if fcflag==0    % no filtering correction & no conversion
            print(f4.fig,'-depsc',[rstpath,'/ANA/','allfam_lf_conv_',YEAR,JDAY,'.eps']);
        else    % filtering correction & no conversion
            print(f4.fig,'-depsc',[rstpath,'/ANA/','allfam_lf_fc_conv_',YEAR,JDAY,'.eps']);  % filtering corrected
        end
        close(f4.fig);
    end
end


















