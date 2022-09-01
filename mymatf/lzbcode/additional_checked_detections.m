%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the script to analyse the inversion results from hypoinverse, set
% a cutoff distance, e.g. 8 km. from the center of lfe fam to their own 
% detections. Perserve the ones within this distance only.
%
% Although the distance cutoff is applied to count and time data, but the 
% saved results should be the same anyway, in other words, because the count
% data does not have the time information and have duplicates from diff fam,
% so only the processed data is worthwhile to save, for the next step,
% merge_fam_afthypo.m 
%   
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/11/25
% Last modified date:   2019/11/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all
clc

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

scrsz=get(0,'ScreenSize');

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/lzbrst');
datapath = '/home/data2/chaosong/matlab/allan/data-no-resp';

FLAG = 'LZB';

winlenhf = 4;

winlenlf = 16;
lofflf = 4;
ccminlf = 0.35;


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
           -123.838500 48.544833 35.6600];

relacont = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),loccont(2,1),loccont(2,2));
relacont(:,1) = dx;
relacont(:,2) = dy;       

% convert absolute loc to relative loc to its own lfe fam
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
        
%% for detections in time
SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));

dcutflag = 1;

if dcutflag 
    hffname = strcat(rstpath, '/evtloc.allfam.dcut.',SUFFIXhf);
    lffname = strcat(rstpath, '/evtloc.allfam.dcut.',SUFFIXlf);
else
    hffname = strcat(rstpath, '/evtloc.allfam.nodcut.',SUFFIXhf);
    lffname = strcat(rstpath, '/evtloc.allfam.nodcut.',SUFFIXlf);
end

hftime = load(hffname);
lftime = load(lffname);
% 25 cols, format is:
%   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date 
%   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum


%% read the alligned seismograms of all detections and get the amplitude info
%%% NOTE: it is impossible to track the order of detections since they have been sorted somewhere
%%% thus the only way to link them is the fam number, day and time
timoffrot= [2003 060;
            2003 061;
            2003 062; % dates of the year, two columns: year Julian_day
            2003 063;
            2003 064;
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

tracepath = '/home/data2/chaosong/matlab/allan/data-no-resp/LZBtrio';


% %% read the original detection results to get the max energy in dtmin per detection window
% %%% NOTE: it is impossible to track the order of detections since they have been sorted somewhere
% %%% thus the only way to link them is the fam number, day and time
% energyhf = zeros(size(hftime,1),2);   % hf, amplitude of detections for all days of all fams
% for i=1: size(nfampool,1)
%     fam = nfampool(i,:);
%     for nd=1:length(timoffrot(:,1))      % num of rows, also num of days
%         year=timoffrot(nd,1);
%         YEAR=int2str(year);
%         jday=timoffrot(nd,2);
%         if jday <= 9
%             JDAY=['00',int2str(jday)];
%         elseif jday<= 99
%             JDAY=['0',int2str(jday)];
%         else
%             JDAY=int2str(jday);
%         end
%         IDENTIF=[YEAR,'.',JDAY,'.',fam,'.lo1.5.cc0.4.22.ms14_1.25-6.5_4s40sps4nsta'];
%         fname = strcat(tracepath, '/MAPS/mapallrst_up',IDENTIF);
%         if isfile(fname)
%             allrst = load(fname);
%             ndet = size(allrst,1);
%             for j = 1: ndet
%                 midtime = allrst(j,1);
%                 objind = find(hftime(:,end)==str2double(fam) & ...
%                               hftime(:,13)==str2double(strcat(YEAR,JDAY)) & ...
%                               abs(hftime(:,15)-midtime)< 0.01);
%                 if isempty(objind) || size(objind,1) > 1
%                     disp(midtime)
%                     break
%                 end
%                 energyhf(objind,1) = allrst(j,6);   % max energy in the dtmin window
%                 energyhf(objind,2) = allrst(j,8);   % max energy in the dtmin window normalized by energy of entire 4s win
%             end
%         else
%             continue
%         end
%     end
% end
% disp('Energy from trio stas linked to hf detections done');
% 
% 
% energylf = zeros(size(lftime,1),2);   % hf, amplitude of detections for all days of all fams
% for i=1: size(nfampool,1)
%     fam = nfampool(i,:);
%     for nd=1:length(timoffrot(:,1))      % num of rows, also num of days
%         year=timoffrot(nd,1);
%         YEAR=int2str(year);
%         jday=timoffrot(nd,2);
%         if jday <= 9
%             JDAY=['00',int2str(jday)];
%         elseif jday<= 99
%             JDAY=['0',int2str(jday)];
%         else
%             JDAY=int2str(jday);
%         end
%         IDENTIF=[YEAR,'.',JDAY,'.',fam,'.lo4.cc0.35.22.ms34_0.5-1.25_16s20sps4nsta'];
%         fname = strcat(tracepath, '/MAPS/mapallrst',IDENTIF);
%         if isfile(fname)
%             allrst = load(fname);
%             ndet = size(allrst,1);
%             for j = 1: ndet
%                 midtime = allrst(j,1);
%                 objind = find(lftime(:,end)==str2double(fam) & ...
%                               lftime(:,13)==str2double(strcat(YEAR,JDAY)) & ...
%                               abs(lftime(:,15)-midtime)< 0.5);
%                 if isempty(objind) || size(objind,1) > 1
%                     disp(midtime)
%                     break
%                 end
%                 energylf(objind,1) = allrst(j,6);   % max energy in the dtmin window
%                 energylf(objind,2) = allrst(j,8);   % max energy in the dtmin window normalized by energy of entire 4s win            end
%             end
%         else
%             continue
%         end
%     end
% end
% disp('Energy from trio stas linked to lf detections done');
% 
% 
% %% save the energy info
% %%% NOTE: this energy file shares the SAME order as 'evtloc.allfam.dcut.' or 'evtloc.allfam.nodcut.'
% %%% the file is for future usage
% ind = find(energyhf(:,1)~=0);
% energynewhf = energyhf(ind,:);
% 
% ind = find(energylf(:,1)~=0);
% energynewlf = energylf(ind,:);
% 
% % format is the same as above
% if dcutflag
%     fidhf = fopen(strcat(rstpath, '/energy.allfam.dcut.',SUFFIXhf),'w+');
%     fidlf = fopen(strcat(rstpath, '/energy.allfam.dcut.',SUFFIXlf),'w+');
% else
%     fidhf = fopen(strcat(rstpath, '/energy.allfam.nodcut.',SUFFIXhf),'w+');
%     fidlf = fopen(strcat(rstpath, '/energy.allfam.nodcut.',SUFFIXlf),'w+');
% end
% 
% fprintf(fidhf,'%10.3e %10.3f \n',energynewhf');
% fclose(fidhf);
% 
% fprintf(fidlf,'%10.3e %10.3f \n',energynewlf');
% fclose(fidlf);


%% read the detections all fam that passed the check from additional stations
%%% in fact, there can be several additional stations, it seems that KLNB and PGC have a better
%%% signal-to-noise ratio, but here get all, since they have been calculated from
%%% plotaddcheck_LZB_allfam.m
stasnew=['PGC  '
         'SSIB '
         'SILB '
         'KLNB '];
nstasnew = size(stasnew,1);   

addrsthf = [];
addrstlf = [];
for ista = 1: nstasnew
    % hf
    for i=1: size(nfampool,1)
        fam = nfampool(i,:);
        
        IDENTIF=[fam,'.up.lo1.5.ccm0.4.22.ms14_1.25-6.5_4s40sps4nsta'];
        fname = strcat(tracepath, '/MAPS/timeadd_',stasnew(ista,:),'_',IDENTIF);
        if isfile(fname)
            addrstfam = load(fname);
            %%% 9 cols
            %%% off124 | off134 | off124sec | off134sec | year-day | timing-of-strongest-arrival |
            %%% timing-of-center-of-detection-window | cc14 | famnum
            addrsthf = [addrsthf; addrstfam];
        else
            continue
        end
    end
    
    addrstnewhf = [];

    for i = 1: size(addrsthf,1)
        objind = find(hftime(:,end)==addrsthf(i,end) & ...
            hftime(:,9)==addrsthf(i,1) & ...
            hftime(:,10)==addrsthf(i,2) & ...
            hftime(:,13)==addrsthf(i,5) & ...
            hftime(:,15)==addrsthf(i,7));
        if ~isempty(objind)
            tmp = [hftime(objind, 1: 16) hftime(objind, 21: 24) addrsthf(i,8:9)];
            addrstnewhf = [addrstnewhf; tmp];
        end
        % 22 cols, format is:
        %   E(043) N(043) E(own) N(own) dist(own) lon lat dep off124 off134 off124sec off134sec date 
        %   main_arrival_time  cent_of_win  avecc ampmax ampmin eng normeng ccadd4 famnum

    end
    
    if dcutflag
        fid = fopen(strcat(rstpath, '/evtloc.',stasnew(ista,:),'check.allfam.dcut.',SUFFIXhf),'w+');
    else
        fid = fopen(strcat(rstpath, '/evtloc.',stasnew(ista,:),'check.allfam.nodcut.',SUFFIXhf),'w+');
    end
    fprintf(fid,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %d %.4f %.1f %.4f %.4e %.4e %10.3e %10.3f %.4f %d \n',...
            addrstnewhf');
    fclose(fid);
    
    disp(['hf detections checked by ',stasnew(ista,:),'linked and saved']);
    
    
    % lf
    for i=1: size(nfampool,1)
        fam = nfampool(i,:);
        
        IDENTIF=[fam,'.lo4.ccm0.35.22.ms34_0.5-1.25_16s20sps4nsta'];
        fname = strcat(tracepath, '/MAPS/timeadd_',stasnew(ista,:),'_',IDENTIF);
        if isfile(fname)
            addrstfam = load(fname);
            %%% 9 cols
            %%% off124 | off134 | off124sec | off134sec | year-day | timing-of-strongest-arrival |
            %%% timing-of-center-of-detection-window | cc14 | famnum
            addrstlf = [addrstlf; addrstfam];
        else
            continue
        end
    end
    
    addrstnewlf = [];
    for i = 1: size(addrstlf,1)
        objind = find(lftime(:,end)==addrstlf(i,end) & ...
                      lftime(:,9)==addrstlf(i,1) & ...
                      lftime(:,10)==addrstlf(i,2) & ...
                      lftime(:,13)==addrstlf(i,5) & ...
                      lftime(:,15)==addrstlf(i,7));
%         tmp = [lftime(objind, 1: 16) lftime(objind, 21: 22) energynewlf(objind,:) addrstlf(i,8:9)];
        if ~isempty(objind)
            tmp = [lftime(objind, 1: 16) lftime(objind, 21: 24) addrstlf(i,8:9)];
            addrstnewlf = [addrstnewlf; tmp];
        end
        % 22 cols, format is:
        %   E(043) N(043) E(own) N(own) dist(own) lon lat dep off124 off134 off124sec off134sec date 
        %   main_arrival_time  cent_of_win  avecc ampmax ampmin eng normeng ccadd4 famnum

    end
    
    if dcutflag
        fid = fopen(strcat(rstpath, '/evtloc.',stasnew(ista,:),'check.allfam.dcut.',SUFFIXlf),'w+');
    else
        fid = fopen(strcat(rstpath, '/evtloc.',stasnew(ista,:),'check.allfam.nodcut.',SUFFIXlf),'w+');
    end
    fprintf(fid,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %d %.4f %.1f %.4f %.4e %.4e %10.3e %10.3f %.4f %d \n',...
            addrstnewlf');
    fclose(fid);
    
    disp(['lf detections checked by ',stasnew(ista,:),'linked and saved']);
    
end


%% read the mean offset of the 4th stations checks
%%% this is different from above, which contains all detections that passed the additonal station
%%% check. Here, we wonder, at each unique detected location, the mean time offsets of multiple 
%%% detections between 4th and 1st station at this location that passed the checked criteria, ie.
%%% the standard derivation, number of detections in the vincinity, so they don't contain time
%%% information, but only the location (12off and 13off), mean and std offset of 4th station, counts
%%%  
%%% in fact, there can be several additional stations, it seems that KLNB and PGC have a better
%%% signal-to-noise ratio, but here get all, since they have been calculated from
%%% plotaddcheck_LZB_allfam.m
stasnew=['PGC  '
         'SSIB '
         'SILB '
         'KLNB '];
nstasnew = size(stasnew,1);  

addrsthf = [];
addrstlf = [];
for ista = 4: nstasnew
    % hf
    for i=1: size(nfampool,1)
        fam = nfampool(i,:);
        
        IDENTIF=[fam,'.up.lo1.5.ccm0.4.22.ms14_1.25-6.5_4s40sps4nsta'];
        fname = strcat(tracepath, '/MAPS/countadd_',stasnew(ista,:),'_',IDENTIF);
        if isfile(fname)
            addrstfam = load(fname);
            %%% 8 cols
            %%% STA124off STA134off STA124offsec STA134offsec STA14meansec STA14stdsec STA14count ...
            %%% famarr
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% whether to convert sta14 offsets from different fams to a uniform frame
            %%% ALGORITHM is confirmed to be true!
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
            if ista == 1    % pgc
                conv = PERMROTS(1,4) - refperm(1,4);
            elseif ista == 2 % ssib
                conv = POLROTS(1,4) - refpol(1,4);
            elseif ista == 3 % silb    
                conv = POLROTS(2,4) - refpol(2,4);
            elseif ista == 4 % klnb
                conv = POLROTS(3,4) - refpol(3,4);
            end
            spshf = 40;
            spsrat = spshf/40;
            tmp = addrstfam(:,5)-conv*spsrat/spshf;   % get the new off based on the ref fam frame
            addrstfam = [addrstfam(:,1:5) tmp addrstfam(:,6:8)];    % become 9 cols
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            addrsthf = [addrsthf; addrstfam];
            
        else
            continue
        end
    end
    
    addrstnewhf = [];
    for i = 1: size(addrsthf,1)
        objind = find(hftime(:,end)==addrsthf(i,end) & ...
                      hftime(:,9)==addrsthf(i,1) & ...
                      hftime(:,10)==addrsthf(i,2));
        aaa = unique(hftime(objind,1));
        if length(aaa) >1
            disp(i)
            disp(objind)
            break
        elseif length(aaa) == 1
            tmp = [hftime(objind(1), 1: 12) addrsthf(i,5:end)];
            addrstnewhf = [addrstnewhf; tmp];
        end
        
        % 17 cols, format is:
        %   E(043) N(043) E(own) N(own) dist(own) lon lat dep off124 off134 off124sec off134sec
        %   mean14sec(own) mean14sec(043) std14sec counts famnum

    end
    
    addrsttmphf = sortrows(addrstnewhf, [1 2]);
    
    addrstnewhf = [];
    [locuniq,~,~] = unique(addrsttmphf(:,1:2),'rows','stable');
    for i = 1: size(locuniq,1)
        [~,ind,~] = intersect(addrsttmphf(:,1:2),locuniq(i,:),'rows','stable');
        [~,ind2] = max(addrsttmphf(ind,16));
        addrstnewhf = [addrstnewhf; addrsttmphf(ind(ind2),:)];
    end
    
    if dcutflag
        fid = fopen(strcat(rstpath, '/addstaoff.',stasnew(ista,:),'check.allfam.dcut.',SUFFIXhf),'w+');
    else
        fid = fopen(strcat(rstpath, '/addstaoff.',stasnew(ista,:),'check.allfam.nodcut.',SUFFIXhf),'w+');
    end
    fprintf(fid,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %d %d \n',...
            addrstnewhf');
    fclose(fid);
    
    disp(['hf mean offset at additional station checked by ',stasnew(ista,:),' saved']);
    
    f1.fig = figure;
    ax=gca;
    scatter(ax,addrstnewhf(:,1), addrstnewhf(:,2), 4, -addrstnewhf(:,14), 'filled','o');    %, 'MarkerEdgeColor', 'w')
    colormap(ax,jet)
    colorbar(ax);
    caxis(ax,[-1.5,1]);
    box on
    axis equal
    xlabel(ax,'N (km)')
    ylabel(ax,'E (km)')
    axis(ax,[-30 50 -30 40]);
    title(ax,strcat(strtrim(stasnew(ista, :))));
    if dcutflag
        print(f1.fig,'-dpdf',strcat(rstpath,'/addstaoff.',stasnew(ista,:),'check.allfam.dcut.',SUFFIXhf,'.pdf'));
    else
        print(f1.fig,'-dpdf',strcat(rstpath,'/addstaoff.',stasnew(ista,:),'check.allfam.nodcut.',SUFFIXhf,'.pdf'));
    end
    
    % lf
    for i=1: size(nfampool,1)
        fam = nfampool(i,:);
        
        IDENTIF=[fam,'.lo4.ccm0.35.22.ms34_0.5-1.25_16s20sps4nsta'];
        fname = strcat(tracepath, '/MAPS/countadd_',stasnew(ista,:),'_',IDENTIF);
        if isfile(fname)
            addrstfam = load(fname);
            %%% 8 cols
            %%% STA124off STA134off STA124offsec STA134offsec STA14meansec STA14stdsec STA14count ...
            %%% famarr
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% whether to convert sta14 offsets from different fams to a uniform frame
            %%% ALGORITHM is confirmed to be true!
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
            if ista == 1    % pgc
                conv = PERMROTS(1,4) - refperm(1,4);
            elseif ista == 2 % ssib
                conv = POLROTS(1,4) - refpol(1,4);
            elseif ista == 3 % silb    
                conv = POLROTS(2,4) - refpol(2,4);
            elseif ista == 4 % klnb
                conv = POLROTS(3,4) - refpol(3,4);
            end
            spslf = 20;
            spsrat = spslf/40;
            tmp = addrstfam(:,5)-conv*spsrat/spslf;   % get the new off based on the ref fam frame
            addrstfam = [addrstfam(:,1:5) tmp addrstfam(:,6:8)];    % become 9 cols
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            addrstlf = [addrstlf; addrstfam];
            
        else
            continue
        end
    end
    
    addrstnewlf = [];
    for i = 1: size(addrstlf,1)
        objind = find(lftime(:,end)==addrstlf(i,end) & ...
                      lftime(:,9)==addrstlf(i,1) & ...
                      lftime(:,10)==addrstlf(i,2), 1, 'first');
                          aaa = unique(hftime(objind,1));
        if length(aaa) >1
            disp(i)
            disp(objind)
            break
        elseif length(aaa) == 1
            tmp = [lftime(objind(1), 1: 12) addrstlf(i,5:end)];
            addrstnewlf = [addrstnewlf; tmp];
        end
        % 17 cols, format is:
        %   E(043) N(043) E(own) N(own) dist(own) lon lat dep off124 off134 off124sec off134sec
        %   mean14sec(own) mean14sec(043) std14sec counts famnum

    end
    
    addrsttmplf = sortrows(addrstnewlf, [1 2]);
    
    addrstnewlf = [];
    [locuniq,~,~] = unique(addrsttmplf(:,1:2),'rows','stable');
    for i = 1: size(locuniq,1)
        [~,ind,~] = intersect(addrsttmplf(:,1:2),locuniq(i,:),'rows','stable');
        [~,ind2] = max(addrsttmplf(ind,16));
        addrstnewlf = [addrstnewlf; addrsttmplf(ind(ind2),:)];
    end
    
    if dcutflag
        fid = fopen(strcat(rstpath, '/addstaoff.',stasnew(ista,:),'check.allfam.dcut.',SUFFIXlf),'w+');
    else
        fid = fopen(strcat(rstpath, '/addstaoff.',stasnew(ista,:),'check.allfam.nodcut.',SUFFIXlf),'w+');
    end
    fprintf(fid,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %d %d \n',...
            addrstnewlf');
    fclose(fid);
    
    disp(['lf mean offset at additional station checked by ',stasnew(ista,:),' saved']);
    

    f2.fig=figure;
    ax=gca;
    scatter(ax,addrstnewlf(:,1), addrstnewlf(:,2), 8, -addrstnewlf(:,14), 'filled','o');    %, 'MarkerEdgeColor', 'w')
    colormap(ax,jet)
    colorbar(ax);
    caxis(ax,[-1.5,1]);
    box on
    axis equal
    xlabel(ax,'N (km)')
    ylabel(ax,'E (km)')
    axis(ax,[-30 50 -30 40]);
    title(ax,strcat(strtrim(stasnew(ista, :))));
    if dcutflag
        print(f2.fig,'-dpdf',strcat(rstpath,'/addstaoff.',stasnew(ista,:),'check.allfam.dcut.',SUFFIXlf,'.pdf'));
    else
        print(f2.fig,'-dpdf',strcat(rstpath,'/addstaoff.',stasnew(ista,:),'check.allfam.nodcut.',SUFFIXlf,'.pdf'));
    end
    
    close all
end        
        
        

        
        
        
        
        
        
        
        
        
        
        