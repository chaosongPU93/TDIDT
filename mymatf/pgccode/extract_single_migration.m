%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is to extract the hf and lf detection results inside a single
% and typical migration, ready for inverting to map locations
%
% Read 'timeori_combined', which is original detected result with 3-sta
% trio, including the times
%
% Feature:
%   1. read lf & hf result file, also should be corrected for filtering
%   effect
%   2. call functions to select data segment according to the timing query
%   3. save it to file
%
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/07/29
% Last modified date:   2019/07/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio

format short e   % Set the format to 5-digit floating point
clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

scrsz=get(0,'ScreenSize');

workpath = getenv('ALLAN');
%%% choose the data path here 
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/PGCtrio');

fam='002';     % family number

stas=['PGC  '         % NOTICE the blank spaces, stas have 5 columns
      'SSIB '
      'SILB '];
  
%% Important parameters same as that during detections    
nsta=size(stas,1);         %  number of stations

%%% IMPORTANT, NEED to change sometimes
%Basics of the cross-correlation:  Window length, number of windows, filter parameters, etc.

%%%% HF %%%%
spshf=40;     % samples per second
winlensechf=4;
winoffsechf=1;        % window offset in sec, which is the step of a moving window
winlenhf=winlensechf*spshf;      % length in smaples
hihf=6.5;
lohf=1.25;
npo=2;     % poles, passes of filters
npa=2;
cyclskiphf = 0;
mshifthf=29+cyclskiphf; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
loopoffmaxhf=2.1; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
xcmaxAVEnminhf=0.44; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?
cncthf = 0.5;

%%%% LF %%%%%
spslf=20;
winlenseclf=12;       % CHECK in detection script for what were used in first-year report
winoffseclf=3;        % window offset in sec, which is the step of a moving window
winlenlf=winlenseclf*spslf;      % length in smaples
hilf=1.25;
lolf=0.5;
npo=2;     % poles, passes of filters
npa=2;
cyclskiplf = 20;
mshiftlf=19+cyclskiplf; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
loopoffmaxlf=2; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
xcmaxAVEnminlf=0.6; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?
cnctlf = 1.1;

% stasnew=['LZB  '
%          'TWKB '
%          'MGCB '];
stasnew=['LZB  '
         'TWKB '
         'MGCB '
         'KLNB '];
nstanew=size(stasnew,1);
daylen = 86400;    % length of sec of a day

PREFIXhf = strcat(fam,'.loff',num2str(loopoffmaxhf),'.ccmin',num2str(xcmaxAVEnminhf),'.nponpa',int2str(npo),int2str(npa),'.ms', int2str(mshifthf));
PREFIXlf = strcat(fam,'.loff',num2str(loopoffmaxlf),'.ccmin',num2str(xcmaxAVEnminlf),'.nponpa',int2str(npo),int2str(npa),'.ms', int2str(mshiftlf));

%% obtain the hf&lf detection inside a migration as query 

%%% read hf detection file & time file  
fname = strcat(rstpath, '/MAPS/timeori_combined_',PREFIXhf,'_',num2str(lohf),'-',num2str(hihf),'_',num2str(winlenhf/spshf),'s',num2str(spshf),'sps', num2str(nstanew), 'newsta');
% STA23time structure: off12 off13 off12sec off13sec year-day timing-of-strongest-arrival timing-of-center-of-detection-window
hftimeori = load(fname);
hftimeori(:,3:4) = vpa(hftimeori(:,3:4),4);  % in case of the precision problem of floats 
hftimeori(:,6) = vpa(hftimeori(:,6),4);

%%% read lf detection file & time file  
fname = strcat(rstpath, '/MAPS/timeori_',PREFIXlf,'_',num2str(lolf),'-',num2str(hilf),'_',num2str(winlenlf/spslf),'s',num2str(spslf),'sps', num2str(nstanew), 'newsta');
% STA23time structure: off12 off13 off12sec off13sec year-day timing-of-strongest-arrival timing-of-center-of-detection-window
lftimeori = load(fname);
lftimeori(:,3:4) = vpa(lftimeori(:,3:4),4);  % in case of the precision problem of floats 
lftimeori(:,6) = vpa(lftimeori(:,6),4); % KEEP in mind that the order of time in hf and lf maybe different

%%% for example, [hfmig, lfmig] = GetDetectionInMigration(fam,date,sttime,edtime,hftimeori,lftimeori)
[hfmig1, lfmig1] = GetDetectionInMigration(fam,2005255,3.45e+4,3.65e+4,hftimeori,lftimeori);


%%% save file 
PREFIX = strcat(fam,'_',num2str(2005255),'_',num2str(3.45e+4),'-',num2str(3.65e+4));

fid = fopen(strcat(rstpath, '/MAPS/migration_hf','_',PREFIX),'w+');
fprintf(fid,'%d %d %.4f %.4f %d %.4f %.1f \n',hfmig1');
fclose(fid);

fid = fopen(strcat(rstpath, '/MAPS/migration_lf','_',PREFIX),'w+');
fprintf(fid,'%d %d %.4f %.4f %d %.4f %.1f \n',lfmig1');
fclose(fid);





















