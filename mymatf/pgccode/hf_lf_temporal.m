% function hf_lf_temporal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to preliminarily analyse the temporal relationship between
% lf and hf detection results after adding several stations to check.
% Read 'timeori_combined', which is original detected result with 3-sta
% trio, including the times.
% You DO NOT have to run this if you don't need contemporaneous detections
%
% Feature:
%   1. read lf & hf result file, also should be corrected for filtering
%   effect
%   2. plot their temporal relationship
%   3. plot and writing files
%
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/06/15
% Last modified date:   2019/07/25
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
winlenseclf=16;       % CHECK in detection script for what were used in first-year report
winoffseclf=8;        % window offset in sec, which is the step of a moving window
winlenlf=winlenseclf*spslf;      % length in smaples
hilf=1.25;
lolf=0.5;
npo=2;     % poles, passes of filters
npa=2;
cyclskiplf = 20;
mshiftlf=19+cyclskiplf; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
loopoffmaxlf=4; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
xcmaxAVEnminlf=0.5; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?
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

%% get the overlapped detections before checking 
%%% read hf time file  
% NOTE! the following files are already corrected for filtering effect
fname = strcat(rstpath, '/MAPS/timeori_',PREFIXhf,'_',num2str(lohf),'-',num2str(hihf),'_',num2str(winlenhf/spshf),'s',num2str(spshf),'sps', num2str(nstanew), 'newsta');
% STA23time structure: off12 off13 off12sec off13sec year-day timing-of-strongest-arrival timing-of-center-of-detection-window
% 'cor' means filtering effect corrected 
hftimeoricor = load(fname);
hftimeoricor(:,3:4) = vpa(hftimeoricor(:,3:4),4);  % in case of the precision problem of floats 
hftimeoricor(:,6) = vpa(hftimeoricor(:,6),4);

%%% read lf time file  
fname = strcat(rstpath, '/MAPS/timeori_',PREFIXlf,'_',num2str(lolf),'-',num2str(hilf),'_',num2str(winlenlf/spslf),'s',num2str(spslf),'sps', num2str(nstanew), 'newsta');
% STA23time structure: off12 off13 off12sec off13sec year-day timing-of-strongest-arrival timing-of-center-of-detection-window
lftimeoricor = load(fname);
lftimeoricor(:,3:4) = vpa(lftimeoricor(:,3:4),4);  % in case of the precision problem of floats 
lftimeoricor(:,6) = vpa(lftimeoricor(:,6),4); % KEEP in mind that the order of time in hf and lf maybe different

% call function DetectOverlap to get overlapped detections in this 2
% frequencies
%   hfalloricor: all hf detections inside the lf detecting window
%   hfallorimedcor: median of the above in one lf window
%   lfalloricor: all lf detections that have contemporal hf part
%   dateoricor: dates that have contemporal part
[hfalloricor,hfallorimedcor,lfalloricor,dateoricor]=DetectOverlap(hftimeoricor,lftimeoricor,winlensechf,winlenseclf); 


%%% load the template filtering effect result
% PREFIX = strcat(fam,'.lf',num2str(lolf),'-',num2str(hilf),'.hf',num2str(lohf),'-',num2str(hihf),...
%                     '.hfwin',num2str(winlensechf),'.lfwin',num2str(winlenseclf),'.resp-corrected');
% fname = strcat(rstpath, '/MAPS/template_filtering_effect',  '_',PREFIX,'_','40sps');
% filtcor = load(fname);
% filtlf = -filtcor(:,2)*spslf/40;   % lf template shift due to filterig effect, but the sign convention is opposite to detection here  
% filthf = -filtcor(:,1)*spshf/40;

PREFIX = strcat(fam,'.lf',num2str(lolf),'-',num2str(hilf),'.hf',num2str(lohf),'-',num2str(hihf),...
    '.hw',num2str(winlensechf),'.lw',num2str(winlenseclf),'.lo',num2str(loopoffmaxlf),'.ccm', ...
    num2str(xcmaxAVEnminlf),'.','80sps');
fname = strcat(rstpath, '/MAPS/tempfeff_PGC_enc','_',PREFIX);
filtcor = load(fname);
filtlf = filtcor(:,2)*spslf/80;   % lf template shift due to filterig effect, sign is the same
filthf = filtcor(:,1)*spshf/80;

%%% get the original offset
lftimeori = lftimeoricor;
if size(filtlf,1) > 2
    lftimeori(:,1) = lftimeori(:,1) + filtlf(2);
    lftimeori(:,2) = lftimeori(:,2) + filtlf(3);
else
    lftimeori(:,1) = lftimeori(:,1) + filtlf(1);
    lftimeori(:,2) = lftimeori(:,2) + filtlf(2);
end
lftimeori(:,3:4) = lftimeori(:,1:2)/spslf;

hftimeori = hftimeoricor;
if size(filthf,1) > 2
    hftimeori(:,1) = hftimeori(:,1) + filthf(2);
    hftimeori(:,2) = hftimeori(:,2) + filthf(3);
else
    hftimeori(:,1) = hftimeori(:,1) + filthf(1);
    hftimeori(:,2) = hftimeori(:,2) + filthf(2);
end
hftimeori(:,3:4) = hftimeori(:,1:2)/spshf;

[hfallori,hfallorimed,lfallori,dateori]=DetectOverlap(hftimeori,lftimeori,winlensechf,winlenseclf); 


% %% individual plots, zoom-in of suspicous migrations, call functions
% 
% %%% 1. migrations that are shown in both lf and hf
% % plot individual hf and lf
% % PlotIndividual(002,2003063,4e+4,4.2e+4,hfallorimed,lfallori);
% % PlotIndividual(002,2004197,4.3e+4,4.4e+4,hfallorimed,lfallori);
% % PlotIndividual(002,2004197,4.63e+4,4.65e+4,hfallorimed,lfallori);
% % PlotIndividual(002,2004197,8.4e+4,8.6e+4,hfallorimed,lfallori);
% % PlotIndividual(002,2004198,1.0e+4,1.4e+4,hfallorimed,lfallori);     % lf seems happen obviously after hf
% % PlotIndividual(002,2004198,8.35e+4,8.6e+4,hfallorimed,lfallori);
% % PlotIndividual(002,2004199,0,1e+4,hfallorimed,lfallori);
% % PlotIndividual(002,2005255,3.45e+4,3.75e+4,hfallorimed,lfallori);   % lf seems happen obviously after hf
% % PlotIndividual(002,2005255,5.85e+4,5.95e+4,hfallorimed,lfallori);   % lf seems happen obviously after hf
% % PlotIndividual(002,2005255,6.15e+4,6.3e+4,hfallorimed,lfallori);
% % PlotIndividual(002,2005255,6.75e+4,6.9e+4,hfallorimed,lfallori);    % lf seems happen obviously after hf
% % PlotIndividual(002,2005255,7.5e+4,7.6e+4,hfallorimed,lfallori);
% % PlotIndividual(002,2005256,0.25e+4,0.5e+4,hfallorimed,lfallori);
% % PlotIndividual(002,2005256,7.7e+4,8.55e+4,hfallorimed,lfallori);   % lf seems happen obviously after hf
% 
% % plot contemporal and average difference hf-lf,
% % Marked 'Y' means useful, seems like a plausible migration
% hflfdif=[];
% [ave12,ave13,hfmig,lfmig] = PlotIndividual2(002,2003063,4e+4,4.2e+4,hfallorimed,lfallori,hfallorimedcor,lfalloricor); 
% hflfdif=[hflfdif; ave12 ave13];
% [ave12,ave13,hfmig,lfmig] = PlotIndividual2(002,2004197,4.3e+4,4.4e+4,hfallorimed,lfallori,hfallorimedcor,lfalloricor);    % s3h,N-->S
% hflfdif=[hflfdif; ave12 ave13];
% [ave12,ave13,hfmig,lfmig] = PlotIndividual2(002,2004197,4.63e+4,4.65e+4,hfallorimed,lfallori,hfallorimedcor,lfalloricor);  % Y
% hflfdif=[hflfdif; ave12 ave13];
% [ave12,ave13,hfmig,lfmig] = PlotIndividual2(002,2004197,8.4e+4,8.6e+4,hfallorimed,lfallori,hfallorimedcor,lfalloricor);
% hflfdif=[hflfdif; ave12 ave13];
% [ave12,ave13,hfmig,lfmig] = PlotIndividual2(002,2004198,1.0e+4,1.4e+4,hfallorimed,lfallori,hfallorimedcor,lfalloricor);     % lf seems happen obviously after hf
% hflfdif=[hflfdif; ave12 ave13];
% [ave12,ave13,hfmig,lfmig] = PlotIndividual2(002,2004198,8.35e+4,8.6e+4,hfallorimed,lfallori,hfallorimedcor,lfalloricor);
% hflfdif=[hflfdif; ave12 ave13];
% [ave12,ave13,hfmig,lfmig] = PlotIndividual2(002,2004199,0,1e+4,hfallorimed,lfallori,hfallorimedcor,lfalloricor);   % Y
% hflfdif=[hflfdif; ave12 ave13];
% [ave12,ave13,hfmig,lfmig] = PlotIndividual2(002,2005255,3.45e+4,3.75e+4,hfallorimed,lfallori,hfallorimedcor,lfalloricor);   % Y, 6g,reversal,lf seems happen obviously after hf
% hflfdif=[hflfdif; ave12 ave13];
% [ave12,ave13,hfmig,lfmig] = PlotIndividual2(002,2005255,5.85e+4,5.95e+4,hfallorimed,lfallori,hfallorimedcor,lfalloricor);   % Y, 6h,reversal,lf seems happen obviously after hf
% hflfdif=[hflfdif; ave12 ave13];
% [ave12,ave13,hfmig,lfmig] = PlotIndividual2(002,2005255,6.15e+4,6.3e+4,hfallorimed,lfallori,hfallorimedcor,lfalloricor);    % Y
% hflfdif=[hflfdif; ave12 ave13];
% [ave12,ave13,hfmig,lfmig] = PlotIndividual2(002,2005255,6.75e+4,6.9e+4,hfallorimed,lfallori,hfallorimedcor,lfalloricor);    % Y, lf seems happen obviously after hf
% hflfdif=[hflfdif; ave12 ave13];
% [ave12,ave13,hfmig,lfmig] = PlotIndividual2(002,2005255,7.5e+4,7.6e+4,hfallorimed,lfallori,hfallorimedcor,lfalloricor);     % Y  
% hflfdif=[hflfdif; ave12 ave13];
% [ave12,ave13,hfmig,lfmig] = PlotIndividual2(002,2005256,0.25e+4,0.5e+4,hfallorimed,lfallori,hfallorimedcor,lfalloricor);    % Y
% hflfdif=[hflfdif; ave12 ave13];
% [ave12,ave13,hfmig,lfmig] = PlotIndividual2(002,2005256,7.7e+4,8.55e+4,hfallorimed,lfallori,hfallorimedcor,lfalloricor);   % Y,lf seems happen obviously after hf
% hflfdif=[hflfdif; ave12 ave13];
% 
% %%% According to the above, define the starting and ending time of those
% %%% plausible migrations
% trange = [2004199,0.20e+4,0.37e+4;
%           2005255,3.42e+4,3.65e+4;
%           2005255,5.80e+4,5.96e+4;
%           2005255,6.15e+4,6.26e+4;
%           2005255,6.70e+4,6.90e+4;
%           2005255,7.50e+4,7.60e+4;
%           2005256,0.36e+4,0.55e+4;
%           2005256,7.70e+4,8.55e+4;];
% seldif = zeros(length(trange),2);         
% for i=1: length(trange)
%     % get the median difference hf-lf & contemporal migration points
%     [seldif(i,1),seldif(i,2),conthfmig,contlfmig] = PlotIndividual2(fam,trange(i,1),trange(i,2),trange(i,3),hfallorimed,lfallori,hfallorimedcor,lfalloricor);    
%     
%     % save the contemporal hfmig and lfmig into files, same length
%     PREFIX = strcat(fam,'_',num2str(trange(i,1)),'_',num2str(trange(i,2)),'-',num2str(trange(i,3)),'.hfwin',num2str(winlensechf),'.lfwin',num2str(winlenseclf),'.',num2str(loopoffmaxlf),'.',num2str(xcmaxAVEnminlf));
%     fid = fopen(strcat(rstpath, '/MAPS/contemp_mig_hf','_',PREFIX),'w+');
%     fprintf(fid,'%d %d %.4f %.4f %d %.4f %.1f \n',conthfmig');
%     fclose(fid);
%     fid = fopen(strcat(rstpath, '/MAPS/contemp_mig_lf','_',PREFIX),'w+');
%     fprintf(fid,'%d %d %.4f %.4f %d %.4f %.1f \n',contlfmig');
%     fclose(fid);
%     
% %     % get the original corresponding migrations to check the general pattern & propagation direction, etc. 
% %     [orihfmig, orilfmig] = GetDetectionInMigration(fam,trange(i,1),trange(i,2),trange(i,3),hftimeoricor,lftimeoricor);
% %     fid = fopen(strcat(rstpath, '/MAPS/ori_mig_hf','_',PREFIX),'w+');
% %     fprintf(fid,'%d %d %.4f %.4f %d %.4f %.1f \n',orihfmig');
% %     fclose(fid);
% %     fid = fopen(strcat(rstpath, '/MAPS/ori_mig_lf','_',PREFIX),'w+');
% %     fprintf(fid,'%d %d %.4f %.4f %d %.4f %.1f \n',orilfmig');
% %     fclose(fid);
% end
%              
% 
% %%% 2. hf migrations shown in Rubin&Armbruster(2013)
% 
% [hfmig1, lfmig1] = GetDetectionInMigration(002,2003062,1.4e+4,3.6e+4,hftimeoricor,lftimeoricor);      % s2a, main front, forward
% [hfmig2, lfmig2] = GetDetectionInMigration(002,2003062,2.2e+4,2.3e+4,hftimeoricor,lftimeoricor);      % s2b, reversal
% [hfmig3, lfmig3] = GetDetectionInMigration(002,2003062,3.2e+4,3.3e+4,hftimeoricor,lftimeoricor);      % s2c, reversal
% [hfmig4, lfmig4] = GetDetectionInMigration(002,2003062,3.8e+4,3.9e+4,hftimeoricor,lftimeoricor);      % s2d, down-dip
% [hfmig5, lfmig5] = GetDetectionInMigration(002,2003062,7.5e+4,7.65e+4,hftimeoricor,lftimeoricor);     % s2e, forward
% [hfmig6, lfmig6] = GetDetectionInMigration(002,2003062,7.7e+4,7.8e+4,hftimeoricor,lftimeoricor);      % s2f, reversal
% [hfmig7, lfmig7] = GetDetectionInMigration(002,2003063,0.6e+4,0.75e+4,hftimeoricor,lftimeoricor);      % s2g, forward
% [hfmig8, lfmig8] = GetDetectionInMigration(002,2003063,6.5e+4,6.55e+4,hftimeoricor,lftimeoricor);      % s2h, N-->S
% [hfmig9, lfmig9] = GetDetectionInMigration(002,2004196,4.0e+4,4.008e+4,hftimeoricor,lftimeoricor);      % s3a, lighting up along the front
% [hfmig10, lfmig10] = GetDetectionInMigration(002,2004196,5.6e+4,5.66e+4,hftimeoricor,lftimeoricor);      % s3b, up-dip
% [hfmig11, lfmig11] = GetDetectionInMigration(002,2004196,7.0e+4,7.08e+4,hftimeoricor,lftimeoricor);      % s3c,reversal
% [hfmig12, lfmig12] = GetDetectionInMigration(002,2004196,7.38e+4,7.45e+4,hftimeoricor,lftimeoricor);      % s3d,reversal
% [hfmig13, lfmig13] = GetDetectionInMigration(002,2004197,3.06e+4,3.168e+4,hftimeoricor,lftimeoricor);      % s3e,reversal
% [hfmig14, lfmig14] = GetDetectionInMigration(002,2004197,3.74e+4,3.85e+4,hftimeoricor,lftimeoricor);      % s3f,reversal
% [hfmig15, lfmig15] = GetDetectionInMigration(002,2004197,3.91e+4,3.96e+4,hftimeoricor,lftimeoricor);      % s3g,forward
% [hfmig16, lfmig16] = GetDetectionInMigration(002,2004197,4.3e+4,4.37e+4,hftimeoricor,lftimeoricor);      % s3h,N-->S
% [hfmig17, lfmig17] = GetDetectionInMigration(002,2005254,2.7e+4,2.8e+4,hftimeoricor,lftimeoricor);      % 6a,reversal
% [hfmig18, lfmig18] = GetDetectionInMigration(002,2005254,3.0e+4,3.2e+4,hftimeoricor,lftimeoricor);      % 6b,reversal
% [hfmig19, lfmig19] = GetDetectionInMigration(002,2005254,3.4e+4,3.5e+4,hftimeoricor,lftimeoricor);      % 6c,N-->S
% [hfmig20, lfmig20] = GetDetectionInMigration(002,2005254,3.5e+4,3.6e+4,hftimeoricor,lftimeoricor);      % 6d,down-dip
% [hfmig21, lfmig21] = GetDetectionInMigration(002,2005254,7.3e+4,7.4e+4,hftimeoricor,lftimeoricor);      % 6e,reversal
% [hfmig22, lfmig22] = GetDetectionInMigration(002,2005255,0.47e+4,0.6e+4,hftimeoricor,lftimeoricor);      % 6f,reversal
% [hfmig23, lfmig23] = GetDetectionInMigration(002,2005255,3.46e+4,3.6e+4,hftimeoricor,lftimeoricor);      % 6g,reversal
% [hfmig24, lfmig24] = GetDetectionInMigration(002,2005255,5.8e+4,6.0e+4,hftimeoricor,lftimeoricor);      % 6h,reversal
% 
% 
% %%% 3. hf and lf migrations worthwhile to check the projection in propagation direction, maybe some
% %%%     of them don't have many contemporaneous detections
% trange3 = [2004197,8.45e+4,8.55e+4;
%            2004198,8.35e+4,8.62e+4;
%            2004199,0.20e+4,0.37e+4;
%            2005255,3.42e+4,3.65e+4;
%            2005255,5.80e+4,5.96e+4;
%            2005255,6.15e+4,6.26e+4;
%            2005255,6.70e+4,6.90e+4;
%            2005255,7.50e+4,7.60e+4;
%            2005256,0.36e+4,0.55e+4;
%            2005256,7.62e+4,8.30e+4;
%            2005256,8.39e+4,8.46e+4;
%            2005256,8.47e+4,8.55e+4;];
% 
% %%% 4. an example, clear hf but almost no lf
% trange4 = [2003062,7.45e+4,7.68e+4;
%            2005255,0.45e+4,0.62e+4;]; 
% 
% %%% 5. an example, clean lf, but messy hf
% trange5 = [2005255,7.2e+4,7.4e+4];

%%
%%% 6. fine selected migrations that undergo linear fit analysis
trange6 = [2005255,3.42e+4,3.65e+4;
          2005255,6.15e+4,6.26e+4;
          2005255,6.70e+4,6.90e+4;
          2005255,7.50e+4,7.60e+4;
          2005256,0.36e+4,0.55e+4;
          2005256,7.62e+4,8.30e+4;];     %2004197,8.45e+4,8.55e+4;
rstpath = strcat(getenv('MHOME'), '/Seisbasics/hypoinverse/forsummary');      

for i = 1: length(trange6)
% i=6;
    date = trange6(i,1);
    sttime = trange6(i,2)/3600;
    edtime = trange6(i,3)/3600;
    datestr = num2str(date);
    msize = 30;
    yr = datestr(1:4);
    jday = datestr(5:end);
    hfday = hftimeori(hftimeori(:,5)==date, :);
    lfday = lftimeori(lftimeori(:,5)==date, :);
    hfday(:,end) = hfday(:,end)/3600;
    lfday(:,end) = lfday(:,end)/3600;
    indhf = find(hfday(:,end)>=sttime & hfday(:,end)<=edtime);
    hfmig = hfday(indhf,:);
    indlf = find(lfday(:,end)>=sttime & lfday(:,end)<=edtime);
    lfmig = lfday(indlf,:);
    spsratio=2;
    
    figure('Position',[scrsz(3)*0.2 scrsz(4)*0.2 scrsz(3)*0.5 scrsz(4)*0.35]);
    
    subplot(1,2,1,'align');
    hold on
    grid on
    box on
    p1=scatter(hfmig(:,end),hfmig(:,1),msize,'filled','ro','MarkerEdgeColor', 'k');   % 12 off hf
    p2=scatter(lfmig(:,end),lfmig(:,1)*spsratio,msize,[153/255 255/255 255/255],'filled','o',...
                'MarkerEdgeColor', 'k');   % 12 off lf    
    % create fit object
    fttpfree = fittype( @(a,b,x) a*x+b);
    [fitobjhf12,gofhf12,~] = fit(hfmig(:,end), hfmig(:,1),fttpfree,'Robust','Bisquare');
    % output fit parameters
    coef12 = coeffvalues(fitobjhf12);
    slopehf12 = coef12(1);
    intcpthf12 = coef12(2);
    fithf12 = feval(fitobjhf12,hfmig(:,end));
    plot(hfmig(:,end),fithf12,'r-','linewidth',2);
%     timeseqhf12 = min(hfmig(:,end)): 0.05: max(hfmig(:,end));
%     fithf12cor = slopehf12*timeseqhf12+intcpthf12;
%     plot(timeseqhf12,fithf12cor,'--','color',[153/255 0/255 0/255],'linewidth',2);
    
    a=slopehf12;
    fttpfix = fittype( @(b,x) a*x+b);
    [fitobjlf12,goflf12,~] = fit(lfmig(:,end), lfmig(:,1)*spsratio,fttpfix,'Robust','Bisquare');
    intcptlf12 = coeffvalues(fitobjlf12);
    fitlf12 = feval(fitobjlf12,lfmig(:,end));
    plot(lfmig(:,end),fitlf12,'-','linewidth',2,'color',[153/255 255/255 255/255]);
    timeseqlf12 = min(lfmig(:,end)): 0.05: max(lfmig(:,end));
    fitlf12cor = slopehf12*timeseqlf12+intcptlf12-(filtlf(1)*spsratio-filthf(1));
    plot(timeseqlf12,fitlf12cor,'--','linewidth',2,'color',[0/255 255/255 255/255]);

    legend([p1,p2],{'Ori. HF off12','Ori. LF off12'},'location','Northeast','Fontsize',10);
    text(0.04,0.93,'a','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
    set(gca, 'xlim',[sttime edtime]);
    set(gca, 'ylim',[-30 30]);
    xlabel(strcat({'Time (hr) of '},yr,jday), 'fontsize',12);
    ylabel('Samples', 'fontsize',12);
%     title(strcat(fam,'.',yr,'.',jday,'\_',num2str(sttime),'-',num2str(edtime),'_hr'), 'fontsize', 14);

    
    subplot(1,2,2,'align');
    hold on
    grid on
    box on
    p1=scatter(hfmig(:,end),hfmig(:,2),msize,'filled','rs','MarkerEdgeColor', 'k');   % 13 off hf
    p2=scatter(lfmig(:,end),lfmig(:,2)*spsratio,msize,[153/255 255/255 255/255],'filled','s',...
                'MarkerEdgeColor', 'k');   % 13 off lf    
    % create fit object
    fttpfree = fittype( @(a,b,x) a*x+b);
%     [fitobjhf13,gofhf13,~] = fit(hfmig(:,end), hfmig(:,2),fttpfree);
    [fitobjhf13,gofhf13,~] = fit(hfmig(:,end), hfmig(:,2),fttpfree,'Robust','Bisquare');
    % output fit parameters
    coef13 = coeffvalues(fitobjhf13);
    slopehf13 = coef13(1);
    intcpthf13 = coef13(2);
    fithf13 = feval(fitobjhf13,hfmig(:,end));
    plot(hfmig(:,end),fithf13,'r-','linewidth',2);
%     timeseqhf13 = min(hfmig(:,end)): 0.05: max(hfmig(:,end));
%     fithf13cor = slopehf13*timeseqhf13+intcpthf13;
%     plot(timeseqhf13,fithf13cor,'--','color',[153/255 0/255 0/255],'linewidth',2);
    
    a=slopehf13;
    fttpfix = fittype( @(b,x) a*x+b);
%     [fitobjlf13,goflf13,~] = fit(lfmig(:,end), lfmig(:,2),fttpfix);
    [fitobjlf13,goflf13,outplf13] = fit(lfmig(:,end), lfmig(:,2)*spsratio,fttpfix,'Robust','Bisquare');
    intcptlf13 = coeffvalues(fitobjlf13);
    fitlf13 = feval(fitobjlf13,lfmig(:,end));
    timeseqlf13 = min(lfmig(:,end)): 0.05: max(lfmig(:,end));
    fitlf13cor = slopehf13*timeseqlf13+intcptlf13-(filtlf(2)*spsratio-filthf(2));
    plot(lfmig(:,end),fitlf13,'-','linewidth',2,'color',[153/255 255/255 255/255]);
    plot(timeseqlf13,fitlf13cor,'--','linewidth',2,'color',[0/255 255/255 255/255]);
    
    legend([p1,p2],{'Ori. HF off13','Ori. LF off13'},'location','Northeast','Fontsize',10);
    text(0.04,0.93,'b','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
    set(gca, 'xlim',[sttime edtime]);
    set(gca, 'ylim',[-30 30]);
    xlabel(strcat({'Time (hr) of '},yr,jday), 'fontsize',12);
    ylabel('Samples', 'fontsize',12);
    
    
    %%% save figure
    print('-depsc2',strcat(rstpath, '/',fam,'offset_statistics.contemp_migration.',num2str(trange6(i,1)),...
    '_',num2str(trange6(i,2)),'-',num2str(trange6(i,3)),'.',num2str(winlensechf),'_',...
    num2str(winlenseclf),'.',num2str(loopoffmaxlf),'.',num2str(xcmaxAVEnminlf),'.eps'));
end
      
      
      
% keyboard


























