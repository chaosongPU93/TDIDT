%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is especially to interpolate the hf detection results after
% adding additional station check, also after avoiding reduntant counting,
% in order to give a reference for lf detection
% 
% Either at the whole rectangle, or only the lf activated points
%
% Feature:
%   1. read result file
%   2. interpolation, using points according to reality
%   3. plot and writing files
%   4. use files 'maporioffcount','mapaddcheck_xxxx_combined', excluding
%   time
%
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/06/13
% Last modified date:   2019/06/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio

format short e   % Set the format to 5-digit floating point
clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

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
spshf=40;     % samples per second

%%% IMPORTANT, NEED to change sometimes
%Basics of the cross-correlation:  Window length, number of windows, filter parameters, etc.
%%%% HF %%%%
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

% stasnew=['LZB  '
%          'TWKB '
%          'MGCB '];
stasnew=['LZB  '
         'TWKB '
         'MGCB '
         'KLNB '];     
nstanew=size(stasnew,1);
mshiftnew=ceil(2.2*mshifthf);
scrsz=get(0,'ScreenSize');

PREFIXhf = strcat(fam,'.loff',num2str(loopoffmaxhf),'.ccmin',num2str(xcmaxAVEnminhf),'.nponpa',int2str(npo),int2str(npa),'.ms', int2str(mshifthf));
PREFIXlf = strcat(fam,'.loff',num2str(loopoffmaxlf),'.ccmin',num2str(xcmaxAVEnminlf),'.nponpa',int2str(npo),int2str(npa),'.ms', int2str(mshiftlf));

%%% interpolation to denser region, supscript 4 means the additional
%%% station
STA124offrange = -40:1:40;
STA134offrange = -40:1:40;

% meshgrid and then reshape is just a way to get every point pair of
% interpolation region
[regionx,regiony] = meshgrid(STA124offrange,STA134offrange);
STA124offitp = reshape(regionx', [], 1);
STA134offitp = reshape(regiony', [], 1);
npts = length(STA124offitp);
STA14meansecitp = zeros(npts, nstanew);

% read lf detection file
fname=strcat(rstpath, '/MAPS/maporioffcount_',PREFIXlf,'_',num2str(lolf),'-',num2str(hilf),'_',num2str(winlenlf/spslf),'s',num2str(spslf),'sps', num2str(nstanew), 'newsta');
savefile23 = load(fname);
%savefile23 = [STA12off STA13off STA12offsec STA13offsec STA23count];
STA12offsec_lf = savefile23(:,3);
STA13offsec_lf = savefile23(:,4);
STA23meansecitp = zeros(length(STA12offsec_lf), nstanew);


figure('Position',[scrsz(3)/10 scrsz(4)/10 3*scrsz(3)/5 4*scrsz(4)/5]);

for istanew = 1: nstanew
    
    % read hf detection file
    fname= strcat(rstpath, '/MAPS/mapaddcheck_', stasnew(istanew,:), '_combined_',PREFIXhf,'_',num2str(lohf),'-',num2str(hihf),'_',num2str(winlenhf/spshf),'s',num2str(spshf),'sps', num2str(nstanew), 'newsta');
    savefile14 = load(fname);
    % savefile14 = [STA124off STA134off STA124offsec STA134offsec STA14meansec STA14stdsec STA14count];
    STA124off_hf = savefile14(:,1);
    STA134off_hf = savefile14(:,2);
    STA124offsec_hf = savefile14(:,3);
    STA134offsec_hf = savefile14(:,4);
    STA14meansec_hf = savefile14(:,5);
    STA14stdsec_hf = savefile14(:,6);
    
    % Here the results ready meet the standard deviation criteria,
    % to prevent potential outliers to contaminate the results, use the
    % stameansec to that has mild value, set a tolerance as well
    tolmean = 10;
    STA14meansec_hf = STA14meansec_hf(abs(STA14meansec_hf) <= tolmean);
    STA124off_hf = STA124off_hf(abs(STA14meansec_hf) <= tolmean);
    STA134off_hf = STA134off_hf(abs(STA14meansec_hf) <= tolmean);
          
    % interpolation with 'scatteredInterpolant' class at a rectangle region
    F = scatteredInterpolant(STA124off_hf,STA134off_hf,STA14meansec_hf,'natural','linear');
    STA14meansecitp(:,istanew) = F(STA124offitp,STA134offitp);
    
    % interpolation only at lf detection points
    F = scatteredInterpolant(STA124offsec_hf,STA134offsec_hf,STA14meansec_hf,'natural','linear');
    STA23meansecitp(:,istanew) = F(STA12offsec_lf, STA13offsec_lf);
    
    %%% plot to check the interpolation result
    msize=15;
    
    ax1=subplot(nstanew, 3, 3*istanew-2,'align');
    tempind = find(abs(STA14meansec_hf)>0.1);
    scatter(STA124off_hf(tempind)/spshf, STA134off_hf(tempind)/spshf, msize, STA14meansec_hf(tempind), 'filled','s');    %, 'MarkerEdgeColor', 'w')
    hold on
    tempind = find(abs(STA14meansec_hf)<=0.1);
    scatter(STA124off_hf(tempind)/spshf, STA134off_hf(tempind)/spshf, msize, STA14meansec_hf(tempind), 'filled','s','MarkerEdgeColor', 'k');    colormap(ax1,jet)
    colorbar
    caxis([-1.5, 1.5]);
    box on
    grid on
    axis equal
    xlabel('STA12 offset (s)')
    ylabel('STA13 offset (s)')
    axis([-1 1.5 -1.5 1]);
    title(strcat(num2str(lohf),'-',num2str(hihf),'\_hf\_ori\_', strtrim(stasnew(istanew, :))));

    ax2=subplot(nstanew, 3, 3*istanew-1,'align');
    scatter(STA124offitp/spshf, STA134offitp/spshf, msize, STA14meansecitp(:,istanew), 'filled','s');    %, 'MarkerEdgeColor', 'w')
    colormap(ax2,jet)
    colorbar
    caxis([-1.5, 1.5]);
    box on
    grid on
    axis equal
    xlabel('STA12 offset (s)')
    ylabel('STA13 offset (s)')
    axis([-1 1.5 -1.5 1]);
    title(strcat(num2str(lohf),'-',num2str(hihf),'\_hf\_itp\_rect\_', strtrim(stasnew(istanew, :))));
    
    ax3=subplot(nstanew, 3, 3*istanew,'align');
    tempind = find(abs(STA23meansecitp(:,istanew))>0.1);
    scatter(STA12offsec_lf(tempind), STA13offsec_lf(tempind), msize, STA23meansecitp(tempind,istanew), 'filled','s');    %, 'MarkerEdgeColor', 'w')
    hold on
    tempind = find(abs(STA23meansecitp(:,istanew))<=0.1);
    scatter(STA12offsec_lf(tempind), STA13offsec_lf(tempind), msize, STA23meansecitp(tempind,istanew), 'filled','s','MarkerEdgeColor', 'k');
    colormap(ax3,jet)
    colorbar
    caxis([-1.5, 1.5]);
    box on
    grid on
    axis equal
    xlabel('STA12 offset (s)')
    ylabel('STA13 offset (s)')
    axis([-1 1.5 -1.5 1]);
    title(strcat(num2str(lohf),'-',num2str(hihf),'\_hf\_itp\_lf\_', strtrim(stasnew(istanew, :))));
    
    
end

% construct the interpolation matrix for the reference in the following
STA14ref = [STA124offitp STA134offitp STA124offitp/spshf STA134offitp/spshf STA14meansecitp];