%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is to analyse the relationship between the hf and low
% frequency detection, find the potential overlapping/union, etc.
%
% Read 'mapaddcheck_xxx_combined'
%
% Feature:
%   1. read lf & hf result file, already corrected for filtering effect
%   2. since the hf result is denser, so extend the lf area roughly 2
%   times, based on each point, 1-->9, then average the overlapping edge
%   3. subtract to get the intersection part (difference if they have a
%   common activatited point
%   4. plot and writing files
%
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/06/15
% Last modified date:   2019/06/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio

format short e   % Set the format to 5-digit floating point
clear
% close all

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
winlenseclf=24;       % CHECK in detection script for what were used in first-year report
winoffseclf=12;        % window offset in sec, which is the step of a moving window
winlenlf=winlenseclf*spslf;      % length in smaples
hilf=1.25;
lolf=0.5;
npo=2;     % poles, passes of filters
npa=2;
cyclskiplf = 20;
mshiftlf=19+cyclskiplf; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
loopoffmaxlf=2; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
xcmaxAVEnminlf=0.45; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?


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

%%% plot
figure('Position',[scrsz(3)/5 scrsz(4)/5 2*scrsz(3)/5 4*scrsz(4)/5]);
ax=NaN(nstanew,1);

nrow = nstanew;
ncol = 1;
msize = 15;

%% get offset difference between hf and lf results
for istanew = 1: nstanew
%     istanew =1;

    %%% READ FILES
    % read hf detection file
    fname= strcat(rstpath, '/MAPS/mapaddcheck_', stasnew(istanew,:), '_combined_',PREFIXhf,'_',num2str(lohf),'-',num2str(hihf),'_',num2str(winlenhf/spshf),'s',num2str(spshf),'sps', num2str(nstanew), 'newsta');
    hfrst = load(fname);
    hfrst(3:6) = vpa(hfrst(3:6),4);
    % savefile14 = [STA124off STA134off STA124offsec STA134offsec STA14meansec STA14stdsec STA14count];
    STA124off_hf = hfrst(:,1);
    STA134off_hf = hfrst(:,2);
    STA124offsec_hf = hfrst(:,3);
    STA134offsec_hf = hfrst(:,4);
    STA14meansec_hf = hfrst(:,5);
    STA14stdsec_hf = hfrst(:,6);

    % read lf detection file
    fname= strcat(rstpath, '/MAPS/mapaddcheck_', stasnew(istanew,:), '_',PREFIXlf,'_',num2str(lolf),'-',num2str(hilf),'_',num2str(winlenlf/spslf),'s',num2str(spslf),'sps', num2str(nstanew), 'newsta');
    lfrst = load(fname);
    lfrst(3:6) = vpa(lfrst(3:6),4);
    % savefile14 = [STA124off STA134off STA124offsec STA134offsec STA14meansec STA14stdsec STA14count];
    STA124off_lf = lfrst(:,1);
    STA134off_lf = lfrst(:,2);
    STA124offsec_lf = lfrst(:,3);
    STA134offsec_lf = lfrst(:,4);
    STA14meansec_lf = lfrst(:,5);
    STA14stdsec_lf = lfrst(:,6);
    
    %%% construct the extended lf result
    extrstraw = [];
    for i = 1: size(lfrst,1)
% i=1
        rangex = lfrst(i,3)-1/spslf/2: 1/spslf/2 :lfrst(i,3)+1/spslf/2;
        rangey = lfrst(i,4)-1/spslf/2: 1/spslf/2 :lfrst(i,4)+1/spslf/2;
        [grdx, grdy] = meshgrid(rangex, rangey);
        duppts = [reshape(grdx,[],1) reshape(grdy,[],1)];   % duplicate points
        % use reptmat to repeat rows
        duprst = [floor(duppts*spshf) duppts repmat(lfrst(i,5:end),[size(duppts,1) 1])]; % duplicate results
        extrstraw = [extrstraw; duprst];
    end
    
    %%% the extended result has repeated locations, unique it according to
    %%% std, and use average of mean offset
    extrstsort = sortrows(extrstraw, [1,2]);
    [ptsuni, isort, ~] = unique(extrstsort(:,1:2), 'rows','stable');
    extrst=[];  % extrst is the extended lf result
    for i = 1: length(isort)
%         ind = find(extrstsort(isort,1:4) == extrstsort(min(isort+1,size(extrstsort,1)), 1:4));
        if extrstsort(isort(i),1:2) == extrstsort(min(isort(i)+1,size(extrstsort,1)), 1:2)
            if extrstsort(isort(i),6) < extrstsort(min(isort(i)+1,size(extrstsort,1)),6)  % preserve the one with smaller std
                extrst = [extrst; extrstsort(isort(i), :)];
            elseif extrstsort(isort(i),6) == extrstsort(min(isort(i)+1,size(extrstsort,1)),6)
                temp = (extrstsort(isort(i),:) + extrstsort(min(isort(i)+1,size(extrstsort,1)),:))/2;
                extrst = [extrst; temp];
            else
                extrst = [extrst; extrstsort(min(isort(i)+1,size(extrstsort,1)),:)];
            end
        else
            extrst = [extrst; extrstsort(isort(i), :)];
        end
        
    end
    
    %%% get the intersect part between extended lf and hf results
    [~, ilf, ihf] = intersect(extrst(:,1:2), hfrst(:,1:2), 'row', 'stable');
    intsectlf = extrst(ilf,:);
    intsecthf = hfrst(ihf,:);
    
    %%% in filtering correction, the 3rd col is hf-lf, be consistent here
    diffoff = [intsectlf(:, 1:4) intsecthf(:,5)-intsectlf(:,5)];    
    

    %%% plot
    % additional sta 1, mean offset
    ax(istanew)=subplot(nrow, ncol, istanew,'align');
    scatter(diffoff(:,3), diffoff(:,4), msize, diffoff(:,5), 'filled','s');    %, 'MarkerEdgeColor', 'w')
    colormap(ax(istanew),jet)
    colorbar
    caxis([-0.5, 0.5]);
    box on
    axis equal
    xlabel('STA12 offset (s)')
    ylabel('STA13 offset (s)')
    axis([-1 1.5 -1.5 1]);
    title(strcat('meanoff\_difference\_in\_hf-lf\_', strtrim(stasnew(istanew,:))));

end

PREFIX = strcat(fam,'_with_extendlf');
print('-depsc',strcat(rstpath,'/FIGS/',PREFIX,'_','SAT14_meanoff_diff_lf_hf','_',num2str(nstanew), 'newsta.eps'));


























