%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This is to driving test script to test the pgc 002 templates from using
%   the different LFE catalogs, as the OLD and NEW catalog have different days
%   and load in different rotatation parameters obtained from the different 
%   different catalog
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2021/07/31
% Last modified date:   2021/07/31
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
format short e   % Set the format to 5-digit floating point
clear
clc
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);
% scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]
% wid=scrsz(3);       % width
% hite=scrsz(4);       % height

workpath = getenv('ALLAN');
%%% choose the data path here
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');

% CATA = 'old';   % use the rot param from the OLD catalog (dates and so on)
CATA = 'new';   % use the rot param from the NEW catalog (dates and so on)

disp(CATA);

temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');

nfampool = ['002';
            ];
nfam = size(nfampool,1);

stas=['PGC  '
    'SSIB '
    'SILB '];

% number of used stations
nsta=size(stas,1);         %  number of stations


%% make template or read template
remake = 1;  % re-make the template or not, 1/0
dstack = [];
ccstack = [];
sps = 100;
templensec = 60;
% templensec = 32*4;
for ifam = 1: nfam
%     ifam=2;
    fam = nfampool(ifam, :);
    disp(fam);
    
    if remake   % if requested templates do not exist, recompute them
        ccmethod = 1; % only using lfes passed the cc thresholds (ccmethod=1); using all lfes (ccmethod=2)  
        plflag = 1;
        [dstack, ccstack] = testmk_bbtemp_PGC(fam,CATA,sps,templensec,ccmethod,plflag);
        
        %write into files
        for ista = 1: nsta
            fid = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
                num2str(sps), 'sps_', num2str(templensec), 's_', ...
                'BBDS_', 'opt_Nof_Non_Chao_cat',CATA,'_m1cc0.5-6.5'), 'w+');  % bandpassed direct stack, no filter, no norm
            
            fprintf(fid, '%f \n', dstack(ista, :)');
            fclose(fid);
            
            fid = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
                num2str(sps), 'sps_', num2str(templensec), 's_', ...
                'BBCCS_', 'opt_Nof_Non_Chao_cat',CATA,'_m1cc0.5-6.5'), 'w+');  % bandpassed cc stack, no filter, no norm
            fprintf(fid, '%f \n', ccstack(ista, :)');
            fclose(fid);
        end

        keyboard
                
    else    % if exist, load them directly
        for ista = 1: nsta
            fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
                num2str(templensec), 's_', 'BBDS_', 'opt_Nof_Non_Chao');
            dstack(ista,:) = load(fname);
            fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
                num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao');
            ccstack(ista,:) = load(fname);
        end
        stack = ccstack;
    end
    
end