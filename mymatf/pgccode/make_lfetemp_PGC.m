%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to make the lfe templates for PGC trio for fam 002
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/11/13
% Last modified date:   2019/11/13

%% Initialization
format short e   % Set the format to 5-digit floating point
% clear
% close all

% set(0,'DefaultFigureVisible','on');
set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
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

ifam = 1;
fam = nfampool(ifam,:);
stas=['PGC  '
      'SSIB '
      'SILB '];
nsta=size(stas,1);         %  number of stations
templensec = 60;

sps = 80;   % sampling rate you want
ccmethod = 2;  % only using lfes passed the cc thresholds (ccmethod=1); using all lfes (ccmethod=2)
plflag = 1;  % whether to plot (1/0)

%% make broadband stacks,
%%% only using lfes passed the cc thresholds (ccmethod=1)
%%% using all lfes (ccmethod=2)
CATA = 'fixed';
ccbp = [2 8];
[dstack, ccstack] = mk_bbtemp_PGC(fam,CATA,sps,templensec,ccmethod,ccbp,plflag);
% write into files
for ista = 1: nsta
    fid = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
        num2str(sps), 'sps_', num2str(templensec), 's_', ...
        'BBDS_', 'opt_Nof_Non_Chao'), 'w+');  % broadband direct stack, no filter, no norm
        
    fprintf(fid, '%f \n', dstack(ista, :)');
    fclose(fid);

    fid = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
        num2str(sps), 'sps_', num2str(templensec), 's_', ...
        'BBCCS_', 'opt_Nof_Non_Chao'), 'w+');  % broadband cc stack, no filter, no norm
    fprintf(fid, '%f \n', ccstack(ista, :)');
    fclose(fid);
end


%% make bandpassed stacks,
%%% only using lfes passed the cc thresholds (ccmethod=1)
%%% using all lfes (ccmethod=2)
bplo = 0.1;
bphi = 15;
[dstack, ccstack] = mk_bptemp_PGC(fam,sps,bplo,bphi,ccmethod,plflag);
% write into files
for ista = 1: nsta
    fid = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
        num2str(sps), 'sps_', num2str(templensec), 's_', ...
        'BPDS_', 'opt_Nof_Non_Chao'), 'w+');  % bandpassed direct stack, no filter, no norm
        
    fprintf(fid, '%f \n', dstack(ista, :)');
    fclose(fid);
    
    fid = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
        num2str(sps), 'sps_', num2str(templensec), 's_', ...
        'BPCCS_', 'opt_Nof_Non_Chao'), 'w+');  % bandpassed cc stack, no filter, no norm
    fprintf(fid, '%f \n', ccstack(ista, :)');
    fclose(fid);
end




















