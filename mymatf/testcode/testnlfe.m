function [timoffrot,nLFE,lfe] = testnlfe(fam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  [timoffrot,nLFE,lfe] = testnlfe(fam)
%
%   small function to see the number of lfe detections in 
%   Bostock's new catalog of fam queried
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/12
% Last modified date:   2020/02/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

workpath = getenv('ALLAN');
temppath = strcat(workpath, '/templates/');
if ~exist(temppath, 'dir')
    mkdir(temppath);
end
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/PGCtrio');

% old LFE family pool, from Yajun's selections
ofampool = ['002';  % 002,246, choose 002
            '013';  % 043,152, choose 043
            '025';  % 141,
            '028';  % 047
            '056';  % 010
            '084';  % 144,102,121,058,257, choose 144
            '099';  % 099,006, choose 099
            '115';  % 068
            '125';  % 125
            '147';  % 147,017, choose 147
            '149']; % 017,047, choose 017 
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
fam
timoffrot = Readbostock(fam);
% get the all catalog LFEs in that family
bostname = ('/BOSTOCK/total_mag_detect_0000_cull_NEW.txt');
%%% 6 cols: 
%%%     family ID, date (ymmdd), hour, sec, magnitude, number of stations that include that event
%%%     in their template or magnitude determination or something
catalog = load(strcat(workpath, bostname));
famnum = str2double(fam);
dateall = catalog(famnum == catalog(:, 1), :);
nday = size(timoffrot, 1)
nLFE = 0;
lfe = [];
for id = 1: nday
    %     nd = 1;
%     close all
    year = timoffrot(id,1);
    YEAR = int2str(year);
    jday = timoffrot(id,2);
    if year == 2003
        date = jday-62+30303;
    elseif year == 2004
        date = jday-196+40714;
    elseif year == 2005
        date = jday-254+50911;
    end
    bostocks = dateall(date == dateall(:, 2), :);
    bostsec = 3600*(bostocks(:,3)-1)+bostocks(:,4);
    bostsamp = round(bostsec*40);
    nlfe1day = size(bostocks, 1);
    nLFE = nLFE + size(bostsamp, 1);
    
    lfe = [lfe; bostocks];
end
nLFE
    