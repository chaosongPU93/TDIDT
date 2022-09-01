%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to call the function 'calc_rots_noresp' to compute the 
% rotation parameters for all fams.
%       
%      
%
% NOTES:
%   1. there are fams that the stacked templates are pretty low-amplitude,
%   then the maximum offset allowed in CC would be necessary, if the result
%   reaches the boundary, then it is necessary to decrease the threshold.
%   However, check the N and E component to see if the splitting should be
%   small or large. For example, fam 010 should allow splitthres <=8, but
%   others could works fine when splitthres <=12, and usually this offset
%   ranges from 0~9.
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/11/10
% Last modified date:   2019/11/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]
wid=scrsz(3);       % width
hite=scrsz(4);       % height

% set path
workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath,'/split_chao/lzbtrio/');

% old LFE family pool, from Yajun's selections
ofampool = ['002';  % 002,246
            '013';  % 043,152
            '025';  % 141
            '028';  % 047
            '056';  % 010
            '084';  % 144,102,121,58,257
            '099';  % 099,006
            '115';  % 068
            '125';  % 125
            '147';  % 147,17
            '149']; % 017,047

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
nfam = length(nfampool);
        
%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
    'LZB'];
POLSTA=['SSIB '           % polaris station names
    'SILB '
    'KLNB '
    'MGCB '
    'TWKB '];

stas=['PGC  '
    'SSIB '
    'SILB '
    'LZB  '
    'TWKB '
    'MGCB '
    'KLNB '];
% number of used stations
nsta=size(stas,1);         %  number of stations


%% calculation
ROTSALL = zeros(nfam*nsta, 4);

remake = 1;  % whether to calculate the rots again. 0/1

lo = 0.5;
hi = 6.5;
sps = 40;
bef = 25;
aft = 35;

% [ROTS]=calc_rots_noresp(fam,lo,hi,sps,splitoff,staoff,mainsta)

if remake   % re-calculate all results
%     %%% fam 002
%     fam = '002';
%     splitoff = 12;
%     staoff = 8;
%     mainsta = 4; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
%     
%     
%     %%% fam 043
%     fam = '043';
%     splitoff = 12;
%     staoff = 8;
%     mainsta = 4; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((2-1)*nsta+1:2*nsta, :) = ROTS;
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
% 
%     
%     %%% fam 141
%     fam = '141';
%     splitoff = 12;
%     staoff = 8;
%     mainsta = 4; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((3-1)*nsta+1:3*nsta, :) = ROTS;
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
%     
%     
%     %%% fam 047
%     fam = '047';
%     splitoff = 12;
%     staoff = 8;
%     mainsta = 4; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((4-1)*nsta+1:4*nsta, :) = ROTS;
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
%     
%     
%     %%% fam 010
%     fam = '010';
%     splitoff = 8;
%     staoff = 8;
%     mainsta = 4; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((5-1)*nsta+1:5*nsta, :) = ROTS;
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
% 
%     
%     %%% fam 144
%     fam = '144';
%     splitoff = 12;
%     staoff = 8;
%     mainsta = 5; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((6-1)*nsta+1:6*nsta, :) = ROTS;
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
%     
%     
%     %%% fam 099
%     fam = '099';
%     splitoff = 12;
%     staoff = 8;
%     mainsta = 4; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((7-1)*nsta+1:7*nsta, :) = ROTS;
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
    
    
    %%% fam 068
    fam = '068';
    splitoff = 12;
    staoff = 8;
    mainsta = 5; 
    [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
    % save all in one file
    ROTSALL((8-1)*nsta+1:8*nsta, :) = ROTS;
    % save results for each family
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
    fprintf(fid,'%d %d %d %d \n',ROTS');
    fclose(fid);
    
    
    %%% fam 125
    fam = '125';
    splitoff = 12;
    staoff = 8;
    mainsta = 5; 
    [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
    % save all in one file
    ROTSALL((9-1)*nsta+1:9*nsta, :) = ROTS;
    % save results for each family
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
    fprintf(fid,'%d %d %d %d \n',ROTS');
    fclose(fid);
    
    
    %%% fam 147
    fam = '147';
    splitoff = 12;
    staoff = 8;
    mainsta = 5; 
    [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
    % save all in one file
    ROTSALL((10-1)*nsta+1:10*nsta, :) = ROTS;
    % save results for each family
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
    fprintf(fid,'%d %d %d %d \n',ROTS');
    fclose(fid);
    
    
    %%% fam 017
    fam = '017';
    splitoff = 12;
    staoff = 8;
    mainsta = 4; 
    [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
    % save all in one file
    ROTSALL((11-1)*nsta+1:11*nsta, :) = ROTS;
    % save results for each family
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
    fprintf(fid,'%d %d %d %d \n',ROTS');
    fclose(fid);
    
    
else    % directly read results
    
    %%% fam 002
    fam = '002';
    mainsta = 4; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    % save all in one file
    ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
    
    %%% fam 043
    fam = '043';
    mainsta = 4; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    % save all in one file
    ROTSALL((2-1)*nsta+1:2*nsta, :) = ROTS;
    
    %%% fam 141
    fam = '141';
    mainsta = 4; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    % save all in one file
    ROTSALL((3-1)*nsta+1:3*nsta, :) = ROTS;
    
    %%% fam 047
    fam = '047';
    mainsta = 4; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    % save all in one file
    ROTSALL((4-1)*nsta+1:4*nsta, :) = ROTS;
    
    %%% fam 010
    fam = '010';
    mainsta = 4; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    % save all in one file
    ROTSALL((5-1)*nsta+1:5*nsta, :) = ROTS;
    
    %%% fam 144
    fam = '144';
    mainsta = 5; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    ROTS(:,4) = ROTS(:,4)-ROTS(4,4);
    % save all in one file
    ROTSALL((6-1)*nsta+1:6*nsta, :) = ROTS;
    
    %%% fam 099
    fam = '099';
    mainsta = 4; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    % save all in one file
    ROTSALL((7-1)*nsta+1:7*nsta, :) = ROTS;
    
    %%% fam 068
    fam = '068';
    mainsta = 5; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    ROTS(:,4) = ROTS(:,4)-ROTS(4,4);
    % save all in one file
    ROTSALL((8-1)*nsta+1:8*nsta, :) = ROTS;
    
    %%% fam 125
    fam = '125';
    mainsta = 5; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas')); 
    ROTS(:,4) = ROTS(:,4)-ROTS(4,4);
    % save all in one file
    ROTSALL((9-1)*nsta+1:9*nsta, :) = ROTS;
    
    %%% fam 147
    fam = '147';
    mainsta = 5; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    ROTS(:,4) = ROTS(:,4)-ROTS(4,4);
    % save all in one file
    ROTSALL((10-1)*nsta+1:10*nsta, :) = ROTS;
    
    %%% fam 017
    fam = '017';
    mainsta = 4; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    % save all in one file
    ROTSALL((11-1)*nsta+1:11*nsta, :) = ROTS;
    
end

%% save resuls of all fam
PREFIX = strcat('rot_para_11fam',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),'LZB');
fid = fopen(strcat(rstpath,PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
fprintf(fid,'%d %d %d %d \n',ROTSALL');
fclose(fid);












