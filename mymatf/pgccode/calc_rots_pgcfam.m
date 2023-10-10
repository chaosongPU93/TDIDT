%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to call the function calc_rots_noresp to compute the 
% rotation parameters for fam '243' with PGC trio
%       
%   reference to the calc_rots_allfam.m  in */lzbcode    
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
%   2021/06/22
%   The ROT param for fam 002 at PGC trio was obtained by a older version of
%   dates and a older catalog of LFEs which only contains fam 002,
%   in this function, we aim to
%   generically and automatically get the dates and lfe timings for any fams
%   at either PGC or TWKB trio according to the updated catalog of LFEs
%   for stacking ('/BOSTOCK/total_mag_detect_0000_cull_NEW.txt'). Correspondingly,
%   the dates for detections would be either similar to 23 dates in TKB trio, or
%   to 9 dates for fam 022 in PGC trio, depending on the number of detections
%   on that date 
%
%   NOTES:
%   '243';   %pretty
%   '253';   %sta TWKB is poor, rest is good
%   '036';   %pretty good
%   '034';   %prefer a bit 036 over 034, pgc trio and klnb is slightly better at 034, but sta TWKB
%             has a diff polarity!, sta SSIB's splitting reached bound, constraints could fix it
%   '061';   %sta TWKB seems to have a diff polarity!, sta LZB is poor
%   '023';   %the twkb trio are poor, pgc trio and klnb are fine
%   '055';   %should prefer 023 over 055, the twkb trio are poor, pgc trio and klnb are fine
%   '026';   %glitches at sta twkb, diff polarity at mgcb, generally poor at twkb trio
%   '251';   %should prefer 251 over 026 and 162, twkb trio is poor, klnb and pgc trio is fine
%   '162';   %need only 1 among 026,251,162, twkb trio is poor, klnb and pgc trio is fine
%   '240';   %twkb trio is poor, twkb has diff polarity, klnb and pgc trio is fine
%   '250';   %use 240 not 250, sta other than pgc trio are unrecognizable, ssib has diff polarity
%   '255';   %12,255,65 is a bit closer to twkb trio, lzb is poor, klnb and pgc has diff polarity
%   '012';   %all are fine, but pgc,ssib,twkb seem to have diff polarity, keep 012 and 255 both
%   '065';   %split time at lzb is too large reaching bound, at ssib and mgcb is big, ssib, lzb,
%             mgcb are poor

%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/12
% Last modified date:   2021/06/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
% close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]
wid=scrsz(3);       % width
hite=scrsz(4);       % height

% set path
workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath,'/split_chao/pgctrio/');

% trial fams, not necessarily the final used ones
fampool = [
%            '002';
%            '243'; 
%            '253';
%            '036';
%            '034'; 
%            '061';
%            '023';
%            '055';
%            '026';
%            '251';
%            '162';
%            '240';
%            '250';
%            '255'; 
%            '012'; 
%            '065'; 
%            '047';
           '246';
           ];

nfam = length(fampool);
        
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
%     splitoff = 12;  % the twkb trio are poor, pgc trio and klnb are fine
%     staoff = 8;
%     mainsta = 1;
% %     lo = 1.25;
% %     hi = 6.5;
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas_newcat'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
%     keyboard

%     %%% fam 243
%     fam = '243';
%     splitoff = 12;
%     staoff = 8;
%     mainsta = 1; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);

%     %%% fam 253
%     fam = '253';
%     splitoff = 12;
%     staoff = 8;
%     mainsta = 1; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);

%     %%% fam 36
%     fam = '036';
%     splitoff = 12;
%     staoff = 8;
%     mainsta = 1; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);

%     %%% fam 34
%     fam = '034';    % sta TWKB has a wrong polarity!, sta SSIB's splitting reached bound!
%     splitoff = 3;
%     staoff = 9;
%     mainsta = 1; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
    
%     %%% fam 61
%     fam = '061';    % sta TWKB has a wrong polarity!, sta LZB is poor
%     splitoff = 12;
%     staoff = 8;
%     mainsta = 1; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
    
%     %%% fam 23
%     fam = '023';
%     splitoff = 12;  % the twkb trio are poor, pgc trio and klnb are fine
%     staoff = 15;
%     mainsta = 1; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
%     keyboard
    
%     %%% fam 55
%     fam = '055';   % should prefer 023 over 055, the twkb trio are poor, pgc trio and klnb are fine
%     splitoff = 12;
%     staoff = 15;
%     mainsta = 1; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
    
%     %%% fam 26
%     fam = '026';    %glitches at sta twkb, wrong polarity at mgcb
%     splitoff = 12;
%     staoff = 8;
%     mainsta = 1; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
    
%     %%% fam 251
%     fam = '251';    % twkb trio is poor, klnb and pgc trio is fine
%     splitoff = 10;
%     staoff = 8;
%     mainsta = 1; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
        
%     %%% fam 162
%     fam = '162';
%     splitoff = 12;  % twkb trio is poor, klnb and pgc trio is fine
%     staoff = 15;
%     mainsta = 1; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
    
%     %%% fam 240
%     fam = '240';    % twkb trio is poor, klnb and pgc trio is fine
%     splitoff = 12;  
%     staoff = 8;
%     mainsta = 1; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);

%     %%% fam 250
%     fam = '250';    % prefer 240 over 250, sta other than pgc trio are poor, ssib has diff polarity
%     splitoff = 12;
%     staoff = 15;
%     mainsta = 1; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
    
%     %%% fam 255
%     fam = '255';    % twkb trio is poor, klnb and pgc has diff polarity
%     splitoff = 12;
%     staoff = 15;
%     mainsta = 1; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);

%     %%% fam 12
%     fam = '012';    % all are fine, but pgc,ssib,twkb seem to have diff polarity
%     splitoff = 8;
%     staoff = 15;
%     mainsta = 1; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
    
%     %%% fam 65
%     fam = '065';
%     splitoff = 13;  % split time at lzb is too large reaching bound, at ssib and mgcb is big
%     staoff = 8;
%     mainsta = 1; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
%     keyboard
    
%     %%% fam 047
%     fam = '047';
%     splitoff = 12;  % the twkb trio are poor, pgc trio and klnb are fine
%     staoff = 8;
%     mainsta = 1;
% %     lo = 1.25;
% %     hi = 6.5;
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
    
    %%% fam 246
    fam = '246';
    splitoff = 12;  % the twkb trio are poor, pgc trio and klnb are fine
    staoff = 8;
    mainsta = 1;
%     lo = 1.25;
%     hi = 6.5;
    [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
    % save all in one file
    ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
    
    % save results for each family
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas_newcat'),'w+');
    fprintf(fid,'%d %d %d %d \n',ROTS');
    fclose(fid);
    keyboard

else   
    %%% fam 243
    fam = '243';
    splitoff = 12;
    staoff = 8;
    mainsta = 1; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    % save all in one file
    ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
    
end







