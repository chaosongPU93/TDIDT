function [timoffrot,stas,bostname,tempoffs] = GetDays(fam,freqflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to get the days and corresponding days and name of
% Bostock's templates according to LFE family name for original PGC trio
%    
% 
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/06/25
% Last modified date:   2019/06/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isequal(fam,'002')  %DON'T FORGET additional family-specific delcarations around line 156

    %%% timeoffrot is the year and date of the data to read
    timoffrot= [
                2003 062; %Could leave out 064 and 257.       % dates of the year, two columns: year Julian_day
                2003 063;
                2004 196;
                2004 197;
                2004 198;
                2004 199;
                2005 254;
                2005 255;
                2005 256];
    
    %%% bostname is the directory of template from Michael Bostock
    bostname=['BOSTOCK/NEW/002-246_2003.062';                  % directory name, path
              'BOSTOCK/NEW/002-246_2003.063';
              'BOSTOCK/NEW/002-246_2004.196';
              'BOSTOCK/NEW/002-246_2004.197';
              'BOSTOCK/NEW/002-246_2004.198';
              'BOSTOCK/NEW/002-246_2004.199';
              'BOSTOCK/NEW/002-246_2005.254';
              'BOSTOCK/NEW/002-246_2005.255';
              'BOSTOCK/NEW/002-246_2005.256'];
    
%     %%% stas = used names of 3-station pair
%     stas=['PGC  '         % NOTICE the blank spaces, stas have 5 columns
%           'SSIB '
%           'SILB '];
        
    if isequal(freqflag,'hf')
        tempoffs=[1210 1210 1210];  % 1.25-6.5hz, 20sps, PGC trio
    elseif isequal(freqflag,'lf')
        tempoffs=[606 607 607];  % 0.5-1.25hz, 20sps, PGC trio
    end
    
    stas=['PGC  '
        'SSIB '
        'SILB '
        'LZB  '
        'TWKB '
        'MGCB '
        'KLNB '];
        
    if isequal(freqflag,'hf')
        tempoffs=[1197 1284 1218 1208 1195 1175 1194];  % 1.25-6.5hz, 40sps, PGC trio
    elseif isequal(freqflag,'lf')
        tempoffs=round([1197 1284 1218 1208 1195 1175 1194]/2);  % 0.5-1.25hz, 20sps, PGC trio
    end      

    
elseif isequal(fam,'068')
    
    timoffrot=[
               2004 198;
               2004 199;
               2004 200;
               2004 201;
               2005 256;
               2005 257;
               2005 258;
               2005 259;
               2005 260;
               2005 261];
    bostname=['BOSTOCK/NEW/068_2004.198';       % what is the data structure of these file??
              'BOSTOCK/NEW/068_2004.199';
              'BOSTOCK/NEW/068_2004.200';
              'BOSTOCK/NEW/068_2004.201';
              'BOSTOCK/NEW/068_2005.256';
              'BOSTOCK/NEW/068_2005.257';
              'BOSTOCK/NEW/068_2005.258';
              'BOSTOCK/NEW/068_2005.259';
              'BOSTOCK/NEW/068_2005.260';
              'BOSTOCK/NEW/068_2005.261'];
          
    stas=['TWKB '
          'LZB  '
          'MGCB '];
    
    if isequal(freqflag,'hf')
        tempoffs=[902 894 901]; %these are the zero crossings for 068: TWKB,LZB,MGCB.
    elseif isequal(freqflag,'lf')
        tempoffs=[902 894 901]; %these are the zero crossings for 068: TWKB,LZB,MGCB.
    end

        
end
