function timeoffrot = Readbostock(fam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This is a callable function to read Michael Bostock's newest catalog,
% according to the family ID, to exreact the date and other information. 
% 
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/04/19
% Last modified date:   2019/04/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set default value
defval('fam','002')
defval('minnum',40)

% set path and name
catapath= ('/home/data2/chaosong/matlab/allan/BOSTOCK/');
cataname = ('total_mag_detect_0000_cull_NEW.txt');

%%% file format of catalog
% 1. family ID
% 2. date: 2-digit year, 2-digit month, 2-digit day
% 3. hour of thay day, minus 1 to real hour
% 4. secs (<3600)
% 5. magnitude (? NOT sure)
% 6. depth (? NOT sure)
catalog = load(strcat(catapath, cataname));

% get all unique dates
% fam = '002';
famnum = str2double(fam);   
famdet = catalog(famnum == catalog(:, 1), 2);   % all detections in that fam
dateall = unique(famdet);   % unique active dates in that fam
ndate = length(dateall);
icount = 0;
for i = 1: ndate
    ind = find(famdet==dateall(i));
    nlfeday = length(ind);
    if nlfeday >= minnum
        icount = icount+1;
        dateava(icount) = dateall(i);   % available dates, with number of LFE >40
    end
end
dateava = dateava';

% year 2003
date1 = dateava(dateava>=30301 & dateava<=30331);
tmp1 = date1- 30303 + 62;  % 30303 = 2003.062
year1 = zeros(length(date1),1);
jday1 = zeros(length(date1),1);
for i = 1: length(date1)
    year1(i) = 2003;
    jday1(i) = tmp1(i);
    if jday1(i) < 60 || jday1(i) > 63  % LZB has the least data from 060~063
        year1(i) = 9999;
        jday1(i) = 9999;
    end
end

% year 2004
date2 = dateava(dateava>=40701 & dateava<=40731); 
tmp2 = date2- 40714 + 196;    % 40714 =2004.196
year2 = zeros(length(date2),1);
jday2 = zeros(length(date2),1);
for i = 1: length(date2)
    year2(i) = 2004;
    jday2(i) = tmp2(i);
    if jday2(i) < 183 || jday2(i) > 208 % LZB has the least data from 060~063
        year2(i) = 9999;
        jday2(i) = 9999;
    end
end

% year 2005
date3 = dateava(dateava>=50901 & dateava<=50930); 
tmp3 = date3- 50911 + 254;    % 50911 =2005.254
year3 = zeros(length(date3),1);
jday3 = zeros(length(date3),1);
for i = 1: length(date3)
    year3(i) = 2005;
    jday3(i) = tmp3(i);
    if jday3(i) < 244 || jday3(i) > 273  || jday3(i) == 257 % LZB has the least data from 060~063
        year3(i) = 9999;
        jday3(i) = 9999;
    end
end

% combine them
year = [year1; year2; year3];
jday = [jday1; jday2; jday3];

% retain only days that have data
year = year((year ~= 9999));
jday = jday((jday ~= 9999));

timeoffrot = [year jday];
















