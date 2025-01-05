function [yr,mo,dy]=yrjday2datestr(yrjday)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the function to convert a date number in the format of 'yyyyddd'
% to the date strings in the format of year, month, and day of the month.
% In 'yyyyddd', yyyy is the 4-digit year, and ddd is the 3-digit julian date.
%  
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/12/18
% Last modified date:   2024/12/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

year = floor(yrjday/1000);
jday = floor(yrjday-year*1000);
a = jul2dat(year,jday);
if a(1) == 9
  mo = 'Sep.';
elseif a(1) == 7
  mo = 'Jul.';
else
  mo = 'Mar.';
end
dy = num2str(a(2));
yr = num2str(a(3));
