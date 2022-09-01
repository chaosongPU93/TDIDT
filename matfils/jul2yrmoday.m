function yrmodom=jul2yrmodom(jday,year)
% written by Allan Rubin

daysinmonth=[31 28 31 30 31 30 31 31 30 31 30 31];
if abs(year-2004)<=1.e-6
    daysinmonth=[31 29 31 30 31 30 31 31 30 31 30 31];
end
%MONTHS=char('JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC');
cumdays=cumsum(daysinmonth);
mo=find(jday<=cumdays,1);
dom=jday-cumdays(mo-1);
MONTHS(find(jday<=cumsum(daysinmonth),1),:);
