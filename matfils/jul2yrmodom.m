function yrmodom=jul2yrmodom(jday,year)

daysinmonth=[31 28 31 30 31 30 31 31 30 31 30 31];
if year==4
    daysinmonth=[31 29 31 30 31 30 31 31 30 31 30 31];
end
cumdays=cumsum(daysinmonth);
mo=find(jday<=cumdays,1);
dom=jday-cumdays(mo-1);
yrmodom=year*10000+mo*100+dom;