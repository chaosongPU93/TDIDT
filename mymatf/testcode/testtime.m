for i = 1: length(timsSTAperm)-1
    difftime(i)=  timsSTAperm(i+1)-timsSTAperm(i);
end

ind = find(difftime<0);
ind
if ind
    disp('the time is not increasing monotonically')
else
    disp('the time is increasing monotonically')
end