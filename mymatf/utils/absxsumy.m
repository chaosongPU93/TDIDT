function [absx, sumy] = absxsumy(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [absx, mergey] = absxmergey(x,y)
% This function is to first get the absolute value of x, for those y
% which shares the same abs(x), sum them. Technically, either abs(x)
% is unique, or shared by at most 2 data points
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2024/03/05
% Last modified date:   2024/03/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=[abs(x) y];
absx = unique(data(:,1));
sumy = zeros(length(absx),1);
for i = 1: length(absx)
  sumy(i)=sum(data(data(:,1)==absx(i),2));
end