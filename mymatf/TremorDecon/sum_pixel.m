function [sumz1d,indices] = sum_pixel(x,y,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [sumz1d,indices] = sum_pixel(x,y,z)
%
% Beyond 'density_pixel', sometimes you not only want the cumulative count
% of data (x,y) at each pixel, but also you want the sum of z at those 
% unique points. The sum of z and indices of unique pixels are returned.
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2023/05/27
% Last modified date:   2023/05/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[density1d,indices] = density_pixel(x,y);
npts = size(density1d,1);
sumz1d = density1d;
for i = 1: npts
  sumz1d(i,3) = sum(z(indices{i}));
end