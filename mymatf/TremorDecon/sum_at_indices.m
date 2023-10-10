function sumz = sum_at_indices(z,indices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sumz = sum_at_indices(z,indices)
%
% Beyond 'density_pixel', sometimes you not only want the cumulative count
% of data z at each pixel, but also you want the sum of z at those 
% unique points. The sum of z is returned.
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2023/05/27
% Last modified date:   2023/05/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

npts = size(indices,1);
sumz = zeros(npts,1);
for i = 1: npts
  sumz(i,1) = sum(z(indices{i}));
end