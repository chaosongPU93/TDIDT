function [medz, meanz] = median_at_indices(z,indices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [medz, meanz] = median_at_indices(z,indices)
%
% Similar to 'sum_at_indices', sometimes you also want the median or mean 
% of data z at each pixel.
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2023/09/12
% Last modified date:   2023/09/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  npts = size(indices,1);
  medz = zeros(npts,1);
  meanz = zeros(npts,1);
  for i = 1: npts
    medz(i,1) = median(z(indices{i}));
    meanz(i,1) = mean(z(indices{i}));
  end