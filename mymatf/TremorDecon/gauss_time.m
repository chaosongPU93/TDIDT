function G = gauss_time(N,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function G = gauss_time(N,sigma)
% Create a symmetric N-point gaussian window in time with 
% a standard deviation of 'sigma', all in samples. See also 'gauss_zerophase' 
%
% INPUT:
%   N: length of window
%   sigma: standard deviation of the gaussian, default is 1
%
% OUTPUT:
%   G: Gaussian window
%
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/11/22
% Last modified date:   2021/11/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('sigma',1);

G = zeros(N,1);   % predefine the matrix
if mod(N,2) == 0  % even number
  x = (1:1:N/2)';
  temp = exp(-x .^ 2 / (2 * sigma ^ 2));
  G(N/2+1: end) = temp;
  G(1: N/2) = flipud(temp);
else  % odd numer
  x = (-floor(N/2):1:floor(N/2))';
  G = exp(-x .^ 2 / (2 * sigma ^ 2));
end

%normalize to make the area be 1
G = G / sum(G);
