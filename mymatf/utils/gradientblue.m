function c = gradientblue(m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GRADIENTBLUE    gradient from light to drak blue color map
%   GRADIENTBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with light blue, and then progrades to dark blue.
%   GRADIENTBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one. 
%   Code is from Allan Rubin found online.
%
%   For example, to reset the colormap of the current figure:
%
%       colormap(gradientblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.
%
% Chao Song, chaosong@princeton.edu
% 
% First created date:   2024/04/25
% Last modified date:   2024/04/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1, m = size(get(gcf,'colormap'),1); end

lightblue=[0.356863,0.811764,0.956863];
darkblue=[0.019608,0.074510,0.670588];

c = zeros(m,3);
for i=1:m
  c(i,:) = lightblue+(darkblue-lightblue)*((i-1)/(m-1));
end

%     solidred=[1,0,0];
%     solidblue=[0,0,1];
%     red2blue=@(i,nbins) solidblue+(solidred-solidblue)*((i-1)/(nbins-1));