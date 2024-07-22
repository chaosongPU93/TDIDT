function covar = mvncov(sigmax,sigmay,rotang)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% covar = mvncov(sigmax,sigmay,rotang)
%
% This function is used to compute the covariance matrix for 2D Gaussian 
% distribution that is required as 'SIGMA' in functions like 'mvnpdf', 
% 'mvnrnd', etc. Suppose you want an 2D Gaussian PDF such that the sigma in
% x direction is 'sigmax' and the sigma in y direction is 'sigmay', and the 
% whole PDF is rotated counterclockwise by 'rotang' degrees. This function
% will return the corresponding covariance matrix between the final 2D 
% Gaussian PDF.
%
% INPUT:  
%   sigmax: sigma along x direction
%   sigmay: sigma along y direction
%   rotang: rotation angle counter-clockwise. in degree, add '-' to incate 
%           clockwise rotation; It controls the off-diagnol values
% 
% OUTPUT:
%   covar:  covariance matrix that can be used to generate 2D Gaussian PDF
%           symmetric, positive semi-definite
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/03/14
% Last modified date:   2024/03/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%scale matrix determined by sigx, sigy
S = [sigmax^2 0; 0 sigmay^2];
%for rotation angle, positive means counterclockwise
rotrad = deg2rad(rotang);
%rotation matrix 
R = [cos(rotrad) -sin(rotrad); sin(rotrad) cos(rotrad)];
%covariance, symmetric, positive semi-definite
covar = R*S*R';

%check if the rotation is succesful using eigenvalue decomposition
[~,eigval] = eig(covar);
sigmaeig = [sqrt(eigval(1,1)) sqrt(eigval(2,2))];
if abs(mean(minmax([sigmax sigmay])-sigmaeig)) > 1e-5
  error('Covariance cannot be computed');
end




