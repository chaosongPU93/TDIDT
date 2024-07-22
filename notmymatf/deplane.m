function [Hxdp,zhat,dist]=deplane(Hx,NyNx,mode)
% Hxdp=DEPLANE(Hx,NyNx,mode)
%
% Calculates the plane that spans the rows of the data space (mode 3) or that 
% best fits the data (mode 1 or 2), and returns the data with the plane removed,
% the plane itself, and the signed L2-norm of the original data to the plane
%
% INPUT:
%
% Hx        Real-valued column vector of spatial-domain topographic data; for
%           current synthetic purposes, we have regularly spaced simulations so
%           specifying elevation and NyNx values are sufficient; however, if we
%           have sparse or irregular data, supply vectorized components in 
%           [Z,Y,X] order
% NyNx      Number of samples in y and x directions; only needed if we do not
%           supply Y and X
% mode      Mode of deplaning
%       1   (fastest for test data) Least squares approach with full calculation
% 			of pseudo-inverse; L.S. design matrix :: [x(:) y(:) ones(:)]; note
% 			that use of OLS implies that we know x and y very well and we assume
%           z is inherently uncertain and z alone can leverage the solution
%       2   Least squares approach: use MATLAB's pinv of L.S. design matrix
%       3   SVD approach: map last right singular vector of data matrix by 
%			data matrix :: [x(:) y(:) z(:)]
%
% OUTPUT:
%
% Hxdp      Real-valued column vector of spatial-domain topographic data; only
%           returns Z component of data
% zhat      Column vector of planar surface fit by LS (mode 1 or 2) or the plane 
%           spanned by first n-1 right singular vectors (mode 3); only returns Z
%           component
% dist      Distance from original data to the plane; calculated as signed 
%           L2 norm
%
% EXAMPLES:
%   rvals = 1:NyNx(1); Hx = randn(NyNx) + rvals; Hxdp = deplane(Hx(:),NyNx,3);
%	figure;mesh(reshape(Hx,NyNx));hold on;surf(reshape(Hxdp,NyNx));hold off
%
% profiling:
%   profile on
%   hxdp1=deplane(Hx,NyNx,1);
%   profile viewer
%   profsave
%   profile off
%
% [Hxdp,zhat,dist]=deplane(Hx,NyNx,3);
% figure();histogram(dist);hold on;histogram(Hxdp)
%
% % If providing non-gridded data Hx with components in [Z,Y,X] order:
% [Hxdp,zhat,dist]=deplane(Hx,[],3);
%
% Last modified 03.28.23 olwalbert-at-princeton.edu

if size(Hx,2)==3
    % note: this hasn't been well tested
    x = Hx(:,3);
    y = Hx(:,2);
    Hx = Hx(:,1);
elseif size(Hx,2)==1
    % Generate a NyNx grid of x and y values
    xv = 1:NyNx(1);
    yv = 1:NyNx(2);
    [X,Y] = meshgrid(xv,yv);
    x = X(:);
    y = Y(:);
else
    error('confirm that data in Hx is provided as (a) column(s)')
end

% Remove mean from grid and topography data
x = x - mean(x);
y = y - mean(y);
z = Hx - mean(Hx);

if mode == 1
	% least squares by mapping design matrix through dual space
	% total time 0.005s for 256 by 256
	G = [x y ones(size(z))]; %design matrix for a plane
	Gd = inv(G'*G)*G';
	b = Gd*z; %the coefficients of the plane, slope in x, slope in y, vertical offset
	zhatp = G*b; %this is the best fit plane of the demeaned data;
    % distance to plane can be found via L2 norm
	Hxdp = z-zhatp;
    zhat = zhatp+mean(Hx);
    dist = sign(Hx-zhat).*sqrt(sum((Hx-zhat).^2,2));
elseif mode == 2
	% now try matlabs pinv
	% total time 0.006s for 256 by 256 test case
	G = [x y ones(size(z))];
	Gp = pinv(G);
	b = Gp*z;
	zhatp = G*b;
	Hxdp = z-zhatp;
    zhat = zhatp+mean(Hx);
    dist = sign(Hx-zhat).*sqrt(sum((Hx-zhat).^2,2));
elseif mode == 3
	A = [x y z];
	% full svd is slow in testing, maybe because rows of A are a power of 2 
    % [~,~,v] = svd(A);

	% Right singular vector that corresponds to smallest singular value  
    % is the normal of an n-dim surface for n-dim data; 
    % for 3-dim data, the surface we seek to fit is planar; the closer the
    % smallest singular value is to 0, the more plane-like the surface will be.

    % We map our data to the last rsv, v, effectively removing the contributions
    % of the planar subspace spanned by the preceding rsvs
    [~,~,v]=svds(A,1,"smallest"); %total time 0.018s for 256 by 256
	Hxdp = A*v;
    % If we are additionally interested in the plane itself or the distance of
    % our data to the plane spanned by the preceding rsvs,
    zhat = z-Hxdp+mean(Hx);
    % or if we had returned addtional vectors, 
    % zhat = @(x,y) x*v(:,1)'+y*v(:,2)'+mean(Hx)
    dist = sign(Hx-zhat).*sqrt(sum((Hx-zhat).^2,2));
end
keyboard