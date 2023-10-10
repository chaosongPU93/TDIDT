function [density1d,indgrid,xgrid,ygrid,densitygrid] = density_matrix(x,y,xran,yran,dx,dy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [density1d,xloc2d,yloc2d,density2d] = density_matrix(x,y,xran,yran,dx,dy)
% 
% This function is to obtain the density (number of points;
% repetition times) of a 2-D data set (x,y) in each
% rectangluar grid. Suitable to plot the data density
% in terms of binning based on the grid, instead of pixels.
% Use 'density_pixel' if you want to bin upon pixels. 
%       
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/12
% Last modified date:   2022/03/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nx = (xran(end) - xran(1))/ dx;
ny = (yran(end) - yran(1))/ dy;

% xloc1d = zeros(nx*ny,1);
% yloc1d = zeros(nx*ny,1);
density1d = zeros(nx*ny,3);   %vector form of the count inside each grid cell
indgrid = cell(nx*ny,1);  %indices of data inside each grid cell 
xgrid = nan(nx,ny); %grid length x
ygrid = nan(nx,ny); %grid length y
densitygrid = nan(nx,ny); % grid form of the count inside each grid cell
k=1;
for i = 1: nx
    for j = 1: ny
        density1d(k,1) = xran(1)+ (i-1+0.5)*dx;
        density1d(k,2) = yran(1)+ (j-1+0.5)*dy;
        indgrid{k} = find(x>= xran(1)+ (i-1)*dx & x< xran(1)+ i*dx &...
                           y>= yran(1)+ (j-1)*dy & y< yran(1)+ j*dy);                 
%         density1d(k,3) = sum(x>= xran(1)+ (i-1)*dx & x< xran(1)+ i*dx &...
%                            y>= yran(1)+ (j-1)*dy & y< yran(1)+ j*dy);
        density1d(k,3) = length(indgrid{k});
        k = k+1;
        
%         xloc2d(i,j) = xran(1)+ (i-1+0.5)*dx;
%         yloc2d(i,j) = yran(1)+ (j-1+0.5)*dy;
%         tmp = sum(x>= xran(1)+ (i-1)*dx & x< xran(1)+ i*dx &...
%                   y>= yran(1)+ (j-1)*dy & y< yran(1)+ j*dy);
%         density2d(i,j) = tmp;
        
        xgrid(j,i) = xran(1)+ (i-1+0.5)*dx;
        ygrid(j,i) = yran(1)+ (j-1+0.5)*dy;
        tmp = sum(x>= xran(1)+ (i-1)*dx & x< xran(1)+ i*dx &...
                  y>= yran(1)+ (j-1)*dy & y< yran(1)+ j*dy);
        densitygrid(j,i) = tmp;
    end
end

ygrid = flipud(ygrid);
densitygrid = flipud(densitygrid);



