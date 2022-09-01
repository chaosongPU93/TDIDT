%%
% Make up some data
x = rand(1,100);
y = rand(1,100);
z = x.^2 + y.^2;
w = x + y;
% Make a contour plot with the w-data.
% [xg,yg]=meshgrid(x,y) creates a mesh grid that represenst the xy 
% coordinates of all points in the range, which is rectangle or square. 
% Draw a plot on paper would help.
[xg, yg] = meshgrid(0 : 0.1 : 1);   

% vq = griddata(x,y,v,xq,yq) can interpolate 2-D or 3-D scattered data
% based on some scattered data (x,y,w) to some other scattered data or a
% grid array.
wg = griddata(x,y,w,xg,yg);
contour(xg,yg,wg);
% Make a scatter3 plot with the z-data, and also use the z-data to color the markers. The fourth input is empty ( [] ) because we want the markers to all be the same size.
hold on
scatter3(x, y, z, [], z);
% Just for good measure
colorbar

%%
% mesh(X,Y,Z) can create a mesh surface plot in 3-d view, (X,Y,Z) have to
% be grid arrays
[X,Y] = meshgrid(-5:.5:5);
Z = Y.*sin(X) - X.*cos(Y);
figure
s = mesh(X,Y,Z,'FaceAlpha','0.5');
s.FaceColor = 'flat';

%%
x = rand(1,100);
y = rand(1,100);
z = x.^2 + y.^2;
figure
% contour(X,Y,Z) creates a contour plot containing the isolines of matrix
% Z, Z must be a 2-D matrix of value at each point with coordinates 
contour(x,y,z);

% contourf(X,Y,Z) creates a filled contour plot, other syntaxes are pretty
% the same