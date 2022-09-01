% test gaussian fitting
% close all

sigma = 4;
hsize = 2*ceil(2.*sigma)+1;
% hsize = [2*ceil(2.*sigma)+1 2*ceil(5.*sigma)+1];
h = fspecial('gaussian',hsize,sigma);
f=figure;
% f.fig.Renderer = 'painters';
ax = gca;
surf(ax,h)
colorbar
xlabel(ax,'x')
ylabel(ax,'y')
zlabel(ax,'z')

