% the boundary is more generalized version which can indicate the shrink factor s, whose default is
% 0.5, if s=0, it is equal to convex hull.
[k1,v1] = boundary(obara(:,6),obara(:,7),0.5);
figure
hold on
scatter(obara(:,6),obara(:,7),6,[0.4 0.4 0.4],'filled','o');
scatter(obara(k1,6),obara(k1,7),20,'b','filled','o');
p1=plot(obara(k1,6),obara(k1,7),'b');

% convex hull gives the biggest boundary, boundary is preferred
[k2,v2] = convhull(obara(:,6),obara(:,7));
scatter(obara(k2,6),obara(k2,7),20,'r','filled','o');
p2=plot(obara(k2,6),obara(k2,7),'r');

% creat the polygon defined by the boundary points
tic
pgon = polyshape(obara(k1,6),obara(k1,7));
v3 = polyarea(obara(k1,6),obara(k1,7));
x=[136, 136.4];
y=[34, 33.8];
p3=scatter(x,y,40,'k','filled','o');
% check if the points are inside the polygon
is = isinterior(pgon,x,y);
toc

% directly query if the points are inside or on a boundary
% inpolygon is faster than isinterior
tic
[is2,ion2] = inpolygon(x,y,obara(k1,6),obara(k1,7));
toc

legend([p1,p2,p3],{'function boundary','function convhull','query points'},'location','northwest');
