clear;
clc;
close all;

figure;
ax = gca;

t=15;
f=1:0.01:25;
Qmin=200;
ymin=exp(-pi*f*t/Qmin);
plot(ax,f,ymin);
hold on

Qmax=500;
ymax=exp(-pi*f*t/Qmax);
plot(ax,f,ymax);

patarea = [f' ymin';
           flipud(f') flipud(ymax');
           f(1) ymin(1)];
patch(ax,patarea(:,1),patarea(:,2),'y','Facealpha',0.3,'edgecolor','none');

for i = 1: length(f)-1
  dymaxdf(i) = (log10(ymax(i+1))-log10(ymax(i))) / (log10(f(i+1))-log10(f(i)));
end
[slope, ind] = min(abs(dymaxdf+1));
dymaxdf(ind)
plot(ax,f(ind),ymax(ind),'k.','markersize',20);

for i = 1: length(f)-1
  dymindf(i) = (log10(ymin(i+1))-log10(ymin(i))) / (log10(f(i+1))-log10(f(i)));
end
[slope, ind] = min(abs(dymindf+1));
dymindf(ind)
plot(ax,f(ind),ymin(ind),'k.','markersize',20);

axis(ax,[1 1e2 1e-2 1]);
ax.XScale = 'log';
ax.YScale = 'log';
% axis equal;

