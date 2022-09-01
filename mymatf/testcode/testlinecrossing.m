% this is to test the expression of a straight line, crossing of lines, etc

          
figure
hold on
box on

x1 = -1;
y1 = 2;
theta1 = 245;
x = -10:0.1:10;
% y = 1/tan(deg2rad(theta1)).*(x-x1) + y1;
a1 = 1/tan(deg2rad(theta1));
b1 = y1-a1*x1;
y = linefcn(x,a1,b1);
scatter(x1,y1,'k');
plot(x, y, 'b-');

x2 = 0;
y2 = 4;
theta2 = 180-22.5;
x = -10:0.1:10;
% y = 1/tan(deg2rad(theta2)).*(x-x2) + y2;
a2 = 1/tan(deg2rad(theta2));
b2 = y2-a2*x2;
y = linefcn(x,a2,b2);
scatter(x2, y2,'k');
plot(x, y, 'r-');

% x0 = (b2-b1)/(a1-a2);
% y0 = (b2*a1-b1*a2)/(a1-a2);
% y0 = linefcn(x0,a2,b2);
[x0,y0] = linecrossing(a1,b1,a2,b2);

scatter(x0, y0,'k','filled');

print('-dpdf',strcat(rstpath,'/test1.pdf'));