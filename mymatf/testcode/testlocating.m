% test locating algorithm
close all
figure
hold on
[x1, y1] = circle_chao(1,2,2,0.1);
[x2, y2] = circle_chao(3,5,2.5,0.1);
[x3, y3] = circle_chao(2,3,sqrt((2-0.7233)^2+(3.967-3)^2),0.1);
plot(x1, y1,'color','b');
plot(x2, y2,'color','b');
plot(x3, y3,'color','b');

[x1, y1] = circle_chao(1,2,2+1,0.1);
[x2, y2] = circle_chao(3,5,2.5+1,0.1);
[x3, y3] = circle_chao(2,3,sqrt((2-0.7233)^2+(3.967-3)^2)+1,0.1);
plot(x1, y1,'color','r');
plot(x2, y2,'color','r');
plot(x3, y3,'color','r');

axis equal


