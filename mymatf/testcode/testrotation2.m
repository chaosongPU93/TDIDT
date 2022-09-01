clc
close all
clear

%%
ns = [-3 3];
ew = [-6 2];
x0 = mean(ew);
y0 = mean(ns);
wid = 8;
hgt = 6;
%just a rectangle, no rotation
[xrec,yrec] = rectangle_chao(x0,y0,wid,hgt,0.01);

% rotate it c-clock by 30 deg, relative to 0,0
[xr1,yr1] = complex_rot(xrec,yrec,30,[0 0]);

% rotate c-clock by 30 deg, relative to -6, -3
[xr2,yr2] = complex_rot(xrec,yrec,30,[-6 -3]);

figure
plot(xrec,yrec,'k-','linew',1); hold on
scatter(0,0,15,'ko','filled'); % origin in X-Y-Z
plot(xr1,yr1,'r-','linew',1);  % projected object in X-Y-Z
plot(xr2,yr2,'b-','linew',1);  % projected object in X-Y-Z
axis equal

%%
ns = [-3 3];
ew = [-6 2];
x0 = mean(ew);
y0 = mean(ns);
wid = 8;
hgt = 6;
%just a rectangle, no rotation
[xrec,yrec] = rectangle_chao(x0,y0,wid,hgt,0.01);

%a rectangle, and rotate c-clock by 30 deg, relative to 0,0
[xr1,yr1] = rectangle_chao(x0,y0,wid,hgt,0.01,30,[0 0]);

%a rectangle, and rotate c-clock by 30 deg, relative to -6, -3
[xr2,yr2] = rectangle_chao(x0,y0,wid,hgt,0.01,30,[-6 -3]);

figure
plot(xrec,yrec,'k-','linew',1); hold on
scatter(0,0,15,'ko','filled'); % origin in X-Y-Z
plot(xr1,yr1,'r-','linew',1);  % projected object in X-Y-Z
plot(xr2,yr2,'b-','linew',1);  % projected object in X-Y-Z
axis equal

%%
semia = 3;
semib = 2;
x0 = 1;
y0 = 2;
%just an ellipse, no rotation
[xell, yell] = ellipse_chao(x0,y0,semia,semib,0.01);

%an ellipse, and rotate c-clock by 30 deg, relative to x0,y0
[xr1,yr1] = ellipse_chao(x0,y0,semia,semib,0.01,30,[x0,y0]);
  
%an ellipse, and rotate c-clock by 30 deg, relative to 0,0
[xr2,yr2] = ellipse_chao(x0,y0,semia,semib,0.01,30,[0 0]);

%an ellipse, and rotate c-clock by 30 deg, relative to 0,0
[xr3,yr3] = ellipse_chao(x0,y0,semia,semib,0.01,30,[-4 2]);

figure
plot(xell,yell,'k-','linew',1); hold on
scatter(0,0,15,'ko','filled'); % origin in X-Y-Z
plot(xr1,yr1,'r-','linew',1);  % projected object in X-Y-Z
plot(xr2,yr2,'b-','linew',1);  % projected object in X-Y-Z
plot(xr3,yr3,'c-','linew',1);  % projected object in X-Y-Z
axis equal
  
  
  
  
  
  
  
  
  
  
  