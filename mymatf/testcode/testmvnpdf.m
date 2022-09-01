% testmvnpdf.m
% This script serves as a test and reminder of 1-D or higher-dimension normal (Gaussian) 
% distribution and how to draw random numbers from that PDF.
%
% --For 1-D case, FJS has a function 'y=gauss(x,mu,sigma)' (modified by myself) to compute its PDF 
% from equation. After checking, now it is equivalent to matlab's built-in function 
% 'y = normpdf(x,mu,sigma)'. Furthermore, it is also equal to my own function 
% 'G = gauss_time(N,sigma)' in principle.
% --For 2-D case, FJS has a function 'Z=gauss2(X,Y,stdn,stdm)' to compute its PDF from a
% multiplication of two 1d gaussian, requiring that X and Y are independent and come out of
% meshgrid. Matlab has a built-in to deal with dimensions higher than 1, 'mvnpdf', which is able to 
% deal with a covariance matrix, meaning that all dimensions do not have to independent from each
% other. Accordingly, 'mvnrnd' is to draw ramdom numbers from the PDF, 'mvncdf' is to compute the 
% CDF. But do notice that in 'y = mvnpdf(X,mu,Sigma)', 'Sigma' is the covariance, not the standard
% deviation.
% --Recommend to use 'mvnpdf', 'mvncdf' and 'mvnrnd'.
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/04/05
% Last modified date:   2022/04/05


clc;
clear;
close all;

%% Design a 2D normal PDF with 'mvnpdf' that has a particular shape

mu = [0 0]; 
% sigma = eye(2).*4;  % 'eye' creates an identity matrix
% sigma = [4 0; 0 1]; % different variance but i.i.d.
% sigma = [1 0.8; 0.8 1]; % off-diagonal is not 0 

cov = 0.4:0.001:0.8;  %possible covariance
for i = 1: length(cov)
  sigma = [1 cov(i); cov(i) 1]; % off-diagonal is not 0 
  
  %obtain eigenvalue and eigenvector, where the sqrt of eigenvalues are the standard deviation of 
  %the principle axes, and the eigenvectors point to the major and minor axis
  [eigvec,eigval] = eig(sigma);  % covariance matrix 'sigma' has to be positive definite!

  %the std of along the principle axes 
  sigmaeig = [sqrt(eigval(1,1)) sqrt(eigval(2,2))];
  ratio(i) = sqrt(eigval(2,2))/sqrt(eigval(1,1));
end

%the ratio closet to desired determines the best covariance, eigenvalue, eigenvector, etc
[~,ind] = min(abs(ratio-2));
ratbst=ratio(ind);
covbst=cov(ind);
sigma = [1 covbst; covbst 1]; % off-diagonal is not 0 
[eigvec,eigval] = eig(sigma);
sigmaeig = [sqrt(eigval(1,1)) sqrt(eigval(2,2))];

for i = 1:2
  angle(i) = atan2d(eigvec(2,i),eigvec(1,i));
end
180+angle(2)

%A matrix is positive definite if it’s symmetric and all its eigenvalues are positive.
%A matrix is positive definite if it’s symmetric and all its pivots are positive.
if sum(diag(eigval)<0)>0
  disp('Sigma is not positive definite!');  
end

%create meshgird
x1 = -4:0.01:4;
x2 = -4:0.01:4;
[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];

% use the designed covariance matrix to create the PDF of the 2D normal distribution
y = mvnpdf(X,mu,sigma);
y = reshape(y,length(x2),length(x1));

figure
% surf(x1,x2,y,'EdgeColor','flat'); 
% surf(X1,X2,y,'EdgeColor','flat'); % looks like for x and y, you can either use meshgrid or a vector
imagesc(x1,x2,y); hold on;% substitute for surf, but need to make 'y' axis direction as normal
ax = gca; ax.YDir = 'normal';
% contour(x1,x2,y,'k-');
% caxis(ax,[min(y(:))-0.5*range(y(:)),max(y(:))])
axis equal tight
% xlim(minmax(x1));
% ylim(minmax(x2));
colormap('jet');
xlabel('x1');
ylabel('x2');
c=colorbar;
c.Label.String = 'Probability Density';

%plot 1-sigma, 2-sigma and 3-sigma ellipse
x0 = 0; 
y0 = 0;
semia = sigmaeig(2);  %1-sigma
semib = sigmaeig(1);
angrot = 45;
[x, y] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);
plot(x,y,'k-');
semia = 2*sigmaeig(2);  %2-sigma
semib = 2*sigmaeig(1);
[x, y] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);
plot(x,y,'k-');
semia = 3*sigmaeig(2);  %3-sigma
semib = 3*sigmaeig(1);
[x, y] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);
plot(x,y,'k-');

%plot the eigenvector, ie., major and minor axes direction
plot([0,eigvec(1,1)],[0,eigvec(2,1)],'k-','linew',2);
plot([0,eigvec(1,2)],[0,eigvec(2,2)],'k-','linew',2);

%check if it matches the desired shape of ellipse 
semia = 2;
semib = 1;
[x, y] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);
plot(x,y,'w-','linew',2);


%% Draw random numbers from the designed 2D normal PDF with 'mvnrnd'
%%% Note that the drawn numbers has a very high 'resolution' the same as your matlab current
%%% settings, so if you want to bin them, do NOT bin upon pixel, but upon grid point instead.
rng('default')  % For reproducibility
rdnum = mvnrnd(mu,sigma,5000);
scatter(rdnum(:,1), rdnum(:,2),10,[.5 .5 .5],'+');

%% Justify if drawn numbers comply to PDF? 1. plot the density map to see if similar to the data PDF
%%% 1. plot the density map of drawn numbers to see if it is similar to the data PDF
%Obtain the density (count), bin by each pixel (intrinsic resolution)
% density = density_pixel(R(:,1), R(:,2));

%Obtain the density (count), bin by each desired grid point
dx = 0.1;
dy = 0.1;
[density] = density_matrix(rdnum(:,1),rdnum(:,2),minmax(x1),minmax(x2),dx,dy);
emppdf = density;
emppdf(:,3) = emppdf(:,3)/sum(emppdf(:,3))/(dx*dy); % normalize to PDF, ==counts/N_total/area_of_bin
% dx = 0.1;
% dy = 0.1;
% [~,xloc2d,yloc2d,density2d] = density_matrix(R(:,1),R(:,2),minmax(x1),minmax(x2),dx,dy);
% emppdf = density2d;
% emppdf = emppdf/sum(emppdf(:))/(dx*dy); % normalize to PDF, ==counts/N_total/area_of_bin

figure;
hold on
scatter(emppdf(:,1),emppdf(:,2),25,emppdf(:,3),'s','filled'); % use scatter if the density is 1D
% imagesc(minmax(x1),minmax(x2),emppdf); hold on; % use 'imagesc' if the density is 2D matrix
colormap(jet)
c=colorbar;
c.Label.String = 'Probability Density';
caxis([0 max(emppdf(:,3))]);
xlim(minmax(x1));
ylim(minmax(x2));
xlabel('x1');
ylabel('x2');
axis equal tight
box on

%%% An alternative to visualize the PDF is to use 'histogram2', theretically gives the same result 
figure;
hold on
h = histogram2(rdnum(:,1), rdnum(:,2),'BinWidth',[dx dy],'FaceColor','flat','Normalization','pdf');
colormap(jet);
c=colorbar;
c.Label.String = 'Probability Density';
% view(3);
h.DisplayStyle = 'tile';
h.ShowEmptyBins = 'on';
% xlim(minmax(x1));
% ylim(minmax(x2));
xlabel('x1');
ylabel('x2');
box on; grid on;

%% Justify if drawn numbers comply to PDF? 2. obtain the histograms (maybe normalized to PDF estimate)
%%% 2. obtain the histograms (maybe normalized to PDF estimate) along the major and minor axes,
%%% estimate its properties and compare with those of the original PDF
%get the coordinates of the drawn numbers along the major and minor axes, ie. the distance to origin
%in the coordinate system defined by the major and minor axes
[rdnum(:,3), rdnum(:,4)] = coordinate_rot(rdnum(:,1), rdnum(:,2),angle(2));
[muHata,sigmaHata] = normfit(rdnum(:,3));   % along major axis a
synpdfa = normpdf(-5:dx:5,muHata,sigmaHata);
[muHatb,sigmaHatb] = normfit(rdnum(:,4));   % along minor axis b
synpdfb = normpdf(-5:dy:5,muHatb,sigmaHatb);

figure
subplot(121)
histogram(rdnum(:,3),'BinWidth',dx,'Normalization','pdf'); hold on
plot(-5:dx:5,synpdfa,'r-','linew',2);
text(0.1,0.9,sprintf('%.2f',muHata),'Units','normalized','HorizontalAlignment','left',...
  'color','r');
text(0.1,0.8,sprintf('%.2f',sigmaHata),'Units','normalized','HorizontalAlignment','left',...
  'color','r');
text(0.3,0.9,sprintf('%.2f',mu(2)),'Units','normalized','HorizontalAlignment','left');
text(0.3,0.8,sprintf('%.2f',sigmaeig(2)),'Units','normalized','HorizontalAlignment','left');
xlim([-5 5]);
xlabel('Major axis a');
ylabel('Probability Density');

subplot(122)
histogram(rdnum(:,4),'BinWidth',dy,'Normalization','pdf'); hold on
plot(-5:dy:5,synpdfb,'r-','linew',2);
text(0.1,0.9,sprintf('%.2f',muHatb),'Units','normalized','HorizontalAlignment','left',...
  'color','r');
text(0.1,0.8,sprintf('%.2f',sigmaHatb),'Units','normalized','HorizontalAlignment','left',...
  'color','r');
text(0.3,0.9,sprintf('%.2f',mu(1)),'Units','normalized','HorizontalAlignment','left');
text(0.3,0.8,sprintf('%.2f',sigmaeig(1)),'Units','normalized','HorizontalAlignment','left');
xlim([-5 5]);
xlabel('Minor axis b');
ylabel('Probability Density');



