% testrand_custompdf_p1.m
%
% This script serves as a test and reminder of obtaining the custom PDF/CDF from a dataset, not
% necessarily a Gaussian (normal) distribution. Given the PDF/CDF, we try to explore how to draw 
% pseudo-random numbers from it.
% 
% Part 1, a custom PDF that is explicitly expressed as an equation; continuous; analytic, etc.
%
% --For 1-D case, if you have a data set 'data' and weights 'wt' and you want to draw random
%   numbers with the same dimension size of the summed weights, by N times (realizations),
%   from its underlying PDF. You can choose 'randgen_pdf', 'randgen_ksdensity', or 
%   'randgen_norm', depending on your needs and if you have the pdf already. The function 
%   'randgen_custom' summarizes the 3 methods.
% --Another option is here: 
%   https://www.mathworks.com/matlabcentral/answers/650778-generate-random-number-from-custom-pdf-and-cdf
%   which creates a PiecewiseLinear distribution or empirical distribution using 'makedist' and then
%   draw random numbers from it using 'random' provided in the statistics and ML toolbox.
%
% --For 2-D case, 'pinky(x1,x2,y)' works fine which has already been tested on a 2D normal 
%   distribution in which x and y are not necessarily independent from each other. You need to know
%   the PDF value at x1,x2, then pinky can draw random numbers from this PDF. The PDF can be
%   written in equation, or only in discrete form. If you don't know about the PDF of a discrete
%   data set at all, ie. you only have x1 and x2 vectors (but from them you have a feeling of 
%   counts), you have a few options. 
%   Option 1 is to bin the x1 and x2 upon pixel (ie. find unique entries) or to bin upon grid point
%   with a size. This is similar to getting the density distribution from x1 and x2. Normalizing the
%   density (counts), you can obtain PDF. Then you can feed this PDF to 'pinky' to draw samples from
%   it.
%   Similarly to option 1, option 2 is to use 'histogram2' to bin for you, specify the binedges
%   carefully, and normalize it to PDF, then you can get the PDF as the height and histogram. From
%   plotting, you could see that they are almost identical.
%   Option 3 is '[f,xi] = ksdensity(x,pts)' where x is the 2-colomn vector [x1 x2], pts is the
%   location points that you want to estimate the PDF (output f) on. Then you can feed this discrete 
%   PDF to 'pinky' to draw samples from it. (help ksdensity and see the last example). But note that
%   the 'function' can NOT deal with 'icdf', meaning that if can NOT draw 2D random numbers using
%   this funtion directly.
%   
% 
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/04/06
% Last modified date:   2022/04/07

clear;
clc;
close all

%% design a custom 2D normal PDF using 'mvnpdf'
x1 = -4:0.01:4;
x2 = -4:0.01:4;
[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];

% use the designed covariance matrix to create the PDF of the 2D normal distribution
mu = [0 0]; 
sigma = [1 0.6; 0.6 1]; % off-diagonal is not 0 
[eigvec,eigval] = eig(sigma);
sigmaeig = [sqrt(eigval(1,1)) sqrt(eigval(2,2))];

y = mvnpdf(X,mu,sigma);
y = reshape(y,length(x2),length(x1));

%% Draw random numbers from the designed 2D normal PDF with 'pinky'
rng('default')  % For reproducibility
N=5000;
rdnum=zeros(N,2);
for i=1:N
  [rdnum(i,1),rdnum(i,2)]=pinky(x1,x2,y);
end

figure;
hold on
imagesc(x1,x2,y)
colormap(jet)
c=colorbar;
c.Label.String = 'Probability Density';
scatter(rdnum(:,1),rdnum(:,2),10,[.5 .5 .5],'+');
xlabel('x1')
ylabel('x2')
axis equal tight
box on

%% Justify if drawn numbers comply to PDF? 1. plot the density map to see if similar to the data PDF
%%% 1. plot the density map of drawn numbers to see if it is similar to the data PDF
%Obtain the density (count), bin by each pixel (intrinsic resolution)
% density = density_pixel(R(:,1), R(:,2));

%Obtain the density (count), bin by each desired grid point
dx = 0.1;
dy = 0.1;
[density,xgrid,ygrid,densitygrid] = density_matrix(rdnum(:,1),rdnum(:,2),minmax(x1),minmax(x2),dx,dy);
emppdf = density;
emppdf(:,3) = emppdf(:,3)/sum(emppdf(:,3))/(dx*dy); % normalize to PDF, ==counts/N_total/area_of_bin
emppdfgrid = densitygrid/sum(sum(densitygrid))/(dx*dy);

%visualize the empirical PDF
figure;
hold on
scatter(emppdf(:,1),emppdf(:,2),25,emppdf(:,3),'s','filled'); % use scatter if the density is 1D
% imagesc(xgrid(1,:)',ygrid(:,1),emppdfgrid); hold on; % use 'imagesc' if the density is 2D matrix
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
[rdnum(:,3), rdnum(:,4)] = coordinate_rot(rdnum(:,1), rdnum(:,2),45);
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
















