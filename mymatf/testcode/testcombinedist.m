% testcombinedist
%
% Imagine you have a bunch of sources coming a from a region, the source location
% follows a 2D Gaussian distribution (starting point), and the origin time has a
% unit distribution within a time period. Now we want to know what the distribution
% like at different stations. Say the reference station is PGC.
%   
% 
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/04/11
% Last modified date:   2022/04/11

clear;
clc;
close all

%% design a custom 2D normal PDF using 'mvnpdf'
x1 = -10:0.1:10;
x2 = -10:0.1:10;
[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];

% use the designed covariance matrix to create the PDF of the 2D normal distribution
mu = [0 0]; 
sigma = [4 2; 2 4]; % off-diagonal is not 0 
[eigvec,eigval] = eig(sigma);
sigmaeig = [sqrt(eigval(1,1)) sqrt(eigval(2,2))];

y = mvnpdf(X,mu,sigma);
y = reshape(y,length(x2),length(x1));

figure;
hold on
imagesc(x1,x2,y)
colormap(jet)
c=colorbar;
c.Label.String = 'Probability Density';
xlabel('x1')
ylabel('x2')
axis equal tight
box on

%% Draw random samples in locations from the designed 2D normal PDF with 'pinky'
rng('default')  % For reproducibility
N=10000;
rdnum=zeros(N,2);
for i=1:N
  [rdnum(i,1),rdnum(i,2)]=pinky(x1,x2,y);
end

scatter(rdnum(:,1),rdnum(:,2),10,[.5 .5 .5],'+');

%% Draw random samples in time from uniform distribution
tori = 10*rand(N,1);
figure
histogram(tori,'Normalization','pdf');
title('Origin time');
xlabel('Time (s)');
ylabel('PDF');

%%
gridpath = strcat(getenv('MHOME'), '/Seisbasics/hypoinverse/forsummary');
gridfile = 'evtloc.offset_002_rectgrid_13_20.mit.ref5';
offgrido = load(fullfile(gridpath, gridfile)); % 2-sample spacing at 40 Hz
mm=27; nn=27;
lon=reshape(offgrido(:,1),mm,nn); %degrees
lat=reshape(offgrido(:,2),mm,nn); %degrees
dep=reshape(offgrido(:,3),mm,nn); %km
% tori=reshape(offgrido(:,4),mm,nn); %origin time, + travel time = 30 s
off12o=reshape(offgrido(:,7),mm,nn); %PGSS is uniform in columns; changes across rows.
off13o=reshape(offgrido(:,8),mm,nn); %PGSI is uniform in rows; changes down columns.

%center, location of 002
lon0 = interp2(off12o,off13o,lon,0,0,'spline');  % interpolate as well
lat0 = interp2(off12o,off13o,lat,0,0,'spline');


[dx,dy] = absloc2relaloc(offgrido(:,1),offgrido(:,2),lon0,lat0); % convert to rela loc around 002
xo = reshape(dx,mm,nn);
yo = reshape(dy,mm,nn);

off12n = griddata(dx,dy,offgrido(:,7),rdnum(:,1),rdnum(:,2),'cubic');
off13n = griddata(dx,dy,offgrido(:,8),rdnum(:,1),rdnum(:,2),'cubic');
torin = griddata(dx,dy,offgrido(:,4),rdnum(:,1),rdnum(:,2),'cubic');

tarrs1 = tori+(30-torin);
tarrs2 = tarrs1-off12n/20+offgrido(1,5)/40;
tarrs3 = tarrs1-off13n/20+offgrido(1,6)/40;

figure
subplot(131)
histogram(tarrs1,'Normalization','pdf');
title('Arrival at sta 1');
xlabel('Time (s)');
ylabel('PDF');
subplot(132)
histogram(tarrs2,'Normalization','pdf');
title('Arrival at sta 2');
xlabel('Time (s)');
subplot(133)
histogram(tarrs3,'Normalization','pdf');
title('Arrival at sta 3');
xlabel('Time (s)');

%%
figure
subplot(131)
histogram(off12n,'Normalization','pdf');
title('Offset 12');
xlabel('Samples at 20 sps');
ylabel('PDF');
subplot(132)
histogram(off13n,'Normalization','pdf');
title('Offset 13');
xlabel('Samples at 20 sps');
subplot(133)
histogram(30-torin,'Normalization','pdf');
title('Travel time');
xlabel('Samples at 20 sps');












