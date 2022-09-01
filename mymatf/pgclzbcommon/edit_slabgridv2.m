function edit_slabgridv2
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% version 2 of editting the slab1.0 grid file to make sure that it
% is compatible to hypoinverse. 'edit_slabgrid' was able to edit the grid
% so that it has the same spatial resolution as the old slab file used by
% John Armbruster and write the result to 'slab1.0.grid'. Now we want to
% go further. In the end, we create a grid of 0.002 deg * 0.002 deg.
% 
% 
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2022/03/07
% Last modified date:   2022/03/07
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load grid file
% close all
new = load('/home/data2/chaosong/Slab1.0_usgs/DepthGrid/cas_slab1.0_clip.xyz');
old = load('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/xyz.grid');

%This defines the maximum possbile range of the resulting grid without extrapolation
% lon = 121.92:0.002:125.76;
% nlon = length(lon);
% lat = 46.5:0.002:49.5;
% nlat = length(lat);
lon = 125.76:-0.002:121.92;
nlon = length(lon);
lat = 49.5:-0.002:46.5;
nlat = length(lat);

%generate a mesh grid based on the desired resolution
[longrd, latgrd] = meshgrid(lon,lat);

%cut the same range from the slab1.0
newcut = new(new(:,1)>=121.92 & new(:,1)<=125.76 & new(:,2)>=46.5 & new(:,2)<=49.5, :);
%generate a mesh grid based on cutted slab1.0
newcut = sortrows(newcut, [1 2],'descend');
nrow = length(unique(newcut(:,2))); % number of unique lat
ncol = length(unique(newcut(:,1))); % number of unique lon
nclongrd = reshape(newcut(:,1),nrow,ncol); 
nclatgrd = reshape(newcut(:,2),nrow,ncol); 
ncdepgrd = reshape(newcut(:,3),nrow,ncol); 

%2-D interpolation from old grid to new grid
depgrd = interp2(nclongrd,nclatgrd,ncdepgrd,longrd,latgrd,'spline');

%reshape for plotting and saving
lonmesh = reshape(longrd,[],1);
latmesh = reshape(latgrd,[],1);
depmesh = reshape(depgrd,[],1);

% %%% check the results
% figure
% subplot(3,1,1)
% scatter(new(:,1),new(:,2), 4, new(:,3), 'filled','o');   hold on  %, 'MarkerEdgeColor', 'w')
% colormap(jet)
% colorbar
% caxis([round(min(new(:,3))), round(max(new(:,3)))]);
% plot([125.76 125.76 121.92 121.92 125.76], [46.5 49.5 49.5 46.5 46.5],'r-');
% axis([121 126 46 50]);
% title('Slab 1.0, 0.02*0.02 deg');
% 
% subplot(3,1,2)
% scatter(lonmesh, latmesh, 4, depmesh, 'filled','o');   hold on  %, 'MarkerEdgeColor', 'w')
% colormap(jet)
% colorbar
% caxis([round(min(new(:,3))), round(max(new(:,3)))]);
% plot([125.76 125.76 121.92 121.92 125.76], [46.5 49.5 49.5 46.5 46.5],'r-');
% axis([121 126 46 50]);
% title('Slab 1.0, interpolated, 0.002*0.002 deg');
% 
% subplot(3,1,3)
% scatter(old(:,1),old(:,2), 4, old(:,3), 'filled','o');   hold on  %, 'MarkerEdgeColor', 'w')
% colormap(jet)
% colorbar
% caxis([round(min(new(:,3))), round(max(new(:,3)))]);
% plot([125.76 125.76 121.92 121.92 125.76], [46.5 49.5 49.5 46.5 46.5],'r-');
% axis([121 126 46 50]);
% title('Old slab, 0.01*0.01 deg');

%%% save file
ngrd = [lonmesh latmesh depmesh];
ngrd = sortrows(ngrd,[2 1],'descend');
fid = fopen('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/slab1.0_2e-3.grid','w+');
fprintf(fid,'%.4f %.4f %.4f \n',ngrd');
fclose(fid);


keyboard