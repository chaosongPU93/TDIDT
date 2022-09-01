% function gen_time2loc_rectgrid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to generate an even, square grid of times in samples at
% a specific sampling rate, so that it can be used in hypoinverse to tranform
% the offset times between station 12 and 13, into relative locations. 
% 
% --Generate several grids at different sps, like 40, 80, 100, can be used for
%   further interpolation if higher resolution is needed
% --This is also useful to test the spatial resolution of hypoinverse. If the
%   sampling in time is too fine so that the space transformation is becoming
%   non-linear and unstable, we'd better use interpolation from a lower-
%   resolution grid.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/02/10
% Last modified date:   2022/02/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clc
clear
% close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);

% set path
workpath = getenv('MHOME');
lzbpath = strcat(workpath, '/Seisbasics/hypoinverse/forsummary');

% fams used
nfampool = [
            '002';
           ];
nfam = size(nfampool, 1);
         
% relative arrival time to main station for each fam at 40 sps.
contpgcoff = [
              86 20;
              ];         
         
% set params
sps = 20;

spsscale = sps/40;
% offmax = (4*1.5)*spsscale;
offmax = 12;

%% generation of rectangle grid
arr = [];
for ifam = 1: nfam
  i = 1;
  tmp = zeros((2*offmax+1)^2, 4);
  for off12 = -offmax: 1: offmax
    for off13 = -offmax: 1: offmax
      tmp(i,1:2) = contpgcoff(ifam,:);
      tmp(i,3:4) = [off12 off13];
      i = i+1;
    end
  end
  arr = [arr; tmp];        

end

figure
ax = gca;
scatter(ax,tmp(:,3),tmp(:,4),10,'r','filled');
box on; grid on
axis equal
axis(ax,[-offmax offmax -offmax offmax])

fid = fopen(fullfile(lzbpath, strcat('offset_002_rectgrid_',num2str(offmax),'_',num2str(sps))),'w+');
fprintf(fid,'%d %d %f %f \n',arr');
fclose(fid);


% %% generate the same grid as 'p085p020'
% arr = [];
% for ifam = 1: nfam
%   i = 1;
%   tmp = zeros(25^2, 4);
%   for off12 = -23:2:25
%     for off13 = -24:2:24
%       tmp(i,1:2) = contpgcoff(ifam,:);
%       tmp(i,3:4) = [off12 off13];
%       i = i+1;
%     end
%   end
%   arr = [arr; tmp];        
% 
% end
% 
% figure
% ax = gca;
% scatter(ax,tmp(:,3),tmp(:,4),10,'r','filled');
% box on; grid on
% axis equal
% axis(ax,[-offmax offmax -offmax offmax])
% 
% fid = fopen(fullfile(lzbpath, strcat('offset_002_rectgrid_','armb','_',num2str(sps))),'w+');
% fprintf(fid,'%d %d %f %f \n',arr');
% fclose(fid);


  
  
  
  