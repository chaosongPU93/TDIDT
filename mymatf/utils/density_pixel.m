function density = density_pixel(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% density1d = density_pixel(x,y)
%
% This function is to obtain the density (number of points;
% repetition times) of a 2-D data set (x,y). Different from
% 'density_matrix' which bins data based on a bin of 
% rectangle grid, this function bins data pixel wise. It 
% simply counts the times that the same (x,y) repeats. 
% Therefore, the input data set should have identical entries.
% The output 'density1d' has 3 columns, first 2 of which are 
% x, y locations while the 3rd columns is the counts at this
% location.
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/03/14
% Last modified date:   2022/03/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = reshape(x,[],1);  % force to be a column vector
y = reshape(y,[],1);
data = [x y]; 

%find the unique entries
dataunq = unique(data,'rows','stable');
density = zeros(size(dataunq,1), 3);
for i = 1: size(dataunq,1)
  ind = find(ismember(data ,dataunq(i,:), 'rows'));
  density(i,1:2) = dataunq(i,:);
  density(i,3) = length(ind);
end




