function medval = median_pixel(x,y,val,method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% medval = median_pixel(x,y,val)
%
% Imagine you have a set of data with 'val' that is distributed at locations
% x,y. Given that a unique location (x,y) might have several data points, you
% want to have a stable estimate (median/mean) of the value of each unique
% location.
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/10/05
% Last modified date:   2022/10/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('method','median');

x = reshape(x,[],1);  % force to be a column vector
y = reshape(y,[],1);
val = reshape(val,[],1);
loc = [x y ]; 

%find the unique entries
locunq = unique(loc,'rows','stable');
medval = zeros(size(locunq,1), 3);
for i = 1: size(locunq,1)
  ind = find(ismember(loc ,locunq(i,:), 'rows'));
  medval(i,1:2) = locunq(i,:);
  if isequal(method, 'median')
    medval(i,3) = median(val(ind));
  elseif isequal(method, 'mean')
    medval(i,3) = mean(val(ind));
  end
    
end

% keyboard