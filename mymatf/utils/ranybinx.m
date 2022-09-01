function [xcnt,ycnt,y1sig] = ranybinx(x,y,avgmethod,nbin,binw)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [xcnt,ycnt,y1sig] = ranybinx(x,y,avgmethod,nbin,binw)
%
% Imagine you have a dataset which has properties 'x' and 'y', between which
% you want to know if there is any relation. The simplest way is just to plot
% the scatter x and y. This function does one more step. It tries to bin the 
% property 'x' based on the number of bins 'nbin', OR, bin width 'binw'. For
% the data points whose property 'x' fall into that bin, obtain the median or 
% mean of property 'y', and its +-1\sigma range. 
% 
% --Note that you should specify either 'nbin' or 'binw', but not both of them,
%   the other one is computed accordingly.
% 
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/03/28
% Last modified date:   2022/03/28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('avgmethod','median');
defval('nbin',10);
defval('binw',[]);

%force to a column vector
x = reshape(x,[],1);  
y = reshape(y,[],1);

%if 'binw' is specified
if isempty(nbin) && ~isempty(binw)
  nbin = ceil((max(x)-min(x))/binw);
%if 'nbin' is specified
elseif ~isempty(nbin) && isempty(binw)
  binw = (max(x)-min(x))/nbin;
%otherwise return error 
else
  error("You can't specify both 'nbin' and 'binw' \n");
end

binEdges = min(x): binw: min(x)+nbin*binw;
xcnt = zeros(nbin,1);
ycnt = zeros(nbin,1);
y1sig = zeros(nbin,1);

for i = 1: length(binEdges)-1
  ybin = y(x>=binEdges(i) & x<binEdges(i+1)); % y fall into the bin
  xcnt(i) = (binEdges(i)+binEdges(i+1))/2;
  if isequal(avgmethod, 'median')
    ycnt(i) = median(ybin);
  elseif isequal(avgmethod, 'mean')
    ycnt(i) = mean(ybin);
  end
  y1sig(i) = std(ybin);  
end

% keyboard