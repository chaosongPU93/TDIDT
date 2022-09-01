function [irun,runcorr,runpval] = Runningcorr(x,y,mwlen,overlap,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [irun,runcorr,runpval] = Runningcorr(x,y,mwlen,type)
% This is a function to calculate the running statistical correlation between 
% quantities x and y. Several 'type' are available, like 'Pearson', 'Spearman'
% and 'Kendall'. Will call the matlab built-in function 'corr' to actually 
% compute them. You can also specify the overlap in the same unit as moving
% window length 'mwlen', whose default is mwlen-1, meaning consecutive windows.
%
% 'runcorr', running correlation coefficient, depending on 'type'.
% 'runpval', running p values, for testing the hypothesis of no correlation
% against the alternative that there is a non-zero correlation. 
%    
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/03/31
% Last modified date:   2022/03/31
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('overlap',mwlen-1);
defval('type','Pearson');

% length of quantity x, should be the same as y
x = reshape(x,[],1);
y = reshape(y,[],1);
lenx = length(x);
leny = length(y);
if lenx ~= leny
  error('The two quantities should have the same size');
end

% index as in the original trace of the center of the moving window
irun = (mwlen/2+1: mwlen-overlap: lenx-mwlen/2)'; 

runcorr = zeros(length(irun), 1);
runpval = zeros(length(irun), 1);
for i = 1: length(irun)
  icnt = irun(i);
  ind = max(icnt-mwlen/2, 1): min(icnt+mwlen/2-1, lenx);
  %the correlation between x and y based on the requested type
  [runcorr(i,1),runpval(i,1)] = corr(x(ind),y(ind),'Type',type);
end

