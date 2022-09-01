function [irun,runamp,runenv] = Runningampenv(trace,mwlen,overlap,avgmethod)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [irun,runamp,runenv] = Runningampenv(trace,mwlen,overlap,avgmethod)
% This is a function to calculate the running absolute amplitude and envelope
% of a trace inside a moving rectangle window with a size of 'mwlen' samples. 
% Options are available for 'avgmethod' as 'median' or 'mean'. You can also 
% specify the overlap in the same unit as moving window length 'mwlen', 
% whose default is mwlen-1, meaning consecutive windows.
% 
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/03/23
% Last modified date:   2022/03/31
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('overlap',mwlen-1);
defval('avgmethod','median');

% length of signal
lsig = length(trace);

% index as in the original trace of the center of the moving window
irun = (mwlen/2+1: mwlen-overlap: lsig-mwlen/2)'; 

%envelope of the detrended full trace
[envup,~] = envelope(detrend(trace));

%absolute amplitude of the detrended full trace
absamp = abs(detrend(trace));

runamp = zeros(length(irun), 1);
runenv = zeros(length(irun), 1);
if isequal(avgmethod, 'median')
  for i = 1: length(irun)
    icnt = irun(i);
    ind = max(icnt-mwlen/2, 1): min(icnt+mwlen/2-1, lsig);
    runamp(i,1) = median(absamp(ind));
    runenv(i,1) = median(envup(ind));
  end
  
elseif isequal(avgmethod, 'mean')
  for i = 1: length(irun)
    icnt = irun(i);
    ind = max(icnt-mwlen/2, 1): min(icnt+mwlen/2-1, lsig);
    runamp(i,1) = mean(absamp(ind));
    runenv(i,1) = mean(envup(ind));
  end
  
end
