function [irunccn,runccn] = RunningCC(trace1, trace2, cclen)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [irunccn,runccn] = RunningCC(trace1, trace2, cclen)
% This is a function to calculate the running normalized cross-correlation 
% between 2 traces within a running (step of 1 samples) window (cclen) that
% are pre-aligned based on the max CC of the whole length.
% The length of 'runccn' (lrcc) should be shorter than the length of trace
% (lsig) by 'cclen'.
% For the moving window where either trace segment is 0, the return would be
% 'nan', i think that is fair, as in such cases, rcc should be undefined.
% 
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/01/04
% Last modified date:   2022/01/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% length of signal
lsig = length(trace1);

%auto dot products
tr1auto = trace1.*trace1;
tr2auto = trace2.*trace2;
csum1auto = cumsum(tr1auto);
csum2auto = cumsum(tr2auto);

%cross dot products
tr12crs = trace1.*trace2;
csum12crs = cumsum(tr12crs);

%numerator and denorminator
if cclen == lsig
  tr12num = sum(tr12crs);
  tr1den = sum(tr1auto);
  tr2den = sum(tr2auto);
  runccn = tr12num./realsqrt(tr1den.*tr2den);
  irunccn = round(cclen/2);
else
  tr12num = csum12crs(cclen+1:lsig)-csum12crs(1:lsig-cclen);
  tr1den = csum1auto(cclen+1:lsig)-csum1auto(1:lsig-cclen);
  tr2den = csum2auto(cclen+1:lsig)-csum2auto(1:lsig-cclen);
  runccn = tr12num./realsqrt(tr1den.*tr2den); % normalized, running CC, --> [-1,1]
  runccn(isnan(runccn)) = 0;   % force the undefined portion to be 0 in case of denominator being 0
  irunccn = (cclen/2+1: lsig-cclen/2)'; % index as in the original trace
end
