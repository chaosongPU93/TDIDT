function [irunccn,runccn] = RunningCC(trace1, trace2, cclen, method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [irunccn,runccn] = RunningCC(trace1, trace2, cclen, method)
% This is a function to calculate the running normalized cross-correlation 
% between 2 traces within a running (step of 1 samples) window (cclen) that
% are pre-aligned based on the max CC of the whole length.
% The length of 'runccn' (lrcc) should be shorter than the length of trace
% (lsig) by 'cclen'.
% For the moving window where either trace segment is 0, the return would be
% 'nan', i think that is fair, as in such cases, rcc should be undefined.
% Other than taking advantage of 'cumsum', you can also choose 'xcorr' in 
% case of extremely small values like 0/0.
% 
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/01/04
% Last modified date:   2022/01/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('method','cumsum');

% length of signal
lsig = length(trace1);
trace1 = detrend(trace1); %just in case
trace2 = detrend(trace2);


if strcmp(method,'cumsum')  
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
    runccn(isnan(runccn) | isinf(runccn)) = 0;   % force the undefined portion to be 0 in case of denominator being 0
%     runccn(abs(tr12num)<1e-7 | abs(tr1den)<1e-7 | abs(tr2den)<1e-7) = 0;
    irunccn = (cclen/2+1: lsig-cclen/2)'; % index as in the original trace
  end
  
  % figure
  % subplot(411);
  % plot(tr12num);
  %
  % subplot(412);
  % plot(tr1den);
  %
  % subplot(413);
  % plot(tr2den);
  %
  % subplot(414);
  % plot(runccn);
  %
  % keyboard
  
elseif strcmp(method,'xcorr')
  % index as in the original trace of the center of the moving window
  irunccn = (cclen/2+1: lsig-cclen/2)';
  runccn = zeros(length(irunccn), 1);
  for i = 1: length(irunccn)
    icnt = irunccn(i);
    ind = max(icnt-cclen/2, 1): min(icnt+cclen/2-1, lsig);
    tr1 = detrend(trace1(ind)); %remove mean
    tr2 = detrend(trace2(ind)); %remove mean
    runccn(i) = xcorr(tr1,tr2,0,'normalized');  %0-lag maximum cc based on current alignment
  end
  runccn(isnan(runccn) | isinf(runccn)) = 0;   % force the undefined portion to be 0 in case of denominator being 0  

end









