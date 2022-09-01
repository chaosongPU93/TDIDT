%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Similar to psd_lfe_tremor2.m, this is also to compare the psd estimate of lfe
% stacks and tremor records, using proper methods determined from trial and
% testing.
%
% Difference from psd_lfe_tremor2.m:
%   1. To simulate the noise in lfe stacks, we should cut the each 20-s data 
%       shifted by 10s before each lfe timing, get the spectrum of them, then
%       get the median of all spectrum as the estimate for noise.
%   2. It is possible that for some sampled frequencies, noise amplitude is 
%       higher than signal, so that subtraction would lead to negative value.
%       When converted to dB, DO NOT plot them.
%   3. To simualte the nosie for tremor, we should choose tremor free days, 
%       several days, several periods the same as migrations. Also, we need to
%       get the spectrum of each 20-s window, then obtain the median/mean, or 
%       say, 90 percentile as the error estimate. Finally, do NOT plot the 
%       negative value when subtract the nosie from the signal
%   
% 
%   NO implementation here, but it is not difficult to do it. The main reason
%   not writing here is that the result in this version should not be too 
%   distinct than version 2.
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/12/20
% Last modified date:   2019/12/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

