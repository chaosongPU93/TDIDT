function [s_f] = filter_freq(s, F, dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[s_f] = filter_freq(s, F, dt)
% filter the time-domain signal s with a filter with response F in the 
% frequency domain via multiplication
%
% INPUT:
%    s: Signal given in time domain.
%    F: Filter's amplitude response (i.e. frequency domain)
%    dt: Sampling interval [s].
%
% OUTPU:T
%   s_f: Filtered signal in time domain
%
% NOTE:
% --the code is almost
%   translated into Matlab from a Python code by Peter Makus (makus@gfz-potsdam.de)
%   https://github.com/PeterMakus/PyGLImER/blob/1b17da8409388702ac0afc21a651000c304d16e4/src/
%     pyglimer/utils/signalproc.py#L121
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/11/16
% Last modified date:   2021/11/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nfft = length(F);

S = fft(s, nfft);
S_f = S.* F *dt;
s_f = real(ifft(S_f, nfft));

