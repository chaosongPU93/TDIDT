function [s_s] = shift_freq(s, N, dt, tshift)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[s_s] = shift_freq(s, N, dt, tshift)
% Shift the time-domain signal by a given time shift in the frequency domain 
%
% INPUT:
%    s: Signal given in time domain.
%    N: number of FFT, usually power of 2
%    dt: Sampling interval [s].
%    tshift: time shift [s], postive means shift to right/would arrive later
%
% OUTPU:T
%   s_s: shifted signal in time domain
%
% NOTE:
% --the code is almost
%   translated into Matlab from a Python code by Peter Makus (makus@gfz-potsdam.de)
%   https://github.com/PeterMakus/PyGLImER/blob/1b17da8409388702ac0afc21a651000c304d16e4/src/
%     pyglimer/utils/signalproc.py#L215
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/11/16
% Last modified date:   2021/11/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = fft(s, N);  % signal in freq domain
k = round(tshift/dt);   % time shift in samples
p = (2*pi*(1:N)*k/N)';  % phase shift, make it as colomun vector too
% p = 2*pi*k/N;
S = S.*(cos(p)-1i.*sin(p));  % dot product between S and phase shift

s_s = real(ifft(S, N)) / cos(2*pi*k/N);  % go back to time-domain and correct for the scaling 
