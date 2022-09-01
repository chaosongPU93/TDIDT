function [cc] = corr_freq(u, v, nf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cross-correlation in frequency domain, note that this is the plain CC without
% normalization; the order of u, v does not matter; one-side compared to built-in
% 'xcorr' with the option 'none', 
% ie. corrf(u, v, nf) = the right side of  xcorr(u, v, nf, 'none')
%
% INPUT:
%   u : flipped vector.
%   v : Vector that the correlation is measured to. here u and v have the same length
%   nf : Array length in frequency domain (use next power of 2).
%
% OUTPU:T
%   cc: cross-correlation
%
% NOTE:
% --the code is almost
%   translated into Matlab from a Python code by Peter Makus (makus@gfz-potsdam.de)
%   https://github.com/PeterMakus/PyGLImER/blob/1b17da8409388702ac0afc21a651000c304d16e4/src/
%     pyglimer/utils/signalproc.py#L91
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/11/16
% Last modified date:   2021/11/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = conj(fft(v,nf));
U = fft(u,nf);
CC = U.*V;
cc = real(ifft(CC, nf));