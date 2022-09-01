function [G] = gauss_zerophase(N, dt, width)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [G] = gauss_zerophase(N, dt, width)
%
% Create a zero-phase Gaussian function (i.e., a low-pass filter). In
% particular meant to be convolved with the impulse
% response that is the output of an iterative deconvolution. Note this
% function is meant to be used in frequency domain, see also 'gauss_time'.
% Note also that this should be the same as 
% '/notmymatf/deconvolutions/GLImER/gaussian.m'.
%
%
% INPUT:
%   N: length of array
%   dt: sampling interval
%   width: width parameter of the Gaussian function
%
% OUTPUT:
%   G: Gaussian window
%
%
% NOTE:
% --the code is almost
%   translated into Matlab from a Python code by Peter Makus (makus@gfz-potsdam.de)
%   https://github.com/PeterMakus/PyGLImER/blob/1b17da8409388702ac0afc21a651000c304d16e4/
%    src/pyglimer/utils/signalproc.py#L121
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/11/14
% Last modified date:   2021/11/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df = 1/(N*dt);  % frequency step
f = df*(0:1:round(0.5*N));  % frequency array
f = reshape(f,[],1);
w = 2*pi.*f;  % angular frequency omega

G = zeros(N,1);  % predefine
temp = exp(-w.^2./(4*width^2))/dt;
G(1: round(0.5*N)+1) = temp;
temp_ud = flipud(temp);
G(round(0.5*N)+2: N) = temp_ud(length(temp_ud)-length(round(0.5*N)+2: N)+1: length(temp_ud));
