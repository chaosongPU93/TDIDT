function y=gauss(x,mu,sigma)
% z=GAUSS(x,mu,sigma)
%
% Retuns a gaussian defined on 'x' with a mean of 'mu', and a standard
% deviation of 'sigma'.
%
% Example:
%
% x=-10:0.1:10; 
% for index=1:4 ; plot(x,gauss(x,index)) ; hold on ; end ; grid
% lin2col(gca)
% Written by FJS, August 6th 1998
% Modified by Chao, 2022/04/05

y=exp(-(x-mu).^2/(2*sigma^2))/(sigma*sqrt(2*pi));