function [b0, b1, xest, yest] = dregr(x,y,delta)

% DREGR Deming regression.
%    [B0, B1] = DREGR(X,Y,DELTA) calculates straight line coefficients
%       Y = B0 + B1*X
%    for N data points {X,Y} where the measurement errors in X and Y are
%    independent and the ratio of their variances is
%       DELTA = VAR{error in Y}/VAR{error in X}
%    When the errors in X and Y are equally likely, DELTA=1 (default).
%
%    [B0, B1, XEST, YEST] = DREGR() provides estimates of the "true" values
%    of X and Y.

% Joe Henning - Fall 2011

if (nargin < 3)
   delta = 1;
end

if (length(x) ~= length(y))
   fprintf('   Error ==> length(x) must be equal to length(y)');
   b0 = NaN;
   b1 = NaN;
   return
end

n = length(x);
xmean = 0;
ymean = 0;
for i = 1:n
   xmean = xmean + x(i)/n;
   ymean = ymean + y(i)/n;
end
sxx = 0;
sxy = 0;
syy = 0;
for i = 1:n
   sxx = sxx + (x(i)-xmean)*(x(i)-xmean)/(n-1);
   sxy = sxy + (x(i)-xmean)*(y(i)-ymean)/(n-1);
   syy = syy + (y(i)-ymean)*(y(i)-ymean)/(n-1);
end

b1 = (syy - delta*sxx + sqrt((syy - delta*sxx)*(syy - delta*sxx) + 4*delta*sxy*sxy))/(2*sxy);
b0 = ymean - b1*xmean;
xest = [];
for i = 1:n
   xest(i) = x(i) + b1/(b1*b1 + delta)*(y(i) - b0 - b1*x(1));
   yest(i) = b0 + b1*xest(i);
end
