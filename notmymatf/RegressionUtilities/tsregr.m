function [b0, b1] = tsregr(x,y)

% TSREGR Theil-Sen estimator.
%    [B0, B1] = TSREGR(X,Y) calculates straight line coefficients
%       Y = B0 + B1*X
%    for N data points {X,Y} using the Theil-Sen estimator.

% Joe Henning - Fall 2011

if (length(x) ~= length(y))
   fprintf('   Error ==> length(x) must be equal to length(y)');
   b0 = NaN;
   b1 = NaN;
   return
end

n = length(x);

m = [];
for i = 1:n
   for j = i:n
      if (i ~= j && x(i) ~= x(j))
         slope = (y(j)-y(i))/(x(j)-x(i));
         m = [m slope];
      end
   end
end

xm = median(x);
ym = median(y);

b1 = median(m);
%b0 = median(y) - b1*median(x);
b0 = median(y - b1*x);
