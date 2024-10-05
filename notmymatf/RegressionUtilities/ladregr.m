function [p, fval, exitflag, output] = ladregr(x, y, sigy, c)

% LADREGR Fit LAD polynomial to data.
%    P = LADREGR(X,Y,SIGY,C) finds the coefficients of a polynomial
%    P(X) of degree C that fits the data Y best in a least absolute
%    deviation sense.  SIGY are the uncertainties associated with
%    measurement Y.  If SIGY are unknown (omitted), they are assumed
%    to be unity.  The Nelder-Mead simplex (direct search) method of
%    FMINSEARCH is used to calculate the regression coefficients.
%
%    C specifies the order of the basis functions (powers of x).  It
%    defaults to 1 (linear regression).
%
%    [P,FVAL,EXITFLAG,OUTPUT] = LADREGR() also returns FVAL, the value
%    of the absolute deviation objective function; EXITFLAG, the exit
%    condition of FMINSEARCH,; and OUTPUT, a structure with
%    FMINSEARCH output parameters.
%
%    See also fminsearch.

% Joe Henning - Fall 2012

if (nargin < 3)
   sigy = ones(size(y));
   c = 1;    
elseif (nargin < 4)
   c = 1;
end

if (isempty(sigy))
   sigy = ones(size(y));
end

if (length(x) ~= length(y))
   fprintf('   Error ==> length(x) must be equal to length(y)\n');
   p = [];
   return
end

if (length(y) ~= length(sigy))
   fprintf('   Error ==> length(y) must be equal to length(sigy)\n');
   p = [];
   return
end

n = length(x);

if (c >= n)
   fprintf('   Warning: polynomial degree >= number of data points; reducting\n');
   c = n-1;
end

% seed the algorithm with the least-squares solution
p0 = svdregr(x,y,sigy,c);

global X Y SIGY C;
X = x;
Y = y;
SIGY = sigy;
C = c;

s = sumdev(p0);

[p1,fval,exitflag,output] = fminsearch(@sumdev,p0);

p = p1;


function s = sumdev(p)
% calculate the absolute values of the residuals
global X Y SIGY C;
s = 0;
n = length(X);
for i = 1:n
   fx = 0;
   for j = 0:C
      fx = fx + p(j+1)*X(i)^(j);
   end
   s = s + abs(Y(i) - fx)/SIGY(i);
end
