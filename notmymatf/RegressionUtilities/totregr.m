function [p, fval, exitflag, output] = totregr(x, y, sigx, sigy, c)

% TOTREGR Fit TLS polynomial to data
%    P = TOTREGR(X,Y,SIGX,SIGY,C) calculates the Total Least Squares
%    polynomial P(X) of N data points (X,Y) using the Delta Method.  SIGX
%    are the uncertainties associated with measurements X, and SIGY are the
%    uncertainties associated with measurements Y.  If SIGY are unknown
%    (omitted), they are assumed to be unity.  The Nelder-Mead simplex
%    (direct search) method of FMINSEARCH is used to calculate the
%    regression coefficients.
%
%    C specifies the order of the basis functions (powers of x).  It
%    defaults to 1 (linear regression).
%
%    [P,FVAL,EXITFLAG,OUTPUT] = TOTREGR() also returns FVAL, the value of
%    the total least squares objective function; EXITFLAG, the exit
%    condition of FMINSEARCH,; and OUTPUT, a structure with FMINSEARCH
%    output parameters.
%
%    See also fminsearch.

% Joe Henning - Fall 2013

if (nargin < 4)
   sigy = ones(size(y));
   c = 1;
elseif (nargin < 5)
   c = 1;
end

if (isempty(sigx))
   sigx = zeros(size(x));
end

if (isempty(sigy))
   sigy = ones(size(y));
end

if (length(x) ~= length(y))
   fprintf('   Error ==> length(x) must be equal to length(y)\n');
   p = [];
   return
end

if (length(x) ~= length(sigx))
   fprintf('   Error ==> length(x) must be equal to length(sigx)\n');
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

% Seed the algorithm with the least-squares solution
p0 = svdregr(x, y, sigy, c);

global X Y SIGX SIGY C;
X = x;
Y = y;
SIGX = sigx;
SIGY = sigy;
C = c;

s = sumtot(p0);

[p1,fval,exitflag,output] = fminsearch(@sumtot,p0);

p = p1;


function s = sumtot(p)
global X Y SIGX SIGY C;
% calculate first derivative of p
pd = 0;
for i = 2:length(p)
   pd = pd + (i-1)*p(i)*X.^(i-2);
end
% use the Delta method to approximate the transformed variance
u = pd.*pd.*SIGX.*SIGX;
% calculate total least squares of the residuals
s = 0;
n = length(X);
for i = 1:n
   fx = 0;
   for j = 0:C
      fx = fx + p(j+1)*X(i)^(j);
   end
   s = s + (Y(i) - fx)*(Y(i) - fx)/(SIGY(i)*SIGY(i) + u(i));
end
