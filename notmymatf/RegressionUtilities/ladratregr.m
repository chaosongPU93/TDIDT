function [pn, pd, fval, exitflag, output] = ladratregr(x, y, sigy, p, q, monflag)

% LADRATREGR Fit LAD rational polynomial to data.
%    [PN,PD] = LADRATREGR(X,Y,SIGY,P,Q) finds the coefficients for the
%    numerator PN(X) of degre P and denominator PD(X) of degree Q of a
%    rational polynomial that fits the data Y best in a least absolute
%    deviation sense.  SIGY are the uncertainties associated with
%    measurement Y.  If SIGY are unknown (omitted), they are assumed to be
%    unity.  The Nelder-Mead simplex (direct search) method of FMINSEARCH
%    is used to calculate the regression coefficients.
%
%    P specifies the order of the of the numerator basis functions (powers
%    of x).  It defaults to 1 (linear regression).  Q specifies the order
%    of the denominator basis functions.  It defaults to 0.
%
%    [PN,PD,FVAL,EXITFLAG,OUTPUT] = LADRATREGR() also returns FVAL, the
%    value of the absolute deviation objective function; EXITFLAG, the exit
%    condition of FMINSEARCH,; and OUTPUT, a structure with FMINSEARCH
%    output parameters.
%
%    [PN,PD] = RATREGR(...,'mon') makes the denominator monic.
%
%    See also RATREGR, fminsearch.

% Joe Henning - Jan 2014

if (nargin < 3)
   sigy = ones(size(y));
   p = 1;    
   q = 0;
elseif (nargin < 4)
   p = 1;
   q = 0;
elseif (nargin < 5)
   q = 0;
end

if (isempty(sigy))
   sigy = ones(size(y));
end

if (length(x) ~= length(y))
   fprintf('   Error ==> length(x) must be equal to length(y)\n');
   pn = [];
   pd = [];
   return
end

if (length(y) ~= length(sigy))
   fprintf('   Error ==> length(y) must be equal to length(sigy)\n');
   pn = [];
   pd = [];
   return
end

n = length(x);

% seed the algorithm with the least-squares solution
[pn0, pd0] = ratregr(x,y,sigy,p,q);
p0 = [pn0; pd0];

global X Y SIGY P Q;
X = x;
Y = y;
SIGY = sigy;
P = p;
Q = q;

s = sumdev(p0);

[p1,fval,exitflag,output] = fminsearch(@sumdev,p0);

pn = p1(1:p+1);
pd = p1(p+2:end);

if nargin == 6
   % make the denominator monic
   tol = eps*max([max(abs(p1)) 1]);
   for i = q+1:-1:1
      d = pd(i);
      if (abs(d) > tol)
         break
      else
         pd = pd(1:i-1);
      end
   end

   pn = pn/d;
   pd = pd/d;
end


function s = sumdev(p)
% calculate the absolute values of the residuals
global X Y SIGY P Q;
pn = p(1:P+1);
pd = p(P+2:end);
s = 0;
n = length(X);
for i = 1:n
   fx = 0;
   denom = 0;
   for j = 0:Q
      denom = denom + pd(j+1)*X(i)^(j);
   end
   for j = 0:P
      fx = fx + pn(j+1)*X(i)^(j)/denom;
   end
   s = s + abs(Y(i) - fx)/SIGY(i);
end
