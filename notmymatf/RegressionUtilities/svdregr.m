function [p, varp, covp, chi2, Q] = svdregr(x, y, sigy, c, basis)

% SVDREGR Generalized linear least squares regression using Singular Value Decomposition.
%    SVDREGR(X,Y,SIGY,C,BASIS) calculates regression coefficients for N
%    data points {X,Y} where SIGY are the uncertainties associated with
%    each measurement Y.  If SIGY are unknown (omitted), they are assumed
%    to be unity.
%
%    C specifies the order of the basis functions (powers of x).  It
%    defaults to 1 (linear regression).
%
%    BASIS specifies the type of basis functions used in the regression.
%    It can be:
%        'standard' - powers of x basis
%        'cheby1'   - Chebyshev polynomials of the first kind basis
%        'cheby2'   - Chebyshev polynomials of the second kind basis
%        'hermite'  - Hermite polynomials basis
%        'laguerre' - Laguerre polynomials basis
%        'legendre' - Legendre polynomials basis
%
%    [P,VARP,COVP] = SVDREGR() also returns VARP, the variances in the
%    estimates of P, and COVP, the covariance of P.
%
%    [P,VARP,COVP,CHI2,Q] = SVDREGR() also returns the chi-square merit
%    CHI2 and probability Q that a value of chi-square as poor as CHI2
%    should occur by chance.
%
%    See also POLYFIT.

% Joe Henning - Fall 2011

q = 0;
if (nargin < 3)
   sigy = ones(size(y));
   q = 1;
   c = 1;
   basis = 'standard';
elseif (nargin < 4)
   c = 1;
   basis = 'standard';
elseif (nargin < 5)
   basis = 'standard';
end

basis_options = {'standard','cheby1','cheby2','hermite','laguerre','legendre'};
indx = find(strncmpi(basis,basis_options,length(basis)));
if isempty(indx)
   fprintf('   Error ==> basis must be a string as defined in the help section\n');
   p = [];
   varp = [];
   covp = [];
   chi2 = NaN;
   Q = NaN;
   return
end
basis = basis_options{indx};

if (isempty(sigy))
   sigy = ones(size(y));
   q = 1;
end

if (length(x) ~= length(y))
   fprintf('   Error ==> length(x) must be equal to length(y)\n');
   p = [];
   varp = [];
   covp = [];
   chi2 = NaN;
   Q = NaN;
   return
end

if (length(y) ~= length(sigy))
   fprintf('   Error ==> length(y) must be equal to length(sigy)\n');
   p = [];
   varp = [];
   covp = [];
   chi2 = NaN;
   Q = NaN;
   return
end

n = length(x);

if (c >= n)
   fprintf('   Warning: polynomial degree >= number of data points; reducing\n');
   c = n-1;
end

A = [];
for i = 0:c
   for j = 1:n
      switch basis
         case 'standard'
            A(j,i+1) = polyb(x(j),i)/sigy(j);
         case 'cheby1'
            A(j,i+1) = cheby1b(x(j),i)/sigy(j);
         case 'cheby2'
            A(j,i+1) = cheby2b(x(j),i)/sigy(j);
         case 'hermite'
            A(j,i+1) = hermb(x(j),i)/sigy(j);
         case 'laguerre'
            A(j,i+1) = lagb(x(j),i)/sigy(j);
         case 'legendre'
            A(j,i+1) = legendb(x(j),i)/sigy(j);
      end
   end
end

[N,M] = size(A);

b = [];
for i = 1:n
   b(i,1) = y(i)/sigy(i);
end

[u,s,v] = svd(A);

tol = eps;

p = zeros(M,1);
for i = 1:M
   % If s(i,i) is zero, its reciprocal is set to zero, not infinity
   if (abs(s(i,i)) > tol*max([1 abs(s(i,i))]))
      p = p + sum(u(:,i).*b)/s(i,i)*v(:,i);
   end
end

varp = zeros(M,1);
for i = 1:M
   for j = 1:M
      if (abs(s(j,j)) > tol*max([1 abs(s(j,j))]))
         varp(i) = varp(i) + (v(i,j)/s(j,j))^2;
      end
   end
end

chi2 = sum((A*p - b).^2);

Q = gammainc((n-2)/2,chi2/2);

if (q == 1)
   % A good fit was assumed, so there is no independent fit
   % probability Q
   varp = varp*chi2/(n-2);
   Q = NaN;
end

covp = zeros(M,M);
for i = 1:M
   for j = 1:M
      if (i ~= j)
         for k = 1:M
            if (abs(s(k,k)) > tol*max([1 abs(s(k,k))]))
               covp(i,j) = covp(i,j) + v(i,k)*v(j,k)/s(k,k)/s(k,k);
            end
         end
      else
         covp(i,j) = varp(i);
      end
   end
end


function yi = polyb(xi, c)
% powers of x basis
yi = xi^(c);


function yi = cheby1b(xi, c)
% Chebyshev polynomials of the first kind basis
t(1) = 1;
t(2) = xi;
if (c == 0)
   yi = t(1);
   return
elseif (c == 1)
   yi = t(2);
   return
end
for k = 3:c+1
   t(k) = 2*xi*t(k-1) - t(k-2);
end
yi = t(c+1);


function yi = cheby2b(xi, c)
% Chebyshev polynomials of the second kind basis
u(1) = 1;
u(2) = 2*xi;
if (c == 0)
   yi = u(1);
   return
elseif (c == 1)
   yi = u(2);
   return
end
for k = 3:c+1
   u(k) = 2*xi*u(k-1) - u(k-2);
end
yi = u(c+1);


function yi = hermb(xi, c)
% Hermite polynomials basis
h(1) = 1;
h(2) = xi;
if (c == 0)
   yi = h(1);
   return
elseif (c == 1)
   yi = h(2);
   return
end
for k = 3:c+1
   h(k) = xi*h(k-1) - (k-2)*h(k-2);
end
yi = h(c+1);


function yi = lagb(xi, c)
% Laguerre polynomials basis
l(1) = 1;
l(2) = 1-xi;
if (c == 0)
   yi = l(1);
   return
elseif (c == 1)
   yi = l(2);
   return
end
for k = 3:c+1
   l(k) = ((2*(k-2) + 1 - xi)*l(k-1) - (k-2)*l(k-2))/(k-1);
end
yi = l(c+1);


function yi = legendb(xi, c)
% Legendre polynomials basis
l(1) = 1;
l(2) = xi;
if (c == 0)
   yi = l(1);
   return
elseif (c == 1)
   yi = l(2);
   return
end
for k = 3:c+1
   l(k) = ((2*(k-2) + 1)*xi*l(k-1) - (k-2)*l(k-2))/(k-1);
end
yi = l(c+1);
