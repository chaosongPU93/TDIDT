function [pn, pd, errflag, chi2, Q] = ratregr(x, y, sigy, p, q, basis, monflag)

% RATREGR Generalized rational function regression.
%    [PN,PD] = RATREGR(X,Y,SIGY,P,Q,BASIS) calculates numerator, PN, and
%    denominator, PD, regression coefficients for N data points (X,Y) where
%    SIGY are the uncertainties associated with each measurement Y.  If
%    SIGY are unknown (omitted), they are assumed to be unity.
%
%    P specifies the order of the of the numerator basis functions (powers
%    of x).  It defaults to 1 (linear regression).  Q specifies the order
%    of the denominator basis functions.  It defaults to 0.
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
%    [PN,PD,ERRFLAG] = RATREGR() also returns an error flag associated
%    with the regression.  Possible values of ERRFLAG and the corresponding
%    error conditions are:
%       -1 : RATREGR terminated before a solution could be calculated
%        0 : RATREGR converged to a valid solution
%        1 : The kernel of the A matrix has no values (the solution may
%            not be sufficient)
%        2 : The kernel of the A matrix has a dimension larger than 1 (the
%            solution may not be sufficient)
%        3 : The numerator of the solution is zero
%        4 : The denominator of the solution is zero
%
%    [PN,PD,ERRFLAG,CHI2,Q] = RATREGR() also returns the chi-square merit
%    CHI2 and probability Q that a value of chi-square as poor as CHI2
%    should occur by chance.
%
%    [PN,PD] = RATREGR(...,'mon') makes the denominator monic.
%
%    See also SVDREGR, POLYFIT.

% Joe Henning - Jan 2014

x = x(:);
y = y(:);
r = 0;

if (nargin < 3)
   sigy = ones(size(y));
   r = 1;
end

if (nargin < 4)
   fprintf('??? Bad p, q input to ratregr ==> p and q should be specified\n');
   p = 1;
end

if (nargin < 5)
   q = 0;
end

if (nargin < 6)
   basis = 'standard';
end

basis_options = {'standard','cheby1','cheby2','hermite','laguerre','legendre'};
indx = find(strncmpi(basis,basis_options,length(basis)));
if isempty(indx)
   fprintf('   Error ==> basis must be a string as defined in the help section\n');
   pn = [];
   pd = [];
   errflag = -1;
   return
end
basis = basis_options{indx};

if (p < 0)
   fprintf('??? Bad p input to ratregr ==> p>=0\n');
   pn = [];
   pd = [];
   errflag = -1;
   return
end

if (q < 0)
   fprintf('??? Bad p input to ratregr ==> q>=0\n');
   pn = [];
   pd = [];
   errflag = -1;
   return
end

if (isempty(sigy))
   sigy = ones(size(y));
   r = 1;
end

if (length(x) ~= length(y))
   fprintf('   Error ==> length(x) must be equal to length(y)\n');
   pn = [];
   pd = [];
   errflag = -1;
   return
end

if (length(y) ~= length(sigy))
   fprintf('   Error ==> length(y) must be equal to length(sigy)\n');
   pn = [];
   pd = [];
   errflag = -1;
   return
end

n = length(x);

A = [];
for i = 0:p
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
for i = 0:q
   for j = 1:n
      switch basis
         case 'standard'
            A(j,p+2+i) = -y(j)*polyb(x(j),i)/sigy(j);
         case 'cheby1'
            A(j,p+2+i) = -y(j)*cheby1b(x(j),i)/sigy(j);
         case 'cheby2'
            A(j,p+2+i) = -y(j)*cheby2b(x(j),i)/sigy(j);
         case 'hermite'
            A(j,p+2+i) = -y(j)*hermb(x(j),i)/sigy(j);
         case 'laguerre'
            A(j,p+2+i) = -y(j)*lagb(x(j),i)/sigy(j);
         case 'legendre'
            A(j,p+2+i) = -y(j)*legendb(x(j),i)/sigy(j);
      end
   end
end

% calculate the kernel of A to find the least-squares solution
%u = null(A)
[u, errflag] = kernel(A);

pn = u(1:p+1);
pd = u(p+2:length(u));

if( length(pn) == 1 && pn(1) == 0)
   pn = [];
   pd = [];
   errflag = 3;
   return
elseif( length(pd) == 1 && pd(1) == 0)
   pn = [];
   pd = [];
   errflag = 4;
   return
end

if (nargin == 7)
   % make the denominator monic
   tol = eps*max([max(abs(u)) 1]);
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

yhat = [];
for j = 1:n
   yhat(j) = 0;
   denom = 0;
   for i = 0:q
      switch basis
         case 'standard'
            denom = denom + pd(i+1)*polyb(x(j),i);
         case 'cheby1'
            denom = denom + pd(i+1)*cheby1b(x(j),i);
         case 'cheby2'
            denom = denom + pd(i+1)*cheby2b(x(j),i);
         case 'hermite'
            denom = denom + pd(i+1)*hermb(x(j),i);
         case 'laguerre'
            denom = denom + pd(i+1)*lagb(x(j),i);
         case 'legendre'
            denom = denom + pd(i+1)*legendb(x(j),i);
      end
   end
   for i = 0:p
      switch basis
         case 'standard'
            yhat(j) = yhat(j) + pn(i+1)*polyb(x(j),i)/denom;
         case 'cheby1'
            yhat(j) = yhat(j) + pn(i+1)*cheby1b(x(j),i)/denom;
         case 'cheby2'
            yhat(j) = yhat(j) + pn(i+1)*cheby2b(x(j),i)/denom;
         case 'hermite'
            yhat(j) = yhat(j) + pn(i+1)*hermb(x(j),i)/denom;
         case 'laguerre'
            yhat(j) = yhat(j) + pn(i+1)*lagb(x(j),i)/denom;
         case 'legendre'
            yhat(j) = yhat(j) + pn(i+1)*legendb(x(j),i)/denom;
      end
   end
end

chi2 = sum((yhat(:)-y).^2);

Q = gammainc((n-2)/2,chi2/2);

if (r == 1)
   Q = NaN;
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


function [Y, err] = kernel(X)
% Calculate the kernel or null space of matrix X
[m,n] = size(X);
minmn = min(m,n);
maxmn = max(m,n);
[U,S,V] = svd(X);
s = [];
for i = 1:minmn
   s = [s; S(i,i)];
end
%tol = eps*max([ones(size(s)) abs(s)],[],2);
tol = eps*max([maxmn max(abs(s))]);
r = sum(s > tol);
Y = V(:,r+1:n);
% check the dimension of the kernel and report if a
% "valid" solution is not found
err = 0;
if (size(Y,2) < 1)
   err = 1;
elseif (size(Y,2) > 1)
   err = 2;
end
% approximate the kernel by the nullest vector
Y = V(:,n);
