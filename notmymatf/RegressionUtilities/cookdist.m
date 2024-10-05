function [d] = cookdist(x, y, sigy, c, basis, plotflag)

% COOKDIST Cook's distance
%    COOKDIST(X,Y,SIGY,C,BASIS) calculates Cook's distance for N data
%    points {X,Y} where SIGY are the uncertainties associated with each
%    measurement Y.  If SIGY are unknown (omitted), they are assumed to be
%    unity.
%
%    Cook's distance measures the effect of deleting a given observation
%    from a least squares regression and can be used to estimate the
%    influence of a data point.  Highly influential data points with large
%    distances can distort the regression.  A typical operational guideline
%    is to declare influential points whose distance > 1.  Another
%    guideline thresholds distances at 4/N or 4/(N-k-1) where N is the
%    number of observations an k is the number of explanatory variables.
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
%    [] = COOKDIST(..,'plot') plots a figure which depicts Cook's distances
%    associated with the input parameters.
%
%    See also SVDREGR.
%
%    Ref:
%    Detection of Influential Observation in Linear Regression
%    R. Dennis Cook
%    Technometrics, Vol. 19, No. 1, February 1977
%    pgs. 15-18

% Joe Henning - Jan 2014

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
   d = [];
   return
end
basis = basis_options{indx};

if (isempty(sigy))
   sigy = ones(size(y));
end

if (length(x) ~= length(y))
   fprintf('   Error ==> length(x) must be equal to length(y)\n');
   d = [];
   return
end

if (length(y) ~= length(sigy))
   fprintf('   Error ==> length(y) must be equal to length(sigy)\n');
   d = [];
   return
end

n = length(x);

% Calculate the full regression prediction and related parameters
p = svdregr(x, y, sigy, c, basis);

yhat = [];
for j = 1:n
   yhat(j) = 0;
   for i = 0:c
      switch basis
         case 'standard'
            yhat(j) = yhat(j) + p(i+1)*polyb(x(j),i);
         case 'cheby1'
            yhat(j) = yhat(j) + p(i+1)*cheby1b(x(j),i);
         case 'cheby2'
            yhat(j) = yhat(j) + p(i+1)*cheby2b(x(j),i);
         case 'hermite'
            yhat(j) = yhat(j) + p(i+1)*hermb(x(j),i);
         case 'laguerre'
            yhat(j) = yhat(j) + p(i+1)*lagb(x(j),i);
         case 'legendre'
            yhat(j) = yhat(j) + p(i+1)*legendb(x(j),i);
      end
   end
end

mse = sum((yhat-y).^2)/n;

% Calculate Cook's distances
nhat = 1:n;
d = [];
for k = 1:n
   ind = find(nhat ~= k);
   pi = svdregr(x(ind), y(ind), sigy(ind), c, basis);
   yhatk = [];
   for j = 1:n
      yhatk(j) = 0;
      for i = 0:c
         switch basis
            case 'standard'
               yhatk(j) = yhatk(j) + pi(i+1)*polyb(x(j),i);
         case 'cheby1'
               yhatk(j) = yhatk(j) + pi(i+1)*cheby1b(x(j),i);
         case 'cheby2'
               yhatk(j) = yhatk(j) + pi(i+1)*cheby2b(x(j),i);
         case 'hermite'
               yhatk(j) = yhatk(j) + pi(i+1)*hermb(x(j),i);
         case 'laguerre'
               yhatk(j) = yhatk(j) + pi(i+1)*lagb(x(j),i);
         case 'legendre'
               yhatk(j) = yhatk(j) + pi(i+1)*legendb(x(j),i);
         end
      end
   end
   d(k) = sum((yhat-yhatk).^2)/((c+1)*mse);
end

if nargin == 6
   figure;
   stem(nhat,d);
   grid;
   xlabel('Observation #');
   ylabel('Cook''s distance');
   title('Cook''s distance plot');
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
