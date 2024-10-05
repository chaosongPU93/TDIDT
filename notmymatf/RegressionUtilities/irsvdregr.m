function [p, varp, covp, chi2, Q] = irsvdregr(x, y, c, k)

% IRSVDREGR IRLS regression using Singular Value Decomposition
%    IRSVDREGR(X,Y,C,K) calculates iterative reweighted least squares
%    regression coefficients for N data points {X,Y}.  C specifies the
%    order of the basis functions (powers of x).  It defaults to 1 (linear
%    regression).  K specifies the order of the number of parameters in the
%    unknown variance model.  It defaults to 1 (linear variance).
%
%    [P,VARP,COVP] = IRSVDREGR() also returns VARP, the variances in the
%    estimates of P, and COVP, the covariance of P.
%
%    [P,VARP,COVP,CHI2,Q] = IRSVDREGR() also returns the chi-square merit
%    CHI2 and probability Q that a value of chi-square as poor as CHI2
%    should occur by chance.

% Joe Henning - Fall 2012

if (nargin < 3)
   c = 1;
   k = 1;
elseif (nargin < 4)
   k = 1;
end

n = length(x);

% Start with uniform sigmas
[p, varp, covp, chi2, Q] = svdregr(x, y, [], c);

MAXITER = 100;
EPS = 1e-8;

for m = 1:MAXITER
   chi2_old = chi2;

   yp = [];
   for i = 1:n
      yp(i) = 0;
      for j = 0:c
         yp(i) = yp(i) + p(j+1)*x(i)^j; 
      end
   end

   resid = (y - yp).^2;

   p2 = svdregr(y, resid, [], k);
   
   vary = [];
   for i = 1:n
      vary(i) = 0;
      for j = 0:k
         vary(i) = vary(i) + p2(j+1)*y(i)^j;
      end
   end
   
   sigy = sqrt(abs(vary));

   [p, varp, covp, chi2, Q] = svdregr(x, y, sigy, c);

   % calculate relative error
   dy = (chi2_old-chi2)/chi2_old;
   if (abs(dy) < EPS)
      break
   end
end

if (m >= MAXITER)
   fprintf('??? Max iterations exceeded.\n');
   p = NaN;
   varp = NaN;
   covp = NaN;
   chi2 = NaN;
   Q = NaN;
end
