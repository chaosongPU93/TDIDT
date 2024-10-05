function [b0, b1, varb0, varb1, covb, r, chi2, Q] = lregr(x, y, sigy)

% LREGR Simple linear regression.
%    [B0, B1] = LREGR(X,Y,SIGY) calculates straight line coefficients
%       Y = B0 + B1*X
%    for N data points {X,Y} where SIGY are the uncertainties associated
%    with each measurement Y.  If SIGY is unknown (omitted), they are
%    assumed to be unity.
%
%    [B0,B1,VARB0,VARB1,COVB] = LREGR() also returns the variances in the
%    estimates of B0 and B1, respectively, and the covariance of B0 and B1.
%
%    [B0,B1,VARB0,VARB1,COVB,R] = LREGR() also returns the correlation
%    coefficient R between the uncertainty in B0 and the uncertainty in B1.
%
%    [B0,B1,VARB0,VARB1,COVB,R,CHI2,Q] = LREGR() also returns the
%    chi-square merit CHI2 and probability Q that a value of chi-square as
%    poor as CHI2 should occur by chance.

% Joe Henning - Fall 2011

c = 0;
if (nargin < 3)
   sigy = ones(size(y));
   c = 1;
end

if (isempty(sigy))
   sigy = ones(size(y));
end

if (length(x) ~= length(y))
   fprintf('   Error ==> length(x) must be equal to length(y)\n');
   b0 = NaN;
   b1 = NaN;
   varb0 = NaN;
   varb1 = NaN;
   r = NaN;
   chi2 = NaN;
   Q = NaN;
   return
end

if (length(y) ~= length(sigy))
   fprintf('   Error ==> length(y) must be equal to length(sigy)\n');
   b0 = NaN;
   b1 = NaN;
   varb0 = NaN;
   varb1 = NaN;
   r = NaN;
   chi2 = NaN;
   Q = NaN;
   return
end

n = length(x);
S = 0;
Sx = 0;
Sy = 0;
Sxx = 0;
Sxy = 0;
for i = 1:n
   S = S + 1/(sigy(i)*sigy(i));
   Sx = Sx + x(i)/(sigy(i)*sigy(i));
   Sy = Sy + y(i)/(sigy(i)*sigy(i));
   Sxx = Sxx + x(i)*x(i)/(sigy(i)*sigy(i));
   Sxy = Sxy + x(i)*y(i)/(sigy(i)*sigy(i));
end

delta = S*Sxx - Sx*Sx;
%b0 = (Sxx*Sy - Sx*Sxy)/delta;
b1 = (S*Sxy - Sx*Sy)/delta;
b0 = (Sy - Sx*b1)/S;

varb0 = Sxx/delta;
varb1 = S/delta;

covb = -Sx/delta;

r = -Sx/sqrt(S*Sxx);

chi2 = 0;
for i = 1:n
   chi2 = chi2 + (y(i) - b0 - b1*x(i))*(y(i) - b0 - b1*x(i))/(sigy(i)*sigy(i));
end

Q = gammainc((n-2)/2,chi2/2);

if (c == 1)
   varb0 = varb0*chi2/(n-2);
   varb1 = varb1*chi2/(n-2);
   Q = NaN;
end
