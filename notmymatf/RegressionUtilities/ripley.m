function [p, varp, chi2, Q] = ripley(x, y, sigx, sigy, inter)

% RIPLEY Ripley's TLS Linear Fit
%    RIPLEY(X,Y,SIGX,SIGY,INTER) calculates the Total Least Squares linear
%    polyomial of N data points (X,Y) using Ripley's algorithm.  SIGX are
%    the uncertainties associated with measurements X, and SIGY are the
%    uncertainties associated with measurements Y.  If SIGY are unknown
%   (omitted), they are assumed to be unity.
%
%    If INTER is non-zero, the y-intercept is not confined to pass through
%    the origin.  It defaults to 1.
%
%    [P,VARP,CHI2,Q] = RIPLEYR() also returns VARP, the variances in the
%    estimates of P, the chi-square merit CHI2, and probability Q that a
%    value of chi-square as poor as CHI2 should occur by chance.
%
%    Ref:
%    Regression techniques for the detection of analytical bias
%    Ripley, B. D. and Thompson, M.
%    Analyst, 112, 1987
%    pgs. 177-183

% Joe Henning - Fall 2013

if (nargin < 4)
   sigy = ones(size(y));
   inter = 1;
if (nargin < 5)
   inter = 1;
end

if (isempty(sigx))
   sigx = zeros(size(x)) + 1E-12;
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

% seed the algorithm with least-squares slope
beta0 = sum((x-mean(x)).*(y-mean(y)))/sum((x-mean(x)).^2);

global X Y SIGX SIGY INTER ALPHA WT;
X = x;
Y = y;
SIGX = sigx;
SIGY = sigy;
INTER = inter;

s = wsos(beta0);

[beta1,fval,exitflag,output] = fminsearch(@wsos,beta0);

alpha = ALPHA;
beta = beta1;
wt = WT;

p = [alpha; beta];

a1 = sum(wt.*x.*x);
a2 = sum(wt.*(x - sum(wt.*x)/sum(wt)).^2);
if (inter)
   sealpha = sqrt(a1/(a2*sum(wt)));
   sebeta = sqrt(1/a2);
else
   sealpha = 0;
   sebeta = sqrt(1/a1);
end

varp = [sealpha; sebeta].^2;

resid = (y - alpha - beta*x).*sqrt(wt);
chi2 = sum(resid.^2);

Q = gammainc((n-2)/2,chi2/2);


function wtsum = wsos(beta)
% calculate the weighted sum of squares given beta
global X Y SIGX SIGY INTER ALPHA WT;
varx = SIGX.*SIGX;
vary = SIGY.*SIGY;
a = 1./(vary + beta*beta*varx);
WT = a;
ALPHA = 0;
if (INTER)
   ALPHA = (sum(a.*Y) - beta*sum(a.*X))/sum(a);
end
u = (vary.*X + beta*varx.*(Y - ALPHA)).*a;
b = Y - ALPHA - beta*u;
wtsum = sum((X-u).^2./varx + b.*b./vary);
