clear;

%http://www.stats.uwo.ca/faculty/bellhouse/Generalized%20Least%20Squares_DaeroKim.pdf

x = [0.345 0.287 0.251 0.225 0.207 0.186 0.161 0.132 0.084 0.060];

y = [367 311 295 268 253 239 220 213 193 192];

sigy = [17 9 9 7 7 6 6 6 5 5];

fprintf('Ordinary Least Squares:\n');
fprintf('Residuals:\n');
fprintf('    Min      1Q  Median      3Q     Max\n');
fprintf('-14.773  -9.319  -2.829   5.571  19.817\n');
fprintf('Coefficients:\n');
fprintf('            Estimate Std. Error t value Pr(>|t|)\n');
fprintf('(Intercept)   135.00      10.08    13.4 9.21e-07 ***\n');
fprintf('x             619.71      47.68    13.0 1.16e-06 ***\n\n');
fprintf('Residual standard error: 12.69 on 8 degrees of freedom\n');
fprintf('Multiple R-squared: 0.9548,     Adjusted R-squared: 0.9491\n');
fprintf('F-statistic: 168.9 on 1 and 8 DF, p-value: 1.165e-06\n');

[p, varp, covp, chi2, Q] = svdregr(x, y)

fprintf('Ordinary Weighted Least Squares:\n');
fprintf('Residuals:\n');
fprintf('    Min      1Q  Median      3Q     Max\n');
fprintf('-2.3230 -0.8842  0.0000  1.3900  2.3353\n');
fprintf('Coefficients:\n');
fprintf('            Estimate Std. Error t value Pr(>|t|)\n');
fprintf('(Intercept)  148.473      8.079   18.38 7.91e-08 ***\n');
fprintf('x            530.835     47.550   11.16 3.71e-06 ***\n\n');
fprintf('Residual standard error: 1.657 on 8 degrees of freedom\n');
fprintf('Multiple R-squared: 0.9397,     Adjusted R-squared: 0.9321\n');
fprintf('F-statistic: 124.6 on 1 and 8 DF, p-value: 3.71e-06\n');

[p2, varp, covp, chi2, Q] = svdregr(x, y, sigy)

fprintf('Iteratively Reweighted Least Squares:\n');
fprintf('Residuals:\n');
fprintf('     Min       1Q   Median       3Q      Max\n');
fprintf('-1.44014 -0.00437  0.37332  0.76038  1.39141\n');
fprintf('Coefficients:\n');
fprintf('            Estimate Std. Error t value Pr(>|t|)\n');
fprintf('(Intercept)  142.892      9.246   15.45 3.06e-07 ***\n');
fprintf('x            531.344     62.002    8.57 2.65e-05 ***\n\n');
fprintf('Residual standard error: 0.9362 on 8 degrees of freedom\n');
fprintf('Multiple R-squared: 0.9018,     Adjusted R-squared: 0.8895\n');
fprintf('F-statistic: 73.44 on 1 and 8 DF, p-value: 2.652e-05\n');

[p3, varp, covp, chi2, Q] = irsvdregr(x, y)

a1 = p(1) + p(2)*x;
a2 = p2(1) + p2(2)*x;
a3 = p3(1) + p3(2)*x;

figure;
plot(x,y,'ko',x,a1,'k--',x,a2,'k',x,a3,'r');
grid;
axis tight;
legend('orig','unweighted','WLS','IRWS','location','NorthWest');
