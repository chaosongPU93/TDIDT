function testchi2gof

bins = 0:5;
obsCounts = [6 16 10 12 4 2];
obsCounts1 = [5.5 16.3 9.5 11.8 4.1 2.5];

n = sum(obsCounts);
pd = fitdist(bins','Poisson','Frequency',obsCounts');
expCounts = n * pdf(pd,bins);
[h,p,st] = chi2gof(bins,'Ctrs',bins,...
                        'Frequency',obsCounts1, ...
                        'Expected',expCounts,...
                        'NParams',1,...
                        'EMin',0)
                    
[testrst,prob,stats] = chi2gof_wt(obsCounts1, expCounts, 1, 0, 0.05)


keyboard