function [f] = plt_srcdtdistCDF(dist,dift,sps)      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_srcdtdistCDF(dist,dift,sps)    
%
% A simple plot of the cdf of the distance and differential time between the
% Nth and (N-n)th sources into one panel, where n depends. It also does not 
% distinguish what kind of differential time 'dift' is here. 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/07/15
% Last modified date:   2022/07/15 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

figure
hold on; box on;
yyaxis left
[cdfval,x] = ecdf(dift/sps); %between Nth and (N-1)th source
plot(cdfval,x,'b','linew',1);
plot([0.5 0.5],[0 median(dift/sps)],'b--');
plot([0 0.5],[median(dift/sps) median(dift/sps)],'b--');
ylabel('Diff. relative origin time (s)');
xlabel('Empirical CDF');
yyaxis right
[cdfval,x] = ecdf(dist); %between Nth and (N-1)th source
plot(cdfval,x,'r','linew',1);
plot([0.5 0.5],[0 median(dist)],'r--');
plot([0.5 1],[median(dist) median(dist)],'r--');
ylabel('Dist. between consecutive sources (km)');
