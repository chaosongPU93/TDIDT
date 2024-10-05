function prob = prob_geq(dataA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   prob = prob_geq(dataA)
% This function will sort 'dataA', and compute the empirical cumulative 
% probability of 'dataA' is greater equal to a particular value that is 
% in 'dataA'
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/08/16
% Last modified date:   2024/08/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aa=sort(dataA);
% bb=(length(aa):-1:1)/length(aa);
[unival, ~, ind] = unique(aa);
counts = accumarray(ind, 1);
countsgeq=zeros(length(counts),1);  %greater equal to 
for i=1:length(counts)
  countsgeq(i,1)=sum(counts(i:end));
end
% Combine the unique values and their counts into a two-column matrix
bb = [unival countsgeq/length(aa)];
prob = bb;