function testttest2_wt

% compare the function 'ttest2' in matlab and 'ttest2_wt'

clc
% %% EXAMPLE 1, from matlab
% load examgrades
% x = grades(:,1);
% y = grades(:,2);
% % use built in, assuming equal variance
% [h,p,ci,stats] = ttest2(x,y,'Vartype','equal')
% 
% % use my own function, use frequency weights, and equal variance 
% [testrst,prob,stats] = ttest2_wt_v2(x,ones(length(x),1),y,ones(length(y),1),0.05, 2, 1)
% 
% 
% p - prob


% %% EXAMPLE 2, random numbers from gaussian of different mean
% x = normrnd(0,1,[100,1]);
% y = normrnd(0.5,1,[150,1]);
% % use built in, assuming equal variance
% [h,p,ci,stats] = ttest2(x,y,'Vartype','equal')
% 
% % use my own function, use frequency weights, and equal variance 
% [testrst,prob,stats] = ttest2_wt_v2(x,ones(length(x),1),y,ones(length(y),1),0.05, 2, 1)
% 
% 
% p - prob


% %% EXAMPLE 3, from matlab, but without assuming equal variance
% load examgrades
% x = grades(:,1);
% y = grades(:,2);
% % use built in, assuming equal variance
% [h,p,ci,stats] = ttest2(x,y,'Vartype','unequal')
% 
% % use my own function, use frequency weights, and unequal variance 
% [testrst,prob,stats] = ttest2_wt_v2(x,ones(length(x),1),y,ones(length(y),1),0.05, 2, 0)
% 
% 
% p - prob


% %% EXAMPLE 4, random numbers from gaussian of different (same) mean and different (same) variance
% x = normrnd(0,5,[100,1]);
% y = normrnd(0,1.5,[150,1]);
% % use built in, assuming equal variance
% [h,p,ci,stats] = ttest2(x,y,'Vartype','unequal')
% 
% % use my own function, use frequency weights, and equal variance 
% [testrst,prob,stats] = ttest2_wt_v2(x,ones(length(x),1),y,ones(length(y),1),0.05, 2, 0)
% 
% 
% p - prob


%% EXAMPLE 5, random numbers from gaussian of different mean and different variance with weights
x = normrnd(0,5,[100,1]);
y = normrnd(1.5,1.5,[150,1]);

% use my own function, use frequency weights (2), and unequal variance (0) 
[testrst,prob1,stats] = ttest2_wt_v2(x,ones(length(x),1),y,ones(length(y),1),0.05, 2, 0)

% use my own function, use reliability weights (1), and unequal variance (0) 
[testrst,prob2,stats] = ttest2_wt_v2(x,ones(length(x),1),y,ones(length(y),1),0.05, 1, 0)

prob1 - prob2

% use my own function, use reliability weights (1), and unequal variance (0) 
[testrst,prob2,stats] = ttest2_wt_v2(x,4*ones(length(x),1),y,4*ones(length(y),1),0.05, 1, 0)

% use my own function, use frequency weights (2), and unequal variance (0) 
[testrst,prob2,stats] = ttest2_wt_v2(x,0.5*ones(length(x),1),y,0.5*ones(length(y),1),0.05, 2, 0)











