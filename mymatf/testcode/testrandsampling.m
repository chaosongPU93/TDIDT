%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is to test to generate random numbers from any kind of distribution
% (here use normal distribtuion as an example, but technically any kind is
% similar).
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/08/26
% Last modified date:   2020/08/26


%%% WAY 1
% so we could use the following codes to create a normal/Gaussian distribution
% with 'makedist' to create a prob. distri. object by specifying param. values
pd = makedist('Normal','mu',0,'sigma',1);
rng('default');
r1 = random(pd);

% Alternatively, you could directly generate random numbers from normal distribution 
rng('default');
r2 = random('Normal', 0, 1); % 1000 is number of samples   


%%% WAY 2
% OR, we could also directly use distribution-specific function 'normpdf' to 
% create the same object with same param. 
x = -3:0.1:3;
xpdf = normpdf(x,0,1); % return the pdf of norm distri. at locations in x


%%%
% we could use 'normrnd' to generate normal random numbers 
rng('default');  % save the current state of random number generator for reproducibility
r3 = normrnd(0,1,size(x));

s = rng; 
r4 = normrnd(0,1,[1,5]);

% use the saved the state of random number generator to shoudl give the same numbers as
% above
rng(s);
r5 = normrnd(0,1,[1,5]);




















