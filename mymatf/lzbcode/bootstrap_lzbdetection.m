%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to bootstrap the current detection results from LZB trio.
% i.e., 
%   1. preturb the detection in samples of time assuming there is error in
%       detection;
%   2. preturb the filtering correction in samples of time (or magnitude and
%       azimuth?? i am not sure yet) assuming there is also error in passband
%       correction;
%   *** the above 2 preturbation assumes that there is error in each detection
%       but does not change the number of detection.
%
%   3. preturb the current selection of detections (and/or the selection of 
%       migrations), i.e., randomly choose a subset of detections or migrations
%       to carry out the same analysis to obtain the summary results, say 1/3,
%       assuming there is potentially bias in selecting detections, or bias in 
%       choosing thresholds in detection. Actually this 3rd perturbation is
%       essential because the thresholds are totally impirical (subjective), and
%       bootstrapping could mitigate this effect.
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




















