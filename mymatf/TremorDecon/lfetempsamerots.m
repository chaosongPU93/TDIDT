% lfetempsamerots.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given the result of the time correction to Michael Bostock's LFE catalog
% time from 'bostcattimecorrect.m', in terms of zero-crossing, max positive
% peak, max absolute peak, etc, we now can safely tackle the correct segment
% of seismogram for each single LFE, and also for stacking them to get the 
% templates at different stations (but for different stations other than PGC,
% you need to account for the travel time difference due to the station 
% location).
% For this script, we try to get the type of templates for all fams, where
% the parameters for shear-wave splitting, rotation to optimal particle motion
% are the same as those obtained from LFEs at fam 002. Then through the 
% comparison of the amplitude of templates, we would know 
% exactly 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/03/31
% Last modified date:   2022/03/31
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%