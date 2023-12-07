function loc = abs_loc_on_axis(axlim,ratio)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loc = abs_loc_on_axis(axlim,ratio)
%
% This simple function is to obtain the ABSOLUTE location on the axis of the 
% plot, based on the current axis range, and a relative ratio from 0 to 1.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/11/04
% Last modified date:   2023/11/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axran = range(axlim);
loc = axlim(1)+ratio*axran;
