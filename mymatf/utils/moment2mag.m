function [M_w1, M_w2, M_w3] = moment2mag(M_0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [M_w1, M_w2] = MOMENT2MAG(M_0)
% This script converts scalar moment in unit 'Nm' to moment magnitude. 
% 
% There are several definitions for moment-magnitude conversion
% 1. Note the equation comes from Hanks and Kanamori 1979
%   Mw = (2/3)*log10(M0*1.0e7)-10.7
% 2. Another equation comes from Kanamori 1977, after derivation, it gives
%   Mw = (2/3)*(log10(M0*1.0e7)-16.1)
% 3. In GCMT project, page bottom of https://www.globalcmt.org/CMTsearch.html
% "Note on calculation of moment magnitude: The moment 
% magnitude is calculated by this software using the formula of Kanamori (1977),
% MW = (2/3)*(log M0 - 16.1), where M0 is given in units of dyne-cm. Prior to 
% February 1, 2006, the quantity (2/3)*16.1 was rounded to the value 10.73. 
% For a small number of earthquakes, searches conducted after 2006/02/01 will 
% give values for MW that differ by 0.1 magnitude unit from values given by 
% searches prior to 2006/02/01."
% So essentially they are using below:
%   Mw = (2/3)*log10(M0*1.0e7)-10.73
%   
% The 2nd and 3rd way seems to be more consistent with Bostock et al. (2015)
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/08/08
% Last modified date:   2024/09/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if following EQ 1
M_w1 = (2/3)*log10(M_0*1.0e7)-10.7;

%if following EQ 2
M_w2 = (2/3)*(log10(M_0*1.0e7)-16.1);

%if following EQ 3
M_w3 = (2/3)*log10(M_0*1.0e7)-10.73;
