function [dtime,dlocxy,eucdist]=srcdistNtoNm(timevec, locxy, m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [dtime, dist]=srcdistNtoNm(timevec, locxy, m)
%
% Function just to ease the computation of differential time and distance between
% sources N and N-1 until N and N-m. The result would be save to cell array for 
% each value in m. 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/02/20
% Last modified date:   2023/02/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dtime = cell(m, 1); %differential time, since it is sorted by time, it is always positive
dlocxy = cell(m, 1); %differential location
eucdist = cell(m, 1); %absolute Euclidean distance

for i = 1: m
  %%%between Nth and (N-i)th source
  dtime{i} = diffcustom(timevec, i,'forward');
  dlocxy{i} = [diffcustom(locxy(:,1), i,'forward') diffcustom(locxy(:,2), i,'forward')];
  eucdist{i} = sqrt(diffcustom(locxy(:,1),i,'forward').^2 + diffcustom(locxy(:,2),i,'forward').^2 );
end