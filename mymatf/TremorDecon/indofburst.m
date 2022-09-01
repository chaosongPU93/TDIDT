function [iets,idate,ibst,idate1,ibst1] = indofburst(trange,globalind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [iets,idate,ibst] = indofburst(trange,globalind)
%
% This function is to identify which ets (iets), which date of the ets 
% (idate), and which burst of all in the same date of the tremor bursts
% (ibst) in query based on its global index 'globalind' in the 'trange'.
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/08/10
% Last modified date:   2022/08/10 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dates = unique(trange(:,1));
years = unique(floor(dates/1000));
nets = length(years);

count = 0;

indtrack = zeros(size(trange,1), 4);

for i = 1: nets
  % dates in each ets
  year = years(i);
  datesets = dates(floor(dates/1000)==year);
  for j = 1: length(datesets)
    rangetemp = trange(trange(:,1)==datesets(j), :);
    for k = 1: size(rangetemp,1)
      count = count + 1;
      indtrack(count,:) = [i j k count];
    end
  end
end

[~,IA,~] = intersect(indtrack(:,end),globalind);
iets = indtrack(IA, 1);
idate = indtrack(IA, 2);
ibst = indtrack(IA, 3);

count1 = 0;
indtrack1 = zeros(size(trange,1), 3);
for i = 1: length(dates)
  rangetemp = trange(trange(:,1)==dates(i), :);
  for j = 1: size(rangetemp,1)
    count1 = count1 + 1;
    indtrack1(count1,:) = [i j count1];
  end
end

[~,IA,~] = intersect(indtrack1(:,end),globalind);
idate1 = indtrack1(IA, 1);
ibst1 = indtrack1(IA, 2);

% keyboard
    