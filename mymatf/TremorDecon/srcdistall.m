function [dtime,dlocxy,eucdist] = srcdistall(timevec,locxy,septran)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [dtime,dlocxy,eucdist] = srcdistall(timevec,locxy,septran)
%
% Obtain the distance between each LFE and all other LFEs, in the sequential
% order of origin time or arrival time. Can specify 
% the range of separation in time in terms of samples, say you only care 
% about distance between sources whose time separation is within [0 1] s, then 
% you can specify 'septran' to be [0 1]*sps  
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/03/08 
% Last modified date:   2023/03/08 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

defval('septran',[]);

nsrc = size(locxy,1);

dtime = cell(nsrc, 1); %differential time, since it is sorted by time, it is always positive
dlocxy = cell(nsrc, 1); %differential location
eucdist = cell(nsrc, 1); %absolute Euclidean distance

%loop for every source
for i = 1: nsrc-1
  othersrc = i+1:nsrc;
  temptime = timevec(othersrc);
  temploc = locxy(othersrc,:);
  if ~isempty(septran)
    ind = find(abs(temptime-timevec(i))>=septran(1) & abs(temptime-timevec(i))<=septran(2));   %closer in time, if not all
  else
    ind = 1:length(othersrc);
  end
  
  dtime{i} = abs(temptime(ind)-timevec(i));  
  dlocxy{i} = [temploc(ind,1)-locxy(i,1) temploc(ind,2)-locxy(i,2)];
  eucdist{i} = sqrt((temploc(ind,1)-locxy(i,1)).^2 + (temploc(ind,2)-locxy(i,2)).^2);
  
end

dtime = cat(1,dtime{:});
dlocxy = cat(1,dlocxy{:});
eucdist = cat(1,eucdist{:});

