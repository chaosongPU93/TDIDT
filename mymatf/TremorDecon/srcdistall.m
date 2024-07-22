function [dtime,dlocxy,eucdist] = srcdistall(timevec,locxy,septran)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [dtime,dlocxy,eucdist] = srcdistall(timevec,locxy,septran)
%
% Obtain the distance between each LFE and all other LFEs, in the sequential
% order of origin time or arrival time. Can specify 
% the range of separation in time in terms of samples, say you only care 
% about distance between sources whose time separation is within [0 1] s, then 
% you can specify 'septran' to be [0 1]*sps
% --For each source, first find all other srcs whose diff time is within the
% range 'septran', then within these sources
% --2024/02/26, the current algorithm is actually smart. For each source i, you
% only need to consider later sources whose diff time is within the range.
% Because if j is selected, then when you look at j, i will be selected
% whatsoever due to symmetry. And the duplicates need to be removed for the
% final set to return as the distance from each unique pair of sources within
% time range.
% --With the current algorithm, diff time is always positive, no matter if you
% are using 'abs' or not.  If you want to focus on events within a certain 
% length of window X sec, meaning that the first and last events are differed
% by at most X sec, essentially the 'septran' should be [0 X/2] sec.
% --2024/04/12, for 'dtime', do not use 'abs' anymore!
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
  othersrc = i+1:nsrc;  %you only need to look later sources, otherwise the duplicates need be removed anyways
  temptime = timevec(othersrc);
  temploc = locxy(othersrc,:);
  if ~isempty(septran)
    ind = find(abs(temptime-timevec(i))>=septran(1) & abs(temptime-timevec(i))<=septran(2));   %closer in time, if not all
  else
    ind = 1:length(othersrc);
  end
  
  dtime{i} = temptime(ind)-timevec(i);  %not using abs anymore!
  dlocxy{i} = [temploc(ind,1)-locxy(i,1) temploc(ind,2)-locxy(i,2)];
  eucdist{i} = sqrt((temploc(ind,1)-locxy(i,1)).^2 + (temploc(ind,2)-locxy(i,2)).^2);
  
end

dtime = cat(1,dtime{:});
dlocxy = cat(1,dlocxy{:});
eucdist = cat(1,eucdist{:});



% %%%%%%%%%%%%%%%% LEGACY %%%%%%%%%%%%%%%%%
% %loop for every source
% tags=[];
% for i = 1: nsrc
%   if ~isempty(septran)
%     ind = find(abs(timevec-timevec(i))>=septran(1) & abs(timevec-timevec(i))<=septran(2));   %closer in time, if not all
%   else
%     ind = 1: nsrc;
%   end
%   ind=reshape(setdiff(ind,i), [],1);
%   
%   tag = [ones(length(ind),1)*i ind];
%   tags=[tags; tag];
% end
%   
% tags2 = [tags(:,2) tags(:,1)];
% [c,ia,ib]=intersect(tags,tags2,'rows');
% indsave = setdiff(1:length(tags), ia);
% 
% dtime{i} = timevec(tag(indsave,1)) - timevec(tag(indsave,2));
% dlocxy{i} = [locxy(tag(indsave,1), 1)-locxy(tag(indsave,2),1) ...
%   locxy(tag(indsave,1), 2)-locxy(tag(indsave,2),2)];
% eucdist{i} = sqrt((locxy(tag(indsave,1), 1)-locxy(tag(indsave,2),1)).^2 + ...
%   (locxy(tag(indsave,1), 2)-locxy(tag(indsave,2),2)).^2 );
% %%%%%%%%%%%%%%%% LEGACY %%%%%%%%%%%%%%%%%

  
  
  










