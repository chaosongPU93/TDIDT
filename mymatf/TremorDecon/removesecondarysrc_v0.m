function [newsrc,indremove,pkinds,mindis,nearpk] = removesecondarysrc_v0(oldsrc,sigsta,maxpkdis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [newsrc,indremove,pkinds,mindis,nearpk] = removesecondarysrc_v0(oldsrc,sigsta,maxpkdis)
%
% version 0, I save it in case it can be inspiring in the future, but it 
% doesn't seem to be a nice algorithm
%
% This is the function to remove the secondary, minor sources that are 'close'
% in arrival time to a major source but smaller in amplitude, from the grouped
% sources in terms of triplets of arrival times at each station.
%
%
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/05/16
% Last modified date:   2022/05/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indremove = [];
for ii = 1: size(oldsrc,1)
  if ismember(ii, indremove)
    continue
  end
  srccomp = oldsrc(ii+1:end,:);
  minrat = 1; % min. amp ratio
  %from experiment, when 2 BP templates closer than 20 samples, the decon would recognize an
  %averaged larger one in the middle, meaning a big impulse that is surrounded by small
  %impulses in the +-10 vicinty could be ignored, since the truth could be 2 sources separated
  %by 20 samples
%   maxsep = 10;
  for ista = 1: nsta
    [pkhgt, pk] = findpeaks(sigsta(:,ista));
    medpksep(ista) = median(diff(pk));
  end
  maxsep = mean(medpksep);
  
  maxoffdiff = 2;   % maximum allowed difference in off12 and off13, ie., the location
%   ind = find(sum(srccomp(:,[4 8 12]),2) < sum(imppairf(ii,[4 8 12]),2)/minrat ...
%     & abs(srccomp(:,1)-imppairf(ii,1)) < maxsep ...
%     & abs((srccomp(:,1)-srccomp(:,5))-(imppairf(ii,1)-imppairf(:,5))) <= 1 ...
%     & abs((srccomp(:,1)-srccomp(:,9))-(imppairf(ii,1)-imppairf(:,9))) <= 1);
%   ind = find(sum(srccomp(:,[4 8 12]),2) < sum(imppairf(ii,[4 8 12]),2)/minrat ...
%     & abs(srccomp(:,1)-imppairf(ii,1)) < maxsep );
%   ind = find(sum(srccomp(:,[2 4 6]),2) < sum(impindep(ii,[2 4 6]),2)/minrat & ...
%     abs(srccomp(:,1)-impindep(ii,1)) < maxsep & ...
%     abs((srccomp(:,1)-srccomp(:,3))-(impindep(ii,1)-impindep(ii,3))) <= minoffdiff & ...
%     abs((srccomp(:,1)-srccomp(:,5))-(impindep(ii,1)-impindep(ii,5))) <= minoffdiff);
  ind = find(sum(srccomp(:,[2 4 6]),2) < sum(oldsrc(ii,[2 4 6]),2)/minrat ...
              & abs(srccomp(:,1)-oldsrc(ii,1)) < maxsep );
  ind = ind+ii; % convert to the index of the original group
  if ~isempty(ind)
    indremove = [indremove; ind];
  end
end

impindep(indremove, :) = [];
