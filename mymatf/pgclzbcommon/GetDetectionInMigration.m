function [hfmig, lfmig] = GetDetectionInMigration(fam,date,sttime,edtime,hfall,lfall)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to return the detections inside the time window of
% a single migration
% 
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/07/29
% Last modified date:   2019/07/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fam=num2str(fam);
datestr = num2str(date);
yr = datestr(1:4);
jday = datestr(5:end);
hfday = hfall(hfall(:,5)==date, :);
lfday = lfall(lfall(:,5)==date, :);

indhf = find(hfday(:,end)>=sttime & hfday(:,end)<=edtime);
hfmig = hfday(indhf,:);
hfmig = sortrows(hfmig,[7,6,1,2]);

indlf = find(lfday(:,end)>=sttime & lfday(:,end)<=edtime);
lfmig = lfday(indlf,:);
lfmig = sortrows(lfmig,[7,6,1,2]);
