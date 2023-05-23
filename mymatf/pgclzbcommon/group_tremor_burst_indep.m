function [burst,ninb] = group_tremor_burst_indep(intert,ttol,ntol)
% 
% This function aims to group the indice of the tremor bursts
% accroding to the threshold on the inter-detection time and
% number of detections in the burst. For the time range defined
% by the detections whose inter-event time (interthf) is smaller than
% the threshold (ttol), the number of detections inside the range 
% for one data set (hf) and the other data set (lf) must both
% exceed their own minimum allowed (ntolhf,ntollf)
% 
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/03/14
% Last modified date:   2022/03/14

%% group all detections into bursts based on 'ttol', number of detections for each burst should be >=1
N = 10000;
burstraw = cell(N,1);   %a large enough cell storing the bursts, each of which is indice of detections
j = 1;  %counting eligible bursts
tmp = [];
for i = 1: size(intert,1)
  if intert(i,2) <= ttol  % 'intert' has 2 cols, time and time to preceding event
    %keep adding them into the same group
    tmp = [tmp; i];   % the globally 1st detection will be included if its inter time is set to 0
  elseif intert(i,2) > ttol  %if the index examined is distant from its preceding
    if ~isempty(tmp)
      if tmp(1)>1
        tmp = [tmp(1)-1; tmp];  % you always miss including the 1st of each group
      end
      n = length(tmp);  % number detections
      tmin = intert(min(tmp),1);  % min and max range of time defined
      tmax = intert(max(tmp),1);
      burstraw{j} = tmp;
      tmp2(j) = n;
      j = j + 1; % move on to the next burst
      tmp = []; % reset the group
      %         i = i-1;
    end
    
    %if its following detection is also distant, then it is isolated and should be in a group that only has itself
    if (i+1 <= size(intert,1) && intert(i+1,2) > ttol) || i == size(intert,1)
      n = 1;  % only has itself
      burstraw{j} = i;
      tmp2(j) = n;
      j = j + 1; % move on to the next burst
      tmp = []; % reset the group
    end
    
    if j > N
      disp('array size is too small');
      break
    end

  end
  
  if i == size(intert,1) && ~isempty(tmp)
    if tmp(1)>1
      tmp = [tmp(1)-1; tmp];  % you always miss including the 1st of each group
    end
    n = length(tmp);  % number detections
    tmin = intert(min(tmp),1);  % min and max range of time defined
    tmax = intert(max(tmp),1);
    burstraw{j} = tmp;
    tmp2(j) = n;
  end
    
%     elseif isempty(tmp)   %if the index examined is distant from its preceding and current group is empty
%         %if its following detection is also distant, then it is isolated and should be in a group that only has itself
%         if (i+1 <= size(intert,1) && intert(i+1,2) > ttol) || i == size(intert,1) 
%           n = 1;  % only has itself
%           burstraw{j} = i;
%           tmp2(j) = n;
%           j = j + 1; % move on to the next burst
%           tmp = []; % reset the group
%         %if it is close to the following, then it should be combined into the next group
%         elseif (i+1 <= size(intert,1) && intert(i+1,2) <= ttol)
% %           j = j + 1; % move on to the next burst
%           tmp = []; % reset the group
%           continue
%         end
%     end
    
%     %current group is complete, move on to the next
%     j = j + 1; % move on to the next burst
%     tmp = []; % reset the group

%     if n >= ntol  % require number to exceed the threshold
%       burstraw{j} = tmp;
%       tmp2(j) = n;
%       j = j + 1;
%     end
%     tmp = [];
      %         continue
%   end
  
%   if j > N
%     disp('array size is too small');
%     break
%   end
  
end

%% check if all detections have been grouped
indgrp = cell2mat(burstraw);
indall = (1: size(intert,1))';
indmiss = setdiff(indall,indgrp);
if ~isempty(indmiss)
  disp('Some detections are missed in grouping:');
  disp(indmiss);
end  
  
%% Apply the constraints of 'ntol'
burst = [];   % tremor bursts
ninb = [];    % number of detections in the bursts
j = 1;
for i = 1: N
  if ~isempty(burstraw{i}) && tmp2(i) >= ntol   % only save the eligible ones
    burst{j} = burstraw{i}; 
    ninb(j) = tmp2(i);
    j = j + 1;
  end
end

if ~isempty(burst)
  burst = burst';
  ninb = ninb';
end

% keyboard


    
