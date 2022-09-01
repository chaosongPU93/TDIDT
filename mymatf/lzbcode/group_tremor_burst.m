function [burst,nlfinb] = group_tremor_burst(interthf,tlf,ttol,ntolhf,ntollf)
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
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/12/27
% Last modified date:   2020/12/27

N = 10000;
burstraw = cell(1,N);
j = 1;
tmp = [];
for i = 1: size(interthf,1)
    if interthf(i,2) <= ttol
        tmp = [tmp; i];
    else
        if ~isempty(tmp)
            nhf = length(tmp);  % number of HF detections
            tmin = interthf(min(tmp),1);  % min and max range of time defined
            tmax = interthf(max(tmp),1);
            nlf = sum(tlf >= tmin & tlf <= tmax); % number of LF detections in the same time range
        else
            nhf = 0;
            nlf = 0;
        end
        
        if nhf >= ntolhf && nlf >= ntollf  % require both numbers to exceed their own thresholds
            burstraw{j} = tmp;
            tmp2(j) = nlf;
            j = j + 1;
        end
        tmp = [];
        continue
    end
    if j > N
        disp('array size is too small');
        break
    end
end

burst = [];   % tremor bursts
nlfinb = [];  % number of LF detections in the bursts
j = 1;
for i = 1: N
    if ~isempty(burstraw{i})
        burst{j} = burstraw{i}; % only save the non-empty ones
        nlfinb(j) = tmp2(i);
        j = j + 1;
    end
end

if ~isempty(burst)
    burst = burst';
end
    
