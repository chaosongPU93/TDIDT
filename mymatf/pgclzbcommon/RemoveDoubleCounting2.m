function [new,isave,idis,dup,dupsel,dist] = RemoveDoubleCounting2(old,dtmin,indext,col)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the function to remove the double counting of events that occurs
% in a short time window, which is considered to be the maximum resolution
% to regard several detections as independent.
% 
% Use both the time AND distance constaint:
%   Assume there are several detections inside the dtmin, accept the one that
%   is closest to its own lfe family, based on the assumption that detections
%   are more reliable around the lfe with the its own rotation parameter. 

% Input:
%   old: old matrix with duplicates of the SAME day!
%   dtmin: min duration
%   indext: index extension, the index range that needs to be checked
%   col: col numbers that are timing and cc coefs
%
% Output:
%   new: new matrix without duplites
%   isave: index that saved
%   ireduce: index that discarded
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2019/11/26
% Last modified date:   2020/05/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% for testing %%%%%%%%
% old=hfday;
% dtmin=dtminhf;
% indext=indexthf;
% col=colnum;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% avoid double counting the same event by keep the highest CC coef window
%%% the min offset between most energetic individual arrivals depends on the duration of dipole in that freq. band
colmaint = col(1);     % col number of main arrival time
colcc = col(2:end);
nrow = size(old,1);

%%% construct a more complete result matrix of all nin windows
new = [];
nskip = 0;    % number of skipped indexes 
idis = [];    % index that is discarded
isave = [];  % index that is saved

dup = [];     % all detections that have duplicates in dtmin, i.e. more than 1
dupsel = [];   % the selected detections (e.g. max CC) among all that have duplicates in dtmin
dist = [];  % distance
indfocus =1;
while indfocus <= nrow
    indfocus = indfocus+ nskip;
    indtemp = [];
    if indfocus == nrow
        new = [new; old(indfocus, :)];
        isave = [isave; indfocus];
        break
    elseif indfocus > nrow
        break
    end    
    indrange=indfocus+1: min(indfocus+indext, nrow);
    ind1=find(abs(old(indrange, colmaint)-old(indfocus,colmaint)) <= dtmin);
    %%% to convert to the original index in old, do 'indfocus+ind1',or
    %%% indrange(ind1), equivalent
    if isempty(ind1)
        new = [new; old(indfocus, :)];
        nskip = 1;
        isave = [isave; indfocus];
    else
        %             i
        indtemp = [indfocus; indfocus+ind1];
        [~,ind2] = max(sum(old(indtemp, colcc), 2)/length(colcc));  % ind2 always has a size of 1
        new = [new; old(indtemp(ind2), :)];
        isave = [isave; indtemp(ind2)];
        indtemp2 = indtemp;
        indtemp2(ind2)=[];
        idis = [idis; indtemp2];
        nskip = length(ind1)+1;
        
        if size(ind1,1) >1
            dup = [dup; old(indtemp(ind1), :)];     % dup contains dupsel 
            dupsel = [dupsel; old(indtemp(ind2), :)];   % always size of 1
            inddiff = setdiff(ind1,ind2);
            disttmp = sqrt((old(indtemp(inddiff),1)-old(indtemp(ind2),1)).^2+...
                      (old(indtemp(inddiff),2)-old(indtemp(ind2),2)).^2);
            dist = [dist; disttmp];  % so the size of dist = size of dup - size of dupsel
        end
    end

end

% keyboard