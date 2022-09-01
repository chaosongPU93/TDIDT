function [new,isave,idis] = RemoveDoubleCounting3(old,dtmin,indext,col,winlensec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the function to remove the double counting of events that occurs
% in a short time window, which is considered to be the maximum resolution
% to regard several detections as independent.
%
% different from 'RemoveDoubleCounting.m', which is the 
% The previous way to consider the duplicates, i.e. the same arrival seen
% by multiple overlapping windows is to view arrivals who are separated by
% less than a threshold as the same arrival, and choose the one with the 
% largest CC coef.
% But now we will use a different way, we choose one that are closer to the
% center of the window, i.e., fairly far away from the start and the end of
% the window. This issue can be seen by 'identify.m' later because some 
% arrivals survived due to the largest CC are located at the end of the window
% so that the col. 6 and 8 would be doubtful as the strongest arrival window
% involved into the computation may not be complete.
%
% Update 2021/06/21:
%   Still preserve the one with higher CC, unless it is too close to the start or
%   or end of the window
%
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
% First created date:   2021/03/22
% Last modified date:   2021/03/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% avoid double counting the same event by keeping the window whose strongest
%%% arrival is the closest to the center
%%% the min offset between most energetic individual arrivals depends on the duration of dipole in that freq. band
colmaint = col(1);     % col number of main arrival time
colcent = col(2);     % col number of the center of the window
colcc = col(3:end);     % col number of CC coef
nrow = size(old,1);

%%% construct a more complete result matrix of all nin windows
new = [];
nskip = 0;    % number of skipped indexes 
idis = [];    % index that is discarded
isave = [];  % index that is saved

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
        sttime = old(indtemp(ind2),colcent)-winlensec/2;
        edtime = old(indtemp(ind2),colcent)-winlensec/2;
        %if the main arrival time is too close to the edges of the window, choose the next largest
        %cc one
        if old(indtemp(ind2),colmaint)-sttime< dtmin || edtime-old(indtemp(ind2),colmaint)< dtmin
            indtemp2 = indtemp;
            indtemp(ind2) = [];
            [~,ind2] = max(sum(old(indtemp, colcc), 2)/length(colcc));
            new = [new; old(indtemp(ind2), :)];
            isave = [isave; indtemp(ind2)];
            idis = [idis; setdiff(indtemp2, indtemp(ind2))];
            nskip = length(ind1)+1;
        else
            new = [new; old(indtemp(ind2), :)];
            isave = [isave; indtemp(ind2)];
            indtemp2 = indtemp;
            indtemp2(ind2)=[];
            idis = [idis; indtemp2];
            nskip = length(ind1)+1;
        end
        
    end

end

%%%%%%%%%%%%% OLD algorithm, commented out 2021/06/21 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% indfocus =1;
% while indfocus <= nrow
%     indfocus = indfocus+ nskip;
%     if indfocus == nrow
%         new = [new; old(indfocus, :)];
%         isave = [isave; indfocus];
%         break
%     elseif indfocus > nrow
%         break
%     end    
%     indrange=indfocus+1: min(indfocus+indext, nrow);
%     ind1=find(abs(old(indrange, colmaint)-old(indfocus,colmaint)) <= dtmin);
%     %%% to convert to the original index in old, do 'indfocus+ind1',or
%     %%% indrange(ind1), equivalent
%     if isempty(ind1)
%         new = [new; old(indfocus, :)];
%         nskip = 1;
%         isave = [isave; indfocus];
%     else
%         %             i
%         indtemp = [indfocus; indfocus+ind1];
%         %[~,ind2] = max(sum(old(indtemp, colcc), 2)/length(colcc));
%         [~,ind] = min(abs(old(indtemp,colcent)-old(indtemp,colmaint)));
%         if length(ind) > 1
%             cc = old(indtemp(ind),colcc);
%             [~,aaa] = max(cc);
%             ind2 = ind(aaa);
%         else
%             ind2 = ind;
%         end
%         new = [new; old(indtemp(ind2), :)];
%         isave = [isave; indtemp(ind2)];
%         indtemp(ind2)=[];
%         idis = [idis; indtemp];
%         nskip = length(ind1)+1;
%     end
% 
% end


