function [newsrc,indremove,pkinds,mindis,nearpk] = removesecondarysrc(oldsrc,sigsta,minsrcsep)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [newsrc,indremove,pkinds,mindis,nearpk] = removesecondarysrc(oldsrc,sigsta,minsrcsep)
%
% This is the function to remove the secondary, minor sources that are 'close'
% in arrival time to a major source but smaller in amplitude, from the grouped
% sources in terms of triplets of arrival times at each station.
%
% --The idea is determine if each of the deconvolved triplet in terms of the 
%   positive peak like the template actually match with the nearst waveform 
%   peak at at least 1 station. If not, discard it; if so, then find if there 
%   are any later-deconvolved triplets that are fairly close to the target. 
%   For those close enough, examine if they are associated with the same 
%   nearest peak at all stations. Consider these as the secondary, minor
%   sources and remove them from the input sources
%
% --Revision 20220523, no need to require the deconvolved positive peak to be 
%   close enough (<= maxpkdis) to the nearest signal waveform peak. Because 
%   rather than to the signal, it HAS to be close to some extent to the peak
%   on the residual, which the algorithm iteratively targets. 
%
% --Revision 200220607, if at least at one station, the deconvolved positive 
%   peaks are pointing to the same peak, throw the secondary ones. In other
%   words, a different, and preserved source after this process would point to
%   DIFFERENT closest waveform peaks at ALL stations.
%
% --Add an option to ask if sources themselves are too close. 
%   This requires a hard threshold where sources close within that
%   range, you throw away the secondary sources.
%
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/05/16
% Last modified date:   2022/06/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('minsrcsep',[]);  % as default, do NOT require a min separation between sources


newsrc = oldsrc;
indremove = [];   % store indices of sources that will be removed

%at each station, find the all peaks of the waveform
nsta = size(sigsta,2);
pkinds = cell(nsta,1);  % indices of peaks of the waveform
medpksep = zeros(nsta,1); % median separation of waveform peaks
ranpksep = zeros(nsta,2); % median separation of waveform peaks 
for ista = 1: nsta
  [~, pkinds{ista,1}] = findpeaks(sigsta(:,ista));
  medpksep(ista) = median(diff(pkinds{ista,1}));
  ranpksep(ista,1:2) = minmax(diff(pkinds{ista,1})');
end
maxsep = mean(medpksep);  % mean separation of all stations

% %if 'maxpkdis' is not specified in the input, then use the min of waveform peak separation 
% if isempty(maxpkdis)
%   maxpkdis = ceil(ranpksep(:,1)'/2);
% end

%for each deconvolved triplets in terms of the arrival time of the positive peak of the template,
%find the nearest peak of the seimogram at each station
nsrc = size(oldsrc,1);  
mindis = zeros(nsrc, nsta);
nearpk = zeros(nsrc, nsta);
for ii = 1: nsrc
  for ista = 1: nsta
    pk = pkinds{ista,1};
    %for any of the triplets, find the nearest peak at each station, store the distance in samples
    %and index of that nearst peak
    [mindis(ii,ista), nearpk(ii,ista)] = min(abs(pk-oldsrc(ii,(ista-1)*2+1)));
  end
end

for ii = 1: nsrc
  if ismember(ii, indremove)
    continue
  end
  
%   %does this nearest peak match with the deconvolved peak at any station? i.e., the distance between
%   %them is no bigger than some allowed maximum, at least at one station
%   match = sum(mindis(ii,:) <= maxpkdis);
%   if match > 0    % at least 1 station is okay

    pkcomp = oldsrc(ii+1:end,:);  % get the other triplets in later iterations
    
    %find all the later deconvolved triplets that are fairly close to the target triplet, using the
    %previous computed mean waveform peak separation of all stations as the threshold
    ind = find(abs(pkcomp(:,1)-oldsrc(ii,1)) <= maxsep | ...
               abs(pkcomp(:,3)-oldsrc(ii,3)) <= maxsep | ...
               abs(pkcomp(:,5)-oldsrc(ii,5)) <= maxsep );
    ind = ind+ii; % convert to the index of the original group
    
    %test if the found close, later deconvolved triplets is associated with the same peak as the
    %target
    if ~isempty(ind)
      for jj = 1: length(ind)
        indtest = ind(jj);
        
        %does it match with the same peak as the target triplet does?
        match = sum(nearpk(indtest,:) == nearpk(ii,:)); %the number of matches

        %if at ANY station, the matched peaks are the same, AND the later deconvolved one is indeed
        %smaller in amplitude, then we will discard it in the end
        %In other words, if at >=1 station the nearest peak is the same one, we say it is
        %secondary source and discard it! So primary and preserved sources need to point
        %to DIFFERENT sources at ALL stations
        if isempty(minsrcsep)   % if sources themselves are not required to be well separated
          if match>=1 && sum(oldsrc(indtest,[2 4 6]),2) < sum(oldsrc(ii,[2 4 6]),2)
            indremove = [indremove; indtest];
          %just in case the later deconvolved actually has a higher amp  
          elseif match>=1 && sum(oldsrc(indtest,[2 4 6]),2) >= sum(oldsrc(ii,[2 4 6]),2)
            indremove = [indremove; ii];  
          end
        else  % if sources themselves cannot be too close
          %are sources too close? ie., their separation is smaller than a threshold?
          tooclose = (abs(oldsrc(indtest,1)-oldsrc(ii,1)) < minsrcsep) + ...
            (abs(oldsrc(indtest,3)-oldsrc(ii,3)) < minsrcsep) + ...
            (abs(oldsrc(indtest,5)-oldsrc(ii,5)) < minsrcsep);
%           if tooclose == 0
%             disp(jj,ii)
%           end
          if (match>=1 || tooclose>=1) && sum(oldsrc(indtest,[2 4 6]),2) < sum(oldsrc(ii,[2 4 6]),2)
            indremove = [indremove; indtest];
          %just in case the later deconvolved actually has a higher amp
          elseif (match>=1 || tooclose>=1) && sum(oldsrc(indtest,[2 4 6]),2) >= sum(oldsrc(ii,[2 4 6]),2)
            indremove = [indremove; ii];
          end
        end
        
      end
    end
%   else
%     %discard the deconvolved peak if it does NOT match with the nearest waveform peak at all stations
%     indremove = [indremove; ii];  
%   end  
end

%discard the sources that are determined to be too close and secondary compared to a major source
%NOTE that 'indremove' may contain duplicates
newsrc(indremove, :) = [];

% keyboard


%%%%%%%%%%%%%%%%%%%%%%%% 1st version, looks more complicated, but give the same result
% indremove = [];
% for ista = 1: nsta
%   [pkhgt, pkinds{ista,1}] = findpeaks(sigsta(:,ista));
%   medpksep(ista) = median(diff(pkinds{ista,1}));
% end
% maxsep = mean(medpksep);
% 
% %turn the sources in terms of the arrival time of zero-crossing to positve peaks
% pkindep = impindep;
% for ista = 1: nsta
%   pkindep(:,(ista-1)*2+1) = pkindep(:,(ista-1)*2+1)+ppeaks(ista)-zcrosses(ista);
% end
% for ii = 1: size(pkindep,1)
%   if ismember(ii, indremove)
%     continue
%   end
%   
%   match = 0;
%   for ista = 1: nsta
%     %for any of the triplets, find the nearest peak at each station
%     pk = pkinds{ista,1};
%     [mindis, nearpk(ii,ista)] = min(abs(pk-pkindep(ii,(ista-1)*2+1)));
%     
%     %does this nearest peak match with the deconvolved peak at any station?
%     match = match+(mindis<=6);
%   end
%   
%   if match > 0
%     pkcomp = pkindep(ii+1:end,:);
%     
%     %find the rest of deconvoled triplets that are fairly close to the target triplet
%     ind = find(abs(pkcomp(:,1)-pkindep(ii,1)) <= maxsep | ...
%       abs(pkcomp(:,3)-pkindep(ii,3)) <= maxsep | ...
%       abs(pkcomp(:,5)-pkindep(ii,5)) <= maxsep );
%     ind = ind+ii; % convert to the index of the original group
%     
%     %what are the nearest peak for each of the found
%     if ~isempty(ind)
%       for jj = 1: length(ind)
%         match = 0;
%         for ista = 1: nsta
%           %for any of the triplets, find the nearest peak at each station
%           pk = pkinds{ista,1};
%           indtest = ind(jj);
%           [mindis,nearpkrst] = min(abs(pk-pkindep(indtest,(ista-1)*2+1)));
%           
%           %is it matching with the same peak as the target triplet does?
%           match = match+(nearpkrst == nearpk(ii,ista));
%         end
%         
%         %if at all stations, the matched peak are the same, and the later deconvolved one is
%         %smaller in amplitude, then we will discard it in the end
%         if match == 3 && sum(pkindep(indtest,[2 4 6]),2) < sum(pkindep(ii,[2 4 6]),2)
%           indremove = [indremove; indtest];
%         end
%         
%       end
%     end
%   else
%     indremove = [indremove; ii];
%   end
%   
% end
%%%%%%%%%%%%%%%%%%%%%%%% 1st version, looks more complicated, but give the same result







