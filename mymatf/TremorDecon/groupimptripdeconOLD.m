function [impindep,imptripf,indtrip] = groupimptripdeconOLD(sigdecon,ampit,rcc,loff_max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [impindep,imppairf,indpair] = groupimptripdeconOLD(sigdecon,ampit,rcc,loff_max)
%
% ABANDONED 20220623, ONLY FOR TESTING, use 'groupimptripdecon.m'
%
% This is another function similar to 'groupimptriplets.m', in order to 
% group the nearest impulses that are deconvolved  independently from 3 
% different stations into triplets. But differently, this function would
% at each station, try to group the impulse that was deconvolved earliest
% among all impulses that are close enough to the target.
%
% --2022/06/20, in 'groupimptriplets.m', we now add the flag for 'method',
%   where 'wtamp' would consider the weighted amplitude only, ie. group 
%   the combo that has the largest summed weighted amp; 'mix' would 
%   consider both the distance and amp, if there is no combo has an 
%   amplitude >= 120% of the rest, group the
%   closest in arrival time.
%
% --In this function, we could still choose one reference station randomly.
%   Rank the impulses at the reference in the deconvolution order and start
%   to group them one by one. Therefore, for any target impulse at the
%   reference, any other impulses close to the target were definitely
%   deconvolved, so there won't be another decision making. At the other 2 
%   stations, we find all candidates are close enough. Among them, choose 
%   the one that were deconvolved earliest. Finally, find the triplets that 
%   has a successful match at all stations.
%
% --However, there is a caveat that might occur, similar to using the 'wtamp'
%   option in 'groupimptriplets.m' that grouped triplets might reverse their
%   order in arrival time at different stations, before and/or after removing 
%   the secondary sources. In such scenario, our eyes tend to group differently
%   so that impulses close enough to each other would be grouped, that is why
%   we added the option 'mix'. This function can NOT prevent such scenario,
%   so we must check the number of cases out of all grouped sources, if the
%   fraction is small, we don't have to worry about it.
%
%   
%
%
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/06/20
% Last modified date:   2022/06/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%different stations has different number of non-zero impulses
npair = sum(sigdecon(:,1)>0); % choose station 1 as the reference

%directly compute the moving window length
mwlen = size(sigdecon,1)-size(rcc,1);

%info of impulses at sta 1
imp1 = [];
imp1(:,1) = find(sigdecon(:,1)>0);  % index of impulse, ascending order of arrival 
imp1(:,2) = sigdecon(sigdecon(:,1)>0, 1);  % amp of impulse
imp1(:,3) = rcc(imp1(:,1)-mwlen/2); % rcc value at the impulse
imp1(:,4) = imp1(:,3).*imp1(:,2);  % rcc*amp
decon1 = ampit{1};   % stored deconvolved index and amp in each iteration
% [~,~,inddecon1] = intersect(decon1(:,1),imp1(:,1),'stable');  % deconvolution order of impulses
[~,indimp,inddecon1] = intersect(decon1(:,1),imp1(:,1),'stable');  % deconvolution order of impulses
% isequal(imp1(indimp,1), decon(inddecon,1))
% isequal(sort(indimp), (1:npair)')

%sort them based on the order in deconvolution
imp1 = imp1(inddecon1,:);  

imptrip = zeros(npair, 4*3);
imptrip(:,1:4) = imp1(:,1:4);

%info of impulses at sta 2
imp2 = [];
imp2(:,1) = find(sigdecon(:,2)>0);
imp2(:,2) = sigdecon(sigdecon(:,2)>0, 2);
imp2(:,3) = rcc(imp2(:,1)-mwlen/2);
imp2(:,4) = imp2(:,3).*imp2(:,2);
decon2 = ampit{2};   % stored deconvolved index and amp in each iteration
[~,~,inddecon2] = intersect(imp2(:,1),decon2(:,1),'stable');  % deconvolution order of impulses
  
%info of impulses at sta 3
imp3 = [];
imp3(:,1) = find(sigdecon(:,3)>0);
imp3(:,2) = sigdecon(sigdecon(:,3)>0, 3);
imp3(:,3) = rcc(imp3(:,1)-mwlen/2);
imp3(:,4) = imp3(:,3).*imp3(:,2);
decon3 = ampit{3};   % stored deconvolved index and amp in each iteration
[~,~,inddecon3] = intersect(imp3(:,1),decon3(:,1),'stable');  % deconvolution order of impulses

for ip = 1: npair
  off12 = imp1(ip,1)-imp2(:,1);   % time offset between sta 2 and 1, be consistent with the sign
  off13 = imp1(ip,1)-imp3(:,1);
  ind2 = find(abs(off12)<=loff_max);  % find the impulse close enough separately for sta 2 and 3
  ind3 = find(abs(off13)<=loff_max);
  off12c = off12(ind2);
  off13c = off13(ind3);
  imp2c = imp2(ind2,:);
  imp3c = imp3(ind3,:);
  off23c = zeros(length(ind2), length(ind3)); % derived time offset23, off12-off13+off23=0, dependent!
  for ii = 1: length(ind2)
    for jj = 1: length(ind3)
      off23c(ii,jj) = off13c(jj)-off12c(ii);  % this meets off12-off13+off23=0 all time
    end
  end
  %if we require the dependent offset 23 be close enough as well
  ind23 = find(abs(off23c)<=loff_max);
  if isempty(ind23)
    continue
  else
    [sub2, sub3] = ind2sub(size(off23c),ind23); % convert linear indices to matrix subscripts
  end
  
  %choose the pair whose sum of decon order is smallest
  %NOTE: you can't do this independently, otherwise 'off23' is not guaranteed to be <=loff_max 
  [~,ind] = min(inddecon2(ind2(sub2))+inddecon3(ind3(sub3))); 
  imptrip(ip,(2-1)*4+1:2*4) = imp2(ind2(sub2(ind)),1:4);
  imp2(ind2(sub2(ind)),:) = [];   % mute the used ones
  inddecon2(ind2(sub2(ind))) = [];
  imptrip(ip,(3-1)*4+1:3*4) = imp3(ind3(sub3(ind)),1:4);
  imp3(ind3(sub3(ind)),:) = [];   % mute the used ones
  inddecon3(ind3(sub3(ind))) = [];
  
  if abs(imptrip(ip,5)-imptrip(ip,9))>loff_max
    disp(ip)
    break
  end
    
end

indtrip = find(sum(imptrip==0,2)==0); %pairs that have corresponding peaks at all stations
imptripf = imptrip(indtrip,:);  % found pairs
%re-organize the information
impindep = imptripf(:, [1 2 5 6 9 10]); % info of impulse index and amp
%adding the time offset 12, 13 and 23, indicating location, off23 == off13-off12 == tarvl2-tarvl3
impindep(:,7:9) = [impindep(:,1)-impindep(:,3) impindep(:,1)-impindep(:,5) ...
                   impindep(:,3)-impindep(:,5)];  % again, be consistent in sign!

% keyboard                 
                 
                 
                 
                 
                 
                 