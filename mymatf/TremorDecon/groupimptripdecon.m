function [impindep,imptripf,indtrip] = groupimptripdecon(sigdecon,ampit,rcc,loff_max,refsta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [impindep,imppairf,indpair] = groupimptripdecon(sigdecon,ampit,rcc,loff_max)
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
% --2022/06/23, added the option of selecting different stations as the 
%   reference, not necessarily the 1st in the station list, eg., PGC  
%
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/06/20
% Last modified date:   2022/06/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('refsta',1);   % default reference station is sta 1, ie. PGC

%different stations has different number of non-zero impulses
npair = sum(sigdecon(:,refsta)>0); 

%directly compute the moving window length
mwlen = size(sigdecon,1)-size(rcc,1);

nsta = size(sigdecon,2);

%cell matrix to store info at diff. sta
impsta = cell(nsta,1);
ideconsta = cell(nsta,1);

%info of impulses at the reference sta
impref = [];
impref(:,1) = find(sigdecon(:,refsta)>0);  % index of impulse, ascending order of arrival 
impref(:,2) = sigdecon(sigdecon(:,refsta)>0, refsta);  % amp of impulse
impref(:,3) = rcc(impref(:,1)-mwlen/2); % rcc value at the impulse
impref(:,4) = impref(:,3).*impref(:,2);  % rcc*amp
deconref = ampit{refsta};   % stored deconvolved index and amp in each iteration
% [~,~,inddecon1] = intersect(decon1(:,1),imp1(:,1),'stable');  % deconvolution order of impulses
[~,indimp,ideconref] = intersect(deconref(:,1),impref(:,1),'stable');  % deconvolution order of impulses
% isequal(imp1(indimp,1), decon(inddecon,1))
% isequal(sort(indimp), (1:npair)')

%sort them based on the order in deconvolution at the ref station
impref = impref(ideconref,:);

%store them in the same order as the station list, ie., PGC, SSIB, SILB
impsta{refsta} = impref;
ideconsta{refsta} = ideconref;

%assign the sorted impulses at ref sta to the according columns of impulse triplets
imptrip = zeros(npair, 4*nsta);
imptrip(:,(refsta-1)*4+1:refsta*4) = impref(:,1:4);

%obtain the impulse info for the other 2 stations, in the same order as the station list, ie., PGC, SSIB, SILB 
for ista = 1: nsta
  if ista == refsta
    continue
  else
    imp = [];
    imp(:,1) = find(sigdecon(:,ista)>0);
    imp(:,2) = sigdecon(sigdecon(:,ista)>0, ista);
    imp(:,3) = rcc(imp(:,1)-mwlen/2);
    imp(:,4) = imp(:,3).*imp(:,2);
    decon = ampit{ista};   % stored deconvolved index and amp in each iteration
    [~,~,idecon] = intersect(imp(:,1),decon(:,1),'stable');  % deconvolution order of impulses
    impsta{ista} = imp;
    ideconsta{ista} = idecon;  
  end
end

%define the 2nd and 3rd station, given which sta is the ref
if refsta == 1      % if sta PGC is ref, 2nd and 3rd station is the same as the original order, SSIB and SILB
  sta2 = 2;
  sta3 = 3;
elseif refsta == 2  % if sta SSIB is ref, 2nd and 3rd station is sta PGC and SILB, respectively
  sta2 = 1;
  sta3 = 3;
else                % if sta SILB is ref, 2nd and 3rd station is sta PGC and SSIB, respectively
  sta2 = 1;
  sta3 = 2;
end
%impulses at the 2nd and 3rd stations are assigned accordingly
imp2 = impsta{sta2};
inddecon2 = ideconsta{sta2};
imp3 = impsta{sta3};
inddecon3 = ideconsta{sta3};

for ip = 1: npair
  %here 'off12', 'off13' mean the offset between ref (1st) sta and 2nd sta, 3rd sta respectively
  off12 = impref(ip,1)-imp2(:,1);
  off13 = impref(ip,1)-imp3(:,1);
  ind2 = find(abs(off12)<=loff_max);  % find the impulse close enough separately for sta 2 and 3
  ind3 = find(abs(off13)<=loff_max);
  off12c = off12(ind2);
  off13c = off13(ind3);
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
  
  %assign the grouped impulses at the 2nd and 3rd sta accordingly
  imptrip(ip,(sta2-1)*4+1: sta2*4) = imp2(ind2(sub2(ind)),1:4);
  imp2(ind2(sub2(ind)),:) = [];   % mute the used ones
  inddecon2(ind2(sub2(ind))) = [];
  imptrip(ip,(sta3-1)*4+1: sta3*4) = imp3(ind3(sub3(ind)),1:4);
  imp3(ind3(sub3(ind)),:) = [];   % mute the used ones
  inddecon3(ind3(sub3(ind))) = [];
  
  %a check for offset 23
  if abs(imptrip(ip,(sta2-1)*4+1)-imptrip(ip,(sta3-1)*4+1))>loff_max
    disp(ip)
    break
  end
  
end
  
indtrip = find(sum(imptrip==0,2)==0); %triplets that have corresponding peaks at all stations
imptripf = imptrip(indtrip,:);  % found triplets
%re-organize the information
impindep = imptripf(:, [1 2 5 6 9 10]); % info of impulse index and amp
%adding the time offset 12, 13 and 23, indicating location, off23 == off13-off12 == tarvl2-tarvl3
impindep(:,7:9) = [impindep(:,1)-impindep(:,3) impindep(:,1)-impindep(:,5) ...
                   impindep(:,3)-impindep(:,5)];  % again, be consistent in sign!

% keyboard                 
                 
                 
                 
                 
                 
                 