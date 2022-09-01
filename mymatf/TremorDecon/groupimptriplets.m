function [impindep,imptripf,indtrip] = groupimptriplets(sigdecon,rcc,loff_max,method,refsta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [impindep,imppairf,indpair] = groupimptriplets(sigdecon,rcc,loff_max,method)
%
% This is the function to group the nearest impulses that are deconvolved 
% independently from 3 different stations into triplets, in the descending
% order of 'significance' of the impulses at a reference station, where the
% significance is quantified by the product of the amplitude of the impulse
% and the rcc value at the impulse.
%
% --The idea is that, for every ranked impulse at the reference station 
% according to the significance (weighted rcc), find all pairs at the other
% 2 stations that are close enough. Among them, choose the one with highest
% weighted rcc. Finally, find the triplets that has a successful match at 
% all stations.
%
% --The way I did in the very begining is to for every ranked impulse at 
% the reference station, find the closest impulse at the other station 
% separately, no matter what is tha amplitude is.
%
% --You have 3 independent impulses, you could only get 2 independent time
% offset off12 and off13, and off23 is derived and dependent, therefore,
% off12-off13+off23=0 is always true!!
%
% --The prealignment between all stations has a rule for the sign. 
%   if off12 >0, sta 1 arrives later than sta 2, so you need to move 2 to
%   the right --> to align with 1. 'ampit' contains the indices of the 
%   impulses at each station and the offsets between each triplets. The
%   sign of the offset is set so if off12 >0, it means sta 1 arrives later
%   than sta 2, so you need to move 2 to the right --> to align with 1.
%   This means, it has the SAME sign rule as the prealignment between the
%   signals at 3 stations. Therefore, the net offset relative to the
%   origin (0,0) is the sum of the 2 parts.
%
% --2022/06/14, we realize the grouping can be imperfect, when there are two
%   triplets are close in arrival time, but their arrival time order could
%   be mixture at diff stations. Check the note 2022-06-10 for a plot of
%   example. So the caveat is, your eyes might say if you can actually 
%   mutually group them differently so that they both are very close in 
%   arrival time, but the algorithm is deciding otherwise due to the 
%   amplitude. So can we combine the distance and amplitude? If there is no
%   combination of triplets that has a substantially larger amplitude, we 
%   group based on the distance instead?
%
% --2022/06/20, add the flag for 'method', where 'wtamp' would consider the 
%   the weighted amplitude only, ie. group the combo that has the largest 
%   summed weighted amp; 'mix' would consider both the distance and amp,
%   If there is no combo has an amplitude >= 120% of the rest, group the
%   closest in arrival time.
%
% --2022/06/23, added the option of selecting different stations as the 
%   reference, not necessarily the 1st in the station list, eg., PGC  
%   
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/05/16
% Last modified date:   2022/06/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('method','wtamp');   % default is to group based on the max summed weighted amp
defval('refsta',1);   % default reference station is sta 1, ie. PGC

%different stations has different number of non-zero impulses
npair = sum(sigdecon(:,refsta)>0); 

%directly compute the moving window length
mwlen = size(sigdecon,1)-size(rcc,1);

nsta = size(sigdecon,2);

%cell matrix to store info at diff. sta
impsta = cell(nsta,1);

%info of impulses at sta 1
impref = [];
impref(:,1) = find(sigdecon(:,refsta)>0);  % index of impulse
impref(:,2) = sigdecon(sigdecon(:,refsta)>0, refsta);  % amp of impulse
impref(:,3) = rcc(impref(:,1)-mwlen/2); % rcc value at the impulse
impref(:,4) = impref(:,3).*impref(:,2);  % rcc*amp

%sort them based on the product of amplitude and rcc
impref = sortrows(impref,4,'descend');

%store them in the same order as the station list, ie., PGC, SSIB, SILB
impsta{refsta} = impref;

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
    impsta{ista} = imp;
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
imp3 = impsta{sta3};

for ip = 1: npair
  %here 'off12', 'off13' mean the offset between ref (1st) sta and 2nd sta, 3rd sta respectively
  off12 = impref(ip,1)-imp2(:,1);   % time offset between sta 2 and 1, be consistent with the sign
  off13 = impref(ip,1)-imp3(:,1);
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
  
  sumwtcoef = imp2c(sub2,4)+imp3c(sub3,4); % obtain sum of the weighted coefs at these pairs

  if isequal(method,'wtamp')
    [msumwtcoef,ind] = max(sumwtcoef);   % choose the pair has the highest weighted coef sum
    %assign the grouped impulses at the 2nd and 3rd sta accordingly
    imptrip(ip,(sta2-1)*4+1: sta2*4) = imp2(ind2(sub2(ind)),1:4);
    imp2(ind2(sub2(ind)),:) = [];   % mute the used ones
    imptrip(ip,(sta3-1)*4+1: sta3*4) = imp3(ind3(sub3(ind)),1:4);
    imp3(ind3(sub3(ind)),:) = [];   % mute the used ones
    
  elseif isequal(method,'mix')
    [sumwtcoefst,indst] = sort(sumwtcoef,'descend');
    %if the combo that has the max summed amp is at least 20% larger than the rest, then it is
    %significant
    if length(sumwtcoef)==1 || (length(sumwtcoef)>1 && sumwtcoefst(1) >= 1.2*sumwtcoefst(2))
      ind = indst(1);
      %assign the grouped impulses at the 2nd and 3rd sta accordingly
      imptrip(ip,(sta2-1)*4+1: sta2*4) = imp2(ind2(sub2(ind)),1:4);
      imp2(ind2(sub2(ind)),:) = [];   % mute the used ones
      imptrip(ip,(sta3-1)*4+1: sta3*4) = imp3(ind3(sub3(ind)),1:4);
      imp3(ind3(sub3(ind)),:) = [];   % mute the used ones
      
    %if there is no combo that has signifcantly high amp, then independently find the closest one at
    %each station. Technically, you can't do this independently, otherwise 'off23' is not guaranteed
    %to be <=loff_max, but for most cases, off23<=loff_max is met
    else
      [~,ind2c] = min(abs(off12c)); % find the closest impulse
      imptrip(ip,(sta2-1)*4+1: sta2*4) = imp2(ind2(ind2c),1:4);
      imp2(ind2(ind2c),:) = [];   % mute the used ones
      [~,ind3c] = min(abs(off13c)); % find the closest impulse
      imptrip(ip,(sta3-1)*4+1: sta3*4) = imp3(ind3(ind3c),1:4);
      imp3(ind3(ind3c),:) = [];   % mute the used ones
    end
    
  end
  
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
                 
                 
                 
                 
                 
                 