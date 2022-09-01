function [impindep,imptripf,indtrip] = groupimptripdecon_refOLD(sigdecon,ampit,irccran,rcccat,off1i,off1iw,loff_max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [impindep,imppairf,indpair] = groupimptripdecon_refOLD(sigdecon,ampit,irccws,rcccat,off1i,off1iw,loff_max)
%
% ABANDONED 20220623, ONLY FOR TESTING, use 'groupimptripdecon_ref.m'
%
% Similar to 'groupimptriplets.m', this function is to group the nearest 
% impulses that are deconvolved 
% independently from 3 different stations into triplets, in the descending
% order of 'significance' of the impulses at a reference station, where the
% significance is quantified by the product of the amplitude of the impulse
% and the rcc (concatenated of multiple subwins) value at the impulse.
%
% This is another function similar to 'groupimptripdecon.m', in order to 
% group the nearest impulses that are deconvolved  independently from 3 
% different stations into triplets. But differently, this function would
% at each station, try to group the impulse that was deconvolved earliest
% among all impulses that are close enough to the target.
%
% --The difference from 'groupimptripdecon.m' is, it will try to group the
% impulses based on a fixed allowable searching boundary 'loff_max' 
% relative to a best alignment (offset) that is predetermined for 
% overlapping sub-windows decomposed from the entire window. If the best 
% alignment changes for consecutive subwins, the grouping range is actually
% 'moving' as well. This may permit the grouping of triplets that were 
% missed due to the single invariant grouping range used in 'groupimptriplets.m'.
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
% --You might use a slightly larger 'loff_max' here than 'groupimptriplets.m'
% and compare the grouped results
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
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/06/20
% Last modified date:   2022/06/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%different stations has different number of non-zero impulses
npair = sum(sigdecon(:,1)>0); % choose station 1 as the reference

%directly compute the moving window length
mwlen = size(sigdecon,1)-size(rcccat,1);

%info of impulses at sta 1
imp1 = [];
imp1(:,1) = find(sigdecon(:,1)>0);  % index of impulse arrival, after some alignment 'off1i'
imp1(:,2) = sigdecon(sigdecon(:,1)>0, 1);  % amp of impulse
imp1(:,3) = rcccat(imp1(:,1)-mwlen/2); % rcc value at the impulse
imp1(:,4) = imp1(:,3).*imp1(:,2);  % rcc*amp
decon1 = ampit{1};   % stored deconvolved index and amp in each iteration
% [~,indimp,inddecon1] = intersect(imp1(:,1),decon1(:,1),'stable');  % deconvolution order of impulses
[~,indimp,inddecon1] = intersect(decon1(:,1),imp1(:,1),'stable');  % deconvolution order of impulses
% isequal(decon1(indimp,1), imp1(inddecon1,1))
% isequal(sort(indimp), (1:npair)')

%sort them based on the order in deconvolution
imp1 = imp1(inddecon1,:);  

%impulse triplets
imptrip = zeros(npair, 4*3);
imptrip(:,1:4) = imp1(:,1:4);

%info of impulses at sta 2
imp2 = [];
imp2(:,1) = find(sigdecon(:,2)>0);  % index of impulse arrival, after some alignment 'off1i'
imp2(:,2) = sigdecon(sigdecon(:,2)>0, 2);
imp2(:,3) = rcccat(imp2(:,1)-mwlen/2);
imp2(:,4) = imp2(:,3).*imp2(:,2);
decon2 = ampit{2};   % stored deconvolved index and amp in each iteration
[~,~,inddecon2] = intersect(imp2(:,1),decon2(:,1),'stable');  % deconvolution order of impulses
%info of impulses at sta 3
imp3 = [];
imp3(:,1) = find(sigdecon(:,3)>0);  % index of impulse arrival, after some alignment 'off1i'
imp3(:,2) = sigdecon(sigdecon(:,3)>0, 3);
imp3(:,3) = rcccat(imp3(:,1)-mwlen/2);
imp3(:,4) = imp3(:,3).*imp3(:,2);
decon3 = ampit{3};   % stored deconvolved index and amp in each iteration
[~,~,inddecon3] = intersect(imp3(:,1),decon3(:,1),'stable');  % deconvolution order of impulses

for ip = 1: npair
  %find which subwin this ref impulse belongs to
  iwin = findwhichrange(imp1(ip,1),irccran);
  %use that subwin to compute the difference in alignment bewteen subwin and entire win
  offdiff = off1iw(iwin,:)-off1i;
  %adjust the impulse location by this alignment difference
  off12 = imp1(ip,1)-imp2(:,1);  % time offset between sta 2 and 1, be consistent with the sign
  off13 = imp1(ip,1)-imp3(:,1);
  ind2 = find(abs(off12 - offdiff(2)) <=loff_max);  % find the impulse close enough separately for sta 2 and 3
  ind3 = find(abs(off13 - offdiff(3)) <=loff_max);
  off12c = off12(ind2);
  off13c = off13(ind3);
  imp2c = imp2(ind2,:);
  imp3c = imp3(ind3,:);
  off23c = zeros(length(ind2), length(ind3)); % derived time offset23, off12-off13+off23=0, dependent!
  for ii = 1: length(ind2)
    for jj = 1: length(ind3)
      off23c(ii,jj) = off13c(jj)-off12c(ii);  % this meets off12-off13+off23=0 all time
%       off23ci(ii,jj) = imp2c(ii,1)-imp3c(jj,1);
    end
  end
%   isequaln(off23c, off23ci)
  
  %if we require the dependent offset 23 be close enough as well, 'offdiff(3)-offdiff(2)' looks a
  %bit counter-intuitive, but would be wrong otherwise 
  ind23 = find(abs(off23c - (offdiff(3)-offdiff(2)) ) <=loff_max);
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


end

indtrip = find(sum(imptrip==0,2)==0); %triplets that have corresponding peaks at all stations
imptripf = imptrip(indtrip,:);  % found triplets
%re-organize the information
impindep = imptripf(:, [1 2 5 6 9 10]); % info of impulse index and amp
%adding the time offset 12, 13 and 23, indicating location, off23 == off13-off12 == tarvl2-tarvl3
impindep(:,7:9) = [impindep(:,1)-impindep(:,3) impindep(:,1)-impindep(:,5) ...
                   impindep(:,3)-impindep(:,5)];  % again, be consistent in sign!


% keyboard                 
                 
                 
                 
                 
                 
                 