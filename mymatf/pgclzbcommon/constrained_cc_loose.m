function [off12,off13,cc,iloff] = constrained_cc_loose(trace,mid,wlen,mshift,ccmin,iup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to 'constrained_cc_interp.m', this function is to carry out a 
% constrained cross-correlation between waveforms at 3 stations to best 
% align them. Unlike 'constrained_cc_interp.m' that finds the global max
% independently at each station and then reject all trials if the circuit
% of offsets from 3 indepmax CC does not meet the thresholds, here we try
% to align trace to our best. In detail, rather than looking at global max
% only, we track all ranked peaks in CC function and look at the triplets
% that meet the offset circuit and find the one with the hightest CC. 
% --This is likely to guarantee a solution (alignment) no worse than using
% a smaller max allowable 'mshift' when you try to increase 'mshift', in case
% that a larger 'mshift' leads to a higher gloabl max of CC, but the causing
% time shift is big.
%
% --The 1st version was to check every local peak, see if there is one best 
%   triplet that has the highest CC and the circuit of offsets is smaller than
%   'loffmax'; then within its proximity, find the triplet that whose circuit
%   is exactly 0. The shortcoming is, the final solution may not necessarily
%   has the highest CC if there is another solution around other peaks that 
%   were abondoned in the 1st step. Therefore, when you allow a larger 'mshift',
%   the solution could be inconsistent.
% --The 2nd version is to check EVERY possible triplets of offsets, and find
%   the one whose the circuit is 0, and has the highest CC. The offset could be
%   located at the boundary. In this manner, 'loffmax' became meaningless.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/09/19
% Last modified date:   2022/09/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stackauto = trace.*trace;
templen = size(trace,2);
lenx = templen-2*mshift;
stack12x = zeros(lenx, 2*mshift+1);
stack13x = zeros(lenx, 2*mshift+1);    % stas: 1->PGC, 2->SSIB, 3->SILB
stack32x = zeros(lenx, 2*mshift+1);

for n=-mshift:mshift
    % PGC corr SSIB, 1+mshift:tracelen-mshift == 1+mshift-n:tracelen-mshift-n == lenx
    stack12x(:,n+mshift+1)=trace(1,1+mshift:templen-mshift).* ...
        trace(2,1+mshift-n:templen-mshift-n);
    % PGC corr SILB
    stack13x(:,n+mshift+1)=trace(1,1+mshift:templen-mshift).* ...
        trace(3,1+mshift-n:templen-mshift-n);
    % SILB corr SSIB
    stack32x(:,n+mshift+1)=trace(3,1+mshift:templen-mshift).* ...
        trace(2,1+mshift-n:templen-mshift-n);
end

sumstack12=zeros(1,2*mshift+1);   % PGC-SSIB
sumstack13=zeros(1,2*mshift+1);   % PGC-SILB
sumstack32=zeros(1,2*mshift+1);   % SILB-SSIB
sumstack1sq=zeros(1,2*mshift+1);  % "sq" are the running cumulative sum, to be used later by differencing the window edpoints
sumstack2sq=zeros(1,2*mshift+1);
sumstack3sq=zeros(1,2*mshift+1);
sumstack3Bsq=zeros(1,2*mshift+1);

istart = mid-round(wlen/2)+1;
iend = istart+wlen;
sumstack12(1,:) = sum(stack12x(istart-mshift: iend-mshift, :));
sumstack13(1,:) = sum(stack13x(istart-mshift: iend-mshift, :));
sumstack32(1,:) = sum(stack32x(istart-mshift: iend-mshift, :));
sumstack1sq(1,:) = sum(stackauto(1, istart: iend));
sumstack3Bsq(1,:) = sum(stackauto(3, istart: iend));
for m = -mshift:mshift
    sumstack2sq(1,m+mshift+1)=sum(stackauto(2, istart-m: iend-m)); %+m??? (yes).
    sumstack3sq(1,m+mshift+1)=sum(stackauto(3, istart-m: iend-m));
end

%An attempt to bypass glitches in data.  Min value of good data typically ~10^{-2}
glitches=1.e-7;
sumstack1sq=max(sumstack1sq,glitches);
sumstack2sq=max(sumstack2sq,glitches);
sumstack3sq=max(sumstack3sq,glitches);    % return maximum between A and B

denomstack12n=realsqrt(sumstack1sq.*sumstack2sq);    % Real square root, An error is produced if X is negative
denomstack13n=realsqrt(sumstack1sq.*sumstack3sq);
denomstack32n=realsqrt(sumstack3Bsq.*sumstack2sq);

sumstack12n=sumstack12./denomstack12n;   % suffix 'n' means normalized
sumstack13n=sumstack13./denomstack13n;
sumstack32n=sumstack32./denomstack32n;


%find all peaks, ie, local max
[ppk12(:,2), ppk12(:,1)] = findpeaks(sumstack12n);
ppk12(:,1) = ppk12(:,1)-mshift-1;
ppk12 = sortrows(ppk12,2,'descend'); % sort the peak based on height
[ppk13(:,2), ppk13(:,1)] = findpeaks(sumstack13n);
ppk13(:,1) = ppk13(:,1)-mshift-1;
ppk13 = sortrows(ppk13,2,'descend'); 
[ppk32(:,2), ppk32(:,1)] = findpeaks(sumstack32n);
ppk32(:,1) = ppk32(:,1)-mshift-1;
ppk32 = sortrows(ppk32,2,'descend');
% figure
% subplot(311); hold on
% plot(-mshift:mshift,sumstack12n); 
% scatter(ppk12(:,1),ppk12(:,2),10);
% subplot(312); hold on
% plot(-mshift:mshift,sumstack13n);
% scatter(ppk13(:,1),ppk13(:,2),10);
% subplot(313); hold on
% plot(-mshift:mshift,sumstack32n);
% scatter(ppk32(:,1),ppk32(:,2),10);


%interpolation if needed
interpstack12n=interp(sumstack12n,iup,3);
interpstack13n=interp(sumstack13n,iup,3);
interpstack32n=interp(sumstack32n,iup,3);
xcmaxconprev=-99999.;  %used to be 0; not good with glitches
% %predefine below in case none qualified are found, but replacable otherwise
iSTA12bang=iup*(2*mshift+1)+1;  %a value outside the possible range
iSTA13bang=iup*(2*mshift+1)+1;
for iSTA12 = 1: length(interpstack12n)
  %3 samples from peak;%intentionally wider than acceptable;%iup-1 are extrapolated points
  for iSTA13 = 1: length(interpstack13n)
    ibangon = (iup*mshift+1)-iSTA13+iSTA12;
    if ibangon >= 1 && ibangon<=iup*(2*mshift+1)
      xcmaxcon=interpstack12n(iSTA12)+interpstack13n(iSTA13)+interpstack32n(ibangon);
      if xcmaxcon > xcmaxconprev
        xcmaxconprev=xcmaxcon;
        iSTA12bang=iSTA12;
        iSTA13bang=iSTA13;
      end
      %             else
    end
  end
end
%%% will result in the max xcmaxcon and corresponding iSTA12,
%%% iSTA13, and save them into xcmaxconprev, iSTA12bang and iSTA13bang
iSTA32bang=(iup*mshift+1)-iSTA13bang+iSTA12bang;
if ~isequal(iSTA12bang, iup*(2*mshift+1)+1) && ~isequal(iSTA13bang, iup*(2*mshift+1)+1) && ...
    iSTA32bang >= 1 && iSTA32bang <= iup*(2*mshift+1) && ...
    interpstack12n(iSTA12bang)+interpstack13n(iSTA13bang)+interpstack32n(iSTA32bang) >= ...
    3*ccmin
  
  off12=(iSTA12bang-(iup*mshift+1))/iup;
  off13=(iSTA13bang-(iup*mshift+1))/iup;
  off32=(iSTA32bang-(iup*mshift+1))/iup;
  
  %%% xcmaxAVEnbang is added by Chao, to distinguish from
  %%% xcmaxAVEn, because it is max average CC coef
  cc=(interpstack12n(iSTA12bang)+interpstack13n(iSTA13bang)+interpstack32n(iSTA32bang))/3;
  iloff = 0;
else
  disp('WRONG! No enclosed alignment could be found');
  off12=mshift+1; off13=mshift+1; off32=mshift+1; %dummy them, if these criteria are met
  cc = (sumstack12n(mshift+1)+sumstack13n(mshift+1)+sumstack32n(mshift+1))/3;
  iloff = [];
  
end

% 
% %if there is no peak, global max would be located at the boundary for sure
% if isempty(ppk12) || isempty(ppk13) || isempty(ppk32)
%   % return the BAD offset and loopoff, and cc at offset (0,0)
%   off12=mshift+1; off13=mshift+1; off32=mshift+1; %dummy them, if these criteria are met
%   cc = (sumstack12n(mshift+1)+sumstack13n(mshift+1)+sumstack32n(mshift+1))/3;
%   iloff = [];
%   disp('WRONG! The best solution is at the boundary');
% else
% 
%   %%%check all possible peaks
%   off12pk = nan(size(ppk12,1), 1); % you can't go more than availble peaks
%   off13pk = nan(size(ppk13,1), 1);
%   off32pk = nan(size(ppk32,1), 1);
%   iloffpk = nan(size(ppk12,1), size(ppk13,1), size(ppk32,1)); % you can't go more than availble peaks
%   for ii = 1: size(ppk12,1)  % sta pair 12, so each dimension can be different, but only the combination matters
%     for jj = 1: size(ppk13,1)  % sta pair 13
%       for kk = 1: size(ppk32,1)  % sta pair 32
%         off12pk(ii) = ppk12(ii,1);
%         off13pk(jj) = ppk13(jj,1);
%         off32pk(kk) = ppk32(kk,1);
%         iloffpk(ii,jj,kk) = off13pk(jj)-off12pk(ii)+off32pk(kk);
%       end
%     end
%   end
%   %find all pairs that are close enough, choose the one with highest weighted CC
%   indclose = find(abs(iloffpk)<=loffmax);
%   if isempty(indclose)
%     % return the BAD offset and loopoff, and cc at offset (0,0)
%     off12=mshift+1; off13=mshift+1; off32=mshift+1; %dummy them, if these criteria are met
%     cc = (sumstack12n(mshift+1)+sumstack13n(mshift+1)+sumstack32n(mshift+1))/3;
%     iloff = [];
%     disp('WRONG! No peaks found to be close enough');
% 
%   else
%     [sub12, sub13, sub32] = ind2sub(size(iloffpk),indclose); % convert linear indices to matrix subscripts
%     xcAVE = (ppk12(sub12,2)+ppk13(sub13,2)+ppk32(sub32,2))/3;
%     [cc,ind] = max(xcAVE);   % choose the pair has the highest CC
%     imaxstack12 = off12pk(sub12(ind));  % again, offset between a station pair has the sign
%     imaxstack13 = off13pk(sub13(ind));
%     imaxstack32 = off32pk(sub32(ind));
%     iloff = iloffpk(sub12(ind),sub13(ind),sub32(ind));  % again, average offset among all stations is an absolute value
%     
%   end
%   
%   %%%for the found triplet of peaks, check the surrounding if there is an exact solution
%   if cc<ccmin
%     % return the BAD offset and loopoff, and cc at offset (0,0)
%     off12=mshift+1; off13=mshift+1; off32=mshift+1; %dummy them, if these criteria are met
%     cc = (sumstack12n(mshift+1)+sumstack13n(mshift+1)+sumstack32n(mshift+1))/3;
%     iloff = [];
%     disp('WRONG! CC is too small');
%     
%   else
%     interpstack12n=interp(sumstack12n,iup,3);
%     interpstack13n=interp(sumstack13n,iup,3);
%     interpstack32n=interp(sumstack32n,iup,3);
%     xcmaxconprev=-99999.;  %used to be 0; not good with glitches
%     %predefine below in case none qualified are found, but replacable otherwise
%     iSTA12bang=mshift+1;
%     iSTA13bang=mshift+1;
%     for iSTA12=max(1,(imaxstack12+mshift+1)*iup-3*iup): ...
%         min((imaxstack12+mshift+1)*iup+3*iup,iup*(2*mshift+1)-(iup-1))
%       %3 samples from peak;%intentionally wider than acceptable;%iup-1 are extrapolated points
%       for iSTA13=max(1,(imaxstack13+mshift+1)*iup-3*iup): ...
%           min((imaxstack13+mshift+1)*iup+3*iup,iup*(2*mshift+1)-(iup-1))
%         ibangon = (iup*mshift+1)-iSTA13+iSTA12;
%         if ibangon >= 1 && ibangon<=iup*(2*mshift+1)
%           xcmaxcon=interpstack12n(iSTA12)+interpstack13n(iSTA13)+interpstack32n(ibangon);
%           if xcmaxcon > xcmaxconprev
%             xcmaxconprev=xcmaxcon;
%             iSTA12bang=iSTA12;
%             iSTA13bang=iSTA13;
%           end
%           %             else
%         end
%       end
%     end
%     %%% will result in the max xcmaxcon and corresponding iSTA12,
%     %%% iSTA13, and save them into xcmaxconprev, iSTA12bang and iSTA13bang
%     iSTA32bang=(iup*mshift+1)-iSTA13bang+iSTA12bang;
%     if abs(iSTA12bang-(imaxstack12+mshift+1)*iup) <= loffmax*iup && ...
%         abs(iSTA13bang-(imaxstack13+mshift+1)*iup) <= loffmax*iup && ...
%         abs(iSTA32bang-(imaxstack32+mshift+1)*iup) <= loffmax*iup && ...
%         interpstack12n(iSTA12bang)+interpstack13n(iSTA13bang)+interpstack32n(iSTA32bang) >= ...
%         3*ccmin
%       
%       off12=(iSTA12bang-(iup*mshift+1))/iup;
%       off13=(iSTA13bang-(iup*mshift+1))/iup;
%       off32=(iSTA32bang-(iup*mshift+1))/iup;
%       
%       %%% xcmaxAVEnbang is added by Chao, to distinguish from
%       %%% xcmaxAVEn, because it is max average CC coef
%       cc=(interpstack12n(iSTA12bang)+interpstack13n(iSTA13bang)+interpstack32n(iSTA32bang))/3;
%       %     end
%       %     if off13-off12+off32 ~=0
%     else
%       disp('WRONG! Loopoff is not enclosed');
%       off12=mshift+1; off13=mshift+1; off32=mshift+1; %dummy them, if these criteria are met
%       cc = (sumstack12n(mshift+1)+sumstack13n(mshift+1)+sumstack32n(mshift+1))/3;
%       iloff = [];
%       
%     end
%     
%   end
% 
% end

% figure
% subplot(311); hold on
% plot(-mshift:mshift,sumstack12n); 
% scatter(ppk12(:,1),ppk12(:,2),10);
% if ~isequal(off12,mshift+1)
%   scatter(off12, sumstack12n(off12+mshift+1),10,'b','filled');
% end
% subplot(312); hold on
% plot(-mshift:mshift,sumstack13n);
% scatter(ppk13(:,1),ppk13(:,2),10);
% if ~isequal(off13,mshift+1)
%   scatter(off13, sumstack13n(off13+mshift+1),10,'b','filled');
% end
% subplot(313); hold on
% plot(-mshift:mshift,sumstack32n);
% scatter(ppk32(:,1),ppk32(:,2),10);
% if ~isequal(off32,mshift+1)
%   scatter(off32, sumstack32n(off32+mshift+1),10,'b','filled');
% end
% 
% keyboard
% 







