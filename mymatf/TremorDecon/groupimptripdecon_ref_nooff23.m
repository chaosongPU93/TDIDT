function [impindep,imptripf,indtrip,sharp] = groupimptripdecon_ref_nooff23(sigdecon,ampit,irccran,rcccat,...
  off1i,off1iw,loff_max,refsta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [impindep,imptripf,indtrip,sharp] = groupimptripdecon_ref_nooff23(sigdecon,ampit,irccran,rcccat,...
%   off1i,off1iw,loff_max,refsta)
%
% All purposes the same as in 'groupimptripdecon_ref.m', the only difference here
% is to remove the constraint on 'off23'.
% 
% --We know that limiting 'off23' to be within the 'loff_max' would lead to 
% the weird zircon shape where the upper-left and lower-right corner in the 
% offset space would be truncated. This inevitably gives me an elongation which
% can be mixed with any potential elongation in deconvolution to data. 
% --To separate out two effects, this version of grouping would get rid of the 
% the constraint on 'off23'. As long as 'off12' and 'off13' is within 'loff_max',
% sources would be accepted, and there is no check on off23.
% 
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2024/04/17
% Last modified date:   2024/04/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('refsta',1);   % default reference station is sta 1, ie. PGC

%different stations has different number of non-zero impulses
npair = sum(sigdecon(:,refsta)>0); 

%directly compute the moving window length
mwlen = size(sigdecon,1)-size(rcccat,1);

nsta = size(sigdecon,2);

%cell matrix to store info at diff. sta
impsta = cell(nsta,1);
ideconsta = cell(nsta,1);

%info of impulses at sta 1
impref = [];
impref(:,1) = find(sigdecon(:,refsta)>0);  % index of impulse arrival, after some alignment 'off1i'
impref(:,2) = sigdecon(sigdecon(:,refsta)>0, refsta);  % amp of impulse
impref(:,3) = rcccat(impref(:,1)-mwlen/2); % rcc value at the impulse
impref(:,4) = impref(:,3).*impref(:,2);  % rcc*amp
deconref = ampit{refsta};   % stored deconvolved index and amp in each iteration
% [~,indimp,inddecon1] = intersect(imp1(:,1),decon1(:,1),'stable');  % deconvolution order of impulses
[~,indimp,ideconref] = intersect(deconref(:,1),impref(:,1),'stable');  % deconvolution order of impulses
% isequal(decon1(indimp,1), imp1(inddecon1,1))
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
    imp(:,3) = rcccat(imp(:,1)-mwlen/2);
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

%adjust the best alignment for each subwin and whole win, based on the ref sta, so that now the 2nd
%and 3rd col is the alignment between 2nd and ref, 3rd and ref sta, using the specified ref 
off1iwn = zeros(size(off1iw));
off1in = zeros(size(off1i));
if refsta == 1      % if sta PGC is ref, 2nd and 3rd station is the same as the original order, SSIB and SILB
  off1iwn = off1iw;
  off1in = off1i;
elseif refsta == 2  % if sta SSIB is ref, 2nd and 3rd station is sta PGC and SILB, respectively
  off1iwn(:,2) = -off1iw(:,2);
  off1iwn(:,3) = off1iw(:,3)-off1iw(:,2); % be careful about the sign here
  off1in(:,2) = -off1i(:,2);
  off1in(:,3) = off1i(:,3)-off1i(:,2);
else                % if sta SILB is ref, 2nd and 3rd station is sta PGC and SSIB, respectively
  off1iwn(:,2) = -off1iw(:,3);
  off1iwn(:,3) = off1iw(:,2)-off1iw(:,3); % be careful about the sign here
  off1in(:,2) = -off1i(:,3);
  off1in(:,3) = off1i(:,2)-off1i(:,3);
end

for ip = 1: npair
  %find which subwin this ref impulse belongs to
  iwin = findwhichrange(impref(ip,1),irccran);
  %use that subwin to compute the difference in alignment bewteen subwin and entire win
%   disp(ip)
  offdiff = off1iwn(iwin,:)-off1in;
  %here 'off12', 'off13' mean the offset between ref (1st) sta and 2nd sta, 3rd sta respectively
  off12 = impref(ip,1)-imp2(:,1);  % time offset between sta 2 and 1, be consistent with the sign
  off13 = impref(ip,1)-imp3(:,1);
  %adjust the impulse location by this alignment difference
  ind2 = find(abs(off12 - offdiff(2)) <=loff_max);  % find the impulse close enough separately for sta 2 and 3
  ind3 = find(abs(off13 - offdiff(3)) <=loff_max);
  off12c = off12(ind2);
  off13c = off13(ind3);
  off23c = zeros(length(ind2), length(ind3)); % derived time offset23, off12-off13+off23=0, dependent!
  for ii = 1: length(ind2)
    for jj = 1: length(ind3)
      off23c(ii,jj) = off13c(jj)-off12c(ii);  % this meets off12-off13+off23=0 all time
%       off23ci(ii,jj) = imp2c(ii,1)-imp3c(jj,1);
    end
  end
%   isequaln(off23c, off23ci)
  
  %In this version, we do NOT constrain the dependent offset 23! Note the diff with 'groupimptripdecon_ref'
  ind23 = find(abs(off23c - (offdiff(3)-offdiff(2)) ) <=inf); %set to infinity basically equals to no constraint on offset23, so that the actual grouping range is a rectangle 
  if isempty(ind23)
    continue
  else
    [sub2, sub3] = ind2sub(size(off23c),ind23); % convert linear indices to matrix subscripts
  end
  
  %choose the pair whose sum of decon order is smallest
  %NOTE: in 'groupimptripdecon', this 'ind' is jointly determined because otherwise 'off23'
  %is not guaranteed to be <=loff_max.
  % [~,ind] = min(inddecon2(ind2(sub2))+inddecon3(ind3(sub3)));
  % ind2min = ind;
  % ind3min = ind;
  %Since now there is no constraint 'off23', you can independently find the earliest-order impulse
  [~,ind2min] = min(inddecon2(ind2(sub2)));
  [~,ind3min] = min(inddecon3(ind3(sub3)));
  
  %assign the grouped impulses at the 2nd and 3rd sta accordingly
  imptrip(ip,(sta2-1)*4+1: sta2*4) = imp2(ind2(sub2(ind2min)),1:4);
  imp2(ind2(sub2(ind2min)),:) = [];   % mute the used ones
  inddecon2(ind2(sub2(ind2min))) = [];
  imptrip(ip,(sta3-1)*4+1: sta3*4) = imp3(ind3(sub3(ind3min)),1:4);
  imp3(ind3(sub3(ind3min)),:) = [];   % mute the used ones
  inddecon3(ind3(sub3(ind3min))) = [];

end

indtrip = find(sum(imptrip==0,2)==0); %triplets that have corresponding peaks at all stations
imptripf = imptrip(indtrip,:);  % found triplets
%re-organize the information
impindep = imptripf(:, [1 2 5 6 9 10]); % info of impulse index and amp
%adding the time offset 12, 13 and 23, indicating location, off23 == off13-off12 == tarvl2-tarvl3
impindep(:,7:9) = [impindep(:,1)-impindep(:,3) impindep(:,1)-impindep(:,5) ...
                   impindep(:,3)-impindep(:,5)];  % again, be consistent in sign!

%output the sharpness of the grouped triplets
sharp = zeros(size(impindep,1),3);
tmp = ampit{1,1};
if size(tmp,2) > 6
  for j = 1: 3
    tmp = ampit{1,j};
    for i = 1: size(impindep,1)
      [~,~,ish] = intersect(impindep(i,(j-1)*2+1),tmp(:,1),'stable');
      if length(ish)>1
        sharp(i,j) = mean(tmp(ish,7));  %if there are multiples, use mean
      else
        sharp(i,j) = tmp(ish,7);
      end
    end
  end
end

% keyboard                 
                 
                 
                 
                 
                 
                 