function zcind = get_zero_crossing(signal,icent,num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to find the approximate zero crossings in the signal. For example, in cc,
%  you want to find the several zero-crossings close to the max cc. 
%
% INPUT:  
%   signal: signal with n points;
%   icent: center index;
%   num:    number of zero crossings on either side
% 
% OUTPUT:
%   zcind:  indexes of zero crossings, size of 2*num
%
% Chao Song, chaosong@princeton.edu
% First created date:   2019/10/01
% Last modified date:   2019/10/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zcf = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
if num == inf   % means get all zero crossing, ignore the center 
    zcind = zcf(signal);
    zcind = reshape(zcind,[],1);
    
else    % means get only the zc along the center
    zcind = zcf(signal);
    zccent = zcind-icent;
    
    izcpos = find(zccent>0,num,'first');
    zcpos = zccent(izcpos)+icent;
    zcpos = reshape(zcpos,[],1);
    
    izcneg = find(zccent<0,num,'last');
    zcneg = zccent(izcneg)+icent;
    zcneg = reshape(zcneg,[],1);
    
    zcind = [zcpos; zcneg];
    
end
    