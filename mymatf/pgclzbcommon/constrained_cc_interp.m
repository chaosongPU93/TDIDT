function [off12,off13,cc,iloff,loff] = constrained_cc_interp(trace,mid,wlen,mshift,loffmax,ccmin,iup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to the constrained cross-correlation between the waveforms
% at the 3 stations to make sure off13-off12+off32=0, note that this is
% different from the independent cross-correlaton 'xcorr', which only align
% separately between two records. This is necessary to ensure only two offsets
% among the 3 station pairs are independent to represent detections
% 
% similar to 'constrained_cc.m', but this allows for interpolation, i.e. 
% upsampling to obtain a higher precision. More important for higher-frequency
% detections. If the 'iup' is set to 1, theoretically should return the same
% result as 'constrained_cc.m'.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/06/17
% Last modified date:   2021/06/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%in case it has non-zero mean
for i = 1:size(trace,1)
  trace(i,:) = detrend(trace(i,:));   
end
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

[xcmaxstack12n,imaxstack12]=max(sumstack12n,[],2);   %Integer-offset max cross-correlation
[xcmaxstack13n,imaxstack13]=max(sumstack13n,[],2);   % along row, max cc val and index in each window
[xcmaxstack32n,imaxstack32]=max(sumstack32n,[],2);

%Parabolic fit:
[xmaxstack12n,ymaxstack12n,aSTA12]=parabol(1,mshift,sumstack12n,imaxstack12); %Interpolated max cross-correlation
[xmaxstack13n,ymaxstack13n,aSTA13]=parabol(1,mshift,sumstack13n,imaxstack13);
[xmaxstack32n,ymaxstack32n,aSTA32]=parabol(1,mshift,sumstack32n,imaxstack32);

%Center them
imaxstack12cent=imaxstack12-mshift-1;  % "cent" is "centered"; imaxSTA12 is original 1: 2*mshift+1, corresponds to -mshift: mshift
imaxstack13cent=imaxstack13-mshift-1;
imaxstack32cent=imaxstack32-mshift-1;
%%% NOTICE: the right order of a closed 3-sta pair is +13, -12, +32, where 13 means 1-->3
iloff=imaxstack13cent-imaxstack12cent+imaxstack32cent; %How well does the integer loop close?
%
xmaxstack12n=xmaxstack12n-mshift-1;
xmaxstack13n=xmaxstack13n-mshift-1;
xmaxstack32n=xmaxstack32n-mshift-1;
loff=xmaxstack13n-xmaxstack12n+xmaxstack32n; %How well does the interpolated loop close?
xcmaxAVEn=(xcmaxstack12n+xcmaxstack13n+xcmaxstack32n)/3;
off12=xmaxstack12n;    % tmp == temporary
off13=xmaxstack13n;
off32=xmaxstack32n;


if xcmaxAVEn<ccmin || abs(loff)>loffmax || isequal(abs(imaxstack12cent),mshift)...
        || isequal(abs(imaxstack13cent),mshift) || isequal(abs(imaxstack32cent),mshift)
    % return the BAD offset and loopoff, and cc at offset (0,0)
    off12=mshift+1; off13=mshift+1; off32=mshift+1; %dummy them, if these criteria are met
%     cc = (sumstack12n(mshift+1)+sumstack13n(mshift+1)+sumstack32n(mshift+1))/3;
    cc = xcmaxAVEn;
    disp('WRONG! The basic criteria is not met');
else
    interpstack12n=interp(sumstack12n,iup,3);
    interpstack13n=interp(sumstack13n,iup,3);
    interpstack32n=interp(sumstack32n,iup,3);
    leninterp=length(interpstack12n);
    [xcmaxinterpstack12n,imaxinterpstack12]=max(interpstack12n(1:leninterp-(iup-1)));
    [xcmaxinterpstack13n,imaxinterpstack13]=max(interpstack13n(1:leninterp-(iup-1)));
    [xcmaxinterpstack32n,imaxinterpstack32]=max(interpstack32n(1:leninterp-(iup-1)));
    xcmaxconprev=-99999.;  %used to be 0; not good with glitches
    %predefine below in case none qualified are found, but replacable otherwise
    iSTA12bang=mshift+1;
    iSTA13bang=mshift+1;
    for iSTA12=max(1,imaxinterpstack12-3*iup):min(imaxinterpstack12+3*iup,iup*(2*mshift+1)-...
            (iup-1))
        %3 samples from peak;%intentionally wider than acceptable;%iup-1 are extrapolated points
        for iSTA13=max(1,imaxinterpstack13-3*iup):min(imaxinterpstack13+3*iup,iup*...
                (2*mshift+1)-(iup-1))
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
    if abs(iSTA12bang-imaxinterpstack12) <= loffmax*iup && ...
       abs(iSTA13bang-imaxinterpstack13) <= loffmax*iup && ...
       abs(iSTA32bang-imaxinterpstack32) <= loffmax*iup && ...
       interpstack12n(iSTA12bang)+interpstack13n(iSTA13bang)+interpstack32n(iSTA32bang) >= ...
       3*ccmin
   
        off12=(iSTA12bang-(iup*mshift+1))/iup;
        off13=(iSTA13bang-(iup*mshift+1))/iup;
        off32=(iSTA32bang-(iup*mshift+1))/iup;
                
        %%% xcmaxAVEnbang is added by Chao, to distinguish from
        %%% xcmaxAVEn, because it is max average CC coef
        cc=(interpstack12n(iSTA12bang)+interpstack13n(iSTA13bang)+interpstack32n(iSTA32bang))/3;
%     end
%     if off13-off12+off32 ~=0
    else
        disp('WRONG! Loopoff is not enclosed');
        off12=mshift+1; off13=mshift+1; off32=mshift+1; %dummy them, if these criteria are met
%         cc = (sumstack12n(mshift+1)+sumstack13n(mshift+1)+sumstack32n(mshift+1))/3;
        cc = xcmaxAVEn;

    end
end


% keyboard








