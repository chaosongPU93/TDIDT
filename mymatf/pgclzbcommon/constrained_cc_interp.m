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
staauto = trace.*trace;
templen = size(trace,2);
lenx = templen-2*mshift;
sta12x = zeros(lenx, 2*mshift+1);
sta13x = zeros(lenx, 2*mshift+1);    % stas: 1->PGC, 2->SSIB, 3->SILB
sta32x = zeros(lenx, 2*mshift+1);

for n=-mshift:mshift
    % PGC corr SSIB, 1+mshift:tracelen-mshift == 1+mshift-n:tracelen-mshift-n == lenx
    sta12x(:,n+mshift+1)=trace(1,1+mshift:templen-mshift).* ...
        trace(2,1+mshift-n:templen-mshift-n);
    % PGC corr SILB
    sta13x(:,n+mshift+1)=trace(1,1+mshift:templen-mshift).* ...
        trace(3,1+mshift-n:templen-mshift-n);
    % SILB corr SSIB
    sta32x(:,n+mshift+1)=trace(3,1+mshift:templen-mshift).* ...
        trace(2,1+mshift-n:templen-mshift-n);
end

sumsta12=zeros(1,2*mshift+1);   % PGC-SSIB
sumsta13=zeros(1,2*mshift+1);   % PGC-SILB
sumsta32=zeros(1,2*mshift+1);   % SILB-SSIB
sumsta1sq=zeros(1,2*mshift+1);  % "sq" are the running cumulative sum, to be used later by differencing the window edpoints
sumsta2sq=zeros(1,2*mshift+1);
sumsta3sq=zeros(1,2*mshift+1);
sumsta3Bsq=zeros(1,2*mshift+1);

istart = mid-round(wlen/2)+1;
iend = istart+wlen;
sumsta12(1,:) = sum(sta12x(istart-mshift: iend-mshift, :));
sumsta13(1,:) = sum(sta13x(istart-mshift: iend-mshift, :));
sumsta32(1,:) = sum(sta32x(istart-mshift: iend-mshift, :));
sumsta1sq(1,:) = sum(staauto(1, istart: iend));
sumsta3Bsq(1,:) = sum(staauto(3, istart: iend));
for m = -mshift:mshift
    sumsta2sq(1,m+mshift+1)=sum(staauto(2, istart-m: iend-m)); %+m??? (yes).
    sumsta3sq(1,m+mshift+1)=sum(staauto(3, istart-m: iend-m));
end

%An attempt to bypass glitches in data.  Min value of good data typically ~10^{-2}
glitches=1.e-7;
sumsta1sq=max(sumsta1sq,glitches);
sumsta2sq=max(sumsta2sq,glitches);
sumsta3sq=max(sumsta3sq,glitches);    % return maximum between A and B

denomsta12n=realsqrt(sumsta1sq.*sumsta2sq);    % Real square root, An error is produced if X is negative
denomsta13n=realsqrt(sumsta1sq.*sumsta3sq);
denomsta32n=realsqrt(sumsta3Bsq.*sumsta2sq);

sumsta12n=sumsta12./denomsta12n;   % suffix 'n' means normalized
sumsta13n=sumsta13./denomsta13n;
sumsta32n=sumsta32./denomsta32n;

[xcmaxsta12n,imaxsta12]=max(sumsta12n,[],2);   %Integer-offset max cross-correlation
[xcmaxsta13n,imaxsta13]=max(sumsta13n,[],2);   % along row, max cc val and index in each window
[xcmaxsta32n,imaxsta32]=max(sumsta32n,[],2);

%Parabolic fit:
[xmaxsta12n,ymaxsta12n,aSTA12]=parabol(1,mshift,sumsta12n,imaxsta12); %Interpolated max cross-correlation
[xmaxsta13n,ymaxsta13n,aSTA13]=parabol(1,mshift,sumsta13n,imaxsta13);
[xmaxsta32n,ymaxsta32n,aSTA32]=parabol(1,mshift,sumsta32n,imaxsta32);

%Center them
imaxsta12cent=imaxsta12-mshift-1;  % "cent" is "centered"; imaxSTA12 is original 1: 2*mshift+1, corresponds to -mshift: mshift
imaxsta13cent=imaxsta13-mshift-1;
imaxsta32cent=imaxsta32-mshift-1;
%%% NOTICE: the right order of a closed 3-sta pair is +13, -12, +32, where 13 means 1-->3
iloff=imaxsta13cent-imaxsta12cent+imaxsta32cent; %How well does the integer loop close?
%
xmaxsta12n=xmaxsta12n-mshift-1;
xmaxsta13n=xmaxsta13n-mshift-1;
xmaxsta32n=xmaxsta32n-mshift-1;
loff=xmaxsta13n-xmaxsta12n+xmaxsta32n; %How well does the interpolated loop close?
xcmaxAVEn=(xcmaxsta12n+xcmaxsta13n+xcmaxsta32n)/3;
off12=xmaxsta12n;    % tmp == temporary
off13=xmaxsta13n;
off32=xmaxsta32n;


if xcmaxAVEn<ccmin || abs(loff)>loffmax || isequal(abs(imaxsta12cent),mshift)...
        || isequal(abs(imaxsta13cent),mshift) || isequal(abs(imaxsta32cent),mshift)
    % return the BAD offset and loopoff, and cc at offset (0,0)
    off12=mshift+1; off13=mshift+1; off32=mshift+1; %dummy them, if these criteria are met
%     cc = (sumstack12n(mshift+1)+sumstack13n(mshift+1)+sumstack32n(mshift+1))/3;
    cc = xcmaxAVEn;
    disp('WRONG! The basic criteria is not met');
else
    interpsta12n=interp(sumsta12n,iup,3);
    interpsta13n=interp(sumsta13n,iup,3);
    interpsta32n=interp(sumsta32n,iup,3);
    leninterp=length(interpsta12n);
    [xcmaxinterpsta12n,imaxinterpsta12]=max(interpsta12n(1:leninterp-(iup-1)));
    [xcmaxinterpsta13n,imaxinterpsta13]=max(interpsta13n(1:leninterp-(iup-1)));
    [xcmaxinterpsta32n,imaxinterpsta32]=max(interpsta32n(1:leninterp-(iup-1)));
    xcmaxconprev=-99999.;  %used to be 0; not good with glitches
    %predefine below in case none qualified are found, but replacable otherwise
    iSTA12bang=mshift+1;
    iSTA13bang=mshift+1;
    for iSTA12=max(1,imaxinterpsta12-3*iup):min(imaxinterpsta12+3*iup,iup*(2*mshift+1)-...
            (iup-1))
        %3 samples from peak;%intentionally wider than acceptable;%iup-1 are extrapolated points
        for iSTA13=max(1,imaxinterpsta13-3*iup):min(imaxinterpsta13+3*iup,iup*...
                (2*mshift+1)-(iup-1))
            ibangon = (iup*mshift+1)-iSTA13+iSTA12;
            if ibangon >= 1 && ibangon<=iup*(2*mshift+1)
                xcmaxcon=interpsta12n(iSTA12)+interpsta13n(iSTA13)+interpsta32n(ibangon);
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
    if abs(iSTA12bang-imaxinterpsta12) <= loffmax*iup && ...
       abs(iSTA13bang-imaxinterpsta13) <= loffmax*iup && ...
       abs(iSTA32bang-imaxinterpsta32) <= loffmax*iup && ...
       interpsta12n(iSTA12bang)+interpsta13n(iSTA13bang)+interpsta32n(iSTA32bang) >= ...
       3*ccmin
   
        off12=(iSTA12bang-(iup*mshift+1))/iup;
        off13=(iSTA13bang-(iup*mshift+1))/iup;
        off32=(iSTA32bang-(iup*mshift+1))/iup;
                
        %%% xcmaxAVEnbang is added by Chao, to distinguish from
        %%% xcmaxAVEn, because it is max average CC coef
        cc=(interpsta12n(iSTA12bang)+interpsta13n(iSTA13bang)+interpsta32n(iSTA32bang))/3;
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








