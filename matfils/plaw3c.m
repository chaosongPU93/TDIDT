function [synths,mommax,sources]=plaw3b(writes,winlen,skiplen,synth,Greens,b,xygrid,fracelsew,seed)
%   From plaw3b but this one shifts but up to a certain distance, with offsets taken from John Armbruster's grid.
%    Format is offPGSS, offPGSI, dx, dy.
%    That's at 40 sps, so multiply by 2.5 to get offsets in 100 sps.
beta=1.+b/1.5;
betam1=beta-1;
Greenlenm1=size(Greens,2)-1;
nsta=size(Greens,1);
%rng('default');
rng(seed);
inwrites=0;
maxmom=1.0;
synths=zeros(nsta,length(writes),winlen-skiplen+1);
rn=rand(writes(end),2); %random numbers are seismogram index, fault patch index
size(rn)
xygrid(:,1:2)=round(2.5*xygrid(:,1:2)); % *2.5 to get to 100 sps from 40.
nspots=size(xygrid,1)
ispot=round(nspots*rn(:,2)+0.499999999); %ispot has length rn.
%ispot(1:10)
%xygrid(245:249,:)
ioff=-xygrid(ispot,1:2); %1 is PGSS; 2 is PGSI %THIS MINUS SIGN IS EMPIRICAL!!!
%size(ioff)
%ioff(1:5,:)
ind=ceil(rn(:,1)*winlen);
sources=[ind,ioff];
%sources(1:10,:)
rnelsew=rand(floor(fracelsew*writes(end)),nsta);
for n=1:writes(end)
    %Amplitude.  So far just using uniform for this subroutine.
    if beta > 100
        momomin=1;
    else
        momomin=(1-rand(1)).^(-1/betam1); %a different random number if moments vary.  Inefficient?  
        if(momomin > maxmom)
            n
            maxmom=momomin
        end
    end
    %Where it goes, +/- a shift
    synth(1,ind(n):ind(n)+Greenlenm1)=synth(1,ind(n):ind(n)+Greenlenm1)+Greens(1,:)*momomin; %Greens(whichG,:) %(for interpolated)
    if nspots == 1 %that is, if everything cpomes from the same location...
        for ista=2:nsta
            synth(ista,ind(n):ind(n)+Greenlenm1)=synth(ista,ind(n):ind(n)+Greenlenm1)+Greens(ista,:)*momomin; %Greens(whichG,:) %(for interpolated)
        end
    else
%       ioff=xygrid(ispot(n),1:2);
%       if n<4
%           ispot(n)
%           ioff
%           size(ioff)
%       end
        for ista=2:nsta %This assumes nsta=3
            if ind(n)+ioff(n,ista-1) <= 0
                trunc=-(ind(n)+ioff(n,ista-1));
                synth(ista,1:Greenlenm1-trunc)=synth(ista,1:Greenlenm1-trunc)+Greens(ista,2+trunc:end)*momomin;
            elseif ind(n)+ioff(n,ista-1)+Greenlenm1 > size(synth,2)
                trunc=ind(n)+ioff(n,ista-1)+Greenlenm1 - size(synth,2)
                synth(ista,ind(n)+ioff(n,ista-1):ind(n)+ioff(n,ista-1)+Greenlenm1-trunc)= ...
                synth(ista,ind(n)+ioff(n,ista-1):ind(n)+ioff(n,ista-1)+Greenlenm1-trunc)+Greens(ista,end-trunc)*momomin; 
            else
                synth(ista,ind(n)+ioff(n,ista-1):ind(n)+ioff(n,ista-1)+Greenlenm1)=synth(ista,ind(n)+ioff(n,ista-1):ind(n)+ioff(n,ista-1)+Greenlenm1)+Greens(ista,:)*momomin; 
            end
        end
    end
    %If time to write it, add elsewhere events first. Assumes momomin=1.
    if(ismember(n,writes))
        for i=1:floor(fracelsew*n)
            for ista=1:nsta
                ind(n)=ceil(rnelsew(i,ista)*winlen);
                synth(ista,ind(n):ind(n)+Greenlenm1)=synth(ista,ind(n):ind(n)+Greenlenm1)+Greens(ista,:);%*momomin;
            end
        end
        inwrites=inwrites+1
        n
        for ista=1:nsta
            synths(ista,inwrites,:)=synth(ista,skiplen:winlen);
        end
        mommax(inwrites)=maxmom;
    end
end
