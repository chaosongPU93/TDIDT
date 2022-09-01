function [synths,mommax,sources]=plaw3c(writes,winlen,skiplen,synth,Greens,b,xygrid,fracelsew,seed)
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
rn=rand(writes(end),2); %random numbers are seismogram index, fault patch index; this generates uniform randoms in [0,1]
size(rn)
xygrid(:,1:2)=round(2.5*xygrid(:,1:2)); % *2.5 to get to 100 sps from 40.
nspots=size(xygrid,1)
ispot=round(nspots*rn(:,2)+0.499999999); %ispot has length rn. I guess +0.499999999 is to ensure >=1 
%ispot(1:10)
%xygrid(245:249,:)
ioff=-xygrid(ispot,1:2); %1 is PGSS; 2 is PGSI %THIS MINUS SIGN IS EMPIRICAL!!!
%size(ioff)
%ioff(1:5,:)
ind=ceil(rn(:,1)*winlen); % this indicates the arrival time at sta 1 is uniform and random
sources=[ind,ioff]; % return the source
%sources(1:10,:)
rnelsew=rand(floor(fracelsew*writes(end)),nsta);  % this generates uniform randoms in [0,1]
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
    %this just add templates arrived at sta 1 that are uniformly randomly in time
    synth(1,ind(n):ind(n)+Greenlenm1)=synth(1,ind(n):ind(n)+Greenlenm1)+Greens(1,:)*momomin; %Greens(whichG,:) %(for interpolated)
    if nspots == 1 %that is, if everything comes from the same location...
        for ista=2:nsta
            %if from same spot, there is no arrival time offset between stations
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
            % ind +ioff is the arrival time in index at sta 2 and 3
            if ind(n)+ioff(n,ista-1) <= 0   %if arrival index is smaller than 1 
                trunc=-(ind(n)+ioff(n,ista-1));
                synth(ista,1:Greenlenm1-trunc)=synth(ista,1:Greenlenm1-trunc)+Greens(ista,2+trunc:end)*momomin;
            elseif ind(n)+ioff(n,ista-1)+Greenlenm1 > size(synth,2) %if arrival index is larger than length
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
                %note that below would change the arrival index at sta 1 when n reaches the number
                %to write, but should be okay given that 'source' is predefined
                ind(n)=ceil(rnelsew(i,ista)*winlen);
                %note that below just treat noise elsewhere as the same-amp template arrive at the
                %same time at all stations as at sta 1
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
