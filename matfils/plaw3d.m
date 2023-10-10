function [synths,mommax,sources]=plaw3d(writes,winlen,skiplen,synth,Greens,b,xygrid,fracelsew,seed)
%   From plaw3b but this one shifts but up to a certain distance, with offsets taken from John Armbruster's grid.
%    Format is offPGSS, offPGSI, dx, dy.
%    That's at 40 sps, so multiply by 2.5 to get offsets in 100 sps.
beta=1.+b/1.5;
betam1=beta-1;
Greenlenm1=size(Greens,2)-1;
nsta=size(Greens,1);
rng(seed);
maxmom=1.0;
synths=zeros(nsta,length(writes),winlen-skiplen+1);
xygrid(:,1:2)=round(2.5*xygrid(:,1:2)); % *2.5 to get to 100 sps from 40.
nspots=size(xygrid,1)
sources=zeros(writes(end),2,length(writes));
rnelsew=rand(floor(fracelsew*writes(end)),nsta);
for n0=1:length(writes)
    rn=rand(writes(n0),2); %random numbers are seismogram index, fault patch index
    sizern=size(rn,1)
    ispot=round(nspots*rn(:,2)+0.499999999); %ispot has length rn.
    rn(:,2)=ispot;
    ind=ceil(rn(:,1)*winlen);
    rn(:,1)=ind; %rn is now full of integers - seismogram index (1) and source location index (2)
    %
    rnsor=sortrows(rn,[2,1]);
    indx=1;
    thisspot=rnsor(indx,2);
    while (thisspot <= nspots) && (indx <= sizern-1)
        while isequal (rnsor(indx+1,2),thisspot) && (indx <= sizern-2) %-2 is so,ething of a kludge.  Last entry might not be distinct.
            indx=indx+1;
            if rnsor(indx,1)-rnsor(indx-1,1) < 32 %if separated by < 1 template duration
                %rnsor(indx,2)=round(nspots*rand+0.499999999);
                rnsor(indx,1)=rnsor(indx-1,1)+32; %This could conceivably extend seismogram beyond end, esp. at high saturation
            end
        end
        thisspot=thisspot+1;
        indx=indx+1;
    end
    %
    rnsor(rnsor(:,1)>winlen,:)=[];
    sizernsor=size(rnsor,1)
    ind=rnsor(:,1); %ind is now back to the original definition (seismogram index), but in order given by rnsor
%     indsor=sortrows(rnsor,1); %Just checking...
%     indsor(end-10:end,:)
%     Greenlenm1
%     synthsize=size(synth)
%     winlen
    ioff=-xygrid(rnsor(:,2),1:2); %1 is PGSS; 2 is PGSI %THIS MINUS SIGN IS EMPIRICAL!!!
    sources(1:sizernsor,:,n0)=rnsor;
    rnsor(1:20,:)
    ioff(1:20,:)
    for n=1:sizernsor %writes(n0)
        %Amplitude.  So far just using uniform for this subroutine.
%         if beta > 100
            momomin=1;
%         else
%             momomin=(1-rand(1)).^(-1/betam1); %a different random number if moments vary.  Inefficient?
%             if(momomin > maxmom)
%                 n
%                 maxmom=momomin
%             end
%         end
        %Where it goes, +/- a shift
        synth(1,ind(n):ind(n)+Greenlenm1)=synth(1,ind(n):ind(n)+Greenlenm1)+Greens(1,:)*momomin; %Greens(whichG,:) %(for interpolated)
        if nspots == 1 %that is, if everything comes from the same location...
            for ista=2:nsta
                synth(ista,ind(n):ind(n)+Greenlenm1)=synth(ista,ind(n):ind(n)+Greenlenm1)+Greens(ista,:)*momomin; %Greens(whichG,:) %(for interpolated)
            end
        else
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
    end
    for i=1:floor(fracelsew*n)
        for ista=1:nsta
            ind(n)=ceil(rnelsew(i,ista)*winlen);
            synth(ista,ind(n):ind(n)+Greenlenm1)=synth(ista,ind(n):ind(n)+Greenlenm1)+Greens(ista,:);%*momomin;
        end
    end
    for ista=1:nsta
        synths(ista,n0,:)=synth(ista,skiplen:winlen);
    end
    mommax(n0)=maxmom;
end
