function [synths,mommax,sources,greensts]=plaw3d(writes,winlen,skiplen,synth,Greens,b,...
  xygrid,sps,fracelsew,seed,tdura)
%   From plaw3b but this one shifts but up to a certain distance, with offsets taken from John Armbruster's grid.
%    Format is offPGSS, offPGSI, dx, dy.
%    That's at 40 sps, so multiply by 2.5 to get offsets in 100 sps.
% The newer version uses plaw3d.m; the difference between that and plaw3c.m, used by synthshift.m,
% is in diff_plaw.  The structure of plaw has changed because of the need to remove 2 events from
% the same spot separated by less than the template duration.  So instead of just adding events as
% the seismograms become more saturated, for each saturation I needed to start making the catalog
% from scratch (zero events).
%
% Return:
%   synths: synthetics, signal length * nsta * n sat level
%   mommax: max moment of sources
%   sources: source info, seismogram index (1); source location index (2); sat level (3)
%
% Coded bY Allan Rubin, modified by Chao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta=1.+b/1.5;
betam1=beta-1;
Greenlenm1=size(Greens,1)-1;
nsta=size(Greens,2);
rng(seed);
maxmom=1.0;
synths=zeros(winlen-skiplen+1, nsta,length(writes));
% xygrid(:,1:2)=round(2.5*xygrid(:,1:2)); % *2.5 to get to 100 sps from 40.
xygrid(:,1:2)=round(sps/40*xygrid(:,1:2)); % *4 to get to 160 sps from 40.
nspots=size(xygrid,1)
sources=zeros(writes(end),2,length(writes));  %seismogram index (1); source location index (2); sat level (3)
rnelsew=rand(floor(fracelsew*writes(end)),nsta);
for n0=1:length(writes)
    rn=rand(writes(n0),2); %uniform randoms in [0,1], random numbers are seismogram index, fault patch index
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
            if rnsor(indx,1)-rnsor(indx-1,1) < tdura*sps %if separated by < 1 template duration
                %rnsor(indx,2)=round(nspots*rand+0.499999999);
                rnsor(indx,1)=rnsor(indx-1,1)+tdura*sps; %This could conceivably extend seismogram beyond end, esp. at high saturation
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
    greenst(1:sizernsor,1)=ind;
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
        synth(ind(n):ind(n)+Greenlenm1,1)=synth(ind(n):ind(n)+Greenlenm1,1)+Greens(:,1)*momomin; %Greens(whichG,:) %(for interpolated)
        if nspots == 1 %that is, if everything comes from the same location...
            for ista=2:nsta
                synth(ista,ind(n):ind(n)+Greenlenm1, ista)=synth(ind(n):ind(n)+Greenlenm1, ista)+...
                  Greens(:, ista)*momomin; %Greens(whichG,:) %(for interpolated)
            end
        else
            for ista=2:nsta %This assumes nsta=3
                if ind(n)+ioff(n,ista-1) <= 0
                    trunc=-(ind(n)+ioff(n,ista-1));
                    synth(1:Greenlenm1-trunc, ista)=synth(1:Greenlenm1-trunc, ista)+...
                      Greens(2+trunc:end, ista)*momomin;
                elseif ind(n)+ioff(n,ista-1)+Greenlenm1 > size(synth,1)
                    trunc=ind(n)+ioff(n,ista-1)+Greenlenm1 - size(synth,1)
                    synth(ind(n)+ioff(n,ista-1):ind(n)+ioff(n,ista-1)+Greenlenm1-trunc, ista)= ...
                        synth(ind(n)+ioff(n,ista-1):ind(n)+ioff(n,ista-1)+Greenlenm1-trunc, ista)+...
                        Greens(end-trunc, ista)*momomin;
                else
                    synth(ind(n)+ioff(n,ista-1):ind(n)+ioff(n,ista-1)+Greenlenm1, ista)=...
                      synth(ind(n)+ioff(n,ista-1):ind(n)+ioff(n,ista-1)+Greenlenm1, ista)+...
                      Greens(:, ista)*momomin;
                end
            end
        end
    end
    %add elsewhere events. Assumes momomin=1.
    for i=1:floor(fracelsew*n)
        for ista=1:nsta
            % ind(n)=ceil(rnelsew(i,ista)*winlen);
            % synth(ista,ind(n):ind(n)+Greenlenm1)=synth(ista,ind(n):ind(n)+Greenlenm1)+Greens(ista,:);%*momomin;
            indelse=ceil(rnelsew(i,ista)*winlen);
            greenstelse(i,1) = indelse;
            if indelse+Greenlenm1 > size(synth,1) %if arrival index is larger than length
              trunc = indelse+Greenlenm1 - size(synth,1);
              synth(indelse: indelse+Greenlenm1-trunc, ista) = ...
                synth(indelse: indelse+Greenlenm1-trunc, ista) + ...
                Greens(end-trunc, ista);
            else
              synth(indelse:indelse+Greenlenm1, ista)=synth(indelse:indelse+Greenlenm1, ista)+...
                Greens(:, ista);%*momomin;
            end
        end
    end
    for ista=1:nsta
        synths(:,ista,n0)=synth(skiplen:winlen,ista);
    end
    greensts{n0}{1}=greenst; % Chao, 2022/01/10, the starting index of each added template, context of full length
    if fracelsew~=0
      greensts{n0}{2}=greenstelse(1:floor(fracelsew*n), 1);
    end
    mommax(n0)=maxmom;
end
