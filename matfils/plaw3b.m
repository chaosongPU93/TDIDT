function [synths,mommax,greensts]=plaw3b(writes,winlen,skiplen,synth,Greens,b,jig,fracelsew,seed)
%UNTITLED2 Summary of this function goes here
%   for 3 stations/templates, shifting by up to +/-mshift (jig)
%   plaw3b accepts an additional fraction, not counted in N, that comes from "elsewhere" (3 stations randomly placed)
beta=1.+b/1.5;
betam1=beta-1;
Greenlenm1=size(Greens,2)-1;
nsta=size(Greens,1);
%rng('default');
rng(seed);
inwrites=0;
maxmom=1.0;
synths=zeros(nsta,length(writes),winlen-skiplen+1);
rn=rand(writes(end),nsta); %random numbers are index, jig2, jig3; this generates uniform randoms in [0,1]
size(rn)
rnelsew=rand(floor(fracelsew*writes(end)),nsta);  % this generates uniform randoms in [0,1]
for n=1:writes(end)
    %Amplitude.  So far just using uniform for this subroutine.
    if beta > 100
        momomin=1; % "momomin" = "moment (of each event) over minimum moment"
    else
        momomin=(1-rand(1)).^(-1/betam1); %a different random number if moments vary.  Inefficient?  
        if(momomin > maxmom)
            n
            maxmom=momomin
        end
    end
    %Where it goes, +/- a shift
    ind=ceil(rn(n,1)*winlen); % this indicates the arrival time at sta 1 is uniform and random
%     %%%%% added by Chao, 2022/01/10 %%%%%%
%     while ind<skiplen	   % to ensure the preserved trace doesn't have contribution from skipped part
%       ind=ceil(rn(n,1)*winlen);
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    greenst(n,1)=ind;  % Chao, 2022/01/10, the starting index of each added template, context of full length
    %this just add templates arrived at sta 1 that are uniformly randomly in time
    synth(1,ind:ind+Greenlenm1)=synth(1,ind:ind+Greenlenm1)+Greens(1,:)*momomin; %Greens(whichG,:) %(for interpolated)
    if jig == 0
        for ista=2:nsta
            synth(ista,ind:ind+Greenlenm1)=synth(ista,ind:ind+Greenlenm1)+Greens(ista,:)*momomin; %Greens(whichG,:) %(for interpolated)
        end
    else
        for ista=2:nsta
            % -0.5 to make the randoms range be [0.5, 0.5] 
            ioff=round(2*jig*(rn(n,ista)-0.5)); %Shift the times at stations 2 and 3 w.r.t. 1, by a random number of samples (at 100 sps) up to mjig.
            if ind+ioff <= 0
                trunc=-(ind+ioff);  % need to truncate some parts of template then add to the start
                synth(ista,1:Greenlenm1-trunc)=synth(ista,1:Greenlenm1-trunc)+Greens(ista,2+trunc:end)*momomin;
            else
                synth(ista,ind+ioff:ind+ioff+Greenlenm1)=synth(ista,ind+ioff:ind+ioff+Greenlenm1)+Greens(ista,:)*momomin; 
            end
        end
    end
    %If time to write it, add elsewhere events first. Assumes momomin=1.
    if(ismember(n,writes))
        for i=1:floor(fracelsew*n)
            for ista=1:nsta
                ind=ceil(rnelsew(i,ista)*winlen); % note that rnelse is also uniform randoms
                %note that below just treat noise elsewhere as the same-amp template arrive at the
                %same time at all stations as at sta 1
                synth(ista,ind:ind+Greenlenm1)=synth(ista,ind:ind+Greenlenm1)+Greens(ista,:);%*momomin;
            end
        end
        inwrites=inwrites+1
        n
        for ista=1:nsta
            synths(ista,inwrites,:)=synth(ista,skiplen:winlen);
        end
        greensts{inwrites}=greenst(1:n,1); % Chao, 2022/01/10, the starting index of each added template, context of full length
        mommax(inwrites)=maxmom;
    end
end

% keyboard