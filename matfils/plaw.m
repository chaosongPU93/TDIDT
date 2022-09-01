function [synths, mommax]=plaw(writes,winlen,skiplen,synth,Greens,b)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
beta=1.+b/1.5;
betam1=beta-1;
Greenlenm1=length(Greens)-1;
rng('default');
rng(1);
inwrites=0;
maxmom=1.0;
synths=zeros(length(writes),winlen-skiplen+1);
for n=1:writes(end)
    %R=rand(1,3);
    ind=ceil(rand(1)*winlen);
    %whichG=ceil(rand(1)*100); For interpolated Green's function
    if beta > 100
        momomin=1;
    else
        momomin=(1-rand(1)).^(-1/betam1);
        if(momomin > maxmom)
            n
            maxmom=momomin
        end
    end
    synth(ind:ind+Greenlenm1)=synth(ind:ind+Greenlenm1)+Greens*momomin; %Greens(whichG,:) %(for interpolated)
    if(ismember(n,writes))
        inwrites=inwrites+1
        synths(inwrites,:)=synth(skiplen:winlen);
        mommax(unwrites)=maxmom;
    end
end



