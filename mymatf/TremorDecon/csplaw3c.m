function [synths,mommax,sources,greensts]=csplaw3c(writes,winlen,skiplen,synth,Greens,b,...
  xygrid,sps,fracelsew,seed,timetype,ftrans,stas)
%   From plaw3b but this one shifts but up to a certain distance, with offsets taken from John Armbruster's grid.
%    Format is offPGSS, offPGSI, dx, dy.
%    That's at 40 sps, so multiply by 2.5 to get offsets in 100 sps.
%
% --The orginal version written by Allan assumes a uniform, random distribution in arrival
% time, this may not the most ideal case if you want to obtain the synthetics at 4th stas.
% Instead, we can assume the same distribution for the origin time. With tools like Taup or
% hypoinverse, travel time between a specific source and station can be predicted with a
% velocity model. Therefore, arrival time can be obtained by adding travel time to origin
% time at any station theoretically, including 4th stations. Thus the synthetics at 4th sta
% can be generated using the same way by adding template at that sta to each arrival time.
% --It should be expected that the bulk of the code is very similar to 'synthgen'.
% --2023/06/14 I changed his sign convention of 'ioff' to be the same as my routines.
%
%
% Return:
%   synths: synthetics, signal length * nsta * n sat level
%   mommax: max moment of sources
%   sources: source info, seismogram index (1); source location index (2); sat level (3)
%
% Coded bY Allan Rubin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('timetype','tarvl');
defval('ftrans','interpchao');
defval('stas',[]);

beta=1.+b/1.5;
betam1=beta-1;
Greenlenm1=size(Greens,1)-1;
nsta=size(Greens,2);
%rng('default');
rng(seed);
inwrites=0;
maxmom=1.0;
synths = zeros(winlen-skiplen+1,nsta,length(writes));
rn=rand(writes(end),2); %random numbers are seismogram index, fault patch index; this generates uniform randoms in [0,1]
size(rn)
% xygrid(:,1:2)=round(2.5*xygrid(:,1:2)); % *2.5 to get to 100 sps from 40.
xygrid(:,1:2)=round(sps/40*xygrid(:,1:2)); % *4 to get to 160 sps from 40.
nspots=size(xygrid,1)
ispot=round(nspots*rn(:,2)+0.499999999); %ispot has length rn. I guess +0.499999999 is to ensure >=1
rn(:,2)=ispot;
%ispot(1:10)
%xygrid(245:249,:)
%use 'off2space002' to interpolate for the travel time to sta 1
[locuse, ~] = off2space002(xygrid(:,1:2),sps,ftrans,0);
% --2023/06/14 I changed his sign convention of 'ioff' to be the same as my routines.
ioff=xygrid(ispot,1:2); %1 is PGSS; 2 is PGSI %THIS MINUS SIGN IS EMPIRICAL!!!
%size(ioff)
%ioff(1:5,:)

if strcmp(timetype,'tarvl')
  sources=zeros(writes(end),2,length(writes));  %seismogram index (1); source location index (2); sat level (3)
elseif strcmp(timetype,'tori')
  sources=zeros(writes(end),4,length(writes)); % 2 more cols for tori and ttrvl
end

ind=ceil(rn(:,1)*winlen); % this indicates the arrival time at sta 1 is uniform and random
if strcmp(timetype,'tarvl')
  %info of synthetic sources
  source=[ind ispot]; % return the source  %seismogram index (1); source location index (2); sat level (3)
elseif strcmp(timetype,'tori')
  %corresponding travel time according to src loc indices
  ttrvl=locuse(ispot,6);
  indttrvl = round(ttrvl*sps);  % in samples
  
  % %corresponding offset 12 and 13 according to src loc indices
  % rnoff=xygrid(ispot,1:2);
  
  % %use 'off2space002' to interpolate for the travel time to sta 1
  % [locuse, ~] = off2space002(rnoff,sps,ftrans,0);
  % ttrvl = locuse(:,6);
  % indttrvl = round(ttrvl*sps);  % in samples
  
  %the arrival time for each source at sta 1 is: origin time + travel time
  indtori = ind;
  indtarvl = indtori + indttrvl;
  
  %info of synthetic sources
  source=[indtarvl ispot indtori indttrvl]; % 4 cols
end
%sources(1:10,:)

ind = source(:,1); %now back to use 'ind' as arrival time to consistent with below

%if also wants synthetics at 4th stas
if nsta > 3
  for ista = 4: nsta
    %predict 'off14' at a 4th station using the plane fit model
    [~,off14] = pred_tarvl_at4thsta(stas(ista,:),ioff(:,1),ioff(:,2),ind);
    ioff = [ioff off14];
  end
end

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
  greenst(n,1)=ind(n);  % although this index is returned in 'sources', but i want to make the correspondance clear
  
  %only add src if its arrival index <= signal length
  if ind(n) <= winlen
    %Where it goes, +/- a shift
    if ind(n) <= 0   %if arrival index is smaller than 1
      trunc = -ind(n);
      synth(1: Greenlenm1-trunc, 1) = synth(1: Greenlenm1-trunc, 1) + ...
        Greens(2+trunc: end, 1)*momomin;
    elseif ind(n)+Greenlenm1 > size(synth,1) %if arrival index is larger than length
      trunc = ind(n)+Greenlenm1 - size(synth,1);
      synth(ind(n): ind(n)+Greenlenm1-trunc, 1) = ...
        synth(ind(n): ind(n)+Greenlenm1-trunc, 1) + ...
        Greens(1: end-trunc, 1)*momomin;
    else
      synth(ind(n): ind(n)+Greenlenm1, 1) = ...
      synth(ind(n): ind(n)+Greenlenm1, 1) + Greens(:,1)*momomin; %Greens(whichG,:) %(for interpolated)
    end
    if nspots == 1 && isequaln(xygrid(1,1:2), [0 0])  %that is, if all sources come from [0,0]
      for ista=2:nsta
        %if from same spot, there is no arrival time offset between stations
        synth(ind(n):ind(n)+Greenlenm1, ista)=synth(ind(n):ind(n)+Greenlenm1, ista)+...
          Greens(:, ista)*momomin; %Greens(whichG,:) %(for interpolated)
      end
    else
      %       ioff=xygrid(ispot(n),1:2);
      %       if n<4
      %           ispot(n)
      %           ioff
      %           size(ioff)
      %       end
      for ista=2:nsta %This assumes nsta=3
        % ind - ioff is the arrival time in index at sta 2 and 3
        if ind(n)-ioff(n,ista-1) <= 0   %if arrival index is smaller than 1
          trunc=-(ind(n)-ioff(n,ista-1));
          synth(1: Greenlenm1-trunc, ista) = synth(1: Greenlenm1-trunc, ista) + ...
            Greens(2+trunc: end, ista)*momomin;
        elseif ind(n)-ioff(n,ista-1)+Greenlenm1 > size(synth,1) %if arrival index is larger than length
          trunc=ind(n)-ioff(n,ista-1)+Greenlenm1 - size(synth,1);
          synth(ind(n)-ioff(n,ista-1):ind(n)-ioff(n,ista-1)+Greenlenm1-trunc, ista)= ...
            synth(ind(n)-ioff(n,ista-1):ind(n)-ioff(n,ista-1)+Greenlenm1-trunc, ista)+...
            Greens(1: end-trunc, ista)*momomin;
        else
          synth(ind(n)-ioff(n,ista-1):ind(n)-ioff(n,ista-1)+Greenlenm1, ista)=...
            synth(ind(n)-ioff(n,ista-1):ind(n)-ioff(n,ista-1)+Greenlenm1, ista)+...
            Greens(:, ista)*momomin;
        end
      end
    end
  end
  
  %If time to write it, add elsewhere events first. Assumes momomin=1.
  if(ismember(n,writes))
    for i=1:floor(fracelsew*n)
      for ista=1:nsta
        %note that below would change the arrival index at sta 1 when n reaches the number
        %to write, but should be okay given that 'source' is predefined
        % ind(n)=ceil(rnelsew(i,ista)*winlen);
        indelse=ceil(rnelsew(i,ista)*winlen);
        greenstelse(i,1) = indelse;
        %note that below just treat noise elsewhere as the same-amp template arrive at the
        %same time at all stations as at sta 1
        if indelse+Greenlenm1 > size(synth,1) %if arrival index is larger than length
          trunc = indelse+Greenlenm1 - size(synth,1);
          synth(indelse: indelse+Greenlenm1-trunc, ista) = ...
            synth(indelse: indelse+Greenlenm1-trunc, ista) + ...
            Greens(1: end-trunc, ista);
        else
          % synth(ista,ind(n):ind(n)+Greenlenm1)=synth(ista,ind(n):ind(n)+Greenlenm1)+Greens(ista,:);%*momomin;
          synth(indelse:indelse+Greenlenm1, ista)=synth(indelse:indelse+Greenlenm1, ista)+...
            Greens(:, ista);%*momomin;
        end
      end
    end
    inwrites=inwrites+1
    n
    for ista=1:nsta
      synths(:,ista,inwrites)=synth(skiplen:winlen,ista);
    end
    a=greenst(1:n,1);
    a(a>winlen)=[]; %discard those arrivals > signal length
    size(a,1)
    greensts{inwrites}{1}=a; % Chao, 2022/01/10, the starting index of each added template, context of full length
    if fracelsew~=0
      greensts{inwrites}{2}=greenstelse(1:floor(fracelsew*n), 1);
    end
    mommax(inwrites)=maxmom;
    b=source(1:n,:);
    b(b(:,1)>winlen,:)=[]; %make sure to discard related sources
    size(b,1)
    sources(1:size(b,1),:,inwrites) = b;
  end
end

