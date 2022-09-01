function [synths,mommax,sources,greensts]=synthgen(writes,winlen,skiplen,synth,Greens,b,xgrid,...
  ygrid,pdfgrid,sps,fracelsew,seed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [data2dpad,xmesh,ymesh] = mat2dzeropad(xvec,yvec,data2d)
%
% This function takes in the templates (or green's function) 'green' and the
% a pre-noised time series 'synth' and output the syhthetic seimograms 
% composed by the stacking of many templates arrive at different times. 
%
% --Assuming there are multiple (defined by
% 'writes') sources randomly coming from a region which are subject to a 
% custom PDF that is defined by 'xgrid', 'ygrid' and 'pdfgrid' of offset in 
% samples between station pair 12 and 13 at any sampling rate. But do make 
% sure that the sampling rate is consistent with the 'synth' and 'green'.
%
% --Their origin time (onset time) has a uniform distribution that spans 
% within the input window length 'winlen'. We have learned from the testing
% code 'testcombinedist' that the distribution of arrival time at any station
% is not uniform as the origin time given a custom distribution of sources
% unless it is also uniform. Because if you bin the arrival time, a few very
% small and very large arrival bins seem to follow the same distribution of 
% the sources. However, if on the seismogram, the very begining or the end is
% never your interest, you can roughly view the arrival time at sta 1 as 
% uniform. But in this code, we aim to take origin time uniform and compute
% the arrival time at sta 1 and use the location in offsets for each source
% to obtain the arrival time at sta 2 and 3. 
%
% --To get the arrival time, we need the travel time estimate. John 
% Armbruster's original catalog does not give travel time or origin time, but
% hypoinverse does report the origin time given the arrival time at each 
% station. 
% One way to do it is, take John's existing solution's depth as
% the fixed depth in hypo, and redo the inversion for a few rounds, ideally
% the result should not change too much relative to the input. But the new 
% result has the travel time information that can be used to interpolate at
% the current grid.
% Another way is, use my own inverted grid which has the travel time estimate
% but less regular than John's and using a newer slab model to do the 
% transformation from samples to locations in 'locinterp002_4s.m', then cut 
% out the ellipse to get the detections inside. From them, you get the PDF in
% offsets and stick all the following procedure to this new PDF. However, this
% would require you to re-run 'tremorbursts002_4s.m', 'seisbursts002_4s.m', as
% the detections inside the boundary are changed, maybe the resulted tremor 
% burst windows, and of course the seismograms to plot.
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/04/13
% Last modified date:   2022/04/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% below contain some params the same as in 'plaw3c.m' from Allan
beta = 1.+b/1.5;
betam1 = beta-1;
Greenlen = size(Greens,1);
Greenlenm1 = Greenlen-1;
nsta = size(Greens,2);
%rng('default');
rng(seed);
inwrites = 0;
maxmom = 1.0;
synths = zeros(winlen-skiplen+1,nsta,length(writes));

%generate random samples for the origin time of sources inside the restricted region
rntori = rand(writes(end),1); 
size(rntori)
indtori = ceil(rntori*winlen); % origin time is uniform within the span of window length

%generate random samples for the location of the source
rnoff = zeros(writes(end),2);
for i = 1: writes(end)
  [rnoff(i,1),rnoff(i,2)] = pinky(xgrid(1,:)',ygrid(:,1),pdfgrid);
  
end

%use 'off2space002' to interpolate for the travel time to sta 1
[locuse, ~] = off2space002(rnoff,sps,'interpArmbreloc',0);
ttrvl = locuse(:,6);
indttrvl = round(ttrvl*sps);  % in samples

%the arrival time for each source at sta 1 is: origin time + travel time
indtarvl = indtori + indttrvl;

%info of synthetic sources
sources=[indtori indttrvl indtarvl rnoff]; % return the source
  
%signal from elsewhere  
rnelsew=rand(floor(fracelsew*writes(end)), nsta);  % this generates uniform randoms in [0,1]
% rnelsew = 0;

%% generate synthetics by convolution with green and their zero-crossing times
% for i = 1: length(writes)
%   
%   nsrc = writes(i);   % number of sources
%     
%   if beta > 100
%     momomin = 1*ones(nsrc,1);
%   else
%     momomin = (1-rand(nsrc,1)).^(-1/betam1); %a different random number if moments vary.  Inefficient?
%     if(momomin > maxmom)
%       maxmom = momomin
%     end
%   end
% 
%   greenst = indtarvl(1: nsrc);  % arrival time index, ie., the starting index of each added template, context of full length
%   %you don't need all impulses, only some of them contribute to the length of truncated record
%   indstsig = 1;
%   indedsig = winlen;
%   indsave = find(greenst>=skiplen-Greenlen+indstsig & greenst<=indedsig+skiplen-1);
%   greenst = greenst(indsave);     
% %   indzcs = sort(indzc)';  
%   greenzc(:,1) = greenst + tgreen(1) ;  % zero-crossing time index -1
%   
%   for ista = 2: nsta %This assumes nsta=3
%     % ind + rnoff is the arrival time in index at sta 2 and 3
%     greenst = indtarvl(1: nsrc) + rnoff(1: nsrc, ista-1);
%     greenst = greenst(greenst>=skiplen-greenlen+indstsig & greenst<=indedsig+skiplen-1);
%     greenzc(:,ista) = greenst + tgreen(ista-1) ;  % zero-crossing time index -1
%   end
%   
%   %make an impulse array
%   impamp = zeros(max(greenzc(:)), nsta);
%   for ista = 1: nsta 
%     for j = 1: size(greenzc,1)
%       impamp(greenzc(j), ista) = impamp(greenzc(j), ista) +1;
%     end
%   end
%   
%   sigconv = conv(Greens(:, ista)',impamp,'full');
% 
% end  


%% generate synthetics by stacking multiple green's arrivals at different times
moment = zeros(writes(end), 1);
for n = 1: writes(end)
  %Amplitude.  So far just using uniform for this subroutine.
  if beta > 100
    momomin = 1;
  else
    momomin = (1-rand(1)).^(-1/betam1); %a different random number if moments vary.  Inefficient?
    if momomin > maxmom
      n;
      maxmom = momomin;
    end
  end
  moment(n,1) = momomin;
  
  %Where it goes, +/- a shift
  %this just add templates arrived at sta 1
  ind = indtarvl(n);  % arrival time index
  greenst(n,1)=ind;  % although this index is returned in 'sources', but i want to make the correspondance clear
  if ind <= 0   %if arrival index is smaller than 1
    trunc = -ind;
    synth(1: Greenlenm1-trunc, 1) = synth(1: Greenlenm1-trunc, 1) + ...
      Greens(2+trunc: end, 1)*momomin;
  elseif ind+Greenlenm1 > size(synth,1) %if arrival index is larger than length
    trunc = ind+Greenlenm1 - size(synth,1);
    synth(ind: ind+Greenlenm1-trunc, 1) = ...
      synth(ind: ind+Greenlenm1-trunc, 1) + ...
      Greens(end-trunc, 1)*momomin;
  else
    synth(ind: ind+Greenlenm1, 1) = synth(ind: ind+Greenlenm1, 1) + Greens(:,1)*momomin; %Greens(whichG,:) %(for interpolated)

  end
%   if nspots == 1 %that is, if everything comes from the same location...
%     for ista=2:nsta
%       %if from same spot, there is no arrival time offset between stations
%       synth(ista,ind(n):ind(n)+Greenlenm1)=synth(ista,ind(n):ind(n)+Greenlenm1)+Greens(ista,:)*momomin; %Greens(whichG,:) %(for interpolated)
%     end
%   else
    for ista = 2: nsta %This assumes nsta=3
      %ind - rnoff is the arrival time in index at sta 2 and 3
      %note the sign here, if off12 >0, move 2 to the right to align with 1, meaning that 2 is
      %earlier than 1, ie., tarvl2 < tarvl1. Be consistent all the time
      if ind-rnoff(n,ista-1) <= 0   %if arrival index is smaller than 1
        trunc = -(ind-rnoff(n,ista-1));
        synth(1: Greenlenm1-trunc, ista) = synth(1: Greenlenm1-trunc, ista) + ...
          Greens(2+trunc: end, ista)*momomin;
      elseif ind-rnoff(n,ista-1)+Greenlenm1 > size(synth,1) %if arrival index is larger than length
        trunc = ind-rnoff(n,ista-1)+Greenlenm1 - size(synth,1);
        synth(ind-rnoff(n,ista-1): ind-rnoff(n,ista-1)+Greenlenm1-trunc, ista) = ...
          synth(ind-rnoff(n,ista-1): ind-rnoff(n,ista-1)+Greenlenm1-trunc, ista) + ...
          Greens(end-trunc, ista)*momomin;
      else
        synth(ind-rnoff(n,ista-1): ind-rnoff(n,ista-1)+Greenlenm1, ista) = ...
          synth(ind-rnoff(n,ista-1): ind-rnoff(n,ista-1)+Greenlenm1, ista) + ...
          Greens(:, ista)*momomin;
      end
    end
%   end

  %If time to write it, add elsewhere events first. Assumes momomin=1.
  if ismember(n,writes)
    for i = 1: floor(fracelsew*n)
      for ista = 1: nsta
        indelse = ceil(rnelsew(i,ista)*winlen);
        greenstelse(i,1) = indelse;
        %note that below just treat noise elsewhere as the same-amp template arrive at the
        %same time at all stations as at sta 1
        if indelse+Greenlenm1 > size(synth,1) %if arrival index is larger than length
          trunc = indelse+Greenlenm1 - size(synth,1);
          synth(indelse: indelse+Greenlenm1-trunc, ista) = ...
            synth(indelse: indelse+Greenlenm1-trunc, ista) + ...
            Greens(end-trunc, ista);
        else
          synth(indelse: indelse+Greenlenm1, ista) = ...
            synth(indelse: indelse+Greenlenm1, ista) + Greens(:, ista);%*momomin;
        end
      end
    end
    inwrites = inwrites+1
    n
    for ista = 1: nsta
      synths(:, ista, inwrites) = synth(skiplen: winlen, ista);
    end
    greensts{inwrites}{1}=greenst(1:n,1); % Chao, 2022/01/10, the starting index of each added template, context of full length
    if fracelsew~=0
      greensts{inwrites}{2}=greenstelse(1:floor(fracelsew*n), 1);
    end
    mommax(inwrites) = maxmom;
  end
end

sources=[sources moment]; % adding the colomn of the moment of each source


% keyboard





