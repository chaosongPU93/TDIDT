function [synths,mommax,sources,greensts]=csplaw3d(writes,winlen,skiplen,synth0,Greens,b,...
  xygrid,sps,fracelsew,seed,tdura,timetype,ftrans,stas)
%   From plaw3b but this one shifts but up to a certain distance, with offsets taken from John Armbruster's grid.
%    Format is offPGSS, offPGSI, dx, dy.
%    That's at 40 sps, so multiply by 2.5 to get offsets in 100 sps.
% The newer version uses plaw3d.m; the difference between that and plaw3c.m, used by synthshift.m,
% is in diff_plaw.  The structure of plaw has changed because of the need to remove 2 events from
% the same spot separated by less than the template duration.  So instead of just adding events as
% the seismograms become more saturated, for each saturation I needed to start making the catalog
% from scratch (zero events).
%
% --The orginal version written by Allan assumes a uniform, random distribution in arrival
% time, this may not the most ideal case if you want to obtain the synthetics at 4th stas.
% Instead, we can assume the same distribution for the origin time. With tools like Taup or
% hypoinverse, travel time between a specific source and station can be predicted with a
% velocity model. Therefore, arrival time can be obtained by adding travel time to origin
% time at any station theoretically, including 4th stations. Thus the synthetics at 4th sta
% can be generated using the same way by adding template at that sta to each arrival time.
% --It should be expected that the bulk of the code is similar to 'synthgen', except that
% this version of 'plaw3?' forces the separation of arrival from the same spot to be larger
% than 1 template duration.
% --2023/06/14 I changed his sign convention of 'ioff' to be the same as my routines.
%
%
% Return:
%   synths: synthetics, signal length * nsta * n sat level
%   mommax: max moment of sources
%   sources: source info, seismogram arrival index (1); source loc index (2)，at diff. sat
%       level, if instead set origin time as random, then tori and ttrvl are added
%   greensts: starting indices of added greens functions
%
% Coded bY Allan Rubin, modified by Chao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('timetype','tarvl');
defval('ftrans','interpchao');
defval('stas',[]);

beta=1.+b/1.5;
betam1=beta-1;
Greenlenm1=size(Greens,1)-1;
nsta=size(Greens,2);
rng(seed);
maxmom=1.0;
synths=zeros(winlen-skiplen+1, nsta,length(writes));
greensts = cell(length(writes),2);
% xygrid(:,1:2)=round(2.5*xygrid(:,1:2)); % *2.5 to get to 100 sps from 40.
xygrid(:,1:2)=round(sps/40*xygrid(:,1:2)); % *4 to get to 160 sps from 40.
nspots=size(xygrid,1)
%use 'off2space002' to interpolate for the travel time to sta 1
%loc has 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
[locuse, ~] = off2space002(xygrid(:,1:2),sps,ftrans,0);

if strcmp(timetype,'tarvl')
  sources=zeros(writes(end),2,length(writes));  %seismogram index (1); source location index (2); sat level (3)
elseif strcmp(timetype,'tori')
  sources=zeros(writes(end),4,length(writes)); % 2 more cols for tori and ttrvl
end
rnelsew=rand(floor(fracelsew*writes(end)),nsta);
for n0=1:length(writes)
  synth = synth0; %for each saturation level you NEED to start from scratch
  rn=rand(writes(n0),2); %uniform randoms in [0,1], random numbers are seismogram index, fault patch index
  sizern=size(rn,1)
  ispot=round(nspots*rn(:,2)+0.499999999); %ispot has length rn.
  if ispot == 0
    ispot = 1;  %just in case
  end
  rn(:,2)=ispot;
  ind=ceil(rn(:,1)*winlen);
  rn(:,1)=ind; %rn is now full of integers - seismogram index (1) and source location index (2)
  if strcmp(timetype,'tori')
    %corresponding travel time according to src loc indices
    ttrvl=locuse(rn(:,2),6);
    indttrvl = round(ttrvl*sps);  % in samples
    
    % %corresponding offset 12 and 13 according to src loc indices
    % rnoff=xygrid(rn(:,2),1:2);
    
    % %use 'off2space002' to interpolate for the travel time to sta 1
    % [locuse, ~] = off2space002(rnoff,sps,ftrans,0);
    % ttrvl = locuse(:,6);
    % indttrvl = round(ttrvl*sps);  % in samples
    
    %the arrival time for each source at sta 1 is: origin time + travel time
    indtori = rn(:,1);
    indtarvl = indtori + indttrvl;
    rn(:,1) = indtarvl;
    rn(:,3:4) = [indtori indttrvl];
  end
  if strcmp(timetype,'tarvl')
    rnsor=sortrows(rn,[2,1]);
  elseif strcmp(timetype,'tori')
    rnsor=sortrows(rn,[2,3]);
  end
  indx=1;
  thisspot=rnsor(indx,2);
  while (thisspot <= nspots) && (indx <= sizern-1)
    while isequal (rnsor(indx+1,2),thisspot) && (indx <= sizern-2) %-2 is something of a kludge.  Last entry might not be distinct.
      indx=indx+1;
      
      if strcmp(timetype,'tarvl')
        if rnsor(indx,1)-rnsor(indx-1,1) < tdura*sps %if separated by < 1 template duration
          %rnsor(indx,2)=round(nspots*rand+0.499999999);
          rnsor(indx,1)=rnsor(indx-1,1)+tdura*sps; %This could conceivably extend seismogram beyond end, esp. at high saturation
        end
      elseif strcmp(timetype,'tori')
        if rnsor(indx,3)-rnsor(indx-1,3) < tdura*sps %if separated by < 1 template duration
          %rnsor(indx,2)=round(nspots*rand+0.499999999);
          rnsor(indx,3)=rnsor(indx-1,3)+tdura*sps; %This could conceivably extend seismogram beyond end, esp. at high saturation
          rnsor(indx,1)=rnsor(indx-1,1)+tdura*sps;
        end
      end
    end
    thisspot=thisspot+1;
    indx=indx+1;
  end
  %limit arrivals to be within signal length
  rnsor=sortrows(rnsor,1);
  idis = find(rnsor(:,1)>winlen);
  rnsor(rnsor(:,1)>winlen,:)=[];
  sizernsor=size(rnsor,1)
  ind=rnsor(:,1); %ind is now back to the original definition (seismogram index), but in order given by rnsor
  %     indsor=sortrows(rnsor,1); %Just checking...
  %     indsor(end-10:end,:)
  %     Greenlenm1
  %     synthsize=size(synth)
  %     winlen
  %2023/06/14 I changed his sign convention of 'ioff' to be the same as my routines.
  ioff=xygrid(rnsor(:,2),1:2); %1 is PGSS; 2 is PGSI %THIS MINUS SIGN IS EMPIRICAL!!!
  % rnsor(1:20,:)
  % ioff(1:20,:)
  %info of synthetic sources
  sources(1:sizernsor,:,n0)=rnsor; % 2 or 4 cols depending on timetype
  
  % keyboard
  %if also wants synthetics at 4th stas
  if nsta > 3
    for ista = 4: nsta
      %predict 'off14' at a 4th station using the plane fit model
      [~,off14] = pred_tarvl_at4thsta(stas(ista,:),ioff(:,1),ioff(:,2),ind);
      ioff = [ioff off14];
    end
  end
  
  greenst = ind;
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
    if ind(n) <= 0   %if arrival index is smaller than 1
      trunc = -ind(n);
      synth(1: Greenlenm1-trunc, 1) = synth(1: Greenlenm1-trunc, 1) + ...
        Greens(2+trunc: end, 1)*momomin;
    elseif ind(n)+Greenlenm1 > size(synth,1) %if arrival index is larger than length
      trunc = ind(n)+Greenlenm1 - size(synth,1);
      synth(ind(n): ind(n)+Greenlenm1-trunc, 1) = ...
        synth(ind(n): ind(n)+Greenlenm1-trunc, 1) + ...
        Greens(end-trunc, 1)*momomin;
    else
      synth(ind(n): ind(n)+Greenlenm1, 1) = synth(ind(n): ind(n)+Greenlenm1, 1) + Greens(:,1)*momomin; %Greens(whichG,:) %(for interpolated)
    end
    
    for ista=2:nsta %now nsta can be more than 3, as 'ioff' has cols for 4th stas too
      if ind(n)-ioff(n,ista-1) <= 0
        trunc=-(ind(n)-ioff(n,ista-1));
        synth(1:Greenlenm1-trunc, ista)=synth(1:Greenlenm1-trunc, ista)+...
          Greens(2+trunc:end, ista)*momomin;
      elseif ind(n)-ioff(n,ista-1)+Greenlenm1 > size(synth,1)
        trunc=ind(n)-ioff(n,ista-1)+Greenlenm1 - size(synth,1);
        synth(ind(n)-ioff(n,ista-1):ind(n)-ioff(n,ista-1)+Greenlenm1-trunc, ista)= ...
          synth(ind(n)-ioff(n,ista-1):ind(n)-ioff(n,ista-1)+Greenlenm1-trunc, ista)+...
          Greens(end-trunc, ista)*momomin;
      else
        synth(ind(n)-ioff(n,ista-1):ind(n)-ioff(n,ista-1)+Greenlenm1, ista)=...
          synth(ind(n)-ioff(n,ista-1):ind(n)-ioff(n,ista-1)+Greenlenm1, ista)+...
          Greens(:, ista)*momomin;
      end
    end
    
    %     f2.fig = figure(100); clf(f2.fig);
    %     f2.fig.Renderer = 'painters';
    %     widin = 12;  % maximum width allowed is 8.5 inches
    %     htin = 3;   % maximum height allowed is 11 inches
    %     [scrsz, resol] = pixelperinch(1);
    %     set(f2.fig,'Position',[2*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);
    %     nrow = 1;
    %     ncol = 1;
    %     for isub = 1:nrow*ncol
    %       f2.ax(isub) = subplot(nrow,ncol,isub);
    %     end
    %     ax=f2.ax(1);
    %     tmp = synth(:,1)-noi(:,1);
    %     plot(ax,tmp,'k'); hold on
    %     scatter(ax,ind(n)+769,0,20,'b');
    %     xlim(ax,[0+skiplen 31*sps+skiplen]);
    %
    %     keyboard
    
  end
  %add elsewhere events. Assumes momomin=1.
  for i=1:floor(fracelsew*writes(n0))
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
    greensts{n0}{2}=greenstelse(1:floor(fracelsew*writes(n0)), 1);
  end
  mommax(n0)=maxmom;
end

%     f2.fig = figure(200); clf(f2.fig);
%     f2.fig.Renderer = 'painters';
%     widin = 12;  % maximum width allowed is 8.5 inches
%     htin = 3;   % maximum height allowed is 11 inches
%     [scrsz, resol] = pixelperinch(1);
%     set(f2.fig,'Position',[2*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);
%     nrow = 1;
%     ncol = 1;
%     for isub = 1:nrow*ncol
%       f2.ax(isub) = subplot(nrow,ncol,isub);
%     end
%     ax=f2.ax(1);
%     tmp = synths(:,1,2)-noi(skiplen:winlen,1);
%     plot(ax,tmp,'k'); hold on
%     scatter(ax,greensts{2}{1}+769-skiplen,zeros(size(greensts{2}{1})),20,'b');
%     xlim(ax,[0 31*sps]);
%
% keyboard