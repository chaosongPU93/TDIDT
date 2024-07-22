function srcvalidate(insat,sps,zcrosses,writes,distrloc,sources,xygrid,...
  greensts,stas,skiplen,greenlen,synths,synth,winlen,green)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the function extract and validate the added impulses of template.
% The impulses used to generate the synthetics in the algorithm was just to keep
% adding a template at a specific time with a scaled amplitude. Technically,
% this is equivalent to a one-time convolution between the full source time 
% function (a series of impulses with arrival time and amp) and the LFE
% templates actually used in generation. 
% Given two ways should be the same, this code therefore can uses convolution
% to test if the generation was done as planned.
%
%
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/07/20
% Last modified date:   2024/07/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmpgrid = xygrid;
tmpgrid(:,1:2)=round(sps/40*tmpgrid(:,1:2)); % *4 to get to 160 sps from 40.

nsta = length(zcrosses);

wlensec = 30; %how long the window to test
bufsec = 1; %need some buffer window for later CC alignment
buffer = bufsec*sps;
wlensecb = wlensec+bufsec;

sigstab = [];  %target window of synthetic signals at all stations
sigpnstab = [];  %target window of synthetic signals plus the noise
sigconvstab = [];  %target window of reproduced synthetic signals using convolution
noistab = [];  %target window of noise

srcpairstab = cell(nsta,1);  %sources whose zero-crossing with time range, [indtarvl rnoff12 rnoff13 amp]

for ista = 1: nsta
  % ista = 3;
  tgreen = zcrosses(ista);  % anchor time of templates, here choose as zero-crossing time
  indstsig = 1;  % starting index of the simulated signal to test
  indedsig = wlensecb*sps+indstsig-1; % ending index of the simulated signal to test
  %   source = squeeze(sources(1:size(greensts{insat}{1},1), :, insat));  % [indtori indttrvl indtarvl rnoff12 rnoff13 amp];
  n=writes(insat);
  if strcmp(distrloc,'uniform')
    %     if forcesep
    a = squeeze(sources(1:n,:,insat));
    b = a(any(a,2),:);
    %     else
    %       a = sources(1:n,:);
    %       b = a(any(a,2),:);
    %     end
    source=b;
    off = tmpgrid(source(:,2),1:2); %note that 'tmpgrid' has the desired sps
  elseif strcmp(distrloc,'custompdf')
    a = sources(1:n,:);
    b = a(any(a,2),:);
    source=b;
    off = source(:,1:2);
  end
  if ista == 1
    greenst = greensts{insat}{1}; % the starting index of each added template, context of full length
  elseif ista <=3
    %ind - rnoff is the arrival time in index at sta 2 and 3
    %note the sign here, if off12 >0, move 2 to the right to align with 1, meaning that 2 is
    %earlier than 1, ie., tarvl2 < tarvl1. Be consistent all the time
    greenst = greensts{insat}{1}-off(:, ista-1);
  else
    greenst = pred_tarvl_at4thsta(stas(ista,:),off(:,1),off(:,2),greensts{insat}{1});
  end
  
  %%%you don't need all impulses, only some of them contribute to the length of truncated record
  %%%cut out the green indice and sources that contribute
  %%%Note that some zero-crossings might appear later than 'indedsig', but the corresponding start
  %%%index of the green is actually smaller than 'indedsig', ie. their 'greenst' <= 'indedsig'
  induse = find(greenst>=skiplen-greenlen+indstsig & greenst<=indedsig+skiplen-1);
  greenst = greenst(induse);
  source = source(induse,:);
  off = off(induse,:);
  source = [source(:,1) off];
  
  greenzc = greenst+tgreen; % index of approximate zero-crossing
  source(:, 1) = source(:, 1)+zcrosses(1);  % now 'greenzc' should be the same as 1st col of 'source'
  impamp = zeros(max(greenzc)+20,1);
  for i = 1: length(greenzc)
    impamp(greenzc(i)) = impamp(greenzc(i))+1;  %assuming every source has an amp of 1
  end
  
  %note that some length of the full simulation 'skiplen' was skipped in 'synthgen.m'
  sigpntemp = squeeze(synths(:,ista,insat));  % simulated signal with white noise
  noitemp = synth(skiplen:winlen,ista);   % account for the skipped part
  sigtemp = sigpntemp-noitemp; % subtract the white noise
  sigpnstab(:,ista) = sigpntemp(indstsig:indedsig); % now focus on the part of interest
  sigstab(:,ista) = sigtemp(indstsig:indedsig);
  noistab(:,ista) = noitemp(indstsig:indedsig);
  
  sigconvtmp = conv(green(:,ista),impamp,'full');
  indtrc = tgreen+skiplen;  % starting index for truncation
  sigconvstab(:,ista) = sigconvtmp(indtrc+indstsig-1:indedsig+indtrc-1);  % cut accordingly
  
  % figure
  % plot(sig,'k'); hold on
  % plot(sigconv+3,'b');
  % plot(sigconv-sig+1,'r-')
  
  %%%Can we reproduce the synthetics with truncated convolution? YES
  figure
  subplot(411)
  plot(green(:,ista),'r');
  xlim([0 greenlen]);
  legend('Template (unfiltered)');
  title('Reproduce synthetics with truncated convolution');
  subplot(412)
  imptemp = find(impamp>0);
  p1=stem(imptemp-skiplen, impamp(imptemp),'b','MarkerSize',4); hold on;
  % p1=stem((1:size(impamp))-skiplen, impamp,'b','MarkerSize',4); hold on;
  ax = gca;
  plot([indstsig indstsig],ax.YLim,'--','color',[.5 .5 .5]);
  plot([indedsig indedsig],ax.YLim,'--','color',[.5 .5 .5]);
  legend(p1,'Synthetic random impulses');
  xran1 = ax.XLim; axpos1 = ax.Position(1);
  subplot(413);
  plot(sigstab(:,ista),'b'); hold on
  ax=gca; xran2 = ax.XLim;
  shrink(ax,xran1/xran2,1);
  ax.Position(1)=axpos1;
  plot(sigconvstab(:,ista),'k');
  legend('Truncated synthetic signal','Truncated signal from convolution');
  subplot(414);
  plot(sigstab(:,ista)-sigconvstab(:,ista),'k'); hold on
  legend('Difference');
  ax=gca; xran3 = ax.XLim;
  shrink(ax,xran1/xran3,1);
  ax.Position(1)=axpos1;
  %   keyboard
  
  %now 'source' contain impulse triplets that contribute to the trace segment of interest, but we
  %want to get ones whose zero-crossing falls into the window.
  source(:,1) = source(:,1)-skiplen;  % index after skipping
  srcpairb = source(source(:,1)>=indstsig & source(:,1)<=indedsig, :);
  srcpairb = sortrows(srcpairb,1,'ascend');
  srcpairstab{ista} = srcpairb; %ideally for each station this should be the same, but coincidence is possible
  
end

%notify if sources are the same at diff stations
if ~isequaln(srcpairstab{1,1},srcpairstab{2,1}) || ~isequaln(srcpairstab{1,1},srcpairstab{3,1})
  disp('Extracted sources at different stations are not the same, re-check!');
end

