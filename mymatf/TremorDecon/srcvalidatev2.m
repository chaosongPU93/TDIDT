function srcvalidatev2(sps,zcrosses,synsrc,synsrcstind,stas,skiplen,greenlen,...
  STAopt,synth,winlen,green)
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

nsta = length(zcrosses);

wlensec = 30; %how long the window to test
bufsec = 1; %need some buffer window for later CC alignment
buffer = bufsec*sps;
wlensecb = wlensec+bufsec;

sigstagt = [];  %target window of synthetic signals at all stations
sigpnstagt = [];  %target window of synthetic signals plus the noise
sigconvstagt = [];  %target window of reproduced synthetic signals using convolution
noistagt = [];  %target window of noise
synsrcgtsta = cell(nsta,1);  %sources whose zero-crossing with time range, [indtarvl rnoff12 rnoff13 amp]
for ista = 1: nsta
  tgreen = zcrosses(ista);  % anchor time of templates, here choose as zero-crossing time
  indst = 1;  % starting index of the simulated signal to test
  inded = wlensecb*sps+indst-1; % ending index of the simulated signal to test
  source = synsrc;
  if ista == 1
    greenst = synsrcstind; % the starting index of each added template, context of full length
  elseif ista <=3
    %ind - rnoff is the arrival time in index at sta 2 and 3
    %note the sign here, if off12 >0, move 2 to the right to align with 1, meaning that 2 is
    %earlier than 1, ie., tarvl2 < tarvl1. Be consistent all the time
    greenst = synsrcstind-source(:, ista);
  else
    greenst = pred_tarvl_at4thsta(stas(ista,:),source(:,2),source(:,3),source(:,1));
  end
  %%%you don't need all impulses, only some of them contribute to the length of truncated record
  %%%cut out the green indice and sources that contribute
  %%%Note that some zero-crossings might appear later than 'indedsig', but the corresponding start
  %%%index of the green is actually smaller than 'indedsig', ie. their 'greenst' <= 'indedsig'
  induse = find(greenst>=skiplen-greenlen+indst & greenst<=inded+skiplen-1);
  greenst = greenst(induse);
  source = source(induse,:);
  %     keyboard
  
  greenzc = greenst+tgreen; % index of approximate zero-crossing
  source(:, 1) = source(:, 1)+zcrosses(1);  % now 'greenzc' should be the same as 1st col of 'source'
  impamp = zeros(max(greenzc)+20,1);  %20 is a bit arbitary
  for i = 1: length(greenzc)
    impamp(greenzc(i)) = impamp(greenzc(i))+source(i, 6);
  end
  
  %note that some length of the full simulation 'skiplen' was skipped in 'synthgen.m'
  sigpntemp = STAopt(:,ista);  % simulated signal with white noise
  noitemp = synth(skiplen:winlen,ista);   % account for the skipped part
  sigtemp = sigpntemp-noitemp; % subtract the white noise
  sigpnstagt(:,ista) = sigpntemp(indst:inded); % now focus on the part of interest
  sigstagt(:,ista) = sigtemp(indst:inded);
  noistagt(:,ista) = noitemp(indst:inded);
  
  %direct convolution
  sigconvtmp = conv(green(:,ista),impamp,'full');
  indtrc = tgreen+skiplen;  % starting index for truncation
  sigconvstagt(:,ista) = sigconvtmp(indtrc+indst-1:inded+indtrc-1);  % cut accordingly
  
  %%%Can we reproduce the synthetics with truncated convolution? --YES
  figure
  subplot(411)
  plot(green(:,ista),'r');
  xlim([0 greenlen]);
  legend('Template (unfiltered)');
  text(0.05,0.9,stas(ista,:),'HorizontalAlignment','left','Units','normalized');
  title('Reproduce synthetics with truncated convolution');
  subplot(412)
  imptemp = find(impamp>0);
  p1=stem(imptemp-skiplen, impamp(imptemp),'b','MarkerSize',4); hold on;
  % p1=stem((1:size(impamp))-skiplen, impamp,'b','MarkerSize',4); hold on;
  ax = gca;
  plot([indst indst],ax.YLim,'--','color',[.5 .5 .5]);
  plot([inded inded],ax.YLim,'--','color',[.5 .5 .5]);
  legend(p1,'Synthetic random impulses');
  xran1 = ax.XLim; axpos1 = ax.Position(1);
  subplot(413);
  plot(sigstagt(:,ista),'b'); hold on
  ax=gca; yran = ax.YLim; xran2 = ax.XLim;
  shrink(ax,xran1/xran2,1);
  ax.Position(1)=axpos1;
  plot(sigconvstagt(:,ista),'k');
  legend('Truncated synthetic signal','Truncated signal from convolution');
  subplot(414);
  plot(sigstagt(:,ista)-sigconvstagt(:,ista),'k'); hold on
  legend('Difference'); ylim(yran)
  
  %now 'source' contain impulse triplets that contribute to the trace segment of interest, but we
  %want to get ones whose zero-crossing falls into the window.
  source(:,1) = source(:,1)-skiplen;  % index after skipping
  tmp = source(source(:,1)>=indst & source(:,1)<=inded, :);
  tmp = sortrows(tmp,1,'ascend');
  synsrcgtsta{ista} = tmp; %ideally for each station this should be the same, but coincidence is possible
  %     keyboard
end
%notify if sources are the same at diff stations
if ~isequaln(synsrcgtsta{1,1},synsrcgtsta{2,1}) || ~isequaln(synsrcgtsta{1,1},synsrcgtsta{3,1})
  disp('Extracted sources at different stations are not the same, re-check!');
end
%   keyboard

