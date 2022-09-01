function [sigdecon,pred,res,dresit,mfitit,ampit,nit,fighdl] = ...
  iterdecon_joint(sig,wlet,noi,dt,twlet,loff_max,dres_min,mfit_min,nit_max,...
  nimp_max,fpltit,fpltend,fpltchk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [sigdecon,pred,res,dresit,mfitit,nit,rf] = ...
%   iterdecon_joint(sig,wlet,noi,dt,twlet,loff_max,dres_min,mfit_min,nit_max,...
%   nimp_max,fpltit,fpltend,fpltchk)
%
% Similar to 'iterdecon', this function is to carry out the iterative deconvolution after 
% Ligorria & Ammon (1999), Essentially, wlet will be deconvolved from sig. The difference
% is to do the joint deconvolution to several stations simultaneously. This means, after
% the signals at different stations are prealigned, presuming now the events to be 
% deconvolved are restricted in space (ideally from the same spot), we allow a small shift
% in samples between stations when deconvolving the 'same' corresponding event during each
% iteration. Therefore, the final impulses from all stations are self-consistent and easier
% to compare with each other.
%
% INPUT:
%   sig: signal / enumerator (impulse response * source wavelet estimation)
%   wlet: The source wavelet estimation / denominator / template
%   dt: sampling interval [s]
%   twlet: Time shift of main arrival of wavelet [s]. The default is 0.
%   width:  Gaussian width parameter. The default is 2.5, e.g., 1.5
%   dres_min:  Convergence control parameter (percentage change in residual
%               per iteration). The default is 0.5.
%   mfit_min:  Convergence control parameter (norm of the residual/misfit).
%               The default is norm(sig). NOT necessary
%   nit_max: Maximal number of iterations before the algorithm interrupts,
%               defaults to 400
%   fpltit: flag for plotting (1) or not (0) for each iteration
%   fpltend: flag for plotting (1) or not (0) for the final iteration
%   fcheck: flag for plotting (1) or not (0) for intermediate computations
%
% OUTPUT:
%   sigdecon: deconvolved impulses, ie., the estimation of impulse responses
%   pred: predicted signal
%   res:  residual, aka signal misfit
%   dresit: percentage of relative change in residual of each iteration
%   mfitit: norm of the residual/misfit of each iteration
%   ampit:  amplitude of decolved impulse of each iteration
%   nit:  actual number of iterations until algorithm converges/stops
%   fighdl: handles of figures plotted
%   rf: receiver function if needed, MAY BE ABANDONED LATER
%
%
% NOTE:
% --Technically, this function is only necessary if the independently-determined 
%   impulses from different stations are not corresponding to each other, i.e. some
%   impulses are too far offseted wrt others.
% --Need to specify the maximum offsets allowed between impulses in each iteration.
%   Given that target signals are prealigned and a high running CC between them, 
%   all events inside should come from a compact area.
% --Evaluate the rcc value at the prominent peaks in master CCs at stations 
%   together. Rank them, if the same-rank peaks are deviated more than max allowed
%   shifts, go to the next rank, compare the cumulative peaks. The number of 
%   combinations of peaks to consider is (n_rank)^n_sta
% --Compute the running CC also iteratively for each iteration, based on the updated
%   residuals at each station.
% --There are multiple choices of finding peaks, you could find all peaks at each
%   stations, or find peaks that are higher than the max amp of the CC between noise
%   and template. But I guess if use the former, you'd better stop the iterations
%   when the deconvloved impulse amplitude at all station are lower than the max amp 
%   of the CC between noise and template. The difference is that, add constraint in 
%   finding peaks exclude small peaks at all staions always, but the stopping criterion
%   can allow some staions (as long as not all of them) has a smaller peak than noise
%   level.
% --There are also different ways of find the 'best' pair of peaks. You can either find
%   the one where the offsets between stations are smallest and choose the one with the
%   highest weighted CC if multiples are closest. Or you can first find all pairs
%   close enough (offset smaller than allowed), then choose the one with the highest 
%   weighted CC. The two emphasize different things, the former wants the events to be
%   as compact in space as possible, while the latter wants the events to be as coherent
%   and 'likely' as possbile. Under some occasions, however, the two can lead to the same
%   results.
% --Again, the number of iterations (nit) >= number of non-zero impulses (nimp) 
%
%
% Things to be added:
% --How to deal with the case when the signal is also a single arrival dominated, but narrower in
%   width compared to the wavelet, how to reflect that in the impulses
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/01/24
% Last modified date:   2022/01/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%% default value for easy debugging or no input
% defval('shift','0');
% defval('width','2.5');
% defval('dfit_min','0.5');
% defval('nit_max','400');
if isempty(twlet)
  twlet = zeros(1,size(sig,2));
end
if isempty(dres_min)
  dres_min = 0.05;
end
if isempty(mfit_min)
  mfit_min = norm(median(sig,2));
end
if isempty(nit_max)
  nit_max = 400;
end
if isempty(nimp_max)
  nimp_max = 400;
end
if isempty(noi)
  noi = zeros(size(sig));
end
if isempty(fpltit)
  fpltit = 1;
end
if isempty(fpltend)
  fpltend = 1;
end
if isempty(fpltchk)
  fpltchk = 0;
end

dres_min = dres_min/100;  % from percent to decimal
itwlet = round(twlet/dt); % time shift of main arrival of wavelet in samples

fighdl = cell(4,1);


%% Preparation, iteration 0
nsta = size(sig, 2); % Note that this means each col is the signal at one station  
lsig = size(sig, 1);   % length of signal
lwlet = size(wlet, 1);  % length of wavelet, usually shorter than signal
nfft = lsig;  % no padding needed

if lwlet>lsig   % if the wavelet is longer than signal, use way 1 to do CC
  CCOPT = 1;
else
  CCOPT = 2; % otherwise, use way 2 to do CC for efficiency
end

%%% iteration 0
nit = 0;  % iteration number
res = sig;  % residual
pred = zeros(nfft,nsta);  % predefine the prediction array
sigdecon = zeros(nfft,nsta);  % predefine the response array
dres = ones(1,nsta);   % decrease of residual
mfit = max(1, norm(sig));   % misfit between data and prediction
ampit = [];   % amplitude of impulse relative to the max amp of the template
dresit = [];  % relative change in residual for each iteration
mfitit = [];  % misfit for each iteration, norm of the residual 
nimp = 0; % total number of impulses 
predchg = sig;   % change in prediction 
amp = ones(1,nsta);   % intial amplitude of impulse
mwtcoef = 1;  % intial amplitude of max master CC weighted by rcc

for i = 1: nsta
  [coefnoi, lagnoi] = xcorr(noi(:,i), wlet(:,i), size(noi, 1), 'none'); % unnormalized cc
  coefnoieff = coefnoi(lagnoi+itwlet(i) >= 1 & lagnoi+itwlet(i) <= size(noi, 1));
  [pkhgt, pkind] = findpeaks(coefnoieff);
  ampnoi(i) = median(pkhgt);
%   [ampnoi(i), ~] = max(coefnoi);   % a measure of noise level
%   figure
%   plot(coefnoi); hold on
%   plot(coefnoieff);
%   scatter(pkind,pkhgt,20,'k');
end
%   noiwletrat = rms(noi)./max(wlet,[],1);  % ratio of amplitude of noise relative to the max amp of wavelet
ampnoin = ampnoi./sum(wlet.^2, 1);   % normalized amp of CC between noise and wavelet 

%average running CC between signals 
cclen=round(1/dt)/2;
[irccs,rccs12] = RunningCC(sig(:,1), sig(:,2), cclen);
[~,rccs13] = RunningCC(sig(:,1), sig(:,3), cclen);
[~,rccs23] = RunningCC(sig(:,2), sig(:,3), cclen);
rccs = (rccs12+rccs13+rccs23)/3;

ampnoiwt = ampnoi.*median(rccs); % amp of CC between noise and wavelet weighted by median rcc
msumwtcoef = 1;

%% iteration 
%%%%%%%%%%%%%%% Different while loop conditions for iterative deconvolution
%%%Some of them are connected tp each other, e.g., The change in residual depends on the prediction
%%%change which is determined by the convolution of the amplitude of the deconvolved impulse
%%%and the template. And it is a open question on how to quantify the noise and residual level, etc.
% while sum(mad(res)>mad(noi))>0 && nimp<nimp_max
% while sum(max(abs(predchg))>mad(noi))>0 && nimp<nimp_max
% while nit<nit_max
% while sum(dres>dres_min)>0 && nit<nit_max
% while sum(dres>dres_min)+sum(mfit>mfit_min)>0 && nit<nit_max
% while sum(mad(res)>mad(noi))>0 && nit<nit_max % compare the amp of residual with amp of noise 
% while sum(max(abs(predchg))>mad(noi))>0 && nit<nit_max % compare the max. amp of change in prediction with amp of noise 
while sum(amp>ampnoin)>0 && nit<nit_max % compare the newest max CC of res and wlet with CC of noise and wlet, all normalized by wlet 
% while msumwtcoef>sum(ampnoiwt) && nit<nit_max % compare the newest max CC of res and wlet weighted by rcc, with median CC of noise and wlet, weighted bt median rcc 
  
  nit = nit+1;

  %compute the running CC between residuals at all stations, at iter 0, residuals equal to signals
  [ircc,rcc12] = RunningCC(res(:,1), res(:,2), cclen);
  [~,rcc13] = RunningCC(res(:,1), res(:,3), cclen);
  [~,rcc23] = RunningCC(res(:,2), res(:,3), cclen);
  rcc = (rcc12+rcc13+rcc23)/3;
  
  ampnoiwtn = ampnoi.*median(rcc)./sum(wlet.^2, 1); % normalized amp of CC between noise and wavelet weighted by median rcc
  ampnoiwt = ampnoi.*median(rcc); % amp of CC between noise and wavelet weighted by median rcc

  %%%%%%%%%%%%%% Different ways to do CC %%%%%%%%%%%%%%
  if CCOPT == 1
    %%% Way 1: redo CC for each iteration, START %%%
    for i = 1: nsta
      %cross-correlation between the residual and wavelet
      [coef(:,i), lag] = xcorr(res(:,i), wlet(:,i), nfft, 'none'); % raw CC
    end
    %%% Way 1, END %%%
  elseif CCOPT == 2
    %%% Way 2: do the CC once, then update it, efficient, START %%%
    for i = 1: nsta
      if nit == 1
        %cross-correlation between the signal and wavelet
        [coef(:,i), lag] = xcorr(sig(:,i), wlet(:,i), nfft, 'none'); % master raw CC
      else
        st = max(idximp(i)-itwlet(i)+1, 1);  % start index of predication change, in case < 1
        ed = min(idximp(i)-itwlet(i)+lwlet, nfft);  % end index of predication change, in case > nfft
        %cross-correlation between the predication change and wavelet
        [coefchg, ~] = xcorr(predchg(st: ed, i), wlet(:,i), lwlet, 'none'); % raw CC
        [~, ind] = max(coefchg);   % max of raw cc

        cnt = lagsamp(i)+nfft+1;  % center index of master raw cc to be updated
        st = 1+(cnt-ind);   % start index of master raw cc to be updated
        ed = 2*lwlet+1+(cnt-ind);   % end index of master raw cc to be updated
        if st < 1 && ed <= 2*nfft+1 % in case out of bounds
          st = 1;
          coefchg = coefchg(2*lwlet+1-(ed-st): 2*lwlet+1);  % drop a few points at the front
        elseif st >= 1 && ed > 2*nfft+1  % in case out of bounds
          ed = 2*nfft+1;
          coefchg = coefchg(1: (ed-st)+1); % drop a few points at the back
        end
        indchg = st: ed;  % indexes of master raw cc to be updated
        %update the master raw CC vector by accounting for the change
        coef(indchg, i) = coef(indchg, i)-coefchg; % note the neg. sign here
      end
    end
    %%% Way 2, END %%%
  end
  %%%%%%%%%%%%%% Different ways to do CC %%%%%%%%%%%%%%

  %%%%%%%%%%%%%% Weight master CC by running CC %%%%%%%%%%%%%%
  %%% Scheme 2: rank all peaks in master CC by the dot product between the rcc and master CC at each
  %%% peak, and choose the one with the highest dot product
  %%% this would perserve the true peaks in master CC without shifting, and no negative
  %%% amplitude would come out
  for i = 1: nsta
    %Effective raw CC that corresponds to the index of the overlapping portion between signal (or 
    %sigdecon) and running CC, same length as rcc
    coefeff(:, i) = coef(lag+itwlet(i) >= 1+round(cclen/2) & ...
                   lag+itwlet(i) <= nfft-round(cclen/2), i);          
    %%%Find all peaks in the effective master raw CC, return the peak value and index in effective  
    %%%raw CC, context of overlapping portion
    pks = [];
    %There will be 2 choices of finding peaks here, seems like both ways are sensible to some
    %extent, but note that limiting the height at each station is more strict, leading to fewer
    %impulses.
    [pks(:,2), pks(:,1)] = findpeaks(coefeff(:, i));
%     [pks(:,2), pks(:,1)] = findpeaks(coefeff(:, i),'MinPeakHeight',ampnoi(i));
    if isempty(pks)   % force stop if there are no more eligible peaks
      pkssta{nit}{i} = pks;
      continue  %
    else
      pks(:,3) = rcc(pks(:,1)) .* pks(:,2); % weight the master cc by rcc value at each peak 
      pkssort = sortrows(pks,3,'descend'); % sort the peak based on rcc value
      pkssta{nit}{i} = pkssort; % pks, [ind_of_peak, peak_height(master cc value), weighted cc by rcc]
    end
  end
  %This is a stupid way to deal with 3 stations only, but don't have to write explicitly like this
  pkss1 = pkssta{nit}{1}; npkss1 = size(pkss1,1);
  pkss2 = pkssta{nit}{2}; npkss2 = size(pkss2,1);
  pkss3 = pkssta{nit}{3}; npkss3 = size(pkss3,1);
  
  %This is useful only when you set limitations to the peaks, otherwise there will always be peaks
  if isempty(pkss1) || isempty(pkss2) || isempty(pkss3)
    nit = nit-1;  % -1 to recognize the current iteration does not proceed to the end
    fprintf('Forced stop at iteration %d \n', nit);
    fprintf('Because no eligible peaks are found (depending on the searching criterion) \n'); 
    fprintf('at the next iteration \n');
    break
  end
  
  %Start from the first pair of sorted peaks
  ipk = 1;
  off12 = pkss1(1,1)-pkss2(1,1);  % sign included
  off13 = pkss1(1,1)-pkss3(1,1);
  off23 = pkss2(1,1)-pkss3(1,1);
  off = (abs(off12)+abs(off13)+abs(off23))/3;   % note that this is the absolute average
  off12b = off12;
  off13b = off13;
  off23b = off23;
  offb = off;
  bsub = [1 1 1];
  
  while (abs(off12b)>loff_max || abs(off13b)>loff_max || abs(off23b)>loff_max) && ...
      ipk<max([npkss1 npkss2 npkss3]) 
    ipk = ipk+1;  % go down to the next peak
    off12 = nan(min(ipk,npkss1), min(ipk,npkss2), min(ipk,npkss3)); % you can't go more than availble peaks
    off13 = off12; 
    off23 = off12;
    for ii = 1: min(ipk,npkss1)  % at station 1, so each dimension can be different, but only the combination matters
      for jj = 1: min(ipk,npkss2)  % at station 2
        for kk = 1: min(ipk,npkss3)  % at station 3
          off12(ii,jj,kk) = pkss1(ii,1)-pkss2(jj,1);
          off13(ii,jj,kk) = pkss1(ii,1)-pkss3(kk,1);
          off23(ii,jj,kk) = pkss2(jj,1)-pkss3(kk,1);
        end
      end
    end
    off = (abs(off12)+abs(off13)+abs(off23))/3;   % note that this is the absolute average
    %%%%%%%%%%%% Different ways to determine the 'best' peak %%%%%%%%%%%%%%%%
    %%%% Under some circumstances, they could produce the same result, but ideas behind them
    %%%% are different. Way2 seems more logical than way1. 
%     %%% Way1: Sequential go down the sorted peaks to find the pair that is closet to each other 
%     indbest = find(off==min(off(:)));
%     [sub1, sub2, sub3] = ind2sub(size(off),indbest); % convert linear indices to matrix subscripts
%     if length(indbest)>1  % if there are multiple combinations whose offsets are equally close
%       sumwtcoef = pkss1(sub1,3)+pkss2(sub2,3)+pkss3(sub3,3); % obtain sum of the weighted coefs at these pairs
%       [~,ind] = max(sumwtcoef);   % choose the pair has the highest weighted coef sum
%       bsub = [sub1(ind), sub2(ind), sub3(ind)];
%     else
%       bsub = [sub1, sub2, sub3];
%     end
%     %%% Way1: Sequential go down the sorted peaks to find the pair that is closet to each other
    
    %%% Way2: find all pairs that are close enough, choose the one with highest weighted CC
    indclose = find(abs(off12)<=loff_max & abs(off13)<=loff_max & abs(off23)<=loff_max);
    if isempty(indclose)
      continue
    else
      [sub1, sub2, sub3] = ind2sub(size(off12),indclose); % convert linear indices to matrix subscripts
      sumwtcoef = pkss1(sub1,3)+pkss2(sub2,3)+pkss3(sub3,3); % obtain sum of the weighted coefs at these pairs
      [msumwtcoef,ind] = max(sumwtcoef);   % choose the pair has the highest weighted coef sum
%       sumind = pkss1(sub1,1)+pkss2(sub2,1)+pkss3(sub3,1); % obtain the sum of the indices
%       [msumind,ind] = min(sumind);   % choose the pair that is the earliest in time
      bsub = [sub1(ind), sub2(ind), sub3(ind)];
    end
    %%% Way2: find all pairs that are close enough, choose the one with highest weighted CC
    %%%%%%%%%%%% Different ways to determine the 'best' peak %%%%%%%%%%%%%%%%
    
    off12b = off12(bsub(1),bsub(2),bsub(3));  % again, offset between a station pair has the sign
    off13b = off13(bsub(1),bsub(2),bsub(3));
    off23b = off23(bsub(1),bsub(2),bsub(3));
    offb = off(bsub(1),bsub(2),bsub(3));  % again, average offset among all stations is an absolute value
  end
  
  %Stop if offset between stations are large than allowed after going through all availble peaks
  if abs(off12b)>loff_max || abs(off13b)>loff_max || abs(off23b)>loff_max || offb>loff_max
    nit = nit-1;  % -1 to recognize the current iteration does not proceed to the end
    fprintf('Forced stop at iteration %d \n', nit);
    fprintf('The offsets of best solution are still larger than allowed \n');
    fprintf('at the next iteration \n');
    fprintf('off12= %d, off13= %d, off23= %d \n',off12b, off13b, off23b);
    break
  else
    fprintf('off12= %d, off13= %d, off23= %d \n',off12b, off13b, off23b);
  end
  %%% Scheme 2, END %%%
  %%%%%%%%%%%%%% Different weighting schemes of CC %%%%%%%%%%%%%%
  
  %obtain the target index, and compute the amplitude from master raw CC
  for i = 1: nsta
    pkstmp = pkssta{nit}{i};
    idximp(i) = pkstmp(bsub(i),1)+round(cclen/2); % convert it to index of the signal (or sigdecon)
    idxcoef = idximp(i)-itwlet(i)+1+nfft;  % index of the raw CC that gives the amp
    lagsamp(i) = lag(idxcoef); % offset in samples, in the context of raw CC
    amp(i) = coef(idxcoef, i)/sum(wlet(:, i).^2); % convert max raw CC to amplitude
%     amp(i) = pkstmp(bsub(i),2)/sum(wlet(:, i).^2); % this should give the same result as above
  end

  %%%A plot to check the weighted CC of found eligible peaks 
  if fpltchk
    figure
    p1=plot(ircc,coefeff(:,1),'r-'); hold on
    p2=plot(ircc,coefeff(:,2),'b-');
    p3=plot(ircc,coefeff(:,3),'k-');
%     scatter(pkss1(:,1)+round(cclen/2),pkss1(:,3),20,'r');   % loc of weighted peaks at each station
%     scatter(pkss2(:,1)+round(cclen/2),pkss2(:,3),20,'b');
%     scatter(pkss3(:,1)+round(cclen/2),pkss3(:,3),20,'k');
    p4=plot(ircc,rcc,'o','color',[.7 .7 .7],'markersize',2);
    scatter(idximp(1),pkss1(bsub(1),3),40,'r','filled');  % chosen peaks
    scatter(idximp(2),pkss2(bsub(2),3),40,'b','filled');
    scatter(idximp(3),pkss3(bsub(3),3),40,'k','filled');
    ax = gca;
    plot(ax.XLim, [ampnoi(1) ampnoi(1)],'r--');
    plot(ax.XLim, [ampnoi(2) ampnoi(2)],'b--');
    plot(ax.XLim, [ampnoi(3) ampnoi(3)],'k--');
    legend([p1,p2,p3,p4],'PGC res-wlet CC','SSIB res-wlet CC','SILB res-wlet CC',...
      '3-station RCC');
    box on; grid on
    xlabel(sprintf('Samples at %d Hz',round(1/dt)),'FontSize',12);
    ylabel('Amplitude','FontSize',12);
    title(sprintf('Iteration: %d',nit));
  end
  
  %%%notify if index is out of bounds
  if sum(idximp<1)+sum(idximp>nfft)>0   
    disp(idximp)
    nit = nit-1;  % -1 to recognize the current iteration does not proceed to the end
    fprintf('Forced stop at iteration %d \n', nit);
    fprintf('Because the best index found is outside the time range \n'); 
    fprintf('at the next iteration \n');
    break
%     keyboard
  end

  for i = 1: nsta
    %update the deconvolved impulse response array
    sigdecon(idximp(i), i) = sigdecon(idximp(i), i)+amp(i);
    ampit(nit, (i-1)*2+1:i*2) = [idximp(i) amp(i)];
    impchg = zeros(nfft,1);   % the change of impulse resulting from the current iteration
    impchg(idximp(i)) = amp(i);
    %change in predicted signal by convolving the change of impulse with wavelet
    tmp = conv(wlet(:,i), impchg, 'full'); 
    predchg(:,i) = tmp(itwlet(i):nfft+itwlet(i)-1);  % cut accordingly
    pred(:,i) = pred(:,i) + predchg(:,i);  % update the prediction
    res_new(:,i) = sig(:,i)-pred(:,i);   % residual
  
  %   dres = abs((sum(res.^2)-sum(res_new.^2))./sum(res.^2));  % target evaluation objective
    dres(i) = abs(norm(res(:,i)).^2-norm(res_new(:,i)).^2)./norm(res(:,i)).^2;  % this is equivalent to the above
  %   dres = sum((res-res_new).^2)./sum(res.^2);  % alternative evaluation objective
  %   mfit = sum(res_new.^2)./sum(sig.^2);  % misfit objective, we want the residual to decrease towards 0
    mfit(i) = norm(res_new(:,i));  % alternative misfit objective
  end
  
  aa = [off12b off13b];
  bb = [idximp(1)-idximp(2) idximp(1)-idximp(3)];
  if ~isequaln(aa,bb)
    keyboard
  end
  ampit(nit, nsta*2+1:nsta*2+3) = [off12b off13b off23b];
%   ampit(nit, nsta*2+1:nsta*2+2) = [idximp(1)-idximp(2) idximp(1)-idximp(3)];
    
  %%%may have a different choice here, do we sum all stations here?
  dresit(nit) = sum(dres*100)/nsta;   % store them, make the relative change in percent
  mfitit(nit) = sum(mfit)/nsta;

  %AIC criterion, eq 28 of Ohta&Ide2017JGRSE
  %this can be used only when there will be a global minimum and pre-compute for a large iterations
  dof = nsta*nit;  % number of free params, for each iteration nit, only the index time idximp is free
  if var(noi(:,1))>0
    AIC(nit) = sum(mfit.^2)/sum(var(noi))+2*dof; % variance of background noise 
%     AIC(nit) = sum(mfit.^2/(lsig-1))/sum(var(noi))+2*dof; % variance of background noise 
  else
    AIC(nit) = 0;
  end
    
  nimp = sum(sigdecon>0, 1); % number of non-zero impulses

  %%%Updating plot for each iteration
  if fpltit
    f1.fig = figure(100); clf(f1.fig)
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 8;   % maximum height allowed is 11 inches
    [scrsz, resol] = pixelperinch(1);
    set(f1.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);
    
    subplot(711); box on; grid on; hold on
    p1=plot(wlet(:,1),'r-'); plot(wlet(:,2),'b-'); plot(wlet(:,3),'k-');
    text(0.98,0.85,'Template','Units','normalized','HorizontalAlignment','right','FontSize',10);
    xlim([0,lwlet]);
    title(sprintf('Iteration: %d',nit)); hold off
    
    subplot(712); box on; grid on; hold on
    yyaxis left; p1=plot(sig(:,1), 'r-'); plot(sig(:,2), 'b-'); plot(sig(:,3), 'k-'); 
    ax=gca; yran = [-max(max(abs(sig))) max(max(abs(sig)))];
    text(0.98,0.85,'Signal','Units','normalized','HorizontalAlignment','right','FontSize',10);
    xlim([0,nfft]); ylim(yran); 
    yyaxis(ax,'right');
    plot(ax,irccs,rccs,'o','color',[.8 .8 .8],'markersize',2);
    ylim(ax,[-1 1]);
    hold off
    
    subplot(713); box on; grid on; hold on
    yyaxis left; p1=plot(res(:,1), 'r-');
    plot(res(:,2), 'b-'); plot(res(:,3), 'k-'); ylim(yran);
    yyaxis right;
    p2=plot(ircc,rcc,'o','color',[.8 .8 .8],'markersize',2);
    ylim([-1 1]);
    text(0.98,0.85,'Residual bf. current iter.','Units','normalized','HorizontalAlignment',...
      'right','FontSize',10);
    xlim([0,nfft]);  hold off
        
    subplot(714); box on; grid on; hold on;
    imptemp = find(sigdecon(:,1)>0); p1=stem(imptemp,sigdecon(imptemp,1), 'r','MarkerSize',4);
    text(0.05,0.85,num2str(length(imptemp)),'color','r','FontSize',10,'Units','normalized');
    imptemp = find(sigdecon(:,2)>0); stem(imptemp,sigdecon(imptemp,2), 'b','MarkerSize',4);
    text(0.1,0.85,num2str(length(imptemp)),'color','b','FontSize',10,'Units','normalized');
    imptemp = find(sigdecon(:,3)>0); stem(imptemp,sigdecon(imptemp,3), 'k','MarkerSize',4);
    text(0.15,0.85,num2str(length(imptemp)),'color','k','FontSize',10,'Units','normalized');
    text(0.05,0.65,num2str(idximp(1)),'color','r','FontSize',10,'Units','normalized');
    text(0.1,0.65,num2str(idximp(1)-idximp(2)),'color','b','FontSize',10,'Units','normalized');
    text(0.15,0.65,num2str(idximp(1)-idximp(3)),'color','k','FontSize',10,'Units','normalized');
    plot([0 nfft],[ampnoin(1) ampnoin(1)],'r--');
    plot([0 nfft],[ampnoin(2) ampnoin(2)],'b--');
    plot([0 nfft],[ampnoin(3) ampnoin(3)],'k--');
    text(0.98,0.85,'Impulses','Units','normalized','HorizontalAlignment',...
      'right','FontSize',10); xlim([0,nfft]); hold off
    
    subplot(715); box on; grid on; hold on   
    p1=plot(pred(:,1), 'r-'); plot(pred(:,2), 'b-'); plot(pred(:,3), 'k-'); 
    text(0.98,0.85,'Prediction','Units','normalized','HorizontalAlignment',...
      'right','FontSize',10);
    xlim([0,nfft]); ylim(yran); hold off
    
    subplot(716); box on; grid on; hold on;   
    yyaxis left; p1=plot(res_new(:,1), 'r-'); 
    plot(res_new(:,2), 'b-'); plot(res_new(:,3), 'k-'); ylim(yran);
    [irccn,rccn12] = RunningCC(res_new(:,1), res_new(:,2), cclen);
    [~,rccn13] = RunningCC(res_new(:,1), res_new(:,3), cclen);
    [~,rccn23] = RunningCC(res_new(:,2), res_new(:,3), cclen);
    rccn = (rccn12+rccn13+rccn23)/3;
    yyaxis right; plot(irccn,rccn,'o','color',[.8 .8 .8],'markersize',2);
    ylim([-1 1]);
    text(0.98,0.85,'Residual af. current iter.','Units','normalized','HorizontalAlignment',...
      'right','FontSize',10);
    xlim([0,nfft]); hold off
    
    subplot(717); box on; grid on; hold on;   
    p1=plot(noi(:,1), 'r-'); plot(noi(:,2), 'b-'); plot(noi(:,3), 'k-');
    text(0.98,0.85,'Noise','Units','normalized','HorizontalAlignment',...
      'right','FontSize',10);
    xlim([0, max(nfft, size(noi,1))]); ylim(yran); hold off
    xlabel(sprintf('Samples at %d Hz',round(1/dt)),'FontSize',10);
    ylabel('Amplitude','FontSize',10);
    
    print(f1.fig,'-djpeg','-r300',strcat('/home/chaosong/Pictures/JointIter',num2str(nit),'.jpg'));

    fighdl{1} = f1;
  end

  res = res_new; % update the residual

end  % iteration stops


%% Post-processing
%Disp information if some stopping criteria have been met 
if sum(idximp<1)+sum(idximp>nfft)==0   % notify if index is out of bounds
  if sum(dres>dres_min)==0
    fprintf('Relative change in L2-norm of residual at all stations has reached dres_min=%f%% \n',...
      dres_min*100);
  end
  if sum(mfit>mfit_min)==0
    fprintf('L2-norm of residual at all stations has reached mfit_min=%e \n',mfit_min);
  end
  if sum(mad(res)>mad(noi))==0
    fprintf('MAD of residual has been no larger than that of noise at all stations \n');
  end
  if sum(max(abs(predchg))>mad(noi))==0
    fprintf('Max of change in prediction has been no larger than MAD of noise at all stations \n');
  end
  if sum(amp>ampnoin)==0
    fprintf('Ampltidue of new impulses has been no larger than noise level at all stations \n');
  end
  if nit>=nit_max
    fprintf('Number of iterations has reached nit_max=%d \n',nit_max);
  end
  if nimp>=nimp_max
    fprintf('Number of impulses has reached nimp_max=%d \n',nimp_max);
  end
end

toc


%% plot for the final iteration
if fpltend
  f2.fig = figure(200);clf(f2.fig);
  widin = 8;  % maximum width allowed is 8.5 inches
  htin = 8;   % maximum height allowed is 11 inches
  [scrsz, resol] = pixelperinch(1);
  set(f2.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);
  nrow = 7;
  ncol = 1;
  for isub = 1:nrow*ncol
    f2.ax(isub) = subplot(nrow,ncol,isub);
  end
  
  ax=f2.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); xlim(ax,[0,lwlet]);
  p1=plot(ax,wlet(:,1),'r-'); 
  plot(ax,wlet(:,2),'b-'); plot(ax,wlet(:,3),'k-');
  %   legend(ax,p1,'Wavelet');
  text(ax,0.98,0.85,'Template','Units','normalized','HorizontalAlignment','right','FontSize',10);
  title(ax,sprintf('Total Iterations: %d',nit)); hold(ax,'off');
    
  ax=f2.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); xlim(ax,[0,lsig]);
  plot(ax,[0 nfft],[ampnoin(1) ampnoin(1)],'r--');
  plot(ax,[0 nfft],[ampnoin(2) ampnoin(2)],'b--');
  plot(ax,[0 nfft],[ampnoin(3) ampnoin(3)],'k--');
  imptemp = find(sigdecon(:,1)>0); p1=stem(ax,imptemp,sigdecon(imptemp,1), 'r','MarkerSize',4); 
  text(ax,0.05,0.85,num2str(length(imptemp)),'color','r','FontSize',10,'Units','normalized');
  imptemp = find(sigdecon(:,2)>0); stem(ax,imptemp,sigdecon(imptemp,2), 'b','MarkerSize',4);
  text(ax,0.1,0.85,num2str(length(imptemp)),'color','b','FontSize',10,'Units','normalized');
  imptemp = find(sigdecon(:,3)>0); stem(ax,imptemp,sigdecon(imptemp,3), 'k','MarkerSize',4);
  text(ax,0.15,0.85,num2str(length(imptemp)),'color','k','FontSize',10,'Units','normalized');
  text(ax, 0.98,0.85,'Impulses','Units','normalized','HorizontalAlignment',...
    'right','FontSize',10);  hold(ax,'off');
%   legend(ax,p1,'Impulses'); 

  yran = [-max(max(abs(sig))) max(max(abs(sig)))];
  ax=f2.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); xlim(ax,[0,lsig]); ylim(ax,yran); 
  yyaxis(ax,'left');
  p1=plot(ax,sig(:,1),'r-'); plot(ax,sig(:,2),'b-'); plot(ax,sig(:,3),'k-');
  yyaxis(ax,'right');
  plot(ax,irccs,rccs,'o','color',[.8 .8 .8],'markersize',2);
  ylim(ax,[-1 1]);
  text(ax, 0.98,0.85,'Signal','Units','normalized','HorizontalAlignment',...
    'right','FontSize',10); hold(ax,'off');

  for i = 1: nsta
    ax=f2.ax(3+i); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); xlim(ax,[0,lsig]);
    yyaxis(ax,'left');
    p1=plot(ax,sig(:,i), 'k-'); ylim(ax,yran);
    p2=plot(ax,pred(:,i), 'b-');
    yyaxis(ax,'right');
    [irccsp, rccsp] = RunningCC(sig(:,i), pred(:,i), cclen);
    plot(ax,irccsp,rccsp,'o','color',[.8 .8 .8],'markersize',2);
    ylim(ax,[-1 1]);
    legend(ax,[p1,p2],'Signal','Prediction','location','southeast'); hold(ax,'off');
  end  
  
  ax=f2.ax(4+nsta); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); xlim(ax,[0,lsig]); ylim(ax,yran); 
  yyaxis(ax,'left');
  p1=plot(ax,res(:,1),'r-'); plot(ax,res(:,2),'b-'); plot(ax,res(:,3),'k-');
%   p2=plot(ax,noi,'-','color',[.7 .7 .7]);
%   legend(ax,[p1,p2],'Residual','Noise');
  ylabel(ax,'Amplitude','FontSize',10);
  xlabel(ax,sprintf('Samples at %d Hz',round(1/dt)),'FontSize',10);
  yyaxis(ax,'right');
  [irccr,rccr12] = RunningCC(res(:,1), res(:,2), cclen);
  [~,rccr13] = RunningCC(res(:,1), res(:,3), cclen);
  [~,rccr23] = RunningCC(res(:,2), res(:,3), cclen);
  rccr = (rccr12+rccr13+rccr23)/3;
  plot(ax,irccr,rccr,'o','color',[.8 .8 .8],'markersize',2);
  ylim(ax,[-1 1]);
  text(ax, 0.98,0.85,'Residual','Units','normalized','HorizontalAlignment',...
    'right','FontSize',10); 
  hold(ax,'off');
  fighdl{2} = f2;

  f3.fig = figure(300); clf(f3.fig);
  ax = gca; box on; grid on;
  hold(ax,'on');
  yyaxis(ax,'left');
  p1=plot(ax,1:nit,dresit,'b-');
  scatter(ax,nit,dresit(end),20,'b','filled');
  text(ax,round(nit*4/5), dresit(end)+0.15, num2str(dresit(end)));
  ylabel(ax,'Average relative change (%)');
  yyaxis(ax,'right');
  p2=plot(ax,1:nit,mfitit,'r-');
  scatter(ax,nit,mfitit(end),20,'r','filled');
  text(ax,round(nit*4/5), mfitit(end)+0.15, num2str(mfitit(end)));
  ylabel(ax,'L2-norm');
  legend(ax,[p1,p2],'Average relative change in L2-norm of residual (dres)','Average L2-norm of residual, or data misfit (mfit)');
  xlabel(ax,'Iteration number');
  hold(ax,'off');
  fighdl{3} = f3;
  
  f4.fig = figure(400); clf(f4.fig);
  plot(log(AIC)); hold on;
  ind = find(AIC==min(AIC)); scatter(ind,log(AIC(ind)),'ro');
  text(0.4,0.5,sprintf('Optimal num. of iterations (min. AIC): %d', ind),'Units','normalized',...
    'HorizontalAlignment','center');
  xlabel('Iteration number');
  ylabel('Log_{10}(AIC)');
  fighdl{4} = f4;
end  


% keyboard
