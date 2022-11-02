function [sigdecon,pred,res,dresit,mfitit,ampit,fighdl] = ...
  iterdecon_4thsta(sig,wlet,rcc,dt,twlet,impindep,stas,off1ia,loff_max,fpltit,fpltend,fpltchk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [sigdecon,pred,res,dresit,mfitit,nit,rf] = ...
%   iterdecon(sig,wlet,dt,twlet,width,dres_min,mfit_min,nit_max,fpltit,fpltend)
%
% For deconvolution theory, go to 'iterdecon.m'. This code is trying to do a 
% similar 'deconvolution' at a 4th station using sources deconvolved from
% trio stations and possibly 2ndary sources have been removed. Utilizing the
% expectation of arrival time at a 4th station for a given source (from plane
% fitting, etc), we try to search for the closest and most important peak 
% of the res-wlet CC. 
% --Different from trio station decon, the process starts from the early
%   sources to later ones, and this 'iterative' process doesn't stop until 
%   some criteria are met, eg., residual change between iterations, misfit,
%   max iteration number, it ONLY depends on the input sources. 
% --If within an allowable range of offset in the vicinity of the expected
%   arrival time, peaks could not be found for some sources, likely they 
%   should be thrown away. Therefore, this decon at 4th stations could help
%   decrease the number of sources even further, which is different from
%   directly removing sources from 4th stations using arrival predictions.
%    
%
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/11/01
% Last modified date:   2022/11/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% default value for easy debugging or no input
% defval('shift','0');
% defval('width','2.5');
% defval('dfit_min','0.5');
% defval('nit_max','400');
if isempty(twlet)
  twlet = 0;
end
if isempty(rcc)
  rcc = ones(size(sig));
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

itwlet = round(twlet/dt); % time shift of main arrival of wavelet in samples
sps = round(1/dt);

fighdl = cell(4,1);

%% Preparation, iteration 0
sig = reshape(sig, [],1); % reshape to a vector if not
lsig = size(sig, 1);   % length of signal
wlet = reshape(wlet, [],1); % reshape to a vector if not
lwlet = size(wlet, 1);  % length of wavelet, usually shorter than signal
nfft = lsig;  % otherwise, the length to that of the signal is enough, no padding needed
if lwlet>lsig   % if the wavelet is longer than signal, use way 1 to do CC
  CCOPT = 1;
else
  CCOPT = 2; % otherwise, use way 2 to do CC for efficiency
end

%%% iteration 0
% nit = 0;  % iteration number
res = sig;  % residual
pred = zeros(nfft,1);  % predefine the prediction array
sigdecon = zeros(nfft,1);  % predefine the response array
dres = 1;   % decrease of residual
mfit = max(1, norm(sig));   % misfit between data and prediction
ampit = [];   % amplitude of impulse relative to the max amp of the template
dresit = [];  % relative change in residual for each iteration
mfitit = [];  % misfit for each iteration, norm of the residual 
nimp = 0; % total number of impulses 
predchg = sig;   % change in prediction 
amp = 1;   % intial amplitude of impulse


%% iteration >= 1 
for isrc = 1: size(impindep)
%   nit = nit+1;
  
  %%%%%%%%%%%%%% Different ways to do CC %%%%%%%%%%%%%%
  if CCOPT == 1
    %%% Way 1: redo CC for each iteration, START %%%
    %cross-correlation between the residual and wavelet
    [coef, lag] = xcorr(res, wlet, nfft, 'none'); % unnormalized master raw CC
    %%% Way 1, END %%%
  elseif CCOPT == 2
    %%% Way 2: do the CC once, then update it, efficient, START %%%
    if isrc == 1
      %cross-correlation between the signal and wavelet
      [coef, lag] = xcorr(sig, wlet, nfft, 'none'); % unnormalized master raw CC
    else
      st = max(idximp-itwlet+1, 1);  % start index of predication change, in case < 1
      ed = min(idximp-itwlet+lwlet, nfft);  % end index of predication change, in case > nfft
      %cross-correlation between the predication change and wavelet
      [coefchg, ~] = xcorr(predchg(st: ed), wlet, lwlet, 'none'); % raw CC
      [~, ind] = max(coefchg);   % max of raw cc
      
      cnt = lagsamp+nfft+1;  % center index of master raw cc to be updated
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
      coef(indchg) = coef(indchg)-coefchg; % note the neg. sign here
    end
    %%% Way 2, END %%%
  end
  %%%%%%%%%%%%%% Different ways to do CC %%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%% Different weighting schemes of CC %%%%%%%%%%%%%%
  %%% Scheme 2: rank all peaks in master CC by the dot product between the rcc and master CC at each
  %%% peak, and choose the one with the highest dot product
  %%% this would perserve the true peaks in master CC without shifting, and no negative
  %%% amplitude would come out
  lrcc = length(rcc); % length of running CC
  ldiff = nfft-lrcc;  % difference in length
  %effective raw CC that corresponds to the index of the overlapping portion between signal (or 
  %sigdecon) and running CC, same length as rcc
  coefeff = coef(lag+itwlet >= 1+round(ldiff/2) & ...
                 lag+itwlet <= nfft-round(ldiff/2));          
  
  %%%Here there are multiple ways to find the eligible peaks, with/without considering the noise level             
  %find all peaks in the effective master raw CC, return the peak height and index in effective raw 
  %CC, context of overlapping portion
  [pkhgt, pkind] = findpeaks(coefeff); 
%   [pkhgt, pkind] = findpeaks(coefeff,'MinPeakHeight',ampnoi); 
%   [pkhgt, pkind] = findpeaks(coefeff,'MinPeakProminence',ampnoi); % promience is good, but do we have enough information on how large it should be?
%   [pkhgt, pkind] = findpeaks(coefeff,'MinPeakHeight',ampnoi,'MinPeakProminence',ampnoi);
  if isempty(pkind)
    ampit(isrc, 1:6) = zeros(1,6);
    dresit(isrc) = -1;
    mfitit(isrc) = -1;
    continue  % go to next source
  end
  %rcc serves the weight as the peak height, aka the master raw cc value at the peak
  wtcoef = [pkind rcc(pkind).* pkhgt]; % weighted master raw cc
  %%% Scheme 2, END %%%
  %%%%%%%%%%%%%% Different weighting schemes of CC %%%%%%%%%%%%%%

  %%%%%%%%%%%% Different ways to determine the 'best' peak %%%%%%%%%%%%%%%%
  %%% Way2: find all pairs that are close enough, choose the one with highest weighted CC
  %source arrival prediction from plane fitting
  idxpred = pred_tarvl_at4thsta(stas,impindep(isrc,7),impindep(isrc,8),impindep(isrc,1),off1ia);
  indclose = find(abs(wtcoef(:,1)-(idxpred-round(ldiff/2)))<=loff_max);
  if isempty(indclose)
    ampit(isrc, 1:6) = zeros(1,6);
    dresit(isrc) = -1;
    mfitit(isrc) = -1;
    continue
  else
    [mwtcoef,ind] = max(wtcoef(indclose,2));  %choose the one with highest weighted CC
  end
  idximp = wtcoef(indclose(ind),1)+round(ldiff/2);  % convert it to index of the signal (or sigdecon)
  idxcoef = idximp-itwlet+1+nfft;  % index of the raw CC that gives the amp
  lagsamp = lag(idxcoef); % offset in samples, in the context of raw CC
  amp = coef(idxcoef)/sum(wlet.^2); % convert max raw CC to amplitude
  rccimp = rcc(wtcoef(indclose(ind),1)); % rcc value at the index of the selected impulse 
  %%% Way2: find all pairs that are close enough, choose the one with highest weighted CC
  %%%%%%%%%%%% Different ways to determine the 'best' peak %%%%%%%%%%%%%%%%

  %update the deconvolved impulse response array
  sigdecon(idximp) = sigdecon(idximp)+amp;
  ampit(isrc, 1:2) = [idximp amp];   %zero-crossing index and amp of decon impulse at this iter  
  [pkhres, pkires] = findpeaks(res); 
  [~,tempind] = min(abs(pkires-idximp));
  ampit(isrc, 3:4) = [pkires(tempind) pkhres(tempind)];  %closest residual waveform positive peak
  [pkhres, pkires] = findpeaks(-res); 
  [~,tempind] = min(abs(pkires-idximp));
  ampit(isrc, 5:6) = [pkires(tempind) -pkhres(tempind)];  %closest residual waveform negative peak  
  impchg = zeros(nfft,1);   % the change of impulse resulting from the current iteration
  impchg(idximp) = amp;
  %change in predicted signal by convolving the change of impulse with wavelet
  predchg = conv(wlet, impchg, 'full'); 
  predchg = predchg(itwlet:nfft+itwlet-1);  % cut accordingly
  pred = pred + predchg;  % update the prediction
  
  res_new = sig-pred;   % residual
%   dres = abs((sum(res.^2)-sum(res_new.^2))./sum(res.^2));  % target evaluation objective
  dres = abs(norm(res).^2-norm(res_new).^2)./norm(res).^2;  % this is equivalent to the above
%   dres = sum((res-res_new).^2)./sum(res.^2);  % alternative evaluation objective
%   mfit = sum(res_new.^2)./sum(sig.^2);  % misfit objective, we want the residual to decrease to 0
  mfit = norm(res_new);  % alternative misfit objective
  dresit(isrc) = dres*100;   % store them, make the relative change in percent
  mfitit(isrc) = mfit;
  
  nimp = sum(sigdecon>0); % number of non-zero impulses

  res = res_new; % update the residual

end  % iteration stops







