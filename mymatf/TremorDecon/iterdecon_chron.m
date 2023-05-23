function [sigdecon,pred,res,dresit,mfitit,ampit,nit,fighdl,rf] = ...
  iterdecon_chron(sig,wlet,rcc,noi,fixthr,dt,twlet,width,dres_min,mfit_min,nit_max,nimp_max,fpltit,fpltend,fpltchk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [sigdecon,pred,res,dresit,mfitit,nit,rf] = ...
%   iterdecon_chron(sig,wlet,dt,twlet,width,dres_min,mfit_min,nit_max,fpltit,fpltend)
%
% This function aims to do iterative deconvolution as well. But different
% from carrying it out following the descending order of ranking the product 
% of RCC and unnormalized res-wlet CC, we try to follow the chronological 
% order in time. In other words, we deconvolve the peak in res-wlet CC that 
% is above the threshold if any and the earliest in time. If the arrival is
% high in amp, logically it is right to resolve it early; if not, the related
% coda from this source won't impact later sources too much. This is one of
% the several ways designed to change the deconvolution order, for example,
% modifying RCC by shifting some stations, using different station pairs, 
% or not using RCC, etc. 
% 
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
%   dresit: percentage of relative change in residual of each iteration, in terms of variance
%           reduction
%   mfitit: norm of the residual/misfit of each iteration
%   ampit:  amplitude of decolved impulse of each iteration
%   nit:  actual number of iterations until algorithm converges/stops
%   fighdl: handles of figures plotted
%   rf: receiver function if needed, MAY BE ABANDONED LATER
%
%
% NOTE:

%
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/02/06
% Last modified date:   2023/02/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%% default value for easy debugging or no input
% defval('shift','0');
% defval('width','2.5');
% defval('dfit_min','0.5');
% defval('nit_max','400');
if isempty(twlet)
  twlet = 0;
end
if isempty(width)
  width = 2.5;
end
if isempty(dres_min)
  dres_min = 0.5; % 0.5 percent
end
if isempty(mfit_min)
  mfit_min = norm(sig);
end
if isempty(nit_max)
  nit_max = 400;
end
if isempty(nimp_max)
  nimp_max = 400;
end
if isempty(rcc)
  rcc = ones(size(sig));
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

itwlet = round(twlet/dt); % time shift of main arrival of wavelet in samples
sps = round(1/dt);

fighdl = cell(4,1);

%% Preparation, iteration 0
sig = reshape(sig, [],1); % reshape to a vector if not
lsig = size(sig, 1);   % length of signal
wlet = reshape(wlet, [],1); % reshape to a vector if not
lwlet = size(wlet, 1);  % length of wavelet, usually shorter than signal
if nargout == 9   % meanning rf is also among the outputs
  nfft = pow2(nextpow2(lsig));  % length that is easy for fft, etc. NOT necessary if no fft is needed
  %pad input with zeros to the same length of next power of 2
  sig = [sig; zeros(nfft-lsig, 1)];
%   wlet = [wlet; zeros(nfft-lwlet, 1)];
else
  nfft = lsig;  % otherwise, the length to that of the signal is enough, no padding needed
end

if lwlet>lsig   % if the wavelet is longer than signal, use way 1 to do CC
  CCOPT = 1;
else
  CCOPT = 2; % otherwise, use way 2 to do CC for efficiency
end

%%% iteration 0
nit = 0;  % iteration number
res = sig;  % residual
pred = zeros(nfft,1);  % predefine the prediction array
sigdecon = zeros(nfft,1);  % predefine the response array
varred = 1;   % decrease of residual, in terms of variance reduction
mfit = max(1, norm(sig));   % misfit between data and prediction
ampit = [];   % amplitude of impulse relative to the max amp of the template
dresit = [];  % relative change in residual for each iteration, in terms of variance reduction
mfitit = [];  % misfit for each iteration, norm of the residual 
nimp = 0; % total number of impulses 
predchg = sig;   % change in prediction 
amp = 1;   % intial amplitude of impulse

%%%%%%%%%%%%%%%%%%
%get the master CC between sig and wlet, to know what is range of weighted master CC
[coef, lag] = xcorr(sig, wlet, nfft, 'none'); % unnormalized master raw CC
lrcc = length(rcc); % length of running CC
ldiff = nfft-lrcc;  % difference in length
%effective raw CC that corresponds to the index of the overlapping portion between signal (or
%sigdecon) and running CC, same length as rcc
coefeff = coef(lag+itwlet >= 1+round(ldiff/2) & ...
  lag+itwlet <= nfft-round(ldiff/2));
%find all peaks in the effective master raw CC
[pkhgt, pkind] = findpeaks(coefeff);
%rcc serves the weight as the peak height, aka the master raw cc value at the peak
wtcoef = rcc(pkind).* pkhgt;
if isempty(fixthr)
  medwtcoef = median(wtcoef);  % median of the weighted master CC, could be percentile?
% medwtcoef = prctile(wtcoef,25);
% medwtcoef = 1.2352e-01; % for test purpose if the records are like noise, sta 1
% medwtcoef = 7.3467e-02; % for test purpose if the records are like noise, sta 2
% medwtcoef = 1.0873e-01; % for test purpose if the records are like noise, sta 3

% medwtcoef = 2.3910e-01; % for test purpose if the records are like noise, sta 1, concatenated case
% medwtcoef = 1.4577e-01; % for test purpose if the records are like noise, sta 2
% medwtcoef = 1.9032e-01; % for test purpose if the records are like noise, sta 3
else
  medwtcoef = fixthr;
end
medwtcoef

mwtcoef = 100;  % intial amplitude of max master CC weighted by rcc, value doesn't really matter 
medrccpeak = median(rcc(pkind));
rccimp = 1;
%%%%%%%%%%%%%%%%%

[coefnoi, lagnoi] = xcorr(noi, wlet, size(noi, 1), 'none'); % unnormalized cc
coefnoieff = coefnoi(lagnoi+itwlet >= 1 & lagnoi+itwlet <= size(noi, 1));
[noipkhgt, noipkind] = findpeaks(coefnoieff);
ampnoi = median(noipkhgt);
% [ampnoi, ~] = max(coefnoi);   % a measure of noise level
% figure
% plot(coefnoi); hold on
% plot(coefnoieff);
% scatter(noipkind,noipkhgt,20,'k');

% noiwletrat = rms(noi)/max(wlet);  % ratio of amplitude of noise relative to the max amp of wavelet
ampnoin = ampnoi/sum(wlet.^2);   % normalized amp of CC between noise and wavelet
ampnoiwtn = ampnoi*median(rcc)/sum(wlet.^2); % normalized amp of CC between noise and wavelet weighted by median rcc
ampnoiwt = ampnoi*median(rcc); % amp of CC between noise and wavelet weighted by median rcc

%% iteration 
while mwtcoef>1.05*medwtcoef   %if discard hard limit on max num of iters

  nit = nit+1;

  %%%%%%%%%%%%%% Different ways to do CC %%%%%%%%%%%%%%
  if CCOPT == 1
    %%% Way 1: redo CC for each iteration, START %%%
    %cross-correlation between the residual and wavelet
    [coef, lag] = xcorr(res, wlet, nfft, 'none'); % unnormalized master raw CC
    %%% Way 1, END %%%
  elseif CCOPT == 2
    %%% Way 2: do the CC once, then update it, efficient, START %%%
    if nit == 1
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

  %%% Scheme 3: still obtain the dot product between the rcc and master CC at each
  %%% peak, and choose the one with higher than threshold and the earliest
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
  if isempty(pkind)   % force stop if there are no more eligible peaks
    nit = nit-1;  % -1 to recognize the current iteration does not proceed to the end
    fprintf('Forced stop at iteration %d \n', nit);
    fprintf('Because no eligible peaks are found (depending on the searching criterion) \n'); 
    fprintf('at the next iteration \n');
    break
  end
  %rcc serves the weight as the peak height, aka the master raw cc value at the peak
  ind = find(rcc(pkind).* pkhgt>=1.05*medwtcoef); %require peaks beyond threshold
  if isempty(ind)
    nit = nit-1;  % -1 to recognize the current iteration does not proceed to the end
    fprintf('Forced stop at iteration %d \n', nit);
    fprintf('Because no peaks beyond threshold are found (depending on the searching criterion) \n'); 
    fprintf('at the next iteration \n');
    break
  end

  %use chronologically the first one
  mwtcoefidx = ind(1);
  mwtcoef = rcc(pkind(mwtcoefidx)) .* pkhgt(mwtcoefidx);

  idximp = pkind(mwtcoefidx)+round(ldiff/2); % convert it to index of the signal (or sigdecon)
  idxcoef = idximp-itwlet+1+nfft;  % index of the raw CC that gives the amp
  lagsamp = lag(idxcoef); % offset in samples, in the context of raw CC
  amp = coef(idxcoef)/sum(wlet.^2); % convert max raw CC to amplitude
  rccimp = rcc(pkind(mwtcoefidx)); % rcc value at the index of the selected impulse  
  %%% Scheme 2, END %%%
  %%%%%%%%%%%%%% Different weighting schemes of CC %%%%%%%%%%%%%%

  %%%A plot to check the weighted CC of found eligible peaks 
  if fpltchk
    figure
    ax = gca;
    hold(ax,'on');
    ax.Box='on'; grid(ax,'on');
    yyaxis(ax,'right');
    p2=plot(ax,(1+round(ldiff/2)): (nfft-round(ldiff/2)),rcc,'o','color',[.7 .7 .7],'markersize',2);  % rcc
    scatter(ax,pkind+round(ldiff/2),rcc(pkind),20,'r');  % rcc value at peaks
    ylim(ax,[-1 1]);
    ylabel('RCC','FontSize',12);
    
    yyaxis(ax,'left');
    p1=plot(ax,(1+round(ldiff/2)): (nfft-round(ldiff/2)),coefeff,'-','color',[.3 .3 .3]); hold on % master cc
    scatter(ax,pkind+round(ldiff/2),pkhgt,20,'r'); % peaks
%     scatter(idximp,pkhgt(mwtcoefidx),10,'b','filled');  % chosen peak
%     scatter(idximp,rcc(pkind(mwtcoefidx)),10,'b','filled'); % chosen peak
    p3=scatter(ax,pkind+round(ldiff/2),rcc(pkind).* pkhgt,20,'b');  % peak heights weighted by rcc
    scatter(ax,idximp,rcc(pkind(mwtcoefidx)).*pkhgt(mwtcoefidx),20,'b','filled'); % chosen peak
    p4=plot(ax,ax.XLim, [medwtcoef medwtcoef],'k--'); % stopping threshold, subject to actual practice
%     p4=plot(ax,ax.XLim, [ampnoiwt ampnoiwt],'k--'); % noise level
%     p4=plot(ax.XLim, [ampnoi ampnoi],'k--'); % noise level
%     plot(ax.XLim,[mwtcoef mwtcoef],'b--');  % indeed the max?
    text(ax,0.05,0.95,sprintf('%d; %.2f',idximp,rcc(pkind(mwtcoefidx)).*pkhgt(mwtcoefidx)),'Units',...
      'normalized','HorizontalAlignment','left');
    ylabel('Amplitude','FontSize',12);
    xlabel(ax,sprintf('Samples at %d Hz',sps),'FontSize',12);
    title(sprintf('Iteration: %d',nit));
    xlim(ax,[max([0,idximp-500]) min([lsig,idximp+500])]);
    legend(ax,[p1,p2,p3,p4],'Residual-template CC','3-station RCC','Weighted peaks','Stopping threshold');
    hold(ax,'on');
  end

  %%%notify if index is out of bounds
  if idximp<1 || idximp>nfft 
    disp(idximp)
    nit = nit-1;  % -1 to recognize the current iteration does not proceed to the end
    fprintf('Forced stop at iteration %d \n', nit);
    fprintf('Because the best index found is outside the time range \n'); 
    fprintf('at the next iteration \n');
    break
%     keyboard
  end

  %update the deconvolved impulse response array
  sigdecon(idximp) = sigdecon(idximp)+amp;
  ampit(nit, 1:2) = [idximp amp];   %zero-crossing index and amp of decon impulse at this iter  
  [pkhres, pkires] = findpeaks(res); 
  [~,tempind] = min(abs(pkires-idximp));
  ampit(nit, 3:4) = [pkires(tempind) pkhres(tempind)];  %closest residual waveform positive peak
  [pkhres, pkires] = findpeaks(-res); 
  [~,tempind] = min(abs(pkires-idximp));
  ampit(nit, 5:6) = [pkires(tempind) pkhres(tempind)];  %closest residual waveform negative peak  
  impchg = zeros(nfft,1);   % the change of impulse resulting from the current iteration
  impchg(idximp) = amp;
  %change in predicted signal by convolving the change of impulse with wavelet
  predchg = conv(wlet, impchg, 'full'); 
  predchg = predchg(itwlet:nfft+itwlet-1);  % cut accordingly
  pred = pred + predchg;  % update the prediction
  
%   %fit a parabola around the peak, to quantify the 'sharpness' of the found peak
%   nptsfit = 5;
%   x = (pkind(mwtcoefidx)-floor(nptsfit/2): pkind(mwtcoefidx)+floor(nptsfit/2))';
% %   autoccres = xcorr(res); %unnormalized autocorrelation
%   normlizer = sqrt(max(xcorr(res))*max(xcorr(wlet))); %use unnormalized autocorrelations as normlizer
%   y = coefeff(x)/normlizer; %normlization in order for a fair comparison if wavefrom amp is diff.
%   poly = polyfit(x,y,2);  %1st param is focal length that determines how flat or skinny parabola is
%   ampit(nit, 7) = poly(1);  %save 1st param to represent the sharpness    
% %   yfit = polyval(poly,x);

  res_new = sig-pred;   % residual
%   dres = abs((sum(res.^2)-sum(res_new.^2))./sum(res.^2));  % target evaluation objective
%   dres = abs(norm(res).^2-norm(res_new).^2)./norm(res).^2;  % this is equivalent to the above
  varred = abs(var(res)-var(res_new))./var(res);  % again equivalent to the above, but called variance reduction
%   dres = sum((res-res_new).^2)./sum(res.^2);  % alternative evaluation objective
%   mfit = sum(res_new.^2)./sum(sig.^2);  % misfit objective, we want the residual to decrease to 0
  mfit = norm(res_new);  % alternative misfit objective
  dresit(nit) = varred*100;   % store them, make the relative change in percent
  mfitit(nit) = mfit;
  
%   %For benchmark, the following 3 values should be identical according to theory in the reference 
%   sum(res_new.^2) 
%   mfit^2
%   sum(sig.^2)-amp^2*sum(wlet.^2)

  %AIC criterion, eq 28 of Ohta&Ide2017JGRSE
  %this can be used only when there will be a global minimum and pre-compute for a large iterations
  dof = 1*nit;  % number of free params, for each iteration nit, only the index time idximp is free
  if var(noi)>0
    AIC(nit) = mfit^2/var(noi)+2*dof; % variance of background noise
  else
    AIC(nit) = 0;
  end

  nimp = sum(sigdecon>0); % number of non-zero impulses

  %%%Updating plot for each iteration
  if fpltit
    f1.fig = figure(100); clf(f1.fig)
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 9;   % maximum height allowed is 11 inches
    [scrsz, resol] = pixelperinch(1);
    set(f1.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);
    
    subplot(811);
    p1=plot(wlet,'r-'); box on; grid on;
    legend(p1,'Template'); xlim([0,lwlet]); %xlim([0,nfft])
    title(sprintf('Iteration: %d',nit));
    longticks(gca,2);
    if lwlet < lsig
      ax = gca;
      loc = getpos(ax);
      shrink(ax,lsig/lwlet,1);
      ax.Position(1)=loc(1);
%       nolabels(ax,1);
    end

    subplot(812); 
    yyaxis left; p1=plot(sig, 'k-'); box on; grid on; ax=gca; yran = ax.YLim; xlim([0,nfft]);
    yyaxis right;
    p2=plot((1+round(ldiff/2)): (nfft-round(ldiff/2)),rcc,'o','color',[.8 .8 .8],'markersize',2);
    ylim([-1 1]);
    legend([p1,p2],'Signal','3-sta RCC');
%     nolabels(gca,1);
    longticks(gca,2);
    
    subplot(813);
    p1=plot(res, 'c-'); box on; grid on; xlim([0,nfft]); ylim(yran);
    legend(p1,'Residual bf. current iter.'); 
%     nolabels(gca,1);
    longticks(gca,2);
    
    subplot(814);
    yyaxis left;
    p1=plot((1+round(ldiff/2)): (nfft-round(ldiff/2)),coefeff,'-','color',[.2 .2 .2]);
    box on; grid on; xlim([0,nfft]); hold on;
    ax = gca;
    plot(ax.XLim,[ampnoi ampnoi],'--','color',[.5 .5 .5]);
    scatter(pkind+round(ldiff/2),pkhgt,20,'r');
    plot([idximp idximp],ax.YLim,'k--');
    yyaxis right;
    p2=plot((1+round(ldiff/2)): (nfft-round(ldiff/2)),rcc,'o','color',[.8 .8 .8],'markersize',2);
    ylim([-1 1]);
    legend([p1,p2],'Res-temp CC','3-sta RCC');
%     nolabels(gca,1);
    longticks(gca,2); hold off
    
    subplot(815);
%     p1=stem(sigdecon, 'b','MarkerSize',4);
    imptemp = find(sigdecon>0); p1=stem(imptemp,sigdecon(imptemp),'b','MarkerSize',4);
    box on; grid on; hold on; xlim([0,nfft]);
    plot([0 nfft],[ampnoin ampnoin],'--','color',[.5 .5 .5]);
    text(0.05,0.85,num2str(nimp),'FontSize',10,'Units','normalized');
    text(0.05,0.65,num2str(idximp),'FontSize',10,'Units','normalized');
    legend(p1,'Impulses'); 
%     nolabels(gca,1);
    longticks(gca,2);
    
    subplot(816);    
    p1=plot(pred, 'b-');  box on; grid on; xlim([0,nfft]); ylim(yran);
    legend(p1,'Prediction');
%     nolabels(gca,1);
    longticks(gca,2);
    
    subplot(817);   
    p1=plot(res_new, 'm-'); box on; grid on; hold on; xlim([0,nfft]); ylim(yran);
    patarea = [max(idximp-itwlet+1, 1) yran(1);
               min(idximp-itwlet+lwlet, nfft) yran(1);
               min(idximp-itwlet+lwlet, nfft) yran(2);
               max(idximp-itwlet+1, 1) yran(2);
               max(idximp-itwlet+1, 1) yran(1);
               ];
    patch(patarea(:,1),patarea(:,2),'k','Facealpha',0.1,'edgecolor','none');
    plot([idximp idximp],yran,'k--');
    p2=plot(noi,'-','color',[.5 .5 .5]);
    legend([p1,p2],'Residual af. current iter.','Noise');
%     nolabels(gca,1);
    longticks(gca,2); hold off
    
    subplot(818);
    p1=plot(res_new, 'm-'); box on; grid on; hold on; ylim(yran); 
    patch(patarea(:,1),patarea(:,2),'k','Facealpha',0.1,'edgecolor','none');
    plot([idximp idximp],yran,'k--'); 
    legend(p1,'Updated residual-zoom'); 
    xlim([max(idximp-itwlet+1, 1), min(idximp-itwlet+lwlet, nfft)]);  % only zoom at portion changed
    xlabel(sprintf('Samples at %d Hz',sps),'FontSize',10);
    ylabel('Amplitude','FontSize',10);
    longticks(gca,2); hold off
    
%     print(f1.fig,'-djpeg','-r300',strcat('/home/chaosong/Pictures/IndepIter',num2str(nit),'.jpg'));

    fighdl{1} = f1;
  end
  
  res = res_new; % update the residual

end  % iteration stops    

% %% do you want to enforce the chronological order?
% chronoflag = 1;
% 
% if chronoflag
%   [~,ikeep,idiscard] = vectorchron(ampit(:,1));
%   sigdecon(ampit(idiscard,1)) = 0;
%   ampit = ampit(ikeep,:);
%   dresit = dresit(ikeep);
%   mfitit = mfitit(ikeep);
% end

%% Post-processing
sigdecon = sigdecon(1:lsig);  %cut deconvolved signal

%make output a column vector
dresit = dresit(:);
mfitit = mfitit(:);
  
%Disp information if some stopping criteria have been met 
if idximp>=1 && idximp<=nfft   % notify if index is out of bounds
  if varred*100<=dres_min
    fprintf('Relative change in L2-norm of residual has reached dres_min=%f%% \n',dres_min);
  end
  if mfit<=mfit_min
    fprintf('L2-norm of residual has reached mfit_min=%e \n',mfit_min);
  end
  if mad(res)<=mad(noi)
    fprintf('MAD of residual has been no larger than that of noise \n');
  end
  if max(abs(predchg))<=mad(noi)
    fprintf('Max of change in prediction has been no larger than MAD of noise \n');
  end
  if amp<=ampnoin
    fprintf('Ampltidue of new impulses has been no larger than noise level \n');
  end
  if nit>=nit_max
    fprintf('Number of iterations has reached nit_max=%d \n',nit_max);
  end
  if nimp>=nimp_max
    fprintf('Number of impulses has reached nimp_max=%d \n',nimp_max);
  end
end

% if required, create receiver function rf, i.e., impulses with a finite width
if nargout == 9   % meanning rf is also among the outputs
  if width   % if the width is not 0
    %     synpdf = normpdf(xx,mu,sigma);
    Gauss = gauss_zerophase(nfft,dt,width);
    rf = filter_freq(sigdecon,Gauss,dt);
  else
    rf = sigdecon;
  end

  % shift and truncate rf
%   if shift  % if the shift is not o
%     rf = shift_freq(rf, nfft, dt, shift);
%   end
  rf = rf(1:lsig);
end

toc

%% plot for the final iteration
if fpltend
  %%%Signal vs. prediction
  f2.fig = figure(200);clf(f2.fig);
  f2.fig.Renderer = 'painters';
  widin = 8;  % maximum width allowed is 8.5 inches
  htin = 6;   % maximum height allowed is 11 inches
  [scrsz, resol] = pixelperinch(1);
  set(f2.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);
  nrow = 4;
  ncol = 1;
  for isub = 1:nrow*ncol
    f2.ax(isub) = subplot(nrow,ncol,isub);
  end
  
  pltxran = [0.1 0.9]; pltyran = [0.1 0.9];
  pltxsep = 0.02; pltysep = 0.03; 
  %get the locations for each axis
  axpos = optaxpos(f2,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

  ax=f2.ax(1); 
  p1=plot(ax,wlet,'r-'); legend(ax,p1,'Template'); 
  xlim(ax,[0,lwlet]);
  title(ax,sprintf('Total Iterations: %d',nit));
  axsym(ax,2);
  axranexp(ax,6,10);
  if lwlet < lsig
    shrink(ax,lsig/lwlet,1);
    ax.Position(1)=axpos(1,1);
    nolabels(ax,1); 
%   else
%     longticks(ax,2); 
  end
  ax.Box='on'; grid(ax,'on');longticks(ax,2); 

  ax=f2.ax(2);
  yyaxis(ax,'left');
%   p1=stem(ax,sigdecon, 'b','MarkerSize',4);
  imptemp = find(sigdecon>0); p1=stem(ax,imptemp,sigdecon(imptemp),'b','MarkerSize',4);
  hold(ax,'on');
  plot(ax,[0 nfft],[ampnoin ampnoin],'--','color',[.3 .3 .3]);
  text(ax,0.05,0.9,num2str(nimp),'FontSize',12,'Units','normalized');
  axranexp(ax,2,10);  ax.Box='on'; grid(ax,'on');
  yyaxis(ax,'right');
  p2=plot(ax,(1+round(ldiff/2)): (nfft-round(ldiff/2)),rcc,'o','color',[.6 .6 .6],'markersize',2);
  ylim(ax,[-1 1]);
  legend(ax,[p1,p2],'Impulses','3-sta RCC'); xlim(ax,[0,lsig]);
  longticks(ax,2); nolabels(ax,1); hold(ax,'off');
  
  ax=f2.ax(3); 
  yyaxis(ax,'left');
  p1=plot(ax,sig, 'k-'); xlim(ax,[0,lsig]); hold(ax,'on'); 
  p2=plot(ax,pred, 'b-'); yran=1.2*[-max(abs([sig; pred])) max(abs([sig; pred]))]; 
  ylim(ax,yran); ax.Box='on'; grid(ax,'on');
  yyaxis(ax,'right');
  rccmwlen = sps/2;
  [irccsp, rccsp] = RunningCC(sig, pred, rccmwlen);
  p3=plot(ax,irccsp,rccsp,'o','color',[.6 .6 .6],'markersize',2);
  ylim(ax,[-1 1]);
  legend(ax,[p1,p2,p3],'Signal','Prediction','Sig-pred RCC','location','southeast');
  longticks(ax,2); nolabels(ax,1); hold(ax,'off');
  
  ax=f2.ax(4); 
  p1=plot(ax,res, 'c-'); xlim(ax,[0,lsig]); ylim(ax,yran); hold(ax,'on');
  p2=plot(ax,noi,'-','color',[.7 .7 .7]);
  legend(ax,[p1,p2],'Residual','Noise');
  xlabel(sprintf('Samples at %d Hz',sps),'FontSize',12);
  ylabel('Amplitude','FontSize',12);
  longticks(ax,2); ax.Box='on'; grid(ax,'on'); hold(ax,'off'); 
  f2.ax(1).XTick = ax.XTick;
  fighdl{2} = f2;

  
  %%%Iteration performance
  f3.fig = figure(300); clf(f3.fig);
  f3.fig.Renderer = 'painters';
  ax = gca; box on; grid on;
  hold(ax,'on');
  yyaxis(ax,'left');
  p1=plot(ax,1:nit,dresit,'b-');
  scatter(ax,nit,dresit(end),20,'b','filled');
  text(ax,round(nit*4/5), dresit(end)+0.15, num2str(dresit(end)));
  ylabel(ax,'Variance reduction (%)');
  yyaxis(ax,'right');
  p2=plot(ax,1:nit,mfitit,'r-');
  scatter(ax,nit,mfitit(end),20,'r','filled');
  text(ax,round(nit*4/5), mfitit(end)+0.15, num2str(mfitit(end)));
  ylabel(ax,'L2-norm');
  legend(ax,[p1,p2],'Variance reduction of residual (dres)','L2-norm of residual, or data misfit (mfit)');
  xlabel(ax,'Iteration number');
  hold(ax,'off');
  fighdl{3} = f3;
  
  %%%AIC variation
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
  