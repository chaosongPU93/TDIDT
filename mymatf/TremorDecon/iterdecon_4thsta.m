function [sigdecon,pred,res,dresit,mfitit,ampit,fighdl] = ...
  iterdecon_4thsta(sig,wlet,irccran,rcc,fixthr,dt,twlet,impindep,stas,off1i,off1iw,offmax,fpltit,fpltend,fpltchk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [sigdecon,pred,res,dresit,mfitit,ampit,fighdl] = ...
%   iterdecon_4thsta(sig,wlet,rcc,dt,twlet,impindep,stas,off1ia,loff_max,fpltit,fpltend,fpltchk)
%
% For deconvolution theory, go to 'iterdecon.m'. This code is trying to do a 
% similar 'deconvolution' at a 4th station using sources deconvolved from
% trio stations and possibly 2ndary sources have been removed. Utilizing the
% expectation of arrival time at a 4th station for a given source (from plane
% fitting, etc), we try to search for the closest and most important peak 
% of the res-wlet CC. NOTE that in default we use the 1st sta PGC as the
% reference sta.
%
% --Different from trio station decon, the process starts from the early
%   sources to later ones, and this 'iterative' process doesn't stop until 
%   some criteria are met, eg., residual change between iterations, misfit,
%   max iteration number, it ONLY depends on the input sources. 
% --If within an allowable range of offset in the vicinity of the expected
%   arrival time, peaks could not be found for some sources, likely they 
%   should be thrown away. Therefore, this decon at 4th stations could help
%   decrease the number of sources even further, which is different from
%   directly removing sources from 4th stations using arrival predictions.
% --The initial algorithm finds all peaks close enough to expected arrival
%   time from plane fitting, and choose one with the highest weighted height
%   in res-wlet cc. There might be a caveat where the found peak is NOT
%   actually an important peak at KLNB independently.
% --Therefore, revisions are to rank all weighted peaks in res-wlet cc in  
%   each iteration that are above some threshold height, eg,
%   median(rcc* sig-wlet cc), among these fairly important peaks, are there 
%   any within the allowable range from the airrval prediction? If none, 
%   probably the source in question is suspicious; if yes but more than 1,
%   their importances have been determined by ranking already. 
% --Potentially this would eliminate more sources, leading to less residual
%   reduction at both 4th stations and trio stations. The advantage, however,
%   is the guarantee of importance of found peaks independently at 4th station.
% --Update 2022/11/11, if, RCC is concatenated from shorter windows, then there
%   might be diff in overall alignment and sub-win alignment between 4th and
%   ref sta (1st). And predicted arrival is always computed as if there is a
%   source from (0,0) so that only the station location diff is corrected. That
%   being said if the net alignment between 4th and ref sta is NOT 0, predicted
%   arrival needs be corrected before matching the peaks in res-wlet CC. Of
%   course if there is no net alignment or RCC is computed upon the whole win,
%   simply set relavant inputs to 0. Also note that as it focuses on one station
%   per run, 'off1i' has a size of 1*1 and 'off1iw' has a size of nsubwin*1.
%    
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/11/01
% Last modified date:   2022/11/02
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

fighdl = cell(3,1);


%% Preparation, iteration 0
sig = reshape(sig, [],1); % reshape to a vector if not
lsig = size(sig, 1);   % length of signal
wlet = reshape(wlet, [],1); % reshape to a vector if not
lwlet = size(wlet, 1);  % length of wavelet, usually shorter than signal
nfft = lsig;  % otherwise, the length to that of the signal is enough, no padding needed
% if lwlet>lsig   % if the wavelet is longer than signal, use way 1 to do CC
%   CCOPT = 1;
% else
%   CCOPT = 2; % otherwise, use way 2 to do CC for efficiency
% end
%use opt 1 to do master raw CC in case for some sources, no eligible peaks at 4th sta are found 
CCOPT = 1; 

%%% iteration 0
res = sig;  % residual
pred = zeros(nfft,1);  % predefine the prediction array
sigdecon = zeros(nfft,1);  % predefine the response array
ampit = [];   % amplitude of impulse relative to the max amp of the template
dresit = [];  % relative change in residual for each iteration
mfitit = [];  % misfit for each iteration, norm of the residual 
predchg = sig;   % change in prediction 

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
% medpksep = round(median(diff(pkind)));
% fprintf('Med sep of peaks in sig-wlet CC is %d spls \n', medpksep);
%rcc serves the weight as the peak height, aka the master raw cc value at the peak
wtcoef = rcc(pkind).* pkhgt;
if isempty(fixthr)
  medwtcoef = median(wtcoef);  % median of the weighted master CC, could be percentile?
else
  medwtcoef = fixthr;
end
% medwtcoef
%%%%%%%%%%%%%%%%%%

%% iteration >= 1 
nsrc = size(impindep,1);
for isrc = 1: size(impindep,1)
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
    ampit(isrc, 7) = 0.5; %otherwise difference should be integer samples 
    dresit(isrc) = -1;
    mfitit(isrc) = -1;
    continue  % go to next source
  end
  %rcc serves the weight as the peak height, aka the master raw cc value at the peak
  wtcoef = [pkind rcc(pkind).* pkhgt]; % weighted master raw cc
  indabvthr = find(wtcoef(:,2)>=medwtcoef); %find weighted peaks above threshold
  if isempty(indabvthr)
    ampit(isrc, 1:6) = zeros(1,6);
    ampit(isrc, 7) = 0.5; %otherwise difference should be integer samples 
    dresit(isrc) = -1;
    mfitit(isrc) = -1;
    continue  % go to next source
  end
%   wtcoef = sortrows(wtcoef(indabvthr,:), 2, 'descend'); %focus on 'important' peaks only!
  wtcoef = wtcoef(indabvthr,:); %focus on 'important' peaks only!
  %%% Scheme 2, END %%%
  %%%%%%%%%%%%%% Different weighting schemes of CC %%%%%%%%%%%%%%

  %%%%%%%%%%%% Different ways to determine the 'best' peak %%%%%%%%%%%%%%%%
  %%% Way2: find all pairs that are close enough, choose the one with highest weighted CC
  %find which subwin this ref impulse belongs to
  if ~isempty(off1iw)
    iwin = findwhichrange(impindep(isrc,1),irccran);
    %use that subwin to compute the difference in alignment bewteen subwin and entire win
    offdiff = off1iw(iwin)-off1i; %just a single value
  else
    offdiff = 0;
  end
  
%   if ~isempty(off1iw)
%     iwin = findwhichrange(impindep(isrc,1),irccran);
%     %use that subwin to compute the difference in alignment bewteen subwin and entire win
%     offdiff = off1iw(iwin); %just a single value
%   else
%     offdiff = off1i;
%   end
  
  %source arrival prediction from plane fitting, including calibrating the arrival prediction in the
  %context of prealigned signal, note sign is '+' 
  [idxpred, off14_plfit] = pred_tarvl_at4thsta(stas,impindep(isrc,7),impindep(isrc,8),impindep(isrc,1),offdiff);
  %find peaks within allowable range, and height is positive??
  indclose = find(abs(wtcoef(:,1)-(idxpred-round(ldiff/2)))<=offmax & pkhgt(indabvthr)>=0); %context of indices of 'coefeff'
  if isempty(indclose)
    ampit(isrc, 1:6) = zeros(1,6);
    ampit(isrc, 7:8) = [0.5 0.5]; %otherwise difference should be integer samples 
    dresit(isrc) = -1;
    mfitit(isrc) = -1;
    continue
  else
    [~,imax] = max(wtcoef(indclose,2));  %choose the one with highest weighted CC
  end
  idximp = wtcoef(indclose(imax),1)+round(ldiff/2);  % convert it to index of the signal (or sigdecon)
  idxcoef = idximp-itwlet+1+nfft;  % index of the raw CC that gives the amp
  lagsamp = lag(idxcoef); % offset in samples, in the context of raw CC
  amp = coef(idxcoef)/sum(wlet.^2); % convert max raw CC to amplitude
  rccimp = rcc(wtcoef(indclose(imax),1)); % rcc value at the index of the selected impulse 
  %%% Way2: find all pairs that are close enough, choose the one with highest weighted CC
  %%%%%%%%%%%% Different ways to determine the 'best' peak %%%%%%%%%%%%%%%%

  %%%A plot to check the weighted CC of found eligible peaks
  if fpltchk
    figure
    ax = gca;
    hold(ax,'on');
    ax.Box='on'; grid(ax,'on');
    yyaxis(ax,'right');
    p2=plot(ax,((1+round(ldiff/2)): (nfft-round(ldiff/2)))/sps,rcc,'-','color',[.7 .7 .7],'markersize',2);  % rcc
%     scatter(ax,pkind+round(ldiff/2),rcc(pkind),20,'r');  % rcc value at peaks
    ylim(ax,[-1 1]);
    ylabel(ax,'RCC','FontSize',12);
    
    yyaxis(ax,'left');
    p1=plot(ax,((1+round(ldiff/2)): (nfft-round(ldiff/2)))/sps,coefeff,'-','color',[.3 .3 .3]); hold on % master cc
    scatter(ax,(pkind+round(ldiff/2))/sps,pkhgt,20,'r'); % peaks
%     scatter(idximp,pkhgt(mwtcoefidx),10,'b','filled');  % chosen peak
%     scatter(idximp,rcc(pkind(mwtcoefidx)),10,'b','filled'); % chosen peak
    p3=scatter(ax,(wtcoef(:,1)+round(ldiff/2))/sps,wtcoef(:,2),20,'b');  % peak heights weighted by rcc
    scatter(ax,idximp/sps,wtcoef(indclose(imax),2),20,'b','filled'); % chosen peak
    axsym(ax,2);
    p4=plot(ax,[idxpred idxpred]/sps,ax.YLim,'k-.','linew',1); % arrival prediction
    p5=plot(ax,ax.XLim, [medwtcoef medwtcoef],'k--'); % stopping threshold, subject to actual practice
%     p4=plot(ax,ax.XLim, [ampnoiwt ampnoiwt],'k--'); % noise level
%     p4=plot(ax.XLim, [ampnoi ampnoi],'k--'); % noise level
%     plot(ax.XLim,[mwtcoef mwtcoef],'b--');  % indeed the max?
    ylabel(ax,'Amplitude','FontSize',12);
%     xlabel(ax,sprintf('Samples at %d Hz',sps),'FontSize',12);
    xlabel(ax,'Time (s)','FontSize',12);
    title(sprintf('Source: %d',isrc));
%     xlim(ax,[1.5e4 2e4]);
    legend(ax,[p1,p2,p3,p4,p5],'Res-temp CC','RCC14','Wgted peaks above thresh','Arrival pred',...
      'threshold');
    hold(ax,'on');
  end

  %update the deconvolved impulse response array
  sigdecon(idximp) = sigdecon(idximp)+amp;
  ampit(isrc, 1:2) = [idximp amp];   %zero-crossing index and amp of decon impulse at this iter  
  [pkhres, pkires] = findpeaks(res); 
  [~,tempind] = min(abs(pkires-idximp));
  ampit(isrc, 3:4) = [pkires(tempind) pkhres(tempind)];  %closest residual waveform positive peak
  [pkhres, pkires] = findpeaks(-res); 
  [~,tempind] = min(abs(pkires-idximp));
  ampit(isrc, 5:6) = [pkires(tempind) pkhres(tempind)];  %closest residual waveform negative peak  
  ampit(isrc, 7) = idximp-(idxpred-round(ldiff/2));  %difference in samples between predicted arrival and found peak  
  [~,ampit(isrc, 8)] = pred_tarvl_at4thsta(stas,impindep(isrc,7),impindep(isrc,8),...
    impindep(isrc,1),offdiff);  %predicted off14
  impchg = zeros(nfft,1);   % the change of impulse resulting from the current iteration
  impchg(idximp) = amp;
  %change in predicted signal by convolving the change of impulse with wavelet
  predchg = conv(wlet, impchg, 'full'); 
  predchg = predchg(itwlet:nfft+itwlet-1);  % cut accordingly
  pred = pred + predchg;  % update the prediction
  
  res_new = sig-pred;   % residual
%   dres = abs((sum(res.^2)-sum(res_new.^2))./sum(res.^2));  % target evaluation objective
%   dres = abs(norm(res).^2-norm(res_new).^2)./norm(res).^2;  % this is equivalent to the above
  varred = abs(var(res)-var(res_new))./var(res);  % again equivalent to the above, but called variance reduction
%   dres = sum((res-res_new).^2)./sum(res.^2);  % alternative evaluation objective
%   mfit = sum(res_new.^2)./sum(sig.^2);  % misfit objective, we want the residual to decrease to 0
  mfit = norm(res_new);  % alternative misfit objective
  dresit(isrc) = varred*100;   % store them, make the relative change in percent, in terms of variance reduction
  mfitit(isrc) = mfit;
  
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
    title(sprintf('Source: %d',isrc));
    longticks(gca,2);
    if lwlet < lsig
      ax = gca;
      loc = getpos(ax);
      shrink(ax,lsig/lwlet,1);
      ax.Position(1)=loc(1);
      nolabels(ax,1);
    end

    subplot(812); 
    yyaxis left; p1=plot(sig, 'k-'); box on; grid on; ax=gca; yran = ax.YLim; xlim([0,nfft]);
    yyaxis right;
    p2=plot((1+round(ldiff/2)): (nfft-round(ldiff/2)),rcc,'-','color',[.8 .8 .8]);
    ylim([-1 1]);
    legend([p1,p2],'Signal','RCC14');
    nolabels(gca,1);
    longticks(gca,2);
    
    subplot(813);
    p1=plot(res, 'c-'); box on; grid on; xlim([0,nfft]); ylim(yran);
    legend(p1,'Residual bf. current iter.'); 
    nolabels(gca,1);
    longticks(gca,2);
    
    subplot(814);
    yyaxis left;
    p1=plot((1+round(ldiff/2)): (nfft-round(ldiff/2)),coefeff,'-','color',[.2 .2 .2]);
    box on; grid on; xlim([0,nfft]); hold on;
    ax = gca;
    scatter(pkind+round(ldiff/2),pkhgt,20,'r');
    plot([idximp idximp],ax.YLim,'k--');
    yyaxis right;
    p2=plot((1+round(ldiff/2)): (nfft-round(ldiff/2)),rcc,'-','color',[.8 .8 .8]);
    ylim([-1 1]);
    legend([p1,p2],'Res-temp CC','RCC14');
    nolabels(gca,1);
    longticks(gca,2); hold off
    
    subplot(815);
%     p1=stem(sigdecon, 'b','MarkerSize',4);
    imptemp = find(sigdecon>0); p1=stem(imptemp,sigdecon(imptemp),'b','MarkerSize',4);
    box on; grid on; hold on; xlim([0,nfft]);
    text(0.05,0.85,num2str(nimp),'FontSize',10,'Units','normalized');
    text(0.05,0.65,num2str(idximp),'FontSize',10,'Units','normalized');
    legend(p1,'Found impulses'); 
    nolabels(gca,1);
    longticks(gca,2);
    
    subplot(816);    
    p1=plot(pred, 'b-');  box on; grid on; xlim([0,nfft]); ylim(yran);
    legend(p1,'Prediction');
    nolabels(gca,1);
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
    legend(p1,'Residual af. current iter.');
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

%% Post-processing
sigdecon = sigdecon(1:lsig);  %cut deconvolved signal

%make output a column vector
dresit = dresit(:);
mfitit = mfitit(:);
  
%Disp information if some stopping criteria have been met 
if ~isempty(ampit)
  fprintf('%d of %d sources have been checked \n',sum(sum(ampit(:,1:6),2)~=0),nsrc);
else
  fprintf('%d of %d sources have been checked \n',nan,nan);
end

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
  title(ax,sprintf('%d of %d sources passed the check',sum(sum(ampit(:,1:6),2)~=0),nsrc));
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
  text(ax,0.05,0.9,num2str(nimp),'FontSize',12,'Units','normalized');
  axranexp(ax,2,10);  ax.Box='on'; grid(ax,'on');
  yyaxis(ax,'right');
  p2=plot(ax,(1+round(ldiff/2)): (nfft-round(ldiff/2)),rcc,'-','color',[.6 .6 .6]);
  ylim(ax,[-1 1]);
  legend(ax,[p1,p2],'Found impulses','RCC14'); xlim(ax,[0,lsig]);
  longticks(ax,2); nolabels(ax,1); hold(ax,'off');
  
  ax=f2.ax(3); 
  yyaxis(ax,'left');
  p1=plot(ax,sig, 'k-'); xlim(ax,[0,lsig]); hold(ax,'on'); 
  p2=plot(ax,pred, 'b-'); yran=1.2*[-max(abs([sig; pred])) max(abs([sig; pred]))]; 
  ylim(ax,yran); ax.Box='on'; grid(ax,'on');
  yyaxis(ax,'right');
  mwlen = sps/2;
  [irccsp, rccsp] = RunningCC(sig, pred, mwlen);
  p3=plot(ax,irccsp,rccsp,'-','color',[.6 .6 .6]);
  ylim(ax,[-1 1]);
  legend(ax,[p1,p2,p3],'Signal','Prediction','Sig-pred RCC','location','southeast');
  longticks(ax,2); nolabels(ax,1); hold(ax,'off');
  
  ax=f2.ax(4); 
  p1=plot(ax,res, 'c-'); xlim(ax,[0,lsig]); ylim(ax,yran); hold(ax,'on');
  legend(ax,[p1],'Residual');
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
  ind = find(dresit~=-1 & mfitit~=-1);
  p1=scatter(ax,ind,dresit(ind),20,'b');
  scatter(ax,ind(end),dresit(ind(end)),20,'b','filled');
  text(ax,round(ind(end)*4/5), dresit(ind(end))+0.15, num2str(dresit(ind(end))));
  ylabel(ax,'Variance reduction (%)');
  yyaxis(ax,'right');
  p2=scatter(ax,ind,mfitit(ind),20,'r');
  scatter(ax,ind(end),mfitit(ind(end)),20,'r','filled');
  text(ax,round(ind(end)*4/5), mfitit(ind(end))+0.15, num2str(mfitit(ind(end))));
  ylabel(ax,'L2-norm');
  legend(ax,[p1,p2],'Variance reduction of residual (dres)','L2-norm of residual, or data misfit (mfit)');
  xlabel(ax,'Source number');
  hold(ax,'off');
  fighdl{3} = f3;
  
end  

% keyboard



