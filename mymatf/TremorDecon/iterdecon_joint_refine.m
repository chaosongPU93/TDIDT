function [sigdecon,pred,res,dresit,mfitit,ampit,nit,fighdl] = ...
  iterdecon_joint_refine(sig,wlet,noi,impraw,dt,twlet,loff_max,dres_min,mfit_min,nit_max,...
  nimp_max,fpltit,fpltend,fpltchk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [sigdecon,pred,res,dresit,mfitit,nit,rf] = ...
%   iterdecon_joint_refine(sig,wlet,dt,twlet,width,dres_min,mfit_min,nit_max,fpltit,fpltend)
%
% Similar to 'iterdecon_joint', this function is to carry out the iterative deconvolution after 
% Ligorria & Ammon (1999), Essentially, wlet will be deconvolved from sig. The difference
% is that, in case the current order of carrying out the deconvolution, i.e., which impulse to 
% deconvolve first, that is determined by the master CC between the residual and template, 
% which is then weighted by the running CC between all stations, can lead to a bias to the 
% results, we now carry out the deconvolution starting from the first arrival in the result
% of 'iterdecon_joint'. This issue might exist when your template is not symmetric, and the 
% left-hand side of the main dipole is likely to be noise, while the right-hand side has coda.
% Therefore, the coda of earlier impulses would have impact on the arrival that should be 
% deconvolved based on the weighted master CC
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
% Chao Song, chaosong@princeton.edu
% First created date:   2022/01/27
% Last modified date:   2022/01/27
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
offbit = [];  % offset between impulses at different stations for each pair/iteration
nimp = 0; % total number of impulses 
predchg = sig;   % change in prediction 
amp = ones(1,nsta);   % intial amplitude of impulse

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
ampnoinorm = ampnoi./sum(wlet.^2, 1);   % normalized amp of CC between noise and wavelet 

%average running CC between signals 
cclen=round(1/dt)/2;
[irccs,rccs12] = RunningCC(sig(:,1), sig(:,2), cclen);
[~,rccs13] = RunningCC(sig(:,1), sig(:,3), cclen);
[~,rccs23] = RunningCC(sig(:,2), sig(:,3), cclen);
rccs = (rccs12+rccs13+rccs23)/3;

%% iteration
%%%This just loops over all deconvolved impulses resulted from 'iterdecon_joint_refine', imagine 
%%%now you have a series of raw impulses with their indices and amplitudes\
impraw = sortrows(impraw,1);  % sort the impulse chronologically 
imprawidx = median(impraw(:,[1 3 5]),2);
for iimp = 1: size(imprawidx,1)
  
  nit = nit+1;

  %compute the running CC between residuals at all stations, at iter 0, use signals
  cclen=round(1/dt)/2;
  [ircc,rcc12] = RunningCC(res(:,1), res(:,2), cclen);
  [~,rcc13] = RunningCC(res(:,1), res(:,3), cclen);
  [~,rcc23] = RunningCC(res(:,2), res(:,3), cclen);
  rcc = (rcc12+rcc13+rcc23)/3;

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

  %%% For each predetermined raw impulse, we try to find the closest peak in the master raw CC
  %%% between the residual and wavelet at each station as the location to deconvolve the residual
  %%% first, so that theoretically we should find the triplets of peaks 
  for i = 1: nsta
    %effective raw CC that corresponds to the index of the overlapping portion between signal (or 
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
      pkssta{nit}{i} = pks; % pks, [ind_of_peak, peak_height(master cc value)]
      if size(impraw,2)>1
        off = pks(:,1)+round(cclen/2)-impraw(iimp,(i-1)*2+1);
      else
        off = pks(:,1)+round(cclen/2)-impraw(iimp,1);
      end
      indclose = find(abs(off)==min(abs(off(:))));
      if length(indclose)>1  % if there are multiple combinations whose offsets are equally close
        [~,ind] = max(pks(indclose,2));   % choose the pair has the highest peak height, ie, master cc 
        bsub(i) = indclose(ind);
      else
        bsub(i) = indclose;
      end
      idximp(i) = pks(bsub(i),1)+round(cclen/2);
      idxcoef = idximp(i)-itwlet(i)+1+nfft;
      lagsamp(i) = lag(idxcoef);
      amp(i) = coef(idxcoef, i)/sum(wlet(:, i).^2);
%       amp2(i) = pks(bsub(i),2)/sum(wlet(:, i).^2); % this should give the same result as above      
    end    
  end
  
  off12b = idximp(1)-idximp(2);
  off13b = idximp(1)-idximp(3);
  off23b = idximp(2)-idximp(3);
  if abs(off12b)>loff_max || abs(off13b)>loff_max || abs(off23b)>loff_max || ...
     (abs(off12b)+abs(off13b)+abs(off23b))/3>loff_max
%     nit = nit-1;  % -1 to recognize the current iteration does not proceed to the end
%     fprintf('Forced stop at iteration %d \n', nit);
%     fprintf('The offsets of best solution are still larger than allowed \n');
%     fprintf('at the next iteration \n');
    fprintf('off12= %d, off13= %d, off23= %d \n',off12b, off13b, off23b);
%     break
  else
    fprintf('off12= %d, off13= %d, off23= %d \n',off12b, off13b, off23b);
  end
   
  %%%A plot to check the weighted CC of found eligible peaks 
  if fpltchk
    figure
    plot(ircc,coefeff(:,1),'r-'); hold on
    plot(ircc,coefeff(:,2),'b-');
    plot(ircc,coefeff(:,3),'k-');
    ax = gca;
    if size(impraw,2)>1
      plot([impraw(iimp,1) impraw(iimp,1)],ax.YLim,'r--');
      plot([impraw(iimp,3) impraw(iimp,3)],ax.YLim,'b--');
      plot([impraw(iimp,5) impraw(iimp,5)],ax.YLim,'k--');
    else
      plot([impraw(iimp,1) impraw(iimp,1)],ax.YLim,'--','color',[.5 .5 .5]);
    end
    pks = pkssta{nit}{1}; scatter(pks(:,1)+round(cclen/2),pks(:,2),20,'r');
    scatter(idximp(1),pks(bsub(1),2),10,'r','filled');
    pks = pkssta{nit}{2}; scatter(pks(:,1)+round(cclen/2),pks(:,2),20,'b');
    scatter(idximp(2),pks(bsub(2),2),10,'b','filled');
    pks = pkssta{nit}{3}; scatter(pks(:,1)+round(cclen/2),pks(:,2),20,'k');
    scatter(idximp(3),pks(bsub(3),2),10,'k','filled');
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
  %   mfit = sum(res_new.^2)./sum(sig.^2);  % misfit objective, we want the residual to decrease to 0
    mfit(i) = norm(res_new(:,i));  % alternative misfit objective
  end  
  ampit(nit, nsta*2+1:nsta*2+3) = [off12b off13b off23b];
  
  %%%may have a different choice here, do we sum all stations here?
  dresit(nit) = sum(dres*100)/nsta;   % store them, make the relative change in percent
  mfitit(nit) = sum(mfit)/nsta;

  %AIC criterion, eq 28 of Ohta&Ide2017JGRSE
  %this can be used only when there will be a global minimum and pre-compute for a large iterations
  dof = nsta*nit;  % number of free params, for each iteration nit, only the index time idximp is free
  if var(noi(:,1))>0
    AIC(nit) = sum(mfit.^2)/sum(var(noi))+2*dof; % variance of background noise 
  else
    AIC(nit) = 0;
  end
    
  nimp = sum(sigdecon>0, 1); % number of non-zero impulses
  
  %%%Updating plot for each iteration
  if fpltit   % if the flag for plotting the result for each iteration is on 
    f1.fig = figure(101); clf(f1.fig)
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 9;   % maximum height allowed is 11 inches
    [scrsz, resol] = pixelperinch(1);
    set(f1.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);
    
    subplot(811); box on; grid on; hold on
    p1=plot(wlet(:,1),'r-'); plot(wlet(:,2),'b-'); plot(wlet(:,3),'k-');
    text(0.95,0.85,'Template','Units','normalized','HorizontalAlignment','right','FontSize',10);
    xlim([0,lwlet]);
    title(sprintf('Iteration: %d',nit)); hold off
    
    subplot(812); box on; grid on; hold on
    p1=plot(sig(:,1), 'r-'); plot(sig(:,2), 'b-'); plot(sig(:,3), 'k-'); 
    ax=gca; yran = [-max(max(abs(sig))) max(max(abs(sig)))];
    text(0.95,0.85,'Signal','Units','normalized','HorizontalAlignment','right','FontSize',10);
    xlim([0,nfft]); ylim(yran); hold off
    
    subplot(813); box on; grid on; hold on
    yyaxis left; p1=plot(res(:,1), 'r-');
    plot(res(:,2), 'b-'); plot(res(:,3), 'k-'); ylim(yran);
    yyaxis right;
    p2=plot(ircc,rcc,'o','color',[.8 .8 .8],'markersize',2);
    ylim([-1 1]);
    text(0.95,0.85,'Residual bf. current iter.','Units','normalized','HorizontalAlignment',...
      'right','FontSize',10);
    xlim([0,nfft]);  hold off
        
    subplot(814); box on; grid on; hold on;
    stem(impraw(1:iimp,1),impraw(1:iimp,2), 'r','MarkerSize',4);
    stem(impraw(1:iimp,3),impraw(1:iimp,4), 'b','MarkerSize',4);
    stem(impraw(1:iimp,5),impraw(1:iimp,6), 'k','MarkerSize',4);    
    text(0.05,0.85,num2str(length(unique(impraw(1:iimp,1)))),'color','r','FontSize',10,'Units','normalized');
    text(0.1,0.85,num2str(length(unique(impraw(1:iimp,3)))),'color','b','FontSize',10,'Units','normalized');
    text(0.15,0.85,num2str(length(unique(impraw(1:iimp,5)))),'color','k','FontSize',10,'Units','normalized');
    text(0.05,0.65,num2str(impraw(iimp,1)),'color','r','FontSize',10,'Units','normalized');
    text(0.1,0.65,num2str(impraw(iimp,7)),'color','b','FontSize',10,'Units','normalized');
    text(0.15,0.65,num2str(impraw(iimp,8)),'color','k','FontSize',10,'Units','normalized');
    text(0.95,0.85,'Old impulses','Units','normalized','HorizontalAlignment',...
      'right','FontSize',10); xlim([0,nfft]); hold off
    
    subplot(815); box on; grid on; hold on;
    imptemp = find(sigdecon(:,1)>0); p1=stem(imptemp,sigdecon(imptemp,1), 'r','MarkerSize',4);
    text(0.05,0.9,num2str(length(imptemp)),'color','r','FontSize',10,'Units','normalized');
    imptemp = find(sigdecon(:,2)>0); stem(imptemp,sigdecon(imptemp,2), 'b','MarkerSize',4);
    text(0.1,0.9,num2str(length(imptemp)),'color','b','FontSize',10,'Units','normalized');
    imptemp = find(sigdecon(:,3)>0); stem(imptemp,sigdecon(imptemp,3), 'k','MarkerSize',4);
    text(0.15,0.9,num2str(length(imptemp)),'color','k','FontSize',10,'Units','normalized');
    text(0.05,0.7,num2str(idximp(1)),'color','r','FontSize',10,'Units','normalized');
    text(0.1,0.7,num2str(idximp(1)-idximp(2)),'color','b','FontSize',10,'Units','normalized');
    text(0.15,0.7,num2str(idximp(1)-idximp(3)),'color','k','FontSize',10,'Units','normalized');
    plot([0 nfft],[ampnoinorm(1) ampnoinorm(1)],'r--');
    plot([0 nfft],[ampnoinorm(2) ampnoinorm(2)],'b--');
    plot([0 nfft],[ampnoinorm(3) ampnoinorm(3)],'k--');
    text(0.95,0.85,'New impulses','Units','normalized','HorizontalAlignment',...
      'right','FontSize',10); xlim([0,nfft]); hold off
    
    subplot(816); box on; grid on; hold on   
    p1=plot(pred(:,1), 'r-'); plot(pred(:,2), 'b-'); plot(pred(:,3), 'k-'); 
    text(0.95,0.85,'Prediction','Units','normalized','HorizontalAlignment',...
      'right','FontSize',10);
    xlim([0,nfft]); ylim(yran); hold off
    
    subplot(817); box on; grid on; hold on;   
    yyaxis left; p1=plot(res_new(:,1), 'r-'); 
    plot(res_new(:,2), 'b-'); plot(res_new(:,3), 'k-'); ylim(yran);
    [irccn,rccn12] = RunningCC(res_new(:,1), res_new(:,2), cclen);
    [~,rccn13] = RunningCC(res_new(:,1), res_new(:,3), cclen);
    [~,rccn23] = RunningCC(res_new(:,2), res_new(:,3), cclen);
    rccn = (rccn12+rccn13+rccn23)/3;
    yyaxis right; plot(irccn,rccn,'o','color',[.8 .8 .8],'markersize',2);
    ylim([-1 1]);
    text(0.95,0.85,'Residual af. current iter.','Units','normalized','HorizontalAlignment',...
      'right','FontSize',10);
    xlim([0,nfft]); hold off
    
    subplot(818); box on; grid on; hold on;   
    p1=plot(noi(:,1), 'r-'); plot(noi(:,2), 'b-'); plot(noi(:,3), 'k-');
    text(0.95,0.85,'Noise','Units','normalized','HorizontalAlignment',...
      'right','FontSize',10);
    xlim([0, max(nfft, size(noi,1))]); ylim(yran); hold off
    
    print(f1.fig,'-djpeg','-r300',strcat('/home/chaosong/Pictures/RefineIter',num2str(nit),'.jpg'));

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
  if sum(amp>ampnoinorm)==0
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
  f2.fig = figure(201);clf(f2.fig);
  widin = 10;  % maximum width allowed is 8.5 inches
  htin = 9;   % maximum height allowed is 11 inches
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
  text(ax,0.95,0.85,'Template','Units','normalized','HorizontalAlignment','right','FontSize',10);
  title(ax,sprintf('Total Iterations: %d',nit)); hold(ax,'off');
    
  ax=f2.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); xlim(ax,[0,lsig]);
  plot(ax,[0 nfft],[ampnoinorm(1) ampnoinorm(1)],'r--');
  plot(ax,[0 nfft],[ampnoinorm(2) ampnoinorm(2)],'b--');
  plot(ax,[0 nfft],[ampnoinorm(3) ampnoinorm(3)],'k--');
  imptemp = find(sigdecon(:,1)>0); p1=stem(ax,imptemp,sigdecon(imptemp,1), 'r','MarkerSize',4); 
  text(ax,0.05,0.85,num2str(length(imptemp)),'color','r','FontSize',10,'Units','normalized');
  imptemp = find(sigdecon(:,2)>0); stem(ax,imptemp,sigdecon(imptemp,2), 'b','MarkerSize',4);
  text(ax,0.1,0.85,num2str(length(imptemp)),'color','b','FontSize',10,'Units','normalized');
  imptemp = find(sigdecon(:,3)>0); stem(ax,imptemp,sigdecon(imptemp,3), 'k','MarkerSize',4);
  text(ax,0.15,0.85,num2str(length(imptemp)),'color','k','FontSize',10,'Units','normalized');
  text(ax, 0.95,0.85,'Impulses','Units','normalized','HorizontalAlignment',...
    'right','FontSize',10);  hold(ax,'off');

  yran = [-max(max(abs(sig))) max(max(abs(sig)))];
  ax=f2.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); xlim(ax,[0,lsig]); ylim(ax,yran); 
  yyaxis(ax,'left');
  p1=plot(ax,sig(:,1),'r-'); plot(ax,sig(:,2),'b-'); plot(ax,sig(:,3),'k-');
  yyaxis(ax,'right');
  [irccs,rccs12] = RunningCC(sig(:,1), sig(:,2), cclen);
  [~,rccs13] = RunningCC(sig(:,1), sig(:,3), cclen);
  [~,rccs23] = RunningCC(sig(:,2), sig(:,3), cclen);
  rccs = (rccs12+rccs13+rccs23)/3;
  plot(ax,irccs,rccs,'o','color',[.8 .8 .8],'markersize',2);
  ylim(ax,[-1 1]);
  text(ax, 0.95,0.85,'Signal','Units','normalized','HorizontalAlignment',...
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
  yyaxis(ax,'right');
  [irccr,rccr12] = RunningCC(res(:,1), res(:,2), cclen);
  [~,rccr13] = RunningCC(res(:,1), res(:,3), cclen);
  [~,rccr23] = RunningCC(res(:,2), res(:,3), cclen);
  rccr = (rccr12+rccr13+rccr23)/3;
  plot(ax,irccr,rccr,'o','color',[.8 .8 .8],'markersize',2);
  ylim(ax,[-1 1]);
  text(ax, 0.95,0.85,'Residual','Units','normalized','HorizontalAlignment',...
    'right','FontSize',10); 
  xlabel(ax,'Samples','FontSize',12);  hold(ax,'off');
  fighdl{2} = f2;

  f3.fig = figure(301); clf(f3.fig);
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
  
  f4.fig = figure(401); clf(f4.fig);
  plot(log(AIC)); hold on;
  ind = find(AIC==min(AIC)); scatter(ind,log(AIC(ind)),'ro');
  text(0.4,0.5,sprintf('Optimal num. of iterations (min. AIC): %d', ind),'Units','normalized',...
    'HorizontalAlignment','center');
  xlabel('Iteration number');
  ylabel('Log_{10}(AIC)');
  fighdl{4} = f4;
end  


% keyboard









