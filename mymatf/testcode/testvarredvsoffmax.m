function [nsrctest,nsrcnorm,rr,vr] = ...
  testvarredvsoffmax(ista,sigsta,optdat,greenf,sps,zcrosses,planefit,irccran,...
  rcc1icat,impindep,stas,off1i,off1iw,overshoot,mswccpksep,pltflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [propang,locxyprop,fitobj,gof,output] = propadirection(timevec,locxy,angleopt)
%
% This function is to first compute the propagation direction events that
% are represented by 'timevec' and 'locxy'. Locations are projected along
% a series of trial directions. A Bisquare robust linear fit is implemented
% to the projected locations. The direction along which, either the RMSE is
% the minimum, or the slope is the maximum, is chosen as the propagation
% direction. Option is made by 'angleopt'. Projected locations and the fit
% statistics along the propagation direction is also returned.
% See 'srcprojdistNtoNm.m' as well.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/01/18
% Last modified date:   2024/01/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('pltflag',1);

wlet = greenf(:,ista);  %template here is best aligned, tapered, linear trend removed, filtered
sig = sigsta(:,ista); %best aligned, filtered, tapered

dt = 1/sps;  % sampling interval
twlet = zcrosses(ista)*dt;
fpltit = 0;  % plot flag for each iteration
fpltend = 0;  % plot flag for the final iteration
fpltchk = 0; % plot flag for intermediate computations
rmse = planefit.gof{ista-3}.rmse;
%+-1-sigma ~0.68, 1.5-sigma ~0.87, 2-sigma ~0.95, 3-sigma ~0.997
%         sigma1 = round(rmse);
%         sigma2 = round(2*rmse);
offmaxtry = 4:1:40;
for ioff = 1: length(offmaxtry)
  offmax = offmaxtry(ioff);
  if isempty(off1iw)
    [~,~,~,~,~,ampiti,~] = iterdecon_4thsta(sig,wlet,...
      irccran,rcc1icat(:,ista-3),[],dt,twlet,impindep,stas(ista,:),...
      off1i(ista),[],offmax,fpltit,fpltend,fpltchk);
  else
    [~,~,~,~,~,ampiti,~] = iterdecon_4thsta(sig,wlet,...
      irccran,rcc1icat(:,ista-3),[],dt,twlet,impindep,stas(ista,:),...
      off1i(ista),off1iw(:,ista),offmax,fpltit,fpltend,fpltchk);
  end
    
  imptemp = impindep;
  imptemp(:,9+(ista-4)*2+1) = ampiti(:,1);
  imptemp(:,9+(ista-3)*2) = ampiti(:,2);
  
  %%% further eliminate sources that fail the check at 4th stations
  trust4th = 7; % trust KLNB the most among all 4th stations
  indremove = find(imptemp(:,9+(trust4th-4)*2+1)==0 & imptemp(:,9+(trust4th-3)*2)==0);
  imptemp(indremove,:) = [];
  imptempst = sortrows(imptemp,1);
  nsrctest(ioff,1) = size(imptempst,1);
  nsrcnorm(ioff,1) = size(imptempst,1)/ size(impindep,1)* 100;  %normalized by # after removing 2ndary srcs
  
  %%%final prediction via convolution between grouped impulses and template at each station
  [~,~,~,~,~,l2normred1,varred1]=...
    predsig_conv_imptemp(sigsta,optdat,imptempst,...
    greenf,zcrosses,overshoot,stas,0,sps);
  rr(ioff,1) = l2normred1(trust4th,3);  %residual reduction at KLNB
  rr(ioff,2) = mean(l2normred1(1:3,3)); %mean residual reduction at trio stas
  vr(ioff,1) = varred1(trust4th,3); %variance reduction at KLNB
  vr(ioff,2) = mean(varred1(1:3,3));%mean variance reduction at trio stas
end

%% plotting
if pltflag
  nrow = 1; % rows and cols of subplots in each figure
  ncol = 1;
  widin = 4; % size of each figure
  htin = 4.5;
  pltxran = [0.12 0.88]; pltyran = [0.12 0.96];
  pltxsep = 0.04; pltysep = 0.03;
  f = initfig(widin,htin,nrow,ncol);
  optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
  ax = f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  yyaxis(ax,'left');
  p0=plot(ax,offmaxtry/rmse,nsrcnorm,'b-','linew',1.5);
  p1=plot(ax,[2 2], ax.YLim,'k--','linew',1.5);
  mseppks=mswccpksep(ista);
  p2=plot(ax,[mseppks mseppks]/2/rmse, ax.YLim,'k-.','linew',1);
  xlabel(ax,sprintf('Allowed difference from predicted arrival (x RMSE)'));
  ylabel(ax,'Preserved sources (%)');
  yyaxis(ax,'right');
  p3=plot(ax,offmaxtry/rmse,vr(:,1),'r-','linew',1.5);
  p4=plot(ax,offmaxtry/rmse,vr(:,2),':','linew',1,'Color',[1 .45 0]); %orange
  text(ax,0.98,0.85,sprintf('%s',strtrim(stas(ista,:))),...
    'Units','normalized','HorizontalAlignment','right','fontsize',12);
  xticks(ax,0:0.5:6);
  ylabel(ax,'Variance reduction (%)');
  legend(ax,[p1,p2,p4],'Chosen threshold',...
    sprintf('Half of median separation \nof peaks in sig-temp CC'),...
    'Mean VR at trio stations','Location','southeast','fontsize',8);
  %       ylim(ax,[16 34]);

%   keyboard
  if isempty(off1iw)
    fname = 'varredvsoffmax1win_no23.pdf';
  else
    fname = 'varredvsoffmax_no23.pdf';
  end
  print(f.fig,'-dpdf',...
    strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));

  %       ax = subplot(1,2,2);
  %       hold(ax,'on');
  %       ax.Box='on'; grid(ax,'on');
  %       scatter(ax,nsrcnorm(k,:),vr(k,:,1),25,'k','filled');
  %       xlabel('# sources left (%)');
  %       ylabel('Change in L-2 norm (%)');
  %       ylim(ax,[16 34]);
end




