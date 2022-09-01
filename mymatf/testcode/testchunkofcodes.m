%%%%%%%%%%%%%%%%%% this part is for testing  %%%%%%%
tst = 35810; % this range is for testing
ted = 35850;
ist = floor(tst*sps+1);
ied = floor(ted*sps);
hfobj = hf(hf(:,seccol)>=tst & hf(:,seccol)<=ted, :); % object hf tremor detections

PREALIGN = 'median';
if strcmp(PREALIGN, 'median')
  off12med = median(hfobj(:, 11))*sps; % median offset from detections
  off13med = median(hfobj(:, 12))*sps;
  off12 = round(off12med);
  off13 = round(off13med);
%   off12 = 1;
%   off13 = -1;
  
  %should account for time offsets between stations here
  optdat = [];
  ortdat = [];
  optdat(:, 1) = STAopt(max(ist,1): min(ied,86400*sps), 2);
  ortdat(:, 1) = STAort(max(ist,1): min(ied,86400*sps), 2);
  optdat(:, 2) = STAopt(max(ist-off12,1): min(ied-off12,86400*sps), 3);
  ortdat(:, 2) = STAort(max(ist-off12,1): min(ied-off12,86400*sps), 3);
  optdat(:, 3) = STAopt(max(ist-off13,1): min(ied-off13,86400*sps), 4);
  ortdat(:, 3) = STAort(max(ist-off13,1): min(ied-off13,86400*sps), 4);

elseif strcmp(PREALIGN, 'constrained')
  %maybe a median offset for a long migration is too rough, better to separate into shorter segments
  %and plot, align them separately, and compute the running CC
  %also, the median offset may not be enclosed
  %%%%%%%%%%% independent enclosed alignment %%%%%%%%%%%%%%%%%%
  %first align the broadband trace
  mshiftadd = sps/8;    % maximum allowed shift between 2 traces
  buffer = 2*mshiftadd;
  opttmp = STAopt(max(ist-buffer,1): min(ied+buffer,86400*sps), 2:end); %some buffer for shifting
  orttmp = STAort(max(ist-buffer,1): min(ied+buffer,86400*sps), 2:end); %some buffer for shifting
  mid = ceil(size(opttmp,1)/2);
  fixlen = ied-ist+1;
  loffmax = 5;
  ccmin = 0.01;  % depending on the length of trace, cc could be very low
  iup = 1;    % times of upsampling
  [off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(opttmp',mid,...
    fixlen,mshiftadd,loffmax,ccmin,iup);
  off12 = round(off12con);
  off13 = round(off13con);
  
  % align the records
  istart = buffer+1;
  iend = istart+ied-ist;
  optdat = [];
  ortdat = [];
  optdat(:, 1) = opttmp(istart: iend, 1);
  ortdat(:, 1) = orttmp(istart: iend, 1);
  optdat(:, 2) = opttmp(istart-round(off12): iend-round(off12), 2);
  ortdat(:, 2) = orttmp(istart-round(off12): iend-round(off12), 2);
  optdat(:, 3) = opttmp(istart-round(off13): iend-round(off13), 3);
  ortdat(:, 3) = orttmp(istart-round(off13): iend-round(off13), 3);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

stasig = [];
stawlet = [];
for ista = 1: 3

  sigbb = optdat(:,ista);
  % sigbb = sigbb(1:50*sps);
  lsig = length(sigbb);
  
%   figure
%   plot(tracetemp(:, ista));
  
  lwlet = pow2(9);
  zerocs = 1198;
  %according analysis to the template, the zero-crossing to the end of coda is safe to set as 8 s
  wletbb = tracetemp(zerocs+8*sps-lwlet+1: zerocs+8*sps, ista);
  
  %% remove mean and linear trend and taper
  %of signal
  sigbbd = detrend(sigbb);
  %and taper with tukeywin, which is actually a tapered cosine window
  %tapered length is adaptative to frequency, maybe at least longer than one full period length of
  %the lowest frequency
  fractap = sps*4/size(sigbbd,1); % if fractap is >=1, n-point von Hann window is returned
  ptstap = fractap/2*size(sigbbd,1); % points tapered at the start/end 
  w = tukeywin(size(sigbbd,1),fractap);
  sigbbdt = w.* sigbbd;
  %detrend again for caution
  sigbbdt=detrend(sigbbdt);

  %romve mean, linear trend of template
  wletbbd = detrend(wletbb);
  %and taper with tukeywin, which is actually a tapered cosine window
  w = tukeywin(size(wletbbd,1),fractap);
  wletbbdt = w.* wletbbd;
  %detrend again for caution
  wletbbdt=detrend(wletbbdt);
  
  %% filter the signal and template to reach a similar spectra shape
  %%% 1.8-4.5 for signal and 1.8-6.0 for template seem like a good pair
  %%% OR, 1.8-5.4 for signal and 1.8-18 for template to keep the most of the template at the high end
  %%% OR, use 1.8-4.5 for both, but use more poles (eg., 3) for signal to achieve a faster decay
  %filter the signal
  hisig=5.4;
  losig=1.8;
  sig = Bandpass(sigbbdt, sps, losig, hisig, npo, npa, 'butter');
  %detrend again for caution
  sig=detrend(sig);
  
  %filter the template
  hiwlet=18;
  lowlet=1.8;
  wlet = Bandpass(wletbbdt, sps, lowlet, hiwlet, npo, npa, 'butter');
  %detrend again for caution
  wlet=detrend(wlet);

  %store them for each station
  stasig(:,ista) = sig;
  stawlet(:,ista) = wlet;

end

%% get the running average CC between 2 stations
cclen=sps*4;
[ircc,rcc12] = RunningCC(stasig(:,1), stasig(:,2), cclen);
[~,rcc13] = RunningCC(stasig(:,1), stasig(:,3), cclen);
[~,rcc23] = RunningCC(stasig(:,2), stasig(:,3), cclen);
rcc = (rcc12+rcc13+rcc23)/3;

%% convolve the original signal with template of other station, and get the running average of it
dt = 1/sps;  % sampling interval
twlet = [190/sps; 191/sps; 191/sps];  % Time shift of main arrival of wavelet
s1w2 = conv(stasig(:,1),stawlet(:,2),'full'); % the order doesn't matter
s1w2 = s1w2(1+round(twlet(2)/dt):lsig+round(twlet(2)/dt));  % cut accordingly
s1w3 = conv(stasig(:,1),stawlet(:,3),'full');
s1w3 = s1w3(1+round(twlet(3)/dt):lsig+round(twlet(3)/dt));
s2w1 = conv(stasig(:,2),stawlet(:,1),'full');
s2w1 = s2w1(1+round(twlet(1)/dt):lsig+round(twlet(1)/dt));  % cut accordingly
s2w3 = conv(stasig(:,2),stawlet(:,3),'full');
s2w3 = s2w3(1+round(twlet(3)/dt):lsig+round(twlet(3)/dt));
s3w1 = conv(stasig(:,3),stawlet(:,1),'full');
s3w1 = s3w1(1+round(twlet(1)/dt):lsig+round(twlet(1)/dt));  % cut accordingly
s3w2 = conv(stasig(:,3),stawlet(:,2),'full');
s3w2 = s3w2(1+round(twlet(2)/dt):lsig+round(twlet(2)/dt));  % cut accordingly
cclen=sps*4;
[ircccon,rcc12c12] = RunningCC(s1w2, s2w1, cclen);  % running cc of 12 between conv of 1 and 2
[~,rcc13c13] = RunningCC(s1w3, s3w1, cclen);
[~,rcc23c23] = RunningCC(s2w3, s3w2, cclen);
rcccon = (rcc12c12+rcc13c13+rcc23c23)/3;

%%%If you want to compare between 3 stations, you have to convolve the record at one station with 
%%%both templates at the other 2 stations
s1w23 = conv(s1w2,stawlet(:,3),'full'); 
s1w23 = s1w23(1+round(twlet(3)/dt):lsig+round(twlet(3)/dt));  % cut accordingly
s2w13 = conv(s2w1,stawlet(:,3),'full');
s2w13 = s2w13(1+round(twlet(3)/dt):lsig+round(twlet(3)/dt));  % cut accordingly
s3w12 = conv(s3w1,stawlet(:,2),'full');
s3w12 = s3w12(1+round(twlet(2)/dt):lsig+round(twlet(2)/dt));  % cut accordingly
cclen=sps*4;
[irccc123,rcc12c123] = RunningCC(s1w23, s2w13, cclen);
[~,rcc13c123] = RunningCC(s1w23, s3w12, cclen);
[~,rcc23c123] = RunningCC(s2w13, s3w12, cclen);
rccc123 = (rcc12c123+rcc13c123+rcc23c123)/3;


%% plot the processed signal in time domain in detail
f.fig = figure;
f.fig.Renderer = 'painters';
widin = 18;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);

ncol = 1;
nrow = 2*3;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
%     grid(f.ax(isub),'on');
end

pltxran = [0.03 0.95]; pltyran = [0.05 0.95];
pltxsep = 0.02; pltysep = 0.04; 
axpos = optaxposition(nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
for isub = 1:nrow*ncol
  set(f.ax(isub), 'position', axpos(isub,:));
end

timemax = (hfobj(:, seccol) -tst)*sps +1; % time of max power rate of half sec
timecnt = (hfobj(:, seccol-1) -tst)*sps +1; % time of center of the 4-s detecting window
off12all = round(hfobj(:, 11)*sps); % offsets of 4-s tremor detections
off13all = round(hfobj(:, 12)*sps);

timearm = round((armobj(:,4) -tst)*sps +1); % time of armbruster's tremor detections
timebost = round((bostobj(:,5) -tst)*sps +1); % time of armbruster's tremor detections

for ista = 1: 3
  isub = (ista-1)*2+1;
  ax = f.ax(isub);
  hold(ax,'on');
  yyaxis(ax,'left');
  if ista == 1
    plot(ax,stasig(:,2),'b-');
    plot(ax,stasig(:,3),'k-'); %,'linew',0.5
  elseif ista == 2
    plot(ax,stasig(:,1),'r-');
    plot(ax,stasig(:,3),'k-'); %,'linew',0.5
  else
    plot(ax,stasig(:,1),'r-');
    plot(ax,stasig(:,2),'b-');
  end
  text(ax,0.05,0.95,'Original','unit','normalized',...
    'HorizontalAlignment','left','FontSize',9);
  xran = [0 ted-tst]*sps;
  yran(2) = max(max(stasig));
  yran(1) = min(min(stasig));
  yloc = yran(1)+(yran(2)-yran(1))*14/15;
  ind = find(timemax>=xran(1) & timemax<=xran(2));
  for i = 1: length(ind)
    barst = timecnt(ind(i))-sps*2+1;
    bared = timecnt(ind(i))+sps*2;
    plot(ax,[barst bared],[yloc yloc],'c-','linew',3);
  end
  scatter(ax,timemax(ind), yloc*ones(size(timemax(ind))),20,'y','filled','MarkerEdgeColor','k');
  yloc12 = yran(1)+(yran(2)-yran(1))*1/6;
  yloc13 = yran(1)+(yran(2)-yran(1))*1/12;
  text(ax,timemax(ind),yloc12*ones(size(timemax(ind))),num2str(off12all(ind)),...
    'HorizontalAlignment','right','FontSize',7,'color','b');
%   text(ax,timemax(ind),yloc13*ones(size(timemax(ind))),num2str(off13all(ind)),...
%     'HorizontalAlignment','right','FontSize',7,'color','k');
  ind = find(timearm>=xran(1) & timearm<=xran(2));
  if ~isempty(ind)
    yloc = yran(1)+(yran(2)-yran(1))*13/15;
    scatter(ax,timearm(ind), yloc*ones(size(timearm(ind))),15,'g','filled','MarkerEdgeColor','k');
  end
  ind = find(timebost>=xran(1) & timebost<=xran(2));
  if ~isempty(ind)
    yloc = yran(1)+(yran(2)-yran(1))*12/15;
    scatter(ax,timebost(ind), yloc*ones(size(timebost(ind))),15,'m','filled','MarkerEdgeColor','k');
  end
  
  xlim(ax,xran);
  ylim(ax,yran);
  xlabel(ax,'Samples','FontSize',10);
  ylabel(ax,'Amplitude','FontSize',10);
  
  yyaxis(ax,'right');
  if ista == 1
    plot(ax,ircc,rcc23,'o','color',[.8 .8 .8],'markersize',2);
    plot(ax,ircccon,rcc23c23,'o','color',[.5 .5 .5],'markersize',2);
  elseif ista == 2
    plot(ax,ircc,rcc13,'o','color',[.8 .8 .8],'markersize',2);
    plot(ax,ircccon,rcc13c13,'o','color',[.5 .5 .5],'markersize',2);
  else
    plot(ax,ircc,rcc12,'o','color',[.8 .8 .8],'markersize',2);
    plot(ax,ircccon,rcc12c12,'o','color',[.5 .5 .5],'markersize',2);
  end
  ylabel(ax,'Running CC','FontSize',10);
  ylim(ax,[-1 1]);
  hold(ax,'off');


  isub = ista*2;
  ax = f.ax(isub);
  hold(ax,'on');
  yyaxis(ax,'left');
  if ista == 1
    plot(ax,s2w3,'b-');
    plot(ax,s3w2,'k-'); %,'linew',0.5
  elseif ista == 2
    plot(ax,s1w3,'r-');
    plot(ax,s3w1,'k-');
  else
    plot(ax,s1w2,'r-');
    plot(ax,s2w1,'b-');
  end
  text(ax,0.05,0.95,'Convolved','unit','normalized',...
    'HorizontalAlignment','left','FontSize',9);
  xran = [0 ted-tst]*sps;
  yran(2) = max(max(stasig));
  yran(1) = min(min(stasig));
  xlim(ax,xran);
  ylim(ax,yran);
  xlabel(ax,'Samples','FontSize',10);
  ylabel(ax,'Amplitude','FontSize',10);
  
  yyaxis(ax,'right');
  if ista == 1
    plot(ax,ircc,rcc23,'o','color',[.8 .8 .8],'markersize',2);
    plot(ax,ircccon,rcc23c23,'o','color',[.5 .5 .5],'markersize',2);
  elseif ista == 2
    plot(ax,ircc,rcc13,'o','color',[.8 .8 .8],'markersize',2);
    plot(ax,ircccon,rcc13c13,'o','color',[.5 .5 .5],'markersize',2);
  else
    plot(ax,ircc,rcc12,'o','color',[.8 .8 .8],'markersize',2);
    plot(ax,ircccon,rcc12c12,'o','color',[.5 .5 .5],'markersize',2);
  end
  ylabel(ax,'Running CC','FontSize',10);
  ylim(ax,[-1 1]);

  hold(ax,'off');

end

%%
f.fig = figure;
f.fig.Renderer = 'painters';
widin = 18;  % maximum width allowed is 8.5 inches
htin = 4;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);

ncol = 1;
nrow = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
%     grid(f.ax(isub),'on');
end

pltxran = [0.03 0.95]; pltyran = [0.05 0.95];
pltxsep = 0.02; pltysep = 0.04; 
axpos = optaxposition(nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
for isub = 1:nrow*ncol
  set(f.ax(isub), 'position', axpos(isub,:));
end

%   isub = (ista-1)*2+1;
  ax = f.ax(1);
  hold(ax,'on');
  yyaxis(ax,'left');
  plot(ax,stasig(:,1),'r-');
  plot(ax,stasig(:,2),'b-');
  plot(ax,stasig(:,3),'k-'); %,'linew',0.5
  text(ax,0.05,0.95,'Original','unit','normalized',...
    'HorizontalAlignment','left','FontSize',9);
  xran = [0 ted-tst]*sps;
  yran(2) = max(max(stasig));
  yran(1) = min(min(stasig));
  yloc = yran(1)+(yran(2)-yran(1))*14/15;
  ind = find(timemax>=xran(1) & timemax<=xran(2));
  for i = 1: length(ind)
    barst = timecnt(ind(i))-sps*2+1;
    bared = timecnt(ind(i))+sps*2;
    plot(ax,[barst bared],[yloc yloc],'c-','linew',3);
  end
  scatter(ax,timemax(ind), yloc*ones(size(timemax(ind))),20,'y','filled','MarkerEdgeColor','k');
  yloc12 = yran(1)+(yran(2)-yran(1))*1/6;
  yloc13 = yran(1)+(yran(2)-yran(1))*1/12;
  text(ax,timemax(ind),yloc12*ones(size(timemax(ind))),num2str(off12all(ind)),...
    'HorizontalAlignment','right','FontSize',7,'color','b');
%   text(ax,timemax(ind),yloc13*ones(size(timemax(ind))),num2str(off13all(ind)),...
%     'HorizontalAlignment','right','FontSize',7,'color','k');
  ind = find(timearm>=xran(1) & timearm<=xran(2));
  if ~isempty(ind)
    yloc = yran(1)+(yran(2)-yran(1))*13/15;
    scatter(ax,timearm(ind), yloc*ones(size(timearm(ind))),15,'g','filled','MarkerEdgeColor','k');
  end
  ind = find(timebost>=xran(1) & timebost<=xran(2));
  if ~isempty(ind)
    yloc = yran(1)+(yran(2)-yran(1))*12/15;
    scatter(ax,timebost(ind), yloc*ones(size(timebost(ind))),15,'m','filled','MarkerEdgeColor','k');
  end
  
  xlim(ax,xran);
  ylim(ax,yran);
  xlabel(ax,'Samples','FontSize',10);
  ylabel(ax,'Amplitude','FontSize',10);
  
  yyaxis(ax,'right');
  plot(ax,ircc,rcc,'o','color',[.8 .8 .8],'markersize',2);
  plot(ax,irccc123,rccc123,'o','color',[.5 .5 .5],'markersize',2);
  ylabel(ax,'Running CC','FontSize',10);
  ylim(ax,[-1 1]);
  hold(ax,'off');


%   isub = ista*2;
  ax = f.ax(2);
  hold(ax,'on');
  yyaxis(ax,'left');
  plot(ax,s1w23,'r-');
  plot(ax,s2w13,'b-');
  plot(ax,s3w12,'k-'); %,'linew',0.5
  text(ax,0.05,0.95,'Convolved','unit','normalized',...
    'HorizontalAlignment','left','FontSize',9);
  xran = [0 ted-tst]*sps;
  yran(2) = max(max(stasig));
  yran(1) = min(min(stasig));
  xlim(ax,xran);
  ylim(ax,yran);
  xlabel(ax,'Samples','FontSize',10);
  ylabel(ax,'Amplitude','FontSize',10);
  
  yyaxis(ax,'right');
  plot(ax,ircc,rcc,'o','color',[.8 .8 .8],'markersize',2);
  plot(ax,irccc123,rccc123,'o','color',[.5 .5 .5],'markersize',2);
  ylabel(ax,'Running CC','FontSize',10);
  ylim(ax,[-1 1]);

  hold(ax,'off');

%%
ista = 3 % which station are we focusing on?
sig = stasig(:,ista);
wlet = stawlet(:,ista);

dt = 1/sps;  % sampling interval
if ista == 1
  twlet = 190/sps;  % Time shift of main arrival of wavelet
else
  twlet = 191/sps;  % Time shift of main arrival of wavelet
end
width = 2.5;  % width for Gaussian filter
dres_min = 0.05;  % tolerance, percentage change in residual per iteration
mfit_min = 5e-1;  % tolerance, norm of the residual/misfit
nit_max = 200;  % max numer of iterations 
nit_max = 10000;  % max numer of iterations 
tdura = 0.4;
% nit_max = 2*round((ted-tst)/tdura);  % assume the max distinguishable number allows half overlapping
fpltit = 0;  % plot flag for each iteration
fpltend = 1;  % plot flag for the final iteration
[sigdecon(:,ista),pred,res,dresit,mfitit,ampit,nit,fighdl] = ...
  iterdecon(sig,wlet,rcc,noise,dt,twlet,width,dres_min,mfit_min,nit_max,fpltit,fpltend);

%% Is there a good way to smooth it? Tried envelope, gaussian filter, not ideal
%%% maybe check out the paper on sliding average
close all
for ista = 1: 3
gf = gauss_time(nit_max,5);
sigdecongf(:,ista) = conv(sigdecon(:,ista),gf,'same');
figure
stem(sigdecon(:,ista)); hold on
plot(sigdecongf(:,ista));
[sigdeconenv(:,ista),~] = envelope(sigdecon(:,ista),nit_max);
% plot(sigdeconenv(:,ista));
end

siddeconuse = sigdecongf;
figure
plot(siddeconuse(:,1)/max(siddeconuse(:,1)),'r'); hold on
plot(siddeconuse(:,2)/max(siddeconuse(:,2)),'b');
plot(siddeconuse(:,3)/max(siddeconuse(:,3)),'k');

[coef12,lag12] = xcorr(siddeconuse(:,1), siddeconuse(:,2), sps/2, 'coeff');
[mcoef12, idx] = max(coef12);   % max of master raw cc
lagsamp12 = lag12(idx);   % offset in samples
[coef13,lag13] = xcorr(siddeconuse(:,1), siddeconuse(:,3), sps/2, 'coeff');
[mcoef13, idx] = max(coef13);   % max of master raw cc
lagsamp13 = lag13(idx);   % offset in samples
[coef23,lag23] = xcorr(siddeconuse(:,2), siddeconuse(:,3), sps/2, 'coeff');
[mcoef23, idx] = max(coef23);   % max of master raw cc
lagsamp23 = lag23(idx);   % offset in samples
text(0.1,0.9,sprintf('lag12: %d, max coef12: %.2f', lagsamp12,mcoef12),'Units','normalized',...
  'HorizontalAlignment','left');
text(0.1,0.8,sprintf('lag13: %d, max coef13: %.2f', lagsamp13,mcoef13),'Units','normalized',...
  'HorizontalAlignment','left');
text(0.1,0.7,sprintf('lag23: %d, max coef23: %.2f', lagsamp23,mcoef23),'Units','normalized',...
  'HorizontalAlignment','left');















