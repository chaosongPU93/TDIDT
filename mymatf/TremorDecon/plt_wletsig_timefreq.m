function [f] = plt_wletsig_timefreq(wletbb,sigbb,wlet,sig,noibb,noi,sps,bpwlet,bpsig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_wletsig_timefreq(wletbb,sigbb,wlet,sig,sps,bpwlet,bpsig)
%
% This function is to plot the waveform and spectrum of template and signal
% that will be fed to deconvolution routines before and after bandpass 
% filtering. From the plot we need to make sure that the signal and template
% have a similar spectral shape after filtering, which has to be the case if 
% you want to carry out deconvolution between the two.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/01/13
% Last modified date:   2022/01/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

f.fig = figure;
f.fig.Renderer = 'painters';
widin = 5;  % maximum width allowed is 8.5 inches
if isempty(noi)
  htin = 8.5;   % maximum height allowed is 11 inches
  nrow = 6;
else
  htin = 10.5;
  nrow = 8;
end
% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, resol] = pixelperinch(1);
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);

lwlet = length(wletbb);
lsig = length(sigbb);
if ~isempty(noi)
  lnoi = length(noibb);
else
  lnoi = [];
end
lmax = max([lwlet; lsig; lnoi]);

f.ax(1)=subplot(nrow,1,1);
hold on
box on
grid on
plot(1:lwlet,wletbb,'r','linew',1);
xlim([0 lwlet]);
ylim([-max(abs(wletbb)) max(abs(wletbb))]);
legend('BB temp.','location','southeast');
%   text(0.9,0.9,stas(ista,:),'unit','normalized');
% xlabel('Samples');
% ylabel('Amplitude','FontSize',10);
nolabels(gca,1); longticks(gca,2); axsym(gca,2); axranexp(gca,6);

f.ax(2)=subplot(nrow,1,2);
hold on
box on
grid on
plot(1:lwlet,wlet,'r','linew',1);
xlim([0 lwlet]);
ylim([-max(abs(wlet)) max(abs(wlet))]);
legend(sprintf('BP temp., %.1f-%.1f Hz',bpwlet(1),bpwlet(2)),'location','southeast');
%   text(0.9,0.9,stas(ista,:),'unit','normalized');
% xlabel('Samples');
ylabel('Amplitude','FontSize',10);
longticks(gca,2); axsym(gca,2); axranexp(gca,6);

f.ax(3)=subplot(nrow,1,3);
hold on
box on
grid on
plot(1:lsig,sigbb,'k','linew',1);
xlim([0 lsig]);
ylim([-max(abs(sigbb)) max(abs(sigbb))]);
legend(sprintf('BB sig.'),'location','southeast');
%   text(0.9,0.9,stas(ista,:),'unit','normalized');
% xlabel('Samples');
% ylabel('Amplitude','FontSize',10);
nolabels(gca,1); longticks(gca,2); axsym(gca,2); axranexp(gca,6);

f.ax(4)=subplot(nrow,1,4);
hold on
box on
grid on
plot(1:lsig,sig,'k','linew',1);
xlim([0 lsig]);
ylim([-max(abs(sig)) max(abs(sig))]);
legend(sprintf('BP sig., %.1f-%.1f Hz',bpsig(1),bpsig(2)),'location','southeast');
%   text(0.9,0.9,stas(ista,:),'unit','normalized');
% xlabel('Samples');
ylabel('Amplitude','FontSize',10);
longticks(gca,2); axsym(gca,2); axranexp(gca,6);

if isempty(noi)
  f.ax(5)=subplot(nrow,1,[5 6]);
else
  f.ax(7)=subplot(nrow,1,[7 8]);
end
%%% plot the spectrum, using 'pchave', which is capable of time-averaging a longer segment
wlen = lwlet;
polap = 75;
nfft = wlen;
Fs = sps;
%note that i request for only 1 window, so there is no robust estimate for the template, it is
%pretty much the same as 'periodogram', CHECK IT
[~,pcfw,~,~,~,~,pcwletbb]=pchave(wletbb,wlen,polap,nfft,Fs,'MAD','dpss');
[~,~,~,~,~,~,pcwlet]=pchave(wlet,wlen,polap,nfft,Fs,'MAD','dpss');
%for the signal and noise, i want a robust estimate from an average of overlapping windows
if lsig<lwlet
  [pcsigbb,pcfs]=pchave(sigbb,pow2(nextpow2(lsig)-1),polap,pow2(nextpow2(lsig)-1),Fs,'MAD','dpss');
  [pcsig,~]=pchave(sig,pow2(nextpow2(lsig)-1),polap,pow2(nextpow2(lsig)-1),Fs,'MAD','dpss');
else
  [pcsigbb,pcfs]=pchave(sigbb,wlen,polap,nfft,Fs,'MAD','dpss');
  [pcsig,~]=pchave(sig,wlen,polap,nfft,Fs,'MAD','dpss');
end
if ~isempty(noi)
  if lnoi<lwlet
    [pcnoibb,pcfn]=pchave(noibb,pow2(nextpow2(lnoi)-1),75,pow2(nextpow2(lnoi)-1),Fs,'MAD','dpss');
    [pcnoi,~]=pchave(noi,pow2(nextpow2(lnoi)-1),75,pow2(nextpow2(lnoi)-1),Fs,'MAD','dpss');
  else
    [pcnoibb,pcfn]=pchave(noibb,wlen,polap,nfft,Fs,'MAD','dpss');
    [pcnoi,~]=pchave(noi,wlen,polap,nfft,Fs,'MAD','dpss');
  end
end
ax = gca;
loglog(ax,pcfw,sqrt(pcwletbb),'r--','linew',1);
hold(ax,'on');
% grid(ax,'on');
loglog(ax,pcfs,sqrt(pcsigbb),'k--','linew',1);
if ~isempty(noi)
  loglog(ax,pcfn,sqrt(pcnoibb),'--','color',[.7 .7 .7],'linew',1);
  loglog(ax,pcfn,sqrt(pcnoi),'color',[.7 .7 .7],'linew',1.5);
end

loglog(ax,pcfw,sqrt(pcwlet),'color','r','linew',1.5);
loglog(ax,pcfs,sqrt(pcsig),'color','k','linew',1.5);

amprat = mean(sqrt(pcsig(pcfs>=bpsig(1) & pcfs<=bpsig(2))))/ ...
         mean(sqrt(pcwlet(pcfw>=bpsig(1) & pcfw<=bpsig(2))));
nsqrtpcsig = sqrt(pcsig)/amprat;
loglog(ax,pcfs,nsqrtpcsig,'color',[.4 .4 .4],'linew',1.5);


yran = [1e-5 1e1];
plot(ax,[1 1],yran,':','color',[0.7 0.7 0.7]);
plot(ax,[2 2],yran,':','color',[0.7 0.7 0.7]);
plot(ax,[4 4],yran,':','color',[0.7 0.7 0.7]);
plot(ax,[8 8],yran,':','color',[0.7 0.7 0.7]);
plot(ax,[bpsig(1) bpsig(1)],yran,'-.','color',[0.5 0.5 0.5]);
plot(ax,[bpsig(2) bpsig(2)],yran,'-.','color',[0.5 0.5 0.5]);
plot(ax,[bpwlet(1) bpwlet(1)],yran,'-.','color',[1 0.65 0]);
plot(ax,[bpwlet(2) bpwlet(2)],yran,'-.','color',[1 0.65 0]);
xticks(ax,[0.1 1 10]);
xlim(ax,[0.1 10]);
ylim(ax,yran);
%     ax.GridLineStyle = '-';
if ~isempty(noi)
  lgd = legend(ax,'BB temp.',sprintf('BB sig.'),...
    sprintf('BB noi.'),sprintf('BP noi., %.1f-%.1f Hz',bpsig(1),bpsig(2)),...
    sprintf('BP temp., %.1f-%.1f Hz',bpwlet(1),bpwlet(2)),...
    sprintf('BP sig., %.1f-%.1f Hz',bpsig(1),bpsig(2)),...
    'Norm. sig.',...
    'location','northwest','fontsize',7);
else
  lgd = legend(ax,'BB temp.',sprintf('BB sig.'),...
    sprintf('BP temp., %.1f-%.1f Hz',bpwlet(1),bpwlet(2)),...
    sprintf('BP sig., %.1f-%.1f Hz',bpsig(1),bpsig(2)),...
    'Norm. sig.','location','northwest','fontsize',7);% make background transparent
end
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));
xlabel(ax,'Frequency (Hz)','FontSize',10);
ylabel(ax,'Amplitude','FontSize',10);
longticks(ax,2);
hold(ax,'off');
  
if ~isempty(noi)
  f.ax(5)=subplot(nrow,1,5);
  hold on
  box on
  grid on
  plot(1:lnoi,noibb,'color',[.7 .7 .7],'linew',1);
  xlim([0 lnoi]);
  ylim([-max(abs(sigbb)) max(abs(sigbb))]);
  legend(sprintf('BB noi.'),'location','southeast');
  %   text(0.9,0.9,stas(ista,:),'unit','normalized');
%   xlabel('Samples');
%   ylabel('Amplitude','FontSize',10);
  nolabels(gca,1); longticks(gca,2); axsym(gca,2); axranexp(gca,6);

  f.ax(6)=subplot(nrow,1,6);
  hold on
  box on
  grid on
  plot(1:lnoi,noi,'color',[.7 .7 .7],'linew',1);
  xlim([0 lnoi]);
  ylim([-max(abs(sig)) max(abs(sig))]);
  legend(sprintf('BP noi., %.1f-%.1f Hz',bpsig(1),bpsig(2)),'location','southeast');
  %   text(0.9,0.9,stas(ista,:),'unit','normalized');
  xlabel(f.ax(6),sprintf('Samples at %d Hz',sps),'FontSize',10);
  ylabel('Amplitude','FontSize',10);
  longticks(gca,2); axsym(gca,2); axranexp(gca,6);
else
  xlabel(f.ax(4),sprintf('Samples at %d Hz',sps),'FontSize',10);
end      

  
  