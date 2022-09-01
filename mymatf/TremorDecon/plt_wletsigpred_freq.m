function [f] = plt_wletsigpred_freq(wlet,sig,pred,sps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_wletsigpred_freq(wletbb,sigbb,wlet,sig,sps,bpwlet,bpsig)
%
% This function is to plot the spectrum of template, signal and the prediction
% from modeling. The expectation could be, if the processed template and signal
% have a similar spectral shape so that you can deconvolve the template from 
% the signal. The predicted waveform from convolving the deconvolved impulses
% and signal should also have a similar spectral shape to those. 
% 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/01/20
% Last modified date:   2022/01/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

f.fig = figure;
widin = 5;  % maximum width allowed is 8.5 inches
htin = 6;   % maximum height allowed is 11 inches
nrow = 4;
% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, resol] = pixelperinch(1);
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);

f.ax(1)=subplot(nrow,1,1);
hold on
box on
grid on
lwlet = length(wlet);
plot(1:lwlet,wlet,'r','linew',1);
xlim([0 lwlet]);
ylim([-max(abs(wlet)) max(abs(wlet))]);
legend('Wavelet','location','southeast');
%   text(0.9,0.9,stas(ista,:),'unit','normalized');
xlabel('Samples');
ylabel('Amplitude');

f.ax(2)=subplot(nrow,1,2);
hold on
box on
grid on
lsig = length(sig);
plot(1:lsig,sig,'k','linew',1);
plot(1:lsig,pred,'b','linew',1);
xlim([0 lsig]);
ylim([-max(abs(sig)) max(abs(sig))]);
legend('Signal','Prediction','location','southeast');
%   text(0.9,0.9,stas(ista,:),'unit','normalized');
xlabel('Samples');
ylabel('Amplitude');

f.ax(3)=subplot(nrow,1,[3 4]);
%%% plot the spectrum, using 'pchave', which is capable of time-averaging a longer segment
wlen = lwlet;
polap = 75;
nfft = wlen;
Fs = sps;
%note that i request for only 1 window, so there is no robust estimate for the template, it is
%pretty much the same as 'periodogram', CHECK IT
[~,pcfw,~,~,~,~,pcwlet]=pchave(wlet,wlen,polap,nfft,Fs,'MAD','dpss');
%for the signal and noise, i want a robust estimate from an average of overlapping windows
if lsig<lwlet
  [pcsig,pcfs]=pchave(sig,pow2(nextpow2(lsig)-1),polap,pow2(nextpow2(lsig)-1),Fs,'MAD','dpss');
  [pcpred,~]=pchave(pred,pow2(nextpow2(lsig)-1),polap,pow2(nextpow2(lsig)-1),Fs,'MAD','dpss');
else
  [pcsig,pcfs]=pchave(sig,wlen,polap,nfft,Fs,'MAD','dpss');
  [pcpred,~]=pchave(pred,wlen,polap,nfft,Fs,'MAD','dpss');
end
ax = gca;
loglog(ax,pcfw,sqrt(pcwlet),'r-','linew',1.5);
hold(ax,'on');
% grid(ax,'on');
loglog(ax,pcfs,sqrt(pcsig),'k-','linew',1.5);
loglog(ax,pcfs,sqrt(pcpred),'b-','linew',1.5);

yran = [1e-5 1e1];
plot(ax,[1 1],yran,':','color',[0.7 0.7 0.7]);
plot(ax,[2 2],yran,':','color',[0.7 0.7 0.7]);
plot(ax,[4 4],yran,':','color',[0.7 0.7 0.7]);
plot(ax,[8 8],yran,':','color',[0.7 0.7 0.7]);
xticks(ax,[0.1 1 10]);
xlim(ax,[0.1 10]);
ylim(ax,yran);
%     ax.GridLineStyle = '-';
lgd = legend(ax,'Wavelet','Signal','Prediction','location','northwest','fontsize',7);
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));
xlabel(ax,'Frequency');
ylabel(ax,'Amplitude');
hold(ax,'off');
















