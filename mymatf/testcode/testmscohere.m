%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is to test the mocohere and cpsd function
%
% when 2*pi*fmin*dt<pi --> dt<1/2fmin --> nshift<Fs/2fmin, the line of phase would pass through (0,0)
% so it is the case if you use unwrap. So if fmin=0, thoretically any nshift would guarantee passing
% through (0,0), but of course 0 is a singularity itself. The conversion equations to lag samples
% would be the same in this case, i.e. using direct slope, is the same as, using the derivative.
% when it is larger than that, there will be a periodic effect, so it may not pass through (0,0),
% depending on how many cycles; and meanwhile, the converted lag at the start of frequency would not
% equal to the true lag, as it has gone through cycles. In that case, unwrap can be used, but won't
% solve the problem.   
% In any case, the slope of the phase is always accurate to get the real lag by
% slope = nshift/(Fs/2), accordingly, using the derivative is okay to convert to the lag samples as
% well since it has the same logic as the using the slope in normalized lag plot 
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/03/25
% Last modified date:   2021/05/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

% rng('default');
% s = rng;
% rng(s);

%% make a synthetic waveform by stacking various frequencies
sps = 80;
Fs = sps;
df = Fs/512;
% df = 0.5;
freq = 0: df: Fs/2-df;
wlensec = 16;
wlen = wlensec*Fs;
t = 0:1/Fs:wlensec-1/Fs;
NOISE = 'no';
% NOISE = 'yes';
if strcmp(NOISE, 'yes')
    dt = rand(length(freq),1)*3/Fs;     % Gaussian noise
    amp = rand(length(freq),1);
elseif strcmp(NOISE, 'no')
    dt = zeros(length(freq),1);     % no noise
    amp = ones(length(freq),1);
end
for i = 1: length(freq)
    x(:, i) = amp(i)*sin(2*pi*freq(i).*(t + dt(i)));
end
xori = sum(x,2);

% when 2*pi*fmin*dt<pi --> dt<1/2fmin --> nshift<Fs/2fmin, the line of phase would pass through (0,0)
% so it is the case if you use unwrap. So if fmin=0, thoretically any nshift would guarantee passing
% through (0,0), but of course 0 is a singularity itself. The conversion equations to lag samples
% would be the same in this case, i.e. using direct slope, is the same as, using the derivative.
% when it is larger than that, there will be a periodic effect, so it may not pass through (0,0),
% depending on how many cycles; and meanwhile, the converted lag at the start of frequency would not
% equal to the true lag, as it has gone through cycles. In that case, unwrap can be used, but won't
% solve the problem.   
% In any case, the slope of the phase is always accurate to get the real lag by
% slope = nshift/(Fs/2), accordingly, using the derivative is okay to convert to the lag samples as
% well since it has the same logic as the using the slope in normalized lag plot 
% 

nshift = 10;     
dtconst = nshift/Fs;
frange = [5 30];
[~,ist] = min(abs(freq-frange(1)));
[~,ied] = min(abs(freq-frange(2))); 
for i = ist: ied
    x(:, i) = amp(i)*sin(2*pi*freq(i).*(t + dt(i)+dtconst));
end
xsft = sum(x,2);

%% plot two traces
f.fig = figure;
f.fig.Renderer = 'painters';
color=['r','b','k'];
ax = gca;
hold(ax, 'on');
plot(ax,t,xori, 'linewidth', 1,'color',...
    color(1)); hold on
plot(ax,t,xsft, 'linewidth', 1,'color',...
    color(2));
text(ax,0.95,0.9,'Broadband','fontsize',10,'unit','normalized','Horizontalalignment','right');
text(ax,0.95,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
    'right','fontsize',10);
text(ax,0.03,0.92,'a','FontSize',9,'unit','normalized','EdgeColor','k','Margin',2);
% text(ax,0.03,0.1,fam,'fontsize',12,'unit','normalized','Horizontalalignment','left');
% text(ax,0.95,0.1,FLAG,'fontsize',12,'unit','normalized','Horizontalalignment','right');
xlabel(ax,'Time (s)','fontsize',10);
ylabel(ax,'Amplitude','fontsize',10);
xlim(ax,[0 wlensec]);
%     ylim(ax,[0 4]);
ax.Box='on';
grid(ax, 'on');
ax.GridLineStyle = '--';


% %% plot the PSD
% 
% %     PSDfunc = 'pmtm';
% PSDfunc = 'periodogram';
% 
% % pmtm param
% nfft = pow2(nextpow2(size(x,1))-1);
% window = hann(size(x,1));
% nw = 4;
% 
% f.fig = figure;
% f.fig.Renderer = 'painters';
% color=['r','b','k'];
% ax = gca;
% hold(ax, 'on');
% xlim(ax,[1e-2, 1e2]);
% % xlim(ax,[0, Fs/2]);
% ylim(ax,[-100,10]);
% 
% if strcmp(PSDfunc, 'pmtm')
%     [psdx,pft] = pmtm(xori,nw,nfft,Fs);
% elseif strcmp(PSDfunc, 'periodogram')
%     [psdx,pft] = periodogram(xori,window,nfft,Fs);
% end
% p1=plot(ax,pft,pow2db(psdx),'o-','linewidth', 1,'color',color(1),'markers',1.5);
% [~,ind] = max(psdx(10:end));
% plot(ax,[pft(ind+9) pft(ind+9)], ax.YLim, '--','color',color(1));
% 
% if strcmp(PSDfunc, 'pmtm')
%     [psdx,pft] = pmtm(xsft,nw,nfft,Fs);
% elseif strcmp(PSDfunc, 'periodogram')
%     [psdx,pft] = periodogram(xsft,window,nfft,Fs);
% end
% p2=plot(ax,pft,pow2db(psdx),'o-','linewidth', 1,'color',color(2),'markers',1.5);
% [~,ind] = max(psdx(10:end));
% plot(ax,[pft(ind+9) pft(ind+9)], ax.YLim, '--','color',color(2));
%     
% ax.XScale = 'log';
% text(ax,0.05,0.9,PSDfunc,'unit','normalized','fontsize',10);
% text(ax,0.95,0.9,'Broadband','unit','normalized','horizontalalignment','right','fontsize',10);
% text(ax,0.95,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
%     'right','fontsize',10);
% xlabel(ax,'Frequency (Hz)');
% ylabel(ax,'PSD (dB/Hz)');
% lgd = legend(ax,[p1, p2],{'Original','Shifted'},'location','southwest','fontsize',8);
% grid on
% box on
% hold(ax, 'off');
% 
% 
%%% obtain magnitude-squared coherence and cross-PSD of templates
% coherence, 0 to 1
%     nfft = pow2(nextpow2(size(x,1))-1);
nfft = 512;
window = hann(nfft);
noverlap = nfft*2/4;
Fs = sps;
% filter param
lolf = 0.5;
lohf = 1.25;
hihf = 6.5;
% 
% f.fig = figure;
% f.fig.Renderer = 'painters';
% color=['r','b','k'];
% ax = gca;
% hold(ax, 'on');
% %     xlim(ax,[1e-2, 1e2]);
% xlim(ax,[0, Fs/2]);
% ylim(ax,[0,1.2]);
% [Cxy,ft] = mscohere(xori,xsft,window,noverlap,nfft,Fs);
% p(1)=plot(ax,ft,Cxy,'o-','linewidth', 1,'color',color(1),'markers',1.5);
% %         [~,ind] = max(psdx(10:end));
% %         plot(ax,[ft(ind+9) ft(ind+9)], ax.YLim, '--','color',color(ista));
% plot(ax,[lolf lolf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
% plot(ax,[lohf lohf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
% plot(ax,[hihf hihf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
% %     ax.XScale = 'log';
% text(ax,0.05,0.9,'Broadband','unit','normalized','horizontalalignment','left','fontsize',10);
% text(ax,0.05,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
%     'left','fontsize',10);
% text(ax,0.05,0.8,strcat({'nfft: '},num2str(nfft)),'unit','normalized','horizontalalignment',...
%     'left','fontsize',10);
% text(ax,0.05,0.75,strcat({'noverlap: '},num2str(noverlap)),'unit','normalized','horizontalalignment',...
%     'left','fontsize',10);
% text(ax,0.05,0.7,strcat({'Hann'}),'unit','normalized','horizontalalignment',...
%     'left','fontsize',10);
% xlabel(ax,'Frequency (Hz)');
% ylabel(ax,'Magnitude-squared coherence');
% %     lgd = legend(ax,p,{'TWKB-LZB','TWKB-MGCB','LZB-MGCB'},'location','southwest','fontsize',8);
% grid on
% box on
% hold(ax, 'off');


%%% phase of cross spectrum in normalized phase lag
%%% the slope of the normalized phase lag = nshift / F_Nyquist, it is most reliable, even though
%%% sometimes the best-fit line does not pass through 0 hz and 0 lag
f.fig = figure;
f.fig.Renderer = 'painters';
color=['r','b','k'];
ax = gca;
hold(ax, 'on');
%     xlim(ax,[1e-2, 1e2]);
xlim(ax,[0, Fs/2]);
%     ylim(ax,[-1.2,1.2]);
[Pxy,ft] = cpsd(xori,xsft,window,noverlap,nfft,Fs);
% p(1)=plot(ax,ft,unwrap(angle(Pxy))/pi,'o-','linewidth', 1,'color',...
%           color(1),'markers',1.5);
p(1)=plot(ax,ft,angle(Pxy)/pi,'o-','linewidth', 1,'color',...
          color(1),'markers',1.5);
[~,ind] = min(abs(ft-frange(1)));
disp(angle(Pxy(ind+1))/pi/2*Fs/ft(ind+1));      
plot(ax,[lolf lolf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
plot(ax,[lohf lohf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
plot(ax,[hihf hihf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
%     ax.XScale = 'log';
text(ax,0.05,0.9,'Broadband','unit','normalized','horizontalalignment','left','fontsize',10);
text(ax,0.05,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
    'left','fontsize',10);
text(ax,0.05,0.8,strcat({'nfft: '},num2str(nfft)),'unit','normalized','horizontalalignment',...
    'left','fontsize',10);
text(ax,0.05,0.75,strcat({'noverlap: '},num2str(noverlap)),'unit','normalized','horizontalalignment',...
    'left','fontsize',10);
text(ax,0.05,0.7,strcat({'Hann'}),'unit','normalized','horizontalalignment',...
    'left','fontsize',10);
xlabel(ax,'Frequency (Hz)');
ylabel(ax,'Normalized Phase lag');
%     lgd = legend(ax,p,{'TWKB-LZB','TWKB-MGCB','LZB-MGCB'},'location','southwest','fontsize',8);
grid on
box on
hold(ax, 'off');

%%% phase of cross spectrum in lag samples
%%% the lag samples would converge to the correct value, but it is not reliable all the time, for
%%% example, if nshift = 15 samples, Fs=80, it would fail
f.fig = figure;
f.fig.Renderer = 'painters';
color=['r','b','k'];
ax = gca;
hold(ax, 'on');
%     xlim(ax,[1e-2, 1e2]);
xlim(ax,[0, Fs/2]);
%     ylim(ax,[-1.5*nshift,1.5*nshift]);
[Pxy,ft] = cpsd(xori,xsft,window,noverlap,nfft,Fs);
% p(1)=plot(ax,ft,unwrap(angle(Pxy))/2/pi*Fs./ft,'o-','linewidth', 1,'color',...
%     color(1),'markers',1.5);
lag = diff(unwrap(angle(Pxy))/pi)/2*Fs./(ft(2)-ft(1));
p(1)=plot(ax,ft,[0; lag],'o-','linewidth', 1,'color',...
    color(1),'markers',1.5);
%     p(ista)=plot(ax,ft,angle(Pxy)/2/pi*Fs./ft,'o-','linewidth', 1,'color',...
%                      color(ista),'markers',1.5);
plot(ax,[lolf lolf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
plot(ax,[lohf lohf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
plot(ax,[hihf hihf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
%     ax.XScale = 'log';
text(ax,0.05,0.9,'Broadband','unit','normalized','horizontalalignment','left','fontsize',10);
text(ax,0.05,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
    'left','fontsize',10);
text(ax,0.05,0.8,strcat({'nfft: '},num2str(nfft)),'unit','normalized','horizontalalignment',...
    'left','fontsize',10);
text(ax,0.05,0.75,strcat({'noverlap: '},num2str(noverlap)),'unit','normalized','horizontalalignment',...
    'left','fontsize',10);
text(ax,0.05,0.7,strcat({'Hann'}),'unit','normalized','horizontalalignment',...
    'left','fontsize',10);
xlabel(ax,'Frequency (Hz)');
ylabel(ax,'Lag (sample)');
%     lgd = legend(ax,p,{'TWKB-LZB','TWKB-MGCB','LZB-MGCB'},'location','southwest','fontsize',8);
grid on
box on
hold(ax, 'off');











