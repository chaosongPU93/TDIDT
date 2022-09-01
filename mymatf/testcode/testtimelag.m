% test the relationship between the time lag from min to max of a stack
% of sinuoids of different frequencies plus noise and the amplitude of
% each contributing sinuoids.

%%
clear
clc
close all

Fs = 40;
df = 0.5;
freq = 0.5: df: 6.5;
wlensec = 2*1/min(freq);
wlen = wlensec*Fs;
t = -wlensec/2:1/Fs:wlensec/2-1/Fs;
% dt = randn(length(freq),1)*2/Fs;     % Gaussian noise
dt = zeros(length(freq),1);     % no noise
% NOISE = 'no';
NOISE = 'yes';


ampflat = ones(length(freq),1);
for i = 1: length(freq)
    xflat(:, i) = ampflat(i)*sin(2*pi*freq(i).*(t + dt(i)));
end
noiscale = 0.5;
if strcmp(NOISE, 'yes')
    noigau = randn(wlen,1)*noiscale;
elseif strcmp(NOISE, 'no')
    noigau = zeros(wlen,1);
end
xsflat = sum(xflat,2)+noigau;
seglt = xsflat(wlen/4+1:wlen/2);
[~,indmin] = min(seglt);
indmin = indmin+wlen/4;  % convert to global index
segrt = xsflat(wlen/2+1:wlen/4*3);
[~,indmax] = max(segrt);
indmax = indmax+wlen/2;  % convert to global index
tlagflat = (indmax-indmin)/Fs;


% ampinc = logspace(-1,2,length(freq));
ampinc = linspace(0.1,100,length(freq));
for i = 1: length(freq)
    xinc(:, i) = ampinc(i)*sin(2*pi*freq(i).*(t + dt(i)));
end
noiscale = 5;
if strcmp(NOISE, 'yes')
    noigau = randn(wlen,1)*noiscale;
elseif strcmp(NOISE, 'no')
    noigau = zeros(wlen,1);
end
xsinc = sum(xinc,2)+noigau;
seglt = xsinc(wlen/4+1:wlen/2);
[~,indmin] = min(seglt);
indmin = indmin+wlen/4;  % convert to global index
segrt = xsinc(wlen/2+1:wlen/4*3);
[~,indmax] = max(segrt);
indmax = indmax+wlen/2;  % convert to global index
tlaginc = (indmax-indmin)/Fs;


% ampdec = logspace(2,-1,length(freq));
ampdec = linspace(100,0.1,length(freq));
for i = 1: length(freq)
    xdec(:, i) = ampdec(i)*sin(2*pi*freq(i).*(t + dt(i)));
end
noiscale = 5;
if strcmp(NOISE, 'yes')
    noigau = randn(wlen,1)*noiscale;
elseif strcmp(NOISE, 'no')
    noigau = zeros(wlen,1);
end
xsdec = sum(xdec,2)+noigau;
seglt = xsdec(wlen/4+1:wlen/2);
[~,indmin] = min(seglt);
indmin = indmin+wlen/4;  % convert to global index
segrt = xsdec(wlen/2+1:wlen/4*3);
[~,indmax] = max(segrt);
indmax = indmax+wlen/2;  % convert to global index
tlagdec = (indmax-indmin)/Fs;


%% plot two traces
f.fig = figure;
f.fig.Renderer = 'painters';
color=['k','b','r'];
subplot(311)
ax = gca;
hold(ax, 'on');
for i = 1: length(freq)
    plot(ax,t,xflat(:,i));
end
plot(ax,t,xsflat, 'linewidth', 1,'color',color(1));
% plot(ax,t,xsft, 'linewidth', 1,'color',...
%     color(2));
% text(ax,0.95,0.9,'Broadband','fontsize',10,'unit','normalized','Horizontalalignment','right');
text(ax,0.6,0.85,sprintf('%.2f s',tlagflat),'fontsize',8,'unit',...
     'normalized','Horizontalalignment','left');
% text(ax,0.03,0.92,'a','FontSize',9,'unit','normalized','EdgeColor','k','Margin',2);
% text(ax,0.03,0.1,fam,'fontsize',12,'unit','normalized','Horizontalalignment','left');
% text(ax,0.95,0.1,FLAG,'fontsize',12,'unit','normalized','Horizontalalignment','right');
xlabel(ax,'Time (s)','fontsize',10);
ylabel(ax,'Amplitude','fontsize',10);
% xlim(ax,[0 wlensec]);
%     ylim(ax,[0 4]);
ax.Box='on';
grid(ax, 'on');
ax.GridLineStyle = '--';

subplot(312)
ax = gca;
hold(ax, 'on');
for i = 1: length(freq)
    plot(ax,t,xinc(:,i));
end
plot(ax,t,xsinc, 'linewidth', 1,'color',color(1));
text(ax,0.6,0.85,sprintf('%.2f s',tlaginc),'fontsize',8,'unit',...
     'normalized','Horizontalalignment','left');
xlabel(ax,'Time (s)','fontsize',10);
ylabel(ax,'Amplitude','fontsize',10);
ax.Box='on';
grid(ax, 'on');
ax.GridLineStyle = '--';

subplot(313)
ax = gca;
hold(ax, 'on');
for i = 1: length(freq)
    plot(ax,t,xdec(:,i));
end
plot(ax,t,xsdec, 'linewidth', 1,'color',color(1));
text(ax,0.6,0.85,sprintf('%.2f s',tlagdec),'fontsize',8,'unit',...
     'normalized','Horizontalalignment','left');
xlabel(ax,'Time (s)','fontsize',10);
ylabel(ax,'Amplitude','fontsize',10);
ax.Box='on';
grid(ax, 'on');
ax.GridLineStyle = '--';











