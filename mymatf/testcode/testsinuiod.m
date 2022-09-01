% test the composition and decomposition of sine wave

%% initialization
clc;
clear;
close all;

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);

rstpath = '/home/chaosong/Desktop/meetinglfhf';

%%
dx = 0.001*pi;
x = -1*pi:dx:1*pi;
nx = length(x);

dfreq = 0.05;
freq = 0.5: dfreq: 6.5;     % frequency
nfreq = length(freq);

a = 1;  % amplitude

x0 = 1/4*pi;

y = zeros(nx, nfreq);
for i = 1: nfreq
    y(:, i) = a.*sin(2*pi*freq(i).*(x-x0));
end

ysum = sum(y,2);

% hf bandpass
lohf = 1.25; 
hihf = 6.5;
npo = 2;
npa = 2;
fs = 1/dx;
ysumbphf = Bandpass(ysum, fs, lohf, hihf, npo, npa, 'butter');
ind = find(freq==lohf);
ysumhf = sum(y(:,ind:end), 2);

% lf bandpass
lolf = 0.5; 
hilf = 1.25;
npo = 2;
npa = 2;
fs = 1/dx;
ysumbplf = Bandpass(ysum, fs, lolf, hilf, npo, npa, 'butter');
ysumlf = sum(y(:,1:ind), 2);


%% plot
%%% plot 1
f1.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 6;   % maximum height allowed is 11 inches
set(f1.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

ax = gca;
hold(ax,'on');
ax.Box = 'on';
grid(ax,'on');
plot(ax,[x0 x0],[-1000 1000],'--','color',[.5 .5 .5],'linew',1);
for i = 1: 10: nfreq
    plot(ax, x, y(:,i));    %,'-','linew',0.5
end

p1 = plot(ax, x, y(:,end),'r-','linew',1);
p2 = plot(ax, x, y(:,1),'b-','linew',1);
p3 = plot(ax, x, ysum,'k-','linew',1.5);

legend(ax,[p1, p2, p3],{'Max freq: 6.5 Hz', 'Min freq: 0.5 Hz', 'Stack of all freq'},...
       'location','southwest');
xlim(ax, [min(x) max(x)]);
ylim(ax, [-nfreq nfreq]);
ax.XTick = min(x):0.5*pi:max(x);
% ax.XTickLabel = ['-2pi';'-1.5pi';'-1pi';'-0.5pi';'0';'0.5pi';'1pi';'1.5pi';'2pi'];

text(ax,0.15,0.85,strcat({'num. of freq: '}, num2str(nfreq)),'FontSize',12,'unit','normalized');
text(ax,0.15,0.75,strcat({'interval of freq: '}, num2str(dfreq),{' Hz'}),'FontSize',12,...
     'unit','normalized');
title(ax,'In-phase stack of sine waves with different freq.')
hold(ax,'off');

print(f1.fig,'-dpdf',strcat(rstpath,'/inphase_sine_stack.pdf'));


%%% plot 2
f2.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 11;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 3;
ncol = 1;

for isub = 1:nrow*ncol
    f2.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
% set(f2.ax(1), 'position', [ 0.1, 0.1, 0.35, 0.75]);
% set(f2.ax(2), 'position', [ 0.55, 0.1, 0.35, 0.75]);


ax = f2.ax(1);
hold(ax,'on');
ax.Box = 'on';
grid(ax,'on');
plot(ax,[x0 x0],[-1000 1000],'--','color',[.5 .5 .5],'linew',1);
for i = 1: 10: nfreq
    plot(ax, x, y(:,i));    %,'-','linew',0.5
end
p1 = plot(ax, x, y(:,end),'r-','linew',1);
p2 = plot(ax, x, y(:,1),'b-','linew',1);
p3 = plot(ax, x, ysum,'k-','linew',1.5);
legend(ax,[p1, p2, p3],{'Max freq: 6.5 Hz', 'Min freq: 0.5 Hz', 'Stack of all freq'},...
       'location','southwest');
xlim(ax, [min(x) max(x)]);
ylim(ax, [-nfreq nfreq]);
ax.XTick = min(x):0.5*pi:max(x);
% ax.XTickLabel = ['-2pi';'-1.5pi';'-1pi';'-0.5pi';'0';'0.5pi';'1pi';'1.5pi';'2pi'];

text(ax,0.15,0.85,strcat({'num. of freq: '}, num2str(nfreq)),'FontSize',12,'unit','normalized');
text(ax,0.15,0.75,strcat({'interval of freq: '}, num2str(dfreq),{' Hz'}),'FontSize',12,...
     'unit','normalized');
title(ax,'In-phase stack of sine waves with different freq.')
hold(ax,'off');


ax = f2.ax(2);
hold(ax,'on');
ax.Box = 'on';
grid(ax,'on');
plot(ax,[x0 x0],[-1000 1000],'--','color',[.5 .5 .5],'linew',1);
p1 = plot(ax, x, ysumhf,'r-','linew',1.5);
p2 = plot(ax, x, ysumbphf,'b-','linew',1.5);
legend(ax,[p1, p2],{'Stack of freq in 1.25-6.5 Hz', 'Bandpassed stack using 1.25-6.5 Hz'},...
       'location','southwest');
xlim(ax, [min(x) max(x)]);
ylim(ax, [-nfreq nfreq]);
ax.XTick = min(x):0.5*pi:max(x);
% ax.XTickLabel = ['-2pi';'-1.5pi';'-1pi';'-0.5pi';'0';'0.5pi';'1pi';'1.5pi';'2pi'];

% text(ax,0.15,0.85,strcat({'num. of freq: '}, num2str(nfreq)),'FontSize',12,'unit','normalized');
% text(ax,0.15,0.75,strcat({'interval of freq: '}, num2str(dfreq),{' Hz'}),'FontSize',12,...
%      'unit','normalized');
title(ax,'Comparison between direct stack of certain freq range and bandpassed stack with sam freq')
hold(ax,'off');


ax = f2.ax(3);
hold(ax,'on');
ax.Box = 'on';
grid(ax,'on');
plot(ax,[x0 x0],[-1000 1000],'--','color',[.5 .5 .5],'linew',1);
p1 = plot(ax, x, ysumlf,'r-','linew',1.5);
p2 = plot(ax, x, ysumbplf,'b-','linew',1.5);
legend(ax,[p1, p2],{'Stack of freq in 0.5-1.25 Hz', 'Bandpassed stack using 0.5-1.25 Hz'},...
       'location','southwest');
xlim(ax, [min(x) max(x)]);
ylim(ax, [-nfreq nfreq]);
ax.XTick = min(x):0.5*pi:max(x);
% ax.XTickLabel = ['-2pi';'-1.5pi';'-1pi';'-0.5pi';'0';'0.5pi';'1pi';'1.5pi';'2pi'];

% text(ax,0.15,0.85,strcat({'num. of freq: '}, num2str(nfreq)),'FontSize',12,'unit','normalized');
% text(ax,0.15,0.75,strcat({'interval of freq: '}, num2str(dfreq),{' Hz'}),'FontSize',12,...
%      'unit','normalized');
title(ax,'Comparison between direct stack of certain freq range and bandpassed stack with sam freq')
hold(ax,'off');

print(f2.fig,'-dpdf',strcat(rstpath,'/bandpass_comparison.pdf'));
