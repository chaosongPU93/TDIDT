% plt_moment_duration_axes.m


clear;
clc;
close all;

htin = 6;   % maximum height allowed is 11 inches
widin = 6;  % maximum width allowed is 8.5 inches
nrow = 1;
ncol = 1;
f = initfig(widin,htin,nrow,ncol); %initialize fig

xran1 = [10 21];
yran1 = [-2 9];

t = tiledlayout(f.fig,1,1);
ax1 = axes(t);
hold(ax1,'on'); ax1.Box = 'on'; grid(ax1,'off'); axis(ax1,'equal');
xlim(ax1,xran1);
xticks(ax1,xran1(1):1:xran1(2));
ylim(ax1,yran1);
xlabel(ax1,'log_{10}{Seismic moment (Nm)}','FontSize',12);
ylabel(ax1,'log_{10}{Characteristic duration (s)}','FontSize',12);
plot(ax1,xran1,[log10(3600) log10(3600)],'k--','linew',0.5);
text(ax1,10.5,log10(3600)+0.3,'1 hour','FontSize',12);
plot(ax1,xran1,[log10(24*3600) log10(24*3600)],'k--','linew',0.5);
text(ax1,10.5,log10(24*3600)+0.3,'1 day','FontSize',12);
plot(ax1,xran1,[log10(24*3600*30) log10(24*3600*30)],'k--','linew',0.5);
text(ax1,10.5,log10(24*3600*30)+0.3,'1 month','FontSize',12);
plot(ax1,xran1,[log10(24*3600*30*12) log10(24*3600*30*12)],'k--','linew',0.5);
text(ax1,10.5,log10(24*3600*30*12)+0.3,'1 year','FontSize',12);
longticks(ax1,2);

ax2 = axes(t);
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';

[~,~,xran2] = moment2mag(10.^(xran1));
xlim(ax2,xran2);
xticks(ax2,1:1:8);
ylim(ax2,yran1);
% yticks(ax2,yran1);
yticklabels(ax2,[]);
% noticks(ax2,2);
xlabel(ax2,'M_{w}','FontSize',12);

longticks(ax2,2);



fname = strcat('m-t_axes.pdf');
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));



