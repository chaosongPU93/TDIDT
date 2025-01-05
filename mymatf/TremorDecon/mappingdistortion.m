%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is to simply convert the a few testing locations in offset 
% space to map space, in order to gauge the mapping distortion caused by the
% the velocity model and station distribution.
% The testing locations are designed to be align off12==off13, and 
% off12+off13==0. 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/11/20
% Last modified date:   2024/11/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e   % Set the format to 5-digit floating point
clear
clc
close all

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, resol] = pixelperinch(1);

fam = '002';   % family number

ftrans = 'interpchao';

sps = 160;

%% testing locations in offset
off12 = (-20: 2: 20)';
off13 = (-20: 2: 20)';
line1 = [off12 off13];
line2 = [off12 -off13];
line3 = [zeros(length(off13),1) -off13];
line4 = [off12 zeros(length(off12),1)];

locline1 = off2space002(line1,sps,ftrans,0);
locline2 = off2space002(line2,sps,ftrans,0);
locline3 = off2space002(line3,sps,ftrans,0);
locline4 = off2space002(line4,sps,ftrans,0);

%% plot
f=initfig(8,4,1,2);

ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on'); axis(ax, 'equal'); 
plot(ax,line1(:,1)/sps,line1(:,2)/sps,'r-','marker','o','markersize',4,'linew',1);
plot(ax,line2(:,1)/sps,line2(:,2)/sps,'b-','marker','o','markersize',4,'linew',1);
plot(ax,line3(:,1)/sps,line3(:,2)/sps,'k-','marker','o','markersize',4,'linew',1);
plot(ax,line4(:,1)/sps,line4(:,2)/sps,'g-','marker','o','markersize',4,'linew',1);
scatter(ax,line1(1,1)/sps,line1(1,2)/sps,50,'r^','filled');  %start
scatter(ax,line1(end,1)/sps,line1(end,2)/sps,50,'rv','filled');  %end
scatter(ax,line2(1,1)/sps,line2(1,2)/sps,50,'b^','filled');
scatter(ax,line2(end,1)/sps,line2(end,2)/sps,50,'bv','filled');
scatter(ax,line3(1,1)/sps,line3(1,2)/sps,50,'k^','filled');
scatter(ax,line3(end,1)/sps,line3(end,2)/sps,50,'kv','filled');
scatter(ax,line4(1,1)/sps,line4(1,2)/sps,50,'g^','filled');
scatter(ax,line4(end,1)/sps,line4(end,2)/sps,50,'gv','filled');
text(ax,0.98,0.75,sprintf('%d%c',45,char(176)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.98,0.25,sprintf('%d%c',90+45,char(176)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
xlabel(ax,sprintf('\\Delta{t}_{12} (s)'));
ylabel(ax,sprintf('\\Delta{t}_{13} (s)'));
xlim(ax,[-0.15 0.15]);
ylim(ax,[-0.15 0.15]);
xticks(ax,-0.15: 0.05 :0.15);
yticks(ax,-0.15: 0.05 :0.15);

ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on'); axis(ax, 'equal'); 
plot(ax,locline1(:,1),locline1(:,2),'r-','marker','o','markersize',4,'linew',1);
plot(ax,locline2(:,1),locline2(:,2),'b-','marker','o','markersize',4,'linew',1);
plot(ax,locline3(:,1),locline3(:,2),'k-','marker','o','markersize',4,'linew',1);
plot(ax,locline4(:,1),locline4(:,2),'g-','marker','o','markersize',4,'linew',1);
scatter(ax,locline1(1,1),locline1(1,2),50,'r^','filled');  %start
scatter(ax,locline1(end,1),locline1(end,2),50,'rv','filled');  %end
scatter(ax,locline2(1,1),locline2(1,2),50,'b^','filled');
scatter(ax,locline2(end,1),locline2(end,2),50,'bv','filled');
scatter(ax,locline3(1,1),locline3(1,2),50,'k^','filled');
scatter(ax,locline3(end,1),locline3(end,2),50,'kv','filled');
scatter(ax,locline3(1,1),locline3(1,2),50,'k^','filled');
scatter(ax,locline3(end,1),locline3(end,2),50,'kv','filled');
scatter(ax,locline4(1,1),locline4(1,2),50,'g^','filled');
scatter(ax,locline4(end,1),locline4(end,2),50,'gv','filled');
[anggeo,anggeooppo]=angatan2d2geo(round(atan2d(2,3)));
text(ax,0.98,0.75,sprintf('%d%c',anggeo,char(176)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.98,0.25,sprintf('%d%c',180,char(176)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
xlim(ax,[-4 4]);
ylim(ax,[-4 4]);
xticks(ax,-4: 1 :4);
yticks(ax,-4: 1 :4);


