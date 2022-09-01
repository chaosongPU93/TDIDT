% testpltlinesegment
clc
close all
clear

figure
ax = gca;
hold on;
xpos = 0;
ypos = 0;
len = 5;
e1 = errorbar(ax,xpos,ypos,len,len,'horizontal','o','markersize',2,'color',...
         'k','linewidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','CapSize',8);
xpos = xpos+len;     
e2 = errorbar(ax,xpos,ypos,len,len,'horizontal','o','markersize',2,'color',...
         'r','linewidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','CapSize',8);
     
xlim(ax,[-10 20]);

print('-dpdf','/home/chaosong/Desktop/linesegment.pdf');