
close all

iflag = 1;  % to control whether to plot diagonal

x1var = linspace(3,5,10);
x2var = linspace(5,3,10);
h = figure(1);
i=1;
x = [1 x1var(i) 7 x2var(i) 1];
y = [1 1 1+2*sqrt(3) 1+2*sqrt(3) 1];
hold on;
plot(x,y,'.-','linewidth',1,'color','b','markersize',14);
if iflag == 1
    plot([1 7],[1 1+2*sqrt(3)],'--','color','b','linewidth',0.8);
    plot([x1var(i) x2var(i)],[1 1+2*sqrt(3)],'--','color','b','linewidth',0.8);
end
plot([1 x1var(i)],[1 1],'-','color','r','linewidth',1);
plot([1 x2var(i)],[1 1+2*sqrt(3)],'-','color','r','linewidth',1);
text(1,0.6,'(1,1)','fontsize',12,'Interpreter','latex');
text(3,0.6,'(3,1)','fontsize',12,'Interpreter','latex');
text(6,4.8,'(7,$1+2\sqrt{3}$)','fontsize',12,'Interpreter','latex');
text(3.5,4.8,'(5,$1+2\sqrt{3}$)','fontsize',12,'Interpreter','latex');
text(2,3,'$l=2\sqrt{7}$','fontsize',14,'Interpreter','latex');
text(2,1.2,'$l=2$','fontsize',14,'Interpreter','latex');
box on;
grid on;
axis equal
axis([0 8 0 6]);
drawnow;
pause(2)
hold off
if iflag == 0
    print('-djpeg','firstone_nodiag.jpg','-r300');
else
    print('-djpeg','firstone.jpg','-r300');
end
frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);   % cm is colormap

clf;
if iflag == 0
    filename = 'dimond_nodiag.gif';
else
    filename = 'dimond.gif';
end

imwrite(imind,cm,filename,'gif','DelayTime',2,'Loopcount',inf);

for i = 2:10
    hold on
    x = [1 x1var(i) 7 x2var(i) 1];
    y = [1 1 1+2*sqrt(3) 1+2*sqrt(3) 1];
    plot(x,y,'.-','linewidth',1,'color','b','markersize',14);
    if iflag == 1
        plot([1 7],[1 1+2*sqrt(3)],'--','color','b','linewidth',0.8);
        plot([x1var(i) x2var(i)],[1 1+2*sqrt(3)],'--','color','b','linewidth',0.8);
    end
    axis equal
    axis([0 8 0 6]);
    box on;
    grid on;
    drawnow    
    pause(1)
    hold off
    
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    clf;
    imwrite(imind,cm,filename,'gif','DelayTime',0.5,'WriteMode','append');
end

i=10;
x = [1 x1var(i) 7 x2var(i) 1];
y = [1 1 1+2*sqrt(3) 1+2*sqrt(3) 1];
plot(x,y,'.-','linewidth',1,'markersize',14);
hold on;
plot(x,y,'.-','color','b','linewidth',1,'markersize',14);
if iflag == 1
    plot([1 7],[1 1+2*sqrt(3)],'--','color','b','linewidth',0.8);
    plot([x1var(i) x2var(i)],[1 1+2*sqrt(3)],'--','color','b','linewidth',0.8);
end
plot([1 x1var(i)],[1 1],'-','color','r','linewidth',1);
plot([1 x2var(i)],[1 1+2*sqrt(3)],'-','color','r','linewidth',1);
text(1,0.6,'(1,1)','fontsize',12,'Interpreter','latex');
text(5,0.6,'(5,1)','fontsize',12,'Interpreter','latex');
text(6,4.8,'(7,$1+2\sqrt{3}$)','fontsize',12,'Interpreter','latex');
text(1.5,4.8,'(3,$1+2\sqrt{3}$)','fontsize',12,'Interpreter','latex');
text(1,3,'$l=4$','fontsize',14,'Interpreter','latex');
text(3,1.2,'$l=4$','fontsize',14,'Interpreter','latex');
axis equal
axis([0 8 0 6]);
box on;
grid on;
drawnow
pause(2)
if iflag == 0
    print('-djpeg','lastone_nodiag.jpg','-r300');
else
    print('-djpeg','lastone.jpg','-r300');
end


frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);

imwrite(imind,cm,filename,'gif','DelayTime',3,'WriteMode','append');