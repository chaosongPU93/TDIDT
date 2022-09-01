%Sits at one spot (rotations; delays) and looks for the cross-correlation.
clear all
close all
format short e

% IDENTIF4='map2005.259.36.75.80.365.385_2-6-4s';
% IDENTIF128='map2005.259.36.75.80.365.385_2-6-128s';
% ALL128='all.36.75.128s';
% OLDER=0;
% tees(1,:)=[44000 45800];
% tees(2,:)=[45800 47800];
% tees(3,:)=[47800 50000];
% tees(4,:)=[50000 50800];
% tees(5,:)=[50800 53000];
% tees(6,:)=[54000 55000];
% tees(7,:)=[56000 57000];
% tees(8,:)=[57000 58500];
% tees(9,:)=[60000 61500];
% tees(10,:)=[67000 69000];
% tees(11,:)=[69000 71000];
% tees(12,:)=[44000 86000];
% ntees=12;
% % IDENTIF4='map2005.260.36.75.80.365.385_2-6-4s';
% % IDENTIF128='map2005.260.36.75.80.365.385_2-6-128s';
% % ALL128='all.36.75.128s';
% % OLDER128='map2005.259.36.75.80.365.385_2-6-128s';
% % OLDER=1;
% % tees(1,:)=[0 3000];
% % tees(2,:)=[21000 24000];
% % tees(3,:)=[0 24000];
% % ntees=3;

IDENTIF4='map2005.259.45.72.65.355.370_2-6-4s';
IDENTIF128='map2005.259.45.72.65.355.370_2-6-128s';
ALL128='map2005.259.45.72.65.355.370_2-6-128s';
OLDER=0;
tees(1,:)=[8000 10000];
tees(2,:)=[10000 12000];
tees(3,:)=[12000 14100];
tees(4,:)=[14100 16000];
tees(5,:)=[16000 18000];
tees(6,:)=[20000 21000];
tees(7,:)=[22000 25000];
tees(8,:)=[26000 27000];
tees(9,:)=[28000 29000];
tees(10,:)=[29500 32500];
tees(11,:)=[33500 37000];
tees(12,:)=[40500 41700];
tees(13,:)=[44000 45500];
tees(14,:)=[45500 48000];
tees(15,:)=[48000 54000];
tees(16,:)=[58000 59000];
tees(17,:)=[62000 65000];
tees(18,:)=[66000 67500];
tees(19,:)=[68000 74000];
tees(20,:)=[2000 74000];
ntees=20;

locs4=load(IDENTIF4);
locs128=load(IDENTIF128);
locsAll=load(ALL128);
if OLDER==1
    olderlocs=load(OLDER128);
end
nin4=length(locs4);
nin128=length(locs128);

scrsz=get(0,'ScreenSize');
for k=ntees:-1:1
    if locs4(1,1) >= tees(k,1)-1
        istart4=1;
    end
    for i=2:nin4-1
        if locs4(i-1,1) < tees(k,1)-1 && locs4(i,1) >= tees(k,1)-1
            istart4=i;
        end
        if locs4(i+1,1) > tees(k,2)+1 && locs4(i,1) <= tees(k,2)+1
            iend4=i;
        end
    end
    if locs4(nin4,1) <= tees(k,2)-1
        iend4=nin4-1;
    end
    if locs128(1,1) >= tees(k,1)-1
        istart128=1;
    end
    for i=2:nin128-1
        if locs128(i-1,1) < tees(k,1)-1 && locs128(i,1) >= tees(k,1)-1
            istart128=i;
        end
        if locs128(i+1,1) > tees(k,2)+1 && locs128(i,1) <= tees(k,2)+1
            iend128=i;
        end
    end
    if locs128(nin128,1) <= tees(k,2)-1
        iend128=nin128;
    end
    if locs128(1,1) >= tees(k,2)-1
        iend128=istart128;
    end
    if locs128(nin128,1) <= tees(k,1)-1
        istart128=iend128;
    end
  
    h=figure('Position',[scrsz(3)/3.1 1 scrsz(3)/3.1 scrsz(4)]);
    subplot(2,1,1,'align')
    colormap(jet) %PGSI=x, PGSS=y
    plot(locsAll(:,2),locsAll(:,3),'ko','MarkerSize',2,'linewidth',0.2)
    hold on
    plot(locs128(1:istart128,2),locs128(1:istart128,3),'ro','MarkerSize',2,'linewidth',0.2)
    if OLDER==1
        plot(olderlocs(:,2),olderlocs(:,3),'ro','MarkerSize',2,'linewidth',0.2)
    end
    scatter(locs128(istart128:iend128,2),locs128(istart128:iend128,3),90,locs128(istart128:iend128,1),'linewidth',1)
    % First on top
    scatter3(locs4(istart4:iend4,2),locs4(istart4:iend4,3),86400-locs4(istart4:iend4,1),25,locs4(istart4:iend4,1),'filled')
    axis equal
    axis([-11 11 -11 11])
    xlabel('PGSI')
    ylabel('PGSS')
    caxis([tees(k,1) tees(k,2)])
    colorbar
    title(IDENTIF4)
    box on

    subplot(2,1,2,'align')
    colormap(jet) %PGSI=x, PGSS=y
    plot(locsAll(:,2),locsAll(:,3),'ko','MarkerSize',2,'linewidth',0.2)
    hold on
    plot(locs128(1:istart128,2),locs128(1:istart128,3),'ro','MarkerSize',2,'linewidth',0.2)
    if OLDER==1
        plot(olderlocs(:,2),olderlocs(:,3),'ro','MarkerSize',2,'linewidth',0.2)
    end
    scatter3(locs4(istart4:iend4,2),locs4(istart4:iend4,3),86400-locs4(istart4:iend4,1),40,log(locs4(istart4:iend4,7)),'filled')
    axis equal
    axis([-11 11 -11 11])
    xlabel('PGSI')
    ylabel('PGSS')
    caxis([-5 1])
    colorbar
    box on
    
    set(h,'PaperPosition',[0.25 0.25 8 10.5])
    print(h,'-depsc',['tmp',int2str(k),'.eps'])
    
end

%plot whole day (or last window)
k=ntees;
if locs4(1,1) >= tees(k,1)-1
    istart4=1;
end
for i=2:nin4-1
    if locs4(i-1,1) < tees(k,1)-1 && locs4(i,1) >= tees(k,1)-1
        istart4=i;
    end
    if locs4(i+1,1) > tees(k,2)+1 && locs4(i,1) <= tees(k,2)+1
        iend4=i;
    end
end
if locs4(nin4,1) <= tees(k,2)-1
    iend4=nin4;
end
if locs128(1,1) >= tees(k,1)-1
    istart128=1;
end
for i=2:nin128-1
    if locs128(i-1,1) < tees(k,1)-1 && locs128(i,1) >= tees(k,1)-1
        istart128=i;
    end
    if locs128(i+1,1) > tees(k,2)+1 && locs128(i,1) <= tees(k,2)+1
        iend128=i;
    end
end
if locs128(nin128,1) <= tees(k,2)-1
    iend128=nin128;
end
if locs128(1,1) >= tees(k,2)-1
    iend128=istart128;
end
if locs128(nin128,1) <= tees(k,1)-1
    istart128=iend128;
end

h=figure('Position',[1 1 scrsz(3)/3.1 scrsz(4)]);
subplot(2,1,1,'align')
colormap(jet) 
plot(locsAll(:,2),locsAll(:,3),'ko','MarkerSize',2,'linewidth',0.2)
hold on
plot(locs128(1:istart128,2),locs128(1:istart128,3),'ro','MarkerSize',2,'linewidth',0.2)
if OLDER==1
    plot(olderlocs(:,2),olderlocs(:,3),'ro','MarkerSize',2,'linewidth',0.2)
end
% Last on top
scatter3(locs4(istart4:iend4,2),locs4(istart4:iend4,3),86400-locs4(istart4:iend4,1),25,locs4(istart4:iend4,1),'filled')
axis equal
axis([-11 11 -11 11])
xlabel('PGSI')
ylabel('PGSS')
caxis([tees(k,1) tees(k,2)])
colorbar
title(IDENTIF4)
box on
subplot(2,1,2,'align')
colormap(jet) 
plot(locsAll(:,2),locsAll(:,3),'ko','MarkerSize',2,'linewidth',0.2)
hold on
plot(locs128(1:istart128,2),locs128(1:istart128,3),'ro','MarkerSize',2,'linewidth',0.2)
if OLDER==1
    plot(olderlocs(:,2),olderlocs(:,3),'ro','MarkerSize',2,'linewidth',0.2)
end
scatter3(locs4(istart4:iend4,2),locs4(istart4:iend4,3),86400-locs4(istart4:iend4,1),25,log(locs4(istart4:iend4,7)),'filled')
axis equal
axis([-11 11 -11 11])
xlabel('PGSI')
ylabel('PGSS')
caxis([-5 1])
colorbar
box on
set(h,'PaperPosition',[0.25 0.25 8 10.5])
print(gcf,'-depsc',['tmp',int2str(k+1),'.eps'],'-zbuffer')

%Reverse order of same
h=figure('Position',[2*scrsz(3)/3 1 scrsz(3)/3 scrsz(4)]);
subplot(2,1,1,'align')
colormap(jet) 
plot(locsAll(:,2),locsAll(:,3),'ko','MarkerSize',2,'linewidth',0.2)
hold on
plot(locs128(1:istart128,2),locs128(1:istart128,3),'ro','MarkerSize',2,'linewidth',0.2)
if OLDER==1
    plot(olderlocs(:,2),olderlocs(:,3),'ro','MarkerSize',2,'linewidth',0.2)
end
% Last on top
scatter3(locs4(istart4:iend4,2),locs4(istart4:iend4,3),locs4(istart4:iend4,1),25,locs4(istart4:iend4,1),'filled')
axis equal
axis([-11 11 -11 11])
xlabel('PGSI')
ylabel('PGSS')
caxis([tees(k,1) tees(k,2)])
colorbar
title(IDENTIF4)
box on
subplot(2,1,2,'align')
colormap(jet) 
plot(locsAll(:,2),locsAll(:,3),'ko','MarkerSize',2,'linewidth',0.2)
hold on
plot(locs128(1:istart128,2),locs128(1:istart128,3),'ro','MarkerSize',2,'linewidth',0.2)
if OLDER==1
    plot(olderlocs(:,2),olderlocs(:,3),'ro','MarkerSize',2,'linewidth',0.2)
end
scatter3(locs4(istart4:iend4,2),locs4(istart4:iend4,3),locs4(istart4:iend4,1),25,log(locs4(istart4:iend4,7)),'filled')
axis equal
axis([-11 11 -11 11])
xlabel('PGSI')
ylabel('PGSS')
caxis([-5 1])
colorbar
box on
set(h,'PaperPosition',[0.25 0.25 8 10.5])
print(gcf,'-depsc',['tmp',int2str(k+2),'.eps'],'-zbuffer')
%print(h,'-depsc', ['tmp',int2str(k+1),'.eps'])
%print(h,'-djpeg','-r600',['tmp',int2str(k+1),'.jpg'])
%print(h,'-dpdf', ['tmp',int2str(k+1),'.pdf'])
     
