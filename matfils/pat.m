%Sits at one spot (rotations; delays) and looks for the cross-correlation.
clear all
close all
format short e

rotang=-(28./360.)*pi;
IDENTIF4='map2003.063.85.20.80.115.55_2-8-4s';
IDENTIF128='map2003.063.85.20.80.115.55_2-6-128s';
ALL128='all.85.20.128s.fort.13.sor';
OLDER128='map2003.061-062.85.20.80.115.55_2-6-128s';
OLDER=1;
tees(1,:)=[6215 7355];
tees(2,:)=[8901 9091];
tees(3,:)=[10397 10475];
tees(4,:)=[0 86400];
ntees=4;
% % IDENTIF4='map2003.064.85.20.80.115.55_2-6-128s';
% % IDENTIF128='map2003.064.85.20.80.115.55_2-6-128s';
% % ALL128='all.85.20.128s.fort.13.sor';
% % OLDER128='map2003.061-063.85.20.80.115.55_2-6-128s';
% % OLDER=1;
% % tees(1,:)=[0 86400];
% % ntees=1;

% locs4=load(IDENTIF4);
% locs128=load(IDENTIF128);
% locsAll=load(ALL128);
% if OLDER==1
%     olderlocs=load(OLDER128);
% end
% nin4=length(locs4);
% nin128=length(locs128);
orig4=load(IDENTIF4);
orig128=load(IDENTIF128);
origAll=load(ALL128);
if OLDER==1
    origolder=load(OLDER128);
end
nin4=length(orig4);
nin128=length(orig128);

locs4=0.33*orig4(:,2)+0.6i*orig4(:,3);
locs4=locs4*exp(1i*rotang);
locs128=0.33*orig128(:,2)+0.6i*orig128(:,3);
locs128=locs128*exp(1i*rotang);
locsAll=0.33*origAll(:,1)+0.6i*origAll(:,2);
locsAll=locsAll*exp(1i*rotang);
if OLDER==1
    olderlocs=0.33*origolder(:,2)+0.6i*origolder(:,3);
    olderlocs=olderlocs*exp(1i*rotang);
end

for k=1:ntees
    if orig4(1,1) >= tees(k,1)-1
        istart4=1;
    end
    for i=2:nin4-1
        if orig4(i-1,1) < tees(k,1)-1 && orig4(i,1) >= tees(k,1)-1
            istart4=i;
        end
        if orig4(i+1,1) > tees(k,2)+1 && orig4(i,1) <= tees(k,2)+1
            iend4=i;
        end
    end
    if orig4(nin4,1) <= tees(k,2)-1
        iend4=nin4;
    end
    if orig128(1,1) >= tees(k,1)-1
        istart128=1;
    end
    for i=2:nin128-1
        if orig128(i-1,1) < tees(k,1)-1 && orig128(i,1) >= tees(k,1)-1
            istart128=i;
        end
        if orig128(i+1,1) > tees(k,2)+1 && orig128(i,1) <= tees(k,2)+1
            iend128=i;
        end
    end
    if orig128(nin128,1) <= tees(k,2)-1
        iend128=nin128;
    end
    if orig128(1,1) >= tees(k,2)-1
        iend128=istart128;
    end
    if orig128(nin128,1) <= tees(k,1)-1
        istart128=iend128;
    end

    h=figure;
    subplot(2,1,1,'align')
    colormap(jet) 
    plot(real(locsAll),imag(locsAll),'ko','MarkerSize',2,'linewidth',0.2)
    hold on
    plot(real(locs128(1:istart128)),imag(locs128(1:istart128)),'ro','MarkerSize',2,'linewidth',0.2)
    if OLDER==1
        plot(real(olderlocs),imag(olderlocs),'ro','MarkerSize',2,'linewidth',0.2)
    end
    scatter(real(locs128(istart128:iend128)),imag(locs128(istart128:iend128)),90,orig128(istart128:iend128,1)...
        ,'linewidth',1)
    % First on top
    scatter3(real(locs4(istart4:iend4)),imag(locs4(istart4:iend4)),86400-orig4(istart4:iend4,1),40,orig4(istart4:iend4,1),'filled')
    % Last on top
    %scatter3(real(locs4(istart4:iend4)),imag(locs4(istart4:iend4)),orig4(istart4:iend4,1),40,orig4(istart4:iend4,1),'filled')
    axis equal
    axis([-6 5 -6 5])
    xlabel('km E')
    ylabel('km N')
    caxis([tees(k,1) tees(k,2)])
    colorbar
    title(IDENTIF4)
    box on
    
    subplot(2,1,2,'align')
    colormap(jet) 
    plot(real(locsAll),imag(locsAll),'ko','MarkerSize',2,'linewidth',0.2)
    hold on
    plot(real(locs128(1:istart128)),imag(locs128(1:istart128)),'ro','MarkerSize',2,'linewidth',0.2)
    if OLDER==1
        plot(real(olderlocs),imag(olderlocs),'ro','MarkerSize',2,'linewidth',0.2)
    end
    scatter3(real(locs4(istart4:iend4)),imag(locs4(istart4:iend4)),86400-orig4(istart4:iend4,1),40,log(orig4(istart4:iend4,7)),'filled')
    axis equal
    axis([-6 5 -6 5])
    xlabel('km E')
    ylabel('km N')
    caxis([-4.5 1.5])
    %caxis([min(log(orig4(istart4:iend4,7))) max(log(orig4(istart4:iend4,7)))])
    colorbar
    box on
    
    set(h,'PaperPosition',[0.25 0.25 8 10.5])
    print(h,'-depsc',['tmp',int2str(k),'.eps'])
      
end

%Reverse order of last plot (sort-of assuming 0-86400)
h=figure;
subplot(2,1,1,'align')
colormap(jet) 
plot(real(locsAll),imag(locsAll),'ko','MarkerSize',2,'linewidth',0.2)
hold on
plot(real(locs128(1:istart128)),imag(locs128(1:istart128)),'ro','MarkerSize',2,'linewidth',0.2)
if OLDER==1
    plot(real(olderlocs),imag(olderlocs),'ro','MarkerSize',2,'linewidth',0.2)
end
scatter(real(locs128(istart128:iend128)),imag(locs128(istart128:iend128)),90,orig128(istart128:iend128,1)...
    ,'linewidth',1)
% Last on top
scatter3(real(locs4(istart4:iend4)),imag(locs4(istart4:iend4)),orig4(istart4:iend4,1),40,orig4(istart4:iend4,1),'filled')
axis equal
axis([-6 5 -6 5])
xlabel('km E')
ylabel('km N')
caxis([tees(k,1) tees(k,2)])
colorbar
title(IDENTIF4)
box on

subplot(2,1,2,'align')
colormap(jet) 
plot(real(locsAll),imag(locsAll),'ko','MarkerSize',2,'linewidth',0.2)
hold on
plot(real(locs128(1:istart128)),imag(locs128(1:istart128)),'ro','MarkerSize',2,'linewidth',0.2)
if OLDER==1
    plot(real(olderlocs),imag(olderlocs),'ro','MarkerSize',2,'linewidth',0.2)
end
scatter3(real(locs4(istart4:iend4)),imag(locs4(istart4:iend4)),orig4(istart4:iend4,1),40,log(orig4(istart4:iend4,7)),'filled')
axis equal
axis([-6 5 -6 5])
xlabel('km E')
ylabel('km N')
caxis([-4.5 1.5])
colorbar
box on

set(h,'PaperPosition',[0.25 0.25 8 10.5])
print(h,'-depsc',['tmp',int2str(k+1),'.eps'])
      

