%Sits at one spot (rotations; delays) and looks for the cross-correlation.
clear all
close all
format short e

IDENTIF128='fort.13.85.20.128s.sor';
locs128=load(IDENTIF128);

    figure
    colormap(jet) %PGSI=x, PGSS=y
    scatter(locs128(:,1),locs128(:,2),10,log(locs128(:,3))/2.303,'filled')
    %plot(0.33*locs128(1:istart128,2),0.6*locs128(1:istart128,3),'ko','MarkerSize',2,'linewidth',0.2)
    axis([-16 16 -16 16])
    xlabel('PGSI')
    ylabel('PGSS')
    %caxis([0 70])
    set(gca,'XTick',[-16:1:16]')
    set(gca,'XTickLabel',['   ';'-15';'   ';'   ';'   ';'   ' ...
                               ;'-10';'   ';'   ';'   ';'   ' ...
                               ;' -5';'   ';'   ';'   ';'   ' ...
                               ;' 0 ';'   ';'   ';'   ';'   ' ...
                               ;' 5 ';'   ';'   ';'   ';'   ' ...
                               ;' 10';'   ';'   ';'   ';'   ' ...
                               ;' 15';'   '])
                          
    set(gca,'YTick',[-16:1:16]')
    set(gca,'YTickLabel',['   ';'-15';'   ';'   ';'   ';'   ' ...
                               ;'-10';'   ';'   ';'   ';'   ' ...
                               ;' -5';'   ';'   ';'   ';'   ' ...
                               ;' 0 ';'   ';'   ';'   ';'   ' ...
                               ;' 5 ';'   ';'   ';'   ';'   ' ...
                               ;' 10';'   ';'   ';'   ';'   ' ...
                               ;' 15';'   '])
                          
    colorbar
    %title(IDENTIF4)
    box on
    print('-depsc','count.eps')
    %print('-depsc',['count',int2str(k),'.eps'])
    
box on
