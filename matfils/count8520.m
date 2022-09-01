%Sits at one spot (rotations; delays) and looks for the cross-correlation.
clear all
close all
format short e

rotang=-(28./360.)*pi;
IDENTIF128='fort.13.85.20.128s.sor';
locs=load(IDENTIF128);
xloc=locs(:,1);
yloc=locs(:,2);
locs128=0.33*xloc+0.6i*yloc;
locs128=locs128*exp(1i*rotang);
    figure
    colormap(jet) %PGSI=x, PGSS=y
    scatter(real(locs128),imag(locs128),10,log(locs(:,3))/2.303,'filled')
    %plot(0.33*locs128(1:istart128,2),0.6*locs128(1:istart128,3),'ko','MarkerSize',2,'linewidth',0.2)
    axis([-6 5 -6 5])
%     xlabel('PGSI')
%     ylabel('PGSS')
%     set(gca,'XTick',[-16:1:16]')
%     set(gca,'XTickLabel',['   ';'-15';'   ';'   ';'   ';'   ' ...
%                                ;'-10';'   ';'   ';'   ';'   ' ...
%                                ;' -5';'   ';'   ';'   ';'   ' ...
%                                ;' 0 ';'   ';'   ';'   ';'   ' ...
%                                ;' 5 ';'   ';'   ';'   ';'   ' ...
%                                ;' 10';'   ';'   ';'   ';'   ' ...
%                                ;' 15';'   '])
%                           
%     set(gca,'YTick',[-16:1:16]')
%     set(gca,'YTickLabel',['   ';'-15';'   ';'   ';'   ';'   ' ...
%                                ;'-10';'   ';'   ';'   ';'   ' ...
%                                ;' -5';'   ';'   ';'   ';'   ' ...
%                                ;' 0 ';'   ';'   ';'   ';'   ' ...
%                                ;' 5 ';'   ';'   ';'   ';'   ' ...
%                                ;' 10';'   ';'   ';'   ';'   ' ...
%                                ;' 15';'   '])
                          
    colorbar
    %title(IDENTIF4)
    box on
    print('-depsc','count8520.eps')
    %print('-depsc',['count',int2str(k),'.eps'])
    
box on
