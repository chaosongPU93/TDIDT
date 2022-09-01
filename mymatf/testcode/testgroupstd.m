clc
clear
close all

aa = 1:1:50;
na = length(aa);

bb = 1:1:8;
nb = length(bb);
cc = 9:1:20;
nc = length(cc);
dd = 21:1:50;
nd = length(dd);
isequal(na,nb+nc+nd)

meana = mean(aa);
stda = std(aa);

meanb = mean(bb);
stdb = std(bb);
meanc = mean(cc);
stdc = std(cc);
meand = mean(dd);
stdd = std(dd);

meanbcd = (nb*meanb+nc*meanc+nd*meand)/(nb+nc+nd);
stdbcd = sqrt(((nb-1)*stdb^2+nb*(meanb-meanbcd)^2+ ...
               (nc-1)*stdc^2+nc*(meanc-meanbcd)^2+ ...
               (nd-1)*stdd^2+nd*(meand-meanbcd)^2 )/(nb+nc+nd-1));
     
isequal(meana,meanbcd)

isequal(stda,stdbcd)




