%
format short e
clear all
close all
set(0,'DefaultFigureVisible','on');

days=[2003.062;
      2003.063;
      2003.064;
      2004.196;
      2004.197;
      2004.198;
      2004.199;
      2005.254;
      2005.255;
      2005.256;
      2005.257];
ndays=length(days);

% %001 030228 6 3147.900 1.726 6
% %001 030302 17 1466.150 1.694 10
% %001 030302 23 1018.875 1.611 8
  % bostocks=load('BOSTOCK/NEW/total_mag_detect_0000.txt');
  % bostsec=3600*(bostocks(:,3)-1)+bostocks(:,4)+4*0.025; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MO=day2month(jday,year);

scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load an Armbruster file.
%   3   2  23  20  43   8   48.2330  123.7413   31.5900   19.0587  -29.6890   54.8633
%   3   2  23  23  39  28   48.8600  122.8453   61.9200   85.0759   40.0304   54.9857
%   3   2  24  13  31  43   48.0175  123.2828   35.9500   52.8408  -53.6516   55.5637
Armb=load('/data1/arubin.old/Armbruster/40sps.pgsssi');
totin=0; %sum over all days
for iday=1:ndays
    yr=floor(days(iday)-2000);
    jul=round(1000*(days(iday)-floor(days(iday))));
    if yr==4
        jul=jul+364; %Should be 365 but I think 364 is what works;
    elseif yr==5
        jul=jul+365+366;
    end
    wins=windows(days(iday));
    totfrac=sum(wins(:,2)-wins(:,1))/86400
    nwins=size(wins,1);
    i=1;
    while(Armb(i,1) < yr)
        i=i+1;
    end
    while(floor(Armb(i,12)) < jul)
        i=i+1;
    end
    evin=0; %sum over all windows for this day
    for iwin=1:nwins
        seconds=86400*(Armb(i,12)-floor(Armb(i,12)));
        while(seconds < wins(iwin,1))
            i=i+1;
            seconds=86400*(Armb(i,12)-floor(Armb(i,12)));
        end
        while(floor(Armb(i,12))==jul && seconds < wins(iwin,2))
            evin=evin+1;
            Armbin(evin,:)=[Armb(i,:),seconds];
            i=i+1;
            seconds=86400*(Armb(i,12)-floor(Armb(i,12)));            
        end
    end
    if evin > 0
        ArmAllin(totin+1:totin+evin,:)=Armbin;
        totin=totin+evin;
        h=figure('Position',[wid/3 1 1.5*wid hite]); 
        scatter(Armbin(:,10),Armbin(:,11),3,Armbin(:,13))
        title(num2str(days(iday)))
        colorbar
        axis equal
    end
    clear Armbin %Need to do this since it's not initialized and size may vary with day.
end
fid = fopen('CATALOGS/ArmbrusterIns','w');
fprintf(fid,'%4i %4i %4i %4i %4i %4i %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %11.4f\n', ArmAllin');
fclose(fid);
%    set(h,'PaperPosition',[0.25 0.25 8 10.5])
%    print(h,'-depsc',['ARMMAP/WIGS/',YEAR,'.',JDAY,'.',int2str(timsPGC(istart)),'_',num2str(lo),'-',num2str(hi),'eye.eps'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Yajun/s catalog
ytime=[2003,60,67;
       2004,194,203;
       2005,254,261];
iyold=0;
totin=0; %sum over all days
for iday=1:ndays
    yr=floor(days(iday)-2000);
    jul=round(1000*(days(iday)-floor(days(iday))));
    wins=windows(days(iday));
    nwins=size(wins,1);
    iy=yr-2;
    wins=wins+(jul-ytime(iy,2))*86400;  %Yajun lists his times as seconds since "day 1"
    if iy ~= iyold
        load(['/data2/arubin/YAJUN/fcna_',num2str(ytime(iy,1)),'_',num2str(ytime(iy,2)),'_',num2str(ytime(iy,3)),'_8km_3.mat'])
        i=1;  %initilize only when reading in a new file
    end
    iyold=iy;
    evin=0;  %sum over all windows for this day
    for iwin=1:nwins
        while(fcna(i,7) < wins(iwin,1))
            i=i+1;
        end
        while(fcna(i,7) < wins(iwin,2))
            evin=evin+1;
            Pengin(evin,:)=fcna(i,:);
            i=i+1;
        end
    end
    if evin > 0
        PenAllin(totin+1:totin+evin,:)=Pengin;
        totin=totin+evin;
        h=figure('Position',[wid/3 1 1.5*wid hite]); 
        scatter(Pengin(:,22),Pengin(:,21),3,mod(Pengin(:,7),86400))
        title(num2str(days(iday)))
        colorbar
        axis equal
    end
    clear Pengin %Need to do this since it's not initialized and size may vary with day.
end
PenAllinDec(:,1)=PenAllin(:,1);  %Family
PenAllinDec(:,2)=PenAllin(:,22);  %kmE
PenAllinDec(:,3)=PenAllin(:,21);  %kmN
PenAllinDec(:,4)=PenAllin(:,7);  %Yajun's times
PenAllinDec(:,5)=mod(PenAllin(:,7),86400);  %daily times
PenAllinDec(:,6)=PenAllin(:,16);  %amplitudes
fid = fopen('CATALOGS/PengIns','w');
fprintf(fid,'%5i %10.4f %10.4f %11.4e\n', PenAllinDec');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Michael's catalog
% 267 030303 1 16.050 1.493 13
% 019 030303 1 134.325 1.689 9
% 019 030303 1 144.175 1.649 16
Bost=load('/data2/arubin/CNDC/BOSTOCK/NEW/total_mag_detect_0000.2003-2005');
totin=0; %sum over all days
for iday=1:ndays
    yr=floor(days(iday)-2000);
    jul=round(1000*(days(iday)-floor(days(iday))));
    yrmodom=jul2yrmodom(jul,yr);
    wins=windows(days(iday));
    nwins=size(wins,1);
    i=1;
    while(Bost(i,2) < yrmodom)
        i=i+1;
    end
    evin=0; %sum over all windows for this day
    for iwin=1:nwins
        seconds=3600*(Bost(i,3)-1)+Bost(i,4);
        while(seconds < wins(iwin,1))
            i=i+1;
            seconds=3600*(Bost(i,3)-1)+Bost(i,4);
        end
        while(seconds < wins(iwin,2))
            evin=evin+1;
            Bostin(evin,:)=[Bost(i,1:5),seconds];
            i=i+1;
            seconds=3600*(Bost(i,3)-1)+Bost(i,4);
        end
    end
    if evin > 0
        BosAllin(totin+1:totin+evin,:)=Bostin;
        totin=totin+evin;
%         h=figure('Position',[wid/3 1 1.5*wid hite]); 
%         scatter(Armbin(:,10),Armbin(:,11),3,Armbin(:,13))
%         title(num2str(days(iday)))
%         colorbar
%         axis equal
    end
    clear Bostin %Need to do this since it's not initialized and size may vary with day.
end
fid = fopen('CATALOGS/BostockIns','w');
fprintf(fid,'%4i %7i %3i %10.3f %7.3f %11.3f\n', BosAllin');
fclose(fid);
