%Sits at one spot (rotations; delays) and looks for the cross-correlation.
%This is for refined views of propagating bursts, after viewing mapgram.
clear all
close all
format short e

rotmap=(-28./180.)*pi; %this is for map view
rotprop=(22./180.)*pi; %this is for projections
dip=[0 -6; 0 6];
rotdip=((68-28)/180.)*pi;

% IDENTIF4='map2003.061.85.20.80.115.55_2-8-4s';
% IDENTIF128='map2003.061.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER=0;
% tees(1,:)=[0 12000];
% tees(2,:)=[12000 21000];
% tees(3,:)=[21000 25000];
% tees(4,:)=[25000 32000];
% tees(5,:)=[32000 53000];
% tees(6,:)=[53000 61800];
% tees(7,:)=[61800 72000];
% tees(8,:)=[72000 82000];
% tees(9,:)=[82000 86400];
% tees(10,:)=[21000 86400];
% ntees=10;
% % IDENTIF4='map2003.062.85.20.80.115.55_2-8-4s';
% % IDENTIF128='map2003.062.85.20.80.115.55_2-6-128s';
% % ALL128='all.85.20.128s.fort.13.sor';
% % OLDER128='map2003.061.85.20.80.115.55_2-6-128s';
% % OLDER=1;
% % tees(1,:)=[16099 16492];
% % tees(2,:)=[16650 16750];
% % tees(3,:)=[17830 17910];
% % tees(4,:)=[18787 19195];
% % tees(5,:)=[22147 22595];
% % tees(6,:)=[23000 23200];
% % tees(7,:)=[25091 25500];
% % tees(8,:)=[31360 31430];
% % tees(9,:)=[32491 32900];
% % tees(10,:)=[34131 34459];
% % %tees(18,:)=[35691 36227]; %diff
% % tees(11,:)=[38339 38891];
% % tees(12,:)=[40597 41243];
% % tees(13,:)=[41300 41500];
% % tees(14,:)=[43129 43450];
% % tees(15,:)=[46347 46763];
% % tees(16,:)=[52250 52590];
% % tees(17,:)=[54920 55115];
% % tees(18,:)=[57610 58020];
% % tees(19,:)=[61595 61753];
% % tees(20,:)=[68403 69011];
% % tees(21,:)=[72147 72443];
% % tees(22,:)=[75060 76400];
% % tees(23,:)=[76971 77683];
% % tees(24,:)=[81943 82776];
% % tees(25,:)=[84400 85200];
% % ntees=25;
% IDENTIF4='map2003.063.85.20.80.115.55_2-8-4s';
% IDENTIF128='map2003.063.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER128='map2003.061-062.85.20.80.115.55_2-6-128s';
% OLDER=1;
% tees(1,:)=[6100 7400];
% tees(2,:)=[12300 12550];
% tees(3,:)=[21403 22576];
% tees(4,:)=[31853 32411];
% tees(5,:)=[32025 33600];
% tees(6,:)=[40600 42166];
% tees(7,:)=[42850 43000];
% tees(8,:)=[64935 65534];
% tees(9,:)=[65590 65760];
% tees(10,:)=[79980 81280];
% ntees=10;

% IDENTIF4='map2004.196.85.20.80.115.55_2-6-4s';
% IDENTIF128='map2004.196.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER=0;
% tees(1,:)=[40899 41395];
% tees(2,:)=[41645 42450];
% tees(3,:)=[42900 43187];
% tees(4,:)=[44227 44539];
% tees(5,:)=[47080 47140];
% tees(6,:)=[48375 49034];
% tees(7,:)=[51035 51600];
% tees(8,:)=[56075 56643];
% tees(9,:)=[56819 57603];
% tees(10,:)=[62576 63683];
% tees(11,:)=[64807 65280];
% tees(12,:)=[70395 70854];
% tees(13,:)=[73855 74566];
% tees(14,:)=[85841 86379];
% ntees=14;
% % IDENTIF4='map2004.197.85.20.80.115.55_2-6-4s';
% % IDENTIF128='map2004.197.85.20.80.115.55_2-6-128s';
% % ALL128='all.85.20.128s.fort.13.sor';
% % OLDER128='map2004.196.85.20.80.115.55_2-6-128s';
% % OLDER=1;
% % tees(1,:)=[15883 16475];
% % tees(2,:)=[24043 24427];
% % tees(3,:)=[25227 25900];
% % tees(4,:)=[27431 28673];
% % tees(5,:)=[30700 31550]; 
% % tees(6,:)=[36850 38587]; 
% % tees(7,:)=[39195 39600];
% % tees(8,:)=[42959 43600];
% % tees(9,:)=[45940 46020];
% % tees(10,:)=[46150 46450];
% % tees(11,:)=[69603 70600];
% % tees(12,:)=[78295 78850]; %79330];
% % tees(13,:)=[84379 85059]; %possible split?
% % ntees=13;
% IDENTIF4='map2004.198.85.20.80.115.55_2-6-4s';
% IDENTIF128='map2004.198.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER128='map2004.196-197.85.20.80.115.55_2-6-128s';
% OLDER=1;
% tees(1,:)=[9200 11000];
% tees(2,:)=[11000 12000];
% tees(3,:)=[0 86400];   
% ntees=3;
% % IDENTIF4='map2004.199.85.20.80.115.55_2-6-4s';
% % IDENTIF128='map2004.198.85.20.80.115.55_2-6-128s';
% % ALL128='all.85.20.128s.fort.13.sor';
% % OLDER128='map2004.196-197.85.20.80.115.55_2-6-128s';
% % OLDER=1;
% % tees(1,:)=[0 2000];
% % tees(2,:)=[2000 4000];
% % tees(3,:)=[4200 4300];   
% % tees(4,:)=[4750 5000];   
% % tees(5,:)=[6150 6300];   
% % tees(6,:)=[6300 8000];   
% % tees(7,:)=[8550 10000];   
% % tees(8,:)=[0 10000];   
% % ntees=8;

% % IDENTIF4='map2005.253.85.20.80.115.55_2-6-4s';
% % IDENTIF128='map2005.253.85.20.80.115.55_2-6-128s';
% % ALL128='all.85.20.128s.fort.13.sor';
% % OLDER=0;
% % tees(1,:)=[0 86400];
% % ntees=1;
% IDENTIF4='map2005.254.85.20.80.115.55_2-6-4s';
% IDENTIF128='map2005.254.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER=0;
% tees(1,:)=[14000 18000];
% tees(2,:)=[18090 18261];
% tees(3,:)=[21106 21348];
% tees(4,:)=[24330 24596];
% tees(5,:)=[27426 28044];
% tees(6,:)=[30342 31676];
% tees(7,:)=[30342 31140];
% tees(8,:)=[31434 31676];
% tees(9,:)=[34498 34942];
% tees(10,:)=[35202 35916];
% tees(11,:)=[48754 49932];
% tees(12,:)=[54226 54972];
% tees(13,:)=[73372 73788];
% tees(14,:)=[79970 81020];
% ntees=14;
IDENTIF4='map2005.255.85.20.80.115.55_2-6-4s';
IDENTIF128='map2005.255.85.20.80.115.55_2-6-128s';
ALL128='all.85.20.128s.fort.13.sor';
OLDER128='map2005.254.85.20.80.115.55_2-6-128s';
OLDER=1;
tees(1,:)=[4835 6060];
tees(2,:)=[22739 23395];
tees(3,:)=[34451 36408];
tees(4,:)=[43013 43867];
tees(5,:)=[45971 46649];
tees(6,:)=[50523 51796];
tees(7,:)=[54109 54979];
tees(8,:)=[58019 59571];
tees(9,:)=[61859 62587];
tees(10,:)=[63669 63979];
tees(11,:)=[67633 68595];
ntees=11;
% IDENTIF4='map2005.256.85.20.80.115.55_2-6-4s';
% IDENTIF128='map2005.256.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER128='map2005.255.85.20.80.115.55_2-6-128s';
% OLDER=1;
% tees(1,:)=[3000 4500];
% tees(2,:)=[3000 5000];
% tees(3,:)=[12000 16000];
% tees(4,:)=[16000 17200];
% tees(5,:)=[79000 83000];
% tees(6,:)=[80100 80200];
% tees(7,:)=[84000 85900];
% tees(8,:)=[85000 85200];
% tees(9,:)=[0 86400];
% ntees=9;
% % IDENTIF4='map2005.257.85.20.80.115.55_2-6-4s';
% % IDENTIF128='map2005.257.85.20.80.115.55_2-6-128s';
% % ALL128='all.85.20.128s.fort.13.sor';
% % OLDER128='map2005.255.85.20.80.115.55_2-6-128s';
% % OLDER=1;
% % tees(1,:)=[2650 2850];
% % tees(2,:)=[82000 86000];
% % tees(3,:)=[0 86400];
% % ntees=3;

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
locs4=locs4*exp(1i*rotmap);
locs128=0.33*orig128(:,2)+0.6i*orig128(:,3);
locs128=locs128*exp(1i*rotmap);
locsAll=0.33*origAll(:,1)+0.6i*origAll(:,2);
locsAll=locsAll*exp(1i*rotmap);
if OLDER==1
    olderlocs=0.33*origolder(:,2)+0.6i*origolder(:,3);
    olderlocs=olderlocs*exp(1i*rotmap);
end
prop4=0.33*orig4(:,2)+0.6i*orig4(:,3);
prop4=prop4*exp(1i*rotprop);
realprop4=real(prop4);
imagprop4=imag(prop4);
prop128=0.33*orig128(:,2)+0.6i*orig128(:,3);
prop128=prop128*exp(1i*rotprop);
propAll=0.33*origAll(:,1)+0.6i*origAll(:,2);
propAll=propAll*exp(1i*rotprop);
if OLDER==1
    olderprop=0.33*origolder(:,2)+0.6i*origolder(:,3);
    olderprop=olderprop*exp(1i*rotprop);
end
dip=dip*exp(1i*rotdip);

clo=-5;
chi=2;
dotsz=15;
scrsz=get(0,'ScreenSize');
wid=scrsz(3)/2.1;
hite=scrsz(4);
scrat=wid/hite;
fidprop = fopen(['PROP/',IDENTIF4(4:length(IDENTIF4))],'w');
fprintf(fidprop,'k  time  medsize  s(1)d(2)f(0)  vel(m/s)  dstrike1 dstrike2 ddip1 ddip2  tau \n');
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

    h=figure('Position',[wid 1 wid hite]);
    %subplot(2,1,1,'align')
    axes('Units','normalized','position',[0.1 0.55 0.5 0.5])
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
    scatter3(real(locs4(istart4:iend4)),imag(locs4(istart4:iend4)),86400-orig4(istart4:iend4,1),25,orig4(istart4:iend4,1),'filled')
    % Last on top
    %scatter3(real(locs4(istart4:iend4)),imag(locs4(istart4:iend4)),orig4(istart4:iend4,1),40,orig4(istart4:iend4,1),'filled')
    plot(dip,'k-')
    axis equal
    axis([-6 5 -6 5])
    xlabel('km E')
    ylabel('km N')
    caxis([tees(k,1) tees(k,2)])
    colorbar
    title(IDENTIF4)
    box on
    
    %subplot(2,1,2,'align')
    axes('position',[0.1 0.2 0.5 0.5])
    colormap(jet) 
    plot(real(locsAll),imag(locsAll),'ko','MarkerSize',2,'linewidth',0.2)
    hold on
    plot(real(locs128(1:istart128)),imag(locs128(1:istart128)),'ro','MarkerSize',2,'linewidth',0.2)
    if OLDER==1
        plot(real(olderlocs),imag(olderlocs),'ro','MarkerSize',2,'linewidth',0.2)
    end
    scatter3(real(locs4(istart4:iend4)),imag(locs4(istart4:iend4)),86400-orig4(istart4:iend4,1),25,log(orig4(istart4:iend4,7)),'filled')
    axis equal
    axis([-6 5 -6 5])
    xlabel('km E')
    ylabel('km N')
    caxis([clo chi])
    %caxis([-4.5 2.5])
    %caxis([min(log(orig4(istart4:iend4,7))) max(log(orig4(istart4:iend4,7)))])
    colorbar
    box on
    
%     axes('position',[0.1 0.15 0.8 0.08])
%     hold on
%     hrf = plotreflinesr(gca,orig4(:,1),'x','r');
%     plot(timsPGC,real(PGCrot),'g');
%     plot(timsPGC,real(SSIBrot),'b');
%     plot(timsPGC,real(SILBrot),'k');
%     axis([tees(k,1) tees(k,2) ltmin ltmax]);
%     %xlabel(WAVES)
%     box on
    
    axes('position',[0.68 0.6 0.3 0.3*scrat])
    colormap(jet) 
    hold on
    scatter3(orig4(istart4:iend4,1),realprop4(istart4:iend4),86400-orig4(istart4:iend4,1),dotsz,log(orig4(istart4:iend4,7)), ...
        'filled')
    medsiz=median(log(orig4(istart4:iend4,7)));
    title(num2str(medsiz))
    %text(0.5*(tees(k,2)+tees(k,1)), 4, num2str(median(log(orig4(istart4:iend4,7)))),'fontsize',8);
    %axis([tees(k,1) tees(k,2) -6 5])
    ylabel('km along strike')
    caxis([clo chi])
    %colorbar
    box on
    hold on
    [tstrike,xstrike,bstrike]=ginput; %do this left to right.  "Enter" after 2 or 0 points
    if length(tstrike)==2
        plot(tstrike,xstrike)
    end
    [dummy,dstrike]=ginput; %do this low to high, regardless
    if length(dstrike)==2
        line1=[dummy(1) dstrike(1); dummy(1)+(tees(k,2)-tees(k,1)) dstrike(1)];
        line2=[dummy(1) dstrike(2); dummy(1)+(tees(k,2)-tees(k,1)) dstrike(2)];
        plot(line1(:,1),line1(:,2))
        plot(line2(:,1),line2(:,2))
    end

    axes('position',[0.68 0.25 0.3 0.3*scrat])
    colormap(jet) 
    hold on
    scatter3(orig4(istart4:iend4,1),imagprop4(istart4:iend4),86400-orig4(istart4:iend4,1),dotsz,log(orig4(istart4:iend4,7)), ...
        'filled')
    %axis([tees(k,1) tees(k,2) -6 5])
    ylabel('km along dip')
    caxis([clo chi])
    %colorbar
    box on
    hold on
    [tdip,xdip,bdip]=ginput; %do this left to right.  "Enter" after 2 or 0 points
    if length(tdip)==2
        plot(tdip,xdip)
    end
    [dummy,ddip]=ginput; %do this low to high, regardless
    if length(ddip)==2
        line1=[dummy(1) ddip(1); dummy(1)+(tees(k,2)-tees(k,1)) ddip(1)];
        line2=[dummy(1) ddip(2); dummy(1)+(tees(k,2)-tees(k,1)) ddip(2)];
        plot(line1(:,1),line1(:,2))
        plot(line2(:,1),line2(:,2))
    end
    
    datfile(1)=k;
    tim=0.5*(tees(k,1)+tees(k,2));
    datfile(2)=tim;
    datfile(3)=medsiz;
    if length(tstrike)==2
        datfile(4)=1;
        vel=(xstrike(1)-xstrike(2))/(tstrike(2)-tstrike(1)); %Should be negative for reversals
    elseif length(tdip)==2
        datfile(4)=2;
        vel=(xdip(1)-xdip(2))/(tdip(2)-tdip(1)); %Should be negative for reversals
    else
        datfile(4)=0;
        vel=0;
    end
    datfile(5)=vel*1000; %(in m/s)
%     strikelen=abs(dstrike(2)-dstrike(1));
%     diplen=abs(ddip(2)-ddip(1));
    datfile(6)=dstrike(1);
    datfile(7)=dstrike(2);
    datfile(8)=ddip(1);
    datfile(9)=ddip(2);
    
    TIDES='TIDES/p085p020';
    tidals=load(TIDES);
    year=str2double(IDENTIF4(4:7));
    dayan=str2double(IDENTIF4(9:11));
    day2=datenum(year,1,1);
    day1=datenum(2000,1,1);
    day=(day2-day1)+dayan-1;
    tidaltimes=tim/86400
    periods=[0.517525  0.500000  0.527429  1.075804  0.997271];
%     tau12p4=tidals(1)*cos(2*pi*(day+tidaltimes)/periods(1)-tidals(2)*pi/180.); %Shear, 12.4 hrs.
    tautot=tidals(1)*cos(2*pi*(day+tidaltimes)/periods(1)-tidals(2)*pi/180.) ...
           +tidals(11)*cos(2*pi*(day+tidaltimes)/periods(2)-tidals(12)*pi/180.) ...
           +tidals(21)*cos(2*pi*(day+tidaltimes)/periods(3)-tidals(22)*pi/180.) ...
           +tidals(31)*cos(2*pi*(day+tidaltimes)/periods(4)-tidals(32)*pi/180.) ...
           +tidals(41)*cos(2*pi*(day+tidaltimes)/periods(5)-tidals(42)*pi/180.); %Shear, all 5, I hope.
%     CFStot=tidals(7)*cos(2*pi*(day+tidaltimes)/periods(1)-tidals(8)*pi/180.) ...
%            +tidals(17)*cos(2*pi*(day+tidaltimes)/periods(2)-tidals(18)*pi/180.) ...
%            +tidals(27)*cos(2*pi*(day+tidaltimes)/periods(3)-tidals(28)*pi/180.) ...
%            +tidals(37)*cos(2*pi*(day+tidaltimes)/periods(4)-tidals(38)*pi/180.) ...
%            +tidals(47)*cos(2*pi*(day+tidaltimes)/periods(5)-tidals(48)*pi/180.); %Shear, all 5, I hope.
    datfile(10)=tautot;

    fprintf(fidprop,'%5i %8.1f %9.3f %3i %11.3e %6.2f %6.2f %6.2f %6.2f %11.3e \n',datfile);

    set(h,'PaperPosition',[0.25 0.25 8 10.5])
    print(h,'-depsc',['PROP/',IDENTIF4,'_',int2str(k),'.eps'])
      
end

