%Sits at one spot (rotations; delays) and looks for the cross-correlation.
% "extract" looks only at that active strip running down the middle
clear all
close all
format short e
rotline=(0./180.)*pi;

IDENTIF4='map2005.259.36.75.80.365.385_2-6-4s';
IDENTIF128='map2005.259.36.75.80.365.385_2-6-128s';
ALL128='all.36.75.128s';

TIDES='TIDES/p036p075';
tidals=load(TIDES);
year=str2double(IDENTIF4(4:7))
dayan=str2double(IDENTIF4(9:11))
day2=datenum(year,1,1);
day1=datenum(2000,1,1);
day=(day2-day1)+dayan-1;
% an=2000;
% day=0;
% while an < year
%     if(mod(an,4)==0)
%         day=day+366;
%     else
%         day=day+365;
%     end
%     an=an+1;
% end
% day=day+dayan-1
tidaltimes=0:0.01:2;
periods=[0.517525  0.500000  0.527429  1.075804  0.997271];
tau12p4=tidals(1)*cos(2*pi*(day+tidaltimes)/periods(1)-tidals(2)*pi/180.); %Shear, 12.4 hrs.
%tau12p4=tidals(7)*cos(2*pi*(day+tidaltimes)/periods(1)-tidals(8)*pi/180.); %Coulomb, 12.4 hrs, I think.
tautot=tidals(1)*cos(2*pi*(day+tidaltimes)/periods(1)-tidals(2)*pi/180.) ...
       +tidals(11)*cos(2*pi*(day+tidaltimes)/periods(2)-tidals(12)*pi/180.) ...
       +tidals(21)*cos(2*pi*(day+tidaltimes)/periods(3)-tidals(22)*pi/180.) ...
       +tidals(31)*cos(2*pi*(day+tidaltimes)/periods(4)-tidals(32)*pi/180.) ...
       +tidals(41)*cos(2*pi*(day+tidaltimes)/periods(5)-tidals(42)*pi/180.); %Shear, all 5, I hope.
CFStot=tidals(7)*cos(2*pi*(day+tidaltimes)/periods(1)-tidals(8)*pi/180.) ...
       +tidals(17)*cos(2*pi*(day+tidaltimes)/periods(2)-tidals(18)*pi/180.) ...
       +tidals(27)*cos(2*pi*(day+tidaltimes)/periods(3)-tidals(28)*pi/180.) ...
       +tidals(37)*cos(2*pi*(day+tidaltimes)/periods(4)-tidals(38)*pi/180.) ...
       +tidals(47)*cos(2*pi*(day+tidaltimes)/periods(5)-tidals(48)*pi/180.); %Shear, all 5, I hope.
% maxtau=max([tau12p4 tautot CFStot]);
% mintau=min([tau12p4 tautot CFStot]);
maxtau=max(tautot);
mintau=min(tautot);
maxtau=max(maxtau,-mintau);

orig4=load(IDENTIF4);
orig128=load(IDENTIF128);
origAll=load(ALL128);
% if OLDER==1
%     origolder=load(OLDER128);
% end
% orig4b=load(IDENTIF4b);
% orig128b=load(IDENTIF128b);
% origAllb=load(ALL128b);
% if OLDERb==1
%     origolderb=load(OLDER128b);
% end

locs4=orig4(:,2)+1i*orig4(:,3);
locs4=locs4*exp(1i*rotline);
orig4(:,2)=real(locs4);
orig4(:,3)=imag(locs4);
locs128=orig128(:,2)+1i*orig128(:,3);
locs128=locs128*exp(1i*rotline);
locsAll=origAll(:,1)+1i*origAll(:,2);
locsAll=locsAll*exp(1i*rotline);
origAll(:,1)=real(locsAll);
origAll(:,2)=imag(locsAll);

ntot=length(orig4);
nin=0;
for n = 1:ntot
    if(abs(orig4(n,3))<=2.)
        nin=nin+1;
        orig4in(nin,:)=orig4(n,:);
    end
end

%plot(orig4(:,2),orig4(:,3),'ko','MarkerSize',2)
scatter3(orig4in(:,2),orig4in(:,3),86400-orig4in(:,1),25,orig4in(:,1),'filled')
axis equal
axis([-16 16 -10 10])

scrsz=get(0,'ScreenSize');
clo=-5;
chi=1;
dotsz=8;

h=figure('Position',[scrsz(3)/10 1 8*scrsz(3)/10 scrsz(4)]);
colormap(jet) 
hold on
%     plot(real(locsAll),imag(locsAll),'ko','MarkerSize',2,'linewidth',0.2)
%     plot(real(locs128(1:istart128)),imag(locs128(1:istart128)),'ro','MarkerSize',2,'linewidth',0.2)
%     if OLDER==1
%         plot(real(olderlocs),imag(olderlocs),'ro','MarkerSize',2,'linewidth',0.2)
%     end
%     scatter(real(locs128(istart128:iend128)),imag(locs128(istart128:iend128)),90,orig128(istart128:iend128,1)...
%         ,'linewidth',1)
%     % First on top
%     scatter3(real(locs4(istart4:iend4)),imag(locs4(istart4:iend4)),86400-orig4(istart4:iend4,1),25,orig4(istart4:iend4,1),'filled')
%     % Last on top
%     %scatter3(real(locs4(istart4:iend4)),imag(locs4(istart4:iend4)),orig4(istart4:iend4,1),40,orig4(istart4:iend4,1),'filled')
scatter3(orig4in(:,1),orig4in(:,2),86400-orig4in(:,1),dotsz,log(orig4in(:,7)),'filled')
plot(tidaltimes*86400,2*tautot/maxtau-10,'k--')
axis([43200 86400 -14 14])
%xlabel('time (s)')
ylabel('km along strike')
caxis([clo chi])
grid minor
%set(axes_handle,'YGrid','on')
colorbar
title(IDENTIF4)
box on
% subplot(4,1,2,'align')
% colormap(jet) 
% hold on
% scatter3(orig4(:,1)/3600,imag(locs4),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
% plot(tidaltimes*86400/3600,tautot*(4/maxtau),'k--')
% plot(tidaltimes*86400/3600,CFStot*(4/maxtau),'k:')
% plot(tidaltimes*86400/3600,tau12p4*(4/maxtau),'r--')
% axis([0 24 -6 4])
% set(gca,'xtick',0:3:24);
% %xlabel('time (s)')
% ylabel('km along dip')
% caxis([clo chi])
% colorbar
% %title(IDENTIF4)
% box on
% 
% subplot(4,1,3,'align')
% colormap(jet) 
% hold on
% scatter3(orig4b(:,1),real(locs4b),86400-orig4b(:,1),dotsz,log(orig4b(:,7)),'filled')
% axis([0 86400 -5 5])
% %xlabel('time (s)')
% ylabel('km along strike')
% caxis([clo chi])
% colorbar
% title(IDENTIF4b)
% box on
% subplot(4,1,4,'align')
% colormap(jet) 
% hold on
% scatter3(orig4b(:,1)/3600,imag(locs4b),86400-orig4b(:,1),dotsz,log(orig4b(:,7)),'filled')
% plot((tidaltimes-1)*24,tautot*(4/maxtau),'k--')
% plot((tidaltimes-1)*24,CFStot*(4/maxtau),'k:')
% plot((tidaltimes-1)*24,tau12p4*(4/maxtau),'r--')
% axis([0 24 -6 4])
% set(gca,'xtick',0:3:24);
% xlabel('time (hrs)')
% ylabel('km along dip')
% caxis([clo chi])
% colorbar
% %title(IDENTIF4b)
% box on
% 
set(h,'PaperPosition',[0.25 0.25 8 10.5])
orient landscape
print(gcf,'-depsc',['AmpvTimeEx',IDENTIF4,'.eps'])

