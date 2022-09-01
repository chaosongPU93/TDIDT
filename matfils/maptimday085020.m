%Sits at one spot (rotations; delays) and looks for the cross-correlation.
clear all
close all
format short e

rotang=(22./180.)*pi;

% % IDENTIF4='map2003.061.85.20.80.115.55_2-8-4s';
% % IDENTIF128='map2003.061.85.20.80.115.55_2-6-128s';
% % ALL128='all.85.20.128s.fort.13.sor';
% % OLDER=0;
%
% IDENTIF4='map2003.062.85.20.80.115.55_2-8-4s';
% IDENTIF128='map2003.062.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER128='map2003.061.85.20.80.115.55_2-6-128s';
% OLDER=1;
% 
% IDENTIF4='map2003.063.85.20.80.115.55_2-8-4s';
% IDENTIF128='map2003.063.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER128='map2003.061-062.85.20.80.115.55_2-6-128s';
% OLDER=1;

% IDENTIF4='map2004.196.85.20.80.115.55_2-6-4s';
% IDENTIF128='map2004.196.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER=0;
% 
% IDENTIF4='map2004.197.85.20.80.115.55_2-6-4s';
% IDENTIF128='map2004.197.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER128='map2004.196.85.20.80.115.55_2-6-128s';
% OLDER=1;
% %
% IDENTIF4='map2004.198.85.20.80.115.55_2-6-4s';
% IDENTIF128='map2004.198.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER128='map2004.196-197.85.20.80.115.55_2-6-128s';
% OLDER=1;
% 
% IDENTIF4='map2004.199.85.20.80.115.55_2-6-4s';
% IDENTIF128='map2004.197.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER128='map2004.196.85.20.80.115.55_2-6-128s';
% OLDER=1;

% % % IDENTIF4='map2005.253.85.20.80.115.55_2-6-4s';
% % % IDENTIF128='map2005.253.85.20.80.115.55_2-6-128s';
% % % ALL128='all.85.20.128s.fort.13.sor';
% % % OLDER=0;
%
IDENTIF4='map2005.254.85.20.80.115.55_2-6-4s';
IDENTIF128='map2005.254.85.20.80.115.55_2-6-128s';
ALL128='all.85.20.128s.fort.13.sor';
OLDER=0;
% 
% IDENTIF4='map2005.255.85.20.80.115.55_2-6-4s';
% IDENTIF128='map2005.255.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER128='map2005.254.85.20.80.115.55_2-6-128s';
% OLDER=1;
% %
% IDENTIF4='map2005.256.85.20.80.115.55_2-6-4s';
% IDENTIF128='map2005.256.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER128='map2005.255.85.20.80.115.55_2-6-128s';
% OLDER=1;
% 
% IDENTIF4='map2005.257.85.20.80.115.55_2-6-4s';
% IDENTIF128='map2005.256.85.20.80.115.55_2-6-128s';
% ALL128='all.85.20.128s.fort.13.sor';
% OLDER128='map2005.255.85.20.80.115.55_2-6-128s';
% OLDER=1;

% TIDES='TIDES/p085p020';
% tidals=load(TIDES);
% year=str2double(IDENTIF4(4:7));
% dayan=str2double(IDENTIF4(9:11));
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
% day=day+dayan-1;
% tidaltimes=0:0.01:1;
% periods=[0.517525  0.500000  0.527429  1.075804  0.997271];
% tau12p4=tidals(1)*cos(2*pi*(day+tidaltimes)/periods(1)-tidals(2)*pi/180.);
% maxtau=max(tau12p4);
% mintau=min(tau12p4);
% maxtau=max(maxtau,-mintau);

orig4=load(IDENTIF4);
orig128=load(IDENTIF128);
origAll=load(ALL128);
if OLDER==1
    origolder=load(OLDER128);
end

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

scrsz=get(0,'ScreenSize');
clo=-5;
chi=2;
dotsz=8;
tees=0:10800:86400;
extra=1000;

h=figure('Position',[scrsz(3)/10 1 8*scrsz(3)/10 scrsz(4)]);
subplot(4,1,1,'align')
colormap(jet) 
hold on
scatter3(orig4(:,1),real(locs4),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
axis([tees(1)-extra tees(2)+extra -5 5])
ylabel('km along strike')
caxis([clo chi])
colorbar
title(IDENTIF4)
box on
subplot(4,1,2,'align')
colormap(jet) 
hold on
scatter3(orig4(:,1),imag(locs4),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
axis([tees(1)-extra tees(2)+extra -6 4])
ylabel('km along dip')
caxis([clo chi])
colorbar
box on
subplot(4,1,3,'align')
colormap(jet) 
hold on
scatter3(orig4(:,1),real(locs4),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
axis([tees(2)-extra tees(3)+extra -5 5])
ylabel('km along strike')
caxis([clo chi])
colorbar
box on
subplot(4,1,4,'align')
colormap(jet) 
hold on
scatter3(orig4(:,1),imag(locs4),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
axis([tees(2)-extra tees(3)+extra -6 4])
xlabel('time (s)')
ylabel('km along dip')
caxis([clo chi])
colorbar
box on
set(h,'PaperPosition',[0.25 0.25 8 10.5])
orient landscape
print(gcf,'-depsc',['time',IDENTIF4,'a.eps'])

h=figure('Position',[scrsz(3)/10 1 8*scrsz(3)/10 scrsz(4)]);
subplot(4,1,1,'align')
colormap(jet) 
hold on
scatter3(orig4(:,1),real(locs4),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
axis([tees(3)-extra tees(4)+extra -5 5])
ylabel('km along strike')
caxis([clo chi])
colorbar
title(IDENTIF4)
box on
subplot(4,1,2,'align')
colormap(jet) 
hold on
scatter3(orig4(:,1),imag(locs4),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
axis([tees(3)-extra tees(4)+extra -6 4])
ylabel('km along dip')
caxis([clo chi])
colorbar
box on
subplot(4,1,3,'align')
colormap(jet) 
hold on
scatter3(orig4(:,1),real(locs4),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
axis([tees(4)-extra tees(5)+extra -5 5])
ylabel('km along strike')
caxis([clo chi])
colorbar
box on
subplot(4,1,4,'align')
colormap(jet) 
hold on
scatter3(orig4(:,1),imag(locs4),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
axis([tees(4)-extra tees(5)+extra -6 4])
xlabel('time (s)')
ylabel('km along dip')
caxis([clo chi])
colorbar
box on
set(h,'PaperPosition',[0.25 0.25 8 10.5])
orient landscape
print(gcf,'-depsc',['time',IDENTIF4,'b.eps'])

h=figure('Position',[scrsz(3)/10 1 8*scrsz(3)/10 scrsz(4)]);
subplot(4,1,1,'align')
colormap(jet) 
hold on
scatter3(orig4(:,1),real(locs4),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
axis([tees(5)-extra tees(6)+extra -5 5])
ylabel('km along strike')
caxis([clo chi])
colorbar
title(IDENTIF4)
box on
subplot(4,1,2,'align')
colormap(jet) 
hold on
scatter3(orig4(:,1),imag(locs4),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
axis([tees(5)-extra tees(6)+extra -6 4])
ylabel('km along dip')
caxis([clo chi])
colorbar
box on
subplot(4,1,3,'align')
colormap(jet) 
hold on
scatter3(orig4(:,1),real(locs4),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
axis([tees(6)-extra tees(7)+extra -5 5])
ylabel('km along strike')
caxis([clo chi])
colorbar
box on
subplot(4,1,4,'align')
colormap(jet) 
hold on
scatter3(orig4(:,1),imag(locs4),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
axis([tees(6)-extra tees(7)+extra -6 4])
xlabel('time (s)')
ylabel('km along dip')
caxis([clo chi])
colorbar
box on
set(h,'PaperPosition',[0.25 0.25 8 10.5])
orient landscape
print(gcf,'-depsc',['time',IDENTIF4,'c.eps'])

h=figure('Position',[scrsz(3)/10 1 8*scrsz(3)/10 scrsz(4)]);
subplot(4,1,1,'align')
colormap(jet) 
hold on
scatter3(orig4(:,1),real(locs4),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
axis([tees(7)-extra tees(8)+extra -5 5])
ylabel('km along strike')
caxis([clo chi])
colorbar
title(IDENTIF4)
box on
subplot(4,1,2,'align')
colormap(jet) 
hold on
scatter3(orig4(:,1),imag(locs4),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
axis([tees(7)-extra tees(8)+extra -6 4])
ylabel('km along dip')
caxis([clo chi])
colorbar
box on
subplot(4,1,3,'align')
colormap(jet) 
hold on
scatter3(orig4(:,1),real(locs4),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
axis([tees(8)-extra tees(9)+extra -5 5])
ylabel('km along strike')
caxis([clo chi])
colorbar
box on
subplot(4,1,4,'align')
colormap(jet) 
hold on
scatter3(orig4(:,1),imag(locs4),86400-orig4(:,1),dotsz,log(orig4(:,7)),'filled')
axis([tees(8)-extra tees(9)+extra -6 4])
xlabel('time (s)')
ylabel('km along dip')
caxis([clo chi])
colorbar
box on
set(h,'PaperPosition',[0.25 0.25 8 10.5])
orient landscape
print(gcf,'-depsc',['time',IDENTIF4,'d.eps'])

