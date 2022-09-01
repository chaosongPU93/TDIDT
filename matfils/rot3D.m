%This version takes home-made template stack, rather than Bostock's.
%For a single station.
%Rotation matrix from IRIS website.  inc is incidence angle measured from
%vertical.  baz is back-azimuth measured clockwise from N.
%Searches for inc and baz that (1) maximizes P; (2) maximizes P/S.  Go for
%the latter.
clear all
close all
datfile='SSIB_002_0.5-16Hz_2pass';
dat=load(datfile);
sta=datfile(1:4);
sps=100;
%rots=[0 90 200];  %PGC , Yajun's "0" changed to 90. (32)
scrsz=get(0,'ScreenSize');
zzero=[0 0;
       60*sps 0];
   
STAe=dat(:,1);
STAn=dat(:,2);
STAz=dat(:,3);
% STA=STAe+1i*STAn;
%     %Rotate & split-correct
%     STAfastslow=STA*exp(-1i*rots(ista,2));
%     STAslow=real(STAfastslow);
%     STAfast=imag(STAfastslow);
%     len=length(STA);
%     offmax=10;
%     STAslow(offmax:len-offmax)=STAslow(offmax+rots(ista,1):len-offmax+rots(ista,1));
%     STAsplitcorrected=(STAslow+1i*STAfast)*exp(1i*rots(ista,2));
%     STAe=real(STAsplitcorrected);
%     STAn=imag(STAsplitcorrected);
h=figure('Position',[scrsz(3)/3.1 1 scrsz(3)/3.1 scrsz(4)]); %center
subplot(6,1,1,'align')
hold on
plot(STAe,'b')
plot(STAn,'r')
plot(-STAz,'k')
mx=max([max(abs(STAe)) max(abs(STAn)) max(abs(STAz))]);
title([sta,' E, N, & -Z'])
xlim([15*sps 30*sps])
ylim([-0.4*mx 0.4*mx])
[xp,~]=ginput(2);
[xs,~]=ginput(2);
box on

subplot(6,1,2,'align')
hold on
plot(STAe,'b')
plot(STAn,'r')
plot(-STAz,'k')
title('E, N, & -Z  (P)')
xlim([xp(1) xp(2)])
box on

subplot(6,1,3,'align')
hold on
plot(STAe,'b')
plot(STAn,'r')
plot(-STAz,'k')
title('E, N, & -Z  (S)')
xlim([xs(1) xs(2)])
box on

ndiv=36;
inc=(0:pi/(4*ndiv):pi/2);
baz=(0:pi/ndiv:2*pi);
%baz=rots(ista,3);
isP=round(xp(1));
ieP=round(xp(2));
isS=round(xs(1));
ieS=round(xs(2));
ampPmax=0;
ratmax=0;
for i=1:2*ndiv+1
    for j=1:2*ndiv+1
        m1=[cos(inc(i)) -sin(inc(i))*sin(baz(j)) -sin(inc(i))*cos(baz(j))];
        L=m1(1)*STAz+m1(2)*STAe+m1(3)*STAn;
%         sqrP=sqrt(L(isP:ieP).*L(isP:ieP));
%         sqrS=sqrt(L(isS:ieS).*L(isS:ieS));
        sqrP=L(isP:ieP).*L(isP:ieP);
        sqrS=L(isS:ieS).*L(isS:ieS);
        ampP(j,i)=sum(sqrP);
        ampS(j,i)=sum(sqrS);
        if ampP(j,i)>ampPmax
            ampPmax=ampP(j,i);
            incPmax=inc(i);
            bazPmax=baz(j);
        end
        if ampP(j,i)/ampS(j,i)>ratmax
            ratmax=ampP(j,i)/ampS(j,i);
            incRmax=inc(i);
            bazRmax=baz(j);
        end
    end
end
incPm=incPmax*180/pi
bazPm=bazPmax*180/pi
incRm=incRmax*180/pi
bazRm=bazRmax*180/pi
subplot(6,1,4,'align')
m1=[cos(incPmax) -sin(incPmax)*sin(bazPmax) -sin(incPmax)*cos(bazPmax)];
m2=[sin(incPmax) cos(incPmax)*sin(bazPmax) cos(incPmax)*cos(bazPmax)];
m3=[0 -cos(bazPmax) sin(bazPmax)];
L=m1(1)*STAz+m1(2)*STAe+m1(3)*STAn;
Q=m2(1)*STAz+m2(2)*STAe+m2(3)*STAn;
T=m3(1)*STAz+m3(2)*STAe+m3(3)*STAn;
hold on
plot(Q,'g')
plot(T,'b')
plot(L,'r')
title('L(r), T(b), & Q(g)')
xlim([15*sps 30*sps])
ylim([-0.4*mx 0.4*mx])
box on

subplot(6,1,5,'align')
hold on
plot(STAz,'k')
plot(L,'g')
xlim([15*sps 30*sps])
title('Z(k) & L(g; maximize P)')
mx=max(abs(STAz));
ylim([-mx mx])
text(15.6*sps, 0.7*mx, ['baz: ',num2str(bazPm)],'fontsize',8);
text(15.6*sps, 0.4*mx, ['inc: ',num2str(incPm)],'fontsize',8);
box on

subplot(6,1,6,'align')
m1=[cos(incRmax) -sin(incRmax)*sin(bazRmax) -sin(incRmax)*cos(bazRmax)];
m2=[sin(incRmax) cos(incRmax)*sin(bazRmax) cos(incRmax)*cos(bazRmax)];
m3=[0 -cos(bazRmax) sin(bazRmax)];
L=m1(1)*STAz+m1(2)*STAe+m1(3)*STAn;
Q=m2(1)*STAz+m2(2)*STAe+m2(3)*STAn;
T=m3(1)*STAz+m3(2)*STAe+m3(3)*STAn;
hold on
plot(STAz,'k')
plot(L,'r')
title('Z(k) & L(r; maximize P/S)')
xlim([15*sps 30*sps])
ylim([-mx mx])
text(15.6*sps, 0.7*mx, ['baz: ',num2str(bazRm)],'fontsize',8);
text(15.6*sps, 0.4*mx, ['inc: ',num2str(incRm)],'fontsize',8);
box on

