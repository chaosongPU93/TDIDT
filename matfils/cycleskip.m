close all
clear all
scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite;
h=figure('Position',[wid/3 1 1.5*wid hite]); %center
seism=zeros(1,1000);
x=(0:49);
lenx=length(x);
seism(500:549)=-sin(2*pi*x/lenx);
seism(550:599)=-0.3*sin(2*pi*x/lenx);
seism(600:649)=-0.2*sin(2*pi*x/lenx);
composite=seism;
subplot(3,1,1,'align')
plot(seism)
box on
for n=1:99
    composite(150:850)=composite(150:850)+exp(-n/20)*seism(150+n:850+n);
    composite(150:850)=composite(150:850)+exp(-n/20)*seism(150-n:850-n);
end
subplot(3,1,2,'align')
hold on
plot(composite)
disp=cumsum(composite);
disp=min(composite)*disp/min(disp);
plot(disp,'r')
box on
composite=seism;
for n=1:49
    composite(150:850)=composite(150:850)+exp(-n/10)*seism(150+n:850+n);
    composite(150:850)=composite(150:850)+exp(-n/10)*seism(150-n:850-n);
end
for n=1:49
    composite(150:850)=composite(150:850)+0.1*exp(-n/10)*seism(200+n:900+n);
    composite(150:850)=composite(150:850)+0.1*exp(-n/10)*seism(200-n:900-n);
end
for n=1:49
    composite(150:850)=composite(150:850)+0.02*exp(-n/10)*seism(250+n:950+n);
    composite(150:850)=composite(150:850)+0.02*exp(-n/10)*seism(250-n:950-n);
end
%composite(100:1900)=composite(100:1900)+1*seism(200:2000);
subplot(3,1,3,'align')
hold on
plot(composite)
disp=cumsum(composite);
disp=min(composite)*disp/min(disp);
plot(disp,'r')
box on
set(h,'PaperPosition',[0.25 0.25 8 10.5])
print(h,'-depsc','misalign.eps')