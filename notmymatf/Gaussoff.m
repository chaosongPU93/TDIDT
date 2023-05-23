close all
clear all
scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);

x = (-4:0.1:4);
lenx=length(x);

h=figure('Units','inches','Position',[1 1 8.5 11]);
subplot(3,2,1)
sigma=0.5;
off=0.5*sigma; 
distr1=gaussmf(x, [sigma,off]);
distr2=gaussmf(x, [sigma,-off]);
plot(x,distr1)
hold on
plot(x,distr2)
distrsum=distr1+distr2;
plot(x,distrsum)
ilast=0;
for i=1:lenx
    num=round(100*distrsum(i));
    if num > 0
        weights(ilast+1:ilast+num)=x(i);
        ilast=ilast+num;
    end
end
pd=fitdist(weights','Normal')
y=pdf(pd,x);
yscale=distrsum(ceil(lenx/2))/y(ceil(lenx/2));
plot(x,yscale*y,'c','linewidth',2)
text(1.,1.3,['\sigma = ','0.560'])
xlim([-2 2])

subplot(3,2,3)
sigma=0.5;
off=0.75*sigma; 
distr1=gaussmf(x, [sigma,off]);
distr2=gaussmf(x, [sigma,-off]);
plot(x,distr1)
hold on
plot(x,distr2)
distrsum=distr1+distr2;
plot(x,distrsum)
ilast=0;
for i=1:lenx
    num=round(100*distrsum(i));
    if num > 0
        weights(ilast+1:ilast+num)=x(i);
        ilast=ilast+num;
    end
end
pd=fitdist(weights','Normal')
y=pdf(pd,x);
yscale=distrsum(ceil(lenx/2))/y(ceil(lenx/2));
plot(x,yscale*y,'c','linewidth',2)
text(1.,1.3,['\sigma = ','0.628'])
xlim([-2 2])

subplot(3,2,5)
sigma=0.5;
off=1.*sigma; 
distr1=gaussmf(x, [sigma,off]);
distr2=gaussmf(x, [sigma,-off]);
plot(x,distr1)
hold on
plot(x,distr2)
distrsum=distr1+distr2;
plot(x,distrsum)
ilast=0;
for i=1:lenx
    num=round(100*distrsum(i));
    if num > 0
        weights(ilast+1:ilast+num)=x(i);
        ilast=ilast+num;
    end
end
pd=fitdist(weights','Normal')
y=pdf(pd,x);
yscale=distrsum(ceil(lenx/2))/y(ceil(lenx/2));
plot(x,yscale*y,'c','linewidth',2)
text(1.,1.3,['\sigma = ','0.708'])
xlim([-2 2])
set(h,'PaperPosition',[0.25 0.25 8 10.5])
print('-depsc','outGauss.eps')

off=1; 
subplot(3,2,2)
sigma=0.94/0.675;
distr1=gaussmf(x, [sigma,off]);
distr2=gaussmf(x, [sigma,-off]);
plot(x,distr1)
hold on
plot(x,distr2)
distrsum=distr1+distr2;
plot(x,distrsum)
ilast=0;
for i=1:lenx
    num=round(100*distrsum(i));
    if num > 0
        weights(ilast+1:ilast+num)=x(i);
        ilast=ilast+num;
    end
end
pd=fitdist(weights','Normal')
y=pdf(pd,x);
yscale=distrsum(ceil(lenx/2))/y(ceil(lenx/2));
plot(x,yscale*y,'c','linewidth',2)
text(1.,1.3,['\sigma = ','1.45'])
xlim([-4 4])

subplot(3,2,4)
sigma=0.94;
distr1=gaussmf(x, [sigma,off]);
distr2=gaussmf(x, [sigma,-off]);
plot(x,distr1)
hold on
plot(x,distr2)
distrsum=distr1+distr2;
plot(x,distrsum)
ilast=0;
for i=1:lenx
    num=round(100*distrsum(i));
    if num > 0
        weights(ilast+1:ilast+num)=x(i);
        ilast=ilast+num;
    end
end
pd=fitdist(weights','Normal')
y=pdf(pd,x);
yscale=distrsum(ceil(lenx/2))/y(ceil(lenx/2));
plot(x,yscale*y,'c','linewidth',2)
text(1.,1.3,['\sigma = ','1.41'])
xlim([-4 4])

% subplot(3,2,5)
% sigma=0.5;
% off=1.*sigma; 
% distr1=gaussmf(x, [sigma,off]);
% distr2=gaussmf(x, [sigma,-off]);
% plot(x,distr1)
% hold on
% plot(x,distr2)
% distrsum=distr1+distr2;
% plot(x,distrsum)
% ilast=0;
% for i=1:lenx
%     num=round(100*distrsum(i));
%     if num > 0
%         weights(ilast+1:ilast+num)=x(i);
%         ilast=ilast+num;
%     end
% end
% pd=fitdist(weights','Normal')
% y=pdf(pd,x);
% yscale=distrsum(ceil(lenx/2))/y(ceil(lenx/2));
% plot(x,yscale*y,'c','linewidth',2)
% text(1.,1.3,['\sigma = ','0.708'])
% xlim([-2 2])
% set(h,'PaperPosition',[0.25 0.25 8 10.5])
% print('-depsc','outGauss.eps')
% 
