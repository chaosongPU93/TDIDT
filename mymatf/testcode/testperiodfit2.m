clc
close all
clear


Fs = 40;
wlensec = 3.2;
wlen = wlensec*Fs;
t = 0:1/Fs:wlensec;
t = reshape(t,[],1);
f = 2;
T = 1/f;
amp = 1;

x = amp*sin(2*pi*f.*t);
x = reshape(x,[],1);

figure
plot(t,x,'k'); hold on

perctap = 0.2;
%T*Fs == Fs/f   % samples of one full cycle
xtap = x.*tukeywin(length(x),perctap);
plot(t,xtap,'b');

lo = 1;
hi = 3;
xbp = Bandpass(x,Fs,lo,hi,2,2,'butter');
plot(t,xbp,'r');

xtapbp = Bandpass(xtap,Fs,lo,hi,2,2,'butter');
plot(t,xbp,'c');    % basically cannot distinguish between the taper--bandpass and direct bandpass 

%%
%although the bandpass doesnot contain the f of data, you still get the wave with same f but with
%greatly reduced amplitude, which is expected due to the feature of the filter, the freq outside
%corners are suppressed but not 0. Meanwhile, you won't see any waves having the frequency of your
%requested range just because it does not have it.
figure
plot(t,xtap,'k'); hold on

lo = 1;
hi = 2;
xbp1 = Bandpass(xtap,Fs,lo,hi,2,2,'butter');
plot(t,xbp1,'r');

lo = 0.1;
hi = 1.5;
xbp2 = Bandpass(xtap,Fs,lo,hi,2,2,'butter');
plot(t,xbp2,'b');


%%
%make little difference when using the untapered data, but the resulting wave with same freq is
%smaller in amp when using a passband does not contain the target freq
figure
plot(t,x,'k'); hold on

lo = 1;
hi = 2;
xbp1 = Bandpass(x,Fs,lo,hi,2,2,'butter');
plot(t,xbp1,'r');

lo = 0.1;
hi = 1;
xbp2 = Bandpass(x,Fs,lo,hi,2,2,'butter');
plot(t,xbp2,'b');


%% periodfit
figure
plot(t,xtap,'k'); hold on

% dT can not be too small, otherwis the 'inv' would be unstable
% Matrix is close to singular or badly scaled
dT = 0.1;
perdx = 0.1: dT: 1.5;
% perdx = 10;
% comp = 'parallel';
comp = 'series';
yofx = xtap;
[vard,resd,pred,contr,perdb]=periodfit(yofx,perdx,t,comp);
% the inverted best prediction has to periodic sine waves, but amp can be smaller if the input
% signal is already tapered, otherwise the amplitude would be exactly the same 
plot(t,pred,'r');

%if you want to invert for some period way outside the target frequency, surely you will get
%something of your requested period but the amplitude would be very small, again, this is due to
%aliasing from finite length of record
perdx1 = 2;
[~,~,pred1,~,~]=periodfit(yofx,perdx1,t,comp);
plot(t,pred1,'b');

%this plot gives the variance reduce of the fitting residual, so the larger the better. it seems
%from this plot that because of the finite length of the record, whether it is tapered or not, there
%will be aliasing in frequency (period) domain, i.e., leaking of energy 
figure
plot(perdx,vard);









