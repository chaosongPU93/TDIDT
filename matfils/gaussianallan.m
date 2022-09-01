format short e
clear all
close all

dt=0.0025;                            %sampling time
Fs=1/dt;                            %sampling frequency
t=-10:dt:10;
L=length(t);
    
%figure
for np=1:1 %5
    Tp=0.1*2^(np-1)                        %predominant period
    
    %disp=(2^(np-1))^2*gaussmf(t,[Tp/4 0]);     %Gaussian

%     ileft=round((-Tp/2-t(1))/dt+1)                          %Boxcar
%     iright=round(L-(t(L)-Tp/2)/dt)                          %Boxcar
%     disp(1:ileft-1)=0.;                                     %Boxcar
%     disp(ileft:iright-1)=0.5*(2^(np-1))^2;                  %Boxcar
%     disp(iright:L)=0.;                                      %Boxcar

    ileft=round((-Tp/2-t(1))/dt+1)                                   %triangle
    iright=round(L-(t(L)-Tp/2)/dt)                                   %triangle
    disp(1:ileft-1)=0.;                                              %triangle
    disp(ileft:(L+1)/2-1)=(2^(np-1))^2*(0:dt/(Tp/2):1-dt/(Tp/2));    %triangle
    disp((L+1)/2:iright-1)=(2^(np-1))^2*(1:-dt/(Tp/2):dt/(Tp/2));    %triangle
    disp(iright:L)=0.;                                               %triangle

    co=25.;
    npol=2;
    npas=1;
    [dispf]=lowpass(disp,Fs,co,npol,npas,'butter','linear');
    figure
    subplot(3,1,1)
    plot(t,dispf,'k')
    xlim([-2.5 2.5]);

    vel=diff(dispf);
    tp=t(1:length(t)-1)+0.5*dt;
    subplot(3,1,2)
    plot(tp,vel)
    xlim([-2.5 2.5]);
    hold on
    lo=1.5;
    hi=6;
    xs=[1/lo 1/hi];
    ys=[0 0];
    plot(xs-2.5,ys,'r','linewidth',6)
    
    npo=2; npa=1;
    [velf]=bandpass(vel,Fs,lo,hi,npo,npa,'butter');
    subplot(3,1,3)
    plot(tp,velf,'r')
    xlim([-2.5 2.5]);

%     NFFT = 2^nextpow2(L);
%     f = Fs/2*linspace(0,1,NFFT/2+1);
%     ffdisp=fft(dispf,NFFT)/L;
%     figure
%     loglog(f,2*abs(ffdisp(1:NFFT/2+1))) 
%     xLimits = [1e-2 1e2];
%     yLimits = [1e-6 1e1];
%     logScale = diff(yLimits)/diff(xLimits);  %# Scale between the x and y ranges
%     powerScale = diff(log10(yLimits))/...    %# Scale between the x and y powers
%                   diff(log10(xLimits));
%     set(gca,'Xlim',xLimits,'YLim',yLimits,...              %# Set the limits and the
%            'DataAspectRatio',[1 logScale/powerScale 1]);  %#   data aspect ratio
%     hold on
%     ffvel=fft(velf,NFFT)/L;
%     loglog(f,2*abs(ffvel(1:NFFT/2+1)),'r') 
    
end

