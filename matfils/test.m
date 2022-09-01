%Reads in seismograms; filters and rotates them, etc.  Can plot coherence 
%on optimal component, or amplitude spectrum on both optimal and orthogonal components,
%and their ratio, for SILB and SSIB on top of one another for specified 
%windows during one day.  
%Avrages many 1024-sample windows (at 100 sps)
format short e
clear all
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

POLSTA='SILB'; 
nsta=size(POLSTA,1);
POLROTS=[0 27.5 210 20 -150];  %SILB.  Last entry is P-wave delay (for SILB, -200+20)
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.; %2 and 3 are now incidence angle and back-azimuth
POLROTS(:,1)=round(POLROTS(:,1)*(100/40)); %Polaris stations at 100 sps; table has split time at 40 sps.
load BOSTOCK/NEW/WFPV/wfpv_002.mat;
stas='SILB';
lensta=size(stas,1);
for ista=1:lensta
    [LIA,idx]=ismember(stas(ista,:),ordlst,'rows');
    STAE=ecomp(idx,:);
    STAN=ncomp(idx,:);
    STAZ=zcomp(idx,:);
    hi=48;
    lo=0.1;
    sps=100;
    tracelen=length(STAE);
    npo=2;
    npa=2;
        [STAEf]=bandpass(STAE,sps,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
        [STANf]=bandpass(STAN,sps,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
        [STAZf]=bandpass(STAZ,sps,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
    clear STAE
    clear STAN
    clear STAZ
%         STAEfd = resample(STAEf,2,5);
%         STANfd = resample(STANf,2,5);
%         STAZfd = resample(STAZf,2,5);
%         sps=100;
%         %Need to rotate E/N by 90 degrees to be consistent with this backaz and incidence angle.  But no, it just gets rotated back again.
%         STA=STAEf+1i*STANf;
%         STAfastslow=STA*exp(-1i*pi/2);
%         STAslow=real(STAfastslow); %(N)
%         STAfast=imag(STAfastslow); %(E)
%         STAsplitcorrected=(STAslow+1i*STAfast)*exp(1i*pi/2);
%         STAscrot=STAsplitcorrected*exp(-1i*POLROTS(ista,3));
    m1=[cos(POLROTS(ista,2)) -sin(POLROTS(ista,2))*sin(POLROTS(ista,3)) -sin(POLROTS(ista,2))*cos(POLROTS(ista,3))];
    m2=[sin(POLROTS(ista,2)) cos(POLROTS(ista,2))*sin(POLROTS(ista,3)) cos(POLROTS(ista,2))*cos(POLROTS(ista,3))];
    m3=[0 -cos(POLROTS(ista,3)) sin(POLROTS(ista,3))];
    L=m1(1)*STAZf+m1(2)*STAEf+m1(3)*STANf;
    Q=m2(1)*STAZf+m2(2)*STAEf+m2(3)*STANf;
    T=m3(1)*STAZf+m3(2)*STAEf+m3(3)*STANf;
    clear STAEf
    clear STANf
    clear STAZf
    len=tracelen'

    STAsoff=round(POLROTS(ista,4)*(100/40))
    if STAsoff > -1
        Q(1:len-STAsoff)=Q(STAsoff+1:len);
        Q(len-STAsoff+1:len)=0;
        T(1:len-STAsoff)=T(STAsoff+1:len);
        T(len-STAsoff+1:len)=0;
    else
        Q(-STAsoff+1:len)=Q(1:len+STAsoff);
        Q(1:-STAsoff)=0;
        T(-STAsoff+1:len)=T(1:len+STAsoff);
        T(1:-STAsoff)=0;
    end
    STApoff=round(POLROTS(ista,5))
    L(-STApoff+1:len)=L(1:len+STApoff); %For P, just assume offset < 0
    L(1:-STAsoff)=0;
    STAopt=Q;
    STAort=T;
    STAvrt=L;  
    clear Q
    clear T
    clear L
    STAoptfilt=bandpass(STAopt,sps,2.5,6.5,npo,npa,'butter');
    STAoptenv=abs(hilbert(STAoptfilt));
    figure
	plot(STAopt)
	hold on
	plot(STAort,'g')
	plot(STAvrt,'r')
end

