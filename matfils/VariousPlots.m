%$Various 2-D tremor plots
format short e
clearvars
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');
scrsz=get(0,'ScreenSize');
wid=0.8*scrsz(3);
hite=0.8*scrsz(4);
scrat=wid/hite; 

%Mo=[0.001 0.01 0.1 1]*1e12;
LogL=1:0.01:3;
L=10.^LogL; %lfe source dimension, m
lenL=length(L);
Lo100=L/100;  %use 100 m as reference     
Lo100sq=Lo100.^2;
Lo100cu=Lo100.^3;
LogMo=9.0:0.01:12; 
Mo=10.^LogMo; %Characteristic LFE moment
Mo6e10=Mo/6e10; %use 6e10 as reference
lenMo=length(Mo);
LogSatAx=0:0.01:2;
SatAx=10.^LogSatAx; %saturation times in amplitude
lenSat=length(SatAx);
LogAmpoMin=0:0.01:1.2;
AmpoMin=10.^LogAmpoMin;
lenAmpoMin=length(AmpoMin);
T=[0.125 0.25]; %LFE duration in s
Cs=[1800 3600]; %shear wave speed in m/s
Vpropocs=[0.4 0.8]; %Vprop as a fraction of Cs
VpropBurst=5; %Burst prop speed in m/s，chao's LFE cat's Vprop can be 6 m/s
W=1000; %Burst pulse width in m, therefore, in the below expression for 'VSlipAve', you can just use 'W' and 'VpropBurst'  
VSlipAve=1/(1000/5); %Average slip speed in mm/s (multiply by [SlipTot/1 mm][Vprop/5 ms^-1]/[PulseWidth/1000 m])
SlipDur=W/VpropBurst;
SlipperBurst=[1 2]; %Burst slip in mm
VavBurst=SlipperBurst/SlipDur; %Average burst slip speed in mm/s
Lmax=max(T)*max(Cs)*max(Vpropocs);  %max lfe source dimension
Lmin=min(T)*min(Cs)*min(Vpropocs);  %min 
% Momax=1.e11;
Momax=5.e10;    %Bostock Mo_{min} of 5x10^{10}Nm       

D=zeros(lenMo,lenL); %slip
tau=zeros(lenMo,lenL);  %stress drop
Veq=zeros(lenMo,lenL);  %lfe slip speed
VeqoVdyn=zeros(lenMo,lenL); %lfe slip speed / Vdyn
LoLb=zeros(lenMo,lenL); %L / min. nucleation size 2L_b
DoDmin=zeros(lenMo,lenL);   %Slip / min. slip
SatN=1.2*ones(lenMo,lenL);  %Average saturation level (arrivals/period)
AreaMin=zeros(lenSat,lenL); %Minimum tremor source area (km^2)
SatB=ones(lenMo,lenAmpoMin);    %Instantaneous saturation
%Area=8e6*ones(lenMo,lenL);
%Nposs=zeros(lenMo,lenL);
for iL=1:lenL
    for iM=1:lenMo
        D(iM,iL)=0.42*Mo6e10(iM)/Lo100sq(iL); %in mm
        tau(iM,iL)=350*Mo6e10(iM)/Lo100cu(iL); %in kPa
        Veq(iM,iL)=1.7*Mo6e10(iM)/Lo100sq(iL); %in mm/s
        SatN(iM,iL)=SatN(iM,iL)/Mo6e10(iM);
        VeqoVdyn(iM,iL)=2.83*Mo6e10(iM)/Lo100sq(iL);
        DoDmin(iM,iL)=21.2*Mo6e10(iM)/Lo100sq(iL);
        LoLb(iM,iL)=0.012*L(iL);
    end
    for iS=1:lenSat
        AreaMin(iS,iL)=SatAx(iS)*(L(iL)/1000)^2;
    end
end
for iA=1:lenAmpoMin
    for iM=1:lenMo
        SatB(iM,iA)=(Mo6e10(iM)^(-2))*AmpoMin(iA)^2;
    end
end
Nrpts=ones(lenMo,lenL)./D;  %# repeats to amount to 1 mm total slip
FracActive=Nrpts*max(T)/SlipDur;    %active fraction of the pulse duration   
%Nposs=Area./L.^2;
%VeqoVavBurst=Veq/max(VavBurst);

h=figure('Units','inches','Position',[1 1 8.5 11]);

subplot(3,2,1)
contour(log10(L),log10(Mo),D,[1e-4 1e-3 1e-2 1e-1 1e-0 1e1 1e2],'-k','ShowText','on')
hold on
plot([log10(Lmax) log10(Lmax)],[log10(min(Mo)) log10(max(Mo))],'k--')
plot([log10(Lmin) log10(Lmin)],[log10(min(Mo)) log10(max(Mo))],'k--')
plot([log10(min(L)) log10(max(L))],[log10(Momax) log10(Momax)],'k--')
title('Slip (mm)','FontWeight','normal')
xticks=([1 1.5 2 2.5 3]);
xticklabels({'10 m',' ','100 m',' ','1 km'});
ax=gca;
ax.XTickLabelRotation = 0.;
xlabel('LFE source dimension')
yticks=([9.0 9.5 10 10.5 11 11.5 12]);
yticklabels({'10^{9}',' ','10^{10}',' ','10^{11}',' ','10^{12}'});
ylabel('Characteristic LFE moment (Nm)')

% keyboard
subplot(3,2,2)
%labels=['0.001' '0.01' '0.1' '1' '10' '100' '1e3' '1e4' '1e5'];
% contour(log10(L),log10(Mo),tau,[0.001 0.01 0.1 1 10 100 1.e3 1.e4 1.e5'],'-k','ShowText','on','LabelFormat','%.3g')
contour(log10(L),log10(Mo),tau,[0.001 0.01 0.1 1 10 100 1.e3 1.e4 1.e5'],'-k','ShowText','on')
hold on
plot([log10(Lmax) log10(Lmax)],[log10(min(Mo)) log10(max(Mo))],'k--')
plot([log10(Lmin) log10(Lmin)],[log10(min(Mo)) log10(max(Mo))],'k--')
plot([log10(min(L)) log10(max(L))],[log10(Momax) log10(Momax)],'k--')
title('''Stress drop'' (kPa) (isolated events)','FontWeight','normal')
xticks=([1 1.5 2 2.5 3]);
xticklabels({'10 m',' ','100 m',' ','1 km'});
ax=gca;
ax.XTickLabelRotation = 0.;
xlabel('LFE source dimension')
yticks=([9.0 9.5 10 10.5 11 11.5 12]);
yticklabels({'10^{9}',' ','10^{10}',' ','10^{11}',' ','10^{12}'});
%ylabel('Characteristic LFE moment (Nm)')

subplot(3,2,5)
contour(log10(L),log10(Mo),Nrpts,[1e4 1e3 1e2 1e1 1e0 1e-1 1e-2],'-k','ShowText','on')
hold on
plot([log10(Lmax) log10(Lmax)],[log10(min(Mo)) log10(max(Mo))],'k--')
plot([log10(Lmin) log10(Lmin)],[log10(min(Mo)) log10(max(Mo))],'k--')
plot([log10(min(L)) log10(max(L))],[log10(Momax) log10(Momax)],'k--')
title(["# repeats",...
    "(multiply by: [^{Slip_{tot}}/_{1 mm}][^{4-Hz slip}/_{Slip_{tot}}] )"],'FontWeight','normal')
xticks=([1 1.5 2 2.5 3]);
xticklabels({'10 m',' ','100 m',' ','1 km'});
ax=gca;
ax.XTickLabelRotation = 0.;
xlabel('LFE source dimension')
yticks=([9.0 9.5 10 10.5 11 11.5 12]);
yticklabels({'10^{9}',' ','10^{10}',' ','10^{11}',' ','10^{12}'});
ylabel('Characteristic LFE moment (Nm)')

subplot(3,2,4)
% contour(log10(L),log10(Mo),Veq/VSlipAve,[1e4 1e3 100 10 1e1 1e0 1e-1],'-k','ShowText','on','LabelFormat','%.3g')
contour(log10(L),log10(Mo),Veq/VSlipAve,[1e4 1e3 100 10 1e1 1e0 1e-1],'-k','ShowText','on')
%contour(log10(L),log10(Mo),FracActive,[0.001 0.01 0.1 1 10 100 1e-4 1e-5],'-k','ShowText','on')
hold on
plot([log10(Lmax) log10(Lmax)],[log10(min(Mo)) log10(max(Mo))],'k--')
plot([log10(Lmin) log10(Lmin)],[log10(min(Mo)) log10(max(Mo))],'k--')
plot([log10(min(L)) log10(max(L))],[log10(Momax) log10(Momax)],'k--')
title(["LFE slip speed / Ave. slip speed in pulse", ...
    "(divide by: [^{Slip_{tot}}/_{1 mm}][^{1 km}/_{W}][^{V_{prop}}/_{5 m/s}] )"],'FontWeight','normal')
    %"(duration = 2 km/5 ms^{-1} = 400 s)"],'FontWeight','normal')
xticks=([1 1.5 2 2.5 3]);
xticklabels({'10 m',' ','100 m',' ','1 km'});
ax=gca;
ax.XTickLabelRotation = 0.;
xlabel('LFE source dimension')
yticks=([9.0 9.5 10 10.5 11 11.5 12]);
yticklabels({'10^{9}',' ','10^{10}',' ','10^{11}',' ','10^{12}'});
%ylabel('Characteristic LFE moment (Nm)')

subplot(3,2,6)
contour(log10(L),log10(Mo),SatN,[1e1 1e0 1e-1 1e-2 1e-3 1e-4],'-k','ShowText','on')
hold on
plot([log10(Lmax) log10(Lmax)],[log10(min(Mo)) log10(max(Mo))],'k--')
plot([log10(Lmin) log10(Lmin)],[log10(min(Mo)) log10(max(Mo))],'k--')
plot([log10(min(L)) log10(max(L))],[log10(Momax) log10(Momax)],'k--')
title(["Average saturation level (arrivals/period)",...
       "(multiply by: [^{Slip_{tot}}/_{1 mm}][^{W_{perp}}/_{4 km}][^{A_{seism}}/_{A}][^{V_{prop}}/_{5 m/s}] )"],'FontWeight','normal')
xticks=([1 1.5 2 2.5 3]);
xticklabels({'10 m',' ','100 m',' ','1 km'});
ax=gca;
ax.XTickLabelRotation = 0.;
xlabel('LFE source dimension')
yticks=([9.0 9.5 10 10.5 11 11.5 12]);
yticklabels({'10^{9}',' ','10^{10}',' ','10^{11}',' ','10^{12}'});
%ylabel('Characteristic LFE moment (Nm)')

subplot(3,2,3)
contour(log10(L),log10(Mo),Veq,[1e2 1e1 1e0 1e-1 1e-2 1e-3 1e-4],'-k','ShowText','on')
hold on
plot([log10(Lmax) log10(Lmax)],[log10(min(Mo)) log10(max(Mo))],'k--')
plot([log10(Lmin) log10(Lmin)],[log10(min(Mo)) log10(max(Mo))],'k--')
plot([log10(min(L)) log10(max(L))],[log10(Momax) log10(Momax)],'k--')
title("LFE slip speed (mm/s)",'FontWeight','normal')
xticks=([1 1.5 2 2.5 3]);
xticklabels({'10 m',' ','100 m',' ','1 km'});
ax=gca;
ax.XTickLabelRotation = 0.;
xlabel('LFE source dimension')
ylabel('Characteristic LFE moment (Nm)')
yticks=([9.0 9.5 10 10.5 11 11.5 12]);
yticklabels({'10^{9}',' ','10^{10}',' ','10^{11}',' ','10^{12}'});
%ylabel('Characteristic LFE moment (Nm)')

keyboard
set(h,'PaperPosition',[0.25 0.25 8 10.5])
%orient landscape
% print('-dpdf','-painters','KinematicConstraints.pdf')



g=figure('Units','inches','Position',[1 1 8.5 11]);

subplot(3,2,1)
contour(log10(L),log10(Mo),VeqoVdyn,[1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2],'-k','ShowText','on')
hold on
plot([log10(Lmax) log10(Lmax)],[log10(min(Mo)) log10(max(Mo))],'k--')
plot([log10(Lmin) log10(Lmin)],[log10(min(Mo)) log10(max(Mo))],'k--')
plot([log10(min(L)) log10(max(L))],[log10(Momax) log10(Momax)],'k--')
title(["V_{eq} / V_{dyn} (multiply by:", ...
    "[^{0.01}/_{a}][^{1 MPa}/_{\sigma}][^{3.6 km/s}/_{c_s}][^{0.25 s}/_{T}] )"],'FontWeight','normal')
xticks=([1 1.5 2 2.5 3]);
xticklabels({'10 m',' ','100 m',' ','1 km'});
ax=gca;
ax.XTickLabelRotation = 0.;
xlabel('LFE source dimension')
yticks=([9.0 9.5 10 10.5 11 11.5 12]);
yticklabels({'10^{9}',' ','10^{10}',' ','10^{11}',' ','10^{12}'});
ylabel('Characteristic LFE moment (Nm)')

% subplot(3,2,2)
% contour(log10(L),log10(Mo),VeqoVdyn,[1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2],'-k','ShowText','on')
% hold on
% plot([log10(Lmax) log10(Lmax)],[log10(min(Mo)) log10(max(Mo))],'k--')
% plot([log10(Lmin) log10(Lmin)],[log10(min(Mo)) log10(max(Mo))],'k--')
% plot([log10(min(L)) log10(max(L))],[log10(Momax) log10(Momax)],'k--')
% title(["V_{eq} / V_{dyn} (multiply by:", ...
%     "[^{0.01}/_{a}][^{1 MPa}/_{\sigma}][^{3.6 km/s}/_{c_s}][^{0.25 s}/_{T}] )"],'FontWeight','normal')
% xticks=([1 1.5 2 2.5 3]);
% xticklabels({'10 m',' ','100 m',' ','1 km'});
% ax=gca;
% ax.XTickLabelRotation = 0.;
% xlabel('LFE source dimension')
% yticks=([9.5 10 10.5 11 11.5 12]);
% yticklabels({' ','10^{10}',' ','10^{11}',' ','10^{12}'});

subplot(3,2,2)
contour(log10(L),log10(Mo),DoDmin,[1e4 1e3 1e2 1e1 1e0 1e-1 1e-2],'-k','ShowText','on')
hold on
plot([log10(Lmax) log10(Lmax)],[log10(min(Mo)) log10(max(Mo))],'k--')
plot([log10(Lmin) log10(Lmin)],[log10(min(Mo)) log10(max(Mo))],'k--')
plot([log10(min(L)) log10(max(L))],[log10(Momax) log10(Momax)],'k--')
title(["Slip / min. slip (multiply by:", ...
    "[^{30 GPa}/_{G}][^{10\mum}/_{D_c}][2/{[a/(b-a)]^{1/2}} )"],'FontWeight','normal')
xticks=([1 1.5 2 2.5 3]);
xticklabels({'10 m',' ','100 m',' ','1 km'});
ax=gca;
ax.XTickLabelRotation = 0.;
xlabel('LFE source dimension')
yticks=([9.0 9.5 10 10.5 11 11.5 12]);
yticklabels({'10^{9}',' ','10^{10}',' ','10^{11}',' ','10^{12}'});
%ylabel('Characteristic LFE moment (Nm)')

subplot(3,2,3)
contour(log10(L),log10(Mo),LoLb,[1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2],'-k','ShowText','on')
hold on
plot([log10(Lmax) log10(Lmax)],[log10(min(Mo)) log10(max(Mo))],'k--')
plot([log10(Lmin) log10(Lmin)],[log10(min(Mo)) log10(max(Mo))],'k--')
plot([log10(min(L)) log10(max(L))],[log10(Momax) log10(Momax)],'k--')
title(["L / min. nucleation size 2L_b (multiply by:", ...
    "[^{b}/_{0.01}][^{\sigma}/_{1 MPa}][^{30 GPa}/_{G}][^{10\mum}/_{D_c}] )"],'FontWeight','normal')
xticks=([1 1.5 2 2.5 3]);
xticklabels({'10 m',' ','100 m',' ','1 km'});
ax=gca;
ax.XTickLabelRotation = 0.;
xlabel('LFE source dimension')
yticks=([9.0 9.5 10 10.5 11 11.5 12]);
yticklabels({'10^{9}',' ','10^{10}',' ','10^{11}',' ','10^{12}'});
ylabel('Characteristic LFE moment (Nm)')

subplot(3,2,4)
%note that 'AreaMin' is computed directly with 'SatAx'
contour(log10(L),log10(SatAx),AreaMin,[0.1 1 10],'-k','ShowText','on')
hold on
plot([log10(Lmax) log10(Lmax)],[log10(min(SatAx)) log10(max(SatAx))],'k--')
plot([log10(Lmin) log10(Lmin)],[log10(min(SatAx)) log10(max(SatAx))],'k--')
%plot([log10(min(L)) log10(max(L))],[log10(Momax) log10(Momax)],'k--')
title('Minimum tremor source area (km^2)','FontWeight','normal')
xticks=([1 1.5 2 2.5 3]);
xticklabels({'10 m',' ','100 m',' ','1 km'});
ax=gca;
ax.XTickLabelRotation = 0.;
xlabel('LFE source dimension')
ylabel('Instantaneous Sat. Level')
yticks=([0 0.5 1 1.5 2]);
yticklabels({'1',' ','10',' ','100'});
%ylabel('Characteristic LFE moment (Nm)')

subplot(3,2,5)
contour(log10(AmpoMin),log10(Mo),SatB,[1e0 1e2 1e4],'-k','ShowText','on')
hold on
plot([log10(min(AmpoMin)) log10(max(AmpoMin))],[log10(Momax) log10(Momax)],'k--')
title(["Instantaneous saturation (multiply by:", ...
    "[^{Low-amplitude Mo}/_{Bostock Mo_{min} of 5x10^{10}Nm}] )"],'FontWeight','normal')
ax=gca;
ax.XTickLabelRotation = 0.;
ax.XTick=([0 0.5 1]);
xticklabels({'1',' ','10'});
xlabel('Amplitude/minimum amplitude')
yticks=([9.0 9.5 10 10.5 11 11.5 12]);
yticklabels({'10^{9}',' ','10^{10}',' ','10^{11}',' ','10^{12}'});
ylabel('Characteristic LFE moment (Nm)')

subplot(3,2,6)
%amp, or moment increases as sqrt of event rate, i.e., Instantaneous saturation level
contour(log10(AmpoMin),log10(Mo),sqrt(SatB),[1e0 10 100],'-k','ShowText','on')
hold on
plot([log10(min(AmpoMin)) log10(max(AmpoMin))],[log10(Momax) log10(Momax)],'k--')
title(["Moment/Observed Moment (multiply by:", ...
    "[^{Low-amplitude Mo}/_{Bostock Mo_{min} of 5x10^{10}Nm}] )"],'FontWeight','normal')
ax=gca;
ax.XTickLabelRotation = 0.;
ax.XTick=([0 0.5 1]);
xticklabels({'1',' ','10'});
xlabel('Amplitude/minimum amplitude')
yticks=([9.0 9.5 10 10.5 11 11.5 12]);
yticklabels({'10^{9}',' ','10^{10}',' ','10^{11}',' ','10^{12}'});
%ylabel('Characteristic LFE moment (Nm)')

% subplot(3,2,5)
% contour(log10(L),log10(Mo),SatLevel,[1e1 1e0 1e-1 1e-2 1e-3 1e-4],'-k','ShowText','on')
% hold on
% plot([log10(Lmax) log10(Lmax)],[log10(min(Mo)) log10(max(Mo))],'k--')
% plot([log10(Lmin) log10(Lmin)],[log10(min(Mo)) log10(max(Mo))],'k--')
% plot([log10(min(L)) log10(max(L))],[log10(Momax) log10(Momax)],'k--')
% title('Fraction of slipping time each source is active','FontWeight','normal')

set(g,'PaperPosition',[0.25 0.25 8 10.5])
%orient landscape
% print('-dpdf','-painters','RSFConstraints.pdf')


    % function labels = expo(vals)
    % vals
    % labels=int2str(vals)
    % % labels(vals == 1e3) = "1.e3";
    % % labels(vals == 1e4) = "1.e4";
    % end
