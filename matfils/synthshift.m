%Reads nsta(=3) template stacks. 
%Builds a synthetic seismogram for each from a user-defined window (ginput).
%shifts 2 and 3 w.r.t. 1 over a pre-defined range.
%This verison is very similar to 'chaosynth.m' and 'synthshift2023.m', the main
%difference is that 'chaosynth' uses 'plaw3b' for generating src arrivals, while 
%'synthshift.m' uses 'plaw3c', and 'synthshift2023.m' uses 'plaw3d'
clear
% close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');
scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite;

datfiles=['PGCopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh0.5-15Hz_le14sh';  %the templates, 0.5-15 Hz at 100 Hz sampling rate
          'SSIopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh0.5-15Hz_le14sh';
          'SILopt_002_1-15Hz_2pass_100sps.le14sh0.75-15Hz_le14sh0.5-15Hz_le14sh'];
stas=['PGC ';
      'SILB';
      'SSIB'];
sps=100;
spsm1=1/sps;
nsta=size(stas,1);

%%%The below aligns the templates by x-correlation
for ista=1:nsta
    STA(ista,:)=load(datfiles(ista,:));
end
figure
hold on
plot(STA(1,:),'r')
plot(STA(2,:),'b')
plot(STA(3,:),'k')

[maxses,imaxses]=max(STA,[],2);
[minses,iminses]=min(STA,[],2);
spread=maxses-minses;
zcrosses=round(0.5*(imaxses+iminses));  % rough, assuming symmetry, Chao 2021/07/16
sampbef=4*sps;
sampaft=10*sps;
is=zcrosses-sampbef;
ie=zcrosses+sampaft;
for ista=1:nsta
    STAtmp(ista,:)=STA(ista,is(ista):ie(ista));
end
mshift=10;
tempxc(1,:)=xcorr(STAtmp(2,:),STAtmp(3,:),mshift,'coeff');
tempxc(2,:)=xcorr(STAtmp(1,:),STAtmp(2,:),mshift,'coeff'); %shift STAtmp(3,:) to right for positive values
tempxc(3,:)=xcorr(STAtmp(1,:),STAtmp(3,:),mshift,'coeff'); %shift STAtmp(2,:) to right for positive values
[~,imax]=max(tempxc,[],2);
imax=imax-(mshift+1) %This would produce a slightly different shift, if filtered seisms were used.
for ista=2:nsta
    STAtmp(ista,mshift+1:end-(mshift+1))=STAtmp(ista,mshift+1-imax(ista):end-(mshift+1)-imax(ista));
end
for ista=1:nsta
    STAtmp(ista,:)=STAtmp(ista,:)/spread(ista);
end
figure
hold on
plot(STAtmp(1,:),'r')
plot(STAtmp(2,:),'b')
plot(STAtmp(3,:),'k')
%%%The above aligns the templates by x-correlation

%%%The below determines the template window, and makes sure it has zero mean with tapered ends
% for tapering, could also use 'tukeywin', Chao 2021/07/16
% greenlen is its length
tlim=[237 1136] %just for consistency; could go back to picking
greenlen=tlim(2)-tlim(1)+1;
green=zeros(nsta,greenlen);
smooth=sps/2;
for ista=1:nsta
    green(ista,:)=STAtmp(ista,tlim(1):tlim(2));
    green(ista,:)=green(ista,:)-mean(green(ista,:));    % remove mean
    x=(0:pi/(2*smooth):pi/2-pi/(2*smooth));    
    green(ista,1:smooth)=sin(x).*green(ista,1:smooth);  % left taper
    x=fliplr(x);
    green(ista,end-(smooth-1):end)=sin(x).*green(ista,end-(smooth-1):end);  % right taper
    green(ista,:)=green(ista,:)-mean(green(ista,:));    % remove mean
    green(ista,:)=green(ista,:)/max(abs(green(ista,:)));    % normalize
end
mean(green,2)
subplot(2,1,1)
hold on
plot(green(1,:),'r')
plot(green(2,:),'b')
plot(green(3,:),'k')
mx=max([abs(green(1,:)) abs(green(2,:))]);
xlim([0 greenlen])
ylim([-mx mx])
box on

for ista=1:nsta
%     greenf(ista,:)=Bandpass(green(ista,:),100,1.25,6.5,2,2,'butter');   % change 'bandpass' to 'Bandpass'
    greenf(ista,:)=Bandpass(green(ista,:),100,1.8,18,2,2,'butter');   % change 'bandpass' to 'Bandpass'
end
subplot(2,1,2)
hold on
plot(greenf(1,:),'r')
plot(greenf(2,:),'b')
plot(greenf(3,:),'k')
mx=max([abs(greenf(1,:)) abs(greenf(2,:))]);
%%%The above determines the template window, and makes sure it has zero mean with tapered ends

%%%The below cross-correlates the templates (Green's functions), just FYI
% running CC using a window length of 'cclen', Chao 2021/07/26
    cclen=50;
    samples=1:length(greenf(1,:));
    ST1auto=greenf(1,:).*greenf(1,:);
    ST1sq=cumsum(ST1auto);
    ST2auto=greenf(2,:).*greenf(2,:);
    ST2sq=cumsum(ST2auto);
    ST3auto=greenf(3,:).*greenf(3,:);
    ST3sq=cumsum(ST3auto);
    STA12tr=greenf(1,:).*greenf(2,:);            
    STA13tr=greenf(1,:).*greenf(3,:);
    STA23tr=greenf(2,:).*greenf(3,:);
    cumsumSTA12tr=cumsum(STA12tr);
    cumsumSTA13tr=cumsum(STA13tr);
    cumsumSTA23tr=cumsum(STA23tr);
    ST12num=cumsumSTA12tr(cclen+1:length(greenf(1,:)))-cumsumSTA12tr(1:length(greenf(1,:))-cclen);
    ST13num=cumsumSTA13tr(cclen+1:length(greenf(1,:)))-cumsumSTA13tr(1:length(greenf(1,:))-cclen);
    ST23num=cumsumSTA23tr(cclen+1:length(greenf(1,:)))-cumsumSTA23tr(1:length(greenf(1,:))-cclen);
    ST1den=ST1sq(cclen+1:length(greenf(1,:)))-ST1sq(1:length(greenf(1,:))-cclen);
    ST2den=ST2sq(cclen+1:length(greenf(1,:)))-ST2sq(1:length(greenf(1,:))-cclen);
    ST3den=ST3sq(cclen+1:length(greenf(1,:)))-ST3sq(1:length(greenf(1,:))-cclen);
    ST12n=ST12num./realsqrt(ST1den.*ST2den);
    ST13n=ST13num./realsqrt(ST1den.*ST3den);
    ST23n=ST23num./realsqrt(ST2den.*ST3den);
    alln=(ST12n+ST13n+ST23n)/3;
    %alln(alln<0)=-10^4*yma; %just so they don't plot.
    plot(samples(cclen/2+1:greenlen-cclen/2),mx*alln,'co','markersize',2)
xlim([0 greenlen])
ylim([-mx mx])
box on
xc23=xcorr(greenf(2,:),greenf(3,:),10,'coeff');
xc13=xcorr(greenf(1,:),greenf(3,:),10,'coeff');
xc12=xcorr(greenf(1,:),greenf(2,:),10,'coeff');
[ccmax23,imax23]=max(xc23)
[ccmax13,imax13]=max(xc13)
[ccmax12,imax12]=max(xc12)
%%%The above cross-correlates the template (Green's functions), just FYI

% % tester=green(315:345);
% % subplot(2,1,2,'align')
% % hold on
% % plot(tester,'b')
% % plot(greenf(315:345),'r')
% % testlen=length(tester);
% % autotest=tester.*tester;
% % cumtest=cumsum(autotest);
% % tester2=cumtest(end);
% % offset=20;

keyboard

%%
%%%%%%%%%%%%%%%%%%
%distr='PL'
%distr='LN'
distr='UN'
%distr='DI'
%distr='EX'
%print('-depsc',['b=',int2str(b),'_template.eps'])
%Twin=3*3600+3+2*ceil(greenlen/sps); %Early and late portions will be deleted. Twin includes only the early portion. 
Twin=0.5*3600+3+2*ceil(greenlen/sps); %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
%Twin=4*3600+3+2*ceil(greenlen/sps); %For DI. 
%Twin=6*3600+3+2*ceil(greenlen/sps); %For DI. 
%Twin=1*3600+3+2*ceil(greenlen/sps); %For DI. 
winlen=Twin*sps+1;
skiplen=greenlen;
%04/25/2022, 0.4s as the duration seems a bit off, more like 0.5s from the broadband template  
satn=2.5*Twin %a single peak is ~20 samples wide; maybe a little less (at 100 sps).
    %Twin is window duration in seconds. Events can fall within Twin of
    %the start, but the synthetics will go to Twin*(sample rate)+Greenlen to
    %avoid checking for subscript overrrun.  When "synths" are written from "synth" in the
    %subroutine, Greenlen from the start and end will not be written.
mjig=0; %5; %max number of samples to jigger 2 and 3 w.r.t. 1.  AT 100 SPS! OBSOLETE? But in file name still
fracelsew=0.5; %0.25 %0.5; %0.6; %The ratio of elsewhere events to local events.  Zero for Discrete Ide.
%writes=[1e7]; %Just a single entry for Discrete Ide.  Value not relevant.
nsat=[0.1 0.4 1 2 4 10 100];
writes=round(nsat*satn) %in degree of saturation
rng('default');
%synth=2.e-7*(rand(nsta,winlen+greenlen+2*mjig)-0.5); %+2*mjig a little extra, for jiggering 2nd & 3rd stations.
synth=2.e-7*(rand(nsta,winlen+greenlen+2*10)-0.5); %+2*10 a little extra, for jiggering 2nd & 3rd stations.% for ista=1:nsta
%     synth(ista,:)=1.e-7*sin(1:winlen+greenlen+2*mjig); %+2*mjig a little extra, for jiggering 2nd & 3rd stations.
% end

nouts=length(writes);
seed=round(writes(1)/5e3); %for random number generator
if strcmp(distr,'PL') || strcmp(distr,'UN') %b>150 for uniform
    b=999. %>150 for uniform size distribution
   %[synths,mommax]=plaw3(writes,winlen,skiplen,synth,green,b,mjig,seed); 
   %[synths,mommax]=plaw3b(writes,winlen,skiplen,synth,green,b,mjig,fracelsew,seed); 
    %maxdist=1.0; %in km; define a circular region on the fault to draw sources from
    %maxdist=9; %This is shorthand for 'NW'; used only to name output file in this case..
    xygrid=load('xygridArmb'); %Made from /ARMMAP/MAPS/2021/interpgrid.m; format PGSS, PGSI, dx, dy
    size(xygrid)
    %sourcestyle='ellipse';
    sourcestyle='rectangle';
    if strcmp(sourcestyle,'rectangle') %The following rotates 45˚ and looks for -5 < x' < 1 and 2.25 < y' < 3.5
        angrot=-45*pi/180;
        loc=complex(xygrid(:,3),xygrid(:,4));
        locrot=loc*exp(1i*angrot);
        xygrid(:,5)=real(locrot);
        xygrid(:,6)=imag(locrot);
        ll=[-5 2.25]; ur=[1 3.5];
        xygrid(xygrid(:,5)<ll(1) | xygrid(:,5)>ur(1) | xygrid(:,6)<ll(2) | xygrid(:,6)>ur(2),:)=[]; %This limits xygrid to the desired range.  PGSS, PGSI, x, y, x', y'
    elseif strcmp(sourcestyle,'ellipse') %The following adds 0.5 km to y, rotates 45˚, and looks for an elliptical region with major semi-axis y' < 3 km and minor semi-axis < 2 km
        angrot=-45*pi/180;
        shiftor=[0 -0.5]; %(in km)
        loc=complex(xygrid(:,3)-shiftor(1),xygrid(:,4)-shiftor(2));
        locrot=loc*exp(1i*angrot);
        xygrid(:,5)=real(locrot);
        xygrid(:,6)=imag(locrot);
        xaxis=0.5; %(semi-, in km)
        yaxis=0.5; %(semi-, in km)
        xygrid(sqrt((xygrid(:,5)/xaxis).^2 + (xygrid(:,6)/yaxis).^2) > 1,:)=[]; %This limits xygrid to the desired range.  PGSS, PGSI, x, y, x', y'
    end
    fid = fopen('tmp.grd','w');
    fprintf(fid,'%8.3f %8.3f %7.2f %7.2f %8.3f %8.3f\n',xygrid');
    fclose(fid);
    size(xygrid)
    [synths,mommax,sources]=plaw3c(writes,winlen,skiplen,synth,green,b,xygrid,sps,fracelsew,seed); 
end
if strcmp(distr,'LN')
    sig=1.;
    [synths,mommax]=lognorm(writes,winlen,skiplen,synth,green,sig,mjig,fracelsew,seed); 
end
if strcmp(distr,'DI')
    stepsize=5; %Sets the step size for the random walk
    frac=1; %Sets the s.d. of the variation as frac*sqrt(moment rate)
    desiredmax=3000; %2000; %Used to shift the center of the Gaussian, to limit large moment rates. 
    alpha=1; %Center of the Gaussian is shifted as -sd*(STF(i-1)/desiredmax)^alpha
    av=30*sps; %in samples, the length of time over which the average moment rate is determined
    
    stepsize=5; %Sets the step size for the random walk
    frac=1; %Sets the s.d. of the variation as frac*sqrt(moment rate)
    desiredmax=1e9; %2000; %Used to shift the center of the Gaussian, to limit large moment rates. 
    alpha=1; %Center of the Gaussian is shifted as -sd*(STF(i-1)/desiredmax)^alpha
    av=1; %in samples, the length of time over which the average moment rate is determined
    
    [synths,cumSTF]=DIde(writes,winlen,skiplen,synth,stepsize,desiredmax,frac,alpha,av,green,mjig,seed);
end
% if strcmp(distr,'EX')                                             %Not operating yet.
%     [synths,mommax]=expo2(writes,winlen,skiplen,synth,green,jig); %Not operating yet.
% end                                                               %Not operating yet.
% [AmpEst,fname1,fname2]=scale2(writes,satn,distr,b,Twin);
% fid = fopen([fname1,'_jig',int2str(jig)],'a');
% for inwrites=1:length(writes)
%     n=writes(inwrites);
%     xc=xcorr(synths(1,inwrites,:),synths(2,inwrites,:),jig+1,'coeff');
%     fprintf(fid,'%9i %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n', n, n/satn, ...
%        AmpEst(inwrites),max(xc),mommax(inwrites),max(synths(1,inwrites,:)),mean(abs(synths(1,inwrites,:))),mean(synths(1,inwrites,:))); ...
% end
% fclose(fid);
if strcmp(distr,'DI')
    for inwrites=1:length(writes)
        n=cumSTF;
        fid = fopen(['3STAS/STAS.',distr,'.dt',int2str(stepsize),'.dmax',int2str(desiredmax),'.av',int2str(av),'.alp',num2str(alpha),'.frac',num2str(frac),'.',int2str(sps),'sps.mjig',int2str(mjig),'else',num2str(fracelsew,2),'_N10e',num2str(log10(n),3)],'w');
        towrite=squeeze(synths(:,inwrites,:));
        fprintf(fid,'%13.5e %13.5e %13.5e\n', towrite);
        fclose(fid);
    end
end
if strcmp(distr,'LN')
    for inwrites=1:length(writes)
        n=writes(inwrites);
        fid = fopen(['3STAS/STAS.',distr,'.sig',num2str(sig,2),'.',int2str(sps),'sps.mjig',int2str(mjig),'else',num2str(fracelsew,2),'_N10e',num2str(log10(n),3)],'w');
        towrite=squeeze(synths(:,inwrites,:));
        fprintf(fid,'%13.5e %13.5e %13.5e\n', towrite);
        fclose(fid);
    end
end
if strcmp(distr,'UN')
    for inwrites=1:length(writes)
        n=writes(inwrites);
        if strcmp(sourcestyle,'ellipse')
            fid = fopen(['3STAS/STAS.',distr,'.',int2str(sps),'sps.',sourcestyle(1:3),'_',num2str(xaxis),'-',num2str(yaxis),'.else',num2str(fracelsew,2),'nsat',num2str(nsat(inwrites))],'w');
        elseif strcmp(sourcestyle,'rectangle')
            fid = fopen(['3STAS/STAS.',distr,'.',int2str(sps),'sps.',sourcestyle(1:3),'_',num2str(ll(1)),num2str(ll(2)),num2str(ur(1)),num2str(ur(2)),'.else',num2str(fracelsew,2),'nsat',num2str(nsat(inwrites))],'w');
        end
        towrite=squeeze(synths(:,inwrites,:));
        fprintf(fid,'%13.5e %13.5e %13.5e\n', towrite);
        fclose(fid);
        if strcmp(sourcestyle,'ellipse')
            fid = fopen(['3STAS/STAS.',distr,'.',int2str(sps),'sps.',sourcestyle(1:3),'_',num2str(xaxis),'-',num2str(yaxis),'.else',num2str(fracelsew,2),'nsat',num2str(nsat(inwrites)),'_sources'],'w');
        elseif strcmp(sourcestyle,'rectangle')
            fid = fopen(['3STAS/STAS.',distr,'.',int2str(sps),'sps.',sourcestyle(1:3),'_',num2str(ll(1)),num2str(ll,2),num2str(ur(1)),num2str(ur,2),'.else',num2str(fracelsew,2),'nsat',num2str(nsat(inwrites)),'_sources'],'w');
        end
        truncsources=sortrows(sources(1:n,:),1);
        fprintf(fid,'%13i %5i %4i\n', truncsources');
        fclose(fid);
    end
end

%% a simple plot the synthetics
figure
nrow = length(writes)+1;
ncol = 1;
subplot(nrow,ncol,1)
hold on
tmp = synth;
plot(tmp(1,:),'r')
plot(tmp(2,:),'b')
plot(tmp(3,:),'k')
xlim([0 0.4e4])

for i = 1: length(writes)
subplot(nrow,ncol,i+1)
hold on
tmp = synths(:,i,:);
plot(tmp(1,:),'r')
plot(tmp(2,:),'b')
plot(tmp(3,:),'k')
xlim([0 0.2e4])
end
