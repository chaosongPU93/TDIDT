%Reads in seismograms and Bostock detections. Filters and rotates seisms, etc.  Stacks based on 
%times of B. detections, using shifts of up to 5 samples to maximize cross-correlation 
%of the central half of the window (at SSIB and SILB only). Puts SILB arrival at
%center of 5-sec window.  Isolates that time +/-12 samples.  Plots energy in that 
%0.6-sec window (both coherent and amp. squared) compared to the total in
%window.  Writes to .En files in BOSTOCK/NEW, with columns 7,8,9 being
%coherent energy (1st 2.2 sec, central 0.6 sec, last 2.2 sec), and 10,11,12
%being the same for amp squared.  

format short e
clear all
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');
scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite;
% h=figure('Position',[wid 1 wid hite]); %center

    timoffrot(1,:)=[2003 062 07 42 12 +086 +020  80 120  50]; %120 for 2-8 Hz
    timoffrot(2,:)=[2003 063 07 42 12 +086 +020  80 115  50];
    %timoffrot(3,:)=[2003 064 07 42 12 +086 +020  80 115  50];
    
    timoffrot(3,:)=[2004 196 07 42 12 +086 +020  80 115  50];
    timoffrot(4,:)=[2004 197 07 42 12 +086 +020  80 115  50];
    timoffrot(5,:)=[2004 198 07 42 12 +086 +020  80 115  50];
    timoffrot(6,:)=[2004 199 07 42 12 +086 +020  80 115  50];
   
    timoffrot(7,:)=[2005 254 07 42 12 +086 +020  80 115  50];
    timoffrot(8,:)=[2005 255 07 42 12 +086 +020  80 115  50];
    timoffrot(9,:)=[2005 256 07 42 12 +086 +020  80 115  50];
    timoffrot(10,:)=[2005 257 07 42 12 +086 +020  80 115  50];

%246  30304  1 3284.125 1.632 10
%  2  30304  2 3108.825 1.681 13
%  2  30304  2 3115.600 1.727 15
%  2  30304  2 3127.925 1.697 17
     bostname(1,:)='BOSTOCK/NEW/002_2003.062';
     bostname(2,:)='BOSTOCK/NEW/002_2003.063';
     %bostname(3,:)='BOSTOCK/NEW/002_2003.063';

     bostname(3,:)='BOSTOCK/NEW/002_2004.196';
     bostname(4,:)='BOSTOCK/NEW/002_2004.197';
     bostname(5,:)='BOSTOCK/NEW/002_2004.198';
     bostname(6,:)='BOSTOCK/NEW/002_2004.199';

     bostname(7,:)='BOSTOCK/NEW/002_2005.254';
     bostname(8,:)='BOSTOCK/NEW/002_2005.255';
     bostname(9,:)='BOSTOCK/NEW/002_2005.256';
     bostname(10,:)='BOSTOCK/NEW/002_2005.257';

% hi=6.;
% lo=1.5;
hi=8.;
lo=1;
% hi=12.;
% lo=0.5;
% hi=10.;
% lo=0.7;
% hi=9.;
% lo=0.6;
npo=2;
npa=2;

STAs=['SILB'
      'SSIB'];
%       'KLNB'
%       'TWKB'
%       'MGCB'];
nsta=size(STAs,1);
POLrots=[0 90 39 20;  %SILB
         6 85 33 86;];  %SSIB from Yajun
%          0 90  7 -5;  %KLNB
%          4 70 48 -4;  %MGCB
%          4 75 38 -27]; %TWKB
POLrots(:,2:3)=pi*POLrots(:,2:3)/180.;
%POLrots(:,1)=round(POLrots(:,1)*(100/40)); %Polaris stations at 100 sps; table has split time at 40 sps.
%POLrots(:,4)=round(POLrots(:,4)*(100/40)); %Polaris stations at 100 sps; table has split time at 40 sps.

SSIBtemp=load(['SSIB_',num2str(lo),'-',int2str(hi),'_npass',int2str(npa)]); %Read in the 100 sps templates 
SILBtemp=load(['SILB_',num2str(lo),'-',int2str(hi),'_npass',int2str(npa)]);
SSIBtemp=resample(SSIBtemp,2,5); %Downsample to 40 sps
SILBtemp=resample(SILBtemp,2,5);
maxlag=200;
xc=xcorr(SSIBtemp,SILBtemp,maxlag); %cross-correlate SSIB with SILB
[~,imax]=max(xc);
SSIBtemp(400:2000)=SSIBtemp(400+(imax-(maxlag+1)):2000+(imax-(maxlag+1))); %Shift SSIB to align
is=880; ilen=80;
xc=xcorr(SSIBtemp(is:is+ilen-1),SILBtemp(is:is+ilen-1)); %cross-correlate a shorter window (2 sec) (largest arrival at 914 @ 40sps)
[~,imax]=max(xc);
SSIBtemp(400:2000)=SSIBtemp(400+(imax-(ilen+1)):2000+(imax-(ilen+1))); %Again shift SSIB to align, if needed

[maxSILBtemp,imaxSILB]=max(SILBtemp); %where is the largest positive peak in the SILB template?
[minSILBtemp,iminSILB]=min(SILBtemp); %where is the largest negative peak in the SILB template?
for i=iminSILB:imaxSILB-1
    if SILBtemp(i)*SILBtemp(i+1)<0
        if abs(SILBtemp(i))<abs(SILBtemp(i+1)) %zero-crossing
            izeroSI=i;
        else
            izeroSI=i+1;
        end
    end
end
winlen=4.*40; 
SILBtemp=SILBtemp(izeroSI-winlen/2:izeroSI+winlen/2-1); %zero out the 1st 50 and last 15 samples of the template
SILBtemp(1:50)=0;
SILBtemp(end-15:end)=0;
SSIBtemp=SSIBtemp(izeroSI-winlen/2:izeroSI+winlen/2-1); %The max positive peak is in the middle of the window.
SSIBtemp(1:50)=0;
SSIBtemp(end-15:end)=0;
[maxSSIBtemp,imaxSSIB]=max(SSIBtemp);

minstretch=20;
maxstretch=80;
STAstr=zeros(nsta,maxstretch-minstretch+1,winlen);
for ista=1:nsta %2; %Dont do more stations yet
    STAtemp=load([STAs(ista,:),'_',num2str(lo),'-',int2str(hi),'_npass',int2str(npa)]); %now for the stretching.  This is 100 sps(?)
    STAstr(ista,:,:)=stretchtemps(STAtemp,minstretch,maxstretch,winlen);
end
% figure
% hold on
% dummy1=squeeze(STAstr(1,:,:));
% dummy2=squeeze(STAstr(2,:,:));
% for i=1:61
%     plot(dummy1(i,:),'r')
% end
% for i=1:61
%     plot(dummy2(i,:),'b')
% end
% xlim([60 100])
% pause

sps=40; %eventually
STAopt=zeros(nsta,sps*24*3600);
STAort=zeros(nsta,sps*24*3600);
nttot=0; %over all days
inwidth=0; %those with sufficiently high cc_max values
for nd=1:length(timoffrot(:,1))

    bostocks=load(bostname(nd,:));

    bostsec=3600*(bostocks(:,3)-1)+bostocks(:,4)-22.675+izeroSI*0.025; %-22.675 to correct for already-PGC-"corrected" 002 times!
    boststart=round(bostsec*sps)-winlen/2;
    bostend=boststart+winlen-1;

    year=timoffrot(nd,1);
    YEAR=int2str(year);
    jday=timoffrot(nd,2);
    if jday <= 9
        JDAY=['00',int2str(jday)];
    elseif jday<= 99
        JDAY=['0',int2str(jday)];
    else
        JDAY=int2str(jday);
    end
    MO=day2month(jday,year);
    direc=[YEAR,'/',MO,'/'];
    prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];
%     rotPGC=pi*(timoffrot(nd,8)-90)/180;
%     rotSSIB=pi*(timoffrot(nd,9)-90)/180;
%     rotSILB=pi*(timoffrot(nd,10)-90)/180;
%     SSIBsoff=timoffrot(nd,6);
%     SILBsoff=timoffrot(nd,7);
%     SSIBtoff=SSIBsoff/40;
%     SILBtoff=SILBsoff/40;

    %Read data:
    for ista=1:nsta
        STAEdat=[prename,'.',STAs(ista,:),'..HHE.D.SAC']; %BHE for permstas.
        STANdat=[prename,'.',STAs(ista,:),'..HHN.D.SAC'];
        [STAE,HdrDataSTA,tnuSTA,pobjSTA,timsSTA]=readsac(STAEdat,0,'l');
        [STAN,~,~,~,~]=readsac(STANdat,0,'l');
        tracelen=length(STAE);
        %cosine taper before filtering:
        x=(0:pi/200:pi/2-pi/200)';
        STAE(1:100)=sin(x).*STAE(1:100);
        STAN(1:100)=sin(x).*STAN(1:100);
        x=flipud(x);
        STAE(tracelen-99:tracelen)=sin(x).*STAE(tracelen-99:tracelen);
        STAN(tracelen-99:tracelen)=sin(x).*STAN(tracelen-99:tracelen);
        %Filter data:
        if year==2003 && jday<213
            [STAEf]=20.0e-3*bandpass(STAE,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
            [STANf]=20.0e-3*bandpass(STAN,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
        else
            [STAEf]=4.0e-3*bandpass(STAE,100,lo,hi,npo,npa,'butter'); 
            [STANf]=4.0e-3*bandpass(STAN,100,lo,hi,npo,npa,'butter'); 
        end
        STAEfd = resample(STAEf,2,5);
        STANfd = resample(STANf,2,5);
        STA=STAEfd+1i*STANfd;
        STAfastslow=STA*exp(-1i*POLrots(ista,2));
        STAslow=real(STAfastslow);
        STAfast=imag(STAfastslow);
        tracelen=length(STA);
        STAslow(20:tracelen-20)=STAslow(20+POLrots(ista,1):tracelen-20+POLrots(ista,1));
        STAsplitcorrected=(STAslow+1i*STAfast)*exp(1i*POLrots(ista,2));
        STAsplitcorrectedrot=STAsplitcorrected*exp(-1i*POLrots(ista,3));
        STAscrot=STAsplitcorrectedrot;

        STAsoff=POLrots(ista,4);
        if STAsoff > -1
            STAscrot(1:tracelen-STAsoff)=STAscrot(STAsoff+1:tracelen);
            STAscrot(tracelen-STAsoff+1:tracelen)=0;
        else
            STAscrot(-STAsoff+1:tracelen)=STAscrot(1:tracelen+STAsoff);
            STAscrot(1:-STAsoff)=0;
        end
        STAopt(ista,:)=real(STAscrot);
        STAort(ista,:)=imag(STAscrot);
    end

    %realPGC=real(PGCrot);
    realSILB=STAopt(1,:);
    realSSIB=STAopt(2,:);
    %realKLNB=STAopt(3,:);
    imagSILB=STAort(1,:);
    imagSSIB=STAort(2,:);
    %imagKLNB=STAort(3,:);

    mlag=7;
    lenb=length(bostsec);
    
    SILBtemplot=SILBtemp;
    SSIBtemplot=SSIBtemp;
    %iSSSIoff=0*(1:lenb);
    %maxSISS=0*(1:lenb);
    nt=0; %nt is "local", drops to zero at the end of each day
    nrow=6;
    mcol=4;
    timbef15=15*sps;
    timbef30=30*sps;
    for ifig=1:floor(lenb/(nrow*mcol))+1
        h=figure('Position',[scrsz(3)/10 scrsz(4)/10 4*scrsz(3)/5 9*scrsz(4)/10]);
        for n = 1:nrow
            for m = 1:mcol
                if nt < lenb %number in 1 day
                    nt=nt+1;
                    subplot(nrow,mcol,(n-1)*mcol+m,'align');
                    if bostend(nt) <= length(realSILB)
                        SILBtr=realSILB(boststart(nt):bostend(nt));
                        SSIBtr=realSSIB(boststart(nt):bostend(nt));
                        %KLNBtr=realKLNB(boststart(nt):bostend(nt));
                        SILBortr=imagSILB(boststart(nt):bostend(nt));
                        
                        SILBpre15=realSILB(boststart(nt)+round(winlen/4)-timbef15:boststart(nt)+round(winlen/4));
                        SILBpre30=realSILB(boststart(nt)+round(winlen/4)-timbef30:boststart(nt)+round(winlen/4));
                        SILBamp15(nttot+nt)=median(SILBpre15.*SILBpre15);
                        SILBamp30(nttot+nt)=median(SILBpre30.*SILBpre30);
                        
                        SISScc=xcorr(SILBtr(winlen/4:3*winlen/4),SSIBtr(winlen/4:3*winlen/4),mlag,'coeff');  %x-correlate over a shorter window
                        [maxSISS(nttot+nt),iSSSIoff]=max(SISScc);
                        iSSSIoff=iSSSIoff-(mlag+1); %This makes a logical lag of zero equal to zero
                        SSIBtr=realSSIB(boststart(nt)-iSSSIoff:bostend(nt)-iSSSIoff); %recenter SSIB
                        SSIBortr=imagSSIB(boststart(nt)-iSSSIoff:bostend(nt)-iSSSIoff); %recenter SSIB
                         
                        SSIBpre15=realSSIB(boststart(nt)+round(winlen/4)-iSSSIoff-timbef15:boststart(nt)+round(winlen/4)-iSSSIoff);
                        SSIBpre30=realSSIB(boststart(nt)+round(winlen/4)-iSSSIoff-timbef30:boststart(nt)+round(winlen/4)-iSSSIoff);
                        SSIBamp15(nttot+nt)=median(SSIBpre15.*SSIBpre15);
                        SSIBamp30(nttot+nt)=median(SSIBpre30.*SSIBpre30);
                        
                        %SIKLcc=xcorr(SILBtr(winlen/4:3*winlen/4),KLNBtr(winlen/4:3*winlen/4),mlag,'coeff');  %x-correlate over a shorter window
                        %[maxSIKL(nttot+nt),iKLSIoff]=max(SIKLcc);
                        %iKLSIoff=iKLSIoff-(mlag+1); %This makes a logical lag of zero equal to zero
                        %KLNBtr=realKLNB(boststart(nt)-iKLSIoff:bostend(nt)-iKLSIoff); %recenter KLNB
                        %KLNBortr=imagKLNB(boststart(nt)-iKLSIoff:bostend(nt)-iKLSIoff); %recenter KLNB
                        
                        SILBmax=max(abs(SILBtr));
                        SSIBmax=max(abs(SSIBtr));
                        %KLNBmax=max(abs(KLNBtr));
                        maxmax=max([SILBmax SSIBmax]); %KLNBmax]);
                        hold on
                        xc=xcorr(SILBtemp,SILBtr,10,'coeff');
                        [maxccSItemp(nttot+nt),imaxSI]=max(xc);
                        SILBtemplot(15:end-15)=SILBtemp(15+(imaxSI-(10+1)):end-15+(imaxSI-(10+1)));
                        plot(SILBtemplot*(maxmax/maxSILBtemp),'g--')
                        plot(SILBtr,'b')
                        plot(SSIBtr,'k')
                        %plot(KLNBtr,'r')
                        for i=1:61
                            dum=STAstr(1,i,:); %dummy version of the stretched template.  Why not "squeezed"?
                            dum(15:end-15)=dum(15+(imaxSI-(10+1)):end-15+(imaxSI-(10+1))); %shift that dummy version to the red trace
                            xc=xcorr(dum,SILBtr,8,'coeff'); %search for a modest additional shift
                            [maxccSIdipole(i),imaxSIdipole(i)]=max(xc); %for each stretched template, the ccmax and its offset
                        end 
                        [maxmaxSIdipole(nttot+nt),imaxSIstretch(nttot+nt)]=max(maxccSIdipole); %imaxSIstretch is the id of the best time-stretch
                        imaxmaxSIdipole(nttot+nt)=imaxSIdipole(imaxSIstretch(nttot+nt)); %imaxmaxSIdipole is the shift, from the red trace, of that best stretch
                        SILBtempsplot=squeeze(STAstr(1,imaxSIstretch(nttot+nt),:)); %set up a dummy version of the stretched dipole to plot
                        SILBtempsplot(25:end-25)=SILBtempsplot(25+(imaxSI-(10+1))+(imaxmaxSIdipole(nttot+nt)-(8+1)) ...
                            :end-25+(imaxSI-(10+1))+(imaxmaxSIdipole(nttot+nt)-(8+1)));
                        localmaxSI=max(SILBtr(winlen/2-20:winlen/2+20));
                        scaleSItemp=scaleSTA(SILBtr,SILBtempsplot,localmaxSI/maxSILBtemp); %scale the SILB template to match SILB
                        scaleSIstore(nttot+nt)=scaleSItemp; %store it
                        plot(scaleSItemp*SILBtempsplot,'r--')
                        %other stations below
                        xc=xcorr(SSIBtemp,SSIBtr,10,'coeff');
                        [maxccSStemp(nttot+nt),imaxSS]=max(xc);
                        SSIBtemplot(15:end-15)=SSIBtemp(15+(imaxSS-(10+1)):end-15+(imaxSS-(10+1)));
                        for i=1:61
                            dum=STAstr(2,i,:); %dummy version of the stretched template.  Why not "squeezed"?
                            dum(15:end-15)=dum(15+(imaxSS-(10+1)):end-15+(imaxSS-(10+1))); %shift that dummy version to the red trace
                            xc=xcorr(dum,SSIBtr,8,'coeff'); %search for a modest additional shift
                            [maxccSSdipole(i),imaxSSdipole(i)]=max(xc); %for each stretched template, the ccmax and its offset
                        end 
                        [maxmaxSSdipole(nttot+nt),imaxSSstretch(nttot+nt)]=max(maxccSSdipole); %imaxSSstretch is the id of the best time-stretch
                        imaxmaxSSdipole(nttot+nt)=imaxSSdipole(imaxSSstretch(nttot+nt)); %imaxmaxSSdipole is the shift, from the red trace, of that best stretch
                        SSIBtempsplot=squeeze(STAstr(1,imaxSSstretch(nttot+nt),:)); %set up a dummy version of the stretched dipole to plot
                        SSIBtempsplot(25:end-25)=SSIBtempsplot(25+(imaxSS-(10+1))+(imaxmaxSSdipole(nttot+nt)-(8+1)) ...
                            :end-25+(imaxSS-(10+1))+(imaxmaxSSdipole(nttot+nt)-(8+1)));
                        localmaxSS=max(SSIBtr(winlen/2-20:winlen/2+20));
                        scaleSStemp=scaleSTA(SSIBtr,SSIBtempsplot,localmaxSS/maxSSIBtemp); %scale the SSIB template to match SSIB
                        scaleSSstore(nttot+nt)=scaleSStemp; %store it
                        %plot(scaleSStemp*SSIBtempsplot,'r--')
                        %other stations above
                        bostmom(nttot+nt)=10.^(1.5*(bostocks(nt,5)+6));
                        title([YEAR(3:4),'.',JDAY,' ',int2str(bostsec(nt)),' ',int2str(iSSSIoff),'  M',num2str(bostmom(nttot+nt)), ...
                            '  cc',num2str(maxSISS(nttot+nt))],'fontsize',8)
                        text(100,-0.96*maxmax,[int2str(imaxSI-10+1),'  ',num2str(maxccSItemp(nttot+nt))],'fontsize',8)
                        text(20,-0.96*maxmax,[int2str(imaxmaxSIdipole(nttot+nt)-8+1),'  ',num2str((imaxSIstretch(nttot+nt)+19)/40)], ...
                            'fontsize',8)
                        ylim([-1.1*maxmax 1.1*maxmax])
                        xlim([0 160])
                        set(gca,'xtick',[40 80 120],'fontsize',8);
                        box on
                        if maxSISS(nttot+nt)>0.5 && maxccSItemp(nttot+nt)>0.6  %used to be 0.55 and 0.6
                            inwidth=inwidth+1;
                            width(inwidth,1:7)=[bostmom(nttot+nt) maxccSItemp(nttot+nt) (imaxSIstretch(nttot+nt)+19)/40 ...
                                scaleSIstore(nttot+nt) (imaxSSstretch(nttot+nt)+19)/40 scaleSSstore(nttot+nt) maxSISS(nttot+nt)];
                        end
                        traces(nttot+nt,1:160)=SILBtr;
                        traces(nttot+nt,161:320)=SSIBtr;
                        traces(nttot+nt,321:331)=[year+0.001*jday bostsec(nt) bostmom(nttot+nt) maxSISS(nttot+nt) ...
                            maxccSItemp(nttot+nt) maxccSStemp(nttot+nt) ...
                            (imaxSIstretch(nttot+nt)+19)/40 (imaxSSstretch(nttot+nt)+19)/40 ...
                            scaleSIstore(nttot+nt) scaleSSstore(nttot+nt) maxSISS(nttot+nt)*maxccSItemp(nttot+nt)];
                        traces(nttot+nt,332:491)=SILBortr;
                        traces(nttot+nt,492:651)=SSIBortr;
                        traces(nttot+nt,652:811)=scaleSItemp*SILBtempsplot;
                        %traces(nttot+nt,652:811)=KLNBtr;
                        %traces(nttot+nt,812:971)=KLNBortr;
                        stats(nttot+nt,:)=[year+0.001*jday bostsec(nt) bostmom(nttot+nt) maxSISS(nttot+nt) ...
                            0.5*(maxccSItemp(nttot+nt)+maxccSStemp(nttot+nt)) ...
                            0.5*((imaxSIstretch(nttot+nt)+19)/40+(imaxSSstretch(nttot+nt)+19)/40) ...
                            0.5*(scaleSIstore(nttot+nt)+scaleSSstore(nttot+nt)) ...
                            0.5*maxSISS(nttot+nt)*(maxccSItemp(nttot+nt)+maxccSStemp(nttot+nt)) ...
                            0.5*(SILBamp15(nttot+nt)+SSIBamp15(nttot+nt)) ...
                            0.5*(SILBamp30(nttot+nt)+SSIBamp30(nttot+nt))];
                    end
                end
            end
        end
        set(h,'PaperPosition',[0.25 0.25 8 10.5])
        orient landscape
        if ifig <= 9
            print(h,'-depsc',['BFIGS/',YEAR,'.',JDAY,'_',num2str(lo),'-',num2str(hi),'_',int2str(npa),'pass_',int2str(0),int2str(0),int2str(ifig),'.eps'])
        elseif ifig <= 99
            print(h,'-depsc',['BFIGS/',YEAR,'.',JDAY,'_',num2str(lo),'-',num2str(hi),'_',int2str(npa),'pass_',int2str(0),int2str(ifig),'.eps'])
        else
            print(h,'-depsc',['BFIGS/',YEAR,'.',JDAY,'_',num2str(lo),'-',num2str(hi),'_',int2str(npa),'pass_',int2str(ifig),'.eps'])
        end
    end
    nttot=nttot+nt
end

save('checkmomout.mat','traces')
save('checkmomstats.mat','stats')
dummySISS=0.5*(maxmaxSIdipole+maxmaxSSdipole);

figure
scatter(bostmom,maxccSItemp,12,maxSISS,'filled')
xlabel('Bostock magnitude')
ylabel('SILB CC_{max} with template')
ylim([0 1])
xlim([0.95*min(bostmom) 1.05*max(bostmom)])
caxis([0 1])
hcb=colorbar;
title(hcb,'SILB CC_{max} with SSIB','rotation',90,'position',[8.4 0.5])
box on
print('-depsc',['BFIGS/',num2str(lo),'-',num2str(hi),'Hz_',int2str(npa),'passSUM1.eps',])

figure
subplot(2,1,1,'align')
%scatter(bostmom,log(((imaxstretch+19)/40).*scaleSIstore.^2),12,maxSISS,'filled')
scatter(bostmom,log10(scaleSIstore),12,maxSISS,'filled')
xlabel('Bostock magnitude')
ylabel('log_{10}(dipole height)')
axis tight
caxis([0 1])
hcb=colorbar;
title(hcb,'SILB CC_{max} with SSIB','rotation',90,'position',[8.4 0.5])
box on
subplot(2,1,2,'align')
%scatter(bostmom,log(((imaxstretch+19)/40).*scaleSIstore.^2),12,maxSISS,'filled')
scatter(bostmom,(imaxSIstretch+19)/40,12,maxSISS,'filled')
xlabel('Bostock magnitude')
ylabel('dipole width')
caxis([0 1])
hcb=colorbar;
title(hcb,'SILB CC_{max} with SSIB','rotation',90,'position',[8.4 0.5])
box on
print('-depsc',['BFIGS/',num2str(lo),'-',num2str(hi),'Hz_',int2str(npa),'passSUM2.eps',])

h=figure('Position',[wid/3 1 wid hite]); 
xs=[11 12];
ys=[-0.5 0.5];
subplot(2,2,4,'align')
dummySI=log10(((imaxSIstretch+19)/40).^2.*scaleSIstore);
dummySS=log10(((imaxSSstretch+19)/40).^2.*scaleSSstore);
scatter(log10(bostmom),0.5*((dummySI+dummySS)-(median(dummySI)+median(dummySS))),10,dummySISS,'filled')
hold on
plot(xs,ys,'k-')
xlabel('Bostock moment')
ylabel('log_{10}(dipole height times width^2)')
%axis tight
xlim([10.8 12.4])
ylim([-1.5 1])
caxis([0 1])
hcb=colorbar('NorthOutside');
title(hcb,'CC_{max} with template')%,'rotation',0,'position',[8.4 0.5])
box on
subplot(2,2,3,'align')
dummySI=log10(((imaxSIstretch+19)/40).*scaleSIstore);
dummySS=log10(((imaxSSstretch+19)/40).*scaleSSstore);
scatter(log10(bostmom),0.5*((dummySI+dummySS)-(median(dummySI)+median(dummySS))),10,dummySISS,'filled')
hold on
plot(xs,ys,'k-')
xlabel('Bostock moment')
ylabel('log_{10}(dipole height times width)')
%axis tight
xlim([10.8 12.4])
ylim([-1.5 1])
caxis([0 1])
hcb=colorbar('NorthOutside');
title(hcb,'CC_{max} with template')%,'rotation',90,'position',[8.4 0.5])
box on
subplot(2,2,1,'align')
% dummySI=log10(((imaxSIstretch+19)/40).*scaleSIstore.^2);
% dummySS=log10(((imaxSSstretch+19)/40).*scaleSSstore.^2);
dummySI=log10(scaleSIstore);
dummySS=log10(scaleSSstore);
scatter(log10(bostmom),0.5*((dummySI+dummySS)-(median(dummySI)+median(dummySS))),10,dummySISS,'filled')
hold on
plot(xs,ys,'k-')
xlabel('Bostock moment')
ylabel('log_{10}(dipole height)')
%axis tight
xlim([10.8 12.4])
ylim([-1.5 1])
caxis([0 1])
hcb=colorbar('NorthOutside');
title(hcb,'CC_{max} with template')%,'rotation',90,'position',[8.4 0.5])
box on
subplot(2,2,2,'align')
dummySI=log10((imaxSIstretch+19)/40);
dummySS=log10((imaxSSstretch+19)/40);
scatter(log10(bostmom),0.5*((dummySI+dummySS)-(median(dummySI)+median(dummySS))),10,dummySISS,'filled')
hold on
xlabel('Bostock moment')
ylabel('log_{10}(dipole width)')
%axis tight
xlim([10.8 12.4])
ylim([-1.5 1])
caxis([0 1])
hcb=colorbar('NorthOutside');
title(hcb,'CC_{max} with template')%,'rotation',90,'position',[8.4 0.5])
box on
set(h,'PaperPosition',[0.25 0.25 8 10.5])
print('-depsc',['BFIGS/',num2str(lo),'-',num2str(hi),'Hz_',int2str(npa),'passSUM3.eps',])

figure
subplot(2,1,1,'align')
scatter(width(:,1),log(width(:,3)),18,width(:,2),'filled') %Bostock moment, stretch factor, CCmax SILB/SSIB
xlabel('Bostock magnitude')
ylabel('scaled width')
ylim([-log(2) log(2)])
title('min 4-s SILB-SSIB CC_{max} = 0.55')
cmax=max(width(:,2));
caxis([0.6 cmax])
hcb=colorbar;
title(hcb,'SILB CC_{max} with template','rotation',90,'position',[8.4 0.7])
box on
subplot(2,1,2,'align')
scatter(log10(width(:,4)),log(width(:,3)),18,width(:,2),'filled') %SILB amplitude, stretch factor, CCmax SILB/SSIB
xlabel('log_{10}(SILB amplitude)')
ylabel('scaled width')
ylim([-log(2) log(2)])
title('min 4-s SILB-SSIB CC_{max} = 0.55')
cmax=max(width(:,2));
caxis([0.6 cmax])
hcb=colorbar;
title(hcb,'SILB CC_{max} with template','rotation',90,'position',[8.4 0.7])
box on
print('-depsc',['BFIGS/',num2str(lo),'-',num2str(hi),'Hz_',int2str(npa),'passWIDTH1.eps',])

figure
subplot(2,1,1,'align')
scatter(log10(width(:,3)),log10(width(:,5)),18,width(:,2),'filled') %Bostock moment, stretch factor, CCmax SILB/SSIB
xlabel('Log_{10} SILB width')
ylabel('Log_{10} SSIB width')
axis equal
xlim([-0.3 0.3])
ylim([-0.3 0.3])
title('min 4-s SILB-SSIB CC_{max} = 0.55')
cmax=max(width(:,2));
caxis([0.6 cmax])
hcb=colorbar;
title(hcb,'SILB CC_{max} with template','rotation',90,'position',[8.4 0.7])
box on
subplot(2,1,2,'align')
scatter(log10(width(:,4)),log10(width(:,6)),18,width(:,2),'filled') %Bostock moment, stretch factor, CCmax SILB/SSIB
xlabel('Log_{10} SILB height')
ylabel('Log_{10} SSIB height')
axis equal
xlim([-4 -2])
ylim([-4 -2])
title('min 4-s SILB-SSIB CC_{max} = 0.55')
cmax=max(width(:,2));
caxis([0.6 cmax])
hcb=colorbar;
title(hcb,'SILB CC_{max} with template','rotation',90,'position',[8.4 0.7])
box on
print('-depsc',['BFIGS/',num2str(lo),'-',num2str(hi),'Hz_',int2str(npa),'passWIDTH2.eps',])

figure
subplot(2,1,1,'align')
scatter(log10(width(:,4)),log10(width(:,3)),18,width(:,2),'filled') %Bostock moment, stretch factor, CCmax SILB/SSIB
xlabel('Log_{10} SILB height')
ylabel('Log_{10} SILB width')
axis tight
title('min 4-s SILB-SSIB CC_{max} = 0.55')
cmax=max(width(:,2));
caxis([0.6 cmax])
hcb=colorbar;
title(hcb,'SILB CC_{max} with template','rotation',90,'position',[8.4 0.7])
box on
subplot(2,1,2,'align')
scatter(log10(width(:,6)),log10(width(:,5)),18,width(:,2),'filled') %Bostock moment, stretch factor, CCmax SILB/SSIB
xlabel('Log_{10} SSIB height')
ylabel('Log_{10} SSIB width')
axis tight
title('min 4-s SILB-SSIB CC_{max} = 0.55')
cmax=max(width(:,2));
caxis([0.6 cmax])
hcb=colorbar;
title(hcb,'SILB CC_{max} with template','rotation',90,'position',[8.4 0.7])
box on
print('-depsc',['BFIGS/',num2str(lo),'-',num2str(hi),'Hz_',int2str(npa),'passWIDTH3.eps',])

h=figure('Position',[wid/3 1 wid hite]); 
medSILBh=median(scaleSIstore);
medSILBw=median((imaxSIstretch+19)/40);
moms=10.^(1.5*(width(:,1)+6));
subplot(2,1,1,'align')
pSILBw=polyfit(log10(moms),log10(width(:,3)/medSILBw),1)
x=[10.5 12.5];
y=pSILBw(1)*x+pSILBw(2);
hold on
plot(x,y,'k--')
scatter(log10(moms),log10(width(:,3)/medSILBw),18,width(:,7),'filled') %Bostock moment, SILB stretch factor, maxCC SI-SS
xlabel('Log_{10} Bostock moment')
ylabel('Log_{10} SILB width')
axis equal
xlim([10.5 12.5])
%ylim([-log10(2) log10(2)])
title('min 4-s SILB-SSIB CC_{max} = 0.55')
cmax=max(width(:,7));
caxis([0.55 cmax])
hcb=colorbar;
title(hcb,'CC_{max} SILB-SSIB','rotation',90,'position',[8.4 0.7])
box on
subplot(2,1,2,'align')
pSILBh=polyfit(log10(moms),log10(width(:,4)/medSILBh),1)
x=[10.5 12.5];
y=pSILBh(1)*x+pSILBh(2);
hold on
plot(x,y,'k--')
scatter(log10(moms),log10(width(:,4)/medSILBh),18,width(:,7),'filled') %Bostock moment, SILB height, maxCC SI-SS
xlabel('Log_{10} Bostock moment')
ylabel('Log_{10} SILB amplitude')
axis equal
xlim([10.5 12.5])
% ylim([-])
%title('min 4-s SILB-SSIB CC_{max} = 0.55')
cmax=max(width(:,7));
caxis([0.55 cmax])
hcb=colorbar;
title(hcb,'CC_{max} SILB-SSIB','rotation',90,'position',[8.4 0.7])
box on
print('-depsc',['BFIGS/',num2str(lo),'-',num2str(hi),'Hz_',int2str(npa),'passWIDTH4.eps',])

h=figure('Position',[wid/3 1 wid hite]); 
medSSIBh=median(scaleSSstore);
medSSIBw=median((imaxSSstretch+19)/40);
moms=10.^(1.5*(width(:,1)+6));
subplot(2,1,1,'align')
pSSIBw=polyfit(log10(moms),log10(width(:,5)/medSSIBw),1)
x=[10.5 12.5];
y=pSSIBw(1)*x+pSSIBw(2);
hold on
plot(x,y,'k--')
scatter(log10(moms),log10(width(:,5)/medSSIBw),18,width(:,7),'filled') %Bostock moment, SILB stretch factor, maxCC SI-SS
xlabel('Log_{10} Bostock moment')
ylabel('Log_{10} SSIB width')
axis equal
xlim([10.5 12.5])
%ylim([-log10(2) log10(2)])
title('min 4-s SILB-SSIB CC_{max} = 0.55')
cmax=max(width(:,7));
caxis([0.55 cmax])
hcb=colorbar;
title(hcb,'CC_{max} SILB-SSIB','rotation',90,'position',[8.4 0.7])
box on
subplot(2,1,2,'align')
pSSIBh=polyfit(log10(moms),log10(width(:,6)/medSSIBh),1)
x=[10.5 12.5];
y=pSSIBh(1)*x+pSSIBh(2);
hold on
plot(x,y,'k--')
scatter(log10(moms),log10(width(:,6)/medSSIBh),18,width(:,7),'filled') %Bostock moment, SILB height, maxCC SI-SS
xlabel('Log_{10} Bostock moment')
ylabel('Log_{10} SSIB amplitude')
axis equal
xlim([10.5 12.5])
% ylim([-])
%title('min 4-s SILB-SSIB CC_{max} = 0.55')
cmax=max(width(:,7));
caxis([0.55 cmax])
hcb=colorbar;
title(hcb,'CC_{max} SILB-SSIB','rotation',90,'position',[8.4 0.7])
box on
print('-depsc',['BFIGS/',num2str(lo),'-',num2str(hi),'Hz_',int2str(npa),'passWIDTH5.eps',])



