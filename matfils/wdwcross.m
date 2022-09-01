% Reads in SSIB and SILB seismograms; filters and rotates them, etc.  For specified 
% windows during one day, finds the best offset for windows of a specified length (e.g. 7.5 s).  Timeshifts up to 
% +/- 5 samples gets most central "hot spot". Uses an extra "cushion" to each side of that s.t. if the max (non-normalized)
% cc offset falls in the cushion, that short window does not contribute.  For all other windows, sum the running dot product
% (to emphasize the coherent energy).  Also the sqrt of the sum of the squares of both stations to get all the radiated energy,
% for comparison.  As of 9/29/18 also takes the running sum of the sqrt of the dot product (after zeroing out all negative values).
% In the plots, black is everything, red is coherent, solids are radiated energy, dashed are amplitude squared amplitude squared),
% dashed is just summed abs(amplitude).
% Make sure any scale factors are consistent with wdwBostockDisp.m!
format short e
clear all
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

   %bostocks=load('BOSTOCK/NEW/002-246_2003.062');
   %bostocks=load('BOSTOCK/NEW/002-246_2003.063');
   %bostocks=load('BOSTOCK/NEW/002-246_2003.064');
   %bostocks=load('BOSTOCK/NEW/002-246_2004.196');
   %bostocks=load('BOSTOCK/NEW/002-246_2004.197');
   %bostocks=load('BOSTOCK/NEW/002-246_2004.198');
   %bostocks=load('BOSTOCK/NEW/002-246_2004.199');
   %bostocks=load('BOSTOCK/NEW/002-246_2005.254');
   %bostocks=load('BOSTOCK/NEW/002-246_2005.255');
   %bostocks=load('BOSTOCK/NEW/002-246_2005.256');
   %bostocks=load('BOSTOCK/NEW/002-246_2005.257');
   %bostsec=3600*(bostocks(:,3)-1)+bostocks(:,4)+4*0.025; %(the latter to correct for the diff between Bostock and Rubin filtering.  But why isn't it 8 samples?)

sps=40;
days=[2003.062;
      2003.063;
      %2003.064;
      2004.196;
      2004.197;
      2004.198;
      2004.199;
      2005.254;
      2005.255;
      2005.256;
      2005.257];
ndays=length(days);
  
bads=[62200 62450; 
      0.025 0.025;
      %0.025 0.025;
      2650 2675;
      43650 44100;
      0.025 0.025;
      0.025 0.025;
      0.025 0.025;
      0.025 0.025;
      18875 18925;
      70620 70645];
bads=bads*sps;

infiles=['ARMMAP/MAPS/map2003.062.86.20.80.120.50_2-8-ms12-4s.dx.onepeak0.5.inbox'; %These are closer to the origin?  For 062, 820 instead of 1577.
         'ARMMAP/MAPS/map2003.063.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox';
         %'ARMMAP/MAPS/map2003.064.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox';
         'ARMMAP/MAPS/map2004.196.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox'; 
         'ARMMAP/MAPS/map2004.197.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox';
         'ARMMAP/MAPS/map2004.198.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox';
         'ARMMAP/MAPS/map2004.199.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox';
         'ARMMAP/MAPS/map2005.254.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox';
         'ARMMAP/MAPS/map2005.255.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox';
         'ARMMAP/MAPS/map2005.256.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox';
         'ARMMAP/MAPS/map2005.257.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5.inbox'];

outfiles=['ARMMAP/MAPS/map2003.062.86.20.80.120.50_2-8-ms12-4s.dx.onepeak0.5';
         'ARMMAP/MAPS/map2003.063.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5';
         %'ARMMAP/MAPS/map2003.064.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5';
         'ARMMAP/MAPS/map2004.196.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5'; 
         'ARMMAP/MAPS/map2004.197.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5';
         'ARMMAP/MAPS/map2004.198.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5';
         'ARMMAP/MAPS/map2004.199.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5';
         'ARMMAP/MAPS/map2005.254.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5';
         'ARMMAP/MAPS/map2005.255.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5';
         'ARMMAP/MAPS/map2005.256.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5';
         'ARMMAP/MAPS/map2005.257.86.20.80.115.50_2-6-ms12-4s.dx.onepeak0.5'];
%As of 08.31.18, the lines that use these files (plotreflines) have been commented out.

hi=8.;
lo=1;
npo=2;
npa=2;

STAs=['SILB'
      'SSIB'];
%       'KLNB'
%       'TWKB'
%       'MGCB'];
nsta=size(STAs,1);
POLrots=[0 90 39 20;  %SILB
         6 85 33 86];  %SSIB from Yajun
         %0 90  7 00;  %KLNB
         %4 70 48 00;  %MGCB
         %4 75 38 00]; %TWKB
POLrots(:,2:3)=pi*POLrots(:,2:3)/180.;

cushion=20; %20; %5; %where, if the max cc coefficient is for that offset, you don't add the dot product.
nshift=6+cushion; %max (-) shift of SSIB relative to SILB (NW corner)
mshift=5+cushion; %max shift of SSIB relative to SILB (SE edge)
%cclen=floor(40/lo)+mod(floor(40/lo),2) %length of running cc window

scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite;

tottotcross=0.;
tottotsum=0.;
tottotcohamp=0.;
tottotamp=0.;
for iday=1:ndays
    year=floor(days(iday));
    YEAR=int2str(year);
    jday=round(1000*(days(iday)-year))
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
    realSILB=STAopt(1,:);
    realSSIB=STAopt(2,:);
    tracelenSILB=length(realSILB);
    tracelenSSIB=length(realSSIB);
    small=1.e-8;
    realSSIB(bads(iday,1):bads(iday,2))=small;
    realSILB(bads(iday,1):bads(iday,2))=small;
    
    rubins=load(infiles(iday,:));
    rubouts=load(outfiles(iday,:));
    LIA=ismember(rubouts(:,1),rubins(:,1));
    rubouts(LIA,:)=[];

    offsetstart=100;
    detects=windows(days(iday));
    %wins=outsidewindows(days(iday));
    %wins=[offsetstart 86400-offsetstart];
    wins=windows(days(iday));
    nwins=size(wins,1);
    %shortwin=30*sps;
    shortwin=7.5*sps; %15*sps;
    totlen=0;
    totcross=0;
    totsum=0;
    totamp=0;
    totcohamp=0;
    for iwin=1:nwins %over all the windows in one day
        istart=wins(iwin,1)*sps; %start sample
        iend=wins(iwin,2)*sps; %end sample
        winlen=iend-istart+1;
        amplitude=0.5*(abs(realSILB(istart:iend))+abs(realSSIB(istart:iend)));
        cumamp=cumsum(amplitude);
        SILBauto=realSILB(istart:iend).*realSILB(istart:iend);
        SILB2=cumsum(SILBauto);
        %SSIBauto=realSSIB(istart-nshift:iend+mshift).*realSSIB(istart-nshift:iend+mshift); %Longer b.c. of mshift.
        SSIBauto=realSSIB(istart:iend).*realSSIB(istart:iend); %Not as precise, but should be OK
        SSIB2=cumsum(SSIBauto);
        denom=sqrt(SILBauto.*SSIBauto);
        cumdenom=cumsum(denom);
	SISSx=zeros(winlen,nshift+mshift+1);
	for k=-nshift:mshift; %cross-correlate over all the shifts
            SISSx(:,k+nshift+1)=realSILB(istart:iend).*realSSIB(istart+k:iend+k); %one column of a.*b for each integer offset
        end
        tmpx=SISSx;
        tmpx(SISSx<0)=0;
        rootSISSx=sqrt(tmpx);
        if nshift+mshift==0 %if no jiggering
            cumsumSISS=cumsum(SISSx,1); 
            cumrootSISS=cumsum(rootSISSx,1); 
        else
            cumsumSISS=cumsum(SISSx); %the more general case.  Sums each column of a.*b for each integer offset
            cumrootSISS=cumsum(rootSISSx); %the more general case.  Sums each column of a.*b for each integer offset
        end
        SISSxamalg=0*(1:winlen); %SISSxamalg zeroed for EACH WINDOW!
        rootSISSxamalg=0*(1:winlen); %likewise for rootSISSxamalg 
        tmp=0*(1:winlen); %Likewise for tmp
        ie=0; %need to initialize for when ie doesn't change (winlen<shortlen)
        for j=1:floor(winlen/shortwin) %for each (e.g.) 15-s sub-window
            is=1+(j-1)*shortwin;
            ie=is+shortwin-1;
            [~,kmax]=max(cumsumSISS(ie,:)-cumsumSISS(is,:)); %find the offset with the max cc value
            if kmax<=cushion || kmax>nshift+mshift-cushion+1 %ignore, if the max cc value lies within the cushion
                SISSxamalg(is:ie)=0;
                rootSISSxamalg(is:ie)=0;
            else
                SISSxamalg(is:ie)=SISSx(is:ie,kmax); %amalg gets the contribution from the strongest 15(e.g.)-sec sub-window
                rootSISSxamalg(is:ie)=rootSISSx(is:ie,kmax); %use the offset from the velocity seismogram (squared).
            end
        end %now add the remainder
        is=ie+1;
        ie=winlen;
        [~,kmax]=max(cumsumSISS(ie,:)-cumsumSISS(is,:)); %find the offset with the max cc value.  Looks like not normalized.
        if kmax<=cushion || kmax>nshift+mshift-cushion+1
            SISSxamalg(is:ie)=0;
            rootSISSxamalg(is:ie)=0;
        else
            SISSxamalg(is:ie)=SISSx(is:ie,kmax); %and again
            rootSISSxamalg(is:ie)=rootSISSx(is:ie,kmax); %amalg gets the contribution from the strongest 15(e.g.)-sec sub-window
        end
        windowssquared(totlen+1:totlen+winlen)=cumdenom; %0.5*(SILB2+SSIB2); %appending the squares of the seismograms for today's saw-tooth record
        windowssummed(totlen+1:totlen+winlen)=cumamp; %0.5*(SILB2+SSIB2); %appending the squares of the seismograms for today's saw-tooth record
        windowscross(totlen+1:totlen+winlen)=cumsum(SISSxamalg); %append this cc window to the saw-tooth record for today
        windowscohamp(totlen+1:totlen+winlen)=cumsum(rootSISSxamalg); %sum the coherent amplitude
        totlen=totlen+winlen; %update today's total length
        totcross=totcross+windowscross(end); %keep a running tally for today by adding the cumsum of each window as it is done
        totsum=totsum+windowssquared(end); %keep a running tally for today by adding the cumsum of each window as it is done
        totamp=totamp+cumamp(end); %keep a running tally for today by adding the cumsum of each window as it is done
        totcohamp=totcohamp+windowscohamp(end);
    end
    tottotcross=tottotcross+totcross; %By this time totcross is the daily contribution.  tottotcross keeps a tally over all days
    tottotsum=tottotsum+totsum; %By this time totsum is the daily contribution.  tottotsum keeps a tally over all days
    tottotcohamp=tottotcohamp+totcohamp; %By this time totsum is the daily contribution.  tottotsum keeps a tally over all days
    tottotamp=tottotamp+totamp; %By this time totsum is the daily contribution.  tottotsum keeps a tally over all days
    if mod(iday,3)==1
        h=figure('Position',[wid/3 1 1.5*wid hite/2]); %center
        subplot(3,1,1,'align')
    elseif mod(iday,3)==2
        subplot(3,1,2,'align')
    else
        subplot(3,1,3,'align')
    end
    hold on
    plot(0.4*windowssquared,'k') %The daily report
    plot(0.4*0.1*windowssummed,'k--') %The daily report
    %plot(offsetstart:0.025:86400-offsetstart,cumsquared/cumsquared(end),'k')
    %hrf = plotreflinesr(gca,rubouts,'x','k');
    %hrf = plotreflinesr(gca,rubins,'x','r');
    plot(windowscross,'r') %The daily report
    plot(0.1*windowscohamp,'r--') %The daily report
    %plot(offsetstart:0.025:86400-offsetstart,cumcross/cumcross(end),'r')
    %hrf = plotreflinesr(gca,detects(:,1),'x','g');
    %hrf = plotreflinesr(gca,detects(:,2),'x','b');
    title([YEAR,'.',JDAY,'  ',num2str(lo),'-',num2str(hi),'Hz, npass',int2str(npa),'  jiggering ',int2str(-nshift+cushion),':', ...
        int2str(mshift-cushion),'  cushion = ',int2str(cushion),'  shortwin =',int2str(shortwin/sps)])
    yma=0.6*max(max(windowssquared),0.1*max(windowssummed));
    ymi=-0.15*yma;
    text(0.05*length(windowssquared),0.8*ymi,num2str(totsum))
    text(0.2*length(windowscross),0.8*ymi,num2str(totcross),'color','r')
    text(0.75*length(windowssquared),0.8*ymi,num2str(totamp))
    text(0.9*length(windowscross),0.8*ymi,num2str(totcohamp),'color','r')
    ylim([ymi yma]);
    box on
    if mod(iday,3)==0 || iday==ndays
        set(h,'PaperPosition',[0.25 0.25 8 10.5])
        print('-depsc',['BFIGS/NEW/',YEAR,'.',JDAY,'_',num2str(lo),'-',int2str(hi),'Hz_npass',int2str(npa), ...
            'jig',int2str(-nshift+cushion),':',int2str(mshift-cushion),'cush',int2str(cushion),'short',num2str(shortwin/sps),'.eps',])
    end
    clear windowssquared
    clear windowscross
    clear windowssummed
    clear windowscohamp
end
tottotcross
tottotsum
tottotcohamp
tottotamp
