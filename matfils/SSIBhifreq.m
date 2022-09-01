%Reads in SILB; filters 8-16 Hz (e.g.) and rotates to search for an optimal direction
%(and splitting)?
format short e
clear all
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

   %timoffrot(1,:)=[2003 062 07 42 12 +086 +020  80 120  50]; %120 for 2-8 Hz
   %timoffrot(1,:)=[2003 063];
   %timoffrot(1,:)=[2003 064 07 42 12 +086 +020  80 115  50];
   %timoffrot(1,:)=[2003 065 07 42 12 +086 +020  80 115  50];
%    
%    timoffrot(1,:)=[2005 254];
     timoffrot(1,:)=[2005 255];
%    timoffrot(1,:)=[2005 256 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2005 257 07 42 12 +086 +020  80 115  50];
%   
%    timoffrot(1,:)=[2004 196];
%    timoffrot(1,:)=[2004 197 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2004 198 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2004 199];

% % %    wins=[26045 27297; %2003.062
% % %          27916 28005;
% % %          28839 29051;
% % %          29757 29922;
% % %          38424 38789];   %wins=[26185 26275];
% % %    wins=[21800 22640; %2003.063
% % %          32300 32900;
% % %          32900 33080;
% % %          41382 41652;
% % %          41652 41832];
% % %    wins=[41473 41528]; %41503 41528]; %2003.063
% % 
% % %   wins=[39460 40060; %2004.196
% % %         40360 40540;
% % %         40940 41364;
% % %         41955 42275;
% % %         42910 43032;
% % %         43610 43755;
% % %         44220 44355;
% % %         46260 46575;
% % %         47080 47142;
% % %         48360 48773;
% % %         51080 51500; 
% % %         56230 56624;
% % %         56845 57315;
% % %         62500 63148;
% % %         64935 65270;
% % %         70420 70684;
% % %         73940 74277;
% % %         78885 79124;
% % %         81520 81768;
% % %         85830 86199];
% % %   wins=[42960 42990]; %2004.196
% %     %wins=[51380 51399]; %2004.196
% % %     wins=[4245 4260; %2004.199
% % %           8968 8985;
% % %           9000 9015;
% % %           9055 9065];
% % 
% % %    wins=[14600 15320; %2005.254
% % %          21130 21400; 
% % %          27400 27580;
% % %          27580 27670;
% % %          27670 27730;
% % %          30300 30630;
% % %          30630 31080;
% % %          34500 34980;
% % %          35200 35710;
% % %          35710 35950;
% % %          49140 49920;
% % %          49920 50160;
% % %          54260 54770;
% % %          68760 69090;
% % %          73350 73800; 
% % %          79950 80370];
% % %      wins=[5310 5800; %2005.255
% % %            8064 8149;
% % %            19006 19245;
% % %            22984 23030;
% % %            35688 36362;
% % %            51033 51548;
% % %            54210 54359;
% % %            58935 59423;
% % %            62310 62626;
% % %            68205 68544;
% % %            75545 75638];

%      wins=[6605 7073; %2003.063
%           8830 9034;
%           10340 10456;
%           12440 12521;
%           21795 22579;
%           32350 33058;
%           41375 41692;
%           50200 50431;
%           65000 65489;
%           80425 80995];
%    wins=[14640 15270; %2005.254 
%          15740 15859; 
%          18090 18170; 
%          20035 20115; 
%          20560 20631; 
%          21155 21331; 
%          24335 24525; 
%          27430 27694;
%          30330 31065;
%          34560 34945;
%          35260 35903;
%          49140 50144;
%          54290 54719;
%          68770 69082;
%          73360 73765; 
%          79970 80340
% 	     82915 83015];
   wins=[1000 1029; %2005.255 
	     5305 5800;
	     8060 8149;
	     19000 19245;
	     22980 23030;
         35680 36362;
         51030 51544;
         54205 54359;
         58930 59423;
         62305 62626;
         68200 68544;
         75540 75638];

nd=size(timoffrot,1);

PERMSTA=['PGC'
         'LZB'];
POLSTA=['SSIB'
        'SILB'
        'KLNB'
        'TWKB'
        'MGCB'];

yr=timoffrot(nd,1);
YEAR=int2str(yr);
jday=timoffrot(nd,2);
if jday <= 9
    JDAY=['00',int2str(jday)];
elseif jday<= 99
    JDAY=['0',int2str(jday)];
else
    JDAY=int2str(jday);
end
MO=day2month(jday,yr);
IDENTIF=[YEAR,'.',JDAY]

direc=[YEAR,'/',MO,'/'];
prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];
scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite;
xero=[1.e-5 1;
      1.e5 1];
  
swin=1024; %for small window
mshift100=15;
for ista=2:2; %1:size(POLSTA,1)
    hi=16;
    lo=8;
    hi=6;
    lo=1.5;
    hi=10;
    lo=5;
    sps=100;
    STAEdat=[prename,'.',POLSTA(ista,:),'..HHE.D.SAC']; %BHE for permstas.
    STANdat=[prename,'.',POLSTA(ista,:),'..HHN.D.SAC'];
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
    npo=2;
    npa=1;
    if yr==2003 && jday<213
        [STAEf]=20.0e-3*bandpass(STAE,sps,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
        [STANf]=20.0e-3*bandpass(STAN,sps,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
    else
        [STAEf]=4.0e-3*bandpass(STAE,sps,lo,hi,npo,npa,'butter'); 
        [STANf]=4.0e-3*bandpass(STAN,sps,lo,hi,npo,npa,'butter'); 
    end
    STA=STAEf+1i*STANf;
    for degrot=0:5:175
        rotang=pi*degrot/180.;

        STArot=STA*exp(-1i*rotang);
        STAreal=real(STArot);
        STAimag=imag(STArot);

%        len=length(STA);
%        STAslow(20:len-20)=STAslow(20+POLROTS(ista,1):len-20+POLROTS(ista,1));
%        
%        STAsplitcorrected=(STAslow+1i*STAfast)*exp(1i*POLROTS(ista,2));
%        STAsplitcorrectedrot=STAsplitcorrected*exp(-1i*POLROTS(ista,3));
%        STAscrot=STAsplitcorrectedrot;
%    
%        STAsoff=POLROTS(ista,4)
%        STAtoff=round(STAsoff/40); %really 40.  That's how it's listed.
    
        itot=0;
        for n=1:size(wins,1)
            winstart=wins(n,1)*sps;
            winend=wins(n,2)*sps;
            nswins=2*floor((winend-winstart)/swin);
            for iswin=1:nswins
                itot=itot+1;
                istart=winstart+(iswin-1)*(swin/2);
                iend=istart-1+swin;
                STAxc=xcorr(STAreal(istart:iend),STAimag(istart:iend),mshift100,'coeff');
                [maxxc, imax]=max(STAxc);
                [minxc, imin]=min(STAxc);
                STAreal2=STAreal(istart:iend).*STAreal(istart:iend);
                STAimag2=STAimag(istart:iend).*STAimag(istart:iend);
                cumreal=cumsum(STAreal2);
                cumimag=cumsum(STAimag2);
                xc(itot,1:5)=[imax maxxc imin minxc cumreal(end)/cumimag(end)];
            end
        end
        figure
        subplot(3,1,1,'align')
        plot(xc(:,1)-16,'bo','markersize',2)
        hold on
        plot(xc(:,3)-16,'ro','markersize',2)
        ylim([-10 10])
        subplot(3,1,2,'align')
        plot(xc(:,2),'bo','markersize',2)
        hold on
        plot(-xc(:,4),'ro','markersize',2)
        ylim([0. 0.666])
        subplot(3,1,3,'align')
        semilogy(xc(:,4),'bo','markersize',2)
        ylim([0.1 10])
        degrot
        xcav(degrot/5+1,:)=[degrot median(xc(:,1)) median(xc(:,2)) median(xc(:,3)) median(xc(:,4)) median(xc(:,5))];
        % xcav(:,2) is imax, (:,3) is maxxc, (:,4) is imin, (:,5) is minxc, (:,6) is ratio.
    end
    figure
    subplot(3,1,1,'align')
    plot(xcav(:,1),xcav(:,2)-16,'b')
    hold on
    plot(xcav(:,1),xcav(:,4)-16,'r')
    xlim([0 180])
    ylim([-13 13])
    subplot(3,1,2,'align')
    plot(xcav(:,1),xcav(:,3),'b')
    hold on
    plot(xcav(:,1),-xcav(:,5),'r')
    xlim([0 180])
    ylim([0 0.666])
    subplot(3,1,3,'align')
    semilogy(xcav(:,1),xcav(:,6))
    xlim([0 180])
    ylim([0.1 10])
    
end %not used; different stations

%                 [STAoptxx(itot,1:swin/2+1) f1] = periodogram(STAopt(istart:iend),[],[],sps);
%                 [STAortxx(itot,1:swin/2+1) f1] = periodogram(STAort(istart:iend),[],[],sps);
%                 STAoptxx(itot,swin/2+2) = mean(STAoptenv(istart:iend));
%                 STAortxx(itot,swin/2+2) = mean(STAoptenv(istart:iend));
%                 %plot(log10(f1),log10(STAoptxx))
%                 %loglog(f1,STAoptxx)
%                 %drawnow
%                 %[STAortxx f1] = periodogram(STAort(istart:iend),[],[],sps);
%                 %loglog(f1,STAortxx,'r')
%                 %loglog(f1,STAoptxx./STAortxx,'g')
%                 %loglog(xero(:,1),xero(:,2),'k')
%                 %loglog(xero(:,1),10*xero(:,2),'k')
%                 %loglog(xero(:,1),1.e-4*xero(:,2),'k')
%                 
%             end
% %         title([POLSTA(ista,:),'  ',YEAR,'.',JDAY,'  ',int2str(wins(n,1)),'-',int2str(wins(n,2)),' sec'])
% %         set(h,'PaperPosition',[0.25 5 8 5])
% %         print(h,'-depsc',['spectr',POLSTA(ista,:),YEAR,'.',JDAY,'_',int2str(wins(n,1)),'-',int2str(wins(n,2)),'.eps'])
%     end
%     a=sortrows(STAoptxx,514);
%     aa=cumsum(a(:,1:513));
%     b=sortrows(STAortxx,514);
%     bb=cumsum(b(:,1:513));
%     h=figure('Position',[wid/10 hite/9 1.*wid hite/1.2]);
%     subplot(3,1,1,'align')
%     aaa(1,:)=aa(64,:)/64;
%     for i=2:floor(itot/64)
%         aaa(i,:)=(aa(i*64,:)-aa((i-1)*64,:))/64;
%     end
%     aaa(i+1,:)=(aa(end,:)-aa(i*64,:))/(itot-i*64);
%     loglog(f1,aaa)
%     title([POLSTA(ista,:),'  ',YEAR,'.',JDAY,'  opt'])
%     xlim([0.1 30])
%     ylim([1.e-7 1])
%     subplot(3,1,2,'align')
%     bbb(1,:)=bb(64,:)/64;
%     for i=2:floor(itot/64)
%         bbb(i,:)=(bb(i*64,:)-bb((i-1)*64,:))/64;
%     end
%     bbb(i+1,:)=(bb(end,:)-bb(i*64,:))/(itot-i*64);
%     loglog(f1,bbb)
%     title([POLSTA(ista,:),'  ',YEAR,'.',JDAY,'  ort'])
%     xlim([0.1 30])
%     ylim([1.e-7 1])
%     subplot(3,1,3,'align')
%     ccc=aaa./bbb;
%     loglog(f1,ccc)
%     hold on
%     title([POLSTA(ista,:),'  ',YEAR,'.',JDAY,'  opt/ort'])
%     xlim([0.1 30])
%     ylim([1.e-1 11])
%     loglog([1.e-5 1.e5],[1 1],'k--')
%     set(h,'PaperPosition',[0.25 0.25 8 10])
%     print(h,'-depsc',['spectr',POLSTA(ista,:),YEAR,'.',JDAY,'.eps'])
%     clear STAoptxx
%     clear STAortxx
% end
% 
% % for ista=1:size(PERMSTA,1)
% %     hi=18;
% %     lo=0.01;
% %     sps=40;
% %     STAEdat=[prename,'.',PERMSTA(ista,:),'..BHE.D.SAC']; %BHE for permstas.
% %     STANdat=[prename,'.',PERMSTA(ista,:),'..BHN.D.SAC'];
% %     [STAE,HdrDataSTA,tnuSTA,pobjSTA,timsSTA]=readsac(STAEdat,0,'l');
% %     [STAN,~,~,~,~]=readsac(STANdat,0,'l');
% %     tracelen=length(STAE);
% %     %cosine taper before filtering:
% %     x=(0:pi/80:pi/2-pi/80)';
% %     %Seems to be necessary at the start of each day for PGCE:
% %     STAE(1:80)=0.;
% %     STAE(81:120)=sin(x).*STAE(81:120); %Only at start of day!
% %     STAN(1:40)=sin(x).*STAN(1:40);
% %     x=flipud(x);
% %     STAE(tracelen-39:tracelen)=sin(x).*STAE(tracelen-39:tracelen);
% %     STAN(tracelen-39:tracelen)=sin(x).*STAN(tracelen-39:tracelen);
% %     %Filter data:
% %     npo=2;
% %     npa=1;
% %     [STAEf]=1.6e-4*bandpass(STAE,sps,lo,hi,npo,npa,'butter');
% %     [STANf]=1.6e-4*bandpass(STAN,sps,lo,hi,npo,npa,'butter');
% % 
% %     STA=STAEf+1i*STANf;
% %     if strcmp(PERMSTA(ista,:),'LZB')
% %         STA=10*STA;
% %     end
% %     STA=STAEf+1i*STANf;
% %     STAfastslow=STA*exp(-1i*POLROTS(ista,2));
% %     STAslow=real(STAfastslow);
% %     STAfast=imag(STAfastslow);
% %     len=length(STA);
% %     STAslow(20:len-20)=STAslow(20+POLROTS(ista,1):len-20+POLROTS(ista,1));
% %     
% %     STAsplitcorrected=(STAslow+1i*STAfast)*exp(1i*POLROTS(ista,2));
% %     STAsplitcorrectedrot=STAsplitcorrected*exp(-1i*POLROTS(ista,3));
% %     STAscrot=STAsplitcorrectedrot;
% % 
% %     STAsoff=PERMROTS(ista,4)
% %     STAtoff=round(STAsoff/40); 
% % 
% %     if STAtoff > -1
% %         STAscrot(1:tracelen-STAsoff)=STAscrot(STAsoff+1:tracelen);
% %         STAscrot(tracelen-STAsoff+1:tracelen)=0;
% %     else
% %         STAscrot(-STAsoff+1:tracelen)=STAscrot(1:tracelen+STAsoff);
% %         STAscrot(1:-STAsoff)=0;
% %     end
% %     STAort=real(STAscrot);
% %     STAopt=imag(STAscrot);
% % 
% %     for n=1:size(wins,1)
% %         istart=wins(n,1)*sps;
% %         iend=wins(n,2)*sps;
% %         %[Csiss,f] = mscohere(realSILB(istart:iend),realSSIB(istart:iend),[],[],[],sps);   
% %         [STAoptxx f1] = pwelch(STAopt(istart:iend),[],[],[],sps);
% %         [STAortxx f1] = pwelch(STAort(istart:iend),[],[],[],sps);
% %         h=figure('Position',[wid/10 hite/3 2.5*wid hite/1.7]);
% %         loglog(f1,STAoptxx,'b')
% %         hold on
% %         loglog(f1,STAortxx,'r')
% %         loglog(f1,STAoptxx./STAortxx,'g')
% %         %loglog(f,10.^(4*(Csiss-1)),'c')
% %         loglog(xero(:,1),xero(:,2),'k')
% %         loglog(xero(:,1),10*xero(:,2),'k')
% %         loglog(xero(:,1),1.e-4*xero(:,2),'k')
% %         xlim([0.005 18])
% %         title([PERMSTA(ista,:),'  ',YEAR,'.',JDAY,'  ',int2str(wins(n,1)),'-',int2str(wins(n,2)),' sec'])
% %         set(h,'PaperPosition',[0.25 5 8 5])
% %         print(h,'-depsc',['spectr',PERMSTA(ista,:),YEAR,'.',JDAY,'_',int2str(wins(n,1)),'-',int2str(wins(n,2)),'.eps'])
% %     end
% % end
% 
%     %
