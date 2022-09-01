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

     %timoffrot(1,:)=[2003 062]; %(Monday)
     %timoffrot(1,:)=[2003 063]; %(Tuesday)
     %timoffrot(1,:)=[2004 196]; %(Wednesday)
     %timoffrot(1,:)=[2004 197]; %(Thursday)
%    %timoffrot(1,:)=[2004 199]; %(Saturday)
     %timoffrot(1,:)=[2005 254]; %(Sunday)
      timoffrot(1,:)=[2005 255]; %(Monday)
     %timoffrot(1,:)=[2004 185];
     %timoffrot(1,:)=[2003 069];
     %timoffrot(1,:)=[2005 213];
     noisedays062=[069]; %No good.  Can come at wrong times of day.
     noisedays063=[069]; %No good.  Can come at wrong times of day.
     noisedays196=[187 184 183]; %[Thur Fri Mon] %No good.  Can come at wrong times of day.
     noisedays197=[187 184 183]; %[Thur Fri Mon] %No good.  Can come at wrong times of day.
     noisedays199=[185 192]; %[Thur Fri Mon] 
     noisedays254=[233 226 219]; %[Sun Sun Sun]
     noisedays255=[213 214 215 216 217 220 221]; %[Mon Tu Wed Thur Fri Mon Tu]

     swin=1800; % small window size
     ip=4;
     while 2^ip < swin
        ip=ip+1;
        npts=2^ip;
     end
tremordays=[62 63 196 197 199 254 255];
if ismember(timoffrot(1,2),tremordays)
    noise=0
else
    noise=1
end
if timoffrot(1,2)==62 || ismember(timoffrot(1,2),noisedays062)
      noisedays=noisedays062;
      wins=[22178 22275;
            22311 22473;
            22514 22540;
            22578 22597.4;
            22745 22831;
            22928 23270;
            25106 25285;
            25761 25789;
            26369 26781;
            26872 26903;
            27178 27198;
            27231 27257;
            27280 27295;
            27501 27604;
            27903 28003;
            28459 28505;
            28834 29044;
            29823 29922;
            30299 30500;
            30926 30987;
            31016 31099;
            31198 31244;
            31279 31439;
            31713 31853;
            32498 32740;
            32789 32846;
            34178 34407;
            35777 36041;
            36094 36133;
            38397 38875;
            40846 40885;
            40923 41193;
            41326 41505;
            43310 43420;
            43641 43678;
            46352 46684;
            49388 49490;
            52252 52486;
            54926 55066;
            68410 68474;
            68519 68745;
            75109 76405;
            77013 77995];
elseif timoffrot(1,2)==63 || ismember(timoffrot(1,2),noisedays063)
      noisedays=noisedays063;
      wins=[6605 7073; %2003.063
            8830 9034;
            10340 10456;
            12440 12521;
            21795 22579;
            32350 33058;
            41375 41692;
            50200 50431;
            65000 65489;
            80425 80995];
% % %    wins=[21800 22640; %2003.063
% % %          32300 32900;
% % %          32900 33080;
% % %          41382 41652;
% % %          41652 41832];
% % %    wins=[41473 41528]; %41503 41528]; %2003.063
elseif timoffrot(1,2)==196 || ismember(timoffrot(1,2),noisedays196)
      noisedays=noisedays196;
      wins=[39460 40060; %2004.196
            40360 40540;
            40940 41364;
            41955 42275;
            42910 43032;
            43610 43755;
            44220 44355;
            46260 46575;
            47080 47142;
            48360 48773;
            51080 51500; 
            56230 56624;
            56845 57315;
            62500 63148;
            64935 65270;
            70420 70684;
            73940 74277;
            78885 79124;
            81520 81768;
            85830 86199];
elseif timoffrot(1,2)==197 || ismember(timoffrot(1,2),noisedays197)
      noisedays=noisedays197;
      wins=[1361 1387; %2004.197 (postdates 30-s)
            5315 5620;
            8261 8465;
            13593 14308;
            19920 20257;
            25229 25625;
            27887 27913;
            31094 31511;
            37815 38086;
            39318 39597;
            43288 43388;
            43517 43546;
            45911 46016;
            46155 46259;
            70165 70193;
            70431 70514;
            77370 77723;
            77942 78043;
            78586 78780;
            79344 79378;
            83567 83714;
            84437 84554;
            84703 84799;
            84968 85133];
elseif timoffrot(1,2)==199 || ismember(timoffrot(1,2),noisedays199)
      noisedays=noisedays199;
      wins=[4245 4260; %2004.199
            8968 8985;
            9000 9015;
            9055 9065];
elseif timoffrot(1,2)==254 || ismember(timoffrot(1,2),noisedays254)
     noisedays=noisedays254;
       wins=[14640 15270; %2005.254 
             15740 15859; 
             18090 18170; 
             20035 20115; 
             20560 20631; 
             21155 21331; 
             24335 24525; 
             27430 27694;
             30330 31065;
             34560 34945;
             35260 35903;
             49140 50144;
             54290 54719;
             68770 69082;
             73360 73765; 
             79970 80340
             82915 83015];
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
elseif timoffrot(1,2)==255 || ismember(timoffrot(1,2),noisedays255)
     noisedays=noisedays255;
     wins=[5310 5800; %2005.255
           8064 8149;
           19006 19245;
           22984 23030;
           35688 36362;
           51033 51548;
           54210 54359;
           58935 59423;
           62310 62626;
           68205 68544;
           75545 75638];
end
if ismember(timoffrot(1,2),noisedays062)
    timoffrot(1:length(noisedays),1)=2003;
    timoffrot(:,2)=noisedays';
elseif ismember(timoffrot(1,2),noisedays063)
    timoffrot(1:length(noisedays),1)=2003;
    timoffrot(:,2)=noisedays';
elseif ismember(timoffrot(1,2),noisedays196)
    timoffrot(1:length(noisedays),1)=2004;
    timoffrot(:,2)=noisedays';
elseif ismember(timoffrot(1,2),noisedays197)
    timoffrot(1:length(noisedays),1)=2004;
    timoffrot(:,2)=noisedays';
elseif ismember(timoffrot(1,2),noisedays199)
    timoffrot(1:length(noisedays),1)=2004;
    timoffrot(:,2)=noisedays';
elseif ismember(timoffrot(1,2),noisedays254)
    timoffrot(1:length(noisedays),1)=2005;
    timoffrot(:,2)=noisedays';
elseif ismember(timoffrot(1,2),noisedays255)
    timoffrot(1:length(noisedays),1)=2005;
    timoffrot(:,2)=noisedays';
end

% % %    wins=[26045 27297; %2003.062
% % %          27916 28005;
% % %          28839 29051;
% % %          29757 29922;
% % %          38424 38789];   %wins=[26185 26275];


POLSTA=['SILB']; %SSIB bad above 20 Hz
%       'SSIB'];
%        'KLNB'];
%        'TWKB' %leave out for 2005.255 !!
%        'MGCB'];
nsta=size(POLSTA,1);
POLROTS=[0 27.5 210 20 -150];  %SILB.  Last entry is P-wave delay (for SILB, -170+20)
%           6 85 33 86;  %SSIB from Yajun
%           0 90  7 -4;  %KLNB
% 	      4 75 38 -5;  %TWKB Not sure about TWKB/MGCB (<== what does that mean ??)
%           4 70 48 -26]; %MGCB
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.; %2 and 3 are now incidence angle and back-azimuth
POLROTS(:,1)=round(POLROTS(:,1)*(100/40)); %Polaris stations at 100 sps; table has split time at 40 sps.

nd=size(timoffrot,1);
for id=1:nd
    yr=timoffrot(id,1);
    YEAR=int2str(yr);
    jday=timoffrot(id,2);
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

    aaatot=zeros(7,npts/2+1);
    bbbtot=zeros(7,npts/2+1);
    dddtot=zeros(7,npts/2+1);
    ccctot=zeros(7,npts/2+1);
    eeetot=zeros(7,npts/2+1);
%     for ista=1:size(PERMSTA,1)
%         hi=20;
%         lo=0.01;
%         sps=40;
%         STAEdat=[prename,'.',PERMSTA(ista,:),'..BHE.D.SAC']; %HHE for polstas.
%         STANdat=[prename,'.',PERMSTA(ista,:),'..BHN.D.SAC'];
%         STAZdat=[prename,'.',PERMSTA(ista,:),'..BHZ.D.SAC'];
%         [STAE,HdrDataSTA,tnuSTA,pobjSTA,timsSTA]=readsac(STAEdat,0,'l');
%         [STAN,~,~,~,~]=readsac(STANdat,0,'l');
%         [STAZ,~,~,~,~]=readsac(STAZdat,0,'l');
%         tracelen=length(STAE);
%         %cosine taper before filtering:
%         x=(0:pi/80:pi/2-pi/80)';
%         STAE(1:40)=sin(x).*STAE(1:40);
%         STAN(1:40)=sin(x).*STAN(1:40);
%         STAZ(1:40)=sin(x).*STAZ(1:40);
%         x=flipud(x);
%         STAE(tracelen-39:tracelen)=sin(x).*STAE(tracelen-39:tracelen);
%         STAN(tracelen-39:tracelen)=sin(x).*STAN(tracelen-39:tracelen);
%         STAZ(tracelen-39:tracelen)=sin(x).*STAZ(tracelen-39:tracelen);
%         %Filter data:
%         npo=2;
%         npa=2;
%         if strcmp(PERMSTA(ista,:),'PGC')
%             fact=1.6e-4;
%         elseif strcmp(PERMSTA(ista,:),'LZB')
%             fact=4.e-3;
%         end
%         [STAEf]=fact*bandpass(STAE,sps,lo,hi,npo,npa,'butter'); 
%         [STANf]=fact*bandpass(STAN,sps,lo,hi,npo,npa,'butter'); 
%         [STAZf]=fact*bandpass(STAZ,sps,lo,hi,npo,npa,'butter'); 
%         STA=STAEf+1i*STANf;
%         STAfastslow=STA*exp(-1i*PERMROTS(ista,2));
%         STAslow=real(STAfastslow);
%         STAfast=imag(STAfastslow);
%         len=length(STA);
%         STAslow(10:len-10)=STAslow(10+PERMROTS(ista,1):len-10+PERMROTS(ista,1));
% 
%         STAsplitcorrected=(STAslow+1i*STAfast)*exp(1i*PERMROTS(ista,2));
%         STAsplitcorrectedrot=STAsplitcorrected*exp(-1i*PERMROTS(ista,3));
%         STAscrot=STAsplitcorrectedrot;
% 
%         STAsoff=PERMROTS(ista,4)
%         %STAtoff=round(STAsoff/40); %really 40.  That's how it's listed.
% 
%         if STAsoff > -1
%             STAscrot(1:tracelen-STAsoff)=STAscrot(STAsoff+1:tracelen);
%             STAscrot(tracelen-STAsoff+1:tracelen)=0;
%         else
%             STAscrot(-STAsoff+1:tracelen)=STAscrot(1:tracelen+STAsoff);
%             STAscrot(1:-STAsoff)=0;
%         end
%         STAopt=real(STAscrot);
%         STAort=imag(STAscrot);
%         STAvrt=STAZf;  %should window this better with its own STAsoff.
%         STAoptfilt=bandpass(STAopt,sps,2.5,6.5,npo,npa,'butter');
%         STAoptenv=abs(hilbert(STAoptfilt));
% 
%         itot=0;
%         for n=1:size(wins,1)
%             winstart=wins(n,1)*sps;
%             winend=wins(n,2)*sps;
%             nswins=floor((winend-winstart)/(swin/2))-1;
%             for iswin=1:nswins
%                 itot=itot+1;
%                 istart=winstart+(iswin-1)*(swin/2);
%                 iend=istart-1+swin;
%                 filt='pmtm';
%                 [STAoptxx(itot,1:npts/2+1) f1] = pmtm(STAopt(istart:iend),[],[],sps);
%                  STAoptxx(itot,1:npts/2+1)=sqrt(STAoptxx(itot,1:npts/2+1));
%                 [STAortxx(itot,1:npts/2+1) f1] = pmtm(STAort(istart:iend),[],[],sps);
%                  STAortxx(itot,1:npts/2+1)=sqrt(STAortxx(itot,1:npts/2+1));
%                 [STAvrtxx(itot,1:npts/2+1) f1] = pmtm(STAvrt(istart:iend),[],[],sps);
%                  STAvrtxx(itot,1:npts/2+1)=sqrt(STAvrtxx(itot,1:npts/2+1));
%                 STAoptxx(itot,npts/2+2) = mean(STAoptenv(istart:iend));  %Last row gets estimate of envelope amplitude (optimal orientation)
%                 STAortxx(itot,npts/2+2) = mean(STAoptenv(istart:iend));  %Last row gets estimate of envelope amplitude (optimal orientation)
%                 STAvrtxx(itot,npts/2+2) = mean(STAoptenv(istart:iend));  %Last row gets estimate of envelope amplitude (optimal orientation)
%                 STAoptxx(itot,npts/2+3) = istart/sps;  %VERY last row gets estimate of envelope amplitude (optimal orientation)
%             end
%         end
%         rnperm=randperm(itot);
%         if noise == 1 %Works as long as there's a PERMSTA station first
%             for i = 1:itot
%     %             a(rnperm(i),:)=STAoptxx(i,:); %Each row is an amplitude spectrum. This is to estimate noise. Unfortunately half-overlapping windows.
%     %             b(rnperm(i),:)=STAortxx(i,:);
%     %             d(rnperm(i),:)=STAvrtxx(i,:);
%                 a(i,:)=STAoptxx(i,:); %Each row is an amplitude spectrum. This is to estimate noise. Unfortunately half-overlapping windows.
%                 b(i,:)=STAortxx(i,:);
%                 d(i,:)=STAvrtxx(i,:);
%             end
%         else
%             a=sortrows(STAoptxx,npts/2+2); %Each row is an amplitude spectrum; first row smallest amplitude; last row largest.
%             b=sortrows(STAortxx,npts/2+2);
%             d=sortrows(STAvrtxx,npts/2+2);
%         end
%         aa=cumsum(a(:,1:npts/2+1));  %cumsum cumsums the columns
%         bb=cumsum(b(:,1:npts/2+1));
%         dd=cumsum(d(:,1:npts/2+1));
% 
%         h=figure('Position',[wid/10 hite/9 1.*wid hite/1.2]);
%         inbin=round(itot/6)
%         ends=round(inbin/2);
%         %times(:,ista)=a(size(a,1):-1:size(a,1)-ends,259); %for MB
%         %How many windows in each bin? Divs gives the answer.
%         divs=[ends round(0.2*(itot-ends)) round(0.4*(itot-ends)) round(0.6*(itot-ends)) round(0.8*(itot-ends)) itot-ends itot];
%         ndivs=length(divs);
%         
%         subplot(3,2,1,'align')
%         aaa(1,:)=aa(divs(1),:)/divs(1);  %The smallest bin (aa is a cumsum)
%         for i=2:ndivs
%             aaa(i,:)=(aa(divs(i),:)-aa(divs(i-1),:))/(divs(i)-divs(i-1)); %aaa has one row for each bin.
%         end
%         %if ~strcmp(PERMSTA(ista,:),'TWKB')
%             aaatot=aaa+aaatot; %a running total as stations are added?  Not normalized, it looks like.
%         %end
%         loglog(f1,aaa)
%         hold on
%         if noise == 1
%             noise1=aaa;
%             save(['noise',PERMSTA(ista,:),'-',JDAY,'_win',int2str(swin),'.',filt,'.mat'],'noise1')
%         else
%             for iday=1:length(noisedays)
%                 if noisedays(iday) <= 9
%                     JNDAY=['00',int2str(noisedays(iday))];
%                 elseif noisedays(iday)<= 99
%                     JNDAY=['0',int2str(noisedays(iday))];
%                 else
%                     JNDAY=int2str(noisedays(iday));
%                 end
%                 load (['noise',PERMSTA(ista,:),'-',JNDAY,'_win',int2str(swin),'.',filt,'.mat']);
%                 loglog(f1,noise1,'color',[0 0 0]+0.75)
%             end
%         end
%         loglog(f1,aaa)
%         title([PERMSTA(ista,:),'  ',YEAR,'.',JDAY,'  win',num2str(swin/sps),'  ',filt,'  opt'])
%         xlim([0.1 20])
%         ylim([5e-4 0.8])
% 
%         subplot(3,2,2,'align')
%         bbb(1,:)=bb(divs(1),:)/divs(1);  %The first 64 rows (aa is a cumsum)
%         for i=2:ndivs
%             bbb(i,:)=(bb(divs(i),:)-bb(divs(i-1),:))/(divs(i)-divs(i-1));
%         end
%         %if ~strcmp(PERMSTA(ista,:),'TWKB')
%             bbbtot=bbb+bbbtot;
%         %end
%         loglog(f1,bbb)
%         hold on
%         %loglog(f1,bbbnoise,'k--')
%         title([PERMSTA(ista,:),'  ',YEAR,'.',JDAY,'  win',num2str(swin/sps),'  ',filt,'  ort'])
%         xlim([0.1 20])
%         ylim([5e-4 0.8])
% 
%         subplot(3,2,3,'align')
%         ddd(1,:)=dd(divs(1),:)/divs(1);  %The first 64 rows (aa is a cumsum)
%         for i=2:ndivs
%             ddd(i,:)=(dd(divs(i),:)-dd(divs(i-1),:))/(divs(i)-divs(i-1));
%         end
%             dddtot=ddd+dddtot;
%         loglog(f1,ddd)
%         hold on
%         title([PERMSTA(ista,:),'  ',YEAR,'.',JDAY,'  win',num2str(swin/sps),'  ',filt,'  vert'])
%         xlim([0.1 20])
%         ylim([5e-4 0.8])
% 
%         subplot(3,2,4,'align')
%         ccc=aaa./bbb;
%         %if ~strcmp(PERMSTA(ista,:),'TWKB')
%             ccctot=ccc+ccctot;
%         %end
%         loglog(f1,ccc)
%         hold on
%         title([PERMSTA(ista,:),'  ',YEAR,'.',JDAY,'  win',num2str(swin/sps),'  ',filt,'  opt/ort'])
%         xlim([0.1 20])
%         ylim([0.2 5])
%         loglog([1.e-5 1.e5],[1 1],'k--')
%         set(h,'PaperPosition',[0.25 0.25 8 10])
%         print(h,'-depsc',['2018/spectr',PERMSTA(ista,:),YEAR,'.',JDAY,'_win',int2str(swin),'.',filt,'.eps'])
%         clear STAoptxx
%         clear STAortxx
%         clear STAvrtxx
%     end

    for ista=1:size(POLSTA,1)
        hi=48;
        lo=0.01;
        sps=100;
        STAEdat=[prename,'.',POLSTA(ista,:),'..HHE.D.SAC']; %BHE for permstas.
        STANdat=[prename,'.',POLSTA(ista,:),'..HHN.D.SAC'];
        STAZdat=[prename,'.',POLSTA(ista,:),'..HHZ.D.SAC'];
        [STAE,HdrDataSTA,tnuSTA,pobjSTA,timsSTA]=readsac(STAEdat,0,'l');
        [STAN,~,~,~,~]=readsac(STANdat,0,'l');
        [STAZ,~,~,~,~]=readsac(STAZdat,0,'l');
        tracelen=length(STAE);
        %cosine taper before filtering:
        x=(0:pi/200:pi/2-pi/200)';
        STAE(1:100)=sin(x).*STAE(1:100);
        STAN(1:100)=sin(x).*STAN(1:100);
        STAZ(1:100)=sin(x).*STAZ(1:100);
        x=flipud(x);
        STAE(tracelen-99:tracelen)=sin(x).*STAE(tracelen-99:tracelen);
        STAN(tracelen-99:tracelen)=sin(x).*STAN(tracelen-99:tracelen);
        STAZ(tracelen-99:tracelen)=sin(x).*STAZ(tracelen-99:tracelen);
        %Filter data:
        npo=2;
        npa=2;
        if yr==2003 && jday<213
            [STAEf]=20.0e-3*bandpass(STAE,sps,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
            [STANf]=20.0e-3*bandpass(STAN,sps,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
            [STAZf]=20.0e-3*bandpass(STAZ,sps,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
        else
            [STAEf]=4.0e-3*bandpass(STAE,sps,lo,hi,npo,npa,'butter'); 
            [STANf]=4.0e-3*bandpass(STAN,sps,lo,hi,npo,npa,'butter'); 
            [STAZf]=4.0e-3*bandpass(STAZ,sps,lo,hi,npo,npa,'butter'); 
        end
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
        len=tracelen;

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
        STApoff=round(POLROTS(ista,5)*(100/40))
        L(-STApoff+1:len)=L(1:len+STApoff); %For P, just assume offset < 0
        L(1:-STApoff)=0;
        STAopt=Q;
        STAort=T;
        STAvrt=L;  
        clear Q
        clear T
        clear L
        STAoptfilt=bandpass(STAopt,sps,2.5,6.5,npo,npa,'butter');
        STAoptenv=abs(hilbert(STAoptfilt));

        itot=0;
        for n=1:size(wins,1)
            winstart=wins(n,1)*sps;
            winend=wins(n,2)*sps;
            nswins=floor((winend-winstart)/(swin/2))-1;
            for iswin=1:nswins
                itot=itot+1;
                istart=winstart+(iswin-1)*(swin/2);
                iend=istart-1+swin;
                filt='pmtm';
                [STAoptxx(itot,1:npts/2+1) f1] = pmtm(STAopt(istart:iend),[],[],sps); %periodogram or multi-taper (pmtm) etc.
                 STAoptxx(itot,1:npts/2+1)=sqrt(STAoptxx(itot,1:npts/2+1));
                [STAortxx(itot,1:npts/2+1) f1] = pmtm(STAort(istart:iend),[],[],sps);
                 STAortxx(itot,1:npts/2+1)=sqrt(STAortxx(itot,1:npts/2+1));
                [STAvrtxx(itot,1:npts/2+1) f1] = pmtm(STAvrt(istart:iend),[],[],sps);
                 STAvrtxx(itot,1:npts/2+1)=sqrt(STAvrtxx(itot,1:npts/2+1));
                STAoptxx(itot,npts/2+2) = mean(STAoptenv(istart:iend));  %Last row gets estimate of envelope amplitude (optimal orientation)
                STAortxx(itot,npts/2+2) = mean(STAoptenv(istart:iend));  %Last row gets estimate of envelope amplitude (optimal orientation)
                STAvrtxx(itot,npts/2+2) = mean(STAoptenv(istart:iend));  %Last row gets estimate of envelope amplitude (optimal orientation)
                STAoptxx(itot,npts/2+3) = istart/sps;  %VERY last row gets estimate of envelope amplitude (optimal orientation)
                %plot(log10(f1),log10(STAoptxx))
                %loglog(f1,STAoptxx)
                %drawnow
                %[STAortxx f1] = periodogram(STAort(istart:iend),[],[],sps);
                %loglog(f1,STAortxx,'r')
                %loglog(f1,STAoptxx./STAortxx,'g')
                %loglog(xero(:,1),xero(:,2),'k')
                %loglog(xero(:,1),10*xero(:,2),'k')
                %loglog(xero(:,1),1.e-4*xero(:,2),'k')

            end
    %         title([POLSTA(ista,:),'  ',YEAR,'.',JDAY,'  ',int2str(wins(n,1)),'-',int2str(wins(n,2)),' sec'])
    %         set(h,'PaperPosition',[0.25 5 8 5])
    %         print(h,'-depsc',['spectr',POLSTA(ista,:),YEAR,'.',JDAY,'_',int2str(wins(n,1)),'-',int2str(wins(n,2)),'.eps'])
        end
        if noise == 1 
            for i = 1:itot
                a(i,:)=STAoptxx(i,:); %Each row is an amplitude spectrum. This is to estimate noise. Unfortunately half-overlapping windows.
                b(i,:)=STAortxx(i,:);
                d(i,:)=STAvrtxx(i,:);
            end
        else
            a=sortrows(STAoptxx,npts/2+2); %Each row is an amplitude spectrum; first row smallest amplitude; last row largest.
            b=sortrows(STAortxx,npts/2+2);
            d=sortrows(STAvrtxx,npts/2+2);
        end
        aa=cumsum(a(:,1:npts/2+1));  %cumsum cumsums the columns
        bb=cumsum(b(:,1:npts/2+1));
        dd=cumsum(d(:,1:npts/2+1));

        h=figure('Position',[wid/10 hite/9 1.*wid hite/1.2]);
        subplot(3,2,1,'align')
    %     aaa(1,:)=aa(64,:)/64;  %The first 64 rows (aa is a cumsum)
    %     for i=2:floor(itot/64)
    %         aaa(i,:)=(aa(i*64,:)-aa((i-1)*64,:))/64;
    %     end
    %     aaa(i+1,:)=(aa(end,:)-aa(i*64,:))/(itot-i*64);
        inbin=round(itot/6)
        ends=round(inbin/2);
        %times(:,ista+size(PERMSTA,1))=a(size(a,1):-1:size(a,1)-ends,259);
        divs=[ends round(0.2*(itot-ends)) round(0.4*(itot-ends)) round(0.6*(itot-ends)) round(0.8*(itot-ends)) itot-ends itot];
        ndivs=length(divs);
        aaa(1,:)=aa(divs(1),:)/divs(1);  %The first 64 rows (aa is a cumsum)
        for i=2:ndivs
            aaa(i,:)=(aa(divs(i),:)-aa(divs(i-1),:))/(divs(i)-divs(i-1));
        end
        %if ~strcmp(POLSTA(ista,:),'TWKB')
            aaatot=aaa+aaatot;
        %end
        loglog(f1,aaa)
        hold on
        if noise == 1
            noise1=aaa;
            save(['noise',POLSTA(ista,:),'100-',JDAY,'_win',int2str(swin),'.rot.mat'],'noise1')
        else
            noisekeep=zeros(ndivs*length(noisedays),npts/2+1);
            for iday=1:length(noisedays)
                if noisedays(iday) <= 9
                    JNDAY=['00',int2str(noisedays(iday))];
                elseif noisedays(iday)<= 99
                    JNDAY=['0',int2str(noisedays(iday))];
                else
                    JNDAY=int2str(noisedays(iday));
                end
                load (['noise',POLSTA(ista,:),'100-',JNDAY,'_win',int2str(swin),'.rot.mat']);
                loglog(f1,noise1,'color',[0 0 0]+0.75)
                noisekeep(ndivs*(iday-1)+1:ndivs*iday,:)=noise1;
            end
            noiseprct=prctile(noisekeep,[10 50 90]);
            loglog(f1,noiseprct,'c')
        end
        loglog(f1,aaa)
        title([POLSTA(ista,:),'  ',YEAR,'.',JDAY,'  win',num2str(swin/sps),'  ',filt,'  Sh'])
        xlim([0.1 50])
        ylim([2e-4 0.8])

        subplot(3,2,2,'align')
        bbb(1,:)=bb(divs(1),:)/divs(1);  %The first 64 rows (aa is a cumsum)
        for i=2:ndivs
            bbb(i,:)=(bb(divs(i),:)-bb(divs(i-1),:))/(divs(i)-divs(i-1));
        end
        %if ~strcmp(POLSTA(ista,:),'TWKB')
            bbbtot=bbb+bbbtot;
        %end
        loglog(f1,bbb)
        hold on
        %loglog(f1,bbbnoise,'k--')
        title([POLSTA(ista,:),'  ',YEAR,'.',JDAY,'  win',num2str(swin/sps),'  ',filt,'  Sv'])
        xlim([0.1 50])
        ylim([2e-4 0.8])

        subplot(3,2,3,'align')
        ddd(1,:)=dd(divs(1),:)/divs(1);  %The first 64 rows (aa is a cumsum)
        for i=2:ndivs
            ddd(i,:)=(dd(divs(i),:)-dd(divs(i-1),:))/(divs(i)-divs(i-1));
        end
            dddtot=ddd+dddtot;
        loglog(f1,ddd)
        hold on
        if noise == 1
            noise1Z=ddd;
            save(['noise',POLSTA(ista,:),'.Z100-',JDAY,'_win',int2str(swin),'.rot.mat'],'noise1Z')
        else
           noisekeepZ=zeros(ndivs*iday,npts/2+1);
           for iday=1:length(noisedays)
                if noisedays(iday) <= 9
                    JNDAY=['00',int2str(noisedays(iday))];
                elseif noisedays(iday)<= 99
                    JNDAY=['0',int2str(noisedays(iday))];
                else
                    JNDAY=int2str(noisedays(iday));
                end
                load (['noise',POLSTA(ista,:),'.Z100-',JNDAY,'_win',int2str(swin),'.rot.mat']);
                noisekeepZ(ndivs*(iday-1)+1:ndivs*iday,:)=noise1Z;
                loglog(f1,noise1Z,'color',[0 0 0]+0.75)
            end
            noiseprctZ=prctile(noisekeepZ,[10 50 90]);
            loglog(f1,noiseprctZ,'c')
        end
        loglog(f1,ddd)
        title([POLSTA(ista,:),'  ',YEAR,'.',JDAY,'  win',num2str(swin/sps),'  P'])
        xlim([0.1 50])
        ylim([2e-4 0.8])

        subplot(3,2,4,'align')
        ccc=aaa./bbb;
        %if ~strcmp(POLSTA(ista,:),'TWKB')
            ccctot=ccc+ccctot;
        %end
        loglog(f1,ccc)
        hold on
        title([POLSTA(ista,:),'  ',YEAR,'.',JDAY,'  win',num2str(swin/sps),'   Sh/Sv'])
        xlim([0.1 50])
        ylim([2.e-1 5])
        loglog([1.e-5 1.e5],[1 1],'k--')

        if noise == 0
            subplot(3,2,5,'align')
            fff=noisekeepZ./noisekeep;
    %         %if ~strcmp(POLSTA(ista,:),'TWKB')
    %             ffftot=fff+ffftot;
    %         %end
            loglog(f1,median(noisekeepZ)./median(noisekeep),'m')
            hold on
            loglog(f1,mean(noisekeepZ)./mean(noisekeep),'c')
            title([POLSTA(ista,:),'  ',YEAR,'.',JDAY,'  win',num2str(swin/sps),'  noiseP/noiseSh'])
            xlim([0.1 50])
            ylim([2.e-1 5])
            loglog([1.e-5 1.e5],[1 1],'k--')
        end

        if noise == 0
            subplot(3,2,6,'align')
            eee=ddd./aaa;
            %if ~strcmp(POLSTA(ista,:),'TWKB')
                eeetot=eee+eeetot;
            %end
            loglog(f1,eee)
            hold on
            loglog(f1,mean(noisekeepZ)./mean(noisekeep),'c')
            title([POLSTA(ista,:),'  ',YEAR,'.',JDAY,'  win',num2str(swin/sps),'  P/Sh'])
            xlim([0.1 50])
            ylim([2.e-1 5])
            loglog([1.e-5 1.e5],[1 1],'k--')
        end

        set(h,'PaperPosition',[0.25 0.25 8 10])
        print(h,'-depsc',['2018/spectr',POLSTA(ista,:),YEAR,'.',JDAY,'_win',int2str(swin),'.rot100.eps'])
        clear STAoptxx
        clear STAortxx
        clear STAvrtxx
    end

    %h=figure('Position',[wid/10 hite/9 1.*wid hite/1.2]);
    if nsta == 5
        allsta=[POLSTA(1,:),' ',POLSTA(2,:),' ',POLSTA(3,:),' ',POLSTA(4,:),' ',POLSTA(5,:),' '];
    elseif nsta == 4
        allsta=[POLSTA(1,:),' ',POLSTA(2,:),' ',POLSTA(3,:),' ',POLSTA(4,:),' '];
    elseif nsta == 3
        allsta=[POLSTA(1,:),' ',POLSTA(2,:),' ',POLSTA(3,:),' '];
    elseif nsta == 2
        allsta=[POLSTA(1,:),' ',POLSTA(2,:),' '];
    elseif nsta == 1
        allsta=[POLSTA(1,:),' '];
    end
%     subplot(3,2,1,'align')
%     loglog(f1,aaatot/nsta)
%     hold on
%     title([allsta,'  ',YEAR,'.',JDAY,'  win',num2str(swin/sps),'  ',filt,'  opt'])
%     xlim([0.1 50])
%     ylim([2e-4 0.8])
%     subplot(3,2,2,'align')
%     loglog(f1,bbbtot/nsta)
%     title([YEAR,'.',JDAY,'  ort'])
%     xlim([0.1 50])
%     ylim([2e-4 0.8])
%     subplot(3,2,3,'align')
%     loglog(f1,dddtot/nsta)
%     title([YEAR,'.',JDAY,'  vert'])
%     xlim([0.1 50])
%     ylim([2e-4 0.8])
%     subplot(3,2,4,'align')
%     loglog(f1,ccctot/nsta)  %size(POLSTA,1)
%     hold on
%     loglog([1.e-5 1.e5],[1 1],'k--')
%     title([YEAR,'.',JDAY,'  opt/ort'])
%     xlim([0.1 50])
%     ylim([8.e-2 5])
%     set(h,'PaperPosition',[0.25 0.25 8 10])
%     print(h,'-depsc',['2018/spectr',YEAR,'.',JDAY,'_win',int2str(swin),'.',filt,'100.eps'])

    if noise 
        noise3=aaatot/nsta;
        save(['noise',int2str(nsta),'100-',JDAY,'_win',int2str(swin),'.rot.mat'],'noise3')
        noise3Z=dddtot/nsta;
        save(['noise',int2str(nsta),'.Z100-',JDAY,'_win',int2str(swin),'.rot.mat'],'noise3Z')
    else
        h=figure('Position',[wid/10 hite/9 1.*wid hite/1.2]);
        subplot(2,1,1,'align')
        loglog(f1,aaatot/nsta,'.')
        hold on
        noisekeep=zeros(ndivs*length(noisedays),npts/2+1);
        for iday=1:length(noisedays)
            if noisedays(iday) <= 9
                JNDAY=['00',int2str(noisedays(iday))];
            elseif noisedays(iday)<= 99
                JNDAY=['0',int2str(noisedays(iday))];
            else
                JNDAY=int2str(noisedays(iday));
            end
            load (['noise',int2str(nsta),'100-',JNDAY,'_win',int2str(swin),'.rot.mat']);
            noisekeep(ndivs*(iday-1)+1:ndivs*iday,:)=noise3; %3 because originally 3 stations were averaged; that number is now obsolete
            loglog(f1,noise3,'color',[0 0 0]+0.75)
        end
        noiseprctZ=prctile(noisekeepZ,[10 50 90]);
        loglog(f1,noiseprctZ,'c')
        loglog(f1,aaatot/nsta,'.')
        loglog(f1,aaatot/nsta)
        minfreq=2.4; maxfreq=4.; %minfreq=2.;
        dum=f1-minfreq;
        [~,f1min]=min(abs(dum));
        dum=f1-maxfreq;
        [~,f1max]=min(abs(dum));
        normlz=mean(aaatot(:,f1min:f1max),2);
        normtot=aaatot/nsta; %just for initialization
        for irow = 1:ndivs
            normtot(irow,:)=aaatot(irow,:)/normlz(irow);
        end
        loglog(f1,0.0035*normtot)
        text(8e-2,2e-3,int2str(itot*swin/(2*sps)),'fontsize',8);
        ylim([1e-4 1e0])
        xlim([5e-2 4e1])
        axis square
        title([allsta,'  ',YEAR,'.',JDAY,'  ',num2str(swin/sps),'-s win','  Sh'])

        subplot(2,1,2,'align')
        loglog(f1,dddtot/nsta,'.')
        hold on
        noisekeepZ=zeros(ndivs*iday,npts/2+1);
        for iday=1:length(noisedays)
            if noisedays(iday) <= 9
                JNDAY=['00',int2str(noisedays(iday))];
            elseif noisedays(iday)<= 99
                JNDAY=['0',int2str(noisedays(iday))];
            else
                JNDAY=int2str(noisedays(iday));
            end
            load (['noise',int2str(nsta),'.Z100-',JNDAY,'_win',int2str(swin),'.rot.mat']);
            noisekeepZ(ndivs*(iday-1)+1:ndivs*iday,:)=noise3Z; %3 because originally 3 stations were averaged; that number is now obsolete
            loglog(f1,noise3Z,'color',[0 0 0]+0.75)
        end
        noiseprctZ=prctile(noisekeepZ,[10 50 90]);
        loglog(f1,noiseprctZ,'c')
        loglog(f1,dddtot/nsta,'.')
        loglog(f1,dddtot/nsta)
%       minfreq=2.4; maxfreq=4.; %minfreq=2.;
%       dum=f1-minfreq;
%       [~,f1min]=min(abs(dum));
%       dum=f1-maxfreq;
%       [~,f1max]=min(abs(dum));
%       normlz=mean(dddtot(:,f1min:f1max),2);
        normtot=dddtot/nsta; %just for initialization
        for irow = 1:ndivs
            normtot(irow,:)=dddtot(irow,:)/normlz(irow);
        end
        loglog(f1,0.0035*normtot)
        text(8e-2,2e-3,int2str(itot*swin/(2*sps)),'fontsize',8);
        ylim([1e-4 1e0])
        xlim([5e-2 4e1])
        axis square
        title([allsta,'  ',YEAR,'.',JDAY,'  ',num2str(swin/sps),'-s win','  P'])

        set(h,'PaperPosition',[0.25 0.25 8 10])
        print('-depsc',['2018/spectr',YEAR,'.',JDAY,'_win',int2str(swin),'.rot.square100.eps'])
    end
    
%     h=figure('Position',[wid/10 hite/9 1.*wid hite/1.2]);
%     ngroup=7;
%     % colo=['r';'b';'g';'c';'m';'k';'r';'b';'g';'c';'m';'k'];
%     % col=colo(1:ngroup-1);
%     for i=ngroup:-1:2
%         rats=zeros(i-1,size(aaatot,2));
%         for j=1:i-1
%             rats(j,:)=aaatot(i,:)./aaatot(j,:);
%         end
%         %loglog(f1,rats/(10^(nsta-i)),col(i-1))
%         loglog(f1,rats/(6^(ngroup-i)))
%         hold on
%     end
%     title([allsta,'  ',YEAR,'.',JDAY,'  optrats'])
%     ymin=8*6^(1-ngroup); ymax=12;
%     ylim([ymin ymax])
%     xlim([4e-2 200])
%     vert=[4 ymin;
%           4 ymax];
%     loglog(vert(:,1),vert(:,2),'k--')
%     title([allsta,'  ',YEAR,'.',JDAY,'  optrats'])
%     set(gcf,'PaperPositionMode','auto')
%     print('-depsc',['2018/spectrat',YEAR,'.',JDAY,'.eps'])
end
