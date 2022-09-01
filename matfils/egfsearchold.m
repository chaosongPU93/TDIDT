%Sits at one spot (rotations; delays) and looks for the cross-correlation.
t=cputime; tt=cputime;
%close all
format short e
%Rotations & offsets:
%timoffrot=[227 14 18 28 +058 -074  90 85 50]; %Quiet day
%timoffrot=[227 20 08 20 +028 +004 100 380 150];
%timoffrot=[246 00 61 00 +060 -085  85  65  45];  %Near beginning; %isolated; 100%
%timoffrot=[246 04 21 16 +062 -084  90  70  45];
%timoffrot=[246 04 21 56 +061 -084  90  60  50]; %Also quite isolated
%timoffrot=[246 04 22 20 +060 -085  85  60  40];
timoffrot=[246 04 26 36 +059 -085  90  60  45];
%timoffrot=[246 04 26 52 +059 -085  90  60  45];
%timoffrot=[246 09 45 56 +061 -083  95  65  50];
%timoffrot=[246 11 08 52 +064 -083  90  60  50];
%timoffrot=[246 12 53 24 +062 -083  95  60  50];
%timoffrot=[246 09 47 40 +060 -084  90  60  50];
%timoffrot=[246 09 56 12 +062 -085  95  55  50]; %Also small
%timoffrot=[246 10 06 36 +061 -085  85  60  45]; %36415-36420, e.g.
%timoffrot=[246 11 11 16 +063 -084  90  60  50];
%timoffrot=[246 12 56 12 +062 -085  90  55  50];
%timoffrot=[246 17 17 48 +062 -087  80  60  40];
%timoffrot=[247 04 59 24 +033 -100  85 370  80];
%timoffrot=[247 10 15 24 +051 -078  85  60  55 ];
%timoffrot=[247 11 14 36 +056 -077  90  55  65];
%timoffrot=[247 14 18 28 +058 -074  90 85 50]; %(the original; catalog has 85  75  55)
%timoffrot=[247 14 42 28 +059 -073  85  65  65];  %Small signal!
%timoffrot=[250 05 09 24 +065 -068  80  75  45]; %Near start of a small migration.  V. messy.
%timoffrot=[250 05 23 56 +060 -072  85  65  50];
%timoffrot=[250 05 44 44 +056 -075  85  70  55]; %Possible EGFs here. FIRST LOOK!
     %timoffrot=[250 05 44 44 -104 +085  85  70  55]; %FAKE!
     %timoffrot=[250 05 44 44 +016 -035  85  70  55]; %FAKE!
%timoffrot=[250 09 13 40 +105 -031  85  50  80];
%timoffrot=[250 11 20 36 +105 -030  90  50  80]; %Big; good for looking for "missed" energy
%timoffrot=[250 11 25 00 +093 -035  85  70  90];
%timoffrot=[250 11 28 36 +090 -037  80  75  65];
%timoffrot=[250 11 30 44 +090 -036  75  60  70]; %early in a migration episode
%timoffrot=[250 11 34 20 +086 -032  80 105  65];
     %timoffrot=[250 11 34 20 -074 +128  80 105  65]; %FAKE!
%timoffrot=[250 11 34 20 +086 -032 +024  80 105  65  80];
%timoffrot=[250 11 47 00 +080 -042  75  70  65];
%timoffrot=[250 15 17 24 +104 -031  90  55  80]; %Big!
%timoffrot=[250 06 21 00 +100 -036  85  50  70]; %Also big, not quite as.jjj
%timoffrot=[253 01 50 20 +093 -024  85  75  85];
%timoffrot=[254 06 47 08 +086 +019  80 105  55]; %Dominated by last half of last panel
%timoffrot=[254 07 40 52 +086 +019  80 110  55]; %Consistent with lxtwmg %(whole 2nd panel, e.g.)
%timoffrot=[254 09 57 08 +086 +019  80 115  50]; %2nd panel quite impressive, again.
%timoffrot=[254 13 48 28 +086 +019  80 115  50];
%timoffrot=[254 13 50 04 +086 +019  80 115  50]; %BIG SPIKE, panel 2.  Coherence starts slightly earlier?
%timoffrot=[254 15 11 56 +086 +019  80 120  55]; %Messy; dominated by panel 2
%timoffrot=[254 20 29 00 +086 +019  85 115  55]; %Messy; mix of in-phase/out-of-phase end of panel 2
%timoffrot=[254 22 17 56 +086 +019  75 105  45];
%****************
  %250 11 23 56 +093 -035  80  75  95 23  5.2   4  16.6  
% 250 11 24 04 +093 -035  85  70  95 23  5.2   4  16.6  
% 250 11 24 12 +093 -035  85  70  95 23  5.2   4  16.6  
% 250 11 24 20 +093 -035  85  70  95 23  5.2   4  16.6  
% 250 11 24 28 +093 -035  85  70  95 23  5.2   4  16.6  
% 250 11 24 36 +093 -035  85  70  90 23  5.2   4  16.6  
% 250 11 24 44 +093 -035  85  70  90 23  5.2   4  16.6  
% 250 11 24 52 +093 -035  85  70  90 23  5.2   4  16.6  
% 250 11 25 00 +093 -035  85  70  90 23  5.2   4  16.6  
% 250 11 25 08 +093 -036  80  70  90 23  5.2   4  16.6  
% 250 11 25 16 +093 -036  80  70  90 23  5.2   4  16.6  
% 250 11 25 24 +093 -036  80  70  85 23  5.2   4  16.6  
% 250 11 25 32 +093 -036  80  70  85 23  5.2   4  16.6  
% 250 11 25 40 +093 -036  80  65  85 23  5.2   4  16.6  
% 250 11 25 48 +093 -036  80  70  85 23  5.2   4  16.6  
% 250 11 25 56 +093 -036  80  70  85 23  5.2   4  16.6  
% 250 11 26 04 +093 -036  85  70  85 23  5.2   4  16.6  
% 250 11 26 12 +092 -037  80  70  85 23  5.2   4  16.6  
% 250 11 26 20 +092 -037  80  70  85 23  5.2   4  16.6  
% 250 11 26 28 +092 -037  80  70  85 23  5.2   4  16.6  
% 250 11 26 36 +092 -037  80  70  85 23  5.2   4  16.6  
% 250 11 26 44 +092 -037  80  70  85 23  5.2   4  16.6  
% 250 11 26 52 +092 -037  80  70  85 23  5.2   4  16.6  
% 250 11 27 24 +091 -037  80  75  80  4  5.2   3   6.1   
% 250 11 27 32 +091 -037  80  75  75  4  5.2   3   6.1   
% 250 11 27 40 +091 -037  80  75  75  4  5.2   3   6.1   
% 250 11 27 48 +091 -037  80  75  75  4  5.2   3   6.1   
% 250 11 28 20 +091 -037  80  75  65 15  5.2   3   6.8   
% 250 11 28 28 +090 -037  75  70  65 15  5.2   3   6.8   
% 250 11 28 36 +090 -037  80  75  65 15  5.2   3   6.8   
% 250 11 28 44 +090 -037  80  75  70 15  5.2   3   6.8   
% 250 11 28 52 +090 -037  75  70  65 15  5.2   3   6.8   
% 250 11 29 00 +090 -037  80  70  65 15  5.2   3   6.8   
% 250 11 29 08 +090 -037  80  70  65 15  5.2   3   6.8   
% 250 11 29 16 +090 -037  80  70  65 15  5.2   3   6.8   
% 250 11 29 24 +090 -037  80  70  65 15  5.2   3   6.8   
% 250 11 29 32 +090 -037  80  70  70 15  5.2   3   6.8   
% 250 11 29 40 +089 -037  75  70  65 15  5.2   3   6.8   
% 250 11 29 48 +090 -036  75  65  65 15  5.2   3   6.8   
% 250 11 29 56 +090 -036  75  70  65 15  5.2   3   6.8   
% 250 11 30 04 +089 -036  75  65  65 15  5.2   3   6.8   
% 250 11 30 12 +089 -036  75  65  65 15  5.2   3   6.8   
% 250 11 30 44 +090 -036  75  60  70  1  5.3   3   4.4   
% 250 11 31 00 +090 -035  75  60  75  1  5.3   4   3.1   
% 250 11 31 32 +089 -034  70  60  70  3  5.4   4   2.8   
% 250 11 31 40 +088 -034  70  70  70  3  5.4   4   2.8   
% 250 11 31 48 +088 -033  75  80  70  3  5.4   4   2.8   
% 250 11 32 12 +087 -033  80  95  70 22  5.5   4   7.1   
% 250 11 32 20 +087 -033  75  95  70 22  5.5   4   7.1   
% 250 11 32 28 +087 -032  80 100  70 22  5.5   4   7.1   
% 250 11 32 36 +087 -032  80 105  70 22  5.5   4   7.1   
% 250 11 32 44 +087 -032  80 100  70 22  5.5   4   7.1   
% 250 11 32 52 +087 -032  80 100  70 22  5.5   4   7.1   
% 250 11 33 00 +087 -032  80 100  70 22  5.5   4   7.1   
% 250 11 33 08 +087 -032  80 100  70 22  5.5   4   7.1   
% 250 11 33 16 +087 -032  80 100  65 22  5.5   4   7.1   
% 250 11 33 24 +087 -032  80 100  65 22  5.5   4   7.1   
% 250 11 33 32 +087 -032  85 100  65 22  5.5   4   7.1   
% 250 11 33 40 +087 -032  80 100  65 22  5.5   4   7.1   
% 250 11 33 48 +087 -032  85 100  65 22  5.5   4   7.1
% 250 11 33 56 +087 -032  80  95  65 22  5.5   4   7.1
% 250 11 34 04 +087 -032  80 100  65 22  5.5   4   7.1
% 250 11 34 12 +086 -032  80 105  65 22  5.5   4   7.1
% 250 11 34 20 +086 -032  80 105  65 22  5.5   4   7.1
% 250 11 34 28 +086 -032  80 105  65 22  5.5   4   7.1
% 250 11 34 36 +086 -032  80 105  65 22  5.5   4   7.1
% 250 11 34 44 +086 -032  80 105  65 22  5.5   4   7.1
% 250 11 34 52 +086 -032  80 105  65 22  5.5   4   7.1
% 250 11 35 00 +086 -032  80 105  65 22  5.5   4   7.1
% 250 11 35 40 +085 -033  80 105  60  1  5.6   4   6.3
% 250 11 36 04 +085 -033  80 105  65 25  5.6   4   6.1
% 250 11 36 12 +085 -034  75 100  60 25  5.6   4   6.1
% 250 11 36 20 +085 -034  75  95  60 25  5.6   4   6.1
% 250 11 36 28 +085 -034  75  95  60 25  5.6   4   6.1
% 250 11 36 36 +084 -034  75 100  60 25  5.6   4   6.1
% 250 11 36 44 +084 -034  75 100  60 25  5.6   4   6.1
% 250 11 36 52 +084 -034  75 100  60 25  5.6   4   6.1
% 250 11 37 00 +084 -034  75  95  60 25  5.6   4   6.1
% 250 11 37 08 +084 -034  75  95  60 25  5.6   4   6.1
% 250 11 37 16 +084 -034  75  95  60 25  5.6   4   6.1
% 250 11 37 24 +084 -034  75  90  60 25  5.6   4   6.1
% 250 11 37 32 +084 -035  75  90  60 25  5.6   4   6.1
% 250 11 37 40 +084 -035  75  85  60 25  5.6   4   6.1
% 250 11 37 48 +084 -035  75  85  55 25  5.6   4   6.1
% 250 11 37 56 +084 -035  75  85  60 25  5.6   4   6.1
% 250 11 38 04 +084 -035  75  85  60 25  5.6   4   6.1
% 250 11 38 12 +084 -035  75  85  60 25  5.6   4   6.1
% 250 11 38 20 +084 -035  75  80  55 25  5.6   4   6.1
% 250 11 38 28 +084 -035  75  80  55 25  5.6   4   6.1
% 250 11 38 36 +084 -035  75  80  55 25  5.6   4   6.1
% 250 11 38 44 +083 -035  75  90  55 25  5.6   4   6.1
% 250 11 38 52 +083 -035  75  90  55 25  5.6   4   6.1
% 250 11 39 00 +083 -035  75  90  55 25  5.6   4   6.1
% 250 11 39 08 +083 -035  75  90  55 25  5.6   4   6.1
% 250 11 39 16 +083 -035  75  90  55 25  5.6   4   6.1
% 250 11 43 48 +081 -040  75  70  60  3  5.7 359   4.3
% 250 11 43 56 +081 -040  75  70  60  3  5.7 359   4.3
% 250 11 44 04 +081 -040  75  65  60  3  5.7 359   4.3
% 250 11 44 52 +081 -041  75  70  65 12  5.7 358   9.5
% 250 11 45 00 +081 -041  75  70  65 12  5.7 358   9.5
% 250 11 45 08 +081 -041  75  70  65 12  5.7 358   9.5
% 250 11 45 16 +081 -041  75  70  65 12  5.7 358   9.5
% 250 11 45 24 +081 -041  75  70  65 12  5.7 358   9.5
% 250 11 45 32 +081 -041  75  70  65 12  5.7 358   9.5
% 250 11 45 40 +081 -041  75  70  65 12  5.7 358   9.5
% 250 11 45 48 +081 -041  75  70  65 12  5.7 358   9.5
% 250 11 45 56 +081 -041  75  70  65 12  5.7 358   9.5
% 250 11 46 04 +081 -041  70  70  65 12  5.7 358   9.5
% 250 11 46 12 +080 -042  75  70  65 12  5.7 358   9.5
% 250 11 46 20 +080 -042  75  70  65 12  5.7 358   9.5
% 250 11 46 36 +080 -042  75  70  70 10  5.7 357   9.5
% 250 11 46 44 +080 -042  75  70  65 10  5.7 357   9.5
% 250 11 46 52 +080 -042  75  70  65 10  5.7 357   9.5
% 250 11 47 00 +080 -042  75  70  65 10  5.7 357   9.5
% 250 11 47 08 +080 -042  75  70  70 10  5.7 357   9.5
% 250 11 47 16 +080 -043  75  65  70 10  5.7 357   9.5
% 250 11 47 24 +080 -043  75  70  70 10  5.7 357   9.5
% 250 11 47 32 +080 -043  75  70  70 10  5.7 357   9.5
% 250 11 47 40 +080 -043  75  70  70 10  5.7 357   9.5
% 250 11 47 48 +080 -043  75  65  70 10  5.7 357   9.5
%***********************
jday=timoffrot(1);
JDAY=int2str(jday)
timlook=3600*timoffrot(2)+60*timoffrot(3)+timoffrot(4);
rotPGC=pi*(timoffrot(7)-90)/180;
rotSSIB=pi*(timoffrot(8)-90)/180;
rotSILB=pi*(timoffrot(9)-90)/180;
SSIBsoff=timoffrot(5);
SILBsoff=timoffrot(6);
SSIBtoff=SSIBsoff/40;
SILBtoff=SILBsoff/40;
%
% SSIBsoff=51;
% SILBsoff=-78;
% SSIBsoff=93;
% SILBsoff=-24;

% %Read Armbruster's detections:
% ArmCat=load('/data2/arubin/ARMB/BC2005/list.2005.pgsssiMAJ');
% detects=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
% m=0;
% Detects=0;
% for n=1:length(detects)
%      %if ArmCat(n,5)==SSIBsoff && ArmCat(n,6)==SILBsoff
%      if isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff+1]) 
%          m=m+1;
%          Detects(m)=detects(n);
%      end
% end
% ArmCat=load('/data2/arubin/ARMB/BC2005/list.2005.pgsssiMAJ_new');
% detects2_8=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
% m=0;
% Detects2_8=0;
% for n=1:length(detects2_8)
%      if isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff SILBsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff-1 SILBsoff+1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff-1]) ||...
%         isequal(ArmCat(n,5:6),[SSIBsoff+1 SILBsoff+1]) 
%          m=m+1;
%          Detects2_8(m)=detects2_8(n);
%      end
% end

%Read data:
prename=['2005.',JDAY,'.00.00.00.0000.CN'];
PGCEdat=[prename,'.PGC..BHE.D.SAC'];
PGCNdat=[prename,'.PGC..BHN.D.SAC'];
%PGCZdat=[prename,'.PGC..BHZ.D.SAC'];
SSIBEdat=[prename,'.SSIB..HHE.D.SAC'];
SSIBNdat=[prename,'.SSIB..HHN.D.SAC'];
SILBEdat=[prename,'.SILB..HHE.D.SAC'];
SILBNdat=[prename,'.SILB..HHN.D.SAC'];

[PGCE,HdrDataPGC,tnuPGC,pobjPGC,timsPGC]=readsac(PGCEdat,0,'l');
[PGCN,~,~,~,~]=readsac(PGCNdat,0,'l');
%[PGCZ,~,~,~,~]=readsac(PGCZdat,0,'l');
[SSIBE,HdrDataSSIB,tnuSSIB,pobjSSIB,timsSSIB]=readsac(SSIBEdat,0,'l');
[SSIBN,~,~,~,~]=readsac(SSIBNdat,0,'l');
[SILBE,HdrDataSILB,tnuSILB,pobjSILB,timsSILB]=readsac(SILBEdat,0,'l');
[SILBN,~,~,~,~]=readsac(SILBNdat,0,'l');
% timsPGC(1)
% timsSSIB(1)
% timsSILB(1)
%
rdsac=cputime-t
t=cputime;
    % Truncate for quicker calculation.  100sps stations are longer for
    % later interpolation
    timwin=2*60+1;
    timstart=max(0, timlook-timwin); timend=min(86400, timlook+timwin);
    PGCE=PGCE(timstart*40+1:timend*40)-mean(PGCE(timstart*40+1:timend*40));
    PGCN=PGCN(timstart*40+1:timend*40)-mean(PGCN(timstart*40+1:timend*40)); 
    %PGCZ=PGCZ(timstart*40+1:timend*40);
    SSIBE=SSIBE(timstart*100:timend*100+1)-mean(SSIBE(timstart*100:timend*100+1));
    SSIBN=SSIBN(timstart*100:timend*100+1)-mean(SSIBN(timstart*100:timend*100+1));
    SILBE=SILBE(timstart*100:timend*100+1)-mean(SILBE(timstart*100:timend*100+1));
    SILBN=SILBN(timstart*100:timend*100+1)-mean(SILBN(timstart*100:timend*100+1));
    timsPGC=timsPGC(timstart*40+1:timend*40);
    timsSSIB=timsSSIB(timstart*100:timend*100+1);
    timsSILB=timsSILB(timstart*100:timend*100+1);
%shrten=cputime-t
%t=cputime;
%
tracelenPGC=length(PGCE);
tracelenSSIB=length(SSIBE);
tracelenSILB=length(SILBE);
    timPGCfirst=timsPGC(1);
    timSSIBfirst=timsSSIB(1); %0.01s earlier than PGC
    timSILBfirst=timsSILB(1); %0.01s earlier than PGC
    timPGClast=timsPGC(tracelenPGC);
    timSSIBlast=timsSSIB(tracelenSSIB);
    timSILBlast=timsSILB(tracelenSILB);
    startdiff=timPGCfirst-timSSIBfirst
    enddiff=timSSIBlast-timPGClast
%pause

%cosine taper before filtering:
x=(0:pi/80:pi/2-pi/80)';
% %Seems to be necessary at the start of each day for PGCE:
%     PGCE(1:80)=0.;
%     PGCE(81:120)=sin(x).*PGCE(81:120); %Only at start of day!
%
PGCE(1:40)=sin(x).*PGCE(1:40);
PGCN(1:40)=sin(x).*PGCN(1:40);
%PGCZ(1:40)=sin(x).*PGCZ(1:40);
x=flipud(x);
PGCE(tracelenPGC-39:tracelenPGC)=sin(x).*PGCE(tracelenPGC-39:tracelenPGC);
PGCN(tracelenPGC-39:tracelenPGC)=sin(x).*PGCN(tracelenPGC-39:tracelenPGC);
%PGCZ(tracelenPGC-39:tracelenPGC)=sin(x).*PGCZ(tracelenPGC-39:tracelenPGC);
x=(0:pi/200:pi/2-pi/200)';
SSIBE(1:100)=sin(x).*SSIBE(1:100);
SSIBN(1:100)=sin(x).*SSIBN(1:100);
SILBE(1:100)=sin(x).*SILBE(1:100);
SILBN(1:100)=sin(x).*SILBN(1:100);
x=flipud(x);
SSIBE(tracelenSSIB-99:tracelenSSIB)=sin(x).*SSIBE(tracelenSSIB-99:tracelenSSIB);
SSIBN(tracelenSSIB-99:tracelenSSIB)=sin(x).*SSIBN(tracelenSSIB-99:tracelenSSIB);
SILBE(tracelenSILB-99:tracelenSSIB)=sin(x).*SILBE(tracelenSILB-99:tracelenSILB);
SILBN(tracelenSILB-99:tracelenSSIB)=sin(x).*SILBN(tracelenSILB-99:tracelenSILB);

% if(rem(tracelen,2)==0)
%     SeisData(tracelen)=[];
%     tims(tracelen)=[];
%     tracelen=tracelen-1;
% end
%Looks like filtering the 40-Hz data followed by interpolation introduces 
%less high-frequency stuff that does not appear in the original filtered
%40-Hz data:
%PGCE=spline(tims,SeisData,timsx);
%[PGCE]=bandpass(PGCE,100,1.0,5.0,2,1,'butter');

%Filter data:
hi=6.;
lo=1.5;
npo=2;
npa=1;
[PGCEf]=1.6e-4*bandpass(PGCE,40,lo,hi,npo,npa,'butter');
[PGCNf]=1.6e-4*bandpass(PGCN,40,lo,hi,npo,npa,'butter');
[SSIBEf]=1.4*4.0e-3*bandpass(SSIBE,100,lo,hi,npo,npa,'butter');
[SSIBNf]=1.4*4.0e-3*bandpass(SSIBN,100,lo,hi,npo,npa,'butter');
[SILBEf]=1.2*4.0e-3*bandpass(SILBE,100,lo,hi,npo,npa,'butter');
[SILBNf]=1.2*4.0e-3*bandpass(SILBN,100,lo,hi,npo,npa,'butter');
%fltr=cputime-t
%t=cputime;
%Decimate the 100 sps data:
SSIBEfd = interp1(timsSSIB,SSIBEf,timsPGC,'linear')';
SSIBNfd = interp1(timsSSIB,SSIBNf,timsPGC,'linear')';
SILBEfd = interp1(timsSILB,SILBEf,timsPGC,'linear')';
SILBNfd = interp1(timsSILB,SILBNf,timsPGC,'linear')';
%dmate=cputime-t
%t=cputime;
% SSIBEfd(1:6)
% SSIBNfd(1:6)
% SILBEfd(1:6)
% SILBNfd(1:6)

PGC=PGCEf+1i*PGCNf;
SSIB=SSIBEfd+1i*SSIBNfd;
SILB=SILBEfd+1i*SILBNfd;
PGCrot=PGC*exp(1i*rotPGC);
SSIBrot=SSIB*exp(1i*rotSSIB);
SILBrot=SILB*exp(1i*rotSILB);
    %rtate=cputime-t
    %t=cputime;

timsSSIB=timsPGC-SSIBtoff; %No longer used, I think
timsSILB=timsPGC-SILBtoff; %No longer used, I think
    %tshft=cputime-t
    %t=cputime;

if SSIBtoff > 0
    SSIBrot(1:tracelenPGC-SSIBsoff)=SSIBrot(SSIBsoff+1:tracelenPGC);
    SSIBrot(tracelenPGC-SSIBsoff+1:tracelenPGC)=0;
else
    SSIBrot(-SSIBsoff+1:tracelenPGC)=SSIBrot(1:tracelenPGC+SSIBsoff);
    SSIBrot(1:-SSIBsoff)=0;
end
if SILBtoff > 0
    SILBrot(1:tracelenPGC-SILBsoff)=SILBrot(SILBsoff+1:tracelenPGC);
    SILBrot(tracelenPGC-SILBsoff+1:tracelenPGC)=0;
else
    SILBrot(-SILBsoff+1:tracelenPGC)=SILBrot(1:tracelenPGC+SILBsoff);
    SILBrot(1:-SILBsoff)=0;
end
    %tracesht=cputime-t
    %t=cputime;

plotstart=1; %(timlook-timstart)*40
plotend=(2*timwin)*40; %(timlook-timstart)*40
ltmax=max([real(PGCrot(plotstart:plotend)); ...
          real(SSIBrot(plotstart:plotend)); ...
          real(SILBrot(plotstart:plotend))]);
ltmin=min([real(PGCrot(plotstart:plotend)); ...
          real(SSIBrot(plotstart:plotend)); ...
          real(SILBrot(plotstart:plotend))]);

figure %Superposed traces, 150s window:
subplot(4,1,1); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([timlook-timwin timlook-timwin/2 ltmin ltmax]);
%xlim([timlook-75 timlook-37.5]);
subplot(4,1,2); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([timlook-timwin/2 timlook ltmin ltmax]);
%xlim([timlook-37.5 timlook]);
subplot(4,1,3); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([timlook timlook+timwin/2 ltmin ltmax]);
%xlim([timlook timlook+37.5]);
subplot(4,1,4); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
axis([timlook+timwin/2 timlook+timwin ltmin ltmax]);
%xlim([timlook+37.5 timlook+75]);

%Autocorrelation of stations:
realPGC=real(PGCrot);
realSSIB=real(SSIBrot);
realSILB=real(SILBrot);
PGCauto=realPGC.*realPGC;
PGC2=cumsum(PGCauto);
SSIBauto=realSSIB.*realSSIB;
SSIB2=cumsum(SSIBauto);
SILBauto=realSILB.*realSILB;
SILB2=cumsum(SILBauto);

mshift=11;
lenx=tracelenPGC-2*mshift;
PGSSx=zeros(lenx, 2*mshift+1);
PGSIx=zeros(lenx, 2*mshift+1);
SISSx=zeros(lenx, 2*mshift+1);
% First index is pointwise multiplication of traces; second is shifting offset.
for n=-mshift:mshift;
    PGSSx(:,n+mshift+1)=realPGC(1+mshift:tracelenPGC-mshift).* ...
        realSSIB(1+mshift-n:tracelenPGC-mshift-n);
    PGSIx(:,n+mshift+1)=realPGC(1+mshift:tracelenPGC-mshift).* ...
        realSILB(1+mshift-n:tracelenPGC-mshift-n);
    SISSx(:,n+mshift+1)=realSILB(1+mshift:tracelenPGC-mshift).* ...
        realSSIB(1+mshift-n:tracelenPGC-mshift-n);
end
cumsumPGSS=cumsum(PGSSx);
cumsumPGSI=cumsum(PGSIx);
cumsumSISS=cumsum(SISSx);
winlen=3*40;
winoff=40/2;
nwin=floor((lenx-winlen)/winoff);
timswin=zeros(nwin,1);
sumsPGSS=zeros(nwin,2*mshift+1);
sumsPGSI=zeros(nwin,2*mshift+1);
sumsSISS=zeros(nwin,2*mshift+1);
sumsPGC2=zeros(nwin,2*mshift+1);
sumsSSIB2=zeros(nwin,2*mshift+1);
sumsSILB2=zeros(nwin,2*mshift+1);
sumsSILB2b=zeros(nwin,2*mshift+1);
for n=1:nwin;
    istart=(n-1)*winoff+1;
    iend=(n-1)*winoff+winlen+1;
    timswin(n)=timsPGC(mshift+istart+1); %+mshift??? yes. Got rid of +winlen/2.  +1 for cumsum subtracting pt 1.
    %First index is window #; second is shifts over +-mshift.
    sumsPGSS(n,:)=cumsumPGSS(iend,:)-cumsumPGSS(istart,:); %Summed x-correlations
    sumsPGSI(n,:)=cumsumPGSI(iend,:)-cumsumPGSI(istart,:);
    sumsSISS(n,:)=cumsumSISS(iend,:)-cumsumSISS(istart,:);
    sumsPGC2(n,:)=PGC2(iend+mshift)-PGC2(istart+mshift);  %PGC doesn't shift, so do this here.  PGC2 is cumsummed. Yes, +mshift.
    sumsSILB2b(n,:)=SILB2(iend+mshift)-SILB2(istart+mshift); %Similar, for the SILB-SSIB connection.
    for m=-mshift:mshift;
        sumsSSIB2(n,m+mshift+1)=SSIB2(iend+mshift-m)-SSIB2(istart+mshift-m); %+m??? (yes).
        sumsSILB2(n,m+mshift+1)=SILB2(iend+mshift-m)-SILB2(istart+mshift-m);
    end
end
denomPGSSn=realsqrt(sumsPGC2.*sumsSSIB2);
denomPGSIn=realsqrt(sumsPGC2.*sumsSILB2);
denomSISSn=realsqrt(sumsSILB2b.*sumsSSIB2);
sumsPGSSn=sumsPGSS./denomPGSSn;
sumsPGSIn=sumsPGSI./denomPGSIn;
sumsSISSn=sumsSISS./denomSISSn;
[xcmaxPGSSn,imaxPGSS]=max(sumsPGSSn,[],2);
[xcmaxPGSIn,imaxPGSI]=max(sumsPGSIn,[],2);
[xcmaxSISSn,imaxSISS]=max(sumsSISSn,[],2);
ix=sub2ind(size(denomPGSSn),(1:nwin)',imaxPGSS);
ampPGSS=sqrt(denomPGSSn(ix)); %This makes amplitude linear rather than quadratic with counts.
ix=sub2ind(size(denomPGSIn),(1:nwin)',imaxPGSI);
ampPGSI=sqrt(denomPGSIn(ix));
ix=sub2ind(size(denomSISSn),(1:nwin)',imaxSISS);
ampSISS=sqrt(denomSISSn(ix)); 
%Center them
imaxPGSS=imaxPGSS-mshift-1;
imaxPGSI=imaxPGSI-mshift-1;
imaxSISS=imaxSISS-mshift-1;

templen=4*40;
tempoff=40/2;
tempoff=10*40;
ntemp=floor((lenx-templen)/tempoff);
%ntemp=floor((lenx-templen-10*40)/tempoff);
timstemp=zeros(nwin,1);
PGCtx=zeros(tracelenPGC-templen+1,ntemp); %no +1 so you can subtract ___(j-1)
PGCtn=zeros(tracelenPGC-templen+1,ntemp);
SSIBtx=zeros(tracelenPGC-templen+1,ntemp); %no +1 so you can subtract ___(j-1)
SSIBtn=zeros(tracelenPGC-templen+1,ntemp);
SILBtx=zeros(tracelenPGC-templen+1,ntemp); %no +1 so you can subtract ___(j-1)
SILBtn=zeros(tracelenPGC-templen+1,ntemp);
alltx=zeros(tracelenPGC-templen+1,ntemp);
alltn=zeros(tracelenPGC-templen+1,ntemp);
tottn=zeros(tracelenPGC-templen+1,ntemp);
n05=0; j05=0; f05=0;
for n=1:ntemp;
    istart=(n-1)*tempoff+mshift+1; %+mshift to coincide with ...?
    %istart=(n-1)*tempoff+mshift+1+5*40; %+mshift to coincide with ...?
    iend=istart+templen-1;
    PGCtempl=realPGC(istart:iend);
    SILBtempl=realSILB(istart:iend);
    SSIBtempl=realSSIB(istart:iend);
    PGCt2=dot(PGCtempl,PGCtempl);
    SSIBt2=dot(SSIBtempl,SSIBtempl);
    SILBt2=dot(SILBtempl,SILBtempl);
    for j=2:istart-40 %so you can subtract ___(j-1) 
        PGCtx(j,n)=dot(PGCtempl,realPGC(j:j+templen-1));
        PGCtarg2=PGC2(j+templen-1)-PGC2(j-1);
        PGCtn(j,n)=PGCtx(j,n)/sqrt(PGCt2*PGCtarg2);
        SSIBtx(j,n)=dot(SSIBtempl,realSSIB(j:j+templen-1));
        SSIBtarg2=SSIB2(j+templen-1)-SSIB2(j-1);
        SSIBtn(j,n)=SSIBtx(j,n)/sqrt(SSIBt2*SSIBtarg2);
        SILBtx(j,n)=dot(SILBtempl,realSILB(j:j+templen-1));
        SILBtarg2=SILB2(j+templen-1)-SILB2(j-1);
        SILBtn(j,n)=SILBtx(j,n)/sqrt(SILBt2*SILBtarg2);
        alltx(j,n)=PGCtx(j,n)+SSIBtx(j,n)+SILBtx(j,n);
        alltn(j,n)=alltx(j,n)/sqrt((PGCt2+SSIBt2+SILBt2)*(PGCtarg2+SSIBtarg2+SILBtarg2));
        tottn(j,n)=(PGCtn(j,n)+SSIBtn(j,n)+SILBtn(j,n))/3.0;
        %if alltn(j,n) > 0.52
        if tottn(j,n) > 0.52
            f05=f05+1
            n05=n05+1; j05=j05+1;
            tempfindt(f05,:)=timsPGC(istart:iend);
            targfindt(f05,:)=timsPGC(j:j+templen-1);
            tempfindPGC(:,f05)=PGCtempl;
            tempfindSSIB(:,f05)=SSIBtempl;
            tempfindSILB(:,f05)=SILBtempl;
            targfindPGC(:,f05)=realPGC(j:j+templen-1);
            targfindSSIB(:,f05)=realSSIB(j:j+templen-1);
            targfindSILB(:,f05)=realSILB(j:j+templen-1);
            find05(f05,:)=[n,j];
        end
%         if (max([PGCtn(j,n) SSIBtn(j,n) SILBtn(j,n)]) > 0.55) && ...
%            (min([PGCtn(j,n) SSIBtn(j,n) SILBtn(j,n)]) > 0.25)
%               [n j PGCtn(j,n) SSIBtn(j,n) SILBtn(j,n)]
%         end
    end
    for j=istart+40:tracelenPGC-templen+1 
        PGCtx(j,n)=dot(PGCtempl,realPGC(j:j+templen-1));
        PGCtarg2=PGC2(j+templen-1)-PGC2(j-1);
        PGCtn(j,n)=PGCtx(j,n)/sqrt(PGCt2*PGCtarg2);
        SSIBtx(j,n)=dot(SSIBtempl,realSSIB(j:j+templen-1));
        SSIBtarg2=SSIB2(j+templen-1)-SSIB2(j-1);
        SSIBtn(j,n)=SSIBtx(j,n)/sqrt(SSIBt2*SSIBtarg2);
        SILBtx(j,n)=dot(SILBtempl,realSILB(j:j+templen-1));
        SILBtarg2=SILB2(j+templen-1)-SILB2(j-1);
        SILBtn(j,n)=SILBtx(j,n)/sqrt(SILBt2*SILBtarg2);
        alltx(j,n)=PGCtx(j,n)+SSIBtx(j,n)+SILBtx(j,n);
        alltn(j,n)=alltx(j,n)/sqrt((PGCt2+SSIBt2+SILBt2)*(PGCtarg2+SSIBtarg2+SILBtarg2));
        tottn(j,n)=(PGCtn(j,n)+SSIBtn(j,n)+SILBtn(j,n))/3.0;
        if tottn(j,n) > 0.52
            f05=f05+1 
            n05=n05+1; j05=j05+1;
            tempfindt(f05,:)=timsPGC(istart:iend);
            targfindt(f05,:)=timsPGC(j:j+templen-1);
            tempfindPGC(:,f05)=PGCtempl;
            tempfindSSIB(:,f05)=SSIBtempl;
            tempfindSILB(:,f05)=SILBtempl;
            targfindPGC(:,f05)=realPGC(j:j+templen-1);
            targfindSSIB(:,f05)=realSSIB(j:j+templen-1);
            targfindSILB(:,f05)=realSILB(j:j+templen-1);
            find05(f05,:)=[n,j];
        end
    end
end
%lump things
find05'
out=zeros(f05,1); keeptot=0;
for i=1:f05
    if out(i) < 0.5 %if not previously kicked out
        clear keeptmp
        clear keeptmpx
        keeptot=keeptot+1; %how many will be kept in the end
        inset=1; %how many in the set to compare; at least this one
        keeptmp(inset)=i; 
        keeptmpx(inset)=tottn(find05(i,2),find05(i,1));
        for k=i+1:f05
            if isequal(find05(i,1),find05(k,1)) && find05(k,2)-find05(i,2)<=80
                inset=inset+1;
                keeptmp(inset)=k; out(k)=1;
                keeptmpx(inset)=tottn(find05(k,2),find05(k,1));
            elseif find05(k,1)-find05(i,1)<=85./tempoff && abs(find05(k,2)-find05(i,2))<=85
                inset=inset+1;
                keeptmp(inset)=k; out(k)=1;
                keeptmpx(inset)=tottn(find05(k,2),find05(k,1));
            end
        end
        keeptmp
        keeptmpx
        [maxkeep,ikeep]=max(keeptmpx) %; 
        keeptmp(ikeep) %keeptmp(ikeep) is between 1 and f05
        tempkeept(keeptot,:)=tempfindt(keeptmp(ikeep),:);
        targkeept(keeptot,:)=targfindt(keeptmp(ikeep),:);
        tempkeepPGC(:,keeptot)=tempfindPGC(:,keeptmp(ikeep));
        tempkeepSSIB(:,keeptot)=tempfindSSIB(:,keeptmp(ikeep));
        tempkeepSILB(:,keeptot)=tempfindSILB(:,keeptmp(ikeep));
        targkeepPGC(:,keeptot)=targfindPGC(:,keeptmp(ikeep));
        targkeepSSIB(:,keeptot)=targfindSSIB(:,keeptmp(ikeep));
        targkeepSILB(:,keeptot)=targfindSILB(:,keeptmp(ikeep));
    end
end
keeptot
ncheck=floor((lenx-templen)/templen); %looks for the max in a string of pt-wise cc m'ments
histmat=zeros(ncheck,ntemp);
for n=1:ncheck
    istart=(n-1)*templen+1;
    iend=n*templen;
    histmat(n,:)=max(tottn(istart:iend,:));
end
nlump=floor(ntemp/4);
hist2mat=zeros(ncheck,nlump);
for n=1:nlump
    istart=(n-1)*4+1;
    iend=n*4;
    hist2mat(:,n)=max(histmat(:,istart:iend),[],2);
end
[m,n]=size(hist2mat);
hist2=reshape(hist2mat,1,m*n);
[look,x]=hist(hist2,50);
look2=[x;look];
fid = fopen('look.txt','w');
fprintf(fid,'%8.4f  %8i\n',look2);
fclose(fid);
        
% [totsor,toti]=sort(alltn);
% [m,n]=size(totsor)
% mm10=m-10;
% mm5=m-5;
% totsor(mm5:m,:)

ampmax=max([ampPGSS; ampPGSI; ampSISS]);
figure 
subplot(4,1,1); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,ampSISS,'ko','MarkerSize',5);
plot(timswin,ampPGSS,'bo','MarkerSize',5);
plot(timswin,ampPGSI,'ro','MarkerSize',5);
plot(timswin,xcmaxSISSn*ampmax,'k--');
plot(timswin,xcmaxPGSSn*ampmax,'b--');
plot(timswin,xcmaxPGSIn*ampmax,'r--');
plot(timswin,(imaxPGSI-imaxPGSS+imaxSISS)*ampmax/5,'k+','MarkerSize',7);
plot(timswin,(imaxPGSS+0.15)*ampmax/5,'b+','MarkerSize',7);
plot(timswin,(imaxPGSI-0.15)*ampmax/5,'r+','MarkerSize',7);
axis([timlook-timwin timlook-timwin/2 -ampmax ampmax]);
subplot(4,1,2); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,ampSISS,'ko','MarkerSize',5);
plot(timswin,ampPGSS,'bo','MarkerSize',5);
plot(timswin,ampPGSI,'ro','MarkerSize',5);
plot(timswin,xcmaxSISSn*ampmax,'k--');
plot(timswin,xcmaxPGSSn*ampmax,'b--');
plot(timswin,xcmaxPGSIn*ampmax,'r--');
plot(timswin,(imaxPGSI-imaxPGSS+imaxSISS)*ampmax/5,'k+','MarkerSize',7);
plot(timswin,(imaxPGSS+0.15)*ampmax/5,'b+','MarkerSize',7);
plot(timswin,(imaxPGSI-0.15)*ampmax/5,'r+','MarkerSize',7);
axis([timlook-timwin/2 timlook -ampmax ampmax]);
subplot(4,1,3); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,ampSISS,'ko','MarkerSize',5);
plot(timswin,ampPGSS,'bo','MarkerSize',5);
plot(timswin,ampPGSI,'ro','MarkerSize',5);
plot(timswin,xcmaxSISSn*ampmax,'k--');
plot(timswin,xcmaxPGSSn*ampmax,'b--');
plot(timswin,xcmaxPGSIn*ampmax,'r--');
plot(timswin,(imaxPGSI-imaxPGSS+imaxSISS)*ampmax/5,'k+','MarkerSize',7);
plot(timswin,(imaxPGSS+0.15)*ampmax/5,'b+','MarkerSize',7);
plot(timswin,(imaxPGSI-0.15)*ampmax/5,'r+','MarkerSize',7);
axis([timlook timlook+timwin/2 -ampmax ampmax]);
subplot(4,1,4); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,ampPGSS,'ko','MarkerSize',5);
plot(timswin,ampPGSS,'bo','MarkerSize',5);
plot(timswin,ampPGSI,'ro','MarkerSize',5);
plot(timswin,xcmaxSISSn*ampmax,'k--');
plot(timswin,xcmaxPGSSn*ampmax,'b--');
plot(timswin,xcmaxPGSIn*ampmax,'r--');
plot(timswin,(imaxPGSI-imaxPGSS+imaxSISS)*ampmax/5,'k+','MarkerSize',7);
plot(timswin,(imaxPGSS+0.15)*ampmax/5,'b+','MarkerSize',7);
plot(timswin,(imaxPGSI-0.15)*ampmax/5,'r+','MarkerSize',7);
axis([timlook+timwin/2 timlook+timwin -ampmax ampmax]);

figure %Superposed traces, 150s window:
subplot(4,1,1); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
for n=1:keeptot
    plot(tempkeept(n,:),tempkeepPGC(:,n),'g');
    plot(targkeept(n,:),0.6*ltmax+((n-1)/keeptot)*0.4*ltmax,'r');
end
axis([timlook-timwin timlook-timwin/2 ltmin ltmax]);
subplot(4,1,2); 
hold on
plot(timswin,zeros(nwin,1),'k:');
% plot(timswin,ampSISS,'ko','MarkerSize',5);
% plot(timswin,ampPGSS,'bo','MarkerSize',5);
% plot(timswin,ampPGSI,'ro','MarkerSize',5);
plot(timswin,xcmaxSISSn*ampmax,'k--');
plot(timswin,xcmaxPGSSn*ampmax,'b--');
plot(timswin,xcmaxPGSIn*ampmax,'r--');
plot(timswin,(imaxPGSS+0.15)*ampmax/5,'b+','MarkerSize',7);
plot(timswin,(imaxPGSI-0.15)*ampmax/5,'r+','MarkerSize',7);
plot(timswin,(imaxPGSI-imaxPGSS+imaxSISS)*ampmax/5,'k+','MarkerSize',7);
axis([timlook-timwin timlook-timwin/2 -ampmax ampmax]);
subplot(4,1,3); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
for n=1:keeptot
    plot(tempkeept(n,:),tempkeepPGC(:,n),'g');
    plot(targkeept(n,:),0.6*ltmax+((n-1)/keeptot)*0.4*ltmax,'r');
end
axis([timlook-timwin/2 timlook ltmin ltmax]);
subplot(4,1,4); 
hold on
plot(timswin,zeros(nwin,1),'k:');
% plot(timswin,ampSISS,'ko','MarkerSize',5);
% plot(timswin,ampPGSS,'bo','MarkerSize',5);
% plot(timswin,ampPGSI,'ro','MarkerSize',5);
plot(timswin,xcmaxSISSn*ampmax,'k--');
plot(timswin,xcmaxPGSSn*ampmax,'b--');
plot(timswin,xcmaxPGSIn*ampmax,'r--');
plot(timswin,(imaxPGSS+0.15)*ampmax/5,'b+','MarkerSize',7);
plot(timswin,(imaxPGSI-0.15)*ampmax/5,'r+','MarkerSize',7);
plot(timswin,(imaxPGSI-imaxPGSS+imaxSISS)*ampmax/5,'k+','MarkerSize',7);
axis([timlook-timwin/2 timlook -ampmax ampmax]);
orient landscape
print('-dpdf',[JDAY,'_',int2str(timlook),'_3.pdf'])

figure
subplot(4,1,1); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
for n=1:keeptot
    plot(tempkeept(n,:),tempkeepPGC(:,n),'g');
    plot(targkeept(n,:),0.6*ltmax+((n-1)/keeptot)*0.4*ltmax,'r');
end
axis([timlook timlook+timwin/2 ltmin ltmax]);
subplot(4,1,2); 
hold on
plot(timswin,zeros(nwin,1),'k:');
% plot(timswin,ampSISS,'ko','MarkerSize',5);
% plot(timswin,ampPGSS,'bo','MarkerSize',5);
% plot(timswin,ampPGSI,'ro','MarkerSize',5);
plot(timswin,xcmaxSISSn*ampmax,'k--');
plot(timswin,xcmaxPGSSn*ampmax,'b--');
plot(timswin,xcmaxPGSIn*ampmax,'r--');
plot(timswin,(imaxPGSS+0.15)*ampmax/5,'b+','MarkerSize',7);
plot(timswin,(imaxPGSI-0.15)*ampmax/5,'r+','MarkerSize',7);
plot(timswin,(imaxPGSI-imaxPGSS+imaxSISS)*ampmax/5,'k+','MarkerSize',7);
axis([timlook timlook+timwin/2 -ampmax ampmax]);
subplot(4,1,3); 
hold on
plot(timsPGC,real(PGCrot),'r');
plot(timsPGC,real(SSIBrot),'b');
plot(timsPGC,real(SILBrot),'k');
for n=1:keeptot
    plot(tempkeept(n,:),tempkeepPGC(:,n),'g');
    plot(targkeept(n,:),0.6*ltmax+((n-1)/keeptot)*0.4*ltmax,'r');
end
axis([timlook+timwin/2 timlook+timwin ltmin ltmax]);
subplot(4,1,4); 
hold on
plot(timswin,zeros(nwin,1),'k:');
% plot(timswin,ampSISS,'ko','MarkerSize',5);
% plot(timswin,ampPGSS,'bo','MarkerSize',5);
% plot(timswin,ampPGSI,'ro','MarkerSize',5);
plot(timswin,xcmaxSISSn*ampmax,'k--');
plot(timswin,xcmaxPGSSn*ampmax,'b--');
plot(timswin,xcmaxPGSIn*ampmax,'r--');
plot(timswin,(imaxPGSS+0.15)*ampmax/5,'b+','MarkerSize',7);
plot(timswin,(imaxPGSI-0.15)*ampmax/5,'r+','MarkerSize',7);
plot(timswin,(imaxPGSI-imaxPGSS+imaxSISS)*ampmax/5,'k+','MarkerSize',7);
axis([timlook+timwin/2 timlook+timwin -ampmax ampmax]);
orient landscape
print('-dpdf',[JDAY,'_',int2str(timlook),'_4.pdf'])

figure
subplot(4,6,1); 
hold on
plot(targkeepPGC(:,1),'r');
plot(tempkeepPGC(:,1),'g');
ylim([ltmin ltmax]);
subplot(4,6,2); 
hold on
plot(targkeepSSIB(:,1),'b');
plot(tempkeepSSIB(:,1),'g');
ylim([ltmin ltmax]);
subplot(4,6,3); 
hold on
plot(targkeepSILB(:,1),'b');
plot(tempkeepSILB(:,1),'g');
ylim([ltmin ltmax]);
subplot(4,6,4); 
hold on
plot(targkeepPGC(:,2),'r');
plot(tempkeepPGC(:,2),'g');
ylim([ltmin ltmax]);
subplot(4,6,5); 
hold on
plot(targkeepSSIB(:,2),'b');
plot(tempkeepSSIB(:,2),'g');
ylim([ltmin ltmax]);
subplot(4,6,6); 
hold on
plot(targkeepSILB(:,2),'b');
plot(tempkeepSILB(:,2),'g');
ylim([ltmin ltmax]);
subplot(4,6,7); 
hold on
plot(targkeepPGC(:,3),'r');
plot(tempkeepPGC(:,3),'g');
ylim([ltmin ltmax]);
subplot(4,6,8); 
hold on
plot(targkeepSSIB(:,3),'b');
plot(tempkeepSSIB(:,3),'g');
ylim([ltmin ltmax]);
subplot(4,6,9); 
hold on
plot(targkeepSILB(:,3),'b');
plot(tempkeepSILB(:,3),'g');
ylim([ltmin ltmax]);
subplot(4,6,10); 
hold on
plot(targkeepPGC(:,4),'r');
plot(tempkeepPGC(:,4),'g');
ylim([ltmin ltmax]);
subplot(4,6,11); 
hold on
plot(targkeepSSIB(:,4),'b');
plot(tempkeepSSIB(:,4),'g');
ylim([ltmin ltmax]);
subplot(4,6,12); 
hold on
plot(targkeepSILB(:,4),'b');
plot(tempkeepSILB(:,4),'g');
ylim([ltmin ltmax]);
subplot(4,6,13); 
hold on
plot(targkeepPGC(:,5),'r');
plot(tempkeepPGC(:,5),'g');
ylim([ltmin ltmax]);
subplot(4,6,14); 
hold on
plot(targkeepSSIB(:,5),'b');
plot(tempkeepSSIB(:,5),'g');
ylim([ltmin ltmax]);
subplot(4,6,15); 
hold on
plot(targkeepSILB(:,5),'b');
plot(tempkeepSILB(:,5),'g');
ylim([ltmin ltmax]);
subplot(4,6,16); 
hold on
plot(targkeepPGC(:,6),'r');
plot(tempkeepPGC(:,6),'g');
ylim([ltmin ltmax]);
subplot(4,6,17); 
hold on
plot(targkeepSSIB(:,6),'b');
plot(tempkeepSSIB(:,6),'g');
ylim([ltmin ltmax]);
subplot(4,6,18); 
hold on
plot(targkeepSILB(:,6),'b');
plot(tempkeepSILB(:,6),'g');
ylim([ltmin ltmax]);
subplot(4,6,19); 
hold on
plot(targkeepPGC(:,7),'r');
plot(tempkeepPGC(:,7),'g');
ylim([ltmin ltmax]);
subplot(4,6,20); 
hold on
plot(targkeepSSIB(:,7),'b');
plot(tempkeepSSIB(:,7),'g');
ylim([ltmin ltmax]);
subplot(4,6,21); 
hold on
plot(targkeepSILB(:,7),'b');
plot(tempkeepSILB(:,7),'g');
ylim([ltmin ltmax]);
subplot(4,6,22); 
hold on
plot(targkeepPGC(:,8),'r');
plot(tempkeepPGC(:,8),'g');
ylim([ltmin ltmax]);
subplot(4,6,23); 
hold on
plot(targkeepSSIB(:,8),'b');
plot(tempkeepSSIB(:,8),'g');
ylim([ltmin ltmax]);
subplot(4,6,24); 
hold on
plot(targkeepSILB(:,8),'b');
plot(tempkeepSILB(:,8),'g');
ylim([ltmin ltmax]);

figure
subplot(4,6,1); 
hold on
plot(targkeepPGC(:,9),'r');
plot(tempkeepPGC(:,9),'g');
ylim([ltmin ltmax]);
subplot(4,6,2); 
hold on
plot(targkeepSSIB(:,9),'b');
plot(tempkeepSSIB(:,9),'g');
ylim([ltmin ltmax]);
subplot(4,6,3); 
hold on
plot(targkeepSILB(:,9),'b');
plot(tempkeepSILB(:,9),'g');
ylim([ltmin ltmax]);
subplot(4,6,4); 
hold on
plot(targkeepPGC(:,10),'r');
plot(tempkeepPGC(:,10),'g');
ylim([ltmin ltmax]);
subplot(4,6,5); 
hold on
plot(targkeepSSIB(:,10),'b');
plot(tempkeepSSIB(:,10),'g');
ylim([ltmin ltmax]);
subplot(4,6,6); 
hold on
plot(targkeepSILB(:,10),'b');
plot(tempkeepSILB(:,10),'g');
ylim([ltmin ltmax]);
subplot(4,6,7); 
hold on
plot(targkeepPGC(:,11),'r');
plot(tempkeepPGC(:,11),'g');
ylim([ltmin ltmax]);
subplot(4,6,8); 
hold on
plot(targkeepSSIB(:,11),'b');
plot(tempkeepSSIB(:,11),'g');
ylim([ltmin ltmax]);
subplot(4,6,9); 
hold on
plot(targkeepSILB(:,11),'b');
plot(tempkeepSILB(:,11),'g');
ylim([ltmin ltmax]);
subplot(4,6,10); 
hold on
plot(targkeepPGC(:,12),'r');
plot(tempkeepPGC(:,12),'g');
ylim([ltmin ltmax]);
subplot(4,6,11); 
hold on
plot(targkeepSSIB(:,12),'b');
plot(tempkeepSSIB(:,12),'g');
ylim([ltmin ltmax]);
subplot(4,6,12); 
hold on
plot(targkeepSILB(:,12),'b');
plot(tempkeepSILB(:,12),'g');
ylim([ltmin ltmax]);
subplot(4,6,13); 
hold on
plot(targkeepPGC(:,13),'r');
plot(tempkeepPGC(:,13),'g');
ylim([ltmin ltmax]);
subplot(4,6,14); 
hold on
plot(targkeepSSIB(:,13),'b');
plot(tempkeepSSIB(:,13),'g');
ylim([ltmin ltmax]);
subplot(4,6,15); 
hold on
plot(targkeepSILB(:,13),'b');
plot(tempkeepSILB(:,13),'g');
ylim([ltmin ltmax]);
subplot(4,6,16); 
hold on
plot(targkeepPGC(:,14),'r');
plot(tempkeepPGC(:,14),'g');
ylim([ltmin ltmax]);
subplot(4,6,17); 
hold on
plot(targkeepSSIB(:,14),'b');
plot(tempkeepSSIB(:,14),'g');
ylim([ltmin ltmax]);
subplot(4,6,18); 
hold on
plot(targkeepSILB(:,14),'b');
plot(tempkeepSILB(:,14),'g');
ylim([ltmin ltmax]);
subplot(4,6,19); 
hold on
plot(targkeepPGC(:,15),'r');
plot(tempkeepPGC(:,15),'g');
ylim([ltmin ltmax]);
subplot(4,6,20); 
hold on
plot(targkeepSSIB(:,15),'b');
plot(tempkeepSSIB(:,15),'g');
ylim([ltmin ltmax]);
subplot(4,6,21); 
hold on
plot(targkeepSILB(:,15),'b');
plot(tempkeepSILB(:,15),'g');
ylim([ltmin ltmax]);
subplot(4,6,22); 
hold on
plot(targkeepPGC(:,16),'r');
plot(tempkeepPGC(:,16),'g');
ylim([ltmin ltmax]);
subplot(4,6,23); 
hold on
plot(targkeepSSIB(:,16),'b');
plot(tempkeepSSIB(:,16),'g');
ylim([ltmin ltmax]);
subplot(4,6,24); 
hold on
plot(targkeepSILB(:,16),'b');
plot(tempkeepSILB(:,16),'g');
ylim([ltmin ltmax]);



% figure
% subplot(3,1,1); 
% hold on
% plot(timsPGC,real(PGCrot),'r');
% plot(timsPGC,real(SSIBrot),'b');
% plot(timsPGC,real(SILBrot),'k');
% axis([timlook-timwin timlook-timwin/2 ltmin ltmax]);
% subplot(3,1,2); 
% hold on
% plot(timswin,zeros(nwin,1),'k:');
% plot(timswin,ampSISS,'ko','MarkerSize',5);
% plot(timswin,ampPGSS,'bo','MarkerSize',5);
% plot(timswin,ampPGSI,'ro','MarkerSize',5);
% plot(timswin,xcmaxSISSn*ampmax,'k--');
% plot(timswin,xcmaxPGSSn*ampmax,'b--');
% plot(timswin,xcmaxPGSIn*ampmax,'r--');
% plot(timswin,(imaxPGSS+0.15)*ampmax/5,'b+','MarkerSize',7);
% plot(timswin,(imaxPGSI-0.15)*ampmax/5,'r+','MarkerSize',7);
% plot(timswin,(imaxPGSI-imaxPGSS+imaxSISS)*ampmax/5,'k+','MarkerSize',7);
% axis([timlook-timwin timlook-timwin/2 -ampmax ampmax]);
% subplot(3,1,3); 
% hold on
% plot(timsPGCtx,PGCtn,'ro','MarkerSize',2);
% plot(timsPGCtx,SSIBtn,'bo','MarkerSize',2);
% plot(timsPGCtx,SILBtn,'ko','MarkerSize',2);
% plot(timsPGCtx,tottn,'g+','MarkerSize',4);
% axis([timlook-timwin timlook-timwin/2 -1 1.01]);
% 
% figure
% subplot(3,1,1); 
% hold on
% plot(timsPGC,real(PGCrot),'r');
% plot(timsPGC,real(SSIBrot),'b');
% plot(timsPGC,real(SILBrot),'k');
% axis([timlook-timwin/2 timlook ltmin ltmax]);
% subplot(3,1,2); 
% hold on
% plot(timswin,zeros(nwin,1),'k:');
% plot(timswin,ampSISS,'ko','MarkerSize',5);
% plot(timswin,ampPGSS,'bo','MarkerSize',5);
% plot(timswin,ampPGSI,'ro','MarkerSize',5);
% plot(timswin,xcmaxSISSn*ampmax,'k--');
% plot(timswin,xcmaxPGSSn*ampmax,'b--');
% plot(timswin,xcmaxPGSIn*ampmax,'r--');
% plot(timswin,(imaxPGSS+0.15)*ampmax/5,'b+','MarkerSize',7);
% plot(timswin,(imaxPGSI-0.15)*ampmax/5,'r+','MarkerSize',7);
% plot(timswin,(imaxPGSI-imaxPGSS+imaxSISS)*ampmax/5,'k+','MarkerSize',7);
% axis([timlook-timwin/2 timlook -ampmax ampmax]);
% subplot(3,1,3); 
% hold on
% plot(timsPGCtx,PGCtn,'ro','MarkerSize',2);
% plot(timsPGCtx,SSIBtn,'bo','MarkerSize',2);
% plot(timsPGCtx,SILBtn,'ko','MarkerSize',2);
% plot(timsPGCtx,tottn,'g+','MarkerSize',4);
% axis([timlook-timwin/2 timlook -1 1.01]);
% 
% figure
% subplot(3,1,1); 
% hold on
% plot(timsPGC,real(PGCrot),'r');
% plot(timsPGC,real(SSIBrot),'b');
% plot(timsPGC,real(SILBrot),'k');
% axis([timlook timlook+timwin/2 ltmin ltmax]);
% subplot(3,1,2); 
% hold on
% plot(timswin,zeros(nwin,1),'k:');
% plot(timswin,ampSISS,'ko','MarkerSize',5);
% plot(timswin,ampPGSS,'bo','MarkerSize',5);
% plot(timswin,ampPGSI,'ro','MarkerSize',5);
% plot(timswin,xcmaxSISSn*ampmax,'k--');
% plot(timswin,xcmaxPGSSn*ampmax,'b--');
% plot(timswin,xcmaxPGSIn*ampmax,'r--');
% plot(timswin,(imaxPGSS+0.15)*ampmax/5,'b+','MarkerSize',7);
% plot(timswin,(imaxPGSI-0.15)*ampmax/5,'r+','MarkerSize',7);
% plot(timswin,(imaxPGSI-imaxPGSS+imaxSISS)*ampmax/5,'k+','MarkerSize',7);
% axis([timlook timlook+timwin/2 -ampmax ampmax]);
% subplot(3,1,3); 
% hold on
% plot(timsPGCtx,PGCtn,'ro','MarkerSize',2);
% plot(timsPGCtx,SSIBtn,'bo','MarkerSize',2);
% plot(timsPGCtx,SILBtn,'ko','MarkerSize',2);
% plot(timsPGCtx,tottn,'g+','MarkerSize',4);
% axis([timlook timlook+timwin/2 -1 1.01]);
% 
% figure
% subplot(3,1,1); 
% hold on
% plot(timsPGC,real(PGCrot),'r');
% plot(timsPGC,real(SSIBrot),'b');
% plot(timsPGC,real(SILBrot),'k');
% axis([timlook+timwin/2 timlook+timwin ltmin ltmax]);
% subplot(3,1,2); 
% hold on
% plot(timswin,zeros(nwin,1),'k:');
% plot(timswin,ampSISS,'ko','MarkerSize',5);
% plot(timswin,ampPGSS,'bo','MarkerSize',5);
% plot(timswin,ampPGSI,'ro','MarkerSize',5);
% plot(timswin,xcmaxSISSn*ampmax,'k--');
% plot(timswin,xcmaxPGSSn*ampmax,'b--');
% plot(timswin,xcmaxPGSIn*ampmax,'r--');
% plot(timswin,(imaxPGSS+0.15)*ampmax/5,'b+','MarkerSize',7);
% plot(timswin,(imaxPGSI-0.15)*ampmax/5,'r+','MarkerSize',7);
% plot(timswin,(imaxPGSI-imaxPGSS+imaxSISS)*ampmax/5,'k+','MarkerSize',7);
% axis([timlook+timwin/2 timlook+timwin -ampmax ampmax]);
% subplot(3,1,3); 
% hold on
% plot(timsPGCtx,PGCtn,'ro','MarkerSize',2);
% plot(timsPGCtx,SSIBtn,'bo','MarkerSize',2);
% plot(timsPGCtx,SILBtn,'ko','MarkerSize',2);
% plot(timsPGCtx,tottn,'g+','MarkerSize',4);
% axis([timlook+timwin/2 timlook+timwin -1 1.01]);
% 
% [totrec,toti]=sort(tottn);
% [PGCrec,PGCi]=sort(PGCtn);
% [SSIBrec,SSIBi]=sort(SSIBtn);
% [SILBrec,SILBi]=sort(SILBtn);
% [PGCrec(length(PGCrec)-30:length(PGCrec)),SSIBrec(length(PGCrec)-30:length(PGCrec)), ...
%     SILBrec(length(PGCrec)-30:length(PGCrec)),totrec(length(PGCrec)-30:length(PGCrec))]
% [PGCi(length(PGCrec)-30:length(PGCrec)),SSIBi(length(PGCrec)-30:length(PGCrec)), ...
%     SILBi(length(PGCrec)-30:length(PGCrec)),toti(length(PGCrec)-30:length(PGCrec))]


%cputime-t;
tot=cputime-tt
