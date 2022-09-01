% WORK IN PROGRESS!
%Sits at one spot (rotations; delays) for one day and looks for the cross-correlation.
%Now also plots aligned 4-s windows.  As of July 2012.
format short e
clear all
close all
set(0,'DefaultFigureVisible','on');
%    timoffrot(1,:)=[2005 261 03 50 44 +085 +044  75  85  45];
%    timoffrot(2,:)=[2005 261 03 50 44 +081 +053  75  60  50];
%    timoffrot(3,:)=[2005 261 03 50 44 +082 +064  80  30  65];
%    timoffrot(4,:)=[2005 261 03 50 44 +076 +072  70  10  55];
%    timoffrot(5,:)=[2005 261 03 50 44 +065 +076  65   0  25];
%    timoffrot(6,:)=[2005 261 03 50 44 +062 +052  70 340  35];

%    timoffrot(1,:)=[2005 258 03 50 44 +054 +064  70 0 35];
%    timoffrot(2,:)=[2005 259 03 50 44 +054 +064  70 0 35];
%    timoffrot(3,:)=[2005  03 50 44 +054 +064  70 0 35];
%    
%    timoffrot(4,:)=[2004 199 12 58 36 +054 +064  75 5 30];  %Should have been 70 0 35
%    timoffrot(5,:)=[2004 200 12 58 36 +054 +064  75 5 30];  %Should have been 70 0 35
%    timoffrot(6,:)=[2004 201 12 58 36 +054 +064  75 5 30];  %Should have been 70 0 35
% 
%    timoffrot(7,:)=[2004 199 03 50 44 +045 +072  70 355 20];
%    timoffrot(8,:)=[2004 200 03 50 44 +045 +072  70 355 20];
%    timoffrot(9,:)=[2004 201 03 50 44 +045 +072  70 355 20];
%    
%    timoffrot(10,:)=[2004 199 12 58 36 +036 +075  75 5 30];
%    timoffrot(11,:)=[2004 200 12 58 36 +036 +075  75 5 30];
%    timoffrot(12,:)=[2004 201 12 58 36 +036 +075  75 5 30];

%    timoffrot(1,:)=[2005 259 03 50 44 +045 +072  70 355 20];
%    timoffrot(2,:)=[2005 260 03 50 44 +045 +072  70 355 20];
%    timoffrot(3,:)=[2005 261 03 50 44 +045 +072  70 355 20];
%    
%    timoffrot(4,:)=[2005 259 12 58 36 +036 +075  75 5 30];
%    timoffrot(5,:)=[2005 260 12 58 36 +036 +075  75 5 30];
%    timoffrot(6,:)=[2005 261 12 58 36 +036 +075  75 5 30];
   
%    timoffrot(1,:)=[2005 255 07 42 12 +086 +033  80 120  45];
%    timoffrot(2,:)=[2005 256 07 42 12 +086 +033  80 120  45];
%    timoffrot(3,:)=[2005 257 07 42 12 +086 +033  80 120  45];
% 
%    timoffrot(1,:)=[2005 254 07 42 12 +074 +024  90   5  40];
%    timoffrot(2,:)=[2005 255 07 42 12 +074 +024  90   5  40];
%   timoffrot(3,:)=[2005 256 07 42 12 +074 +024  90   5  40];
%   timoffrot(4,:)=[2005 257 07 42 12 +074 +024  90   5  40];
%    timoffrot(5,:)=[2005 258 07 42 12 +074 +024  90   5  40];
%    timoffrot(6,:)=[2005 259 07 42 12 +074 +024  90   5  40];
%   timoffrot(5,:)=[2005 260 07 42 12 +074 +024  90   5  40];
% 
%    timoffrot(8,:)=[2005 255 07 42 12 +069 +039  80   5  45];
%    timoffrot(9,:)=[2005 256 07 42 12 +069 +039  80   5  45];
%    timoffrot(10,:)=[2005 257 07 42 12 +069 +039  80   5  45];
%    timoffrot(11,:)=[2005 258 07 42 12 +069 +039  80   5  45];
  

   timoffrot(1,:)=[2005 253 07 42 12 +086 +020  80 115  50];
   timoffrot(2,:)=[2005 254 07 42 12 +086 +020  80 115  50];
   timoffrot(3,:)=[2005 255 07 42 12 +086 +020  80 115  50];
   timoffrot(4,:)=[2005 256 07 42 12 +086 +020  80 115  50];
   timoffrot(5,:)=[2005 257 07 42 12 +086 +020  80 115  50];
   timoffrot(6,:)=[2005 258 07 42 12 +086 +020  80 115  50];
   timoffrot(7,:)=[2005 259 07 42 12 +086 +020  80 115  50];
   timoffrot(8,:)=[2005 260 07 42 12 +086 +020  80 115  50];
%   
%    timoffrot(7,:)=[2005 254 07 42 12 +086 +013  80 115  50];
%    timoffrot(8,:)=[2005 255 07 42 12 +086 +013  80 115  50];
%    timoffrot(9,:)=[2005 256 07 42 12 +086 +013  80 115  50];
%    timoffrot(10,:)=[2005 257 07 42 12 +086 +013  80 115  50];
% 
%    timoffrot(11,:)=[2005 254 07 42 12 +086 +027  80 120  45];
%    timoffrot(12,:)=[2005 255 07 42 12 +086 +027  80 120  45];
%    timoffrot(13,:)=[2005 256 07 42 12 +086 +027  80 120  45];
%    timoffrot(14,:)=[2005 257 07 42 12 +086 +027  80 120  45];
%    timoffrot(15,:)=[2005 258 07 42 12 +086 +027  80 120  45];
% 
   %timoffrot(3,:)=[2004 196 07 42 12 +086 +020  80 115  50];
   %timoffrot(1,:)=[2004 197 07 42 12 +086 +020  80 115  50];
%   timoffrot(4,:)=[2004 198 07 42 12 +086 +020  80 115  50];
%    timoffrot(19,:)=[2004 199 07 42 12 +086 +020  80 115  50];
%   
%    timoffrot(20,:)=[2004 195 07 42 12 +086 +013  80 115  50];
%    timoffrot(21,:)=[2004 196 07 42 12 +086 +013  80 115  50];
%    timoffrot(22,:)=[2004 197 07 42 12 +086 +013  80 115  50];
%    timoffrot(23,:)=[2004 198 07 42 12 +086 +013  80 115  50];
% 
%    timoffrot(24,:)=[2004 196 07 42 12 +086 +027  80 120  45];
%    timoffrot(25,:)=[2004 197 07 42 12 +086 +027  80 120  45];
%    timoffrot(26,:)=[2004 198 07 42 12 +086 +027  80 120  45];
%    timoffrot(27,:)=[2004 199 07 42 12 +086 +027  80 120  45];
% 
   %timoffrot(1,:)=[2003 062 07 42 12 +086 +020  80 115  50];
   %timoffrot(1,:)=[2003 063 07 42 12 +086 +020  80 115  50];
   %timoffrot(3,:)=[2003 064 07 42 12 +086 +020  80 115  50];
   %timoffrot(4,:)=[2003 065 07 42 12 +086 +020  80 115  50];
%    
%    timoffrot(32,:)=[2003 061 07 42 12 +086 +013  80 115  50];
%    timoffrot(33,:)=[2003 062 07 42 12 +086 +013  80 115  50];
%    timoffrot(34,:)=[2003 063 07 42 12 +086 +013  80 115  50];
%    timoffrot(35,:)=[2003 064 07 42 12 +086 +013  80 115  50];
% 
%    timoffrot(36,:)=[2003 062 07 42 12 +086 +027  80 120  45];
%    timoffrot(37,:)=[2003 063 07 42 12 +086 +027  80 120  45];
%    timoffrot(38,:)=[2003 064 07 42 12 +086 +027  80 120  45];
%    timoffrot(39,:)=[2003 065 07 42 12 +086 +027  80 120  45];

%    timoffrot(1,:)=[2005 254 07 42 12 +086 +020  80 115  50];
%    timoffrot(2,:)=[2005 255 07 42 12 +086 +020  80 115  50];
%    timoffrot(3,:)=[2005 256 07 42 12 +086 +020  80 115  50];
%    timoffrot(4,:)=[2005 257 07 42 12 +086 +020  80 115  50];
%   
%    timoffrot(1,:)=[2004 196 07 42 12 +086 +020  80 115  50];
%    timoffrot(2,:)=[2004 197 07 42 12 +086 +020  80 115  50];
%    timoffrot(3,:)=[2004 198 07 42 12 +086 +020  80 115  50];
%    timoffrot(4,:)=[2004 199 07 42 12 +086 +020  80 115  50];
%   
%    timoffrot(9,:)=[2003 062 07 42 12 +086 +020  80 115  50];
%    timoffrot(10,:)=[2003 063 07 42 12 +086 +020  80 115  50];
%    timoffrot(11,:)=[2003 064 07 42 12 +086 +020  80 115  50];
   
%    %%timoffrot(1,:)=[2005 253 07 42 12 +086 +020  80 120  50];
%    timoffrot(2,:)=[2005 254 07 42 12 +086 +020  80 120  50];
%    timoffrot(3,:)=[2005 255 07 42 12 +086 +020  80 120  50];
%    timoffrot(4,:)=[2005 256 07 42 12 +086 +020  80 120  50];
%    timoffrot(5,:)=[2005 257 07 42 12 +086 +020  80 120  50];
%    %%timoffrot(6,:)=[2005 258 07 42 12 +086 +020  80 120  50];
%   
%    timoffrot(7,:)=[2005 254 07 42 12 +086 +013  80 115  55];
%    timoffrot(8,:)=[2005 255 07 42 12 +086 +013  80 115  55];
%    timoffrot(9,:)=[2005 256 07 42 12 +086 +013  80 115  55];
%    timoffrot(10,:)=[2005 257 07 42 12 +086 +013  80 115  55];
% 
%    timoffrot(11,:)=[2005 254 07 42 12 +086 +027  80 120  40];
%    timoffrot(12,:)=[2005 255 07 42 12 +086 +027  80 120  40];
%    timoffrot(13,:)=[2005 256 07 42 12 +086 +027  80 120  40];
%    timoffrot(14,:)=[2005 257 07 42 12 +086 +027  80 120  40];
%    %%timoffrot(15,:)=[2005 258 07 42 12 +086 +027  80 120  40];
% 
%    timoffrot(16,:)=[2004 196 07 42 12 +086 +020  80 120  50];
%    timoffrot(17,:)=[2004 197 07 42 12 +086 +020  80 120  50];
%    timoffrot(18,:)=[2004 198 07 42 12 +086 +020  80 120  50];
%    timoffrot(19,:)=[2004 199 07 42 12 +086 +020  80 120  50];
%   
%    timoffrot(21,:)=[2004 196 07 42 12 +086 +013  80 115  55];
%    timoffrot(22,:)=[2004 197 07 42 12 +086 +013  80 115  55];
%    timoffrot(23,:)=[2004 198 07 42 12 +086 +013  80 115  55];
% 
%    timoffrot(24,:)=[2004 196 07 42 12 +086 +027  80 120  40];
%    timoffrot(25,:)=[2004 197 07 42 12 +086 +027  80 120  40];
%    timoffrot(26,:)=[2004 198 07 42 12 +086 +027  80 120  40];
%    timoffrot(27,:)=[2004 199 07 42 12 +086 +027  80 120  40];
% 
%    timoffrot(28,:)=[2003 062 07 42 12 +086 +020  80 120  50];
%    timoffrot(29,:)=[2003 063 07 42 12 +086 +020  80 120  50];
%    timoffrot(30,:)=[2003 064 07 42 12 +086 +020  80 120  50];
%    %%timoffrot(31,:)=[2003 065 07 42 12 +086 +020  80 120  50];
%    
%    timoffrot(32,:)=[2003 061 07 42 12 +086 +013  80 115  55];
%    timoffrot(33,:)=[2003 062 07 42 12 +086 +013  80 115  55];
%    timoffrot(34,:)=[2003 063 07 42 12 +086 +013  80 115  55];
%    timoffrot(35,:)=[2003 064 07 42 12 +086 +013  80 115  55];
% 
%    timoffrot(36,:)=[2003 062 07 42 12 +086 +027  80 120  40];
%    timoffrot(37,:)=[2003 063 07 42 12 +086 +027  80 120  40];
%    timoffrot(38,:)=[2003 064 07 42 12 +086 +027  80 120  40];

%    timoffrot(1,:)=[2003 269 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(2,:)=[2003 270 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(3,:)=[2003 271 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(4,:)=[2003 272 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(5,:)=[2004 117 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(6,:)=[2004 118 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(7,:)=[2004 191 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(8,:)=[2004 192 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(9,:)=[2004 193 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(10,:)=[2004 194 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(11,:)=[2004 195 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(12,:)=[2005 072 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(13,:)=[2005 073 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(14,:)=[2005 075 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(15,:)=[2005 077 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(16,:)=[2006 119 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(17,:)=[2006 120 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(18,:)=[2006 121 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(19,:)=[2005 246 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(20,:)=[2005 247 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(21,:)=[2005 248 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(22,:)=[2005 249 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(23,:)=[2004 362 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  THESE THREE SKIPPED 1ST TIME
%    timoffrot(24,:)=[2004 363 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  
%    timoffrot(25,:)=[2004 364 05 19 48 +059 -073  95  95  50];   %FOR HIGHER FREQUENCY  

%    timoffrot(1,:)=[2003 269 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(2,:)=[2003 270 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(3,:)=[2003 271 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(4,:)=[2003 272 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(5,:)=[2004 117 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(6,:)=[2004 118 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(7,:)=[2004 191 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(8,:)=[2004 192 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(9,:)=[2004 193 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(10,:)=[2004 194 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(11,:)=[2004 195 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(12,:)=[2005 072 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(13,:)=[2005 073 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(14,:)=[2005 075 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(15,:)=[2005 077 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(16,:)=[2006 119 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(17,:)=[2006 120 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(18,:)=[2006 121 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(19,:)=[2005 246 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(20,:)=[2005 247 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(21,:)=[2005 248 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(22,:)=[2005 249 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(23,:)=[2004 362 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY
%    timoffrot(24,:)=[2004 363 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  
%    timoffrot(25,:)=[2004 364 05 19 48 +059 -073  85  65  55];   %FOR LOWER FREQUENCY  

length(timoffrot(:,1))
for nd=1:length(timoffrot(:,1))
t=cputime; tt=cputime;
% loopoffmax=1.5
% xcmaxAVEnmin=0.4 %0.4 for 4s 1.5-6 Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?
% winlen=4*40;
% %winlen=128*40;
% winoff=40/1;
% %winoff=40*8;
year=timoffrot(nd,1);
YEAR=int2str(year);
jday=timoffrot(nd,2);
if jday <= 9
    JDAY=['00',int2str(jday)];
elseif jday<= 99
    JDAY=['0',int2str(jday)];
else
    JDAY=int2str(jday)
end
MO=day2month(jday,year)
timlook=3600*timoffrot(nd,3)+60*timoffrot(nd,4)+timoffrot(nd,5)
rotPGC=pi*(timoffrot(nd,8)-90)/180;
rotSSIB=pi*(timoffrot(nd,9)-90)/180;
rotSILB=pi*(timoffrot(nd,10)-90)/180;
SSIBsoff=timoffrot(nd,6);
SILBsoff=timoffrot(nd,7);
SSIBtoff=SSIBsoff/40;
SILBtoff=SILBsoff/40;
IDENTIF=[YEAR,'.',JDAY,'.',int2str(timoffrot(nd,6)),'.',int2str(timoffrot(nd,7)), ...
    '.',int2str(timoffrot(nd,8)),'.',int2str(timoffrot(nd,9)),'.',int2str(timoffrot(nd,10))]

%%Read Armbruster's detections:
%ArmCat=load(['/data2/arubin/CNDC/',YEAR,'/list.',YEAR,'.pgsssiMAJ']);
%%ArmCat=load(['/data2/arubin/CNDC/',YEAR,'/list.',YEAR,'.pgsssiMIN']);
%detects=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
%vectoff=[ArmCat(:,5)-SSIBsoff ArmCat(:,6)-SILBsoff (ArmCat(:,6)-SILBsoff)-(ArmCat(:,5)-SSIBsoff)];
%distoff=vectoff.*vectoff;
%distoff2=sum(distoff,2);
%gridoff=max(abs(vectoff),[],2); %Try this instead.  Should change w/mshift (defined later, though)
%m0=0;
%Detects=0;
%for n=1:length(detects)
%    %if distoff2(n)<=9
%    if gridoff(n)<=mshift-1
%         m0=m0+1;
%         Detects(m0)=detects(n);
%     end
%end
%%abs(xmaxPGSIn(n)-xmaxPGSSn(n)+xmaxSISSn(n))>loopoffma
%ArmCat=load(['/data2/arubin/CNDC/',YEAR,'/list.',YEAR,'.pgsssiMAJ_new']);
%detects2_8=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
%vectoff=[ArmCat(:,5)-SSIBsoff ArmCat(:,6)-SILBsoff];
%distoff=vectoff.*vectoff;
%distoff2=sum(distoff,2);
%gridoff=max(abs(vectoff),[],2); %Try this instead.  Should change w/mshift (defined later, though)
%m1=0;
%Detects2_8=0;
%for n=1:length(detects2_8)
%    %if distoff2(n)<=9
%    if gridoff(n)<=mshift-1
%         m1=m1+1;
%         Detects2_8(m1)=detects2_8(n);
%     end
%end

%Read data:
direc=[YEAR,'/',MO,'/'];
prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];
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
samplratPGC=round(1/(timsPGC(2)-timsPGC(1)))
samplratSSIB=round(1/(timsSSIB(2)-timsSSIB(1)))
samplratSILB=round(1/(timsSILB(2)-timsSILB(1)))
% timsPGC(1)
% timsSSIB(1)
% timsSILB(1)
%
rdsac=cputime-t;
t=cputime;
%     % Truncate for quicker calculation.  100sps stations are longer for
%     % later interpolation
%     timwin=4*60+1;
%     timstart=max(0, timlook-timwin); timend=min(86400, timlook+timwin);
%     PGCE=PGCE(timstart*40+1:timend*40)-mean(PGCE(timstart*40+1:timend*40));
%     PGCN=PGCN(timstart*40+1:timend*40)-mean(PGCN(timstart*40+1:timend*40)); 
%     %PGCZ=PGCZ(timstart*40+1:timend*40);
%     SSIBE=SSIBE(timstart*100:timend*100+1)-mean(SSIBE(timstart*100:timend*100+1));
%     SSIBN=SSIBN(timstart*100:timend*100+1)-mean(SSIBN(timstart*100:timend*100+1));
%     SILBE=SILBE(timstart*100:timend*100+1)-mean(SILBE(timstart*100:timend*100+1));
%     SILBN=SILBN(timstart*100:timend*100+1)-mean(SILBN(timstart*100:timend*100+1));
%     timsPGC=timsPGC(timstart*40+1:timend*40);
%     timsSSIB=timsSSIB(timstart*100:timend*100+1);
%     timsSILB=timsSILB(timstart*100:timend*100+1);
% %shrten=cputime-t
% %t=cputime;
% %
tracelenPGC=length(PGCE);
tracelenSSIB=length(SSIBE);
tracelenSILB=length(SILBE);

    timPGCfirst=timsPGC(1);
    timSSIBfirst=timsSSIB(1); %0.01s earlier than PGC
    timSILBfirst=timsSILB(1); %0.01s earlier than PGC
    timPGClast=timsPGC(tracelenPGC);
    timSSIBlast=timsSSIB(tracelenSSIB);
    timSILBlast=timsSILB(tracelenSILB);
    timstart=timPGCfirst; timend=timPGClast; 
    timwin=(timPGClast-timPGCfirst)/2;
    startdiff=timPGCfirst-timSSIBfirst
    enddiff=timSSIBlast-timPGClast

%cosine taper before filtering:
x=(0:pi/80:pi/2-pi/80)';
% %Seems to be necessary at the start of each day for PGCE:
    PGCE(1:80)=0.;
    PGCE(81:120)=sin(x).*PGCE(81:120); %Only at start of day!
%PGCE(1:40)=sin(x).*PGCE(1:40);
PGCN(1:40)=sin(x).*PGCN(1:40);
x=flipud(x);
PGCE(tracelenPGC-39:tracelenPGC)=sin(x).*PGCE(tracelenPGC-39:tracelenPGC);
PGCN(tracelenPGC-39:tracelenPGC)=sin(x).*PGCN(tracelenPGC-39:tracelenPGC);
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
hi=0.05;
lo=0.02;
resampPGC=[1 40];
resampSSIB=[1 100];
resampSILB=[1 100];
samplratPGC=samplratPGC*resampPGC(1)/resampPGC(2)
samplratSSIB=samplratSSIB*resampSSIB(1)/resampSSIB(2)
samplratSILB=samplratSILB*resampSILB(1)/resampSILB(2)
PGCEd = resample(PGCE,resampPGC(1),resampPGC(2));
PGCNd = resample(PGCN,resampPGC(1),resampPGC(2));
timsPGCd = (0:1/samplratPGC:86400-1/samplratPGC);
SSIBEd = resample(SSIBE,resampSSIB(1),resampSSIB(2));
SSIBNd = resample(SSIBN,resampSSIB(1),resampSSIB(2));
timsSSIBd = (0:1/samplratSSIB:86400-1/samplratSSIB);
SILBEd = resample(SILBE,resampSILB(1),resampSILB(2));
SILBNd = resample(SILBN,resampSILB(1),resampSILB(2));
timsSILBd = (0:1/samplratSILB:86400-1/samplratSILB);
tracelenPGC=length(PGCEd)
tracelenSSIB=length(SSIBEd)
tracelenSILB=length(SILBEd)
SSIBsoff=round(SSIBsoff*resampPGC(1)/resampPGC(2))
SILBsoff=round(SILBsoff*resampPGC(1)/resampPGC(2))

npo=2;
npa=1;
[PGCEf]=1.6e-4*bandpass(PGCEd,samplratPGC,lo,hi,npo,npa,'butter');
[PGCNf]=1.6e-4*bandpass(PGCNd,samplratPGC,lo,hi,npo,npa,'butter');
%[PGCZf]=1.6e-4*bandpass(PGCZ,40,lo,hi,npo,npa,'butter');
if year==2003 && jday<213
    [SSIBEf]=20.0e-3*bandpass(SSIBEd,samplratSSIB,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
    [SSIBNf]=20.0e-3*bandpass(SSIBNd,samplratSSIB,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
    [SILBEf]=20.0e-3*bandpass(SILBEd,samplratSSIB,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0* 
    [SILBNf]=20.0e-3*bandpass(SILBNd,samplratSSIB,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
else
    [SSIBEf]=4.0e-3*bandpass(SSIBEd,samplratSILB,lo,hi,npo,npa,'butter'); 
    [SSIBNf]=4.0e-3*bandpass(SSIBNd,samplratSILB,lo,hi,npo,npa,'butter'); 
    [SILBEf]=4.0e-3*bandpass(SILBEd,samplratSILB,lo,hi,npo,npa,'butter'); 
    [SILBNf]=4.0e-3*bandpass(SILBNd,samplratSILB,lo,hi,npo,npa,'butter'); 
end
fltr=cputime-t;
t=cputime;

%Decimate the 100 sps data:
% SSIBEfd = interp1(timsSSIB,SSIBEf,timsPGC,'linear')';
% SSIBNfd = interp1(timsSSIB,SSIBNf,timsPGC,'linear')';
% SILBEfd = interp1(timsSILB,SILBEf,timsPGC,'linear')';
% SILBNfd = interp1(timsSILB,SILBNf,timsPGC,'linear')';
% dmate=cputime-t
% t=cputime;
% SSIBEfd = resample(SSIBEf,2,5); %Now done before filtering
% SSIBNfd = resample(SSIBNf,2,5);
% SILBEfd = resample(SILBEf,2,5);
% SILBNfd = resample(SILBNf,2,5);
dmate=cputime-t
t=cputime;

PGC=PGCEf+1i*PGCNf;
SSIB=SSIBEf+1i*SSIBNf;
SILB=SILBEf+1i*SILBNf;
PGCrot=PGC*exp(1i*rotPGC);
SSIBrot=SSIB*exp(1i*rotSSIB);
SILBrot=SILB*exp(1i*rotSILB);

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

realPGC=real(PGCrot);
realSSIB=real(SSIBrot);
realSILB=real(SILBrot);
% hilPGC=abs(hilbert(realPGC));
% hilSSIB=abs(hilbert(realSSIB));
% hilSILB=abs(hilbert(realSILB));

% [PGChilb,locs]=findpeaks(hilPGC);
% PGChilbtims=timsPGC(locs);
% [SSIBhilb,locs]=findpeaks(hilSSIB);
% SSIBhilbtims=timsPGC(locs);
% [SILBhilb,locs]=findpeaks(hilSILB);
% SILBhilbtims=timsPGC(locs);
% 
% PGCpeaks=[PGChilbtims' PGChilb];
% SSIBpeaks=[SSIBhilbtims' SSIBhilb];
% SILBpeaks=[SILBhilbtims' SILBhilb];
% % realPGC=imag(PGCrot); %A quick/dirty way to try the same on the orthogonal components
% % realSSIB=imag(SSIBrot);
% % realSILB=imag(SILBrot);

%timsPGC = resample(timsPGC,2,5); I think for 100 sps.
PGCtrace=[timsPGCd' realPGC];
SSIBtrace=[timsSSIBd' realSSIB];
SILBtrace=[timsSILBd' realSILB];

figure
plot(PGCtrace(:,1),PGCtrace(:,2),'r')
hold on
plot(SSIBtrace(:,1),SSIBtrace(:,2),'b')
plot(SILBtrace(:,1),SILBtrace(:,2),'k')
xlim([8000 8300])  %xlim([8131 8148])  %xlim([79360 79376])  %xlim([32612 32632])  %

fid = fopen(['HILBERTS/PGCTRIO/PGC.',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',int2str(samplratPGC),'sps'],'w');
fprintf(fid,'%13.5f %9.5f\n',PGCtrace(:,:)');
fclose(fid);
fid = fopen(['HILBERTS/PGCTRIO/SSIB.',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',int2str(samplratSSIB),'sps'],'w');
fprintf(fid,'%13.5f %9.5f\n',SSIBtrace(:,:)');
fclose(fid);
fid = fopen(['HILBERTS/PGCTRIO/SILB.',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',int2str(samplratSILB),'sps'],'w');
fprintf(fid,'%13.5f %9.5f\n',SILBtrace(:,:)');
fclose(fid);


end


