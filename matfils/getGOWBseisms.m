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

%sequence is PG SS SI KLNB SILBz GOWBz 

   timoffrot(1,:)=[2005 253 +086 +020 -004 -180 -080  80 115  50  90];
   timoffrot(2,:)=[2005 254 +086 +020 -004 -180 -080  80 115  50  90];
   timoffrot(3,:)=[2005 255 +086 +020 -004 -180 -080  80 115  50  90];
   timoffrot(4,:)=[2005 256 +086 +020 -004 -180 -080  80 115  50  90];
   timoffrot(5,:)=[2005 257 +086 +020 -004 -180 -080  80 115  50  90];
   timoffrot(6,:)=[2005 258 +086 +020 -004 -180 -080  80 115  50  90];
   timoffrot(7,:)=[2005 259 +086 +020 -004 -180 -080  80 115  50  90];
   timoffrot(8,:)=[2005 260 +086 +020 -004 -180 -080  80 115  50  90];
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
%sequence is PG SS SI KLNB SILBz GOWBz 
%timoffrot(1,:)=[2005 253 +086 +020 -004 -180 -080  80 115  50  90];
MO=day2month(jday,year)
rotPGC=pi*(timoffrot(nd,8)-90)/180;
rotSSIB=pi*(timoffrot(nd,9)-90)/180;
rotSILB=pi*(timoffrot(nd,10)-90)/180;
rotKLNB=pi*(timoffrot(nd,11)-90)/180;
SSIBsoff=timoffrot(nd,3);
SILBsoff=timoffrot(nd,4);
KLNBsoff=timoffrot(nd,5);
SILZsoff=timoffrot(nd,6);
GOWZsoff=timoffrot(nd,7);
SSIBtoff=SSIBsoff/40;
SILBtoff=SILBsoff/40;
KLNBtoff=KLNBsoff/40;
SILZtoff=SILZsoff/40;
GOWZtoff=GOWZsoff/40;
IDENTIF=[YEAR,'.',JDAY,'.',int2str(timoffrot(nd,3)),'.',int2str(timoffrot(nd,4)), ...
    '.',int2str(timoffrot(nd,8)),'.',int2str(timoffrot(nd,9)),'.',int2str(timoffrot(nd,10))]

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
SILBZdat=[prename,'.SILB..HHZ.D.SAC'];
KLNBEdat=[prename,'.KLNB..HHE.D.SAC'];
KLNBNdat=[prename,'.KLNB..HHN.D.SAC'];
% GOWBZdat=[prename,'.GOWB..HHZ.D.SAC'];

[PGCE,HdrDataPGC,tnuPGC,pobjPGC,timsPGC]=readsac(PGCEdat,0,'l');
[PGCN,~,~,~,~]=readsac(PGCNdat,0,'l');
%[PGCZ,~,~,~,~]=readsac(PGCZdat,0,'l');
[SSIBE,HdrDataSSIB,tnuSSIB,pobjSSIB,timsSSIB]=readsac(SSIBEdat,0,'l');
[SSIBN,~,~,~,~]=readsac(SSIBNdat,0,'l');
[SILBE,HdrDataSILB,tnuSILB,pobjSILB,timsSILB]=readsac(SILBEdat,0,'l');
[SILBN,~,~,~,~]=readsac(SILBNdat,0,'l');
[SILBZ,~,~,~,~]=readsac(SILBZdat,0,'l');
[KLNBE,HdrDataKLNB,tnuKLNB,pobjKLNB,timsKLNB]=readsac(KLNBEdat,0,'l');
[KLNBN,~,~,~,~]=readsac(KLNBNdat,0,'l');
% [GOWBE,HdrDataGOWB,tnuGOWB,pobjGOWB,timsGOWB]=readsac(GOWBZdat,0,'l');
samplratPGC=round(1/(timsPGC(2)-timsPGC(1)))
samplratSSIB=round(1/(timsSSIB(2)-timsSSIB(1)))
samplratSILB=round(1/(timsSILB(2)-timsSILB(1)))
samplratKLNB=round(1/(timsKLNB(2)-timsKLNB(1)))
% samplratGOWB=round(1/(timsGOWB(2)-timsGOWB(1)))
%
rdsac=cputime-t;
t=cputime;

tracelenPGC=length(PGCE);
tracelenSSIB=length(SSIBE);
tracelenSILB=length(SILBE);
tracelenKLNB=length(KLNBE);
% tracelenGOWB=length(GOWBE);

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
SILBZ(1:100)=sin(x).*SILBZ(1:100);
KLNBE(1:100)=sin(x).*KLNBE(1:100);
KLNBN(1:100)=sin(x).*KLNBN(1:100);
% GOWBZ(1:100)=sin(x).*GOWBZ(1:100);
x=flipud(x);
SSIBE(tracelenSSIB-99:tracelenSSIB)=sin(x).*SSIBE(tracelenSSIB-99:tracelenSSIB);
SSIBN(tracelenSSIB-99:tracelenSSIB)=sin(x).*SSIBN(tracelenSSIB-99:tracelenSSIB);
SILBE(tracelenSILB-99:tracelenSILB)=sin(x).*SILBE(tracelenSILB-99:tracelenSILB);
SILBN(tracelenSILB-99:tracelenSILB)=sin(x).*SILBN(tracelenSILB-99:tracelenSILB);
SILBZ(tracelenSILB-99:tracelenSILB)=sin(x).*SILBZ(tracelenSILB-99:tracelenSILB);
KLNBE(tracelenKLNB-99:tracelenKLNB)=sin(x).*KLNBE(tracelenKLNB-99:tracelenKLNB);
KLNBN(tracelenKLNB-99:tracelenKLNB)=sin(x).*KLNBN(tracelenKLNB-99:tracelenKLNB);
% GOWBZ(tracelenGOWB-99:tracelenGOWB)=sin(x).*GOWBZ(tracelenGOWB-99:tracelenGOWB);

%Filter data:
hi=0.05;
lo=0.02;
resampPGC=[1 10];
resampSSIB=[1 25];
resampSILB=[1 25];
resampKLNB=[1 25];
% resampGOWB=[1 25];
samplratPGC=samplratPGC*resampPGC(1)/resampPGC(2)
samplratSSIB=samplratSSIB*resampSSIB(1)/resampSSIB(2)
samplratSILB=samplratSILB*resampSILB(1)/resampSILB(2);
samplratKLNB=samplratKLNB*resampKLNB(1)/resampKLNB(2);
% samplratGOWB=samplratGOWB*resampGOWB(1)/resampGOWB(2);
PGCEd = resample(PGCE,resampPGC(1),resampPGC(2));
PGCNd = resample(PGCN,resampPGC(1),resampPGC(2));
timsPGCd = (0:1/samplratPGC:86400-1/samplratPGC);
SSIBEd = resample(SSIBE,resampSSIB(1),resampSSIB(2));
SSIBNd = resample(SSIBN,resampSSIB(1),resampSSIB(2));
timsSSIBd = (0:1/samplratSSIB:86400-1/samplratSSIB);
SILBEd = resample(SILBE,resampSILB(1),resampSILB(2));
SILBNd = resample(SILBN,resampSILB(1),resampSILB(2));
SILBZd = resample(SILBZ,resampSILB(1),resampSILB(2));
timsSILBd = (0:1/samplratSILB:86400-1/samplratSILB);
KLNBEd = resample(KLNBE,resampKLNB(1),resampKLNB(2));
KLNBNd = resample(KLNBN,resampKLNB(1),resampKLNB(2));
timsKLNBd = (0:1/samplratKLNB:86400-1/samplratKLNB);
% GOWBZd = resample(GOWBZ,resampGOWB(1),resampGOWB(2));
% timsGOWBd = (0:1/samplratGOWB:86400-1/samplratGOWB);
tracelenPGC=length(PGCEd)
tracelenSSIB=length(SSIBEd)
tracelenSILB=length(SILBEd);
tracelenKLNB=length(KLNBEd);
% tracelenGOWB=length(GOWBZd);
SSIBsoff=round(SSIBsoff*resampPGC(1)/resampPGC(2))
SILBsoff=round(SILBsoff*resampPGC(1)/resampPGC(2));
KLNBsoff=round(KLNBsoff*resampPGC(1)/resampPGC(2));
% GOWsoff=round(GOWBsoff*resampPGC(1)/resampPGC(2));

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
    [SILBZf]=20.0e-3*bandpass(SILBZd,samplratSSIB,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
    [KLNBEf]=20.0e-3*bandpass(KLNBEd,samplratSSIB,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0* 
    [KLNBNf]=20.0e-3*bandpass(KLNBNd,samplratSSIB,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
%     [GOWBZf]=20.0e-3*bandpass(GOWBZd,samplratGOWB,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
else
    [SSIBEf]=4.0e-3*bandpass(SSIBEd,samplratSILB,lo,hi,npo,npa,'butter'); 
    [SSIBNf]=4.0e-3*bandpass(SSIBNd,samplratSILB,lo,hi,npo,npa,'butter'); 
    [SILBEf]=4.0e-3*bandpass(SILBEd,samplratSILB,lo,hi,npo,npa,'butter'); 
    [SILBNf]=4.0e-3*bandpass(SILBNd,samplratSILB,lo,hi,npo,npa,'butter'); 
    [SILBZf]=4.0e-3*bandpass(SILBZd,samplratSILB,lo,hi,npo,npa,'butter'); 
    [KLNBEf]=4.0e-3*bandpass(KLNBEd,samplratSILB,lo,hi,npo,npa,'butter'); 
    [KLNBNf]=4.0e-3*bandpass(KLNBNd,samplratSILB,lo,hi,npo,npa,'butter'); 
%     [GOWBZf]=4.0e-3*bandpass(GOWBZd,samplratSILB,lo,hi,npo,npa,'butter'); 
end
fltr=cputime-t;
t=cputime;

dmate=cputime-t
t=cputime;

PGC=PGCEf+1i*PGCNf;
SSIB=SSIBEf+1i*SSIBNf;
SILB=SILBEf+1i*SILBNf;
SILBZ=SILBZf;
% GOWBZ=GOWBZf;
KLNB=KLNBEf+1i*KLNBNf;
PGCrot=PGC*exp(1i*rotPGC);
SSIBrot=SSIB*exp(1i*rotSSIB);
SILBrot=SILB*exp(1i*rotSILB);
KLNBrot=KLNB*exp(1i*rotKLNB);

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
if KLNBtoff > 0
    KLNBrot(1:tracelenPGC-KLNBsoff)=KLNBrot(KLNBsoff+1:tracelenPGC);
    KLNBrot(tracelenPGC-KLNBsoff+1:tracelenPGC)=0;
else
    KLNBrot(-KLNBsoff+1:tracelenPGC)=KLNBrot(1:tracelenPGC+KLNBsoff);
    KLNBrot(1:-KLNBsoff)=0;
end
if SILZtoff > 0
    SILBZ(1:tracelenPGC-SILZsoff)=SILBZ(SILZsoff+1:tracelenPGC);
    SILBZ(tracelenPGC-SILZsoff+1:tracelenPGC)=0;
else
    SILBZ(-SILZsoff+1:tracelenPGC)=SILBZ(1:tracelenPGC+SILZsoff);
    SILBZ(1:-SILZsoff)=0;
end
% if GOWZtoff > 0
%     GOWBZ(1:tracelenPGC-GOWZsoff)=GOWBZ(GOWZsoff+1:tracelenPGC);
%     GOWBZ(tracelenPGC-GOWZsoff+1:tracelenPGC)=0;
% else
%     GOWBZ(-GOWZsoff+1:tracelenPGC)=GOWBZ(1:tracelenPGC+GOWZsoff);
%     GOWBZ(1:-GOWZsoff)=0;
% end

realPGC=real(PGCrot);
realSSIB=real(SSIBrot);
realSILB=real(SILBrot);
realKLNB=real(KLNBrot);

%timsPGC = resample(timsPGC,2,5); I think for 100 sps.
PGCtrace=[timsPGCd' realPGC];
SSIBtrace=[timsSSIBd' realSSIB];
SILBtrace=[timsSILBd' realSILB];
KLNBtrace=[timsKLNBd' realKLNB];
SILBZtrace=[timsSILBd' SILBZ];
% GOWBZtrace=[timsGOWBd' GOWBZ];

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
fid = fopen(['HILBERTS/PGCTRIO/KLNB.',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',int2str(samplratKLNB),'sps'],'w');
fprintf(fid,'%13.5f %9.5f\n',KLNBtrace(:,:)');
fclose(fid);
fid = fopen(['HILBERTS/PGCTRIO/SILBZ.',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',int2str(samplratSILB),'sps'],'w');
fprintf(fid,'%13.5f %9.5f\n',SILBZtrace(:,:)');
fclose(fid);
% fid = fopen(['HILBERTS/PGCTRIO/GOWBZ.',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',int2str(samplratGOWB),'sps'],'w');
% fprintf(fid,'%13.5f %9.5f\n',GOWBZtrace(:,:)');
% fclose(fid);

end


