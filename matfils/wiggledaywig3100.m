%Sits at one spot (rotations; delays) for one day and looks for the cross-correlation.
%Now also plots aligned 4-s windows.  As of July 2012.
format short e
clear all
close all
set(0,'DefaultFigureVisible','on');
%    timoffrot(1,:)=[2005 258 03 50 44 +054 +064  70 0 35];
%    timoffrot(2,:)=[2005 259 03 50 44 +054 +064  70 0 35];
%    timoffrot(3,:)=[2005 260 03 50 44 +054 +064  70 0 35];
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
%    timoffrot(4,:)=[2005 255 07 42 12 +074 +024  90   5  40];
%    timoffrot(5,:)=[2005 256 07 42 12 +074 +024  90   5  40];
%    timoffrot(6,:)=[2005 257 07 42 12 +074 +024  90   5  40];
%    timoffrot(7,:)=[2005 258 07 42 12 +074 +024  90   5  40];
% 
%    timoffrot(8,:)=[2005 255 07 42 12 +069 +039  80   5  45];
%    timoffrot(9,:)=[2005 256 07 42 12 +069 +039  80   5  45];
%    timoffrot(10,:)=[2005 257 07 42 12 +069 +039  80   5  45];
%    timoffrot(11,:)=[2005 258 07 42 12 +069 +039  80   5  45];
  

%    timoffrot(1,:)=[2005 253 07 42 12 +086 +020  80 115  50];
%   timoffrot(1,:)=[2005 254 07 42 12 +086 +020  80 115  50];
    timoffrot(1,:)=[2005 255 07 42 12 +086 +020  80 115  50];
%    timoffrot(4,:)=[2005 256 07 42 12 +086 +020  80 115  50];
%    timoffrot(5,:)=[2005 257 07 42 12 +086 +020  80 115  50];
%    timoffrot(6,:)=[2005 258 07 42 12 +086 +020  80 115  50];
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
%   timoffrot(2,:)=[2004 196 07 42 12 +086 +020  80 115  50];
%   timoffrot(3,:)=[2004 197 07 42 12 +086 +020  80 115  50];
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
%   timoffrot(5,:)=[2003 062 07 42 12 +086 +020  80 115  50];
%   timoffrot(6,:)=[2003 063 07 42 12 +086 +020  80 115  50];
%   timoffrot(7,:)=[2003 064 07 42 12 +086 +020  80 115  50];
%    timoffrot(31,:)=[2003 065 07 42 12 +086 +020  80 115  50];
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
%    timoffrot(5,:)=[2004 196 07 42 12 +086 +020  80 115  50];
%    timoffrot(6,:)=[2004 197 07 42 12 +086 +020  80 115  50];
%    timoffrot(7,:)=[2004 198 07 42 12 +086 +020  80 115  50];
%    timoffrot(8,:)=[2004 199 07 42 12 +086 +020  80 115  50];
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
%Rotations & offsets:
%    timoffrot=[2005 255 07 42 12 +085 +042  75  85  45];
%    timoffrot=[2005 256 07 42 12 +085 +042  75  85  45];
%    timoffrot=[2005 257 07 42 12 +085 +042  75  85  45];
%    timoffrot=[2005 258 07 42 12 +085 +042  75  85  45];/+045
% 
%    timoffrot=[2005 255 07 42 12 +080 +053  75  50  50];
%    timoffrot=[2005 256 07 42 12 +080 +053  75  50  50];
%    timoffrot=[2005 257 07 42 12 +080 +053  75  50  50];
%    timoffrot=[2005 258 07 42 12 +080 +053  75  50  50];

   %timoffrot=[2005 246 04 25 16 +060 -085  90  60  45]; %agu
   %timoffrot=[2005 254 04 05 00 +087 +018  85 110  45];
   %timoffrot=[2005 255 04 05 00 +087 +018  85 110  45];
   %timoffrot=[2005 254 07 42 12 +086 +016  80 110  60];
   %timoffrot=[2005 262 05 19 48 +085 +020  80 115  55];   
   %timoffrot=[2005 254 05 19 48 +085 +020  80 115  55]; %Started yesterday.  
   %timoffrot=[2005 254 05 19 48 +086 +019  80 105  55]; %Started yesterday.  
   %timoffrot=[2005 256 05 19 48 +085 +020  80 115  55]; %Started yesterday.  
   %timoffrot=[2005 255 01 27 56 +086 +027  75 125  45]; %Started yesterday 
   %timoffrot=[2005 257 05 19 48 +085 +020  80 115  55]; %v. little
   %timoffrot=[2005 260 05 19 48 +085 +020  80 115  55]; % Sept. 18.  Nothing at this spot.
   %timoffrot=[2005 259 12 58 36 +036 +075  80 365 385]; %New spot
   %timoffrot=[2005 260 12 58 36 +036 +075  80 365 385]; %New spot
   %timoffrot=[2005 259 03 50 44 +045 +072  65 355 370]; %New spot
   %timoffrot=[2005 260 03 50 44 +045 +072  65 355 370]; %New spot
   %timoffrot=[2004 193 06 52 36 +105 -031  80  45  80]; %Earlier spot 2004
   %timoffrot=[2004 199 14 15 32 +085 +020  80 115  55];
   %timoffrot=[2004 196 14 15 32 +086 +013  80 115  55];
   %timoffrot=[2005 254 14 15 32 +086 +013  80 115  55];
   %timoffrot=[2004 196 14 15 32 +085 +020  80 110  55];
   %timoffrot=[2004 197 14 15 32 +086 +027  75 120  45]; %downstream
   %timoffrot=[2004 198 14 15 32 +086 +027  75 120  45]; %downstream
   %timoffrot=[2004 198 14 15 32 +085 +020  80 115  55];
   %timoffrot=[2004 196  9 45 48   87   14  75 115  50]; %focus on early, 2004, 2-8Hz catalog
   %timoffrot=[2003 063 14 15 32 +085 +020  80 115  55]; %spot, 2003
   %timoffrot=[2003 062 14 15 32 +085 +020  80 115  55]; %spot, 2003
   %timoffrot=[2003 063 14 15 32 +086 +027  75 125  45]; %downstream, 2003
   %timoffrot=[2003 062 14 15 32 +086 +013  80 110  55]; %upstream, 2003
   %timoffrot=[2003 061 14 15 32 +085 +020  80 115  55]; %faked %COULD BE 80 120 50/55, 2-8Hz
   %timoffrot=[2003 069 14 15 32 +085 +020  80 115  55]; %spot, 2003.  NADA Anywhere! on 065
   %timoffrot=[2005 250 05 19 48 +092 -037  80 70  75];   
   %timoffrot=[2003 062  9 45 48   87   14  75 115  50]; %focus on early, 2003, 2-8Hz catalog
   %timoffrot=[2005 257 05 19 48 +085 +043  75 80  45];   
mshift=19; %16 %11 %14 %12_for_+085+020 %15_for_+059-073
hi=6.;
lo=1.5;
% hi=8.;
% lo=2.;
% hi=6.;
% lo=0.75;
loopoffmax=1.5
xcmaxAVEnmin=0.4 %0.4 for 4s 1.5-6 Hz; 0.36 for 4s2-8 Hz ;  0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?
winlen=4*40;
%winlen=128*40;
winoff=40/1;
%winoff=40*8;
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

%Read Armbruster's detections:
ArmCat=load(['/data2/arubin/CNDC/',YEAR,'/list.',YEAR,'.pgsssiMAJ']);
%ArmCat=load(['/data2/arubin/CNDC/',YEAR,'/list.',YEAR,'.pgsssiMIN']);
detects=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
vectoff=[ArmCat(:,5)-SSIBsoff ArmCat(:,6)-SILBsoff (ArmCat(:,6)-SILBsoff)-(ArmCat(:,5)-SSIBsoff)];
distoff=vectoff.*vectoff;
distoff2=sum(distoff,2);
gridoff=max(abs(vectoff),[],2); %Try this instead.  Should change w/mshift (defined later, though)
m0=0;
Detects=0;
for n=1:length(detects)
    %if distoff2(n)<=9
    if gridoff(n)<=mshift-1
         m0=m0+1;
         Detects(m0)=detects(n);
     end
end
%abs(xmaxPGSIn(n)-xmaxPGSSn(n)+xmaxSISSn(n))>loopoffma
ArmCat=load(['/data2/arubin/CNDC/',YEAR,'/list.',YEAR,'.pgsssiMAJ_new']);
detects2_8=86400*(ArmCat(:,1)-jday)+3600*ArmCat(:,2)+60*ArmCat(:,3)+ArmCat(:,4);
vectoff=[ArmCat(:,5)-SSIBsoff ArmCat(:,6)-SILBsoff];
distoff=vectoff.*vectoff;
distoff2=sum(distoff,2);
gridoff=max(abs(vectoff),[],2); %Try this instead.  Should change w/mshift (defined later, though)
m1=0;
Detects2_8=0;
for n=1:length(detects2_8)
    %if distoff2(n)<=9
    if gridoff(n)<=mshift-1
         m1=m1+1;
         Detects2_8(m1)=detects2_8(n);
     end
end

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
% timsPGC(1)
% timsSSIB(1)
% timsSILB(1)
%
rdsac=cputime-t
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

winbig=2*(tracelenPGC/2-120);
timbig=winbig/(2*40);
igstart=floor(tracelenPGC/2-winbig/2)+1;
nwin=floor((winbig-winlen)/winoff);
spsfactor=1;
nzeros=glitches2(PGCE,nwin,winlen,winoff,igstart,spsfactor);
nzerosPGCE=nzeros;
nzeros=glitches2(PGCN,nwin,winlen,winoff,igstart,spsfactor);
nzerosPGCN=nzeros;
%nzerosPGC=nzerosPGCE+nzerosPGCN;
nzerosPGC=max(nzerosPGCE,nzerosPGCN);
spsfactor=2.5;
nzeros=glitches2(SSIBE,nwin,winlen,winoff,igstart,spsfactor);
nzerosSSIBE=nzeros;
nzeros=glitches2(SSIBN,nwin,winlen,winoff,igstart,spsfactor);
nzerosSSIBN=nzeros;
%nzerosSSIB=nzerosSSIBE+nzerosSSIBN;
nzerosSSIB=max(nzerosSSIBE,nzerosSSIBN);
nzeros=glitches2(SILBE,nwin,winlen,winoff,igstart,spsfactor);
nzerosSILBE=nzeros;
nzeros=glitches2(SILBN,nwin,winlen,winoff,igstart,spsfactor);
nzerosSILBN=nzeros;
%nzerosSILB=nzerosSILBE+nzerosSILBN;
nzerosSILB=max(nzerosSILBE,nzerosSILBN);
%pause

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
npo=2;
npa=1;
[PGCEf]=1.6e-4*bandpass(PGCE,40,lo,hi,npo,npa,'butter');
[PGCNf]=1.6e-4*bandpass(PGCN,40,lo,hi,npo,npa,'butter');
%[PGCZf]=1.6e-4*bandpass(PGCZ,40,lo,hi,npo,npa,'butter');
if year==2003 && jday<213
    [SSIBEf]=20.0e-3*bandpass(SSIBE,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
    [SSIBNf]=20.0e-3*bandpass(SSIBN,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
    [SILBEf]=20.0e-3*bandpass(SILBE,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0* 
    [SILBNf]=20.0e-3*bandpass(SILBN,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
else
    [SSIBEf]=4.0e-3*bandpass(SSIBE,100,lo,hi,npo,npa,'butter'); 
    [SSIBNf]=4.0e-3*bandpass(SSIBN,100,lo,hi,npo,npa,'butter'); 
    [SILBEf]=4.0e-3*bandpass(SILBE,100,lo,hi,npo,npa,'butter'); 
    [SILBNf]=4.0e-3*bandpass(SILBN,100,lo,hi,npo,npa,'butter'); 
end
fltr=cputime-t
t=cputime;

%Decimate the 100 sps data:
% SSIBEfd = interp1(timsSSIB,SSIBEf,timsPGC,'linear')';
% SSIBNfd = interp1(timsSSIB,SSIBNf,timsPGC,'linear')';
% SILBEfd = interp1(timsSILB,SILBEf,timsPGC,'linear')';
% SILBNfd = interp1(timsSILB,SILBNf,timsPGC,'linear')';
% dmate=cputime-t
% t=cputime;
SSIBEfd = resample(SSIBEf,2,5);
SSIBNfd = resample(SSIBNf,2,5);
SILBEfd = resample(SILBEf,2,5);
SILBNfd = resample(SILBNf,2,5);
dmate=cputime-t
t=cputime;

PGC=PGCEf+1i*PGCNf;
SSIB=SSIBEfd+1i*SSIBNfd;
SILB=SILBEfd+1i*SILBNfd;
PGCrot=PGC*exp(1i*rotPGC);
SSIBrot=SSIB*exp(1i*rotSSIB);
SILBrot=SILB*exp(1i*rotSILB);
    %rtate=cputime-t
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

realPGC=real(PGCrot);
realSSIB=real(SSIBrot);
realSILB=real(SILBrot);
% realPGC=imag(PGCrot); %A quick/dirty way to try the same on the orthogonal components
% realSSIB=imag(SSIBrot);
% realSILB=imag(SILBrot);

plotstart=1; %(timlook-timstart)*40
plotend=tracelenPGC; %(timlook-timstart)*40
ltmax=max([real(PGCrot(plotstart:plotend)); ...
          real(SSIBrot(plotstart:plotend)); ...
          real(SILBrot(plotstart:plotend))]);
ltmax=min(ltmax,0.75);
ltmin=min([real(PGCrot(plotstart:plotend)); ...
          real(SSIBrot(plotstart:plotend)); ...
          real(SILBrot(plotstart:plotend))]);
ltmin=max(ltmin,-0.75);

figure %Superposed traces, 150s window:
subplot(4,1,1); 
hold on
plot(timsPGC,realPGC,'r');
plot(timsPGC,realSSIB,'b');
plot(timsPGC,realSILB,'k');
axis([0 timwin/2 ltmin ltmax]);
%xlim([timlook-75 timlook-37.5]);
subplot(4,1,2); 
hold on
plot(timsPGC,realPGC,'r');
plot(timsPGC,realSSIB,'b');
plot(timsPGC,realSILB,'k');
axis([timwin/2 timwin ltmin ltmax]);
%xlim([timlook-37.5 timlook]);
subplot(4,1,3); 
hold on
plot(timsPGC,realPGC,'r');
plot(timsPGC,realSSIB,'b');
plot(timsPGC,realSILB,'k');
axis([timwin 3*timwin/2 ltmin ltmax]);
%xlim([timlook timlook+37.5]);
subplot(4,1,4); 
hold on
plot(timsPGC,realPGC,'r');
plot(timsPGC,realSSIB,'b');
plot(timsPGC,realSILB,'k');
axis([3*timwin/2 2*timwin ltmin ltmax]);
%xlim([timlook+37.5 timlook+75]);
traceplt=cputime-t
t=cputime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Autocorrelation of stations.  Those that end in "2" are the running
%   cumulative sum, to be used later by differncing the window edpoints.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PGCauto=realPGC.*realPGC;
PGC2=cumsum(PGCauto);
SSIBauto=realSSIB.*realSSIB;
SSIB2=cumsum(SSIBauto);
SILBauto=realSILB.*realSILB;
SILB2=cumsum(SILBauto);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Cross-correlation between stations, with small offsets up to +/- mshift.
%  First index is pointwise multiplication of traces; second is shifting offset.
%  lenx is shorter than tracelenPGC by mshift at each end (see notebook sketch)
%  For PGSS and PGSI, SSI and SIL are shifted relative to PGC, by 1 each time through loop.
%  For SISS, SSI is shifted relative to SILB.
%  cumsumPGSS etc. are the running cumulative sum of the x-correlation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mshift=16;
lenx=tracelenPGC-2*mshift;
PGSSx=zeros(lenx, 2*mshift+1);
PGSIx=zeros(lenx, 2*mshift+1);
SISSx=zeros(lenx, 2*mshift+1);
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
xcshifts=cputime-t
t=cputime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  "winbig" is now the whole day, minus 3 sec at each end (apparently).
%  "timbig" is the time of half that.
%  igstart is the index of the starting sample.
%  winlen and winoff refer to the small windows.
%  timswin refers to the central times of those small windows.
%  sumsPGSS (etc.) is the cross-correlation sum over the window.  The first
%    index refers to the window number and the second the shift over +/-mshift.
%  Normalized x-correlation:
%    For PGSS and PGSI, for a given window PGC does not shift but SSI and 
%    SIL do.  So can compute sumsPGC2 (from the running cum. sum PGC2) just
%    once for each window.  Same for sumsSILB2b.  But for the stations that
%    shift, SSI and SIL (for PGC) and SSI (for SIL), must compute sumsSSIB2 
%    and sumsSILB2 upon each shift (actually, this is is easy book-keeping
%    but not efficient).  Again, the first index refers to the window
%    number and the second the shift over +/-mshift.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% winbig=2*(tracelenPGC/2-120);
% timbig=winbig/(2*40);
% igstart=floor(tracelenPGC/2-winbig/2)+1;
% nwin=floor((winbig-winlen)/winoff);
timswin=zeros(nwin,1);
sumsPGSS=zeros(nwin,2*mshift+1);
sumsPGSI=zeros(nwin,2*mshift+1);
sumsSISS=zeros(nwin,2*mshift+1);
sumsPGC2=zeros(nwin,2*mshift+1);
sumsSSIB2=zeros(nwin,2*mshift+1);
sumsSILB2=zeros(nwin,2*mshift+1);
sumsSILB2b=zeros(nwin,2*mshift+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  sumsPGSS is shorter than sumsPGC2 by 2*mshift.  This is why sumsPGC2 etc
%  is shifted by +mshift.  cumsumPGSS(1,:)=cumsum(PGSSx)(1,:) starts mshift
%  to the right of the first data sample.  igstart is how many to the right
%  of that.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:nwin;
    istart=igstart+(n-1)*winoff;
    iend=istart+winlen;
    timswin(n)=timsPGC(istart+winlen/2); 
    sumsPGSS(n,:)=cumsumPGSS(iend,:)-cumsumPGSS(istart-1,:); 
    sumsPGSI(n,:)=cumsumPGSI(iend,:)-cumsumPGSI(istart-1,:);
    sumsSISS(n,:)=cumsumSISS(iend,:)-cumsumSISS(istart-1,:);
    sumsPGC2(n,:)=PGC2(iend+mshift)-PGC2(istart+mshift-1);  %PGC2 is cumsummed. Yes, +mshift.
    sumsSILB2b(n,:)=SILB2(iend+mshift)-SILB2(istart+mshift-1); %Similar, for the SILB-SSIB connection.
    for m=-mshift:mshift;
        sumsSSIB2(n,m+mshift+1)=SSIB2(iend+mshift-m)-SSIB2(istart+mshift-1-m); %+m??? (yes).
        sumsSILB2(n,m+mshift+1)=SILB2(iend+mshift-m)-SILB2(istart+mshift-1-m);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Denominator for the normalization.  A 2D array, nwin by 2*mshift+1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%An attempt to bypass glitches in data.  Min value of good data typically ~10^{-2}
glitches=1.e-7;
sumsPGC2=max(sumsPGC2,glitches);
sumsSSIB2=max(sumsSSIB2,glitches);
sumsSILB2=max(sumsSILB2,glitches);
%
denomPGSSn=realsqrt(sumsPGC2.*sumsSSIB2);
denomPGSIn=realsqrt(sumsPGC2.*sumsSILB2);
denomSISSn=realsqrt(sumsSILB2b.*sumsSSIB2);
%
% denomPGSSn=max(denomPGSSn,eps);
% denomPGSIn=max(denomPGSIn,eps);
% denomSISSn=max(denomSISSn,eps);
%
sumsPGSSn=sumsPGSS./denomPGSSn;
sumsPGSIn=sumsPGSI./denomPGSIn;
sumsSISSn=sumsSISS./denomSISSn;
[xcmaxPGSSn,imaxPGSS]=max(sumsPGSSn,[],2); %Integer-offset max cross-correlation
[xcmaxPGSIn,imaxPGSI]=max(sumsPGSIn,[],2);
[xcmaxSISSn,imaxSISS]=max(sumsSISSn,[],2);
%Parabolic fit:
% [xmaxPGSSn,ymaxPGSSn,xloPGSSn,xhiPGSSn]=parabol(nwin,mshift,sumsPGSSn,imaxPGSS); %These return some poor measure of "confidence".
% [xmaxPGSIn,ymaxPGSIn,xloPGSIn,xhiPGSIn]=parabol(nwin,mshift,sumsPGSIn,imaxPGSI);
% [xmaxSISSn,ymaxSISSn,xloSISSn,xhiSISSn]=parabol(nwin,mshift,sumsSISSn,imaxSISS);
[xmaxPGSSn,ymaxPGSSn,aPGSS]=parabol(nwin,mshift,sumsPGSSn,imaxPGSS); %Interpolated max cross-correlation
[xmaxPGSIn,ymaxPGSIn,aPGSI]=parabol(nwin,mshift,sumsPGSIn,imaxPGSI);
[xmaxSISSn,ymaxSISSn,aSISS]=parabol(nwin,mshift,sumsSISSn,imaxSISS);

ix=sub2ind(size(denomPGSSn),(1:nwin)',imaxPGSS);
ampPGSS=sqrt(denomPGSSn(ix)); %This makes amplitude linear rather than quadratic with counts.
ampPGC2=sumsPGC2(ix); %by construction PGC2 is the same for all shifts
ampSSIB2=sumsSSIB2(ix);
ix=sub2ind(size(denomPGSIn),(1:nwin)',imaxPGSI);
ampPGSI=sqrt(denomPGSIn(ix));
ampSILB2=sumsSILB2(ix);
ix=sub2ind(size(denomSISSn),(1:nwin)',imaxSISS);
ampSISS=sqrt(denomSISSn(ix));
AmpComp(1:4)=0;
AmpComp(5:nwin)=((ampPGC2(5:nwin)+ampSSIB2(5:nwin)+ampSILB2(5:nwin))- ...
                (ampPGC2(1:nwin-4)+ampSSIB2(1:nwin-4)+ampSILB2(1:nwin-4)))./ ...
                ((ampPGC2(5:nwin)+ampSSIB2(5:nwin)+ampSILB2(5:nwin))+ ...
                (ampPGC2(1:nwin-4)+ampSSIB2(1:nwin-4)+ampSILB2(1:nwin-4))) ;
%Center them
imaxPGSS=imaxPGSS-mshift-1;
imaxPGSI=imaxPGSI-mshift-1;
imaxSISS=imaxSISS-mshift-1;
iloopoff=imaxPGSI-imaxPGSS+imaxSISS; %How well does the integer loop close?
xmaxPGSSn=xmaxPGSSn-mshift-1;
xmaxPGSIn=xmaxPGSIn-mshift-1;
xmaxSISSn=xmaxSISSn-mshift-1;
loopoff=xmaxPGSIn-xmaxPGSSn+xmaxSISSn; %How well does the interpolated loop close?
xcmaxAVEn=(xcmaxPGSSn+xcmaxPGSIn+xcmaxSISSn)/3;
xcnshifts=cputime-t
t=cputime;

ampmax=max([ampPGSS; ampPGSI; ampSISS]);
% figure 
% subplot(4,1,1,'align'); 
% hold on
% plot(timswin,zeros(nwin,1),'k:');
% %plot(timswin,3.5+zeros(nwin,1),'k:');
% plot(timswin,xcmaxAVEn*mshift,'g');
% plot(timswin,xmaxPGSSn,'bs','MarkerSize',2);
% plot(timswin,xmaxPGSIn,'ro','MarkerSize',2);
% plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',2);
% %plot(timswin,xmaxSISSn,'k*','MarkerSize',2);
% axis([0 timbig/2 -mshift mshift]);
% ylabel('samples')
% box on
% subplot(4,1,2,'align'); 
% hold on
% plot(timswin,zeros(nwin,1),'k:');
% %plot(timswin,3.5+zeros(nwin,1),'k:');
% plot(timswin,xcmaxAVEn*mshift,'g');
% plot(timswin,xmaxPGSSn,'bs','MarkerSize',2);
% plot(timswin,xmaxPGSIn,'ro','MarkerSize',2);
% plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',2);
% %plot(timswin,xmaxSISSn,'k*','MarkerSize',2);
% axis([timbig/2 timbig -mshift mshift]);
% ylabel('samples')
% box on
% subplot(4,1,3,'align'); 
% hold on
% plot(timswin,zeros(nwin,1),'k:');
% %plot(timswin,3.5+zeros(nwin,1),'k:');
% plot(timswin,xcmaxAVEn*mshift,'g');
% plot(timswin,xmaxPGSSn,'bs','MarkerSize',2);
% plot(timswin,xmaxPGSIn,'ro','MarkerSize',2);
% plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',2);
% %plot(timswin,xmaxSISSn,'k*','MarkerSize',2);
% axis([timbig 3*timbig/2 -mshift mshift]);
% ylabel('samples')
% box on
% subplot(4,1,4,'align'); 
% hold on
% plot(timswin,zeros(nwin,1),'k:');
% %plot(timswin,3.5+zeros(nwin,1),'k:');
% plot(timswin,xcmaxAVEn*mshift,'g');
% plot(timswin,xmaxPGSSn,'bs','MarkerSize',2);
% plot(timswin,xmaxPGSIn,'ro','MarkerSize',2);
% plot(timswin,(xmaxPGSIn-xmaxPGSSn+xmaxSISSn),'k*','MarkerSize',2);
% %plot(timswin,xmaxSISSn,'k*','MarkerSize',2);
% axis([3*timbig/2 2*timbig -mshift mshift]);
% xlabel('sec')
% ylabel('samples')
% box on
% title([IDENTIF,'_{',int2str(lo),'-',int2str(hi),'}'])
% orient landscape
% print('-depsc',['FIGS/',IDENTIF,'-',int2str(winlen/40),'s_',int2str(lo),'-',int2str(hi),'a.eps'])

medxcmaxAVEn=median(xcmaxAVEn)
xmaxPGSSntmp=xmaxPGSSn;
xmaxPGSIntmp=xmaxPGSIn;
xmaxSISSntmp=xmaxSISSn;
%loopoff=xmaxPGSIn-xmaxPGSSn+xmaxSISSn; 

%Go to 100 sps
srate=100;
PGCEfu = resample(PGCEf,5,2);
PGCNfu = resample(PGCNf,5,2);
PGC=PGCEfu+1i*PGCNfu;
SSIB=SSIBEf+1i*SSIBNf;
SILB=SILBEf+1i*SILBNf;
PGCrot=PGC*exp(1i*rotPGC);
SSIBrot=SSIB*exp(1i*rotSSIB);
SILBrot=SILB*exp(1i*rotSILB);
SSIBsoff=round(2.5*SSIBsoff)
SILBsoff=round(2.5*SILBsoff)
if SSIBtoff > 0
    SSIBrot(1:tracelenSSIB-SSIBsoff)=SSIBrot(SSIBsoff+1:tracelenSSIB);
    SSIBrot(tracelenSSIB-SSIBsoff+1:tracelenSSIB)=0;
else
    SSIBrot(-SSIBsoff+1:tracelenSSIB)=SSIBrot(1:tracelenSSIB+SSIBsoff);
    SSIBrot(1:-SSIBsoff)=0;
end
if SILBtoff > 0
    SILBrot(1:tracelenSSIB-SILBsoff)=SILBrot(SILBsoff+1:tracelenSSIB);
    SILBrot(tracelenSSIB-SILBsoff+1:tracelenSSIB)=0;
else
    SILBrot(-SILBsoff+1:tracelenSSIB)=SILBrot(1:tracelenSSIB+SILBsoff);
    SILBrot(1:-SILBsoff)=0;
end
realPGC100=real(PGCrot);
realSSIB100=real(SSIBrot);
realSILB100=real(SILBrot);
igstart100=floor(tracelenSSIB/2-2.5*winbig/2)+1; %This is close to time of igstart at 40 sps
mshift100=ceil(2.5*mshift);
winlen100=round(2.5*winlen);
winoff100=round(2.5*winoff);
loopoffmax100=round(loopoffmax*2.5)

nin=0;
zerosallowed=20*winlen/160;  %keep as winlen, not winlen100; same with mshift and loopoffmax in if statement below
for n=1:nwin
    if xcmaxAVEn(n)<xcmaxAVEnmin || abs(xmaxPGSIn(n)-xmaxPGSSn(n)+xmaxSISSn(n))>loopoffmax ...
            || isequal(abs(imaxPGSS(n)),mshift) || isequal(abs(imaxPGSI(n)),mshift) || isequal(abs(imaxSISS(n)),mshift) ...
            || nzerosPGC(n)>zerosallowed || nzerosSSIB(n)>zerosallowed || nzerosSILB(n)>zerosallowed
        xmaxPGSSntmp(n)=srate; xmaxPGSIntmp(n)=srate; xmaxSISSntmp(n)=srate; 
%         if timswin(n) > 35204
%             [1 n timswin(n) xcmaxAVEn(n) xmaxPGSIn(n)-xmaxPGSSn(n)+xmaxSISSn(n)]
%             pause
%         end
    else
        istart=igstart100+(n-1)*winoff100+mshift100; 
        iend=istart+winlen100-1;
        winPGC=realPGC100(istart:iend); winSSIB=realSSIB100(istart:iend); winSILB=realSILB100(istart:iend);
        xcPGSSn=xcorr(winSSIB,winPGC,mshift100,'coeff');
        xcPGSIn=xcorr(winSILB,winPGC,mshift100,'coeff');
        xcSISSn=xcorr(winSSIB,winSILB,mshift100,'coeff');
%         interpPGSSn=interp(sumsPGSSn(n,:),iup,3);
%         interpPGSIn=interp(sumsPGSIn(n,:),iup,3);
%         interpSISSn=interp(sumsSISSn(n,:),iup,3);
%         leninterp=length(interpPGSSn);
        [xcmaxPGSSn100,imaxPGSS100]=max(xcPGSSn);
        [xcmaxPGSIn100,imaxPGSI100]=max(xcPGSIn);
        [xcmaxSISSn100,imaxSISS100]=max(xcSISSn);
        xcmaxconprev=-99999.;  %used to be 0; not good with glitches
        for iPGSS=max(1,imaxPGSS100-7):min(imaxPGSS100+7,2*mshift100+1) %7 samples from peak; 
            for iPGSI=max(1,imaxPGSI100-7):min(imaxPGSI100+7,2*mshift100+1)
                ibangon = (mshift100+1)-iPGSI+iPGSS;
                if ibangon >= 1 && ibangon<=2*mshift100+1
                    xcmaxcon=xcPGSSn(iPGSS)+xcPGSIn(iPGSI)+xcSISSn(ibangon);
                    if xcmaxcon > xcmaxconprev
                        xcmaxconprev=xcmaxcon;
                        iPGSSbang=iPGSS;
                        iPGSIbang=iPGSI;
                    end
                end
            end
        end
        iSISSbang=mshift100+1-iPGSIbang+iPGSSbang;
%         if abs(timswin(n)-55158) <= 0.4
%             [2 timswin(n) xcmaxAVEn(n) xmaxPGSIn(n)-xmaxPGSSn(n)+xmaxSISSn(n) ...
%             iPGSSbang imaxinterpPGSS xcmaxinterpPGSSn iPGSIbang imaxinterpPGSI xcmaxinterpPGSIn ...
%             iSISSbang imaxinterpSISS xcmaxinterpSISSn]
%             pause
%         end
        xcmaxAVEn100=(xcPGSSn(iPGSSbang)+xcPGSIn(iPGSIbang)+xcSISSn(iSISSbang))/3;
        if abs(iPGSSbang-imaxPGSS100) <= loopoffmax100 && ...
           abs(iPGSIbang-imaxPGSI100) <= loopoffmax100 && ...
           abs(iSISSbang-imaxSISS100) <= loopoffmax100 && ...
           xcmaxAVEn100 >= xcmaxAVEnmin
            xmaxPGSSntmp(n)=-(iPGSSbang-(mshift100+1));  %empirical -, to be consistent w/previoue convention (whatever that was)
            xmaxPGSIntmp(n)=-(iPGSIbang-(mshift100+1));
            xmaxSISSntmp(n)=-(iSISSbang-(mshift100+1));
            
            %for plotting traces
            imaxPGSSwr=round(xmaxPGSSntmp(n));
            imaxPGSIwr=round(xmaxPGSIntmp(n));

            
            imid=round((istart+iend)/2);
            %Check power spectrum
            [PGCxx fp] = pwelch(winPGC,[],[],[],srate); %40 is sps
            PGCxx=PGCxx/max(PGCxx);
            [SSIBxx fp] = pwelch(realSSIB100(istart-imaxPGSSwr:iend-imaxPGSSwr),[],[],[],srate);
            SSIBxx=SSIBxx/max(SSIBxx);
            [SILBxx fp] = pwelch(realSILB100(istart-imaxPGSIwr:iend-imaxPGSIwr),[],[],[],srate);
            SILBxx=SILBxx/max(SILBxx);
            flo=find(fp > lo,1)-1;
            fhi=find(fp > hi,1)+1; %extra 1 for good measure
            belowcut=median([PGCxx(2:flo); SSIBxx(2:flo); SILBxx(2:flo)]);
            ppeaksPGC=findpeaks(PGCxx(flo+1:fhi));
            if length(ppeaksPGC)>=1
                maxppeakPGC=max(ppeaksPGC);
            else
                maxppeakPGC=0.;
            end
            ppeaksSSIB=findpeaks(SSIBxx(flo+1:fhi));
            if length(ppeaksSSIB)>=1
                maxppeakSSIB=max(ppeaksSSIB);
            else
                maxppeakSSIB=0.;
            end
            ppeaksSILB=findpeaks(SILBxx(flo+1:fhi));
            if length(ppeaksSILB)>=1
                maxppeakSILB=max(ppeaksSILB);
            else
                maxppeakSILB=0.;
            end
            abovecut=median([maxppeakPGC maxppeakSSIB maxppeakSILB]);
            if abovecut > 0.9*belowcut %-1
                PGSStr=winPGC.*realSSIB100(istart-imaxPGSSwr:iend-imaxPGSSwr);
                PGSItr=winPGC.*realSILB100(istart-imaxPGSIwr:iend-imaxPGSIwr);
                SISStr=realSILB100(istart-imaxPGSIwr:iend-imaxPGSIwr).*realSSIB100(istart-imaxPGSSwr:iend-imaxPGSSwr);
                cumsumtr=cumsum(PGSStr)+cumsum(PGSItr)+cumsum(SISStr);
                [cumsumtrdiff idiff]=max(cumsumtr(srate+1:winlen100)-cumsumtr(1:winlen100-srate));

                PGCfile(nin*winlen100+1:(nin+1)*winlen100,1:2)=[timsSSIB(istart:iend)' winPGC];
                SSIBfile(nin*winlen100+1:(nin+1)*winlen100,1:2)=[timsSSIB(istart:iend)' realSSIB100(istart-imaxPGSSwr:iend-imaxPGSSwr)];
                SILBfile(nin*winlen100+1:(nin+1)*winlen100,1:2)=[timsSSIB(istart:iend)' realSILB100(istart-imaxPGSIwr:iend-imaxPGSIwr)];
                PGSSfile(nin+1,1:2)=[imaxPGSSwr xcmaxPGSSn(n)];
                PGSIfile(nin+1,1:2)=[imaxPGSIwr xcmaxPGSIn(n)];
                SISSfile(nin+1,1:3)=[cumsumtrdiff/cumsumtr(winlen100) xcmaxSISSn(n) idiff];

                nin=nin+1;
                aPGSSkeep(nin,:)=[timswin(n) aPGSS(n)];
                aPGSIkeep(nin,:)=[timswin(n) aPGSI(n)];
                aSISSkeep(nin,:)=[timswin(n) aSISS(n)];
                loopoffkeep(nin,:)=[timswin(n) loopoff(n)];
                mapfile(nin,:)=[timswin(n) xmaxPGSIntmp(n) xmaxPGSSntmp(n) ...  %this now has the ave CONSISTENT xcorr value
                    xcmaxAVEn100 loopoff(n) AmpComp(n) cumsumtrdiff timswin(n)-2+idiff/srate cumsumtrdiff/cumsumtr(winlen100)];
            else
                xmaxPGSSntmp(n)=srate; xmaxPGSIntmp(n)=srate; xmaxSISSntmp(n)=srate;
            end
%             if timswin(n) > 35204
%                 [2 n timswin(n) xcmaxAVEn100]
%                 pause
%             end
        else
            xmaxPGSSntmp(n)=srate; xmaxPGSIntmp(n)=srate; xmaxSISSntmp(n)=srate; 
%             if timswin(n) > 35204
%                 [3 n timswin(n) xcmaxAVEn100]
%                 pause
%             end
        end
    end
end
%fid = fopen(['ARMMAP/MAPS/map',IDENTIF,'_',int2str(lo),'-',int2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
fid = fopen(['ARMMAP/100/MAPS/map',IDENTIF,'_',int2str(lo),'-',int2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
fprintf(fid,'%9.1f %9.5f %9.5f %9.4f %9.2f %9.3f %10.4f %10.4f %8.3f\n',mapfile(1:nin,:)');
fclose(fid);
figure 
subplot(4,1,1,'align'); 
hold on
plot(timswin,xcmaxAVEnmin*mshift100+zeros(nwin,1),'k:');
plot(timsPGC(winlen:2*winlen),7+zeros(winlen+1,1),'k','linewidth',2);
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*mshift100,'g');
plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
axis([0 timbig/2 -mshift100 mshift100]);
ylabel('samples')
box on
subplot(4,1,2,'align'); 
hold on
plot(timswin,xcmaxAVEnmin*mshift100+zeros(nwin,1),'k:');
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*mshift100,'g');
plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
axis([timbig/2 timbig -mshift100 mshift100]);
ylabel('samples')
box on
subplot(4,1,3,'align'); 
hold on
plot(timswin,xcmaxAVEnmin*mshift100+zeros(nwin,1),'k:');
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*mshift100,'g');
plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
axis([timbig 3*timbig/2 -mshift100 mshift100]);
ylabel('samples')
box on
subplot(4,1,4,'align'); 
hold on
plot(timswin,xcmaxAVEnmin*mshift100+zeros(nwin,1),'k:');
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn*mshift100,'g');
plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
axis([3*timbig/2 2*timbig -mshift100 mshift100]);
xlabel('sec')
ylabel('samples')
box on
title([IDENTIF,'_{',int2str(lo),'-',int2str(hi),'}'])
orient landscape
print('-depsc',['ARMMAP/100/FIGS/',IDENTIF,'-',int2str(winlen/40),'s_',int2str(lo),'-',int2str(hi),'b.eps'])

for n=1:nwin;
    if abs(loopoff(n))>2.0  %I think loopoff is still the old (40sps) version
        loopoff(n)=srate; 
    end
end
figure 
subplot(4,1,1,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
plot(timsPGC(winlen:2*winlen),7+zeros(winlen+1,1),'k','linewidth',2);
axis([0 timbig/2 -mshift100 mshift100]);
title('blue = PGSS ;   red = PGSI ;  black = SISS')
box on
subplot(4,1,2,'align'); 
hold on
hrf = plotreflinesr(gca,detects,'x','k'); 
hrf = plotreflinesr(gca,Detects,'x','r'); 
hrf = plotreflinesr(gca,detects2_8,'x','g'); 
hrf = plotreflinesr(gca,Detects2_8,'x','b'); 
plot(timswin,xcmaxAVEnmin+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn,'k');
plot(timswin,-1+0.1*abs(loopoff),'g.','MarkerSize',1);
plot(timsPGC,min(0,realPGC),'r');
plot(timsPGC,min(0,realSSIB),'b');
plot(timsPGC,min(0,realSILB),'k');
axis([0 timbig/2 -1 1]);
box on
subplot(4,1,3,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
axis([timbig/2 timbig -mshift100 mshift100]);
box on
subplot(4,1,4,'align'); 
hold on
hrf = plotreflinesr(gca,detects,'x','k'); 
hrf = plotreflinesr(gca,Detects,'x','r'); 
hrf = plotreflinesr(gca,detects2_8,'x','g'); 
hrf = plotreflinesr(gca,Detects2_8,'x','b'); 
plot(timswin,xcmaxAVEnmin+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn,'k');
plot(timswin,-1+0.1*abs(loopoff),'g.','MarkerSize',1);
plot(timsPGC,min(0,realPGC),'r');
plot(timsPGC,min(0,realSSIB),'b');
plot(timsPGC,min(0,realSILB),'k');
axis([timbig/2 timbig -1 1]);
box on
title([IDENTIF,'_{',int2str(lo),'-',int2str(hi),'}'])
orient landscape
print('-depsc',['ARMMAP/100/FIGS/',IDENTIF,'-',int2str(winlen/40),'s_',int2str(lo),'-',int2str(hi),'c.eps'])

figure
subplot(4,1,1,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
plot(timsPGC(tracelenPGC/2+winlen:tracelenPGC/2+2*winlen),7+zeros(winlen+1,1),'k','linewidth',2);
axis([timbig 3*timbig/2 -mshift100 mshift100]);
box on
subplot(4,1,2,'align'); 
hold on
hrf = plotreflinesr(gca,detects,'x','k'); 
hrf = plotreflinesr(gca,Detects,'x','r'); 
hrf = plotreflinesr(gca,detects2_8,'x','g'); 
hrf = plotreflinesr(gca,Detects2_8,'x','b'); 
plot(timswin,xcmaxAVEnmin+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn,'k');
plot(timswin,-1+0.1*abs(loopoff),'g.','MarkerSize',1);
plot(timsPGC,min(0,realPGC),'r');
plot(timsPGC,min(0,realSSIB),'b');
plot(timsPGC,min(0,realSILB),'k');
axis([timbig 3*timbig/2 -1 1]);
box on
subplot(4,1,3,'align'); 
hold on
plot(timswin,zeros(nwin,1),'k:');
plot(timswin,xmaxPGSSntmp,'bs','MarkerSize',2);
plot(timswin,xmaxPGSIntmp,'ro','MarkerSize',2);
plot(timswin,xmaxSISSntmp,'k*','MarkerSize',2);
axis([3*timbig/2 2*timbig -mshift100 mshift100]);
box on
subplot(4,1,4,'align'); 
hold on
hrf = plotreflinesr(gca,detects,'x','k'); 
hrf = plotreflinesr(gca,Detects,'x','r'); 
hrf = plotreflinesr(gca,detects2_8,'x','g'); 
hrf = plotreflinesr(gca,Detects2_8,'x','b'); 
plot(timswin,xcmaxAVEnmin+zeros(nwin,1),'k:');
plot(timswin,xcmaxAVEn,'k');
plot(timswin,-1+0.1*abs(loopoff),'g.','MarkerSize',1);
plot(timsPGC,min(0,realPGC),'r');
plot(timsPGC,min(0,realSSIB),'b');
plot(timsPGC,min(0,realSILB),'k');
axis([3*timbig/2 2*timbig -1 1]);
box on
title([IDENTIF,'_{',int2str(lo),'-',int2str(hi),'}'])
orient landscape
print('-depsc',['ARMMAP/100/FIGS/',IDENTIF,'-',int2str(winlen/40),'s_',int2str(lo),'-',int2str(hi),'d.eps'])

figure
colormap(jet)
scatter(xmaxPGSIn-xmaxPGSSn+xmaxSISSn,xcmaxAVEn,3,AmpComp)
hold on 
plot(-5:5,xcmaxAVEnmin+zeros(11,1),'k:');
axis([-5 5 -0.2 1.0])
hrf = plotreflinesr(gca,-1.5,'x','k');colorbar
hrf = plotreflinesr(gca,1.5,'x','k');colorbar
box on
title([IDENTIF,'_{',int2str(lo),'-',int2str(hi),'}'])
print('-depsc',['ARMMAP/100/FIGS/',IDENTIF,'-',int2str(winlen/40),'s_',int2str(lo),'-',int2str(hi),'e.eps'])

if winlen<=200
scrsz=get(0,'ScreenSize');
nt=0;
nrow=4;
mcol=6;
for ifig=1:floor(nin/(nrow*mcol))+1
    figure('Position',[scrsz(3)/10 scrsz(4)/10 4*scrsz(3)/5 9*scrsz(4)/10]);
    for n = 1:nrow
        for m = 1:mcol
            nt=nt+1;
            if nt <= nin
    %             if PGSSfile(nt,1) <= 0 && PGSSfile(nt,1) >= -2 && abs(PGSIfile(nt,1))<=2
                    subplot(3*nrow,mcol,3*(n-1)*mcol+m,'align');
                    plot(PGCfile(winlen100*(nt-1)+1:winlen100*nt,1),PGCfile(winlen100*(nt-1)+1:winlen100*nt,2),'r')
                    hold on
                    plot(SSIBfile(winlen100*(nt-1)+1:winlen100*nt,1),SSIBfile(winlen100*(nt-1)+1:winlen100*nt,2),'b')
                    plot(SILBfile(winlen100*(nt-1)+1:winlen100*nt,1),SILBfile(winlen100*(nt-1)+1:winlen100*nt,2),'k')
                    is = PGCfile(winlen100*(nt-1)+1,1);
                    ien= PGCfile(winlen100*nt,1);
                    %yma=0.4;
                    yma=max(max([PGCfile(winlen100*(nt-1)+1:winlen100*nt,2) SSIBfile(winlen100*(nt-1)+1:winlen100*nt,2) ...
                        SILBfile(winlen100*(nt-1)+1:winlen100*nt,2)]));
                    ymi=min(min([PGCfile(winlen100*(nt-1)+1:winlen100*nt,2) SSIBfile(winlen100*(nt-1)+1:winlen100*nt,2) ...
                        SILBfile(winlen100*(nt-1)+1:winlen100*nt,2)]));
                    xvect=[is is+2*(yma-ymi)*(winlen100/(4*srate))]; %amplitude bar originally scaled for 4-s window
                    yma=2.4*max(yma,-ymi);
                    yvect=[-0.9*yma -0.9*yma];
                    plot(xvect,yvect,'r','linewidth',3)
                    plot([is is+1/lo],[-0.8*yma -0.8*yma],'k','linewidth',3)
                    text(is+0.2, 0.66*yma, int2str(PGSIfile(nt,1)),'fontsize',6);
                    text(ien-0.6, 0.66*yma, int2str(PGSSfile(nt,1)),'fontsize',6);
                    box on
                    axis([is ien -yma yma])
                    set(gca,'XTick',[is is+2],'fontsize',6);

                    subplot(3*nrow,mcol,3*(n-1)*mcol+mcol+m,'align');
%                     [PGCxx f] = pwelch(PGCfile(winlen*(nt-1)+1:winlen*nt,2)-mean(PGCfile(winlen*(nt-1)+1:winlen*nt,2)),[],[],[],40);
%                     [SSIBxx f] = pwelch(SSIBfile(winlen*(nt-1)+1:winlen*nt,2)-mean(SSIBfile(winlen*(nt-1)+1:winlen*nt,2)),[],[],[],40);
%                     [SILBxx f] = pwelch(SILBfile(winlen*(nt-1)+1:winlen*nt,2)-mean(SILBfile(winlen*(nt-1)+1:winlen*nt,2)),[],[],[],40);
                    [PGCxx f] = pwelch(PGCfile(winlen100*(nt-1)+1:winlen100*nt,2),[],[],[],srate);
                    [SSIBxx f] = pwelch(SSIBfile(winlen100*(nt-1)+1:winlen100*nt,2),[],[],[],srate);
                    [SILBxx f] = pwelch(SILBfile(winlen100*(nt-1)+1:winlen100*nt,2),[],[],[],srate);
%                     [Cpgss,f] = mscohere(PGCfile(winlen*(nt-1)+1:winlen*nt,2),SSIBfile(winlen*(nt-1)+1:winlen*nt,2), ...
%                         [],[],[],40);
%                     [Cpgsi,f] = mscohere(PGCfile(winlen*(nt-1)+1:winlen*nt,2),SILBfile(winlen*(nt-1)+1:winlen*nt,2), ...
%                         [],[],[],40);
%                     [Csiss,f] = mscohere(SILBfile(winlen*(nt-1)+1:winlen*nt,2),SSIBfile(winlen*(nt-1)+1:winlen*nt,2), ...
%                         [],[],[],40);
%                     plot(f,Cpgss,'b')
%                     hold on
%                     plot(f,Cpgsi,'r')
%                     plot(f,Csiss,'k')
                    plot(f,PGCxx/max(PGCxx),'r+')
                    hold on
                    plot(f,SSIBxx/max(SSIBxx),'b+')
                    plot(f,SILBxx/max(SILBxx),'k+')
                    PGSStr=PGCfile(winlen100*(nt-1)+1:winlen100*nt,2).*SSIBfile(winlen100*(nt-1)+1:winlen100*nt,2);            
                    PGSItr=PGCfile(winlen100*(nt-1)+1:winlen100*nt,2).*SILBfile(winlen100*(nt-1)+1:winlen100*nt,2);
                    SISStr=SSIBfile(winlen100*(nt-1)+1:winlen100*nt,2).*SILBfile(winlen100*(nt-1)+1:winlen100*nt,2);
                    % find the peaks etc.
                    avedots=(PGSStr+PGSItr+SISStr)/3.;
                    [peaks, locs]=findpeaks(avedots,'minpeakdistance',3);
                    npeaks=length(peaks);
                    [maxpk, imaxpk]=max(peaks);
                    pks=zeros(npeaks,2);
                    pks(:,2)=peaks;
                    pks(:,1)=PGCfile(winlen100*(nt-1)+locs);
                    pksort=sortrows(pks,2);
                    rat14=maxpk/pksort(npeaks-4,2); %ratio of max to 4th largest anywhere in window
                    if imaxpk==1 
                        maxpkp=peaks(2);
                        maxpkm=-9e9
                        pkwid=2*(pks(2,1)-pks(1,1));
                        pksid12=maxpk/maxpkp;
                        pksid13=maxpk/peaks(3);
                    elseif imaxpk==npeaks
                        maxpkp=-9e9
                        maxpkm=peaks(npeaks-1);
                        pkwid=2*(pks(npeaks,1)-pks(npeaks-1,1));
                        pksid12=maxpk/maxpkm;
                        pksid13=maxpk/peaks(npeaks-2);
                    else
                        maxpkp=peaks(imaxpk+1);
                        maxpkm=peaks(imaxpk-1);
                        pkwid=pks(imaxpk+1,1)-pks(imaxpk-1,1);
                        pksid12=maxpk/max(maxpkp,maxpkm);
                        pksid13=maxpk/min(maxpkp,maxpkm);
                    end
                    cumsumPGSStr=cumsum(PGSStr);
                    cumsumPGSItr=cumsum(PGSItr);
                    cumsumSISStr=cumsum(SISStr);
%                     yma=1.1; %yma=max(avedots);
%                     ymi=-0.1; %ymi=min(avedots);
                    %hold on
%                     plot(PGCfile(winlen*(nt-1)+1:winlen*nt,1),PGSStr,'b')
%                     plot(PGCfile(winlen*(nt-1)+1:winlen*nt,1),PGSItr,'r')
%                     plot(PGCfile(winlen*(nt-1)+1:winlen*nt,1),SISStr,'k')
%                     plot(PGCfile(winlen*(nt-1)+1:winlen*nt,1),avedots,'c') %,'linewidth',2)
%                     plot(PGCfile(winlen*(nt-1)+1:winlen*nt,1),zeros(winlen,1),'k')
%                     plot(pks(:,1),pks(:,2),'r+')
                    box on
                    text(is+0.1, ymi+0.82*(yma-ymi), num2str(PGSIfile(nt,2),3),'fontsize',6);
                    text(is+0.1, ymi+0.64*(yma-ymi), num2str(PGSSfile(nt,2),3),'fontsize',6);
                    text(is+0.1, ymi+0.46*(yma-ymi), num2str(SISSfile(nt,2),3),'fontsize',6);
                    text(ien-0.8, ymi+0.04*(yma-ymi), num2str(SISSfile(nt,1),3),'fontsize',6);
                    %axis([lo-1 hi+1 ymi yma])
                    xlim([0 hi+2])
                    set(gca,'XTick',(0:20),'fontsize',6);
                    %axis tight
                    %set(gca,'XTick',[is is+2],'fontsize',6);
                    pkfile(nt,:)=[PGCfile(winlen100*(nt-1)+1+(winlen100/2)) pks(imaxpk,1) maxpk pksid12 pksid13 pkwid rat14];

                    PGCauto=PGCfile(winlen100*(nt-1)+1:winlen100*nt,2).*PGCfile(winlen100*(nt-1)+1:winlen100*nt,2);
                    PGC2=cumsum(PGCauto);
                    SSIBauto=SSIBfile(winlen100*(nt-1)+1:winlen100*nt,2).*SSIBfile(winlen100*(nt-1)+1:winlen100*nt,2);
                    SSIB2=cumsum(SSIBauto);
                    SILBauto=SILBfile(winlen100*(nt-1)+1:winlen100*nt,2).*SILBfile(winlen100*(nt-1)+1:winlen100*nt,2);
                    SILB2=cumsum(SILBauto);
                    PGSSnum=cumsumPGSStr((srate/2+1):winlen100)-cumsumPGSStr(1:winlen100-srate/2);
                    PGSInum=cumsumPGSItr((srate/2+1):winlen100)-cumsumPGSItr(1:winlen100-srate/2);
                    SISSnum=cumsumSISStr((srate/2+1):winlen100)-cumsumSISStr(1:winlen100-srate/2);
                    PGCden=PGC2((srate/2+1):winlen100)-PGC2(1:winlen100-srate/2);
                    SSIBden=SSIB2((srate/2+1):winlen100)-SSIB2(1:winlen100-srate/2);
                    SILBden=SILB2((srate/2+1):winlen100)-SILB2(1:winlen100-srate/2);
                    subplot(3*nrow,mcol,3*(n-1)*mcol+2*mcol+m,'align');
                    PGSSn=PGSSnum./realsqrt(PGCden.*SSIBden);
                    PGSIn=PGSInum./realsqrt(PGCden.*SILBden);
                    SISSn=SISSnum./realsqrt(SSIBden.*SILBden);
                    alln=(PGSSn+PGSIn+SISSn)/3;
                    idiff=SISSfile(nt,3);
                    maxxc=max(alln(idiff:idiff+srate/2)); %endpoints of alln window should be endpoints of 1 sec cumsum interval
                    plot(PGCfile(winlen100*(nt-1)+(srate/4+1):winlen100*nt-(srate/4),1),alln,'c','linewidth',3)
                    hold on
                    plot(PGCfile(winlen100*(nt-1)+(srate/4+1):winlen100*nt-(srate/4),1),PGSSn,'b')
                    plot(PGCfile(winlen100*(nt-1)+(srate/4+1):winlen100*nt-(srate/4),1),PGSIn,'r')
                    plot(PGCfile(winlen100*(nt-1)+(srate/4+1):winlen100*nt-(srate/4),1),SISSn,'k')
                    plot(PGCfile(winlen100*(nt-1)+1:winlen100*nt,1),0.75*maxxc*ones(winlen100,1),'k:');
                    plot(PGCfile(winlen100*(nt-1)+1:winlen100*nt,1),0.65*maxxc*ones(winlen100,1),'k:');
                    axis([is ien -0.5 1])
                    set(gca,'XTick',[is is+2],'fontsize',6);
    %             end
            end
        end
    end
    orient landscape
    if ifig <= 9
        print('-depsc',['ARMMAP/100/WIGS/',IDENTIF,'-',int2str(lo),'-',int2str(hi),'.',int2str(0),int2str(0),int2str(ifig),'.eps'])
    elseif ifig <= 99
        print('-depsc',['ARMMAP/100/WIGS/',IDENTIF,'-',int2str(lo),'-',int2str(hi),'.',int2str(0),int2str(ifig),'.eps'])
    else
        print('-depsc',['ARMMAP/100/WIGS/',IDENTIF,'-',int2str(lo),'-',int2str(hi),'.',int2str(ifig),'.eps'])
    end

end
%fid = fopen(['ARMMAP/MAPS/pks',IDENTIF,'_',int2str(lo),'-',int2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
fid = fopen(['ARMMAP/100/pks',IDENTIF,'_',int2str(lo),'-',int2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
fprintf(fid,'%9.1f %9.3f %10.6f %9.3f %9.3f %9.3f %9.3f \n',pkfile(1:nin,:)');
fclose(fid);
end

medlok=median(abs(loopoffkeep))
medaPGSS=median(aPGSSkeep)
medaPGSI=median(aPGSIkeep)
medaSISS=median(aSISSkeep)


%cputime-t;
tot=cputime-tt
end
