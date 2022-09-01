%Gets Hilbert transforms of several days of TWKB seismograms; writes out sampe every 2 seconds.  ALTERED SINCE LAST RUN - CHECK!
%Also gets local max and min of each seismogram.  NOT YET!
format short e
clear all
%set(0,'DefaultFigureVisible','off');
%      timoffrot(1,:)=[2005 254 07 42 12 +086 +020  80 115  50];
%      timoffrot(2,:)=[2005 255 07 42 12 +086 +020  80 115  50];
%      timoffrot(3,:)=[2005 256 07 42 12 +086 +020  80 115  50];
%      timoffrot(4,:)=[2005 257 07 42 12 +086 +020  80 115  50];
%      timoffrot(5,:)=[2005 258 07 42 12 +086 +020  80 115  50];
%      timoffrot(6,:)=[2005 259 07 42 12 +086 +020  80 115  50];
%      timoffrot(7,:)=[2005 260 07 42 12 +086 +020  80 115  50];
%      timoffrot(8,:)=[2005 261 07 42 12 +086 +020  80 115  50];

%     timoffrot(1,:)=[2004 194 07 42 12 +086 +020  80 115  50];
%     timoffrot(2,:)=[2004 195 07 42 12 +086 +020  80 115  50];
%     timoffrot(3,:)=[2004 196 07 42 12 +086 +020  80 115  50];
%     timoffrot(4,:)=[2004 197 07 42 12 +086 +020  80 115  50];
%     timoffrot(5,:)=[2004 198 07 42 12 +086 +020  80 115  50];
%     timoffrot(6,:)=[2004 199 07 42 12 +086 +020  80 115  50];
%     timoffrot(7,:)=[2004 200 07 42 12 +086 +020  80 115  50];
%     timoffrot(8,:)=[2004 201 07 42 12 +086 +020  80 115  50];
%     timoffrot(9,:)=[2004 202 07 42 12 +086 +020  80 115  50];
%     timoffrot(10,:)=[2004 203 07 42 12 +086 +020  80 115  50];

   timoffrot(1,:)=[2003 060 07 42 12 +086 +020  80 115  50];
   timoffrot(2,:)=[2003 061 07 42 12 +086 +020  80 115  50];
   timoffrot(3,:)=[2003 062 07 42 12 +086 +020  80 115  50];
   timoffrot(4,:)=[2003 063 07 42 12 +086 +020  80 115  50];
   timoffrot(5,:)=[2003 064 07 42 12 +086 +020  80 115  50];
   timoffrot(6,:)=[2003 065 07 42 12 +086 +020  80 115  50];
   timoffrot(7,:)=[2003 066 07 42 12 +086 +020  80 115  50];
   timoffrot(8,:)=[2003 067 07 42 12 +086 +020  80 115  50];

ndays=length(timoffrot(:,1))
samperhr=0.25/3600  %assuming 40 sps; printout in hours
for nd=1:length(timoffrot(:,1))
close all
 hi=6.;
 lo=1.5;
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
%rotPGC=pi*(timoffrot(nd,8)-90)/180;
%rotTWKB=pi*(timoffrot(nd,9)-90)/180;
%rotSILB=pi*(timoffrot(nd,10)-90)/180;
%TWKBsoff=timoffrot(nd,6);
%SILBsoff=timoffrot(nd,7);
%TWKBtoff=TWKBsoff/40;
%SILBtoff=SILBsoff/40;
IDENTIF=[YEAR,'.',JDAY,'.',int2str(timoffrot(nd,6)),'.',int2str(timoffrot(nd,7)), ...
    '.',int2str(timoffrot(nd,8)),'.',int2str(timoffrot(nd,9)),'.',int2str(timoffrot(nd,10))]

%Read data:
direc=[YEAR,'/',MO,'/'];
prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];
%PGCEdat=[prename,'.LZB..BHE.D.SAC'];
%PGCNdat=[prename,'.LZB..BHN.D.SAC'];
TWKBEdat=[prename,'.TWKB..HHE.D.SAC'];
TWKBNdat=[prename,'.TWKB..HHN.D.SAC'];
%SILBEdat=[prename,'.MGCB..HHE.D.SAC'];
%SILBNdat=[prename,'.MGCB..HHN.D.SAC'];

%[PGCE,HdrDataPGC,tnuPGC,pobjPGC,timsPGC]=readsac(PGCEdat,0,'l');
%[PGCN,~,~,~,~]=readsac(PGCNdat,0,'l');
[TWKBE,HdrDataTWKB,tnuTWKB,pobjTWKB,timsTWKB]=readsac(TWKBEdat,0,'l');
[TWKBN,~,~,~,~]=readsac(TWKBNdat,0,'l');
%[SILBE,HdrDataSILB,tnuSILB,pobjSILB,timsSILB]=readsac(SILBEdat,0,'l');
%[SILBN,~,~,~,~]=readsac(SILBNdat,0,'l');
%
%tracelenPGC=length(PGCE);
tracelenTWKB=length(TWKBE);
%tracelenSILB=length(SILBE);

%    timPGCfirst=timsPGC(1);
%    timSILBfirst=timsSILB(1); %0.01s earlier than PGC
%    timPGClast=timsPGC(tracelenPGC);
%    timSILBlast=timsSILB(tracelenSILB);

%%cosine taper before filtering:
%x=(0:pi/80:pi/2-pi/80)';
%% %Seems to be necessary at the start of each day for PGCE:
%    PGCE(1:80)=0.;
%    PGCE(81:120)=sin(x).*PGCE(81:120); %Only at start of day!
%%PGCE(1:40)=sin(x).*PGCE(1:40);
%PGCN(1:40)=sin(x).*PGCN(1:40);
%x=flipud(x);
%PGCE(tracelenPGC-39:tracelenPGC)=sin(x).*PGCE(tracelenPGC-39:tracelenPGC);
%PGCN(tracelenPGC-39:tracelenPGC)=sin(x).*PGCN(tracelenPGC-39:tracelenPGC);
x=(0:pi/200:pi/2-pi/200)';
%PGCE(1:100)=sin(x).*PGCE(1:100);
%PGCN(1:100)=sin(x).*PGCN(1:100);
TWKBE(1:100)=sin(x).*TWKBE(1:100);
TWKBN(1:100)=sin(x).*TWKBN(1:100);
%SILBE(1:100)=sin(x).*SILBE(1:100);
%SILBN(1:100)=sin(x).*SILBN(1:100);
x=flipud(x);
%PGCE(tracelenPGC-99:tracelenPGC)=sin(x).*PGCE(tracelenPGC-99:tracelenPGC);
%PGCN(tracelenPGC-99:tracelenPGC)=sin(x).*PGCN(tracelenPGC-99:tracelenPGC);
TWKBE(tracelenTWKB-99:tracelenTWKB)=sin(x).*TWKBE(tracelenTWKB-99:tracelenTWKB);
TWKBN(tracelenTWKB-99:tracelenTWKB)=sin(x).*TWKBN(tracelenTWKB-99:tracelenTWKB);
%SILBE(tracelenSILB-99:tracelenTWKB)=sin(x).*SILBE(tracelenSILB-99:tracelenSILB);
%SILBN(tracelenSILB-99:tracelenTWKB)=sin(x).*SILBN(tracelenSILB-99:tracelenSILB);

%Looks like filtering the 40-Hz data followed by interpolation introduces 
%less high-frequency stuff that does not appear in the original filtered

%Filter data:
npo=2;
npa=1;
%[PGCEf]=1.6e-4*bandpass(PGCE,40,lo,hi,npo,npa,'butter');
%[PGCNf]=1.6e-4*bandpass(PGCN,40,lo,hi,npo,npa,'butter');
%[PGCZf]=1.6e-4*bandpass(PGCZ,40,lo,hi,npo,npa,'butter');
if year==2003 && jday<213
%    [PGCEf]=20.0e-3*bandpass(PGCE,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
%    [PGCNf]=20.0e-3*bandpass(PGCN,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
    [TWKBEf]=20.0e-3*bandpass(TWKBE,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
    [TWKBNf]=20.0e-3*bandpass(TWKBN,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
%    [SILBEf]=20.0e-3*bandpass(SILBE,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0* 
%    [SILBNf]=20.0e-3*bandpass(SILBN,100,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
else
    [TWKBEf]=4.0e-3*bandpass(TWKBE,100,lo,hi,npo,npa,'butter'); 
    [TWKBNf]=4.0e-3*bandpass(TWKBN,100,lo,hi,npo,npa,'butter'); 
%    [SILBEf]=4.0e-3*bandpass(SILBE,100,lo,hi,npo,npa,'butter'); 
%    [SILBNf]=4.0e-3*bandpass(SILBN,100,lo,hi,npo,npa,'butter'); 
end

%Resample the 100 sps data:
%PGCEfd = resample(PGCEf,2,5);
%PGCNfd = resample(PGCNf,2,5);
TWKBEfd = resample(TWKBEf,2,5);
TWKBNfd = resample(TWKBNf,2,5);
%SILBEfd = resample(SILBEf,2,5);
%SILBNfd = resample(SILBNf,2,5);
%     timTWKBfirst=timsTWKB(1); %0.01s earlier than PGC
%     timTWKBlast=timsTWKB(tracelenTWKB);
    
%PGC=PGCEfd+1i*PGCNfd;
TWKB=TWKBEfd+1i*TWKBNfd;
%SILB=SILBEfd+1i*SILBNfd;
%PGCrot=PGC*exp(1i*rotPGC);
TWKBrot=abs(TWKB);
%SILBrot=SILB*exp(1i*rotSILB);
    %rtate=cputime-t
    %t=cputime;

%realPGC=real(PGCrot);
%realTWKB=real(TWKBrot);
%realSILB=real(SILBrot);
%hilPGC=abs(hilbert(realPGC));
hilTWKB=abs(hilbert(TWKBrot));
hilTWKBd=resample(hilTWKB,1,80);
%hilSILB=abs(hilbert(realSILB));

len=length(hilTWKBd);
TWKBhilalldays(len*(nd-1)+1:len*nd)=hilTWKBd;

a=reshape(TWKBrot,80,43200);
TWKBmax=max(a);
len=length(TWKBmax);
TWKBmaxalldays(len*(nd-1)+1:len*nd)=TWKBmax;
% [TWKBpks,locs]=findpeaks(TWKBrot); %,'minpeakdistance',40);  %OLD
% if nd ==1
%     TWKBpksalldays=TWKBpks;
%     TWKBtimsalldays=locs*samperhr;
% else
%     len=length(TWKBpksalldays)
%     lennew=length(TWKBpks);
%     TWKBpksalldays(len+1:len+lennew)=TWKBpks;
%     TWKBtimsalldays(len+1:end+lennew)=24*(nd-1)+locs*samperhr;
% end

%[PGChilb,locs]=findpeaks(hilPGC);
%PGChilbtims=timsPGC(locs);
%[TWKBhilb,locs]=findpeaks(hilTWKB);
%TWKBhilbtims=timsPGC(locs);
%[SILBhilb,locs]=findpeaks(hilSILB);
%SILBhilbtims=timsPGC(locs);

% PGCpeaks=[PGChilbtims' PGChilb];
% TWKBpeaks=[TWKBhilbtims' TWKBhilb];
% SILBpeaks=[SILBhilbtims' SILBhilb];
% 
% timsPGC = resample(timsPGC,2,5);
% PGCtrace=[timsPGC' realPGC];
% TWKBtrace=[timsPGC' realTWKB];
% SILBtrace=[timsPGC' realSILB];
% 
% fid = fopen(['HILBERTS/LZB/PGC.',IDENTIF,'_',int2str(lo),'-',int2str(hi)],'w');
% fprintf(fid,'%13.5f %9.5f\n',PGCtrace(:,:)');
% fclose(fid);
% fid = fopen(['HILBERTS/LZB/TWKB.',IDENTIF,'_',int2str(lo),'-',int2str(hi)],'w');
% fprintf(fid,'%13.5f %9.5f\n',TWKBtrace(:,:)');
% fclose(fid);
% fid = fopen(['HILBERTS/LZB/SILB.',IDENTIF,'_',int2str(lo),'-',int2str(hi)],'w');
% fprintf(fid,'%13.5f %9.5f\n',SILBtrace(:,:)');
% fclose(fid);
% 
% fid = fopen(['HILBERTS/LZB/PGChp.',IDENTIF,'_',int2str(lo),'-',int2str(hi)],'w');
% fprintf(fid,'%13.5f %9.5f\n',PGCpeaks(:,:)');
% fclose(fid);
% fid = fopen(['HILBERTS/LZB/TWKBhp.',IDENTIF,'_',int2str(lo),'-',int2str(hi)],'w');
% fprintf(fid,'%13.5f %9.5f\n',TWKBpeaks(:,:)');
% fclose(fid);
% fid = fopen(['HILBERTS/LZB/SILBhp.',IDENTIF,'_',int2str(lo),'-',int2str(hi)],'w');
% fprintf(fid,'%13.5f %9.5f\n',SILBpeaks(:,:)');
% fclose(fid);

end
fid = fopen(['HILBERTS/TWKBhil.',YEAR,'_',int2str(timoffrot(1,2)),'-',int2str(timoffrot(ndays,2)),'_',int2str(lo),'-',int2str(hi)],'w');
tims=1:2:2*length(TWKBhilalldays);
fil(:,1)=tims;
fil(:,2)=TWKBhilalldays;
fprintf(fid,'%12i %9.5f\n',fil');
fclose(fid);

fid = fopen(['HILBERTS/TWKBpks.',YEAR,'_',int2str(timoffrot(1,2)),'-',int2str(timoffrot(ndays,2)),'_',int2str(lo),'-',int2str(hi)],'w');
tims=1:2:2*length(TWKBmaxalldays);
fil2(:,1)=tims;
fil2(:,2)=TWKBmaxalldays;
fprintf(fid,'%12i %9.5f\n',fil2');
fclose(fid);
