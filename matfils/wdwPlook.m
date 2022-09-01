% From wdwSTA; this one stretches and scales templates to match "strongest arrival".
% As of 12/20/17:  
% Writes to 2018 directories
% Re-introduces sub-sample time offsets, based on interpolated dot products
%    (NOT interpolated seisms).  
% writes time stretches as factors, not indices that must be remembered
%    when analyzing data.
% Includes an offset buffer of 0.5*concntr in the prior 1.25s data
% wdwP won't write mapfiles.  Just for wigs.  Comnents out stretched templates.
format short e
clear all
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');
% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo? 
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs
scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite; 

fam='065';
if isequal(fam,'002')  %DON'T FORGET additional family-specific delcarations around line 156
    %timoffrot= [2003 062];
    timoffrot= [2003 062; %Could leave out 064 and 257.
                2003 063;
    %           2003 064;
                2004 196;
                2004 197;
                2004 198;
                2004 199;
                2005 254;
                2005 255;
                2005 256];
    %           2005 257];
    %bostname='BOSTOCK/NEW/002-246_2003.062'; 
    bostname=['BOSTOCK/NEW/002-246_2003.062'; 
              'BOSTOCK/NEW/002-246_2003.063';
    %         'BOSTOCK/NEW/002-246_2003.064';
              'BOSTOCK/NEW/002-246_2004.196';
              'BOSTOCK/NEW/002-246_2004.197';
              'BOSTOCK/NEW/002-246_2004.198';
              'BOSTOCK/NEW/002-246_2004.199';
              'BOSTOCK/NEW/002-246_2005.254';
              'BOSTOCK/NEW/002-246_2005.255';
              'BOSTOCK/NEW/002-246_2005.256'];
    %         'BOSTOCK/NEW/002-246_2005.257'];
    PERMROTS=[0 90 32 0;  %PGC , Yajun's "0" changed to 90.  Don't think that matters, if no splitting.
              0 90 54 9]; %LZB
    POLROTS=[6 85 33 86;  %SSIB from Yajun
             0 90 39 20;  %SILB
             0 90  7 -4;  %KLNB
             4 70 48 -26; %MGCB
             4 75 38 -5]; %TWKB
    stas=['PGC  '
         'SSIB '
         'SILB '];  
    scaleseisms=[1.0 0.76 0.95]; %I think this is now important only for plotting purposes.  DIVIDES by scaleseisms.
    hi=6.5; 
    lo=1.25;
    npo=2;
    npa=2;
    mshift=19; %maximum shift for the x-correlations. 19 for 002 Stanford.  4 for 068?
    loopoffmax=1.5; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.
    xcmaxAVEnmin=0.44; %0.4 for 068? %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?    %tempoffs=[911 998 932]; %these are the OLD zero crossings for 002: PGC,SSIB,SILB.  87 or 86 samples for B's templates, it seems. 
    tempoffs=[101 101 101]; %these are the zero crossings for 002: PGC,SSIB,SILB.
    tempzeros=[93 109;
               92 109;
               92 108]; %These, inclusive, are the samples that get zeroed to make a single dipole.
    bostdelay=22.83-22.675; %002; 22.675 b.c. 002_ catalogs already "corrected" for PGC start; from tempseps.f. 22.83 a refinement.
elseif isequal(fam,'068')
    timoffrot=[2004 198; 
               2004 199;
               2004 200;
               2004 201;
               2005 256;
               2005 257;
               2005 258;
               2005 259;
               2005 260;
               2005 261];
    bostname=['BOSTOCK/NEW/068_2004.198'; 
              'BOSTOCK/NEW/068_2004.199';
              'BOSTOCK/NEW/068_2004.200';
              'BOSTOCK/NEW/068_2004.201';
              'BOSTOCK/NEW/068_2005.256';
              'BOSTOCK/NEW/068_2005.257';
              'BOSTOCK/NEW/068_2005.258';
              'BOSTOCK/NEW/068_2005.259';
              'BOSTOCK/NEW/068_2005.260';
              'BOSTOCK/NEW/068_2005.261'];
    PERMROTS=[0  0  0  0;  %PGC
              8 65  6 -8]; %LZB
    POLROTS =[0  0  0  0;  %SSIB 
              0  0  0  0;  %SILB
              2 50 25  0;  %KLNB (erroneous delay w.r.t. TWKB last column)
              5 65 41 -1;  %MGCB
              4 60 14  0]; %TWKB
    stas=['TWKB '
          'LZB  '
          'MGCB '];
    scaleseisms=[1.0 0.6 0.6]; %I think this is now important only for plotting purposes.  DIVIDES by scaleseisms.
    hi=8;
    lo=1.25;
    npo=2;
    npa=2;
    mshift=4; %maximum shift for the x-correlations. 19 for 002 Stanford.  4 for 068?
    loopoffmax=1.5; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.
    xcmaxAVEnmin=0.4; %0.4 for 068? Check it. %0.44 for 002 Stanford %0.45;  
    s%tempoffs=[911 998 932]; %these are the OLD zero crossings for 002: PGC,SSIB,SILB.  87 or 86 samples for B's templates, it seems. 
    tempoffs=[101 101 101]; %these are the zero crossings.
    tempzeros=[91 109;
               91 110;
               92 109]; %The last ones, and all those outside, are the samples that get zeroed to make a single dipole.
    bostdelay=22.625; %068; TWKB comes in at 905.
elseif isequal(fam,'065')
    % timoffrot=[2003 066;
    %            2003 067;
    %            2003 068;
    timoffrot=[2004 200];
%                2004 201;
%                2004 204;
%                2005 259;
%                2005 260];
    % bostname=['BOSTOCK/NEW/065_2003.066';
    %           'BOSTOCK/NEW/065_2003.067';
    %           'BOSTOCK/NEW/065_2003.068';
    bostname=['BOSTOCK/NEW/065_2004.200'];
%               'BOSTOCK/NEW/065_2004.201';
%               'BOSTOCK/NEW/065_2004.204';
%               'BOSTOCK/NEW/065_2005.259';
%               'BOSTOCK/NEW/065_2005.260'];
    PERMROTS=[0  0  0  0;  %PGC  4th column here is shift w.r.t. stas(1,:).  Best set to zero?
              0  0  0  0;  %LZB 
              0  0  0  0]; %LZBz
    POLROTS=[0   0   0   0;  %SSIB
             1 -70  65 103;  %SILB
             3  70  30   0;  %KLNB 
             0   0   0   0;  %MGCB 
             9  65  20 -16;  %TWKB
             0   0   0   0;  %SILBz
             0   0   0   0]; %SSIBz
    %        1  20 155   0;  %SILB
    %        3 160 120   0;  %KLNB 
    %        9 155 110   0;  %TWKB
    stas=['KLNB '
          'TWKB '
          'SILB '];
    stasP='SILB ';
    scaleseisms=[1.0 0.7 0.5]; %I think this is now important only for plotting purposes.  DIVIDES by scaleseisms.
    scalePseisms=[0.5]; %I think this is now important only for plotting purposes.  DIVIDES by scaleseisms.
    POLPROTS=[37.5 275 103 219]; %SILB; 4th col. is amount by which positive P peak preceeds positive S peak (samples)
    nstaP=size(stasP,1);
    hi=6.5;
    lo=1.25;
    npo=2;
    npa=2;
    mshift=5; %maximum shift for the x-correlations. 19 for 002 Stanford.  4 for 068?
    loopoffmax=1.5; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.
    xcmaxAVEnmin=0.4; %0.4 for 068? Check it. %0.44 for 002 Stanford %0.45;  
    %tempoffs=[911 998 932]; %these are the OLD zero crossings for 002: PGC,SSIB,SILB.  87 or 86 samples for B's templates, it seems. 
    tempoffs=[101 100 100]; %these are the zero crossings.
    tempzeros=[92 108;
               90 110;
               91 111]; %The last ones, and all those outside, are the samples that get zeroed to make a single dipole.
    bostdelay=23.250; %065; where KLNB comes in.
end

PERMSTA=['PGC'
         'LZB'];
POLSTA=['SSIB '
        'SILB '
        'KLNB '
        'MGCB '
        'TWKB '];
% Column order is:
%   1-Bostock family ID# 
%   2-fast/slow time offset (at 40 sps)
%   3-slow direction 
%   4-fast/slow cc coefficient 
%   5-corrected polarization angle
%   6-uncorrected polarization angle.  Angles are c-clockwise from East
% Station order is 1.PGC 2.LZB 3.VGZ 4.SSIB 5.SILB 6.TWKB 7.MGCB 8.TWBB 9.TSJB 10.KLNB
% 2.0000  0.0000  0.0000  0.0000  32.0000 22.8366
% 2.0000  0.0000  0.0000  0.0000  55.0000 17.7860    
  
nsta=size(stas,1);
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;
%POLROTS(:,1)=round(POLROTS(:,1)*(100/40)); %Polaris stations at 100 sps; table has split time at 40 sps.
%POLROTS(:,4)=round(POLROTS(:,4)*(100/40)); %Polaris stations at 100 sps; table has split time at 40 sps.
sps=100; %40;
tempwinlen=60*sps;
stack=zeros(nsta,tempwinlen);
stackort=zeros(nsta,tempwinlen);

%Basics of the cross-correlation:  Window length, number of windows, filter parameters, etc.
winlensec=4;
winoffsec=1;
winlen=winlensec*sps;
winoff=winoffsec*sps;
prewinsec=7;
prewin=prewinsec*sps;
wintot=winlen+prewin;
tracelen=86400*sps; %one day of data at 40 sps
winbig=2*(tracelen/2-(2*sps)); %ignore 2 seconds at each end of day
timbig=winbig/(2*sps); %half that time, in seconds
igstart=floor(tracelen/2-winbig/2)+1; %start counting seis data from here
nwin=floor((winbig-winlen)/winoff);
%UPGRADING SINCE MODIFYING READPOLS & READPERMS STOPED HERE
%hi=6.5;  %002 Stanford
%lo=1.25; %002 Stanford
% hi=6;
% lo=1.5;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Get the templates for the stations. Won't be always required.  PGCopt_002_1.25-6.5Hz_2pass_100sps.le14shift.resamp
% for ista=1:nsta
%     if stas(ista,4) == ' '
%         temptemps(ista,:)=load([stas(ista,1:3),'opt_',fam,'_',num2str(lo),'-',num2str(hi), ...
%                          'Hz_',int2str(npa),'pass_',int2str(sps),'sps.input'],'w');
%     else
%         temptemps(ista,:)=load([stas(ista,1:4),'opt_',fam,'_',num2str(lo),'-',num2str(hi), ...
%                          'Hz_',int2str(npa),'pass_',int2str(sps),'sps.input'],'w');
%     end
% end
% % tempbef=70;
% % tempaft=89;
% tempbef=59;
% tempaft=60;
% templen=tempbef+tempaft+1;
% for ista=1:nsta
%     STAtemps(ista,:)=temptemps(ista,tempoffs(ista)-tempbef:tempoffs(ista)+tempaft);
%     %snips(templen*(ista-1)+1:ista*templen)=STAtemps(ista,:);
% end
% %% scalefact scales templates; scaleseisms scales seisms.  Strategy changes with family.
% if isequal(fam,'002')
%     tempoffs=tempoffs-1; %Center of strongest window is 1 or 2 samples ahead of zero crossing (002); make it 1.
%     whichtoplot=2;
% elseif isequal(fam,'068')
%     whichtoplot=1;
% elseif isequal(fam,'065')
%     whichtoplot=1;
% end
% minses=-min(STAtemps,[],2); %STAtemps is (currently) 3 by 120
% maxses= max(STAtemps,[],2);
% plustominus=maxses./minses;
% scalefact=minses*max(plustominus); %This is used to scale templates, just for plotting purposes
% for ista=1:nsta
%     STAtemps(ista,:)=STAtemps(ista,:)/scalefact(ista); %This plots the templates with the largest positive value (of any) at +1
% end
% figure
% plot(STAtemps(1,:),'r')
% hold on
% plot(STAtemps(2,:),'b')
% plot(STAtemps(3,:),'k')
% drawnow
% % A different strategy for stretched templates?
% minstretch=12; %20
% maxstretch=80;
% stretchlen=60; %This is duration of window, in samples; not max-min
% STAstr=zeros(nsta,maxstretch-minstretch+1,stretchlen);
% ttemp=temptemps;
% minses=-min(ttemp,[],2); 
% maxses= max(ttemp,[],2);
% plustominus=maxses./minses;
% scalefact=minses*max(plustominus); %This is used to scale templates, just for plotting purposes
% for ista=1:nsta
%     ttemp(ista,:)=ttemp(ista,:)/scalefact(ista); %This plots the templates with the largest positive value (of any) at +1
% end
% figure
% hold on
% for ista=1:nsta
%     ttemp(ista,1:tempzeros(ista,1))=0;
%     ttemp(ista,tempzeros(ista,2):end)=0;
%     temp=ttemp(ista,tempoffs(ista)-stretchlen+1:tempoffs(ista)+stretchlen); %this has some slop at the ends 
%     plot(temp)
%     STAstr(ista,:,:)=stretchtemps(temp,minstretch,maxstretch,stretchlen,sps);
% end
% figure
% toplot=squeeze(STAstr(:,1,:));
% plot(toplot(1,:),'r')
% hold on
% plot(toplot(2,:),'b')
% plot(toplot(3,:),'k')
% figure
% toplot=squeeze(STAstr(:,end,:));
% plot(toplot(1,:),'r')
% hold on
% plot(toplot(2,:),'b')
% plot(toplot(3,:),'k')
% figure
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cycle over each day:
for nd=1:length(timoffrot(:,1))
    close all
    %Bostock's detections:
    bostocks=load(bostname(nd,:));
    bostsec=3600*(bostocks(:,3)-1)+bostocks(:,4)+bostdelay; 
    bostsamp=round(bostsec*40);
    %Which days of data to read?
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
    IDENTIF=[YEAR,'.',JDAY,'.',fam,'.loff',num2str(loopoffmax),'.ccmin',num2str(xcmaxAVEnmin),'.nponpa',int2str(npo),int2str(npa)]
    direc=[YEAR,'/',MO,'/'];
    prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];

    %Read the data; find glitches (many consecutive zeros) and flag those windows in STAnzeros.
    %Get timsSTA from the permanent stations (last one over-writes):
    STAopt=zeros(nsta,tracelen);
    %STAort=STAopt;
    STAnzeros=zeros(nsta,nwin);
    for ista=1:nsta
        found=0;
        [LIA,idx]=ismember(stas(ista,:),PERMSTA,'rows');
        if LIA
            found=found+LIA;
            if strcmp(PERMSTA(idx,1:3),'PGC')
                fact=1.6e-4;
            elseif strcmp(PERMSTA(idx,1:3),'LZB')
                fact=4.e-3;
            end
            [opt,ort,nzeros,timsSTA]=readperms(prename,PERMSTA,PERMROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
        end
        [LIA,idx]=ismember(stas(ista,:),POLSTA,'rows');
        if LIA
            found=found+LIA; %better be 1
            if year==2003 && jday<213
                fact=20.0e-3;
            else
                fact=4.0e-3; 
            end
            [opt,ort,nzeros]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
            timsSTA=0:1/sps:86400-1/sps;
        end
        found=found
        %factr1(ista)=prctile(abs(opt),90); %Not normalized
        %factr2(ista)=factr1(ista)/factr1(1); %This is what is used; keeps 1st station unchanged but scales the others
        STAopt(ista,:)=opt/scaleseisms(ista); 
        %STAort(ista,:)=ort;
        STAnzeros(ista,:)=nzeros;
    end
    %Now for broader band (bb)
    lobb=0.5;
    hibb=8;
    STAoptbb=zeros(nsta,tracelen);
    for ista=1:nsta
        [LIA,idx]=ismember(stas(ista,:),PERMSTA,'rows');
        if LIA
            if strcmp(PERMSTA(idx,1:3),'PGC')
                fact=1.6e-4;
            elseif strcmp(PERMSTA(idx,1:3),'LZB')
                fact=4.e-3;
            end
            [opt,~,~,~]=readperms(prename,PERMSTA,PERMROTS,idx,sps,lobb,hibb,npo,npa,fact,nwin,winlen,winoff,igstart);
        end
        [LIA,idx]=ismember(stas(ista,:),POLSTA,'rows');
        if LIA
            if year==2003 && jday<213
                fact=20.0e-3;
            else
                fact=4.0e-3; 
            end
            [opt,~,~]=readpols(prename,POLSTA,POLROTS,idx,sps,lobb,hibb,npo,npa,fact,nwin,winlen,winoff,igstart);
        end
        STAoptbb(ista,:)=opt/scaleseisms(ista); 
    end
    %Now for optimal P
    STAoptP=zeros(nstaP,tracelen);
    for ista=1:nstaP
        [LIA,idx]=ismember(stasP(ista,:),PERMSTA,'rows');
        if LIA
            found=found+LIA;
            if strcmp(PERMSTA(idx,1:3),'PGC')
                fact=1.6e-4;
            elseif strcmp(PERMSTA(idx,1:3),'LZB')
                fact=4.e-3;
            end
            [optP]=readpermsP(prename,PERMSTA,PERMPROTS,idx,sps,lo,hi+10,npo,npa,fact,nwin,winlen,winoff,igstart);
        end
        [LIA,idx]=ismember(stasP(ista,:),POLSTA,'rows');
        if LIA
            found=found+LIA; %better be 1
            if year==2003 && jday<213
                fact=20.0e-3;
            else
                fact=4.0e-3; 
            end
            [optP]=readpolsP(prename,POLSTA,POLPROTS,idx,sps,lo,hi+10,npo,npa,fact,nwin,winlen,winoff,igstart);
        end
        STAoptP(ista,:)=optP/scalePseisms(ista); 
    end
    
    figure
    timP=48858;
    plot(STAoptP(timP*sps+1:timP*sps+wintot))

end
