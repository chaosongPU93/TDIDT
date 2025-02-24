% From wdwSYNTH; this one stretches and scales templates to match "strongest arrival".
% As of 12/20/17:  
% Writes to 2018 directories
% Re-introduces sub-sample time offsets, based on interpolated dot products
%    (NOT interpolated seisms).  
% writes time stretches as factors, not indices that must be remembered
%    when analyzing data.
% Includes an offset buffer of 0.5*concntr in the prior 1.25s data
format short e
clear all
close all
stretchy=1
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');
% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo? 
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs
scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite; 

fam='002';
if isequal(fam,'002')  %DON'T FORGET additional family-specific delcarations around line 156 
        datnames=['3STAS/STAS.UN.100sps.mdist1_else0.25_N10e3.26';
                  '3STAS/STAS.UN.100sps.mdist1_else0.25_N10e4.26';
                  '3STAS/STAS.UN.100sps.mdist1_else0.25_N10e5.26';
                  '3STAS/STAS.UN.100sps.mdist1_else0.25_N10e6.26';
                  '3STAS/STAS.UN.100sps.mdist1_else0.25_N10e7.26'];
        datnames=['3STAS/STAS.UN.100sps.ell_3-2.else0nsat0.1'; %This one done later; sources seem different (not part of the same sequance as all the rest)
                  '3STAS/STAS.UN.100sps.ell_3-2.else0nsat0.4';
                  '3STAS/STAS.UN.100sps.ell_3-2.else0nsat1  ';
                  '3STAS/STAS.UN.100sps.ell_3-2.else0nsat2  ';
                  '3STAS/STAS.UN.100sps.ell_3-2.else0nsat4  ';
                  '3STAS/STAS.UN.100sps.ell_3-2.else0nsat10 ';
                  '3STAS/STAS.UN.100sps.ell_3-2.else0nsat100'];
        datnames=['3STAS/STAS.UN.100sps.ell_0.5-0.5.else0.25nsat0.1'; 
                  '3STAS/STAS.UN.100sps.ell_0.5-0.5.else0.25nsat0.4';
                  '3STAS/STAS.UN.100sps.ell_0.5-0.5.else0.25nsat1  ';
                  '3STAS/STAS.UN.100sps.ell_0.5-0.5.else0.25nsat2  ';
                  '3STAS/STAS.UN.100sps.ell_0.5-0.5.else0.25nsat4  ';
                  '3STAS/STAS.UN.100sps.ell_0.5-0.5.else0.25nsat10 ';
                  '3STAS/STAS.UN.100sps.ell_0.5-0.5.else0.25nsat100'];
        datnames=['3STAS/STAS.UN.100sps.rec_-52.2513.5.else0nsat0.1'; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0nsat0.4'; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0nsat1  '; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0nsat2  '; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0nsat4  '; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0nsat10 '; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0nsat100']; 
        datnames=['3STAS/STAS.UN.100sps.rec_-52.2513.5.else0.25nsat0.1'; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0.25nsat0.4'; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0.25nsat1  '; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0.25nsat2  '; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0.25nsat4  '; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0.25nsat10 '; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0.25nsat100'; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0.25nsat500']; 
        datnames=['3STAS/STAS.UN.100sps.rec_-52.2513.5.else0.5nsat0.1'; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0.5nsat0.4'; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0.5nsat1  '; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0.5nsat2  '; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0.5nsat4  '; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0.5nsat10 '; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0.5nsat100'; 
                  '3STAS/STAS.UN.100sps.rec_-52.2513.5.else0.5nsat500']; 
        datnames=['3STAS23/STAS.UN.100sps.ell_3-1.5.else0nsat0.1'; 
                  '3STAS23/STAS.UN.100sps.ell_3-1.5.else0nsat0.4';
                  '3STAS23/STAS.UN.100sps.ell_3-1.5.else0nsat1  ';
                  '3STAS23/STAS.UN.100sps.ell_3-1.5.else0nsat2  ';
                  '3STAS23/STAS.UN.100sps.ell_3-1.5.else0nsat4  ';
                  '3STAS23/STAS.UN.100sps.ell_3-1.5.else0nsat10 ';
                  '3STAS23/STAS.UN.100sps.ell_3-1.5.else0nsat100'];
        datnames=['3STAS23/STAS.UN.100sps.ell_2.75-1.else0nsat0.1'; 
                  '3STAS23/STAS.UN.100sps.ell_2.75-1.else0nsat0.4';
                  '3STAS23/STAS.UN.100sps.ell_2.75-1.else0nsat1  ';
                  '3STAS23/STAS.UN.100sps.ell_2.75-1.else0nsat2  ';
                  '3STAS23/STAS.UN.100sps.ell_2.75-1.else0nsat4  ';
                  '3STAS23/STAS.UN.100sps.ell_2.75-1.else0nsat10 '];
        datnames=['3STAS23/STAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat0.1'; 
                  '3STAS23/STAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat0.4';
                  '3STAS23/STAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat1  ';
                  '3STAS23/STAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat2  ';
                  '3STAS23/STAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat4  ';
                  '3STAS23/STAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat10 ';
                  '3STAS23/STAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat100'];
        datnames=['3STAS23/STAS.UN.100sps.ell_3-2.diam1.else0nsat0.1'; 
                  '3STAS23/STAS.UN.100sps.ell_3-2.diam1.else0nsat0.4';
                  '3STAS23/STAS.UN.100sps.ell_3-2.diam1.else0nsat1  ';
                  '3STAS23/STAS.UN.100sps.ell_3-2.diam1.else0nsat2  ';
                  '3STAS23/STAS.UN.100sps.ell_3-2.diam1.else0nsat4  ';
                  '3STAS23/STAS.UN.100sps.ell_3-2.diam1.else0nsat10 '];
        datnames=['3STAS23/STAS.UN.100sps.ell_3-1.5.diam1.else0nsat0.1'; 
                  '3STAS23/STAS.UN.100sps.ell_3-1.5.diam1.else0nsat0.4';
                  '3STAS23/STAS.UN.100sps.ell_3-1.5.diam1.else0nsat1  ';
                  '3STAS23/STAS.UN.100sps.ell_3-1.5.diam1.else0nsat2  ';
                  '3STAS23/STAS.UN.100sps.ell_3-1.5.diam1.else0nsat4  ';
                  '3STAS23/STAS.UN.100sps.ell_3-1.5.diam1.else0nsat10 '];   datsps=100;
        datnames=['3STAS23/STAS.UN.100sps.ell_2.75-1.25.diam0.3.else0nsat0.1'; 
                  '3STAS23/STAS.UN.100sps.ell_2.75-1.25.diam0.3.else0nsat0.4';
                  '3STAS23/STAS.UN.100sps.ell_2.75-1.25.diam0.3.else0nsat1  ';
                  '3STAS23/STAS.UN.100sps.ell_2.75-1.25.diam0.3.else0nsat2  ';
                  '3STAS23/STAS.UN.100sps.ell_2.75-1.25.diam0.3.else0nsat4  ';
                  '3STAS23/STAS.UN.100sps.ell_2.75-1.25.diam0.3.else0nsat10 ';
                  '3STAS23/STAS.UN.100sps.ell_2.75-1.25.diam0.3.else0nsat100'];
        datnames=['3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0nsat0.1'; 
                  '3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0nsat0.4';
                  '3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0nsat1  ';
                  '3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0nsat2  ';
                  '3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0nsat4  ';
                  '3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0nsat10 ';
                  '3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0nsat100'];
        datnames=['3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0.5nsat0.1'; 
                  '3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0.5nsat0.4';
                  '3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0.5nsat1  ';
                  '3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0.5nsat2  ';
                  '3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0.5nsat4  ';
                  '3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0.5nsat10 ';
                  '3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0.5nsat40 '];
        datnames=['3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0.25nsat0.1'; 
                  '3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0.25nsat0.4';
                  '3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0.25nsat1  ';
                  '3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0.25nsat2  ';
                  '3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0.25nsat4  ';
                  '3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0.25nsat10 ';
                  '3STAS23/STAS.UN.100sps.ell_1.5-0.75.diam0.3.else0.25nsat40 '];
        datnames=['3STAS23/STAS.UN.3621s.ell_1.5-0.75.diam0.3.else0.25nsat0.1'; 
                  '3STAS23/STAS.UN.3621s.ell_1.5-0.75.diam0.3.else0.25nsat0.4';
                  '3STAS23/STAS.UN.3621s.ell_1.5-0.75.diam0.3.else0.25nsat1  ';
                  '3STAS23/STAS.UN.3621s.ell_1.5-0.75.diam0.3.else0.25nsat2  ';
                  '3STAS23/STAS.UN.3621s.ell_1.5-0.75.diam0.3.else0.25nsat4  ';
                  '3STAS23/STAS.UN.3621s.ell_1.5-0.75.diam0.3.else0.25nsat10 ';
                  '3STAS23/STAS.UN.3621s.ell_1.5-0.75.diam0.3.else0.25nsat40 '];
        datnames=['3STAS23/STAS.UN.3621s.ell_2.75-1.25.diam0.3.else0nsat0.1'; 
                  '3STAS23/STAS.UN.3621s.ell_2.75-1.25.diam0.3.else0nsat0.4';
                  '3STAS23/STAS.UN.3621s.ell_2.75-1.25.diam0.3.else0nsat1  ';
                  '3STAS23/STAS.UN.3621s.ell_2.75-1.25.diam0.3.else0nsat2  ';
                  '3STAS23/STAS.UN.3621s.ell_2.75-1.25.diam0.3.else0nsat4  ';
                  '3STAS23/STAS.UN.3621s.ell_2.75-1.25.diam0.3.else0nsat10 ';
                  '3STAS23/STAS.UN.3621s.ell_2.75-1.25.diam0.3.else0nsat40 '];
        datnames=['3STAS23/STAS.UN.1821s.ell_1.5-0.75.diam0.3.else0.25nsat0.1'; 
                  '3STAS23/STAS.UN.1821s.ell_1.5-0.75.diam0.3.else0.25nsat0.4';
                  '3STAS23/STAS.UN.1821s.ell_1.5-0.75.diam0.3.else0.25nsat1  ';
                  '3STAS23/STAS.UN.1821s.ell_1.5-0.75.diam0.3.else0.25nsat2  ';
                  '3STAS23/STAS.UN.1821s.ell_1.5-0.75.diam0.3.else0.25nsat4  ';
                  '3STAS23/STAS.UN.1821s.ell_1.5-0.75.diam0.3.else0.25nsat10 ';
                  '3STAS23/STAS.UN.1821s.ell_1.5-0.75.diam0.3.else0.25nsat40 '];
    stas=['PGC '
          'SSIB'
          'SILB'];  
    %tempoffs=[911 998 932]; %these are the OLD zero crossings for 002: PGC,SSIB,SILB.  87 or 86 samples for B's templates, it seems. 
    tempoffs=[101 101 101]; %these are the zero crossings for 002: PGC,SSIB,SILB.  
    tempzeros=[93 109;
               92 109;
               92 108]; %These, exclusive of the endpoints, are the samples that don't get zeroed to make a single dipole.
elseif isequal(fam,'068')
    stas=['TWKB'
          'LZB '
          'MGCB'];
    tempoffs=[902 894 901]; %these are the (OLD!!!) zero crossings for 068: TWKB,LZB,MGCB.
end
nsta=size(stas,1);
sps=40;
spsm1=1/sps;

%Basics of the cross-correlation:  Window length, number of windows, filter parameters, etc.
winlensec=4;
winoffsec=1;
% winlensec=128;
% winoffsec=64;
winlen=winlensec*sps;
winoff=winoffsec*sps;
%dat=load(datnames(1,:));

hi=6.5;
lo=1.25;
npo=2;
npa=2;
mshift=19; %maximum shift for the x-correlations. 19 for 002 Stanford
loopoffmax=1.5; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.
xcmaxAVEnmin=0.44; %0.11; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the templates for the stations. Won't be always required.  PGCopt_002_1.25-6.5Hz_2pass_100sps.le14shift.resamp
for ista=1:nsta
    if stas(ista,4) == ' '
        temptemps(ista,:)=load([stas(ista,1:3),'opt_',fam,'_',num2str(lo),'-',num2str(hi), ...
                         'Hz_',int2str(npa),'pass_',int2str(sps),'sps.input'],'w');
    else
        temptemps(ista,:)=load([stas(ista,1:4),'opt_',fam,'_',num2str(lo),'-',num2str(hi), ...
                         'Hz_',int2str(npa),'pass_',int2str(sps),'sps.input'],'w');
    end
end
tempbef=59;
tempaft=60;
templen=tempbef+tempaft+1;
for ista=1:nsta
    STAtemps(ista,:)=temptemps(ista,tempoffs(ista)-tempbef:tempoffs(ista)+tempaft);
    %snips(templen*(ista-1)+1:ista*templen)=STAtemps(ista,:);
end
% scalefact scales templates; scaleseisms scales seisms.  Strategy changes with family.
if isequal(fam,'002')
    tempoffs=tempoffs-1; %Center of strongest window is 1 or 2 samples ahead of zero crossing (002); make it 1.
    whichtoplot=2;
    scaleseisms=[1.0 0.76 0.95];  %I think this is now important only for plotting purposes.
elseif isequal(fam,'068')
    whichtoplot=1;
    scaleseisms=[1.0 1.0 1.0];
end
minses=-min(STAtemps,[],2); %STAtemps is (currently) 3 by 120
maxses= max(STAtemps,[],2);
plustominus=maxses./minses;
scalefact=minses*max(plustominus); %This is used to scale templates, just for plotting purposes
for ista=1:nsta
    STAtemps(ista,:)=STAtemps(ista,:)/scalefact(ista); %This plots the templates with the largest positive value (of any) at +1
end
figure
plot(STAtemps(1,:),'r')
hold on
plot(STAtemps(2,:),'b')
plot(STAtemps(3,:),'k')
drawnow
% A different strategy for stretched templates?
minstretch=20;
maxstretch=80;
stretchlen=60; %This is duration of window, in samples; not max-min
STAstr=zeros(nsta,maxstretch-minstretch+1,stretchlen);
ttemp=temptemps;
minses=-min(ttemp,[],2); 
maxses= max(ttemp,[],2);
plustominus=maxses./minses;
scalefact=minses*max(plustominus); %This is used to scale templates, just for plotting purposes
for ista=1:nsta
    ttemp(ista,:)=ttemp(ista,:)/scalefact(ista); %This plots the templates with the largest positive value (of any) at +1
end
figure
hold on
for ista=1:nsta
    ttemp(ista,1:tempzeros(ista,1))=0;
    ttemp(ista,tempzeros(ista,2):end)=0;
    temp=ttemp(ista,tempoffs(ista)-stretchlen+1:tempoffs(ista)+stretchlen); %this has some slop at the ends 
    plot(temp)
    STAstr(ista,:,:)=stretchtemps(temp,minstretch,maxstretch,stretchlen,sps);
end
figure
toplot=squeeze(STAstr(:,1,:));
plot(toplot(1,:),'r')
hold on
plot(toplot(2,:),'b')
plot(toplot(3,:),'k')
figure
toplot=squeeze(STAstr(:,end,:));
plot(toplot(1,:),'r')
hold on
plot(toplot(2,:),'b')
plot(toplot(3,:),'k')
figure

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cycle over each day:
for nd=1:size(datnames,1) %a little kludgy; takes input file names of different lengths and makes them readable and writable at similar lengths.
    for j=1:size(datnames,2)
        if not(strcmp(datnames(nd,j),' '))
            empt=j+1;
        end
    end
    tmpname=datnames(nd,:);
    for ilett=empt:size(datnames,2)
        tmpname(ilett)='.';
    end
    %IDENTIF=[tmpname(7:end),'.loff',num2str(loopoffmax),'.ccmin',num2str(xcmaxAVEnmin),'.nponpa',int2str(npo),int2str(npa)]
    IDENTIF=[tmpname(9:end),'.loff',num2str(loopoffmax),'.ccmin',num2str(xcmaxAVEnmin),'.nponpa',int2str(npo),int2str(npa)]
    tmp=load(datnames(nd,1:empt-1))';
    tmport=load([datnames(nd,1:empt-1),'.orth'])';
    for ista=1:nsta
        [STAop(ista,:)]=bandpass(tmp(ista,:),datsps,lo,hi,npo,npa,'butter');
        STAopt(ista,:)=resample(STAop(ista,:),2,5);
        STAhil(ista,:)=abs(hilbert(STAopt(ista,:)));
        [STAor(ista,:)]=bandpass(tmport(ista,:),datsps,lo,hi,npo,npa,'butter');
        STAort(ista,:)=resample(STAor(ista,:),2,5);
    end
    medhil=median(STAhil); %Median of the 3 hilbert transforms for each sample.    
    tracelen=size(STAopt,2); 
    winbig=2*(tracelen/2-sps); %ignore 1 second at each end of data
    timbig=winbig/(2*sps); %half that time, in seconds
    igstart=floor(tracelen/2-winbig/2)+1; %start counting seis data from here
    nwin=floor((winbig-winlen)/winoff);
    timsSTA=spsm1*(1:tracelen);
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Autocorrelation of stations.  Those that end in "sq" are the running
    %   cumulative sum, to be used later by differncing the window edpoints.
    %   (Used to be PGCauto, PGC2, SSIBauto, SSIB2, etc.)
    %   Station to itself is in a 3 x tracelen array
    %   Cross-station measurements are in their own linear array
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    STAauto=STAopt.*STAopt;
    STAsq=cumsum(STAauto,2);
    STAortauto=STAort.*STAort;
    STAortsq=cumsum(STAortauto,2);
    clear STAortauto
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Cross-correlation between stations, with small offsets up to +/- mshift.
    %  First index is pointwise multiplication of traces; second is shifting offset.
    %  lenx is shorter than tracelen by mshift at each end (see notebook sketch)
    %  For STA12 and PGSI, SSI and SIL are shifted relative to PGC, by 1 each time through loop.
    %  For SISS, SSI is shifted relative to SILB.
    %  cumsumSTA12 etc. are the running cumulative sum of the x-correlation.
    %  PGSSx becomes STA12x, PGSI -> STA13, SISS -> STA32
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lenx=tracelen-2*mshift;
    STA12x=zeros(lenx, 2*mshift+1);
    STA13x=zeros(lenx, 2*mshift+1);
    STA32x=zeros(lenx, 2*mshift+1);
    for n=-mshift:mshift;
        STA12x(:,n+mshift+1)=STAopt(1,1+mshift:tracelen-mshift).* ...
            STAopt(2,1+mshift-n:tracelen-mshift-n);
        STA13x(:,n+mshift+1)=STAopt(1,1+mshift:tracelen-mshift).* ...
            STAopt(3,1+mshift-n:tracelen-mshift-n);
        STA32x(:,n+mshift+1)=STAopt(3,1+mshift:tracelen-mshift).* ...
            STAopt(2,1+mshift-n:tracelen-mshift-n);
    end
    cumsumSTA12=cumsum(STA12x);  %prev cumsumPGSS
    cumsumSTA13=cumsum(STA13x);  %prev cumsumPGSI
    cumsumSTA32=cumsum(STA32x);  %prev cumsumSISS
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  "winbig" is now the whole day, minus 2 sec at each end (apparently).
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
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timswin=zeros(nwin,1);
%     sumsPGSS=zeros(nwin,2*mshift+1);
%     sumsPGSI=zeros(nwin,2*mshift+1);
%     sumsSISS=zeros(nwin,2*mshift+1);
%     sumsPGC2=zeros(nwin,2*mshift+1);
%     sumsSSIB2=zeros(nwin,2*mshift+1);
%     sumsSILB2=zeros(nwin,2*mshift+1);
%     sumsSILB2b=zeros(nwin,2*mshift+1);
    sumsSTA12=zeros(nwin,2*mshift+1);
    sumsSTA13=zeros(nwin,2*mshift+1);
    sumsSTA32=zeros(nwin,2*mshift+1);
    sumsSTA1sq=zeros(nwin,2*mshift+1);
    sumsSTA2sq=zeros(nwin,2*mshift+1);
    sumsSTA3sq=zeros(nwin,2*mshift+1);
    sumsSTA3Bsq=zeros(nwin,2*mshift+1);
%The following commented out 7/12/18 to fix the mshift business.
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     %  sumsPGSS is shorter than sumsPGC2 by 2*mshift.  This is why sumsPGC2 etc
% %     %  is shifted by +mshift.  cumsumPGSS(1,:)=cumsum(PGSSx)(1,:) starts mshift
% %     %  to the right of the first data sample.  igstart is how many to the right
% %     %  of that.
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     for n=1:nwin;
% %         istart=igstart+(n-1)*winoff;
% %         iend=istart+winlen;
% %         timswin(n)=timsSTA(istart+winlen/2); 
% %         sumsSTA12(n,:)=cumsumSTA12(iend,:)-cumsumSTA12(istart-1,:); 
% %         sumsSTA13(n,:)=cumsumSTA13(iend,:)-cumsumSTA13(istart-1,:);
% %         sumsSTA32(n,:)=cumsumSTA32(iend,:)-cumsumSTA32(istart-1,:);
% %         sumsSTA1sq(n,:)=STAsq(1,iend+mshift)-STAsq(1,istart+mshift-1);  %PGC2 is cumsummed. Yes, +mshift.
% %         sumsSTA3Bsq(n,:)=STAsq(3,iend+mshift)-STAsq(3,istart+mshift-1); %Similar, for the SILB-SSIB connection.
% %         for m=-mshift:mshift;
% %             sumsSTA2sq(n,m+mshift+1)=STAsq(2,iend+mshift-m)-STAsq(2,istart+mshift-1-m); %+m??? (yes).
% %             sumsSTA3sq(n,m+mshift+1)=STAsq(3,iend+mshift-m)-STAsq(3,istart+mshift-1-m);
% %         end
% %     end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  sumsPGSS is shorter than sumsPGC2 by 2*mshift.  This is why sumsPGC2 etc
    %  is shifted by +mshift.  cumsumPGSS(1,:)=cumsum(PGSSx)(1,:) starts mshift
    %  to the right of the first data sample.  igstart is how many to the right
    %  of that.
    %  07/06/2018:  I'm pretty sure I want to subtract mshift from every iend or
    %  istart index to the right of the equal signs in the following FOR loop.
    %  The commented version above looks like it has the proper relative shifts but 
    %  the absolute time isn't registered to the global igstart properly.
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for n=1:nwin;
        istart=igstart+(n-1)*winoff;
        iend=istart+winlen;
        timswin(n)=timsSTA(istart+winlen/2); 
        sumsSTA12(n,:)=cumsumSTA12(iend-mshift,:)-cumsumSTA12(istart-mshift-1,:);  %Yes, -mshift (7/06/18)
        sumsSTA13(n,:)=cumsumSTA13(iend-mshift,:)-cumsumSTA13(istart-mshift-1,:);
        sumsSTA32(n,:)=cumsumSTA32(iend-mshift,:)-cumsumSTA32(istart-mshift-1,:);
        sumsSTA1sq(n,:)=STAsq(1,iend)-STAsq(1,istart-1);  %PGC2 is cumsummed. Yes, +mshift.  No, no mshift (7/06/18)
        sumsSTA3Bsq(n,:)=STAsq(3,iend)-STAsq(3,istart-1); %Similar, for the SILB-SSIB connection.
        for m=-mshift:mshift;
            sumsSTA2sq(n,m+mshift+1)=STAsq(2,iend-m)-STAsq(2,istart-1-m); %+m??? (yes).
            sumsSTA3sq(n,m+mshift+1)=STAsq(3,iend-m)-STAsq(3,istart-1-m);
        end
    end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Denominator for the normalization.  A 2D array, nwin by 2*mshift+1.
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %An attempt to bypass glitches in data.  Min value of good data typically ~10^{-2}
%     glitches=1.e-7;
%     sumsSTA1sq=max(sumsSTA1sq,glitches);
%     sumsSTA2sq=max(sumsSTA2sq,glitches);
%     sumsSTA3sq=max(sumsSTA3sq,glitches);
    %
    denomSTA12n=realsqrt(sumsSTA1sq.*sumsSTA2sq);
    denomSTA13n=realsqrt(sumsSTA1sq.*sumsSTA3sq);
    denomSTA32n=realsqrt(sumsSTA3Bsq.*sumsSTA2sq);
    %
    sumsSTA12n=sumsSTA12./denomSTA12n;
    sumsSTA13n=sumsSTA13./denomSTA13n;
    sumsSTA32n=sumsSTA32./denomSTA32n;
    [xcmaxSTA12n,imaxSTA12]=max(sumsSTA12n,[],2); %Integer-offset max cross-correlation
    [xcmaxSTA13n,imaxSTA13]=max(sumsSTA13n,[],2);
    [xcmaxSTA32n,imaxSTA32]=max(sumsSTA32n,[],2);
    %Parabolic fit:
    [xmaxSTA12n,ymaxSTA12n,aSTA12]=parabol(nwin,mshift,sumsSTA12n,imaxSTA12); %Interpolated max cross-correlation
    [xmaxSTA13n,ymaxSTA13n,aSTA13]=parabol(nwin,mshift,sumsSTA13n,imaxSTA13);
    [xmaxSTA32n,ymaxSTA32n,aSTA32]=parabol(nwin,mshift,sumsSTA32n,imaxSTA32);
    
    %h=figure('Position',[0.1*wid 1 2.5*wid hite]); %center

    ix=sub2ind(size(denomSTA12n),(1:nwin)',imaxSTA12); %Find the linear index of the largest denominator
    ampSTA12=sqrt(denomSTA12n(ix)); %This makes amplitude linear rather than quadratic with counts.
    ampSTA1sq=sumsSTA1sq(ix); %by construction PGC2 is the same for all shifts  % sumsPGC2 becomes sumsSTA1sq
    ampSTA2sq=sumsSTA2sq(ix); % sumsSSIB2 becomes sumsSTA2sq
    ix=sub2ind(size(denomSTA13n),(1:nwin)',imaxSTA13);
    ampSTA13=sqrt(denomSTA13n(ix));
    ampSTA3sq=sumsSTA3sq(ix);
    ix=sub2ind(size(denomSTA32n),(1:nwin)',imaxSTA32);
    ampSTA32=sqrt(denomSTA32n(ix));
    AmpComp(1:4)=0;
    %AmpComp seems to be amplitude squared in 4s window minus amp squared in prior window,
    %divided by sum of amp squared in the two windows.  And why?
    AmpComp(5:nwin)=((ampSTA1sq(5:nwin)+ampSTA2sq(5:nwin)+ampSTA3sq(5:nwin))- ...
                    (ampSTA1sq(1:nwin-4)+ampSTA2sq(1:nwin-4)+ampSTA3sq(1:nwin-4)))./ ...
                    ((ampSTA1sq(5:nwin)+ampSTA2sq(5:nwin)+ampSTA3sq(5:nwin))+ ...
                    (ampSTA1sq(1:nwin-4)+ampSTA2sq(1:nwin-4)+ampSTA3sq(1:nwin-4))) ;
    %Center them
    imaxSTA12cent=imaxSTA12-mshift-1;  % "cent" is "centered"; imaxSTA12 is original 1:2*mshift+1
    imaxSTA13cent=imaxSTA13-mshift-1;
    imaxSTA32cent=imaxSTA32-mshift-1;
    iloopoff=imaxSTA13cent-imaxSTA12cent+imaxSTA32cent; %How well does the integer loop close?
    xmaxSTA12n=xmaxSTA12n-mshift-1;
    xmaxSTA13n=xmaxSTA13n-mshift-1;
    xmaxSTA32n=xmaxSTA32n-mshift-1;
    loopoff=xmaxSTA13n-xmaxSTA12n+xmaxSTA32n; %How well does the interpolated loop close?
    xcmaxAVEn=(xcmaxSTA12n+xcmaxSTA13n+xcmaxSTA32n)/3;
    % xcnshifts=cputime-t
    % t=cputime;
    for n=1:nwin;  %Added 12/29/21.  Maybe there's a more efficient way of doing this.
        istart=igstart+(n-1)*winoff;
        iend=istart+winlen; %(-1, says Chao)
        ampopt(n)=STAsq(1,iend)-STAsq(1,istart-1)+...
               STAsq(2,iend-imaxSTA12(n))-STAsq(2,istart-imaxSTA12(n)-1)+...
               STAsq(3,iend-imaxSTA13(n))-STAsq(3,istart-imaxSTA13(n)-1);
        amport(n)=STAortsq(1,iend)-STAortsq(1,istart-1)+...
               STAortsq(2,iend-imaxSTA12(n))-STAortsq(2,istart-imaxSTA12(n)-1)+...
               STAortsq(3,iend-imaxSTA13(n))-STAortsq(3,istart-imaxSTA13(n)-1);
    end
    ampmax=max([ampSTA12; ampSTA13; ampSTA32]);
    medxcmaxAVEn=median(xcmaxAVEn)
    xmaxSTA12ntmp=xmaxSTA12n;
    xmaxSTA13ntmp=xmaxSTA13n;
    xmaxSTA32ntmp=xmaxSTA32n;

    iup=4; 
    nin=0;
    zerosallowed=20*winlen/160;
    concentration=0.55; %in seconds; how concentrated is the coherent energy within the window?
    cncntr=concentration*sps;
    offset=round(0.5*cncntr);
    for n=1:nwin
        if xcmaxAVEn(n)<xcmaxAVEnmin || abs(xmaxSTA13n(n)-xmaxSTA12n(n)+xmaxSTA32n(n))>loopoffmax ...
                || isequal(abs(imaxSTA12cent(n)),mshift) || isequal(abs(imaxSTA13cent(n)),mshift) ...
                || isequal(abs(imaxSTA32cent(n)),mshift) || isequal(timswin(n),timsSTA(igstart+winlen/2))     
            xmaxSTA12ntmp(n)=mshift+1; xmaxSTA13ntmp(n)=mshift+1; xmaxSTA32ntmp(n)=mshift+1; %dummy them, if these criteria are met
        else
            % %%%%%%%%%%%%%%%%%%%%%
            interpSTA12n=interp(sumsSTA12n(n,:),iup,3);
            interpSTA13n=interp(sumsSTA13n(n,:),iup,3);
            interpSTA32n=interp(sumsSTA32n(n,:),iup,3);
            leninterp=length(interpSTA12n);
            [xcmaxinterpSTA12n,imaxinterpSTA12]=max(interpSTA12n(1:leninterp-(iup-1)));
            [xcmaxinterpSTA13n,imaxinterpSTA13]=max(interpSTA13n(1:leninterp-(iup-1)));
            [xcmaxinterpSTA32n,imaxinterpSTA32]=max(interpSTA32n(1:leninterp-(iup-1)));
            xcmaxconprev=-99999.;  %used to be 0; not good with glitches
            for iSTA12=max(1,imaxinterpSTA12-3*iup):min(imaxinterpSTA12+3*iup,iup*(2*mshift+1)-(iup-1)) %3 samples from peak; 
                                                                                     %intentionally wider than acceptable;
                                                                                     %iup-1 are extrapolated points
                for iSTA13=max(1,imaxinterpSTA13-3*iup):min(imaxinterpSTA13+3*iup,iup*(2*mshift+1)-(iup-1))
                    ibangon = (iup*mshift+1)-iSTA13+iSTA12;
                    if ibangon >= 1 && ibangon<=iup*(2*mshift+1)
                        xcmaxcon=interpSTA12n(iSTA12)+interpSTA13n(iSTA13)+interpSTA32n(ibangon);
                        if xcmaxcon > xcmaxconprev
                            xcmaxconprev=xcmaxcon;
                            iSTA12bang=iSTA12;
                            iSTA13bang=iSTA13;
                        end
                    end
                end
            end
            % %%%%%%%%%%%%%%%%%%%%%
%             xcmaxconprev=-99999.;  %used to be 0; not good with glitches
%             imaxSTA12n=imaxSTA12(n); %This "n" for nth window; other "n's" for "normalized".  Unfortunately.
%             imaxSTA13n=imaxSTA13(n);
%             imaxSTA32n=imaxSTA32(n);
%             sumsSTA12nn=sumsSTA12n(n,:);
%             sumsSTA13nn=sumsSTA13n(n,:);
%             sumsSTA32nn=sumsSTA32n(n,:);
%             for iSTA12 =     max(1,imaxSTA12n-floor(loopoffmax+1)):min(imaxSTA12n+floor(loopoffmax+1),2*mshift+1)
%                 for iSTA13 = max(1,imaxSTA13n-floor(loopoffmax+1)):min(imaxSTA13n+floor(loopoffmax+1),2*mshift+1)
%                     ibangon = (mshift+1)-iSTA13+iSTA12;
%                     if ibangon >= 1 && ibangon <= 2*mshift+1
%                         xcmaxcon=sumsSTA12nn(iSTA12)+sumsSTA13nn(iSTA13)+sumsSTA32nn(ibangon);
%                         if xcmaxcon > xcmaxconprev
%                             xcmaxconprev=xcmaxcon;
%                             iSTA12bang=iSTA12;
%                             iSTA13bang=iSTA13;
%                         end
%                     end
%                 end
%             end
            % %%%%%%%%%%%%%%%%%%%%%
            iSTA32bang=(iup*mshift+1)-iSTA13bang+iSTA12bang;
            if abs(iSTA12bang-imaxinterpSTA12) <= loopoffmax*iup && ...
               abs(iSTA13bang-imaxinterpSTA13) <= loopoffmax*iup && ...
               abs(iSTA32bang-imaxinterpSTA32) <= loopoffmax*iup && ...
               interpSTA12n(iSTA12bang)+interpSTA13n(iSTA13bang)+interpSTA32n(iSTA32bang) >= 3*xcmaxAVEnmin
                xmaxSTA12ntmp(n)=(iSTA12bang-(iup*mshift+1))/iup;
                xmaxSTA13ntmp(n)=(iSTA13bang-(iup*mshift+1))/iup;
                xmaxSTA32ntmp(n)=(iSTA32bang-(iup*mshift+1))/iup;
            % %%%%%%%%%%%%%%%%%%%%%
%             iSTA32bang=(mshift+1)-iSTA13bang+iSTA12bang;
%             if abs(iSTA12bang-imaxSTA12n) <= loopoffmax && ...  %not sure if these 3 lines are satisfied automatically ...
%                abs(iSTA13bang-imaxSTA13n) <= loopoffmax && ...
%                abs(iSTA32bang-imaxSTA32n) <= loopoffmax && ...
%                sumsSTA12n(n,iSTA12bang)+sumsSTA13n(n,iSTA13bang)+sumsSTA32n(n,iSTA32bang) >= 3*xcmaxAVEnmin
%                 xmaxSTA12ntmp(n)=iSTA12bang-(mshift+1); %without interpolation this is just centering.
%                 xmaxSTA13ntmp(n)=iSTA13bang-(mshift+1);
%                 xmaxSTA32ntmp(n)=iSTA32bang-(mshift+1);

                %for plotting traces
                imaxSTA12wr=round(xmaxSTA12ntmp(n)); %without interpolation this is not needed.
                imaxSTA13wr=round(xmaxSTA13ntmp(n));
    % 
                istart=igstart+(n-1)*winoff; %+mshift; %a better way might exist?  %ADDED mshift 10/20/12; DELETED IT 1/19/17.
                    %ADDED IT BACK 10/4/2017 to fix bug.  PGC is offset from igstart by mshift before first x-correlation.
                    %Not sure why mshift was added.  It changes STA12tr, STA1file etc. relative to the window that was used
                    %in the original x-correlation.  This will affect the stated time of max energy (through idiff).
                    %GOT RID of the mshift, yet again, 07/06/18, but only after subtracting mshift from all those istarts and
                    %iends starting around line 444.
                iend=istart+winlen-1;
                imid=round((istart+iend)/2);                
                %Check power spectrum for reasonableness
                [STA1xx fp] = pwelch(STAopt(1,istart:iend),[],[],[],40); %40 is sps
                STA1xx=STA1xx/max(STA1xx);
                [STA2xx fp] = pwelch(STAopt(2,istart-imaxSTA12wr:iend-imaxSTA12wr),[],[],[],40);
                STA2xx=STA2xx/max(STA2xx);
                [STA3xx fp] = pwelch(STAopt(3,istart-imaxSTA13wr:iend-imaxSTA13wr),[],[],[],40);
                STA3xx=STA3xx/max(STA3xx);
                flo=find(fp > lo,1)-1;
                fhi=find(fp > hi,1)+1; %extra 1 for good measure
                belowcut=median([STA1xx(2:flo); STA2xx(2:flo); STA3xx(2:flo)]);
                ppeaksSTA1=findpeaks(STA1xx(flo+1:fhi));
                if length(ppeaksSTA1)>=1
                    maxppeakSTA1=max(ppeaksSTA1);
                else
                    maxppeakSTA1=0.;
                end
                ppeaksSTA2=findpeaks(STA2xx(flo+1:fhi));
                if length(ppeaksSTA2)>=1
                    maxppeakSTA2=max(ppeaksSTA2);
                else
                    maxppeakSTA2=0.;
                end
                ppeaksSTA3=findpeaks(STA3xx(flo+1:fhi));
                if length(ppeaksSTA3)>=1
                    maxppeakSTA3=max(ppeaksSTA3);
                else
                    maxppeakSTA3=0.;
                end
                abovecut=median([maxppeakSTA1 maxppeakSTA2 maxppeakSTA3]);
                if abovecut > 0.9*belowcut %-1 %This checks for frequency range; make sure it's not too narrow?
                    STA12tr=STAopt(1,istart:iend).*STAopt(2,istart-imaxSTA12wr:iend-imaxSTA12wr);
                    STA13tr=STAopt(1,istart:iend).*STAopt(3,istart-imaxSTA13wr:iend-imaxSTA13wr);
                    STA32tr=STAopt(3,istart-imaxSTA13wr:iend-imaxSTA13wr).*STAopt(2,istart-imaxSTA12wr:iend-imaxSTA12wr);
                    cumsumtr=cumsum(STA12tr)+cumsum(STA13tr)+cumsum(STA32tr);
                    [cumsumtrdiff idiff]=max(cumsumtr(cncntr+1:winlen)-cumsumtr(1:winlen-cncntr));
                    %What is amp squared in strongest coherent 1/2 sec?
                    isdiff=istart+idiff-1; %Start of strong 0.5s
                    iediff=istart+idiff-1+cncntr;
                    dummy=STAopt(1,isdiff:iediff).*STAopt(1,isdiff:iediff)+ ...
                          STAopt(2,isdiff-imaxSTA12wr:iediff-imaxSTA12wr).*STAopt(2,isdiff-imaxSTA12wr:iediff-imaxSTA12wr)+ ...
                          STAopt(3,isdiff-imaxSTA13wr:iediff-imaxSTA13wr).*STAopt(3,isdiff-imaxSTA13wr:iediff-imaxSTA13wr);
                    dum2=cumsum(dummy);
                    Ampsq(nin+1)=dum2(end);
                    %Energy in prior 1.25 seconds, with offset
                    if isdiff-2.5*cncntr-offset < -5
                        Prev(nin+1)=999.;
                    else
                    dummy=STAopt(1,isdiff-2.5*cncntr-offset:isdiff-offset).*STAopt(1,isdiff-2.5*cncntr-offset:isdiff-offset)+ ...
                          STAopt(2,isdiff-2.5*cncntr-imaxSTA12wr-offset:isdiff-imaxSTA12wr-offset).*STAopt(2,isdiff-2.5*cncntr-imaxSTA12wr-offset:isdiff-imaxSTA12wr-offset)+ ...
                          STAopt(3,isdiff-2.5*cncntr-imaxSTA13wr-offset:isdiff-imaxSTA13wr-offset).*STAopt(3,isdiff-2.5*cncntr-imaxSTA13wr-offset:isdiff-imaxSTA13wr-offset);
                        dum2=cumsum(dummy);
                        Prev(nin+1)=dum2(end);
                    end
                    clear dummy
                    %CC in same window (test)
                    dummy(1,:)=STAopt(1,istart:iend);
                    dummy(2,:)=STAopt(2,istart-imaxSTA12wr:iend-imaxSTA12wr);
                    dummy(3,:)=STAopt(3,istart-imaxSTA13wr:iend-imaxSTA13wr);
                    denoms=dot(dummy,dummy,2);
                    cc(nin+1)=(dot(dummy(1,:),dummy(2,:))/sqrt(denoms(1)*denoms(2))+dot(dummy(2,:),dummy(3,:))/sqrt(denoms(2)*denoms(3))+ ...
                               dot(dummy(3,:),dummy(1,:))/sqrt(denoms(3)*denoms(1)))/3;
                    clear dummy
                    %CC in arrival
                    dummy(1,:)=STAopt(1,isdiff:iediff);
                    dummy(2,:)=STAopt(2,isdiff-imaxSTA12wr:iediff-imaxSTA12wr);
                    dummy(3,:)=STAopt(3,isdiff-imaxSTA13wr:iediff-imaxSTA13wr);
                    denoms=dot(dummy,dummy,2);
                    ccarr(nin+1)=(dot(dummy(1,:),dummy(2,:))/sqrt(denoms(1)*denoms(2))+dot(dummy(2,:),dummy(3,:))/sqrt(denoms(2)*denoms(3))+ ...
                                       dot(dummy(3,:),dummy(1,:))/sqrt(denoms(3)*denoms(1)))/3;
                    clear dummy
                    %CC in prior 1.25 seconds
                    if isdiff-2.5*cncntr-offset < -5
                        ccprior125(nin+1)=999;
                    else
                        dummy(1,:)=STAopt(1,isdiff-2.5*cncntr-offset:isdiff-1-offset);
                        dummy(2,:)=STAopt(2,isdiff-imaxSTA12wr-2.5*cncntr-offset:isdiff-imaxSTA12wr-1-offset);
                        dummy(3,:)=STAopt(3,isdiff-imaxSTA13wr-2.5*cncntr-offset:isdiff-imaxSTA13wr-1-offset);
                        denoms=dot(dummy,dummy,2);
                        ccprior125(nin+1)=(dot(dummy(1,:),dummy(2,:))/sqrt(denoms(1)*denoms(2))+dot(dummy(2,:),dummy(3,:))/sqrt(denoms(2)*denoms(3))+ ...
                                          dot(dummy(3,:),dummy(1,:))/sqrt(denoms(3)*denoms(1)))/3;
                    end
                    clear dummy
                    % %Energy in the 1.25s before that, with offset overlap (assuming 0.5sncntr)
                    if isdiff > 5*cncntr+mshift
                        dummy=STAopt(1,isdiff-5*cncntr:isdiff-2.5*cncntr).*STAopt(1,isdiff-5*cncntr:isdiff-2.5*cncntr)+ ...
                              STAopt(2,isdiff-5*cncntr-imaxSTA12wr:isdiff-2.5*cncntr-imaxSTA12wr).*STAopt(2,isdiff-5*cncntr-imaxSTA12wr:isdiff-2.5*cncntr-imaxSTA12wr)+ ...
                              STAopt(3,isdiff-5*cncntr-imaxSTA13wr:isdiff-2.5*cncntr-imaxSTA13wr).*STAopt(3,isdiff-5*cncntr-imaxSTA13wr:isdiff-2.5*cncntr-imaxSTA13wr);
                    else
                        dummy=STAopt(1,mshift:mshift+2.5*cncntr-1).*STAopt(1,mshift:mshift+2.5*cncntr-1)+ ...
                              STAopt(2,mshift-imaxSTA12wr:mshift-imaxSTA12wr+2.5*cncntr-1).*STAopt(2,mshift-imaxSTA12wr:mshift-imaxSTA12wr+2.5*cncntr-1)+ ...
                              STAopt(3,mshift-imaxSTA13wr:mshift-imaxSTA13wr+2.5*cncntr-1).*STAopt(3,mshift-imaxSTA13wr:mshift-imaxSTA13wr+2.5*cncntr-1);
                    end
                    dum2=cumsum(dummy);
                    Prev2(nin+1)=dum2(end);
                    clear dummy
                    %CC in prior 4 seconds
                    if isdiff > winlen+mshift+offset
                        dummy(1,:)=STAopt(1,isdiff-winlen-offset:isdiff-1-offset);
                        dummy(2,:)=STAopt(2,isdiff-imaxSTA12wr-winlen-offset:isdiff-imaxSTA12wr-1-offset);
                        dummy(3,:)=STAopt(3,isdiff-imaxSTA13wr-winlen-offset:isdiff-imaxSTA13wr-1-offset);
                    else
                        dummy(1,:)=STAopt(1,mshift:isdiff-1-offset);
                        dummy(2,:)=STAopt(2,mshift-imaxSTA12wr:isdiff-imaxSTA12wr-1-offset);
                        dummy(3,:)=STAopt(3,mshift-imaxSTA13wr:isdiff-imaxSTA13wr-1-offset);
                    end
                    denoms=dot(dummy,dummy,2);
                    ccprior(nin+1)=(dot(dummy(1,:),dummy(2,:))/sqrt(denoms(1)*denoms(2))+dot(dummy(2,:),dummy(3,:))/sqrt(denoms(2)*denoms(3))+ ...
                                    dot(dummy(3,:),dummy(1,:))/sqrt(denoms(3)*denoms(1)))/3;
                    clear dummy
                    %CC in following 4 seconds
                    if iediff < winbig-winlen-mshift
                        dummy(1,:)=STAopt(1,iediff+1:iediff+winlen);
                        dummy(2,:)=STAopt(2,iediff+1-imaxSTA12wr:iediff-imaxSTA12wr+winlen);
                        dummy(3,:)=STAopt(3,iediff+1-imaxSTA13wr:iediff-imaxSTA13wr+winlen);
                    else
                        dummy(1,:)=STAopt(1,iediff+1:winbig-mshift);
                        dummy(2,:)=STAopt(2,iediff+1-imaxSTA12wr:winbig-mshift-imaxSTA12wr);
                        dummy(3,:)=STAopt(3,iediff+1-imaxSTA13wr:winbig-mshift-imaxSTA13wr);
                    end
                    denoms=dot(dummy,dummy,2);
                    ccpost(nin+1)=(dot(dummy(1,:),dummy(2,:))/sqrt(denoms(1)*denoms(2))+dot(dummy(2,:),dummy(3,:))/sqrt(denoms(2)*denoms(3))+ ...
                                  dot(dummy(3,:),dummy(1,:))/sqrt(denoms(3)*denoms(1)))/3;
                    clear dummy
                    %Energy in prior 15 seconds
                    if isdiff > 15*sps+mshift
                        dummy=median(STAauto(1,isdiff-15*sps:isdiff))+ ...
                              median(STAauto(2,isdiff-15*sps-imaxSTA12wr:isdiff-imaxSTA12wr))+ ...
                              median(STAauto(3,isdiff-15*sps-imaxSTA13wr:isdiff-imaxSTA13wr));
                    else
                        dummy=median(STAauto(1,1:isdiff))+ ...
                              median(STAauto(2,1:isdiff-imaxSTA12wr))+ ...
                              median(STAauto(3,1:isdiff-imaxSTA13wr));
                    end
                    Prev15(nin+1)=dummy;
                    clear dummy
%                     %Energy in prior 30 seconds
%                     if isdiff > 30*sps+mshift
%                         dummy=median(STAauto(1,isdiff-30*sps:isdiff))+ ...
%                               median(STAauto(2,isdiff-30*sps-imaxSTA12wr:isdiff-imaxSTA12wr))+ ...
%                               median(STAauto(3,isdiff-30*sps-imaxSTA13wr:isdiff-imaxSTA13wr));
%                     else
%                         dummy=median(STAauto(1,1:isdiff))+ ...
%                               median(STAauto(2,1:isdiff-imaxSTA12wr))+ ...
%                               median(STAauto(3,1:isdiff-imaxSTA13wr));
%                     end
%                     Prev30(nin+1)=dummy;
%                     clear dummy
                    dummy=STAopt(1,iediff:iediff+2.5*cncntr).*STAopt(1,iediff:iediff+2.5*cncntr)+ ...
                          STAopt(2,iediff-imaxSTA12wr:iediff+2.5*cncntr-imaxSTA12wr).*STAopt(2,iediff-imaxSTA12wr:iediff+2.5*cncntr-imaxSTA12wr)+ ...
                          STAopt(3,iediff-imaxSTA13wr:iediff+2.5*cncntr-imaxSTA13wr).*STAopt(3,iediff-imaxSTA13wr:iediff+2.5*cncntr-imaxSTA13wr);
                    dum2=cumsum(dummy);
                    Post(nin+1)=dum2(end);
                    % %%%%%%Below is for comparison to template, + and -
                    indx=istart+idiff-1+round(cncntr/2); %indx should be centered ~ on zero-crossing of main arrival
                    traces(1,:)=STAopt(1,indx-tempbef:indx+tempaft);
                    traces(2,:)=STAopt(2,indx-tempbef-imaxSTA12wr:indx+tempaft-imaxSTA12wr);
                    traces(3,:)=STAopt(3,indx-tempbef-imaxSTA13wr:indx+tempaft-imaxSTA13wr);
                    tempxc=zeros(nsta,2*floor(cncntr/2)+1);
                    for ista=1:nsta
                        tempxc(ista,:)=xcorr(traces(ista,:),STAtemps(ista,:),floor(cncntr/2),'coeff');
                    end
                    sumxc=sum(tempxc)/nsta;
                    [match(nin+1,1) ioff]=max(sumxc);
                    ioff=ioff-(floor(cncntr/2)+1); %shift STAtemps by match(nin+1) (shift right for positive values)
                    timstemp(nin*templen+1:(nin+1)*templen)=timsSTA(indx-tempbef+ioff:indx+tempaft+ioff);
                    tempxcneg=zeros(nsta,2*floor(cncntr/2)+1);
                    for ista=1:nsta
                        tempxcneg(ista,:)=xcorr(traces(ista,:),-STAtemps(ista,:),floor(cncntr/2),'coeff');
                    end
                    sumxcneg=sum(tempxcneg)/nsta;
                    [match(nin+1,2) ioff]=max(sumxcneg);
                    ioff=ioff-(floor(cncntr/2)+1); %shift "snips" by match(nin+1) (shift right for positive values)
                    timstempneg(nin*templen+1:(nin+1)*templen)=timsSTA(indx-tempbef+ioff:indx+tempaft+ioff);
                    clear traces
% %%%%%%%%%%%%%%%%%%%%%%%%%%Above is for comparison to template, + and -
                    if stretchy
% %%%%%%%%%%%%%%%%%%%%%%%%%Below is for comparison to STRETCHED template, + only
                        %indx still good
                        traces(1,:)=STAopt(1,indx-stretchlen/2+1:indx+stretchlen/2); %"traces" are ~centered on zero-xing
                        traces(2,:)=STAopt(2,indx-stretchlen/2+1-imaxSTA12wr:indx+stretchlen/2-imaxSTA12wr);
                        traces(3,:)=STAopt(3,indx-stretchlen/2+1-imaxSTA13wr:indx+stretchlen/2-imaxSTA13wr);
                        maxccdipole=zeros(nsta,maxstretch-minstretch+1);
                        imaxdipole=zeros(nsta,maxstretch-minstretch+1);
                        tempxc=zeros(nsta,maxstretch-minstretch+1,2*floor(cncntr/2)+1);
                        for ista=1:nsta
                            for i=1:maxstretch-minstretch+1
                                dum=STAstr(ista,i,:); %dummy version of the stretched template.  Why not "squeezed"?
                                tempxc(ista,i,:)=xcorr(traces(ista,:),squeeze(dum)',floor(cncntr/2),'coeff'); %search for a best shift
    %                           [maxccdipole(ista,i),imaxdipole(ista,i)]=max(tempxc); %for each stretched template, the ccmax and its offset
                            end 
                        end
                        tempmat=max(tempxc,[],2); %this takes the max along the second, stretch dimension, and yields an array of size (nsta,1,cncntr+1)
                        tempmat2=squeeze(tempmat); %This should be a matrix of size (nsta,cncntr+1).
                        sumxc=sum(tempmat2)/nsta; %this should sum along the first, nsta dimension, and yield an array of size (1,cncntr+1)
                        [~,ibestshift]=max(sumxc); %ibestshift is approximate. Next, search for the best stretch and shift near this shift.
                        mibshift=4;
                        for ista=1:nsta 
                            if ibestshift>mibshift && 2*floor(cncntr/2)+1-ibestshift>mibshift
                                %tmpmat3(maxstretch-minstretch+1,2*mibshift+1)=squeeze(STAstr(ista,:,ibestshift-mibshift:ibestshift+mibshift)); %This excises a piece of STAstr around the "best" shift.
                                tmpmat3=squeeze(tempxc(ista,:,ibestshift-mibshift:ibestshift+mibshift)); %This excises a piece of STAstr around the "best" shift.
                                [maxmaxdipole(ista,nin+1),ind]=max(tmpmat3(:)); %tmpmat3(:) must turn tmpmat3 into a 1-D array.  maxmaxdipole is ccmax near the "best" shift.
                                [imaxstretch(ista,nin+1),imaxmaxdipole(ista,nin+1)]=ind2sub(size(tmpmat3),ind); %imaxstretch is the id of the best time-stretch; imaxmaxdipole is the best timeshift of that best stretch
                                imaxmaxdipole(ista,nin+1)=imaxmaxdipole(ista,nin+1)+(ibestshift-mibshift-1)-(floor(cncntr/2)+1); %correct imaxmaxdipole for ibestshift and floor(cncntr/2)
                            elseif ibestshift>mibshift
                                %tmpmat3(maxstretch-minstretch+1,ibestshift-mibshift:size(STAstr,3))=STAstr(ista,:,ibestshift-mibshift:end); %This excises a piece of STAstr around the "best" shift.
                                tmpmat3=squeeze(tempxc(ista,:,ibestshift-mibshift:end)); %This excises a piece of STAstr around the "best" shift.
                                [maxmaxdipole(ista,nin+1),ind]=max(tmpmat3(:)); %tmpmat3(:) must turn tmpmat3 into a 1-D array.  maxmaxdipole is ccmax near the "best" shift.
                                [imaxstretch(ista,nin+1),imaxmaxdipole(ista,nin+1)]=ind2sub(size(tmpmat3),ind); %imaxstretch is the id of the best time-stretch; imaxmaxdipole is the best timeshift of that best stretch
                                imaxmaxdipole(ista,nin+1)=imaxmaxdipole(ista,nin+1)+(ibestshift-mibshift-1)-(floor(cncntr/2)+1); %correct imaxsmaxdipole for ibestshift and floor(cncntr/2)
                            else 
                                %tmpmat3(maxstretch-minstretch+1,ibestshift+mibshift)=STAstr(ista,:,1:ibestshift+mibshift); %This excises a piece of STAstr around the "best" shift.
                                tmpmat3=squeeze(tempxc(ista,:,1:ibestshift+mibshift)); %This excises a piece of STAstr around the "best" shift.
                                [maxmaxdipole(ista,nin+1),ind]=max(tmpmat3(:)); %tmpmat3(:) must turn tmpmat3 into a 1-D array.  maxmaxdipole is ccmax near the "best" shift.
                                [imaxstretch(ista,nin+1),imaxmaxdipole(ista,nin+1)]=ind2sub(size(tmpmat3),ind); %imaxstretch is the id of the best time-stretch; imaxmaxdipole is the best timeshift of that best stretch
                                imaxmaxdipole(ista,nin+1)=imaxmaxdipole(ista,nin+1)-(floor(cncntr/2)+1); %correct imaxmaxdipole for floor(cncntr/2)
                            end 
                            tempsplotstr=squeeze(STAstr(ista,imaxstretch(ista,nin+1),:)); %set up a dummy version of the stretched dipole 
                            buffer=floor(cncntr/2)+1;
                            tempsplotstr(buffer:end-buffer)=tempsplotstr(buffer-imaxmaxdipole(ista,nin+1):end-buffer-imaxmaxdipole(ista,nin+1));
                            timstempstr(ista,nin*stretchlen+1:(nin+1)*stretchlen)= ...
                                timsSTA(indx-stretchlen/2+1+imaxmaxdipole(ista,nin+1):indx+stretchlen/2+imaxmaxdipole(ista,nin+1));
                            localmin=-min(traces(ista,:));
                            tempscaletemp=scaleSTA(traces(ista,:),tempsplotstr,localmin); %scale the template to match the trace
                            scaletemp(ista,nin+1)=tempscaletemp; %store it
                        end
                        stretches(:,nin+1)=(imaxstretch(:,nin+1)+(minstretch-1))/sps; 
                        clear tmpmat3
                        clear traces
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Above is for comparison to STRETCHED template, + only
                    end
                    % STA1file, STA2file, STA3file are back-to-back windows for the day (plus sample times)
                    STA1file(nin*winlen+1:(nin+1)*winlen,1:2)=[timsSTA(istart:iend)' STAopt(1,istart:iend)'];
                    STA2file(nin*winlen+1:(nin+1)*winlen,1:2)=[timsSTA(istart:iend)' STAopt(2,istart-imaxSTA12wr:iend-imaxSTA12wr)'];
                    STA3file(nin*winlen+1:(nin+1)*winlen,1:2)=[timsSTA(istart:iend)' STAopt(3,istart-imaxSTA13wr:iend-imaxSTA13wr)'];
                    STAamp(nin+1,1)=prctile(abs(STAopt(1,istart:iend)),80);
                    STAamp(nin+1,2)=prctile(abs(STAopt(2,istart-imaxSTA12wr:iend-imaxSTA12wr)),80);
                    STAamp(nin+1,3)=prctile(abs(STAopt(3,istart-imaxSTA13wr:iend-imaxSTA13wr)),80);
                    STAamp(nin+1,:)=STAamp(nin+1,:)/STAamp(nin+1,1);
%                     STA1bbfile(nin*winlen+1:(nin+1)*winlen,1:2)=[timsSTA(istart:iend)' STAoptbb(1,istart:iend)'];
%                     STA2bbfile(nin*winlen+1:(nin+1)*winlen,1:2)=[timsSTA(istart:iend)' STAoptbb(2,istart-imaxSTA12wr:iend-imaxSTA12wr)'];
%                     STA3bbfile(nin*winlen+1:(nin+1)*winlen,1:2)=[timsSTA(istart:iend)' STAoptbb(3,istart-imaxSTA13wr:iend-imaxSTA13wr)'];
                    STA12file(nin+1,1:2)=[imaxSTA12wr xcmaxSTA12n(n)];
                    STA13file(nin+1,1:2)=[imaxSTA13wr xcmaxSTA13n(n)];
                    STA32file(nin+1,1:3)=[cumsumtrdiff/cumsumtr(winlen) xcmaxSTA32n(n) idiff];

                    nin=nin+1;
                    istartkeep(nin)=istart; %For adding other stations later
                    aSTA12keep(nin,:)=[timswin(n) aSTA12(n)];
                    aSTA13keep(nin,:)=[timswin(n) aSTA13(n)];
                    aSTA32keep(nin,:)=[timswin(n) aSTA32(n)];
                    loopoffkeep(nin,:)=[timswin(n) loopoff(n)];
                    mapfile(nin,:)=[timswin(n) xmaxSTA13ntmp(n) xmaxSTA12ntmp(n) ...
                        xcmaxAVEn(n) loopoff(n) Ampsq(nin) ccarr(nin) timswin(n)-winlensec/2+idiff/sps cumsumtrdiff/cumsumtr(winlen) ...
                        match(nin,1) match(nin,2) Prev(nin) Post(nin) Prev15(nin) Prev2(nin) ...
                        ccprior125(nin) ccprior(nin) ccpost(nin) STAamp(nin,2) STAamp(nin,3) xcmaxSTA12n(n) xcmaxSTA13n(n) xcmaxSTA32n(n) ];
% %%%%%%%%%%%%%%%%%%%%%%%%%%Below is for comparison to STRETCHED template, + only
                    if stretchy
                        maxshiftdiff=max([abs(imaxmaxdipole(1,nin)-imaxmaxdipole(2,nin)), ...
                                          abs(imaxmaxdipole(1,nin)-imaxmaxdipole(3,nin)), ...
                                          abs(imaxmaxdipole(2,nin)-imaxmaxdipole(3,nin))]);
                        stretchfile(nin,:)=[timswin(n) timswin(n)-winlensec/2+idiff/sps xmaxSTA13ntmp(n) xmaxSTA12ntmp(n) xcmaxAVEn(n) ...
                            Ampsq(nin) cumsumtrdiff/cumsumtr(winlen) Prev(nin) ccarr(nin) ...
                            maxmaxdipole(1,nin) maxmaxdipole(2,nin) maxmaxdipole(3,nin) ... %these are cc coefficients
                            stretches(1,nin) stretches(2,nin) stretches(3,nin) ... %these are time stretches
                            scaletemp(1,nin) scaletemp(2,nin) scaletemp(3,nin) maxshiftdiff]; %these are amplitude scalings
                    end
% %%%%%%%%%%%%%%%%%%%%%%%%%%Above is for comparison to STRETCHED template, + only
                else
                    xmaxSTA12ntmp(n)=20; xmaxSTA13ntmp(n)=20; xmaxSTA32ntmp(n)=20;
                end
            else
                xmaxSTA12ntmp(n)=20; xmaxSTA13ntmp(n)=20; xmaxSTA32ntmp(n)=20; 
            end
        end
    end
    %fid = fopen(['3STAS/MAPS/2018test/map',IDENTIF,'_',num2str(lo),'-',num2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w'); %uncommented/compiled Jun 27  2019 wdwSYNTHstretch.m
    %fid = fopen(['3STAS/MAPS/2021/map',IDENTIF,'_',num2str(lo),'-',num2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
    fid = fopen(['3STAS23/MAPS/2023/map',IDENTIF,'_',num2str(lo),'-',num2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
    fprintf(fid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %8.3f %10.3f %7.3f %7.3f %7.3f %10.3e %10.3e %10.3e %10.3e %6.3f %6.3f %6.3f %5.2f %5.2f %7.3f %7.3f %7.3f\n',mapfile(1:nin,:)');
    fclose(fid);
% %%%%%%%%%%%%%%%%%%%%%%%%%%Below is for comparison to STRETCHED template, + only    
    if stretchy
        %fid = fopen(['3STAS/MAPS/2018test/str',IDENTIF,'_',num2str(lo),'-',num2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w'); %uncommented/compiled Jun 27  2019 wdwSYNTHstretch.m
        %fid = fopen(['3STAS/MAPS/2021/str',IDENTIF,'_',num2str(lo),'-',num2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
        fid = fopen(['3STAS23/MAPS/2023/str',IDENTIF,'_',num2str(lo),'-',num2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
        fprintf(fid,'%10.3f %10.3f %6.2f %6.2f %8.3f %10.3e %7.3f %10.3e %6.3f %6.3f %6.3f %6.3f %5.2f %5.2f %5.2f %7.3f %7.3f %7.3f %5i \n',stretchfile(1:nin,:)');
        fclose(fid);
    end
% %%%%%%%%%%%%%%%%%%%%%%%%%%Above is for comparison to STRETCHED template, + only
    figure 
    subplot(4,1,1,'align'); 
    hold on
    plot(timswin,xcmaxAVEnmin*mshift+zeros(nwin,1),'k:');
    plot(timsSTA(winlen:2*winlen),7+zeros(winlen+1,1),'k','linewidth',2);
    plot(timswin,zeros(nwin,1),'k:');
    plot(timswin,xcmaxAVEn*mshift,'g');
    plot(timswin,xmaxSTA12ntmp,'bs','MarkerSize',2);
    plot(timswin,xmaxSTA13ntmp,'ro','MarkerSize',2);
    plot(timswin,xmaxSTA32ntmp,'k*','MarkerSize',2);
    axis([0 timbig/2 -mshift mshift]);
    ylabel('samples')
    box on
    subplot(4,1,2,'align'); 
    hold on
    plot(timswin,xcmaxAVEnmin*mshift+zeros(nwin,1),'k:');
    plot(timswin,zeros(nwin,1),'k:');
    plot(timswin,xcmaxAVEn*mshift,'g');
    plot(timswin,xmaxSTA12ntmp,'bs','MarkerSize',2);
    plot(timswin,xmaxSTA13ntmp,'ro','MarkerSize',2);
    plot(timswin,xmaxSTA32ntmp,'k*','MarkerSize',2);
    axis([timbig/2 timbig -mshift mshift]);
    ylabel('samples')
    box on
    subplot(4,1,3,'align'); 
    hold on
    plot(timswin,xcmaxAVEnmin*mshift+zeros(nwin,1),'k:');
    plot(timswin,zeros(nwin,1),'k:');
    plot(timswin,xcmaxAVEn*mshift,'g');
    plot(timswin,xmaxSTA12ntmp,'bs','MarkerSize',2);
    plot(timswin,xmaxSTA13ntmp,'ro','MarkerSize',2);
    plot(timswin,xmaxSTA32ntmp,'k*','MarkerSize',2);
    axis([timbig 3*timbig/2 -mshift mshift]);
    ylabel('samples')
    box on
    subplot(4,1,4,'align'); 
    hold on
    plot(timswin,xcmaxAVEnmin*mshift+zeros(nwin,1),'k:');
    plot(timswin,zeros(nwin,1),'k:');
    plot(timswin,xcmaxAVEn*mshift,'g');
    plot(timswin,xmaxSTA12ntmp,'bs','MarkerSize',2);
    plot(timswin,xmaxSTA13ntmp,'ro','MarkerSize',2);
    plot(timswin,xmaxSTA32ntmp,'k*','MarkerSize',2);
    axis([3*timbig/2 2*timbig -mshift mshift]);
    xlabel('sec')
    ylabel('samples')
    box on
    title([IDENTIF,'_{',num2str(lo),'-',num2str(hi),'}'])
    orient landscape
    %print('-depsc',['3STAS/FIGS/2018test/',IDENTIF,'-',int2str(winlen/sps),'s_',num2str(lo),'-',num2str(hi),'b.eps'])
    print('-depsc',['3STAS23/FIGS/2023/',IDENTIF,'-',int2str(winlen/sps),'s_',num2str(lo),'-',num2str(hi),'b.eps'])

    figure
    colormap(jet)
    scatter(xmaxSTA13n-xmaxSTA12n+xmaxSTA32n,xcmaxAVEn,3,log(ampopt./amport))
    hold on 
    plot(-50:50,xcmaxAVEnmin+zeros(101,1),'k:');
    axis([min(-5,-2.5*loopoffmax) max(5,2.5*loopoffmax) -0.0 1.0])
    hrf = plotreflinesr(gca,-loopoffmax,'x','k');colorbar
    hrf = plotreflinesr(gca,loopoffmax,'x','k');colorbar
    box on
    title([IDENTIF,'_{',num2str(lo),'-',num2str(hi),'}'])
    %print('-depsc',['3STAS/FIGS/2018test/',IDENTIF,'-',int2str(winlen/40),'s_',num2str(lo),'-',num2str(hi),'e.eps'])
    print('-depsc',['3STAS23/FIGS/2023/',IDENTIF,'-',int2str(winlen/40),'s_',num2str(lo),'-',num2str(hi),'e.eps'])

    if winlen>500
    scrsz=get(0,'ScreenSize');
    nt=0;
    nset=2; %# of panels per detection
    nrow=6;
    mcol=6; %2 for 20s;
    figmax=12;
    figmax=min(figmax,floor(nin/(nrow*mcol))+1);
    for ifig=1:floor(nin/(nrow*mcol))+1
        figure('Position',[scrsz(3)/10 scrsz(4)/10 4*scrsz(3)/5 9*scrsz(4)/10]);
        fign(ifig)=gcf;
        if ifig > 5
            close(fign(ifig-5))
        end
        for n = 1:nrow
            for m = 1:mcol
                nt=nt+1;
                if nt <= nin
                     %if STA12file(nt,1) >= 10 && STA12file(nt,1) <= 16 && STA13file(nt,1) >= 2 && STA13file(nt,1) <= 8
                        subplot(nset*nrow,mcol,nset*(n-1)*mcol+m,'align');
                        yma=max(max([STA1file(winlen*(nt-1)+1:winlen*nt,2) STA2file(winlen*(nt-1)+1:winlen*nt,2) ...
                            STA3file(winlen*(nt-1)+1:winlen*nt,2)]));
                        ymi=min(min([STA1file(winlen*(nt-1)+1:winlen*nt,2) STA2file(winlen*(nt-1)+1:winlen*nt,2) ...
                            STA3file(winlen*(nt-1)+1:winlen*nt,2)]));
                        yma=2.4*max(yma,-ymi);
% Lines below plot + or - template:
                        if match(nt,1) >= match(nt,2) %if positive cc with template is larger ...
                            plot(timstemp((nt-1)*templen+1:nt*templen),STAtemps(whichtoplot,:)*yma/2.4,'c','linewidth',2)
                        else %if cc with negative template is larger ...
                            plot(timstempneg((nt-1)*templen+1:nt*templen),-STAtemps(whichtoplot,:)*yma/2.4,'m','linewidth',2)
                        end
% Lines above plot + or - template:
% Lines below plot + stretched template:
%                       tempplot=squeeze(STAstr(whichtoplot,imaxstretch(whichtoplot,nt),:))*scaletemp(whichtoplot,nt);
%                       plot(timstempstr(whichtoplot,(nt-1)*stretchlen+1:nt*stretchlen),tempplot,'c','linewidth',1)
% Lines above plot + stretched template:                        
                        hold on
                        plot(STA1file(winlen*(nt-1)+1:winlen*nt,1),STA1file(winlen*(nt-1)+1:winlen*nt,2),'r')
                        plot(STA2file(winlen*(nt-1)+1:winlen*nt,1),STA2file(winlen*(nt-1)+1:winlen*nt,2),'b')
                        plot(STA3file(winlen*(nt-1)+1:winlen*nt,1),STA3file(winlen*(nt-1)+1:winlen*nt,2),'k')
                        is = STA1file(winlen*(nt-1)+1,1);
                        ien= STA1file(winlen*nt,1);
                        issamp=round(is*sps);
                        iensamp=round(ien*sps);
                        plot(STA1file(winlen*(nt-1)+1:winlen*nt,1),medhil(issamp:iensamp),'g') %just to test Hilbert transform
                        axis([is ien -yma yma])
%                         xvect=[is is+2*(yma-ymi)*(winlen/160.)]; %*mean(scalefact)]; %amplitude bar originally scaled for 4-s window.  Not sure "mean" is necessary.
%                         yvect=[-0.9*yma -0.9*yma];
%                         plot(xvect,yvect,'r','linewidth',3)
                        plot([is+1/hi is+1/lo],[-0.9*yma -0.9*yma],'k','linewidth',3)
                        xcvect=[is+STA32file(nt,3)/sps is+(STA32file(nt,3)+cncntr)/sps];
                        ycvect=[0.95*yma 0.95*yma];
                        plot(xcvect,ycvect,'b','linewidth',3)
                        %plot(bostsec,0.93*yma,'ro','MarkerSize',4,'MarkerFaceColor','r')
                        text(is+STA32file(nt,3)/sps,  0.58*yma, num2str(match(nt,1),2),'fontsize',6,'color','b');
                        text(is+STA32file(nt,3)/sps, -0.75*yma, num2str(match(nt,2),2),'fontsize',6,'color','r');
                        text(is+0.2, 0.66*yma, int2str(STA13file(nt,1)),'fontsize',6);
                        text(ien-0.6, 0.66*yma, int2str(STA12file(nt,1)),'fontsize',6);
                        box on
                        set(gca,'XTick',[is (is+ien)/2],'fontsize',6);
                        % %%%%%% Below: Find the max of the hilbert transform within the arrival; then walk out from there +/- 
                        [maxhil loc]=max(medhil(issamp+STA32file(nt,3):issamp+STA32file(nt,3)+cncntr-1));
                        maxloc=issamp+STA32file(nt,3)+loc;
                        icheck=maxloc;
                        while medhil(icheck)>0.5*maxhil && icheck<tracelen-1
                            icheck=icheck+1;
                        end
                        ccdat(nt,13)=icheck-maxloc;
                        while medhil(icheck)>0.25*maxhil && icheck<tracelen-1
                            icheck=icheck+1;
                        end
                        ccdat(nt,14)=icheck-maxloc;
                        while medhil(icheck)>0.125*maxhil && icheck<tracelen-1
                            icheck=icheck+1;
                        end
                        ccdat(nt,15)=icheck-maxloc;
                        icheck=maxloc;
                        while medhil(icheck)>0.5*maxhil && icheck>1
                            icheck=icheck-1;
                        end
                        ccdat(nt,16)=maxloc-icheck;
                        while medhil(icheck)>0.25*maxhil && icheck>1
                            icheck=icheck-1;
                        end
                        ccdat(nt,17)=maxloc-icheck;
                        while medhil(icheck)>0.125*maxhil && icheck>1
                            icheck=icheck-1;
                        end
                        ccdat(nt,18)=maxloc-icheck;
                        % This for longer-term averages:
                        maxhil=median(medhil(issamp+STA32file(nt,3):issamp+STA32file(nt,3)+cncntr-1)); %median of arrival window
                        maxloc=issamp+STA32file(nt,3); %"maxloc" is really start of arrival window
                        icheck=maxloc;
                        while median(medhil(icheck:icheck+cncntr-1))>0.666*maxhil && icheck+cncntr<tracelen && icheck-maxloc < 240*sps
                            icheck=icheck+1;
                        end
                        ccdat(nt,19)=icheck-maxloc;
                        while median(medhil(icheck:icheck+cncntr-1))>0.5*maxhil && icheck+cncntr<tracelen && icheck-maxloc < 240*sps
                            icheck=icheck+1;
                        end
                        ccdat(nt,20)=icheck-maxloc;
                        while median(medhil(icheck:icheck+cncntr-1))>0.25*maxhil && icheck+cncntr<tracelen && icheck-maxloc < 240*sps
                            icheck=icheck+1;
                        end
                        ccdat(nt,21)=icheck-maxloc;
                        maxloc=issamp+STA32file(nt,3)+cncntr-1; %"maxloc" is really end of arrival window
                        icheck=maxloc;
                        while median(medhil(icheck-cncntr+1:icheck))>0.666*maxhil && icheck-cncntr>0 && maxloc-icheck < 240*sps
                            icheck=icheck-1;
                        end
                        ccdat(nt,22)=maxloc-icheck;
                        while median(medhil(icheck-cncntr+1:icheck))>0.5*maxhil && icheck-cncntr>0 && maxloc-icheck < 240*sps
                            icheck=icheck-1;
                        end
                        ccdat(nt,23)=maxloc-icheck;
                        while median(medhil(icheck-cncntr+1:icheck))>0.25*maxhil && icheck-cncntr>0 && maxloc-icheck < 240*sps
                            icheck=icheck-1;
                        end
                        ccdat(nt,24)=maxloc-icheck;
                        % %%%%%% Above: Find the max of the hilbert transform within the arrival; then walk out from there +/- 

%                         subplot(3*nrow,mcol,3*(n-1)*mcol+mcol+m,'align');
% %                         plot(STA1bbfile(winlen*(nt-1)+1:winlen*nt,1),STA1bbfile(winlen*(nt-1)+1:winlen*nt,2),'r')
% %                         hold on
% %                         plot(STA2bbfile(winlen*(nt-1)+1:winlen*nt,1),STA2bbfile(winlen*(nt-1)+1:winlen*nt,2),'b')
% %                         plot(STA3bbfile(winlen*(nt-1)+1:winlen*nt,1),STA3bbfile(winlen*(nt-1)+1:winlen*nt,2),'k')
% %                         is = STA1bbfile(winlen*(nt-1)+1,1);
% %                         ien= STA1bbfile(winlen*nt,1);
% % %                         yma=max(max([STA1bbfile(winlen*(nt-1)+1:winlen*nt,2) STA2bbfile(winlen*(nt-1)+1:winlen*nt,2) ...
% % %                             STA3bbfile(winlen*(nt-1)+1:winlen*nt,2)]));
% % %                         ymi=min(min([STA1bbfile(winlen*(nt-1)+1:winlen*nt,2) STA2bbfile(winlen*(nt-1)+1:winlen*nt,2) ...
% % %                             STA3bbfile(winlen*(nt-1)+1:winlen*nt,2)]));
% % %                         xvect=[is is+2*(yma-ymi)*(winlen/160.)]; %amplitude bar originally scaled for 4-s window
% % %                         yma=2.4*max(yma,-ymi);
% % %                         yvect=[-0.9*yma -0.9*yma];
% % %                         plot(xvect,yvect,'r','linewidth',3)
% %                         plot([is+1/hibb is+1/lobb],[-0.8*yma -0.8*yma],'k','linewidth',3)
%                         text(is+0.2, 0.66*yma, int2str(STA13file(nt,1)),'fontsize',6);
%                         text(ien-0.6, 0.66*yma, int2str(STA12file(nt,1)),'fontsize',6);
%                         box on
%                         axis([is ien -yma yma])
%                         set(gca,'XTick',[is (is+ien)/2],'fontsize',6);

                    subplot(nset*nrow,mcol,nset*(n-1)*mcol+mcol+m,'align');
                    STA12tr=STA1file(winlen*(nt-1)+1:winlen*nt,2).*STA2file(winlen*(nt-1)+1:winlen*nt,2);            
                    STA13tr=STA1file(winlen*(nt-1)+1:winlen*nt,2).*STA3file(winlen*(nt-1)+1:winlen*nt,2);
                    STA23tr=STA2file(winlen*(nt-1)+1:winlen*nt,2).*STA3file(winlen*(nt-1)+1:winlen*nt,2);
%                     % find the peaks etc.
%                     avedots=(STA12tr+STA13tr+STA23tr)/3.;
%                     [peaks, locs]=findpeaks(avedots,'minpeakdistance',3);
%                     npeaks=length(peaks);
%                     [maxpk, imaxpk]=max(peaks);
%                     pks=zeros(npeaks,2);
%                     pks(:,2)=peaks;
%                     pks(:,1)=STA1file(winlen*(nt-1)+locs);
%                     pksort=sortrows(pks,2);
%                     rat14=maxpk/pksort(npeaks-4,2); %ratio of max to 4th largest anywhere in window
%                     if imaxpk==1 
%                         maxpkp=peaks(2);
%                         maxpkm=-9e9;
%                         pkwid=2*(pks(2,1)-pks(1,1));
%                         pksid12=maxpk/maxpkp;
%                         pksid13=maxpk/peaks(3);
%                     elseif imaxpk==npeaks
%                         maxpkp=-9e9;
%                         maxpkm=peaks(npeaks-1);
%                         pkwid=2*(pks(npeaks,1)-pks(npeaks-1,1));
%                         pksid12=maxpk/maxpkm;
%                         pksid13=maxpk/peaks(npeaks-2);
%                     else
%                         maxpkp=peaks(imaxpk+1);
%                         maxpkm=peaks(imaxpk-1);
%                         pkwid=pks(imaxpk+1,1)-pks(imaxpk-1,1);
%                         pksid12=maxpk/max(maxpkp,maxpkm);
%                         pksid13=maxpk/min(maxpkp,maxpkm);
%                     end
                    cumsumSTA12tr=cumsum(STA12tr);
                    cumsumSTA13tr=cumsum(STA13tr);
                    cumsumSTA23tr=cumsum(STA23tr);
%                     yma=1.1; %yma=max(avedots);
%                     ymi=-0.1; %ymi=min(avedots);
% The following for running cc
                    cclen=20; %running cc window length, in samples
                    yma=max(max([cumsumSTA12tr cumsumSTA13tr cumsumSTA23tr]));
                    ymi=min(min([cumsumSTA12tr cumsumSTA13tr cumsumSTA23tr]));
                    ST1auto=STA1file(winlen*(nt-1)+1:winlen*nt,2).*STA1file(winlen*(nt-1)+1:winlen*nt,2);
                    ST1sq=cumsum(ST1auto);
                    ST2auto=STA2file(winlen*(nt-1)+1:winlen*nt,2).*STA2file(winlen*(nt-1)+1:winlen*nt,2);
                    ST2sq=cumsum(ST2auto);
                    ST3auto=STA3file(winlen*(nt-1)+1:winlen*nt,2).*STA3file(winlen*(nt-1)+1:winlen*nt,2);
                    ST3sq=cumsum(ST3auto);
                    ST12num=cumsumSTA12tr(cclen+1:winlen)-cumsumSTA12tr(1:winlen-cclen);
                    ST13num=cumsumSTA13tr(cclen+1:winlen)-cumsumSTA13tr(1:winlen-cclen);
                    ST23num=cumsumSTA23tr(cclen+1:winlen)-cumsumSTA23tr(1:winlen-cclen);
                    ST1den=ST1sq(cclen+1:winlen)-ST1sq(1:winlen-cclen);
                    ST2den=ST2sq(cclen+1:winlen)-ST2sq(1:winlen-cclen);
                    ST3den=ST3sq(cclen+1:winlen)-ST3sq(1:winlen-cclen);
                    ST12n=ST12num./realsqrt(ST1den.*ST2den);
                    ST13n=ST13num./realsqrt(ST1den.*ST3den);
                    ST23n=ST23num./realsqrt(ST2den.*ST3den);
                    alln=(ST12n+ST13n+ST23n)/3;

                    meanalln=mean(alln);
                    ccdat(nt,1)=0;
                    ccdat(nt,2)=0;
                    ccdat(nt,3)=ifig;
                    ccdat(nt,4)=STA1file(winlen*(nt-1)+1,1);
                    ccdat(nt,5)=STA13file(nt,1);
                    ccdat(nt,6)=STA12file(nt,1);
                    ccdat(nt,7)=Prev(nt);
                    ccdat(nt,8)=(STA13file(nt,2)+STA12file(nt,2)+STA32file(nt,2))/3;
                    ccdat(nt,9)=meanalln;
                    ccdat(nt,10)=STA32file(nt,1);
                    ccdat(nt,11)=2.5*Ampsq(nt)/Prev(nt);
                    ccdat(nt,12)=match(nt,1)-match(nt,2);

                    alln(alln<0)=-10^4*yma; %just so they don't plot.
                    % %idiff=STA32file(nt,3);
                    % %maxxc=max(alln(idiff:idiff+cclen)); %endpoints of alln window should be endpoints of 1 sec cumsum interval
                    plot(STA1file(winlen*(nt-1)+cclen/2+1:winlen*nt-cclen/2,1),(yma+ymi)/2+(yma-ymi)*alln/2,'co','markersize',1)
                    hold on
                    %plot(STA1file(winlen*(nt-1)+cclen/2+1:winlen*nt-cclen/2,1),ST12n,'b')
                    %plot(STA1file(winlen*(nt-1)+cclen/2+1:winlen*nt-cclen/2,1),ST13n,'r')
                    %plot(STA1file(winlen*(nt-1)+cclen/2+1:winlen*nt-cclen/2,1),ST23n,'k')
                    %%plot(STA1file(winlen*(nt-1)+1:winlen*nt,1),0.75*maxxc*ones(winlen,1),'k:');
                    %%plot(STA1file(winlen*(nt-1)+1:winlen*nt,1),0.65*maxxc*ones(winlen,1),'k:');
% The above for running cc
% The following for running-sum dot product
                    plot(STA1file(winlen*(nt-1)+1:winlen*nt,1),cumsumSTA12tr,'k') %Use the color of the excluded station
                    hold on
                    plot(STA1file(winlen*(nt-1)+1:winlen*nt,1),cumsumSTA13tr,'b')
                    plot(STA1file(winlen*(nt-1)+1:winlen*nt,1),cumsumSTA23tr,'r')
                    axis([is ien ymi yma])
                    set(gca,'XTick',(0:20),'fontsize',6);
% The above for running-sum dot product
                    box on
                    text(is+0.1, ymi+0.82*(yma-ymi), num2str(STA13file(nt,2),2),'fontsize',6);
                    text(is+0.1, ymi+0.64*(yma-ymi), num2str(STA12file(nt,2),2),'fontsize',6);
                    text(is+0.1, ymi+0.46*(yma-ymi), num2str(STA32file(nt,2),2),'fontsize',6);
                    text(ien-0.6, ymi+0.1*(yma-ymi), num2str(STA32file(nt,1),2),'fontsize',6);
                    text(ien-2.2, ymi+0.1*(yma-ymi), num2str(STAamp(nt,2),2),'fontsize',6);
                    text(ien-1.4, ymi+0.1*(yma-ymi), num2str(STAamp(nt,3),2),'fontsize',6);
                    %axis([lo-1 hi+1 ymi yma])
                    %axis tight
                    %set(gca,'XTick',[is is+2],'fontsize',6);
                    %pkfile(nt,:)=[STA1file(winlen*(nt-1)+1+(winlen/2)) pks(imaxpk,1) maxpk pksid12 pksid13 pkwid rat14];
                end
            end
        end
        drawnow
        orient landscape
%        if ifig <= 9
%            print('-depsc',['3STAS/WIGS/2018test/',IDENTIF,'-',num2str(lo),'-',num2str(hi),'ms',int2str(mshift),'-',int2str(winlen/sps),'s','.',int2str(0),int2str(0),int2str(ifig),'.eps'])
%        elseif ifig <= 15 %99
%            print('-depsc',['3STAS/WIGS/2018test/',IDENTIF,'-',num2str(lo),'-',num2str(hi),'ms',int2str(mshift),'-',int2str(winlen/sps),'s','.',int2str(0),int2str(ifig),'.eps'])
% %         else
% %             print('-depsc',['3STAS/WIGS/2018test/',IDENTIF,'-',num2str(lo),'-',num2str(hi),'ms',int2str(mshift),'-',int2str(winlen/sps),'s','.',int2str(ifig),'.eps'])
%        end
        if ifig <= 9
            %print('-depsc',['3STAS/WIGS/2021/',IDENTIF,'-',num2str(lo),'-',num2str(hi),'ms',int2str(mshift),'-',int2str(winlen/sps),'s','.',int2str(0),int2str(0),int2str(ifig),'.eps'])
            print('-depsc',['3STAS23/WIGS/2023/',IDENTIF,'-',num2str(lo),'-',num2str(hi),'ms',int2str(mshift),'-',int2str(winlen/sps),'s','.',int2str(0),int2str(0),int2str(ifig),'.eps'])
%        elseif ifig <= 15 %99
%            print('-depsc',['3STAS/WIGS/2021/',IDENTIF,'-',num2str(lo),'-',num2str(hi),'ms',int2str(mshift),'-',int2str(winlen/sps),'s','.',int2str(0),int2str(ifig),'.eps'])
        end

    end
    %fid = fopen(['ARMMAP/MAPS/pks',IDENTIF,'_',num2str(lo),'-',num2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
%     fid = fopen(['3STAS/MAPS/2018/pks',IDENTIF,'_',num2str(lo),'-',num2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/sps),'s'],'w');
%     fprintf(fid,'%9.1f %9.3f %10.6f %9.3f %9.3f %9.3f %9.3f \n',pkfile(1:nin,:)');
%     fclose(fid);
    %fid = fopen(['3STAS/MAPS/2018test/ccdat_',IDENTIF,'_',num2str(lo),'-',num2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w'); %uncommented/compiled Jun 27  2019 wdwSYNTHstretch.m
    %fid = fopen(['3STAS/MAPS/2021/ccdat_',IDENTIF,'_',num2str(lo),'-',num2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
    fid = fopen(['3STAS23/MAPS/2023/ccdat_',IDENTIF,'_',num2str(lo),'-',num2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
    fprintf(fid,'%5i %4i %4i %11.4e %4i %4i %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e  %6i %5i %5i  %7i %5i %5i  %8i %6i %6i  %8i %6i %6i \n',ccdat(1:nin,:)');
    fclose(fid);
    clear STA1file; clear STA2file; clear STA3file; 
    end %(if winlen>500)
end

medlok=median(abs(loopoffkeep))
medaSTA12=median(aSTA12keep)
medaSTA13=median(aSTA13keep)
medaSTA32=median(aSTA32keep)
