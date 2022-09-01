%An autocorrelation.
format short e
close all
clear all
set(0,'DefaultFigureVisible','on');
scrsz=get(0,'ScreenSize');
h=figure('Position',[1 2.5*scrsz(4)/4 scrsz(3)/1.1 scrsz(4)/4]);

days=[2005 255];
tims=[%990 1040; %2005.255 (postdates 30-s)
      %5300 5800;
      %8050 8150;
      %35680 36370]; %nwins 1337
      19000 19250;
      22980 23030;
      51030 51550;
      54200 54360;
      58800 59600];
%       62200 62700;
%       68100 68700;
%       75400 75800];
ntims=size(tims,1);

STAs=['PGC ';
      'LZB ';
      %'SNB ';
      %'VGZ ';
      'SSIB';
      'SILB';
      'KLNB';
      'MGCB';
      'TWKB'];
%     'TSJB';
%     'TWBB'];
nsta=size(STAs,1);
STAoffs=[0;
         13;
         %205;
         %-43;
         87;
         20;
         -5;
         -22;
         66];

% 002
% STArots=[0 90 32 00;  %PGC , Yajun's "0" changed to 90.
%          6 85 33 86;  %SSIB from Yajun
%          0 90 39 20];  %SILB
%        0 90 54 00;  %LZB
%        0 90  7 -5;  %KLNB
%        4 70 48 -4;  %MGCB
%        4 75 38 -27]; %TWKB
%STArots(:,2:3)=pi*STArots(:,2:3)/180.;

hi=8.;
lo=3.;
npo=2;
npa=1;
winlen=3*40; %4 for 1.5-6Hz; 8 for 0.5-1.5Hz
winoff=1*40; %1 for 4s; 2 for 8 s; 3 for 12 s; 8 for 128s

t=cputime; tt=cputime;
year=days(1);
YEAR=int2str(year);
jday=days(2);
if jday <= 9
    JDAY=['00',int2str(jday)];
elseif jday<= 99
    JDAY=['0',int2str(jday)];
else
    JDAY=int2str(jday);
end
MO=day2month(jday,year)
%IDENTIF=[YEAR,'.',JDAY,'.',int2str(STArots(2,4)),'.',int2str(STArots(3,4)), ...
%    '.loff',num2str(loopoffmax),'.ccmin',num2str(xcmaxAVEnmin),'.nponpa',int2str(npo),int2str(npa)]

%Get data:
direc=[YEAR,'/',MO,'/'];
prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];
for ista=1:nsta
    if strcmp(STAs(ista,end),' ')
        sps=40;
        STAEdat=[prename,'.',STAs(ista,1:3),'..BHE.D.SAC']; %BHE for permstas.
        STANdat=[prename,'.',STAs(ista,1:3),'..BHN.D.SAC'];
%       spsfactor=1; %for glitches
    else
        sps=100;
        STAEdat=[prename,'.',STAs(ista,:),'..HHE.D.SAC']; %BHE for permstas.
        STANdat=[prename,'.',STAs(ista,:),'..HHN.D.SAC'];
%       spsfactor=2.5; %for glitches
    end
    %%%%%%%%%%%%%%%%
    %Read data:
    [STAE,HdrDataSTA,tnuSTA,pobjSTA,timsSTA]=readsac(STAEdat,0,'l');
    [STAN,~,~,~,~]=readsac(STANdat,0,'l');
    tracelen=length(STAE)
    if ista==1
        tracelenPGC=tracelen;
        timsPGC=timsSTA;
    end
    %%%%%%%%%%%%%%%%
    %Find glitches prior to filtering...
%   nzerosE=glitches2(STAE,nwin,winlen,winoff,igstart,spsfactor); 
%   nzerosN=glitches2(STAN,nwin,winlen,winoff,igstart,spsfactor);
%   nzeros(ista,:)=max(nzerosE,nzerosN);
    %%%%%%%%%%%%%%%%
    %cosine taper before filtering:
    x=(0:pi/(2*sps):pi/2-pi/(2*sps))';   
    if strcmp(STAs(ista,end),' ') %This seems to be necessary at the start of each day for PGCE:
        STAE(1:80)=0.;
        STAE(81:120)=sin(x).*STAE(81:120); %Only at start of day!
    else
        STAE(1:sps)=sin(x).*STAE(1:sps);
    end
    STAN(1:sps)=sin(x).*STAN(1:sps);
    x=flipud(x);
    STAE(tracelen-(sps-1):tracelen)=sin(x).*STAE(tracelen-(sps-1):tracelen);
    STAN(tracelen-(sps-1):tracelen)=sin(x).*STAN(tracelen-(sps-1):tracelen);
    %%%%%%%%%%%%%%%%
    %Filter data:
    if strcmp(STAs(ista,end),' ') 
        [STAE]=1.6e-4*bandpass(STAE,sps,lo,hi,npo,npa,'butter');
        [STAN]=1.6e-4*bandpass(STAN,sps,lo,hi,npo,npa,'butter');
    else
        if year==2003 && jday<213
            [STAE]=20.0e-3*bandpass(STAE,sps,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
            [STAN]=20.0e-3*bandpass(STAN,sps,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
        else
            [STAE]=4.0e-3*bandpass(STAE,sps,lo,hi,npo,npa,'butter'); 
            [STAN]=4.0e-3*bandpass(STAN,sps,lo,hi,npo,npa,'butter'); 
        end
        STAE=resample(STAE,2,5);
        STAN=resample(STAN,2,5);
        tracelen=length(STAE)
    end
%   %%%%%%%%%%%%%%%%
%   %Split-correct:
%   STA=STAEf+1i*STANf;
%   STAfastslow=STA*exp(-1i*STArots(ista,2));
%   STAslow=real(STAfastslow);
%   STAfast=imag(STAfastslow);
%   len=length(STA);
%   STAslow(10:tracelen-10)=STAslow(10+STArots(ista,1):tracelen-10+STArots(ista,1));
%   STAsplitcorrected=(STAslow+1i*STAfast)*exp(1i*STArots(ista,2));
%   STAsplitcorrectedrot=STAsplitcorrected*exp(-1i*STArots(ista,3));
%   STAscrot=STAsplitcorrectedrot;
%   %%%%%%%%%%%%%%%%
    %Center on desired region:
    STAoff=STAoffs(ista);
    if STAoff > 0
        STAE(1:tracelen-STAoff)=STAE(STAoff+1:tracelen);
        STAE(tracelen-STAoff+1:tracelen)=0;
        STAN(1:tracelen-STAoff)=STAN(STAoff+1:tracelen);
        STAN(tracelen-STAoff+1:tracelen)=0;
    elseif STAoff < 0
        STAE(-STAoff+1:tracelen)=STAE(1:tracelen+STAoff);
        STAE(1:-STAoff)=0;
        STAN(-STAoff+1:tracelen)=STAN(1:tracelen+STAoff);
        STAN(1:-STAoff)=0;
    end
%   STAopt(ista,:)=real(STAscrot);
%   STAort(ista,:)=imag(STAscrot);

    nsamp=0*1:ntims;
    itstart=0*1:ntims;
    sofar=0;
    for itim=1:ntims
        nsamp(itim)=40*(tims(itim,2)-tims(itim,1));
        itstart(itim)=40*tims(itim,1);
        itend=itstart(itim)+nsamp(itim)-1;
        timeseriesE(sofar+1:sofar+nsamp(itim))=STAE(itstart(itim):itend);
        timeseriesN(sofar+1:sofar+nsamp(itim))=STAN(itstart(itim):itend);
        sofar=sofar+nsamp(itim);
    end
    tslen=sofar
    tsE=timeseriesE;
    tsN=timeseriesN;
    tot=cputime-tt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Autocorrelation of stations.  Those that end in "2" are the running
    %   cumulative sum, to be used later by differncing the window edpoints.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    autoE=tsE.*tsE;
    autoN=tsN.*tsN;
    E2=cumsum(autoE);
    N2=cumsum(autoN);
    %PGCauto=realPGC.*realPGC;
    %PGC2=cumsum(PGCauto);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Cross-correlation between stations, with small offsets up to +/- mshift.
    %  First index is pointwise multiplication of traces; second is shifting offset.
    %  lenx is shorter than tracelenPGC by mshift at each end (see notebook sketch)
    %  For PGSS and PGSI, SSI and SIL are shifted relative to PGC, by 1 each time through loop.
    %  For SISS, SSI is shifted relative to SILB.
    %  cumsumPGSS etc. are the running cumulative sum of the x-correlation.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nwins=floor((tslen-winlen)/winoff);
    ncomps=tslen-winlen-1;  %-1 unnecessary?
    if ista==1
        ccE=zeros(nwins,ncomps);
        ccN=zeros(nwins,ncomps);
    end
    for iwin=1:nwins
        iwstart=(iwin-1)*winoff+1;
        iwend=iwstart+winlen-1;
        shortE=tsE(iwstart:iwend);
        shortN=tsN(iwstart:iwend);
        if iwin==1
            shortEvar=E2(iwend);  %Let's assume zero mean.
            shortNvar=N2(iwend);
        else
            shortEvar=E2(iwend)-E2(iwstart-1);
            shortNvar=N2(iwend)-N2(iwstart-1);
        end
        for icomp=1:ncomps
            itarstart=icomp;
            itarend=itarstart+winlen-1;
            targetE=tsE(itarstart:itarend);
            targetN=tsN(itarstart:itarend);
            if icomp==1
                targetEvar=E2(itarend);  %Let's again assume zero mean
                targetNvar=N2(itarend);
            else
                targetEvar=E2(itarend)-E2(itarstart-1);
                targetNvar=N2(itarend)-N2(itarstart-1);
            end
            auto=shortE.*targetE;
            cumauto=cumsum(auto);
            ccE(iwin,icomp)=ccE(iwin,icomp)+cumauto(end)/sqrt(shortEvar*targetEvar);
            auto=shortN.*targetN;
            cumauto=cumsum(auto);
            ccN(iwin,icomp)=ccN(iwin,icomp)+cumauto(end)/sqrt(shortNvar*targetNvar);
        end
        %h=figure('Position',[1 2.5*scrsz(4)/4 scrsz(3)/1.1 scrsz(4)/4]);
        plot((ccE(iwin,:)+ccN(iwin,:))/(2*ista))
        ylim([-1.1 1.1])
        drawnow
        %axis=([0 ncomps -1.1 1.1]);
    end
    tot=cputime-tt
end
cc=(ccE+ccN)/(2*nsta);
npks=50;
bbb=zeros(nwins,npks);
aaa=zeros(nwins,npks);
for i=1:nwins
    bb=findpeaks(cc(i,:),'SORTSTR','descend');
    bbb(i,:)=bb(2:npks+1);
    aa=findpeaks(-cc(i,:),'SORTSTR','descend');
    aaa(i,:)=aa(2:npks+1);
end
[bbbb ibbbb]=sort(bbb(:,1),'descend');
[aaaa iaaaa]=sort(aaa(:,1),'descend');
figure
hold on
for n=1:nwins
    plot(bbb(n,:))
end
   
%plot(cc(iwin,:))
tot=cputime-tt