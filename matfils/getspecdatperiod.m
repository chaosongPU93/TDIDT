%Currently takes 6 days of data, 2 for each ETS event, at times when 002 is
%active, and 4 days of noise, 2, 3, or 4 weeks prior (not 2003; no data
%currently for that).  Does data first, to get amplitude scaling for the
%stations, but just stores that and then gets and plots noise.  This
%version uses periodogram (no mean or taper removel yet).  18-s windows; half-overlapping.
format short e
clear all
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite; 

% jdays= [2003 062];
%         %2004 197]; 
jdays= [2003 062; %jdays are for data.
        2003 063;
        2004 196;
        2004 197;
        2005 254;
        2005 255];
noidays= [2004 182; %2 weeks 196
          2004 183; %2 weeks 197
          2005 233; %3 weeks 254
          2005 227];%4 weeks 255 
% % noidays= [2005 219; %Sun (254) %noidays are for noise
% %           2005 226; %Sun (254)
% %           2005 233; %Sun (254)
% %           2005 220; %Mon (255-weekday)
% %           2005 214; %Tue (255-weekday)
% %           2005 215; %Wed (255-weekday)
% %           2005 216; %Thu (255-weekday)
% %           2005 217];%Fri (255-weekday) %Could add 220 Mon 221 Tue
PERMROTS=[0 90 11 0;  %PGC , Yajun's "0" changed to 90.  Don't think that matters, if no splitting.
          0 90 54 9]; %LZB
POLROTS=[6 90 41 87;  %SSIB 
         0 90 41 20;  %SILB
         0 90  7 -4;  %KLNB
         4 70 48 -26; %MGCB
         4 75 38 -5]; %TWKB
%0 0 11 0   %PGC, SSIB, SILB From Chao
%6 90 41 87 
%0 0 41 20 
%From Yajun:
%     PERMROTS=[0 90 32 0;  %PGC , Yajun's "0" changed to 90.  Don't think that matters, if no splitting.
%               0 90 54 9]; %LZB
%     POLROTS=[6 85 33 86;  %SSIB from Yajun
%              0 90 39 20;  %SILB
%              0 90  7 -4;  %KLNB
%              4 70 48 -26; %MGCB
%              4 75 38 -5]; %TWKB
    stas=['PGC  '
         'SSIB '
         'SILB '];  

PERMSTA=['PGC '
         'LZB '];
POLSTA=['SSIB '
        'SILB '
        'KLNB '
        'MGCB '
        'TWKB '];
  
nsta=size(stas,1);
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;
winsec=18; %window length in seconds (not used if not to set sps)
swin=720; %window length in samples
sps=40; %round(swin/winsec)

hi=19;
lo=.1; 
npo=2;
npa=2;

ip=4;
while 2^ip < swin
   ip=ip+1;
   npts=2^ip;
end

%%%%%%%%DATA%%%%%%%%%%%
itotwin=0; 
optdat=zeros(86400*40,nsta+1);
optsnrdat=optdat;
ortdat=optdat;
%cycle over each day:
for nd=1:length(jdays(:,1))
    year=jdays(nd,1);
    YEAR=int2str(year);
    jday=jdays(nd,2);
    if jday <= 9
        JDAY=['00',int2str(jday)];
    elseif jday<= 99
        JDAY=['0',int2str(jday)];
    else
        JDAY=int2str(jday);
    end
    MO=day2month(jday,year);
    direc=[YEAR,'/',MO,'/'];
    prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];
    %Get timsSTA from the permanent stations (last one over-writes):
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
            [opt,ort,timsSTA]=readpermspec(prename,PERMSTA,PERMROTS,idx,sps,lo,hi,npo,npa,fact);
            [optsnr,~,~]=readpermspec(prename,PERMSTA,PERMROTS,idx,sps,2.5,5,npo,npa,fact);
        end
        [LIA,idx]=ismember(stas(ista,:),POLSTA,'rows');
        if LIA
            found=found+LIA; %better be 1
            if year==2003 && jday<213
                fact=20.0e-3;
            else
                fact=4.0e-3; 
            end
            [opt,ort]=readpolspec(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact);
            [optsnr,~]=readpolspec(prename,POLSTA,POLROTS,idx,sps,2.5,5,npo,npa,fact);
        end
        found=found
        STAopt(:,ista)=opt;
        STAoptsnr(:,ista)=optsnr;
        STAort(:,ista)=ort;
    end %finished one station
    %Have the day-long data in hand; record the windows.
    wins=getspecwins(jdays(nd,:));
    for n=1:size(wins,1)
        winstart=wins(n,1)*sps;
        winend=wins(n,2)*sps;
        nswins=floor((winend-winstart)/(swin/2))-1;
        for iswin=1:nswins
            itotwin=itotwin+1;
            istart=winstart+(iswin-1)*(swin/2);
            iend=istart-1+swin;
            optdat((itotwin-1)*swin+1:itotwin*swin,1)=timsSTA(istart:iend); %The first column of all of these are the time; prob. unnecessary.
            optdat((itotwin-1)*swin+1:itotwin*swin,2:4)=STAopt(istart:iend,:);
         optsnrdat((itotwin-1)*swin+1:itotwin*swin,1)=timsSTA(istart:iend);
         optsnrdat((itotwin-1)*swin+1:itotwin*swin,2:4)=STAoptsnr(istart:iend,:);
            ortdat((itotwin-1)*swin+1:itotwin*swin,1)=timsSTA(istart:iend);
            ortdat((itotwin-1)*swin+1:itotwin*swin,2:4)=STAort(istart:iend,:);
        end %finished one window
    end %finished one day
end %finished the days
optdat(itotwin*swin+1:end,:)=[]; %All the optimal windows
ortdat(itotwin*swin+1:end,:)=[]; %All the orthogonal windows
optsnrdat(itotwin*swin+1:end,:)=[]; %All the optimal windows in the good signal-to-noise range
optsnrenv=abs(hilbert(optsnrdat(:,2:4))); %The envelope of those good snr windows.  env lacks the first column of times.
ampsta=median(optsnrenv) %scale the optimal and orthogal windows by the entire record at each station
%maopt=mean(abs(optdat))
%maort=mean(abs(ortdat))
for ista=1:nsta
    optdat(:,ista+1)=optdat(:,ista+1)/ampsta(ista);
    ortdat(:,ista+1)=ortdat(:,ista+1)/ampsta(ista);
end
%maopt=mean(abs(optdat))
%maort=mean(abs(ortdat))
fid = fopen('specdata','w');
fprintf(fid,'%11.3f %10.3e %9.3e %9.3e\n',optdat');
fclose(fid);

nseg=itotwin;
%optxx=zeros(nsta*nseg,npts/2+2);
%ortxx=zeros(nsta*nseg,npts/2+2);
optpmtm=zeros(nsta,npts/2+2);
ortpmtm=zeros(nsta,npts/2+2);
optpmtmave=zeros(nseg,npts/2+2);
ortpmtmave=zeros(nseg,npts/2+2);
optxx=zeros(nsta,npts/2+2);
ortxx=zeros(nsta,npts/2+2);
optxxave=zeros(nseg,npts/2+2);
ortxxave=zeros(nseg,npts/2+2);
for iseg=1:nseg
    amps=mean(optsnrenv((iseg-1)*swin+1:iseg*swin,:)); %optsnrenv has no first column of times
    dtopt=detrend(optdat((iseg-1)*swin+1:iseg*swin,2:4)); %detrendeds have no first column of times
    dtort=detrend(ortdat((iseg-1)*swin+1:iseg*swin,2:4));
    for ista=1:nsta
        [optpmtm(ista,1:npts/2+1) f1] = pmtm(dtopt(:,ista),[],[],sps);
         optpmtm(ista,1:npts/2+1)=sqrt(optpmtm(ista,1:npts/2+1));
        [ortpmtm(ista,1:npts/2+1) ~] = pmtm(dtort(:,ista),[],[],sps);
         ortpmtm(ista,1:npts/2+1)=sqrt(ortpmtm(ista,1:npts/2+1));
        optpmtm(ista,npts/2+2) = amps(ista);  %Last row gets estimate of envelope amplitude (optimal orientation)
        ortpmtm(ista,npts/2+2) = amps(ista);  %Last row gets estimate of envelope amplitude (optimal orientation)
    end
    %%%%%%%%%
    taper=hamming(size(dtopt,1));
    dtopt(:,1)=taper.*dtopt(:,1);
    dtopt(:,2)=taper.*dtopt(:,2);
    dtopt(:,3)=taper.*dtopt(:,3);
    dtort(:,1)=taper.*dtort(:,1);
    dtort(:,2)=taper.*dtort(:,2);
    dtort(:,3)=taper.*dtort(:,3);
    %%%%%%%%%
    for ista=1:nsta
    %    %First the individual spectra at each station, nsta*nseg lines
    %    [optxx(nsta*(iseg-1)+ista,1:npts/2+1) f1] = periodogram(dtopt(:,ista),[],[],sps);
    %     optxx(nsta*(iseg-1)+ista,1:npts/2+1)=sqrt(optxx(nsta*(iseg-1)+ista,1:npts/2+1));
    %    [ortxx(nsta*(iseg-1)+ista,1:npts/2+1) f1] = periodogram(dtort(:,ista),[],[],sps);
    %     ortxx(nsta*(iseg-1)+ista,1:npts/2+1)=sqrt(ortxx(nsta*(iseg-1)+ista,1:npts/2+1));
    %    optxx(nsta*(iseg-1)+ista,npts/2+2) = amps(ista);  %Last row gets estimate of envelope amplitude (optimal orientation)
    %    ortxx(nsta*(iseg-1)+ista,npts/2+2) = amps(ista);  %Last row gets estimate of envelope amplitude (optimal orientation)
    %    %Now for the averages of the 3 stations
    %    optxxave(iseg,1:npts/2+1)=optxxave(iseg,1:npts/2+1)+optxx(nsta*(iseg-1)+ista,1:npts/2+1)/nsta;
        %This version doesn't write the individual stations
        [optxx(ista,1:npts/2+1) f1] = periodogram(dtopt(:,ista),[],[],sps);
         optxx(ista,1:npts/2+1)=sqrt(optxx(ista,1:npts/2+1));
        [ortxx(ista,1:npts/2+1) ~] = periodogram(dtort(:,ista),[],[],sps);
         ortxx(ista,1:npts/2+1)=sqrt(ortxx(ista,1:npts/2+1));
        optxx(ista,npts/2+2) = amps(ista);  %Last row gets estimate of envelope amplitude (optimal orientation)
        ortxx(ista,npts/2+2) = amps(ista);  %Last row gets estimate of envelope amplitude (optimal orientation)
    end
    %Now for the averages of the 3 stations
    optxxave(iseg,:)=mean(optxx); %Used to be mean. %wasmedian
    ortxxave(iseg,:)=mean(ortxx); %Used to be mean. %wasmedian
    optpmtmave(iseg,:)=mean(optpmtm); %Used to be mean. %wasmedian
    ortpmtmave(iseg,:)=mean(ortpmtm); %Used to be mean. %wasmedian
end
optxxsor=sortrows(optxxave,npts/2+2);
ortxxsor=sortrows(ortxxave,npts/2+2);
optpmtmsor=sortrows(optpmtmave,npts/2+2);
ortpmtmsor=sortrows(ortpmtmave,npts/2+2);
nbins=14 %This includes the half-sized bins at the high and low ends.
inbin=round(nseg/(nbins-1)); %-1 to account for the half-sized bins at the ends.
ends=round(inbin/2);
divs=[1,ends:inbin:nseg,nseg];  % division regime of bins, 1st is half the size of bins in the middle, and last size is variable
ndivs=length(divs)-1 %This should be nbins
optxxfile=zeros(nbins,npts/2+1);
optortfile=zeros(nbins,npts/2+1);
optpmtmfile=zeros(nbins,npts/2+1);
optortpmtmfile=zeros(nbins,npts/2+1);
for i=1:nbins
    optxxfile(i,:)=mean(optxxsor(divs(i):divs(i+1)-1,1:npts/2+1)); %Used to be mean. Note here there is a tiny flaw, bins are a little overalpping
    optortfile(i,:)=mean(optxxsor(divs(i):divs(i+1)-1,1:npts/2+1))./mean(ortxxsor(divs(i):divs(i+1),1:npts/2+1)); %Used to be mean.
    optpmtmfile(i,:)=mean(optpmtmsor(divs(i):divs(i+1)-1,1:npts/2+1)); %Used to be mean.
    optortpmtmfile(i,:)=mean(optpmtmsor(divs(i):divs(i+1)-1,1:npts/2+1))./mean(ortpmtmsor(divs(i):divs(i+1),1:npts/2+1)); %Used to be mean.
    %loglog(f1,mean(optxxsor(divs(i):divs(i+1),1:npts/2+1)))
    %hold on
end

%%%%%%NOISE%%%%%%%%
itotnoi=0; 
optnoi=zeros(86400*40,nsta+1); %Don't know how big this will be; one full day is generous
optsnrnoi=optnoi;
ortnoi=optnoi;
%cycle over each day:
for nd=1:length(noidays(:,1))
    year=noidays(nd,1);
    YEAR=int2str(year);
    jday=noidays(nd,2);
    if jday <= 9
        JDAY=['00',int2str(jday)];
    elseif jday<= 99
        JDAY=['0',int2str(jday)];
    else
        JDAY=int2str(jday);
    end
    MO=day2month(jday,year);
    direc=[YEAR,'/',MO,'/'];
    prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];
    %Get timsSTA from the permanent stations (last one over-writes):
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
            [opt,ort,timsSTA]=readpermspec(prename,PERMSTA,PERMROTS,idx,sps,lo,hi,npo,npa,fact);
            [optsnr,~,~]=readpermspec(prename,PERMSTA,PERMROTS,idx,sps,2.5,5,npo,npa,fact);
        end
        [LIA,idx]=ismember(stas(ista,:),POLSTA,'rows');
        if LIA
            found=found+LIA; %better be 1
            if year==2003 && jday<213
                fact=20.0e-3;
            else
                fact=4.0e-3; 
            end
            [opt,ort]=readpolspec(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact);
            [optsnr,~]=readpolspec(prename,POLSTA,POLROTS,idx,sps,2.5,5,npo,npa,fact);
        end
        found=found
        STAopt(:,ista)=opt;
        STAoptsnr(:,ista)=optsnr;
        STAort(:,ista)=ort;
    end %finished one station
    %Have the day-long data in hand; record the windows.
    wins=getspecwins(noidays(nd,:));
    for n=1:size(wins,1)
        winstart=wins(n,1)*sps;
        winend=wins(n,2)*sps;
        nswins=floor((winend-winstart)/(swin/2))-1;
        for iswin=1:nswins
            itotnoi=itotnoi+1;
            istart=winstart+(iswin-1)*(swin/2);
            iend=istart-1+swin;
            optnoi((itotnoi-1)*swin+1:itotnoi*swin,1)=timsSTA(istart:iend); %The first column of all of these are the time; prob. unnecessary.
            optnoi((itotnoi-1)*swin+1:itotnoi*swin,2:4)=STAopt(istart:iend,:);
         optsnrnoi((itotnoi-1)*swin+1:itotnoi*swin,1)=timsSTA(istart:iend);
         optsnrnoi((itotnoi-1)*swin+1:itotnoi*swin,2:4)=STAoptsnr(istart:iend,:);
            ortnoi((itotnoi-1)*swin+1:itotnoi*swin,1)=timsSTA(istart:iend);
            ortnoi((itotnoi-1)*swin+1:itotnoi*swin,2:4)=STAort(istart:iend,:);
        end %finished one window
    end %finished one day
end %finished the days
optnoi(itotnoi*swin+1:end,:)=[]; %All the optimal windows
ortnoi(itotnoi*swin+1:end,:)=[]; %All the orthogonal windows
optsnrnoi(itotnoi*swin+1:end,:)=[]; %All the optimal windows in the good signal-to-noise range
optsnrnoienv=abs(hilbert(optsnrnoi(:,2:4))); %The envelope of those good snr windows.  env lacks the first column of times.
%ampsta=median(optsnrnoienv) %scale the optimal and orthogal windows by the entire record at each station
for ista=1:nsta
    optnoi(:,ista+1)=optnoi(:,ista+1)/ampsta(ista); %These ampsta come from the data, not the noise
    ortnoi(:,ista+1)=ortnoi(:,ista+1)/ampsta(ista);
end

nsegnoi=itotnoi;
optxxnoi=zeros(nsta,npts/2+2);
ortxxnoi=zeros(nsta,npts/2+2);
optxxnoiave=zeros(nsegnoi,npts/2+2);
ortxxnoiave=zeros(nsegnoi,npts/2+2);
optpmtmnoi=zeros(nsta,npts/2+2);
ortpmtmnoi=zeros(nsta,npts/2+2);
optpmtmnoiave=zeros(nsegnoi,npts/2+2);
ortpmtmnoiave=zeros(nsegnoi,npts/2+2);
for iseg=1:nsegnoi
    amps=mean(optsnrnoienv((iseg-1)*swin+1:iseg*swin,:)); %optsnrnoienv has no first column of times
    dtopt=detrend(optnoi((iseg-1)*swin+1:iseg*swin,2:4)); %detrendeds have no first column of times
    dtort=detrend(ortnoi((iseg-1)*swin+1:iseg*swin,2:4));
    for ista=1:nsta
        [optpmtmnoi(ista,1:npts/2+1) f1] = pmtm(dtopt(:,ista),[],[],sps);
         optpmtmnoi(ista,1:npts/2+1)=sqrt(optpmtmnoi(ista,1:npts/2+1));
        [ortpmtmnoi(ista,1:npts/2+1) ~] = pmtm(dtort(:,ista),[],[],sps);
         ortpmtmnoi(ista,1:npts/2+1)=sqrt(ortpmtmnoi(ista,1:npts/2+1));
        optpmtmnoi(ista,npts/2+2) = amps(ista);  %Last row gets estimate of envelope amplitude (optimal orientation)
        ortpmtmnoi(ista,npts/2+2) = amps(ista);  %Last row gets estimate of envelope amplitude (optimal orientation)
    end
    %%%%%%%%%
    taper=hamming(size(dtopt,1));
    dtopt(:,1)=taper.*dtopt(:,1);
    dtopt(:,2)=taper.*dtopt(:,2);
    dtopt(:,3)=taper.*dtopt(:,3);
    dtort(:,1)=taper.*dtort(:,1);
    dtort(:,2)=taper.*dtort(:,2);
    dtort(:,3)=taper.*dtort(:,3);
    %%%%%%%%%
    for ista=1:nsta
        %This version doesn't write the individual stations
        [optxxnoi(ista,1:npts/2+1) f1] = periodogram(dtopt(:,ista),[],[],sps);
         optxxnoi(ista,1:npts/2+1)=sqrt(optxxnoi(ista,1:npts/2+1));
        [ortxxnoi(ista,1:npts/2+1) ~] = periodogram(dtort(:,ista),[],[],sps);
         ortxxnoi(ista,1:npts/2+1)=sqrt(ortxxnoi(ista,1:npts/2+1));
        optxxnoi(ista,npts/2+2) = amps(ista);  %Last row gets estimate of envelope amplitude (optimal orientation)
        ortxxnoi(ista,npts/2+2) = amps(ista);  %Last row gets estimate of envelope amplitude (optimal orientation)
    end
    %Now for the averages of the 3 stations
    optxxnoiave(iseg,:)=mean(optxxnoi); %Used to be mean. %was median
    ortxxnoiave(iseg,:)=mean(ortxxnoi); %Used to be mean. %was median
    optpmtmnoiave(iseg,:)=mean(optpmtmnoi); %Used to be mean. %was median
    ortpmtmnoiave(iseg,:)=mean(ortpmtmnoi); %Used to be mean. %was median
end
medianoptxxnoiave=median(optxxnoiave); %This used to "correct" tremor spectra %was median %Stick w/median of all noise windows.
medianoptpmtmnoiave=median(optpmtmnoiave); %This used to "correct" tremor spectra %was median %Stick w/median of all noise windows.
%Binning of noise follows
optxxnoisor=sortrows(optxxnoiave,npts/2+2);
ortxxnoisor=sortrows(ortxxnoiave,npts/2+2); %This sorts the orthogonal component independently.
optpmtmnoisor=sortrows(optpmtmnoiave,npts/2+2);
ortpmtmnoisor=sortrows(ortpmtmnoiave,npts/2+2); %This sorts the orthogonal component independently.
nbins=14 %This includes the half-sized bins at the high and low ends.
inbin=round(nsegnoi/(nbins-1)); %-1 to account for the half-sized bins at the ends.
ends=round(inbin/2);
divs=[1,ends:inbin:nsegnoi,nsegnoi];    % ?? Chao, 2021/09/10
ndivs=length(divs)-1 %This should be nbins
optxxnoifile=zeros(nbins,npts/2+1);
optortnoifile=zeros(nbins,npts/2+1);
optpmtmnoifile=zeros(nbins,npts/2+1);
optortpmtmnoifile=zeros(nbins,npts/2+1);
for i=1:nbins
    optxxnoifile(i,:)=mean(optxxnoisor(divs(i):divs(i+1),1:npts/2+1)); %Used to be mean.
    optortnoifile(i,:)=mean(optxxnoisor(divs(i):divs(i+1),1:npts/2+1))./mean(ortxxnoisor(divs(i):divs(i+1),1:npts/2+1)); %Used to be mean.
    optpmtmnoifile(i,:)=mean(optpmtmnoisor(divs(i):divs(i+1),1:npts/2+1)); %Used to be mean.
    optortpmtmnoifile(i,:)=mean(optpmtmnoisor(divs(i):divs(i+1),1:npts/2+1))./mean(ortpmtmnoisor(divs(i):divs(i+1),1:npts/2+1)); %Used to be mean.
    %loglog(f1,mean(optxxnoisor(divs(i):divs(i+1),1:npts/2+1)))
    %hold on
end
%Binning of noise is above
%%%NOISE%%%
h=figure('Position',[wid 1 1.5*wid hite]);
subplot(2,3,1)
loglog(f1,optxxnoifile,'color',[0.7 0.7 0.7])
hold on
% subplot(2,3,3)
% %loglog(f1,optortnoifile,'color',[0.7 0.7 0.7])
% semilogx(f1,optortpmtmnoifile,'color',[0.7 0.7 0.7])
% hold on
%%%DATA%%%
subplot(2,3,1)
loglog(f1,optxxfile)
hold on
%loglog(f1,optxxfile(end,:),'ko','markersize',2,'markerfacecolor','k')
%loglog(f1,optxxfile(end,:),'k.')
loglog(f1(1:26),optxxfile(end,1:26),'k.') %26 for swin=720
loglog(f1(1:26),optxxfile(1,1:26),'b.')
ylim([0.003 30])
xlim([0.1 20])
minfreq=2.4; maxfreq=4.; %minfreq=2.;
dum=f1-minfreq;
[~,f1min]=min(abs(dum));
dum=f1-maxfreq;
[~,f1max]=min(abs(dum));
normlzer=mean(optxxfile(:,f1min:f1max),2);
normopt=zeros(ndivs,npts/2+1);
for irow = 1:ndivs
    normopt(irow,:)=optxxfile(irow,:)/normlzer(irow);
end
%hold on
loglog(f1,0.018*normopt)    % 0.018?? Chao, 2021/09/10
loglog(f1(1:26),0.018*normopt(end,(1:26)),'k.')
loglog(f1(1:26),0.018*normopt(1,(1:26)),'b.')
minus1x=[3.3 6.6];
minus1y=[0.009 0.0045];
loglog(0.95*minus1x,3.0*minus1y,'k-','linewidth',2)  % reference line with slope of 1
%loglog(1.1*minus1x,3*minus1y,'k-','linewidth',2)
%%%%%%%%doesn't work (axis equal workaround)
%xLimits = [0.003 30];                    %# Limits for the x axis
%yLimits = [0.1 20];                      %# Limits for the y axis
%logScale = diff(yLimits)/diff(xLimits);  %# Scale between the x and y ranges
%powerScale = diff(log10(yLimits))/diff(log10(xLimits));
%set(gca,'Xlim',xLimits,'YLim',yLimits,'DataAspectRatio',[1 logScale/powerScale 1]);  %#   data aspect ratio
%%%%%%%%
set(gca,'XTick',[0.1 1 10]')
title([int2str(jdays(1,2)),'-',int2str(jdays(end,2)),' ',num2str(swin/sps),'-sec windows'])
subplot(2,3,3)
% loglog(f1,optortfile)
% loglog(f1(1:26),optortfile(end,1:26),'k.')
% loglog(f1(1:26),optortfile(1,1:26),'b.')
semilogx(f1,optortpmtmfile)
hold on
semilogx(f1(1:26),optortpmtmfile(end,1:26),'k.')
semilogx(f1(1:26),optortpmtmfile(1,1:26),'b.')
semilogx(f1,optortpmtmnoifile,'color',[0.7 0.7 0.7])
ylim([0 4])
xlim([0.1 20])
set(gca,'XTick',[0.1 1 10]')

tempdat=load('ChaoPMTM002_0.1-15Hz_2pass_100sps_0.1-15Hz_CC2-8Hz_le14sh');  %This from CNDC/comp.m
%tempdat=load('0.1-15Hz_2pass_100sp_0.1-15Hz_le14shSRR_pmtm');  %This from BOSTOCK/NEW/WFPV/tempspec.m
%tempdat=load('1-15Hz_2pass_100sps_0.1-15Hz_le14sh_pmtm'); %This from BOSTOCK/NEW/WFPV/tempspec.m
ftemp=tempdat(:,1)';
tempspec=tempdat(:,2:4)';
%noisespec=tempdat(:,3:end)';
%pmtmdat=load('62_6days_18sec_14bins.pmtm'); %This from getspecdat.m (pmtm version)
%fpmtm=pmtmdat(:,1)';
%pmtmspec=pmtmdat(:,12:end)'; %This does the larger-amplitude half, I think.
subplot(2,3,2)
%loglog(ftemp,28*noisespec,'color',[0.7 0.7 0.7])
minfreq=2.4; maxfreq=4.; %minfreq=2.;
dum=f1-minfreq;
[~,f1min]=min(abs(dum));
dum=f1-maxfreq;
[~,f1max]=min(abs(dum));
normlzer=mean(optpmtmfile(:,f1min:f1max),2);
normoptpmtm=zeros(ndivs,npts/2+1);
for irow = 1:ndivs
    normoptpmtm(irow,:)=optpmtmfile(irow,:)/normlzer(irow);
end
loglog(ftemp,33*mean(tempspec),'m') %was median, to shift the temp spectrum vertically
hold on
loglog(ftemp(1:21),33*mean(tempspec(:,1:21)),'m.') %was median
%loglog(ftemp,28*(tempspec-median(noisespec)),'c')
loglog(f1,normoptpmtm(end-3:end,:))
loglog(1.05*minus1x,140*minus1y,'k-','linewidth',2)
xlim([0.1 20])
ylim([0.003 30])
set(gca,'XTick',[0.1 1 10]')

subplot(2,3,2)
optpmtmrefined=optpmtmfile;
for i=1:size(optpmtmrefined,1)
    optpmtmrefined(i,:)=optpmtmrefined(i,:)-medianoptpmtmnoiave(1:end-1);
end
% loglog(f1,optxxrefined,'color',[0.7 0.7 0.7])
% hold on
normlzer=mean(optpmtmrefined(:,f1min:f1max),2); %I think this is redundant (from above) but not harmful
normoptref=zeros(ndivs,npts/2+1);
for irow = 1:ndivs
    normoptref(irow,:)=optpmtmrefined(irow,:)/normlzer(irow);
end
% %hold on
minfreq=0.6; maxfreq=1.8; 
dum=f1-minfreq;
[~,f1min]=min(abs(dum));
dum=f1-maxfreq;
[~,f1max]=min(abs(dum));
loglog(f1(f1min:f1max),normoptref(end,f1min:f1max),'k')
% loglog(f1(1:26),0.01*normoptref(end,(1:26)),'k.')
% loglog(f1(1:26),0.01*normoptref(1,(1:26)),'b.')
% loglog(ftemp,0.28*median(tempspec),'c')
% ylim([0.0003 3])
% xlim([0.1 20])
% set(gca,'XTick',[0.1 1 10]')

set(h,'PaperPosition',[0.25 0.25 8 10.5])
print('-depsc',[int2str(jdays(1,2)),'_',int2str(size(jdays,1)),'days_',num2str(swin/sps),'sec','_',int2str(nbins),'bins_period.eps'])
