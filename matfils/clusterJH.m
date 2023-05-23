%plots maps only, on one page.
 clear all
 close all
scrsz=get(0,'ScreenSize');
wid=scrsz(3);
hite=scrsz(4);
scrat=wid/hite; 
set(0,'DefaultFigureVisible','on');

mapinputs=['mapSTAS.UN.100sps.mdist6_else0.5_N10e3.05.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.mdist6_else0.5_N10e3.65.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.mdist6_else0.5_N10e3.95.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.mdist6_else0.5_N10e4.56.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.mdist6_else0.5_N10e5.56.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.mdist6_else0.5_N10e6.56.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5'];
mapinputs=['mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat100.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5'];
mapinputs=['mapSTAS.UN.100sps.ell_3-1.5.else0nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';           
           'mapSTAS.UN.100sps.ell_3-1.5.else0nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.ell_3-1.5.else0nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.ell_3-1.5.else0nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.ell_3-1.5.else0nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.ell_3-1.5.else0nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5'];
           %'mapSTAS.UN.100sps.ell_3-1.5.else0nsat100.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5'];
mapinputs=['mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';           
           'mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
           'mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5'];
           %'mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat100.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5'];
% mapinputs=['mapSTAS.UN.100sps.ell_2.75-1.else0nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%            'mapSTAS.UN.100sps.ell_2.75-1.else0nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%            'mapSTAS.UN.100sps.ell_2.75-1.else0nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%            'mapSTAS.UN.100sps.ell_2.75-1.else0nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%            'mapSTAS.UN.100sps.ell_2.75-1.else0nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%            'mapSTAS.UN.100sps.ell_2.75-1.else0nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5'];
% mapinputs=['mapSTAS.UN.100sps.ell_3-1.25.else0nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';           
%            'mapSTAS.UN.100sps.ell_3-1.25.else0nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%            'mapSTAS.UN.100sps.ell_3-1.25.else0nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%            'mapSTAS.UN.100sps.ell_3-1.25.else0nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%            'mapSTAS.UN.100sps.ell_3-1.25.else0nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%            'mapSTAS.UN.100sps.ell_3-1.25.else0nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5'];
tims=[300 500;
      700 900;
      100 1700];
nwin=size(tims,1);
winlen=tims(:,2)-tims(:,1);
ALL128='~arubin/data2/arubin/CNDC/ARMMAP/MAPS/all_3-day.86.20.80.115.50_2-6-ms19-128s.dx.sorted';
origAll=load(ALL128);
%for a rectangle oriented at 45˚
boxx=[-5 2.25; -5 3.5; 1 3.5; 1 2.25; -5 2.25];
angrot=45*pi/180;
boxtorot=complex(boxx(:,1),boxx(:,2));
boxrot=boxtorot*exp(1i*angrot);
%for an ellipse oriented at 45˚
xaxis=3.;
yaxis=1.5;
xaxis=2.75;
yaxis=1.0;
xaxis=3.;
yaxis=1.25;
xxp=-xaxis:0.01:xaxis;
yyp=sqrt(yaxis.^2*(1 -(xxp/xaxis).^2 ));
%xxp=[xxp,xxp];
%yyp=[yyp,-yyp];
locp=complex(xxp, yyp);
locn=complex(xxp,-yyp);
angrot=45*pi/180;
locprot=locp*exp(1i*angrot);
locnrot=locn*exp(1i*angrot);
shiftor=complex(0,-0.5); %(in km)
ellipsep=locprot+shiftor;
ellipsen=locnrot+shiftor;
%
nfile=size(mapinputs,1);
%nclust=11; %make this odd
%target=9;
%prctarget=100*(target-0.5)/nclust;
%offset=(nclust-1)/2;
maxlag=1800; %in seconds
h=figure('Units','inches','Position',[1 1 11 8.5]);
%h=figure('Position',[wid/20 1 4.5*wid/5 hite]); %center
nevs=zeros(1,nfile);
for ifile = 1:nfile
    mapfil=load(mapinputs(ifile,:));
    nevs(ifile)=size(mapfil,1);
end
maxnev=max(nevs);
dxp=zeros(maxnev*(maxnev-1)/2,nfile)-999;
dyp=dxp;
for ifile = 1:nfile
    mapfil=load(mapinputs(ifile,:));
    nsat=strfind(mapinputs(ifile,:),'nsat');
    elsew=strfind(mapinputs(ifile,:),'else');
    spsloc=strfind(mapinputs(ifile,:),'sps');
    %nevs(ifile)=size(mapfil,1);
    extrac1=mapfil(1:nevs(ifile)-1,:);
    extrac2=mapfil(2:nevs(ifile),:);
    dists=sqrt((extrac2(:,15)-extrac1(:,15)).^2+(extrac2(:,14)-extrac1(:,14)).^2);
    %h=figure('Position',[wid/5 1 3*wid/5 hite]); %center
    %subplot(2,1,1)
    %plot(extrac1(:,1),dists,'r+')
    %%plot(extrac1(:,1),dists2,'bo')
    %ymax=5;
    %ylim([0 ymax])
    %text(0.9*extrac1(end,1), 0.85*ymax, num2str(median(dists)));
    %title(mapinputs(ifile,4:11))
    %drawnow
    xses=mapfil(:,15);
    yses=mapfil(:,14);
    orig=mapfil(:,1);
%     times=zeros(1,nevs-nclust+1);
%     xymeds=zeros(1,nevs-nclust+1);
%     xytarget=zeros(1,nevs-nclust+1);
% needed (below)?
    rotang=(-50./180.)*pi; %Not used unless you want to check "along strike" and "along dip"
    %rotline=(pi/2+rotang);
    complexloc=xses+1i*yses;
    rotlocs=complexloc*exp(1i*rotang);
    xprime=real(rotlocs);
    yprime=imag(rotlocs);
%    xmed=zeros(nevs-nclust+1,1);
%    ymed=zeros(nevs-nclust+1,1);
%    nsmed=zeros(nevs-nclust+1,1);
%    ccav=zeros(nevs-nclust+1,1);
%    cc12=zeros(nevs-nclust+1,1);
%    cc13=zeros(nevs-nclust+1,1);
%    cc32=zeros(nevs-nclust+1,1);
% needed (above)?
    entry=0;
    for i = 1:nevs(ifile)-1
        step=1;
        while (orig(i+step)-orig(i) <= maxlag) && (i+step<nevs(ifile))
	        if orig(i+step,1)-orig(i,1)>=4
                entry=entry+1;
                dxp(entry,ifile)=xprime(i+step)-xprime(i);
                dyp(entry,ifile)=yprime(i+step)-yprime(i);
            end
            step=step+1; 
        end
    end
%     clear dxp
%     clear dyp
%       xmed(i)=median(xses(i:i+nclust-1)); %median x location of detections in temporal cluster.
%       ymed(i)=median(yses(i:i+nclust-1));
%       nsmed(i)=median(mapfil(i:i+nclust-1,17));
%       xdists=xses(i:i+nclust-1)-xmed(i); %x distance of those in cluster from that median value.
%       ydists=yses(i:i+nclust-1)-ymed(i);
%       xydists=sqrt(xdists.^2+ydists.^2); %absolute distance of those in cluster from that median value.
%       xymeds(i)=median(xydists);
%       xymeds(i)=median(xydists);         %median of those distances (red plusses).
%       times(i)=mapfil(i+offset,1);
%       ccav(i)=mapfil(i+offset,4);
%       cc12(i)=mapfil(i+offset,21);
%       cc13(i)=mapfil(i+offset,22);
%       cc32(i)=mapfil(i+offset,23);
%       xytarget(i)=prctile(xydists,prctarget); %chose percentile (blue circles).
    %end
    %subplot(2,1,2)
    %plot(times,xymeds,'r+')
    %hold on
    %plot(times,xytarget,'bo')
    %ymax=1.25;
    %ylim([0 ymax])
    %text(0.9*times(end), 0.85*ymax, num2str(median(xymeds)), 'color','r');
    %text(0.8*times(end), 0.85*ymax, num2str(median(xytarget)), 'color','b');
    %title([mapinputs(ifile,nsat:nsat+6),'  ',mapinputs(ifile,elsew:elsew+7)])

    for iwin = 1:nwin
        if orig(1,1) >= tims(iwin,1)
            istart=1;
        end
        for i=2:nevs(ifile)-1
            if orig(i-1,1) < tims(iwin,1) && orig(i,1) >= tims(iwin,1)
                istart=i;
            end
            if orig(i+1,1) > tims(iwin,2) && orig(i,1) <= tims(iwin,2)
                iend=i;
            end
        end
        if orig(nevs(ifile),1) <= tims(iwin,2)
            iend=nevs(ifile);
        end
        colormap(jet)
        hold on
        %drawnow
        subplot(2,3,ifile,'align')
        plot(origAll(:,2),origAll(:,1),'o','Color',[0.6 0.6 0.6],'MarkerSize',0.3,'MarkerFaceColor',[0.6 0.6 0.6]); %'linewidth',0.2)
        scatter3(xses(istart:iend),yses(istart:iend),orig(istart:iend),13,orig(istart:iend),'filled');
        %scatter(xses(istart:iend),yses(istart:iend),25,orig(istart:iend),'filled')
%         %pause(0.5);
%         %uistack(h(1),'top');
%         %pause(0.5);
%         uistack(h(2),'top')
%         %end
        %Rotated
        axis equal
        axis([-7 4 -5 4])
        set(gca,'XTick',-7:1:4);
        clo=tims(iwin,1);
        chi=tims(iwin,2);
        caxis([clo chi])
        %colorbar
        title([mapinputs(ifile,nsat:nsat+6),'  ',mapinputs(ifile,elsew:elsew+4)],'fontweight','normal')
        box on
        %plot(real(boxrot),imag(boxrot),'k-')
        plot(real(ellipsep),imag(ellipsep),'k','linewidth',1)
        plot(real(ellipsen),imag(ellipsen),'k','linewidth',1)

		%h=figure('Position',[wid/2 1 wid/2 hite]); %center
        %subplot(3,1,1)
        %ymax=5;
        %axis([tims(iwin,1) tims(iwin,2) 0 ymax])
        %hold on
        %plot(extrac1(:,1),dists,'r+') %This plots distances between consecutive detections.
        %%plot(times,(nsmed*(ymax/2)/3)+ymax/2) %I think this puts +/3 km at top and bottom.  Medians.
        %%plot(times,ymax/2,'k','linewidth',2) 
        %%axis([tims(iwin,1)-winlen(iwin)/2 tims(iwin,2)+winlen(iwin)/2 0 ymax])
        %text(0.9*extrac1(end,1), 0.85*ymax, num2str(median(dists)));
        %title(mapinputs(ifile,4:11))

        %subplot(3,1,2)
        %axis([tims(iwin,1) tims(iwin,2) 0 ymax])
        %hold on
        %%hrf=plotreflines3rds(gca,bostsec3,'x','k','b');
        %%hrf=plotreflines3rds(gca,bostsec2,'x','c','m');
        %%hrf=plotreflines3rds(gca,bostsec1,'x','r','t');
        %plot(times,ymax*cc12,'g')
        %plot(times,ymax*cc13,'r')
        %plot(times,ymax*cc32,'b')
        %plot(times,ymax*cc32,'bo') %*ymax2/1,'go') 
        %%plot(PGC(:,1),8*PGC(:,2),'r') %PGC(:,2)/2
        %xpdiff=xprime-circshift(xprime,1);
        %ypdiff=yprime-circshift(yprime,1);
        %plot(orig,abs(xpdiff),'r+')
        %plot(orig,abs(ypdiff),'bd')
        %plot([tims(iwin,1) tims(iwin,2)],[0 0],'k')
        %title([mapinputs(ifile,nsat:nsat+6),'  ',mapinputs(ifile,elsew:elsew+7)])

        %subplot(3,1,3)
        %ymax2=1.2;
        %%plot(PGC(:,1),PGC(:,2),'r') %PGC(:,2)/2
        %hold on
        %%plot(xcav(:,1),xcav(:,2),'g')
        %plot(times,xymeds,'r+') %median of distances from meam location (red plusses).
        %plot(times,xytarget,'bo') %100*(target-0.5)/nclust percentile of distances from meam location (red plusses).
        %%plot(times,(nsmed*(ymax2/2)/3)+ymax2/2) %I think this puts +/3 km at top and bottom.  Medians.
        %%plot(times,ymax2/2,'k--') 
        %plot(times,ccav,'go') %*ymax2/1,'go') 
        %%axis([tims(iwin,1)-winlen(iwin)/2 tims(iwin,2)+winlen(iwin)/2 0 ymax2])
        %axis([tims(iwin,1) tims(iwin,2) 0 ymax2])
        %plot([tims(iwin,1) tims(iwin,2)],[0 0],'k')
    end
end

 %set(h,'PaperPosition',[0.25 0.25 8 10.5])
 %set(h,'PaperPosition',[0.25 0.25 10.5 8])
 orient(h,'landscape')
 print(h,'-depsc','-painters',['maps_',mapinputs(ifile,spsloc+4:elsew+7),'.eps'])
 %print(h,'-dpdf',['maps_',mapinputs(ifile,elsew:elsew+7),'.pdf'])

h=figure('Units','inches','Position',[1 1 11 8.5]);
%h=figure('Position',[wid/20 1 4.5*wid/5 hite]); %center
for ifile = 1:nfile
    subplot(2,3,ifile,'align')
    xcol=dxp(:,ifile);
    ycol=dyp(:,ifile);
    xcol(xcol==-999)=[];
    ycol(ycol==-999)=[];
    ymax=6e4;
    xmed=median(abs(xcol));
    ymed=median(abs(ycol));
    histogram(abs(xcol),'binwidth',0.1)
    hold on
    text(xmed,0.8*ymax,num2str(xmed))
    plot([xmed,xmed],[0,ymax],'b--')
    histogram(abs(ycol),'binwidth',0.1)
    plot([ymed,ymed],[0,ymax],'m--')
    text(ymed,0.9*ymax,num2str(ymed))
    xlim([0 6])
    box on
end
 %set(h,'PaperPosition',[0.25 0.25 8 10.5])
 orient(h,'landscape')
 print(h,'-depsc',['dlocs_',mapinputs(ifile,spsloc+4:elsew+7),'.eps'])
