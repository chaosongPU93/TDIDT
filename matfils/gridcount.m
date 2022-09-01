%
clear all
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

scrsz=get(0,'ScreenSize');
wid=scrsz(3);
hite=scrsz(4);
scrat=wid/hite; 
h=figure('Position',[0.3*wid 1 wid/2 hite]); %center

cutout='ellipse';
%cutout='rectang';
if strcmp(cutout,'ellipse')
    alim=1.5; %semi-major axis
    blim=1.0; %semi-minor axis
    alim=1.75; %semi-major axis
    blim=1.0; %semi-minor axis
    % alim=3; %semi-major axis
    % blim=3; %semi-minor axis
    xel=-alim:0.01:alim;
    yel=blim*sqrt(1-(xel/alim).^2);
    xel=[xel,fliplr(xel)];
    yel=[yel,-fliplr(yel)];
    angrot=45*pi/180;
    loc=complex(xel,yel);
    locrot=loc*exp(1i*angrot);
    xelrot=real(locrot);
    yelrot=imag(locrot);
    shiftor=[-0.1 -0.5]; %(in km)
    outline(:,1)=xelrot+shiftor(1);
    outline(:,2)=yelrot+shiftor(2);
elseif strcmp(cutout,'rectang')
    NS=[-3 3];
    EW=[-6 2];
    xses=EW(1):0.01:EW(2);
    yses=NS(1):0.01:NS(2);
    lenx=length(xses);
    leny=length(yses);
    outline(1:lenx,1)=xses;
    outline(1:lenx,2)=NS(1);
    outline(lenx+1:lenx+leny,1)=EW(2);
    outline(lenx+1:lenx+leny,2)=yses;
    outline(lenx+leny+1:2*lenx+leny,1)=fliplr(xses);
    outline(lenx+leny+1:2*lenx+leny,2)=NS(2);
    outline(2*lenx+leny+1:2*lenx+2*leny,1)=EW(1);
    outline(2*lenx+leny+1:2*lenx+2*leny,2)=fliplr(yses);
end

%     mapinputs=['map2003.062.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk0.5';
%                'map2003.063.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk0.5';
%                'map2004.196.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk0.5';
%                'map2004.197.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk0.5';
%                'map2004.198.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk0.5';
%                'map2004.199.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk0.5';
%                'map2005.254.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk0.5';
%                'map2005.255.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk0.5';
%                'map2005.256.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk0.5'];
fname= 'mapall.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk0.5';
aaa=load(fname);
%subplot(2,2,1,'align')
axes('Units','normalized','position',[0.09 0.6 0.4 0.35])
% [singletons,multiples]=gridcat(fname);
% scatter(singletons(:,1),singletons(:,2),5,log10(singletons(:,3)))
% hold on
% scatter(multiples(:,1),multiples(:,2),5,log10(multiples(:,3)),'filled')
density1d = density_pixel(aaa(:,15),aaa(:,14));
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(dumhf(:,1),dumhf(:,2),5,log10(dumhf(:,3)));
hold on
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(dumhf(:,1),dumhf(:,2),5,log10(dumhf(:,3)),'filled');
axis equal
axis([-8 3 -4 4])
SEbord=-1.35; NWbord=1.6; %SEbord=1.3; NWbord=1.5;
SE=[-10 SEbord;
  10 SEbord];
NW=[-10 NWbord;
  10 NWbord];
angrot=45*pi/180;
SEline=complex(SE(:,1),SE(:,2));
NWline=complex(NW(:,1),NW(:,2));
SEline=SEline*exp(1i*angrot);
NWline=NWline*exp(1i*angrot);
ax=gca;
plot([0 0],ax.YLim,'k--');
plot(ax.XLim,[0 0],'k--');
scatter(shiftor(1),shiftor(2),10,'ko','linew',1);
plot(SEline,'r--','linewidth',2)
plot(NWline,'r--','linewidth',2)
% angs=0:pi/100:2*pi;
% circ=1.0*exp(1i*angs);
% plot(real(circ)+circcent(1,1),imag(circ)+circcent(1,2),'k','linewidth',1)
plot(outline(:,1),outline(:,2),'k-','linewidth',2)
%plot(xelrot+shiftor(1),yelrot+shiftor(2),'k-','linewidth',2)
colormap('jet');
cbarul=colorbar('NorthOutside');
caxis([0 1.4]);
xlabel(cbarul,'Log_{10}(# tremor detections / pixel)');
set(cbarul,'XTick',0.2:0.2:1.2); 
xlabel('km East','fontsize',10);
ylabel('km North','fontsize',10)
text(-7.5,3.5,'a','fontsize', 12)
box on

%This no longer plots the proper elliptical boundary, now that shiftor has been moved to the top.
% axes('Units','normalized','position',[0.55 0.6 0.4 0.35])
% scatter(singletons(:,4),singletons(:,5),5,log10(singletons(:,3)))
% hold on
% scatter(multiples(:,4),multiples(:,5),5,log10(multiples(:,3)),'filled')
% axis equal
% axis([-8 3 -4 4])
% shiftor=[-0.35 -0.25]; %(in km)
% plot(outline(:,1),outline(:,2),'k-','linewidth',2)
% cbarul=colorbar('NorthOutside');
% xlabel(cbarul,'Log_{10}(# tremor detections / pixel)');
% set(cbarul,'XTick',0.2:0.2:1.2); 
% xlabel('km East','fontsize',10);
% ylabel('km North','fontsize',10)
% text(-7.5,3.5,'a','fontsize', 12)
% box on

axes('Units','normalized','position',[0.09 0.22 0.4 0.35])
locfile=load(fname);
excat=locfile;
if strcmp(cutout,'ellipse') 
    cshiftor=complex(shiftor(1),shiftor(2)); %Need to do this below as well
    rotshift=cshiftor*exp(-1i*angrot);
    PseudoE=(locfile(:,16)-real(rotshift))/alim;
    PseudoN=(locfile(:,17)-imag(rotshift))/blim;
    A=PseudoE.^2+PseudoN.^2>1; 
    excat(A,:)=[];
elseif strcmp(cutout,'rectang')
    xloc=locfile(:,15);
    yloc=locfile(:,14);
    A = xloc<EW(1) | xloc>EW(2) | yloc<NS(1) | yloc>NS(2);
    excat(A,:)=[];
end
plot(excat(:,15),excat(:,14),'k.')
axis equal
axis([-8 3 -4 4])
xlabel('km East','fontsize',10);
ylabel('km North','fontsize',10)
text(-7.5,3.5,'a','fontsize', 12)
box on
% if strcmp(cutout,'ellipse')
%     fid=fopen([fname,'_',num2str(alim),'-',num2str(blim)],'w');
% elseif strcmp(cutout,'rectang')
%     fid=fopen([fname,'_',num2str(EW(1)),num2str(EW(2)),num2str(NS(1)),num2str(NS(2))],'w');
% end
% fprintf(fid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %8.3f %10.3f %7.3f %7.3f %7.3f %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %5i %5i %10.3e %6.3f %6.3f %6.3f %7.3f\n',excat');
% fclose(fid);

%%
mapinputs=['map2003.062.002.loff1.5.ccmin0.42.nponpa22_1.25-6.5-ms18-4s.onepk.5';
%            'map2003.063.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk0.5';
%            'map2004.196.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk0.5';
%            'map2004.197.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk0.5';
%            'map2004.198.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk0.5';
%            'map2004.199.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk0.5';
%            'map2005.254.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk0.5';
%            'map2005.255.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk0.5';
%            'map2005.256.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk0.5'
           ];
nfile=size(mapinputs,1);
nwintot=0;
for ifile=1:nfile
    fname=mapinputs(ifile,:);
    locfile=load(fname);
    excat=locfile;
    if strcmp(cutout,'ellipse')
        cshiftor=complex(shiftor(1),shiftor(2)); 
        rotshift=cshiftor*exp(-1i*angrot);
        PseudoE=(locfile(:,16)-real(rotshift))/alim;
        PseudoN=(locfile(:,17)-imag(rotshift))/blim;
        A=PseudoE.^2+PseudoN.^2>1; 
        excat(A,:)=[];
    elseif strcmp(cutout,'rectang')
    end
    tims=excat(:,8);  % start time of strongest 0.5-s win
    tdiffminus=tims(2:end)-tims(1:end-1); % separation time between itself and its preceding 
    tdiffmin=min(tdiffminus(1:end-1),tdiffminus(2:end));  % inter-event time as min sep. between its neighbour
    timcheck=tims(2:end-1); % time of detection that actually have inter-event times
    PGSIs=excat(:,2);
    PGSSs=excat(:,3);
    SSSIs=PGSSs-PGSIs;
    mshift=18;
    halfdaysec=86400/2;
    
    %plot detections and inter-event times
    h=figure('Position',[0.3*wid 1 wid/2.9 hite]); %center
    subplot(8,1,1,'align');
    hold on
    plot(tims,SSSIs,'ks','MarkerSize',2);
    plot(tims,PGSIs,'ro','MarkerSize',2);
    %plot(tims,PGSSs,'k*','MarkerSize',2);
    axis([0 halfdaysec/2 -mshift mshift]);
    ylabel('samples')
    title(mapinputs(ifile,:))
    box on
    subplot(8,1,2,'align');
    hold on
    plot(tims,SSSIs,'ks','MarkerSize',2);
    plot(tims,PGSIs,'ro','MarkerSize',2);
    %plot(tims,PGSSs,'k*','MarkerSize',2);
    axis([halfdaysec/2 halfdaysec -mshift mshift]);
    ylabel('samples')
    box on
    subplot(8,1,3,'align');
    hold on
    plot(tims,SSSIs,'ks','MarkerSize',2);
    plot(tims,PGSIs,'ro','MarkerSize',2);
    %plot(tims,PGSSs,'k*','MarkerSize',2);
    axis([halfdaysec 3*halfdaysec/2 -mshift mshift]);
    ylabel('samples')
    box on
    subplot(8,1,4,'align');
    hold on
    plot(tims,SSSIs,'ks','MarkerSize',2);
    plot(tims,PGSIs,'ro','MarkerSize',2);
    %plot(tims,PGSSs,'k*','MarkerSize',2);
    axis([3*halfdaysec/2 2*halfdaysec -mshift mshift]);
    xlabel('sec')
    ylabel('samples')
    box on
    %
    cutoff=30; %20; %100;
    xcut=[0 86400];
    ycut=[cutoff cutoff];
    subplot(8,1,5,'align');
    semilogy(timcheck,tdiffmin,'k.');
    xlim([0 halfdaysec/2]);
    ylim([0.4 1000]);
    hold on
    plot(xcut,ycut,'r--')
    plot(locfile(:,1),cutoff+100*(10.^(locfile(:,4))-10^0.44),'g') %2.51=10^0.4
    xlabel('sec')
    ylabel('log_{10}(Dt)')
    set(gca,'YTick',[1 10 100])
    box on
    subplot(8,1,6,'align');
    semilogy(timcheck,tdiffmin,'k.');
    xlim([halfdaysec/2 halfdaysec]);
    ylim([0.4 1000]);
    hold on
    plot(xcut,ycut,'r--')
    plot(locfile(:,1),cutoff+100*(10.^(locfile(:,4))-10^0.44),'g')
    xlabel('sec')
    ylabel('log_{10}(Dt)')
    set(gca,'YTick',[1 10 100])
    box on
    subplot(8,1,7,'align');
    semilogy(timcheck,tdiffmin,'k.');
    xlim([halfdaysec 3*halfdaysec/2]);
    ylim([0.4 1000]);
    hold on
    plot(xcut,ycut,'r--')
    plot(locfile(:,1),cutoff+100*(10.^(locfile(:,4))-10^0.44),'g')
    xlabel('sec')
    ylabel('log_{10}(Dt)')
    set(gca,'YTick',[1 10 100])
    box on
    subplot(8,1,8,'align');
    semilogy(timcheck,tdiffmin,'k.');
    xlim([3*halfdaysec/2 2*halfdaysec]);
    ylim([0.4 1000]);
    hold on
    plot(xcut,ycut,'r--')
    plot(locfile(:,1),cutoff+100*(10.^(locfile(:,4))-10^0.44),'g')
    xlabel('sec')
    ylabel('log_{10}(Dt)')
    set(gca,'YTick',[1 10 100])
    box on
    drawnow
    
    keyboard
    %
    nev=size(excat,1);
    times=excat(:,1);
    nwin=0;
    j=0;
    first=1;
    %date=str2double(mapinputs(ifile,4:11)); Original strategy; turn '2003.062' etc. into integers.
    %dat=single(date);
    %bracket=[int32(date) int32(1000*(date-floor(date)))];
    if isequal(ifile,1)
        if strcmp(cutout,'ellipse') 
            fid=fopen(['timwins_',num2str(alim),'-',num2str(blim)],'w');
        elseif strcmp(cutout,'rectang')
            fid=fopen(['timwins_',num2str(EW(1)),num2str(EW(2)),num2str(NS(1)),num2str(NS(2))],'w');
        end
        fprintf(fid,['if isequal(',mapinputs(ifile,4:7),' ',mapinputs(ifile,9:11),')']);
    else
        if strcmp(cutout,'ellipse') 
            fid=fopen(['timwins_',num2str(alim),'-',num2str(blim)],'a');
        elseif strcmp(cutout,'rectang')
            fid=fopen(['timwins_',num2str(EW(1)),num2str(EW(2)),num2str(NS(1)),num2str(NS(2))],'a');
        end
        fprintf(fid,['\n','elseif isequal(''',mapinputs(ifile,4:7),'''',mapinputs(ifile,9:11),')']);
    end
    fprintf(fid,['\n','mapinputs=[''',mapinputs(ifile,4:11),''',suffix];']);
    while j < nev
        j=j+1;
        i=j+1;
        while i < nev && times(i)-times(i-1) < cutoff
            i=i+1; %consecutive events are close enough to check the next
        end
        i=i-1; %the last one didn't make it
        nin=i-j+1;
        if nin > 5
            nwin=nwin+1;
            %bracket(nwin+1,:)=[int32(times(j)) int32(times(i))];
            if isequal(first,1)
                fprintf(fid,['\n','    wins=[',int2str(times(j)),' ',int2str(times(i)),';']);
                first=0;
            else
                fprintf(fid,['\n','          ',int2str(times(j)),' ',int2str(times(i)),';']);
            end
        end
        j=i; %this i didn't make it.  Make j the next one, the next time through the "while" loop.
    end
    %bracketsave(nwintot+1:nwintot+nwin+1,:)=bracket; 
    %nwintot=nwintot+(nwin+1); % +1 b.c. of date line
    %clear bracket
    fprintf(fid,'\n');
    fclose(fid);
    
    if strcmp(cutout,'ellipse')
        fid=fopen([fname,'_',num2str(alim),'-',num2str(blim)],'w');
    elseif strcmp(cutout,'rectang')
        fid=fopen([fname,'_',num2str(EW(1)),num2str(EW(2)),num2str(NS(1)),num2str(NS(2))],'w');
    end
    fprintf(fid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %8.3f %10.3f %7.3f %7.3f %7.3f %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %5i %5i %10.3e %6.3f %6.3f %6.3f %7.3f\n',excat');
    fclose(fid);
    %fprintf(fid,'%8i %8i\n',bracketsave');
end
