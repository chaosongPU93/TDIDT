% test for plotting only one time winow in chaoSTA

close all
figure

%% subplot 1
subplot(3, 1, 1,'align');
nt=1;
%%% STA1file(*, 2) stores the seismogram of each window
yma=max(max([STA1file(winlen*(nt-1)+1:winlen*nt,2) STA2file(winlen*(nt-1)+1:winlen*nt,2) ...
    STA3file(winlen*(nt-1)+1:winlen*nt,2)]));
ymi=min(min([STA1file(winlen*(nt-1)+1:winlen*nt,2) STA2file(winlen*(nt-1)+1:winlen*nt,2) ...
    STA3file(winlen*(nt-1)+1:winlen*nt,2)]));
yma=2.4*max(yma,-ymi);
% % Lines below plot + or - template; unstretched:
%                         if match(nt,1) >= match(nt,2) %if positive cc with template is larger ...
%                             plot(timstemp((nt-1)*templen+1:nt*templen),STAtemps(whichtoplot,:)*yma/2.4,'c','linewidth',2)
%                         else %if cc with negative template is larger ...
%                             plot(timstempneg((nt-1)*templen+1:nt*templen),-STAtemps(whichtoplot,:)*yma/2.4,'m','linewidth',2)
%                         end
% % Lines above plot + or - template; unstretched:
% Lines below plot + stretched template:
%                         tempplot=squeeze(STAstr(whichtoplot,imaxstretch(whichtoplot,nt),:))*scaletemp(whichtoplot,nt);
%                         plot(timstempstr(whichtoplot,(nt-1)*stretchlen+1:nt*stretchlen),tempplot,'c','linewidth',1)
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
xvect=[is is+2*(yma-ymi)*(winlen/160.)]; %*mean(scalefact)]; %amplitude bar originally scaled for 4-s window.  Not sure "mean" is necessary.
yvect=[-0.9*yma -0.9*yma];
plot(xvect,yvect,'r','linewidth',3)
plot([is+1/hi is+1/lo],[-0.8*yma -0.8*yma],'k','linewidth',3)
xcvect=[is+STA32file(nt,3)/sps is+(STA32file(nt,3)+cncntr)/sps];
ycvect=[0.95*yma 0.95*yma];
plot(xcvect,ycvect,'b','linewidth',3)
plot(bostsec,0.93*yma,'ro','MarkerSize',4,'MarkerFaceColor','r')
text(is+STA32file(nt,3)/sps,  0.58*yma, num2str(match(nt,1),2),'fontsize',6,'color','b');
text(is+STA32file(nt,3)/sps, -0.75*yma, num2str(match(nt,2),2),'fontsize',6,'color','r');
%                         text(is+0.2, 0.66*yma, int2str(STA13file(nt,1)),'fontsize',6);
%                         text(ien-0.6, 0.66*yma, int2str(STA12file(nt,1)),'fontsize',6);
if n ==1 && m ==1
    title([int2str(timoffrot(nd,1)),'.',int2str(timoffrot(nd,2)),'  p.',int2str(ifig)])
end
box on
set(gca,'XTick',[is (is+ien)/2],'fontsize',6);
%%%%%%% Below: Find the max of the hilbert transform within the arrival; then walk out from there +/-
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
%%%%%%% Above: Find the max of the hilbert transform within the arrival; then walk out from there +/-

%% subplot 2
%%% plot subfigures at 2nd row (0*3+2), 5th row (1*3+2), (n-1)*3+2
subplot(3, 1, 2, 'align');
plot(STA1bbfile(winlen*(nt-1)+1:winlen*nt,1),STA1bbfile(winlen*(nt-1)+1:winlen*nt,2),'r')
hold on
plot(STA2bbfile(winlen*(nt-1)+1:winlen*nt,1),STA2bbfile(winlen*(nt-1)+1:winlen*nt,2),'b')
plot(STA3bbfile(winlen*(nt-1)+1:winlen*nt,1),STA3bbfile(winlen*(nt-1)+1:winlen*nt,2),'k')
is = STA1bbfile(winlen*(nt-1)+1,1);
ien= STA1bbfile(winlen*nt,1);
%                         yma=max(max([STA1bbfile(winlen*(nt-1)+1:winlen*nt,2) STA2bbfile(winlen*(nt-1)+1:winlen*nt,2) ...
%                             STA3bbfile(winlen*(nt-1)+1:winlen*nt,2)]));
%                         ymi=min(min([STA1bbfile(winlen*(nt-1)+1:winlen*nt,2) STA2bbfile(winlen*(nt-1)+1:winlen*nt,2) ...
%                             STA3bbfile(winlen*(nt-1)+1:winlen*nt,2)]));
%                         xvect=[is is+2*(yma-ymi)*(winlen/160.)]; %amplitude bar originally scaled for 4-s window
%                         yma=2.4*max(yma,-ymi);
%                         yvect=[-0.9*yma -0.9*yma];
%                         plot(xvect,yvect,'r','linewidth',3)
plot([is+1/hibb is+1/lobb],[-0.8*yma -0.8*yma],'k','linewidth',3)
text(is+0.2, 0.66*yma, int2str(STA13file(nt,1)),'fontsize',6);
text(ien-0.6, 0.66*yma, int2str(STA12file(nt,1)),'fontsize',6);
box on
axis([is ien -yma yma])
set(gca,'XTick',[is (is+ien)/2],'fontsize',6);

%% subplot 3
%%% plot subfigures at 3rd row (0*3+3), 6th row (1*3+3), (n-1)*3+3
subplot(3, 1, 3,'align');
STA12tr=STA1file(winlen*(nt-1)+1:winlen*nt,2).*STA2file(winlen*(nt-1)+1:winlen*nt,2);
STA13tr=STA1file(winlen*(nt-1)+1:winlen*nt,2).*STA3file(winlen*(nt-1)+1:winlen*nt,2);
STA23tr=STA2file(winlen*(nt-1)+1:winlen*nt,2).*STA3file(winlen*(nt-1)+1:winlen*nt,2);
%                     % find the peaks etc.               %Commented out 9/11/18
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
%                     ccdat(nt,1)=timoffrot(nd,1);
%                     ccdat(nt,2)=timoffrot(nd,2);
%                     ccdat(nt,3)=ifig;
%                     ccdat(nt,4)=STA1file(winlen*(nt-1)+1,1);
%                     ccdat(nt,5)=STA13file(nt,1);
%                     ccdat(nt,6)=STA12file(nt,1);
%                     ccdat(nt,7)=Prev(nt);
%                     ccdat(nt,8)=(STA13file(nt,2)+STA12file(nt,2)+STA32file(nt,2))/3;
%                     ccdat(nt,9)=meanalln;
%                     ccdat(nt,10)=STA32file(nt,1);
%                     ccdat(nt,11)=2.5*Ampsq(nt)/Prev(nt);
%                     ccdat(nt,12)=match(nt,1)-match(nt,2);

alln(alln<0)=-10^4*yma; %just so they don't plot.
%%idiff=STA32file(nt,3);
%%maxxc=max(alln(idiff:idiff+cclen)); %endpoints of alln window should be endpoints of 1 sec cumsum interval
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
% %                     STA12file(nin+1,1:2)=[imaxSTA12wr xcmaxSTA12n(n)];
% %                     STA13file(nin+1,1:2)=[imaxSTA13wr xcmaxSTA13n(n)];
% %                     STA32file(nin+1,1:3)=[cumsumtrdiff/cumsumtr(winlen) xcmaxSTA32n(n) idiff];
%                     text(is+0.1, ymi+0.82*(yma-ymi), num2str(STA13file(nt,2),2),'fontsize',6);
%                     text(is+0.1, ymi+0.64*(yma-ymi), num2str(STA12file(nt,2),2),'fontsize',6);
%                     text(is+0.1, ymi+0.46*(yma-ymi), num2str(STA32file(nt,2),2),'fontsize',6);
text(is+0.1, ymi+0.82*(yma-ymi), num2str((STA13file(nt,2)+STA12file(nt,2)+STA32file(nt,2))/3,2),'fontsize',6);
text(ien-0.6, ymi+0.1*(yma-ymi), num2str(STA32file(nt,1),2),'fontsize',6);
%                     text(ien-2.2, ymi+0.1*(yma-ymi), num2str(STAamp(nt,2),2),'fontsize',6);
%                     text(ien-1.4, ymi+0.1*(yma-ymi), num2str(STAamp(nt,3),2),'fontsize',6);
%axis([lo-1 hi+1 ymi yma])
%axis tight
%set(gca,'XTick',[is is+2],'fontsize',6);
%                   pkfile(nt,:)=[STA1file(winlen*(nt-1)+1+(winlen/2)) pks(imaxpk,1) maxpk pksid12 pksid13 pkwid rat14];

drawnow
orient TALL
% print('-depsc2','-r600', [getenv('ALLAN'),'/WIGS/',IDENTIF,'-',num2str(lo),'-',num2str(hi),'ms',int2str(mshift),'-',int2str(winlen/sps),'s','.','one-win-example','.eps'])
