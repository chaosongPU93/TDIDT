% test for plotting only one time winow in chaoSTA

close all
figure('Position',[scrsz(3)/5 scrsz(4)/5 1.5*scrsz(3)/5 2.5*scrsz(4)/5]);

%% subplot 1
subplot(3, 1, 1,'align');
nt=1;
%%% STA1file(*, 2) stores the seismogram of each window
yma=max(max([STA1file(winlen*(nt-1)+1:winlen*nt,2) STA2file(winlen*(nt-1)+1:winlen*nt,2) ...
STA3file(winlen*(nt-1)+1:winlen*nt,2)]));
ymi=min(min([STA1file(winlen*(nt-1)+1:winlen*nt,2) STA2file(winlen*(nt-1)+1:winlen*nt,2) ...
STA3file(winlen*(nt-1)+1:winlen*nt,2)]));
yma=2.4*max(yma,-ymi);
ymakeep=yma/2.4;
% Lines below plot + or - template:
%                         if match(nt,1) >= match(nt,2) %if positive cc with template is larger ...
%                             plot(timstemp((nt-1)*templen+1:nt*templen),STAtemps(whichtoplot,:)*yma/2.4,'c','linewidth',2)
%                         else %if cc with negative template is larger ...
%                             plot(timstempneg((nt-1)*templen+1:nt*templen),-STAtemps(whichtoplot,:)*yma/2.4,'m','linewidth',2)
%                         end
% Lines above plot + or - template:
hold on
%%% plot the seismogram at 3 stations of each win
plot(STA1file(winlen*(nt-1)+1:winlen*nt,1),STA1file(winlen*(nt-1)+1:winlen*nt,2),'r')
plot(STA2file(winlen*(nt-1)+1:winlen*nt,1),STA2file(winlen*(nt-1)+1:winlen*nt,2),'b')
plot(STA3file(winlen*(nt-1)+1:winlen*nt,1),STA3file(winlen*(nt-1)+1:winlen*nt,2),'k')
is = STA1file(winlen*(nt-1)+1,1);   % start time of each win
ien= STA1file(winlen*nt,1);         % end time of each win
axis([is ien -yma yma])
% plot amplitude scaling ?
xvect=[is is+2*(yma-ymi)*(winlen/160.)]; %*mean(scalefact)]; %amplitude bar originally scaled for 4-s window.  Not sure "mean" is necessary.
yvect=[-0.9*yma -0.9*yma];
plot(xvect,yvect,'r','linewidth',3)
% plot freq. band
plot([is+1/hi is+1/lo],[-0.8*yma -0.8*yma],'k','linewidth',3)
% STA32file(nt,3), 3rd col. -> idiff, /sps to time instead of samples
xcvect=[is+STA32file(nt,3)/sps is+(STA32file(nt,3)+cncntr)/sps];
ycvect=[0.95*yma 0.95*yma];
plot(xcvect,ycvect,'b','linewidth',3)
% plot Bostock's detections (in sec), red circle
plot(bostsec,0.93*yma,'ro','MarkerSize',4,'MarkerFaceColor','r')
% plot text of max cc coeff val
% match(1/2) -> max cc coeff between positive/negative template and the window with the main arrival, 2s win, 120 samples
text(is+STA32file(nt,3)/sps,  0.58*yma, num2str(match(nt,1),2),'fontsize',6,'color','b');
text(is+STA32file(nt,3)/sps, -0.75*yma, num2str(match(nt,2),2),'fontsize',6,'color','r');
% makes more sense now, because STA13file(nt,1) -> imaxSTA13wr, which stores the max cc shift of
% qualified window only, whose size is smaller than nwin
text(is+0.2, 0.66*yma, int2str(STA13file(nt,1)),'fontsize',6);
text(ien-0.6, 0.66*yma, int2str(STA12file(nt,1)),'fontsize',6);

if n ==1 && m ==1
    title([int2str(timoffrot(nd,1)),'.',int2str(timoffrot(nd,2)),'  p.',int2str(ifig)])
end
box on
set(gca,'XTick',[is (is+ien)/2],'fontsize',6);

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
% plot dot product of 1 and 2, which is already shifted
STA12tr=STA1file(winlen*(nt-1)+1:winlen*nt,2).*STA2file(winlen*(nt-1)+1:winlen*nt,2);
STA13tr=STA1file(winlen*(nt-1)+1:winlen*nt,2).*STA3file(winlen*(nt-1)+1:winlen*nt,2);
STA23tr=STA2file(winlen*(nt-1)+1:winlen*nt,2).*STA3file(winlen*(nt-1)+1:winlen*nt,2);
% find the peaks etc.
avedots=(STA12tr+STA13tr+STA23tr)/3.;
[peaks, locs]=findpeaks(avedots,'minpeakdistance',3);
npeaks=length(peaks);
[maxpk, imaxpk]=max(peaks);
pks=zeros(npeaks,2);
pks(:,2)=peaks;     % store peaks val in 2nd col.
pks(:,1)=STA1file(winlen*(nt-1)+locs);      % store global ind of loc in 1st col.
pksort=sortrows(pks,2);     % sort pks according to peaks val in ascending order
rat14=maxpk/pksort(npeaks-4,2); %ratio of max to 4th largest anywhere in window, 5th instead of 4th???
if imaxpk==1
    maxpkp=peaks(2);
    maxpkm=-9e9;
    pkwid=2*(pks(2,1)-pks(1,1));
    pksid12=maxpk/maxpkp;
    pksid13=maxpk/peaks(3);
elseif imaxpk==npeaks
    maxpkp=-9e9;
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
%%% cumsumSTA12tr is cum sum of dot product
cumsumSTA12tr=cumsum(STA12tr);
cumsumSTA13tr=cumsum(STA13tr);
cumsumSTA23tr=cumsum(STA23tr);
%                     yma=1.1; %yma=max(avedots);
%                     ymi=-0.1; %ymi=min(avedots);
% The following for running cc
cclen=20; %running cc window length, in samples, 0.5s
yma=max(max([cumsumSTA12tr cumsumSTA13tr cumsumSTA23tr]));   % here, yma & ymi are redefined
ymi=min(min([cumsumSTA12tr cumsumSTA13tr cumsumSTA23tr]));
%%% auto dot product, ST1auto can be compared to STA12tr
ST1auto=STA1file(winlen*(nt-1)+1:winlen*nt,2).*STA1file(winlen*(nt-1)+1:winlen*nt,2);
ST1sq=cumsum(ST1auto);  % cum sum of squre
ST2auto=STA2file(winlen*(nt-1)+1:winlen*nt,2).*STA2file(winlen*(nt-1)+1:winlen*nt,2);
ST2sq=cumsum(ST2auto);  % ST2sq can be compared with cumsumSTA12tr
ST3auto=STA3file(winlen*(nt-1)+1:winlen*nt,2).*STA3file(winlen*(nt-1)+1:winlen*nt,2);
ST3sq=cumsum(ST3auto);
ST12num=cumsumSTA12tr(cclen+1:winlen)-cumsumSTA12tr(1:winlen-cclen);  % length of cumsumSTA12tr is (winlen)
ST13num=cumsumSTA13tr(cclen+1:winlen)-cumsumSTA13tr(1:winlen-cclen);  % result in sum of square of cclen long
ST23num=cumsumSTA23tr(cclen+1:winlen)-cumsumSTA23tr(1:winlen-cclen);  % window from [2, cclen+1]
ST1den=ST1sq(cclen+1:winlen)-ST1sq(1:winlen-cclen);
ST2den=ST2sq(cclen+1:winlen)-ST2sq(1:winlen-cclen);
ST3den=ST3sq(cclen+1:winlen)-ST3sq(1:winlen-cclen);
ST12n=ST12num./realsqrt(ST1den.*ST2den);
ST13n=ST13num./realsqrt(ST1den.*ST3den);
ST23n=ST23num./realsqrt(ST2den.*ST3den);  % length of alln, ST12n, ST1den, ST12num are (winlen-cclen)
alln=(ST12n+ST13n+ST23n)/3;   % alln --> [-1,1]
alln(alln<0)=-10*yma;   %just so they don't plot. %% reassign those alln<0 to be larger then yaxis
%%idiff=STA32file(nt,3);
%%maxxc=max(alln(idiff:idiff+cclen)); %endpoints of alln window should be endpoints of 1 sec cumsum interval
plot(STA1file(winlen*(nt-1)+cclen/2+1:winlen*nt-cclen/2,1),(yma+ymi)/2+(yma-ymi)*alln/2,'co','markersize',1)
hold on
plot(STA1file(winlen*(nt-1)+1:winlen*nt,1), (yma+ymi)/2+zeros(winlen,1), 'k--');
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
%                     if ismember(round(STA1file(winlen*(nt-1)+1)),displaytim)
%                         %%%%This is to write stuff
%                         i1=winlen*(nt-1)+1; i2=winlen*nt;
%                         cumsumSTA=cumsumSTA12tr+cumsumSTA13tr+cumsumSTA23tr;
%                         cumsumSTA=cumsumSTA/cumsumSTA(end);
%                         datfile(1:9,winlen*(writeouts-1)+1:winlen*writeouts)=[STA1file(i1:i2,1)';
%                                                                               STA1file(i1:i2,2)';
%                                                                               STA2file(i1:i2,2)';
%                                                                               STA3file(i1:i2,2)';
%                                                                               STA1bbfile(i1:i2,2)';
%                                                                               STA2bbfile(i1:i2,2)';
%                                                                               STA3bbfile(i1:i2,2)';
%                                                                               cumsumSTA';
%                                                                               -99*ones(1,winlen)];
%                         datfile(9,winlen*(writeouts-1)+1+cclen/2:winlen*writeouts-cclen/2)=(ST12n+ST13n+ST23n)/3;
%                         guidefile(writeouts,:)=[writeouts STA1file(i1,1) ymakeep];
%                         %%%%That was to write stuff
%                     end
axis([is ien ymi yma])
set(gca,'XTick',[is (is+ien)/2],'fontsize',6);
% The above for running-sum dot product
box on
%%% text from up to down is max cc coeff of 13, 12, 32
%%% text from right to left is normalized cumsumtrdiff of 32, 80 percentiles of data 3, data 2
text(is+0.1, ymi+0.82*(yma-ymi), num2str(STA13file(nt,2),2),'fontsize',6);
text(is+0.1, ymi+0.64*(yma-ymi), num2str(STA12file(nt,2),2),'fontsize',6);
text(is+0.1, ymi+0.46*(yma-ymi), num2str(STA32file(nt,2),2),'fontsize',6);
text(ien-0.6, ymi+0.1*(yma-ymi), num2str(STA32file(nt,1),2),'fontsize',6);
text(ien-2.2, ymi+0.1*(yma-ymi), num2str(STAamp(nt,2),2),'fontsize',6);
text(ien-1.4, ymi+0.1*(yma-ymi), num2str(STAamp(nt,3),2),'fontsize',6);
%axis([lo-1 hi+1 ymi yma])
%axis tight
%set(gca,'XTick',[is is+2],'fontsize',6);
%                     pkfile(nt,:)=[STA1file(winlen*(nt-1)+1+(winlen/2)) pks(imaxpk,1) maxpk pksid12 pksid13 pkwid rat14];
%%% pkfile:  time   index of maxpk     maxpk val
%%% pksid12=maxpk/max(maxpkp,maxpkm)    pksid13=maxpk/min(maxpkp,maxpkm)
%%% pkwid=pks(imaxpk+1,1)-pks(imaxpk-1,1)   rat14=maxpk/pksort(npeaks-4,2)

drawnow
% orient TALL
print('-depsc', [workpath,'/PGCtrio/WIGS/',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s_',num2str(sps),'sps','_','one-win-example.eps']);
