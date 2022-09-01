%Reads in synthetic seismogram. Does instantaneous amp and frequency stuff.
format short e
clear all 
close all 
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');


hi=6;
lo=1.5;
% hi2nd=hi; %4.5;
% lo2nd=lo;

%Read data:
seism=load('BOSTOCK/NEW/WFPV/synth_beta9999_Twin600_N100000');
seism=seism(30*40:600*40);
tracelen=length(seism);

% %cosine taper before filtering:
% x=(0:pi/80:pi/2-pi/80)';
% PGCN(1:40)=sin(x).*PGCN(1:40);
% x=flipud(x);
% PGCN(tracelenPGC-39:tracelenPGC)=sin(x).*PGCN(tracelenPGC-39:tracelenPGC);
% 
% %Filter data:
% npo=2;
% npa=1;
% [PGCEf]=1.6e-4*bandpass(PGCE,40,lo,hi,npo,npa,'butter');
% [PGCNf]=1.6e-4*bandpass(PGCN,40,lo,hi,npo,npa,'butter');

scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite;

%for n=1:size(wins,1)
    istart=30*40;
    iend=600*40;
% instantaneous amplitude and frequecy using hi2nd, lo2nd filtering.
    hSI=hilbert(seism); %For amplitude information, keep it as before
    SIamp=abs(hSI);
    xprime=0*(1:tracelen);
    yprime=0*(1:tracelen);
    xprime(2:tracelen-1)=real(hSI(3:end))-real(hSI(1:end-2));
    xprime(1)=real(hSI(2))-real(hSI(1));
    xprime(tracelen)=real(hSI(tracelen))-real(hSI(tracelen-1));
    yprime(2:tracelen-1)=imag(hSI(3:end))-imag(hSI(1:end-2));
    yprime(1)=imag(hSI(2))-imag(hSI(1));
    yprime(tracelen)=imag(hSI(tracelen))-imag(hSI(tracelen-1));
    SIfreq=((40/2)/(2*pi))*(seism.*yprime'-xprime'.*imag(hSI))./abs(hSI).^2; %divided by 2 b.c. of centered f.d.
    loglog(SIamp,SIfreq,'bo','markersize',1)
    xlim([.4 40])
    ylim([0.2 20])

%     npages=seclen/30;
%     for j=1:npages
%         %h=figure('Position',[wid/3 1 2.5*wid hite]); %UNCOMMENT
%         istart=plotstart+(j-1)*shortlen+1; %added +1 late
%         iend=istart+shortlen;
%         inormstart=istart-plotstart;
%         inormend=inormstart+shortlen;
%         iscmplx=istart-plotstart;
%         iecmplx=iscmplx+shortlen;
%         ltmax=max([realSSIB(istart-nshift:iend+mshift); realSILB(istart:iend)]);
%         ltmin=min([realSSIB(istart-nshift:iend+mshift); realSILB(istart:iend)]);
%         ltmax=1.5*max(ltmax,-ltmin);
%         [dummy,kmax]=max(cumSISSn(iecmplx,:)-cumSISSn(iscmplx,:));
%         runningcc(iscmplx:iecmplx)=SISSn(iscmplx:iecmplx,kmax); %make the running cc equal to that of the most coherent 30 sec.
%         kmax=kmax-1-nshift %to center
%         runningSSamp(iscmplx:iecmplx)=SSamp(iscmplx+nshift+kmax:iecmplx+nshift+kmax); %to account for shift of SSIB w.r.t SILB
%         runningSSfreq(iscmplx:iecmplx)=SSfreq(iscmplx+nshift+kmax:iecmplx+nshift+kmax); %same as above
% %         for k=-nshift:mshift
% %             subplot(nshift+mshift+2,1,k+nshift+2,'align'); 
% %             hold on
% %             for inr=inrubstart:inrubend  %inefficient
% %                 if k == -nshift
% %                     if round(rubins(inr,2)-rubins(inr,3)) < k
% %                         irstart=(rubins(inr,1)-2)*40-round(rubins(inr,2));
% %                         irend=irstart+4*40-1;
% %                         plot(timsPGC(irstart:irend),0.9*ltmax,'r','linewidth',3)
% %                     end
% %                 end
% %                 if (k == mshift) && (round(rubins(inr,2)-rubins(inr,3)) > k)
% %                     irstart=(rubins(inr,1)-2)*40-round(rubins(inr,2));
% %                     irend=irstart+4*40-1;
% %                     plot(timsPGC(irstart:irend),0.9*ltmax,'r','linewidth',3)
% %                 end
% %                 if round(rubins(inr,2)-rubins(inr,3))==k 
% %                     irstart=(rubins(inr,1)-2)*40-round(rubins(inr,2));
% %                     irend=irstart+4*40-1;
% %                     plot(timsPGC(irstart:irend),realPGC(irstart+round(rubins(inr,2)):irend+round(rubins(inr,2))),'r')
% %                     text(timsPGC(irstart),-0.9*ltmax,int2str(round(rubins(inr,2))),'fontsize',6);
% %                 end
% %             end
% %             if k==kmax
% %                 plot(timsPGC(istart:iend),SIamp(iscmplx:iecmplx),'g');
% %             end
% %             plot(timsPGC(istart:iend),SISSn(inormstart:inormend,k+nshift+1)*ltmax,'co','markersize',1)
% %             if k==0
% %                 for inb=inbosstart:inbosend
% %                     plot(bostsec,0.85*ltmax,'ko','MarkerSize',4,'MarkerFaceColor','k')
% %                 end
% %             end
% %             plot(timsPGC(istart:iend),realSILB(istart:iend),'k');
% %             plot(timsPGC(istart:iend),realSSIB(istart+k:iend+k),'b');
% %             axis([timsPGC(istart) timsPGC(iend) -ltmax ltmax]);
% %             box on
% %        end
% %        subplot(nshift+mshift+2,1,1,'align'); 
% %        hold on
% %        %semilogy(timsPGC(istart:iend),SIfreq(iscmplx:iecmplx),'b'); %DON'T KNOW WHY THIS DOESN'T WORK! 
% %        plot(timsPGC(istart:iend),SIfreq(iscmplx:iecmplx),'b');
% %        plot(timsPGC(istart:iend),SIamp(iscmplx:iecmplx)./SIortamp(iscmplx:iecmplx),'r');
% %        axis([timsPGC(istart) timsPGC(iend) 0.5 10]);
% %        set(gca,'ytick',[1 3 5 10]);
% %        title([IDENTIF,'  ',num2str(lo),'-',num2str(hi),' Hz'])
% %        box on
% %        
% %        set(h,'PaperPosition',[0.25 0.25 8 10.5])
% %        orient landscape
% %        %if j <= 9
% %        print(h,'-depsc',['ARMMAP/CMPWIGS/',YEAR,'.',JDAY,'.',int2str(timsPGC(istart)),'_',num2str(lo),'-',num2str(hi),'eye.eps'])
% %        %else
% %        %    print(h,'-depsc',['ARMMAP/WIGS/',YEAR,'.',JDAY,'.',int2str(wins(n,1)),'_',num2str(lo),'-',num2str(hi),'_',int2str(j),'eye.eps'])
% %        %end
%     end
%     cmplxseisms(ndata+1:ndata+longlen,1)=SIamp;
%     cmplxseisms(ndata+1:ndata+longlen,2)=SIfreq;
%     cmplxseisms(ndata+1:ndata+longlen,3)=SIamp./SIortamp;
%     cmplxseisms(ndata+1:ndata+longlen,4)=runningcc;
%     cmplxseisms(ndata+1:ndata+longlen,5)=runningSSamp;
%     cmplxseisms(ndata+1:ndata+longlen,6)=runningSSfreq;
%     ndata=ndata+longlen;
% %end
% aa=cmplxseisms;
% tf1=cmplxseisms(:,1)<1.e-3;
% tf2=cmplxseisms(:,1)>10;
% tf3=cmplxseisms(:,2)<0.1;
% tf4=cmplxseisms(:,2)>10;
% tf5=cmplxseisms(:,3)<1.5;
% TFall=tf1 | tf2 | tf3 | tf4 | tf5 ;
% aa(tf5,:)=[];
% aa2=aa;
% tf6=aa(:,4)<0.5;
% aa2(tf6,:)=[];
% aa=sortrows(aa,-1);
% aa2=sortrows(aa2,-1);
% %values=hist3(aa(:,1:2),[201 201]);
% %xb = linspace(min(aa(:,1)),max(aa(:,1)),size(values,1)+1);
% %yb = linspace(min(aa(:,2)),max(aa(:,2)),size(values,1)+1);
% %imagesc(xb,yb,values')
% figure
% loglog(cmplxseisms(:,1),cmplxseisms(:,2),'bo','markersize',1)
% hold on
% loglog(aa(:,1),aa(:,2),'ro','markersize',1)
% loglog(aa2(:,1),aa2(:,2),'go','markersize',1)
% in=300; %without cc quality control
% navs=floor(size(aa,1)/in);
% x=0*(1:navs);
% y=zeros(navs,3);
% for i=1:navs
%     x(i)=median(aa((i-1)*in+1:i*in,1));
%     y(i,:)=prctile(aa((i-1)*in+1:i*in,2),[25 50 75]);
% end
% loglog(x,y(:,2),'c','linewidth',2)
% loglog(x,y(:,1),'c','linewidth',1)
% loglog(x,y(:,3),'c','linewidth',1)
% in=200;
% navs=floor(size(aa2,1)/in);
% x=0*(1:navs);
% y=zeros(navs,3);
% for i=1:navs
%     x(i)=median(aa2((i-1)*in+1:i*in,1));
%     y(i,:)=prctile(aa2((i-1)*in+1:i*in,2),[25 50 75]);
% end
% loglog(x,y(:,2),'k','linewidth',2)
% loglog(x,y(:,1),'k','linewidth',1)
% loglog(x,y(:,3),'k','linewidth',1)
% medy=median(y(1:floor(0.2*length(y)),2));  %amps (x's) a;ready sorted high to low.
% med_o_cent=medy/sqrt(hi*lo);
% text(0.2,0.864,['median (x > ',num2str(x(floor(0.2*length(x)))),') = ',num2str(medy)],'fontsize',7);
% text(0.2,0.72,['median/central = ',num2str(medy/sqrt(hi*lo))],'fontsize',7);
% text(0.2,0.6,'[Log(median)-Log(central)]','fontsize',7);
% text(0.2,0.5,['/[Log(hi)-Log(central)] = ',num2str(log10(med_o_cent)/log10(sqrt(hi/lo)))],'fontsize',7);
% loline=[1.e-5 lo2nd;
%         1.e5 lo2nd];
% hiline=[1.e-5 hi2nd;
%         1.e5 hi2nd];
% loglog(loline(:,1),loline(:,2),'k--')
% loglog(hiline(:,1),hiline(:,2),'k--')
% xlim([0.001 2])
% ylim([0.3 30])
% xlabel('log envelope amplitude')
% ylabel('log frequency (Hz)')
% title(['SILB  ',YEAR,'.',JDAY,'  ',num2str(lo),'-',num2str(hi),' Hz'])
% print('-depsc',['Amp_v_FreqSILB_',YEAR,'.',JDAY,'_',num2str(lo),'-',num2str(hi),'.eps'])
% 
% aa=sortrows(aa,-5);
% aa2=sortrows(aa2,-5);
% figure
% loglog(cmplxseisms(:,5),cmplxseisms(:,6),'bo','markersize',1)
% hold on
% loglog(aa(:,5),aa(:,6),'ro','markersize',1)
% loglog(aa2(:,5),aa2(:,6),'go','markersize',1)
% in=300; %without cc quality control
% navs=floor(size(aa,1)/in);
% x=0*(1:navs);
% y=zeros(navs,3);
% for i=1:navs
%     x(i)=median(aa((i-1)*in+1:i*in,5));
%     y(i,:)=prctile(aa((i-1)*in+1:i*in,6),[25 50 75]);
% end
% loglog(x,y(:,2),'c','linewidth',2)
% loglog(x,y(:,1),'c','linewidth',1)
% loglog(x,y(:,3),'c','linewidth',1)
% in=200;
% navs=floor(size(aa2,1)/in);
% x=0*(1:navs);
% y=zeros(navs,3);
% for i=1:navs
%     x(i)=median(aa2((i-1)*in+1:i*in,5));
%     y(i,:)=prctile(aa2((i-1)*in+1:i*in,6),[25 50 75]);
% end
% loglog(x,y(:,2),'k','linewidth',2)
% loglog(x,y(:,1),'k','linewidth',1)
% loglog(x,y(:,3),'k','linewidth',1)
% medy=median(y(1:floor(0.2*length(y)),2));  %amps (x's) a;ready sorted high to low.
% med_o_cent=medy/sqrt(hi*lo);
% text(0.2,0.864,['median (x > ',num2str(x(floor(0.2*length(x)))),') = ',num2str(medy)],'fontsize',7);
% text(0.2,0.72,['median/central = ',num2str(medy/sqrt(hi*lo))],'fontsize',7);
% text(0.2,0.6,'[Log(median)-Log(central)]','fontsize',7);
% text(0.2,0.5,['/[Log(hi)-Log(central)] = ',num2str(log10(med_o_cent)/log10(sqrt(hi/lo)))],'fontsize',7);
% loline=[1.e-5 lo2nd;
%         1.e5 lo2nd];
% hiline=[1.e-5 hi2nd;
%         1.e5 hi2nd];
% loglog(loline(:,1),loline(:,2),'k--')
% loglog(hiline(:,1),hiline(:,2),'k--')
% xlim([0.001 2])
% ylim([0.3 30])
% xlabel('log envelope amplitude')
% ylabel('log frequency (Hz)')
% title(['SSIB  ',YEAR,'.',JDAY,'  ',num2str(lo),'-',num2str(hi),' Hz'])
% print('-depsc',['Amp_v_FreqSSIB_',YEAR,'.',JDAY,'_',num2str(lo),'-',num2str(hi),'.eps'])
