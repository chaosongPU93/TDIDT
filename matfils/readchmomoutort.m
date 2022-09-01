%                 traces(nttot+nt,1:160)=SILBtr;
%                 traces(nttot+nt,161:320)=SSIBtr;
%                 traces(nttot+nt,321:331)=[year+0.001*jday bostsec(nt) bostmom(nttot+nt) maxSISS(nttot+nt) ...
%                     maxccSItemp(nttot+nt) maxccSStemp(nttot+nt) ...
%                     (imaxSIstretch(nttot+nt)+19)/40 (imaxSSstretch(nttot+nt)+19)/40 ...
%                     scaleSIstore(nttot+nt) scaleSSstore(nttot+nt) maxSISS(nttot+nt)*maxccSItemp(nttot+nt)];
%                 traces(nttot+nt,332:491)=SILBortr;
%                 traces(nttot+nt,492:651)=SSIBortr;
%                 traces(nttot+nt,652:811)=scaleSItemp*SILBtempsplot;
%                 %traces(nttot+nt,652:811)=KLNBtr;
%                 %traces(nttot+nt,812:971)=KLNBortr;
close all
clear all
lo=1;
hi=8;
npa=2;
load('checkmomout.mat')
traceplot=sortrows(traces,-331);
ntot=size(traceplot,1)
totSILB=0;
totSSIB=0;
totSISS=0;
intot=0; %the number of Bostock detections in the time windows I've chosen
totinSILB=0;
totinSSIB=0;
totinSISS=0;
nt=0; %nt is "local", drops to zero at the end of each day
nrow=3;
mcol=4;
scrsz=get(0,'ScreenSize');
for ifig=1:floor(ntot/(nrow*mcol))+1
    h=figure('Position',[scrsz(3)/10 scrsz(4)/10 4*scrsz(3)/5 9*scrsz(4)/10]);
    for n = 1:nrow
        for m = 1:mcol
            if nt < ntot 
                nt=nt+1;
                wins=windows(traceplot(nt,321));
                nwins=size(wins,1);
                inwin=0;
                for i=1:nwins
                    if traceplot(nt,322)>=wins(i,1) && traceplot(nt,322)<=wins(i,2)
                        inwin=1;
                        intot=intot+1
                    end
                end
                subplot(2*nrow,mcol,2*(n-1)*mcol+m,'align');
                maxmax=max([abs(traceplot(nt,1:320)) abs(traceplot(nt,652:811))]);
                hold on
                plot(traceplot(nt,1:160),'b')
                plot(traceplot(nt,161:320),'k')
                plot(traceplot(nt,652:811),'r--')
                title([num2str(traceplot(nt,321)),'  ',int2str(traceplot(nt,322)),'  M',num2str(traceplot(nt,323)), ...
                    '  cc',num2str(traceplot(nt,324))],'fontsize',8)
                text(5,0.9*maxmax,int2str(inwin),'fontsize',12)
                text(100,-0.98*maxmax,[num2str(traceplot(nt,325)),' SI'],'fontsize',7)
                text(100,-0.78*maxmax,[num2str(traceplot(nt,326)),' SS'],'fontsize',7)
                text(5,-0.98*maxmax,[num2str(traceplot(nt,327)),'  ',num2str(1000*traceplot(nt,329))],'fontsize',7)
                text(5,-0.78*maxmax,[num2str(traceplot(nt,328)),'  ',num2str(1000*traceplot(nt,330))],'fontsize',7)
                ylim([-1.1*maxmax 1.1*maxmax])
                xlim([0 160])
                set(gca,'xtick',[40 80 120],'fontsize',8);
                box on
                subplot(2*nrow,mcol,2*(n-1)*mcol+mcol+m,'align');
                hold on
                plot(traceplot(nt,332:491),'b')
                plot(traceplot(nt,492:651),'k')
                text(20,-0.9*maxmax,'orthogonal','fontsize',7)
                ylim([-1.1*maxmax 1.1*maxmax])
                xlim([0 160])
                set(gca,'xtick',[40 80 120],'fontsize',8);
                box on
                
                SILBtr=traceplot(nt,1:160);
                %SILBtr2=SILBtr(61:100).*SILBtr(61:100);
                SILBtr2=SILBtr(66:95).*SILBtr(66:95);
                %SILBtr2=SILBtr(71:90).*SILBtr(71:90);
                sumSILBtr2=sum(SILBtr2);
                totSILB=totSILB+sumSILBtr2;
                if inwin == 1
                    totinSILB=totinSILB+sumSILBtr2;
                end
                
                SSIBtr=traceplot(nt,161:320);
                %SSIBtr2=SSIBtr(61:100).*SSIBtr(61:100);
                SSIBtr2=SSIBtr(66:95).*SSIBtr(66:95);
                %SSIBtr2=SSIBtr(71:90).*SSIBtr(71:90);
                sumSSIBtr2=sum(SSIBtr2);
                totSSIB=totSSIB+sumSSIBtr2;
                if inwin == 1
                    totinSSIB=totinSSIB+sumSSIBtr2;
                end
                
                %SISS=SILBtr(61:100).*SSIBtr(61:100);
                SISS=SILBtr(66:95).*SSIBtr(66:95);
                %SISS=SILBtr(71:90).*SSIBtr(71:90);
                sumSISS=sum(SISS);
                totSISS=totSISS+sumSISS;
                if inwin == 1
                    totinSISS=totinSISS+sumSISS;
                end
            end
        end
    end
    set(h,'PaperPosition',[0.25 0.25 8 10.5])
    orient landscape
    if ifig <= 9
        print(h,'-depsc',['BFIGS/',num2str(lo),'-',int2str(hi),'_',int2str(npa),'pass_',int2str(0),int2str(0),int2str(ifig),'.eps'])
    elseif ifig <= 99
        print(h,'-depsc',['BFIGS/',num2str(lo),'-',int2str(hi),'_',int2str(npa),'pass_',int2str(0),int2str(ifig),'.eps'])
    else
        print(h,'-depsc',['BFIGS/',num2str(lo),'-',int2str(hi),'_',int2str(npa),'pass_',int2str(ifig),'.eps'])
    end
end
totSILB
totinSILB
totSSIB
totinSSIB
totSISS
totinSISS