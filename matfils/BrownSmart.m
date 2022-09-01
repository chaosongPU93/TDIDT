%Sits at one spot (rotations; delays) and looks for the cross-correlation.
clear all
close all
format short e
t=cputime; tt=cputime;

IDPGC='PGC2005.254.85.20.80.115.55_2-6-4s';
IDSSIB='SSIB2005.254.85.20.80.115.55_2-6-4s.shift';
IDSILB='SILB2005.254.85.20.80.115.55_2-6-4s.shift';
% IDPGC='PGC2005.254.85.20.80.115.55_1-8-4s';
% IDSSIB='SSIB2005.254.85.20.80.115.55_1-8-4s.shift';
% IDSILB='SILB2005.254.85.20.80.115.55_1-8-4s.shift';
PGC=load(IDPGC);
SSIB=load(IDSSIB);
SILB=load(IDSILB);
timesPGC = PGC(:,1); datPGC = PGC(:,2);
timesSSIB = SSIB(:,1); datSSIB = SSIB(:,2); offsSSIB=SSIB(:,3);
timesSILB = SILB(:,1); datSILB = SILB(:,2); offsSILB=SILB(:,3);
%offSILB(168361:168400)
% [timsSSIB datSSIB]=load(IDSSIB);
% [timsSILB datSILB]=load(IDSILB);
timstart=34460; timend=34980; delt=timend-timstart; %day 254
timstart=34600; timend=34720; delt=timend-timstart; %day 254
% timstart=35044; timend=36412; delt=timend-timstart; %day 255; TOO BIG
% timstart=35044; timend=35880; delt=timend-timstart; %day 255; less big
i=1;
while timesPGC(i) < timstart
    lookstart=i;
    i=i+1;
end
lookstart=lookstart+1;
while timesPGC(i) <= timend
    lookend=i;
    i=i+1;
end
realPGC=datPGC(lookstart:lookend);
realSSIB=datSSIB(lookstart:lookend);
realSILB=datSILB(lookstart:lookend);
timsPGC=timesPGC(lookstart:lookend);
timsSSIB=timesSSIB(lookstart:lookend);
timsSILB=timesSILB(lookstart:lookend);
offSSIB=offsSSIB(lookstart:lookend);
offSILB=offsSILB(lookstart:lookend);
PGCauto=realPGC.*realPGC;
PGC2=cumsum(PGCauto);
SSIBauto=realSSIB.*realSSIB;
SSIB2=cumsum(SSIBauto);
SILBauto=realSILB.*realSILB;
SILB2=cumsum(SILBauto);

lenx=length(realPGC);
%templen=3*40;
templen=4*40;
tempoff=40/2;
ntemp=floor((lenx-templen)/tempoff);
timstemp=zeros(ntemp,1);
PGCtx=zeros(lenx-templen+1,ntemp); %no +1 so you can subtract ___(j-1)
PGCtn=zeros(lenx-templen+1,ntemp);
SSIBtx=zeros(lenx-templen+1,ntemp); %no +1 so you can subtract ___(j-1)
SSIBtn=zeros(lenx-templen+1,ntemp);
SILBtx=zeros(lenx-templen+1,ntemp); %no +1 so you can subtract ___(j-1)
SILBtn=zeros(lenx-templen+1,ntemp);
alltx=zeros(lenx-templen+1,ntemp);
alltn=zeros(lenx-templen+1,ntemp);
tottn=zeros(lenx-templen+1,ntemp);
n05=0; j05=0; f05=0;
quantum=240;
for n=1:ntemp;
    istart=(n-1)*tempoff+1; %+mshift to coincide with ...?
    iend=istart+templen-1;
    if mod((iend-1),quantum) >= templen-1
        n
%         istart
%         iend
        PGCtempl=realPGC(istart:iend);
        SILBtempl=realSILB(istart:iend);
        SSIBtempl=realSSIB(istart:iend);
        PGCt2=dot(PGCtempl,PGCtempl);
        SSIBt2=dot(SSIBtempl,SSIBtempl);
        SILBt2=dot(SILBtempl,SILBtempl);
        for j=istart+40:lenx-templen+1 
            if quantum-mod((j-1),quantum) >= templen
            if abs(timsPGC(istart)-timsPGC(j))>1.5 %In case data string is not monotonically increasing; guards against xc=1
                PGCtx(j,n)=dot(PGCtempl,realPGC(j:j+templen-1));
                PGCtarg2=PGC2(j+templen-1)-PGC2(j-1);
                PGCtn(j,n)=PGCtx(j,n)/sqrt(PGCt2*PGCtarg2);
                SSIBtx(j,n)=dot(SSIBtempl,realSSIB(j:j+templen-1));
                SSIBtarg2=SSIB2(j+templen-1)-SSIB2(j-1);
                SSIBtn(j,n)=SSIBtx(j,n)/sqrt(SSIBt2*SSIBtarg2);
                SILBtx(j,n)=dot(SILBtempl,realSILB(j:j+templen-1));
                SILBtarg2=SILB2(j+templen-1)-SILB2(j-1);
                SILBtn(j,n)=SILBtx(j,n)/sqrt(SILBt2*SILBtarg2);
                alltx(j,n)=PGCtx(j,n)+SSIBtx(j,n)+SILBtx(j,n);
                alltn(j,n)=alltx(j,n)/sqrt((PGCt2+SSIBt2+SILBt2)*(PGCtarg2+SSIBtarg2+SILBtarg2));
                tottn(j,n)=(PGCtn(j,n)+SSIBtn(j,n)+SILBtn(j,n))/3.0;
                if alltn(j,n) > 0.45 && abs(timsPGC(istart)-timsPGC(j))>=2
                    %alltn(j,n)
%                     if alltn(j,n)>=0.8
%                         j
%                         timsPGC(istart)
%                         timsPGC(j)
%                         pause
%                     end
                    f05=f05+1;
                    n05=n05+1; j05=j05+1;
                    tempfindt(f05,:)=timsPGC(istart:iend);
                    targfindt(f05,:)=timsPGC(j:j+templen-1);
                    tempoffSSIB(f05)=offSSIB(istart);
                    tempoffSILB(f05)=offSILB(istart);
                    targoffSSIB(f05)=offSSIB(j);
                    targoffSILB(f05)=offSILB(j);
                    if offSSIB(istart)~=offSSIB(iend)
                        n
                        j
                        1+0
                    end
                    if offSILB(istart)~=offSILB(iend)
                        n
                        j
                        1+0
                    end
                    if offSSIB(j)~=offSSIB(j+templen-1)
                        n
                        j
                        1+2
                    end
                    if offSILB(j)~=offSILB(j+templen-1)
                        n
                        j
                        1+3
                    end
                    tempfindPGC(:,f05)=PGCtempl;
                    tempfindSSIB(:,f05)=SSIBtempl;
                    tempfindSILB(:,f05)=SILBtempl;
                    targfindPGC(:,f05)=realPGC(j:j+templen-1);
                    targfindSSIB(:,f05)=realSSIB(j:j+templen-1);
                    targfindSILB(:,f05)=realSILB(j:j+templen-1);
                    find05(f05,:)=[n,j];
                end
            end %In case data string is not monotonically increasing
            end
        end
    end
end

%lump things
find05'
out=zeros(f05,1); 
keeptot=0;
for i=1:f05
    if out(i) < 0.5 %if not previously kicked out (out=0 instead of 1)
        clear keeptmp
        clear keeptmpx
        keeptot=keeptot+1; %how many will be kept in the end
        inset=1; %how many in the set to compare; at least this one
        keeptmp(inset)=i; 
        keeptmpx(inset)=alltn(find05(i,2),find05(i,1)); %alltn(n,j) of the i_th member of f05
        keeptmpPGCtn(inset)=PGCtn(find05(i,2),find05(i,1));
        keeptmpSSIBtn(inset)=SSIBtn(find05(i,2),find05(i,1));
        keeptmpSILBtn(inset)=SILBtn(find05(i,2),find05(i,1));
        keeptmpnj(inset,:)=[find05(i,1), find05(i,2)];
        for k=i+1:f05
            %if isequal(find05(i,1),find05(k,1)) && find05(k,2)-find05(i,2)<=80 %Same template, targets closer than 80 samples
            %if isequal(find05(i,1),find05(k,1)) && abs(targfindt(i,1)-targfindt(k,1))<=2 %Same temp, targs closer than 2s
            if abs(tempfindt(i,1)-tempfindt(k,1))<=1 && abs(targfindt(i,1)-targfindt(k,1))<=2 %temps & targs closer than 2s
                inset=inset+1;
                keeptmp(inset)=k; out(k)=1;
                keeptmpx(inset)=alltn(find05(k,2),find05(k,1));
                keeptmpPGCtn(inset)=PGCtn(find05(k,2),find05(k,1));
                keeptmpSSIBtn(inset)=SSIBtn(find05(k,2),find05(k,1));
                keeptmpSILBtn(inset)=SILBtn(find05(k,2),find05(k,1));
                keeptmpnj(inset,:)=[find05(k,1), find05(k,2)];
            elseif find05(k,1)-find05(i,1)<=85.0/tempoff && abs(find05(k,2)-find05(i,2))<=85 %Templates within a few (4); targets?
                inset=inset+1;
                keeptmp(inset)=k; out(k)=1;
                keeptmpx(inset)=alltn(find05(k,2),find05(k,1));
                keeptmpPGCtn(inset)=PGCtn(find05(k,2),find05(k,1));
                keeptmpSSIBtn(inset)=SSIBtn(find05(k,2),find05(k,1));
                keeptmpSILBtn(inset)=SILBtn(find05(k,2),find05(k,1));
                keeptmpnj(inset,:)=[find05(k,1), find05(k,2)];
           end
        end
        keeptmp
        keeptmpx
        [maxkeep(keeptot),ikeep]=max(keeptmpx); %keep track of the kept maxima
        ikeep %ikeep varies only from 1 to inset, the number in the set being compared
        keeptmp(ikeep) %keeptmp(ikeep) varies between 1 and f05
        tempkeept(keeptot,:)=tempfindt(keeptmp(ikeep),:);
        targkeept(keeptot,:)=targfindt(keeptmp(ikeep),:);
        tempkeepPGC(:,keeptot)=tempfindPGC(:,keeptmp(ikeep));
        tempkeepSSIB(:,keeptot)=tempfindSSIB(:,keeptmp(ikeep));
        tempkeepSILB(:,keeptot)=tempfindSILB(:,keeptmp(ikeep));
        targkeepPGC(:,keeptot)=targfindPGC(:,keeptmp(ikeep));
        targkeepSSIB(:,keeptot)=targfindSSIB(:,keeptmp(ikeep));
        targkeepSILB(:,keeptot)=targfindSILB(:,keeptmp(ikeep));
        tempkeepoffSSIB(keeptot)=tempoffSSIB(keeptmp(ikeep));
        tempkeepoffSILB(keeptot)=tempoffSILB(keeptmp(ikeep));
        targkeepoffSSIB(keeptot)=targoffSSIB(keeptmp(ikeep));
        targkeepoffSILB(keeptot)=targoffSILB(keeptmp(ikeep));
        
        keepPGCtn(keeptot)=keeptmpPGCtn(ikeep);
        keepSSIBtn(keeptot)=keeptmpSSIBtn(ikeep);
        keepSILBtn(keeptot)=keeptmpSILBtn(ikeep);
        keepnj(keeptot,:)=keeptmpnj(ikeep,:);
    end
end
keeptot
[maxkeepsort,mksi] = sort(maxkeep,'descend');
ncheck=floor((lenx-templen)/templen); %looks for the max in a string of pt-wise cc m'ments
histmat=zeros(ncheck,ntemp);
for n=1:ncheck
    istart=(n-1)*templen+1;
    iend=n*templen;
    histmat(n,:)=max(alltn(istart:iend,:));
end
nlump=floor(ntemp/4);
hist2mat=zeros(ncheck,nlump);
for n=1:nlump
    istart=(n-1)*4+1;
    iend=n*4;
    hist2mat(:,n)=max(histmat(:,istart:iend),[],2);
end
[m,n]=size(hist2mat);
hist2=reshape(hist2mat,1,m*n);
[look,x]=hist(hist2,50);
look2=[x;look];
fid = fopen('look.txt','w');
fprintf(fid,'%8.4f  %8i\n',look2);
fclose(fid);
        
        % orient landscape
        % print('-dpdf',[JDAY,'_',int2str(timlook),'_4.pdf'])tempkeept
        
ltmax=max([realPGC; realSSIB; realSILB]);
ltmin=min([realPGC; realSSIB; realSILB]);
slope=(ltmax-ltmin)/160;
for ifig=1:floor(keeptot/8)
    figure
    for isub=1:8
        icount=(ifig-1)*8+isub;
        crossoff=(2*ltmax+(0:159));
        for ij=1:icount-1  %Get rid of "duplicates" (same small offsets for template and target)
            if abs((tempkeept(mksi(icount),1)-tempkeept(mksi(ij),1)) ...
                  -(targkeept(mksi(icount),1)-targkeept(mksi(ij),1))) <=0.1 && ...
               abs((tempkeept(mksi(icount),1)-tempkeept(mksi(ij),1))) <= 4.0
%                     targkeepPGC(:,mksi(icount))=0.;
%                     targkeepSSIB(:,mksi(icount))=0.;
%                     targkeepSILB(:,mksi(icount))=0.;
                    crossoff=(ltmin+slope*(0:159));
            end
        end
        subplot(4,6,(isub-1)*3+1,'align'); 
        hold on
        plot(targkeepPGC(:,mksi(icount)),'r');
        plot(tempkeepPGC(:,mksi(icount)),'g');
        plot(crossoff,'k--')
        text(10, 0.8*ltmax, int2str(icount));
        text(100, 0.8*ltmax, num2str(maxkeep(mksi(icount))),'fontsize',7);
        text(110, 0.8*ltmin, num2str(keepPGCtn(mksi(icount))),'fontsize',7);
        text(10, 0.8*ltmin, 'PGC','fontsize',7);
        text(60, 0.6*ltmax, int2str(tempkeept(mksi(icount),floor(templen/2))),'fontsize',7); %This writes the central time of window
        text(110, 0.6*ltmax, int2str(targkeept(mksi(icount),floor(templen/2))),'fontsize',7);
        axis([0 160 ltmin ltmax]);
        set(gca,'fontsize',8);
        box on
        subplot(4,6,(isub-1)*3+2,'align'); 
        hold on
        plot(targkeepSSIB(:,mksi(icount)),'b');
        plot(tempkeepSSIB(:,mksi(icount)),'g');
        plot(crossoff,'k--')
        text(10, 0.8*ltmax, int2str(icount));
        text(100, 0.8*ltmax, num2str(maxkeep(mksi(icount))),'fontsize',7);
        text(110, 0.8*ltmin, num2str(keepSSIBtn(mksi(icount))),'fontsize',7);
        text(10, 0.8*ltmin, 'SSIB','fontsize',7);
        text(60, 0.8*ltmin, int2str(tempkeepoffSSIB(mksi(icount))),'fontsize',7);
        text(80, 0.8*ltmin, int2str(targkeepoffSSIB(mksi(icount))),'fontsize',7);
        text(60, 0.6*ltmax, int2str(tempkeept(mksi(icount),floor(templen/2))),'fontsize',7);
        text(110, 0.6*ltmax, int2str(targkeept(mksi(icount),floor(templen/2))),'fontsize',7);
        %text(100, 0.6*ltmax, int2str(keepnj(mksi(icount),:)));
        axis([0 160 ltmin ltmax]);
        set(gca,'fontsize',8);
        box on
        subplot(4,6,(isub-1)*3+3,'align'); 
        hold on
        plot(targkeepSILB(:,mksi(icount)),'k');
        plot(tempkeepSILB(:,mksi(icount)),'g');
        plot(crossoff,'k--')
        text(10, 0.8*ltmax, int2str(icount));
        text(100, 0.8*ltmax, num2str(maxkeep(mksi(icount))),'fontsize',7);
        text(110, 0.8*ltmin, num2str(keepSILBtn(mksi(icount))),'fontsize',7);
        text(10, 0.8*ltmin, 'SILB','fontsize',7);
        text(60, 0.8*ltmin, int2str(tempkeepoffSILB(mksi(icount))),'fontsize',7);
        text(80, 0.8*ltmin, int2str(targkeepoffSILB(mksi(icount))),'fontsize',7);
        text(60, 0.6*ltmax, int2str(tempkeept(mksi(icount),floor(templen/2))),'fontsize',7);
        text(110, 0.6*ltmax, int2str(targkeept(mksi(icount),floor(templen/2))),'fontsize',7);
        %text(100, 0.6*ltmax, int2str(keepnj(mksi(icount),:)));
        axis([0 160 ltmin ltmax]);
        set(gca,'fontsize',8);
        box on
    end
    orient landscape
    print('-dpdf',['255_',int2str(timstart),'_BrownS',int2str(ifig),'.pdf'])
end

% 
ltmax=max([datPGC(lookstart:lookend); datSSIB(lookstart:lookend); datSILB(lookstart:lookend)]);
ltmin=min([datPGC(lookstart:lookend); datSSIB(lookstart:lookend); datSILB(lookstart:lookend)]);
% timbig=4.33333333334*60; 
% winbig=2*timbig*40;
figure 
subplot(4,1,1); 
hold on
plot(timsPGC,realPGC,'r');
plot(timsSSIB,realSSIB,'b');
plot(timsSILB,realSILB,'k');
axis([timstart timstart+delt/4 ltmin ltmax]);
box on
subplot(4,1,2); 
hold on
plot(timsPGC,realPGC,'r');
plot(timsSSIB,realSSIB,'b');
plot(timsSILB,realSILB,'k');
axis([timstart+delt/4 timstart+delt/2 ltmin ltmax]);
box on
subplot(4,1,3); 
hold on
plot(timsPGC,realPGC,'r');
plot(timsSSIB,realSSIB,'b');
plot(timsSILB,realSILB,'k');
axis([timstart+delt/2 timend-delt/4 ltmin ltmax]);
box on
subplot(4,1,4); 
hold on
plot(timsPGC,realPGC,'r');
plot(timsSSIB,realSSIB,'b');
plot(timsSILB,realSILB,'k');
axis([timend-delt/4 timend ltmin ltmax]);
box on
% 
%title([IDENTIF,'_{',int2str(lo),'-',int2str(hi),'}'])
orient landscape
print('-depsc','tmp.eps')

cputime-t;
tot=cputime-tt

% % 
% % nin4=length(locs4);
% % nin128=length(locs128);
% % for k=1:ntees
% %     for i=2:nin4-1
% %         if locs4(i-1,1) < tees(k,1)-1 && locs4(i,1) >= tees(k,1)-1
% %             istart4=i;
% %         end
% %         if locs4(i+1,1) > tees(k,2)+1 && locs4(i,1) <= tees(k,2)+1
% %             iend4=i;
% %         end
% %     end
% %     for i=2:nin128-1
% %         if locs128(i-1,1) < tees(k,1)-1 && locs128(i,1) >= tees(k,1)-1
% %             istart128=i;
% %         end
% %         if locs128(i+1,1) > tees(k,2)+1 && locs128(i,1) <= tees(k,2)+1
% %             iend128=i;
% %         end
% %     end
% % 
% %     figure
% %     colormap(jet) %PGSI=x, PGSS=y
% %     plot(0.33*locs128(1:istart128,2),0.6*locs128(1:istart128,3),'ko','MarkerSize',2,'linewidth',0.2)
% %     hold on 
% %     if OLDER==1
% %         plot(0.33*olderlocs(:,2),0.6*olderlocs(:,3),'ko','MarkerSize',2,'linewidth',0.2)
% %     end
% %     scatter(0.33*locs128(istart128:iend128,2),0.6*locs128(istart128:iend128,3),90,locs128(istart128:iend128,1)...
% %         ,'linewidth',1)
% %     scatter(0.33*locs4(istart4:iend4,2),0.6*locs4(istart4:iend4,3),40,locs4(istart4:iend4,1),'filled')
% %     axis([-5 5 -5 5])
% %     xlabel('PGSI - KM ESE')
% %     ylabel('PGSS - KM NNE')
% %     caxis([tees(k,1) tees(k,2)])
% %     colorbar
% %     title(IDENTIF4)
% %     box on
% %     print('-depsc',['tmp',int2str(k),'.eps'])
% % box on
% % end
% 
