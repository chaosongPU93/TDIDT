% stats(nttot+nt,:)=[year+0.001*jday bostsec(nt) bostmom(nttot+nt) maxSISS(nttot+nt) ...
%                        5   0.5*(maxccSItemp(nttot+nt)+maxccSStemp(nttot+nt)) ...
%                        6   0.5*((imaxSIstretch(nttot+nt)+19)/40+(imaxSSstretch(nttot+nt)+19)/40) ...
%                        7   0.5*(scaleSIstore(nttot+nt)+scaleSSstore(nttot+nt)) ...
%                        8   0.5*maxSISS(nttot+nt)*(maxccSItemp(nttot+nt)+maxccSStemp(nttot+nt)) ...
%                        9   0.5*(SILBamp15(nttot+nt)+SSIBamp15(nttot+nt)) ...
%                       10   0.5*(SILBamp30(nttot+nt)+SSIBamp30(nttot+nt))];
close all
clear all
scrsz=get(0,'ScreenSize');
wid=scrsz(3);
hite=scrsz(4);
scrat=wid/hite;

load('checkmomstats.mat')
%stats(:,9:10)=log10(stats(:,9:10));
ndetects=size(stats,1);
ncol=size(stats,2);

h=figure('Position',[wid/3 1 0.4*wid hite]);

subplot(6,2,1,'align')
plot(log10(stats(:,3)),log10(stats(:,6)),'ro','MarkerSize',2)
%scatter(log10(stats(:,3)),log10(stats(:,6)),12,stats(:,8))
hold on
tmp=sortrows(stats,3);
nbatch=16
inbatch=floor(ndetects/nbatch)
leftover=ndetects-nbatch*inbatch
medians=zeros(nbatch,ncol);
for ibat=1:nbatch-1
    medians(ibat,:)=median(tmp((ibat-1)*inbatch+1:ibat*inbatch,:));
end
ibat=nbatch;
medians(ibat,:)=median(tmp((ibat-1)*inbatch+1:ndetects,:));
plot(log10(medians(:,3)),log10(medians(:,6)),'ko','markerfacecolor','k','markersize',4)
xlabel('log_{10} Bostock moment')
ylabel('log_{10}(average time stretch)')
title('002 data, SSIB & SILB only, fat templates')
%caxis([0 0.7])
%hcb=colorbar;
%ylabel(hcb,'Ave. CC, template & STAs')
axis equal
xlim([10.8 12.6])
%ylim([-0.5 2])
box on

sorstats=sortrows(stats,3);
mstats=sorstats;
TF1 = mstats(:,3) < 1.15e11;
TF2 = mstats(:,3) > 9e11;
TFall = TF1 | TF2;
mstats(TFall,:)=[];
mdetects=size(mstats,1);
inbatch1=50 %163
nbatch=floor(mdetects/inbatch1)
leftover=mdetects-nbatch*inbatch1
subplot(6,2,2,'align')
hold on
for i=1:nbatch-1
    istart=(i-1)*inbatch1+1;
    iend=i*inbatch1;
    %scatter(log10(sorstats(istart:iend,3)),log10(sorstats(istart:iend,6)),12,sorstats(istart:iend,8),'filled')
    mstats(istart:iend,ncol+1)=log10(mstats(istart:iend,7))-median(log10(mstats(istart:iend,7)));
    plot(mstats(istart:iend,ncol+1),log10(mstats(istart:iend,6).^2),'ro','MarkerSize',2)
   % scatter(mstats(istart:iend,ncol+1),log10(mstats(istart:iend,6)),12,mstats(istart:iend,8))
end
i=nbatch;
istart=(i-1)*inbatch1+1;
iend=mdetects;
mstats(istart:iend,ncol+1)=log10(mstats(istart:iend,7))-median(log10(mstats(istart:iend,7)));
plot(mstats(istart:iend,ncol+1),log10(mstats(istart:iend,6).^2),'ro','MarkerSize',2)

tmp=sortrows(mstats,ncol+1);
nbatch=16
inbatch=floor(mdetects/nbatch)
leftover=mdetects-nbatch*inbatch
medians=zeros(nbatch,ncol+1);
for ibat=1:nbatch-1
    medians(ibat,:)=median(tmp((ibat-1)*inbatch+1:ibat*inbatch,:));
end
ibat=nbatch;
medians(ibat,:)=median(tmp((ibat-1)*inbatch+1:mdetects,:));
plot(medians(:,ncol+1),log10(medians(:,6).^2),'ko','markerfacecolor','k','markersize',4)
xlabel('log_{10} (av. amplitude stretch) - median in bin')
ylabel('log_{10}(average time stretch^2)')
title('002 data, binned by Bostock moment')
%caxis([0 0.7])
%hcb=colorbar;
%ylabel(hcb,'Ave. CC, template & STAs')
axis equal
xlim([-1.2 1.2])
%ylim([-0.5 2])
box on

sorstats=sortrows(stats,9);
inbatch1=50 %163
nbatch=floor(ndetects/inbatch1)
leftover=ndetects-nbatch*inbatch1
subplot(3,1,2,'align')
hold on
for i=1:nbatch-1
    istart=(i-1)*inbatch1+1;
    iend=i*inbatch1;
    %scatter(log10(sorstats(istart:iend,3)),log10(sorstats(istart:iend,6)),12,sorstats(istart:iend,8),'filled')
    sorstats(istart:iend,ncol+1)=log10(sorstats(istart:iend,3))-median(log10(sorstats(istart:iend,3)));
    plot(sorstats(istart:iend,ncol+1),log10(sorstats(istart:iend,6)),'ro','MarkerSize',2)
   % scatter(sorstats(istart:iend,ncol+1),log10(sorstats(istart:iend,6)),12,sorstats(istart:iend,8))
end
i=nbatch;
istart=(i-1)*inbatch1+1;
iend=ndetects;
sorstats(istart:iend,ncol+1)=log10(sorstats(istart:iend,3))-median(log10(sorstats(istart:iend,3)));
plot(sorstats(istart:iend,ncol+1),log10(sorstats(istart:iend,6)),'ro','MarkerSize',2)

tmp=sortrows(sorstats,ncol+1);
nbatch=16
inbatch=floor(ndetects/nbatch)
leftover=ndetects-nbatch*inbatch
medians=zeros(nbatch,ncol+1);
for ibat=1:nbatch-1
    medians(ibat,:)=median(tmp((ibat-1)*inbatch+1:ibat*inbatch,:));
end
ibat=nbatch;
medians(ibat,:)=median(tmp((ibat-1)*inbatch+1:ndetects,:));
plot(medians(:,ncol+1),log10(medians(:,6)),'ko','markerfacecolor','k','markersize',4)
xlabel('log_{10} (Bostock moment) - median in that bin')
ylabel('log_{10}(average time stretch)')
title('002 data, binned by amp^2 in prior 15 s')
%caxis([0 0.7])
%hcb=colorbar;
%ylabel(hcb,'Ave. CC, template & STAs')
axis equal
xlim([-1 1])
%ylim([-0.5 2])
box on

subplot(3,1,3,'align')
hold on
for i=1:nbatch-1
    istart=(i-1)*inbatch1+1;
    iend=i*inbatch1;
    %scatter(log10(sorstats(istart:iend,3)),log10(sorstats(istart:iend,6)),12,sorstats(istart:iend,8),'filled')
    sorstats(istart:iend,ncol+1)=log10(sorstats(istart:iend,7).*sorstats(istart:iend,6).^2)- ...
        median(log10(sorstats(istart:iend,7).*sorstats(istart:iend,6).^2));
    plot(sorstats(istart:iend,ncol+1),log10(sorstats(istart:iend,6)),'ro','MarkerSize',2)
    %scatter(sorstats(istart:iend,ncol+1),log10(sorstats(istart:iend,6)),12,sorstats(istart:iend,8))
end
i=nbatch;
istart=(i-1)*inbatch1+1;
iend=ndetects;
sorstats(istart:iend,ncol+1)=log10(sorstats(istart:iend,7).*sorstats(istart:iend,6).^2)- ...
    median(log10(sorstats(istart:iend,7).*sorstats(istart:iend,6).^2));
plot(sorstats(istart:iend,ncol+1),log10(sorstats(istart:iend,6)),'ro','MarkerSize',2)

tmp=sortrows(sorstats,ncol+1);
nbatch=16
inbatch=floor(ndetects/nbatch)
leftover=ndetects-nbatch*inbatch
medians=zeros(nbatch,ncol+1);
for ibat=1:nbatch-1
    medians(ibat,:)=median(tmp((ibat-1)*inbatch+1:ibat*inbatch,:));
end
ibat=nbatch;
medians(ibat,:)=median(tmp((ibat-1)*inbatch+1:ndetects,:));
plot(medians(:,ncol+1),log10(medians(:,6)),'ko','markerfacecolor','k','markersize',4)
xlabel('log_{10} average(amp stretch x time stretch^2) - median in that bin')
ylabel('log_{10}(average time stretch)')
title('002 data, binned by amp^2 in prior 15 s')
%caxis([0 0.7])
%hcb=colorbar;
%ylabel(hcb,'Ave. CC, template & STAs')
axis equal
xlim([-1 1])
%ylim([-0.5 2])
box on

set(h,'PaperPosition',[0.25 0.25 8 10])
print(h,'-depsc','localBostock.eps')

