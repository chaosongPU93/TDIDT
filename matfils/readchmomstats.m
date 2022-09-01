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
stats(:,9:10)=log10(stats(:,9:10));
%sorstats=sortrows(stats,9);
%ntot=size(sorstats,1)

%inbatch=floor(ntot/20);
%for i=1:20
%    h=figure('Position',[wid/3 1 0.6*wid 0.6*hite]); 
%    istart=(i-1)*inbatch+1;
%    iend=i*inbatch;
%    scatter(log10(sorstats(istart:iend,3)),log10(sorstats(istart:iend,6)),12,sorstats(istart:iend,8),'filled')
%end

catlg=stats;
ndetects=size(catlg,1)
range=20; 
imin=3;

%Deal with time of day
for i=1:ndetects-1
    if catlg(i+1,2)<catlg(i,2) 
        catlg(i+1:end,2)=catlg(i+1:end,2)+100000;
    end
end

%Label groups
group=ceil(catlg(:,2)/range);

i=1;
ikept=0;
while i<ndetects
    tmp(1,:)=catlg(i,:);
    j=i+1;
    while group(j) == group(i) && j<ndetects
        tmp(j-i+1,:)=catlg(j,:);
        j=j+1;
    end
    ingroup=j-i;
    if size(tmp,1)>=imin
        y(ikept+1:ikept+ingroup)=log10(tmp(:,6));
        x(ikept+1:ikept+ingroup)=log10(tmp(:,3));
        x(ikept+1:ikept+ingroup)=x(ikept+1:ikept+ingroup)-mean(x(ikept+1:ikept+ingroup));
        catkept(ikept+1:ikept+ingroup,:)=tmp;
        ikept=ikept+ingroup;
    end
    i=j;
    clear tmp
end

%Now for the half-overlapping windows
catlg(:,2)=catlg(:,2)+range/2;
group=ceil(catlg(:,2)/range);

i=1;
ikept=ikept;
while i<ndetects
    tmp(1,:)=catlg(i,:);
    j=i+1;
    while group(j) == group(i) && j<ndetects
        tmp(j-i+1,:)=catlg(j,:);
        j=j+1;
    end
    ingroup=j-i;
    if size(tmp,1)>=imin
        y(ikept+1:ikept+ingroup)=log10(tmp(:,6));
        x(ikept+1:ikept+ingroup)=log10(tmp(:,3));
        x(ikept+1:ikept+ingroup)=x(ikept+1:ikept+ingroup)-mean(x(ikept+1:ikept+ingroup));
        catkept(ikept+1:ikept+ingroup,:)=tmp;
        ikept=ikept+ingroup;
    end
    i=j;
    clear tmp
end
ikept

%w(:,1)=x;
%w(:,2)=y;
%nbatch=9;
%ingroup=floor(ikept/nbatch)
%leftover=ikept-nbatch*ingroup
%fil=sortrows(w,1);
%%fid = fopen([filnam(1:50),'pre'],'w');
%%for i=1:ndetects
%% fprintf(fid,'%9.1f %6.2f %6.2f %8.3f %7.2f %8.3f %8.3f %10.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',fil(i,1:19)); %transpose?
%%end
%%fclose(fid);

%%medians=zeros(nbatch,2);
%for ibat=1:nbatch-1
% medians(ibat,:)=median(fil((ibat-1)*ingroup+1:ibat*ingroup,:));
%end
%ibat=nbatch;
%medians(ibat,:)=median(fil((ibat-1)*ingroup+1:ikept,:));
%%medians(ibat+1,:)=median(fil);
%%medians(ibat+1,19)=9999;
%fid = fopen('w.medians','w');
%for i=1:nbatch
% fprintf(fid,'%8.3f %8.3f \n',medians(i,1:2)); %transpose?
%end
%fclose(fid);

h=figure('Position',[0.1*wid 0.25*hite wid/2.5 hite/2]); %center
scatter(x,y,30,catkept(:,8),'filled')
%plot(x,y,'ro','MarkerSize',2,'MarkerFaceColor','r')
%hold on
%plot(medians(:,1),medians(:,2),'ko','MarkerSize',6,'MarkerFaceColor','k')
%xlim([-1.5 1.5])
%title(['range ',int2str(range),'s;  min #events ',int2str(imin)]);
%text(-1.2,-0.5,int2str(ikept),'fontsize',8);
%axis equal
%box on
%print(h,'-depsc',['local_range',int2str(range),'_min',int2str(imin),'.eps'])
 
% xlabel('amplitude of velocity peak')
% ylabel('x-corr coefficient')
% set(h,'PaperPosition',[0.25 5 4 3])
% print(h,'-depsc','BFIGS/datacc.eps')
% 
% set(h,'PaperPosition',[0.25 0.25 8 10.5])
% orient landscape
% print(h,'-depsc',['BFIGS/',num2str(lo),'-',num2str(hi),'_',int2str(npa),'pass_',int2str(0),int2str(0),int2str(ifig),'.eps'])
