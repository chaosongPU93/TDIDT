function PlotIndividual(fam,date,sttime,edtime,hfall,lfall)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to plot 12 & 13 offset in lf and hf of their common
% detections during the individual suspicous migrations, including 2 
% subplots, one is the respective offset, one is the difference.  
% 
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/06/21
% Last modified date:   2019/06/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fam=num2str(fam);
datestr = num2str(date);
yr = datestr(1:4);
jday = datestr(5:end);
hfday = hfall(hfall(:,5)==date, :);
lfday = lfall(lfall(:,5)==date, :);

spsratio=2;
msize = 4;
scrsz=get(0,'ScreenSize');
figure('Position',[scrsz(3)/5 scrsz(4)/5 2*scrsz(3)/5 2*scrsz(4)/5]);

subplot(2,1,1,'align');  % plot the original migration, hf and lf respectively
hold on
p1=plot(lfday(:,end),hfday(:,1),'bo','MarkerSize',msize);  % 12 off hf
p2=plot(lfday(:,end),hfday(:,2),'bd','MarkerSize',msize);  % 13 off hf
p3=plot(lfday(:,end),lfday(:,1)*spsratio,'ro','MarkerSize',msize);  % 12 off lf
p4=plot(lfday(:,end),lfday(:,2)*spsratio,'rd','MarkerSize',msize);  % 13 off lf
legend([p1,p2,p3,p4],{'ori hf off12','ori hf off13','ori lf off12','ori lf off13'},'location','best','Fontsize',10);
set(gca, 'xlim',[sttime edtime]);
set(gca, 'ylim',[-30 30]);
ylabel('Samples', 'fontsize',12);
title(strcat(fam,'.',yr,'.',jday,'\_',num2str(sttime),'-',num2str(edtime),'sec'), 'fontsize', 14);
box on

subplot(2,1,2,'align');
p1=plot(lfday(:,end),hfday(:,1)-lfday(:,1)*spsratio,'bo','MarkerSize',msize); hold on  % diff 12 off
p2=plot(lfday(:,end),hfday(:,2)-lfday(:,2)*spsratio,'rd','MarkerSize',msize); hold on  % diff 13 off
plot([sttime edtime],[0 0], 'k--','linewidth',1); hold on
legend([p1,p2],{'ori (hf-lf)12','ori (hf-lf)13'},'location','best','Fontsize',10);
set(gca, 'xlim',[sttime edtime]);
set(gca, 'ylim',[-10 10]);
ylabel('Samples', 'fontsize',12);
box on