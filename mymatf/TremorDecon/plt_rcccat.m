function [f] = plt_rcccat(rccpair,sps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_rcccat(rccpaircat)
%
% Plot the running cross-correlation (RCC) of different station pairs and  
% the average, and the empirical CDF of them.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/07/07
% Last modified date:   2022/07/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

lrcc = size(rccpair,1);

f.fig = figure;
subplot(321)
plot((1:lrcc)/sps,rccpair(:,1),'r','linew',1); hold on ;
plot((1:lrcc)/sps,rccpair(:,2),'b','linew',1);
plot((1:lrcc)/sps,rccpair(:,3),'k','linew',1);
% xlabel('Samples');
xlabel('Time (s) since the start');
ylabel('Running CC');

subplot(322)
plot((1:lrcc)/sps,mean(rccpair,2),'Color',[.5 .5 .5],'linew',1); hold on ;
xlabel('Time (s) since the start');
ylabel('Running CC');

subplot(323)
[cdfval,x] = ecdf(rccpair(:,1));
plot(x,cdfval,'r','linew',1); hold on ;
[cdfval,x] = ecdf(rccpair(:,2));
plot(x,cdfval,'b','linew',1);
[cdfval,x] = ecdf(rccpair(:,3));
plot(x,cdfval,'k','linew',1);
text(0.95,0.5,sprintf('med=%.2f',median(rccpair(:,1))),'HorizontalAlignment','right','Color','r');
text(0.95,0.4,sprintf('med=%.2f',median(rccpair(:,2))),'HorizontalAlignment','right','Color','b');
text(0.95,0.3,sprintf('med=%.2f',median(rccpair(:,3))),'HorizontalAlignment','right','Color','k');
ylabel('Empirical CDF');
xlabel('Running CC');
legend('sta1-2','sta1-3','sta2-3','Location','northwest');

subplot(324)
[cdfval,x] = ecdf(mean(rccpair,2));
plot(x,cdfval,'Color',[.5 .5 .5],'linew',1); hold on ;
text(0.95,0.5,sprintf('med=%.2f',median(mean(rccpair,2))),'HorizontalAlignment','right');
ylabel('Empirical CDF');
xlabel('Running CC');
legend('average','Location','northwest');

subplot(325)
[pdfval,x] = ksdensity(rccpair(:,1),-1:0.01:1);
% h=histogram(rccpair(:,1),'binw',0.01,'Normalization','pdf');
plot(x,pdfval,'r','linew',1); hold on ;
[pdfval,x] = ksdensity(rccpair(:,2),-1:0.01:1);
plot(x,pdfval,'b','linew',1); 
[pdfval,x] = ksdensity(rccpair(:,3),-1:0.01:1);
plot(x,pdfval,'k','linew',1);
ylabel('Empirical PDF');
xlabel('Running CC');

subplot(326)
[pdfval,x] = ksdensity(mean(rccpair,2),-1:0.01:1);
plot(x,pdfval,'Color',[.5 .5 .5],'linew',1); hold on ;
ylabel('Empirical PDF');
xlabel('Running CC');

