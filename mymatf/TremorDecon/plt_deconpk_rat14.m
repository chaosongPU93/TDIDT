function f=plt_deconpk_rat14(f,impindepstall,srcamprall,color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_deconpk_rat14(f,srcamprall,color)
%
% This function is to specifically plot the relationship between the amp at
% sta 1 and 4, and their ratio. This plot aims to lump all found sources 
% from all burst windows as contribution to the plot.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/11/17
% Last modified date:   2022/11/17 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trust4th = 7;
ax = f.ax(1);
hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
scatter(ax,log10(impindepstall(:,2)),log10(impindepstall(:,9+(trust4th-3)*2)),20,color);
axis(ax,'equal');
%       axis(ax,[0 1.2 0 1.2]);
plot(ax,[-1.5 1.5],[-1.5 1.5],'--','linew',1.5,'color',[.3 .3 .3]);
xlabel(ax,'log_{10}{Amp at 1}');
ylabel(ax,'log_{10}{Amp at 4}');

ax = f.ax(2);
hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
histogram(ax,log10(srcamprall(:,end)),'Normalization','pdf',...
  'BinWidth',0.05,'FaceColor',color);
plot(ax,[median(log10(srcamprall(:,end))) median(log10(srcamprall(:,end)))],ax.YLim,'--',...
  'color',color,'linew',1);
errorbar(ax,median(log10(srcamprall(:,end))), 0.9*ax.YLim(2), mad(log10(srcamprall(:,end)),1),...
  mad(log10(srcamprall(:,end)),1),'horizontal','o',...
  'color',color,'linewidth',0.8,'CapSize',5,'MarkerSize',0.5);
% quant = log10(srcamprall(:,end));
% [muHat,sigmaHat] = normfit(quant);
% pdfhat = normpdf(-1:0.02:1,muHat,sigmaHat);
% plot(ax,-1:0.02:1,pdfhat,'r-','linew',2);
% plot(ax,[muHat muHat],ax.YLim,'r--','linew',1);
% plot(ax,[median(quant) median(quant)],ax.YLim,'b--','linew',1);
% text(ax,0.05,0.9,sprintf('1*sigma=%.2f',sigmaHat),'Units','normalized');
xlabel(ax,'log_{10}{Amp ratio 1/4}');
ylabel(ax,'PDF');

ax = f.ax(3);
hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
scatter(ax,log10(impindepstall(:,2)),log10(srcamprall(:,end)),20,color);
xlabel(ax,'log_{10}{Amp at 1}');
ylabel(ax,'log_{10}{Amp ratio 1/4}');
axis(ax,'equal');

%   ax = f.ax(3);
%   hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%   scatter(ax,impindepstall(:,2),srcamprall(:,end),20,'k');
%   xlabel(ax,'Amp ratio at 1');
%   ylabel(ax,'Amp ratio 1/4');
%   axis(ax,'equal');


