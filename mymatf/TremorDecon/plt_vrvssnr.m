function f=plt_vrvssnr(greentype,mfit,varsig,mfitort,varsigort,...
          mfitvert,varsigvert,stas,ista,snrf,snrfort,snrfvert,...
          varred,varredort,varredvert,varred4th,varredort4th,varredvert4th,...
          detecttype,saveflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is a simple function to plot the variance reduction vs the 
% Signal-to-noise ratio of templates if using the grouped sources or 
% 4-sta sources. 
% If the direct-stack templates are used as the Green's functions, then
% 1 extra panel of relative variance change with iteration is added.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/07/16
% Last modified date:   2024/07/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if strcmp(greentype, 'dstack')
  
  nit = length(mfit);
  relavarort = mfitort/varsigort*100;
  relavarvert = mfitvert/varsigvert*100;
  relavar = mfit/varsig*100;

  widin = 7.8;  % maximum width allowed is 8.5 inches
  htin = 3;   % maximum height allowed is 11 inches
  nrow = 1;
  ncol = 3;
  f = initfig(widin,htin,nrow,ncol);
  pltxran = [0.08 0.98]; pltyran = [0.15 0.98]; % optimal axis location
  pltxsep = 0.07; pltysep = 0.05;
  optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
  
  ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
  p(1)=plot(ax,1:nit,relavar,'r-','linew',1);
  p(2)=plot(ax,1:nit,relavarort,'b-','linew',1);
  p(3)=plot(ax,1:nit,relavarvert,'k-','linew',1);
  ylabel(ax,'Relative variance of residual (%)');
  xlabel(ax,'Iteration #');
  lgd=legend(ax,p,'Optimal','Orthogonal','Vertical','Location','southwest');
  set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
  title(lgd,strtrim(stas(ista,:)),'fontsize',9,'FontWeight','normal');
%   text(ax,0.99,0.3,strtrim(stas(ista,:)),...
%     'Units','normalized','HorizontalAlignment','right','fontsize',10);
  text(ax,0.98,0.94,'a','FontSize',10,'unit','normalized','HorizontalAlignment',...
    'right','EdgeColor','k','Margin',1,'backgroundcolor','w');
  ylim(ax,[0 100]);

  ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
  pltsta=[1 2 3];
  p=[];
  p(1)=scatter(ax,snrf(pltsta),varred(pltsta,3),20,'ro');
  p(2)=scatter(ax,snrfort(pltsta),varredort(pltsta,3),24,'bs');
  p(3)=scatter(ax,snrfvert(pltsta),varredvert(pltsta,3),20,'k^');
  scatter(ax,snrf(1),varred(1,3),20,'ro','filled');
  scatter(ax,snrfort(1),varredort(1,3),24,'bs','filled');
  scatter(ax,snrfvert(1),varredvert(1,3),20,'k^','filled');
  lgd=legend(ax,p,'Optimal','Orthogonal','Vertical','location','southeast');
  set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
  xlabel(ax,'LFE template SNR');
  ylabel(ax,'Variance reduction (%)');
  %       title(ax,'Right after grouping');
  text(ax,0.01,0.94,'Grouped sources',...
    'Units','normalized','HorizontalAlignment','left','fontsize',10);
  text(ax,0.98,0.94,'b','FontSize',10,'unit','normalized','HorizontalAlignment',...
    'right','EdgeColor','k','Margin',1,'backgroundcolor','w');
  xlim(ax,[0 80]);
  ylim(ax,[0 100]);
%   ylim(ax,[0 70]);
  
  ax=f.ax(3); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
  pltsta=[1 2 3 7];
  scatter(ax,snrf(pltsta),varred4th(pltsta,3),20,'ro');
  scatter(ax,snrfort(pltsta),varredort4th(pltsta,3),24,'bs');
  scatter(ax,snrfvert(pltsta),varredvert4th(pltsta,3),20,'k^');
  scatter(ax,snrf(1),varred4th(1,3),20,'ro','filled');
  scatter(ax,snrfort(1),varredort4th(1,3),24,'bs','filled');
  scatter(ax,snrfvert(1),varredvert4th(1,3),20,'k^','filled');
  xlabel(ax,'LFE template SNR');
  ylabel(ax,'Variance reduction (%)');
  %       title(ax,'4-station srcs');
  text(ax,0.01,0.94,'4-station sources',...
    'Units','normalized','HorizontalAlignment','left','fontsize',10);
  text(ax,0.98,0.94,'c','FontSize',10,'unit','normalized','HorizontalAlignment',...
    'right','EdgeColor','k','Margin',1,'backgroundcolor','w');
  xlim(ax,[0 80]);
  ylim(ax,[0 100]);
%   ylim(ax,[0 70]);
  
  if saveflag
    fname = strcat('vrvssnr',greentype,detecttype,'.pdf');
    print(f.fig,'-dpdf',...
      strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
  end
  
% elseif strcmp(greentype, 'ccstack')
%   widin = 5;  % maximum width allowed is 8.5 inches
%   htin = 3;   % maximum height allowed is 11 inches
%   nrow = 1;
%   ncol = 2;
%   f = initfig(widin,htin,nrow,ncol);
%   pltxran = [0.08 0.98]; pltyran = [0.15 0.98]; % optimal axis location
%   pltxsep = 0.1; pltysep = 0.05;
%   optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
% 
%   ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%   pltsta=[1 2 3];
%   p(1)=scatter(ax,snrf(pltsta),varred(pltsta,3),20,'ro');
%   p(2)=scatter(ax,snrfort(pltsta),varredort(pltsta,3),24,'bs');
%   p(3)=scatter(ax,snrfvert(pltsta),varredvert(pltsta,3),20,'k^');
%   legend(ax,p,'Optimal','Orthogonal','Vertical','location','southeast');
%   xlabel(ax,'Signal-to-noise ratio');
%   ylabel(ax,'Variance reduction (%)');
%   text(ax,0.15,0.94,'Grouped sources',...
%     'Units','normalized','HorizontalAlignment','left','fontsize',10);
%   text(ax,0.02,0.94,'a','FontSize',10,'unit','normalized','HorizontalAlignment',...
%     'left','EdgeColor','k','Margin',1,'backgroundcolor','w');
%   xlim(ax,[0 80]);
%   ylim(ax,[0 70]);
%   
%   pltsta=[1 2 3 7];
%   ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
%   scatter(ax,snrf(pltsta),varred4th(pltsta,3),20,'ro');
%   scatter(ax,snrfort(pltsta),varredort4th(pltsta,3),24,'bs');
%   scatter(ax,snrfvert(pltsta),varredvert4th(pltsta,3),20,'k^');
%   xlabel(ax,'Signal-to-noise ratio');
%   ylabel(ax,'Variance reduction (%)');
%   text(ax,0.15,0.94,'4-station sources',...
%     'Units','normalized','HorizontalAlignment','left','fontsize',10);
%   text(ax,0.02,0.94,'b','FontSize',10,'unit','normalized','HorizontalAlignment',...
%     'left','EdgeColor','k','Margin',1,'backgroundcolor','w');
%   xlim(ax,[0 80]);
%   ylim(ax,[0 70]);
%   
%   if saveflag
%     fname = strcat('vrvssnr',detecttype,'.pdf');
%     print(f.fig,'-dpdf',...
%       strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
%   end
% end







