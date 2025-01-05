function f=plt_vrvssnr_v2(stas,ista,snrf,snrfort,snrfvert,...
          varred,varredort,varredvert,varred4th,varredort4th,varredvert4th)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to 'plt_vrvssnr', but plot the result for using other templates
% for deconvolution other than '002'. Technically you can call that function,
% but if the change of variance reduction (VR) with iteration at one station 
% has been plotted in the same panel showing 002 as a comparison, then you 
% only need the other 2 panels showing the final VR for 3-sta and 4-sta 
% sources. This function is to plot the other 2 panels only.
% 
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/11/13
% Last modified date:   2024/11/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('saveflag',1);

% if strcmp(greentype, 'dstack')
  
  widin = 5.5;  % maximum width allowed is 8.5 inches
  htin = 3;   % maximum height allowed is 11 inches
  nrow = 1;
  ncol = 2;
  f = initfig(widin,htin,nrow,ncol);
  pltxran = [0.1 0.98]; pltyran = [0.15 0.98]; % optimal axis location
  pltxsep = 0.1; pltysep = 0.05;
  optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
  
  ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
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
  text(ax,0.98,0.94,'a','FontSize',10,'unit','normalized','HorizontalAlignment',...
    'right','EdgeColor','k','Margin',1,'backgroundcolor','w');
  xlim(ax,[0 80]);
  % ylim(ax,[0 100]);
  ylim(ax,[-20 100]);
  
  ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
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
  text(ax,0.98,0.94,'b','FontSize',10,'unit','normalized','HorizontalAlignment',...
    'right','EdgeColor','k','Margin',1,'backgroundcolor','w');
  xlim(ax,[0 80]);
  % ylim(ax,[0 100]);
  ylim(ax,[-20 100]);

  
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







