function f=plt_shortwin_demo_noise(imp,timevec,locxy,imp4th,timevec4th,...
  locxy4th,lsig,sps,ttype,saveflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_shortwin_demo_noise(torisplst,implocorist,lsig,...
%   sps,ttype,wt,torisplst4th,implocorist4th,wt4th)
%
% Similar to 'plt_shortwin_summary_noise.m', this script would plot the detections
% BOTH after 2ndary removed AND after 4th-sta checked, but for noise, without
% seismograms
% --Note that this function do not plot the projection along the propagation
% direction as there should be no expected migration from noise. Something might
% be there because you force individual alignment for each short window.
% --This is to plot a summary figure for one burst window for the purpose of
% demonstration, should be called in 'deconv_ref_4s_exp_4thsta_181demo.m'. 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/03/23
% Last modified date:   2024/03/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('sps',160);
defval('ttype','tori');  % default is sort the sources by origin time
defval('saveflag',1);  %default is short-win detections

%the cut-out boundary of 4-s detections
cutout = 'ellipse';
x0 = 0.2; 
y0 = 0.2;
semia = 1.75;
semib = 1.25;
angrot = 45;
[xcut, ycut] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);

widin = 8.3;  % maximum width allowed is 8.5 inches
htin = 4;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f = initfig(widin,htin,nrow,ncol); %initialize fig

%%% 2ndary removed, map view
%mannually set the locations for each axis
set(f.ax(1), 'position', [0.08 0.1 0.40 0.8]);
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
plot(ax,xcut,ycut,'k-','linew',2);
cran = [0 lsig/sps];
xran = [-4 4];
yran = [-4 4];
msize = 35;
if ~isempty(imp)
  wt = median(imp(:,[2 4 6]),2);
  wtmax = prctile(wt,95); %use percentile in case
  refscl = wt./wtmax;
  refscl(refscl>=1) = 1;  %force the larger amp to be plotted as the same size in case of saturation
  scatter(ax,locxy(:,1),locxy(:,2),msize*refscl,timevec/sps,'filled','o',...
    'MarkerEdgeColor',[.5 .5 .5]);
end
text(ax,0.99,0.05,sprintf('%d events',size(locxy,1)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',8);
% colormap(ax,flipud(colormap(ax,'kelicol')));
% colormap(ax,'jet');
colormap(ax,'viridis');
c1=colorbar(ax);
caxis(ax,cran);
if isequal(ttype,'tori')
  c1.Label.String = sprintf('Relative origin time (s)');
elseif isequal(ttype,'tarvl')
  c1.Label.String = sprintf('Arrival time (s)');
end
c1.Label.FontSize = 9;
scatter(ax,xran(1)+0.1*range(xran),yran(2)-0.05*range(yran),msize,'w','filled',...
  'MarkerEdgeColor',[.5 .5 .5],'linew',1);
text(ax,0.02,0.9,strcat({'Amplitude '},'$\geq$',{' 95th prctile'}) ,'Units','normalized',...
  'HorizontalAlignment','left','FontSize',8,'interpreter','latex');
axis(ax,'equal');
xticks(ax,xran(1):1:xran(2));
yticks(ax,yran(1):1:yran(2));
axis(ax,[xran yran]);
% text(ax,0.99,0.95,'Secondary removed','HorizontalAlignment','right','Units',...
%   'normalized','FontSize',8);  
text(ax,0.99,0.95,'3-station','HorizontalAlignment','right','Units',...
  'normalized','FontSize',10);  
text(ax,0.02,0.05,'a','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
longticks(ax,2);

%%% 4th-sta checked, map view
%mannually set the locations for each axis
% cpos=c.Position;
% ax1pos=f.ax(1).Position;
% set(f.ax(2), 'position', [cpos(1)+cpos(3)+0.12 cpos(2) ax1pos(3) cpos(4)]);
set(f.ax(2), 'position', [0.56 0.1 0.40 0.8]);
ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
plot(ax,xcut,ycut,'k-','linew',2);
if ~isempty(imp4th) && ~isempty(imp)
  wt4th = median(imp4th(:,[2 4 6]),2);
  % wt4thmax = prctile(wt4th,95); %use percentile in case
  refscl4th = wt4th./wtmax;
  refscl4th(refscl4th>=1) = 1;  %force the larger amp to be plotted as the same size in case of saturation
  scatter(ax,locxy4th(:,1),locxy4th(:,2),msize*refscl4th,timevec4th/sps,'filled',...
    'o','MarkerEdgeColor',[.5 .5 .5]);
end
text(ax,0.99,0.05,sprintf('%d events',size(locxy4th,1)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',8);
% colormap(ax,flipud(colormap(ax,'kelicol')));
colormap(ax,'viridis');
c2=colorbar(ax);
caxis(ax,cran);
if isequal(ttype,'tori')
  c2.Label.String = sprintf('Relative origin time (s)');
elseif isequal(ttype,'tarvl')
  c2.Label.String = sprintf('Arrival time (s)');
end
c2.Label.FontSize = 9;
% scatter(ax,xran(1)+0.1*range(xran),yran(2)-0.05*range(yran),msize,'w','filled',...
%   'MarkerEdgeColor',[.5 .5 .5]);
% text(ax,0.02,0.9,'Amplitude\geq 95th prctile','Units','normalized',...
%   'HorizontalAlignment','left','FontSize',8);
axis(ax,'equal');
xticks(ax,xran(1):1:xran(2));
yticks(ax,yran(1):1:yran(2));
axis(ax,[xran yran]);
% text(ax,0.99,0.95,'Checked at KLNB','HorizontalAlignment','right','Units',...
%   'normalized','FontSize',8);  
text(ax,0.99,0.95,'4-station','HorizontalAlignment','right','Units',...
  'normalized','FontSize',10);  
text(ax,0.02,0.05,'b','FontSize',10,'unit','normalized','EdgeColor','k','Margin',1);
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
longticks(ax,2);

if saveflag
  fname = strcat('shortwinsum_noi_demo.pdf');
  print(f.fig,'-dpdf',...
    strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
end

% keyboard