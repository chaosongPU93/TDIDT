function [f,den1d,conmat,muopt,sigmaopt,mdistprojopt,dprojxopt,countnopt,...
  muort,sigmaort,mdistprojort,dprojxort,countnort] = ...
  plt_srcdloc(dloc,disttype,msize,cstr,marker,scale,binmethod,xran,yran,dx,dy,...
  smoothsigma,ncont)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_srcdloc(dt,dloc,disttype,timetype,bintype)
%
% Function to create a similar plot to figure 4 of Rubin & Gillard (JGRSE,
% 2000) and figure 4 of Rubin (JGRSE, 2002). For each detection, it gets 
% placed at the origin, and then all the other detections in the adopted
% differential time bin (or overall within some range) get plotted in their
% relative position on the fault plane. Here we talk about the event pairs
% composed by either each detection and all the other (that occur within  
% some time range), or consecutive events 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/03/08
% Last modified date:   2024/01/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('dtran',[]);
defval('disttype','km');
defval('timetype','tarvl');
defval('scale','log10');
defval('binmethod','grid');
defval('dx',[]);
defval('dy',[]);
defval('smoothsigma',[]);

nrow = 1; % rows and cols of subplots in each figure
ncol = 3;
widin = 8.4; % size of each figure
htin = 3.75;
pltxran = [0.07 0.98]; pltyran = [0.18 0.88];
pltxsep = 0.02; pltysep = 0.03;
f = initfig(widin,htin,nrow,ncol);
% axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
axpos = [0.05 0.08 0.3 0.88;
         0.36 0.08 0.3 0.88;
         0.72 0.175 0.27 0.69];

for isub = 1:nrow*ncol
  set(f.ax(isub), 'position', axpos(isub,:));
end

%cumulative density, bin type determined by 'bintype'
ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
if strcmp(binmethod,'pixel')
  den1d = density_pixel(dloc(:,1),dloc(:,2));
elseif strcmp(binmethod,'grid')
  den1d = density_matrix(dloc(:,1),dloc(:,2),xran,yran,dx,dy);
end
den1d = den1d(den1d(:,3)>0, :);
den1d = sortrows(den1d,3);

dum = den1d;
dum(dum(:,3)>1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'linew',0.2);  %, 'MarkerEdgeColor', 'k'

dum = sortrows(den1d,3);
dum(dum(:,3)==1, :) = [];
if strcmp(scale,'log10')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),msize,dum(:,3),marker,'filled','MarkerEdgeColor','none');
text(ax,0.98,0.05,sprintf('%d pairs',size(dloc,1)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
colormap(ax,'plasma');
c=colorbar(ax,'SouthOutside');
ax.CLim(2) = prctile(dum(:,3),99);
if strcmp(scale,'log10')
  c.Label.String = strcat('log_{10}(',cstr{1},')');
elseif strcmp(scale,'linear')
  c.Label.String = cstr{1};
end
%Principal component analysis
[coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dloc);
plot(ax,ellx,elly,'-','linew',1.5,'color','k');
plot(ax,x0+semib*[coeff(1,2),-coeff(1,2)],y0+semib*[coeff(2,2),-coeff(2,2)],...
  '-','linew',2,'color','b');
plot(ax,x0+semia*[coeff(1,1),-coeff(1,1)],y0+semia*[coeff(2,1),-coeff(2,1)],...
  ':','linew',2,'color','b');
% text(ax,0.98,0.16,strcat(num2str(round(anglegeo(2))),'$^{\,\circ}$',{'; '},...
%   num2str(round(anglegeo(1))),'$^{\,\circ}$'),'FontSize',11,...
%   'unit','normalized','interpreter','latex','HorizontalAlignment','right');
text(ax,0.98,0.17,sprintf('%d%c; %d%c',round(anglegeo(2)),char(176),...
  round(anglegeo(1)),char(176)),'Units','normalized',...
  'HorizontalAlignment','right','FontSize',9);
text(ax,0.98,0.11,strcat({'Asp. ratio: '},sprintf('%.1f',semia/semib)),'Units',...
  'normalized','HorizontalAlignment','right','FontSize',9);
text(ax,0.02,0.95,'a','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
if strcmp(disttype,'spl')
  xlabel(ax,'Diff off12 (samples)');
  ylabel(ax,'Diff off13 (samples)');
  xticks(ax,xran(1):10:xran(2));
  yticks(ax,yran(1):10:yran(2));
elseif strcmp(disttype,'km')
  xlabel(ax,'Location difference E (km)');
  ylabel(ax,'Location difference N (km)');
  xticks(ax,xran(1):1:xran(2));
  yticks(ax,yran(1):1:yran(2));
end
% plot(ax,ax.XLim,ax.YLim,'k--','linew',1);
% pos = ax.Position;
% c.Position = [pos(1), pos(2)-0.12, pos(3), 0.02];
% 0.08 0.55 0.4 0.4
% c.Position = [0.10 0.515 0.42 0.02];
c.Position = [0.05 0.11 0.3 0.04];
longticks(ax,2);
% nolabels(ax,3);
hold(ax,'off');

% keyboard
%contours of cumulative density
ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on'); axis(ax, 'equal');
if strcmp(disttype,'spl')
  [~,xgrid,ygrid,zgrid] = ...
    zeropadmat2d(den1d,min(den1d(:,1)):1:max(den1d(:,1)),...
    min(den1d(:,2)):1:max(den1d(:,2)));
elseif strcmp(disttype,'km')
  [~,xgrid,ygrid,zgrid] = ...
    zeropadmat2d(den1d,floor(min(den1d(:,1))):dx:ceil(max(den1d(:,1))),...
    floor(min(den1d(:,2))):dy:ceil(max(den1d(:,2))));
end
if ~isempty(smoothsigma)
  zgridgf = imgaussfilt(zgrid,smoothsigma);  %smooth it a bit
else
  zgridgf = zgrid;
end
if strcmp(scale,'log10')
  zgridgf = log10(zgridgf);
end

[conmat,conobj] = contour(ax,xgrid,ygrid,zgridgf,ncont,'-','linew',1);
if ~isempty(conmat)
  contable = getContourLineCoordinates(conmat);
  conmat=table2array(contable);
end
%a line cross (0,0) in the projection direction
x = reshape(xran(1):0.001:xran(2), [], 1);
yopt = linefcn(x,tand(angle(2)),0);
plot(ax,x,yopt,'-','linew',2,'color','b');
%a line cross (0,0) in the orthogonal direction
yort = linefcn(x,tand(angle(1)),0);
plot(ax,x,yort,':','linew',2,'color','b');
colormap(ax,'plasma');
c=colorbar(ax,'SouthOutside');
%     ax.CLim(2) = prctile(dum(:,3),99);
ax.CLim(2) = max(conmat(:,1));
if strcmp(scale,'log10')
  c.Label.String = strcat('log_{10}(',cstr{1},')');
elseif strcmp(scale,'linear')
  c.Label.String = cstr{1};
end
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
if strcmp(disttype,'spl')
  %     xlabel(ax,'Diff off12 (samples)');
  %     ylabel(ax,'Diff off13 (samples)');
  xticks(ax,xran(1):5:xran(2));
  yticks(ax,yran(1):5:yran(2));
elseif strcmp(disttype,'km')
  %     xlabel(ax,'Diff E loc (km)');
  %     ylabel(ax,'Diff N loc (km)');
  xticks(ax,xran(1):1:xran(2));
  yticks(ax,yran(1):1:yran(2));
end
% pos = ax.Position;
% c.Position = [pos(1), pos(2)-0.12, pos(3), 0.02];
c.Position = [0.36 0.11 0.3 0.04];
longticks(ax,2);
nolabels(ax,3);
text(ax,0.02,0.95,'b','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
hold(ax,'off');

%interpolation to obtain the intersection between PCA directions and contours
%cross-sections of contours along the 2 PCA directions
F = scatteredInterpolant(conmat(:,3),conmat(:,4),conmat(:,1),'linear','none');
normalizer=length(dloc);
ax=f.ax(3); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
[ax,muopt,sigmaopt,mdistprojopt,dprojxopt,countnopt]=plt_dloccrssect(ax,F,x,yopt,anglegeo(2),...
  ['-';'-'],[.5 .5 .5; 0 0 1],xran,'both',normalizer);
[ax,muort,sigmaort,mdistprojort,dprojxort,countnort]=plt_dloccrssect(ax,F,x,yort,anglegeo(1),...
  [':';':'],[.5 .5 .5; 0 0 1],xran,'both',normalizer);
text(ax,0.01,0.65,sprintf('SE \\mu=%.2f;\nSE \\sigma=%.2f;\nmed(|x|)=%.2f',...
  muopt,sigmaopt,mdistprojopt),'Units','normalized',...
  'HorizontalAlignment','left','FontSize',9,'color','b');
text(ax,0.65,0.65,sprintf('NE \\mu=%.2f;\nNE \\sigma=%.2f;\nmed(|x|)=%.2f',...
  muort,sigmaort,mdistprojort),'Units','normalized',...
  'HorizontalAlignment','left','FontSize',9,'color','b');
text(ax,0.02,0.05,'c','FontSize',10,'unit','normalized','EdgeColor','k',...
  'Margin',1,'backgroundcolor','w');
axranexp(ax,2,20);
% xlabel(ax,'');
% ylabel(ax,'Probability');
lgd = legend(ax,'SE data','SE gsfit','NE data','NE gsfit',...
  'Location','north','fontsize',7,'NumColumns',2,'orientation','horizontal');
%make background transparent
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));

hold(ax,'off');



% keyboard






