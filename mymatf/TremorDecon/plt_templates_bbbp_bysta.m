function f=plt_templates_bbbp_bysta(green,greenf,stas,staplt,lowlet,hiwlet,sps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_templates_bbbp_bysta(green,greenf,stas,staplt,lowlet,hiwlet,sps)
%
% this function simply plot the broadband (BB) and bandpassed (BP) filtered
% templates (or green's functions), on one component (optimal, orthogonal, 
% or vertical). BB and BP all in one panel for each station.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/08/29
% Last modified date:   2024/08/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

compname = {'Broadband',sprintf('Band-passed, %.1f-%.1f Hz',lowlet,hiwlet)};

nsta = length(staplt);
ncol = 1;
nrow = 1;
widin = 5;  % maximum width allowed is 8.5 inches
htin = 5.5;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);

pltxran = [0.1 0.98]; pltyran = [0.1 0.98]; % optimal axis location
pltxsep = 0.04; pltysep = 0.015;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

green=green(:,staplt);
greenf=greenf(:,staplt);

mx=1.2*max(max(abs(green(:,:))));

ax = f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
lwlet = size(greenf,1);

vsep = 0.5;
ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
xlabel(ax,'Time (s)','fontsize',10);
ylabel(ax,'Amplitude','fontsize',10);

for i = 1: nsta
  p(1)=plot(ax,(1:lwlet)/sps,green(:,i)+i*vsep,'-','color','k','linew',1);
  p(2)=plot(ax,(1:lwlet)/sps,greenf(:,i)+i*vsep,'-','color','r','linew',1);
  text(ax,0.5,i*vsep+0.15,sprintf('%s',stas(staplt(i),:)),'HorizontalAlignment',...
    'left','fontsize',11);
  xlim(ax,[0 12]);
  ylim(ax,[0 (nsta+1)*vsep]);  
  longticks(ax,6);
  
end
lgd=legend(ax,p,compname,'Location','northeast','Orientation','horizontal');
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
hold(ax,'off');


