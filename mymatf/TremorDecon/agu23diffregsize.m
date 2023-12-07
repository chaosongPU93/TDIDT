function f=agu23diffregsize(xran,yran)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% agu23diffregsize.m
%
% This script is to plot the schematic diagram of how synthetics are generated
% from different size of regions with no added noise. Basically to plot a 
% a series of concentric circles.
% 
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/12/04
% Last modified date:   2023/12/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrow = 1;
ncol = 1;
widin = 4*ncol;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
[f] = initfig(widin,htin,nrow,ncol);

pltxran = [0.08 0.95]; pltyran = [0.15 0.9];
pltxsep = 0.08; pltysep = 0.05; 
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

% xran = [-6 4];
% yran = [-4 4];

ax = f.ax(1);
hold(ax,'on');
ax.Box = 'on';
grid(ax, 'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
semia = 1.75*(0.6:0.2:2.0);
semib = 1.25*(0.6:0.2:2.0);
nreg = length(semia);
for ireg = 1: nreg
  xaxis = semia(ireg);
  yaxis = semib(ireg);
  % xaxis=1.75; %axis length of the same ellipse of my 4-s catalog
  % yaxis=1.25;
  shiftor=[0.2 0.2]; %(in km) %center of the same ellipse of my 4-s catalog
  [xcut,ycut] = ellipse_chao(shiftor(1),shiftor(2),xaxis,yaxis,0.01,45,shiftor);
  if ireg == 3
    scatter(ax,shiftor(1),shiftor(2),10,'k','filled');
    plot(ax,xcut,ycut,'k-','linew',2);
  else
    plot(ax,xcut,ycut,'-','Color',[.5 .5 .5],'linew',1);
  end
   
end
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
% xticks(ax,xran(1):5:xran(2));
% yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
longticks(ax,2);
hold(ax,'off');

