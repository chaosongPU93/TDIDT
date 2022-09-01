function [f] = plt_agu2022abstract(sigsta,resgrp,impindepst,implocst,torisplst,nsep,sps,tmaxi,tmaxo,tbosti,...
  tbosto,ircccat,rcccat,irccrcat,rccrcat,off1iw,off1i,loff_max,tstbuf,tedbuf,dy,mo,yr,k,ftrans,xcut,ycut)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_agu2022abstract(sigsta,resgrp,impindepst,implocst,torisplst,nsep,sps,tmaxi,tmaxo,tbosti,...
%   tbosto,ircccat,rcccat,irccrcat,rccrcat,off1iw,off1i,loff_max,tstbuf,tedbuf,dy,mo,yr,k,ftrans,xcut,ycut)
%
% This function would create a summary figure for each data win, so that we 
% will have a feeling for all migrations. Similar to 'plt_decon_summary.m', 
% but this deal with the refinement by using subwins and a larger grouping 
% boundary in each subwin.
%
% --It contains the comparison between
%   the deconvolved positive peaks and signal waveform peaks at each station;
%   the signal waveform with the 4-s tremor detections and M. Bostock's LFE
%   detections indicated; residual waveform with the arrival time of deconvolved
%   impulses at sta 1 indicated; scatter of sources in terms of offsets, 
%   accounting for prealignment offset; the scatter of sources in terms of rela 
%   locations.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/06/29
% Last modified date:   2022/06/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%compose a summary figure for each data win and save it, so that we will have a feeling for
%%%all migrations
f.fig = figure;
f.fig.Renderer = 'painters';
widin = 8.5;  % maximum width allowed is 8.5 inches
htin = 10;   % maximum height allowed is 11 inches
% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, resol] = pixelperinch(1);
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);
nrow = 7;
ncol = 1;
for isub = 1:nrow*ncol
  f.ax(isub) = subplot(nrow,ncol,isub);
end

%get the locations for each axis
axpos = [0.08 0.87 0.85 0.1;
         0.08 0.755 0.85 0.1;
         0.08 0.34 0.4 0.4;
         0.55 0.34 0.4 0.4;
         0.08 0.09 0.25 0.25;
         0.415 0.09 0.25 0.25;
         0.68 0.09 0.25 0.25;
         ];
for isub = 1:nrow*ncol
  set(f.ax(isub), 'position', axpos(isub,:));
end

color = ['r';'b';'k'];

ym = max(abs(sigsta(:)));
yran=1.2*[-ym ym];

lsig = size(sigsta,1); 
nsta = size(sigsta,2); 

%%%seismograms of signal
ax=f.ax(1);
hold(ax,'on');
yyaxis(ax,'right');
% plot(ax,ircccat/sps,rcccat,'o','color',[.7 .7 .7],'markersize',1);
plot(ax,ircccat/sps,rcccat,'color',[.7 .7 .7]);
ylim(ax,[-1 1]);
yyaxis(ax,'left');
for i = 1: nsta
  plot(ax,(1:lsig)/sps, sigsta(:,i), '-','Color',color(i,:)); 
end
xlim(ax,[0,lsig/sps]);  ylim(ax,yran); 
ax.Box='on'; grid(ax,'on');
%plot the strongest 0.5-s arrival outside ellipse
ind = find(tmaxo>=tstbuf & tmaxo<=tedbuf);
yloc = (yran(1)+range(yran)*0.05);
for ii = 1: length(ind)
  barst = tmaxo(ind(ii))-tstbuf; % arrival of strongest 0.5-s, in sec
  bared = tmaxo(ind(ii))-tstbuf+0.5; % last for 0.5 s
  plot(ax,[barst bared],[yloc yloc],'-','linew',2.5,'color',[.6 .6 .6]);
end
%plot the strongest 0.5-s arrival inside ellipse
ind = find(tmaxi>=tstbuf & tmaxi<=tedbuf);
for ii = 1: length(ind)
  barst = tmaxi(ind(ii))-tstbuf; % arrival of strongest 0.5-s, in sec
  bared = tmaxi(ind(ii))-tstbuf+0.5; % last for 0.5 s
  plot(ax,[barst bared],[yloc yloc],'k-','linew',2.5);
end
%plot the bostock's LFE catalog outside rectangle
ind = find(tbosto>=tstbuf & tbosto<=tedbuf);
scatter(ax,tbosto(ind)-tstbuf, yloc*ones(size(tbosto(ind))),10,'m','linew',1); % bostock
%plot the bostock's LFE catalog inside rectangle
ind = find(tbosti>=tstbuf & tbosti<=tedbuf);
scatter(ax,tbosti(ind)-tstbuf, yloc*ones(size(tbosti(ind))),10,'m','filled');
text(ax,0.98,0.9,'Signal','Units','normalized','HorizontalAlignment','right','FontSize',10);
text(ax,0.01,0.9,'a','FontSize',11,'unit','normalized');
longticks(ax,3); nolabels(ax,1); 
hold(ax,'off');

%%%seismograms of residual
ax=f.ax(2);
hold(ax,'on');
yyaxis(ax,'right');
% [irccr,rccr12] = RunningCC(resgrp(:,1), resgrp(:,2), mwlen);
% [~,rccr13] = RunningCC(resgrp(:,1), resgrp(:,3), mwlen);
% [~,rccr23] = RunningCC(resgrp(:,2), resgrp(:,3), mwlen);
% rccr = (rccr12+rccr13+rccr23)/3;
% plot(ax,irccr/sps,rccr,'o','color',[.7 .7 .7],'markersize',1.5);
% plot(ax,irccr/sps,rccr,'color',[.7 .7 .7]);
plot(ax,irccrcat/sps,rccrcat,'color',[.7 .7 .7]);
ylim(ax,[-1 1]);
ylabel(ax,'Running CC','FontSize',10);
yyaxis(ax,'left');
for i = 1: nsta
  plot(ax,(1:lsig)/sps, resgrp(:,i), '-','Color',color(i,:)); 
end
xlim(ax,[0,lsig/sps]);  ylim(ax,yran); 
ax.Box='on'; grid(ax,'on');
scatter(ax,impindepst(:,1)/sps,yloc*ones(size(impindepst(:,1))),8,'r','filled');
text(ax,0.5,0.9,num2str(size(impindepst,1)),'fontsize',10,...
  'Units','normalized');
ylabel(ax,'Amplitude','FontSize',10);
text(ax,0.98,0.9,'Residual','Units','normalized','HorizontalAlignment','right','FontSize',10);
%       xlabel(sprintf('Samples at %d Hz',round(1/dt)),'FontSize',12);
xlabel(ax,sprintf('Time (s) since %.1f s on %s %s %s',tstbuf,dy,mo,yr),'FontSize',10);
text(ax,0.01,0.9,'b','FontSize',11,'unit','normalized');
longticks(ax,3); 
hold(ax,'off');

%%%plot the scatter of offsets, accounting for prealignment offset, == true offset
ax=f.ax(3);
span = max(range(off1iw(:,2))+2*loff_max, range(off1iw(:,3))+2*loff_max);
xran = [round(mean(minmax(off1iw(:,2)'))-span/2)-1, round(mean(minmax(off1iw(:,2)'))+span/2)+1];
yran = [round(mean(minmax(off1iw(:,3)'))-span/2)-1, round(mean(minmax(off1iw(:,3)'))+span/2)+1];
cran = [0 lsig];
[ax,~,~,~,~,c] = plt_decon_imp_scatter_ref(ax,impindepst,xran,yran,cran,off1iw,loff_max,...
  sps,35,'mean','tori','comb');
scatter(ax,off1i(k,2),off1i(k,3),20,'ks','filled','MarkerEdgeColor','k');
text(ax,0.98,0.04,'c','FontSize',11,'unit','normalized','HorizontalAlignment','right');
xticks(ax,xran(1): 4 : xran(2));
yticks(ax,yran(1): 4 : yran(2));
% pos = ax.Position;
% c.Position = [pos(1)+pos(3)+0.07, pos(2), 0.02, 0.32];
hold(ax,'off');

%%%plot the scatter of sources in terms of rela locations
ax=f.ax(4);
xran = [-5 5];
yran = [-5 5];
cran = [0 lsig/sps];
[ax,~,~,~,c] = plt_decon_imp_scatter_space_ref(ax,impindepst,xran,yran,cran,off1iw,loff_max,...
  sps,35,ftrans,'mean','tori','comb');
text(ax,0.98,0.04,'d','FontSize',11,'unit','normalized','HorizontalAlignment','right');
% pos = ax.Position;
% c.Position = [pos(1)+pos(3)+0.085, pos(2), 0.02, pos(3)];
plot(ax,xcut,ycut,'k-','linew',2);
xticks(ax,xran(1): 1 : xran(2));
yticks(ax,xran(1): 1 : xran(2));

angle = 0:5:355;
slope = zeros(length(angle),1);
rmse = zeros(length(angle),1);
fttpfree = fittype( @(a,b,x) a*x+b);

%%% find best angle, now is mainly to get the variation of se and slope with the trial angle
for iang = 1: length(angle)
  %%% propagation trial
  implocdum = implocst;
  for jj = 1: size(implocst,1)
    x0 = implocdum(jj,1);
    y0 = implocdum(jj,2);
    [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),[0 0]);
    implocdum(jj,1) = newx;
    implocdum(jj,2) = newy;
  end
  % linear robust least square
  [fitobj,gof,~] = fit(torisplst/sps, implocdum(:,1),fttpfree,'Robust','Bisquare',...
    'StartPoint',[1 1]);
  % output fit parameters
  coef = coeffvalues(fitobj);
  slope(iang) = coef(1);
  rmse(iang) = gof.rmse;
end

%%% best angle estimate from hf
ind = find(slope>0);
ind3 = find(rmse(ind)==min(rmse(ind)));     % rmse, Root Mean Squared Error, the smaller, the better
if length(ind3) > 1
  disp(strcat({'multiple angles have same rmse: '},num2str(angle(ind3))));
end
angrmse = angle(ind(ind3(1)));
ind6 = find(slope==max(slope)); % one with the largest slope, i.e., migrating speed
if length(ind6) > 1
  disp(strcat({'multiple angles have same slope:'},num2str(angle(ind6))));
end
angslope = angle(ind6(1));

[rotx, roty] = complex_rot(0,1,-angrmse);
xvect = [0.5-rotx 0.5+rotx];
yvect = [-2.5-roty -2.5+roty];
drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1);
% text(ax,0.63,0.2,strcat(num2str(angrmse),'^{\circ}'),'FontSize',9,...
%   'unit','normalized','horizontalalignment','left');
text(ax,0.63,0.2,strcat(num2str(angrmse),'$^{\circ}$'),'FontSize',9,...
  'unit','normalized','interpreter','latex');
hold(ax,'off');


ax=f.ax(5);
hold(ax,'on');
%%% Actually used is the pre-determined best prop direc to do the fitting
implocdum = implocst;
for jj = 1: size(implocst,1)
  x0 = implocdum(jj,1);
  y0 = implocdum(jj,2);
  [newx,newy] = coordinate_rot(x0,y0,-(angrmse-90),[0 0]);
  implocdum(jj,1) = newx;
  implocdum(jj,2) = newy;
end
p2=scatter(ax,torisplst/sps,implocdum(:,2),15,'k');
p1=scatter(ax,torisplst/sps,implocdum(:,1),15,[.5 .5 .5],'filled','o',...
  'MarkerEdgeColor','k');
% text(ax,0.9,0.9,sprintf('%d',size(implocst,1)),'Units','normalized',...
%   'HorizontalAlignment','right');
% linear robust least square
[fitobjprop,~,~] = fit(torisplst/sps, implocdum(:,1),fttpfree,'Robust','Bisquare',...
  'StartPoint',[1 1]);
% output fit parameters
coefprop = coeffvalues(fitobjprop);
slopeprop = coefprop(1);
intcptprop = coefprop(2);
fitprop = feval(fitobjprop,torisplst/sps);
plot(ax,torisplst/sps,fitprop,'-','linewidth',2,'color','k');
legend(ax,[p1,p2],'Along prop. (min rmse)', 'Along ort. to prop.','Location','south');
xlabel(ax,'Relative origin time (s)','FontSize',10);
ylabel(ax,'Projected location (km)','FontSize',10);
ax.Box='on'; grid(ax,'on');
xlim(ax,[0 ceil(max(torisplst/sps))]);
text(ax,0.04,0.96,'e','FontSize',11,'unit','normalized');
longticks(ax,2);
hold(ax,'off');

ax=f.ax(6);
hold(ax,'on');
dt = diffcustom(torisplst,nsep,'forward');
distort = abs(diffcustom(implocdum(:,2),nsep,'forward'));
p2=scatter(ax,dt/sps,distort,15,'k');
distprop = abs(diffcustom(implocdum(:,1),nsep,'forward'));
p1=scatter(ax,dt/sps,distprop,15,[.5 .5 .5],'filled','o','MarkerEdgeColor','k');
% text(ax,0.9,0.8,strcat({'med of dt\leq1.0s: '},sprintf('%.2f km',median(distprop(dt/sps<=1)))),...
%   'Units','normalized','HorizontalAlignment','right');
% text(ax,0.9,0.7,strcat({'med of dt\leq0.5s: '},sprintf('%.2f km',median(distprop(dt/sps<=0.5)))),...
%   'Units','normalized','HorizontalAlignment','right');
% ylim([0 3.5]);
xlabel(ax,'Diff. relative origin time (s)','FontSize',10);
ylabel(ax,'Projected dist. between sources (km)','FontSize',10);
xlim(ax,[0 2]);
ax.Box='on'; grid(ax,'on');
text(ax,0.04,0.96,'f','FontSize',11,'unit','normalized');
longticks(ax,2);
hold(ax,'off');

ax=f.ax(7);
hold(ax,'on');
p2=histogram(ax,distort,'BinWidth',0.1,'FaceColor','w');
p1=histogram(ax,distprop,'BinWidth',0.1,'FaceColor','k');
text(ax,0.9,0.7,sprintf('med: %.2f km',median(distprop)),'Units','normalized',...
  'HorizontalAlignment','right');
% xlim([0 3.5]);
% xlabel(ax,'Projected dist. between consecutive sources (km)','FontSize',10);
ylabel(ax,'Counts','FontSize',10);
ax.Box='on'; grid(ax,'on');
view(ax,-90, 90);
set(ax, 'ydir', 'reverse');
nolabels(ax,1);
longticks(ax,2);
xlim(ax,f.ax(6).YLim);
text(ax,0.04,0.96,'g','FontSize',11,'unit','normalized');
hold(ax,'off');







