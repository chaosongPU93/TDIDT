function [f] = plt_srcprojdistCDF(implocst,nsepvect,sps,torisplst)      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_srcdist(implocst,impindepstst,nsep,sps,dist,dt,torisplst,ttype,mapview)
%
% Plot the distance between the source N and source N-nsep, in the sequential
% order of origin time or arrival time (specified by 'ttype'). Has the option 
% of plotting map view of sources in 'offset' space (samples at specific sps)
% and 'relaloc' space (N and E in km).
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/06/21
% Last modified date:   2022/06/21 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% defval('ttype','tori');  % default is sort the sources by origin time

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


f.fig = figure;
f.fig.Renderer = 'painters';
widin = 6;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
[scrsz, resol] = pixelperinch(1);
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);

nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
  f.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
% set(f.ax(1), 'position', [ 0.1, 0.1, 0.8, 0.8]);
ax = f.ax(1);
%%% Actually used is the pre-determined best prop direc to do the fitting
implocdum = implocst;
for jj = 1: size(implocst,1)
  x0 = implocdum(jj,1);
  y0 = implocdum(jj,2);
  [newx,newy] = coordinate_rot(x0,y0,-(angrmse-90),[0 0]);
  implocdum(jj,1) = newx;
  implocdum(jj,2) = newy;
end
ncomp = length(nsepvect);
dt = cell(ncomp,1);
dist = cell(ncomp,1);
lgdstr = cell(ncomp,1);
for i = 1: ncomp
  nsep = nsepvect(i);
  dtorinni = diffcustom(torisplst,nsep,'forward');
  distpropnni = abs(diffcustom(implocdum(:,1),nsep,'forward'));
  dt{i} = dtorinni;
  dist{i} = distpropnni;
  lgdstr{i} = sprintf('Nth - (N-%d)th',nsep);
end
color = jet(ncomp);
hold(ax,'on');
ax.Box='on';
grid(ax,'on');
for i = 1: ncomp
  [cdfval,x] = ecdf(dist{i}); %between Nth and (N-1)th source
  plot(ax,x,cdfval,'linew',1,'color',color(i,:));
end
% [cdfval,x] = ecdf(dist{1}); %between Nth and (N-1)th source
% plot(ax,x,cdfval,'r','linew',1);
% [cdfval,x] = ecdf(dist{2}); %between Nth and (N-2)th source
% plot(ax,x,cdfval,'b','linew',1);
% [cdfval,x] = ecdf(dist{3}); %between Nth and (N-3)th source
% plot(ax,x,cdfval,'k','linew',1);
% plot(ax,[0 0.5 0.5 0 0],[0 0 0.7 0.7 0],'Color',[.5 .5 .5],'linew',2);
legend(ax,lgdstr,'Location','southeast');
% legend(ax,sprintf('Nth - (N-%d)th',nsepvect(1)),...
%   sprintf('Nth - (N-%d)th',nsepvect(2)),...
%   sprintf('Nth - (N-%d)th',nsepvect(3)),'Location','southeast');
xlabel(ax,'Proj. dist. along prop. (min rmse) between sources (km)');
ylabel(ax,'Empirical CDF');
hold(ax,'off');

% set(f.ax(2), 'position', [ 0.46, 0.17, 0.4, 0.4]);
% ax = f.ax(2);
% hold(ax,'on');
% ax.Box='on';
% grid(ax,'on');
% [cdfval,x] = ecdf(dist{1});
% plot(ax,x,cdfval,'r','linew',1);
% [cdfval,x] = ecdf(dist{2});
% plot(ax,x,cdfval,'b','linew',1);
% [cdfval,x] = ecdf(dist{3});
% plot(ax,x,cdfval,'k','linew',1);
% xlim(ax,[0 0.5]);
% hold(ax,'off');

ax = f.ax(2);
%%% Actually used is the pre-determined best prop direc to do the fitting
implocdum = implocst;
for jj = 1: size(implocst,1)
  x0 = implocdum(jj,1);
  y0 = implocdum(jj,2);
  [newx,newy] = coordinate_rot(x0,y0,-(angslope-90),[0 0]);
  implocdum(jj,1) = newx;
  implocdum(jj,2) = newy;
end
ncomp = length(nsepvect);
dt = cell(ncomp,1);
dist = cell(ncomp,1);
for i = 1: ncomp
  nsep = nsepvect(i);
  dtorinni = diffcustom(torisplst,nsep,'forward');
  distpropnni = abs(diffcustom(implocdum(:,1),nsep,'forward'));
  dt{i} = dtorinni;
  dist{i} = distpropnni;
end
hold(ax,'on');
ax.Box='on';
grid(ax,'on');
for i = 1: ncomp
  [cdfval,x] = ecdf(dist{i}); %between Nth and (N-1)th source
  plot(ax,x,cdfval,'linew',1,'color',color(i,:));
end
% [cdfval,x] = ecdf(dist{1}); %between Nth and (N-1)th source
% plot(ax,x,cdfval,'r','linew',1);
% [cdfval,x] = ecdf(dist{2}); %between Nth and (N-2)th source
% plot(ax,x,cdfval,'b','linew',1);
% [cdfval,x] = ecdf(dist{3}); %between Nth and (N-3)th source
% plot(ax,x,cdfval,'k','linew',1);
% plot(ax,[0 0.5 0.5 0 0],[0 0 0.7 0.7 0],'Color',[.5 .5 .5],'linew',2);
xlabel(ax,'Proj. dist. along prop. (max speed) between sources (km)');
hold(ax,'off');














