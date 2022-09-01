function [f] = plt_srcdist(implocst,impindepstst,nsep,sps,dist,dt,tsplst,ttype,mapview)
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

defval('ttype','tori');  % default is sort the sources by origin time
defval('mapview','relaloc');  % default is plot map view in N and E in km

f.fig = figure;
f.fig.Renderer = 'painters';

subplot(2,2,1); hold on
if isequal(mapview,'relaloc')
  scatter(implocst(1+nsep:end,1),implocst(1+nsep:end,2),25,dist,'filled',...
    'MarkerEdgeColor',[.5 .5 .5]);
elseif isequal(mapview,'offset')
  scatter(impindepstst(1+nsep:end,7),impindepstst(1+nsep:end,8),25,dist,'filled',...
    'MarkerEdgeColor',[.5 .5 .5]);
end
text(0.9,0.9,sprintf('%d',size(implocst,1)-nsep),'Units','normalized',...
  'HorizontalAlignment','right');
colormap jet
c=colorbar;
c.Label.String = sprintf('Dist. between consecutive sources (km)');
caxis([0 3.5]);
axis equal
if isequal(mapview,'relaloc')
  xlabel('E (km)');
  ylabel('N (km)');
  axis([-4 4 -4 4]);
  xticks(-4:1:4);
  yticks(-4:1:4);
  plot([-4 4],[-4 4],'k--','linew',2);
elseif isequal(mapview,'offset')
  xlabel(sprintf('Offset 12 (samples at %d Hz)',sps),'FontSize',11);
  ylabel(sprintf('Offset 13 (samples at %d Hz)',sps),'FontSize',11);
  %       xran = [round(mean(minmax(off1iw(:,2)'))-span/2)-1, round(mean(minmax(off1iw(:,2)'))+span/2)+1];
  %       yran = [round(mean(minmax(off1iw(:,3)'))-span/2)-1, round(mean(minmax(off1iw(:,3)'))+span/2)+1];
  %       xlim(xran); xticks(xran(1): max(0.5,2) : xran(2));
  %       ylim(yran); yticks(yran(1): max(0.5,2) : yran(2));
  plot([-25 25],[-25 25],'k--','linew',2);
end
grid on; box on;

subplot(2,2,2); hold on
grid on; box on;
yyaxis left
histogram(dist,'BinWidth',0.1);
text(0.9,0.9,sprintf('med: %.2f km',median(dist)),'Units','normalized',...
  'HorizontalAlignment','right');
% xlim([0 3.5]);
xlabel('Dist. between consecutive sources (km)');
ylabel('Counts');
yyaxis right
[cdfval,x] = ecdf(dist); %between Nth and (N-1)th source
plot(x,cdfval,'r','linew',1);
ax=gca;
plot([median(dist) ax.XLim(2)],[0.5 0.5],'r--');
plot([median(dist) median(dist)],[0 0.5],'r--');
ylabel('Empirical CDF');

subplot(2,2,3)
plot(tsplst(1+nsep:end)/sps,dist,'MarkerSize',2);
if isequal(ttype,'tori')
  xlabel('Relative origin time (s)');
elseif isequal(ttype,'tarvl')
  xlabel('Arrival time (s)');
end
ylabel('Dist. between consecutive sources (km)');
grid on; box on;

subplot(2,2,4); hold on
grid on; box on;
yyaxis left
scatter(dt/sps,dist,15);
text(0.9,0.9,strcat({'med of dt\leq1s: '},sprintf('%.2f km',median(dist(dt/sps<=1)))),...
  'Units','normalized','HorizontalAlignment','right');
text(0.9,0.75,strcat({'med of dt\leq0.5s: '},sprintf('%.2f km',median(dist(dt/sps<=0.5)))),...
  'Units','normalized','HorizontalAlignment','right');
plot([0.1 0.1],[0 3.5],'k--');
plot([0 2],[median(dist(dt/sps<=1)) median(dist(dt/sps<=1))],'k--');
% ylim([0 3.5]);
if isequal(ttype,'tori')
  xlabel('Diff. relative origin time (s)');
elseif isequal(ttype,'tarvl')
  xlabel('Diff. arrival time (s)');
end
ylabel('Dist. between consecutive sources (km)');
xlim([0 2]);
yyaxis right
[cdfval,x] = ecdf(dt/sps); %between Nth and (N-1)th source
plot(x,cdfval,'r','linew',1);
ax=gca;
plot([median(dt/sps) ax.XLim(2)],[0.5 0.5],'r--');
plot([median(dt/sps) median(dt/sps)],[0 0.5],'r--');
ylabel('Empirical CDF');






