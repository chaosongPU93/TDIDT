function [f] = plt_abnsrcdist_p1(abnloce,abnlocl,abne,abnl,nsep,sps,abndist,abndt,abntoril,dist,torisplst,mapview)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_abnsrcdist_p1(abnloce,abnlocl,abne,abnl,nsep,sps,abndist,abndt,torisplst,mapview)
%
% Plot the abnormal sources shown in the plot of source distance ordered 
% by origin time that have a very short differential origin time, but a 
% relatively big distance. Show their characteristics in different aspects.
% 
% This is part 1 to show the distance between each abnormal pair. Has the 
% option of plotting map view of sources in 'offset' space (samples at
% specific sps) and 'relaloc' space (N and E in km)
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/06/21
% Last modified date:   2022/06/21 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

defval('mapview','relaloc');  % default is plot map view in N and E in km

f.fig = figure;
f.fig.Renderer = 'painters';

nabn = size(abnloce,1);

subplot(2,2,1); hold on
if isequal(mapview,'relaloc')
  for ii = 1: nabn
    plot([abnloce(ii,1) abnlocl(ii,1)], [abnloce(ii,2) abnlocl(ii,2)],'k-');
    scatter(abnloce(ii,1),abnloce(ii,2),15,'w','filled',...
      'MarkerEdgeColor',[.5 .5 .5]);
    scatter(abnlocl(ii,1),abnlocl(ii,2),15,'r','filled',...
      'MarkerEdgeColor',[.5 .5 .5]);
  end
  axis equal
  xlabel('E (km)');
  ylabel('N (km)');
  axis([-4 4 -4 4]);
  xticks(-4:1:4);
  yticks(-4:1:4);
  plot([-4 4],[-4 4],'b--');
elseif isequal(mapview,'offset')
  for ii = 1: nabn
    plot([abne(ii,7) abnl(ii,7)], [abne(ii,8) abnl(ii,8)],'k-');
    scatter(abne(ii,7),abne(ii,8),15,'w','filled',...
      'MarkerEdgeColor',[.5 .5 .5]);
    scatter(abnl(ii,7),abnl(ii,8),15,'r','filled',...
      'MarkerEdgeColor',[.5 .5 .5]);
  end
  xlabel(sprintf('Offset 12 (samples at %d Hz)',sps),'FontSize',11);
  ylabel(sprintf('Offset 13 (samples at %d Hz)',sps),'FontSize',11);
  axis equal
  %       xran = [round(mean(minmax(off1iw(:,2)'))-span/2)-1, round(mean(minmax(off1iw(:,2)'))+span/2)+1];
  %       yran = [round(mean(minmax(off1iw(:,3)'))-span/2)-1, round(mean(minmax(off1iw(:,3)'))+span/2)+1];
  %       xlim(xran); xticks(xran(1): max(0.5,2) : xran(2));
  %       ylim(yran); yticks(yran(1): max(0.5,2) : yran(2));
  plot([-25 25],[-25 25],'b--');
end
text(0.9,0.9,sprintf('%d',nabn),'Units','normalized',...
  'HorizontalAlignment','right');
grid on; box on;


subplot(2,2,2); hold on
histogram(abndist,'BinWidth',0.1);
xlim([0 3.5]);
xlabel('Distance between consecutive sources (km)');
ylabel('Counts');
grid on; box on;

subplot(2,2,3)
plot(torisplst(1+nsep:end)/sps,dist,'MarkerSize',2); hold on
scatter(abntoril/sps,abndist,15,'r','filled','MarkerEdgeColor',[.5 .5 .5]);
ylim([0 3.5]);
xlabel('Relative origin time (s)');
ylabel('Distance between consecutive sources (km)');
grid on; box on;

subplot(2,2,4); hold on
scatter(abndt/sps,abndist,15);
ylim([0 3.5]);
xlabel('Differential relative origin time (s)');
ylabel('Distance between consecutive sources (km)');
xlim([0 2]);
grid on; box on;




