function [f] = plt_abnsrcdist_p2(abne,abnl,stas,abndtarvl,abndist)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_abnsrcdist_p2(abne,abnl,stas,abndtarvl,abndist)
%
% Plot the abnormal sources shown in the plot of source distance ordered 
% by origin time that have a very short differential origin time, but a 
% relatively big distance. Show their characteristics in different aspects.
% 
% This is part 2 to show the arrival time relationship between each abnormal
% pair at each station. 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/06/21
% Last modified date:   2022/06/21 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

f.fig = figure;
f.fig.Renderer = 'painters';

nabn = size(abne,1);
nsta = 3;

for ista = 1: nsta
  subplot(3,2,(ista-1)*2+1); hold on
  for ii = 1: nabn
    plot([abne(ii,(ista-1)*2+1) abnl(ii,(ista-1)*2+1)], ...
      [abne(ii,ista*2) abnl(ii,ista*2)],'k-');
    scatter(abne(ii,(ista-1)*2+1),abne(ii,ista*2),15,'w','filled','MarkerEdgeColor',...
      [.5 .5 .5]);
    scatter(abnl(ii,(ista-1)*2+1),abnl(ii,ista*2),15,'r','filled','MarkerEdgeColor',...
      [.5 .5 .5]);
  end
  text(0.9,0.9,strtrim(stas(ista,:)),'Units','normalized',...
    'HorizontalAlignment','right');
  xlabel('Arrival time in samples');
  ylabel('Amplitude');
  grid on; box on;
end
subplot(3,2,2); hold on
scatter(abnl(:,1),abndtarvl(:,1),15,'r','filled');
scatter(abnl(:,3),abndtarvl(:,2),15,'b','filled');
scatter(abnl(:,5),abndtarvl(:,3),15,'k','filled');
text(0.9,0.9,sprintf('min(mean(abs))=%.1f',min(mean(abs(abndtarvl),2))),'Units','normalized',...
  'HorizontalAlignment','right');
text(0.9,0.8,sprintf('med(mean(abs))=%.1f',median(mean(abs(abndtarvl),2))),'Units','normalized',...
  'HorizontalAlignment','right');
xlabel('Arrival time in samples');
ylabel('differential arrival time');
grid on; box on;

subplot(3,2,4); hold on
scatter(abnl(:,1),abndist,15,'r','filled');
xlabel('Arrival time in samples');
ylabel('Distance between consecutive sources (km)');
grid on; box on;

subplot(3,2,6); hold on
scatter(mean(abnl(:,[2 4 6])./abne(:,[2 4 6]), 2),abndist,15);
xlabel('Mean amplitude ratio');
ylabel('Distance (km)');
grid on; box on;
      
