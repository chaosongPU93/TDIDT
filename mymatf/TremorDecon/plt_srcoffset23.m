function [f] = plt_srcoffset23(imploc,impindepst,torispl,indsort,mapview)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_srcoffset23(imploc,impindepst,torispl,indsort,mapview)
%
% Plot the offset 23 from different perspectives. For example, 
% 1. map view in the space set by off12 and off13, or relative location N 
%   and E. 
% 2. Histogram.
% 3. offset vs. origin time (you can also use the arrival time)
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/07/01
% Last modified date:   2022/07/01 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

defval('mapview','relaloc');  % default is plot map view in N and E in km

f.fig = figure;
f.fig.Renderer = 'painters';

subplot(2,2,1); hold on
if isequal(mapview,'relaloc')
  scatter(imploc(:,1),imploc(:,2),30,impindepst(:,9),'filled',...
    'MarkerEdgeColor',[.5 .5 .5]);
  plot([-3 3],[-3 3],'k--','linew',2);
  xlabel('E (km)');
  ylabel('N (km)');
elseif isequal(mapview,'offset')
  scatter(impindepst(:,7),impindepst(:,8),30,impindepst(:,9),'filled',...
    'MarkerEdgeColor',[.5 .5 .5]);
  plot([-25 25],[-25 25],'k--','linew',2);
  colormap jet
  xlabel('Offset 12 (samples)');
  ylabel('Offset 13 (samples)');
end
c=colorbar;
c.Label.String = sprintf('Offset 23 (samples)');
axis equal
grid on; box on;

subplot(2,2,2)
histogram(impindepst(:,9));
xlabel('Offset 23 (samples)');
ylabel('Counts');
grid on; box on;

subplot(2,2,[3 4])
plot(torispl(indsort),impindepst(indsort,9),'MarkerSize',2);
ylabel('Offset 23 (samples)');
xlabel('Relative origin time (samples)');
grid on; box on;

      