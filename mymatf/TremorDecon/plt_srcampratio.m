function [f,amprat,ampsclp,ampscln,ampratsclp,ampratscln] = plt_srcampratio(imploc,impindepst,greenf,sps,mapview)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,amprat,ampratscl] = plt_srcampratio(imploc,impindepst,greenf,sps,mapview)
%
% Plot the amplitude ratio between sta 1 and 2, 13, and 23, for the sources
% grouped from independent deconvolution and preserved after removing the 
% secondary sources. Has the option of plotting map view of sources in 
% 'offset' space (samples at specific sps) and 'relaloc' space (N and E in km)
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/06/21
% Last modified date:   2022/06/21 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

defval('mapview','relaloc');  % default is plot map view in N and E in km

f.fig = figure;
f.fig.Renderer = 'painters';

%%% Direct deconvolved impulse amplitude
for ista = 2:3
  subplot(2,3,ista-1); hold on
  if isequal(mapview,'relaloc')
    scatter(imploc(:,1),imploc(:,2),20,log10(impindepst(:,2)./impindepst(:,ista*2)),'filled',...
      'MarkerEdgeColor',[.5 .5 .5]);
    xlabel('E (km)');
    ylabel('N (km)');
    axis equal
    axis([-5 3 -4 4]);
    xticks(-5:1:3);
    yticks(-4:1:4);
  elseif isequal(mapview,'offset')
    scatter(impindepst(:,7),impindepst(:,8),20,log10(impindepst(:,2)./impindepst(:,ista*2)),'filled',...
      'MarkerEdgeColor',[.5 .5 .5]);
    if ista == 2
      xlabel(sprintf('Offset 12 (samples at %d Hz)',sps),'FontSize',11);
      ylabel(sprintf('Offset 13 (samples at %d Hz)',sps),'FontSize',11);
    end
    axis equal
  end
  colormap jet
  c=colorbar;
  c.Label.String = sprintf('log_{10}(Impulse amp ratio 1 / %d)',ista);
  caxis([-0.8 1.2]);
  grid on; box on;
end

subplot(2,3,3); hold on
if isequal(mapview,'relaloc')
  scatter(imploc(:,1),imploc(:,2),20,log10(impindepst(:,4)./impindepst(:,6)),'filled',...
    'MarkerEdgeColor',[.5 .5 .5]);
  axis equal
  axis([-5 3 -4 4]);
  xticks(-5:1:3);
  yticks(-4:1:4);
elseif isequal(mapview,'offset')
  scatter(impindepst(:,7),impindepst(:,8),20,log10(impindepst(:,4)./impindepst(:,6)),'filled',...
    'MarkerEdgeColor',[.5 .5 .5]);
  axis equal
end
colormap jet
c=colorbar;
c.Label.String = sprintf('log_{10}(Impulse amp ratio 2 / 3)');
caxis([-0.8 1.2]);
grid on; box on;

amprat = [impindepst(:,2)./impindepst(:,4) impindepst(:,2)./impindepst(:,6) ...
  impindepst(:,4)./impindepst(:,6)];


%Deconvolved impulse amplitude * amp range of template
for ista = 2:3
  subplot(2,3,ista+2); hold on
  if isequal(mapview,'relaloc')
    scatter(imploc(:,1),imploc(:,2),20,log10(impindepst(:,2)*range(greenf(:,1))./...
      (impindepst(:,ista*2)*max(greenf(:,ista)))),'filled',...
      'MarkerEdgeColor',[.5 .5 .5]);
    xlabel('E (km)');
    ylabel('N (km)');
    axis equal
    axis([-5 3 -4 4]);
    xticks(-5:1:3);
    yticks(-4:1:4);
  elseif isequal(mapview,'offset')
    scatter(impindepst(:,7),impindepst(:,8),20,log10(impindepst(:,2)*range(greenf(:,1))./...
      (impindepst(:,ista*2)*max(greenf(:,ista)))),'filled',...
      'MarkerEdgeColor',[.5 .5 .5]);
    if ista == 2
      xlabel(sprintf('Offset 12 (samples at %d Hz)',sps),'FontSize',11);
      ylabel(sprintf('Offset 13 (samples at %d Hz)',sps),'FontSize',11);
    end
    axis equal
  end
  colormap jet
  c=colorbar;
  c.Label.String = sprintf('log_{10}(Scaled template amp ratio 1 / %d)',ista);
  caxis([-0.8 1.2]);
  grid on; box on;
end

subplot(2,3,6); hold on
if isequal(mapview,'relaloc')
  scatter(imploc(:,1),imploc(:,2),20,log10(impindepst(:,4)*range(greenf(:,2))./...
    (impindepst(:,6)*max(greenf(:,3)))),'filled',...
    'MarkerEdgeColor',[.5 .5 .5]);
  axis equal
  axis([-5 3 -4 4]);
  xticks(-5:1:3);
  yticks(-4:1:4);
elseif isequal(mapview,'offset')
scatter(impindepst(:,7),impindepst(:,8),20,log10(impindepst(:,4)*range(greenf(:,2))./...
    (impindepst(:,6)*max(greenf(:,3)))),'filled',...
    'MarkerEdgeColor',[.5 .5 .5]);
  axis equal
end
colormap jet
c=colorbar;
c.Label.String = sprintf('log_{10}(Scaled template amp ratio 2 / 3)');
caxis([-0.8 1.2]);
grid on; box on;

ampratsclp = [impindepst(:,2)*max(greenf(:,1))./(impindepst(:,4)*max(greenf(:,2))) ...
  impindepst(:,2)*max(greenf(:,1))./(impindepst(:,6)*max(greenf(:,3))) ...
  impindepst(:,4)*max(greenf(:,2))./(impindepst(:,6)*max(greenf(:,3)))];

ampsclp = [impindepst(:,2)*max(greenf(:,1)) impindepst(:,4)*max(greenf(:,2)) ...
  impindepst(:,6)*max(greenf(:,3))];

ampratscln = [impindepst(:,2)*min(greenf(:,1))./(impindepst(:,4)*min(greenf(:,2))) ...
  impindepst(:,2)*min(greenf(:,1))./(impindepst(:,6)*min(greenf(:,3))) ...
  impindepst(:,4)*min(greenf(:,2))./(impindepst(:,6)*min(greenf(:,3)))];

ampscln = [impindepst(:,2)*min(greenf(:,1)) impindepst(:,4)*min(greenf(:,2))...
  impindepst(:,6)*min(greenf(:,3))];

