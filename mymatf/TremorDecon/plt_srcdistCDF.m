function [f] = plt_srcdistCDF(implocst,torisplst,nsepvect)      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_srcdistCDF(implocst,torisplst,nsepvect)  
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

ncomp = length(nsepvect);
dt = cell(ncomp,1);
dist = cell(ncomp,1);
lgdstr = cell(ncomp,1);

for i = 1: ncomp
  nsep = nsepvect(i);
  %between Nth and (N-i)th source
  dtorinni = diffcustom(torisplst,nsep,'forward');
  distnni = sqrt(diffcustom(implocst(:,1),nsep,'forward').^2 + ...
    diffcustom(implocst(:,2),nsep,'forward').^2 );
  dt{i} = dtorinni;
  dist{i} = distnni;
  lgdstr{i} = sprintf('Nth - (N-%d)th',nsep);
end

color = jet(ncomp);

f.fig=figure;
f.fig.Renderer='Painters';
widin = 10;  % maximum width allowed is 8.5 inches
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
% plot(ax,[0 1 1 0 0],[0 0 0.5 0.5 0],'Color',[.5 .5 .5],'linew',2);
legend(ax,lgdstr,'Location','southeast');
% legend(ax,sprintf('Nth - (N-%d)th',nsepvect(1)),...
%   sprintf('Nth - (N-%d)th',nsepvect(2)),...
%   sprintf('Nth - (N-%d)th',nsepvect(3)),'Location','southeast');
xlabel(ax,'Distance between sources (km)');
ylabel(ax,'Empirical CDF');
hold(ax,'off');

ax = f.ax(2);
hold(ax,'on');
ax.Box='on';
grid(ax,'on');
for i = 1: ncomp
  [pdfval,x] = ksdensity(dist{i});
  % h=histogram(rccpair(:,1),'binw',0.01,'Normalization','pdf');
  plot(ax,x,pdfval,'linew',1,'color',color(i,:));
end
xlabel(ax,'Distance between sources (km)');
ylabel(ax,'Empirical PDF');
hold(ax,'off');





