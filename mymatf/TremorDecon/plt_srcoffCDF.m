function [f] = plt_srcoffCDF(impindepstst,torisplst,nsepvect) 
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

ncomp = length(nsepvect);
dt = cell(ncomp,1);
doff = cell(ncomp,1);

for i = 1: ncomp
  nsep = nsepvect(i);
  %between Nth and (N-i)th source
  dtorinni = diffcustom(torisplst,nsep,'forward');
  doffnni = abs(diffcustom(impindepstst(:,7:9),nsep,'forward'));
  dt{i} = dtorinni;
  doff{i} = doffnni;  
end

f.fig=figure;
f.fig.Renderer='Painters';
widin = 10;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
[scrsz, resol] = pixelperinch(1);
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);

nrow = 2;
ncol = ncomp;
for isub = 1:nrow*ncol
  f.ax(isub) = subplot(nrow,ncol,isub);
end

for isub = 1: 1*ncol
  ax = f.ax(isub);
  hold(ax,'on');
  ax.Box='on';
  grid(ax,'on');
  doffnni = doff{isub};
  [cdfval,x] = ecdf(doffnni(:,1)); %between Nth and (N-1)th source
  plot(ax,x,cdfval,'r','linew',1);
  [cdfval,x] = ecdf(doffnni(:,2)); %between Nth and (N-2)th source
  plot(ax,x,cdfval,'b','linew',1);
  [cdfval,x] = ecdf(doffnni(:,3)); %between Nth and (N-3)th source
  plot(ax,x,cdfval,'k','linew',1);
  % plot(ax,[0 1 1 0 0],[0 0 0.5 0.5 0],'Color',[.5 .5 .5],'linew',2);
  text(ax,0.02,0.95,sprintf('Nth - (N-%d)th',nsepvect(isub)),'Units','normalized',...
    'HorizontalAlignment','left','FontSize',10);
  if isub == 1
    text(ax,0.02,0.6,sprintf('[%d, %d]', round(prctile(impindepstst(:,7),2)), ...
      round(prctile(impindepstst(:,7),98))),'Units','normalized',...
      'HorizontalAlignment','left','FontSize',10,'Color','r');
    text(ax,0.02,0.5,sprintf('[%d, %d]', round(prctile(impindepstst(:,8),2)), ...
      round(prctile(impindepstst(:,8),98))),'Units','normalized',...
      'HorizontalAlignment','left','FontSize',10,'Color','b');
    text(ax,0.02,0.4,sprintf('[%d, %d]', round(prctile(impindepstst(:,9),2)), ...
      round(prctile(impindepstst(:,9),98))),'Units','normalized',...
      'HorizontalAlignment','left','FontSize',10,'Color','k');
    legend(ax,'\Delta{off12}','\Delta{off13}','\Delta{off23}','Location','southeast');
    xlabel(ax,'Distance in offset between sources (samples)');
    ylabel(ax,'Empirical CDF');
  end
  longticks(ax,2);
  hold(ax,'off');
end

for isub = 1*ncol+1: nrow*ncol
  ax = f.ax(isub);
  hold(ax,'on');
  ax.Box='on';
  grid(ax,'on');
  doffnni = doff{isub-ncol};
  [pdfval,x] = ksdensity(doffnni(:,1)); %between Nth and (N-1)th source
  plot(ax,x,pdfval,'r','linew',1);
  [pdfval,x] = ksdensity(doffnni(:,2)); %between Nth and (N-2)th source
  plot(ax,x,pdfval,'b','linew',1);
  [pdfval,x] = ksdensity(doffnni(:,3)); %between Nth and (N-3)th source
  plot(ax,x,pdfval,'k','linew',1);
  text(ax,0.02,0.95,sprintf('Nth - (N-%d)th',nsepvect(isub-ncol)),'Units','normalized',...
    'HorizontalAlignment','left','FontSize',10);
  if isub == 1*ncol+1
    xlabel(ax,'Distance in offset between sources (samples)');
    ylabel(ax,'Empirical PDF');
  end
  longticks(ax,2);
  hold(ax,'off');
end




