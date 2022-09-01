function [f] = plt_srcalldistCDF(impindepstst,implocst,torisplst,septime,sps)
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

defval('septime',[]);

nsrc = size(impindepstst,1);
dist = cell(nsrc, 1);
for i = 1: nsrc
  if ~isempty(septime)
    ind = find(abs(torisplst-torisplst(i))<septime);   %closer in time, if not all
  else
    ind = 1:nsrc;
  end
  ind = setdiff(ind,i);
  dist{i} = sqrt((implocst(i,1)-implocst(ind,1)).^2 + (implocst(i,2)-implocst(ind,2)).^2);
%   temp = [];
%   for j = 1: length(ind)
%     temp(j) = sqrt((implocst(i,1)-implocst(ind(j),1))^2 + (implocst(i,2)-implocst(ind(j),2))^2);
%   end
%   dist{i} = temp;
  
end

color = jet(nsrc);

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
for i = 1: nsrc
  if ~isempty(dist{i})
    [cdfval,x] = ecdf(dist{i}); %between Nth and (N-1)th source
  %   plot(ax,x,cdfval,'linew',1,'color',color(i,:));
    plot(ax,x,cdfval,'linew',1,'color',[0 0 0  0.15]);
    xtar = 1;
    [~,ind] = min(abs(x-xtar));
    cdftar(i) = cdfval(ind); 
  end
end
median(cdftar)
% [cdfval,x] = ecdf(dist{1}); %between Nth and (N-1)th source
% plot(ax,x,cdfval,'r','linew',1);
% [cdfval,x] = ecdf(dist{2}); %between Nth and (N-2)th source
% plot(ax,x,cdfval,'b','linew',1);
% [cdfval,x] = ecdf(dist{3}); %between Nth and (N-3)th source
% plot(ax,x,cdfval,'k','linew',1);
% plot(ax,[0 1 1 0 0],[0 0 0.5 0.5 0],'Color',[.5 .5 .5],'linew',2);
% legend(ax,lgdstr,'Location','southeast');
% legend(ax,sprintf('Nth - (N-%d)th',nsepvect(1)),...
%   sprintf('Nth - (N-%d)th',nsepvect(2)),...
%   sprintf('Nth - (N-%d)th',nsepvect(3)),'Location','southeast');
xlabel(ax,'Distance between each source to all others (km)');
ylabel(ax,'Empirical CDF');
if ~isempty(septime)
  text(ax,0.9,0.1,sprintf('Sources within %d s',round(septime/sps)),'Units','normalized',...
    'HorizontalAlignment','right');
else
  text(ax,0.9,0.1,sprintf('All Sources'),'Units','normalized',...
    'HorizontalAlignment','right');
end
hold(ax,'off');

ax = f.ax(2);
hold(ax,'on');
ax.Box='on';
grid(ax,'on');
for i = 1: nsrc
  if ~isempty(dist{i})
    [pdfval,x] = ksdensity(dist{i});
    % h=histogram(rccpair(:,1),'binw',0.01,'Normalization','pdf');
    plot(ax,x,pdfval,'linew',1,'color',[0 0 0  0.15]);
  end
end
xlabel(ax,'Distance between each source to all others (km)');
ylabel(ax,'Empirical PDF');
hold(ax,'off');

% keyboard
