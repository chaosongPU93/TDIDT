function [f,dt,dist,cdftar] = plt_srcalldistCDF(tsplst,implocst,septran,sps,xtar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,dt,dist,cdftar] = plt_srcalldistCDF(tsplst,implocst,maxsept,sps,xtar)
%
% Plot the distance between each LFE and all other LFEs, in the sequential
% order of origin time or arrival time (specified by 'ttype'). Can specify 
% the range of separation in time in terms of samples, say you only care 
% about distance between sources whose time separation is within 1 s. Also
% accepts a target distance reference to know the fraction of sources whose
% distance to others is larger than that reference. Return the handle of 
% figure, distance cell array, and fraction of distance to others beyond ref
% for each source.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/06/21
% Last modified date:   2023/02/15 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

defval('septran',[]);
defval('xtar',1);

nsrc = size(implocst,1);
dist = cell(nsrc, 1);
%loop for every source
for i = 1: nsrc
  if ~isempty(septran)
    ind = find(abs(tsplst-tsplst(i))>=septran(1) & abs(tsplst-tsplst(i))<=septran(2));   %closer in time, if not all
  else
    ind = 1:nsrc;
  end
  ind = setdiff(ind,i);
  %absolute Euclidean distance
  dt{i} = abs(tsplst(ind)-tsplst(i));
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
    % xtar = 1;
    [~,ind] = min(abs(x-xtar));
    cdftar(i) = cdfval(ind); 
  end
end
% median(cdftar)
xlabel(ax,'Distance between each source to all others (km or spl)');
ylabel(ax,'Empirical CDF');
if ~isempty(septran)
  text(ax,0.9,0.1,sprintf('Sources within %.1f-%.1f s',septran(1)/sps,septran(2)/sps),...
    'Units','normalized','HorizontalAlignment','right');
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
xlabel(ax,'Distance between each source to all others (km or spl)');
ylabel(ax,'Empirical PDF');
hold(ax,'off');

% keyboard
