function [fracsrc2all,dfracsrc2all,mmaxzero,f]=frac_uniqevt_incluster(nbst,imp,nsrc,mmax,sps,pltflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [fracsrc2all,dfracsrc2all,f]=frac_uniqevt_incluster(nbst,imp,nsrc,mmax,sps)
%
% This function is to compute the fraction of unique events within the 
% cluster defined by consecutive events within N & N-m source pairs. 
% The diff time of eligible N & N-m source pairs needs be smaller than
% a 'dtcut = 0.25*m+0.125'. Returns the inclusive fractions 'fracsrc2all',
% and exclusive fractions 'dfracsrc2all', and a figure handle. 
% Deal with data and synthetic noise only, NOT for synthetics.
% Deal with a series of m < mmax.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/11/03
% Last modified date:   2023/11/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('pltflag',1);

fracsrc2all = [];
for m = 1:mmax
  
  dtcut = 0.25*m+0.125;
  nsrcsep = nsrc-m;
  nsrcsep(nsrcsep<0) = 0;
  ndcutpair = zeros(nbst, 1);
  ndcutsrc = zeros(nbst, 1);
  ndcutsrc2 = zeros(nbst, 1);
  imppair = cell(nbst, 1);
  imppairuni = cell(nbst, 1);
  impcont = cell(nbst, 1);
  impcontuni = cell(nbst, 1);
  indcont = cell(nbst, 1);
  indcontuni = cell(nbst, 1);
  
  for i = 1: nbst
    if nsrc(i) == 0
      continue
    end
    ist = sum(nsrc(1:i-1))+1;
    ied = ist+nsrc(i)-1;
    impi = imp(ist:ied,:);
    %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
    dtarvl = srcdistNtoNm(impi(:,1),impi(:,7:8),mmax);
    %   dtarvlnnsepall = [dtarvlnnsepall; dtarvl{nsep}];
    if isempty(dtarvl{m})
      continue
    end
    impbf = impi(1:end-m,:);
    impaf = impi(1+m:end,:);
    if ~isequal(size(impbf,1),length(dtarvl{m}))
      disp('Check');
    end
    ind = find(dtarvl{m}/sps <= dtcut);
    ndcutpair(i,1) = length(ind);
    tmp = [];
    for j = 1: length(ind)
      tmp = [tmp; impbf(ind(j),:); impaf(ind(j),:)];
    end
    imppair{i,1} = tmp;
    % impsort = sortrows(imppair,1);
    imppairuni{i,1} = unique(tmp,'rows','stable');
    % impall = [impall; impuni];
    ndcutsrc(i,1) = size(imppairuni{i,1},1);
    
    tmp = [];
    for j = 1: length(ind)
      tmp = [tmp; impi(ind(j):ind(j)+m,:)];
    end
    impcont{i,1} = tmp;
    impcontuni{i,1} = unique(tmp,'rows','stable');
    % impall = [impall; impuni];
    ndcutsrc2(i,1) = size(impcontuni{i,1},1);
    
    tmp = [];
    for j = 1: length(ind)
      tmp = [tmp; reshape(ind(j):ind(j)+m, [], 1)];
    end
    indcont{i,1} = tmp;
    indcontuni{i,1} = unique(tmp,'rows','stable');
  end
  
%   imppairunia = cat(1,imppairuni{:});
%   ndcutsrcall = sum(ndcutsrc);
%   fracsrcall = ndcutsrcall/sum(nsrc)*100;
  
  impcontunia = cat(1,impcontuni{:});
  ndcutsrc2all(m,1) = sum(ndcutsrc2);
  fracsrc2all(m,1) = ndcutsrc2all(m,1)/sum(nsrc)*100;
  
%   ndcutpairall(m,1) = sum(ndcutpair);
%   fracpairall = ndcutpairall(m,1)/sum(nsrcsep)*100;
  
end

mmaxzero = find(ndcutsrc2all==0,1);
if isempty(mmaxzero)  %if the trial mmax is not big enough to exhaust the list
  mmaxzero = mmax;
  warning('mmax is not big enough, better increase it');
end

% for m = 1:mmaxzero
%   dtcut = 0.25*m+0.125;
%   fprintf('%d clusters of %d consecutive events (%d/%d unique) w/i %.3f s \n',...
%     ndcutpairall(m,1), m+1, ndcutsrc2all(m,1), sum(nsrc), dtcut);
% end


fracsrc2all = [100; fracsrc2all];
fprintf('%.3f \n',fracsrc2all);
dfracsrc2all = fracsrc2all(1:end-1) - fracsrc2all(2:end);
fprintf('%.3f \n',dfracsrc2all);

% fracsrc2all = fracsrc2all(1:mmaxzero+1);
% dfracsrc2all = dfracsrc2all(1:mmaxzero);

if pltflag
  f=initfig;
  ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  plot(ax,0:1:mmaxzero,fracsrc2all(1:mmaxzero+1),'ko-','linew',1,'markersize',4);
  plot(ax,1:mmaxzero,dfracsrc2all(1:mmaxzero),'ro-','linew',1,'markersize',4);
  xlabel(ax,'m');
  ylabel(ax,'Frac of srcs');
  title(ax,'Frac of srcs in cluster of N & N-m whose diff time w/i 0.25*m+0.125 s');
  legend(ax,'inclusive','exclusive');
  xlim(ax,[0 mmax]);
else
  f=[];
end

