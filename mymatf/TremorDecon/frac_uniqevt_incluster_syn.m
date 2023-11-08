function [f1,f2,frac,dfrac]=frac_uniqevt_incluster_syn(f1,f2,impplt,mmax,nsat,nrounds,label,sps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [fracsrc2all,dfracsrc2all,f]=frac_uniqevt_incluster_syn(nbst,imp,nsrc,mmax,sps)
%
% This function is to compute the fraction of unique events within the
% cluster composed by consecutive events within N & N-m source pairs.
% The diff time of eligible N & N-m source pairs needs be smaller than
% a 'dtcut = 0.25*m+0.125'. Will create 2 separate plots of inclusive and 
% exclusive fractions. 
% Summarize all saturation levels, and trials of synthetics (either different
% source region size or noise levels)
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/11/03
% Last modified date:   2023/11/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nnsat = length(nsat);
color = jet(nrounds);

%%%loop for saturation level
for insat = 1: nnsat
  disp(nsat(insat));
  
  %%%loop for region size or noise level
  for iround = 1: nrounds
    impi = impplt{insat,iround};
    
    for m = 1:mmax
      
      dcut = 0.25*m+0.125;
      impi = impplt{insat,iround};
      nsrc = size(impi,1);
      nsrcsep = nsrc-m;
      nsrcsep(nsrcsep<0) = 0;
      %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
      dtarvl = srcdistNtoNm(impi(:,1),impi(:,7:8),mmax);
      %   dtarvlnnsepall = [dtarvlnnsepall; dtarvl{m}];
      if isempty(dtarvl{m})
        Nn=[];
        ampbincnt=[];
        fraci=[];
        return
      end
      dtarvlplt = dtarvl{m};
      
      impbf = impi(1:end-m,:);
      impaf = impi(1+m:end,:);
      if ~isequal(size(impbf,1),length(dtarvl{m}))
        disp('Check');
      end
      ind = find(dtarvl{m}/sps <= dcut);
      ndcutpair = length(ind);
      tmp = [];
      for j = 1: length(ind)
        tmp = [tmp; impbf(ind(j),:); impaf(ind(j),:)];
      end
      imppair = tmp;
      % impsort = sortrows(imppair,1);
      imppairuni = unique(tmp,'rows','stable');
      % impall = [impall; impuni];
      ndcutsrc = size(imppairuni,1);
      
      tmp = [];
      for j = 1: length(ind)
        tmp = [tmp; impi(ind(j):ind(j)+m,:)];
      end
      impcont = tmp;
      impcontuni = unique(tmp,'rows','stable');
      % impall = [impall; impuni];
      ndcutsrc2 = size(impcontuni,1);
      
      tmp = [];
      for j = 1: length(ind)
        tmp = [tmp; reshape(ind(j):ind(j)+m, [], 1)];
      end
      indcont = tmp;
      indcontuni = unique(tmp,'rows','stable');
      
      %   imppairunia = cat(1,imppairuni{:});
      %   ndcutsrcall = sum(ndcutsrc);
      %   fracsrcall = ndcutsrcall/sum(nsrc)*100;
      
      %   impcontunia = cat(1,impcontuni{:});
      ndcutsrc2all(m,1) = sum(ndcutsrc2);
      fracsrc2all(m,1) = ndcutsrc2all(m,1)/sum(nsrc)*100;
      
      ndcutpairall(m,1) = sum(ndcutpair);
      fracpairall = ndcutpairall(m,1)/sum(nsrcsep)*100;
      
    end
    
%     for m = 1:mmax
%       dtcut = 0.25*m+0.125;
%       fprintf('%d clusters of %d consecutive events (%d/%d unique) w/i %.3f s \n',...
%         ndcutpairall(m,1), m+1, ndcutsrc2all(m,1), sum(nsrc), dtcut);
%     end
    
    frac(:,iround) = [100; fracsrc2all];
%     fprintf('%.3f \n',frac(:,iround));
    dfrac(:,iround) = frac(1:end-1,iround) - frac(2:end,iround);
%     fprintf('%.3f \n',dfrac(:,iround));
    
    ax=f1.ax(insat);
    hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    p1(iround) = plot(ax,0:1:mmax,frac(:,iround),'o-','linew',1,'markersize',4,'color',color(iround,:));
    ylim(ax,[0 100]);
    % title(ax,'Frac of srcs in cluster of N & N-m whose diff time w/i 0.25*m+0.125 s');
  %   legend(ax,'inclusive','exclusive');

    ax=f2.ax(insat);
    hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    p2(iround) = plot(ax,1:mmax,dfrac(:,iround),'o-','linew',1,'markersize',4,'color',color(iround,:));
    ylim(ax,[0 100]);
  end
  text(f1.ax(insat),0.98,0.3,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','right');
  text(f2.ax(insat),0.98,0.3,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','right');

end
legend(f1.ax(1),p1,label);
xlabel(f1.ax(1),'m');
ylabel(f1.ax(1),'Frac of srcs');
legend(f2.ax(1),p2,label);
xlabel(f2.ax(1),'m');
ylabel(f2.ax(1),'Frac of srcs');


