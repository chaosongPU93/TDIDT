function [f,Nn,frac]=plt_difftime_NNm_syn(f,impplt,mmax,nsat,nrounds,label,sps,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,Nn,frac]=plt_difftime_NNm_syn(f,impplt,mmax,nsat,nrounds,label,sps,m)
%
% Similar purpose as 'plt_fracdifftime_NNm_mmax', This function is for the certain m,
% to plot the diff time distribution of srcs, binned by amp
% first, for each m smaller than mmax, of N & N-m source pairs. 
% Amp here is the median amp of the 
% cluster composed by consecutive events within N & N-m source pairs.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/11/05
% Last modified date:   2023/11/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('m',1);

if m < 3
  xran = [0 2];
else
  xran = [0 2*ceil(0.25*m+0.125)];
end

nnsat = length(nsat);
color = jet(nrounds);
dtcut = 0.25*m+0.125;
 
%%%loop for saturation level
for insat = 1: nnsat
  disp(nsat(insat));
  
  ax=f.ax(insat); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  binwdt = 0.05;
  if m == 1
    yran = [0 0.2];
  elseif m == 2
    yran = [0 0.06];
  end
  nx = round(xran(2)/binwdt)+1;
  edges = xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2;
  cnt = xran(1): binwdt: xran(2);
  
  %%%loop for region size or noise level
  for iround = 1: nrounds
    impi = impplt{insat,iround};
    nsrc = size(impi,1);
    
    %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
    dtarvl = srcdistNtoNm(impi(:,1),impi(:,7:8),mmax);
    %   dtarvlnnsepall = [dtarvlnnsepall; dtarvl{m}];
    if isempty(dtarvl{m})
      continue
    end
    dtarvlplt = dtarvl{m};
    
    N=histcounts(dtarvlplt/sps,edges,'normalization','count');
    Nn = N./length(dtarvlplt);
%     Nn = N./denom;
    frac(insat,iround) = sum(dtarvlplt/sps<=dtcut)/length(dtarvlplt);
    p(iround) = plot(ax,cnt,Nn,'color',color(iround,:),'LineWidth',1);
    % if ~singleflag
    %   label{iround} = sprintf('a/2=%.2f,b/2=%.2f',semia(iround),semib(iround));
    % else
    %   label{iround} = sprintf('noise=%.1f',perctrial(iround));
    % end

  end
  %   text(ax,0.95,0.85,sprintf('Fraction w/i %.3f s: %.2f',dtcut,frac),'HorizontalAlignment','right',...
  %     'Units','normalized','FontSize',12);
  if insat == 1
    ylabel(ax,'Normalized count');
    xlabel(ax,sprintf('Diff. arrival between sources N and N-%d (s)',m));
    legend(ax,p,label);
  end
  text(ax,0.98,0.2,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','right');  
  xlim(ax,xran);
  ylim(ax,yran);
  longticks(ax,2);
  hold(ax,'off');
end