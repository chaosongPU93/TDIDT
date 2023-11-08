function [f,ampbincnt,fraci,dtarvlplt,ampplt] = ...
  plt_fracdifftime_NNm_syn_binamp(f,impplt,mmax,nsat,nrounds,label,sps,m,nbin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,ampbincnt,fraci,dtarvlplt,ampplt] = ...
%   plt_fracdifftime_NNm_syn_binamp(f,impplt,mmax,nsat,nrounds,label,sps,m,nbin)
%
% Similar purpose as 'plt_frac_difftime_NNm_syn_binamp', This function is to
% summarize all saturation levels, and trials of synthetics (either different
% source region size or noise levels), 
% plot only the fraction of diff time distribution of N & N-m source pairs, 
% binned by amp first, w/i some dtcut. 
% Amp here is the median amp of the 
% cluster defined by consecutive events within N & N-m source pairs.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/11/05
% Last modified date:   2023/11/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('m',1);
defval('nbin',5);

dtcut = 0.25*m+0.125;
if m < 3
  xran = [0 2];
else
  xran = [0 2*ceil(dtcut)];
end

nnsat = length(nsat);
color = jet(nrounds);

%%%loop for saturation level
for insat = 1: nnsat
  disp(nsat(insat));
  
  ax=f.ax(insat); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  
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
    
    ampcont = [];
    for i = 1: nsrc-m
      impcont = impi(i: i+m,:);
      ampcont(i,1) = median(mean(impcont(:,[2 4 6]),2));
    end
    ampplt = ampcont;
    
    %%%
    %%%%%%%%%% if bin by amp with a equal number
%     nbin = 5;
    [ampbin,indbin,n] = binxeqnum(ampplt,nbin);
    binwdt = 0.05;
    dtcut = 0.25*m+0.125;
    if m < 3
      xran = [0 2];
    else
      xran = [0 2*ceil(dtcut)];
    end
    nx = round(xran(2)/binwdt)+1;
    Nn = zeros(nx, nbin);
    edges = xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2;
    cnt = xran(1): binwdt: xran(2);
    fraci = [];
    for i = 1: nbin
      ind = indbin{i};
      dtarvli = dtarvlplt(ind);
      ampbincnt(i) = median(ampbin{i});
      N=histcounts(dtarvli/sps,edges,'normalization','count');
      Nn(:,i) = N./mode(n);
      fraci(i,iround) = sum(dtarvli/sps<=dtcut)/mode(n);
    end
    p(iround) = plot(ax,log(ampbincnt),fraci(:,iround),'-','Color',color(iround,:),'linew',1,...
      'marker','o','markersize',4,'markerfacec',color(iround,:));
%     if ~singleflag
%       label{iround} = sprintf('a/2=%.2f,b/2=%.2f',semia(iround),semib(iround));
%     else
%       label{iround} = sprintf('noise=%.1f',perctrial(iround));
%     end
    
  end
  if insat == 1
    xlabel(ax,'Median log_{10}{amp}');
    % ylabel(ax,'Median diff. arrival (s)');
    ylabel(ax,sprintf('Frac. of diff. arrival w/i %.3f s',dtcut));
    legend(ax,p,label);
  end
  text(ax,0.02,0.05,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','left');
  yran=[0 1];
  ylim(ax,yran);
  longticks(ax,2);
  hold(ax,'off');
    
end
