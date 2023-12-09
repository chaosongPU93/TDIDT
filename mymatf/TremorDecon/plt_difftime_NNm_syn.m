function [f,Nn,frac]=plt_difftime_NNm_syn(f,impplt,nsat,nround,label,sps,m,ref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,Nn,frac]=plt_difftime_NNm_syn(f,impplt,mmax,nsat,nrounds,label,sps,m)
%
% This function is to plot the diff time distribution of N & N-m source 
% pairs, for a certain m. It summarizes different saturation, and different
% noise levels or source region sizes of synthetics. 
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/11/05
% Last modified date:   2023/11/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('m',1);
defval('ref',[]);

if m < 3
  xran = [0 2];
else
  xran = [0 2*ceil(0.25*m+0.125)];
end
binwdt = 0.05;
nx = round(xran(2)/binwdt)+1;
edges = xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2;
cnt = xran(1): binwdt: xran(2);
  
nnsat = length(nsat);
color = jet(nround);
dtcut = 0.25*m+0.125;
 
%%%loop for saturation level
for insat = 1: nnsat
  disp(nsat(insat));
  
  ax=f.ax(insat); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  if m == 1
    yran = [0 0.15];
  elseif m == 2
    yran = [0 0.06];
  end
  %%%loop for region size or noise level
  for iround = 1: nround
    impi = impplt{insat,iround};
    nsrc = size(impi,1);
    
    %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
    dtarvl = diffcustom(impi(:,1), m,'forward');
%     dtarvl = srcdistNtoNm(impi(:,1),impi(:,7:8),mmax);
    %   dtarvlnnsepall = [dtarvlnnsepall; dtarvl{m}];
    if isempty(dtarvl)
      continue
    end
    dtarvlplt = dtarvl;
    
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
  if ~isempty(ref)
    Nnd = ref(:,m);
    p(nround+1) = plot(ax,cnt,Nnd(1:length(cnt)),'k--','LineWidth',1);
  end
  %   text(ax,0.95,0.85,sprintf('Fraction w/i %.3f s: %.2f',dtcut,frac),'HorizontalAlignment','right',...
  %     'Units','normalized','FontSize',12);
  text(ax,0.98,0.2,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
    'HorizontalAlignment','right');  
  xlim(ax,xran);
  ylim(ax,yran);
  longticks(ax,2);
  hold(ax,'off');
end

ylabel(f.ax(5),'Normalized count');
xlabel(f.ax(5),sprintf('Arrival time difference N to N-%d (s)',m));
legend(f.ax(1),p,label);



