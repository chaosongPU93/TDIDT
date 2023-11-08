function [f,frac]=plt_frac_difftime_NNm_syn(f,impplt,mmax,nsat,nrounds,label,sps,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,frac]=plt_frac_difftime_NNm_syn(f,impplt,mmax,nsat,nrounds,label,sps,m)
%
% This function is to plot the fraction of diff time distribution of N & N-m
% source pairs, for a certain m. It summarizes different saturation, and
% different noise levels or source region sizes of synthetics. 
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
    
    frac(insat,iround) = sum(dtarvlplt/sps<=dtcut)/length(dtarvlplt);
  end
end

ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%%%loop for region size
for iround = 1: nrounds
  p(iround) = plot(ax,log10(nsat),frac(:,iround),'-o','markersize',4,'color',color(iround,:));
end
legend(ax,p,label);
xlabel(ax,'log_{10}(Saturation)');
ylabel(ax,sprintf('Frac. of diff. arrival w/i %.3f s',dtcut));
yran = [0 1];
ylim(ax,yran);
xlim(ax,[-0.5 2]);
