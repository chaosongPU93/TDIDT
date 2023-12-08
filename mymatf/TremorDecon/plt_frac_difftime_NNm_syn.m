function [f,frac]=plt_frac_difftime_NNm_syn(f,impplt,nsat,nround,label,sps,m,ref)
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
defval('ref',[]);

nnsat = length(nsat);
color = jet(nround);
dtcut = 0.25*m+0.125;
  
%%%loop for saturation level
for insat = 1: nnsat
  disp(nsat(insat));
  
  %%%loop for region size or noise level
  for iround = 1: nround
    impi = impplt{insat,iround};
    nsrc(insat,iround) = size(impi,1);
    
    %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
    dtarvl = diffcustom(impi(:,1), m,'forward');
%     dtarvl = srcdistNtoNm(impi(:,1),impi(:,7:8),mmax);
    %   dtarvlnnsepall = [dtarvlnnsepall; dtarvl{m}];
    if isempty(dtarvl)
      continue
    end
    dtarvlplt = dtarvl;
    
    frac(insat,iround) = sum(dtarvlplt/sps<=dtcut)/length(dtarvlplt);
  end
end

ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
%%%loop for region size
for iround = 1: nround
%     p(iround) = plot(ax,log10(nsat),frac(:,iround),'-o','markersize',4,...
%       'color',color(iround,:),'MarkerEdgeColor','k');
  p(iround) = plot(ax,log10(nsat),frac(:,iround),'-','Color',color(iround,:),'linew',1);
  scatter(ax,log10(nsat),frac(:,iround),nsrc(:,iround)/50,color(iround,:),...
    'filled','MarkerEdgeColor','k');
end
if ~isempty(ref)
  fracd = ref(:,m)*ones(1,nnsat);
  p(nround+1) = plot(ax,log10(nsat),fracd,'k--','LineWidth',1);
end

legend(ax,p,label,'NumColumns',2,'Location','north');
xlabel(ax,'log_{10}(Saturation)');
ylabel(ax,sprintf('Fraction'));
yran = [0 1];
ylim(ax,yran);
xlim(ax,[-0.5 2]);
longticks(ax,2);
hold(ax,'off');

