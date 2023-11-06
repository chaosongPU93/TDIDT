function [f,Nn,ampbincnt,fraci]=plt_frac_difftime_NNm_syn_binamp(f,impplt,mmax,insat,iround,sps,m,nbin)
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
defval('nbin',5);

dtcut = 0.25*m+0.125;
if m < 3
  xran = [0 2];
else
  xran = [0 2*ceil(dtcut)];
end

% nnsat = length(nsat);
% color = jet(nrounds);
  
%%%first bin source by amp, then plot diff arrival time for N and N-1 for each amp bin
%%%loop for saturation level
% for insat = 1: nnsat
  % disp(nsat(insat));
  %%%loop for region size or noise level
  % for iround = 1: nrounds
    impi = impplt{insat,iround};
    nsrc = size(impi,1);
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
    
    ampcont = [];
    for i = 1: nsrc-m
      impcont = impi(i: i+m,:);
      ampcont(i,1) = median(mean(impcont(:,[2 4 6]),2));
    end
    ampplt = ampcont;
    
    %%%
    ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    %%%%%%%%%% if bin by amp with a equal number
%     nbin = 5;
    [ampbin,indbin,n] = binxeqnum(ampplt,nbin);
    color = jet(nbin);
    binwdt = 0.05;
    nx = round(xran(2)/binwdt)+1;
    Nn = zeros(nx, nbin);
    edges = xran(1)-binwdt/2: binwdt: xran(2)+binwdt/2;
    cnt = xran(1): binwdt: xran(2);
    label = []; p = []; fraci = [];
    for i = 1: nbin
      ind = indbin{i};
      dtarvli = dtarvlplt(ind);
      ampbincnt(i) = median(ampbin{i});
      N=histcounts(dtarvli/sps,edges,'normalization','count');
      Nn(:,i) = N./mode(n);
      fraci(i) = sum(dtarvli/sps<=dtcut)/mode(n);
      p(i)=plot(ax,cnt,Nn(:,i),'color',color(i,:),'LineWidth',1);
      label{i} = sprintf('amp of %.1f',ampbincnt(i));
      % label{i} = sprintf('%d/%dth amp',i,nbin);
      % keyboard
    end
    %%%%%%%%%% if bin by amp with a equal number
    Nnm = mean(Nn,2);
    p(nbin+1)=plot(ax,cnt,Nnm,'k-','LineWidth',1.5);
    label{nbin+1} = sprintf('mean');  %median
    legend(ax,p,label);
    ylabel(ax,'Normalized count');
    xlabel(ax,sprintf('Diff. arrival between sources N and N-%d (s)',m));
    xlim(ax,xran);
    if m == 1
      yran = [0 0.2];
    elseif m == 2
      yran = [0 0.06];
    else
      yran = ax.YLim;
    end
    ylim(ax,yran);
    longticks(ax,2);
    hold(ax,'off');
%     title(ax,'Data');
    
    ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    yran=[0 1];
    plot(ax,log(ampbincnt),fraci,'k-','linew',1,'marker','o','markersize',4,'markerfacec','k');
    xlabel(ax,'Median log_{10}{amp}');
    % ylabel(ax,'Median diff. arrival (s)');
    ylabel(ax,sprintf('Frac. of diff. arrival w/i %.3f s',dtcut));
    ylim(ax,yran);
%     title(ax,'Data');
    
  % end
  
% end