function [f,lgd]=plt_deconpk_rat_stat(f,nsat,label,msrcampr,madsrcampr,ref,sybsize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axall=plt_deconpk_rat_stat(axall,stat)
%
% This function is to plot the derived statistics from the amp ratio of
% sources 12, 13, 23 (and 14 if a 4th sta is involved), rather than the
% histogram or scatter like 'plt_deconpk_rat_comb.m'. The purpose is to
% ease the summary of variation in ratio wrt noise level, sat level, etc.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/10/18
% Last modified date:   2023/10/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('ref',[]);
defval('sybsize',[]);

nnsat = size(msrcampr,1);
ntrial = size(msrcampr,2);
% color = jet(ntrial);
color = gradientblue(ntrial);
symbol = ['o';'^';'s';'d'];
for itrial = 1: ntrial
  mip = cat(1,msrcampr{:,itrial});
  madip = cat(1,madsrcampr{:,itrial});
  nsta = size(mip,2);
  for i = 1: nsta
    ax=f.ax(i); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log';
    if i == 1
      if isempty(sybsize)
        p(itrial) = plot(ax,nsat,mip(:,i),'-','marker',symbol(i,:),'markersize',4,...
          'color',color(itrial,:),'MarkerEdgeColor','k');
      else
        p(itrial) = plot(ax,nsat,mip(:,i),'-','color',color(itrial,:),'linew',1);
        scatter(ax,nsat,mip(:,i),sybsize(:,i),color(itrial,:),...
          symbol(i,:),'filled','MarkerEdgeColor','k');
      end
      %       label{itrial} = sprintf('noise=%.1f',trial(itrial));
      ylabel(ax,'Median of amp ratio');
      xlabel(ax,'Saturation');
    else
      if isempty(sybsize)
        plot(ax,nsat,mip(:,i),'-','marker',symbol(i,:),'markersize',4,...
          'color',color(itrial,:),'MarkerEdgeColor','k');
      else
        plot(ax,nsat,mip(:,i),'-','color',color(itrial,:),'linew',1);
        scatter(ax,nsat,mip(:,i),sybsize(:,i),color(itrial,:),...
          symbol(i,:),'filled','MarkerEdgeColor','k');
      end
    end
    if itrial == ntrial
      if ~isempty(ref)
        mamprd = ref{1};
        p(ntrial+1) = plot(ax,nsat,mamprd(i)*ones(nnsat,1),'k--','linew',1.5);
        if length(ref)>2
          mamprn = ref{3};
          p(ntrial+2) = plot(ax,nsat,mamprn(i)*ones(nnsat,1),'r--','linew',1.5);
        end
      end
    end
%     yran = [-0.1 0.1];
    yran = log10([0.7 1.4]);
    ylim(ax,yran);
    ytks = 0.7:0.1:1.4;
    yticks(ax,log10(ytks));
    ytklbls=cell(length(ytks),1);
    for kk=1:size(ytklbls,1)
      ytklbls{kk}=num2str(ytks(kk));
    end
    yticklabels(ax,ytklbls);
    xticks(ax,nsat);
    longticks(ax,2);
  end
  
  for i = nsta+1: nsta+nsta
    ax=f.ax(i); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); ax.XScale='log';
    if isempty(sybsize)
      plot(ax,nsat,madip(:,i-nsta),'-','marker',symbol(i-nsta,:),'markersize',4,...
        'color',color(itrial,:),'MarkerEdgeColor','k');
    else
      plot(ax,nsat,madip(:,i-nsta),'-','color',color(itrial,:),'linew',1);
      scatter(ax,nsat,madip(:,i-nsta),sybsize(:,i-nsta),color(itrial,:),...
        symbol(i-nsta,:),'filled','MarkerEdgeColor','k');
    end
    if i == nsta+1
      ylabel(ax,'STD of amp ratio');
      xlabel(ax,'Saturation');
    end
    if itrial == ntrial
      if ~isempty(ref)
        madamprd = ref{2};
        plot(ax,nsat,madamprd(i-nsta)*ones(nnsat,1),'k--','linew',1.5);
        if length(ref)>2
          madamprn = ref{4};
          plot(ax,nsat,madamprn(i-nsta)*ones(nnsat,1),'r--','linew',1.5);
        end
      end
    end
%     yran = [0 0.4];
    yran = log10([1.2 2.2]);
    ylim(ax,yran);
    ytks = 1.2:0.2:2.2;
    yticks(ax,log10(ytks));
    ytklbls=cell(length(ytks),1);
    for kk=1:size(ytklbls,1)
      ytklbls{kk}=num2str(ytks(kk));
    end
    yticklabels(ax,ytklbls);
    xticks(ax,nsat);
    longticks(ax,2);
  end
  
end
nolabels(f.ax(2:nsta),3);
nolabels(f.ax(nsta+2:nsta*2),3);

lgd(1)=legend(f.ax(1),p,label,'NumColumns',2,'Location','best','fontsize',6);  %'Orientation','vertical'
set(lgd(1).BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
lgd(2)=legend(f.ax(nsta+1),p,label,'NumColumns',2,'Location','best','fontsize',6);  %'Orientation','vertical'
set(lgd(2).BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  %make background transparent
text(f.ax(1),0.98,0.1,'PGC/SSIB','HorizontalAlignment','right',...
  'Units','normalized');
text(f.ax(1+nsta),0.98,0.05,'PGC/SSIB','HorizontalAlignment','right',...
  'Units','normalized');
text(f.ax(2),0.98,0.05,'PGC/SILB','HorizontalAlignment','right',...
  'Units','normalized');
text(f.ax(2+nsta),0.98,0.05,'PGC/SILB','HorizontalAlignment','right',...
  'Units','normalized');
text(f.ax(3),0.98,0.05,'SSIB/SILB','HorizontalAlignment','right',...
  'Units','normalized');
text(f.ax(3+nsta),0.98,0.05,'SSIB/SILB','HorizontalAlignment','right',...
  'Units','normalized');
if nsta == 4
  text(f.ax(4),0.98,0.05,'PGC/KLNB','HorizontalAlignment','right',...
    'Units','normalized');
  text(f.ax(4+nsta),0.98,0.05,'PGC/KLNB','HorizontalAlignment','right',...
    'Units','normalized');
end
