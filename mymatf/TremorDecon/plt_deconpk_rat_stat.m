function f=plt_deconpk_rat_stat(f,nsat,label,msrcampr,madsrcampr,ref,sybsize)
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
color = jet(ntrial);
symbol = ['o';'^';'s';'d'];
for itrial = 1: ntrial
  mip = cat(1,msrcampr{:,itrial});
  madip = cat(1,madsrcampr{:,itrial});  
  nsta = size(mip,2);
  for i = 1: nsta
    ax=f.ax(i); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    if i == 1
      if isempty(sybsize)
        p(itrial) = plot(ax,log10(nsat),mip(:,i),'-','marker',symbol(i,:),'markersize',4,...
          'color',color(itrial,:));
      else
        p(itrial) = plot(ax,log10(nsat),mip(:,i),'-','color',color(itrial,:),'linew',1);
        scatter(ax,log10(nsat),mip(:,i),sybsize(:,i),color(itrial,:),...
          symbol(i,:),'filled');
      end
%       label{itrial} = sprintf('noise=%.1f',trial(itrial));
      ylabel(ax,'Median of log_{10}(amp ratio)');
      xlabel(ax,'log_{10}(Saturation)');
    else
      if isempty(sybsize)
        plot(ax,log10(nsat),mip(:,i),'-','marker',symbol(i,:),'markersize',4,...
          'color',color(itrial,:));
      else
        plot(ax,log10(nsat),mip(:,i),'-','color',color(itrial,:),'linew',1);
        scatter(ax,log10(nsat),mip(:,i),sybsize(:,i),color(itrial,:),...
          symbol(i,:),'filled');
      end  
    end
    if itrial == ntrial
      if ~isempty(ref)
        mamprd = ref{1};
        p(ntrial+1) = plot(ax,log10(nsat),mamprd(i)*ones(nnsat,1),'k--');
      end
    end
    yran = [-0.1 0.1];
    ylim(ax,yran);
    longticks(ax,2);
  end
  
  for i = nsta+1: nsta+nsta
    ax=f.ax(i); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    if isempty(sybsize)
      plot(ax,log10(nsat),madip(:,i-nsta),'-','marker',symbol(i-nsta,:),'markersize',4,...
        'color',color(itrial,:));
    else
      plot(ax,log10(nsat),madip(:,i-nsta),'-','color',color(itrial,:),'linew',1);
      scatter(ax,log10(nsat),madip(:,i-nsta),sybsize(:,i-nsta),color(itrial,:),...
        symbol(i-nsta,:),'filled');
    end
    if i == nsta+1
      ylabel(ax,'MAD of log_{10}(amp ratio)');
      xlabel(ax,'log_{10}(Saturation)');
    end
    if itrial == ntrial
      if ~isempty(ref)
        madamprd = ref{2};
        plot(ax,log10(nsat),madamprd(i-nsta)*ones(nnsat,1),'k--');
      end
    end
    yran = [0 0.3];
    ylim(ax,yran);
    longticks(ax,2);
  end
end
legend(f.ax(1),p,label,'NumColumns',2,'Location','north');
legend(f.ax(4),p,label,'NumColumns',2,'Location','north');
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
