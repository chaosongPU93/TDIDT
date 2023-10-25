function f=plt_deconpk_rat_stat(f,nsat,label,msrcampr,madsrcampr)
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
      p(itrial) = plot(ax,log10(nsat),mip(:,i),'-','marker',symbol(i,:),'markersize',4,...
        'color',color(itrial,:));
%       label{itrial} = sprintf('noise=%.1f',trial(itrial));
      ylabel(ax,'Mean of amp ratio');
      xlabel(ax,'log_{10}(Saturation)');
      title(ax,'PGC/SSIB');
    else
      plot(ax,log10(nsat),mip(:,i),'-','marker',symbol(i,:),'markersize',4,...
        'color',color(itrial,:));
      if i==2
        title(ax,'PGC/SILB');
      elseif i==3
        title(ax,'SSIB/SILB');
      elseif i==4
        title(ax,'PGC/KLNB');
      end
    end
    yran = [-0.1 0.1];
    ylim(ax,yran);
    longticks(ax,2);
  end
  
  for i = nsta+1: nsta+nsta
    ax=f.ax(i); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    plot(ax,log10(nsat),madip(:,i-nsta),'-','marker',symbol(i-nsta,:),'markersize',4,...
      'color',color(itrial,:));
    if i == nsta+1
      ylabel(ax,'MAD of amp ratio');
      xlabel(ax,'log_{10}(Saturation)');
      title(ax,'PGC/SSIB');
    elseif i==nsta+2
      title(ax,'PGC/SILB');
    elseif i==nsta+3
      title(ax,'SSIB/SILB');
    elseif i==nsta+4
      title(ax,'PGC/KLNB');
    end
    yran = [0 0.2];
    ylim(ax,yran);
    longticks(ax,2);
  end
end
legend(f.ax(1),p,label);

