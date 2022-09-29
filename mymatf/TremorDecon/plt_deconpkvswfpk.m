function plt_deconpkvswfpk(f,clppkhtwfrall,psrcamprsall,clnpkhtwfrall,nsrcamprsall,color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plt_deconpkvswfpk(f,clppkhtwfrall,psrcamprsall,clnpkhtwfrall,nsrcamprsall,color)
%
% This function is to plot the scatter between the closest (i.e. associated)
% deconvolved peak amp ratio and waveform peak amp ratio (both pos and neg)
% so see how they fit. The expectation is, they are close to each other. 
% Combine sources from all burst windows together.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/09/28
% Last modified date:   2022/09/28 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1: 3
  ax=f.ax(i);
  hold(ax,'on');  
  grid(ax,'on');
  scatter(ax,log10(clppkhtwfrall(:,i)),log10(psrcamprsall(:,i)),5,color,'o');
  if isequal(color,'k')
    scatter(ax,median(log10(clppkhtwfrall(:,i))),median(log10(psrcamprsall(:,i))),40,'b','^','filled',...
      'markeredgecolor','k','linew',1);
  else
    scatter(ax,median(log10(clppkhtwfrall(:,i))),median(log10(psrcamprsall(:,i))),40,color,'^','filled',...
      'markeredgecolor','k','linew',1);
  end
  axis(ax,'equal');
  axis(ax,[-1.5 1.5 -1.5 1.5]);
  plot(ax,[-1.5 1.5 ],[-1.5 1.5],'r--');
  if i ==1
    ylabel(ax,'log_{10}{closest pos amp ratio 1/2}');
    xlabel(ax,'log_{10}{pos waveform peak height ratio 1/2}');
  elseif i ==2
    ylabel(ax,'log_{10}{closest pos amp ratio 1/3}');
    xlabel(ax,'log_{10}{pos waveform peak height ratio 1/3}');
  else
    ylabel(ax,'log_{10}{closest pos amp ratio 2/3}');
    xlabel(ax,'log_{10}{pos waveform peak height ratio 2/3}');
  end
  
  ax=f.ax(3+i);
  hold(ax,'on');  
  grid(ax,'on');
  scatter(ax,log10(clnpkhtwfrall(:,i)),log10(nsrcamprsall(:,i)),5,color,'o');
  if isequal(color,'k')
    scatter(ax,median(log10(clnpkhtwfrall(:,i))),median(log10(clnpkhtwfrall(:,i))),40,'b','^','filled',...
      'markeredgecolor','k','linew',1);
  else
    scatter(ax,median(log10(clnpkhtwfrall(:,i))),median(log10(nsrcamprsall(:,i))),40,color,'^','filled',...
      'markeredgecolor','k','linew',1);
  end
  axis(ax,'equal');
  axis(ax,[-1.5 1.5 -1.5 1.5]);
  plot(ax,[-1.5 1.5 ],[-1.5 1.5],'r--');
  if i ==1
    ylabel(ax,'log_{10}{closest neg amp ratio 1/2}');
    xlabel(ax,'log_{10}{neg waveform peak height ratio 1/2}');
  elseif i ==2
    ylabel(ax,'log_{10}{closest neg amp ratio 1/3}');
    xlabel(ax,'log_{10}{neg waveform peak height ratio 1/3}');
  else
    ylabel(ax,'log_{10}{closest neg amp ratio 2/3}');
    xlabel(ax,'log_{10}{neg waveform peak height ratio 2/3}');
  end

end


