function plt_deconpkvswfpk(f,clppkhtwfall,psrcampsall,clnpkhtwfall,nsrcampsall,color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plt_deconpkvswfpk(f,clppkhtwfrall,psrcamprsall,clnpkhtwfrall,nsrcamprsall,color)
%
% This function is to plot the scatter between the closest (i.e. associated)
% deconvolved peak amp and waveform peak amp (both pos and neg)
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
  scatter(ax,clppkhtwfall(:,i),psrcampsall(:,i),5,color,'o');
  if isequal(color,'k')
    scatter(ax,median(clppkhtwfall(:,i)),median(psrcampsall(:,i)),40,'b','^','filled',...
      'markeredgecolor','k','linew',1);
  else
    scatter(ax,median(clppkhtwfall(:,i)),median(psrcampsall(:,i)),40,color,'^','filled',...
      'markeredgecolor','k','linew',1);
  end
  axis(ax,'equal');
  axis(ax,[-1.5 1.5 -1.5 1.5]);
  plot(ax,[-1.5 1.5 ],[-1.5 1.5],'r--');
  if i ==1
    ylabel(ax,'closest pos amp at sta 1');
    xlabel(ax,'pos waveform peak height at sta 1');
  elseif i ==2
    ylabel(ax,'closest pos amp at sta 2');
    xlabel(ax,'pos waveform peak height at sta 2');
  else
    ylabel(ax,'closest pos amp at sta 3');
    xlabel(ax,'pos waveform peak height at sta 3');
  end
  
  ax=f.ax(3+i);
  hold(ax,'on');  
  grid(ax,'on');
  scatter(ax,clnpkhtwfall(:,i),nsrcampsall(:,i),5,color,'o');
  if isequal(color,'k')
    scatter(ax,median(clnpkhtwfall(:,i)),median(clnpkhtwfall(:,i)),40,'b','^','filled',...
      'markeredgecolor','k','linew',1);
  else
    scatter(ax,median(clnpkhtwfall(:,i)),median(nsrcampsall(:,i)),40,color,'^','filled',...
      'markeredgecolor','k','linew',1);
  end
  axis(ax,'equal');
  axis(ax,[-1.5 1.5 -1.5 1.5]);
  plot(ax,[-1.5 1.5 ],[-1.5 1.5],'r--');
  if i ==1
    ylabel(ax,'closest neg amp at sta 1');
    xlabel(ax,'neg waveform peak height at sta 1');
  elseif i ==2
    ylabel(ax,'closest neg amp at sta 2');
    xlabel(ax,'neg waveform peak height at sta 2');
  else
    ylabel(ax,'closest neg amp at sta 3');
    xlabel(ax,'neg waveform peak height at sta 3');
  end
end


