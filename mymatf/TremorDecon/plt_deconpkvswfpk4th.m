function f=plt_deconpkvswfpk4th(f,clppkhtwfall,psrcampsall,clnpkhtwfall,nsrcampsall,color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_deconpkvswfpk(f,clppkhtwfrall,psrcamprsall,clnpkhtwfrall,nsrcamprsall,color)
%
% This function is to plot the SCATTER between the closest (i.e. associated)
% deconvolved peak amp and waveform peak amp (both pos and neg) between all 
% station pairs (12,13, and 23) including 14 so see how they fit. The 
% expectation is, they are close to each other, ie., distribute along the 
% line with a slope of 1. Combine sources from all burst windows together.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/11/21
% Last modified date:   2022/11/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1: 4
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
%   axis(ax,[-1.5 1.5 -1.5 1.5]);
  plot(ax,[min([clppkhtwfall(:,i); psrcampsall(:,i)]), max([clppkhtwfall(:,i); psrcampsall(:,i)])],...
    [min([clppkhtwfall(:,i); psrcampsall(:,i)]), max([clppkhtwfall(:,i); psrcampsall(:,i)])],'r--');
  if i ==1
    ylabel(ax,'closest pos amp at sta 1');
    xlabel(ax,'pos waveform peak height at sta 1');
  elseif i ==2
    ylabel(ax,'closest pos amp at sta 2');
    xlabel(ax,'pos waveform peak height at sta 2');
  elseif i ==3 
    ylabel(ax,'closest pos amp at sta 3');
    xlabel(ax,'pos waveform peak height at sta 3');
  elseif i ==4 
    ylabel(ax,'closest pos amp at sta 4');
    xlabel(ax,'pos waveform peak height at sta 4');
  end
  
  ax=f.ax(4+i);
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
%   axis(ax,[-1.5 1.5 -1.5 1.5]);
  plot(ax,[min([clnpkhtwfall(:,i); nsrcampsall(:,i)]), max([clnpkhtwfall(:,i); nsrcampsall(:,i)])],...
    [min([clnpkhtwfall(:,i); nsrcampsall(:,i)]), max([clnpkhtwfall(:,i); nsrcampsall(:,i)])],'r--');
  if i ==1
    ylabel(ax,'closest neg amp at sta 1');
    xlabel(ax,'neg waveform peak height at sta 1');
  elseif i ==2
    ylabel(ax,'closest neg amp at sta 2');
    xlabel(ax,'neg waveform peak height at sta 2');
  elseif i ==3
    ylabel(ax,'closest neg amp at sta 3');
    xlabel(ax,'neg waveform peak height at sta 3');
  elseif i ==4
    ylabel(ax,'closest neg amp at sta 4');
    xlabel(ax,'neg waveform peak height at sta 4');
  end
end


