function plt_deconpkvswfpk_rat(f,clppkhtwfrall,psrcamprsall,clnpkhtwfrall,nsrcamprsall,color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plt_deconpkvswfpk_rat(f,clppkhtwfrall,psrcamprsall,clnpkhtwfrall,nsrcamprsall,color)
%
% This function is to plot the histogram of ratio between the closest 
% (i.e. associated) deconvolved peak amp ratio and waveform peak amp 
% ratio (both pos and neg) so see how they fit. The expectation is, 
% they are close to each other. Combine sources from all burst windows 
% together.
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
  temp = clppkhtwfrall(:,i)./psrcamprsall(:,i);
  temp = temp(temp>0);
  histogram(ax,log10(temp),'FaceColor',color);
  plot(ax,[median(log10(temp)) median(log10(temp))],ax.YLim,'--',...
    'color',color,'linew',1);
%   text(ax,0.95,0.9,sprintf('med=%.2f',median(log10(temp))),'Units','normalized',...
%     'HorizontalAlignment','right');
  if i ==1
    xlabel(ax,'log_{10}{closest pos amp ratio 1/2 / pos waveform peak height ratio 1/2}');
    ylabel('Count');
  elseif i ==2
    xlabel(ax,'log_{10}{closest pos amp ratio 1/3 / pos waveform peak height ratio 1/3}');
  else
    xlabel(ax,'log_{10}{closest pos amp ratio 2/3 / pos waveform peak height ratio 2/3}');
  end
  xlim(ax,[-1 1]);

  ax=f.ax(3+i);
  hold(ax,'on');  
  grid(ax,'on');
  temp = clnpkhtwfrall(:,i)./nsrcamprsall(:,i);
  temp = temp(temp>0);
  histogram(ax,log10(temp),'FaceColor',color);
  plot(ax,[median(log10(temp)) median(log10(temp))],ax.YLim,'--',...
    'color',color,'linew',1);
%   text(ax,0.95,0.9,sprintf('med=%.2f',median(log10(temp))),'Units','normalized',...
%     'HorizontalAlignment','right');
  if i ==1
    xlabel(ax,'log_{10}{closest neg amp ratio 1/2 / neg waveform peak height ratio 1/2}');
    ylabel(ax,'Count');
  elseif i ==2
    xlabel(ax,'log_{10}{closest neg amp ratio 1/3 / neg waveform peak height ratio 1/3}');
  else
    xlabel(ax,'log_{10}{closest neg amp ratio 2/3 / neg waveform peak height ratio 2/3}');
  end
  xlim(ax,[-1 1]);

end


