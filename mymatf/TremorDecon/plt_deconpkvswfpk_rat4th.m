function f=plt_deconpkvswfpk_rat4th(f,clppkhtwfall,psrcampsall,clnpkhtwfall,nsrcampsall,color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_deconpkvswfpk_rat(f,clppkhtwfrall,psrcamprsall,clnpkhtwfrall,nsrcamprsall,color)
%
% Similar to 'plt_deconpkvswfpk_rat',it is to plot the HISTOGRAM of ratio between the closest 
% (i.e. associated) deconvolved peak amp and waveform peak amp 
% (both pos and neg) between all station pairs (12,13, and 23) 
% including 14, so see how they fit. The expectation is, 
% they are close to each other. Combine sources from all burst windows 
% together.
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
  temp = clppkhtwfall(:,i)./psrcampsall(:,i);
  temp = temp(temp>0);
  histogram(ax,log10(temp),'FaceColor',color,'Normalization','pdf');
  plot(ax,[median(log10(temp)) median(log10(temp))],ax.YLim,'--',...
    'color',color,'linew',1);
%   text(ax,0.95,0.9,sprintf('med=%.2f',median(log10(temp))),'Units','normalized',...
%     'HorizontalAlignment','right');
  if i ==1
    xlabel(ax,'log_{10}{closest pos amp / pos waveform peak height at sta 1}');
    ylabel('PDF');
  elseif i ==2
    xlabel(ax,'log_{10}{closest pos amp / pos waveform peak height at sta 2}');
  elseif i ==3
    xlabel(ax,'log_{10}{closest pos amp / pos waveform peak height at sta 3}');
  elseif i ==4
    xlabel(ax,'log_{10}{closest pos amp / pos waveform peak height at sta 4}');
  end
  xlim(ax,[-1 1]);

  ax=f.ax(4+i);
  hold(ax,'on');  
  grid(ax,'on');
  temp = clnpkhtwfall(:,i)./nsrcampsall(:,i);
  temp = temp(temp>0);
  histogram(ax,log10(temp),'FaceColor',color,'Normalization','pdf');
  plot(ax,[median(log10(temp)) median(log10(temp))],ax.YLim,'--',...
    'color',color,'linew',1);
%   text(ax,0.95,0.9,sprintf('med=%.2f',median(log10(temp))),'Units','normalized',...
%     'HorizontalAlignment','right');
  if i ==1
    xlabel(ax,'log_{10}{closest neg amp / neg waveform peak height at sta 1}');
    ylabel('PDF');
  elseif i ==2
    xlabel(ax,'log_{10}{closest neg amp / neg waveform peak height at sta 2}');
  elseif i ==3
    xlabel(ax,'log_{10}{closest neg amp / neg waveform peak height at sta 3}');
  elseif i ==4
    xlabel(ax,'log_{10}{closest neg amp / neg waveform peak height at sta 4}');
  end
  xlim(ax,[-1 1]);

end


