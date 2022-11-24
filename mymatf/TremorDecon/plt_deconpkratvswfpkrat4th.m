function f=plt_deconpkratvswfpkrat4th(f,clppkhtwfall,psrcampsall,clnpkhtwfall,nsrcampsall)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=plt_deconpkratvswfpkrat4th(f,clppkhtwfrall,psrcamprsall,clnpkhtwfrall,nsrcamprsall,color)
%
% Similar to 'plt_deconpkratvswfpkrat', this function is to plot the
% histogram of the closet waveform peak ratio between each station pair
% and the deconvolved source amp ratio between each station pair,
% including pair 14. Note that the ratio is to the same quantity at the 
% same station pair. This is to know if the distribution of waveform peaks
% are wider or narrower than the deconvolved sources, may inspire us on
% if sources with an amp ratio too deviated from the center should be 
% reasonably thrown away, if it turns out that source amp ratio is wider
% in distribution than the waveform.
% 
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/11/21
% Last modified date:   2022/11/21 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ppkr = [clppkhtwfall(:,1)./clppkhtwfall(:,2) clppkhtwfall(:,1)./clppkhtwfall(:,2) ...
  clppkhtwfall(:,2)/clppkhtwfall(:,3) clppkhtwfall(:,1)./clppkhtwfall(:,4)];
npkr = [clnpkhtwfall(:,1)./clnpkhtwfall(:,2) clnpkhtwfall(:,1)./clnpkhtwfall(:,2) ...
  clnpkhtwfall(:,2)/clnpkhtwfall(:,3) clnpkhtwfall(:,1)./clnpkhtwfall(:,4)];
psrcr = [psrcampsall(:,1)./psrcampsall(:,2) psrcampsall(:,1)./psrcampsall(:,2) ...
  psrcampsall(:,2)/psrcampsall(:,3) psrcampsall(:,1)./psrcampsall(:,4)];
nsrcr = [nsrcampsall(:,1)./nsrcampsall(:,2) nsrcampsall(:,1)./nsrcampsall(:,2) ...
  nsrcampsall(:,2)/nsrcampsall(:,3) nsrcampsall(:,1)./nsrcampsall(:,4)];

for i = 1: 4
  ax=f.ax(i);
  hold(ax,'on');  
  grid(ax,'on');
  temp = ppkr(:,i);
  temp = temp(temp>0);
  p(1)=histogram(ax,log10(temp),'FaceColor','k','Normalization','pdf');
  plot(ax,[median(log10(temp)) median(log10(temp))],ax.YLim,'--',...
    'color','k','linew',1);
  temp = psrcr(:,i);
  temp = temp(temp>0);
  p(2)=histogram(ax,log10(temp),'FaceColor','r','Normalization','pdf');
  plot(ax,[median(log10(temp)) median(log10(temp))],ax.YLim,'--',...
    'color','r','linew',1);
  if i ==1
    ylabel('PDF');
  end
  xlim(ax,[-1 1]);
  if i == 1
    legend(ax,[p1 p2],'closest pos waveform peaks','scaled pos deconvolved peaks');
  end
  
  ax=f.ax(4+i);
  hold(ax,'on');  
  grid(ax,'on');
  temp = npkr(:,i);
  temp = temp(temp>0);
  p(1)=histogram(ax,log10(temp),'FaceColor','k','Normalization','pdf');
  plot(ax,[median(log10(temp)) median(log10(temp))],ax.YLim,'--',...
    'color','k','linew',1);
  temp = nsrcr(:,i);
  temp = temp(temp>0);
  p(2)=histogram(ax,log10(temp),'FaceColor','r','Normalization','pdf');
  plot(ax,[median(log10(temp)) median(log10(temp))],ax.YLim,'--',...
    'color','r','linew',1);
  if i ==1
    xlabel(ax,'log_{10}{amp ratio 1/2}');
    ylabel('PDF');
  elseif i ==2
    xlabel(ax,'log_{10}{amp ratio 1/2}');
  elseif i ==3
    xlabel(ax,'log_{10}{amp ratio 2/3}');
  elseif i ==4
    xlabel(ax,'log_{10}{amp ratio 1/4}');
  end
  xlim(ax,[-1 1]);
  if i == 1
    legend(ax,[p1 p2],'closest neg waveform peaks','scaled neg deconvolved peaks');
  end

end