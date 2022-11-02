function f1 = plt_traceindetail(trace,stas,lwin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
% This is a simple function to plot 'trace' at different stations. This trace
% could be jusy signal waveform, or signal envelope, or sig-wlet CC. It only
% aims to plot the trace in a way that you can see the wiggles, and visually
% feel the coherence between them.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/10/27
% Last modified date:   2022/10/27 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

widin = 15;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
lsig = size(trace,1);
nrow = ceil(lsig/lwin); 
ncol = 1;
f1 = initfig(widin,htin,nrow,ncol);

pltxran = [0.06 0.96]; pltyran = [0.06 0.96];
pltxsep = 0.02; pltysep = 0.035;
optaxpos(f1,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);


ym = max(abs(trace(:)));
yran=1.2*[-ym ym];

nsta = size(trace,2);
color = jet(nsta);
for i = 1: nrow
  ax=f1.ax(i);
  hold(ax,'on');
  for j = 1: nsta
    plot(ax, trace(:,j)+(j-1)*ym/3, '-','Color',color(j,:)); 
  end
  xlim(ax,[(i-1)*lwin,i*lwin]); 
%   ylim(ax,yran);
  ax.Box='on'; 
  longticks(ax,8);
  if i~=nrow
%     nolabels(ax,1);
  else
    xlabel(ax,'Samples','fontsize',10);
    ylabel(ax,'Amplitude','fontsize',10);
    legend(ax,stas);
  end
  if i==1
  end

end  
  
  
  
  
  
  
  
  
