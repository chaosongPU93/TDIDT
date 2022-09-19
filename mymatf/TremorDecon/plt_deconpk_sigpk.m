function [f,pk,pkhgt] = plt_deconpk_sigpk(sigsta,pkindep,indremove,fpltremv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_deconpk_sigpk(sigsta,pkindep,indremove,fpltremv)
%
% plot the positive peaks indicated by the grouped triplets, see if they 
% indeed match the peaks of the signal at each station, could be helpful
% to look in detail which triplets are minor
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/06/21
% Last modified date:   2022/06/21 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('indremove',[]);   % default input for indice for removing is empty
defval('fpltremv',0);   % default is do not plot the removed ones if any

f.fig = figure;
f.fig.Renderer = 'painters';
widin = 12;  % maximum width allowed is 8.5 inches
htin = 6;   % maximum height allowed is 11 inches
[scrsz, resol] = pixelperinch(1);
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);
nrow = size(sigsta,2); 
ncol = 1;
for isub = 1:nrow*ncol
  f.ax(isub) = subplot(nrow,ncol,isub);
end
pltxran = [0.1 0.9]; pltyran = [0.15 0.9];
pltxsep = 0.02; pltysep = 0.05;
%get the locations for each axis
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

color = ['r';'b';'k'];

ym = max(abs(sigsta(:)));
yran=1.2*[-ym ym];

lsig = size(sigsta,1); 

pkhgt = cell(size(sigsta,2),1);
pk = cell(size(sigsta,2),1);
for i = 1: nrow
  ax=f.ax(i);
  hold(ax,'on');
  plot(ax, sigsta(:,i), '-','Color',color(i,:)); 
  xlim(ax,[0,lsig]); ylim(ax,yran);
  [pkhgt{i}, pk{i}] = findpeaks(sigsta(:,i));
  scatter(ax,pk{i},pkhgt{i},8,color(i,:));
  median(diff(pk{i}))
  if ~isempty(indremove)
    for ii = 1: size(pkindep,1)
      if ~ismember(ii,indremove)
        plot(ax,[pkindep(ii,(i-1)*2+1) pkindep(ii,(i-1)*2+1)], ax.YLim, '--','Color',color(i,:));
        text(ax,pkindep(ii,(i-1)*2+1),yran(2)-0.1*range(yran),num2str(ii),'FontSize',6,...
          'HorizontalAlignment','center');
      else
        if fpltremv
          plot(ax,[pkindep(ii,(i-1)*2+1) pkindep(ii,(i-1)*2+1)], ax.YLim, '-.','Color',color(i,:));
          text(ax,pkindep(ii,(i-1)*2+1),yran(2)-0.1*range(yran),num2str(ii),'FontSize',6,...
            'HorizontalAlignment','center');
        end
      end
    end

  else
    for ii = 1: size(pkindep,1)
      plot(ax,[pkindep(ii,(i-1)*2+1) pkindep(ii,(i-1)*2+1)], ax.YLim, '--','Color',color(i,:));
      text(ax,pkindep(ii,(i-1)*2+1),yran(2)-0.1*range(yran),num2str(ii),'FontSize',6,...
        'HorizontalAlignment','center');
    end
  end
  ax.Box='on'; grid(ax,'on');
end
xlabel(f.ax(3),'Samples');
ylabel(f.ax(3),'Amplitude');
