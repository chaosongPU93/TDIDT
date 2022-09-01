function [f] = plt_srcrccwithoffset(impindepst,rccpaircat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_srcrccwithoffset(impindepst,rccpaircat)
%
% plot the corresponding running CC at the arrival time of the deconvolved
% impulse at each station, colorcoded by the time offset 12, 13, and 23 
% respectively. If there is some structure in source location on the fault,
% you might see high rcc value along with a nice progression in color 
% (which represents the offset), given time elapses. 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/07/01
% Last modified date:   2022/07/01 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

f.fig = figure;
f.fig.Renderer = 'painters';

nrow = 3;
ncol = 4;
% nsta = 3;
for ii = 1: nrow
  for jj = 1: ncol-1
    subplot(nrow,ncol,(ii-1)*ncol+jj);
    hold on; box on; grid on
    scatter(impindepst(:,(jj-1)*2+1),rccpaircat(impindepst(:,(jj-1)*2+1),jj),25,...
      impindepst(:,6+ii),'filled','MarkerEdgeColor',[.5 .5 .5]);
    colormap jet
    c=colorbar;
    ylim([-1 1]);
    box on;
    if jj == 1
      title('Sta 1-2');
      if ii == 1
        c.Label.String = sprintf('Offset 12 (samples)');
        nolabels(gca,3);
      elseif ii == 2
        c.Label.String = sprintf('Offset 13 (samples)');
        nolabels(gca,3);
      else
        xlabel('Samples');
        ylabel('Running CC');
        c.Label.String = sprintf('Offset 23 (samples)');
      end
    elseif jj == 2
      title('Sta 1-3');
      nolabels(gca,3);
    else
      title('Sta 2-3');
      nolabels(gca,3);
    end
    longticks(gca);
  end
  
  subplot(nrow,ncol,ii*ncol);
  hold on; box on; grid on
  scatter(impindepst(:,1),mean(rccpaircat(impindepst(:,1),:), 2),25,...
    impindepst(:,6+ii),'filled','MarkerEdgeColor',[.5 .5 .5]);
  prc2 = prctile(mean(rccpaircat(impindepst(:,1),:), 2), 2);
  med = median(mean(rccpaircat(impindepst(:,1),:), 2));
  ax = gca;
  plot(ax.XLim,[prc2 prc2],'--','Color',[.5 .5 .5],'linew',1.5);
  plot(ax.XLim,[med med],'k-','linew',1.5);
  colormap jet
  colorbar;
%   xlabel('Samples');
  ylim([-1 1]);
  box on;
  title('Mean');  
  longticks(gca);
  nolabels(gca,3);
end





