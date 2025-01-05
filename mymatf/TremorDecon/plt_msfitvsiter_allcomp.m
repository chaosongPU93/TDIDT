function [f1,f2] = plt_msfitvsiter_allcomp(mfit,varsig,mfitort,varsigort,...
  mfitvert,varsigvert,stas,ista,detecttype,saveflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mfit = msfitvsiter(ampit,sigsta,greenf,zcrosses,ista,sps)
%
% According to the misfit computed from 'msfitvsiter.m' on all components
% (optimal, ort., and vertical), this function simply plots the misfit vs. the
% iteration in 2 ways, one is absolute value of variance of residual, the other
% is relative to the variance of data, ie, before iteration 0.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/06/21
% Last modified date:   2024/06/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nit = length(mfit);

%% absolute value of variance of residual
widin = 6;  % maximum width allowed is 8.5 inches
htin = 3.5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 2;
f1 = initfig(widin,htin,nrow,ncol);

pltxran = [0.08 0.96]; pltyran = [0.12 0.95]; % optimal axis location
pltxsep = 0.1; pltysep = 0.05;
optaxpos(f1,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
ax=f1.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
plot(ax,1:nit,mfit,'r-','linew',1);
ylabel(ax,'Variance of residual');
xlabel(ax,'Iteration #');
legend(ax,'Optimal','Location','east');
text(ax,0.99,0.95,strtrim(stas(ista,:)),...
  'Units','normalized','HorizontalAlignment','right','fontsize',10);
ax=f1.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
plot(ax,1:nit,mfitort,'b-','linew',1);
plot(ax,1:nit,mfitvert,'k-','linew',1);
legend(ax,'Orthogonal','Vertical','Location','best');
ylabel(ax,'Variance of residual');
xlabel(ax,'Iteration #');

if saveflag
  % orient(f.fig,'landscape');
  fname = strcat('resvsiter_',strtrim(stas(ista,:)),detecttype,'.pdf');
  print(f1.fig,'-dpdf',...
    strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
end

%% variance of residual relative to the variance of data
relavarort = mfitort/varsigort*100;
relavarvert = mfitvert/varsigvert*100;
relavar = mfit/varsig*100;

widin = 3.2;  % maximum width allowed is 8.5 inches
htin = 3.5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 1;
f2 = initfig(widin,htin,nrow,ncol);

pltxran = [0.18 0.95]; pltyran = [0.12 0.96]; % optimal axis location
pltxsep = 0.1; pltysep = 0.05;
optaxpos(f2,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);
ax=f2.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax, 'on');
plot(ax,1:nit,relavar,'r-','linew',1);
plot(ax,1:nit,relavarort,'b-','linew',1);
plot(ax,1:nit,relavarvert,'k-','linew',1);
ylabel(ax,'Relative variance of residual (%)');
xlabel(ax,'Iteration #');
legend(ax,'Optimal','Orthogonal','Vertical','Location','best');
text(ax,0.99,0.95,strtrim(stas(ista,:)),...
  'Units','normalized','HorizontalAlignment','right','fontsize',10);
% ylim(ax,[0 120]);

if saveflag
  % orient(f.fig,'landscape');
  fname = strcat('relaresvsiter_',strtrim(stas(ista,:)),detecttype,'.pdf');
  print(f2.fig,'-dpdf',...
    strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2024/figures/',fname));
end





