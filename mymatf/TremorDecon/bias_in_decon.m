% bias_in_decon.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We have sort of seen that the complete procedure of deconvolution would 
% lead to a characteristic temporal resolution, ~0.25 s. For example, no 
% matter how saturated the synthetic sources are, the deconvolution would 
% probably fit many impulses around a dominant one, but the impulse 
% grouping, and secondary source removal likely discard those impulses
% anyway. Therefore, in the end, sources could be on average separated by
% each other by an amount that is comparable to the characteristic period
% of the seismogram, i.e., 0.25 s for the 4-Hz signal.
% -- This code aims to gain insight of this bias to the deconvolution.
% We first generate a bunch of random indice for the arrival times of ground
% truth sources, then put them into sequential time bins that are 0.25-s 
% wide. If one source index falls in one bin; all sources in the same time
% bin are regarded as the same source to mimic the resolution of 
% deconvolution. Then we just need to count the number of filled bins, 
% compared to total number of ground-truth sources.
%  
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/08/05
% Last modified date:   2024/08/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrun = 1000;  %num of realizations
for irun = 1: nrun
  seed = irun;

  tdura = 0.25;
  sps = 160;
  greenlen = pow2(9)*sps/40;
  bufsec = 1;
  msftaddm = bufsec*sps;  %buffer range for later CC alignment, +1 for safety
  rccmwsec = 0.5;
  rccmwlen = rccmwsec*sps;  %window length for computing RCC
  overshoot = rccmwlen/2; %RCC points to the center of computing window, so extra room is needed
  Twin=3*3600+(greenlen+msftaddm*2+overshoot*2-2)/sps; %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
  winlen=Twin*sps+1;
  satn=1/tdura*Twin;
  %variation of saturation level
  % nsat=[0.1 0.4 1 2 4 10 20 40 100];  % times of saturation, # of arrivals per duration 
  nsat=[0.4 1 2 4 10 20 40 100];  % times of saturation, # of arrivals per duration 
  nnsat = length(nsat);
  writes=round(nsat*satn); %how many templates to throw in, under different degrees of saturation
 
  %generate random arrival times
  rn = rand(writes(end),1);

  %set time bins
  tbincnts = 0:0.25:Twin;
  bincnts = tbincnts*sps+1;
  tbinedges = (tbincnts(2:end)+tbincnts(1:end-1))/2;
  tbinedges = [0 tbinedges];
  binedges = tbinedges*sps+1;
  nbins = length(binedges)-1;
  
  for insat = 1: nnsat
    rni = rn(1:writes(insat));
    rni = ceil(rni*winlen);
    ngt(irun,insat) = length(rni);
    N = histcounts(rni,binedges);
    Nfilled = sum(N>0);
    ndet(irun,insat) = Nfilled;
    %ratio of filled bins to all bins
    ratfilled(irun,insat) = Nfilled/nbins;
    %temporal density, related to above
    tden(irun,insat)= Nfilled/Twin;
    %ratio of detection, # of detection / # of ground-truth
    ratdec(irun,insat) = Nfilled/length(rni);
    % keyboard
  end

end

%%
nrow = 2; ncol = 3;
widin = 8.4;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
f = initfig(widin,htin,nrow,ncol);
pltxran = [0.08 0.98]; pltyran = [0.08 0.98]; % optimal axis location
pltxsep = 0.06; pltysep = 0.06;
optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

%num of GT srcs
ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');ax.XScale='log';
for irun = 1: nrun
  plot(ax,nsat,ngt(irun,:),'-',...
    'linew',1,'color','k','marker','o','markersize',3,...
    'markerfacec','k');
end
xticks(ax,nsat);
xlabel(ax,'Saturation');
ylabel(ax,'# of ground-truth sources');

%num of detections
ax=f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');ax.XScale='log';
for irun = 1: nrun
  plot(ax,nsat,ndet(irun,:),'-',...
    'linew',1,'color','k','marker','o','markersize',3,...
    'markerfacec','k');
end
xticks(ax,nsat);
xlabel(ax,'Saturation');
ylabel(ax,'# of detections');

%ratio of filled bins to all bins
ax=f.ax(3); hold(ax,'on'); ax.Box='on'; grid(ax,'on');ax.XScale='log';
for irun = 1: nrun
  plot(ax,nsat,ratfilled(irun,:),'-',...
    'linew',1,'color','k','marker','o','markersize',3,...
    'markerfacec','k');
end
xticks(ax,nsat);
xlabel(ax,'Saturation');
ylabel(ax,'ratio of filled bins to all bins');

%temporal density
ax=f.ax(4); hold(ax,'on'); ax.Box='on'; grid(ax,'on');ax.XScale='log';
for irun = 1: nrun
  plot(ax,nsat,tden(irun,:),'-',...
    'linew',1,'color','k','marker','o','markersize',3,...
    'markerfacec','k');
end
xticks(ax,nsat);
xlabel(ax,'Saturation');
ylabel(ax,'# of detections per sec');

%ratio of detection, # of detection / # of ground-truth
ax=f.ax(5); hold(ax,'on'); ax.Box='on'; grid(ax,'on');ax.XScale='log';
for irun = 1: nrun
  plot(ax,nsat,ratdec(irun,:),'-',...
    'linew',1,'color','k','marker','o','markersize',3,...
    'markerfacec','k');
end
xticks(ax,nsat);
xlabel(ax,'Saturation');
ylabel(ax,'# of detections / # of ground-truth');

delete(f.ax(6));

fname = strcat('aa.pdf');
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));

