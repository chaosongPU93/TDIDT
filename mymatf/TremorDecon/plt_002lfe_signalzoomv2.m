function f=plt_002lfe_signalzoomv2(greenf,sigsta,rcc,impindepst,sps,xzoom,...
    tstbuf,dy,mo,yr,tbostplt,idxbst,yran2l,saveflag)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % This is the function to plot the seismograms of the stacked LFE templates
  % 'greenf' that are bandpassed. And a zoom-in window of the signal seismograms
  % 'signal', this is supposed to be a figure shown in the paper. 
  %
  %
  %
  % Chao Song, chaosong@princeton.edu
  % First created date:   2024/02/08
  % Last modified date:   2024/02/08
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  defval('saveflag',1);  %default is short-win detections
  
  lsig = size(sigsta,1); 
  nsta = 3; 
  lwlet = size(greenf,1);
  % lwlet = 12*sps;
  wletzoom = [2 8];
  lwletzoom = range(wletzoom)*sps;
  lsigzoom = range(xzoom)*sps;
  reflen = 25*sps;  %any length smaller is shrinked
  
  %%%compose a summary figure for each data win and save it, so that we will have a feeling for
  %%%all migrations
  widin = 8.4;  % maximum width allowed is 8.5 inches
  htin = 3;   % maximum height allowed is 11 inches
  nrow = 2;
  ncol = 1;
  f=initfig(widin,htin,nrow,ncol);
  
  % figxran = [0.08 0.94]; figyran = [0.1 0.95];
  % figxsep = 0.05; figysep = 0.04;
  % axpos=optaxpos(f,nrow,ncol,figxran,figyran,figxsep,figysep);
  axpos = [0.06 0.62 0.87*lwletzoom/reflen 0.35;
           % 0.06+0.87*lwletzoom/reflen+0.08 0.62 0.87*lsigzoom/reflen 0.35;];
           0.06 0.13 0.87*lsigzoom/reflen 0.35];
  for i = 1:nrow*ncol
    set(f.ax(i), 'position', axpos(i,:));
  end
  
  color = ['r';'b';'k'];
  
  if isempty(yran2l)
    ym = max(abs(sigsta(:)));
    yran=1.2*[-ym ym];
  else
    yran=yran2l;
  end
    
  %%%LFE templates, bandpassed
  ax=f.ax(1); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on'); 
  yyaxis(ax,'left');
  for i = 1: nsta
    p(i)=plot(ax,(1:lwlet)/sps, greenf(1:lwlet,i), '-','Color',color(i,:),'linew',1); 
  end
  xlim(ax,wletzoom);  
  % ylim(ax,yran/2);
  ylim(ax,[-0.5 0.5]);
  text(ax,0.98,0.75,'1.8-18 Hz','Units','normalized','HorizontalAlignment','right',...
    'FontSize',9);
  text(ax,0.98,0.9,'LFE templates at 002','Units','normalized',...
    'HorizontalAlignment','right','FontSize',9);
  text(ax,0.03,0.85,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1);
  xlabel(ax,'Time (s)','fontsize',10);
  ylabel(ax,'Amplitude','fontsize',10);
  % shrink(ax,range(xzoom)*sps/lwlet,1);
  % ax.Position(1)=axpos(1,1);
  xticks(ax,0: 2: lwletzoom);
  longticks(ax,2); 
%   nolabels(ax,2);
  ax.YColor=ax.XColor;

  yyaxis(ax,'right');
  mwlen = sps/2;
  [irccgf,rccgf12] = RunningCC(greenf(:,1), greenf(:,2), mwlen);
  [~,rccgf13] = RunningCC(greenf(:,1), greenf(:,3), mwlen);
  [~,rccgf23] = RunningCC(greenf(:,2), greenf(:,3), mwlen);
  % irccgf = irccgf+*sps;
  % rccgf = (rccgf12+rccgf13+rccgf23)/3;
  rccgf = (rccgf12+rccgf13)/2;
  % plot(ax,irccgf/sps,yran(2)*rccgf,'o','markersize',1,'Color',[.5 .5 .5]); % scale with data amp.
  plot(ax,irccgf/sps,rccgf,'o','color',[.8 .8 .8],'markersize',1);
  % legend(ax,p,{'PGC','SSIB','SILB'},'Location','southeast','NumColumns',3,...
  %   'fontsize',8);
  ylabel(ax,'Running CC','FontSize',10);
  ylim(ax,[-1 1]); 
  ax.YColor=[.6 .6 .6];
  hold(ax,'off');
     
  %%%seismograms of signal, zoom-in
  ax=f.ax(2); hold(ax,'on'); ax.Box = 'on'; grid(ax,'on'); 
  if isempty(rcc)
    %%%obtain a single best alignment based on the zoom-in segment
    msftadd = 20; 
    optcc = sigsta(xzoom(1)*sps+1: xzoom(2)*sps, :);
    ccmid = ceil(size(optcc,1)/2);
    ccwlen = round(size(optcc,1)-2*(msftadd+1));
    loffmax = 4*sps/40;
    ccmin = 0.01;  % depending on the length of trace, cc could be very low
    iup = 1;    % times of upsampling
    [off12con,off13con,ccali,iloopoff,loopoff] = constrained_cc_interp(optcc',ccmid,...
      ccwlen,msftadd,loffmax,ccmin,iup);
    % if a better alignment cannot be achieved, use 0,0
    if off12con == msftadd+1 && off13con == msftadd+1
      off12con = 0;
      off13con = 0;
      fprintf('This segment cannot be properly aligned, double-check needed \n');
    end
    sigseg(:, 1) = sigsta(xzoom(1)*sps+1: xzoom(2)*sps, 1);
    sigseg(:, 2) = sigsta(xzoom(1)*sps+1-off12con: xzoom(2)*sps-off12con, 2); % sta 2
    sigseg(:, 3) = sigsta(xzoom(1)*sps+1-off13con: xzoom(2)*sps-off13con, 3); % sta 3
  
    mwlen = sps/2;
    [ircc,rcc12] = RunningCC(sigseg(:,1), sigseg(:,2), mwlen);
    [~,rcc13] = RunningCC(sigseg(:,1), sigseg(:,3), mwlen);
    [~,rcc23] = RunningCC(sigseg(:,2), sigseg(:,3), mwlen);
    % ircc = ircc+*sps;
    % rcc = (rcc12+rcc13+rcc23)/3;
    rcc = (rcc12+rcc13)/2;
  else
    sigseg = sigsta(xzoom(1)*sps+1: xzoom(2)*sps, 1:3);
    rcc=rcc(xzoom(1)*sps+1: xzoom(2)*sps);
    ircc=1:length(rcc);
  end
  
  lseg = size(sigseg,1);
  yyaxis(ax,'left');
  for i = 1: nsta
    p(i)=plot(ax,(1:lseg)/sps, sigseg(:,i), '-','Color',color(i,:),'linew',1); 
  end
  xlim(ax,[0 lseg/sps]);  
  xticks(ax,0: 2: lseg/sps);
  % tkval = xzoom(1): 2: xzoom(2);
  % tklbl = cell(length(tkval),1);
  % for jj = 1: length(tkval)
  %   tklbl{jj} = num2str(tkval(jj));
  % end
  % xticklabels(ax,tklbl);
  ylim(ax,yran); 
  text(ax,0.99,0.9,'1.8-6.3 Hz, tremor velocity signal','Units','normalized',...
    'HorizontalAlignment','right','FontSize',10);
%   text(ax,0.01,0.2,'1.8-6.3 Hz','Units','normalized','HorizontalAlignment','left',...
%     'FontSize',10);
%   text(ax,0.01,0.85,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1);
  ylabel(ax,'Amplitude','FontSize',10);
  xlabel(ax,sprintf('Time (s) since %.1f s on %s %s %s',tstbuf+xzoom(1),dy,mo,yr),'FontSize',10);
  % xlabel(ax,sprintf('Time (s) since %.1f s on %s %s %s',tstbuf,dy,mo,yr),...
  %   'FontSize',10);
  longticks(ax,5); 
  % %simply plot the deconvolved sources
  % yloc = (yran(1)+range(yran)*0.05);
  % ind = find(impindepst(:,1)/sps >= xzoom(1) & impindepst(:,1)/sps <= xzoom(2));
  % scatter(ax,impindepst(ind,1)/sps-xzoom(1),yloc*ones(size(impindepst(ind,1))),10,'k');
  %plot deconvolved sources separated by clusters
  [clusibst,nclusibst,tclusibst]=clusters_in_burst(idxbst);
  ind = find(tclusibst(:,1)>=xzoom(1) & tclusibst(:,2)<=xzoom(2));
  clusibstxzoom = clusibst(ind);
  tclusibstxzoom = tclusibst(ind);
  nclus = length(ind);
  for i=1:2:nclus
    yloc = (yran(1)+range(yran)*0.05);
    impclus=clusibstxzoom{i};
    scatter(ax,impclus(:,1)/sps-xzoom(1),yloc*ones(size(impclus(:,1))),10,'^','k');
  end
  for i=2:2:nclus
    yloc = (yran(1)+range(yran)*0.05);
    impclus=clusibstxzoom{i};
    scatter(ax,impclus(:,1)/sps-xzoom(1),yloc*ones(size(impclus(:,1))),10,'v','k');
  end
  impuniclus = unique(cat(1,clusibstxzoom{:}),'rows');
  ind = find(impindepst(:,1)/sps >= xzoom(1) & impindepst(:,1)/sps <= xzoom(2));
  impiso = setdiff(impindepst(ind,:),impuniclus,'rows');
  scatter(ax,impiso(:,1)/sps-xzoom(1),yloc*ones(size(impiso(:,1))),10,'k');
  
  %%%7 events between 27.1 and 28.68 s form a cluster, highlight them
  % indclus = find(impindepst(:,1)/sps >= 27.1 & impindepst(:,1)/sps <= 28.68);
  % scatter(ax,impindepst(indclus,1)/sps-xzoom(1),yloc*ones(size(impindepst(indclus,1))),10,'k','filled');
  % impclus=clusibstxzoom{nclus};
  % scatter(ax,impclus(:,1)/sps-xzoom(1),yloc*ones(size(impclus(:,1))),10,'k','filled');
  %plot bostock's LFE detections
  yloc = (yran(1)+range(yran)*0.11);
  tbost=tbostplt-tstbuf;
  ind = find(tbost >= xzoom(1) & tbost <= xzoom(2));
  scatter(ax,tbost(ind)-xzoom(1),yloc*ones(length(ind),1),15,'g','filled'); %,'MarkerEdgeColor',[.5 .5 .5]
  
  oldc = colormap(ax,'kelicol');
  newc = flipud(oldc);
  colormap(ax,newc);
  caxis(ax,xzoom);
  ax.YColor=ax.XColor;
  
  yyaxis(ax,'right');
  % plot(ax,ircc/sps,yran(2)*rcc,'o','markersize',1,'Color',[.5 .5 .5]); % scale with data amp.
  plot(ax,ircc/sps,rcc,'o','color',[.8 .8 .8],'markersize',1);
  % lgd=legend(ax,p,{'PGC','SSIB','SILB'},'Location','north','NumColumns',3,...
  %   'fontsize',8);
  % set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));
  ylabel(ax,'Running CC','FontSize',10);
  ylim(ax,[-1 1]); 
  ax.YColor=[.6 .6 .6];
  hold(ax,'off');
  
%   keyboard
