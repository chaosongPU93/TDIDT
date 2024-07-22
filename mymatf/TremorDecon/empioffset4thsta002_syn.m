% empioffset4thsta002_syn.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script aims to obtain an empirical estimate of the time offset needed
% for the 4th station among LZB, TWKB, MGCB, KLNB to align with the detected
% sources. The estimate comes from the 4th stational check result to 4-s 
% tremor detections in 'detection002_4s.m'. The estimate will be used to
% determine the offset needed for a deconvolved LFE source based on its
% location, which is the output in codes like 'deconvbursts002_ref_4s_exp.m'.
%
% --4-s tremor detections checked by diff 4th stations are diff. The same 
%   location have duplicates. Ideally, the observational (or theroratical)
%   predication of the time offset between sta 4 and 1 is smoothly varying
%   wrt the source location, which is also your eyes seeing. However, there
%   are many 'spikes' in the data sets on top of a smooth surface.
% --to smooth out these spikes, there is a testing script 
%   'testempioffset4thsta002' that tries a few ways, like 2D 
%   interpolation, gradient estimate, and 2D plane fitting. The test shows
%   that interpolation does not work at all; the gradient might not be very
%   stable; plane fitting seems like the best solution so far in terms of 
%   having a 1st-order estimate for the time offset. Therefore, this script
%   is going to use plane fitting only. 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/07/10
% Last modified date:   2024/07/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e   % Set the format to 5-digit floating point
clear
clc
close all

  %% for easy testing
  defval('normflag',0); %whether to normalize templates
  defval('rccmwsec',0.5); %moving win len in sec for computing RCC
  
  %% Initialization
  %%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
  %%% AND if using the same family, same station trio
  % if pltflag
  %     set(0,'DefaultFigureVisible','on');
  % else
  %     set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
  % end
  
  [scrsz, resol] = pixelperinch(1);
  
  % WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band
  % (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs
  
  workpath = getenv('ALLAN');
  datapath = strcat(workpath,'/data-no-resp');
  temppath = strcat(datapath, '/templates/PGCtrio/');
  rstpath = strcat(datapath, '/PGCtrio');
  
  stas=['PGC  '
    'SSIB '
    'SILB '
    %   'LZB  '
    %   'TWKB '
    %   'MGCB '
    'KLNB '
    ]; % determine the trio and order, here the 1st sta is PGC
  nsta=size(stas,1);         %  number of stations

  %%%specify distribution for source location
  distr='UN';  % uniform distribution
  
  %%%specify distribution for source location
  % distrloc = 'custompdf'; %using a custom PDF function
  distrloc = 'uniform'; %uniformly random in a specified region,
  
  fracelsew=0; %0.25 %0.5; %0.6; %The ratio of elsewhere events to local events.  Zero for Discrete Ide.
  
  %%%specify if considering the physical size of each source
  % physicalsize = 1;
  physicalsize = 0;
  
  %%%diameter of physical size
  if physicalsize
    diam=0.15;% 0.5; %0.6; %
  else
    diam=0;
  end
  
  %%%specify regime for transformation from time offset to map location
  % ftrans = 'interpArmb';
  % ftrans = 'interpArmbreloc';
  ftrans = 'interpchao';
  
  %%%specify if forcing a min speration of arrival time for 2 events from the same spot
  % forcesep = 1;
  forcesep = 0;
  
  %%%specify shape of the source region
  srcregion='ellipse';
  % srcregion='rectangle';
  % srcregion='circle';
  
  %%%flag for plot the src distribution for each region size
  pltsrc = 1;
  % pltsrc = 0;
  
  %%%flag for plot the projection along min-scatter for each region size & saturation level
  % pltproj = 1;
  pltproj = 0;
  
  %variation of source region size
  if strcmp(srcregion,'ellipse')
    semia = 1.75*(0.6:0.2:2.0);
    semib = 1.25*(0.6:0.2:2.0);
    nreg = length(semia);
  end

  tdura = 0.25;  %must be consistent with what synthetics actually used
% tdura = 0.4;

  %times of saturation
  % nsat=[0.1 0.4 1 2 4 10 20 40 100];
  nsat=[0.4 1 2 4 10 20 40 100];
  nnsat = length(nsat);
  
  %length of each simulation
  sps = 160;%Sampling rate of the data will be used, samples per second
  greenlen = pow2(9)*sps/40;
  bufsec = 1;
  msftaddm = bufsec*sps;  %buffer range for later CC alignment, +1 for safety
  rccmwsec = 0.5;
  rccmwlen = rccmwsec*sps;  %window length for computing RCC
  overshoot = rccmwlen/2; %RCC points to the center of computing window, so extra room is needed

  % Twin=0.5*3600+(greenlen+msftaddm*2+overshoot*2-2)/sps; %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
  % nrun = 6;

  Twin=3*3600+(greenlen+msftaddm*2+overshoot*2-2)/sps; %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
  nrun = 1;

  winlensec=4; 
  
  cyclskip = 0;
  %   mshift=24; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
  %   loopoffmax=5; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
  %   xcmaxAVEnmin=0.6; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?
  mshift=26*sps/40; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
  loopoffmax=2.1*sps/40; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
  xcmaxAVEnmin=0.44; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?
  
  ireg = 6;
  xaxis = semia(ireg);
  yaxis = semib(ireg);
  
  insat = 4;  
  
  stasnew=[
    'KLNB '
    ];  % twkb lzb mgcb
  nstanew=size(stasnew,1);

  tmrrst = load([workpath,'/synthetics/tmrall',num2str(winlensec),'.',distr,'.',int2str(sps),...
    'sps.',srcregion(1:3),'_',num2str(xaxis),'-',num2str(yaxis),'.diam',num2str(diam),...
    '.else',num2str(fracelsew,2),'nsat',num2str(nsat(insat)),stasnew(nstanew,:)],'w');
  
%% 2D plane fitting
data = tmrrst;

modparam = [];
modparammap = [];
off14pred = [];
for i = nstanew: nstanew
  colflg = 30+i;  % column num for the flag to indicate if the 4th sta passes the check, 43-46
  colloff = 31+i; % column num for the summed difference in offset needed to align 4th sta and original trio, 47-50
  colioff = 32+i; % column num for the offset needed to align 4th sta and original trio, 51-54
  colcc = 33+i; % column num for ave CC of aligning 4th sta and original trio, 55-58
  ccavethres = 0.8*xcmaxAVEnmin;  %% average cc threshold
  offavethres = 2.1*sps/40;    % average differential offset threshold
  
  ind = find(data(:,colflg)==1);
  loff = data(ind,colloff);  %convert to the same sps 
  ioff = data(ind,colioff);  %convert to the same sps, essentially 'off14'
  cc = data(ind,colcc);
  off12 = data(ind,2);  %note that the location off12 and off13 are already consistent here
  off13 = data(ind,3);
  loce = data(ind,1);
  locn = data(ind,2);
  mloff = median_pixel(off12,off13,loff);
  mioff = median_pixel(off12,off13,ioff);
  mcc = median_pixel(off12,off13,cc);

%   if isequaln(data, hfall)
%     if i==1
%       ran = [-30 20 -75 35];
%     elseif i==2
%       ran = [-30 25 -70 40];
%     elseif i==3
%       ran = [-30 25 -65 40];
%     elseif i==4
%       ran = [-40 45 -70 50];
%     end
%     ind2 = find(off12>=ran(1) & off12<=ran(2) & off13>=ran(3) & off13<=ran(4));
%     off12 = off12(ind2);
%     off13 = off13(ind2);
%     ioff = ioff(ind2);
%     mloff = mloff(mloff(:,1)>=ran(1) & mloff(:,1)<=ran(2) & mloff(:,2)>=ran(3) & mloff(:,2)<=ran(4), :);
%     mioff = mioff(mioff(:,1)>=ran(1) & mioff(:,1)<=ran(2) & mioff(:,2)>=ran(3) & mioff(:,2)<=ran(4), :);
%     mcc = mcc(mcc(:,1)>=ran(1) & mcc(:,1)<=ran(2) & mcc(:,2)>=ran(3) & mcc(:,2)<=ran(4), :);
%   end

  
  %%%plane fitting
%   [xData, yData, zData] = prepareSurfaceData(mioff(:,1),mioff(:,2),mioff(:,3));
  [xData, yData, zData] = prepareSurfaceData(off12/sps,off13/sps,ioff/sps);

  % Set up fittype and options.
  opts = fitoptions( 'Method', 'NonlinearLeastSquares' ); %
  opts.Display = 'Off';
  opts.Robust = 'Bisquare';
  opts.StartPoint = [0.8 0.8 0.1];
  
  % Fit model to data.
  %%%plane fit model
  ftplane = fittype( 'a*x+b*y+c', 'independent', {'x', 'y'}, 'dependent', 'z' );
  [fitresult, gof{i},output{i}] = fit( [xData, yData], zData, ftplane, opts );
  modparam(i,:) = [fitresult.a fitresult.b fitresult.c];
  off14pred = fitresult.a .* xData + fitresult.b .* yData + fitresult.c;
  
  %test if the returned goodness-of-fit metrics are also weighted, answer is YES 
  sseauto = gof{i}.sse;  %sum of squares due to error
  ssecalc = sum((zData-off14pred).^2); %assuming equal weights of 1
  rmseauto = gof{i}.rmse;  %root mean squared error
  rmsecalc = sqrt(sum((zData-off14pred).^2)./(length(zData)-2));  %assuming equal weights of 1
  
%   %%%parabolic fit model
%   ftpara = fittype( '(x-x0)^2/a^2+(y-y0)^2/b^2+z0', 'independent', {'x', 'y'}, 'dependent', 'z' );
%   opts.StartPoint = [0.5 0.5 0.5 0.5 0.5];
%   [fitresult, gof] = fit( [xData, yData], zData, ftpara, opts );
%   modparam(i,:) = [fitresult.a fitresult.b fitresult.x0 fitresult.y0 fitresult.z0];
%   off14pred{i} = (off12-(fitresult.x0)).^2./(fitresult.a)^2+(off13-(fitresult.y0)).^2./(fitresult.b)^2+...
%     (fitresult.z0);

  %%% Plot fit with data.
  widin = 8;  % maximum width allowed is 8.5 inches
  htin = 4.5;   % maximum height allowed is 11 inches
  nrow = 1;
  ncol = 2;
  f=initfig(widin,htin,nrow,ncol);

%   figxran = [0.1 0.95]; figyran = [0.12 0.95];
%   figxsep = 0.1; figysep = 0.07;
%   axloc=optaxpos(f,nrow,ncol,figxran,figyran,figxsep,figysep);

  ax = f.ax(2); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  p1 = scatter3(xData, yData, zData, 15,[.3 .3 .3],'filled','MarkerEdgeColor','w');
  %     aa=feval(fitresult,xData,yData);
  %     surf(ax,xData, yData,aa);
  p2 = plot(fitresult);
  colormap('jet');
  %     c=colorbar;
  %     c.Label.String = sprintf('off14 at %d sps',sps);
  legend([p1 p2], sprintf('\\Delta{t}'), 'Plane fit', 'Location', 'northwest');
  xlabel(sprintf('\\Delta{t}_{12} (s)'));
  ylabel(sprintf('\\Delta{t}_{13} (s)'));
  zlabel(sprintf('\\Delta{t}_{14} (s)'));
  axis('equal');
  zlim([-1.5 1.5]);
  xlim([-1 1]);
  ylim([-1 1]);
  text(0.02,0.10,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1);
  view(243,5.3);
  %     view(2);

  ax = f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  error = zData-off14pred;
  binw=1/sps;
  histogram(ax,error,'Normalization','count','BinWidth',binw,...
    'Facec','k','edgec','none');
  p1=plot(ax,[mean(error) mean(error)],ax.YLim,'r--','linew',2);
  p2=plot(ax,[rmseauto rmseauto],ax.YLim,'k--','linew',1.5);
  plot(ax,[-rmseauto -rmseauto],ax.YLim,'k--','linew',1.5);
  text(ax,0.98,0.95,sprintf('%s',strtrim(stasnew(i,:))),...
    'Units','normalized','HorizontalAlignment','right','fontsize',12);
  xran = [-0.25,0.25];
  %     [muHat,sigmaHat] = normfit(error);
  %     pdfhat = normpdf(xran(1):binw:xran(2),muHat,sigmaHat)*binw;
  %     p3=plot(ax,xran(1):binw:xran(2),pdfhat,'r','linew',2);
  xlabel(ax,sprintf('Deviation from plane fitting (s)'));
  ylabel(ax,'Count');
  xlim(ax,xran);
  xticks(ax,xran(1):0.25:xran(2));
  text(ax,0.02,0.10,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',1);
  legend([p1 p2], 'Mean', '\pm1 RMSE', 'Location', 'northwest');

end  
  
fname = strcat('planefit',strtrim(stasnew(i,:)),'_syn.pdf');
print(f.fig,'-dpdf',strcat('/home/chaosong/Pictures/',fname));
  
  
  
  
  
  
  