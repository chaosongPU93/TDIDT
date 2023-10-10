% function [off1ic,ccali,iloopoff,loopoff,mampsqsum]=locateMBlfe(famnum)
% locateMBlfe.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --This script aims to detect and locate Michael Bostock's LFEs 
% independently. I have the locations of each of his detected LFE fams, 
% his result is not guaranteed to be reliable in our context. So the idea
% is to read in a specific fam number, find all member LFEs from his LFE
% catalog. Choose a time window with a proper length around the timing, then
% do a 3-station detection similar to 4-s tremor detection. I think you can
% still enforce a threshold for the cricuit of time offsets 12-13+23 and mean
% CC between the same station pairs, in order to know which of his LFEs are
% detectable using our method. But the first step is to get the distribution
% of offset circuit, CC, and the resulting off12, off13 for each LFE. 
% --From off12 and off13, you can invert for locations for memeber LFEs, from
% which you would know the rough location for the entire fam given some 
% location uncertainty and noise (the fam location should be some average).
% But inverting from offset to map locations is not necessary, if your purpose
% is not the absolute location.
% --This might be very important to know whether some families are actually
% close in space. If some families are pretty close, then MB's catalog should
% have plenty redundant detection from multiple fams, then you probably need
% to lump them all, identify unique ones and discard others. If they are distant
% otherwise, then in particular cases like estimating the geodetic moment of 
% a certain region, you may include specific families only.
%
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/06/03
% Last modified date:   2023/06/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
% close all

defval('famnum',2);

format short e   % Set the format to 5-digit floating point

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
temppath = strcat(workpath, '/templates/');
%%% choose the data path here
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/PGCtrio');

FLAG = 'PGC';

%%%load MB's catalog
%load the LFE catalog of Michael Bostock, inside and outside the rectangle in 'locinterp002_4s.m'
%his time is approximately corrected to the postive waveform peak 
bosdir = ('/home/data2/chaosong/matlab/allan/BOSTOCK');
lfefnm = ('newlfeloc');
lfeloc = load(fullfile(bosdir, lfefnm));
loc0 = lfeloc(lfeloc(:,1)==2,:); %obtain the location of fam 002, lon0 and lat0
bostcat = ReformBostock(loc0(3),loc0(2),1);
% keyboard
bostcat = bostcat((bostcat(:,2)==2003 & bostcat(:,3)==3) | ...
                  (bostcat(:,2)==2004 & bostcat(:,3)==7) | ...
                  (bostcat(:,2)==2005 & bostcat(:,3)==9), :);
%format: [fam yyyy mm dd sec dx dy lon lat dep magnitude number-of-stations], 12 cols
%time should point to the peak, zero-crossing?

%choose a particular fam
fam = num2zeropadstr(famnum, 3)
bost = bostcat(bostcat(:,1)==famnum,:);

dates = unique(bost(:,2)*1e4 + bost(:,3)*1e2 + bost(:,4));
nd = length(dates);
  
%Use new LFE catalog
CATA = 'new';

% get permanent and polaris station rotation parameters, based on 40-sps data
sft2=0;     % centroid shift of station 2
sft3=0;     % centroid shift of station 3
[PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,sft2,sft3);
if ~isequal(CATA, 'fixed')
% if (famnum==2 && ~isequal(CATA, 'fixed')) || (famnum~=2 && ~isequal(CATA, 'fixed'))
  reftime = PERMROTS(1,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
  PERMROTS(:,4) = PERMROTS(:,4)-reftime;    % to make sure that 1st station is 0
  POLROTS(:,4) = POLROTS(:,4)-reftime;  
end
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;     % convert 2-3 columns to rad
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
  'LZB'];
POLSTA=['SSIB '           % polaris station names
  'SILB '
  'KLNB '
  'MGCB '
  'TWKB '];

stas=['PGC  '
  'SSIB '
  'SILB '];     % determine the trio and order, here the 1st sta is PGC
nsta=size(stas,1);         %  number of stations
 
%Sampling rate of the data will be used
sps=160;     % samples per second
  
hi=6.5;    % frequency band
lo=1.25;
npo=2;     % poles, passes of filters
npa=2;

% if famnum == 2
  msft = 16*sps/40;  %+1 for safety
  ccmin = 0.4;  % depending on the length of trace, cc could be very low
  loffmax = 2*sps/40;
% elseif famnum == 246
%   msft = 30*sps/40;  %+1 for safety
%   ccmin = 0.3;  % depending on the length of trace, cc could be very low
%   loffmax = 4*sps/40;
% 
% end

concentration=0.5; %in seconds; how concentrated is the coherent energy within the window?
cncntr=concentration*sps;   % in samples, 20

pltflag = 0;

%% identify and locate each MB's LFEs 
nlfe = 0;
k = 0;
bostsave = [];  %MB LFEs that actually been relocated
datesave = [];  %dates actually been analyzed

%detection window around MB's LFE timing
wlensec = 2;
wlen = wlensec*sps;
zccor = 0.2;

%cosine taper of the detection win
fractap = 0.25/wlensec;

%loop over days
for id = 1: nd
  fprintf('%d / %d \n', id, nd);
  date = dates(id);
  yr = floor(date/1e4);
  mo = floor((date-yr*1e4)/1e2);
  dy = date-yr*1e4-mo*1e2;

  jday = dat2jul(mo,dy,yr);
  JDAY = num2zeropadstr(jday,3);
  YEAR=int2str(yr);
  MO=day2month(jday,yr);     % EXTERNAL function, day2month, get the month of one particular date

  direc=[datapath, '/arch', YEAR,'/',MO,'/'];     % directory name
  prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];    %  path plus prefix of data file,
  %     disp(prename);

  [STAopt,~,~,fileflag] = rd_daily_bpdata(yr,jday,prename,stas,PERMSTA,POLSTA,...
  PERMROTS,POLROTS,sps,lo,hi,npo,npa,[],[],[],[]);
  STAtime = STAopt(:,1); %time column
  STAopt = STAopt(:,2:end); 

  if fileflag == 0    % means there are missing files
    fprintf('Day %s / %s will be omitted because of missing files. \n', yr, JDAY);
    continue    % continue to the next day
  else
    datesave = [datesave; date];  
  end
  
  %Bostock's LFE catalog on the same date
  bosti = bost(bost(:,2)==yr & bost(:,3)==mo & bost(:,4)==dy,:);
  bosti = sortrows(bosti, 5);
  nbosti = size(bosti,1);  
  nlfe = nlfe+nbosti;
  bostsave = [bostsave; bosti];

  %loop over LFEs of the same day
  for i = 1: nbosti 
    
    k = k+1;
    %chop a segment around MB's detection for computing best alignment
    tbosti = bosti(i,5);
    indst = round((tbosti-wlensec/2-zccor)*sps);
    inded = indst+wlen-1;
    optcc = detrend(STAopt(indst: inded, :));
    taper = tukeywin(wlen,fractap);
    for ista = 1: nsta
      optcc(:,ista) = taper.* optcc(:,ista);
    end
    ccmid = ceil(size(optcc,1)/2);
    ccwlen = round(size(optcc,1)-2*(msft+1));
    iup = 1;    % times of upsampling
    [off12con,off13con,ccali(k),iloopoff(k),loopoff(k)] = constrained_cc_interp(optcc(:,1:3)',ccmid,...
      ccwlen,msft,loffmax,ccmin,iup);
    % if a better alignment cannot be achieved, use 0,0
    if off12con == msft+1 && off13con == msft+1
      % off12con = 0;
      % off13con = 0;
      fprintf('MB LFE %d on the %dth day cannot be properly aligned \n',k,id);
    end
    off1ic(k,1) = 0;
    off1ic(k,2) = round(off12con);
    off1ic(k,3) = round(off13con);

    sigsta = [];
    for ista = 1: nsta
      sigsta(:,ista) = STAopt(indst-off1ic(k,ista):inded-off1ic(k,ista), ista);
    end
    rccmwsec = 0.5; 
    rccmwlen = rccmwsec*sps;  %window length for computing RCC
    %compute running CC between 3 stations
    [ircc,rcc12,rcc13,rcc23] = RunningCC3sta(sigsta,rccmwlen);
    rccpair = [rcc12 rcc13 rcc23];
    rcc = mean(rccpair(:,[1 2]), 2);

    STA12tr=sigsta(:,1).*sigsta(:,2);
    STA13tr=sigsta(:,1).*sigsta(:,3);
    STA32tr=sigsta(:,2).*sigsta(:,3);
    cumsumtr=cumsum(STA12tr)+cumsum(STA13tr)+cumsum(STA32tr); % sum of the cumsum of all traces
    %%% first get the squared sum of each 0.5s window, then get the maximum and start indice
    %%% this gives the proxy of Energy of 0.5s window, i.e. the integral of \dot(E)
    %%% of eqn 1 in Rubin et al. 2013
    %%% idiff<==>win [idiff+1, idiff+20] in regional ind<==>win [istart+idiff, istart+idiff+19] in global index  
    [mampsqsum(k),idiff]=max(cumsumtr(cncntr+1:end)-cumsumtr(1:end-cncntr));
    
    if pltflag
      %plot the overview of both the signal and noise together
      widin = 8;  % maximum width allowed is 8.5 inches
      htin = 4;   % maximum height allowed is 11 inches
      nrow = 1;
      ncol = 1;
      f = initfig(widin,htin,nrow,ncol); %initialize fig
      
      xran = [0.08 0.92]; yran = [0.15 0.90];
      xsep = 0.02; ysep = 0.04;
      optaxpos(f,nrow,ncol,xran,yran,xsep,ysep);
      
      ym = max(abs(sigsta(:)));
      yran = [-1.5*ym 1.5*ym];
      
      ax = f.ax(1); hold(ax,'on');
      yyaxis(ax,'left');
      time = STAtime(indst:inded);
      tstbuf = STAtime(indst);
      plot(ax,time-tstbuf,sigsta(:,1),'r-','linew',1);
      plot(ax,time-tstbuf,sigsta(:,2),'b-','linew',1);
      plot(ax,time-tstbuf,sigsta(:,3),'k-','linew',1); %,'linew',0.5
      xran = [time(1)-tstbuf time(end)-tstbuf]; % x axis range in sec
      xlim(ax,xran);
      ylim(ax,yran);
      longticks(ax,5);
      %plot the bostock's LFE catalog
      yloc = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))*0.06;
      scatter(ax,tbosti-tstbuf, yloc*ones(size(tbosti)),30,'g','filled');
      plot(ax,[tbosti-tstbuf tbosti-tstbuf],ax.YLim,'k--');
      text(ax,0.01,0.2,num2str(off1ic(k,2)),'unit','normalized',...
        'HorizontalAlignment','left','fontsize',10,'color','b'); % offset 12
      text(ax,0.05,0.2,num2str(off1ic(k,3)),'unit','normalized',...
        'HorizontalAlignment','left','fontsize',10,'color','k'); % offset 13
      text(ax,0.09,0.2,sprintf('%.3f',ccali(k)),'unit','normalized',...
        'HorizontalAlignment','left','fontsize',10,'color','k'); % mean CC of 3 pairs
      xlabel(sprintf('Time (s) on %d/%d/%d',dy,mo,yr),'FontSize',10);
      ylabel(ax,'Amplitude','FontSize',10);
      yyaxis(ax,'right');
%       plot(ax,ircc/sps,yran(2)/4*rcc+yran(2)/4*3,'-','Color',[.6 .6 .6],'linew',0.5); % scale with data amp.
      plot(ax,ircc/sps,yran(2)*rcc,'-','Color',[.6 .6 .6],'linew',0.5); % scale with data amp.
      ylim(ax,[-1 1]);
      ylabel(ax,'log_{10}{RCC}','FontSize',12);
      hold(ax,'off');
      
    end
        
  end %loop end for LFEs of the same day

  % keyboard
end %loop end for days

%% convert offsets from different fams to a uniform frame
if famnum ~= 2
  reffam = num2zeropadstr(2, 3);  % use '002' as the reference center
  [refperm, refpol] = GetRotsCommon(FLAG,reffam,CATA,datapath,0,0);
  reftref = refperm(1,4);     % reftime is the reference time at the 1st station, PGC here
  refperm(:,4) = refperm(:,4)-reftref;    % to make sure that 1st station is 0
  refpol(:,4) = refpol(:,4)-reftref;

  % conversion time constants, NOTE that the slow/fast time offset and offset between stas
  % is all calculated based on 40 sps
  conv12 = POLROTS(1,4) - refpol(1,4);   % station 2, SSIB
  conv13 = POLROTS(2,4) - refpol(2,4);    % station 3, SILB

  spsrat = sps/40;
  off1ic(:,2) = off1ic(:,2)-conv12*spsrat;  % get the new off based on the ref fam frame
  off1ic(:,3) = off1ic(:,3)-conv13*spsrat;
end
            
%% statistics for eligible detections
%%%plot the loopoff vs. CC of all MB LFEs
f = initfig(6,5,1,1); 
ax = f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); 
[mampsqsumst,indsort] = sort(mampsqsum,'ascend');
scatter(ax,loopoff(indsort),ccali(indsort),5,mampsqsumst,'o','filled','MarkerEdgeColor','none');
% scatter(ax,loopoff(indsort),ccali(indsort),5,mampsqsumst,'o','filled','MarkerEdgeColor','none');
oldc = colormap(ax, 'kelicol');
newc = flipud(oldc);
colormap(ax, newc);
caxis(ax, [prctile(mampsqsum,5) prctile(mampsqsum,95)]);
% text(ax,0.98,0.95,sprintf('Satur=%.1f',nsat(insat)),'Units','normalized',...
%     'HorizontalAlignment','right');
xran = [min(-5,-2.5*loffmax) max(5,2.5*loffmax)];
yran = [0.2 1];
axis(ax,[xran yran]);
yticks(ax,yran(1): 0.2 : yran(2));
plot(ax,ax.XLim,[ccmin ccmin],'k--');
plot(ax,[loffmax loffmax],ax.YLim,'k--');
plot(ax,[-loffmax -loffmax],ax.YLim,'k--');
c=colorbar(ax);
c.Label.String = 'sum of amp^2 of strongest 0.5s';
ylabel(ax,'Average CC');
xlabel(ax,sprintf('Summed time offsets of sta pairs (%d Hz)',sps));

%%
%%%plot the cumulative density of eligible MB LFEs
npass = sum(~(off1ic(:,2) == msft+1 & off1ic(:,3) == msft+1));
ratpass = npass/nlfe
ngood = sum(~(abs(off1ic(:,2)) <= 30 & abs(off1ic(:,3)) <= 30));
ratgood = ngood/nlfe

% den1d = density_pixel(off1ic(ind,2),off1ic(ind,3));
den1d = density_pixel(off1ic(:,2),off1ic(:,3));
% ftrans = 'interpchao';
% [imploc, ~] = off2space002(den1d(:,1:2),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
% den1dloc = [imploc(:,1:2) den1d(:,3)];
scale = 'linear';
% scale = 'log';

f = initfig(6,6,1,1); 
ax = f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on'); 
dum = den1d;
dum(dum(:,3)>1, :) = [];
if strcmp(scale,'log')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),20,dum(:,3),'o','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dum = sortrows(den1d,3);
dum(dum(:,3)==1, :) = [];
nmulti = sum(dum(1:end-1,3));
ratmulti = nmulti/nlfe
if strcmp(scale,'log')
  dum(:,3) = log10(dum(:,3));
end
scatter(ax,dum(:,1),dum(:,2),20,dum(:,3),'o','filled','MarkerEdgeColor',[.5 .5 .5]);
colormap(ax,newc);
c=colorbar(ax);
caxis(ax, [0 prctile(dum(:,3),95)]);
if strcmp(scale,'log')
  c.Label.String = strcat('log_{10}(# per pixel)');
elseif strcmp(scale,'linear')
  c.Label.String = strcat('# per pixel');  
end
axis(ax,'equal');
xran = [-msft/2 msft/2];
yran = [-msft/2 msft/2];
axis(ax,[xran yran]);
xlabel(ax,sprintf('PGC-SSIB offset (%d Hz)',sps));
ylabel(ax,sprintf('PGC-SILB offset (%d Hz)',sps));
