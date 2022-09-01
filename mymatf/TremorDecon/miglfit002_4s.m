% miglfit002_4s.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to analyse and linear fit several relatively large
% scale migrations imaged by detections around 002 region. These RTMs 
% typically migrated through the high-density ellipse around 002. The 
% purpose is to get the a sense of the migration direction, speed, and a
% the time needed to pass the ellipse given its size.
% 
% These information might be useful for determining the appropriate max.
% threshold for the separation between detections to group them into 
% temporally-cluster tremor bursts. 
% 
% --Through the plots, looks like the direction ranges from 115;130;130;145
%   degrees, while the speed through the 002 ellipse is around 20 km/h.
%   Let's take the direction as the median, 130. The time needed to pass
%   through the ellipse is about 0.1 h = 6 m = 360 s.
% --Some useful scale: If speed is 30 km/h, then time needed for the same
%   distance 2 km is about 4 m = 240 s; Given that the location error of 
%   +-1 sample at 40 sps using PGC trio possibily is around +-0.75 km in
%   in the strike direction that is also close to the propagation direction,
%   then it means, change in 1 sample would need 1.5 m = 90 s if speed is 
%   30 km/h and 135 s if speed is 20 km/h
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/03/29
% Last modified date:   2022/03/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
format short e   % Set the format to 5-digit floating point
clear
clc
close all

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

freqflag='hf';  % flag to indicate whether to do hf or lf;

FLAG = 'PGC'; % detector
  
fam = '002';   % family number

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

sps = 40;

iup = 4;  % upsample 4 times

%load detections
if isequal(fam,'002')
  cyclskip = 0;
  mshift=26+cyclskip; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
  loopoffmax=2.1; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
  xcmaxAVEnmin=0.44; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?
end
hi=6.5;    % frequency band
lo=1.25;
npo=2;     % poles, passes of filters
npa=2;
winlensec=4;
winoffsec=1;        % window offset in sec, which is the step of a moving window
winlen=winlensec*sps;      % length in smaples

PREFIX = strcat(fam,'.up.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),'.',...
                        int2str(npo),int2str(npa),'.ms', int2str(mshift));

%load all detections
hfall = load(strcat(rstpath, '/MAPS/eloc.pgc.tdec_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',...
            num2str(winlen/sps),'s',num2str(sps),'sps4add_'));
daycol = 14;
seccol = 16;  % 2 choices here, one is the center of detecting window, other is the start of strongest .5 s 
hfall = sortrows(hfall, [daycol, seccol]);
%format as follows:
%%% 8+34+4*nstanew cols, if 4 new stas, then will be 58 cols
%%% UPDATED at 2021/06/23
%%%   dx,dy,lon,lat,dep,tori,off12,off13 (integer samples at upsampled sps)
%%%   1:off12(n) 2:off13(n) 3:off12sec 4:off13sec 5:fam 6:date 7:timswin(n)
%%%   8:timswin(n)-winlensec/2+idiff/sps 9:xcmaxAVEnbang(nin) 10:loopoff(n)
%%%   11:cumsumtrdiff  12:cumsumtrdiff/cumsumtr(winlen)
%%%   for the rest, n=n+4;
%%%   9:Ampsq(nin) 10:Prior(nin) 11:Post(nin) 12:Prior4(nin) 13:Post4(nin) 14:Prior8(nin)
%%%   15:Post8(nin) 16:Prior16(nin) 17:Post16(nin) 18:cc(nin) 19:ccprior(nin)
%%%   20:ccpost(nin) 21:ccprior4(nin) 22:ccpost4(nin) 23:STAamp(nin,2) 24:STAamp(nin,3)
%%%   25:sumsSTA12n(n,iSTA12bang) 26:sumsSTA13n(n,iSTA13bang)
%%%   27:sumsSTA32n(n,iSTA32bang) 28:sigma(nin) 29:sigma12(nin) 30:sigma13(nin)
%%%   in(1:nstanew) loff(1:nstanew) ioff(1:nstanew) ccmaxave(1:nstanew)

%%
%Time ranges of RTMs
trange = [
  2005255   1.35*3600    1.70*3600;
  2005255   9.55*3600   10.10*3600;
  2005255  16.10*3600   16.55*3600;
  2003063   1.70*3600    2.05*3600;
  ];
                      
                      
xran = [-7 3];
yran = [-4 4];

angle = 0:5:355;

slopehf = zeros(size(trange,1),length(angle));
rmsehf = zeros(size(trange,1),length(angle));
angrmsehf = zeros(size(trange,1),1);
angslopehf = zeros(size(trange,1),1);                      
                      
% create fit object with free constraints
fttpfree = fittype( @(a,b,x) a*x+b);

x0 = 0.2;
y0 = 0.2;
semia = 1.75;
semib = 1.0;
angrot = 45;
[xell, yell] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);
        
for i = 1: size(trange,1)
  % for i = 1: 211
  %     i=46;
  disp(trange(i,:));
  indhf = find(hfall(:,daycol)==trange(i,1) & hfall(:,seccol)>=trange(i,2) & ...
    hfall(:,seccol)<=trange(i,3));
  mighf = hfall(indhf,:);
  
  %%% find best angle, now is mainly to get the variation of se and slope with the trial angle
  for iang = 1: length(angle)
    %%% propagation trial of hf
    mighfdum = mighf;
    for j = 1: size(mighf,1)
      x0 = mighfdum(j,1);
      y0 = mighfdum(j,2);
      [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),[0 0]);
      mighfdum(j,1) = newx;
      mighfdum(j,2) = newy;
    end
    % linear robust least square
    [fitobj,gof,~] = fit(mighfdum(:,seccol)/3600, mighfdum(:,1),fttpfree,'Robust','Bisquare',...
      'StartPoint',[1 1]);
    % output fit parameters
    coef = coeffvalues(fitobj);
    slopehf(i,iang) = coef(1);
    rmsehf(i,iang) = gof.rmse;
    
  end
                      
  %%% best angle estimate from hf
  ind = find(slopehf(i,:)>0);
  ind3 = find(rmsehf(i,ind)==min(rmsehf(i,ind)));     % rmse, Root Mean Squared Error, the smaller, the better
  if length(ind3) > 1
    disp(strcat({'multiple angles have same rmse: '},num2str(angle(ind3))));
  end
  angrmsehf(i) = angle(ind(ind3(1)));
  
  ind6 = find(slopehf(i,:)==max(slopehf(i,:))); % one with the largest slope, i.e., migrating speed
  if length(ind6) > 1
    disp(strcat({'multiple angles have same slope:'},num2str(angle(ind6))));
  end
  angslopehf(i) = angle(ind6(1));
                    
                    
  %%% define and position the figure frame and axes of each plot
  f.fig=figure;
  widin = 8;  % maximum width allowed is 8.5 inches
  htin = 9;   % maximum height allowed is 11 inches
  set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
  nrow = 2;
  ncol = 2;
  for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
  end
  
  %%% reposition
  set(f.ax(1), 'position', [ 0.08, 0.64, 0.36, 0.32]);
  set(f.ax(2), 'position', [ 0.52, 0.64, 0.36, 0.32]);
  set(f.ax(3), 'position', [ 0.08, 0.35, 0.36, 0.22]);
  set(f.ax(4), 'position', [ 0.52, 0.35, 0.36, 0.22]);
  
  % subplot 1 of figure i
  hold(f.ax(1),'on');
  plot(f.ax(1),[-100 100],[0 0],'k--');
  plot(f.ax(1),[0 0],[-100 100],'k--');
  plot(f.ax(1),xell,yell,'k-','linew',1.5);
  f.ax(1).FontSize = 9;
  mighf = sortrows(mighf,-seccol);
  scatter(f.ax(1),mighf(:,1),mighf(:,2), 20, mighf(:,seccol)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
  colormap(f.ax(1),'jet');
  c=colorbar(f.ax(1),'SouthOutside');
  pos = f.ax(1).Position;
  c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
  %     c.TickLabels=[];
  juldate = num2str(trange(i,1));
  yr = str2double(juldate(1:4));
  date = str2double(juldate(5:end));
  a = jul2dat(yr,date);
  mo = a(1);
  if mo == 9
    mo = {' Sep. '};
  elseif mo == 7
    mo = {' Jul. '};
  else
    mo = {' Mar. '};
  end
  day = num2str(a(2));
  yr = num2str(a(3));
  c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
  %     c.Label.String = strcat(num2str(trange(i,1)),' of HF',' (hr)');
  c.Label.FontSize = 11;
  caxis(f.ax(1),[trange(i,2)/3600 trange(i,3)/3600])
  text(f.ax(1),0.85,0.1,'HF','FontSize',12,'unit','normalized');
  text(f.ax(1),0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
  text(f.ax(1),0.04,0.1,'PGC','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
  f.ax(1).Box = 'on';
  grid(f.ax(1), 'on');
  axis(f.ax(1), 'equal');
  f.ax(1).GridLineStyle = '--';
  f.ax(1).XAxisLocation = 'top';
  medxhf = median(mighf(:,1));
  medyhf = median(mighf(:,2));
  [rotx, roty] = complex_rot(0,1,-angrmsehf(i));
  xvect = [medxhf-rotx medxhf+rotx];
  yvect = [medyhf-roty medyhf+roty];
  drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5);
  [rotx, roty] = complex_rot(0,1,-angslopehf(i));
  xvect = [medxhf-rotx medxhf+rotx];
  yvect = [medyhf-roty medyhf+roty];
  drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
  xticks(f.ax(1),xran(1):1:xran(2));
  yticks(f.ax(1),yran(1):1:yran(2));
  xlabel(f.ax(1),'E (km)','fontsize',11);
  ylabel(f.ax(1),'N (km)','fontsize',11);
  axis(f.ax(1),[xran yran]);
  text(f.ax(1),0.5,0.1,num2str(i),'FontSize',12,'unit','normalized');
  text(f.ax(1),0.84,0.95,strcat(num2str(size(mighf,1)),{' detections'}),'FontSize',8,...
    'unit','normalized','horizontalalignment','center');
  text(f.ax(1),0.84,0.90,strcat({'in '},num2str(trange(i,3)-trange(i,2)),{' s'}),'FontSize',8,...
    'unit','normalized','horizontalalignment','center');
  rate = sprintf('%.3f',size(mighf,1)/(trange(i,3)-trange(i,2)));
  text(f.ax(1),0.84,0.85,strcat({'rate: '},rate),'FontSize',8,...
    'unit','normalized','horizontalalignment','center');
  text(f.ax(1),0.84,0.78,strcat(num2str(angrmsehf(i)),{'{\circ}'}),'FontSize',10,...
    'unit','normalized','horizontalalignment','center');
  hold(f.ax(1),'off');
  
  % subplot 3 of figure i
  hold(f.ax(3),'on');
  %%% Actually used is the pre-determined best prop direc to do the fitting
  mighfdum = mighf;
  for j = 1: size(mighf,1)
    x0 = mighfdum(j,1);
    y0 = mighfdum(j,2);
    [newx,newy] = coordinate_rot(x0,y0,-(angrmsehf(i)-90),[0 0]);
    mighfdum(j,1) = newx;
    mighfdum(j,2) = newy;
  end
  f.ax(3).FontSize = 9;
  scatter(f.ax(3),mighfdum(:,seccol)/3600,mighfdum(:,1),20,[0.6 0.6 0.6],'filled','o',...
    'MarkerEdgeColor','k');
  
  % create fit object  
  [fitobjhfprop,gofhfprop,outphfprop] = fit(mighfdum(:,seccol)/3600, mighfdum(:,1),fttpfree,'Robust',...
    'Bisquare','StartPoint',[1 1]);
  % output fit parameters
  coefprop = coeffvalues(fitobjhfprop);
  slopeprophf(i) = coefprop(1);
  intcptprophf(i) = coefprop(2);
  fitprophf = feval(fitobjhfprop,mighfdum(:,seccol)/3600);
  plot(f.ax(3),mighfdum(:,seccol)/3600,fitprophf,'-','linewidth',2,'color',[0.6 0.6 0.6]);
  
  f.ax(3).Box = 'on';
  grid(f.ax(3), 'on');
  f.ax(3).GridLineStyle = '--';
  xran1 = [trange(i,2)/3600 trange(i,3)/3600];
  %     yran1 = [round(min([mighfdum(:,1);miglfdum(:,1)]))-1 ...
  %              round(max([mighfdum(:,1);miglfdum(:,1)]))+1];
  aa = round(prctile(mighfdum(:,1), 98));
  bb = round(prctile(mighfdum(:,1), 2));
  aap = aa + ceil((aa-bb)/6);
  bbp = bb - ceil((aa-bb)/3);
  yran1 = [bbp aap];
  xlim(f.ax(3),xran1);
  ylim(f.ax(3),yran1);
  text(f.ax(3),0.04,0.91,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
  xlabel(f.ax(3),'Time (hr)','fontsize',11);
  ylabel(f.ax(3),'Dist. along prop. (km)','fontsize',11);
  
  % compute the HF weights in robust linear regression, see NOTES above
  rhf = outphfprop.residuals;   % usual residuals
  x = [mighfdum(:,seccol)/3600];
  hatmat = x*inv(x'*x)*x';
  h = zeros(size(hatmat,1),1);    % leverage of least square
  for jj = 1 : size(hatmat,1)
    h(jj) = hatmat(jj,jj);
  end
  radj = rhf./sqrt(1-h);      % adjusted residuals
  K = 4.685;
  s = mad(rhf,1)/0.6745;
  u = radj/(K*s);
  wthf = zeros(length(u),1);    % rubust weight of next iteration
  for jj = 1 : length(u)
    if abs(u(jj)) < 1
      wthf(jj) = (1-(u(jj))^2)^2;
    else
      wthf(jj) = 0;
    end
  end
  
  % get the standard error of the estimated parameters, may indicate the compare the quality
  % of fitting cross different RTMs, NOTE that rmse, mse are definitely not a good measure
  % the obtained CI is comfirmed to be correct by comparing it with the 'fitobjhfprop'
  slopesehf = gofhfprop.rmse./sqrt(sum((x-mean(x)).^2));
  est = slopeprophf(i);
  slopeprophfCI(i,:) = confidence_interval_general(est,slopesehf,length(x)-2,95);
  interceptsehf = slopesehf.*sqrt(sum(x.^2)./length(x));
  est = intcptprophf(i);
  intcptprophfCI(i,:) = confidence_interval_general(est,interceptsehf,length(x)-2,95);
  sehf(i, :) = [gofhfprop.rmse slopesehf interceptsehf];
  
  x = mighfdum(:,seccol)/3600;
  y = mighfdum(:,1);
  x_bar = wt_mean(x,wthf);
  y_bar = wt_mean(y,wthf);
  x_var = sum(wthf.*(x-x_bar).^2) / sum(wthf);
  y_var = sum(wthf.*(y-y_bar).^2) / sum(wthf);
  xy_cov = sum(wthf.*(x-x_bar).*(y-y_bar)) / sum(wthf);
  pearwthf(i) = xy_cov / sqrt(x_var*y_var);
  
  text(f.ax(3),0.45,0.2,sprintf('Slope: %.1f km/h',slopeprophf(i)),'FontSize',8,...
    'unit','normalized','horizontalalignment','left');
  text(f.ax(3),0.97,0.2,sprintf('SE: %.2f',slopesehf),'FontSize',8,...
    'unit','normalized','horizontalalignment','right');
  text(f.ax(3),0.45,0.13,sprintf('Pearson: %.3f',pearwthf(i)),'FontSize',8,...
    'unit','normalized','horizontalalignment','left');
  hold(f.ax(3),'off');
  
  % subplot 4 of figure i
  hold(f.ax(4),'on');
  f.ax(4).FontSize = 9;
  f.ax(4).Box = 'on';
  yyaxis(f.ax(4),'left');
  plot(f.ax(4),angle,rmsehf(i,:),'o-','linew',1,'color','k','markers',4);
  yranl = f.ax(4).YLim;
    plot(f.ax(4),[angrmsehf(i) angrmsehf(i)],[0 100],'--','linew',0.5,'color','k');
  f.ax(4).YColor = 'k';
  f.ax(4).YLim = yranl;
  ylabel(f.ax(4),'RMSE of HF','fontsize',11);
  
  yyaxis(f.ax(4),'right');
  plot(f.ax(4),angle,slopehf(i,:),'^:','linew',1,'color',[0.12 0.56 1],'markers',4);
  yranr = f.ax(4).YLim;
  plot(f.ax(4),[angslopehf(i) angslopehf(i)],f.ax(4).YLim,'--','linew',0.5,...
    'color',[0.12 0.56 1]);
  f.ax(4).YColor = [0.12 0.56 1]; %[0 0 0.5];
  f.ax(4).YLim = yranr;
  ylabel(f.ax(4),'Slope of HF (km/h)','fontsize',11);
  text(f.ax(4),0.04,0.91,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
  xlim(f.ax(4),[0 360]);
  %     ymax = f.ax(4).YLim(2)+0.1;
  %     ylim(f.ax(4),[0 ymax]);
  %     ylim(f.ax(4),[0 1]);
  xlabel(f.ax(4),'Trial propagation direction (deg)','fontsize',11);
  hold(f.ax(4),'off');

end


%% estimate the time needed to pass through 002 ellipse
%angle between propagation direction and semi-minor axis
theta = median(angrmsehf)-(45+90);
% theta = 0;
a = 1.75;
b = 1.0;
%the length of propagation vector on the ellipse, using elliptical equation: (x/a)^2+(y/b)^2=1;
len = 2* sqrt(1/(sind(theta)^2/a.^2+cosd(theta)^2/b.^2));

%take 20 km/h as the average propagation speed
speed = 20;
time = len/speed *3600;







