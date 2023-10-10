% function synthrccvssrcregsize.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is similar in large part to 'synthrccvssatur.m', especially 
% for the start. The goal is to analysize the variation of CC (median of 
% running CC and overall CC after best alignment), with respect to the change
% of source region size. 
%
% --Different from 'synthrccvssatur.m', we are NOT 
% using the source distribution from the observation as in 'locinterp002_4s.m',
% instead, we will use a 2D Gaussian distribution with different spreadings 
% (varying sigmas). In this case, we will make the sigma_x same as sigma_y, 
% which will give circles of different sizes. We will also truncate the PDF
% at +-3 sigma, then concatenate it with zeros outside. But that requires some
% some adjustment on PDF value inside to ensure the sum is 1.
% 
% --The way to generate a custom shape of 2D Gaussian (normal) PDF can be found
% in 'testmvnpdf.m' and 'testrand_custompdf_p1.m'.
%
% --If you instead use a uniform random distribution for sources, then you 
% should use code 'analyze_synth'. 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/05/05
% Last modified date:   2023/09/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
clear
clc
% close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');

[scrsz, resol] = pixelperinch(1);

%% prepare templates (Green's functions)
ccstack = [];
sps = 160;
templensec = 60;

fam = '002';
disp(fam);

stas=['PGC ';
      'SSIB';
      'SILB';
      ];
nsta=size(stas,1);

for ista = 1: nsta
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
        num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao_catnew');
    ccstack(:,ista) = load(fname);
end
STA = ccstack;

% for ista=1:nsta
%     STA(:,ista)=Bandpass(STA(:,ista),sps,0.1,15,2,2,'butter');   % change 'bandpass' to 'Bandpass'
% end
%plot the raw templates, not filtered, not best aligned
% figure
% hold on
% plot(STA(:,1),'r')
% plot(STA(:,2),'b')
% plot(STA(:,3),'k')

%%%The below aligns the templates by x-correlation
[maxses,imaxses]=max(STA,[],1);
[minses,iminses]=min(STA,[],1);
spread=maxses-minses;
% zcrosses=round(0.5*(imaxses+iminses));  % rough, assuming symmetry, Chao 2021/07/16
%automatically find the zero-crossings
zcrosses = zeros(nsta,1);
for ista = 1:nsta
    seg = STA(iminses(ista): imaxses(ista),ista);  % for zero-crossing timing, only use the main station
    [~,zcrosses(ista)] = min(abs(seg));
    zcrosses(ista) = zcrosses(ista)-1+iminses(ista);  % convert to global index
end
%now we want to cut a segment around the zero-crossing at each station
sampbef=6*sps;
sampaft=10*sps;
is=zcrosses-sampbef;
ie=zcrosses+sampaft;
for ista=1:nsta
    STAtmp(:,ista)=STA(is(ista):ie(ista),ista);  % this means templates are 'aligned' at zero-crossings
end
%x-correlation independently between each station pair 
mshiftadd=10*sps/40;
tempxc(:,1)=xcorr(STAtmp(:,2),STAtmp(:,3),mshiftadd,'coeff');
tempxc(:,2)=xcorr(STAtmp(:,1),STAtmp(:,2),mshiftadd,'coeff'); %shift STAtmp(3,:) to right for positive values
tempxc(:,3)=xcorr(STAtmp(:,1),STAtmp(:,3),mshiftadd,'coeff'); %shift STAtmp(2,:) to right for positive values
[~,imax]=max(tempxc,[],1);
imax=imax-(mshiftadd+1); %This would produce a slightly different shift, if filtered seisms were used.
imax(2)-imax(3)+imax(1)   %enclosed if it equals to 0
for ista=2:nsta
    STAtmp(mshiftadd+1:end-(mshiftadd+1),ista)=STAtmp(mshiftadd+1-imax(ista):end-(mshiftadd+1)-imax(ista),ista);
end
for ista=1:nsta
    STAtmp(:,ista)=STAtmp(:,ista)/spread(ista); % now templates are 'aligned' indeoendently by x-corr wrt. sta 1
end
% figure
% hold on
% plot(STAtmp(:,1),'r')
% plot(STAtmp(:,2),'b')
% plot(STAtmp(:,3),'k')
%%%The above aligns the templates by x-correlation

%%%detrend, taper and bandpass templates
tmpwlet = STAtmp; % no bandpass
tmpwletf = STAtmp;  % bandpassed version
fractap = sps/size(tmpwlet,1);
for ista = 1: nsta
  %romve mean, linear trend of template
  tmpwlet(:,ista) = detrend(tmpwlet(:,ista));
  %and taper with tukeywin, which is actually a tapered cosine window
  w = tukeywin(size(tmpwlet(:,ista),1),fractap);
  tmpwlet(:,ista) = w.* tmpwlet(:,ista);
  %detrend again for caution
  tmpwlet(:,ista)=detrend(tmpwlet(:,ista));
  %filter the template
  hiwlet=18;
  lowlet=1.8;
  tmpwletf(:,ista) = Bandpass(tmpwlet(:,ista), sps, lowlet, hiwlet, 2, 2, 'butter');
  %detrend again for caution
  tmpwletf(:,ista)=detrend(tmpwletf(:,ista));
end

%%%constrained CC, so that only 2 offsets are independent
ccmid = round(size(tmpwletf,1)/2);
ccwlen = 10*sps;
loffmax = 5*sps/40;
ccmin = 0.01;  % depending on the length of trace, cc could be very low
iup = 1;    % times of upsampling
[off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(tmpwletf',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
offwlet1i(1) = 0;
offwlet1i(2) = round(off12con);
offwlet1i(3) = round(off13con);

%%%automatically find the rough zero-crossing time, whose abs. value is closest to 0, whether + or -
[~,imin] = min(tmpwletf(:,1));
[~,imax] = max(tmpwletf(:,1));
[~,zcsta1] = min(abs(tmpwletf(imin:imax,1)));
zcsta1 = zcsta1+imin-1;
greenlen = pow2(9)*sps/40;
green = zeros(greenlen,nsta); % no bandpass
greenf = zeros(greenlen,nsta);  % bandpassed version
for ista = 1: nsta
  %cut according to the zero-crossing and the time shift from the constrained CC
  green(:,ista) = tmpwlet(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenf(:,ista) = tmpwletf(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  %detrend again for caution
  green(:,ista)=detrend(green(:,ista));
  greenf(:,ista)=detrend(greenf(:,ista));
  %normalize by max amp
  green(:,ista)=green(:,ista)/max(abs(green(:,ista)));    % normalize
  greenf(:,ista)=greenf(:,ista)/max(abs(green(:,ista)));    % normalize
  
  %re-find the zero-crossing as the template length has changed
  [~,imin] = min(greenf(:,ista));
  [~,imax] = max(greenf(:,ista));
  [~,zcrosses(ista)] = min(abs(greenf(imin:imax,ista)));
  zcrosses(ista) = zcrosses(ista)+imin-1;
end
%the following is just a check, because now the templates must be best aligned 
ccmid = round(size(greenf,1)/2);
ccwlen = 4*sps;
loffmax = 5*sps/40;
ccmin = 0.01;  % depending on the length of trace, cc could be very low
iup = 1;    % times of upsampling
[off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(greenf',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
if ~(off12con==0 && off13con==0)
  disp('Filtered templates are NOT best aligned');
end

%%%plot the unfiltered and filtered templates
mean(green,1)
figure
subplot(2,1,1)
hold on
plot(green(:,1),'r')
plot(green(:,2),'b')
plot(green(:,3),'k')
text(0.95,0.9,'Raw','Units','normalized','HorizontalAlignment',...
  'right');
mx=max([abs(green(:,1)); abs(green(:,2))]);
xlim([0 greenlen])
ylim([-mx mx])
title("Green's functions");
box on

subplot(2,1,2)
hold on
plot(greenf(:,1),'r')
plot(greenf(:,2),'b')
plot(greenf(:,3),'k')
text(0.95,0.9,sprintf('%.1f-%.1f Hz',lowlet,hiwlet),'Units','normalized','HorizontalAlignment',...
  'right');
mx=max([abs(greenf(:,1)); abs(greenf(:,2))]);

%%%running CC using a window length of 'cclen'
mwlen=sps/2;
[ircc,rcc12] = RunningCC(greenf(:,1), greenf(:,2), mwlen);
[~,rcc13] = RunningCC(greenf(:,1), greenf(:,3), mwlen);
[~,rcc23] = RunningCC(greenf(:,2), greenf(:,3), mwlen);
rcc = (rcc12+rcc13+rcc23)/3;
%alln(alln<0)=-10^4*yma; %just so they don't plot.
% plot(samples(cclen/2+1:greenlen-cclen/2),mx*alln,'co','markersize',2); % scale with data amp.
plot(ircc,mx*rcc,'co','markersize',2); % scale with data amp.
xlim([0 greenlen])
ylim([-mx mx])
box on
xc23=xcorr(greenf(:,2),greenf(:,3),10,'coeff');
xc13=xcorr(greenf(:,1),greenf(:,3),10,'coeff');
xc12=xcorr(greenf(:,1),greenf(:,2),10,'coeff');
[ccmax23,imax23]=max(xc23)
[ccmax13,imax13]=max(xc13)
[ccmax12,imax12]=max(xc12)


%% design a custom 2D normal PDF using 'mvnpdf'
%%%From the true source distribution and the transform from offset to space, we roughly know that
%%%a ellipse with semia=1.75km, semib=1.25km corresponds to +-12 and +-15 samples in off12 and
%%%off13. So we can set a range of max offsets

xvec = -48:1:48;
yvec = -48:1:48;
[X1,X2] = meshgrid(xvec,yvec);
X = [X1(:) X2(:)];

sigmatry = (1:1:15)'.^2;

mrcc = zeros(length(sigmatry), 1);
mcc = zeros(length(sigmatry), 1);
nsrc = zeros(length(sigmatry), 1);
ntot = zeros(length(sigmatry), 1);

for i = 1: length(sigmatry)
  % i = 1;
  %use the designed covariance matrix to create the PDF of the 2D normal distribution
%   mu = [0 0];
%   sigma = sigmatry(i).*[1 0; 0 1]; % off-diagonal=0 --> circular
  mu = [2 -6];
  sigma = sigmatry(i).*[1 0.6; 0.6 1]; % off-diagonal=0.6 --> 45 deg ellipse
  [eigvec,eigval] = eig(sigma);
  sigmaeig = [sqrt(eigval(1,1)) sqrt(eigval(2,2))];
  
  pdffull = mvnpdf(X,mu,sigma);
  
  figure
  imagesc(xvec,yvec,reshape(pdffull,length(yvec),length(xvec))); hold on;% substitute for surf, but need to make 'y' axis direction as normal
  ax = gca; ax.YDir = 'normal';
  % contour(x1,x2,y,'k-');
  % caxis(ax,[min(y(:))-0.5*range(y(:)),max(y(:))])
  axis equal tight
  % xlim(minmax(x1));
  % ylim(minmax(x2));
  colormap('jet');
  xlabel('x1');
  ylabel('x2');
  c=colorbar;
  c.Label.String = 'Probability Density';
  
  %plot the horizontal and vectical boundary similar to limiting offset12, 13 and 23
  loff_max = round(3*sqrt(sigmatry(i)));  % maximum allowed shift is 4 samples at 40 Hz, it is also true from the synthetic distribution
  offxran = [-loff_max+mu(1) loff_max+mu(1)];
  offyran = [-loff_max+mu(2) loff_max+mu(2)];
  off12m = round(range(offxran)/2);
  off13m = round(range(offyran)/2);
  xbndlow = [offxran(1):offxran(2), offxran(2)*ones(1,off13m)];
  xbndupp = [offxran(2)-1:-1:offxran(1), offxran(1)*ones(1,off13m)];
  ybndlow = [offyran(1)*ones(1, off12m), offyran(1):offyran(2)];
  ybndupp = [offyran(2)*ones(1, off12m-1), offyran(2):-1:offyran(1)];
  xbnd = [xbndlow xbndupp]';
  ybnd = [ybndlow ybndupp]';
  plot(xbnd,ybnd,'-','Color','w','linew',2);
  
  %plot 1-sigma, 2-sigma and 3-sigma ellipse
  x0 = mu(1);
  y0 = mu(2);
  for j = 1:2
    angle(j) = atan2d(eigvec(2,j),eigvec(1,j));
  end
%   angrot = 0;
  angrot = min(angle);
  semia = 2*sigmaeig(2);  %2-sigma
  semib = 2*sigmaeig(1);
  [x, y] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);
  plot(x,y,'k-','linew',1);
  semia = 3*sigmaeig(2);  %3-sigma
  semib = 3*sigmaeig(1);
  [xcut, ycut] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);
  plot(xcut,ycut,'k-','linew',2);
  
  %truncate the PDF by 3-sigma
  bnd = [xcut ycut];
  [iin,ion] = inpolygon(X(:,1),X(:,2),bnd(:,1),bnd(:,2));
  isinbnd = iin | ion;
  Xtr = X(isinbnd == 1, :);
  pdftr = pdffull(isinbnd == 1, :);
  
  % %check the truncation result
  % figure
  % scatter(Xtr(:,1),Xtr(:,2),50,pdftr,'s','filled');hold on;
  % axis equal
  % % axis([-10 10 -10 10])
  % xlim(minmax(xvec));
  % ylim(minmax(yvec));
  % colormap('jet');
  % xlabel('x1');
  % ylabel('x2');
  % c=colorbar;
  % c.Label.String = 'Probability Density';
  % plot(xcut,ycut,'k-');
  
  %Compose the PDF
  pdfoff = [Xtr pdftr];
  %because of truncation, the sum is NOT 1, so adjust it
  pdfoff(:,3) = pdfoff(:,3)/sum(pdfoff(:,3));
  
  %%% Note that 'pinky(x1,x2,y)' can only deal with an evenly-spaced vector x1 and x2, y is a grid of
  %%% PDF values on the grid points that are defined by x1 and x2.
  %%% Therefore, 'pinky(x1,x2,y)' can be applied when either your PDF is generated by binning upon an
  %%% evenly-spaced grid (or ksdensity), OR, you bin by pixel in locations but return to the sample
  %%% domain, in which case the grid is even. But the latter case needs zero-padding to grid points
  %%% that is otherwise empty.
  [pdfoffpad,xgrid,ygrid,pdfoffgrid,ind1] = zeropadmat2d(pdfoff,xvec,yvec);
  
  %%%plot the PDF, note the color should be the SAME as density, as it is just a normalization
  figure
  dum = pdfoffpad;
  dum(dum(:,3)~=0, :) = [];
  scatter(dum(:,1),dum(:,2),6,dum(:,3),'linew',0.2);  hold on
  dum = pdfoffpad;
  dum(dum(:,3)==0, :) = [];
  scatter(dum(:,1),dum(:,2),30,dum(:,3),'filled');
  plot(xbnd,ybnd,'-','Color','w','linew',2);
  axis equal
  axis([minmax(xvec) minmax(yvec)]);
  colormap(jet)
  c=colorbar;
  c.Label.String = 'Probability Density';
  caxis([0 max(pdfoffpad(:,3))]);
  xlabel(strcat({'Offset 12 (samples at '},num2str(sps),' sps)'));
  ylabel(strcat({'Offset 13 (samples at '},num2str(sps),' sps)'));
  title('Data PDF in smaple space');
  box on; grid on;
  
  
  %% synthetics generation
  %%%Specify the amplitude-frequency (counts) distribution
  distr='UN'  % uniform distribution
  
  Twin=250+3+2*ceil(greenlen/sps); %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
  winlen=Twin*sps+1;
  skiplen=greenlen;
  tdura = 0.5;  % duration from the broadband template
  % tdura = 0.75; % for bandpassed, can be ~0.3s or 0.7s depending on the definition of 'duration'
  satn=1/tdura*Twin   % if just saturated, how many templates can be fit in? a single peak is ~20 samples wide; maybe a little less (at 100 sps).
  %Twin is window duration in seconds. Events can fall within Twin of
  %the start, but the synthetics will go to Twin*(sample rate)+Greenlen to
  %avoid checking for subscript overrrun.  When "synths" are written from "synth" in the
  %subroutine, Greenlen from the start and end will not be written.
  fracelsew=0; %0.25 %0.5; %0.6; %The ratio of elsewhere events to local events.  Zero for Discrete Ide.
  
  %here, i want to ensure for same length of sysnthetics and region of different sizes, the spatial
  %'density' of sources are roughly the same, so for a larger source region, synthetic waveform will be
  %more saturated
%   nsat(i) = pi*(3*sqrt(sigmatry(i)))^2 *1 ./(pi*(3*sqrt(sigmatry(1)))^2);  % times of saturation
  %Alternatively, use the same saturation level for all region size, so the larger region will have
  %sparser sources.
  nsat(i) = 0.82;
  
  writes=round(nsat(i)*satn) %how many templates to throw in, under different degrees of saturation
  rng('default');
  %%%Here the noise is 'uniform' in time!
  % mampnoi = 1e-3;
  % synth=mampnoi*2*(rand(winlen+greenlen+2*10,nsta)-0.5); %+2*10 a little extra, for jiggering 2nd & 3rd stations.% for ista=1:nsta
  %%%Below use the Gaussian white noise
  noistd = 5e-2;
  synth=noistd*(randn(winlen+greenlen+2*10,nsta)); %+2*10 a little extra, for jiggering 2nd & 3rd stations.% for ista=1:nsta
  
  nouts=length(writes);
  seed=round(writes(1)/5e3); %for random number generator
  
  if strcmp(distr,'PL') || strcmp(distr,'UN') %b>150 for uniform
    b=999; %>150 for uniform size distribution
    [synths,mommax,sources,greensts]=synthgen(writes,winlen,skiplen,synth,greenf,b,xgrid,...
      ygrid,pdfoffgrid,sps,fracelsew,seed);
  end
  
  
  %% median(rcc) or (cc) VS. saturation rate
  % wlensectry = (15:30:300)';
  
  insat = 1;
  %   disp(nsat(insat));
  %   for iwlen = 1: length(wlensectry)
  %     wlensec = wlensectry(iwlen);
  %     iwlen = 1;
  wlensec = 250; %this is the median win length of from burst grouping result in 'tremorbursts002_4s.m'
  bufsec = 1; %need some buffer window for later CC alignment
  buffer = bufsec*sps;
  wlensecb = wlensec+bufsec;
  
  sigstab = [];  %target window of synthetic signals at all stations
  sigpnstab = [];  %target window of synthetic signals plus the noise
  noistab = [];  %target window of noise
  
  srcpairstab = cell(nsta,1);  %sources whose zero-crossing with time range, [indtarvl rnoff12 rnoff13 amp]
  
  for ista = 1: nsta
    % ista = 3;
    
    tgreen = zcrosses(ista);  % anchor time of templates, here choose as zero-crossing time
    indstsig = 1;  % starting index of the simulated signal to test
    indedsig = wlensecb*sps+indstsig-1; % ending index of the simulated signal to test
    source = sources(1:size(greensts{insat}{1},1), :);  % [indtori indttrvl indtarvl rnoff12 rnoff13 amp];
    if ista == 1
      greenst = greensts{insat}{1}; % the starting index of each added template, context of full length
    else
      %ind - rnoff is the arrival time in index at sta 2 and 3
      %note the sign here, if off12 >0, move 2 to the right to align with 1, meaning that 2 is
      %earlier than 1, ie., tarvl2 < tarvl1. Be consistent all the time
      greenst = greensts{insat}{1}-source(:, 2+ista);
    end
    
    %%%you don't need all impulses, only some of them contribute to the length of truncated record
    %%%cut out the green indice and sources that contribute
    %%%Note that some zero-crossings might appear later than 'indedsig', but the corresponding start
    %%%index of the green is actually smaller than 'indedsig', ie. their 'greenst' <= 'indedsig'
    induse = find(greenst>=skiplen-greenlen+indstsig & greenst<=indedsig+skiplen-1);
    greenst = greenst(induse);
    source = source(induse,:);
    
    greenzc = greenst+tgreen; % index of approximate zero-crossing
    source(:, 3) = source(:, 3)+tgreen;  % now 'greenzc' should be the same as 3rd col of 'source'
    impamp = zeros(max(greenzc)+20,1);
    for ii = 1: length(greenzc)
      impamp(greenzc(ii)) = impamp(greenzc(ii))+source(ii, 6);
    end
    
    %note that some length of the full simulation 'skiplen' was skipped in 'synthgen.m'
    sigpntemp = squeeze(synths(:,ista,insat));  % simulated signal with white noise
    noitemp = synth(skiplen:winlen,ista);   % account for the skipped part
    sigtemp = sigpntemp-noitemp; % subtract the white noise
    sigpnstab(:,ista) = sigpntemp(indstsig:indedsig); % now focus on the part of interest
    sigstab(:,ista) = sigtemp(indstsig:indedsig);
    noistab(:,ista) = noitemp(indstsig:indedsig);
    
    %now 'source' contain impulse triplets that contribute to the trace segment of interest, but we
    %want to get ones whose zero-crossing falls into the window.
    source(:,3) = source(:,3)-skiplen;  % index after skipping
    srcpairb = source(source(:,3)>=indstsig & source(:,3)<=indedsig, 3:6);
    srcpairb = sortrows(srcpairb,1,'ascend');
    srcpairstab{ista} = srcpairb; %ideally for each station this should be the same, but coincidence is possible
    
  end
  
  %notify if sources are the same at diff stations
  if ~isequaln(srcpairstab{1,1},srcpairstab{2,1}) || ~isequaln(srcpairstab{1,1},srcpairstab{3,1})
    disp('Extracted sources at different stations are not the same, re-check!');
  end
  
  
  %% Best alignment for the testing window
  ccmid = round(size(sigpnstab,1)/2);
  ccwlen = round(size(sigpnstab,1)*0.8);
  loffmax = 5*sps/40;
  ccmin = 0.01;  % depending on the length of trace, cc could be very low
  iup = 1;    % times of upsampling
  mshiftadd=10*sps/40;
  [off12con,off13con,cc,iloopoff,loopoff] = constrained_cc_interp(sigpnstab',ccmid,...
    ccwlen,mshiftadd,loffmax,ccmin,iup);
  off1i = zeros(nsta,1);
  off1i(2) = round(off12con);
  off1i(3) = round(off13con);
  
  %align the signals, noises, etc
  sigpnsta = zeros(wlensec*sps, nsta);
  noista = zeros(wlensec*sps, nsta);
  sigsta = zeros(wlensec*sps, nsta);
  for ista = 1: nsta
    sigpnsta(:,ista) = sigpnstab(round(buffer/2)+1-off1i(ista): end-round(buffer/2)-off1i(ista), ista);
    noista(:,ista) = noistab(round(buffer/2)+1-off1i(ista): end-round(buffer/2)-off1i(ista), ista);
    sigsta(:,ista) = sigstab(round(buffer/2)+1-off1i(ista): end-round(buffer/2)-off1i(ista), ista); % sta 2
  end
  
  srcpairb = srcpairstab{1,1};
  %store the info of synthetic impulses, same format as in paired deconvolution,
  %9 cols, [ind1 amp1 ind2 amp2 ind3 amp3 off12 off13 off23]
  srcpair = zeros(size(srcpairb,1), 9);
  srcpair(:,1) = srcpairb(:,1)-round(buffer/2);  % cut out the buffer
  srcpair(:,[2 4 6]) = repmat(srcpairb(:,4), 1,3);
  %the indices of synthetic impulses need to be shifted too
  srcpair(:,3) = srcpair(:,1)-srcpairb(:,2)-off1i(2); % note the sign is consistent!
  srcpair(:,5) = srcpair(:,1)-srcpairb(:,3)-off1i(3);
  srcpair(:,7) = srcpairb(:,2);  % no need to shift offset, because they are 'true' offset
  srcpair(:,8) = srcpairb(:,3);
  srcpair(:,9) = srcpair(:,8)-srcpair(:,7);  % off23, == off13 - off12 == tarvl2 - tarvl3
  
  %compute running CC
  mwlen=sps/2;
  % mwlen=sps;
  [ircc,rcc12] = RunningCC(sigpnsta(:,1), sigpnsta(:,2), mwlen);
  [~,rcc13] = RunningCC(sigpnsta(:,1), sigpnsta(:,3), mwlen);
  [~,rcc23] = RunningCC(sigpnsta(:,2), sigpnsta(:,3), mwlen);
  rcc = (rcc12+rcc13+rcc23)/3;  %average
  mrcc(i) = median(rcc); %
  
  %compute the median absolute amplitude and envelope of the same moving window
  %for the moving window at the same station, sensable to use median
  [ir,ramp1,renv1] = Runningampenv(sigpnsta(:,1),mwlen,mwlen-1,'median');
  [~,ramp2,renv2] = Runningampenv(sigpnsta(:,2),mwlen,mwlen-1,'median');
  [~,ramp3,renv3] = Runningampenv(sigpnsta(:,3),mwlen,mwlen-1,'median');
  %looks like using the amplitude and envelope are pretty similar
  %maybe better to mean, otherwise stations PGC will be underrepresented due to a larger amp
  %variation
  ramp = mean([ramp1 ramp2 ramp3], 2);  % use mean or median??
  renv = mean([renv1 renv2 renv3], 2);
  
  %overall average max CC based on the current best alignment
  mcc12 = sum(sigpnsta(:,1).*sigpnsta(:,2))./ ...
    (sqrt(sum(sigpnsta(:,1).^2)).*sqrt(sum(sigpnsta(:,2).^2)));
% mcc12 = xcorr(sigpnsta(:,1), sigpnsta(:,2),0,'normalized');   % same result, might be more efficient
  mcc13 = sum(sigpnsta(:,1).*sigpnsta(:,3))./ ...
    (sqrt(sum(sigpnsta(:,1).^2)).*sqrt(sum(sigpnsta(:,3).^2)));
  mcc23 = sum(sigpnsta(:,2).*sigpnsta(:,3))./ ...
    (sqrt(sum(sigpnsta(:,2).^2)).*sqrt(sum(sigpnsta(:,3).^2)));
  mcc(i) = (mcc12+mcc13+mcc23)/3;
  
  %number of sources actually of interest
  nsrc(i) = size(srcpair,1);
  ntot(i) = writes;
  
  %     %ratio of number of sources to size of region
  %     nrat = nsrc(i) / ( pi*(3*sqrt(sigmatry(i)))^2 );  %./ (pi*(3*sqrt(sigmatry(1)))^2) .* satn
  
  %distribution of synthetic sources, ground truth
  xran = [-loff_max+mu(1)-1 loff_max+mu(1)+1];
  yran = [-loff_max+mu(2)-1 loff_max+mu(2)+1];
  lsig = size(sigpnsta,1);
  cran = [0 lsig];
  f1.fig = figure;
  f1.fig.Renderer = 'painters';
  ax1=gca;
  [ax1] = plt_decon_imp_scatter(ax1,srcpair,xran,yran,cran,offxran,offyran,...
    sps,75,'mean','tori');
  scatter(ax1,off1i(2),off1i(3),20,'ks','filled','MarkerEdgeColor','k');
  title(ax1,sprintf('Synthetic sources, using data of %d Hz',sps));
  %
  %     %%% scatter of rela locations, and account for prealignment offset
  %     xran = [-3 3];
  %     yran = [-3 3];
  %     cran = [0 lsig];
  %     ftrans = 'interpchao';
  %     [f] = plt_decon_imp_scatter_space(srcpair,xran,yran,cran,offxran,offyran,sps,50,ftrans,'mean','tori');
  %     plot(gca,xcut,ycut,'k-','linew',2);
  %     title(gca,sprintf('Synthetic sources, radius of %d samples', 3*sqrt(sigmatry(i))));
  
end

%%
%%%For each window length, plot of median(rcc)/mcc  vs.  saturation rate
close all

i = 1;
figure
color = jet(length(sigmatry));
ax(1) = subplot(221); hold on
plot(3*sqrt(sigmatry),mcc,'o-','markersize',4,'Color',color(i,:));
for i = 1: length(sigmatry)
  text(3*sqrt(sigmatry(i)), mcc(i)*1.02, num2str(nsrc(i)), 'FontSize',8);
end
% legend(num2str(sigmatry),'Location','southeastx');
xlabel('Radius of source region');
ylabel('Overall CC');
box on; grid on

ax(2) = subplot(222); hold on
plot(3*sqrt(sigmatry),mrcc,'o-','markersize',4,'Color',color(i,:));
for i = 1: length(sigmatry)
  text(3*sqrt(sigmatry(i)), mrcc(i)*1.02, num2str(nsrc(i)), 'FontSize',8);
end
xlabel('Radius of source region');
ylabel('Median of running CC');
box on; grid on

ax(3) = subplot(223); hold on
% for i = 1: length(sigmatry)
  plot(nsat,nsrc./(pi*(3*sqrt(sigmatry)).^2),'o-','markersize',4,'Color',color(i,:));
  plot(nsat,ntot./(pi*(3*sqrt(sigmatry)).^2),'o-','markersize',4,'Color',color(end,:));
% end
xlabel('Saturation level');
ylabel('Number of sources inside window / area of source region');
box on; grid on

supertit(ax([1 2]),...
  sprintf('CC vs. region radius for Synthetic record length: %d s / %d s',wlensec, Twin),12);












