% testconvolvegaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
% This function is to test the convolution between the LFE template and an
% Gaussian function / heaviside step function / box car function, etc., to 
% see how the shape would change.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/11/24
% Last modified date:   2022/11/24 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, resol] = pixelperinch(1);

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');
  
fam = '002';   % family number
stas=['PGC  '
  'SSIB '
  'SILB '
  'LZB  '
  'TWKB '
  'MGCB '
  'KLNB ']; % determine the trio and order, here the 1st sta is PGC
nsta=size(stas,1);         %  number of stations



%% prepare templates (Green's functions), from 'lfetemp002_160sps.m'
sps = 160;
templensec = 60;

ccstack = [];
for ista = 1: nsta
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
        num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao_catnew');
    ccstack(:,ista) = load(fname);
end
STA = detrend(ccstack);

ccstackort = [];
for ista = 1: nsta
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
        num2str(templensec), 's_', 'BBCCS_', 'ort_Nof_Non_Chao_catnew');
    ccstackort(:,ista) = load(fname);
end
STAort = detrend(ccstackort);

%flag of normalization
normflag = 0;

% %plot the raw templates, not filtered, not best aligned
% figure
% subplot(211)
% hold on
% for ista = 1: nsta
%   plot(STA(:,ista));
% end
% subplot(212)
% hold on
% for ista = 1: nsta
%   plot(STAort(:,ista));
% end

%%%The below aligns the templates by x-correlation
ist = templensec*sps*4/10;  %not using whole window in case any station has very long-period energy
ied = templensec*sps*6/10;
[maxses,imaxses]=max(STA(ist:ied,:),[],1);
[minses,iminses]=min(STA(ist:ied,:),[],1);
spread=maxses-minses;
imaxses = imaxses+ist-1;  %convert to global indices
iminses = iminses+ist-1;
% zcrosses=round(0.5*(imaxses+iminses));  % rough, assuming symmetry, Chao 2021/07/16
%automatically find the zero-crossings
zcrosses = zeros(nsta,1);
for ista = 1:nsta
    seg = detrend(STA(iminses(ista): imaxses(ista),ista));  % for zero-crossing timing, only use the main station
    [~,zcrosses(ista)] = min(abs(seg));
    zcrosses(ista) = zcrosses(ista)-1+iminses(ista);  % convert to global index
end
%now we want to cut a segment around the zero-crossing at each station
sampbef=6*sps;
sampaft=10*sps;
is=zcrosses-sampbef;
ie=zcrosses+sampaft;
for ista=1:nsta
    STAtmp(:,ista)=detrend(STA(is(ista):ie(ista),ista));  % this means templates are 'aligned' at zero-crossings
    STAtmport(:,ista)=detrend(STAort(is(ista):ie(ista),ista)); 
end
%x-correlation independently between each station pair 
mshiftadd=10*sps/40;
tempxc(:,1)=xcorr(STAtmp(:,2),STAtmp(:,3),mshiftadd,'coeff');
tempxc(:,2)=xcorr(STAtmp(:,1),STAtmp(:,2),mshiftadd,'coeff'); %shift STAtmp(3,:) to right for positive values
tempxc(:,3)=xcorr(STAtmp(:,1),STAtmp(:,3),mshiftadd,'coeff'); %shift STAtmp(2,:) to right for positive values
for ista = 4: nsta
  tempxc(:,ista)=xcorr(STAtmp(:,1),STAtmp(:,ista),mshiftadd,'coeff'); %shift STAtmp(2,:) to right for positive values
end
[~,imax]=max(tempxc,[],1);
imax=imax-(mshiftadd+1); %This would produce a slightly different shift, if filtered seisms were used.
imax(2)-imax(3)+imax(1);   %enclosed if it equals to 0
for ista=2:nsta
    STAtmp(mshiftadd+1:end-(mshiftadd+1),ista)=detrend(STAtmp(mshiftadd+1-imax(ista):end-(mshiftadd+1)-imax(ista),ista));
    STAtmport(mshiftadd+1:end-(mshiftadd+1),ista)=detrend(STAtmport(mshiftadd+1-imax(ista):end-(mshiftadd+1)-imax(ista),ista));
end
%normalization
if normflag 
  for ista=1:nsta
      STAtmp(:,ista)=STAtmp(:,ista)/spread(ista); % now templates are 'aligned' indeoendently by x-corr wrt. sta 1
      STAtmport(:,ista)=STAtmport(:,ista)/spread(ista);
  end
end
% figure
% subplot(211)
% hold on
% for ista = 1: nsta
%   plot(STAtmp(:,ista));
% end
% subplot(212)
% hold on
% for ista = 1: nsta
%   plot(STAtmport(:,ista));
% end
%%%The above aligns the templates by x-correlation

%%%detrend, taper and bandpass templates
tmpwlet = STAtmp; % no bandpass
tmpwletf = STAtmp;  % bandpassed version
fractap = sps/size(tmpwlet,1);
tmpwletort = STAtmport; % no bandpass
tmpwletfort = STAtmport;  % bandpassed version
for ista = 1: nsta
  %romve mean, linear trend of template
  tmpwlet(:,ista) = detrend(tmpwlet(:,ista));
  %and taper with tukeywin, which is actually a tapered cosine window
  w = tukeywin(size(tmpwlet(:,ista),1),fractap);
  tmpwlet(:,ista) = w.* tmpwlet(:,ista);
  %detrend again for caution
  tmpwlet(:,ista) = detrend(tmpwlet(:,ista));
  %filter the template
  hiwlet=18;
  lowlet=1.8;
  tmpwletf(:,ista) = Bandpass(tmpwlet(:,ista), sps, lowlet, hiwlet, 2, 2, 'butter');
  %detrend again for caution
  tmpwletf(:,ista) = detrend(tmpwletf(:,ista));
  
  %same process for orthogonal
  tmpwletort(:,ista) = detrend(tmpwletort(:,ista));
  tmpwletort(:,ista) = w.* tmpwletort(:,ista);
  tmpwletort(:,ista) = detrend(tmpwletort(:,ista));
  tmpwletfort(:,ista) = Bandpass(tmpwletort(:,ista), sps, lowlet, hiwlet, 2, 2, 'butter');
  tmpwletfort(:,ista) = detrend(tmpwletfort(:,ista));
end

%%%constrained CC, so that only 2 offsets are independent
ccmid = round(size(tmpwletf,1)/2);
ccwlen = 10*sps;
loffmax = 5*sps/40;
ccmin = 0.01;  % depending on the length of trace, cc could be very low
iup = 1;    % times of upsampling
[off12con,off13con,cc] = constrained_cc_interp(tmpwletf(:,1:3)',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
offwlet1i(1) = 0;
offwlet1i(2) = round(off12con);
offwlet1i(3) = round(off13con);

for ista = 4: nsta
  [coef,lag] = xcorr(tmpwletf(:,1), tmpwletf(:,ista), mshiftadd, 'coeff');
  [mcoef, idx] = max(coef);   % max of master raw cc
  offwlet1i(ista) = lag(idx);   % offset in samples  
end

%%%automatically find the rough zero-crossing time, whose abs. value is closest to 0, whether + or -
[~,imin] = min(tmpwletf(:,1));
[~,imax] = max(tmpwletf(:,1));
[~,zcsta1] = min(abs(tmpwletf(imin:imax,1)));
zcsta1 = zcsta1+imin-1;
greenlen = pow2(9)*sps/40;
green = zeros(greenlen,nsta); % no bandpass
greenf = zeros(greenlen,nsta);  % bandpassed version
ppeaks = zeros(nsta,1); % positive peaks
npeaks = zeros(nsta,1); % negative peaks
greenort = zeros(greenlen,nsta); % no bandpass
greenfort = zeros(greenlen,nsta);  % bandpassed version
for ista = 1: nsta
  %cut according to the zero-crossing and the time shift from the constrained CC
  green(:,ista) = tmpwlet(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenf(:,ista) = tmpwletf(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  %detrend again for caution
  green(:,ista) = detrend(green(:,ista));
  greenf(:,ista) = detrend(greenf(:,ista));
  if normflag
    %normalize by max amp
    green(:,ista) = green(:,ista)/max(abs(green(:,ista)));    % normalize
    greenf(:,ista) = greenf(:,ista)/max(abs(green(:,ista)));    % normalize
  end
  
  %same process for orthogonal
  greenort(:,ista) = tmpwletort(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenfort(:,ista) = tmpwletfort(zcsta1+8*sps-greenlen+1-offwlet1i(ista): zcsta1+8*sps-offwlet1i(ista), ista);
  greenort(:,ista) = detrend(greenort(:,ista));
  greenfort(:,ista) = detrend(greenfort(:,ista));
  if normflag
    greenort(:,ista) = greenort(:,ista)/max(abs(green(:,ista)));    % normalize
    greenfort(:,ista) = greenfort(:,ista)/max(abs(green(:,ista)));    % normalize
  end
  
  %re-find the zero-crossing as the template length has changed
  [~,imin] = min(greenf(:,ista));
  [~,imax] = max(greenf(:,ista));
  [~,zcrosses(ista)] = min(abs(greenf(imin:imax,ista)));
  zcrosses(ista) = zcrosses(ista)+imin-1;
  ppeaks(ista) = imax;
  npeaks(ista) = imin;
end
%the following is just a check, because now the templates must be best aligned 
ccmid = round(size(greenf,1)/2);
ccwlen = 4*sps;
loffmax = 5*sps/40;
ccmin = 0.01;  % depending on the length of trace, cc could be very low
iup = 1;    % times of upsampling
[off12con,off13con,cc] = constrained_cc_interp(greenf(:,1:3)',ccmid,...
  ccwlen,mshiftadd,loffmax,ccmin,iup);
if ~(off12con==0 && off13con==0)
  disp('Filtered templates are NOT best aligned \n');
end
for ista = 4: nsta
  [coef,lag] = xcorr(greenf(:,1), greenf(:,ista), mshiftadd, 'coeff');
  [mcoef, idx] = max(coef);   % max of master raw cc
  if lag(idx)~=0   % offset in samples
    fprintf('Filtered templates are NOT best aligned at %s \n',stas(ista,:));
  end
end

amprat(1,:) = minmax(greenf(:,1)')./minmax(greenf(:,2)');	% amp ratio between max at sta 3 and 2 or min
amprat(2,:) = minmax(greenf(:,1)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min  
amprat(3,:) = minmax(greenf(:,2)')./minmax(greenf(:,3)');	% amp ratio between max at sta 3 and 1 or min  
spread = range(greenf);   % range of the amp of template

%%%plot the unfiltered and filtered templates
% plt_templates(green,greenf,stas,greenort,greenfort,lowlet,hiwlet,sps);

%just the filtered templates
% plt_templates_bp(greenf,stas,lowlet,hiwlet,sps);

sig = greenf(:,3);
lsig = length(sig);
times = (1:lsig)/sps;

% %%
% %all using the same unit, sec
% timeg = 1/sps: 1/sps: 6; 
% mu = 3;
% sigma = 1;
% gau = gauss(timeg,mu,sigma);
% plot(timeg,gau);%timeg,
% xlabel('Time (s)');
% ylabel('Amplitude');
% text(0.9,0.9,sprintf('\\mu=%.2f, \\sigma=%.2f',mu,sigma),'Units','normalized','HorizontalAlignment','right');

%% define a time-domain Gaussian function
% close all

timeg = 1/sps: 1/sps: 5;
mu = 2*sps;
% sigma = [0.01, 1, 2, 4, 8, 16, 32];
sigma = [8, 16, 32, 64, 128];
nsigma = length(sigma);
figure;
ncol = 3;
for i = 1: nsigma
% i=3;
  isigma = sigma(i);
  gau = gauss(timeg,mu/sps,isigma/sps);
%   denorm = max(gau);  %normalize by max
  denorm = trapz(gau);  %normalize to an area of 1
%   denorm = 1/isigma;  %
%   denorm = 1;
  gau = gau'/ denorm;
  subplot(nsigma,ncol,(i-1)*ncol+1);
  plot(timeg,gau); hold on
  ax=gca;
  plot([mu+1*isigma mu+1*isigma]/sps,ax.YLim,'k--','linew',1); 
  plot([mu-1*isigma mu-1*isigma]/sps,ax.YLim,'k--','linew',1); 
  text(0.9,0.9,sprintf('(\\mu=%d, \\sigma=%d) /160',mu,isigma),'Units','normalized','HorizontalAlignment','right');
  if i== nsigma
    xlabel('Time (s)');
    ylabel('Amplitude');    
  end
  
  filter = gau;
  [~,itfilter] = max(gau);
  
  [pf,ff] = periodogram(filter,hann(length(filter)),pow2(nextpow2(length(filter))),sps);
  [ps,fs] = periodogram(sig,hann(length(sig)),pow2(nextpow2(length(sig))),sps);
  subplot(nsigma,ncol,(i-1)*ncol+2);
  loglog(ff,pf,'r'); hold on
  loglog(fs,ps,'k');
%   axis equal
  xlim([0.1 sps/2]);
  ylim([1e-40 1]);
  ax=gca;
  loglog([1/(1*isigma/sps) 1/(1*isigma/sps)],ax.YLim,'k--','linew',1);  % a reference line with a slope of -1
  loglog([1/(2*isigma/sps) 1/(2*isigma/sps)],ax.YLim,'k-.','linew',1);  % a reference line with a slope of -1
  loglog([1/(3*isigma/sps) 1/(3*isigma/sps)],ax.YLim,'k:','linew',1);  % a reference line with a slope of -1
%   text(0.9,0.9,sprintf('\\mu=%.2f, \\sigma=%.4f',mu,isigma),'Units','normalized','HorizontalAlignment','right');
  if i== nsigma
    xlabel('Frequency (Hz)');
    ylabel('PSD');    
  end
  
%   loglog(ax.XLim,[sqrt(min(pxx)*max(pxx)) sqrt(min(pxx)*max(pxx))],'r--','linew',0.5);
  
  
  sconvgfull = conv(sig, filter, 'full');
  sconvg = sconvgfull(itfilter:lsig+itfilter-1);  % cut accordingly
  sconvgs = sconvg*range(sig)/range(sconvg);

  subplot(nsigma,ncol,i*ncol);
  plot(times,sig,'k'); hold on %
  plot((1:length(sconvgfull))/sps, sconvgfull+0.2,'b');
  plot(times,sconvgs,'r');
  ylim([-0.2 0.4]);
  if i== nsigma
    xlabel('Time (s)');
    ylabel('Amplitude');    
  end  
end

keyboard

%%
timeg = 1/sps: 1/sps: 5;
mu = 2.5*sps;
isigma = 32;
gau = gauss(timeg,mu/sps,isigma/sps);
%   denorm = max(gau);  %normalize by max
denorm = trapz(gau);  %normalize to an area of 1
%   denorm = 1/isigma;  %
%   denorm = 1;
gau = gau'/ denorm;

[~,itfilter] = max(gau);

nfft = pow2(nextpow2(max(length(gau),length(sig))));
[gauF,f,amp1,pha1,power,psd] = fftspectrum(gau, nfft, sps, 'twosided');
[sigF,f,amp2,pha2,power,psd] = fftspectrum(sig, nfft, sps, 'twosided');
sconvgF = gauF.*sigF;
sconvgfull2 = real(ifft(sconvgF,nfft));
sconvgF2 = (amp1.*nfft.*exp(1i.*pha1)).*(amp2.*nfft.*exp(1i.*(pha2)));
sconvgfull22 = real(ifft(sconvgF2,nfft));

% sconvgfull = sconvgfull*range(sig)/range(sconvgfull);
figure
plot(times,sig,'k'); hold on %
plot((1:length(sconvgfull))/sps, sconvgfull*50+0.2,'r');
plot((1:length(sconvgfull2))/sps, sconvgfull2*50+0.2,'b');
plot((1:length(sconvgfull22))/sps, sconvgfull22*50+0.2,'g');
ax=gca;
plot([itfilter itfilter]/sps,ax.YLim,'b--');
plot([lsig+itfilter-1 lsig+itfilter-1]/sps,ax.YLim,'b--');
% plot(times,sconvgs,'r');

keyboard

%%
figure
plot(times,sig,'k'); hold on %
plot(times,sconvgfull2*10,'r');


%% define a heaviside step function
close all

timeh = -2.5: 1/sps: 2.5;
hea = heaviside(timeh);
denorm = 1;
% denorm = sum(hea);  %normalize to an area of 1
hea = hea/ denorm;

% plot(hea);

filter = hea;
itfilter = find(hea==0.5/denorm);
sconvh = conv(sig, filter, 'full');
sconvh = sconvh(itfilter:lsig+itfilter-1);  % cut accordingly

figure;
plot(times,5*sig,'k'); hold on
plot(times,sconvh,'r');





