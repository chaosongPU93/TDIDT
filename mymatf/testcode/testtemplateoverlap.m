% testtemplateoverlap.m
%
% this script is to test the stacking of 2 templates due to the different overlapping in time, or in
% other words, how close their zero-crossing are. Hopefully we can get some feeling about how close
% the 2 templates are can it result in the construction or destruction of the peaks/troughs.
%
% The way to do it is to set one template not moving, and move the zero-crossing of the other by
% different amounts in samples relative to that of the unmoving template. or check the resulting
% statistics.

%% Initialization
clear
clc
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');

[scrsz, res] = pixelperinch(1);

%% prepare templates, inherited from 'synthshift_chao.m'
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

for ista=1:nsta
    STA(:,ista)=Bandpass(STA(:,ista),sps,0.1,15,2,2,'butter');   % change 'bandpass' to 'Bandpass'
end

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
  greenf(:,ista)=greenf(:,ista)/max(abs(greenf(:,ista)));    % normalize
  
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

%% choose one station, keep one template not moving, move the other by some amount
ista = 1;
temp = greenf(:,1);
tgreen = zcrosses(ista);

inddelay = (-50:1:50)';
ntry = length(inddelay);
rst = zeros(2*greenlen,ntry);
indst = round(greenlen/2)+1;
inded = indst+greenlen-1;
rst(indst:inded, :) = repmat(temp,1,ntry);

Greenlenm1 = greenlen-1;

% figure
% plot(temp,'k');

figure
hold on
grid on; box on
for i = 1:5:ntry
  plot(10*rst(:,i)+(i-1)*2,'k');
end
ax=gca;
plot([tgreen+indst-1 tgreen+indst-1],ax.YLim,'--','linew',2);
axis([1500 2500 -10 210]);
xlabel('Samples at 160 sps');

indzc = zeros(ntry,1); 
for i = 1:ntry
  ind = indst+inddelay(i);
  if ind <= 0   %if arrival index is smaller than 1
    trunc = -ind;
    rst(1: Greenlenm1-trunc, i) = rst(1: Greenlenm1-trunc, i) + ...
      temp(2+trunc: end);
  elseif ind+Greenlenm1 > size(rst,1) %if arrival index is larger than length
    trunc = ind+Greenlenm1 - size(rst,1);
    rst(ind: ind+Greenlenm1-trunc, i) = ...
      rst(ind: ind+Greenlenm1-trunc, i) + ...
      temp(end-trunc);
  else
    rst(ind: ind+Greenlenm1, i) = rst(ind: ind+Greenlenm1, i) + temp; %Greens(whichG,:) %(for interpolated)

  end
  indzc(i) = ind+tgreen-1;
  
end

figure
hold on
grid on; box on
for i = 1:5:ntry
  plot(10*rst(:,i)+(i-1)*2,'k');
  text(2400,(i-1)*2-1,num2str(inddelay(i)),'Color','b','FontSize',12);
  scatter(indzc(i),10*rst(indzc(i),i)+(i-1)*2, 20,'r');
end
plot([tgreen+indst-1 tgreen+indst-1],ax.YLim,'--','linew',2);
axis([1500 2500 -10 210]);
xlabel('Samples at 160 sps');

figure
hold on
grid on; box on
for i = 31:1:36
  plot(10*rst(:,i)+(i-31)*3,'k');
  text(2400,(i-31)*3-1,num2str(inddelay(i)),'Color','b','FontSize',12);
  scatter(indzc(i),10*rst(indzc(i),i)+(i-31)*3, 20,'r');
end
for i = 66:1:71
  plot(10*rst(:,i)+(i-58)*3,'k');
  text(2400,(i-58)*3-1,num2str(inddelay(i)),'Color','b','FontSize',12);
  scatter(indzc(i),10*rst(indzc(i),i)+(i-58)*3, 20,'r');
end
plot([tgreen+indst-1 tgreen+indst-1],ax.YLim,'--','linew',2);
axis([1500 2500 -10 50]);
xlabel('Samples at 160 sps');


%% feed the stacked seimogram to 'iterdecon' see how can deconvolution separate them
for i = 1:1:ntry
  wlet = temp;
  lwlet = length(wlet);
  sig = rst(:,i);
  noi = [];
  
  %detrend and taper
  sig = detrend(sig);
  fractap = 0.05; % if fractap is >=1, n-point von Hann window is returned
  ptstap = fractap/2*size(sig,1); % if fractap is >=1, n-point von Hann window is returned
  w = tukeywin(size(sig,1),fractap);
  % sig = w.* sig;
  %detrend again for caution
  sig = detrend(sig);
  lsig = length(sig);
  
  dt = 1/sps;  % sampling interval
  twlet = tgreen*dt;
  width = 2.5;  % width for Gaussian filter
  dres_min = 0.5;  % tolerance, percentage change in residual per iteration
  mfit_min = 5e-1;  % tolerance, norm of the residual/misfit
  nit_max = 4;  % max numer of iterations
  npul_max = 4;%a single peak is ~20 samples wide; maybe a little less (at 100 sps). ~0.4s, 1/0.4=2.5
  fpltit = 0;  % plot flag for each iteration
  fpltend = 0;  % plot flag for the final iteration
  fcheck = 0; % plot flag for intermediate computations
  % rcc = [];  % running CC between diff. stations
  rcc = ones(length(sig) ,1); % for testing purpose
  [sigdecon(:,i),pred,res,dresit,mfitit,ampit,nit,fighdl] = ...
    iterdecon(sig,wlet,rcc,noi,dt,twlet,width,dres_min,mfit_min,nit_max,npul_max,fpltit,fpltend,fcheck);
  
  nit
  imp(:,i) = ampit(:,1);
  impamp(:,i) = ampit(:,2);
end

%%
figure
hold on
grid on; box on
plot(inddelay,inddelay,'r-','linew',1);
plot(inddelay,-25:0.5:25,'b-','linew',1);
for i = 1:1:ntry 
%   [~,ind] = sort(impamp(:,i),'descend');
%   hi2(:,1) = imp(ind(1:2), i);
%   hi2(:,2) = impamp(ind(1:2), i);
%   [~,ind1] = min(abs(hi2(:,1)-(tgreen+indst-1)));
%   ind2 = setdiff(1:2, ind1);
%   delay(i,1) = hi2(ind2,1)-hi2(ind1,1);
%   scatter(inddelay(i),delay(i,1),20,'ko');

  scatter(inddelay(i)*ones(4,1),imp(:,i)-(tgreen+indst-1),impamp(:,i)*15,'ko','filled');

end
axis equal
axis([-60 60 -60 60]);
xlabel('Delay from ground truth (samples at 160 sps)');
ylabel('Delay from deconvolved impulses (samples at 160 sps)');

% print('-dpdf','-r300',strcat('/home/chaosong/Pictures/IndepIter',num2str(nit),'.jpg'));






