%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Similar to psd_lfe_tremor.m, this is also to compare the psd estimate of lfe
% stacks and tremor records, using proper methods determined from trial and
% testing.
%
% FEATURES   
%   1. Compare the different lfe stacks: broadband direct stack; broadband CC 
%       stack; bandpassed (0.1-15 hz) direct stack; bandpassed (0.1-15 hz) CC
%       stack
%   2. Include noise psd as well, take the real psd as the (Signal_psd-
%       Noise_psd), with a proper proxy for representing noise; to both lfe and
%       tremor;
%   3. For psd estimate function: choose mtm/periodgram (no averaging) for lfe;
%       choose pchave for tremor (averaging, taper and robust algorithm written
%       by Frederik J. Simons)
%
%
% NOTES:
%   1. To simulate the noise in lfe stacks, we should cut the each 20-s data 
%       shifted by 10s before each lfe timing, get the spectrum of them, then
%       get the median of all spectrum as the estimate for noise.
%   2. It is possible that for some sampled frequencies, noise amplitude is 
%       higher than signal, so that subtraction would lead to negative value.
%       When converted to dB, DO NOT plot them.
%   3. To simualte the nosie for tremor, we should choose tremor free days, 
%       several days, several periods the same as migrations. Also, we need to
%       get the spectrum of each 20-s window, then obtain the median/mean, or 
%       say, 90 percentile as the error estimate. Finally, do NOT plot the 
%       negative value when subtract the nosie from the signal
%       
% NOTICE:
%   set(gcf,'renderer','Painters')   to force the figure to be vector image
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/11/25
% Last modified date:   2019/12/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
% scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]
scrsz=get(0,'MonitorPositions'); % return the position of all monitors

% set path
workpath = getenv('ALLAN');
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
temppath1 = strcat(datapath, '/templates/PGCtrio/');
temppath2 = strcat(datapath, '/templates/LZBtrio/');

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
    'LZB'];
POLSTA=['SSIB '           % polaris station names
    'SILB '
    'KLNB '
    'MGCB '
    'TWKB '];
trioflag = 1;

if trioflag == 1    % use PGC trio
    fam = '002';

    % get rots
    FLAG = 'PGC';
    CATA = 'fixed';
    [PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,0,0);
    
    stas=['PGC  '
          'SSIB '
          'SILB '];
      
    temppath = temppath1;
      
elseif trioflag == 2    % use LZB trio
    
    nfampool = ['002';
    '043';
    '141';
    '047';
    '010';
    '144';
    '099';
    '068';
    '125';
    '147';
    '017'];

    ifam = 1;
    fam = nfampool(ifam, :);

    % get permanent and polaris station rotation parameters
    sft2=0;     % centroid shift of station 2
    sft3=0;     % centroid shift of station 3
    FLAG = 'LZB';
    CATA = 'NEW';
    [PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,sft2,sft3);
    reftime = POLROTS(5,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
    PERMROTS(:,4) = PERMROTS(:,4)-reftime;    % to make sure that 1st station is 0
    POLROTS(:,4) = POLROTS(:,4)-reftime;
    
    stas=['TWKB '
          'LZB  '
          'MGCB '];
      
    temppath = temppath2;  
      
end

% number of used stations
nsta=size(stas,1);         %  number of stations

% convert angles to rads
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;     % convert 2-3 columns to rad
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;


%% load lfe stacks
templensec = 60;
sps = 40;
templen = templensec * sps;
lfe = zeros(nsta, templen);
for ista = 1:nsta
    
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
        num2str(sps), 'sps_', num2str(templensec), 's_', ...
        'BBDS_', 'opt_Nof_Non_Chao');        
    bbdslfe(ista,:) = load(fname);
    
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
        num2str(sps), 'sps_', num2str(templensec), 's_', ...
        'BBCCS_', 'opt_Nof_Non_Chao');
    bbcclfe(ista,:) = load(fname);
    
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
        num2str(sps), 'sps_', num2str(templensec), 's_', ...
        'BPDS_', 'opt_Nof_Non_Chao');
    bpdslfe(ista,:) = load(fname);
    
    strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
        num2str(sps), 'sps_', num2str(templensec), 's_', ...
        'BPCCS_', 'opt_Nof_Non_Chao');
    bpcclfe(ista,:) = load(fname);
    
end

f.fig = figure;
set(f.fig,'Position',[scrsz(2,1)*0.9 0 scrsz(2,3)*0.6 scrsz(2,4)*0.5]);
subplot(4,1,1)
for ista = 1:nsta
    plot((1:templen)/sps,bbdslfe(ista,:),'linew',1); hold on
end
title('broadband direct lfe stacks')
legend({'PGC','SSIB','SILB'},'location','southeast');

subplot(4,1,2)
for ista = 1:nsta
    plot((1:templen)/sps,bbcclfe(ista,:),'linew',1); hold on
end
title('broadband cc lfe stacks')
legend({'PGC','SSIB','SILB'},'location','southeast');

subplot(4,1,3)
for ista = 1:nsta
    plot((1:templen)/sps,bpdslfe(ista,:),'linew',1); hold on
end
title('bandpassed (0.1-15 hz) direct lfe stacks')
legend({'PGC','SSIB','SILB'},'location','southeast');

subplot(4,1,4)
for ista = 1:nsta
    plot((1:templen)/sps,bpcclfe(ista,:),'linew',1); hold on
end
title('bandpassed (0.1-15 hz) cc lfe stacks')
legend({'PGC','SSIB','SILB'},'location','southeast');


%% load one day of tremor records, take 2005 255 as an example
timoffrot= [2005 255];
bostname=['BOSTOCK/NEW/002-246_2005.255'];
% read one day of day using rotation para without filtering, broadband
[tremor,timeperm] = rdnofiltPGC(fam,timoffrot,datapath,stas,sps,PERMSTA,PERMROTS,POLSTA,...
                                POLROTS);
f.fig = figure;
set(f.fig,'Position',[scrsz(2,1)*0.9 0 scrsz(2,3)*0.7 scrsz(2,4)*0.3]);
for ista = 1:nsta
    plot(timeperm,tremor(ista,:)+(ista-2)*5,'linew',1); hold on
end
ylim([-15 15]);
title('tremor records')
legend({'PGC','SSIB','SILB'},'location','southeast');

trange = [
          2005255,3.42e+4,3.65e+4;
          2005255,5.80e+4,5.96e+4;          
          2005255,6.70e+4,6.90e+4];

tremoruse=[];
for ista=1:nsta
    tmp1=[];
    for j = 1: length(trange)
         tmp2 = tremor(ista,timeperm>=trange(j,2)&timeperm<=trange(j,3));
         tmp1 = [tmp1 tmp2];
    end
    tremoruse(ista,:) = tmp1;
end


%% load one tremor-free day, to simulate the noise
timoffrot= [2005 244];
% read one day of day using rotation para without filtering, broadband
[tremor,timeperm] = rdnofiltPGC(fam,timoffrot,datapath,stas,sps,PERMSTA,PERMROTS,POLSTA,...
                                POLROTS);
f.fig = figure;
set(f.fig,'Position',[scrsz(2,1)*0.9 0 scrsz(2,3)*0.7 scrsz(2,4)*0.3]);
for ista = 1:nsta
    plot(timeperm,tremor(ista,:)+(ista-2)*5,'linew',1); hold on
end
ylim([-15 15]);
title('tremor noise')
legend({'PGC','SSIB','SILB'},'location','southeast');

trange2 = [
          2005244,3.42e+4,3.65e+4;
          2005244,5.80e+4,5.96e+4;          
          2005244,6.70e+4,6.90e+4];

tremornoi=[];
for ista=1:nsta
    tmp1=[];
    for j = 1: length(trange2)
         tmp2 = tremor(ista,timeperm>=trange2(j,2)&timeperm<=trange2(j,3));
         tmp1 = [tmp1 tmp2];
    end
    tremornoi(ista,:) = tmp1;
end

%% use periodgram to get psd estimate for lfe stacks
f.fig = figure;
set(f.fig,'Position',[scrsz(2,1)*0.9 0 scrsz(2,3)*0.6 scrsz(2,4)*0.6]);
subplot(4,1,1)
cpool1 = [0.5 0.5 1; 1 0.5 0.5; 0.5 1 0.5];
cpool2 = [0 0 1;1 0 0; 0 1 0];
for ista = 1:nsta
    x = bbdslfe(ista,:);
    n = length(x);
    sig = x(n/3+1:n*2/3);
    noi = x(1:n/3);
    nfft = pow2(nextpow2(length(sig)-1));   % no zero-padding
    window = hamming(length(sig));
    Fs = sps;
    [psig,ft] = periodogram(sig,window,nfft,Fs);   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(psig),'--','color',cpool1(ista,:),'linew',1); hold on
    
    [pnoi,ft] = periodogram(noi,window,nfft,Fs);
    semilogx(ft, 10*log10(pnoi),':','color',[0.5 0.5 0.5],'linew',1);
    
    p(ista) = semilogx(ft, 10*log10(psig-pnoi),'o-','markers',1.5,'color',cpool2(ista,:),'linew',1);
    
end
xlim([1e-2,Fs]);
ylim([-150,-20]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('BBDS lfe, periodgram, hamming, nfft=nextpow2-1')
legend([p(1) p(2) p(3)],{'PGC','SSIB','SILB'},'location','southwest');

subplot(4,1,2)
for ista = 1:nsta
    x = bbcclfe(ista,:);
    n = length(x);
    sig = x(n/3+1:n*2/3);
    noi = 0.5*(x(1:n/3) + x(n*2/3+1:n));
    nfft = pow2(nextpow2(length(sig)-1));   % no zero-padding
    window = hamming(length(sig));
    Fs = sps;
    [psig,ft] = periodogram(sig,window,nfft,Fs);   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(psig),'--','color',cpool1(ista,:),'linew',1); hold on
    
    [pnoi,ft] = periodogram(noi,window,nfft,Fs);
    semilogx(ft, 10*log10(pnoi),':','color',[0.5 0.5 0.5],'linew',1);
    
    p(ista) = semilogx(ft, 10*log10(psig-pnoi),'o-','markers',1.5,'color',cpool2(ista,:),'linew',1);
    
end
xlim([1e-2,Fs]);
ylim([-150,-20]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('BBCC lfe, periodgram, hamming, nfft=nextpow2-1')
legend([p(1) p(2) p(3)],{'PGC','SSIB','SILB'},'location','southwest');

subplot(4,1,3)
for ista = 1:nsta
    x = bpdslfe(ista,:);
    n = length(x);
    sig = x(n/3+1:n*2/3);
    noi = 0.5*(x(1:n/3) + x(n*2/3+1:n));
    nfft = pow2(nextpow2(length(sig)-1));   % no zero-padding
    window = hamming(length(sig));
    Fs = sps;
    [psig,ft] = periodogram(sig,window,nfft,Fs);   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(psig),'--','color',cpool1(ista,:),'linew',1); hold on
    
    [pnoi,ft] = periodogram(noi,window,nfft,Fs);
    semilogx(ft, 10*log10(pnoi),':','color',[0.5 0.5 0.5],'linew',1);
    
    p(ista) = semilogx(ft, 10*log10(psig-pnoi),'o-','markers',1.5,'color',cpool2(ista,:),'linew',1);
    
end
xlim([1e-2,Fs]);
ylim([-150,-20]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('BPDS lfe, periodgram, hamming, nfft=nextpow2-1')
legend([p(1) p(2) p(3)],{'PGC','SSIB','SILB'},'location','southwest');

subplot(4,1,4)
for ista = 1:nsta
    x = bpcclfe(ista,:);
    n = length(x);
    sig = x(n/3+1:n*2/3);
    noi = 0.5*(x(1:n/3) + x(n*2/3+1:n));
    nfft = pow2(nextpow2(length(sig)-1));   % no zero-padding
    window = hamming(length(sig));
    Fs = sps;
    [psig,ft] = periodogram(sig,window,nfft,Fs);   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(psig),'--','color',cpool1(ista,:),'linew',1); hold on
    
    [pnoi,ft] = periodogram(noi,window,nfft,Fs);
    semilogx(ft, 10*log10(pnoi),':','color',[0.5 0.5 0.5],'linew',1);
    
    p(ista) = semilogx(ft, 10*log10(psig-pnoi),'o-','markers',1.5,'color',cpool2(ista,:),'linew',1);
    
end
xlim([1e-2,Fs]);
ylim([-150,-20]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('BPCC lfe, periodgram, hamming, nfft=nextpow2-1')
legend([p(1) p(2) p(3)],{'PGC','SSIB','SILB'},'location','southwest');



%% use pmtm to get psd estimate for lfe stacks
f.fig = figure;
set(f.fig,'Position',[scrsz(2,1)*0.9 0 scrsz(2,3)*0.6 scrsz(2,4)*0.6]);
subplot(4,1,1)
cpool1 = [0.5 0.5 1; 1 0.5 0.5; 0.5 1 0.5];
cpool2 = [0 0 1;1 0 0; 0 1 0];
for ista = 1:nsta
    x = bbdslfe(ista,:);
    n = length(x);
    sig = x(n/3+1:n*2/3);
    noi = 0.5*(x(1:n/3) + x(n*2/3+1:n));
    nfft = pow2(nextpow2(length(sig)-1));   % no zero-padding
    nw = 4;  % For multi-taper spectral analysis, typical choices for NW are 2, 5/2, 3, 7/2, or 4.
    Fs = sps;
    [psig,ft] = pmtm(sig,nw,nfft,Fs);   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(psig),'--','color',cpool1(ista,:),'linew',1); hold on
    
    [pnoi,ft] = pmtm(noi,nw,nfft,Fs);
    semilogx(ft, 10*log10(pnoi),':','color',[0.5 0.5 0.5],'linew',1);
    
    p(ista) = semilogx(ft, 10*log10(psig-pnoi),'o-','markers',1.5,'color',cpool2(ista,:),'linew',1);
    
end
xlim([1e-2,Fs]);
ylim([-150,-20]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('BBDS lfe, pmtm, nw=4, nfft=nextpow2-1')
legend([p(1) p(2) p(3)],{'PGC','SSIB','SILB'},'location','southwest');

subplot(4,1,2)
for ista = 1:nsta
    x = bbcclfe(ista,:);
    n = length(x);
    sig = x(n/3+1:n*2/3);
    noi = 0.5*(x(1:n/3) + x(n*2/3+1:n));
    nfft = pow2(nextpow2(length(sig)-1));   % no zero-padding
    nw = 4;  % For multi-taper spectral analysis, typical choices for NW are 2, 5/2, 3, 7/2, or 4.
    Fs = sps;
    [psig,ft] = pmtm(sig,nw,nfft,Fs);   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(psig),'--','color',cpool1(ista,:),'linew',1); hold on
    
    [pnoi,ft] = pmtm(noi,nw,nfft,Fs);
    semilogx(ft, 10*log10(pnoi),':','color',[0.5 0.5 0.5],'linew',1);
    
    p(ista) = semilogx(ft, 10*log10(psig-pnoi),'o-','markers',1.5,'color',cpool2(ista,:),'linew',1);
    
end
xlim([1e-2,Fs]);
ylim([-150,-20]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('BBCC lfe, pmtm, nw=4, nfft=nextpow2-1')
legend([p(1) p(2) p(3)],{'PGC','SSIB','SILB'},'location','southwest');

subplot(4,1,3)
for ista = 1:nsta
    x = bpdslfe(ista,:);
    n = length(x);
    sig = x(n/3+1:n*2/3);
    noi = 0.5*(x(1:n/3) + x(n*2/3+1:n));
    nfft = pow2(nextpow2(length(sig)-1));   % no zero-padding
    nw = 4;  % For multi-taper spectral analysis, typical choices for NW are 2, 5/2, 3, 7/2, or 4.
    Fs = sps;
    [psig,ft] = pmtm(sig,nw,nfft,Fs);   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(psig),'--','color',cpool1(ista,:),'linew',1); hold on
    
    [pnoi,ft] = pmtm(noi,nw,nfft,Fs);
    semilogx(ft, 10*log10(pnoi),':','color',[0.5 0.5 0.5],'linew',1);
    
    p(ista) = semilogx(ft, 10*log10(psig-pnoi),'o-','markers',1.5,'color',cpool2(ista,:),'linew',1);
    
end
xlim([1e-2,Fs]);
ylim([-150,-20]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('BPDS lfe, pmtm, nw=4, nfft=nextpow2-1')
legend([p(1) p(2) p(3)],{'PGC','SSIB','SILB'},'location','southwest');

subplot(4,1,4)
for ista = 1:nsta
    x = bpcclfe(ista,:);
    n = length(x);
    sig = x(n/3+1:n*2/3);
    noi = 0.5*(x(1:n/3) + x(n*2/3+1:n));
    nfft = pow2(nextpow2(length(sig)-1));   % no zero-padding
    nw = 4;  % For multi-taper spectral analysis, typical choices for NW are 2, 5/2, 3, 7/2, or 4.
    Fs = sps;
    [psig,ft] = pmtm(sig,nw,nfft,Fs);   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(psig),'--','color',cpool1(ista,:),'linew',1); hold on
    
    [pnoi,ft] = pmtm(noi,nw,nfft,Fs);
    semilogx(ft, 10*log10(pnoi),':','color',[0.5 0.5 0.5],'linew',1);
    
    p(ista) = semilogx(ft, 10*log10(psig-pnoi),'o-','markers',1.5,'color',cpool2(ista,:),'linew',1);
    
end
xlim([1e-2,Fs]);
ylim([-150,-20]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('BPCC lfe, pmtm, nw=4, nfft=nextpow2-1')
legend([p(1) p(2) p(3)],{'PGC','SSIB','SILB'},'location','southwest');



%% use pchave to get psd estimate for tremor detections
% break continous records into segments
wlen = 20*sps;
nwin = floor(size(tremoruse,2)/wlen);
tresig1=zeros(nwin,wlen);
tresig2=zeros(nwin,wlen);
tresig3=zeros(nwin,wlen);
for i=1:nwin
    tresig1(i,:) = tremoruse(1,(i-1)*wlen+1:i*wlen);
    tresig2(i,:) = tremoruse(2,(i-1)*wlen+1:i*wlen);
    tresig3(i,:) = tremoruse(3,(i-1)*wlen+1:i*wlen);
end

trenoi1=zeros(nwin,wlen);
trenoi2=zeros(nwin,wlen);
trenoi3=zeros(nwin,wlen);
for i=1:nwin
    trenoi1(i,:) = tremornoi(1,(i-1)*wlen+1:i*wlen);
    trenoi2(i,:) = tremornoi(2,(i-1)*wlen+1:i*wlen);
    trenoi3(i,:) = tremornoi(3,(i-1)*wlen+1:i*wlen);
end
trenoi(1,:) = sum(trenoi1,1)/nwin;
trenoi(2,:) = sum(trenoi2,1)/nwin;
trenoi(3,:) = sum(trenoi3,1)/nwin;


pxxmat1=[];
pxxmat2=[];
pxxmat3=[];


f.fig = figure;
set(f.fig,'Position',[scrsz(2,1)*0.9 0 scrsz(2,3)*0.6 scrsz(2,4)*0.45]);
f.ax(1)=subplot(3,1,1);
f.ax(2)=subplot(3,1,2);
f.ax(3)=subplot(3,1,3);

cpool1 = [0.7 0.7 1; 1 0.7 0.7; 0.7 1 0.7];
cpool2 = [0.3 0.3 1; 1 0.3 0.3; 0.3 1 0.3];
cpool3 = [0 0 1;1 0 0; 0 1 0];

nfft = 256;
lwin = 256;
olap = 50;
Fs = sps;
pnoi = zeros(nfft/2+1,nsta);
for i = 1: nsta
    x=trenoi(i,:);
    [pnoi(:,i),ft]=pchave(x,lwin,olap,nfft,Fs,'MAD','dpss');
    semilogx(f.ax(i),ft,10*log10(pnoi(:,i)),':','color',[0.5 0.5 0.5]); hold(f.ax(i),'on');
    xlim(f.ax(i),[5e-2,50]);
    ylim(f.ax(i),[-110,0]);
    xlabel(f.ax(i),'Frequency (Hz)');
    ylabel(f.ax(i),'PSD (dB/Hz)');
end

for i=1:nwin
    x=tresig1(i,:);
    nfft = 256;
    lwin = 256;
    olap = 50;    
    Fs = sps;
    [pxx,ft]=pchave(x,lwin,olap,nfft,Fs,'MAD','dpss');
    semilogx(f.ax(1),ft,10*log10(pxx),'--','color',cpool1(1,:),'linew',2); hold(f.ax(1),'on');
%     p1=semilogx(f.ax(1),ft, 10*log10(pxx),'color',[0.5 0.5 1]); hold(f.ax(1),'on');
%     p1.Color(4) = 0.5;
    x = ft';
    y = real(10*log10(pxx-pnoi(:,1))');
    % Now make the 'line' (actually a surface)...
    z = zeros(size(x));
    surface(f.ax(1),[x;x],[y;y],[z;z],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',2,...
            'edgealpha',.2,...
            'edgecolor',cpool2(1,:)); hold(f.ax(1),'on');
    set(f.ax(1),'Xscale','log');
    
    pxxmat1(i,:) = (pxx-pnoi(:,1))';
           
    x=tresig2(i,:);
    [pxx,ft]=pchave(x,lwin,olap,nfft,Fs,'MAD','dpss');
    semilogx(f.ax(2),ft,10*log10(pxx),'--','color',cpool1(2,:),'linew',1); hold(f.ax(2),'on');
%     p2=semilogx(ft, 10*log10(pxx),'color','b');
%     p2.Color(4) = 0.5;
    x = ft';
    y = real(10*log10(pxx-pnoi(:,2))');
    % Now make the 'line' (actually a surface)...
    z = zeros(size(x));
    surface(f.ax(2),[x;x],[y;y],[z;z],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',2,...
            'edgealpha',.2,...
            'edgecolor',cpool2(2,:)); hold(f.ax(2),'on');
    set(f.ax(2),'Xscale','log');

    pxxmat2(i,:) = (pxx-pnoi(:,2))';
    
    x=tresig3(i,:);
    [pxx,ft]=pchave(x,lwin,olap,nfft,Fs,'MAD','dpss');
    semilogx(f.ax(3),ft,10*log10(pxx),'--','color',cpool1(3,:),'linew',1); hold(f.ax(3),'on');
%     p3=semilogx(ft, 10*log10(pxx),'color','r');
%     p3.Color(4) = 0.5;
    x = ft';
    y = real(10*log10(pxx-pnoi(:,3))');
    % Now make the 'line' (actually a surface)...
    z = zeros(size(x));
    surface(f.ax(3),[x;x],[y;y],[z;z],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',2,...
            'edgealpha',.2,...
            'edgecolor',cpool2(3,:)); hold(f.ax(3),'on');
    set(f.ax(3),'Xscale','log');

    pxxmat3(i,:) = (pxx-pnoi(:,3))';
    
end

pxxm1 = median(pxxmat1,1);
ln1=semilogx(f.ax(1),ft, real(10*log10(pxxm1)),'color',cpool3(1,:),'linew',1.5,'marker','o','markers',2);
pxxm2 = median(pxxmat2,1);
ln2=semilogx(f.ax(2),ft, real(10*log10(pxxm2)),'color',cpool3(2,:),'linew',1.5,'marker','o','markers',2);
pxxm3 = median(pxxmat3,1);
ln3=semilogx(f.ax(3),ft, real(10*log10(pxxm3)),'color',cpool3(3,:),'linew',1.5,'marker','o','markers',2);

title(f.ax(1),'20-s tremor records,pchave,dpss,nfft=256,nw=4,MAD')
legend(f.ax(1),[ln1,ln2,ln3],{'PGC','SSIB','SILB'},'location','south');
title(f.ax(2),strcat('tremor:',num2str(trange(1,1)),'  noise:',num2str(trange2(1,1))));




 


































