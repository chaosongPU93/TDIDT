%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to compare the psd estimate of lfe stacks and tremor 
% records, using different implementation methods in matlab
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/11/02
% Last modified date:   2019/11/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
% clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
% scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]
scrsz=get(0,'MonitorPositions'); % return the position of all monitors

% set path
workpath = getenv('ALLAN');
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');

fam = '002';
% get days, 3-sta trio name, bostock's template name, zero crossing index
% of templates

% get rots
FLAG = 'PGC';
CATA = 'fixed';
[PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,0,0);
% POLROTS(1:4)=0;

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
    'LZB'];
POLSTA=['SSIB '           % polaris station names
    'SILB '
    'KLNB '
    'MGCB '
    'TWKB '];

% stas=['PGC  '
%     'SSIB '
%     'SILB '
%     'LZB  '
%     'TWKB '
%     'MGCB '
%     'KLNB '];

stas=['PGC  '
      'SSIB '
      'SILB '];

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
        num2str(sps), 'sps_', num2str(templensec), 's_','DS_', 'opt_Nof_Non_Chao');        
    lfe(ista,:) = load(fname);
end
figure
for ista = 1:nsta
    plot((1:templen)/sps,lfe(ista,:),'linew',1); hold on
end
title('lfe stacks')
legend({'PGC','SSIB','SILB'},'location','south');


%% load one day of tremor records, take 2005 255 as an example
timoffrot= [2005 255];
bostname=['BOSTOCK/NEW/002-246_2005.255'];
% read one day of day using rotation para without filtering, broadband
[tremor,timeperm] = rdnofiltPGC(fam,timoffrot,bostname,datapath,stas,sps,PERMSTA,PERMROTS,POLSTA,...
                                POLROTS);
figure
for ista = 1:nsta
    plot(timeperm,tremor(ista,:)+(ista-2)*5,'linew',1); hold on
end
ylim([-20 20]);
title('tremor records')
legend({'PGC','SSIB','SILB'},'location','south');

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

%% load one day of tremor detections, take 2005 255 hf detection as an example
datestr = num2str(2005255);
yr = datestr(1:4);
day = datestr(5:end);
winlenhf = 4;
winlen=winlenhf*sps;      % length in smaples
hi=6.5;
lo=1.25;
loff = 2.1;
ccmin = 0.44;
npo = 2;
npa = 2;
mshift = 29;
nstanew = 4;
IDENTIF = strcat(yr,'.',day,'.',fam,'.loff',num2str(loff),'.ccmin',num2str(ccmin),'.nponpa', ...
                     num2str(npo),num2str(npa),'.ms',num2str(mshift));
fname = strcat(datapath,'/PGCtrio/MAPS/seistraceall_',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                    num2str(winlen/sps),'s',num2str(sps),'sps',num2str(nstanew), 'newsta');
tmp = load(fname);    % detected windows
tremordet = tmp(:,1:4);

nwin = 50;
for i=1:nwin
    tredet1(i,:)=tremordet((i-1)*winlen+1:i*winlen,2);
    tredet2(i,:)=tremordet((i-1)*winlen+1:i*winlen,3);
    tredet3(i,:)=tremordet((i-1)*winlen+1:i*winlen,4);
end

%% use different psd estimate methods to lfe, tremor&detections
% close all
%% 1 periodgram
%%%% 1.1 to lfe
figure
for ista = 1:nsta
    x = lfe(ista,:);
    nfft = length(x);   % means no zero-padding
    window = hamming(nfft);
    Fs = sps;
    [pxx,ft] = periodogram(x,window,nfft,Fs,'psd');   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(pxx)); hold on
end
xlim([1e-2,20]);
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('lfe, periodgram, hamming, nfft=n')
legend({'PGC','SSIB','SILB'},'location','south');


%%%% 1.2 to tremor
figure
for ista = 1:nsta
    x = tremor(ista,:);
    nfft = length(x);   % means no zero-padding
    window = hamming(nfft);
    Fs = sps;
    [pxx,ft] = periodogram(x,window,nfft,Fs,'psd');   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(pxx)); hold on
end
xlim([1e-2,20]);
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('tremor, periodgram, hamming, nfft=n')
legend({'PGC','SSIB','SILB'},'location','south');


%% 2 pwelch
%%%% 2.1 to lfe
figure
for ista = 1:nsta
    x = lfe(ista,:);
    nfft = length(x)/4;   % means no zero-padding
    window = hamming(nfft);
    noverlap = nfft/2;    
    Fs = sps;
    [pxx,ft]=pwelch(x,window,noverlap,nfft,Fs,'psd');   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(pxx)); hold on
end
xlim([1e-2,20]);
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('lfe, pwelch, hamming, nfft=n/4')
legend({'PGC','SSIB','SILB'},'location','south');
    
%%%% 2.1.1 the previous way to lfe, with default options
figure
for ista = 1:nsta
    x = lfe(ista,:); 
    Fs = sps;
    [pxx,ft]=pwelch(x,[],[],[],Fs,'psd');   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(pxx)); hold on
end
xlim([1e-2,20]);
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('lfe, pwelch, defaults')
legend({'PGC','SSIB','SILB'},'location','south');

%%%% 2.1.2 to lfe, shorter window, more averging
figure
for ista = 1:nsta
    x = lfe(ista,:);
    nfft = length(x)/6;   % means no zero-padding
    window = hamming(nfft);
    noverlap = nfft/2;    
    Fs = sps;
    [pxx,ft]=pwelch(x,window,noverlap,nfft,Fs,'psd');   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(pxx)); hold on
end
xlim([1e-2,20]);
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('lfe, pwelch, hamming, nfft=n/6')
legend({'PGC','SSIB','SILB'},'location','south');

%%%% 2.1.3 to lfe, more and more shorter window, more averging
figure
for ista = 1:nsta
    x = lfe(ista,:);
    nfft = length(x)/8;   % means no zero-padding
    window = hamming(nfft);
    noverlap = nfft/2;    
    Fs = sps;
    [pxx,ft]=pwelch(x,window,noverlap,nfft,Fs,'psd');   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(pxx)); hold on
end
xlim([1e-2,20]);
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('lfe, pwelch, hamming, nfft=n/8')
legend({'PGC','SSIB','SILB'},'location','south');

%%%% 2.2 to tremor
figure
for ista = 1:nsta
    x = tremor(ista,:);
    nfft = length(x)/4;   % means no zero-padding
    window = hamming(nfft);
    noverlap = nfft/2;    
    Fs = sps;
    [pxx,ft]=pwelch(x,window,noverlap,nfft,Fs,'psd');   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(pxx)); hold on
end
xlim([1e-2,20]);
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('tremor, pwelch, hamming, nfft=n/4')
legend({'PGC','SSIB','SILB'},'location','south');

%%%% 2.2.1 to tremor, more averging
figure
for ista = 1:nsta
    x = tremor(ista,:);
    nfft = 600;   % means no zero-padding
    window = hamming(nfft);
    noverlap = nfft/2;    
    Fs = sps;
    [pxx,ft]=pwelch(x,window,noverlap,nfft,Fs,'psd');   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(pxx)); hold on
end
xlim([1e-2,20]);
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('tremor, pwelch, hamming, nfft=600')
legend({'PGC','SSIB','SILB'},'location','south');


%%%% 2.2.2 to tremor, less averging
figure
for ista = 1:nsta
    x = tremor(ista,:);
    nfft = pow2(nextpow2(4000));   % means no zero-padding
    window = hamming(nfft);
    noverlap = nfft/2;    
    Fs = sps;
    [pxx,ft]=pwelch(x,window,noverlap,nfft,Fs,'psd');   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(pxx)); hold on
end
xlim([1e-2,20]);
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('tremor, pwelch, hamming, nfft=4096')
legend({'PGC','SSIB','SILB'},'location','south');

%%%% 2.2.3 to tremor, more and more averging
figure
for ista = 1:nsta
    x = tremor(ista,:);
    nfft = 256;   % means no zero-padding
    window = hamming(nfft);
    noverlap = nfft/2;    
    Fs = sps;
    [pxx,ft]=pwelch(x,window,noverlap,nfft,Fs,'psd');   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(pxx)); hold on
%     [pxx,ft]=pwelch(x,window,noverlap,nfft,Fs,'power');   % pxx and ft has the length of nfft/2+1, if nfft is even
%     loglog(ft, pxx); hold on
end
xlim([1e-2,20]);
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('tremor, pwelch, hamming, nfft=256')
legend({'PGC','SSIB','SILB'},'location','south');


%%%% 2.2.4 to tremor, more and more averging, just power spectrum, i.e. amplitude
figure
for ista = 1:nsta
    x = tremor(ista,:);
    nfft = 256;   % means no zero-padding
    window = hamming(nfft);
    noverlap = nfft/2;    
    Fs = sps;
    [pxx,ft]=pwelch(x,window,noverlap,nfft,Fs,'power');   % pxx and ft has the length of nfft/2+1, if nfft is even
    loglog(ft, pxx); hold on
end
xlim([1e-2,20]);
xlabel('frequency (Hz)');
ylabel('Amplitude ');
title('tremor, pwelch, hamming, nfft=256')
legend({'PGC','SSIB','SILB'},'location','south');


%%%% 2.2.5 to tremor, try to use a appropriate nfft
figure
for ista = 1:nsta
    x = tremor(ista,:);
    nfft = 1024;   % means no zero-padding
    window = hamming(nfft);
    noverlap = nfft/2;    
    Fs = sps;
    [pxx,ft]=pwelch(x,window,noverlap,nfft,Fs,'psd');   % pxx and ft has the length of nfft/2+1, if nfft is even
    semilogx(ft, 10*log10(pxx)); hold on
%     [pxx,ft]=pwelch(x,window,noverlap,nfft,Fs,'power');   % pxx and ft has the length of nfft/2+1, if nfft is even
%     loglog(ft, pxx); hold on
end
xlim([1e-2,20]);
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('tremor, pwelch, hamming, nfft=256')
legend({'PGC','SSIB','SILB'},'location','south');


%% 3 pmtm, can't divide data into several windows
%%%% 3.1 to lfe, 
figure2
for ista = 1:nsta
    x = lfe(ista,:);
    nfft = pow2(nextpow2(length(x))-1);   % means no zero-padding
    nw = 4;  % For multi-taper spectral analysis, typical choices for NW are 2, 5/2, 3, 7/2, or 4. 
    Fs = sps;
    [pxx,ft]=pmtm(x,nw,nfft,Fs);
    semilogx(ft, 10*log10(pxx)); hold on
end
xlim([1e-2,20]);
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('lfe, pmtm, hamming, nfft=2^i(<=n), nw=4')
legend({'PGC','SSIB','SILB'},'location','south');


%%% NOTE: pmtm requires huge computation when dealing with large data sets, and will incease with
%%% the number of nw, it seems that it is not suitable for the entire records
% %%%% 3.2 to entire tremor 
% figure2
% for ista = 1:nsta
%     x = tremor(ista,:);
%     nfft = length(x);   % means no zero-padding
%     nw = 4;
%     Fs = sps;
%     [pxx,ft]=pmtm(x,nw,nfft,Fs);
%     semilogx(ft, 10*log10(pxx)); hold on
% end
% xlim([1e-2,20]);
% xlabel('frequency (Hz)');
% ylabel('PSD (dB/Hz)');
% title('tremor, pmtm, hamming, nfft=n, nw=4')
% legend({'PGC','SSIB','SILB'},'location','south');

%% 3.2 to tremor detections
figure2
pxxmat1=[];
pxxmat2=[];
pxxmat3=[];
for i=1:nwin
    x=tredet1(i,:);
    nfft = pow2(nextpow2(length(x))-1);
    nw=4;   
    Fs = sps;
    [pxx,ft]=pmtm(x,nw,nfft,Fs);
    semilogx(ft, 10*log10(pxx),'color',[0.7 0.7 0.7]); hold on  % hold option must be after semilog
%     p1=semilogx(ft, 10*log10(pxx),'color','k');
%     p1.Color(4) = 0.5;
    pxxmat1(i,:) = pxx';
    
    x=tredet2(i,:);
    [pxx,ft]=pmtm(x,nw,nfft,Fs);
    semilogx(ft, 10*log10(pxx),'color',[0.7 1 1]);
%     p2=semilogx(ft, 10*log10(pxx),'color','b');
%     p2.Color(4) = 0.5;
    pxxmat2(i,:) = pxx';
    
    x=tredet3(i,:);
    [pxx,ft]=pmtm(x,nw,nfft,Fs);
    semilogx(ft, 10*log10(pxx),'color',[1 0.7 0.7]);
%     p3=semilogx(ft, 10*log10(pxx),'color','r');
%     p3.Color(4) = 0.5;
    pxxmat3(i,:) = pxx';
end
pxxm1 = median(pxxmat1,1);
ln1=semilogx(f.ax(1),ft, 10*log10(pxxm1),'k','linew',1.5,'marker','o','markers',2,'markerf','k');

pxxm2 = median(pxxmat2,1);
ln2=semilogx(f.ax(2),ft, 10*log10(pxxm2),'b','linew',1.5,'marker','o','markers',2,'markerf','b');

pxxm3 = median(pxxmat3,1);
ln3=semilogx(f.ax(3),ft, 10*log10(pxxm3),'r','linew',1.5,'marker','o','markers',2,'markerf','r');
xlim(f.ax(1),[5e-2,50]);
xlim(f.ax(2),[5e-2,50]);
xlim(f.ax(3),[5e-2,50]);

xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('tremor detections,pmtm,nfft=2^i(<=n),nw=4')
legend([ln1,ln2,ln3],{'PGC','SSIB','SILB'},'Location','south');

%% 3.3 to tremor records
% break continous records into segments
wlen = 18*sps;
nw = floor(size(tremoruse,2)/wlen);
tre1=zeros(nw,wlen);
tre2=zeros(nw,wlen);
tre3=zeros(nw,wlen);
for i=1:nw
    tre1(i,:) = tremoruse(1,(i-1)*wlen+1:i*wlen);
    tre2(i,:) = tremoruse(2,(i-1)*wlen+1:i*wlen);
    tre3(i,:) = tremoruse(3,(i-1)*wlen+1:i*wlen);
end
pxxmat1=[];
pxxmat2=[];
pxxmat3=[];

f.fig=figure;
set(f.fig,'Position',[scrsz(2,1)*0.9 0 scrsz(2,3)*0.6 scrsz(2,4)*0.8]);
f.ax(1)=subplot(3,1,1);
f.ax(2)=subplot(3,1,2);
f.ax(3)=subplot(3,1,3);
for i=1:nw
    x=tre1(i,:);
    nfft = pow2(nextpow2(length(x))-1);
    nw=4;   
    Fs = sps;
    [pxx,ft]=pmtm(x,nw,nfft,Fs);
    semilogx(f.ax(1),ft,10*log10(pxx),'color',[0.7 0.7 0.7]); hold(f.ax(1),'on');
%     p1=semilogx(ft, 10*log10(pxx),'color','k');
%     p1.Color(4) = 0.5;
    pxxmat1(i,:) = pxx';
    
    x=tre2(i,:);
    [pxx,ft]=pmtm(x,nw,nfft,Fs);
    semilogx(f.ax(2),ft,10*log10(pxx),'color',[0.7 1 1]); hold(f.ax(2),'on');
%     p2=semilogx(ft, 10*log10(pxx),'color','b');
%     p2.Color(4) = 0.5;
    pxxmat2(i,:) = pxx';
    
    x=tre3(i,:);
    [pxx,ft]=pmtm(x,nw,nfft,Fs);
    semilogx(f.ax(3),ft,10*log10(pxx),'color',[1 0.7 0.7]); hold(f.ax(3),'on');
%     p3=semilogx(ft, 10*log10(pxx),'color','r');
%     p3.Color(4) = 0.5;
    pxxmat3(i,:) = pxx';
end
pxxm1 = median(pxxmat1,1);
ln1=semilogx(f.ax(1),ft, 10*log10(pxxm1),'k','linew',1.5,'marker','o','markers',2,'markerf','k');

pxxm2 = median(pxxmat2,1);
ln2=semilogx(f.ax(2),ft, 10*log10(pxxm2),'b','linew',1.5,'marker','o','markers',2,'markerf','b');

pxxm3 = median(pxxmat3,1);
ln3=semilogx(f.ax(3),ft, 10*log10(pxxm3),'r','linew',1.5,'marker','o','markers',2,'markerf','r');
xlim(f.ax(1),[5e-2,50]);
xlim(f.ax(2),[5e-2,50]);
xlim(f.ax(3),[5e-2,50]);

xlabel(f.ax(1),'frequency (Hz)');
ylabel(f.ax(1),'PSD (dB/Hz)');
title(f.ax(1),'18-s tremor records,pmtm,nfft=2^i(<=n),nw=4')
legend(f.ax(1),[ln1,ln2,ln3],{'PGC','SSIB','SILB'},'location','south');


%% 4 pchave, compared to ptmt, it has averaging similar as pwelch
%%%% 4.1 to lfe
figure
for ista = 1:nsta
    x = lfe(ista,:);    
    nfft = 256;
    lwin = 256;
    olap = 50;    
    Fs = sps;
    [pxx,ft,~,~,~,~,pxxnon]=pchave(x,lwin,olap,nfft,Fs,'MAD','dpss');
    semilogx(ft, 10*log10(pxx)); hold on
%     semilogx(ft, 10*log10(pxxnon),'--'); hold on
end
xlim([1e-2,20]);
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('lfe, pchave, dpss, nfft=n, nw=4, MAD')
legend({'PGC','SSIB','SILB'},'location','south');


%%%% 4.1 to lfe
figure
for ista = 1:nsta
    x = lfe(ista,:);    
    nfft = 256;
    lwin = 256;
    olap = 50;    
    Fs = sps;
    [pxx,ft]=pchave(x,lwin,olap,nfft,Fs,'IQ','dpss');
    semilogx(ft, 10*log10(pxx)); hold on
end
xlim([1e-2,20]);
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('lfe, pchave, dpss, nfft=n, nw=4, IQ')
legend({'PGC','SSIB','SILB'},'location','south');



%%%% 4.2 to tremor detections
figure
pxxmat1=[];
pxxmat2=[];
pxxmat3=[];
for i=1:nwin
    x=tredet1(i,:);
    nfft = 64;
    lwin = 64;
    olap = 50;    
    Fs = sps;
    [pxx,ft]=pchave(x,lwin,olap,nfft,Fs,'MAD','dpss');
    semilogx(ft, 10*log10(pxx),'color',[0.7 0.7 0.7]); hold on  % hold option must be after semilog
%     p1=semilogx(ft, 10*log10(pxx),'color','k');
%     p1.Color(4) = 0.5;
    pxxmat1(i,:) = pxx';
    
    x=tredet2(i,:);
    [pxx,ft]=pchave(x,lwin,olap,nfft,Fs,'MAD','dpss');
    semilogx(ft, 10*log10(pxx),'color',[0.7 1 1]);
%     p2=semilogx(ft, 10*log10(pxx),'color','b');
%     p2.Color(4) = 0.5;
    pxxmat2(i,:) = pxx';
    
    x=tredet3(i,:);
    [pxx,ft]=pchave(x,lwin,olap,nfft,Fs,'MAD','dpss');
    semilogx(ft, 10*log10(pxx),'color',[1 0.7 0.7]);
%     p3=semilogx(ft, 10*log10(pxx),'color','r');
%     p3.Color(4) = 0.5;
    pxxmat3(i,:) = pxx';
end
pxxm1 = median(pxxmat1,1);
ln1=semilogx(f.ax(1),ft, 10*log10(pxxm1),'k','linew',1.5,'marker','o','markers',2,'markerf','k');

pxxm2 = median(pxxmat2,1);
ln2=semilogx(f.ax(2),ft, 10*log10(pxxm2),'b','linew',1.5,'marker','o','markers',2,'markerf','b');

pxxm3 = median(pxxmat3,1);
ln3=semilogx(f.ax(3),ft, 10*log10(pxxm3),'r','linew',1.5,'marker','o','markers',2,'markerf','r');
xlim(f.ax(1),[5e-2,50]);
xlim(f.ax(2),[5e-2,50]);
xlim(f.ax(3),[5e-2,50]);

xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('4-s tremor detections,pchave,dpss,nfft=32,nw=4,MAD')
legend([ln1,ln2,ln3],{'PGC','SSIB','SILB'},'Location','south');


%%
%%%% 4.3 to tremor records

% break continous records into segments
wlen = 18*sps;
nw = floor(size(tremoruse,2)/wlen);
tre1=zeros(nw,wlen);
tre2=zeros(nw,wlen);
tre3=zeros(nw,wlen);
for i=1:nw
    tre1(i,:) = tremoruse(1,(i-1)*wlen+1:i*wlen);
    tre2(i,:) = tremoruse(2,(i-1)*wlen+1:i*wlen);
    tre3(i,:) = tremoruse(3,(i-1)*wlen+1:i*wlen);
end
pxxmat1=[];
pxxmat2=[];
pxxmat3=[];

f.fig=figure;
set(f.fig,'Position',[scrsz(2,1)*0.9 0 scrsz(2,3)*0.6 scrsz(2,4)*0.8]);
f.ax(1)=subplot(3,1,1);
f.ax(2)=subplot(3,1,2);
f.ax(3)=subplot(3,1,3);
for i=1:nw
    x=tre1(i,:);
    nfft = 256;
    lwin = 256;
    olap = 50;    
    Fs = sps;
    [pxx,ft]=pchave(x,lwin,olap,nfft,Fs,'MAD','dpss');
    semilogx(f.ax(1),ft,10*log10(pxx),'color',[0.7 0.7 0.7]); hold(f.ax(1),'on');
%     p1=semilogx(ft, 10*log10(pxx),'color','k');
%     p1.Color(4) = 0.5;
    pxxmat1(i,:) = pxx';
    
    x=tre2(i,:);
    [pxx,ft]=pchave(x,lwin,olap,nfft,Fs,'MAD','dpss');
    semilogx(f.ax(2),ft,10*log10(pxx),'color',[0.7 1 1]); hold(f.ax(2),'on');
%     p2=semilogx(ft, 10*log10(pxx),'color','b');
%     p2.Color(4) = 0.5;
    pxxmat2(i,:) = pxx';
    
    x=tre3(i,:);
    [pxx,ft]=pchave(x,lwin,olap,nfft,Fs,'MAD','dpss');
    semilogx(f.ax(3),ft,10*log10(pxx),'color',[1 0.7 0.7]); hold(f.ax(3),'on');
%     p3=semilogx(ft, 10*log10(pxx),'color','r');
%     p3.Color(4) = 0.5;
    pxxmat3(i,:) = pxx';
end
pxxm1 = median(pxxmat1,1);
ln1=semilogx(f.ax(1),ft, 10*log10(pxxm1),'k','linew',1.5,'marker','o','markers',2,'markerf','k');

pxxm2 = median(pxxmat2,1);
ln2=semilogx(f.ax(2),ft, 10*log10(pxxm2),'b','linew',1.5,'marker','o','markers',2,'markerf','b');

pxxm3 = median(pxxmat3,1);
ln3=semilogx(f.ax(3),ft, 10*log10(pxxm3),'r','linew',1.5,'marker','o','markers',2,'markerf','r');
xlim(f.ax(1),[5e-2,50]);
xlim(f.ax(2),[5e-2,50]);
xlim(f.ax(3),[5e-2,50]);

xlabel(f.ax(1),'frequency (Hz)');
ylabel(f.ax(1),'PSD (dB/Hz)');
title(f.ax(1),'18-s tremor records,pchave,dpss,nfft=256,nw=4,MAD')
legend(f.ax(1),[ln1,ln2,ln3],{'PGC','SSIB','SILB'},'location','south');


