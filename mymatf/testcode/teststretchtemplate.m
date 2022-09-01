% teststretchtemplate
% try to stretch the template and see how that would vary the spectra
% note that this stretching would change not only the 'signal' but also the 'noise'

format short e   % Set the format to 5-digit floating point
clc
clear
close all

[scrsz, res] = pixelperinch(1);


%%% this part of codes are similar to that in 'testcutdipole.m' with fair modifications
temppath = '/Users/admin/Documents/MATLAB/';

dstack = [];
ccstack = [];
sps = 100;
templensec = 60;

fam = '002';
disp(fam);

stas=['PGC '
    'SSIB'
    'SILB'];     % determine the trio and order
nsta=size(stas,1);         %  number of stations

%     CATA = 'old';
CATA = 'new';

% YSHAPE = 'psd';
YSHAPE = 'amp';

TAPER = 'tukey';
fractap = 0.2;

for ista = 1: nsta
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
        num2str(templensec), 's_', 'BBDS_', 'opt_Nof_Non_Chao_cat',CATA,'_m1');
    dstack(ista,:) = load(fname);
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
        num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao_cat',CATA,'_m1');
    ccstack(ista,:) = load(fname);
end
stack = ccstack;

templen = size(stack,2);
mid = templen/2;

wlensec = 20.48;
wlen = wlensec*sps;
tracebb = stack(:, mid-wlen/2-sps: mid+wlen/2+sps-1)';
%     tracebb = stack(:, mid-wlen/4-sps: mid+wlen/4*3+sps-1)';
mean(tracebb)

%%% detrend data, == remove mean, remove linear trend
tracebb = detrend(tracebb);
mean(tracebb)

lohf = 0.1;
hihf = 15;
npo = 2;
npa = 2;

%%% taper data
if strcmp(TAPER, 'tukey')
    % taper with tukeywin, which is actually a tapered cosine window
    %tapered length is adaptative to frequency, maybe at least longer than one full period length of
    %the lowest frequency
    %     fractap = round(1/lohf*sps)/size(tracebb,1)*2;
    %     fractap = max(fractap,0.1); % if fractap is >=1, n-point von Hann window is returned
    w = tukeywin(size(tracebb,1),fractap);
    tracebbtap = w.* tracebb;
end

%%% bandpass the broadband (unfiltered) data
tracebpadd = tracebbtap;
for ista = 1: nsta
    tracebpadd(:,ista) = Bandpass(tracebbtap(:,ista), sps, lohf, hihf, npo, npa, 'butter');
end

[maxses,imaxses]=max(tracebpadd,[],1);
[minses,iminses]=min(tracebpadd,[],1);
spread=maxses-minses;

for ista=1:nsta
    tracebpadd(:, ista)=tracebpadd(:, ista)/spread(ista);
end

figure
ax = gca;
hold on
%%% plot the seismogram at 3 stations of each win
plot(ax,tracebpadd(:, 1),'r','linew',0.5);
plot(ax,tracebpadd(:, 2),'b','linew',0.5);
plot(ax,tracebpadd(:, 3),'k','linew',0.5);

% re-align the filtered traces, because the alignment was in LF, which may not be true for
% higher frequency
mshiftadd = sps/5;    % maximum allowed shift between 2 traces
mid = ceil(size(tracebpadd,1)/2);
fixlen = 2*sps;
loffmax = 4;
ccmin = 0.3;  % 0.4/0.35
iup = 1;    % times of upsampling
[off12add,off13add,cc,iloopoff,loopoff] = constrained_cc_interp(tracebpadd',mid,...
    fixlen,mshiftadd,loffmax,ccmin,iup);

% align the records
istart = sps+1;
iend = wlen+sps;
tracebp = [];
tracebp(:, 1) = tracebpadd(istart: iend, 1);
tracebp(:, 2) = tracebpadd(istart-round(off12add): iend-round(off12add), 2);
tracebp(:, 3) = tracebpadd(istart-round(off13add): iend-round(off13add), 3);
mean(tracebp)

%%% since it is possibily realigned, for safety, detrend and taper the aligned traces again!
tracebp = detrend(tracebp);
mean(tracebp)

if strcmp(TAPER, 'tukey')
    % taper with tukeywin, which is actually a tapered cosine window
    %tapered length is adaptative to frequency, maybe at least longer than one full period length of
    %the lowest frequency
    w = tukeywin(size(tracebp,1),fractap);
    tracebptap = w.* tracebp;
end
mean(tracebptap)

%in case there is other unpredicted noise, use only the short segment around the dipole
[maxamppos,indmax] = max(tracebp(max(wlen/2-0.5*fixlen,1): ...
    min(wlen/2+0.5*fixlen-1,wlen),:), [],1);
[minampneg,indmin] = min(tracebp(max(wlen/2-0.5*fixlen,1): ...
    min(wlen/2+0.5*fixlen-1,wlen),:), [],1);


%% stretch the template
% sps = 20;
scale = 2;
len = size(tracebp,1);
samp = 1:1:len;
sampstr = 1/scale:1/scale:len;
tracebpstr = interp1(samp, tracebp, sampstr, 'spline') /sqrt(scale);
figure
for ista = 1: nsta
    subplot(3,1,ista)
    hold on
    plot(samp-len/2, tracebp(:,ista), 'k','linew',2);
    plot(sampstr*scale-len*scale/2, tracebpstr(:,ista), 'b','linew',2);
end

%% cross-correlation
offsetmax = sps/4;
figure
for ista = 1: nsta
    subplot(3,2,(ista-1)*2+1)
    hold on
    x1 = tracebp(len/2-len/4+1: len/2+len/4, ista);
    x2 = tracebpstr(len*scale/2-len/4+1: len*scale/2+len/4, ista);
    plot(1:1:len/2, detrend(x1), 'k','linew',2);
    plot(1:1:len/2, detrend(x2), 'b','linew',2);
    
    [coef, lag] = xcorr(detrend(x1), detrend(x2), offsetmax, 'coeff');
    [maxcoef(ista), idx] = max(coef);
    lagsamp(ista) = lag(idx);
    x1ali = tracebp(len/2-len/4+1+lagsamp(ista): len/2+len/4+lagsamp(ista), ista);
    subplot(3,2,(ista-1)*2+2)
    hold on
    plot(1:1:len/2, detrend(x1ali), 'k','linew',2);
    plot(1:1:len/2, detrend(x2), 'b','linew',2);
end
maxcoef
lagsamp


%% plot spectra
%%%%%%%%%%%%%%%%% subplot, power spectral density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs = sps;

color=['r','b','k'];

jj = 1;

if strcmp(YSHAPE, 'amp')
    figure
    for ista = 1: nsta
        subplot(3,1,ista)
        ax = gca;
        % use the untapered version, cut a segment of different length centered at the main dipole
        x = tracebp(:, ista);
        nfft = size(x,1);
        nw = 4;
        
        [psdx,pft] = pmtm(x,nw,nfft,Fs);
        p(1)=loglog(ax,pft,sqrt(psdx),'-','linewidth',2,...
            'color','k','markers',1.5);
        hold(ax, 'on');
        grid(ax, 'on');
        xlim(ax,[1e-1, 1.2*Fs/2]);
        ylim(ax,[5e-4, 1e-1]);
        
        x = tracebpstr(:, ista);
        nfft = size(x,1);
        nw = 4;
        
        [psdx,pft] = pmtm(x,nw,nfft,Fs);
        p(2)=loglog(ax,pft,sqrt(psdx),'-','linewidth',2,...
            'color','b','markers',1.5);
        p(3)=loglog(ax,pft*scale,sqrt(psdx),'-','linewidth',2,...
            'color','r','markers',1.5);
%         p(4)=loglog(ax,pft*scale,sqrt(psdx/scale),'-','linewidth',2,...
%             'color','c','markers',1.5);
        plot(ax, [2 2],ax.YLim,'r--','linewidth',1 );
        plot(ax, [5 5],ax.YLim,'r--','linewidth',1 );
        plot(ax, [2 2]/scale,ax.YLim,'b--','linewidth',1 );
        plot(ax, [5 5]/scale,ax.YLim,'b--','linewidth',1 );
        title(ax,sprintf('%.1f-%d Hz, pmtm, sps: %d',...
            lohf,hihf,sps));
        xlabel(ax,'Frequency (Hz)');
        ylabel(ax,'Amp /Hz');
%         lgd = legend(ax,p,{'original','stretched','freq*scale','amp/scale'},...
%             'location','southwest','fontsize',10);
        lgd = legend(ax,p,{'original','stretched','freq*scale'},...
            'location','southwest','fontsize',10);
        hold(ax, 'off');
        
    end
 
    figure
    for ista = 1: nsta
        subplot(3,1,ista)
        ax = gca;
        % use the tapered version, cut a segment of different length centered at the main dipole
        %     x = tracebptap(mid-psdlen(jj)/2+1: mid+psdlen(jj)/2, :);
        x = tracebp(:, ista);
        nfft = size(x,1);
        window = hanning(size(x,1));
        [psdx,pft] = periodogram(x,window,nfft,Fs);
        p(1)=loglog(ax,pft,sqrt(psdx),'-','linewidth',2,...
            'color','k','markers',1.5);
        hold(ax, 'on');
        grid(ax, 'on');
        
        x = tracebpstr(:, ista);
        nfft = size(x,1);
        window = hanning(size(x,1));
        [psdx,pft] = periodogram(x,window,nfft,Fs);
        p(2)=loglog(ax,pft,sqrt(psdx),'-','linewidth',2,...
            'color','b','markers',1.5);
        p(3)=loglog(ax,pft*scale,sqrt(psdx),'-','linewidth',2,...
            'color','r','markers',1.5);
%         p(4)=loglog(ax,pft*scale,sqrt(psdx/scale),'-','linewidth',2,...
%             'color','c','markers',1.5);

        xlim(ax,[1e-1, 1.2*Fs/2]);
        ylim(ax,[5e-4 1e-1]);
        title(ax,sprintf('%.1f-%d Hz, periodogram, sps: %d, Hann',...
            lohf,hihf,sps));
        xlabel(ax,'Frequency (Hz)');
        ylabel(ax,'Amp /Hz');
%         lgd = legend(ax,p,{'original','stretched','freq*scale','amp/scale'},...
%             'location','southwest','fontsize',10);
        lgd = legend(ax,p,{'original','stretched','freq*scale'},...
            'location','southwest','fontsize',10);
        hold(ax, 'off');
    
    end
    
elseif strcmp(YSHAPE, 'psd')
    ax = f.ax(2);
    hold(ax, 'on');
    % use the untapered version, cut a segment of different length centered at the main dipole
    x = tracebp(mid-psdlen(jj)/2+1: mid+psdlen(jj)/2, :);
    % nfft = pow2(nextpow2(size(x,1))-1);
    % nfft = 1024;
    nfft = size(x,1);
    window = hann(size(x,1));
    nw = 4;
    
    [psdx,pft] = pmtm(x,nw,nfft,Fs);
    psdxmean = mean(psdx,2);
    psdxmed = median(psdx,2);
    
    for ista = 1: nsta
        p(ista)=plot(ax,pft,pow2db(psdx(:,ista)),linesty{jj},'linewidth',1,...
            'color',color(ista),'markers',1.5);
    end
    
    %     p(1)=loglog(ax,pft,pow2db(psdxmed),linesty{jj},'linewidth',1,...
    %             'color',color(1),'markers',1.5);
    %     hold(ax, 'on');
    %     p(2)=loglog(ax,pft,pow2db(psdxmean),linesty{jj},'linewidth',1,...
    %             'color',color(2),'markers',1.5);
    
    ax.XScale = 'log';
    xlim(ax,[1e-1, 1.2*Fs/2]);
    ylim(ax,[-150,-20]);
    title(ax,sprintf('%.1f-%d Hz, pmtm, sps: %d, nfft: %d, Hann',...
        lohf,hihf,sps,nfft));
    xlabel(ax,'Frequency (Hz)');
    ylabel(ax,'PSD (dB/Hz)');
    %     lgd = legend(ax,p,{'median','mean'},'location','southwest','fontsize',8);
    lgd = legend(ax,p,stas,'location','southwest','fontsize',8);
    hold(ax, 'off');
    
    ax = f.ax(3);
    hold(ax, 'on');
    % use the tapered version, cut a segment of different length centered at the main dipole
    %     x = tracebptap(mid-psdlen(jj)/2+1: mid+psdlen(jj)/2, :);
    x = tracebp(mid-psdlen(jj)/2+1: mid+psdlen(jj)/2, :);
    
    [psdx,pft] = periodogram(x,window,nfft,Fs);
    psdxmean = mean(psdx,2);
    psdxmed = median(psdx,2);
    
    for ista = 1: nsta
        p(ista)=plot(ax,pft,pow2db(psdx(:,ista)),linesty{jj},'linewidth',1,...
            'color',color(ista),'markers',1.5);
    end
    
    %     p(1)=loglog(ax,pft,pow2db(psdxmed),linesty{jj},'linewidth',1,...
    %             'color',color(1),'markers',1.5);
    %     hold(ax, 'on');
    %     p(2)=loglog(ax,pft,pow2db(psdxmean),linesty{jj},'linewidth',1,...
    %             'color',color(2),'markers',1.5);
    
    ax.XScale = 'log';
    xlim(ax,[1e-1, 1.2*Fs/2]);
    ylim(ax,[-150,-20]);
    title(ax,sprintf('%.1f-%d Hz, periodogram, sps: %d, nfft: %d, Hann',...
        lohf,hihf,sps,nfft));
    xlabel(ax,'Frequency (Hz)');
    ylabel(ax,'PSD (dB/Hz)');
    %     lgd = legend(ax,p,{'median','mean'},'location','southwest','fontsize',8);
    lgd = legend(ax,p,stas,'location','southwest','fontsize',8);
    hold(ax, 'off');
    
    
end















