%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare the spectra of the templates obtained by Chao and Allan, as it seems like Allan's spectra
% have an obvious peak around the 1/t_dur, but mine are rather flattish.
% try to see why this happens, test for the param used to get the spectra, assuming the different
% templates are obtained in their own way correctly.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/07/30
% Last modified date:   2021/07/30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format short e   % Set the format to 5-digit floating point
clc
clear
close all

[scrsz, res] = pixelperinch(1);

%% choose templates from Chao or Allan
TEMPFLAG = 'Chao';
% TEMPFLAG = 'Allan';

% YSHAPE = 'psd';
YSHAPE = 'amp';

TAPER = 'tukey';
fractap = 0.2;

%%% load Chao's templates 
if strcmp(TEMPFLAG, 'Chao')
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
    
    for ista = 1: nsta
%         fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
%             num2str(templensec), 's_', 'BBDS_', 'opt_Nof_Non_Chao');
%         dstack(ista,:) = load(fname);
%         fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
%             num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao');
%         ccstack(ista,:) = load(fname);
        
        fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
            num2str(templensec), 's_', 'BBDS_', 'opt_Nof_Non_Chao_cat',CATA,'_m1cc0.5-6.5');
        dstack(ista,:) = load(fname);
        fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
            num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao_cat',CATA,'_m1cc0.5-6.5');
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

    
%%% load Chao's templates
elseif strcmp(TEMPFLAG, 'Allan')
    
    temppath = '/Users/admin/Documents/MATLAB/allan/matfils/';
    datfiles=['PGCopt_002_1.25-6.5Hz_2pass_100sps_0.1-15Hz_le14shSRR';  %the templates, 0.5-15 Hz at 100 Hz sampling rate
              'SSIopt_002_1.25-6.5Hz_2pass_100sps_0.1-15Hz_le14shSRR';
              'SILopt_002_1.25-6.5Hz_2pass_100sps_0.1-15Hz_le14shSRR'];
    
    lohf = 0.1;
    hihf = 15;
    
    stas=['PGC '
          'SSIB'
          'SILB'];     % determine the trio and order
    nsta=size(stas,1);
    for ista=1:nsta
        fname = strcat(temppath, datfiles(ista,:));
        tmp = load(fname);
        stack(:, ista) = tmp(:,1);
    end
    sps = 100;
    templen = size(stack,1);
    
    mid = [8273; 8496; 8326];
    wlensec = 20.48;
    wlen = wlensec*sps;
    for ista=1:nsta
        tracebpadd(:, ista)=stack(mid(ista)-wlen/2-sps: mid(ista)+wlen/2+sps-1, ista);
    end
    mean(tracebpadd)
    
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
        
end

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


%% plot the trace and spectra
f.fig = figure;
f.fig.Renderer = 'painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 11;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/15 scrsz(4)/10 widin*res htin*res]);
nrow = 3;
ncol = 1;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).FontSize = 8;
    f.ax(isub).Box = 'on';
    grid(f.ax(isub), 'on');
    f.ax(isub).GridLineStyle = '--';
end

% % reposition
% set(f.ax(2), 'position', [ 0.1, 0.07, 0.85, 0.48]);     % for spectrum
% set(f.ax(1), 'position', [ 0.1, 0.62, 0.85, 0.33]);     % for LF

%%%%%%%% broader band, instead of broadband
ax = f.ax(1);
hold(ax,'on');

time = (0:1/sps:wlensec-1/sps)';
is = round(time(1));  % start time of each win
ien = round(time(end));   % end time of each win

ymami = max([maxamppos; -minampneg], [], 1);
pltscale = 1.2;
ym = pltscale*max(ymami);
axis(ax,[time(1) time(end) -ym ym]);

%annotate the window to find the max/min
plot(ax,[is+(wlen/2-0.5*fixlen)/sps is+(wlen/2-0.5*fixlen)/sps],ax.YLim,'k--',...
    'linew',1.5);
plot(ax,[is+(wlen/2+0.5*fixlen-1)/sps is+(wlen/2+0.5*fixlen-1)/sps],ax.YLim,'k--',...
    'linew',1.5);

%%% plot the seismogram at 3 stations of each win
plot(ax,time,tracebp(:, 1),'r','linew',0.5);
plot(ax,time,tracebp(:, 2),'b','linew',0.5);
plot(ax,time,tracebp(:, 3),'k','linew',0.5);

title(ax,strcat(TEMPFLAG));
text(ax,0.02,0.9,strcat(num2str(lohf),'-',num2str(hihf),{' Hz'}),'fontsize',8,'unit',...
    'normalized','Horizontalalignment','left');
text(ax,0.98,0.9,sprintf('%.4f',ymami(1)),'fontsize',8,'unit',...
    'normalized','Horizontalalignment','right','color','r');
text(ax,0.98,0.75,sprintf('%.4f',ymami(2)),'fontsize',8,'unit',...
    'normalized','Horizontalalignment','right','color','b');
text(ax,0.98,0.6,sprintf('%.4f',ymami(3)),'fontsize',8,'unit',...
    'normalized','Horizontalalignment','right','color','k');

% set(ax,'XTick',[is is+(ien-is)/4 is+(ien-is)*2/4 is+(ien-is)*3/4 ien],'fontsize',8);
ylabel(ax,'Amplitude','fontsize',9);
xlabel(ax,'Time (s)','fontsize',9);

% ax.XTickLabel = [];
%         ax.YTick = [];
hold(ax,'off');


%%%%%%%%%%%%%%%%% subplot, power spectral density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mid = ceil(size(tracebp,1)/2);

psdlen = wlen;
linesty = {'-'};

Fs = sps;

color=['r','b','k'];

jj = 1;

if strcmp(YSHAPE, 'psd')
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


elseif strcmp(YSHAPE, 'amp')
    ax = f.ax(2);
    % use the untapered version, cut a segment of different length centered at the main dipole
    x = tracebp(mid-psdlen(jj)/2+1: mid+psdlen(jj)/2, :);
    % nfft = pow2(nextpow2(size(x,1))-1);
    % nfft = 1024;
    nfft = size(x,1);
    window = hamming(size(x,1));
    nw = 4;

    [psdx,pft] = pmtm(x,nw,nfft,Fs);
    psdxmean = mean(psdx,2);
    psdxmed = median(psdx,2);
    p(1)=loglog(ax,pft,sqrt(psdxmed),linesty{jj},'linewidth',1,...
            'color',color(1),'markers',1.5);
    hold(ax, 'on');
    grid(ax, 'on');
    p(2)=loglog(ax,pft,sqrt(psdxmean),linesty{jj},'linewidth',1,...
            'color',color(2),'markers',1.5);    
    xlim(ax,[1e-1, 1.2*Fs/2]);
    ylim(ax,[5e-4, 1e-1]);
    title(ax,sprintf('%.1f-%d Hz, pmtm, sps: %d, nfft: %d',...
        lohf,hihf,sps,nfft));
    xlabel(ax,'Frequency (Hz)');
    ylabel(ax,'Amp/Hz');
    lgd = legend(ax,p,{'median','mean'},'location','southwest','fontsize',8);
    hold(ax, 'off');

    ax = f.ax(3);
    % use the tapered version, cut a segment of different length centered at the main dipole
%     x = tracebptap(mid-psdlen(jj)/2+1: mid+psdlen(jj)/2, :);
    x = tracebp(mid-psdlen(jj)/2+1: mid+psdlen(jj)/2, :);
    [psdx,pft] = periodogram(x,window,nfft,Fs);
    psdxmean = mean(psdx,2);
    psdxmed = median(psdx,2);
    p(1)=loglog(ax,pft,sqrt(psdxmed),linesty{jj},'linewidth',1,...
            'color',color(1),'markers',1.5);
    hold(ax, 'on');    
    grid(ax, 'on');
    p(2)=loglog(ax,pft,sqrt(psdxmean),linesty{jj},'linewidth',1,...
            'color',color(2),'markers',1.5);
    xlim(ax,[1e-1, 1.2*Fs/2]);
    ylim(ax,[5e-4 1e-1]);
    title(ax,sprintf('%.1f-%d Hz, periodogram, sps: %d, nfft: %d, Hamming',...
    lohf,hihf,sps,nfft));
    xlabel(ax,'Frequency (Hz)');
    ylabel(ax,'Amp/Hz');
    lgd = legend(ax,p,{'median','mean'},'location','southwest','fontsize',8);
    hold(ax, 'off');

end




       




















