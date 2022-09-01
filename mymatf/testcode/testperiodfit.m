% Test the usage of 'periodfit', and compare the result with the 'Bandpass'
%
% Whenever you want to filter a trace, always taper first, this is to avoid
% the artifact from non-periodicity as filtering in Fourier domain assumes 
% the incoming signal is periodic which ends at where it starts, and taper
% would make sure both ends go to 0.
% -- Ideally, I would expect if I bandpass the data with a decent-width filter,
% the filtered signal should be very close to the prediction by inverting
% the cos/sin contribution in the same frequency range; and it should be 
% more accurate if i use the shorter segment around the dipole rather than
% a longer trace.
% -- However, the direct experiment shows that it is more accurate when the
% signal is long and period spacing used in 'periodfit' is smaller enough,
% and of course if the spacing is too small, the inversion matrix would be
% ill-conditioned (the length effect i think is caused by the alising effect 
% of the finite signal rather than the infinite long one). The recovery of
% bandpassed signal is also the best for the main arrival, which is exactly
% what we want.
% -- Also, i am surprised that if the feeding signal for periodfit is short,
% the main part you want to recover almost always fail, even the inverison is
% stable. So it seems that if you want to recover a short sement of signal,
% you need to feed a much longer window, also, long the window, smaller the 
% period spacing
% -- I think maybe 'periodfit' is more useful when the filter bandpass is too
% narrow so that filtering would be numerically unstable, but you still want
% to know the contribution from a specific period, then you can directly invert
% for that period for that small segment of the signal
% -- Finally, there is catch fundamentally in this comparison which makes it
% a little unfair. That is, bandpass designs a filter which is NOT completely
% flat inside the corners, meaning NOT all freqs within corners get 100% gain,
% and the freqs outside the corners would receive a strong decay, but NOT zero
% gain, so there is definitely leaking of energy at other freqs.  The periodfit
% asks for the contribution inside several periods, but you could never include
% all possible periods in the range, and when the spacing is too small, inversion
% would fail, even through here the corners can be exact. In terms of accuracy
% in the exact content of a specific frequency range, maybe 'periodfit' is better
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/05/04
% Last modified date:   2021/05/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
clc
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);
% scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]
% wid=scrsz(3);       % width
% hite=scrsz(4);       % height

workpath = getenv('ALLAN');
%%% choose the data path here
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');

% FLAG = 'PGC';
FLAG = 'TWKB';

disp(FLAG);

if strcmp(FLAG, 'TWKB')
    temppath = strcat(datapath, '/templates/LZBtrio/');
    rstpath = strcat(datapath, '/LZBtrio');
    
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
            '017';
            '006';
            '001';
            ];   
    nfam = length(nfampool);

    stas=['TWKB '
          'LZB  '
          'MGCB '];     % determine the trio and order

elseif strcmp(FLAG, 'PGC')    
    temppath = strcat(datapath, '/templates/PGCtrio/');
    rstpath = strcat(datapath, '/PGCtrio');
    
    nfampool = ['002';
               ];
    nfam = length(nfampool);

    stas=['PGC  '
          'SSIB '
          'SILB '];

end
    
% number of used stations
nsta=size(stas,1);         %  number of stations


%% make template or read template
remake = 0;  % re-make the template or not, 1/0
dstack = [];
ccstack = [];
sps = 80;
templensec = 60;
% for ifam = 1: size(nfampool,1)
    ifam=6;
    fam = nfampool(ifam, :);
    disp(fam);
    
    if remake
        ccmethod = 2; % only using lfes passed the cc thresholds (ccmethod=1); using all lfes (ccmethod=2)  
        plflag = 0;
        [dstack, ccstack] = mk_bbtemp_LZB(fam,sps,ccmethod,plflag);
        % write into files
        for ista = 1: nsta
            fid = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
                num2str(sps), 'sps_', num2str(templensec), 's_', ...
                'BBDS_', 'opt_Nof_Non_Chao'), 'w+');  % bandpassed direct stack, no filter, no norm
            
            fprintf(fid, '%f \n', dstack(ista, :)');
            fclose(fid);
            
            fid = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
                num2str(sps), 'sps_', num2str(templensec), 's_', ...
                'BBCCS_', 'opt_Nof_Non_Chao'), 'w+');  % bandpassed cc stack, no filter, no norm
            fprintf(fid, '%f \n', ccstack(ista, :)');
            fclose(fid);
        end
        ind = [5 4 6];
        %     stack = dstack(ind, :);
        stack = ccstack(ind, :);
        
    else
        for ista = 1: nsta
            fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
                num2str(templensec), 's_', 'BBDS_', 'opt_Nof_Non_Chao');
            dstack(ista,:) = load(fname);
            fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
                num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao');
            ccstack(ista,:) = load(fname);
        end
        stack = ccstack;
    end
    
    templen = size(stack,2);
    mid = templen/2;
    
    wlensec = 60;
    wlen = wlensec*sps;   
    tracebb = stack(:, mid-wlen/2+1: mid+wlen/2)';
    
    % taper with tukeywin, which is actually a tapered cosine window
    w = tukeywin(size(tracebb,1),0.1);
    for ista = 1:nsta
        tracebb(:,ista) = tracebb(:,ista).*w;        
    end
%     
%     figure
%     plot(w)
    
%     %% plot the templates
%     %%% plot the 1-step template
%     %%% figure 3
%     f.fig = figure;
%     f.fig.Renderer = 'painters';
%     color=['r','b','k'];
%     ax = gca;
%     hold(ax, 'on');
%     for ista = 1: nsta
%         plot(ax,(1:wlen)/sps,tracebb(:,ista), 'linewidth', 1,'color',...
%             color(ista)); hold on
%         %     plot(stack(ista, mid-5*sps: mid+5*sps) + 1* ista, 'linewidth', 1,'color',...
%         %          color(ista)); hold on
% %         text(ax,1,0,stas(ista,:),'fontsize',9);
%     end
%     text(ax,0.95,0.9,'Broadband','fontsize',10,'unit','normalized','Horizontalalignment','right');
%     text(ax,0.95,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
%          'right','fontsize',10);
%     text(ax,0.03,0.92,'a','FontSize',9,'unit','normalized','EdgeColor','k','Margin',2);
%     text(ax,0.03,0.1,fam,'fontsize',12,'unit','normalized','Horizontalalignment','left');
%     text(ax,0.95,0.1,FLAG,'fontsize',12,'unit','normalized','Horizontalalignment','right');
%     xlabel(ax,'Time (s)','fontsize',10);
%     ylabel(ax,'Amplitude','fontsize',10);
%     xlim(ax,[0 wlensec]);
% %     ylim(ax,[0 4]);
%     ax.Box='on';
%     grid(ax, 'on');
%     ax.GridLineStyle = '--';
%     
%     %     wavehil = abs(hilbert(tracefile(winlen*(ndet-1)+1:winlen*ndet,2:4)));
%     [envup, envlo] = envelope(tracebb);
%     
%     medenvup = median(envup,2);
%     medenvlo = median(envlo,2);
%     
%     plot(ax,(1:wlen)/sps,medenvup,'color',[0.7 0.7 0.7],'linew',1.5);
%     plot(ax,(1:wlen)/sps,medenvlo,'color',[0.7 0.7 0.7],'linew',1.5);
%     hold(ax, 'off');


%% check the contribution of specific frequencies by inversion, periodfit
    
    f.fig = figure;
    f.fig.Renderer = 'painters';
    widin = 6;  % maximum width allowed is 8.5 inches
    htin = 9;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
    nrow = 3;
    ncol = 1;
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
        f.ax(isub).FontSize = 8;
        f.ax(isub).Box = 'on';
        grid(f.ax(isub), 'on');
        f.ax(isub).GridLineStyle = '--';
    end
    color=['r','b','k'];

    %%%%%%%% broadband
    tracef = tracebb;
    ax = f.ax(1);
    hold(ax,'on');
    
    maxamppos = max(tracef, [],1);
    minampneg = min(tracef, [],1);
    ymami = max([maxamppos; -minampneg], [], 1);
    ym = 1.1*max(ymami);
    
    is = 0;   % start time of each win
    ien= round(wlen/sps);         % end time of each win
    axis(ax,[is ien -ym ym]);
    
    %%% plot the seismogram at 3 stations of each win
    plot(ax,(1:wlen)/sps,tracef(:, 1),'r','linew',0.5);
    plot(ax,(1:wlen)/sps,tracef(:, 2),'b','linew',0.5);
    plot(ax,(1:wlen)/sps,tracef(:, 3),'k','linew',0.5);

    text(ax,0.05,0.9,strcat('Broadband'),'fontsize',9,'unit',...
        'normalized','Horizontalalignment','left');
    text(ax,0.95,0.9,sprintf('%.4f',ymami(1)),'fontsize',8,'unit',...
        'normalized','Horizontalalignment','right','color','r');
    text(ax,0.95,0.75,sprintf('%.4f',ymami(2)),'fontsize',8,'unit',...
        'normalized','Horizontalalignment','right','color','b');
    text(ax,0.95,0.6,sprintf('%.4f',ymami(3)),'fontsize',8,'unit',...
        'normalized','Horizontalalignment','right','color','k');
    
    set(ax,'XTick',is:4:ien,'fontsize',8);
    %         ylabel(ax,'Amplitude','fontsize',9);
    ax.XTickLabel = [];
    ax.YTick = [];
    hold(ax,'off');
    
    octavef = 2.^(0:-1:-1);
    for ifreq = 1: length(octavef)-1
        
        ax = f.ax(ifreq+1);
        hold(ax,'on');
        
        lo = octavef(ifreq+1);
        hi = octavef(ifreq);
        npo = 2;    % npo being 2 or 3 is proper
        npa = 2;
        %             if diff(2*[lohf hihf]/sps)>=0.01  % the bandwidth is moderate
        for ista = 1: nsta
            tracef(:,ista) = Bandpass(tracebb(:,ista), sps, lo, hi, npo, npa, 'butter');
        end
        %             else   % the bandwidth is too narrow, thus the filtering may be less reliable
        %                 for ista = 1: nsta
        %                     nsps = 5;
        %                     [num, denom] = rat(nsps/sps);
        %                     tracenofds(:,ista+1) = resample(tracenof(:,ista+1),num,denom);   % times = num/denom
        %                     tracehfds(:,ista+1) = Bandpass(tracenofds(:,ista+1), nsps, lohf, hihf, npo, npa, 'butter');
        %                     [num, denom] = rat(sps/nsps);
        %                     tracehf(:,ista+1) = resample(tracehfds(:,ista+1),num,denom);
        %                 end
        %             end
        
        %%%%%%%%%%%%%%%%
        
        maxamppos = max(tracef, [],1);
        minampneg = min(tracef, [],1);
        ymami = max([maxamppos; -minampneg], [], 1);
        ym = 1.1*max(ymami);

        axis(ax,[is ien -ym ym]);
        
        %%% plot the seismogram at 3 stations of each win
        plot(ax,(1:wlen)/sps,tracef(:, 1),'r','linew',0.5);
        plot(ax,(1:wlen)/sps,tracef(:, 2),'b','linew',0.5);
        plot(ax,(1:wlen)/sps,tracef(:, 3),'k','linew',0.5);
        
        %%% plot freq. band
        %         plot(ax,[is+1/hihf is+1/lohf],[-0.8*yma -0.8*yma],'k','linew',3)
        if lo<1 || hi <1
            lop = round(1/lo);
            hip = round(1/hi);
            text(ax,0.05,0.9,strcat(num2str(hip),'-',num2str(lop),{' s'}),'fontsize',9,'unit',...
                'normalized','Horizontalalignment','left');
        else
            text(ax,0.05,0.9,strcat(num2str(lop),'-',num2str(hip),{' Hz'}),'fontsize',9,'unit',...
                'normalized','Horizontalalignment','left');
        end
        text(ax,0.95,0.9,sprintf('%.4f',ymami(1)),'fontsize',8,'unit',...
            'normalized','Horizontalalignment','right','color','r');
        text(ax,0.95,0.75,sprintf('%.4f',ymami(2)),'fontsize',8,'unit',...
            'normalized','Horizontalalignment','right','color','b');
        text(ax,0.95,0.6,sprintf('%.4f',ymami(3)),'fontsize',8,'unit',...
            'normalized','Horizontalalignment','right','color','k');
        
        set(ax,'XTick',is:4:ien,'fontsize',8);
        
        if ifreq ~= length(octavef)-1
            ax.XTickLabel = [];
        else
            %                 ylabel(ax,'Amplitude','fontsize',9);
            xlabel(ax,strcat({'Time (s)'}),'fontsize',9);
        end
%         ax.YTick = [];
        hold(ax,'off');
        
        ax = f.ax(ifreq+2);
        hold(ax,'on');
        Fs = sps;
        nfft = 2048;
        df = Fs/nfft;
%         df = 0.02;
        x = reshape((1:wlen)/sps,[],1);
%         perdx = linspace(1,2,df);
        % dT can not be too small, otherwis the 'inv' would be unstable 
        % Matrix is close to singular or badly scaled
        dT = 0.03;  
        perdx = 1: dT: 2;
%         perdx = 10;
        comp = 'parallel';
%         comp = 'series';
        for ista = 1: nsta
            yofx = tracebb(:,ista);
            [vard,resd,pred,contr,perdb]=periodfit(yofx,perdx,x,comp);
            plot(ax,x,pred,color(ista),'linew',0.5);
        end
        if lo<1 || hi <1
            lop = round(1/lo);
            hip = round(1/hi);
            text(ax,0.05,0.9,strcat(num2str(hip),'-',num2str(lop),{' s'}),'fontsize',9,'unit',...
                'normalized','Horizontalalignment','left');
        else
            text(ax,0.05,0.9,strcat(num2str(lop),'-',num2str(hip),{' Hz'}),'fontsize',9,'unit',...
                'normalized','Horizontalalignment','left');
        end
        text(ax,0.15,0.9,strcat({'Periodfit'}),'fontsize',9,'unit',...
             'normalized','Horizontalalignment','left');
        axis(ax,[is ien -ym ym]);
        set(ax,'XTick',is:4:ien,'fontsize',8);
        hold(ax,'off');
        
        
%         figure
%         plot(vard);
    end

% end    
    
    
    
    
    
    
    
    
    
    
    
    
    