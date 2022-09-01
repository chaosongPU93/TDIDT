% function mig_spectrogram_LZB_exp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to obtain the spectrogram of the seismogram during the
% time of the migrations in LZB trio. Read in the time of the migration and
% then read the corresponding data at three stations, then use the built-in 
% matlab function spectrogram to calculate
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/02/17
% Last modified date:   2021/02/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all
clc

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);

rstpath = strcat(getenv('MHOME'), '/Seisbasics/hypoinverse/lzbrst');

datapath = strcat(getenv('ALLAN'),'/data-no-resp');

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
%             '158';      % 158, 20200916,testing purpose
            ];     % 2020/09/07, i am adding a new family 006 to the pool, final version
        
nfam = size(nfampool,1);
disp(nfam); 

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
         'LZB'];
POLSTA=['SSIB '           % polaris station names
        'SILB '
        'KLNB '
        'MGCB '
        'TWKB '];
    
stas=['TWKB '
      'LZB  '
      'MGCB '];     % determine the trio and order         

  
% load files
winlenhf = 4;
% ccminhf = 0.4;
ccminhf = 0.01;

winlenlf = 16;
lofflf = 4;
ccminlf = 0.35;

SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'.',num2str(ccminhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname1 = strcat(rstpath, '/evtloc.allfam.exp.nodcutnodou.',SUFFIXhf);
fname2 = strcat(rstpath, '/evtloc.allfam.exp.eq8kmdcutnodou.',SUFFIXhf);
fname3 = strcat(rstpath, '/evtloc.allfam.exp.eq12kmdcutnodou.',SUFFIXhf);

fname = fname3;
hftime = load(fname);
% 25 cols, format is:
%   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date
%   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum
%

%%% NOTE 
SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname1 = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.nodcutnodou.',SUFFIXlf);
fname2 = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.eq8kmdcutnodou.',SUFFIXlf);
fname3 = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.eq12kmdcutnodou.',SUFFIXlf);

fname = fname3;
lftime = load(fname);

% this is inverted from (0,0) of all fams, same order, location of control points
loccont = [-123.492667 48.451500 38.1400; 
           -123.772167 48.493000 35.5900; 
           -123.863167 48.528167 35.2100;
           -123.603333 48.440167 36.7100;
           -123.800167 48.408833 34.5200;
           -123.893333 48.536500 35.0700;
           -123.864500 48.498667 34.8800;
           -123.753333 48.525667 36.2000;
           -123.703667 48.502667 36.4100;
           -123.814333 48.538667 35.7900;
           -123.838500 48.544833 35.6600;
           -123.908000 48.494167 34.5100;       % 006
           -123.879667 48.446167 34.2600;       % 001
%            -123.967000 48.458667 33.7800;       % 158, 20200916,testing purpose
           ];
       

relacont = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),loccont(2,1),loccont(2,2));
relacont(:,1) = dx;
relacont(:,2) = dy;

%%% load the location resolution points, +-1 & 2 samples
reslzb1 = load(fullfile(rstpath,'evtloc.13fam_locres_hf1spl'));
reslzb2 = load(fullfile(rstpath,'evtloc.13fam_locres_hf2spl'));

% sort according to day, sec
hftime = sortrows(hftime, [13, 15]);
lftime = sortrows(lftime, [13, 15]);

% this is new time ranges that come from 'identify_RTMs_v3.m'
trange = [
    2004198   8.3520e+04   8.6400e+04   % combine, y
    ];

% hardcopy of angbest from either min. rmse or max. slope of accepted ones ntranacc
angbest = [
   240
   ];


trange = [
    2004198     84222	86400   % #13
    2005256     76381   77940   % #38
	2005259     5232	7714    % #49   
    ];

% extend the time range by -+0.2 hr
trange1(:,1) = trange(:,1);
trange1(:,2) = trange(:,2)-0.2*3600;
trange1(:,3) = trange(:,3)+0.2*3600;
    
angbest = [
    245
    250
    255
    ];

%%

hfhi = 6.5;
hflo = 1.25;
lfhi = 1.25;
lflo = 0.5;
npa = 2;
npo = 2;

for i = 1: size(trange1,1)
    
    indhf = find(hftime(:,13)==trange1(i,1) & hftime(:,15)>=trange1(i,2) & ...
                 hftime(:,15)<=trange1(i,3));
    mighf = hftime(indhf,:);
    
    indlf = find(lftime(:,13)==trange1(i,1) & lftime(:,15)>=trange1(i,2) & ...
                 lftime(:,15)<=trange1(i,3));
    miglf = lftime(indlf,:); 
    
    % use the fam that most hf detections belong to
    % knowing the family is important bc it is related to rots para to read the data
    famhf = mode(mighf(:,end));
    
    if famhf < 10 
        fam = strcat('00',num2str(famhf));
    elseif famhf < 100
        fam = strcat('0',num2str(famhf));
    else
        fam = num2str(famhf);
    end
    
    % get permanent and polaris station rotation parameters
    sft2=0;     % centroid shift of station 2
    sft3=0;     % centroid shift of station 3
    FLAG = 'LZB';
    CATA = 'new';
    [PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,sft2,sft3);
    reftime = POLROTS(5,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
    PERMROTS(:,4) = PERMROTS(:,4)-reftime;    % to make sure that 1st station is 0
    POLROTS(:,4) = POLROTS(:,4)-reftime;
    
    nsta=size(stas,1);         %  number of stations
    PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;     % convert 2-3 columns to rad 
    POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;
    sps=40;
    
    %%% 'scalefact' scales templates; 'scaleseisms' scales seisms.  Strategy changes with family.
    if isequal(fam,'002')
        scaleseisms=[1.0 0.76 0.95];       % scaleseisms scales seismograms
    elseif isequal(fam,'043')
        scaleseisms=[1.0 1.0 1.0];
    else
        scaleseisms=[1.0 1.0 1.0];
    end

    date = num2str(trange1(i,1));
    YEAR = date(1:4);
    year = str2double(YEAR);
    JDAY = date(5:7);
    jday = str2double(JDAY);
    MO=day2month(jday,year);     % EXTERNAL function, day2month, get the month of one particular date
    fprintf('Migration %d, using rots para from fam %s \n',i,fam);
    direc=[datapath, '/', YEAR,'/',MO,'/'];     % directory name
    
    prename = [direc, YEAR,'.',JDAY,'.00.00.00.0000.CN'];
    
    if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
        POLSTA(3,:)='KELB ';
    else
        POLSTA(3,:)='KLNB ';  % remember to change it back
    end
    
    %Read the data; find glitches (many consecutive zeros) and flag those windows in STAnzeros.
    %Get timsSTA from the permanent stations (last one over-writes):
    lenofday = 24*3600;
    STAopt=zeros(nsta,sps * lenofday);
    STAort=zeros(nsta,sps * lenofday);
    Topt=zeros(nsta,sps * lenofday);
    
    STAopthf=zeros(nsta,sps * lenofday);
    STAorthf=zeros(nsta,sps * lenofday);
    
    STAoptlf=zeros(nsta,sps * lenofday);
    STAortlf=zeros(nsta,sps * lenofday);
    
    fileflag = 1;   % 1 means all files exist, the file status is normal
    for ista=1:nsta
        found=0;
        
        [LIA,idx]=ismember(stas(ista,:),PERMSTA,'rows');     % ismember, to determine whether each row of stas is contained in PERMSTA, return logical value 1/0 and index
        if LIA      % if station is a permanent one
            found=found+LIA;
            if strcmp(PERMSTA(idx,1:3),'PGC')     % string compare
                %%% WAHT is the meanning of 'fact', similar to instrument
                %%% response
                fact=1.0e-3;
            elseif strcmp(PERMSTA(idx,1:3),'LZB')
                fact=1.8e-3;
            end
            %%% readperms is an EXTERNAL FUNCTION
            % opt: optimal seismogram after rotations
            % ort: orthogonal seismogram after rotations
            % nzeros: number of zeros in the trace
            % timsSTA: time sequence
            fname = strcat(prename,'.',PERMSTA(idx,1:3),'..BHE.D.SAC');
            if isfile(fname)    % if have the data file
%                 % 1. this is for data without removing station response
%                 [opt,ort,nzeros,timsSTA]=readperms(prename,PERMSTA,PERMROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
                % 2. this is for data with no response
                [opt, ort, time] = readperm_nofilterv2(prename, ...
                                                       PERMSTA, PERMROTS, idx, sps, fact);                                                       
                [opthf,orthf,~] = readpermsv3(prename,PERMSTA,PERMROTS,idx,sps,hflo,hfhi,npo,...
                                                npa,fact);
                [optlf,ortlf,~] = readpermsv3(prename,PERMSTA,PERMROTS,idx,sps,lflo,lfhi,npo,...
                                                npa,fact);                            
            else
                fileflag = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for station %s in day %s %s, this day will be omitted. \n',...
                        PERMSTA(idx,1:3), YEAR, JDAY);
                break   % break the entire station loop
            end                
        end
        
        [LIA,idx]=ismember(stas(ista,:),POLSTA,'rows');    
        if LIA        % if are in POLSTA
            found=found+LIA; %better be 1
            if year==2003 && jday<213        % should result from some criteria
                fact=7.5e-3;
            else
                fact=1.5e-3; 
            end
            
            %%% readpols is an EXTERNAL FUNCTION
            % opt:
            % ort:
            % nzeros:
            fname = strcat(prename,'.',POLSTA(idx,1:4),'..HHE.D.SAC');
            if isfile(fname)    % if have the data file
%                 % 1. this is for data without removing station response
%                 [opt,ort,nzeros]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
                % 2. this is for data with no response
                [opt, ort, time] = readpola_nofilterv2(prename, ...
                                                       POLSTA, POLROTS, idx, sps, fact);
                [opthf,orthf]=readpolsv3(prename,POLSTA,POLROTS,idx,sps,hflo,hfhi,npo,npa,fact);
                [optlf,ortlf]=readpolsv3(prename,POLSTA,POLROTS,idx,sps,lflo,lfhi,npo,npa,fact);                       
            else
                fileflag = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for station %s in day %s / %s. \n',POLSTA(idx,1:4), YEAR, JDAY);
                break   % break the entire station loop
            end                  
        end
%         found=found         % could be a benchmark
        %factr1(ista)=prctile(abs(opt),90); %Not normalized
        %factr2(ista)=factr1(ista)/factr1(1); %This is what is used; keeps 1st station unchanged but scales the others
                
        STAopt(ista,:)=opt./scaleseisms(ista);   
        STAort(ista,:)=ort./scaleseisms(ista);
        Topt(ista,:) = time;
        
        STAopthf(ista,:)=opthf./scaleseisms(ista);   
        STAorthf(ista,:)=orthf./scaleseisms(ista);
        STAoptlf(ista,:)=optlf./scaleseisms(ista);   
        STAortlf(ista,:)=ortlf./scaleseisms(ista);
        
    end

    if fileflag == 0    % means there are missing files
        fprintf('Day %s / %s will be omitted because of missing files. \n', YEAR, JDAY);
        continue    % continue to the next day
    end
    
    sttime = max(trange1(i,2), 0);
    edtime = min(trange1(i,3), 86400);
    
    sttimelong = max(sttime-0.2*3600,0);
    edtimelong = min(edtime+0.2*3600,86400);
    stindlong = find(abs(sttimelong-Topt(1,:))<= 1/sps,1,'first');    % origin time index
    edindlong = find(abs(edtimelong-Topt(1,:))<= 1/sps,1,'first');    % origin time index
    STAopthflong = STAopthf(:, stindlong+1: edindlong);
    STAoptlflong = STAoptlf(:, stindlong+1: edindlong);
    
    
    stind = find(abs(sttime-Topt(1,:))<= 1/sps,1,'first');    % origin time index
    edind = find(abs(edtime-Topt(1,:))<= 1/sps,1,'first');    % origin time index
    Topt = Topt(:, stind+1: edind);
    STAopt = STAopt(:, stind+1: edind);
    STAort = STAort(:, stind+1: edind);
    
    STAopthf = STAopthf(:, stind+1: edind);
    STAorthf = STAorthf(:, stind+1: edind);
    STAoptlf = STAoptlf(:, stind+1: edind);
    STAortlf = STAortlf(:, stind+1: edind);
    
    juldate = num2str(trange1(i,1));
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

%%    
    f1.fig = figure;
%     f1.fig.Renderer='Painters';
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 9;   % maximum height allowed is 11 inches
    set(f1.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
%     nrow = 3;
    nrow = 9;
    ncol = 1;    
    for isub = 1:nrow*ncol
        f1.ax(isub) = subplot(nrow,ncol,isub);
    end
    
    %%% reposition
    set(f1.ax(3), 'position', [ 0.08, 0.12, 0.75, 0.22]);
    set(f1.ax(2), 'position', [ 0.08, 0.44, 0.75, 0.22]);
    set(f1.ax(1), 'position', [ 0.08, 0.76, 0.75, 0.22]);
    
    set(f1.ax(6), 'position', [ 0.2, 0.27, 0.5, 0.06]);
    set(f1.ax(5), 'position', [ 0.2, 0.59, 0.5, 0.06]);
    set(f1.ax(4), 'position', [ 0.2, 0.91, 0.5, 0.06]);
    
    set(f1.ax(9), 'position', [ 0.08, 0.07, 0.75, 0.04]);
    set(f1.ax(8), 'position', [ 0.08, 0.39, 0.75, 0.04]);
    set(f1.ax(7), 'position', [ 0.08, 0.71, 0.75, 0.04]);
    
    
    f2.fig = figure;
%     f1.fig.Renderer='Painters';
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 9;   % maximum height allowed is 11 inches
    set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
%     nrow = 3;
    nrow = 9;
    ncol = 1;    
    for isub = 1:nrow*ncol
        f2.ax(isub) = subplot(nrow,ncol,isub);
    end
    
%     %%% reposition
%     set(f2.ax(3), 'position', [ 0.08, 0.12, 0.8, 0.09]);
%     set(f2.ax(2), 'position', [ 0.08, 0.44, 0.8, 0.09]);
%     set(f2.ax(1), 'position', [ 0.08, 0.76, 0.8, 0.09]);
%     
%     set(f2.ax(6), 'position', [ 0.08, 0.27, 0.8, 0.09]);
%     set(f2.ax(5), 'position', [ 0.08, 0.59, 0.8, 0.09]);
%     set(f2.ax(4), 'position', [ 0.08, 0.91, 0.8, 0.09]);
%     
%     set(f2.ax(9), 'position', [ 0.08, 0.07, 0.8, 0.09]);
%     set(f2.ax(8), 'position', [ 0.08, 0.39, 0.8, 0.09]);
%     set(f2.ax(7), 'position', [ 0.08, 0.71, 0.8, 0.09]);
    
    for ista = 1: nsta
        tmpopt = STAopt(ista, :);
        % each overlapping window length
%         seglensec = 16;
%         seglen = seglensec*sps;
        seglen = 512;
        window = hann(seglen);
        nfft = pow2(nextpow2(seglen));
        Fs = sps;
        nov = seglen*0.75; 

%         %%%%%%%% the following is using the default plot settings %%%%%%%%
%         %%% from the benchmarking, we know that it is plotting the time t,
%         %%% frequency f, power spectral density PSD using pwelch method
%         %%% the default ploting function is 'surf', and it seems that the 
%         %%% colormap cannot be changed although i don't know why, something
%         %%% weird would happen if you want to change caxis as well.
%         f1 = figure;
%         spectrogram(tmpopt,window,nov,nfft,Fs,'yaxis');
%         ax = gca;
%         ax.YScale = 'log';
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         %%%%%%%%%%  this following codes could reproduce the default plot setting
%         %%% NOTE that 'spectrogram' internally uses PWELCH to compute the PSD
%         %%% estimate for each segment, which might not be recommended.
%         [spec, freq, t, psd, fc, tc] = spectrogram(tmpopt,window,nov,nfft,Fs);
% 
%         ax = f1.ax(ista);
%         ax.Box = 'on';
%         hold(ax, 'off');
% %         imagesc(ax,t, f, 10*log10(psd));
%         thr = (t+trange1(i,2))/3600;
%         surf(ax, thr,freq,10*log10(psd),'EdgeColor','none');
% %         colorlim = get(ax,'clim');
% %         newlim = [(colorlim(1)*0.8),colorlim(2)];
% %         set(ax,'clim',(newlim));
%             
%         ylim(ax,[min(freq) max(freq)]);
%         xlim(ax,[min(thr) max(thr)]);
%         ax.YScale = 'log';
%         xlabel(ax, strcat({'Time (hr) on '}, day, mo, yr));
%         ylabel(ax, 'Frequency (Hz)');
% %         colormap(ax,'jet');
%         c=colorbar(ax);
%         c.Label.String = 'Power/frequency (dB/Hz)';
% %         caxis(ax,[-160,0]);
%         view(ax,0,90);
%         ax.TickDir = 'out';
%         text(ax,0.8,0.9,stas(ista,:),'FontSize',12,'unit','normalized');
%         text(ax,0.1,0.9,num2str(i),'FontSize',12,'unit','normalized');
%         hold(ax, 'off');

        %%% Recommended this way %%%%%%%%%%%%       
        %%%%%%%%%%%%% the following is to test the function written by FJS
        %%% Well, it seems that 'spectrogram2' basically gives the similar result
        %%% with the same parameter setting.
        [psd, freq, t, psddB]=spectrogram2(tmpopt,nfft,Fs,seglen,nov,'s');
        
        %%% now get the sqaured amplitude sum ratio of LF to HF;
        tmpopthf = STAopthf(ista, :);
        tmpoptlf = STAoptlf(ista, :);
        enghf = zeros(size(t));
        englf = zeros(size(t));
        engrat = zeros(size(t));
        for j = 1: length(t)
            indst = t(j)*sps+1;
            inded = t(j)*sps+seglen;
            tmp = tmpopthf(indst: inded);
            enghf(j) = sum(tmp.*tmp);
            
            tmp = tmpoptlf(indst: inded);
            englf(j) = sum(tmp.*tmp);
            
            engrat(j) = englf(j)./enghf(j);
        end
        meanrat = mean(engrat);
        
        tmpopthflong = STAopthflong(ista, :);
        tmpoptlflong = STAoptlflong(ista, :);
        len = floor((edtimelong-sttimelong)/3.2)-3;
        
        tlong = t(1): 3.2: t(1)+3.2*(len-1);
        tlong = tlong';
        enghflong = zeros(size(tlong));
        englflong = zeros(size(tlong));
        engratlong = zeros(size(tlong));
        for j = 1: length(tlong)
            indst = floor(tlong(j)*sps+1);
            inded = floor(tlong(j)*sps+seglen);
            tmp = tmpopthflong(indst: inded);
            enghflong(j) = sum(tmp.*tmp);
            
            tmp = tmpoptlflong(indst: inded);
            englflong(j) = sum(tmp.*tmp);
            
            engratlong(j) = englflong(j)./enghflong(j);
        end
        
        ax = f1.ax(ista);
        ax.Box = 'on';
        hold(ax, 'on');
        
        yyaxis(ax, 'left');
%         imagesc(ax,t, f, 10*log10(psd));
        thr = (t+trange1(i,2))/3600;
        surf(ax, thr,freq,10*log10(psd),'EdgeColor','none');
%         colorlim = get(ax,'clim');
%         newlim = [(colorlim(1)*0.8),colorlim(2)];
%         set(ax,'clim',(newlim));
        plot(ax, [thr(1), thr(end)], [hfhi, hfhi], 'k--');
        plot(ax, [thr(1), thr(end)], [hflo, hflo], 'b--');
        plot(ax, [thr(1), thr(end)], [lflo, lflo], 'k--');
        
%         % mark the HF and LF detections on the plot as well
%         xhf = [mighf(:,15)/3600' mighf(:,15)/3600'];
%         yhf = [hfhi*ones(size(mighf(:,15))) max(freq)*ones(size(mighf(:,15)))];
%         for j = 1: size(mighf(:,15),1)
%             plot(ax, xhf(j,:), yhf(j,:), 'w -');
%         end
%         
%         xlf = [miglf(:,15)/3600' miglf(:,15)/3600'];
%         ylf = [0.08*ones(size(miglf(:,15))) lflo*ones(size(miglf(:,15)))];
%         for j = 1: size(miglf(:,15),1)
%             plot(ax, xlf(j,:), ylf(j,:), 'w-');
%         end
        
        ylim(ax,[min(freq) max(freq)]);
        xlim(ax,[min(thr) max(thr)]);
        ax.YScale = 'log';
%         xlabel(ax, strcat({'Time (hr) on '}, day, mo, yr));
        ylabel(ax, 'Frequency (Hz)','fontsize',10);
%         colormap(ax,'jet');
        c=colorbar(ax,'East');
        pos = ax.Position;
        c.Position = [pos(1)+pos(3)+0.08, pos(2), 0.02, pos(4)];
        c.Label.String = 'Power/frequency (dB/Hz)';
%         caxis(ax,[-160,0]);
        view(ax,0,90);
        ax.TickDir = 'out';
        text(ax,0.9,0.9,stas(ista,:),'FontSize',12,'unit','normalized');
        text(ax,0.05,0.9,num2str(i),'FontSize',12,'unit','normalized');
        

        yyaxis(ax, 'right');
        % plot the ratio at those detections only 
        hft = mighf(:,15)- 4/2;
        lft = miglf(:,15)- 16/2;
%         hft = mighf(:,14);
%         lft = miglf(:,14);
        hfengrat = zeros(length(hft),1);
        for jj = 1: length(hft)
            [~,ind] = min(abs(t+trange1(i,2)-hft(jj)));
            hfengrat(jj) = engrat(ind);
        end
        lfengrat = zeros(length(lft),1);
        for jj = 1: length(lft)
            [~,ind] = min(abs(t+trange1(i,2)-lft(jj)));
            lfengrat(jj) = engrat(ind);
        end
        
        lfthr = (lft)/3600;
        plot(ax, lfthr,lfengrat,'.-','color','w','linew',1);
        
        hfthr = (hft)/3600;
%         plot(ax, hfthr,hfengrat+max(lfengrat),'.-','color','w','linew',1);        
        plot(ax, hfthr,hfengrat, '.-','color',[0.5 0.5 0.5],'linew',1);
        ylabel(ax, 'Energy ratio (LF/HF)','fontsize',10);
        
        ax.XTickLabel = [];
        hold(ax, 'off');
        
        % plot energy ratio of a longer window
        ax = f1.ax(ista+3);
        ax.Box = 'off';
        hold(ax, 'on');
        ax.Color = 'none';
        thrlong = (tlong+sttimelong)/3600;
        plot(ax, thrlong,engratlong, 'k-','linew',1);
        plot(ax, [thr(1) thr(1)], [ax.YLim(1) ax.YLim(2)], 'r-','linew',2);
        plot(ax, [thr(end) thr(end)], [ax.YLim(1) ax.YLim(2)], 'r-','linew',2);
        plot(ax, [trange(i,2)/3600 trange(i,2)/3600], [ax.YLim(1) ax.YLim(2)], 'b-','linew',2);
        plot(ax, [trange(i,3)/3600 trange(i,3)/3600], [ax.YLim(1) ax.YLim(2)], 'b-','linew',2);

        xlim(ax,[min(thrlong) max(thrlong)]);

        hold(ax, 'off');
        
        % plot seismograms
        ax = f1.ax(ista+6);
        ax.Box = 'on';
        hold(ax, 'on');
        plot(ax, Topt(ista,:)/3600,STAopt(ista,:), 'k-','linew',1);
        plot(ax, [trange(i,2)/3600 trange(i,2)/3600], [ax.YLim(1) ax.YLim(2)], 'b-','linew',2);
        plot(ax, [trange(i,3)/3600 trange(i,3)/3600], [ax.YLim(1) ax.YLim(2)], 'b-','linew',2);

        xlim(ax,[min(thr) max(thr)]);
        hold(ax, 'off');
        
        
        %%% figure 2
        indhf = find(hftime(:,13)==trange1(i,1) & hftime(:,15)>=sttimelong & ...
                 hftime(:,15)<=edtimelong);
        mighf = hftime(indhf,:);
    
        indlf = find(lftime(:,13)==trange1(i,1) & lftime(:,15)>=sttimelong & ...
                 lftime(:,15)<=edtimelong);
        miglf = lftime(indlf,:); 
    
        % plot energy ratio of a longer window
        ax = f2.ax((ista-1)*3+1);
        ax.Box = 'on';
        hold(ax, 'on');
        ax.Color = 'none';
        thrlong = (tlong+sttimelong)/3600;
        plot(ax, thrlong,engratlong, 'k-'); %,'linew',1
        plot(ax, [thr(1) thr(1)], [ax.YLim(1) ax.YLim(2)], 'r-','linew',2);
        plot(ax, [thr(end) thr(end)], [ax.YLim(1) ax.YLim(2)], 'r-','linew',2);
        plot(ax, [trange(i,2)/3600 trange(i,2)/3600], [ax.YLim(1) ax.YLim(2)], 'b-','linew',2);
        plot(ax, [trange(i,3)/3600 trange(i,3)/3600], [ax.YLim(1) ax.YLim(2)], 'b-','linew',2);
        text(ax,0.9,0.9,stas(ista,:),'FontSize',8,'unit','normalized');
        text(ax,0.05,0.9,num2str(i),'FontSize',8,'unit','normalized');
%         text(ax,0.5,0.9,'Energy ratio (LF/HF)','FontSize',8,'unit','normalized');
        ylabel(ax,'Energy ratio (LF/HF)','FontSize',8);
        xlim(ax,[min(thrlong) max(thrlong)]);
        plot(ax, ax.XLim, [1 1], '--','linew',1,'color',[.5 .5 .5]);        
        hold(ax, 'off');

        % plot lf energy ratio of a longer window
        ax = f2.ax((ista-1)*3+2);
        ax.Box = 'on';
        hold(ax, 'on');
        ax.Color = 'none';
        thrlong = (tlong+sttimelong)/3600;
        plot(ax, thrlong,englflong, 'k-');
        plot(ax, [thr(1) thr(1)], [ax.YLim(1) ax.YLim(2)], 'r-','linew',2);
        plot(ax, [thr(end) thr(end)], [ax.YLim(1) ax.YLim(2)], 'r-','linew',2);
        plot(ax, [trange(i,2)/3600 trange(i,2)/3600], [ax.YLim(1) ax.YLim(2)], 'b-','linew',2);
        plot(ax, [trange(i,3)/3600 trange(i,3)/3600], [ax.YLim(1) ax.YLim(2)], 'b-','linew',2);
        scatter(ax, miglf(:,14)/3600,ax.YLim(2)*2/3*ones(1,length(miglf(:,14))), 6, 'o');
%         text(ax,0.9,0.9,stas(ista,:),'FontSize',8,'unit','normalized');
%         text(ax,0.05,0.9,num2str(i),'FontSize',8,'unit','normalized');
%         text(ax,0.5,0.9,'Energy LF','FontSize',8,'unit','normalized');
        ylabel(ax,'Energy LF','FontSize',8);
        xlim(ax,[min(thrlong) max(thrlong)]);
        hold(ax, 'off');
        
        % plot hf energy ratio of a longer window
        ax = f2.ax((ista-1)*3+3);
        ax.Box = 'on';
        hold(ax, 'on');
        ax.Color = 'none';
        thrlong = (tlong+sttimelong)/3600;
        plot(ax, thrlong,enghflong, 'k-');
        ax.YLim = f2.ax((ista-1)*3+2).YLim*1;
        plot(ax, [thr(1) thr(1)], [ax.YLim(1) ax.YLim(2)], 'r-','linew',2);
        plot(ax, [thr(end) thr(end)], [ax.YLim(1) ax.YLim(2)], 'r-','linew',2);
        plot(ax, [trange(i,2)/3600 trange(i,2)/3600], [ax.YLim(1) ax.YLim(2)], 'b-','linew',2);
        plot(ax, [trange(i,3)/3600 trange(i,3)/3600], [ax.YLim(1) ax.YLim(2)], 'b-','linew',2);
        scatter(ax, mighf(:,14)/3600,ax.YLim(2)*2/3*ones(1,length(mighf(:,14))), 6, '^');
%         text(ax,0.9,0.9,stas(ista,:),'FontSize',8,'unit','normalized');
%         text(ax,0.05,0.9,num2str(i),'FontSize',8,'unit','normalized');
%         text(ax,0.5,0.9,'Energy HF','FontSize',8,'unit','normalized');
        ylabel(ax,'Energy HF','FontSize',8);
        xlim(ax,[min(thrlong) max(thrlong)]);
        hold(ax, 'off');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
    
    
    xlabel(f1.ax(9),strcat({'Time (hr) on '}, day, mo, yr),'fontsize',10);
    title(f1.ax(1),strcat({'Mig. '},num2str(i),{': '},num2str(trange1(i,1)),{' '},...
          num2str(trange1(i,2)),{'-'},num2str(trange1(i,3)),'s'));
      
    % save figure
    print(f1.fig,'-dpdf',strcat(rstpath,'/LZB.mig.autortm.spectrogram',num2str(trange1(i,1)),...
          '_',num2str(trange1(i,2)),'-',num2str(trange1(i,3)),'.',num2str(winlenhf),'_',...
          num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));
    
      
    xlabel(f2.ax(9),strcat({'Time (hr) on '}, day, mo, yr),'fontsize',10);        
      
    print(f2.fig,'-dpdf',strcat(rstpath,'/LZB.mig.autortm.energyratio',num2str(trange1(i,1)),...
          '_',num2str(trange1(i,2)),'-',num2str(trange1(i,3)),'.',num2str(winlenhf),'_',...
          num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));  
end
       
       
          
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       