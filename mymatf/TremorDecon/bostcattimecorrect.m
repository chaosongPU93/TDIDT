% bostcattimecorrect.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It has been noticed that the time of the LFE catalog from Michael Bostock
% does not necessarily point to the peak or zero-crossing of the LFEs. This
% is partially solved when I tried to create LFE templates for several fams
% using TWKB/PGC trio in 'mk_bbtemp_PGC.m' and 'mk_bbtemp_LZB.m' by 
% introducing 'offsec'. This script tries to systematically estimate this 
% time correction for all families.
%
% --Instead of creating templates for each possible fam at different stations,
%   in this scriptm we are only interested in the offset time between the true 
%   peak/zero-crossing and catalog time at station PGC.
% --So, we don't have to compute the rotation parameters or anything. Just 
%   simply stack the N and E component data based on the catalog timing, and 
%   see how far away the peak/zero-crossing is from the start. Use the median
%   time of N and E.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/03/30
% Last modified date:   2022/03/30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio
format short e   % Set the format to 5-digit floating point
clear
clc
close all

pltflag = 1;
if pltflag == 1
  set(0,'DefaultFigureVisible','on');
else
  set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
end

[scrsz, res] = pixelperinch(1);

% set path
workpath = getenv('ALLAN');

% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');

sps=40;     % samples per second

%locations of the origin, say it's the offset (0,0)
loc0 = off2space002([0 0],sps,'interpArmb',0); % 8 cols, format: dx,dy,lon,lat,dep,tori,off12,off13
%get the location of lfe fams that have detections
lfeloc = LfelocBostock(loc0(3),loc0(4)); % format: [fam dx dy lon lat dep], 6 cols
%note that fam have detections and fam have locations are NOT the same
famall = unique(lfeloc(:,1));
nfam = size(famall,1);
for i = 1: nfam
  if famall(i)<10
    fampool(i,:) = strcat('00',num2str(famall(i)));
  elseif famall(i)<100
    fampool(i,:) = strcat('0',num2str(famall(i)));
  else
    fampool(i,:) = num2str(famall(i));
  end  
end

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
    'LZB'];
POLSTA=['SSIB '           % polaris station names
    'SILB '
    'KLNB '
    'MGCB '
    'TWKB '];

stas=['PGC  '
     ];
% number of used stations
nsta=size(stas,1);         %  number of stations

% get the all catalog LFEs in that family
bostname = ('/BOSTOCK/total_mag_detect_0000_cull_NEW.txt');
catalog = load(strcat(workpath, bostname));

%%%This is posterior knowledge after checking every family
%there are a few suspious fams that need to be checked carefully again
susfam = ['006';
          '070';
          '010';
          '045';
          '074';
          '132';
          '231';
          '232';
          '233';
          ];
[~,~,susind] = intersect(susfam,fampool,'rows','stable');       

remake = 0;

%% stack at bostock's detection timings on N and E components 
if remake
  for i = 1: nfam
    % for i = 1:30
    % for i = 1:length(susind)
    %   if ~ismember(i,susind)
    %     continue
    %   end
    fam = fampool(i,:)
    %   fam = fampool(susind(i),:)
    %   i = find(famall==240);
    
    [timoffrot,~] = GetDays4Stack(fam);
    nday = size(timoffrot, 1);
    fprintf('%d days \n', nday);
    
    if nday == 0
      disp("No day with detections found in this family of BOSTOCK's catalog");
      return
    end
    
    %%% data parameters
    lenofday = 24*3600;
    STAE=zeros(sps*lenofday, nsta);
    STAN=zeros(sps*lenofday, nsta);
    %%% desired template parameters
    templensec = 50;
    templen = templensec * sps;
    extra = 2*sps;    % use extra more samples than the original window.
    stackexE = zeros(templen+ 2*extra,nsta);
    stackexN = zeros(templen+ 2*extra,nsta);
    
    catfam = catalog(famall(i) == catalog(:, 1), :);
    
    nstack = 0;
    nLFE = 0;
    
    %%% loop for day
    for id = 1: nday
      %     fprintf('%d / %d \n',id, nday);
      
      year = timoffrot(id,1);
      YEAR = int2str(year);
      jday = timoffrot(id,2);
      if year == 2003
        date = jday-62+30303;
      elseif year == 2004
        date = jday-196+40714;
      elseif year == 2005
        date = jday-254+50911;
      end
      bostocks = catfam(date == catfam(:, 2), :);
      bostsec = 3600*(bostocks(:,3)-1)+bostocks(:,4);
      bostsamp = round(bostsec*sps);
      nlfe1day = size(bostsamp, 1);
      nLFE = nLFE + nlfe1day;
      
      if jday <= 9
        JDAY=['00',int2str(jday)];
      elseif jday<= 99
        JDAY=['0',int2str(jday)];
      else
        JDAY=int2str(jday);
      end
      %     MO = 'SEP';
      MO=day2month(jday,year);     % EXTERNAL function, day2month, get the month of one particular date
      direc=[datapath, '/arch', YEAR,'/',MO,'/'];     % directory name
      %     fprintf('%s \n',direc);
      datafnm = [direc, YEAR,'.',JDAY,'.00.00.00.0000.CN'];
      %     fprintf('Processing date: %s / %s \n',YEAR, JDAY);
      
      if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
        POLSTA(3,:)='KELB ';
      else
        POLSTA(3,:)='KLNB ';  % remember to change it back
      end
      
      fileflag = 1;   % 1 means all files exist, the file status is normal
      
      %%% loop for each station for reading data
      for ista = 1: nsta
        found = 0;
        [LIA, idx] = ismember(stas(ista, :), PERMSTA, 'rows');
        
        %%% if station is a permanent one
        if LIA
          found = found+ LIA;
          if strcmp(PERMSTA(idx, 1:3),'PGC')
            fact = 1.0e-3;
          elseif strcmp(PERMSTA(idx, 1:3),'LZB')
            fact = 1.8e-3;
          end
          fname = strcat(datafnm,'.',PERMSTA(idx,1:3),'..BHE.D.SAC');
          if isfile(fname)    % if have the data file
            % 1. this is for data without removing station response
            %   [opt,ort,timeperm]=readperm_1daynofilterv2(datafnm, PERMSTA, PERMROTS, idx, sps, fact);
            % 2. this is for data with no response
            [dataE, dataN, ~] = readperms_norotsv2(datafnm, PERMSTA, idx, sps, fact);
          else
            fileflag = 0;   % change the file flag to 0, meaning abnormal
            fprintf('No data for station %s in day %s %s, this day will be omitted. \n',...
              PERMSTA(idx,1:3), YEAR, JDAY);
            break   % break the entire station loop
          end
          
        end
        
        %%% if station is a polaris one
        [LIA, idx] = ismember(stas(ista, :), POLSTA, 'rows');
        if LIA
          found = found+ LIA;
          if year == 2003 && jday < 213
            fact = 7.5e-3;
          else
            fact = 1.5e-3;
          end
          fname = strcat(datafnm,'.',POLSTA(idx,1:4),'..HHE.D.SAC');
          if isfile(fname)    % if have the data file
            %                 % 1. this is for data without removing station response
            %                 [opt,ort,nzeros]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
            % 2. this is for data with no response
            [dataE, dataN, ~] = readpols_norotsv2(datafnm, POLSTA, idx, sps, fact);
          else
            fileflag = 0;   % change the file flag to 0, meaning abnormal
            fprintf('No data for station %s in day %s / %s. \n',POLSTA(idx,1:4), YEAR, JDAY);
            break   % break the entire station loop
          end
          
        end
        
        STAE(:, ista) = dataE;
        STAN(:, ista) = dataN;
      end
      
      if fileflag == 0    % means there are missing files
        fprintf('Day %s / %s will be omitted because of missing files. \n', YEAR, JDAY);
        continue    % continue to the next day
      end
      
      len = size(STAE,1);
      for n=1: nlfe1day
        % set the timing be at the center of the window
        if (bostsamp(n) +templen -1 +extra <= len) && ...
            (bostsamp(n) -extra >= 1)
          
          %%% it is better to use a longer window than templen to make
          %%% sure the shifting in step2 would preserve the frequency
          %%% content as much as possible
          windata = STAE(bostsamp(n) -extra: bostsamp(n) +templen -1 +extra, :);
          % remove mean and linear trend
          tmp = detrend(windata);
          stackexE = stackexE+ tmp;
          
          windata2 = STAN(bostsamp(n) -extra: bostsamp(n) +templen -1 +extra, :);
          tmp2 = detrend(windata2);
          stackexN = stackexN+ tmp2;
          
          nstack = nstack + 1;
          
          % datamat is 3D, for storing all qualified wins for stacking
          % size is (templen, nsta, nstack)
          dmatE(:, :, nstack) = STAE(bostsamp(n) -extra: bostsamp(n) +templen -1 +extra, :);
          dmatN(:, :, nstack) = STAN(bostsamp(n) -extra: bostsamp(n) +templen -1 +extra, :);
          
        end
      end
      
    end     % loop for days
    
    %detrend and average by the number of stacked LFEs
    stackexE = detrend(stackexE/nstack);
    stackexN = detrend(stackexN/nstack);
    
    %chop the extra
    stackE = detrend(stackexE(1+extra: end-extra, :));
    stackN = detrend(stackexN(1+extra: end-extra, :));
    
    %filter the trace with buffering
    lo = 0.5;
    hi = 6.5;
    npo = 2;
    npa = 2;
    for ista =1: nsta
      stackexEf(:,ista) = Bandpass(stackexE(:,ista), sps, lo, hi, npo, npa, 'butter');
      stackexNf(:,ista) = Bandpass(stackexN(:,ista), sps, lo, hi, npo, npa, 'butter');
    end
    stackEf =detrend(stackexEf(1+extra: end-extra, :));
    stackNf =detrend(stackexNf(1+extra: end-extra, :));
    
    
    %plot the stack
    figure  % plot original stack, fig 1
    ax = subplot(2,1,1);  % plot E component
    hold(ax,'on');
    for ista = 1:nsta
      plot(ax,(1:templen)/sps, stackE(:,ista)*10 + 1* ista,'b');
      plot(ax,(1:templen)/sps, stackN(:,ista)*10 - 1* ista,'r');
      text(ax,50/sps,1*ista,strcat(stas(ista,:),{' E'}));
      text(ax,50/sps,-1*ista,strcat(stas(ista,:),{' N'}));
      text(ax,0.05,0.9,fam,'Units','normalized','HorizontalAlignment','left','FontSize',10);
      text(ax,0.95,0.9,'Unfiltered','Units','normalized','HorizontalAlignment','right','FontSize',10);
    end
    ax.Box = 'on';
    hold(ax,'off');
    
    ax = subplot(2,1,2);  % plot N component
    hold(ax,'on');
    for ista = 1:nsta
      plot(ax,(1:templen)/sps, stackEf(:,ista)*10 + 1* ista,'b');
      plot(ax,(1:templen)/sps, stackNf(:,ista)*10 - 1* ista,'r');
      text(ax,50/sps,1*ista,strcat(stas(ista,:),{' E'}));
      text(ax,50/sps,-1*ista,strcat(stas(ista,:),{' N'}));
      text(ax,0.05,0.9,fam,'Units','normalized','HorizontalAlignment','left','FontSize',10);
      text(ax,0.95,0.9,sprintf('BP %.2f-%.2f Hz',lo,hi),'Units','normalized','HorizontalAlignment','right',...
        'FontSize',10);
    end
    
    %   vertical_cursors;
    
    %looks like the filtered traces have a better signal noise ratio to identify the zero-crossing
    %time
    ista = 1;
    [maxe,imaxe] = max(stackEf(:,ista));
    [mine,imine] = min(stackEf(:,ista));
    min2maxe = imaxe-imine; % min to max time in samples
    if min2maxe < 0
      disp('The polarity at E is reversed, needs check');
    end
    if abs(min2maxe) > 2*sps % usually this shouldn't be too big, unless SNR is low and arrival is ambiguious
      disp('The dipole at E might be unreliable, needs check');
    end
    
    %target the short segment between min and max
    [~, zcsta1e] = min(abs(stackEf(min(imine,imaxe): max(imine,imaxe), ista)));
    zcsta1e = zcsta1e+min(imine,imaxe)-1;
    
    %%%similarly to N component
    [maxn,imaxn] = max(stackNf(:,ista));
    [minn,iminn] = min(stackNf(:,ista));
    min2maxn = imaxn-iminn; % min to max time in samples
    if min2maxn < 0
      disp('The polarity at N is reversed, needs check');
    end
    if abs(min2maxn) > 2*sps % usually this shouldn't be too big, unless SNR is low and arrival is ambiguious
      disp('The dipole at N might be unreliable, needs check');
    end
    
    %target the short segment between min and max
    [~, zcsta1n] = min(abs(stackNf(min(iminn,imaxn): max(iminn,imaxn), ista)));
    zcsta1n = zcsta1n+min(iminn,imaxn)-1;
    
    zcsta1 = [zcsta1e, zcsta1n];
    %The zero-crossing is an average if both components are OK
    if abs(zcsta1e-zcsta1n) < 0.375*sps
      zcsta1 = round((zcsta1e+zcsta1n)/2)/sps;
    else
      %use the one which has a higher amplitude of the dipole
      [~,ind] = max([maxe-mine, maxn-minn]);
      zcsta1 = round(zcsta1(ind))/sps;
    end
    zc(i,1) = zcsta1;
    
    %also obtain the location of the positive peak
    pospk1 = [imaxe, imaxn];
    %The positive peak is an average if both components are OK
    if abs(imaxe-imaxn) < 0.375*sps
      pospk1 = round((imaxe+imaxn)/2)/sps;
    else
      %use the one which has a higher amplitude of the dipole
      [~,ind] = max([maxe-mine, maxn-minn]);
      pospk1 = round(pospk1(ind))/sps;
    end
    ppk(i,1) = pospk1;
    
    %absolute peak
    [maxae,imaxae] = max(abs(stackEf(:,ista)));
    [maxan,imaxan] = max(abs(stackNf(:,ista)));
    abspk1 = [imaxae, imaxan];
    %The absolute peak is an average if both components are OK
    if abs(imaxae-imaxan) < 0.375*sps
      abspk1 = round((imaxae+imaxan)/2)/sps;
    else
      %use the one which has a higher amplitude of the dipole
      [~,ind] = max([maxae, maxan]);
      abspk1 = round(abspk1(ind))/sps;
    end
    apk(i,1) = abspk1;
    
    plot(ax,[zcsta1e/sps zcsta1e/sps],ax.YLim,'b--');
    plot(ax,[zcsta1n/sps zcsta1n/sps],ax.YLim,'r--');
    plot(ax,[zcsta1 zcsta1],ax.YLim,'k--');
    plot(ax,[pospk1 pospk1],ax.YLim,'k:');
    plot(ax,[abspk1 abspk1],ax.YLim,'k-.');
    ax.Box = 'on';
    hold(ax,'off');
    
    
    %   keyboard
  end
  
  %%% manual corrected some families
  zc(5) = 25.62;
  ppk(5) = 25.75;
  apk(5) = 25.5;
  
  ppk(51) = 25.38;
  
  zc(68) = 31.02;
  ppk(68) = 31.15;
  apk(68) = 30.85;
  
  zc(69) = 32.02;
  ppk(69) = 32.15;
  apk(69) = 31.85;

  %save to file
  savemat = [famall zc ppk apk];
  fid = fopen(strcat(workpath,'/BOSTOCK/bostcattime_correction'),'w+');
  fprintf(fid,'%d %.4f %.4f %.4f \n', savemat');
  fclose(fid);

%Do NOT remake it, read from file instead
else
  timecor = load(strcat(workpath,'/BOSTOCK/bostcattime_correction'));
  
end

%% Make a summary plot of the time correction in map
%recall that we have the locations of the same LFE families
%'lfeloc': format: [fam dx dy lon lat dep], 6 cols

%options to colorcode the fams with time correction based on different targets
% timecol = 2;  %using zero-crossing
% timecol = 3;  %using max postive peak
timecol = 4;  %using max absolute peak

nrow = 1; % rows and cols of subplots in each figure
ncol = 1; 
widin = 6; % size of each figure
htin = 6;
f = initfig(widin,htin,nrow,ncol);
ax = f.ax(1);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
scatter(ax,lfeloc(:,2),lfeloc(:,3),20,timecor(:,timecol),'filled'); %,'MarkerEdgeColor','k'
colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
if timecol == 2
  corstr = 'zero-crossing';
elseif timecol == 3
  corstr = 'max postive peak';
elseif timecol == 4
  corstr = 'max absolute peak';
end
c.Label.String = strcat({'Time (s) correction based on '},corstr);
c.Label.FontSize = 8;
text(ax,lfeloc(:,2)+1,lfeloc(:,3)+1,num2str(lfeloc(:,1)),'FontSize',6);
xran = [-70 60];
yran = [-50 70];
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xran(1):10:xran(2));
yticks(ax,yran(1):10:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);

%get the locations of stations, here PGC-trio is enough
dtdir = ('/home/data2/chaosong/matlab/allan/2004/JUL/');
flist= ('datalist');
outf = ('staloc.txt');
stainfo = getstainfo(dtdir, flist, outf);
ind=[5,8,6];
pgc = stainfo(ind,:);
for i = 1: length(ind)
    [dxpgc(i),dypgc(i)] = absloc2relaloc(str2num(pgc(i,3)),str2num(pgc(i,2)),loc0(3),loc0(4));
end
plot(ax,dxpgc,dypgc,'^','MarkerSize',8,'MarkerEdgeColor','k','markerf','r');
text(ax,dxpgc+1.5,dypgc+1.5,pgc(:,1),'FontSize',8);

bostname = ('/BOSTOCK/total_mag_detect_0000_cull_NEW.txt');
catalog = load(strcat(workpath, bostname));
if timecol == 2
  suffix = 'zc';
elseif timecol == 3
  suffix = 'ppk';
elseif timecol == 4
  suffix = 'apk';
end
print(f.fig,'-dpdf',strcat(workpath,'/BOSTOCK/bostcattimecor_',suffix,'.pdf'));













