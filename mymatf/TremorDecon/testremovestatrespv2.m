%testremovestatrespv2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a follow-up script to 'testremovestatresp.m', in order to compare 
% the processed data at the same station that used different instrument 
% response files instrument response removal (IRR).
% 
% --In 2019, I did IRR using 2019 SAC_PZs, which has been confirmed to be 
% almost identical to IRR using 2019 RESPs in 'testremovestatresp.m'. In case
% these RESP files are not correct, I moved the original IRR data to folders
% like '/home/data2/chaosong/matlab/allan/data-no-resp/arch2003/'.
%
% --In 2022/06/10, CNDC has confirmed several nominal RESP files for the CN
% (permanent) and PO (polaris) stations, according to the table file that 
% records the time of instrument change. These RESP files are stored in
% '/home/data2/chaosong/matlab/allan/data-no-resp/stdresp/CN/ and PO/'.
% Therefore, I have done the IRR once
% again using these standard RESP files carefully and store the processed
% data in folders like '/home/data2/chaosong/matlab/allan/data-no-resp/2003/'.
%
% Check '/home/data2/chaosong/Notes/20220608.md' for
% how the RESPs are obtained.
% 
% --After test, in 2004/2005, the raw IRR data for sta PGC/LZB/SSIB/SILB/KLNB/...
%   using the 2019 PZ files or 2022 RESP files are almost the SAME; 
%   In 2003, the raw IRR19 for PGC is almost SAME as IRR22; for SSIB/SILB/...
%   IRR19 *5 is almost the same as IRR22; for KLNB,  IRR19 *4 is almost the 
%   same as IRR22.
% --In 2003, for LZB, IRR19 using pz is dramatically different from
%   IRR22 using RESP. see details in 'testremovestatrespv3.m'. It looks like 
%   the RESP files in 2003 is not reliable. Instead, we could redo LZB using
%   pz2022, or just copy from IRR using pz2019, two of which are very close.
% 
%
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/06/16
% Last modified date:   2022/06/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format short e   % Set the format to 5-digit floating point
clear
clc
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
temppath = strcat(workpath, '/templates/');
%%% choose the data path here
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/PGCtrio');

FLAG = 'PGC';
% FLAG = 'LZB';

fam = '002'

freqflag='hf';  % flag to indicate whether to do hf or lf;

% get days, 3-sta trio name, bostock's template name, zero crossing index
% of templates
timoffrot= [
            2003 061;
%             2003 062;
%             2003 063;
            2004 196;
%             2004 197;
%             2004 198;
%             2004 199;
            2005 254;
%             2005 255;
%             2005 256;
            ];
% [timoffrot,~] = GetDays4Stack(fam);
nday = size(timoffrot, 1);

%Use new LFE catalog
CATA = 'new';

% get permanent and polaris station rotation parameters, based on 40-sps data
sft2=0;     % centroid shift of station 2
sft3=0;     % centroid shift of station 3
[PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,sft2,sft3);
if ~isequal(CATA, 'fixed')
  reftime = PERMROTS(1,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
  PERMROTS(:,4) = PERMROTS(:,4)-reftime;    % to make sure that 1st station is 0
  POLROTS(:,4) = POLROTS(:,4)-reftime;
end
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;     % convert 2-3 columns to rad
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
  'LZB'];
POLSTA=['SSIB '           % polaris station names
  'SILB '
  'KLNB '
  'MGCB '
  'TWKB '];

stas=['PGC  '
      'SSIB '
      'SILB '
      'KLNB '
      'LZB  '
      'TWKB '
      'MGCB '
      ];
nsta=size(stas,1);         %  number of stations

%Sampling rate of the data will be used
sps=40;     % samples per second

hi=6.5;    % frequency band
lo=1.25;

npo=2;     % poles, passes of filters
npa=2;

%%% 'scalefact' scales templates; 'scaleseisms' scales seisms.  Strategy changes with family.
if isequal(fam,'002')
  whichtoplot=2;          % flag to which to plot
  scaleseisms=[1.0 0.76 0.95];       % scaleseisms scales seismograms
elseif isequal(fam,'068')
  whichtoplot=1;
  scaleseisms=[1.0 0.6 0.6];
else
  scaleseisms=[1.0 0.76 0.95];
end

%% START TO LOOP FOR every day
%cycle over each day:
for nd=1: nday      % num of rows, also num of days
  %     nd=1;
  
  year=timoffrot(nd,1);
  jday=timoffrot(nd,2);
  
  YEAR=int2str(year);
  JDAY = num2zeropadstr(jday,3);
  MO=day2month(jday,year);     % EXTERNAL function, day2month, get the month of one particular date

  %%%% data processed using the SAC_PZs in 2019
  direc19=[datapath, '/arch', YEAR,'/',MO,'/'];     % directory name
  prename19=[direc19,YEAR,'.',JDAY,'.00.00.00.0000.CN'];    %  path plus prefix of data file,
  disp(prename19);  
  %read horizontal optimal and orthogonal components
  [STAopt19,~,~,fileflag] = rd_daily_bpdata(year,jday,prename19,stas,PERMSTA,POLSTA,...
    PERMROTS,POLROTS,sps,lo,hi,npo,npa,[],[],[],[]);
  %rearrange the output
  timsSTA19 = STAopt19(:,1);
  STAopt19 = STAopt19(:,2:end);
  %     %scale seismogram relatively, from historical codes, maybe not ncecessary
  %     for ista = 1: nsta
  %       STAopt(:,ista+1) = STAopt(:,ista+1)/scaleseisms(ista);
  %     end
  
  %%%% data processed using the latest RESP files in 2022
  direc22=[datapath, '/', YEAR,'/',MO,'/'];     % directory name
  prename22=[direc22,YEAR,'.',JDAY,'.00.00.00.0000.CN'];    %  path plus prefix of data file,
  disp(prename22);
  %read horizontal optimal and orthogonal components
  [STAopt22,~,~,fileflag] = rd_daily_bpdata2022(year,jday,prename22,stas,PERMSTA,POLSTA,...
    PERMROTS,POLROTS,sps,lo,hi,npo,npa,[],[],[],[]);
  %rearrange the output
  timsSTA22 = STAopt22(:,1);
  STAopt22 = STAopt22(:,2:end);

  
  %%
  f.fig = figure;
  f.fig.Renderer = 'painters';
  widin = 10;  % maximum width allowed is 8.5 inches
  htin = 8;   % maximum height allowed is 11 inches
  % get the scrsz in pixels and number of pixels per inch of monitor 1
  [scrsz, res] = pixelperinch(1);
  set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
  
  nrow = nsta;
  ncol = 1;
  
  for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
  end
  
  pltxran = [0.08 0.96]; pltyran = [0.08 0.96];
  pltxsep = 0.05; pltysep = 0.05;
  axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

  for ista = 1: nsta
    ax = f.ax(ista);  hold(ax,'on');
    plot(ax,timsSTA19,STAopt19(:,ista)+1,'r');
    plot(ax,timsSTA22,STAopt22(:,ista)-1,'k'); 
    
    text(ax,0.98,0.9,stas(ista,:),'unit','normalized','HorizontalAlignment','right');
    xlabel(ax,'Time (s)');
    ylabel(ax,'Amplitude');
    ylim(ax,[-3 3]);
    if ista ==1 
      plot(ax,timsSTA22,STAopt22(:,ista)-STAopt19(:,ista),'b'); 
      text(ax,0.2,0.5,'black-red','unit','normalized','HorizontalAlignment','left','color','b');
      text(ax,0.2,0.9,'raw*1.0e-3','unit','normalized','HorizontalAlignment','left','color','r');
      text(ax,0.2,0.1,'raw*1.0e-3','unit','normalized','HorizontalAlignment','left','color','k');
      title(ax,strcat(YEAR,JDAY));
    elseif ista ==5
      plot(ax,timsSTA22,STAopt22(:,ista)*1.8-STAopt19(:,ista),'b'); 
      text(ax,0.2,0.5,'black*1.8-red','unit','normalized','HorizontalAlignment','left','color','b');
      text(ax,0.2,0.9,'raw*1.8e-3','unit','normalized','HorizontalAlignment','left','color','r');
      text(ax,0.2,0.1,'raw*1.0e-3','unit','normalized','HorizontalAlignment','left','color','k');
    elseif ista == 4
      if year==2003
        plot(ax,timsSTA22,STAopt22(:,ista)*1.9-STAopt19(:,ista),'b');   
        text(ax,0.2,0.5,'black*1.9-red','unit','normalized','HorizontalAlignment','left','color','b');
      else
        plot(ax,timsSTA22,STAopt22(:,ista)*1.5-STAopt19(:,ista),'b');  
        text(ax,0.2,0.5,'black*1.5-red','unit','normalized','HorizontalAlignment','left','color','b');
      end
      if year==2003 && jday<213 
        text(ax,0.2,0.9,'raw*7.5e-3','unit','normalized','HorizontalAlignment','left','color','r');
      else
        text(ax,0.2,0.9,'raw*1.5e-3','unit','normalized','HorizontalAlignment','left','color','r');
      end
      text(ax,0.2,0.1,'raw*1.0e-3','unit','normalized','HorizontalAlignment','left','color','k');

    else
      plot(ax,timsSTA22,STAopt22(:,ista)*1.5-STAopt19(:,ista),'b');  
      text(ax,0.2,0.5,'black*1.5-red','unit','normalized','HorizontalAlignment','left','color','b');
      if year==2003 && jday<213 
        text(ax,0.2,0.9,'raw*7.5e-3','unit','normalized','HorizontalAlignment','left','color','r');
      else
        text(ax,0.2,0.9,'raw*1.5e-3','unit','normalized','HorizontalAlignment','left','color','r');
      end
      text(ax,0.2,0.1,'raw*1.0e-3','unit','normalized','HorizontalAlignment','left','color','k');
    end
    
    hold(ax,'off');
  end
  
end
  
 


