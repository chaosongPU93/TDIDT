% rm002246overlap.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Through 'locateMBlfe.m' by locating all MB LFEs in the newest, most reliable
% catalog 'total_mag_detect_0000_cull_NEW.txt' which has the updated moment
% estimate, we realize that 002 and 246 are very close in space. To avoid 
% double counting and include all unique LFEs, we need to lump 2 fams, if 
% we want to use the geodetic constraint for total moment, etc.
% Here we firts plot out LFEs from both fams, and see which ones are unique,
% then come up a way to combine them.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/03/22
% Last modified date:   2023/03/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

famnum = [2; 246];

format short   % Set the format to 5-digit floating point

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

%%%load MB's catalog
%load the LFE catalog of Michael Bostock, inside and outside the rectangle in 'locinterp002_4s.m'
%his time is approximately corrected to the postive waveform peak 
bosdir = ('/home/data2/chaosong/matlab/allan/BOSTOCK');
lfefnm = ('newlfeloc');
lfeloc = load(fullfile(bosdir, lfefnm));
loc0 = lfeloc(lfeloc(:,1)==2,:); %obtain the location of fam 002, lon0 and lat0

%this specifies which MB LFE catalog to use, 'AMR' is updated within 2-4 Hz
%while 'NEW' is within 0.5-1 Hz
% flag = 'NEW';
flag = 'AMR';
bostname = strcat('/BOSTOCK/update20230916/total_mag_detect_0000_cull_',flag,'.txt');
bostcat = ReformBostock(loc0(3),loc0(2),1,bostname);
% keyboard
bostcat = bostcat((bostcat(:,2)==2003 & bostcat(:,3)==3) | ...
                  (bostcat(:,2)==2004 & bostcat(:,3)==7) | ...
                  (bostcat(:,2)==2005 & bostcat(:,3)==9), :);
%format: [fam yyyy mm dd sec dx dy lon lat dep magnitude number-of-stations], 12 cols
%time should point to the peak, zero-crossing?

%choose a particular fam
fam1 = num2zeropadstr(famnum(1), 3)
bost1 = bostcat(bostcat(:,1)==famnum(1),:);
dates1 = unique(bost1(:,2)*1e4 + bost1(:,3)*1e2 + bost1(:,4));
nd1 = length(dates1);

fam2 = num2zeropadstr(famnum(2), 3)
bost2 = bostcat(bostcat(:,1)==famnum(2),:);
dates2 = unique(bost2(:,2)*1e4 + bost2(:,3)*1e2 + bost2(:,4));
nd2 = length(dates2);

dates = unique([dates1; dates2]);
nd = length(dates);
% isequaln(dates1,da)
 
%Sampling rate of the data will be used
sps=160;     % samples per second
  

%% identify and locate each MB's LFEs 
nlfe1 = 0;
nlfe2 = 0;

% bostuni1 = [];
% bostuni2 = [];
% bostdup = [];
bostuni = [];

pltflag = 0;

maxdtime = 1;
maxdmag = 0.1;

%loop over days
for id = 1: nd
  fprintf('%d / %d \n', id, nd);
  date = dates(id);
  yr = floor(date/1e4);
  mo = floor((date-yr*1e4)/1e2);
  dy = date-yr*1e4-mo*1e2;

  jday = dat2jul(mo,dy,yr);
  JDAY = num2zeropadstr(jday,3);
  YEAR=int2str(yr);
  MO=day2month(jday,yr);     % EXTERNAL function, day2month, get the month of one particular date
  
  %Bostock's LFE catalog on the same date
  bosti1 = bost1(bost1(:,2)==yr & bost1(:,3)==mo & bost1(:,4)==dy,:);
  bosti1 = sortrows(bosti1, 5);
  nbosti1 = size(bosti1,1);  
  nlfe1 = nlfe1+nbosti1;
  tbosti1 = bosti1(:,5);
  mbosti1 = bosti1(:,end-1);
  
  bosti2 = bost2(bost2(:,2)==yr & bost2(:,3)==mo & bost2(:,4)==dy,:);
  bosti2 = sortrows(bosti2, 5);
  nbosti2 = size(bosti2,1);  
  nlfe2 = nlfe2+nbosti2;
  tbosti2 = bosti2(:,5);
  mbosti2 = bosti2(:,end-1);
  
  maxamp(id,1) = max([mbosti1; mbosti2]);
  
  idup1 = [];
  idup2 = [];
  iuni1 = [];
  iuni2 = [];
  for j = 1:nbosti1
%     ind = find(abs(tbosti2-tbosti1(j))<=maxdtime & abs(mbosti2-mbosti1(j))<=maxdmag);
    ind = find(abs(tbosti2-tbosti1(j))<=maxdtime);
    if length(ind)>1  %usually just 1 duplicate
      fprintf('Multiple duplicates found, check %dth LFE of fam 002',j);
    end
    if ~isempty(ind)
      idup1 = [idup1; j];
      idup2 = [idup2; ind];
    else
      iuni1 = [iuni1; j];
    end
    
  end
  iuni2 = setdiff(1:nbosti2,idup2);
  bostuni1 = bosti1(iuni1,:); %unique sources from 002
  bostuni2 = bosti2(iuni2,:); %unique sources from 246
  bostdup = bosti1(idup1,:);  %in this case, always save the duplicate from 002
  tmp = [bostuni1; bostuni2; bostdup];
  tmp = sortrows(tmp,5);
  bostuni = [bostuni; tmp];
  
  if pltflag
    f = initfig(15,9,6,1); %initialize fig
    xran = [0.05 0.95]; yran = [0.08 0.92];
    xsep = 0.02; ysep = 0.04;
    optaxpos(f,6,1,xran,yran,xsep,ysep);
    for i = 1: 6
      ax=f.ax(i); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
      if ~isempty(bosti1)
        p1=scatter(ax,tbosti1,mbosti1,15,'b','filled');  %fam 002
      end
      if ~isempty(bosti2)
        p2=scatter(ax,tbosti2,mbosti2,15,'r','filled');  %fam 246  
      end
      xlim(ax,[(i-1)*86400/6 i*86400/6]);
      ylim(ax,[0.5 2.5]);
    end
    xlabel(ax,sprintf('Time (s) on %d/%d/%d',dy,mo,yr));
    ylabel(ax,'Magnitude');
    text(f.ax(1),0.9,0.9,sprintf('%.2f',maxamp(id)),'Units','normalized',...
      'HorizontalAlignment','right');
    legend(ax,[p1,p2],fam1,fam2);
  end

end
  
%% save lumped catalog
fid = fopen(strcat(bosdir,'/002-246_lumped.2003-2005_cull_',flag,'_chao'),'w');
%format: [fam yyyy mm dd sec dx dy lon lat dep magnitude number-of-stations]
fprintf(fid,'%d %d %d %d %e %f %f %f %f %f %f %d \n',bostuni');
fclose(fid);









  