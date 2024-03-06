% tremorbursts002_4s.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The next script to run after 'locinterp002_4s.m', is to read in the tremor
% detections inside a cutoff region output by 'locinterp002_4s.m' that has
% the highest density. Then to this catalog, obtain the inter-event time.
% By setting a max separation time and min # of detections inside, we extract
% the time of tremor bursts that are clustered in time. The time windows are
% the targets for the future deconvolution.
%
%
% --I shall reiterate the order of the entire flow of detection here: 
%   'detection002_4s.m' --> 'mergeday002_4s.m' --> 'locinterp002_4s.m'
%   --> 'tremorbursts002_4s.m'
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/03/15
% Last modified date:   2022/03/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio

format short e   % Set the format to 5-digit floating point
clear
clc
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

scrsz=get(0,'ScreenSize');

%WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band 
%(4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
%%% choose the data path here 
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/PGCtrio');

freqflag='hf';  % flag to indicate whether to do hf or lf;

FLAG = 'PGC'; % detector
  
fam = '002';   % family number

sps = 40;

iup = 4;  % upsample 4 times

cutout = 'ellipse';

%load detections inside the cutout boundary of interest, an output of 'locinterp002_4s.m'
if isequal(fam,'002')
  cyclskip = 0;
  mshift=26+cyclskip; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
  loopoffmax=2.1; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
  xcmaxAVEnmin=0.44; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?
end
hi=6.5;    % frequency band
lo=1.25;
npo=2;     % poles, passes of filters
npa=2;
winlensec=4;
winoffsec=1;        % window offset in sec, which is the step of a moving window
winlen=winlensec*sps;      % length in smaples
PREFIX = strcat(fam,'.up.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),'.',...
                        int2str(npo),int2str(npa),'.ms', int2str(mshift));
hfbnd = load(strcat(rstpath, '/MAPS/eloc.pgc.tdec_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',...
            num2str(winlen/sps),'s',num2str(sps),'sps4add_',cutout(1:4)));
daycol = 14;
seccol = 16;  % 2 choices here, one is the center of detecting window, other is the start of strongest .5 s 
hfbnd = sortrows(hfbnd, [daycol, seccol]);

%load all detections, an output of 'locinterp002_4s.m'
hfall = load(strcat(rstpath, '/MAPS/eloc.pgc.tdec_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',...
            num2str(winlen/sps),'s',num2str(sps),'sps4add_'));
hfall = sortrows(hfall, [daycol, seccol]);
%get ones outside your boundary
hfout = setdiff(hfall,hfbnd,'rows');           
hfout = sortrows(hfout, [daycol, seccol]);

%%
hfanalyse = hfbnd;
% hfanalyse = hfout;

%% compute relative time in days
%sort according to day, sec of the strongest arrival of the window
daycol = 14;
seccol = 16;
hfanalyse = sortrows(hfanalyse, [daycol, seccol]);
ncol = size(hfanalyse,2);

%obtain the relative time in days
i03 = find(hfanalyse(:,daycol) < 2004*1000);
hfanalyse(i03,ncol+1) = (hfanalyse(i03,daycol)-2003060)+hfanalyse(i03,seccol)./(3600.*24);
i04 = find(hfanalyse(:,daycol) < 2005*1000 & hfanalyse(:,daycol) > 2004*1000);
hfanalyse(i04,ncol+1) = (hfanalyse(i04,daycol)-2004194)+hfanalyse(i04,seccol)./(3600.*24);
i05 = find(hfanalyse(:,daycol) > 2005*1000);
hfanalyse(i05,ncol+1) = (hfanalyse(i05,daycol)-2005254)+hfanalyse(i05,seccol)./(3600.*24);
          
%% separation in time between itself and its preceding detection
%for the detections, obtain the separation in time between itself and its preceding
%detection for the catalog of each ETS.
%1st col: occurence time
%2nd col: inter-detection time to its preceding detection
hfinter03 = interevt_time(hfanalyse(i03,end));
hfinter04 = interevt_time(hfanalyse(i04,end));
hfinter05 = interevt_time(hfanalyse(i05,end));
hfinter = [hfinter03; hfinter04; hfinter05];
          
%% plot the inter-detection time for all 3 ETS, in days
%%% HF catalog
%set a threshold of inter-detection time
ttol1 = 1e-3*ones(3,1); %in days
% ttol2 = 5e-4*ones(3,1);
ttol2 = [median(hfinter03(:,2)); median(hfinter04(:,2)); median(hfinter05(:,2))];
ymax = 1.2*max([ttol2; ttol1]);
[f] = plt_interevt_time_diffets(hfinter03,hfinter04,hfinter05,ymax,ttol1,ttol2);
% text(f.ax(1),0.04,0.9,'HF','FontSize',11,'unit','normalized','horizontalalignment','left',...
%     'EdgeColor','k','Margin',2);
% print(f.fig,'-dpdf',strcat(figpath,'/inter-detection_time_HF',FLAG,bndflag,'.pdf')); 
          
%% group the tremor bursts 
% ttol = 1e-3*ones(3,1); 
% ntol = 10;
ttol = 35/86400*ones(3,1);  % this is about 4*median
ntol = 3;
% ttol = 35/86400*ones(3,1);  % this is about 4*median
% ntol = 3;
[bursthf03, n03] = group_tremor_burst_indep(hfinter03,ttol(1),ntol);
[bursthf04, n04] = group_tremor_burst_indep(hfinter04,ttol(2),ntol);
[bursthf05, n05] = group_tremor_burst_indep(hfinter05,ttol(3),ntol);

% recover the index of the detections to occurence times with the same format as trange
if ~isempty(bursthf03)
  [thf03,tranhf03,ntot03] = burst_range(bursthf03,hfinter03,2003060);
  perchf03 = ntot03/length(i03)*100;
else
  thf03 = []; tranhf03 = []; ntot03 = 0; perchf03 = 0;
end
if ~isempty(bursthf04)
  [thf04,tranhf04,ntot04] = burst_range(bursthf04,hfinter04,2004194);
  perchf04 = ntot04/length(i04)*100;
else
  thf04 = []; tranhf04 = []; ntot04 = 0; perchf04 = 0; 
end
if ~isempty(bursthf05)
  [thf05,tranhf05,ntot05] = burst_range(bursthf05,hfinter05,2005254);
  perchf05 = ntot05/length(i05)*100;
else
  thf05 = []; tranhf05 = []; ntot05 = 0; perchf05 = 0;
end

%% obtain the separation time between each burst and the length in time of each burst
tsep03 = zeros(size(thf03,1),1);   %what is the separation in time between each window?
tlen03 = (thf03(:,2)-thf03(:,1))*86400;  %what is the length of each window?
for i = 2: size(thf03,1)
  tsep03(i) = (thf03(i,1)-thf03(i-1,2))*86400;
end
tsep04 = zeros(size(thf04,1),1);
tlen04 = (thf04(:,2)-thf04(:,1))*86400;
for i = 2: size(thf04,1)
  tsep04(i) = (thf04(i,1)-thf04(i-1,2))*86400;
end
tsep05 = zeros(size(thf05,1),1);
tlen05 = (thf05(:,2)-thf05(:,1))*86400;
for i = 2: size(thf05,1)
  tsep05(i) = (thf05(i,1)-thf05(i-1,2))*86400;
end
tsep = [tsep03; tsep04; tsep05];
tlen = [tlen03; tlen04; tlen05];

figure
subplot(221)
scatter(tlen03,log10(tsep03),15,'r','filled'); hold on
xlabel('Length in time (s) of each burst ');
ylabel('log_{10}(Separation in time (s) between bursts)');
scatter(tlen04,log10(tsep04),15,'b','filled');
scatter(tlen05,log10(tsep05),15,'k','filled');
ax=gca; plot(ax.XLim,[log10(ttol(1)*86400) log10(ttol(1)*86400)], 'k--');
box on; grid on;
legend('2003','2004','2005','T_{max}','Location','northeast');

subplot(222)
scatter(tlen,log10(tsep),15,[.4 .4 .4],'filled'); hold on
xlabel('Length in time (s) of each burst ');
ylabel('log_{10}(Separation in time (s) between bursts)');
ax=gca; plot(ax.XLim,[log10(ttol(1)*86400) log10(ttol(1)*86400)], 'k--');
box on; grid on;
legend('All','T_{max}','Location','northeast');

subplot(223)
scatter(tlen03,tsep03,15,'r','filled'); hold on
xlabel('Length in time (s) of each burst ');
ylabel('Separation in time (s) between bursts');
scatter(tlen04,tsep04,15,'b','filled');
scatter(tlen05,tsep05,15,'k','filled');
axis([-10 400 0 100]);
ax=gca; plot(ax.XLim,[ttol(1)*86400 ttol(1)*86400], 'k--');
box on; grid on;
legend('2003','2004','2005','T_{max}','Location','northeast');

subplot(224)
scatter(tlen,tsep,15,[.4 .4 .4],'filled'); hold on
xlabel('Length in time (s) of each burst ');
ylabel('Separation in time (s) between bursts');
axis([-10 400 0 100]);
ax=gca; plot(ax.XLim,[ttol(1)*86400 ttol(1)*86400], 'k--');
box on; grid on;
legend('All','T_{max}','Location','northeast');
text(0.5,0.4,sprintf('<=40 s: %d/%d',sum(tsep<=40)-3,size(tsep,1)),'Units','normalized');
text(0.5,0.5,sprintf('<=50 s: %d/%d',sum(tsep<=50)-3,size(tsep,1)),'Units','normalized');
text(0.5,0.6,sprintf('<=60 s: %d/%d',sum(tsep<=60)-3,size(tsep,1)),'Units','normalized');
text(0.5,0.7,sprintf('<=80 s: %d/%d',sum(tsep<=80)-3,size(tsep,1)),'Units','normalized');


%% plot the grouped bursts for all 3 ETS, in days
ymax = 1e3;
[f] = plt_burst_diffets(ymax,thf03,thf04,thf05,hfinter03,hfinter04,hfinter05);
text(f.ax(1),0.06,0.5,strcat(num2str(round(perchf03)),'%'),'FontSize',12,'unit','normalized',...
    'horizontalalignment','left');
text(f.ax(2),0.06,0.5,strcat(num2str(round(perchf04)),'%'),'FontSize',12,'unit','normalized',...
    'horizontalalignment','left');
text(f.ax(3),0.06,0.5,strcat(num2str(round(perchf05)),'%'),'FontSize',12,'unit','normalized',...
    'horizontalalignment','left');

%% summarize the inter-event time and grouped bursts
ymax = 1e-3;
[f] = plt_intertimeandbursts(hfinter03,hfinter04,hfinter05,thf03,thf04,thf05,ymax,ttol);

keyboard

%% plot the inter-detection time for each ETS in more detail, in sec
ttol1 = 1e-3*ones(3,1)*86400; 
ttol2 = ttol*86400; 
ymax = 1.2*max([ttol2; ttol1]);

[f] = plt_interevt_time_burst(thf03,hfinter03,2003060,ymax,ttol1,ttol2);

[f] = plt_interevt_time_burst(thf04,hfinter04,2004194,ymax,ttol1,ttol2);

[f] = plt_interevt_time_burst(thf05,hfinter05,2005254,ymax,ttol1,ttol2);

%% some simple statistics of the bursts
%times of all tremor bursts
trange = [tranhf03;tranhf04;tranhf05];

%what is the typical length of each window? Useful for plotting all of them later
tlen = trange(:,3)-trange(:,2);

figure
histogram(tlen,'BinWidth',10); hold on;
text(0.5,0.9,sprintf('Med: %.1f s',median(tlen)),'Units','normalized');
xlabel('Length in time (s) of each burst ');
ylabel('Counts');

% get the new catalog composed only by the bursts
hfnew = refine_catalog_by_trange(hfanalyse,daycol,seccol,trange);
hfnew = hfnew(:,1:end-1); % ignore the last column of relative time in days

perchfnew = size(hfnew,1)/size(hfanalyse,1)

% figure
% histogram(log10(tsep),'BinWidth',1); hold on;
% text(0.1,0.9,sprintf('Med: %d s',median(log10(tsep))),'Units','normalized');
% xlabel('log_{10}(Separation in time (s) between bursts)');
% ylabel('Counts');

keyboard

%% output results to file
%format as follows:
%%% 8+34+4*nstanew cols, if 4 new stas, then will be 58 cols
%%% UPDATED at 2021/06/23
%%%   dx,dy,lon,lat,dep,tori,off12,off13 (integer samples at upsampled sps)
%%%   1:off12(n) 2:off13(n) 3:off12sec 4:off13sec 5:fam 6:date 7:timswin(n)
%%%   8:timswin(n)-winlensec/2+idiff/sps 9:xcmaxAVEnbang(nin) 10:loopoff(n)
%%%   11:cumsumtrdiff  12:cumsumtrdiff/cumsumtr(winlen)
%%%   for the rest, n=n+4;
%%%   9:Ampsq(nin) 10:Prior(nin) 11:Post(nin) 12:Prior4(nin) 13:Post4(nin) 14:Prior8(nin)
%%%   15:Post8(nin) 16:Prior16(nin) 17:Post16(nin) 18:cc(nin) 19:ccprior(nin)
%%%   20:ccpost(nin) 21:ccprior4(nin) 22:ccpost4(nin) 23:STAamp(nin,2) 24:STAamp(nin,3)
%%%   25:sumsSTA12n(n,iSTA12bang) 26:sumsSTA13n(n,iSTA13bang)
%%%   27:sumsSTA32n(n,iSTA32bang) 28:sigma(nin) 29:sigma12(nin) 30:sigma13(nin)
%%%   in(1:nstanew) loff(1:nstanew) ioff(1:nstanew) ccmaxave(1:nstanew)

% fid = fopen(strcat(rstpath, '/MAPS/eloc.pgc.tdec.bst',num2str(round(ttol(1)*86400)),'s_',PREFIX,'_',...
%   num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps4add_',cutout(1:4)),'w+');
% fprintf(fid,'%.4f %.4f %.4f %.4f %.4f %.4f %d %d %6.2f %6.2f %.4f %.4f %d %d %9.1f %10.3f %8.3f %7.2f %10.3e %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f \n',...
%         hfnew');
  
if isequaln(hfanalyse(:,1:end-1),hfbnd)
  fid = fopen(strcat(rstpath, '/MAPS/tdec.bstran',num2str(round(ttol(1)*86400)),'s.',...
    num2str(ntol),'.pgc002.',cutout(1:4)),'w+');
  % fid = fopen(strcat(rstpath, '/MAPS/tdec.bstran',num2str(round(ttol(1)*86400)),'s.pgc002.',cutout(1:4)),'w+');

elseif isequaln(hfanalyse(:,1:end-1),hfout)
  fid = fopen(strcat(rstpath, '/MAPS/tdec.bstran',num2str(round(ttol(1)*86400)),'s.',...
    num2str(ntol),'.pgcout002.',cutout(1:4)),'w+');

end

fprintf(fid,'%d %9.2f %9.2f \n',trange');
fclose(fid);


