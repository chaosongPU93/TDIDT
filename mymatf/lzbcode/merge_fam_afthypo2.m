%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the script to merge the detections from all families IF you want 
% to set a distance cutoff first, then merge fams to remove duplicates base
% on time and CC coef, so run this AFTER hypoinverse.
% The input is already been removed off the double counting. Since for each
% fam, it wiped put the same data, so the detections should still have a lot of
% duplicates, preserve the ones within the one duration that have the highest 
% CC value. In the mean time, detections in each family have a different
% reference (0,0) point, so the merging can also convert all to the same
% frame (similar to changing centroid).
%
% V2:
%   Try to use a new way to remove the double counted detections, ABONDONED!
%
%
% NOTE :
%   1*. This should be ran AFTER hypoinverse, before do migration analyzing!!!
%       
%        
%   2. run with fcflag=0 (DEFAULT) would save the no-duplicate results for 
%       plotaddcheck_LZB_allfam and figures; run with fcflag=1 would correct 
%       the filtering effect for no-duplicate results to generate and save new
%        figures, but not save files
%   3. run with convflag=0 (DEFAULT) would keep the results in its own fam frame;
%       run with convflag=1 would convert results of all fams to the same frame
%       of the reference fam, which is now 043, can be changed easily.
%   4. IN DEFAULT, in order to write rst files for next-step plotaddcheck scripts, set
%       fcflag=0 and convflag=0; other options ONLY save plots to eyeball analysis
%   
%        
% Conversion algorithm:
%       (off12')_fam - (off12')_ref = (off12)_fam - (off12)_ref
%   where ' means the detection in that family, no ' means the original
%   reference point in each family, i.e. the location of that family 
%   (considered as 0,0 in each fam). Since we already have the travel time
%   difference between sta 12 & 13 in each family (4th col in rots para), 
%   what we need to do here is to convert them to a uniform frame.
%   Thus, relative to a selected fam, all detections in other fam can be
%   written in the new frame as:
%       (off12')_ref = (off12')_fam - [(off12)_fam - (off12)_ref]
%
% Features:
%   1. read detection offsets and trace files
%   2. option flag to correct filtering effect (1/0)
%   3. option flag to convert all fams to the same frame (1/0)
%   4. remove double counting, duplicates, NECESSARY
%   5. re-write those files (only when fcflag=0 and convflag=0)
%   6. re-plot the offsets variation w/ time after removing duplicates (always save plots)
%   
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/11/25
% Last modified date:   2019/11/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
clc
% close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/lzbrst');

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
        
winlenhf = 4;

winlenlf = 16;
lofflf = 4;
ccminlf = 0.35;


%% load detections with or without a distance cutoff
SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));

dcutflag = 0;

if dcutflag 
    hffname = strcat(rstpath, '/evtloc.allfam.dcut.',SUFFIXhf);
    lffname = strcat(rstpath, '/evtloc.allfam.dcut.',SUFFIXlf);
else
    hffname = strcat(rstpath, '/evtloc.allfam.nodcut.',SUFFIXhf);
    lffname = strcat(rstpath, '/evtloc.allfam.nodcut.',SUFFIXlf);
end

hfold = load(hffname);
lfold = load(lffname);
% 23 cols, format is:
%   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date 
%   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin famnum
%



%% sort according to time, remove double counting
disp('Start removing double counts ...');

%%% for hf
%%% sort according to date and main arrival time, could add offset
hfsort = sortrows(hfold, [13,14]);
dateall = unique(hfsort(:,13));
hfnew = [];
hfdup = [];
hfdupsel = [];
hfdist = [];
for id = 1: length(dateall)
    hfday = hfsort(hfsort(:,13) == dateall(id), :);
    dtminhf = 0.5;      % min time during which only 1 detection is retained
    colnum = [14 16];
    indexthf = length(nfampool);    % the extreme case is that every fam contains a duplicate
    [hfdaynew,indsave,inddis,dupday,dupselday,distday] = ...
            RemoveDoubleCounting2(hfday,dtminhf,indexthf,colnum);
    hfnew = [hfnew; hfdaynew];    % hfdaynew is the file of all fam for one day  
    hfdup = [hfdup; dupday];
    hfdupsel = [hfdupsel; dupselday];
    hfdist = [hfdist; distday];
end

%%% for lf
%%% sort according to date and main arrival time, could add offset
lfsort = sortrows(lfold, [13,14]);
dateall = unique(lfsort(:,13));
lfnew = [];
lfdup = [];
lfdupsel = [];
lfdist = [];
for id = 1: length(dateall)
    lfday = lfsort(lfsort(:,13) == dateall(id), :);
    dtminlf = 1.1;      % min time during which only 1 detection is retained
    colnum = [14 16];
    indextlf = length(nfampool);    % the extreme case is that every fam contains a duplicate
    [lfdaynew,indsave,inddis,dupday,dupselday,distday] = ...
            RemoveDoubleCounting2(lfday,dtminlf,indextlf,colnum);
    lfnew = [lfnew; lfdaynew];    % hfdaynew is the file of all fam for one day  
    lfdup = [lfdup; dupday];
    lfdupsel = [lfdupsel; dupselday];
    lfdist = [lfdist; distday];
end

%%
f.fig=figure;
f.fig.Renderer='Painters';
widin = 4;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
f.ax = gca;
hold(f.ax, 'on');
f.ax.Box = 'on';
grid(f.ax, 'on');
f.ax.GridLineStyle = '--';
% patarea = [0 0;
%            2.5 0;
%            2.5 1;
%            0 1;
%            0 0];
% patch(f5.ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.1,'edgecolor','none');
h = histogram(f.ax,hfdist,'binwidth',0.25,'normalization','probability',...
              'facecolor','k','facea',0.5);
perc = sum(abs(hfdist)<=2)/size(hfdist,1);
perc90 = prctile(hfdist,90);
% for i = 1: length(h.Values)
%     text(f2.ax,h.BinEdges(i)+0.1,h.Values(i),sprintf('%.2f',h.Values(i)));
% end
text(f.ax,0.5,0.95,num2str(size(hfdist,1)),'fontsize',13,'unit','normalized');
text(f.ax,0.8,0.95,'HF','fontsize',13,'unit','normalized');
plot([2, 2], [0 1],'b--','linew',1);
text(f.ax,2,0.15,sprintf('%.2f',perc),'fontsize',13);
plot([perc90, perc90], [0 1],'r--','linew',1);
text(f.ax,perc90,0.15,'90 percentile','fontsize',13);
xlim(f.ax,[0 12]);
ylim(f.ax,[0 0.2]);
xlabel(f.ax,'Distance (km)');
ylabel(f.ax,'Fraction');
hold(f.ax, 'off');

print(f.fig,'-dpdf',strcat(rstpath,'/hf_duplicate_distance.pdf'));



f2.fig=figure;
f2.fig.Renderer='Painters';
widin = 4;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
f2.ax = gca;
hold(f2.ax, 'on');
f2.ax.Box = 'on';
grid(f2.ax, 'on');
f2.ax.GridLineStyle = '--';
% patarea = [0 0;
%            2.5 0;
%            2.5 1;
%            0 1;
%            0 0];
% patch(f5.ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.1,'edgecolor','none');
h = histogram(f2.ax,lfdist,'binwidth',0.5,'normalization','probability',...
              'facecolor','k','facea',0.5);
perc = sum(abs(lfdist)<=4)/size(lfdist,1);
perc90 = prctile(lfdist,90);
% for i = 1: length(h.Values)
%     text(f2.ax,h.BinEdges(i)+0.1,h.Values(i),sprintf('%.2f',h.Values(i)));
% end
text(f2.ax,0.5,0.95,num2str(size(lfdist,1)),'fontsize',13,'unit','normalized');
text(f2.ax,0.8,0.95,'LF','fontsize',13,'unit','normalized');
plot([2, 2], [0 1],'b--','linew',1);
text(f2.ax,4,0.15,sprintf('%.2f',perc),'fontsize',13);
plot([perc90, perc90], [0 1],'r--','linew',1);
text(f2.ax,perc90,0.15,'90 percentile','fontsize',13);
xlim(f2.ax,[0 12]);
ylim(f2.ax,[0 0.2]);
xlabel(f2.ax,'Distance (km)');
ylabel(f2.ax,'Fraction');
hold(f2.ax, 'off');

print(f2.fig,'-dpdf',strcat(rstpath,'/lf_duplicate_distance.pdf'));







%%
%%% save the results here
% 23 cols, format is:
%   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date 
%   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin famnum
%
if dcutflag 
    fidhf = fopen(strcat(rstpath, '/evtloc.allfam.dcutnodou.',SUFFIXhf),'w+');
    fidlf = fopen(strcat(rstpath, '/evtloc.allfam.dcutnodou.',SUFFIXlf),'w+');
else
    fidhf = fopen(strcat(rstpath, '/evtloc.allfam.nodcutnodou.',SUFFIXhf),'w+');
    fidlf = fopen(strcat(rstpath, '/evtloc.allfam.nodcutnodou.',SUFFIXlf),'w+');
end

% if dcutflag 
%     fidhf = fopen(strcat(rstpath, '/aaaa.dcutnodou.',SUFFIXhf),'w+');
%     fidlf = fopen(strcat(rstpath, '/bbbb.dcutnodou.',SUFFIXlf),'w+');
% else
%     fidhf = fopen(strcat(rstpath, '/aaaa.nodcutnodou.',SUFFIXhf),'w+');
%     fidlf = fopen(strcat(rstpath, '/bbbb.nodcutnodou.',SUFFIXlf),'w+');
% end

fprintf(fidhf,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %d %.4f %.1f %.4f %.4f %.4f %.4f %.4f %.4e %.4e %d \n',...
        hfnew');
fclose(fidhf);

fprintf(fidlf,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %d %.4f %.1f %.4f %.4f %.4f %.4f %.4f %.4e %.4e %d \n',...
        lfnew');
fclose(fidlf);











