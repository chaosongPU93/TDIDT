%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This is the script to plot the location of the stations in your data and
% the LFE families with their names
% 
% Features:
%   
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/04/24
% Last modified date:   2019/04/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

dtdir = ('/home/data2/chaosong/matlab/allan/2004/JUL/');
flist= ('datalist');
outf = ('staloc.txt');

stainfo = getstainfo(dtdir, flist, outf);
ind=[3,11,4];
lzbtrio = stainfo(ind,:);

ind=[5,8,6];
pgctrio = stainfo(ind,:);

bosdir = ('/home/data2/chaosong/matlab/allan/BOSTOCK');
lfefnm = ('newlfeloc');
lfeloc = load(fullfile(bosdir, lfefnm));    
nlfe = size(lfeloc,1);

f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 10;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
hold on
%%% lfeloc is loctions of all fams, but not necessarily each fam has detections
% scatter(lfeloc(:,3), lfeloc(:,2), 12, 'r', 'filled', 'MarkerEdgeColor', 'k'); hold on
for jj=1:3
    % p1 is lzb trio stations
    p1 = plot(str2num(lzbtrio(jj,3)),str2num(lzbtrio(jj,2)),'^','MarkerSize',8,'MarkerEdgeColor',...
        'k','markerf','k');
    text(str2num(lzbtrio(jj,3))+0.01,str2num(lzbtrio(jj,2))+0.01,lzbtrio(jj,1));
    % p2 is pgc trio stations
    p2 = plot(str2num(pgctrio(jj,3)),str2num(pgctrio(jj,2)),'^','MarkerSize',8,'MarkerEdgeColor',...
        'k','markerf','r');
    text(str2num(pgctrio(jj,3))+0.01,str2num(pgctrio(jj,2))+0.01,lzbtrio(jj,1));    
end
box on

%%% find id of families that Yajun used
famused = [2,13,25,28,56,84,99,115,125,147,149];
%%% see if all families's locations can be found 
[~,~,ind] = intersect(famused, lfeloc(:, 1));
if size(ind,1) == length(famused)
    disp('Locations of all fams used by Yajun can be found');
else
    disp('WRONG! Locations of some fams used by Yajun cannot be found'); 
end
% p3 are the fams used by Yajun
p3 = scatter(lfeloc(ind,3), lfeloc(ind,2), 40, ' b', 'filled', 'MarkerEdgeColor', 'k');
text(lfeloc(ind,3),  lfeloc(ind,2), num2str(lfeloc(ind,1)), 'fontsize',12,'color','k');


%%% set path and name
catapath= ('/home/data2/chaosong/matlab/allan/BOSTOCK/');
cataname = ('total_mag_detect_0000_cull_NEW.txt');
catalog = load(strcat(catapath, cataname));
%%% find the families that have LFE detections
famall = unique(catalog(:, 1));

%%% find the intersection between Yajun's used LFE fam and fam have
%%% detections, this the evidence that most fam are changed, so that old and
%%% new fam pool are different
[~,~,ind] = intersect(famused,famall);
famava =  famall(ind);
if length(famava) == length(famused)
    disp('All fams used by Yajun have detections');
else
    disp('WRONG! Some fams used by Yajun do not have detections');
end

%%% find the intersection between fam have detections and fam have locations
[~,~,ind] = intersect(famall, lfeloc(:, 1));
% p4, fams that have locations and detections
p4 = scatter(lfeloc(ind,3), lfeloc(ind,2), 12, 'r', 'filled', 'MarkerEdgeColor', 'k'); hold on
text(lfeloc(ind,3),  lfeloc(ind,2), num2str(lfeloc(ind,1)), 'fontsize',12,'color','g');

%%% plot the location of new fam pool
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
[~,~,ind] = intersect(str2num(nfampool(:,:)), lfeloc(:, 1),'stable');
% p5, fams used by Chao from the compromise to both close locations and # of detections
p5 = scatter(lfeloc(ind,3), lfeloc(ind,2), 60, 'y', 'filled', 'MarkerEdgeColor', 'k'); hold on
text(lfeloc(ind,3),  lfeloc(ind,2), num2str(lfeloc(ind,1)), 'fontsize',12,'color','c');

famloc = lfeloc(ind,:);


%%% according to the hit density plot from fam 002 with PGC trio, there is a spatial difference
%%% close to (123.4,48.525); (123.53,48.517), but close to the detection edge of 002, so try to
%%% use another fam to resolve this area
region = [-123.55 -123.35 -123.35 -123.55 -123.55;
          48.515 48.515 48.53 48.53 48.515];
plot(region(1,:), region(2,:), 'k','linew',1);

legend([p1,p2,p3,p4,p5],{'LZB trio','PGC trio','Fams by Yajun','Fams have detections',...
       'Fams by Chao'},'location','southwest');


workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/lzbrst');

print('-dpdf',strcat(rstpath,'/entire_sta_fam_distri.pdf'));












