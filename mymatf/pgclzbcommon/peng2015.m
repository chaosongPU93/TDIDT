%%% read and analyse the detections from Peng et al. 2015
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/11/23
% Last modified date:   2019/11/23


format short e   % Set the format to 5-digit floating point
clear
% close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/lzbrst');

fname = strcat(workpath,'/PAPERS/Peng_etal/Pengetal2015-sup-0002-supinfo.txt');
det2005 = load(fname);

fname = strcat(workpath,'/PAPERS/Peng_etal/Pengetal2015-sup-0003-supinfo.txt');
det2004 = load(fname);

fname = strcat(workpath,'/PAPERS/Peng_etal/Pengetal2015-sup-0004-supinfo.txt');
det2003 = load(fname);

jdmin = 254;
jdmax = 258;
refyr = 2005;
daymax = 15;
tmp = jul2dat(refyr,jdmin);   % return month, day, year
refday = tmp(2);  
det1 = det2005(det2005(:,2)>=jdmin & det2005(:,2)<=jdmax, :);
timeint = det1(:,2)-jdmin+refday;
timedec = det1(:,3)/3600/24;
time = timeint+timedec;
det1(:,15) = time;

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
           -123.838500 48.544833 35.6600];

relacont = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),loccont(2,1),loccont(2,2));
relacont(:,1) = dx;
relacont(:,2) = dy;

[dx, dy] = absloc2relaloc(det1(:,5),det1(:,4),loccont(2,1),loccont(2,2));
det1(:,6) = dx;
det1(:,7) = dy;

%%% define and position the figure frame and axes of each plot 
f.fig = figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end    

%%% reposition
for isub = 1:nrow*ncol
    if mod(isub,ncol)==0
        icol = ncol;
        irow = floor(isub/ncol);
    else
        icol = mod(isub,ncol);
        irow = floor(isub/ncol)+1;
    end
    set(f.ax(isub),'Position', [(icol-0.9)/ncol, 1-irow/nrow*0.95, 1/ncol*0.8 1/nrow*0.8]);
end


%%% deal with each subplot
cran = [refday daymax];

% marker size
msizehf = 4;
msizelf = 8;

% subplot 1 of figure 1
hold(f.ax(1),'on');
dumdet = sortrows(det1,-15);
scatter(f.ax(1),dumdet(:,6),dumdet(:,7), msizehf, dumdet(:,end), 'filled','o');
scatter(f.ax(1),relacont(:,1),relacont(:,2),20,'k','filled','o');
colormap(f.ax(1),'jet');
c=colorbar(f.ax(1),'NorthOutside');
c.Label.String = 'Day in Sept. 2005';
c.Label.FontSize = 12;
caxis(f.ax(1),cran);
plot(f.ax(1),[-100 100],[0 0],'k--');
plot(f.ax(1),[0 0],[-100 100],'k--');
text(f.ax(1),0.85,0.9,'HF','FontSize',15,'unit','normalized');
text(f.ax(1),0.05,0.9,'(a)','FontSize',12,'unit','normalized');
f.ax(1).Box = 'on';
grid(f.ax(1), 'on');
f.ax(1).GridLineStyle = '--';
axis(f.ax(1), 'equal');
axis(f.ax(1),[-20 35 -25 20]);
f.ax(1).XLabel.String = 'E (km)';
f.ax(1).YLabel.String = 'N (km)';
hold(f.ax(1),'off');

% subplot 2 of figure 1
hold(f.ax(2),'on');
scatter(f.ax(2),det1(:,6),det1(:,7), msizehf, det1(:,end), 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(f.ax(2),relacont(:,1),relacont(:,2),20,'k','filled','o');
colormap(f.ax(2),'jet');
c=colorbar(f.ax(2),'NorthOutside');
c.Label.String = 'Day in Sept. 2005';
c.Label.FontSize = 12;
caxis(f.ax(2),cran);
plot(f.ax(2),[-100 100],[0 0],'k--');
plot(f.ax(2),[0 0],[-100 100],'k--');
text(f.ax(2),0.85,0.9,'HF','FontSize',15,'unit','normalized');
text(f.ax(2),0.05,0.9,'(b)','FontSize',12,'unit','normalized');
f.ax(2).Box = 'on';
grid(f.ax(2), 'on');
f.ax(2).GridLineStyle = '--';
axis(f.ax(2), 'equal');
axis(f.ax(2),[-20 35 -25 20]);
f.ax(2).XLabel.String = 'E (km)';
f.ax(2).YLabel.String = 'N (km)';
hold(f.ax(2),'off');


print(f.fig,'-dpdf',strcat(rstpath,'/yajun.',num2str(refyr),'.timerela.pdf'));
  
  
% det2 = det2014(det2014(:,2)==195,:);

% figure
% dum2 = sortrows(det2,-3);
% scatter(dum2(:,6),dum2(:,7), 7, dum2(:,3)/3600, 'filled','o');
% caxis([0 24]);