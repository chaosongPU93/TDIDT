% function plot_all_rela
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to plot the locations of all detections inverted by 
% hypoinverse in km relative to offset (0,0)
% 
% Looks like figure 3b in Rubin&Armbruster (2013)
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/08/20
% Last modified date:   2019/08/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/forsummary');

fam = '002';
winlenhf = 4;

winlenlf = 16;
lofflf = 4;
ccminlf = 0.5;


SUFFIXhf = strcat('hf.alldays.count.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/eventloc.',fam,'.',SUFFIXhf);
hfmaplocall = load(fname);
% 8 cols, format is:
%   lon lat dep off12 off13 off12sec off13sec count

SUFFIXlf = strcat('lf.alldays.count.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/eventloc.',fam,'.',SUFFIXlf);
lfmaplocall = load(fname);

[ind,~] = ind2sub(size(hfmaplocall), find(hfmaplocall(:,4)==0 & hfmaplocall(:,5)==0, 1, 'first'));
lon0 = hfmaplocall(ind,1);
lat0 = hfmaplocall(ind,2);
if isempty(lon0)
    lon0 = -123.5850;
    lat0 = 48.4367;
end

hfrelalocall = hfmaplocall;
[dx, dy] = absloc2relaloc(hfmaplocall(:,1),hfmaplocall(:,2),lon0,lat0);

hfrelalocall(:,1) = dx;
hfrelalocall(:,2) = dy;

lfrelalocall = lfmaplocall;
[dx, dy] = absloc2relaloc(lfmaplocall(:,1),lfmaplocall(:,2),lon0,lat0);

lfrelalocall(:,1) = dx;
lfrelalocall(:,2) = dy;

%% plot the accumulative hit count
%%% figure 1, plot the overall hf and lf detections
% define and position the figure frame and axes of each plot 
f.fig = figure(1);
widin = 8;  % maximum width allowed is 8.5 inches
htin = 4.5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end 

% reposition
for isub = 1:nrow*ncol
    if mod(isub,ncol)==0
        icol = ncol;
        irow = floor(isub/ncol);
    else
        icol = mod(isub,ncol);
        irow = floor(isub/ncol)+1;
    end
    set(f.ax(isub),'Position', [(icol-0.85)/ncol, 1-irow/nrow*0.9, 1/ncol*0.8 1/nrow*0.8]);
    f.ax(isub).Box = 'on';
    grid(f.ax(isub), 'on');
    f.ax(isub).GridLineStyle = '--';
end

% marker size
msizehf = 6;
msizelf = 12;

% subplot 1 of figure 1
hold(f.ax(1),'on');
dumhf = hfrelalocall;
dumhf(dumhf(:,8)>1, :) = [];
scatter(f.ax(1),dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,8)),'o');  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(hfrelalocall,8);
dumhf(dumhf(:,8)==1, :) = [];
scatter(f.ax(1),dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,8)), 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(f.ax(1),0,0,20,'k','o','LineWidth',1);
plot(f.ax(1),[-100 100],[0 0],'k--');
plot(f.ax(1),[0 0],[-100 100],'k--');
text(f.ax(1),0.05,0.93,'a','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f.ax(1),0.91,0.93,'PGC','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f.ax(1),0.91,0.1,'HF','FontSize',15,'unit','normalized','horizontalalignment','center');
colormap(f.ax(1),'jet');
c=colorbar(f.ax(1),'NorthOutside');
c.Label.String = 'log_{10}(hits)';
c.Label.FontSize = 12;
caxis(f.ax(1),[0 2]);
axis(f.ax(1), 'equal');
axis(f.ax(1),[-30 20 -30 25]);
f.ax(1).XLabel.String = 'E (km)';
f.ax(1).YLabel.String = 'N (km)';
hold(f.ax(1),'off');

% subplot 2 of figure 1
hold(f.ax(2),'on');
dumlf = lfrelalocall;
dumlf(dumlf(:,8)>1, :) = [];
scatter(f.ax(2),dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,8)),'o');   %, 'MarkerEdgeColor', 'w')
dumlf = sortrows(lfrelalocall,8);
dumlf(dumlf(:,8)==1, :) = [];
scatter(f.ax(2),dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,8)), 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(f.ax(2),0,0,20,'k','o','LineWidth',1);
plot(f.ax(2),[-100 100],[0 0],'k--');
plot(f.ax(2),[0 0],[-100 100],'k--');
text(f.ax(2),0.05,0.93,'b','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f.ax(2),0.91,0.93,'PGC','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f.ax(2),0.91,0.1,'LF','FontSize',15,'unit','normalized','horizontalalignment','center');
colormap(f.ax(2),'jet');
c=colorbar(f.ax(2),'NorthOutside');
c.Label.String = 'log_{10}(hits)';
c.Label.FontSize = 12;
caxis(f.ax(2),[0 1]);
axis(f.ax(2), 'equal');
axis(f.ax(2),[-30 20 -30 25]);
f.ax(2).XLabel.String = 'E (km)';
f.ax(2).YLabel.String = 'N (km)';
hold(f.ax(2),'off');


%%% save figure to file
print(f.fig,'-dpdf',strcat(rstpath,'/',fam,'.hflf.alld.accuhit',num2str(winlenhf),'_',num2str(winlenlf),...
      '.',num2str(lofflf),'.',num2str(ccminlf),'.rela.pdf'));
% keyboard



%% plot the time evolution 
SUFFIXhf = strcat('hf.alldays.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/eventloc.',fam,'.',SUFFIXhf);
hfmaploctime = load(fname);
% 10 cols, format is:
%   lon lat dep off12 off13 off12sec off13sec date strongest_arr_time center_of_win


SUFFIXlf = strcat('lf.alldays.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/eventloc.',fam,'.',SUFFIXlf);
lfmaploctime = load(fname);

% convert absolute loc to relative loc
[ind,~] = ind2sub(size(hfmaploctime), find(hfmaploctime(:,4)==0 & hfmaploctime(:,5)==0, 1, 'first'));
lon0 = hfmaploctime(ind,1);
lat0 = hfmaploctime(ind,2);
if isempty(lon0)
    lon0 = -123.5850;
    lat0 = 48.4367;
end

hfrelaloctime = hfmaploctime;
[dx, dy] = absloc2relaloc(hfmaploctime(:,1),hfmaploctime(:,2),lon0,lat0);
hfrelaloctime(:,1) = dx;
hfrelaloctime(:,2) = dy;

lfrelaloctime = lfmaploctime;
[dx, dy] = absloc2relaloc(lfmaploctime(:,1),lfmaploctime(:,2),lon0,lat0);
lfrelaloctime(:,1) = dx;
lfrelaloctime(:,2) = dy;


%%% select a year
refyr = 2005;
refjday = 254;
tmp = jul2dat(refyr,refjday);   % return month, day, year
refday = tmp(2);    
daymax = 14;
selhf = hfrelaloctime(hfrelaloctime(:,8)>=refyr*1000,:);
timeint = selhf(:,8)-(refyr*1000+refjday)+refday;
timedec = selhf(:,9)/3600/24;
time = timeint+timedec;
selhf(:,11) = time;

sellf = lfrelaloctime(lfrelaloctime(:,8)>=refyr*1000,:);
timeint = sellf(:,8)-(refyr*1000+refjday)+refday;
timedec = sellf(:,9)/3600/24;
time = timeint+timedec;
sellf(:,11) = time;

%%% define and position the figure frame and axes of each plot 
f2.fig = figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f2.ax(isub) = subplot(nrow,ncol,isub);
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
    set(f2.ax(isub),'Position', [(icol-0.9)/ncol, 1-irow/nrow*0.95, 1/ncol*0.8 1/nrow*0.8]);
end

% marker size
msizehf = 4;
msizelf = 8;

%%% deal with each subplot
% subplot 1 of figure 1
hold(f2.ax(1),'on');
dumhf = sortrows(selhf,-9);
scatter(f2.ax(1),dumhf(:,1),dumhf(:,2), msizehf, dumhf(:,end), 'filled','o');
scatter(f2.ax(1),0,0,20,'k','o','LineWidth',1);
colormap(f2.ax(1),'jet');
c=colorbar(f2.ax(1),'NorthOutside');
c.Label.String = 'Day in Sept. 2005';
c.Label.FontSize = 12;
caxis(f2.ax(1),[refday daymax]);
plot(f2.ax(1),[-100 100],[0 0],'k--');
plot(f2.ax(1),[0 0],[-100 100],'k--');
text(f2.ax(1),0.05,0.93,'a','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f2.ax(1),0.91,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f2.ax(1),0.91,0.1,'HF','FontSize',15,'unit','normalized','horizontalalignment','center');
f2.ax(1).Box = 'on';
grid(f2.ax(1), 'on');
f2.ax(1).GridLineStyle = '--';
axis(f2.ax(1), 'equal');
axis(f2.ax(1),[-30 20 -30 20]);
f2.ax(1).XLabel.String = 'E (km)';
f2.ax(1).YLabel.String = 'N (km)';
hold(f2.ax(1),'off');

% subplot 2 of figure 1
hold(f2.ax(2),'on');
scatter(f2.ax(2),selhf(:,1),selhf(:,2), msizehf, selhf(:,end), 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(f2.ax(2),0,0,20,'k','o','LineWidth',1);
colormap(f2.ax(2),'jet');
c=colorbar(f2.ax(2),'NorthOutside');
c.Label.String = 'Day in Sept. 2005';
c.Label.FontSize = 12;
caxis(f2.ax(2),[refday daymax]);
plot(f2.ax(2),[-100 100],[0 0],'k--');
plot(f2.ax(2),[0 0],[-100 100],'k--');
text(f2.ax(2),0.05,0.93,'b','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f2.ax(2),0.91,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f2.ax(2),0.91,0.1,'HF','FontSize',15,'unit','normalized','horizontalalignment','center');
f2.ax(2).Box = 'on';
grid(f2.ax(2), 'on');
f2.ax(2).GridLineStyle = '--';
axis(f2.ax(2), 'equal');
axis(f2.ax(2),[-30 20 -30 20]);
f2.ax(2).XLabel.String = 'E (km)';
f2.ax(2).YLabel.String = 'N (km)';
hold(f2.ax(2),'off');

% subplot 3 of figure 1
hold(f2.ax(3),'on');
dumlf = sortrows(sellf,-9);
scatter(f2.ax(3),dumlf(:,1),dumlf(:,2), msizelf, dumlf(:,end), 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(f2.ax(3),0,0,20,'k','o','LineWidth',1);
colormap(f2.ax(3),'jet');
colorbar(f2.ax(3),'NorthOutside');
caxis(f2.ax(3),[refday daymax]);
plot(f2.ax(3),[-100 100],[0 0],'k--');
plot(f2.ax(3),[0 0],[-100 100],'k--');
text(f2.ax(3),0.05,0.93,'c','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f2.ax(3),0.91,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f2.ax(3),0.91,0.1,'LF','FontSize',15,'unit','normalized','horizontalalignment','center');
f2.ax(3).Box = 'on';
grid(f2.ax(3), 'on');
f2.ax(3).GridLineStyle = '--';
axis(f2.ax(3), 'equal');
axis(f2.ax(3),[-30 20 -30 20]);
f2.ax(3).XLabel.String = 'E (km)';
f2.ax(3).YLabel.String = 'N (km)';
hold(f2.ax(3),'off');

% subplot 4 of figure 1
hold(f2.ax(4),'on');
scatter(f2.ax(4),sellf(:,1),sellf(:,2), msizelf, sellf(:,end), 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(f2.ax(4),0,0,20,'k','o','LineWidth',1);
colormap(f2.ax(4),'jet');
colorbar(f2.ax(4),'NorthOutside');
caxis(f2.ax(4),[refday daymax]);
plot(f2.ax(4),[-100 100],[0 0],'k--');
plot(f2.ax(4),[0 0],[-100 100],'k--');
text(f2.ax(4),0.05,0.93,'d','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f2.ax(4),0.91,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f2.ax(4),0.91,0.1,'LF','FontSize',15,'unit','normalized','horizontalalignment','center');
f2.ax(4).Box = 'on';
grid(f2.ax(4), 'on');
f2.ax(4).GridLineStyle = '--';
axis(f2.ax(4), 'equal');
axis(f2.ax(4),[-30 20 -30 20]);
f2.ax(4).XLabel.String = 'E (km)';
f2.ax(4).YLabel.String = 'N (km)';
hold(f2.ax(4),'off');

    
%%% save figure
print(f2.fig,'-dpdf',strcat(rstpath,'/',fam,'hflf.alld.tevo',num2str(refyr),'.',num2str(winlenhf),'_',...
      num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.rela.pdf'));
% close(f2.fig)





















