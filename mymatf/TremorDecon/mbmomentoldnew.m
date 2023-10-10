% lfecatmoment.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Micheal Bostock has updated his moment estimate for the LFE catalog. This 
% code aims to compare the moment of the new and old catalog.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/09/16
% Last modified date:   2023/09/16
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

[scrsz, resol] = pixelperinch(1);

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');

bosdir = ('/home/data2/chaosong/matlab/allan/BOSTOCK');
lfefnm = ('newlfeloc');
lfeloc = load(fullfile(bosdir, lfefnm));
loc0 = lfeloc(lfeloc(:,1)==2,:);
%%%if still use the catalog separated by each fam
bostname = '/BOSTOCK/update20230916/total_mag_detect_0000_cull_NEW.txt'
bostcat = load(strcat(workpath,bostname));
% bostcat = ReformBostock(loc0(3),loc0(2),1,bostname);
% bostcat = bostcat(bostcat(:,2)==2003 | bostcat(:,2)==2004 | bostcat(:,2)==2005, :);
% bostcat = bostcat(bostcat(:,1)~=2 & bostcat(:,1)~=47 & bostcat(:,1)~=246, :);
bostcat(:,end+1) = mag2moment(bostcat(:,5));

bostnamenew = '/BOSTOCK/update20230916/total_mag_detect_0000_cull_AMR.txt'
bostcatnew = load(strcat(workpath,bostnamenew));
% bostcatnew = ReformBostock(loc0(3),loc0(2),1,bostnamenew);
% bostcatnew = bostcatnew(bostcatnew(:,2)==2003 | bostcatnew(:,2)==2004 | bostcatnew(:,2)==2005, :);
% bostcatnew = bostcatnew(bostcatnew(:,1)~=2 & bostcatnew(:,1)~=47 & bostcatnew(:,1)~=246, :);
bostcatnew(:,end+1) = mag2moment(bostcatnew(:,5));


%% cross plot of old and new moments
widin = 5;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
nrow = 1;
ncol = 1;
f = initfig(widin,htin,nrow,ncol); %initialize fig
  
fttpfree = fittype( @(a,b,x) a*x+b);

ax=f.ax(1);
hold(ax,'on'); ax.Box = 'on'; grid(ax,'on');
x = log10(bostcat(:,end));
y = log10(bostcatnew(:,end));
scatter(ax,x,y,15,'r','filled');
%linear robust least square
[fitobj,gof,output] = fit(x,y,fttpfree,'Robust','Bisquare','StartPoint',[1 1]);
% [fitobj,gof,output] = fit(x,y,'poly1');
% stats = statsofrobustlnfit(fitobj,gof,output,x,y);
coef = coeffvalues(fitobj);
slope = coef(1);
intcpt = coef(2);
fitx = linspace(9.5,13,100);
fity = feval(fitobj,fitx);
plot(ax,fitx,fity,'-','linewidth',2,'color','k');
yeqx = fitx;
plot(ax,fitx,yeqx,'--','linewidth',1,'color','k');
text(ax,0.95,0.95,sprintf('y = %.2f x + %.2f',slope,intcpt),'HorizontalAlignment',...
  'right','Units','normalized');
% text(ax,0.95,0.90,sprintf('weighted pearson = %.2f',stats.pearwt),'HorizontalAlignment',...
%   'right','Units','normalized');
% text(ax,0.95,0.25,sprintf('In burst: %d; %.2f; %.2e',nbostia,rbostia,bostiams),'HorizontalAlignment',...
%   'right','Units','normalized');
% text(ax,0.95,0.2,sprintf('Outside: %d; %.2f; %.2e',nbostimis,rbostimis,bostimisms),'HorizontalAlignment',...
%   'right','Units','normalized');
% text(ax,0.95,0.05,sprintf('%d/%d in common',length(x),size(bosta,1)),'HorizontalAlignment',...
%   'right','Units','normalized');
xlabel(ax,"log_{10}{moment from NEW.txt (NM)}");
ylabel(ax,"log_{10}{moment from AMR.txt (NM)}");
axis(ax,'equal');
axis(ax,[9.5 13 9.5 13]);


