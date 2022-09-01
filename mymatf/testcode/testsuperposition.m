%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test the forward modeling of the superposition of several dipoles exactly from the templates that
% are separated in time, to see if the example we have seen that has a short time lag from the min
% to max is short compared to the full duration.
%
% But first we need to examine the duration of the templates using the way we do to the best
% example, to decide the spacing of the component templates.
% Second, we need to examine the min-to-max time lag of the templates, and compare it to the best
% example, to see if the templates need to shrink/stretch to match the time lag of the example.
%
% We can view this test as the forward modeling of the deconvolution. If the result is not intuitive
% enough, we may have to do the deconvolution directly.
%
% -- This part of the code is the following code of 'testcutmaindipole.m'. It will read in the
%   seismogram, do stretching orshrinking if necessary, and then try to superimpose several cycles 
%   of templates that are separated in time in some manners, to see if the example detection could
%   be simulated by this kind of superposition.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/07/15
% Last modified date:   2021/07/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e   % Set the format to 5-digit floating point
clc
clear
% close all


%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
%%% choose the data path here
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');
figpath = strcat(workpath, '/project2021/PGCtrio');
hypopath = strcat(getenv('MHOME'), '/Seisbasics/hypoinverse/forproj21');

fam = '243';

FLAG = 'PGC';

iup = 1;

timoffrot= [
    2003 061;
    2003 062; %Could leave out 064 and 257.       % dates of the year, two columns: year Julian_day
    2003 063;
    2004 196;
    2004 197;
    2004 198;
    2004 199;
    2005 254;
    2005 255;
    2005 256;
    ];
nday = size(timoffrot, 1);

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
    'LZB'];
POLSTA=['SSIB '           % polaris station names
    'SILB '
    'KLNB '
    'MGCB '
    'TWKB '];

stas=['PGC '
    'SSIB'
    'SILB'];     % determine the trio and order
nsta=size(stas,1);         %  number of stations

%load the saved file of template and best example
date = 2004199;
nwinlen = 32*20;
fname = strcat(rstpath, '/MAPS/bestexample_',fam,'_',num2str(date),'_',num2str(nwinlen/20),...
    's',num2str(20),'sps');
tracebest = load(fname);

wlen = 32*40;
fname = strcat(rstpath, '/MAPS/templateseg_',fam,'_',num2str(wlen/40),...
    's',num2str(40),'sps');
tracetemp = load(fname);

%resample the template to the same sps, 20
[num, denom] = rat(20/40);
tracetemp = resample(tracetemp,num,denom);   % times = num/denom

bst1 = 310;
bed1 = 327;
bst2 = 301;
bed2 = 330;
bzc = 320;
segbest1 = tracebest(bst1:bed1, :);
segbest2 = tracebest(bst2:bed2, :);

tst1 = 312;
ted1 = 326;
tst2 = 309;
ted2 = 330;
tzc = 320;
segtemp1 = tracetemp(tst1:ted1, :);
segtemp2 = tracetemp(tst2:ted2, :);

%%
figure
subplot(211)
ax = gca;
hold(ax,'on');
box on
grid on
plot(ax,[0 700],[0 0],'k--','linew',1);
plot(ax,tracebest(:, 1),'r','linew',1);
plot(ax,tracebest(:, 2),'b','linew',1);
plot(ax,tracebest(:, 3),'k','linew',1);
plot(ax,[bst1 bst1],ax.YLim,'k--','linew',0.5);
plot(ax,[bed1 bed1],ax.YLim,'k--','linew',0.5);
plot(ax,[bst2 bst2],ax.YLim,'k--','linew',0.5);
plot(ax,[bed2 bed2],ax.YLim,'k--','linew',0.5);
ylabel(ax,'Amplitude','fontsize',9);
hold(ax,'off');


subplot(212)
ax = gca;
hold(ax,'on');
box on
grid on
plot(ax,[0 700],[0 0],'k--','linew',1);
plot(ax,tracetemp(:, 1),'r','linew',1);
plot(ax,tracetemp(:, 2),'b','linew',1);
plot(ax,tracetemp(:, 3),'k','linew',1);
plot(ax,[tst1 tst1],ax.YLim,'k--','linew',0.5);
plot(ax,[ted1 ted1],ax.YLim,'k--','linew',0.5);
plot(ax,[tst2 tst2],ax.YLim,'k--','linew',0.5);
plot(ax,[ted2 ted2],ax.YLim,'k--','linew',0.5);
ylabel(ax,'Amplitude','fontsize',9);
hold(ax,'off');

vertical_cursors

%% try different ways to compose the synthetics
%%%%%% way 1, use the original trace %%%%%%%
%%% stretch/shrink the record if necessary
sps = 20;
scale = 2;
len = 4*sps;
test = tracetemp(size(tracetemp,1)/2-len/2+1: size(tracetemp,1)/2+len/2, 1);
samp = 1:1:len;
sampstr = 1/scale:1/scale:len;
teststr = interp1(samp, test, sampstr, 'spline');
figure
hold on
plot(samp-len/2, test, 'k','linew',2);
plot(sampstr*scale-len*scale/2, teststr, 'b','linew',1);

vel = test;
acc = diff(vel);
figure
hold on
plot(vel, 'k','linew',2);
plot(acc, 'b','linew',1);


% % shift the trace, if the shifted amount is positive, then is shifting to the right
% st = 10;
% ed = len-10;
% sft = 4;
% test1 = test(st: ed);
% test1sft = test(st-sft: ed-sft);
% figure
% hold on
% plot(test1, 'k','linew',2);
% plot(test1sft, 'b','linew',2);

% %according to the automated estimate of the min-to-max time lag, we need to strech the template
% sps = 20;
% scale = 0.2/0.15;
% len = size(tracetemp,1);
% samp = 1:1:len;
% sampstr = 1/scale:1/scale:len;
% tracetempstr = interp1(samp, tracetemp, sampstr, 'spline');

% shift the template, if the shifted amount is negative, then is shifting to the right
sps = 20;
sft = tst2 - bst2;
st = 2*sps+1;
ed = size(tracetemp,1)-2*sps;
tempcn = tracetemp(st: ed, 3);
temprt = tracetemp(st-sft: ed-sft, 3);
templt = tracetemp(st+sft: ed+sft, 3);

durab = bed2-bst2+1;
durat = ted2-tst2+1;
duradif = durab - durat;

%if the spacing is equal
% SPACING = 'eqampeqspa';
% SPACING = 'uneqampeqspa';
% SPACING = 'eqampuneqspa';
% SPACING = 'uneqampuneqspa';
SPACING = 'n3';
% SPACING = 'n5';
% SPACING = 'n7';


if isequal(SPACING, 'eqampeqspa')
    ntemp = 9;
    istf = round(duradif/(ntemp-1));
    itemp = [];
    ampsca = ones(ntemp,1);
    for i = 1: ntemp
        itemp(:,i) = ampsca(i).*tracetemp(st+sft-istf*(i-1): ed+sft-istf*(i-1), 3);
    end
    syn = sum(itemp,2);

elseif isequal(SPACING, 'uneqampeqspa')
    ntemp = 9;
    istf = round(duradif/(ntemp-1));
    itemp = [];
    ampsca = [1:1:5 fliplr(1:1:4)];
    for i = 1: ntemp
        itemp(:,i) = ampsca(i).*tracetemp(st+sft-istf*(i-1): ed+sft-istf*(i-1), 3);
    end
    syn = sum(itemp,2);
    
elseif isequal(SPACING, 'eqampuneqspa')
    ntemp = 7;
    ampsca = ones(ntemp,1);
    itemp(:,1) = ampsca(1).*tracetemp(st+sft-0: ed+sft-0, 3);
    itemp(:,2) = ampsca(2).*tracetemp(st+sft-2: ed+sft-2, 3);
    itemp(:,3) = ampsca(3).*tracetemp(st+sft-3: ed+sft-3, 3);
    itemp(:,4) = ampsca(4).*tracetemp(st+sft-4: ed+sft-4, 3);
    itemp(:,5) = ampsca(5).*tracetemp(st+sft-(duradif-3): ed+sft-(duradif-3), 3);
    itemp(:,6) = ampsca(6).*tracetemp(st+sft-(duradif-2): ed+sft-(duradif-2), 3);
    itemp(:,7) = ampsca(7).*tracetemp(st+sft-duradif: ed+sft-duradif, 3);
    syn = sum(itemp,2);
    
elseif isequal(SPACING, 'uneqampuneqspa')
    ntemp = 7;
    ampsca = [1:1:4 fliplr(1:1:3)];
    itemp(:,1) = ampsca(1).*tracetemp(st+sft-0: ed+sft-0, 3);
    itemp(:,2) = ampsca(2).*tracetemp(st+sft-2: ed+sft-2, 3);
    itemp(:,3) = ampsca(3).*tracetemp(st+sft-3: ed+sft-3, 3);
    itemp(:,4) = ampsca(4).*tracetemp(st+sft-4: ed+sft-4, 3);
    itemp(:,5) = ampsca(5).*tracetemp(st+sft-(duradif-3): ed+sft-(duradif-3), 3);
    itemp(:,6) = ampsca(6).*tracetemp(st+sft-(duradif-2): ed+sft-(duradif-2), 3);
    itemp(:,7) = ampsca(7).*tracetemp(st+sft-duradif: ed+sft-duradif, 3);
    syn = sum(itemp,2);

elseif isequal(SPACING, 'n7')
    ntemp = 7;
%     ampsca = [0.1 0.4 1.6 6.4 1.6 0.4 0.1];
    ampsca = [0.1 0.4 1.6 6.4 1.6 0.4 0.1]/1.5;
%     ampsca = [0.1 0.5 2.5 12.5 2.5 0.5 0.1];
    itemp(:,1) = ampsca(1).*tracetemp(st+sft-0: ed+sft-0, 3);
    itemp(:,2) = ampsca(2).*tracetemp(st+sft-2: ed+sft-2, 3);
    itemp(:,3) = ampsca(3).*tracetemp(st+sft-3: ed+sft-3, 3);
    itemp(:,4) = ampsca(4).*tracetemp(st+sft-4: ed+sft-4, 3);
    itemp(:,5) = ampsca(5).*tracetemp(st+sft-(duradif-3): ed+sft-(duradif-3), 3);
    itemp(:,6) = ampsca(6).*tracetemp(st+sft-(duradif-2): ed+sft-(duradif-2), 3);
    itemp(:,7) = ampsca(7).*tracetemp(st+sft-duradif: ed+sft-duradif, 3);
    syn = sum(itemp,2);
    
elseif isequal(SPACING, 'n5')
    ntemp = 5;
%     ampsca = [0.1 0.4 1.6 0.4 0.1];
    ampsca = [0.1 0.4 1.6 0.4 0.1]*3;
%     ampsca = [0.1 0.5 2.5 0.5 0.1]*3;
    itemp(:,1) = ampsca(1).*tracetemp(st+sft-0: ed+sft-0, 3);
    itemp(:,2) = ampsca(2).*tracetemp(st+sft-3: ed+sft-3, 3);
    itemp(:,3) = ampsca(3).*tracetemp(st+sft-4: ed+sft-4, 3);
    itemp(:,4) = ampsca(4).*tracetemp(st+sft-(duradif-3): ed+sft-(duradif-3), 3);
    itemp(:,5) = ampsca(5).*tracetemp(st+sft-duradif: ed+sft-duradif, 3);
    syn = sum(itemp,2);
    
elseif isequal(SPACING, 'n3')
    ntemp = 3;
    ampsca = [0.5 5 0.5];
    itemp(:,1) = ampsca(1).*tracetemp(st+sft-0: ed+sft-0, 3);
    itemp(:,2) = ampsca(2).*tracetemp(st+sft-4: ed+sft-4, 3);
    itemp(:,3) = ampsca(3).*tracetemp(st+sft-duradif: ed+sft-duradif, 3);
    syn = sum(itemp,2);
    
end



best = tracebest(st: ed, 3);

figure
subplot(211)
hold on
box on
grid on
plot(best,'k','linew',2);

subplot(212)
hold on
box on
grid on
% plot(tempcn, 'k','linew',2);
% plot(temprt, 'b','linew',2);
% plot(templt, 'r','linew',2);
for i = 1: ntemp
    plot(itemp(:,i),'linew',1);
end
plot(syn,'k','linew',2);

vertical_cursors


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% way 2, use the trace segment only %%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% way 2, use the trace segment only %%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














