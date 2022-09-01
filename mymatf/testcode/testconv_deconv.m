% testconv_deconv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is to test the algorithm of convolution, to see if the stacking of
% templates directly could be achieved by convolution as well
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/07/16
% Last modified date:   2021/07/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e   % Set the format to 5-digit floating point
clc
clear
close all


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
segtemp1 = tracebest(tst1:ted1, :);
segtemp2 = tracebest(tst2:ted2, :);

%%
sps = 20;
scale = 2;
len = 4*sps;
temp = tracetemp(size(tracetemp,1)/2-len/2+1: size(tracetemp,1)/2+len/2, 3);

amp = zeros(size(temp,1),1);
ind = [len*1/4+1, len/2+1,len*3.5/4+1];
% ind = [len/2+1];
% amp(ind) = 0.5;
amp(ind) = [0.6 0.5 0.7];

rst = conv(temp,amp,'same');

figure
subplot(411)
hold on
box on
grid on
plot(1:size(temp,1),temp,'k','linew',2);
xlim([0 len]);

subplot(412)
hold on
box on
grid on
stem((1:size(temp,1))-len/2-1,amp,'b','linew',1);
xlim([-len/2 len/2]);

subplot(413)
hold on
box on
grid on
plot(1:size(rst,1),rst,'k','linew',2);
% ylim([-1,1]);
xlim([0 len]);

%%%%%%%% use online code with spectral division, works, but not perfect %%%%%%%%
[t,ampinv] = deconvolution(rst,temp);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(414)
hold on
box on
grid on
stem(t,ampinv,'b','linew',1);
xlim([-len/2 len/2]);


%%
temp = tracetemp(:, 3);
best = tracebest(:,3);
len = size(temp,1);

bst2 = 301;
bed2 = 330;
durab = bed2-bst2+1;
durat = ted2-tst2+1;
duradif = durab - durat;

amp = zeros(len,1);
% SPACING = 'eqampeqspa';
SPACING = 'n7';
if isequal(SPACING, 'eqampeqspa')
    % ind = [len/2+1];
    % ind = [bst2]+(len/2+1-tst2);
    ind = (0:1:8)+(bst2+len/2+1-tst2);   % index of non-zero
    ampind = ones(length(ind),1);
    for i = 1: length(ind)
        amp(ind(i)) = ampind(i);
    end
elseif isequal(SPACING, 'n7')
    ind = [0 2 3 4 5 6 8]+(bst2+len/2+1-tst2);
    ampind = [0.1 0.4 1.6 6.4 1.6 0.4 0.1]/1.5;
    for i = 1: length(ind)
        amp(ind(i)) = ampind(i);
    end
end

rst = conv(temp,amp,'same');
% rst = conv(amp,temp);

figure
subplot(411)
hold on
box on
grid on
plot(1:length(best),best,'r','linew',2);
plot(1:length(temp),temp,'k','linew',2);
xlim([0 len]);

subplot(412)
hold on
box on
grid on
stem((1:length(amp))-len/2-1,amp,'b','linew',1);
xlim([-len/2 len/2]);

subplot(413)
hold on
box on
grid on
plot(1:length(rst),rst,'k','linew',2);
xlim([0 len]);

% vertical_cursors


%%%%%%%% use matlab built-in 'deconv', NOT working %%%%%%%%
% ampinv = zeros(len,1);
% ampdec = deconv(rst,temp);
% 
% [ampdec, res] = deconv(best,temp);
%     for i = 1: length(ind)
%         ampinv(ind(i)) = ampdec(i);
%     end
% ampinv = deconv(rst,temp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%% use matlab built-in 'deconvblind', NOT working %%%%%%%%
% %%% the problem is the inverted amplitude is offsetted in time, idk why should it be
% [ampinv,tempr] = deconvblind(rst,temp,20);
% subplot(414)
% hold on
% box on
% grid on
% stem((1:length(ampinv)),ampinv,'b','linew',1);
% stem((1:length(tempr)),tempr,'r','linew',1);
% xlim([0 len]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% use online code with spectral division, works, but not perfect %%%%%%%%
[t,ampinv] = deconvolution(best,temp);  
subplot(414)
hold on
box on
grid on
stem(t,ampinv,'b','linew',1);
xlim([-len/2 len/2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









