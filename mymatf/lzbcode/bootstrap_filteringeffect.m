%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to bootstrap the filtering effect correction results 
% from LZB trio.
% i.e., 
%   1. preturb the detection in samples of time assuming there is error in
%       detection;
%   2. preturb the filtering correction in samples of time (or magnitude and
%       azimuth?? i am not sure yet) assuming there is also error in passband
%       correction;
%   *** the above 2 preturbation assumes that there is error in each detection
%       but does not change the number of detection.
%
%   3. preturb the current selection of detections (and/or the selection of 
%       migrations), i.e., randomly choose a subset of detections or migrations
%       to carry out the same analysis to obtain the summary results, say 1/3,
%       assuming there is potentially bias in selecting detections, or bias in 
%       choosing thresholds in detection. Actually this 3rd perturbation is
%       essential because the thresholds are totally impirical (subjective), and
%       bootstrapping could mitigate this effect.
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/08/27
% Last modified date:   2020/08/27

%% Initialization
format short e   % Set the format to 5-digit floating point
clc
clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not


%% PGC trio
%%% fc is inverted with (lf-hf_12, lf-hf_13) relative to fam 002
fcpgc = [-123.592000 48.436500 36.7900];

%%% this is inverted from (0,0) relative to fam 002, location of control points
contpgc = [-123.585000 48.436667 36.8800];

%%% convert to km
[dx, dy] = absloc2relaloc(fcpgc(1),fcpgc(2),contpgc(1),contpgc(2));
fcvecpgc = [dx dy];
fcmagpgc = sqrt(dx.^2+dy.^2);


%% LZB trio

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
nfam = length(nfampool);

%%% fc is inverted with (lf-hf_12, lf-hf_13) relative to all fams
%%% the same content as /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/LZBallfamfc, which is
%%% inverted by hypoinverse from /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/chat00p.fc
fclzb = [-123.509500 48.457000 38.0200;
    -123.801000 48.561333 36.2500;
    -123.896333 48.594000 35.7800;
    -123.638833 48.474000 36.7200;
    -123.797167 48.440333 34.8100;
    -123.925000 48.599333 35.5600;
    -123.898667 48.555833 35.2500;
    -123.772167 48.575333 36.7300;
    -123.734833 48.562667 36.8900;
    -123.837500 48.587500 36.3200;
    -123.867833 48.590000 36.0100];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VERSION 1 for the locations of lfe families, by inversion from hypoinverse
%%% this is inverted from (0,0) of all fams, same order, location of control points
%%% inverted from /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/chat00p.reffam
contlzb = [-123.492667 48.451500 38.1400;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% convert each filterring correction to km relative to its own fam
fcveclzb = zeros(size(fclzb,1),2);
fcmaglzb = zeros(size(fclzb,1),1);
for i = 1: size(fclzb,1)
    [dx, dy] = absloc2relaloc(fclzb(i,1),fclzb(i,2),contlzb(i,1),contlzb(i,2));
    fcveclzb(i,:) = [dx dy];
    fcmaglzb(i) = sqrt(dx.^2+dy.^2);
end

%%% convert lfe location to km relative to 043
loclfe = zeros(size(contlzb,1),2);
for i = 1: size(contlzb,1)
    [dx, dy] = absloc2relaloc(contlzb(i,1),contlzb(i,2),contlzb(2,1),contlzb(2,2));
    loclfe(i,:) = [dx dy];
end

%%% get the mean and standard deviation of current correction
% exclude those farther and inconsistent fams, 010, 002, 047
indexc = [1; 4; 5];

tempvec = fcveclzb;
tempmag = fcmaglzb;

tempvec(indexc,:)=[];
tempmag(indexc,:)=[];

% use the correction in km as estimate of std in x and y direction 
vecmu = mean(tempvec,1);
vecstd = std(tempvec,1);

magmu = mean(tempmag);
magstd = std(tempmag);


% Take the filtering effect correction of each fam as the 'mean' or mu of the true distribution of
% the true value, and std of all fams as the common 'std' or sigma of the true distribution of the
% true value. Then we resample this distribution multiple times to obtain new filtering correction.

% 

%%% WAY 1
% so we could use the following codes to create a normal/Gaussian distribution
% with 'makedist' to create a prob. distri. object by specifying param. values
pd = makedist('Normal','mu',0,'sigma',1);
rng('default');
r1 = random(pd, [100,1]);

%%% WAY 2
% Alternatively, you could directly generate random numbers from normal distribution 
rng('default');
r2 = random('Normal', 0, 1, [100, 1]); % 1000 is number of samples   

%%% WAY 3
% we could use 'normrnd' to generate normal random numbers 
rng('default');  % save the current state of random number generator for reproducibility
r3 = normrnd(0,1,size(x));

s = rng; 
r4 = normrnd(0,1,[1,5]);

% use the saved the state of random number generator to shoudl give the same numbers as
% above
rng(s);
r5 = normrnd(0,1,[1,5]);
























