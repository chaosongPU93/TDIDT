% function gen_loc_errorest_gauss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to generate a random gaussian distribution of the
% +/-2 samples in the detection, in order to see how does that translate
% into space (shape). Previously we use 'gen_loc_resolution_time' to 
% generate a fixed +/-2-sample in off12 and 13 to see how does the sample diff
% translate to space, but that give us a pseudo error ellipse. For this code,
% we like to generate a dense grid to see more clearly the variation of this 
% error in space with the azimuth.
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/08/24
% Last modified date:   2021/08/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clc
clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);

% set path
workpath = getenv('MHOME');
lzbpath = strcat(workpath, '/Seisbasics/hypoinverse/lzbrst');


%% LZB trio
% fams used
nfampool = [
%             '002';
            '043';
%             '141';
%             '047';
%             '010';
%             '144';
%             '099';
%             '068';
%             '125';
%             '147';
%             '017';
%             '006';
%             '001';
            ];
nfam = size(nfampool, 1);

% relative arrival time to main station for each fam
contlzboff = [
%               14 -22;
              -14 -1;
%               -23 7;
%               2 -15;
%               -25 -3;
%               -26 10;
%               -26 6;
%               -8 -1;
%               -4 -6;
%               -15 4;
%               -18 6;
%               -32 9;
%               -32 5;
              ];

% +2/-2 samples in 40 sps, equivalent to +1/-1 samples in 20 sps, ie, LF resolution
rng('default');

DIST = 'gaussian';
% DIST = 'uniform';
%%%%%% Gaussian distribution
if isequal(DIST, 'gaussian')
    mu = 0;
    sigma = 1;
    off12raw = mu+ sigma*randn(5000,1);    % mean of 0, std of 1;
    off13raw = mu+ sigma*randn(5000,1);    % mean of 0, std of 1;
    j = 0;
%     for i = 1: length(off12raw)
%         off23raw = off13raw(i)-off12raw(i);     % constraint
%         if abs(off13raw(i))<=2*sigma && abs(off12raw(i))<=2*sigma && abs(off23raw)<=2*sigma
%             j = j+1;
%             off12(j) = off12raw(i);
%             off13(j) = off13raw(i);
%         end
%     end
    for i = 1: length(off12raw)
        if abs(off13raw(i))<=2*sigma && abs(off12raw(i))<=2*sigma
            j = j+1;
            off12(j) = off12raw(i);
            off13(j) = off13raw(i);
        end
    end
    
    
%%%%%% Uniform distribution
elseif isequal(DIST, 'uniform')
    mu = 0;
    fact = 4;
    off12raw = mu+ fact*(rand(5000,1)-0.5);    % mean of 0, std of [-2,2];
    off13raw = mu+ fact*(rand(5000,1)-0.5);    % mean of 0, std of [-2,2];
    j = 0;
    for i = 1: length(off12raw)
        off23raw = off13raw(i)-off12raw(i);     % constraint
        if abs(off13raw(i))<=fact/2 && abs(off12raw(i))<=fact/2 && abs(off23raw)<=fact/2
            j = j+1;
            off12(j) = off12raw(i);
            off13(j) = off13raw(i);
        end
    end
end

%%% fc is inverted with (lf-hf_12, lf-hf_13) relative to all fams
%%% the same content as /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/LZBallfamfc, which is
%%% inverted by hypoinverse from /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/chat00p.fc
fclzb = [
%     -123.509500 48.457000 38.0200;
    -123.801000 48.561333 36.2500;
%     -123.896333 48.594000 35.7800;
%     -123.638833 48.474000 36.7200;
%     -123.797167 48.440333 34.8100;
%     -123.925000 48.599333 35.5600;
%     -123.898667 48.555833 35.2500;
%     -123.772167 48.575333 36.7300;
%     -123.734833 48.562667 36.8900;
%     -123.837500 48.587500 36.3200;
%     -123.867833 48.590000 36.0100;
%     -123.930667 48.545167 34.8600;       % 006
%     -123.892500 48.467333 34.3600;       % 001
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VERSION 1 for the locations of lfe families, by inversion from hypoinverse
%%% this is inverted from (0,0) of all fams, same order, location of control points
%%% inverted from /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/chat00p.reffam
contlzb = [
%     -123.492667 48.451500 38.1400;
    -123.772167 48.493000 35.5900;
%     -123.863167 48.528167 35.2100;
%     -123.603333 48.440167 36.7100;
%     -123.800167 48.408833 34.5200;
%     -123.893333 48.536500 35.0700;
%     -123.864500 48.498667 34.8800;
%     -123.753333 48.525667 36.2000;
%     -123.703667 48.502667 36.4100;
%     -123.814333 48.538667 35.7900;
%     -123.838500 48.544833 35.6600;
%     -123.908000 48.494167 34.5100;       % 006
%     -123.879667 48.446167 34.2600;       % 001
    ];

arr = [];
for ifam = 1: nfam
    for i = 1: length(off12)
        tmp(i,1:2) = contlzboff(ifam,:);
        tmp(i,3:4) = [off12(i) off13(i)];
    end
    arr = [arr; tmp];        
end

figure
ax = gca;
scatter(ax,off12,off13,10,'k','filled');
axis(ax,[-2.5 2.5 -2.5 2.5])
axis equal

if isequal(DIST, 'gaussian')
    fid = fopen(strcat(lzbpath, '/offset_043_locerr2_40_gaus_sq'),'w+');
    fprintf(fid,'%d %d %f %f \n',arr');
    fclose(fid);
elseif isequal(DIST, 'uniform')
    fid = fopen(strcat(lzbpath, '/offset_043_locerr2_40_unif'),'w+');
    fprintf(fid,'%d %d %f %f \n',arr');
    fclose(fid);
end








