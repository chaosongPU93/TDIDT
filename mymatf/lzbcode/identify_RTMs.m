% function identify_RTMs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to find the best propagation direction for any possible
% periods of tremor bursts from 'identify_tremor_bursts_visually.m'
% Search over a series of trial directions, project HF directions onto this
% direciton, then apply a robust linear regression to projections vs. time,
% obtain the standard error of the slope and the Pearson correlation coefficient
% Try to find a threshold to discard some bad ones. The rest are viewed as
% RTMs that have a unidirectional propagation direction
%
%   
% NOTES:
%   2020/09/08, i am adding a new fam 006, so that i want to check if original
%               migrations are still reasonable
%   2020/09/10, i am adding a fam 001 back again for last try, 13 fam in total now
% 
% Chao Song, chaosong@princeton.edu
% First created date:   2020/12/24
% Last modified date:   2020/12/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
% format short e   % Set the format to 5-digit floating point
clear
close all
clc

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

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
            '017';
            '006';
            '001';
            ];     % 2020/09/07, i am adding a new family 006 to the pool, final version

nfam = size(nfampool,1);
disp(nfam); 

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

  
% load files
winlenhf = 4;

winlenlf = 16;
lofflf = 4;
ccminlf = 0.35;

SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.dcutnodou.',SUFFIXhf);
hftime = load(fname);
% 25 cols, format is:
%   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date
%   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum
%


SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.dcutnodou.',SUFFIXlf);
lftime = load(fname);

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
           -123.838500 48.544833 35.6600;
           -123.908000 48.494167 34.5100;       % 006
           -123.879667 48.446167 34.2600;       % 001
           ];


relacont = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),loccont(2,1),loccont(2,2));
relacont(:,1) = dx;
relacont(:,2) = dy;


% sort according to day, sec
hftime = sortrows(hftime, [13, 15]);
lftime = sortrows(lftime, [13, 15]);


%%
% periods totally obtained visually, i.e., time-space clustered
% commented times are periods within which the number of detections in hf/lf is lower than thres.
ntran04vis = [
%    2004195  2.8512e+04   2.9600e+04;
   2004196  3.2400e+04   3.7800e+04;
   2004196  4.5220e+04   4.8800e+04;
   2004196  7.2000e+04   7.6032e+04;
   2004197  7.7760e+03   8.6400e+03;
   2004197  1.1880e+04   1.2960e+04;
   2004197  1.8540e+04   1.9872e+04;
   2004197  2.4192e+04   2.6784e+04;
   2004197  3.0060e+04   3.2832e+04;
   2004197  3.3910e+04   3.5424e+04;
   2004197  3.6000e+04   3.6900e+04;
   2004197  4.0608e+04   4.4064e+04;
   2004197  4.4496e+04   4.5792e+04;
   2004197  5.1120e+04   5.2704e+04;
   2004197  5.4432e+04   5.7024e+04;
   2004197  6.9984e+04   7.3440e+04;
   2004197  7.8480e+04   7.9488e+04;
   2004197  7.9920e+04   8.2080e+04;
   2004198  9.9000e+03   1.0980e+04;
   2004198  2.1780e+04   2.3688e+04;
   2004198  3.1968e+04   3.4560e+04;
   2004198  4.3200e+04   4.3848e+04;
   2004198  6.2640e+04   6.4440e+04;
   2004199  2.5920e+03   4.3200e+03;
   2004199  6.6240e+03   7.5600e+03;
   2004199  1.7640e+04   1.8720e+04;
   2004199  3.2400e+04   3.3840e+04;
   2004199  3.5424e+04   3.7152e+04;
   2004199  4.2840e+04   4.3380e+04;
   2004199  5.7024e+04   6.1344e+04;
%    2004199  6.6744e+04   6.7320e+04;
   2004199  8.2080e+04   8.3520e+04;
   2004200  1.0800e+04   1.2096e+04;
   2004200  2.8260e+04   2.9376e+04;
   2004200  4.1472e+04   4.5792e+04;
   2004200  4.7520e+04   5.2704e+04;
   2004201           0   1.6200e+03;
   2004201  4.3200e+03   9.5040e+03;
   2004201  2.4192e+04   2.6784e+04;
%    2004202  1.4688e+04   1.6560e+04;
   2004202  4.6440e+04   4.8384e+04;
%    2004202  5.0544e+04   5.2200e+04;
   2004203  5.9400e+04   6.3072e+04;
   2004203  6.5880e+04   6.7320e+04;
   ];

% periods totally obtained visually, i.e., time-space clustered
% commented times are periods within which the number of detections in hf/lf is lower than thres.
ntran05vis = [
%    2005254  1.9872e+04   2.1600e+04;
   2005254  5.1840e+04   5.4432e+04;
   2005254  6.7392e+04   7.0560e+04;
%    2005255  4.3200e+03   6.0480e+03;
   2005255  1.8180e+04   1.8900e+04;
%    2005255  2.2680e+04   2.3328e+04;
%    2005255  4.3200e+04   4.5792e+04;
   2005255  5.0976e+04   5.1950e+04;
   2005255  5.8320e+04   5.9400e+04;
   2005255  5.9616e+04   6.1920e+04;
   2005255  7.1712e+04   7.2360e+04;
   2005256  1.1880e+04   1.6416e+04;
   2005256  1.8144e+04   1.9872e+04;
   2005256  1.9080e+04   3.1968e+04;
   2005256  3.2400e+04   3.3408e+04;
   2005256  3.5424e+04   3.7152e+04;
   2005256  4.5360e+04   4.6656e+04;
%    2005256  6.9552e+04   7.0416e+04;
   2005256  7.7760e+04   7.9920e+04;
   2005256  8.0784e+04   8.4672e+04;
   2005257  2.8800e+03   7.5600e+03;
   2005257  7.5600e+03   1.7280e+04;
%    2005257  2.5488e+04   2.6352e+04;
   2005258  2.4624e+04   2.6280e+04;
   2005260  1.2600e+03   2.5920e+03;
%    2005260  8.6400e+03   9.5040e+03;
   2005260  4.1040e+04   4.3200e+04;
   2005261  9.9360e+03   1.2960e+04;
   ];

% bursts resulting from automated methods in 'identify_tremor_bursts_intertime.m' 
% with ttol=5e-4
ntran04 = [
    2004196   7.8820e+04   7.9181e+04;
    2004196   8.1252e+04   8.1821e+04;
    2004196   8.5875e+04   8.6292e+04;
    2004197   8.3030e+03   8.6740e+03;
    2004197   1.8285e+04   1.9175e+04;
    2004197   2.6409e+04   2.7106e+04;
    2004197   2.9788e+04   3.0480e+04;
    2004197   3.0986e+04   3.1493e+04;
    2004197   3.6529e+04   3.7063e+04;
    2004197   4.1348e+04   4.1759e+04;
    2004197   4.1751e+04   4.2198e+04;
    2004197   4.2453e+04   4.2863e+04;
    2004197   4.2971e+04   4.3692e+04;
    2004197   4.4603e+04   4.5221e+04;
    2004197   4.6255e+04   4.6846e+04;
    2004197   5.4674e+04   5.5597e+04;
    2004197   5.6210e+04   5.6754e+04;
    2004197   7.0349e+04   7.1040e+04;
    2004197   7.7824e+04   7.9047e+04;
    2004197   8.4602e+04   8.5666e+04;
    2004197   8.5620e+04   8.6255e+04;
    2004198   5.8930e+03   8.3360e+03;
    2004198   8.2850e+03   9.3250e+03;
    2004198   2.0175e+04   2.1221e+04;
    2004198   2.6696e+04   2.7321e+04;
    2004198   3.2740e+04   3.4163e+04;
    2004198   4.2572e+04   4.3612e+04;
    2004198   5.3959e+04   5.5213e+04;
    2004198   5.5283e+04   5.6047e+04;
    2004198   5.5978e+04   5.7886e+04;
    2004198   5.9650e+04   6.0424e+04;
    2004198   6.3137e+04   6.3502e+04;
    2004198   7.3572e+04   7.5560e+04;
    2004198   8.4757e+04   8.5300e+04;
    2004198   8.5405e+04   8.5982e+04;
    2004198   8.6138e+04   8.6984e+04;
    2004199   5.2040e+03   6.0940e+03;
    2004199   3.3120e+04   3.3751e+04;
    2004199   3.6033e+04   3.6406e+04;
    2004199   4.2785e+04   4.3238e+04;
    2004199   4.5726e+04   4.6072e+04;
    2004199   4.6280e+04   4.7096e+04;
    2004199   4.7167e+04   4.7693e+04;
    2004199   4.7842e+04   4.8471e+04;
    2004199   4.8424e+04   4.9087e+04;
    2004199   4.9684e+04   5.0674e+04;
    2004199   6.6532e+04   6.7428e+04;
    2004199   8.0920e+04   8.2747e+04;
    2004199   8.2803e+04   8.3492e+04;
    2004200   1.2202e+04   1.2852e+04;
    2004200   1.3039e+04   1.3679e+04;
    2004200   1.3840e+04   1.4921e+04;
    2004200   1.4847e+04   1.5390e+04;
    2004200   1.6274e+04   1.7020e+04;
    2004200   1.9045e+04   1.9938e+04;
    2004200   4.8049e+04   4.8636e+04;
    2004201   6.9880e+03   7.4390e+03;
    2004203   1.6544e+04   1.7367e+04;
    2004203   1.8567e+04   2.0621e+04;
    2004203   6.6477e+04   6.6932e+04;
   ];

% bursts resulting from automated methods in 'identify_tremor_bursts_intertime.m' 
% with ttol=5e-4
ntran05 = [
    2005255   3.4544e+04   3.5264e+04;
    2005255   4.6043e+04   4.6664e+04;
    2005255   5.1221e+04   5.2005e+04;
    2005255   5.8831e+04   5.9150e+04;
    2005255   6.7604e+04   6.8595e+04;
    2005255   7.1573e+04   7.2079e+04;
    2005255   7.3802e+04   7.4175e+04;
    2005255   7.5190e+04   7.5652e+04;
    2005255   8.4452e+04   8.4931e+04;
    2005255   8.5070e+04   8.5607e+04;
    2005256   3.5600e+03   4.7570e+03;
    2005256   8.4960e+03   9.1590e+03;
    2005256   1.3662e+04   1.4406e+04;
    2005256   1.4405e+04   1.4918e+04;
    2005256   1.9369e+04   2.0099e+04;
    2005256   2.1448e+04   2.2185e+04;
    2005256   2.8275e+04   2.9232e+04;
    2005256   2.9655e+04   3.0166e+04;
    2005256   3.0118e+04   3.0856e+04;
    2005256   3.2691e+04   3.3392e+04;
    2005256   3.4756e+04   3.5298e+04;
    2005256   3.7056e+04   3.7804e+04;
    2005256   3.7741e+04   3.8384e+04;
    2005256   3.8704e+04   3.9587e+04;
    2005256   4.1848e+04   4.3176e+04;
    2005256   4.3116e+04   4.4808e+04;
    2005256   4.5256e+04   4.6138e+04;
    2005256   4.6147e+04   4.7961e+04;
    2005256   4.8992e+04   4.9804e+04;
    2005256   5.1664e+04   5.2925e+04;
    2005256   5.5560e+04   5.6630e+04;
    2005256   6.0663e+04   6.1414e+04;
    2005256   6.1734e+04   6.2684e+04;
    2005256   6.2753e+04   6.3398e+04;
    2005256   7.1457e+04   7.3974e+04;
    2005256   7.6321e+04   7.6954e+04;
    2005256   7.7265e+04   7.7905e+04;
    2005256   8.4849e+04   8.5324e+04;
    2005257   5.9820e+03   6.5190e+03;
    2005257   6.4500e+03   7.6100e+03;
    2005257   1.0305e+04   1.1498e+04;
    2005257   1.2486e+04   1.3826e+04;
    2005257   1.5913e+04   1.7852e+04;
    2005257   2.1010e+04   2.2196e+04;
    2005257   2.2124e+04   2.4116e+04;
    2005257   2.4045e+04   2.4636e+04;
    2005257   2.5725e+04   2.7224e+04;
    2005257   3.5290e+04   3.7587e+04;
    2005257   3.8537e+04   4.3992e+04;
    2005257   4.3927e+04   4.4577e+04;
    2005257   4.7110e+04   4.8678e+04;
    2005257   6.1772e+04   6.2608e+04;
    2005257   6.9150e+04   6.9952e+04;
    2005257   7.3657e+04   7.7914e+04;
    2005257   7.8604e+04   8.0071e+04;
    2005257   8.0004e+04   8.1009e+04;
    2005257   8.1108e+04   8.2051e+04;
    2005257   8.2151e+04   8.2872e+04;
    2005258   2.4757e+04   2.6277e+04;
    2005258   3.4940e+04   3.6505e+04;
    2005258   3.6756e+04   3.7190e+04;
    2005258   3.7158e+04   3.8247e+04;
    2005258   3.8178e+04   3.9933e+04;
    2005259   9.4000e+01   1.6940e+03;
    2005259   2.1990e+03   3.8540e+03;
    2005259   5.2490e+03   6.3820e+03;
    2005259   6.5320e+03   7.4710e+03;
    2005259   3.7921e+04   3.9524e+04;
    2005259   3.9968e+04   4.0433e+04;
    2005259   4.0366e+04   4.1279e+04;
    2005259   7.4469e+04   7.4783e+04;
    2005260   2.6830e+03   3.1550e+03;
    2005260   3.5970e+03   4.5700e+03;
    2005260   4.5030e+03   6.6730e+03;
    2005260   6.7930e+03   7.0910e+03;
    2005260   9.5230e+03   1.0210e+04;
    2005260   5.6362e+04   5.6817e+04;
    2005260   5.6894e+04   5.7240e+04;
    2005261   1.1290e+04   1.2163e+04;
   ];
   
% bursts resulting from automated methods in 'identify_tremor_bursts_intertime.m' 
% with ttol=1e-3;
ntran03 = [
    2003064    1.7589e+04   1.9399e+04 
    ];

ntran04 = [
    2004195   4.1784e+04   4.2666e+04
    2004195   4.2745e+04   4.3392e+04
    2004196   3.5325e+04   3.7173e+04
    2004196   4.6700e+04   4.9040e+04
    2004196   6.9421e+04   7.1726e+04
    2004196   7.6128e+04   7.6577e+04
    2004196   7.8850e+04   7.9202e+04
    2004196   8.1282e+04   8.1791e+04
    2004196   8.5905e+04   8.6262e+04
    2004197   5.3810e+03   5.7990e+03
    2004197   6.7010e+03   7.2020e+03
    2004197   8.3330e+03   9.3110e+03
    2004197   1.3635e+04   1.4281e+04
    2004197   1.7984e+04   2.0411e+04
    2004197   2.3906e+04   2.5930e+04
    2004197   2.6439e+04   2.7195e+04
    2004197   2.9084e+04   3.1862e+04
    2004197   3.3831e+04   3.5500e+04
    2004197   3.6559e+04   3.7198e+04
    2004197   3.7241e+04   3.7761e+04
    2004197   4.0944e+04   4.2168e+04
    2004197   4.2292e+04   4.2971e+04
    2004197   4.3001e+04   4.3662e+04
    2004197   4.4633e+04   4.5752e+04
    2004197   4.6227e+04   4.6816e+04
    2004197   5.1066e+04   5.2277e+04
    2004197   5.4704e+04   5.5897e+04
    2004197   5.6240e+04   5.7257e+04
    2004197   6.9920e+04   7.1142e+04
    2004197   7.7484e+04   7.9841e+04
    2004197   7.9919e+04   8.1483e+04
    2004197   8.1517e+04   8.2240e+04
    2004197   8.2282e+04   8.3073e+04
    2004197   8.3766e+04   8.6638e+04
    2004198   5.9230e+03   9.7670e+03
    2004198   9.8900e+03   1.0628e+04
    2004198   1.0790e+04   1.1471e+04
    2004198   1.9998e+04   2.1646e+04
    2004198   2.6610e+04   2.7291e+04
    2004198   3.2495e+04   3.5376e+04
    2004198   4.2602e+04   4.3751e+04
    2004198   5.3476e+04   5.8080e+04
    2004198   5.9624e+04   6.0467e+04
    2004198   6.2936e+04   6.3915e+04
    2004198   7.3510e+04   7.6059e+04
    2004198   8.4192e+04   8.5952e+04
    2004198   8.6118e+04   8.7278e+04
    2004199   1.2710e+03   2.0930e+03
    2004199   2.1580e+03   3.1600e+03
    2004199   3.4790e+03   3.8300e+03
    2004199   4.1370e+03   6.1680e+03
    2004199   6.3370e+03   6.8880e+03
    2004199   3.2419e+04   3.3857e+04
    2004199   3.5955e+04   3.6376e+04
    2004199   4.1691e+04   4.2144e+04
    2004199   4.2815e+04   4.3333e+04
    2004199   4.5756e+04   4.6042e+04
    2004199   4.6310e+04   4.7066e+04
    2004199   4.7197e+04   4.7663e+04
    2004199   4.7777e+04   4.9156e+04
    2004199   4.9714e+04   5.0644e+04
    2004199   5.9481e+04   6.0765e+04
    2004199   6.6239e+04   6.7398e+04
    2004199   8.0831e+04   8.3462e+04
    2004200   1.1676e+04   1.5829e+04
    2004200   1.5884e+04   1.7155e+04
    2004200   1.9075e+04   2.0140e+04
    2004200   2.8286e+04   2.9498e+04
    2004200   4.8079e+04   4.8696e+04
    2004200   4.9633e+04   5.0437e+04
    2004201   8.3400e+02   1.6590e+03
    2004201   6.8970e+03   7.7990e+03
    2004203   1.6574e+04   1.7649e+04
    2004203   1.8492e+04   2.0806e+04
    2004203   6.1190e+04   6.1630e+04
    2004203   6.4548e+04   6.5397e+04
    2004203   6.6308e+04   6.6902e+04
    ];

% bursts resulting from automated methods in 'identify_tremor_bursts_intertime.m' 
% with ttol=1e-3;
ntran05 = [
    2005254   5.6351e+04   5.7024e+04
    2005254   7.0830e+04   7.1642e+04
    2005255   3.4370e+04   3.5642e+04
    2005255   4.6069e+04   4.6979e+04
    2005255   5.1099e+04   5.1975e+04
    2005255   5.3299e+04   5.3837e+04
    2005255   5.8020e+04   5.9341e+04
    2005255   6.0351e+04   6.1587e+04
    2005255   6.1648e+04   6.2517e+04
    2005255   6.7238e+04   6.8565e+04
    2005255   7.1603e+04   7.2049e+04
    2005255   7.3832e+04   7.4266e+04
    2005255   7.4760e+04   7.5622e+04
    2005255   8.4266e+04   8.4901e+04
    2005255   8.5100e+04   8.5664e+04
    2005256   3.5900e+03   4.9870e+03
    2005256   8.0030e+03   1.0384e+04
    2005256   1.2890e+04   1.3299e+04
    2005256   1.3382e+04   1.4376e+04
    2005256   1.4435e+04   1.4992e+04
    2005256   1.5396e+04   1.6040e+04
    2005256   1.9039e+04   2.0233e+04
    2005256   2.0407e+04   2.2324e+04
    2005256   2.7689e+04   2.9202e+04
    2005256   2.9685e+04   3.1233e+04
    2005256   3.2606e+04   3.4264e+04
    2005256   3.4786e+04   3.5268e+04
    2005256   3.5597e+04   4.1142e+04
    2005256   4.1174e+04   4.5219e+04
    2005256   4.5286e+04   4.8260e+04
    2005256   4.8818e+04   5.0215e+04
    2005256   5.1364e+04   5.4083e+04
    2005256   5.5237e+04   5.7052e+04
    2005256   6.0693e+04   6.4846e+04
    2005256   7.1030e+04   7.3944e+04
    2005256   7.6351e+04   7.6924e+04
    2005256   7.6984e+04   7.8694e+04
    2005256   7.8931e+04   8.0564e+04
    2005256   8.1726e+04   8.2092e+04
    2005256   8.4791e+04   8.5294e+04
    2005257   4.9140e+03   7.7760e+03
    2005257   8.9780e+03   1.3883e+04
    2005257   1.5943e+04   1.8685e+04
    2005257   2.1040e+04   2.5395e+04
    2005257   2.5465e+04   2.8085e+04
    2005257   3.5320e+04   3.7628e+04
    2005257   3.8483e+04   4.5322e+04
    2005257   4.6123e+04   4.8648e+04
    2005257   6.1419e+04   6.3758e+04
    2005257   6.9110e+04   6.9922e+04
    2005257   7.3563e+04   8.2073e+04
    2005257   8.2129e+04   8.3410e+04
    2005257   8.3484e+04   8.4261e+04
    2005258   2.4787e+04   2.6247e+04
    2005258   2.6329e+04   2.8992e+04
    2005258   3.4970e+04   3.9993e+04
    2005259   1.2400e+02   1.6640e+03
    2005259   1.9280e+03   4.0760e+03
    2005259   5.2020e+03   7.7440e+03
    2005259   3.7267e+04   4.1249e+04
    2005259   7.4225e+04   7.4753e+04
    2005260   1.8820e+03   2.5110e+03
    2005260   2.6210e+03   3.4680e+03
    2005260   3.5740e+03   6.7080e+03
    2005260   6.7390e+03   7.0610e+03
    2005260   7.0970e+03   7.5000e+03
    2005260   7.6440e+03   8.0480e+03
    2005260   9.3790e+03   1.0387e+04
    2005260   5.1863e+04   5.2397e+04
    2005260   5.6392e+04   5.7210e+04
    2005261   8.5110e+03   9.6510e+03
    2005261   9.7170e+03   1.0343e+04
    2005261   1.0653e+04   1.1148e+04
    2005261   1.1320e+04   1.2133e+04
   ];

% visual modifications to original time range, and need some further check
ntrantest = [
    2004197   8.3766e+04   8.4888e+04   % old 6, new 34, ACCEPTED
    2004197   8.4888e+04   8.5248e+04   % old 7, new 34, too large SE, DISCARDED
    2004197   8.5320e+04   8.6638e+04   % old 8, new 34, ACCEPTED
    2005255	  34370        35642        % new 3, old 21, too small Pearson, DISCARDED
    2005255   6.7238e+04   6.8565e+04   % new 10, old 22, questionable direction, 155-->65
    2005256   3590         4987         % new 16, old 25, questionable direction, 175-->245
    2005256   8.0030e+03   9.8000e+03   % new 17, old 26, change ending time as old, ACCEPTED
    2005256   2.8260e+04   2.9202e+04   % new 24, old 28, change start time as 7.85, ACCEPTED
    2005256   5.2020e+04   5.3208e+04   % new 32, old 30, try a new range, ACCEPTED
    2005256   5.5548e+04   5.7052e+04   % new 33, change start time as 15.43, ACCEPTED
    2005256   7.1496e+04   7.3944e+04   % new 35, old 31, change start time as 19.86, ACCEPTED
    2005257   6120         7776         % new 41, old 33, change start time as 1.7, ACCEPTED 
];

%%
trange = ntran05;

angbestl1 = zeros(size(trange,1),1);
angbestl2 = zeros(size(trange,1),1);
angbestl3 = zeros(size(trange,1),1);
angbestl4 = zeros(size(trange,1),1);
angbestl5 = zeros(size(trange,1),1);

xran = [-15 25];
yran = [-20 20];

resprophf = nan(size(trange,1)+1,200);
resproplf = nan(size(trange,1)+1,50);

% create fit object with free constraints
fttpfree = fittype( @(a,b,x) a*x+b);

indcheck = 1: size(trange,1);
indplt = indcheck;

for i = 1: length(indplt)
% for i = 41: 63
    disp(trange(indplt(i),:));
    indhf = find(hftime(:,13)==trange(indplt(i),1) & hftime(:,15)>=trange(indplt(i),2) & ...
                 hftime(:,15)<=trange(indplt(i),3));
    mighf = hftime(indhf,:);

    indlf = find(lftime(:,13)==trange(indplt(i),1) & lftime(:,15)>=trange(indplt(i),2) & ...
                 lftime(:,15)<=trange(indplt(i),3));
    miglf = lftime(indlf,:);
    
    if size(mighf,1)< 15 || size(miglf,1) < 15
        disp(strcat(num2str(indplt(i)),' not enough detections'));
    end
% end
    angle = 0:5:360;
    
%     l1normhf = zeros(length(angle),1);
%     l2normhf = zeros(length(angle),1);
    slopehf = zeros(length(angle),1);
%     ssehf = zeros(length(angle),1);
    rmsehf = zeros(length(angle),1);
    rsquarehf = zeros(length(angle),1);
    
%     l1normlf = zeros(length(angle),1);
%     l2normlf = zeros(length(angle),1);
    slopelf = zeros(length(angle),1);
%     sself = zeros(length(angle),1);
    rmself = zeros(length(angle),1);
    rsquarelf = zeros(length(angle),1);
    
    for iang = 1: length(angle)
        %%% propagation trial of hf
        mighfdum = mighf;
        for j = 1: size(mighf,1)
            x0 = mighfdum(j,1);
            y0 = mighfdum(j,2);
            [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),0,0);
            mighfdum(j,1) = newx;
            mighfdum(j,2) = newy;
        end        
        % linear robust least square
        [fitobj,gof,~] = fit(mighfdum(:,15)/3600, mighfdum(:,1),fttpfree,'Robust','Bisquare',...
                                  'StartPoint',[1 1]);
        % output fit parameters
        coef = coeffvalues(fitobj);
        slopehf(iang) = coef(1);
        fitprophf = feval(fitobj,mighfdum(:,15)/3600);
%         l1normhf(iang) = sum(abs(mighfdum(:,1)-fitprophf))/(length(mighfdum(:,1)));
%         l2normhf(iang) = sum((mighfdum(:,1)-fitprophf).^2)/(length(mighfdum(:,1)));
%         ssehf(iang) = gof.sse;
        rmsehf(iang) = gof.rmse;
        rsquarehf(iang) = gof.rsquare;    % r-square, defined as ratio of the sum of squares of the regression and the total sum of squares
    
        %%% propagation trial of lf
        miglfdum = miglf;
        for j = 1: size(miglf,1)
            x0 = miglfdum(j,1);
            y0 = miglfdum(j,2);
            [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),0,0);
            miglfdum(j,1) = newx;
            miglfdum(j,2) = newy;
        end
        % linear robust least square
        [fitobj,gof,~] = fit(miglfdum(:,15)/3600, miglfdum(:,1),fttpfree,'Robust','Bisquare',...
                                  'StartPoint',[0.5 5]);
        % output fit parameters
        coef = coeffvalues(fitobj);
        slopelf(iang) = coef(1);
        fitproplf = feval(fitobj,miglfdum(:,15)/3600);
%         l1normlf(iang) = sum(abs(miglfdum(:,1)-fitproplf))/(length(miglfdum(:,1)));
%         l2normlf(iang) = sum((miglfdum(:,1)-fitproplf).^2)/(length(miglfdum(:,1)));
%         sself(iang) = gof.sse;
        rmself(iang) = gof.rmse;
        rsquarelf(iang) = gof.rsquare;    % r-square, defined as ratio of the sum of squares of the regression and the total sum of squares
        
    end

    %%% best angle estimate from hf
    ind = find(slopehf>0);
%     ind1 = find(l1normhf(ind)==min(l1normhf(ind)));
%     angbestl1(indplt(i)) = angle(ind(ind1(1)));
%     ind2 = find(l2normhf(ind)==min(l2normhf(ind)));
%     angbestl2(indplt(i)) = angle(ind(ind2(1)));
    ind3 = find(rmsehf(ind)==min(rmsehf(ind)));     % rmse, Root Mean Squared Error, the smaller, the better
    angbestl3(indplt(i)) = angle(ind(ind3(1)));
%     ind4 = find(rsquarehf(ind)==max(rsquarehf(ind)));   % R-square, 0-1, the larger the better
%     angbestl4(indplt(i)) = angle(ind(ind4(1)));
%     ind5 = find(ssehf(ind)==min(ssehf(ind)));   % sum of square due to error, the smaller, the better
%     angbestl5(indplt(i)) = angle(ind(ind5(1)));
%     disp([angbestl3(indplt(i)),angbestl4(indplt(i)),angbestl5(indplt(i))])
%     angbesthf(indplt(i)) = median([angbestl3(indplt(i)),angbestl4(indplt(i)),angbestl5(indplt(i))]);
    angbesthf(indplt(i)) = angbestl3(indplt(i));
    
    
    %%% best angle estimate from hf
    ind = find(slopelf>0);
%     ind1 = find(l1normlf(ind)==min(l1normlf(ind)));
%     angbestl1(indplt(i)) = angle(ind(ind1(1)));
%     ind2 = find(l2normlf(ind)==min(l2normlf(ind)));
%     angbestl2(indplt(i)) = angle(ind(ind2(1)));
    ind3 = find(rmself(ind)==min(rmself(ind)));     % rmse, Root Mean Squared Error, the smaller, the better
    angbestl3(indplt(i)) = angle(ind(ind3(1)));
%     ind4 = find(rsquarelf(ind)==max(rsquarelf(ind)));   % R-square, 0-1, the larger the better
%     angbestl4(indplt(i)) = angle(ind(ind4(1)));
%     ind5 = find(sself(ind)==min(sself(ind)));   % sum of square due to error, the smaller, the better
%     angbestl5(indplt(i)) = angle(ind(ind5(1)));
%     disp([angbestl3(indplt(i)),angbestl4(indplt(i)),angbestl5(indplt(i))])
%     angbestlf(indplt(i)) = median([angbestl3(indplt(i)),angbestl4(indplt(i)),angbestl5(indplt(i))]);
    angbestlf(indplt(i)) = angbestl3(indplt(i));
    
    % determine if more inspection is needed
    if abs(angbesthf(indplt(i))-angbestlf(indplt(i))) > 20
        disp('Difference in propagation direction between hf & lf is too large!');
    end
    
end

angbesthf = angbesthf';
angbestlf = angbestlf';

% angbest = [angbesthf angbestlf];
angbest = angbesthf;

%%
for i = 1: size(trange,1)
% for i = 1: 211
%     i=5;
    disp(trange(i,:));
    indhf = find(hftime(:,13)==trange(i,1) & hftime(:,15)>=trange(i,2) & ...
                 hftime(:,15)<=trange(i,3));
    mighf = hftime(indhf,:);
    
    mighfdum = mighf;
    for j = 1: size(mighf,1)
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        mighfdum(j,1) = newx;
        mighfdum(j,2) = newy;
    end
    
    indlf = find(lftime(:,13)==trange(i,1) & lftime(:,15)>=trange(i,2) & ...
                 lftime(:,15)<=trange(i,3));
    miglf = lftime(indlf,:); 
    
    miglfdum = miglf;
    for j = 1: size(miglf,1)
        x0 = miglfdum(j,1);
        y0 = miglfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        miglfdum(j,1) = newx;
        miglfdum(j,2) = newy;
    end
    
    %%% define and position the figure frame and axes of each plot
    f.fig=figure;
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 9;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    nrow = 2;
    ncol = 2;    
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
    end

    %%% reposition
    set(f.ax(1), 'position', [ 0.1, 0.64, 0.32, 0.32]);
    set(f.ax(2), 'position', [ 0.5, 0.64, 0.32, 0.32]);
    set(f.ax(3), 'position', [ 0.1, 0.35, 0.32, 0.22]);
    set(f.ax(4), 'position', [ 0.5, 0.35, 0.32, 0.22]);
%     set(f.ax(5), 'position', [ 0.1, 0.078, 0.32, 0.22]);
%     set(f.ax(6), 'position', [ 0.5, 0.078, 0.32, 0.22]);
    
    % subplot 1 of figure i
    hold(f.ax(1),'on');
    plot(f.ax(1),[-100 100],[0 0],'k--');
    plot(f.ax(1),[0 0],[-100 100],'k--');
    f.ax(1).FontSize = 9;
    mighf = sortrows(mighf,-15);
    scatter(f.ax(1),mighf(:,1),mighf(:,2), 20, mighf(:,15)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(1),'jet');
    c=colorbar(f.ax(1),'SouthOutside');
    pos = f.ax(1).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    juldate = num2str(trange(i,1));
    yr = str2double(juldate(1:4));
    date = str2double(juldate(5:end));
    a = jul2dat(yr,date);
    mo = a(1);
    if mo == 9 
        mo = {' Sep. '};
    elseif mo == 7
        mo = {' Jul. '};
    else
        mo = {' Mar. '};
    end
    day = num2str(a(2));
    yr = num2str(a(3));
    c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
%     c.Label.String = strcat(num2str(trange(i,1)),' of HF',' (hr)');
    c.Label.FontSize = 11;
    caxis(f.ax(1),[trange(i,2)/3600 trange(i,3)/3600])
    text(f.ax(1),0.85,0.1,'HF','FontSize',12,'unit','normalized');
    text(f.ax(1),0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(1),0.04,0.1,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(1).Box = 'on';
    grid(f.ax(1), 'on');
    axis(f.ax(1), 'equal');
    f.ax(1).GridLineStyle = '--';
    f.ax(1).XAxisLocation = 'top';
    medxhf = median(mighf(:,1));
    medyhf = median(mighf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medxhf-rotx medxhf+rotx];
    yvect = [medyhf-roty medyhf+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medxhf-rotx medxhf+rotx];
    yvect = [medyhf-roty medyhf+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
    xticks(f.ax(1),xran(1):5:xran(2));
    yticks(f.ax(1),yran(1):5:yran(2));    
    xlabel(f.ax(1),'E (km)','fontsize',11);
    ylabel(f.ax(1),'N (km)','fontsize',11);
    text(f.ax(1),0.5,0.1,num2str(i),'FontSize',12,'unit','normalized');
    text(f.ax(1),0.83,0.95,strcat(num2str(size(mighf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(1),0.83,0.90,strcat({'in '},num2str(trange(i,3)-trange(i,2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(mighf,1)/(trange(i,3)-trange(i,2)));
    text(f.ax(1),0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    hold(f.ax(1),'off');

    % subplot 2 of figure i
    hold(f.ax(2),'on');
    plot(f.ax(2),[-100 100],[0 0],'k--');
    plot(f.ax(2),[0 0],[-100 100],'k--');
    f.ax(2).FontSize = 9;
    miglf = sortrows(miglf,-15);
    scatter(f.ax(2),miglf(:,1),miglf(:,2), 20, miglf(:,15)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(2),'jet');
    c=colorbar(f.ax(2),'SouthOutside');
    pos = f.ax(2).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
%     c.Label.String = strcat(num2str(trange(i,1)),' of LF',' (hr)');
    c.Label.FontSize = 11;
    caxis(f.ax(2),[trange(i,2)/3600 trange(i,3)/3600])
    text(f.ax(2),0.85,0.1,'LF','FontSize',12,'unit','normalized');
    text(f.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(2),0.04,0.1,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(2).Box = 'on';
    grid(f.ax(2), 'on');
    axis(f.ax(2), 'equal');
    f.ax(2).GridLineStyle = '--';
    f.ax(2).XAxisLocation = 'top';
    medxlf = median(miglf(:,1));
    medylf = median(miglf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
    xticks(f.ax(2),xran(1):5:xran(2));
    yticks(f.ax(2),yran(1):5:yran(2));    
    xlabel(f.ax(2),'E (km)','fontsize',11);
    ylabel(f.ax(2),'N (km)','fontsize',11);
    text(f.ax(2),0.5,0.1,num2str(i),'FontSize',12,'unit','normalized');
    text(f.ax(2),0.83,0.95,strcat(num2str(size(miglf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(2),0.83,0.9,strcat({'in '},num2str(trange(i,3)-trange(i,2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(miglf,1)/(trange(i,3)-trange(i,2)));
    text(f.ax(2),0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    hold(f.ax(2),'off');
    
    % subplot 3 of figure i
    hold(f.ax(3),'on');
    f.ax(3).FontSize = 9;
%     scatter(f.ax(3),mighfdum(:,15)/3600,mighfdum(:,1),20,[1 0 0],'filled','o',...
%             'MarkerEdgeColor','k'); 
    scatter(f.ax(3),mighfdum(:,15)/3600,mighfdum(:,1),20,[0.6 0.6 0.6],'filled','o',...
            'MarkerEdgeColor','k');    
    scatter(f.ax(3),miglfdum(:,15)/3600,miglfdum(:,1),20,[0.6 1 1],'filled','o',...
            'MarkerEdgeColor','k');
    
%     indpool = [12;13;16;19;21;22;29;32;36;40];
%     indpool = [indpool; 3;5;7;17;20;23;25;26;30];
%     if ismember(i,indpool)
        mintime = max(min(mighfdum(:,15)), min(miglfdum(:,15)));
        maxtime = min(max(mighfdum(:,15)), max(miglfdum(:,15)));
        mighfdum2 = mighfdum(mighfdum(:,15)>=mintime & mighfdum(:,15)<=maxtime, :);
        miglfdum2 = miglfdum(miglfdum(:,15)>=mintime & miglfdum(:,15)<=maxtime, :);
%     else
%         mighfdum2 = mighfdum;
%         miglfdum2 = miglfdum;
%     end

    [fitobjhfprop,gofhfprop,outphfprop] = fit(mighfdum2(:,15)/3600, mighfdum2(:,1),fttpfree,'Robust',...
                                              'Bisquare','StartPoint',[1 1]);
    % output fit parameters
    coefprop = coeffvalues(fitobjhfprop);
    slopeprophf(i) = coefprop(1);
    intcptprophf(i) = coefprop(2);
    fitprophf = feval(fitobjhfprop,mighfdum2(:,15)/3600);
    plot(f.ax(3),mighfdum2(:,15)/3600,fitprophf,'-','linewidth',2,'color',[0.6 0.6 0.6]);
        
%     intcptproplf = ones(size(miglfdum(:,end)))\(miglfdum(:,1)-slopeprophf*miglfdum(:,end)/3600); % direct LS fit
%     fitproplf = slopeprophf*miglfdum(:,end)/3600+intcptproplf;
    a=slopeprophf(i);
    fttpfix = fittype( @(b,x) a*x+b);
    [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum2(:,15)/3600, miglfdum2(:,1),fttpfix,'Robust',...
                                              'Bisquare','StartPoint',1);
%     [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum(:,10)/3600, miglfdum(:,1),fttpfix,'Robust',...
%                                               'Bisquare','Weights',miglfdum(:,11));
    fitproplf = feval(fitobjlfprop,miglfdum2(:,15)/3600);
    intcptproplf(i) = coeffvalues(fitobjlfprop);
    offset(i) = intcptproplf(i)-intcptprophf(i);
    plot(f.ax(3),miglfdum2(:,15)/3600,fitproplf,'-','linewidth',2,'color',[0.6 1 1]);
    f.ax(3).Box = 'on';
    grid(f.ax(3), 'on');
    f.ax(3).GridLineStyle = '--';
    xran1 = [trange(i,2)/3600 trange(i,3)/3600];
    yran1 = [round(min([mighfdum(:,1);miglfdum(:,1)]))-1 ...
             round(max([mighfdum(:,1);miglfdum(:,1)]))+1];
    xlim(f.ax(3),xran1);
    ylim(f.ax(3),yran1);
    text(f.ax(3),0.04,0.91,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    xlabel(f.ax(3),'Time (hr)','fontsize',11);
    ylabel(f.ax(3),'Dist. along prop. (km)','fontsize',11); 
    hold(f.ax(3),'off');
    
    % compute the HF weights in robust linear regression, see NOTES above
    rhf = outphfprop.residuals;   % usual residuals
    x = [mighfdum2(:,15)/3600];
    hatmat = x*inv(x'*x)*x';
    h = zeros(size(hatmat,1),1);    % leverage of least square
    for jj = 1 : size(hatmat,1)
        h(jj) = hatmat(jj,jj);
    end    
    radj = rhf./sqrt(1-h);      % adjusted residuals
    K = 4.685;
    s = mad(rhf,1)/0.6745;
    u = radj/(K*s);
    wthf = zeros(length(u),1);    % rubust weight of next iteration
    for jj = 1 : length(u)
        if abs(u(jj)) < 1
            wthf(jj) = (1-(u(jj))^2)^2;
        else
            wthf(jj) = 0;
        end
    end
    
    % get the standard error of the estimated parameters, may indicate the compare the quality
    % of fitting cross different RTMs, NOTE that rmse, mse are definitely not a good measure
    % the obtained CI is comfirmed to be correct by comparing it with the 'fitobjhfprop'
    slopese = gofhfprop.rmse./sqrt(sum((x-mean(x)).^2));
    est = slopeprophf(i);
    slopeprophfCI(i,:) = confidence_interval_general(est,slopese,length(x)-2,95);
    interceptse = slopese.*sqrt(sum(x.^2)./length(x));
    est = intcptprophf(i);
    intcptprophfCI(i,:) = confidence_interval_general(est,interceptse,length(x)-2,95);
    sehf(i, :) = [slopese interceptse]; 
    
    x = mighfdum2(:,15)/3600;
    y = mighfdum2(:,1);
    x_bar = wt_mean(x,wthf);
    y_bar = wt_mean(y,wthf);
    x_var = sum(wthf.*(x-x_bar).^2) / sum(wthf);
    y_var = sum(wthf.*(y-y_bar).^2) / sum(wthf);
    xy_cov = sum(wthf.*(x-x_bar).*(y-y_bar)) / sum(wthf);
    pearwthf(i) = xy_cov / sqrt(x_var*y_var);
    
    % compute the LF weights in robust linear regression, see NOTES above
    rlf = outplfprop.residuals;   % usual residuals
    x = [miglfdum2(:,15)/3600];
    hatmat = x*inv(x'*x)*x';
    h = zeros(size(hatmat,1),1);    % leverage of least square
    for jj = 1 : size(hatmat,1)
        h(jj) = hatmat(jj,jj);
    end    
    radj = rlf./sqrt(1-h);      % adjusted residuals
    K = 4.685;
    s = mad(rlf,1)/0.6745;
    u = radj/(K*s);
    wtlf = zeros(length(u),1);    % rubust weight of next iteration
    for jj = 1 : length(u)
        if abs(u(jj)) < 1
            wtlf(jj) = (1-(u(jj))^2)^2;
        else
            wtlf(jj) = 0;
        end
    end
    
    % subplot 4 of figure i
    hold(f.ax(4),'on');
    f.ax(4).FontSize = 9;
    numhf(i) = length(mighfdum2(:,1));
    numlf(i) = length(miglfdum2(:,1));
    resprophf(i,1:numhf(i)) = outphfprop.residuals;
    resproplf(i,1:numlf(i)) = outplfprop.residuals+(intcptproplf(i)-intcptprophf(i));
    [xlochf, a, b, pdfvalhf] = weighted_bar_pdf(resprophf(i,1:numhf(i)), wthf, 0.5, 'int');
    hfHdl = bar(f.ax(4),xlochf, pdfvalhf,1,'stacked');
    hfHdl(1).FaceColor = [0.6 0.6 0.6];
    [xloclf, ~, ~, pdfvallf] = weighted_bar_pdf(resproplf(i,1:numlf(i)), wtlf, 0.5, 'int');
    lfHdl = bar(f.ax(4),xloclf, pdfvallf,1,'stacked','facea',0.6);
    lfHdl(1).FaceColor = [0.6 1 1];
    text(f.ax(4),0.04,0.91,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(4),0.4,0.92,sprintf('LF-HF = %.2f km',offset(i)),'FontSize',11,'unit','normalized');
    f.ax(4).Box = 'on';
    grid(f.ax(4), 'on');
    f.ax(4).GridLineStyle = '--';
    ymax = f.ax(4).YLim(2)+0.1;
    ylim(f.ax(4),[0 ymax]);
%     ylim(f.ax(4),[0 1]);
    xlabel(f.ax(4),'Residual in prop. (km)','fontsize',11);
    ylabel(f.ax(4),'PDF estimate','fontsize',11);
    hold(f.ax(4),'off');

%     if pearwthf(i) >= 0.5 && sehf(i,1)<3
%         print(f.fig,'-dpdf',strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2020/Figures/nrtm03_',...
%             num2str(i),'.pdf'));
%     end

end

[tmp,tmpind] = sort(sehf(:,1),'descend');
sortsehf = [tmpind tmp];
[tmp,tmpind] = sort(pearwthf(:),'ascend');
sortpearwthf = [tmpind tmp];

%%
tmp1 = sortpearwthf(sortpearwthf(:,2)>=0.5, 1);
tmp2 = sortsehf(sortsehf(:,2)<3, 1);

if issame(trange,ntran04)
    ind04 = intersect(tmp1,tmp2);
    acctran04 = ntran04(ind04,:);
elseif issame(trange,ntran05)
    ind05 = intersect(tmp1,tmp2);
    acctran05 = ntran05(ind05,:);
else
    ind03 = intersect(tmp1,tmp2);
    acctran03 = ntran03(ind03,:);
end


%% hardcopy of the accepted RTMs passed SE and Pearson thresholds
acctran04 = [
2004195   4.1784e+04   4.2666e+04 % 2004195,4.15e+4,4.60e+4;
2004197   3.3831e+04   3.5500e+04
2004197   4.0944e+04   4.2168e+04
2004197   5.4704e+04   5.5897e+04
2004197   5.6240e+04   5.7257e+04
2004197   6.9920e+04   7.1142e+04
2004197   7.7484e+04   7.9841e+04
2004197   7.9919e+04   8.1483e+04
2004197   8.3766e+04   8.6638e+04
2004198   5.9230e+03   9.7670e+03
2004198   1.9998e+04   2.1646e+04
2004198   4.2602e+04   4.3751e+04
2004198   5.3476e+04   5.8080e+04
2004198   5.9624e+04   6.0467e+04
2004198   6.2936e+04   6.3915e+04
2004198   7.3510e+04   7.6059e+04
2004198   8.4192e+04   8.5952e+04
2004199   4.1370e+03   6.1680e+03
2004199   4.2815e+04   4.3333e+04
2004199   4.7777e+04   4.9156e+04
2004199   4.9714e+04   5.0644e+04
2004199   8.0831e+04   8.3462e+04
2004200   1.1676e+04   1.5829e+04
2004200   1.5884e+04   1.7155e+04
2004200   1.9075e+04   2.0140e+04
2004200   4.8079e+04   4.8696e+04
2004201   8.3400e+02   1.6590e+03
2004203   1.8492e+04   2.0806e+04
];

ind04 = [
     1
    18
    21
    27
    28
    29
    30
    31
    34
    35
    38
    41
    42
    43
    44
    45
    46
    51
    56
    60
    61
    64
    65
    66
    67
    69
    71
    74
    ];

angbest04 = [
   260
   155
   165
   220
   185
   195
   195
   195
   260
   235
   205
   230
    80
    85
    65
   235
   250
   310
   250
   135
   195
   140
   170
    65
   215
   240
   260
   120
];


acctran05 = [
    2005255   4.6069e+04   4.6979e+04
    2005255   5.8020e+04   5.9341e+04
    2005255   6.0351e+04   6.1587e+04
    2005255   6.7238e+04   6.8565e+04
    2005256   8.0030e+03   1.0384e+04
    2005256   1.3382e+04   1.4376e+04
    2005256   2.0407e+04   2.2324e+04
    2005256   2.7689e+04   2.9202e+04
    2005256   2.9685e+04   3.1233e+04
    2005256   3.4786e+04   3.5268e+04
    2005256   5.5237e+04   5.7052e+04
    2005256   6.0693e+04   6.4846e+04
    2005256   7.1030e+04   7.3944e+04
    2005256   7.6351e+04   7.6924e+04
    2005256   7.6984e+04   7.8694e+04
    2005257   1.5943e+04   1.8685e+04
    2005257   2.1040e+04   2.5395e+04
    2005257   2.5465e+04   2.8085e+04
    2005257   3.5320e+04   3.7628e+04
    2005257   3.8483e+04   4.5322e+04
    2005257   6.1419e+04   6.3758e+04
    2005257   7.3563e+04   8.2073e+04
    2005257   8.2129e+04   8.3410e+04
    2005258   2.4787e+04   2.6247e+04
    2005258   3.4970e+04   3.9993e+04
    2005259   1.9280e+03   4.0760e+03
    2005259   5.2020e+03   7.7440e+03
    2005259   3.7267e+04   4.1249e+04
    2005259   7.4225e+04   7.4753e+04
    2005260   2.6210e+03   3.4680e+03
    2005260   3.5740e+03   6.7080e+03
    2005260   9.3790e+03   1.0387e+04
];

ind05 = [
     4
     7
     8
    10
    17
    19
    23
    24
    25
    27
    33
    34
    35
    36
    37
    43
    44
    45
    46
    47
    49
    51
    52
    54
    56
    58
    59
    60
    61
    63
    64
    68
    ];

angbest05 = [
   345
   165
   110
   155
   215
    80
   170
   140
   125
   260
   240
   175
   185
   260
   195
   170
   250
   165
   335
   165
   245
   230
    10
   250
   130
   140
   210
   115
   245
   120
   160
    20
    ];


%%
% modified, accepted by further visual checking
acctran04vis = [
    2004195   4.1784e+04   4.2666e+04 % 2004195,4.15e+4,4.60e+4;
    2004197   3.3831e+04   3.5500e+04
    2004197   4.0944e+04   4.2168e+04
    2004197   5.4704e+04   5.5897e+04
    2004197   5.6240e+04   5.7257e+04
    2004197   6.9920e+04   7.1142e+04
    2004197   7.7484e+04   7.9841e+04
    2004197   7.9919e+04   8.1483e+04
    2004197   8.3766e+04   8.4888e+04  % old 6, new 34, ACCEPTED
    2004197   8.5320e+04   8.6638e+04  % old 8, new 34, ACCEPTED
    2004198   5.9230e+03   9.7670e+03
    2004198   1.9998e+04   2.1646e+04
    2004198   4.2602e+04   4.3751e+04
    2004198   5.3476e+04   5.8080e+04
    2004198   6.2936e+04   6.3915e+04
    2004198   7.3510e+04   7.6059e+04
    2004198   8.4192e+04   8.5952e+04
    2004199   4.1370e+03   6.1680e+03
    2004199   4.2815e+04   4.3333e+04
    2004199   4.7777e+04   4.9156e+04
    2004199   4.9714e+04   5.0644e+04
    2004199   8.0831e+04   8.3462e+04
    2004200   1.1676e+04   1.5829e+04
    2004200   1.5884e+04   1.7155e+04
    2004200   1.9075e+04   2.0140e+04
    2004200   4.8079e+04   4.8696e+04
    2004201   8.3400e+02   1.6590e+03
    2004203   1.8492e+04   2.0806e+04
    ];


angbest04vis = [
   260
   155
   165
   220
   185
   195
   195
   195
   230
   270   
   235
   205
   230
    80
    65
   235
   250
   310
   250
   135
   195
   140
   170
    65
   215
   240
   260
   120
];

% modified, accepted by further visual checking
acctran05vis = [
    2005255   4.6069e+04   4.6979e+04
    2005255   5.8020e+04   5.9341e+04
    2005255   6.0351e+04   6.1587e+04
    2005255   6.7238e+04   6.8565e+04   % new 10, old 22, questionable direction, 155-->65
    2005256   3590         4987         % new 16, old 25, questionable direction, 175-->245
    2005256   8.0030e+03   9.8000e+03   % new 17, old 26, change ending time as old, ACCEPTED
    2005256   1.3382e+04   1.4376e+04
    2005256   2.0407e+04   2.2324e+04
    2005256   2.8260e+04   2.9202e+04   % new 24, old 28, change start time as 7.85, ACCEPTED
    2005256   2.9685e+04   3.1233e+04
    2005256   3.4786e+04   3.5268e+04
    2005256   5.2020e+04   5.3208e+04   % new 32, old 30, try a new range, ACCEPTED
    2005256   5.5548e+04   5.7052e+04   % new 33, change start time as 15.43, ACCEPTED
    2005256   6.0693e+04   6.4846e+04
    2005256   7.1496e+04   7.3944e+04   % new 35, old 31, change start time as 19.86, ACCEPTED
    2005256   7.6351e+04   7.6924e+04
    2005256   7.6984e+04   7.8694e+04
    2005257   6120         7776         % new 41, old 33, change start time as 1.7, ACCEPTED 
    2005257   1.5943e+04   1.8685e+04
    2005257   2.1040e+04   2.5395e+04
    2005257   2.5465e+04   2.8085e+04
    2005257   3.5320e+04   3.7628e+04
    2005257   3.8483e+04   4.5322e+04
    2005257   6.1419e+04   6.3758e+04
    2005257   7.3563e+04   8.2073e+04
    2005257   8.2129e+04   8.3410e+04
    2005258   2.4787e+04   2.6247e+04
    2005258   3.4970e+04   3.9993e+04
    2005259   1.9280e+03   4.0760e+03
    2005259   5.2020e+03   7.7440e+03
    2005259   3.7267e+04   4.1249e+04
    2005259   7.4225e+04   7.4753e+04
    2005260   2.6210e+03   3.4680e+03
    2005260   3.5740e+03   6.7080e+03
    2005260   9.3790e+03   1.0387e+04
    ];

angbest05vis = [
   345
   165
   110
   155
   175
   215
    80
   170
   105
   125
   260
   175
   255
   175
   180
   260
   195
   195
   170
   250
   165
   335
   165
   245
   230
    10
   250
   130
   140
   210
   115
   245
   120
   160
    20
   ];












