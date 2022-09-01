%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS script is used to get the locations (map coordinates) of detections
% by interpolate existing locations grid through hypoinverse from John
% Armbruster.
%
% Modified from grid_cat.m written by Prof. Rubin  
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2019/05/23
% Last modified date:   2019/05/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all
rads=pi/180.;
erad=6372.028;  % radius of earth
lat0=48.0+26.32/60.;    % lat of the event inverted from (0,0), off12,off13
lon0=123.0+35.07/60.;   % lon of the event inverted from (0,0), off12,off13
srad=erad*cos(lat0*rads);   % radius of main station

fam='002'
A=load('p085p020');
mm=25; nn=25;
dlat=reshape(A(:,3),mm,nn); %degrees
mlat=reshape(A(:,4),mm,nn); %minutes
lat=dlat+mlat/60.;
%%% Vq = interp2(V,K) returns the interpolated values on a refined grid 
%%% formed by repeatedly halving the intervals K times in each dimension.
%%% This results in 2^K-1 interpolated points between sample values.
lat2=interp2(lat,3);    % 24*(2^3-1)+25=193
dlon=reshape(A(:,5),mm,nn); %degrees
mlon=reshape(A(:,6),mm,nn); %minutes
lon=dlon+mlon/60.;
lon2=interp2(lon,3);
dep=reshape(A(:,7),mm,nn);

%%%%%%%%% algorithm %%%%%%%%%%%%%%%%
% arr2 = arr1 - off12 + 86(centPGSS)
% arr3 = arr1 - off13 + 20(centPGSI)
%  ||
% arr2 - arr1 = 86(centPGSS) - off12 == 1st col in A
% arr3 - arr1 = 20(centPGSI) - off13 == 2nd col in A
%  ||
% off12 = 86(centPGSS) - (arr2 - arr1)
% off13 = 20(centPGSI) - (arr3 - arr1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
centPGSS=86
centPGSI=20
offPGSS=reshape(centPGSS-A(:,1),mm,nn); %PGSS is uniform in columns; changes across rows.
offPGSI=reshape(centPGSI-A(:,2),mm,nn); %PGSI is uniform in rows; changes down columns.
offPGSS2=interp2(offPGSS,3); %scalar is number of times grid is halved
offPGSI2=interp2(offPGSI,3); %scalar is number of times grid is halved
offPGSS1=offPGSS2(1,:); %a 1-d array of offsets at 1/4 sample
offPGSI1=offPGSI2(:,1);
mPGSS=length(offPGSS1); %
nPGSI=length(offPGSS2); %This and the previous line look iffy!  Should be offPGSS2 and offPGSI2?  Probably works only b.c. grid is square.
% In p085p020, order is PGSS PISI.  PGSI changes more rapidly.
% in output of wiggledaywig2, order is [timswin(n) xmaxPGSIntmp(n) xmaxPGSSntmp(n) ...

workpath = getenv('ALLAN');
datapath = strcat(workpath, '/PGCtrio/MAPS');
cd(datapath);
system('rm -f mapall.002.loff2.1.ccmin0.44.nponpa22.ms19_1.25-6.5_4s40sps');
system('cat map*.002.loff2.1.ccmin0.44.nponpa22.ms19_1.25-6.5_4s40sps > mapall.002.loff2.1.ccmin0.44.nponpa22.ms19_1.25-6.5_4s40sps');

locfile = load('mapall.002.loff2.1.ccmin0.44.nponpa22.ms19_1.25-6.5_4s40sps');

ncol=size(locfile,2);
PGSIs=locfile(:,2);
PGSSs=locfile(:,3);
nevt=length(PGSIs);
nhits=zeros(nPGSI,mPGSS);
evtlon = zeros(nevt,1);
evtlat = zeros(nevt,1);
hitcount = zeros(nevt,1);
for i=1:nevt
    if PGSIs(i)>=min(offPGSI1) && PGSIs(i)<=max(offPGSI1) && PGSSs(i)>=min(offPGSS1) && PGSSs(i)<=max(offPGSS1) %if w/in Armbruster's grid
        for j=1:nPGSI
            if (abs(PGSIs(i)-offPGSI1(j))) <= 1.e-6
                PGSImatch=j;                            %This just assigns PGSI offset
            end
        end
        for j=1:mPGSS
            if (abs(PGSSs(i)-offPGSS1(j))) <= 1.e-6
                PGSSmatch=j;                            %This just assigns PGSS offset
            end
        end
        evtlat(i) = lat2(PGSImatch,PGSSmatch);
        evtlon(i) = -lon2(PGSImatch,PGSSmatch);
        nhits(PGSImatch,PGSSmatch)=nhits(PGSImatch,PGSSmatch)+1;
        locfile(i,14)=PGSImatch;
        locfile(i,15)=PGSSmatch;
    else
        evtlat(i) = 91;
        evtlon(i) = 181;
    end
    
end

for i=1:nevt
    hitcount(i) = nhits(locfile(i,14),locfile(i,15)); 
end

evtinfo_tmp = [hitcount evtlat evtlon];
evtinfo_tmp = sortrows(evtinfo_tmp, [2,3]);
evtinfo = unique(evtinfo_tmp, 'rows');

fid = fopen('/home/data2/chaosong/Seisbasics/hypoinverse/Chaotest/eventloc.002.interp', 'w+');
fprintf(fid,'%d %.6f %.6f \n',evtinfo');
fclose(fid);




















