%One more column (mean of running average CC) added 7/6/21
%
%This one cleans up some glitches 2/22/22 - eliminates detections occuring w/in
%offset of kept detections EVEN if they occur > offset from the first in a
%series of overlapping detections. Keeps the last entry. Fixes a bug in
%treating a sequence at the end of the daily catalog.
%
%Also keeps the largest CC detection, even one at the margins of a window.
%For when locations are more important than statistics of the "arrival".
%
 clear all
 close all
 scrsz=get(0,'ScreenSize');
 wid=scrsz(3);
 hite=scrsz(4);
 scrat=wid/hite; 
 rads=pi/180.;
 erad=6372.028;
 lat0=48.0+26.32/60.;
 lon0=123.0+35.07/60.;
 srad=erad*cos(lat0*rads)
 
 fam='002';
 A=load('p085p020');
 mm=25; nn=25;
 dlat=reshape(A(:,3),mm,nn); %degrees
 mlat=reshape(A(:,4),mm,nn); %minutes
 lat=dlat+mlat/60.;
 lat2=interp2(lat,3); 
 dlon=reshape(A(:,5),mm,nn); %degrees
 mlon=reshape(A(:,6),mm,nn); %minutes
 lon=dlon+mlon/60.;
 lon2=interp2(lon,3); 
 dep=reshape(A(:,7),mm,nn);

 centPGSS=86
 centPGSI=20
 offPGSS=reshape(centPGSS-A(:,1),mm,nn); %PGSS is uniform in columns; changes across rows.
 offPGSI=reshape(centPGSI-A(:,2),mm,nn); %PGSI is uniform in rows; changes down columns.
 offPGSS2=interp2(offPGSS,3); %scalar is number of times grid is halved
 offPGSI2=interp2(offPGSI,3); %scalar is number of times grid is halved
 offPGSS1=offPGSS2(1,:);
 offPGSI1=offPGSI2(:,1);
 mPGSS=length(offPGSS1); %
 nPGSI=length(offPGSS2); %This and the previous line look iffy!  Should be offPGSS2 and offPGSI2?  Probably works only b.c. grid is square.
% In p085p020, order is PGSS PISI.  PGSI changes more rapidly.
% in output of wiggledaywig2, order is [timswin(n) xmaxPGSIntmp(n) xmaxPGSSntmp(n) ...

offset=0.5;
suff=num2str(offset);
fam='002';
if isequal(fam,'002')
%    mapinputs='map2003.063.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s';
    mapinputs=['map2003.062.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s';
%                'map2003.063.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s';
%                'map2004.196.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s';
%                'map2004.197.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s';
%                'map2004.198.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s';
%                'map2004.199.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s';
%                'map2005.254.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s';
%                'map2005.255.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s';
%                'map2005.256.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s'
               ];
    mapoutputs=mapinputs;
    for n=1:size(mapoutputs,1)
        mapoutputs(n,size(mapinputs,2)+1:size(mapinputs,2)+6+length(suff))=['.onepk',suff];
    end
%     mapoutputs='map2003.063.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk.5'
%     mapoutputs=['map2003.062.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk.5';
%                 'map2003.063.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk.5';
%                 'map2004.196.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk.5';
%                 'map2004.197.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk.5';
%                 'map2004.198.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk.5';
%                 'map2004.199.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk.5';
%                 'map2005.254.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk.5';
%                 'map2005.255.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk.5';
%                 'map2005.256.002.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms18-4s.onepk.5'];
    sta=['PGC ';
         'SSIB';
         'SILB'];
elseif isequal(fam,'068')
    mapinputs=['map2004.198.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s';
               'map2004.199.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s';
               'map2004.200.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s';
               'map2004.201.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s';
               'map2005.256.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s';
               'map2005.257.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s';
               'map2005.258.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s';
               'map2005.259.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s';
               'map2005.260.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s';
               'map2005.261.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s'];
    mapoutputs=['map2004.198.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s.onepk.5';
                'map2004.199.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s.onepk.5';
                'map2004.200.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s.onepk.5';
                'map2004.201.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s.onepk.5';
                'map2005.256.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s.onepk.5';
                'map2005.257.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s.onepk.5';
                'map2005.258.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s.onepk.5';
                'map2005.259.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s.onepk.5';
                'map2005.260.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s.onepk.5';
                'map2005.261.068.loff1.5.ccmin0.4.nponpa22_1.25-8-ms4-4s.onepk.5'];
    sta=['TWKB '
         'LZB  '
         'MGCB '];
end
nfile=size(mapinputs,1);

%offset=0.5;
ntotkept=0;
for ifile=1:nfile
    mapfilwin=load(mapinputs(ifile,:));
    mapfil=sortrows(mapfilwin,8); %B.C. sometimes the strongest 1/2 sec is not chronological...
    ncol=size(mapfil,2);
    mapfil(:,16:19)=0;
%   mapfil(:,13)=mapfil(:,15);  %(:,15) was optimal/orthogonal;  %obsolete; wdwSTA now writes this into columm 13.
    PGSIs=mapfil(:,2);
    PGSSs=mapfil(:,3);
    nevs=length(PGSIs);
    %nhits=zeros(nPGSI,mPGSS);
    for i=1:nevs
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
        %nhits(PGSImatch,PGSSmatch)=nhits(PGSImatch,PGSSmatch)+1;
        dy=rads*(lat2(PGSImatch,PGSSmatch)-lat0)*erad;
        dx=-rads*(lon2(PGSImatch,PGSSmatch)-lon0)*srad; %(minus bcause 123 is -123)
        mapfil(i,14)=dy;                               %This assigns north coord
        mapfil(i,15)=dx;                               %This assigns east coord
        angrot=-45*pi/180;
        spot=complex(dx,dy);
        spotrot=spot*exp(1i*angrot);
        mapfil(i,16)=real(spotrot);                    %This assigns pseudo-east coord
        mapfil(i,17)=imag(spotrot);                    %This assigns pseudo-north coord        
        %mapfil(i,18)=lat2(PGSImatch,PGSSmatch);        %This is not used
        %mapfil(i,19)=lon2(PGSImatch,PGSSmatch);        %This is not used
        mapfil(i,18)=PGSImatch;
        mapfil(i,19)=PGSSmatch;
    %   if locfile(i,16) >= -6 && locfile(i,16) <= 3.2 && locfile(i,17) >= -3.5 && locfile(i,17) <= 4.5
    %   fprintf(fid,'%9.1f %6.2f %6.2f %8.3f %7.2f %8.3f %8.3f %10.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',locfile(i,:)); %transpose?
%   20.5  -0.50   4.25    0.491    1.35  1.516e-01  9.961e-02     19.075   0.294   0.356   0.377  5.178e-01  3.306e-01  1.794e-03  1.644e-03   0.154   0.023   0.202  0.79  0.71   0.433   0.464   0.575
%   43.5   4.75   1.00    0.442    0.70  1.078e-01  5.525e-02     44.200   0.350   0.214   0.307  1.307e-01  2.339e-01  9.596e-04  1.388e-03   0.523   0.334   0.277  0.43  0.68   0.385   0.398   0.545
%   72.5  -0.75  -1.00    0.559   -1.07  2.572e-01  2.087e-01     71.200   0.478   0.490   0.399  1.427e-01  3.509e-01  2.357e-03  2.133e-03  -0.300  -0.097   0.583  0.65  0.87   0.556   0.565   0.555
%   1      2      3       4        5     6          7             8        9       10      11     12         13         14         15          16      17      18     19    20
%  9.1f    6.2f    6.2f    8.3f    7.2f    10.3e    8.3f          10.3f    7.3f    7.3f    7.3f    10.3e    10.3e        10.3e      10.3e      7.3f    7.3f    7.3f    5.2f 5.2f
    end
%   mapfil(:,2:3)=2*mapfil(:,2:3);
    ntot=size(mapfil,1); %I think this must be the same as nevs
    mapfid = fopen(mapoutputs(ifile,:),'w');
    ndupes=0;
    i=0; j=0; kept=0; l=0;
    tim=mapfil(:,1);
    timpk=mapfil(:,8);
    xc=mapfil(:,4);
    while j < ntot%-1  %Changed 2/21/22
        i=i+1; %i is incremented
        j=i+1; %j is one more than i
        while abs(timpk(j)-timpk(i)) < offset && j < ntot %keep incrementing j until it is more than offset after i. IT IS POSSIBLE FOR timpk(j) < timpk(i)!!!
            j=j+1;
        end
        if j-i==1 && (tim(i)-timpk(i) > 1.97 || tim(i)-timpk(i) < -1.49) %I think this just counts isolated detections at the margins of a window; increments l if so. 
            l=l+1;
        end
        if j-i==1 && abs(timpk(j)-timpk(i)) > offset  %If isolated, write it.  The second half guards against the first being true b.c. j=ntot
            fprintf(mapfid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %8.3f %10.3f %7.3f %7.3f %7.3f %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %5i %5i %10.3e %6.3f %6.3f %6.3f %7.3f\n',mapfil(i,:));
            kept=kept+1; 
            mapfilkeep(kept,:)=mapfil(i,:);
        else 
            ndupes=ndupes+1;
            xcmax=0.;
            if abs(timpk(j)-timpk(i)) < offset  %if this is true, j=ntot
                jlast=j
                tim(j)
            else
                jlast=j-1;
            end
            for n=i:jlast
                %Use the first below if doing statistics on arrival.  Use the second
                %if you just want the maxx cc value of a string of detections.
                %if xc(n)> xcmax && tim(n)-timpk(n) < 1.97 && tim(n)-timpk(n) > -1.49 %search for the largest cc value that isn't at the margins.
                if xc(n)> xcmax %search for the largest cc value.
                    xcmax=xc(n);
                    nkeep=n;
%                     searchformaxcc=n
%                     searchformaxcc=tim(n)
                end
            end
            if xcmax > 1.e-6 %if larger than the zero initial value, write it
                fprintf(mapfid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %8.3f %10.3f %7.3f %7.3f %7.3f %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %5i %5i %10.3e %6.3f %6.3f %6.3f %7.3f\n',mapfil(nkeep,:)');
                kept=kept+1;
                mapfilkeep(kept,:)=mapfil(nkeep,:);
            end
        end
        i=j-1; %if the detection was a singleton, i is now that event; if it was part of a cluster, i is now 
               %the last of that cluster. It will be incremented by one at the start of the loop above.
        while abs(timpk(i+1)-mapfilkeep(kept,8)) < offset && i < ntot-1 %keep incrementing j until it is more than offset after i. IT IS POSSIBLE FOR timpk(j) < timpk(i)!!!
            i=i+1;
        end
        if j==ntot %check to see if the last entry is a singleton.  Could be more efficient to do this after the next "end".
            if abs(timpk(j)-timpk(j-1)) > offset
                fprintf(mapfid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %8.3f %10.3f %7.3f %7.3f %7.3f %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %5i %5i %10.3e %6.3f %6.3f %6.3f %7.3f\n',mapfil(j,:));
                kept=kept+1; 
                mapfilkeep(kept,:)=mapfil(j,:);
            end
        end
    end
    ntot
    ndupes
    kept
    nmargins=l
    mapfullfil(ntotkept+1:ntotkept+kept,:)=mapfilkeep(1:kept,:);
    ntotkept=ntotkept+kept;
end

