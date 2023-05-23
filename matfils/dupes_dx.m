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

fam='002';
if isequal(fam,'002')
%   mapinputs=['mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat100.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s'];
%   mapoutputs=['mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.ell_1.5-1.else0.25nsat100.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5'];
%   mapinputs=['mapSTAS.UN.100sps.rec_-52.2513.5.else0nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0nsat100.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s'];
%   mapoutputs=['mapSTAS.UN.100sps.rec_-52.2513.5.else0nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0nsat100.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5'];
%   mapinputs=['mapSTAS.UN.100sps.rec_-52.2513.5.else0.25nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0.25nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0.25nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0.25nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0.25nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0.25nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0.25nsat100.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0.25nsat500.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s'];
%   mapoutputs=['mapSTAS.UN.100sps.rec_-52.2513.5.else0.25nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0.25nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0.25nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0.25nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0.25nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0.25nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0.25nsat100.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0.25nsat500.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5'];
%   mapinputs=['mapSTAS.UN.100sps.rec_-52.2513.5.else0.5nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0.5nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0.5nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0.5nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0.5nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0.5nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0.5nsat100.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
%              'mapSTAS.UN.100sps.rec_-52.2513.5.else0.5nsat500.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s'];
%   mapoutputs=['mapSTAS.UN.100sps.rec_-52.2513.5.else0.5nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0.5nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0.5nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0.5nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0.5nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0.5nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0.5nsat100.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
%               'mapSTAS.UN.100sps.rec_-52.2513.5.else0.5nsat500.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5'];
    mapinputs=['mapSTAS.UN.100sps.ell_3-1.5.else0nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_3-1.5.else0nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_3-1.5.else0nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_3-1.5.else0nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_3-1.5.else0nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_3-1.5.else0nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_3-1.5.else0nsat100.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s'];
    mapoutputs=['mapSTAS.UN.100sps.ell_3-1.5.else0nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_3-1.5.else0nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_3-1.5.else0nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_3-1.5.else0nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_3-1.5.else0nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_3-1.5.else0nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_3-1.5.else0nsat100.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5'];
    mapinputs=['mapSTAS.UN.100sps.ell_2.75-1.else0nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_2.75-1.else0nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_2.75-1.else0nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_2.75-1.else0nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_2.75-1.else0nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_2.75-1.else0nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s'];
    mapoutputs=['mapSTAS.UN.100sps.ell_2.75-1.else0nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_2.75-1.else0nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_2.75-1.else0nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_2.75-1.else0nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_2.75-1.else0nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_2.75-1.else0nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5'];
    mapinputs=['mapSTAS.UN.100sps.ell_3-1.25.else0nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_3-1.25.else0nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_3-1.25.else0nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_3-1.25.else0nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_3-1.25.else0nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_3-1.25.else0nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s'];
    mapoutputs=['mapSTAS.UN.100sps.ell_3-1.25.else0nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_3-1.25.else0nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_3-1.25.else0nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_3-1.25.else0nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_3-1.25.else0nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_3-1.25.else0nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5'];
    mapinputs=['mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s';
               'mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat100.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s'];
    mapoutputs=['mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat0.1.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat0.4.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat1...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat2...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat4...loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat10..loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5';
               'mapSTAS.UN.100sps.ell_3-1.25.diam0.3.else0nsat100.loff1.5.ccmin0.44.nponpa22_1.25-6.5-ms19-4s.onepk.5'];

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

offset=0.5;
ntotkept=0;
for ifile=1:nfile
    mapfil=load(mapinputs(ifile,:));
    ncol=size(mapfil,2);
    mapfil(:,16:19)=0;
    PGSIs=mapfil(:,2); %KLUDGY MINUS SIGN! Deleted after changing the sign of the offset in the synthetic seismogram.
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
%   10.0  -2.25   0.25    0.758   -0.11  1.231e+01    0.865      9.075   0.861   0.449   0.640  6.276e-01  2.031e+00  6.058e-03  8.832e-01  0.309  0.785  0.252  0.77  1.19   0.694   0.794   0.786
%   11.0  -2.25   0.25    0.753   -0.06  1.231e+01    0.865      9.075   0.899   0.449   0.640  6.276e-01  2.031e+00  6.058e-03  8.832e-01  0.309  0.785  0.252  0.69  1.12   0.703   0.789   0.766
%   15.0  -1.50   0.50    0.929    0.22  8.302e+00    0.968     16.450   1.011   0.858   0.567  1.422e-01  7.920e+00  4.701e-03  1.822e-01 -0.031  0.123  0.618  1.02  1.16   0.932   0.925   0.929

    end
%   mapfil(:,2:3)=2*mapfil(:,2:3);
    ntot=size(mapfil,1); %I think this must be the same as nevs
    mapfid = fopen(mapoutputs(ifile,:),'w');
    ndupes=0;
    i=0; j=0; kept=0; l=0;
    tim=mapfil(:,1);
    timpk=mapfil(:,8);
    xc=mapfil(:,4);
    while j < ntot-1
        i=i+1; %i is incremented
        j=i+1; %j is one more than i
        while timpk(j)-timpk(i) < offset && j < ntot %keep incrementing j until it is more than offset after i.
            j=j+1;
        end
        if j-i==1 && (tim(i)-timpk(i) > 5.97 || tim(i)-timpk(i) < -4.49) %I think this just counts isolated detections at the margins of a window; increments l if so. Now for 12s window.
            l=l+1;
        end
        if j-i==1 && tim(i)-timpk(i) < 5.97 && tim(i)-timpk(i) > -4.49  %If isolated and not at the margins, write it.
            fprintf(mapfid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %8.3f %10.3f %7.3f %7.3f %7.3f %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %5.2f %5.2f %8.3f %7.3f %7.3f\n',mapfil(i,:));
            kept=kept+1; 
            mapfilkeep(kept,:)=mapfil(i,:);
        else 
            ndupes=ndupes+1;
            xcmax=0.;
            for n=i:j-1
                if xc(n)> xcmax && tim(n)-timpk(n) < 5.97 && tim(n)-timpk(n) > -4.49 %search for the largest cc value that isn't at the margins.
                    xcmax=xc(n);
                    nkeep=n;
                end
            end
            if xcmax > 1.e-6 %if larger than the zero initial value, write it
                fprintf(mapfid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %8.3f %10.3f %7.3f %7.3f %7.3f %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %5.2f %5.2f %8.3f %7.3f %7.3f\n',mapfil(nkeep,:)');
                kept=kept+1;
                mapfilkeep(kept,:)=mapfil(nkeep,:);
            end
        end
        if j < ntot
            i=j-1;
        else
%         If last entries are a string of duplicates, they might appear twice.
%         check and take the highest x-correlation value.
            xcmax=0.;
            for n=i:j
                if xc(n) > xcmax && tim(n)-timpk(n) < 5.97 && tim(n)-timpk(n) > -4.49
                xcmax=xc(n);
                nkeep=n;
                end
            end
            fprintf(mapfid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %8.3f %10.3f %7.3f %7.3f %7.3f %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %5.2f %5.2f %8.3f %7.3f %7.3f\n',mapfil(nkeep,:));
            kept=kept+1;
            mapfilkeep(kept,:)=mapfil(nkeep,:);
        end
    end
    ntot
    ndupes
    kept
    nmargins=l
    mapfullfil(ntotkept+1:ntotkept+kept,:)=mapfilkeep(1:kept,:);
    ntotkept=ntotkept+kept;
end

