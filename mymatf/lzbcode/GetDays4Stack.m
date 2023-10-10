function [timoffrot,tempoffs] = GetDays4Stack(fam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to get the days and corresponding approximate zero
% crossings for stacking (all 7 stations), based on 40 sps, station order 
% is ['PGC 'SSIB ''SILB ' 'LZB  ' 'TWKB ' 'MGCB ' 'KLNB '];
% 
% when the plot shows that the dipole is later than 400, then plus the same
% amount to the tempoffs.
%
%
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/06/25
% Last modified date:   2019/06/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isequal(fam,'002')  %DON'T FORGET additional family-specific delcarations around line 156
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
%     % tempoffs are the zero-crossings at E component before rotation
    tempoffs=[1211 1297 1231 1221 1207 1186 1207]; %these are the zero crossings


elseif isequal(fam, '043')  % use 043 instead of 013, but with error
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1290 1374 1349 1195 1210 1208 1254]; %these are the zero crossings
    
elseif isequal(fam, '141')  % use 141 instead of 025, but with error
        
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1318 1396 1391 1182 1206 1212 1269]; %these are the zero crossings
    
elseif isequal(fam, '047')  % use 047 instead of 028
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1240 1331 1274 1214 1213 1197 1227]; %these are the zero crossings
    
elseif isequal(fam, '010')  % use 010 instead of 056
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[10];
    timoffrot(ind, :)=[];
    
    tempoffs=[1288 1391 1296 1179 1207 1200 1203]; %these are the zero crossings
    
elseif isequal(fam, '144')  % use 144 instead of 084
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1342 1411 1416 1191 1218 1226 1285]; %these are the zero crossings
    
    
elseif isequal(fam, '099')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[8];
    timoffrot(ind, :)=[];
    
    tempoffs=[1315 1398 1378 1179 1206 1211 1268]; %these are the zero crossings
    
elseif isequal(fam, '068')  % use 068 instead of 115   
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1281 1347 1352 1194 1202 1202 1242]; %these are the zero crossings
    
elseif isequal(fam, '125')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1262 1340 1316 1200 1205 1198 1235]; %these are the zero crossings for 002: PGC,SSIB,SILB.
    
elseif isequal(fam, '147')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1309 1379 1383 1191 1207 1212 1262]; %these are the zero crossings for 002: PGC,SSIB,SILB.
    
elseif isequal(fam, '017')    % use 017 instead of 149
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[2 7];
    timoffrot(ind, :)=[];
    
    tempoffs=[1311 1379 1389 1186 1204 1210 1262]; %these are the zero crossings for 002: PGC,SSIB,SILB.

    
% added 2020/02/12, for fam 243 with PGC trio 
elseif isequal(fam, '243')
    
    % generate unique dates matrix that in that family
%     timoffrot= [2003 062; 
%                 2003 063;
%                 2004 196;
%                 2004 197;
%                 2004 198;
%                 2004 199;
%                 2005 254;
%                 2005 255;
%                 2005 256];
    timoffrot = Readbostock(fam);
    
    tempoffs=[1241 1300 1270 1278 1262 1242 1238]; %these are the zero crossings
    

elseif isequal(fam, '001')    % use 017 instead of 149
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1267 1393 1366 1171 1204 1208 1266]; %these are the zero crossings 
    
    
elseif isequal(fam, '019')  % use 047 instead of 028
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1258 1344 1293 1194 1205 1196 1227]; %these are the zero crossings
    
elseif isequal(fam, '021')  % use 047 instead of 028
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1260 1349 1300 1191 1204 1196 1231]; %these are the zero crossings 1253
    
elseif isequal(fam, '045')  % use 010 instead of 056
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1300 1404 1355 1179 1218 1217 1289]; %these are the zero crossings    
    
elseif isequal(fam, '076')  % use 010 instead of 056
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1300 1403 1365 1168 1203 1196 1258]; %these are the zero crossings  

elseif isequal(fam, '176')  % use 010 instead of 056
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1295 1400 1353 1166 1200 1204 1254]; %these are the zero crossings
    
elseif isequal(fam, '015')  % use 141 instead of 025, but with error
        
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1337 1412 1415 1170 1206 1218 1284]; %these are the zero crossings

elseif isequal(fam, '158')  % use 141 instead of 025, but with error
        
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1331 1406 1410 1168 1207 1219 1290]; %these are the zero crossings
    
elseif isequal(fam, '234')  % use 141 instead of 025, but with error
        
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1331 1406 1410 1168 1207 1219 1290]; %these are the zero crossings    
    tempoffs=tempoffs+148;

elseif isequal(fam, '231')  % use 141 instead of 025, but with error
        
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1477 1552 1540 1297 1341 1356 1424]; %these are the zero crossings
    
elseif isequal(fam, '006')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1325 1399 1383 1174 1206 1213 1274]; %these are the zero crossings
    
elseif isequal(fam, '253')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1306 1393 1285 1366 1318 1319 1323]; %these are the zero crossings

elseif isequal(fam, '036')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1171 1263 1134 1230 1215 1188 1194]; %these are the zero crossings

elseif isequal(fam, '034')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1181 1284 1145 1223 1214 1188 1201]; %these are the zero crossings

elseif isequal(fam, '061')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1129 1203 1048 1238 1210 1182 1174]; %these are the zero crossings

elseif isequal(fam, '023')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1099 1158 1013 1235 1207 1175 1155]; %these are the zero crossings

elseif isequal(fam, '055')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1096 1156 1012 1220 1195 1170 1146]; %these are the zero crossings

elseif isequal(fam, '026')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1125 1188 1055 1260 1217 1187 1172]; %these are the zero crossings

elseif isequal(fam, '251')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1315 1365 1242 1463 1417 1387 1363]; %these are the zero crossings

elseif isequal(fam, '162')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1087 1136 1024 1235 1200 1160 1135]; %these are the zero crossings

elseif isequal(fam, '240')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1228 1267 1219 1341 1315 1289 1253]; %these are the zero crossings

elseif isequal(fam, '250')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1251 1288 1223 1366 1345 1326 1280]; %these are the zero crossings

elseif isequal(fam, '255')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1245 1247 1310 1252 1224 1227 1218]; %these are the zero crossings

elseif isequal(fam, '012')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1240 1222 1319 1227 1200 1208 1206]; %these are the zero crossings

elseif isequal(fam, '065')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1260 1304 1331 1234 1226 1216 1230]; %these are the zero crossings

elseif isequal(fam, '246')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1260 1346 1280 1270 1258 1238 1256]; %these are the zero crossings
    
else
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    
    tempoffs=[1315 1398 1378 1179 1206 1211 1268]; %these are the zero crossings

end

% 
% stas=['PGC  '
%     'SSIB '
%     'SILB '
%     'LZB  '
%     'TWKB '
%     'MGCB '
%     'KLNB '];















