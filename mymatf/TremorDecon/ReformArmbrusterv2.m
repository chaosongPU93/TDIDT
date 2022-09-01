function newcat = ReformArmbrusterv2(lon0,lat0,remake)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This is a function with default input to rearrange the tremor catalog from
% John Armbruster, format: [yyyy mm dd sec dx dy lon lat dep], 9 cols
% 
% Different from 'ReformArmbruster.m', this catalog comes from web page:
% 'https://www.ldeo.columbia.edu/~armb/TREM/3_station.txt'.
% The catalog is indeed what is claimed in Armbruster et al. 2014 JGRSE paper.
% The detection windows are 150-s long offset every 8 s. Check the page 3 of 
% that paper. 
% it looks like the time in that file is the start of the center of the 
% detecting window at station PGC, in pratice, maybe more useful.
% 
% --Comparing the 2 catalog via finding the detection that is closest to fam 
%   002, the time difference between the same detection (ind=2822) in the two 
%   catalogs is about 11.93 s, which is roughly the travel time needed for
%   a source at 002 received at station PGC. So yes, the guess is right, 
%   catalog in this version denotes 'arrival time', while the catalog in earlier
%   version denotes 'origin time' which is corrected for travel time, but still
%   not sure if it is start of the center of a detection window.
%
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/03/29
% Last modified date:   2022/03/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('lon0',-123.585000);   % default value for loc of fam 002 is from direct Hypoinverse
defval('lat0',48.436667);
defval('remake',0); % if read the original catalog and compile it again

workpath = getenv('ALLAN');

if remake
  fid = fopen(fullfile(workpath,'/BOSTOCK/JohnArmb_3station_40sps.txt'),'r');
  
  format='%s %s %d %d %f %f %d %d %d %d %f \n';
  % format='%s \n';
  datacell = textscan(fid,format,'HeaderLines',3);
  
  datestr = datacell{1,1};  % cell array
  timestr = datacell{1,2};
  lat = datacell{1,5};  % double array
  lon = datacell{1,6};
  
  time = [];
  for i = 1: size(datestr,1)
    str = datestr{i};
    yyyy = str2double(str(1:2))+2000;
    jday = str2double(str(4:end));
    a = jul2dat(yyyy,jday);
    mm = a(1);
    dd = a(2);
    
    str = timestr{i};
    sec = str2double(str(1:2))*3600 + str2double(str(4:5))*60 + str2double(str(7:8));
    time = [time; yyyy mm dd sec];
  end
  
  %interpolation to get the depth according to the old slab file
  olds = load('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssi/xyz.grid');
  F = scatteredInterpolant(-olds(:,1),olds(:,2),olds(:,3),'linear','linear');
  dep = F(lon,lat);
  
  %convert to cart locations relative to loc of fam 002
  [dx, dy] = absloc2relaloc(lon,lat,lon0,lat0);
  
  %rearrange the order to format: [yyyy mm dd sec dx dy lon lat dep], 9 cols
  newcat = [time dx dy lon lat dep];

  %save to file
  fid = fopen(strcat(workpath,'/BOSTOCK/JohnArmb_3station_40sps_compiled.txt'),'w+');
  fprintf(fid,'%d %d %d %e %f %f %f %f %f \n', newcat'); 
  fclose(fid);

else
  newcat = load(strcat(workpath,'/BOSTOCK/JohnArmb_3station_40sps_compiled.txt'));
    
end  

% keyboard







