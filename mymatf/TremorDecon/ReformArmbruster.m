function newcat = ReformArmbruster(lon0,lat0,remake)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This is a function with default input to rearrange the tremor catalog from
% John Armbruster, format: [yyyy mm dd sec dx dy lon lat dep], 9 cols
% 
% Validating the catalog by plotting the differential time of events on the
% same date, you will see that the minimum is about 8 s, and the separations
% roughly are several times of 8 s. This proves that the catalog is indeed 
% what is claimed in Armbruster et al. 2014 JGRSE paper. The detection windows
% are 150-s long offset every 8 s. Check the page 3 of that paper. 
% 
% After comparing '40sps.pgsssi.hyp' with John's catalog on this web page:
% 'https://www.ldeo.columbia.edu/~armb/TREM/3_station.txt', it looks like 
% the time in '40sps.pgsssi.hyp' is the event origin time, but in pratice,
% maybe the arrival time at PGC is more useful.
%
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2021/12/22
% Last modified date:   2021/12/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('lon0',-123.585000);   % default value for loc of fam 002 is from direct Hypoinverse
defval('lat0',48.436667);
defval('remake',0); % if read the original catalog and compile it again

workpath = getenv('ALLAN');

if remake
  johncat = load(strcat(workpath,'/BOSTOCK/40sps.pgsssi.hyp'));

  yr = floor(johncat(:,1)/10000);
  yyyy = 2000+yr;
  mm = floor((johncat(:,1)-yr*10000)/100);
  dd = johncat(:,1)-yr*10000-mm*100;

  hr = floor(johncat(:,2)/100);
  mn = johncat(:,2)-hr*100;
  time = hr*3600+mn*60+johncat(:,3);

  lat = johncat(:,4)+johncat(:,5)/60;
  lon = -(johncat(:,6)+johncat(:,7)/60);
  dep = johncat(:,8);

  %convert to cart locations relative to loc of fam 002 
  [dx, dy] = absloc2relaloc(lon,lat,lon0,lat0);

  %rearrange the order to format: [yyyy mm dd sec dx dy lon lat dep], 9 cols
  newcat = [yyyy mm dd time dx dy lon lat dep];
  
  %save to file
  fid = fopen(strcat(workpath,'/BOSTOCK/40sps.pgsssi.hyp.compiled'),'w+');
  fprintf(fid,'%d %d %d %e %f %f %f %f %f \n', newcat'); 
  fclose(fid);

else
  newcat = load(strcat(workpath,'/BOSTOCK/40sps.pgsssi.hyp.compiled'));
    
end  
  
% keyboard  
  
  
  
  
  
  
  
  
  
  
