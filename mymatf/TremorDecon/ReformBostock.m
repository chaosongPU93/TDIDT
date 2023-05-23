function [newcat,newcatfull] = ReformBostock(lon0,lat0,remake,bostname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This is a function with default input to rearrange the LFE catalog from
% Michael Bostock, format: [fam yyyy mm dd sec dx dy lon lat dep magnitude
% number-of-stations], 12 cols. the first output is the LFEs that have 
% information of location, the second contains all
% 
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2021/12/22
% Last modified date:   2021/12/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('lon0',-123.585000);   % default value for loc of fam 002 is from direct Hypoinverse
defval('lat0',48.436667);
defval('remake',0); % if read the original catalog and compile it again
defval('bostname', '/BOSTOCK/total_mag_detect_0000_cull_NEW.txt');

workpath = getenv('ALLAN');

if remake

%   bostname = ('/BOSTOCK/total_mag_detect_0000_cull_NEW.txt');
  
  %%% 6 cols: 
  %%%     family ID, date (ymmdd), hour, sec, magnitude, number of stations that include that event
  %%%     in their template or magnitude determination or something
  % 1. family ID
  % 2. date: 2-digit year, 2-digit month, 2-digit day
  % 3. hour of thay day, minus 1 to real hour
  % 4. secs (<3600)
  % 5. magnitude 
  % 6. number of stations that include that event in their template or magnitude determination or something
  catalog = load(strcat(workpath, bostname));

  yy = floor(catalog(:,2)/10000);
  mm = floor((catalog(:,2)-yy*10000)/100);
  dd = catalog(:,2)-yy*10000-mm*100;
  yyyy = yy+2000;
  hh = catalog(:,3)-1;
  sec = hh*3600+catalog(:,4);

  temp = [catalog(:,1) yyyy mm dd sec catalog(:,5) catalog(:,6)]; % 7 cols

  %get the location of lfe fams that have detections
  lfeloc = LfelocBostock(lon0,lat0); % format: [fam dx dy lon lat dep], 6 cols
  %note that fam have detections and fam have locations are NOT the same, so some columns could be 0
  fam = unique(lfeloc(:,1));
  
  %load the time correction for Michael's catalog time
  timecor = load(strcat(workpath,'/BOSTOCK/bostcattime_correction'));
  timecol = 2;  %using zero-crossing
  
  for i=1:length(fam)
    ind = find(temp(:,1)==fam(i));
    temp(ind,8:12) = repmat(lfeloc(lfeloc(:,1)==fam(i),2:end), [length(ind),1]);

    %correction to the time in sec according to fam ID, should point to the rough zero-crossing of LFE
    %although the order of fam are the same, but it is now done for safety
    offsec = timecor(timecor(:,1)==fam(i),timecol);
    
    temp(ind,5)=temp(ind,5)+offsec; % apply time correction ONLY for fams that have locations

  end
  
  %turn zeros to nans for those fams that do not have locations
  indnan = find(sum(temp(:,8:12),2)==0);
  temp(indnan,8:12) = nan(length(indnan),5);

  %rearrange the order to format: [fam yyyy mm dd sec dx dy lon lat dep magnitude
  % number-of-stations], 12 cols
  newcatfull = [temp(:,1:5) temp(:,8:12) temp(:,6:7)];

  %get rid of fams without locations
  ind = setdiff(1:size(temp,1), indnan);
  newcat = newcatfull(ind,:);

%   %save to file
%   fid = fopen(strcat(workpath,'/BOSTOCK/total_mag_detect_0000_cull_NEW_compiled.txt'),'w+');
%   fprintf(fid,'%d %d %d %d %e %f %f %f %f %f %f %d \n', newcat'); 
%   fclose(fid);
%   fid = fopen(strcat(workpath,'/BOSTOCK/total_mag_detect_0000_cull_NEW_compiledfull.txt'),'w+');
%   fprintf(fid,'%d %d %d %d %e %f %f %f %f %f %f %d \n', newcatfull'); 
%   fclose(fid);
  
else
  newcat = load(strcat(workpath,'/BOSTOCK/total_mag_detect_0000_cull_NEW_compiled.txt'));
  newcatfull = load(strcat(workpath,'/BOSTOCK/total_mag_detect_0000_cull_NEW_compiledfull.txt'));
  
end
  
% keyboard






