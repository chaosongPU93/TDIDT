function newloc = LfelocBostock(lon0,lat0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This is a function with default input to output the locations in lat and lon
% of the intersection between fam have detections and fam have locations,
% format: [fam dx dy lon lat dep], 6 cols
% 
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/01/04
% Last modified date:   2022/01/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('lon0',-123.585000);   % default value for loc of fam 002 is from direct Hypoinverse
defval('lat0',48.436667);

bosdir = ('/home/data2/chaosong/matlab/allan/BOSTOCK');
lfefnm = ('newlfeloc');
lfeloc = load(fullfile(bosdir, lfefnm));    
nlfe = size(lfeloc,1);

%%% set path and name
catapath= ('/home/data2/chaosong/matlab/allan/BOSTOCK/');
cataname = ('total_mag_detect_0000_cull_NEW.txt');
catalog = load(strcat(catapath, cataname));
%%% find the families that have LFE detections
famall = unique(catalog(:, 1));

%%% find the intersection between fam have detections and fam have locations
[~,~,ind] = intersect(famall, lfeloc(:, 1));
lfeuse = lfeloc(ind,:); % format: [fam lat lon dep], 4 cols

%convert to cart locations relative to loc of fam 002 
[dx,dy] = absloc2relaloc(lfeuse(:,3),lfeuse(:,2),lon0,lat0);

%rearrange the order to format: [fam dx dy lon lat dep], 6 cols
newloc = [lfeuse(:,1) dx dy lfeuse(:,3) lfeuse(:,2) lfeuse(:,4)];