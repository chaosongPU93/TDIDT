% test the stability of templates at different stations when trying to gett
% the travel time difference between different stations using different
% main stations in cross-correlation

% set path
workpath = getenv('ALLAN');
temppath = strcat(workpath, '/templates/');
if ~exist(temppath, 'dir')
    mkdir(temppath);
end
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/PGCtrio');


mgcb = load(strcat(datapath, '/split_chao/','rot_para_11fam0.5-6.5_25_35MGCB_40sps7stas'));

twkb = load(strcat(datapath, '/split_chao/','rot_para_11fam0.5-6.5_25_35TWKB_40sps7stas'));

lzb = load(strcat(datapath, '/split_chao/','rot_para_11fam0.5-6.5_25_35LZB_40sps7stas'));

% convert them to a same reference station, for example, LZB
for i=1:11
    twkb((i-1)*7+1:i*7, 4) = twkb((i-1)*7+1:i*7, 4)-twkb((i-1)*7+4, 4);
    mgcb((i-1)*7+1:i*7, 4) = mgcb((i-1)*7+1:i*7, 4)-mgcb((i-1)*7+4, 4);
end

diffmgcb = mgcb-lzb;
difftwkb = twkb-lzb;
