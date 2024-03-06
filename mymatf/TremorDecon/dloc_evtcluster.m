function [dloc,dt,k,dloc_spl]=dloc_evtcluster(impcluster,implocclus,sps,ftrans,m,n,timetype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [impcluster,impuni,ncluster,impclusterk,impunii]=evtcluster(nbst,imp,nsrc,m,dtcut,sps,timetype)
%
% This function is to obtain the location and time difference between events
% in the cluster that is defined by the first (N-m) and last (N) event in the 
% cluster.  Each cluster contains m+1 events, and you can ask the difference
% between event pair N and N-n inside each cluster, where n<=m. For example, 
% you can target all consecutive events in doublets, triplets, etc., or N and
% N-2 in triplets, quaduplets, quintuplets, etc.
%
% 2024/02/06, mapping from samples to relation location already contains error,
% so when you try to obtain the location difference, or distance, start with 
% the difference in samples first. For example, if you ask the difference in 
% certain directions, first obtain the difference in off12 and off13, then map
% the diff to location, then you can project along any direction to get the
% loc difference along that direction, or the distance (abs of loc diff) in
% that direction, etc. In all, it is NOT recommended to use output 'dloc' 
% unless a lof of data points would average out the error. 
% 
% 2024/02/15, looks like if you need the loc difference in map-view, you have
% to first invert time offsets to map locations, then obtain the difference.
% Despite the possible location error propagation, you need it. Since if you 
% start from difference in samples, then transform them to map locations,
% you are essentially assuming one of the two sources is located at the 
% origin, which is not the case in reality, even if the actual difference in
% the two ways may be small.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/01/23
% Last modified date:   2024/01/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('timetype','tarvl');

ncluster = round(size(impcluster,1)/(m+1));

k=0;
dt=[];
dloc=[];
dloc_spl=[];
for i = 1: ncluster
  ist=(i-1)*(m+1)+1;
  ied=i*(m+1);
  implocplt=implocclus(ist:ied,:);
  impplt=impcluster(ist:ied,:);
  if strcmp(timetype,'tori')
    %return the origin time, sorted srcs, and indices
    [tevt,impplt,indsort]=tarvl2tori(impplt,sps,ftrans,1);
    implocplt = implocplt(indsort, :);
  else
    tevt=impplt(:,1);
  end

%   keyboard
  dti = diffcustom(tevt,n,'forward');
  dloci = diffcustom(implocplt(:,1:2),n,'forward'); %relative loc in km
  dloci_spl = diffcustom(impplt(:,7:8),n,'forward'); %time offset in samples
  dt = [dt; dti];
  dloc = [dloc; dloci];
  dloc_spl = [dloc_spl; dloci_spl];
  k=k+size(dloci,1);

  % for j = 1:m
  %   k=k+1;
  %   dt(k,1)=tevt(j+1,1)-tevt(j,1);
  %   dloc(k,:)=implocplt(j+1,1:2)-implocplt(j,1:2);
  %   dprojloc(k,:)=impprojlocplt(j+1,1:2)-impprojlocplt(j,1:2);
  % end
end