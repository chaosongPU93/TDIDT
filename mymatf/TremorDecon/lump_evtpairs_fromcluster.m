function [imppairst,imppaired,mampcont,amppair,dt,dloc,dloc_spl,imppairuevt]=...
  lump_evtpairs_fromcluster(catclus,n,lumpst,lumped,sps,ftrans)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [catclus,catimp,catuimp,catind]=...
%   evtcluster_ex(nbst,imp,nsrc,mmax,sps,timetype)
%
% Given clusters that are categorized, lump all unique event pairs of interest
% into one from all types of clusters ('catclus'). For example, you want to lump all 
% consecutive event pairs ('n'=1) from doublets, triplets, and so on, so 'lumpst'
% is the same as 'n' and 'lumped' is the size of 'catclus', which is also 'mmax'
% 
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/04/25
% Last modified date:   2024/04/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('lumpst',n);
defval('lumped',size(catclus,1));
defval('sps',160);
defval('ftrans','interpchao');

mmax = size(catclus,1); %largest cluster size

catimppairst = cell(mmax,1);  %for a certain category, start event of the pair
catimppaired = cell(mmax,1);  %for a certain category, end event of the pair
catmampcont = cell(mmax,1);  %for a certain category, median amp of ALL evts between event pair
catdt=cell(mmax,1); %for a certain category, diff time between each event pair
catdloc=cell(mmax,1); %for a certain category, diff loc between each event pair
catdlocspl=cell(mmax,1); %for a certain category, diff loc in spls between each event pair

%%%'m' decides what type of cluster start to look at, eg., if look at all N &
%%%N-1 in doublets, triplets and above, then n=1, and m starts from 1 to mmax
for m=n:mmax
  catclususe = catclus{m};
  nclus = size(catclususe,1); %num of clusters, consecu. clusters may share events!
  imppairst = []; %start event of the pair for ALL clusters, but same category
  imppaired = []; %end event of the pair for ALL clusters, but same category
  mampcont = []; %median amp of all continous evts between the start and end of the pair for ALL clusters, but same category
  for i = 1: nclus  %loop over all clusters
    impi = catclususe{i};
    impipairst = impi(1:end-n, :);  %list of starting src of a PAIR
    impipaired = impi(1+n:end, :);  %list of starting src of a PAIR
    imppairst = [imppairst; impipairst];
    imppaired = [imppaired; impipaired];    
    %median amp of all continous events between 'imppairst' and 'imppaired'
    mampicont = zeros(size(impipairst,1), 1);
    for j = 1: size(impipairst,1)
      impcont = impi(j: j+n,:);
      mampicont(j,1) = median(mean(impcont(:,[2 4 6]),2));
    end
    mampcont = [mampcont; mampicont];
  end
  if size(imppairst,1) ~= nclus*(m+1-n)
    error('number does not match');
  end
  [~,iunist] = unique(imppairst,'rows','stable');
  idupst = setdiff((1:size(imppairst,1))', iunist); %indices for duplicates of start event
  [~,iunied] = unique(imppaired,'rows','stable');
  iduped = setdiff((1:size(imppaired,1))', iunied); %indices for duplicates of end event
  %%%if a index shows in both 'idupst' and 'iduped', then the src PAIR is a duplicate
  idup = intersect(idupst,iduped);
  %%%make sure for certain category of cluster, no duplicate 
  imppairst(idup,:) = []; %remove the duplicates from both lists
  imppaired(idup,:) = [];
  mampcont(idup) = [];
  
  %store src pairs for a certain category of cluster
  catimppairst{m} = imppairst;
  catimppaired{m} = imppaired;
  catmampcont{m} = mampcont;
  
  %for found UNIQUE event PAIRS
  if ~isempty(imppaired)
    dt = imppaired(:,1) - imppairst(:,1);
    dloc_spl = imppaired(:,7:8) - imppairst(:,7:8);
    %for diff in loc in map view, do the mapping first, then diff
    implocpairst = off2space002(imppairst(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
    implocpaired = off2space002(imppaired(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
    dloc = implocpaired(:,1:2) - implocpairst(:,1:2);
  else
    dt=[];
    dloc=[];
    dloc_spl=[];
  end
  catdt{m}=dt;
  catdloc{m}=dloc;
  catdlocspl{m}=dloc_spl;  
end
% keyboard

%%%'lumpst' determine from which category of clusters to lump together, must be >= n
% lumpst = n; %e.g., n=1,lumpst=n, means lump all N and N-1 from doublets and above
% lumped = mmax;
imppairst=cat(1,catimppairst{lumpst:lumped});
imppaired=cat(1,catimppaired{lumpst:lumped});
mampcont=cat(1,catmampcont{lumpst:lumped});

%%%in case of duplicate event PAIRS between categories, find them and remove
[~,iunist] = unique(imppairst,'rows','stable');
idupst = setdiff((1:size(imppairst,1))', iunist);
[~,iunied] = unique(imppaired,'rows','stable');
iduped = setdiff((1:size(imppaired,1))', iunied);
%if a index shows in both 'idupst' and 'iduped', then the src PAIR is a duplicate
idup = intersect(idupst,iduped);
%make sure for certain category of cluster, no duplicate
imppairst(idup,:) = []; %remove the duplicates from both lists
imppaired(idup,:) = [];
mampcont(idup,:) = [];

%mean amp of the pair
amppair = mean([mean(imppairst(:,[2 4 6]),2) mean(imppaired(:,[2 4 6]),2)], 2); 

%obtain the diff location, etc. of the LUMPPED src PAIRS
dt = imppaired(:,1) - imppairst(:,1);
dloc_spl = imppaired(:,7:8) - imppairst(:,7:8);
implocpairst = off2space002(imppairst(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
implocpaired = off2space002(imppaired(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
dloc = implocpaired(:,1:2) - implocpairst(:,1:2);

%%%unlike in 'lfedloc_incluster', essentially, we want to know the UNIQUE EVENTS
%%%associated with these event pairs
imppairevt = [imppairst; imppaired];
imppairuevt = unique(imppairevt,'rows','stable');








