function [impcluster,impuni,ncluster,impclusterk,impunii]=evtcluster(nbst,imp,nsrc,m,dtcut,sps,timetype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [impcluster,impuni,ncluster,impclusterk,impunii]=evtcluster(nbst,imp,nsrc,m,dtcut,sps,timetype)
%
% This function is to obtain the events in the cluster that is defined by
% the differential time of the first (N-m) and last (N) event in the cluster 
% is smaller than 'dtcut' in sec.  All events in the eligible cluster are 
% returned, including duplicates. The unique events from the clusters lumped
% for each burst are also included.  
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/01/23
% Last modified date:   2024/01/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,dtplt,indimpdtcut] = med_amp_incluster(nbst,imp,nsrc,m,dtcut,sps,timetype);
ncluster=size(indimpdtcut,1); %num of clusters
impclusterk=cell(ncluster,1); %events in each cluster 
impunii=cell(nbst,1); %unique events in clusters for each burst
k = 0;  %count of cluster
for i = 1: nbst
%   i=181;
  ist = sum(nsrc(1:i-1))+1;
  ied = ist+nsrc(i)-1;
  impi = imp(ist:ied,:);
  % tevt = impi(:,1);

  %which type of time to use
  if strcmp(timetype,'tori')
    [~,impi]=tarvl2tori(impi,sps,ftrans,1);  %return the origin time 
  end

  indibst = find(dtplt(:,2)==i);  %ind of starting evt of the eligible cluster   
  if ~isempty(indibst)
    idx2 = cat(1,indimpdtcut{indibst}); %cat all ind from the same burst
    idx3 = unique(idx2);  %only save the unique ones
    impunii{i} = impi(idx3,:); %recover the events 
    for j = 1: length(indibst)
      k=k+1;
      impclusterk{k} = impi(indimpdtcut{indibst(j)},:);
    end
  end
end
impuni = cat(1,impunii{:}); %unique events included in clusters
impcluster = cat(1,impclusterk{:}); %events included in clusters, with duplicates
