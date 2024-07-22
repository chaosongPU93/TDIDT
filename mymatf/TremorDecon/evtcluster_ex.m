function [catclus,catclusbst,catimp,catuimp,catmedamp,catdtnnm,catind]=...
  evtcluster_ex(nbst,imp,nsrc,mmax,sps,timetype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [catclus,catimp,catuimp,catind]=...
%   evtcluster_ex(nbst,imp,nsrc,mmax,sps,timetype)
%
% Function 'evtcluster' tries to find clusters of events whose first and last
% event is separated by less than 'dcut'. However, the events included in that 
% cluster is NOT exclusive. For example, two events are included in doublets, 
% they can also appear in triplets. It is fine when you want ALL doublets, but 
% it contains duplicates (events or event pairs) if you ONLY want you doublets 
% are NOT part of any triplets or quadupltets. This function is trying to 
% satisfy that requirement -- to obtain mutually exclusive catagories of
% doublets, triplets, etc. 
% --To do that, start from the mmax, the max m where the cluster of m+1
% events within 'dcut' exists. 
% --From mmax, first get the cluster of mmax+1 events, identify them, put all 
% events aside to the category.
% --Remove events from 2 to mmax from the main catalog, in the remaining cat,
% find the next clusters of mmax events.
% --Continue till doublets.
% --So essentially you always start from the mmax, instead of 1. And this
% function would classify all EXCLUSIVE categories at one time. 
% --Note that the final event pairs from these EXCLUSIVE categories should be
% unique, but since you are preserving the firts and last events in each round
% just in case you would miss a cluster in the next round, the events in each
% cluster might be double-counted across different clusters. This might affect
% when you try to make a statement like 'the fraction of the catalog that only
% appears as doublets', in which case you will exclude events from doublet
% category that are also shown in a longer cluster!!
%
% --RETURN:
%     catclus = cell(mmax,1); %each category, all srcs included in EACH cluster  
%     catclusbst = cell(mmax,1); %each category, all srcs included in EACH cluster, separated by burst 
%     catimp = cell(mmax,1);  %each category, all srcs included in ALL clusters lumped 
%     catuimp = cell(mmax,1); %each category, UNIQUE srcs included in ALL clusters lumped 
%     catmedamp = cell(mmax,1); %each category, the median amp of srcs included in EACH cluster
%     catdtnnm = cell(mmax,1); %each category,  the diff time of between N and N-m of EACH cluster
%     catind = cell(mmax,1);  %each category, each burst, indices of srcs included in EACH cluster
%
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/03/15
% Last modified date:   2024/03/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('sps',160);
defval('timetype','tarvl');

%store different categories of clusters, from 1 to mmax
catclus = cell(mmax,1); %each category, all srcs included in EACH cluster  
catclusbst = cell(mmax,1); %each category, all srcs included in EACH cluster, separated by burst 
catimp = cell(mmax,1);  %each category, all srcs included in ALL clusters lumped 
catuimp = cell(mmax,1); %each category, UNIQUE srcs included in ALL clusters lumped 
catind = cell(mmax,1);  %each category, each burst, indices of srcs included in EACH cluster
catmedamp = cell(mmax,1); %each category, the median amp of srcs included in EACH cluster
catdtnnm = cell(mmax,1); %each category,  the diff time of between N and N-m of EACH cluster

%CURRENT catalog where srcs will be kept being removed
impcur = imp;
nsrccur = nsrc;

for m = mmax:-1:1
%   m = mmax;
  dtcut = 0.25*m+0.125;
  
  %find eligible clusters, return their indices in terms of the CURRENT catalog
  [medamp,dtnnm,indimpdtcut,inddtcutbst] = med_amp_incluster(nbst,impcur,nsrccur,m,dtcut,sps,timetype);
  
  %note the indices of srcs in the CURRENT catalog included in the clusters for each burst
  catind{m} = inddtcutbst; %if you keep removing srcs from each burst, then src ind is invalid 
  
  %median amplitude of the cluster
  catmedamp{m} = medamp;
  
  %diff time between N and N-m of the cluster
  catdtnnm{m} = dtnnm;
  
  ncluster=size(indimpdtcut,1); %num of clusters, in total
  impclusteribst=cell(nbst,1); %events in each cluster
  impclusterkclus=cell(ncluster,1); %events in each cluster
  impunii=cell(nbst,1); %unique events in clusters for each burst
  k = 0;  %count of clusters in total (all bursts combined)
  
  impdum = [];  %dummy catalog to keep removing srcs from it
  nsrcdum = zeros(nbst,1);  %num of srcs in each burst of the updated dummy catalog
  for i = 1: nbst
%     i
    %for the CURRENT catalog, identify srcs in one burst 
    ist = sum(nsrccur(1:i-1))+1;
    ied = ist+nsrccur(i)-1;
    impi = impcur(ist:ied,:);
    impidum = impi;
    % tevt = impi(:,1);
    
    %which type of time to use
    if strcmp(timetype,'tori')
      [~,impi]=tarvl2tori(impi,sps,ftrans,1);  %return sorted src by origin time
    end
    
    %correspond src indices to actual impulse information
    inddtcuti = inddtcutbst{i};
    if ~isempty(inddtcuti)
      idx2 = cat(1,inddtcuti{:}); %cat all ind from the same burst
      idx3 = unique(idx2);  %indices of unique srcs
      %for smaller m, there can be unique clusters, but the first or last events
      %can be shared by other clusters in the same burst
      impunii{i} = impi(idx3,:); %recover the events
      nclusteri = size(inddtcuti,1); %num of clusters in burst i
      impclusterj = cell(nclusteri,1);
      indrem = [];
      for j = 1: nclusteri
        k=k+1;
        ind = inddtcuti{j};
        impclusterkclus{k} = impi(ind,:);
        impclusterj{j} = impi(ind,:);
%         locator(k,:) = [i j ind(1) ind(end)];
%         locator{m} = aa;
        
        %removing all events in the middle except the first and last, eg. 3 from
        %1,2,3
        indrem = [indrem; reshape(ind(2:end-1),[],1)];
      end
      impclusteribst{i} = impclusterj;
      
      impidum(indrem, :) = [];
    end
    
    nsrcdum(i) = size(impidum,1);
    impdum = [impdum; impidum];
    
  end
  
  impuni = cat(1,impunii{:}); %unique events included in clusters
  catuimp{m} = impuni;
  catclusbst{m} = impclusteribst;
  catclus{m} = impclusterkclus;
  catimp{m} = cat(1,impclusterkclus{:});
  
%   %check if their are duplicates  
%   aa = categoryind{m};
%   bb = find(~cellfun(@isempty,aa));
%   for ichk = 1: length(bb)
%     ibst = bb(ichk); 
%     ist = sum(nsrcdum(1:ibst-1))+1;
%     ied = ist+nsrcdum(ibst)-1;
%     impdumchk = impdum(ist:ied,:);
% 
%     ist = sum(nsrcuse(1:ibst-1))+1;
%     ied = ist+nsrcuse(ibst)-1;
%     impusechk = impuse(ist:ied,:);
%     
%     cc=aa{ibst};
%     for jchk = 1: size(cc,1)
%       dd = cc{jchk};
%       dum = impusechk;
%       dum(dd(2:end-1),:) = [];
%       if ~isequaln(dum,impdumchk)
%         fprintf('Something wrong at m=%d, bst=%d (ichk=%d), %d th cluster',m,ibst,ichk,jchk);
%       end  
%     end
% 
%   end
    
  %refresh the remaining list
  nsrccur = nsrcdum;
  impcur = impdum;
  
end















