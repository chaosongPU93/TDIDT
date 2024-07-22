function [fracuimp,nuimp,fracuimpsum]=frac_uniqevt_incluster2(catuimp,catclus,imp,nsrc,mmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the fraction of the catalog in terms of unique events in different
% type of clusters. Categorization of clusters have been done and is taken
% as the input 'catclus'. 'catuimp' essentially stores the unique events in 
% different clusters. 
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/04/04
% Last modified date:   2024/04/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuimp = zeros(mmax,1);
for m=1: mmax
  nuimp(m) = size(catuimp{m},1);
  dtcut = 0.25*m+0.125;
  % fprintf('%d clusters of %d consecutive events (%d/%d unique) w/i %.3f s \n',...
  %   size(catclus{m},1), m+1, nuimp(m), sum(nsrc), dtcut);
end
niso = length(imp)-sum(nuimp);
fraciso = niso/length(imp)*100;
nuimp = [niso; nuimp];
fracuimp = nuimp/length(imp)*100;

for i = 1:length(fracuimp)
  fracuimpsum(i) = sum(fracuimp(i:end));
end
