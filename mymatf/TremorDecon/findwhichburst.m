function ibst = findwhichburst(impi,imp,nsrc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ibst = findwhichburst(impi,imp,nsrc)
%
% Imagine you have a few particular source impulses 'impi', and you want to
% know which bursts they belong to, given you have all sources 'imp' and the
% num of sources for each burst 'nsrc'.
%
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2024/11/06
% Last modified date:   2024/11/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbst = length(nsrc);
for i = 1: nbst
  ist = sum(nsrc(1:i-1))+1;
  ied = ist+nsrc(i)-1;
  indrange(i,:) = [ist ied];
end   
  
[~, ind] = ismember(impi,imp,'rows');

ibst = findwhichrange(ind,indrange);