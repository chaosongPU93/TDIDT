function mmax=getmmaxcluster(nbst,imp,nsrc,sps,timetype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the maximum of m for a cluster of m+1 events to exit for a catalog 
% 'imp' and number of src for each burst 'nsrc'. Essentiall to increase 
% 'm' until there is no such cluster any more.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/04/04
% Last modified date:   2024/04/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmax=1;
dtcut = 0.25*mmax+0.125;
%it is fine to use this function to obtain the mmax, but not recommended to
%obtain the actual exclusive clusters!
[~,~,indimpdtcut] = med_amp_incluster(nbst,imp,nsrc,mmax,dtcut,sps,timetype);
while ~isempty(indimpdtcut)
  mmax=mmax+1;
  dtcut = 0.25*mmax+0.125;
  [~,~,indimpdtcut] = med_amp_incluster(nbst,imp,nsrc,mmax,dtcut,sps,timetype);
end
mmax=mmax-1;
