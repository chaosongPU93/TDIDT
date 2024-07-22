function [distnn1,dist2all,pca2vec,projang1,dprojxy1nn1,dprojxy12all,...
  projang2,dprojxy2nn1,dprojxy22all]=dist_evtcustom(impiwin,implociwin,sps,ftrans,timetype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [mdistnn1,mdist2all,pca2vec,projang1,mprojx1nn1,mprojx12all,projang2,mprojx2nn1,mprojx22all]=...
%   dist_evtcustom(impiwin,implociwin,sps,ftrans,timetype)
%
% Similar to 'dist_evtcluster', this function is able to compute the various
% type of distance between events inside a customized set of impulses
% 'impiwin'. To define 'consecutive' events, an event time, either the
% arrival time, or relative origin time needs to be assigned by 'timetype'.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/01/23
% Last modified date:   2024/01/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('timetype','tori');

if size(impiwin,1) ==1
  distnn1=[];
  dist2all=[];
  pca2vec=[];
  projang1=[];
  dprojxy1nn1=[];
  dprojxy12all=[];
  projang2=[];
  dprojxy2nn1=[];
  dprojxy22all=[];
  return
else
  
  %which type of time to use
  if strcmp(timetype,'tori')
    %return the origin time, sorted srcs, and indices
    [tevt,~,indsort]=tarvl2tori(impiwin,sps,ftrans,1);
    implociwin = implociwin(indsort, :);
  else
    tevt=impiwin(:,1);
  end
  
  %abs distance between consecutive events
  [~,~,eucdist] = srcdistNtoNm(tevt,implociwin(:,1:2),1);
  distnn1=eucdist{1};
  % mdistnn1 = median(distnn1);
  
  %abs distance from each to all others
  [~,~,dist2all] = srcdistall(tevt,implociwin(:,1:2));
  % mdist2all = median(dist2all);
  
  if size(impiwin,1) > 2
    %Principal component analysis
    [coeff,~,angle,anglegeo]=pcaellipse(implociwin(:,1:2));
    pca2vec=[coeff(1,2) coeff(2,2)];
    
    %if using the 2nd PC
    projang1 = angle(2);
    
    %distance between consecutive events along the 2nd pca direction
    [~,~,projxy] = customprojection(implociwin(:,1:2),projang1);
    dprojxy1nn1 = diffcustom(projxy,1,'forward');
    % mdprojx1nn1 = median(abs(dprojxy1nn1(:,1)));
    
    %distance from each to all others along the pca direction
    [~,dprojxy12all] = srcdistall(tevt,projxy);
    % mdprojx12all = median(abs(dprojxy12all(:,1)));
    
    
    %if using the propagation direction, can be rough though
    projang2 = propadirection(tevt/sps,implociwin(:,1:2),0);
    if isempty(projang2) && size(implociwin,1)==2  %if there happens to be 2 points
      %atan2d(dY,dX)
      dum=atan2d(implociwin(2,2)-implociwin(1,2), implociwin(2,1)-implociwin(1,1));
      [anggeo,anggeooppo]=angatan2d2geo(dum);
      projang2=min([anggeo,anggeooppo]);  %clockwise from N
    end
    
    %distance between consecutive events along the propagation direction
    [~,~,projxy] = customprojection(implociwin(:,1:2),projang2);
    dprojxy2nn1 = diffcustom(projxy,1,'forward');
    % mdprojx2nn1 = median(abs(dprojxy2nn1(:,1)));
    
    %distance from each to all others along the propagation direction
    [~,dprojxy22all] = srcdistall(tevt,projxy);
    % mdprojx22all = median(abs(dprojxy22all(:,1)));
    
  else
    pca2vec=[];
    projang1=[];
    dprojxy1nn1=[];
    dprojxy12all=[];
    projang2=[];
    dprojxy2nn1=[];
    dprojxy22all=[];
  end
end
