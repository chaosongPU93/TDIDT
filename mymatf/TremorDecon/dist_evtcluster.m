function [mdistnn1,mdist2all,pca2vec,projang1,mdprojx1nn1,mdprojx12all,...
  projang2,mdprojx2nn1,mdprojx22all,distnn1cat,dist2allcat,...
  dprojxy1nn1cat,dprojxy12allcat,dprojxy2nn1cat,dprojxy22allcat]=...
  dist_evtcluster(catclusm,sps,ftrans,timetype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [mdistnn1,mdist2all,pca2vec,projang1,mdprojx1nn1,mdprojx12all,...
%   projang2,mdprojx2nn1,mdprojx22all,distnn1cat,dist2allcat,...
%   dprojxy1nn1cat,dprojxy12allcat,dprojxy2nn1cat,dprojxy22allcat]=...
%   dist_evtcluster(catclusm,sps,ftrans,timetype)
%
% This function is to compute the various type of distance between events
% inside EACH cluster 'impcluster'. Types include the absolute distance
% between consecutive events, and from each to all others. They can be 
% absolute euclidean distance, or abs location difference projected along
% some direction. To define 'consecutive' events, an event time, either the
% arrival time, or relative origin time needs to be assigned by 'timetype'.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/01/23
% Last modified date:   2024/01/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('timetype','tori');

ncluster = size(catclusm,1); %num of clusters, consecu. clusters may share events!

mdistnn1 = zeros(ncluster,1);
mdist2all = zeros(ncluster,1);
pca2vec = zeros(ncluster,2);
projang1 = zeros(ncluster,1);
mdprojx1nn1 = zeros(ncluster,1);
mdprojx12all = zeros(ncluster,1);
projang2 = zeros(ncluster,1);
mdprojx2nn1 = zeros(ncluster,1);
mdprojx22all = zeros(ncluster,1);

distnn1cat=[];
dist2allcat=[];
dprojxy1nn1cat=[];
dprojxy12allcat=[];
dprojxy2nn1cat=[];
dprojxy22allcat=[];

for i = 1: ncluster
  imp = catclusm{i};
  tevt=imp(:,1);
  imploc = off2space002(imp(:,7:8),sps,ftrans,0); % 8 cols, format: dx,dy,lon,lat,dep,ttrvl,off12,off13
  %which type of time to use
  if strcmp(timetype,'tori')
    [torisplst,~,indsort]=tarvl2tori(imp,sps,ftrans,1);  %return the origin time 
    imploc = imploc(indsort, :);
    tevt = torisplst;
  end
  imploc = imploc(:,1:2);
  
  %Principal component analysis
  [coeff,~,angle,anglegeo]=pcaellipse(imploc);
  pca2vec(i,:)=[coeff(1,2) coeff(2,2)];

  %abs distance between consecutive events
  [~,~,eucdist] = srcdistNtoNm(tevt,imploc,1);
  distnn1=eucdist{1};
  mdistnn1(i,1) = median(distnn1);
  distnn1cat=[distnn1cat; distnn1];

  %abs distance from each to all others
  [~,~,dist2all] = srcdistall(tevt,imploc);
  mdist2all(i,1) = median(dist2all);
  dist2allcat=[dist2allcat; dist2all];

  %if using the 2nd PC  
  projang1(i,1) = anglegeo(2); 

  %distance between consecutive events along the pca direction
  [~,~,projxy] = customprojection(imploc,projang1(i,1));
  dprojxy1nn1 = diffcustom(projxy,1,'forward');
  mdprojx1nn1(i,1) = median(abs(dprojxy1nn1(:,1)));
  dprojxy1nn1cat = [dprojxy1nn1cat; dprojxy1nn1];

  %distance from each to all others along the pca direction
  [~,dprojxy12all] = srcdistall(tevt,projxy);
  mdprojx12all(i,1) = median(abs(dprojxy12all(:,1)));
  dprojxy12allcat = [dprojxy12allcat; dprojxy12all];

  %if using the propagation direction, can be rough though
  % projang2(i,1) = propadirection(tevt/sps,imploc,0);

  %if using the 1st PCA direction
  projang2(i,1) = anglegeo(1); 

  %if using the othoogonal direction to the domiant propagation, ie, the
  %elongation direction, NE
  projang2(i,1) = 45; 
  
  %distance between consecutive events along the projected direction
  [~,~,projxy] = customprojection(imploc,projang2(i,1));
  dprojxy2nn1 = diffcustom(projxy,1,'forward');
  mdprojx2nn1(i,1) = median(abs(dprojxy2nn1(:,1)));
  dprojxy2nn1cat = [dprojxy2nn1cat; dprojxy2nn1];

  %distance from each to all others along the projected direction
  [~,dprojxy22all] = srcdistall(tevt,projxy);
  mdprojx22all(i,1) = median(abs(dprojxy22all(:,1)));
  dprojxy22allcat = [dprojxy22allcat; dprojxy22all];
end




