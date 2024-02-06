function [mdistnn1,mdist2all,pca2vec,projang1,mprojx1nn1,mprojx12all,...
  projang2,mprojx2nn1,mprojx22all,distnn1cat,dist2allcat,...
  dprojxy1nn1cat,dprojxy12allcat,dprojxy2nn1cat,dprojxy22allcat]=...
  dist_evtcluster(impcluster,implocclus,sps,ftrans,m,timetype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [mdistnn1,mdist2all,projang1,mprojx1nn1,mprojx12all,projang2,mprojx2nn1,mprojx22all]=...
%   dist_evtcluster(impcluster,implocclus,sps,ftrans,m,timetype)
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

ncluster = round(size(impcluster,1)/(m+1));

mdistnn1 = zeros(ncluster,1);
mdist2all = zeros(ncluster,1);
pca2vec = zeros(ncluster,2);
projang1 = zeros(ncluster,1);
mprojx1nn1 = zeros(ncluster,1);
mprojx12all = zeros(ncluster,1);
projang2 = zeros(ncluster,1);
mprojx2nn1 = zeros(ncluster,1);
mprojx22all = zeros(ncluster,1);

distnn1cat=[];
dist2allcat=[];
dprojxy1nn1cat=[];
dprojxy12allcat=[];
dprojxy2nn1cat=[];
dprojxy22allcat=[];

for i = 1: ncluster
  ist=(i-1)*(m+1)+1;
  ied=i*(m+1);
  implocplt=implocclus(ist:ied,:);
  impplt=impcluster(ist:ied,:);
  tevt=impplt(:,1);

  %which type of time to use
  if strcmp(timetype,'tori')
    [torisplst,~,indsort]=tarvl2tori(impplt,sps,ftrans,1);  %return the origin time 
    implocplt = implocplt(indsort, :);
    tevt = torisplst;
  end

  %Principal component analysis
  [coeff,score,latent,tsquared,explained] = pca(implocplt(:,1:2));
  %each column in coeff represent the unit vector of each principal component
%   x0=mean(implocplt(:,1));
%   y0=mean(implocplt(:,2));
%   angle1=atan2d(coeff(2,1),coeff(1,1)); %angle of principal axis 1
%   [anggeo,anggeooppo]=angatan2d2geo(angle1);
%   ang1=min([anggeo,anggeooppo]);  %clockwise from N
  % plot(ax,x0+[0,coeff(1,1)],y0+[0,coeff(2,1)],'k-','linew',1);
  angle2=atan2d(coeff(2,2),coeff(1,2)); %angle of principal axis 2
  pca2vec(i,:)=[coeff(1,2) coeff(2,2)];
  [anggeo,anggeooppo]=angatan2d2geo(angle2);
  ang2=min([anggeo,anggeooppo]);  %clockwise from N
%   %Compute 95% CI using percentile method, in case not Gaussian
%   %if Gaussian, use 2*std~95% is fine
%   p=95; % CI level
%   CIx=prctile(score(:,1), [(100-p)/2; p+(100-p)/2]); % x CI [left, right]
%   CIy=prctile(score(:,2), [(100-p)/2; p+(100-p)/2]); % y CI [lower, upper]
%   semia=range(CIx)/2;
%   semib=range(CIy)/2;
%   angrot = angle1;
%   %this would create an ellipse roughly circumventing 95% of data
%   [x, y] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);
  % plot(ax,x,y,'k-');
  % keyboard

  %abs distance between consecutive events
  [~,~,eucdist] = srcdistNtoNm(tevt,implocplt,1);
  distnn1=eucdist{1};
  mdistnn1(i,1) = median(distnn1);
  distnn1cat=[distnn1cat; distnn1];

  %abs distance from each to all others
  [~,~,dist2all] = srcdistall(tevt,implocplt);
  mdist2all(i,1) = median(dist2all);
  dist2allcat=[dist2allcat; dist2all];

  %if using the 2nd PC  
  projang1(i,1) = ang2; 

  %distance between consecutive events along the pca direction
  [~,~,projxy] = customprojection(implocplt,projang1(i,1));
  dprojxy1nn1 = diffcustom(projxy,1,'forward');
  mdprojx1nn1 = median(abs(dprojxy1nn1(:,1)));
  dprojxy1nn1cat = [dprojxy1nn1cat; dprojxy1nn1];

  %distance from each to all others along the pca direction
  [~,dprojxy12all] = srcdistall(tevt,projxy);
  mprojx12all(i,1) = median(abs(dprojxy12all(:,1)));
  dprojxy12allcat = [dprojxy12allcat; dprojxy12all];

  %if using the propagation direction, can be rough though
  projang2(i,1) = propadirection(tevt/sps,implocplt,0);

  %distance between consecutive events along the propagation direction
  [~,~,projxy] = customprojection(implocplt,projang2(i,1));
  dprojxy2nn1 = diffcustom(projxy,1,'forward');
  mdprojx2nn1 = median(abs(dprojxy2nn1(:,1)));
  dprojxy2nn1cat = [dprojxy2nn1cat; dprojxy2nn1];

  %distance from each to all others along the propagation direction
  [~,dprojxy22all] = srcdistall(tevt,projxy);
  mprojx22all(i,1) = median(abs(dprojxy22all(:,1)));
  dprojxy22allcat = [dprojxy22allcat; dprojxy22all];
end
