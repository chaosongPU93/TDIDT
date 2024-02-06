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

%which type of time to use
if strcmp(timetype,'tori')
  [tevt,~,indsort]=tarvl2tori(impiwin,sps,ftrans,1);  %return the origin time
  implociwin = implociwin(indsort, :);
else
  tevt=impiwin(:,1);
end

%Principal component analysis
[coeff,score,latent,tsquared,explained] = pca(implociwin(:,1:2));
%each column in coeff represent the unit vector of each principal component
%   x0=mean(implociwin(:,1));
%   y0=mean(implociwin(:,2));
%   angle1=atan2d(coeff(2,1),coeff(1,1)); %angle of principal axis 1
%   [anggeo,anggeooppo]=angatan2d2geo(angle1);
%   ang1=min([anggeo,anggeooppo]);  %clockwise from N
% plot(ax,x0+[0,coeff(1,1)],y0+[0,coeff(2,1)],'k-','linew',1);
angle2=atan2d(coeff(2,2),coeff(1,2)); %angle of principal axis 2
pca2vec=[coeff(1,2) coeff(2,2)];
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
[~,~,eucdist] = srcdistNtoNm(tevt,implociwin,1);
distnn1=eucdist{1};
% mdistnn1 = median(distnn1);

%abs distance from each to all others
[~,~,dist2all] = srcdistall(tevt,implociwin);
% mdist2all = median(dist2all);

%if using the 2nd PC
projang1 = ang2;

%distance between consecutive events along the pca direction
[~,~,projxy] = customprojection(implociwin,projang1);
dprojxy1nn1 = diffcustom(projxy,1,'forward');
% mdprojx1nn1 = median(abs(dprojxy1nn1(:,1)));

%distance from each to all others along the pca direction
[~,dprojxy12all] = srcdistall(tevt,projxy);
% mdprojx12all = median(abs(dprojxy12all(:,1)));

%if using the propagation direction, can be rough though
projang2 = propadirection(tevt/sps,implociwin,0);

%distance between consecutive events along the propagation direction
[~,~,projxy] = customprojection(implociwin,projang2);
dprojxy2nn1 = diffcustom(projxy,1,'forward');
% mdprojx2nn1 = median(abs(dprojxy2nn1(:,1)));

%distance from each to all others along the propagation direction
[~,dprojxy22all] = srcdistall(tevt,projxy);
% mdprojx22all = median(abs(dprojxy22all(:,1)));
