function [coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dataxy,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [coeff,score,angle,anglegeo,x0,y0,semia,semib,ellx,elly]=pcaellipse(dataxy,p)
%
% This function is to carry out a Principal Component Analysis to 'dataxy'.
% Would return the unit vector and direction of each principal component, 
% and the projected location along that direction. According to the PCA,
% it also returns an ellipse that cover the data to the confidence level
% of 'p'.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/02/05
% Last modified date:   2024/02/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('p',95);  % CI level

%apply a PCA to the data to infer the structure
%Principal component analysis
dataxy = dataxy(:,1:2);
[coeff,score] = pca(dataxy);
%each column in coeff represent the unit vector of each principal component
x0=mean(dataxy(:,1));
y0=mean(dataxy(:,2));
angle1=atan2d(coeff(2,1),coeff(1,1)); %angle of principal axis 1
[anggeo,anggeooppo]=angatan2d2geo(angle1);
ang1=min([anggeo,anggeooppo]);  %clockwise from N
% plot(ax,x0+[0,coeff(1,1)],y0+[0,coeff(2,1)],'k-','linew',1);
if size(coeff) > 1
  angle2=atan2d(coeff(2,2),coeff(1,2)); %angle of principal axis 2
  [anggeo,anggeooppo]=angatan2d2geo(angle2);
  ang2=min([anggeo,anggeooppo]);  %clockwise from N
else
  angle2=[];
  ang2=[];
end
angle = [angle1 angle2]; %PCA direction, four-quadrant convention, ie, from E cclockwise to 180, cclockwise to -180
anglegeo = [ang1 ang2]; %PCA direction, geographic convention, ie, measure from N clockwise

if size(coeff) > 1
  %Compute 95% CI using percentile method, in case not Gaussian
  %if Gaussian, use 2*std~95% is fine
  % p=95; % CI level
  CIx=prctile(score(:,1), [(100-p)/2; p+(100-p)/2]); % x CI [left, right]
  CIy=prctile(score(:,2), [(100-p)/2; p+(100-p)/2]); % y CI [lower, upper]
  semia=range(CIx)/2;
  semib=range(CIy)/2;
  angrot = angle1;
  %this would create an ellipse roughly circumventing 95% of data
  [ellx, elly] = ellipse_chao(x0,y0,semia,semib,0.01,angrot,[x0,y0]);
  % keyboard
else
  semia=[]; semib=[]; ellx=[]; elly=[]; 
end


