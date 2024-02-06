function [propang,locxyprop,fitobj,gof,output] = ...
  propadirection(timevec,locxy,robustflg,angleopt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [propang,locxyprop,fitobj,gof,output] = propadirection(timevec,locxy,angleopt)
%
% This function is to first compute the propagation direction events that
% are represented by 'timevec' and 'locxy'. Locations are projected along
% a series of trial directions. A Bisquare robust linear fit is implemented
% to the projected locations. The direction along which, either the RMSE is
% the minimum, or the slope is the maximum, is chosen as the propagation
% direction. Option is made by 'angleopt'. Projected locations and the fit
% statistics along the propagation direction is also returned.
% See 'srcprojdistNtoNm.m' as well.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/01/18
% Last modified date:   2024/01/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('robustflg',1);
defval('angleopt','minrmse');

%%%number of sources has to be larger than 2, otherwise there can be only 1 possible line to fit
if size(locxy,1) <= 2
  propang=[];
  locxyprop = [];
  fitobj = [];
  gof=[];
  output=[];
  return
else
  
  dangle = 5; %increment of trial angles
  angle = 0: dangle: 360-dangle;
  slope = zeros(length(angle),1);
  rmse = zeros(length(angle),1);
  fttpfree = fittype( @(a,b,x) a*x+b);
  
  %%% find best angle, now is mainly to get the variation of se and slope with the trial angle
  for iang = 1: length(angle)
    %%% propagation trial
    locxyprop = locxy;
    for jj = 1: size(locxy,1)
      x0 = locxyprop(jj,1);
      y0 = locxyprop(jj,2);
      [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),[0 0]);
      locxyprop(jj,1) = newx;
      locxyprop(jj,2) = newy;
    end
    % linear robust least square
    if robustflg
      [fitobj,gof,~] = fit(timevec, locxyprop(:,1),fttpfree,'Robust','Bisquare',...
        'StartPoint',[1 1]);
    else
      [fitobj,gof,~] = fit(timevec, locxyprop(:,1),fttpfree,'StartPoint',[1 1]);
    end
    % output fit parameters
    coef = coeffvalues(fitobj);
    slope(iang) = coef(1);
    rmse(iang) = gof.rmse;
  end
  
  %%% best angle estimate
  ind = find(slope>0);
  ind3 = find(rmse(ind)==min(rmse(ind)));     % rmse, Root Mean Squared Error, the smaller, the better
  if length(ind3) > 1
    disp(strcat({'multiple angles have same rmse: '},num2str(angle(ind3))));
  end
  angrmse = angle(ind(ind3(1)));
  ind6 = find(slope==max(slope)); % one with the largest slope, i.e., migrating speed
  if length(ind6) > 1
    disp(strcat({'multiple angles have same slope:'},num2str(angle(ind6))));
  end
  angslope = angle(ind6(1));
  
  if strcmp(angleopt,'minrmse')
    propang = angrmse;
  elseif strcmp(angleopt,'maxslope')
    propang = angslope;
  end
  
  %%%use the min-rmse direction as the projection direction
  locxyprop = locxy;
  for jj = 1: size(locxy,1)
    x0 = locxyprop(jj,1);
    y0 = locxyprop(jj,2);
    [newx,newy] = coordinate_rot(x0,y0,-(propang-90),[0 0]);
    locxyprop(jj,1) = newx;
    locxyprop(jj,2) = newy;
  end
  % linear robust least square
  if robustflg
    [fitobj,gof,output] = fit(timevec, locxyprop(:,1),fttpfree,'Robust','Bisquare',...
      'StartPoint',[1 1]);
  else
    [fitobj,gof,output] = fit(timevec, locxyprop(:,1),fttpfree,'StartPoint',[1 1]);
  end
  
end
