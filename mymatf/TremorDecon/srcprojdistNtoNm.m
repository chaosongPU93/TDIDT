function [locxyproj,dlocxyproj,stats] = srcprojdistNtoNm(timevec,locxy,m,sps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [locxyproj,dlocxyproj,stats] = srcprojdistNtoNm(timevec,locxy,m,sps)
%
% Compute the projected distance between the sources N and N-1 until N and
% N-m, in the sequential order of origin time or arrival time,
% along the proposed propagation direction (direction that has the min RMSE
% of the robust linear regression). 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/03/08
% Last modified date:   2023/03/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('m',1);
defval('sps',160);

%%%number of sources has to be larger than 2, otherwise there can be only 1 possible line to fit 
if size(locxy,1) <= 2
  locxyproj = [];
  dlocxyproj = [];
  stats=[];
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
    locxyproj = locxy;
    for jj = 1: size(locxy,1)
      x0 = locxyproj(jj,1);
      y0 = locxyproj(jj,2);
      [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),[0 0]);
      locxyproj(jj,1) = newx;
      locxyproj(jj,2) = newy;
    end
    % linear robust least square
    [fitobj,gof,~] = fit(timevec/sps, locxyproj(:,1),fttpfree,'Robust','Bisquare',...
      'StartPoint',[1 1]);
    % output fit parameters
    coef = coeffvalues(fitobj);
    slope(iang) = coef(1);
    rmse(iang) = gof.rmse;
  end
  
  %%% best angle estimate from hf
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
  
  angproj = angrmse;
  
  %%%use the min-rmse direction as the projection direction
  locxyproj = locxy;
  for jj = 1: size(locxy,1)
    x0 = locxyproj(jj,1);
    y0 = locxyproj(jj,2);
    [newx,newy] = coordinate_rot(x0,y0,-(angproj-90),[0 0]);
    locxyproj(jj,1) = newx;
    locxyproj(jj,2) = newy;
  end
  % linear robust least square
  [fitobjproj,gof,output] = fit(timevec/sps, locxyproj(:,1),fttpfree,'Robust','Bisquare',...
    'StartPoint',[1 1]);
  
  %%%Some statistics
  x = timevec/sps;
  y = locxyproj(:,1);
  stats = statsofrobustlnfit(fitobjproj,gof,output,x,y);
  %add more statistics into the structure
  stats.angrmse = angrmse;
  stats.angslope = angslope;
  stats.angproj = angproj;
  
  dlocxyproj = cell(m, 1); %differential location along projected direction and its orthogonal
  for i = 1: m
    %%%between Nth and (N-i)th source
    dlocxyproj{i} = [diffcustom(locxyproj(:,1), i,'forward') diffcustom(locxyproj(:,2), i,'forward')];
  end
  
end


