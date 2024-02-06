function [arearat,area,normarea] = area_of_overlap_circs_grid(cnt,rad,gridsize,pltflag)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % [arearat,area,normarea] = area_of_overlap_circs_grid(cnt,rad,nrand,nrepeat,pltflag)
  % 
  % This is the function to roughly estimate the total covered area by 
  % multiple (possibly) overlapping circles. Each circles are defined by the
  % center and radius. Radi are the same for all circles if the input is a 
  % single value instead of the same size as the center. These N circles have
  % a nominal area of N*pi*r^2 in total. The function would also return the 
  % ratio between 2 areas, denoting the extent of their overlapping. 
  % -- the algorithm is to generate a very fine grid, determined by users, 
  % then ask which grid points are inside ANY circle, unique then count. 
  % --note that there is a lower limit for the total covered area, when all
  % circles are perfectly overlapping. When they have the same size, the 
  % limit is the area of one circle. Therefore, the area ratio has a minimum
  % of 1/ncirc.
  %
  %
  %
  % Chao Song, chaosong@princeton.edu
  % First created date:   2024/01/04
  % Last modified date:   2024/01/04
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  defval('gridsize',0.01);
  defval('pltflag',0);
  
  ncirc = size(cnt,1);
  if length(rad) == 1
    rad = rad*ones(ncirc,1);
  end
  
  %generate circles
  circ=cell(ncirc,1);
  for i = 1: ncirc
    radi=rad(i);
    circi=[];
    [circi(:,1),circi(:,2)] = circle_chao(cnt(i,1),cnt(i,2),radi,0.1);
    circ{i}=circi;
  end

  %min rectangle that outlines circles
  tmp=cat(1,circ{:});  %concatenate all circles' coordinates
  xlim=minmax(tmp(:,1)');  %xlim and ylim define a rectangle
  ylim=minmax(tmp(:,2)');
  xran=range(tmp(:,1)); %length and width of rectangle
  yran=range(tmp(:,2));
  recarea=xran*yran;

  %generate grid points inside the rectangle
  xvec=xlim(1): gridsize: xlim(2);
  yvec=ylim(1): gridsize: ylim(2);
  [xmesh,ymesh] = meshgrid(xvec,yvec);
  gridpts = [xmesh(:) ymesh(:)];

  %which and how many points inside ANY circle
  tmp1 = [];
  for j = 1: ncirc
    bnd = circ{j};
    [iin,ion] = inpolygon(gridpts(:,1),gridpts(:,2),bnd(:,1),bnd(:,2));
    isinbnd = iin | ion;
    indj = find(isinbnd==1);
    tmp1 = [tmp1; reshape(indj,[],1)];  
  end
  %use counts to approximate area
  normarea=length(tmp1);  %norminal area, sum of all circles  
  tmp2=unique(tmp1);  
  area=length(tmp2);
  arearat=area/normarea;
  % if arearat > 1
  %   arearat=1;
  % end
  % overlaprat=1-arearat;    
  
  if pltflag
    f=initfig;
    ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
    for i = 1: ncirc
      circi=circ{i};
      plot(ax,circi(:,1),circi(:,2),'k-');
    end
    scatter(ax,gridpts(tmp2,1),gridpts(tmp2,2),2,'r','filled');   
    text(ax,0.98,0.05,sprintf('tot area %d',area),'Units','normalized',...
      'HorizontalAlignment','right');
    text(ax,0.98,0.1,sprintf('norm area %d',normarea),...
      'Units','normalized','HorizontalAlignment','right');
    text(ax,0.98,0.15,sprintf('arearat %.2f',arearat),...
      'Units','normalized','HorizontalAlignment','right');
    % text(ax,0.98,0.2,sprintf('overlaprat %.2f',overlaprat),...
    %   'Units','normalized','HorizontalAlignment','right');  
    axis equal
  end 
  
  % keyboard