function [overlapratm,arearatm,aream] = area_of_overlap_circs_MC(cnt,rad,nrand,nrepeat,pltflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [area,overlaprat,arearat] = area_of_overlap_circs_MC(cnt,rad,nrand,nrepeat,pltflag)
% 
% This is the function to roughly estimate the total covered area by 
% multiple (possibly) overlapping circles. Each circles are defined by the
% center and radius. Radi are the same for all circles if the input is a 
% single value instead of the same size as the center. These N circles have
% a nominal area of N*pi*r^2 in total. The function would also return the 
% ratio between 2 areas, denoting the extent of their overlapping. 
% -- the algorithm is Monte-Carlo sampling. A large enough number points are
% uniformly randomly sampled in the rectange area. For each point, ask if it
% inside ANY circle. The ratio of points inside any circle with respect to 
% the total sampled is the same as the ratio between the total covered area
% of these circles and area of the rectangle.   
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/01/03
% Last modified date:   2024/01/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('nrand',1e4);
defval('nrepeat',10);
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
  normareai(i)=polyarea(circi(:,1),circi(:,2)); %approximation, NOT analytic form
  normareaai(i)=pi*radi^2;
end
normarea=sum(normareai);  %norminal area, sum of all circles
normareaa=sum(normareaai);  %analytic norminal area, sum of all circles

%generate 2D random points inside the a rectangle
tmp=cat(1,circ{:});  %concatenate all circles' coordinates
xlim=minmax(tmp(:,1)');  %xlim and ylim define a rectangle
ylim=minmax(tmp(:,2)');
xran=range(tmp(:,1)); %length and width of rectangle
yran=range(tmp(:,2));
recarea=xran*yran;

for ii = 1: nrepeat
  % rng('default');
  tmp = rand(nrand,2);  % uniformly distributed random numbers within [0,1];
  rnd(:,1) = xran.*tmp(:,1)+xlim(1);
  rnd(:,2) = yran.*tmp(:,2)+ylim(1);
  % k=0;  %count of points inside ANY circle
  tmp1=[];

  for j = 1: ncirc
    bnd = circ{j};
    [iin,ion] = inpolygon(rnd(:,1),rnd(:,2),bnd(:,1),bnd(:,2));
    isinbnd = iin | ion;
    indj = find(isinbnd==1);
    tmp1 = [tmp1; reshape(indj,[],1)];  
  end
  tmp1 = unique(tmp1);
  countin(ii)=length(tmp1);

  % for i = 1: nrand
  %   isinbnd=0;  %reset the flag which indicates if a point is inside or on a polygon
  %   for j = 1: ncirc
  %     bnd = circ{j};
  %     [iin,ion] = inpolygon(rnd(i,1),rnd(i,2),bnd(:,1),bnd(:,2));
  %     isinbnd = iin | ion;
  %     if isinbnd == 1 
  %       break  %continue to next rand point if a point is inside ANY circle,
  %     end
  %   end
  %   if isinbnd==1
  %     k=k+1;
  %     tmp1(k)=i;
  %   end
  % end
  % countin(ii)=k;

  tmp2=setdiff(1:nrand,tmp1);
  indin{ii}=tmp1;
  indout{ii}=tmp2;

  inrat(ii)=countin(ii)/nrand;  %ratio of points inside == ratio of area
  area(ii)=inrat(ii)*recarea;
  arearat(ii)=area(ii)/normarea;
  if arearat(ii) > 1
    arearat(ii)=1;
  end
  overlaprat(ii)=1-arearat(ii);

end
overlapratm = median(overlaprat);
arearatm = median(arearat);
aream = median(area);

if pltflag
  f=initfig;
  ax=f.ax(1); hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  for i = 1: ncirc
    circi=circ{i};
    plot(ax,circi(:,1),circi(:,2),'k-');
  end
  scatter(ax,rnd(indin{nrepeat},1),rnd(indin{nrepeat},2),5,'r','filled');
  scatter(ax,rnd(indout{nrepeat},1),rnd(indout{nrepeat},2),5,'c','filled');
  text(ax,0.98,0.05,sprintf('%d in',countin(nrepeat)),'Units','normalized',...
    'HorizontalAlignment','right');
  text(ax,0.98,0.1,sprintf('rec %.2f*%.2f=%.2f',xran,yran,recarea),'Units','normalized',...
    'HorizontalAlignment','right');
  text(ax,0.98,0.15,sprintf('norm area acc %.2f',normarea/normareaa),...
    'Units','normalized','HorizontalAlignment','right');
  text(ax,0.98,0.2,sprintf('norm area %.2f',normarea),...
    'Units','normalized','HorizontalAlignment','right');    
  text(ax,0.98,0.25,sprintf('tot area %.2f',aream),'Units','normalized',...
    'HorizontalAlignment','right');
  text(ax,0.98,0.3,sprintf('arearat %.2f',arearatm),...
    'Units','normalized','HorizontalAlignment','right');
  text(ax,0.98,0.35,sprintf('overlaprat %.2f',overlapratm),...
    'Units','normalized','HorizontalAlignment','right');  
  axis equal
end 

% keyboard