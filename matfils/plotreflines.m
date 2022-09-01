function [h] = plotreflines(pin,vls,xy)

% plots dotted lines at the locations specified
%
% INPUT
% 
% pin      axis (or axes) to plot on
% vls      values at which to plot lines
% xy       whether the lines should be vertical 'x' 
%              horizontal 'y' (default)
%            BELOW: OFTEN DOESN'T WORK
%            or a number that is the slope of just one line to be plotted
%              on the axes type that exists
%            or two numbers that are another point on the line
%            if xy is a number or coordinate point, then vls should be a
%            point somewhere on the line
%
% OUTPUT
%
% h        handles to the lines plotted

defval('xy','y')
pin = pin(:);

if isfloat(xy) & length(xy)>1
  for k=1:length(xy)
    h(k,:) = plotreflines(pin,vls,xy(k));
  end
end

for n = 1:length(pin)
  p = pin(n);
  
  xlims = get(p,'xlim');
  ylims = get(p,'ylim');
  xlg = isequal(get(p,'xscale'),'log');
  ylg = isequal(get(p,'yscale'),'log');

  np = get(p,'nextplot');
  set(p,'nextplot','add')

  if isstr(xy)
    switch xy
     case 'y'
      h(:,n) = plot(p,repmat(xlims,[length(vls(:)),1]).',...
		    (vls(:)*[1 1]).');
     case 'x'
      h(:,n) = plot(p,(vls(:)*[1 1]).',...
		    (repmat(ylims,[length(vls(:)),1])).');
    end
  elseif length(xy)==1
    if xlg
      vls(1) = log(vls(1));
      xlims = log(xlims);
    end
    if ylg 
      vls(2) = log(vls(2));
      ylims = log(ylims);
    end
    yepts = vls(2)+xy*(xlims-vls(1));
    xepts = vls(1)+(ylims-vls(2))/xy;
    yinlims = yepts <= ylims(2) & yepts >= ylims(1);
    xinlims = xepts <= xlims(2) & xepts >= xlims(1);
    pts = [xlims(yinlims)' yepts(yinlims)';
	   xepts(xinlims)' ylims(xinlims)'];
    pts = repfun(pts,'remrep',1);
    if length(pts)==0
      if xlg xlims = exp(xlims); end
      if ylg ylims = exp(ylims); end
      h = plot(p,xlims(1),ylims(1));
      set(h,'visible','off')
    else
      if xlg pts(:,1) = exp(pts(:,1)); end
      if ylg pts(:,2) = exp(pts(:,2)); end
      h = plot(p,pts(:,1),pts(:,2));
    end
  elseif size(xy,2)==2

  end
  set(h,'linestyle',':','color','k')
  
  set(p,'nextplot',np)

end
